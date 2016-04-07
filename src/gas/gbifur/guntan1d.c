/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*
*				guntan1d.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Resolves untangles of interior curves.
*
*/

#if defined(ONED)

#include <gdecs/gdecs.h>


	/* LOCAL Function Declarations */
LOCAL	boolean    consistent_params1d(INTERFACE*);
LOCAL	int	untangle_interior_cross1d(Front*,CROSS*);
LOCAL	void	eliminate_point(POINT*);
LOCAL	void	oned_riemann_problem(Locstate,Locstate,Locstate,Locstate,
				     RIEMANN_SOLVER_WAVE_TYPE*,
				     RIEMANN_SOLVER_WAVE_TYPE*,double*);
LOCAL	void	propagate_point_past_dirichlet_boundary(POINT*,POINT*,SIDE);
LOCAL	void	reflect_point_at_neumann_boundary(POINT*,POINT*,SIDE,Front*);
LOCAL	void	untangle_boundary_interior_cross1d(Front*,CROSS*);

/*ARGSUSED*/
EXPORT int g_untangle_front1d(
	Front *front,
	CROSS **cross,
	int   flag)
{
	CROSS *cr;
	int status;

	if (cross==NULL || *cross==NULL)
	    return CURVES_UNTANGLED;

	for (cr = *cross; cr; cr = cr->next)
	{
	    status = untangle_interior_cross1d(front,cr);
	    if (status != CURVES_UNTANGLED)
		return status;
	    delete_from_cross_list(cr);
	    if (cr == *cross)
	        *cross = cr->next;
	}
	*cross = NULL;
	if (!consistent_params1d(front->interf))
	{
	    screen("ERROR in g_untangle_front1d(), inconsistent params\n");
	    (void) print_interface(front->interf);
	    clean_up(ERROR);
	    return ERROR_IN_UNTANGLE;
	}
	return CURVES_UNTANGLED;
}		/*end g_untangle_front1d*/

LOCAL	int	untangle_interior_cross1d(
	Front *front,
	CROSS *cr)
{
	COMPONENT                ncomp[7], pcomp[7];
	INTERFACE                *intfc = front->interf;
	RECT_GRID                *gr = front->rect_grid;
	POINT                    **p, *new_pt, *ptl, *ptr;
	POINT                    *before, *after;
	double                    new_coords[7][3], coords[3];
	double                    speed[17];
	double                    w_speed[7];
	double                    max_speed;
	double                    dt = front->dt;
	double                    dx = gr->h[0];
	double                    dt_frac;
	static double             nor[3] = {1.0, 0.0, 0.0};
	boolean                     track_new_wave[7];
	boolean                     reduce_time_step;
	int		         n_new_waves;
	int                      np, i;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	static const int w_type[7] = { BACKWARD_SOUND_WAVE_LE,
				       BACKWARD_SHOCK_WAVE,
				       BACKWARD_SOUND_WAVE_TE,
				       CONTACT,
				       FORWARD_SOUND_WAVE_TE,
				       FORWARD_SHOCK_WAVE,
				       FORWARD_SOUND_WAVE_LE };
	static const int w_class[7] = { RAREF_LEADING_EDGE,
					SHOCK_WAVE,
					RAREF_TRAILING_EDGE,
					CONTACT_WAVE,
					RAREF_TRAILING_EDGE,
					SHOCK_WAVE,
					RAREF_LEADING_EDGE };
	static Locstate  sl = NULL, sr = NULL, sl_mid = NULL, sr_mid = NULL;
	static Locstate  left[7], right[7], ahead[7], behind[7];

	debug_print("guntan","Entered untangle_interior_cross1d()\n");
	if (debugging("guntan"))
	{
	    (void) printf("Number of interaction points = cr->npt = %d\n",
			  cr->npt);
	    for (i = 0; i < cr->npt; ++i)
	    {
		(void) printf("Interacting point[%d]\n",i);
		print_point(cr->pt[i]);
	    }
	    (void) printf("Tangled interface\n");
	    print_interface(intfc);
	}
	 
	if (sl == NULL)
	{
	    alloc_state(intfc,&sl,front->sizest);
	    alloc_state(intfc,&sl_mid,front->sizest);
	    alloc_state(intfc,&sr_mid,front->sizest);
	    alloc_state(intfc,&sr,front->sizest);
	    left[0] = sl;      right[0] = sl;
	    left[1] = sl;      right[1] = sl_mid;
	    left[2] = sl_mid;  right[2] = sl_mid;
	    left[3] = sl_mid;  right[3] = sr_mid;
	    left[4] = sr_mid;  right[4] = sr_mid;
	    left[5] = sr_mid;  right[5] = sr;
	    left[6] = sr;      right[6] = sr;
	    ahead[0] = sl;     behind[0] = sl_mid;
	    ahead[1] = sl;     behind[1] = sl_mid;
	    ahead[2] = sl;     behind[2] = sl_mid;
	    ahead[3] = sl_mid; behind[3] = sr_mid;
	    ahead[4] = sr;     behind[4] = sr_mid;
	    ahead[5] = sr;     behind[5] = sr_mid;
	    ahead[6] = sr;     behind[6] = sr_mid;
	}

	ptl = cr->pt[0];
	ptr = cr->pt[cr->npt-1];
	if (debugging("guntan"))
	{
	    (void) printf("ptl\n");
	    print_point(ptl);
	    (void) printf("ptr\n");
	    print_point(ptr);
	}

	/* Remove intermediate interior points */
	for (i = 1; i < cr->npt-1; ++i)
	{
	    if (debugging("guntan"))
	    {
		(void) printf("Deleting intermediate interaction point\n");
		print_point(cr->pt[i]);
	    }
	    eliminate_point(cr->pt[i]);
	}

	p = intfc->points;
	np = intfc->num_points;
	for (before = NULL, i = 0; i < np-1; ++i)
	{
	    if ((p[i+1] == ptl) || (p[i+1] == ptr))
	    {
		before = p[i];
		break;
	    }
	}
	for (after = NULL, i = np-1; i > 0; i--)
	{
	    if ((p[i-1] == ptl) || (p[i-1] == ptr))
	    {
		after = p[i];
		break;
	    }
	}
	if ((before == NULL) || (after == NULL))
	{
	    screen("ERROR in untangle_interior_cross1d(), "
		   "crossing points not on interface\n");
	    clean_up(ERROR);
	}

	if (debugging("guntan"))
	{
	    (void) printf("Point to left of interaction\n");
	    print_point(before);
	    (void) printf("left incoming point\n");
	    print_point(ptl);
	    (void) printf("right incoming point\n");
	    print_point(ptr);
	    (void) printf("Point to right of interaction\n");
	    print_point(after);
	}

	coords[0] = 0.5*(Coords(ptl)[0] + Coords(ptr)[0]);
	if ((Coords(ptr)[0] < gr->L[0]) &&
		rect_boundary_type(intfc,0,0) == REFLECTION_BOUNDARY)
	{
	    copy_state(sr,right_state(ptr));
	    copy_state(sl,right_state(ptr));
	    nor[0] = 1.0;
	    reflect_state(sl,intfc,Coords(ptr),gr->L,nor);
	}
	else if ((Coords(ptl)[0] > gr->U[0]) &&
		rect_boundary_type(intfc,0,1) == REFLECTION_BOUNDARY)
	{
	    copy_state(sl,left_state(ptl));
	    copy_state(sr,left_state(ptl));
	    nor[0] = -1.0;
	    reflect_state(sr,intfc,Coords(ptl),gr->U,nor);
	}
	else
	{
	    copy_state(sl,left_state(ptl));
	    copy_state(sr,right_state(ptr));
	}
	oned_riemann_problem(sl,sr,sl_mid,sr_mid,&l_wave,&r_wave,speed);
	if (debugging("guntan"))
	{
	    (void) printf("Riemann Problem for Interaction\n");
	    (void) printf("Left/Right Data States\n");
	    verbose_print_state("sl",sl);
	    verbose_print_state("sr",sr);
	    (void) printf("Left/Right Answer Mid States\n");
	    verbose_print_state("sl_mid",sl_mid);
	    verbose_print_state("sr_mid",sr_mid);
	    (void) printf("l_wave = %s\n",rsoln_wave_name(l_wave));
	    (void) printf("r_wave = %s\n",rsoln_wave_name(r_wave));
	    (void) printf("dt = %g, dx = %g\n",dt,dx);
	    for (i = 0; i < 17; ++i)
		(void) printf("    speed[%d] = %g, speed[%d]*dt/dx = %g\n",
			      i,speed[i],i,speed[i]*dt/dx);
	}
	w_speed[0] = speed[3];
	w_speed[1] = speed[3];
	w_speed[2] = speed[4];
	w_speed[3] = speed[8];
	w_speed[4] = speed[12];
	w_speed[5] = speed[13];
	w_speed[6] = speed[13];
	dt_frac = 0.75;
	reduce_time_step = NO;
	for (max_speed = 0, i = 0; i < 17; ++i)
	{
	    double spd = fabs(speed[i]);
	    if (spd*dt > dx)
	    {
		dt_frac = min(dt_frac,dx/(spd*dt));
		reduce_time_step = YES;
	    }
	    if (max_speed < spd)
		max_speed = spd;
	}
	for (i = 0; i < 7; ++i)
	{
	    double ds;
	    new_coords[i][0] = coords[0] + dt*w_speed[i];
	    new_coords[i][1] = 0.0;
	    new_coords[i][2] = 0.0;
	    if (new_coords[i][0] <= Coords(before)[0])
	    {
		reduce_time_step = YES;
		ds = coords[0] - Coords(before)[0];
		dt_frac = min(dt_frac,fabs(ds/(dt*w_speed[i])));
	    }
	    if (Coords(after)[0] <= new_coords[i][0])
	    {
		reduce_time_step = YES;
		ds = Coords(after)[0] - coords[0];
		dt_frac = min(dt_frac,fabs(ds/(dt*w_speed[i])));
	    }
	}
	if (reduce_time_step)
	{
	    *front->dt_frac = min(dt_frac,*front->dt_frac);
	    debug_print("guntan","Left untangle_interior_cross1d()\n");
	    return MODIFY_TIME_STEP_TO_UNTANGLE;
	}

	set_max_front_speed(0,max_speed,return_obst_state(),coords,
			    front);
	set_max_front_speed(1,max_speed,return_obst_state(),coords,
			    front);

	n_new_waves = 0;
	for (i = 0; i < 7; ++i)
	    track_new_wave[i] = NO;
	if (!ComponentIsFlowSpecified(positive_component(before),front))
	{
	    if (l_wave == SHOCK)
	    {
	        if (tracked_oned_scattered_wave(w_class[1],ahead[1],
					        behind[1],front))
	        {
	            track_new_wave[1] = YES;
	            ++n_new_waves;
	        }
	    }
	    else if (l_wave == RAREFACTION)
	    {
	        if (tracked_oned_scattered_wave(w_class[0],ahead[0],
					        behind[0],front))
	        {
	            track_new_wave[0] = YES;
	            ++n_new_waves;
	        }
	        if (tracked_oned_scattered_wave(w_class[2],ahead[2],
					        behind[2],front))
	        {
	            track_new_wave[2] = YES;
	            ++n_new_waves;
	        }
	    }
	}
	if (tracked_oned_scattered_wave(w_class[3],ahead[3],
					behind[3],front))
	{
	    track_new_wave[3] = YES;
	    ++n_new_waves;
	}
	if (!ComponentIsFlowSpecified(negative_component(after),front))
	{
	    if (r_wave == SHOCK)
	    {
	        if (tracked_oned_scattered_wave(w_class[5],ahead[5],
					        behind[5],front))
	        {
	            track_new_wave[5] = YES;
	            ++n_new_waves;
	        }
	    }
	    else if (r_wave == RAREFACTION)
	    {
	        if (tracked_oned_scattered_wave(w_class[4],ahead[4],
					        behind[4],front))
	        {
	            track_new_wave[4] = YES;
	            ++n_new_waves;
	        }
	        if (tracked_oned_scattered_wave(w_class[6],ahead[6],
					        behind[6],front))
	        {
	            track_new_wave[6] = YES;
		    ++n_new_waves;
	        }
	    }
	}

	if (debugging("guntan"))
	    (void) printf("number of new tracked waves = %d\n",n_new_waves);
	if (n_new_waves == 0)
	{
	    if (debugging("guntan"))
		(void) printf("Untracking all waves from the interaction\n");
	    eliminate_point(ptl);
	    eliminate_point(ptr);
	    if (negative_component(after) != positive_component(before))
		set_equivalent_comps(negative_component(after),
			             positive_component(before),intfc);
	    negative_component(after) = positive_component(before);
	    debug_print("guntan","Left untangle_interior_cross1d()\n");
	    return CURVES_UNTANGLED;
	}
	else
	{
	    COMPONENT comp[8], new_comp;
	    int       j, k;
	    comp[0] = positive_component(before);
	    comp[7] = negative_component(after);
	    for (i = 0; i < 7; ++i)
	    {
		if (track_new_wave[i] == YES)
		{
		    for (j = i+1; j < 7; ++j)
			if (track_new_wave[j] == YES)
			    break;
		    if (j == 7) /* No further tracked waves ahead */
		    {
			for (k = i+1; k < 7; ++k)
			    comp[k] = comp[7];
			break;
		    }
		    else
		    {
			new_comp = new_component(NEW_COMP);
			for (k = i+1; k <= j; ++k)
			    comp[k] = new_comp;
			i = j - 1;
		    }
		}
		else
		    comp[i+1] = comp[i];
	    }
	    for (i = 0; i < 7; ++i)
	    {
		ncomp[i] = comp[i];
		pcomp[i] = comp[i+1];
	    }
	}
	if (debugging("guntan"))
	{
	    static const char *mesg[] = {
		"Tracking new backward rarefaction leading edge",
		"Tracking new backward shock wave",
		"Tracking new backward rarefaction trailing edge",
		"Tracking new contact",
		"Tracking new forward rarefaction trailing edge",
		"Tracking new forward shock wave",
		"Tracking new forward rarefaction leading edge" };
	    (void) printf("pcomp(before) = %d\n",positive_component(before));
	    (void) printf("ncomp[0] = %d\n",ncomp[0]);
	    for (i = 0; i < 7; ++i)
	    {
	        if (track_new_wave[i] == YES)
		    (void) printf("%s, ncomp = %d, pcomp = %d\n",
				  mesg[i],ncomp[i],pcomp[i]);
	    }
	    (void) printf("pcomp[6] = %d\n",pcomp[6]);
	    (void) printf("ncomp(after) = %d\n",negative_component(after));
	}
	for (i = 0; i < 7; ++i)
	{
	    if (track_new_wave[i] == YES)
	    {
		new_pt = make_point(new_coords[i],ncomp[i],pcomp[i]);
		wave_type(new_pt) = w_type[i];
		copy_state(left_state(new_pt),left[i]);
		copy_state(right_state(new_pt),right[i]);
	    }
	}
	eliminate_point(ptl);
	eliminate_point(ptr);

	if (debugging("guntan"))
	{
	    (void) printf("Untangled interface\n");
	    print_interface(intfc);
	}
	 
	debug_print("guntan","Left untangle_interior_cross1d()\n");
	return CURVES_UNTANGLED;
}		/*end untangle_interior_cross1d*/

/*
*			g_bdry_untangle1d():
*
*	Main driver for resolving one dimensional interactions between
*	boundaries and tracked interior waves.
*/

/*ARGSUSED*/
EXPORT int g_bdry_untangle1d(
	Front    *front,
	CROSS    **cross,
	RPROBLEM *rp,
	NODE     *n_shift,
	int      flag)
{
	CROSS *cr;
	int i;

	debug_print("guntan","Entered g_bdry_untangle1d()\n");
	if (cross==NULL || *cross==NULL)
	{
	    if (debugging("guntan"))
		(void) printf("No crossings nothing to do\n");
	    debug_print("guntan","Left g_bdry_untangle1d()\n");
	    return CURVES_UNTANGLED;
	}

	for (cr = *cross; cr; cr = cr->next)
	{
	    for (i = 0; i < cr->npt; ++i)
	    {
		if (is_passive_boundary(cr->pt[i]))
		{
	            if (debugging("guntan"))
		        (void) printf("Deleting cross with passive boundary\n");
		    delete_from_cross_list(cr);
		    if (cr == *cross)
			*cross = cr->next;
		    break;
		}
	    }
	}
	for (cr = *cross; cr; cr = cr->next)
	{
	    for (i = 0; i < cr->npt; ++i)
	    {
		if (wave_type(cr->pt[i]) < FIRST_PHYSICS_WAVE_TYPE)
		    break;
	    }
	    if (i == cr->npt)
	    {
	        if (debugging("guntan"))
		    (void) printf("Skiping non-boundary cross\n");
		continue;
	    }
	    untangle_boundary_interior_cross1d(front,cr);
	    delete_from_cross_list(cr);
	    if (cr == *cross)
	        *cross = cr->next;
	}

	debug_print("guntan","Left g_bdry_untangle1d()\n");
	return CURVES_UNTANGLED;
}		/*end g_bdry_untangle1d*/

/*
*		untangle_boundary_interior_cross1d():
*
*	Resolves the interactions between a set of interior waves and a
*	boundary.  The present version assumes only one boundary point
*	occurs in the interaction.	
*/

LOCAL	void	untangle_boundary_interior_cross1d(
	Front *front,
	CROSS *cr)
{
	POINT *bdry_pt;
	POINT *interior_pt;
	SIDE  bdry_side;
	int   i;

	debug_print("guntan","Entered untangle_boundary_interior_cross1d()\n");
	if (wave_type(cr->pt[0]) < FIRST_PHYSICS_WAVE_TYPE)
	{
	    bdry_pt = cr->pt[0];
	    interior_pt = cr->pt[cr->npt-1];
	    bdry_side = NEGATIVE_SIDE;
	}
	else if (wave_type(cr->pt[cr->npt-1]) < FIRST_PHYSICS_WAVE_TYPE)
	{
	    interior_pt = cr->pt[0];
	    bdry_pt = cr->pt[cr->npt-1];
	    bdry_side = POSITIVE_SIDE;
	}

	/* Remove intermediate interior points */
	for (i = 1; i < cr->npt-1; ++i)
	    eliminate_point(cr->pt[i]);
	if (bdry_side == NEGATIVE_SIDE)
	    positive_component(bdry_pt) = negative_component(interior_pt);
	else
	    negative_component(bdry_pt) = positive_component(interior_pt);
	
	switch (wave_type(bdry_pt))
	{
	case SUBDOMAIN_BOUNDARY:
	case REFLECTION_BOUNDARY:
	    /* These interactions are handled by scatter_front. */
	    return;
	case DIRICHLET_BOUNDARY:
	    propagate_point_past_dirichlet_boundary(bdry_pt,interior_pt,
						    bdry_side);
	    return;
	case NEUMANN_BOUNDARY:
	    reflect_point_at_neumann_boundary(bdry_pt,interior_pt,
					      bdry_side,front);
	    return;
	default:
	    screen("ERROR in untangle_boundary_interior_cross1d(), "
		   "unknown or unsupported boundary type\n");
	    clean_up(ERROR);
	}
	debug_print("guntan","Left untangle_boundary_interior_cross1d()\n");

}		/*end untangle_boundary_interior_cross1d*/

LOCAL	void propagate_point_past_dirichlet_boundary(
	POINT *bdry_pt,
	POINT *interior_pt,
	SIDE  bdry_side)
{
	Locstate  interior_state;
	debug_print("guntan","Entered propagate_point_past_dirichlet_boundary()\n");
	switch (bdry_side)
	{
	case NEGATIVE_SIDE:
	    positive_component(bdry_pt) = positive_component(interior_pt);
	    interior_state = right_state(interior_pt);
	    copy_state(right_state(bdry_pt),interior_state);
	    break;

	case POSITIVE_SIDE:
	    negative_component(bdry_pt) = negative_component(interior_pt);
	    interior_state = left_state(interior_pt);
	    copy_state(left_state(bdry_pt),interior_state);
	    break;

	default:
	    screen("ERROR in propagate_point_past_dirichlet_boundary(), "
		   "unknown value for bdry_side\n");
	    clean_up(ERROR);
	}
	eliminate_point(interior_pt);
	debug_print("guntan","Left propagate_point_past_dirichlet_boundary()\n");
}		/*end propagate_point_past_dirichlet_boundary*/


/*
*		reflect_point_at_neumann_boundary():
*
*	Resolve the interaction of a wave with a refecting boundary.
*	NOTE: The current code assumes the Neumann boundary is stationary.
*       TODO: Add support for moving boundaries.
*/

LOCAL	void reflect_point_at_neumann_boundary(
	POINT *bdry_pt,
	POINT *interior_pt,
	SIDE  bdry_side,
	Front *front)
{
	INTERFACE                *intfc = front->interf;
	Locstate                 sl, sr;
	double                    coords[3];
	double                    speed[17];
	double                    max_speed;
	double                    dt = front->dt;
	int                      i;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	static double             nor[3] = {1.0, 0.0, 0.0};
	static Locstate          scratch = NULL, sl_mid = NULL, sr_mid = NULL;

	debug_print("guntan","Entered reflect_point_at_neumann_boundary()\n");
	if (debugging("guntan"))
	{
	    (void) printf("Input data\n");
	    print_general_vector("bdry_pt = ",Coords(bdry_pt),
	                         front->rect_grid->dim,"\n");
	    print_general_vector("interior_pt = ",Coords(interior_pt),
	                         front->rect_grid->dim,"\n");
	    (void) printf("bdry_side = %s\n",side_name(bdry_side));
	}
	if (scratch == NULL)
	{
	    alloc_state(intfc,&scratch,front->sizest);
	    alloc_state(intfc,&sl_mid,front->sizest);
	    alloc_state(intfc,&sr_mid,front->sizest);
	}

	switch (bdry_side)
	{
	case NEGATIVE_SIDE:
	    sl = scratch;
	    sr = right_state(interior_pt);
	    coords[0] = Coords(bdry_pt)[0] + front->rect_grid->h[0];
	    copy_state(sl,sr);
	    reflect_state(sl,intfc,coords,Coords(bdry_pt),nor);
	    break;

	case POSITIVE_SIDE:
	    sl = left_state(interior_pt);
	    sr = scratch;
	    coords[0] = Coords(bdry_pt)[0] - front->rect_grid->h[0];
	    copy_state(sr,sl);
	    reflect_state(sr,intfc,coords,Coords(bdry_pt),nor);
	    break;

	default:
	    screen("ERROR in reflect_point_at_neumann_boundary(), "
		   "unknown value for bdry_side\n");
	    clean_up(ERROR);
	}
	oned_riemann_problem(sl,sr,sl_mid,sr_mid,&l_wave,&r_wave,speed);
	for (max_speed = 0, i = 0; i < 17; ++i)
	{
	    if (max_speed < fabs(speed[i]))
		max_speed = fabs(speed[i]);
	}
	set_max_front_speed(0,max_speed,return_obst_state(),Coords(bdry_pt),
			    front);
	set_max_front_speed(1,max_speed,return_obst_state(),Coords(bdry_pt),
			    front);

	if (bdry_side == NEGATIVE_SIDE)
	{
	    copy_state(right_state(bdry_pt),sr_mid);
	    if ((r_wave == SHOCK) && (tracked_oned_scattered_wave(SHOCK_WAVE,
							           sr,sr_mid,
								   front)))
	    {
	        Coords(interior_pt)[0] = Coords(bdry_pt)[0] + dt*speed[12];
	        copy_state(left_state(interior_pt),sr_mid);
	        copy_state(right_state(interior_pt),sr);
	        wave_type(interior_pt) = FORWARD_SHOCK_WAVE;
	    }
	    else
	    {
		positive_component(bdry_pt) = positive_component(interior_pt);
	        eliminate_point(interior_pt);
	    }
	}
	else /* bdry_side == POSITIVE_SIDE by above switch */
	{
	    copy_state(left_state(bdry_pt),sl_mid);
	    if ((l_wave == SHOCK) && (tracked_oned_scattered_wave(SHOCK_WAVE,
					                           sl,sl_mid,
								   front)))
	    {
	        Coords(interior_pt)[0] = Coords(bdry_pt)[0] + dt*speed[3];
	        copy_state(left_state(interior_pt),sl);
	        copy_state(right_state(interior_pt),sl_mid);
	        wave_type(interior_pt) = BACKWARD_SHOCK_WAVE;
	    }
	    else
	    {
		negative_component(bdry_pt) = negative_component(interior_pt);
	        eliminate_point(interior_pt);
	    }
	}
	if (debugging("guntan"))
	{
	    (void) printf("Output data\n");
	    print_general_vector("bdry_pt = ",Coords(bdry_pt),
	                         front->rect_grid->dim,"\n");
	    print_general_vector("interior_pt = ",Coords(interior_pt),
	                         front->rect_grid->dim,"\n");
	}
	debug_print("guntan","Left reflect_point_at_neumann_boundary()\n");
}		/*end reflect_point_at_neumann_boundary*/

/*
*			oned_riemann_problem():
*
*	Solves a one dimensional Riemann problem with data sl and sr.
*
*	Input:
*	       sl - left state
*	       sr - right state
*	Output: 
*	       sl_mid - left mid state
*	       sr_mid - right mid state
*	       l_wave - left wave family
*	       r_wave - right wave family
*/

LOCAL	void oned_riemann_problem(
	Locstate                 sl,
	Locstate                 sr,
	Locstate                 sl_mid,
	Locstate                 sr_mid,
	RIEMANN_SOLVER_WAVE_TYPE *l_wave,
	RIEMANN_SOLVER_WAVE_TYPE *r_wave,
	double                    *speed)
{
	double pjump = 0.0;
	double pmidl, pmidr;
	double umidl, umidr;
	double ml, mr;
	double ul, cl, uml, cml, umr, cmr, ur, cr;

	(void) find_mid_state(sl,sr,pjump,&pmidl,&pmidr,&umidl,&umidr,&ml,
			      &mr,l_wave,r_wave);
	Dens(sl_mid) = (*l_wave == SHOCK) ? dens_Hugoniot(pmidl,sl) :
					    dens_rarefaction(pmidl,sl);
	Vel(sl_mid)[0] = umidl;
	Press(sl_mid) = pmidl;
	Set_params(sl_mid,sl);
	set_type_of_state(sl_mid,TGAS_STATE);
	set_state(sl_mid,state_type(sl),sl_mid);

	Dens(sr_mid) = (*r_wave == SHOCK) ? dens_Hugoniot(pmidr,sr) :
					    dens_rarefaction(pmidr,sr);
	Vel(sr_mid)[0] = umidr;
	Press(sr_mid) = pmidr;
	Set_params(sr_mid,sr);
	set_type_of_state(sr_mid,TGAS_STATE);
	set_state(sr_mid,state_type(sr),sr_mid);

	ul  = vel(0,sl);      cl = sound_speed(sl);
	uml = vel(0,sl_mid); cml = sound_speed(sl_mid);
	umr = vel(0,sr_mid); cmr = sound_speed(sr_mid);
	ur  = vel(0,sr);      cr = sound_speed(sr);

	speed[0] = ul - cl; speed[1] = ul; speed[2] = ul + cl;
	if (*l_wave == SHOCK)
	{
	    speed[3] = speed[4] = ul - ml/Dens(sl);
	}
	else
	{
	    speed[3] = speed[0];
	    speed[4] = uml - cml;
	}
	speed[5] = uml - cml; speed[6] = uml; speed[7] = uml + cml;
	speed[8] = 0.5*(uml+umr);
	speed[9] = umr - cmr; speed[10] = umr; speed[11] = umr + cmr;
	if (*r_wave == SHOCK)
	{
	    speed[12] = speed[13] = ur + mr/Dens(sr);
	}
	else
	{
            speed[12] = speed[11];
	    speed[13] = ur + cr;
	}
	speed[14] = ur - cr; speed[15] = ur; speed[16] = ur + cr;
}		/*end oned_riemann_problem*/

LOCAL	void	eliminate_point(
	POINT* pt)
{
	INTERFACE *intfc = pt->interface;
	positive_component(pt) = NO_COMP;
	negative_component(pt) = NO_COMP;
	reset_intfc_components(intfc);
	(void) delete_point(pt);
}		/*end eliminate_point*/

LOCAL boolean consistent_params1d(
	INTERFACE *intfc)
{
	POINT **p;
	int   i, np;

	p = intfc->points;
	if (p == NULL)
	    return YES;
	np = intfc->num_points;
	for (i = 1; i < np; ++i)
	{
	    if ( !is_subdomain_boundary(Hyper_surf(p[i-1])) &&
	         !is_subdomain_boundary(Hyper_surf(p[i])) &&
	         (Params(right_state(p[i-1])) != Params(left_state(p[i]))) )
	        return NO;
	}
	return YES;
}		/*end consistent_components1d */

#endif /* defined(ONED) */
