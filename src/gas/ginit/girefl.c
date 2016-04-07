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
*				girefl.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for initializing and analyzing shock
*	reflection or refraction problems.
*/


#if defined(TWOD)
#include <ginit/ginit.h>


	/* LOCAL Function Declarations */
LOCAL	boolean	compute_anom(double,double*,POINTER);
LOCAL	boolean	compute_detachment(double,double*,POINTER);
LOCAL	boolean	compute_mech_eq(double,double*,POINTER);
LOCAL	boolean	compute_sonic(double,double*,POINTER);
LOCAL	boolean	compute_vnr(double,double*,POINTER);
LOCAL	void	compute_mach_angles(double,RP_DATA*);
LOCAL	void	free_bubble_comp_type(COMP_TYPE*);
LOCAL	void	get_state_bubble(double*,Locstate,COMP_TYPE*,
				 HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	init_mach_reflection(int,CURVE*,Bubble*,Front*);
LOCAL	void	init_regular_reflection(double*,int,CURVE*,Bubble*,Front*);
LOCAL	void	set_bubble_comp_type(COMP_TYPE*,Front*);

/*
*			init_ramp_reflection():
*
*	This function drives the initialization of reflections created by
*	a normal shock passing a sharp compressive corner in a wall.  The
*	result can be either a regular reflection, or a mach reflection,
*	depending on the ramp angle and incident shock strength.  The
*	geometry can either be NORMAL, meaning the incident is normal to the
*	grid, and the grid runs along the refl_wall.  We can also have
*	OBLIQUE geometry where the incident runs obliquely through the grid,
*	letting the grid align with the bow_wall.  We also allow for the
*	shock to be ahead or behind the corner (before or after bifurcation).
*	The NORMAL case allows for a second corner in the ramp, without
*	allowing the shock to initially pass this point.
*/

EXPORT void init_ramp_reflection(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	Bubble		*bubble;
	CURVE		*ramp = NULL;
	CURVE		*incident;
	CURVE		**curves;
	Gas_param	*params;
	NODE		*ns, *ne;
	POINT		*p;
	RECT_GRID	*rect_grid = front->rect_grid;
	char		s[Gets_BUF_SIZE];
	double		*L = rect_grid->L;
	double		*U = rect_grid->U;
	double		angle1;			/* angle of first corner */
	double		angle2;			/* angle of second corner */
	double		thickness;		/* height of second corner */
	double		dist_from_left;		/* x posn of first corner */
	double		shock_corner_dist;	/* if inc not past corner */
	double		rho0, p0; 		/* ahead data */
	double		M0;			/* inc shock data */
	double		inc_ang;		/* from pos x axis */
	double		inc_n[MAXD];		/* inc normal to ahead */
	double		ramp_height;
	double		inc_speed;
	double		node_v[MAXD];
	double		corner_posn[MAXD];
	double		ramp_top[MAXD];
	double		coords0[MAXD];
	double		coords1[MAXD];
	int		constant_ahead, constant_behind;
	int		shock_past_corner;
	int		prob_type;		/* OBLIQUE or NORMAL */
	int		i, dim = rect_grid->dim;
	boolean		sav_scss = interpolate_states_at_split_curve_node();

	screen("Enter the incident orientation, choices are:\n");
	screen("\toblique (O), normal (N), or double (D): ");
	(void) Gets(s);
	switch (s[0])
	{
	case 'N':
	case 'n':
	    prob_type = NORMAL;
	    screen("Enter the ramp angles (degrees),\n");
	    screen("\tthickness, and distance to inlet: ");
	    (void) Scanf("%f %f %f %f\n",&angle1,&angle2,
			 &thickness,&dist_from_left);
	    inc_ang = PI / 2.;
	    angle1 = radians(angle1);	angle2 = radians(angle2);

	    make_ramp(prob_type,angle1,angle2,thickness,dist_from_left,
		      COMPOBST,COMPA,rect_grid,&ramp);

	    corner_posn[0] = Coords(ramp->end->posn)[0];
	    corner_posn[1] = Coords(ramp->end->posn)[1];
	    break;
	case 'O':
	case 'o':
	    prob_type = OBLIQUE;
	    screen("Enter the incident angle and the height of the ramp: ");
	    (void) Scanf("%f %f\n",&inc_ang,&ramp_height);
	    inc_ang = radians(inc_ang);

	    ramp_top[0] = L[0];
	    ramp_top[1] = L[1] + ramp_height;

	    corner_posn[0] = L[0] + ramp_height * tan(inc_ang);
	    corner_posn[1] = L[1];
		
			/* make ramp */

	    ns = make_node(Point(corner_posn));
	    ne = make_node(Point(ramp_top));
	    node_type(ne) = FIXED_NODE;
	    ramp = make_curve(COMPOBST,COMPB,ns,ne);
	    wave_type(ramp) = NEUMANN_BOUNDARY;
	    start_status(ramp) = FIXED;
	    end_status(ramp) = FIXED;
	    break;
	case 'D':
	case 'd':
	    prob_type = DOUBLE_REFL;
	    screen("Enter the ramp angles and distance to inlet: ");
	    (void) Scanf("%f %f %f\n",&angle1,&angle2,&dist_from_left);
	    inc_ang = PI / 2.;
	    angle1 = radians(angle1);	angle2 = radians(angle2);

	    make_ramp(DOUBLE_REFL,angle1,angle2,thickness,dist_from_left,
		      COMPOBST,COMPA,rect_grid,&ramp);

		/*Note:  the corner here is not the same as above*/
	    corner_posn[0] = Coords(ramp->start->posn)[0];
	    corner_posn[1] = Coords(ramp->start->posn)[1];

	    screen("ERROR in init_ramp_reflection(), "
	           "DOUBLE_REFL geometry - CODE NEEDED\n");
	    clean_up(ERROR);
	    break;
        case 'A':
	case 'a':
	    prob_type = MACH_ANGLES;
	    break;	       
	default:
	    screen("ERROR in init_ramp_reflection(), Invalid choice %s\n",s);
	    clean_up(ERROR);
	}
	ramp_reflection_corner_posn(corner_posn,YES,dim);

	screen("Type 'y' if the incident shock has passed the corner: ");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	{
	    shock_past_corner = YES;
	    set_bubble_comp_type(comp_type(COMPBOW),front);
	    bubble = (Bubble*)comp_type(COMPBOW)->extra;

	    bubble->comp_bow = COMPBOW;
	    bubble->is_node_at_corner = YES;

	    screen("Enter the distance travelled past the corner: ");
	    (void) Scanf("%f\n",&bubble->refl_length);
	}
	else
	{
	    screen("Enter the distance from the shock to the corner: ");
	    (void) Scanf("%f\n",&shock_corner_dist);
	    shock_past_corner = NO;
	}

	params = init_eos_params(init,ip,"",YES);
	screen("Enter the ahead state -- rho, pressure: ");
	(void) Scanf("%f %f\n",&rho0,&p0);
	screen("Enter the incident Mach number: ");
	(void) Scanf("%f\n",&M0);

	set_ambient_comp_type(comp_type(COMPA),front);

	screen("Is the flow ahead of the wave constant (dflt = no): ");
	(void) Gets(s);
	constant_ahead = (s[0] == 'y' || s[0] == 'Y') ? YES : NO;

	set_ambient_comp_type(comp_type(COMPB),front);

	screen("Is the flow behind the wave constant (dflt = no): ");
	(void) Gets(s);
	constant_behind = (s[0] == 'y' || s[0] == 'Y') ? YES : NO;

	set_obstacle_comp_type(comp_type(COMPOBST),front);

	Init_params(Ambient(comp_type(COMPA)),params);
	Dens(Ambient(comp_type(COMPA))) = rho0;
	Vel(Ambient(comp_type(COMPA)))[0] = 0.;
	Vel(Ambient(comp_type(COMPA)))[1] = 0.;
	Press(Ambient(comp_type(COMPA))) = p0;
	set_type_of_state(Ambient(comp_type(COMPA)),TGAS_STATE);
	set_state(Ambient(comp_type(COMPA)),GAS_STATE,
		  Ambient(comp_type(COMPA)));

	if (constant_ahead)
	    (void)SetConstantFlowRegion(COMPA,Ambient(comp_type(COMPA)),
					front->interf);
	if (prob_type == NORMAL)
	{
	    inc_n[0] = 1.0;
	    inc_n[1] = 0.0;
	}
	else
	{
	    inc_n[0] = cos(inc_ang - PI/2.);
	    inc_n[1] = sin(inc_ang - PI/2.);
	}
	(void) s_polar_4(SHOCK_MACH_NUMBER,M0,&inc_speed,inc_n,
		         Ambient(comp_type(COMPA)),Ambient(comp_type(COMPB)),
			 GAS_STATE);
	if (constant_behind)
	    (void)SetConstantFlowRegion(COMPB,Ambient(comp_type(COMPB)),
					front->interf);
	if (shock_past_corner)
	{
	    bubble->RP->ang_dir = COUNTER_CLOCK;
	    bubble->RP->ang[1] = inc_ang;
	    for (i = 0; i < dim; i++)
		bubble->cor_v[i] = 0.0;

	    ft_assign(bubble->RP->state[0],Ambient(comp_type(COMPA)),
	    	front->sizest);
	    ft_assign(bubble->RP->state[1],Ambient(comp_type(COMPB)),
	    	front->sizest);

	    if (prob_type == MACH_ANGLES)
	    	compute_mach_angles(inc_speed,bubble->RP);

	    for (i = 0; i < dim; i++)
	    	bubble->corner_posn[i] = corner_posn[i];

	    if (prob_type == NORMAL)
	    {
	    	bubble->inc_t[0] = 0.0;
	    	bubble->inc_t[1] = 1.0;

	    	bubble->aw_t[0] = cos(angle1);
	    	bubble->aw_t[1] = sin(angle1);

	    	bubble->bw_t[0] = -1.0;
	    	bubble->bw_t[1] = 0.0;

	    	node_v[0] = inc_speed;
	    	node_v[1] = inc_speed * bubble->aw_t[1] / bubble->aw_t[0];
	    }
	    else
	    {
	    	bubble->inc_t[0] = cos(inc_ang);
	    	bubble->inc_t[1] = sin(inc_ang);

	    	bubble->aw_t[0] = 1.0;
	    	bubble->aw_t[1] = 0.0;

	    	bubble->bw_t[0] = cos(inc_ang + PI/2.0);
	    	bubble->bw_t[1] = sin(inc_ang + PI/2.0);

	    	node_v[0] = inc_speed / bubble->inc_t[1];
	    	node_v[1] = inc_speed * 0.0;
	    }

	    if (is_regular_reflection(node_v,front,bubble->RP))
	    {
	    	if (prob_type == NORMAL)
	    	{
	    	    bubble->RP->ang[0] = angle1;
	    	    bubble->RP->ang[3] = angle1 + PI;
	    	}
		else
		{
		    bubble->RP->ang[0] = 0.0;
		    bubble->RP->ang[3] = inc_ang + PI/2.;
	    	}
	    	init_regular_reflection(node_v,prob_type,ramp,bubble,front);
	    }
	    else
	    {
	    	/* Note:  the angles are labelled differently for
	    	 * reg vs mach nodes.  Compare B_reflect_node_propagate()
	    	 * with Mach_node_propatate().
		 */

	    	bubble->RP->ang[0] = bubble->RP->ang[1];
	    	init_mach_reflection(prob_type,ramp,bubble,front);
	    }
	}
	else
	{
	    if (prob_type == NORMAL)
	    {
	    	coords0[0] = coords1[0] = 
	    		L[0] + dist_from_left - shock_corner_dist;
	    	coords0[1] = L[1];	coords1[1] = U[1];
	    	ns = make_node(Point(coords0));
	    	ne = make_node(Point(coords1));
	    	incident = make_curve(COMPB,COMPA,ns,ne);
	    }
	    else
	    {
	    	coords0[0] = corner_posn[0] +
	    		shock_corner_dist * cos(inc_ang + PI/2.);
	    	coords0[1] = corner_posn[1] +
	    		shock_corner_dist * sin(inc_ang + PI/2.);

	    	coords1[0] = coords0[0] + (U[1] - coords0[1]) / tan(inc_ang);
	    	if (coords1[0] > U[0])
	    	{	/* end inc on right wall, not upper wall */
	    		coords1[0] = U[0];
	    		coords1[1] = coords0[1] +
	    			(U[0] - coords0[0]) * tan(inc_ang);
	    	}
	    	else
	    		coords1[1] = U[1];

	    	p = Point(coords0);
	    	set_interpolate_states_at_split_curve_node(NO);
	    	curves = split_curve(p,ramp->last,ramp,
	    		COMPOBST,COMPA,COMPOBST,COMPB);
	    	set_interpolate_states_at_split_curve_node(sav_scss);
	    	end_status(curves[0]) = FIXED;
	    	start_status(curves[1]) = FIXED;
	    	ns = curves[0]->end;
	    	ne = make_node(Point(coords1));
	    	incident = make_curve(COMPB,COMPA,ns,ne);
	    }
	    wave_type(incident) = FORWARD_SHOCK_WAVE;
	    start_status(incident) = INCIDENT;
	    end_status(incident) = INCIDENT;

	    if (debugging("iramp"))
	    {
	        (void) printf("not past corner in init_ramp_reflection()\n");
	        verbose_print_state("ahead state",Ambient(comp_type(COMPA)));
		verbose_print_state("behind state",Ambient(comp_type(COMPB)));
	    }
	}

	/* circle_D_extend() causes problems if the incident shock passes
	*  a corner at the end opposite from the reflection. For example
	*  in the oblique geometry, it can pass the upper right corner, and
	*  will no longer necessarily be straight.
	*/
	set_use_circle_D_extend(NO);

	if (debugging("iramp"))
	{
	   print_interface(front->interf);
	}
	if (debugging("Mach_stat"))
	{
	    RP_DATA *RP = bubble->RP;

	    (void) printf("\nStatic Mach node state stats:\n");
	    (void) printf("theta_w = %g\n",degrees(angle1));
	    (void) printf("M_0 = %g\n",M0);
	    (void) printf("corner[0] = %g, corner[1] = %g\n",
		          corner_posn[0],corner_posn[1]);
	    verbose_print_state("RP->state[0]",RP->state[0]);
	    verbose_print_state("RP->state[1]",RP->state[1]);
	    (void) printf("\n");
	}
}		/*end init_ramp_reflection*/


LOCAL	void init_regular_reflection(
	double		*node_v,
	int		prob_type,
	CURVE		*ramp,
	Bubble		*bubble,
	Front		*front)
{
	CURVE		*bow;
	CURVE		*incident;
	CURVE		**curves;
	NODE		*ns, *ne;
	POINT		*p;
	double		len;			/* refl to bow */
	double		end_inc[MAXD];
	double		mid_posn[MAXD];		/* between refl and bow */
	double		refl_ang = bubble->RP->ang[2];
	double		*refl_length = &bubble->refl_length;
	double		*bow_length = &bubble->bow_length;
	double		*refl_posn = bubble->refl_posn;
	double		*bow_posn = bubble->bow_posn;
	double		*corner_posn = bubble->corner_posn;
	double		*inc_t = bubble->inc_t;
	double		*aw_t = bubble->aw_t;
	double		*bw_t = bubble->bw_t;
	double		*U = front->rect_grid->U;
	int		*gmax = front->rect_grid->gmax;
	boolean		sav_scss = interpolate_states_at_split_curve_node();

	init_bubble(node_v,REFL_WALL_LENGTH,bubble->refl_length,bubble,front);
	
	if (prob_type == NORMAL)
	{
		refl_posn[0] = corner_posn[0] + *refl_length*aw_t[0];
		refl_posn[1] = corner_posn[1] + *refl_length*aw_t[1];

		end_inc[0] = refl_posn[0];
		end_inc[1] = U[1];

		bow_posn[0] = corner_posn[0] - *bow_length;
		bow_posn[1] = corner_posn[1];
	}
	else
	{
		refl_posn[0] = corner_posn[0] + *refl_length;
		refl_posn[1] = corner_posn[1];

		end_inc[0] = refl_posn[0] + 
			(U[1] - refl_posn[1]) * inc_t[0] / inc_t[1];
		if (end_inc[0] > U[0])
		{		/* end_inc on right wall, not upper wall */
			end_inc[0] = U[0];
			end_inc[1] = refl_posn[1] + 
				(U[0] - refl_posn[0]) * inc_t[1] / inc_t[0];
		}
		else
			end_inc[1] = U[1];

		bow_posn[0] = corner_posn[0] + *bow_length*bw_t[0];
		bow_posn[1] = corner_posn[1] + *bow_length*bw_t[1];
	}

		/* Make incident */

	if (prob_type == NORMAL)
	{
		p = Point(bubble->refl_posn);
		set_interpolate_states_at_split_curve_node(NO);
		curves = split_curve(p,ramp->last,ramp,
			COMPOBST,COMPA,COMPOBST,COMPBOW);
		set_interpolate_states_at_split_curve_node(sav_scss);
		end_status(curves[0]) = start_status(curves[1]) = FIXED;
		
		ns = curves[1]->start;
	}
	else
	{
		ns = make_node(Point(refl_posn));
	}
	node_type(ns) = B_REFLECT_NODE;
	copy_RP_DATA_structure(Rp_data(ns),bubble->RP);
	ne = make_node(Point(end_inc));
	incident = make_curve(COMPB,COMPA,ns,ne);
	wave_type(incident) = FORWARD_SHOCK_WAVE;
	start_status(incident) = end_status(incident) = INCIDENT;

		/* Make reflected shock */

	if (prob_type == NORMAL)
	{
		if (bubble->is_attached)
		{
			ne = curves[1]->end;
			node_type(ne) = ATTACHED_B_NODE;
		}
		else
		{
			ne = make_node(Point(bow_posn));
			/* no node type since boundary not created yet */
		}
	}
	else
	{
		if (bubble->is_attached)
		{
			ne = ramp->start;
			node_type(ne) = ATTACHED_B_NODE;
		}
		else
		{
			p = Point(bow_posn);
			set_interpolate_states_at_split_curve_node(NO);
			curves = split_curve(p,ramp->last,ramp,
					     COMPOBST,COMPBOW,COMPOBST,COMPB);
			set_interpolate_states_at_split_curve_node(sav_scss);
			end_status(curves[0]) = 
				start_status(curves[1]) = FIXED;

			ne = curves[0]->end;
			node_type(ne) = NEUMANN_NODE;
		}
	}
	len = hypot(refl_posn[0] - bow_posn[0],refl_posn[1] - bow_posn[1]);
	mid_posn[0] = refl_posn[0] + .5 * len * cos(refl_ang);
	mid_posn[1] = refl_posn[1] + .5 * len * sin(refl_ang);
	bow = make_half_bow_curve(COMPBOW,COMPB,
		2*max(gmax[0],gmax[1]),POSITIVE_ORIENTATION,
		ns,ne,mid_posn,front->rect_grid,bubble);
	wave_type(bow) = FORWARD_SHOCK_WAVE;
	start_status(bow) = REFLECTED;
	end_status(bow) = INCIDENT;

}		/*end init_regular_reflection*/


LOCAL void init_mach_reflection(
	int		prob_type,
	CURVE		*ramp,
	Bubble		*bubble,
	Front		*front)
{
	CURVE		*bow;
	CURVE		*incident;
	CURVE		*cur;
	CURVE		**curves;
	NODE		*ns, *ne;
	POINT		*p;
	double		surf_ten;
	double		mid_posn[MAXD];		/* between pt and bow */
	double		end_inc[MAXD];
	double		*bow_length = &bubble->bow_length;
	double		*mach_height = &bubble->mach_height;
	double		*refl_length = &bubble->refl_length;
	double		slip_length;		/* corner to base of contact */
	double		*refl_ang = &bubble->RP->ang[1];
	double		*contact_ang = &bubble->RP->ang[2];
	double		*mach_ang = &bubble->RP->ang[3];
	double		*corner_posn = bubble->corner_posn;
	double		*refl_posn = bubble->refl_posn;
	double		*bow_posn = bubble->bow_posn;
	double		*slip_posn = bubble->slip_posn;
	double		*mach_posn = bubble->mach_posn;
	double		*inc_t = bubble->inc_t;
	double		*aw_t = bubble->aw_t;
	double		*bw_t = bubble->bw_t;
	double		*L = front->rect_grid->L;
	double		*U = front->rect_grid->U;
	int		*gmax = front->rect_grid->gmax;
	boolean		sav_scss = interpolate_states_at_split_curve_node();

	screen("In Mach reflection case,\n\t");
	surf_ten = prompt_for_surface_tension(CONTACT,
			"for the slip line ");

	bubble->comp_mach = COMPMACH;
	comp_type(COMPMACH)->extra = (POINTER) bubble;
	set_bubble_comp_type(comp_type(COMPMACH),front);

	if (!init_mach_bubble(REFL_WALL_LENGTH,bubble->refl_length,
				 bubble,NORMAL_TO_MACH_REFLECTION,front))
	{
		screen("ERROR in init_mach_reflection(), ");
		screen("init_mach_bubble() failed\n");
		clean_up(ERROR);
	}

	if (prob_type == NORMAL)
	{
		mach_posn[0] = corner_posn[0] + *refl_length*aw_t[0];
		mach_posn[1] = corner_posn[1] + *refl_length*aw_t[1];

		refl_posn[0] = mach_posn[0] - *mach_height*aw_t[1];
		refl_posn[1] = mach_posn[1] + *mach_height*aw_t[0];

		bow_posn[0] = corner_posn[0] - *bow_length;
		bow_posn[1] = corner_posn[1];

		slip_length = *refl_length - 
			*mach_height*tan(*mach_ang - *contact_ang);
		slip_posn[0] = corner_posn[0] + slip_length*aw_t[0];
		slip_posn[1] = corner_posn[1] + slip_length*aw_t[1];

		end_inc[0] = refl_posn[0];
		end_inc[1] = U[1];
	}
	else
	{
		mach_posn[0] = corner_posn[0] + *refl_length;
		mach_posn[1] = L[1];

		refl_posn[0] = mach_posn[0];
		refl_posn[1] = mach_posn[1] + *mach_height;

		bow_posn[0] = corner_posn[0] + *bow_length*bw_t[0];
		bow_posn[1] = corner_posn[1] + *bow_length*bw_t[1];

		slip_length = *refl_length - 
			*mach_height*tan(*mach_ang - *contact_ang);
		slip_posn[0] = corner_posn[0] + slip_length;
		slip_posn[1] = L[1];

		end_inc[0] = refl_posn[0] +
			(U[1] - refl_posn[1]) * inc_t[0] / inc_t[1];
		if (end_inc[0] > U[0])
		{		/* end inc on right wall, not upper wall */
			end_inc[0] = U[0];
			end_inc[1] = refl_posn[1] +
				(U[0] - refl_posn[0]) * inc_t[1] / inc_t[0];
		}
		else
			end_inc[1] = U[1];
	}

	mid_posn[0] = refl_posn[0] + 0.5 * *mach_height * cos(*refl_ang);
	mid_posn[1] = refl_posn[1] + 0.5 * *mach_height * sin(*refl_ang);

		/* Make incident shock */

	ns = make_node(Point(refl_posn));
	node_type(ns) = MACH_NODE;
	copy_RP_DATA_structure(Rp_data(ns),bubble->RP);
	ne = make_node(Point(end_inc));
	incident = make_curve(COMPB,COMPA,ns,ne);
	wave_type(incident) = FORWARD_SHOCK_WAVE;
	start_status(incident) = INCIDENT;
	end_status(incident) = INCIDENT;

		/* Make reflected shock */
	if (bubble->is_attached)
	{
		ne = ramp->start;
		node_type(ne) = ATTACHED_B_NODE;
	}
	else if (prob_type == NORMAL)
	{
		ne = make_node(Point(bow_posn));
		/* no node type since boundary curves not created yet */
	}
	else
	{
		p = Point(bubble->bow_posn);
		set_interpolate_states_at_split_curve_node(NO);
		curves = split_curve(p,ramp->last,ramp,
				     COMPOBST,COMPBOW,COMPOBST,COMPB);
		set_interpolate_states_at_split_curve_node(sav_scss);
		end_status(curves[0]) = start_status(curves[1]) = FIXED;
		ramp = curves[1];
		ne = curves[0]->end;
		node_type(ne) = NEUMANN_NODE;
	}
	bow = make_half_bow_curve(COMPBOW,COMPB,2*max(gmax[0],gmax[1]),
		POSITIVE_ORIENTATION,ns,ne,mid_posn,front->rect_grid,bubble);
	wave_type(bow) = FORWARD_SHOCK_WAVE;
	start_status(bow) = REFLECTED;
	end_status(bow) = INCIDENT;

		/* Make mach stem */

	ns = incident->start;
	if (prob_type == NORMAL)
	{
		p = Point(bubble->mach_posn);
		set_interpolate_states_at_split_curve_node(NO);
		curves = split_curve(p,ramp->last,ramp,
			COMPOBST,COMPA,COMPOBST,COMPMACH);
		set_interpolate_states_at_split_curve_node(sav_scss);
		end_status(curves[0]) = start_status(curves[1]) = FIXED;
		ramp = curves[1];
		ne = curves[1]->start;
	}
	else
	{
		ne = make_node(Point(mach_posn));
	}
	cur = make_curve(COMPA,COMPMACH,ns,ne);
	wave_type(cur) = BACKWARD_SHOCK_WAVE;
	start_status(cur) = MACH_STEM;
	end_status(cur) = INCIDENT;

		/* Make contact */

	ns = incident->start;
	if (prob_type == NORMAL)
	{
		p = Point(bubble->slip_posn);
		set_interpolate_states_at_split_curve_node(NO);
		curves = split_curve(p,ramp->last,ramp,
			COMPOBST,COMPMACH,COMPOBST,COMPBOW);
		set_interpolate_states_at_split_curve_node(sav_scss);
		end_status(curves[0]) = start_status(curves[1]) = FIXED;
		ramp = curves[1];
		ne = curves[1]->start;
	}
	else
	{
		ne = make_node(Point(slip_posn));
	}
	cur = make_curve(COMPMACH,COMPBOW,ns,ne);
	wave_type(cur) = CONTACT;
	start_status(cur) = SLIP;
	end_status(cur) = INCIDENT;
	surface_tension(cur) = surf_ten;

}		/*end init_mach_reflection*/

typedef struct {
	double		inc_speed;
	double		a, b, c;
	double		p_detach;
	RP_DATA		*RP;
} MACH_ANG_PARAMS;

/*
*			compute_mach_angles():
*
*	This is a driver for the following series of functions.  Note that
*	use is made of the ordering of work done in this function, even
*	though each step doesn't technically follow from the previous.
*
*	anom_ang -- the transition point to anomolous reflection.  p1 has
*		reached the sonic point on the incident polar, so the three
*		shock solution ceases to exist here.
*	detach_ang -- the theoretical limit of regular reflection.  At this
*		point, the reflected polar is tangent to the vertical axis
*		in the pressure-turning angle shock polar diagram.
*	mech_eq_ang -- the mechanical equilibrium point is where the 
*		reflected polar intersects the vertical axis at the max
*		pressure on the incident polar.
*	sonic_ang -- this where the reflected polar intersects the vertical
*		axis at its sonic point.
*
*	In the second part of the routine, theoretical sets of various
*	quantities are computed for the three shock configuration.
*	The following are all in the steady frame, where the node is
*	stationary.
**	th01, th12, th13 -- the turning angles across the incident, reflecte
*		and Mach shocks respectively
*	q0 -- the ahead state speed in the steady frame
*	beta1 -- the angle between the reflected and its behind streamline
*	betas -- the angle between the Mach stem and its behind streamline
*	omega0 -- the angle between the incident and the incoming streamline
*	omega1 -- the angle between the reflected shock and its incoming
*		streamline
*	Omega1 -- the angle between the incident and reflected shocks
*	omegas -- the angle between the Mach stem and the incoming streamline
*	omega1p -- the angle between the reflected shock and the incoming
*		streamline ahead of the incident
*	chi -- the node trajectory angle relative to the wall
*	thetaw -- the ramp angle
*	a, b, c -- the coefficients in a numerically determined quadratic
*		relation between thetaw and chi -- input by user
*/


LOCAL void compute_mach_angles(
	double		inc_speed,
	RP_DATA		*RP)
{
	MACH_ANG_PARAMS mach_ang_params;
	double		epsilon, delta;
	double		min_ang, max_ang;
	double		anom_ang;
	double		detach_ang;
	double		mech_eq_ang;
	double		sonic_ang;
	double		vnr_ang;
	double		a, b, c;
	double		step;
	double		p3, minp3, maxp3;
	double		q0, q1, node_v[MAXD];
	double		m01, m12, m03;		/* mass fluxes */
	double		th01, th12, th03;	/* turn angles */
	double		omega0, omegas;
	double		omega1, Omega1, omega1p;
	double		beta1, betas;
	double		thetaw, chi;
	double		Mnsq, p1;

	mach_ang_params.inc_speed = inc_speed;
	mach_ang_params.RP = RP;
	min_ang = PI / 72.0; /*should be okay except for VERY strong shocks*/
	max_ang = PI / 2.0;
	delta = epsilon = 0.00001;

	if (find_root(compute_anom,(POINTER)&mach_ang_params,0.,&anom_ang,
		      min_ang,max_ang,epsilon,delta) == FUNCTION_FAILED)
	{
		screen("ERROR in compute_mach_angles(), ");
		(void) printf("find_root() failed\n");
		clean_up(ERROR);
	}
	(void) printf("\n\n\n");
	print_angle("anom_angle = ",anom_ang,"\n");
	max_ang = anom_ang;

	if (find_root(compute_detachment,(POINTER)&mach_ang_params,0.,
		&detach_ang,min_ang,max_ang,epsilon,delta) == FUNCTION_FAILED)
	{
		screen("ERROR in compute_mach_angles(), ");
		(void) printf("bisection_find_root() failed\n");
		clean_up(ERROR);
	}
	print_angle("detachment_angle = ",detach_ang,"\n");
	mach_ang_params.p_detach = pressure(RP->state[2]);

	if (find_root(compute_sonic,(POINTER)&mach_ang_params,0.,&sonic_ang,
		      min_ang,detach_ang,epsilon,delta) == FUNCTION_FAILED)
	{
		screen("ERROR in compute_mach_angles(), ");
		(void) printf("find_root() failed\n");
		clean_up(ERROR);
	}
	print_angle("sonic_angle = ",sonic_ang,"\n");

	if (find_root(compute_mech_eq,(POINTER)&mach_ang_params,0.,&mech_eq_ang,
		      min_ang,max_ang,epsilon,delta) == FUNCTION_FAILED)
	{
		screen("ERROR in compute_mach_angles(), ");
		(void) printf("find_root() failed\n");
		clean_up(ERROR);
	}
	print_angle("mech_eq_angle = ",mech_eq_ang,"\n");
	maxp3 = pressure(RP->state[2]) - EPS;

	if (find_root(compute_vnr,(POINTER)&mach_ang_params,0.,&vnr_ang,
		      mech_eq_ang,max_ang,epsilon,delta) == FUNCTION_FAILED)
	{
		screen("ERROR in compute_mach_angles(), ");
		(void) printf("find_root() failed\n");
		clean_up(ERROR);
	}
	print_angle("vnr_angle = ",vnr_ang,"\n");
	minp3 = pressure(RP->state[3]) + EPS;

	screen("Input the coefficients a, b, and c of the\n");
	screen("\tquadratic for chi in terms of thetaw: ");
	(void) Scanf("%f %f %f\n",&a,&b,&c);

	(void) printf("\n\n%10s %10s %10s %10s %10s %10s %10s\n",
	       "omega0","thetaw","chi","omega1p","Omega1","beta1","betas");
	step = detach_ang/50.0;
	chi = beta1 = betas = 0.0;	/* not relevant for regular refl */
	RP->ang[1] = PI/2.0;
	for (omega0 = step; omega0 <= detach_ang; omega0 += step)
	{
		node_v[0] = inc_speed;
		node_v[1] = inc_speed / tan(omega0);

		Mnsq = inc_speed*inc_speed/sound_speed_squared(RP->state[0]);
		p1 = max_behind_shock_pr(Mnsq,RP->state[0]);

		Check_return(
		    s_polar_3(RP->state[0],YES,p1,YES,NO,node_v,RP->state[1],
			      &RP->ang[0],&th01),
		    compute_mach_angles)

		th12 = -th01;
		if (!s_polar_2(RP->state[1],YES,YES,
			  th12,node_v,RP->state[2],&RP->ang[2]))
		{
			(void) printf("s_polar_2() failed");
			(void) printf(" for omega0 = %g\n",omega0);
			continue;
		}

		thetaw = PI/2.0 - omega0;
		Omega1 = RP->ang[2] - RP->ang[1];
		omega1p = PI - Omega1 - omega0;

		print_line_of_floats(7,omega0,thetaw,chi,omega1p,
				     Omega1,beta1,betas);
	}

	step = (maxp3 - minp3)/50.;
	for (p3 = minp3; p3 <= maxp3; p3 += step)
	{
		if (!find_steady_ahead_speed(CLOCKWISE,RP->state[0],
			RP->state[1],RP->state[3],p3,&q0,&th01,&th12,&th03))
		{
			(void) printf("find_steady_ahead_speed() failed");
			(void) printf(" for p3 = %g\n",p3);
			continue;
		}

		node_v[0] = inc_speed;
		node_v[1] = sqrt(q0*q0 - inc_speed*inc_speed);

		Check_return(
		    s_polar_3(RP->state[0],YES,pressure(RP->state[1]),YES,NO,
			      node_v,RP->state[1],&RP->ang[0],&th01),
		    compute_mach_angles)
 
		Check_return(
		    s_polar_3(RP->state[1],YES,p3,NO,YES,node_v,RP->state[2],
			      &RP->ang[1],&th12),
		    compute_mach_angles)
 
		Check_return(
		    s_polar_3(RP->state[0],YES,p3,YES,YES,node_v,RP->state[3],
			      &RP->ang[3],&th03),
		    compute_mach_angles)

		m01 = mass_flux(pressure(RP->state[1]),RP->state[0]);
		m12 = mass_flux(p3,RP->state[1]);
		m03 = mass_flux(p3,RP->state[0]);

		q1 = sqrt(q0*q0 - 
			  (pressure(RP->state[1]) - pressure(RP->state[0]))*
			  (1./Dens(RP->state[1]) + 1./Dens(RP->state[0])));

		omega0 = asin(m01/(Dens(RP->state[0])*q0));
		omega1 = asin(m12/(Dens(RP->state[1])*q1));
		omegas = asin(m03/(Dens(RP->state[0])*q0));

		beta1 = omega1 + th12;		/* th12 < 0 */
		betas = omegas - th03;
		omega1p = omega1 - th01;
		Omega1 = PI - omega1 - (omega0 - th01);

		thetaw = (-(b+1) +
			  sqrt((b+1)*(b+1) - 4*a*(c-PI/2.+omega0)))/(2.*a);
		chi = a*thetaw*thetaw + b*thetaw + c;

		print_line_of_floats(7,omega0,thetaw,chi,omega1p,
				     Omega1,beta1,betas);
	}

	clean_up(0);
}		/*end compute_mach_angles*/


/*
*			compute_anom():
*
*	This function is used to find the angle at which the solution for
*	the three shock configuration ceases to exist.  This is the transition
*	point to anomolous reflection, and corresponds to the point where
*	the pressure behind the incident has reached the sonic point on
*	the incident polar.
*/

LOCAL boolean compute_anom(
	double		anom_ang,
	double		*diff,
	POINTER		parameters)
{
	RP_DATA		*RP = ((MACH_ANG_PARAMS *)parameters)->RP;
	double		inc_speed = ((MACH_ANG_PARAMS *)parameters)->inc_speed;
	double		node_v[MAXD], rv0[MAXD];
	double		M0sq, Mnsq;
	double		p1;

	node_v[0] = inc_speed;
	node_v[1] = inc_speed / tan(anom_ang);

	M0sq = mach_number_squared(RP->state[0],node_v,rv0);
	Mnsq = inc_speed*inc_speed/sound_speed_squared(RP->state[0]);
	p1 = max_behind_shock_pr(Mnsq,RP->state[0]);

	*diff = p1 - pressure_at_sonic_point(M0sq,RP->state[0]);
	return FUNCTION_SUCCEEDED;
}		/*end compute_anom*/


/*
*			compute_detachment():
*
*	This function is simply a shell for the function s_polar_2() for
*	use in the root finder.  If s_polar2() succeeds, then a regular
*	reflection is possible.  The point where it begins to fail is 
*	the point where the reflected polar is tangent to the verticle
*	axis in the shock polar diagram.  This is the so called detachment
*	angle.  The iteration proceeds on the angle between the incident
*	and the incoming streamline.  The bisection proceeds by creating a
*	step function.  If s_polar_2() succeeds, we return +1, if not -1.  
*	Thus the point of jump is the detachment angle.  The pressure
*	at the detachment angle is stored for use in the following function.
*/

LOCAL boolean compute_detachment(
	double		detach_ang,
	double		*diff,
	POINTER		parameters)
{
	RP_DATA		*RP = ((MACH_ANG_PARAMS *)parameters)->RP;
	double		p1, turn_angle;
	double		node_v[MAXD], Mnsq;
	double		inc_speed = ((MACH_ANG_PARAMS *)parameters)->inc_speed;
	int		is_back_weak = YES;

	node_v[0] = inc_speed;
	node_v[1] = inc_speed / tan(detach_ang);

	Mnsq = inc_speed*inc_speed/sound_speed_squared(RP->state[0]);
	p1 = max_behind_shock_pr(Mnsq,RP->state[0]);

	if (!s_polar_3(RP->state[0],YES,p1,YES,NO,node_v,
			       RP->state[1],&RP->ang[0],&turn_angle))
	{
		return FUNCTION_FAILED;
	}

	if (!s_polar_2(RP->state[1],YES,is_back_weak,
			turn_angle,node_v,RP->state[2],&RP->ang[2]))
		*diff = 1.0;
	else
		*diff = -1.0;
	return FUNCTION_SUCCEEDED;
}		/*end compute_detachment*/

/*
*			compute_sonic():
*
*	This function is very similar to the preceeding.  Here we
*	want to compute the sonic angle, which is defined as the angle
*	at which the reflected polar crosses the verticle axis at the
*	sonic point.  Here, we do not create a step function, but just
*	use the difference between the computed pressure in state 2
*	and the pressure at the sonic point on the reflected polar.  The
*	variable is again the angle between the incoming streamline and
*	the incident.  This angle needs to be increased (decreased) as
*	(p2 - psonic) is greater (less) than zero.  In the case of no
*	regular reflection solution (cannot compute p2), the pressure at
*	the detachment angle is substituted for p2 in the difference.
*/

LOCAL boolean compute_sonic(
	double		sonic_ang,
	double		*diff,
	POINTER		parameters)
{
	RP_DATA		*RP = ((MACH_ANG_PARAMS *)parameters)->RP;
	int		is_given_ahead = YES, is_back_weak = YES;
	double		turn_angle;
	double		p1;
	double		rv1[MAXD], M1sq, Mnsq;
	double		node_v[MAXD];
	double		inc_speed = ((MACH_ANG_PARAMS *)parameters)->inc_speed;

	node_v[0] = inc_speed;
	node_v[1] = inc_speed / tan(sonic_ang);

	Mnsq = inc_speed*inc_speed/sound_speed_squared(RP->state[0]);
	p1 = max_behind_shock_pr(Mnsq,RP->state[0]);

	if (!s_polar_3(RP->state[0],YES,p1,YES,NO,node_v,
			       RP->state[1],&RP->ang[0],&turn_angle))
	{
		return FUNCTION_FAILED;
	}

	M1sq = mach_number_squared(RP->state[1],node_v,rv1);

	if (!s_polar_2(RP->state[1],is_given_ahead,is_back_weak,
			turn_angle,node_v,RP->state[2],&RP->ang[2]))
	{
		*diff = ((MACH_ANG_PARAMS *)parameters)->p_detach - 
			pressure_at_sonic_point(M1sq,RP->state[1]);
		return FUNCTION_SUCCEEDED;
        }

	*diff = pressure(RP->state[2]) -
		pressure_at_sonic_point(M1sq,RP->state[1]);
	return FUNCTION_SUCCEEDED;
}		/*end compute_sonic*/

/*
*			compute_mech_eq():
*
*	This function computes the mechanical equilibrium angle, and is
*	very similar to compute_sonic().  The trial angle should be 
*	increased (decreased) as (p2 - max pressure on incident polar) is
*	greater (less) than zero.   In the case of no
*	regular reflection solution (cannot compute p2), the pressure at
*	the detachment angle is substituted for p2 in the difference.
*/

LOCAL boolean compute_mech_eq(
	double		mech_eq_ang,
	double		*diff,
	POINTER		parameters)
{
	RP_DATA		*RP = ((MACH_ANG_PARAMS *)parameters)->RP;
	double		turn_angle;
	double		p1, M0sq, Mnsq;
	double		node_v[MAXD], rv0[MAXD];
	double		inc_speed = ((MACH_ANG_PARAMS *)parameters)->inc_speed;
	int		is_given_ahead = YES, is_back_weak = YES;

	node_v[0] = inc_speed;
	node_v[1] = inc_speed / tan(mech_eq_ang);

	M0sq = mach_number_squared(RP->state[0],node_v,rv0);
	Mnsq = inc_speed*inc_speed/sound_speed_squared(RP->state[0]);
	p1 = max_behind_shock_pr(Mnsq,RP->state[0]);

	if (!s_polar_3(RP->state[0],YES,p1,YES,NO,node_v,
			       RP->state[1],&RP->ang[0],&turn_angle))
	{
		return FUNCTION_FAILED;
	}

	if (!s_polar_2(RP->state[1],is_given_ahead,is_back_weak,
			turn_angle,node_v,RP->state[2],&RP->ang[2]))
	{
		*diff = ((MACH_ANG_PARAMS *)parameters)->p_detach -
			max_behind_shock_pr(M0sq,RP->state[0]);
		return FUNCTION_SUCCEEDED;
        }
	*diff = pressure(RP->state[2]) - 
		max_behind_shock_pr(M0sq,RP->state[0]);
	return FUNCTION_SUCCEEDED;
}		/*end compute_mech_eq*/


/*
*			compute_vnr():
*
*	This function is used by the root finder to compute the angle
*	of incidence corresponding to the onset of von Neumann reflection.
*	The trial should be increased (decreased) as p3 (the pressure
*	at the interesection of the incident and reflected polars) is 
*	greater (less) than the max pressure on the reflected polar.  The
*	case of angles that are too large, for which there is no solution
*	for a Mach configuration, is handled by setting the upper bound
*	computed in the previous function.
*/

LOCAL boolean compute_vnr(
	double		vnr_ang,
	double		*diff,
	POINTER		parameters)
{
	RP_DATA		*RP = ((MACH_ANG_PARAMS *)parameters)->RP;
	double		inc_speed = ((MACH_ANG_PARAMS *)parameters)->inc_speed;
	double		turn_angle;
	double		p1, p2, Mnsq, M1sq;
	double		node_v[MAXD], rv1[MAXD];
	int		is_given_ahead = YES;

	node_v[0] = inc_speed;
	node_v[1] = inc_speed / tan(vnr_ang);

	Mnsq = inc_speed*inc_speed/sound_speed_squared(RP->state[0]);

	p1 = max_behind_shock_pr(Mnsq,RP->state[0]);
	if (!s_polar_3(RP->state[0],YES,p1,YES,NO,node_v,
			       RP->state[1],&RP->ang[0],&turn_angle))
	{
		return FUNCTION_FAILED;
	}

	M1sq = mach_number_squared(RP->state[1],node_v,rv1);
	p2 = max_behind_shock_pr(M1sq,RP->state[1]);

	if (!s_polar_2(RP->state[0],is_given_ahead,NO,
			       turn_angle,node_v,RP->state[3],&RP->ang[3]))
	{
		return FUNCTION_FAILED;
	}

	*diff = pressure(RP->state[3]) - p2;
	return FUNCTION_SUCCEEDED;
}		/*end compute_vnr*/

LOCAL	void	set_bubble_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	comp_type->type = BUBBLE;
	if (comp_type->extra == NULL)
	    comp_type->extra=(POINTER)allocate_bubble(NO,front,NULL,GAS_STATE);
	comp_type->_get_state = get_state_bubble;
	comp_type->free_comp_type_extra = free_bubble_comp_type;
}		/*end set_bubble_comp_type*/

LOCAL	void	free_bubble_comp_type(
	COMP_TYPE	*comp_type)
{
	Bubble	*bubble;

	if (comp_type->type != BUBBLE)
		return;

	bubble = (Bubble*)comp_type->extra;
	if (bubble != NULL)
		free_bubble(bubble,YES);
	comp_type->extra = NULL;
}		/*end free_bubble_comp_type*/

/*ARGSUSED*/
LOCAL	void	get_state_bubble(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	debug_print("init_states","Entered get_state_bubble()\n");
	find_bubble_state(s,coords,ct->comp,intfc,stype);
}		/*end get_state_bubble*/

#endif /* defined(TWOD) */
