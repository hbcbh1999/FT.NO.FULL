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
*				gbifur.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Performs bifurcations of two dimensional interfaces
*
*/


#include <gdecs/gdecs.h>

#if defined(TWOD)
	/* LOCAL Function Declarations */
LOCAL	double	inc_angle(double*,BOND*,ORIENTATION,int);
LOCAL	void	create_mach_contact(O_CURVE*,O_CURVE*,Bubble*,Front*);
LOCAL	void	create_mach_stem(O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,BOND*,
				 Bubble*,Front*,Wave*,double,NODE_FLAG);
LOCAL	void	init_bow(COMPONENT,O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,BOND*,
			 Locstate*,NODE*,Bubble*,Front*,int,boolean);
LOCAL	void	insert_normal_contact_at_wall(Front*,O_CURVE*,O_CURVE*,
					      O_CURVE*,O_CURVE*);
LOCAL	void	normal_to_regular_reflection(double,double*,BOND*,O_CURVE*,
					     O_CURVE*,O_CURVE*,COMPONENT,
					     Bubble*,Front*);
LOCAL	void	normal_to_mach_reflection(double,BOND*,O_CURVE*,O_CURVE*,
					  O_CURVE*,O_CURVE*,O_CURVE*,
					  COMPONENT,Bubble*,Front*,Wave*,
					  double,NODE_FLAG);
LOCAL	int	untrack_B_reflect_node(Front*,Wave*,O_CURVE*,O_CURVE*,
				       O_CURVE*,O_CURVE*,ANGLE_DIRECTION,
				       double,double*,RPROBLEM*,NODE_FLAG);
#endif /* defined(TWOD) */

#if defined(ONED) || defined(TWOD)

/*
*			track_scattered_wave():
*
*	This function decides which scattered waves to track when a
*	bifurcation or tangle occurs.  Input includes the node type,
*	wave class (shock, contact or rarefaction edge), and the curve
*	status to be considered.  An ahead state (s0) and behind state (s1)
*	are also included.  A wave strength is then computed and compared
*	to the cutoff stored in scatter_wave_tolerance(fr,ntype,wclass,nstate).
*	For shocks, the strength is the pressure ratio (behind/ahead) and for
*	contacts the density ratio.  For rarefactions, the strength is again
*	the pressure ratio but the behind state should be that behind the
*	entire rarefaction, not just the edge under consideration, and this
*	strength is ahead/behind since pressure decreases across a rarefaction.
*
*	Note: this function is not really appropriate for use in checking
*	every time step whether to track a curve or not.  That would require
*	a point by point check along the curve and some sort of criterion to
*	decide if some aggregate wave strength should still be considered
*	large enough.
*/

EXPORT	boolean track_scattered_wave(
	int		n_type,
	int		w_class,
	int		status,
	Locstate	s0,
	Locstate	s1,
	Front		*fr)
{
	SCAT_WV_TOL *swt;
	double	    strength, tol;
	const double eps = 10.0*MACH_EPS;/*TOLERANCE*/

	if (!Same_params(s0,s1))
	    return YES;

	swt = scattered_wave_tolerance(fr,n_type,w_class,status);

	if (w_class == CONTACT)
	{
	    if (Different_params(s0,s1))
	        return YES;
	}
	switch (swt->tol_type)
	{
	case NEVER_TRACK:
	    return NO;
	case ALWAYS_TRACK:
	    switch (w_class)
	    {
	    case SHOCK_WAVE:
	        strength = pressure(s1)/pressure(s0) - 1.0;
		return (strength > eps) ? YES : NO;
	    case CONTACT_WAVE:
	        strength = fabs(Dens(s1) - Dens(s0))/(Dens(s0) + Dens(s1));
		return (strength > eps) ? YES : NO;
	    case RAREF_LEADING_EDGE:
	    case RAREF_TRAILING_EDGE:
	        strength = pressure(s0)/pressure(s1) - 1.0;
		return (strength > eps) ? YES : NO;
	    default:
	        return YES;
	    }
	case PRESSURE_RATIO:
	    if (w_class == SHOCK_WAVE)
	        strength = pressure(s1)/pressure(s0) - 1.0;
	    else
	        strength = pressure(s0)/pressure(s1) - 1.0;
	    break;
	case ATWOOD_NUMBER:
	    switch (w_class)
	    {
	    case SHOCK_WAVE:
	        strength = (Dens(s1) - Dens(s0))/(Dens(s0) + Dens(s1));
	        break;
	    case CONTACT_WAVE:
	        strength = fabs(Dens(s1) - Dens(s0))/(Dens(s0) + Dens(s1));
	        break;
	    case RAREF_LEADING_EDGE:
	    case RAREF_TRAILING_EDGE:
	        strength = (Dens(s0) - Dens(s1))/(Dens(s0) + Dens(s1));
	        break;
	    }
	    break;
	case MACH_NUMBER_AHEAD:
	    if (w_class != SHOCK_WAVE)
		return NO;
	    strength = mass_flux(pressure(s1),s0)/acoustic_impedance(s0) - 1.0;
	}
	tol = swt->wave_tol;
	return (strength > tol) ? YES : NO;
}		/*end track_scattered_wave*/

#endif /* defined(ONED) || defined(TWOD) */

#if defined(TWOD)

#define ERROR_ANGLE (.4166667*PI)	/* TOLERANCE -- 75 degrees */

/*
*			g_B_node_bifurcation():
*
*	A B_node can bifurcate into a regular reflection, Mach stem or
*	attached bow shock. This routine determines if such a bifurcation
*	should take place and depending on the variable flag, carries out
*	the bifurcation.
*	It is assumed that the boundary is a NEUMANN_BOUNDARY, that the
*	computational frame is such that the boundary is at rest and
*	that the physical curve is a shock. Then the ahead side of the
*	boundary is the side towards which the physical curve is moving,
*	while the forward_facing side is the lower pressure side. When
*	these two sides agree, we set forward_facing_is_ahead = YES.
*	Otherwise forward_facing_is_ahead = NO.  The incident angle is defined
*	as the angle between the ahead side and the physical curve,
*	normalized to be in the interval [0,PI).  If forward_facing_is_ahead,
*       then incident angle < PI/2 means a transition to a regular reflection
*       or (simple, complex, double, ...) mach stem, while incident angle >
*       PI/2 means the physical curve simply moves past an
*       expansion corner which may also be a node.  Otherwise (forward_facing_
*       is_ahead == NO), means a bow shock, which may be attached (leading
*       to a bifurcation and change in topology) or detached (and no
*       bifurcation).
*/

EXPORT int g_B_node_bifurcation(
	Front		*fr,
	POINTER		p2wave,
	O_CURVE		*oldcphys,
	O_CURVE		*newcphys,
	O_CURVE		*oldca,
	O_CURVE		*newca,
	O_CURVE		*oldcaprop,
	O_CURVE		*newcaprop,
	O_CURVE		*oldcb,
	O_CURVE		*newcb,
	O_CURVE		*oldcbprop,
	O_CURVE		*newcbprop,
	POINT		*oldposn,
	Locstate	sl,
	Locstate	sr,
	ANGLE_DIRECTION	i_to_prop_dir,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	BOND		*bcorner;
	COMPONENT	compbow_ext;
	CURVE		*tmpc;
	NODE		*interact_node[5];
	POINT		*newposn = Node_of_o_curve(newcphys)->posn;
	POINT		*corner;
	RP_DATA		*RP;
	Wave		*wave = (Wave*)p2wave;
	double		initial_time_elapsed;
	double		node_v[MAXD];
	double		inc_t[MAXD];
	double		inc_ang;		/* ahead wall to inc */
	double		tol_angle;
	int		forward_facing_is_ahead = 1;
	int		i, dim = fr->interf->dim;
	ORIENTATION	tmpc_orient;
	Bubble		*bubble;

	debug_print("gBnode","Entered g_B_node_bifurcation\n");

	if (is_scalar_wave(wave_type(newcphys->curve)) &&
	    (wave_type(newcaprop->curve) == NEUMANN_BOUNDARY))
	{
	    insert_normal_contact_at_wall(fr,oldcphys,newcphys,oldcaprop,
					  oldcb);
	    debug_print("gBnode","Left g_B_node_bifurcation\n");
	    return GOOD_NODE;
	}
	if ((wave_type(newcphys->curve) < FIRST_VECTOR_PHYSICS_WAVE_TYPE) ||
		(is_rarefaction_wave(wave_type(newcphys->curve))) ||
		(wave_type(newcaprop->curve) == DIRICHLET_BOUNDARY) ||
		is_subdomain_boundary(Hyper_surf(newcaprop->curve)))
	{
	    debug_print("gBnode","Left g_B_node_bifurcation\n");
	    return GOOD_NODE;
	}
	if (debugging("gBnode"))
	{
	    (void) printf("oldposn = <%g, %g>\n",Coords(oldposn)[0],
	    	      Coords(oldposn)[1]);
	    (void) printf("oldcphys - ");   print_o_curve(oldcphys);
	    (void) printf("newcphys - ");   print_o_curve(newcphys);
	    (void) printf("oldca - ");  print_o_curve(oldca);
	    (void) printf("newca - ");  print_o_curve(newca);
	    (void) printf("oldcaprop - ");  print_o_curve(oldcaprop);
	    (void) printf("newcaprop - ");  print_o_curve(newcaprop);
	    (void) printf("oldcb - "); print_o_curve(oldcb);
	    (void) printf("newcb - "); print_o_curve(newcb);
	    (void) printf("oldcbprop - ");  print_o_curve(oldcbprop);
	    (void) printf("newcbprop - ");  print_o_curve(newcbprop);
	    verbose_print_state("sl",sl);
	    verbose_print_state("sr",sr);
	}

	if (i_to_prop_dir == COUNTER_CLOCK)
	    forward_facing_is_ahead *= -1;
	if (	(wave_type(newcphys->curve) == BACKWARD_SHOCK_WAVE &&
	         newcphys->orient == POSITIVE_ORIENTATION)
		 ||
		(wave_type(newcphys->curve) == FORWARD_SHOCK_WAVE &&
		 newcphys->orient == NEGATIVE_ORIENTATION)) 
	{
	    forward_facing_is_ahead *= -1;
	}
	forward_facing_is_ahead =(forward_facing_is_ahead == 1) ? YES : NO;

			/* No Transition? */

	find_tangent_to_propagated_curve(newposn,
		                         Bond_at_node_of_o_curve(newcphys),
					 oldcphys,newcphys,inc_t,fr,p2wave,dt);
	inc_ang = inc_angle(inc_t,Bond_at_node_of_o_curve(newcaprop),
		            newcaprop->orient,dim);
	if (continue_past_fixed_node(flag) == YES)
	    tol_angle = 0.5*PI;
	else
	    tol_angle = (forward_facing_is_ahead) ?
			0.5*PI - ERROR_ANGLE : 0.5*PI + ERROR_ANGLE;

	if (debugging("gBnode"))
	{
	    (void) printf("Bonds: phys newcaprop newcb:\n");
	    print_bond(Bond_at_node_of_o_curve(newcphys));
	    print_bond(Bond_at_node_of_o_curve(newcaprop));
	    print_bond(Bond_at_node_of_o_curve(newcb));
	    print_angle_direction("i_to_prop_dir = ",i_to_prop_dir,"\n");
	    (void) printf("forward_facing_is_ahead = %s\n",
	    	      (forward_facing_is_ahead) ? "YES" : "NO");
	}

	if (forward_facing_is_ahead)
	{
	    if (inc_ang >= tol_angle)
	    {
	    	if (debugging("gBnode"))
	    	{
	    	    (void) printf("No transition, ");
	    	    (void) printf("forward_facing_is_ahead = YES\n");
	    	    print_angle("inc_ang =",inc_ang," ");
	    	    print_angle(">= tol_angle =",tol_angle,"\n");
	    	}
	    	debug_print("gBnode","Left g_B_node_bifurcation\n");
	    	return GOOD_NODE;
	    }
	}
	else if (inc_ang <= tol_angle)
	{
	    if (debugging("gBnode"))
	    {
		(void) printf("No transition, forward_facing_is_ahead = NO\n");
		print_angle("inc_ang =",inc_ang," ");
		print_angle("<= tol_angle =",tol_angle,"\n");
	    }
	    debug_print("gBnode","Left g_B_node_bifurcation\n");
	    return GOOD_NODE;
	}

	if (is_subdomain_boundary(Hyper_surf(newcb->curve)))
	{
	    /* Cannot install bifurcation across virtual boundary */

	    (void) printf("Warning in g_B_node_bifurcation(), ");
	    (void) printf("cannot install bifurcation across subdomain\n");
	    debug_print("gBnode","Left g_B_node_bifurcation\n");
	    return MODIFY_TIME_STEP_NODE;
	}

			/* Initialize RP */

	RP = Rp_data(Node_of_o_curve(newcphys));
	if (wave_type(newcphys->curve) == FORWARD_SHOCK_WAVE)
	{
	    ft_assign(RP->state[0],sr,fr->sizest);
	    ft_assign(RP->state[1],sl,fr->sizest);
	}
	else
	{
	    ft_assign(RP->state[0],sl,fr->sizest);
	    ft_assign(RP->state[1],sr,fr->sizest);
	}
	RP->ang_dir = (forward_facing_is_ahead) ?
		Opposite_ang_dir(i_to_prop_dir) : i_to_prop_dir;
	RP->ang[1] = angle(inc_t[0],inc_t[1]);

			/* Find corner data */

	set_corner_for_bifurcation(oldcphys,newcphys,oldcaprop,newcaprop,
		                   oldcb,newcb,oldposn,newposn,
				   fr,wave,&corner,&bcorner,
				   forward_facing_is_ahead,node_v,
				   RP,dt,dt_frac);

			/* Initialize bubble */

	bubble = allocate_bubble(YES,fr,RP,NO_COMP);
	bubble->is_node_at_corner =
	    (Node_of_o_curve(newcaprop) == Node_of_o_curve(newcb)) ?
	    NO : YES;
	for (i = 0; i < dim; i++)
	    bubble->inc_t[i] = inc_t[i];
	find_tangent_to_curve(newposn,Bond_at_node_of_o_curve(newcaprop),
		              newcaprop->curve,newcaprop->orient,
			      bubble->aw_t,fr);
	find_tangent_to_curve(Node_of_o_curve(newcb)->posn,
		              Bond_at_node_of_o_curve(newcb),
			      newcb->curve,newcb->orient,
			      bubble->bw_t,fr);

	if (debugging("gBnode"))
	{
	    (void) printf("Bond: bcorner\n");
	    print_bond(bcorner);
	    (void) printf("is_node_at_corner = %s node_v = %g %g dt_frac %g\n",
			  (bubble->is_node_at_corner)?"YES":"NO",
			  node_v[0],node_v[1],*dt_frac);
	    (void) printf("behind_w_ang = %g d corner = <%g, %g>\n",
			  degrees(angle(bubble->bw_t[0],bubble->bw_t[1])),
			  Coords(corner)[0],Coords(corner)[1]);
	}

			/* Is bowshock detached? */

	if (!forward_facing_is_ahead)
	{
	    double	attached_shock_ang;

	    if (!bubble->is_node_at_corner)
	    {
	    	debug_print("gBnode","Left g_B_node_bifurcation\n");
	    	return GOOD_NODE;
	    }
	    if (!is_bow_shock_attached(RP->state[0],
	    	                          angle(bubble->bw_t[0],
						bubble->bw_t[1]),
			                  RP->state[1],&attached_shock_ang))
	    {

	    	/* TODO: move node back to detached position */

	    	screen("ERROR in g_B_node_bifurcation(), "
		       "code to re-detach shock needed\n");
	    	clean_up(ERROR);
	    	debug_print("gBnode","Left g_B_node_bifurcation\n");
	    	return GOOD_NODE;
	    }
	}
	    	/* Initialize rp */

	if (continue_past_bifurcation(flag) != YES)
	{
	    debug_print("gBnode","Initializing rp\n");
	    interact_node[0] = Node_of_o_curve(newcphys);
	    interact_node[1] = Node_of_o_curve(oldcb);
	    interact_node[2] = NULL;

	    if (bubble->is_node_at_corner)
	    {
	        interact_node[2] = Node_of_o_curve(newcb);
	        Check_return(next_boundary(oldcb->curve,oldcb->orient,
				           &tmpc,&tmpc_orient),
			     g_B_node_bifurcation)
		interact_node[3] = Node_of(tmpc,Opposite_orient(tmpc_orient));
		interact_node[4] = NULL;
	    }
	    /*TMP*/
	    printf("Test augment_rproblem_list() position 1\n");
	    augment_rproblem_list(rp,interact_node,dt,*dt_frac,
			          oldcb->curve->interface,
			          newcb->curve->interface,fr,p2wave);
	    debug_print("gBnode","Completed initialization of rp\n");
	    debug_print("gBnode","Left g_B_node_bifurcation\n");
	    return BIFURCATION_NODE;
	}

	if (!forward_facing_is_ahead)
	{
	    /* Transition to an attached bow shock */

	    node_type(Node_of_o_curve(newcphys)) = ATTACHED_B_NODE;
	    if (wave_type(newcphys->curve) == FORWARD_SHOCK_WAVE)
	    	ft_assign(sl,RP->state[1],fr->sizest);
	    else
	    	ft_assign(sr,RP->state[1],fr->sizest);
	    /* TODO: Delete extra (null) bdry curve;	*
	     *   move node to corner;			*
	     *  assign states				*/
	    screen("ERROR in g_B_node_bifurcation(), "
		   "init attached bow shock code needed\n");
	    clean_up(ERROR);
	    debug_print("gBnode","Left g_B_node_bifurcation\n");
	    return GOOD_NODE;
	}

	bubble->comp_bow = new_component(NEW_COMP);
	compbow_ext = (wave_type(newcphys->curve) == FORWARD_SHOCK_WAVE) ?
	    negative_component(newcphys->curve) :
	    positive_component(newcphys->curve);
	for (i = 0; i < dim; i++)
	    bubble->corner_posn[i] = Coords(corner)[i];
	for (i = 0; i < dim; i++) bubble->cor_v[i] = 0.0;
	initial_time_elapsed = dt * (1.0 - *dt_frac);

	debug_print("gBnode","Calling is_regular_reflection\n");
	if (is_regular_reflection(node_v,fr,RP))
	{
	    if (!track_scattered_wave(B_REFLECT_NODE,SHOCK_WAVE,
					 REFLECTED,RP->state[1],
					 RP->state[2],fr))
	    {
	    	debug_print("gBnode","Left_g_B_node_bifurcation()\n");
	    	return GOOD_NODE;
	    }

	    normal_to_regular_reflection(initial_time_elapsed,node_v,
			                 bcorner,newcphys,newcb,newcaprop,
			                 compbow_ext,bubble,fr);
	}
	else
	{			
	    /* Cannot check here for tracking Mach stem because state
	     * behind Mach stem hasn't been computed yet. */

	    normal_to_mach_reflection(initial_time_elapsed,bcorner,oldcphys,
				      newcphys,newcb,newcaprop,newcbprop,
			              compbow_ext,bubble,fr,wave,dt,flag);
	}

	free_bubble(bubble,NO);
	debug_print("gBnode","Left g_B_node_bifurcation\n");
	return GOOD_NODE;
}		/*end g_B_node_bifurcation*/


/*
*			normal_to_regular_reflection():
*
*	This function is responsible for the installatation and change
*	of topology from a normal shock (B_NODE) to a regular reflection
*	configuration (B_REFLECT_NODE).
*/

LOCAL void normal_to_regular_reflection(
	double		initial_time_elapsed,
	double		*node_v,
	BOND		*bcorner,
	O_CURVE		*newcphys,
	O_CURVE		*newcb,
	O_CURVE		*newcaprop,
	COMPONENT	compbow_ext,
	Bubble		*bubble,
	Front		*fr)
{
	NODE		*rnode;
	O_CURVE		Bow_shock;
	O_CURVE		Bow_wall;
	O_CURVE		Refl_wall;		/* should == newcbprop */
	Locstate	behind = bubble->RP->state[1];
	Locstate	refl_st = bubble->RP->state[2];
	Locstate	exterior_bow_state;
	ANGLE_DIRECTION	i_to_a_dir;
	int		i, dim = fr->interf->dim;

	debug_print("gBnode","Creating normal_to_regular_reflection bifurcation\n");

	i_to_a_dir = Opposite_ang_dir(bubble->RP->ang_dir);
	bubble->RP->ang[0] = angle(bubble->aw_t[0],bubble->aw_t[1]);
	bubble->RP->ang[3] = angle(bubble->bw_t[0],bubble->bw_t[1]);

	init_bubble(node_v,TIME_ELAPSED,initial_time_elapsed,bubble,fr);

	if (bubble->is_attached)
	{
	    screen("ERROR in normal_to_regular_reflection(), "
	           "bow shock is attached (reg refl). CODE NEEDED\n");
	    clean_up(ERROR);
	}
	if (debugging("gBnode"))
	{
	    (void) printf("real reflected wall length = %g\n",
	    	          bubble->refl_length);
	    (void) printf("real bow wall length = %g\n",bubble->bow_length);
	}

		/* Create bow shock */

	Bow_shock.orient = POSITIVE_ORIENTATION;
	rnode = Node_of_o_curve(newcphys);
	node_type(rnode) = B_REFLECT_NODE;
	if (is_bdry(newcaprop->curve))
	    set_is_bdry(rnode);
	else
	    set_not_bdry(rnode);

	init_bow(compbow_ext,newcb,&Bow_shock,&Bow_wall,&Refl_wall,
		 bcorner,&exterior_bow_state,rnode,bubble,fr,NO,
		 track_scattered_wave(B_REFLECT_NODE,SHOCK_WAVE,REFLECTED,
				      bubble->RP->state[1],
				      bubble->RP->state[2],fr));
	if (Bow_shock.curve == NULL)
	{
	    screen("ERROR in normal_to_regular_reflection(), "
	           "untracked refl shock. CODE NEEDED\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < dim; i++)
	    bubble->refl_posn[i] = Coords(rnode->posn)[i];

	    	/* Initialize bubble interface states */

	if (curve_ang_oriented_l_to_r(i_to_a_dir,Bow_shock.orient))
	{
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,NEGATIVE_SIDE,
					refl_st,bubble->bow,fr->sizest);
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,POSITIVE_SIDE,
					behind,exterior_bow_state,fr->sizest);
	}
	else
	{
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,POSITIVE_SIDE,
					refl_st,bubble->bow,fr->sizest);
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,NEGATIVE_SIDE,
					behind,exterior_bow_state,fr->sizest);
	}

	if (bubble->is_node_at_corner)
	{
	    if (curve_ang_oriented_l_to_r(i_to_a_dir,Bow_wall.orient))
	    {
	    	set_states_by_interpolation(Bow_wall.curve,NULL,NULL,
					    POSITIVE_SIDE,bubble->refl_corner,
					    bubble->bow,fr->sizest);
		set_states_by_interpolation(Refl_wall.curve,NULL,NULL,
					    POSITIVE_SIDE,refl_st,
					    bubble->refl_corner,fr->sizest);
	    }
	    else
	    {
	    	set_states_by_interpolation(Bow_wall.curve,NULL,NULL,
					    NEGATIVE_SIDE,bubble->bow,
					    bubble->refl_corner,fr->sizest);
	    	set_states_by_interpolation(Refl_wall.curve,NULL,NULL,
					    NEGATIVE_SIDE,bubble->refl_corner,
					    refl_st,fr->sizest);
	    }
	}
	else
	{
	    if (curve_ang_oriented_l_to_r(i_to_a_dir,Bow_wall.orient))
	    {
	    	set_states_by_interpolation(Bow_wall.curve,NULL,NULL,
					    POSITIVE_SIDE,refl_st,bubble->bow,
					    fr->sizest);
	    }
	    else
	    {
	    	set_states_by_interpolation(Bow_wall.curve,NULL,NULL,
					    NEGATIVE_SIDE,bubble->bow,
					    refl_st,fr->sizest);
	    }
	}
}		/*end normal_to_regular_reflection*/


/*
*			normal_to_mach_reflection():
*
*	This function is responsible for the installatation and change
*	of topology from a normal shock (B_NODE) to a Mach reflection
*	configuration (MACH_NODE).
*/

LOCAL void normal_to_mach_reflection(
	double		initial_time_elapsed,
	BOND		*bcorner,
	O_CURVE		*oldcphys,
	O_CURVE		*newcphys,
	O_CURVE		*newcb,
	O_CURVE		*newcaprop,
	O_CURVE		*newcbprop,
	COMPONENT	compbow_ext,
	Bubble		*bubble,
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE_FLAG	flag)
{
	Locstate	behind = bubble->RP->state[1];
	Locstate	refl_st = bubble->RP->state[2];
	Locstate	exterior_bow_state;
	O_CURVE		Bow_shock;
	O_CURVE		Bow_wall;
	O_CURVE		Refl_wall;		/* should == newcbprop */
	ANGLE_DIRECTION	i_to_a_dir;

	debug_print("gBnode","Creating Mach reflection bifurcation\n");

	    /* Angles labeled differently for mach node */
	bubble->RP->ang[0] = bubble->RP->ang[1];
	i_to_a_dir = Opposite_ang_dir(bubble->RP->ang_dir);

	if (debugging("gBnode"))
	{
	    (void) printf("Calling init_mach_bubble() ");
	    (void) printf("in g_B_node_bifurcation()\n");
	}

	if (!init_mach_bubble(TIME_ELAPSED,initial_time_elapsed,
				 bubble,NORMAL_TO_MACH_REFLECTION,fr)
    		 ||
	    !track_scattered_wave(MACH_NODE,SHOCK_WAVE,MACH_STEM,
				     bubble->RP->state[0],
				     bubble->RP->state[3],fr))
	{
	    node_type(Node_of_o_curve(newcphys)) = NEUMANN_NODE;
	    return;
	}

	bubble->comp_mach = new_component(NEW_COMP);

	if (bubble->is_attached)
	{
	    screen("ERROR in normal_to_mach_reflection(), "
	           "bow shock is attached (mach refl). CODE NEEDED\n");
	    clean_up(ERROR);
	}

	create_mach_stem(oldcphys,newcphys,newcaprop,newcbprop,
			 Bond_at_node_of_o_curve(newcphys),bubble,
			 fr,wave,dt,flag);
	
		/* Create bow shock */

	if (track_scattered_wave(MACH_NODE,SHOCK_WAVE,REFLECTED,
				 bubble->RP->state[1],bubble->RP->state[2],fr))
	{
	    Bow_shock.orient = POSITIVE_ORIENTATION;
	    init_bow(compbow_ext,newcb,&Bow_shock,&Bow_wall,
	    	     &Refl_wall,bcorner,&exterior_bow_state,
		     Node_of_o_curve(newcphys),bubble,fr,YES,YES);
	}

			/* Initialize bubble interface states */

	if (curve_ang_oriented_l_to_r(i_to_a_dir,Bow_shock.orient))
	{
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,NEGATIVE_SIDE,
					refl_st,bubble->bow,fr->sizest);
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,
					POSITIVE_SIDE,behind,exterior_bow_state,
					fr->sizest);
	}
	else
	{
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,POSITIVE_SIDE,
					refl_st,bubble->bow,fr->sizest);
	    set_states_by_interpolation(Bow_shock.curve,NULL,NULL,NEGATIVE_SIDE,
					behind,exterior_bow_state,fr->sizest);
	}

	if (track_scattered_wave(MACH_NODE,CONTACT_WAVE,SLIP,
				 bubble->RP->state[2],bubble->RP->state[3],fr))
	    create_mach_contact(newcphys,newcbprop,bubble,fr);
}		/*end normal_to_mach_reflection*/


/*
*			g_reflect_node_bifurcation():
*
*	This routine gives the bifurcation of a regular reflection node to
*	a Mach triple point.
*/

/*ARGSUSED*/
EXPORT void g_reflect_node_bifurcation(
	Front		*fr,
	Wave		*wave,
	O_CURVE		*oldcinc,
	O_CURVE		*cinc,
	O_CURVE		*oldcrefl,
	O_CURVE		*crefl,
	O_CURVE		*oldca,
	O_CURVE		*cahead,
	O_CURVE		*oldcb,
	O_CURVE		*cbehind,
	O_CURVE		*oldcaprop,
	O_CURVE		*caprop,
	O_CURVE		*oldcbprop,
	O_CURVE		*cbprop,
	POINT		*oldposn,
	POINT		*posn,
	POINT		*corner,
	BOND		*crossbinc,
	BOND		*crossbprop,
	RP_DATA		*RP,
	RPROBLEM	*rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	Bubble		*bubble;
	Locstate	left_st, right_st;
	Locstate	behind = RP->state[1];
	Locstate	refl_bow = RP->state[2];
	double		initial_time_elapsed;
	double		i_n[MAXD];		/* inc normal to ahead */
	double		*inc_t;
	double		*aw_t;
	double		*refl_ang = &RP->ang[2];/* pos x axis to refl shock */
	boolean		track_slip;
	ANGLE_DIRECTION	i_to_a_dir;
	int		i, dim = fr->rect_grid->dim;

	debug_print("grBnode","Entered g_reflect_node_bifurcation\n");

	bubble = allocate_bubble(YES,fr,RP,NO_COMP);
	inc_t = bubble->inc_t;
	aw_t = bubble->aw_t;

	    /* Find tangents and normals */

	find_tangent_to_propagated_curve(posn,crossbinc,oldcinc,cinc,inc_t,
					 fr,(POINTER)wave,dt);
	find_tangent_to_propagated_curve(posn,crossbprop,oldcaprop,caprop,aw_t,
					 fr,(POINTER)wave,dt);
	if (RP->ang_dir == CLOCKWISE)
	{
	    i_to_a_dir = COUNTER_CLOCK;
	    i_n[0] = -inc_t[1]; i_n[1] = inc_t[0];
	}
	else
	{
	    i_to_a_dir = CLOCKWISE;
	    i_n[0] = inc_t[1]; i_n[1] = -inc_t[0];
	}
	if (aw_t[0]*i_n[0] + aw_t[1]*i_n[1] <= 0.)
	{
	    screen("ERROR in g_reflect_node_bifurcation(), "
	           "Propagated shock faces corner. New code needed.\n");
	    clean_up(ERROR);
	}
	if (debugging("grBnode"))
	{
	    print_general_vector("ahead wall tangent = ",aw_t,dim,"\n");
	    print_general_vector("incident normal = ",i_n,dim,"\n");
	}

	    	/* Initialize bubble parameters */

	if (crefl->curve != NULL)
	{
	    bubble->comp_bow = 
			(curve_ang_oriented_l_to_r(i_to_a_dir,crefl->orient)) ?
				negative_component(crefl->curve) :
				positive_component(crefl->curve);
	}
	else
	{
	    if (debugging("grBnode"))
	    	(void) printf("reflected wave untracked\n");

	    bubble->comp_bow = 
			(curve_ang_oriented_l_to_r(i_to_a_dir,cinc->orient)) ?
				negative_component(cinc->curve) :
				positive_component(cinc->curve);
	}

	initial_time_elapsed = dt * (1.0 - *dt_frac);
	bubble->is_node_at_corner =
	    (Node_of_o_curve(caprop) == Node_of_o_curve(cbprop)) ? NO : YES;
	for (i = 0; i < dim; i++)
	    bubble->cor_v[i] = 0.0;
	if (bubble->is_node_at_corner)
	{
	    find_tangent_to_curve(Node_of_o_curve(cbehind)->posn,
			          Bond_at_node_of_o_curve(cbehind),
			          cbehind->curve,POSITIVE_ORIENTATION,
				  bubble->bw_t,fr);
	}
	else
	{
	    for (i = 0; i < dim; i++)
		bubble->bw_t[i] = -bubble->aw_t[i];
	}

	    /* angles labeled differently for mach node */
	RP->ang[0] = RP->ang[1];

	for (i = 0; i < dim; i++)
	    bubble->corner_posn[i] = Coords(corner)[i];

	if (debugging("grBnode"))
	{
	    (void) printf("Calling init_mach_bubble\n");
	    (void) printf("initial_time_elapsed = %g dt = %g dt_frac = %g\n",
			  initial_time_elapsed,dt,*dt_frac);
	}
	if (!init_mach_bubble(TIME_ELAPSED,initial_time_elapsed,
				 bubble,REGULAR_TO_MACH_REFLECTION,fr))
	{
	    if (untrack_node(fr,B_REFLECT_NODE))
	    {
	    	if (!untrack_B_reflect_node(fr,wave,oldcinc,cinc,oldcrefl,
					       crefl,i_to_a_dir,dt,dt_frac,rp,
					       flag))
		{
		    screen("ERROR in g_reflect_node_bifurcation(), "
		           "unable to untrack node\n");
		    clean_up(ERROR);
		}
		debug_print("grBnode","Left g_reflect_node_bifurcation\n");
		return;
	    }
	    else
	    {
	    	screen("ERROR in g_reflect_node_bifurcation(), "
	    	       "init_mach_bubble() failed\n");
	    	clean_up(ERROR);
	    }
	}

	if (bubble->is_attached)
	{
	    screen("ERROR in g_reflect_node_bifurcation(), "
	           "bow shock is attached. CODE NEEDED\n");
	    clean_up(ERROR);
	}

	track_slip = track_scattered_wave(MACH_NODE,CONTACT_WAVE,SLIP,
					  bubble->RP->state[2],
					  bubble->RP->state[3],fr);

	bubble->comp_mach = (track_slip) ? new_component(NEW_COMP) :
					   bubble->comp_bow;

	    /* Propagate curves before change of topology */

	left_st = Left_state_at_node_of_o_curve(cbprop);
	right_st = Right_state_at_node_of_o_curve(cbprop);
	ft_assign(left_st,Left_state_at_node_of_o_curve(oldcbprop),fr->sizest);
	ft_assign(right_st,Right_state_at_node_of_o_curve(oldcbprop),fr->sizest);

	create_mach_stem(oldcinc,cinc,caprop,cbprop,
			 crossbinc,bubble,fr,wave,dt,flag);

	if (crefl->curve != NULL)
	{
	    /* Transfer reflected curve to Mach triple point */

	    change_node_of_curve(crefl->curve,crefl->orient,
	    		         Node_of_o_curve(cinc));
	    if (curve_ang_oriented_l_to_r(i_to_a_dir,crefl->orient))
	    {
	    	left_st = refl_bow;
	    	right_st = behind;
	    }
	    else
	    {
	    	left_st = behind;
	    	right_st = refl_bow;
	    }
	    (void) propagate_curve_near_node(Node_of_o_curve(oldcrefl),
					     Node_of_o_curve(crefl),
					     oldcrefl,crefl,
					     Bond_at_node_of_o_curve(crefl),
					     left_st,right_st,YES,*refl_ang,fr,
					     wave,dt,flag);
	    for (i = 0; i < dim; i++)
	    	bubble->bow_posn[i] =
		    Coords(Opp_node_of_o_curve(crefl)->posn)[i];


	    /* For bifurcation with sharp corner, we need to make sure
	     * the point that was added to enforce the refl_ang is not
	     * removed during redistribution, or else the refl curve
	     * could end up running through the ramp itself. */
	    if (bubble->is_node_at_corner)
	    	do_not_redistribute(crefl->curve) = YES;
	}

	if (debugging("grBnode"))
	{
	    if (crefl->curve != NULL)
	    {
	    	(void) printf("After transfer of crefl\n");
	    	(void) printf("Triple point node on cinc:\n");
	    	print_node(Node_of_o_curve(cinc));
	    	(void) printf("Triple point node on crefl:\n");
	    	print_node(Node_of_o_curve(crefl));
	    	(void) printf("New cinc:\n");
	    	print_o_curve(cinc);
	    	(void) printf("New crefl:\n");
	    	print_o_curve(crefl);
	    }
	    else
	    {
	    	(void) printf("Triple point node on cinc:\n");
	    	print_node(Node_of_o_curve(cinc));
	    	(void) printf("New cinc:\n");
	    	print_o_curve(cinc);
	    }
	}

	if (track_slip)
	    create_mach_contact(cinc,cbprop,bubble,fr);

	free_bubble(bubble,NO);
	debug_print("grBnode","Left g_reflect_node_bifurcation\n");
}		/*end g_reflect_node_bifurcation*/


LOCAL int untrack_B_reflect_node(
	Front		*fr,
	Wave		*wave,
	O_CURVE		*oldcinc,
	O_CURVE		*cinc,
	O_CURVE		*oldcrefl,
	O_CURVE		*crefl,
	ANGLE_DIRECTION	i_to_a_dir,
	double		dt,
	double		*dt_frac,
	RPROBLEM	*rp,
	NODE_FLAG	flag)
{
	COMPONENT	newcomp;
	NODE		*oldn = Node_of_o_curve(oldcinc);
	NODE		*newn = Node_of_o_curve(cinc);

	newcomp = (curve_ang_oriented_l_to_r(i_to_a_dir,cinc->orient)) ?
		negative_component(cinc->curve) :
		positive_component(cinc->curve);

	if (crefl->curve != NULL)
	{
	    UNTRACK_FLAG  uflag;
	    boolean       oppn_ss;
	    NODE	  *oppn = Opp_node_of_o_curve(crefl);
	    oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	    set_untrack_flag(uflag,crefl->orient,YES,oppn_ss,YES,YES,YES);
	    (void) untrack_curve(crefl,oldcrefl,newcomp,dt,fr,
				 (POINTER)wave,rp,uflag);
	    oldcrefl->curve = crefl->curve = NULL;
	}

	node_type(newn) = NEUMANN_NODE;
	if (!B_node_propagate(fr,(POINTER)wave,oldn,newn,
				 &rp,dt,dt_frac,flag))
	{
	    screen("ERROR in untrack_B_reflect_node(), "
	           "B_node_propagate() failed\n");
	    clean_up(ERROR);
	}

	return GOOD_NODE;
}		/*end untrack_B_reflect_node*/


/*
*			create_mach_stem():
*
*	This function is responsible for changing the topology of the
*	incident curve to create the mach stem.  We first shift the base
*	of the incident to the base position of the mach stem, and ft_assign
*	the new states at the wall.  We then rotate the points near the
*	base, thus enforcing the new angle for the mach stem.  We finally
*	install a point in the incident corresponding to the new triple
*	point, and split the curve at this point.  
*/

LOCAL void create_mach_stem(
	O_CURVE		*oldcinc,
	O_CURVE		*newcinc,
	O_CURVE		*cahead,
	O_CURVE		*cbehind,
	BOND		*crossbinc,
	Bubble		*bubble,
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE_FLAG	flag)
{
	BOND		*bshift;
	CURVE		*mach_stem;
	CURVE		**curves;
	INTERFACE	*intfc = newcinc->curve->interface;
	HYPER_SURF_ELEMENT	*hseshift;
	Locstate	right_st1, left_st1;
	Locstate	right_st2, left_st2;
	Locstate	ahead = bubble->RP->state[0];
	Locstate	behind = bubble->RP->state[1];
	Locstate	refl_mach = bubble->RP->state[3];
	POINT		*mach_base;
	POINT		*mach_tp;
	POINT		*pend;
	double		aw_n[MAXD];		/* ahead wall normal */
	double		mach_ang = bubble->RP->ang[3];
						/* pos x axis to mach stem */
	double		*aw_t = bubble->aw_t;
	double		coords[MAXD];
	double		normal_dist;
	double		dr[MAXD];
	double		len;
	ANGLE_DIRECTION	i_to_a_dir = Opposite_ang_dir(bubble->RP->ang_dir);
	int		i, dim = fr->rect_grid->dim;

		/* Set states on wall at foot of Mach stem */
		/* TODO: Correct if prop-ahead != forward facing */

	if (curve_ang_oriented_l_to_r(i_to_a_dir,cahead->orient))
	{
	    left_st1 = ahead;
	    right_st1 = return_obst_state();
	}
	else
	{
	    left_st1 = return_obst_state();
	    right_st1 = ahead;
	}
	if (curve_ang_oriented_l_to_r(i_to_a_dir,cbehind->orient))
	{
	    left_st2 = return_obst_state();
	    right_st2 = bubble->mach;
	}
	else
	{
	    left_st2 = bubble->mach;
	    right_st2 = return_obst_state();
	}
	if (debugging("grBnode"))
	{
	    (void) printf("Left Right st1:\n");
	    (*fr->print_state)(left_st1);
	    (*fr->print_state)(right_st1);
	    (void) printf("Left Right st2:\n");
	    (*fr->print_state)(left_st2);
	    (*fr->print_state)(right_st2);
	}

	    	/* Shift newcinc base to mach stem base */

	coords[0] = 0.0;	coords[1] = 0.0;
	mach_base = Point(coords);
	if (bubble->is_node_at_corner)
	{
	    /* Note: base of incident has already been shifted past corner,
	     * so that much must be subtracted from bubble->refl_length. */

	    len = bubble->refl_length -
		  distance_between_positions(
				       Coords(Node_of_o_curve(cahead)->posn),
				       bubble->corner_posn,dim);
	}
	else
	{
	    /* Note: base of incident has NOT been shifted, so we need
	     * to add the length to the corner to bubble->refl_length. */

	    len = bubble->refl_length +
		  distance_between_positions(
				       Coords(Node_of_o_curve(cahead)->posn),
				       bubble->corner_posn,dim);
	}
	states_at_distance_along_curve(Node_of_o_curve(cahead)->posn,
		                       Bond_at_node_of_o_curve(cahead),
		                       cahead->curve,cahead->orient,len,1,
				       NULL,NULL,NULL,&hseshift,
				       NULL,&mach_base,fr);
	bshift = Bond_of_hse(hseshift);
	for (i = 0; i < dim; i++)
	    bubble->mach_posn[i] = Coords(mach_base)[i];
	shift_node_past(mach_base,bshift,cahead->curve,cahead->orient,
		        cbehind->curve,cbehind->orient,i_to_a_dir,
		        Node_of_o_curve(cbehind),fr,flag,
		        left_st1,right_st1,left_st2,right_st2);

		/* Modify newcinc at mach stem base */

	if (curve_ang_oriented_l_to_r(i_to_a_dir,newcinc->orient))
	{
	    left_st1 = bubble->mach;
	    right_st1 = ahead;
	}
	else
	{
	    left_st1 = ahead;
	    right_st1 = bubble->mach;
	}
	if (debugging("grBnode"))
	{
	    (void) printf("Before call to modify -- crossbinc:\n");
	    print_bond(crossbinc);
	}
	(void) propagate_curve_near_node(Node_of_o_curve(oldcinc),
		                         Node_of_o_curve(newcinc),
					 oldcinc,newcinc,crossbinc,
		                         left_st1,right_st1,NO,
					 normalized_angle(mach_ang + PI),
		                         fr,wave,dt,flag);

		/* Split incident shock at Mach triple point */

	if (i_to_a_dir == CLOCKWISE)
	{
	    aw_n[0] = -aw_t[1]; aw_n[1] = aw_t[0];
	}
	else
	{
	    aw_n[0] = aw_t[1]; aw_n[1] = -aw_t[0];
	}

	for (i = 0; i < dim; i++)
	{
	    bubble->refl_posn[i] = Coords(mach_base)[i] +
				   aw_n[i] * bubble->mach_height;
	}
	mach_tp = Point(bubble->refl_posn);
	for (bshift = Bond_at_node_of_o_curve(newcinc); bshift != NULL;
	     bshift = Following_bond(bshift,newcinc->orient))
	{
	    pend = Point_of_bond(bshift,Opposite_orient(newcinc->orient));

	    for (i = 0; i < dim; i++)
	    	dr[i] = Coords(pend)[i] - Coords(mach_base)[i];
	    normal_dist = scalar_product(dr,aw_n,dim);
	    if (normal_dist >= bubble->mach_height)
	    {
	        if (insert_point_in_bond(mach_tp,bshift,newcinc->curve) !=
		    FUNCTION_SUCCEEDED)
	        {
	            screen("ERROR in create_mach_stem(), "
		           "insert_point_in_bond() failed\n");
	            clean_up(ERROR);
	        }
	        break;
	    }
	    else
	    {
	        for (i = 0; i < dim; i++)
	            Coords(pend)[i] -= scalar_product(dr,aw_t,dim) * aw_t[i];
	    }
	}
	if (!bshift)
	{
	    screen("ERROR in create_mach_stem(), "
		   "inc too short to construct Mach triple pt\n");
	    clean_up(ERROR);
	}
	if (debugging("grBnode"))
	{
	    (void) printf("posn of Mach triple point = <%g, %g>\n",
			  Coords(mach_tp)[0],Coords(mach_tp)[1]);
	    (void) printf("bshift:\n");
	    print_bond(bshift);
	    (void) printf("newcinc before split\n");
	    (void) printf("newcinc\n");
	    print_o_curve(newcinc);
	}
	interpolate_intfc_states(intfc) = NO;
	set_interpolate_states_at_split_curve_node(NO);
	curves = split_curve(mach_tp,bshift,newcinc->curve,
			     negative_component(newcinc->curve),
			     positive_component(newcinc->curve),
			     negative_component(newcinc->curve),
			     positive_component(newcinc->curve));
	interpolate_intfc_states(intfc) = YES;

		/* Set curve and node data */

	if (newcinc->orient == POSITIVE_ORIENTATION)
	{
	    mach_stem = curves[0];
	    newcinc->curve = curves[1]; 
	    end_status(mach_stem) = MACH_STEM;
	    start_status(newcinc->curve) = INCIDENT;
	    for (i = 0; i < dim; i++)
	    {
	    	Node_vel(mach_stem->end)[i] = bubble->node_v[i];
	    	Node_vel(mach_stem->start)[i] = bubble->mach_speed * 
	    		                        bubble->aw_t[i];
	    }
	    if (i_to_a_dir == CLOCKWISE)
	    {
	    	negative_component(mach_stem) = bubble->comp_mach;
	    	set_states_by_interpolation(mach_stem,NULL,NULL,NEGATIVE_SIDE,
				            bubble->mach,refl_mach,fr->sizest);
		ft_assign(left_start_state(newcinc->curve),behind,fr->sizest);
	    }
	    else
	    {
	    	positive_component(mach_stem) = bubble->comp_mach;
	    	set_states_by_interpolation(mach_stem,NULL,NULL,POSITIVE_SIDE,
				            bubble->mach,refl_mach,fr->sizest);
		ft_assign(right_start_state(newcinc->curve),behind,fr->sizest);
	    }
	}
	else
	{
	    newcinc->curve = curves[0]; 
	    mach_stem = curves[1];
	    end_status(newcinc->curve) = INCIDENT;
	    start_status(mach_stem) = MACH_STEM;
	    for (i = 0; i < dim; i++)
	    {
	    	Node_vel(mach_stem->start)[i] = bubble->node_v[i];
	    	Node_vel(mach_stem->end)[i] = bubble->mach_speed * 
	    		                      bubble->aw_t[i];
	    }
	    if (i_to_a_dir == CLOCKWISE)
	    {
	    	positive_component(mach_stem) = bubble->comp_mach;
	    	set_states_by_interpolation(mach_stem,NULL,NULL,POSITIVE_SIDE,
				            refl_mach,bubble->mach,fr->sizest);
		ft_assign(right_end_state(newcinc->curve),behind,fr->sizest);
	    }
	    else
	    {
	    	negative_component(mach_stem) = bubble->comp_mach;
	    	set_states_by_interpolation(mach_stem,NULL,NULL,NEGATIVE_SIDE,
				            refl_mach,bubble->mach,fr->sizest);
		ft_assign(left_end_state(newcinc->curve),behind,fr->sizest);
	    }
	}
	zero_corr_of_hyper_surf(Hyper_surf(mach_stem));
	node_type(curves[0]->end) = MACH_NODE;
	copy_RP_DATA_structure(Rp_data(curves[0]->end),bubble->RP);
	bubble->RP = Rp_data(curves[0]->end);
}		/*end create_mach_stem*/


/*
*			create_mach_contact():
*
*	This function creates the slip line for a new mach configuration.
*	We first split the behind wall at the base of the contact, and 
*	then simply make_curve() between the triple point and the base.
*
*	TODO: what is the node velocity at the base?
*/

LOCAL void create_mach_contact(
	O_CURVE		*newcinc,
	O_CURVE		*refl_wall,
	Bubble		*bubble,
	Front		*fr)
{
	BOND		*bshift;
	CURVE		*corner_cont_wall;
	CURVE		*cont_stem_wall;
	CURVE		*contact;
	CURVE		**curves;
	INTERFACE	*intfc = newcinc->curve->interface;
	HYPER_SURF_ELEMENT	*hseshift;
	RP_DATA		*RP = bubble->RP;
	Locstate	refl_mach = RP->state[3];
	Locstate	refl_bow = RP->state[2];
	NODE		*node_of_contact_base;
	POINT		*cont_base;
	ANGLE_DIRECTION	i_to_a_dir = Opposite_ang_dir(RP->ang_dir);
	int		i, dim = fr->rect_grid->dim;
	boolean		sav_scss = interpolate_states_at_split_curve_node();
	double		cont_ang = RP->ang[2];	/* pos x axis to contact */
	double		mach_ang = RP->ang[3];	/* pos x axis to mach stem */
	double		mach_cont_dist;		/* corner to base of contact */
	double		coords[MAXD];

		/* Determine location of new nodes on wall */

	mach_cont_dist = bubble->mach_height *
		fabs(tan(cont_ang - mach_ang));
	if (debugging("grBnode"))
	{
	    (void) printf("dist corner to mach = %g mach_cont_dist = %g\n",
			  bubble->refl_length,mach_cont_dist); 
	}

	    /* Split wall at contact base point */

	coords[0] = 0.0;	coords[1] = 0.0;
	cont_base = Point(coords);
	states_at_distance_along_curve(Node_of_o_curve(refl_wall)->posn,
		                       Bond_at_node_of_o_curve(refl_wall),
		                       refl_wall->curve,refl_wall->orient,
				       mach_cont_dist,1,NULL,NULL,NULL,
				       &hseshift,NULL,&cont_base,fr);
	bshift = Bond_of_hse(hseshift);
	for (i = 0; i < dim; i++)
	    bubble->slip_posn[i] = Coords(cont_base)[i];
	if (insert_point_in_bond(cont_base,bshift,refl_wall->curve) !=
	    FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in create_mach_contact(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	interpolate_intfc_states(intfc) = NO;
	set_interpolate_states_at_split_curve_node(NO);
	curves = split_curve(cont_base,bshift,refl_wall->curve,
			     negative_component(refl_wall->curve),
			     positive_component(refl_wall->curve),
			     negative_component(refl_wall->curve),
			     positive_component(refl_wall->curve));
	set_interpolate_states_at_split_curve_node(sav_scss);
	interpolate_intfc_states(intfc) = YES;
	zero_corr_of_hyper_surf(Hyper_surf(curves[0]));
	zero_corr_of_hyper_surf(Hyper_surf(curves[1]));

	    /* Set curve and node data */

	if (refl_wall->orient == POSITIVE_ORIENTATION)
	{
	    cont_stem_wall = curves[0];
	    corner_cont_wall = curves[1];
	}
	else
	{
	    corner_cont_wall = curves[0];
	    cont_stem_wall = curves[1];
	}

	if (i_to_a_dir == CLOCKWISE)
	{
	    if (refl_wall->orient == POSITIVE_ORIENTATION)
	    {
	        positive_component(cont_stem_wall) = bubble->comp_mach;
	        positive_component(corner_cont_wall) = bubble->comp_bow;
	        ft_assign(right_start_state(corner_cont_wall),
	    	       bubble->contact_bow,fr->sizest);
	        if (bubble->is_node_at_corner)
	    	    ft_assign(right_end_state(corner_cont_wall),
	    		   bubble->refl_corner,fr->sizest);
	        ft_assign(right_end_state(cont_stem_wall),
	    	       bubble->contact_mach,fr->sizest);
	    }
	    else
	    {
	        negative_component(cont_stem_wall) = bubble->comp_mach;
	        negative_component(corner_cont_wall) = bubble->comp_bow;
	        ft_assign(left_end_state(corner_cont_wall),
	    	       bubble->contact_bow,fr->sizest);
		if (bubble->is_node_at_corner)
		    ft_assign(left_start_state(corner_cont_wall),
			   bubble->refl_corner,fr->sizest);
		ft_assign(left_start_state(cont_stem_wall),
		       bubble->contact_mach,fr->sizest);
	    }
	}
	else
	{
	    if (refl_wall->orient == POSITIVE_ORIENTATION)
	    {
	        negative_component(cont_stem_wall) = bubble->comp_mach;
	        negative_component(corner_cont_wall) = bubble->comp_bow;
	        ft_assign(left_start_state(corner_cont_wall),
	    	       bubble->contact_bow,fr->sizest);
	        if (bubble->is_node_at_corner)
	    	    ft_assign(left_end_state(corner_cont_wall),
	    		   bubble->refl_corner,fr->sizest);
	        ft_assign(left_end_state(cont_stem_wall),
	    	       bubble->contact_mach,fr->sizest);
	    }
	    else
	    {
	        positive_component(cont_stem_wall) = bubble->comp_mach;
	        positive_component(corner_cont_wall) = bubble->comp_bow;
	        ft_assign(right_end_state(corner_cont_wall),
	    	       bubble->contact_bow,fr->sizest);
	        if (bubble->is_node_at_corner)
	    	    ft_assign(right_start_state(corner_cont_wall),
	    		   bubble->refl_corner,fr->sizest);
	        ft_assign(right_start_state(cont_stem_wall),
	    	       bubble->contact_mach,fr->sizest);
	    }
	}

	end_status(curves[0]) = FIXED;
	start_status(curves[1]) = FIXED;
	node_type(cont_stem_wall->start) = NEUMANN_NODE;
	node_type(cont_stem_wall->end) = NEUMANN_NODE;
	node_of_contact_base = curves[0]->end;

		/* Make contact */

	if (i_to_a_dir == CLOCKWISE)
	{
	    contact = make_curve(bubble->comp_mach,bubble->comp_bow,
	    		         Node_of_o_curve(newcinc),node_of_contact_base);
	    set_states_by_interpolation(contact,NULL,NULL,NEGATIVE_SIDE,
					refl_mach,bubble->contact_mach,
					fr->sizest);
	    set_states_by_interpolation(contact,NULL,NULL,POSITIVE_SIDE,
					refl_bow,bubble->contact_bow,
					fr->sizest);
	}
	else
	{
	    contact = make_curve(bubble->comp_bow,bubble->comp_mach,
	    	                 Node_of_o_curve(newcinc),node_of_contact_base);
	    set_states_by_interpolation(contact,NULL,NULL,NEGATIVE_SIDE,
					refl_bow,bubble->contact_bow,
					fr->sizest);
	    set_states_by_interpolation(contact,NULL,NULL,POSITIVE_SIDE,
					refl_mach,bubble->contact_mach,
					fr->sizest);
	}
	wave_type(contact) = CONTACT;
	start_status(contact) = SLIP;
	end_status(contact) = INCIDENT;
	refl_wall->curve = corner_cont_wall;

	for (i = 0; i < dim; i++)
	    Node_vel(contact->end)[i] = bubble->contact_speed*bubble->aw_t[i];

}		/*end create_mach_contact*/


/*
*			set_corner_for_bifurcation():
*
*	This function is used to compute the node velocity as a reflected
*	shock bifurcates to a more complicated structure.  Computed quantities
*	are corner, bcorner, node_v, and dt_frac.
*
*	Note: dt_frac is the fraction of dt *already used* for the incident
*	to move from it old position to the corner, ie
*		oldposn[] + dt_frac*dt*inc_speed*shock_nor[] = corner_posn[].
*	It is not the fraction of dt since the bifurcation.
*/

EXPORT void set_corner_for_bifurcation(
	O_CURVE		*oldcphys,
	O_CURVE		*newcphys,
	O_CURVE		*oldca,
	O_CURVE		*newca,
	O_CURVE		*oldcb,
	O_CURVE		*newcb,
	POINT		*oldposn,
	POINT		*newposn,
	Front		*fr,
	Wave		*wave,
	POINT		**corner,
	BOND		**bcorner,
	int		forward_facing_is_ahead,
	double		*node_v,
	RP_DATA		*RP,
	double		dt,
	double		*dt_frac)
{
	double		a_t[MAXD], b_t[MAXD];
	double		s_t[MAXD], speed;
	double		corner_angle;		/* ahead to behind walls */
	double		ahead_wall_ang;		/* pos x axis to ahead wall */
	double		behind_wall_ang;
	double		dp[MAXD], dc[MAXD];
	double		num, den;
	NODE		*n = Node_of_o_curve(newcphys);
	CURVE		*nextbc;
	BOND		*bwall, *bfollow;
	POINT		*next_corner;
	int		i, dim = fr->interf->dim;
	ORIENTATION	nextbc_orient;

	    /* Find corner and corner angle */

	debug_print("gBnode","Entered set_corner_for_bifurcation\n");
	if (debugging("gBnode"))
	{
	    (void) printf("oldposn = <%g, %g>, newposn = <%g, %g>\n",
	    	          Coords(oldposn)[0],Coords(oldposn)[1],
	    	          Coords(newposn)[0],Coords(newposn)[1]);
	}
	find_tangent_to_propagated_curve(newposn,
		                         Bond_at_node_of_o_curve(newcphys),
					 oldcphys,newcphys,s_t,fr,
					 (POINTER)wave,dt);
	if (debugging("gBnode"))
	    (void) printf("s_t = <%g, %g>\n",s_t[0],s_t[1]);
	for (i = 0; i < dim; i++)
	    dp[i] = Coords(newposn)[i] - Coords(oldposn)[i];
	if (Node_of_o_curve(newca) != Node_of_o_curve(newcb))
	{
	        /* newcphys has crossed fixed node = corner */

	    debug_print("gBnode","Case of fixed node = corner\n");
	    Check_return(
		next_boundary(newcb->curve,newcb->orient,&nextbc,
			      &nextbc_orient),
		set_corner_for_bifurcation)
	    corner_angle = angle_from_c1_to_c2_at_common_node(
	            newca->curve,newca->orient,newcb->curve,
	            newcb->orient,fr);
	    find_tangent_to_propagated_curve(Node_of_o_curve(newca)->posn,
	                                     Bond_at_node_of_o_curve(newca),
	                                     oldca,newca,a_t,fr,
					     (POINTER)wave,dt);
	    find_tangent_to_propagated_curve(Node_of_o_curve(newcb)->posn,
	                                     Bond_at_node_of_o_curve(newcb),
					     oldcb,newcb,b_t,fr,
					     (POINTER)wave,dt);
	    ahead_wall_ang = angle(a_t[0],a_t[1]);
	    behind_wall_ang = angle(b_t[0],b_t[1]);
	    if (node_type(n) == NEUMANN_NODE)
	    {
	        debug_print("gBnode","Neumann node case\n");
	        *corner = Node_of_o_curve(newcb)->posn;
	        *bcorner = Bond_at_node_of_o_curve(newcb);
	    }
	    else /* REGULAR_REFLECTION */
	    {
	        debug_print("gBnode","Regular reflection node case\n");
	        *corner = 
	        Node_of(nextbc,Opposite_orient(nextbc_orient))->posn;
	        *bcorner = Bond_at_node(nextbc,
	                                Opposite_orient(nextbc_orient));
	    }
	    if (debugging("gBnode"))
	    {
		(void) printf("corner = <%g, %g>\n",
			      Coords((*corner))[0],Coords((*corner))[1]);
		(void) printf("a_t = <%g, %g>, b_t = <%g, %g>\n",a_t[0],a_t[1],
			      b_t[0],b_t[1]);
		print_angle("ahead_wall_ang =",ahead_wall_ang,", ");
		print_angle("behind_wall_ang =",behind_wall_ang,"\n");
	    }
	}
	else
	{
	    /* No intermediate node in bifurcation */

	    debug_print("gBnode","Case of no fixed node at corner\n");
	    find_tangent_to_propagated_curve(Node_of_o_curve(newcb)->posn,
	                                     Bond_at_node_of_o_curve(newcb),
	                                     oldcb,newcb,b_t,fr,
					     (POINTER)wave,dt);
	    debug_print("gBnode","b_t = <%g, %g>\n",b_t[0],b_t[1]);
	    behind_wall_ang = angle(b_t[0],b_t[1]);
	    switch (node_type(n))
	    {
	    case NEUMANN_NODE:
	        debug_print("gBnode","Start bwall loop\n");
	        find_tangent_to_propagated_curve(Node_of_o_curve(newca)->posn,
	                                     Bond_at_node_of_o_curve(newca),
	                                     oldca,newca,a_t,fr,
					     (POINTER)wave,dt);
	       corner_angle = angle_from_c1_to_c2_at_common_node(newca->curve,
							     newca->orient,
	                                                     oldcb->curve,
							     newcb->orient,
							     fr);
	        ahead_wall_ang = angle(a_t[0],a_t[1]);

	        /* Set corner and bcorner */

	        for (bwall = Bond_at_node_of_o_curve(newcb);
	             Following_bond(bwall,newcb->orient);
	             bwall = Following_bond(bwall,newcb->orient))
	        {
	            *corner = (newcb->orient == POSITIVE_ORIENTATION) ?
			      bwall->start : bwall->end;
	            debug_print("gBnode","corner = <%g, %g>\n",
	                  Coords((*corner))[0],Coords((*corner))[1]);
	            *bcorner = bwall;
	            if (inc_angle(s_t,Following_bond(bwall,newcb->orient),
	                          Opposite_orient(newcb->orient),dim)
	                >= 0.5*PI - ERROR_ANGLE) 
	                break;
	            next_corner = (newcb->orient == POSITIVE_ORIENTATION) ?
	                          bwall->end : bwall->start;
		    for (i = 0; i < dim; i++)
		        dc[i] = Coords(next_corner)[i] - Coords(oldposn)[i];
		    if (scalar_product(dc,b_t,dim)*scalar_product(dp,b_t,dim)
			<=0.0)
	                break;
	        }
	        debug_print("gBnode","End of bwall loop\n");
	        if (Following_bond(bwall,newcb->orient))
	            *bcorner = Following_bond(bwall,newcb->orient);
	        break;
	    case B_REFLECT_NODE:
	        bwall = *bcorner;
	        find_tangent_to_propagated_curve(newposn,bwall,oldca,newca,
					         a_t,fr,(POINTER)wave,dt);
	        ahead_wall_ang = angle(a_t[0],a_t[1]);
	        corner_angle = normalized_angle(angle(b_t[0],b_t[1]) -
	                       angle(a_t[0],a_t[1]));
	        if (corner_angle > PI) 
	            corner_angle = 2.0 * PI - corner_angle;
	        debug_print("gBnode","Start bwall loop\n");

	            /* set corner and bcorner */
	        do
	        {
	            /* NOTE: Shift node and cut curve haven't been
	             * called, so the loop runs backward on the 
	             * ahead curve */

		    bfollow = Following_bond(bwall,
			                     Opposite_orient(newca->orient));
	    	    *corner = (newca->orient == POSITIVE_ORIENTATION) ?
			      bwall->start : bwall->end;
		    for (i = 0; i < dim; i++)
		        dc[i] = Coords(*corner)[i] - Coords(oldposn)[i];
		    (void) vector_product(s_t,dc,&num,dim);
		    (void) vector_product(s_t,dp,&den,dim);
		    *dt_frac = num/den;
	            if (*dt_frac < 0.9)
		    {
		        for (i = 0; i < dim; i++)
		        {
	                    node_v[i] = (Coords(newposn)[i] - Coords(*corner)[i]) /
	                    	        (dt * (1.0 - *dt_frac));
		        }
	            }
	            else
		    {
		        for (i = 0; i < dim; i++)
		        {
	                    node_v[i] = (Coords(*corner)[i] - Coords(oldposn)[i]) /
	                                (dt * *dt_frac);
		        }
	                speed = mag_vector(node_v,dim) *
			        fabs(cos(angle(-s_t[1],s_t[0]) - behind_wall_ang))/
	                        fabs(cos(angle(-s_t[1],s_t[0]) - ahead_wall_ang));
		        for (i = 0; i < dim; i++)
	                    node_v[i] = a_t[i] * speed;
	            }
	            if(is_regular_reflection(node_v,fr,RP))
		        break;
		    if (bfollow)
		    {
		        next_corner = (newca->orient == POSITIVE_ORIENTATION) ?
				      bfollow->start : bfollow->end;
		        for (i = 0; i < dim; i++)
		        {
		    	    dc[i] = Coords(next_corner)[i] - Coords(oldposn)[i];
		        }
		        if (scalar_product(dc,b_t,dim)/scalar_product(dp,b_t,dim)
			    <= 0.0)
		    	    break;
		        bwall = bfollow;
		    }
	        }
	        while (bfollow);		/*end do loop*/

	        debug_print("gBnode","End of bwall loop\n");
	        if (Following_bond(bwall,Opposite_orient(newca->orient)))
	            *bcorner = Following_bond(bwall,
					      Opposite_orient(newca->orient));
	        break;
	    default:
	        screen("ERROR in set_corner_for_bifurcation(), "
		       "unexpected case\n");
	        clean_up(ERROR);
	    }
	}
	for (i = 0; i < dim; i++)
	    dc[i] = Coords(*corner)[i] - Coords(oldposn)[i];
	(void) vector_product(s_t,dc,&num,dim);
	(void) vector_product(s_t,dp,&den,dim);
	*dt_frac = num/den;
	if (*dt_frac < 0.9)
	{
	    for (i = 0; i < dim; i++)
	        node_v[i] = (Coords(newposn)[i] - Coords(*corner)[i]) /
				(dt * (1.0 - *dt_frac));
	}
	else
	{
	    for (i = 0; i < dim; i++)
	        node_v[i] = dc[i] / (dt * *dt_frac);
	    speed = mag_vector(node_v,dim) *
	        fabs(cos(angle(-s_t[1],s_t[0])-behind_wall_ang)) /
	            fabs(cos(angle(-s_t[1],s_t[0]) - ahead_wall_ang));
	    for (i = 0; i < dim; i++)
	    	node_v[i] = a_t[i] * speed;
	}
		/* Is ahead = forward_facing */

	if (!forward_facing_is_ahead)
	    behind_wall_ang = ahead_wall_ang;
	if (debugging("gBnode"))
	{
	    print_angle("Corner angle =",corner_angle,"\n");
	    (void) printf("dt_frac %g, node_v <%g, %g>\n",*dt_frac,node_v[0],
			  node_v[1]);
	}
	debug_print("gBnode","Left set_corner_for_bifurcation\n");
}		/*end set_corner_for_bifurcation*/


/*
*			g_identify_physical_node():
*
*	This routine is called in scalar_unravel. Its purpose is to set the
*	node_type and curve start/end status of new nodes (especially CC_NODES)
*	created within this routine.
*/

EXPORT void g_identify_physical_node(
	NODE		*n)
{
	CURVE		**c;
	int		num_curves = 0;
	int             num_vector = 0;
	int		num_contact = 0;
	int		num_thinflame = 0;
	int             num_neumann = 0;
	int             num_dirichlet = 0;
	int             num_passive = 0;
	int             num_subdomain = 0;
	int             num_phys;
	int             wtype;

	debug_print("identify","Entered g_identify_physical_node()\n");
	if (debugging("identify")) print_node(n);
	for (c = n->in_curves; c && *c; c++)
	{
	    wtype = wave_type(*c);
	    num_curves++;
	    if (is_scalar_wave(wtype))
	    	num_contact++;
	    if (is_thinflame_wave(wtype))
	    	num_thinflame++;
	    if (is_vector_wave(wtype))
		num_vector++;
	    if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	    {
		switch (wtype)
		{
		case NEUMANN_BOUNDARY:
		    num_neumann++;
		    break;
		case DIRICHLET_BOUNDARY:
		    num_dirichlet++;
		    break;
		case PASSIVE_BOUNDARY:
		    num_passive++;
		    break;
		case SUBDOMAIN_BOUNDARY:
		case REFLECTION_BOUNDARY:
		    num_subdomain++;
		    break;
		}
	    }
	}
	for (c = n->out_curves; c && *c; c++)
	{
	    wtype = wave_type(*c);
	    num_curves++;
	    if (is_scalar_wave(wtype))
	    	num_contact++;
	    if (is_vector_wave(wtype))
		num_vector++;
	    if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	    {
		switch (wtype)
		{
		case NEUMANN_BOUNDARY:
		    num_neumann++;
		    break;
		case DIRICHLET_BOUNDARY:
		    num_dirichlet++;
		    break;
		case PASSIVE_BOUNDARY:
		    num_passive++;
		    break;
		case SUBDOMAIN_BOUNDARY:
		case REFLECTION_BOUNDARY:
		    num_subdomain++;
		    break;
		}
	    }
	}
	num_phys = num_contact + num_thinflame + num_vector;
	if (num_subdomain > 0)
	{
	    node_type(n) = SUBDOMAIN_NODE;
	    for (c = n->in_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            end_status(*c) = PASSIVE;
		else if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	            end_status(*c) = FIXED;
		else
	            end_status(*c) = INCIDENT;
	    }
	    for (c = n->out_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            start_status(*c) = PASSIVE;
		else if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	            start_status(*c) = FIXED;
		else
	            start_status(*c) = INCIDENT;
	    }
	}
	else if (num_phys == 0)
	{
	    if ((num_neumann > 0) || (num_dirichlet > 0))
	        node_type(n) = FIXED_NODE;
	    else
	        node_type(n) = PASSIVE_NODE;
	    for (c = n->in_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            end_status(*c) = PASSIVE;
		else
	            end_status(*c) = FIXED;
	    }
	    for (c = n->out_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            start_status(*c) = PASSIVE;
		else
	            start_status(*c) = FIXED;
	    }
	}
	else if ((num_neumann > 0) && (num_dirichlet > 0))
	{
	    node_type(n) = ATTACHED_B_NODE;
	    for (c = n->in_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            end_status(*c) = PASSIVE;
		else if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	            end_status(*c) = FIXED;
		else
	            end_status(*c) = INCIDENT;
	    }
	    for (c = n->out_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            start_status(*c) = PASSIVE;
		else if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	            start_status(*c) = FIXED;
		else
	            start_status(*c) = INCIDENT;
	    }
	}
	else if (((num_neumann > 0) || (num_dirichlet > 0)) && (num_phys == 1))
	{
	    node_type(n) = (num_neumann > 0) ? NEUMANN_NODE : DIRICHLET_NODE;
	    for (c = n->in_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            end_status(*c) = PASSIVE;
		else if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	            end_status(*c) = FIXED;
		else
	            end_status(*c) = INCIDENT;
	    }
	    for (c = n->out_curves; c && *c; c++)
	    {
	        wtype = wave_type(*c);
		if (wtype == PASSIVE_BOUNDARY)
	            start_status(*c) = PASSIVE;
		else if (wtype < FIRST_PHYSICS_WAVE_TYPE)
	            start_status(*c) = FIXED;
		else
	            start_status(*c) = INCIDENT;
	    }
	}
	else if ((num_curves == 3) && (num_contact == 3))
	{
	    node_type(n) = CC_NODE;
	    /* TODO: reset with two SLIP and one INCIDENT */
	    for (c = n->in_curves; c && *c; c++)
	        end_status(*c) = SLIP;
	    for (c = n->out_curves; c && *c; c++)
	        start_status(*c) = SLIP;
	}
	else
	{
	    node_type(n) = ERROR;
	    debug_print("identify","Left g_identify_physical_node\n");
	    return;
	}
	debug_print("identify","Left g_identify_physical_node\n");
}		/*end g_identify_physical_node*/



LOCAL	void insert_normal_contact_at_wall(
	Front		*fr,
	O_CURVE		*oldcphys,
	O_CURVE		*newcphys,
	O_CURVE		*oldca,
	O_CURVE		*oldcb)
{
	INTERFACE	*intfc;
	CWNP		*cwnp;
	POINT		*p, *pnd;
	double		*h = fr->rect_grid->h;
	double		W[MAXD], t[MAXD], nor_a[MAXD], nor_b[MAXD], nor[MAXD];
	double		len, alpha;
	double		pjump, curvature;
	boolean		sav_intrp;
	int		i, dim = fr->rect_grid->dim;

	cwnp = contact_wall_node_params(newcphys->curve->interface);
	if ((cwnp->adjust != YES) ||
			(newcphys->curve->num_points < 3) ||
			(fr->step < cwnp->first_adjust_step) ||
			(fr->time < cwnp->first_adjust_time))
		return;
	find_tangent_to_curve(Node_of_o_curve(oldcphys)->posn,
		              Bond_at_node_of_o_curve(oldcphys),oldcphys->curve,
		              oldcphys->orient,t,fr);
	normal(Node_of_o_curve(oldca)->posn,
	       Hyper_surf_element(Bond_at_node_of_o_curve(oldca)),
	       Hyper_surf(oldca->curve),nor_a,fr);
	normal(Node_of_o_curve(oldcb)->posn,
	       Hyper_surf_element(Bond_at_node_of_o_curve(oldcb)),
	       Hyper_surf(oldcb->curve),nor_b,fr);
	for (i = 0; i < dim; i++)
	    nor[i] = 0.5*(nor_a[i]+nor_b[i]);
	len = mag_vector(nor,dim);
	for (i = 0; i < dim; i++)
	    nor[i] /= len;
	if (scalar_product(t,nor,dim) < 0.0)
	    for (i = 0; i < dim; i++)
		nor[i] = -nor[i];
	pnd = Node_of_o_curve(newcphys)->posn;
	p = Point_adjacent_to_node(newcphys->curve,newcphys->orient);
	len = scaled_separation(p,pnd,h,dim);
	if (len < 1.5*cwnp->wall_bond_len)
	{
	    alpha = len/scaled_hypot(nor,h,dim);
	    for (i = 0; i < dim; i++)
	    	Coords(p)[i] = Coords(pnd)[i] + alpha*nor[i];
	}
	else
	{
	    intfc = newcphys->curve->interface;
	    sav_intrp = interpolate_intfc_states(intfc);
	    interpolate_intfc_states(intfc) = YES;
	    alpha = cwnp->wall_bond_len/scaled_hypot(nor,h,dim);
	    p = Point(Coords(pnd));
	    for (i = 0; i < dim; i++)
	    	Coords(p)[i] += alpha*nor[i];
	    insert_point_adjacent_to_node(p,newcphys->curve,
	    			          newcphys->orient);
	    interpolate_intfc_states(intfc) = sav_intrp;
	}
	if (newcphys->curve->num_points < 5)
	    return;
	normal(p,Hyper_surf_element(Bond_at_node_of_o_curve(newcphys)),
	       Hyper_surf(newcphys->curve),nor,fr);
	if (fabs(surface_tension(newcphys->curve)) >= MACH_EPS)
	{
	    pjump = 2.0*(p->curvature)*surface_tension(newcphys->curve);
	    /*
	    double tension = surface_tension(newcphys->curve);
	    switch (dim)
	    {
	    
	    case 2:
		curvature = mean_curvature_at_point(p,
			Hyper_surf_element(Bond_at_node_of_o_curve(newcphys)),
			Hyper_surf(newcphys->curve),fr);
	        pjump = tension*curvature;
		break;
            case 3:
	        pjump = 2.0*tension*(p->curvature);
		break;
	    }*/
	}
	else
	    pjump = 0.0;

	w_speed(Coords(p),left_state(p),right_state(p),
		left_state(p),right_state(p),
		W,pjump,nor,wave_type(newcphys->curve),fr);

}		/*end insert_normal_contact_at_wall*/

/*
*			inc_angle():
*
*	Computes the incident angle, defined to be the angle between the
*	bond (with orientation) bwall and the vector s_tx, s_ty. The answer
*	is normalized to lie in [0,PI).
*/

LOCAL double inc_angle(
	double		*s_t,
	BOND		*bwall,
	ORIENTATION	bwall_orient,
	int		dim)
{
	double		x_to_i_angle,x_to_w_angle,ans;
	double		w_t[MAXD], len;
	int		i;

	debug_print("gBnode","Entered inc_angle()\n");
	for (i = 0; i < dim; i++)
	    w_t[i] = Coords(bwall->end)[i] - Coords(bwall->start)[i];
	len = separation(bwall->start,bwall->end,dim);
	for (i = 0; i < dim; i++)
	    w_t[i] /= len;
	if (bwall_orient == NEGATIVE_ORIENTATION)
	{
	    for (i = 0; i < dim; i++)
	    	w_t[i] *= -1.;
	}
	x_to_i_angle = angle(s_t[0],s_t[1]);
	x_to_w_angle = angle(w_t[0],w_t[1]);
	ans = normalized_angle(x_to_i_angle - x_to_w_angle);
	if (ans >= PI) ans = 2.0 * PI - ans;
	if (debugging("gBnode"))
	{
	    (void) printf("s_t = <%g, %g>\n",s_t[0],s_t[1]);
	    (void) printf("w_t = <%g, %g>\n",w_t[0],w_t[1]);
	    print_angle("x_to_i_angle =",x_to_i_angle,", ");
	    print_angle("x_to_w_angle =",x_to_w_angle,"\n");
	    print_angle("ans =",ans,"\n");
	}
	debug_print("gBnode","Left inc_angle()\n");
	return ans;
}		/*end inc_angle*/


/*
*			init_bow():
*
*	This function creates the reflected or bow shock when initializing
*	a reflection (reg or mach).
*/

LOCAL void init_bow(
	COMPONENT	compbow_ext,
	O_CURVE		*cbehind,
	O_CURVE		*bow_shock,
	O_CURVE		*bow_wall,
	O_CURVE		*refl_wall,
	BOND		*bcorner,
	Locstate	*exterior_bow_state,
	NODE		*ref_node,
	Bubble		*bubble,
	Front		*fr,
	int		is_mach_refl,
	boolean		track_rs)
{
	COMPONENT	comp_bow = bubble->comp_bow;
	INTERFACE	*intfc = cbehind->curve->interface;
	RECT_GRID	*gr = fr->rect_grid;
	NODE		*ne;
	CURVE		**curves;
	BOND		*bbow;
	POINT		*pbow, *pmid;
	int		w_type = wave_type(cbehind->curve);
	int		i, dim = fr->rect_grid->dim;
	boolean		sav_intrp;
	boolean		sav_copy;
	ANGLE_DIRECTION	i_to_a_dir = Opposite_ang_dir(bubble->RP->ang_dir);
	boolean		sav_scss = interpolate_states_at_split_curve_node();
	double		len, bbow_len;
	double		coords[MAXD];
	double		refl_angle;

	debug_print("gBnode","Entered init_bow\n");

	    /* Find crossing of cbehind with bow wave */

	len = bubble->bow_length;
	bbow = bcorner;
	bbow_len = separation(bbow->start,bbow->end,gr->dim);
	while (len > bbow_len)
	{
	    if (Following_bond(bbow,cbehind->orient) == NULL)
		break;
	    len -= bbow_len;
	    bbow = Following_bond(bbow,cbehind->orient);
	    bbow_len = separation(bbow->start,bbow->end,gr->dim);
	}
	len = min(len,bbow_len);
	if (cbehind->orient == POSITIVE_ORIENTATION)
	{
	    for (i = 0; i < dim; i++)
	    	bubble->bow_posn[i] =
		    Coords(bbow->start)[i] + (len/bbow_len)*
		    (Coords(bbow->end)[i] - Coords(bbow->start)[i]);
	}
	else
	{
	    for (i = 0; i < dim; i++)
	    	bubble->bow_posn[i] =
	    	    Coords(bbow->end)[i] + (len/bbow_len)*
		    (Coords(bbow->start)[i] - Coords(bbow->end)[i]);
	}
	pbow = Point(bubble->bow_posn);

	    /* Split cbehind at crossing with bow wave */

	sav_intrp = interpolate_intfc_states(intfc);
	sav_copy = copy_intfc_states();
	interpolate_intfc_states(intfc) = YES;
	if (insert_point_in_bond(pbow,bbow,cbehind->curve)!=FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in init_bow(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	bow_wall->orient = cbehind->orient;
	if (!track_rs)
	{
	    bow_wall->curve = cbehind->curve;
	    *exterior_bow_state =
		    (curve_ang_oriented_l_to_r(i_to_a_dir,bow_wall->orient)) ?
		    right_state(pbow) : left_state(pbow);
	}
	else
	{
	    set_copy_intfc_states(YES);
	    interpolate_intfc_states(intfc) = NO;
	    set_interpolate_states_at_split_curve_node(NO);
	    curves = split_curve(pbow,bbow,cbehind->curve,
	    		         negative_component(cbehind->curve),
				 positive_component(cbehind->curve),
				 negative_component(cbehind->curve),
				 positive_component(cbehind->curve));
	    set_interpolate_states_at_split_curve_node(sav_scss);
	    interpolate_intfc_states(intfc) = sav_intrp;
	    set_copy_intfc_states(sav_copy);
	    ne = curves[0]->end;
	    switch (w_type)
	    {
	    case NEUMANN_BOUNDARY:
	        node_type(ne) = NEUMANN_NODE;
	    	break;
	    case DIRICHLET_BOUNDARY:
	    	node_type(ne) = DIRICHLET_NODE;
	    	break;
	    case SUBDOMAIN_BOUNDARY:
	    	node_type(ne) = SUBDOMAIN_NODE;
	    	break;
	    default:
	    	screen("ERROR in init_bow(), invalid wave type\n");
	    	clean_up(ERROR);
	    	break;
	    }
	    end_status(curves[0]) = start_status(curves[1]) = FIXED;
	    debug_print("gBnode","Setting bow_wall data\n");
	    if (cbehind->orient == POSITIVE_ORIENTATION)
	    {
	    	bow_wall->curve = curves[0];
	    	zero_corr_of_hyper_surf(Hyper_surf(curves[0]));
	    	if (i_to_a_dir == CLOCKWISE)
	    	{
	    	    positive_component(bow_wall->curve) = comp_bow;
	    	    *exterior_bow_state = right_start_state(curves[1]);
	        }
		else
		{
		    negative_component(bow_wall->curve) = comp_bow;
		    *exterior_bow_state = left_start_state(curves[1]);
		}
	    }
	    else
	    {
	    	bow_wall->curve = curves[1];
	    	zero_corr_of_hyper_surf(Hyper_surf(curves[1]));
	    	if (i_to_a_dir == CLOCKWISE)
	    	{
	    	      negative_component(bow_wall->curve) = comp_bow;
	    	      *exterior_bow_state = left_end_state(curves[0]);
	    	}
	    	else
	    	{
	    	      positive_component(bow_wall->curve) = comp_bow;
	    	      *exterior_bow_state = right_end_state(curves[0]);
	    	}
	    }

	    	/* Initialize bow wave */

	    len = separation(ref_node->posn,pbow,gr->dim);
	    refl_angle = (is_mach_refl) ?
	    	bubble->RP->ang[1] : bubble->RP->ang[2];
	    coords[0] = Coords(ref_node->posn)[0] + 0.5*len*cos(refl_angle);
	    coords[1] = Coords(ref_node->posn)[1] + 0.5*len*sin(refl_angle);
	    pmid = Point(coords);
	    if (curve_ang_oriented_l_to_r(i_to_a_dir,bow_shock->orient))
	    {
	    	bow_shock->curve =
		    make_half_bow_curve(comp_bow,compbow_ext,
					2*max(gr->gmax[0],gr->gmax[1]),
				        bow_shock->orient,ref_node,ne,
				        Coords(pmid),gr,bubble);
		wave_type(bow_shock->curve) = FORWARD_SHOCK_WAVE;
	    }
	    else
	    {
	    	bow_shock->curve =
		    make_half_bow_curve(compbow_ext,comp_bow,
					2*max(gr->gmax[0],gr->gmax[1]),
				        bow_shock->orient,ref_node,ne,
				        Coords(pmid),gr,bubble);
		wave_type(bow_shock->curve) = BACKWARD_SHOCK_WAVE;
	    }
	    bow_shock->orient = POSITIVE_ORIENTATION;
	    set_status_at_node(bow_shock->curve,bow_shock->orient,REFLECTED);
	    set_status_at_node(bow_shock->curve,
			       Opposite_orient(bow_shock->orient),INCIDENT);
	    set_no_tan_propagate(bow_shock->curve);
	}
	if (bubble->is_node_at_corner)
	{
	    Check_return(next_boundary(bow_wall->curve,bow_wall->orient,
			               &refl_wall->curve,&refl_wall->orient),
		         init_bow)
	    zero_corr_of_hyper_surf(Hyper_surf(refl_wall->curve));
	    if (curve_ang_oriented_l_to_r(i_to_a_dir,refl_wall->orient))
	    	negative_component(refl_wall->curve) = comp_bow;
	    else
	    	positive_component(refl_wall->curve) = comp_bow;
	}
	else
	{
	    refl_wall->curve = bow_wall->curve;
	    refl_wall->orient = bow_wall->orient;
	}
	debug_print("gBnode","Left init_bow()\n");
}		/*end init_bow*/


/*
*			allocate_bubble():
*
*	Allocates a Bubble structure.  Usually, the RP will already exist and
*	thus be supplied, but it can be allocated here by passing a NULL
*	argument and appropriate state type.
*/

EXPORT	Bubble *allocate_bubble(
	boolean		dynamic,
	Front		*front,
	RP_DATA		*RP,
	int		stype)
{
	Bubble		*bubble;
	size_t		sizest = front->sizest;

	if (dynamic == YES)
	{
	    bubble = (Bubble *) store(sizeof(Bubble));
	    bubble->bow = alloc_intfc_state(front->interf,sizest);
	    bubble->refl_corner = alloc_intfc_state(front->interf,sizest);
	    bubble->bow_corner = alloc_intfc_state(front->interf,sizest);

	    bubble->mach = alloc_intfc_state(front->interf,sizest);
	    bubble->contact_bow = alloc_intfc_state(front->interf,sizest);
	    bubble->contact_mach = alloc_intfc_state(front->interf,sizest);
	}
	else
	{
	    scalar(&bubble,sizeof(Bubble));
	    alloc_state(front->interf,&bubble->bow,sizest);
	    alloc_state(front->interf,&bubble->refl_corner,sizest);
	    alloc_state(front->interf,&bubble->bow_corner,sizest);

	    alloc_state(front->interf,&bubble->mach,sizest);
	    alloc_state(front->interf,&bubble->contact_bow,sizest);
	    alloc_state(front->interf,&bubble->contact_mach,sizest);
	}
	bubble->intfc_table_storage = dynamic;

	if (RP == NULL)
	    bubble->RP = allocate_RP_DATA_structure(sizest,dynamic,stype);
	else
	    bubble->RP = RP;

	return bubble;
}		/*end allocate_bubble*/

/*
*			free_bubble():
*
*	Frees a the storage for a Bubble structure. If free_rp is YES then
*	the bubble->RP structure will also be freed.
*/

EXPORT	void free_bubble(
	Bubble		*bubble,
	boolean		free_rp)
{
	if (bubble->intfc_table_storage == YES)
	    return;

	if (free_rp)
	    free_RP_DATA_structure(bubble->RP);

	free(bubble->bow);
	free(bubble->refl_corner);
	free(bubble->bow_corner);

	free(bubble->mach);
	free(bubble->contact_bow);
	free(bubble->contact_mach);

	free(bubble);
}		/*end free_bubble*/
#endif /* defined(TWOD) */
