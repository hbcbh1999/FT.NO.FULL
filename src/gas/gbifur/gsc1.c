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
*				gsc1.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Solves two dimensional shock contact Riemann problems.
*
*/

#if defined(TWOD) && defined(FULL_PHYSICS)

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	int	angles_consistent_at_diff_node(RP_DATA*,int);
LOCAL	int	identify_curves_at_sc_rp(RPROBLEM*,O_CURVE**,O_CURVE**,
					 O_CURVE**,O_CURVE**,ANGLE_DIRECTION*);
LOCAL	int	n_install_dfrctn_to_bdry(CURVE**,ORIENTATION*,RPROBLEM*,Front*,
					 O_CURVE*,O_CURVE*,ANGLE_DIRECTION);
LOCAL	int	regular_shock_contact_rp(Front*,Wave*,RPROBLEM*,O_CURVE*,
					 O_CURVE*,O_CURVE*,O_CURVE*,
					 ANGLE_DIRECTION);
LOCAL	void	reset_wall_node(POINT*,Locstate,Locstate,BOND*,
				CURVE*,ORIENTATION,
				CURVE*,ORIENTATION,
				ANGLE_DIRECTION,Front*);
LOCAL	void	set_curve_at_shock_diffraction(CURVE**,ORIENTATION,COMPONENT,
					       COMPONENT,Locstate,Locstate,
					       int,int,NODE*,NODE*);


/*
*			shock_contact_rp():
*
*	This is the main driver for the shock contact interaction
*	Riemann problem.  The interacting curves are identified
*	and the configuration of a diffraction node is installed.
*	The curves at the diffraction node are numbered as follows
*
*	0.  Incident shock 
*	1.  Leading edge of reflected rarefaction wave (if it exists)
*	2.  Reflected shock (if it exists)
*	3.  Trailing edge of reflected rarefaction wave (if it exists)
*	4.  Contact behind the reflected wave
*	5.  Transmitted shock.
*	6.  Contact in front of incident shock.
*/



EXPORT int shock_contact_rp(
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp)
{
	O_CURVE		*shock, *oldsh, *contact, *oldct;
	ANGLE_DIRECTION	i_to_f_dir;
	int		status;

	debug_print("shock_contact","Entered shock_contact_rp()\n");
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    NODE **n;

	    (void) printf("Interfaces into shock_contact_rp()\n");
	    (void) printf("Old Interface\n");
	    print_interface(rp->old_intfc);
	    (void) printf("New Interface\n");
	    print_interface(rp->new_intfc);
	    for (n = rp->new_intfc->nodes; n && *n; n++)
	    	print_propagation_status(*n);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	if (!identify_curves_at_sc_rp(rp,&shock,&oldsh,&contact,&oldct,
					 &i_to_f_dir))
	{
	    (void) printf("WARNING in shock_contact_rp(), "
	                  "identify_curves_at_sc_rp() failed\n");
	    return ERROR_IN_STEP;
	}

	status = regular_shock_contact_rp(front,wave,rp,shock,oldsh,contact,
					  oldct,i_to_f_dir);
	if (status != GOOD_STEP)
	{
	    (void) printf("WARNING in shock_contact_rp(), "
	                  "regular_shock_contact_rp() failed\n"
	                  "CODE NEEDED FOR IRREGULAR SHOCK DIFFRACTION\n");
	    debug_print("shock_contact",
		  "Left shock_contact_rp(), status = %s\n",
		  time_step_status_as_string(status));
	    return status;
	}
	debug_print("shock_contact","Left shock_contact_rp(), status = GOOD_STEP\n");
	return GOOD_STEP;
}		/*end shock_contact_rp*/

LOCAL	int identify_curves_at_sc_rp(
	RPROBLEM	*rp,
	O_CURVE		**shock,
	O_CURVE		**oldsh,
	O_CURVE		**contact,
	O_CURVE		**oldct,
	ANGLE_DIRECTION	*i_to_f_dir)
{
	O_CURVE		*oc;

	debug_print("shock_contact","Entered identify_curves_at_sc_rp()\n");
	*shock = NULL;
	*oldsh = NULL;
	*contact = NULL;
	*oldct = NULL;
	if (!rp) 
	{
	    debug_print("shock_contact",
		  "Left identify_curves_at_sc_rp(), ans = YES\n");
	    return YES;
	}


	if (!(rp->ang_ordered_curves && rp->ang_ordered_curves->first &&
		 rp->ang_ordered_curves->first->next))
	{
	    (void) printf("WARNING in identify_curves_at_sc_rp(), "
	                  "Angle ordered curve list not set up\n");
	    debug_print("shock_contact",
		  "Left identify_curves_at_sc_rp(), ans = NO\n");
	    return NO;
	}

	oc = rp->ang_ordered_curves->first;

	if (is_scalar_wave(wave_type(oc->next->curve)))
	{
	    *i_to_f_dir = CLOCKWISE;
	    *contact = oc->next;
	    *oldct = rp->old_ang_ordered_curves->first->next;
	    *shock = oc->next->next;
	    *oldsh = (*oldct)->next;
	}
	else if (is_shock_wave(wave_type(oc->next->curve)))
	{
	    *i_to_f_dir = COUNTER_CLOCK;
	    *shock = oc->next;
	    *oldsh = rp->old_ang_ordered_curves->first->next;
	    *contact = oc->next->next;
	    *oldct = (*oldsh)->next;
	}
	else
	{
	    (void) printf("WARNING in identify_curves_at_sc_rp(), "
	                  "Shock or contact not found\n");
	    debug_print("shock_contact",
	    	  "Left identify_curves_at_sc_rp(), ans = NO\n");
	    return NO;
	}
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("Old incident shock\n");
	    print_o_curve(*oldsh);
	    (void) printf("Old front contact\n");
	    print_o_curve(*oldct);

	    if (debugging("curve_states"))
	    {
	    	(void) printf("States on old incident shock\n");
	    	verbose_print_curve_states((*oldsh)->curve);
	    	(void) printf("States on old incident contact\n");
	    	verbose_print_curve_states((*oldct)->curve);
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	debug_print("shock_contact","Left identify_curves_at_sc_rp(), ans = YES\n");
	return YES;
}		/*end identify_curves_at_sc_rp*/



/*
*			regular_shock_contact_rp():
*
*	This function performs the bifurcation which occurs when
*	an incident shock collides with a contact discontinuity.
* 	It is assumed that the sound speed of the gas in which
*	the incident shock propagates has a faster sound speed
*	than the gas on the other side of the contact.
*	The result of this interaction is a diffraction node.
*
*	The curves at the diffraction node are numbered as follows
*
*	0.  Contact in front of incident shock.
*	1.  Incident shock 
*	2.  Reflected shock (if it exists)
*	3.  Contact behind the reflected wave
*	4.  Transmitted shock.
*/


LOCAL	int regular_shock_contact_rp(
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp,
	O_CURVE		*shock,
	O_CURVE		*oldsh,
	O_CURVE		*contact,
	O_CURVE		*oldct,
	ANGLE_DIRECTION	i_to_f_dir)
{
	BOND		*crossb_sh, *crossb_ct;
	COMPONENT	l_comp[7], r_comp[7];
	COMPONENT	newcomp[7];
	CURVE		**curves;
	CURVE		*dncur[7];
	INTERFACE	*intfc = rp->new_intfc;
	INTERFACE	*save_intfc;
	Locstate	l_st[7], r_st[7];
	NODE		*diffn, *oppn;
	ORIENTATION	c_orient[7];
	NODE_FLAG	flag;
	POINT		*pc;
	RP_DATA		*RP;
	boolean		sav_interp;
	boolean		sav_scss = interpolate_states_at_split_curve_node();
	double		dt = rp->dt;
	double		tcr_sh, tcr_ct;
	double		nod_v[MAXD];
	boolean		is_reflected_shock;
	int		j;
	boolean		track[7];
	int		bdry_type;
	int		w_type[7];
	static	int	status[7] = { INCIDENT,
				      REFLECTED,
				      REFLECTED,
				      REFLECTED,
				      SLIP,
				      TRANSMITTED,
				      CONTACT_TARGET};
	static	double	**t = NULL;

	debug_print("shock_contact","Entered regular_shock_contact_rp()\n");

	/* Allocate storage */

	if (t == NULL)
	    bi_array(&t,2,2,FLOAT);

	if (rp->bdry_type1 != rp->bdry_type2) 
	{
	    (void) printf("WARNING: in regular_shock_contact_rp(), "
	                  "Code needed for different boundary types\n");
	    debug_print("shock_contact",
	          "Left regular_shock_contact_rp(), ans = ERROR_IN_STEP\n");
	    return ERROR_IN_STEP;
	}
	else
	    bdry_type = rp->bdry_type1;

	for (j = 0; j < 7; j++)
	{
	    dncur[j] = NULL;
	    w_type[j] = ERROR;
	    track[j] = NO;
	    l_comp[j] = ERROR;
	    r_comp[j] = ERROR;
	    newcomp[j] = ERROR;
	    l_st[j] = NULL;
	    r_st[j] = NULL;
	}

	save_intfc = current_interface();
	set_current_interface(intfc);
	sav_interp = interpolate_intfc_states(intfc);
	

	RP = allocate_RP_DATA_structure(front->sizest,YES,GAS_STATE);
	RP->ang_dir = Opposite_ang_dir(i_to_f_dir);

	set_orients_and_w_type_for_shock_diffraction(wave_type(shock->curve),
						     shock->orient,
						     wave_type(contact->curve),
						     contact->orient,
						     c_orient,w_type);


	/* Diffraction node position is intersection of propagated curves */

	clear_node_flag(flag);
	single_extend_to_cross(flag) = YES;
	double_extend_to_cross(flag) = YES;
	if (intersection_of_two_o_curves(shock,oldsh,contact,oldct,&crossb_sh,
					 &crossb_ct,&pc,&tcr_sh,&tcr_ct,front,
					 (POINTER)wave,dt,flag) == NO) 
	{
	    set_current_interface(save_intfc);
	    interpolate_intfc_states(intfc) = sav_interp;
	    (void) printf("WARNING in shock_contact_rp(), "
	                  "No cross of propagated shock and contact\n");
	    debug_print("shock_contact",
	          "Left regular_shock_contact_rp(), ans = ERROR_IN_STEP\n");
	    return ERROR_IN_STEP;
	}
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("Propagated curves intersect at (%g, %g)\n",
	    	          Coords(pc)[0],Coords(pc)[1]);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	/* Read off states from the propagated incident curve */

	assign_ahead_states_at_shock_diffraction(tcr_sh,crossb_sh,shock->curve,
						 shock->orient,tcr_ct,crossb_ct,
						 contact->curve,contact->orient,
						 newcomp,RP);

	/* Approximate node velocity */

	find_node_vel_at_rp(pc,tcr_sh,crossb_sh,oldsh,shock,t[0],tcr_ct,
			    crossb_ct,oldct,contact,t[1],RP->ang_dir,
		            DIFFRACTION_NODE,front,wave,dt,nod_v);

	/* Find states around the diffraction node */

	set_to_next_node_only(flag);
	switch (is_regular_diffraction_node(Coords(pc),nod_v,NULL,t,RP,NULL,
					    &is_reflected_shock,front,
					    DIFFRACTION_NODE,flag)) 
	{
	case REGULAR_DIFFRACTION:
	    break;

	case ANOMALOUS_REFLECTION:
	case REGULAR_TO_MACH_DIFFRACTION:
	case ERROR_DIFFRACTION:
	default:
	    set_current_interface(save_intfc);
	    interpolate_intfc_states(intfc) = sav_interp;
	    (void) printf("WARNING in regular_shock_contact_rp(), "
	                  "is_regular_diffraction_node() failed\n"
	                  "possible bifurcation,  CODE NEEDED\n");
	    rp->dt_frac = min(rp->dt_frac,0.5);
	    debug_print("shock_contact",
	          "Left regular_shock_contact_rp(), ans = MODIFY_TIME_STEP\n");
	    return MODIFY_TIME_STEP;
	}

	/* Check for consistent angles */
	
	if (!angles_consistent_at_diff_node(RP,is_reflected_shock))
	{
	    rp->dt_frac = min(rp->dt_frac,
	    		  Min_time_step_modification_factor(front));
	    set_current_interface(save_intfc);
	    interpolate_intfc_states(intfc) = sav_interp;
	    (void) printf("WARNING in regular_shock_contact_rp(), "
	                  "is_regular_diffraction_node() returns "
	                  "inconsistent angles\n"
	                  "possible bifurcation,  CODE NEEDED\n");
	    debug_print("shock_contact",
	          "Left regular_shock_contact_rp(), ans = MODIFY_TIME_STEP\n");
	    return MODIFY_TIME_STEP;
	}

	/* Determine components and states of new curves */

	set_track_and_newcomp_list_for_shock_diffraction(track,newcomp,w_type,
							 RP,is_reflected_shock,
							 front);

	/* Set left and right states and components at diffraction node */

	set_states_and_comps_about_shock_diffraction(newcomp,i_to_f_dir,
						     c_orient,l_comp,r_comp,
						     l_st,r_st,RP);

	/* Split contact into front and back parts */

	interpolate_intfc_states(intfc) = NO;
	set_copy_intfc_states(YES);
	if (insert_point_in_bond(pc,crossb_ct,contact->curve) !=
	    FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in regular_shock_contact_rp(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}

	set_interpolate_states_at_split_curve_node(NO);
	if (contact->orient == POSITIVE_ORIENTATION) 
	{
	    curves = split_curve(pc,crossb_ct,contact->curve,
				 l_comp[4],r_comp[4],l_comp[6],r_comp[6]);
	    dncur[4] = curves[0];
	    dncur[6] = curves[1];
	}
	else 
	{
	    curves = split_curve(pc,crossb_ct,contact->curve,
				 l_comp[6],r_comp[6],l_comp[4],r_comp[4]);
	    dncur[4] = curves[1];
	    dncur[6] = curves[0];
	}
	set_interpolate_states_at_split_curve_node(sav_scss);
	roclists_after_split(rp,contact->curve,curves,YES);
	delete_interior_points_of_curve(front,dncur[4]);
	diffn = Node_of(dncur[6],c_orient[6]);
	Node_vel(diffn)[0] = nod_v[0];
	Node_vel(diffn)[1] = nod_v[1];
	node_type(diffn) = DIFFRACTION_NODE;
	copy_RP_DATA_structure(Rp_data(diffn),RP);
	RP = Rp_data(diffn);

	/* Set states at node of ahead contact */

	copy_state(Left_state_at_node(dncur[6],c_orient[6]),l_st[6]);
	copy_state(Right_state_at_node(dncur[6],c_orient[6]),r_st[6]);
	set_status_at_node(dncur[6],c_orient[6],status[6]);

	/* Install incident shock at new node */

	dncur[0] = shock->curve;
	change_node_of_curve(dncur[0],c_orient[0],diffn);
	cut_curve(pc,crossb_sh,dncur[0],c_orient[0],front,l_st[0],r_st[0]);
	set_status_at_node(dncur[0],c_orient[0],status[0]);

#if defined(DEBUG_SHOCK_CONTACT)
	    if (debugging("shock_contact"))
	    {
	    	(void) printf("Interface before "
	    	              "delete_null_boundary_curves()\n");
	    	print_interface(rp->new_intfc);
	    }
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	    delete_null_boundary_curves(rp,front,(POINTER)wave);
#if defined(DEBUG_SHOCK_CONTACT)
	    if (debugging("shock_contact"))
	    {
	    	(void) printf("Interface after "
	    	              "delete_null_boundary_curves()\n");
	    	print_interface(rp->new_intfc);
	    }
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	/*
	*  Make reflected and transmitted curves, set status at node
	*  and assign states at new diffraction node 
	*/

	oppn = Node_of(dncur[4],Opposite_orient(c_orient[4]));
	for (j = 1; j < 6; j++)
	{
	    if (!track[j])
		continue;
	    set_curve_at_shock_diffraction(&dncur[j],c_orient[j],l_comp[j],
					   r_comp[j],l_st[j],r_st[j],w_type[j],
					   status[j],diffn,oppn);
	}


	interpolate_intfc_states(intfc) = YES;
	if (!install_dfrctn_to_bdry(dncur,c_orient,i_to_f_dir,bdry_type,
				       newcomp,rp,RP,front))
	{
	    set_current_interface(save_intfc);
	    interpolate_intfc_states(intfc) = sav_interp;
	    rp->dt_frac = min(rp->dt_frac,0.5);
	    (void) printf("WARNING in regular_shock_contact_rp() "
	                  "install to boundary failed\n");
	    debug_print("shock_contact",
	    	  "Left regular_shock_contact_rp(), ans = MODIFY_TIME_STEP\n");
	    return MODIFY_TIME_STEP;
	}

	/* Install curves correctly onto boundaries */
	
	if (Apply_CFL_at_nodes(front) == YES)
	{
	    RP_NODE *rpn;
	    double sep;
	    double max_sep = Max_new_node_separation(front);
	    double *h = front->rect_grid->h;
	    int dim = front->rect_grid->dim;

	    for (rpn = rp->first_rp_node; rpn != NULL; rpn = rpn->next)
	    {
	    	if (rpn->old_node == NULL)
		    continue;
	    	sep = scaled_separation(diffn->posn,rpn->old_node->posn,h,dim);
	    	if (sep > max_sep)
	    		return MODIFY_TIME_STEP;
	    }
	}
	switch (bdry_type) 
	{
	case DIRICHLET_BOUNDARY:
	    break;

	case NEUMANN_BOUNDARY:
	    if (!n_install_dfrctn_to_bdry(dncur,c_orient,rp,front,oldsh,
					     oldct,i_to_f_dir))
	    {
	    	set_current_interface(save_intfc);
	    	interpolate_intfc_states(intfc) = sav_interp;
	    	(void) printf("WARNING in regular_shock_contact_rp(), "
	    	              "n_install_dfrctn_to_bdry() failed\n");
	    	debug_print("shock_contact",
	    	      "Left regular_shock_contact_rp(), ans = ERROR_IN_STEP\n");
	    	return ERROR_IN_STEP;
	    }
	    break;

	default:
	    set_current_interface(save_intfc);
	    interpolate_intfc_states(intfc) = sav_interp;
	    (void) printf("WARNING in regular_shock_contact_rp(), "
	                  "Unknown boundary type\n");
	    debug_print("shock_contact",
	    	  "Left regular_shock_contact_rp(), ans = ERROR_IN_STEP\n");
	    return ERROR_IN_STEP;
	}
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("Interface at end of regular_shock_contact_rp()\n");
	    print_interface(intfc);
	    print_correspond_hyper_surf_list(intfc);
	    if (debugging("states"))
	    	show_intfc_states(intfc);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	set_current_interface(save_intfc);
	interpolate_intfc_states(intfc) = sav_interp;

	/* free storage */


	debug_print("shock_contact","Left regular_shock_contact_rp(), "
			      "ans = GOOD_STEP\n");
	return GOOD_STEP;
}		/*end regular_shock_contact_rp*/

EXPORT	void set_orients_and_w_type_for_shock_diffraction(
	int		w_type_c0,
	ORIENTATION	c0_orient,
	int		w_type_c6,
	ORIENTATION	c6_orient,
	ORIENTATION	*c_orient,
	int		*w_type)
{
	c_orient[0] = c0_orient;
	c_orient[1] = Opposite_orient(c0_orient);
	c_orient[2] = Opposite_orient(c0_orient);
	c_orient[3] = Opposite_orient(c0_orient);
	c_orient[4] = Opposite_orient(c6_orient);
	c_orient[5] = Opposite_orient(c0_orient);
	c_orient[6] = c6_orient;

	w_type[0] = w_type[5] = w_type_c0;
	w_type[4] = w_type[6] = w_type_c6;

	if (is_forward_wave(w_type_c0))
	{
	    w_type[1] = BACKWARD_SOUND_WAVE_LE;
	    w_type[2] = BACKWARD_SHOCK_WAVE;
	    w_type[3] = BACKWARD_SOUND_WAVE_TE;
	}
	else
	{
	    w_type[1] = FORWARD_SOUND_WAVE_LE;
	    w_type[2] = FORWARD_SHOCK_WAVE;
	    w_type[3] = FORWARD_SOUND_WAVE_TE;
	}
}		/*end set_orients_and_w_type_for_shock_diffraction*/

EXPORT	void find_node_vel_at_rp(
	POINT		*pc,
	double		tcr1,
	BOND		*crossb1,
	O_CURVE		*oldc1,
	O_CURVE		*newc1,
	double		*t1,
	double		tcr2,
	BOND		*crossb2,
	O_CURVE		*oldc2,
	O_CURVE		*newc2,
	double		*t2,
	ANGLE_DIRECTION	c2_to_c1_dir,
	int		n_type,
	Front		*front,
	Wave		*wave,
	double		dt,
	double		*node_v)
{
	double		v1[MAXD], v2[MAXD];
	double		vm[MAXD], dv[MAXD];
	double		n1[MAXD], n2[MAXD];
	double		den, num1, num2;
	int		dim = front->interf->dim;
	int		i;
	Locstate	st1a, st1b, st2a, st2b;
	static	Locstate s1l = NULL, s2l = NULL, ansl = NULL, s1r = NULL,
	                 s2r = NULL, ansr = NULL, st0 = NULL;

	debug_print("node_vel","Entered find_node_vel_at_rp()\n");
	if (s1l == NULL)
	{
	    alloc_state(front->interf,&s1l,front->sizest);
	    alloc_state(front->interf,&s1r,front->sizest);
	    alloc_state(front->interf,&s2l,front->sizest);
	    alloc_state(front->interf,&s2r,front->sizest);
	    alloc_state(front->interf,&ansl,front->sizest);
	    alloc_state(front->interf,&ansr,front->sizest);
	    alloc_state(front->interf,&st0,front->sizest);
	}

	left_state_along_bond(tcr1,crossb1,newc1->curve,s1l);
	right_state_along_bond(tcr1,crossb1,newc1->curve,s1r);
	find_tangent_to_propagated_curve(pc,crossb1,oldc1,newc1,t1,front,
					 (POINTER)wave,dt);

	left_state_along_bond(tcr2,crossb2,newc2->curve,s2l);
	right_state_along_bond(tcr2,crossb2,newc2->curve,s2r);
	find_tangent_to_propagated_curve(pc,crossb2,oldc2,newc2,t2,front,
					 (POINTER)wave,dt);

	find_states_at_node_of_interaction(newc1->orient,newc2->orient,
					   c2_to_c1_dir,s1l,s1r,s2l,s2r,&st1a,
					   &st1b,&st2a,&st2b);
	interpolate_states(front,0.5,0.5,
		           Coords(Node_of_o_curve(newc1)->posn),st1a,
		           Coords(Node_of_o_curve(newc2)->posn),st2a,st0);
	if (compute_node_velocity(wave_type(newc1->curve),
				  wave_type(newc2->curve),st0,st1b,st2b,
				  t1,t2,node_v,dim,n_type,c2_to_c1_dir))
	    return;

	if (newc1->orient == POSITIVE_ORIENTATION)
	{
	    n1[0] =  t1[1];
	    n1[1] = -t1[0];
	}
	else
	{
	    n1[0] = -t1[1];
	    n1[1] =  t1[0];
	}
	w_speed(Coords(pc),s1l,s1r,ansl,ansr,v1,0.0,n1,
		wave_type(newc1->curve),front);

	if (newc2->orient == POSITIVE_ORIENTATION)
	{
	    n2[0] =  t2[1];
	    n2[1] = -t2[0];
	}
	else
	{
	    n2[0] = -t2[1];
	    n2[1] =  t2[0];
	}
	w_speed(Coords(pc),s2l,s2r,ansl,ansr,v2,0.0,n2,
		wave_type(newc2->curve),front);

	if (debugging("node_vel"))
	{
	    (void) printf("t1 = %g, %g, n1 = %g, %g\n",t1[0],t1[1],n1[0],n1[1]);
	    (void) printf("\tv1 = %g, %g\n",v1[0],v1[1]);
	    (void) printf("t2 = %g, %g, n2 = %g, %g\n",t2[0],t2[1],n2[0],n2[1]);
	    (void) printf("\tv2 = %g, %g\n",v2[0],v2[1]);
	}
	(void) vector_product(t2,t1,&den,dim);
	if (fabs(den) < EPSILON)
	{
	    screen("ERROR in find_node_vel_at_rp(), "
	           "Propagated bonds are parallel\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < dim; i++)
	{
	    vm[i] = 0.5*(v1[i] + v2[i]);
	    dv[i] = 0.5*(v1[i] - v2[i]);
	}
	(void) vector_product(dv,t2,&num1,dim);
	(void) vector_product(dv,t1,&num2,dim);
	for (i = 0; i < dim; i++)
	    node_v[i] = vm[i] + (num1*(t1[i]) + num2*(t2[i]))/den;
	debug_print("node_vel","Left find_node_vel_at_rp()\n");
}		/*end find_node_vel_at_rp*/

LOCAL	int angles_consistent_at_diff_node(
	RP_DATA		*RP,
	int		is_reflected_shock)
{
	double		b[5];
	int		num_angs, i;

	debug_print("shock_contact","Entered angles_consistent_at_diff_node()\n");
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("is_reflected_shock = %s, ",
			 (is_reflected_shock) ? "YES" : "NO");
	    print_angle_direction("RP->ang_dir =",RP->ang_dir,"\n");
	    (void) printf("ang[0] = %g\n",RP->ang[0]);
	    if (is_reflected_shock)
	    	(void) printf("ang[2] = %g\n",RP->ang[2]);
	    else
	    {
	    	(void) printf("rarefaction_angle0 = %g\n",RP->ang[1]);
	    	(void) printf("rarefaction_angle1 = %g\n",RP->ang[3]);
	    }
	    (void) printf("ang[4] = %g\n",RP->ang[4]);
	    (void) printf("ang[5] = %g\n",RP->ang[5]);
	    (void) printf("ang[6] = %g\n",RP->ang[6]);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	if (is_reflected_shock)
	{
	    num_angs = 4;
	    b[0] = normalized_angle(RP->ang[2] - RP->ang[0]);
	    b[1] = normalized_angle(RP->ang[4] - RP->ang[0]);
	    b[2] = normalized_angle(RP->ang[5] - RP->ang[0]);
	    b[3] = normalized_angle(RP->ang[6] - RP->ang[0]);
	}
	else
	{
	    num_angs = 5;
	    b[0] = normalized_angle(RP->ang[1] - RP->ang[0]);
	    b[1] = normalized_angle(RP->ang[3] - RP->ang[0]);
	    b[2] = normalized_angle(RP->ang[4] - RP->ang[0]);
	    b[3] = normalized_angle(RP->ang[5] - RP->ang[0]);
	    b[4] = normalized_angle(RP->ang[6] - RP->ang[0]);
	}
	if (RP->ang_dir == CLOCKWISE)
	{
	    for (i = 0; i < num_angs; i++)
	    	b[i] = 2.0 * PI - b[i];
	}
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    for (i = 0; i < num_angs; i++)
	    	(void) printf("b[%d] = %g\n",i,b[i]);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	for (i = 1; i < num_angs; i++)
	{
	    if (b[i] < b[i-1])
		return NO;
	}
	debug_print("shock_contact","Left angles_consistent_at_diff_node()\n");
	return YES;
}		/*end angles_consistent_at_diff_node*/



EXPORT	void assign_ahead_states_at_shock_diffraction(
	double		tcr_sh,
	BOND		*crossb_sh,
	CURVE		*shock,
	ORIENTATION	sh_or,
	double		tcr_ct,
	BOND		*crossb_ct,
	CURVE		*contact,
	ORIENTATION	ct_or,
	COMPONENT	*comp,
	RP_DATA		*RP)
{
	F_USER_INTERFACE *fuh;
	void		(*sav_interpolator)(double,double,double*,Locstate,double*,
					    Locstate,RECT_GRID*,Locstate);
	boolean		(*sav_tri_interpolator)(double,double,double,double*,
						Locstate,double*,Locstate,
						double*,Locstate,RECT_GRID*,
						Locstate);
	INTERFACE	*intfc = contact->interface;

	if (curve_ang_oriented_l_to_r(RP->ang_dir,sh_or))
	{
	    if (crossb_sh != NULL)
	    {
	    	right_state_along_bond(tcr_sh,crossb_sh,shock,RP->state[1]);
	    }
	    else
	    {
	    	copy_state(RP->state[1],Right_state_at_node(shock,sh_or));
	    }
	    if (comp != NULL)
	    {
	    	comp[0] = negative_component(shock);
	    	comp[1] = positive_component(shock);
	    }
	}
	else 
	{
	    if (crossb_sh != NULL)
	    {
	    	left_state_along_bond(tcr_sh,crossb_sh,shock,RP->state[1]);
	    }
	    else
	    {
	    	copy_state(RP->state[1],Left_state_at_node(shock,sh_or));
	    }
	    if (comp != NULL)
	    {
	    	comp[0] = positive_component(shock);
	    	comp[1] = negative_component(shock);
	    }
	}

	/* Read off states from propagated front contact */

	fuh = &f_user_interface(intfc);
	sav_interpolator = fuh->_bi_interpolate_intfc_states;
	sav_tri_interpolator = fuh->_tri_interpolate_intfc_states;
	fuh->_bi_interpolate_intfc_states = gt_lin_comb_states;
	fuh->_tri_interpolate_intfc_states = gt_tri_lin_comb_states;
	if (curve_ang_oriented_l_to_r(RP->ang_dir,ct_or))
	{
	    if (crossb_ct != NULL)
	    {
	    	right_state_along_bond(tcr_ct,crossb_ct,contact,RP->state[0]);
		left_state_along_bond(tcr_ct,crossb_ct,contact,RP->state[6]);
	    }
	    else
	    {
	    	copy_state(RP->state[0],Right_state_at_node(contact,ct_or));
		copy_state(RP->state[6],Left_state_at_node(contact,ct_or));
	    }
	    if (comp != NULL)
	    	comp[6] = negative_component(contact);
	}
	else 
	{
	    if (crossb_ct != NULL)
	    {
	    	left_state_along_bond(tcr_ct,crossb_ct,contact,RP->state[0]);
		right_state_along_bond(tcr_ct,crossb_ct,contact,RP->state[6]);
	    }
	    else
	    {
	    	copy_state(RP->state[0],Left_state_at_node(contact,ct_or));
		copy_state(RP->state[6],Right_state_at_node(contact,ct_or));
	    }
	    if (comp != NULL)
	    	comp[6] = positive_component(contact);
	}
	fuh->_bi_interpolate_intfc_states = sav_interpolator;
	fuh->_tri_interpolate_intfc_states = sav_tri_interpolator;
}		/*end assign_ahead_states_at_shock_diffraction*/

/*
*		set_track_and_newcomp_list_for_shock_diffraction():
*
*	Assumes that newcomp[0,1,6] are already set.
*/

EXPORT	void set_track_and_newcomp_list_for_shock_diffraction(
	boolean		*track,
	COMPONENT	*newcomp,
	int		*w_type,
	RP_DATA		*RP,
	boolean		is_reflected_shock,
	Front		*fr)
{
	int		i;

	debug_print("shock_contact",
	    "Entered set_track_and_newcomp_list_for_shock_diffraction()\n");
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	    (void) printf("is_reflected_shock = %s\n",
	                  y_or_n(is_reflected_shock));
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	track[0] = track[4] = track[6] = YES;
	for (i = 2; i < 6; i++)
	    newcomp[i] = NO_COMP;

	if (is_shock_wave(w_type[0]) &&
	    track_scattered_wave(DIFFRACTION_NODE,SHOCK_WAVE,TRANSMITTED,
	    		         RP->state[6],RP->state[5],fr))
	{
	    newcomp[5] = new_component(NEW_COMP);
	    track[5] = YES;
	}
	else
	{
	    newcomp[5] = newcomp[6];
	    track[5] = NO;
	}

	if (is_reflected_shock)
	{
	    track[1] = track[3] = NO;
	    if (is_shock_wave(w_type[0]) &&
	        track_scattered_wave(DIFFRACTION_NODE,SHOCK_WAVE,REFLECTED,
	    		             RP->state[1],RP->state[4],fr))
	    {
		track[2] = YES;
		newcomp[4] = new_component(NEW_COMP);
	    }
	    else
	    {
	    	track[2] = NO;
	    	newcomp[4] = newcomp[1];
	    }
	    newcomp[2] = newcomp[1];
	    newcomp[3] = newcomp[4];
	}
	else
	{
	    track[2] = NO;
	    if (is_shock_wave(w_type[0]) &&
	        track_scattered_wave(DIFFRACTION_NODE,RAREF_LEADING_EDGE,
				     REFLECTED,RP->state[1],RP->state[4],fr))
	    {
	    	track[1] = track[3] = YES;
	    	newcomp[3] = newcomp[2] = new_component(NEW_COMP);
	    	newcomp[4] = new_component(NEW_COMP);
	    }
	    else
	    {
	    	track[1] = track[3] = NO;
	    	newcomp[4] = newcomp[3] = newcomp[2] = newcomp[1];
	    }
	}
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    int i;

	    for (i = 0; i < 7; i++)
	    {
	    	(void) printf("track[%d] = %s, newcomp[%d] = %d\n",i,
	    		      (track[i]) ? "YES" : "NO",i,newcomp[i]);
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	debug_print("shock_contact",
	      "Left set_track_and_newcomp_list_for_shock_diffraction()\n");
}		/*end set_track_and_newcomp_list_for_shock_diffraction*/

EXPORT	void set_states_and_comps_about_shock_diffraction(
	COMPONENT	*newcomp,
	ANGLE_DIRECTION	i_to_f_dir,
	ORIENTATION	*c_orient,
	COMPONENT	*l_comp,
	COMPONENT	*r_comp,
	Locstate	*l_st,
	Locstate	*r_st,
	RP_DATA		*RP)
{
	int		j;

	debug_print("shock_contact",
	      "Entered set_states_and_comps_about_shock_diffraction()\n");
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    char s[80];

	    print_angle_direction("i_to_f_dir = ",i_to_f_dir,"\n");
	    for (j = 0; j < 7; j++)
	    {
	    	(void) sprintf(s,"c_orient[%d] = ",j);
	    	print_orientation(s,c_orient[j],"\n");
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	for (j = 0; j < 7; j++) 
	{
	    if (curve_ang_oriented_l_to_r(i_to_f_dir,c_orient[j]))
	    {
	    	l_st[j] = RP->state[(j+1)%7];
	    	l_comp[j] = newcomp[(j+1)%7];
	    	r_st[j] = RP->state[j];
	    	r_comp[j] = newcomp[j];
	    }
	    else 
	    {
	    	l_st[j] = RP->state[j];
	    	l_comp[j] = newcomp[j];
	    	r_st[j] = RP->state[(j+1)%7];
	    	r_comp[j] = newcomp[(j+1)%7];
	    }
	}
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    for (j = 0; j < 7; j++)
	    {
	    	(void) printf("l_comp[%d] = %d, r_comp[%d] = %d\n",
	    		      j,l_comp[j],j,r_comp[j]);
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	debug_print("shock_contact",
	      "Left set_states_and_comps_about_shock_diffraction()\n");
}		/*end set_states_and_comps_about_shock_diffraction*/

LOCAL	void set_curve_at_shock_diffraction(
	CURVE		**curve,
	ORIENTATION	orient,
	COMPONENT	l_comp,
	COMPONENT	r_comp,
	Locstate	l_st,
	Locstate	r_st,
	int		w_type,
	int		n_status,
	NODE		*diffn,
	NODE		*oppn)
{
	ORIENTATION	opp_orient = Opposite_orient(orient);

	debug_print("shock_contact","Entered set_curve_at_shock_diffraction()\n");
	if (*curve == NULL)
	{
	    if (orient == POSITIVE_ORIENTATION)
	    	*curve = make_curve(l_comp,r_comp,diffn,oppn);
	    else
	    	*curve = make_curve(l_comp,r_comp,oppn,diffn);
	}
	else
	{
	    negative_component((*curve)) = l_comp;
	    positive_component((*curve)) = r_comp;
	}
	wave_type(*curve) = w_type;
	copy_state(Left_state_at_node(*curve,orient),l_st);
	copy_state(Right_state_at_node(*curve,orient),r_st);

	copy_state(Left_state_at_node(*curve,opp_orient),l_st);
	copy_state(Right_state_at_node(*curve,opp_orient),r_st);

	/* Set node status */

	set_status_at_node(*curve,orient,n_status);
	set_status_at_node(*curve,opp_orient,INCIDENT);

	debug_print("shock_contact","Left set_curve_at_shock_diffraction()\n");
}		/*end set_curve_at_shock_diffraction*/


LOCAL	int n_install_dfrctn_to_bdry(
	CURVE		**dncur,
	ORIENTATION	*c_orient,
	RPROBLEM	*rp,
	Front		*front,
	O_CURVE		*oldsh,
	O_CURVE		*oldct,
	ANGLE_DIRECTION	i_to_f_dir)
{
	INTERFACE *intfc = rp->new_intfc;
	CURVE		*tbc, *rbc;
	POINT		*pt[5];
	Locstate	sr, sl;
	double		coords[MAXD], partial_dt;
	double		t[MAXD], n[MAXD];
	ORIENTATION	copp_or[7];
	ORIENTATION	tbc_orient, rbc_orient;
	int		dim = front->rect_grid->dim;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	int		i;
	static int	w_type[5] = { BACKWARD_SOUND_WAVE_LE,
				      BACKWARD_SHOCK_WAVE,
				      BACKWARD_SOUND_WAVE_TE,
				      CONTACT,
				      FORWARD_SHOCK_WAVE};
	static double	pjump[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};

	debug_print("nbdry","Enter n_install_dfrctn_to_bdry()\n");
	for (i = 1; i < 6; i++)
	    if (dncur[i] != NULL)
		break;
	if (i == 6) /* No tracked curves at boundary */
	{
	    debug_print("nbdry","Left n_install_dfrctn_to_bdry()\n");
	    return YES;
	}
	if (debugging("nbdry"))
	{
	    char mesg[80];

	    print_angle_direction("i_to_f_dir =",i_to_f_dir,"\n");
	    for (i = 0; i < 7; i++)
	    {
	    	if (dncur[i] == NULL)
		    continue;
	    	(void) sprintf(mesg,"c_orient[%d] =",i);
	    	print_orientation(mesg,c_orient[i],"\n");
	    }
	    (void) printf("Interface before n_install_dfrctn_to_bdry()\n");
	    print_interface(intfc);
	}

	identify_bdry_curves_at_shock_diffraction(&tbc,&tbc_orient,
						  &rbc,&rbc_orient,
						  i_to_f_dir,rp);
	find_tangent_to_curve(Node_of(tbc,tbc_orient)->posn,
			     Bond_at_node(tbc,tbc_orient),
			     tbc,tbc_orient,t,front);
	normal(Node_of(tbc,tbc_orient)->posn,
	       Hyper_surf_element(Bond_at_node(tbc,tbc_orient)),
	       Hyper_surf(tbc),n,front);
	if (debugging("nbdry"))
	{
	    char ltitle[80], rtitle[80];
	    (void) printf("Boundary nodes before reset\n");
	    for (i = 1; i < 6; i++)
	    {
	    	if (dncur[i] == NULL)
		    continue;
	    	print_node(Node_of(dncur[i],Opposite_orient(c_orient[i])));
	    }
	    (void) printf("Wall states before reset\n");
	    for (i = 1; i < 6; i++)
	    {
	    	if (dncur[i] == NULL)
		    continue;
	    	(void) sprintf(ltitle,"Left state dncur[%d]",i);
	    	(void) sprintf(rtitle,"Right state dncur[%d]",i);
	    	sl = Left_state_at_node(dncur[i],Opposite_orient(c_orient[i]));
	    	sr = Right_state_at_node(dncur[i],Opposite_orient(c_orient[i]));
	    	if (curve_ang_oriented_l_to_r(i_to_f_dir,c_orient[i]))
	    	{
	    	    verbose_print_state(rtitle,sr);
	    	    verbose_print_state(ltitle,sl);
	    	}
	    	else
	    	{
	    	    verbose_print_state(ltitle,sl);
	    	    verbose_print_state(rtitle,sr);
	    	}
	    }
	    (void) printf("Transmitted wall tangent = <%g, %g>\n",
	    	      t[0],t[1]);
	    (void) printf("Transmitted wall normal  = <%g, %g>\n",
	    	      n[0],n[1]);
	    (void) printf("Transmitted side boundary, ");
	    print_orientation("orient =",tbc_orient,"\n");
	    print_curve(tbc);
	    (void) printf("Reflected side boundary, ");
	    print_orientation("orient =",rbc_orient,"\n");
	    print_curve(rbc);
	}
	if (is_obstacle_state(Right_state_at_node(tbc,tbc_orient)))
	{
	    for (i = 0; i < dim; i++)
		n[i] *= -1.0;
	}
	sr = (curve_ang_oriented_l_to_r(i_to_f_dir,oldct->orient)) ?
		Right_state_at_node_of_o_curve(oldct) :
		Left_state_at_node_of_o_curve(oldct);
	sl = (curve_ang_oriented_l_to_r(i_to_f_dir,oldsh->orient)) ?
		Left_state_at_node_of_o_curve(oldsh) :
		Right_state_at_node_of_o_curve(oldsh);
	pjump[3] = set_pjump_at_wave(Node_of_o_curve(oldct)->posn,
			Hyper_surf_element(Bond_at_node_of_o_curve(oldct)),
			Hyper_surf(oldct->curve),front,NULL);

	for (i = 0; i < 7; i++)
	    copp_or[i] = Opposite_orient(c_orient[i]);
	interpolate_intfc_states(intfc) = YES;

	partial_dt = find_position_and_dt_of_intersection(oldsh->curve,
				                          oldsh->orient,
							  oldct->curve,
							  oldct->orient,
							  coords,front,rp->dt);

	for (i = 0; i < 5; i++)
	{
	    if (debugging("nbdry"))
	    {
	    	if (dncur[i+1] != NULL)
	    	    (void) printf("Attaching curve %d\n",i+1);
	    }
	    attach_phys_curve_to_n_bdry(dncur[i+1],copp_or[i+1],
				        w_type[i],pt+i,sl,sr,coords,
				        t,n,front,pjump[i],partial_dt,front);
	}
	if (debugging("nbdry"))
	{
	    char ltitle[80], rtitle[80];

	    verbose_print_state("left state into w_speed",sl);
	    verbose_print_state("right state into w_speed",sr);
	    (void) printf("New boundary node positions\n");
	    for (i = 0; i < 5; i++)
	    {
	    	if (pt[i] == NULL)
		    continue;
	    	(void) printf("\tpt[%d] = <%g, %g>\n",i,
	    		      Coords(pt[i])[0],Coords(pt[i])[1]);
	    }
	    (void) printf("Wall states from w_speed\n");
	    for (i = 1; i < 6; i++)
	    {
	    	if (dncur[i] == NULL)
		    continue;
	    	(void) sprintf(ltitle,"Left state pt[%d]",i-1);
	    	(void) sprintf(rtitle,"Right state pt[%d]",i-1);
	    	sl = left_state(pt[i-1]);
	    	sr = right_state(pt[i-1]);
	    	verbose_print_state(ltitle,sl);
	    	verbose_print_state(rtitle,sr);
	    }
	}

	if (!realign_phys_curves_at_nbdry(dncur,copp_or,rbc,rbc_orient,tbc,
					     tbc_orient,i_to_f_dir,pt,front))
	{
	    (void) printf("WARNING in n_install_dfrctn_to_bdry(), "
	                  "unable to realign boundary curves\n");
	    debug_print("nbdry","Left n_install_dfrctn_to_bdry()\n");
	    return NO;
	}


	interpolate_intfc_states(intfc) = sav_intrp;
	if (debugging("nbdry"))
	{
	    char ltitle[80], rtitle[80];
	    (void) printf("Wall states after reset\n");
	    for (i = 1; i < 6; i++)
	    {
	    	if (dncur[i] == NULL)
		    continue;
	    	(void) sprintf(ltitle,"Left state dncur[%d]",i);
	    	(void) sprintf(rtitle,"Right state dncur[%d]",i);
	    	sl = Left_state_at_node(dncur[i],Opposite_orient(c_orient[i]));
		sr = Right_state_at_node(dncur[i],Opposite_orient(c_orient[i]));
		if (curve_ang_oriented_l_to_r(i_to_f_dir,c_orient[i]))
		{
		    verbose_print_state(rtitle,sr);
		    verbose_print_state(ltitle,sl);
		}
		else
		{
		    verbose_print_state(ltitle,sl);
		    verbose_print_state(rtitle,sr);
		}
	    }
	    (void) printf("Interface after n_install_dfrctn_to_bdry()\n");
	    print_interface(intfc);
	    show_intfc_states(intfc);
	}
	debug_print("nbdry","Left n_install_dfrctn_to_bdry()\n");
	return YES;
}		/*end n_install_dfrctn_to_bdry*/


EXPORT	int realign_phys_curves_at_nbdry(
	CURVE		**dncur,
	ORIENTATION	*c_orient,
	CURVE		*rbc,
	ORIENTATION	rbc_orient,
	CURVE		*tbc,
	ORIENTATION	tbc_orient,
	ANGLE_DIRECTION	i_to_f_dir,
	POINT		**pt,
	Front		*front)
{
	CURVE		*c, *next_bc;
	ORIENTATION	orient;
	ORIENTATION	next_bc_orient;
	int		i;


	/* Realign reflected boundary */
	for (i = 0; i < 5; i++)
	    if (pt[i] != NULL)
		break;
	reset_wall_node(pt[i],left_state(pt[i]),right_state(pt[i]),NULL,
			rbc,rbc_orient,dncur[i+1],c_orient[i+1],
			i_to_f_dir,front);

	/* Realign transmitted boundary */
	for (i = 4; i >= 0; i--)
	    if (pt[i] != NULL)
		break;
	reset_wall_node(pt[i],right_state(pt[i]),left_state(pt[i]),NULL,
			tbc,tbc_orient,dncur[i+1],c_orient[i+1],
			i_to_f_dir,front);

	/* Realign intermediate boundaries */
	for (i = 0; i < 5; i++)
	    if (pt[i] != NULL)
		break;
	Check_return(next_boundary(rbc,rbc_orient,&c,&orient),
		     realign_phys_curves_at_nbdry)
	orient = Opposite_orient(orient);
	for (i++; c != tbc && i < 5; i++)
	{
	    Check_return(next_boundary(c,orient,&next_bc,&next_bc_orient),
	    	     realign_phys_curves_at_nbdry)
	    delete_interior_points_of_curve(front,c);
	    if (next_bc == tbc)
		break;
	    while (i < 5 && pt[i] == NULL) i++;
	    if (i == 5)
	    {
	        (void) printf("WARNING in realign_phys_curves_at_nbdry(), "
	                      "inconsistent boundary curves\n");
	        return NO;
	    }
	    reset_wall_node(pt[i],left_state(pt[i]),right_state(pt[i]),
			    c->first,c,orient,dncur[i+1],c_orient[i+1],
			    i_to_f_dir,front);
	    orient = Opposite_orient(next_bc_orient);
	    c = next_bc;
	}

	return YES;
}		/*end realign_phys_curves_at_nbdry*/

EXPORT	void interpolate_state_next_to_node(
	CURVE		*c,
	ORIENTATION	c_orient,
	Front		*front)
{
	BOND		*b, *bf;
	double		alpha, beta;
	double		*crds0, *crds2;
	Locstate	s0, s1, s2;
	ORIENTATION	opp_orient;
	
	if (c == NULL)
	    return;
	opp_orient = Opposite_orient(c_orient);

	b = Bond_at_node(c,c_orient);
	bf = Following_bond(b,c_orient);
	if (bf == NULL)
	    return;
	alpha = bond_length(bf)/(bond_length(b) + bond_length(bf));
	beta = bond_length(b)/(bond_length(b) + bond_length(bf));

	crds0 = Coords(Node_of(c,c_orient)->posn);
	crds2 = Coords(Point_of_bond(bf,opp_orient));
	s0 = Left_state_at_node(c,c_orient);
	s1 = left_state_at_point_on_curve(Point_of_bond(bf,c_orient),bf,c);
	s2 = left_state_at_point_on_curve(Point_of_bond(bf,opp_orient),bf,c);
	interpolate_states(front,alpha,beta,crds0,s0,crds2,s2,s1);

	s0 = Right_state_at_node(c,c_orient);
	s1 = right_state_at_point_on_curve(Point_of_bond(bf,c_orient),bf,c);
	s2 = right_state_at_point_on_curve(Point_of_bond(bf,opp_orient),bf,c);
	interpolate_states(front,alpha,beta,crds0,s0,crds2,s2,s1);
}		/*end interpolate_state_next_to_node*/
	
EXPORT void attach_phys_curve_to_n_bdry(
	CURVE		*c,
	ORIENTATION	c_orient,
	int		w_type,
	POINT		**pt,
	Locstate	sl,
	Locstate	sr,
	double		*avg_coords,
	double		*t,
	double		*n,
	Front		*front,
	double		pjump,
	double		dt,
	Front		*fr)
{
	POINT		*p;
	BOND		*b;
	Locstate	stmp, ansl, ansr;
	double		W[MAXD];
	double		sep_max, *h = front->rect_grid->h;
	boolean		sav_intrp;
	int		dim = front->rect_grid->dim;
	ORIENTATION	opp_or = Opposite_orient(c_orient);
	int		j;
	static Locstate templ = NULL, tempr = NULL, obst = NULL;

	debug_print("nbdry","Entered attach_phys_curve_to_n_bdry()\n");

	if (templ == NULL)
	{
	    alloc_state(front->interf,&templ,front->sizest);
	    alloc_state(front->interf,&tempr,front->sizest);
	    obst = return_obst_state();
	}
	    	/*TOLERANCE*/
	sep_max = sqrt(front->Redist.spacing[0] * front->Redist.spacing[1]);

	if (c == NULL)
	{
	    *pt = NULL;
	    debug_print("nbdry","Left attach_phys_curve_to_n_bdry()\n");
	    return;
	}
	    /* Impose boundary conditions on sl and sr */
	*pt = Point(avg_coords);
	ansr = right_state(*pt);
	ansl = left_state(*pt);
	w_speed(avg_coords,obst,sl,ansl,templ,W,0.0,n,NEUMANN_BOUNDARY,front);
	w_speed(avg_coords,obst,sr,ansl,tempr,W,0.0,n,NEUMANN_BOUNDARY,front);
	w_speed(avg_coords,templ,tempr,ansl,ansr,W,pjump,t,w_type,front);
	for (j = 0; j < dim; j++)
	    Coords(*pt)[j] += W[j]*dt;

	if (debugging("nbdry"))
	{
	    (void) printf("Calculated wall node position = (%g, %g)\n",
	    	          Coords(*pt)[0],Coords(*pt)[1]);
	    verbose_print_state("ansr",ansr);
	    verbose_print_state("ansl",ansl);
	}

	/* Project onto boundary */
	if (is_bdry(Node_of(c,c_orient)))
	    nearest_boundary_point(Coords(*pt),Coords(*pt),front->rect_grid);

	p = Point(Coords(*pt));

	if ((!intersect_ray_with_curve(*pt,n,NULL,NULL,c,opp_or,&b,p)) ||
	    (scaled_separation(*pt,p,h,dim) > sep_max))
	{
	    (void) printf("WARNING in attach_phys_curve_to_n_bdry(), ");
	    (void) printf("can't install new wave normal to wall\n");
	    debug_print("nbdry","Left attach_phys_curve_to_n_bdry()\n");
	    return;
	}
	sav_intrp = interpolate_intfc_states(c->interface);
	interpolate_intfc_states(c->interface) = NO;
	copy_state(left_state(p),left_state(*pt));
	copy_state(right_state(p),right_state(*pt));
	stmp = left_state_at_point_on_curve(b->start,b,c);
	if (is_obstacle_state(stmp))
	    copy_state(stmp,left_state(*pt));
	stmp = right_state_at_point_on_curve(b->start,b,c);
	if (is_obstacle_state(stmp))
	    copy_state(stmp,right_state(*pt));
	stmp = left_state_at_point_on_curve(b->end,b,c);
	if (is_obstacle_state(stmp))
	    copy_state(stmp,left_state(*pt));
	stmp = right_state_at_point_on_curve(b->end,b,c);
	if (is_obstacle_state(stmp))
	    copy_state(stmp,right_state(*pt));
	if (insert_point_in_bond(p,b,c) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in attach_phys_curve_to_n_bdry(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	interpolate_intfc_states(c->interface) = sav_intrp;
	while (p != Point_adjacent_to_node(c,c_orient))
	    (void) delete_point_adjacent_to_node(fr,c,c_orient);
	debug_print("nbdry","Left attach_phys_curve_to_n_bdry()\n");
	return;
}		/*end attach_phys_curve_to_n_bdry*/

LOCAL	void	reset_wall_node(
	POINT		*newp,
	Locstate	st0,
	Locstate	st1,
	BOND		*b,
	CURVE		*cb0,
	ORIENTATION	cb0_orient,
	CURVE		*cp,
	ORIENTATION	cp_orient,
	ANGLE_DIRECTION	i_to_f_dir,
	Front		*front)
{
	POINT		*p;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF	*hs;
	CURVE		*cb1;
	ORIENTATION	cb1_orient;
	int		i, dim = front->rect_grid->dim;
	double		s;
	double		coords[MAXD];

	debug_print("nbdry","Entered reset_wall_node()\n");

	Check_return(next_boundary(cb0,cb0_orient,&cb1,&cb1_orient),
		     reset_wall_node)

	if (debugging("nbdry"))
	{
	    verbose_print_state("st0",st0);
	    verbose_print_state("st1",st1);
	    (void) printf("newp = <%g, %g>\n",Coords(newp)[0],Coords(newp)[1]);
	    (void) printf("cb0, ");
	    print_orientation("cb0_orient =",cb0_orient,"\n");
	    print_curve(cb0);
	    (void) printf("cb1, ");
	    print_orientation("cb1_orient =",cb1_orient,"\n");
	    print_curve(cb1);
	}

	if (b != NULL)
	{
	    p = Point_of_bond(b,Opposite_orient(cb0_orient));
	}
	else
	{
	    if (long_nearest_interface_point(Coords(newp),NO_COMP,
					     cb0->interface,INCLUDE_BOUNDARIES,
					     Hyper_surf(cb0),coords,&s,
					     &hse,&hs) != YES)
	    {
		screen("ERROR in reset_wall_node(), "
		       "long_nearest_interface_point() failed\n");
		clean_up(ERROR);
	    }
	    p = Point_of_bond(Bond_of_hse(hse),Opposite_orient(cb0_orient));
	}
	while (p != Point_adjacent_to_node(cb0,cb0_orient))
	    (void) delete_point_adjacent_to_node(front,cb0,cb0_orient);

	if (is_obstacle_state(Right_state_at_node(cb0,cb0_orient)))
	{
	    copy_state(Left_state_at_node(cb0,cb0_orient),st0);
	}
	else
	{
	    copy_state(Right_state_at_node(cb0,cb0_orient),st0);
	}
	if (is_obstacle_state(Right_state_at_node(cb1,cb1_orient)))
	{
	    copy_state(Left_state_at_node(cb1,cb1_orient),st1);
	}
	else
	{
	    copy_state(Right_state_at_node(cb1,cb1_orient),st1);
	}
	for (i = 0; i < dim; i++)
	    Coords(Node_of(cb0,cb0_orient)->posn)[i] = Coords(newp)[i];
	set_bond_length(Bond_at_node(cb0,cb0_orient),dim);

	if (curve_ang_oriented_l_to_r(i_to_f_dir,cp_orient))
	{
	    copy_state(Left_state_at_node(cp,cp_orient),left_state(newp));
	    copy_state(Right_state_at_node(cp,cp_orient),right_state(newp));
	}
	else
	{
	    copy_state(Left_state_at_node(cp,cp_orient),right_state(newp));
	    copy_state(Right_state_at_node(cp,cp_orient),left_state(newp));
	}
	debug_print("nbdry","Left reset_wall_node()\n");
}		/*end reset_wall_node*/
#endif /* defined(TWOD) && defined(FULL_PHYSICS) */
