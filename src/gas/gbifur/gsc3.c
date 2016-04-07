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
*				gsc3.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Solves two dimensional shock contact Riemann problems.
*
*/

#if defined(FULL_PHYSICS) && defined(TWOD)

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	int	identify_curves_after_n_bdry_cross(CURVE**,ORIENTATION*,
						   CURVE**,ORIENTATION*,
						   CURVE**,ORIENTATION*,
						   RPROBLEM*,ANGLE_DIRECTION);
LOCAL	int	leaving_shock_contact(const char*,int);
LOCAL	int	new_cross_bifurcate(NODE*,O_CURVE*,O_CURVE*,Front*,RPROBLEM*);
LOCAL	int	te_overtakes_inc_shock(O_CURVE**,O_CURVE**,O_CURVE**,
				       O_CURVE**,Front*,Wave*,RPROBLEM*);



EXPORT int diffracted_shock_exits_contact(
	Front       *front,
	Front       *newfront,
	Wave        *wave,
	RPROBLEM    *rp)
{
	NODE		*new_node;
	CURVE		*joined_c;
	O_CURVE		*oc, *oldoc;
	O_CURVE		*newc[7][2], *oldc[7][2];
	RP_NODE		*rpn[2];
	COMPONENT	comp[7][2];
	int		status;
	int		w_type, n_status;
	ANGLE_DIRECTION	i_to_f_dir[2];
	int		i, j, k;
	int		jmax, jmin;
	int		closed_curve;
	static	int	add_newn[7] = { NO, NO, YES, NO, NO, YES, NO};

	debug_print("shock_contact","Entering diffracted_shock_exits_contact()\n");
#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    if (!debugging("2drp"))
		print_rproblem(rp);
	    (void) printf("Old interface\n");
	    print_interface(rp->old_intfc);
	    (void) printf("New interface\n");
	    print_interface(rp->new_intfc);
	    print_correspond_hyper_surf_list(rp->old_intfc);
	    print_correspond_hyper_surf_list(rp->new_intfc);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	rpn[0] = rp->first_rp_node;	rpn[1] = rp->last_rp_node;

	for (i = 0; i < 7; ++i) 
	{
	    newc[i][0] = newc[i][1] = NULL;
	    oldc[i][0] = newc[i][1] = NULL;
	}

	for (oc = rp->ang_ordered_curves->first,
	     oldoc = rp->old_ang_ordered_curves->first;
	    oc && oldoc; oc = oc->next, oldoc = oldoc->next)
	{
	    j = (Node_of_o_curve(oldoc) == rpn[0]->old_node) ?  0 : 1;
	    n_status = status_at_node(oldoc->curve,oldoc->orient);
	    switch (n_status)
	    {
	    case INCIDENT:
#if defined(DEBUG_SHOCK_CONTACT)
	    	if (debugging("shock_contact"))
	    	    (void) printf("Incident curve found\n");
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	    	i = 0;
	    	break;
	    case REFLECTED:
#if defined(DEBUG_SHOCK_CONTACT)
	    	if (debugging("shock_contact"))
	    	    (void) printf("Reflected curve found\n");
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	    	w_type = wave_type(oldoc->curve);
	    	i = (is_rarefaction_leading_edge(w_type)) ? 1 :
	    	    	(is_shock_wave(w_type)) ?  2 : 3;
	    	break;
	    case SLIP:
#if defined(DEBUG_SHOCK_CONTACT)
	    	if (debugging("shock_contact"))
	    	    (void) printf("Contact curve found\n");
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	    	i = 4;
	    	break;
	    case TRANSMITTED:
#if defined(DEBUG_SHOCK_CONTACT)
	    	if (debugging("shock_contact"))
	    	    (void) printf("Transmitted curve found\n");
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	    	i = 5;
	    	break;
	    default:
	    	continue;
	    }
	    newc[i][j] = oc;	oldc[i][j] = oldoc;
	    if (oc == rp->ang_ordered_curves->last ||
	    	    oldoc == rp->old_ang_ordered_curves->last)
	    	break;
	}

	/* Evaluate local components and angle orientations */

	for (j = 0; j < 2; ++j)
	{
	    CURVE	*cinc;
	    ORIENTATION	c_or;
	    int	        l_to_r;

	    find_curve_with_status(rpn[j]->node,&cinc,&c_or,INCIDENT);
	    i_to_f_dir[j] = incident_shock_orientation(wave_type(cinc),c_or);

	    l_to_r = curve_ang_oriented_l_to_r(i_to_f_dir[j],c_or);
	    comp[0][j] = (l_to_r) ? positive_component(cinc) :
	    	    	            negative_component(cinc);
	    comp[1][j] = (l_to_r) ? negative_component(cinc) :
	    			    positive_component(cinc);
	}

	for (i = 1; i < 6; ++i)
	{
	    if (newc[i][0] == NULL && newc[i][1] == NULL)
	    {
	    	comp[i+1][0] = comp[i][0]; 
	    	comp[i+1][1] = comp[i][1]; 
	    	continue;
	    }
	    if (newc[i][0] == NULL || newc[i][1] == NULL)
	    {
	    	screen("ERROR in diffracted_shock_exits_contact(), "
	    	       "UNEXPECTED CASE, CODE NEEDED\n");
	    	return leaving_shock_contact("diffracted_shock_exits_contact",
					     ERROR_IN_STEP);
	    }
	    for (j = 0; j < 2; ++j)
	    {
	    	comp[i+1][j] = (curve_ang_oriented_l_to_r(i_to_f_dir[j],
	    					          newc[i][j]->orient))
		              ? negative_component(newc[i][j]->curve) : 
				positive_component(newc[i][j]->curve);
	    }
	}

#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("COMPONENTS before reset\n");
	    for (i = 0; i < 7; ++i)
	    {
	    	(void) printf("comp[%d][0] = %d, comp[%d][1] = %d\n",
	    		      i,comp[i][0],i,comp[i][1]);
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	for (i = 0; i < 7; ++i)
	    set_equivalent_comps(comp[i][0],comp[i][1],newfront->interf);

	/* Ensure consistent components of new loops */

	for (i = 1; i < 6; ++i)
	{
	    if (newc[i-1][0] == NULL)
		continue;
	    if (comp[i][0] != comp[i][1])
	    {
	    	jmin = (comp[i][0] < comp[i][1]) ? 0 : 1;
	    	jmax = (jmin+1)%2;
	    	reset_component_of_loop(newc[i-1][jmax]->curve,
				       newc[i-1][jmax]->orient,i_to_f_dir[jmax],
				       comp[i][jmin],front);
	    	for (k = i; k < 6; ++k)
	    	{
	    	    comp[k][jmax] = comp[k][jmin];
	    	    if (newc[k][jmin] != NULL)
			break;
	    	}
	    }
	}

#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("COMPONENTS after reset\n");
	    for (i = 0; i < 7; ++i)
	    {
	    	(void) printf("comp[%d][0] = %d, comp[%d][1] = %d\n",
	    	              i,comp[i][0],i,comp[i][1]);
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	/* Delete null curves */

	delete_null_physical_curves(rp);

	/* Join interacting curves */

	if (newc[0][0] && newc[3][0])
	{
	    status = te_overtakes_inc_shock(newc[0],oldc[0],
				            newc[3],oldc[3],front,wave,rp);

	    if (status != GOOD_STEP)
	    {
	    	(void) printf("WARNING: in "
	    	              "diffracted_shock_exits_contact()\n"
	    	              "te_overtakes_inc_shock() failed for "
	    	              "curves %d\n",i);
	    	return leaving_shock_contact("diffracted_shock_exits_contact",
					     status);
	    }

	}

	for (i = 0; i < 6; ++i)
	{
	    if (newc[i][0] == NULL) continue;

	/* Install new curve configurations */

	    closed_curve = (newc[i][0]->curve == newc[i][1]->curve) ? YES : NO;

#if defined(DEBUG_SHOCK_CONTACT)
	    if (debugging("shock_contact"))
	    {
	        (void) printf("processing curves %d, closed_curve = %s\n",
	    		      i,(closed_curve) ? "YES" : "NO");
	    }
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	    status = join_propagated_curves(&new_node,&joined_c,newc[i],oldc[i],
				            ((add_newn[i] || closed_curve) ?
					     YES : NO),
				            front,newfront,(POINTER)wave,rp);

	    if (status != GOOD_STEP)
	    {
	    	(void) printf("WARNING: in "
	    	              "diffracted_shock_exits_contact()\n"
			      "join_propagated_curves() failed for "
		              "curves %d\n",i);
	    	return leaving_shock_contact("diffracted_shock_exits_contact",
					     status);
	    }

	    /* Now check for possible bifurctions in transmitted shocks*/

	    if (closed_curve)
	    {
	    	CURVE *c = newc[i][0]->curve;

	    	set_incident_status(c,POSITIVE_ORIENTATION);
	    	set_incident_status(c,NEGATIVE_ORIENTATION);
	    }

	    if (!add_newn[i])
	    {
	    	/* Do not allow redistribute to create a tangle. */

	    	do_not_redistribute(joined_c) = YES;
	    	continue;
	    }

	    if (!new_cross_bifurcate(new_node,newc[i][0],newc[i][1],
					front,rp))
	    {
	       (void) printf("ERROR in diffracted_shock_exits_contact(), "
	                     "new_cross_bifurcate() failed\n");
	       return leaving_shock_contact("diffracted_shock_exits_contact",
					    ERROR_IN_STEP);
	    }
	}

	(void) delete_node(rpn[0]->node);	rpn[0]->node = NULL;
	(void) delete_node(rpn[1]->node);	rpn[1]->node = NULL;

#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("New interface after "
			  "diffracted_shock_exits_contact()\n");
	    print_interface(rp->new_intfc);
    	    if (debugging("states"))
	    {
	        (void) printf("\n\tSTATES AFTER SHOCK EXITS CONTACT\n\n");
	        show_intfc_states(rp->new_intfc);
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	return leaving_shock_contact("diffracted_shock_exits_contact",
		                     GOOD_STEP);
}		/*end diffracted_shock_exits_contact*/

EXPORT int n_diffracted_shock_exits_contact(
	Front		*front,
	Front		*newfront,
	Wave		*wave,
	RPROBLEM	*rp)
{
	ANGLE_DIRECTION	i_to_f_dir;
	CURVE		*dncur[7];
	CURVE		*tbc, *rbc;
	INTERFACE	*intfc = rp->new_intfc;
	Locstate	sl, sr;
	NODE_FLAG	flag;
	O_CURVE		Osh, Oct;
	O_CURVE		*doc[7];
	ORIENTATION	c_orient[7];
	ORIENTATION	tbc_orient, rbc_orient;
	POINT		*pt[5];
	RP_NODE		*rpn;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	double		t[MAXD], n[MAXD];
	double		intsct_coords[MAXD], partial_dt;
	int		dim = front->rect_grid->dim;
	int		status;
	int		i;
	static double	pjump[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};
	static int	w_type[5] = { BACKWARD_SOUND_WAVE_LE,	
				      BACKWARD_SHOCK_WAVE,
				      BACKWARD_SOUND_WAVE_TE,
				      CONTACT,
				      FORWARD_SHOCK_WAVE};

	debug_print("n_d_shock","Entered n_diffracted_shock_exits_contact()\n");


	for (rpn = rp->first_rp_node; rpn != NULL; rpn = rpn->next)
	{
	    if (node_type(rpn->node) == DIFFRACTION_NODE)
	    	break;
	}
	for (i = 0; i < 7; ++i) scalar(&doc[i],sizeof(O_CURVE));
	if (curves_at_shock_diffraction(rpn->node,doc,&i_to_f_dir,NO) != YES)
	{
	    (void) printf("WARNING in n_diffracted_shock_exits_contact(), "
	                  "Unable to find curves\n");
	    for (i = 0; i < 7; ++i)
		free(doc[i]);
	    debug_print("n_d_shock","Left n_diffracted_shock_exits_contact()\n");
	    return ERROR_IN_STEP;
	}
	for (i = 0; i < 7; ++i)
	{
	    dncur[i] = doc[i]->curve;
	    c_orient[i] = (dncur[i] != NULL) ? doc[i]->orient :
					       ORIENTATION_NOT_SET;
	}
	for (i = 0; i < 7; ++i)
	{
	    free(doc[i]);
	    doc[i] = NULL;
	}
	delete_interior_points_of_curve(front,dncur[0]);
	delete_interior_points_of_curve(front,dncur[6]);
	find_curve_with_status(rpn->old_node,&Osh.curve,&Osh.orient,INCIDENT);
	find_curve_with_status(rpn->old_node,&Oct.curve,&Oct.orient,
			       CONTACT_TARGET);
	i_to_f_dir = incident_shock_orientation(wave_type(Osh.curve),
						Osh.orient);

	partial_dt = find_position_and_dt_of_intersection(Osh.curve,
			                         Opposite_orient(Osh.orient),
						 Oct.curve,
			                         Opposite_orient(Oct.orient),
						 intsct_coords,front,rp->dt);

	if (i_to_f_dir == CLOCKWISE)
	{
	    rbc = rp->bdry_curves->first->curve;
	    rbc_orient = rp->bdry_curves->first->orient;
	    tbc = rp->bdry_curves->last->curve;
	    tbc_orient = rp->bdry_curves->last->orient;
	}
	else
	{
	    tbc = rp->bdry_curves->first->curve;
	    tbc_orient = rp->bdry_curves->first->orient;
	    rbc = rp->bdry_curves->last->curve;
	    rbc_orient = rp->bdry_curves->last->orient;
	}
	find_tangent_to_curve(Node_of(tbc,tbc_orient)->posn,
	    	              Bond_at_node(tbc,tbc_orient),tbc,tbc_orient,t,
			      front);
	normal(Node_of(tbc,tbc_orient)->posn,
	       Hyper_surf_element(Bond_at_node(tbc,tbc_orient)),
	       Hyper_surf(tbc),n,front);
	if (is_obstacle_state(Right_state_at_node(tbc,tbc_orient)))
	{
	    for (i = 0; i < dim; ++i) n[i] *= -1.0;
	}
	/*
	*  As written this will take sr and sl to be the
	*  states next to the diffraction node along the respective curves.
	*  Take sr and sl to be the wall states since they at
	*  least satisfy the boundary conditions.
	*/
	sr = (curve_ang_oriented_l_to_r(i_to_f_dir,Oct.orient)) ?
	                        Right_state_at_node_of_o_curve(&Oct) :
				Left_state_at_node_of_o_curve(&Oct);
	sl = (curve_ang_oriented_l_to_r(i_to_f_dir,Osh.orient)) ?
				Left_state_at_node_of_o_curve(&Osh) :
				Right_state_at_node_of_o_curve(&Osh);

	pjump[3] = set_pjump_at_wave(Node_of_o_curve(&Oct)->posn,
			      Hyper_surf_element(Bond_at_node_of_o_curve(&Oct)),
			      Hyper_surf(Oct.curve),front,NULL);

	clear_node_flag(flag);
	dont_correct_angles_at_node(flag) = YES;
	status = phys_node_crosses_bdry(front,newfront,(POINTER)wave,rp,flag);

	if (status != GOOD_STEP)
	{
	    (void) printf("WARNING in n_diffracted_shock_exits_contact()\n");
	    (void) printf("phys_node_crosses_bdry() returns ");
	    print_time_step_status("time step status = ",status,"\n");
	    debug_print("n_d_shock","Left n_diffracted_shock_exits_contact()\n");
	    return status;
	}
	/* Reset rbc and tbc after boundary cross */

	if (!identify_curves_after_n_bdry_cross(&rbc,&rbc_orient,
			                           &tbc,&tbc_orient,
						   dncur,c_orient,
						   rp,i_to_f_dir))
	{
	    (void) printf("WARNING in n_diffracted_shock_exits_contact(), "
	                  "can't identify new phys curves after bdry cross\n");
	    debug_print("n_d_shock","Left n_diffracted_shock_exits_contact()\n");
	    return ERROR_IN_STEP;
	}

	interpolate_intfc_states(intfc) = YES;
	for (i = 0; i < 5; ++i)
	{
	    attach_phys_curve_to_n_bdry(dncur[i+1],c_orient[i+1],
					w_type[i],pt+i,sl,sr,
					intsct_coords,t,n,front,
					pjump[i],partial_dt,front);
	}

	if (!realign_phys_curves_at_nbdry(dncur,c_orient,rbc,rbc_orient,
			                     tbc,tbc_orient,
					     Opposite_ang_dir(i_to_f_dir),
					     pt,front))
	{
	    (void) printf("WARNING in n_diffracted_shock_exits_contact(), ");
	    (void) printf("unable to realign boundary curves\n");
	    debug_print("n_d_shock","Left n_diffracted_shock_exits_contact()\n");
	    return ERROR_IN_STEP;
	}
	for (i = 0; i < 5; ++i)
	{
	    if (dncur[i+1] == NULL)
		continue;
	    node_type(Node_of(dncur[i+1],c_orient[i+1])) = NEUMANN_NODE;
	    do_not_redistribute(dncur[i+1]) = YES;
	}
	interpolate_intfc_states(intfc) = sav_intrp;

	if (debugging("n_d_shock"))
	{
	    CURVE	*ctmp;
	    ORIENTATION	ctmp_orient;

	    (void) printf("States on curves after "
	                  "n_diffracted_shock_exits_contact()\n"
	                  "\nReflected Side Boundary\n");
	    verbose_print_curve_states(rbc);
	    Check_return(next_boundary(rbc,rbc_orient,&ctmp,&ctmp_orient),
			 n_diffracted_shock_exits_contact)
	    i = 1;
	    while (ctmp != tbc)
	    {
	    	(void) printf("\nIntermediate boundary %d\n",i);
	    	verbose_print_curve_states(ctmp);
	    	ctmp_orient = Opposite_orient(ctmp_orient);
	    	Check_return(next_boundary(ctmp,ctmp_orient,&ctmp,&ctmp_orient),
			     n_diffracted_shock_exits_contact)
	    	++i;
	    }
	    (void) printf("Transmitted Side Boundary\n");
	    verbose_print_curve_states(tbc);
	    for (i = 0; i < 5; ++i)
	    {
	    	if (dncur[i+1] == NULL) continue;
	    	(void) printf("\nPhysical Curve %d\n",i+1);
	    	verbose_print_curve_states(dncur[i+1]);
	    }
	}

	debug_print("n_d_shock","Left n_diffracted_shock_exits_contact()\n");

	return GOOD_STEP;
}		/*end n_diffracted_shock_exits_contact*/

LOCAL	int identify_curves_after_n_bdry_cross(
	CURVE		**rbc,
	ORIENTATION	*rbc_orient,
	CURVE		**tbc,
	ORIENTATION	*tbc_orient,
	CURVE		**dncur,
	ORIENTATION	*c_orient,
	RPROBLEM	*rp,
	ANGLE_DIRECTION	i_to_f_dir)
{
	CURVE		*c;
	int		i;
	ORIENTATION	c_or;

	if (i_to_f_dir == CLOCKWISE)
	{
	    *rbc = rp->bdry_curves->first->curve;
	    *rbc_orient = rp->bdry_curves->first->orient;
	    *tbc = rp->bdry_curves->last->curve;
	    *tbc_orient = rp->bdry_curves->last->orient;
	}
	else
	{
	    *tbc = rp->bdry_curves->first->curve;
	    *tbc_orient = rp->bdry_curves->first->orient;
	    *rbc = rp->bdry_curves->last->curve;
	    *rbc_orient = rp->bdry_curves->last->orient;
	}

	dncur[0] = dncur[6] = NULL;
	c = *tbc;	c_or = *tbc_orient;
	if (dncur[5] != NULL)
	{
	    if ((dncur[5] = find_physical_curve_at_node(Node_of(c,c_or),
						        c_orient+5)) == NULL)
	    	return NO;
	    Check_return(next_boundary(c,c_or,&c,&c_or),
			 identify_curves_after_n_bdry_cross)
	    c_or = Opposite_orient(c_or);
	}
	if (dncur[4] != NULL)
	{
	    if ((dncur[4] = find_physical_curve_at_node(Node_of(c,c_or),
						        c_orient+4)) == NULL)
	    	return NO;
	}

	c = *rbc;	c_or = *rbc_orient;
	for (i = 1; i < 4; ++i)
	{
	    if (dncur[i] != NULL)
	    {
	    	if ((dncur[i] = find_physical_curve_at_node(Node_of(c,c_or),
							    c_orient+i))==NULL)
	    	    return NO;
	    	Check_return(next_boundary(c,c_or,&c,&c_or),
			     identify_curves_after_n_bdry_cross)
	    	c_or = Opposite_orient(c_or);
	    }
	}
	return YES;
}		/*end identify_curves_after_n_bdry_cross*/


LOCAL int leaving_shock_contact(
	const char	*message,
	int		status)
{
	debug_print("shock_contact", "Leaving %s(), ",message);
	if (debugging("shock_contact"))
	    print_time_step_status("time step status = ",status,"\n");
	return status;
}		/*end leaving_shock_contact*/

LOCAL	int te_overtakes_inc_shock(
	O_CURVE		**nis,
	O_CURVE		**ois,
	O_CURVE		**nte,
	O_CURVE		**ote,
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp)
{
	ANGLE_DIRECTION	i_to_f_dir[2];
	BOND		*newbis, *newbte;
	COMPONENT	l_comp, r_comp;
	CURVE		*newc;
	Locstate	lst, rst, bst;
	POINT		*pcrs;
	NODE		*newn[2], *ns, *ne;
	NODE_FLAG       flag;
	RPROBLEM	*rptmp;
	double		tis, tte, dt_frac;
	int		i;
	int		status;
	size_t		sizest = front->sizest;
	boolean		c_ext[2];
	static	double	origin[MAXD];

	debug_print("shock_contact","Entered te_overtakes_inc_shock()\n");

	set_to_next_node_only(flag);
	for (i = 0; i < 2; ++i)
	{
	    pcrs = Point(origin);
	    lst = left_state(pcrs);	rst = right_state(pcrs);
	    dt_frac = 1.0;
	    status = cross_or_extend_to_cross_two_propagated_curves(ois[i],
	                 nis[i],ote[i],nte[i],&pcrs,&newbis,&newbte,&tis,&tte,
	                 front,(POINTER)wave,&rptmp,rp->dt,&dt_frac,
			 flag,c_ext);

	    if (status != GOOD_NODE)
	    {
	    	free_rp(rptmp);
	    	return status;
	    }
	    if (c_ext[0])
	    {
	    	ft_assign(lst,Left_state_at_node_of_o_curve(nis[i]),sizest);
	    	ft_assign(rst,Right_state_at_node_of_o_curve(nis[i]),sizest);
	    }
	    else
	    {
	    	left_state_along_bond(tis,newbis,nis[i]->curve,lst);
	    	right_state_along_bond(tis,newbis,nis[i]->curve,rst);
	    }
	    if (is_forward_wave(wave_type(nis[i]->curve)))
	    {
	    	bst = lst;
	    	i_to_f_dir[i] = (nis[i]->orient == POSITIVE_ORIENTATION) ?
				CLOCKWISE : COUNTER_CLOCK;
	    }
	    else
	    {
	    	bst = rst;
	    	i_to_f_dir[i] = (nis[i]->orient == NEGATIVE_ORIENTATION) ?
				CLOCKWISE : COUNTER_CLOCK;
	    }

	    newn[i] = make_node(pcrs);
	    node_type(newn[i]) = OVERTAKE_NODE;
	    change_node_of_curve(nis[i]->curve,nis[i]->orient,newn[i]);
	    cut_curve(pcrs,newbis,nis[i]->curve,nis[i]->orient,front,lst,rst);
	    set_status_at_node(nis[i]->curve,nis[i]->orient,OVERTOOK);
	    change_node_of_curve(nte[i]->curve,nte[i]->orient,newn[i]);
	    cut_curve(pcrs,newbte,nte[i]->curve,nte[i]->orient,front,bst,bst);
	    set_status_at_node(nte[i]->curve,nte[i]->orient,INCIDENT);
	}

	if (nis[0]->orient == NEGATIVE_ORIENTATION)
	{
	    ns = newn[0];	ne = newn[1];
	}
	else
	{
	    ne = newn[0];	ns = newn[1];
	}
	if (is_forward_wave(wave_type(nis[0]->curve)))
	{
	    r_comp = positive_component(nis[0]->curve);
	    l_comp = ((i_to_f_dir[0] == COUNTER_CLOCK &&
	    	      nte[0]->orient == POSITIVE_ORIENTATION)
		 ||
		      (i_to_f_dir[0] == CLOCKWISE &&
		       nte[0]->orient == NEGATIVE_ORIENTATION)) ?
	    	  positive_component(nte[0]->curve) :
	    	  negative_component(nte[0]->curve);
	}
	else
	{
	    l_comp = negative_component(nis[0]->curve);
	    r_comp = ((i_to_f_dir[0] == COUNTER_CLOCK &&
	    	      nte[0]->orient == POSITIVE_ORIENTATION)
	    	 ||
	    	     (i_to_f_dir[0] == CLOCKWISE &&
	    	      nte[0]->orient == NEGATIVE_ORIENTATION)) ?
	    	  positive_component(nte[0]->curve) :
	    	  negative_component(nte[0]->curve);
	}

	newc = make_curve(l_comp,r_comp,ns,ne);
	start_status(newc) = end_status(newc) = TRANSMITTED;
	wave_type(newc) = wave_type(nis[0]->curve);
	ft_assign(Left_state_at_node(newc,Opposite_orient(nis[0]->orient)),
	       Left_state_at_node_of_o_curve(nis[0]),sizest);
	ft_assign(Right_state_at_node(newc,Opposite_orient(nis[0]->orient)),
	       Right_state_at_node_of_o_curve(nis[0]),sizest);
	if (nis[1]->orient == nis[0]->orient)
	{
	    ft_assign(Left_state_at_node(newc,nis[0]->orient),
	    	   Right_state_at_node_of_o_curve(nis[1]),sizest);
	    ft_assign(Right_state_at_node(newc,nis[1]->orient),
	    	   Left_state_at_node_of_o_curve(nis[1]),sizest);
	}
	else
	{
	    ft_assign(Left_state_at_node(newc,nis[0]->orient),
		   Left_state_at_node_of_o_curve(nis[1]),sizest);
	    ft_assign(Right_state_at_node(newc,nis[0]->orient),
		   Right_state_at_node_of_o_curve(nis[1]),sizest);
	}

	if (debugging("shock_contact"))
	{
	    (void) printf("new transmitted shock\n");
	    verbose_print_curve_states(newc);
	}
	nis[0] = nis[1] = NULL;
	nte[0] = nte[1] = NULL;
	return GOOD_NODE;
}		/*end te_overtakes_inc_shock*/
	
EXPORT	int overtake_diffraction_collision(
	Front		*front,
	Front		*newfront,
	Wave		*wave,
	RPROBLEM	*rp)
{
	RPROBLEM	*newrp = NULL;
	O_CURVE		*oc, *oldoc, *oldc[4], *newc[4];
	RP_NODE		*rpn, *otake_rpn, *diff_rpn;
	double		dt_frac, dt = rp->dt;
	int		status = ERROR_IN_STEP;
	int		n_status, n_type;
	int		i;
	boolean		sav_intrp = interpolate_intfc_states(newfront->interf);
	NODE_FLAG	flag;
	UNTRACK_FLAG	uflag;

	debug_print("shock_contact","Entering overtake_diffraction_collision()\n");

	set_to_next_node_only(flag);
	if (rp->num_nod != 2)
	{
	    status = ERROR_IN_STEP;
	    (void) printf("WARNING in overtake_diffraction_collision(), "
	                  "wrong number of nodes (%d)\n",rp->num_nod);
	    goto leave;
	}
	otake_rpn = diff_rpn = NULL;
	for (rpn = rp->first_rp_node; rpn != NULL; rpn = rpn->next)
	{
	    if (node_type(rpn->old_node) == DIFFRACTION_NODE)
	    	diff_rpn = rpn;
	    if (node_type(rpn->old_node) == OVERTAKE_NODE)
	    	otake_rpn = rpn;
	}
	if (diff_rpn == NULL || otake_rpn == NULL)
	{
	    status = ERROR_IN_STEP;
	    (void) printf("WARNING in overtake_diffraction_collision(), "
	                  "Unable to identify nodes\n");
	    goto leave;
	}

	/* Identify curves */

	interpolate_intfc_states(newfront->interf) = YES;
	for (i = 0; i < 4; ++i)
	{
	    newc[i] = oldc[i] = NULL;
	}
	for (oc = rp->ang_ordered_curves->first,
	     oldoc = rp->old_ang_ordered_curves->first;
	     oc && oldoc; oc = oc->next, oldoc = oldoc->next)
	{
	    n_status = status_at_node(oldoc->curve,oldoc->orient);
	    n_type = node_type(Node_of_o_curve(oldoc));
	    if (n_status==SLIP && n_type==DIFFRACTION_NODE)
	    {
	    	oldc[1] = oldoc; newc[1] = oc;
	    }
	    else if (n_status==CONTACT_TARGET && n_type==DIFFRACTION_NODE)
	    {
		oldc[3] = oldoc; newc[3] = oc; 
	    }
	    else if (n_status==TRANSMITTED && n_type==DIFFRACTION_NODE)
	    {
	    	oldc[2] = oldoc; newc[2] = oc;
	    }
	    else if (n_status==TRANSMITTED && n_type==OVERTAKE_NODE)
	    {
	    	oldc[0] = oldoc; newc[0] = oc;
	    }
	    else
	    {
	    	set_states_at_node_by_propagate(front,(POINTER)wave,
						oldoc,oc,dt);
	    	set_untrack_flag(uflag,oc->orient,YES,NO,NO,YES,NO);
	    	(void) untrack_curve(oc,oldoc,negative_component(oc->curve),
	    		             dt,front,(POINTER)wave,rp,uflag);
	    }
	    if (oc == rp->ang_ordered_curves->last ||
	    		oldoc == rp->old_ang_ordered_curves->last)
	    	break;
	}
	if (debugging("shock_contact"))
	{
	    (void) printf("After identification of curves\n");
	    for (i = 0; i < 4; ++i)
	    {
	    	if (newc[i] == NULL)
	    	    (void) printf("newc[%d] = NULL\n",i);
	    	else
	    	{
	    	    (void) printf("newc[%d] = %llu\n",
	    	    	          i,curve_number(newc[i]->curve));
	    	    print_o_curve(newc[i]);
	    	}
	    }
	}

	if (newc[0] != NULL && newc[0]->curve != NULL)
	{
	    /* Transmitted wave at overtake tracked */

	    debug_print("shock_contact","New incident curve is tracked\n");

	    change_node_of_curve(newc[0]->curve,newc[0]->orient,diff_rpn->node);
	    set_status_at_node(newc[0]->curve,newc[0]->orient,INCIDENT);
	    delete_null_physical_curves(rp);
	    debug_print("shock_contact","Calling node_propagate\n");
	    dt_frac = 1.0;
	    if (propagation_status(diff_rpn->node) != PROPAGATED_NODE)
	    {
	    	propagation_status(diff_rpn->node) = VEL_COMPUTED_NODE;
	    	status = (*front->node_propagate)(front,(POINTER)wave,
						  diff_rpn->old_node,
						  diff_rpn->node,&newrp,dt,
						  &dt_frac,flag,
						  (POINTER)otake_rpn->old_node);
	    }
	    if (status != GOOD_NODE) 
	    {
	    	free_rp(newrp);
		screen("ERROR in overtake_diffraction_collision()\n"
		       "node propagation failed\n");
	    	status = ERROR_IN_STEP;
	    }
	    else
	    	status = GOOD_STEP;
	}
	else
	{
	    /* Transmitted wave at overtake untracked */

	    debug_print("shock_contact","New incident curve is untracked\n");

	    delete_null_physical_curves(rp);

	    /*
	    *  Since the incident curve is untracked, turn off
	    *  tracking on transmitted curve at diffraction node
	    *  if it is tracked.
	    */

	    set_states_at_node_by_propagate(front,(POINTER)wave,
					    oldc[1],newc[1],dt);
	    set_states_at_node_by_propagate(front,(POINTER)wave,
					    oldc[3],newc[3],dt);
	    if (newc[2] != NULL && newc[2]->curve != NULL)
	    {
	    	/* Turn off tracking of transmitted curve */

	    	set_untrack_flag(uflag,newc[2]->orient,YES,NO,NO,NO,NO);
	    	(void) untrack_curve(newc[2],oldc[2],
			             negative_component(newc[2]->curve),
			             dt,front,(POINTER)wave,rp,uflag);
	    	newc[2] = NULL;
	    }

	    	/* Join contacts into a single curve */

	    debug_print("shock_contact","Deleting redundant node\n");
	    (void) delete_redundant_node(diff_rpn->node,NULL,rp,front);
	}
	if (delete_node(otake_rpn->node) == FUNCTION_FAILED)
	{
	    status = ERROR_IN_STEP;
	    (void) printf("WARNING overtake_diffraction_collision(), "
	                  "Unable to delete overtake node\n");
	}
	else
	    status = GOOD_STEP;

leave:
	interpolate_intfc_states(newfront->interf) = sav_intrp;
	if (debugging("shock_contact"))
	{
	    (void) printf("Interface after overtake_diffraction_collision()\n");
	    print_interface(rp->new_intfc);
	}
	return leaving_shock_contact("overtake_diffraction_collision",status);
}		/*end overtake_diffraction_collision*/


EXPORT	void install_bdry_cross(
	CROSS		*cr,
	O_CURVE		*newc,
	RP_NODE		*rpn,
	RPROBLEM	*rp,
	Front		*front,
	SIDE		*int_side)
{
	INTERFACE	*intfc = rp->new_intfc;
	CURVE		**curs;
	COMPONENT	bdryncomp[2], bdrypcomp[2];
	Locstate	st[2];
	size_t		sizest = front->sizest;
	boolean		save_intrp;
	boolean		sav_scss = interpolate_states_at_split_curve_node();
	SIDE		side;
	Locstate        sl, sr;
	int		i;
	static	ORIENTATION	orient[2] = { NEGATIVE_ORIENTATION,
					      POSITIVE_ORIENTATION};

	save_intrp = interpolate_intfc_states(intfc);
	interpolate_intfc_states(intfc) = YES;
	if (cr->prev)
	    cr->prev->next = cr->next;
	if (cr->next)
	    cr->next->prev = cr->prev;
	if (insert_point_in_bond(cr->p,cr->b1,cr->c1) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in install_bdry_cross(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	rcl_after_insert_point(cr,cr->p,cr->b1);
	if (newc->orient == POSITIVE_ORIENTATION)
	    cr->b1 = cr->b1->next;
	cut_curve(cr->p,cr->b1,cr->c1,newc->orient,front,
	          left_state(cr->p),right_state(cr->p));
	cr->p = Node_of(cr->c1,newc->orient)->posn;
	nearest_boundary_point(Coords(cr->p),Coords(cr->p),front->rect_grid);
	side = physical_side_of_bdry_curve(cr->c2);
	if ((side == POSITIVE_SIDE && newc->orient == POSITIVE_ORIENTATION) ||
	    (side == NEGATIVE_SIDE && newc->orient == NEGATIVE_ORIENTATION))
	{
	    st[0] = Right_state_at_node_of_o_curve(newc);
	    st[1] = Left_state_at_node_of_o_curve(newc);
	    if (side == POSITIVE_SIDE)
	    {
	    	bdrypcomp[0] = positive_component(newc->curve);
	    	bdrypcomp[1] = negative_component(newc->curve);
	    	bdryncomp[0] = bdryncomp[1] = negative_component(cr->c2);
	    }
	    else
	    {
	    	bdryncomp[0] = positive_component(newc->curve);
	    	bdryncomp[1] = negative_component(newc->curve);
	    	bdrypcomp[0] = bdrypcomp[1] = negative_component(cr->c2);
	    }
	}
	else
	{
	    st[0] = Left_state_at_node_of_o_curve(newc);
	    st[1] = Right_state_at_node_of_o_curve(newc);
	    if (side == POSITIVE_SIDE)
	    {
	    	bdrypcomp[0] = negative_component(newc->curve);
	    	bdrypcomp[1] = positive_component(newc->curve);
	    	bdryncomp[0] = bdryncomp[1] = negative_component(cr->c2);
	    }
	    else
	    {
	    	bdryncomp[0] = negative_component(newc->curve);
	    	bdryncomp[1] = positive_component(newc->curve);
	    	bdrypcomp[0] = bdrypcomp[1] = negative_component(cr->c2);
	    }
	}
	interpolate_intfc_states(intfc) = NO;
	if (insert_point_in_bond(cr->p,cr->b2,cr->c2) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in install_bdry_cross(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	rcl_after_insert_point(cr,cr->p,cr->b2);
	set_interpolate_states_at_split_curve_node(NO);
	curs = split_curve(cr->p,cr->b2,cr->c2,bdryncomp[0],bdrypcomp[0],
	    		   bdryncomp[1],bdrypcomp[1]);
	set_interpolate_states_at_split_curve_node(sav_scss);
	roclists_after_split(rp,cr->c2,curs,YES);
	rcl_after_split(cr,cr->p,cr->b2,cr->c2,curs);
	set_status_at_node(curs[0],orient[0],FIXED);
	set_status_at_node(curs[1],orient[1],FIXED);
	switch (wave_type(curs[0]))
	{
	case SUBDOMAIN_BOUNDARY:
	    for (i = 0; i < 2; ++i)
	    {
	    	obstacle_state(front->interf,
			       Left_state_at_node(curs[i],orient[i]),sizest);
	    	obstacle_state(front->interf,
			       Right_state_at_node(curs[i],orient[i]),sizest);
	    }
	    break;
	case DIRICHLET_BOUNDARY:
	case NEUMANN_BOUNDARY:
	    for (i = 0; i < 2; ++i)
	    {
	    	sl = Left_state_at_node(curs[i],orient[i]);
	    	sr = Right_state_at_node(curs[i],orient[i]);
	    	if (int_side != NULL && *int_side == POSITIVE_SIDE)
	    	{
	    	    ft_assign(sr,st[i],sizest);
	    	    obstacle_state(front->interf,sl,sizest);
	    	}
	    	else if (int_side != NULL && *int_side == NEGATIVE_SIDE)
	    	{
	    	    ft_assign(sl,st[i],sizest);
	    	    obstacle_state(front->interf,sr,sizest);
	    	}
	    	else if (is_excluded_comp(positive_component(curs[i]),intfc))
	    	{
	    	    ft_assign(sl,st[i],sizest);
	    	    obstacle_state(front->interf,sr,sizest);
	    	}
	    	else
	    	{
	    	    ft_assign(sr,st[i],sizest);
	    	    obstacle_state(front->interf,sl,sizest);
	    	}
	    }
	    break;
	}
	rpn->states_assigned_at_node = YES;
	interpolate_intfc_states(intfc) = save_intrp;
}		/*end install_bdry_cross*/

LOCAL	int new_cross_bifurcate(
	NODE		*newn,
	O_CURVE		*oc0,
	O_CURVE		*oc4,
	Front		*front,
	RPROBLEM	*rp)
{
	CURVE		*c0 = oc0->curve, *c4 = oc4->curve;
	ORIENTATION	c0_orient = oc0->orient, c4_orient = oc4->orient;
	CURVE		*cur;

	if (newn == NULL) return YES;
	if (c0 == NULL && c4 == NULL) return YES;
	set_incident_status(c0,c0_orient);
	set_incident_status(c4,c4_orient);

	/*
	*	The present version joins the two incident curves at the
	*  	cross node.  Other resolution might include untracking
	*	the shocks or the bifurcation to mach reflection.
	*/
	

	interpolate_intfc_states(c0->interface) = YES;
	if (c0_orient == c4_orient)
	{
	    invert_curve(c0);
	    c0_orient = Opposite_orient(c0_orient);
	}
	if (c0 == c4)
	{
	    size_t 	    sizest = front->sizest;
	    static Locstate	st = NULL;

	    if (st == NULL)
	    	alloc_state(front->interf,&st,sizest);
	    cur = c0;
	    node_type(newn) = CLOSED_NODE;
	    interpolate_states(front,0.5,0.5,Coords(newn->posn),
			       left_start_state(cur),Coords(newn->posn),
			       left_end_state(cur),st);
	    ft_assign(left_start_state(cur),st,sizest);
	    ft_assign(left_end_state(cur),st,sizest);
	    interpolate_states(front,0.5,0.5,Coords(newn->posn),
			       right_start_state(cur),Coords(newn->posn),
			       right_end_state(cur),st);
	    ft_assign(right_start_state(cur),st,sizest);
	    ft_assign(right_end_state(cur),st,sizest);
	}
	else if (c4_orient == POSITIVE_ORIENTATION)
	{
	    cur = join_curves(c0,c4,negative_component(c4),
	    	              positive_component(c4),(BOND **)NULL);

	    /* Do not allow redistribute to create a tangle. */
	    do_not_redistribute(cur) = YES;

	    roclists_after_join(rp,c0,oc0,c4,oc4,cur);
	    (void) delete_node(newn);
	}
	else
	{
	    cur = join_curves(c4,c0,negative_component(c4),
			      positive_component(c4),(BOND **)NULL);

	    /* Do not allow redistribute to create a tangle. */
	    do_not_redistribute(cur) = YES;

	    roclists_after_join(rp,c4,oc4,c0,oc0,cur);
	    (void) delete_node(newn);
	}
#if defined(DEBUG_SHOCK_CONTACT)
	debug_print("shock_contact","Curves %d and %d merged into curve %d\n",
	      c0,c4,cur);
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	return YES;
}		/*end new_cross_bifurcate*/


/*
*			estimate_node_vel():
*
*	This function provides an estimate of the node velocity
*	of a given node.   This estimate is found by calculating
*	the intersection of the lines through the propagated
*	bonds at the node of two given curves and using this
*	intersection as the new node position.  Since this
*	function may be called for curves on the new interface
*	before the interior propagation step, the second order
*	corrections to the Riemann problem at the points to
*	be propagated can not be included.  The states adjacent
*	to these points are just set equal to the states at
*	the points.
*
*/

EXPORT	int estimate_node_vel(
	CURVE		*c1,
	ORIENTATION	c1_orient,
	double		*t1,
	CURVE		*c2,
	ORIENTATION	c2_orient,
	double		*t2,
	ANGLE_DIRECTION	c2_to_c1_dir,
	int		n_type,
	double		*node_v,
	Front		*front)
{
	int		dim = front->interf->dim;
	Locstate	st1a, st1b, st2a, st2b;
	static Locstate st0 = NULL;

	debug_print("node_vel","Entering estimate_node_vel()\n");

	if (st0 == NULL)
	{
	    alloc_state(front->interf,&st0,front->sizest);
	}

	find_tangent_to_curve(Node_of(c1,c1_orient)->posn,
	    	              Bond_at_node(c1,c1_orient),
			      c1,c1_orient,t1,front);
	find_tangent_to_curve(Node_of(c2,c2_orient)->posn,
			      Bond_at_node(c2,c2_orient),
			      c2,c2_orient,t2,front);
	find_states_at_node_of_interaction(c1_orient,c2_orient,c2_to_c1_dir,
	    Left_state_at_node(c1,c1_orient),Right_state_at_node(c1,c1_orient),
	    Left_state_at_node(c2,c2_orient),Right_state_at_node(c2,c2_orient),
	    &st1a,&st1b,&st2a,&st2b);

	interpolate_states(front,0.5,0.5,
		           Coords(Node_of(c1,c1_orient)->posn),st1a,
		           Coords(Node_of(c2,c2_orient)->posn),st2a,st0);

	if (compute_node_velocity(wave_type(c1),wave_type(c2),
	                          st0,st1b,st2b,t1,t2,node_v,dim,n_type,
			          c2_to_c1_dir))
	    return YES;
	
	return NO;
}		/*end estimate_node_vel*/



/*
*			compute_node_velocity():
*
*	Computes the node velocity for a series of nodes.  st0 is the state
*	between the interacting curves (labelled 1 and 2), with st1b and
*	st2b the corresponding behind states.
*
*	DIFFRACTION_NODE:
*	B_REFLECT_NODE:
*	TOT_INT_REFL_NODE:  The mass flux across the incident is computed,
*	yielding q0 = m/(rho0*sinb) where sinb is the angle between the 
*	incident shock and contact.  The steady state velocity is then
*	parallel to the contact, so node_v = v0 + q*t where t is the tangent
*	to the contact.  No minus is needed because t points in the opposite
*	direction from the steady state velocity.
*	
*	CROSS_NODE: For a cross node we write 
*			(node_v - vel0) = alpha1*t1 + alpha2*t2.
*       Thus we are writing the steady ahead velocity as a linear combination
*	of the tangents t1 and t2.  Crossing this equation with t1 yields
*	alpha2 since (steady vel)xt1 = q0*sinb1 = m1/rho0 where q0 is the
*	steady ahead speed and b1 is the incident angle of the steady ahead
*	flow on wave1.  A similar process yields alpha2, and thus node_v.
*
*	OVERTAKE_NODE: This is very similar to a cross node.  We write the
*	same linear combination as in that case.  However, the mass flux
*	across wave 1 (and thus alpha1) needs an extra minus sign.  This is
*	due to the unintuitive labelling of the states in this case.  st0
*	is still between the two waves even though it is not really the ahead
*	state for an overtake configuration.  Thus the computation of the
*	mass flux across wave 1 is going in the wrong direction (from behind
*	to ahead).
*
*	CC_NODE: Here we consider the interaction of two contacts, so the
*	node is simply being convected with the ahead flow.  Thus the node
*	velocity is the same as the ahead state velocity.
*/

EXPORT	int compute_node_velocity(
	int		wt1,	/* wave type of incident1 */
	int		wt2,	/* wave type of incident2 */
	Locstate	st0,
	Locstate	st1b,
	Locstate	st2b,
	double		*t1,	/* incident1 tangent */
	double		*t2,	/* incident2 tangent */
	double		*node_v,
	int		dim,
	int		n_type,
	ANGLE_DIRECTION	c2_to_c1_dir)
{
	double		sinb, cosb;		/* b is angle of t1 and t2 */
	double		m, m1, m2;
	double		q0, rho0, c0;
	double		*t;
	int		i;
	static const double VEL_TOL = 1.0e-12;	/*TOLERANCE*/

	debug_print("node_vel","Entered compute_node_velocity()\n");
	cosb = scalar_product(t2,t1,dim);
	(void) vector_product(t2,t1,&sinb,dim);
	if (debugging("node_vel"))
	{
	    double b;
	    b = angle(cosb,sinb);
	    (void) printf("t1 = <%g, %g>, mag = %g\n",t1[0],t1[1],
	    	          mag_vector(t1,dim));
	    (void) printf("t2 = <%g, %g>, mag = %g\n",t2[0],t2[1],
			  mag_vector(t2,dim));
	    (void) printf("cos(b)  = %g,  sin(b) = %g, ",cosb,sinb);
	    print_angle("b =",b,"\n");
	    verbose_print_state("st0",st0);
	    verbose_print_state("st1b",st1b);
	    verbose_print_state("st2b",st2b);
	}
	if (cosb < 0.0 && n_type != CROSS_NODE)
	{
	    (void) printf("WARNING in compute_node_velocity(), "
	                  "Angle between interacting curves > PI/2\n");
	    debug_print("node_vel","Left compute_node_velocity(), ans = NO\n");
	    return NO;
	}
	else if ((n_type == CROSS_NODE) &&
		  ((c2_to_c1_dir == CLOCKWISE && sinb > 0.0) ||
		  (c2_to_c1_dir == COUNTER_CLOCK && sinb < 0.0)))
	{
	    (void) printf("WARNING in compute_node_velocity(), "
		          "Cross node with angle between "
	                  "interacting curves > PI\n");
	    debug_print("node_vel","Left compute_node_velocity(), ans = NO\n");
	    return NO;
	}
	sinb = fabs(sinb);
	rho0 = Dens(st0);
	c0 = sound_speed(st0);
	switch (n_type)
	{
	case B_REFLECT_NODE:
	case DIFFRACTION_NODE:
	    if (is_vector_wave(wt2))
	    {
	        m = mass_flux(pressure(st2b),st0);
	        t = t1;
	    }
	    else if (is_vector_wave(wt1))
	    {
	        m = mass_flux(pressure(st1b),st0);
	        t = t2;
	    }
	    else
	    {
	        (void) printf("WARNING in estimate_node_vel(), "
	                      "Shock interaction with not vector wave\n");
	        debug_print("node_vel",
	    	      "Left compute_node_velocity(), ans = NO\n");
	        return NO;
	    }
	    if (rho0*c0*sinb <= VEL_TOL*m)
	    {
	    	(void) printf("WARNING in estimate_node_vel(), "
	    	              "Node velocity too large\n"
	       	              "m = %g, rho0 = %g, c0 = %g, sinb = %g\n",
		              m,rho0,c0,sinb);
	    	(void) printf("rho0*c0*sinb = %g, VEL_TOL*m = %g\n",
			      rho0*c0*sinb, VEL_TOL*m);
	    	debug_print("node_vel","Left compute_node_velocity(), ans = NO\n");
	    	return NO;
	    }
	    q0 = m/(rho0*sinb);
	    for (i = 0; i < dim; ++i)
	    	node_v[i] = vel(i,st0) + q0*t[i];
	    debug_print("node_vel","Left compute_node_velocity(), ans = YES\n");
	    return YES;
	case CROSS_NODE:
	    m1 = mass_flux(pressure(st1b),st0);
	    m2 = mass_flux(pressure(st2b),st0);
	    for (i = 0; i < dim; ++i)
	    {
	    	node_v[i] = vel(i,st0) + (m2*t1[i] + m1*t2[i])/(rho0*sinb);
	    }
	    debug_print("node_vel","Left compute_node_velocity(), ans = YES\n");
	    return YES;
	case OVERTAKE_NODE:
	    m1 = mass_flux(pressure(st1b),st0);
	    m2 = mass_flux(pressure(st2b),st0);
	    for (i = 0; i < dim; ++i)
	    {
	    	node_v[i] = vel(i,st0) + (m2*t1[i] - m1*t2[i])/(rho0*sinb);
	    }
	    debug_print("node_vel","Left compute_node_velocity(), ans = YES\n");
	    return YES;
	case CC_NODE:
	    for (i = 0; i < dim; ++i)
	    	node_v[i] = vel(i,st0);
	    debug_print("node_vel","Left compute_node_velocity(), ans = YES\n");
	    return YES;
	}
	debug_print("node_vel","Left compute_node_velocity(), ans = NO\n");
	return NO;
}		/*end compute_node_velocity*/

EXPORT	void find_states_at_node_of_interaction(
	ORIENTATION	c1_orient,
	ORIENTATION	c2_orient,
	ANGLE_DIRECTION	c2_to_c1_dir,
	Locstate	s1l,
	Locstate	s1r,
	Locstate	s2l,
	Locstate	s2r,
	Locstate	*st1a,
	Locstate	*st1b,
	Locstate	*st2a,
	Locstate	*st2b)
{
	if (c2_to_c1_dir == COUNTER_CLOCK)
	{
	    if (c2_orient == POSITIVE_ORIENTATION)
	    {
	    	*st2a = s2l;
	    	*st2b = s2r;
	    }
	    else
	    {
	    	*st2b = s2l;
	    	*st2a = s2r;
	    }
	    if (c1_orient == POSITIVE_ORIENTATION)
	    {
	    	*st1a = s1r;
	    	*st1b = s1l;
	    }
	    else
	    {
	    	*st1b = s1r;
	    	*st1a = s1l;
	    }
	}
	else
	{
	    if (c2_orient == NEGATIVE_ORIENTATION)
	    {
	    	*st2a = s2l;
	    	*st2b = s2r;
	    }
	    else
	    {
	    	*st2b = s2l;
	    	*st2a = s2r;
	    }
	    if (c1_orient == NEGATIVE_ORIENTATION)
	    {
	    	*st1a = s1r;
	    	*st1b = s1l;
	    }
	    else
	    {
	    	*st1b = s1r;
	    	*st1a = s1l;
	    }
	}
}		/*end find_states_at_node_of_interaction*/
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */
