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
*				grp.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Solves two dimensional Riemann problems.
*
*/

#if defined(TWOD)

#define DEBUG_STRING "2drp"
#include <gdecs/gdecs.h>
	/* LOCAL Function Declarations */
LOCAL	boolean	is_refl_curve_crosses_corner(RPROBLEM*);
LOCAL	boolean	prepare_curves_for_bdry_cross(O_CURVE*,O_CURVE*,O_CURVE*,
					      ANGLE_DIRECTION,RP_NODE**,
					      Front*,Wave*,RPROBLEM*);
LOCAL	int	connect_interacting_curves(Front*,Wave*,RPROBLEM*);
LOCAL	int	g_interior_rp(Front*,Front*,Wave*,RPROBLEM*);
LOCAL   int     incident_curves_exit_at_bdry(Front*,Front*,Wave*,RPROBLEM*);
LOCAL	int	num_curves_with_wave_type_at_node(NODE*,int);
LOCAL	int	refl_curve_crosses_corner(Front*,Front*,Wave*,RPROBLEM*);
LOCAL	int	remove_oc_duplicates(O_CURVE_FAMILY**);
LOCAL	void	g_delete_curve_from_rp_node(CURVE*,RP_NODE*,RPROBLEM*);
LOCAL	void	g_init_rp_nodes(RPROBLEM*);
LOCAL	void	g_set_phys_ocurves_to_null(RP_NODE*,RPROBLEM*);
LOCAL	void	g_set_rp_statistics(RPROBLEM*);
LOCAL	void	g_user_free_rp_node(RP_NODE*,RPROBLEM*);
LOCAL	void	g_user_print_rp_node(RP_NODE*,RPROBLEM*);
LOCAL	void	g_user_print_rproblem(RPROBLEM*);
LOCAL	void	init_physical_ocurves(RPROBLEM*rp);
LOCAL	void	print_rp_statistical_data(RPROBLEM*);
LOCAL	void	replace_null_curves_in_family(O_CURVE_FAMILY**,RPROBLEM*);
LOCAL	void	substitute_contact_continuation(O_CURVE**,O_CURVE_FAMILY**);
LOCAL	void	substitute_shock_continuation(O_CURVE**,O_CURVE_FAMILY**);

#if defined(FULL_PHYSICS)
LOCAL	boolean	is_deprecursor(RPROBLEM*);
LOCAL	boolean	is_diffracted_shock_exits_contact(RPROBLEM*);
LOCAL	boolean	is_n_diffracted_shock_exits_contact(RPROBLEM*);
LOCAL	boolean	is_overtake_diffraction_collision(RPROBLEM*);
LOCAL	boolean	is_overtake_overtake_collision(RPROBLEM*);
LOCAL	int	deprecurse_shock_refraction(Front*,Wave*,RPROBLEM*);
LOCAL	int	is_premature_deprecursion(Front*,Wave*,RPROBLEM*);
LOCAL	int	overtake_overtake_collision(Front*,Front*,Wave*,RPROBLEM*);
LOCAL	int	repropagate_node(Front*,Wave*,RPROBLEM*,int,NODE_FLAG);
LOCAL	int	untrack_precursor(Front*,Wave*,double,RPROBLEM*);
#endif /* defined(FULL_PHYSICS) */

#if !defined(__INTEL_COMPILER)
#pragma	noinline	is_refl_curve_crosses_corner
#pragma	noinline	prepare_curves_for_bdry_cross
#pragma	noinline	connect_interacting_curves
#pragma	noinline	g_interior_rp
#pragma	noinline	num_curves_with_wave_type_at_node
#pragma	noinline	refl_curve_crosses_corner
#pragma	noinline	remove_oc_duplicates
#pragma	noinline	find_curves_with_status
#pragma	noinline	g_delete_curve_from_rp_node
#pragma	noinline	g_init_rp_nodes
#pragma	noinline	g_set_phys_ocurves_to_null
#pragma	noinline	g_set_rp_statistics
#pragma	noinline	g_user_free_rp_node
#pragma	noinline	g_user_print_rp_node
#pragma	noinline	g_user_print_rproblem
#pragma	noinline	init_physical_ocurves
#pragma	noinline	replace_null_curves_in_family
#pragma	noinline	substitute_contact_continuation
#pragma	noinline	substitute_shock_continuation

#if defined(FULL_PHYSICS)
#pragma	noinline	is_deprecursor
#pragma	noinline	is_diffracted_shock_exits_contact
#pragma	noinline	is_n_diffracted_shock_exits_contact
#pragma	noinline	is_overtake_diffraction_collision
#pragma	noinline	is_overtake_overtake_collision
#pragma	noinline	deprecurse_shock_refraction
#pragma	noinline	is_premature_deprecursion
#pragma	noinline	overtake_overtake_collision
#pragma	noinline	repropagate_node
#pragma	noinline	untrack_precursor
#endif /* defined(FULL_PHYSICS) */
#endif /*!defined(__INTEL_COMPILER)*/


EXPORT	void	set_g_rproblem_hooks(void)
{
	F_USER_RPROBLEM *rpuh;

	DEBUG_ENTER(set_g_rproblem_hooks)
	rpuh = rp_user_hook();
	rpuh->size_rproblem = sizeof(G_RPROBLEM);
	rpuh->size_rp_node = sizeof(G_RP_NODE);
	rpuh->_init_rp_nodes = g_init_rp_nodes;
	rpuh->_delete_curve_from_rp_node = g_delete_curve_from_rp_node;
	rpuh->_user_free_rp_node = g_user_free_rp_node;
	rpuh->_user_print_rp_node = g_user_print_rp_node;
	rpuh->_user_print_rproblem = g_user_print_rproblem;
	rpuh->_set_rp_statistics = g_set_rp_statistics;
	rpuh->_set_phys_ocurves_to_null = g_set_phys_ocurves_to_null;
	DEBUG_LEAVE(set_g_rproblem_hooks)
}		/*end set_g_rproblem_hooks*/


/*
*			g_init_2drproblem():
*/

EXPORT void g_init_2drproblem(
	RPROBLEM	*rp,
	Front		*front)
{
	DEBUG_ENTER(g_init_2drproblem)
	init_rp_nodes(rp);
	init_ocurve_lists(rp,front);
	init_physical_ocurves(rp);
	set_rp_statistics(rp);
	if (DEBUG)
	{
	    (void) printf("Initialized rproblem %p\n",(POINTER)rp);
	    print_rproblem(rp);
	}
	DEBUG_LEAVE(g_init_2drproblem)
}		/*end g_init_2drproblem*/



/*
*			g_2drproblem():
*
*	Recognized return statuses are
*		ERROR_IN_STEP, GOOD_STEP, MODIFY_TIME_STEP
*/

EXPORT int g_2drproblem(
	Front		*front,
	Front		*newfront,
	POINTER		p2wave,
	RPROBLEM	**prp)
{
	Wave		*wave = (Wave*)p2wave;
	RPROBLEM	*rp = *prp;
	int		status;

	DEBUG_ENTER(g_2drproblem)
	if (DEBUG)
	{
	    (void) printf("Rproblem into g_2drproblem()\n");
	    print_rproblem(rp);
	}
	(*front->init_2drproblem)(rp,front);
	if (is_bdry_type(rp)) 
	{
	    if (DEBUG)
	    {
	    	(void) printf("\tBdry_type Riemann Problem:\n");
	    	print_rp_statistical_data(rp);
	    }
	    if ((rp->bdry_type1 == SUBDOMAIN_BOUNDARY) ||
	    	rp->bdry_type2 == SUBDOMAIN_BOUNDARY)
	    {
	    	RP_NODE *rp_n;

	    	if (  (rp->bdry_type1 == DIRICHLET_BOUNDARY ||
	    	         rp->bdry_type1 == NEUMANN_BOUNDARY)
		     ||
		    ( (rp->bdry_type2 == DIRICHLET_BOUNDARY ||
	    	       rp->bdry_type2 == NEUMANN_BOUNDARY)))
	    	{
	    	    status = pp_curve_exits_at_bdry(front,newfront,wave,prp);
	    	}
	    	else
	    	{
	    	    if (DEBUG)
	    	    {
	    	      (void) printf("WARNING in g_2drproblem(), "
	    	                    "interaction involving subdomain "
	    	                    "boundaries\n\tCalling "
	    	                    "set_node_states_and_continue()\n");
	    	    }
	    	    for (rp_n = rp->first_rp_node;
	    	    	 rp_n != NULL; rp_n = rp_n->next)
	    	    {
	    	        if (!set_node_states_and_continue(rp_n->old_node,
							     rp_n->node,front))
	    	        {
	    	    	    DEBUG_LEAVE(g_2drproblem)
	    	    	    return ERROR_IN_STEP;
	    	        }
	    	        rp_n->states_assigned_at_node = YES;
	    	    }
	    	    status = GOOD_STEP;
	        }
	    }
	    else if (rp->num_phys == 1 && 
	    	     rp->bdry_type1 == DIRICHLET_BOUNDARY && 
		     rp->bdry_type2 == DIRICHLET_BOUNDARY)
	    {
		NODE_FLAG flag;
		clear_node_flag(flag);
	    	status = phys_node_crosses_bdry(front,newfront,p2wave,rp,flag);
	    }
#if defined(FULL_PHYSICS)
	    else if (is_n_diffracted_shock_exits_contact(rp) == YES)
	    {
	    	status = n_diffracted_shock_exits_contact(front,newfront,wave,
							  rp);
	    }
	    else if (rp_num_inc_shock(rp) == 1 && rp_num_contact(rp) == 1)
	    {
	    	status = shock_contact_rp(front,wave,rp);
	    }
#endif /* defined(FULL_PHYSICS) */
	    else if (rp_num_incident(rp) == 0 && rp->num_fxd == 2)
	    	status = curve_exits_parallel_to_bdry(front,p2wave,rp);
	    else if (is_refl_curve_crosses_corner(rp) == YES)
	    	status = refl_curve_crosses_corner(front,newfront,wave,rp);
	    else if (rp_num_incident(rp) == 0)
	    {
	    	status = incident_curves_exit_at_bdry(front,newfront,wave,rp);
	    }
	    else if (rp_num_incident(rp) == 1 && rp->num_fxd == 0)
	    	status = bdry_rp_1i_0f(front,wave,rp);
	    else if (rp_num_incident(rp) >= 1 && rp->num_fxd == 1)
	    {
	    	status = incident_curve_crosses_fixed_node(front,newfront,
							   p2wave,rp);
	    }
	    else if (rp_num_incident(rp) == 2)
	    	status = bdry_rp_2i(front,newfront,wave,rp);
	    else 
	    {
	        (void) printf("WARNING: unexpected case in g_2drproblem\n");
	    	(void) print_rp_statistical_data(rp);
		if (DEBUG)
	            print_rproblem(rp);
	        DEBUG_LEAVE(g_2drproblem)
	        return ERROR_IN_STEP;
	    }
	}
	else
	    status = g_interior_rp(front,newfront,wave,rp);
	DEBUG_LEAVE(g_2drproblem)
	return status;
}		/*end g_2drproblem*/

/*
*			init_physical_ocurves():
*
*	Initializes the physical curves associated with a given RPROBLEM.
*/

LOCAL void init_physical_ocurves(
	RPROBLEM	*rp)
{
	RP_NODE		*rp_node;
	CURVE		*c1,*c2;
	ORIENTATION	orient1,orient2;

	DEBUG_ENTER(init_physical_ocurves)

		/* Initialize physical curves */

	for (rp_node = rp->first_rp_node; rp_node; rp_node = rp_node->next) 
	{
	    free_o_curve_family(rpn_contact1(rp_node));
	    rpn_contact1(rp_node) = NULL;
	    free_o_curve_family(rpn_contact2(rp_node));
	    rpn_contact2(rp_node) = NULL;

	    free_o_curve_family(rpn_inc_shock1(rp_node)); 
	    rpn_inc_shock1(rp_node) = NULL;
	    free_o_curve_family(rpn_inc_shock2(rp_node)); 
	    rpn_inc_shock2(rp_node) = NULL;

	    free_o_curve_family(rpn_reflected1(rp_node)); 
	    rpn_reflected1(rp_node) = NULL;
	    free_o_curve_family(rpn_reflected2(rp_node)); 
	    rpn_reflected2(rp_node) = NULL;

	    free_o_curve_family(rpn_transmitted(rp_node)); 
	    rpn_transmitted(rp_node) = NULL;

	    find_curves_with_wave_type(rp_node->node,&c1,&orient1,&c2,&orient2,
				       CONTACT);
	    if (c1)
		init_cfamily(&rpn_contact1(rp_node),c1,orient1);
	    if (c2)
		init_cfamily(&rpn_contact2(rp_node),c2,orient2);

	    find_curves_with_status(rp_node->node,&c1,&orient1,
	    			    &c2,&orient2,INCIDENT);
	    if (c1 && wave_type(c1) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
	    	init_cfamily(&rpn_inc_shock1(rp_node),c1,orient1);
	    if (c2 && wave_type(c2) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
	    	init_cfamily(&rpn_inc_shock2(rp_node),c2,orient2);

	    find_curves_with_status(rp_node->node,&c1,&orient1,
	    			    &c2,&orient2,REFLECTED);
	    if (c1)
		init_cfamily(&rpn_reflected1(rp_node),c1,orient1);
	    if (c2)
		init_cfamily(&rpn_reflected2(rp_node),c2,orient2);

	    find_curves_with_status(rp_node->node,&c1,&orient1,
				    &c2,&orient2,TRANSMITTED);
	    if (c1)
		init_cfamily(&rpn_transmitted(rp_node),c1,orient1);
	}

		/* Replace null curves by continuations */

	for (rp_node = rp->first_rp_node; rp_node; rp_node = rp_node->next) 
	{
	    replace_null_curves_in_family(&rpn_inc_shock1(rp_node),rp);
	    replace_null_curves_in_family(&rpn_inc_shock2(rp_node),rp);
	    replace_null_curves_in_family(&rpn_reflected1(rp_node),rp);
	    replace_null_curves_in_family(&rpn_reflected2(rp_node),rp);
	    replace_null_curves_in_family(&rpn_transmitted(rp_node),rp);
	    replace_null_curves_in_family(&rpn_contact1(rp_node),rp);
	    replace_null_curves_in_family(&rpn_contact2(rp_node),rp);

	    	/* Null curves and families prefered in location 2 */

	    relocate_null_pointer((POINTER *)&rpn_inc_shock1(rp_node),
	    	                  (POINTER *)&rpn_inc_shock2(rp_node));
	    relocate_null_pointer((POINTER *)&rpn_reflected1(rp_node),
			          (POINTER *)&rpn_reflected2(rp_node));
	    relocate_null_pointer((POINTER *)&rpn_contact1(rp_node),
			          (POINTER *)&rpn_contact2(rp_node));
	}
	DEBUG_LEAVE(init_physical_ocurves)
}		/*end init_physical_ocurves*/


/*
*			replace_null_curves_in_family():
*
*	All null curves in an O_CURVE_FAMILY associated with a RPROBLEM
*	are replaced by their continuation. If the continuation is
*	nonunique, this generates new O_CURVES within the same family.
*	The continuation of a reflected wave through a node is not defined,
*	so the algorithm cannot handle such cases.
*/

LOCAL void replace_null_curves_in_family(
	O_CURVE_FAMILY	**cfamily,
	RPROBLEM	*rp)
{
	O_CURVE		*oc;
	O_CURVE		Baseoc;

	DEBUG_ENTER(replace_null_curves_in_family)
		/* Test for null curves */
	
	if (!*cfamily)
	{
	    DEBUG_LEAVE(replace_null_curves_in_family)
	    return;
	}

	Baseoc.curve = NULL;
testnull:
	for (oc = (*cfamily)->first; oc; oc = oc->next)
	    if (is_null_curve(oc->curve,rp)) break;

	if (!oc)
	{
	    DEBUG_LEAVE(replace_null_curves_in_family)
	    return;
	}

	if (DEBUG)
	{
	    (void) printf("Null curve found\n");
	    print_curve(oc->curve);
	}

	if (Baseoc.curve == NULL) Baseoc = *oc;

	    /* Replace null curves */

	oc->orient = Opposite_orient(oc->orient);
	switch(node_type(Node_of(oc->curve,oc->orient))) 
	{
	case PASSIVE_NODE:
	case DIRICHLET_NODE:
	case NEUMANN_NODE:
	case FIXED_NODE:
	case B_REFLECT_NODE:
	case ATTACHED_B_NODE:
	case SUBDOMAIN_NODE:
#if defined(FULL_PHYSICS)
	case WAVE_END_NODE:
#endif /* defined(FULL_PHYSICS) */
	    delete_oc_curve_from_family(&oc,cfamily);
	    break;
#if defined(FULL_PHYSICS)
	case OVERTAKE_NODE:
	    if (is_scalar_wave(wave_type(oc->curve)))
	    	substitute_contact_continuation(&oc,cfamily);

	    /* allow rarefactions as incidents at overtake nodes */
	    if (is_rarefaction_wave(wave_type(oc->curve)) ||
	    	        is_shock_wave(wave_type(oc->curve)))
	    	substitute_shock_continuation(&oc,cfamily);

	    /* Check for possible circular loop of null curves */

	    if (oc && oc->curve==Baseoc.curve && oc->orient==Baseoc.orient)
	    {
	    	delete_oc_curve_from_family(&oc,cfamily);
	    	Baseoc.curve = NULL;
	    }
	    break;
	case DIFFRACTION_NODE:
	case TRANSMISSION_NODE:
	case CROSS_NODE:
#endif /* defined(FULL_PHYSICS) */
	case MACH_NODE:
	    if (is_scalar_wave(wave_type(oc->curve)))
	    	substitute_contact_continuation(&oc,cfamily);
	    if (is_shock_wave(wave_type(oc->curve)))
	    	substitute_shock_continuation(&oc,cfamily);

	    /* Check for possible circular loop of null curves */

	    if (oc && oc->curve==Baseoc.curve && oc->orient==Baseoc.orient)
	    {
	    	delete_oc_curve_from_family(&oc,cfamily);
	    	Baseoc.curve = NULL;
	    }
	    break;
#if defined(FULL_PHYSICS)
	case CC_NODE:
	    substitute_contact_continuation(&oc,cfamily);

	    /* Check for possible circular loop of null curves */

	    if (oc && oc->curve==Baseoc.curve && oc->orient==Baseoc.orient)
	    {
	    	delete_oc_curve_from_family(&oc,cfamily);
	    	Baseoc.curve = NULL;
	    }
	    break;
#endif /* defined(FULL_PHYSICS) */
	}
	if (*cfamily) goto testnull;
	DEBUG_LEAVE(replace_null_curves_in_family)
}		/*end replace_null_curves_in_family*/

/*
*		substitute_contact_continuation():
*/

LOCAL void substitute_contact_continuation(
	O_CURVE		**poc,
	O_CURVE_FAMILY	**cfamily)
{
	O_CURVE		*oc;
#if defined(FULL_PHYSICS)
	CURVE		**c;
#endif /* defined(FULL_PHYSICS) */

	DEBUG_ENTER(substitute_contact_continuation)

	if (poc == NULL)
	{
	    DEBUG_LEAVE(substitute_contact_continuation)
	    return;
	}
	oc = *poc;
	switch(node_type(Node_of(oc->curve,oc->orient))) 
	{
#if defined(FULL_PHYSICS)
	case CROSS_NODE:
	case OVERTAKE_NODE:
	case WAVE_END_NODE:
#endif /* defined(FULL_PHYSICS) */
	case MACH_NODE:
	case ATTACHED_B_NODE:
	case DIRICHLET_NODE:
	case NEUMANN_NODE:
	    delete_oc_curve_from_family(poc,cfamily);
	    break;
	case PASSIVE_NODE:
	case FIXED_NODE:
	case B_REFLECT_NODE:
	    screen("ERROR in substitute_contact_continuation(), "
	           "invalid node type\n");
	    clean_up(ERROR);
	    break;
#if defined(FULL_PHYSICS)
	case DIFFRACTION_NODE:
	case TRANSMISSION_NODE:
	    for (c = Node_of(oc->curve,oc->orient)->in_curves; c && *c; c++)
	    {
	    	if (is_scalar_wave(wave_type(*c)) &&
		        ((oc->orient == POSITIVE_ORIENTATION) || 
		    (oc->curve != *c)))
	    	{
	    	    oc->curve = *c;
	    	    oc->orient = NEGATIVE_ORIENTATION;
	    	    DEBUG_LEAVE(substitute_contact_continuation)
	    	    return;
	    	}
	    }
	    for (c = Node_of(oc->curve,oc->orient)->out_curves; c && *c; c++)
	    {
	    	if (is_scalar_wave(wave_type(*c)) &&
	    	        ((oc->orient == NEGATIVE_ORIENTATION) || 
	    	     (oc->curve != *c)))
	    	{
	    	    oc->curve = *c;
	    	    oc->orient = POSITIVE_ORIENTATION;
	    	    DEBUG_LEAVE(substitute_contact_continuation)
	    	    return;
	    	}
	    }
	    break;
	case CC_NODE:
	    /* Nothing is done for CC_NODE in 
	    	    substitute_contact_continuation */
	    break;
#endif /* defined(FULL_PHYSICS) */
	}
	DEBUG_LEAVE(substitute_contact_continuation)
}		/*end substitute_contact_continuation*/


/*
*		substitute_shock_continuation():
*/

LOCAL void substitute_shock_continuation(
	O_CURVE		**poc,
	O_CURVE_FAMILY	**cfamily)
{
	O_CURVE		*oc;

	DEBUG_ENTER(substitute_shock_continuation)

	if (poc == NULL)
	{
	    DEBUG_LEAVE(substitute_shock_continuation)
	    return;
	}
	oc = *poc;

	switch(node_type(Node_of(oc->curve,oc->orient))) 
	{
	case DIRICHLET_NODE:
	case NEUMANN_NODE:
	case B_REFLECT_NODE:
	case ATTACHED_B_NODE:
#if defined(FULL_PHYSICS)
	case WAVE_END_NODE:
#endif /* defined(FULL_PHYSICS) */
	    delete_oc_curve_from_family(poc,cfamily);
	    break;
	case PASSIVE_NODE:
	case FIXED_NODE:
#if defined(FULL_PHYSICS)
	case CC_NODE:
#endif /* defined(FULL_PHYSICS) */
	    screen("ERROR in substitute_shock_continuation(), ");
	    screen("invalid node type\n");
	    clean_up(ERROR);
	    break;
#if defined(FULL_PHYSICS)
	case DIFFRACTION_NODE:
	case TRANSMISSION_NODE:
	    delete_oc_curve_from_family(poc,cfamily);
	    break;
	case CROSS_NODE:
	case OVERTAKE_NODE:
	    delete_oc_curve_from_family(poc,cfamily);
	    break;
#endif /* defined(FULL_PHYSICS) */
	case MACH_NODE:
	    break;
	default:
	    screen("ERROR in substitute_shock_continuation(), ");
	    screen("invalid node type\n");
	    clean_up(ERROR);
	    break;
	}
	DEBUG_LEAVE(substitute_shock_continuation)
}		/*end substitute_shock_continuation*/


/*
*			incident_curves_exit_at_bdry():
*
*                                      |
*	     ------------            A *
*	    /            \             |\      
*	---*--------------*---         | \     
*          A              B            |  \    
*                                      |   \      
*				       *----*---
*                                           B
*
*	This function deals with the case were an incident shock tries to
*	pass through the boundary, so that nodes A and B interact (left
*	picture).  On the right is a similar case, but the curve is passing
*	through a corner, so there is a fixed node between A and B.  It is
*	possible to have reflected waves at A and B (on a Neumann boundary),
*	which then interact (case 2 below).  On a Dirichlet boundary, there
*	are no reflected waves (case 0 below).  In both cases, the null
*	curves are deleted, at which time the states behind the exiting
*	curve are copied onto the boundary.  In the case on the right, the
*	process can be thought of as extending the behind boundaries to the
*	fixed node, so that the states at the node come from behind the 
*	physical curve.
*/

LOCAL int incident_curves_exit_at_bdry(
	Front		*front,
	Front		*newfront,
	Wave		*wave,
	RPROBLEM	*rp)
{
	NODE		*newn;
	CURVE		*joined_c;
	O_CURVE		*newoc[2], *oldoc[2];
	COMPONENT	comp[3][2];
	ANGLE_DIRECTION	i_to_a_dir[2];
	int		i;
	int		status = ERROR_IN_STEP;

	DEBUG_ENTER(incident_curves_exit_at_bdry)

	switch (rp_num_refl(rp))
	{
	case 2:
	    if (find_curves_at_rp_with_status(oldoc,newoc,2,rp,REFLECTED) ==
							      FUNCTION_FAILED)
	    {
	       (void) printf("WARNING in incident_curves_exit_at_bdry(), "
	                     "can't find reflected curves\n");
	       DEBUG_LEAVE(incident_curves_exit_at_bdry)
	       return ERROR_IN_STEP;
	    }
	    for (i = 0; i < 2; i++)
	    {
	    	CURVE	    *cinc;
	    	ORIENTATION c_orient;
	    	int	    l_to_r;

	    	find_curve_with_status(Node_of_o_curve(newoc[i]),
				       &cinc,&c_orient,INCIDENT);
	    	i_to_a_dir[i] = incident_shock_orientation(wave_type(cinc),
							   c_orient);
	    	l_to_r = curve_ang_oriented_l_to_r(i_to_a_dir[i],c_orient);
	    	comp[0][i] = (l_to_r) ? positive_component(cinc) :
	    				negative_component(cinc);
	    	comp[1][i] = (l_to_r) ? negative_component(cinc) :
	    				positive_component(cinc);
	    	l_to_r = curve_ang_oriented_l_to_r(i_to_a_dir[i],
					           newoc[i]->orient);
	    	comp[2][i] = (l_to_r) ? negative_component(newoc[i]->curve) :
					positive_component(newoc[i]->curve);
	    }
	    for (i = 1; i < 3; i++)
	    {
	    	int jmin, jmax;

	    	if (comp[i][0] == comp[i][1]) continue;
	    	jmin = (comp[i][0] < comp[i][1]) ? 0 : 1;
	    	jmax = (jmin+1)%2;
	    	reset_component_of_loop(newoc[jmax]->curve,
				        newoc[jmax]->orient,i_to_a_dir[jmax],
				        comp[i][jmin],front);
	    }


	    delete_null_physical_curves(rp);
	    status = join_propagated_curves(&newn,&joined_c,newoc,oldoc,NO,
			                    front,newfront,(POINTER)wave,rp);
	    break;
	case 0:
	    status = GOOD_STEP;
	    break;
	default:
	       print_rproblem(rp);
	       (void) printf("WARNING in incident_curves_exit_at_bdry(), "
	                     "code needed for other than 2 reflected curves\n");
	       DEBUG_LEAVE(incident_curves_exit_at_bdry)
	       return ERROR_IN_STEP;
	}

	delete_null_physical_curves(rp);
	delete_null_boundary_curves(rp,front,(POINTER)wave);

	set_rp_statistics(rp);
	DEBUG_LEAVE(incident_curves_exit_at_bdry)
	return status;
}		/*end incident_curves_exit_at_bdry*/


LOCAL	boolean is_refl_curve_crosses_corner(
	RPROBLEM	*rp)
{

	DEBUG_ENTER(is_refl_curve_crosses_corner)
	if (rp == NULL)
	{
	    DEBUG_LEAVE(is_refl_curve_crosses_corner)
	    return NO;
	}

	if (rp_num_incident(rp) == 0 && rp_num_refl(rp) == 1 &&
			rp->num_nod == 3 && 
			rp->num_bdry_nod == 3 && rp->num_fxd == 1)
	{
	    DEBUG_LEAVE(is_refl_curve_crosses_corner)
	    return YES;
	}

	DEBUG_LEAVE(is_refl_curve_crosses_corner)
	return NO;
}		/*end is_refl_curves_crosses_corner*/

LOCAL	int  refl_curve_crosses_corner(
	Front		*front,
	Front		*newfront,
	Wave		*wave,
	RPROBLEM	*rp)
{
	CROSS		*cross;
	RP_NODE		*rpn;
	RP_NODE		*newrpn;
	O_CURVE		*ca;
	CURVE		*nbc;
	Locstate	st, sr, sl;
	ORIENTATION	nbc_orient;
	ANGLE_DIRECTION	i_to_f_dir;
	SIDE		int_side;
	O_CURVE		Cp, *cp = &Cp;
	O_CURVE		OldCp, *oldcp = &OldCp;

	DEBUG_ENTER(refl_curve_crosses_corner)

	delete_null_physical_curves(rp);
	delete_null_boundary_curves(rp,newfront,(POINTER)wave);

	find_rpn_with_physical_node(&rpn,rp,NO);
	find_curve_with_status(rpn->node,&cp->curve,&cp->orient,REFLECTED);
	Check_return(
	    find_correspond_of_oriented_curve(cp,oldcp,rpn->old_node,
					      front,rp->old_intfc),
	    refl_curve_crosses_corner)
	i_to_f_dir = incident_shock_orientation(wave_type(cp->curve),
						cp->orient);
	ca = (i_to_f_dir == COUNTER_CLOCK) ? rp->ang_ordered_curves->last :
			                     rp->ang_ordered_curves->first;

	if (prepare_curves_for_bdry_cross(cp,oldcp,ca,i_to_f_dir,&newrpn,front,
					  wave,rp) == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in refl_curve_crosses_corner(), "
	                  "can't find position past corner\n");
	    DEBUG_LEAVE(refl_curve_crosses_corner)
	    return ERROR_IN_STEP;
	}

	int_side = (curve_ang_oriented_l_to_r(i_to_f_dir,ca->orient)) ?
		   NEGATIVE_SIDE : POSITIVE_SIDE;
	if (!generate_boundary_cross_list(&cross,rp,front,(POINTER)wave))
	{
	    (void) printf("WARNING in refl_curve_crosses_corner(), "
	                  "generate_boundary_cross_list() failed\n");
	    DEBUG_LEAVE(refl_curve_crosses_corner)
	    return MODIFY_TIME_STEP;
	}
	install_bdry_cross(cross,cp,rpn,rp,front,&int_side);
	delete_null_boundary_curves(rp,newfront,(POINTER)wave);
	Check_return(next_boundary(ca->curve,ca->orient,&nbc,&nbc_orient),
		     refl_curve_crosses_corner)
	st = (curve_ang_oriented_l_to_r(i_to_f_dir,cp->orient)) ?
		Right_state_at_node_of_o_curve(cp) :
		Left_state_at_node_of_o_curve(cp);
	switch (wave_type(nbc))
	{
	case SUBDOMAIN_BOUNDARY:
	    sl = return_obst_state();
	    sr = return_obst_state();
	    break;
	case DIRICHLET_BOUNDARY:
	    sl = sr = st;
	    break;
	case NEUMANN_BOUNDARY:
	    if (int_side == POSITIVE_SIDE)
	    {
	    	sr = st;
	    	sl = return_obst_state();
	    }
	    else
	    {
	    	sl = st;
	    	sr = return_obst_state();
	    }
	    break;
	default:
	    (void) printf("WARNING in refl_curve_crosses_corner(), "
	                  "unknown wave type\n");
	    DEBUG_LEAVE(refl_curve_crosses_corner)
	    return ERROR_IN_STEP;
	}
	assign_interacting_states(Node_of(nbc,nbc_orient)->posn,
	    		          nbc,Opposite_orient(nbc_orient),front,sl,sr);

	DEBUG_LEAVE(refl_curve_crosses_corner)
	return GOOD_STEP;
}		/*end refl_curve_crosses_corner*/

LOCAL	boolean prepare_curves_for_bdry_cross(
	O_CURVE		*cp,
	O_CURVE		*oldcp,
	O_CURVE		*ca,
	ANGLE_DIRECTION	i_to_f_dir,
	RP_NODE		**newrpn,
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp)
{
	NODE		*newn;
	INTERFACE	*intfc = rp->new_intfc;
	RPROBLEM	*newrp = NULL;
	O_CURVE		*oldca;
	BOND		*pbcr, *abcr;
	Locstate	l_st, r_st;
	double		sp, sa;
	double		dt_frac;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	NODE_FLAG	flag;
	boolean		c_ext[2];
	POINT		*pc;

	DEBUG_ENTER(prepare_curves_for_bdry_cross)

	set_to_next_node_only(flag);
	pc = Point(Coords(Node_of_o_curve(cp)->posn));
	oldca = (i_to_f_dir == COUNTER_CLOCK) ?
				rp->old_ang_ordered_curves->last :
				rp->old_ang_ordered_curves->first;
	set_status_at_node(cp->curve,cp->orient,INCIDENT);

	if (!cross_or_extend_to_cross_two_propagated_curves(oldcp,cp,oldca,
							       ca,&pc,&pbcr,
							       &abcr,&sp,&sa,
							       front,
							       (POINTER)wave,
			                                       &newrp,rp->dt,
							       &dt_frac,flag,
							       c_ext))
	{
	    DEBUG_LEAVE(prepare_curves_for_bdry_cross)
	    return FUNCTION_FAILED;
	}
	newn = make_node(pc);
	*newrpn = add_to_rp_node_list(rp,newn,NULL);
	change_node_of_curve(cp->curve,cp->orient,newn);
	switch (wave_type(ca->curve))
	{
	case NEUMANN_BOUNDARY:
	    node_type(newn) = NEUMANN_NODE;
	    break;
	case DIRICHLET_BOUNDARY:
	    node_type(newn) = DIRICHLET_NODE;
	    break;
	default:
	    (void) printf("WARNING in prepare_curves_for_bdry_cross(), "
	                  "invalid wave type\n");
	    DEBUG_LEAVE(prepare_curves_for_bdry_cross)
	    return FUNCTION_FAILED;
	}
	set_status_at_node(cp->curve,cp->orient,INCIDENT);

	interpolate_intfc_states(intfc) = YES;
	l_st = left_state(pc);	r_st = right_state(pc);
	left_state_along_bond(sp,pbcr,cp->curve,l_st);
	right_state_along_bond(sp,pbcr,cp->curve,r_st);
	cut_curve(pc,pbcr,cp->curve,cp->orient,front,l_st,r_st);

	if (node_type(newn) == NEUMANN_NODE)
	{
	    POINT    *pt;
	    RP_NODE  *rpn;
	    Locstate  sl, sr, st;
	    double     t[MAXD], n[MAXD], v[MAXD];
	    int	      i, dim = front->rect_grid->dim;
	    size_t    sizest = front->sizest;

	    alloc_state(front->interf,&st,sizest);
	    alloc_state(front->interf,&sl,sizest);
	    alloc_state(front->interf,&sr,sizest);
	    find_tangent_to_curve(Node_of_o_curve(oldca)->posn,
	    		          Bond_at_node_of_o_curve(oldca),
				  oldca->curve,oldca->orient,t,front);
	    normal(Node_of_o_curve(oldca)->posn,
		   Hyper_surf_element(Bond_at_node_of_o_curve(oldca)),
		   Hyper_surf(oldca->curve),n,front);
	    if (is_obstacle_state(Right_state_at_node_of_o_curve(oldca)))
	    {
	    	for (i = 0; i < dim; i++) n[i] *= -1.0;
	    }
	    if (curve_ang_oriented_l_to_r(i_to_f_dir,oldcp->orient))
	    {
	    	ft_assign(st,Right_state_at_node_of_o_curve(oldcp),sizest);
	    }
	    else
	    {
	    	ft_assign(st,Left_state_at_node_of_o_curve(oldcp),sizest);
	    }
	    zero_normal_velocity(st,n,dim);
	    for (i = 0; i < dim; i++) v[i] = -vel(i,st);
	    ft_assign(sr,st,sizest);	ft_assign(sl,st,sizest);
	    if (curve_ang_oriented_l_to_r(i_to_f_dir,cp->orient))
	    {
	    	zero_state_velocity(sl,dim);
	    	add_velocity_to_state(sl,v);
	    }
	    else
	    {
	    	zero_state_velocity(sr,dim);
	    	add_velocity_to_state(sr,v);
	    }
	    for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	    {
	    	if (is_fixed_node(rpn->old_node))
		    break;
	    }
	    attach_phys_curve_to_n_bdry(cp->curve,
	    			        Opposite_orient(cp->orient),
					wave_type(cp->curve),&pt,sl,sr,
					Coords(rpn->old_node->posn),
					t,n,front,0.0,rp->dt,front);
	    ft_assign(Left_state_at_node_of_o_curve(cp),left_state(pt),sizest);
	    ft_assign(Right_state_at_node_of_o_curve(cp),right_state(pt),sizest);

	    interpolate_state_next_to_node(cp->curve,cp->orient,front);


	    free_these(3,st,sl,sr);
	}

	interpolate_intfc_states(intfc) = sav_intrp;
	DEBUG_LEAVE(prepare_curves_for_bdry_cross)
	return FUNCTION_SUCCEEDED;
}		/*end prepare_curves_for_bdry_cross*/


/*
*			g_interior_rp():
*
*	Solves Riemann problems which result from interior interactions
*/

LOCAL	int g_interior_rp(
	Front		*front,
	Front		*newfront,
	Wave		*wave,
	RPROBLEM	*rp)
{
	int		status = ERROR_IN_STEP;

	DEBUG_ENTER(g_interior_rp)
	if (DEBUG)
	{
	    (void) printf("interior type\n");
	    (void) printf("Rproblem into g_interior_rp()\n");
	    print_rproblem(rp);
	}
	
#if defined(FULL_PHYSICS)
	if (is_diffracted_shock_exits_contact(rp) == YES)
	    status = diffracted_shock_exits_contact(front,newfront,wave,rp);
	else if (is_overtake_diffraction_collision(rp) == YES)
	    status = overtake_diffraction_collision(front,newfront,wave,rp);
	else if (rp->num_ang_ordered_curves == 2)
	    status = connect_interacting_curves(front,wave,rp);
	else if (is_overtake_overtake_collision(rp) == YES)
	    status = overtake_overtake_collision(front,newfront,wave,rp);
	else if (is_deprecursor(rp) == YES)
	    status = deprecurse_shock_refraction(front,wave,rp);
	else
	{
	    status = MODIFY_TIME_STEP;
	    (void) printf("WARNING in g_interior_rp(), "
	                  "code needed for g_2drproblem\n"
	                  "Reducing time step\n");
	    if (DEBUG) print_rproblem(rp);
	}
#endif /* defined(FULL_PHYSICS) */
	if (status == ERROR_IN_STEP)
	{
	    (void) printf("interior type rp\n");
	    (void) printf("WARNING in g_interior_rp(), "
			  "unable to resolve interior rp\n");
	    print_rproblem(rp);
	}
	DEBUG_LEAVE(g_interior_rp)
	return status;
}		/*end g_interior_rp*/

#if defined(FULL_PHYSICS)
LOCAL	boolean is_deprecursor(
	RPROBLEM	*rp)
{
	RP_NODE		*rpn;
	int		is_irnode, is_tnode, is_cnode;

	DEBUG_ENTER(is_deprecursor)
	is_irnode = is_tnode = is_cnode = NO;
	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (!rpn->old_node)
	    {
	    	DEBUG_LEAVE(is_deprecursor)
	    	return NO;
	    }
	    if (node_type(rpn->old_node) == TRANSMISSION_NODE) 
	    	is_tnode = YES;
	    if (node_type(rpn->old_node) == TOT_INT_REFL_NODE) 
	    	is_irnode = YES;
	    if (node_type(rpn->old_node) == CROSS_NODE) 
	    	is_cnode = YES;
	}
	DEBUG_LEAVE(is_deprecursor)
	return (is_tnode && is_irnode && is_cnode) ? YES : NO;
}		/*end is_deprecursor*/

LOCAL	boolean is_diffracted_shock_exits_contact(
	RPROBLEM	*rp)
{
	RP_NODE		*rpn;

	DEBUG_ENTER(is_diffracted_shock_exits_contact)

	if (rp->num_nod != 2)
	{
	    DEBUG_LEAVE(is_diffracted_shock_exits_contact)
	    return NO;
	}

	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (!rpn->old_node)
	    {
	    	DEBUG_LEAVE(is_diffracted_shock_exits_contact)
	    	return NO;
	    }
	    if (node_type(rpn->old_node) != DIFFRACTION_NODE)
	    {
	    	DEBUG_LEAVE(is_diffracted_shock_exits_contact)
	    	return NO;
	    }
	}
	DEBUG_LEAVE(is_diffracted_shock_exits_contact)
	return YES;
}		/*end is_diffracted_shock_exits_contact*/

LOCAL	boolean is_n_diffracted_shock_exits_contact(
	RPROBLEM	*rp)
{
	RP_NODE		*rpn;

	DEBUG_ENTER(is_n_diffracted_shock_exits_contact)

	if (rp->num_phys != 1 || rp->bdry_type1 != NEUMANN_BOUNDARY || 
		    		 rp->bdry_type2 != NEUMANN_BOUNDARY)
	{
	    DEBUG_LEAVE(is_n_diffracted_shock_exits_contact)
	    return NO;
	}
	if (rp->num_nod != 3 || rp_num_incident(rp) != 0)
	{
	    DEBUG_LEAVE(is_n_diffracted_shock_exits_contact)
	    return NO;
	}

	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (!rpn->old_node)
	    {
	    	DEBUG_LEAVE(is_n_diffracted_shock_exits_contact)
	    	return NO;
	    }
	    if (node_type(rpn->old_node) == NEUMANN_NODE) continue;
	    if (node_type(rpn->old_node) != DIFFRACTION_NODE)
	    {
	    	DEBUG_LEAVE(is_n_diffracted_shock_exits_contact)
	    	return NO;
	    }
	}
	DEBUG_LEAVE(is_n_diffracted_shock_exits_contact)
	return YES;
}		/*end is_n_diffracted_shock_exits_contact*/

LOCAL	boolean is_overtake_diffraction_collision(
	RPROBLEM	*rp)
{
	RP_NODE		*rpn, *rpn1;
	int		nt, nt1;

	DEBUG_ENTER(is_overtake_diffraction_collision)

	if (rp->num_nod != 2)
	{
	    DEBUG_LEAVE(is_overtake_diffraction_collision)
	    return NO;
	}

	rpn = rp->first_rp_node;	rpn1 = rpn->next;
	if (rpn->old_node == NULL || rpn1->old_node == NULL)
	{
	    DEBUG_LEAVE(is_overtake_diffraction_collision)
	    return NO;
	}
	nt = node_type(rpn->old_node);	nt1 = node_type(rpn1->old_node);
	if ((nt == DIFFRACTION_NODE && nt1 == OVERTAKE_NODE) ||
	    (nt == OVERTAKE_NODE && nt1 == DIFFRACTION_NODE))
	{
	    DEBUG_LEAVE(is_overtake_diffraction_collision)
	    return YES;
	}
	DEBUG_LEAVE(is_overtake_diffraction_collision)
	return NO;
}		/*end is_overtake_diffraction_collision*/

LOCAL	boolean is_overtake_overtake_collision(
	RPROBLEM	*rp)
{
	RP_NODE		*rpn;

	DEBUG_ENTER(is_overtake_overtake_collision)

	if (rp->num_nod != 2)
	{
	    DEBUG_LEAVE(is_overtake_overtake_collision)
	    return NO;
	}

	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (!rpn->old_node)
	    {
	    	DEBUG_LEAVE(is_overtake_overtake_collision)
	    	return NO;
	    }
	    if (node_type(rpn->old_node) != OVERTAKE_NODE)
	    {
	    	DEBUG_LEAVE(is_overtake_overtake_collision)
	    	return NO;
	    }
	}
	DEBUG_LEAVE(is_overtake_overtake_collision)
	return YES;
}		/*end is_overtake_overtake_collision*/


/*
*			deprecurse_shock_refraction():
*
*	This function is responsible for the change in topology from
*	a precursor configuration back to a single diffraction node.  This
*	involves deleting the cross, overtake and transmission nodes, and
*	the relevant curves, and propagating the diffraction node.
*/

LOCAL	int deprecurse_shock_refraction(
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp)
{
	CURVE		*new_trans_c, *c;
	NODE		*newn = NULL;
	O_CURVE		*oc;
	RP_NODE		*rpn;
	RPROBLEM	*newrp = NULL;
	double		dt_frac;
	NODE_FLAG 	flag;
	ORIENTATION	ntc_or, orient;
	int		status;

	DEBUG_ENTER(deprecurse_shock_refraction)

	/* Note: we don't correct angles to avoid tangles.  After the new
	   configuration is installed, it doesn't matter because we can
	   reduce the timestep. */
	clear_node_flag(flag);
	to_next_node_only(flag) = YES;
	use_subsonic_state(flag) = YES;
	node_vel_by_angle(flag) = YES;
	dont_correct_angles_at_node(flag) = YES;

	status = is_premature_deprecursion(front,wave,rp);

	if (status == GOOD_NODE)
	{
	    if (DEBUG)
	    	(void) printf("premature precursion detected\n");
	    DEBUG_LEAVE(deprecurse_shock_refraction)
	    return GOOD_STEP;
	}
	else if (status == MODIFY_TIME_STEP_NODE)
	{
	    (void) printf("WARNING in deprecurse_shock_refraction(), "
	                  "is_premature_deprecursion() returns ");
	    print_node_status("status = ",status,"\n");
	    DEBUG_LEAVE(deprecurse_shock_refraction)
	    return MODIFY_TIME_STEP;
	}
	else if (untrack_node(front,PRECURSOR_RR_DIFFRACTION) == YES)
	{
	    status = untrack_precursor(front,wave,front->dt,rp);
	    DEBUG_LEAVE(deprecurse_shock_refraction)
	    return status;
	}


	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (node_type(rpn->node) == TRANSMISSION_NODE) 
	    {
	    	newn = rpn->node;
	    	break;
	    }
	}
	if (newn == NULL)
	{
	    (void) printf("WARNING in deprecurse_shock_refraction(), "
	                  "transmission node not found\n");
	    DEBUG_LEAVE(deprecurse_shock_refraction)
	    return ERROR_IN_STEP;
	}
	delete_null_physical_curves(rp);

	    /* Find incident at old transmission node */
	find_curve_with_status(newn,&new_trans_c,&ntc_or,INCIDENT);
	if (new_trans_c == NULL)
	{
	    (void) printf("WARNING in deprecurse_shock_refraction(), "
	                  "Unable to find new transmitted curve\n");
	    DEBUG_LEAVE(deprecurse_shock_refraction)
	    return ERROR_IN_STEP;
	}

	for (oc = rp->ang_ordered_curves->first; oc != NULL; oc = oc->next)
	{
	    if (oc->curve && (Node_of_o_curve(oc) != newn))
	    	change_node_of_curve(oc->curve,oc->orient,newn);
	    if (oc == rp->ang_ordered_curves->last) break;
	}
	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (rpn->node == newn) continue;
	    if (delete_node(rpn->node) == FUNCTION_FAILED)
	    {
	        (void) printf("WARNING in deprecurse_shock_refraction(), "
	                      "Unable to delete node\n");
	        print_node(rpn->node);
	        DEBUG_LEAVE(deprecurse_shock_refraction)
	        return ERROR_IN_STEP;
	    }
	    rpn->node = NULL;
	}

	/* Reset status of transmitted wave at old overtake node */
	find_curve_with_status(newn,&c,&orient,TRANSMITTED);
	if (c != NULL) set_status_at_node(c,orient,REFLECTED);

	    /* Reset status of incident at old transmission node */
	set_status_at_node(new_trans_c,ntc_or,TRANSMITTED);

	node_type(newn) = DIFFRACTION_NODE;
	if (DEBUG)
	    (void) printf("Calling node_propagate\n");
	dt_frac = 1.0;
	status = (*front->node_propagate)(front,(POINTER)wave,NULL,newn,&newrp,
	    			          rp->dt,&dt_frac,flag,NULL);
	if (status != GOOD_NODE) 
	{
	    /* Try again */
	    clear_node_flag(flag);
	    to_next_node_only(flag) = YES;
	    use_subsonic_state(flag) = YES;
	    free_rp(newrp);
	    newrp = NULL;
	    dt_frac = 1.0;
	    status = (*front->node_propagate)(front,(POINTER)wave,NULL,
					      newn,&newrp,rp->dt,&dt_frac,
					      flag,NULL);
	    if (status != GOOD_NODE) 
	    {
	        free_rp(newrp);
	        (void) printf("WARNING in deprecurse_shock_refraction(), "
	                      "node propagation failed\n");
	        DEBUG_LEAVE(deprecurse_shock_refraction)
	        return MODIFY_TIME_STEP;
	    }
	}
	DEBUG_LEAVE(deprecurse_shock_refraction)
	return GOOD_STEP;

}		/*end deprecurse_shock_refraction*/


/*
*			untrack_precursor():
*
*	This function provides an untracking algorithm for a precursor for
*	use in the case where the incident angle continues to increase, and
*	the cross node should eventually bifurcate to a double Mach
*	configuration.  The cross node is identified, and drives the untracking.
*	The behind curves for the cross node are identified and untracked.
*	The two incidents are given the same wave type, and will be joined
*	automatically by the untracking function.  All that will be left is
*	the transmission node.
*/

LOCAL int untrack_precursor(
	Front		*front,
	Wave		*wave,
	double		dt,
	RPROBLEM	*rp)
{
	ANGLE_DIRECTION	i0_to_i4_dir;
	COMPONENT	newcomp;
	NODE		*cnode = NULL;
	NODE		*tnode = NULL, *old_tnode = NULL;
	RP_NODE		*rpn;
	RPROBLEM	*newrp = NULL;
	UNTRACK_FLAG	uflag;
	NODE_FLAG	flag;
	double		dt_frac;
	int		i;
	int		status;
	static O_CURVE	*oldoc = NULL;
	static O_CURVE	*oc[5];

	DEBUG_ENTER(untrack_precursor)

	if (oldoc == NULL)
	{
	    for (i = 0; i < 5; i++)
	    	scalar(&oc[i],sizeof(O_CURVE));
	    scalar(&oldoc,sizeof(O_CURVE));
	}

	for (i = 0; i < 5; i++)
	    oc[i]->curve = NULL;

	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (node_type(rpn->node) == CROSS_NODE)
	    {
	    	cnode = rpn->node;
	    	break;
	    }
	    if (node_type(rpn->node) == TRANSMISSION_NODE)
	    {
	    	tnode = rpn->node;
	    	old_tnode = rpn->old_node;
	    	break;
	    }
	}

	if ((cnode == NULL) || (tnode == NULL))
	{
	    DEBUG_LEAVE(untrack_precursor)
	    return ERROR_IN_STEP;
	}

	identify_curves_with_status(cnode,oc[0],oc[4],INCIDENT);
	find_curve_with_status(cnode,&oc[2]->curve,&oc[2]->orient,SLIP);
	identify_curves_with_status(cnode,oc[1],oc[3],REFLECTED);

	if ((oc[0]->orient == POSITIVE_ORIENTATION &&
	    	is_forward_wave(wave_type(oc[0]->curve)))
	    	    ||
	    (oc[0]->orient == NEGATIVE_ORIENTATION &&
	    	is_backward_wave(wave_type(oc[0]->curve))))
	    i0_to_i4_dir = CLOCKWISE;
	else
	    i0_to_i4_dir = COUNTER_CLOCK;

	newcomp = (curve_ang_oriented_l_to_r(i0_to_i4_dir,oc[0]->orient)) ?
	    	negative_component(oc[0]->curve) :
	    	positive_component(oc[0]->curve);

	if (wave_type(oc[0]->curve) != wave_type(oc[4]->curve))
	    invert_curve(oc[4]->curve);

	for (i = 1; i < 4; i++)
	{
	    boolean       oppn_ss;
	    NODE	  *oppn;
	    /* Make sure curve was tracked initially, and not untracked
	     * recursively earlier in the for loop.
	     */
	    if ((oc[i]->curve == NULL) ||
	        (oc[i]->curve->interface == NULL))
	    	continue;
	    oppn = Opp_node_of_o_curve(oc[i]);
	    oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	    set_untrack_flag(uflag,oc[i]->orient,YES,oppn_ss,YES,YES,YES);
	    find_corr_cur_in_rp(oc[i],oldoc,front,rp);
	    (void) untrack_curve(oc[i],oldoc,newcomp,dt,front,
				 (POINTER)wave,rp,uflag);
	}
	dt_frac = 1.0;
	set_to_next_node_only(flag);
	status = (*front->node_propagate)(front,(POINTER)wave,old_tnode,
					  tnode,&newrp,rp->dt,&dt_frac,
					  flag,NULL);
	if (status != GOOD_NODE)
	{
	    free_rp(newrp);
	    (void) printf("WARNING in untrack_precursor(), "
	                  "node propagation failed\n");
	    DEBUG_LEAVE(untrack_precursor)
	    return MODIFY_TIME_STEP;
	}

	DEBUG_LEAVE(untrack_precursor)
	return GOOD_STEP;
}		/*end untrack_precursor*/


/*
*			is_premature_deprecursion():
*
*	This function checks to make sure we actually want to deprecurse.
*	Since we have a number of nodes very close together, the propagation
*	often fails merely because the nodes were not propagated in the
*	correct order.  This function attempts such an ordering, and only
*	when this fails do we signal a TRUE deprecursion.  We allow for the
*	case of an untracked overtake node.
*/

LOCAL	int is_premature_deprecursion(
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp)
{
	NODE_FLAG flag;
	int	  status;

	DEBUG_ENTER(is_premature_deprecursion)

	clear_node_flag(flag);
	to_next_node_only(flag) = YES;
	set_virtuals_by_adjacent_bond(flag) = YES;
	status = repropagate_node(front,wave,rp,TRANSMISSION_NODE,flag);
	if (status != GOOD_NODE)
	{
	    if (DEBUG)
		print_node_status("status = ",status,"\n");
	    DEBUG_LEAVE(is_premature_deprecursion)
	    return status;
	}
	set_virtuals_by_adjacent_bond(flag) = NO;

	status = repropagate_node(front,wave,rp,CROSS_NODE,flag);
	if (status != GOOD_NODE)
	{
	    if (DEBUG)
		print_node_status("status = ",status,"\n");
	    DEBUG_LEAVE(is_premature_deprecursion)
	    return status;
	}

	status = repropagate_node(front,wave,rp,TOT_INT_REFL_NODE,flag);
	if (status != GOOD_NODE)
	{
	    if (DEBUG)
		print_node_status("status = ",status,"\n");
	    DEBUG_LEAVE(is_premature_deprecursion)
	    return status;
	}

	status = repropagate_node(front,wave,rp,OVERTAKE_NODE,flag);
	if (status == NO_STORAGE_NODE)
	    status = GOOD_NODE;		/* overtake node not tracked */

	if (DEBUG)
	    print_node_status("status = ",status,"\n");
	DEBUG_LEAVE(is_premature_deprecursion)
	return status;
}		/*end is_premature_deprecursion*/


LOCAL int repropagate_node(
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp,
	int		n_type,
	NODE_FLAG	flag)
{
	RP_NODE		*rpn;
	NODE		*newn, *oldn;
	RPROBLEM	*newrp = NULL;
	double		dt_frac = 1.0;
	int		status;

	DEBUG_ENTER(repropagate_node)

	newn = oldn = NULL;
	for (rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    if (node_type(rpn->node) == n_type) 
	    {
	    	newn = rpn->node;	oldn = rpn->old_node;
	    	break;
	    }
	}
	if (newn == NULL)
	{
	    (void) printf("WARNING in repropagate_node(), ");
	    print_node_type("no such node, ",n_type,"\n",front->interf);
	    DEBUG_LEAVE(repropagate_node)
	    return NO_STORAGE_NODE;
	}
	propagation_status(newn) = UNPROPAGATED_NODE;
	status = (*front->node_propagate)(front,(POINTER)wave,oldn,newn,
					  &newrp,rp->dt,&dt_frac,flag,NULL);

	if (status != GOOD_NODE) free_rp_list(&newrp);

	DEBUG_LEAVE(repropagate_node)
	return status;
}		/*end repropagate_node*/

/*
*			overtake_overtake_collision():
*
*	This function handles the case of one wave completely overtaking the
*	another.  Thus the INCIDENT and OVERTOOK curves disappear, leaving
*	only the reflected, contact and transmitted waves (if tracked).  We
*	do not handle the case where the overtake reversed direction, and the
*	OVERTOOK curve pulls ahead of the INCIDENT, thus causing disappearance
*	of the overtake configuration.  We assume that both incidents and the
*	transmitted waves are tracked.
*/

LOCAL	int overtake_overtake_collision(
	Front		*front,
	Front		*newfront,
	Wave		*wave,
	RPROBLEM	*rp)
{
	COMPONENT	comp[7][2];
	CURVE		*joined_c;
	NODE		*new_node;
	O_CURVE		*newc[7][2], *oldc[7][2];
	O_CURVE		*oc, *oldoc;
	RP_NODE		*rpn[2];
	int		i, j, k;
	int		w_type, n_status;
	int		err = NO;
	ANGLE_DIRECTION	i0_to_i1_dir[2];
	int		jmin, jmax;
	int		closed_curve;
	int		status;

	DEBUG_ENTER(overtake_overtake_collision)

	rpn[0] = rp->first_rp_node;	rpn[1] = rp->last_rp_node;

	for (i = 0; i < 7; i++) 
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
	    case OVERTOOK:
	    	i = 0;
	    	break;
	    case INCIDENT:
	    	i = 1;
	    	break;
	    case REFLECTED:
	    	w_type = wave_type(oldoc->curve);
	    	i = (is_rarefaction_leading_edge(w_type)) ? 2 :
	    	    	(is_shock_wave(w_type)) ?  3 : 4;
	    	break;
	    case SLIP:
	    	i = 5;
	    	break;
	    case TRANSMITTED:
	    	i = 6;
	    	break;
	    default:
	    	continue;
	    }
	    newc[i][j] = oc;	oldc[i][j] = oldoc;
	    if (oc == rp->ang_ordered_curves->last ||
	    	    oldoc == rp->old_ang_ordered_curves->last)
	    	break;
	}

	for (j = 0; j < 2; j++)
	{
	    if (newc[0][j] == NULL) err = YES;
	    if (newc[1][j] == NULL) err = YES;
	    if (newc[6][j] == NULL) err = YES;
	}
	if (newc[1][0] != newc[1][1]) err = YES;

	if (err)
	{
	    if (DEBUG)
	    {
	        (void) printf("WARNING in overtake_overtake_collision(), "
	                      "unexpected topology in rproblem\n");
	        (void) printf("newc[0][0] = %p, newc[0][1] = %p\n",
	    	              (POINTER)newc[0][0],(POINTER)newc[0][1]);
	        (void) printf("newc[1][0] = %p, newc[1][1] = %p\n",
			      (POINTER)newc[1][0],(POINTER)newc[1][1]);
	        (void) printf("newc[6][0] = %p, newc[6][1] = %p\n",
				  (POINTER)newc[6][0],(POINTER)newc[6][1]);
	    }
	    DEBUG_LEAVE(overtake_overtake_collision)
	    return ERROR_IN_STEP;
	}

	for (j = 0; j < 2; j++)
	{
	    i0_to_i1_dir[j] = (((newc[0][j]->orient == POSITIVE_ORIENTATION) &&
	    	     is_forward_wave(wave_type(newc[0][j]->curve)))
			    		 ||
	    	    ((newc[0][j]->orient == NEGATIVE_ORIENTATION) &&
	    	     is_backward_wave(wave_type(newc[0][j]->curve)))) ?
	    	         COUNTER_CLOCK : CLOCKWISE;

	    if (curve_ang_oriented_l_to_r(i0_to_i1_dir[j],newc[0][j]->orient))
	    {
	        comp[0][j] = negative_component(newc[0][j]->curve);
	        comp[1][j] = positive_component(newc[0][j]->curve);
	    }
	    else
	    {
	        comp[0][j] = positive_component(newc[0][j]->curve);
	        comp[1][j] = negative_component(newc[0][j]->curve);
	    }
	}

	for (i = 1; i < 7; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
	    	comp[i+1][j] = (curve_ang_oriented_l_to_r(
	    	    	i0_to_i1_dir[j],newc[i][j]->orient)) ? 
	    	    positive_component(newc[i][j]->curve) : 
	    	    negative_component(newc[i][j]->curve);
	    }
	}

	/* Ensure consistent components of new loops */
	/* mismatch of comp[2][j] is assumed impossible */

	for (i = 3; i < 7; i++)
	{
	    if (newc[i-1][0] == NULL) continue;
	    if (comp[i][0] != comp[i][1])
	    {
	    	jmin = (comp[i][0] < comp[i][1]) ? 0 : 1;
	    	jmax = (jmin+1)%2;
	    	reset_component_of_loop(newc[i-1][jmax]->curve,
				        newc[i-1][jmax]->orient,
					i0_to_i1_dir[jmax],
				        comp[i][jmin],front);
	    	for (k = i; k < 6; k++)
	    	{
	    	    comp[k][jmax] = comp[k][jmin];
	    	    if (newc[k][jmin] != NULL) break;
	    	}
	    }
	}

	delete_null_physical_curves(rp);

	for (i = 2; i < 7; i++)
	{
	    if (newc[i][0] == NULL) continue;

	    /* Install new curve configurations */

	    closed_curve = (newc[i][0]->curve == newc[i][1]->curve) ? YES : NO;

	    status = join_propagated_curves(&new_node,&joined_c,newc[i],oldc[i],
					    closed_curve ? YES : NO,
					    front,newfront,(POINTER)wave,rp);

	    if (status != GOOD_STEP)
	    {
	        (void) printf("WARNING: in overtake_overtake_collision(), \n"
	                      "join_propagated_curves() failed for "
	                      "curves %d\n",i);
	        if (DEBUG)
	            (void) printf("status = %d\n",status);
		DEBUG_LEAVE(overtake_overtake_collision)
		return status;
	    }

	    /* Now check for possible bifurctions in transmitted shocks */

	    if (closed_curve)
	    {
		CURVE *c = newc[i][0]->curve;

		set_incident_status(c,POSITIVE_ORIENTATION);
		set_incident_status(c,NEGATIVE_ORIENTATION);
	    }
	}

	(void) delete_node(rpn[0]->node);	rpn[0]->node = NULL;
	(void) delete_node(rpn[1]->node);	rpn[1]->node = NULL;

	DEBUG_LEAVE(overtake_overtake_collision)
	return GOOD_STEP;
}		/*end overtake_overtake_collision*/



#endif /* defined(FULL_PHYSICS) */


LOCAL void g_set_rp_statistics(
	RPROBLEM	*rp)
{
	RP_NODE		*rp_node;
	O_CURVE		*oc;
	O_CURVE_FAMILY	*contacts = NULL;

	DEBUG_ENTER(g_user_set_rp_statistics)

	if (!rp)
	{
	    DEBUG_LEAVE(g_user_set_rp_statistics)
	    return;
	}
	f_set_rp_statistics(rp);
	rp->num_ang_ordered_curves = 0;
	for (oc = rp->ang_ordered_curves->first; oc != NULL; oc = oc->next)
	{
	    rp->num_ang_ordered_curves++;
	    if (oc == rp->ang_ordered_curves->last) break;
	}
	rp_num_incident(rp) = 0;
	rp_num_refl(rp) = 0;
	rp_num_transm(rp) = 0;
	rp_num_contact(rp) = 0;
	rp_num_inc_shock(rp) = 0;
	for (rp_node = rp->first_rp_node; rp_node; rp_node = rp_node->next) 
	{
	    if (rpn_inc_shock1(rp_node)) 
	    {
	    	rp_num_inc_shock(rp)++;
	    	rp_num_incident(rp)++;
	    }
	    if (rpn_inc_shock2(rp_node))
	    {
	    	rp_num_inc_shock(rp)++;
	    	rp_num_incident(rp)++;
	    }
	    if (rpn_reflected1(rp_node))
	    	rp_num_refl(rp)++;
	    if (rpn_reflected2(rp_node))
	    	rp_num_refl(rp)++;
	    if (rpn_transmitted(rp_node))
	    	rp_num_transm(rp)++;
	    if (rpn_contact1(rp_node))
	    {
	    	join_cfamilies(&contacts,rpn_contact1(rp_node));
	    	for (oc = rpn_contact1(rp_node)->first; oc; oc = oc->next)
	    	    if (status_at_node(oc->curve,oc->orient) == INCIDENT)
	    	    	rp_num_incident(rp)++;
	    }
	    if (rpn_contact2(rp_node))
	    {
	    	join_cfamilies(&contacts,rpn_contact2(rp_node));
	    	for (oc = rpn_contact2(rp_node)->first; oc; oc = oc->next)
	    	    if (status_at_node(oc->curve,oc->orient) == INCIDENT)
	    	    	rp_num_incident(rp)++;
	    }
	}

	/* Remove duplicates in contact cfamily */

	rp_num_contact(rp) = remove_oc_duplicates(&contacts);
	free_o_curve_family(contacts);

	if (DEBUG)
	    print_rp_statistical_data(rp);
	DEBUG_LEAVE(g_user_set_rp_statistics)
}		/*end g_user_set_rp_statistics*/

LOCAL	int	remove_oc_duplicates(
		O_CURVE_FAMILY	**ocf)
{
	int 	num_oc;
	O_CURVE	*oc, *oc1;
	if (ocf == NULL || *ocf == NULL)
		return 0;
start:
	for (oc = (*ocf)->first; oc != NULL; oc = oc->next)
	{
	    for (oc1 = oc->next; oc1; oc1 = oc1->next)
	    {
	    	if (oc1->curve == oc->curve)
	    	{
	    	    delete_oc_curve_from_family(&oc1,ocf);
	    	    if (*ocf == NULL)
	    	    	return 0;
	    	    goto start;
	    	}
	    }
	}

	for (num_oc = 0, oc = (*ocf)->first; oc; oc = oc->next)
	    num_oc++;
	return num_oc;

}		/*end remove_oc_duplicates*/

LOCAL	void print_rp_statistical_data(
	RPROBLEM	*rp)
{
	DEBUG_ENTER(print_rp_statistical_data)

	(void) printf("\nStatistical Data for rproblem %p\n",(POINTER)rp);
	(void) printf("number of incident curves = %d\n",rp_num_incident(rp));
	(void) printf("number of incident shocks = %d\n", rp_num_inc_shock(rp));
	(void) printf("number of contact curves = %d\n",rp_num_contact(rp));
	(void) printf("number of reflected curves = %d\n",rp_num_refl(rp));
	(void) printf("number of transmitted curves = %d\n",rp_num_transm(rp));
	(void) printf("number of Neumann and Dirichlet = %d\n",rp->num_nd);
	(void) printf("number of fixed nodes = %d\n",rp->num_fxd);
	(void) printf("number of source nodes = %d\n",rp->num_srce);
	(void) printf("number of boundary nodes = %d\n",rp->num_bdry_nod);
	(void) printf("number of physical nodes = %d\n",rp->num_phys);
	(void) printf("total number of nodes = %d\n",rp->num_nod);
	print_wave_type("bdry_type1 = ",rp->bdry_type1,"\n",rp->new_intfc);
	print_wave_type("bdry_type2 = ",rp->bdry_type2,"\n",rp->new_intfc);
	(void) printf("End of Statistical Data for rproblem %p\n\n",
		      (POINTER)rp);
	DEBUG_LEAVE(print_rp_statistical_data)
}		/*end print_rp_statistical_data*/

EXPORT 	double find_position_and_dt_of_intersection(
	CURVE		*c0,	/* interacting curve 0 */
	ORIENTATION	c0_or,	/* c0 orient */
	CURVE		*c1,	/* interacting curve 1 */
	ORIENTATION	c1_or,	/* c1 orient */
	double		*coords,/* calculated estimate of interaction coords */
	Front		*fr,
	double		dt)
{
	BOND		*b0, *b1;
	double		wsp0[MAXD],wsp1[MAXD],dp[MAXD],dv[MAXD];
	double		v0,v1,normdp;
	double		dtf;
	double		partial_dt;
	int		i,dim;
	POINT		*p0, *p1;
	static POINT	*newp0 = NULL, *newp1 = NULL;

	DEBUG_ENTER(find_position_and_dt_of_intersection)

	if (newp0 == NULL)
	{
	    newp0 = Static_point(fr->interf);
	    newp1 = Static_point(fr->interf);
	}

	p0 = Node_of(c0,c0_or)->posn;
	b0 = Bond_at_node(c0,c0_or);

	g_pt_prop_by_w_speed(fr,p0,newp0,
		             Hyper_surf_element(b0),Hyper_surf(c0),dt,wsp0);

	p1 = Node_of(c1,c1_or)->posn;
	b1 = Bond_at_node(c1,c1_or);

	g_pt_prop_by_w_speed(fr,p1,newp1,
		             Hyper_surf_element(b1),Hyper_surf(c1),dt,wsp1);

	dim = fr->rect_grid->dim;
	for (i = 0; i < dim; i++)
	{
	    dp[i] = Coords(p1)[i] - Coords(p0)[i];
	    dv[i] = wsp0[i] - wsp1[i];
	}

	normdp = scalar_product(dp,dp,dim);
	dtf = normdp/scalar_product(dv,dp,dim);
	if (dtf < 0.0) dtf = 0.0;
	if (dtf > dt) dtf = dt;
	partial_dt = dt - dtf;

	/* Estimate position of interaction.  Make sure the position
	*  lies on segment connecting p0 & p1 */

	v0 = scalar_product(wsp0,dp,dim)/normdp;
	v1 = scalar_product(wsp1,dp,dim)/normdp;

	for (i = 0; i < dim; i++)
	    coords[i] = 0.5*(Coords(p0)[i]+Coords(p1)[i]+(v0+v1)*dp[i]*dtf);

	DEBUG_LEAVE(find_position_and_dt_of_intersection)
	return partial_dt;
}		/*end find_position_and_dt_of_intersection*/

LOCAL	int connect_interacting_curves(
	Front		*front,
	Wave		*wave,
	RPROBLEM	*rp)
{
	double		x, y;
	double		V[MAXD];
	NODE		*nd[2];
	POINT		*newp;
	O_CURVE		*oc[2];
	O_CURVE		*old_oc[2];
	CURVE		*newc;
	COMPONENT	left, right;
	int		i;
	boolean		sav_interp;

	DEBUG_ENTER(connect_interacting_curves)

	oc[0] = rp->ang_ordered_curves->first;
	old_oc[0] = rp->old_ang_ordered_curves->first;
	oc[1] = rp->ang_ordered_curves->last;
	old_oc[1] = rp->old_ang_ordered_curves->last;

	if (DEBUG)
	{
	    (void) printf("Interacting curves\n");
	    (void) printf("oc[0],\t");	print_o_curve(oc[0]);
	    (void) printf("old_oc[0],\t");	print_o_curve(old_oc[0]);
	    (void) printf("oc[1],\t");	print_o_curve(oc[1]);
	    (void) printf("old_oc[1],\t");	print_o_curve(old_oc[1]);
	}
	x = y = 0.0;
	for (i = 0; i < 2; i++)
	{
	    nd[i] = Node_of_o_curve(oc[i]);
	    newp = nd[i]->posn;
	    point_propagate(front,(POINTER)wave,
			    Node_of_o_curve(old_oc[i])->posn,
			    newp,Bond_at_node_of_o_curve(old_oc[i]),
			    old_oc[i]->curve,rp->dt,V);
	    ft_assign(Left_state_at_node_of_o_curve(oc[i]),left_state(newp),
		   front->sizest);
	    obstacle_state(front->interf,left_state(newp),front->sizest);
	    ft_assign(Right_state_at_node_of_o_curve(oc[i]),
	    	   right_state(newp),front->sizest);
	    obstacle_state(front->interf,right_state(newp),front->sizest);
	    x += 0.5 * Coords(newp)[0];
	    y += 0.5 * Coords(newp)[1];
	}
	delete_null_physical_curves(rp);
	Coords(nd[0]->posn)[0] = x;
	Coords(nd[0]->posn)[1] = y;
	change_node_of_curve(oc[1]->curve,oc[1]->orient,nd[0]);
	if (delete_node(nd[1]) == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in conect_interacting_curves(), "
	                  "unexpected case,  delete_node() failed\n");
	    DEBUG_LEAVE(connect_interacting_curves)
	    return ERROR_IN_STEP;
	}
	sav_interp = interpolate_intfc_states(rp->new_intfc);
	interpolate_intfc_states(rp->new_intfc) = YES;
	if (oc[0]->orient == oc[1]->orient)
	{
	    invert_curve(oc[1]->curve);
	    oc[1]->orient = Opposite_orient(oc[1]->orient);
	}

	left = negative_component(oc[0]->curve);
	right = positive_component(oc[0]->curve);
	if (oc[0]->orient == NEGATIVE_ORIENTATION)
	{
	    newc = join_curves(oc[0]->curve,oc[1]->curve,left,right,NULL);
	    roclists_after_join(rp,oc[0]->curve,NULL,oc[1]->curve,NULL,newc);
	}
	else
	{
	    newc = join_curves(oc[1]->curve,oc[0]->curve,left,right,NULL);
	    roclists_after_join(rp,oc[1]->curve,NULL,oc[0]->curve,NULL,newc);
	}
	if (delete_node(nd[0]) == FUNCTION_FAILED)
	{
	    (void) printf("WARNING in conect_interacting_curves(), "
	                  "unexpected case,  delete_node() failed\n");
	    DEBUG_LEAVE(connect_interacting_curves)
	    return ERROR_IN_STEP;
	}
	if (DEBUG)
	{
	    (void) printf("Newc,\t");
	    show_curve_states(newc);
	}
	interpolate_intfc_states(rp->new_intfc) = sav_interp;
	DEBUG_LEAVE(connect_interacting_curves)
	return GOOD_STEP;
}		/*end connect_interacting_curves*/

/*
*			g_untrack_curve():
*
*	This function is a driver for the function f_untrack_curve.  After
*	deleting the curve using f_untrack_curve, we then check to see if 
*	we are untracking an incident curve.  If we are and if the relevant 
*	recurse flag is set, then we recursively call g_untrack_curve to 
*	untrack all physical vector curves at that node.  This recursion 
*	check is done first at the given node of oc, and then at the oppn.
*/

EXPORT boolean g_untrack_curve(
	O_CURVE		*oc,
	O_CURVE		*oldoc,
	COMPONENT	newcomp,
	double		dt,
	Front		*fr,
	POINTER		wave,
	RPROBLEM	*rp,
	UNTRACK_FLAG	flag)
{
	COMPONENT	comp;
	CURVE		**pc;
	INTERFACE	*intfc = oc->curve->interface;	/*new interface*/
	NODE		*node = Node_of_o_curve(oc);
	NODE		*oppn = Opp_node_of_o_curve(oc);
	NODE		*old_node, *old_oppn, *on;
	O_CURVE		*old_phys_oc;
	O_CURVE		Phys_oc, Old_phys_oc;
	int		status, opp_status;
	ORIENTATION	orient = oc->orient, opp_or = Opposite_orient(orient);
	boolean		sts_set = states_set_at_node(flag,orient);
	boolean		opp_sts_set = states_set_at_node(flag,opp_or);
	boolean         untrack_status;
	boolean         n_ss, on_ss;
	int		num_scal;
	int		i;

	debug_print("untrack","Entered g_untrack_curve()\n");

	if (debugging("untrack"))
	{
	    (void) printf("Untracking curve:  \n");
	    print_curve(oc->curve);
	    (void) printf("newcomp = %d\n",newcomp);
	}

	status = status_at_node(oc->curve,oc->orient);
	opp_status = status_at_node(oc->curve,Opposite_orient(oc->orient));

	untrack_status = f_untrack_curve(oc,oldoc,newcomp,dt,fr,wave,rp,flag);

#define retain_curve_at_node(c)						\
	((wave_type(c) < FIRST_SCALAR_PHYSICS_WAVE_TYPE) ||		\
	(wave_type(c) < FIRST_VECTOR_PHYSICS_WAVE_TYPE && (num_scal > 1)))

	if (((status == INCIDENT) || (status == OVERTOOK)) &&
		    untrack_recurse(flag,orient))
	{
	    num_scal = num_curves_with_wave_type_at_node(node,CONTACT);

	    for (i = 0; i < 2; i++)
	    {
		if (i == 0)
		{
		    pc = node->in_curves;
		    Phys_oc.orient = NEGATIVE_ORIENTATION;
		    flag.end_states_set = sts_set;
		}
		else
		{
		    pc = node->out_curves;
		    Phys_oc.orient = POSITIVE_ORIENTATION;
		    flag.start_states_set = sts_set;
		}
		while (pc && *pc)
		{
		    if (retain_curve_at_node(*pc))
			pc++;
		    else
		    {
		        Phys_oc.curve = *pc;
		        comp = (status_at_node(Phys_oc.curve,Phys_oc.orient) 
				== TRANSMITTED) ?
			    positive_component(Phys_oc.curve) : newcomp;
		        if (oldoc != NULL)
		        {
			    old_node = Node_of_o_curve(oldoc);
			    old_phys_oc = &Old_phys_oc;
			    if (!find_correspond_of_oriented_curve(&Phys_oc,
				    old_phys_oc,old_node,fr,
				    old_node->interface))
			    {
			        screen("ERROR in g_untrack_curve(), "
				       "find_correspond_of_oriented_curve() "
				       "failed\n");
			        clean_up(ERROR);
			    }
		        }
		        else
			    old_phys_oc = NULL;
			on = Opp_node_of_o_curve(&Phys_oc);
			on_ss = (propagation_status(on)==PROPAGATED_NODE) ?
			    YES : NO;
	                set_states_set_at_node_flag(flag,
			    Opposite_orient(Phys_oc.orient),on_ss);
		        g_untrack_curve(&Phys_oc,old_phys_oc,comp,
				        dt,fr,wave,rp,flag);
		        pc = (i == 0) ? node->in_curves : node->out_curves;
		    }
		}
	    }
	}

	if (((opp_status == INCIDENT) || (opp_status == OVERTOOK))
	    		 &&
			untrack_recurse(flag,opp_or))
	{
	    num_scal = num_curves_with_wave_type_at_node(oppn,CONTACT);
	    for (i = 0; i < 2; i++)
	    {
		if (i == 0)
		{
		    pc = oppn->in_curves;
		    Phys_oc.orient = NEGATIVE_ORIENTATION;
		    flag.end_states_set = opp_sts_set;
		}
		else
		{
		    pc = oppn->out_curves;
		    Phys_oc.orient = POSITIVE_ORIENTATION;
		    flag.start_states_set = opp_sts_set;
		}
		while (pc && *pc)
		{
		    if (retain_curve_at_node(*pc))
			pc++;
		    else
		    {
		        Phys_oc.curve = *pc;
		        comp = (status_at_node(Phys_oc.curve,Phys_oc.orient) 
				    == TRANSMITTED) ?
			    positive_component(Phys_oc.curve) : newcomp;
		        if (oldoc != NULL)
		        {
			    old_oppn = Opp_node_of_o_curve(oldoc);
			    old_phys_oc = &Old_phys_oc;
			    if (!find_correspond_of_oriented_curve(&Phys_oc,
				    &Old_phys_oc,old_oppn,fr,
				    old_oppn->interface))
			    {
			        screen("ERROR in g_untrack_curve(), "
				       "find_correspond_of_oriented_curve() "
				       "failed\n");
			        clean_up(ERROR);
			    }
		        }
		        else
			    old_phys_oc = NULL;
			on = Opp_node_of_o_curve(&Phys_oc);
			on_ss = (propagation_status(on)==PROPAGATED_NODE) ?
			    YES : NO;
	                set_states_set_at_node_flag(flag,
			    Opposite_orient(Phys_oc.orient),on_ss);
		        g_untrack_curve(&Phys_oc,old_phys_oc,comp,
				        dt,fr,wave,rp,flag);
		        pc = (i == 0) ? oppn->in_curves : oppn->out_curves;
		    }
		}
	    }
	}

#undef retain_curve_at_node

	if (mono_comp_curves(intfc) == YES && untrack_mono_comp_recurse(flag))
	{
	    UNTRACK_FLAG mono_uflag;

	    pc = intfc->curves;
	    while (pc && *pc)
	    {
	    	if (is_mono_comp_curve(*pc))
		{
	    	    Phys_oc.curve = *pc;
	    	    Phys_oc.orient = POSITIVE_ORIENTATION;
	            n_ss =
		        (propagation_status((*pc)->start)==PROPAGATED_NODE) ?
			    YES : NO;
	            on_ss = (propagation_status((*pc)->end)==PROPAGATED_NODE) ?
			    YES : NO;
	    	    set_untrack_flag(mono_uflag,Phys_oc.orient,
	    	    		     n_ss,on_ss,YES,YES,YES);
	    	    g_untrack_curve(&Phys_oc,NULL,
	    		            positive_component(Phys_oc.curve),dt,fr,
				    wave,rp,mono_uflag);
	    	    pc = intfc->curves;
		}
		else
		    pc++;
	    }
	}

	debug_print("untrack","Left g_untrack_curve()\n");
	return untrack_status;
}		/*end g_untrack_curve*/


LOCAL	int num_curves_with_wave_type_at_node(
	NODE		*node,
	int		w_type)
{
	CURVE		**pc;
	int		i, num_cur = 0;

	DEBUG_ENTER(num_curves_with_wave_type_at_node)
	for (i = 0; i < 2; i++)
	{
	    for (pc = (i == 0) ? node->in_curves : node->out_curves;
							pc && *pc; pc++)
	    if (wave_type(*pc) == w_type) num_cur++;
	}
	DEBUG_LEAVE(num_curves_with_wave_type_at_node)
	return num_cur;
}		/*end num_curves_with_wave_type_at_node*/


LOCAL	void	g_init_rp_nodes(
	RPROBLEM	*rp)
{
	RP_NODE		*rpn;

	DEBUG_ENTER(g_init_rp_nodes)
	f_init_rp_nodes(rp);
	for(rpn = rp->first_rp_node; rpn; rpn = rpn->next)
	{
	    rpn_inc_shock1(rpn)  = rpn_inc_shock2(rpn)  = NULL;
	    rpn_reflected1(rpn) = rpn_reflected2(rpn) = NULL;
	    rpn_contact1(rpn) = rpn_contact2(rpn) = NULL;
	    rpn_transmitted(rpn) = NULL;
	}
	DEBUG_LEAVE(g_init_rp_nodes)
}		/*end g_init_rp_nodes*/

LOCAL	void g_delete_curve_from_rp_node(
	CURVE		*curve,
	RP_NODE		*rpn,
	RPROBLEM	*rp)
{

	DEBUG_ENTER(g_delete_curve_from_rp_node)
	f_delete_curve_from_rp_node(curve,rpn,rp);
	if (!is_bdry(curve))
	{
	    delete_curve_from_o_curve_family(curve,&rpn_inc_shock1(rpn));
	    delete_curve_from_o_curve_family(curve,&rpn_inc_shock2(rpn));
	    delete_curve_from_o_curve_family(curve,&rpn_reflected1(rpn));
	    delete_curve_from_o_curve_family(curve,&rpn_reflected2(rpn));
	    delete_curve_from_o_curve_family(curve,&rpn_transmitted(rpn));
	    delete_curve_from_o_curve_family(curve,&rpn_contact1(rpn));
	    delete_curve_from_o_curve_family(curve,&rpn_contact2(rpn));
	}
	DEBUG_LEAVE(g_delete_curve_from_rp_node)
}		/*end g_delete_curve_from_rp_node*/

/*ARGSUSED*/
LOCAL	void	g_user_free_rp_node(
	RP_NODE		*rpn,
	RPROBLEM	*rp)
{

	DEBUG_ENTER(g_user_free_rp_node)
	free_o_curve_family(rpn_inc_shock1(rpn));
	free_o_curve_family(rpn_inc_shock2(rpn));
	free_o_curve_family(rpn_reflected1(rpn));
	free_o_curve_family(rpn_reflected2(rpn));
	free_o_curve_family(rpn_transmitted(rpn));
	free_o_curve_family(rpn_contact1(rpn));
	free_o_curve_family(rpn_contact2(rpn));
	DEBUG_LEAVE(g_user_free_rp_node)
}		/*end g_user_free_rp_node*/


/*ARGSUSED*/
LOCAL	void	g_user_print_rp_node(
	RP_NODE		*rpn,
	RPROBLEM	*rp)
{

	DEBUG_ENTER(g_user_print_rp_node)
	(void) printf("Incident curve families:\n");
	print_o_curve_family(rpn_inc_shock1(rpn));
	print_o_curve_family(rpn_inc_shock2(rpn));
	(void) printf("Reflected curve families:\n");
	print_o_curve_family(rpn_reflected1(rpn));
	print_o_curve_family(rpn_reflected2(rpn));
	(void) printf("Transmitted curve family:\n");
	print_o_curve_family(rpn_transmitted(rpn));
	(void) printf("Contact curve families:\n");
	print_o_curve_family(rpn_contact1(rpn));
	print_o_curve_family(rpn_contact2(rpn));
	DEBUG_LEAVE(g_user_print_rp_node)
}		/*end g_user_print_rp_node*/

LOCAL	void	g_user_print_rproblem(
	RPROBLEM	*rp)
{

	DEBUG_ENTER(g_user_print_rproblem)
	(void) printf("Number of incident curves = %d\n",rp_num_incident(rp));
	(void) printf("Number of reflected curves = %d\n",rp_num_refl(rp));
	(void) printf("Number of transmitted curves = %d\n",rp_num_transm(rp));
	(void) printf("Number of incident shocks = %d\n",rp_num_inc_shock(rp));
	(void) printf("Number of contacts = %d\n",rp_num_contact(rp));
	DEBUG_LEAVE(g_user_print_rproblem)
}		/*end g_user_print_rproblem*/

/*ARGSUSED*/
LOCAL	void	g_set_phys_ocurves_to_null(
	RP_NODE		*rpn,
	RPROBLEM	*rp)
{

	DEBUG_ENTER(g_set_phys_ocurves_to_null)
	rpn_inc_shock1(rpn) = NULL;
	rpn_reflected1(rpn) = NULL;
	rpn_transmitted(rpn) = NULL;
	rpn_inc_shock2(rpn) = NULL;
	rpn_reflected2(rpn) = NULL;
	rpn_contact1(rpn) = NULL;
	rpn_contact2(rpn) = NULL;
	DEBUG_LEAVE(g_set_phys_ocurves_to_null)
}		/*end g_set_phys_ocurves_to_null*/
#endif /* defined(TWOD) */
