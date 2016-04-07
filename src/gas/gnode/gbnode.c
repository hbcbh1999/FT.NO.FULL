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
*				gbnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the functions for the first advance of nodes corresponding
*	to physical-boundary curve interactions.
*	
*
*/

#if defined(TWOD)

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	void	modify_regular_reflection_node(POINT*,BOND*,BOND*,O_CURVE*,
					       O_CURVE*,O_CURVE*,O_CURVE*,
					       O_CURVE*,O_CURVE*,
					       ANGLE_DIRECTION,ANGLE_DIRECTION,
					       Front*,Wave*,NODE_FLAG,double,
					       RP_DATA*);

/*
*			B_reflect_node_propagate():
*
*
*            reflected (ang 2)    RP->state[1]    incident (ang 1)
*                             \                  /
*                              \                /
*                               \              /
*                                \            /
*                                 \          /
*          RP->state[2]            \        /             RP->state[0]
*                                   \      /
*                                    \    /
*                                     \  /
*  behind wall (ang 3)_________________\/___________________ahead wall (ang 0)
*
*
*	A B_reflect_node is a node at which four curves meet, with two
*	having a shock wave type and the others having NEUMANN boundary
*	wave types.
*
*	The forward facing (low pressure) side of the incident shock defines 
*	an ahead component, and hence an ahead and behind boundary and an
*	incident to ahead angular orientation of the node. Propagation
*	of the incident shock also defines a component into which
*	the node propagates, and hence a (possibly distinct) angular 
*	orientation.  If the node propagates into the the physical domain,
*	the physical curve is extended to meet the boundary using 
*	H_extend_crossing_of_two_propagated_curves().  This case is not the
*	normal case, and is presumably accompanied by a bifurcation to 
*	a Mach triple point configuration.  
*
*	Otherwise the intersection of the propagated physical curve
*	is found with the propagated ahead boundary.  The the ahead
*	and behind boundaries are updated by calling shift_node(),
*	and the physical curve is updated using cut_curve().
*
*/

EXPORT	int B_reflect_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	BOND		*crossbinc;		/* intersecting bonds */
	BOND		*crossbahead;		/*  on newcinc, newcahead */
	BOND		*bcorner;
	COMPONENT	ahead_comp;
	COMPONENT	propagation_comp;
	NODE		*interact_nodes[3];
	O_CURVE		Oldcinc, Newcinc;	/* incident curve */
	O_CURVE		Oldcref, Newcref;	/* reflected curve */
	O_CURVE		Oldcahead, Newcahead;	/* ahead boundary */
	O_CURVE		Oldcbehind, Newcbehind;	/* behind boundary */
	O_CURVE		Oldcaprop, Newcaprop;	/* prop dir ahead bdry */
	O_CURVE		Oldcbprop, Newcbprop;	/* prop dir behind bdry */
	POINT		Pc;			/* cross point */
	POINT		*corner;
	RP_DATA		*RP;
	double		tcr_inc;	/* fractional distances to cross */
	double		tcr_ahead;
	double		ta[MAXD];		/* ahead wall tangent */
	double		tb[MAXD];		/* behind wall tangent */
	double		t[MAXD];		/* inc tangent */
	double		node_v[MAXD];
	int		status;
	ANGLE_DIRECTION	i_to_a_dir;	/* dir of angle cinc to  bdry faced by
					 *  low pressure side of cinc  */
	ANGLE_DIRECTION	i_to_prop_dir;	/* dir of angle cinc to bdry in
					    direction of node propagation */
	ORIENTATION	tmp_orient;
	SIDE		inc_side;
	SIDE		propagation_side;
	int		forward_facing_is_ahead;

	debug_print("B_reflect_node","Entered B_reflect_node_propagate()\n");
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("B_reflect_node")) 
	{
	    (void) printf("\n\tOLD NODE:\n");
	    print_node(oldn);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

		/* Identify curves and components */

	find_curve_with_status(oldn,&Oldcinc.curve,&Oldcinc.orient,INCIDENT);
	Check_return(find_correspond_of_oriented_curve(&Oldcinc,&Newcinc,
					               newn,fr,newn->interface),
	             B_reflect_node_propagate)

		/* In the case of untracked refl, it is assumed that the
		 * following two functions set the curve fields to NULL.
		 */

	find_curve_with_status(oldn,&Oldcref.curve,&Oldcref.orient,REFLECTED);
	(void) find_correspond_of_oriented_curve(&Oldcref,&Newcref,
						 newn,fr,newn->interface);

	find_propagation_orientation(fr,(POINTER)wave,oldn,newn,NULL,
	                             &Oldcinc,dt,&i_to_prop_dir,
				     &Oldcahead,&Newcahead,
				     &Oldcbehind,&Newcbehind,&inc_side,
				     &propagation_side,&ahead_comp,
				     &propagation_comp);
	find_adjacent_curves(&Oldcinc,&i_to_a_dir,&Oldcahead,&Oldcbehind,
		             &ahead_comp);
	if (Oldcref.curve != NULL)
	{
	    Oldcbehind.curve = adjacent_curve(Oldcref.curve,Oldcref.orient,
			                      Opposite_ang_dir(i_to_a_dir),
					      &Oldcbehind.orient);
	}
	else
	{
	    Oldcbehind.curve = adjacent_curve(Oldcinc.curve,Oldcinc.orient,
			                      Opposite_ang_dir(i_to_a_dir),
					      &Oldcbehind.orient);
	}
	Check_return(find_correspond_of_oriented_curve(&Oldcbehind,&Newcbehind,
						       newn,fr,newn->interface),
	             B_reflect_node_propagate)

	forward_facing_is_ahead = (i_to_prop_dir == i_to_a_dir); 
	if (forward_facing_is_ahead)
	{
	    copy_o_curve(&Oldcaprop,&Oldcahead);
	    copy_o_curve(&Newcaprop,&Newcahead);
	    copy_o_curve(&Oldcbprop,&Oldcbehind);
	    copy_o_curve(&Newcbprop,&Newcbehind);
	}
	else 
	{
	    copy_o_curve(&Oldcaprop,&Oldcbehind);
	    copy_o_curve(&Newcaprop,&Newcbehind);
	    copy_o_curve(&Oldcbprop,&Oldcahead);
	    copy_o_curve(&Newcbprop,&Newcahead);
	}
	if ((continue_past_fixed_node(flag) == YES) && 
	    is_fixed_node(Opp_node_of_o_curve(&Oldcaprop)))
	{
	    Check_return(next_boundary(Oldcaprop.curve,
				       Opposite_orient(Oldcaprop.orient),
			               &Oldcaprop.curve,&tmp_orient),
		         B_reflect_node_propagate)
	    if (Oldcaprop.orient != tmp_orient) 
	    {
	    	inc_side = Opposite_side(inc_side);
	    	propagation_side = Opposite_side(propagation_side);
	    }
	    Oldcaprop.orient =  tmp_orient;
	    Check_return(next_boundary(Newcaprop.curve,
			               Opposite_orient(Newcaprop.orient),
			               &Newcaprop.curve,&Newcaprop.orient),
		         B_reflect_node_propagate)
	}

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("B_reflect_node")) 
	{
	    (void) printf("\t\tOLD INCIDENT CURVE:\n");
	    print_o_curve(&Oldcinc);
	    (void) printf("\t\tOLD REFLECTED CURVE:\n");
	    print_o_curve(&Oldcref);
	    (void) printf("\t\tOLD AHEAD BOUNDARY:\n");
	    print_o_curve(&Oldcahead);
	    (void) printf("\t\tOLD BEHIND BOUNDARY:\n");
	    print_o_curve(&Oldcbehind);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */


		/* Identify new position of node */

	if (propagation_side == inc_side) 
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("B_reflect_node","oldn propagates into ahead comp\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = H_extend_crossing_of_two_propagated_curves(&Oldcaprop,
				&Newcaprop,&Oldcinc,&Newcinc,ahead_comp,
				propagation_comp,&Pc,&crossbahead,&crossbinc,
				&tcr_ahead,&tcr_inc,fr,(POINTER)wave,rp,
				dt,dt_frac,flag);
	}
	else 
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("B_reflect_node","oldn propagates out of ahead comp\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = crossing_of_two_propagated_curves(&Oldcinc,&Newcinc,
				&Oldcaprop,&Newcaprop,&Pc,&crossbinc,
				&crossbahead,&tcr_inc,&tcr_ahead,
				fr,(POINTER)wave,rp,dt,dt_frac,flag);

	}
	if (status != GOOD_NODE)
	{
	    debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
	    return status;
	}

		/* Keep cross point on outer boundary */

	if (is_bdry(newn) && (to_next_node_only(flag) == YES))
	    nearest_boundary_point(Coords(&Pc),Coords(&Pc),fr->rect_grid);


		/* Read off states from the propagated incident curve */

	RP = Rp_data(newn);
	RP->ang_dir = Opposite_ang_dir(i_to_a_dir);
	find_tangent_to_propagated_curve(&Pc,crossbahead,
			&Oldcaprop,&Newcaprop,ta,fr,(POINTER)wave,dt);
	find_tangent_to_propagated_curve(&Pc,
			Bond_at_node_of_o_curve(&Newcbprop),
			&Oldcbprop,&Newcbprop,tb,fr,(POINTER)wave,dt);
	if (forward_facing_is_ahead)
	{
	    RP->ang[0] = angle(ta[0],ta[1]);
	    RP->ang[3] = angle(tb[0],tb[1]);
	}
	else
	{
	    RP->ang[3] = angle(ta[0],ta[1]);
	    RP->ang[0] = angle(tb[0],tb[1]);
	}
	find_tangent_to_propagated_curve(&Pc,crossbinc,&Oldcinc,&Newcinc,t,fr,
					 (POINTER)wave,dt);
	RP->ang[1] = angle(t[0],t[1]);

	if (curve_ang_oriented_l_to_r(i_to_a_dir,Newcinc.orient))
	{
	    left_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,RP->state[1]);
	    right_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,
			           RP->state[0]);
	}
	else 
	{
	    left_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,RP->state[0]);
	    right_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,
			           RP->state[1]);
	}

		/* Calculate the new state according to shock polars */

	if (continue_past_fixed_node(flag) == YES) 
	{
	    bcorner = crossbahead;
	    set_corner_for_bifurcation(&Oldcinc,&Newcinc,&Oldcaprop,&Newcaprop,
				       &Oldcbprop,&Newcbprop,oldn->posn,&Pc,fr,
				       wave,&corner,&bcorner,
				       forward_facing_is_ahead,
				       node_v,RP,dt,dt_frac);
	    Node_vel(newn)[0] = node_v[0];
	    Node_vel(newn)[1] = node_v[1];
	    status = velocity_satisfies_CFL(newn,dt,dt_frac,fr);
	    if (status != GOOD_NODE)
	    {
	    	debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
		return status;
	    }
	}
	else 
	{
	    node_v[0] = Node_vel(newn)[0];
	    node_v[1] = Node_vel(newn)[1];
	}

	if (is_regular_reflection(node_v,fr,RP)) 
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("B_reflect_node","\t\tREGULAR REFLECTION:\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    modify_regular_reflection_node(&Pc,crossbinc,crossbahead,
			                   &Oldcinc,&Newcinc,&Oldcref,&Newcref,
					   &Newcaprop,&Newcbprop,i_to_a_dir,
					   i_to_prop_dir,fr,wave,flag,dt,RP);
		status = GOOD_NODE;
	}
	else if (continue_past_bifurcation(flag) != YES) 
	{
	    interact_nodes[0] = newn;
	    interact_nodes[1] = oldn;
	    interact_nodes[2] = NULL;
	    if (continue_past_fixed_node(flag) != YES)
	    	set_corner_for_bifurcation(&Oldcinc,&Newcinc,&Oldcaprop,
					   &Newcaprop,&Oldcbprop,&Newcbprop,
				           oldn->posn,&Pc,fr,wave,&corner,
					   &crossbahead,forward_facing_is_ahead,
					   node_v,RP,dt,dt_frac);
	    augment_rproblem_list(rp,interact_nodes,dt,*dt_frac,oldn->interface,
				  newn->interface,fr,(POINTER)wave);
	    propagation_status(newn) = VEL_COMPUTED_NODE;
	    status = BIFURCATION_NODE;
	    debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
	    return status;
	}
	else 
	{
	    bcorner = crossbahead;
	    if (continue_past_fixed_node(flag) != YES)
	    	set_corner_for_bifurcation(&Oldcinc,&Newcinc,&Oldcaprop,
					   &Newcaprop,&Oldcbprop,&Newcbprop,
					   oldn->posn,&Pc,fr,wave,&corner,
					   &bcorner,forward_facing_is_ahead,
					   node_v,RP,dt,dt_frac);
	    g_reflect_node_bifurcation(fr,wave,&Oldcinc,&Newcinc,&Oldcref,
				       &Newcref,&Oldcahead,&Newcahead,
				       &Oldcbehind,&Newcbehind,&Oldcaprop,
				       &Newcaprop,&Oldcbprop,&Newcbprop,
				       oldn->posn,&Pc,corner,crossbinc,
				       crossbahead,RP,*rp,dt,dt_frac,flag);
	    propagation_status(newn) = PROPAGATED_NODE;
	    debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
	    return GOOD_NODE;
	}
	propagation_status(newn) = PROPAGATED_NODE;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("B_reflect_node")) 
	{
	    (void) printf("\n\tNEW NODE:\n");
	    print_node(newn);
	    (void) printf("\t\tNEW INCIDENT CURVE:\n");
	    print_o_curve(&Newcinc);
	    (void) printf("\t\tNEW REFLECTED CURVE:\n");
	    print_o_curve(&Newcref);
	    (void) printf("\t\tNEW AHEAD BOUNDARY:\n");
	    print_o_curve(&Newcahead);
	    (void) printf("\t\tNEW BEHIND BOUNDARY:\n");
	    print_o_curve(&Newcbehind);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	debug_print("B_reflect_node","Left B_reflect_node_propagate()\n");
	return status;
}		/*end B_reflect_node_propagate*/


/*
*               modify_regular_reflection_node():
*
*	Given that a regular reflection node remains such a node under
*	propagation, this routine sets up the local environment of the
*	propagated node.  The node newn on the incident curve newcinc
*	is shifted along the curves newcehind and newcahead to the position
*	given by the cross-point pc on bond newbahead.  The node newn is
*	also propagated with respect to the reflected curve newcref, becoming
*	the last "interacting" point of that curve with respect to the node;
*	a new point is inserted to maintain the correct reflection angle.
*	The ahead component is ahead_comp, and the state date is specified
*	by the RP_DATA structure *RP.
*/

LOCAL void modify_regular_reflection_node(
	POINT		*pc,
	BOND		*newbinc,
	BOND		*newbahead,
	O_CURVE		*oldcinc,
	O_CURVE		*newcinc,
	O_CURVE		*oldcref,
	O_CURVE		*newcref,
	O_CURVE		*newcahead,
	O_CURVE		*newcbehind,
	ANGLE_DIRECTION	i_to_a_dir,
	ANGLE_DIRECTION	i_to_prop_dir,
	Front		*fr,
	Wave		*wave,
	NODE_FLAG	flag,
	double		dt,
	RP_DATA		*RP)
{
	NODE		*oldn,*newn;
	Locstate	st_left_inc, st_right_inc;
	Locstate	st_left_ref, st_right_ref;
	Locstate	st_left_ahead, st_right_ahead;
	Locstate	st_left_behind, st_right_behind;

	debug_print("B_reflect_node","Entered modify_regular_reflection_node()\n");

		/* Identify left and right states on curves */

#define left_is_toward_incident(i_to_a_dir,c_orient)			\
	((i_to_a_dir) == CLOCKWISE ?					\
		((c_orient) == POSITIVE_ORIENTATION ? YES : NO) :	\
		((c_orient) == NEGATIVE_ORIENTATION ? YES : NO))

	if (left_is_toward_incident(i_to_a_dir,newcahead->orient)) 
	{
		if (i_to_a_dir == i_to_prop_dir)
			st_left_ahead = RP->state[0];
		else	st_left_ahead = RP->state[2];
		st_right_ahead = return_obst_state();
	}
	else 
	{
		st_left_ahead = return_obst_state();
		if (i_to_a_dir == i_to_prop_dir)
			st_right_ahead = RP->state[0];
		else st_right_ahead = RP->state[2];
	}
	if (left_is_toward_incident(i_to_a_dir,newcbehind->orient)) 
	{
		st_left_behind = return_obst_state();
		if (i_to_a_dir == i_to_prop_dir)
			st_right_behind = RP->state[2];
		else st_right_behind = RP->state[0];
	}
	else 
	{
		if (i_to_a_dir == i_to_prop_dir)
			st_left_behind = RP->state[2];
		else st_left_behind = RP->state[0];
		st_right_behind = return_obst_state();
	}
	if (left_is_toward_incident(i_to_a_dir,newcinc->orient)) 
	{
		st_left_inc  = RP->state[1];
		st_right_inc = RP->state[0];
	}
	else 
	{
		st_left_inc  = RP->state[0];
		st_right_inc = RP->state[1];
	}

		/* Shift node on wall and assign wall states */

	newn = Node_of_o_curve(newcinc);
	debug_print("B_reflect_node","Calling shift_node_past\n");
	shift_node_past(pc,newbahead,newcahead->curve,newcahead->orient,
		newcbehind->curve,newcbehind->orient,i_to_a_dir,newn,fr,flag,
		st_left_ahead,st_right_ahead,st_left_behind,st_right_behind);

	oldn = Node_of_o_curve(oldcinc);

		/* Cut extra part of the incident curve */

	(void) propagate_curve_near_node(oldn,newn,oldcinc,newcinc,
					 newbinc,st_left_inc,st_right_inc,
					 NO,RP->ang[1],fr,wave,dt,flag);

		/* Propagate old reflected curve and insert a new bond */

	if ((newcref != NULL) && (newcref->curve != NULL))
	{
		if (left_is_toward_incident(i_to_a_dir,newcref->orient)) 
		{
			st_left_ref = RP->state[2];
			st_right_ref = RP->state[1];
		}
		else 
		{
			st_left_ref = RP->state[1];
			st_right_ref = RP->state[2];
		}
		(void) propagate_curve_near_node(oldn,newn,oldcref,newcref,
			Bond_at_node_of_o_curve(newcref),st_left_ref,
			st_right_ref,YES,RP->ang[2],fr,wave,dt,flag);
	}
	debug_print("B_reflect_node","Left modify_regular_reflection_node()\n");
}		/*end modify_regular_reflection_node*/


/*
*			attached_b_node_propagate():
*
*	This function handles the case of a curve that is fixed at one
*	end at a corner where multiple boundary-like curves meet.  The
*	propagation is relatively simple.  A call to is made to
*	is_bow_shock_attached() to determine if the curve can remain
*	attached.  If so, the new states returned are ft_assigned.  If not,
*	a bifurcation is signalled.
*/

/* ARGSUSED */
EXPORT int attached_b_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	Locstate	ahead, behind;
	Locstate	new_ahead, new_behind;
	NODE		*interact_node[3];
	O_CURVE		Ci, NewCi;
	O_CURVE		Ca, NewCa;
	O_CURVE		Cb, NewCb;
	POINT		*oldp;
	double		b_ang;		/* pos x axis to wall behind */
	double		b_t[MAXD];	/* direction vector for b_ang */
	double		shock_angle;	/* pos x axis to attached shock */
	double		dir[MAXD];	/* direction vector for shock_angle */
	double		V[MAXD];
	double		len;
	size_t		sizest = fr->sizest;
	int		status;
	static POINT	*newp = NULL;
	
	debug_print("attached_b_node","Entered attached_b_node_propagate()\n");

	if (newp == NULL)
	    newp = Static_point(fr->interf);

	Ci.curve = find_physical_curve_at_node(oldn,&Ci.orient);
	NewCi.curve = find_physical_curve_at_node(newn,&NewCi.orient);
	oldp = Node_of_o_curve(&Ci)->posn;
	if (wave_type(Ci.curve) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE) 
	{
	    point_propagate(fr,(POINTER)wave,oldp,newp,
			    Bond_at_node_of_o_curve(&Ci),Ci.curve,dt,V);
	    if (wave_type(Ci.curve) == FORWARD_SHOCK_WAVE) 
	    {
	        ahead = right_state(newp);
	    	behind = left_state(newp);
	    	if (Ci.orient == POSITIVE_ORIENTATION)
	    	{
	    	    Ca.curve = adjacent_curve(Ci.curve,Ci.orient,
					      CLOCKWISE,&Ca.orient);
	    	    Cb.curve = adjacent_curve(Ci.curve,Ci.orient,
					      COUNTER_CLOCK,&Cb.orient);
	    	}
	    	else
	    	{
	    	    Ca.curve = adjacent_curve(Ci.curve,Ci.orient,
					      COUNTER_CLOCK,&Ca.orient);
	    	    Cb.curve = adjacent_curve(Ci.curve,Ci.orient,
					      CLOCKWISE,&Cb.orient);
	    	}
	    }
	    else 
	    {
	    	ahead = left_state(newp);
	    	behind = right_state(newp);
	    	if (Ci.orient == POSITIVE_ORIENTATION)
	    	{
	    	    Ca.curve = adjacent_curve(Ci.curve,Ci.orient,
					      COUNTER_CLOCK,&Ca.orient);
	       	    Cb.curve = adjacent_curve(Ci.curve,Ci.orient,
					      CLOCKWISE,&Cb.orient);
	    	}
	    	else
	    	{
	    	    Ca.curve = adjacent_curve(Ci.curve,Ci.orient,
					      CLOCKWISE,&Ca.orient);
	    	    Cb.curve = adjacent_curve(Ci.curve,Ci.orient,
					      COUNTER_CLOCK,&Cb.orient);
	    	}
	    }
	    if (wave_type(NewCi.curve) == FORWARD_SHOCK_WAVE) 
	    {
	        new_ahead = Right_state_at_node_of_o_curve(&NewCi);
	        new_behind = Left_state_at_node_of_o_curve(&NewCi);
	        if (NewCi.orient == POSITIVE_ORIENTATION)
	        {
	    	    NewCa.curve = adjacent_curve(NewCi.curve,NewCi.orient,
						 CLOCKWISE,&NewCa.orient);
	    	    NewCb.curve = adjacent_curve(NewCi.curve,NewCi.orient,
						 COUNTER_CLOCK,&NewCb.orient);
	        }
	        else
	        {
	    	    NewCa.curve = adjacent_curve(NewCi.curve,NewCi.orient,
						 COUNTER_CLOCK,&NewCa.orient);
	    	    NewCb.curve = adjacent_curve(NewCi.curve,NewCi.orient,
						 CLOCKWISE,&NewCb.orient);
	        }
	    }
	    else
	    {
	        new_ahead = Left_state_at_node_of_o_curve(&NewCi);
	        new_behind = Right_state_at_node_of_o_curve(&NewCi);
	        if (NewCi.orient == POSITIVE_ORIENTATION)
	        {
	    	    NewCa.curve = adjacent_curve(NewCi.curve,NewCi.orient,
						 COUNTER_CLOCK,&NewCa.orient);
	    	    NewCb.curve = adjacent_curve(Ci.curve,NewCi.orient,
						 CLOCKWISE,&NewCb.orient);
	        }
	        else
	        {
	    	    NewCa.curve = adjacent_curve(NewCi.curve,NewCi.orient,
						 CLOCKWISE,&NewCa.orient);
	    	    NewCb.curve = adjacent_curve(NewCi.curve,NewCi.orient,
						 COUNTER_CLOCK,&NewCb.orient);
	        }
	    }
	    find_tangent_to_curve(Node_of_o_curve(&Cb)->posn,
			          Bond_at_node_of_o_curve(&Cb),Cb.curve,
			          Cb.orient,b_t,fr);
	    b_ang = angle(b_t[0],b_t[1]);
	    if (is_bow_shock_attached(ahead,b_ang,behind,&shock_angle)) 
	    {
	        ft_assign(new_ahead,ahead,sizest);
	        ft_assign(new_behind,behind,sizest);

	        len = bond_length(Bond_at_node_of_o_curve(&Ci));
	        dir[0] = cos(shock_angle);	dir[1] = sin(shock_angle);

	        if (wave_type(NewCi.curve) == FORWARD_SHOCK_WAVE)
	        {
	    	    (void) adjust_angle_at_node(newn,&Ci,&NewCi,new_ahead,
						new_behind,dir,dt,len,fr,wave);

	    	    if (NewCi.orient == POSITIVE_ORIENTATION)
	    	    {
	    	        if (NewCb.orient == POSITIVE_ORIENTATION)
	    		    ft_assign(right_start_state(NewCb.curve),new_behind,
			           sizest);
	    	        else
	    		    ft_assign(left_end_state(NewCb.curve),new_behind,
				   sizest);

	    	        if (NewCa.orient == POSITIVE_ORIENTATION)
	    		    ft_assign(left_start_state(NewCa.curve),new_ahead,
				   sizest);
	    	        else
	    		    ft_assign(right_end_state(NewCa.curve),new_ahead,
				   sizest);
	    	    }
	    	    else
	    	    {
	    	        if (NewCb.orient == POSITIVE_ORIENTATION)
	    		    ft_assign(left_start_state(NewCb.curve),new_behind,
				   sizest);
	    	        else
	    		    ft_assign(right_end_state(NewCb.curve),new_behind,
				   sizest);
    
	    	        if (NewCa.orient == POSITIVE_ORIENTATION)
	    		    ft_assign(right_start_state(NewCa.curve),new_ahead,
				   sizest);
	    	        else
	    		    ft_assign(left_end_state(NewCa.curve),new_ahead,
				   sizest);
	    	    }
	        }
	        else
		{
		    (void) adjust_angle_at_node(newn,&Ci,&NewCi,new_behind,
						new_ahead,dir,dt,len,fr,wave);

		    if (NewCi.orient == POSITIVE_ORIENTATION)
		    {
		        if (NewCb.orient == POSITIVE_ORIENTATION)
		    	    ft_assign(left_start_state(NewCb.curve),new_behind,
				   sizest);
		        else
		    	    ft_assign(right_end_state(NewCb.curve),new_behind,
				   sizest);

		        if (NewCa.orient == POSITIVE_ORIENTATION)
		    	    ft_assign(right_start_state(NewCa.curve),new_ahead,
				   sizest);
		        else
		    	    ft_assign(left_end_state(NewCa.curve),new_ahead,
				   sizest);
		    }
		    else
		    {
		        if (NewCb.orient == POSITIVE_ORIENTATION)
		    	    ft_assign(right_start_state(NewCb.curve),new_behind,
				   sizest);
		        else
		    	    ft_assign(left_end_state(NewCb.curve),new_behind,
				   sizest);

		        if (NewCa.orient == POSITIVE_ORIENTATION)
		    	    ft_assign(left_start_state(NewCa.curve),new_ahead,
				   sizest);
		        else
		    	    ft_assign(right_end_state(NewCa.curve),new_ahead,
				   sizest);
		    }
		}
	        debug_print("attached_b_node","Left attached_b_node_propagate()\n");
		return GOOD_NODE;
	    }
	    if (continue_past_bifurcation(flag) != YES) 
	    {
	    	interact_node[0] = newn;
	    	interact_node[1] = oldn;
	    	interact_node[2] = NULL;
	    	augment_rproblem_list(rp,interact_node,dt,*dt_frac,
				      Ci.curve->interface,
				      NewCi.curve->interface,fr,(POINTER)wave);
	        debug_print("attached_b_node","Left attached_b_node_propagate()\n");
	    	return BIFURCATION_NODE;
	    }
	    else 
	    {
	    	/* TODO: Carry out bifurcation */
	        debug_print("attached_b_node","Left attached_b_node_propagate()\n");
	    	return BIFURCATION_NODE;
	    }
	}
	else 
	{
			/* TODO: Improve code for attached contact */
	    status = fixed_node_propagate(fr,(POINTER)wave,oldn,newn,dt);
	    debug_print("attached_b_node","Left attached_b_node_propagate()\n");
	    return status;
	}
}		/*end attached_b_node_propagate*/


			/* SHOCK POLAR ANALYSIS */

/*
*			is_regular_reflection():
*
*       This is shock boundary regular reflection.   The input and output
*	states st_0 and st_1 in RP are given in the computational frame.
*	st_0 is the state ahead and st_1 is the state behind the
*       incident shock, while st_2 is the state behind reflected shock.
*       Returns 1 if there is a regular reflection according to the von Neumann
*	criteria, and 0 otherwise. If there is a regular reflection, then the
*	state st_2 behind the reflected shock is computed and placed in
*	RP->state[2]. This state is given in the computational frame.
*	It is assumed that the wall is stationary in the computational frame.
*	TODO: Remove this assumption.
*/

EXPORT int is_regular_reflection(
	double		*nod_v,
	Front		*fr,
	RP_DATA		*RP)
{
	int		is_given_ahead = YES, is_back_weak = YES;
	double		turn_angle;
	double		M0sq, p2;
	double		rv0[MAXD];
	int		debug_sp2 = NO;
	
	debug_print("B_reflect_node","Entered is_regular_reflection()\n");

	if (debugging("B_reflect_node"))
	{
		(void) printf("nod_v = %g, %g\n",nod_v[0],nod_v[1]);
		(void) printf("rest frame state0 vel = %g %g\n",
			vel(0,RP->state[0]) - nod_v[0],
			vel(1,RP->state[0]) - nod_v[1]);
		(void) printf("rest frame state1 vel = %g %g\n",
			vel(0,RP->state[1]) - nod_v[0],
			vel(1,RP->state[1]) - nod_v[1]);
	}

	/*	Find -1.0 * turn_angle from rest frame velocities */

	M0sq = mach_number_squared(RP->state[0],nod_v,rv0);
	turn_angle = angle(rv0[0],rv0[1]) -
			angle(vel(0,RP->state[1]) - nod_v[0],
				vel(1,RP->state[1]) - nod_v[1]);

	if (debugging("B_reflect_node"))
	{
		(void) printf("state0:");
		(*fr->print_state)(RP->state[0]);
		(void) printf("state1:");
		(*fr->print_state)(RP->state[1]);
		print_angle("turn_angle =",-turn_angle,"\n");
		if (debugging("s_polar_2"))	debug_sp2 = YES;
		else				add_to_debug("spolar2");
	}

	if (!s_polar_2(RP->state[1],is_given_ahead,is_back_weak,
			turn_angle,nod_v,RP->state[2],&RP->ang[2]))
	{
		/* detachment transition criterion */
		debug_print("B_reflect_node","Left is_regular_reflection()\n");
		return NO;
	}

	p2 = pressure(RP->state[2]);
	if (p2 > max_behind_shock_pr(M0sq,RP->state[0]))
	{
		/* Mechanical equilibrium transition criteria */
		debug_print("B_reflect_node","Left is_regular_reflection()\n");
		return NO;
	}

	if (debugging("B_reflect_node"))
	{
		if (!debug_sp2)	remove_from_debug("spolar2");
		print_RP_node_states("States after is_regular_reflection()",
			nod_v,RP,B_REFLECT_NODE);
	}
	debug_print("B_reflect_node","Left is_regular_reflection()\n");
	return YES;
}		/*end is_regular_reflection*/

/*
*			is_bow_shock_attached():
*
*	Given the ahead state and the absolute wedge angle,
*	this routine checks whether it is possible to form an
*	attached bow shock.  
*
*	Returns YES and the behind_state, attached shock angle if there is an
*	attached bow shock; otherwise return NO.
*/

EXPORT int is_bow_shock_attached(
	Locstate	ahead_state,
	double		wedge_angle, /* pos x axis to wall behind attach bow */
	Locstate	answer_state,
	double		*shock_angle)/* pos x axis to shock attached at bow */
{
	double		turn_angle;
	int		is_given_ahead = YES;
	int		is_behind_weak = YES;
	static double	node_v[MAXD];	/* statics initialized to zero*/

	debug_print("is_bow_shock_attached","Entered is_bow_shock_attached()\n");

	turn_angle = normalized_angle(wedge_angle - 
			angle(Mom(ahead_state)[0],Mom(ahead_state)[1]));
	if (debugging("is_bow_shock_attached"))
	{
		print_angle("wedge_angle =",wedge_angle,"\n");
		print_angle("turn_angle =",turn_angle,"\n");
	}
	if (turn_angle >= PI) turn_angle = 2.0*PI - turn_angle;
	if (!s_polar_2(ahead_state,is_given_ahead,is_behind_weak,
			turn_angle,node_v,answer_state,shock_angle))
	{
		debug_print("is_bow_shock_attached",
			"Left is_bow_shock_attached(), ans = NO\n");
		return NO;
	}
	debug_print("is_bow_shock_attached",
		"Left is_bow_shock_attached(), ans = YES\n");

	return YES;
}		/*end is_bow_shock_attached*/
#endif /* defined(TWOD) */
