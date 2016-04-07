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
*				gssnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Contains the routine used in the first advance of nodes
*	corresponding to shock-shock wave interactions.
*/

#if defined(TWOD)

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	int	correct_overtake_for_subsonic(Front*,Wave*,double,RP_DATA*,
					      NODE*,O_CURVE**,O_CURVE**,
					      BOND**,POINT*,double*,
					      boolean*,boolean);
LOCAL	int	find_curves_at_shock_crossing(NODE*,NODE*,O_CURVE**,O_CURVE**,
					      Front*,ANGLE_DIRECTION*,double**,
					      double*,double*);
LOCAL	int	new_cross_position(O_CURVE**,O_CURVE**,BOND**,POINT**,
				   RPROBLEM**,Front*,Wave*,RP_DATA*,double*,
				   double,double*,double*,double*,boolean*,
				   int*,ANGLE_DIRECTION,NODE_FLAG);
LOCAL	int	set_states_about_shock_crossing(NODE*,O_CURVE**,O_CURVE**,
						BOND**,RP_DATA*,Locstate*,
						Locstate*,Front*,Wave*,
						POINT**,double**,double*,double*,
						double,ANGLE_DIRECTION,
						int,boolean*,int*,NODE_FLAG);
LOCAL	void	check_for_incomplete_cnode_deprecursor(int,NODE*,NODE*,
						       O_CURVE**,O_CURVE**,
						       double,double,Front*,
						       Wave*,RPROBLEM**);
LOCAL	void	degenerate_cross_node_prop(O_CURVE**,POINT**,ANGLE_DIRECTION,
					   RP_DATA*,double*,Front*,
					   Wave*,double);
LOCAL	void	modify_Mach_node(POINT*,BOND*,BOND*,NODE*,O_CURVE*,O_CURVE*,
				 O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
				 O_CURVE*,int,Front*,Wave*,double,
				 RP_DATA*,NODE_FLAG);
LOCAL	void	normal_at_degenerate_node(CURVE*,ORIENTATION,
					  CURVE*,ORIENTATION,
					  double*,Front*);
LOCAL	void	untrack_Mach_node(O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
				  O_CURVE*,O_CURVE*,ANGLE_DIRECTION,Front*,
				  Wave*,double,RPROBLEM*);
LOCAL	void	untrack_overtake_node(NODE*,O_CURVE**,O_CURVE**,ANGLE_DIRECTION,
				      Front*,Wave*,double,RPROBLEM*);

#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void	print_cross_node_data(NODE*,NODE*,O_CURVE**,O_CURVE**,
				      ANGLE_DIRECTION,double**,double*,double*);
LOCAL	void	print_Mach_node_curves(NODE*,O_CURVE*,O_CURVE*,O_CURVE*,
				       O_CURVE*,int,int,const char*);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

/*
*			Mach_node_propagate():
*
*					|cinc
*	        behind_comp; state1	|
*					|
*			____________ 	|
*		      /		     \  |
*	            /		      \ |
*	          /		       \|
*	    crefl		      /	\
*		    refl_comp;       /	 \      ahead_comp; state0;
*		    state2          /	  \
*			           /	   \
*				cslip	    \
*					     \ cmach
*				  stem_comp;
*				  state3
*/

EXPORT int Mach_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	BOND		*crossbinc, *crossbmach;/* intersecting two bonds */
						/* on Newcinc, Newcmach */
	COMPONENT	ahead_comp, prop_comp;
	POINT		Pc;			/* cross point */
	O_CURVE		Oldcinc, Newcinc;	/* the incident curve */
	O_CURVE		Oldcrefl, Newcrefl;	/* the reflected curve */
	O_CURVE		Oldcmach, Newcmach;	/* the Mach stem */
	O_CURVE		Oldcslip, Newcslip;	/* the slip curve */
	RP_DATA		*RP;
	int		status;
	ANGLE_DIRECTION	i_to_a_dir;		/* angular dir cinc to cmach */
	SIDE		inc_side, prop_side;
	ANGLE_DIRECTION	i_to_prop_dir;
	int		tracked_refl = YES;	/* is refl wave tracked? */

	debug_print("Mach_node","Entered Mach_node_propagate()\n");
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("Mach_node")) 
	{
	    (void) printf("\n\tOLD NODE:\n");
	    print_node(oldn);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    /* Identify curves */

	find_curve_with_status(oldn,&Oldcinc.curve,&Oldcinc.orient,INCIDENT);
	find_propagation_orientation(fr,(POINTER)wave,oldn,newn,NULL,
	                             &Oldcinc,dt,
	                             &i_to_prop_dir,&Oldcmach,&Newcmach,
				     &Oldcrefl,&Newcrefl,&inc_side,&prop_side,
				     &ahead_comp,&prop_comp);
	find_adjacent_curves(&Oldcinc,&i_to_a_dir,&Oldcmach,&Oldcrefl,
		             &ahead_comp);
	if (i_to_prop_dir != i_to_a_dir) 
	{
	    /*TODO: Correct routine */
	    screen("ERROR in mach_node_propagate(), "
		   "i_to_prop_dir != i_to_a_dir, CODE NEEDED\n");
	    clean_up(ERROR);
	}
	find_curve_with_status(oldn,&Oldcslip.curve,&Oldcslip.orient,SLIP);
	Check_return(
	    find_correspond_of_oriented_curve(&Oldcinc,&Newcinc,
					      newn,fr,newn->interface),
	    Mach_node_propagate)
	if ((Oldcrefl.curve == NULL)
	    	 ||
	    ((start_status(Oldcrefl.curve) != REFLECTED) &&
	     (end_status(Oldcrefl.curve) != REFLECTED)))
	{
	    /* If the reflected wave is untracked, the adjacent curve
	     * will still be found.  Above test checks for this */

	    Oldcrefl.curve = Newcrefl.curve = NULL;
	    tracked_refl = NO;
	}
	else
	{
	    if (!find_correspond_of_oriented_curve(&Oldcrefl,&Newcrefl,
	    				              newn,fr,newn->interface))
	    	tracked_refl = NO;
	}
	Check_return(
	    find_correspond_of_oriented_curve(&Oldcmach,&Newcmach,newn,fr,
					      newn->interface),
	    Mach_node_propagate)
	if (Oldcslip.curve != NULL)
	{
	    Check_return(
		find_correspond_of_oriented_curve(&Oldcslip,&Newcslip,newn,fr,
						  newn->interface),
		Mach_node_propagate)
	}

	RP = allocate_RP_DATA_structure(fr->sizest,YES,GAS_STATE);
	Rp_data(newn) = RP;
	RP->ang_dir = Opposite_ang_dir(i_to_a_dir);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("Mach_node")) 
	{
	    (void) printf("ahead_comp = %d\n",ahead_comp);
	    print_Mach_node_curves(oldn,&Oldcinc,&Oldcrefl,&Oldcmach,
	    	                   &Oldcslip,tracked_refl,
	    	                   (Oldcslip.curve != NULL &&
				   Newcslip.curve != NULL) ? YES : NO,"OLD");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */


	    /* Calculate the node speed and the states */

	status = find_Mach_node_states(oldn,&Oldcinc,&Newcinc,&Oldcmach,
				       &Newcmach,&Oldcslip,&crossbinc,
				       &crossbmach,&Pc,wave,fr,rp,dt,
				       dt_frac,RP,flag);
	if (status == MODIFY_TIME_STEP_NODE) return status;
	if (status != GOOD_NODE)
	{
	    (void) printf("WARNING, find_Mach_node_states() failed\n");
	    if (untrack_node(fr,MACH_NODE) == YES)
	    {
	    	(void) printf("\tuntracking Mach node\n");
	    	untrack_Mach_node(&Newcinc,&Newcmach,&Oldcrefl,&Newcrefl,
	    			  &Oldcslip,&Newcslip,i_to_a_dir,fr,wave,
	    			  dt,*rp);
	    	return GOOD_NODE;
	    }
	    else
	    	return status;
	}

	    /* Modify the curves near the node and assign the new states */

	modify_Mach_node(&Pc,crossbinc,crossbmach,newn,&Oldcinc,&Newcinc,
		         &Oldcrefl,&Newcrefl,&Oldcslip,&Newcslip,&Oldcmach,
			 &Newcmach,tracked_refl,fr,wave,dt,RP,flag);

	propagation_status(newn) = PROPAGATED_NODE;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("Mach_node")) 
	{
	    print_Mach_node_curves(newn,&Newcinc,&Newcrefl,&Newcmach,
			           &Newcslip,tracked_refl,
			           (Oldcslip.curve != NULL &&
				   Newcslip.curve != NULL) ? YES : NO,"NEW");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	debug_print("Mach_node","Left Mach_node_propagate()\n");
	return status;
}		/*end Mach_node_propagate*/


/*
*			       	modify_Mach_node():
*
*	   This routine finds the position and states in the neighborhood of a
*	Mach node.
*	   It is given the new position of the node (at the point PC), the in-
*	cident bond and the Mach bond together with the updated states and
*	angles of the reflected shock and slip line. After shifting the node
*	position to PC and ft_assigning the states to the new Mach bond and the
*	incident bond (an extra piece of which needs to be cut off first), we
*	turn to the old reflected shock and slipline. Since their node has been
* 	shifted, we propagate these curves near the node, insert a new bond
*	to maintain the correct reflection angle and assign their states.
*/

LOCAL	void  modify_Mach_node(
	POINT		*pc,
	BOND		*newbinc,
	BOND		*newbmach,
	NODE		*newn,
	O_CURVE		*oldcinc,
	O_CURVE		*newcinc,
	O_CURVE		*oldcrefl,
	O_CURVE		*newcrefl,
	O_CURVE		*oldcslip,
	O_CURVE		*newcslip,
	O_CURVE		*oldcmach,
	O_CURVE		*newcmach,
	int		tracked_refl,
	Front		*fr,
	Wave		*wave,
	double		dt,
	RP_DATA		*RP,
	NODE_FLAG	flag)
{
	NODE		*oldn;
	Locstate	st_left_inc, st_right_inc;
	Locstate	st_left_mach, st_right_mach;
	Locstate	st_left_refl, st_right_refl;
	Locstate	st_left_slip, st_right_slip;
	ANGLE_DIRECTION	i_to_a_dir = Opposite_ang_dir(RP->ang_dir);

	debug_print("Mach_node","Entered modify_Mach_node()\n");

	    /* Find states on incident curve */

	if (curve_ang_oriented_l_to_r(i_to_a_dir,newcinc->orient))
	{
	    st_left_inc = RP->state[1];
	    st_right_inc = RP->state[0];
	}
	else 
	{
	    st_left_inc = RP->state[0];
	    st_right_inc = RP->state[1];
	}

	    /* Find states on reflected curve */

	if (tracked_refl)
	{
	    if (curve_ang_oriented_l_to_r(i_to_a_dir,newcrefl->orient))
	    {
	    	st_left_refl = RP->state[2];
	    	st_right_refl = RP->state[1];
	    }
	    else 
	    {
	    	st_left_refl = RP->state[1];
	    	st_right_refl = RP->state[2];
	    }
	}

	    /* Find states on slip */

	if (curve_ang_oriented_l_to_r(i_to_a_dir,newcslip->orient))
	{
	    st_left_slip = RP->state[3];
	    st_right_slip = RP->state[2];
	}
	else 
	{
	    st_left_slip = RP->state[2];
	    st_right_slip = RP->state[3];
	}

	    /* Find states on mach stem */

	if (curve_ang_oriented_l_to_r(i_to_a_dir,newcmach->orient))
	{
	    st_left_mach = RP->state[0];
	    st_right_mach = RP->state[3];
	}
	else
	{
	    st_left_mach = RP->state[3];
	    st_right_mach = RP->state[0];
	}

	    /* Shift the node position */

	Coords(newn->posn)[0] = Coords(pc)[0];
	Coords(newn->posn)[1] = Coords(pc)[1];

	    /* Add points using Mach angle and states for these points */

	oldn = Node_of_o_curve(oldcmach);

	    /* Modify Mach curve */

	(void) propagate_curve_near_node(oldn,newn,oldcmach,newcmach,newbmach,
	                                 st_left_mach,st_right_mach,YES,
					 RP->ang[3],fr,wave,dt,flag);

	    /* Cut extra part of the incident curve and assign states */

	(void) propagate_curve_near_node(oldn,newn,oldcinc,newcinc,newbinc,
		                         st_left_inc,st_right_inc,NO,
					 RP->ang[0],fr,wave,dt,flag);

	    /* Propagate old reflected curve and insert a new bond */

	if (tracked_refl)
	    (void) propagate_curve_near_node(oldn,newn,oldcrefl,newcrefl,
	    	                             Bond_at_node_of_o_curve(newcrefl),
					     st_left_refl,st_right_refl,YES,
					     RP->ang[1],fr,wave,dt,flag);

	    /* Propagate old slipline and insert a new bond */

	if (oldcslip->curve != NULL)
	    (void) propagate_curve_near_node(oldn,newn,oldcslip,newcslip,
	    	                             Bond_at_node_of_o_curve(newcslip),
					     st_left_slip,st_right_slip,YES,
					     RP->ang[2],fr,wave,dt,flag);

	debug_print("Mach_node","Left modify_Mach_node()");
}		/*end modify_Mach_node*/

/*
*			ramp_reflection_corner_posn():
*
*	Stores the corner position of the ramp for future reference.
*	set_posn flags storing/retrieving:  YES to store a value, NO
*	to retrieve.  We return ERROR_FLOAT if a retrieval is requested
*	and the position has never been set.  This is in case this function
*	gets called in a problem other than a ramp reflection.
*/

EXPORT	void	ramp_reflection_corner_posn(
	double		*cposn,
	int		set_posn,
	int		dim)
{
	int		i;
	static int	is_set = NO;
	static double	corner_posn[MAXD];

	if (set_posn)
	{
	    is_set = YES;
	    for (i = 0; i < dim; i++)
		corner_posn[i] = cposn[i];
	}
	else if (is_set)
	{
	    for (i = 0; i < dim; i++)
		cposn[i] = corner_posn[i];
	}
	else
	    corner_posn[0] = ERROR_FLOAT;

}		/*end ramp_reflection_corner_posn*/


/*
*			untrack_Mach_node():
*
*	This algorithm provides some robustness by automatically untracking
*	a Mach node when the shock polar analysis fails.  Often this will be
*	caused by an increase in the incident angle, in which case the reflected
*	shock and slip line are becoming increasingly weak and untracking is
*	the most reasonable thing to do.
*
*	The reflected shock and slip line are untracked, with the incident
*	and Mach stem joined into a new curve.  The joining is done automatically
*	in the untracking code once there are only two curves remaining at the
*	node.  Recursive untracking is applied.
*/

LOCAL void untrack_Mach_node(
	O_CURVE		*newcinc,
	O_CURVE		*newcmach,
	O_CURVE		*oldcrefl,
	O_CURVE		*newcrefl,
	O_CURVE		*oldcslip,
	O_CURVE		*newcslip,
	ANGLE_DIRECTION	i_to_a_dir,
	Front		*fr,
	Wave		*wave,
	double		dt,
	RPROBLEM	*rp)
{
	COMPONENT    newcomp;
	NODE	     *oppn;
	UNTRACK_FLAG flag;
	boolean      oppn_ss;

	newcomp = (curve_ang_oriented_l_to_r(i_to_a_dir,newcinc->orient)) ?
	          negative_component(newcinc->curve) :
		  positive_component(newcinc->curve);

	if (newcrefl->curve != NULL)
	{
	    oppn = Opp_node_of_o_curve(newcrefl);
	    oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	    set_untrack_flag(flag,newcrefl->orient,YES,oppn_ss,YES,YES,YES);
	    (void) untrack_curve(newcrefl,oldcrefl,newcomp,dt,
			         fr,(POINTER)wave,rp,flag);
	    oldcrefl->curve = newcrefl->curve = NULL;
	}

	if (wave_type(newcmach->curve) != wave_type(newcinc->curve))
	    invert_curve(newcmach->curve);

	if (newcslip->curve != NULL)
	{
	    oppn = Opp_node_of_o_curve(newcslip);
	    oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	    set_untrack_flag(flag,newcslip->orient,NO,oppn_ss,YES,YES,YES);
	    (void) untrack_curve(newcslip,oldcslip,newcomp,dt,
			         fr,(POINTER)wave,rp,flag);
	    oldcslip->curve = newcslip->curve = NULL;
	}
}		/*end untrack_Mach_node*/


#if defined(DEBUG_NODE_PROPAGATE)
LOCAL void print_Mach_node_curves(
	NODE		*node,
	O_CURVE		*inc,
	O_CURVE		*refl,
	O_CURVE		*mach,
	O_CURVE		*slip,
	int		tracked_refl,
	int		tracked_slip,
	const char	*string)
{
	(void) printf("\nPRINTOUT OF CURVES AT %s MACH NODE\n",string);
	print_node(node);
	(void) printf("\t\t%s INCIDENT CURVE:\n",string);
	print_o_curve(inc);
	if (tracked_refl)
	{
	    (void) printf("\t\t%s REFLECTED CURVE:\n",string);
	    print_o_curve(refl);
	}
	else
	{
	    (void) printf("\n\t\t%s REFLECTED wave untracked\n",string);
	}
	(void) printf("\t\t%s MACH STEM\n",string);
	print_o_curve(mach);
	if (tracked_slip)
	{
	    (void) printf("\t\t%s SLIP CURVE\n",string);
	    print_o_curve(slip);
	}
	else
	{
	    (void) printf("\n\t\t%s SLIP CURVE untracked\n",string);
	}
	(void) printf("\nEND PRINTOUT OF CURVES AT %s MACH NODE\n",string);
}		/*end print_Mach_node_curves*/
#endif /* defined(DEBUG_NODE_PROPAGATE) */


#if defined(FULL_PHYSICS)
/*
*			  cross_node_propagate():
*
*
*
*               reflected wave \  state1   / 
*                   curve 1     \         /incident shock
*		state 2		 \	 /
*                                 \     /
*                                  \   /
*                                   \ /
*              ----------------------/     	state 0
* 	          contact           / \
*                   curve 2        /   \
*                                 / 	\
*                  state 3       /	 \
*                               /state 4  \
*	        reflected shock/ 	   \incident shock
*		   curve 3			curve 4
*
*	A cross node is a node at which two incident shocks collide
*	producing two reflected waves. 
*	The forward facing (low pressure) side of the two incident shocks
*	defines the front side of the wave.
*
*	Front and back refer to the low and high pressure sides of the 
*	incident shock, while ahead and behind refer to the direction
*	of node propagation.
*/

#define cnp_failed(name,status)						\
{									\
	if (debugging("cross_node"))					\
	{								\
	    (void) printf("WARNING in cross_node_propagate(), ");	\
	    (void) printf("%s() failed\n",name);			\
	    print_node_status("status = ",status,"\n");			\
	}								\
	debug_print("cross_node","Left cross_node_propagate()\n");		\
}


struct _POINT_LIST {
	POINT *p;
	struct _POINT_LIST *next;
	struct _POINT_LIST *prev;
};

typedef struct _POINT_LIST POINT_LIST;

/* ARGSUSED */
EXPORT int cross_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	BOND		*newb[5];
	Locstate	l_st_phys[2], r_st_phys[2];
	RP_DATA		*RP;
	double		ang0[5], ang4[5];
	double		tcr[2];		 /* fractional distances to cross */
	double		shift[MAXD];
	int		status;
	ANGLE_DIRECTION	i0_to_i4_dir;	/* dir of angle from incident curve 0*/
					/* to incident curve 4 */
	boolean		is_refl_raref;
	int		i;
	size_t		sizest = fr->sizest;
	int		num_curves;
	int		degenerate_cross_node = NO;
	static O_CURVE	*oldc[5], *newc[5];	/* Curves at node */
	static POINT	**pc = NULL;			/* cross point */
	static char	name[20];
	static double	**t = NULL;
	static boolean	correct_angle_at_node[5] = { NO, YES, YES, YES, NO };

	debug_print("cross_node","Entered cross_node_propagate()\n");


	    /* Allocate storage for RP_DATA */

	if (pc == NULL) 
	{
	    (void) strcpy(name,"RP->state[0]");
	    uni_array(&pc,2,sizeof(POINT *));
	    bi_array(&t,5,MAXD,FLOAT);
	    for (i = 0; i < 2; i++)
	    	pc[i] = Static_point(fr->interf);
	    for (i = 0; i < 5; i++)
	    {
	    	scalar(&oldc[i],sizeof(O_CURVE));
	    	scalar(&newc[i],sizeof(O_CURVE));
	    }
	}

	RP = allocate_RP_DATA_structure(sizest,YES,GAS_STATE);
	num_curves = find_curves_at_shock_crossing(oldn,newn,oldc,newc,fr,
						   &i0_to_i4_dir,t,ang0,ang4);
	RP->ang_dir = Opposite_ang_dir(i0_to_i4_dir);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("cross_node")) 
		print_cross_node_data(oldn,newn,oldc,newc,
				      i0_to_i4_dir,t,ang0,ang4);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

		/* Identify new position of node */

	status = new_cross_position(oldc,newc,newb,pc,rp,fr,wave,RP,ang0,
				    dt,dt_frac,shift,tcr,&is_refl_raref,
				    &degenerate_cross_node,i0_to_i4_dir,flag);

	if (status != GOOD_NODE)
	{
	    cnp_failed("new_cross_position",status);
	    propagation_status(newn) = VEL_COMPUTED_NODE;
	    check_for_incomplete_cnode_deprecursor(status,oldn,newn,
						   oldc,newc,dt,
						   *dt_frac,fr,wave,rp);
	    return status;
	}

	if (propagation_status(newn) != PROPAGATED_NODE)
	{
	    status = set_states_about_shock_crossing(newn,oldc,newc,newb,RP,
						     l_st_phys,r_st_phys,fr,
						     wave,pc,t,tcr,shift,dt,
						     i0_to_i4_dir,num_curves,
						     &is_refl_raref,
						     &degenerate_cross_node,
						     flag);

	    if (status != GOOD_NODE)
	    {
	    	cnp_failed("set_states_about_shock_crossing",status);
	    	return status;
	    }


	    if (!modify_curves_at_node(pc[0],newb,newn,oldc,newc,5,
	    			          correct_angle_at_node,fr,wave,
					  dt,RP,flag))
	    {
	    	status = ERROR_NODE;
	    	cnp_failed("modify_curves_at_node",status);
	    	return status;
	    }

	    propagation_status(newn) = PROPAGATED_NODE;
	    status = GOOD_NODE;
	}

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("cross_node")) 
	{
	    (void) printf("CURVES AT NEW CROSS NODE\n");
	    print_cross_node(newn,newc[0],newc[1],newc[2],
	    	             newc[3],newc[4]);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	debug_print("cross_node","Left cross_node_propagate(), ");
	if (debugging("cross_node"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end cross_node_propagate*/


LOCAL	int set_states_about_shock_crossing(
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	RP_DATA		*RP,
	Locstate	*l_st_phys,
	Locstate	*r_st_phys,
	Front		*fr,
	Wave		*wave,
	POINT		**pc,
	double		**t,
	double		*tcr,
	double		*shift,
	double		dt,
	ANGLE_DIRECTION	i0_to_i4_dir,
	int		num_curves,
	boolean		*is_refl_raref,
	int		*degenerate_cross_node,
	NODE_FLAG	flag)
{
	double		u[MAXD];
	double		cp;
	int		dim = fr->rect_grid->dim;
	boolean		is_plus_or;
	int		status;
	int		i;
	size_t		sizest = fr->sizest;
	static Locstate st01 = NULL, st04 = NULL;

	if (st04 == NULL)
	{
	    alloc_state(fr->interf,&st01,sizest);
	    alloc_state(fr->interf,&st04,sizest);
	}

	if (*degenerate_cross_node) return GOOD_NODE;

	if (curve_ang_oriented_l_to_r(i0_to_i4_dir,newc[0]->orient))
	{
	    right_state_along_bond(tcr[0],newb[0],newc[0]->curve,st01);
	    left_state_along_bond(tcr[0],newb[0],newc[0]->curve,RP->state[1]);
	    r_st_phys[0] = RP->state[0];
	    l_st_phys[0] = RP->state[1];
	}
	else 
	{
	    left_state_along_bond(tcr[0],newb[0],newc[0]->curve,st01);
	    right_state_along_bond(tcr[0],newb[0],newc[0]->curve,RP->state[1]);
	    l_st_phys[0] = RP->state[0];
	    r_st_phys[0] = RP->state[1];
	}

	if (curve_ang_oriented_l_to_r(i0_to_i4_dir,newc[4]->orient))
	{
	    left_state_along_bond(tcr[1],newb[4],newc[4]->curve,st04);
	    right_state_along_bond(tcr[1],newb[4],newc[4]->curve,RP->state[4]);
	    l_st_phys[1] = RP->state[0];
	    r_st_phys[1] = RP->state[4];
	}
	else 
	{
	    right_state_along_bond(tcr[1],newb[4],newc[4]->curve,st04);
	    left_state_along_bond(tcr[1],newb[4],newc[4]->curve,RP->state[4]);
	    r_st_phys[1] = RP->state[0];
	    l_st_phys[1] = RP->state[4];
	}
	interpolate_states(fr,0.5,0.5,Coords(newn->posn),st01,
			   Coords(newn->posn),st04,RP->state[0]);
	for (i = 0; i < dim; i++)
	    u[i] = vel(i,RP->state[0]) - Node_vel(newn)[i];
	(void) vector_product(u,t[0],&cp,dim);
	is_plus_or = ( cp > 0.0) ? NO : YES;
	
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("cross_node"))
	{
	    double M;
	    double *nod_v = Node_vel(newn);
	    static char name[20] = "RP->state[0]";
	    
	    (void) printf("node velocity = <%g, %g>\n",nod_v[0],nod_v[1]);
	    (void) printf("States in RP before ");
	    (void) printf("find_cross_node_states()\n");
	    for (i = 0; i < 5; i++)
	    {
	    	if (i == 2 || i == 3) continue;
	    	name[10] = '0' + i;
	    	M = mach_number(RP->state[i],nod_v);
	    	verbose_print_state(name,RP->state[i]);
	    	if (M >= 1.)
	    	{
	    	    (void) printf("RP->state[%d] is ",i);
	    	    (void) printf("supersonic, M%d = %g\n",i,M);
	    	}
		else
	    	{
	    	    (void) printf("RP->state[%d] is subsonic, ",i);
	    	    (void) printf("M%d = %g\n",i,M);
	    	}
	    }
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	if (!find_cross_node_states(Node_vel(newn),RP,
			               is_refl_raref,is_plus_or)) 
	{
	    *degenerate_cross_node = YES;
	
	/* Check for Degenerate cross node */

#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("cross_node"))
	    {
	    	double cos_ang;
	    	cos_ang = -scalar_product(t[0],t[4],dim);
	    	(void) printf("num_curves = %d, cos_ang = %g\n",
	    		      num_curves,cos_ang);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    if (num_curves != 2)
	    {
	    	status = ERROR_NODE;
	    	(void) printf("WARNING in cross_node_propagate(), "
	    	              "find_cross_node_states() failed\n"
	    	              "possible bifurcation, CODE NEEDED\n");
	    	return status;
	    }

	    if (end_of_curve(tcr[0],newb[0],newc[0]->curve,newc[0]->orient) &&
	        end_of_curve(tcr[1],newb[4],newc[4]->curve,newc[4]->orient) &&
	        (phys_virtuals_preset(flag) != YES))
	    {
	    	degenerate_cross_node_prop(oldc,pc,i0_to_i4_dir,
	    				   RP,shift,fr,wave,dt);
	    }
	    ft_assign(RP->state[2],RP->state[1],sizest);
	    ft_assign(RP->state[3],RP->state[1],sizest);
	    ft_assign(RP->state[4],RP->state[1],sizest);
	}
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("cross_node"))
	{
	    static char name[20] = "RP->state[0]";
	    (void) printf("States in RP after find_cross_node_states()\n");
	    for (i = 0; i < 5; i++)
	    {
	    	name[10] = '0' + i;
	    	verbose_print_state(name,RP->state[i]);
	    }
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (oldc[1]->curve != NULL && *is_refl_raref) 
	{
	    status = ERROR_NODE;
	    (void) printf("WARNING: in cross_node_propagate()\n"
	                  "The reflected shock has disappeared, it "
	                  "was present on last time step.\n"
	                  "Bifurcation has occured, code needed\n");
	    return status;
	}
	return GOOD_NODE;
}		/*end set_states_about_shock_crossing*/

/*
*		check_for_incomplete_cnode_deprecursor();
*
*	This function checks to make sure all relevant nodes are included
*	in the rp for a deprecursion from precursor with reflected 
*	rarefaction back to a single diffraction node.  We allow for the
*	case of untracked overtake node.  The cross in the name means
*	the deprecursion is being triggered by a cross node.  Things are
*	complicated by the fact that the numbering of the curves at a cross
*	node is somewhat arbitrary.
*/

LOCAL	void	check_for_incomplete_cnode_deprecursor(
	int		status,
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	double		dt,
	double		dt_frac,
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp)
{
	NODE		*interact_nodes[9];
	RECT_GRID	*gr;
	RP_NODE		*rpn;
	int		tir_node, c_node, t_node, o_node;
	int		t_index, tir_index, o_index;

	if (status != CROSS_NODE_NODE) return; 
	if (rp == NULL || *rp == NULL) return;
	tir_node = c_node = t_node = o_node = NO;
	for (rpn = (*rp)->first_rp_node; rpn != NULL; rpn = rpn->next)
	{
	    if (node_type(rpn->node) == TOT_INT_REFL_NODE)
	    	tir_node = YES;
	    if (node_type(rpn->node) == CROSS_NODE)
	    	c_node = YES;
	    if (node_type(rpn->node) == TRANSMISSION_NODE)
	    	t_node = YES;
	    if (node_type(rpn->node) == OVERTAKE_NODE)
	    	o_node = YES;
	}
	if (tir_node && c_node && t_node && o_node) return;

	gr = fr->rect_grid;
	if (c_node == NO) return;

	t_index = tir_index = o_index = ERROR;
	if (node_type(Opp_node_of_o_curve(oldc[0])) == TRANSMISSION_NODE)
	{
	    t_index = 0;
	    tir_index = 1;
	    o_index = 3;
	}
	else if (node_type(Opp_node_of_o_curve(oldc[4])) == TRANSMISSION_NODE)
	{
	    t_index = 4;
	    tir_index = 3;
	    o_index = 1;
	}
	if (t_index == ERROR) return;

	t_node = t_node || is_short_curve(oldc[t_index]->curve,
					  oldc[t_index]->orient,gr,1.0);

	o_node = o_node ||
	         (oldc[o_index]->curve      /* allow for untracked rarefaction*/
	     			 &&
	         (node_type(Opp_node_of_o_curve(oldc[o_index]))==OVERTAKE_NODE)
		 		 &&
	         is_short_curve(oldc[o_index]->curve,
			        oldc[o_index]->orient,gr,1.0));

	tir_node = tir_node ||
		   ((node_type(Opp_node_of_o_curve(oldc[tir_index])) ==
		       TOT_INT_REFL_NODE)
		 		 &&
		    is_short_curve(oldc[tir_index]->curve,
				   oldc[tir_index]->orient,gr,1.0));

	if (t_node && tir_node)
	{
	    interact_nodes[0] = Opp_node_of_o_curve(newc[tir_index]);
	    interact_nodes[1] = Opp_node_of_o_curve(oldc[tir_index]);
	    interact_nodes[2] = newn;
	    interact_nodes[3] = oldn;
	    interact_nodes[4] = Opp_node_of_o_curve(newc[t_index]);
	    interact_nodes[5] = Opp_node_of_o_curve(oldc[t_index]);
	    if (o_node)
	    {
	    	interact_nodes[6] = Opp_node_of_o_curve(newc[o_index]);
	    	interact_nodes[7] = Opp_node_of_o_curve(oldc[o_index]);
	    	interact_nodes[8] = NULL;
	    }
	    else
	    	interact_nodes[6] = NULL;

	    augment_rproblem_list(rp, interact_nodes,dt,dt_frac,
	    	                  oldn->interface,newn->interface,fr,
				  (POINTER)wave);
	}
}		/*end check_for_incomplete_cnode_deprecursor*/


#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void print_cross_node_data(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	ANGLE_DIRECTION	i0_to_i4_dir,
	double		**t,
	double		*ang0,
	double		*ang4)
{
	int		dim = newn->interface->dim;

	print_angle_direction("i0_to_i4_dir =",i0_to_i4_dir,"\n");
	(void) printf("CURVES AT OLD CROSS NODE\n");
	print_cross_node(oldn,oldc[0],oldc[1],oldc[2],
			 oldc[3],oldc[4]);
	(void) printf("CURVES AT NEW CROSS NODE\n");
	print_cross_node(newn,newc[0],newc[1],newc[2],
			 newc[3],newc[4]);
	(void) printf("\nTangents and angles at old cross node\n");
	(void) printf("OLD INCIDENT CURVE 0:\n");
	print_general_vector("t[0] = ",t[0],dim,"");
	(void) printf("ang0[0] = %g, ang4[0] = %g\n",ang0[0],ang4[0]);
	if (oldc[1]->curve) 
	{
	    (void) printf("OLD REFLECTED CURVE 1:\n");
	    print_general_vector("t[1] = ",t[1],dim,"");
	    (void) printf("ang0[1] = %g, ang4[1] = %g\n",ang0[1],ang4[1]);
	}
	if (oldc[2]->curve) 
	{
	    (void) printf("OLD CONTACT CURVE 2:\n");
	    print_general_vector("t[2] = ",t[2],dim,"");
	    (void) printf("ang0[2] = %g, ang4[2] = %g\n",ang0[2],ang4[2]);
	}
	if (oldc[3]->curve) 
	{
	    (void) printf("OLD REFLECTED CURVE 3:\n");
	    (void) printf("t[3] = %g, %g ",t[3][0],t[3][1]);
	    print_general_vector("t[3] = ",t[3],dim,"");
	    (void) printf("ang0[3] = %g, ang4[3] = %g\n",ang0[3],ang4[3]);
	}
	(void) printf("\t\tOLD INCIDENT CURVE 4:\n");
	(void) printf("t[4] = %g, %g ",t[4][0],t[4][1]);
	print_general_vector("t[4] = ",t[4],dim,"");
	(void) printf("ang0[4] = %g, ang4[4] = %g\n",ang0[4],ang4[4]);
}		/*end print_cross_node_data*/
#endif /* defined(DEBUG_NODE_PROPAGATE) */

LOCAL	int new_cross_position(
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		**pc,
	RPROBLEM	**rp,
	Front		*fr,
	Wave		*wave,
	RP_DATA		*RP,
	double		*ang0,
	double		dt,
	double		*dt_frac,
	double		*shift,
	double		*tcr,
	boolean		*is_refl_raref,
	int		*degenerate_cross_node,
	ANGLE_DIRECTION	i0_to_i4_dir,
	NODE_FLAG	flag)
{
	int		status;
	int		i;
	static const double ERR_ANG = 0.45; /*TOLERANCE*/

	if ((phys_virtuals_preset(flag) == YES) &&
	    ((i0_to_i4_dir == COUNTER_CLOCK && ang0[4] < PI - ERR_ANG) ||
	    (i0_to_i4_dir == CLOCKWISE && ang0[4] > PI + ERR_ANG)))
	{
	    *is_refl_raref = NO;
	    for (i = 0; i < 2; i++)
	    {
	    	newb[4*i] = Following_bond(Bond_at_node_of_o_curve(newc[4*i]),
					   newc[4*i]->orient);
	    	if (newc[4*i]->orient == POSITIVE_ORIENTATION)
	    	{
	    	    tcr[i] = 0.0;
	    	    Coords(pc[i])[0] = Coords(newb[4*i]->start)[0];
	    	    Coords(pc[i])[1] = Coords(newb[4*i]->start)[1];
	    	}
	    	else
	    	{
	    	    tcr[i] = 1.0;
	    	    Coords(pc[i])[0] = Coords(newb[4*i]->end)[0];
	    	    Coords(pc[i])[1] = Coords(newb[4*i]->end)[1];
	    	}
	    }
	    status = GOOD_NODE;
	}
	else if ((i0_to_i4_dir == COUNTER_CLOCK && ang0[4] < PI - ERR_ANG) ||
		(i0_to_i4_dir == CLOCKWISE && ang0[4] > PI + ERR_ANG))
	{

	    status = cross_or_extend_to_cross_two_propagated_curves(oldc[0],
				newc[0],oldc[4],newc[4],pc,&newb[0],&newb[4],
				&tcr[0],&tcr[1],fr,(POINTER)wave,rp,dt,dt_frac,
				flag,NULL);

	}
	else
	{
	    *degenerate_cross_node = YES;
	    degenerate_cross_node_prop(oldc,pc,i0_to_i4_dir,RP,
			               shift,fr,wave,dt);
	    for (i = 0; i < 2; i++)
	    	newb[4*i] = Bond_at_node_of_o_curve(newc[4*i]);
	    status = GOOD_NODE;
	}

	if (status != GOOD_NODE)
	    return status;
	    
	for (i = 1; i < 4; i++)
	{
	    newb[i] = (newc[i] && newc[i]->curve) ?
				Bond_at_node_of_o_curve(newc[i]) : NULL;
	}
	return status;
}		/*end new_cross_position*/



LOCAL	int find_curves_at_shock_crossing(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	Front		*fr,
	ANGLE_DIRECTION	*i0_to_i4_dir,
	double		**t,
	double		*ang0,
	double		*ang4)
{
	CURVE		*ctmp;
	CURVE		**c;
	ORIENTATION	ctmp_orient;
	double		sin0[5], cos0[5];
	double		sin4[5], cos4[5];
	double		t_tmp;
	int		num_curves;
	int		i;
	int		dim = fr->rect_grid->dim;

	for (i = 0; i < 5; i++)
	{
	    zero_scalar(newc[i],sizeof(O_CURVE));
	    zero_scalar(oldc[i],sizeof(O_CURVE));
	}

	identify_curves_with_status(newn,newc[0],newc[4],INCIDENT);
	Check_return(
	    find_correspond_of_oriented_curve(newc[0],oldc[0],oldn,fr,
					      oldn->interface),
	    find_curves_at_shock_crossing)
	Check_return(
	    find_correspond_of_oriented_curve(newc[4],oldc[4],oldn,fr,
					      oldn->interface),
	    find_curves_at_shock_crossing)
	
	if ((oldc[0]->orient == POSITIVE_ORIENTATION &&
			is_forward_wave(wave_type(oldc[0]->curve)))
			 ||
	    (oldc[0]->orient == NEGATIVE_ORIENTATION &&
			is_backward_wave(wave_type(oldc[0]->curve))))
	    *i0_to_i4_dir = CLOCKWISE;
	else
	    *i0_to_i4_dir = COUNTER_CLOCK;



		/* Identify other physical curves */

	num_curves = 0;
	for (c = oldn->in_curves; c && *c; c++)
	{
	    if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
	    	num_curves++;
	}
	for (c = oldn->out_curves; c && *c; c++)
	{
	    if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
	    	num_curves++;
	}
	

	/* Is the contact discontinuity tracked? */

	find_curve_with_status(oldn,&oldc[2]->curve,&oldc[2]->orient,SLIP);
	if (oldc[2]->curve != NULL)
	{
	    Check_return(
		find_correspond_of_oriented_curve(oldc[2],newc[2],newn,fr,
						  newn->interface),
		find_curves_at_shock_crossing)
	}

	/* Are any reflected shocks tracked? */

	identify_curves_with_status(oldn,oldc[1],oldc[3],REFLECTED);
	Check_return(
	    find_correspond_of_oriented_curve(oldc[1],newc[1],newn,fr,
		                              newn->interface),
	    find_curves_at_shock_crossing)
	Check_return(
	    find_correspond_of_oriented_curve(oldc[3],newc[3],newn,fr,
			                      newn->interface),
	    find_curves_at_shock_crossing)


	for (i = 0; i < 5; i++)
	{
	    if (oldc[i]->curve)
	    {
	    	double spacing = (i == 2) ? Front_spacing(fr,GENERAL_WAVE) :
				           Front_spacing(fr,VECTOR_WAVE);
	    	find_secant_to_curve(Node_of_o_curve(oldc[i])->posn,
				     Bond_at_node_of_o_curve(oldc[i]),
				     oldc[i]->curve,oldc[i]->orient,t[i],fr,
				     spacing);
	    }		
	}

	/* The reflected curves may need to be switched */

	for (i = 0; i < 5; i++)
	{
	    if (oldc[i]->curve)
	    {
	    	(void) vector_product(t[0],t[i],&sin0[i],dim);
	    	cos0[i] = scalar_product(t[0],t[i],dim);
	    	ang0[i] = angle(cos0[i],sin0[i]);
	    	(void) vector_product(t[4],t[i],&sin4[i],dim);
	    	cos4[i] = scalar_product(t[4],t[i],dim);
	    	ang4[i] = angle(cos4[i],sin4[i]);
	    }
	}
	if (oldc[1]->curve && oldc[3]->curve)
	{
	    if ((*i0_to_i4_dir == CLOCKWISE && ang0[1] > ang0[3]) ||
	    	(*i0_to_i4_dir == COUNTER_CLOCK && ang0[1] < ang0[3]))
	    {
	    	ctmp = oldc[1]->curve;
	    	ctmp_orient = oldc[1]->orient;
	    	oldc[1]->curve = oldc[3]->curve;
	    	oldc[1]->orient = oldc[3]->orient;
	    	oldc[3]->curve = ctmp;		
	    	oldc[3]->orient = ctmp_orient;
	    	for (i = 0; i < dim; i++)
	    	{
	    	    t_tmp = t[1][i]; t[1][i] = t[3][i];
	    	    t[3][i] = t_tmp;
	    	}
	    	t_tmp = ang0[1];
	    	ang0[1] = ang0[3];	
	    	ang0[3] = t_tmp;
	    	t_tmp = ang4[1];
	    	ang4[1] = ang4[3];	
	    	ang4[3] = t_tmp;
	    	ctmp = newc[1]->curve;
	    	ctmp_orient = newc[1]->orient;
	    	newc[1]->curve = newc[3]->curve;
	    	newc[1]->orient = newc[3]->orient;
	    	newc[3]->curve = ctmp;		
	    	newc[3]->orient = ctmp_orient;
	    }
	}
	else if (oldc[1]->curve && oldc[2]->curve)
	{
	    if ((*i0_to_i4_dir == CLOCKWISE && ang0[1] > ang0[2]) ||
	    	(*i0_to_i4_dir == COUNTER_CLOCK && ang0[1] < ang0[2]))
	    {
	    	oldc[3]->curve = oldc[1]->curve;
	    	oldc[3]->orient = oldc[1]->orient;
	    	newc[3]->curve = newc[1]->curve;
	    	newc[3]->orient = newc[1]->orient;
	    	for (i = 0; i < dim; i++)
	    	    t[3][i] = t[1][i];	
	    	ang0[3] = ang0[1];
	    	ang4[3] = ang4[1];
	    	oldc[1]->curve = NULL;
	    	newc[1]->curve = NULL;
	    }
	}
	else if (oldc[1]->curve)
	{
	    if ((*i0_to_i4_dir == CLOCKWISE && ang0[1] > 2.*PI - ang4[1])
	        		 ||
	    (*i0_to_i4_dir == COUNTER_CLOCK && ang4[1] < 2.*PI - ang0[1]))
	    {
	    	oldc[3]->curve = oldc[1]->curve;
	    	oldc[3]->orient = oldc[1]->orient;
	    	newc[3]->curve = newc[1]->curve;
	    	newc[3]->orient = newc[1]->orient;
	    	for (i = 0; i < dim; i++)
	    	    t[3][i] = t[1][i];	
	    	ang0[3] = ang0[1];
	    	ang4[3] = ang4[1];
	    	oldc[1]->curve = NULL;
	    	newc[1]->curve = NULL;
	    }
	}
	return num_curves;
}		/*end find_curves_at_shock_crossing*/


LOCAL void degenerate_cross_node_prop(
	O_CURVE		**oldc,
	POINT		**pc,
	ANGLE_DIRECTION	i0_to_i4_dir,
	RP_DATA		*RP,
	double		*shift,
	Front		*fr,
	Wave		*wave,
	double		dt)
{
	double		*coords0, *coords4;
	double		*nor, dn;
	double		coordsl0[MAXD], coordsr0[MAXD];
	double		coordsl4[MAXD], coordsr4[MAXD];
	double		V[MAXD];
	Locstate	sl0, sr0, sl4, sr4, sl, sr, ansl, ansr;
	int		i, j, dim = fr->interf->dim;
	static WSSten	*wssten = NULL;
	static Locstate sll0 = NULL, sll4 = NULL, srr0 = NULL, srr4 = NULL;

	debug_print("cross_node","DEGENERATE CROSS NODE\n");
	if (wssten == NULL) 
	{
	    wssten = AllocDefaultWSSten(fr);
	    alloc_state(fr->interf,&sll0,fr->sizest);
	    alloc_state(fr->interf,&sll4,fr->sizest);
	    alloc_state(fr->interf,&srr0,fr->sizest);
	    alloc_state(fr->interf,&srr4,fr->sizest);
	}
	else
	    ClearWSStenData(wssten);
	wssten->front = fr;
	wssten->wave = wave;
	wssten->ncomp = negative_component(oldc[0]->curve),
	wssten->pcomp = positive_component(oldc[0]->curve),
	wssten->pjump = 0.0;
	wssten->dt = dt;
	wssten->w_type = wave_type(oldc[0]->curve);
	wssten->hs = Hyper_surf(oldc[0]->curve);
	nor = wssten->nor;

	wssten->coords = coords0 = Coords(Node_of_o_curve(oldc[0])->posn);
	coords4 = Coords(Node_of_o_curve(oldc[4])->posn);
	normal_at_degenerate_node(oldc[0]->curve,oldc[0]->orient,
	                          oldc[4]->curve,oldc[4]->orient,nor,fr);
	wssten->dn = dn = grid_size_in_direction(nor,fr->rect_grid->h,dim);

	for (i = 0; i < dim; i++)
	{
	}
	sl0 = Left_state_at_node_of_o_curve(oldc[0]);
	sr0 = Right_state_at_node_of_o_curve(oldc[0]);
	sl4 = Left_state_at_node_of_o_curve(oldc[4]);
	sr4 = Right_state_at_node_of_o_curve(oldc[4]);
	sl = wssten->sl[0];
	sr = wssten->sr[0];
	if (oldc[0]->orient == oldc[4]->orient)
	{
	    interpolate_states(fr,0.5,0.5,coords0,sl0,coords4,sr4,sl);
	    interpolate_states(fr,0.5,0.5,coords0,sr0,coords4,sl4,sr);
	}
	else
	{
	    interpolate_states(fr,0.5,0.5,coords0,sl0,coords4,sl4,sl);
	    interpolate_states(fr,0.5,0.5,coords0,sr0,coords4,sr4,sr);

	}
	for (i = 1; i < wssten->nsts; i++)
	{
	    for (j = 0; j < dim; j++)
	    {
	        coordsl0[j] = wssten->lcrds[i][j] = coords0[j] - i*nor[j]*dn;
	        coordsr0[j] = wssten->rcrds[i][j] = coords0[j] + i*nor[j]*dn;
	        coordsl4[j] = coords4[j] - i*nor[j]*dn;
	        coordsr4[j] = coords4[j] + i*nor[j]*dn;
	    }
	    hyp_solution(coordsl0,negative_component(oldc[0]->curve),
	                 Hyper_surf(oldc[0]->curve),NEGATIVE_SIDE,fr,
			 wave,sll0,sl);
	    hyp_solution(coordsr0,positive_component(oldc[0]->curve),
	                 Hyper_surf(oldc[0]->curve),POSITIVE_SIDE,
			 fr,wave,srr0,sr);
	    hyp_solution(coordsl4,negative_component(oldc[4]->curve),
	                 Hyper_surf(oldc[4]->curve),NEGATIVE_SIDE,fr,
			 wave,sll4,sl);
	    hyp_solution(coordsr4,positive_component(oldc[4]->curve),
	                 Hyper_surf(oldc[4]->curve),POSITIVE_SIDE,
			 fr,wave,srr4,sr);

	    if (oldc[0]->orient == oldc[4]->orient)
	    {
	        interpolate_states(fr,0.5,0.5,coordsl0,sll0,coordsr4,srr4,
				   wssten->sl[i]);
	        interpolate_states(fr,0.5,0.5,coordsr0,srr0,coordsl4,sll4,
				   wssten->sr[i]);
	    }
	    else
	    {
	        interpolate_states(fr,0.5,0.5,coordsl0,sll0,coordsl4,sll4,
				   wssten->sl[i]);
	        interpolate_states(fr,0.5,0.5,coordsr0,srr0,coordsr4,srr4,
				   wssten->sr[i]);
	    }
	}
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("cross_node"))
	{
	    (void) printf("States in npt_w_speed() ");
	    (void) printf("in degenerate cross node case\n");
	    PrintWSSten(wssten);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	if (curve_ang_oriented_l_to_r(i0_to_i4_dir,oldc[0]->orient))
	{
	    ansr = RP->state[0];
	    ansl = RP->state[1];
	}
	else 
	{
	    ansl = RP->state[0];
	    ansr = RP->state[1];
	}
	npt_w_speed(wssten,ansl,ansr,V);
	for (i = 0; i < dim; i++)
	    Coords(pc[0])[i] = 0.5*(coords0[i]+coords4[i]+shift[i]) + V[i]*dt;
}		/*end degenerate_cross_node_prop*/

LOCAL void normal_at_degenerate_node(
	CURVE		*c1,
	ORIENTATION	c1_orient,
	CURVE		*c2,
	ORIENTATION	c2_orient,
	double		*nor,
	Front		*fr)
{
	double		tgnt[MAXD];

	tangent_at_degenerate_node(c1,c1_orient,c2,c2_orient,tgnt,fr);
	nor[0] = tgnt[1];
	nor[1] = -tgnt[0];
}		/*end normal_at_degenerate_node*/

EXPORT void tangent_at_degenerate_node(
	CURVE		*c1,
	ORIENTATION	c1_orient,
	CURVE		*c2,
	ORIENTATION	c2_orient,
	double		*tgnt,
	Front		*fr)
{
	double		*p1, *p2;
	double		*pt1, *pt2;
	double		t1[MAXD], t2[MAXD];
	double		len;
	int		i, dim = fr->rect_grid->dim;
	static	BOND	*b1 = NULL, *b2 = NULL;

	if (b1 == NULL)
	{
	    scalar(&b1,sizeof(BOND));
	    b1->start = Static_point(fr->interf);
	    b1->end = Static_point(fr->interf);
	    scalar(&b2,sizeof(BOND));
	    b2->start = Static_point(fr->interf);
	    b2->end = Static_point(fr->interf);
	}

	bond_tangent_to_curve(Node_of(c1,c1_orient)->posn,
		              Bond_at_node(c1,c1_orient),c1,c1_orient,b1,fr);
	pt1 = (c1_orient == POSITIVE_ORIENTATION) ? Coords(b1->end) :
						    Coords(b1->start);
	bond_tangent_to_curve(Node_of(c2,c2_orient)->posn,
	                      Bond_at_node(c2,c2_orient),c2,c2_orient,b2,fr);
	pt2 = (c2_orient == POSITIVE_ORIENTATION) ? Coords(b2->end) :
						    Coords(b2->start);
	p1 = Coords(Node_of(c1,c1_orient)->posn);
	p2 = Coords(Node_of(c2,c2_orient)->posn);
	for (i = 0; i < dim; i++)
	{
	    t1[i] = pt1[i] - p1[i];
	    t2[i] = pt2[i] - p2[i];
	}
	if (c1_orient == POSITIVE_ORIENTATION)
	{
	    for (i = 0; i < dim; i++)
		tgnt[i] = t1[i] - t2[i];
	}
	else
	{
	    for (i = 0; i < dim; i++)
		tgnt[i] = t2[i] - t1[i];
	}
	len = mag_vector(tgnt,dim);
	for (i = 0; i < dim; i++)
	    tgnt[i] /= len;
}		/*end tangent_at_degenerate_node*/

/*
*			  overtake_node_propagate():
*
*
*	reflected wave
* (curve 3)  reflected shock or
* (curves 2 and 4) leading and 
* trailing edges of reflected 
* rarefaction            \                       /
*		          \	incident shock1 /
*			   \         |         /
*			    \state 2 |state 1 /
*                            \       |       / 
*                             \      |      / incident shock0
*			       \     |     /
*				\    |    /
*		state 5		 \   |   /
*                                 \  |  /
*                                  \ | /
*                                   \|/
*              ----------------------/     	state 0
* 	          contact5          /
*                                  /
*                                 /
*                  state 6       /
*                               /
*	    transmitted shock6 / 
*
*	A overtake node is a node at which one incident shock overtakes
*	another from behind, producing a transmitted shock and a
*	reflected wave.
*	The forward facing (low pressure) side of the two incident shocks
*	defines the front side of the wave.
*
*	Front and back refer to the low and high pressure sides of the 
*	incident shock, while ahead and behind refer to the direction
*	of node propagation.
*/


/* ARGSUSED */
EXPORT int overtake_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	BOND             *newb[7];
	O_CURVE          Ctmp0, Ctmp1;
	double            tcr[2];        /* fractional distances to overtake */
	int              status;
	ANGLE_DIRECTION  i0_to_i1_dir;    /* ang dir from inc 0 to inc 1 */
	boolean          is_plus_or;
	boolean          is_refl_raref;
	int              i;
	static Locstate  tmpst = NULL;
	static O_CURVE   **oldc = NULL, **newc = NULL; /* Curves at node */
	static POINT     *pc = NULL;            /* overtake point */
	static RP_DATA   *RP = NULL;
	static boolean   correct_angle_at_node[7] =
	                     { NO, NO, YES, YES, YES, YES, YES };
#if defined(DEBUG_NODE_PROPAGATE)
	static char    header[7][100] = {"AHEAD INCIDENT CURVE",
	                  "OVERTAKING INCIDENT CURVE",
	                  "REFLECTED RAREFACTION LEADING EDGE",
	                  "REFLECTED SHOCK CURVE",
	                  "REFLECTED RAREFACTION TRAILING EDGE",
	                  "BACK CONTACT CURVE",
	                  "TRANSMITTED CURVE"};
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	debug_print("overtake","Entered overtake_node_propagate()\n");

	if (tmpst == NULL) 
	{
	    pc = Static_point(fr->interf);
	    uni_array(&oldc,7,sizeof(O_CURVE *));
	    uni_array(&newc,7,sizeof(O_CURVE *));
	    alloc_state(fr->interf,&tmpst,fr->sizest);
	    RP = allocate_RP_DATA_structure(fr->sizest,NO,GAS_STATE);
	    for (i = 0; i < 7; i++)
	    {
	        scalar(&oldc[i],sizeof(O_CURVE));
	        scalar(&newc[i],sizeof(O_CURVE));
	    }
	}

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("overtake")) 
	{
	    (void) printf("\n\tOLD NODE:\n");
	    print_node(oldn);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    /* Identify incident curves */

	find_curve_with_status(oldn,&oldc[0]->curve,&oldc[0]->orient,OVERTOOK);
	find_curve_with_status(oldn,&oldc[1]->curve,&oldc[1]->orient,INCIDENT);

	i0_to_i1_dir = ((oldc[0]->orient == POSITIVE_ORIENTATION &&
	                    is_forward_wave(wave_type(oldc[0]->curve)))
	                        ||
	                (oldc[0]->orient == NEGATIVE_ORIENTATION &&
	                    is_backward_wave(wave_type(oldc[0]->curve)))) ?
	               COUNTER_CLOCK : CLOCKWISE;
	RP->ang_dir = i0_to_i1_dir;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("overtake"))
	    print_angle_direction("i0_to_i1_dir = ",i0_to_i1_dir,"\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    /* Identify other curves */

	find_curve_with_status(oldn,&oldc[5]->curve,&oldc[5]->orient,SLIP);
	find_curve_with_status(oldn,&oldc[6]->curve,&oldc[6]->orient,
	                       TRANSMITTED);

	identify_curves_with_status(oldn,&Ctmp0,&Ctmp1,REFLECTED);
	if (Ctmp0.curve == NULL)
	{
	    /* Untracked reflected wave */
	    oldc[2]->curve = oldc[3]->curve = oldc[4]->curve = NULL;
	}
	else if (is_shock_wave(wave_type(Ctmp0.curve)))
	{
	    /* Tracked reflected shock */
	    *oldc[3] = Ctmp0;
	    oldc[2]->curve = oldc[4]->curve = NULL;
	}
	else 
	{
	    /* Tracked reflected rarefaction */
	    oldc[3]->curve = NULL;
	    if (is_rarefaction_leading_edge(wave_type(Ctmp0.curve)))
	    {
	        *oldc[2] = Ctmp0;
	        *oldc[4] = Ctmp1;
	    }
	    else
	    {
	        *oldc[2] = Ctmp1;
	        *oldc[4] = Ctmp0;
	    }
	}

	for (i = 0; i < 7; i++)
	{
	    Check_return(
		find_correspond_of_oriented_curve(oldc[i],newc[i],newn,
						  fr,newn->interface),
		overtake_node_propagate)
	}

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("overtake")) 
	{
	    for (i = 0; i < 7; i++)
	    {
	        if (oldc[i]->curve == NULL)
		    continue;
	        (void) printf("\t\tOLD %s:\n",header[i]);
	        if (debugging("states"))
	            verbose_print_curve_states(oldc[i]->curve);
	        else
	            print_o_curve(oldc[i]);
	    }
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    /* Identify new position of node */

	status = cross_or_extend_to_cross_two_propagated_curves(oldc[0],newc[0],
								oldc[1],newc[1],
								&pc,&newb[0],
								&newb[1],
								&tcr[0],
								&tcr[1],
								fr,
								(POINTER)wave,
								rp,dt,dt_frac,
								flag,NULL);
	
	if (status != GOOD_NODE) 
	{
	    debug_print("overtake","Left overtake_node_propagate(), status = %s\n",
		  node_status_as_string(status));
	    return status;
	}
	
	for (i = 2; i < 7; i++)
	    newb[i] = (newc[i]->curve) ? Bond_at_node_of_o_curve(newc[i]) :
					 NULL;


	if (is_backward_wave(wave_type(newc[0]->curve)))
	    left_state_along_bond(tcr[0],newb[0],newc[0]->curve,RP->state[0]);
	else 
	    right_state_along_bond(tcr[0],newb[0],newc[0]->curve,RP->state[0]);

	if (is_backward_wave(wave_type(newc[1]->curve)))
	{
	    left_state_along_bond(tcr[1],newb[1],newc[1]->curve,RP->state[1]);
	    right_state_along_bond(tcr[1],newb[1],newc[1]->curve,RP->state[2]);
	}
	else 
	{
	    right_state_along_bond(tcr[1],newb[1],newc[1]->curve,RP->state[1]);
	    left_state_along_bond(tcr[1],newb[1],newc[1]->curve,RP->state[2]);
	}


#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("overtake"))
	{
	    double M;
	    int i;
	    char mesg[80];

	    (void) printf("States on crossing bonds\n");
	    (void) printf("tcr[0] = %g\n",tcr[0]);
	    verbose_print_bond_states("newb[0]",newb[0],newc[0]->curve);
	    (void) printf("tcr[1] = %g\n",tcr[1]);
	    verbose_print_bond_states("newb[1]",newb[1],newc[1]->curve);
	    (void) printf("Node velocity = <%g, %g>\n",
	        Node_vel(newn)[0],Node_vel(newn)[1]);
	    (void) printf("\n\nInterpolated states\n");
	    for ( i = 0; i < 3; i++)
	    {
	        M = mach_number(RP->state[i],Node_vel(newn));
	        (void) sprintf(mesg,"RP->state[%d] is %s, M%d = %g",
	            i,(M >= 1.0) ? "supersonic" : "subsonic",i,M);
	        verbose_print_state(mesg,RP->state[i]);
	    }
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	clear_state(fr->interf,RP->state[3],fr->sizest);
	clear_state(fr->interf,RP->state[4],fr->sizest);
	ft_assign(tmpst,RP->state[1],fr->sizest);

	is_plus_or = (i0_to_i1_dir == CLOCKWISE) ? NO : YES;

	if (!find_overtake_node_states(Node_vel(newn),RP,
	          &is_refl_raref,is_plus_or))
	{
	    status = ERROR_NODE;
	    if (mach_number(RP->state[1],Node_vel(newn)) < 1.0)
	    {
	        status = correct_overtake_for_subsonic(fr,wave,dt,RP,
	                                               newn,oldc,newc,
						       newb,pc,tcr,
						       &is_refl_raref,
						       is_plus_or);

	    }

	    if (status != GOOD_NODE)
	    {
	        (void) printf("WARNING in overtake_node_propagate(), "
	                      "find_overtake_node_states() failed, "
	                      "possible bifurcation.\n");

	        /* Don't untrack if part of a precursor configuration.
	         * This test will probably eventually be inadequate for
	         * identifying this case.
	         */

	        if ((untrack_node(fr,OVERTAKE_NODE) == YES)
	                &&
	        (node_type(Opp_node_of_o_curve(newc[1])) != TOT_INT_REFL_NODE))
	        {
	            (void) printf("\tUntracking overtake node.\n");
	            untrack_overtake_node(newn,oldc,newc,i0_to_i1_dir,fr,
					  wave,dt,*rp);
	            status = GOOD_NODE;
	        }
	        debug_print("overtake","Left overtake_node_propagate(), ");
	        if (debugging("overtake"))
		    print_node_status("status = ",status,"\n");
	        return status;
	    }
	}

	if ((oldc[3]->curve != NULL && is_refl_raref) ||
	    (!is_refl_raref && (oldc[2]->curve || oldc[4]->curve)))
	{
	    (void) printf("WARNING in overtake_node_propagate(), "
	                  "Bifurcation has occured, code needed\n");
	    status = ERROR_NODE;
	    debug_print("overtake","Left overtake_node_propagate(), ");
	    if (debugging("overtake"))
		print_node_status("status = ",status,"\n");
	    return status;
	}

	if (!modify_curves_at_node(pc,newb,newn,oldc,newc,7,
	                              correct_angle_at_node,fr,wave,dt,RP,flag))
	{
	    (void) printf("WARNING in overtake_node_propagate(), "
	                  "modify_curves_at_node failed\n");
	    status = ERROR_NODE;
	    debug_print("overtake","Left overtake_node_propagate(), ");
	    if (debugging("overtake"))
		print_node_status("status = ",status,"\n");
	    return status;
	}

	status = GOOD_NODE;

	propagation_status(newn) = PROPAGATED_NODE;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("overtake")) 
	{
	    for (i = 0; i < 7; i++)
	    {
	        if (newc[i]->curve == NULL)
		    continue;
	        (void) printf("\t\tNEW %s:\n",header[i]);
	        if (debugging("states"))
	            verbose_print_curve_states(newc[i]->curve);
	        else
	            print_o_curve(newc[i]);
	    }
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	debug_print("overtake","Left overtake_node_propagate(), ");
	if (debugging("overtake"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end overtake_node_propagate*/


LOCAL	int correct_overtake_for_subsonic(
	Front		*fr,
	Wave		*wave,
	double		dt,
	RP_DATA		*RP,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		*pc,
	double		*tcr,
	boolean		*is_refl_raref,
	boolean		is_plus_or)
{
	POINT		*ps, *pe;
	double		tau;
	double		p1, c1, M1;
	double		V[MAXD], dx, dy;
	double		sign, len, cst;
	double		ang, theta;
	int		i, j, dim = fr->rect_grid->dim;
	int		status = ERROR_NODE;
	static const int	MAX_NUM_ITER = 5;/*TOLERANCE*/
	static	POINT	*ptmps = NULL, *ptmpe = NULL;

	debug_print("overtake","Entered correct_overtake_for_subsonic()\n");
	if (ptmpe == NULL)
	{
	    ptmps = Static_point(fr->interf);
	    ptmpe = Static_point(fr->interf);
	}

	if (newb[0] == newc[0]->curve->first)
	{
	    ps = ptmps;
	    point_propagate(fr,(POINTER)wave,Node_of_o_curve(oldc[0])->posn,ps,
			    Bond_at_node_of_o_curve(oldc[0]),
			    oldc[0]->curve,dt,V);
	}
	else
	    ps = newb[0]->start;
	if (newb[0] == newc[0]->curve->last)
	{
	    pe = ptmpe;
	    point_propagate(fr,(POINTER)wave,Opp_node_of_o_curve(oldc[0])->posn,
			    pe,Bond_at_opp_node_of_o_curve(oldc[0]),
			    oldc[0]->curve,dt,V);
	}
	else
	    pe = newb[0]->end;

	for (i = 0; i < MAX_NUM_ITER; i++)
	{
	    p1 = pressure(RP->state[1]);
	    c1 = sound_speed(RP->state[1]);
	    for (j = 0; j < dim; j++)
	    {
	        V[j] = (vel(j,RP->state[1]) - Node_vel(newn)[j])/c1;
	    }
	    M1 = mag_vector(V,dim); V[0] /= M1;	V[1] /= M1;
	    dx = Coords(pe)[0] - Coords(ps)[0];
	    dy = Coords(pe)[1] - Coords(ps)[1];
	    len = hypot(dx,dy);
	    cst = (V[0]*dx + V[1]*dy)/len;
	    sign = (cst > 0.0) ? -1.0 : 1.0;
	    len /= (c1*dt);
	    tau = sign*(sqrt(sqr(M1*cst)+1.0- sqr(M1)) - M1*fabs(cst))/len;
	    tcr[0] += tau;
	    Coords(pc)[0] += tau*dx;	Coords(pc)[1] += tau*dy;
	    Node_vel(newn)[0] += tau*dx/dt;	Node_vel(newn)[1] += tau*dy/dt;
	    if (debugging("overtake"))
	    {
	    	(void) printf("Sonic overtake correction factor %d = %g\n",
			      i,tau);
	    }
	    if (!s_polar_3(RP->state[0],YES,p1,is_plus_or,NO,Node_vel(newn),
			      RP->state[1],&ang,&theta))
	    {
	    	(void) printf("WARNING in correct_overtake_for_subsonic(), "
	    	              "s_polar_3() failed\n");
	    	return ERROR_NODE;
	    }
	}
	if (debugging("overtake"))
	{
	    (void) printf("Final sonic overtake correction factor = %g\n"
		          "Mach number RP->state[1] = %g\n",
			  tau,mach_number(RP->state[1],Node_vel(newn)));
	}
	if (is_rarefaction_wave(wave_type(newc[1]->curve)))
	    ft_assign(RP->state[2],RP->state[1],fr->sizest);
	if (find_overtake_node_states(Node_vel(newn),RP,
	                              is_refl_raref,is_plus_or))
	    status = GOOD_NODE;
	else
	{
	    status = ERROR_NODE;
	    (void) printf("WARNING in correct_overtake_for_subsonic() "
	                  "find_overtake_node_states() failed\n"
	                  "possible bifurcation, CODE NEEDED\n");
	}
	debug_print("overtake","Left correct_overtake_for_subsonic(), ");
	if (debugging("overtake"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end correct_overtake_for_subsonic*/


/*
*			untrack_overtake_node():
*
*	This routine provides some robustness to the code by untracking an
*	overtake node when the shock polar analysis fails.  All behind wave
*	and incident1 are untracked.  incident0 and the transmitted shock
*	are joined if possible, otherwise all waves at the node are untracked.
*/

LOCAL void untrack_overtake_node(
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	ANGLE_DIRECTION	i0_to_i1_dir,
	Front		*fr,
	Wave		*wave,
	double		dt,
	RPROBLEM	*rp)
{
	COMPONENT    newcomp;
	NODE	     *oppn;
	UNTRACK_FLAG flag;
	int	     i, num_in, num_out;
	boolean      oppn_ss;

	newcomp = (curve_ang_oriented_l_to_r(i0_to_i1_dir,newc[0]->orient)) ?
			positive_component(newc[0]->curve) :
			negative_component(newc[0]->curve);

	    /* untrack reflected waves and slip */

	for (i = 2; i < 6; i++)
	{
	    if (newc[i]->curve == NULL) continue;
	    oppn = Opp_node_of_o_curve(newc[i]);
	    oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	    set_untrack_flag(flag,newc[i]->orient,YES,oppn_ss,YES,YES,YES);
	    (void) untrack_curve(newc[i],oldc[i],newcomp,dt,fr,
				 (POINTER)wave,rp,flag);
	    oldc[i]->curve = newc[i]->curve = NULL;
	}

	if (newc[1]->curve->interface == NULL)
	    return;		/* i1 was recursively deleted */

	/* If there are still three curves at the node, both incidents and the 
	 * transmitted are still tracked.  Invert the transmitted shock for
	 * automatic join with i0 when i1 is untracked.  If join is not
	 * possible, all remaining waves will be untracked automatically
	 * by untracking i1.
	 */

	if ((num_curves_at_node(newn,&num_in,&num_out) == 3)
	    		 &&
	    (wave_type(newc[0]->curve) != wave_type(newc[6]->curve)))
	    invert_curve(newc[6]->curve);

	oppn = Opp_node_of_o_curve(newc[1]);
	oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	set_untrack_flag(flag,newc[1]->orient,YES,oppn_ss,YES,YES,YES);
	(void) untrack_curve(newc[1],oldc[1],newcomp,dt,fr,
			     (POINTER)wave,rp,flag);
	oldc[1]->curve = newc[1]->curve = NULL;
	return;
}		/*end untrack_overtake_node*/
#endif /* defined(FULL_PHYSICS) */
#endif /* defined(TWOD) */
