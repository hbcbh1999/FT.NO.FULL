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


#if defined(TWOD)
/*
*
*				gnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Contains the main driver used in the first advance of the nodes:
*
*	g_node_propagate() propagates a node and some of the
*		neighboring points on the curves associated with the node.
*
*	and several of its support subroutines.
*	All access to this functions is through function pointers
*	in the Front structure that are set by init_physics() in ginit.c.
*
*	Curves are classified according to their wave_type and their start/
*	end status. The wave type is further divided into two general
*	categories, namely boundary curves and physical curves.
*	Nodes are classified into types, according to the types
*	and status of the curves at the node. We make this classification
*	precise in the following table.
*
*	NODE_TYPE	WAVE_TYPE		STATUS		COMMENTS
*
*	PASSIVE_NODE	PASSIVE_BOUNDARY	PASSIVE	out of comp. region
*
*	FIXED_NODE	boundary curve		FIXED	no change of type
*			boundary curve		FIXED		or angle
*			PASSIVE_BOUNDARY	PASSIVE	any number of these
*
*	CLOSED_NODE	physical curve		INCIDENT node of closed curve
*			physical curve		INCIDENT
*
*	NEUMANN_NODE	NEUMANN_BOUNDARY	FIXED	no change of angle
*			NEUMANN_BOUNDARY	FIXED
*			PASSIVE_BOUNDARY	PASSIVE	any number of these
*			physical curve		INCIDENT
*
*	DIRICHLET_NODE	DIRICHLET_BOUNDARY	FIXED	no change of angle
*			DIRICHLET_BOUNDARY	FIXED
*			PASSIVE_BOUNDARY	PASSIVE	any number of these
*			physical curve		INCIDENT
*
*	SUBDOMAIN_NODE	SUBDOMAIN_BOUNDARY	FIXED	no change of type
*			SUBDOMAIN_BOUNDARY	FIXED	or angle
*			PASSIVE_BOUNDARY	PASSIVE	any number of these
*			physical curve		INCIDENT
*
*	ATTACHED_B_NODE	boundary curve		FIXED	change of type or
*			boundary curve		FIXED	angle required
*			PASSIVE_BOUNDARY	PASSIVE	any number of these
*			physical curve		INCIDENT second phys curve is
*			physical curve		INCIDENT optional
*
*	B_REFLECT_NODE	boundary curve		FIXED	no change of type or
*			boundary curve		FIXED	angle allowed
*			PASSIVE_BOUNDARY	PASSIVE	any number of these
*			shock			INCIDENT
*			shock			REFLECTED
*
*	MACH_NODE	shock			INCIDENT
*			shock			REFLECTED
*			shock			MACH_STEM
*			CONTACT			SLIP
*
*	CROSS_NODE	shock			INCIDENT
*			shock			INCIDENT
*			shock			REFLECTED
*			shock			REFLECTED
*
*	OVERTAKE_NODE	shock			INCIDENT
*			shock			INCIDENT
*			shock			TRANSMITTED
*
*	DIFFRACTION_NODEshock			INCIDENT
*			shock			TRANSMITTED
*			shock			REFLECTED
*			rarefaction		REFLECTED
*			rarefaction		REFLECTED
*			CONTACT			SLIP
*			CONTACT			SLIP
*
*	TOT_INT_REFL_NODEshock			INCIDENT
*			rarefaction		REFLECTED
*			rarefaction		REFLECTED
*			CONTACT			SLIP
*			CONTACT			SLIP
*
*	TRANSMISSION_NODEshock			INCIDENT
*			shock			TRANSMITTED
*			CONTACT			SLIP
*			CONTACT			SLIP
*
*	CC_NODE		CONTACT			SLIP or INCIDENT
*			CONTACT			SLIP
*			CONTACT			SLIP
*
*	WAVE_END_NODE	any physics type 	INCIDENT
*	WAVE_END_NODE	any bdry type	FIXED or PASSIVE according to wave type
*/


#include <gdecs/gdecs.h>



/*
*                      g_node_propagate():
*
*       This is the high level routine to propagate a gas node.
*	It allows topology of the interface to be changed during propagation.
*
*	The routine calls different subroutines for different kinds of nodes,
*	e.g. closed_node_propagate() for node_type(oldn) == CLOSED_NODE.
*	In those subroutines, the node might change its type and the local
*	geometry around the node might also change because of the interaction.
*	The new	local states and local geometry around a node are calculated
*	there by the  Rankine-Hugoniot conditions, that is, from shock polar
*	analysis.  These local states record the waves produced due to
*	the interaction at the node.  These waves are different from
*	those waves scattered from a curve.  The calculated local states will
*	be stored in the states on the curves.
*/

/* ARGSUSED */
EXPORT int g_node_propagate(
	Front		*fr,
	POINTER		p2wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag,
	POINTER		user)
{
	Wave		*wave = (Wave*)p2wave;
	int		status = ERROR_NODE;
	int		i, dim = fr->rect_grid->dim;

	debug_print("node_propagate","Entered g_node_propagate()\n");

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("node_propagate")) print_node(oldn);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (oldn != NULL && oldn->in_curves==NULL && oldn->out_curves==NULL) 
	{
	    status = GOOD_NODE;
	    goto leave;
	}

		/* Propagate node according to its type */

	switch (node_type(newn)) 
	{
	case PASSIVE_NODE:
	case FIXED_NODE:
	case SUBDOMAIN_NODE:
	case CLOSED_NODE:
	case DIRICHLET_NODE:
	    status = f_node_propagate(fr,p2wave,oldn,newn,rp,dt,dt_frac,
	    			flag,user);
	    break;
	case NEUMANN_NODE: 
	{
	    void (*save_impose_bc)(POINT*,BOND*,CURVE*,double*,Front*,
				   boolean,boolean);
	    CURVE *cphys;
	    ORIENTATION cphys_orient;
	    cphys = find_physical_curve_at_node(oldn,&cphys_orient);
	    save_impose_bc = fr->impose_bc;

	    if (is_shock_wave(wave_type(cphys)))
	        fr->impose_bc = f_impose_bc;

	    status = f_node_propagate(fr,p2wave,oldn,newn,rp,dt,dt_frac,
	    			flag,user);
	    fr->impose_bc = save_impose_bc;
	    break;
	}

	case B_REFLECT_NODE:
	    status = B_reflect_node_propagate(fr,wave,oldn,newn,rp,
					      dt,dt_frac,flag);
	    break;

	case MACH_NODE:
	    status = Mach_node_propagate(fr,wave,oldn,newn,rp,dt,dt_frac,flag);
	    break;

	case ATTACHED_B_NODE:
	    status = attached_b_node_propagate(fr,wave,oldn,newn,rp,
					       dt,dt_frac,flag);
	    break;

#if defined(FULL_PHYSICS)
	case CROSS_NODE:
	    status = cross_node_propagate(fr,wave,oldn,newn,rp,dt,dt_frac,flag);
	    break;

	case OVERTAKE_NODE:
	    status = overtake_node_propagate(fr,wave,oldn,newn,rp,
					     dt,dt_frac,flag);
	    break;

	case DIFFRACTION_NODE:
	case TOT_INT_REFL_NODE:
	    status = diffraction_node_propagate(fr,wave,oldn,newn,
					        rp,dt,dt_frac,flag,user);
	    break;
	case TRANSMISSION_NODE:
	    status = transmission_node_propagate(fr,wave,oldn,newn,
					         rp,dt,dt_frac,flag);
	    break;

	case CC_NODE:
	    status = cc_node_propagate(fr,wave,oldn,newn,rp,dt,dt_frac,flag);
	    break;

	case WAVE_END_NODE:
	    status = wave_end_propagate(fr,oldn,newn,rp,dt,dt_frac,flag);
	    break;
#endif /* defined(FULL_PHYSICS) */

	default:
	    screen("ERROR in g_node_propagate(), "
	           "unknown node type %d\n",node_type(newn));
	    clean_up(ERROR);
	    break;
	}

leave:
#if defined(CHECK_FOR_BAD_STATES)
	if ((status == GOOD_NODE) &&
	    (debugging("bad_state") &&
	     (is_bad_state_at_node("g_node_propagate",newn))))
	{
	    screen("ERROR in g_node_propagate(), bad state at newn\n");
	    print_node(newn);
	    print_interface(newn->interface);
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
	if ((status == GOOD_NODE) &&
		is_physical_node(newn) &&
		(propagation_status(newn) == PROPAGATED_NODE) &&
		(Apply_CFL_at_nodes(fr) == YES))
	{
	    double *h = fr->rect_grid->h;
	    int i, dim = fr->rect_grid->dim;

	    for (i = 0; i < dim; i++)
	    {
	    	set_max_front_speed(i,fabs(Node_vel(newn)[i]),
				    return_obst_state(),
				    Coords(newn->posn),fr);
	    }
	    set_max_front_speed(dim,scaled_hypot(Node_vel(newn),h,dim),
				return_obst_state(),Coords(newn->posn),fr);
	}
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("node_propagate"))
	    print_node_status("status = ",status,"\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	debug_print("node_propagate","Left g_node_propagate(), \n");
	return status;
}	/*end g_node_propagate*/

EXPORT	boolean	is_bad_state_at_node(
	const char *function,
	NODE	*node)
{
	CURVE    **c;

	for (c = node->in_curves; c && *c; c++)
	{
	    if (is_bad_state(left_end_state(*c),YES,function))
	    {
		(void) printf("WARNING in is_bad_state_at_node(), "
			      "left_end_state of curve %llu is bad\n",
			      curve_number(*c));
	        if (debugging("bad_state"))
		    fprint_raw_gas_data(stdout,left_end_state(*c),
					node->interface->dim);
		return YES;
	    }
	    if (is_bad_state(right_end_state(*c),YES,function))
	    {
		(void) printf("WARNING in is_bad_state_at_node(), "
			      "right_end_state of curve %llu is bad\n",
			      curve_number(*c));
	        if (debugging("bad_state"))
		    fprint_raw_gas_data(stdout,right_end_state(*c),
					node->interface->dim);
		return YES;
	    }
	}
	for (c = node->out_curves; c && *c; c++)
	{
	    if (is_bad_state(left_start_state(*c),YES,function))
	    {
		(void) printf("WARNING in is_bad_state_at_node(), "
			      "left_start_state of curve %llu is bad\n",
			      curve_number(*c));
	        if (debugging("bad_state"))
		    fprint_raw_gas_data(stdout,left_start_state(*c),
					node->interface->dim);
		return YES;
	    }
	    if (is_bad_state(right_start_state(*c),YES,function))
	    {
		(void) printf("WARNING in is_bad_state_at_node(), "
			      "right_start_state of curve %llu is bad\n",
			      curve_number(*c));
	        if (debugging("bad_state"))
		    fprint_raw_gas_data(stdout,right_start_state(*c),
					node->interface->dim);
		return YES;
	    }
	}
	return NO;
}		/*end is_bad_state_at_node*/

EXPORT	void	node_warning(
	const char *fname,
	const char *mesg,
	const char *end,
	NODE_FLAG  flag)
{
	if (node_warnings_off(flag) == YES)
	    return;
	(void) printf("WARNING in %s(), %s%s",fname,mesg,end);
}		/*end node_warning*/
#endif /* defined(TWOD) */
