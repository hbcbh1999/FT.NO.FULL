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
*				gtypes.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains functions for setting, testing,  and copying
*	integer type flags,  including,  wave type, node type, state type, and
*	start and end status.
*
*	
*/


#include <gdecs/gdecs.h>

#if defined(TWOD)
	/* LOCAL function prototypes */
LOCAL	Gas_param	*find_params_on_curve(CURVE*,SIDE);
#endif /* defined(TWOD) */

EXPORT	int opposite_wave_type(
	int		w_type)
{
	switch (w_type)
	{
	case FORWARD_SHOCK_WAVE:
	    return BACKWARD_SHOCK_WAVE;
	case BACKWARD_SHOCK_WAVE:
	    return FORWARD_SHOCK_WAVE;
	case FORWARD_SOUND_WAVE_LE:
	    return BACKWARD_SOUND_WAVE_LE;
	case BACKWARD_SOUND_WAVE_LE:
	    return FORWARD_SOUND_WAVE_LE;
	case FORWARD_SOUND_WAVE_TE:
	    return BACKWARD_SOUND_WAVE_TE;
	case BACKWARD_SOUND_WAVE_TE:
	    return FORWARD_SOUND_WAVE_TE;
	default:
	    return w_type;
	}
}		/*end opposite_wave_type*/



#if defined(TWOD)
/*
*			set_status_at_node():
*/

EXPORT void set_status_at_node(
	CURVE		*c,
	ORIENTATION	c_orient,
	int		status)
{
	if (c == NULL) return;
	if (c_orient == POSITIVE_ORIENTATION)
		start_status(c) = status;
	else
		end_status(c) = status;
}		/*end set_status_at_node*/


/*
*			set_fixed_status():
*/

EXPORT void set_fixed_status(
	CURVE		*c,
	ORIENTATION	c_orient)
{
	if (c_orient == POSITIVE_ORIENTATION)
		start_status(c) = FIXED;
	else	end_status(c) = FIXED;
}		/*end set_fixed_status*/

/*
*			set_incident_status():
*/

EXPORT void set_incident_status(
	CURVE		*c,
	ORIENTATION	c_orient)
{
	if (c_orient == POSITIVE_ORIENTATION)
		start_status(c) = INCIDENT;
	else	end_status(c) = INCIDENT;
}		/*end set_incident_status*/

/*
*			set_start_end_status():
*
*	Sets the start and end status of curves in an interface from the
*	knowledge of the wave types and node types.  Only boundary type
*	statuses need to be set.  Physical status types are set by
*	the copy curve in zoom_interface() and should not be reset here.
*/

EXPORT void set_start_end_status(
	Front		*fr)
{
	INTERFACE	*intfc = fr->interf;
	CURVE		**c;
	NODE            **n;

	for (c = intfc->curves; *c; ++c)
	{
	    if (start_status(*c) == UNKNOWN_CURVE_STATUS)
	    {
	    	if (wave_type(*c) == PASSIVE_BOUNDARY)
	    	    start_status(*c) = PASSIVE;
	    	else if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE)
	    	    start_status(*c) = FIXED;
	    	else
	    	    start_status(*c) = INCIDENT;
	    }
	    if (end_status(*c) == UNKNOWN_CURVE_STATUS)
	    {
	    	if (wave_type(*c) == PASSIVE_BOUNDARY)
	    	    end_status(*c) = PASSIVE;
	    	else if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE)
	    	    end_status(*c) = FIXED;
	    	else
	    	    end_status(*c) = INCIDENT;
	    }
	}
        for (n = intfc->nodes; n && *n; ++n)
        {
            ORIENTATION orient;

            if (node_type(*n) == ATTACHED_B_NODE &&
                find_physical_curve_at_node(*n,&orient) == NULL)
            {
		node_type(*n) = UNKNOWN_NODE_TYPE;
		if (((*n)->in_curves == NULL) && ((*n)->out_curves == NULL))
		{
		    delete_node(*n);
		    n = intfc->nodes - 1;
		}
		else
		{
		    boolean subdomain = NO, neumann = NO, dirichlet = NO;
                    for (c = (*n)->in_curves; c && *c; ++c)
		    {
			switch (wave_type(*c))
			{
			case NEUMANN_BOUNDARY:
			    neumann = YES;
			    break;
			case DIRICHLET_BOUNDARY:
			    dirichlet = YES;
			    break;
			case SUBDOMAIN_BOUNDARY:
			    subdomain = YES;
			    break;
			default:
			    break;
			}
		    }
                    for (c = (*n)->out_curves; c && *c; ++c)
		    {
			switch (wave_type(*c))
			{
			case NEUMANN_BOUNDARY:
			    neumann = YES;
			    break;
			case DIRICHLET_BOUNDARY:
			    dirichlet = YES;
			    break;
			case SUBDOMAIN_BOUNDARY:
			    subdomain = YES;
			    break;
			default:
			    break;
			}
		    }
		    if (subdomain == YES)
			node_type(*n) = SUBDOMAIN_NODE;
		    else if (neumann == YES)
			node_type(*n) = NEUMANN_NODE;
		    else if (dirichlet == YES)
			node_type(*n) = DIRICHLET_NODE;
		}
            }
        }
}		/*end set_start_end_status*/

LOCAL	Gas_param *find_params_on_curve(
	CURVE *c,
	SIDE  side)
{
	BOND	  *b;
	INTERFACE *intfc = c->interface;
	Gas_param *param;
	Gas_param **plist = gas_params_list(intfc);
	int       i, np = num_gas_params(intfc);

	param = NULL;
	if (side == NEGATIVE_SIDE)
	{
	    if (is_excluded_comp(negative_component(c),intfc) == YES)
		return NULL;/*Obstacle state*/
	    if (!is_obstacle_state(left_start_state(c)))
	    	param = Params(left_start_state(c));
	    else if (!is_obstacle_state(left_end_state(c)))
	    	param = Params(left_end_state(c));
	    else
	    {
	    	for (b = c->first; b != c->last; b = b->next)
	    	{
	    	    if (!is_obstacle_state(left_state(b->end)))
	    	    {
	    	    	param = Params(left_state(b->end));
	    	    	break;
	    	    }
	    	}
	    }
	}
	else if (side == POSITIVE_SIDE)
	{
	    if (is_excluded_comp(positive_component(c),intfc) == YES)
		return NULL;/*Obstacle state*/
	    if (!is_obstacle_state(right_start_state(c)))
	    	param = Params(right_start_state(c));
	    else if (!is_obstacle_state(right_end_state(c)))
	    	param = Params(right_end_state(c));
	    else
	    {
	    	for (b = c->first; b != c->last; b = b->next)
	    	{
	    	    if (!is_obstacle_state(right_state(b->end)))
	    	    {
	    	        param = Params(right_state(b->end));
	    	    	break;
	    	    }
	    	}
	    }
	}
	for (i = 0; i < np; i++)
	    if (param == plist[i])
		break;
	if (i == np)
	{
	    screen("ERROR in find_params_on_curve(), param %llu not in "
		   "params list\n",gas_param_number(param));
	    print_curve(c);
	    print_interface(intfc);
	    show_intfc_states(intfc);
	    clean_up(ERROR);
	}
	return param;
}		/*end find_params_on_curve*/

#endif /* defined(TWOD) */

/*
*			g_is_correspondence_possible():
*
*	Returns YES if c1 and c2 match in wave_type, match in start and
*	end status at the given nodes, and have the same params on both
*	sides.  Returns NO otherwise.
*/

/*ARGSUSED*/
EXPORT boolean g_is_correspondence_possible(
	HYPER_SURF	*hs1,
	HYPER_SURF	*hs2,
	HYPER_SURF_BDRY **p_hsb,
	HYPER_SURF_BDRY **n_hsb)
{
	if (wave_type(hs1) != wave_type(hs2))
	    return NO;

#if defined(TWOD)
	if (hs1->interface->dim == 2)
	{
	    CURVE     *c1 = Curve_of_hs(hs1);
	    CURVE     *c2 = Curve_of_hs(hs2);
	    NODE      *ns = ((p_hsb != NULL) && (*p_hsb != NULL)) ?
				Node_of_hsb(*p_hsb) : NULL;
	    NODE      *ne = ((n_hsb != NULL) && (*n_hsb != NULL)) ?
				Node_of_hsb(*n_hsb) : NULL;
	    Gas_param *param1 = NULL, *param2 = NULL;

	    if ((ns != NULL) && (start_status(c1) != start_status(c2)))
	    	return NO;
	    if ((ne != NULL) && (end_status(c1) != end_status(c2)))
	    	return NO;

	    param1 = find_params_on_curve(c1,NEGATIVE_SIDE);
	    param2 = find_params_on_curve(c2,NEGATIVE_SIDE);
	    if (param1 != param2)
	    	return NO;

	    param1 = find_params_on_curve(c1,POSITIVE_SIDE);
	    param2 = find_params_on_curve(c2,POSITIVE_SIDE);
	    if (param1 != param2)
	    	return NO;
	}
#endif /* defined(TWOD) */
	return YES;

}		/*end g_is_correspondence_possible*/
