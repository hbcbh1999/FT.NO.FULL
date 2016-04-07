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
*				gmdnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the routines used to handle a bifurcation from a regular
*	to a mach diffraction.
*
*/


#if defined(TWOD)
#include <gdecs/gdecs.h>

/*
*			reg_to_mach_diff_reconfigure():
*
*	This is a simple minded function to deal with the regular to mach
*	diffraction.  We simply untrack the reflected wave, and change the
*	difraction node into a transmission node.  
*/



EXPORT	int reg_to_mach_diff_reconfigure(
	O_CURVE		**newc,
	O_CURVE		**oldc,
	NODE		*oldn)
{
	int		i;

	debug_print("mach_diff","Entering reg_to_mach_diff_reconfigure()\n");
	if (debugging("mach_diff"))
	{
	    (void) printf("WARNING in diffraction_node_propagate(), "
	                  "bifurcation from regular diffraction "
	                  "to mach difraction.\n");
	}

	if (newc[5]->orient == POSITIVE_ORIENTATION)
	{
	    start_status(oldc[5]->curve) = INCIDENT;
	    start_status(newc[5]->curve) = INCIDENT;
	}
	else
	{
	    end_status(oldc[5]->curve) = INCIDENT;
	    end_status(newc[5]->curve) = INCIDENT;
	}
	if (newc[0]->orient == POSITIVE_ORIENTATION)
	{
	    start_status(oldc[0]->curve) = TRANSMITTED;
	    start_status(newc[0]->curve) = TRANSMITTED;
	}
	else
	{
	    end_status(oldc[0]->curve) = TRANSMITTED;
	    end_status(newc[0]->curve) = TRANSMITTED;
	}

		/* untrack any reflected waves */

	for (i = 1; i < 4; i++)
	    if (newc[i]->curve != NULL)
	    	untracked_hyper_surf(newc[i]->curve) = YES;

	node_type(oldn) = TRANSMISSION_NODE;

	debug_print("mach_diff","Leaving reg_to_mach_diff_reconfigure()\n");
	return REPEAT_TIME_STEP_NODE;
}		/*end reg_to_mach_diff_reconfigure*/
#endif /* defined(TWOD) */
