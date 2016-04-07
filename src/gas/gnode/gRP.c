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
*				gRP.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*
*	Contains functions for the allocation and maintenance of the
*	RP_DATA structure.
*/


#if defined(TWOD)
#include <gdecs/gdecs.h>

EXPORT	RP_DATA *allocate_RP_DATA_structure(
	size_t		sizest,
	boolean		dynamic,
	int		stype)
{
	RP_DATA		*RP;
	int		i;

	if (sizest == 0)
	    return NULL;
	if (dynamic == YES)
	{
	    RP = (RP_DATA *) store(sizeof(RP_DATA));
	    if (RP == NULL)
		return NULL;

	    for (i = 0; i < MAX_N_CURVES; i++)
	    {
	        if ((RP->state[i] = g_alloc_intfc_state(sizest)) == NULL)
	    	    return NULL;
	    }
	}
	else
	{
	    scalar(&RP,sizeof(RP_DATA));
	    if (RP == NULL)
		return NULL;

	    for (i = 0; i < MAX_N_CURVES; i++)
	    {
	    	g_alloc_state(&RP->state[i],sizest);
	    	if (RP->state[i]  == NULL)
	    	    return NULL;
	    }
	}
	RP->intfc_table_storage = dynamic;
	RP->stype = stype;
	RP->ang_dir = ANGLE_DIRECTION_NOT_SET;

	return RP;
}		/*end allocate_RP_DATA_structure*/

EXPORT	void	free_RP_DATA_structure(
	RP_DATA	*RP)
{
	int	i;
	if (RP->intfc_table_storage == YES)
		return;

	for (i = 0; i < MAX_N_CURVES; i++)
	{
		if (RP->state[i] != NULL)
			free(RP->state[i]);
	}
	free(RP);
}		/*end free_RP_DATA_structure*/

EXPORT	void copy_RP_DATA_structure(
	RP_DATA		*cRP,
	RP_DATA		*RP)
{
	int		i;

	if (cRP == NULL) return;
	if (RP == NULL) /* Clear data in cRP */
	{
		cRP->ang_dir = ANGLE_DIRECTION_NOT_SET;
		cRP->stype = UNKNOWN_STATE;
		for (i = 0; i < MAX_N_CURVES; i++)
		{
			cRP->ang[i] = 0.0;
			g_clear_state(cRP->state[i],g_sizest());
			cRP->M[i] = 0.0;
			cRP->theta[i] = 0.0;
		}
		return;
	}

	cRP->ang_dir = RP->ang_dir;
	cRP->stype = RP->stype;
	for (i = 0; i < MAX_N_CURVES; i++)
	{
		cRP->ang[i] = RP->ang[i];
		set_state(cRP->state[i],cRP->stype,RP->state[i]);
		cRP->M[i] = RP->M[i];
		cRP->theta[i] = RP->theta[i];
	}
}		/*end copy_RP_DATA_structure*/
#endif /* defined(TWOD) */
