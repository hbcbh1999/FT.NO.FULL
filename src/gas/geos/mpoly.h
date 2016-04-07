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
*				mpoly.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if defined(MULTI_COMPONENT)
#if !defined(_MPOLY_H)
#define _MPOLY_H

#include <geos/geosdecs.h>

	/* Structures for accessing equation of state information */

struct _MPOLY_EOS {
	EOS	eos;

	/* Add private EOS parameters below */

	double _gamma[MAX_NUM_GAS_COMPS];/* Molecular weights of components */
	double _M[MAX_NUM_GAS_COMPS];	/* Molecular weights of components */
	double R;			/* Ideal gas constant PV = R T */
};
typedef struct _MPOLY_EOS MPOLY_EOS;

#define	MPOLY_Eos(state)	((MPOLY_EOS *)Params(state)->eos)

	/* Macros */
#define	gamma(state)	(MPOLY_Eos(state))->_gamma
#define	M(state)	(MPOLY_Eos(state))->_M

#endif /* !defined(_MPOLY_H) */
#endif /* defined(MULTI_COMPONENT) */
