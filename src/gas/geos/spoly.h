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
*				spoly.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if !defined(_SPOLY_H)
#define _SPOLY_H

#include <geos/poly.h>

	/* Structures for accessing equation of state information */

struct _SPOLY_EOS {
	POLY_EOS Poly;
	double	 pinf;
	double	 einf;
	double    rhoinf;
	double	 et;
};
typedef struct _SPOLY_EOS SPOLY_EOS;

#define	SPOLY_Eos(state)	((SPOLY_EOS *)Params(state)->eos)

	/* Macros */

#define Pinf(state)		(SPOLY_Eos(state)->pinf)
#define Einf(state)		(SPOLY_Eos(state)->einf)
#define Et(state)		(SPOLY_Eos(state)->et)
#define	stiff_pressure(state) 	(pressure(state) + Pinf(state))

#endif /* !defined(_SPOLY_H) */
