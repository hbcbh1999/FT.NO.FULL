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
*				gentest.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if !defined(_GENTEST_H)
#define _GENTEST_H

#include <geos/poly.h>

	/* Structures for accessing equation of state information */

struct _GENTEST_EOS {
	POLY_EOS	Poly;
	double		pinf;
};
typedef struct _GENTEST_EOS GENTEST_EOS;

#define	GENTEST_Eos(state)	((GENTEST_EOS *)Params(state)->eos)

	/* Macros */

#define Pinf(state)		(GENTEST_Eos(state)->pinf)
#define	stiff_pressure(state) 	(pressure(state) + Pinf(state))

#endif /* !defined(_GENTEST_H) */
