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
*				geosdecs.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains definitions and structures for equation of 
*	state calculations.
*
*	THIS FILE MAY ONLY BE VALIDLY INCLUDED IN GEOS FILES.  INCLUSION
*	OF THIS FILE IN FILES IN OTHER DIRECTORIES WILL DESTROY THE EQUATION 
*	OF STATE MODULARITY OF THE CODE.
*/

#if !defined(_GEOSDECS_H)
#define _GEOSDECS_H

#include <gdecs/gdecs.h>

		/* possible equations of state */

enum {
	POLYTROPIC	      = FIRST_PHYSICAL_EOS_TYPE +  1,
	POLYTROPIC_ZND,
	STIFFENED_POLYTROPIC,
	BKW,
	LJD,
	HOM_GAS,
	HOM_SOLID,
	SESAME,
	MULTI_COMP_POLYTROPIC,
	JWL,
	MIE_GRUNEISEN,
	ISENTROPIC_TWO_PHASE,
	STELLAR,
	GENTEST
};

#if defined(LOG10)
#define plog(x) log10(x)
#define pexp(x) pow(10.0,x)
#else /* defined(LOG10) */
#define plog(x) log(x)
#define pexp(x) exp(x)
#endif  /* defined(LOG10) */

	/* Geos limited function prototypes */

/* geosutils.c */
IMPORT	void	limit_pressure(double*,double*,int);


#endif /* !defined(_GEOSDECS_H) */
