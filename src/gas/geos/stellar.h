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
*				stellar.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if !defined(_STELLAR_H)
#define _STELLAR_H

#include <geos/geosdecs.h>

	/* Structures for accessing equation of state information */

#define		MAX_NUM_ISOTOPS		20

struct _STELLAR_EOS {
	EOS	eos;

	/* Add private EOS parameters below */
	int	_ionmax;		/* number of isotops */
	double	_xmass[MAX_NUM_ISOTOPS]; /* mass fractions of isotops */
	double	_aion[MAX_NUM_ISOTOPS]; /* atomic weight of isotops */
	double	_zion[MAX_NUM_ISOTOPS]; /* atomic number of isotops */

};
typedef struct _STELLAR_EOS STELLAR_EOS;

#define	STELLAR_Eos(state)	((STELLAR_EOS *)Params(state)->eos)

	/* Macros */

#define ionmax(state)    	(STELLAR_Eos(state))->_ionmax
#define xmass(state)    	(STELLAR_Eos(state))->_xmass
#define aion(state)    		(STELLAR_Eos(state))->_aion
#define zion(state)    		(STELLAR_Eos(state))->_zion

#define Cv(state)         ((Gas *) (state))->cv
#define GAMMA(state)      (((Gas *) (state))->gamma - 1.0)
#define Sqrtgm(state)     sqrt(Local_gamma(state))
#define Coef1(state)      (0.5*(Local_gamma(state) + 1.0))
#define Coef2(state)      (0.5*GAMMA(state))
#define Coef3(state)      (0.5*GAMMA(state)/Local_gamma(state))
#define Coef4(state)      (GAMMA(state)/(Local_gamma(state) + 1.0))
#define Coef5(state)      (2.0/(Local_gamma(state) + 1.0))
#define Coef6(state)      (1.0/GAMMA(state))
#define Coef7(state)      (Local_gamma(state)/GAMMA(state))
#define Coef8(state)      (2.0*Local_gamma(state)/(Local_gamma(state) + 1.0))
#define Mu(state)         (sqrt(GAMMA(state)/(Local_gamma(state) + 1.0)))
/*#define R(state)          (Cv(state)*GAMMA(state)) */
#define R(state)         ((Gas *) (state))->R

#endif /* !defined(_STELLAR_H) */
