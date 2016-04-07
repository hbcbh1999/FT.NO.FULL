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
*				jwl.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if !defined(_JWL_H)
#define _JWL_H

#include <geos/geosdecs.h>

	/* Structures for accessing equation of state information */

struct _JWL_EOS {
	EOS	eos;

	/* Basic independent parameters defining material */
	double	_Rho0;	/* Reference density */
	double	_P1, _P2;	/* Units of pressure */
	double	_R1, _R2;	/* dimensionless */
	double	_W;	/* Gruneisen coefficient, dimensionless */
	double	_R;	/* Gas constant, defines temperature units */
	double	_dH;	/* heat of detonation */

	/* Useful combinations of EOS parameters */

	double	_Rho1;	/* Rho0*R1 */
	double	_Rho2;	/* Rho0*R2 */
	double	_W_1;	/* W/r1 */
	double	_W_2;	/* W/r2 */
	double	_Wp1;	/* W + 1*/
	double	_Pref0;	/* p1*exp(-R1) + p2*exp(-R2)*/
};
typedef struct _JWL_EOS JWL_EOS;

#define	JWL_Eos(state)	((JWL_EOS *)Params(state)->eos)
#define	Rho0(state)	(JWL_Eos(state))->_Rho0
#define	P1(state)	(JWL_Eos(state))->_P1
#define	P2(state)	(JWL_Eos(state))->_P2
#define	R1(state)	(JWL_Eos(state))->_R1
#define	R2(state)	(JWL_Eos(state))->_R2
#define	W(state)	(JWL_Eos(state))->_W
#define	Rho1(state)	(JWL_Eos(state))->_Rho1
#define	Rho2(state)	(JWL_Eos(state))->_Rho2
#define	W1(state)	(JWL_Eos(state))->_W_1
#define	W2(state)	(JWL_Eos(state))->_W_2
#define	Wp1(state)	(JWL_Eos(state))->_Wp1
#define	R(state)	(JWL_Eos(state))->_R
#define	Pref0(state)	(JWL_Eos(state))->_Pref0
#define	dH(state)	(JWL_Eos(state))->_dH

#endif /* !defined(_JWL_H) */
