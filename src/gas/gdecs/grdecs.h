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
*				grdecs.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GRDECS_H)
#define _GRDECS_H


		/* adjustable parameters */
#define EPS		0.000001 /* TOLERANCE */
				 /* Governs convergence of iterations in */
				 /* Riemann problem and is used to test if */
				 /* certain fabs(a - b) < EPS */

		/* Sonic transition tolerances */
#define	SONIC_TOL 	0.01	/* TOLERANCE */
#define SONIC_MINUS	0.99	/* TOLERANCE 1.0 - SONIC_TOL */
#define SONIC_MINUS_SQR	0.0001	/* TOLERANCE (1.0 - SONIC_TOL)^2 */
#define	SONIC_PLUS	1.01	/* TOLERANCE 1.0 + SONIC_TOL */
#define	SONIC_PLUS_SQR	1.0201	/* TOLERANCE (.10 + SONIC_TOL)^2 */

		/* Types of Riemann solvers */

enum _RIEMANN_SOLVER_TYPE {
	EXACT_RIEMANN_SOLVER	     = 1,
	LINEAR_APPROX_RIEMANN_SOLVER = 2,
	CG_APPROX_RIEMANN_SOLVER     = 3  /*Colella & Glaz's Riemann Solver*/
};
typedef enum _RIEMANN_SOLVER_TYPE RIEMANN_SOLVER_TYPE;

		/* Solution methods for exact Riemann Solver */

enum _RIEMANN_SOLVER_METHOD {
	GODUNOV_RS = 1,
	SECANT_RS  = 2
};
typedef enum _RIEMANN_SOLVER_METHOD RIEMANN_SOLVER_METHOD;

		/* possible types of coordinate systems */

enum	_GEOMETRY	{RECTANGULAR, CYLINDRICAL, SPHERICAL};
typedef	enum _GEOMETRY GEOMETRY;

		/* wave family values */

enum _WAVE_FAMILY {
	LEFT_FAMILY  = 1,
	RIGHT_FAMILY = 2
};
typedef enum _WAVE_FAMILY WAVE_FAMILY;

		/* parameter choices for s_polar_4() */

enum _SHOCK_STRENGTH_PARAMETER {
	BEHIND_PRESSURE	= 1,
	BEHIND_VELOCITY,
	SHOCK_SPEED,
	SHOCK_MACH_NUMBER
};
typedef enum _SHOCK_STRENGTH_PARAMETER SHOCK_STRENGTH_PARAMETER;

		/* choice of roots for the intersection of two shock polars */

enum {
	WEAK = 1,
	STRONG
};

#if defined(FULL_PHYSICS)
		/* Values returned by is_regular_diffraction_node() */

enum _DIFFRACTION_STATUS {
	ERROR_DIFFRACTION		      = 0,
	REGULAR_DIFFRACTION,
	ANOMALOUS_REFLECTION,
	REGULAR_TO_MACH_DIFFRACTION,
	PRECURSOR_WITH_REFLECTED_RAREFACTION,
	PRECURSOR_WITH_REFLECTED_SHOCK
};
typedef enum _DIFFRACTION_STATUS DIFFRACTION_STATUS;

		/* Values returned by find_transmission_node_states() */

enum _TRANSMISSION_STATUS {
	ERROR_TRANSMISSION = 0,
	REGULAR_TRANSMISSION,
	BIFURCATION_TRANSMISSION
};
typedef enum _TRANSMISSION_STATUS TRANSMISSION_STATUS;

#endif /* defined(FULL_PHYSICS) */


#endif /* !defined(_GRDECS_H) */
