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
*				poly.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if !defined(_POLY_H)
#define _POLY_H

#include <geos/geosdecs.h>

	/* Structures for accessing equation of state information */

struct _POLY_EOS {
	EOS	eos;

	/* Add private EOS parameters below */

	double gamma; /* polytropic gas constant */
	double R;     /* Ideal gas constant PV = R T */
	double GAMMA; /* Gruneisent exponent = gamma - 1 */
	double cv;    /* Specific heat at constant volume, cv = R/GAMMA */

	/* useful combs of gamma */

	double sqrtgm;	/* sqrt of gamma */

	double mu;	/* mu = sqrt((gamma - 1)/(gamma + 1)) */

			/*                                    2             */
	double coef1;	/* coef1 = (gamma + 1)/2 = 1 / (1 - mu )            */

                        /*                           2          2           */
	double coef2;	/* coef2 = (gamma - 1)/2 = mu  / (1 - mu )          */

			/*                                   2          2   */
	double coef3;	/* coef3 = (gamma - 1)/(2*gamma) = mu  / (1 + mu )  */

        		/*                                     2            */
	double coef4;	/* coef4 = (gamma - 1)/(gamma + 1) = mu             */

        		/*                               2         2        */
	double coef5;	/* coef5 = 2/(gamma + 1) = 1 - mu = (1 + mu )/gamma */

			/*                                2       2         */
	double coef6;	/* coef6 = 1/(gamma - 1) = (1 - mu )/(2*mu )        */

			/*                                    2       2     */
	double coef7;	/* coef7 = gamma/(gamma - 1) = (1 + mu )/(2*mu )    */

        		/*                                       2          */
	double coef8;	/* coef8 = 2 * gamma/(gamma + 1) = 1 + mu           */

	/* The followings are from "Fluid Mechanics" by Landau and Lifshitz,
	 * page 45
	 *
	 * sigma_{ij} = -p delta_{ij} +
	 *             eta (e_{ij}+e_{ji})+(zeta-2/3*eta) (div v) delta_{ij}
	 * e_{ij} = \partial v_i/\partial x_j
	 * eta = dynamic viscosity = ordinary viscosity = shear viscosity
	 * zeta = bulk viscosity = second viscosity = expansion viscosity =
	 *        dilatational viscosity
	 *
	 * Comments from "Physics of Shock Waves and High-Temperature
	 *        Hydrodynamic Phenomena"
	 * by Zeldovich and Raizer, vol 1, page 73
	 *
	 * Sometimes, (zeta-2/3*eta) = lambda is called as second viscosity,
	 * so lambda=-2/3*eta means zeta=0
	 */

        double  shear_visc;    /* dynamic/shear viscosity for fluid */
        double  bulk_visc;    /* bulk viscosity for fluid */
        double  heat_coeff;    /* N-S thermal conduction coefficient */
};
typedef struct _POLY_EOS POLY_EOS;

#define	POLY_Eos(state)	((POLY_EOS *)Params(state)->eos)

	/* Macros */

#define Gamma(state)	  (POLY_Eos(state)->gamma)
#define R(state)	  (POLY_Eos(state)->R)
#define GAMMA(state)	  (POLY_Eos(state)->GAMMA)
#define Cv(state)	  (POLY_Eos(state)->cv)
#define	Sqrtgm(state) 	  (POLY_Eos(state)->sqrtgm)
#define	Mu(state) 	  (POLY_Eos(state)->mu)
#define	Coef1(state) 	  (POLY_Eos(state)->coef1)
#define	Coef2(state) 	  (POLY_Eos(state)->coef2)
#define	Coef3(state) 	  (POLY_Eos(state)->coef3)
#define	Coef4(state) 	  (POLY_Eos(state)->coef4)
#define	Coef5(state) 	  (POLY_Eos(state)->coef5)
#define	Coef6(state) 	  (POLY_Eos(state)->coef6)
#define	Coef7(state) 	  (POLY_Eos(state)->coef7)
#define	Coef8(state) 	  (POLY_Eos(state)->coef8)
#define	Shear_visc(state) (POLY_Eos(state)->shear_visc)
#define	Bulk_visc(state)  (POLY_Eos(state)->bulk_visc)
#define Heat_coeff(state)	(POLY_Eos(state)->heat_coeff)

	/* Polytropic specific functions */
IMPORT	void	set_POLY_coefs(Gas_param*);

#endif /* !defined(_POLY_H) */
