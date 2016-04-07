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
*				s2phase.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*
*/

#if !defined(_S2PHASE_H)
#define _S2PHASE_H

#include <geos/geosdecs.h>

struct _Spline {
    double *x; /* Independent variable array, monotone nondecreasing */
    double *y;   /* Dependent variable array */
    int    m;    /* Number of data points */

    int    k;    /* Degree of spline polynomial */
    int    n;    /* Number of knots in the spline approximation */
    double *t;   /* Knot array */
    double *c;   /* B-spline coefficients */
    int   nest; /* Allocated length of the arrays t and c */

    double     atol;
    double     rtol;
    double     terr;
};
typedef struct _Spline Spline;


	/* Structures for accessing equation of state information */

struct _S2PHASE_EOS {
  EOS    eos;
  double	 a_sat_l;     /* sound speed squared of saturated liquid */ 
  double	 a_sat_v;     /* sound speed squared of saturated vapor */ 
  double	 p_sat_l;     /* pressure of saturated liquid */
  double  p_inf;       /* pressure P_infinity for liquid */ 
  double	 p_sat_v;     /* pressure of saturated vapor */
  double  p_vl;        /* mixed domain coefficient */
  double	 rho_sat_l;   /* density of saturated liquid */
  double	 rho_sat_v;   /* density of saturated vapor */
  double  e_sat_l;     /* specific internal energy of saturated liquid */
  double  e_sat_v;     /* specific internal energy of saturated vapor */
  double	 delta_e;     /* vapor - liquid energy jump (related to the heat of vaporization)*/
  double  e_inf;       /* energy E_infinity for liquid */ 
  double  S_0;         /* entropy (remains constant) */
  double  R_l;         /* gas constant of saturated liquid */
  double  R_v;         /* gas constant of saturated vapor */
  double	 t_sat_v;     /* temperature of saturated vapor */
  double	 t_sat_l;     /* temperature of saturated liquid */
  double  cv_l;        /* C_V of saturated liquid */
  double  gamma_l;     /* adiabatic exponent of liquid */
  double  gamma_v;     /* adiabatic exponent of vapor */
  double  eta_l;       /* eos coef. for liquid */
  double  eta_v;       /* eos coef. for vapor */
  
  /* The followings are from "Fluid Mechanics" by Landau and Lifshitz, page 45
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

  double  shear_visc_l;    /* dynamic viscosity for liquid */
  double  shear_visc_v;    /* dynamic viscosity for vapor */
  double  bulk_visc_l;    /* bulk viscosity for liquid */
  double  bulk_visc_v;    /* bulk viscosity for vapor */

  Spline *int_c_drho_over_rho_mix_spline;
};
typedef struct _S2PHASE_EOS S2PHASE_EOS;

#define	S2PHASE_Eos(state)	((S2PHASE_EOS *)Params(state)->eos)

	/* Macros */

#define A_sat_l(state)		(S2PHASE_Eos(state)->a_sat_l)
#define A_sat_v(state)		(S2PHASE_Eos(state)->a_sat_v)
#define P_sat_l(state)		(S2PHASE_Eos(state)->p_sat_l)
#define P_inf(state)		(S2PHASE_Eos(state)->p_inf)
#define P_sat_v(state)		(S2PHASE_Eos(state)->p_sat_v)
#define P_vl(state)		(S2PHASE_Eos(state)->p_vl)
#define Rho_sat_l(state)	(S2PHASE_Eos(state)->rho_sat_l)
#define Rho_sat_v(state)	(S2PHASE_Eos(state)->rho_sat_v)
#define Delta_e(state)	        (S2PHASE_Eos(state)->delta_e)
#define E_inf(state)	        (S2PHASE_Eos(state)->e_inf)
#define E_sat_l(state)	        (S2PHASE_Eos(state)->e_sat_l)
#define E_sat_v(state)	        (S2PHASE_Eos(state)->e_sat_v)
#define S_0(state)	        (S2PHASE_Eos(state)->S_0)
#define R_l(state)            	(S2PHASE_Eos(state)->R_l)
#define R_v(state)            	(S2PHASE_Eos(state)->R_v)
#define Gamma_l(state)	        (S2PHASE_Eos(state)->gamma_l)
#define Gamma_v(state)	        (S2PHASE_Eos(state)->gamma_v)
#define Eta_l(state)	        (S2PHASE_Eos(state)->eta_l)
#define Eta_v(state)	        (S2PHASE_Eos(state)->eta_v)
#define T_sat_l(state)	        (S2PHASE_Eos(state)->t_sat_l)
#define T_sat_v(state)	        (S2PHASE_Eos(state)->t_sat_v)

#define Shear_visc_l(state)	        (S2PHASE_Eos(state)->shear_visc_l)
#define Shear_visc_v(state)	        (S2PHASE_Eos(state)->shear_visc_v)
#define Bulk_visc_l(state)	        (S2PHASE_Eos(state)->bulk_visc_l)
#define Bulk_visc_v(state)	        (S2PHASE_Eos(state)->bulk_visc_v)
#define cv_v(state)           		(S2PHASE_Eos(state)->cv_v)

#define Int_c_drho_over_rho_mix_spline(state)				\
	(S2PHASE_Eos(state)->int_c_drho_over_rho_mix_spline)

#endif /* !defined(_S2PHASE_H) */
