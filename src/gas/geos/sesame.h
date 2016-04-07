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
*				sesame.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains definitions and structures for the SESAME
*	equation of state.
*/

#if !defined(_SESAME_H)
#define _SESAME_H

#include <geos/geosdecs.h>

	/* Structures for accessing equation of state information */

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
    typedef void (SGET)(int*,int*,double*,int*,int*,int*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

struct _SESAME_EOS {
	EOS	     eos;

	/* Add private EOS parameters below */
	double	     _rmin, _rmax, _rref;
	double	     _tmin, _tmax, _tref;
	double	     _emin, _emax, _eref;
	double	     _pmin, _pmax, _pref;
	double	     _smin, _smax;
	double	     _cref;
	double	     pmin;
	double	     _rs_entropy_scale, _rs_entropy_shift;
	double	     _ps_entropy_scale, _ps_entropy_shift;
	double	     _gru_gam_abs_max;
	double	     eps; /* Measures discreteness of sesame table */
	CHART	     **root;
	Front	     **fr;
	Wave	     **wave;
	int	     num_ses_tabs;
	SGET	     *_sget;
	void	     (*_set_RT_entropy)(Front*,Wave*,struct _SESAME_EOS*);

		/* Information for sesame table lookup */

	struct _SESTAB {
	    struct _S2DIR {
	    	int	lcmx; /* Dimension of tbls array */
	    	int	nrs;  /* number or regions in problem */
	    	int	lcfw[10]; /*array used for directory to tbls*/
	    } *s2dir;
	    double *tbls; /* Original SESAME tables */
	} sestab;

	struct _DEPARAMS {
	    double rmid, tmid;
	    double rlgmid, tlgmid;
	    double eta, xi;
	    double temperature, rho;
#if defined(PHASE_CODE)
	    double dlogr, dlogt;
	    double dxi, deta;
	    double rlgmin, tlgmin;
	    double t2min, t2max;
	    double r2min, r2max;
	    double *R, *T, *E, *P, *slopeE, *slopeP;
	    int	  len;
	    double *pb_R, *pb_T, *pb_E, *pb_P;
	    double *pb_slopeR, *pb_slopeE, *pb_slopeP;
	    int   pb_len;
#endif /* defined(PHASE_CODE) */
	} de_params;

	double      abser0; /* Absolute truncation error for computing entropy*/
	double      reler0; /* Relative truncation error for computing entropy*/
	double      terr0; /* density/temp tolerance for computing entropy*/

		/* Sesame data file name and material id */

	char       seslib[120]; /* Sesame data file name */
	int        ids2; /* material id */

		/* File name with restart Sesame data */
	char	      ses_restart_file_name[120];
	const IO_TYPE *restart_io_type;

		/* Arrays used in intializing the sesame  */
		/* hyp solution functions.  This storage   */
		/* is freed and these pointers are set to */
		/* null after initialization.		  */

#if defined(PHASE_CODE)
	double	   *S_on_ref_curve;	/* Entropy along reference curve  */
	int	   n_pts_ref_curve;	/*number of points on the ref curve*/
#endif /* defined(PHASE_CODE) */
	double	   *slpe; /* Spline slopes for energy on isotherms */
	double	   *slpp; /* Spline slopes for pressure on isotherms */
	double	   Smax; /* Maximum entropy on vapor dome; used in RP solver */
	double	   BIS_EPS;	/*TOLERANCE*/
	double	   ABS_SES_EPS;	/*TOLERANCE*/
	const char **_ses_table_name;
};
typedef struct _SESAME_EOS SESAME_EOS;
#if defined(__cplusplus)
typedef struct _SESAME_EOS::_SESTAB SESTAB;
typedef struct _SESAME_EOS::_SESTAB::_S2DIR S2DIR;
typedef struct _SESAME_EOS::_DEPARAMS DEPARAMS;
#else /* defined(__cplusplus) */
typedef struct _SESTAB SESTAB;
typedef struct _S2DIR S2DIR;
typedef struct _DEPARAMS DEPARAMS;
#endif /* defined(__cplusplus) */

struct _Ses_Front {
	Front		front;
	SESAME_EOS	*seos;
};
typedef struct _Ses_Front Ses_Front;
#define Ses_front_seos(front)	((Ses_Front*)(front))->seos

struct _Ses_Wave {
	Wave		wave;
	SESAME_EOS	*seos;
};
typedef struct _Ses_Wave Ses_Wave;
#define Ses_wave_seos(wave)	((Ses_Wave*)(wave))->seos

#define SESAME_Eos(state) ((SESAME_EOS *)Params(state)->eos)
#define Nrho_hyp(seos)	  ((seos)->fr[SESAME_RHO_TEMP]->rect_grid->gmax[0])
#define Ntemp_hyp(seos)	  ((seos)->fr[SESAME_RHO_TEMP]->rect_grid->gmax[1])
#define Rho_min(seos)	  ((seos)->_rmin)
#define Rho_max(seos)	  ((seos)->_rmax)
#define Rho_ref(seos)	  ((seos)->_rref)

#define Temp_min(seos)    ((seos)->_tmin)
#define Temp_max(seos)	  ((seos)->_tmax)
#define Temp_ref(seos)	  ((seos)->_tref)

#define Reference_sound_speed(seos)	((seos)->_cref)

#define Pressure_min(seos)	((seos)->_pmax)
#define Pressure_max(seos)	((seos)->_pmin)
#define Pressure_ref(seos)	((seos)->_pref)

#define Energy_min(seos)	((seos)->_emax)
#define Energy_max(seos)	((seos)->_emin)
#define Energy_ref(seos)	((seos)->_eref)

#define Entropy_min(seos)	((seos)->_smin)
#define Entropy_max(seos)	((seos)->_smax)

#define Gru_gam_abs_max(seos)	((seos)->_gru_gam_abs_max)

#define RS_entropy_scale(seos)	((seos)->_rs_entropy_scale)
#define RS_entropy_shift(seos)	((seos)->_rs_entropy_shift)
#define PS_entropy_scale(seos)	((seos)->_ps_entropy_scale)
#define PS_entropy_shift(seos)	((seos)->_ps_entropy_shift)
#define GetSesameTable(seos,ir,ids2,tbls,lcnt,lu,ifl)			\
    (*(seos)->_sget)(ir,ids2,tbls,lcnt,lu,ifl)
#define ses_table_name(seos,i)	     ((seos)->_ses_table_name)[i]
#define set_RT_entropy(fr,wave,seos) (*(seos)->_set_RT_entropy)(fr,wave,seos)

	/* Declarations used by the SESAME EOS */

	/*   ordering of groups of tables      */

enum _SESAME_TABLE_TYPE {
	SESAME_RHO_TEMP=0,
	SESAME_RHO_ENERGY,
	SESAME_RHO_ENTROPY,
	SESAME_PRESS_ENTROPY,
	SESAME_VOLUME_PRESSURE,
	NUMBER_SESAME_TABLES
};
typedef enum _SESAME_TABLE_TYPE SESAME_TABLE_TYPE;

#define Ses_params(state)	((SESAME_EOS *)Params(state)->eos)

#define Sestab(state)		(Ses_params(state)->sestab)

#define TBLS(state)		(Sestab(state).tbls)

#define LCMX(state)		(Sestab(state).s2dir->lcmx)

#define NRS(state)		(Sestab(state).s2dir->nrs)

#define LCFW(state)		(Sestab(state).s2dir->lcfw)

#define Eps(state)		(Ses_params(state)->eps)

#define Rmin(state)		(rmin(Ses_params(state)))

#define Rmax(state)		(rmax(Ses_params(state)))

#define Tmin(state)		(tmin(Ses_params(state)))

#define Tmax(state)		(tmax(Ses_params(state)))

#define Pmin(state)		(Ses_params(state)->pmin)

enum {
	COMP_PURE_PHASE = MIN_INTERIOR_COMP,
	COMP_MIXED_PHASE
};

	/* Flags for the return of data from tri solution */

enum {
	EVALUATE_DENSITY     =	0x01,
	EVALUATE_PRESSURE    =	0x02,
	EVALUATE_ENERGY	     =	0x04,
	EVALUATE_ENTROPY     =	0x08,
	EVALUATE_TEMPERATURE =	0x10,
	EVALUATE_ADB_GAMMA   =	0x20,
	EVALUATE_GRU_GAMMA   =	0x40,
	EVALUATE_IDPOCR	     =	0x80
};

#define evaluate_density(flag)		((flag) & EVALUATE_DENSITY)
#define evaluate_pressure(flag)		((flag) & EVALUATE_PRESSURE)
#define evaluate_energy(flag)		((flag) & EVALUATE_ENERGY)
#define evaluate_entropy(flag)		((flag) & EVALUATE_ENTROPY)
#define evaluate_temperature(flag)	((flag) & EVALUATE_TEMPERATURE)
#define evaluate_adb_gamma(flag)	((flag) & EVALUATE_ADB_GAMMA)
#define evaluate_gru_gamma(flag)	((flag) & EVALUATE_GRU_GAMMA)
#define evaluate_idpocr(flag)		((flag) & EVALUATE_IDPOCR)


#if defined(OLD_HUGONIOTS)

struct _PR_HG_SESAME_PARAMS {
	double p_avg, p1, p1_grid, tau0;
	Locstate state;
};
typedef struct _PR_HG_SESAME_PARAMS PR_HG_SESAME_PARAMS;

struct _DENS_HG_SESAME_PARAMS {
	double htd, rho1, rho1_grid, p0;
	Locstate state;
};
typedef struct _DENS_HG_SESAME_PARAMS DENS_HG_SESAME_PARAMS;

#else /* defined(OLD_HUGONIOTS) */

struct _PR_HG_SESAME_PARAMS {
	double p_avg_inv, e0, tau0;
	Locstate state;
};
typedef struct _PR_HG_SESAME_PARAMS PR_HG_SESAME_PARAMS;

struct _DENS_HG_SESAME_PARAMS {
	double rho1, rho1_grid, p0, e0;
	Locstate state;
};
typedef struct _DENS_HG_SESAME_PARAMS DENS_HG_SESAME_PARAMS;

#endif /* defined(OLD_HUGONIOTS) */

	/* State structure and macros for density - temperature table */

struct _SES_RT_STATE {
	double cold_p;
	double cold_e;
	double reduced_p;
	double reduced_e;
	double S;
	double adb_gam;
	double gru_gam;
	double riv;
};
typedef struct _SES_RT_STATE SES_RT_STATE;

#define ses_rt_state(state)	((SES_RT_STATE *) state)
#define ses_rt_coldp(state)	ses_rt_state(state)->cold_p
#define ses_rt_colde(state)	ses_rt_state(state)->cold_e
#define ses_rt_redp(state)	ses_rt_state(state)->reduced_p
#define ses_rt_rede(state)	ses_rt_state(state)->reduced_e
#define ses_rt_S(state)		ses_rt_state(state)->S
#define ses_rt_adb_gam(state)	ses_rt_state(state)->adb_gam
#define ses_rt_gru_gam(state)	ses_rt_state(state)->gru_gam
#define ses_rt_riv(state)	ses_rt_state(state)->riv

enum {RT_CP=0,RT_CE,RT_RP,RT_RE,RT_S,RT_AG,RT_GG,RT_RF};

#define SESAME_LOGS

#if defined(SESAME_LOGS)
#   define ses_rt_grid_from_rho(rho,seos)	log((rho)/Rho_ref(seos))
#   define ses_rt_rho_from_grid(rho_grid,seos)	(Rho_ref(seos)*exp(rho_grid))
#   define dses_rt_rho_from_grid(rho,seos)	(rho)
#   define dlog_rho_drho_grid(rho,seos)	(1.0)
#   define ses_rt_grid_from_temp(T,seos)        log1p((T)/Temp_ref(seos))
#   define ses_rt_temp_from_grid(T_grid,seos)	(Temp_ref(seos)*expm1(T_grid))
#   define dses_rt_temp_from_grid(T,seos)	((T) + Temp_ref(seos))
#   define dT_grid_dlogT(T,seos)		((T)/((T) + Temp_ref(seos)))
#else /* defined(SESAME_LOGS) */
#   define ses_rt_grid_from_rho(rho,seos)	((rho)/Rho_ref(seos))
#   define ses_rt_rho_from_grid(rho_grid,seos)	(Rho_ref(seos)*(rho_grid))
#   define dses_rt_rho_from_grid(rho,seos)	Rho_ref(seos)
#   define dlog_rho_drho_grid(rho,seos)	(Rho_ref(seos)/(rho))
#   define ses_rt_grid_from_temp(T,seos)	((T)/Temp_ref(seos))
#   define ses_rt_temp_from_grid(T_grid,seos)	(Temp_ref(seos)*(T_grid))
#   define dses_rt_temp_from_grid(T,seos)	Temp_ref(seos)
#   define dT_grid_dlogT(T,seos)		(Temp_ref(seos)/(T))
#endif /* defined(SESAME_LOGS) */

	/* State structure and macros for density - energy table */

struct _SES_RE_STATE {
	double cold_p;
	double reduced_p;
	double T_var;
	double S;
	double adb_gam;
	double gru_gam;
	double riv;
};
typedef struct _SES_RE_STATE SES_RE_STATE;

#define ses_re_state(state)	((SES_RE_STATE *) state)
#define ses_re_coldp(state)	ses_re_state(state)->cold_p
#define ses_re_redp(state)	ses_re_state(state)->reduced_p
#define ses_re_Tvar(state)	ses_re_state(state)->T_var
#define ses_re_S(state)		ses_re_state(state)->S
#define ses_re_adb_gam(state)	ses_re_state(state)->adb_gam
#define ses_re_gru_gam(state)	ses_re_state(state)->gru_gam
#define ses_re_riv(state)	ses_re_state(state)->riv

enum {RE_CP=0,RE_RP,RE_T,RE_S,RE_AG,RE_GG,RE_RF};

#if defined(SESAME_LOGS)
#   define ses_re_grid_from_engy(e,seos)	log1p((e)/Energy_ref(seos))
#   define ses_re_grid_from_rho(rho,seos)	log((rho)/Rho_ref(seos))
#   define ses_re_rho_from_grid(x,seos)	(Rho_ref(seos)*exp(x))
#   define ses_re_engy_from_grid(x,seos)	(Energy_ref(seos)*expm1(x))
#else /* defined(SESAME_LOGS) */
#   define ses_re_grid_from_rho(x,seos)	((x)/Rho_ref(seos))
#   define ses_re_rho_from_grid(x,seos)	(Rho_ref(seos)*(x))
#   define ses_re_grid_from_engy(x,seos)	((x)/Energy_ref(seos))
#   define ses_re_engy_from_grid(x,seos)	(Energy_ref(seos)*(x))
#endif /* defined(SESAME_LOGS) */
#   define ses_re_temp_from_var(x,seos)	ses_rt_temp_from_grid(x,seos)

	/* State structure and macros for density - entropy table */

struct _SES_RS_STATE {
	double cold_p;
	double cold_e;
	double reduced_p;
	double reduced_e;
	double T_var;
	double adb_gam;
	double gru_gam;
	double riv;
};
typedef struct _SES_RS_STATE SES_RS_STATE;

#define ses_rs_state(state)	((SES_RS_STATE *) state)
#define ses_rs_coldp(state)	ses_rs_state(state)->cold_p
#define ses_rs_colde(state)	ses_rs_state(state)->cold_e
#define ses_rs_redp(state)	ses_rs_state(state)->reduced_p
#define ses_rs_rede(state)	ses_rs_state(state)->reduced_e
#define ses_rs_Tvar(state)	ses_rs_state(state)->T_var
#define ses_rs_adb_gam(state)	ses_rs_state(state)->adb_gam
#define ses_rs_gru_gam(state)	ses_rs_state(state)->gru_gam
#define ses_rs_riv(state)	ses_rs_state(state)->riv

enum {RS_CP=0,RS_CE,RS_RP,RS_RE,RS_T,RS_AG,RS_GG,RS_RF,NUM_SES_VAR};

#if defined(SESAME_LOGS)
#   define ses_rs_grid_from_rho(x,seos)	log((x)/Rho_ref(seos))
#   define ses_rs_rho_from_grid(x,seos)	(Rho_ref(seos)*exp(x))
#   define ses_rs_grid_from_entpy(x,seos)			\
			(RS_entropy_scale(seos)*(x) + RS_entropy_shift(seos))
#   define ses_rs_entpy_from_grid(x,seos)			\
			(((x) - RS_entropy_shift(seos))/RS_entropy_scale(seos))
#else /* defined(SESAME_LOGS) */
#   define ses_rs_grid_from_rho(x,seos)	((x)/Rho_ref(seos))
#   define ses_rs_rho_from_grid(x,seos)	(Rho_ref(seos)*(x))
#   define ses_rs_grid_from_entpy(x,seos)	(x)
#   define ses_rs_entpy_from_grid(x,seos)	(x)
#endif /* defined(SESAME_LOGS) */
#   define ses_rs_temp_from_var(x,seos)	ses_rt_temp_from_grid(x,seos)

	/* State structure and macros for pressure - entropy table */

struct _SES_PS_STATE {
	double cold_e;
	double reduced_e;
	double T_var;
	double rho_var;
	double adb_gam;
	double gru_gam;
	double riv;
};
typedef struct _SES_PS_STATE SES_PS_STATE;

#define ses_ps_state(state)	((SES_PS_STATE *) state)
#define ses_ps_colde(state)	ses_ps_state(state)->cold_e
#define ses_ps_rede(state)	ses_ps_state(state)->reduced_e
#define ses_ps_Tvar(state)	ses_ps_state(state)->T_var
#define ses_ps_rho_var(state)	ses_ps_state(state)->rho_var
#define ses_ps_adb_gam(state)	ses_ps_state(state)->adb_gam
#define ses_ps_gru_gam(state)	ses_ps_state(state)->gru_gam
#define ses_ps_riv(state)	ses_ps_state(state)->riv

enum {PS_CE=0,PS_RE,PS_T,PS_RHO,PS_AG,PS_GG,PS_RF};

#if defined(SESAME_LOGS)
#   define ses_ps_grid_from_press(x,seos)	log1p((x)/Pressure_ref(seos))
#   define ses_ps_press_from_grid(x,seos)	(Pressure_ref(seos)*expm1(x))
#   define ses_ps_grid_from_entpy(x,seos)				\
			(PS_entropy_scale(seos)*(x) + PS_entropy_shift(seos))
#   define ses_ps_entpy_from_grid(x,seos)				\
			(((x) - PS_entropy_shift(seos))/PS_entropy_scale(seos))
#else /* defined(SESAME_LOGS) */
#   define ses_ps_grid_from_press(x,seos)	((x)/Pressure_ref(seos))
#   define ses_ps_press_from_grid(x,seos)	(Pressure_ref(seos)*(x))
#   define ses_ps_grid_from_entpy(x,seos)	(x)
#   define ses_ps_entpy_from_grid(x,seos)	(x)
#endif /* defined(SESAME_LOGS) */
#   define ses_ps_temp_from_var(x,seos)	ses_rt_temp_from_grid(x,seos)
#   define ses_ps_rho_from_var(x,seos)	ses_rt_rho_from_grid(x,seos)

	/* State structure and macros for volume - pressure table */

struct _SES_VP_STATE {
	double cold_e;
	double reduced_e;
	double T_var;
	double S;
	double adb_gam;
	double gru_gam;
	double riv;
};
typedef struct _SES_VP_STATE SES_VP_STATE;

#define ses_vp_state(state)	((SES_VP_STATE *) state)
#define ses_vp_colde(state)	ses_vp_state(state)->cold_e
#define ses_vp_rede(state)	ses_vp_state(state)->reduced_e
#define ses_vp_Tvar(state)	ses_vp_state(state)->T_var
#define ses_vp_S(state)		ses_vp_state(state)->S
#define ses_vp_adb_gam(state)	ses_vp_state(state)->adb_gam
#define ses_vp_gru_gam(state)	ses_vp_state(state)->gru_gam
#define ses_vp_riv(state)	ses_vp_state(state)->riv

enum {VP_CE=0,VP_RE,VP_T,VP_S,VP_AG,VP_GG,VP_RF};

#if defined(SESAME_LOGS)
#   define ses_vp_grid_from_vol(x,seos)	  log(Rho_ref(seos)*(x))
#   define ses_vp_vol_from_grid(x,seos)	  (exp(x)/Rho_ref(seos))
#   define ses_vp_grid_from_press(x,seos) log1p((x)/Pressure_ref(seos))
#   define ses_vp_press_from_grid(x,seos) (Pressure_ref(seos)*expm1(x))
#else /* defined(SESAME_LOGS) */
#   define ses_vp_grid_from_vol(x,seos)	  (Rho_ref(seos)*(x))
#   define ses_vp_vol_from_grid(x,seos)	  ((x)/Rho_ref(seos))
#   define ses_vp_grid_from_press(x,seos) ((x)/Pressure_ref(seos))
#   define ses_vp_press_from_grid(x,seos) (Pressure_ref(seos)*(x))
#endif /* defined(SESAME_LOGS) */
#   define ses_vp_temp_from_var(x,seos)	ses_rt_temp_from_grid(x,seos)

		/* Macros to access solution functions */

#define ses_solution(coords,comp,hs,side,front,wave,state)		\
    (*(front)->_hyp_solution)(coords,comp,hs,side,front,wave,state,NULL)

#define ses_grad_solution(coords,comp,hs,side,front,wave,grd_st)	\
    (*(wave)->hyp_grad_soln)(coords,comp,hs,side,front,wave,grd_st)

enum {
	ISOTHERM_BOUNDARY = FIRST_BOUNDARY_EOS_WAVE_TYPE,
	ISOCHOR_BOUNDARY,
	PHASE_BOUNDARY    = FIRST_INTERIOR_EOS_WAVE_TYPE
};

enum {
	EOS_BOUNDARY_NODE = FIRST_EOS_NODE_TYPE,
	PHASE_BDRY_NODE
};

#if defined(PHASE_CODE)

/* Data on the phase boundary */           


	struct _PHASE_BDRY {
	    int dim;  /* Number of raw data points on phase boundary */
	    int place[2]; /*Positions for crossing of phase bound at tmin */
	    int n_pts_dome;
	    double *rvar; /* Transformed density along phase boundary */
	    double *Tvar, *slphT; /* Transformed temp along phase boundary */
	    double *rp, *slphrp; /* Reduced pressure along phase boundary */
	    double *re, *slphe; /* Reduced energy along phase boundary */
	    double *cp; /* Cold pressure along phase boundary */
	    double *ce; /* Cold energy along phase boundary */
	    double *S;
	};
	typedef struct _PHASE_BDRY PHASE_BDRY;

	/* Data and spline information on the cold curve */

	struct _COLD_CURVE {
	    double *rvar;
	    double *cp, *slcp;
	    double *ce, *slce;
	};
	typedef struct _COLD_CURVE COLD_CURVE;

#endif /* defined(PHASE_CODE) */

#include <geos/gsesprotos.h>

#endif /* !defined(_SESAME_H) */
