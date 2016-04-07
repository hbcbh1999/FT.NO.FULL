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
*				gstate.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the structures for gas dynamics states.
*/

#if !defined(_GSTATE_H)
#define _GSTATE_H

#include <driver/ddecs.h>

#if !defined(SMAXD)
enum {SMAXD = MAXD};
#endif /* !defined(SMAXD) */

	/* Artificial viscosity and heat conduction parameters */

#define	use_linear_artificial_viscosity(avisc)	((avisc).use_lin_av)
#define	use_lapidus_artificial_viscosity(avisc)	((avisc).use_lapidus_av)
#define	use_upwind_artificial_viscosity(avisc)	((avisc).use_upwind_av)
#define	use_muscl_slope_flattening(avisc)	((avisc).use_msf)
#define lapidus_stability_factor(visc_nl)				\
	( 1.0 / (sqrt(1.0 + 0.25*sqr(visc_nl)) - 0.5*(visc_nl)) )

	/* Prompts used for printing and restart of AVISC params */

#define	Use_nlav	"Use nonlinear artificial viscosity = "
#define	Coef_nlav	"Coefficient of nonlinear artificial viscosity = "
#define	Use_lav		"Use linear artificial viscosity = "
#define	Coef_lav	"Coefficient of linear artificial viscosity = "
#define	Use_uwav	"Use upwind artificial viscosity = "
#define	Coef_uwav	"Coefficient of upwind artificial viscosity = "
#define	Use_msf		"Use MUSCL slope flattening = "
#define	Coef_msf_ieta	"Muscl slope flattening parameter eta inverse = "
#define	Coef_msf_ms	"Muscl slope flattening minimum shock strength = "
#define	Coef_msf_msvj	"Muscl slope flattening minimum shock specific volume jump = "
#define	Coef_char_speed_cutoff						\
			"Muscl slope flattening charateristic speed cutoff  = "
#define	Coef_hc		"Coefficient of artificial heat conduction = "
#define	Coef_dst	"Coefficient of dynamic surface tension = "
#define	Coef_sp		"Artificial viscosity stability coefficient = "
#define Coef_contact_detector                                           \
                        "Contact detector coefficient = "

typedef	struct
{
	/* Flag for type of artificial viscosity */

	const char *hyp_method;
	boolean	   use_lin_av;
	boolean	   use_lapidus_av;
	boolean	   use_upwind_av;
	boolean	   use_msf;
	double	   linear_visc_coef;	/* coef. of lin. art. visc. */
	double	   lapidus_visc_coef;	/* coef. of Lapidus art. visc. */
	double	   upwind_visc_coef;	/* coef. of upwind art. visc */
	double	   msf_ieta;		/* Slope flattening coef.	*/
	double	   min_shock_jump;	/* Captured shock threshold */
	double	   min_sp_vol_jump;	/* Tolerance for sp. vol. jump */
	double	   heat_cond;		/* coef. of art. heat cond. */
	double	   sp_coef;		/* art. visc. stability factor */
	double      char_speed_cutoff;   /* No msf unless
					   lambdaR*lambdaL/cR*cL <= cutoff */
	double      dynamic_st;             /* Dynamic surface tension coef.
					   See set_pjump_at_wave() for desc. */
        double      contact_detector;    /* Dimensionless contact cutoff    */
                                        /* Corresponds to K0 in PPM paper  */
        double      viscosity_coef;      /* N-S viscosity coeff */
	boolean       use_stokes_vis;      /* Use Stokes viscosity*/
	double      diffusivity_coef;    /* N-S mass diffusion coeff */
	double      conductivity_coef;   /* N-S thermal conductivity coeff */
	boolean       subgrid_vis;         /* Subgrid model for viscosity*/
	boolean       subgrid_md;         /* Subgrid model for mass diffusion*/
        boolean       subgrid_con;
        boolean       aver_planar;
        boolean       aver_sphere;
        boolean       turbulence;
} AVISC;

		/* various gas parameters */

struct _Gas_param {
	struct _EOS *eos;
	const char  *prt_label;
	int         dim;
	size_t      sizest;
	int         n_comps;
	AVISC	    avisc;
	boolean        (*_invalid_state)(const char*,Locstate,boolean);
	void        (*_fprint_state)(FILE*,Locstate);
	void        (*_verbose_fprint_state)(FILE*,const char*,Locstate);
	void        (*_fprint_Gas_param)(FILE*,struct _Gas_param*);
	void        (*_alloc_state)(Locstate*,size_t);
	Locstate    (*_alloc_intfc_state)(size_t);
	void        (*_clear_state)(Locstate,size_t);
	void        (*_set_state)(Locstate,int,Locstate);
	void        (*_set_params)(Locstate,Locstate);

	/* Tolerances for general flows */

	double       min_energy;
	double       min_pressure;
	double       vacuum_dens;
	double       raref_press;

#if defined(COMBUSTION_CODE)

	/* Tolerances for combustion */
	double       tol_alpha;
	double       tol_press;
#endif /* defined(COMBUSTION_CODE) */


#if defined(COMBUSTION_CODE)
	int               composition_type;
	double             critical_temperature;
	int               burned;	/* BURNED or UNBURNED */
	double             q;	/* heat of combustion */
		        /* = 0 if burned = BURNED */
		        /* >= 0 if burned = UNBURNED */
	double             rate_mult;	/* multiplier for reaction rate */
	double             activ_en;	/* activation energy */
	char              code;         /* Flame velocity character code */
	double             flame_coeff[3];   /* Srabasti: coeffs of flame velo*/
	int               ncoeffs;      /* Number of used coeffs */
	double             (*_flame_velocity)(Locstate); 
				/* User provided flame velocity function */
	struct _Gas_param *other_params;
		/* "other" gas parameters corresponding to the alternate
		* phase of gas.  For ZND and TWO_CONSTITUENT_REACTIVE,
		* params == unburned params and
		* other_params == burned params */ 
	double             speed_corr, p_corr, u_corr, rho_corr;

#endif /* defined(COMBUSTION_CODE) */
};
typedef struct _Gas_param Gas_param;

#define Num_gas_components(st)						\
    Params(st)->n_comps
#define fprint_gas_state(file,st)					\
    (Params(st) != NULL) ?						\
        (*Params(st)->_fprint_state)((file),(st)) : g_fprint_state((file),(st))
#define print_gas_state(st)						\
    fprint_gas_state(stdout,st)
#define verbose_fprint_state(file,name,st)				\
    (Params(st) != NULL) ?	\
        (*Params(st)->_verbose_fprint_state)((file),(name),(st)) : \
	g_verbose_fprint_state((file),(name),(st))
#define verbose_print_state(name,st)					\
    verbose_fprint_state(stdout,name,st)
#define fprint_Gas_param(file,params)	((params) != NULL) ?		\
    (*(params)->_fprint_Gas_param)((file),(params)) :			\
    g_fprint_Gas_param((file),(params))
#define print_Gas_param(params)						\
    fprint_Gas_param(stdout,params)
#define print_avisc_structure(avisc,interactive)			\
    fprint_avisc_structure(stdout,avisc,interactive)
#define set_state(st1,type_st1,st2)					\
    (Params(st2) != NULL) ?						\
        (*(Params(st2))->_set_state)(st1,type_st1,st2) :		\
	g_set_state(st1,type_st1,st2)
#define copy_state(st1,st2)						\
    set_state(st1,state_type(st2),st2)
#define invalid_state(fnc,st,print_warning)				\
    ((Params(st) != NULL) && (Params(st)->_invalid_state != NULL)) ? \
	(*Params(st)->_invalid_state)(fnc,st,print_warning) : NO

		/*
		*  Possible representations of states
		*
		*  Below:
		*		rho    = density
		*		E      = total energy/(unit volume)
		*		e      = specific internal energy
		*		P      = pressure
		*		T      = temperature
		*		V      = 1/rho = specific volume
		*		S      = specific entropy
		*		lambda = reaction progress
		*/

					/* State representation */
enum _GAS_STATE_TYPE {
	OBSTACLE_STATE = 0xff,	/* state behind wall */
	UNKNOWN_STATE  = 0x00,	/* NONE */
	GAS_STATE      = 0x01,	/* rho, E, momentum */
	EGAS_STATE     = 0x02,	/* rho, e, velocity */
	FGAS_STATE     = 0x03,	/* rho, T, velocity */
	TGAS_STATE     = 0x04,	/* rho, P, velocity */
	VGAS_STATE     = 0x05,	/* rho, P, velocity + redundant vars */
	ZGAS_STATE     = 0x06,	/* rho, E, momentum, rho*lambda */
	CGAS_STATE     = 0x07	/* rho, E, momentum, rho1 */
};
typedef enum _GAS_STATE_TYPE GAS_STATE_TYPE;

		/* possible values for the fluid composition type */

enum _COMPOSITION_TYPE {
	PURE_NON_REACTIVE       = 0,
	MULTI_COMP_NON_REACTIVE,
	PTFLAME,
	ZND,
	TWO_CONSTITUENT_REACTIVE,
	THINFLAME
};
typedef enum _COMPOSITION_TYPE COMPOSITION_TYPE;

		/* choice of thermodynamic phases for the gas */

enum {
	BURNED	 = 1,
	UNBURNED = 2
};

		/* possible values for interpolation coefficient weights */

enum _INTERPOLATION_WEIGHT {
	UNWEIGHTED_INTERP_COEFS		  = 1,
	VOLUME_WEIGHTED_INTERP_COEFS	  = 2,
	RGAS_WEIGHTED_INTERP_COEFS	  = 3,
	RGAS_VOLUME_WEIGHTED_INTERP_COEFS = 4
};
typedef enum _INTERPOLATION_WEIGHT INTERPOLATION_WEIGHT;


		/* macros */

	/* Macro for accessing equation of state data for all states */
#define Params(state)		((Gas *) (state))->params
#define Set_params(s1,s2)						\
	if (Params(s2) == NULL) 					\
	    Params(s1) = NULL;						\
	else								\
	    (Params(s2)->_set_params)(s1,s2)
#define Init_params(state,params)	Params(state) = params
#define Same_params(s1,s2)	(Params(s1) == Params(s2))
#define Different_params(s1,s2)	(Params(s1) != Params(s2))
#if defined(COMBUSTION_CODE)
#define Set_other_params(s1,s2)	Params(s1) = Params(s2)->other_params
#define	Unburned(state)		(Params(state)->burned == UNBURNED)
#define	Burned(state)		(Params(state)->burned == BURNED)
#define	Heat_release(state)		(Params(state)->q)
#define	Rate_multiplier(state)	Params(state)->rate_mult
#define	Other_rate_multiplier(state)	Params(state)->other_params->rate_mult
#endif /* defined(COMBUSTION_CODE) */
#define	Composition_type(s)	Params(s)->composition_type

#if !defined(UNRESTRICTED_THERMODYNAMICS)
#define floor_pressure(p,pmin)	max((p),(pmin))
#else /* !defined(UNRESTRICTED_THERMODYNAMICS) */
#define floor_pressure(p,pmin)	(p)
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	/* Macro for accessing density of all states */

#define Dens(state)	((Gas *) (state))->rho
#if defined(SUBGRID)
#define C(state)     ((Gas *) (state))->c
#define CI(state)     ((Gas *) (state))->ci
#define SCT(state)     ((Gas *) (state))->sct 
#define SCT0(state)     ((Gas *) (state))->sct0
#define PRT(state)     ((Gas *) (state))->prt
#define STR(state)     ((Gas *) (state))->str
#define Tau(state)     ((Gas *) (state))->tau
#define Qh(state)     ((Gas *) (state))->qh
#define Qc(state)     ((Gas *) (state))->qc
#define DCSn(state)     ((Gas *) (state))->dcsn
#define DCIn(state)     ((Gas *) (state))->dcin
#define CSn(state)     ((Gas *) (state))->csn
#define CIn(state)     ((Gas *) (state))->cin
#define CSt1(state)     ((Gas *) (state))->cst1
#define DCSt1(state)     ((Gas *) (state))->dcst1
#define Svisx(state)     ((Gas *) (state))->svisx
#define Svisy(state)     ((Gas *) (state))->svisy
#define Svisz(state)     ((Gas *) (state))->svisz
#define Sethermal(state)     ((Gas *) (state))->sethermal
#define Scmd(state)     ((Gas *) (state))->scmd
#endif /* defined SUBGRID */
#define Dvnnn(state)     ((Gas *) (state))->dvnnn
#define Dvnn(state)     ((Gas *) (state))->dvnn
#define Dvt1nn(state)     ((Gas *) (state))->dvt1nn
#define Dvt2nn(state)     ((Gas *) (state))->dvt2nn
#define Dvt1n(state)     ((Gas *) (state))->dvt1n
#define Dvt2n(state)     ((Gas *) (state))->dvt2n
#define Mvisx(state)     ((Gas *) (state))->mvisx
#define Mvisy(state)     ((Gas *) (state))->mvisy
#define Mvisz(state)     ((Gas *) (state))->mvisz
#define Evis(state)     ((Gas *) (state))->evis
#define Ethermal(state)     ((Gas *) (state))->ethermal
#define Cmd(state)     ((Gas *) (state))->cmd
#define Emd(state)     ((Gas *) (state))->emd
#define Dvt1t1(state)     ((Gas *) (state))->dvt1t1
#define Dvt2t1(state)     ((Gas *) (state))->dvt2t1

	/* Macros for accessing total energy and momentum of Gas states */
	
#define Energy(state)	((Gas *) (state))->e
#define Mom(state)	(((Gas *) (state))->m)

	/* Macros for accessing pressure and velocity of TGas or VGas states */

#define Press(state)	((Gas *) (state))->e
#define Vel(state)	(((Gas *) (state))->m)
#define reset_gamma(state)	((Gas *) (state))->gamma_set = NO
#define Local_gamma(state)      ((Gas *) (state))->gamma
#define Local_gamma_set(state)  ((Gas *) (state))->gamma_set

	/* State type representation */

#define	set_type_of_state(state,type)					\
	((((Gas *)(state))->_type = (type)),(((Gas *)(state))->_failed = NO))
#define state_type(state)	((unsigned int)((Gas *)(state))->_type)
#define material_failure(state) ((unsigned int)((Gas *)(state))->_failed)
#define	set_material_failure(state,y_n)	 (((Gas *)(state))->_failed = (y_n))

	/* Macros for accessing temperature of FGAS_STATE type state */

#define Temperature(state)	((Gas *) (state))->e

	/* Macros for accessing multicomponent state information */

#define pdens(state)	(&((MGas *) (state))->rho0)

	/* Macros for accessing additional VGas state information */

#define	Int_en(state)		((VGas *) (state))->e
#define	Entropy(state)		((VGas *) (state))->S
#define	Sound_speed(state)	((VGas *) (state))->c
#if defined(VERBOSE_GAS_PLUS)
#define	Temp(state)		((VGas *) (state))->T
#define	Enthalpy(state)		((VGas *) (state))->i
#endif /* defined(VERBOSE_GAS_PLUS) */
#if defined(PHASE_CODE)
#define	Wave_curve(state)	((VGas *) (state))->wave_cur
#endif /* defined(PHASE_CODE) */

	/* Macros for accessing ZGas state information */

#define React(state)	pdens(state)[0]
#define Prod(state)	pdens(state)[0]

	/* Macros for accessing binary mixture information */

#define	MassFraction(state)	pdens(state)[0]	/* Reaction variable lambda */
#define	DensityFraction(state)	pdens(state)[0]	/* rho*lambda above */

	/* Macros for accessing CGas state information */

#define Dens1(state)	pdens(state)[1]

#define is_obstacle_state(state)	(Params(state)==NULL)

		/* ordinary gas */

struct _Gas {
	double		rho;	   /* density */
	double		e;	   /* energy density 		- GAS_STATE
				   *  pressure 			- TGAS_STATE
				   *  specific internal  energy - EGAS_STATE
				   */
	double		m[SMAXD];  /* momentum density 		- GAS_STATE
				   *  velocity			- TGAS_STATE
				   *  velocity			- EGAS_STATE
				   */
	Gas_param	*params;   /* Structure defined EOS */
	double		gamma;
	double		cv;
	double		R;
	boolean		gamma_set;
	unsigned int	_type:8;   /* state type representation flag */
	unsigned int	_failed:8; /* flag for detection of material failure*/
	/* May 7 2004: Myoung-Nyoun Kim: gradient of temperature */
        double           dTdn;

        double           dvnnn;
        double           dvnn;
        double           dvt1n;
        double           dvt1nn;
        double           dvt2nn;
        double           dvt2n;
        double           mvisx;
        double           mvisy;
        double           mvisz;
        double           evis;
        double           ethermal;
        double           cmd[2];
        double           emd;
        double           dvt1t1;
        double           dvt2t1;
#if defined(SUBGRID)
        double           c;
        double           ci;
        double           sct;
        double           sct0;
        double           prt;
        double           str;
        double           tau[SMAXD][SMAXD];   /*turbulent stress */
        double           qh[SMAXD];
        double           qc[SMAXD*2];
        double           dcsn;
        double           dcin;
        double           csn;
        double           cin;
        double           cst1;
        double           dcst1;
        double           svisx;
        double           svisy;
        double           svisz;
        double           sethermal;
        double           scmd[2];
#endif /* defined SUBGRID */
};
typedef struct _Gas Gas;


enum { MAX_NUM_GAS_COMPS = 10 };

/* May 7 2004: Myoung-Nyoun Kim: gradient of temperature */
#define dTdn(state) ((Gas *)(state))->dTdn

		/* multiple component gas */

typedef struct {
	Gas   gas;
	double rho0;	/* head of the partial density array */
} MGas;

		/* Gas structure for uni_arrayization */

struct _Vec_Gas {
	double       *rho;	/* Mass Density */
	double       *e;         /* Specific internal energy */
	double       *re;        /* Internal energy density */
	double       *en_den;    /* Total Energy density */
	double       *p;         /* Thermodynamic pressure */
	double       *c;         /* Thermodynamic sound speed */
	double       *GAM;	/* Gruneisen coefficient */
	double       *c2;	/* Square of sound speed */
	double       *FD;	/* Fundamental derivative */
	double       **m;	/* momentum density 0 = normal */
	double       **v;	/* velocity 0 = normal */
	double       **rho0;     /* rho0[i] = array of fractional densities */
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double       *min_pressure;
	double       *min_energy;
	double       *vacuum_dens;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int         nprms;
	int         *prms_jmp;
	Locstate    *state;
	double       **coords;
	const double* const *Q;	/* coordinate transformation */
	struct {
	    boolean rho;
	    boolean e; 
	    boolean re; 
	    boolean en_den;
	    boolean p;
	    boolean c;
	    boolean GAM;
	    boolean c2;
	    boolean FD;
	    boolean m;
	    boolean v;
	    boolean rho0;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    boolean min_pressure;
	    boolean min_energy;
	    boolean vacuum_dens;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    boolean prms_jmp;
	    boolean state;
	    boolean coords;
	    boolean vst;
	} alloc;
	struct {
	    boolean rho;
	    boolean e; 
	    boolean re; 
	    boolean en_den;
	    boolean p;
	    boolean c;
	    boolean GAM;
	    boolean c2;
	    boolean FD;
	    boolean m;
	    boolean v;
	    boolean rho0;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    boolean min_pressure;
	    boolean min_energy;
	    boolean vacuum_dens;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    boolean prms_jmp;
	    boolean state;
	    boolean coords;
	} set;
};
typedef struct _Vec_Gas Vec_Gas;

#define Vec_Gas_field_allocated(vst,ptr)	(vst)->alloc.ptr
#define Vec_Gas_field_set(vst,ptr)		(vst)->set.ptr


		/* Source uni_arrays */

struct _Vec_Src {
	double    *mass;   /* mass[j] = mass source in cell j                */
	double    *energy; /* energy[j] = total energy source in cell j      */
	double    **mom;   /* mom[i][j] = ith component of momentum source   */
	                  /*             in cell j                          */

	double    *radii, rmin;  /* For calculating cylindrical source terms */
	double    **rho0;  /* rho0[i][j] = source for fractional density i   */
	                  /* in cell j                                      */

	struct {
	    boolean mass;
	    boolean energy;
	    boolean mom;
	    boolean radii;
	    boolean rho0;
	    boolean src;
	} alloc;
};
typedef struct _Vec_Src Vec_Src;

enum {
	MAX_NUM_WAVES  = 5,
	NO_PTS_ON_COMP = 100
};

/* Wave curve information for the Riemann solver */

	/* possible types of waves in Riemann solver */

enum _RIEMANN_SOLVER_WAVE_TYPE {
	UNSET_RIEMANN_SOLVER_WAVE_TYPE = -1,
	RAREFACTION = 1,
	SHOCK,
	COMPOSITE,
	STRONG_DET,
	CJ_DET,
	LAST_GAS_RIEMANN_SOLVER_WAVE_TYPE
};
typedef enum _RIEMANN_SOLVER_WAVE_TYPE RIEMANN_SOLVER_WAVE_TYPE;

#if defined(PHASE_CODE)
typedef struct {
	Locstate st[MAX_NUM_WAVES];
	double ml[MAX_NUM_WAVES], mr[MAX_NUM_WAVES];
	int w_type[MAX_NUM_WAVES];
	int num_waves;
	double pstart[NO_PTS_ON_COMP], pend[NO_PTS_ON_COMP];
	double rstart[NO_PTS_ON_COMP], rend[NO_PTS_ON_COMP];
	double ustart[NO_PTS_ON_COMP], uend[NO_PTS_ON_COMP];
	double rhoc[NO_PTS_ON_COMP];
	int special;
} WAVE_CURVE;
#endif /* defined(PHASE_CODE) */

	/*	 verbose gas, verbose version of TGas  */

typedef struct {
	Gas gas;

	double rho_list[MAX_NUM_GAS_COMPS];	/* density of components */

		/* Additional information */

	double e;	/* specific internal energy */
	double S;	/* specific entropy */
	double c;	/* sound speed */
#if defined(VERBOSE_GAS_PLUS)
	double T;	/* temperature */
	double i;	/* specific enthalpy */
#endif /* defined(VERBOSE_GAS_PLUS) */
#if defined(PHASE_CODE)
	WAVE_CURVE *wave_cur;
#endif /* defined(PHASE_CODE) */
} VGas;

typedef struct {
	int state_type;
	Gas_param **params;
	double *rho;
	double *p;       /* pressure */
	double *e;       /* specific internal energy */
	double *S;       /* specific entropy */
	double *c;       /* sound speed */
	double *T;       /* temperature */
	double *i;       /* specific enthalpy */
	double *gam;
	double *GAM;
} Vec_Thermo;

struct _MODE {
	double amplitude,phase;
			/* dim-1 vector quantity */
	double wave_number[2];

	double growth_rate;
	double propagation_speed;	/* For KH random surface	*/
};
typedef struct _MODE MODE;

enum _CONSISTENT_PARAMS {
	PARAMS_CONSISTENT,
	PARAMS_ALL_OBSTACLE_STATES,
	PARAMS_INCONSISTENT,
	OBSTACLE_STATE_FOUND
};
typedef enum _CONSISTENT_PARAMS CONSISTENT_PARAMS;

	/* Time dependent pressure data */
struct _FD_DATA {
	double   tr;	/* rise time */
	double   tp;	/* peak time */
	double   ts;	/* shut-off time */
	double   pr_a;	/* ambient pressure */
	double   pr_p;	/* peak pressure */
	Locstate state; /* reference state */
};
typedef struct _FD_DATA FD_DATA;

struct _TD_BSTATE {
        int    n;/*number of elements*/
        double  *t;
        double  *p;
        double  *rho;
        double  *v;/*velocity in normal direction from exterior to interior*/
	char   fn[1029];
        Gas_param *params;
};
typedef struct _TD_BSTATE TD_BSTATE;

#include <gdecs/geos.h>

#endif /* !defined(_GSTATE_H) */
