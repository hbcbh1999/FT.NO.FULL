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
*				gdecs.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the structures for gas dynamics.
*/

#if !defined(_GDECS_H)
#define _GDECS_H

#include <driver/ddecs.h>
#include <gdecs/guserint.h>
#include <gdecs/grdecs.h>

enum _SCAT_WV_TOLERANCE_TYPE {
	NEVER_TRACK,
	ALWAYS_TRACK,
	PRESSURE_RATIO,
	ATWOOD_NUMBER,
	MACH_NUMBER_AHEAD
};
typedef enum _SCAT_WV_TOLERANCE_TYPE SCAT_WV_TOLERANCE_TYPE;

struct _SCAT_WV_TOL {
	const char             *wave_name;
	const char             *tol_type_name;
	double	               wave_tol;
	SCAT_WV_TOLERANCE_TYPE tol_type;
};
typedef struct _SCAT_WV_TOL SCAT_WV_TOL;

struct _UNTRACK_NODE_FLAG {
	const char *_node_name;
	boolean	   _untrack;
};
typedef struct _UNTRACK_NODE_FLAG UNTRACK_NODE_FLAG;

struct _EXTREME_VALUES {
	double  rho_max,  coords_rho_max[3],  rho_min,  coords_rho_min[3];
	double    p_max,    coords_p_max[3],    p_min,    coords_p_min[3];
	double    e_max,    coords_e_max[3],    e_min,    coords_e_min[3];
	double   ke_max,   coords_ke_max[3],   ke_min,   coords_ke_min[3];
	double    E_max,    coords_E_max[3],    E_min,    coords_E_min[3];
	double m_max[3], coords_m_max[3][3], m_min[3], coords_m_min[3][3];
	double v_max[3], coords_v_max[3][3], v_min[3], coords_v_min[3][3];
};
typedef struct _EXTREME_VALUES EXTREME_VALUES;

struct _G_Front {
	D_Front	dfront;

	struct _G_FRONT_EXTENSION {
	    /* Minimum wave strengths for tracking of scattered waves
	     *  during bifurcation. It is indexed by node type, wave class,
	     *  and curve status. */
	    SCAT_WV_TOL *_scat_wv_tol_list;

	    void (*_operator_split_point_propagate)(Front*,POINTER,POINT*,
	                                            POINT*,HYPER_SURF_ELEMENT*,
						    HYPER_SURF*,double,double*);
	    void (*_unsplit_point_propagate)(Front*,POINTER,POINT*,POINT*,
	                                     HYPER_SURF_ELEMENT*,HYPER_SURF*,
					     double,double*);
	    void (*_set_USWSSten2d)(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			            HYPER_SURF*,Front*,Wave*,double);

	    /* one dimensional scheme */
	    void (*_oned_tangential_scheme)(int,int*,int*,Wave*,Wave*,Front*,
		                            Front*,Stencil*,Tan_stencil*,
					    int,int,Vec_Gas*,Vec_Src*,double,
					    double,int);

	    WSSten *(*_AllocWSSten)(int,int,Front*);
	    int	_nsts_WSSten;
	    int	_stype_WSSten;

#if defined(FULL_PHYSICS)
	    /*
	     * Yes/no flag for use of untracking when node propagation
	     * algorithm fails.  Indexed by node type.
	     */
	    UNTRACK_NODE_FLAG *_untrack_node_list;
#endif /* defined(FULL_PHYSICS) */
	    boolean _include_wall_normal_velocity;
	    double _curvature_factor;

	    EXTREME_VALUES _extreme_front_vals;
	} g_extension;
};
typedef struct _G_Front G_Front;

#define g_front(front)			((G_Front *)(front))
#define	GasFrontExtension(front)	g_front(front)->g_extension

#define scat_wv_tol_list(front)						\
		GasFrontExtension(front)._scat_wv_tol_list
#define extreme_front_vals(front)					\
		GasFrontExtension(front)._extreme_front_vals
#define	g_front_oned_tangential_scheme(front)				\
		GasFrontExtension(front)._oned_tangential_scheme
#define scat_wv_index(nt,wc,ns)						\
    (ns + NUM_NODE_STATUS*(wc + NUM_WAVE_CLASSES*(nt-FIRST_PHYSICS_NODE_TYPE)))
#define	scattered_wave_tolerance(front,nt,wc,ns)			\
    (scat_wv_tol_list(front)+scat_wv_index(nt,wc,ns))
#define tracked_oned_scattered_wave(family,s0,s1,front)			\
	track_scattered_wave(ONED_INTERACTION,family,ONED_WAVE,s0,s1,front)

#define	OperatorSplitPointPropagate(front)				\
    GasFrontExtension(front)._operator_split_point_propagate
#define	UnsplitPointPropagate(front)					\
    GasFrontExtension(front)._unsplit_point_propagate
#define	WSStenAllocator(front)						\
    GasFrontExtension(front)._AllocWSSten
#define	SetUSWSStencil(front)						\
    GasFrontExtension(front)._set_USWSSten2d
#define untrack_node_list(front)					\
    GasFrontExtension(front)._untrack_node_list

#define	untrack_node_index(nt)						\
    (nt-FIRST_PHYSICS_NODE_TYPE)
#define	untrack_node(front,nt)						\
    untrack_node_list(front)[untrack_node_index(nt)]._untrack
#define	untrack_node_name(front,nt)					\
    untrack_node_list(front)[untrack_node_index(nt)]._node_name
#define include_wall_normal_velocity(front)				\
    GasFrontExtension(front)._include_wall_normal_velocity
#define curvature_factor(front)						\
    GasFrontExtension(front)._curvature_factor

#define	operator_split_point_propagate(fr,wave,oldp,newp,oldhse,oldhs,dt,V) \
	(*OperatorSplitPointPropagate(fr))(fr,wave,oldp,newp,oldhse,oldhs,dt,V)

#define	unsplit_point_propagate(fr,wave,oldp,newp,oldhse,oldhs,dt,V)	\
	(*UnsplitPointPropagate(fr))(fr,wave,oldp,newp,oldhse,oldhs,dt,V)

#define	oned_tangential_scheme(fr,ts,of,vs,vst,src,dt,dn,dim)		\
    (*g_front_oned_tangential_scheme(fr))(ERROR,NULL,NULL,NULL,NULL,	\
                                          fr,NULL,NULL,ts,of,vs,vst,src,\
					  dt,dn,dim)

#define	nsts_WSSten(fr)							\
	GasFrontExtension(fr)._nsts_WSSten

#define	stype_WSSten(fr)						\
	GasFrontExtension(fr)._stype_WSSten

#define	AllocWSSten(rad,stype,fr)					\
	(*WSStenAllocator(fr))(rad,stype,fr)

#define	AllocDefaultWSSten(fr)						\
	(*WSStenAllocator(fr))(nsts_WSSten(fr),stype_WSSten(fr),fr)

#define	set_USWSSten2d(uswssten,p,hse,hs,fr,wave,dt)		\
	(*SetUSWSStencil(fr))(uswssten,p,hse,hs,fr,wave,dt)

/* values returned by load_state_vectors */

enum _TIME_DYNAMICS {
	CONSTANT_IN_TIME = 0,
	DYNAMIC_IN_TIME
};
typedef enum _TIME_DYNAMICS TIME_DYNAMICS;


struct _G_Wave {
	D_Wave	dwave;

	struct _G_WAVE_EXTENSION {

	    /* gather data from Wave into Vec_Gas/Vec_Src structures */
	    struct _G_USER_WAVE {
	        TIME_DYNAMICS (*_load_state_vectors)(int,int*,Vec_Gas*,int,int,
				                     Wave*,Wave*,int*,int);

	    /* one dimensional scheme */
	        void (*_oned_interior_scheme)(int,int*,int*,Wave*,Wave*,
			                      Front*,Front*,Stencil*,
					      Tan_stencil*,int,int,Vec_Gas*,
					      Vec_Src*,double,double,int);

	    /* scatter data from Vec_Gas into Wave structure */

	        void (*_assign_wave_state_vectors)(int,int*,Wave*,Wave*,
		                                   Vec_Gas*,int,int,int*,int);
	        EXTREME_VALUES _extreme_wave_vals;
	    } g_user_wave;

	    struct _G_WAVE_POINTERS {
	    	Vec_Gas *vgas;
	    	Vec_Src *vsrc;
	    } g_wave_pointers;
	} g_extension;
};
typedef struct _G_Wave G_Wave;
#if defined(__cplusplus)
typedef struct _G_Wave::_G_WAVE_EXTENSION::_G_USER_WAVE G_USER_WAVE;
typedef struct _G_Wave::_G_WAVE_EXTENSION::_G_WAVE_POINTERS G_WAVE_POINTERS;
#else /* defined(__cplusplus) */
typedef struct _G_USER_WAVE G_USER_WAVE;
typedef struct _G_WAVE_POINTERS G_WAVE_POINTERS;
#endif /* defined(__cplusplus) */

#define	g_wave(wave)		((G_Wave *)(wave))
#define	GasWaveExtension(wave)	g_wave(wave)->g_extension
#define	g_user_wave(wave)	GasWaveExtension(wave).g_user_wave
#define	g_wave_load_state_vectors(wave)				\
				g_user_wave(wave)._load_state_vectors
#define	g_wave_oned_interior_scheme(wave)			\
				g_user_wave(wave)._oned_interior_scheme
#define	g_wave_assign_wave_state_vectors(wave)			\
				g_user_wave(wave)._assign_wave_state_vectors
#define extreme_wave_vals(wave)					\
		g_user_wave(wave)._extreme_wave_vals
#define	g_wave_pointers(wave)	GasWaveExtension(wave).g_wave_pointers
#define	g_wave_vgas(wave)	g_wave_pointers(wave).vgas
#define	g_wave_vsrc(wave)	g_wave_pointers(wave).vsrc

#define	load_state_vectors(sn,ip,vst,imn,imx,wv,newwv,ic,pb)		\
	(*g_wave_load_state_vectors(wv))(sn,ip,vst,imn,imx,wv,newwv,ic,pb)

#define	oned_interior_scheme(sn,ip,ic,wv,nwv,fr,nfr,stn,of,vs,vst,src,dt,dn,dim)\
    (*g_wave_oned_interior_scheme(wv))(sn,ip,ic,wv,nwv,fr,nfr,stn,NULL,	\
                                       of,vs,vst,src,dt,dn,dim)

#define	assign_wave_state_vectors(sn,ip,wv,nwv,vst,imn,imx,ic,pb)	\
	(*g_wave_assign_wave_state_vectors(wv))(sn,ip,wv,nwv,vst,imn,imx,ic,pb)

/*
*	Initialization structure for gas library
*/

struct _G_INIT_PHYSICS {
	INIT_PHYSICS    Iphys;
	boolean            _supports_riemann_problem_waves;
	Gas_param*	(*_prompt_for_eos_params)(INIT_DATA*,INIT_PHYSICS*,
						  boolean,const char*);
	Gas_param**     (*_prompt_for_eos_params_list)(INIT_DATA*,INIT_PHYSICS*,
	                                               boolean,int*);
	Gas_param*	(*_init_eos_params)(INIT_DATA*,INIT_PHYSICS*,
					    const char*,boolean);
	int    _problem;
	int    (*_init_composition_type)(INIT_PHYSICS*,INIT_DATA*,size_t*,int*);
	void   (*_init_printing_and_plotting)(INIT_PHYSICS*,INIT_DATA*,
					      Printplot*,int);
	PRINTING_LIST	*(*_set_printing_list)(INIT_PHYSICS*,INIT_DATA*,
					       Printplot*);
	INPUT_SOLN	**(*_set_input_solution)(Printplot*);
	OUTPUT_VALUE	*(*_phys_solution)(OUTPUT_SOLN*,double*,int*);
	void	(*_phys_intfc_solution)(OUTPUT_SOLN*,POINT*,HYPER_SURF_ELEMENT*,
					HYPER_SURF*,OUTPUT_VALUE*,
					OUTPUT_VALUE*);
	int	(*_prompt_for_boundary_state)(int,const char*,double*,
					      COMPONENT,int,HYPER_SURF*,
					      INIT_DATA*,INIT_PHYSICS*);
	int	(*_prompt_for_wave_type)(const char*,INTERFACE*,INIT_PHYSICS*);
	void	(*_init_hyperbolic_method)(INIT_DATA*,INIT_PHYSICS*);
	void	(*_init_problem_type)(INIT_PHYSICS*);
	void	(*_prompt_for_equation_of_state)(INIT_DATA*,Gas_param**,
						 const char*,
						 const char*,INIT_PHYSICS*);
	void	(*_prompt_for_problem_specific_data)(INIT_DATA*,INIT_PHYSICS*);
	void	(*_set_basic_phys_parameters)(INIT_DATA*,INIT_PHYSICS*);
};
typedef struct _G_INIT_PHYSICS G_INIT_PHYSICS;

#define	g_iphys(ip)	((G_INIT_PHYSICS *) (ip))

#define	supports_riemann_problem_waves(ip)				\
	g_iphys(ip)->_supports_riemann_problem_waves

#define	prompt_for_eos_params(init,ip,prompt_for_visc,mesg)		\
	(*g_iphys(ip)->_prompt_for_eos_params)(init,ip,prompt_for_visc,mesg)

#define prompt_for_eos_params_list(init,ip,pfv,num_eos)     \
        (*g_iphys(ip)->_prompt_for_eos_params_list)(init,ip,pfv,num_eos)

#define	init_eos_params(init,ip,mesg,prompt_for_visc)			\
	(*g_iphys(ip)->_init_eos_params)(init,ip,mesg,prompt_for_visc)

#define	problem_type(ip)	g_iphys(ip)->_problem

#define	init_composition_type(ip,init,psizest,pnfloats)			\
	(*g_iphys(ip)->_init_composition_type)(ip,init,psizest,pnfloats)

#define init_printing_and_plotting(ip,init,prt,nfloats)			\
	(*g_iphys(ip)->_init_printing_and_plotting)(ip,init,prt,nfloats)

#define	set_printing_list(ip,init,prt)					\
	(*g_iphys(ip)->_set_printing_list)(ip,init,prt)

#define	prompt_for_boundary_state(w_type,name,coords,comp,index,hs,init,ip)\
	(*g_iphys(ip)->_prompt_for_boundary_state)(w_type,name,coords,comp,\
						   index,hs,init,ip)

#define	init_hyperbolic_method(init,ip)					\
	(*g_iphys(ip)->_init_hyperbolic_method)(init,ip)

#define init_problem_type(ip)						\
	(*g_iphys(ip)->_init_problem_type)(ip)

#define prompt_for_equation_of_state(init,params,m1,m2,ip)	\
	(*g_iphys(ip)->_prompt_for_equation_of_state)(init,params,m1,m2,ip)

#define	prompt_for_problem_specific_data(init,ip)			\
	(*g_iphys(ip)->_prompt_for_problem_specific_data)(init,ip)

#define	set_basic_phys_parameters(init,ip)				\
	(*g_iphys(ip)->_set_basic_phys_parameters)(init,ip)

#define	prompt_for_wave_type(mesg,intfc,ip)				\
	(*g_iphys(ip)->_prompt_for_wave_type)(mesg,intfc,ip)

/* possible problem types */
enum {
	UNKNOWN_PROBLEM_TYPE = -3,
        UNSPECIFIED	     =  1,
	PLANE_FRONT,
	SHOCKED_THERMAL_LAYER,
	KELVIN_HELMHOLTZ,
	RAYLEIGH_TAYLOR,
	RANDOM_SURFACE,
#if defined(FULL_PHYSICS)
	MESHKOV,
	RICHTMYER_MESHKOV,
#endif /* defined(FULL_PHYSICS) */
	ASTROPHYSICAL_JET,
	BOWSHOCK,
	RAMP_REFLECTION,
	EXPANDING_SHOCK,
#if defined(FULL_PHYSICS)
	SHOCK_DIFFRACTION,
	SHOCK_TRANSMISSION,
#endif /* defined(FULL_PHYSICS) */
	TRI_GRID_TEST,
	RIEMANN2D,
#if defined(FULL_PHYSICS)
	CC_NODE_TEST,
#endif /* defined(FULL_PHYSICS) */
	AMBIENT_STATE_TEST,
	BUBBLES_DROPS,
	EXPANDING_SHELLS,
	SHOCK_JET,
	RICHTMYER_LINEAR_THEORY,
	SUPERNOVA,
	IMPLOSION,
	ONED_TEST,
	INJECTION_INLET,
	RADIAL_RAYLEIGH_TAYLOR,
	NE_BOOSTER,
	VELOCITY_TEST,
	FLUID_RIGID_BODY,
	DENSITY_STEP,
	LAST_GAS_PROBLEM_TYPE
};

#if defined(FULL_PHYSICS)

		/* physics values for propagation flag */

enum _G_NODE_FLAG_INDEX {
	_USE_SUBSONIC_STATE_INDEX	    = _FIRST_PHYSICS_NODE_FLAG_INDEX,
	_NODE_VEL_BY_ANGLE_INDEX,
	_DONT_CORRECT_ANGLES_AT_NODE_INDEX,
	_UNTRACK_SCATTERED_WAVES_AT_REFRACTION_PRECURSOR_INDEX,
	_NODE_WARNINGS_OFF,
	_FIRST_NON_GAS_NODE_FLAG_INDEX
};

#define use_subsonic_state(flag)					\
    (flag)._node_flags[_USE_SUBSONIC_STATE_INDEX]
#define node_vel_by_angle(flag)						\
    (flag)._node_flags[_NODE_VEL_BY_ANGLE_INDEX]
#define dont_correct_angles_at_node(flag)				\
    (flag)._node_flags[_DONT_CORRECT_ANGLES_AT_NODE_INDEX]
#define untrack_scattered_waves_at_refraction_precursor(flag)		\
    (flag)._node_flags[_UNTRACK_SCATTERED_WAVES_AT_REFRACTION_PRECURSOR_INDEX]
#define node_warnings_off(flag)						\
    (flag)._node_flags[_NODE_WARNINGS_OFF]

		/* values for tran_node_parameter_choice */

enum {
	USE_SLIP	   = 1,
	USE_INCIDENT_ANGLE = 2
};
#endif /* defined(FULL_PHYSICS) */

		/* values for init_bubble_parameter_choice */

enum {
	BOW_LENGTH	 = 1,
	TIME_ELAPSED	 = 2,
	REFL_WALL_LENGTH = 3
};

		/* choice of roots for Mach speed in Mach_node_speed() */

enum {
	SLOW = 1,
	FAST = 2
};

		/* Values for flag in bifurcation of reflections */

enum {
	NORMAL_TO_MACH_REFLECTION    = 1,
	NORMAL_TO_REGULAR_REFLECTION = 2,
	REGULAR_TO_MACH_REFLECTION   = 3
};


	/* Initialization info for the bubble in mach or regular reflections.
	   The states and angles around the reflection point are assumed to be
	   in a companion RP data structure. */

struct _Bubble {
		/* These states are the same for reg or mach. */
	Locstate	bow,		/* behind refl shock at wall */
			refl_corner,	/* need two for convex corner */
			bow_corner;	/* for concave, should be same */

		/* These states don't exist for reg bubbles. */
	Locstate	mach,			/* behind mach stem at wall */
			contact_bow,		/* behind slip for bow */
			contact_mach;		/* ahead of slip for mach */

	RP_DATA		*RP;			/* states and angles around */
						/* reflection point */
	COMPONENT       comp_bow,
			comp_mach;		/* doesn't exist for reg refl */

	double	        refl_posn[MAXD],	/* location of refl pt */
		        bow_posn[MAXD],		/* base of bow shock */
		        slip_posn[MAXD],	/* base of contact */
		        mach_posn[MAXD],	/* base of mach stem */
			corner_posn[MAXD];

	double		inc_t[MAXD],		/* away from node */
			aw_t[MAXD],		/* ahead wall tangent */
			bw_t[MAXD];		/* behind wall tangent */

	double		refl_length,		/* corner to base of mach */
			bow_length,		/* corner to bow */
			mach_height;

	double		node_v[MAXD];		/* node velocity */
	double		cor_v[MAXD];		/* corner velocity */
	double		mach_speed;		/* for base of Mach stem */
	double		contact_speed;		/* for base of contact */

	int             is_node_at_corner;

	int		is_attached;		/* is bow shock attached */
	double		bow_base_ang;		/* for attached bow */
	boolean		intfc_table_storage;
};
typedef struct _Bubble Bubble;

struct _RANDOM_STATE {

	Locstate mu;	/* mean state */
	Locstate sigma;	/* variance state */

	/* Specification of correlation ellipsoid */
	double *Q[3];/* Rotation matrix for the ellisoid axes */
	double lambda[3];/* correlation lengths */
	double *A[3];/* Bilinear form defining the correlation ellipsoid */

	int N;		/* number of evaluations per tau or lambda */
	double tau;	/* correlation time */
	double tlast;	/* time of last correlation computation */
	double delta_t;	/* time between correlation computations */

	int	M_ell;	/* number grid points in correlation ellipsoid */

	int correlated;

	/* Specification of rectangular arrays of states */
	RECT_GRID grid; /* Spatial grid for indep. and correlated states */

	double	time_of_save;

	/* Random number generator seeds */
	unsigned short int xsubi[3];
	unsigned short int seed_after_save[3];
	unsigned short int seed_after_old[3];

	/* Storage for random states grids */
	Locstate  ***indep_st;	   /* array of indepedent states */
	byte	  *indep_st_store; /* storage for array of indepedent states */
	Locstate  ***corr_st[3];   /* array of correlated states */
	byte      *corr_st_store;  /* storage for array of correlated states */
	Locstate  ***save_st, ***old_st, ***new_st;

	double	Qstore[3][3], Astore[3][3];

	double _RadialVelocityDecayExponent;
	double _RadialVelocityDecayScale;

};
typedef struct _RANDOM_STATE RANDOM_STATE;

#define	Mean(rstate)	((rstate)->mu)
#define	Sigma(rstate)	((rstate)->sigma)
#define	uncorrelated_random_state(state,rstate)				\
	(*(rstate)->_uncorrelated_random_state)(state,rstate)
#define RadialVelocityDecayExponent(rstate)				\
	((rstate)->_RadialVelocityDecayExponent)
#define RadialVelocityDecayScale(rstate)				\
	((rstate)->_RadialVelocityDecayScale)

struct _GraphUnits {
	double	_time_scale,     _initial_tmin, _initial_tmax;
	double	_length_scale,   _initial_lmin, _initial_lmax;
	double	_velocity_scale, _initial_vmin, _initial_vmax;
};
typedef struct _GraphUnits GraphUnits;


	/* possible initializations for interface instability problems */

enum _PERT_TYPE {
	UNSET_PERT_TYPE =       0,      /* Unset value            */
	COMPRESSIBLE =		1,	/* compressible theory	  */
	INCOMPRESSIBLE =	2,	/* incompressible theory  */
	VELOCITY_FIELD = 	3,	/* initial velocity field */
	RT =			4,	/* Rayleigh-Taylor	  */
	KELHELM =		7,	/* Kelvin-Helmholtz	  */
	RT_AMB =		8	/* R-T no linear analys	  */
};
typedef enum _PERT_TYPE PERT_TYPE;


enum _PERT_BDRY_TYPE {
	PERIODIC,
	SYMMETRIC,
	UNMODIFIED
};
typedef enum _PERT_BDRY_TYPE PERT_BDRY_TYPE;

struct _SINE_PERT {
	double		z_intfc, z_bdry;
	Locstate	amb_st;
	PERT_TYPE	init_type;
	PERT_BDRY_TYPE  pert_bdry_type[3];

		/* sigma occurs in the growth rate, as a multiple of t */
	double	stream_velocity[2];/* For KH_.. only	*/
	double	delta_v[2];	/* For KH_.. only	*/
	double	surf_ten;
	int	number_modes;
	MODE	*mode;

                /* The following added 24/2/2004 by egeorge.
		   To aid code readability, this addendum along 
		   with some of the original SINE_PERT members should 
		   be ft_assigned a separate struct called, for example,
		   PERT_FROM_FILE. */
        double    *x_coord;
        double    *y_coord;
        double    **amplitudes;
        int      Nx;
        int      Ny;
        int      read_amplitudes_from_file; /* This entry might no longer 
					       be needed if the preceding 
					       PERT_FROM_FILE structure is 
					       defined. */
};
typedef struct _SINE_PERT SINE_PERT;


enum _Gravity_type {
	NO_GRAVITY	       = 0, /* No gravity */
	CONSTANT_GRAVITY       = 1, /* constant in space and time */
	TIME_DEPENDENT_GRAVITY = 2, /* Spatially constant but time dependent */
	ASTROPHYSICAL_GRAVITY  = 3, /* g = 1/r^2 */
	RADIAL_GRAVITY         = 4, /* g = g_0 r,  r = radial vector */
	GENERALIZED_ASTROPHYSICAL_GRAVITY = 5,
	USER_DEFINED           = 6  /* user defined gravity */
};
typedef enum _Gravity_type Gravity_type;

struct _GRAVITY {
	Gravity_type	type;

	/* Parameters for constant gravity*/
	double	g[3];

	/* Parameters for time dependent gravity */
	double	**g_of_t;
	int	num_time_points;
	int	dim;

	/* Parameters for Astrophysical gravity */
	double	center[3];	/* center of gravity */
	double	G;		/* Gravity constant  */
	double	M;		/* mass of gravity center */
	int     N;              /* number of shells dividing the domain */
	double	R;		/* Maximum range of shells */
	double	dR;		/* width of of shells */
        double   *Mp;            /* partial mass of each sector */
	double   ***g_of_grid;   /* gravity at each center of grid block */
	CHART   ** roots;       /* Charts */
        int     nroots;         /* Number of charts */
	INIT_PHYSICS* ip;       /* Physics */
	int gmax[3];            /* Number of grid blocks */
        double GL[3];            /* Lower corner of global grid */
        double GU[3];            /* Upper corner of global grid */
  
	/* Arbitary user defined gravity */
	const double    *(*_user_defined_gravity)(const double*,const double,
						 struct _GRAVITY*);
};
typedef struct _GRAVITY GRAVITY;

IMPORT void g_set_gravity_charts(CHART**, int);
IMPORT boolean eval_partial_mass(double*,int,int);
		
#define	user_defined_gravity(coords,time,grav_data)			\
	(*(grav_data)->_user_defined_gravity)(coords,time,grav_data)

typedef enum {
	MOC_TYPE_UNSET   = -1,
	RIEMANN          =  1,
	MOC_PLUS_RIEMANN =  2,
	MOC_PLUS_RH      =  3
}	MOC_TYPE;

enum _GEOM_SOURCE_METHOD {
	MODIFIED_EULER,
	ANALYTIC_SOLUTION,
	BACKWARD_EULER
};
typedef enum _GEOM_SOURCE_METHOD GEOM_SOURCE_METHOD;


struct _NptWSpeedOpts {
	double    Mach_tol;
	double    A_tol;
	double    Wall_limiter;
	MOC_TYPE vector_moc;
	MOC_TYPE scalar_moc;
	boolean	 _scalar_filter_outgoing_waves_at_contact;
	boolean	 use_cheap_ahead_state_moc;
	boolean	 use_neumann_cheap_moc;
	void	 (*vector_ahead_state_moc)(double*,Locstate,Locstate,Locstate,
					   Locstate,Locstate,double,double,
					   double,double,double*,double*,int,
					   double,Front*);
	void	 (*neumann_moc)(double*,Locstate,double,double,Locstate,SIDE,
			        Locstate,double,double*,Front*);
	GEOMETRY_REMAP	   remap;
	GEOM_SOURCE_METHOD geom_source_method;
};
typedef struct _NptWSpeedOpts NptWSpeedOpts;

#include <gdecs/gbifur.h>
#include <gdecs/gprt.h>
#include <gdecs/guserrp.h>

/* Function prototypes */

#include <gbifur/gbifurprotos.h>
#include <geos/geosprotos.h>
#include <gintfc/gintfcprotos.h>
#include <gnode/gnodeprotos.h>
#include <gprop/gpropprotos.h>
#include <gprt/gprtprotos.h>
#include <gstate/gstateprotos.h>

#endif /* !defined(_GDECS_H) */
