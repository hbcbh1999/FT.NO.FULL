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
*				ghyp.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GHYP_H)
#define _GHYP_H

#include <gdecs/gdecs.h>

		/* eigenvalue and left/right eigenuni_arrays */

typedef struct {
	int	negn;
	double	**lambda; /* eigenvalues  */
	double	**sgnl;   /* 1 for left moving characteric 0 otherwise */
	double	**sgnr;   /* 1 for right moving characteric 0 otherwise */
	double	**sgnm;   /* algegbraic sign of average characteristic speed */
	double	***l;     /* left eigenuni_arrays */
	double	***r;     /* right engenuni_arrays */

	struct {
	    boolean lambda;
	    boolean sgnl;   /* 1 for left moving characteric 0 otherwise */
	    boolean sgnr;   /* 1 for right moving characteric 0 otherwise */
	    boolean sgnm;   /* algegbraic sign of average characteristic speed */
	    boolean l;     /* left eigenuni_arrays */
	    boolean r;     /* right engenuni_arrays */
	    boolean vegn;
	} alloc;
} Vec_Eigen;

typedef struct {
	double	**g;		/* artificial viscosity coefficients */
	double	*cs_ave;	/* Average sound speed squared */
	double	*c_ave;		/* Average sound speed */
	double	*vn_ave;	/* Average normal velocity */
	double	**b;		/* Wave speed differences times
				 * viscosity parameter */
	double	*visc;		/* Vector of viscosity parameters */
	double	**uconM;	/* Conservative form of mid state */
	double	*mdlambda;	/* max(lambda[i] - lambda[i+1]) */
	boolean	use_lapidus;	/* use Lapidus artificial viscosity */
	boolean	use_linear;	/* use linear artificial viscosity */
	boolean	use_upwind;	/* use upwind weighted artificial viscosity */

	struct {
	    boolean g;
	    boolean cs_ave;
	    boolean c_ave;
	    boolean vn_ave;
	    boolean b;
	    boolean visc;
	    boolean uconM;
	    boolean mdlambda;
	    boolean avisc;
	} alloc;
} Vec_Avisc;

typedef struct {
	double *chi;
	struct {
	    boolean chi;
	    boolean msf;
	} alloc;
} Vec_MSF;

/* CG_PARAMS contains different 'fudge' factors and parameters */

struct _CG_PARAMS {
	int pv_iterations;    /* Number of times to iterate in
	                       * Pressure-Velocity plane */
	double min_p_jump;
	double min_v_jump;    /* velocity convergence criteria */
	double sw_tol;
	double min_mass_flux;
	double sqr_min_mass_flux;
};
typedef struct _CG_PARAMS CG_PARAMS;

struct _MUSCL_FLUX {
	double	**F;	/* Volume flux vector F[k][j] is the volume flux of */
	                /* the kth state variable (0 = density, 1 = total   */
			/* energy density, 2 = sweep direction momemtum,    */
			/* etc,  see the definition of ucon above) across   */
			/* the cell boundary between cells j-1 and cell j   */
	double   **H;    /* Line flux uni_arrays, see the documentation above   */
};
typedef struct _MUSCL_FLUX MUSCL_FLUX;


/* User defined options for MUSCL method */
struct _Muscl_Opts {
	int	kmax;	/* Number of independent state variables */
	int     worksp_len;
	int     nfloats;
	boolean    monotone_reconstruction;
	boolean    link_reconstructions;

	/* Load data for MUSCL solvers */
	struct _Vec_Muscl* (*_load_state_data)(struct _Muscl_Opts*,int,int*,
                                               Front*,Wave*,Stencil*,
                                               Tan_stencil*,Vec_Gas*,Vec_Src*,
					       int,int,int,double,double);

	/* compute eigenuni_arrays for df/du */
	void (*_compute_eigens)(int,int,struct _Vec_Muscl*);

	/* compute artificial viscosity */
	void (*_compute_art_visc_coefs)(int,int,struct _Vec_Muscl*);

	/* State reconstruction function */
	void (*_reconstructor)(int,int,struct _Vec_Muscl*);

	/* Half-step propagation */
	void (*_half_step)(int,int,double,double,struct _Vec_Muscl*);
	void (*_strong_wave_half_step)(int,int,double,double,struct _Vec_Muscl*);

	/* Riemann solver */
	void (*_rsolver)(int,int,double**,Vec_Gas*,double**,Vec_Gas*,
                         double**,Vec_Gas*,MUSCL_FLUX*,struct _Vec_Muscl*);

	void (*_rmidstate)(int,int,Vec_Gas*,Vec_Gas*,
	                   double*,double*,struct _Vec_Muscl*);

	/* Method of characteristics */
	void (*_characteristic_solve)(int,int,struct _Vec_Muscl*);

	/* Compute flux uni_arrays */
	void (*_flux_vectors)(int,int,double**,double*,MUSCL_FLUX*,
                              struct _Vec_Muscl*);

	/* addition of artificial viscosity before half step */
	void (*_add_art_visc1)(int,int,double**,struct _Vec_Muscl*);

	/* addition of artificial viscosity after half step */
	void (*_add_art_visc2)(int,int,struct _Vec_Muscl*);

	/* Transforms src terms from noncons. to cons. form */
	void (*_cons_src)(int,int,int,int*,Tan_stencil*,struct _Vec_Muscl*);

	void (*_print_internal_energy)(const char*,double**,
				       struct _Vec_Muscl*,int,int);

	/* Stencil operators */

	void (*_npt_solver)(double,double,Locstate,const double*,
	                    int,int*,int*,Stencil*);
	void (*_one_side_npt_tang_solver)(double,double,Tan_stencil*,Locstate,
					  Front*);
	int               _npts_tan_sten;
	void (*_npt_tang_solver)(double,double,Tan_stencil*,
	                         Locstate,Locstate,Front*);
        void (*_npt_parab_tan_solver2d)(double,double,Tan_stencil*,
                                Locstate,Locstate,struct _Front*);
        void (*_npt_parab_tan_solver3d)(struct _Front*,const Tparams*,
                                const Tparams*,POINT*);

	/* Storage allocators */

	void (*_alloc_phys_vecs)(Wave*,int);
	void (*_free_phys_vecs)(Wave*);

	/* State Loaders */

	TIME_DYNAMICS (*_load_state_vectors)(int,int*,Vec_Gas*,int,int,
				             Wave*,Wave*,int*,int);
	void (*_assign_wave_state_vectors)(int,int*,Wave*,Wave*,
		                           Vec_Gas*,int,int,int*,int);

	CG_PARAMS	_Cg_params;

	struct _MUSCL_PromptType_Reconstructor        *Sintrp;
	struct _MUSCL_PromptType_Rsolver              *Rsolver;
	struct _MUSCL_PromptType_characteristic_solve *Moc;
};
typedef struct _Muscl_Opts Muscl_Opts;


/*
*    The Vec_Muscl structure assumes a system of conservation laws of
*    the form
*
*	U_t + d(A*F)/dV(x) + dH/dx = S
*
*    For gas dynamics the vector U is the vector of conservative
*    quantities U = (rho, E, m) where rho is the mass density,
*    E is the total energy (kinetic + internal) density,  and m is the
*    momentum density uni_array.
*
*    A(x) is the cross sectional area of the section at x.  For example
*    in rectangular geometry A(x) = 1,  while A(x) = x for cylindrical
*    geometry and A(x) = x^2 for spherical geometry.  V is the volume
*    dV = A dx.   S is the body source term.
*    
*    Depending on the solver the definitions of F, H, and S may vary.
*    Note that all of the conservation laws below are equivalent to
*    the one above.
*
*	U_t + d(F+H)/dx = S - F*(dA/dx)/A
*
*	U_t + d(A*(F+H))/dV = S + H*dA/dV = S + H*(dA/dx)/A
*
*    For the oned_MUSCL solver we use
*
*	F = (m, m*m/rho + P, m*E + m*P/rho), H = 0,
*       S = ( -m*(dA/dx)/A,
*             rho*g - (m*m/rho + P)*(dA/dx)/A,
*             <m,g> - (m*E + m*P/rho)*(dA/dx)/A ).
*    
*    For the oned_PLM solver we use
*
*	F = ((m, m*m/rho, m*E + m*P/rho), H = (0,P,0), S = (0,rho*g,<m,g>)
*/

struct _Vec_Muscl_Index {	/* Indices of states stored in u             */
	int  density;
	int  energy;
	int  internal_energy_density;
	int  v[3];
	int  pressure;
	int  sound_speed;
	int  GAM;
	int  FD;
	int  grav;
        int  prho[MAX_NUM_GAS_COMPS];  /* fractional mass index */
};
typedef struct _Vec_Muscl_Index Vec_Muscl_Index;

enum _SWEEP_TYPE {
	UNSET_SWEEP_TYPE = 0,
	REGULAR_INTERIOR_SWEEP,
	IRREGULAR_INTERIOR_SWEEP,
	TANGENTIAL_SWEEP
};
typedef enum _SWEEP_TYPE SWEEP_TYPE;

struct _Vec_Muscl {
	Vec_Gas	*vst;		/* Initial data in vec gas form */
	double   **CellEdge;	/* Coordinates of cell edges    */
	double   **coords;       /* Coordinates of cell centers  */

	double	**u;		/* Initial state data, unconservative form */
	             /* u[index.density] = rho = density                    */
		     /* u[index.energy] = e = specific internal energy      */
		     /* u[index.v[0]] = sweep direction velocity component  */
		     /* u[index.v[i]], 0 < i < dim, velocity components     */
		     /*                ...,    orthogonal                   */
		     /* For PLM the additional values are set               */
		     /* u[index.pressure] = pressure                        */
		     /* u[index.sound_speed] = sound speed                  */
		     /* u[index.GAM] = Gruneisen exponent                   */
		     /* u[index.FD] = Fundamental derivative                */
		     /* u[index.grav] = sweep direction gravity component   */

	int             nvar_u; /* Number of variables stored in the u array */

	Vec_Muscl_Index index;

	double	**ucon;		/* Initial state data, conservative form   */
	                        /* ucon[0] = rho = density                 */
				/* ucon[1] = total energy density =        */
				/*           rho*(0.5*v*v + e)             */
				/* ucon[2] = sweep direction momentum      */
				/*           density                       */
				/* ucon[3],    momentum density components */
				/* ...,        orthogonal                  */
				/* ucon[1+dim] to sweep direction          */

				/* Note: the arrays stored in src and sourc */
				/* can be in either conservative or         */
				/* non-conservative form depending on the   */
				/* the section of the algorithm             */
	Vec_Src	*src;		/* Source terms */
	double	**source;	/* Array form of src */

	Vec_Eigen	Vegn;	/* Eigenuni_array matricies for initial data */

		/* State reconstruction data */
	double	**central;	/* Central differences */
	double	**forward;	/* Forward differences */
	double	**backward;	/* Backward differences */
	double   **dulim;        /* Limited forward and backward differences */
	double   **duf;		/* Van Leer differences */
	double   **sgn;          /* Signs of central differences */
	double	**du;	  	/* limited state slope */

	double	**q;	  	/* state in eigen coordinate */
	double	**dq;	        /* Slopes from linear reconstruction,   */
	                        /* in eigenvalue coordinates            */

		/* Half step data structure */
	Vec_Gas	VLst;	 /* Vec Gas form of left half step state      */
	                 /* or u+c family characteristic trace back   */
	double	**uL;	 /* Half step evolved left state              */
	double	*pL;	 /* Half step evolved left state pressure     */
	double   *gL;     /* Gravity value at left position            */

	Vec_Gas	VMst;	 /* Vec Gas form of mid state from Riemann solution */
	                 /* or contact family characteristic trace back     */
	double	**uM;	 /* Mid state from Riemann solution                 */
	double	*pM;	 /* pressure at mid state                           */
	double   *gM;     /* Gravity value at mid position                   */

	Vec_Gas	VRst;	 /* Vec Gas form of right half step state           */
	                 /* or u-c family characteristic trace back         */
	double	**uR;	 /* Half step evolved right state                   */
	double	*pR;	 /* Half step evolved right state pressure          */
	double   *gR;     /* Gravity value at left position            */

	double	**right_wvs;  /* Contributions from right moving waves */
	double	**left_wvs;   /* Contributions from left moving waves */
	double	**source_sum; /* Source term contributions from 1/2 step */

	MUSCL_FLUX Flux;

	double   **uM0;  /* Cell edge state at start of time step */
	Vec_Gas	VM0st;	/* Vec Gas form of right half step state */
	double	*pM0;	/* pressure at mid state */
	double   *gM0;	/* gravity at mid state */

	double   **uM1;       /* Cell edge state at end of time step   */
	Vec_Gas	VM1st;	/* Vec Gas form of right half step state */
	double	*pM1;	/* pressure at mid state */
	double   *gM1;	/* gravity at mid state */

	MUSCL_FLUX Flux0; /* Fluxes at begining of time step       */
	MUSCL_FLUX Flux1; /* Fluxes at begining of time step       */

	double	**awv;	/* Eigen coordinates of jump at cell boundaries */

	/* Grid description */
	double    *A;      /* A[j] = area of cell face between cells j-1     */
	                  /*        and cell j                              */
	double    *dV;     /* Volume of cell j                               */
	double    alpha;   /* Geometry factor */

	Muscl_Opts Opts;

	Vec_Avisc	*avisc; /* artificial viscosity paramters */
	Vec_MSF		*msf; /* Muscl slope flattening parameters */
	boolean            monotone_reconstruction;
	boolean            link_reconstructions;

	int		offset;	/* Starting index of Vec_Gas arrarys */
	int		vsize;	/* Size of arrays in Vec_Gas structures */
	int		max_vsize;	/* Allocated length of arrays */
	boolean		is_src;	/* YES if source terms are present */
	int		dim;	/* Spatial dimension of flow */
	size_t		sizest;	/* Size of Locstate */
	int		sten_rad;/* Radius of finite difference stencil */
	int		idir;	/* Sweep direction */
	const double* const *Q;	/* rotation vector */
	double		dt;	/* time step */
	double		dn;	/* spatial step */
	double		*grav;	/* array of cell centered projected gravities */

	Front		*front;
	Wave		*wave;

	Stencil         *sten;
	Tan_stencil     *tsten;

	double		***worksp;

	struct {
	    boolean CellEdge;
	    boolean u;
	    boolean ucon;
	    boolean source;
	    boolean central;
	    boolean forward;
	    boolean backward;
	    boolean dulim;
	    boolean duf;
	    boolean sgn;
	    boolean du;
	    boolean q;
	    boolean dq;
	    boolean uL;
	    boolean pL;
	    boolean uM;
	    boolean pM;
	    boolean uR;
	    boolean pR;
	    boolean right_wvs;
	    boolean left_wvs;
	    boolean source_sum;
	    struct {
	        boolean F;
	        boolean H;
	    } Flux;
	    boolean uM0;
	    boolean uM1;
	    struct {
	        boolean F;
	        boolean H;
	    } Flux0;
	    struct {
	        boolean F;
	        boolean H;
	    } Flux1;
	    boolean awv;
	    boolean grav;
	    boolean worksp;
	    boolean A;
	    boolean dV;
	    boolean vmuscl;
	} alloc;

	SWEEP_TYPE _MusclSweepType;
};
typedef struct _Vec_Muscl Vec_Muscl;

#define MusclSweepType(vmuscl) vmuscl->_MusclSweepType

struct _MUSCL_PromptType_Reconstructor {
	const char    *prompt, *select;
	void (*reconstructor)(int,int,Vec_Muscl*);
	void (*half_step)(int,int,double,double,struct _Vec_Muscl*);
	void (*strong_wave_half_step)(int,int,double,double,struct _Vec_Muscl*);
};
typedef struct _MUSCL_PromptType_Reconstructor MUSCL_PromptType_Reconstructor;

struct _MUSCL_PromptType_Rsolver {
	const char *prompt, *select;
	void (*rsolver)(int,int,double**,Vec_Gas*,double**,Vec_Gas*,
                        double**,Vec_Gas*,MUSCL_FLUX*,Vec_Muscl*);
	void (*rmidstate)(int,int,Vec_Gas*,Vec_Gas*,
	                  double*,double*,Vec_Muscl*);
};
typedef struct _MUSCL_PromptType_Rsolver MUSCL_PromptType_Rsolver;

struct _MUSCL_PromptType_characteristic_solve {
	const char    *prompt, *select;
	void (*characteristic_solve)(int,int,struct _Vec_Muscl*);
};
typedef struct _MUSCL_PromptType_characteristic_solve MUSCL_PromptType_characteristic_solve;

/* Access macros for Muscl_Opts functions */
	
#define compute_eigens(start,end,vmuscl)				\
    if ((vmuscl)->Opts._compute_eigens)					\
        (*(vmuscl)->Opts._compute_eigens)(start,end,vmuscl)

#define compute_art_visc_coefs(start,end,vmuscl)			\
    if ((vmuscl)->Opts._compute_art_visc_coefs)				\
        (*(vmuscl)->Opts._compute_art_visc_coefs)(start,end,vmuscl)

#define reconstructor(start,end,vmuscl)					\
    (*(vmuscl)->Opts._reconstructor)(start,end,vmuscl)

#define half_step(start,end,dt,dn,vmuscl)				\
    (*(vmuscl)->Opts._half_step)(start,end,dt,dn,vmuscl)

#define rsolver(start,end,uL,vlst,uR,vrst,uM,vmst,flux,vmuscl)		\
    (*(vmuscl)->Opts._rsolver)(start,end,uL,vlst,uR,vrst,uM,vmst,flux,vmuscl)

#define rmidstate(start,end,vlst,vrst,pm,vm,vmuscl)			\
    (*(vmuscl)->Opts._rmidstate)(start,end,vlst,vrst,pm,vm,vmuscl)

#define characteristic_solve(start,end,vmuscl)				\
    (*(vmuscl)->Opts._characteristic_solve)(start,end,vmuscl)

#define flux_vectors(start,end,uM,p,Flux,vmuscl)			\
    (*(vmuscl)->Opts._flux_vectors)(start,end,uM,p,Flux,vmuscl)

#define add_art_visc1(start,end,uM,vmuscl)				\
    if ((vmuscl)->Opts._add_art_visc1)					\
        (*(vmuscl)->Opts._add_art_visc1)(start,end,uM,vmuscl)

#define add_art_visc2(start,end,vmuscl)					\
    if ((vmuscl)->Opts._add_art_visc2)					\
        (*(vmuscl)->Opts._add_art_visc2)(start,end,vmuscl)

#define cons_src(start,end,swp_num,iperm,tsten,vmuscl)		\
    (*(vmuscl)->Opts._cons_src)(start,end,swp_num,iperm,tsten,vmuscl)

#define print_internal_energy(mesg,ucon,vmuscl,start,end)		\
    (*(vmuscl)->Opts._print_internal_energy)(mesg,ucon,vmuscl,start,end)

#define	Cg_params(mopts)	((mopts)->_Cg_params)

#include <ghyp/ghypprotos.h>

#endif /* !defined(_GHYP_H) */
