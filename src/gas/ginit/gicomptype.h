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
*				gicomptype.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the structures for comp_type initialization.
*/


#if !defined(_GICOMPTYPE_H)
#define	_GICOMPTYPE_H

#include <ginit/ginit.h>

		/* possible component types */

enum _COMP_TYPE_TYPE {
	UNSET_COMP_TYPE          =  0,
	EXTERIOR,
	OBSTACLE,
	AMBIENT,
	ELLIPTICAL,
	KH_SINE_PERTURBED,      	/* obsolete */
	NWAVE,
	RANDOM_SURFACE_PERTURBED,	/* obsolete */
	BUBBLE,
	PRANDTL_MEYER_WAVE,
	UNTRACKED_SHOCK,
	RESTART,
	TAYLOR_WAVE,
	RT_KH,
	RT_PERTURBED,
	KH_PERTURBED,
	RAREFACTION_WAVE_1D,
	TRANS_LAYER,
	STRETCHING,
	ONE_DIMENSIONAL_OVERLAY,
	RANDOM_REGION,
	ONE_DIMENSIONAL_TABLE,
	TABULATED_REGION
};
typedef enum _COMP_TYPE_TYPE COMP_TYPE_TYPE;

struct _COMP_TYPE {
	COMP_TYPE_TYPE 	type;
	COMPONENT	comp;
	Gas_param 	*params;
	void 		(*_get_state)(double*,Locstate,struct _COMP_TYPE*,
				      HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
	void		(*free_comp_type_extra)(struct _COMP_TYPE*);
	POINTER 	extra;
};
typedef struct _COMP_TYPE COMP_TYPE;


	/* thermodynamic properties */

enum _STRATIFICATION_TYPE {
    CONSTANT	     = 0,
    ISOTHERMAL	     = 1,
    ADIABATIC	     = 2,
    CONSTANT_DENSITY = 3,
    HYDRO_STATIC     = 4
};
typedef enum _STRATIFICATION_TYPE STRATIFICATION_TYPE;


	/* macros for the comp_type extra field */

#define	Get_state(coords,state,comp_type,hs,intfc,init,type)		\
	(*((comp_type)->_get_state))(coords,state,comp_type,hs,intfc,init,type)
#define	Get_tgas_state(coords,state,comp_type,intfc,init)		\
	Get_state(coords,state,comp_type,NULL,intfc,init,TGAS_STATE)

#define	Ambient(ct)	          ((Locstate)(ct)->extra)
#define	Untracked_shock(ct)	  ((UT_SHOCK *)	 (ct)->extra)
#define	Taylor_wave(ct)	          ((_TAYLOR_WAVE *)(ct)->extra)
#define	Elliptical(ct)	          ((_ELLIPTICAL *)(ct)->extra)
#define	Wall_edge_equilibrium(ct) ((_WAll_EDGE_EQUILIBRIUM *)(ct)->extra)
#define	Rt_kh(ct)	          ((_RT_KH *)(ct)->extra)
#define	Rt_perturbed(ct)	  ((_RT_PERTURBED *)(ct)->extra)
#define	Rarefaction_wave_1d(ct)	  ((_RAREFACTION_WAVE_1D *) (ct)->extra)
#define	Trans_layer(ct)           ((_TRANS_LAYER *)(ct)->extra)
#define	Stretching(ct)	          ((_STRETCHING *)(ct)->extra)
#define Riemann_2d(ct)  	  ((_RIEMANN_2D *)(ct)->extra )
#define	One_d_overlay(ct)	  ((ONED_OVERLAY *)(ct)->extra)
#define	Random_state(ct)	  ((RANDOM_STATE *)(ct)->extra)
#define	Sine_pert(ct)	          ((SINE_PERT *)(ct)->extra)
#define Tabulated_region_data(ct) ((TABULATED_REGION_DATA*)(ct)->extra)
#define One_d_table(comp_type)    ((ONED_TABLE *)(comp_type)->extra )

struct _ONED_TABLE{
        int             nx;
        double           *d;
        double           *p;
        double           *v;
        double           *pt;
};
typedef struct  _ONED_TABLE ONED_TABLE;

struct _TABULATED_REGION_DATA {
	double     GL[3], GU[3];
	double     h[3];
	double     *rho;
	double     *p;
	double     *v[3];
	int       n[3], len;
	boolean      tab_region_set;
	COMP_TYPE *allocated_from;
};
typedef struct _TABULATED_REGION_DATA TABULATED_REGION_DATA;

enum _OVERLAY_TYPE {
	OVERLAY_TYPE_UNSET  = -1,
	RADIAL_OVERLAY      =  0,
	CYLINDRICAL_OVERLAY,
	RECTANGULAR_OVERLAY
};
typedef enum _OVERLAY_TYPE OVERLAY_TYPE;

struct	_ONED_OVERLAY {
	OVERLAY_TYPE	overlay_type;
	double		origin[3],	/* The 1D data is extended by */
			direction[3];	/* by orthogonal extension */
					/* from the directed line with   */
					/* the given origin and direction*/
	INTERFACE	*intfc1d;
	Front		*front;
	Printplot       *prt;
	INPUT_SOLN	**is;		/* Stores oned data */
};
typedef	struct	_ONED_OVERLAY ONED_OVERLAY;

		/* Initialization info for untracked shocks */

struct _UT_SHOCK {
	Locstate	state0, state1;	/* state0 ahead, state1 behind */
	double		nor[MAXD];	/* high to low pressure normal */
	double		posn[MAXD];	/* Given point on shock wave   */
	double		width;		/* thickness of wave           */
	int		_wave_type;
	COMP_TYPE	*ctype0, *ctype1;
	boolean		free_with_comp_type;
};
typedef struct _UT_SHOCK UT_SHOCK;

struct __TAYLOR_WAVE {
	Locstate st0, st1, stas;
	double z0, z1, zs;
	double sgn;
	double vz0, vz1, va[MAXD];
	double zbar, tbar;
	size_t sizest;
	int l_or_r;
};		/* obsolete */
typedef struct __TAYLOR_WAVE _TAYLOR_WAVE;

struct __RAREFACTION_WAVE_1D {
	WAVE_FAMILY	l_or_r;	      /* left or right family */
	double		zbar, tbar;   /* center of wave */
	double		zl,  zt;      /* positions of edges */
	Locstate	stl, stt;     /* TGAS_STATE at l & t edges */
	double           zmin, zmax;   /* range of wave */
	double           spl, spt;     /* wave speed at edges */
	Locstate        stmin, stmax; /* states at edges of wave range*/
	ELLIPSOID	*el_lead;     /* ellipsoid at leading edge */
	ELLIPSOID	*el_trail;    /* ellipsoid at trailing edge */
	LAYER_SURF      *lead;        /* layer surface at leading edge */
	LAYER_SURF      *trail;       /* layer surface at trailing edge */
};			/* 1d centered rarefaction wave */
typedef struct __RAREFACTION_WAVE_1D _RAREFACTION_WAVE_1D;

struct __ELLIPTICAL {
	ELLIPSOID            *ellipsoid;
	RANDOM_STATE	     *rstate;
	Locstate	     state;
	Locstate	     wkstate[2];
	double                weight[MAXD];
	double		     r0;
	_RAREFACTION_WAVE_1D *rw1d;
	STRATIFICATION_TYPE  stratification_type;
};
typedef struct __ELLIPTICAL _ELLIPTICAL;
#define	RadialVelocity(ellip)	Vel(ellip->state)[0]
						 
struct __RT_KH {
	STRATIFICATION_TYPE stratification_type;
	double		    ref_coords[MAXD];
	Locstate	    ref_state;
};
typedef struct __RT_KH _RT_KH;

struct __RT_PERTURBED {
	int		layer_label;
	LAYER_SURF      *lower_surf,    *upper_surf;
	_RT_KH		*rt_kh;
	int		num_modes;
	int		lin_pert_intvl;
	NORMAL_MODE	**normal_mode;
};
typedef struct __RT_PERTURBED _RT_PERTURBED;

struct __TRANS_LAYER {
	LAYER_SURF      *lower_surf,    *upper_surf;
	Locstate	lower_st,	upper_st;
};
typedef struct __TRANS_LAYER _TRANS_LAYER;

struct __STRETCHING {
	double		v[8][3];
	double		L[3], U[3];
	Locstate	ambient;
};
typedef struct __STRETCHING _STRETCHING;

struct __RIEMANN_2D {
	int		num_sections;
	double		*start_angle;
	Locstate	*states;
};
typedef struct __RIEMANN_2D _RIEMANN_2D;

#endif /* !defined(_GICOMPTYPE_H) */
