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
*				guserint.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*			User Supplied Structure Components
*/

#if !defined(_GUSERINT_H)
#define _GUSERINT_H

#include <gdecs/gstate.h>

		/* structure used for interacting tracked waves */

enum {MAX_N_CURVES = 7 };

struct _RP_DATA {
	ANGLE_DIRECTION	ang_dir;
	double		ang[MAX_N_CURVES];
	Locstate	state[MAX_N_CURVES];
	int		stype;
	double		M[MAX_N_CURVES];
	double		theta[MAX_N_CURVES];
	boolean		intfc_table_storage;
};
typedef struct _RP_DATA RP_DATA;

/*	Gas extenstions to Interface data structures */

struct _G_HYPER_SURF {
	F_HYPER_SURF	f_hs;
	int             _layer_index;
	boolean            _no_slip;
	double		_adherence_coeff;
	double           _create_time;
};
typedef struct _G_HYPER_SURF G_HYPER_SURF;

	/* G_HYPER_SURF access macros */
#define g_hyper_surf(hs)        ((G_HYPER_SURF *) Hyper_surf(hs))
#define	layer_index(hs)		(g_hyper_surf(hs)->_layer_index)
#define	no_slip(hs)		(g_hyper_surf(hs)->_no_slip)
#define	adherence_coeff(hs)	(g_hyper_surf(hs)->_adherence_coeff)
#define create_time(hs)         (g_hyper_surf(hs)->_create_time)

struct _G_CURVE {
	F_CURVE	fcurve;
};
typedef struct _G_CURVE G_CURVE;

	/* G_CURVE access macros */
#define g_curve(curve)                  ((G_CURVE *) (curve))

typedef enum { DONT_ADJUST_ANGLE, ADJUST_ANGLE} ADJUST_ANGLE_VALUE;

struct _G_NODE {
	F_NODE	fnode;
#if defined(TWOD)
	RP_DATA *RP;
	struct _ADJUST_ANGLE_DATA {
	    double	_adjust_len;
	    ADJUST_ANGLE_VALUE (*_adjust_angle_len)(NODE*,O_CURVE*,O_CURVE*,
						    Locstate,Locstate,
						    double*,double,double*,
						    double*,Front*,Wave*);
	} _adjust_angle;
#endif /* defined(TWOD) */
};
typedef struct _G_NODE G_NODE;


	/* G_NODE access macros */
#define g_node(node)	((G_NODE *) (node))

#if defined(TWOD)

#if defined(__cplusplus)
typedef struct _G_NODE::_ADJUST_ANGLE_DATA ADJUST_ANGLE_DATA;
#else /* defined(__cplusplus) */
typedef struct _ADJUST_ANGLE_DATA ADJUST_ANGLE_DATA;
#endif /* defined(__cplusplus) */

#define Rp_data(node)	(g_node(node)->RP)
#define	adjust_angle(node)	g_node(node)->_adjust_angle
#define adjust_len(node)	adjust_angle(node)._adjust_len
#define adjust_angle_len(newn,oldc,newc,st_l,st_r,dir,dt,adjust_len,sonic_rad,fr,wave) \
	(*adjust_angle(newn)._adjust_angle_len)(newn,oldc,newc,   \
						st_l,st_r,dir,dt, \
					   	adjust_len,sonic_rad,fr,wave)

#endif /* defined(TWOD) */


struct _WSSten {
	double		   *coords;
	double		   **lcrds;/*Position of left states at start of step*/
	double		   **rcrds;/*Position of right states at start of step*/
	double              *V;     /*Velocity of moving frame for point */
	Locstate	   *sl;	   /*Array of left states*/
	Locstate	   *tsl;   /*Array of left states same types as slopes*/
	Locstate	   *dsl;   /*Left state slopes*/
	Locstate	   *sr;	   /*Array of right states*/
	Locstate	   *dsr;   /*Right state slopes*/
	Locstate	   *tsr;   /*Array of left states same types as slopes*/
	int		   nsts;   /*Number of states on left and right*/
	COMPONENT	   ncomp;
	COMPONENT	   pcomp;
	int		   w_type;
	int		   stype;
	POINT		   *p;     /*Point being propagated */
	HYPER_SURF	   *hs;	   /*Hypersurface being propagated*/
	HYPER_SURF_ELEMENT *hse;   /*Hypersurface element containing p*/
	double		   *nor;   /*Normal vector to front*/
	double		   dn;	   /*Spatial grid spacing in normal direction*/
	double		   dt;	   /*Time step*/
	double		   pjump;  /*Pressure jump allowed across a contact*/
	Front		   *front;
	Wave		   *wave;

		/*Function hook for interpolation along stencil*/
	void		   (*_ws_interpolate)(Locstate,double,
		                              SIDE,int,struct _WSSten*);
		/*Function hook to set slopes for interpolation*/
	void		   (*_set_ws_slopes)(struct _WSSten*);

		/*Clear an existing WSSten structure for reuse*/
	void		   (*_ClearWSStenData)(struct _WSSten*);

		/*Print the stencil data*/
	void		   (*_FPrintWSSten)(FILE*,struct _WSSten*);

	/* Internal storage for states and positions */
	byte		   *sl_store, *dsl_store, *tsl_store;
	byte		   *sr_store, *dsr_store, *tsr_store;
	double		   **lcrds_store, **rcrds_store;
	double		   coords_store[3];
	double		   nor_store[3];
};
typedef struct _WSSten WSSten;

#define	ws_interpolate(ans,x,side,stype,wssten)    			\
	(*(wssten)->_ws_interpolate)(ans,x,side,stype,wssten)

#define	set_ws_slopes(wssten)					\
	(*(wssten)->_set_ws_slopes)(wssten)

#define	ClearWSStenData(wssten)	(*(wssten)->_ClearWSStenData)(wssten)

#define	FPrintWSSten(file,wssten)	(*(wssten)->_FPrintWSSten)(file,wssten)

#define	PrintWSSten(wssten)	FPrintWSSten(stdout,wssten)

struct	_USWSSten2d {
	Tan_stencil	**_tan_fr;	/* Data for tang. updates */
	WSSten		**_nor_fr;	/* Data for normal updates */
	int		nor_rad;	/* normal radius of 2d stencil */
	int		tan_rad;	/* tangential radius of 2d stencil */
	Tan_stencil	*nor_ans;	/* Ans. from normal updates */
	WSSten		*tan_ans;	/* Ans. from tang. updates */
	double		ds, dn, dt;
	double		tngt[MAXD];
	double           **_nor_vec;	/* Normal vector for each WSSten */
	Front		*fr;
	Wave		*wave;
	HYPER_SURF	*hs;
};
typedef struct  _USWSSten2d USWSSten2d;

#define	tan_fr(sten)	(sten)->_tan_fr
#define	nor_fr(sten)	(sten)->_nor_fr
#define nor_vec(sten)	(sten)->_nor_vec

struct _CWNP {
	double   wall_bond_len;
	double   first_adjust_time;
	int     first_adjust_step;
	boolean adjust;
};
typedef struct _CWNP CWNP;

struct _G_USER_INTERFACE {
	int        _intfc_type;
	int        num_params;
	Gas_param  **params_list;
	void       (*_stratified_state)(Locstate,double,double,Locstate);
	const char *stratified_state_name;
	void	   (*_w_speed)(double*,Locstate,Locstate,Locstate,Locstate,
			       double*,double,double*,int,Front*);
	void	   (*_npt_w_speed)(WSSten*,Locstate,Locstate,double*);
	void	   (*_unsplit_w_speed2d)(USWSSten2d*,Locstate,Locstate,double*);
	CWNP       _ContactWallNodeParams;
};
typedef struct _G_USER_INTERFACE G_USER_INTERFACE;

struct _G_INTERFACE {
	F_INTERFACE f_intfc;
	G_USER_INTERFACE g_user_intfc;
	int	_num_layers;
};
typedef struct _G_INTERFACE G_INTERFACE;

	/* G_INTERFACE access macros */
#define g_interface(intfc)	((G_INTERFACE *) (intfc))
#define g_user_interface(intfc)	(g_interface(intfc)->g_user_intfc)
#define interface_type(intfc)	(g_user_interface(intfc)._intfc_type)
#define gas_params_list(intfc)	(g_user_interface(intfc).params_list)
#define num_gas_params(intfc)	(g_user_interface(intfc).num_params)
#define	set_stratified_state(f)	g_set_stratified_state(f,#f)
#define	num_layers(intfc)	(g_interface(intfc)->_num_layers)
#define	ContactWallNodeParams(intfc)					\
    (g_user_interface(intfc)._ContactWallNodeParams)
#define	contact_wall_node_params(intfc)					\
    &ContactWallNodeParams(intfc)


/*	End of Gas extenstions to Interface data structures */

/* 
*	Possible types of interfaces in gas dynamics simulations.
*	The default value set in make_interface is PHYSICAL_INTERFACE
*	which is ordinary interface in the physical x-y space.
*	Other types must be set directly.  These include EOS_INTERFACE
*	which correspond to interfaces in the thermodyamic phase spaces.
*/

enum _G_INTERFACE_TYPE {
	PHYSICAL_INTERFACE = 1,
	EOS_INTERFACE
};
typedef enum _G_INTERFACE_TYPE G_INTERFACE_TYPE;


	/* possible values for the wave_type of a CURVE */

enum {
	TIME_DIRICHLET_BOUNDARY	     = FIRST_PHYSICS_WAVE_TYPE,
	VELOCITY_SPECIFIED,
	NO_SLIP_NEUMANN_BOUNDARY,
	FIRST_BOUNDARY_EOS_WAVE_TYPE,
	CONTACT	                     = FIRST_SCALAR_PHYSICS_WAVE_TYPE,
	THIN_FLAME,
	FIRST_INTERIOR_EOS_WAVE_TYPE,
	BACKWARD_SHOCK_WAVE          = FIRST_VECTOR_PHYSICS_WAVE_TYPE,
	FORWARD_SHOCK_WAVE,
	FORWARD_SOUND_WAVE_LE,
	FORWARD_SOUND_WAVE_TE,
	BACKWARD_SOUND_WAVE_LE,
	BACKWARD_SOUND_WAVE_TE,
	MARKED_CURVE                  = (FIRST_VECTOR_PHYSICS_WAVE_TYPE + 100),
	RIEMANN_PROBLEM_WAVE
};

enum _RAREFACTION_EDGE_TYPE {
	LEADING_EDGE,
	TRAILING_EDGE
};
typedef enum _RAREFACTION_EDGE_TYPE RAREFACTION_EDGE_TYPE;

#define is_scalar_wave(w_type)	((w_type) == CONTACT)
#define is_thinflame_wave(w_type)	((w_type) == THIN_FLAME)

#define is_vector_wave(w_type)						\
(		((w_type) ==     FORWARD_SHOCK_WAVE ) ||		\
		((w_type) ==    BACKWARD_SHOCK_WAVE ) ||		\
		((w_type) ==  FORWARD_SOUND_WAVE_LE ) ||		\
		((w_type) == BACKWARD_SOUND_WAVE_LE ) ||		\
		((w_type) ==  FORWARD_SOUND_WAVE_TE ) ||		\
		((w_type) == BACKWARD_SOUND_WAVE_TE )			)

#define is_shock_wave(w_type)						\
	(((w_type) == FORWARD_SHOCK_WAVE) || ((w_type) == BACKWARD_SHOCK_WAVE))

#define is_rarefaction_wave(w_type)					\
	(   ((w_type) == FORWARD_SOUND_WAVE_LE) ||			\
		((w_type) == BACKWARD_SOUND_WAVE_LE) ||			\
		((w_type) == FORWARD_SOUND_WAVE_TE) ||			\
		((w_type) == BACKWARD_SOUND_WAVE_TE)    )

#define is_rarefaction_leading_edge(w_type)				\
	((w_type)==FORWARD_SOUND_WAVE_LE || (w_type)==BACKWARD_SOUND_WAVE_LE)

#define is_rarefaction_trailing_edge(w_type)				\
	((w_type)==FORWARD_SOUND_WAVE_TE || (w_type)==BACKWARD_SOUND_WAVE_TE)

#define is_forward_wave(w_type)						\
	(((w_type) == FORWARD_SHOCK_WAVE) || 				\
	((w_type) == FORWARD_SOUND_WAVE_LE) || 				\
	((w_type) == FORWARD_SOUND_WAVE_TE))

#define is_backward_wave(w_type)					\
	(((w_type) == BACKWARD_SHOCK_WAVE) || 				\
	((w_type) == BACKWARD_SOUND_WAVE_LE) || 			\
	((w_type) == BACKWARD_SOUND_WAVE_TE))

	/* possible values for the class of waves */
/*
* IMPORTANT NOTE:  the values below are used as an index in an array 
* to look up user supplied information on whether to track waves at 
* a node when it is produced by a dynamic bifurcation.  Thus the particular 
* values of each symbolic name is important and should not be changed.
*/

enum {
	CONTACT_WAVE	    = 0,
	SHOCK_WAVE,
	RAREF_LEADING_EDGE,
	RAREF_TRAILING_EDGE,
	NUM_WAVE_CLASSES
};

	/* possible values for the start/end_status of a CURVE */
/*
* IMPORTANT NOTE:  the value of the status at node below
* is used as an index in an array to look up user supplied information
* on whether to track waves at a node when it is produced by a dynamic 
* bifurcation.  Thus the particular values of each status is important 
* and should not be changed.
*/

enum {
	MACH_STEM =  	FIRST_PHYSICS_CURVE_STATUS,
	SLIP,
	TRANSMITTED,
	OVERTOOK,
	CONTACT_TARGET,
	ONED_WAVE,
	NUM_NODE_STATUS
};

enum {

	/* possible values for the type of a _HYPER_SURF_BDRY */

/*
* IMPORTANT NOTE:  the value node_type - FIRST_PHYSICS_HSBDRY_TYPE
* is used as an index in an array to look up user supplied information
* on whether to track waves at a node of the given type when it is
* produced by a dynamic bifurcation.  Thus the particular values of
* each node type are important and should not be changed.
*/

	B_REFLECT_HSBDRY  = FIRST_PHYSICS_HSBDRY_TYPE,
	ATTACHED_B_HSBDRY,
	MACH_HSBDRY,
	CROSS_HSBDRY,
	OVERTAKE_HSBDRY,
	DIFFRACTION_HSBDRY,
	TRANSMISSION_HSBDRY,
	CC_HSBDRY,
	WAVE_END_HSBDRY,
	TOT_INT_REFL_HSBDRY,

	/* possible values for the type of a NODE */
/*
* IMPORTANT NOTE:  the value node_type - FIRST_PHYSICS_NODE_TYPE
* is used as an index in an array to look up user supplied information
* on whether to track waves at a node of the given type when it is
* produced by a dynamic bifurcation.  Thus the particular values of
* each node type are important and should not be changed.
*/

	B_REFLECT_NODE	  = B_REFLECT_HSBDRY,
				/* 1 incident, 1 reflected, 2 boundaries */
	ATTACHED_B_NODE	  = ATTACHED_B_HSBDRY,
				/* 1 shock or contact, >= 2 boundaries */
	MACH_NODE	  = MACH_HSBDRY,
				/* 1 incident, 1 reflected, */
				/* 1 contact, 1 Mach */
	CROSS_NODE	  = CROSS_HSBDRY,
				/* shock-shock: 4 shocks, 1 contact */
	OVERTAKE_NODE	  = OVERTAKE_HSBDRY,
				/* shock-shock: 3 shocks, 1 contact */
	DIFFRACTION_NODE  = DIFFRACTION_HSBDRY,
				/* shock-contact: 3 shocks, 2 contacts */
	TRANSMISSION_NODE = TRANSMISSION_HSBDRY,
				/* shock-contact: 2 shocks, 2 contacts */
	CC_NODE		  = CC_HSBDRY,
				/* contact-contact */
	WAVE_END_NODE	  = WAVE_END_HSBDRY,
				/* end of tracked wave */
	TOT_INT_REFL_NODE  = TOT_INT_REFL_HSBDRY,
				/* total internal reflection node */

/* The other "node types" are used only for the node untrack flags. */
	PRECURSOR_RR_DIFFRACTION, /* cluster of nodes, transmission, cross,
				   * total internal reflection, and possibly
				   * overtake nodes.
				   */
	ONED_INTERACTION,         /*one dimensional "node"*/
	NUM_PHYS_NODE_TYPES	= ONED_INTERACTION - B_REFLECT_NODE + 1,

	FIRST_EOS_NODE_TYPE = FIRST_PHYSICS_NODE_TYPE + 100,
        /* possible values for the type of a BOUNDARY CURVE */
        ATTACHED_B_CURVE          = ATTACHED_B_HSBDRY
};

enum {
	ST_PARAMS_ID = FIRST_PHYSICS_MESSAGE_ID,
	NCOMPS_ID,
	COMPS_ID,
	RM_DATA_ID
};

#endif /* !defined(_GUSERINT_H) */
