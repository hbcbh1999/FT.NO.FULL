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
*				ginit.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the structures for gas dynamics initialization.
*/

#if !defined(_GINIT_H)
#define _GINIT_H

#include <gdecs/gdecs.h>

	/* types of ramp problems */
enum {
	NORMAL	    = 1,
	OBLIQUE,
	DOUBLE_REFL,
	MACH_ANGLES
};

	/* component ft_assignments */

enum {
	COMPOBST = MIN_INTERIOR_COMP,
	FIRST_DYNAMIC_COMPONENT,

	COMPL	      = FIRST_DYNAMIC_COMPONENT,
	COMPR	      = FIRST_DYNAMIC_COMPONENT + 1,

	COMPA	      = FIRST_DYNAMIC_COMPONENT,
	COMPB,
	COMPC,
	COMPMID	      = COMPC,
	COMPBOW	      = COMPC,
	COMPMACH,
	COMPU,
	COMPD,

	COMPAL	      = FIRST_DYNAMIC_COMPONENT,
	COMPAR,
	COMPBL,
	COMPBR,

	COMPBS	      = FIRST_DYNAMIC_COMPONENT,
	COMPAS,
	COMPBC,
	COMPIP,
	COMPBP,

		/* Components for diffraction node */

	COMPIF	      = FIRST_DYNAMIC_COMPONENT,
	COMPIB,
	COMPRM,
	COMPRB,
	COMPTB,
	COMPTF,
	COMPBB	      = FIRST_DYNAMIC_COMPONENT + 7
};

                 /* various properties of a multilayer system */

struct _LAYER_FLAG {
	boolean _MOVING_FRAME;
	boolean _ALL_CONTACT;
	boolean _HAS_ELLIPSOID;
	boolean _ELLIPTICAL_REGION;
	boolean _RT_LINEAR_THEORY;
	boolean _INNERMOST_ELLIPSOID;
	boolean _ALLOW_SECTORED_ELLIPSOIDS;
};
typedef struct _LAYER_FLAG LAYER_FLAG;

#define CLEAR_LAYER_FLAG(FLAG)	        zero_scalar(&FLAG,sizeof(LAYER_FLAG));
#define clear_layer_flag(flag)	        zero_scalar(flag,sizeof(LAYER_FLAG));
#define	MOVING_FRAME(FLAG)	        (FLAG)._MOVING_FRAME
#define	moving_frame(flag)	        (flag)->_MOVING_FRAME
#define	ALL_CONTACT(FLAG)	        (FLAG)._ALL_CONTACT
#define	all_contact(flag)	        (flag)->_ALL_CONTACT
#define	HAS_ELLIPSOID(FLAG)	        (FLAG)._HAS_ELLIPSOID
#define	has_ellipsoid(flag)	        (flag)->_HAS_ELLIPSOID
#define	ELLIPTICAL_REGION(FLAG)	        (FLAG)._ELLIPTICAL_REGION
#define	elliptical_region(flag)         (flag)->_ELLIPTICAL_REGION
#define	RT_LINEAR_THEORY(FLAG)	        (FLAG)._RT_LINEAR_THEORY
#define	rt_linear_theory(flag)	        (flag)->_RT_LINEAR_THEORY
#define	INNERMOST_ELLIPSOID(FLAG)       (FLAG)._INNERMOST_ELLIPSOID
#define	innermost_ellipsoid(flag)       (flag)->_INNERMOST_ELLIPSOID
#define	ALLOW_SECTORED_ELLIPSOIDS(FLAG) (FLAG)._ALLOW_SECTORED_ELLIPSOIDS
#define	allow_sectored_ellipsoids(flag) (flag)->_ALLOW_SECTORED_ELLIPSOIDS
 

		/* these symbols are used in range checking */

enum {
	GE_AND_LE = 1,
	LE_OR_GE,
	GE_,
	_LE
};
 

		/* possible types of sines for fronts */

enum {
	SINE		  = 1,	/* horizontal sine front */
	SINE_IN_VORTICITY	/* for KH */
};


typedef	struct	{ double	nu,epsilon; }	SIN_SQR_PERT_PARAMS;


typedef struct {
	double y0, y1;
	Locstate st0, st1;
} Lin_comb_comp;

typedef struct {
	double	      center_of_mass[3];/* Center of mass for rigid body */
	double	      total_mass;	/* Total mass of rigid body */
	double	      mom_of_inertial;	/* Moment of inertial of rigid body */
	boolean	      rotation_only;	/* Yes if for rotation only */
	boolean	      vertical_motion_only;		/* fixing the x axis */
} RIGID_BODY_PARAMS;

/*
*			Ellipsoid
*
*	Describes a possibly perturbed ellipsoid
*
*	A point on the ellipsoid is given by
*
*	p[i] = Q[i][j]*D[j][k]*r[k] (assuming summing convention)
*
*	where
*			Q is a rotation bi_array,
*			D[j][k] = delta[j][k]*(rad[j] + pert(r)[j])
*			r is a unit vector in S^(dim-1)
*	and
*			pert(r) is a Fourier polynomial in r.
*	The possible null data structre fpoly contains the Fourier
*	coefficients needed to compute pert(r).  If fpoly = NULL
*	then pert(r) is identically zero and we have a unperturbed
*	ellipsoid.
*
*/

struct _ELLIPSOID {
	double 	      cen[3];		/* Center of Ellipsoid */
	double 	      rad[3];		/* Lengths of radii */
	double 	      **Q;		/* Rotation matrix */
	double	      *Qrows[3],
		      Qstore[3][3];
	double 	      ThetaS[2],
		      ThetaE[2];	/* Spherical coords of start and end */
	double         scale;		/* Scale radius */
	boolean          closed;		/* Closed ellipsoid if YES */
	ORIENTATION   nor_orient;	/* Specifies inward or outward normal */
	COMPONENT     compin, compout;  /* left and right component */
	FOURIER_POLY  *fpoly;		/* Fourier Perturbation factors */
	LEGENDRE_POLY *lpoly;		/* Legendre Perturbation factors */
	HYPER_SURF    *hs;		/* Hypersurface of ellipsoid */
	double 	      surf_tension;	/* Surface tension */
	RIGID_BODY_PARAMS rgb_params;
	int 	      wv_type;		/* Wave type of ellipsoid */
	int 	      dim;		/* Dimension of embedding space */
	boolean	      untracked;	/* a flag for tracking the ellip */
	int	      layer_index;	/* Identifies specific layers */ 
	boolean          n2o_enforced;
	struct _ELLIPSOID *_inner, *_outer;
	HYPER_SURF    *(*_make_ellipsoid)(struct _ELLIPSOID*,
					  COMPONENT,COMPONENT,Front*);
	void          (*_make_ellip_region_boundaries)(struct _ELLIPSOID*,
						       Front*);
	int           btype[4];
	COMPONENT     bcomp[4];
	boolean          rbdry[4];
	int           obtype[4];
	double 	      odir[4][3];
	COMPONENT     obcomp[4];
	boolean          reset_position;
	double         vr[3];
	struct _ELLIPSOID *rpfronts[7];
};
typedef struct _ELLIPSOID ELLIPSOID;

#define	make_ellipsoid(ellip,compin,compout,front)			\
	(*(ellip)->_make_ellipsoid)(ellip,compin,compout,front)

#define make_ellip_region_boundaries(ellip,front)			\
	(*(ellip)->_make_ellip_region_boundaries)(ellip,front)

#define inner_ellipsoid(ellip)	(ellip)->_inner
#define outer_ellipsoid(ellip)	(ellip)->_outer

#if defined(COMBUSTION_CODE)
		/* initialization info for reacting flows */

typedef	struct {
	int num_of_steps;
	double reaction_zone_width,x_intfc,x_bndry,wv_speed;
	double progress_data[300],step;
	Locstate initial_state;	/* TGas state behind front */
} Reaction;
#endif /* defined(COMBUSTION_CODE) */

struct _LAYER_SURF {
	COMPONENT		l_comp, r_comp;
	double			pbar[3];
	double			s_max, s_min;
	double			nor[3];
	FOURIER_POLY		*fpoly;
	int			wv_type;
	int			dim;
	double			surf_ten;
	int			layer_index;
	boolean                    created;
	boolean                    reset_position;
	double                   velocity[3];
	struct _UT_SHOCK	*untracked;
};
typedef struct _LAYER_SURF LAYER_SURF;
	
struct _LAYER {
	int		layer_label;
	COMPONENT	comp;
	int		num_ellips;
	ELLIPSOID	**ellip;
	LAYER_SURF	*lower_surf,    *upper_surf;
	struct	_LAYER	*next, *prev;
};
typedef	struct _LAYER LAYER;

typedef	struct {
	int		num_layers;
	LAYER		**layer;
	Front		*front;
	double           dt;
	LAYER_FLAG	flag;
} LAYER_SYS;

typedef struct {
	double		*wv_num, phase, amp_max;
	double		sigma_r, sigma_i;
	double		ksq, **a, **b;
} NORMAL_MODE;

typedef struct {
	int		intfc_label, tot_pts;		/* for RT */
	NORMAL_MODE	*normal_mode;
	LAYER_SYS	*layer_sys;
	POINTER		pointer;
	INIT_DATA	*init;
	double		g_z;
} LIN_PERT;

typedef struct {
	Prompt_type	ptype;
	void		(*ppsd)(INIT_DATA*,INIT_PHYSICS*);
} Prob_type;

struct _Scalar_Plot_choice {
	Plot_choice	Choice;
	ScalarPlotItem	Item;
	boolean	use_for_dim[4];
};
typedef	struct _Scalar_Plot_choice Scalar_Plot_choice;

#define	scalar_plot_choice(pc)	((Scalar_Plot_choice*)(pc))
#define	scalar_plot_item(pc)	(&scalar_plot_choice(pc)->Item)

#include <ginit/gicomptype.h>
#include <ghyp/ghyp.h>

struct _G_PT_PROP_OPTS {
	boolean		use_unsplit_pt_prop;
	NptWSpeedOpts	npt_w_speed_opts;
};
typedef struct _G_PT_PROP_OPTS G_PT_PROP_OPTS;

enum _INTERPOLATION_OPTION {
    CONSERVATIVE_VARIABLES = 1,
    THERMODYNAMIC_VARIABLES = 2
};
typedef enum _INTERPOLATION_OPTION INTERPOLATION_OPTION;

struct _G_INIT_DATA {
	D_INIT_DATA	D_init_data;

	const char	*_ex_name;
	void	(*_prompt_for_composition_type)(INIT_DATA*);
	int	_material_composition_type;		/*PROMPTED*/

	void	(*_prompt_for_printing_and_plotting)(INIT_DATA*);

	Prompt_type	*_optional_printing_prompt_types;
	char	*_optional_printing_variables;		/*PROMPTED*/

	void	(*_prompt_for_hyperbolic_method)(INIT_DATA*);

	Muscl_Opts	_MusclOptions;

	AVISC	_global_artifical_viscosity_parameters;

	COMP_TYPE	*(*_get_comp_type_type)(LAYER_SYS*,int,int,
						INIT_PHYSICS*,INIT_DATA*);

	void	(*_prompt_for_ref_state)(const char*,Locstate,int,Gas_param*,
					 INIT_DATA*);

	void    (*_prompt_for_constant_flow_region)(COMPONENT,Locstate,
	                                            INTERFACE*);

	GRAVITY	*_gravity_data;

	void	(*_prompt_for_tracked_bifurcations)(INIT_DATA*);
	SCAT_WV_TOL *_scat_wv_tol_list;
	UNTRACK_NODE_FLAG *_untrack_node_list;

	void	(*_prompt_for_wave_capture_options)(INIT_DATA*);
	G_WAVE_CAPTURE	*_wave_capture_data;

	void	(*_prompt_for_point_propagation_options)(INIT_DATA*);
	G_PT_PROP_OPTS	_pt_propagate_opts;

	void	(*_prompt_for_interpolation_options)(INIT_DATA*);
	INTERPOLATION_OPTION	_interpolation_option;
	INTERPOLATION_WEIGHT	_rotational_symmetry_interpolation_flag;

	void	(*_prompt_for_maximum_number_of_components)(void);

	int     (*_prompt_for_bdry_wave_type)(INIT_DATA*,const char*,
					      const Prompt_type*);

	const Prompt_type *(*_bdry_wave_type_prompt_type)(void);
	int     _promptForBoundaryFlagsOffset;

	boolean	_reflect_small_loop_shocks;

	double _RL_eff;

	int   _type_of_state;
};
typedef struct _G_INIT_DATA G_INIT_DATA;
#define	g_init_data(init)	((G_INIT_DATA*)(init))

#define ex_name(init)	g_init_data(init)->_ex_name

#define	prompt_for_composition_type(init)				\
	(*g_init_data(init)->_prompt_for_composition_type)(init)
#define	material_composition_type(init)					\
	g_init_data(init)->_material_composition_type

#define	prompt_for_printing_and_plotting(init)				\
	(*g_init_data(init)->_prompt_for_printing_and_plotting)(init)

#define	optional_printing_prompt_types(init)				\
	g_init_data(init)->_optional_printing_prompt_types
#define	optional_printing_variables(init)				\
	g_init_data(init)->_optional_printing_variables

#define	prompt_for_hyperbolic_method(init)				\
	(*g_init_data(init)->_prompt_for_hyperbolic_method)(init)

#define	MusclOptions(init)		g_init_data(init)->_MusclOptions

#define	global_artifical_viscosity_parameters(init)			\
	g_init_data(init)->_global_artifical_viscosity_parameters

#define	get_comp_type_type(ls,ll,rl,ip,init)	\
	(*g_init_data(init)->_get_comp_type_type)(ls,ll,rl,ip,init)


#define prompt_for_constant_flow_region(comp,st,intfc,init)             \
        (*g_init_data(init)->_prompt_for_constant_flow_region)(comp,st,intfc)

#define	prompt_for_ref_state(msg,st,type,params,init)			\
	(*g_init_data(init)->_prompt_for_ref_state)(msg,st,type,params,init)

#define gravity_data(init)						\
	g_init_data(init)->_gravity_data

#define	prompt_for_tracked_bifurcations(init)				\
	(*g_init_data(init)->_prompt_for_tracked_bifurcations)(init)

#define	scattered_wave_tolerance_list(init)				\
	g_init_data(init)->_scat_wv_tol_list
	
#define	untrack_node_options_list(init)					\
	g_init_data(init)->_untrack_node_list

#define	prompt_for_wave_capture_options(init)				\
	(*g_init_data(init)->_prompt_for_wave_capture_options)(init)

#define	wave_capture_data(init)						\
	g_init_data(init)->_wave_capture_data

#define	f_wave_capture_data(init)					\
	(&g_init_data(init)->_wave_capture_data->F_wave_capture)

#define	prompt_for_point_propagation_options(init)			\
	(*g_init_data(init)->_prompt_for_point_propagation_options)(init)

#define	pt_propagate_opts(init)						\
	g_init_data(init)->_pt_propagate_opts

#define	prompt_for_interpolation_options(init)				\
	(*g_init_data(init)->_prompt_for_interpolation_options)(init)

#define	interpolation_option(init)					\
	g_init_data(init)->_interpolation_option

#define	rotational_symmetry_interpolation_flag(init)			\
	g_init_data(init)->_rotational_symmetry_interpolation_flag

#define	prompt_for_maximum_number_of_components(init)			\
	(*g_init_data(init)->_prompt_for_maximum_number_of_components)()

#define prompt_for_bdry_wave_type(init,mesg,ptypes)                     \
	(*g_init_data(init)->_prompt_for_bdry_wave_type)(init,mesg,ptypes)

#define bdry_wave_type_prompt_type(init)                                \
	(*g_init_data(init)->_bdry_wave_type_prompt_type)()

#define	reflect_small_loop_shocks(init)					\
	g_init_data(init)->_reflect_small_loop_shocks

#define	RL_eff(init)							\
	(g_init_data(init)->_RL_eff)

#define	type_of_state(init)						\
	(g_init_data(init)->_type_of_state)

#define	print_rarefaction_wave_1d(rw1d)	fprint_rarefaction_wave_1d(stdout,rw1d)

#include <ginit/ginitprotos.h>

#endif /* !defined(_GINIT_H) */
