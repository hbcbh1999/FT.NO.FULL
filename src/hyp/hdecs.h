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
*			 	hdecs.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	External definitions and declarations for hyperbolic code.
*/

#if !defined(_HDECS_H)
#define _HDECS_H

#include <tri/trigrid.h>

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

	/* Types of Point Sources */

enum {
	SOURCE	      = 1,
	SINK,
	NOT_SPECIFIED
};


	/* Macros */

	/*
	*	state pointer, x- and y-components of position, and
	*	component corresponding to the center of a grid square
	*/

#define Rect_coords(icoords,wave)					\
	(Coords(Regular_grid_node((icoords),wave_tri_soln(wave)->tri_grid)))

#define Rect_comp(icoords,wave)						\
	Regular_grid_comp((icoords),wave_tri_soln(wave)->tri_grid)

#define Rect_state(icoords,wave)					\
	Regular_grid_state((icoords),wave_tri_soln(wave)->tri_grid)

#define	Rect_blk_el0(icoords,wave)					\
	Regular_blk_el0((icoords),wave_tri_soln(wave)->tri_grid)

#define	Rect_blk_el1(icoords,wave)					\
	Regular_blk_el1((icoords),wave_tri_soln(wave)->tri_grid)

#define Rect_crossing(icoords,dir,wave)				\
	(*(wave)->_nearest_crossing)((icoords),dir,			\
				     wave_tri_soln(wave)->tri_grid)

	/*
	*	line integral of flux out of and area of the part
	*	of a grid square belonging to a given component
	*/

#define Rect_flux(icoords,comp,wave,flux_integral,front_state)		\
		(*wave_tri_soln(wave)->flux)((icoords),comp,		\
			wave_tri_soln(wave),flux_integral,front_state)

#define Rect_source(icoords,dt,wave,newwave,front,newfront,source,index,comp)\
		source_contributions(wave_tri_soln(wave),dt,wave,newwave,\
			front,newfront,					\
			Regular_grid_state((icoords),			\
				wave_tri_soln(wave)->tri_grid),		\
			source,Rect_coords(icoords,wave),icoords,index,comp)

#if defined(ONED)
#define comp_index_1d(comp,icoords,gmax)				\
	((comp)*(gmax)[0] + (icoords)[0])
#endif /* defined(ONED) */
#if defined(TWOD)
#define comp_index_2d(comp,icoords,gmax)				\
	(((comp)*(gmax)[1]+(icoords)[1])*(gmax)[0] + (icoords)[0])
#endif /* defined(TWOD) */
#if defined(THREED)
#define comp_index_3d(comp,icoords,gmax)				\
	((((comp)*(gmax)[2]+(icoords)[2])*(gmax)[1] +(icoords)[1])*(gmax)[0]\
	+(icoords)[0])
#endif /* defined(THREED) */
#if defined(ONED) && defined(TWOD) && defined(THREED)
#define comp_index(comp,icoords,gmax,dim)				\
		(((dim) == 2) ?	comp_index_2d(comp,icoords,gmax) :	\
		 ((dim) == 3) ?	comp_index_3d(comp,icoords,gmax) :	\
				comp_index_1d(comp,icoords,gmax))
#elif defined(TWOD) && defined(THREED)
#define comp_index(comp,icoords,gmax,dim)				\
		(((dim) == 2) ?	comp_index_2d(comp,icoords,gmax) :	\
		comp_index_3d(comp,icoords,gmax))
#elif defined(ONED) && defined(THREED)
#define comp_index(comp,icoords,gmax,dim)				\
		(((dim) == 1) ?	comp_index_1d(comp,icoords,gmax) :	\
				comp_index_3d(comp,icoords,gmax))
#elif defined(ONED) && defined(TWOD)
#define comp_index(comp,icoords,gmax,dim)				\
		(((dim) == 2) ?	comp_index_2d(comp,icoords,gmax) :	\
				comp_index_1d(comp,icoords,gmax))
#elif defined(ONED)
#define comp_index(comp,icoords,gmax,dim) comp_index_1d(comp,icoords,gmax)
#elif defined(TWOD)
#define comp_index(comp,icoords,gmax,dim) comp_index_2d(comp,icoords,gmax)
#elif defined(THREED)
#define comp_index(comp,icoords,gmax,dim) comp_index_3d(comp,icoords,gmax)
#endif /* defined(ONED) && defined(TWOD) && defined(THREED) */

#define Rect_area(ic,comp,wv,gmax,dim)					\
	wave_areas(wv)[comp_index((comp)-wave_min_comp(wv),ic,gmax,dim)]


	/*
	*	area integral of state over the part of a dual grid square
	*	belonging to a given component, divided by area of square
	*/

#define Rect_state_integral(icoords,comp,wave,area_integral,user)	\
		(*wave_tri_soln(wave)->integral)(icoords,comp,		\
			wave_tri_soln(wave),area_integral,(POINTER)user)


	/* Used to provide information on stencil positions */
	/*		for (irregular) boundaries	    */
	/*	     in different hyperbolic stencils	    */

typedef struct {
	boolean		reg_stencil, prev_reg_stencil;
	int		npts;
	int		**icoords, **icoords_store;
	int		*nc, *ncstore;
	COMPONENT	newcomp;
	COMPONENT	*comp, *compstore;
	HYPER_SURF	**hs, **hsstore;	/* TODO REMOVE */
	CRXING          ***crx, ***crxstore;
	POINT		**p, **pstore;
	Locstate	*st, *ststore;
	Locstate	*worksp, *worksp_store;
	byte		*worksp_st_store;
	Front		*fr, *newfr;
	struct _Wave	*wave, *newwave;
} Stencil;



struct _CLIP {
	struct _CLIP *prev, *next;
	COMPONENT comp;
#if defined(TWOD)
	int       ix, iy, num_dn;
	int       is_02, is_12, is_22,
		  is_01,        is_21,
	    	  is_00, is_10, is_20;
#endif /* defined(TWOD) */
    	double     incr_mass, mb_fact;
};
typedef struct _CLIP CLIP;

/* FOR AMR,  Overture init params */
/*Overture init params */
struct _OVERPARAM {
        double   errorThreshold;
        double   efficiency;
        int     numberOfSmooths;
        int     period;
        int     baseLevel;
        int     numberOfRefinementLevels;
        int     refinementRatio;
        int     extra_buf;
};
typedef struct _OVERPARAM Overparam;

struct _Wave {

	/* TIME DYNAMIC STRUCTURES */
	struct _Wave *old_wave;
	struct _WAVE_POINTERS {
	    struct _TRI_SOLN {

		    /* input variables */

	        INTERFACE	*intfc;	/* TODO: remove this */
	        TRI_GRID	*tri_grid;
	        size_t		sizest;

		TRI_GRID_HOOKS Tri_grid_hooks;

	        struct _INTERPOLATORS {
	            boolean (*linear_cell)(double*,LINEAR_ELEMENT*,
				           struct _TRI_SOLN*,Locstate);
	            boolean (*grad_linear_cell)(double**,LINEAR_ELEMENT*,
					        struct _TRI_SOLN*,Locstate*);
	            boolean (*least_square)(LEAST_SQR_CLUSTER*,struct _TRI_SOLN*,
		    			     double*,Locstate);
	            void    (*bilinear_cell)(double*,BILINEAR_ELEMENT*,
					     struct _TRI_SOLN*,Locstate);
	            void    (*grad_bilinear_cell)(double*,double*,
					          BILINEAR_ELEMENT*,
					          struct _TRI_SOLN*,Locstate*);
	            void (*grad_bond)(double,BOND*,Locstate,Locstate,Locstate*);
	        } interpolator;

	        struct _UNSPLIT {
	            void (*flux)(struct _TRI_SOLN*,int,Locstate,Locstate,
				 double,Locstate);
	            void (*flux_obl)(struct _TRI_SOLN*,Locstate,Locstate,
				     double*,Locstate);
	            void (*sources)(double,struct _Wave*,struct _Wave*,
				    Front*,Front*,Locstate,Locstate,
				    double*,int*,int,COMPONENT);
	        } unsplit;

	        struct _EL_INTEGRALS {
	            double (*linear_cell)(LINEAR_ELEMENT*,struct _TRI_SOLN*,
					 Locstate,POINTER);
	            double (*bilinear_cell)(BILINEAR_ELEMENT*,
					   struct _TRI_SOLN*,Locstate,POINTER);
	        } el_integral;


	            /* output variables */

	        boolean	(*_tri_solution)(double*,COMPONENT,struct _TRI_SOLN*,
				         Locstate,Locstate);
	        boolean	(*_grad_tri_solution)(double*,COMPONENT,
					      struct _TRI_SOLN*,Locstate*);
	        void	(*flux)(int*,COMPONENT,struct _TRI_SOLN*,Locstate,
				Locstate);
	        double	(*integral)(int*,COMPONENT,struct _TRI_SOLN*,
				    Locstate,POINTER);
#if defined(USE_OVERTURE)
            POINTER                  cg_over; /* pointer of CompositeGrid */
            POINTER                  cg_over_function; /* pointer of doubleCompositeGridFunction */
            int                      patch_number;
            int                      patch_level;
            int                      NumberOfLevels;
            int                      use_overture_state;
            int                      overture_init_step;
            COMPONENT                patch_component; /* it's for patches without interface */
#endif /* defined(USE_OVERTURE) */
	    } *tri_soln;		/* triangulated solution structure */

	    double *areas;		/* component areas of grid squares */
	    COMPONENT min_comp;		/* for use by area_of_component() */
	} wave_pointers;

	/* END TIME DYNAMIC STRUCTURES */


	/* TIME STATIC PARAMETERS */

		/* State Variable Specifications */

	size_t sizest;		/* size of Locstate in bytes */

		/* Optional storage saving flag */

	boolean min_storage;

		/* Real time */

	double time;

		/* Functions for copying and freeing Wave structures */

	struct _Wave*	(*_copy_wave)(struct _Wave*);
	void	(*_copy_into_wave)(struct _Wave*,struct _Wave*);
	void	(*_assign_wave_parameters)(struct _Wave*,struct _Wave*);
	void	(*_assign_wave_pointers)(struct _Wave*,struct _Wave*);
	void	(*_assign_copy_wave_pointers)(struct _Wave*,struct _Wave*);
	void	(*_clear_wave_pointers)(struct _Wave*);
	void	(*_free_wave)(struct _Wave*);
	void	(*_free_wave_pointers)(struct _Wave*);
	void	(*_free_copy_wave_pointers)(struct _Wave*);

	boolean	(*_detect_and_load_mix_state)(int,Stencil*,int);

		/* Variables for Debugging */

	int nfloats;		/* number of floats in a Locstate */
	/* prints a Locstate */
	void (*print_state)(Locstate);
	boolean (*_bad_state_data)(const char*,Front*,struct _Wave*);

	/* prints the wave Locstate's */
	void (*show_wave_states)(struct _Wave*);

	/* prints wave Locstates in tri_soln */
	void (*show_tri_soln)(Front*,struct _Wave*);

	/* allows direct plotting of hyp soln */
	void (*plot_hyp_soln)(Front*,struct _Wave*,int);

	/* parallel communication */
	boolean    (*_scatter_states)(struct _Wave*,Front*,int*,int);


		/* Variables to Define the Hyperbolic Method */
	
	const char *method;	/* name of method */
	int npts_sten;		/* stencil has npts_sten points */
	int sten_rad;		/* maximum stencil radius */
	int npts_vsten;		/* vector stencil has npts_sten points */

	RECT_GRID *rect_grid;		/* hyperbolic rectangular grid */

	PP_GRID* pp_grid;

                /* Mass Conservation Support Structures and Functions */

	/* Tri grid generation */

	TRI_GRID_HOOKS Tri_grid_hooks;

	/* split scheme stencil operators */
	void (*_npt_solver)(double,double,Locstate,const double*,
	                    int,int*,int*,Stencil*);

	/* vector scheme stencil operator */
	void (*_vec_solver)(int,int*,double*,struct _Wave*,struct _Wave*,
			    Front*,Front*,int*,int,int,double,double);

	/* alloc/free physics state uni_arrays */
	void (*_alloc_phys_vecs)(struct _Wave*,int);
	void (*_free_phys_vecs)(struct _Wave*);

#if defined(__cplusplus)

	/* unsplit scheme flux functions */
	struct _Wave::_WAVE_POINTERS::_TRI_SOLN::_UNSPLIT       unsplit;

	/* state interpolators */
	struct _Wave::_WAVE_POINTERS::_TRI_SOLN::_INTERPOLATORS interpolator;

	/* Integration over elements */
	struct _Wave::_WAVE_POINTERS::_TRI_SOLN::_EL_INTEGRALS  el_integral;

#else /* defined(__cplusplus) */

	/* unsplit scheme flux functions */
	struct _UNSPLIT	unsplit;

	/* state interpolators */
	struct _INTERPOLATORS	interpolator;

	/* Integration over elements */
	struct _EL_INTEGRALS	el_integral;

#endif /* defined(__cplusplus) */

	/* returns max speed of a state */
	double (*max_wave_speed)(Locstate);

	/* Computes hyp maximum time step */
	double (*max_hyp_time_step)(struct _Wave*,double*);

	/* Bundles state info */
	void (*bundle_states)(int*,int*,struct _Wave*,byte*);

	/* Unbundles state info */
	void (*unbundle_states)(int*,int*,struct _Wave*,byte*);

	/* For reaction term */
	void (*_react)(double,double*,int*,Locstate,COMPONENT);

	/* For mass conservation in interior */
	void (*clip_interior_soln)(int*,int,Locstate,POINTER,COMPONENT,
				   struct _Wave*,struct _Wave*,INTERFACE*,
				   double*,double*,double*,CLIP*);
	void (*_assign_excess_mass)(CLIP*,struct _Wave*,double*);

		/* Source Data */

			/* Point Sources */

	int num_point_sources;		/* number of sources or sinks */
	int *source_type;		/* type of source: SOURCE or SINK */
	int **isrcv;                    /* row, column influence flags */
	int **isrc_min, **isrc_max;     /* mesh columns influenced */
	double **srcpt;			/* coordinates of sources */
	double *strength;		/* strength, positive for source */
	double *prod_rate;		/* rate of production/injection */
	double **pt_src_diam;            /* rectangle influenced by src/sinks */
	void (*stat_prod_rate)(int,struct _Wave*,Front*,double*);
	Locstate *composition;		/* flow composition at source */

		/*Locstate allocation and clearing*/
	void	(*_alloc_state)(Locstate*,size_t);
	void	(*_clear_state)(Locstate,size_t);
	void	(*_obstacle_state)(Locstate,size_t);

		/* Grid crossing information */

	CRXING* (*_nearest_crossing)(int*,GRID_DIRECTION,TRI_GRID*);

		/* Internal Variables */

		/* wave speed accumlator */

	struct _h_MaxWaveSpeed *_MaxWaveSpeed;
	struct _h_MaxWaveSpeed *(*_alloc_MaxWaveSpeed)(struct _h_MaxWaveSpeed*,
		                                       struct _Wave*);
	struct _h_MaxViscosity *_MaxViscosity;
	struct _h_MaxViscosity *(*_alloc_MaxViscosity)(struct _h_MaxViscosity*,
		                                       struct _Wave*);

#if defined(USE_OVERTURE)
        POINTER                   cg_over;  /* pointer of CompositeGrid */
        POINTER                   cg_over_function; /* pointer of doubleCompositeGridFunction */
        int                       patch_number;
        int                       totalNumberOfPatches;
        int                       patch_level;
        int                       NumberOfLevels;
        int                       use_overture_state;
        int                       overture_init_step;
        COMPONENT                 patch_component; /* it's for patches without interface */
    /*  Neighbor                  *neighbor;  */
        Patch_bdry_flag           *pd_flag;
        void  (*overture_assign_wave_params)(Locstate,Locstate);
        void  (*overture_assign_wave_st_type)(Locstate,Locstate);
        void  (*overture_to_ft_st)(Locstate,POINTER,int,int*);
        void  (*ft_to_overture_st)(Locstate,POINTER,int,int*);


        void  (*trans_wv_st_to_overfunc)(struct _Wave*,Front*,POINTER,int,int);
        void  (*fill_root_extr_overture_cell_st_from_wv)(struct _Wave*,Front*,
                  POINTER,int,int);
        boolean  (*overture_init_interpolation_coarse_to_fine)(struct _Wave**,Front**);
        boolean  (*overture_interpolation_fine_to_coarse)(Wv_on_pc**,
                  struct _Wave**,Front**,struct _Wave***,Front***,int,int);
        boolean  (*overture_undistribute_interpolation_fine_to_coarse)(
                  struct _Wave**,Front**);
        boolean  (*scatter_patch_states)(Overparam*,Wv_on_pc**,struct _Wave**,
                  Front**,int*);
        boolean  (*scatter_patch_states_in_sweep_dir)(Overparam*,Wv_on_pc**,
                  struct _Wave**,Front**,int*,int);
        boolean  (*overture_fill_patch_amr_buffer)(struct _Wave**,Front**);
        boolean  (*overture_fill_patch_amr_buffer_pt)(int*,int,struct _Wave**,Front**);
        boolean  (*overture_injection_after_repatch)(struct _Wave**,Front**,
                                  struct _Wave**,Front**);
#endif /* defined(USE_OVERTURE) */

};
typedef struct _Wave Wave;
#if defined(__cplusplus)

typedef struct _Wave::_WAVE_POINTERS WAVE_POINTERS;
typedef struct _Wave::_WAVE_POINTERS::_TRI_SOLN TRI_SOLN;
typedef struct _Wave::_WAVE_POINTERS::_TRI_SOLN::_UNSPLIT UNSPLIT;
typedef struct _Wave::_WAVE_POINTERS::_TRI_SOLN::_INTERPOLATORS INTERPOLATORS;
typedef struct _Wave::_WAVE_POINTERS::_TRI_SOLN::_EL_INTEGRALS EL_INTEGRALS;

#else /* defined(__cplusplus) */

typedef struct _WAVE_POINTERS WAVE_POINTERS;
typedef struct _TRI_SOLN TRI_SOLN;
typedef struct _UNSPLIT UNSPLIT;
typedef struct _INTERPOLATORS INTERPOLATORS;
typedef struct _EL_INTEGRALS EL_INTEGRALS;

#endif /* defined(__cplusplus) */

#define npt_solver(dh,dt,ans,dir,swp_num,iperm,index,sten,wave)		\
	(*(wave)->_npt_solver)(dh,dt,ans,dir,swp_num,iperm,index,sten)
#define	scatter_states(wave,front,iperm,swp)				\
	(*(wave)->_scatter_states)(wave,front,iperm,swp)
#define	h_wave(wave)	((Wave*)wave)
#define bad_state_data(fnc,fr,wv)					\
    ((wv)->_bad_state_data != NULL) ? (*(wv)->_bad_state_data)(fnc,fr,wv) : NO
#define detect_and_load_mix_state(idir,sten,indx)                       \
        (*(sten->wave)->_detect_and_load_mix_state)(idir,sten,indx)

	/* Structure for tracking maximum hyperbolic wave speed */

struct	_h_MaxWaveSpeed {
	double		_maxsp[MAXD];	/* Max wave speeds in coord dirs    */
	double		**_coords;	/* Location of maximum wave speeds  */
	Locstate	*_mxspst;	/* Copy of state which set the      */
					/* maxium wave speed.		    */
	size_t		_sizest;
	struct _h_MaxWaveSpeedOperators {
		void	(*_set)(int,double,Locstate,double*,Wave*);
		void	(*_include)(struct _h_MaxWaveSpeed*,Wave*);
		void	(*_initialize)(Wave*);
		void	(*_print)(FILE*,Wave*);
		boolean	(*_read_print)(INIT_DATA*,const IO_TYPE*,Wave*);
		struct _h_MaxWaveSpeed *(*_copy)(Wave*);
		void    (*_destroy)(Wave*);
	} operators;
};
typedef struct _h_MaxWaveSpeed h_MaxWaveSpeed;

#if defined(__cplusplus)

typedef struct _h_MaxWaveSpeed::_h_MaxWaveSpeedOperators h_MaxWaveSpeedOperators;

#else /* defined(__cplusplus) */

typedef struct _h_MaxWaveSpeedOperators h_MaxWaveSpeedOperators;

#endif /* defined(__cplusplus) */

	/* Structure for tracking maximum viscosity */

struct	_h_MaxViscosity {
        double		_maxvisc;	/* Max viscosity */
	double		*_coords;	/* Location of maximum viscosity  */
	Locstate	*_mxviscst;	/* Copy of state which set the      */
					/* maxium viscosity.		    */
	size_t		_sizest;
	struct _h_MaxViscosityOperators {
		void	(*_set)(double,Locstate,double*,Wave*);
		void	(*_include)(struct _h_MaxViscosity*,Wave*);
		void	(*_initialize)(Wave*);
		void	(*_print)(FILE*,Wave*);
		boolean	(*_read_print)(INIT_DATA*,const IO_TYPE*,Wave*);
		struct _h_MaxViscosity *(*_copy)(Wave*);
		void    (*_destroy)(Wave*);
	} operators;
};
typedef struct _h_MaxViscosity h_MaxViscosity;

#if defined(__cplusplus)

typedef struct _h_MaxViscosity::_h_MaxViscosityOperators h_MaxViscosityOperators;

#else /* defined(__cplusplus) */

typedef struct _h_MaxViscosityOperators h_MaxViscosityOperators;

#endif /* defined(__cplusplus) */

#define alloc_MaxWaveSpeed(mxsp,wv)	(*(wv)->_alloc_MaxWaveSpeed)(mxsp,wv)
#define	MaxWaveSpeed(wv)	(wv)->_MaxWaveSpeed
#define	Maxsp(wv)		MaxWaveSpeed(wv)->_maxsp
#define	MaxWaveSpeedState(wv)	MaxWaveSpeed(wv)->_mxspst
#define	MaxWaveSpeedCoords(wv)	MaxWaveSpeed(wv)->_coords

#define MaxViscosity(wv)        (wv)->_MaxViscosity
#define Maxvisc(wv)             MaxViscosity(wv)->_maxvisc
#define MaxViscosityState(wv)   MaxViscosity(wv)->_mxviscst
#define MaxViscosityCoords(wv)  MaxViscosity(wv)->_coords

#define MaxWaveSpeedOperators(wv)	MaxWaveSpeed(wv)->operators
#define	SetMaxWaveSpeed(wv)		MaxWaveSpeedOperators(wv)._set
#define	IncludeMaxWaveSpeedInfo(wv)	MaxWaveSpeedOperators(wv)._include
#define	InitializeMaxWaveSpeed(wv)	MaxWaveSpeedOperators(wv)._initialize
#define	PrintMaxWaveSpeedInfo(wv)	MaxWaveSpeedOperators(wv)._print
#define	ReadPrintMaxWaveSpeedInfo(wv)	MaxWaveSpeedOperators(wv)._read_print
#define CopyMaxWaveSpeed(wv)		MaxWaveSpeedOperators(wv)._copy
#define DestroyMaxWaveSpeed(wv)		MaxWaveSpeedOperators(wv)._destroy

#define MaxViscosityOperators(wv)       MaxViscosity(wv)->operators
#define SetMaxViscosity(wv)             MaxViscosityOperators(wv)._set
#define IncludeMaxViscosityInfo(wv)     MaxViscosityOperators(wv)._include
#define InitializeMaxViscosity(wv)      MaxViscosityOperators(wv)._initialize
#define PrintMaxViscosityInfo(wv)       MaxViscosityOperators(wv)._print
#define ReadPrintMaxViscosityInfo(wv)   MaxViscosityOperators(wv)._read_print
#define CopyMaxViscosity(wv)            MaxViscosityOperators(wv)._copy
#define DestroyMaxViscosity(wv)         MaxViscosityOperators(wv)._destroy

#define	set_max_wave_speed(i,spd,st,crds,wv)				\
			(*SetMaxWaveSpeed(wv))(i,spd,st,crds,wv)
#define	include_max_wave_speed_info(mxsp,wv)				\
			(*IncludeMaxWaveSpeedInfo(wv))(mxsp,wv)
#define	initialize_max_wave_speed(wv)					\
			(*InitializeMaxWaveSpeed(wv))(wv)
#define	print_max_wave_speed_info(file,wv)				\
			(*PrintMaxWaveSpeedInfo(wv))(file,wv)
#define	read_print_max_wave_speed_info(init,io_type,wv)			\
			(*ReadPrintMaxWaveSpeedInfo(wv))(init,io_type,wv)
#define copy_max_wave_speed(wv)						\
			(*CopyMaxWaveSpeed(wv)._copy)(wv)
#define destroy_max_wave_speed(wv)					\
			(*DestroyMaxWaveSpeed(wv))(wv)
#define	set_max_viscosity(mu,st,crds,wv)				\
			(*SetMaxViscosity(wv))(mu,st,crds,wv)
#define	include_max_viscosity_info(mxvisc,wv)				\
			(*IncludeMaxViscosityInfo(wv))(mxvisc,wv)
#define	initialize_max_viscosity(wv)					\
			(*InitializeMaxViscosity(wv))(wv)
#define	print_max_viscosity_info(file,wv)				\
			(*PrintMaxViscosityInfo(wv))(file,wv)
#define	read_print_max_viscosity_info(file,wv)				\
			(*ReadPrintMaxViscosityInfo(wv))(file,wv)
#define copy_max_viscosity(wv)						\
			(*CopyMaxViscosity(wv)._copy)(wv)
#define destroy_max_viscosity(wv)					\
			(*DestroyMaxViscosity(wv))(wv)

	/* Macros to access wave copy/free functions */

#define	copy_wave(wave)			(*(wave)->_copy_wave)(wave)
#define	copy_into_wave(newwv,wv)	(*(wv)->_copy_into_wave)(newwv,wv)

#define	assign_wave_parameters(newwave,wave)				\
		(*(wave)->_assign_wave_parameters)(newwave,wave)

#define	assign_wave_pointers(newwave,wave)				\
		(*(wave)->_assign_wave_pointers)(newwave,wave)

#define	assign_copy_wave_pointers(newwave,wave)				\
		(*(wave)->_assign_copy_wave_pointers)(newwave,wave)

#define	clear_wave_pointers(wave)					\
		(*(wave)->_clear_wave_pointers)(wave)

#define	free_wave(wave)		      (*(wave)->_free_wave)(wave)

#define	free_wave_pointers(wave)      (*(wave)->_free_wave_pointers)(wave)

#define	free_copy_wave_pointers(wave) (*(wave)->_free_copy_wave_pointers)(wave)

struct _H_Front {
	Front	front;
	struct {
		Wave	*wave;
	} h_extension;
};
typedef struct _H_Front H_Front;
#define h_front(front)		((H_Front *)(front))
#define HypFrontExtension(fr)	h_front(fr)->h_extension
#define wave_of_front(fr)	HypFrontExtension(fr).wave

#define	wave_pointers(wave)	((wave)->wave_pointers)	/* NOT A POINTER */
#define	wave_tri_soln(wave)	(wave)->wave_pointers.tri_soln
#define	wave_areas(wave)	(wave)->wave_pointers.areas
#define	wave_min_comp(wave)	(wave)->wave_pointers.min_comp

#define	vec_solver(sn,ip,dr,wv,nwv,fr,nfr,ic,imn,imx,dt,dh)		\
	(*(wv)->_vec_solver)(sn,ip,dr,wv,nwv,fr,nfr,ic,imn,imx,dt,dh)


typedef struct {
	Prompt_type ptype;
	int         npts_sten;
	int         npts_vsten;
	int         sten_rad;
	int         (*hyp_driver)(double,double*,Wave*,Front*);
	void        (*_prompt_for_hyp_method_options)(INIT_DATA*);
} Hyp_method;

#define	prompt_for_hyp_method_options(init)				\
    if (hyperbolic_method(init)->_prompt_for_hyp_method_options != NULL) \
        (*hyperbolic_method(init)->_prompt_for_hyp_method_options)(init)

	/* Macros to access interpolator functions */

#define	Least_square_interpolate(eq,soln,crds,ans)			\
	(*(soln)->interpolator.least_square)(eq,soln,crds,ans)

#define	Bilinear_cell_interpolate(f,eq,soln,ans)			\
	(*(soln)->interpolator.bilinear_cell)(f,eq,soln,ans)

#define	Linear_cell_interpolate(f,t,soln,ans)				\
	(*(soln)->interpolator.linear_cell)(f,t,soln,ans)

#define	Grad_bilinear_cell_interpolate(f,d,eq,soln,ans)			\
	(*(soln)->interpolator.grad_bilinear_cell)(f,d,eq,soln,ans)

#define	Grad_linear_cell_interpolate(f,et,soln,ans)			\
	(*(soln)->interpolator.grad_linear_cell)(f,et,soln,ans)

#define	Grad_bond_interpolate(soln,tt,b,s1,s2,ans)			\
	(*(soln)->interpolator.grad_bond)(tt,b,s1,s2,ans)

	/* Macros to access integrator functions */

#define	Bilinear_cell_integrate(eq,soln,ans,user)			\
	(*(soln)->el_integral.bilinear_cell)(eq,soln,ans,user)

#define	Linear_cell_integrate(et,soln,ans,user)				\
	(*(soln)->el_integral.linear_cell)(et,soln,ans,user)

	/* Macros to access flux functions */

#define	flux_across_grid_segment(tri_soln,idir,s1,s2,dh,ans)		\
	if ((tri_soln)->unsplit.flux != NULL)				\
		(*(tri_soln)->unsplit.flux)(tri_soln,idir,s1,s2,dh,ans)

#define	flux_across_line_segment(tri_soln,s1,s2,dir,ans)		\
	if ((tri_soln)->unsplit.flux_obl != NULL)			\
		(*(tri_soln)->unsplit.flux_obl)(tri_soln,s1,s2,dir,ans)

#define source_contributions(tri_soln,dt,wv,newwv,fr,newfr,state,ans,	\
			     coords,icoords,index,comp)			\
	if ((tri_soln)->unsplit.sources != NULL)			\
		(*(tri_soln)->unsplit.sources)(dt,wv,newwv,fr,newfr,	\
				state,ans,coords,icoords,index,comp)

		/* Macros to access solution functions */

#define	tri_solution(coords,comp,tri_soln,state,dflt_state)		\
	(*(tri_soln)->_tri_solution)(coords,comp,tri_soln,state,dflt_state)

#define	grad_tri_solution(coords,comp,tri_soln,state)			\
	(*(tri_soln)->_grad_tri_solution)(coords,comp,tri_soln,state)

#define maximum_wave_speed(wave,state)					\
	( ((wave)->max_wave_speed != NULL) ?				\
		(*(wave)->max_wave_speed)(state) : 0.0 )

#define is_vector_method(wave)		((wave)->npts_vsten != 0)
#define stencil_radius(wave)		((wave)->npts_sten/2)
#define vsten_radius(wave)		((wave)->npts_vsten/2)

	/*
	*	Macros for allocation and free of state uni_arrays 
	*	for uni_arrayized schemes
	*/

#define alloc_phys_vecs(wv,vs)	(*(wv)->_alloc_phys_vecs)(wv,vs)

#define free_phys_vecs(wv)	(*(wv)->_free_phys_vecs)(wv)

/* Unsplit stencil structure */

struct _UnsplitStencilOptions {
	boolean use_hyp_solution;
};
typedef struct _UnsplitStencilOptions UnsplitStencilOptions;

struct _UnsplitStencil {
	COMPONENT new_comp, max_comp;

	/* Stencil geometry */
	double     *dirs[3], dir_store[9];
	double     dt, dh[3];
	int       npts, rad; /* Stencil is npts in each direction
			      * with radius rad*/
	int       nrad[3];
	int       *iperm;
	int       *idirs[3], idir_store[9];

	/* Indexed lists */
	COMPONENT *comp, *comp_list;
	Locstate  *state, *state_list; byte *state_store;
	double     **coords, **coords_list, *coords_store;
	int       **ic, **ic_list, *ic_store;
	boolean   *state_set, *state_set_list;

	Stencil   *sten;
	Front     *front, *newfront, *infront[3], *outfront[3];
	Wave      *wave,  *inwave[3],  *outwave[3];
	Wave      *wk_wv1, *wk_wv2, *tmpwave;
};
typedef struct _UnsplitStencil UnsplitStencil;

/*
*	Initialization structure for hyp library
*/

struct _H_INIT_DATA {
	TRI_INIT_DATA	Tri_init_data;		/*PROMPTED*/
	boolean	_use_mv_states;			/*PROMPTED*/
	UnsplitStencilOptions _USopts;		/*PROMPTED*/

	Hyp_method	*_hyperbolic_method;	/*PROMPTED*/
	Hyp_method	*_available_hyperbolic_methods;
	Hyp_method	*_default_hyperbolic_method;
	void	(*_setup_available_hyperbolic_methods_list)(INIT_DATA*);
	void    (*_init_wave_mv_state_list)(INIT_DATA*);
	void    (*_alloc_wave_mv_state_list)(INIT_DATA*,Wave*);
	void    (*_init_mv_list_states)(Wave*,Front*,INIT_DATA*,
					void (*)(double*,COMPONENT,Locstate,
						 INTERFACE*,INIT_DATA*));
	boolean _use_blocking_pp_comp;/*COMMAND LINE*/
};
typedef struct _H_INIT_DATA H_INIT_DATA;
#define	h_init_data(init)	((H_INIT_DATA*)(init))

#define	use_mv_states(init)	h_init_data(init)->_use_mv_states

#define	hyperbolic_method(init)	h_init_data(init)->_hyperbolic_method
#define	hyperbolic_method_name(init)					\
	hyperbolic_method(init)->ptype.type.ctype
#define	available_hyperbolic_methods(init)				\
	h_init_data(init)->_available_hyperbolic_methods
#define	default_hyperbolic_method(init)					\
	h_init_data(init)->_default_hyperbolic_method
#define	setup_available_hyperbolic_methods_list(init)			\
	(*h_init_data(init)->_setup_available_hyperbolic_methods_list)(init)
#define	USopts(init) h_init_data(init)->_USopts
#define init_wave_mv_state_list(init) 					\
	if (h_init_data(init)->_init_wave_mv_state_list != NULL)	\
	    (*h_init_data(init)->_init_wave_mv_state_list)(init)
#define alloc_wave_mv_state_list(init,wv)				\
	if (h_init_data(init)->_alloc_wave_mv_state_list != NULL)	\
	    (*h_init_data(init)->_alloc_wave_mv_state_list)(init,wv)
#define init_mv_list_states(wv,fr,init,initializer)			\
	if (h_init_data(init)->_init_mv_list_states != NULL)		\
	    (*h_init_data(init)->_init_mv_list_states)(wv,fr,init,initializer)
#define use_blocking_pp_comp(init)					\
	h_init_data(init)->_use_blocking_pp_comp

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#include <hyp/hprotos.h>
#endif /* !defined(_HDECS_H) */
