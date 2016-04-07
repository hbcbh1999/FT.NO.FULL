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
*				ddecs.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains External Declarations for the Driver Code.
*/

#if !defined(_DDECS_H)
#define _DDECS_H

#include <hyp/hdecs.h>


#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif


struct _Grid {
			/* Space Grid */

	RECT_GRID	*rect_grid;	/* All rect grid parameters, XL etc. */

	PP_GRID *pp_grid;

			/* Time Grid */

	double	dt, dt_last;	/* Time Interval */
	double	_dt_lim;	/* Upper bound on time step */
	double	time;		/* Current Time */
	int	step;		/* Current Time Step */
	void	(*user_dt_control)(struct _Grid*,double*);

			/* Time Step Reduction Control */
	int	_max_num_mts;	/* Max. no. time step modifications per step */
	int	num_mts;	/* Current number of time step modifications */

	double	lscale;		/* Length/time scales for conversion */
	double	tscale;		/* to dimensioned units		     */


	struct _STOP	*_stop;
	struct _PAUSE	*_pause;

	INIT_DATA	*init;

			/* Miscellaneous Variables */
	void	(*_print_Grid_structure)(struct _Grid*);
	void	(*_print_remap_values)(void);

	char	*temporary_input_file;

	int initialization_complete;

			/* Information for repartition in output */
	boolean	repart_at_end_of_run;
	PP_GRID *new_pp_grid;
};
typedef struct _Grid Grid;

#define print_Grid_structure(grid)					\
	(*(grid)->_print_Grid_structure)(grid)
#define	dt_lim(s)			(s)->_dt_lim
#define	max_num_time_step_mods(s)	(s)->_max_num_mts

/* Macros for accessing Remap data */

#include <driver/damr.h>
#include <driver/dprt.h>

struct _STOP {
			/* Stopping Criteria */
	PrtMode _stop_time_mode;	/* EXACT_TIME or CONSTANT_TIME */
	double	_stop_time;		/* Max Total Time */
	int	_stop_step;		/* Max Number of Time Steps */
	int 	_stop_mode;          	/* Alternate Stopping Criterion */
	int	(*_stop_run)(Grid*,Front*,Wave*);
	void	(*_print_stopping_criteria)(Grid*);
};
typedef struct _STOP STOP;

#define	stop(grid)	((grid)->_stop)
#define	stop_time(grid)	(stop(grid)->_stop_time)
#define	stop_step(grid)	(stop(grid)->_stop_step)
#define	stop_mode(grid)	(stop(grid)->_stop_mode)
#define	stop_time_mode(grid)	(stop(grid)->_stop_time_mode)
#define stop_run(grid,front,wave)					\
	(*stop(grid)->_stop_run)(grid,front,wave)
#define print_stopping_criteria(grid)					\
	(*stop(grid)->_print_stopping_criteria)(grid)

struct _PAUSE {		/* Pause Criteria */
	PrtMode	_pause_mode;	/* Real vs. mesh time for pauses */
	double	_pause_time;	/* Time at which program pauses */
	int	_pause_step;	/* Mesh Time at which program pauses */
	void (*_await)(Grid*,Front*,Wave*,Printplot*,boolean);
};
typedef struct _PAUSE PAUSE;

#define	pause(grid)		((grid)->_pause)
#define	pause_mode(grid)	(pause(grid)->_pause_mode)
#define	pause_time(grid)	(pause(grid)->_pause_time)
#define	pause_step(grid)	(pause(grid)->_pause_step)
#define await(grid,front,wave,prt,got_intfc_from_file)			\
	(*pause(grid)->_await)(grid,front,wave,prt,got_intfc_from_file)

struct _INIT_PHYSICS {
	CHART		*root;
	Printplot	*prt;
	RECT_GRID	top_grid;
	double		initial_dt;
	void (*init_run)(int*,char***,INIT_DATA*,struct _INIT_PHYSICS*);
	void (*_init_physical_units)(INIT_DATA*);
	void (*_init_remap_and_rect_grid)(RECT_GRID*);
	void (*_init_physics)(INIT_DATA*,struct _INIT_PHYSICS*);
	void (*init_interface)(INIT_DATA*,struct _INIT_PHYSICS*);
	void (*init_parabolic)(Front*);
	void (*init_hyp_bdry)(Front*,void (*)(int,Locstate*,int,int));
	void (*init_bdry_state)(int,Locstate*,int,int);
	void (*init_cauchy_data_pointers)(struct _INIT_PHYSICS*,boolean);
	void (*cauchy_deposition)(Wave*,Front*);
	void (*initializer)(double*,COMPONENT,Locstate,INTERFACE*,INIT_DATA*);
	void (*intfc_initializer)(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				  Locstate,Locstate,INIT_DATA*);
	void (*pt_source_initializer)(Wave*,int,Locstate*);
	void (*restart_initializer)(int*,int,COMPONENT,
				    INPUT_SOLN**,Locstate,INIT_DATA*);
	void (*restart_intfc_initializer)(POINT*,HYPER_SURF_ELEMENT*,
					 HYPER_SURF*,Locstate,Locstate,
					 INIT_DATA*);
	void (*restart_pt_source_initializer)(Wave*,int,Locstate*);
	void (*init_point_sources)(Wave*);
	void (*set_Dirichlet_boundary_states)(Front*);

        void (*set_gravity_charts)(CHART**,int);
        boolean (*stratify_density_for_gravity)(void);
			
	void (*_init_stopping_criteria)(INIT_DATA*,struct _INIT_PHYSICS*);
	void (*_init_pause_criteria)(INIT_DATA*,struct _INIT_PHYSICS*);
#if defined(USE_AMR)
        void (*free_comp_types)(void);
#endif /* if defined(USE_AMR) */
};
typedef struct _INIT_PHYSICS INIT_PHYSICS;
/*
*	Initialization structure for driver library
*/

struct _D_INIT_DATA {
	H_INIT_DATA    H_init_data;
	double	       _initial_dt_lim;			/*PROMPTED*/
	double	       _initial_time;			/*PROMPTED*/
	STOP	       *_initial_stop;
	PAUSE	       *_initial_pause;
	int	       _max_num_mts;			/*PROMPTED*/
	char	       _output_filename[1024];		/*COMMAND LINE*/
	boolean	       _compress_output;		/*COMMAND LINE*/

	void (*_set_interface_hooks)(int,INIT_DATA*);

	PRINT_OPTIONS	_Default_print_options;

	void (*_prompt_for_printing_options)(INIT_DATA*);
	PRINT_OPTIONS	_prt_opts[NUM_OUTPUT_FORMATS];	/*PROMPTED*/

	PRINT_OPTIONS	*_geomview_opts;		/*PROMPTED*/

	void (*_prompt_for_physics_options)(INIT_DATA*,struct _INIT_PHYSICS*);

	Plot_choice	*(*_set_plot_choices)(INIT_DATA*);
	Plot_choice	*_plot_choices;			/*PROMPTED*/

	char	*_selected_tri_plots;			/*PROMPTED*/
	char	*_selected_gd_plots;			/*PROMPTED*/
	char	*_selected_PROSTAR_plots;			/*PROMPTED*/

#if defined(USE_HDF)
	HDF_PRINT_OPTIONS	_HDF_prt_options;	/*PROMPTED*/
	HDF_FRAME_OPTS		*_HDF_frame_options;	/*PROMPTED*/
	HDF_PRINT_OPTIONS	_SDS_prt_options;	/*PROMPTED*/
	HDF_FRAME_OPTS		*_SDS_frame_options;	/*PROMPTED*/
#endif /* defined(USE_HDF) */
       /* needed for VTK */
       VTK_PRINT_OPTIONS       _VTK_prt_options;       /*PROMPTED*/
       VTK_FRAME_OPTS          *_VTK_frame_options;    /*PROMPTED*/
       /* end needed for VTK */
#if defined(__GD__)
       /* needed for GD */
       GD_PRINT_OPTIONS       _GD_prt_options;       /*PROMPTED*/
       GD_FRAME_OPTS          *_GD_frame_options;    /*PROMPTED*/
       /* end needed for GD */
#endif /* defined(__GD__) */

        PROSTAR_PRINT_OPTIONS   _PROSTAR_prt_options;   /*PROMPTED*/
        PROSTAR_FRAME_OPTS	*_PROSTAR_frame_options;	/*PROMPTED*/
	void	(*_prompt_for_initial_intfc_options)(INIT_DATA*);
	void	(*_restart_clean_up)(INIT_DATA*);

	PRINT_OPTIONS	*_wall_time_dump_opts;
};
typedef struct _D_INIT_DATA D_INIT_DATA;
#define	d_init_data(init)	((D_INIT_DATA*)(init))
#define	initial_time(init)	d_init_data(init)->_initial_time
#define	initial_dt_lim(init)	d_init_data(init)->_initial_dt_lim
#define	initial_stop(init)	d_init_data(init)->_initial_stop
#define	initial_pause(init)	d_init_data(init)->_initial_pause
#define	initial_max_num_time_step_mods(init)	d_init_data(init)->_max_num_mts
#define	output_filename(init)	d_init_data(init)->_output_filename
#define	compress_output(init)	d_init_data(init)->_compress_output
#define	geomview_opts(init)	d_init_data(init)->_geomview_opts
#define initial_intfc_options(init) d_init_data(init)->_initial_intfc_options

#define	set_interface_hooks(dim,init)					\
    (*d_init_data(init)->_set_interface_hooks)(dim,init)

#define	Default_print_options(init)					\
    d_init_data(init)->_Default_print_options

#define	prompt_for_printing_options(init)				\
    (*d_init_data(init)->_prompt_for_printing_options)(init)
#define	prt_opts(init)	          d_init_data(init)->_prt_opts

#define	prompt_for_physics_options(init,ip)				\
    (*d_init_data(init)->_prompt_for_physics_options)(init,ip)

#define	prompt_for_initial_intfc_options(init)				\
    (*d_init_data(init)->_prompt_for_initial_intfc_options)(init)

#define set_plot_choices(init)						\
    (*d_init_data(init)->_set_plot_choices)(init)
#define	plot_choices(init)	d_init_data(init)->_plot_choices

#define	selected_tri_plots(init)	d_init_data(init)->_selected_tri_plots
#define	selected_gd_plots(init)		d_init_data(init)->_selected_gd_plots
#define	selected_PROSTAR_plots(init)	d_init_data(init)->_selected_PROSTAR_plots

/* needed for VTK */
#define VTK_prt_options(init)           d_init_data(init)->_VTK_prt_options
#define VTK_frame_options(init)         d_init_data(init)->_VTK_frame_options	
#if defined(USE_HDF)
/* end needed for VTK */
#define	HDF_prt_options(init)		d_init_data(init)->_HDF_prt_options
#define	HDF_frame_options(init)		d_init_data(init)->_HDF_frame_options
#define	SDS_prt_options(init)		d_init_data(init)->_SDS_prt_options
#define	SDS_frame_options(init)		d_init_data(init)->_SDS_frame_options
#endif /* defined(USE_HDF) */
#define	PROSTAR_prt_options(init)		d_init_data(init)->_PROSTAR_prt_options
#define	PROSTAR_frame_options(init)		d_init_data(init)->_PROSTAR_frame_options
#define	GD_prt_options(init)		d_init_data(init)->_GD_prt_options
#define	GD_frame_options(init)		d_init_data(init)->_GD_frame_options

#define	wall_time_dump_opts(init)					\
    d_init_data(init)->_wall_time_dump_opts

#define	restart_clean_up(init)				\
    (*d_init_data(init)->_restart_clean_up)(init)


#define	init_remap_and_rect_grid(grid,ip)				\
				(*(ip)->_init_remap_and_rect_grid)(grid)
#define	init_physical_units(init,ip)					\
	if ((ip)->_init_physical_units) (*(ip)->_init_physical_units)(init)

#define	init_stopping_criteria(init,ip)					\
				(*(ip)->_init_stopping_criteria)(init,ip)
#define	init_pause_criteria(init,ip)					\
				(*(ip)->_init_pause_criteria)(init,ip)

#define	init_physics(init,ip)						\
	if ((ip)->_init_physics != NULL) (*(ip)->_init_physics)(init,ip)

struct _D_Front {
	H_Front	hfront;
	struct {
		CHART	*chart;
	} d_extension;
};
typedef struct _D_Front D_Front;

#define d_front(front)			((D_Front *)(front))
#define DriverFrontExtension(fr)	d_front(fr)->d_extension
#define chart_of_front(fr)		DriverFrontExtension(fr).chart

struct _D_Wave {
	Wave	wave;

	struct {
		/* Variables for plotting */

		void (*_plot_states)(FILE*,Front*,Wave*,Printplot*);

		CHART   *chart;
	} d_extension;
};
typedef struct _D_Wave D_Wave;
#define d_wave(wave)			((D_Wave *)(wave))
#define	DriverWaveExtension(wave)	d_wave(wave)->d_extension
#define	plot_states_function(wave)	DriverWaveExtension(wave)._plot_states
#define chart_of_wave(wv)		DriverWaveExtension(wv).chart

		/* Macros to access plotting functions */
#define plot_states(file,front,wave,prt)				\
	if (plot_states_function(wave) != NULL)				\
		(*plot_states_function(wave))(file,front,wave,prt)

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#include <driver/dprotos.h>
#endif /* !defined(_DDECS_H) */
