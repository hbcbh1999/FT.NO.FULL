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
*			dprotos.h
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if !defined(_DPROTOS_H)
#define _DPROTOS_H

#include <driver/ddecs.h>

	/* Driver EXPORTED Function Prototypes*/

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

#if defined(USE_OVERTURE)
        /* doverturepatch.c */
IMPORT  void    read_overture_option(INIT_DATA*, INIT_PHYSICS*);
IMPORT  void    overture_init_amr(CHART*);
IMPORT  int     d_overture_amr_main(INIT_DATA*, INIT_PHYSICS*);
IMPORT  double   find_amr_time_step(CHART**,int,Printplot*,boolean);
IMPORT  void    overture_dump_st_printout(CHART*,Printplot*);
        /* doverturepatch2.c */
IMPORT  int     pseudo_advance_amr_fronts2d(Front***,double,double*,
                     Overparam*,Wave**,Front**,int,Wv_on_pc**,int);
IMPORT  void    save_overture_show_files(Printplot*,Wave*,int);
#endif /* if defined(USE_OVERTURE) */

	/* dpatchmesh.c*/
IMPORT	CHART	*current_chart(void);
IMPORT	CHART	*read_print_CHART(const IO_TYPE*,CHART*,int*);
IMPORT	LEVEL	*read_print_LEVEL(const IO_TYPE*);
IMPORT	int	amr(CHART*,LEVEL*,Printplot*,int,double,double);
IMPORT	int	hyp_amr(double,double*,Wave*,Front*);
IMPORT	void	amr_hyp_grad_solution(double*,COMPONENT,HYPER_SURF*,SIDE,
				      Front*,POINTER,Locstate*);
IMPORT	void	amr_hyp_solution(double*,COMPONENT,HYPER_SURF*,SIDE,Front*,
				 POINTER,Locstate,Locstate);
IMPORT	void	init_amr(CHART*);
IMPORT	void	interpolate_from_finer_grids_to_level(LEVEL*);
IMPORT	void	print_CHART(FILE*,CHART*);
IMPORT	void	print_LEVEL(FILE*,LEVEL*);
IMPORT	void	reinit_amr_hyp_solution_functions(LEVEL*);
IMPORT	void	set_current_chart(CHART*);

	/* dinit.c*/
IMPORT	INTERFACE	*read_print_intfc_and_grids(INIT_DATA*,const IO_TYPE*,
                                                    RECT_GRID*,RECT_GRID*,
						    RECT_GRID*);
IMPORT	boolean	position_file_at_read_time(FILE*,int,double*,double*);
IMPORT	void	d_init_run(int*,char***,INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	d_init_stopping_criteria(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	d_prompt_for_initial_intfc_options(INIT_DATA*);
IMPORT	void	d_restart_clean_up(INIT_DATA*);
IMPORT	void	d_set_interface_hooks(int,INIT_DATA*);
IMPORT	void	read_next_dt(INIT_DATA*,const IO_TYPE*);
IMPORT	void	set_driver_hooks(INIT_DATA*,INIT_PHYSICS*);

	/* dinout.c*/
IMPORT	boolean	read_state_variables(const IO_TYPE*,int,INTERFACE*,
                                     INPUT_SOLN**,int);
IMPORT	double	is_state(INPUT_SOLN*,int*);
IMPORT	void	determine_read_version(FILE*);
IMPORT	void	d_print_states(FILE*,Printplot*);
IMPORT	void	d_print_states1d(FILE*,Printplot*);
IMPORT	void	free_input_soln(INPUT_SOLN*);
IMPORT	void	print_INPUT_SOLN_structure(INPUT_SOLN*);
IMPORT	void	print_OUTPUT_SOLN_structure(OUTPUT_SOLN*);
IMPORT	void	record_print_version(FILE*);
IMPORT	void	set_print_version(int);

	/* diprt.c */
IMPORT	FILE	*open_data_file(Front*,const char*,boolean,boolean,
				const char*,char**,const char*,char**);
IMPORT	char*	read_plotting_choices(char*,INIT_DATA*);
IMPORT	void	add_user_output_function(void(*)(Grid*,Wave*,Front*,
						 Printplot*,OUTPUT_DATA*,
						 boolean),
					 OUTPUT_DATA*,Printplot*);
IMPORT	void	copy_output_data(OUTPUT_DATA*,OUTPUT_DATA*,boolean,char*);
IMPORT	void	create_output_data_file(OUTPUT_DATA*);
IMPORT	void	d_prompt_for_printing_options(INIT_DATA*);
IMPORT	void	d_init_printplot(INIT_DATA*,Grid*,Front*,Printplot*);
IMPORT	void	init_peep_time(INIT_DATA*,Grid*,Printplot*);
IMPORT	void	init_output_data(INIT_DATA*,OUTPUT_DATA*,Grid*,Printplot*,
				 const char*,boolean,boolean,boolean);
IMPORT	void	print_plotting_choices(const char*,Plot_choice*);
IMPORT	void	prompt_for_output_filenames(PRINT_OPTIONS*,boolean,const char*);
IMPORT	void	prompt_for_tri_plots(INIT_DATA*);
IMPORT	void	prompt_for_gd_plots(INIT_DATA*);
IMPORT	void	prompt_for_prostar_plots(INIT_DATA*);
IMPORT	void	set_defaults_for_print_options(PRINT_OPTIONS*,INIT_DATA*);
IMPORT	void	set_next_print_time(Grid*,PRINT_CONTROL*);
IMPORT void	set_output_data(PRINT_OPTIONS*,OUTPUT_DATA*,Grid*,Printplot*,
			    boolean,boolean);
IMPORT	void	init_cross_sections(INIT_DATA*,Front*,Grid*,Printplot*);
IMPORT	void	init_tri_plots(INIT_DATA*,Printplot*);
IMPORT	void	init_gd_movie_plots(INIT_DATA*,Printplot*);
#if defined(USE_HDF)
IMPORT	void	init_HDF_plots(HDF_DATA_TYPE,INIT_DATA*,Printplot*);
IMPORT	void	prompt_for_hdf_plots(HDF_DATA_TYPE,INIT_DATA*);
IMPORT	void	set_hdf_append(Printplot*);
#endif /* defined(USE_HDF) */
/* needed for VTK */
IMPORT  void    prompt_for_vtk_plots(INIT_DATA*);
IMPORT  void    init_vtk_plots(INIT_DATA*,Printplot*);
/* end needed for VTK */
IMPORT	void	init_PROSTAR_plots(INIT_DATA*,Printplot*);
IMPORT  void    init_gd_movie_plots(INIT_DATA*,Printplot*);

	/* dmain.c*/
IMPORT	int	d_stop_run(Grid*,Front*,Wave*);
IMPORT	int	dmain(int,char**,INIT_DATA*,INIT_PHYSICS*);
IMPORT	int	time_step(CHART*,double,double*,int,double);
IMPORT	void	d_clean_up_printout(int);
IMPORT	void	d_clean_up(void);
IMPORT	void	perform_initialization(int,char**,INIT_DATA*,INIT_PHYSICS*);

	/* dprint.c*/
IMPORT	OUTPUT_DATA	**d_alloc_output_datas(int);
IMPORT	boolean	is_print_step(int,PRINT_CONTROL*);
IMPORT	boolean	is_print_time(double,PRINT_CONTROL*);
IMPORT	double	find_time_step(Grid*,Wave*,Front*,Printplot*,boolean);
IMPORT	double	nonphysics_timestep_reductions(Grid*,Wave*,Front*,
					       Printplot*,double);
IMPORT	void	adjoin_node_number(char*);
IMPORT	void	d_await(Grid*,Front*,Wave*,Printplot*,boolean);
IMPORT	void	d_print_D_Front_structure(Front*);
IMPORT	void	d_print_initial_data(FILE*,CHART*,Printplot*);
IMPORT	void	d_printout(CHART*,Printplot*,boolean,int);
IMPORT	void	d_print_stopping_criteria(Grid*);
IMPORT	void	give_peep(Grid*,Wave*,Front*,Printplot*,OUTPUT_DATA*,boolean);
IMPORT  void    interactive_printing_control(Printplot*,Grid*,boolean*,boolean*);
IMPORT	void	plot_states1d(FILE*,Front*,Wave*,Printplot*);
IMPORT	void	plot_states2d(FILE*,Front*,Wave*,Printplot*);
IMPORT	void	plot_states3d(FILE*,Front*,Wave*,Printplot*);
IMPORT	void	set_output_file_name(int,char*,const char*,const int,const int);

	/* dstat.c*/
IMPORT	void	d_dual_grid_ident_statistics3d(Grid_Stats_data*,
					       Grid*,Wave*,Front*,int);
IMPORT	void	d_init_grid_statistics(INIT_DATA*,Front*,Grid*,Wave*,Printplot*,
				       void(*)(Grid*,Wave*,Front*,Printplot*,
					       OUTPUT_DATA*,boolean));
IMPORT	void	d_reg_grid_cyl_statistics2d(Grid_Stats_data*,
					    Grid*,Wave*,Front*,int);
IMPORT	void	d_reg_grid_ident_statistics2d(Grid_Stats_data*,
					      Grid*,Wave*,Front*,int);

	/* dsub.c*/
IMPORT	const char *output_format_name(int);
IMPORT	boolean	regrid_front_and_wave(Front*,Wave*,RECT_GRID*,RECT_GRID*,boolean);
IMPORT	void	d_assign_wave_parameters(Wave*,Wave*);
IMPORT	void	d_copy_into_front(Front*,Front*);
IMPORT	void	d_copy_into_wave(Wave*,Wave*);
IMPORT	void	d_print_Grid_structure(Grid*);
IMPORT	void	d_set_default_front_parameters(INIT_DATA*,Front*);
IMPORT	void	d_set_default_wave_parameters(INIT_DATA*,Wave*);
IMPORT	void	delete_untracked_hyper_surfaces(Front*,Wave*);
IMPORT	void	print_Printplot_structure(Printplot*);
IMPORT	void	print_cross_sectional_graphs(Grid*,Wave*,Front*,Printplot*,
					     OUTPUT_DATA*,boolean);

	/* ioshock.c*/
IMPORT	void	cdumps(float,float,int,int,int,float*,float*,int,float,
		       float**,float*,int,int);
IMPORT	void	setdmp(char*,int);

#endif /* !defined(_DPROTOS_H) */
/*gstglobs.c*/
IMPORT void print_gravity(INTERFACE *);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif
