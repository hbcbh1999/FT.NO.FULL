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
*				ginitprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Function prototypes for the ginit library.
*/

#if !defined(_GINITPROTOS_H)
#define _GINITPROTOS_H

/*	gibifur.c */
IMPORT	void	g_prompt_for_tracked_bifurcations(INIT_DATA*);

/*	gicc.c */
IMPORT	void	init_CC_interaction(INIT_DATA*,INIT_PHYSICS*);

/*	gielrp.c */
IMPORT	ELLIPSOID	*allocate_ellipsoid(ELLIPSOID*,int);
IMPORT	ELLIPSOID	*prompt_for_ellipsoid(Front*,const char*,double*,
					      double*,const LAYER_FLAG*,
					      INIT_PHYSICS*);
IMPORT	double	max_radii(ELLIPSOID*);
IMPORT	void	init_el_riem_prob(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	print_ellipsoid(ELLIPSOID*,INTERFACE*);
IMPORT	void	set_default_ellipsoid_structure(ELLIPSOID*,int);
IMPORT	void	set_elliptical_comp_type(COMP_TYPE*,INIT_PHYSICS*);

/*	giglobs.c */
IMPORT	COMP_TYPE	*comp_type(COMPONENT);
IMPORT	const COMPONENT *comps_with_type(COMP_TYPE_TYPE,int*);
IMPORT	int	g_init_composition_type(INIT_PHYSICS*,INIT_DATA*,size_t*,int*);
IMPORT	int	max_num_comps(void);
IMPORT	void	check_float_input(const char*,double,double,double,int);
IMPORT	void	check_int_input(const char*,int,int,int,int);
IMPORT	void	g_compute_sizest(int,size_t*,int*,int);
IMPORT	void	g_prompt_for_composition_type(INIT_DATA*);
IMPORT	void	g_prompt_for_maximum_number_of_components(void);
IMPORT	void	get_state_ambient(double*,Locstate,COMP_TYPE*,
				  HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
IMPORT	void	get_state_2d_riemann(double*,Locstate,COMP_TYPE*,
				  HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
IMPORT	void	free_comp_types(void);
IMPORT	void	prompt_for_gravity(INIT_DATA*,INIT_PHYSICS*);
IMPORT	boolean	overlays_exist(void);
IMPORT	void	get_state_1d_overlay(double*,Locstate,COMP_TYPE*,
				     HYPER_SURF*,INTERFACE*,INIT_DATA*,int);

/*	gictype.c */
IMPORT	COMP_TYPE	*g_get_comp_type_type(LAYER_SYS*,int,int,
				     INIT_PHYSICS*,INIT_DATA*);
IMPORT	void	prompt_for_comp_type_type(COMP_TYPE*,const char*,INIT_PHYSICS*);
IMPORT	void	set_rarefaction_wave_1d_comp_type(COMP_TYPE*,Front*);

/*	girt.c */
IMPORT	STRATIFICATION_TYPE prompt_for_stratification(const char*);
IMPORT	void	get_state_in_stratified_region(STRATIFICATION_TYPE,Locstate,
					       double,double,Locstate);
IMPORT	void	get_state_rt_kh_perturbed(double*,Locstate,COMP_TYPE*,
					  HYPER_SURF*,INTERFACE*,INIT_DATA*,
					  int);
IMPORT	void	free_rt_kh_comp_type(COMP_TYPE*);
IMPORT	void	set_rt_kh_comp_type(COMP_TYPE*,Front*);

/*	gilayer.c */
IMPORT	LAYER_SURF	*alloc_layer_surf(void);
IMPORT  void    g_prompt_for_constant_flow_region(COMPONENT,Locstate,
                                                  INTERFACE*);
IMPORT	boolean	rarefaction_edge_at_coords(double*,HYPER_SURF*,
					   RAREFACTION_EDGE_TYPE);
IMPORT	double	get_surf_height(double*,LAYER_SURF*);
IMPORT	void	init_multi_layer(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	set_stretching_comp_type(COMP_TYPE*,Front*);
IMPORT	void	set_2d_riemann_comp_type(COMP_TYPE*,Gas_param*,INIT_DATA*);
IMPORT	void	prompt_for_1d_overlay(INIT_DATA*,COMP_TYPE*,INIT_PHYSICS*,
				      const char*,OVERLAY_TYPE);
IMPORT	void	set_1d_overlay_comp_type(COMP_TYPE*);
IMPORT	LAYER   *prompt_for_each_layer(int,int,const LAYER_FLAG*,Front*,
				      LAYER**,INIT_PHYSICS*);
IMPORT	void    make_layer_surf(Front*,LAYER_SURF*,INIT_DATA*);
IMPORT	void    init_comp_type(INIT_DATA*,INIT_PHYSICS*,LAYER_SYS*);
IMPORT	void	prompt_for_rigid_body_params(RIGID_BODY_PARAMS*,
				 	const char*,int);

/*	girgb.c */
IMPORT	void	init_fluid_rigid_body(INIT_DATA*,INIT_PHYSICS*);

/*	girpregion.c */
IMPORT	void	set_up_riemann_problem_region(int,int,int,int,double*,
                                              double*,SIDE,Locstate,Locstate,
					      LAYER_SYS*,INIT_PHYSICS*,
					      INIT_DATA*);

/*	gimkcur.c */
IMPORT	HYPER_SURF	*g_make_ellipse(ELLIPSOID*,COMPONENT,COMPONENT,Front*);
IMPORT	void	coords_on_pert_ellipsoid(double*,double*,ELLIPSOID*);
IMPORT	void	g_make_ellip_region_boundaries2d(ELLIPSOID*,Front*);
IMPORT	void	g_make_fourier_curve(int,int,double,double,FOURIER_POLY*,
				     COMPONENT,COMPONENT,double);
IMPORT	void	make_neumann_curve(INTERFACE*,NODE*,NODE*,double,double,double,
				   double,int,int,COMPONENT,COMPONENT,CURVE**);
IMPORT	void	make_rotated_curve(double*,double,double*,double*,
				   double (*f) (double,POINTER),
				   double (*fprime) (double,POINTER),
				   POINTER,NODE*,COMPONENT,COMPONENT,
				   CURVE**,ORIENTATION,
				   int,int,int,double,Front*);
IMPORT	void	make_vertical_axis_curve(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	make_ramp(int,double,double,double,double,COMPONENT,COMPONENT,
			  RECT_GRID*,CURVE**);
IMPORT	void	set_tangent_to_origin(TANGENT_FUNCTION*);
IMPORT	double	sin_sqr_pert(double,POINTER);
IMPORT	double	sin_sqr_pert_prime(double,POINTER);
IMPORT	void	set_normal_to_origin(NORMAL_FUNCTION*);

/*	gimksurf.c */
IMPORT	HYPER_SURF *g_make_ellipsoid(ELLIPSOID*,COMPONENT,COMPONENT,Front*);
IMPORT	void	   g_make_ellip_region_boundaries3d(ELLIPSOID*,Front*);
IMPORT	void	   make_random_surface(int,Front*,SINE_PERT*,COMPONENT,
				       COMPONENT,double,int);
/*	gimkbub3d.c  */
IMPORT  void    g_prompt_for_dynamic_bubble_insertion3d(CHART*);

/*	ginitintfc.c */
IMPORT	const char *comp_type_name(COMP_TYPE_TYPE);
IMPORT  int     g_prompt_for_bdry_wave_type(INIT_DATA*,
					    const char*,const Prompt_type*);
IMPORT	int	g_prompt_for_boundary_state(int,const char*,double*,
					    COMPONENT,int,HYPER_SURF*,
					    INIT_DATA*,INIT_PHYSICS*);
IMPORT	int	g_prompt_for_wave_type(const char*,INTERFACE*,INIT_PHYSICS*);
IMPORT	void	g_init_cauchy_data_pointers(INIT_PHYSICS*,boolean);
IMPORT	void	g_init_interface(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	g_init_problem_type(INIT_PHYSICS*);
IMPORT  void    g_init_parabolic(Front*);
IMPORT	void	gi_set_interface_hooks(int,INIT_DATA*);
IMPORT	void	prompt_for_ambient_state(COMP_TYPE*,Gas_param*,const char*,
					 Front*,INIT_DATA*);
IMPORT	void	prompt_for_problem_type(INIT_PHYSICS*,Prob_type*);
IMPORT	void	g_prompt_for_ref_state(const char*,Locstate,int,
				       Gas_param*,INIT_DATA*);
IMPORT  void    read_time_dependent_data(TD_BSTATE*);
IMPORT	void	set_ambient_boundary_state(Locstate,double*,COMPONENT,
					   INTERFACE*,INIT_DATA*);
IMPORT	void	set_ambient_comp_type(COMP_TYPE*,Front*);
IMPORT	void	set_obstacle_comp_type(COMP_TYPE*,Front*);

/*      giphysprompt.c */
IMPORT  void    g_prompt_for_physics_options(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	g_set_prompting_hooks(INIT_DATA*);
IMPORT	void	init_artificial_parameters(INIT_DATA*);

/*	ginitphys.c */
IMPORT	void	g_init_physics(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	g_copy_into_front(Front*,Front*);
IMPORT	void	g_copy_into_wave(Wave*,Wave*);
IMPORT	void	g_set_basic_phys_parameters(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	read_print_stratified_state_function(const IO_TYPE*,INTERFACE*);
IMPORT	void	set_gas_hooks(INIT_DATA*,INIT_PHYSICS*);

/*	giparams.c */
IMPORT	Gas_param	*g_init_eos_params(INIT_DATA*,INIT_PHYSICS*,
					   const char*,boolean);
IMPORT	Gas_param	*g_prompt_for_eos_params(INIT_DATA*,INIT_PHYSICS*,
						 boolean,const char*);
IMPORT  Gas_param       **g_prompt_for_eos_params_list(INIT_DATA*,INIT_PHYSICS*,
                                                       boolean,int*);
IMPORT	void	prompt_for_artificial_viscosity_and_heat_conduction(INIT_DATA*,
								    const char*,
								    const char*,
								    boolean,
								    AVISC*);
IMPORT	void	read_print_flame_velocity_params(INIT_DATA*,const IO_TYPE*,Gas_param*);

/*	gipert.c */
IMPORT	double	get_rho_prime(Locstate,COMP_TYPE*,double);
IMPORT	void	init_kelvin_helmholtz(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	init_random_surface(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	lin_pert_state(Locstate,double*,Front*,Wave*,int);
IMPORT	void	rt_ml_linear_pert(LAYER_SYS*,INIT_DATA*);

/*	gijet.c */
IMPORT	void	init_fuel_injection_jet(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	init_injection_inlet_jet(INIT_DATA*,INIT_PHYSICS*);

/*      gi3comp.c */
IMPORT	void	init_3comp_jet3d(INIT_DATA*,INIT_PHYSICS*);

/*      giboone.c */
IMPORT  void    init_neutrino_booster_detector(INIT_DATA*,INIT_PHYSICS*);
IMPORT  void    set_state_cracked_pmt_bottom(double*,COMPONENT,Locstate);

/*	gipppert.c */
IMPORT	int	pp_get_next_mode(int*,int,LIN_PERT*);
IMPORT	int	pp_list_modes(NORMAL_MODE**,LIN_PERT*,int,int,int);
IMPORT	void	brdcst_info(LIN_PERT*,NORMAL_MODE**,int,int*);

/*	giprt.c */
IMPORT	INPUT_SOLN	**g_set_input_solution(Printplot*);
IMPORT	PRINTING_LIST	*g_set_printing_list(INIT_PHYSICS*,INIT_DATA*,
					     Printplot*);
IMPORT	PROBE	*init_probe(double*,double*,int*,const double*,
                            INIT_DATA*,INIT_PHYSICS*);
IMPORT	Plot_choice	*g_set_plot_choices(INIT_DATA*);
IMPORT	void	fprint_rarefaction_wave_1d(FILE*,_RAREFACTION_WAVE_1D*);
IMPORT	void	g_init_printing_and_plotting(INIT_PHYSICS*,INIT_DATA*,
					     Printplot*,int);
IMPORT	void	g_init_statistics(Front*,Grid*,Wave*,Printplot*,INIT_DATA*);
IMPORT	void	g_prompt_for_optional_printing_variables(INIT_DATA*);
IMPORT	void	g_prompt_for_printing_and_plotting(INIT_DATA*);

/*	gireadstate.c */
IMPORT	Locstate g_read_print_state_data(INIT_DATA*,const IO_TYPE*,
                                         Locstate,INTERFACE*);
IMPORT	int	read_Gas_param_list(INIT_DATA*,const IO_TYPE*,uint64_t**,
				    Gas_param***,int,size_t,boolean);
IMPORT	void	g_free_restart_params_list(void);
IMPORT	void	g_read_print_ContactWallNodeParams(const IO_TYPE*,INTERFACE*);
IMPORT	void	g_read_print_Dirichlet_bdry_states(INIT_DATA*,const IO_TYPE*,
                                                   INTERFACE*);
IMPORT	void	g_read_print_Gas_param_list(INIT_DATA*,const IO_TYPE*,
                                            INTERFACE*);
IMPORT	void	g_read_print_RP_DATA_at_nodes(INIT_DATA*,const IO_TYPE*,
                                              INTERFACE*);
IMPORT  void    g_read_print_boundary_state_data(INIT_DATA*,const IO_TYPE*,
                                                 INTERFACE*,int);
IMPORT	void	read_print_avisc_params(INIT_DATA*,const IO_TYPE*,int,AVISC*);
IMPORT	void	read_print_gas_data(const IO_TYPE*,Locstate*,int,int,size_t,int,
				    uint64_t*,Gas_param**,int);
IMPORT	void	read_print_thermodynamic_restrictions(Gas_param*,
                                                      const IO_TYPE*);
IMPORT	void	reset_artificial_viscosity_and_heat_conduction(INIT_DATA*);
IMPORT	void	set_restart_params(Gas_param**,uint64_t,int,
				   uint64_t*,Gas_param**);

/*	girefl.c */
IMPORT	void	init_ramp_reflection(INIT_DATA*,INIT_PHYSICS*);

/*	girestrt.c */
IMPORT	void	g_restart_initializer(int*,int,COMPONENT,
				      INPUT_SOLN**,Locstate,INIT_DATA*);
IMPORT	void	g_restart_intfc_initializer(POINT*,HYPER_SURF_ELEMENT*,
					    HYPER_SURF*,Locstate,Locstate,
					    INIT_DATA*);
IMPORT	void	g_restart_set_intfc_states(double*,double*,int,POINT*,
					   HYPER_SURF_ELEMENT*,HYPER_SURF*);

/*	girstate.c.c */
IMPORT	RANDOM_STATE	*allocate_random_state_structure(INTERFACE*);
IMPORT	RANDOM_STATE	*copy_random_state_structure(RANDOM_STATE*,INTERFACE*);
IMPORT	int	prompt_for_random_flow_inlet(double*,COMPONENT,HYPER_SURF*,
					     INTERFACE*,int,INIT_DATA*);
IMPORT	void	free_random_state_structure(RANDOM_STATE*);
IMPORT	void	init_random_state_region(Locstate,RANDOM_STATE*,Front*);
IMPORT	void	random_region_state(Locstate,double*,RANDOM_STATE*);
IMPORT	void	read_print_random_velocity_inlet_data(INIT_DATA*,const IO_TYPE*,
						      BOUNDARY_STATE*,
						      INTERFACE*);
IMPORT	void	set_random_region_comp_type(COMP_TYPE*,Front*);

/*	gisc.c */
IMPORT	void	set_untracked_shock_wave_comp_type(COMP_TYPE*,UT_SHOCK*,Front*);
IMPORT	void	init_meshkov(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	init_shock_diffraction(INIT_DATA*,INIT_PHYSICS*);
IMPORT	void	init_shock_transmission(INIT_DATA*,INIT_PHYSICS*);

/*	gistate.c */
IMPORT	_RAREFACTION_WAVE_1D *allocate_RAREFACTION_WAVE_1D(Front*);
IMPORT	int	init_state_type(int);
IMPORT	void	init_shock_states(double,double,double,double*,Gas_param*,
				  Locstate,Locstate);
IMPORT	void	prompt_for_behind_contact_state(Locstate,Locstate,
						Gas_param*,boolean,double*,int,
						INIT_DATA*);
IMPORT	void	prompt_for_behind_shock_state(Locstate,Locstate,boolean,double*,
					      int,boolean,INIT_DATA*);
#if defined(COMBUSTION_CODE)
IMPORT	void	prompt_for_burning(Gas_param**,const char*);
#endif /* defined(COMBUSTION_CODE) */

/*    girm_linear.c */
IMPORT	void	rm_linear(INIT_DATA*,INIT_PHYSICS*,LAYER_SYS*);

/*    gihypinit.c */
IMPORT	void	g_set_hyp_solvers(INIT_DATA*,Wave*,Front*);
IMPORT	void	g_setup_available_hyperbolic_methods_list(INIT_DATA*);

/*    gitabregion.c */
IMPORT	void	set_tabulated_region_comp_type(COMP_TYPE*,Front*);
IMPORT	void	set_up_read_state_from_file_region(COMP_TYPE*,const char*,
                                                   Front*);

/*    spolars.c */
IMPORT	int	sp_main(int,char**,INIT_DATA*,INIT_PHYSICS*);

/*    testsolver.c */
IMPORT	int	gsolv_main(int,char**,INIT_DATA*,INIT_PHYSICS*);

/*    gcauchy.c */
IMPORT  void    g_cauchy_deposition(Wave*,Front*);

#endif /* !defined(_GINITPROTOS_H) */
