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
*				gprtprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GPRTPROTOS_H)
#define _GPRTPROTOS_H

	/* gprt EXPORTED Function Declarations */

/*	gprint.c */
IMPORT	Gas_param  *gas_params_for_comp(COMPONENT,INTERFACE*);
IMPORT	const char *avisc_print_style(void);
IMPORT	const char *print_wave_family(FILE*,const char*,WAVE_FAMILY,
				      const char*);
IMPORT	void	evaluate_probe(Grid*,Wave*,Front*,Printplot*,OUTPUT_DATA*,boolean);
IMPORT	void	fprint_avisc_structure(FILE*,AVISC*,boolean);
IMPORT	void	fprint_thermodynamic_restrictions(FILE*,Gas_param*);
IMPORT	void	g_FPrintWSSten(FILE*,WSSten*);
IMPORT	void	g_PrintUSWSSten2d(USWSSten2d*);
IMPORT	void	g_fprint_Dirichlet_bdry_states(FILE*,INTERFACE*);
IMPORT	void	g_fprint_Gas_param(FILE*,Gas_param*);
IMPORT	void	g_fprint_Gas_param_list(FILE*,INTERFACE*);
IMPORT  void    g_fprint_front(Front*,FILE*);
IMPORT	void	g_fprint_random_velocity_inlet_data(FILE*,INTERFACE*,
						    BOUNDARY_STATE*);
IMPORT	void	g_show_wave_states(Wave*);
IMPORT	void	print_WSStenData(WSSten*);
IMPORT	void	set_default_data_file_names(const char*,const char*,
					    char**,char**,INIT_DATA*);
IMPORT	void	g_print_time_dependent_boundary_state(FILE*,INTERFACE*,
                                                      BOUNDARY_STATE*);
IMPORT	void	g_print_extreme_values(FILE*,CHART*,Printplot*);
IMPORT	void	print_graph_footer(FILE*,const char*,boolean);
IMPORT	void	print_graph_header(FILE*,const char*,const char*,
				   boolean,const Front*);
IMPORT  Gas_param*  param_for_comp(COMPONENT);
IMPORT  void        communicate_comp_params(Front*, int);
/*#bjet2 */
IMPORT	void test_print_Tan_stencil(Tan_stencil*);

/*	gprt/gprstate.c */
IMPORT	const char *g_wave_type_as_string(int);
IMPORT	const char *rsoln_wave_name(RIEMANN_SOLVER_WAVE_TYPE);
IMPORT	const char *state_type_name(int);
IMPORT	int	g_read_hsbdry_type_from_string(const char*,INTERFACE*);
IMPORT	int	g_read_state_type_from_string(const char*);
IMPORT	int	g_read_wave_type_from_string(const char*);
IMPORT	void	fprint_gas_data(FILE*,Locstate);
IMPORT	void	fprint_raw_gas_data(FILE*,Locstate,int);
IMPORT	void	fprint_state_type(FILE*,const char*,int);
IMPORT	void	g_fprint_hsbdry_type(FILE*,const char*,int,const char*,
				     INTERFACE*);
IMPORT	void	g_fprint_intfc_state(FILE*,Locstate,INTERFACE*);
IMPORT	void	g_fprint_state(FILE*,Locstate);
IMPORT	void	g_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	g_print_state(Locstate);
IMPORT	void	g_verbose_fprint_intfc_state(FILE*,Locstate,INTERFACE*);
IMPORT	void	g_verbose_print_state(Locstate);
IMPORT	void	print_RP_node_states(const char*,double*,RP_DATA*,int);
IMPORT	void	fprint_curve_status(FILE*,const char*,int);
IMPORT  void    g_fprint_tdp_boundary_state_data(FILE*,INTERFACE*,
				BOUNDARY_STATE*);
#if defined(DEBUG_NODE_PROPAGATE)
IMPORT	void	print_diffraction_status(const char*,int);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
IMPORT	void	print_rsoln_wave(const char*,RIEMANN_SOLVER_WAVE_TYPE,
                                 const char*);
IMPORT	void	print_state_type(const char*,int);
IMPORT	void	g_verbose_fprint_state(FILE*,const char*,Locstate);
IMPORT	int	read_curve_status_from_string(const char*);
IMPORT	boolean	out_spectral_analysis(char*,Wave*,Front*);
IMPORT	boolean	output_spectral_in_time(char*,Wave*,Front*);

/*	grmdata.c */
IMPORT	GraphUnits	*RmGraphUnits(OUTPUT_DATA**);
IMPORT	OUTPUT_DATA	**rm_alloc_output_datas(int);
IMPORT	void	record_rm_amp_and_vel_data(Grid*,Wave*,Front*,Printplot*,
					   OUTPUT_DATA*,boolean);
IMPORT	void	record_radial_rm_amp_and_vel_data(Grid*,Wave*,Front*,Printplot*,
						  OUTPUT_DATA*,boolean);
IMPORT	void	set_rm_layer_indices(OUTPUT_DATA**,int);

/*	gdriverdstat.c */
IMPORT	double	g_quad_integral(BILINEAR_ELEMENT*,TRI_SOLN*,Locstate,POINTER);
IMPORT	double	g_tri_integral(LINEAR_ELEMENT*,TRI_SOLN*,Locstate,POINTER);
IMPORT	void	g_init_grid_statistics(INIT_DATA*,Front*,
				       Grid*,Wave*,Printplot*);

/*	g2dprint.c */
IMPORT	void	fprint_RP_DATA(FILE*,RP_DATA*,double*);
IMPORT	void	g_fprint_ContactWallNodeParams(FILE*,INTERFACE*);
IMPORT	void	g_fprint_RP_DATA_at_nodes(FILE*,INTERFACE*);
IMPORT	void	g_fgraph_front_states(FILE*,Front*);
IMPORT	void	multi_bubble_velocities(Grid*,Wave*,Front*,Printplot*,
					OUTPUT_DATA*,boolean);
IMPORT	void	print_cross_node(NODE*,O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
				 O_CURVE*);
IMPORT	void	print_ramp_refl_stats(Grid*,Wave*,Front*,Printplot*,
				      OUTPUT_DATA*,boolean);
IMPORT	void	record_jet_velocity(Grid*,Wave*,Front*,Printplot*,
				    OUTPUT_DATA*,boolean);
IMPORT	void	set_initial_time_elapsed(double);
IMPORT	void	show_contact_states(Grid*,Wave*,Front*,Printplot*,
				    OUTPUT_DATA*,boolean);
IMPORT	void	show_contact_and_dirichlet_bdry_states(Grid*,Wave*,Front*,
						       Printplot*,OUTPUT_DATA*,
						       boolean);
IMPORT	void	show_front_states_along_lower_wall(Grid*,Wave*,Front*,
						   Printplot*,OUTPUT_DATA*,
						   boolean);
IMPORT	void	show_front_states_for_bowshocks(Grid*,Wave*,Front*,Printplot*,
						OUTPUT_DATA*,boolean);
IMPORT	void	show_front_states_for_expanding_shocks(Grid*,Wave*,Front*,
						       Printplot*,
						       OUTPUT_DATA*,boolean);
IMPORT	void	show_front_states_for_ramp_reflections(Grid*,Wave*,Front*,
						       Printplot*,
						       OUTPUT_DATA*,boolean);
IMPORT	void	show_front_states_for_rm_problem(Grid*,Wave*,Front*,Printplot*,
						 OUTPUT_DATA*,boolean);
/*	gprt/gprcur.c */
IMPORT	void	g_fgraph_curve_states(FILE*,CURVE*,Front*,double*);
IMPORT	void	print_header_for_self_similar_states_along_curve(FILE*,char*,
								 int);
IMPORT	void	g_fprint_header_for_graph_curve_states(FILE*,Front*,
						       const char*);
IMPORT	void	print_self_similar_front_states_along_curve(FILE*,CURVE*,
							    Front*,double*,
							    double,double*);
#if defined(DEBUG_NODE_PROPAGATE)
IMPORT	void	verbose_print_bond_states(const char*,BOND*,CURVE*);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
IMPORT	void	verbose_print_curve_states(CURVE*);


/*      gintext.c */
IMPORT  void    init_intfc_extrema(INIT_DATA*,Front*,Grid*,Printplot*);
IMPORT  double   get_mean_initial_interface_height(int,double*);

/*	glayeravg.c */
IMPORT	void	init_layer_stats(INIT_DATA*,Front*,Grid*,Printplot*);

/*	gprt/glpdiff.c */
IMPORT  void    Lp_diff(Grid*,Wave*,Front*,Printplot*,OUTPUT_DATA*,boolean);

/*	gintstat.c */
IMPORT	void	init_intfc_stats(INIT_DATA*,Front*,Grid*,Printplot*);

/*	grectstat.c */
IMPORT	void	init_rect_state_stats(INIT_DATA*,Front*,Grid*,Printplot*);

/*	gprtst.c */
IMPORT  void  tecplot_interior_states(char *, Wave *);
IMPORT  void  test_print_Stencil(Stencil*);
IMPORT  void  print_curve_pointers(SURFACE *surf);
IMPORT  void  fprint_curve_states(FILE *fp, CURVE **c);
IMPORT  void  test_print_curve_states(INTERFACE  *intfc);
IMPORT  void  test_print_closed_curve(INTERFACE  *);
IMPORT  void  fprint_surface_states(FILE *fp, SURFACE **s);
IMPORT  void  fprint_intfc_states(FILE *, INTERFACE *);
IMPORT  void  print_intfc_in_box(INTERFACE *intfc);
IMPORT  void  out_surf(char *fname, INTERFACE *intfc);
IMPORT  void  tecplot_surface_states(const char *, FILE*,SURFACE*);
IMPORT  void  tecplot_interface_states(const char*, INTERFACE*);
IMPORT  void  check_print_nodes(INTERFACE *);

/*	gdiagnostic.c */
IMPORT  void    summarize_front_states3d(Front*);
IMPORT  void    summary_of_curve_states3d(CURVE*);

/*	grgbdata.c */
IMPORT	void	init_moving_body_data(INIT_DATA*,Front*,Grid*,Printplot*);
#endif /* !defined(_GPRTPROTOS_H) */
