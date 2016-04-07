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
*			hprotos.h
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if !defined(_HPROTOS_H)
#define _HPROTOS_H

#include <hyp/hdecs.h>

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif


	/* Hyp EXPORTED Function Declarations*/
	/* hdriver.c*/
IMPORT	int	hyp_scalar_driver(double,double*,Wave*,Front*);
IMPORT	int	hyp_vector_driver(double,double*,Wave*,Front*);

	/* hinit.c*/
IMPORT	boolean	h_read_print_max_wave_speed_info(INIT_DATA*,const IO_TYPE*,
                                                 Wave*);
IMPORT  boolean    h_read_print_max_viscosity_info(INIT_DATA*,const IO_TYPE*,Wave*);
IMPORT	void	h_set_interface_hooks(int,INIT_DATA*);
IMPORT	void	init_point_source_composition(Wave*,
					      void (*)(Wave*,int,Locstate*));
IMPORT	void	 prompt_for_hyp_method(INIT_DATA*);
IMPORT	void	 prompt_for_unsplit_options(INIT_DATA*);
IMPORT	void     set_hyp_hooks(INIT_DATA*);
IMPORT	void	 set_hyperbolic_method(INIT_DATA*,Wave*,
				       int (**)(double,double*,Wave*,Front*));

	/* hnpt.c*/
IMPORT	void	hyp_npt(int,int*,double,double,Wave*,Wave*,Front*,
			Front*,COMPONENT);
IMPORT	void	set_hyp_npt_globals(Wave*);
IMPORT  COMPONENT       Find_rect_comp(int,int*,Wave*);

	/* hprint.c*/
IMPORT	void	h_fprint_max_wave_speed_info(FILE*,Wave*);

/* For parabolic step */
IMPORT	void	h_fprint_max_viscosity_info(FILE*,Wave*);

IMPORT	void	h_print_H_Front_structure(Front*);
IMPORT	void	print_Stencil(Stencil*);
IMPORT	void	print_Wave_structure(Wave*);
IMPORT	void	print_pt_sources(Wave*);
#if defined(TWOD)
IMPORT	void	print_tri_soln(FILE*,Front*,Wave*,TRI_SOLN*,size_t,
			       double (**)(double*,Front*,POINTER,
					 COMPONENT,Locstate));
#endif /* defined(TWOD) */

	/* hpseudo.c */
IMPORT	int            pseudo_unsplit_driver(double,double*,Wave*,Front*);
IMPORT	void	       h_set_unsplit_options(UnsplitStencilOptions*);


	/* hscatter.c*/
IMPORT  boolean    h_iscatter_states(Wave*,Front*,int*,int);
IMPORT  boolean    h_scatter_states(Wave*,Front*,int*,int);
IMPORT	void	SetDefaultHypPPBlockSize(int);
IMPORT	void	SetHypPPBlockSize(size_t);
IMPORT  void    reflect_states_across_domain(int*,int,int,Front*,Wave*);
IMPORT  void    pp_receive_large_data(byte*,size_t,int);
IMPORT  void    pp_send_large_data(byte*,size_t,int);

	/* hsoln.c*/
IMPORT	double	measure_of_linear_element(LINEAR_ELEMENT*,TRI_SOLN*);
IMPORT	int	copy_hyp_solution_function(Wave*,Wave*);
IMPORT	int	hyp_tri_grid_driver(Front*,Wave*,RECT_GRID*);
IMPORT	int	init_hyp_solution_function(Wave*,Front*);
IMPORT	void	evaluate_dirichlet_boundary_state(double*,HYPER_SURF*,Front*,
						  Wave*,Locstate);
IMPORT	void	h_hyp_solution(double*,COMPONENT,HYPER_SURF*,SIDE,Front*,
			       POINTER,Locstate,Locstate);
IMPORT	void	reinit_hyp_solution_function(Wave*,Front*);
IMPORT	void	states_on_bilinear_element(Locstate*,BILINEAR_ELEMENT*,
					   TRI_GRID*);
IMPORT  void    set_tri_soln_struct(TRI_SOLN*,INTERFACE*,TRI_GRID*,size_t,
                                     INTERPOLATORS*,EL_INTEGRALS*,UNSPLIT*);
	/* hsub.c*/
IMPORT	Stencil	*alloc_stencil(int,Front*);
IMPORT	double	h_max_hyp_time_step(Wave*,double*);

/* For parabolic step */
IMPORT	double	h_max_parab_time_step(Wave*,double*);

IMPORT	int	*set_iperm(int,int);
IMPORT	int	is_source_block(Wave*,INTERFACE*,COMPONENT,int*);
IMPORT  int	reflect_pt_about_Nbdry(double*,double*,double*,COMPONENT,
				       HYPER_SURF*,Front*);
IMPORT	void	free_pt_source_interior_vectors(Wave*);
IMPORT	void	h_copy_into_front(Front*,Front*);
IMPORT	void	h_set_default_front_parameters(INIT_DATA*,Front*);
IMPORT	void	set_pt_source_interior_vectors(Wave*);
IMPORT	void	set_sweep_limits(Wave*,int,int*,int*,int*);
IMPORT	void	set_limits_for_open_bdry(Wave*,Front*,int,int*,int*,int*);

	/* hvec.c*/
IMPORT	void	hyp_reg_vec(int,int*,double,double,Wave*,Wave*,
			    Front*,Front*,COMPONENT);

	/* hwave.c*/
IMPORT	h_MaxWaveSpeed	*h_alloc_MaxWaveSpeed(h_MaxWaveSpeed*,Wave*);
IMPORT	void	h_assign_wave_parameters(Wave*,Wave*);
IMPORT	void	h_copy_into_wave(Wave*,Wave*);
IMPORT	void	h_free_wave(Wave*);
IMPORT	void	h_initialize_max_wave_speed(Wave*);
IMPORT	void	h_set_default_wave_parameters(INIT_DATA*,Wave*);
IMPORT	void	h_set_max_wave_speed(int,double,Locstate,double*,Wave*);
IMPORT  void    h_set_max_viscosity(double,Locstate,double*,Wave*);
IMPORT  h_MaxViscosity  *h_alloc_MaxViscosity(h_MaxViscosity*,Wave*);

#if defined(USE_OVERTURE)
       /* hoverture_dricer.c */
IMPORT  int     hyp_amr_split_driver(double,double*,Wave**,Front**,
                   Front**,int,Wv_on_pc**,int,CompositeGrid*,
                   doubleCompositeGridFunction*,
                   doubleCompositeGridFunction*,Overparam*,
                   void(*)(int,int*,double,double,
                    Wave*,Wave*,Front*,Front*,COMPONENT));
IMPORT  int     hyp_patch_split_driver(double,double*,Wave*,Front*,Front*,
                   CompositeGrid*, doubleCompositeGridFunction*,
                   doubleCompositeGridFunction*,
                   void (*)(int,int*,double,double,Wave*,Wave*,
                    Front*,Front*,COMPONENT));
IMPORT  int     reinstore_undistribute_patch(Wave***,Front***,Wave**,
                   Front**,Wv_on_pc**,int,int);
IMPORT  int     reinstore_mini_undistribute_patch(Wave***,Front***,int*,
                   Wave**,Front**,Wv_on_pc**,int,int,int);
IMPORT  TRI_SOLN  *copy_AMR_tri_soln_storage(TRI_SOLN*,Wave*);
IMPORT  boolean    h_scatter_patch_states(Wave**,Front**,Overparam*,
                   Wv_on_pc**,int,int*,int);
IMPORT  void    patch_wave_trans(Wave**,Wave**,Wv_on_pc**,int,int,int,int);
IMPORT  void    free_copy_patches(Wave**,Front**,Wv_on_pc**,int);
#endif /* if defined(USE_OVERTURE) */

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif


#endif /* !defined(_HPROTOS_H) */
