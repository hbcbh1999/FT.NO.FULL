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
*				ghypprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GHYPPROTOS_H)
#define _GHYPPROTOS_H

#if defined(cray)
#define artvsc ARTVSC
#define lwvec1 LWVEC1
#define lwvec2 LWVEC2
#endif /* defined(cray) */

	/* ghyp EXPORTED Function Declarations */
#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
/*	gflw.F */
FORTRAN	void FORTRAN_NAME(artvsc)(double*,double*,double*,double*,double*,double*,
				  double*,double*,double*,double*,double*,double*,
				  double*,int*);
FORTRAN	void	FORTRAN_NAME(lwvec1)(double*,double*,double*,double*,double*,double*,
				     double*,double*,double*,double*,double*,double*,
				     double*,double*,double*,double*,double*,double*,
				     double*,double*,boolean*,int*);
FORTRAN	void	FORTRAN_NAME(lwvec2)(double*,double*,double*,double*,double*,double*,
				     double*,double*,double*,double*,double*,double*,
				     double*,double*,double*,double*,double*,double*,
				     int*,double*,double*,boolean*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

/*	ggodunov.c */
IMPORT	void	godunov(double,double,Locstate,const double*,int,
                        int*,int*,Stencil*);
IMPORT	void	godunovobl(double,double,Tan_stencil*,Locstate,Front*);

/*	glf.c */
IMPORT	void	LF(double,double,Locstate,const double*,int,int*,int*,Stencil*);
IMPORT	void	LFoblique(double,double,Tan_stencil*,Locstate,Front*);

/*	ghyperbolic.c */
IMPORT	TIME_DYNAMICS g_load_state_vectors(int,int*,Vec_Gas*,int,int,
				           Wave*,Wave*,int*,int);
IMPORT	boolean	check_ans(const char*,double,double,Locstate,
			  COMPONENT,Stencil*,int);
IMPORT	boolean	check_gas(const char*,Locstate*,Locstate,Tan_stencil*,
	                  int,Front*);
IMPORT	void	g_two_side_npt_tang_solver(double,double,Tan_stencil*,
                                           Locstate,Locstate,Front*);
IMPORT	void	point_FD(double,double,Locstate,const double*,
                         int,int*,int*,Stencil*);
IMPORT	void	set_rotation(double**,const double*,int);
IMPORT  boolean    g_detect_and_load_mix_state(int,Stencil*,int);

/*	ghypvec.c */
IMPORT	void	vector_FD(int,int*,double*,Wave*,Wave*,Front*,Front*,
                         int*,int,int,double,double);
IMPORT	void	copy_states_to_new_time_level(int*,Wave*,Wave*,
					      int,int,int*,int);
IMPORT  void    check_and_correct_bad_state(int,int*,int*,Front*,Front*,Wave*,
                                Wave*,double,double,const double*,int,int,int);

/*	glw.c */
IMPORT	void	LW_alloc_phys_vecs(Wave*,int);
IMPORT	void	LW_free_phys_vecs(Wave*);
IMPORT	void    LW(double,double,Locstate,const double*,int,int*,int*,Stencil*);
IMPORT	void	LWoblique(double,double,Tan_stencil*,Locstate,Front*);
IMPORT	void	LW_vec(int,int*,double*,Wave*,Wave*,Front*,Front*,
                       int*,int,int,double,double);

/*	gtvd.c */
IMPORT	void	oned_TVD(int,int*,int*,Wave*,Wave*,Front*,Front*,Stencil*,
                           Tan_stencil*,int,int,Vec_Gas*,Vec_Src*,
			   double,double,int);

/*	gmuscl.c */
IMPORT	void	oned_MUSCL(int,int*,int*,Wave*,Wave*,Front*,Front*,Stencil*,
                           Tan_stencil*,int,int,Vec_Gas*,Vec_Src*,
			   double,double,int);
IMPORT	void	set_muscl_default_hooks(Muscl_Opts*,Hyp_method*,int);
IMPORT	void	set_muscl_default_opts(Muscl_Opts*,Hyp_method*,int);

/*	gplm.c */
IMPORT	void	oned_PLM(int,int*,int*,Wave*,Wave*,Front*,Front*,Stencil*,
                         Tan_stencil*,int,int,Vec_Gas*,Vec_Src*,
			 double,double,int);
IMPORT	void	set_plm_default_opts(Muscl_Opts*,Hyp_method*,int);

/*	ghypprt.c */
IMPORT	void	g_print_internal_energy(const char*,double**,Vec_Muscl*,int,int);
IMPORT  void    g_print_state_vectors(int,int,const char**,
                                      double*,double*,double**,int);
IMPORT	void	g_printout_vec_data(const char*,double*,double*,double**,
				    int,int,int,const char*);

/*	grsolve.c */
IMPORT	void	g_exact_Riemann_midstate(int,int,Vec_Gas*,Vec_Gas*,
                                         double*,double*,Vec_Muscl*);
IMPORT	void	g_gamma_law_fit_Riemann_midstate(int,int,Vec_Gas*,Vec_Gas*,
                                                 double*,double*,Vec_Muscl*);
IMPORT	void	g_gamma_law_fit_rsolver(int,int,double**,Vec_Gas*,double**,
                                        Vec_Gas*,double**,Vec_Gas*,
				        MUSCL_FLUX*,Vec_Muscl*);
IMPORT	void	g_linear_us_up_Riemann_midstate(int,int,Vec_Gas*,Vec_Gas*,
                                                double*,double*,Vec_Muscl*);
IMPORT	void	g_linear_us_up_rsolver(int,int,double**,Vec_Gas*,double**,
                                       Vec_Gas*,double**,Vec_Gas*,
				       MUSCL_FLUX*,Vec_Muscl*);
IMPORT	void	g_muscl_exact_rsolver(int,int,double**,Vec_Gas*,double**,Vec_Gas*,
                                      double**,Vec_Gas*,MUSCL_FLUX*,Vec_Muscl*);

/*	gcgrsolve.c */
IMPORT	void	g_cg_Riemann_midstate(int,int,Vec_Gas*,Vec_Gas*,double*,double*,
				      Vec_Muscl*);
IMPORT	void	cg_rsolve(int,int,double**,Vec_Gas*,double**,Vec_Gas*,
                          double**,Vec_Gas*,MUSCL_FLUX*,Vec_Muscl*);
IMPORT	void	set_cg_params(CG_PARAMS*);

/*	ghypsub.c */
IMPORT	Muscl_Opts *muscl_options(void);
IMPORT	Vec_Gas	   *g_alloc_vgas(Vec_Gas*,int,int);
IMPORT	Vec_Muscl  *alloc_Vec_Muscl(Vec_Muscl*);
IMPORT	Vec_Muscl  *g_muscl_free_wk_space(Vec_Muscl*);
IMPORT	Vec_Muscl  *load_state_data(int,int*,Front*,Wave*,Stencil*,
                                    Tan_stencil*,Vec_Gas*,Vec_Src*,
				    int,int,int,double,double);
IMPORT	Vec_Src    *muscl_alloc_vsrc(Vec_Src*,int,Wave*);
IMPORT	void	oblique_FD(double,double,Tan_stencil*,Locstate,Front*);
IMPORT	void	g_assign_wave_state_vectors(int,int*,Wave*,Wave*,Vec_Gas*,
                                            int,int,int*,int);
IMPORT	void	g_muscl_alloc_phys_vecs(Wave*,int);
IMPORT	void	g_muscl_free_phys_vecs(Wave*);
IMPORT	void	clear_Vec_Gas_set_flags(Vec_Gas*);
IMPORT	void	compute_slope_limiting_coeffs(int,int,Vec_Muscl*);
IMPORT	void	copy_vec_state_params(Vec_Gas*,Vec_Gas*,int,int);
IMPORT	void	g_free_vgas(Vec_Gas*);
IMPORT	void	g_free_vsrc(Vec_Src*);
IMPORT	void	g_load_VGas_state_vectors(int,int,Vec_Gas*,int);
IMPORT	void	g_load_muscl_flux(int,int,double**,Vec_Gas*,
                                  MUSCL_FLUX*,Vec_Muscl*);

/*	gmoc.c */
IMPORT	void	  g_implicit_characteristic_solve(int,int,Vec_Muscl*);
IMPORT	void	  g_riemann_characteristic_solve(int,int,Vec_Muscl*);
IMPORT	void	  g_first_order_direct_characteristic_solve(int,int,Vec_Muscl*);

/*      gvisc.c */
IMPORT	int	parab_driver(double,double*,Wave*,Front*);
IMPORT  void    parab_npt(double,Front*,Wave*,Wave*);
IMPORT  void    parab_tan_solver2d(double,double,Tan_stencil*,
                                           Locstate,Locstate,Front*);
IMPORT  void    parab_tan_solver3d(Front*,const Tparams*,const Tparams*,POINT*);
#if defined(SUBGRID)
IMPORT  void    SGS(double, Front*, Wave*);
#endif /* defined(SUBGRID) */

#endif /* !defined(_GHYPPROTOS_H) */
