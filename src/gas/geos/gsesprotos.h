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
*				gsesprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
 */

#if !defined(_GSESPROTOS_H)
#define _GSESPROTOS_H

#if defined(SESAME_CODE) && defined(TWOD)

#if defined(cray)
#define	cllib		CLLIB
#define	oplib_f		OPLIB_F
#define	oplib_uf	OPLIB_UF
#define	oplib_f		OPLIB_F
#define	oplib_uf	OPLIB_UF
#define	s4get		S4GET
#define	s2get		S2GET
#define	s2eos		S2EOS
#endif /* defined(cray) */

#include <geos/sesame.h>

	/* geos EXPORTED Function Declarations */

/*	geos/gseshyp.c */
IMPORT	const SESAME_TABLE_TYPE *sesame_table_numbers(void);
IMPORT	const char **sesame_table_names(void);
IMPORT	const char **sesame_headers(void);
IMPORT	const size_t *sesame_size_of_states(void);
IMPORT	int	init_ses_hyp_soln_func(Wave*,Front*);
IMPORT	void	init_sesame_hyp_tri_solns(SESAME_EOS*,POINTER,POINTER);
IMPORT	void	set_default_ses_wave_and_front(INIT_DATA*,Front*,Wave*,
					       size_t,boolean);
IMPORT	void	set_ses_hooks(Front*,SESAME_TABLE_TYPE);
IMPORT	void	set_ses_inv_hyp_solns(SESAME_EOS*);
IMPORT	void	set_user_hooks_for_sesame(void);

/*	geos/gsesintrp.c */
IMPORT	boolean	ses_tri_interpolator(double*,LINEAR_ELEMENT*,TRI_SOLN*,
				     Locstate);
IMPORT	boolean	ses_tri_lin_comb_states(double,double,double,double*,Locstate,
					double*,Locstate,double*,
					Locstate,RECT_GRID*,Locstate);
IMPORT	void	ses_lin_comb_states(double,double,double*,Locstate,double*,
				    Locstate,RECT_GRID*,Locstate);
IMPORT	void	ses_quad_interpolator(double*,BILINEAR_ELEMENT*,TRI_SOLN*,
				      Locstate);
IMPORT	void	set_ses_intrp_flag(int,SESAME_TABLE_TYPE);
IMPORT	void	set_ses_intrp_flag_all(SESAME_TABLE_TYPE);
#if defined(PHASE_CODE)

/*	geos/gsesphase.c */
IMPORT	int	get_phase_hyp_state(double,SESAME_EOS*,double*,double*,double*);
IMPORT	int	get_phase_state(double,SESAME_EOS*,PHASE_BDRY*,
				double*,double*,double*);
IMPORT	void	cold_PE_spline(SESAME_EOS*,COLD_CURVE*);
IMPORT	void	get_phase_temp_state(double,Front*,double*,
				     double*,double*,int*,SESAME_EOS*);
IMPORT	void	init_RT_interior_states(Wave*,Front*,SESAME_EOS*,
					COLD_CURVE*,PHASE_BDRY*);
IMPORT	void	init_new_phase_bound(Wave*,Front*,SESAME_EOS*,COLD_CURVE*,
				     PHASE_BDRY*);
IMPORT	void	lookspl(SESAME_EOS*,PHASE_BDRY*,double,int,int,int,int,
			double*,double*,double*,double*);
IMPORT	void	phase_spline(SESAME_EOS*,COLD_CURVE*,PHASE_BDRY*);
IMPORT	void	ses_phase_states(Front*,SESAME_EOS*,PHASE_BDRY*);
IMPORT	void	set_boundary_states(Front*,SESAME_EOS*,PHASE_BDRY*,
				    COLD_CURVE*);

/*	geos/gsesspline.c */
IMPORT	void	get_temp_state(double,int,SESAME_EOS*,PHASE_BDRY*,double*,double*,
			       double*,int*);
IMPORT	void	init_PE_phase_grid(SESAME_EOS*,double,PHASE_BDRY*,int*);
IMPORT	void	init_PE_spline(SESAME_EOS*,PHASE_BDRY*);
IMPORT	void	set_cross_states(Front*,SESAME_EOS*,PHASE_BDRY*,double*,
				 double*,double*,double*,double*);
#endif /* defined(PHASE_CODE) */

/*	geos/gsesprint.c */
IMPORT	void	fprint_SESAME_params(FILE*,SESAME_EOS*);
IMPORT	void	print_ps_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_re_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_rt_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_rs_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_title_for_sesame(FILE*,SESAME_EOS*);
IMPORT	void	print_vp_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	read_print_SESAME_params(SESAME_EOS*,const IO_TYPE*);
IMPORT	void	ses_rt_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_re_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_rs_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_ps_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_vp_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	verbose_ses_show_intfc_states(INTERFACE*,SESAME_TABLE_TYPE,
					      SESAME_EOS*);

/*	geos/gsesinout.c */
IMPORT	boolean	restart_sesame(INIT_DATA*,SESAME_EOS*);
IMPORT	void	init_ses_prt(Printplot*,Front*,int);
IMPORT	void	ses_printout(CHART*,Printplot*,boolean,int);

/*	geos/sesinv.c */
IMPORT	boolean	zero_temperature_cold_curve(SESAME_EOS*);
IMPORT	double	sesame_rt_sound_speed_squared(double,double,double*,double*);
IMPORT	void	s2eos_lookup(double,double,double*,double*,SESAME_EOS*);
IMPORT	void	set_RT_entropy_from_cold_curve(Front*,Wave*,SESAME_EOS*);
IMPORT	void	set_RT_entropy_from_mid_point(Front*,Wave*,SESAME_EOS*);
IMPORT	void	setrt(double,double,Locstate,SESAME_EOS*);

/*	geos/sesspln.c */
IMPORT	int	splcomp(double*,double*,double*,int,double*,double*,
				      double*,double*,double*);
IMPORT	int	splcomp2(double*,double*,double*,int,double*,double*,double*,double*,
			 double*,double,double);
IMPORT	void	spline(double,double,double,double,double,double,double,double*,double*);

/*	geos/sesstate.c */
IMPORT	void	phbnd(double*,double,double,double,SESAME_EOS*);
IMPORT	void	sets(double,double*,int,int*,SESAME_EOS*);
IMPORT	void	setspb(double,double*,int,double,SESAME_EOS*);

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
/*	geos/gsestoc.F */
    FORTRAN	void	FORTRAN_NAME(cllib)(int*);
#   if defined(cray)
    FORTRAN	void	SFORTRAN_NAME(oplib_f)(int*,_fcd);
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,_fcd);
#   else /* defined(cray) */
    FORTRAN	void	SFORTRAN_NAME(oplib_f)(int*,const char*,int);
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,const char*,int);
#   endif /* defined(cray) */

/*	geos/sesadd.F */
    FORTRAN	void	FORTRAN_NAME(s4get)(int*,int*,double*,int*,int*,int*);

/*	geos/sesame.F */
    FORTRAN	void	FORTRAN_NAME(s2get)(int*,int*,double*,int*,int*,int*);
    FORTRAN	void	FORTRAN_NAME(s2eos)(int*,double*,double*,double*,
                                            double*,double*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

#endif /* defined(SESAME_CODE) && defined(TWOD) */
#endif /* !defined(_GSESPROTOS_H) */
