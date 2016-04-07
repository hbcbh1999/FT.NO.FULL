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
*				geosprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GEOSPROTOS_H)
#define _GEOSPROTOS_H

#if defined(cray)
#define khstate KHSTATE
#define speed SPEED
#endif /* defined(cray) */

	/* geos EXPORTED Function Declarations*/

/*	geos/eosplot.c */
IMPORT	int	eosp_main(int,char**,INIT_DATA*,INIT_PHYSICS*);

/*	geos/generic-eos.c */
IMPORT	EOS	*set_GENERIC_eos(EOS*);

/*	geos/mpoly-eos.c */
IMPORT	EOS	*set_MPOLY_eos(EOS*);

/*	geos/poly-eos.c */
IMPORT	EOS	*set_POLY_eos(EOS*);

/*	geos/sesame-eos.c */
IMPORT	EOS	*set_SESAME_eos(INIT_DATA*,EOS*);

/*	geos/spoly-eos.c */
IMPORT	EOS	*set_SPOLY_eos(EOS*);
IMPORT	void	set_SPOLY_coefs(Gas_param*);

/*	geos/stellar-eos.c */
IMPORT	EOS	*set_STELLAR_eos(EOS*,char *);
IMPORT	void	set_STELLAR_coefs(Gas_param*);

/*    geos/gentest-eos.c */
IMPORT	EOS	*set_GENTEST_eos(EOS*);

/*	geos/jwl-eos.c */
IMPORT	EOS	*set_JWL_eos(EOS*);

/*	geos/mg-eos.c */
IMPORT	EOS	*set_MG_eos(EOS*);

/*	geos/s2phase-eos.c */
IMPORT	EOS	*set_S2PHASE_eos(EOS*);

/*	geos/geosutils.c */
IMPORT	void	load_pressure(Vec_Gas*,int,int);
IMPORT	void	load_pressure_and_gammas(Vec_Gas*,int,int);
IMPORT	void	load_pressure_and_sound_speed(Vec_Gas*,int,int);
IMPORT	void	load_specific_internal_energy(Vec_Gas*,int,int);
IMPORT	void	load_internal_energy_density(Vec_Gas*,int,int);
IMPORT	void	load_sound_speed(Vec_Gas*,int,int);
IMPORT	void	set_params_jumps(Vec_Gas*,int,int);

/*	geos/giniteos.c */
IMPORT	Gas_param *read_print_EOS_data(INIT_DATA*,const IO_TYPE*,Gas_param*);
IMPORT	int	  prompt_for_equation_of_state_type(const char*,const char*,
						    char**,char**,
						    Prompt_type*,int);
IMPORT	void	g_prompt_for_equation_of_state(INIT_DATA*,struct _Gas_param**,
					       const char*,const char*,
					       INIT_PHYSICS*);

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
/*	gpertsub.F */
    FORTRAN	void	FORTRAN_NAME(khstate)(double*,double*,double*,double*,
                                              double*,double*,double*,double*,
					      double*,double*,double*,double*,
					      double*,double*,double*,double*,
					      double*,double*);
    FORTRAN	void	FORTRAN_NAME(speed)(double*,double*,double*,double*,double*,
                                            double*,double*,double*,double*,double*,
					    double*,double*,double*,double*);
/*    helm_eos.f   */
    FORTRAN     void    FORTRAN_NAME(helmstate)(double*,double*,double*,double*,
                                                double*,double*,double*,double*,
                                                double*,double*,double*,double*,
                                                double*,int*,double*,double*,
                                                double*,int*,char*,int*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */


#if defined(TWOD) && defined(PHASE_CODE)
/*	geos/gphriem.c */
IMPORT	int	intrsct_wv_crv_wth_phs_bdry(Locstate,Locstate,double*,double*,
					    double*,int*,double*,double*,double*,
					    int*,RIEMANN_SOLVER_WAVE_TYPE,
					    RIEMANN_SOLVER_WAVE_TYPE);
IMPORT	int	is_retrograde_bndry(double,Locstate);
IMPORT	void	get_ph_sound_spd(double*,double*,Locstate);
IMPORT	void	make_wave_crv(int,Locstate,double*,double*,double*,int,
			      WAVE_CURVE*);
IMPORT	void	state_on_comp(double,double*,double*,double*,double*,WAVE_CURVE*);
#endif /* defined(TWOD) && defined(PHASE_CODE) */

#endif /* !defined(_GEOSPROTOS_H) */
