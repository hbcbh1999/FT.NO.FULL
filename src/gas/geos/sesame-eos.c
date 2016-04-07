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
*
*				sesame-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*	NOTE:  The data contained in the table tbls,  returned from
*	a call to s2get() is arranged in the following way.
*
*	tbls[0] = material index
*	tbls[1] = reference density
*	nr = tbls[2] = number of density points in the table
*	nt = tbls[3] = number of temperature points in the table
*	tbls[4]...tbls[3+nr] = density values
*	tbls[4+nr]...tbls[3+nr+nt] = temperature values
*	tbls[4+nr+nt]...tbls[3+nr+nt+nr*nt] = pressure values
*	tbls[4+nr+nt+nr*nt]...tbls[3+nr+nt+2*nr*nt] = energy values
*
*	The pressure and energy at index points i, j where
*
*	rho = tbls[4+i], 0 <= i < nr, T = tbls[4+nr+j], 0 <= j < nt
*
*	are given by
*
*	p = tbls[4+nr+nt       + nr*j+i],
*	e = tbls[4+nr+nt+nr*nt + nr*j+i]
*/

#if defined(SESAME_CODE) && defined(TWOD)

#define DEBUG_STRING    "ses_hyp"

#include <geos/sesame.h>

	/* LOCAL function prototypes */
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL	double	SESAME_internal_energy(Locstate);
LOCAL	double	SESAME_pressure(Locstate);
LOCAL	double	SESAME_sound_speed_squared(Locstate);
LOCAL	double	SESAME_acoustic_impedance_squared(Locstate);
LOCAL	double	SESAME_specific_internal_energy(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL	double	SESAME_temperature(Locstate);
LOCAL	double	SESAME_entropy(Locstate);
LOCAL	double	SESAME_adiabatic_gamma(Locstate);
LOCAL	double	SESAME_gruneisen_gamma(Locstate);
LOCAL	double	SESAME_C_V(Locstate);

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
LOCAL	void	SESAME_single_eos_load_pressure_and_sound_speed2(Vec_Gas*,
								 int,int);
LOCAL	void	SESAME_single_eos_load_pressure_and_gammas(Vec_Gas*,int,int);
LOCAL	void	SESAME_single_eos_load_pressure(Vec_Gas*,int,int);
LOCAL	void	SESAME_single_eos_load_sound_speed2(Vec_Gas*,int,int);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */

	/* Velocity Related Hugoniot Functions */

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL	double	SESAME_dens_rarefaction(double,Locstate);
LOCAL	double	SESAME_pressure_rarefaction(double,Locstate);
LOCAL	void	SESAME_state_on_adiabat_with_pr(Locstate,double,Locstate,int);
LOCAL	void	SESAME_state_on_adiabat_with_dens(Locstate,double,Locstate,int);

	/* General Wave Curve Functions */
LOCAL	double	SESAME_mass_flux(double,Locstate);
LOCAL	double	SESAME_mass_flux_squared(double,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
LOCAL	double	SESAME_oned_fan_state(double,Locstate,Locstate,Locstate,
				      int,boolean*);

	/* Functions to Compute Riemann Solutions */
LOCAL	double	SESAME_riemann_wave_curve(Locstate,double);
LOCAL	double	SESAME_eps_for_Godunov(Locstate,double,double);
LOCAL	void	SESAME_initialize_riemann_solver(Locstate,Locstate,double*,
						 double*,double,double*,double*,
						 boolean(*)(Locstate,Locstate,
							 double,double*,double*,
							 double*,double*,double*,
							 double*,
							 RIEMANN_SOLVER_WAVE_TYPE*,
							 RIEMANN_SOLVER_WAVE_TYPE*));


	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	SESAME_prompt_for_state(Locstate,int,Gas_param*,const char*);
LOCAL	void	SESAME_prompt_for_thermodynamics(Locstate,Gas_param*,
						 const char*);
LOCAL	void	SESAME_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	SESAME_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                             Gas_param*);
LOCAL	EOS*	SESAME_free_EOS_params(EOS*);
LOCAL	void	SESAME_prompt_for_EOS_params(INIT_DATA*,
                                             Gas_param*,const char*,
					     const char*);

	/* Problem Type Specific Initialization Functions */
LOCAL	void	set_eos_function_hooks(EOS*);

	/* Local support function prototypes */
LOCAL	boolean	invert_pe(double,double*,POINTER);
LOCAL	boolean	initialize_Sesame_tables(SESAME_EOS*);
LOCAL	double	int_dp_over_c_rho(double,Locstate);
#if defined(PHASE_CODE)
LOCAL	double	mass_flux_for_comp(double,Locstate);
LOCAL	double	mass_flux_for_phase_change(double,Locstate);
#endif /* defined(PHASE_CODE) */
#if  defined(RAREFACTION_STATE_BY_PRESSURE)
LOCAL	boolean	ses_p_cprf(double,double*,POINTER);
#else /* defined(RAREFACTION_STATE_BY_PRESSURE) */
LOCAL	boolean	ses_r_cprf(double,double*,POINTER);
#endif /* defined(RAREFACTION_STATE_BY_PRESSURE) */
LOCAL	void	sesinv(SESAME_EOS*);
LOCAL	void	open_fortran_device(char*,int);
LOCAL	void	print_list_of_sesame_materials(FILE*,char*);
LOCAL	void	prompt_for_window_params(SESAME_EOS*);

EXPORT	EOS	*set_SESAME_eos(
	INIT_DATA *init,
	EOS	  *eos)
{
	int i;
	if (eos == NULL)
	{
	    SESAME_EOS *seos;
	    const char **ses_table_names;
	    const size_t *sizests;

	    ses_table_names = sesame_table_names();
	    sizests = sesame_size_of_states();

	    scalar(&seos,sizeof(SESAME_EOS));
	    seos->_ses_table_name = ses_table_names;
	    seos->num_ses_tabs = NUMBER_SESAME_TABLES;
	    uni_array(&seos->root,seos->num_ses_tabs,sizeof(CHART *));
	    uni_array(&seos->fr,seos->num_ses_tabs,sizeof(Front *));
	    uni_array(&seos->wave,seos->num_ses_tabs,sizeof(Wave *));
	    for (i = 0; i < seos->num_ses_tabs; ++i)
	    {
	    	RECT_GRID	*gr;
	    	Front		*fr;
	    	Wave		*wave;
	    	CHART		*root;

	    	scalar(&fr,sizeof(Ses_Front));
	    	scalar(&wave,sizeof(Ses_Wave));
	    	scalar(&root,sizeof(CHART));
	    	scalar(&gr,sizeof(RECT_GRID));

	    	Ses_front_seos(fr) = Ses_wave_seos(seos) = seos;
	    	root->front = fr;
	    	root->wave = wave;
	    	fr->rect_grid = wave->rect_grid = gr;
	    	seos->root[i] = root;
	    	seos->fr[i] = fr;
	    	seos->wave[i] = wave;
	    	set_default_ses_wave_and_front(init,fr,wave,sizests[i],YES);
	    }
	    seos->BIS_EPS =	0.1; /* TOLERANCE */
	    seos->ABS_SES_EPS = 1.0e-07; /* TOLERANCE */
	    eos = &seos->eos;
	}
	(void) set_GENERIC_eos(eos);
	set_eos_function_hooks(eos);
	return eos;
}		/*end set_SESAME_eos*/

LOCAL	void	set_eos_function_hooks(
	EOS		*eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_internal_energy = SESAME_internal_energy;
	eos->_pressure = SESAME_pressure;
	eos->_sound_speed_squared = SESAME_sound_speed_squared;
	eos->_acoustic_impedance_squared = SESAME_acoustic_impedance_squared;
	eos->_specific_internal_energy = SESAME_specific_internal_energy;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_temperature = SESAME_temperature;
	eos->_entropy = SESAME_entropy;
	eos->_adiabatic_gamma = SESAME_adiabatic_gamma;
	eos->_gruneisen_gamma = SESAME_gruneisen_gamma;
	eos->_C_V = SESAME_C_V;

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
	eos->_single_eos_load_pressure_and_sound_speed2 =
		SESAME_single_eos_load_pressure_and_sound_speed2;
	eos->_single_eos_load_pressure_and_gammas =
		SESAME_single_eos_load_pressure_and_gammas;
	eos->_single_eos_load_pressure = SESAME_single_eos_load_pressure;
	eos->_single_eos_load_sound_speed2 =
	    SESAME_single_eos_load_sound_speed2;

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */

	/* Velocity Related Hugoniot Functions */

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	eos->_dens_rarefaction = SESAME_dens_rarefaction;
	eos->_pressure_rarefaction = SESAME_pressure_rarefaction;
	eos->_state_on_adiabat_with_pr = SESAME_state_on_adiabat_with_pr;
	eos->_state_on_adiabat_with_dens = SESAME_state_on_adiabat_with_dens;

	/* General Wave Curve Functions */
	eos->_mass_flux = SESAME_mass_flux;
	eos->_mass_flux_squared = SESAME_mass_flux_squared;

	/* Functions for the Evaluation of Riemann Solutions */
	eos->_oned_fan_state = SESAME_oned_fan_state;

	/* Functions to Compute Riemann Solutions */
	eos->_riemann_wave_curve = SESAME_riemann_wave_curve;
	eos->_eps_for_Godunov = SESAME_eps_for_Godunov;
	eos->_initialize_riemann_solver = SESAME_initialize_riemann_solver;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_prompt_for_state = SESAME_prompt_for_state;
	eos->_prompt_for_thermodynamics = SESAME_prompt_for_thermodynamics;
	eos->_fprint_EOS_params = SESAME_fprint_EOS_params;
	eos->_read_print_EOS_params = SESAME_read_print_EOS_params;
	eos->_free_EOS_params = SESAME_free_EOS_params;
	eos->_prompt_for_EOS_params = SESAME_prompt_for_EOS_params;
}		/*end set_eos_function_hooks*/


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*			SESAME_internal_energy():
*
*	Returns the internal energy per unit volume of a state.
*/

LOCAL	double	SESAME_internal_energy(
	Locstate state)
{
	Front *fr;
	Wave *wv;
	COMPONENT comp;
	double var[NUM_SES_VAR];
	double coords[MAXD];
	double rho;
	double T;

	switch (state_type(state)) 
	{
	case GAS_STATE:
		return	Energy(state) - kinetic_energy(state);

	case EGAS_STATE:
		return	Energy(state)*Dens(state);

	case VGAS_STATE:
		return Dens(state) * Int_en(state);
		
	case TGAS_STATE:
		return Dens(state)*specific_internal_energy(state);

	case FGAS_STATE:
		fr = Ses_params(state)->fr[SESAME_RHO_TEMP];
		wv = Ses_params(state)->wave[SESAME_RHO_TEMP];
		rho = Dens(state);
		T = Temperature(state);

		set_ses_intrp_flag(EVALUATE_ENERGY,SESAME_RHO_TEMP);
		coords[0] = ses_rt_grid_from_rho(rho,Ses_params(state));
		coords[1] = ses_rt_grid_from_temp(T,Ses_params(state));
		comp = nearest_interior_comp(multiphase_eos(Eos(state)),
				COMP_PURE_PHASE,coords,fr->interf);
		ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
		return rho*(var[RT_CE] + T*var[RT_RE]);

	default:
		screen("ERROR: in SESAME_internal_energy(), "
		       "no such state type\n");
		clean_up(ERROR);
	}
	return ERROR_FLOAT;
}		/*end SESAME_internal_energy*/


/*
*			SESAME_pressure():
*
*	Returns the thermodynamic pressure of a state.
*
*				     dE  |
*			     P = -  ---- |
*		                     dV  |S
*
*	Where E = specific internal energy,  V = specific volume,  and
*	S = specific entropy.
*/

LOCAL	double	SESAME_pressure(
	Locstate state)
{
	Front *fr;
	Wave *wv;
	double  e, T, coords[MAXD], var[NUM_SES_VAR];
	COMPONENT comp;
	double rho, pr;

	if (is_obstacle_state(state))
		return HUGE_VAL;
	rho = Dens(state);
	switch (state_type(state)) 
	{
	case	GAS_STATE:
#if !defined(UNRESTRICTED_THERMODYNAMICS)
		if (rho < Vacuum_dens(state))
			return Min_pressure(state);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
		fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
		wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];

		e = specific_internal_energy(state);
		coords[0] = ses_re_grid_from_rho(rho,Ses_params(state));
		coords[1] = ses_re_grid_from_engy(e,Ses_params(state));
		set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_ENERGY);
		comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
		ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
		T = ses_re_temp_from_var(var[RE_T],Ses_params(state));
		pr = rho*T*var[RE_RP] + var[RE_CP];
		break;

	case	EGAS_STATE:
		fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
		wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
		e = Energy(state);
		coords[0] = ses_re_grid_from_rho(rho,Ses_params(state));
		coords[1] = ses_re_grid_from_engy(e,Ses_params(state));
		set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_ENERGY);
		comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
		ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
		T = ses_re_temp_from_var(var[RE_T],Ses_params(state));
		pr = rho*T*var[RE_RP] + var[RE_CP];
		break;

	case	FGAS_STATE:
		fr = Ses_params(state)->fr[SESAME_RHO_TEMP];
		wv = Ses_params(state)->wave[SESAME_RHO_TEMP];
		T = Temperature(state);
		coords[0] = ses_rt_grid_from_rho(rho,Ses_params(state));
		coords[1] = ses_rt_grid_from_temp(T,Ses_params(state));
		set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_TEMP);
		comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
		ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
		pr = rho*T*var[RT_RP] + var[RT_CP];
		break;

	case	TGAS_STATE:
	case	VGAS_STATE:
		pr = Press(state);
		break;

	default:
		screen("ERROR in SESAME_pressure(), no such state type\n");
		clean_up(ERROR);
		break;
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(state));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end SESAME_pressure*/


/*
*			SESAME_sound_speed_squared():
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

LOCAL	double	SESAME_sound_speed_squared(
	Locstate state)
{
	Front     *fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
	Wave      *wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
	double     rho, en, coords[MAXD];
	double     pr, T, var[NUM_SES_VAR];
	double     c2;
	COMPONENT comp;
		
	rho = Dens(state);
	en = specific_internal_energy(state);
	coords[0] = ses_re_grid_from_rho(rho,Ses_params(state));
	coords[1] = ses_re_grid_from_engy(en,Ses_params(state));
	set_ses_intrp_flag(EVALUATE_PRESSURE | EVALUATE_ADB_GAMMA,
			SESAME_RHO_ENERGY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	T = ses_re_temp_from_var(var[RE_T],Ses_params(state));
	pr = rho*T*var[RE_RP] + var[RE_CP];
	c2 = var[RE_AG]*pr/rho;
	return c2;
}		/*end SESAME_sound_speed_squared*/


/*
*		SESAME_acoustic_impedance_squared():
*
*	Returns the square of the local acoustic impedance of the state.
*
*                        2     dP  |
*			i = - ---- |
*                              dV  |S
*/

LOCAL	double	SESAME_acoustic_impedance_squared(
	Locstate state)
{
	Front *fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
	Wave *wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
	double  rho, en, coords[MAXD];
	double  pr, T, var[NUM_SES_VAR];
	double  i2;
	COMPONENT comp;
		
	rho = Dens(state);
	en = specific_internal_energy(state);
	coords[0] = ses_re_grid_from_rho(rho,Ses_params(state));
	coords[1] = ses_re_grid_from_engy(en,Ses_params(state));
	set_ses_intrp_flag(EVALUATE_PRESSURE | EVALUATE_ADB_GAMMA,
			SESAME_RHO_ENERGY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	T = ses_re_temp_from_var(var[RE_T],Ses_params(state));
	pr = rho*T*var[RE_RP] + var[RE_CP];
	i2 = var[RE_AG]*pr*rho;
	return i2;
}		/*end SESAME_acoustic_impedance_squared*/

/*
*			SESAME_specific_internal_energy():
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/

typedef struct {
	Front		*fr;
	Wave		*wv;
	double		rho;
	double		rho_grid;
	boolean		phase;
	SESAME_EOS	*seos;
} SES_INVERT_PE_PARAMS;

LOCAL	boolean	invert_pe(
	double		e_grid,
	double		*p,
	POINTER		params)
{
	double		T, var[NUM_SES_VAR];

	Front		*fr = ((SES_INVERT_PE_PARAMS *) params)->fr;
	Wave		*wv = ((SES_INVERT_PE_PARAMS *) params)->wv;
	double		rho = ((SES_INVERT_PE_PARAMS *) params)->rho;
	boolean		phase = ((SES_INVERT_PE_PARAMS *) params)->phase;
	SESAME_EOS	*seos = ((SES_INVERT_PE_PARAMS *) params)->seos;
	double		coords[MAXD];
	COMPONENT	comp;

	coords[0] = ((SES_INVERT_PE_PARAMS *) params)->rho_grid;
	coords[1] = e_grid;
	comp = nearest_interior_comp(phase,COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	T = ses_re_temp_from_var(var[RE_T],seos);
	*p = var[RE_CP] + rho*T*var[RE_RP];
	return FUNCTION_SUCCEEDED;
	
}		/*end invert_pe*/

LOCAL	double	SESAME_specific_internal_energy(
	Locstate state)
{
	SES_INVERT_PE_PARAMS	inv_pe_params;
	COMPONENT		comp;
	double			T, rho, pr, e, e_grid, emin_grid, emax_grid;
	double			epsilon, delta;
	double			rel_err;
	double			coords[2];
	double			var[NUM_SES_VAR];
	Front			*fr;
	Wave			*wv;

	switch (state_type(state))
	{

	case	GAS_STATE:
	    return (Energy(state) - kinetic_energy(state))/Dens(state);

	case	EGAS_STATE:
	    return Energy(state);

	case	TGAS_STATE:
	    rho = Dens(state);
	    rel_err = Ses_params(state)->reler0;
	    fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
		
	    pr = pressure(state);
	    inv_pe_params.seos = Ses_params(state);
	    inv_pe_params.phase = multiphase_eos(Eos(state));
	    inv_pe_params.fr = fr;
	    inv_pe_params.wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
	    inv_pe_params.rho = rho;
	    inv_pe_params.rho_grid =
		ses_re_grid_from_rho(rho,Ses_params(state));
		
	    emin_grid = fr->rect_grid->L[1];
	    emax_grid = fr->rect_grid->U[1];
	    epsilon = rel_err*pr;
	    delta = 0.5*rel_err*(emin_grid + emax_grid);
	    set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_ENERGY);
	    if (find_root(invert_pe,(POINTER)&inv_pe_params,pr,&e_grid,
				     emin_grid,emax_grid,epsilon,delta) ==
		FUNCTION_FAILED) 
	    {
	        screen("ERROR in SESAME_specific_internal_energy(), "
	               "Unable to invert SESAME function.\n");
	        (void) printf("rho = %g, p = %g,emin = %g,emax = %g\n",rho,pr,
			ses_re_engy_from_grid(emin_grid,Ses_params(state)),
	        	ses_re_engy_from_grid(emax_grid,Ses_params(state)));
	        (void) printf("epsilon = %g, delta = %g\n",epsilon,delta);
	        clean_up (ERROR);
	    }
	    e = ses_re_engy_from_grid(e_grid,Ses_params(state));
	    return e;

	case	FGAS_STATE:
	    fr = Ses_params(state)->fr[SESAME_RHO_TEMP];
	    wv = Ses_params(state)->wave[SESAME_RHO_TEMP];
	    rho = Dens(state);
	    T = Temperature(state);

	    set_ses_intrp_flag(EVALUATE_ENERGY,SESAME_RHO_TEMP);
	    coords[0] = ses_rt_grid_from_rho(rho,Ses_params(state));
	    coords[1] = ses_rt_grid_from_temp(T,Ses_params(state));
	    comp = nearest_interior_comp(multiphase_eos(Eos(state)),
	    		COMP_PURE_PHASE,coords,fr->interf);
	    ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	    return var[RT_CE] + T*var[RT_RE];
	
	case	VGAS_STATE:
	    return Int_en(state);

	default:
	    screen("ERROR in SESAME_specific_internal_energy(), "
	           "no such state type\n");
	    clean_up(ERROR);
	    break;
	}
	return ERROR_FLOAT;
}		/*end SESAME_specific_internal_energy*/

/***************END PRIMARY THERMODYNAMIC FUNCTIONS ************************/
/***************SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS ***********/


/*
*			SESAME_temperature():
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	double	SESAME_temperature(
	Locstate state)
{
	Front *fr;
	Wave *wv;
	double  en, coords[MAXD], T, var[NUM_SES_VAR];
	COMPONENT comp;

	if (state_type(state) == FGAS_STATE) return Temperature(state);

		
	fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
	wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
	coords[0] = ses_re_grid_from_rho(Dens(state),Ses_params(state));
	en = specific_internal_energy(state);
	coords[1] = ses_re_grid_from_engy(en,Ses_params(state));
	set_ses_intrp_flag(EVALUATE_TEMPERATURE,SESAME_RHO_ENERGY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	T = ses_re_temp_from_var(var[RE_T],Ses_params(state));
	return T;
}		/*end SESAME_temperature*/

/*
*			SESAME_entropy():
*
*	Returns the specific entropy of a state.
*/

LOCAL	double	SESAME_entropy(
	Locstate state)
{
	Front *fr;
	Wave *wv;
	double  coords[MAXD], var[NUM_SES_VAR];
	COMPONENT comp;

	if (state_type(state) == VGAS_STATE)
	    return Entropy(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(state) < Vacuum_dens(state))
	    return 0.0;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
	wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
		
	coords[0] = ses_re_grid_from_rho(Dens(state),Ses_params(state));
	coords[1] = ses_re_grid_from_engy(
			specific_internal_energy(state),Ses_params(state));
	set_ses_intrp_flag(EVALUATE_ENTROPY,SESAME_RHO_ENERGY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state)),
				COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	return var[RE_S];
}		/*end SESAME_entropy*/

/*
*			SESAME_adiabatic_gamma():
*
*	Returns the dimensionless sound speed
*
*		gamma = - d(log P)/d(log V) | .
*					     S
*	As usual P = thermodynamic pressure,  V = specific volume
*	and S = specific entropy.
*/

LOCAL	double	SESAME_adiabatic_gamma(
	Locstate state)
{
	Front *fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
	Wave *wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
	double  rho, en, coords[MAXD], var[NUM_SES_VAR];
	COMPONENT comp;
		
	rho = Dens(state);
	en = specific_internal_energy(state);
	coords[0] = ses_re_grid_from_rho(rho,Ses_params(state));
	coords[1] = ses_re_grid_from_engy(en,Ses_params(state));
	set_ses_intrp_flag(EVALUATE_ADB_GAMMA,SESAME_RHO_ENERGY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	return var[RE_AG];
}		/*end SESAME_adiabatic_gamma*/


/*
*			SESAME_gruneisen_gamma():
*
*	Returns the dimensionless Gruneisen exponent
*
*
*                                                 dP/dE |
*		GAMMA = - d(log T)/d(log V) |  =  -----  V
*                                            S     rho
*
*	As usual P = thermodynamic pressure,  V = specific volume
*	rho = density, E = specific internal energy,
*	and  S = specific entropy.
*
*
*/

LOCAL	double	SESAME_gruneisen_gamma(
	Locstate state)
{
	Front *fr = Ses_params(state)->fr[SESAME_RHO_ENERGY];
	Wave *wv = Ses_params(state)->wave[SESAME_RHO_ENERGY];
	double  rho, en, coords[MAXD], var[NUM_SES_VAR];
	COMPONENT comp;
		
	rho = Dens(state);
	en = specific_internal_energy(state);
	coords[0] = ses_re_grid_from_rho(rho,Ses_params(state));
	coords[1] = ses_re_grid_from_engy(en,Ses_params(state));
	set_ses_intrp_flag(EVALUATE_GRU_GAMMA,SESAME_RHO_ENERGY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	return var[RE_GG];
}		/*end SESAME_gruneisen_gamma*/


/*
*			SESAME_C_V():
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*/

LOCAL	double	SESAME_C_V(
	Locstate state)
{
	double T = temperature(state);
	double rho = Dens(state);
	double pvec[3], evec[3];

	s2eos_lookup(rho,T,pvec,evec,Ses_params(state));
	return evec[2];
}	/* end SESAME_C_V */


/***************END SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS *******/
/***************VECTORIZED THERMODYNAMIC FUNCTIONS *************************/

/*
*		SESAME_single_eos_load_pressure_and_sound_speed2():
*
*	Loads a vector of pressures and sound speeds into the
*	appropriate fields of the Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	SESAME_single_eos_load_pressure_and_sound_speed2(
	Vec_Gas *vst,
	int     offset,
	int     vsize)
{
	Gas_param  *params = Params(vst->state[offset]);
	SESAME_EOS *seos;
	Front      *fr;
	Wave       *wv;
	double      coords[MAXD];
	double      T, var[NUM_SES_VAR];
	double      *rho = vst->rho + offset;
	double      *p = vst->p + offset, *c2 = vst->c2 + offset;
	double      *e = vst->e + offset;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double      *min_pressure = vst->min_pressure + offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int        interp_flag;
	int        k;

	seos = (SESAME_EOS *)params->eos;
	fr = seos->fr[SESAME_RHO_ENERGY];
	wv = seos->wave[SESAME_RHO_ENERGY];
	interp_flag = EVALUATE_PRESSURE | EVALUATE_ADB_GAMMA;
	set_ses_intrp_flag(interp_flag,SESAME_RHO_ENERGY);
		
	for (k = 0; k < vsize; ++k)
	{
	    coords[1] = ses_re_grid_from_engy(e[k],seos);
	    coords[0] = ses_re_grid_from_rho(rho[k],seos);
	    ses_solution(coords,COMP_PURE_PHASE,NULL,POSITIVE_SIDE,fr,wv,var);
	    T = ses_re_temp_from_var(var[RE_T],seos);
	    p[k] = rho[k]*T*var[RE_RP] + var[RE_CP];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (p[k] < min_pressure[k])
		p[k] = min_pressure[k];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    c2[k] = var[RE_AG]*p[k]/rho[k]; 
	}
	if (vst->FD != NULL)
	{
	    Locstate state = vst->state[offset];
	    double    *FD = vst->FD + offset;
	    static   Locstate tmpst = NULL;

	    if (tmpst == NULL)
	    {
	        (*Params(state)->_alloc_state)(&tmpst,Params(state)->sizest);
		set_type_of_state(tmpst,EGAS_STATE);
	    }
	    for (k = 0; k < vsize; ++k)
	    {
	        Dens(tmpst) = rho[k];
		Energy(tmpst) = e[k];
		Set_params(tmpst,state);
	        FD[k] = fundamental_derivative(tmpst);
	    }
	}
}		/*end SESAME_single_eos_load_pressure_and_sound_speed2*/


/*
*		SESAME_single_eos_load_pressure_and_gammas():
*
*	Loads the pressure, adiabatic exponent, and Gruneisen
*	coefficient uni_arrays of the Vec_Gas state vst.
*	This function assumes that the specific internal energy
*	uni_array vst->e is already loaded.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	SESAME_single_eos_load_pressure_and_gammas(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	Gas_param  *params = Params(vst->state[offset]);
	SESAME_EOS *seos;
	Front      *fr;
	Wave       *wv;
	double      coords[MAXD];
	double      T, var[NUM_SES_VAR];
	double      *rho = vst->rho + offset;
	double      *p = vst->p + offset, *e = vst->e + offset;
	double      *c2 = vst->c2 + offset, *GAM = vst->GAM + offset;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double      *min_pressure = vst->min_pressure + offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int        interp_flag;
	int        k;

	seos = (SESAME_EOS *)params->eos;
	fr = seos->fr[SESAME_RHO_ENERGY];
	wv = seos->wave[SESAME_RHO_ENERGY];
	interp_flag = EVALUATE_PRESSURE | EVALUATE_ADB_GAMMA |
				EVALUATE_GRU_GAMMA;
	set_ses_intrp_flag(interp_flag,SESAME_RHO_ENERGY);
		
	for (k = 0; k < vsize; ++k)
	{
	    coords[1] = ses_re_grid_from_engy(e[k],seos);
	    coords[0] = ses_re_grid_from_rho(rho[k],seos);
	    ses_solution(coords,COMP_PURE_PHASE,NULL,POSITIVE_SIDE,fr,wv,var);
	    T = ses_re_temp_from_var(var[RE_T],seos);
	    p[k] = rho[k]*T*var[RE_RP] + var[RE_CP];
	    GAM[k] = var[RE_GG];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (p[k] < min_pressure[k])
	    	p[k] = min_pressure[k];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    c2[k] = var[RE_AG]*p[k]/rho[k];
	}
	if (vst->FD != NULL)
	{
	    Locstate state = vst->state[offset];
	    double    *FD = vst->FD + offset;
	    static   Locstate tmpst = NULL;

	    if (tmpst == NULL)
	    {
	        (*Params(state)->_alloc_state)(&tmpst,Params(state)->sizest);
		set_type_of_state(tmpst,EGAS_STATE);
	    }
	    for (k = 0; k < vsize; ++k)
	    {
	        Dens(tmpst) = rho[k];
		Energy(tmpst) = e[k];
		Set_params(tmpst,state);
	        FD[k] = fundamental_derivative(tmpst);
	    }
	}
}		/*end SESAME_single_eos_load_pressure_and_gammas*/

/*
*			SESAME_single_eos_load_pressure():
*
*	Loads a vector of pressures into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	SESAME_single_eos_load_pressure(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	Gas_param *params = Params(vst->state[offset]);
	Front *fr;
	SESAME_EOS *seos;
	Wave *wv;
	double  coords[MAXD], T, var[NUM_SES_VAR];
	double *p = vst->p + offset;
	double *rho = vst->rho + offset;
	double *e = vst->e + offset;
	int k;

	seos = (SESAME_EOS *)params->eos;
	fr = seos->fr[SESAME_RHO_ENERGY];
	wv = seos->wave[SESAME_RHO_ENERGY];
	for (k = 0; k < vsize; ++k)
	{
	    coords[0] = ses_re_grid_from_rho(rho[k],seos);
	    coords[1] = ses_re_grid_from_engy(e[k],seos);
	    set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_ENERGY);
	    ses_solution(coords,COMP_PURE_PHASE,NULL,POSITIVE_SIDE,fr,wv,var);
	    T = ses_re_temp_from_var(var[RE_T],seos);
	    p[k] = rho[k]*T*var[RE_RP] + var[RE_CP];
	}

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
}		/*end SESAME_single_eos_load_pressure*/

/*
*			SESAME_single_eos_load_sound_speed2():
*
*	Loads a vector of sound speeds into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure.
*	Assumes that the pressure field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	SESAME_single_eos_load_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	Gas_param *params = Params(vst->state[offset]);
	Front *fr;
	SESAME_EOS *seos;
	Wave *wv;
	double coords[MAXD], T, var[NUM_SES_VAR];
	double p;
	double *rho = vst->rho + offset;
	double *e = vst->e + offset;
	double *c2 = vst->c2 + offset;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double *min_pressure = vst->min_pressure + offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int   k;

	seos = (SESAME_EOS *)params->eos;
	fr = seos->fr[SESAME_RHO_ENERGY];
	wv = seos->wave[SESAME_RHO_ENERGY];
	set_ses_intrp_flag(EVALUATE_PRESSURE | EVALUATE_ADB_GAMMA,
	                   SESAME_RHO_ENERGY);
		
	for (k = 0; k < vsize; ++k)
	{
	    coords[1] = ses_re_grid_from_engy(e[k],seos);
	    coords[0] = ses_re_grid_from_rho(rho[k],seos);
	    ses_solution(coords,COMP_PURE_PHASE,NULL,POSITIVE_SIDE,fr,wv,var);
	    T = ses_re_temp_from_var(var[RE_T],seos);
	    p = rho[k]*T*var[RE_RP] + var[RE_CP];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (p < min_pressure[k])
		p = min_pressure[k];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    c2[k] = var[RE_AG]*p/rho[k]; 
	}
	if (vst->FD != NULL)
	{
	    Locstate state = vst->state[offset];
	    double    *FD = vst->FD + offset;
	    static   Locstate tmpst = NULL;

	    if (tmpst == NULL)
	    {
	        (*Params(state)->_alloc_state)(&tmpst,Params(state)->sizest);
		set_type_of_state(tmpst,EGAS_STATE);
	    }
	    for (k = 0; k < vsize; ++k)
	    {
	        Dens(tmpst) = rho[k];
		Energy(tmpst) = e[k];
		Set_params(tmpst,state);
	        FD[k] = fundamental_derivative(tmpst);
	    }
	}
}		/*end SESAME_single_eos_load_sound_speed2*/

/***************END VECTORIZED THERMODYNAMIC FUNCTIONS *********************/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/

/***************Purely Thermodynamic Hugoniot Functions*********************/

/***************End Purely Thermodynamic Hugoniot Functions*****************/
/***************Velocity Related Hugoniot Functions*************************/

/***************End Velocity Related Hugoniot Functions*********************/
/***************Purely Thermodynamic Adiabatic Wave Curve Functions*********/

/*	
*			SESAME_dens_rarefaction():
*
*	Given the state state0 and the pressure on the other side of
*	a simple wave in steady irrotational flow, this
* 	function returns the density on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equation
*
*		dh/dP = V,  h(p0) = h0;
*
*	where h is the specific enthalpy,  and the derivatives are taken
*	at constant entropy.
*/

LOCAL	double	SESAME_dens_rarefaction(
	double p1,
	Locstate state0)
{
	Front		*fr;
	Wave		*wv;
	double		coords[MAXD], var[NUM_SES_VAR];
	SESAME_EOS	*seos = Ses_params(state0);
	COMPONENT	comp;

	fr = seos->fr[SESAME_PRESS_ENTROPY];
	wv = seos->wave[SESAME_PRESS_ENTROPY];
	coords[0] = ses_ps_grid_from_press(p1,seos);
	coords[1] = ses_ps_grid_from_entpy(entropy(state0),seos);
	set_ses_intrp_flag(EVALUATE_DENSITY,SESAME_PRESS_ENTROPY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state0)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	return ses_ps_rho_from_var(var[PS_RHO],seos);
}		/*end SESAME_dens_rarefaction*/

/*	
*			SESAME_pressure_rarefaction():
*
*	Given the state state0 and the density on the other side of
*	a simple wave in steady irrotational flow, this
* 	function returns the pressure on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equation
*
*		de/dV = -P,  e(V0) = e0;
*
*	where e is the specific internal energy,  and the derivatives are taken
*	at constant entropy.
*/

LOCAL	double	SESAME_pressure_rarefaction(
	double rho1,
	Locstate state0)
{
	Front		*fr;
	Wave		*wv;
	double		coords[MAXD];
	double		T, var[NUM_SES_VAR];
	SESAME_EOS	*seos = Ses_params(state0);
	COMPONENT	comp;

	fr = seos->fr[SESAME_RHO_ENTROPY];
	wv = seos->wave[SESAME_RHO_ENTROPY];
	coords[0] = ses_rs_grid_from_rho(rho1,seos);
	coords[1] = ses_rs_grid_from_entpy(entropy(state0),seos);
	set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_ENTROPY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state0)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	T = ses_rs_temp_from_var(var[RS_T],seos);
	return rho1*T*var[RS_RP] + var[RS_CP];
}		/*end SESAME_pressure_rarefaction*/


/*	
*			SESAME_state_on_adiabat_with_pr():
*
*	Given the state state0 and the pressure on the other side of
*	a simple wave in steady irrotational flow, this function returns
*	the thermodynamic variable on the other side.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in state1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*/

LOCAL	void	SESAME_state_on_adiabat_with_pr(
	Locstate state0,
	double p1,
	Locstate state1,
	int stype1)
{
	Front		*fr;
	Wave		*wv;
	double		r1, e1, coords[MAXD];
	double		T, var[NUM_SES_VAR];
	int		interp_flag;
	COMPONENT	comp;
	SESAME_EOS	*seos = Ses_params(state0);
		
	zero_state_velocity(state1,Params(state0)->dim);
	Set_params(state1,state0);
	set_type_of_state(state1,stype1);
		
	if (stype1 == TGAS_STATE)
	{
		Dens(state1) = dens_rarefaction(p1,state0);
		Press(state1) = p1;
		return;
	}
	coords[0] = ses_ps_grid_from_press(p1,seos);
	coords[1] = ses_ps_grid_from_entpy(entropy(state0),seos);
	fr = seos->fr[SESAME_PRESS_ENTROPY];
	wv = seos->wave[SESAME_PRESS_ENTROPY];
	interp_flag = EVALUATE_DENSITY | EVALUATE_ENERGY;
	set_ses_intrp_flag(interp_flag,SESAME_PRESS_ENTROPY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state0)),
			COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	r1 = ses_ps_rho_from_var(var[PS_RHO],seos);
	Dens(state1) = r1;
        T = ses_ps_temp_from_var(var[PS_T],seos);
	e1 = T*var[PS_RE] + var[PS_CE];

	switch(stype1)
	{
	case GAS_STATE:
		Energy(state1) = r1 * e1;
		break;
	case EGAS_STATE:
		Energy(state1) = e1;
		break;
	case VGAS_STATE:
		Energy(state1) = e1;
		set_type_of_state(state1,EGAS_STATE);
		set_state(state1,VGAS_STATE,state1);
		break;
	default:
		screen("ERROR in state_on_adiabat_with_pr()\n"
		       "Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end SESAME_state_on_adiabat_with_pr*/

/*	
*			SESAME_state_on_adiabat_with_dens():
*
*	Given the state state0 and the density on the other side of
*	a simple wave in steady irrotational flow, this	function returns
*	the pressure and internal energy on the other side.
*
*	IMPORTANT NOTES:
*		1.  If stype1 == GAS_STATE the energy in state1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*		2.  Dens(st1) cannot be set to rho1 before the evaluation of
*		the pressure of st0.  This allows this function to work
*		even in the case were st0 = st1 (ie they both point to the
*		same area in storage).
*/

LOCAL	void	SESAME_state_on_adiabat_with_dens(
	Locstate state0,
	double rho1,
	Locstate state1,
	int stype1)
{
	Front		*fr;
	Wave		*wv;
	double		coords[MAXD];
	double		T, var[NUM_SES_VAR];
	COMPONENT	comp;
	SESAME_EOS	*seos = Ses_params(state0);
	SIDE		side;

	Set_params(state1,state0);
	zero_state_velocity(state1,Params(state0)->dim);
	set_type_of_state(state1,stype1);
	if (stype1 == TGAS_STATE)
	{
		Press(state1) = pressure_rarefaction(rho1,state0);
		Dens(state1) = rho1;
		return;
	}
	fr = seos->fr[SESAME_RHO_ENTROPY];
	wv = seos->wave[SESAME_RHO_ENTROPY];
	coords[0] = ses_rs_grid_from_rho(rho1,seos);
	coords[1] = ses_rs_grid_from_entpy(entropy(state0),seos);
	set_ses_intrp_flag(EVALUATE_ENERGY,SESAME_RHO_ENTROPY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state0)),
				COMP_PURE_PHASE,coords,fr->interf);
	side = (multiphase_eos(Eos(state0)) && comp == COMP_PURE_PHASE) ?
	       NEGATIVE_SIDE : POSITIVE_SIDE;
	ses_solution(coords,comp,NULL,side,fr,wv,var);
	T = ses_rs_temp_from_var(var[RS_T],seos);
	Dens(state1) = rho1;
	switch(stype1)
	{
	case GAS_STATE:
		Energy(state1) = rho1*(T*var[RS_RE]+var[RS_CE]);
		break;
	case EGAS_STATE:
		Energy(state1) = T*var[RS_RE] + var[RS_CE];
		break;
	case VGAS_STATE:
		Energy(state1) = T*var[RS_RE] + var[RS_CE];
		set_type_of_state(state1,EGAS_STATE);
		set_state(state1,VGAS_STATE,state1);
		break;
	default:
		screen("ERROR in state_on_adiabat_with_dens(), "
		       "Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end SESAME_state_on_adiabat_with_dens*/




/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/
/***************General Wave Curve Functions********************************/

/*
*			SESAME_mass_flux():
*
*	Returns the mass flux across a wave.
*
*				
*		     | (P - P0) |
*		m  = | -------  |
*		     | (U - U0) |
*
*	Where 
*		P0 = pressure ahead of the shock
*		U0 = velocity ahead of the shock
*		P = pressure behind the shock
*		U = velocity behind the shock
*
*/

LOCAL	double	SESAME_mass_flux(
	double p,
	Locstate state0)
{
	double		p0, rho0, V, V0, denom;
	double		mf;

	p0 = pressure(state0);
	rho0 = Dens(state0);

	mf = ERROR_FLOAT;

	if (fabs(p0 - p) < EPS*p0) /*TOLERANCE*/
		return acoustic_impedance(state0);
#if defined(PHASE_CODE)
	if (multiphase_eos(Eos(state0)))
	{
		WAVE_CURVE	*wave_cur;
		wave_cur = (state_type(state0) == VGAS_STATE) ?
				Wave_curve(state0) : NULL;
		if (wave_cur != NULL && wave_cur->special)
			return mass_flux_for_phase_change(p,state0);
	}
#endif /* defined(PHASE_CODE) */
	if (p < p0)
	{
		denom = int_dp_over_c_rho(p,state0);
		mf = (p0 - p)/denom;
	}
	else
	{
		V = 1.0/dens_Hugoniot(p,state0);
		V0 = 1/rho0;
		mf = sqrt((p - p0) / (V0 - V));
	}
	return mf;
}		/*end SESAME_mass_flux*/

/*
*			SESAME_mass_flux_squared():
*
*	Returns the square of the mass flux across a wave.
*
*				 2
*		 2   | (P - P0) |
*		m  = | -------  |
*		     | (U - U0) |
*
*	Where 
*		P0 = pressure ahead of the shock
*		U0 = velocity ahead of the shock
*		P = pressure behind the shock
*		U = velocity behind the shock
*
*/

LOCAL	double	SESAME_mass_flux_squared(
	double p,
	Locstate state0)
{
	double		p0, rho0, V, V0, denom;
	double		m;

	p0 = pressure(state0);
	rho0 = Dens(state0);

	if (fabs(p0 - p) < EPS*p0) /*TOLERANCE*/
		return acoustic_impedance_squared(state0);

#if defined(PHASE_CODE)
	if (multiphase_eos(Eos(state0)))
	{
		WAVE_CURVE	*wave_cur;
		wave_cur = (state_type(state0) == VGAS_STATE) ?
				Wave_curve(state0) : NULL;
		if (wave_cur != NULL && wave_cur->special)
		{
			m = mass_flux_for_phase_change(p,state0);
			return m*m;
		}
	}
#endif /* defined(PHASE_CODE) */
	if (p < p0)
	{
		denom = int_dp_over_c_rho(p,state0);
		m = (p0 - p)/denom;
		return m*m;
	}
	else
	{
		V = 1.0/dens_Hugoniot(p,state0);
		V0 = 1/rho0;
		return (p - p0) / (V0 - V);
	}
}		/*end SESAME_mass_flux_squared*/


/***************End General Wave Curve Functions****************************/
/***************Functions for the Evaluation of Riemann Solutions***********/

/*
*				SESAME_oned_fan_state():
*
*	This is a utility function provided for the evaluation of states
*	in a simple wave.   Given sta, it solves for stm using the
*	equation:
*
*	                     / p_m        |            / c_m        |
*	                    /             |           /             |
*	                    \       dP    |           \        dc   |
*	    w = c_m - c_a +  \    -----   |         =  \     ------ |
*	                      \   rho c   |             \     mu^2  |
*	                      /           |             /           |
*	                     /p_a         | S = S_a    / c_a        | S = S_a
*
*	here c is the sound speed,  rho the density,  S the specific entropy,
*	p the pressure,  and mu^2 = (G - 1)/G,  where G is the fundamental
*	derivative of gas dynamics.  The returned state1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	state1 can be used to evaluate the pressure, density, and sound speed
*	of the state inside the rarefaction fan.
*	
*	Input data:
*		w = value of w as defined above
*		sta = state ahead of fan
*		stb = state behind fan
*
*	Output data:
*		stm = state inside fan
*		vacuum = 1 if stm is a vacuum,  0 otherwise
*
*	Returns the sound speed of the answer state stm.
*/

typedef struct {
	Locstate sta, stm;
	int stype_m;
	double c_m;
	double S0_grid, icor0;
	Front *fr;
	Wave *wv;
	boolean phase;
} RFUNC_PARAMS;


LOCAL	double	SESAME_oned_fan_state(
	double    w,
	Locstate sta,
	Locstate stb,
	Locstate stm,
	int      stype_m,
	boolean  *vacuum)
{
	COMPONENT	comp_r;
	Front		*fr;
	Locstate	str; /* sta or stb if sta is vacuum */
	RFUNC_PARAMS	Rfpr;
	Wave		*wv;
	double		cref;
	double		psi;
	double		c_m, c_r, S_r, crds_r[MAXD];
	double		epsilon, delta;
	double		var[NUM_SES_VAR];
	double		psi_b, psi_a;
	int		interp_flag;
	SIDE		side;
#if defined(PHASE_CODE)
	COMPONENT	comp_b;
	double		crds_b[MAXD];
#endif /* defined(PHASE_CODE) */
#if !defined(RAREFACTION_STATE_BY_PRESSURE)
	double		rho_r;
	double		rho_b_grid, rho_a_grid, rho_m_grid;
	static const double PSI_FAC = 0.1; /* TOLERANCE */

	Set_params(stm,sta);
	zero_state_velocity(stm,Params(sta)->dim);
	*vacuum = NO;

	str = sta;
	c_r = sound_speed(str);
	psi = w + c_r;

	cref = Reference_sound_speed(Ses_params(str));
	rho_r = Dens(str);
	S_r = entropy(str);
	crds_r[0] = ses_rs_grid_from_rho(rho_r,Ses_params(str));
	epsilon = Eps(str)*fabs(psi);
	delta = Eps(str)*fabs(crds_r[0]);
	if (delta <= 0.0) delta = Eps(str);

	Rfpr.sta = str;
	Rfpr.stm = stm;
	set_type_of_state(stm,TGAS_STATE);
	Rfpr.stype_m = state_type(stm);
	Rfpr.phase = multiphase_eos(Eos(str));
	Rfpr.S0_grid = crds_r[1] = ses_rs_grid_from_entpy(S_r,Ses_params(str));
	Rfpr.fr= fr = Ses_params(str)->fr[SESAME_RHO_ENTROPY];
	Rfpr.wv = wv = Ses_params(str)->wave[SESAME_RHO_ENTROPY];
	comp_r = nearest_interior_comp(Rfpr.phase,
				COMP_PURE_PHASE,crds_r,fr->interf);
#if defined(PHASE_CODE)
	crds_b[0] = ses_rs_grid_from_rho(Dens(stb),Ses_params(stb));
	crds_b[1] = ses_rs_grid_from_entpy(entropy(stb),Ses_params(stb));
	comp_b = nearest_interior_comp(Rfpr.phase,
			COMP_PURE_PHASE,crds_b,fr->interf);
	if (Rfpr.phase)
	{
		if (comp_r==COMP_PURE_PHASE && comp_b==COMP_PURE_PHASE)
			side = NEGATIVE_SIDE;
		else if (comp_r==COMP_MIXED_PHASE && comp_b==COMP_MIXED_PHASE)
			side = POSITIVE_SIDE;
		else if (comp_r==COMP_PURE_PHASE && comp_b==COMP_MIXED_PHASE)
			side = NEGATIVE_SIDE;
		else if (comp_r==COMP_MIXED_PHASE && comp_b==COMP_PURE_PHASE)
			side = POSITIVE_SIDE;
		else
		{
			screen("ERROR in ses_oned_fan_state(), "
			       "THREE PHASE CODE NEEDED\n");
			clean_up(ERROR);
		}
	}
	else
#else /* defined(PHASE_CODE) */
		side = POSITIVE_SIDE;
#endif /* defined(PHASE_CODE) */

	interp_flag = EVALUATE_ADB_GAMMA | EVALUATE_IDPOCR;
	set_ses_intrp_flag(interp_flag,SESAME_RHO_ENTROPY);
	ses_solution(crds_r,comp_r,NULL,side,fr,wv,var);
	Rfpr.icor0 = 0.5*(c_r+cref)*var[RS_RF];
	set_ses_intrp_flag(EVALUATE_PRESSURE | EVALUATE_ADB_GAMMA |
			EVALUATE_IDPOCR, SESAME_RHO_ENTROPY);
	if (find_root(ses_r_cprf,(POINTER)&Rfpr,
			psi,&rho_m_grid,
			ses_rs_grid_from_rho(Dens(stb),Ses_params(stb)),
			ses_rs_grid_from_rho(Dens(sta),Ses_params(sta)),
			epsilon,delta) == FUNCTION_FAILED)
	{
		rho_b_grid = ses_rs_grid_from_rho(Dens(stb),Ses_params(stb));
		rho_a_grid = ses_rs_grid_from_rho(Dens(sta),Ses_params(sta));
		(void) ses_r_cprf(rho_b_grid,&psi_b,(POINTER)&Rfpr);
		(void) ses_r_cprf(rho_a_grid,&psi_a,(POINTER)&Rfpr);
		if(fabs(psi_a - psi) <= PSI_FAC*fabs(psi))
			rho_m_grid = rho_a_grid;
		else if(fabs(psi_b - psi) <= PSI_FAC*fabs(psi))
			rho_m_grid = rho_b_grid;
		else
		{
		    screen("ERROR in ses_oned_fan_state(),"
		           "SESAME bisection failed.\n");
		    verbose_print_state("sta",sta);
		    verbose_print_state("stb",stb);
		    (void) printf("psi = %g, epsilon = %g, delta = %g\n",
		    	psi,epsilon,delta);
		    (void) printf("rho_b = %g, psi_b = %g\n",
		    	Dens(stb),psi_b);
		    (void) printf("rho_a = %g, psi_a = %g\n",Dens(sta),psi_a);
		    clean_up (ERROR);
		}
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Dens(Rfpr.stm) = max(Dens(Rfpr.stm),Vacuum_dens(sta));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	set_state(stm,stype_m,Rfpr.stm);
	c_m = Rfpr.c_m;

#else /* !defined(RAREFACTION_STATE_BY_PRESSURE) */
	double p_r;
	double p_b_grid, p_a_grid;
	double p_m_grid, p_b, p_a;

	Set_params(stm,sta);
	zero_state_velocity(stm,dim);

	p_b = pressure(stb);
	p_a = pressure(sta);
	str = sta;
	c_r = sound_speed(str);
	p_r = p_a;
	psi = w + c_r;

	Rfpr.sta = str;
	Rfpr.stm = stm;
	Rfpr.stype_m = TGAS_STATE;
	S_r = entropy(str);
	crds_r[0] = ses_ps_grid_from_press(p_r,Ses_params(str));
	epsilon = Eps(str)*fabs(psi);
	delta = Eps(str)*fabs(crds_r[0]);
	Rfpr.phase = multiphase_eos(Eos(str));
	Rfpr.S0_grid = crds_r[1] = ses_ps_grid_from_entpy(S_r,Ses_params(str));
	Rfpr.fr= fr = Ses_params(str)->fr[SESAME_PRESS_ENTROPY];
	Rfpr.wv = wv = Ses_params(str)->wave[SESAME_PRESS_ENTROPY];
	comp_r = nearest_interior_comp(Rfpr.phase,
				COMP_PURE_PHASE,crds_r,fr->interf);
#if defined(PHASE_CODE)
	crds_b[0] = ses_ps_grid_from_press(pressure(stb),Ses_params(stb));
	crds_b[1] = ses_ps_grid_from_entpy(entropy(stb),Ses_params(stb));
	comp_b = nearest_interior_comp(Rfpr.phase,
				COMP_PURE_PHASE,crds_b,fr->interf);
	if (Rfpr.phase)
	{
		if (comp_r==COMP_PURE_PHASE && comp_b==COMP_PURE_PHASE)
			side = NEGATIVE_SIDE;
		else if (comp_r==COMP_MIXED_PHASE && comp_b==COMP_MIXED_PHASE)
			side = POSITIVE_SIDE;
		else if (comp_r==COMP_PURE_PHASE && comp_b==COMP_MIXED_PHASE)
			side = POSITIVE_SIDE;
		else if (comp_r==COMP_MIXED_PHASE && comp_b==COMP_PURE_PHASE)
			side = NEGATIVE_SIDE;
		else
		{
			screen("ERROR in ses_oned_fan_state(), "
			       "THREE PHASE CODE NEEDED\n");
			clean_up(ERROR);
		}
	}
	else
#else /* defined(PHASE_CODE) */
		side = POSITIVE_SIDE;
#endif /* defined(PHASE_CODE) */

	interp_flag = EVALUATE_ADB_GAMMA | EVALUATE_IDPOCR;
	set_ses_intrp_flag(interp_flag,SESAME_PRESS_ENTROPY);
	ses_solution(crds_r,comp_r,NULL,side,fr,wv,var);
	Rfpr.icor0 = c_r*var[PS_RF];
	set_ses_intrp_flag(EVALUATE_DENSITY | EVALUATE_ADB_GAMMA |
			EVALUATE_IDPOCR, SESAME_PRESS_ENTROPY);
	if (find_root(ses_p_cprf,(POINTER)&Rfpr,
		psi,&p_m_grid,ses_ps_grid_from_press(p_b,Ses_params(stb)),
		ses_ps_grid_from_press(p_a,Ses_params(sta)),
		epsilon,delta) == FUNCTION_FAILED)
	{
		p_b_grid = ses_ps_grid_from_press(p_b,Ses_params(stb));
		p_a_grid = ses_ps_grid_from_press(p_a,Ses_params(sta));
		(void) ses_p_cprf(p_b,&psi_b,(POINTER)&Rfpr);
		(void) ses_p_cprf(p_a_grid,&psi_a,(POINTER)&Rfpr);
		if (fabs(psi_a - psi) <= PSI_FAC*fabs(psi))
			p_m_grid = p_a_grid;
		else if (fabs(psi_b - psi) <= PSI_FAC*fabs(psi))
			p_m_grid = p_b_grid;
		else
		{
		    screen("ERROR in ses_oned_fan_state(),"
		           "no root found.\n");
		    verbose_print_state("sta",sta);
		    verbose_print_state("stb",stb);
		    (void) printf("psi = %g, epsilon = %g, delta = %g\n",
		    	psi,epsilon,delta);
		    (void) printf("rmid = %g, psi_b = %g\n",
		    	Dens(stb),psi_b);
		    (void) printf("rho0 = %g, psi_a = %g\n",Dens(sta),psi_a);
		    	clean_up (ERROR);
		}
	}
	set_state(stm,stype_m,Rfpr.stm);
	c_m = Rfpr.c_m;
#endif /* !defined(RAREFACTION_STATE_BY_PRESSURE) */

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Press(stm) < Min_pressure(sta))
	{
		state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
		c_m = 0.0;
		*vacuum = YES;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return c_m;
}		/* end oned_fan_state*/

#if !defined(RAREFACTION_STATE_BY_PRESSURE)
/*
*			ses_r_cprf():
*
*	Computes c + integral from rho0 to rho of  c/rho drho
*	for the sesame hyp solution equation of state.
*
*	IMPORT:  Evaluation flag must be set before a call to this
*	function.
*
*	set_ses_intrp_flag(EVALUATE_PRESSURE | EVALUATE_ADB_GAMMA |
*					EVALUATE_IDPOCR,SESAME_RHO_ENTROPY);
*/

LOCAL	boolean	ses_r_cprf(
	double		rho_grid,
	double		*frho,
	POINTER		parameters)
{
	RFUNC_PARAMS	*prms = (RFUNC_PARAMS *) parameters;
	SESAME_EOS	*seos;
	double		coords[MAXD], rho, c;
	double		pr, T, var[NUM_SES_VAR];
	double		cref;
	Front		*fr;
	Wave		*wv;
	boolean		phase;
	SIDE		side;
	COMPONENT	comp;
	
	coords[0] = rho_grid;
	coords[1] = prms->S0_grid;
	fr = prms->fr;
	wv = prms->wv;
	seos = Ses_front_seos(fr);
	*frho = prms->icor0;
	rho = ses_rs_rho_from_grid(rho_grid,seos);
	cref = Reference_sound_speed(seos);

	phase = prms->phase;
	comp = nearest_interior_comp(phase,COMP_PURE_PHASE,coords,fr->interf);
#if defined(PHASE_CODE)
	side = (phase && comp == COMP_PURE_PHASE) ?
		NEGATIVE_SIDE : POSITIVE_SIDE;
#else /* defined(PHASE_CODE) */
	side = POSITIVE_SIDE;
#endif /* defined(PHASE_CODE) */
	ses_solution(coords,comp,NULL,side,fr,wv,var);
	T = ses_rs_temp_from_var(var[RS_T],seos);
	pr = rho*T*var[RS_RP] + var[RS_CP];
	prms->c_m = c = sqrt(pr*var[RS_AG]/rho);
	*frho = c - *frho + 0.5*(c+cref)*var[RS_RF];
	Press(prms->stm) = pr;	/*Assumes stm is TGAS_STATE*/
	Dens(prms->stm) = rho;
	return FUNCTION_SUCCEEDED;
}		/*end ses_r_cprf*/

#else /* !defined(RAREFACTION_STATE_BY_PRESSURE) */



/*
*			ses_p_cprf():
*
*	Computes c + integral from p0 to p of  c/rho drho
*	for the sesame hyp solution equation of state.
*
*	IMPORT:  Evaluation flag must be set before a call to this
*	function.
*
*	set_ses_intrp_flag(EVALUATE_DENSITY | EVALUATE_ADB_GAMMA |
*					EVALUATE_IDPOCR,SESAME_PRESS_ENTROPY);
*/

LOCAL	boolean	ses_p_cprf(
	double		pr_grid,
	double		*frho,
	POINTER		parameters)
{
	RFUNC_PARAMS	*prms = (RFUNC_PARAMS *) parameters;
	double		coords[MAXD], rho, c;
	double		pr, var[NUM_SES_VAR];
	double		cref;
	Front		*fr;
	Wave		*wv;
	SESAME_EOS	*seos;
	COMPONENT	comp;
	boolean		phase;
	
	coords[0] = pr_grid;
	coords[1] = prms->S0_grid;
	fr = prms->fr;
	wv = prms->wv;
	seos = Ses_front_seos(fr);
	cref = Reference_sound_speed(seos);
	*frho = prms->icor0;
	pr = ses_ps_press_from_grid(pr_grid,seos);
	phase = prms->phase;
	comp = nearest_interior_comp(phase,COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	rho = ses_ps_rho_from_var(var[PS_RHO],seos);
	c = prms->c_m = sqrt(var[PS_AG]*pr/rho);
	*frho = c - *frho + 0.5*(c+cref)*var[PS_RF];
	Dens(prms->stm) = rho;	/* Assumes stm is a TGAS_STATE*/
	Press(prms->stm) = pr;
	return FUNCTION_SUCCEEDED;
}		/*end ses_p_cprf*/
#endif /* !defined(RAREFACTION_STATE_BY_PRESSURE) */
/***************End Functions for the Evaluation of Riemann Solutions********/



/***************Functions to Compute Riemann Solutions**********************/


/*
*			SESAME_riemann_wave_curve():
*
*	Evalutes the forward wave family wave curve defined by
*
*		 _
*		|
*		|
*		|                                1/2
*               |   [ (Pstar  -  P0) * ( V0 - V) ]     if Pstar > P0
*		|
*		|
*	        / 
*	       /
*              \
*		\		
*		|
*               |        / Pstar     |
*               |       /            |
*               |       \      dP    |
*               |        \   ------  |		       if Pstar < P0
*               |         \   rho c  |
*               |         /          |
*               |        / P0        | S
*               |_
*
*/

LOCAL	double	SESAME_riemann_wave_curve(
	Locstate state0,
	double pstar)
{
	double		rho0, rho_star, p0;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (pstar < Min_pressure(state0)) pstar = Min_pressure(state0);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	
	p0 = pressure(state0);
	if (pstar >= p0)
	{
		rho0 = Dens(state0);
		rho_star = dens_Hugoniot(pstar,state0);
		return (rho_star <= rho0) ? 0.0 :
			sqrt((pstar-p0)*(rho_star-rho0)/(rho0*rho_star));
	}
	else
		return -int_dp_over_c_rho(pstar,state0);
}		/*end SESAME_riemann_wave_curve*/


/*
*			SESAME_eps_for_Godunov():
*
*	Returns a tolerance to be used to determine convergence of the
*	of the Riemann solver.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	double	SESAME_eps_for_Godunov(
	Locstate state,
	double pstar,
	double r_eps)
{
	double pr_r, delta, sesame_eps;

	if (state == NULL || is_obstacle_state(state)) return 0.0;

	sesame_eps = Eps(state);
	pr_r = pressure(state);
	delta = min(1.0,fabs((pstar - pr_r)/pr_r));
	return (delta<=0.001) ? sesame_eps : sesame_eps/delta;
}		/*end SESAME_eps_for_Godunov*/

/*
*			SESAME_initialize_riemann_solver()
*
*	Computes the epsilons and the initial guess for the pressure
*	in the secant iteration of find_mid_state.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	void	SESAME_initialize_riemann_solver(
	Locstate Tsl,
	Locstate Tsr,
	double *pstar,
	double *p_min,
	double eps,
	double *eps_u,
	double *eps_p,
	boolean (*fd_md_st)(Locstate,Locstate,double,double*,double*,
			 double*,double*,double*,double*,
			 RIEMANN_SOLVER_WAVE_TYPE*,RIEMANN_SOLVER_WAVE_TYPE*))
{
	double	     pl, pr, ptmp;
	static const double EPS_FAC = 10.0;	/* TOLERANCE */

	pl = pressure(Tsl), pr = pressure(Tsr);
	if (Eos(Tsl) != Eos(Tsr))
	{
#if defined(UNRESTRICTED_THERMODYNAMICS)
	    *p_min = -HUGE_VAL;
#else /* defined(UNRESTRICTED_THERMODYNAMICS) */
	    *p_min = max(Min_pressure(Tsl),Min_pressure(Tsr));
#endif /* defined(UNRESTRICTED_THERMODYNAMICS) */
	    *eps_u = *eps_p = eps;
	    ptmp = ses_ps_press_from_grid(
		Ses_params(Tsl)->fr[SESAME_PRESS_ENTROPY]->rect_grid->L[0],
		Ses_params(Tsl));
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    *p_min = max(*p_min,ptmp);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    *eps_p = max(Eps(Tsl),*eps_p);
	    ptmp = ses_ps_press_from_grid(
		Ses_params(Tsr)->fr[SESAME_PRESS_ENTROPY]->rect_grid->L[0],
		Ses_params(Tsr));
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    *p_min = max(*p_min,ptmp);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    *eps_p = max(Eps(Tsr),*eps_p);
	    *eps_u = EPS_FAC *(*eps_p);
	    *pstar = 0.5*(pl + pr);
	    *pstar = max(*pstar,*p_min);
	    return;
	}

	*eps_p = Eps(Tsl);
	*eps_u = EPS_FAC * Eps(Tsl);

	*p_min = max(Pmin(Tsl),Pmin(Tsr));
	*pstar = 0.5*(max(pl,*p_min) + max(pr,*p_min));
}		/*end SESAME_initialize_riemann_solver*/


/***************End Functions to Compute Riemann Solutions******************/
/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/

/***************TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*********/

/***************END TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*****/


/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			SESAME_prompt_for_state():
*
*	Prompts for a hydrodynamical state.  The form of
*	the states depends of the Eos. 	The type of the state
*	is returned.
*/

LOCAL	void	SESAME_prompt_for_state(
	Locstate   state,
	int        stype,
	Gas_param  *params,
	const char *mesg)
{
	int i, dim;
	static  char velmesg[3][11] = {"x velocity","y velocity","z velocity"};

	if (params == NULL)
	{
		g_obstacle_state(state,g_sizest());
		return;
	}
	dim = params->dim;
	set_type_of_state(state,TGAS_STATE);
	Params(state) = params;
	screen("Enter the density, pressure");
	for (i = 0; i < dim; ++i)
	{
		screen(", ");
		if (i == (dim - 1)) screen("and ");
		screen("%s",velmesg[i]);
	}
	screen("%s: ",mesg);
	(void) Scanf("%f %f",&Dens(state),&Press(state));
	for (i = 0; i < dim; ++i)
		(void) Scanf("%f",&Vel(state)[i]);
	(void) getc(stdin); /*read trailing newline*/

	set_state(state,stype,state);
}		/* end SESAME_prompt_for_state */

/*
*			SESAME_prompt_for_thermodynamics():
*
*	Prompts for the thermodynamic variables in a state.  Returns
*	a state with the appropriate thermodynamic state and zero velocity.
*	The return status gives the state type representation of the state.
*/

LOCAL	void	SESAME_prompt_for_thermodynamics(
	Locstate   state,
	Gas_param  *params,
	const char *mesg)
{
	if (params == NULL)
	{
		g_obstacle_state(state,g_sizest());
		return;
	}
	set_type_of_state(state,TGAS_STATE);
	zero_state_velocity(state,MAXD);
	Params(state) = params;
	screen("Enter the density and pressure %s: ",mesg);
	(void) Scanf("%f %f\n",&Dens(state),&Press(state));
}		/* end SESAME_prompt_for_thermodynamics */

/*
*			SESAME_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

LOCAL	void	SESAME_fprint_EOS_params(
	FILE *file,
	Gas_param *params)
{
	SESAME_EOS	*seos = (SESAME_EOS*)params->eos;

	fprint_SESAME_params(file,seos);

	(void) fprintf(file,"Sesame restart library = %s\n",
		            (strlen(seos->ses_restart_file_name) != 0) ?
		            seos->ses_restart_file_name : "NONE");


}		/*end SESAME_fprint_EOS_params */

/*
*			SESAME_read_print_EOS_params():
*
*	Reads the equation of state data as printed by SESAME_fprint_EOS_params.
*	This is restart function.
*/

/*ARGSUSED*/
LOCAL	void	SESAME_read_print_EOS_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	FILE       *file = io_type->file;
	SESAME_EOS *seos = (SESAME_EOS*)params->eos;

	read_print_SESAME_params(seos,io_type);

	(void) fgetstring(file,"Sesame restart library = ");
	(void) fscanf(file,"%s",seos->ses_restart_file_name);

	if (strcmp(seos->ses_restart_file_name,"NONE") == 0)
	{
	    (void) strcpy(seos->ses_restart_file_name,"");
	    if (initialize_Sesame_tables(seos) == FUNCTION_FAILED)
	    {
	        screen("ERROR in SESAME_read_print_EOS_params(), "
	               "can't initialize tables\n");
	        clean_up(ERROR);
	    }
	    sesinv(seos);
	}
	else
	{
	    IO_TYPE IO_type;
	    FILE *file;
	    if ((file = fopen(seos->ses_restart_file_name,"r")) == NULL)
	    {
	    	screen("ERROR in SESAME_read_print_EOS_params(), "
	    	       "Unable to open %s\n",seos->ses_restart_file_name);
	    	clean_up(ERROR);
	    }
	    determine_io_type(file,&IO_type);
	    seos->restart_io_type = &IO_type;
	    if (restart_sesame(init,seos) == FUNCTION_FAILED)
	    {
	    	screen("ERROR in SESAME_read_print_EOS_params(), "
	    	       "restart_sesame() failed\n");
	    	clean_up(ERROR);
	    }
	    (void) fclose(IO_type.file);
	    seos->restart_io_type = NULL;
	}
}		/*end SESAME_read_print_EOS_params*/

/*
*			SESAME_free_EOS_params():
*
*	Frees the storage allocated for an equation of state parameter
*	function.
*/

LOCAL	EOS*	SESAME_free_EOS_params(
	EOS	*eos)
{
	SESAME_EOS	*seos = (SESAME_EOS*)eos;
	int i;

	for (i = 0; i < seos->num_ses_tabs; ++i)
	{
		free_wave_pointers(seos->wave[i]);
		free_front(seos->fr[i]);
		free(seos->fr[i]->rect_grid);
		free(seos->wave[i]);
		free(seos->fr[i]);
		free(seos->root[i]);
	}
	free_these(3,seos->root,seos->fr,seos->wave);
	free(seos);
	return NULL;
}		/*end SESAME_free_EOS_params*/

/*
*			SESAME_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

LOCAL	void	SESAME_prompt_for_EOS_params(
	INIT_DATA     *init,
	Gas_param     *params,
	const char    *message1,
	const char    *message2)
{
	SESAME_EOS	*seos = (SESAME_EOS*)params->eos;
	FILE            *file;
	static const char 	*DEFAULT_SES_LIBRARY =
			"../../databases/gas/sesdata/m101394.bin";
	int i;
	char s[Gets_BUF_SIZE];

	if (seos == NULL)
	{
	    screen("ERROR in prompt_for_sesame(), "
	           "Unable to allocate more storage\n");
	    clean_up(ERROR);
	}
 	screen("If available, enter the name of a file with "
	       "the SESAME interfaces?\n\t(dflt=no file available): ");
	(void) Gets(seos->ses_restart_file_name);
 	if (seos->ses_restart_file_name[0] != '\0') 
 	{
	    IO_TYPE IO_type;
	    FILE    *file;
	    if ((file = fopen(seos->ses_restart_file_name,"r")) == NULL)
 	    {
 	    	screen("ERROR in SESAME_prompt_for_EOS_params(), "
 	    	       "Unable to open %s\n",s);
 	    	clean_up(ERROR);
 	    }
	    start_clock("read_interfaces");
	    determine_io_type(file,&IO_type);
	    seos->restart_io_type = &IO_type;
 	    if (restart_sesame(init,seos) == FUNCTION_FAILED)
 	    {
 	    	screen("ERROR in SESAME_prompt_for_EOS_params(), "
	    	       "restart_sesame() failed\n");
 	    	clean_up(ERROR);
 	    }
 	    (void) fclose(IO_type.file);
	    seos->restart_io_type = NULL;
	    stop_clock("read_interfaces");
 	}
	else
	{
	    (void) strcpy(seos->seslib,DEFAULT_SES_LIBRARY);
	    screen("Enter the file name of the sesame library (default=%s): ",
		    seos->seslib);
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) strcpy(seos->seslib,s);
	    if ((file = fopen(seos->seslib,"r")) == NULL)
	    {
	        screen("ERROR in SESAME_prompt_for_EOS_params(), "
		       "no such file %s\n",seos->seslib);
	        clean_up(ERROR);
	    }
	    (void) fclose(file);
	    file = NULL;
	    for (i = 0; i < 2; ++i)
	    {
	        screen("Enter the material index for the%s gas%s.\n"
	               "Default or P prints list.\n"
	               "\tEnter index here: ",
	        	message1,message2);
	        (void) Gets(s);
	        if (s[0] == '\0' || s[0] == 'P' || s[0] == 'p')
	        {
	            IMPORT	boolean	suppress_prompts;
	            print_list_of_sesame_materials(stdout,seos->seslib);
	            if (suppress_prompts == NO)
	                print_list_of_sesame_materials(stderr,seos->seslib);
	        }
	        else
	        {
	            (void) sscanf(s,"%d",&seos->ids2);
	            break;
	        }
	    }
	    if (i >= 2)
	    {
	        screen("ERROR in prompt_for_sesame(), "
	               "What's the matter aren't 10 chances enough for you\n"
	               "to find the material you want?\n");
	        clean_up(ERROR);
	    }

	    multiphase_eos(seos) = NO;
	    seos->_sget = FORTRAN_NAME(s2get);
#if defined(PHASE_CODE)
	    screen("Do you wish phase boundaries to be included? [%s]: ",
	    	( multiphase_eos(seos)==YES)?"y(dflt),n":"y,n(dflt)");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	        multiphase_eos(seos) = YES;
	        seos->_sget = FORTRAN_NAME(s4get);
	    }
	    if (s[0] == 'n' || s[0] == 'N')
	    {
	        multiphase_eos(seos) = NO;
	        seos->_sget = FORTRAN_NAME(s2get);
	    }
#endif /* defined(PHASE_CODE) */

	    if (initialize_Sesame_tables(seos) == FUNCTION_FAILED)
	    {
	        screen("ERROR in SESAME_prompt_for_EOS_params(), "
	               "can't initialize tables\n");
	        clean_up(ERROR);
	    }
	    prompt_for_window_params(seos);
	    sesinv(seos);
	}
}		/*end SESAME_prompt_for_EOS_params*/

LOCAL	void	prompt_for_window_params(
	SESAME_EOS *seos)
{
	char	           s[Gets_BUF_SIZE];
	RECT_GRID          Gr;
	REMAP              Remap;
	double	           *tbls = seos->sestab.tbls;
	const double	   eps_mach = 100.0*MACH_EPS;/*TOLERANCE*/
	double	           *rho, *T;
	double	           dlogr, dlogT;
	double	           p[3], e[3];
	double	           csq;
	int	           nr = irint(tbls[2]), nt = irint(tbls[3]);
	int	           ir_min, ir_max, jt_min, jt_max;
	int	           i, j;
	static const double MIN_REL_ERROR = 1.0e-07; /*TOLERANCE*/

	rho = tbls+4;
	T = tbls+4+nr;

	set_remap(2,IDENTITY_REMAP,&Remap);
	dlogr = 0.0;
	for (i = 1; i < nr; ++i)
	    dlogr += (rho[i] - rho[i-1])/rho[i];
	dlogr /= tbls[2]-1.0;
	dlogT = 0.0;
	for (j = 1; j < nt; ++j)
	    dlogT += (T[j] - T[j-1])/T[j];
	dlogT /= tbls[3] - 1.0;

	seos->eps = 0.5*(dlogr + dlogT);
	screen("Enter a small positive constant to measure the\n"
	       "\tdiscreteness of the sesame table (default=%g): ",seos->eps);
	(void) Gets(s);
	(void) sscan_float(s,&seos->eps);

	screen("Choose a density and temperature window for "
	       "the local inversions.\n");
	screen("The %d density points in the table\n",nr);
	for (i = 0; i < nr; ++i)
	    screen("%s%12g%s", (i%5 == 0) ? "    " : "",rho[i],
		            ((i%5 == 4) || (i == (nr-1))) ? "\n"   : " ");
	screen("The %d temperature points in the table\n",nt);
	for (j = 0; j < nt; ++j)
	    screen("%s%12g%s",(j%5 == 0) ? "    " : "",T[j],
			    ((j%5 == 4) || (j == (nt - 1))) ? "\n" : " ");

	for (ir_min = 0; ((rho[ir_min] <= 0.0) && (ir_min < nr)); ++ir_min);
	ir_max = nr - 1;
	for (jt_min = 0; ((T[jt_min] <= 0.0) && (jt_min < nr)); ++jt_min);
	jt_max = nt - 1;

	while ((ir_min < ir_max) && (jt_min < jt_max))
	{
	    double r, t, cv, cp;
	    for (i = ir_min; i <= ir_max; ++i)
	    {
		r = rho[i];

		t = T[jt_min];
	        s2eos_lookup(r,t,p,e,seos);
		cv = e[2];
	        csq = sesame_rt_sound_speed_squared(r,t,p,e);
	        if ((cv <= 0.0) || (csq <= 0.0) || (p[1] == 0))
		{
	            ++jt_min;
		    break;
		}
		cp = cv + (p[0]/(r*r) - e[1])*p[2]/p[1];
		if (cp < cv)
		{
	            ++jt_min;
		    break;
		}

		t = T[jt_max];
	        s2eos_lookup(r,t,p,e,seos);
		cv = e[2];
	        csq = sesame_rt_sound_speed_squared(r,t,p,e);
	        if ((cv <= 0.0) || (csq <= 0.0) || (p[1] == 0))
		{
	            jt_max--;
		    break;
		}
		cp = cv + (p[0]/(r*r) - e[1])*p[2]/p[1];
		if (cp < cv)
		{
	            jt_max--;
		    break;
		}
	    }

	    for (j = jt_min; j <= jt_max; ++j)
	    {
		t = T[j];

		r = rho[ir_min];
	        s2eos_lookup(r,t,p,e,seos);
		cv = e[2];
	        csq = sesame_rt_sound_speed_squared(r,t,p,e);
		if ((cv <= 0.0) || (csq <= 0.0) || (p[1] == 0))
	        {
	            ++ir_min;
		    break;
	        }
		cp = cv + (p[0]/(r*r) - e[1])*p[2]/p[1];
		if (cp < cv)
	        {
	            ++ir_min;
		    break;
	        }

		r = rho[ir_max];
	        s2eos_lookup(r,t,p,e,seos);
		cv = e[2];
	        csq = sesame_rt_sound_speed_squared(r,t,p,e);
		if ((cv <= 0.0) || (csq <= 0.0) || (p[1] == 0))
	        {
	            --ir_max;
		    break;
	        }
		cp = cv + (p[0]/(r*r) - e[1])*p[2]/p[1];
		if (cp < cv)
	        {
	            --ir_max;
		    break;
	        }
	    }
	    if ((i > ir_max) && (j > jt_max))
		break;
	}
	if ((ir_min >= ir_max) || (jt_min >= jt_max))
	{
	    screen("ERROR in prompt_for_window_params(), corrupt "
		   "sesame table\n");
	    clean_up(ERROR);
	}

	/*DEFAULT values for density limits*/
	Rho_min(seos) = rho[ir_min];
	Rho_max(seos) = rho[ir_max];
	Rho_ref(seos) = tbls[1];

	/*DEFAULT values for temperature limits*/
	Temp_min(seos) = T[jt_min];
	Temp_max(seos) = T[jt_max];
	Temp_ref(seos) = 300.0;

	screen("Enter the minimum density "
	       "(default = %g gram/cc): ",Rho_min(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Rho_min(seos));
	for (i = nr-1; i >= 0; i--)
	    if (rho[i] <= Rho_min(seos))
		break;
	if (i < ir_min)
	    (void) printf("WARNING in prompt_for_window_params(), "
			  "unstable thermodynamics in window range\n");
	ir_min = max(i,0);

 	screen("Enter the maximum density (default = %g gram/cc): ",
	       Rho_max(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Rho_max(seos));
	for (i = 0; i < nr; ++i)
	    if (rho[i] >= Rho_max(seos))
		break;
	if (i > ir_max)
	    (void) printf("WARNING in prompt_for_window_params(), "
			  "unstable thermodynamics in window range\n");
	ir_max = min(i,nr-1);
	if ((Rho_max(seos) < Rho_min(seos)) || (ir_max < ir_min))
	{
	    screen("ERROR in prompt_for_window_params(), "
		   "invalid density range\n");
	    clean_up(ERROR);
	}

 	screen("Enter the reference density (default = %g gram/cc): ",
	       Rho_ref(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Rho_ref(seos));

	screen("Enter the minimum temperature(Kelvin) "
	       "(default = %g degrees Kelvin): ",Temp_min(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Temp_min(seos));
	for (j = nt-1; j >= 0; j--)
	    if (T[j] <= Temp_min(seos))
		break;
	if (j < jt_min)
	    (void) printf("WARNING in prompt_for_window_params(), "
			  "unstable thermodynamics in window range\n");
	jt_min = max(j,0);

	screen("Enter the maximum temperature(Kelvin) "
	       "(default = %g degrees Kelvin): ",Temp_max(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Temp_max(seos));
	for (j = 0; j < nt; ++j)
	    if (T[j] >= Temp_max(seos))
		break;
	if (j > jt_max)
	    (void) printf("WARNING in prompt_for_window_params(), "
			  "unstable thermodynamics in window range\n");
	jt_max = min(j,nt-1);
	if ((Temp_max(seos) < Temp_min(seos)) || (jt_max < jt_min))
	{
	    screen("ERROR in prompt_for_window_params(), "
		   "invalid temperature range\n");
	    clean_up(ERROR);
	}

	screen("Enter the reference temperature(Kelvin) "
	       "(default = %g degrees Kelvin): ",Temp_ref(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Temp_ref(seos));

	s2eos_lookup(Rho_ref(seos),Temp_ref(seos),p,e,seos);
	csq = sesame_rt_sound_speed_squared(Rho_ref(seos),Temp_ref(seos),p,e);
	Reference_sound_speed(seos) = sqrt(csq);

	Pressure_min(seos) = HUGE_VAL;
	Pressure_max(seos) = -HUGE_VAL;
	Pressure_ref(seos) = p[0];

	Energy_min(seos) = HUGE_VAL;
	Energy_max(seos) = -HUGE_VAL;
	Energy_ref(seos) = e[0];

	Entropy_min(seos) = HUGE_VAL;
	Entropy_max(seos) = -HUGE_VAL;

	RS_entropy_scale(seos) = 1.0;
	PS_entropy_scale(seos) = 1.0;

	Gr.L[0] = ses_rt_grid_from_rho(Rho_min(seos),seos);
	Gr.U[0] = ses_rt_grid_from_rho(Rho_max(seos),seos);
	Gr.gmax[0] = 10*(ir_max - ir_min + 1);/*TOLERANCE*/
	Gr.L[1] = ses_rt_grid_from_temp(Temp_min(seos),seos);
	Gr.U[1] = ses_rt_grid_from_temp(Temp_max(seos),seos);
	Gr.gmax[1] = 10*(jt_max - jt_min + 1);/*TOLERANCE*/
	Gr.dim = 2;
	set_rect_grid(Gr.L,Gr.U,Gr.L,Gr.U,NOBUF,NOBUF,Gr.gmax,2,&Remap,&Gr);
	(void) adjust_top_grid_for_square(&Gr,&Gr);

	Nrho_hyp(seos) = Gr.gmax[0]; /*DEFAULT*/
	screen("Enter the number of density mesh points desired "
	       "(default = %d): ",Nrho_hyp(seos));
	(void) Gets(s);
	if (s[0] != '\0')
		(void) sscanf(s,"%d",&Nrho_hyp(seos));

	Ntemp_hyp(seos) = Gr.gmax[1]; /*DEFAULT*/
	screen("Enter the number of temperature mesh points desired "
	       "(default = %d): ",Ntemp_hyp(seos));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&Ntemp_hyp(seos));

	seos->abser0 = eps_mach; /*DEFAULT*/
	screen("Enter the absolute truncation error desired for computing "
	       "entropy\n\t(default= %g): ",seos->abser0);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&seos->abser0);
	if(seos->abser0 < eps_mach) 
	    seos->abser0 = eps_mach;

	seos->reler0 = MIN_REL_ERROR;
	screen("Enter the relative truncation error desired for computing "
	       "entropy\n\t(default= %g): ",seos->reler0);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&seos->reler0);
	if(seos->reler0 < MIN_REL_ERROR)
	    seos->reler0 = MIN_REL_ERROR;

	seos->terr0 = HUGE_VAL;

	seos->_set_RT_entropy = (zero_temperature_cold_curve(seos) == YES) ?
		set_RT_entropy_from_cold_curve : set_RT_entropy_from_mid_point;
	screen("Two methods are available to compute entropy.\n"
	       "\tIntegration from mid point of density/temperature window "
	       "(midpoint)\n"
	       "\tIntegration from cold curve (cold curve)\n"
	       "Enter choice (default = %s): ",
		(seos->_set_RT_entropy==set_RT_entropy_from_cold_curve) ?
		"cold curve" : "midpoint");
	(void) Gets(s);
	if (s[0] == 'm' || s[0] == 'M')
	    seos->_set_RT_entropy = set_RT_entropy_from_mid_point;
	else if (s[0] == 'c' || s[0] == 'C')
	    seos->_set_RT_entropy = set_RT_entropy_from_cold_curve;
}		/*end prompt_for_window_params*/

LOCAL	void	print_list_of_sesame_materials(
	FILE	*file,
	char	*filename)
{
	IO_TYPE    IO_Type;
	FILE	   *fp = fopen(filename,"r");
	char	   *sbuff = NULL;
	const char *month;
	char       *s;
	double      dtmp;
	double      *ids;
	double      date;
	double      *tabtype = NULL;
	int        year, day;
	int        idate;
	int        sbuff_size = 0;
	int        *tablen = NULL;
	int        max_num_tabs = 0;
	int        num_tabtypes;
	int        mat, i, num_materials;
	static const char *months[] = { "January",
				        "February",
				        "March",
				        "April",
				        "May",
				        "June",
				        "July",
				        "August",
				        "September",
				        "October",
				        "November",
				        "December"};

#define	skip_end_of_record(fp)	(void)fseek(fp,sizeof(double),SEEK_CUR)
#define	end_of_material(fp)	(void)fseek(fp,5,SEEK_CUR)
#define	start_fortran_read(fp)	(void)fseek(fp,INT,SEEK_SET)

	if (fp == NULL)
	{
	    screen("ERROR in print_list_of_sesame_materials(), "
	           "can't open %s\n",filename);
	    clean_up(ERROR);
	}
	determine_io_type(fp,&IO_Type);
	IO_Type.read_float_size = sizeof(double);
	IO_Type.read_endian = FT_BIG_ENDIAN;
	if (ft_endian_type() == FT_LITTLE_ENDIAN)
	{
	    IO_Type.reverse_endian = YES;
	}
	start_fortran_read(fp);
	read_binary_real_array(&dtmp,1,&IO_Type);
	num_materials = irint(dtmp);
	(void) fprintf(file,"Number of materials = %d\n",num_materials);
	read_binary_real_array(&date,1,&IO_Type);
	idate = (int)(date + 0.5);
	year = idate%100;
	year += (year < 70) ? 2000 : 1900;
	i = idate/100%100-1;
	if ((i < 0) || (i > 11))
	{
	    (void) fprintf(file,"Invalid date for table\n");
	}
	else
	{
	    month = months[i];
	    day = ((int)(date))/10000;
	    (void) fprintf(file,"Binary table produced on %s %d, %d\n",
	    	           month,day,year);
	}
	(void) fseek(fp,FLOAT,SEEK_CUR);

	skip_end_of_record(fp);

	uni_array(&ids,num_materials,FLOAT);
	read_binary_real_array(ids,num_materials,&IO_Type);
	(void)fseek(fp,(long)(2*FLOAT*num_materials),SEEK_CUR);

	skip_end_of_record(fp);
	end_of_material(fp);
	skip_end_of_record(fp);

	for (mat = 0; mat < num_materials; ++mat)
	{
	    (void)fseek(fp,FLOAT*4,SEEK_CUR);

	    read_binary_real_array(&dtmp,1,&IO_Type);
	    num_tabtypes = irint(dtmp);

	    if (num_tabtypes > max_num_tabs)
	    {
	    	if (tabtype != NULL)
	    	    free(tabtype);
	    	if (tablen != NULL)
	    	    free(tablen);
	    	max_num_tabs = num_tabtypes;
	    	uni_array(&tabtype,max_num_tabs,FLOAT);
	    	uni_array(&tablen,max_num_tabs,sizeof(int));
	    }
	    read_binary_real_array(tabtype,num_tabtypes,&IO_Type);

	    for (i = 0; i < num_tabtypes; ++i)
	    {
	        read_binary_real_array(&dtmp,1,&IO_Type);
	    	tablen[i] = irint(dtmp);
	    }

	    skip_end_of_record(fp);

	    if (tabtype[0] == 101.0)
	    {
	    	if ((tablen[0]*FLOAT+1) > sbuff_size)
	    	{
	    	    if (sbuff != NULL)
	    	    	free(sbuff);
	    	    sbuff_size = (int)(tablen[0]*FLOAT+1);
	    	    uni_array(&sbuff,sbuff_size,CHAR);
	    	}
	    	(void)fread((void*)sbuff,CHAR,tablen[0]*FLOAT,fp);
	    	sbuff[tablen[0]*FLOAT] = '\0';
	    	if ((s = strstr(sbuff,"/source")) != NULL)
	   	    s[0] = '\0';
	        (void)fprintf(file,"Id %g, %s\n",ids[mat],sbuff);
	    }
	    else
	    	(void)fseek(fp,(long)(FLOAT*tablen[0]),SEEK_CUR);

	    (void) fprintf(file,"Table types:");
	    for (i = 0; i < num_tabtypes; ++i)
	        (void) fprintf(file," %g",tabtype[i]);
	    (void) fprintf(file,"\n\n");

	    skip_end_of_record(fp);

	    for (i = 1; i < num_tabtypes; ++i)
	    {
	    	(void)fseek(fp,(long)(FLOAT*tablen[i]),SEEK_CUR);
	    	skip_end_of_record(fp);
	    }

	    end_of_material(fp);
	    skip_end_of_record(fp);
	}
	free_these(3,ids,tabtype,tablen,sbuff);
	(void) fclose(fp);
#undef	skip_end_of_record
#undef	end_of_material
#undef	start_fortran_read
}		/*end print_list_of_sesame_materials*/

/*
*		     sesinv():
*
*	Inverts the original sesame density-temperature tables.
*
*	ids2  (input) = Sesame material ID number
*	tbls (output) = array containing original SESAME tables
*	lcnt (output) = contains location of last entry in tbls.
*	lu    (input) = unit number of input file
*	ifl   (input) = error flag
*/

LOCAL	void sesinv(
	SESAME_EOS    *seos)
{
	POINTER pbdry = NULL;
	POINTER ccur = NULL;
	double	rgmax = ses_rt_grid_from_rho(Rho_max(seos),seos);
	double	rgmin = ses_rt_grid_from_rho(Rho_min(seos),seos);
	double	tgmax = ses_rt_grid_from_temp(Temp_max(seos),seos);
	double	tgmin = ses_rt_grid_from_temp(Temp_min(seos),seos);

	DEBUG_ENTER(sesinv)

	g_preserve_user_hooks(2,SAVE_HOOKS);/*no return before matching call */
	set_user_hooks_for_sesame();

	seos->de_params.rlgmid = 0.5*(rgmax+rgmin);
	seos->de_params.tlgmid = 0.5*(tgmax+tgmin);
	seos->de_params.rmid =
	    ses_rt_rho_from_grid(seos->de_params.rlgmid,seos);
	seos->de_params.tmid =
	    ses_rt_temp_from_grid(seos->de_params.tlgmid,seos);

#if defined(PHASE_CODE)
	if(multiphase_eos(seos) == YES)
	{
	    PHASE_BDRY *phase_bound;
	    COLD_CURVE *cold_curve;
	    double	*tbls = seos->sestab.tbls;
	    int nrho, ntemp;
	    int n_pts_dome;
	    int n_pts_cold;
	    int n_pts_gr;
	    int ntv;
	    int i;

	    seos->n_pts_ref_curve = 2*(Nrho_hyp(seos))+2;
	    uni_array(&seos->S_on_ref_curve,seos->n_pts_ref_curve,FLOAT);

	    nrho = (int) tbls[2];		ntemp = (int) tbls[3];
	    ntv = (int) tbls[4+nrho+ntemp+3*ntemp*nrho];

	    scalar(&phase_bound,sizeof(PHASE_BDRY));
	    phase_bound->n_pts_dome = n_pts_dome = 2*ntv;
	    debug_print("sesinv","In sesinv(), n_pts_dome = %d\n",n_pts_dome);

	    uni_array(&phase_bound->S,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->slphT,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->slphrp,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->slphe,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->re,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->rp,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->ce,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->cp,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->rvar,n_pts_dome,FLOAT);
	    uni_array(&phase_bound->Tvar,n_pts_dome,FLOAT);

	    for (i = 0; i < n_pts_dome; ++i)
	    {
	    	phase_bound->S[i] = 0.0;
	    	phase_bound->slphT[i] = 0.0;
	    	phase_bound->slphrp[i] = 0.0;
	    	phase_bound->slphe[i] = 0.0;
	    	phase_bound->re[i] = 0.0;
	    	phase_bound->rp[i] = 0.0;
	    	phase_bound->ce[i] = 0.0;
	    	phase_bound->cp[i] = 0.0;
	    	phase_bound->rvar[i] = 0.0;
	    	phase_bound->Tvar[i] = 0.0;
	    }

	    pbdry = (POINTER) phase_bound;

	    scalar(&cold_curve,sizeof(COLD_CURVE));

	    n_pts_cold = nrho + 4;
	    debug_print("sesinv","In sesinv(), n_pts_cold = %d\n",n_pts_cold);

	    uni_array(&cold_curve->slcp,n_pts_cold,FLOAT);
	    uni_array(&cold_curve->cp,n_pts_cold,FLOAT);
	    uni_array(&cold_curve->ce,n_pts_cold,FLOAT);
	    uni_array(&cold_curve->slce,n_pts_cold,FLOAT);
	    uni_array(&cold_curve->rvar,n_pts_cold,FLOAT);

	    for (i = 0; i < n_pts_cold; ++i)
	    {
	    	cold_curve->slcp[i] = 0.0;
	    	cold_curve->cp[i] = 0.0;
	    	cold_curve->ce[i] = 0.0;
	    	cold_curve->slce[i] = 0.0;
	    	cold_curve->rvar[i] = 0.0;
	    }

	    ccur = (POINTER) cold_curve;

	    n_pts_gr = 2*(nrho + 4)*ntemp;
	    debug_print("sesinv","In sesinv(), n_pts_gr = %d\n",n_pts_gr);

	    uni_array(&seos->slpe,n_pts_gr,FLOAT);
	    uni_array(&seos->slpp,n_pts_gr,FLOAT);

	    for (i = 0; i < n_pts_gr; ++i)
	    {
	    	seos->slpe[i] = 0.0;
	    	seos->slpp[i] = 0.0;
	    }

	    cold_PE_spline(seos,cold_curve);

	    phase_spline(seos,cold_curve,phase_bound);
	}
	else
#endif /* defined(PHASE_CODE) */
	{
	    double	p[3], e[3];

	    /* sets pmin = P(rhomin,Tmin) */
	    s2eos_lookup(Rho_min(seos),Temp_min(seos),p,e,seos);
	    seos->pmin = p[0];
	}

	init_sesame_hyp_tri_solns(seos,pbdry,ccur);

#if defined(PHASE_CODE)
	if (multiphase_eos(seos) == YES)
	{
	    free_these(3,seos->S_on_ref_curve,seos->slpe,seos->slpp);
	    seos->n_pts_ref_curve = 0;
	    seos->S_on_ref_curve = NULL;
	    seos->slpe = NULL;
	    seos->slpp = NULL;
	}

	if (pbdry != NULL)
	{
	    PHASE_BDRY *phase_bound = (PHASE_BDRY *) pbdry;

	    free_these(9,phase_bound->slphT,phase_bound->slphrp,
	    	       phase_bound->slphe,phase_bound->re,
		       phase_bound->rp,phase_bound->ce,
		       phase_bound->cp,phase_bound->rvar,
		       phase_bound->Tvar);
	    free(phase_bound);
	    pbdry = NULL;
	}
	if (ccur != NULL)
	{
	    COLD_CURVE *cold_curve = (COLD_CURVE *)ccur;

	    free_these(5,cold_curve->slcp,cold_curve->cp,
	    	       cold_curve->ce,cold_curve->slce,
	    	       cold_curve->rvar);
	    free(cold_curve);
	    ccur = NULL;
	}
#endif /* defined(PHASE_CODE) */
	g_preserve_user_hooks(2,RESTORE_HOOKS);
	DEBUG_LEAVE(sesinv)
}		/*end sesinv*/



LOCAL void open_fortran_device(
	char *filename,
	int lu)
{
	char s[120];

	(void) sprintf(s,"%s",filename);
#if defined(cray)
	SFORTRAN_NAME(oplib_uf)(&lu,_cptofcd(s,strlen(s)));
#else /* defined(cray) */
	/*SFORTRAN_NAME(oplib_uf)(&lu,s,(int)strlen(s)); */
#endif /* defined(cray) */
}		/*end open_fortran_device*/

/***************END INITIALIZATION UTILITY FUNCTIONS************************/


/*
*			int_dp_over_c_rho();
*
*	Evaluates the integral from p to p0 of dp/c*rho with constant entropy.
*/


LOCAL	double	int_dp_over_c_rho(
	double		p,
	Locstate	state0)
{
	double		p0 = pressure(state0);
	double		rho0 = Dens(state0);
	double		ans;
	double		coords[MAXD], S0, c, var[NUM_SES_VAR], rho1;
	double		cref;
	int		interp_flag;
	Front		*fr;
	Wave		*wv;
	COMPONENT	comp;
	SESAME_EOS	*seos = Ses_params(state0);
		
	cref = Reference_sound_speed(seos);
	S0 = entropy(state0);
	coords[1] = ses_ps_grid_from_entpy(S0,seos);
	coords[0] = ses_ps_grid_from_press(p0,seos);
	fr = seos->fr[SESAME_PRESS_ENTROPY];
	wv = seos->wave[SESAME_PRESS_ENTROPY];
	interp_flag = EVALUATE_ADB_GAMMA | EVALUATE_IDPOCR;
	set_ses_intrp_flag(interp_flag,SESAME_PRESS_ENTROPY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state0)),
		COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	c = 0.5*(sqrt(var[PS_AG]*p0/rho0)+cref);
	ans = c*var[PS_RF];


	coords[0] = ses_ps_grid_from_press(p,seos);
	interp_flag = EVALUATE_DENSITY | EVALUATE_ADB_GAMMA | EVALUATE_IDPOCR;
	set_ses_intrp_flag(interp_flag,SESAME_PRESS_ENTROPY);
	comp = nearest_interior_comp(multiphase_eos(Eos(state0)),
		COMP_PURE_PHASE,coords,fr->interf);
	ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,wv,var);
	rho1 = ses_ps_rho_from_var(var[PS_RHO],seos);
	c = 0.5*(sqrt(var[PS_AG]*p/rho1)+cref);
	return ans - c*var[PS_RF];
}		/*end int_dp_over_c_rho*/


#if defined(PHASE_CODE)
LOCAL	double	mass_flux_for_phase_change(
	double		p,	/* ratio of pressures across sound waves */
	Locstate	st0)
{
	WAVE_CURVE	*wave_cur = Wave_curve(st0);
	int		num = wave_cur->num_waves;
	int		i;
	double		p1, p2, rho1, V, V1, V0, ms;
	double		p0, rho0, rho;
	double		denom;

	rho0 = Dens(st0);
	p0 = pressure(st0);

	if (p < Press(wave_cur->st[0]))
	{
		switch (wave_cur->w_type[0])
		{
		case COMPOSITE:
			return mass_flux_for_comp(p,st0);

		case RAREFACTION:
			if (wave_cur->w_type[1] == RAREFACTION)
			{
				/* Split rarefaction */
				p1 = Press(wave_cur->st[0]);
				denom = int_dp_over_c_rho(p,wave_cur->st[1]);
				return (p0 - p)/denom;
			}
			else
			{
				/* Single rarefaction */
				denom = int_dp_over_c_rho(p,wave_cur->st[0]);
				return (p0 - p)/denom;
			}
		}
	}
	for (i = 0; i < num-1; ++i)
	{
		p1 = Press(wave_cur->st[i]);
		p2 = Press(wave_cur->st[i+1]);
		if (Between(p,p1,p2))
		{
			switch (wave_cur->w_type[i])
			{
			case RAREFACTION:
				denom = int_dp_over_c_rho(p,st0);
				return (p0 - p)/denom;
			case COMPOSITE:
				switch (wave_cur->w_type[i+1])
				{
				case COMPOSITE:
					/* Composite */
					return mass_flux_for_comp(p,st0);

				case RAREFACTION:
					/* Split rarefaction */
					denom = int_dp_over_c_rho(p,st0);
					return (p0 - p)/denom;
				}
				break;
			case SHOCK:
				/* Check if split shock */
				if (wave_cur->w_type[i-1] != SHOCK)
				{
					/* Simple shock */
					V = 1.0/dens_Hugoniot(p,st0);
					V0 = 1.0/rho0;
					return ((p-p0)< p0*EPS) ?
						acoustic_impedance(st0) :
						sqrt((p - p0) / (V0 - V));
				}
				else
				{
					/* Split shock */
					rho = dens_Hugoniot(p,wave_cur->st[i]);
					V = 1.0/rho;
					rho1 = Dens(wave_cur->st[i]);
					V1 = 1.0/rho1;
					ms = (p - p1)/sqrt((p1 - p)/(V - V1));
					V0 = 1.0/rho0;
					ms = ms + sqrt((p1 - p0)*(V1 - V0));
				}
			}
		}
	}
	if ( p > Press(wave_cur->st[num-1]))
	{
		if (wave_cur->w_type[num - 2] == SHOCK && num-3 >= 0
			 && wave_cur->w_type[num-3] != SHOCK)
		{
			/* Split shock */
			/* TODO: This may not be robust */	
			rho = dens_Hugoniot(p,wave_cur->st[num-1]);
			V = 1.0/rho;
			p1 = Press(wave_cur->st[num-1]);
			rho1 = Dens(wave_cur->st[num-1]);
			V1 = 1.0/rho1;
			ms = (p - p1)/sqrt((p1 - p)/(V - V1));
			V0 = 1.0/rho0;
			ms = ms + sqrt((p1 - p0)*(V1 - V0));
			return (p - p0)/ms;

		}
		else
		{
				/* Simple shock */
			V  = 1.0/dens_Hugoniot(p,st0);
			V0 = 1.0/rho;
			return ((p-p0)< p0*EPS) ?
				acoustic_impedance(st0) :
				sqrt((p - p0) / (V0 - V));
		}
	}
	return ERROR_FLOAT;
}		/*end mass_flux_for_phase_change*/

LOCAL	double	mass_flux_for_comp(
	double		p,
	Locstate	st0)
{
	WAVE_CURVE	*wave_cur;
	int		i, wtp, num;
	double		ms, p0, slope;
	double		pend, pstart, rhoc;
	double		pend1, pstart1, rhoc1;

	wave_cur = Wave_curve(st0);
	ms = 0.0;
	p0 = pressure(st0);
	num = wave_cur->num_waves;
	num = num - 1;
	for (i = 0; i < num; ++i)
	{
	    if (p < Press(wave_cur->st[i]))
	    {
	        wtp = wave_cur->w_type[i+1];
	        if (wtp == COMPOSITE)
	        {
	    	    ms = int_dp_over_c_rho(p,wave_cur->st[i]);
	    	    pend = wave_cur->pend[NO_PTS_ON_COMP-1];
	    	    pstart = wave_cur->pstart[NO_PTS_ON_COMP-1];
	    	    rhoc = wave_cur->rhoc[NO_PTS_ON_COMP-1];
	    	    ms = ms + (pend - pstart)/rhoc; 
	    	    ms = ms + int_dp_over_c_rho(pstart,st0);
	    	    return (p - p0)/ms;
	        }
	    }
	}
	for (i = 0; i < NO_PTS_ON_COMP - 1; ++i)
	{
	/* Find mass flux on composite portion through linear interpolation */	
		pend = wave_cur->pend[i];
		pend1 = wave_cur->pend[i+1];
		if (Between(p,pend1,pend))
		{
			rhoc = wave_cur->rhoc[i];
			rhoc1 = wave_cur->rhoc[i+1];
			pstart = wave_cur->pstart[i];
			pstart1 = wave_cur->pstart[i+1];
			slope = (p - pend)/(pend1 - pend);
			rhoc = (rhoc1 - rhoc)*slope + rhoc;
			pstart = (pstart1 - pstart)*slope + pstart;
			break;
		}
	}
	ms = (p - pstart)/rhoc;
	ms = ms + int_dp_over_c_rho(pstart,st0);
	return  (p - p0)/ms;
}		/*end mass_flux_for_comp*/
#endif /* defined(PHASE_CODE) */

EXPORT	S2DIR	FORTRAN_NAME(s2dir);

LOCAL	boolean	initialize_Sesame_tables(
	SESAME_EOS *seos)
{
	IMPORT	S2DIR	FORTRAN_NAME(s2dir);
	boolean	status;
	int	ids2 = seos->ids2;
	int	ifl, ir = 1;
	int	lcnt;
	int	lu = 9;
	static	S2DIR	S2dflt = {
				    10000,/*default dimension*/
				    1,
				    0,0,0,0,0,0,0,0,0,0
				 };

	debug_print("init","Initializing a new sesame table\n");

	FORTRAN_NAME(s2dir) = S2dflt;
	seos->sestab.s2dir = &FORTRAN_NAME(s2dir);
	uni_array(&seos->sestab.tbls,FORTRAN_NAME(s2dir).lcmx,FLOAT);
	if (seos->sestab.tbls == NULL)
	{
	    screen("ERROR: in initialize_Sesame_tables(), "
	           "Unable to allocate table storage\n");
	    clean_up(ERROR);
	}

	lcnt = 1;
	open_fortran_device(seos->seslib,lu);

	GetSesameTable(seos,&ir,&ids2,seos->sestab.tbls,&lcnt,&lu,&ifl);

	if (ifl < 0)
	{
	    debug_print("init","Not enough storage in tables, reallocating\n");
	    free(seos->sestab.tbls);
	    FORTRAN_NAME(s2dir) = S2dflt;
	    FORTRAN_NAME(s2dir).lcmx -= ifl;
	    uni_array(&seos->sestab.tbls,FORTRAN_NAME(s2dir).lcmx,FLOAT);
	    if (seos->sestab.tbls == NULL)
	    {
	        screen("ERROR: in initialize_Sesame_tables(), "
	               "Unable to allocate more storage\n");
	        clean_up(ERROR);
	    }
	    GetSesameTable(seos,&ir,&ids2,seos->sestab.tbls,&lcnt,&lu,&ifl);
	}
	if (ifl == 0)
	{
	    screen("Material %d not found.\n",ids2);
	    status = FUNCTION_FAILED;
	}
	else
	    status = FUNCTION_SUCCEEDED;

	FORTRAN_NAME(cllib)(&lu);
	return status;
}		/*end initialize_Sesame_tables*/

#endif /* defined(SESAME_CODE) && defined(TWOD) */
