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
*				generic-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Generic file for the implementation of an equation of state model
*
*/

#define	DEBUG_STRING	"generic"
#include <geos/geosdecs.h>

/* Generic EOS prototypes */

	/* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL	double	GENERIC_internal_energy(Locstate);
LOCAL	double	GENERIC_pressure(Locstate);
LOCAL   double   GENERIC_density(Locstate);
LOCAL	double	GENERIC_sound_speed_squared(Locstate);
LOCAL	double	GENERIC_acoustic_impedance_squared(Locstate);
LOCAL	double	GENERIC_specific_internal_energy(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL	double	GENERIC_specific_enthalpy(Locstate);
LOCAL	double	GENERIC_temperature(Locstate);
LOCAL	double	GENERIC_entropy(Locstate);
LOCAL	double	GENERIC_adiabatic_gamma(Locstate);
LOCAL	double	GENERIC_gruneisen_gamma(Locstate);
LOCAL	double	GENERIC_fundamental_derivative(Locstate);
LOCAL	double	GENERIC_C_V(Locstate);
LOCAL	double	GENERIC_C_P(Locstate);
LOCAL	double	GENERIC_K_S(Locstate);
LOCAL	double	GENERIC_K_T(Locstate);

	/* MATERIAL PROPERTY FUNCTIONS */
LOCAL	double	GENERIC_bulk_viscosity(Locstate);
LOCAL	double	GENERIC_shear_viscosity(Locstate);
LOCAL 	double   GENERIC_heat_coeff(Locstate);

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
LOCAL	void    GENERIC_single_eos_load_pressure_and_sound_speed2(Vec_Gas*,
								  int,int);
LOCAL	void	GENERIC_single_eos_load_pressure_and_gammas(Vec_Gas*,int,int);
LOCAL	void	GENERIC_single_eos_load_pressure(Vec_Gas*,int,int);
LOCAL	void	GENERIC_single_eos_load_sound_speed2(Vec_Gas*,int,int);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
LOCAL	double	GENERIC_dens_Hugoniot(double,Locstate);
LOCAL	void	GENERIC_state_w_pr_on_Hugoniot(Locstate,double,Locstate,int);
LOCAL	boolean	GENERIC_state_w_mf_sqr_on_Hugoniot(Locstate,double,Locstate,int);

	/* Velocity Related Hugoniot Functions */
LOCAL	double	GENERIC_pr_normal_vel_wave_curve(double,Locstate);

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL	double	GENERIC_dens_rarefaction(double,Locstate);
LOCAL	double	GENERIC_pressure_rarefaction(double,Locstate);
LOCAL	void	GENERIC_state_on_adiabat_with_pr(Locstate,double,Locstate,int);
LOCAL	void	GENERIC_state_on_adiabat_with_dens(Locstate,double,Locstate,int);

	/* General Wave Curve Functions */
LOCAL	double	GENERIC_mass_flux(double,Locstate);
LOCAL	double	GENERIC_mass_flux_squared(double,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
LOCAL	double	GENERIC_oned_fan_state(double,Locstate,Locstate,Locstate,
				       int,boolean*);

	/* Functions to Compute Riemann Solutions */
LOCAL	double	GENERIC_riemann_wave_curve(Locstate,double);
LOCAL	void	GENERIC_set_state_for_find_mid_state(Locstate,Locstate);
LOCAL	double	GENERIC_eps_for_Godunov(Locstate,double,double);
LOCAL	void	GENERIC_initialize_riemann_solver(Locstate,Locstate,double*,
						  double*,double,double*,double*,
						  boolean(*)(Locstate,Locstate,
							  double,double*,double*,
							  double*,double*,double*,
							  double*,
							  RIEMANN_SOLVER_WAVE_TYPE*,
							  RIEMANN_SOLVER_WAVE_TYPE*));


	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
LOCAL	boolean	GENERIC_steady_state_wave_curve(double,double,double*,Locstate);
LOCAL	double	GENERIC_pressure_at_sonic_point(double,Locstate);
LOCAL	boolean	GENERIC_pr_at_max_turn_angle(double*,double,Locstate);
LOCAL	double	GENERIC_state_in_prandtl_meyer_wave(double,double,Locstate,double,
						    Locstate,Locstate,int);

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
LOCAL	double	GENERIC_CJ_state(Locstate,int,Locstate,int,int);
LOCAL	void	GENERIC_progress_state(double,Locstate,Locstate,double);
LOCAL	void	GENERIC_fprint_combustion_params(FILE*,Gas_param*);
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
LOCAL	void	GENERIC_neumann_riem_inv_moc(double*,Locstate,double,double,
	                                     Locstate,SIDE,Locstate,double,
	                                     double*,Front*);
LOCAL	void	GENERIC_shock_ahead_state_riem_inv_moc(double*,Locstate,
						       Locstate,Locstate,
						       Locstate,Locstate,
						       double,double,double,double,
						       double*,double*,int,
						       double,Front*);
LOCAL	boolean	GENERIC_shock_moc_plus_rh(double*,Locstate,Locstate,Locstate,
					  Locstate,double,double*,double*,int,
					  Front*);

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	GENERIC_prompt_for_state(Locstate,int,Gas_param*,const char*);
LOCAL	void	GENERIC_prompt_for_thermodynamics(Locstate,Gas_param*,
						  const char*);
LOCAL	void	GENERIC_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	GENERIC_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                              Gas_param*);
LOCAL	EOS*	GENERIC_free_EOS_params(EOS*);
LOCAL	void	GENERIC_prompt_for_EOS_params(INIT_DATA*,Gas_param*,
					      const char*,const char*);

	/* Problem Type Specific Initialization Functions */
LOCAL	double	GENERIC_RT_RS_f(double,Locstate,double,double,double);
LOCAL	void	GENERIC_RT_single_mode_perturbation_state(Locstate,
							  double*,double,
							  Locstate,double,
							  double,MODE*,double);
LOCAL	void	GENERIC_KH_single_mode_state(Locstate,double*,double,Locstate,
					     double,double,double,MODE*);
LOCAL	void	GENERIC_compute_isothermal_stratified_state(Locstate,double,
							    double,Locstate);
LOCAL	void	GENERIC_compute_isentropic_stratified_state(Locstate,double,
							    double,Locstate);
LOCAL	void	GENERIC_compute_constant_density_stratified_state(Locstate,
							          double,double,
								  Locstate);

	/* Equation of state domain functions */
LOCAL	double	GENERIC_Min_energy(Locstate);
LOCAL	double	GENERIC_Min_pressure(Locstate);
LOCAL	double	GENERIC_Vacuum_dens(Locstate);
LOCAL	double	GENERIC_Raref_press(Locstate);
#if defined(COMBUSTION_CODE)
LOCAL	double	GENERIC_Tol_alpha(Locstate);
LOCAL	double	GENERIC_Tol_pressure(Locstate);
#endif /* defined(COMBUSTION_CODE) */

	/* LOCAL function prototypes */
LOCAL	boolean	Mach_number_squared_behind_oblique_shock(double,double*,Locstate);
LOCAL	boolean	mta_aux(double,double*,Locstate);
LOCAL	boolean	oned_fan_aux(double,double*,POINTER);
LOCAL	boolean	pmsr(double,double*,double*,int,Locstate);
LOCAL	boolean	pr_normal_vel_wave_curve_aux(double,double*,Locstate);
LOCAL	boolean	prandtl_meyer_speed_rinv(double*,double,Locstate,double*,double,
					 Locstate,double*);
LOCAL	boolean	pmw_aux(double,double*,Locstate);
LOCAL	boolean	rndr(double,double*,double*,int,Locstate);
LOCAL	boolean	rnpr(double,double*,double*,int,Locstate);
LOCAL	double	int_dp_over_rho_c(double,Locstate,Locstate);
LOCAL	void	set_eos_function_hooks(EOS*);
LOCAL   double   gaussian_int(int,double,double,Locstate);
LOCAL   void    legendre_init(int,double*,double*);
LOCAL   void    linear_transform(int,double*,double*,double,double);
LOCAL   void    exp_transform(int,double*,double*);
LOCAL   double   one_over_rho_c(double,Locstate);

EXPORT	EOS	*set_GENERIC_eos(
	EOS	*eos)
{
	if (eos == NULL)
	    scalar(&eos,sizeof(EOS));
	set_eos_function_hooks(eos);
	return eos;
}		/*end set_GENERIC_eos*/

/*
*		BASIC EQUATION OF STATE FUNCTIONS
*		NO DEFAULT IMPLEMENTATION AVAILABLE
*/

/*
*			GENERIC_pressure():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
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

LOCAL	double	GENERIC_pressure(
	Locstate state)
{
	double pr;
	if (is_obstacle_state(state))
	    return HUGE_VAL;
	switch (state_type(state)) 
	{
	case	GAS_STATE:
	case	EGAS_STATE:
	case	FGAS_STATE:
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (Dens(state) < Vacuum_dens(state))
	    	return Min_pressure(state);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    pr = ERROR_FLOAT;
	    screen("ERROR in GENERIC_pressure(), function unimplemented\n");
	    clean_up(ERROR);
	    break;

	case	TGAS_STATE:
	case	VGAS_STATE:
	    pr = Press(state);
	    break;

	default:
	    screen("ERROR in GENERIC_pressure(), no such state type\n");
	    clean_up(ERROR);
	    break;
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(state));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end GENERIC_pressure*/

LOCAL   double   GENERIC_density(
        Locstate state)
{
        switch (state_type(state))
	{
	case	GAS_STATE:
	case	EGAS_STATE:
	case	FGAS_STATE:
	case	TGAS_STATE:
	case	VGAS_STATE:
	        return Dens(state);	
	default:
	        screen("ERROR in GENERIC_density(), function unimplemented\n");
		clean_up(ERROR);
		break;
	}
	return ERROR_FLOAT;
}               /*end GENERIC_density*/

/*
*			GENERIC_sound_speed_squared():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

/*ARGSUSED*/
LOCAL	double	GENERIC_sound_speed_squared(
	Locstate state)
{
	screen("ERROR in GENERIC_sound_speed_squared(), "
	       "function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end GENERIC_sound_speed_squared*/

/*
*			GENERIC_specific_internal_energy():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/

LOCAL	double	GENERIC_specific_internal_energy(
	Locstate state)
{
	double rho = Dens(state);

	switch (state_type(state))
	{

	case	GAS_STATE:
	    return (Energy(state) - kinetic_energy(state))/rho;

	case	EGAS_STATE:
	    return Energy(state);

	case	TGAS_STATE:
	case	FGAS_STATE:
	    screen("ERROR in GENERIC_specific_internal_energy(), "
	           "function unimplemented\n");
	    clean_up(ERROR);
	    return ERROR_FLOAT;
	
	case	VGAS_STATE:
	    return Int_en(state);

	default:
	    screen("ERROR in GENERIC_specific_internal_energy(), "
	           "no such state type\n");
	    clean_up(ERROR);
	    break;
	}
	return ERROR_FLOAT;
}		/*end GENERIC_specific_internal_energy*/

/*
*			GENERIC_gruneisen_gamma():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
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

/*ARGSUSED*/
LOCAL	double	GENERIC_gruneisen_gamma(
	Locstate state)
{
	screen("ERROR in GENERIC_gruneisen_gamma(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end GENERIC_gruneisen_gamma*/

/*
*			GENERIC_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_fprint_EOS_params(
	FILE *file,
	Gas_param *params)
{
	screen("ERROR in GENERIC_fprint_EOS_params(), no generic printing "
	       "of EOS parameters");
	clean_up(ERROR);
}		/*end GENERIC_fprint_EOS_params */

/*
*			GENERIC_read_print_EOS_params():
*
*	Reads the equation of state data as printed by GENERIC_fprint_EOS_params.
*	This is restart function.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_read_print_EOS_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	screen("ERROR in GENERIC_read_print_EOS_params(), no generic reading "
	       "of EOS parameters");
	clean_up(ERROR);
}		/*end GENERIC_read_print_EOS_params*/

/*
*			GENERIC_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_prompt_for_EOS_params(
	INIT_DATA  *init,
	Gas_param  *params,
	const char *message1,
	const char *message2)
{
	screen("ERROR in GENERIC_prompt_for_EOS_params(), no generic "
	       "prompting for EOS parameters\n");
	clean_up(ERROR);
}		/*end GENERIC_prompt_for_EOS_params*/
/*		END OF BASIC EQUATION OF STATE FUNCTIONS		*/

/*
*		BASIC TEMPERATURE DEPENDENT FUNCTIONS
*		NO DEFAULT IMPLEMENTATION AVAILABLE,
*		BUT MY NOT BE NEEDED FOR ALL APPLICATIONS
*/

/*
*			GENERIC_temperature():
*
*		BASIC TEMPERATURE DEPENDENT FUNCTION
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	double	GENERIC_temperature(
	Locstate state)
{
	if (state_type(state) == FGAS_STATE)
	    return Temperature(state);

	screen("ERROR in GENERIC_temperature(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end GENERIC_temperature*/

/*
*			GENERIC_entropy():
*
*		BASIC TEMPERATURE DEPENDENT FUNCTION
*
*	Returns the specific entropy of a state.
*/

LOCAL	double	GENERIC_entropy(
	Locstate state)
{
	if (state_type(state) == VGAS_STATE)
	    return Entropy(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(state) < Vacuum_dens(state))
	    return 0.0;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	screen("ERROR in GENERIC_entropy(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end GENERIC_entropy*/

/*
*			GENERIC_C_V():
*
*		BASIC TEMPERATURE DEPENDENT FUNCTION
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*/

/*ARGSUSED*/
LOCAL	double	GENERIC_C_V(
	Locstate state)
{
	screen("ERROR in GENERIC_C_V(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}	/* end GENERIC_C_V */

/*		END OF BASIC TEMPERATURE DEPENDENT FUNCTIONS		*/

LOCAL	void	set_eos_function_hooks(
	EOS *eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_internal_energy = GENERIC_internal_energy;
	eos->_pressure = GENERIC_pressure;
	eos->_density = GENERIC_density;
	eos->_sound_speed_squared = GENERIC_sound_speed_squared;
	eos->_acoustic_impedance_squared = GENERIC_acoustic_impedance_squared;
	eos->_specific_internal_energy = GENERIC_specific_internal_energy;
	
	/* MATERIAL PROPERTY FUNCTIONS */
	eos->_bulk_viscosity = GENERIC_bulk_viscosity;
	eos->_shear_viscosity = GENERIC_shear_viscosity;
	eos->_heat_coeff = GENERIC_heat_coeff;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_specific_enthalpy = GENERIC_specific_enthalpy;
	eos->_temperature = GENERIC_temperature;
	eos->_entropy = GENERIC_entropy;
	eos->_adiabatic_gamma = GENERIC_adiabatic_gamma;
	eos->_gruneisen_gamma = GENERIC_gruneisen_gamma;
	eos->_fundamental_derivative = GENERIC_fundamental_derivative;
	eos->_C_V = GENERIC_C_V;
	eos->_C_P = GENERIC_C_P;
	eos->_K_S = GENERIC_K_S;
	eos->_K_T = GENERIC_K_T;

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
	eos->_single_eos_load_pressure_and_sound_speed2 =
	    GENERIC_single_eos_load_pressure_and_sound_speed2;
	eos->_single_eos_load_pressure_and_gammas =
	    GENERIC_single_eos_load_pressure_and_gammas;
	eos->_single_eos_load_pressure = GENERIC_single_eos_load_pressure;
	eos->_single_eos_load_sound_speed2 =
	    GENERIC_single_eos_load_sound_speed2;

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
	eos->_dens_Hugoniot = GENERIC_dens_Hugoniot;
	eos->_state_w_pr_on_Hugoniot = GENERIC_state_w_pr_on_Hugoniot;
	eos->_state_w_mf_sqr_on_Hugoniot = GENERIC_state_w_mf_sqr_on_Hugoniot;

	/* Velocity Related Hugoniot Functions */
	eos->_pr_normal_vel_wave_curve = GENERIC_pr_normal_vel_wave_curve;

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	eos->_dens_rarefaction = GENERIC_dens_rarefaction;
	eos->_pressure_rarefaction = GENERIC_pressure_rarefaction;
	eos->_state_on_adiabat_with_pr = GENERIC_state_on_adiabat_with_pr;
	eos->_state_on_adiabat_with_dens = GENERIC_state_on_adiabat_with_dens;

	/* General Wave Curve Functions */
	eos->_mass_flux = GENERIC_mass_flux;
	eos->_mass_flux_squared = GENERIC_mass_flux_squared;

	/* Functions for the Evaluation of Riemann Solutions */
	eos->_oned_fan_state = GENERIC_oned_fan_state;

	/* Functions to Compute Riemann Solutions */
	eos->_riemann_wave_curve = GENERIC_riemann_wave_curve;
	eos->_set_state_for_find_mid_state =
	    GENERIC_set_state_for_find_mid_state;
	eos->_eps_for_Godunov = GENERIC_eps_for_Godunov;
	eos->_initialize_riemann_solver = GENERIC_initialize_riemann_solver;

	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
	eos->_steady_state_wave_curve = GENERIC_steady_state_wave_curve;
	eos->_pressure_at_sonic_point = GENERIC_pressure_at_sonic_point;
	eos->_pr_at_max_turn_angle = GENERIC_pr_at_max_turn_angle;
	eos->_state_in_prandtl_meyer_wave = GENERIC_state_in_prandtl_meyer_wave;

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
	eos->_CJ_state = GENERIC_CJ_state;
	eos->_progress_state = GENERIC_progress_state;
	eos->_fprint_combustion_params = GENERIC_fprint_combustion_params;
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
	eos->_neumann_riem_inv_moc = GENERIC_neumann_riem_inv_moc;
	eos->_shock_ahead_state_riem_inv_moc =
	    GENERIC_shock_ahead_state_riem_inv_moc;
	eos->_shock_moc_plus_rh = GENERIC_shock_moc_plus_rh;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_prompt_for_state = GENERIC_prompt_for_state;
	eos->_prompt_for_thermodynamics = GENERIC_prompt_for_thermodynamics;
	eos->_fprint_EOS_params = GENERIC_fprint_EOS_params;
	eos->_read_print_EOS_params = GENERIC_read_print_EOS_params;
	eos->_free_EOS_params = GENERIC_free_EOS_params;
	eos->_prompt_for_EOS_params = GENERIC_prompt_for_EOS_params;

	/* Problem Type Specific Initialization Functions */
	eos->_RT_RS_f = GENERIC_RT_RS_f;
	eos->_RT_single_mode_perturbation_state =
	    GENERIC_RT_single_mode_perturbation_state;
	eos->_KH_single_mode_state = GENERIC_KH_single_mode_state;
	eos->_compute_isothermal_stratified_state =
	    GENERIC_compute_isothermal_stratified_state;
	eos->_compute_isentropic_stratified_state =
	    GENERIC_compute_isentropic_stratified_state;
	eos->_compute_constant_density_stratified_state = GENERIC_compute_constant_density_stratified_state;
	eos->_multiphase_eos = NO;
	eos->_compute_ns_terms = NO;

	/* Equation of state domain functions */
	eos->_Min_energy = GENERIC_Min_energy;
	eos->_Min_pressure = GENERIC_Min_pressure;
	eos->_Vacuum_dens = GENERIC_Vacuum_dens;
	eos->_Raref_press = GENERIC_Raref_press;
#if defined(COMBUSTION_CODE)
	eos->_Tol_alpha = GENERIC_Tol_alpha;
	eos->_Tol_pressure = GENERIC_Tol_pressure;
#endif /* defined(COMBUSTION_CODE) */
}


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*			GENERIC_internal_energy():
*
*	Returns the internal energy per unit volume of a state.
*/

LOCAL	double	GENERIC_internal_energy(
	Locstate state)
{
	switch (state_type(state)) 
	{
	case GAS_STATE:
	    return	Energy(state) - kinetic_energy(state);

	case EGAS_STATE:
	    return	Energy(state)*Dens(state);

	case VGAS_STATE:
	    return Dens(state) * Int_en(state);

	case FGAS_STATE:
	case TGAS_STATE:
	    return Dens(state)*specific_internal_energy(state);

	default:
	    screen("ERROR: in GENERIC_internal_energy(), no such state type\n");
	    clean_up(ERROR);
	}
	return ERROR_FLOAT;
}		/*end GENERIC_internal_energy*/




/*
*		GENERIC_acoustic_impedance_squared():
*
*	Returns the square of the local acoustic impedance of the state.
*
*                        2     dP  |
*			i = - ---- |
*                              dV  |S
*/

LOCAL	double	GENERIC_acoustic_impedance_squared(
	Locstate state)
{
	return sqr(Dens(state))*sound_speed_squared(state);
}		/*end GENERIC_acoustic_impedance_squared*/



/***************END PRIMARY THERMODYNAMIC FUNCTIONS ************************/
/***************SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS ***********/

/*
*			GENERIC_specific_enthalpy():
*
*	This function computes the specific enthalpy of the given state.
*
*			H = E + P*V
*
*	E = specific internal energy, P = pressure, V = specific volume.
*
*/

LOCAL	double	GENERIC_specific_enthalpy(
	Locstate state)
{
#if defined(VERBOSE_PLUS_GAS)
	if (state_type(state) == VGAS_STATE)
	    return Enthalpy(state);
#endif /* defined(VERBOSE_PLUS_GAS) */

	return specific_internal_energy(state) + pressure(state)/Dens(state);
}		/*end GENERIC_specific_enthalpy*/


/*
*			GENERIC_adiabatic_gamma():
*
*	Returns the dimensionless sound speed
*
*		gamma = - d(log P)/d(log V) | .
*					     S
*	As usual P = thermodynamic pressure,  V = specific volume
*	and S = specific entropy.
*/

LOCAL	double	GENERIC_adiabatic_gamma(
	Locstate state)
{
	return sound_speed_squared(state)*Dens(state)/pressure(state);
}		/*end GENERIC_adiabatic_gamma*/


/*
*			GENERIC_fundamental_derivative():
*
*	Returns the fundamental derivative of gas dynamics for the state.
*	This quantity is defined by the formula
*
*			    2      2
*		           d P / dV  |
*                                    |S
*             G = -0.5 V -----------------
*                          dP / dV |
*                                  |S
*
*	Where P is the thermodynamic pressure,  V is the specific volume
*	and S is the specific entropy.  Both derivatives are taken at
*	constant S.
*/

LOCAL	double	GENERIC_fundamental_derivative(
	Locstate state)
{
	double rho = Dens(state);
	double cp2, cm2;
	double p = pressure(state);
	double dp = 0.01*p;/*TOLERANCE*/
	static Locstate tmpst = NULL;

	if (tmpst == NULL)
	    (*Params(state)->_alloc_state)(&tmpst,Params(state)->sizest);

	state_on_adiabat_with_pr(state,p+dp,tmpst,TGAS_STATE);
	cp2 = sound_speed_squared(tmpst);
	state_on_adiabat_with_pr(state,p-dp,tmpst,TGAS_STATE);
	cm2 = sound_speed_squared(tmpst);
	return 1.0 + 0.25*rho*(cp2 - cm2)/dp;
}		/*end GENERIC_fundamental_derivative*/

/*
*			GENERIC_C_P():
*
*	Specific heat at constant pressure.
*
*
*                        dS  |
*		C_P = T ---- |
*                        dT  | P
*/

LOCAL	double	GENERIC_C_P(
	Locstate state)
{
	double c2 = sound_speed_squared(state);
	double GAM = gruneisen_gamma(state);
	double T = temperature(state);
	double c_v = C_V(state);

	return c2*c_v/(c2 - GAM*GAM*T*c_v);
}	/* end GENERIC_C_P */

/*
*			GENERIC_K_S():
*
*	Isentropic compressibility.
*
*                        1   dV  |
*		K_S = - --- ---- |
*                        V   dP  | S
*/

LOCAL	double	GENERIC_K_S(
	Locstate state)
{
	return Dens(state)/acoustic_impedance_squared(state);
}	/* end GENERIC_K_S */

/*
*			GENERIC_K_T():
*
*	Isothermal compressibility.
*
*                        1   dV  |
*		K_T = - --- ---- |
*                        V   dP  | T
*/

LOCAL	double	GENERIC_K_T(
	Locstate state)
{
	double c2 = sound_speed_squared(state);
	double GAM = gruneisen_gamma(state);
	double T = temperature(state);
	double rho = Dens(state);
	double c_v = C_V(state);

	return 1.0/(rho*(c2 - GAM*GAM*T*c_v));
}	/*end GENERIC_K_T */



/************** END SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS *******/

/************** MATERIAL PROPERTY FUNCTIONS ********************************/
/*ARGSUSED*/
LOCAL	double	GENERIC_shear_viscosity(
	Locstate state)
{
	return 0.0;
}	/*end GENERIC_shear_viscosity */

/*ARGSUSED*/
LOCAL	double	GENERIC_bulk_viscosity(
	Locstate state)
{
	return 0.0;
}	/*end GENERIC_bulk_viscosity */

/*ARGSUSED*/
LOCAL   double   GENERIC_heat_coeff(
        Locstate state)
{
	return 0.0;
}               /*end GENERIC_heat_coeff*/
/************** END MATERIAL PROPERTY FUNCTIONS ****************************/

/************** VECTORIZED THERMODYNAMIC FUNCTIONS *************************/

/*
*		GENERIC_single_eos_load_pressure_and_sound_speed2():
*
*	Loads a vector of pressures and sound speeds into the
*	appropriate fields of the Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	GENERIC_single_eos_load_pressure_and_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	static Locstate tmpst = NULL;
	Locstate        *state = vst->state + offset;
	double           *rho = vst->rho + offset;
	double           *p = vst->p + offset, *c2 = vst->c2 + offset;
	double           *e = vst->e + offset;
	double           *FD = NULL;
	int             k;

	if (vst->FD != NULL)
	    FD = vst->FD + offset;

	if (tmpst == NULL)
	{
	    (*Params(state[0])->_alloc_state)(&tmpst,Params(state[0])->sizest);
	    zero_scalar(tmpst,Params(state[0])->sizest);
	    set_type_of_state(tmpst,EGAS_STATE);
	}
	for (k = 0; k < vsize; ++k)
	{
	    Dens(tmpst) = rho[k];
	    Energy(tmpst) = e[k];
	    Set_params(tmpst,state[k]);
	    p[k] = pressure(tmpst);
	    c2[k] = sound_speed_squared(tmpst);
	    if (FD != NULL)
	        FD[k] = fundamental_derivative(tmpst);
	}
}		/*end GENERIC_single_eos_load_pressure_and_sound_speed2*/


/*
*		GENERIC_single_eos_load_pressure_and_gammas():
*
*	Loads the pressure, adiabatic exponent, and Gruneisen
*	coefficient uni_arrays of the Vec_Gas state vst.
*	This function assumes that the specific internal energy
*	uni_array vst->e is already loaded.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_gammas.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	GENERIC_single_eos_load_pressure_and_gammas(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	static Locstate tmpst = NULL;
	Locstate        state = vst->state[offset];
	double           *rho = vst->rho + offset;
	double           *p = vst->p + offset;
	double           *e = vst->e + offset;
	double           *c2 = vst->c2 + offset, *GAM = vst->GAM + offset;
	double           *FD = NULL;
	int             k;

	if (vst->FD != NULL)
	    FD = vst->FD + offset;

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
	    p[k] = pressure(tmpst);
	    c2[k] = sound_speed_squared(tmpst);
	    GAM[k] = gruneisen_gamma(tmpst);
	    if (FD != NULL)
	        FD[k] = fundamental_derivative(tmpst);
	}
}		/*end GENERIC_single_eos_load_pressure_and_gammas*/

/*
*			GENERIC_single_eos_load_pressure():
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

/*ARGSUSED*/
LOCAL	void	GENERIC_single_eos_load_pressure(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	static Locstate tmpst = NULL;
	Locstate *state = vst->state + offset;
	double *rho = vst->rho + offset;
	double *p = vst->p + offset;
	double *e = vst->e + offset;
	int k;

	if (tmpst == NULL)
	{
	    (*Params(state[0])->_alloc_state)(&tmpst,Params(state[0])->sizest);
	    zero_scalar(tmpst,Params(state[0])->sizest);
	    set_type_of_state(tmpst,EGAS_STATE);
	}
	for (k = 0; k < vsize; ++k)
	{
	    Dens(tmpst) = rho[k];
	    Energy(tmpst) = e[k];
	    Set_params(tmpst,state[k]);
	    p[k] = pressure(tmpst);
	}
}		/*end GENERIC_single_eos_load_pressure*/

/*
*		GENERIC_single_eos_load_sound_speed2():
*
*	Loads a vector of sound speeds into the	appropriate fields of the
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	GENERIC_single_eos_load_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	static Locstate tmpst = NULL;
	Locstate        *state = vst->state + offset;
	double           *rho = vst->rho + offset;
	double           *c2 = vst->c2 + offset;
	double           *e = vst->e + offset;
	double           *FD = NULL;
	int             k;

	if (vst->FD != NULL)
	    FD = vst->FD + offset;

	if (tmpst == NULL)
	{
	    (*Params(state[0])->_alloc_state)(&tmpst,Params(state[0])->sizest);
	    zero_scalar(tmpst,Params(state[0])->sizest);
	    set_type_of_state(tmpst,EGAS_STATE);
	}
	for (k = 0; k < vsize; ++k)
	{
	    Dens(tmpst) = rho[k];
	    Energy(tmpst) = e[k];
	    Set_params(tmpst,state[k]);
	    c2[k] = sound_speed_squared(tmpst);
	    if (FD != NULL)
	        FD[k] = fundamental_derivative(tmpst);
	}
}		/*end GENERIC_single_eos_load_sound_speed2*/

/***************END VECTORIZED THERMODYNAMIC FUNCTIONS *********************/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/

/***************Purely Thermodynamic Hugoniot Functions*********************/

/*
*			GENERIC_dens_Hugoniot():
*
*	Given the state state0 on one side of an oblique shock and the pressure
*	p1 on the other side, this function returns the density rho1 of the
*	state with pressure p1.  Rho1 is found by solving the Hugoniot relation
*
*		(p1 + p0)*(1/rho0 - 1/rho1) = 2*(e1 - e0)
*
*	where e0 and e1 are the specific internal energies of the two
*	respective states.  For a given equation of state the specific
*	internal energy can be expressed as a function of the
*	pressure and density.  Thus the above equation can be solved to
*	give rho1 as a function of state0 and p1.
*
*
*	Reference: Courant and Friedrichs page 302 ff.
*/


LOCAL	double	GENERIC_dens_Hugoniot(
	double    p1,
	Locstate state0)
{
	static const int MAX_NUM_ITER_GENERIC = 20;/*TOLERANCE*/
	static const int MAX_NUM_ITER_BISECTION  = 40;/*TOLERANCE*/
	static const int UNSTABLE		  =  2;/*TOLERANCE*/

	int n, unstable = 0;
	double pbar, p0;
	double e0, emin, emax, emax0;
	double rho, rho0;
	double fmin, fmax, fmax0, f, dfde;
	double x;
	double p, GAM;
	double i2, i02;
	static Locstate ans = NULL;

	if (ans == NULL)
	{
	    (*Params(state0)->_alloc_state)(&ans,Params(state0)->sizest);
	    set_type_of_state(ans,EGAS_STATE);
	    zero_state_velocity(ans,MAXD);
	}
	Set_params(ans,state0);
	p0 = pressure(state0);
	e0 = specific_internal_energy(state0);
	i02 = acoustic_impedance_squared(state0);
	GAM = gruneisen_gamma(state0);
	rho0 = Dens(state0);
	pbar = 0.5*(p1 + p0);
	/*
	 * This initial guess at the behind state energy is obtained by
	 * fitting an SPOLY EOS to the general EOS at the ahead state and
	 * solving the SPOLY Hugoniot eqn.  The relevant formulae are:
	 * 
	 * a) in the sentence containing eqn (7.31) in Menikoff and Plohr,
	 * Rev. Mod. Phys., 1989.
	 *
	 * b) c_0^2 = gam*p0/rho0.
	 *
	 * c) Eqn. (6) in B. Plohr, "Shockless acceleration of thin plates ...",
	 * AIAA, 1988.
	 *
	 * NOTE: Roughly a page of algebra is also required.
	 */
	Energy(ans) = e0 + pbar*(p1-p0)/(i02 + 0.5*(GAM+2.0)*rho0*(p1-p0));

	emin = e0;
	emax = emax0 = e0 + pbar/rho0;

	fmin = p1 - p0;
	fmax = fmax0 = -HUGE_VAL;
	if ( Energy(ans) < emin || Energy(ans) > emax)
	    Energy(ans) = 0.5*(emin+emax);

	/* Solve Hugoniot by Newton iteration */
	for ( n = 0; n < MAX_NUM_ITER_GENERIC; ++n)
	{
	    /* Solve hug. for density */
	    Dens(ans) = rho = pbar/(emax0 - Energy(ans));
	    p = pressure(ans);
	    i2 = acoustic_impedance_squared(state0);
	    GAM = gruneisen_gamma(ans);
	    f = p1 - p;
	    if (fabs(f) < EPSILON*p1)
	    {
	        /*
	         * For a vanishingly weak wave it may occur 
	         * that the convergence criterion is satisfied 
	         * but Dens(ans) < rho0.
	         * In this case return rho0.
	         */
	        return (Dens(ans) > rho0) ? Dens(ans) : rho0;
	    }
	        /* reset energy bounds */
	    if( f > 0.0 )
	    {
	        fmin = f;
	        emin = Energy(ans);
	    }
	    else
	    {
	        fmax = f;
	        emax = Energy(ans);
	    }
	        /*
	         * Derivative for Newton interation.
	         * Remember to use formula above relating density to energy.
	         * Eqs. 2.46, 2.47 from Menikoff & Plohr are also required.
	         */
	    dfde = (GAM*p*rho - i2)/pbar - rho*GAM;
	    Energy(ans) -= f/dfde;

	    /* exit after UNSTABLE failures to satisfy bounds */ 
	    if( Energy(ans) < emin)
	    {
	        Energy(ans) = emin + 0.1*(emax-emin);
	        if (unstable++ > UNSTABLE)
	            break;
	    }
	    else if (Energy(ans) > emax)
	    {
	        Energy(ans) = emax - 0.1*(emax-emin);
	        if (unstable++ > UNSTABLE)
	            break;
	    }
	}
	    
	    /* Bisection method for robust alternative */

	for (n = 0; n < MAX_NUM_ITER_BISECTION; ++n)
	{
	    x = (fmax == fmax0) ? 0.5 : fmin/(fmin-fmax);
	    if (x < 0.1)
	        x = 0.2;
	    else if(x > 0.9)
	        x = 0.8;

	    Energy(ans) = emin + x*(emax - emin);
	    Dens(ans) = pbar/(emax0 - Energy(ans));
	    f = p1 - pressure(ans);
	    if (fabs(f) < EPSILON*p1)
	    {
	        return Dens(ans);
	    }
	    if( f > 0.0 )
	    {
	        fmin = f;
	        emin = Energy(ans);
	    }
	    else
	    {
	        fmax = f;
	        emax = Energy(ans);
	    }
	}
	(void) printf("WARNING in GENERIC_dens_Hugoniot(), "
	              "No convergence to solution\n");
	(void) printf("Dens(ans) = %g, p1 = %g\n",Dens(ans),p1);
	verbose_print_state("state0",state0);
	return Dens(ans);
}	    /*end GENERIC_dens_Hugoniot*/


/*
*			GENERIC_state_w_pr_on_Hugoniot():
*
*	Given the state state0 on one side of an oblique shock and the pressure
*	p1 on the other side, this function returns the thermodynamic variables
*	of the state with pressure p1 (density and internal energy for a
*	GAS_STATE, pressure and density for a TGAS_STATE).  Rho1 is found by
*	solving the Hugoniot relation
*
*		e1 - e0 - 0.5*(p1 + p0)*(V0 - V1) = 0
*
*	where e0 and e1 are the specific internal energies of the two
*	respective states.  For a given equation of state the specific
*	internal energy can be expressed as a function of the
*	pressure and density.  Thus the above equation can be solved to
*	give rho1 and e1 as a function of state0 and p1.  The internal
*	energy is then given by E1 = r1 * e1.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in state1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*	Reference: Courant and Friedrichs page 302 ff.
*/

LOCAL	void	GENERIC_state_w_pr_on_Hugoniot(
	Locstate state0,
	double p1,
	Locstate state1,
	int stype1)
{
	Dens(state1) = dens_Hugoniot(p1,state0);
	Press(state1) = p1;
	zero_state_velocity(state1,MAXD);
	Set_params(state1,state0);
	set_type_of_state(state1,TGAS_STATE);
	set_state(state1,stype1,state1);
}		/*end GENERIC_state_w_pr_on_Hugoniot*/

/*
*			GENERIC_state_w_mf_sqr_on_Hugoniot():
*
*	Given the state state0 on one side of an oblique shock and the square
*	of the mass flux across the shock, this function returns the
*	thermodynamic variables of the state on the opposite side of the shock.
*
*	By definition the square of the mass flux across a shock is given by
*
*			mf_sqr = (p1 - p0) / (V0 - V1)
*
*	where pi and Vi denote the pressure and specific volume on
*	of the two states on either side of the shock.
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

/*ARGSUSED*/
LOCAL	boolean	GENERIC_state_w_mf_sqr_on_Hugoniot(
	Locstate state0,
	double m2,
	Locstate ans,
	int stype1)
{
	static const int MAX_NUM_ITER_GENERIC = 20;/*TOLERANCE*/
	static const int MAX_NUM_ITER_BISECTION  = 40;/*TOLERANCE*/
	static const int UNSTABLE		  =  2;/*TOLERANCE*/

	int n, unstable = 0;
	double b;
	double p0, p1;
	double e0, emin, emax;
	double rho0, V, V0, V1;
	double fmin, fmax, fmax0, f, dfde;
	double x;
	double GAM, i2;

	p0 = pressure(state0);
	emin = e0 = specific_internal_energy(state0);
	rho0 = Dens(state0);
	V0 = 1.0/rho0;
	i2 = acoustic_impedance_squared(state0);
	GAM = gruneisen_gamma(state0);
	set_type_of_state(ans,EGAS_STATE);
	Dens(ans) = (GAM+2.0)*m2/(GAM*m2 + 2.0*i2);
	V1 = 1.0/Dens(ans);
	p1 = p0 + 2.0*(m2 - i2)/(rho0*(GAM+1.0));
	Energy(ans) = e0 + 0.5*(p1 + p0)*(V0 - V1);
	Set_params(ans,state0);
	set_type_of_state(ans,GAS_STATE);
	zero_state_velocity(ans,MAXD);

	emax = e0 + p0*V0 + 0.5*m2*V0*V0;

	fmin = HUGE_VAL;
	fmax = fmax0 = -HUGE_VAL;

	if ( Energy(ans) < emin || Energy(ans) > emax)
	    Energy(ans) = 0.5*(emin+emax);

	for ( n = 0; n < MAX_NUM_ITER_GENERIC; ++n)
	{
	    b = sqrt(p0*p0 + 2.0*m2*(Energy(ans)-e0));
	    V = V0 - (b - p0)/m2;
	    Dens(ans) = 1.0/V;
	    p1 = pressure(ans);
	    f = b - p1;

	    if (fabs(f) < EPSILON*m2)
	    {
	    	set_state(ans,stype1,ans);
	    	return FUNCTION_SUCCEEDED;
	    }

	    if( f > 0.0 )
	    {
	    	fmin = f;
	    	emin = Energy(ans);
	    }
	    else
	    {
	    	fmax = f;
	    	emax = Energy(ans);
	    }
	    i2 = acoustic_impedance_squared(ans);
	    GAM = gruneisen_gamma(ans);

	    dfde = Dens(ans)*GAM + (i2 - m2 - GAM*p1*Dens(ans))/b;
	    if (fabs(dfde) <= f*EPSILON)
	        break;

	    Energy(ans) -= f/dfde;
	    if( Energy(ans) < emin)
	    {
	    	Energy(ans) = emin + 0.1*(emax-emin);
	    	if (unstable++ > UNSTABLE)
	            break;
	    }
	    else if (Energy(ans) > emax)
	    {
	    	Energy(ans) = emax - 0.1*(emax-emin);
	    	if (unstable++ > UNSTABLE)
	            break;
	    }
	}

	/* Bisection method for robust alternative */

	for (n = 0; n < MAX_NUM_ITER_BISECTION; ++n)
	{
	    x = (fmax == fmax0) ? 0.5 : fmin/(fmin-fmax);
	    if (x < 0.1)
	        x = 0.2;
	    else if(x > 0.9)
	        x = 0.8;

	    Energy(ans) = emin + x*(emax - emin);
	    b = sqrt(p0*p0 + 2.0*m2*(Energy(ans)-e0));
	    V = V0 - (b - p0)/m2;
	    Dens(ans) = 1.0/V;
	    p1 = pressure(ans);
	    f = b - p1;
	    if (fabs(f) < EPSILON*m2)
	    {
	    	set_state(ans,stype1,ans);
	    	return FUNCTION_SUCCEEDED;
	    }

	    if( f > 0.0 )
	    {
	    	fmin = f;
	    	emin = Energy(ans);
	    }
	    else
	    {
	    	fmax = f;
	    	emax = Energy(ans);
	    }
	}
	set_state(ans,stype1,ans);
	(void) printf("WARNING in GENERIC_state_w_mf_sqr_on_Hugoniot(), ");
	(void) printf("No convergence to solution\n");
	return FUNCTION_FAILED;
}		/*end GENERIC_state_w_mf_sqr_on_Hugoniot*/


/***************End Purely Thermodynamic Hugoniot Functions*****************/
/***************Velocity Related Hugoniot Functions*************************/

/*
*			GENERIC_pr_normal_vel_wave_curve():
*
*	Computes the pressure on the forward Riemann wave curve given the
*	velocity difference across the wave.
*
*	If du > 0 returns the solution to the system:
*
*                 2
*		du   = (p - p0)*(V0 - V)
*		de   = 0.5*(p + p0)*(V0 - V)
*
*	if du < 0 returns the solution to the system:
*
*		     /p    dP   |
*		du = \    ----  |
*		     /p0  rho c |S
*/

LOCAL	double	GENERIC_pr_normal_vel_wave_curve(
	double du,	/* normal velocity change across shock = (u1 - u0)*/
	Locstate state0)
{
	double p;
	double p_max, p_min;
	double c0 = sound_speed(state0);

	if (du < 0.0)
	{
	    p_max = pressure(state0);
	    p_min = p_max + acoustic_impedance(state0)*du;
	    while ((p_min > Min_pressure(state0)) &&
	    	(riemann_wave_curve(state0,p_min) > du))
	    {
	    	p_max = p_min;
	    	p_min = 0.5*(p_min + Min_pressure(state0));
	    }
	    if (p_min < Min_pressure(state0))
	    	p_min = Min_pressure(state0);
	}
	else
	{
	    p_min = pressure(state0);
	    p_max = p_min + acoustic_impedance(state0)*du;
	    while (riemann_wave_curve(state0,p_max) < du)
	    {
	    	p_min = p_max;
	    	p_max *= 2.0;
	    }
	}

	if (find_root(pr_normal_vel_wave_curve_aux,state0,
		      du,&p,p_min,p_max,EPSILON*c0,EPSILON*(p_max-p_min)) ==
							FUNCTION_FAILED) 
	{
	    screen("ERROR in GENERIC_pr_normal_vel_wave_curve(), "
	           "can't find root\n");
	    clean_up(ERROR);
	}
	return p;
}		/*end GENERIC_pr_normal_vel_wave_curve*/

LOCAL	boolean	pr_normal_vel_wave_curve_aux(
	double		p,
	double		*du,
	Locstate	state0)
{
	*du = riemann_wave_curve(state0,p);
	return FUNCTION_SUCCEEDED;
}		/*end pr_normal_vel_wave_curve_aux*/


/***************End Velocity Related Hugoniot Functions*********************/
/***************Purely Thermodynamic Adiabatic Wave Curve Functions*********/


/*	
*			GENERIC_dens_rarefaction():
*
*	Given the state state0 and the pressure on the other side  
*	of a simple wave in steady irrotational flow, this function 
* 	returns the density on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equations 
*
*                 dh    |
*	   	------- |    = V,  h(p0) = h0;
*		  dp    | s
*       and     
*
*		d (rho)     |           1 
*            ------------   |    =   -------,  rho(p0) = rho0; 
*               d p         | s        c^2
*
*       where h is specific enthalpy,  p is pressure, V is specific 
*	volume (V = 1/rho), rho is density, s is entropy, and c is 
*	sound speed. The derivatives are taken at constant entropy. 
*	
*	Two equations are necessary, because c^2 depends on p which is
*	computed by the EOS-specific pressure fcn using e (energy) which is
*	obtained from h, which is the solution of the first equation.
*/

/*ARGSUSED*/
LOCAL	boolean	rndr(
	double p,
	double *y,
	double *f,
	int n,
	Locstate st)
{
	if (y[1] == 0.0)
	{
	    (void) printf("WARNING in rndr(): divided by 0"); 
	    return FUNCTION_FAILED; 
	}
	Dens(st) = y[1];
	Energy(st) = y[0] - p/y[1];
	f[0] = 1.0/y[1];
	f[1] = 1.0/sound_speed_squared(st);

	return FUNCTION_SUCCEEDED;
}		/* end rndr*/

LOCAL	double	GENERIC_dens_rarefaction(
	double p1,
	Locstate state0)
{
	static Locstate tmpst = NULL;
	double p0, y0[2], y1[2], H;
	double tolLocal; 

	if (tmpst == NULL)
	{
	    (*Params(state0)->_alloc_state)(&tmpst,Params(state0)->sizest);
	    set_type_of_state(tmpst,EGAS_STATE);
	    zero_state_velocity(tmpst,MAXD);
	}
	Set_params(tmpst,state0);
	p0 = pressure(state0);
	y0[0] = specific_enthalpy(state0);
	y0[1] = Dens(state0);
	H = 0.5*(p1 + p0);

	tolLocal = max(0.0001*EPSILON,EPSILON*p0); 
	if (!runga_kutta(p0,y0,p1,y1,&H,2,rndr,tolLocal,tmpst))
	{
	    screen("ERROR in GENERIC_dens_rarefaction(), "
	           "Runga Kutta integration failed\n");
	    clean_up(ERROR);
	}
	return y1[1];
}		/*end GENERIC_dens_rarefaction*/

/*	
*			GENERIC_pressure_rarefaction():
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

/*ARGSUSED*/
LOCAL	boolean	rnpr(
	double rho,
	double *y,
	double *f,
	int n,
	Locstate st)
{
	Dens(st) = rho;
	Energy(st) = y[0];
	f[0] = pressure(st)/sqr(rho);
	return FUNCTION_SUCCEEDED;
}		/* end rnpr*/

LOCAL	double	GENERIC_pressure_rarefaction(
	double rho1,
	Locstate state0)
{
	static Locstate tmpst = NULL;
	double rho0, e0, e1, H;
	if (tmpst == NULL)
	{
	    (*Params(state0)->_alloc_state)(&tmpst,Params(state0)->sizest);
	    set_type_of_state(tmpst,EGAS_STATE);
	    zero_state_velocity(tmpst,MAXD);
	}
	Set_params(tmpst,state0);
	rho0 = Dens(state0);
	e0 = specific_internal_energy(state0);
	H = 0.5*(rho1 + rho0);
	if (!runga_kutta(rho0,&e0,rho1,&e1,&H,2,rnpr,EPSILON*rho0,tmpst))
	{
	    screen("ERROR in GENERIC_pressure_rarefaction(), "
	           "Runga Kutta integration failed\n");
	    clean_up(ERROR);
	}
	return pressure(tmpst);
}		/*end GENERIC_pressure_rarefaction*/


/*	
*			GENERIC_state_on_adiabat_with_pr():
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

LOCAL	void	GENERIC_state_on_adiabat_with_pr(
	Locstate state0,
	double p1,
	Locstate state1,
	int stype1)
{
	Dens(state1) = dens_rarefaction(p1,state0);
	Press(state1) = p1;
	Set_params(state1,state0);
	set_type_of_state(state1,TGAS_STATE);
	zero_state_velocity(state1,MAXD);
	set_state(state1,stype1,state1);
}		/*end GENERIC_state_on_adiabat_with_pr*/

/*	
*			GENERIC_state_on_adiabat_with_dens():
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

LOCAL	void	GENERIC_state_on_adiabat_with_dens(
	Locstate state0,
	double rho1,
	Locstate state1,
	int stype1)
{
	Press(state1) = pressure_rarefaction(rho1,state0);
	Dens(state1) = rho1;
	Set_params(state1,state0);
	set_type_of_state(state1,TGAS_STATE);
	zero_state_velocity(state1,MAXD);
	set_state(state1,stype1,state1);
}		/*end GENERIC_state_on_adiabat_with_dens*/




/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/
/***************General Wave Curve Functions********************************/

/*
*			GENERIC_mass_flux():
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

LOCAL	double	GENERIC_mass_flux(
	double p,
	Locstate state0)
{
	return sqrt(GENERIC_mass_flux_squared(p,state0));
}		/*end GENERIC_mass_flux*/

/*
*			GENERIC_mass_flux_squared():
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

LOCAL	double	GENERIC_mass_flux_squared(
	double p,
	Locstate state0)
{
	double p0 = pressure(state0);

	if (fabs(p - p0) < p0*EPSILON)
	    return acoustic_impedance_squared(state0);
	if (p < p0)
	{
	    double du = riemann_wave_curve(state0,p);

	    return sqr((p-p0)/du);
	}
	else
	    return (p-p0)/(1.0/Dens(state0)-1.0/dens_Hugoniot(p,state0));
}		/*end GENERIC_mass_flux_squared*/


/***************End General Wave Curve Functions****************************/
/***************Functions for the Evaluation of Riemann Solutions***********/

/*
*			GENERIC_oned_fan_state():
*
*	This is a utility function provided for the evaluation of states
*	in a simple wave.   Given sta, it solves for stm using the
*	equation:
*
*	                     / p_m        |            / p_m        |
*	                    /             |           /             |
*	                    \       dP    |           \       G dP  |
*	    w = c_m - c_a +  \    -----   |         =  \     ------ |
*	                      \   rho c   |             \    rho c  |
*	                      /           |             /           |
*	                     /p_a         | S = S_a    / p_a        | S = S_a
*
*	                                               / c_m        |
*	                                              /             |
*	                                              \        dc   |
*	                                            =  \     ------ |
*	                                                \     mu^2  |
*	                                                /           |
*	                                               / c_a        | S = S_a
*
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

struct FAN_AUX	{
	Locstate sta;
	Locstate stm;
	double	ca, cm;
};

LOCAL	double	GENERIC_oned_fan_state(
	double    w,
	Locstate sta,
	Locstate stb,
	Locstate stm,
	int      stype_m,
	boolean     *vacuum)
{
	double	pa = pressure(sta);
	double	pb = pressure(stb);
	double	pm, ca;
	double   wb;
	struct FAN_AUX	Fprms;

	*vacuum = NO;
	Fprms.sta = sta;
	Fprms.stm = stm;
	Fprms.ca = ca = sound_speed(sta);
	Set_params(stm,sta);

	if (!oned_fan_aux(pb,&wb,(POINTER)&Fprms))
	{
	    screen("ERROR in GENERIC_oned_fan_state(), can't evalute "
		   "state behind rarefaction fan\n");
	    clean_up(ERROR);
	    return Fprms.cm;
	}
	if (((w <= wb) && (wb <= 0.0)) || ((0.0 <= wb) && (wb <= w)))
	{
	    set_state(stm,stype_m,stb);
            Fprms.cm = sound_speed(stm);
	}
	else if (((wb <= 0.0) && (0.0 <= w)) || ((w <= 0.0) && (0.0 <= wb)))
	{
	    set_state(stm,stype_m,sta);
            Fprms.cm = sound_speed(stm);
	}
	else if (find_root(oned_fan_aux,(POINTER)&Fprms,w,&pm,pa,pb,
		      EPSILON*ca,EPSILON*(pb-pa))) 
	{
	    set_state(stm,stype_m,stm);
	}
	else
	{
	    *vacuum = YES;
	    screen("ERROR in GENERIC_oned_fan_state(), can't find root\n");
	    clean_up(ERROR);
	    return Fprms.cm;
	}

	return Fprms.cm;
}		/* end GENERIC_oned_fan_state*/

LOCAL	boolean	oned_fan_aux(
	double	p,
	double	*w,
	POINTER	prms)
{
	struct FAN_AUX	*fprms = (struct FAN_AUX  *)prms;
	Locstate sta = fprms->sta;
	Locstate stm = fprms->stm;
	double    cm;

	*w = int_dp_over_rho_c(p,sta,stm);

	/* The use of gaussian_int() may have non-monotonicity */
	/* need to be checked for the function */
	/* *w = gaussian_int(2,pressure(sta),p,stm); */

	fprms->cm = cm = sound_speed(stm);
	*w += cm - fprms->ca;
	return FUNCTION_SUCCEEDED;
}		/*end oned_fan_aux*/
/***************End Functions for the Evaluation of Riemann Solutions********/



/***************Functions to Compute Riemann Solutions**********************/


/*
*			GENERIC_riemann_wave_curve():
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

LOCAL	double	GENERIC_riemann_wave_curve(
	Locstate state0,
	double pstar)
{
	double p0 = pressure(state0);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (pstar < Min_pressure(state0))
	    pstar = Min_pressure(state0);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	if (pstar > p0)
	{
	    double V0 = 1.0/Dens(state0);
	    double Vstar = 1.0/dens_Hugoniot(pstar,state0);
	    return sqrt((pstar-p0)*(V0 - Vstar));
	}
	else
	    return int_dp_over_rho_c(pstar,state0,NULL);
}		/*end GENERIC_riemann_wave_curve*/


/*
*		GENERIC_set_state_for_find_mid_state():
*
*	Copies the Gas state st into the thermodynamic
*	version Tst, for some EOSs a VGas state is set.
*
*	Technical function added for enhanced performance.
*/

LOCAL	void	GENERIC_set_state_for_find_mid_state(
	Locstate Tst,
	Locstate st)
{
	set_state(Tst,VGAS_STATE,st);
}		/*end GENERIC_set_state_for_find_mid_state*/

/*
*			GENERIC_eps_for_Godunov():
*
*	Returns a tolerance to be used to determine convergence of the
*	of the Riemann solver.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	double	GENERIC_eps_for_Godunov(
	Locstate state,
	double pstar,
	double r_eps)
{
	return r_eps;
}		/*end GENERIC_eps_for_Godunov*/

/*
*			GENERIC_initialize_riemann_solver()
*
*	Computes the epsilons and the initial guess for the pressure
*	in the secant iteration of find_mid_state.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_initialize_riemann_solver(
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
	double	pl = pressure(Tsl),		pr = pressure(Tsr);
	double	ul = vel(0,Tsl),		ur = vel(0,Tsr);
	double	il = acoustic_impedance(Tsl),	ir = acoustic_impedance(Tsr);
	*pstar = (il*pr + ir*pl + ir*il*(ul - ur))/(ir + il);
	*p_min = max(Min_pressure(Tsr),Min_pressure(Tsl));
	*eps_u = *eps_p = eps;
}		/*end GENERIC_initialize_riemann_solver*/


/***************End Functions to Compute Riemann Solutions******************/
/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/

/***************TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*********/

/*
*			GENERIC_steady_state_wave_curve():
*
*	Calculates the steady state wave curve through the state
*	state0 with steady state flow speed q0,  parameterized
*	by pressure.  In general the value returned is given by
*
*                       __                          --                  -- --
*                       |                           |      2             |  |
*                       |     p1/p0 - 1             |    M0              |  |
*                       |                           |                    |  |
*       theta =  arctan |-------------------  * sqrt| ---------    -    1|  |
*                       |         2                 |      2             |  |
*                       |gamma0*M0  - (p1/p0 - 1)   |    Mn              |  |
*			|                           |                    |  |
*                       --                          --                  -- --
*							        for p1 > p0
*
*	where gamma0 = adiabatic_gamma(state0), Mn = m/(rho0*c0) =
*	shock normal Mach number, and m = mass_flux across the shock.
*	and
*
*
*			/ p1
*		       /                 2
*	theta =        \     sqrt(1 - (M) ) *  dp		for p1 < p0
*			\		     ------
*			 \                        2
*			 /		     M * c * rho
*		        / p0
*
*	
*	Returns FUNCTION_SUCCEEDED on success,  FUNCTION_FAILED on failure.
*/

LOCAL	boolean	GENERIC_steady_state_wave_curve(
	double p1,
	double M0sq,
	double *theta,
	Locstate state0)
{
	double		rho0 = Dens(state0);
	double		p0 = pressure(state0);
	double		dp, c0sq;
	double		Mnsq, cotb;
	double		A0, A1;
	double		tan_theta;
	double		gam0;

	c0sq = sound_speed_squared(state0);
	if (p1 >= p0)		/* Shock */
	{
	    dp = (p1 - p0) / p0;
	    gam0 = rho0 * c0sq / p0;
	    tan_theta = dp / (gam0*M0sq - dp);
	    Mnsq = mass_flux_squared(p1,state0) / (sqr(rho0)*c0sq);
	    cotb = M0sq/Mnsq - 1.0;
	    cotb = max(0.0,cotb);
	    cotb = sqrt(cotb);
	    tan_theta *= cotb;
	    *theta = (tan_theta > 0.0) ? atan(tan_theta) : 0.0;
	    return FUNCTION_SUCCEEDED;
	}
	else		/* Rarefaction */
	{
	    double M1sq;
	    static	Locstate st1 = NULL;

	    if (st1 == NULL)
	    	(*Params(state0)->_alloc_state)(&st1,Params(state0)->sizest);

	    state_on_adiabat_with_pr(state0,p1,st1,state_type(state0));
	    M1sq = (M0sq*c0sq + 2.0*(specific_enthalpy(state0) -
				     specific_enthalpy(st1))) /
	           sound_speed_squared(st1);
	    if (prandtl_meyer_speed_rinv(&A0,M0sq,state0,&A1,M1sq,
	    			         st1,theta) == FUNCTION_FAILED)
	    {
	    	(void) printf("WARNING in GENERIC_steady_state_wave_curve(), "
		              "prandtl_meyer_speed_rinv() failed\n");
	        return FUNCTION_FAILED;
	    }
	    *theta = -*theta;
	    return FUNCTION_SUCCEEDED;
	}
}		/*end GENERIC_steady_state_wave_curve*/


/*
*			GENERIC_pressure_at_sonic_point():
*
*	Returns the pressure at the sonic point of the shock polar
*	through the state state0 with steady state Mach number M0.
*/

LOCAL	double	GENERIC_pressure_at_sonic_point(
	double M0sq,
	Locstate state0)
{
	double	psonic, pmin, pmax;
	VGas St0;
	Locstate st0 = (Locstate) &St0;

	set_state(st0,VGAS_STATE,state0);
	Vel(st0)[0] = M0sq*sound_speed_squared(st0) +
		      2.0*specific_enthalpy(st0);

	pmin = pressure(st0);
	pmax = max_behind_shock_pr(M0sq,st0);
	if (find_root(Mach_number_squared_behind_oblique_shock,st0,
		      1.0,&psonic,pmin,pmax,SONIC_TOL,EPSILON*(pmax-pmin)) ==
							FUNCTION_FAILED) 
	{
	    screen("ERROR in GENERIC_pressure_at_sonic_point(), "
	           "root not found\n");
	    clean_up(ERROR);
	}

	return psonic;
}		/*end GENERIC_pressure_at_sonic_point*/

LOCAL	boolean	Mach_number_squared_behind_oblique_shock(
	double		p,
	double		*M2,
	Locstate	st0)
{
	static Locstate st1 = NULL;

	if (st1 == NULL)
	    (*Params(st0)->_alloc_state)(&st1,Params(st0)->sizest);
	state_w_pr_on_Hugoniot(st0,p,st1,EGAS_STATE);
	*M2 = (Vel(st0)[0] - 2.0*specific_enthalpy(st1))/
				sound_speed_squared(st1);
	return FUNCTION_SUCCEEDED;
}

/*
*			GENERIC_pr_at_max_turn_angle():
*
*	Given state0 and the Mach number (squared) of state0 in the frame
*	of a shock, this function calculates the pressure at the point of
*	maximum turning angle on the turning angle pressure shock polar
*	through state0.
*
*	Returns FUNCTION_SUCCEEDED if sucessful,  FUNCTION_FAILED otherwise.
*/

LOCAL	boolean	GENERIC_pr_at_max_turn_angle(
	double *prm,
	double M0sq,	/* Mach number of state0 in the frame of the shock */
	Locstate state0)
{
	double	pmin, pmax;
	VGas St0;
	Locstate st0 = (Locstate) &St0;

	set_state(st0,VGAS_STATE,state0);
	Vel(st0)[0] = M0sq;

	pmin = pressure(st0);
	pmax = max_behind_shock_pr(M0sq,st0);
	if (find_root(mta_aux,st0,1.0,prm,pmin,
		      pmax,SONIC_TOL,EPSILON*(pmax-pmin)) == FUNCTION_FAILED) 
	{
	    (void) printf("WARNING in GENERIC_pr_at_max_turn_angle(), "
		          "no root found\n");
	    return FUNCTION_FAILED;
	}

	return FUNCTION_SUCCEEDED;
}		/*end GENERIC_pr_at_max_turn_angle*/

LOCAL	boolean	mta_aux(
	double		p,
	double		*y,
	Locstate	st0)
{
	double M02, p0, i02, i2, m2, rho0, rho;
	double GAM;
	static Locstate st = NULL;

	if (st == NULL)
	    (*Params(st0)->_alloc_state)(&st,Params(st0)->sizest);
	M02 = Vel(st0)[0];
	state_w_pr_on_Hugoniot(st0,p,st,EGAS_STATE);
	p0 = pressure(st0);	p = pressure(st);
	rho0 = Dens(st0);	rho = Dens(st);
	i02 = acoustic_impedance_squared(st0);
	i2 = acoustic_impedance_squared(st);
	m2 = (p - p0)/(1.0/rho0 - 1.0/rho);
	GAM = gruneisen_gamma(st);
	*y = 0.5*(M02*i02/m2 + rho0/rho)*(i2 - m2)/(i2 - 0.5*GAM*(p - p0)*rho);
	return FUNCTION_SUCCEEDED;
}		/*end mta_aux*/


/*
*		GENERIC_state_in_prandtl_meyer_wave():
*
*	This is a utility function provided for the evaluation of states
*	in a Prandtl-Meyer wave.   Given sta, it solves for stm using the
*	equation:
*
*	                / p_m        |             /A_m                |
*	               /             |            /            2       |
*	               \   cos(a) dP |            \       csc A  dA    |
*	w = A_m - A_a + \  --------  |             \  ------------     |
*	                 \  rho c q  |             /          2    2   |
*	                 /           |B=B_a       /    (1 + mu  cot A) |B=B_a
*	                /p_a         |S=S_a      /A_a                  |S=S_a
*
*	The integrals are evaluted at constant entropy and constant
*	B = 0.5*q*q + i, where i is the specific enthalpy.  Here
*	c is the sound speed, q is the flow speed, sin(A) = c/q is the
*	              2
*	Mach angle, mu  = (G - 1)/G, and G is the fundamental derivative
*	G = 1 + rho c dc/dp|S. Note that mu may be complex for some
*	equations of state.
*
*	The returned state1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	state1 can be used to evaluate the pressure, density, and sound speed
*	of the state inside the rarefaction fan.
*	
*	Input data:
*		w = value of w as defined above
*		sta = state ahead of fan
*		A_a = Positive Mach angle of sta, sin(A_a) = c_a/q_a
*		A_b = Positive Mach angle of stb, sin(A_b) = c_b/q_b
*		stype_m = state type of stm
*
*	Output data:
*		stm = state inside fan
*
*	Returns the Mach angle of stm.
*/


/*ARGSUSED*/
LOCAL	double	GENERIC_state_in_prandtl_meyer_wave(
	double w,
	double A_a,
	Locstate sta,
	double A_b,
	Locstate stb,
	Locstate stm,
	int stype_m)
{
	double pa, pb;
	double pm, M0;
	VGas St[2];
	Locstate st0 = (Locstate) St;
	Locstate stmid = (Locstate) (St+1);

	pa = pressure(sta);
	pb = pressure(stb);
	set_state(st0,VGAS_STATE,sta);
	M0 = 1.0/sin(A_a);
	Vel(st0)[0] = M0*M0;

	if (find_root(pmw_aux,st0,w,&pm,pa,pb,SONIC_TOL,
		      EPSILON*fabs(pb-pa)) == FUNCTION_FAILED) 
	{
	    screen("ERROR in GENERIC_state_in_prandtl_meyer_wave(), "
	           "can't find mid state\n");
	    clean_up(ERROR);
	}
	set_state(stm,stype_m,stmid);
	return Vel(stm)[0];
}		/* end GENERIC_state_in_prandtl_meyer_wave*/

LOCAL	boolean	pmw_aux(
	double		p,
	double		*y,
	Locstate	st0)
{
	Locstate	stm = (Locstate) (((VGas*)st0)+1);
	double		A0, A1;
	double		c0sq, M0sq = Vel(st0)[0];
	double		M1sq;

	c0sq = sound_speed_squared(st0);
	state_on_adiabat_with_pr(st0,p,stm,EGAS_STATE);
	M1sq = (M0sq*c0sq+2.0*(specific_enthalpy(st0)-specific_enthalpy(stm)))/
	       sound_speed_squared(stm);
	if (prandtl_meyer_speed_rinv(&A0,M0sq,st0,&A1,M1sq,stm,y) ==
							FUNCTION_FAILED)
	{
	    (void) printf("WARNING in pmw_aux(), "
	                  "prandtl_meyer_speed_rinv() failed\n");
	    return FUNCTION_FAILED;
	}
	*y = A1 - A0 - *y;
	return FUNCTION_SUCCEEDED;
}		/*end pmw_aux*/

/***************END TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*****/

#if defined(COMBUSTION_CODE)
/***************DETONATION SPECIFIC UTILITY FUNCTIONS*********************/

/*
*			GENERIC_CJ_state():
*
* 	This routine finds the state behind a CJ_detonation.
*	The inputs are the initial state "start"
*	and the side (l_or_r, -1 or 1) we are on.
*/

/*ARGSUSED*/
LOCAL	double	GENERIC_CJ_state(
	Locstate CJ,
	int st_type_CJ,
	Locstate start,
	int l_or_r,
	int avail)
{
	screen("ERROR in GENERIC_CJ_state(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}	/* end GENERIC_CJ_state*/


/*
*	 		GENERIC_progress_state(): 
*
*	Finds the gas state as a function of reaction progress
*	in the steady frame.  
*/
	
/*ARGSUSED*/
LOCAL	void	GENERIC_progress_state(
	double prog,		/* reaction progress */
	Locstate init,		/* init = state behind front */
	Locstate ans,		/* TGas states, init = state behind front */
	double max_vol)		/* maximum allowed volume of reacted state */
{
	screen("ERROR in GENERIC_progress_state(), function unimplemented\n");
	clean_up(ERROR);
}	/* end GENERIC_progress_state*/

/*
*			GENERIC_fprint_combustion_params():
*
*	Prints combustion related parameters.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_fprint_combustion_params(
	FILE *file,
	Gas_param *params)
{
	screen("ERROR in GENERIC_fprint_combustion_params(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}
/***************END DETONATION SPECIFIC UTILITY FUNCTIONS*****************/
#endif /* defined(COMBUSTION_CODE) */


/***************METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**************/

/*
*			GENERIC_neumann_riem_inv_moc():
*
*	Uses the characteristic equations for the k-Riemann invariants
*	to update the states along a Neumann boundary.
*
*	This function integrates the characteristic form of Euler's
*	equations:
*
*	         1     dP       dU             acU
*	       -----  ----  +  ----  =  g  -  -----       (+)
*	       rho c   dl       dl              r
*	                 +        +
*		     
*	         1     dP       dU             acU
*	       -----  ----  -  ----  = -g  -  -----       (-)
*	       rho c   dl       dl              r
*	                 -        -
*
*	               dS
*	              ---- = 0                            (0S)
*	               dl
*	                 0
*		     
*	               dV
*	              ---- = 0                            (0V)
*	               dl
*	                 0
*
*	Here:
*		rho = density
*		P   = pressure
*		S   = specific entropy
*		c   = sound speed
*		U   = component of velocity in the normal direction
*		V   = component of velocity orthogonal to the normal
*		g   = gravitational acceleration
*		a   = geometric coefficient = 0 for rectangular geometry
*		                            = nor[0] for cylindrical symmetry
*		                            = 2 for spherical symmetry
*
*	Basic geometry:
*
*		side = POSITIVE_SIDE  if the flow is to the right of the wall
*		                      (-) family is directed towards the wall
*
*			/|ans      flow region
*			/|   \
*		    wall/|    \
*			/|     \(-) characteristic
*			/|      \
*			/|_______\________________
*		       pt[0]   pt[0]+dn
*		        st0     st1
*
*		side = NEGATIVE_SIDE if the flow is to the left of the wall
*		                     (+) family is directed towards the wall
*
*			        flow region    ans|\
*			                      /   |\
*			                     /    |\
*		         (+) characteristic /     |\
*			                   /      |\
*			  ________________/_______|\
*				     pt[0]+dn   pt[0]
*		                        st1      st0
*		                flow region
*
*	Basic algorithm:
*	In this function we use equations (0S), (0V), and the slip
*	boundary condition at the wall to compute the updated entropy and
*	velocity at the new state.  Note that we allow a possibly nonzero wall
*	velocity u0.  Thus the entropy at the new state equals the entropy
*	at state st0,  the tangential velocity of the new state equals that
*	of st0,  and the normal component of velocity of the new state
*	equals u0.  The pressure of the new state is found by integrating
*	either equation (+) or (-) depending on the characteristic family
*	directed towards the wall.  The input variable side determines
*	the side of the computational region with respect to the wall.
*
*	NOTE:  some applications may include an artificial heat
*	conduction.  This can be implemented in a variety of ways.
*	One method is to allow an entropy change between st0 and ans
*	that is proportional to the quantity (T1 - T0)/T0  where
*	T0 and T1 are the termperatures of states st0 and st1 respectively.
*
*	Input variables:
*		pt	coordinates of the wall
*		st0	state at the wall at start of time step
*		u0	wall normal velocity
*		c0	sound speed at wall
*		st1	state at position pt + nor*dn at start of time step
*		side	side of flow relative to the wall
*		dn	grid spacing in wall normal direction
*		nor	wall normal
*		front   Front structure
*	Output variables:
*		ans	time updated state at the wall
*		
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_neumann_riem_inv_moc(
	double     *pt,
	Locstate  st0,
	double     u0,
	double     c0,
	Locstate  st1,
	SIDE      side,
	Locstate  ans,
	double     dn,
	double     *nor,
	Front     *front)
{
	RECT_GRID *gr = front->rect_grid;
	double     p1, pa, pb;
	double     u1, du, vnor;
	double     heat_cond;
	double     i1, ia;
	double     sgn = (side == POSITIVE_SIDE) ? -1.0 : 1.0;
	double     *v0;
	double     time = front->time;
	double     dt = front->dt;
	double     gans[3], g1[3], g;
	double     ptans[3], pt1[3];
	int       number_of_interations = 0;
	int       MAX_NUM_INTERATIONS = 4;/*TOLERANCE*/
	int       i, dim = gr->dim;
	double       alpha = nor[0]*rotational_symmetry();
	double       ra, r, rmin;
	boolean	    is_rot_symmetry;

	is_rot_symmetry = (gr->Remap.remap == IDENTITY_REMAP) ? NO : YES;
	for (i = 0; i < dim; ++i)
	{
	    pt1[i] = pt[i]+dn*nor[i];
	    ptans[i] = pt[i] + u0*nor[i]*dt;
	}
	eval_gravity(pt1,time,g1);
	eval_gravity(ptans,time+dt,gans);
	for (g = 0.0, i = 0; i < dim; ++i)
	    g += 0.5*(g1[i]+gans[i])*nor[i];

	set_state(ans,EGAS_STATE,st0);
	if ((heat_cond = Params(st0)->avisc.heat_cond) > 0.0)
	{
	    double T0, T1, dS;
	    double rho0 = Dens(st0);
	    double c02 = sound_speed_squared(st0);
	    double p0 = pressure(st0);
	    double GAM = gruneisen_gamma(st0);

	    T0 = temperature(st0);
	    T1 = temperature(st1);
	    dS = heat_cond*(T1/T0 - 1.0);

	    Dens(ans) -= Dens(st0)*T0*GAM*dS/c02;
	    Energy(ans) += (1.0 - GAM*p0/(rho0*c02))*T0;
	}

	p1 = pressure(st1);
	i1 = acoustic_impedance(st1);
	u1 = scalar_product(VelocityVector(st1,NULL),nor,dim);

	state_on_adiabat_with_pr(ans,p1,ans,EGAS_STATE);

	du = sgn*(g*dt + u1 - u0);
	if (is_rot_symmetry)
	{
	    rmin = fabs(pos_radius(0.0,gr));
	    r = pos_radius(pt1[0],gr);
	    ra = pos_radius(ptans[0],gr);
	    if ((fabs(r) > rmin) && (fabs(ra) > rmin))
	    	du -= 0.5*alpha*sound_speed(st1)*u1*dt/r;
	}

	pa = p1;
	do
	{
	    double dp;
	    pb = pa;
	    ia = acoustic_impedance(ans);
	    dp = du;
	    if (is_rot_symmetry)
	    {
	    	if ((fabs(r) > rmin) && (fabs(ra) > rmin))
	            dp -= 0.5*alpha*sound_speed(ans)*u0*dt/ra;
	    }
	    dp *= 2.0*ia*i1/(i1 + ia);
	    pa = p1 + dp;
	    if (pa < Min_pressure(st0))
	      pa = Min_pressure(st0);
	    state_on_adiabat_with_pr(ans,pa,ans,EGAS_STATE);
	} while ((fabs(pb - pa) > EPSILON*p1) &&
		 (++number_of_interations < MAX_NUM_INTERATIONS));

	v0 = VelocityVector(st0,NULL);
	vnor = scalar_product(nor,v0,dim);
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = u0*nor[i] + v0[i] - vnor * nor[i];

	set_state(ans,TGAS_STATE,ans);
}		/*end GENERIC_neumann_riem_inv_moc*/

/*
*		GENERIC_shock_ahead_state_riem_inv_moc():
*
*	Uses the characteristic equations for the k-Riemann invariants
*	to update the state ahead of a shock wave.
*
*	This function integrates the characteristic form of Euler's
*	equations:
*
*	         1     dP       dU             acU
*	       -----  ----  +  ----  =  g  -  -----       (+)
*	       rho c   dl       dl              r
*	                 +        +
*		     
*	         1     dP       dU             acU
*	       -----  ----  -  ----  = -g  -  -----       (-)
*	       rho c   dl       dl              r
*	                 -        -
*
*	               dS
*	              ---- = 0                            (0S)
*	               dl
*	                 0
*		     
*	               dV
*	              ---- = 0                            (0V)
*	               dl
*	                 0
*
*	Here:
*		rho = density
*		P   = pressure
*		S   = specific entropy
*		c   = sound speed
*		U   = component of velocity in the normal direction
*		V   = component of velocity orthogonal to the normal
*		g   = gravitational acceleration
*		a   = geometric coefficient = 0 for rectangular geometry
*		                            = 1 for cylindrical symmetry
*		                            = 2 for spherical symmetry
*
*	Basic geometry:
*
*		side = POSITIVE_SIDE
*		                    /
*		                   /ans(position = pt + W*dt)
*		                  /  + 0 -   
*		                 /    +  0  -   
*		                /      +   0   -   
*		          shock/        +    0    -   
*		         front/          +     0     -   
*		             /            +      0      -   
*			____/st0__________st3____st2____st1____
*		           pt
*
*		side = NEGATIVE_SIDE
*	                                  \
*			                ans\(position = pt + W*dt)
*			             + 0 -  \
*			          +  0  -    \shock
*			       +   0   -      \front
*			    +    0    -        \
*		         +     0     -          \
*                     +      0      -            \
*		____st3____st2____st1__________st0\______
*	                                         pt
*
*		+ = forward characteristic (velocity = U + c)
*		0 = middle characteristic (velocity = U)
*		- = backward characteristic (velocity = U - c)
*
*	Basic algorithm:
*	The entropy and tangential component of velocity of the state
*	ans are found by integrating the (0S) and (0V) characteristic
*	equations.  State st2 is the state at the foot of these
*	characteristics,  so the entropy and tangential velocity of
*	ans is the same as that of st2.  The pressure and normal
*	component of velocity of ans are then obtained by integrating
*	the characteristic equations (+) and (-).
*
*	Input variables:
*		st0	state at foot of shock at start of time step
*		st1	state at foot of (-) or (U-c) charateristic 
*		st2	state at foot of (0) or (U) characteristic
*		st3	state at foot of (+) or (U+c) characteristic
*		pt0	coordinates of shock at start of time step
*		side	ahead side of shock
*		add_source	include gravitational and geometric sources
*		dn	grid spacing in wall normal direction
*	        f1      location of st1 = pt0 + f1*dn*nor
*	        f2      location of st2 = pt0 + f2*dn*nor
*	        f3      location of st3 = pt0 + f3*dn*nor
*		nor	wall normal
*		W	predicted shock front velocity
*		front   Front structure
*	Ouput variables:
*		ans	time updated state ahead of the shock
*
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_shock_ahead_state_riem_inv_moc(
	double     *pt0,
	Locstate  st0,
	Locstate  st1,
	Locstate  st2,
	Locstate  st3,
	Locstate  ans,
	double     dn,
	double     f1,
	double     f2,
	double     f3,
	double     *nor,
	double     *W,
	int       add_source,
	double     dt,
	Front     *front)
{
	RECT_GRID *gr = front->rect_grid;
	double     time = front->time;
	double	  p, plast;
	double	  u;
	double	  A, B, C;
	double	  u1, u3, ubar, du;
	double	  i, i1, i3, ibar, di;
	double	  i13bar;
	double	  p1, p3, pbar, dp;
	double	  *v2, vnor;
	double     ptans[3], pt1[3], pt3[3];
	double	  gdt = 0.0, dgdt = 0.0, g13dt = 0.0, dg13dt = 0.0;
	double     chi13 = 0.0, dchi13 = 0.0;
	double     eta = 0.0, eta13 = 0.0, deta13 = 0.0;
	double	  c1 = 0.0, c3 = 0.0;
	double	  alpha = dt*nor[0]*rotational_symmetry();
	double	  r = 0.0;
	boolean	  is_rot_symmetry;
	boolean   	  add_cyl_source = NO;
	int	  number_of_interations = 0;
	int	  MAX_NUM_INTERATIONS = 4;/*TOLERANCE*/
	int	  j, dim = gr->dim;

	DEBUG_ENTER(GENERIC_shock_ahead_state_riem_inv_moc)
	if (DEBUG)
	{
	    print_general_vector("pt0 = ",pt0,dim,", ");
	    print_general_vector("nor = ",nor,dim,"\n");
	    (void) printf("dn = %g, dt = %g\n",dn,dt);
	    (void) printf("add_source = %d\n",add_source);
	    verbose_print_state("st0",st0);
	    verbose_print_state("st1",st1);
	    verbose_print_state("st2",st2);
	    verbose_print_state("st3",st3);
	}

	is_rot_symmetry = (gr->Remap.remap == IDENTITY_REMAP) ? NO : YES;
	for (j = 0; j < dim; ++j)
	{
	    pt1[j] = pt0[j] + f1*dn*nor[j];
	    pt3[j] = pt0[j] + f3*dn*nor[j];
	    ptans[j] = pt0[j] + W[j]*dt;
	}

	u1 = scalar_product(VelocityVector(st1,NULL),nor,dim);
	p1 = pressure(st1);
	c1 = sound_speed(st1);
	i1 = Dens(st1)*c1;

	u3 = scalar_product(VelocityVector(st3,NULL),nor,dim);
	p3 = pressure(st3);
	c3 = sound_speed(st3);
	i3 = Dens(st3)*c3;

	ubar = 0.5*(u1+u3);
	du = u1 - u3;

	pbar = 0.5*(p1 + p3);
	dp = p1 - p3;

	ibar = i13bar = 0.5*(i1 + i3);
	di = i1 - i3;

	if (add_source)
	{
	    double gvans[3], gv1[3], gv3[3];
	    double gans, g1, g3;
	    double g1_bar, g3_bar;

	    eval_gravity(pt1,time,gv1);
	    eval_gravity(pt3,time,gv3);
	    eval_gravity(ptans,time+dt,gvans);
	    for (gans = 0.0, g1 = 0.0, g3 = 0.0, j = 0; j < dim; ++j)
	    {
	        g1 += gv1[j]*nor[j];
	        g3 += gv3[j]*nor[j];
		gans += gvans[j]*nor[j];
	    }
	    g13dt = 0.5*(g1 + g3)*dt;
	    dg13dt = (g1 - g3)*dt;
	    g1_bar = 0.5*(gans + g1);
	    g3_bar = 0.5*(gans + g3);
	    gdt = 0.5*(g1_bar + g3_bar)*dt;
	    dgdt = (g1_bar - g3_bar)*dt;
	    if (is_rot_symmetry)
	    {
	    	if (alpha != 0.0)
	    	{
		    double chi1, chi3;
		    double eta1, eta3;
	            double r1 = pos_radius(pt1[0],gr);
	            double r3 = pos_radius(pt3[0],gr);
	            double rmin = fabs(pos_radius(0.0,gr));

	            r = pos_radius(ptans[0],gr);
	            if ((fabs(r) > rmin) && (fabs(r1) > rmin) && 
		        (fabs(r3) > rmin))
		    {
		        add_cyl_source = YES;

		        eta1 = 0.5*alpha*c1/r1;
		        chi1 = alpha*c1*u1/r1;

		        eta3 = 0.5*alpha*c3/r3;
		        chi3 = alpha*c3*u3/r3;

		        chi13 = 0.5*(chi1 + chi3);
		        dchi13 = chi1 - chi3;

		        eta13 = 0.5*(eta1 + eta3);
		        deta13 = eta1 - eta3;
		    }
	    	}
	    }
	}

	A = ibar;
	B = 0.25*di*dp + (ibar*ibar - 0.25*di*di)*(chi13 + 0.5*(du + dg13dt));
	C = 0.25*di*du +
	    ibar*(g13dt + 0.5*dchi13) + 0.5*di*(chi13 + 0.5*dg13dt) - 0.5*dp;

	(*Params(st0)->_clear_state)(ans,Params(st0)->sizest);
	p = pbar - B/A;
	u = ubar + C/A;
	state_on_adiabat_with_pr(st2,p,ans,EGAS_STATE);

	do {
	    plast = p;
	    i = acoustic_impedance(ans);
	    ibar = 0.5*(i + i13bar);
	    if (is_rot_symmetry)
	    {
	    	if (add_cyl_source == YES)
	    	{
	            double c = sound_speed(ans);
		    eta = 0.5*alpha*c/r;
	    	}
	    }
	    
	    A = ibar - 0.25*eta*di;
	    B = (ibar*ibar - 0.0625*di*di)*(
		ubar*(eta+eta13) + 0.5*du*(1.0 + 0.5*deta13) + 0.5*dgdt +
		eta13*(0.5*ubar*deta13 + 0.5*du*eta13 + gdt));
	    C = 0.5*du*(0.25*di*(1.0 + 0.5*deta13) + ibar*eta13) +
		ubar*(0.5*ibar*deta13 + 0.25*(eta + eta13)*di) +
		ibar*gdt + 0.125*dgdt*di;

	    p = pbar - B/A;
	    u = ubar + C/A;
	    state_on_adiabat_with_pr(ans,p,ans,EGAS_STATE);
	} while ((fabs(p - plast) > EPSILON*pbar) &&
	    	 (++number_of_interations < MAX_NUM_INTERATIONS));

	v2 = VelocityVector(st2,NULL);
	vnor = scalar_product(nor,v2,dim);
	for (j = 0; j < dim; ++j)
	    Vel(ans)[j] = u*nor[j] + v2[j] - vnor * nor[j];

	set_state(ans,TGAS_STATE,ans);
	if (DEBUG)
	    verbose_print_state("ans",ans);
	DEBUG_LEAVE(GENERIC_shock_ahead_state_riem_inv_moc)
}		/*end GENERIC_shock_ahead_state_riem_inv_moc*/


/*
*			GENERIC_shock_moc_plus_rh():
*
*	Given the updated state st0 on the ahead shock, and the state
*	st1 at the foot of the characteristic through the behind state
*	this function uses the method of characteristics and the
*	Rankine-Hugoniot conditions to solve for the updated state
*	behind the shock.
*
*	This function integrates the characteristic form of Euler's
*	equations:
*
*	         1     dP       dU             acU
*	       -----  ----  +  ----  =  g  -  -----       (+)
*	       rho c   dl       dl              r
*	                 +        +
*		     
*	         1     dP       dU             acU
*	       -----  ----  -  ----  = -g  -  -----       (-)
*	       rho c   dl       dl              r
*	                 -        -
*
*	               dS
*	              ---- = 0                            (0S)
*	               dl
*	                 0
*		     
*	               dV
*	              ---- = 0                            (0V)
*	               dl
*	                 0
*
*	together with the Hugoniot conditions across a shock
*
*		      rho  *(U  - s) = rho  * (U  - s) = m
*                        0    0           1     1
*
*                             2                       2
*		rho  *(U  - s)  + P  = rho1 * (U  - s)  + P1
*		   0    0          0      1     1
*
*	                                1
*		             E  - E  = --- (P + P ) * (V  - V )
*		              1    0    2    1   0      0    1 
*
*				  V  = V
*	                           0    1
*
*			
*
*	Here:
*		rho = density
*		P   = pressure
*		S   = specific entropy
*		c   = sound speed
*		s   = shock normal velocity
*		U   = component of velocity in the normal direction
*		V   = component of velocity orthogonal to the normal
*		g   = gravitational acceleration
*		a   = geometric coefficient = 0 for rectangular geometry
*		                            = 1 for cylindrical symmetry
*		                            = 2 for spherical symmetry
*
*		The subscripts on the Hugoniot equations represent
*		the ahead (0) and behind (1) shock states respectively.
*
*	Basic geometry:
*
*		side = POSITIVE_SIDE
*		                    /
*		                ans/sta(position = pt + W*dt)
*		                + /  + 0 -   
*		             +   /    +  0  -   
*		          +     /      +   0   -   
*		       +       /shock   +    0    -   
*		    +         /front     +     0     -   
*	         +           /            +      0      -   
*	    __stb________stm/______________+_______0_______-___
*		           pt
*
*		side = NEGATIVE_SIDE
*	                                  \
*			                sta\ans(position = pt + W*dt)
*			             + 0 -  \ -
*			          +  0  -    \    -    
*			       +   0   -      \       - 
*			    +    0    -   shock\         -
*		         +     0     -     front\            -
*                     +      0      -            \              -
*		___+_______0_______-______________\stm___________stb__
*	                                          pt
*
*		+ = forward characteristic (velocity = U + c)
*		0 = middle characteristic (velocity = U)
*		- = backward characteristic (velocity = U - c)
*
*	Basic algorithm:
*	The state behind the shock and the shock velocity is determined from
*	the state ahead of the shock and one additional piece of information
*	which is obtained by integrating the incoming behind shock
*	characteristic.  Basically the discretized integration of the
*	behind shock incoming characteristic and the Rankine-Hugoniot 
*	equations across the shock provide a complete set of equations
*	to determine the time updated state behind the shock,  and the
*	time updated shock velocity.  It is common in practice to take
*	the net shock velocity and the average of the wave velocity
*	computed from the Riemann solution at the start of the time step
*	and the valued computed from the above set of equations for the
*	velocity at the end of the time step.
*
*	Input variables:
*		sta	state ahead of time updated shock
*		stm	state at foot of shock at start of time step
*		stb	state behind shock at foot of incoming chacteristic
*		pt	coordinates of shock at start of time step
*		dn	grid spacing in wall normal direction
*		nor	shock normal
*		W	first prediction of shock front velocity
*		front   Front structure
*	Output variables:
*		ans	times updated state behind shock
*		W	updated shock speed
*
*		The answer state ans is returned in Gas format.
*/

/*ARGSUSED*/
LOCAL	boolean	GENERIC_shock_moc_plus_rh(
	double     *pt,
	Locstate  sta,
	Locstate  stm,
	Locstate  stb,
	Locstate  ans,
	double     dn,
	double     *nor,
	double     *W,
	int       w_type,
	Front     *front)
{
	RECT_GRID *gr = front->rect_grid;
	boolean      status;
	double     time = front->time;
	double     dt = front->dt;
	double     g, ga[3], gb[3];
	double	  pa, pm, pb, dp, p, plast, pbar;
	double	  ib, ibar, i, m;
	double	  u, ua, ub, du, ca, Ws;
	double	  sgn;
	double	  num, den;
	double	  gdt;
	double	  *va, vnor;
	double     pta[3], ptb[3];
	int	  j, dim = gr->dim;
	int	  number_of_interations = 0;
	int	  MAX_NUM_INTERATIONS = 4;/*TOLERANCE*/
	double	  cb = 0.0;
	double	  r = 0.0, rb = 0.0, rmin = fabs(pos_radius(0.0,gr));
	double	  alpha = dt*nor[0]*rotational_symmetry();
	boolean	  is_rot_symmetry;

	DEBUG_ENTER(GENERIC_shock_moc_plus_rh)
	if (DEBUG)
	{
	    print_general_vector("pt = ",pt,dim,"\n");
	    print_general_vector("nor = ",nor,dim,"\n");
	    (void) printf("dn = %g, dt = %g\n",dn,dt);
	    print_wave_type("w_type = ",w_type,"\n",current_interface());
	    verbose_print_state("sta",sta);
	    verbose_print_state("stm",stm);
	    verbose_print_state("stb",stb);
	}
	sgn = (is_forward_wave(w_type)) ? 1.0 : -1.0;
	is_rot_symmetry = (gr->Remap.remap == IDENTITY_REMAP) ? NO : YES;

	ua = scalar_product(VelocityVector(sta,NULL),nor,dim);
	if (is_rarefaction_wave(w_type))
	{
	    set_state(ans,GAS_STATE,sta);
	    ca = sound_speed(sta);
	    Ws = ua + sgn*ca;
	}
	else
	{
	    for (j = 0; j < dim; ++j)
	    {
		pta[j] = pt[j] + W[j]*dt;
		ptb[j] = pt[j] - dn*nor[j];
	    }
	    eval_gravity(pta,time + dt,ga);
	    eval_gravity(ptb,time,gb);
	    g = 0.5*(scalar_product(ga,nor,dim) + scalar_product(gb,nor,dim));

	    gdt = g*dt;
	    pa = pressure(sta);
	    pb = pressure(stb);
	    dp = pb - pa;
	    pbar = 0.5*(pb + pa);
	    ib = acoustic_impedance(stb);
	    ub = scalar_product(VelocityVector(stb,NULL),nor,dim);
	    du = ub - ua;
	    set_state(ans,EGAS_STATE,stm);
	    p = pressure(ans);
	    u = scalar_product(VelocityVector(ans,NULL),nor,dim);
	    if (is_rot_symmetry)
	    {
	    	if (alpha > 0.0)
	    	{
	            cb = sound_speed(stb);
	    	}
	    }
	    do
	    {
	    	if (DEBUG)
	    	    (void) printf("p = %g\n",p);
	    	plast = p;
	    	i = acoustic_impedance(ans);
	    	ibar = 0.5*(i + ib);
	    	m = mass_flux(p,sta);
	    	Ws = ua + sgn*m/Dens(sta);
	    	num = dp + sgn*ibar*(du + gdt);
	    	den = 1.0 + ibar/m;
	    	if (is_rot_symmetry)
		{
		    if (alpha > 0.0)
		    {
		    	double c = sound_speed(ans);
		    	r = pos_radius(pta[0],gr);
		    	rb = pos_radius(ptb[0],gr);
		    	if ((r > rmin) && (rb > rmin))
		    	{
		            num -= 0.5*alpha*(i*c*ua/r + ib*cb*ub/rb);
		            den += 0.5*sgn*alpha*i*c/(r*m);
		    	}
		    }
		}
		p = pa + num/den;
		u = ua + sgn*(p - pa)/m;
		state_w_pr_on_Hugoniot(sta,p,ans,EGAS_STATE);
	    } while ((fabs(p - plast) > EPSILON*pbar) &&
			(++number_of_interations < MAX_NUM_INTERATIONS));
	    va = VelocityVector(sta,NULL);
	    vnor = scalar_product(nor,va,dim);
	    for (j = 0; j < dim; ++j)
	    	Vel(ans)[j] = u*nor[j] + va[j] - vnor * nor[j];
	}
	for (j = 0; j < dim; ++j)
	    W[j] = 0.5 * (W[j] + Ws * nor[j]);

	set_state(ans,GAS_STATE,ans);
	pm = pressure(stm);
	status = FUNCTION_SUCCEEDED;
	if (fabs((p - pm)/pm) > 0.2) /*TOLERANCE*/
	{
	    if (debugging("smocprh"))
	    {
	        (void) printf("WARNING in GENERIC_shock_moc_plus_rh(), "
	                      "change in solution too strong\n"
	                      "Incoming pressure = %g, updated pressure = %g, "
			      "relative change = %g\n",pm,p,(p-pm)/pm);
	    }
	    status = FUNCTION_FAILED;
	}
	if (DEBUG)
	    verbose_print_state("ans",ans);
	DEBUG_LEAVE(GENERIC_shock_moc_plus_rh)
	return status;
}		/*end GENERIC_shock_moc_plus_rh*/

/***************END METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**********/

/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			GENERIC_prompt_for_state():
*
*	Prompts for a hydrodynamical state.  The form of
*	the states depends of the Eos. 	The type of the state
*	is returned.
*/

LOCAL	void	GENERIC_prompt_for_state(
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
	    if (i == (dim - 1))
		screen("and ");
	    screen("%s",velmesg[i]);
	}
	screen("%s: ",mesg);
	(void) Scanf("%f %f",&Dens(state),&Press(state));
	for (i = 0; i < dim; ++i)
	    (void) Scanf("%f",&Vel(state)[i]);
	(void) getc(stdin); /*read trailing newline*/

	set_state(state,stype,state);
}		/* end GENERIC_prompt_for_state */

/*
*			GENERIC_prompt_for_thermodynamics():
*
*	Prompts for the thermodynamic variables in a state.  Returns
*	a state with the appropriate thermodynamic state and zero velocity.
*	The return status gives the state type representation of the state.
*/

LOCAL	void	GENERIC_prompt_for_thermodynamics(
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
	screen("Enter the density and pressure%s: ",mesg);
	(void) Scanf("%f %f\n",&Dens(state),&Press(state));
}		/* end GENERIC_prompt_for_thermodynamics */


/*
*			GENERIC_free_EOS_params():
*
*	Frees the storage allocated for an equation of state parameter
*	function.
*/

LOCAL	EOS*	GENERIC_free_EOS_params(
	EOS *params)
{
	free(params);
	return NULL;
}		/*end GENERIC_free_EOS_params*/



/***************Problem Type Specific Initialization Functions**************/

/*
*			GENERIC_RT_RS_f():
*
*	Support function for the computation of a solution to the linearized
*	Rayleigh-Taylor flow.
*
*	NEEDED:  More complete documentation
*/

/*ARGSUSED*/
LOCAL	double	GENERIC_RT_RS_f(
	double		s,
	Locstate	amb_st,
	double		dz,
	double		k_sqr,
	double		g_z)
{
	screen("ERROR in GENERIC_RT_RS_f(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end GENERIC_RT_RS_f*/


/*
*			GENERIC_RT_single_mode_perturbation_state():
*
*	Computes the perturbation term for the solution to the linearized
*	Euler equations in the Rayleigh-Taylor problem.  See the appendix to
*
*	``The Dynamics of Bubble Growth for
*				Rayleigh-Taylor Unstable Interfaces''
*	C. L. Gardner, J. Glimm, O. McBryan, R. Menikoff, D. H. Sharp,
*	and Q. Zhang, Phys. Fluids 31 (3), 447-465 (1988).
*
*	for an explanation of the formulas for the case where both fluids
*	are stiffened polytropic gases.
*
*	Note that ans is only the perturbation to the ambient state.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_RT_single_mode_perturbation_state(
	Locstate	ans,
	double		*coords,
	double		t,
	Locstate	amb_st,
	double		z_intfc,
	double		z_bdry,
	MODE		*mode,
	double		g_z)
{
	screen("ERROR in GENERIC_RT_single_mode_perturbation_state(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}		/*end GENERIC_RT_single_mode_perturbation_state*/

/*
*		KH_single_mode_state():
*
*	Computes the state at location coords and time t for the solution of
*	the linearized Euler equations for a single mode Kelvin-Helmholtz
*	perturbation.
*	
*	See I. G. Currie, Fundamental Mechanics of Fluids, Chapter 6.,
*	or see Lamb's Hydrodynamics for the incompressible analysis;
*	for the compressible analysis, use Crocco's equation in place
*	of Bernoulli's equation and the wave equation in place of
*	Laplace's equation.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_KH_single_mode_state(
	Locstate	ans,
	double		*coords,
	double		t,
	Locstate	amb_st,
	double		stream_velocity,
	double		z_intfc,
	double		z_bdry,
	MODE		*mode)
{
	screen("ERROR in GENERIC_KH_single_mode_state(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}		/*end GENERIC_KH_single_mode_state*/


/*
*		GENERIC_compute_isothermal_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isothermal one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*	    P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, T = T_r.
*
*	The solution is given implicitly by
*
*	    f = f_R + g*dz
*	    T = T_r
*
*	where f = e + P/rho - T*S.
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_compute_isothermal_stratified_state(
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	screen("ERROR in GENERIC_compute_isothermal_stratified_state(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}		/*end GENERIC_compute_isothermal_stratified_state */

/*
*		GENERIC_compute_isentropic_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isentropic one dimensional steady flow.
*
*	The solution is computed by solving the differential equation:
*
*	    P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, S = S_r.
*
*	The solution is given implicitly by
*
*	    h =  h_r + g*gz
*	    S = S_r
*
*	where h = e + P/rho is the specific enthalpy.
*
*/

/*ARGSUSED*/
LOCAL	void	GENERIC_compute_isentropic_stratified_state(
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	double h_ref, h_ans, dh;
	double rho_ref;
	double dP, P, P_ref;
	int   i;
	const double eps = 1.0e3*MACH_EPS;/*TOLERANCE*/
	static const int MAX_ITER = 10;/*TOLERANCE*/

	rho_ref = Dens(ref_st);
	h_ref = specific_enthalpy(ref_st);
	P_ref = pressure(ref_st);
	h_ans = h_ref + gz*dz;

	/* Solve h(P,S) = h_ans using Newton's formula */

	P = P_ref + rho_ref*gz*dz;
	state_on_adiabat_with_pr(ref_st,P,ans,TGAS_STATE);
	for (i = 0; i < MAX_ITER; ++i)
	{
	    dh = h_ans - specific_enthalpy(ans);
	    if (fabs(dh) < eps*h_ref)
		break;

	    dP = Dens(ans)*dh;
	    if (fabs(dP) < eps*P_ref)
		break;
	    P += dP;
	    state_on_adiabat_with_pr(ref_st,P,ans,TGAS_STATE);
	}
	if (i == MAX_ITER)
	{
	    (void) printf("WARNING in "
			  "GENERIC_compute_isentropic_stratified_state(), "
			  "Newton iteration did not converge\n");
	}
	zero_state_velocity(ans,Params(ref_st)->dim);
}	/*end GENERIC_compute_isentropic_stratified_state*/

/*
*		GENERIC_compute_constant_density_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in a constant density one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_R
*/

LOCAL	void	GENERIC_compute_constant_density_stratified_state(
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	double rho, pr;
	Set_params(ans,ref_st);
	rho = Dens(ref_st);
	pr = pressure(ref_st);
	Dens(ans) = rho;
	Press(ans) = pr + rho*gz*dz;
	set_type_of_state(ans,TGAS_STATE);
	zero_state_velocity(ans,Params(ref_st)->dim);
}	/*end GENERIC_compute_constant_density_stratified_state*/

/***************End Problem Type Specific Initialization Functions**********/
/***************END INITIALIZATION UTILITY FUNCTIONS************************/

/***************EQUATION OF STATE DOMAIN FUNCTIONS**************************/

LOCAL	double	GENERIC_Min_energy(
	Locstate	state)
{
	return Params(state)->min_energy;
}	/*end GENERIC_Min_energy*/

LOCAL	double	GENERIC_Min_pressure(
	Locstate	state)
{
	return Params(state)->min_pressure;
}	/*end GENERIC_Min_pressure*/

LOCAL	double	GENERIC_Vacuum_dens(
	Locstate	state)
{
	return Params(state)->vacuum_dens;
}	/*end GENERIC_Vacuum_dens*/

LOCAL	double	GENERIC_Raref_press(
	Locstate	state)
{
	return Params(state)->raref_press;
}	/*end GENERIC_Raref_press*/

#if defined(COMBUSTION_CODE)
LOCAL	double	GENERIC_Tol_alpha(
	Locstate	state)
{
	return Params(state)->tol_alpha;
}	/*end GENERIC_Tol_alpha*/

LOCAL	double	GENERIC_Tol_pressure(
	Locstate	state)
{
	return Params(state)->tol_press;
}	/*end GENERIC_tol_pressure*/
#endif /* defined(COMBUSTION_CODE) */

/***************END EQUATION OF STATE DOMAIN FUNCTIONS*************************/

/*
*			prandtl_meyer_speed_rinv():
*
*	This function computes the purely speed dependent part of
*	the Riemann invariant for a steady irrotational two dimensional
*	supersonic flow.  The full riemann invariants are given by
*
*		theta - rinv for the GAMMA+ family
*
*	and
*
*		theta + rinv for the GAMMA- family.
*
*	This function returns rinv as the value of the integral
*
*	/ p0					/ q1
*      /                   2		       /                   2
*      \     sqrt(1 - (c/q)  ) * dp      =     \     sqrt(1 - (c/q)  ) * dq
*	\		        ----            \                       ----
*	/		        q*c*rho         /			 c
*      / p1				       / q0
*
*	Where the sound speed c and the flow speed q are related by
*	the Bernoulli relation
*
*			0.5*q*q + i = 0.5*sqr(q0) + i0
*
*	where i is the specific internal enthalpy and qhat is the
*	constant critical speed.
*
*	It is assumed that state1 is has correctly ft_assigned thermodynamic
*	data so that state1 and state0 have the same entropy and that
*	q0, state0, q1 and state1 satisfy the Bernoulli relation above.
*
*	In addition this function computes and returns the Mach angles 
*	A0, and A1 for the flow given by
*
*			sin(Ai) = qi/ci, i = 0,1.
*
*	The flow speeds of the two states may be passed either directly
*	or as Mach numbers.  If the flag mach_numbers_given is YES
*	the quantities x0sq and x1sq are assumed to be the steady
*	state Mach numbers of the two states respectively.  Otherwise
*	these numbers are assumed to be the steady state flow speeds.
*/

LOCAL	boolean	prandtl_meyer_speed_rinv(
	double		*A0,
	double		x0sq,
	Locstate	state0,
	double		*A1,
	double		x1sq,
	Locstate	state1,
	double		*rinv)
{
	double		M0sq, M1sq;
	double		cotA0, cotA1;
	double		p0, p1, H, drinv;
	double		zero = 0.0;
	static		Locstate tmpst = NULL;

	M0sq = x0sq;	M1sq = x1sq;
	if (M0sq < SONIC_MINUS_SQR || M1sq < SONIC_MINUS_SQR) 
	{
	    (void) printf("WARNING in prandtl_meyer_speed_rinv(), "
	                  "Subsonic state in simple wave\n"
	                  "M0sq = %g, M1sq = %g\n",M0sq,M1sq);
	    return FUNCTION_FAILED;
	}
	M0sq = max(1.0,M0sq);		M1sq = max(1.0,M1sq);
	cotA0 = sqrt(M0sq - 1.0);	cotA1 = sqrt(M1sq - 1.0);
	*A0 = atan2(1.0,cotA0);		*A1 = atan2(1.0,cotA1);

	if (tmpst == NULL)
	    (*Params(state0)->_alloc_state)(&tmpst,sizeof(VGas));

	set_state(tmpst,VGAS_STATE,state0);
	zero_state_velocity(tmpst,Params(tmpst)->dim);
	Vel(tmpst)[0] = sin(*A0);
	p0 = pressure(state0);
	p1 = pressure(state1);
	H = 0.5*(p1 + p0);
	if (!runga_kutta(p1,&zero,p0,&drinv,&H,1,pmsr,EPSILON,tmpst))
	{
	    screen("ERROR in prandtl_meyer_speed_rinv(), "
	           "Runga Kutta integration failed\n");
	    clean_up(ERROR);
	    return FUNCTION_FAILED;
	}
	*rinv = drinv;
	return FUNCTION_SUCCEEDED;
}		/*end prandtl_meyer_speed_rinv*/


/*ARGSUSED*/
LOCAL	boolean	pmsr(
	double p,
	double *y,
	double *f,
	int n,
	Locstate state0)
{
	double sinA, sinA0;
	double c, c0;
	double H, H0;
	static Locstate st = NULL;

	if (st == NULL)
	    (*Params(state0)->_alloc_state)(&st,Params(state0)->sizest);

	sinA0 = Vel(state0)[0];
	c0 = sound_speed(state0);
	H0 = specific_enthalpy(state0);
	state_on_adiabat_with_pr(state0,p,st,state_type(state0));
	c = sound_speed(st);
	H = specific_enthalpy(st);
	sinA = c*sinA0/(sqrt(c0*c0 + 2.0*(H0-H)*sqr(sinA0)));
	*y = sinA*sqrt(1.0 - sqr(sinA))/(Dens(st)*c*c);
	return FUNCTION_SUCCEEDED;
}		/*end pmsr*/


/*
*			int_dp_over_rho_c():
*
*
*	Returns the value of the integral:
*
*		 
*                        / p1        |
*                       /            |
*                       \      dP    |
*                        \   ------  |
*                         \   rho c  |
*                         /          |
*                        / p0        | S
*
*	Also returns the state on the adiabat containing state0 with
*	pressure P.
*
*/

LOCAL	double	int_dp_over_rho_c(
	double		p1,
	Locstate	state0,
	Locstate	state1)
{
	double	ans;
	double	p, dp, p0 = pressure(state0);
	double	V;
	double	sgn;
	double	eps = 1.0e-2;/*TOLERANCE*/
	int	i, n;

	DEBUG_ENTER(int_dp_over_rho_c)
	
	if (p1 == -HUGE_VAL) /*Reduce to finite interval*/
	    p1 = pressure_rarefaction(0.001*Dens(state0),state0);/*TOLERANCE*/

	if (p0 != 0.0) 
	    dp = (p1 - p0)/p0;
	else if (p1 != 0.0) 
	    dp = (p1 - p0)/p1;
	else 
	    return 0.0;   
	sgn = (dp < 0.0) ? -1.0 : 1.0;
	dp *= sgn;
	n = irint(sqrt(dp*dp*dp/eps));
	n = max(1,n);
	if (n > 100)
	{
	    /* This maximum number of integral segments is set at 100 
	     * after numerical tests. Please change, if accuracy is  
	     * found to be deteriorated due to insufficient number of 
	     * segments
	     */ 

	    if (DEBUG)
	        (void) printf("WARNING in generic int_dp_over_rho_c(), "
			      "number of intervals reset from %d to %d\n",
			      n,100);
	    n = 100; 
	}
	dp = (p1 - p0)/n;
	sgn = (dp < 0.0) ? -1.0 : 1.0;

	if (DEBUG)
	{
	    (void) printf("p1 = %g, n = %d, dp = %g\n",p1,n,dp);
	    verbose_print_state("state0",state0);
	}
	if (state1 == NULL)
	{
	    static Locstate tmpst = NULL;

	    if (tmpst == NULL)
	    	(*Params(state0)->_alloc_state)(&tmpst,Params(state0)->sizest);
	    state1 = tmpst;
	}

	set_state(state1,EGAS_STATE,state0);
	for (p = p0, ans = 0.0, i = 0; i < n; ++i)
	{
	    p += dp;   
	    V = 1.0/Dens(state1);
	    state_on_adiabat_with_pr(state1,p,state1,EGAS_STATE);
	    if(DEBUG && ( dp*(V - 1.0/Dens(state1)) < 0 ) ) 
	    {
	        (void) printf("\nWARNING in int_dp_over_rho_c(), "
			      "dp*dV = %g", dp*(V - 1.0/Dens(state1)) ); 
	    }
	    ans += sqrt(fabs(dp*(V - 1.0/Dens(state1))));  
		/* the sign was recovered at the return */ 
	}
	ans *= sgn; 

	DEBUG_LEAVE(int_dp_over_rho_c)
	return ans;
	
}		/*end int_dp_over_rho_c*/

LOCAL   double           gaussian_int(
	int             np,
	double           a,
	double           b,
	Locstate	state)
{
        double  pt[6], wt[6], ans=0.0;
	int    i,j,nn;
	double  a0,al,dt;

	if (a <= 0 || b <= 0)
	{
	    screen("ERROR in gaussian_int(), "
	           "nonpositive numbers : a = %g, b = %g\n",a,b);
	    clean_up(ERROR);
	}

	a0 = log(a);
	al = log(b);

	if(debugging("fast_gaussian"))
	    nn = 1;
	else
	{
	    dt = fabs(al - a0)/10;
	    for(nn = 1; nn < dt; nn++);
	}
	dt=(al - a0)/nn;

	for(j = 0; j < nn; j++)
	{
	    (void) legendre_init(np,pt,wt);
	    (void) linear_transform(np,pt,wt,a0+j*dt,al-(nn-1-j)*dt);
	    (void) exp_transform(np,pt,wt);
	    for(i = 0; i < np; i++)
		ans += wt[i]*one_over_rho_c(pt[i],state);
	}

	state_on_adiabat_with_pr(state,b,state,EGAS_STATE); /* Set to b */

	return ans;
}		/*end gaussian_int*/

LOCAL   void    legendre_init(
	int     np,
	double*  point,
	double*  weight)
{
        int i;
        switch(np)
	{
	case 1:
	       point[0] = 0;
	       weight[0] =2;
	       break;
	case 2:
	       point[1] = .577350269189626;
	       weight[1] = 1;
	       break;
	case 3:
	       point[1] = 0;
	       point[2] = .774596669241483;
	       weight[1] =.8888888888888889;
	       weight[2] =.5555555555555556;
	       break;
	case 4:
	       point[2] = .339981043584856;
	       point[3] = .861136311594053;
	       weight[2] = .652145154862546;
	       weight[3] = .347854845137454;
	       break;
	case 5:
	       point[2] = 0;
	       point[3] = .538469310105683;
	       point[4] = .906179845938664;
	       weight[2] = .5688888888888889;
	       weight[3] = .478628670499366;
	       weight[4] = .236926885056189;
	       break;
	case 6:
	       point[3] = .238619186083197;
	       point[4] = .661209386466265;
	       point[5] = .932469514203152;
	       weight[3] = .467913934572691;
	       weight[4] = .360761573048139;
	       weight[5] = .171324492379170;
	       break;
	default:
	       screen("ERROR in legendre_init(), "
	              "given number %d is greater than 6", np);
	       clean_up(ERROR);
	       break;
	}
	for (i = 0; i < np/2; i++)
	{
	       point[i] = -point[np-i-1];
	       weight[i] = weight[np-i-1];
	}
}		/*end legendre_init*/

LOCAL   void    linear_transform(
	int     np,
	double*  point,
	double*  weight,
	double   a,
	double   b)
{
        int i;
	for (i = 0; i < np; i++)
	{
	    point[i] = ((b-a)*point[i] + (b+a))/2;
	    weight[i] = weight[i]*(b-a)/2;
	}
}		/*end linear_transform*/       

LOCAL   void    exp_transform(
	int     np,
	double*  point,
	double*  weight)
{
        int i;

	for (i = 0; i < np; i++)
	{
	    point[i] = exp(point[i]);
	    weight[i] = point[i]*weight[i];
	}
}		/*end exp_transform*/       

LOCAL   double   one_over_rho_c(
	double		p,
	Locstate	state)
{
        double rho, c, rhoc;

	state_on_adiabat_with_pr(state,p,state,EGAS_STATE);
	rho = Dens(state);
	c = sound_speed(state);
	rhoc = rho*c;
	if (rhoc <= 0.0)
	{
	    screen("ERROR in one_over_rho_c(), negative values : "
		   "rho = %g, c = %g, rho*c = %g\n",rho,c,rhoc);
	    clean_up(ERROR);
#if defined(USE_OVERTURE)
            return -HUGE_VAL;
#else
            return NAN;
#endif /* if defined(USE_OVERTURE) */
	}
	return 1/rhoc;
}		/*end one_over_rho_c*/
