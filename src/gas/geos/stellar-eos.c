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
*				stellar-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <geos/stellar.h>

	/* LOCAL Function Prototypes */
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL	double	STELLAR_internal_energy(Locstate);
LOCAL	double	STELLAR_pressure(Locstate);
LOCAL	double	STELLAR_sound_speed_squared(Locstate);
LOCAL	double	STELLAR_acoustic_impedance_squared(Locstate);
LOCAL	double	STELLAR_specific_internal_energy(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL	double	STELLAR_specific_enthalpy(Locstate);
LOCAL	double	STELLAR_temperature(Locstate);
LOCAL	double	STELLAR_entropy(Locstate);
LOCAL	double	STELLAR_adiabatic_gamma(Locstate);
LOCAL	double	STELLAR_gruneisen_gamma(Locstate);
LOCAL	double	STELLAR_fundamental_derivative(Locstate);
LOCAL	double	STELLAR_C_V(Locstate);
LOCAL	double	STELLAR_C_P(Locstate);
LOCAL	double	STELLAR_K_S(Locstate);
LOCAL	double	STELLAR_K_T(Locstate);

	/* MATERIAL PROPERTY FUNCTIONS */
LOCAL	double	STELLAR_bulk_viscosity(Locstate);
LOCAL	double	STELLAR_shear_viscosity(Locstate);
LOCAL   double   STELLAR_density(Locstate);
LOCAL 	double   STELLAR_heat_coeff(Locstate);

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
LOCAL	void	STELLAR_single_eos_load_pressure_and_sound_speed2(Vec_Gas*,
                                                               int,int);
LOCAL	void	STELLAR_single_eos_load_pressure_and_gammas(Vec_Gas*,int,int);
LOCAL	void	STELLAR_single_eos_load_pressure(Vec_Gas*,int,int);
LOCAL	void	STELLAR_single_eos_load_sound_speed2(Vec_Gas*,int,int);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
LOCAL	double	STELLAR_dens_Hugoniot(double,Locstate);
LOCAL	void	STELLAR_state_w_pr_on_Hugoniot(Locstate,double,Locstate,int);
LOCAL	boolean	STELLAR_state_w_mf_sqr_on_Hugoniot(Locstate,double,Locstate,int);

	/* Velocity Related Hugoniot Functions */
LOCAL	double	STELLAR_pr_normal_vel_wave_curve(double,Locstate);

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL	double	STELLAR_dens_rarefaction(double,Locstate);
LOCAL	double	STELLAR_pressure_rarefaction(double,Locstate);
LOCAL	void	STELLAR_state_on_adiabat_with_pr(Locstate,double,Locstate,int);
LOCAL	void	STELLAR_state_on_adiabat_with_dens(Locstate,double,Locstate,int);

	/* General Wave Curve Functions */
LOCAL	double	STELLAR_mass_flux(double,Locstate);
LOCAL	double	STELLAR_mass_flux_squared(double,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
LOCAL	double	STELLAR_oned_fan_state(double,Locstate,Locstate,Locstate,int,
				    boolean*);

	/* Functions to Compute Riemann Solutions */
LOCAL	double	STELLAR_riemann_wave_curve(Locstate,double);
LOCAL	void	STELLAR_set_state_for_find_mid_state(Locstate,Locstate);
LOCAL	double	STELLAR_eps_for_Godunov(Locstate,double,double);
LOCAL	void	STELLAR_initialize_riemann_solver(Locstate,Locstate,double*,double*,
	                                       double,double*,double*,
	                                       boolean(*)(Locstate,Locstate,
						       double,double*,double*,
						       double*,double*,double*,
						       double*,
						       RIEMANN_SOLVER_WAVE_TYPE*,
						       RIEMANN_SOLVER_WAVE_TYPE*));


	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
LOCAL	boolean	STELLAR_steady_state_wave_curve(double,double,double*,Locstate);
LOCAL	double	STELLAR_pressure_at_sonic_point(double,Locstate);
LOCAL	boolean	STELLAR_pr_at_max_turn_angle(double*,double,Locstate);
LOCAL	double	STELLAR_state_in_prandtl_meyer_wave(double,double,Locstate,double,
	                                         Locstate,Locstate,int);

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
LOCAL	double	STELLAR_CJ_state(Locstate,int,Locstate,int,int);
LOCAL	void	STELLAR_progress_state(double,Locstate,Locstate,double);
LOCAL	void	STELLAR_fprint_combustion_params(FILE*,Gas_param*);
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
LOCAL	void	STELLAR_neumann_riem_inv_moc(double*,Locstate,double,double,Locstate,
	                                  SIDE,Locstate,double,double*,Front*);
LOCAL	void	STELLAR_shock_ahead_state_riem_inv_moc(double*,Locstate,Locstate,
	                                            Locstate,Locstate,
	                                            Locstate,double,double,
	                                            double,double,double*,double*,
						    int,double,Front*);
LOCAL	boolean	STELLAR_shock_moc_plus_rh(double*,Locstate,Locstate,Locstate,
	                               Locstate,double,double*,double*,int,
	                               Front*);

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	STELLAR_prompt_for_state(Locstate,int,Gas_param*,const char*);
LOCAL	void	STELLAR_prompt_for_thermodynamics(Locstate,Gas_param*,
					       const char*);
LOCAL	void	STELLAR_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	STELLAR_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                           Gas_param*);
LOCAL	EOS*	STELLAR_free_EOS_params(EOS*);
LOCAL	void	STELLAR_prompt_for_EOS_params(INIT_DATA*,Gas_param*,
					   const char*,const char*);

	/* Problem Type Specific Initialization Functions */
LOCAL	double	STELLAR_RT_RS_f(double,Locstate,double,double,double);
LOCAL	void	STELLAR_RT_single_mode_perturbation_state(Locstate,double*,double,
	                                               Locstate,double,double,
	                                               MODE*,double);
LOCAL	void	STELLAR_KH_single_mode_state(Locstate,double*,double,Locstate,
	                                  double,double,double,MODE*);
LOCAL	void	STELLAR_compute_isothermal_stratified_state(Locstate,double,double,
	                                                 Locstate);
LOCAL	void	STELLAR_compute_isentropic_stratified_state(Locstate,double,double,
	                                                 Locstate);

LOCAL	void	set_eos_function_hooks(EOS*);

typedef struct {
	double ca, cb, dS, alpha, beta, gam, c2, c4, psi, eta;
	double ua, rad, delta;
	boolean cyl_coord;
} STELLAR_S_MPRH;

	/* Polytropic specific utility functions */
LOCAL	boolean	stellar_fS(double,double*,POINTER);
LOCAL	void	set_local_gamma(Locstate);
LOCAL   char    gaspath[256];

EXPORT	EOS	*set_STELLAR_eos(
	EOS	*eos,
	char	*path)
{
	if (eos == NULL)
	    scalar(&eos,sizeof(STELLAR_EOS));
	strcpy(gaspath,path);
	(void) set_GENERIC_eos(eos);
	set_eos_function_hooks(eos);
	return eos;
}	/*end set_STELLAR_eos*/

LOCAL	void	set_eos_function_hooks(
	EOS *eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_internal_energy = STELLAR_internal_energy;
	eos->_pressure = STELLAR_pressure;
	eos->_sound_speed_squared = STELLAR_sound_speed_squared;
	eos->_acoustic_impedance_squared = STELLAR_acoustic_impedance_squared;
	eos->_specific_internal_energy = STELLAR_specific_internal_energy;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_specific_enthalpy = STELLAR_specific_enthalpy;
	eos->_temperature = STELLAR_temperature;
	eos->_entropy = STELLAR_entropy;
	eos->_adiabatic_gamma = STELLAR_adiabatic_gamma;
	eos->_gruneisen_gamma = STELLAR_gruneisen_gamma;
	eos->_fundamental_derivative = STELLAR_fundamental_derivative;
	eos->_C_V = STELLAR_C_V;
	eos->_C_P = STELLAR_C_P;
	eos->_K_S = STELLAR_K_S;
	eos->_K_T = STELLAR_K_T;

	/* MATERIAL PROPERTY FUNCTIONS */
	eos->_bulk_viscosity = STELLAR_bulk_viscosity;
	eos->_shear_viscosity = STELLAR_shear_viscosity;
	eos->_density = STELLAR_density;
	eos->_heat_coeff = STELLAR_heat_coeff;

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
	eos->_single_eos_load_pressure_and_sound_speed2 =
	    STELLAR_single_eos_load_pressure_and_sound_speed2;
	eos->_single_eos_load_pressure_and_gammas =
	    STELLAR_single_eos_load_pressure_and_gammas;
	eos->_single_eos_load_pressure = STELLAR_single_eos_load_pressure;
	eos->_single_eos_load_sound_speed2 = STELLAR_single_eos_load_sound_speed2;

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
	eos->_dens_Hugoniot = STELLAR_dens_Hugoniot;
	eos->_state_w_pr_on_Hugoniot = STELLAR_state_w_pr_on_Hugoniot;
	eos->_state_w_mf_sqr_on_Hugoniot = STELLAR_state_w_mf_sqr_on_Hugoniot;

	/* Velocity Related Hugoniot Functions */
	eos->_pr_normal_vel_wave_curve = STELLAR_pr_normal_vel_wave_curve;

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	eos->_dens_rarefaction = STELLAR_dens_rarefaction;
	eos->_pressure_rarefaction = STELLAR_pressure_rarefaction;
	eos->_state_on_adiabat_with_pr = STELLAR_state_on_adiabat_with_pr;
	eos->_state_on_adiabat_with_dens = STELLAR_state_on_adiabat_with_dens;

	/* General Wave Curve Functions */
	eos->_mass_flux = STELLAR_mass_flux;
	eos->_mass_flux_squared = STELLAR_mass_flux_squared;

	/* Functions for the Evaluation of Riemann Solutions */
	eos->_oned_fan_state = STELLAR_oned_fan_state;

	/* Functions to Compute Riemann Solutions */
	eos->_riemann_wave_curve = STELLAR_riemann_wave_curve;
	eos->_set_state_for_find_mid_state = STELLAR_set_state_for_find_mid_state;
	eos->_eps_for_Godunov = STELLAR_eps_for_Godunov;
	eos->_initialize_riemann_solver = STELLAR_initialize_riemann_solver;

	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
	eos->_steady_state_wave_curve = STELLAR_steady_state_wave_curve;
	eos->_pressure_at_sonic_point = STELLAR_pressure_at_sonic_point;
	eos->_pr_at_max_turn_angle = STELLAR_pr_at_max_turn_angle;
	eos->_state_in_prandtl_meyer_wave = STELLAR_state_in_prandtl_meyer_wave;

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
	eos->_CJ_state = STELLAR_CJ_state;
	eos->_progress_state = STELLAR_progress_state;
	eos->_fprint_combustion_params = STELLAR_fprint_combustion_params;
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
	eos->_neumann_riem_inv_moc = STELLAR_neumann_riem_inv_moc;
	eos->_shock_ahead_state_riem_inv_moc =
	    STELLAR_shock_ahead_state_riem_inv_moc;
	eos->_shock_moc_plus_rh = STELLAR_shock_moc_plus_rh;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_prompt_for_state = STELLAR_prompt_for_state;
	eos->_prompt_for_thermodynamics = STELLAR_prompt_for_thermodynamics;
	eos->_fprint_EOS_params = STELLAR_fprint_EOS_params;
	eos->_read_print_EOS_params = STELLAR_read_print_EOS_params;
	eos->_free_EOS_params = STELLAR_free_EOS_params;
	eos->_prompt_for_EOS_params = STELLAR_prompt_for_EOS_params;

	/* Problem Type Specific Initialization Functions */
	eos->_RT_RS_f = STELLAR_RT_RS_f;
	eos->_RT_single_mode_perturbation_state =
	    STELLAR_RT_single_mode_perturbation_state;
	eos->_KH_single_mode_state = STELLAR_KH_single_mode_state;
	eos->_compute_isothermal_stratified_state = STELLAR_compute_isothermal_stratified_state;
	eos->_compute_isentropic_stratified_state = STELLAR_compute_isentropic_stratified_state;
}


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*			STELLAR_internal_energy():
*
*	Returns the internal energy per unit volume of a state.
*/

LOCAL	double	STELLAR_internal_energy(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	switch (state_type(state)) 
	{
	case GAS_STATE:
	    return Energy(state) - kinetic_energy(state);

	case EGAS_STATE:
	    return Energy(state)*Dens(state);

	case VGAS_STATE:
	    return Dens(state) * Int_en(state);

	case TGAS_STATE:
	    return Press(state) / GAMMA(state);

	case FGAS_STATE:
	    return Dens(state)*Cv(state)*Temperature(state);

	default:
	    screen("ERROR: in STELLAR_internal_energy(), no such state type\n");
	    clean_up(ERROR);
	}
	return ERROR_FLOAT;
}		/*end STELLAR_internal_energy*/


/*
*			STELLAR_pressure():
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

LOCAL	double	STELLAR_pressure(
	Locstate state)
{
	double pr, rho;

	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	if (is_obstacle_state(state))
	    return HUGE_VAL;
	rho = Dens(state);
	switch (state_type(state)) 
	{

	case	GAS_STATE:
	    pr = GAMMA(state) * (Energy(state) - kinetic_energy(state));
	    break;

	case	EGAS_STATE:
	    pr = GAMMA(state) * Energy(state) * rho;
	    break;

	case	FGAS_STATE:
	    pr = R(state)*Temperature(state)*Dens(state);
	    break;

	case	TGAS_STATE:
	case	VGAS_STATE:
	    pr = Press(state);
	    break;

	default:
	    screen("ERROR in STELLAR_pressure(), no such state type\n");
	    clean_up(ERROR);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(state));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end STELLAR_pressure*/


/*
*			STELLAR_sound_speed_squared():
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

LOCAL	double	STELLAR_sound_speed_squared(
	Locstate state)
{
	double c2;
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	if (state_type(state) == VGAS_STATE)
	{
	    double c = Sound_speed(state);
	    return c*c;
	}
	c2 = Local_gamma(state)*pressure(state)/Dens(state);
	return c2;
}		/*end STELLAR_sound_speed_squared*/


/*
*		STELLAR_acoustic_impedance_squared():
*
*	Returns the square of the local acoustic impedance of the state.
*
*                        2     dP  |
*			i = - ---- |
*                              dV  |S
*/

LOCAL	double	STELLAR_acoustic_impedance_squared(
	Locstate state)
{
	double i2;
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	if (state_type(state) == VGAS_STATE)
	{
	    double i = Dens(state)*Sound_speed(state);
	    return i*i;
	}
	i2 =  Local_gamma(state)*pressure(state)*Dens(state);
	return i2;
}		/*end STELLAR_acoustic_impedance_squared*/

/*
*			STELLAR_specific_internal_energy():
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/

LOCAL	double	STELLAR_specific_internal_energy(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	switch (state_type(state)) {

	case	GAS_STATE:
	    return (Energy(state) - kinetic_energy(state))/Dens(state);

	case	EGAS_STATE:
	    return Energy(state);

	case	TGAS_STATE:
	    return Coef6(state)* Press(state) / Dens(state);

	case	FGAS_STATE:
	    return Cv(state)*Temperature(state);

	case	VGAS_STATE:
	    return Int_en(state);

	default:
	    screen("ERROR in STELLAR_specific_internal_energy(), "
	           "no such state type\n");
	    clean_up(ERROR);
	    break;
	}
	return ERROR_FLOAT;

}		/*end STELLAR_specific_internal_energy*/


/***************END PRIMARY THERMODYNAMIC FUNCTIONS ************************/
/***************SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS ***********/

/*
*			STELLAR_specific_enthalpy():
*
*	This function computes the specific enthalpy of the given state.
*
*			H = E + P*V
*
*	E = specific internal energy, P = pressure, V = specific volume.
*
*/

LOCAL	double	STELLAR_specific_enthalpy(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
#if defined(VERBOSE_PLUS_GAS)
	if (state_type(state) == VGAS_STATE) return Enthalpy(state);
#endif /* defined(VERBOSE_PLUS_GAS) */

	return Coef7(state) * pressure(state) / Dens(state);
}		/*end STELLAR_specific_enthalpy*/


/*
*			STELLAR_temperature():
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	double	STELLAR_temperature(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	if (state_type(state) == FGAS_STATE)
	    return Temperature(state);

	return (pressure(state)/Dens(state))/R(state);
}		/*end STELLAR_temperature*/

/*
*			STELLAR_entropy():
*
*	Returns the specific entropy of a state.
*/

LOCAL	double	STELLAR_entropy(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	if (state_type(state) == VGAS_STATE)
	    return Entropy(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(state) < Vacuum_dens(state))
	    return 0.0;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	return (log(pressure(state)) - Local_gamma(state)*log(Dens(state))) *
	       Coef6(state)*R(state);
}		/*end STELLAR_entropy*/

/*
*			STELLAR_adiabatic_gamma():
*
*	Returns the dimensionless sound speed
*
*		gamma = - d(log P)/d(log V) | .
*					     S
*	As usual P = thermodynamic pressure,  V = specific volume
*	and S = specific entropy.
*/

LOCAL	double	STELLAR_adiabatic_gamma(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	return Local_gamma(state);
}		/*end STELLAR_adiabatic_gamma*/


/*
*			STELLAR_gruneisen_gamma():
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

LOCAL	double	STELLAR_gruneisen_gamma(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	return GAMMA(state);
}		/*end STELLAR_gruneisen_gamma*/

/*
*			STELLAR_fundamental_derivative():
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

LOCAL	double	STELLAR_fundamental_derivative(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	return 0.5*(Local_gamma(state) + 1.0);
}		/*end STELLAR_fundamental_derivative*/

/*
*			STELLAR_C_V():
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*/

LOCAL	double	STELLAR_C_V(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	return Cv(state);
}	/* end STELLAR_C_V */

/*
*			STELLAR_C_P():
*
*	Specific heat at constant pressure.
*
*
*                        dS  |
*		C_P = T ---- |
*                        dT  | P
*/

LOCAL	double	STELLAR_C_P(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	return Local_gamma(state)*Cv(state);
}	/* end STELLAR_C_P */

/*
*			K_S():
*
*	Isentropic compressibility.
*
*                        1   dV  |
*		K_S = - --- ---- |
*                        V   dP  | S
*/

LOCAL	double	STELLAR_K_S(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	return 1.0/(Local_gamma(state)*pressure(state));
}	/* end STELLAR_K_S */

/*
*			STELLAR_K_T():
*
*	Isothermal compressibility.
*
*                        1   dV  |
*		K_T = - --- ---- |
*                        V   dP  | T
*/

LOCAL	double	STELLAR_K_T(
	Locstate state)
{
	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	return 1.0/pressure(state);
}	/* end STELLAR_K_T */



/***************END SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS *******/

/************** MATERIAL PROPERTY FUNCTIONS ********************************/

LOCAL	double	STELLAR_bulk_viscosity(
	Locstate state)
{
	return 0.0;
}	/*end STELLAR_bulk_viscosity */

LOCAL	double	STELLAR_shear_viscosity(
	Locstate state)
{
	return 0.0;
}	/*end STELLAR_shear_viscosity */

/*
*			STELLAR_density():
*
*	Returns the density of a state using other variables.
*
*/

LOCAL	double	STELLAR_density(
	Locstate state)
{
	double pr, rho;

	if (Local_gamma_set(state) == NO)
	    set_local_gamma(state);
	if (is_obstacle_state(state))
	    return HUGE_VAL;
	pr = Press(state);
	switch (state_type(state)) 
	{

	case	EGAS_STATE:
	    rho = pr / GAMMA(state) / Energy(state);
	    /*	from    pr = GAMMA(state) * Energy(state) * rho; */
	    break;

	case	FGAS_STATE:
	    rho = pr / R(state) / Temperature(state);
	    /* from   pr = R(state)*Temperature(state)*Dens(state); */
	    break;

	case	GAS_STATE:
	case	TGAS_STATE:
	case	VGAS_STATE:
	    rho = Dens(state);
	    break;

	default:
	    screen("ERROR in STELLAR_density(), no such state type\n");
	    clean_up(ERROR);
	}
	return rho;
}		/*end STELLAR_density*/

/*
*               STELLAR_heat_coeff():
*/
LOCAL   double   STELLAR_heat_coeff(
        Locstate state)
{
	return 0.0;

}               /*end STELLAR_heat_coeff*/

/************** END MATERIAL PROPERTY FUNCTIONS ****************************/

/***************VECTORIZED THERMODYNAMIC FUNCTIONS *************************/

/*
*		STELLAR_single_eos_load_pressure_and_sound_speed2():
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

LOCAL	void	STELLAR_single_eos_load_pressure_and_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	double *rho = vst->rho + offset;
	double *p = vst->p + offset, *c2 = vst->c2 + offset;
	Locstate *state = vst->state + offset;
	double gm;
	int   k;

	if (Vec_Gas_field_set(vst,re))
	{
	    double *re = vst->re + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        p[k] = (gm-1.0)*re[k];
	    }
	}
	else
	{
	    double *e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        p[k] = (gm-1.0)*rho[k]*e[k];
	    }
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	for (k = 0; k < vsize; ++k)
	{
	    if (Local_gamma_set(state[k]) == NO)
	        set_local_gamma(state[k]);
	    gm = Local_gamma(state[k]);
	    c2[k] = gm*p[k]/rho[k];
	}
	if (vst->FD != NULL)
	{
	    double *FD = vst->FD + offset;
	    for (k = 0; k < vsize; ++k)
	    {
	    	if (Local_gamma_set(state[k]) == NO)
	            set_local_gamma(state[k]);
	        gm = Local_gamma(state[k]);
	        FD[k] = 0.5*(gm+1.0);
	    }
	}
}		/*end STELLAR_single_eos_load_pressure_and_sound_speed2*/


/*
*		STELLAR_single_eos_load_pressure_and_gammas():
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

LOCAL	void	STELLAR_single_eos_load_pressure_and_gammas(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	double *rho = vst->rho + offset;
	double *p = vst->p + offset;
	double *c2 = vst->c2 + offset, *GAM = vst->GAM + offset;
	Locstate *state = vst->state + offset;
	double gm;
	int   k;

	if (Vec_Gas_field_set(vst,re))
	{
	    double *re = vst->re + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        GAM[k] = gm - 1.0;
	        p[k] = (gm-1.0)*re[k];
	    }
	}
	else
	{
	    double *e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        GAM[k] = gm - 1.0;
	        p[k] = (gm-1.0)*rho[k]*e[k];
	    }
	}
	if (vst->FD != NULL)
	{
	    double *FD = vst->FD + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        FD[k] = 0.5*(gm+1.0);
	    }
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	for (k = 0; k < vsize; ++k)
	{
	    if (Local_gamma_set(state[k]) == NO)
	        set_local_gamma(state[k]);
	    gm = Local_gamma(state[k]);
	    c2[k] = gm*p[k]/rho[k];
	}
}		/*end STELLAR_single_eos_load_pressure_and_gammas*/

/*
*			STELLAR_single_eos_load_pressure():
*
*	Loads a vector of pressures into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	STELLAR_single_eos_load_pressure(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	double *p = vst->p + offset;
	double *rho = vst->rho + offset;
	Locstate *state = vst->state + offset;
	double gm;
	int   k;

	if (Vec_Gas_field_set(vst,re))
	{
	    double *re = vst->re + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        p[k] = (gm-1.0)*re[k];
	    }
	}
	else
	{
	    double *e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        p[k] = (gm-1.0)*rho[k]*e[k];
	    }
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
}		/*end STELLAR_single_eos_load_pressure*/

/*
*		STELLAR_single_eos_load_sound_speed2():
*
*	Loads a vector of sound speeds into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	STELLAR_single_eos_load_sound_speed2(
	Vec_Gas *vst,
	int     offset,
	int     vsize)
{
	double *e = vst->e + offset;
	double *c2 = vst->c2 + offset;
	Locstate *state = vst->state + offset;
	double gm;
	int   k;

	for (k = 0; k < vsize; ++k)
	{
	    if (Local_gamma_set(state[k]) == NO)
	    	set_local_gamma(state[k]);
	    gm = Local_gamma(state[k]);
	    c2[k] = gm*(gm-1.0)*e[k];
	}
	if (vst->FD != NULL)
	{
	    double *FD = vst->FD + offset;
	    for (k = 0; k < vsize; ++k)
	    {
		if (Local_gamma_set(state[k]) == NO)
	    	    set_local_gamma(state[k]);
		gm = Local_gamma(state[k]);
	        FD[k] = 0.5*(gm+1.0);
	    }
	}
}		/*end STELLAR_single_eos_load_sound_speed2*/


/***************END VECTORIZED THERMODYNAMIC FUNCTIONS *********************/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/

/***************Purely Thermodynamic Hugoniot Functions*********************/

/*
*			STELLAR_dens_Hugoniot():
*
*	Given the state st0 on one side of an oblique shock and the pressure
*	p1 on the other side, this function returns the density rho1 of the
*	state with pressure p1.  Rho1 is found by solving the Hugoniot relation
*
*		(p1 + p0)*(1/rho0 - 1/rho1) = 2*(e1 - e0)
*
*	where e0 and e1 are the specific internal energies of the two
*	respective states.  For a given equation of state the specific
*	internal energy can be expressed as a function of the
*	pressure and density.  Thus the above equation can be solved to
*	give rho1 as a function of st0 and p1.
*
*
*	Reference: Courant and Friedrichs page 302 ff.
*/


LOCAL	double	STELLAR_dens_Hugoniot(
	double p1,
	Locstate st0)
{
	double	p0, c4;

	p0 = pressure(st0);
	c4 = Coef4(st0);
	return Dens(st0)*(p1 + p0*c4)/ (p0 + p1*c4);
}		/*end STELLAR_dens_Hugoniot*/


/*
*			STELLAR_state_w_pr_on_Hugoniot():
*
*	Given the state st0 on one side of an oblique shock and the pressure
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
*	give rho1 and e1 as a function of st0 and p1.  The internal
*	energy is then given by E1 = r1 * e1.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in st1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*	Reference: Courant and Friedrichs page 302 ff.
*/

LOCAL	void	STELLAR_state_w_pr_on_Hugoniot(
	Locstate st0,
	double p1,
	Locstate st1,
	int stype1)
{
	double	p0;
	double	rho0;
	double	c4;

	p0 = pressure(st0);
	rho0 = Dens(st0);
	c4 = Coef4(st0);
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,TGAS_STATE);
	Press(st1) = p1;
	Dens(st1) = rho0*(p1 + p0*c4) / (p0 + p1*c4);
	set_local_gamma(st1);/*because gamma is used in the following, we have to set it now. */
	set_type_of_state(st1,stype1);
	switch(stype1)
	{
	case TGAS_STATE:
	    Press(st1) = p1;
	    break;
	case GAS_STATE:
	    Energy(st1) = p1*Coef6(st1);
	    reset_gamma(st1);
	    break;
	case EGAS_STATE:
	    Energy(st1) = p1*Coef6(st1)/Dens(st1);
	    reset_gamma(st1);
	    break;
	case FGAS_STATE:
	    Temperature(st1) = p1/(Dens(st1)*R(st1));
	    reset_gamma(st1);
	    break;
	case VGAS_STATE:
	    Press(st1) = p1;
	    set_type_of_state(st1,TGAS_STATE);
	    set_state(st1,VGAS_STATE,st1);
	    break;
	default:
	    screen("ERROR in state_w_pr_on_Hugoniot(), "
	           "Unknown state type %d\n",stype1);
	    clean_up(ERROR);
	}
}		/*end STELLAR_state_w_pr_on_Hugoniot*/

/*
*			STELLAR_state_w_mf_sqr_on_Hugoniot():
*
*	Given the state st0 on one side of an oblique shock and the square
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
*		If stype1 == GAS_STATE the energy in st1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*/

LOCAL	boolean	STELLAR_state_w_mf_sqr_on_Hugoniot(
	Locstate st0,
	double m2,
	Locstate st1,
	int stype1)
{
	double	p0, p1;
	double 	r0;
	double	c4;

	p0 = pressure(st0);
	r0 = Dens(st0);
	c4 = Coef4(st0);
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,TGAS_STATE);
	p1 = Coef5(st0)*m2/r0 - c4*p0;
	Press(st1) = p1;
	Dens(st1) = Dens(st0)*(p1 + p0*c4)/(p0 + p1*c4);
	set_local_gamma(st1);/*because gamma is used in the following, we have to set it now */
	set_type_of_state(st1,stype1);
	switch(stype1)
	{
	case TGAS_STATE:
	    Press(st1) = p1;
	    break;
	case GAS_STATE:
	    Energy(st1) = p1 * Coef6(st1);
	    reset_gamma(st1);
	    break;
	case EGAS_STATE:
	    Energy(st1) = p1 * Coef6(st1)/Dens(st1);
	    reset_gamma(st1);
	    break;
	case FGAS_STATE:
	    Temperature(st1) = p1/(Dens(st1)*R(st1));
	    reset_gamma(st1);
	    break;
	case VGAS_STATE:
	    Press(st1) = p1;
	    set_type_of_state(st1,TGAS_STATE);
	    set_state(st1,VGAS_STATE,st1);
	    break;
	default:
	    screen("ERROR in state_w_mf_sqr_on_Hugoniot(), "
	           "Unknown state type %d\n",stype1);
	    clean_up(ERROR);
	}
	return FUNCTION_SUCCEEDED;
}		/*end STELLAR_state_w_mf_sqr_on_Hugoniot*/


/***************End Purely Thermodynamic Hugoniot Functions*****************/
/***************Velocity Related Hugoniot Functions*************************/

/*
*			STELLAR_pr_normal_vel_wave_curve():
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

LOCAL	double	STELLAR_pr_normal_vel_wave_curve(
	double du,	/* normal velocity change across shock = (u1 - u0)*/
	Locstate st0)
{
	double b, c, disc, p0;

	p0 = pressure(st0);

	if (du > 0.0)
	{
	    c = Local_gamma(st0)*Dens(st0)*sqr(du) / p0;
	    b = c / (1.0 + Coef4(st0));
	    disc = (b*b + 4.0*c);
	    return p0 * (1.0 + 0.5*(b + sqrt(disc)));
	}
	else if (du < 0.0)
	{
	    double y, p1;

	    y = 1.0 + Coef2(st0)*du/sound_speed(st0);
	    if (y <= 0.0)
	    	return Min_pressure(st0);
	    p1 = pow(y,1.0/Coef3(st0));
	    if (p1 < Min_pressure(st0))
	    	p1 = Min_pressure(st0);
	    return p1;
	}
	else
	    return p0;
}		/*end STELLAR_pr_normal_vel_wave_curve*/


/***************End Velocity Related Hugoniot Functions*********************/
/***************Purely Thermodynamic Adiabatic Wave Curve Functions*********/

/*	
*			STELLAR_dens_rarefaction():
*
*	Given the state st0 and the pressure on the other side of
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

LOCAL	double	STELLAR_dens_rarefaction(
	double p1,
	Locstate st0)
{
	if (Local_gamma_set(st0) == NO)
	    set_local_gamma(st0);
	return Dens(st0) * pow(p1/pressure(st0),1.0/Local_gamma(st0));
}		/*end STELLAR_dens_rarefaction*/

/*	
*			STELLAR_pressure_rarefaction():
*
*	Given the state st0 and the density on the other side of
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

LOCAL	double	STELLAR_pressure_rarefaction(
	double rho1,
	Locstate st0)
{
	if (Local_gamma_set(st0) == NO)
	    set_local_gamma(st0);
	return pressure(st0) * pow(rho1/Dens(st0),Local_gamma(st0));
}		/*end STELLAR_pressure_rarefaction*/


/*	
*			STELLAR_state_on_adiabat_with_pr():
*
*	Given the state st0 and the pressure on the other side of
*	a simple wave in steady irrotational flow, this function returns
*	the thermodynamic variable on the other side.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in st1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*/

LOCAL	void	STELLAR_state_on_adiabat_with_pr(
	Locstate st0,
	double p1,
	Locstate st1,
	int stype1)
{
	if (Local_gamma_set(st0) == NO)
	    set_local_gamma(st0);
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,TGAS_STATE);
	Press(st1) = p1;
	Dens(st1) = Dens(st0)*pow(p1/pressure(st0),1.0/Local_gamma(st0));
	set_local_gamma(st1);/*because gamma is used in the following, we have to set it now. */
	set_type_of_state(st1,stype1);
	switch(stype1)
	{
	case TGAS_STATE:
	    Press(st1) = p1;
	    break;
	case GAS_STATE:
	    Energy(st1) = p1 * Coef6(st1);
	    reset_gamma(st1);
	    break;
	case EGAS_STATE:
	    Energy(st1) = p1*Coef6(st1)/Dens(st1);
	    reset_gamma(st1);
	    break;
	case FGAS_STATE:
	    Temperature(st1) = p1/(Dens(st1)*R(st1));
	    reset_gamma(st1);
	    break;
	case VGAS_STATE:
	    Press(st1) = p1;
	    set_type_of_state(st1,TGAS_STATE);
	    set_state(st1,VGAS_STATE,st1);
	    break;
	default:
	    screen("ERROR in state_on_adiabat_with_pr(), "
	           "Unknown state type %d\n",stype1);
	    clean_up(ERROR);
	}
}		/*end STELLAR_state_on_adiabat_with_pr*/

/*	
*			STELLAR_state_on_adiabat_with_dens():
*
*	Given the state st0 and the density on the other side of
*	a simple wave in steady irrotational flow, this	function returns
*	the pressure and internal energy on the other side.
*
*	IMPORTANT NOTES:
*		1.  If stype1 == GAS_STATE the energy in st1 is
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

LOCAL	void	STELLAR_state_on_adiabat_with_dens(
	Locstate st0,
	double rho1,
	Locstate st1,
	int stype1)
{
	double p1; 	/* pressure of the answer state */

	if (Local_gamma_set(st0) == NO)
	    set_local_gamma(st0);
	Set_params(st1,st0);
	zero_state_velocity(st1,Params(st0)->dim);
	set_type_of_state(st1,TGAS_STATE);
	p1 =  pressure(st0) * pow(rho1/Dens(st0),Local_gamma(st0));
	Press(st1) = p1;
	Dens(st1) = rho1;
	set_local_gamma(st1);/*because gamma is used in the following, we have to set it now. */
	set_type_of_state(st1,stype1);
	switch(stype1)
	{
	case TGAS_STATE:
	    Press(st1) = p1;
	    break;
	case GAS_STATE:
	    Energy(st1) = p1 * Coef6(st1);
	    reset_gamma(st1);
	    break;
	case EGAS_STATE:
	    Energy(st1) = p1*Coef6(st1)/rho1;
	    reset_gamma(st1);
	    break;
	case FGAS_STATE:
	    Temperature(st1) = p1/(rho1*R(st1));
	    reset_gamma(st1);
	    break;
	case VGAS_STATE:
	    Press(st1) = p1;
	    set_type_of_state(st1,TGAS_STATE);
	    set_state(st1,VGAS_STATE,st1);
	    break;
	default:
	    screen("ERROR in state_on_adiabat_with_dens(), "
	           "Unknown state type %d\n",stype1);
	    clean_up(ERROR);
	}
}		/*end STELLAR_state_on_adiabat_with_dens*/

/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/
/***************General Wave Curve Functions********************************/

/*
*			STELLAR_mass_flux():
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

LOCAL	double	STELLAR_mass_flux(
	double p,
	Locstate st0)
{
	double p0, rho0;
	double xi, m, i0;

	p0 = pressure(st0);
	rho0 = Dens(st0);
	if (p < p0)
	{
	    i0 = acoustic_impedance(st0);
	    xi = (p > 0.0) ? p/p0 : 0.0;
	    if ( (1.0 - xi) < EPS)
	    	return i0;
	    m = Coef3(st0)*(1.0-xi)/(1.0-pow(xi,Coef3(st0)));
	    return i0*m;
	}
	else
	    return sqrt(rho0*(Coef1(st0)*p + Coef2(st0)*p0));
}		/*end STELLAR_mass_flux*/

/*
*			STELLAR_mass_flux_squared():
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

LOCAL	double	STELLAR_mass_flux_squared(
	double p,
	Locstate st0)
{
	double p0, rho0;
	double xi, m, i02;

	p0 = pressure(st0);
	rho0 = Dens(st0);
	if (p < p0)
	{
	    i02 = acoustic_impedance_squared(st0);
	    xi = p/p0;
	    if ((1.0 - xi) < EPS)
	        return i02;
	    m = Coef3(st0)*(1.0-xi)/(1.0-pow(xi,Coef3(st0)));
	    return i02*m*m;
	}
	else
	    return rho0*(Coef1(st0)*p + Coef2(st0)*p0);
}		/*end STELLAR_mass_flux_squared*/


/***************End General Wave Curve Functions****************************/
/***************Functions for the Evaluation of Riemann Solutions***********/

/*
*				STELLAR_oned_fan_state():
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
*	derivative of gas dynamics.  The returned st1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	st1 can be used to evaluate the pressure, density, and sound speed
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

/*ARGSUSED*/
LOCAL	double	STELLAR_oned_fan_state(
	double    w,
	Locstate sta,
	Locstate stb,
	Locstate stm,
	int      stype_m,
	boolean  *vacuum)
{
	double	c_a, c_m, p_a;
	double	c2, c3, c4;

	zero_state_velocity(stm,Params(sta)->dim);
	*vacuum = NO;

	p_a = pressure(sta);
	set_type_of_state(stm,TGAS_STATE);
	c_a = sound_speed(sta);
	c2 = 1.0 / Coef2(sta);
	c3 = 1.0 / Coef3(sta);
	c4 = Coef4(sta);
	c_m = c_a + c4*w;
	if (c_m <= 0.0)
	{
	    /* rarefaction to vacuum */
	    state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
	    c_m = 0.0;
	    *vacuum = YES;
	}
	else
	{
	    Set_params(stm,sta);
	    Dens(stm) = Dens(sta)*pow(c_m/c_a,c2);
	    Press(stm) = p_a*pow(c_m/c_a,c3);
	    reset_gamma(stm);
	}

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Press(stm) < Min_pressure(sta))
	{
	    state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
	    c_m = 0.0;
	    *vacuum = YES;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	set_state(stm,stype_m,stm);
	return c_m;
}		/* end oned_fan_state*/
/***************End Functions for the Evaluation of Riemann Solutions********/



/***************Functions to Compute Riemann Solutions**********************/


/*
*			STELLAR_riemann_wave_curve():
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

LOCAL	double	STELLAR_riemann_wave_curve(
	Locstate st0,
	double pstar)
{
	double rho0 = Dens(st0), p0 = pressure(st0);
	double c1, c2, c3;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (pstar < Min_pressure(st0))
	    pstar = Min_pressure(st0);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	c1 = Coef1(st0);
	c2 = Coef2(st0);
	c3 = Coef3(st0);

	return (pstar < p0) ?
	        sound_speed(st0)*(pow(pstar/p0,c3) - 1.0)/ c2 :
	        (pstar-p0)/sqrt(rho0*(c1*pstar+c2*p0));
}		/*end STELLAR_riemann_wave_curve*/


/*
*		STELLAR_set_state_for_find_mid_state():
*
*	Copies the Gas state st into the thermodynamic
*	version Tst, for some EOSs a VGas state is set.
*
*	Technical function added for enhanced performance.
*/

LOCAL	void	STELLAR_set_state_for_find_mid_state(
	Locstate Tst,
	Locstate st)
{
	set_state(Tst,TGAS_STATE,st);
}		/*end STELLAR_set_state_for_find_mid_state*/

/*
*			STELLAR_eps_for_Godunov():
*
*	Returns a tolerance to be used to determine convergence of the
*	of the Riemann solver.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	double	STELLAR_eps_for_Godunov(
	Locstate state,
	double pstar,
	double r_eps)
{
	return r_eps;
}		/*end STELLAR_eps_for_Godunov*/

/*
*			STELLAR_initialize_riemann_solver()
*
*	Computes the epsilons and the initial guess for the pressure
*	in the secant iteration of find_mid_state.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	void	STELLAR_initialize_riemann_solver(
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
	double pl, pr;
	double cl, cr, ul_tdl, ur_tdl, z;
	double c2l, c2r, c3l, c3r;
	double vl = vel(0,Tsl), vr = vel(0,Tsr);
	double dutdl;

	*eps_u = *eps_p = eps;
	pl = pressure(Tsl), pr = pressure(Tsr);
	if (Eos(Tsl) != Eos(Tsr))
	{
#if defined(UNRESTRICTED_THERMODYNAMICS)
	    *p_min = -HUGE_VAL;
#else /* defined(UNRESTRICTED_THERMODYNAMICS) */
	    *p_min = max(Min_pressure(Tsl),Min_pressure(Tsr));
#endif /* defined(UNRESTRICTED_THERMODYNAMICS) */

	    /*
	     * Setting pstar = 0.5*(pl + pr) at an SPOLY--POLY contact can
	     * yield pstar < 0 for situations in which cavitation DOES NOT
	     * occur (that is, vacuum DOES NOT form at contact).
	     * In this case, the (negative) initial guess for pstar may be
	     * passed to STELLAR_mass_flux(), from secant_find_mid_state() or
	     * godunov_find_mid_state(), which should never happen.
	     */
	    *pstar = 0.5*(pl + pr);
	    *pstar = max(*pstar,*p_min);
	    return;
	}

	c2l = Coef2(Tsl);
	c3l = Coef3(Tsl);
	c2r = Coef2(Tsr);
	c3r = Coef3(Tsr);
	cl = sound_speed(Tsl);
	cr = sound_speed(Tsr);
	ul_tdl = vl + cl/c2l;
	ur_tdl = vr - cr/c2r;
	dutdl = ul_tdl - ur_tdl;
	if (pl >= pr)
	{
	    z = (c2l*cr/(c2r*cl))*pow(pl/pr,c3l);
	    *pstar = (dutdl > 0.0) ? pl*pow(c2l*dutdl/((1.0 + z)*cl),1.0/c3l) :
				     0.5*min(pl,pr);
	}
	else
	{
	    z = (c2r*cl/(c2l*cr))*pow(pr/pl,c3r);
	    *pstar = (dutdl > 0.0) ? pr*pow(c2r*dutdl/((1.0 + z)*cr),1.0/c3r) :
				     0.5*min(pl,pr);
	}
#if defined(UNRESTRICTED_THERMODYNAMICS)
	*p_min = -HUGE_VAL;
#else /* defined(UNRESTRICTED_THERMODYNAMICS) */
	*pstar = max(*pstar,Min_pressure(Tsl));
	*pstar = max(*pstar,Min_pressure(Tsr));
	*p_min = Min_pressure(Tsl);
#endif /* defined(UNRESTRICTED_THERMODYNAMICS) */
	*pstar = max(*pstar,*p_min);
}		/*end STELLAR_initialize_riemann_solver*/


/***************End Functions to Compute Riemann Solutions******************/
/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/

/***************TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*********/

/*
*			STELLAR_steady_state_wave_curve():
*
*	Calculates the steady state wave curve through the state
*	st0 with steady state flow speed q0,  parameterized
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
*	where gamma0 = adiabatic_gamma(st0), Mn = m/(rho0*c0) =
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

LOCAL	boolean	STELLAR_steady_state_wave_curve(
	double p1,
	double M0sq,
	double *theta,
	Locstate st0)
{
	double	p0 = pressure(st0);
	double	dp, tmp;
	double	c4, cf8;
	double	Mnsq, cotb;
	double	A0, A1;
	double	tan_theta;
	double	gam0;

	if (p1 >= p0)		/* shock */
	{
	    dp = (p1 - p0) / p0;
	    gam0 = Local_gamma(st0);
	    tan_theta = dp / (gam0*M0sq - dp);
	    c4 = Coef4(st0);
	    cf8 = Coef8(st0);
	    Mnsq = (p1/p0 + c4) / cf8;
	    cotb = M0sq/Mnsq - 1.0;
	    cotb = max(0.0,cotb);
	    cotb = sqrt(cotb);
	    tan_theta *= cotb;
	    *theta = (tan_theta > 0.0) ? atan(tan_theta) : 0.0;
	    return FUNCTION_SUCCEEDED;
	}
	else		/* rarefaction */
	{
	    double mu;
	    double GC3 = Coef3(st0);

	    c4 = Coef4(st0);
	    mu = Mu(st0);
	    if (M0sq < SONIC_MINUS_SQR)
	    {
	    	(void) printf("WARNING in STELLAR_steady_state_wave_curve(), "
	                      "Subsonic state in rarefaction\n");
	        (void) printf("Mach number = %g, (squared %g)\n",
	                      sqrt(M0sq),M0sq);
	        return FUNCTION_FAILED;
	    }
	    else
	    {
	    	M0sq = max(1.0,M0sq);
	    	A0 = asin(sqrt(1.0/M0sq));
	    	A1 = atan(mu*pow(p1/p0,GC3)/
	             sqrt(1.0 + c4*(M0sq - 1.0) - pow(p1/p0,2.0*GC3)));
	        tmp = angle(tan(A1),mu);
	        *theta = -(A1 - A0 + (tmp - atan(mu/tan(A0)))/mu);
	        return FUNCTION_SUCCEEDED;
	    }
	}
}		/*end STELLAR_steady_state_wave_curve*/


/*
*			STELLAR_pressure_at_sonic_point():
*
*	Returns the pressure at the sonic point of the shock polar
*	through the state st0 with steady state Mach number M0.
*/

LOCAL	double	STELLAR_pressure_at_sonic_point(
	double M0sq,
	Locstate st0)
{
	double x;

	x = 0.5*(M0sq - 1.0);
	return pressure(st0) * (x + sqrt(1.0 + 2.0*Coef4(st0)*x + x*x));
}		/*end STELLAR_pressure_at_sonic_point*/


/*
*			STELLAR_pr_at_max_turn_angle():
*
*	Given st0 and the Mach number (squared) of st0 in the frame
*	of a shock, this function calculates the pressure at the point of
*	maximum turning angle on the turning angle pressure shock polar
*	through st0.
*
*	Returns FUNCTION_SUCCEEDED if sucessful,  FUNCTION_FAILED otherwise.
*/

LOCAL	boolean	STELLAR_pr_at_max_turn_angle(
	double *prm,
	double M0sq,	/* Mach number of st0 in the frame of the shock */
	Locstate st0)
{
	double xi;
	double c4;

	if (M0sq < SONIC_MINUS_SQR)
	{
	    (void) printf("WARNING in STELLAR_pr_at_max_turn_angle(), "
	                  "subsonic ahead state\n");
	    return FUNCTION_FAILED;
	}

	M0sq = max(1.0,M0sq);
	c4 = Coef4(st0);

	xi = 0.5 * (M0sq - 4.0) + 
	sqrt(0.25 * sqr(M0sq-4.0) + 2.0*(1.0+c4)*(M0sq-1.0) );
	*prm = pressure(st0)*(xi+1.0);
	return FUNCTION_SUCCEEDED;
}		/*end STELLAR_pr_at_max_turn_angle*/

/*
*		STELLAR_state_in_prandtl_meyer_wave():
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
*	The returned st1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	st1 can be used to evaluate the pressure, density, and sound speed
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
LOCAL	double	STELLAR_state_in_prandtl_meyer_wave(
	double w,
	double A_a,
	Locstate sta,
	double A_b,
	Locstate stb,
	Locstate stm,
	int stype_m)
{
	double	c_a, c_m, p_a, A_m;
	double	cos_muw, sin_muw, cot_A_a;
	double	c2, c3, mu;

	zero_state_velocity(stm,Params(sta)->dim);

	set_type_of_state(stm,TGAS_STATE);
	p_a = pressure(sta);
	c_a = sound_speed(sta);
	c2 = 1.0 / Coef2(sta);
	c3 = 1.0 / Coef3(sta);
	mu = Mu(sta);
	cos_muw = cos(mu*w);
	sin_muw = sin(mu*w);
	cot_A_a = 1.0/tan(A_a);
	c_m = cos_muw + mu*cot_A_a*sin_muw;
	A_m = mu*c_m/(mu*cos_muw*cot_A_a - sin_muw);
	if (c_m <= 0.0)
	{
	    /* rarefaction to vacuum */
	    state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
	    c_m = 0.0;
	    A_m = 0.0;
	}
	else
	{
	    Set_params(stm,sta);
	    Dens(stm) = Dens(sta)*pow(c_m/c_a,c2);
	    Press(stm) = p_a*pow(c_m/c_a,c3);
	    reset_gamma(stm);
	}

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Press(stm) < Min_pressure(sta))
	{
	    state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
	    A_m = 0.0;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	set_state(stm,stype_m,stm);
	return A_m;
}		/* end STELLAR_state_in_prandtl_meyer_wave*/

/***************END TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*****/

#if defined(COMBUSTION_CODE)
/***************DETONATION SPECIFIC UTILITY FUNCTIONS*********************/

/*
*			STELLAR_CJ_state():
*
* 	This routine finds the state behind a CJ_detonation.
*	The inputs are the initial state "start"
*	and the side (l_or_r, -1 or 1) we are on.
*/

/*ARGSUSED*/
LOCAL	double	STELLAR_CJ_state(
	Locstate CJ,
	int st_type_CJ,
	Locstate start,
	int l_or_r,
	int avail)
{
	screen("ERROR in STELLAR_CJ_state(), nonreactive EOS\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}	/* end STELLAR_CJ_state*/


/*
*	 		STELLAR_progress_state(): 
*
*	Finds the gas state as a function of reaction progress
*	in the steady frame.  
*/

/*ARGSUSED*/
LOCAL	void	STELLAR_progress_state(
	double prog,		/* reaction progress */
	Locstate init,		/* init = state behind front */
	Locstate ans,		/* TGas states, init = state behind front */
	double max_vol)		/* maximum allowed volume of reacted state */
{
	screen("ERROR in STELLAR_progress_state(), nonreactive EOS\n");
	clean_up(ERROR);
}	/* end STELLAR_progress_state*/

/*
*			STELLAR_fprint_combustion_params():
*
*	Prints combustion related parameters.
*/

/*ARGSUSED*/
LOCAL	void	STELLAR_fprint_combustion_params(
	FILE *file,
	Gas_param *params)
{
	(void) fprintf(file,"\tcomposition_type = %d %s\n",
	               PURE_NON_REACTIVE,"PURE_NON_REACTIVE");
}
/***************END DETONATION SPECIFIC UTILITY FUNCTIONS*****************/
#endif /* defined(COMBUSTION_CODE) */


/***************METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**************/

/*
*			STELLAR_neumann_riem_inv_moc():
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
*		                            = 1 for cylindrical symmetry
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

LOCAL	void	STELLAR_neumann_riem_inv_moc(
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
	screen("Entering STELLAR_neumann_riem_inv_moc()\n");
	screen("Code needed\n");
	clean_up(ERROR);
}		/*end STELLAR_neumann_riem_inv_moc*/

/*
*		STELLAR_shock_ahead_state_riem_inv_moc():
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
*		                   /ans(position = pt0 + W*dt)
*		                  /  + 0 -   
*		                 /    +  0  -   
*		                /      +   0   -   
*		          shock/        +    0    -   
*		         front/          +     0     -   
*		             /            +      0      -   
*			____/st0__________st3____st2____st1____
*		           pt0
*
*		side = NEGATIVE_SIDE
*	                                  \
*			                ans\(position = pt0 + W*dt)
*			             + 0 -  \
*			          +  0  -    \shock
*			       +   0   -      \front
*			    +    0    -        \
*		         +     0     -          \
*                     +      0      -            \
*		____st3____st2____st1__________st0\______
*	                                         pt0
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
LOCAL	void	STELLAR_shock_ahead_state_riem_inv_moc(
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
	double     u1, u2, u3;
	double     S1, S2, S3;
	double     u;
	double     c, c1, c3, c1sqr, c2sqr, c3sqr;
	double     pr, p2, rho2;
	double     alpha1, alpha3;
	double     gam, GAM, cv, cf6;
	double     ptans[3], pt1[3], pt3[3];
	double     gans[3], g1[3], g3[3], g1_bar, g3_bar;
	double     A1, A3;
	double     ua1, ua3;
	int       i, dim;
	double     radans = 0.0;
	double     alpha = nor[0]*rotational_symmetry();

	dim = gr->dim;
	u1 = 0.0;	u2 = 0.0;	u3 = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    u1 += nor[i] * Vel(st1)[i];
	    u2 += nor[i] * Vel(st2)[i];
	    u3 += nor[i] * Vel(st3)[i];
	    pt1[i] = pt0[i] + f1*dn*nor[i];
	    pt3[i] = pt0[i] + f3*dn*nor[i];
	    ptans[i] = pt0[i] + W[i]*dt;
	}
	if (add_source)
	{
	    eval_gravity(pt1,time,g1);
	    eval_gravity(pt3,time,g3);
	    eval_gravity(ptans,time+dt,gans);
	    for (g1_bar = 0.0, g3_bar = 0.0, i = 0; i < dim; ++i)
	    {
	        g1_bar += 0.5*(gans[i] + g1[i])*nor[i];
	        g3_bar += 0.5*(gans[i] + g3[i])*nor[i];
	    }
	}
	else
	{
	    g1_bar = 0.0;
	    g3_bar = 0.0;
	}

	c1sqr = sound_speed_squared(st1);
	S1 = entropy(st1);
	c1 = sqrt(c1sqr);

	rho2 = Dens(st2);
	p2 = pressure(st2);
	c2sqr = sound_speed_squared(st2);
	S2 = entropy(st2);

	c3sqr = sound_speed_squared(st3);
	S3 = entropy(st3);
	c3 = sqrt(c3sqr);

	Set_params(ans,st0);
	gam = Local_gamma(st0);
	GAM = GAMMA(st0);
	cv = Cv(st0);
	cf6 = Coef6(st0);

	alpha1 = (S2 - S1)/(4.0*gam*cv);
	alpha3 = (S2 - S3)/(4.0*gam*cv);
	A1 = 0.5*GAM*(u1 + g1_bar*dt) - c1*(1.0 + alpha1);
	A3 = 0.5*GAM*(u3 + g3_bar*dt) + c3*(1.0 + alpha3);

	if (is_rotational_symmetry() && add_source && alpha != 0.0)
	{
	    double rmin = fabs(pos_radius(0.0,gr));
	    double ra = pos_radius(ptans[0],gr);
	    double r1 = pos_radius(pt1[0],gr);
	    double r3 = pos_radius(pt3[0],gr);

	    radans = (fabs(ra) > rmin) ? 0.5*alpha*dt/ra : 0.0;

	    if (fabs(r1) > rmin)
	        A1 += 0.25*GAM*alpha*c1*u1*dt/r1;
	    if (fabs(r3) > rmin)
	        A3 -= 0.25*GAM*alpha*c3*u3*dt/r3;
	}

	if (is_rotational_symmetry() && radans != 0.0)
	{
	    double A, B, C, D;

	    A = radans*(alpha1 - alpha3);
	    B = 2.0 - (alpha1 + alpha3) +  radans*(A1 + A3);
	    C = A3 - A1;
	    D = 1.0 - 4.0*A*C/(B*B);
	    if (D < 0.0)
	    {
		screen("ERROR in POLY_shock_ahead_state_riem_inv_moc(), "
		       "negative discriminate\n");
		clean_up(ERROR);
	    }
	    c = (C/B)/(0.5*(1.0 + sqrt(D)));
	}
	else
	    c = (A3 - A1)/(2.0 - (alpha1 + alpha3));

	Dens(ans) = rho2 * pow(c*c/c2sqr,cf6);
	pr = p2 * pow(Dens(ans)/rho2,gam);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(ans));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	Press(ans) = pr;
	set_type_of_state(ans,TGAS_STATE);
	ua1 = (2.0/GAM)*(A1 + c*(1.0 - alpha1));
	ua3 = (2.0/GAM)*(A3 - c*(1.0 - alpha3));
	if (is_rotational_symmetry() && radans != 0.0)
	{
	    ua1 /= 1.0 - radans*c;
	    ua3 /= 1.0 + radans*c;
	}
	u = 0.5*(ua1 + ua3);
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = Vel(st2)[i] + (u - u2)*nor[i];
	reset_gamma(ans);
}		/*end STELLAR_shock_ahead_state_riem_inv_moc*/


/*
*			STELLAR_shock_moc_plus_rh():
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
*		nor	wall normal
*		W	first prediction of shock front velocity
*		front   Front structure
*	Output variables:
*		ans	times updated state behind shock
*		W	updated shock speed
*
*		The answer state ans is returned in Gas format.
*/


LOCAL	boolean	STELLAR_shock_moc_plus_rh(
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
	RECT_GRID   *gr = front->rect_grid;
	double       time = front->time;
	double       dt = front->dt;
	double       g, ga[3], gb[3];
	double	    hld_va[MAXD], hld_vb[MAXD];
	double	    ua, ub;
	double	    Ws, u_ans;
	double	    nr = nor[0];
	double       ca;
	double       pa, ra, pb, pm, rb, xa, xl, xr, pi;
	double       sgn;
	double       p, r, u;
	double       f0;
	double       epsilon, delta;
	const double eps = 10.0*MACH_EPS;/*TOLERANCE*/
	double       c4, gam, eta, psi;
	double       pta[3], ptb[3];
	double       fpi;
	STELLAR_S_MPRH Smprh;
	double       d;
	double       cf2, cf3;
	boolean        status = FUNCTION_SUCCEEDED;
	int	    i, dim;
	static const int   mnth = 100; /*TOLERANCE*/
	double alpha = rotational_symmetry();

	sgn = (is_forward_wave(w_type)) ? 1.0 : -1.0;
	dim = gr->dim;

	if (is_rarefaction_wave(w_type))
	{
	    set_state(ans,GAS_STATE,sta);
	    for (ua = 0.0, i = 0; i < dim; ++i)
	    	ua += nor[i] * Vel(sta)[i];
	    ca = sound_speed(sta);
	    Ws = ua + sgn*ca;
	    for (i = 0; i < dim; ++i)
	    	W[i] = 0.5 * (W[i] + Ws * nor[i]);
	    return status;
	}

	ua = ub = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    hld_va[i] = Vel(sta)[i];
	    hld_vb[i] = Vel(stb)[i];
	    ua += nor[i] * Vel(sta)[i];
	    ub += nor[i] * Vel(stb)[i];
	    pta[i] = pt[i] + W[i]*dt;
	    ptb[i] = pt[i] - dn*nor[i];
	}
	eval_gravity(pta,time + dt,ga);
	eval_gravity(ptb,time,gb);
	g = 0.5*(scalar_product(ga,nor,dim) + scalar_product(gb,nor,dim));

	Vel(sta)[0] = ua;
	Vel(stb)[0] = ub;


	ra = Dens(sta);		pa = Press(sta);	ua = Vel(sta)[0];
	rb = Dens(stb);		pb = Press(stb);	ub = Vel(stb)[0];
	pm = Press(stm);

	gam = Local_gamma(sta);
	cf2 = Coef2(sta);
	cf3 = Coef3(sta);
	c4 = Coef4(sta);

	Smprh.beta = 0.25/gam;
	Smprh.gam = gam;
	Smprh.c2 = cf2;
	Smprh.c4 = c4;
	Smprh.ca = ca = sound_speed(sta);
	Smprh.cb = sound_speed(stb);
	Smprh.psi = psi = cf3;		 /* (gam - 1)/(2*gam) */
	Smprh.eta = eta = 1.0 - psi;	 /* (gam + 1)/(2*gam) */
	Smprh.alpha = ca*psi;
	Smprh.dS = Smprh.beta * (log(pa/pb) - gam*log(ra/rb));

	f0 = ub - ua;
	if (is_gravity() == YES)
	    f0 += g*dt;
	f0 *= sgn;

	if (is_rotational_symmetry() && alpha > 0.0)
	{
	    double rmin, rada, radb;

	    Smprh.ua = ua;
	    Smprh.delta = sgn/cf2;
	    rmin = fabs(pos_radius(0.0,gr));

	    Ws = scalar_product(W,nor,dim);
	    rada = pos_radius(pta[0],gr);
	    radb = pos_radius(ptb[0],gr);

	    if ((fabs(rada) > rmin) && (fabs(radb) > rmin))
	    {
	        double cb;
	        Smprh.cyl_coord = YES;
	        cb = sound_speed(stb);
	        f0 -= 0.5*nr*cb*ub*dt*alpha/radb;
	        Smprh.rad = 0.5*nr*cf2*dt*alpha/rada;
	    }
	    else
	    {
	        Smprh.cyl_coord = NO;
	        Smprh.rad = 0.0;
	    }
	}
	else
	    Smprh.cyl_coord = NO;

	xa = pm/pa;
	xa = max(1.0,xa);

	f0 *= cf2;
	delta = max(xa*EPS, eps);
	epsilon = max(f0*EPS, eps);
	xl = max(1.0 - delta,0.5*xa);/*TOLERANCE*/
	xr = 1.5*xa;/*TOLERANCE*/

	if ((find_root(stellar_fS,(POINTER) &Smprh,f0,&pi,
		       xl,xr,epsilon,delta) == FUNCTION_FAILED) 
	    &&
	    (search_harder_for_root(stellar_fS,(POINTER) &Smprh,f0,&pi,
				    xl,xr,&xl,&xr,1.0,HUGE_VAL,
				    mnth,epsilon,delta) == FUNCTION_FAILED))
	{
	    double rpi;
	    if (find_root(stellar_fS,(POINTER) &Smprh,f0,&rpi,
		          0.0,1.0,epsilon,delta) == FUNCTION_FAILED)
	    {
	        (void) printf("WARNING in STELLAR_shock_moc_plus_rh(), "
	                      "for %s shock solution\n"
	                      "no root found for stellar_fS()\n",
		              (is_backward_wave(w_type)) ?
			      "backward" : "forward");
	        print_general_vector("pt =",pt,dim,"\n");
	        verbose_print_state("sta",sta);
	        verbose_print_state("stm",stm);
	        verbose_print_state("stb",stb);
	        (void) printf("dn = %g, ",dn);
	        print_general_vector("nor =",nor,dim,"\n");
	        (void) printf("g = %g, dt = %g\n",g,dt);
	        print_wave_type("w_type = ",w_type,"\n",current_interface());
	        (void) printf("interval = [%g, %g]\n",xl,xr);
	        (void) printf("epsilon = %g, delta = %g\n",epsilon,delta);
	        (void) stellar_fS(pi,&fpi,(POINTER) &Smprh);
	        (void) printf("f0 = %g, stellar_fS(pi) = %g, "
			      "stellar_fS(pi) - f0 = %g\n",f0,fpi,fpi-f0);
	        print_function_values(stellar_fS,(POINTER) &Smprh,
	    	                      f0,0.0,xr,100,"stellar_fS",stdout);
		status = FUNCTION_FAILED;
	    }
	    else
	    {
	        double frpi;
	        (void) printf("WARNING in STELLAR_shock_moc_plus_rh(), "
			      "shock reduced to zero strength\n");
		if (stellar_fS(pi,&fpi,(POINTER)&Smprh) &&
		    stellar_fS(rpi,&frpi,(POINTER)&Smprh))
		{
		    if (fabs(frpi - f0) < fabs(fpi - f0))
		        pi = 1.0;
		}
		else
		{
		    status = FUNCTION_FAILED;
	            pi = 1.0;
	        }
	    }
	}
	if (debugging("smocprh"))
	{
	    (void) printf("Data for STELLAR_shock_moc_plus_rh(), "
	                      "for %s shock solution\n",
		              (is_backward_wave(w_type)) ?
			      "backward" : "forward");
	    print_general_vector("pt =",pt,dim,"\n");
	    verbose_print_state("sta",sta);
	    verbose_print_state("stm",stm);
	    verbose_print_state("stb",stb);
	    (void) printf("dn = %g, ",dn);
	    print_general_vector("nor =",nor,dim,"\n");
	    (void) printf("g = %g, dt = %g\n",g,dt);
	    print_wave_type("w_type = ",w_type,"\n",current_interface());
	    (void) printf("interval = [%g, %g]\n",xl,xr);
	    (void) printf("epsilon = %g, delta = %g\n",epsilon,delta);
	    (void) stellar_fS(pi,&fpi,(POINTER) &Smprh);
	    (void) printf("f0 = %g, stellar_fS(pi) = %g, "
			  "stellar_fS(pi) - f0 = %g\n",f0,fpi,fpi-f0);
	    print_function_values(stellar_fS,(POINTER) &Smprh,
	    	                  f0,0.0,xr,100,"stellar_fS",stdout);
	}
	if (pi < (1.0 - delta)) 
	{
	    (void) printf("WARNING in STELLAR_shock_moc_plus_rh(), "
	                  "Non-physical solution for R-H conditions\n");
	    (void) printf("pi = %g\n",pi);
	    print_general_vector("pt =",pt,dim,"\n");
	    verbose_print_state("sta",sta);
	    verbose_print_state("stm",stm);
	    verbose_print_state("stb",stb);
	    (void) printf("dn = %g, ",dn);
	    print_general_vector("nor =",nor,dim,"\n");
	    (void) printf("g = %g, dt = %g\n",g,dt);
	    print_wave_type("w_type = ",w_type,"\n",current_interface());
	    (void) printf("interval = [%g, %g]\n",xl,xr);
	    (void) printf("epsilon = %g, delta = %g\n",epsilon,delta);
	    (void) printf("f0 = %g\n",f0);
	    print_function_values(stellar_fS,(POINTER) &Smprh,f0,xl,xr,
	                          100,"stellar_fS",stdout);
	    status = FUNCTION_FAILED;
	}
	if (pi < 1.0)
	    pi = 1.0;
	p = pa * pi;
	r = ra * (pi + c4)/(1.0 + c4*pi);
	d = sgn * ca / sqrt(psi + eta*pi);
	u = ua + d * (pi - 1.0) / gam;

	Dens(ans) = r;
	Vel(ans)[0] = u;
	Press(ans) = p;
	reset_gamma(ans);
	Set_params(ans,sta);
	set_type_of_state(ans,TGAS_STATE);
	Ws = ua + d * (pi + c4) / (1.0 + c4);

	for (i = 0; i < dim; ++i)
	{
	    Vel(sta)[i] = hld_va[i];
	    Vel(stb)[i] = hld_vb[i];
	}

	u_ans = Vel(ans)[0];
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = Vel(sta)[i] + (u_ans - ua) * nor[i];
	set_type_of_state(ans,TGAS_STATE);

	set_state(ans,GAS_STATE,ans);

	for (i = 0; i < dim; ++i)
	    W[i] = 0.5 * (W[i] + Ws * nor[i]);
	if (fabs((p - pm)/pm) > 0.2) /*TOLERANCE*/
	{
	    if (debugging("smocprh"))
	    {
	        (void) printf("WARNING in STELLAR_shock_moc_plus_rh(), "
	                      "change in solution too strong\n"
	                      "Incoming pressure = %g, updated pressure = %g, "
			      "relative change = %g\n",pm,p,(p-pm)/pm);
	    }
	    status = FUNCTION_FAILED;
	}
	return status;
}		/*end STELLAR_shock_moc_plus_rh*/


/***************END METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**********/

/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			STELLAR_prompt_for_state():
*
*	Prompts for a hydrodynamical state.  The form of
*	the states depends of the Eos. 	The type of the state
*	is returned.
*/

LOCAL	void	STELLAR_prompt_for_state(
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
}		/* end STELLAR_prompt_for_state */

/*
*			STELLAR_prompt_for_thermodynamics():
*
*	Prompts for the thermodynamic variables in a state.  Returns
*	a state with the appropriate thermodynamic state and zero velocity.
*	The return status gives the state type representation of the state.
*/

LOCAL	void	STELLAR_prompt_for_thermodynamics(
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
	screen("Enter the density and pressure");
	screen("%s: ",mesg);
	(void) Scanf("%f %f\n",&Dens(state),&Press(state));
}		/* end STELLAR_prompt_for_thermodynamics */


/*
*			STELLAR_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

LOCAL	void	STELLAR_fprint_EOS_params(
	FILE *file,
	Gas_param *params)
{
	int	i, nc = ((STELLAR_EOS *) params->eos)->_ionmax;
	double	*xmass = ((STELLAR_EOS *) params->eos)->_xmass;
	double	*aion = ((STELLAR_EOS *) params->eos)->_aion;
	double	*zion = ((STELLAR_EOS *) params->eos)->_zion;

	(void) fprintf(file,"\tEquation of state = %d STELLAR\n",STELLAR);
	(void) fprintf(file,"\tnumber of isotops = %d\n",nc);
	fprint_float(file,"\tmass fraction = ",xmass[0],"");
	for (i = 1; i < nc; i++)
	    fprint_float(file,"",xmass[i],"");
	fprint_float(file,"\tatomic weight = ",aion[0],"");
	for (i = 1; i < nc; i++)
	    fprint_float(file,"",aion[i],"");
	fprint_float(file,"\tatomic number = ",zion[0],"");
	for (i = 1; i < nc; i++)
	    fprint_float(file,"",zion[i],"");
	(void) fprintf(file,"\n");
}		/*end STELLAR_fprint_EOS_params */

/*
*			STELLAR_read_print_EOS_params():
*
*	Reads the equation of state data as printed by STELLAR_fprint_EOS_params.
*	This is restart function.
*/

LOCAL	void	STELLAR_read_print_EOS_params(
	INIT_DATA     *init,
        const IO_TYPE *io_type,
	Gas_param	*params)
{
	FILE		*file = io_type->file;
	STELLAR_EOS *eos = (STELLAR_EOS *) params->eos;
	int	c, i, nc;
	double	*xmass = eos->_xmass;
	double	*aion = eos->_aion;
	double	*zion = eos->_zion;

	if (fgetstring(file,"number of isotops = "))
	    fscanf(file,"%d",&nc);
	else
	{
	    screen("ERROR in STELLAR_read_print_EOS_params(),"
		   "can't find stellar printout\n");
	}
	eos->_ionmax = nc;
	

	if (fgetstring(file,"mass fraction = "))
	    for (i = 0; i < nc; i++)
		xmass[i] = fread_float(NULL,io_type);
	if (fgetstring(file,"\tatomic weight = "))
	    for (i = 0; i < nc; i++)
		aion[i] = fread_float(NULL,io_type);
	if (fgetstring(file,"\tatomic number = "))
	    for (i = 0; i < nc; i++)
		zion[i] = fread_float(NULL,io_type);
}		/*end STELLAR_read_print_EOS_params*/

/*
*			STELLAR_free_EOS_params():
*
*	Frees the storage allocated for an equation of state parameter
*	function.
*/

LOCAL	EOS*	STELLAR_free_EOS_params(
	EOS *eos)
{
        printf("Entering STELLAR_free_EOS_params(), please check!\n");
	clean_up(ERROR);
		
	
	free(eos);
	return NULL;
}		/*end STELLAR_free_EOS_params*/

/*
*			STELLAR_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

LOCAL	void	STELLAR_prompt_for_EOS_params(
	INIT_DATA  	*init,
	Gas_param	*params,
	const char	*message1,
	const char	*message2)
{
	STELLAR_EOS     *speos = (STELLAR_EOS*) params->eos;
	int             n,nc,len;
	double		*xmass = speos->_xmass;
	double		*aion = speos->_aion;
	double		*zion = speos->_zion;
	char		s[100];

	nc = 1;
	screen("Enter the number of isotops (<< %d)",MAX_NUM_ISOTOPS);
	screen("%s",(len > 25) ? "\n\t" : " ");
        screen("for the%s gas%s (default = %d): ",message1,message2,nc);
        (void) Gets(s);
        if (s[0] != '\0')
                (void) sscanf(s,"%d\n",&nc);
	if (nc <= 0 || nc > MAX_NUM_GAS_COMPS)
        {
	    screen("ERROR in STELLAR_prompt_for_EOS_params(), "
                      "invalid number of isotops %d\n"
                      "Value must be in the range 1 <= n <= %d\n",
                      nc,MAX_NUM_ISOTOPS);
	    clean_up(ERROR);
	}
	speos->_ionmax = nc;
	for (n = 0; n < nc; ++n)
        {
	    screen("Enter the mass fraction for %d%s isotop: ",
	    		n+1,ordinal_suffix(n+1));
	    (void) Scanf("%f\n",xmass+n);
	    screen("Enter the atomic weight for %d%s isotop: ",
	    		n+1,ordinal_suffix(n+1));
	    (void) Scanf("%f\n",aion+n);
	    screen("Enter the atomic number for %d%s isotop: ",
	    		n+1,ordinal_suffix(n+1));
	    (void) Scanf("%f\n",zion+n);
	}
}		/*end STELLAR_prompt_for_EOS_params*/



/***************Problem Type Specific Initialization Functions**************/

/*
*			STELLAR_RT_RS_f():
*
*	Support function for the computation of a solution to the linearized
*	Rayleigh-Taylor flow.
*
*	NEEDED:  More complete documentation
*/

LOCAL	double	STELLAR_RT_RS_f(
	double		s,
	Locstate	amb_st,
	double		dz,
	double		k_sqr,
	double		g_z)
{
	double gam, GAM;
	double h_sqr;
	double csqr = sound_speed_squared(amb_st);
	double rho = Dens(amb_st);
	double alpha1, alpha2;
	double D, D1, D2, N;
	double arg1, arg2;
	double beta;

	if (Local_gamma_set(amb_st) == NO)
	    set_local_gamma(amb_st);
	gam = Local_gamma(amb_st);
	GAM = GAMMA(amb_st);
	beta = gam * g_z / csqr;
	h_sqr = 0.25*sqr(beta) + s/csqr + k_sqr + GAM*sqr(g_z)*k_sqr/(s*csqr);
	if (h_sqr < 0.0)
	{
	    screen("ERROR in RT_RS_f(), h_sqr = %g < 0\n",h_sqr);
	    screen("s = %g, beta = %g, csqr = %g\n",s,beta,csqr);
	    screen("k_sqr = %g, gam = %g, g_z = %g\n",k_sqr,gam,g_z);
	    clean_up(ERROR);
	}
	alpha1 = -sqrt(h_sqr) - 0.5*beta;
	alpha2 =  sqrt(h_sqr) - 0.5*beta;
	if (alpha1*dz >= alpha2*dz)
	{
	    arg1 = alpha1*dz;
	    arg2 = alpha2*dz;
	    D = 1.0 - exp(arg2 - arg1);
	}
	else
	{
	    arg1 = alpha2*dz;
	    arg2 = alpha1*dz;
	    D = exp(arg2 - arg1) - 1.0;
	}
	N = GAM*sqr(g_z) + s*csqr;
	D1 = GAM*g_z + alpha1*csqr;
	D2 = GAM*g_z + alpha2*csqr;
	return  rho*N*(exp(alpha1*dz-arg1)/D1 -
	        exp(alpha2*dz-arg1)/D2)/D - g_z*rho;
}		/*end STELLAR_RT_RS_f*/


/*
*			STELLAR_RT_single_mode_perturbation_state():
*
*	Computes the perturbation term for the solution to the linearized
*	Euler equations in the Rayleigh-Taylor problem.  See the appendix to
*
*	``The Dynamics of Bubble Growth for
*				Rayleigh-Taylor Unstable Interfaces''
*	C. L. Gardner, J. Glimm, O. McBryan, R. Menikoff, D. H. Sharp,
*	and Q. Zhang, Phys. Fluids 31 (3), 447-465 (1988).
*
*	for an explanation of the formulas.
*
*       To avoid overflows in some exponential terms, the solution is divided
*       into the light and heavy fluid cases, as determined by the sign of
*       zh = z_bdry - z_intfc.  For the case where zh > 0 (light fluid), the
*       equation for the pressure perturbation dP is multiplied by
*       exp(alpham*zh)/exp(alpham*zh), and for the case where zh < 0 (heavy
*       fluid), it's multiplied by exp(alphap*zh)/exp(alphap*zh).  All other
*       formulas are exactly as they appear in Gardner et al, with the
*       exception of a sign difference in the density perturbation which
*       corrects a typo in the paper.
*
*       Note that ans is only the perturbation to the ambient state.
*/

/*ARGSUSED*/
LOCAL	void	STELLAR_RT_single_mode_perturbation_state(
	Locstate	ans,
	double		*coords,
	double		t,
	Locstate	amb_st,
	double		z_intfc,
	double		z_bdry,
	MODE		*mode,
	double		g_z)
{
	int   j, dim = Params(amb_st)->dim;

	double gam;
	double GAM;
	double csqr;
	double beta;

	double A = mode->amplitude;
	double sigma = mode->growth_rate;
	double timefac = exp(sigma*t);
	double *k = mode->wave_number;
	double phi = scalar_product(k,coords,dim-1) - mode->phase;

	double h, alphap, alpham;
	double z, zh, N, D, Dp, Dm;
	double prefac, expp, expm;
	double dP, dPdz, drho, dv_z;

	if (Local_gamma_set(amb_st) == NO)
	    set_local_gamma(amb_st);
	gam = Local_gamma(amb_st);
	GAM = GAMMA(amb_st);
	csqr = sound_speed_squared(amb_st);
	beta = gam * g_z / csqr;

	h = 0.25*sqr(beta) + sqr(sigma)/csqr +
	  scalar_product(k,k,dim-1)*(1.0+GAM*sqr(g_z)/(sqr(sigma)*csqr));
	h = sqrt(h);

	alphap = -0.5*beta + h;
	alpham = -0.5*beta - h;

	z = coords[dim-1];
	zh = z_bdry - z_intfc;
	N = GAM * sqr(g_z) + sqr(sigma) * csqr;
	Dp = GAM * g_z + alphap * csqr;
	Dm = GAM * g_z + alpham * csqr;

	/* To avoid overflows in the exponential terms, we have to
	   separate the light (above the interface) and heavy (below
	   the interface) fluid cases, as determined by the sign of
	   zh. */

	if (zh > 0.0)      /* light fluid */
	{
	    D = 1.0 - exp(-2.0*h*zh);
	    expm = exp((alpham+beta) * (z-z_intfc)) / Dm;
	    expp = exp((alpham+beta) * zh +
	           (alphap+beta) * (z-z_bdry)) / Dp;
	}
	else               /* heavy fluid */
	{
	    D = exp(2.0*h*zh) - 1.0;
	    expm = exp((alpham+beta) * (z-z_bdry) +
	           (alphap+beta) * zh) / Dm;
	    expp = exp((alphap+beta) * (z-z_intfc)) / Dp;
	}
	prefac = N * Dens(amb_st) * A / D;

	dP = -prefac * (expm - expp);
	dPdz = -prefac * ((alpham+beta) * expm - (alphap+beta) * expp);
	drho = (sqr(sigma) * dP + GAM * g_z * dPdz) / N;
	dv_z = (g_z * drho - dPdz) / (sigma * Dens(amb_st));

	set_type_of_state(ans,TGAS_STATE);
	Set_params(ans,amb_st);

	Press(ans) = dP*timefac*sin(phi);
	Dens(ans) = drho*timefac*sin(phi);
	Vel(ans)[dim-1] = dv_z*timefac*sin(phi);
	reset_gamma(ans);

	for (j = 0; j < dim-1; ++j)
	    Vel(ans)[j] = -k[j]*dP*cos(phi)/(sigma*Dens(amb_st));

	if (debugging("RT_RS_all")) 
	{
	    (void) printf("\nIn STELLAR_RT_single_mode_perturbation_state(),\n");
	    (void) printf("A = %g; sigma = %g; timefac = %g; "
	                  "k[0] = %g; xfac = %g\n",
	                  A,sigma,timefac,k[0],sin(phi));
	    (void) printf("gam = %g; csqr = %g; beta = %g\n",gam,csqr,beta);
	    (void) printf("h = %g; alphap = %g; alpham = %g\n",h,alphap,alpham);
	    (void) printf("zh = %g; N = %g; Dp = %g; Dm = %g\n",zh,N,Dp,Dm);
	    (void) printf("D = %g; expm = %g; expp = %g; prefac = %g\n",
	                  D,expm,expp,prefac);
	    (void) printf("dP = %g; dPdz = %g; drho = %g; dv_z = %g\n",
	                  dP,dPdz,drho,dv_z);
	    (void) printf("x = %g; z = %g; t = %g; z_intfc = %g; z_bdry = %g\n",
	                  coords[0],coords[1],t,z_intfc,z_bdry);
	    (void) printf("P = %g; d = %g; Vx = %g; Vz = %g\n",
	                  Press(amb_st)+Press(ans),Dens(amb_st)+Dens(ans),
	                  Vel(amb_st)[0]+Vel(ans)[0],
	                  Vel(amb_st)[dim-1]+Vel(ans)[dim-1]);
	}	
}		/*end STELLAR_RT_single_mode_perturbation_state*/


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

LOCAL	void	STELLAR_KH_single_mode_state(
	Locstate	ans,
	double		*coords,
	double		t,
	Locstate	amb_st,
	double		stream_velocity,
	double		z_intfc,
	double		z_bdry,
	MODE		*mode)
{
	double gam;
	double amb_press = pressure(amb_st);

	Set_params(ans,amb_st);
	set_type_of_state(ans,TGAS_STATE);
	gam = Local_gamma(amb_st);
	                /*TODO 3D*/
	FORTRAN_NAME(khstate)(coords,coords+1,&t,&Dens(amb_st),&z_intfc,
	                      &z_bdry,&mode->amplitude,mode->wave_number,
	                      &mode->growth_rate,&mode->propagation_speed,
	                      &mode->phase,Vel(ans),Vel(ans)+1,
	                      &Dens(ans),&Press(ans),&amb_press,
	                      &stream_velocity,&gam);
}		/*end STELLAR_KH_single_mode_state*/


/*
*		STELLAR_compute_isothermal_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isothermal one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, T = T_r.
*/

LOCAL	void	STELLAR_compute_isothermal_stratified_state(
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	double gam, tmp;
	double pr, rho;
	int dim = Params(ref_st)->dim;

	if (Local_gamma_set(ref_st) == NO)
	    set_local_gamma(ref_st);
	gam = Local_gamma(ref_st);
	Set_params(ans,ref_st);
	tmp = exp((gam*gz/sound_speed_squared(ref_st))*dz);
	rho = Dens(ref_st);
	pr = pressure(ref_st);
	Dens(ans) = rho*tmp;
	Press(ans) = pr*tmp;
	set_type_of_state(ans,TGAS_STATE);
	zero_state_velocity(ans,dim);
	reset_gamma(ans);
}		/*end STELLAR_compute_isothermal_stratified_state */


/*
*		STELLAR_compute_isentropic_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isentropic one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, S = S_r.
*/

LOCAL	void	STELLAR_compute_isentropic_stratified_state(
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	double rho_r, p_r, gam, GAM;
	double cr2;
	double csq_ratio;

	if (Local_gamma_set(ref_st) == NO)
	    set_local_gamma(ref_st);
	Set_params(ans,ref_st);
	gam = Local_gamma(ref_st);
	GAM = GAMMA(ref_st);
	cr2 = sound_speed_squared(ref_st);
	csq_ratio = 1.0 + GAM*gz*dz/cr2;
	if (csq_ratio <= 0.0)
	{
	    Dens(ans) = Vacuum_dens(ref_st);
	    Press(ans) = Min_pressure(ref_st);
	}
	else
	{
	    rho_r = Dens(ref_st);
	    p_r = pressure(ref_st);
	    Dens(ans) = rho_r*pow(csq_ratio,1.0/GAM);
	    Press(ans) = p_r*pow(csq_ratio,gam/GAM);
	}
	reset_gamma(ans);
	set_type_of_state(ans,TGAS_STATE);
	zero_state_velocity(ans,Params(ref_st)->dim);
}	/*end STELLAR_compute_isentropic_stratified_state*/


/***************End Problem Type Specific Initialization Functions**********/
/***************END INITIALIZATION UTILITY FUNCTIONS************************/

LOCAL	boolean	stellar_fS(double x, double *fans, POINTER prm)
{
	double c, rr, v;
	double ca, cb, dS, alpha, beta, c4, eta, psi;

	ca = ((STELLAR_S_MPRH *) prm)->ca;
	cb = ((STELLAR_S_MPRH *) prm)->cb;
	dS = ((STELLAR_S_MPRH *) prm)->dS;
	alpha = ((STELLAR_S_MPRH *) prm)->alpha;
	beta = ((STELLAR_S_MPRH *) prm)->beta;
	eta = ((STELLAR_S_MPRH *) prm)->eta;
	psi = ((STELLAR_S_MPRH *) prm)->psi;
	c4 = ((STELLAR_S_MPRH *) prm)->c4;

	if (x >= 1.0)
	{
	    rr = (x + c4)/(1.0 + c4*x);
	    v = alpha * (x - 1.0)/sqrt(psi + eta*x);
	    c = ca * sqrt(x/rr);
	    dS += beta*log(x) - 0.25*log(rr);
	}
	else
	{
	    double c2, gam;
	    gam = ((STELLAR_S_MPRH *) prm)->gam;
	    c2 = ((STELLAR_S_MPRH *) prm)->c2;
	    rr = pow(x,1.0/gam);
	    c = ca * pow(x,(gam-1.0)/gam);
	    v = (c - ca)/c2;
	}

	*fans = c - cb + v - (c + cb)*dS;
	if (is_rotational_symmetry() && ((STELLAR_S_MPRH *) prm)->cyl_coord == YES)
	{
	    double ua, rad, delta;

	    ua = ((STELLAR_S_MPRH *) prm)->ua;
	    rad = ((STELLAR_S_MPRH *) prm)->rad;
	    delta = ((STELLAR_S_MPRH *) prm)->delta;

	    *fans = *fans + rad*c*(ua + delta*v);
	}
	return FUNCTION_SUCCEEDED;
}		/*end stellar_fS*/

/*TMP */
LOCAL	double max_gamma = 0;
LOCAL	double min_gamma = 2;
LOCAL	void set_local_gamma(
	Locstate state)
{
        int     ionmax, input;
        double   *xmass, *aion, *zion;
        double   den, temp, ener, pres,entropy, csond; 
	double   dpd, dpt, c_v, c_p, gammac, gam2, gam3;
	int	len = strlen(gaspath);

        ionmax = ionmax(state);
        xmass = xmass(state);
        aion = aion(state);
        zion = zion(state);

	den = Dens(state);
	switch (state_type(state))
	{
	case GAS_STATE:
	    input = 2;
	    ener = (Energy(state) - kinetic_energy(state))/Dens(state);
	    FORTRAN_NAME(helmstate)(&den,&temp,&ener,&pres,&entropy, &csond,
                                 &dpd, &dpt, &c_v, &c_p, &gammac, &gam2,
                                 &gam3, &ionmax, xmass, aion, zion, &input,
				 gaspath,&len);
	    break;

	case EGAS_STATE:
	    input = 2;
	    ener = Energy(state);
	    FORTRAN_NAME(helmstate)(&den,&temp,&ener,&pres,&entropy, &csond,
                                 &dpd, &dpt, &c_v, &c_p, &gammac, &gam2,
                                 &gam3, &ionmax, xmass, aion, zion, &input,
				 gaspath,&len);
	    break;
	       
	case FGAS_STATE:
	    input = 1;
	    temp = Temperature(state);
	    FORTRAN_NAME(helmstate)(&den,&temp,&ener,&pres,&entropy, &csond,
                                 &dpd, &dpt, &c_v, &c_p, &gammac, &gam2,
                                 &gam3, &ionmax, xmass, aion, zion, &input,
				 gaspath,&len);
	    break;
	       
	case TGAS_STATE:
	case VGAS_STATE:
	    input = 3;
	    pres = Press(state);
	    FORTRAN_NAME(helmstate)(&den,&temp,&ener,&pres,&entropy, &csond,
                                 &dpd, &dpt, &c_v, &c_p, &gammac, &gam2,
                                 &gam3, &ionmax, xmass, aion, zion, &input,
				 gaspath,&len);
	    break;
	default:
	    screen("ERROR in Local_gamma(), no such state type\n");
	    clean_up(ERROR);
	}
	Local_gamma_set(state) = YES;
	Local_gamma(state) = gammac;
	Cv(state) = c_v;
	Cv(state) = 1/(gammac-1);
	R(state) = (gammac-1)*ener/temp;
}	/* end set_local_gamma */

