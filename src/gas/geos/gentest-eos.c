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
*				gentest-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#include <geos/gentest.h>

/* Local prototypes */

LOCAL	double	GENTEST_pressure(Locstate);
LOCAL	double	GENTEST_sound_speed_squared(Locstate);
LOCAL	double	GENTEST_specific_internal_energy(Locstate);

LOCAL	double	GENTEST_temperature(Locstate);
LOCAL	double	GENTEST_entropy(Locstate);
LOCAL	double	GENTEST_gruneisen_gamma(Locstate);

LOCAL	double	GENTEST_C_V(Locstate);

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	GENTEST_prompt_for_state(Locstate,int,Gas_param*,const char*);
LOCAL	void	GENTEST_prompt_for_thermodynamics(Locstate,Gas_param*,
						  const char*);
LOCAL	void	GENTEST_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	GENTEST_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                              Gas_param*);
LOCAL	EOS*	GENTEST_free_EOS_params(EOS*);
LOCAL	void	GENTEST_prompt_for_EOS_params(INIT_DATA*,Gas_param*,
					      const char*,const char*);

	/* Equation of state domain functions */
LOCAL	double	GENTEST_Min_energy(Locstate);
LOCAL	double	GENTEST_Min_pressure(Locstate);

LOCAL	void	set_eos_function_hooks(EOS*);
LOCAL	void	set_GENTEST_coefs(Gas_param*);


/* LOCAL structs/variables */

typedef struct {
	double	ca, cb, dS, alpha, beta, c4, psi, eta;
} GENTEST_S_MPRH;

EXPORT	EOS	*set_GENTEST_eos(
	EOS	*eos)
{
	if (eos == NULL)
	{
		scalar(&eos,sizeof(GENTEST_EOS));
	}
	(void) set_GENERIC_eos(eos);
	set_eos_function_hooks(eos);
	return eos;
}	/*end set_GENTEST_eos*/

LOCAL	void	set_GENTEST_coefs(
	Gas_param *params)
{
	double	pinf = ((GENTEST_EOS *) params->eos)->pinf;
	double	gam = ((POLY_EOS *) params->eos)->gamma;
	double	coef6;

	set_POLY_coefs(params);
	coef6 = ((POLY_EOS *) params->eos)->coef6;
	params->min_pressure = -0.95*pinf;
	params->min_energy = coef6*(params->min_pressure+gam*pinf);
}	/*end set_GENTEST_coefs*/

LOCAL	void	set_eos_function_hooks(
	EOS *eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_pressure = GENTEST_pressure;
	eos->_sound_speed_squared = GENTEST_sound_speed_squared;
	eos->_specific_internal_energy = GENTEST_specific_internal_energy;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_temperature = GENTEST_temperature;
	eos->_entropy = GENTEST_entropy;
	eos->_gruneisen_gamma = GENTEST_gruneisen_gamma;
	eos->_C_V = GENTEST_C_V;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_prompt_for_state = GENTEST_prompt_for_state;
	eos->_prompt_for_thermodynamics = GENTEST_prompt_for_thermodynamics;
	eos->_fprint_EOS_params = GENTEST_fprint_EOS_params;
	eos->_read_print_EOS_params = GENTEST_read_print_EOS_params;
	eos->_free_EOS_params = GENTEST_free_EOS_params;
	eos->_prompt_for_EOS_params = GENTEST_prompt_for_EOS_params;

	/* Equation of state domain functions */
	eos->_Min_energy = GENTEST_Min_energy;
	eos->_Min_pressure = GENTEST_Min_pressure;
}	/*end set_eos_function_hooks*/


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*			GENTEST_pressure():
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

LOCAL	double	GENTEST_pressure(
	Locstate state)
{
	double pr, ke, gam;
	if (is_obstacle_state(state))
		return HUGE_VAL;
	switch (state_type(state)) 
	{
	case	GAS_STATE:
		gam = Gamma(state);
		ke = kinetic_energy(state);
		pr = (gam-1.0)*(Energy(state) - ke) - gam*Pinf(state);
		break;

	case	EGAS_STATE:
		gam = Gamma(state);
		pr = (gam-1.0)*Energy(state)*Dens(state) - gam*Pinf(state);
		break;

	case    FGAS_STATE:
		pr = R(state)*Temperature(state)*Dens(state) - Pinf(state);
		break;

	case	TGAS_STATE:
	case	VGAS_STATE:
		pr = Press(state);
		break;

	default:
		screen("ERROR in GENTEST_pressure(), no such state type\n");
		clean_up(ERROR);
		break;
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(state));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end GENTEST_pressure*/

/*
*			GENTEST_sound_speed_squared():
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

LOCAL	double	GENTEST_sound_speed_squared(
	Locstate state)
{
	double c2 = Gamma(state)*stiff_pressure(state)/Dens(state);
	return c2;
}		/*end GENTEST_sound_speed_squared*/


/*
*			GENTEST_specific_internal_energy():
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/

LOCAL	double	GENTEST_specific_internal_energy(
	Locstate state)
{
	switch (state_type(state)) {

	case	GAS_STATE:
		return	(Energy(state) - kinetic_energy(state))/Dens(state);

	case	EGAS_STATE:
		return	Energy(state);

	case	TGAS_STATE:
		return (Press(state) + Gamma(state)*Pinf(state))/
					(Dens(state)*(Gamma(state) - 1.));
	
	case	FGAS_STATE:
		return Pinf(state)/Dens(state) + R(state)*Temperature(state)*Coef6(state);

	case	VGAS_STATE:
		return Int_en(state);

	default:
		screen("ERROR in GENTEST_specific_internal_energy(), ");
		screen("no such state type\n");
		clean_up(ERROR);
		break;
	}
	return ERROR_FLOAT;
}		/*end GENTEST_specific_internal_energy*/


/***************SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS ***********/

/*
*			GENTEST_temperature():
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	double	GENTEST_temperature(
	Locstate state)
{
	if (state_type(state) == FGAS_STATE) return Temperature(state);
	return stiff_pressure(state)/(R(state)*Dens(state));
}		/*end GENTEST_temperature*/

/*
*			GENTEST_entropy():
*
*	Returns the specific entropy of a state.
*/

LOCAL	double	GENTEST_entropy(
	Locstate state)
{
	if (state_type(state) == VGAS_STATE)
	    return Entropy(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(state) < Vacuum_dens(state))
	    return 0.0;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	return (log(stiff_pressure(state)) - Gamma(state)*log(Dens(state))) *
			Coef6(state)*R(state);
}		/*end GENTEST_entropy*/

/*
*			GENTEST_gruneisen_gamma():
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

LOCAL	double	GENTEST_gruneisen_gamma(
	Locstate state)
{
	return Gamma(state)-1.0;
}		/*end GENTEST_gruneisen_gamma*/

/*
*			GENTEST_C_V():
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*/

LOCAL	double	GENTEST_C_V(
	Locstate state)
{
	return R(state)/(Gamma(state) - 1.0);
}	/* end GENTEST_C_V */


/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			GENTEST_prompt_for_state():
*
*	Prompts for a hydrodynamical state.  The form of
*	the states depends of the Eos. 	The type of the state
*	is returned.
*/

LOCAL	void	GENTEST_prompt_for_state(
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
}		/* end GENTEST_prompt_for_state */

/*
*			GENTEST_prompt_for_thermodynamics():
*
*	Prompts for the thermodynamic variables in a state.  Returns
*	a state with the appropriate thermodynamic state and zero velocity.
*	The return status gives the state type representation of the state.
*/

LOCAL	void	GENTEST_prompt_for_thermodynamics(
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
}		/* end GENTEST_prompt_for_thermodynamics */

/*
*			GENTEST_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

LOCAL	void	GENTEST_fprint_EOS_params(
	FILE *file,
	Gas_param *params)
{
	double gam = ((POLY_EOS *) params->eos)->gamma;
	double R = ((POLY_EOS *) params->eos)->R;
	double pinf = ((GENTEST_EOS *) params->eos)->pinf;

	(void) fprintf(file,"\tEquation of state = %d GENTEST\n",
		GENTEST);
	(void) fprintf(file,"\tgamma = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &gam,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,gam);
	(void) fprintf(file,", pinfinity = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &pinf,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,pinf);
	(void) fprintf(file,", R = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &R,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,R);
	(void) fprintf(file,"\n");
}		/*end GENTEST_fprint_EOS_params */

/*
*			GENTEST_read_print_EOS_params():
*
*	Reads the equation of state data as printed by GENTEST_fprint_EOS_params.
*	This is restart function.
*/

/*ARGSUSED*/
LOCAL	void	GENTEST_read_print_EOS_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	FILE      *file = io_type->file;
	GENTEST_EOS *speos = (GENTEST_EOS *)params->eos;
	POLY_EOS *peos = (POLY_EOS *)params->eos;
	int c;

	peos->gamma = fread_float("gamma = ",io_type);
	speos->pinf = fread_float("pinfinity = ",io_type);
	if ((c = getc(file)) == ',')
	{
	    peos->R = fread_float("R = ",io_type);
	}
	else
	{
	    (void) ungetc(c,file);
	    peos->R = 1.0;
	}
	set_GENTEST_coefs(params);
}		/*end GENTEST_read_print_EOS_params*/

/*
*			GENTEST_free_EOS_params():
*
*	Frees the storage allocated for an equation of state parameter
*	function.
*/

LOCAL	EOS*	GENTEST_free_EOS_params(
	EOS *eos)
{
	free(eos);
	return NULL;
}		/*end GENTEST_free_EOS_params*/

/*
*			GENTEST_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

/*ARGSUSED*/
LOCAL	void	GENTEST_prompt_for_EOS_params(
	INIT_DATA  *init,
	Gas_param  *params,
	const char *message1,
	const char *message2)
{
	GENTEST_EOS *speos = (GENTEST_EOS*) params->eos;
	POLY_EOS *peos = (POLY_EOS*) params->eos;
	char s[2048];
	long offset;
	int n;
	static const char *fmt = "%lf %lf %lf";

	peos->R = 1;
	screen("Enter the ratio of specific heats (gamma), \n");
	screen("\tthe stiffened gas constant p infinity, and\n");
	screen("\tthe ideal gas constant (R, (P+P_inf)V = RT, ");
	screen("default for R = %g)\n",peos->R);
	screen("\tfor the%s gas%s: ",message1,message2);
	offset = ftell(stdout);
	(void) Gets(s);
	n = sscanf(s,fmt,&peos->gamma,&speos->pinf,&peos->R);
	if (n == 1)
	{
	    /* Compatibility mode for old style input files */
	    (void) fseek(stdout,offset,SEEK_SET);
	    (void) fscan_float(stdin,&speos->pinf);
	    (void) fgets(s,2046,stdin);/*Grad trailing newline*/
	    screen("%g %g %g\n",peos->gamma,speos->pinf,peos->R);
	}
	set_GENTEST_coefs(params);
}		/*end GENTEST_prompt_for_EOS_params*/


/***************EQUATION OF STATE DOMAIN FUNCTIONS**************************/

LOCAL	double	GENTEST_Min_energy(
	Locstate	state)
{
#if defined(UNRESTRICTED_THERMODYNAMICS)
	return -HUGE_VAL;
#else /* defined(UNRESTRICTED_THERMODYNAMICS) */
	double	emin = Coef6(state)*
			(Min_pressure(state)+Gamma(state)*Pinf(state));
	return max(emin,Params(state)->min_energy);
#endif /* defined(UNRESTRICTED_THERMODYNAMICS) */
}	/*end GENTEST_Min_energy*/

LOCAL	double	GENTEST_Min_pressure(
	Locstate	state)
{
#if defined(UNRESTRICTED_THERMODYNAMICS)
	return -HUGE_VAL;
#else /* defined(UNRESTRICTED_THERMODYNAMICS) */
	return max(-Pinf(state),Params(state)->min_pressure);
#endif /* defined(UNRESTRICTED_THERMODYNAMICS) */
}	/*end GENTEST_Min_energy*/

