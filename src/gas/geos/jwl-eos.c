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
*				jwl-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#define	DEBUG_STRING	"jwl"
#include <geos/jwl.h>

	/* LOCAL Function Prototypes */
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL	double	JWL_internal_energy(Locstate);
LOCAL	double	JWL_pressure(Locstate);
LOCAL	double	JWL_sound_speed_squared(Locstate);
LOCAL	double	JWL_acoustic_impedance_squared(Locstate);
LOCAL	double	JWL_specific_internal_energy(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL	double	JWL_specific_enthalpy(Locstate);
LOCAL	double	JWL_temperature(Locstate);
LOCAL	double	JWL_entropy(Locstate);
LOCAL	double	JWL_adiabatic_gamma(Locstate);
LOCAL	double	JWL_gruneisen_gamma(Locstate);
LOCAL	double	JWL_fundamental_derivative(Locstate);
LOCAL	double	JWL_C_V(Locstate);
LOCAL	double	JWL_C_P(Locstate);
LOCAL	double	JWL_K_T(Locstate);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
LOCAL	double	JWL_dens_Hugoniot(double,Locstate);
LOCAL	void	JWL_state_w_pr_on_Hugoniot(Locstate,double,Locstate,int);

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL	double	JWL_dens_rarefaction(double,Locstate);
LOCAL	double	JWL_pressure_rarefaction(double,Locstate);
LOCAL	void	JWL_state_on_adiabat_with_pr(Locstate,double,Locstate,int);
LOCAL	void	JWL_state_on_adiabat_with_dens(Locstate,double,Locstate,int);

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	JWL_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	JWL_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                          Gas_param*);
LOCAL	void	JWL_prompt_for_EOS_params(INIT_DATA*,Gas_param*,
					  const char*,const char*);

	/* LOCAL Function prototypes*/
LOCAL	double	pc(double,JWL_EOS*);
LOCAL	double	pc0(double,JWL_EOS*);
LOCAL	double	pinfinity(Locstate);
LOCAL	double	dpinfinity(Locstate);
LOCAL	double	rg0(double,JWL_EOS*);
LOCAL	void	set_eos_function_hooks(EOS*);
LOCAL	void	set_JWL_coefs(JWL_EOS*);

#define	Pc(state)	pc(Dens(state),JWL_Eos(state))
#define	Pc0(state)	pc0(Dens(state),JWL_Eos(state))
#define	Rg0(state)	rg0(Dens(state),JWL_Eos(state))
#define	stiff_pressure(state)	(pressure(state) + pinfinity(state))

EXPORT	EOS	*set_JWL_eos(
	EOS	*eos)
{
	if (eos == NULL)
		scalar(&eos,sizeof(JWL_EOS));
	(void) set_GENERIC_eos(eos);
	set_eos_function_hooks(eos);
	return eos;
}

LOCAL	void	set_eos_function_hooks(
	EOS *eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_internal_energy = JWL_internal_energy;
	eos->_pressure = JWL_pressure;
	eos->_sound_speed_squared = JWL_sound_speed_squared;
	eos->_acoustic_impedance_squared = JWL_acoustic_impedance_squared;
	eos->_specific_internal_energy = JWL_specific_internal_energy;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_specific_enthalpy = JWL_specific_enthalpy;
	eos->_temperature = JWL_temperature;
	eos->_entropy = JWL_entropy;
	eos->_adiabatic_gamma = JWL_adiabatic_gamma;
	eos->_gruneisen_gamma = JWL_gruneisen_gamma;
	eos->_fundamental_derivative = JWL_fundamental_derivative;
	eos->_C_V = JWL_C_V;
	eos->_C_P = JWL_C_P;
	eos->_K_T = JWL_K_T;

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
	eos->_dens_Hugoniot = JWL_dens_Hugoniot;
	eos->_state_w_pr_on_Hugoniot = JWL_state_w_pr_on_Hugoniot;

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	eos->_dens_rarefaction = JWL_dens_rarefaction;
	eos->_pressure_rarefaction = JWL_pressure_rarefaction;
	eos->_state_on_adiabat_with_pr = JWL_state_on_adiabat_with_pr;
	eos->_state_on_adiabat_with_dens = JWL_state_on_adiabat_with_dens;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_fprint_EOS_params = JWL_fprint_EOS_params;
	eos->_read_print_EOS_params = JWL_read_print_EOS_params;
	eos->_prompt_for_EOS_params = JWL_prompt_for_EOS_params;
}


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*			JWL_internal_energy():
*
*	Returns the internal energy per unit volume of a state.
*/

LOCAL	double	JWL_internal_energy(
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
		
	case TGAS_STATE:
		return	(Press(state)-Pc(state))/W(state);

	case FGAS_STATE:
		return Dens(state)*specific_internal_energy(state);

	default:
		screen("ERROR: in JWL_internal_energy(), ");
		screen("no such state type\n");
		clean_up(ERROR);
	}
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end JWL_internal_energy*/


/*
*			JWL_pressure():
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

LOCAL	double	JWL_pressure(
	Locstate state)
{
	double pr, rho, E;

	if (is_obstacle_state(state))
		return HUGE_VAL;
	rho = Dens(state);
	switch (state_type(state)) 
	{

	case	GAS_STATE:
		E = Energy(state) - kinetic_energy(state);
		pr = W(state)*E + Pc(state);
		break;

	case	EGAS_STATE:
		E = Energy(state)*rho;
		pr = W(state)*E + Pc(state);
		break;

	case	FGAS_STATE:
		pr = R(state)*Temperature(state)*Dens(state) - pinfinity(state);
		break;

	case	TGAS_STATE:
	case	VGAS_STATE:
		pr = Press(state);
		break;

	default:
		screen("ERROR in JWL_pressure(), no such state type\n");
		clean_up(ERROR);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(state));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end JWL_pressure*/


/*
*			JWL_sound_speed_squared():
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

LOCAL	double	JWL_sound_speed_squared(
	Locstate state)
{
	double	c2, rho;
	double	p1, p2, rho1, rho2, e1, e2;

	if (state_type(state) == VGAS_STATE)
	{
	    double c = Sound_speed(state);
	    return c*c;
	}

	rho = Dens(state);
	p1 = P1(state);
	p2 = P2(state);
	rho1 = Rho1(state);
	rho2 = Rho2(state);
	e1 = exp(-rho1/rho);
	e2 = exp(-rho2/rho);
	c2 = (p1 * rho1 * e1 + p2 * rho2 * e2)/ (rho * rho)
		+ W(state) * Wp1(state) * (specific_internal_energy(state)
		+ dH(state) - (p1 * e1 / rho1 + p2 * e2 / rho2));
	return c2;
}		/*end JWL_sound_speed_squared*/


/*
*		JWL_acoustic_impedance_squared():
*
*	Returns the square of the local acoustic impedance of the state.
*
*                        2     dP  |
*			i = - ---- |
*                              dV  |S
*/

LOCAL	double	JWL_acoustic_impedance_squared(
	Locstate state)
{
	double	i2;
	double	rho = Dens(state);

	if (state_type(state) == VGAS_STATE)
	{
	    double i = Dens(state)*Sound_speed(state);
	    return i*i;
	}
	i2 = rho*rho*sound_speed_squared(state);
	return i2;
}		/*end JWL_acoustic_impedance_squared*/

/*
*			JWL_specific_internal_energy():
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/

LOCAL	double	JWL_specific_internal_energy(
	Locstate state)
{
	switch (state_type(state))
	{

	case	GAS_STATE:
		return	(Energy(state) - kinetic_energy(state))/Dens(state);

	case	EGAS_STATE:
		return	Energy(state);

	case	TGAS_STATE:
		return	(Press(state) - Pc(state))/(W(state)*Dens(state));

	case	FGAS_STATE:
		return	(pressure(state) - Pc(state))/(W(state)*Dens(state));
	
	case	VGAS_STATE:
		return Int_en(state);

	default:
		screen("ERROR in JWL_specific_internal_energy(), "
		       "no such state type\n");
		clean_up(ERROR);
		break;
	}
	return ERROR_FLOAT;
}		/*end JWL_specific_internal_energy*/


/***************END PRIMARY THERMODYNAMIC FUNCTIONS ************************/
/***************SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS ***********/

/*
*			JWL_specific_enthalpy():
*
*	This function computes the specific enthalpy of the given state.
*
*			H = E + P*V
*
*	E = specific internal energy, P = pressure, V = specific volume.
*
*/

LOCAL	double	JWL_specific_enthalpy(
	Locstate state)
{
#if defined(VERBOSE_PLUS_GAS)
	if (state_type(state) == VGAS_STATE)
	    return Enthalpy(state);
#endif /* defined(VERBOSE_PLUS_GAS) */

	return ((1.0 + 1.0/W(state))*pressure(state) - Pc(state)/W(state))/Dens(state);
}		/*end JWL_specific_enthalpy*/


/*
*			JWL_temperature():
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	double	JWL_temperature(
	Locstate state)
{
	if (state_type(state) == FGAS_STATE) return Temperature(state);

	return stiff_pressure(state)/(Dens(state)*R(state));
}		/*end JWL_temperature*/

/*
*			JWL_entropy():
*
*	Returns the specific entropy of a state.
*/

LOCAL	double	JWL_entropy(
	Locstate state)
{
	if (state_type(state) == VGAS_STATE)
	    return Entropy(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(state) < Vacuum_dens(state))
	    return 0.0;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	return R(state)*(log(stiff_pressure(state)) -
				Wp1(state)*log(Dens(state)))/W(state);
}		/*end JWL_entropy*/

/*
*			JWL_adiabatic_gamma():
*
*	Returns the dimensionless sound speed
*
*		gamma = - d(log P)/d(log V) | .
*					     S
*	As usual P = thermodynamic pressure,  V = specific volume
*	and S = specific entropy.
*/

LOCAL	double	JWL_adiabatic_gamma(
	Locstate state)
{
	return sound_speed_squared(state) * Dens(state)/pressure(state);
}		/*end JWL_adiabatic_gamma*/


/*
*			JWL_gruneisen_gamma():
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

LOCAL	double	JWL_gruneisen_gamma(
	Locstate state)
{
	return W(state);
}		/*end JWL_gruneisen_gamma*/

/*
*			JWL_fundamental_derivative():
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

LOCAL	double	JWL_fundamental_derivative(
	Locstate state)
{
	double	W = W(state);
	double	p = pressure(state);
	double	V = 1.0/Dens(state);
	double	pr, prp, prpp;
	double	p_1 = P1(state);
	double	p_2 = P2(state);
	double	rho_1 = Rho1(state);
	double	rho_2 = Rho2(state);
	double	exp1, exp2;
	double	dp;

	exp1 = p_1*exp(-rho_1*V);
	exp2 = p_2*exp(-rho_2*V);
	pr = exp1 + exp2;
	prp = -rho_1*exp1 - rho_2*exp2;
	prpp = rho_1*rho_1*exp1 + rho_2*rho_2*exp2;
	dp = (W+1.0)*(p-pr);

	return 0.5*(W + 2.0)*(dp + V*V*prpp)/(dp - V*prp);
}		/*end JWL_fundamental_derivative*/

/*
*			JWL_C_V():
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*/

/*ARGSUSED*/
LOCAL	double	JWL_C_V(
	Locstate state)
{
	return R(state)/W(state);
}	/* end JWL_C_V */

/*
*			JWL_C_P():
*
*	Specific heat at constant pressure.
*
*
*                        dS  |
*		C_P = T ---- |
*                        dT  | P
*/

/*ARGSUSED*/
LOCAL	double	JWL_C_P(
	Locstate state)
{
	double	x = dpinfinity(state)/stiff_pressure(state);

	return (R(state)/W(state))*(Wp1(state) + x)/(1.0 + x);
}	/* end JWL_C_P */

/*
*			JWL_K_T():
*
*	Isothermal compressibility.
*
*                        1   dV  |
*		K_T = - --- ---- |
*                        V   dP  | T
*/

/*ARGSUSED*/
LOCAL	double	JWL_K_T(
	Locstate state)
{
	double	x = dpinfinity(state)/stiff_pressure(state);
	return (Wp1(state)+x)/((1.0+x)*adiabatic_gamma(state)*pressure(state));
}	/* end JWL_K_T */



/***************END SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS *******/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/

/***************Purely Thermodynamic Hugoniot Functions*********************/

/*
*			JWL_dens_Hugoniot():
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


LOCAL	double	JWL_dens_Hugoniot(
	double p1,
	Locstate st0)
{
	JWL_EOS	*jwleos = JWL_Eos(st0);
	double	A, B;
	double	p0 = pressure(st0), pbar;
	double	rho, rhon, rho0 = Dens(st0);
	double	rho_min, rho_max;
	double	fr, dfr;
	double	x, p_c;
	int	n;
	static const	int	max_num_iter = 10; /*TOLERANCE*/

	DEBUG_ENTER(JWL_dens_Hugoniot)
	if (DEBUG)
	{
		(void) printf("p1 = %g\n",p1);
		verbose_print_state("st0",st0);
	}
	if (p1 < 0.0)
	{
		screen("ERROR in JWL_dens_Hugoniot(), ");
		screen("negative pressure\n");
		clean_up(ERROR);
	}
	pbar = 0.5*(p1 + p0);
	x = specific_internal_energy(st0) + dH(st0) + pbar/Dens(st0);
	A = (p1/W(st0) + pbar)/x;
	B = 1.0/(W(st0)*x);

	/*
	*	Solve the equation
	*
	*	rho = A - B*pc0(rho,jwleos)
	*/

	rho_min = rho0;
	rho_max = rho0*(W(st0)+2.0)/W(st0);
	rho = A - B*pc0(rho0,jwleos);
	rho = min(rho,(1.0-EPSILON)*rho_max);
	rho = max(rho,rho_min);
	for (n = 0; n < max_num_iter; ++n)
	{
		p_c = pc0(rho,jwleos);
		fr = rho + B*p_c - A;
		if (DEBUG)
		{
			(void) printf("rho = %g, fr = %g\n",rho,fr);
			(void) printf("rho_min = %g, rho_max = %g\n",
				rho_min,rho_max);
		}
		if (fabs(fr) < rho0*EPSILON)
		{
			rhon = rho;
			break;
		}
		if (fr < 0.0)
			rho_min = max(rho,rho_min);
		else
			rho_max = min(rho,rho_max);
		dfr = 1.0 + B*p_c*rg0(rho,jwleos)/rho;
		rhon = rho - fr/dfr;
		if ((rhon < rho_min) || (rho_max < rhon))
		{
			if (DEBUG)
				(void) printf("Newton unstable\n");
			rho = 0.5*(rho_max+rho_min);
			p_c = pc0(rho,jwleos);
			fr = rho + B*p_c - A;
			if (fr < 0.0)
				rho_min = max(rho,rho_min);
			else
				rho_max = min(rho,rho_max);
			if (fabs(rho_max - rho_min) < EPSILON*rho0)
				break;
		}
		else if (fabs(rhon - rho) < EPSILON*rho0)
			break;
		rho = rhon;
	}
	if (n == max_num_iter)
	{
		screen("ERROR in JWL_dens_Hugoniot(), ");
		screen("Newton iteration did not converge\n");
		clean_up(ERROR);
	}
	if (DEBUG)
		(void) printf("rho = %g, number of iterations = %d\n",rho,n);
	DEBUG_LEAVE(JWL_dens_Hugoniot)
	return rho;
}		/*end JWL_dens_Hugoniot*/


/*
*			JWL_state_w_pr_on_Hugoniot():
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

LOCAL	void	JWL_state_w_pr_on_Hugoniot(
	Locstate st0,
	double p1,
	Locstate st1,
	int stype1)
{
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,stype1);
	Dens(st1) = dens_Hugoniot(p1,st0);
	switch(stype1)
	{
	case TGAS_STATE:
		Press(st1) = p1;
		break;
	case GAS_STATE:
		Energy(st1) = (p1 - Pc(st1))/W(st1);
		break;
	case EGAS_STATE:
		Energy(st1) = (p1 - Pc(st1))/(Dens(st1)*W(st1));
		break;
	case FGAS_STATE:
		Temperature(st1) = (p1 + pinfinity(st1))/(R(st1)*Dens(st1));
		break;
	case VGAS_STATE:
		Press(st1) = p1;
		set_type_of_state(st1,TGAS_STATE);
		set_state(st1,VGAS_STATE,st1);
		break;
	default:
		screen("ERROR in JWL_state_on_adiabat_with_pr()\n");
		screen("Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end JWL_state_w_pr_on_Hugoniot*/

/***************End Purely Thermodynamic Hugoniot Functions*****************/

/***************Purely Thermodynamic Adiabatic Wave Curve Functions*********/

/*	
*			JWL_dens_rarefaction():
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

LOCAL	double	JWL_dens_rarefaction(
	double p1,
	Locstate st0)
{
	double	k;
	double	p_1 = P1(st0), p_2 = P2(st0);
	double	rho_1 = Rho1(st0), rho_2 = Rho2(st0);
	double	p0 = pressure(st0);
	double	rho0 = Dens(st0), V0 = 1.0/rho0;
	double	rho;
	double	exp1, exp2;
	int	n;
	static const	int	max_num_iter = 10; /*TOLERANCE*/

	DEBUG_ENTER(JWL_dens_rarefaction)

	k = pow(V0,Wp1(st0))*(p0 - p_1*exp(-rho_1*V0) - p_2*exp(-rho_2*V0));

	if (DEBUG)
	{
		(void) printf("p1 = %g\n",p1);
		verbose_print_state("st0",st0);
		(void) printf("k = %g\n",k);
	}


	if (p1 > p0) /*Compression wave*/
	{
		/*
		 *	Solve the equation (for V)
		 *
		 *	 w+1
		 *	v    (p1 - p_1*exp(-rho_1*v) - p_2*exp(-rho_2*v)) = k
		 */
		double	Vn, V0 = 1.0/Dens(st0);
		double	V, Vmin, Vmax;
		double	pv, dpv;
		double	VWp1;

		if (DEBUG)
			(void) printf("Compression wave\n");

		Vmin = 0.0;
		Vmax = V0;
		V = V0;
		for (n = 0; n < max_num_iter; ++n)
		{
			VWp1 = pow(V,Wp1(st0));
			exp1 = p_1*exp(-rho_1*V);
			exp2 = p_2*exp(-rho_2*V);
			pv = VWp1*(p1 - exp1 - exp2);
			if (fabs(k - pv) <= EPSILON*k)
			{
				Vn = V;
				break;
			}
			if (k <= pv)
				Vmax = V;
			else
				Vmin = V;
			dpv = Wp1(st0)*pv/V + VWp1*(rho_1*exp1 + rho_2*exp2);
			Vn = V + (k - pv)/dpv;
			if ((Vn < Vmin) || (Vmax < Vn))
			{
				if (DEBUG)
					(void) printf("Newton unstable\n");
				V = 0.5*(Vmax+Vmin);
				VWp1 = pow(V,Wp1(st0));
				exp1 = p_1*exp(-rho_1*V);
				exp2 = p_2*exp(-rho_2*V);
				pv = VWp1*(p1 - exp1 - exp2);
				if (k <= pv)
					Vmax = V;
				else
					Vmin = V;
				Vn = V;
				if (fabs(Vmax - Vmin) < EPSILON*V0)
					break;
			}
			else if (fabs(Vn - V) < EPSILON*V0)
				break;
			V = Vn;
		}
		rho = 1.0/Vn;
	}
	else	/*rarefaction wave*/
	{
		/*
		*	Solve the equation (for r)
		*
		*     w+1
		*  k r    - (p1 - p_1*exp(-rho_1/r) - p_2*exp(-rho_2/r)) = 0
		*/
		double	rmin, rmax;
		double	r, rn;
		double	fr, dfr;
		double	krWp1;

		if (DEBUG)
			(void) printf("Rarefaction wave\n");

		rmin = 0.0;
		rmax = rho0;
		r = rho0;
		for (n = 0; n < max_num_iter; ++n)
		{
			krWp1 = k*pow(r,Wp1(st0));
			exp1 = p_1*exp(-rho_1/r);
			exp2 = p_2*exp(-rho_2/r);
			fr = krWp1 + exp1 + exp2 - p1;
			if (fabs(fr) <= EPSILON*p0)
			{
				rn = r;
				break;
			}
			if (fr <= 0.0)
				rmin = r;
			else
				rmax = r;
			dfr = Wp1(st0)*krWp1/r + (rho_1*exp1+rho_2*exp2)/(r*r);
			rn = r - fr/dfr;
			if ((rn < rmin) || (rmax < rn))
			{
				if (DEBUG)
					(void) printf("Newton unstable\n");
				r = 0.5*(rmin + rmax);
				krWp1 = k*pow(r,Wp1(st0));
				exp1 = p_1*exp(-rho_1/r);
				exp2 = p_2*exp(-rho_2/r);
				fr = krWp1 + exp1 + exp2 - p1;
				if (fr <= 0.0)
					rmin = r;
				else
					rmax = r;
				rn = r;
				if (fabs(rmax - rmin) < EPSILON*rho0)
					break;
			}
			else if (fabs(rn - r) < EPSILON*rho0)
				break;
			r = rn;
		}
		rho = rn;
	}
	if (DEBUG)
	{
		(void) printf("Number of iteratations = %d\n",n);
		(void) printf("rho = %g\n",rho);
	}
	if (n == max_num_iter)
	{
		screen("ERROR in JWL_dens_rarefaction(), ");
		screen("Newton iteration did not converge\n");
		clean_up(ERROR);
	}
	DEBUG_LEAVE(JWL_dens_rarefaction)
	return rho;
}		/*end JWL_dens_rarefaction*/

/*	
*			JWL_pressure_rarefaction():
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

LOCAL	double	JWL_pressure_rarefaction(
	double rho1,
	Locstate st0)
{
	double k;
	double	p_1 = P1(st0), p_2 = P2(st0);
	double	rho_1 = Rho1(st0), rho_2 = Rho2(st0);
	double	V = 1.0/rho1, V0 = 1.0/Dens(st0);
	double	p0 = pressure(st0);

	k = pow(V0,Wp1(st0))*(p0 - p_1*exp(-rho_1*V0) - p_2*exp(-rho_2*V0));
	return p_1*exp(-rho_1*V) + p_2*exp(-rho_2*V) + k*pow(rho_1,Wp1(st0));

}		/*end JWL_pressure_rarefaction*/


/*	
*			JWL_state_on_adiabat_with_pr():
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

LOCAL	void	JWL_state_on_adiabat_with_pr(
	Locstate st0,
	double p1,
	Locstate st1,
	int stype1)
{
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,stype1);
	Dens(st1) = dens_rarefaction(p1,st0);
	switch(stype1)
	{
	case TGAS_STATE:
		Press(st1) = p1;
		break;
	case GAS_STATE:
		Energy(st1) = (p1 - Pc(st1))/W(st1);
		break;
	case EGAS_STATE:
		Energy(st1) = (p1 - Pc(st1))/(Dens(st1)*W(st1));
		break;
	case FGAS_STATE:
		Temperature(st1) = (p1 + pinfinity(st1))/(R(st1)*Dens(st1));
		break;
	case VGAS_STATE:
		Press(st1) = p1;
		set_type_of_state(st1,TGAS_STATE);
		set_state(st1,VGAS_STATE,st1);
		break;
	default:
		screen("ERROR in JWL_state_on_adiabat_with_pr()\n");
		screen("Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end JWL_state_on_adiabat_with_pr*/

/*	
*			JWL_state_on_adiabat_with_dens():
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

LOCAL	void	JWL_state_on_adiabat_with_dens(
	Locstate st0,
	double rho1,
	Locstate st1,
	int stype1)
{
	double	p1;
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,stype1);
	p1 = pressure_rarefaction(rho1,st0);
	Dens(st1) = rho1;
	switch(stype1)
	{
	case TGAS_STATE:
		Press(st1) = p1;
		break;
	case GAS_STATE:
		Energy(st1) = (p1 - Pc(st1))/W(st1);
		break;
	case EGAS_STATE:
		Energy(st1) = (p1 - Pc(st1))/(Dens(st1)*W(st1));
		break;
	case FGAS_STATE:
		Temperature(st1) = (p1 + pinfinity(st1))/(R(st1)*Dens(st1));
		break;
	case VGAS_STATE:
		Press(st1) = p1;
		set_type_of_state(st1,TGAS_STATE);
		set_state(st1,VGAS_STATE,st1);
		break;
	default:
		screen("ERROR in JWL_state_on_adiabat_with_pr()\n");
		screen("Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end JWL_state_on_adiabat_with_dens*/

/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/

/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/

/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			JWL_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

LOCAL	void	JWL_fprint_EOS_params(
	FILE *file,
	Gas_param *params)
{
	JWL_EOS *jwleos = (JWL_EOS*) params->eos;

	(void) fprintf(file,"\tEquation of state = %d JWL\n",JWL);
	(void) fprintf(file,"\tRho0 = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &jwleos->_Rho0,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,jwleos->_Rho0);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"\tA, B = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",2);
	    (void) fwrite((const void *) &jwleos->_P1,FLOAT,2,file);
	}
	else
	    (void) fprintf(file,"%"FFMT" %"FFMT,jwleos->_P1,jwleos->_P2);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"\tR1, R2 = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",2);
	    (void) fwrite((const void *) &jwleos->_R1,FLOAT,2,file);
	}
	else
	    (void) fprintf(file,"%"FFMT" %"FFMT,jwleos->_R1,jwleos->_R2);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"\tGruneisen coefficient = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &jwleos->_W,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,jwleos->_W);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"\tR = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &jwleos->_R,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,jwleos->_R);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"\tdH = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &jwleos->_dH,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,jwleos->_dH);
	(void) fprintf(file,"\n");
}		/*end JWL_fprint_EOS_params */

/*
*			JWL_read_print_EOS_params():
*
*	Reads the equation of state data as printed by JWL_fprint_EOS_params.
*	This is restart function.
*/

/*ARGSUSED*/
LOCAL	void	JWL_read_print_EOS_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	JWL_EOS *jwleos = (JWL_EOS*) params->eos;

	jwleos->_Rho0 = fread_float("Rho0 = ",io_type);
	jwleos->_P1 = fread_float("A, B = ",io_type);
	jwleos->_P2 = fread_float(NULL,io_type);
	jwleos->_R1 = fread_float("R1, R2 = ",io_type);
	jwleos->_R2 = fread_float(NULL,io_type);
	jwleos->_W = fread_float("Gruneisen coefficient = ",io_type);
	jwleos->_R = fread_float("R = ",io_type);
	jwleos->_dH = fread_float("dH = ",io_type);
	set_JWL_coefs(jwleos);
}		/*end JWL_read_print_EOS_params*/

/*
*			JWL_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

/*ARGSUSED*/
LOCAL	void	JWL_prompt_for_EOS_params(
	INIT_DATA  *init,
	Gas_param  *params,
	const char *message1,
	const char *message2)
{
	char s[Gets_BUF_SIZE];
	JWL_EOS *jwleos = (JWL_EOS*) params->eos;
	screen("Input JWL EOS parameters for the%s gas%s\n",message1,message2);
	screen("Enter the reference density Rho0: ");
	(void) Scanf("%f\n",&jwleos->_Rho0);
	screen("Enter the two reference pressures A and B: ");
	(void) Scanf("%f %f\n",&jwleos->_P1,&jwleos->_P2);
	screen("Enter the parameters R1 and R2: ");
	(void) Scanf("%f %f\n",&jwleos->_R1,&jwleos->_R2);
	screen("Enter the Gruneisen coefficient: ");
	(void) Scanf("%f\n",&jwleos->_W);
	jwleos->_R = 1.0;
	screen("Enter the gas constant R (RT = (P+P_inf)*V) (dflt = %g): ",
		jwleos->_R);
	(void) Gets(s);
	if (s[0] != '\0')
		(void) sscan_float(s,&jwleos->_R);
	jwleos->_dH = 0.0;
	screen("Enter the heat of detonation (dflt = %g): ",jwleos->_dH);
	(void) Gets(s);
	if (s[0] != '\0')
		(void) sscan_float(s,&jwleos->_dH);

	set_JWL_coefs(jwleos);
}		/*end JWL_prompt_for_EOS_params*/
/***************END INITIALIZATION UTILITY FUNCTIONS************************/

/*
*	Auxilary support functions.
*/

LOCAL	void	set_JWL_coefs(
	JWL_EOS *jwleos)
{
	jwleos->_Rho1 = jwleos->_Rho0*jwleos->_R1;
	jwleos->_Rho2 = jwleos->_Rho0*jwleos->_R2;
	jwleos->_W_1 = jwleos->_W/jwleos->_Rho1;
	jwleos->_W_2 = jwleos->_W/jwleos->_Rho2;
	jwleos->_Wp1 = jwleos->_W + 1.0;
	jwleos->_Pref0 = jwleos->_P1*exp(-jwleos->_R1) +
			 jwleos->_P2*exp(-jwleos->_R2);
}		/*end set_JWL_coefs*/


LOCAL	double	pc0(
	double	rho,
	JWL_EOS	*jwleos)
{
	double	p_1 = jwleos->_P1;
	double	p_2 = jwleos->_P2;
	double	W_1 = jwleos->_W_1;
	double	W_2 = jwleos->_W_2;
	double	rho_1 = jwleos->_Rho1;
	double	rho_2 = jwleos->_Rho2;

	return p_1*(1.0 - W_1*rho)*exp(-rho_1/rho) +
				p_2*(1.0 - W_2*rho)*exp(-rho_2/rho);
}	/*end pc0*/

LOCAL	double	pc(
	double	rho,
	JWL_EOS	*jwleos)
{
	return	pc0(rho,jwleos) + jwleos->_W*rho*jwleos->_dH;
}	/*end pc*/

/*
*			rg0():
*
*
*	rg0 = - dlog(pc0)/dlog(v).
*
*/

LOCAL	double	rg0(
	double	rho,
	JWL_EOS	*jwleos)
{
	double	p_1 = jwleos->_P1;
	double	p_2 = jwleos->_P2;
	double	W_1 = jwleos->_W_1;
	double	W_2 = jwleos->_W_2;
	double	rho_1 = jwleos->_Rho1;
	double	rho_2 = jwleos->_Rho2;
	double	p_c, e1, e2;

	e1 = exp(-rho_1/rho);
	e2 = exp(-rho_2/rho);
	p_c = p_1*(1.0 - W_1*rho)*e1 + p_2*(1.0 - W_2*rho)*e2;
	return -(1.0/(rho*p_c))*
		(  (p_1*W_1*rho*rho - p_1*rho_1*(1.0 - W_1*rho))*e1 +
		   (p_2*W_2*rho*rho - p_2*rho_2*(1.0 - W_2*rho))*e2 );

}	/*end rg0*/


LOCAL	double	pinfinity(
	Locstate	state)
{
	return -(P1(state)*exp(-Rho1(state)/Dens(state)) +
		 P2(state)*exp(-Rho2(state)/Dens(state))) +
		pow(Dens(state)/Rho0(state),Wp1(state))*Pref0(state);
}	/*end pinfinity*/

LOCAL	double	dpinfinity(
	Locstate	state)
{
	return -(Pc0(state)*Dens(state)/W(state))*(1.0 - Rg0(state))
		- Dens(state)*Wp1(state)*pinfinity(state);
}	/*end dpinfinity*/
