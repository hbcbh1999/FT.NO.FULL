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
*				gsesintrp.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the interpolators for Sesame hyp solutions.
*
*
*	Organization of state storage.
*
*	The state variables in the sesame hyp solution are arrays of
*	floats, corresponding to the dependent variables in the various
*	representations of the equation of states.  Depending on the
*	thermodynamic variables regarded as the independent variables,
*	these states are stored as follows.
*
*
*	DENSITY and TEMPERATURE 
*			state[0] = cold pressure = P(T=0)
*			state[1] = cold internal energy = e(T=0)
*			state[2] = reduced pressure = (P - P(cold))/(rho*T)
*			state[3] = reduced internal energy = (e - e(cold))/T
*			state[4] = specific entropy
*			state[5] = adiabatic gamma = P*c**2/rho
*			state[6] = gruneisen gamma
*			state[7] = integral dp/(c*rho)
*
*	DENSITY and ENERGY
*			state[0] = cold pressure = P(T=0)
*			state[1] = reduced pressure = (P - P(cold))/(rho*T)
*			state[2] = ses_rt_temp_from_grid(T,seos)
*			state[3] = specific entropy
*			state[4] = adiabatic gamma = P*c**2/rho
*			state[5] = gruneisen gamma
*			state[6] = integral dp/(c*rho)
*
*	DENSITY and ENTROPY
*			state[0] = cold pressure = P(T=0)
*			state[1] = cold internal energy = e(T=0)
*			state[2] = reduced pressure = (P - P(cold))/(rho*T)
*			state[3] = reduced internal energy = (e - e(cold))/T
*			state[4] = ses_rt_temp_from_grid(T,seos)
*			state[5] = adiabatic gamma = P*c**2/rho
*			state[6] = gruneisen gamma
*			state[7] = integral dp/(c*rho)
*
*	PRESSURE and ENTROPY
*			state[0] = cold internal energy = e(T=0)
*			state[1] = reduced internal energy = (e - e(cold))/T
*			state[2] = ses_rt_temp_from_grid(T,seos)
*			state[3] = density or ses_rt_grid_from_rho(rho,seos)
*			state[4] = adiabatic gamma = P*c**2/rho
*			state[5] = gruneisen gamma
*			state[6] = integral dp/(c*rho)
*
*/

#if defined(SESAME_CODE) && defined(TWOD)
#include <geos/sesame.h>

LOCAL	int	var[NUM_SES_VAR], num_vars;

EXPORT	void set_ses_intrp_flag(
	int		  var_flag,
	SESAME_TABLE_TYPE eos_type)
{
	num_vars = 0;
	switch (eos_type)
	{
	case SESAME_RHO_TEMP:
	    if (evaluate_pressure(var_flag))
	    	var[num_vars++] = 0;
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 1;
	    if (evaluate_pressure(var_flag))
	    	var[num_vars++] = 2;
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 3;
	    if (evaluate_entropy(var_flag))
	    	var[num_vars++] = 4;
	    if (evaluate_adb_gamma(var_flag))
	    	var[num_vars++] = 5;
	    if (evaluate_gru_gamma(var_flag))
	    	var[num_vars++] = 6;
	    if (evaluate_idpocr(var_flag))
	    	var[num_vars++] = 7;
		break;

	case SESAME_RHO_ENERGY:
	    if (evaluate_pressure(var_flag))
	    	var[num_vars++] = 0;
	    if (evaluate_pressure(var_flag))
	    	var[num_vars++] = 1;
	    if (evaluate_pressure(var_flag) || evaluate_temperature(var_flag))
		var[num_vars++] = 2;
	    if (evaluate_entropy(var_flag))
	    	var[num_vars++] = 3;
	    if (evaluate_adb_gamma(var_flag))
	    	var[num_vars++] = 4;
	    if (evaluate_gru_gamma(var_flag))
	    	var[num_vars++] = 5;
	    if (evaluate_idpocr(var_flag))
	    	var[num_vars++] = 6;
	    break;

	case SESAME_PRESS_ENTROPY:
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 0;
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 1;
	    if (evaluate_temperature(var_flag) || evaluate_energy(var_flag))
	    	var[num_vars++] = 2;
	    if(evaluate_density(var_flag))
	    	var[num_vars++] = 3;
	    if (evaluate_adb_gamma(var_flag))
	    	var[num_vars++] = 4;
	    if (evaluate_gru_gamma(var_flag))
	    	var[num_vars++] = 5;
	    if (evaluate_idpocr(var_flag))
	    	var[num_vars++] = 6;
	    break;

	case SESAME_VOLUME_PRESSURE:
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 0;
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 1;
	    if (evaluate_temperature(var_flag) || evaluate_energy(var_flag))
	    	var[num_vars++] = 2;
	    if(evaluate_entropy(var_flag))
	    	var[num_vars++] = 3;
	    if (evaluate_adb_gamma(var_flag))
	    	var[num_vars++] = 4;
	    if (evaluate_gru_gamma(var_flag))
	    	var[num_vars++] = 5;
	    if (evaluate_idpocr(var_flag))
	    	var[num_vars++] = 6;
	    break;

	case SESAME_RHO_ENTROPY:
	    if (evaluate_pressure(var_flag))
	    	var[num_vars++] = 0;
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 1;
	    if (evaluate_pressure(var_flag))
	    	var[num_vars++] = 2;
	    if (evaluate_energy(var_flag))
	    	var[num_vars++] = 3;
	    if (evaluate_pressure(var_flag) || evaluate_energy(var_flag) ||
	    	evaluate_temperature(var_flag))
	        var[num_vars++] = 4;
	    if (evaluate_adb_gamma(var_flag))
	    	var[num_vars++] = 5;
	    if (evaluate_gru_gamma(var_flag))
	    	var[num_vars++] = 6;
	    if (evaluate_idpocr(var_flag))
	    	var[num_vars++] = 7;
	    break;

	default:
	    screen("Unknown equation of state type "
	           "in set_ses_intrp_flag()\n");
	    clean_up(ERROR);
	}
}		/*end set_ses_intrp_flag*/

EXPORT	void set_ses_intrp_flag_all(
	SESAME_TABLE_TYPE table)
{
	int		i;

	switch (table)
	{
	case SESAME_RHO_TEMP:
	case SESAME_RHO_ENTROPY:
	    num_vars = 8;
	    break;

	case SESAME_RHO_ENERGY:
	case SESAME_PRESS_ENTROPY:
	case SESAME_VOLUME_PRESSURE:
	    num_vars = 7;
	    break;
	}
	for (i = 0; i < num_vars; i++)
	    var[i] = i;
}		/*end set_ses_intrp_flag_all*/

/*
*			ses_lin_comb_states():
*
*/

/*ARGSUSED*/
EXPORT void ses_lin_comb_states(
	double		alpha,
	double		beta,
	double		*crds1,
	Locstate	s1,
	double		*crds2,
	Locstate	s2,
	RECT_GRID	*gr,
	Locstate	ans)
{
	double		*Ans, *S1, *S2;
	int		i, j;

	Ans = (double *) ans; S1 = (double *) s1; S2 = (double *) s2;
	for (i = 0; i < num_vars; i++)
	{
		j = var[i];
		Ans[j] = alpha * S1[j] + beta * S2[j];
	}

}		/*end ses_lin_comb_states*/

/*
*			ses_tri_lin_comb_states():
*
*/

/*ARGSUSED*/
EXPORT boolean ses_tri_lin_comb_states(
	double		f0,
	double		f1,
	double		f2,
	double		*crds0,
	Locstate	s0,
	double		*crds1,
	Locstate	s1,
	double		*crds2,
	Locstate	s2,
	RECT_GRID	*gr,
	Locstate	ans)
{
	double		*Ans, *S0, *S1, *S2;
	int		i, j;

	Ans = (double *) ans;
	S0 = (double *) s0; S1 = (double *) s1; S2 = (double *) s2;

	for (i = 0; i < num_vars; i++)
	{
		j = var[i];
		Ans[j] = f0 * S0[j] + f1 * S1[j] + f2 * S2[j];
	}

	return FUNCTION_SUCCEEDED;
}		/*end ses_tri_lin_comb_states*/


/*
*			ses_tri_interpolator():
*
*	Forms a linear combination of three local states.  This function
*	is generally used to calculate an interpolated state from the three
*	states at the vertices of a triangle.
*
*/


/*ARGSUSED*/
EXPORT boolean ses_tri_interpolator(
	double	       *f,
	LINEAR_ELEMENT *t,
	TRI_SOLN       *soln,
	Locstate       ans)
{
	double		f0 = f[0], f1 = f[1], f2 = f[2];
	Locstate	s0 = t->s[0], s1 = t->s[1], s2 = t->s[2];

	return ses_tri_lin_comb_states(f0,f1,f2,Coords(t->p[0]),s0,
			               Coords(t->p[1]),s1,
				       Coords(t->p[2]),s2,
				       &soln->tri_grid->comp_grid,ans);

}		/*end ses_tri_interpolator*/



/*
*			ses_quad_interpolator():
*
*	Forms a bilinear combination of four local states.  This function
*	is generally used to calculate an interpolated state from the four 
*	states at the vertices of a quadrangle.
*
*/


/*ARGSUSED*/
EXPORT void ses_quad_interpolator(
	double		*f,
	BILINEAR_ELEMENT *q,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	Locstate	s[(1<<MAXD)];
	double		*Ans, *S0, *S1, *S2, *S3;
	double		fx = f[0], fy = f[1];
	double		f0,f1,f2,f3;
	int		i, j;

	states_on_bilinear_element(s,q,soln->tri_grid);

	Ans = (double *) ans;
	S0 = (double *) s[0];	S1 = (double *) s[1];
	S2 = (double *) s[2];	S3 = (double *) s[3];

	f1 = fx * (1.0 - fy);
	f2 = fx * fy;
	f3 = fy * (1.0 - fx);
	f0 = 1.0 - f1 - f2 - f3;


	for (i = 0; i < num_vars; i++)
	{
		j = var[i];
		Ans[j] = f0*S0[j] + f1*S1[j] + f2*S2[j] + f3*S3[j];
	}
}		/*end ses_quad_interpolator*/
#endif /* defined(SESAME_CODE) && defined(TWOD) */
