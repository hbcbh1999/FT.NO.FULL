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
*				gistate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains initialization routines for gas states.
*/

#include <ginit/ginit.h>



/*	
*		init_shock_states():
*
*	Computes the state behind a shock moving into stagnant gas 
*	of density rho and pressure pr given M, the mach number of the
*	shock relative to the sound speed ahead.
*/



EXPORT void init_shock_states(
	double		rho,
	double		pr,
	double		M,
	double		*s_n,
	Gas_param	*params,
	Locstate	ahead,
	Locstate	behind)
{
	int		i, dim;
	double		shock_speed;

		/* Set ahead state */

	Init_params(ahead,params);
	dim = params->dim;
#if defined(COMBUSTION_CODE)
	prompt_for_burning(&Params(ahead)," ahead");
#endif /* defined(COMBUSTION_CODE) */
	Dens(ahead) = rho;
	Press(ahead) = pr;
	for (i = 0; i < dim; i++)
		Vel(ahead)[i] = 0.;
	set_type_of_state(ahead,TGAS_STATE);
	set_state(ahead,GAS_STATE,ahead);

		/* Set behind state */

#if defined(COMBUSTION_CODE)
	if (params->composition_type == PTFLAME)
	{
		Locstate CJ;
		size_t sizest = params->sizest;

		(*params->_alloc_state)(&CJ,sizest);
		shock_speed = CJ_det(CJ,TGAS_STATE,ahead,RIGHT_FAMILY);
		Set_params(Params(behind),Params(ahead)->other_params);
		Dens(behind) = Dens(CJ);
		Press(behind) = Press(CJ);
		for (i = 0; i < dim; i++)
			Vel(behind)[i] = Vel(CJ)[i] * s_n[i];
		set_type_of_state(behind,TGAS_STATE);
		set_state(behind,GAS_STATE,behind);
		free(CJ);
	}
	else
#endif /* defined(COMBUSTION_CODE) */
		(void) s_polar_4(SHOCK_MACH_NUMBER,M,
			&shock_speed,s_n,ahead,behind,GAS_STATE);


	if (debugging("init"))
	{
		verbose_print_state("ahead",ahead);
		verbose_print_state("behind",behind);
	}
}		/*end init_shock_states*/

EXPORT	_RAREFACTION_WAVE_1D *allocate_RAREFACTION_WAVE_1D(
	Front *front)
{
	_RAREFACTION_WAVE_1D *rw1d;
	scalar(&rw1d,sizeof(_RAREFACTION_WAVE_1D));
	alloc_state(front->interf,&rw1d->stl,front->sizest);
	alloc_state(front->interf,&rw1d->stt,front->sizest);
	return rw1d;
}		/*end allocate_RAREFACTION_WAVE_1D*/


#if defined(COMBUSTION_CODE)
EXPORT void prompt_for_burning(
	Gas_param	**params,
	const char	*message)
{
	if ((*params)->composition_type == PTFLAME ||
	    (*params)->composition_type == THINFLAME)
	{
	    char s[Gets_BUF_SIZE];

	    screen("Type 'y' if the gas%s has burned: ",message);
	    (void) Gets(s);
	    if (s[0] == 'y')
	    {
	    	if ((*params)->burned == UNBURNED)
	    		*params = (*params)->other_params;
	    }
	    else
	    {
	    	if ((*params)->burned == BURNED)
	    		*params = (*params)->other_params;
	    }
	}
}		/*end prompt_for_burning*/
#endif /* defined(COMBUSTION_CODE) */

/*
*		prompt_for_behind_contact_state():
*
*	Prompts for a state (behind) that is connected to the ahead
*	state by a contact with normal nor.
*/

/*ARGSUSED*/
EXPORT void prompt_for_behind_contact_state(
	Locstate	ahead,
	Locstate	behind,
	Gas_param	*params,
	boolean            allow_arbitrary_behind_state,
	double		*nor,
	int		state_type,
	INIT_DATA	*init)
{
	char 	s[Gets_BUF_SIZE];
	double 	dens;
	double	shear[MAXD];
	int 	i, dim = params->dim;
#if defined(TWOD) || defined(THREED)
	double	shear_nor;
#endif /* defined(TWOD) || defined(THREED) */
 
	if (debugging("init_prompt"))
		verbose_print_state("ahead contact state",ahead);

	if (allow_arbitrary_behind_state)
	{
	    screen("Type 'y' to input an independent behind state: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	        prompt_for_ref_state(" behind the contact",behind,
			             TGAS_STATE,params,init);
	        return;
	    }
	}
	screen("\nIn addition to the ahead state, two more parameters\n"
	       "\tare needed to specify the contact configuration.\n");

	screen("Enter the density behind the contact: ");
	(void) Scanf("%f\n",&dens);
	for (i = 0; i < dim; i++)
	    shear[i] = 0.0;
	switch (dim)
	{
#if defined(ONED)
	case 1:
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    shear_nor = 0.0;
	    screen("Enter the velocity jump (shear) across "
	           "the contact (dflt = %g): ",shear_nor);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&shear_nor);
	    	shear[0] = -nor[1]*shear_nor;
	    	shear[1] =  nor[0]*shear_nor;
	    }
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    screen("You can add a optional velocity shear across the contact\n"
	           "\tThis can be any vector orthogonal to the "
	           "normal vector <%g, %g, %g>\n",nor[0],nor[1],nor[2]);
	    screen("Enter the the velocity jump (shear) across "
	           "the contact (dflt = (%g, %g, %g)): ",
		   shear[0],shear[1],shear[2]);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	int   n;
	    	double tmp[3];
	    	static const char *fmt = "%lf %lf %lf";

	    	n = sscanf(s,fmt,tmp,tmp+1,tmp+2);
	    	for (i = 0; i < n; i++)
	    	    shear[i] = tmp[i];
	    }
	    shear_nor = scalar_product(shear,nor,dim);
	    if (shear_nor > EPSILON)
	    {
	    	(void) printf("WARNING in prompt_for_behind_contact_state(),"
	    	              " velocity shear not orthogonal to surface\n"
	    	              "Canceling normal component of shear\n");
	    }
	    for (i = 0; i < dim; i++)
	    	shear[i] -= shear_nor*nor[i];
	    break;
#endif /* defined(THREED) */
	}

	Init_params(behind,params);
	Dens(behind) = dens;
	Press(behind) = pressure(ahead);
	for (i = 0; i < dim; i++)
	    Vel(behind)[i] = vel(i,ahead) + shear[i];
	set_type_of_state(behind,TGAS_STATE);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            Gas_param *params = Params(ahead);
            if(params->n_comps != 1)
                return;
        }
	set_state(behind,state_type,behind);
	if (debugging("init_prompt"))
	    verbose_print_state("behind contact state",behind);
}		/*end prompt_for_behind_contact_state*/


EXPORT 	void 	prompt_for_behind_shock_state(
	Locstate	ahead,
	Locstate	behind,
	boolean            allow_arbitrary_behind_state,
	double		*nor,
	int		state_type,
	boolean		isforward,
	INIT_DATA	*init)
{
	char 		s[Gets_BUF_SIZE];
	double 		parameter, shock_speed, nor_b2a[MAXD];
	int		sign = isforward ? 1 : -1;
	int 		which_parameter, i;
	size_t 		sizest = Params(ahead)->sizest;
	int 		dim = Params(ahead)->dim;
	Locstate 	gahead;

	debug_print("behind_shock","Entered prompt_for_behind_shock_state()\n");

	if (allow_arbitrary_behind_state)
	{
	    screen("Type 'y' to input an independent behind state: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	        prompt_for_ref_state(" across the shock",behind,
			             TGAS_STATE,Params(ahead),init);
	        return;
	    }
	}
	(*Params(ahead)->_alloc_state)(&gahead,sizest);
 
	if (debugging("behind_shock"))
	    verbose_print_state("ahead shock state",ahead);
	screen("\nIn addition to the ahead state, one more parameter\n"
	       "is needed to specify the shock wave configuration.\n"
	       "The choices are\n"
	       "\tThe pressure behind the shock (P)\n"
	       "\tThe magnitude of the normal component of the gas velocity "
	           "behind the shock (V)\n"
	       "\tThe normal component of the shock speed (S)\n"
	       "\tThe normal shock Mach number (M)\n"
	       "Enter choice here: ");
	(void) Gets(s);

	if (s[0] == 'p' || s[0] == 'P')
	{
	    which_parameter = BEHIND_PRESSURE;
	    screen("Enter the pressure behind the shock: ");
	    (void) Scanf("%f\n",&parameter);
	}
	else if (s[0] == 'v' || s[0] == 'V')
	{
	    which_parameter = BEHIND_VELOCITY;
	    screen("Enter the magnitude of the normal component of the gas "
		   "velocity behind the shock: ");
	    (void) Scanf("%f\n",&parameter);
	}
	else if (s[0] == 's' || s[0] == 'S')
	{
	    which_parameter = SHOCK_SPEED;
	    screen("Enter the normal component of the speed of the shock: ");
	    (void) Scanf("%f\n",&parameter);
	}
	else if (s[0] == 'm' || s[0] == 'M')
	{
	    which_parameter = SHOCK_MACH_NUMBER;
	    screen("Enter the normal shock Mach number of the wave: ");
	    (void) Scanf("%f\n",&parameter);
	}
	else
	{
	    screen("No such choice\n");
	    clean_up(ERROR);
	}

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* some options in s_polar_4() involves
             * computing of sound speed and so on,
             * which needs partial density.
             */
            Gas_param *params = Params(ahead);
            if(params->n_comps != 1)
            {
                set_state(gahead, TGAS_STATE, ahead);
                pdens(gahead)[0] = Dens(gahead);
                set_state(gahead, GAS_STATE, gahead);
            }
            else
                set_state(gahead, GAS_STATE, ahead);
        }
        else
	    set_state(gahead, GAS_STATE, ahead);
	for ( i = 0; i < dim; i++ )
	    nor_b2a[i] = sign*nor[i];
	(void) s_polar_4(which_parameter,parameter,&shock_speed,
			 nor_b2a,gahead,behind,state_type);

	if ( entropy(ahead) - entropy(behind) > 0.0 )
	{
	    screen("ERROR in prompt_for_behind_shock_state(), "
		   "Entropy of the behind state is smaller than "
		   "that of the ahead state!");
	    clean_up(ERROR);
	}

	if (debugging("behind_shock"))
	{
	    double	abs_v[MAXD];

	    for (i = 0; i < dim; i++)
	    	abs_v[i] = shock_speed*nor_b2a[i];
	    verbose_print_state("behind shock state",behind);
	    (void) printf("Shock Mach number = %g\n",
	    	          mach_number(ahead,abs_v));
	    (void) printf("Shock speed = %g\n",shock_speed);
	}

	free(gahead);
	debug_print("behind_shock","Left prompt_for_behind_shock_state()\n");
}		/*end prompt_for_behind_shock_state*/

EXPORT	int	init_state_type(
	int	stype)
{
	char		s[Gets_BUF_SIZE];

	screen("Enter the desired state type representation, choices are\n");
	screen("\tDensity, total energy, momentum, (GAS_STATE)\n");
	screen("\tDensity, pressure, velocity, (TGAS_STATE)\n");
	screen("\tDensity, specific internal energy, velocity, (EGAS_STATE)\n");
	screen("\tDensity, Temperature, velocity, (FGAS_STATE)\n");
	screen("Enter choice (dflt = %s): ",state_type_name(stype));
	(void) Gets(s);
	if (s[0] != '\0')
	    stype = g_read_state_type_from_string(s);
	return stype;
}		/*end init_state_type*/

