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
*				giparams.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains high level drivers for the initialization of printout
*	of gas params structures.
*/

#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	boolean	prompt_for_artificial_parameter(double*,double,int,const char*,
						const char*,const char*);
LOCAL	void	prompt_for_muscl_artificial_parameters(AVISC*,int,
						       const char*,const char*);
LOCAL	void	prompt_for_thermodynamic_restrictions(Gas_param*);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
LOCAL	void	prompt_for_thermodynamic_restrictions_values(Gas_param*,
							     const char*);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

#if defined(COMBUSTION_CODE)
LOCAL   double   timmes_woosley_flame_velocity(Locstate);
LOCAL   double   experimental_flame_velocity(Locstate);
LOCAL   double   dynamic_flame_velocity(Locstate);
#endif /* defined(COMBUSTION_CODE) */


/*	
*			g_prompt_for_eos_params():
*/

EXPORT  Gas_param       *g_prompt_for_eos_params(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	boolean		prompt_for_visc,
	const char	*mesg)
{
	static  int             num_eos = 0;
	static  Gas_param       **paramlist = NULL;
	int                     j = 0;

	if ((paramlist == NULL) && (num_eos == 0))
	{
	    paramlist = prompt_for_eos_params_list(init,ip,prompt_for_visc,
	                			   &num_eos);
	    return NULL;
	}
	else if (paramlist == NULL)
	{
	    screen("ERROR: in g_prompt_for_eos_params()\n"
	           "paramlist is null.\n");
	    clean_up(ERROR);
	}

 	if (num_eos > 1)
	{
	    char	s[2048];
	    char	prompt[256];

	    (void) sprintf(prompt,"%s%d%s%s: ",
	   	           "\nInput the EOS model (0 <= an integer <= ",
	                   num_eos-1,", p prints available options) ",mesg);
	    screen("\n\t%s",prompt);
	    (void) fgets(s,2046,stdin);
	    while (s[0] == 'p' || s[0] == 'P')
	    {
	        IMPORT	boolean suppress_prompts;
	        if (suppress_prompts == NO)
	        {
	            for (j = 0; j < num_eos; ++j)
	            {
	                (void) fprintf(stderr,
	                	       "Equation of state number %d\n",j);
	                fprint_Gas_param(stderr,paramlist[j]);
	            }
	            screen("%s",prompt);
	        }
	        (void) fgets(s,2046,stdin);
	    }
	    (void) fprintf(stdout,"%s\n",s);
	    (void) sscanf(s,"%d\n",&j);
	    check_int_input("for the EOS model",j,0,num_eos-1,GE_AND_LE);
	}
	return  paramlist[j];
}		/*end g_prompt_for_eos_params*/


/*
*		g_prompt_for_eos_params_list():
*
*	Initializes an array of Gas_param structures.
*/

EXPORT	Gas_param **g_prompt_for_eos_params_list(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	boolean		prompt_for_visc,
	int		*num_eos)
{
	Gas_param **prms;
	int i;
	char mesg[80];

	screen("\nYou will now be prompted for the number\n"
	       "\tof different equations of state models,  "
	       "followed by prompts\n"
	       "\tfor the parameters of each EOS.  The various equations\n"
	       "\tof state will then be referred to by the integer\n"
	       "\tthat corresponds to the order in which they are prompted.\n"
	       "Enter the number of EOS models to be used: ");
	(void) Scanf("%d\n",num_eos);
	check_int_input("for the number of EOS models", *num_eos, 1, 1, GE_);

	prms = NULL;
	uni_array(&prms,*num_eos,sizeof(Gas_param *));

	for (i = 0; i < *num_eos; ++i)
	{
	    screen("\n");
	    (void) sprintf(mesg," with index %d",i);
	    prms[i] = init_eos_params(init,ip,mesg,prompt_for_visc);
	}
	return prms;
}		/*end g_prompt_for_eos_params_list*/

/*
*			g_init_eos_params():
*/

EXPORT Gas_param *g_init_eos_params(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	const char	*message,
	boolean		prompt_for_visc)
{
	Gas_param	*params;
	boolean		bo;
	char		s[Gets_BUF_SIZE];

#if !defined(COMBUSTION_CODE)
	prompt_for_equation_of_state(init,&params,"",message,ip);
	if (params == NULL) 
	{
	    (void) printf("Initialized params\n");
	    print_Gas_param(params);
	    return NULL; /* Obstacle EOS */
	}
	prompt_for_artificial_viscosity_and_heat_conduction(init,
	    "special ","for this EOS model, ",prompt_for_visc,
	    (params != NULL) ? &params->avisc : (AVISC *)NULL);
#else /* !defined(COMBUSTION_CODE) */
	if (g_composition_type() == PURE_NON_REACTIVE)
	{
	    prompt_for_equation_of_state(init,&params,"",message,ip);
	    if (params == NULL) 
	    {
	    	(void) printf("Initialized params\n");
	    	print_Gas_param(params);
	    	return NULL; /* Obstacle EOS */
	    }
	    prompt_for_artificial_viscosity_and_heat_conduction(init,
	    	"special ","for this EOS model, ",prompt_for_visc,
	        (params != NULL) ? &params->avisc : (AVISC *)NULL);
	    params->composition_type = PURE_NON_REACTIVE;
	    params->burned = UNBURNED;
	    params->q = 0.;
	    params->critical_temperature = 0.;
	    params->activ_en = 0.;
	    params->rate_mult = 0.;
	    params->speed_corr = 0.;
	    params->p_corr = 0.;
	    params->u_corr = 0.;
	    params->rho_corr = 0.;
	    params->other_params = NULL;
	}
	else
	{
	    int		composition_type = g_composition_type();
	    Gas_param	*unburned_params;
	    Gas_param	*burned_params;

	    /*
	    *	A combustible fluid has a mutually linked pair
	    *	of Gas_param structures, corresponding to the
	    * 	burned and unburned phases.  To flip from one
	    *	phase to the other, simply use
	    *
	    *		params = params->other_params
	    *
	    *	The call params = init_eos_params(...) initializes
	    *	params to correspond to the unburned phase, and
	    *	initializes params->other_params to correspond to the
	    *	burned phase.  To get the burned parameters,
	    *	flip params after the call to init_eos_params().
	    */

	    prompt_for_equation_of_state(init,&burned_params,
					 " burned",message,ip);
	    if (params == NULL) 
	    {
	    	(void) printf("Initialized params\n");
	    	print_Gas_param(params);
	    	return NULL; /* Obstacle EOS */
	    }
	    prompt_for_artificial_viscosity_and_heat_conduction(init,"special ",
	        "for this burned EOS model, ",prompt_for_visc,
	        (burned_params != NULL) ? &burned_params->avisc :
		    NULL);
	    burned_params->composition_type = composition_type;
	    burned_params->burned = BURNED;

	    prompt_for_equation_of_state(init,&unburned_params," unburned",
	                                 message,ip);
	    if (params == NULL) 
	    {
	    	(void) printf("Initialized params\n");
	    	print_Gas_param(params);
	    	return NULL; /* Obstacle EOS */
	    }
	    prompt_for_artificial_viscosity_and_heat_conduction(init,
	                "special ","for this unburned EOS model, ",
	                prompt_for_visc,
	                (unburned_params != NULL) ? &unburned_params->avisc :
	                    NULL);
	    burned_params->other_params = unburned_params;
	    unburned_params->composition_type = composition_type;
	    unburned_params->burned = UNBURNED;
	    unburned_params->burned = UNBURNED;

	    /* I don't know why */
	    burned_params->other_params = unburned_params;
	    unburned_params->other_params = burned_params;

	    screen("Enter heat released upon combustion: ");
	    (void) Scanf("%f\n",&burned_params->q);
	    unburned_params->q = burned_params->q;

	    debug_print("init","heat released%s initialized to %g\n",
	    	message,unburned_params->q);

	    screen("Enter the critical temperature for burning: ");
	    (void) Scanf("%f\n",&burned_params->critical_temperature);
	    unburned_params->critical_temperature =
	                burned_params->critical_temperature;

	    debug_print("init","critical temperature%s initialized to %g\n",
	    	message,unburned_params->critical_temperature);

	    if (unburned_params->composition_type == PTFLAME)
	    {
	    	screen("Enter curvature correction coefficients \n");
	    	screen("for det_speed, pressure, velocity, density: ");
	    	(void) Scanf("%f %f %f %f\n",&burned_params->speed_corr,
	    		     &burned_params->p_corr,&burned_params->u_corr,
	                     &burned_params->rho_corr);
	        unburned_params->speed_corr = burned_params->speed_corr;
	        unburned_params->p_corr = burned_params->p_corr;
	        unburned_params->u_corr = burned_params->u_corr;
	        unburned_params->rho_corr = burned_params->rho_corr;
	    }
	    else
	    {
	    	burned_params->speed_corr = 0.;
	    	burned_params->p_corr = 0.;
	    	burned_params->u_corr = 0.;
	    	burned_params->rho_corr = 0.;
	    	unburned_params->speed_corr = 0.;
	    	unburned_params->p_corr = 0.;
	    	unburned_params->u_corr = 0.;
	    	unburned_params->rho_corr = 0.;
	    }
	    if (debugging("init"))
	    {
	    	(void) printf("curvature corrections for speed,p,u,");
	    	(void) printf("density = %g %g %g %g\n",
	    	              burned_params->speed_corr,burned_params->p_corr,
	                      burned_params->u_corr,burned_params->rho_corr);
	    }

	    if (composition_type == ZND)
	    {
	    	screen("Enter the activation energy: ");
	    	(void) Scanf("%f\n",&burned_params->activ_en);
	    	unburned_params->activ_en = burned_params->activ_en;

	    	debug_print("init","activation energy%s initialized to %g\n",
	    		message,unburned_params->activ_en);
	    }

	    if (composition_type == THINFLAME)
            {
		screen("Enter method to calculate flame speed "
		        "w.r.t unburned fuel.\n"
		        "The choices are\n"
                        "\t Timmes-Woosley's formula (T)\n"
			"\t Experimental formula (E)\n"
			"\t Dynamic formula (D)\n"
			"\t Constant value (C)\n" 
			"Enter choice here: "); 
		(void) Gets(s);
		if((s[0] == 't') || (s[0] == 'T'))
		{
			unburned_params->_flame_velocity = 
				timmes_woosley_flame_velocity;
			unburned_params->code = 'T';
			
			unburned_params->ncoeffs = 0;
		}
		else if((s[0] == 'e') || (s[0] == 'E'))
		{	
			screen("Enter K and Q: ");
			(void) Scanf("%lf %lf\n", 
				     &unburned_params->flame_coeff[0],
				     &unburned_params->flame_coeff[1]);	
			unburned_params->_flame_velocity = 
				experimental_flame_velocity;
			unburned_params->code = 'E';
			unburned_params->ncoeffs = 2;
		}
		else if((s[0] == 'd') || (s[0] == 'D'))
		{
			unburned_params->_flame_velocity = 
				dynamic_flame_velocity;
			unburned_params->code = 'D';
			unburned_params->ncoeffs = 0;
		}
		else if((s[0] == 'c') || (s[0] == 'C'))
		{	
			screen("Enter flame speed w.r.t unburned fuel: ");
			(void) Scanf("%f\n", &unburned_params->flame_coeff[0]);
			unburned_params->code = 'C';
			unburned_params->ncoeffs = 1;
		}

		burned_params->_flame_velocity = 0;
		burned_params->code = 'C';
		burned_params->flame_coeff[0] = 0;
		burned_params->ncoeffs = 1;
	    }
	    /* flame */
	    params = unburned_params;
	}
#endif /* !defined(COMBUSTION_CODE) */
	prompt_for_thermodynamic_restrictions(params);
	bo = is_binary_output();
	set_binary_output(NO);
	(void) printf("Initialized params\n");
	print_Gas_param(params);
	set_binary_output(bo);
	screen("\n");
	return params;
}		/*end g_init_eos_params*/

/*ARGSUSED*/
LOCAL	void	prompt_for_thermodynamic_restrictions(
	Gas_param	*params)
{
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	char	s[Gets_BUF_SIZE];

	screen("Use defaults for thermodynamic restrictions (dflt = yes): ");
	(void) Gets(s);
	if ((s[0] != 'n') && (s[0] != 'N'))
	    return;

#if !defined(COMBUSTION_CODE)
	prompt_for_thermodynamic_restrictions_values(params,"");
#else /* !defined(COMBUSTION_CODE) */
	if (params->composition_type == PURE_NON_REACTIVE)
	    prompt_for_thermodynamic_restrictions_values(params,"");
	else
	{
	    prompt_for_thermodynamic_restrictions_values(params,"unburned ");
	    prompt_for_thermodynamic_restrictions_values(params->other_params,
	                				 "burned ");
	}
#endif /* !defined(COMBUSTION_CODE) */

#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
}

#if !defined(UNRESTRICTED_THERMODYNAMICS)
LOCAL	void	prompt_for_thermodynamic_restrictions_values(
	Gas_param	*params,
	const char	*mesg)
{
	char	s[Gets_BUF_SIZE];
	const double   eps = MACH_EPS;

	if (params == NULL)
	    return;

	screen("Enter the minimum allowed %spressure "
	       "(dflt = %g, mach for machine epsilon): ",mesg,
		params->min_pressure);
	(void) Gets(s);
	if (strncasecmp(s,"m",1) == 0)
	    params->min_pressure = eps;
	else if (s[0] != '\0')
	    (void) sscan_float(s,&params->min_pressure);
	screen("Enter the minimum allowed %sinternal energy "
	       "(dflt = %g, mach for machine epsilon): ",
		mesg,params->min_energy);
	(void) Gets(s);
	if (strncasecmp(s,"m",1) == 0)
	    params->min_energy = eps;
	else if (s[0] != '\0')
	    (void) sscan_float(s,&params->min_energy);
	screen("Enter the %svacuum density cutoff "
	       "(dflt = %g, mach for machine epsilon): ",mesg,
	       params->vacuum_dens);
	(void) Gets(s);
	if (strncasecmp(s,"m",1) == 0)
	    params->vacuum_dens = eps;
	else if (s[0] != '\0')
	    (void) sscan_float(s,&params->vacuum_dens);
	screen("Enter the %srarefaction pressure tolerance "
	       "(dflt = %g, mach for 1 - machine epsilon): ",
		mesg,params->raref_press);
	(void) Gets(s);
	if (strncasecmp(s,"m",1) == 0)
	    params->raref_press = 1.0 - eps;
	else if (s[0] != '\0')
	    (void) sscan_float(s,&params->raref_press);

#if defined(COMBUSTION_CODE)
	if (params->composition_type == PURE_NON_REACTIVE)
	    return;
	
	screen("Enter the %scombustion minimum rarefaction pressure ratio "
	       "(dflt = %g, mach for machine epsilon): ",
		mesg,params->tol_alpha);
	(void) Gets(s);
	if (strncasecmp(s,"m",1) == 0)
	    params->tol_alpha = eps;
	else if (s[0] != '\0')
	    (void) sscan_float(s,&params->tol_alpha);

	screen("Enter the %scombustion "
	       "rarefaction pressure tolerance "
	       "(dflt = %g, mach for machine epsilon): ",
		mesg,params->tol_press);
	(void) Gets(s);
	if (strncasecmp(s,"m",1) == 0)
	    params->tol_press = eps;
	else if (s[0] != '\0')
	    (void) sscan_float(s,&params->tol_press);
#endif /* defined(COMBUSTION_CODE) */
}		/*end prompt_for_thermodynamic_restrictions*/
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

EXPORT	void prompt_for_artificial_viscosity_and_heat_conduction(
	INIT_DATA	*init,
	const char	*mesg1,
	const char	*mesg2,
	boolean		prompt,
	AVISC		*avisc)
{
	IMPORT		boolean	suppress_prompts;
	boolean		interactive;
	char		s[Gets_BUF_SIZE];
	AVISC		Avisc;
	int		restart = (avisc == NULL) ? YES : NO;
	static boolean	first = YES;

	if (!prompt)
	    return;
	(void) printf("\n");

	interactive = (suppress_prompts) ? NO : YES;
	default_artificial_viscosity(&Avisc);
	if (init != NULL)
	    Avisc.hyp_method = hyperbolic_method_name(init);

	if (first == YES)
	{
	    boolean bin_io = is_binary_output();

	    first = NO;
	    screen("Current defaults for artificial viscosity parameters---\n");
	    set_binary_output(NO);
	    print_avisc_structure(&Avisc,interactive);
	    set_binary_output(bin_io);
	}

	screen("Use current defaults for artificial viscosity parameters\n");
	screen("\t(dflt = y,  type p to print defaults): ");
	(void) Gets(s);
	if (s[0] == 'p' || s[0] == 'P')
	{
	    boolean	bin_io = is_binary_output();

	    set_binary_output(NO);
	    print_avisc_structure(&Avisc,interactive);
	    set_binary_output(bin_io);
	    screen("Use current defaults for artificial ");
	    screen("viscosity parameters (dflt = y): ");
	    (void) Gets(s);
	}
	if (s[0] != 'n' && s[0] != 'N')
	{
	    if (avisc != NULL)
	    	*avisc = Avisc;
	    (void) printf("\n");
	    return;
	}

	if ((strstr(Avisc.hyp_method,"VECTOR_MUSCL") != NULL) ||
	   (strstr(Avisc.hyp_method,"PLM") != NULL))
	{
	    prompt_for_muscl_artificial_parameters(&Avisc,restart,mesg1,mesg2);
	}

	use_lapidus_artificial_viscosity(Avisc) =
	    prompt_for_artificial_parameter(&Avisc.lapidus_visc_coef,0.5,
	                                    restart,
	                		    "Lapidus nonlinear artificial "
	                			"viscosity",
	                                    mesg1,mesg2);
	if (use_lapidus_artificial_viscosity(Avisc))
	    Avisc.sp_coef = lapidus_stability_factor(Avisc.lapidus_visc_coef);

	use_linear_artificial_viscosity(Avisc) =
	    prompt_for_artificial_parameter(&Avisc.linear_visc_coef,0.05,
	                                    restart,
	                		    "linear artificial viscosity",
	                		    mesg1,mesg2);

	use_upwind_artificial_viscosity(Avisc) =
	    prompt_for_artificial_parameter(&Avisc.upwind_visc_coef,0.1,
	                                   restart,
	                		   "upwind artificial viscosity",
	                		   mesg1,mesg2);

	(void) prompt_for_artificial_parameter(&Avisc.heat_cond,0.05,
	                                       restart,
	                		       "artificial heat conduction",
	                		       mesg1,mesg2);

	/* see set_pjump_at_wave() for description of dynamic surface tension */
	screen("Dynamic surface tension is used to stabilize contacts\n"
	       "against shear instabilities.  "
	       "The value for this coefficient\n"
	       "should be of the same magnitude "
	       "as the wavelength (in zones)\n"
	       "of the instabilities you wish to suppress, i.e. 2-4 zones.\n");

	(void) prompt_for_artificial_parameter(&Avisc.dynamic_st,0.0,
	                                       restart,
	                		       "dynamic surface tension",
	                		       mesg1,mesg2);

        screen("The contact detection coefficient is used to detect flow\n"
               "regions where contact discontinuities occur,\n"
               "large values of this parameter will increase the chance that\n"
               "a region will be labeled as having a contact discontinuity\n");

        (void) prompt_for_artificial_parameter(&Avisc.contact_detector,0.0,
                                               restart,
                                               "contact dectection",
                                               mesg1,mesg2);

	if (avisc != NULL)
	    *avisc = Avisc;
	(void) printf("\n");
}		/*end prompt_for_artificial_viscosity_and_head_conduction*/

LOCAL	boolean prompt_for_artificial_parameter(
	double		*prm,
	double		nzdflt,
	int		restart,
	const char	*pname,
	const char	*mesg1,
	const char	*mesg2)
{
	char		s[Gets_BUF_SIZE];
	double		ftmp;

	if (mesg1 == NULL)
	    mesg1 = "";
	if (mesg2 == NULL)
	    mesg2 = "";

	screen("To have a %scoefficient of %s\n",mesg1,pname);
	if (strcmp(mesg2,"") != 0)
	    screen("\t%s\n",mesg2);
	screen("\tenter the coefficient ");
	if (restart == NO)
	{
	    screen("(default = %g, type d to use %g): ",*prm,nzdflt);
	}
	else
	    screen("(ignored for restart): ");
	(void) Gets(s);
	if (s[0] == 'd' || s[0] == 'D')
	    *prm = nzdflt;
	else if (s[0] != '\0')
	{
	    (void) sscan_float(s,&ftmp);
	    if (restart == NO)
		*prm = ftmp;
	}
	return (*prm != 0.0) ? YES : NO;
}		/*end prompt_for_artificial_parameter*/

LOCAL	void prompt_for_muscl_artificial_parameters(
	AVISC		*avisc,
	int		restart,
	const char	*mesg1,
	const char	*mesg2)
{
	static const double DEFAULT_MIN_SHOCK_JUMP   = 0.25;   /*TOLERANCE*/
	static const double DEFAULT_MIN_SP_VOL_JUMP  = 1.0e-6; /*TOLERANCE*/
	static const double DEFAULT_CHAR_SPD_CUTOFF  = 0.0;    /*TOLERANCE*/
	double	    min_shock_jump = DEFAULT_MIN_SHOCK_JUMP;
	double	    min_sp_vol_jump = DEFAULT_MIN_SP_VOL_JUMP;
	double	    char_speed_cutoff = DEFAULT_CHAR_SPD_CUTOFF;
	double       eta;
	char	    s[Gets_BUF_SIZE];


	screen("Do you wish to use slope flattening at strong waves\n"
	       "\t(dflt %s, d = global defaults)): ",
		(use_muscl_slope_flattening(*avisc)) ? "yes" : "no");
	(void) Gets(s);
	if (s[0] == 'd' || s[0] == 'D')
	    return;
	if (use_muscl_slope_flattening(*avisc))
	{
	    /* Default = YES */
	    if (s[0] == 'n' || s[0] == 'N')
	    {
	    	use_muscl_slope_flattening(*avisc) = NO;
	    	return;
	    }
	}
	else
	{
	    /* Default = NO */
	    if (s[0] != 'y' && s[0] != 'Y')
	    {
	    	use_muscl_slope_flattening(*avisc) = NO;
	    	return;
	    }
	    use_muscl_slope_flattening(*avisc) = YES;
	}


	eta = (avisc->msf_ieta > 1.0) ? 1.0/avisc->msf_ieta : 1.0;
	(void) prompt_for_artificial_parameter(&eta,eta,restart,
					       "the wave speed weight eta, "
					       "(0(max limiting) < eta < 1)",
					       mesg1,mesg2);
	if (eta > 1.0)
	    avisc->msf_ieta = 1.0;
	else if (eta > 0.0)
	    avisc->msf_ieta = 1.0/eta;
	else
	    avisc->msf_ieta = 1.0/eta;

	if (avisc->min_shock_jump != 0.0)
	    min_shock_jump = avisc->min_shock_jump;
	(void) prompt_for_artificial_parameter(&min_shock_jump,
	        			       DEFAULT_MIN_SHOCK_JUMP,restart,
					     "the minimum shock jump tolerance",
					       mesg1,mesg2);
	avisc->min_shock_jump = min_shock_jump;

	if (avisc->min_sp_vol_jump != 0.0)
	    min_sp_vol_jump = avisc->min_sp_vol_jump;
	(void) prompt_for_artificial_parameter(&min_sp_vol_jump,
	                                       DEFAULT_MIN_SP_VOL_JUMP,
		                               restart,
			       "\n\tthe minimum specific volume jump tolerance",
		                               mesg1,mesg2);
	avisc->min_sp_vol_jump = min_sp_vol_jump;
	if (avisc->char_speed_cutoff != 0.0)
	    char_speed_cutoff = avisc->char_speed_cutoff;
	(void) prompt_for_artificial_parameter(&char_speed_cutoff,
	                                       DEFAULT_CHAR_SPD_CUTOFF,restart,
				   "\n\tthe scaled characteristic speed cutoff",
		                               mesg1,mesg2);
	avisc->char_speed_cutoff = char_speed_cutoff;
}		/*end prompt_for_muscl_artificial_parameters*/

/* flame */
#if defined(COMBUSTION_CODE)
/*
*		timmes_woosley_flame_velocity():
*
*	From "The conductive propagation of nuclear flames - Timmes&Woosley
*	APJ-1992, Vol-396, pp-649-667, Equation-43.
*/

LOCAL  double   timmes_woosley_flame_velocity(
       Locstate state)
{
	double Xc = 0.5;			/* Molar fraction of carbon */
	double d = density(state);	/* Density of the unburned fuel */

	/* Ensure that it is cm/s */
	/*double v = 92e5*pow(d/(2e9), 0.805)*pow(Xc/.5, 0.889); */
	double v = 0.092*pow(d/2, 0.805)*pow(Xc*2, 0.889);
	/* Ensure that it is cm/s */

	return v;
}	/* end timmes_woosley_flame_velocity */

LOCAL double experimental_flame_velocity(
      Locstate state)
{
	double v;
	double p = pressure(state);
	double d = density(state);

	v = Params(state)->flame_coeff[0]*pow(p/d, 
			Params(state)->flame_coeff[1]);
	return v;
}	/* experimental_flame_velocity */

LOCAL double dynamic_flame_velocity(
      Locstate state)
{
	double v;
	double gr  = adiabatic_gamma(state);
	double mu2 = (gr-1)/(gr+1);
	double p   = pressure(state);
	double d   = density(state);
	double T   = p/d;
	double c   = sound_speed(state);
	double q   = Heat_release(state);

	double a = 1 + 4*q*mu2/(T*(1-mu2*mu2));
	v = c*sqrt(a - sqrt(a*a-1));
	return v;
}	/* end dynamic_flame_velocity */

EXPORT void read_print_flame_velocity_params(
	    INIT_DATA *init,
	    const IO_TYPE *io_type,
	    Gas_param *params)
{
	int i;
	FILE    *file = io_type->file;

	(void) fgetstring(file,"flame_velocity = ");
	(void) fscanf(file,"%c",&params->code);
	switch(params->code)
	{
	case 'T':
		params->_flame_velocity = timmes_woosley_flame_velocity;
	break;

	case 'E':
		params->_flame_velocity = experimental_flame_velocity;
	break;

	case 'D':
		params->_flame_velocity = dynamic_flame_velocity;
	break;

	case 'C':
		params->_flame_velocity = 0;
	break;

	default:
		screen("ERROR in read_print_flame_velocity_params(), "
		       "unable to read flame velocity params\n");
		clean_up(ERROR);
	break;
	}
	(void) fscanf(file,"%d",&params->ncoeffs);
	for(i=0; i<params->ncoeffs; i++)
	{
	    if (fgetstring(file," "))
	        params->flame_coeff[i] = fread_float(NULL,io_type);
	}
}	/* end read_print_flame_velocity_params */
#endif /* defined(COMBUSTION_CODE) */
/* flame */
