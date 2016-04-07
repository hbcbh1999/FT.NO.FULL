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
*			eosplot.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if defined(TWOD)
#include <geos/geosdecs.h>
#include <ginit/ginit.h>
#include <time.h>

enum {
	HUGONIOTS =		1,		/*   type of plot	   */
	ADIABATS =		2,
	HUGONIOTS_ADIABATS =	3,
	ISOTHERMS =		4,
	ISOCHORES =		5,
	ISOBARS =		6,
	CONST_ENERGY_PLOTS =	7,
	SES_TRI =		8,
	SES_PHS =		9
};

LOCAL	size_t sizest;

	/* LOCAL Function Declarations */
LOCAL	int	prompt_for_num_points(const char*,int);
LOCAL	void	const_energy_plots(Gas_param*);
LOCAL	void    enter_state(const char*,Locstate,int,Gas_param*);
LOCAL	void	eosplot_clean_up(void);
LOCAL	void	get_thermo_data(Locstate,double*,double*,double*,double*,
				double*,double*);
LOCAL	void	hugoniots_adiabats(int,Gas_param*);
LOCAL	void	isobars(Gas_param*);
LOCAL	void	isochores(Gas_param*);
LOCAL	void	isotherms(Gas_param*);
LOCAL	void	start_up(int,char**,INIT_DATA*);

#if defined(SESAME_CODE)
LOCAL	void	ses_print_tri_soln(Gas_param*);
#if defined(PHASE_CODE)
LOCAL	void	print_phase_bound(Gas_param*);
LOCAL	void	init_phase_bdry(Locstate);
#endif /* defined(PHASE_CODE) */
#endif /* defined(SESAME_CODE) */

LOCAL char *temporary_input_file = NULL;
LOCAL const char *ffmt = "%-"FFMT" %-"FFMT" %-"FFMT" %-"FFMT" %-"FFMT"\n";
LOCAL const char *sfmt = "%-"SFMT" %-"SFMT" %-"SFMT" %-"SFMT" %-"SFMT"\n";

/*ARGSUSED*/
EXPORT	int eosp_main(
	int		argc,
	char		**argv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	IMPORT	boolean	suppress_prompts;
	CHART		*root = ip->root;
	Front		*front = root->front;
	Grid		*grid = root->grid;
	Wave		*wave = root->wave;
	char		s[Gets_BUF_SIZE];
	int		plot_type;
	Gas_param	*params;
	REMAP           Remap;
	static	double	L[3] = {0.0, 0.0, 0.0};
	static	double	U[3] = {1.0, 1.0, 1.0};
	static	int	gmax[3] = {1, 1, 1};
	static	RECT_GRID	Rgr;

	sizest = g_sizest();
	suppress_prompts = NO;
	start_up(argc,argv,init);

	set_interface_hooks(2,init);
	i_intfc(init) = make_interface(2);
	set_remap(2,IDENTITY_REMAP,&Remap);
	set_rect_grid(L,U,L,U,NOBUF,NOBUF,gmax,2,&Remap,&Rgr);
	grid->rect_grid = wave->rect_grid = front->rect_grid = &Rgr;
	screen("Request binary/non-binary output %s",
	    (is_binary_output() == YES) ? "[b(dflt),n]: " : "[b,n(dflt)]: ");
	(void) Gets(s);
	if (s[0] == 'n' || s[0] == 'N')
	    set_binary_output(NO);
	else if (s[0] == 'b' || s[0] == 'B')
	    set_binary_output(YES);
	
	params = init_eos_params(init,ip,"",NO);
        for (;;)
	{
	    screen("Enter the type of plot you would like; choices are\n");
	    screen("\t%d = Hugoniots only.\n",HUGONIOTS);
	    screen("\t%d = Adiabats only.\n",ADIABATS);
	    screen("\t%d = Hugoniots and adiabats.\n",HUGONIOTS_ADIABATS);
	    screen("\t%d = Isotherms.\n",ISOTHERMS);
	    screen("\t%d = Isochores (constant density curves).\n",ISOCHORES);
	    screen("\t%d = Isobars (constant pressure curves).\n",ISOBARS);
	    screen("\t%d = Constant spec. int. energy curves.\n",
		   CONST_ENERGY_PLOTS);
#if defined(SESAME_CODE)
	    screen("\t%d = Contour plots.\n",SES_TRI);
#if defined(PHASE_CODE)
	    screen("\t%d = Phase boundaries.\n",SES_PHS);
#endif /* defined(PHASE_CODE) */
#endif /* defined(SESAME_CODE) */
	    screen("\tDefault = quit\n");
	    screen("Enter choice here: ");
	    (void) Gets(s);
	    if (s[0] == '\0' || s[0] == 'q' || s[0] == 'Q')
		break;
	    (void) sscanf(s,"%d\n", &plot_type);
		
	    switch (plot_type)
	    {
	    case HUGONIOTS_ADIABATS:
	    case HUGONIOTS:
	    case ADIABATS:
	    	hugoniots_adiabats(plot_type,params);
	    	break;
	    case ISOTHERMS:
		isotherms(params);
		break;
	    case ISOCHORES:
	    	isochores(params);
	    	break;
	    case ISOBARS:
	    	isobars(params);
	    	break;
	    case CONST_ENERGY_PLOTS:
	    	const_energy_plots(params);
	    	break;
#if defined(SESAME_CODE)
	    case SES_TRI:
	    	ses_print_tri_soln(params);
	    	break;
#if defined(PHASE_CODE)
	    case SES_PHS:
	    	print_phase_bound(params);
	    	break;
#endif /* defined(PHASE_CODE) */
#endif /* defined(SESAME_CODE) */
	    default:
	    	screen("Unknown plot type\n");
	    	break;
	    }
	    screen("Do you wish to produce "
	           "another type of plot? (default=no or quit): ");
	    (void) Gets(s);
	    if (s[0] != 'y' && s[0] != 'Y')
		break;
	}
	clean_up(0);
	return 0;
}		/*end eosp_main*/

LOCAL void hugoniots_adiabats(
	int		plot_type,
	Gas_param	*params)
{
	boolean		keep_prompting = YES;
	char		s[Gets_BUF_SIZE];
	int		number_hug, number_adbt;
	int		i, logs, var, MORE = 0;
	double		rho0, p0, e0, T0, rho1, p1, e1, T1, dpr, drho, pmax;
	double		rhomin, pmin, csq, csq0, S;
	double		log_pmax, log_p1, log_rho0, log_rho1;
	int             stype;
	static Locstate state = NULL, state1 = NULL;

	if (state == NULL)
	{
	    g_alloc_state(&state,sizest);
	    g_alloc_state(&state1,sizest);
	}
	stype = TGAS_STATE;
	stype = init_state_type(stype);
	while (keep_prompting == YES)
	{
	    if (!MORE)
	    {
	        screen("Do you wish to use log scales (default = no): ");
		(void) Gets(s);
		logs = (s[0] == 'y' || s[0] == 'Y') ? YES : NO;
	    }
	    enter_state("initial",state,stype,params);
	    verbose_print_state("Initial state",state);
	    (void) printf("Raw state data\n");
	    fprint_raw_gas_data(stdout,state,Params(state)->dim);
	    get_thermo_data(state,&rho0,&p0,&T0,&e0,&csq0,&S);
	    if (plot_type != ADIABATS)
	    {
	        if (!MORE)
		    number_hug = prompt_for_num_points("Hugoniot",100);
	        screen("Enter the pressure on the other side"
	               " of the shock ( > p0 = %g): ",p0);
	        (void) Scanf("%f\n", &pmax);
	        if (plot_type != HUGONIOTS)
	        {
	            if (!MORE)
		        number_adbt = prompt_for_num_points("adiabatic",100);
	            screen("Enter the pressure on the "
	                   "other side of the adiabat (< p0 = %g): ",p0);
	            (void) Scanf("%f\n", &pmin);
		    var = -1;
	        }
	        if (logs != 1)
	        {
	            dpr = (pmax - p0)/(number_hug - 1.0);
		    (void) output();
	            (void) printf(sfmt,"DENSITY","TEMPERATURE","PRESSURE",
				  "INT_ENERGY","SOUND_SP_SQR");
	            for (p1 = pmax; p1 >= p0 + 0.5*dpr; p1 -= dpr)
	            {
			state_w_pr_on_Hugoniot(state,p1,state1,
					       state_type(state));
			get_thermo_data(state1,&rho1,&p1,&T1,&e1,&csq,&S);
			(void) printf(ffmt,rho1,T1,p1,e1,sqrt(csq));
	            }
	            (void) printf(ffmt,rho0,T0,p0,e0,csq0);
	        }
	        else
	        {
	            log_pmax = plog(pmax);
	            dpr = (log_pmax - plog(p0))/(number_hug - 1.0);
		    (void) output();
	            (void) printf(sfmt,
	                          "LOG_DENSITY","LOG_PRESSURE","LOG_INT_ENERGY",
	                          "LOG_TEMPERATURE","SOUND_SP_SQR");
#if defined(PHASE_CODE)
		    init_phase_bdry(state);
#endif /* defined(PHASE_CODE) */
	            log_p1 = log_pmax;
	            for (i = 0; i < number_hug; ++i)
	            {
	    	        p1 = pexp(log_p1);
	    	        state_w_pr_on_Hugoniot(state,p1,state1,
					       state_type(state));
	    	        get_thermo_data(state1,&rho1,&p1,&T1,&e1,&csq,&S);
	    	        (void) printf(ffmt,
	    		              plog(rho1),log_p1,plog(e1+1.0),
	    		              plog(T1),csq);
	    	        log_p1 -= dpr;
	            }
	            (void) printf(sfmt,plog(rho0),plog(p0),plog(e0+1.0),
	                          plog(T0),csq0);
	        }
	    }
	    if (plot_type != HUGONIOTS)
	    {
	        if (plot_type == ADIABATS && !MORE)
	        {    
		    number_adbt = prompt_for_num_points("adiabatic",100);
	            screen("Enter 0 to give the density on the"
	                   " other side of the adiabat, 1 to give"
	                   " pressure: ");
	            (void) Scanf("%d\n", &var);
	        }
	        if ((var == 1) || (var == -1))
	        {
	            if ((plot_type == ADIABATS) ||
			(plot_type == HUGONIOTS_ADIABATS))
	            {    
		        if (var == 1)
			{
	                    screen("Enter the pressure on the "
	                           "other side of the adiabat: ");
	                    (void) Scanf("%f\n", &pmin);
			}
		        (void) output();
	                (void) printf(sfmt,
	                              "DENSITY","TEMPERATURE","PRESSURE",
	                              "INT_ENERGY","SOUND_SP_SQ");
	                (void) printf(ffmt,rho0,T0,p0,e0,csq0);
	            }
		    if (pmin <= p0)
		    {
	                dpr = (p0 - pmin)/(number_adbt - 1.0);
	                for (p1 = p0 - dpr; p1 >= pmin - 0.5*dpr; p1 -= dpr)
	                {
			    state_on_adiabat_with_pr(state,p1,state1,
						     state_type(state));
			    get_thermo_data(state1,&rho1,&p1,&T1,&e1,&csq,&S);
			    (void) printf(ffmt,rho1,T1,p1,e1,csq);
	                }
		    }
		    else
		    {
			pmax = pmin;
	                dpr = (pmax - p0)/(number_adbt - 1.0);
	                for (p1 = pmax; p1 >= p0 + 0.5*dpr; p1 -= dpr)
	                {
			    state_on_adiabat_with_pr(state,p1,state1,
						     state_type(state));
			    get_thermo_data(state1,&rho1,&p1,&T1,&e1,&csq,&S);
			    (void) printf(ffmt,rho1,T1,p1,e1,csq);
	                }
		    }
	        }
	        else
	        {
	            if (plot_type == ADIABATS)
	            {    
	                screen("Enter the density on the other side of the "
			       "adiabat ( < rho0 = %g): ",rho0);
	                (void) Scanf("%f\n",&rhomin);
	                if (logs != 1)
	                {
			    (void) output();
	                    (void) printf(sfmt,
	                   		  "DENSITY","TEMPERATURE","PRESSURE",
					  "INT_ENERGY","SOUND_SP_SQ");
	                    (void) printf(ffmt,rho0,T0,p0,e0,csq0);
	                }
	                else
	                {
			    (void) output();
	                    (void) printf(sfmt,
	                   		  "LOG_DENSITY","LOG_PRESSURE",
	                   		  "LOG_INT_ENERGY","LOG_TEMPERATURE",
	                		  "SOUND_SP_SQ");
	                    (void) printf(ffmt,plog(rho0),plog(p0),plog(e0),
	                    		  plog(T0),csq0);
	                }
	            }
	            if (logs != 1)
	            {
	                drho = (rho0 - rhomin)/(number_adbt - 1.0);
			for (rho1 = rho0 - drho; 
			    rho1 >= rhomin - 0.5*drho; 
			    rho1 = rho1 - drho)
			{
			    state_on_adiabat_with_dens(state,rho1,state1,
						       state_type(state));
			    get_thermo_data(state1,&rho1,&p1,&T1,&e1,&csq,&S);
			    (void) printf(ffmt,rho1,T1,p1,e1,csq);
	                }
	            }
	            else
	            {
			log_rho0 = plog(rho0);
			drho = (log_rho0 - plog(rhomin))/(number_adbt - 1.0);
			log_rho1 = log_rho0;
			for (i = 0; i < number_adbt; ++i)
			{
			    log_rho1 -= drho;
			    rho1 = pexp(log_rho1);
			    state_on_adiabat_with_dens(state,rho1,state1,
						       state_type(state));
			    get_thermo_data(state1,&rho1,&p1,&T1,&e1,&csq,&S);
			    (void) printf(ffmt,plog(rho1),plog(p1),
				          plog(e1+1.0),plog(T1),csq);
	                }
	            }
	        }
	    }
	    (void) printf("\n");
	    screen("Enter 1 for another plot, 0 or default to quit: ");
	    (void) Gets(s);
	    MORE = 0;
	    if (s[0] != '\0')
	        (void) sscanf(s,"%d", &MORE);
	    if (MORE != 1)
		break;
	}
}		/*end hugoniots_adiabats*/

LOCAL	int prompt_for_num_points(
	const char	*mesg,
	int		dflt)
{
	char		s[Gets_BUF_SIZE];

	if (mesg == NULL)
	    mesg = " ";
	screen("Enter the number of points on the %s"
	       "\n\tcurve to be plotted (default = %d): ",mesg,dflt);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d\n", &dflt);
	if (dflt <= 1)
	    dflt = 2;
	return dflt;
}		/*end prompt_for_num_points*/

LOCAL void isotherms(
	Gas_param	*params)
{
	boolean		keep_prompting = YES;
	int		i, num_points, first = YES, logs = NO;
	double		rho, p, e, T, drho, rhomin, rhomax, csq, S;
	double		log_rho, log_p, log_e, dlog_rho, log_rhomin, log_rhomax;
	char		s[Gets_BUF_SIZE];
	static Locstate state = NULL;

	if (state == NULL)
	{
	    g_alloc_state(&state,sizest);
	}

	Init_params(state,params);
	set_type_of_state(state,FGAS_STATE);
	while (keep_prompting == YES)
	{
	    screen("Enter the constant temperature, "
		   "negative or default stops: ");
	    (void) Gets(s);
	    T = -HUGE_VAL;
	    if (s[0] != '\0')
	        (void) sscanf(s,"%lf", &T);
	    if (T < 0.0)
		break;
	    Temperature(state) = T;
	    if (first)
	    {
	    	screen("Do you want log scales?: ");
	    	(void) Gets(s);
	    	if (s[0] == 'y' || s[0] == 'Y')
		    logs = YES;
	    	screen("Enter the minimum and maximum densities: ");
	    	(void) Scanf("%f %f\n", &rhomin, &rhomax);
	    	screen("Enter the number of points to be plotted: ");
	    	(void) Scanf("%d\n", &num_points);
	    }
	    if (logs)
	    {
	    	(void) output();
	    	(void) printf(sfmt,"LOG_DENSITY", "LOG_PRESSURE","LOG_INT_ENGY",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	log_rhomax = plog(rhomax);
	    	log_rhomin = plog(rhomin);
	    	log_rho = log_rhomin;
	    	dlog_rho = (log_rhomax - log_rhomin)/(num_points - 1.0);
	    }
	    else
	    {
	    	(void) output();
	    	(void) printf(sfmt,
			      "DENSITY", "PRESSURE","INT_ENERGY",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	drho = (rhomax - rhomin)/(num_points - 1.0);
	    }
	    rho = rhomin;
	    for (i = 1; i <= num_points; ++i)
	    {
	    	Dens(state) = rho;
	    	get_thermo_data(state,&rho,&p,&T,&e,&csq,&S);
	    	if (logs)
	    	{
	    	    log_p = plog(p);
	    	    log_e = plog(e + 1.0);
	    	    (void) printf(ffmt,log_rho,log_p,log_e,csq,S);
	    	    log_rho += dlog_rho;
	    	    rho = pexp(log_rho);
	    	}
	    	else
	    	{
	    	    (void) printf(ffmt,rho,p,e,csq,S);
	    	    rho += drho;
	    	}
	    }
	    (void) printf("\n");
	    first = NO;
	}
}		/*end isotherms*/

LOCAL void isochores(
	Gas_param	*params)
{
	boolean		keep_prompting = YES;
	double		rho, p, e, T, dvar;
	double		indep_var, indep_min, indep_max, csq, S;
	int		i, num_points, var, first = YES, logs = NO;
	double		log_T, log_p, log_e, dlog_var;
	double		log_indmin, log_indmax, log_var;
	char		s[Gets_BUF_SIZE];
	static Locstate state = NULL;

	if (state == NULL)
	{
	    g_alloc_state(&state,sizest);
	}

	Init_params(state,params);
	while (keep_prompting == YES)
	{
	    screen("Enter the constant density, negative or default stops: ");
	    rho = -HUGE_VAL;
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,"%lf", &rho);
	    if (rho < 0.0)
		break;
	    Dens(state) = rho;
	    if (first)
	    {
	    	screen("Do you want log scales?: ");
	    	(void) Gets(s);
	    	if (s[0] == 'y' || s[0] == 'Y')
		    logs = YES;
	    	screen("Choose the independent variable, choices are:\n"
	    	       "\t1 = pressure, 2 = Temperature, "
		       "3 = sp. int. energy\n"
	    	       "Enter choice here: ");
	    	(void) Scanf("%d\n", &var);
	    	screen("Enter the minimum and maximum values "
	    	       "of the indep. variable: ");
	    	(void) Scanf("%f %f\n", &indep_min, &indep_max);
	    	screen("Enter the number of points to be plotted: ");
	    	(void) Scanf("%d\n", &num_points);
	    }
	    if (logs)
	    {
	    	(void) output();
	    	(void) printf(sfmt,"LOG_TEMP","LOG_PRESSURE","LOG_INT_ENGY",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	log_indmax = plog(indep_max);
	    	log_indmin = plog(indep_min);
	    	log_var = log_indmin;
	    	dlog_var = (log_indmax - log_indmin)/(num_points - 1.0);
	    }
	    else
	    {
	    	(void) output();
	    	(void) printf(sfmt,
	    		      "TEMPERATURE","PRESSURE","INT_ENERGY",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	dvar = (indep_max - indep_min)/(num_points - 1.0);
	    }
	    indep_var = indep_min;

	    for (i = 1; i <= num_points; ++i)
	    {
	    	switch (var)
	    	{
	    	case 1:       /*     pressure    */
	    	    Press(state) = indep_var;
	    	    set_type_of_state(state,TGAS_STATE);
	    	    get_thermo_data(state,&rho,&p,&T,&e,&csq,&S);
	    	    if (logs)
	    	    {
	    	    	log_T = plog(T);
	    	    	log_p = log_var;
	    	    	log_e = plog(e);
	    	    	(void) printf(ffmt,log_T,log_p,log_e,csq,S);
	    	    	log_var += dlog_var;
	    	    	indep_var = pexp(log_var);
	    	    }
	    	    else
	    	    {
	    	    	(void) printf(ffmt,T,p,e,csq,S);
	    	    	indep_var += dvar;
	    	    }
	    	    break;
	    	case 2:       /*     temperature    */
	    	    Temperature(state) = indep_var;
	    	    set_type_of_state(state,FGAS_STATE);
	    	    get_thermo_data(state,&rho,&p,&T,&e,&csq,&S);
	    	    if (logs)
	    	    {
	    	    	log_T = log_var;
	    	    	log_p = plog(p);
	    	    	log_e = plog(e);
	    	    	(void) printf("%-"FFMT" %-"FFMT" %-"FFMT" ",
				      log_T,log_p,log_e);
	    		(void) printf("%-"FFMT" %-"FFMT" \n",
				      csq,S);
	    		log_var += dlog_var;
	    		indep_var = pexp(log_var);
	    	    }
	    	    else
	    	    {
	    	    	(void) printf("%-"FFMT" %-"FFMT" %-"FFMT" ",T,p,e);
	    	    	(void) printf("%-"FFMT" %-"FFMT" \n",
				      csq,S);
	    		indep_var += dvar;
	    	    }
	    	    break;
	    	case 3:       /*     sp. int. energy    */
	    	    Energy(state) = indep_var;
	    	    set_type_of_state(state,EGAS_STATE);
	    	    get_thermo_data(state,&rho,&p,&T,&e,&csq,&S);
	    	    if (logs)
	    	    {
	    	    	log_T = plog(T);
	    	    	log_p = plog(p);
	    	    	log_e = log_var;
	    	    	(void) printf("%-"FFMT" %-"FFMT" %-"FFMT" ",
				      log_T,log_p,log_e);
	    	    	(void) printf("%-"FFMT" %-"FFMT" %-"FFMT"\n",
				      csq,S);
	    	    	log_var += dlog_var;
	    	    	indep_var = pexp(log_var);
	    	    }
	    	    else
	    	    {
	    	    	(void) printf("%-"FFMT" %-"FFMT" %-"FFMT" ",
				      T,p,e);
	    	    	(void) printf("%-"FFMT" %-"FFMT" \n",
				      csq,S);
	    		indep_var += dvar;
	    	    }
	    	    break;
	    	}
	    }
	    (void) printf("\n");
	    first = NO;
	}
}		/*end isochores*/

LOCAL void isobars(
	Gas_param	*params)
{
	boolean		keep_prompting = YES;
	int		i, num_points, first = YES, logs = NO;
	double		rho, p, e, T, drho, rhomin, rhomax, csq, S;
	double		log_rho, log_T, log_e, dlog_rho, log_rhomin, log_rhomax;
	char		s[Gets_BUF_SIZE];

	static Locstate state = NULL;

	if (state == NULL)
	{
	    g_alloc_state(&state,sizest);
	}

	Init_params(state,params);
	set_type_of_state(state,TGAS_STATE);

	while (keep_prompting == YES)
	{
	    screen("Enter the constant pressure, negative or default stops: ");
	    p = -HUGE_VAL;
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,"%lf",&p);
	    if (p < 0.0)
		break;
	    Press(state) = p;
	    if (first)
	    {
	    	screen("Do you want log scales?: ");
	    	(void) Gets(s);
	    	if (s[0] == 'y' || s[0] == 'Y')
		    logs = YES;
	    	screen("Enter the minimum and maximum densities: ");
	    	(void) Scanf("%f %f\n", &rhomin, &rhomax);
	    	screen("Enter the number of points to be plotted: ");
	    	(void) Scanf("%d\n", &num_points);
	    }
	    if (logs)
	    {
	    	(void) output();
	    	(void) printf(sfmt,"LOG_DENSITY","LOG_TEMP"," LOG_INT_ENGY",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	log_rhomax = plog(rhomax);
	    	log_rhomin = plog(rhomin);
	    	log_rho = log_rhomin;
	    	dlog_rho = (log_rhomax - log_rhomin)/(num_points - 1.0);
	    }
	    else
	    {
	    	(void) output();
	    	(void) printf(sfmt,
	    		      "DENSITY","TEMPERATURE","INT_ENERGY",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	drho = (rhomax - rhomin)/(num_points - 1.0);
	    }
	    rho = rhomin;
	    for (i = 1; i <= num_points; ++i)
	    {
	    	Dens(state) = rho;
	    	get_thermo_data(state,&rho,&p,&T,&e,&csq,&S);
	    	if (logs)
	    	{
	    	    log_T = plog(T);
	    	    log_e = plog(e);
	    	    (void) printf(ffmt,log_rho,log_T,log_e,csq,S);
	    	    log_rho += dlog_rho;
	    	    rho = pexp(log_rho);
	    	}
	    	else
	    	{
	    	    (void) printf(ffmt,rho,T,e,csq,S);
	    	    rho += drho;
	    	}
	    }
	    (void) printf("\n");
	    first = NO;
	}
}		/*end isobars*/

LOCAL void const_energy_plots(
	Gas_param	*params)
{
	boolean		keep_prompting = YES;
	int		i, num_points, first = YES, logs = NO;
	double		rho, p, e, T, drho, rhomin, rhomax, csq, S;
	double		log_rho, log_p, log_T, dlog_rho, log_rhomin, log_rhomax;
	char		s[Gets_BUF_SIZE];

	static Locstate state = NULL;

	if (state == NULL)
	{
	    g_alloc_state(&state,sizest);
	}

	Init_params(state,params);
	set_type_of_state(state,EGAS_STATE);

	while (keep_prompting == YES)
	{
	    screen("Enter the constant specific internal energy, "
	           "negative or default stops: ");
	    e = -HUGE_VAL;
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscanf(s,"%lf",&e);
	    if (e < 0.0)
		break;
	    Energy(state) = e;
	    if (first)
	    {
	    	screen("Do you want log scales?: ");
	    	(void) Gets(s);
	    	if (s[0] == 'y' || s[0] == 'Y')
		    logs = YES;
	    	screen("Enter the minimum and maximum densities: ");
	    	(void) Scanf("%f %f\n", &rhomin, &rhomax);
	    	screen("Enter the number of points to be plotted: ");
	    	(void) Scanf("%d\n", &num_points);
	    }
	    if (logs)
	    {
	    	(void) output();
	    	(void) printf(sfmt,"LOG_DENSITY","LOG_TEMP"," LOG_PRESSURE",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	log_rhomax = plog(rhomax);
	    	log_rhomin = plog(rhomin);
	    	log_rho = log_rhomin;
	    	dlog_rho = (log_rhomax - log_rhomin)/(num_points - 1.0);
	    }
	    else
	    {
	    	(void) output();
	    	(void) printf(sfmt,"DENSITY","TEMPERATURE","PRESSURE",
	    		      "SOUND_SP_SQ","ENTROPY");
	    	drho = (rhomax - rhomin)/(num_points - 1.0);
	    }
	    rho = rhomin;
	    for (i = 1; i <= num_points; ++i)
	    {
	    	Dens(state) = rho;
	    	get_thermo_data(state,&rho,&p,&T,&e,&csq,&S);
	    	if (logs)
	    	{
	    	    log_T = plog(T);
	    	    log_p = plog(p);
	    	    (void) printf(ffmt,log_rho,log_T,log_p,csq,S);
	    	    log_rho += dlog_rho;
	    	    rho = pexp(log_rho);
	    	}
	    	else
	    	{
	    	    (void) printf(ffmt,rho,T,p,csq,S);
	    	    rho += drho;
	    	}
	    }
	    (void) printf("\n");
	    first = NO;
	}

}		/*end const_energy_plots*/


/*
*				start_up():
*
*	Calls Routines to handle system error messages, and prints
*	Current messages at top of output or on the screen.
*
*/

#define SHIFT (--argc,++argv)

LOCAL void start_up(
	int	     argc,
	char	     **argv,
	INIT_DATA    *init)
{
	IMPORT	boolean	suppress_prompts;
	setbuf(stdin,NULL);
	init_clean_up(eosplot_clean_up,NULL);

	SHIFT;
	while (argc && argv[0][0]=='-')
	{
	    if (strncmp(argv[0],"-i",2) == 0)
	    {
	    	static char infile[1024];
	    	SHIFT;
		temporary_input_file = infile;
	    	stripcomm(infile,argv[0]);
	    	if (freopen(infile,"r",stdin) == NULL)
	    	{
	    	    screen("ERROR in start_up(), "
	    	           "can't reopen %s to %s\n","stdin",infile);
	    	    clean_up(ERROR);
	    	}
	    	suppress_prompts = YES;
	    }
	    else if (strncmp(argv[0],"-o",2) == 0)
	    {
	    	SHIFT;
	    	if (freopen(argv[0],"w",stdout) == NULL)
	    	{
	    	    screen("ERROR in start_up(), "
	    	           "can't reopen %s to %s\n","stdout",argv[0]);
	    	    clean_up(ERROR);
	    	}
	    }
	    else if (strncmp(argv[0],"-e",2) == 0)
	    {
	    	SHIFT;
	    	if (freopen(argv[0],"w",stderr) == NULL)
	    	{
	    	    screen("ERROR in start_up(), ");
	    	    screen("can't reopen %s to %s\n","stderr",argv[0]);
	    	    clean_up(ERROR);
	    	}
	    }
	    SHIFT;
	}
	set_error_immediate(stdout);
	record_print_version(stdout);
	print_title(stdout,title(init));
	init_prompting_and_debugging(init);
}		/*end start_up*/

LOCAL	void eosplot_clean_up(void)
{
	if ((temporary_input_file != NULL) && (is_io_node(pp_mynode())))
	    (void) unlink(temporary_input_file);
}		/*end eosplot_clean_up*/

LOCAL void get_thermo_data(
	Locstate	state,
	double		*rho,
	double		*p,
	double		*T,
	double		*e,
	double		*csq,
	double		*S)
{
	*rho = Dens(state);
	*p = pressure(state);
	*T = temperature(state);
	*e = specific_internal_energy(state);
	*S = entropy(state);
	*csq = sound_speed_squared(state);
	/**csq = 0.0; */
	return;
}		/*end get_thermo_data*/

LOCAL	void enter_state(
	const char *mesg,
	Locstate  state,
	int       stype,
	Gas_param *params)
{
	double rho, p, e, T;

	zero_scalar(state,params->sizest);
	set_type_of_state(state,stype);
	Init_params(state,params);
	screen("Enter the %s state density: ",mesg);
	(void) Scanf("%f\n", &rho);
	Dens(state) = rho;

	switch (stype)
	{
	case GAS_STATE:
	    screen("Enter the %s state total energy: ",mesg);
	    (void) Scanf("%f\n", &e);
	    Energy(state) = e;
	    break;
	case TGAS_STATE:
	    screen("Enter the %s state pressure: ",mesg);
	    (void) Scanf("%f\n", &p);
	    Press(state) = p;
	    break;
	case EGAS_STATE:
	    screen("Enter the %s state specific internal energy: ",mesg);
	    (void) Scanf("%f\n", &e);
	    Energy(state) = e;
	    break;
	case FGAS_STATE:
	    screen("Enter the %s state temperature: ",mesg);
	    (void) Scanf("%f\n", &T);
	    Temperature(state) = T;
	    break;
	default:
	    screen("ERROR in enter_state(), unsupported state type %s\n",
		    state_type_name(stype));
	    clean_up(ERROR);
	    break;
	}
}		/*end enter_state*/

#if defined(SESAME_CODE)
#include <geos/sesame.h>

LOCAL	void ses_print_tri_soln(
	Gas_param	*params)
{
	FILE		  *file;
	char		  *c, s[Gets_BUF_SIZE];
	SESAME_TABLE_TYPE itable;
	SESAME_EOS	  *seos = (SESAME_EOS *)params->eos;

 	screen("Enter the Table to be printed,  choices are\n"
	       "\tDENSITY TEMPERATURE\n"
	       "\tDENSITY ENERGY\n"
	       "\tDENSITY ENTROPY\n"
	       "\tPRESSURE ENTROPY\n");
	if (seos->fr[SESAME_VOLUME_PRESSURE]->interf != NULL)
	    screen("\tVOLUME PRESSURE\n");
	screen("Enter choice: ");
	(void) Gets(s);
	for (c = s; *c != '\0'; ++c)
	    *c = tolower(*c);
	if ((strncmp(s,"1",1) == 0) || (strstr(s,"temperature") != NULL))
	    itable = SESAME_RHO_TEMP;
	else if ((strncmp(s,"2",1) == 0) || (strstr(s,"energy") != NULL))
	    itable = SESAME_RHO_ENERGY;
	else if ((strncmp(s,"3",1) == 0) || ((strstr(s,"density") != NULL) &&
					     (strstr(s,"entropy") != NULL)))
	    itable = SESAME_RHO_ENTROPY;
	else if ((strncmp(s,"4",1) == 0) || ((strstr(s,"pressure") != NULL) &&
					     (strstr(s,"entropy") != NULL)))
	    itable = SESAME_PRESS_ENTROPY;
	else if ((seos->fr[SESAME_VOLUME_PRESSURE]->interf == NULL) &&
		 (strncmp(s,"5",1) == 0) || (strstr(s,"volume") != NULL))
	    itable = SESAME_VOLUME_PRESSURE;
	else
	{
	    screen("NO SUCH CHOICE %s\n",s);
	    return;
	}
	file = stdout;
	screen("Enter an optional filename for the tri solution output: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((file = fopen(s,"w")) == NULL)
	    {
	        screen("ERROR in ses_print_tri_soln(), can't open %s\n",s);
	        clean_up(ERROR);
	    }
	}
 
	switch(itable)
	{
	case SESAME_RHO_TEMP:
	    print_rt_tri_soln(file,seos);
	    break;

	case SESAME_RHO_ENERGY:
	    print_re_tri_soln(file,seos);
	    break;

	case SESAME_RHO_ENTROPY:
	    print_rs_tri_soln(file,seos);
	    break;

	case SESAME_PRESS_ENTROPY:
	    print_ps_tri_soln(file,seos);
	    break;

	case SESAME_VOLUME_PRESSURE:
	    if (seos->fr[SESAME_VOLUME_PRESSURE]->interf != NULL)
	    	print_vp_tri_soln(file,seos);
	    break;
	
	default:
	    break;
	}

	if (file != stdout)
	{
	    trace_foutput(file);
	    (void) fclose(file);
	}
}		/*end ses_print_tri_soln*/

#if defined(PHASE_CODE)
LOCAL 	void print_phase_bound(
	Gas_param	*params)
{
	SESAME_EOS	*seos = (SESAME_EOS *)params->eos;
	CURVE		**cur, *phsbdry;
	BOND		*b;
	Front		*fr = seos->fr[SESAME_RHO_TEMP];
	INTERFACE	*intfc = fr->interf;
	Locstate	lstate;
	double		p, rho, temp, S;

	for (cur = intfc->curves; cur && *cur; ++cur)
	{
	    if (wave_type(*cur) == PHASE_BOUNDARY)
	    {
	    	phsbdry = *cur;
	    }
	}

	(void) output();
	(void) printf(sfmt,"LOG_DENSITY", "LOG_PRESSURE","LOG_TEMP",
		           "ENTROPY","ENTROPY");
	for (b = (phsbdry)->first; b != NULL; b = b->next)
	{
	    lstate = left_state_at_point_on_curve(b->start,b,phsbdry);
	    rho = ses_rt_rho_from_grid(Coords(b->start)[0],seos);
	    temp = ses_rt_temp_from_grid(Coords(b->start)[1],seos);
	    p = ses_rt_coldp(lstate) + ses_rt_redp(lstate)*rho*temp;
	    S = ses_rt_S(lstate);
		
	    (void) printf(ffmt,plog(rho),plog(p),plog(temp),S,S);

	    lstate = left_state_at_point_on_curve(b->end,b,phsbdry);
	    rho = ses_rt_rho_from_grid(Coords(b->end)[0],seos);
	    temp = ses_rt_temp_from_grid(Coords(b->end)[1],seos);
	    p = ses_rt_coldp(lstate) + ses_rt_redp(lstate)*rho*temp;
	    S = ses_rt_S(lstate);
		
	    (void) printf(ffmt,plog(rho),plog(p),plog(temp),S,S);
		
	}

}		/*end print_phase_bound*/


LOCAL	void	init_phase_bdry(
	Locstate state)
{
	double pi[2],ri[2],ui[2];
	int   ni[2];
	static      Locstate Ts = NULL;

	if (!multiphase_eos(Eos(state)))
	    return;

	if (Ts == NULL)
	    g_alloc_state(&Ts,Params(state)->sizest);


	(void) printf("Calling intrsct_wv_crv\n");
	set_state_for_find_mid_state(Ts,state);
	if (intrsct_wv_crv_wth_phs_bdry(Ts,Ts,pi,ui,ri,ni,pi,ui,ri,ni,
		SHOCK,SHOCK))
	{
	    (void) printf("%-"FFMT" %-"FFMT" %-"FFMT" "
			  "%-"FFMT" %-"FFMT"\n",
			  plog(ri[0]),plog(pi[0]),0.0,0.0,0.0);
	}
}		/*end init_phase_bdry*/
#endif /* defined(PHASE_CODE) */
#endif /* defined(SESAME_CODE) */
#endif /* defined(TWOD) */
