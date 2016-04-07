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
*				giphysprompt.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <ginit/ginit.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>


	/* LOCAL Function Declarations */
LOCAL	void	g_copy_redistribution_values(INIT_DATA*,Front*);
LOCAL	void	g_no_comps(void);
LOCAL	void	g_prompt_for_hyperbolic_method(INIT_DATA*);
LOCAL	void	g_prompt_for_initial_intfc_options(INIT_DATA*);
LOCAL	void	g_prompt_for_point_propagation_options(INIT_DATA*);
LOCAL	void	g_prompt_for_interpolation_options(INIT_DATA*);
LOCAL	void	g_prompt_for_redistribute(INIT_DATA*);
LOCAL	void	g_set_redistribution_defaults(INIT_DATA*);
LOCAL	void	print_point_prop_parameters(G_PT_PROP_OPTS*);
#if defined(FULL_PHYSICS) && defined(TWOD)
LOCAL	void	g_prompt_for_wave_capture_options(INIT_DATA*);
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */

EXPORT	void	g_set_prompting_hooks(
	INIT_DATA	*init)
{
	f_init_data(init)->_set_redistribution_defaults =
	    g_set_redistribution_defaults;
	f_init_data(init)->_copy_redistribution_values =
	    g_copy_redistribution_values;
	f_init_data(init)->_prompt_for_redistribute =
	    g_prompt_for_redistribute;
	d_init_data(init)->_prompt_for_physics_options =
	    g_prompt_for_physics_options;
	g_init_data(init)->_prompt_for_composition_type =
	    g_prompt_for_composition_type;
	g_init_data(init)->_prompt_for_printing_and_plotting =
	    g_prompt_for_printing_and_plotting;
	d_init_data(init)->_set_plot_choices = g_set_plot_choices;
	d_init_data(init)->_prompt_for_initial_intfc_options =
	    g_prompt_for_initial_intfc_options;
	g_init_data(init)->_prompt_for_hyperbolic_method =
	    g_prompt_for_hyperbolic_method;
	g_init_data(init)->_prompt_for_tracked_bifurcations =
	    g_prompt_for_tracked_bifurcations;
	g_init_data(init)->_get_comp_type_type = g_get_comp_type_type;
	g_init_data(init)->_prompt_for_constant_flow_region =
	    g_prompt_for_constant_flow_region;
	g_init_data(init)->_prompt_for_ref_state = g_prompt_for_ref_state;

#if defined(FULL_PHYSICS) && defined(TWOD)
	g_init_data(init)->_prompt_for_wave_capture_options =
	    g_prompt_for_wave_capture_options;
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */

	g_init_data(init)->_prompt_for_interpolation_options =
	    g_prompt_for_interpolation_options;
	interpolation_option(init) = CONSERVATIVE_VARIABLES;

	g_init_data(init)->_prompt_for_point_propagation_options =
	    g_prompt_for_point_propagation_options;
	if (strcmp(ex_name(init),"spolars") == 0)
	{
	    g_init_data(init)->_prompt_for_maximum_number_of_components =
		g_no_comps;
	}
	else
	{
	    g_init_data(init)->_prompt_for_maximum_number_of_components =
	        g_prompt_for_maximum_number_of_components;
	}
	h_init_data(init)->_setup_available_hyperbolic_methods_list =
	    g_setup_available_hyperbolic_methods_list;
	
}		/*end g_set_prompting_hooks*/

/*
*                       g_prompt_for_physics_options():
*
*       Prompts for user selected physics options.
*/

EXPORT  void    g_prompt_for_physics_options(
	INIT_DATA       *init,
	INIT_PHYSICS    *ip)
{
	prompt_for_composition_type(init);
	prompt_for_printing_and_plotting(init);
	prompt_for_hyperbolic_method(init);

	prompt_for_gravity(init,ip);
	set_gravity(gravity_data(init)); /* Required here so that the
					  * user interface hooks are set
					  * properly for restart
					  */

	prompt_for_tracked_bifurcations(init);
	prompt_for_wave_capture_options(init);
	prompt_for_point_propagation_options(init);
	prompt_for_interpolation_options(init);
	prompt_for_maximum_number_of_components(init);
}               /*end g_prompt_for_physics_options*/


LOCAL	void	g_prompt_for_hyperbolic_method(
	INIT_DATA	*init)
{
	setup_available_hyperbolic_methods_list(init);
	if (!debugging("dflt_hyp"))
	    default_hyperbolic_method(init) = NULL;

	prompt_for_hyp_method(init);

	init_artificial_parameters(init);
}		/*end g_prompt_for_hyperbolic_method*/


LOCAL	void	g_prompt_for_initial_intfc_options(
	INIT_DATA	*init)
{
	d_prompt_for_initial_intfc_options(init);
	if (restart_io_type(init) != NULL)
	    reset_artificial_viscosity_and_heat_conduction(init);
}		/*end g_prompt_for_initial_intfc_options*/


/*
*		init_artificial_parameters():
*/


EXPORT	void init_artificial_parameters(
	INIT_DATA		*init)
{
	AVISC	*avisc = &global_artifical_viscosity_parameters(init);

	/*
	*			IMPORTANT !!!!!!
	*
	*   NOTE: default sp_coef = 1.0/(sqrt(1.0 + 0.25*sqr(chi)) - 0.5*(chi))
	*   where chi = lapidus_visc_coef above.  If you change this number
	*   you are responsible for reading this message and changing the
	*   default value of sp_coef.
	*/

	zero_scalar(avisc,sizeof(AVISC));

	avisc->hyp_method = hyperbolic_method_name(init);
	avisc->use_lin_av =		NO;
	avisc->use_upwind_av = 		NO;
	avisc->linear_visc_coef = 	0.0;
	avisc->upwind_visc_coef =	0.0;
	avisc->heat_cond =		0.0;
	avisc->char_speed_cutoff =	0.0;
	avisc->dynamic_st =              0.0;
        avisc->contact_detector =       0.0;
	if (strstr(avisc->hyp_method,"LAX_WENDROFF") != NULL)
	{
	    avisc->use_lapidus_av =	YES;
	    avisc->use_msf = 		NO;
	    avisc->lapidus_visc_coef =	0.5;
	    avisc->msf_ieta =		0.0;
	    avisc->min_shock_jump =	0.0;
	    avisc->min_sp_vol_jump =	0.0;
	    avisc->sp_coef =		1.2808;
	}
	else if (strstr(avisc->hyp_method,"VECTOR_MUSCL") != NULL)
	{
	    avisc->use_lapidus_av =	NO;
	    avisc->use_msf = 		YES;
	    avisc->lapidus_visc_coef =	0.0;
	    avisc->msf_ieta =		2.0;
	    avisc->min_shock_jump =	0.25;
	    avisc->min_sp_vol_jump =	1e-06;
	    avisc->sp_coef =		1.0;
            avisc->contact_detector =   0.1;
	}
	else if (strstr(avisc->hyp_method,"PLM") != NULL)
	{
	    avisc->use_lapidus_av =	NO;
	    avisc->use_msf = 		YES;
	    avisc->lapidus_visc_coef =	0.0;
	    avisc->msf_ieta =		2.0;
	    avisc->min_shock_jump =	0.25;
	    avisc->min_sp_vol_jump =	1e-06;
	    avisc->sp_coef =		1.0;
            avisc->contact_detector =   0.1;
	}

	set_default_artificial_viscosity(avisc);
	prompt_for_artificial_viscosity_and_heat_conduction(init,
		"global default ","",YES,avisc);
	set_default_artificial_viscosity(avisc);
	debug_print("init","heat conductivity initialized to %g\n",avisc->heat_cond);
}		/*end init_artificial_parameters*/

LOCAL	void g_prompt_for_interpolation_options(
	INIT_DATA	*init)
{
	char		s[Gets_BUF_SIZE];
	RECT_GRID	*gr = &Comp_grid(init);
	GEOMETRY_REMAP	remap = gr->Remap.remap;

	rotational_symmetry_interpolation_flag(init) = UNWEIGHTED_INTERP_COEFS;
	if (remap == INVALID_REMAP)
	{
	    screen("ERROR in g_prompt_for_interpolation_options(), "
		   "invalid remap value\n");
	    clean_up(ERROR);
	}
	if (remap != IDENTITY_REMAP)
	    rotational_symmetry_interpolation_flag(init) =
		VOLUME_WEIGHTED_INTERP_COEFS;

	screen("The current defaults for the "
	       "linear interpolation options are\n");
	switch(interpolation_option(init))
	{
	case CONSERVATIVE_VARIABLES:
	    screen("\tLinear interpolation based on conserved variables\n");
	    break;
	case THERMODYNAMIC_VARIABLES:
	    screen("\tLinear interpolation based on thermodynamic variables\n");
	    break;
	default:
	    screen("ERROR in g_prompt_for_interpolation_options(), "
		   "invalid interpolation option %d\n",
		   interpolation_option(init));
	    clean_up(ERROR);
	}
	if (remap != IDENTITY_REMAP)
	{
	    switch (rotational_symmetry_interpolation_flag(init))
	    {
	    case UNWEIGHTED_INTERP_COEFS:
	        screen("\tUnweighted interpolation coefficients\n");
		break;
	    case VOLUME_WEIGHTED_INTERP_COEFS:
	        screen("\tVolume interpolation coefficients\n");
		break;
	    case RGAS_WEIGHTED_INTERP_COEFS:
	        screen("\tLinear interpolation using volume averaged states\n");
		break;
	    case RGAS_VOLUME_WEIGHTED_INTERP_COEFS:
		screen("\tVolume weighted/volume averages interpolation\n");
		break;
	    default:
	        screen("ERROR in g_prompt_for_interpolation_options(), "
		       "invalid rotational symmetry interpolation flag %d\n",
		       rotational_symmetry_interpolation_flag(init));
	        clean_up(ERROR);
	    }
	}

	screen("Use current defaults for linear interpolation options\n");
	screen("\t(default = y): ");
	(void) Gets(s);

	if (s[0] != 'n' && s[0] != 'N')
	    return;

	screen_print_long_string("State interpolation is based on "
	                         "linear interpolation of the conserved "
	                         "variables [C] or on "
	                         "the logarithms of the density and "
	                         "temperature [T].\n");
	screen("Enter choice for state interpolators ");
	switch(interpolation_option(init))
	{
	case CONSERVATIVE_VARIABLES:
	    screen("[C]: ");
	    break;
	case THERMODYNAMIC_VARIABLES:
	    screen("[T]: ");
	    break;
	default:
	    screen("ERROR in g_prompt_for_interpolation_options(), "
		   "invalid interpolation option %d\n",
		   interpolation_option(init));
	    clean_up(ERROR);
	}
	(void) Gets(s);
	switch (s[0])
	{
	case 'T':
	case 't':
	    interpolation_option(init) = THERMODYNAMIC_VARIABLES;
	    break;
	case 'C':
	case 'c':
	default:
	    interpolation_option(init) = CONSERVATIVE_VARIABLES;
	    break;
	}

	if (remap != IDENTITY_REMAP)
	{
	    screen("\n");
	    screen_print_long_string("Geometric weighting can be applied to "
				     "the interpolation coefficients using "
				     "one of four different algorithms.\n");
	    screen("\t(%d) Linear interpolation using weights "
		   "based on the volume\n"
	           "\t\tseparating the geometric locations of the states.\n",
	           VOLUME_WEIGHTED_INTERP_COEFS);
	    screen("\t(%d) Linear interpolation using volume averaged states\n",
	           RGAS_WEIGHTED_INTERP_COEFS);
	    screen("\t(%d) A combination of methods (1) and (2).\n",
	           RGAS_VOLUME_WEIGHTED_INTERP_COEFS);
	    screen("\t(%d) No geometric weighting\n",UNWEIGHTED_INTERP_COEFS);
	    screen("Enter choice  [%d]: ",
	           rotational_symmetry_interpolation_flag(init));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        int  iflag;
	        (void) sscanf(s,"%d",&iflag);
	        switch (iflag)
	        {
	        case UNWEIGHTED_INTERP_COEFS:
	            rotational_symmetry_interpolation_flag(init) =
		        UNWEIGHTED_INTERP_COEFS;
		    break;
	        case VOLUME_WEIGHTED_INTERP_COEFS:
	            rotational_symmetry_interpolation_flag(init) =
		        VOLUME_WEIGHTED_INTERP_COEFS;
		    break;
	        case RGAS_WEIGHTED_INTERP_COEFS:
	            rotational_symmetry_interpolation_flag(init) =
		        RGAS_WEIGHTED_INTERP_COEFS;
		    break;
	        case RGAS_VOLUME_WEIGHTED_INTERP_COEFS:
	            rotational_symmetry_interpolation_flag(init) =
		        RGAS_VOLUME_WEIGHTED_INTERP_COEFS;
		    break;
	        default:
		    screen("ERROR in g_prompt_for_interpolation_options(), "
		           "invalid interpolation weigth flag %d\n",iflag);
		    clean_up(ERROR);
		    break;
	        }
	    }
	}
}		/*end g_prompt_for_interpolation_options*/


LOCAL	void	g_prompt_for_point_propagation_options(
	INIT_DATA	*init)
{
	G_PT_PROP_OPTS	*opts = &pt_propagate_opts(init);
	NptWSpeedOpts	*wsopts = &opts->npt_w_speed_opts;
	char		s[Gets_BUF_SIZE];
	int		dim = i_intfc(init)->dim;
	RECT_GRID	*gr = &Comp_grid(init);
	GEOMETRY_REMAP	remap = gr->Remap.remap;

	/*Set defaults*/
	opts->use_unsplit_pt_prop = NO;
	wsopts->Mach_tol = 0.25;
	wsopts->Wall_limiter = 1.0;
	wsopts->A_tol = 0.25;
	wsopts->vector_moc = MOC_PLUS_RH;
	wsopts->scalar_moc = RIEMANN;
	wsopts->_scalar_filter_outgoing_waves_at_contact = NO;
	wsopts->use_cheap_ahead_state_moc = NO;
	wsopts->use_neumann_cheap_moc = YES;	/* putting this as default */
	wsopts->vector_ahead_state_moc = NULL;
	wsopts->neumann_moc = NULL;
	wsopts->remap = remap;
	wsopts->geom_source_method = BACKWARD_EULER;
	print_point_prop_parameters(opts);

	screen("Use defaults for point propagation operators (dflt = y): ");
	(void) Gets(s);
	if (s[0] != 'n' && s[0] != 'N')
	    return;

	if (is_rotational_symmetry() &&
	    (rotational_symmetry() > 0) && (dim == 2))
	    opts->use_unsplit_pt_prop = YES;/*Default for rotational symmetry*/

#if defined(TWOD)
	if (dim == 2)
	{
	    screen("Use unsplit point propagate operator (dflt = %s): ",
	    	   (opts->use_unsplit_pt_prop == YES) ? "y" : "n");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    	opts->use_unsplit_pt_prop = YES;
	    else if (s[0] == 'n' || s[0] == 'N')
	    	opts->use_unsplit_pt_prop = NO;
	}
#endif /* defined(TWOD) */

	screen("Use defaults for three point wave speed (dflt = y): ");
	(void) Gets(s);
	if (s[0] != 'n' && s[0] != 'N')
	    return;

	screen("Enter the Mach number tolerance for strong waves "
	       "(dflt = %g): ",wsopts->Mach_tol);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if (strcmp(s,"inf")==0 || strcmp(s,"Inf")==0 || strcmp(s,"INF")==0)
	        wsopts->Mach_tol = HUGE_VAL;
	    else
	        (void) sscan_float(s,&wsopts->Mach_tol);
	}
	if (wsopts->Mach_tol < 0.0)
	{
	    screen("ERROR in g_prompt_for_point_propagation_options(), "
	           "negative Mach_tol entered\n");
	    clean_up(ERROR);
	}

	screen("Enter the Wall limter for strong waves at Neumann boundaries "
	       "(dflt = %g): ",wsopts->Wall_limiter);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if (strcmp(s,"inf")==0 || strcmp(s,"Inf")==0 || strcmp(s,"INF")==0)
	        wsopts->Wall_limiter = HUGE_VAL;
	    else
	        (void) sscan_float(s,&wsopts->Wall_limiter);
	}
	if (wsopts->Wall_limiter < 0.0)
	{
	    screen("ERROR in g_prompt_for_point_propagation_options(), "
	           "negative Wall_limiter entered\n");
	    clean_up(ERROR);
	}

	screen("Enter the Atwood number tolerance for strong waves ");
	screen("(dflt = %g): ",wsopts->A_tol);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if (strcmp(s,"inf")==0 || strcmp(s,"Inf")==0 || strcmp(s,"INF")==0)
	        wsopts->A_tol = HUGE_VAL;
	    else
	        (void) sscan_float(s,&wsopts->A_tol);
	}
	if (wsopts->A_tol < 0.0)
	{
	    screen("ERROR in g_prompt_for_point_propagation_options(), "
	           "negative A_tol entered\n");
	    clean_up(ERROR);
	}
	screen("Use the cheap method of characteristics "
	       "for ahead states? (dflt = %s): ",
	       (wsopts->use_cheap_ahead_state_moc == YES) ? "yes" : "no");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	    wsopts->use_cheap_ahead_state_moc = YES;
	else if ((s[0] == 'n') || (s[0] == 'N'))
	    wsopts->use_cheap_ahead_state_moc = NO;

	screen("Use the cheap method of characteristics "
	       "at Neumann boundaries? (dflt = %s): ",
	    (wsopts->use_neumann_cheap_moc == YES) ? "yes" : "no");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	     wsopts->use_neumann_cheap_moc = YES;
	else if ((s[0] == 'n') || (s[0] == 'N'))
	     wsopts->use_neumann_cheap_moc = NO;

	screen("Enter choice for method of characteristics type ");
	screen("across vector waves\n");
	screen("\tChoices are\n");
	screen("\t\tRiemann problems (R)\n");
	screen("\t\tMethod of Characteristics plus Riemann problems (M)\n");
	screen("\t\tMethod of characteristics plus Rankine-Hugoniot (H)\n");
	screen("Enter choice");
	switch (wsopts->vector_moc)
	{
	case MOC_PLUS_RIEMANN:
	     screen(" (dflt = M)");
	     break;
	case RIEMANN:
	     screen(" (dflt = R)");
	     break;
	case MOC_PLUS_RH:
	     screen(" (dflt = H)");
	     break;
	default:
	    break;
	}
	screen(": ");
	(void) Gets(s);
	switch(s[0])
	{
	case 'M':
	case 'm':
	     wsopts->vector_moc = MOC_PLUS_RIEMANN;
	     break;
	case 'R':
	case 'r':
	     wsopts->vector_moc = RIEMANN;
	     break;
	case 'H':
	case 'h':
	     wsopts->vector_moc = MOC_PLUS_RH;
	     break;
	}

	screen("Enter choice for method of characteristics type ");
	screen("across scalar waves\n");
	screen("\tChoices are\n");
	screen("\t\tRiemann problems (R)\n");
	screen("\t\tMethod of Characteristics plus Riemann problems (M)\n");
	screen("\t\tMethod of characteristics plus Rankine-Hugoniot (H)\n");
	screen("Enter choice");
	switch (wsopts->scalar_moc)
	{
	case MOC_PLUS_RIEMANN:
	     screen(" (dflt = M)");
	     break;
	case RIEMANN:
	     screen(" (dflt = R)");
	     break;
	case MOC_PLUS_RH:
	     screen(" (dflt = H)");
	     break;
	default:
	    break;
	}
	screen(": ");
	(void) Gets(s);
	switch(s[0])
	{
	case 'M':
	case 'm':
	     wsopts->scalar_moc = MOC_PLUS_RIEMANN;
	     break;
	case 'R':
	case 'r':
	     wsopts->scalar_moc = RIEMANN;
	     break;
	case 'H':
	case 'h':
	     wsopts->scalar_moc = MOC_PLUS_RH;
	     break;
	}

	if (wsopts->scalar_moc == RIEMANN)
	{
	     screen("Use filtering for outgoing waves at contact? "
	            "(dflt = %s): ",
	             (wsopts->_scalar_filter_outgoing_waves_at_contact==YES) ?
	             "yes" : "no");
	    (void) Gets(s);
	    if ((s[0] == 'y') || (s[0] == 'Y'))
	        wsopts->_scalar_filter_outgoing_waves_at_contact = YES;
	    else if ((s[0] == 'n') || (s[0] == 'N'))
	        wsopts->_scalar_filter_outgoing_waves_at_contact = NO;
	}

	if (remap != IDENTITY_REMAP)
	{
	    screen("Specify the method for integrating radial source terms\n");
	    screen("\tChoices are\n");
	    screen("\t\tModified Euler (M%s)\n",
	        (wsopts->geom_source_method==MODIFIED_EULER)?", default":"");
	    screen("\t\tAnalytic solution (A%s)\n",
	        (wsopts->geom_source_method==ANALYTIC_SOLUTION)?", default":"");
	    screen("\t\tBackward Euler (B%s)\n",
	        (wsopts->geom_source_method==BACKWARD_EULER)?", default":"");
	    screen("Enter choice: ");
	    (void) Gets(s);
	    switch(s[0])
	    {
	    case 'M':
	    case 'm':
	         wsopts->geom_source_method = MODIFIED_EULER;
	         break;
	    case 'A':
	    case 'a':
	         wsopts->geom_source_method = ANALYTIC_SOLUTION;
	         break;
	    case 'B':
	    case 'b':
	         wsopts->geom_source_method = BACKWARD_EULER;
	         break;
	    }
	}
	print_point_prop_parameters(opts);
}		/*end g_prompt_for_point_propagation_options*/

LOCAL	void	print_point_prop_parameters(
	G_PT_PROP_OPTS	*opts)
{
	screen("\nPoint propagation parameters have "
	       "the current default values\n");
	if (opts->use_unsplit_pt_prop == YES)
	    screen("\tUnsplit normal/tangential update\n");
	else
	    screen("\tOperator split normal/tangential update\n");
	print_g_npt_w_speed_opts(&opts->npt_w_speed_opts);
	screen("\n");
}		/*end print_point_prop_parameters*/

#if defined(FULL_PHYSICS) && defined(TWOD)
/*
*			g_prompt_for_wave_capture_options():
*/

LOCAL void	g_prompt_for_wave_capture_options(
	INIT_DATA	*init)
{
	G_WAVE_CAPTURE	*gwc;
	char		s[Gets_BUF_SIZE];

	wave_capture_data(init) = NULL;
	/*
	*  Currently wave capture is only implemented in 2D
	*/
	if (i_intfc(init)->dim != 2)
	    return;

	screen("Type yes to request automatic wave capture: ");
	(void) Gets(s);
	if (s[0] != 'y' && s[0] != 'Y')
	    return;

	scalar(&gwc,sizeof(G_WAVE_CAPTURE));
	wave_capture_data(init) = gwc;
	gwc->F_wave_capture._wave_capture = g_shock_capture;

	screen("Enter dimensionless thresholds for shock capture\n"
	       "\tinitiation, extension and contraction: ");
	(void) Scanf("%f %f %f\n",&gwc->_init_threshold,
				  &gwc->_expand_threshold,
				  &gwc->_contract_threshold);
}		/*end g_prompt_for_wave_capture_options*/
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */

LOCAL void    g_prompt_for_redistribute(
	INIT_DATA   *init)
{
	char	s[Gets_BUF_SIZE];
	int     dim = Comp_grid(init).dim;

	f_prompt_for_redistribute(init);

	if (dim == 2)
	{
	    screen("\n\t\tSmall loop control\n\n");
	    screen("Reflect small loop shocks (dflt = %s): ",
		   y_or_n(reflect_small_loop_shocks(init)));
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
		reflect_small_loop_shocks(init) = YES;        
	}

}   /*end to g_prompt_for_redistribute()*/

LOCAL	void	g_set_redistribution_defaults(
	INIT_DATA	*init)
{
	f_set_redistribution_defaults(init);
	reflect_small_loop_shocks(init) = NO;        
}		/*end g_set_redistribution_defaults*/

LOCAL	void	g_copy_redistribution_values(
	INIT_DATA *init,
	Front	  *front)
{
	f_copy_redistribution_values(init,front);
	if (Delete_small_loops_function(front) == g_delete_small_loops)
	    reflect_small_loop_shocks(init) = YES;        
	else
	    reflect_small_loop_shocks(init) = NO;        
}		/*end g_set_redistribution_defaults*/

LOCAL void g_no_comps(void)
{
	return;
}               /*end g_no_comps*/
