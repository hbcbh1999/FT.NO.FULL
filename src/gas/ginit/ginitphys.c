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
*				ginitphys.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	Front*	g_copy_front(Front*);
#if defined(USE_OVERTURE)
LOCAL   Front*  g_deep_copy_front(Front*);
#endif /* if defined(USE_OVERTURE) */
LOCAL	MAX_FRONT_SPEED	*g_alloc_MaxFrontSpeed(MAX_FRONT_SPEED*,INTERFACE*,
                                               size_t);
LOCAL	h_MaxWaveSpeed	*g_alloc_MaxWaveSpeed(h_MaxWaveSpeed*,Wave*);
LOCAL	Wave*	g_copy_wave(Wave*);
LOCAL	const Prompt_type *g_bdry_wave_type_prompt_type(void);
LOCAL	int	gi_user_read_print_interface(INIT_DATA*,const IO_TYPE*,
                                             INTERFACE*,boolean);
LOCAL	int	select_hyper_surfs(INTERFACE*,HYPER_SURF***,int);
LOCAL	int	sort_by_components(INTERFACE*,int,HYPER_SURF***);
LOCAL	int	sort_by_wave_type(INTERFACE*,int,HYPER_SURF***,int,int);
LOCAL	void	delete_tracked_waves_from_interface(Grid*,Front*,Wave*);
LOCAL	void	g_assign_wave_parameters(Wave*,Wave*);
LOCAL	void	g_init_hyperbolic_method(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	g_init_run(int*,char***,INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	g_no_hyp_methods(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	g_phys_await(Grid*,Front*,Wave*,Printplot*,boolean);
LOCAL	void	g_read_print_front(INIT_DATA*,Front*);
LOCAL	void	g_restart_clean_up(INIT_DATA*);
LOCAL	void	g_set_default_front_parameters(INIT_DATA*,Front*);
LOCAL	void	g_set_default_wave_parameters(INIT_DATA*,Wave*);
LOCAL	void	gi_prompt_for_front_options(INIT_DATA*,Front*);
LOCAL	void	gi_read_print_front_options(INIT_DATA*,Front*);
LOCAL	void	gi_set_normal_function(const char*,NORMAL_FUNCTION*,INTERFACE*);
LOCAL	void	gi_set_tangent_function(const char*,TANGENT_FUNCTION*,
                                        INTERFACE*);
LOCAL	void	init_interpolators(INIT_DATA*,Front*,Wave*);
LOCAL	void	init_point_propagation_opts(INIT_DATA*,Front*);
LOCAL	void	reset_component_numbers(INTERFACE*);
LOCAL	void	reset_flow_specified_comps(Front*);
LOCAL	void	reset_surface_tension(Front*);
LOCAL	void	unknown_problem_type(INIT_DATA*,INIT_PHYSICS*);
LOCAL	boolean	gi_user_read_print_curve(CURVE*,const IO_TYPE*,boolean);
LOCAL	void	gi_user_read_curve(CURVE*);
LOCAL	void	gi_user_read_print_node(NODE*,const IO_TYPE*,boolean);
LOCAL	void	gi_user_read_print_point(POINT*,const IO_TYPE*,boolean);
LOCAL	int	sort_by_status(int,HYPER_SURF***,int);
LOCAL	void	g_init_redistribute(INIT_DATA*,Front*);
LOCAL	void	gi_user_read_print_surface(SURFACE*,const IO_TYPE*,boolean);

EXPORT	void	set_gas_hooks(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	set_driver_hooks(init,ip);
	ip->init_run = g_init_run;
	ip->prt->_init_statistics = g_init_statistics;
	ip->initial_dt = 0.0;
	ip->_init_physical_units = NULL;
	ip->_init_physics = g_init_physics;
	ip->init_hyp_bdry = NULL;
	ip->init_point_sources = NULL;
	ip->set_Dirichlet_boundary_states = NULL;
	problem_type(ip) = UNKNOWN_PROBLEM_TYPE;
	g_iphys(ip)->_init_composition_type = g_init_composition_type;
	g_iphys(ip)->_set_basic_phys_parameters = g_set_basic_phys_parameters;
	g_iphys(ip)->_init_printing_and_plotting = g_init_printing_and_plotting;
	g_iphys(ip)->_set_printing_list = g_set_printing_list;
	g_iphys(ip)->_phys_solution = g_solution;
	g_iphys(ip)->_phys_intfc_solution = g_intfc_solution;
	g_iphys(ip)->_set_input_solution = g_set_input_solution;
	g_iphys(ip)->_init_hyperbolic_method = g_init_hyperbolic_method;
	g_iphys(ip)->_init_problem_type = g_init_problem_type;
	g_iphys(ip)->_prompt_for_problem_specific_data = unknown_problem_type;
	g_iphys(ip)->_prompt_for_boundary_state = g_prompt_for_boundary_state;
	g_iphys(ip)->_prompt_for_eos_params = g_prompt_for_eos_params;
	g_iphys(ip)->_prompt_for_eos_params_list = g_prompt_for_eos_params_list;
	g_iphys(ip)->_init_eos_params = g_init_eos_params;
	g_iphys(ip)->_prompt_for_equation_of_state =
	    g_prompt_for_equation_of_state;
	d_init_data(init)->_restart_clean_up = g_restart_clean_up;
	g_iphys(ip)->_prompt_for_wave_type = g_prompt_for_wave_type;
 	g_init_data(init)->_prompt_for_bdry_wave_type =
	    g_prompt_for_bdry_wave_type;
	g_init_data(init)->_promptForBoundaryFlagsOffset = 4;
	g_init_data(init)->_bdry_wave_type_prompt_type =
	    g_bdry_wave_type_prompt_type;
	if (strcmp(ex_name(init),"spolars") == 0)
	    g_iphys(ip)->_init_hyperbolic_method = g_no_hyp_methods;

	type_of_state(init) = GAS_STATE;
	set_g_rproblem_hooks();
	d_init_data(init)->_set_interface_hooks = gi_set_interface_hooks;
	f_init_data(init)->_read_print_front_options =
	    gi_read_print_front_options;
	f_init_data(init)->_prompt_for_front_options =
	    gi_prompt_for_front_options;
	g_set_prompting_hooks(init);
}		/*end set_gas_hooks*/

LOCAL	const Prompt_type *g_bdry_wave_type_prompt_type(void)
{
	static const Prompt_type ptype[] = {
	    {"Unknown",    "UN", 1, {UNKNOWN_BOUNDARY_TYPE} },
	    {"Periodic",   "PE", 2, {SUBDOMAIN_BOUNDARY}    },
	    {"Reflecting", "RE", 1, {REFLECTION_BOUNDARY}   },
	    {"Mixed",      "M",  1, {MIXED_TYPE_BOUNDARY}   },
	    {"Neumann",    "NE", 2, {NEUMANN_BOUNDARY}      },
	    {"No Slip Neumann",    "NO", 2, {NO_SLIP_NEUMANN_BOUNDARY}      },
	    {"Dirichlet",  "DI", 1, {DIRICHLET_BOUNDARY}    },
	    {"Passive",    "PA", 2, {PASSIVE_BOUNDARY}      },
	    {NULL,         NULL, 0, {UNKNOWN_BOUNDARY_TYPE} }
	};
	return ptype;
}		/*end g_bdry_wave_type_prompt_type*/


LOCAL void g_init_run(
	int		*pargc,
	char		***pargv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	d_init_run(pargc,pargv,init,ip);
	
	/*#bjet2 */
	communicate_comp_params(ip->root->front, max_num_comps());

	/*Clear up storage from comp_types*/
	free_comp_types();
}		/*g_init_run*/

/*ARGSUSED*/
LOCAL	void	unknown_problem_type(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	screen("ERROR in unknown_problem_type, initialization unspecifed\n"
	       "function to specify problem dependent data must be given\n");
	clean_up(ERROR);
}		/*end unknown_problem_type*/

/*
*			g_init_physics():
*
*	Initializes the physics dependent functions.
*	Called by init() in dinit.c.
*/

EXPORT void g_init_physics(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	CHART		*root = ip->root;
	Printplot	*prt = ip->prt;
	Front		*fr = root->front;
	Wave		*wave = root->wave;
	size_t		sizest;
	int		nfloats;	/* number of stat, prt, and */
					/* show_state vars in a Gas state. */

	set_composition_type(init_composition_type(ip,init,&sizest,&nfloats));
	g_set_sizeof_state(root,sizest,nfloats);

	set_basic_phys_parameters(init,ip);

	init_printing_and_plotting(ip,init,prt,nfloats);

	init_hyperbolic_method(init,ip);

	scat_wv_tol_list(fr) = scattered_wave_tolerance_list(init);
	untrack_node_list(fr) = untrack_node_options_list(init);

#if defined(FULL_PHYSICS)
	if (fr->rect_grid->dim == 2)	/* TODO THREED */
	{
	    fr->_replace_unphys_loop = g_replace_unphys_loop;
	    f_wave_capture(fr) = f_wave_capture_data(init);
	}
#endif /* defined(FULL_PHYSICS) */

	ip->set_gravity_charts = g_set_gravity_charts;

	init_point_propagation_opts(init,fr);
	init_interpolators(init,fr,wave);
}		/*end g_init_physics*/

EXPORT void g_set_basic_phys_parameters(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	CHART		*root = ip->root;
	Grid		*grid = root->grid;
	Front		*fr = root->front;
	Wave		*wave = root->wave;
	RECT_GRID	*gr = fr->rect_grid;
	int		dim = gr->dim;

	set_coord_sys(gr->Remap.remap,dim);

	if (pause(grid) != NULL)
	    pause(grid)->_await = g_phys_await;

	g_set_default_wave_parameters(init,wave);

	/* Set front function pointers */

	g_set_default_front_parameters(init,fr);

	ip->init_interface = g_init_interface;
	ip->init_parabolic = g_init_parabolic;
	ip->init_bdry_state = NULL;
	ip->init_cauchy_data_pointers = g_init_cauchy_data_pointers;
	ip->cauchy_deposition = g_cauchy_deposition;
	supports_riemann_problem_waves(ip) = NO;
}		/*end g_set_basic_phys_parameters*/


/*ARGSUSED*/
LOCAL	void	g_no_hyp_methods(
	INIT_DATA	*init,
	INIT_PHYSICS    *ip)
{
	return;
}		/*end g_no_hyp_methods*/

LOCAL	void	g_init_hyperbolic_method(
	INIT_DATA	*init,
	INIT_PHYSICS    *ip)
{
	CHART		*root = ip->root;
	Front		*fr = root->front;
	Wave		*wave = root->wave;
	
	set_hyperbolic_method(init,wave,&root->hyp_solver);
	g_set_hyp_solvers(init,wave,fr);
}		/*end g_init_hyperbolic_method*/

/*ARGSUSED*/
LOCAL	void init_interpolators(
	INIT_DATA	*init,
	Front		*fr,
	Wave		*wave)
{
	int		dim = fr->rect_grid->dim;

	switch (interpolation_option(init))
	{
	case THERMODYNAMIC_VARIABLES:
	    fr->_state_interpolator = gt_lin_comb_states;
	    switch (dim)
	    {
	    case 2:
	        fr->_tri_state_interpolator = gt_tri_lin_comb_states;
	    	wave->interpolator.linear_cell = gt_tri_interpolator;
	    	wave->interpolator.bilinear_cell = gt_quad_interpolator;
	    	break;
	    case 3:
	        fr->_tri_state_interpolator = gt_tri_lin_comb_states;
	    	wave->interpolator.linear_cell = gt_tetra_interpolator;
	    	wave->interpolator.bilinear_cell = gt_cube_interpolator;
	    	wave->interpolator.least_square = g_least_sqr_interpolator;
	    	break;
	    }
	    break;
	case CONSERVATIVE_VARIABLES:
	default:
	    fr->_state_interpolator = g_lin_comb_states;
	    switch (dim)
	    {
	    case 2:
	        fr->_tri_state_interpolator = g_tri_lin_comb_states;
	    	wave->interpolator.linear_cell = g_tri_interpolator;
	    	wave->interpolator.bilinear_cell = g_quad_interpolator;
	    	break;
	    case 3:
	        fr->_tri_state_interpolator = g_tri_lin_comb_states;
	    	wave->interpolator.linear_cell = g_tetra_interpolator;
	    	wave->interpolator.bilinear_cell = g_cube_interpolator;
	    	wave->interpolator.least_square = g_least_sqr_interpolator;
	    	break;
	    }
	    break;
	}

	switch (dim)
	{
	case 2:
	    wave->interpolator.grad_linear_cell = g_grad_tri_interpolator;
	    wave->interpolator.grad_bond = NULL;
	    wave->interpolator.grad_bilinear_cell = g_grad_quad_interpolator;
	    break;
	case 3:
	    /* TODO */
	    break;
	}
	if (is_rotational_symmetry())
	    set_interpolation_weights(
		rotational_symmetry_interpolation_flag(init));
}		/*end init_interpolators*/

LOCAL	void g_phys_await(
	Grid		*grid,
	Front		*front,
	Wave		*wave,
	Printplot	*prt,
	boolean		got_intfc_from_file)
{
	IMPORT boolean  suppress_prompts;
	INIT_DATA	*init = grid->init;
	boolean		modify_intfc = NO;
	boolean		sav_suppress_prompts;
	char		s[Gets_BUF_SIZE];

	sav_suppress_prompts = suppress_prompts;
	suppress_prompts = (interactive_prompting(init) == YES) ? NO : YES;

	if (got_intfc_from_file == YES)		/* Restart */
	{
	    screen("Do you wish to modify the interface on restart "
	           "(dflt = no): ");
	    (void) Gets(s);
	    if ((s[0] == 'y') || (s[0] == 'Y'))
		modify_intfc = YES;
	}
	else					/* Pause time */
	{
	    d_await(grid,front,wave,prt,got_intfc_from_file);

	    /* TODO: allow modification of user statistics printing.
	     * This would require upgrading init_statistics_for_problem()
	     * and resetting the last print times.  Hopefully, 
	     * the statistics controls are flexible enough that this 
	     * won't be necessary. */

	    screen("Type 'm' to modify the interface: ");
	    (void) Gets(s);
	    if ((s[0] == 'm') || (s[0] == 'M'))
		modify_intfc = YES;
	}

	if (!modify_intfc)
	{
	    suppress_prompts = sav_suppress_prompts;
	    return;
	}

	screen("Do you want to turn off tracking of "
	       "certain curves? (no dflt): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    delete_tracked_waves_from_interface(grid,front,wave);

	screen("Do you wish to reset any constant comp values? (no dflt): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    reset_flow_specified_comps(front);

	screen("Do you wish to reset the surface tension "
	       "of any curve (no dflt): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    reset_surface_tension(front);

	screen("Do you wish to reset the component numbers of any curves "
	       "(no dflt): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    reset_component_numbers(front->interf);
	suppress_prompts = sav_suppress_prompts;
}		/*end g_phys_await*/

LOCAL	void reset_surface_tension(
	Front		*front)
{
	char	   s[Gets_BUF_SIZE];
	int	   n_hss, i, nhs;
	int	   dim = front->rect_grid->dim;
	HYPER_SURF **hslist, **hs;
	double	   surf_ten;
	static const char *hsname[] = {"points", "curves", "surfaces"};

	n_hss = 0;
	screen("Enter the number of families of %s ",hsname[dim]);
	screen("to be reset (dflt = %d): ",n_hss);
	(void) Gets(s);
	if (s[0] != '\0') (void) sscanf(s,"%d",&n_hss);

	for (i = 0; i < n_hss; ++i)
	{
	    nhs = select_hyper_surfs(front->interf,&hslist,CONTACT_WAVE);
	    if (nhs != 0)
	    {
		for (hs = hslist; hs && *hs; ++hs)
		{
		    screen("Enter the new value of surface tension "
		           "(dflt = %g): ",surface_tension(*hs));
		    (void) Gets(s);
		    if (s[0] != '\0')
		    {
			(void) sscan_float(s,&surf_ten);
			surface_tension(*hs) = surf_ten;
		    }
		}
	    }
	}
}		/*end reset_surface_tension*/


LOCAL	void	reset_component_numbers(
	INTERFACE	*intfc)
{
	COMPONENT	oldcomp, newcomp;
	HYPER_SURF	**hs;
	int		ncomps, i;
	char		s[Gets_BUF_SIZE];

	ncomps = 0;
	screen("How many components do you wish to reset (dflt = %d): ",ncomps);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&ncomps);
	if (ncomps == 0)
	    return;
	for (i = 1; i <= ncomps; ++i)
	{
	    oldcomp = NO_COMP;
	    screen("Enter the %d%s old component number (dflt = continue): ",
	           i,ordinal_suffix(i));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d",&oldcomp);
	    if (oldcomp == NO_COMP)
		continue;
	    newcomp = NO_COMP;
	    screen("Enter the %d%s new component number (dflt = continue): ",
	           i,ordinal_suffix(i));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d",&newcomp);
	    if (newcomp == NO_COMP) continue;
	    for (hs = intfc->hss; hs && *hs; ++hs)
	    {
	    	if (negative_component(*hs) == oldcomp)
	    	{
	    	    negative_component(*hs) = newcomp;
		    intfc->modified = YES;
		}
		if (positive_component(*hs) == oldcomp)
		{
		    positive_component(*hs) = newcomp;
		    intfc->modified = YES;
		}
	    }
	}
}		/*end reset_component_numbers*/

LOCAL	void reset_flow_specified_comps(
	Front	*front)
{
	COMPONENT	comp;
	char		s[Gets_BUF_SIZE];
	char		*c;
	int		n;

	n = 0;
	screen("Enter a list of component numbers which are ");
	screen("to be reset from flow specified to computed: ");
	(void) Gets(s);
	if (s[0] == '\0')
		return;

	for (c = s; c[0] != '\0'; ++c)
	{
		n = c[0];
		if (!isspace(n))
		{
			(void) sscanf(c,"%d",&comp);
			for (; c[0] != '\0'; ++c)
				if (isspace(n)) break;
			if (c[0] == '\0')
				return;
			SetActiveFlowComponent(comp,front);
		}
	}
}		/*end reset_flow_specified_comps*/


/*ARGSUSED*/
LOCAL	void delete_tracked_waves_from_interface(
	Grid		*grid,
	Front		*front,
	Wave		*wave)
{
	Front		*tempfront;
	HYPER_SURF	**hs;
	char		s[Gets_BUF_SIZE];
	boolean		sav_interp = interpolate_intfc_states(front->interf);
	int		n_waves;
	int		dim = front->rect_grid->dim;
	COMPONENT	newcomp;
	HYPER_SURF	**hslist;
	int		i, nc;
	O_CURVE		Oc;
	UNTRACK_FLAG	flag;

	debug_print("untrack","Entered delete_tracked_waves_from_interface()\n");

	if (debugging("untrack"))
	{
	    (void) printf("Hypersurface list of front\n");
	    for (hs = front->interf->hss; hs && *hs; ++hs)
	        (void) printf("Hypersurface %llu\n",hypersurface_number(*hs));
	    print_interface(front->interf);
	}
	tempfront = copy_front(front);
	set_size_of_intfc_state(front->sizest);
	set_copy_intfc_states(YES);
	set_add_to_correspond_list(YES);
	tempfront->interf = copy_interface(front->interf);
	if (tempfront->interf == NULL)
	{
	    screen("ERROR in delete_tracked_waves_from_interface(), "
		   "unable to copy interface\n");
	    free_front(tempfront);
	    clean_up(ERROR);
	}
	interpolate_intfc_states(tempfront->interf) = YES;

	screen("Delete all vector waves? (dflt = no): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    if (debugging("untrack"))
	    {
		(void) printf("Removing all vector waves from interface %llu\n",
			      interface_number(tempfront->interf));
		(void) printf("Hypersurface list\n");
		for (hs = tempfront->interf->hss; hs && *hs; ++hs)
		    (void) printf("Hypersurface %llu\n",
				  hypersurface_number(*hs));
		(void) printf("Selected candidates\n");
		for (hs = tempfront->interf->hss; hs && *hs; ++hs)
		{
	            if (wave_type(*hs) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
		    {
			(void) printf("Hypersurface %llu\n",
				      hypersurface_number(*hs));
		    }
		}
		print_interface(tempfront->interf);
	    }
	    switch (dim)
	    {
	    case 1:
	        hs = tempfront->interf->hss;
		while (hs && *hs)
	        {
	            if (wave_type(*hs) < FIRST_VECTOR_PHYSICS_WAVE_TYPE)
		    {
			++hs;
	                continue;
		    }
	            if (is_forward_wave(wave_type(*hs)))
	                newcomp = negative_component((*hs));
	            else if (is_backward_wave((wave_type(*hs))))
	                newcomp = positive_component((*hs));
	            else
	                newcomp = negative_component((*hs));
		    if (untrack_point(Point_of_hs(*hs),newcomp,tempfront)==YES)
	                hs = tempfront->interf->hss;
		}
		break;
	    case 2:
restart_loop:
	        for (hs = tempfront->interf->hss; hs && *hs; ++hs)
	        {
	            if (wave_type(*hs) < FIRST_VECTOR_PHYSICS_WAVE_TYPE)
	                continue;
	            if (is_forward_wave(wave_type(*hs)))
	                newcomp = negative_component((*hs));
	            else if (is_backward_wave((wave_type(*hs))))
	                newcomp = positive_component((*hs));
	            else
	                newcomp = negative_component((*hs));
	            Oc.curve = Curve_of_hs(*hs);
	            Oc.orient = POSITIVE_ORIENTATION;
		    set_untrack_flag(flag,Oc.orient,YES,YES,YES,YES,NO);
	            if (untrack_curve(&Oc,NULL,newcomp,grid->dt,tempfront,
				      wave,NULL,flag) == YES)
	                goto restart_loop;
	        }
		break;
	    case 3: /* TODO */
		screen("ERROR in delete_tracked_waves_from_interface(), "
		       "untracking of surfaces not implemented\n");
		clean_up(ERROR);
		break;
	    }
	}

	n_waves = 0;
	screen("Enter the number of additional wave families "
	       "to be deleted (dflt = %d): ",n_waves);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&n_waves);

	screen("\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i < n_waves; ++i)
	    {
	    	nc = select_hyper_surfs(tempfront->interf,&hslist,
					UNKNOWN_WAVE_TYPE);
		if (nc == 0) continue;

		for (hs = hslist; hs && *hs; ++hs)
		{
	            if (is_forward_wave(wave_type(*hs)))
			newcomp = negative_component((*hs));
		    else if (is_backward_wave((wave_type(*hs))))
			newcomp = positive_component((*hs));
		    else
		        newcomp = negative_component((*hs));
		    screen("Enter the component of the new "
		           "loop (dflt = %d): ", newcomp);
		    (void) Gets(s);
		    if (s[0] != '\0')
			(void) sscanf(s,"%d",&newcomp);
		    (void) untrack_point(Point_of_hs(*hs),newcomp,tempfront);
		}
		screen("\n");
	    }
	    break;
	case 2:
	    for (i = 0; i < n_waves; ++i)
	    {
	    	nc = select_hyper_surfs(tempfront->interf,&hslist,
					UNKNOWN_WAVE_TYPE);
		if (nc == 0)
		    continue;

		for (hs = hslist; hs && *hs; ++hs)
		{
	            if (is_forward_wave(wave_type(*hs)))
			newcomp = negative_component((*hs));
		    else if (is_backward_wave((wave_type(*hs))))
			newcomp = positive_component((*hs));
		    else
		        newcomp = negative_component((*hs));
		    screen("Enter the component of the new "
		           "loop (dflt = %d): ", newcomp);
		    (void) Gets(s);
		    if (s[0] != '\0')
			(void) sscanf(s,"%d",&newcomp);
		    Oc.curve = Curve_of_hs(*hs);
		    Oc.orient = POSITIVE_ORIENTATION;
		    set_untrack_flag(flag,Oc.orient,YES,YES,YES,YES,NO);
		    (void) untrack_curve(&Oc,NULL,newcomp,grid->dt,tempfront,
			                 wave,NULL,flag);
		}
		screen("\n");
	    }
	    break;
	case 3: /* TODO */
	    screen("ERROR in delete_tracked_waves_from_interface(), "
	           "untracking of surfaces not implemented\n");
	    clean_up(ERROR);
	    break;
	}

	assign_interface_and_free_front(front,tempfront);
	interpolate_intfc_states(front->interf) = sav_interp;

	reinit_hyp_solution_function(wave,front);

	debug_print("untrack","Left delete_tracked_waves_from_interface()\n");
}		/*end delete_tracked_waves_from_interface*/

LOCAL	int select_hyper_surfs(
	INTERFACE	*intfc,
	HYPER_SURF	***hslist,
	int		w_type)
{
	char	   s[Gets_BUF_SIZE];
	int	   nhs;
	int	   nfam;
	int	   dim = intfc->dim;
	static const char *hsname[] = { "point", "curve", "surface"};

	*hslist = NULL;
	nhs = 0;

	nfam = 1;
	screen("Enter the number of %ss in this family "
	       "(dflt = 1, a = all): ",hsname[dim-1]);
	(void) Gets(s);
	if (s[0] == 'a' || s[0] == 'A')
	{
	    nfam = -1;
	}
	else if (s[0] != '\0')
	{
	    (void) sscanf(s,"%d",&nfam);
	    if (nfam <= 0)
	    {
	    	*hslist = NULL;
	    	return 0;
	    }
	}

	nhs = sort_by_wave_type(intfc,nhs,hslist,w_type,nfam);
	if (nfam > 0 && nhs < nfam)
	{
	    screen("Sorry but I can't find %d %s%s of this type\n"
	           "so I will do nothing.\n",nfam,hsname[dim],
					     (nfam == 1) ? "" : "s");
	    *hslist = NULL;
	    return 0;
	}
	if (dim == 2)
	{
	    if (nfam > 0 && nhs != nfam)
	    	nhs = sort_by_status(nhs,hslist,nfam);
	    if (nfam > 0 && nhs < nfam)
	    {
	    	screen("Sorry but I can't find %d %s of this type\n"
	    	       "so I will do nothing.\n",
	    	       nfam, (nfam == 1) ? "curve" : "curves");
	    	*hslist = NULL;
	    	return 0;
	    }
	}
	if (nfam > 0 && nhs != nfam)
	    nhs = sort_by_components(intfc,nhs,hslist);
	if (nfam > 0 && nhs < nfam)
	{
	    screen("Sorry but I can't find %d %s%s of this type\n"
	           "so I will do nothing.\n",nfam,
		   hsname[dim], (nfam == 1) ? "" : "s");
		*hslist = NULL;
		return 0;
	}
	else if (nfam > 0 && nhs != nfam)
	{
	    screen("I still have more than %d candidate%s "
		   "to select.\n",nfam,(nfam == 1) ? "" : "s");
	    screen("Do you want to select them all? (no dflt): ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
		return nhs;
	    *hslist = NULL;
	    return 0;
	}
	return nhs;
}		/*end select_hyper_surfs*/

LOCAL	int sort_by_wave_type(
	INTERFACE	*intfc,
	int		nhs,
	HYPER_SURF	***hslist,
	int		w_type,
	int		nfam)
{
	HYPER_SURF **hs;
	char	   s[Gets_BUF_SIZE];
	int	   dim = intfc->dim;
	static const char *hsname[] = { "point", "curve", "surface"};

	if (w_type == UNKNOWN_WAVE_TYPE)
	{
	    screen("Enter the wave type of the %ss to be selected,\n",
	    	hsname[dim]);
	    screen("\tChoices are\n");
	    screen("\t\tSHOCK (S)\n");
	    screen("\t\tRAREFACTION LEADING EDGE (RL)\n");
	    screen("\t\tRAREFACTION TRAILING EDGE (RT)\n");
	    screen("\t\tCONTACT (C)\n");
	    screen("Enter choice: ");
	    (void) Gets(s);
	    if (strncasecmp(s,"S",1) == 0)
	    	w_type = SHOCK_WAVE;
	    else if (strncasecmp(s,"C",1) == 0)
	    	w_type = CONTACT_WAVE;
	    else if (strncasecmp(s,"RL",2) == 0)
	    	w_type = RAREF_LEADING_EDGE;
	    else if (strncasecmp(s,"RT",2) == 0)
	    	w_type = RAREF_TRAILING_EDGE;
	    else
	    	return 0;
	}

	switch(w_type)
	{
	case SHOCK_WAVE:
	    for (hs = intfc->hss; hs && *hs; ++hs)
	    {
	    	if (is_shock_wave(wave_type(*hs)))
	    	{
	    	    ++nhs;
	    	    if (!add_to_pointers(*hs,hslist))
	    	    {
	    	        screen("ERROR in sort_by_wave_type(), "
	    	               "add_to_pointers() failed\n");
	    	        clean_up(ERROR);
	    	    }
	    	}
	    }
	    break;
	case RAREF_LEADING_EDGE:
	    for (hs = intfc->hss; hs && *hs; ++hs)
	    {
	    	if (is_rarefaction_leading_edge(wave_type(*hs)))
	    	{
	    	    ++nhs;
	    	    if (!add_to_pointers(*hs,hslist))
	    	    {
	    	        screen("ERROR in sort_by_wave_type(), "
	    	               "add_to_pointers() failed\n");
	    	        clean_up(ERROR);
	    	    }
	    	}
	    }
	    break;
	case RAREF_TRAILING_EDGE:
	    for (hs = intfc->hss; hs && *hs; ++hs)
	    {
	    	if (is_rarefaction_trailing_edge(wave_type(*hs)))
	    	{
	    	    ++nhs;
	    	    if (!add_to_pointers(*hs,hslist))
	    	    {
	    	        screen("ERROR in sort_by_wave_type(), "
	    	               "add_to_pointers() failed\n");
	    	        clean_up(ERROR);
	    	    }
	    	}
	    }
		break;
	case CONTACT_WAVE:
	    for (hs = intfc->hss; hs && *hs; ++hs)
	    {
	    	if (wave_type(*hs) == CONTACT)
	    	{
	    	    ++nhs;
	    	    if (!add_to_pointers(*hs,hslist))
	    	    {
	    	        screen("ERROR in sort_by_wave_type(), "
	    	               "add_to_pointers() failed\n");
	    		clean_up(ERROR);
	    	    }
	    	}
	    }
	    break;
	}
	if (nhs == 0) 
	{
	    screen("Sorry but no wave of this wave type found.\n");
	    return 0;
	}
	if ((nfam > 0) && (nhs != nfam) && (w_type == SHOCK_WAVE || 
	        		w_type == RAREF_LEADING_EDGE || 
	        		w_type == RAREF_TRAILING_EDGE))
	{
	    screen("Is this a forward (F) or a backward (B) wave?: ");
	    (void) Gets(s);
	    if (s[0] == 'f' || s[0] == 'F')
	    {
	        for (hs = *hslist; hs && *hs; ++hs)
	        {
	            if (is_backward_wave(wave_type(*hs)))
	            {
	                if (!delete_from_pointers(*hs,hslist))
			{
			    screen("ERROR in sort_by_wave_type(), "
			           "delete_from_pointers() failed\n");
			    clean_up(ERROR);
			}
			--hs;
			--nhs;
	            }
	        }
	    }
	    if (s[0] == 'b' || s[0] == 'B')
	    {
	        for (hs = *hslist; hs && *hs; ++hs)
	        {
	            if (is_forward_wave(wave_type(*hs)))
		    {
			if (!delete_from_pointers(*hs,hslist))
			{
			    screen("ERROR in sort_by_wave_type(), "
			           "delete_from_pointers() failed\n");
			    clean_up(ERROR);
			}
			--hs;
			--nhs;
		    }
	        }
	    }
	}
	return nhs;
}		/*end sort_by_wave_type*/

LOCAL	int sort_by_components(
	INTERFACE	*intfc,
	int		nhs,
	HYPER_SURF	***hslist)
{
	COMPONENT  neg_c, pos_c;
	HYPER_SURF **hs;
	char	   s[Gets_BUF_SIZE];
	int	   dim = intfc->dim;
	static const char *hsname[] = { "point", "curve", "surface"};
	static const char *ncname[] = { "left", "left", "negative" };
	static const char *pcname[] = { "right", "right", "positive" };

	screen("Enter the %s and %s component of\n",ncname[dim],pcname[dim]);
	screen("\tthe %s to be selected, (p prints interface): ",hsname[dim]);
	(void) Gets(s);
	if (s[0] == 'p' || s[0] == 'P')
	{
	    print_interface(intfc);
	    screen("Enter the %s and %s component of\n",
	    	   ncname[dim],pcname[dim]);
	    screen("\tthe %s to be selected, (p prints interface): ",
	    	   hsname[dim]);
	    (void) Gets(s);
	}
	(void) sscanf(s,"%d %d",&neg_c,&pos_c);
	for (hs = *hslist; hs && *hs; ++hs)
	{
	    if (negative_component((*hs)) != neg_c ||
		positive_component((*hs)) != pos_c)
	    {
	    	if (!delete_from_pointers(*hs,hslist))
	    	{
	    	    screen("ERROR in sort_by_components(), "
	    	           "delete_from_pointers() failed\n");
	    	    clean_up(ERROR);
	    	}
	    	--hs;
	    	--nhs;
	    }
	 }
	 return nhs;
}		/*end sort_by_components*/

LOCAL	int sort_by_status(
	int		nc,
	HYPER_SURF	***hslist,
	int		nfam)
{
	CURVE		*c;
	HYPER_SURF	**hs;
	char		s[Gets_BUF_SIZE];
	int		status;
	int		start, end;

	screen("Enter the status at one of the nodes of the curve: ");
	(void) Gets(s);
	status = read_curve_status_from_string(s);
	for (hs = *hslist; hs && *hs; ++hs)
	{
	    c = Curve_of_hs(*hs);
	    if (start_status(c) != status && end_status(c) != status)
	    {
	    	if (!delete_from_pointers(*hs,hslist))
	    	{
	    	    screen("ERROR in sort_by_status(), "
	    	           "delete_from_pointers() failed\n");
	    	    clean_up(ERROR);
	    	}
	    	--hs;
	    	--nc;
	    }
	}
	if (nc <= nfam) return nc;

	screen("Enter the start status of the curve: ");
	(void) Gets(s);
	start = read_curve_status_from_string(s);
	screen("Enter the end status of the curve: ");
	(void) Gets(s);
	end = read_curve_status_from_string(s);

	for (hs = *hslist; hs && *hs; ++hs)
	{
	    c = Curve_of_hs(*hs);
	    if (start_status(c) != start || end_status(c) != end)
	    {
	    	if (!delete_from_pointers(*hs,hslist))
	    	{
	    	    screen("ERROR in sort_by_status(), "
	    	           "delete_from_pointers() failed\n");
	    	    clean_up(ERROR);
	    	}
	    	--hs;
	    	--nc;
	    }
	 }
	 return nc;
}		/*end sort_by_status*/

/* July 13 2004: Myoung-Nyoun Kim: local Lax-Friedrichs for small curve */
LOCAL   void g_tan_curve_propagate(
        Front           *fr,
        Front           *newfr,
        INTERFACE       *tempintfc,
        CURVE           *tempc,
        CURVE           *newc,
        double           dt)
{
        boolean catch_all = NO;
        void (*one_side_npt_tang_solver)(double,double,Tan_stencil*,Locstate,
                                  struct _Front*);
        if ( debugging("local_LF") )
        {
            COMPONENT ncp = negative_component(tempc),pcp = positive_component(tempc);
            boolean      change = NO;

            if ( ( ncp == 2 && pcp == 3 ) || ( ncp == 3 && pcp == 2 ) ||
		 ( ncp == 2 && pcp == 4 ) || ( ncp == 4 && pcp == 2 ) ) change = YES;
            if ( change )
            {
                BOND     *bb;
                Locstate lst, rst;
                double    Tl, Tr;
                if  ( tempc->num_points < 30 )
                    catch_all = YES;
                else
                {
                   for ( bb = tempc->first; bb != NULL; bb = bb->next )
                   {
                       lst = left_state_at_point_on_curve(bb->end,bb,tempc);
                       rst = right_state_at_point_on_curve(bb->end,bb,tempc);
                       Tl  = temperature(lst);
                       Tr  = temperature(rst);
                       if ( ( (pcp == 3 || pcp == 4) && Tr > Tl * 1.2) ||
                            ( (ncp == 3 || ncp == 4) && Tl > Tr * 1.2) )
                       {
                           catch_all = YES;
                           break;
                       }
                   }
                }
            }
            if ( catch_all )
            {
                one_side_npt_tang_solver = fr->_one_side_npt_tang_solver;
                fr->_one_side_npt_tang_solver = LFoblique;

                if ( debugging("check_local_LF") )
                {
                    printf("small curve detected\n");
                    print_curve(tempc);
                }
            }
        }

        f_tan_curve_propagate(fr,newfr,tempintfc,tempc,newc,dt);
        if (catch_all)
        {
            fr->_one_side_npt_tang_solver = one_side_npt_tang_solver;
        }
}

LOCAL	void init_point_propagation_opts(
	INIT_DATA	*init,
	Front		*fr)
{
	G_PT_PROP_OPTS	*opts = &pt_propagate_opts(init);
	int		dim = fr->rect_grid->dim;

	switch (dim)
	{
	case 1:
	    fr->snd_node_propagate = NULL;
	    fr->tan_curve_propagate = NULL;
	    fr->curve_propagate = NULL;
	    break;
	case 2:
	    /* July 13 2004: Myoung-Nyoun Kim: local Lax-Friedrichs for small curve */
            if ( debugging("local_LF") )
                fr->tan_curve_propagate = g_tan_curve_propagate;
            else
	        fr->tan_curve_propagate = f_tan_curve_propagate;
	    fr->snd_node_propagate = g_snd_node_propagate;
	    fr->curve_propagate = g_curve_propagate2d;
	    fr->_compute_force_and_torque = g_compute_force_and_torque;
	    break;
	case 3:
	    fr->snd_node_propagate = NULL;/* TODO */
	    fr->tan_curve_propagate = NULL;/* TODO */
	    fr->curve_propagate = g_curve_propagate_3d;
	    break;
	}

	set_point_propagate(fr,opts->use_unsplit_pt_prop);
	set_npt_wspeed_options(&opts->npt_w_speed_opts);
}		/*end init_point_propagation_opts*/

EXPORT  void gi_set_interface_hooks(
	int		dim,
	INIT_DATA	*init)
{
	I_USER_INTERFACE *iuh = i_user_hook(dim);
	F_USER_INTERFACE *fuh = f_user_hook(dim);

	g_set_interface_hooks(dim,init);

	/* Set hooks for interface initialization */

	iuh->_user_read_print_interface = gi_user_read_print_interface;

	iuh->_user_read_curve = gi_user_read_curve;
	iuh->_user_read_print_curve = gi_user_read_print_curve;
	iuh->_user_read_print_node = gi_user_read_print_node;
	iuh->_user_read_print_surface = gi_user_read_print_surface;

	if (dim == 1)
	    iuh->_user_read_print_point = gi_user_read_print_point;

	fuh->_read_wave_type_from_string = g_read_wave_type_from_string;
	iuh->_read_boundary_type_from_string = g_read_wave_type_from_string;
	fuh->_reflect_state = g_reflect_state;
	fuh->_alloc_state = g_alloc_state;
	fuh->_alloc_intfc_state = g_alloc_intfc_state;
	fuh->_clear_state = g_clear_state;
	fuh->_obstacle_state = g_obstacle_state;
	fuh->_read_print_state_data = g_read_print_state_data;
	fuh->_read_print_boundary_state_data = g_read_print_boundary_state_data;
	fuh->_alloc_MaxFrontSpeed = g_alloc_MaxFrontSpeed;
	fuh->_set_tangent_function = gi_set_tangent_function;
	fuh->_set_normal_function = gi_set_normal_function;
	switch (interpolation_option(init))
	{
	case THERMODYNAMIC_VARIABLES:
	    fuh->_bi_interpolate_intfc_states = gt_lin_comb_states;
	    fuh->_tri_interpolate_intfc_states = gt_tri_lin_comb_states;
	    break;
	case CONSERVATIVE_VARIABLES:
	default:
	    fuh->_bi_interpolate_intfc_states = g_lin_comb_states;
	    fuh->_tri_interpolate_intfc_states = g_tri_lin_comb_states;
	    break;
	}
	if (debugging("verbose"))
	    fuh->_fprint_intfc_state = g_verbose_fprint_intfc_state;
	else
	    fuh->_fprint_intfc_state = g_fprint_intfc_state;
	fuh->_read_hsbdry_type_from_string = g_read_hsbdry_type_from_string;
	if (dim == 2)
	{
	    fuh->_form_subintfc_via_communication =
		g_form_subintfc_via_communication2d;
#if defined(USE_OVERTURE)
            fuh->_form_patch_subintfc_via_cut = g_form_patch_subintfc_via_cut2d;
            fuh->_form_patch_subintfc = g_form_patch_subintfc_2d;
            fuh->_assembly_fine_patch_fronts_to_one = g_assembly_fine_patch_fronts_to_one;
#endif /* if defined(USE_OVERTURE) */
	}
}		/*end gi_set_interface_hooks*/

LOCAL	void gi_user_read_curve(
	CURVE		*curve)
{
	f_user_read_curve(curve);
	if (curve->interface->dim == 2)
	{
	    char type[120];
	    screen("Enter start status: ");
	    (void) Scanf("%s\n",type);
	    start_status(curve) = read_curve_status_from_string(type);
	    screen("Enter end status: ");
	    (void) Scanf("%s\n",type);
	    end_status(curve) = read_curve_status_from_string(type);
	    surface_tension(curve) =
	    	prompt_for_surface_tension(wave_type(curve),"for the curve ");
	}
}		/*end gi_user_read_curve*/

LOCAL	boolean gi_user_read_print_curve(
	CURVE	      *curve,
	const IO_TYPE *io_type,
	boolean          overlay)
{
	if (!f_user_read_print_curve(curve,io_type,overlay))
	{
	    (void) printf("WARNING in gi_user_read_print_curve(), "
	                  "f_user_read_print_curve() failed\n");
	    return NO;
	}
	if (curve->interface->dim == 2)
	{
	    FILE      *file = io_type->file;
	    INTERFACE *intfc = curve->interface;
	    char type[120];

	    if (!fgetstring(file,"curve->start_status = "))
	    {
	        (void) printf("WARNING in gi_user_read_print_curve(), "
	                      "can't find start_status\n");
	        return NO;
	    }
	    if (fscanf(file,"%s",type) != 1)
	    {
	        (void) printf("WARNING in gi_user_read_print_curve(), "
	                      "can't read start_status type\n");
	        return NO;
	    }
	    start_status(curve) = read_curve_status_from_string(type);
	    if (!fgetstring(file,"curve->end_status = "))
	    {
	        (void) printf("WARNING in gi_user_read_print_curve(), "
	                      "can't find end_status\n");
	        return NO;
	    }
	    if (fscanf(file,"%s",type) != 1)
	    {
	        (void) printf("WARNING in gi_user_read_print_curve(), "
	                      "can't read end type\n");
	        return NO;
	    }
	    end_status(curve) = read_curve_status_from_string(type);
	    surface_tension(curve) = read_print_float(
	                                  "curve->surface_tension = ",0.0,
					  io_type);
	    if (fgetstring(file,"curve->layer_index = "))
	    {
	        if (fscanf(file,"%d",&layer_index(Hyper_surf(curve))) != 1)
		{
	            (void) printf("WARNING in gi_user_read_print_curve(), "
	                          "can't read layer index\n");
	            return NO;
		}
	        num_layers(intfc) = max(num_layers(intfc),
	    	                               layer_index(Hyper_surf(curve)));
	    }
	    if (wave_type(curve) == NEUMANN_BOUNDARY)
	        if (fgetstring(file,"no slip = "))
		{
		    no_slip(Hyper_surf(curve)) = fread_boolean(file);
		    if (no_slip(Hyper_surf(curve)) == YES)
		        adherence_coeff(Hyper_surf(curve)) =
		    	    read_print_float("adherence coefficient = ",1.0,
					io_type);
		}
	}
	return YES;
}		/*end gi_user_read_print_curve*/

LOCAL	void gi_user_read_print_node(
	NODE	      *node,
	const IO_TYPE *io_type,
	boolean          overlay)
{
	FILE	*file = io_type->file;
	char	s[80];
	long	current;
	int	i, c, dim = node->interface->dim;

	f_user_read_print_node(node,io_type,overlay);
	if (node->interface &&
			interface_type(node->interface) != PHYSICAL_INTERFACE)
		return;

	current = ftell(file);
	(void) fscanf(file,"%s",s);
	(void) fseek(file,current,SEEK_SET);
	if (strncmp(s,"End",3) == 0)
	    return;

	(void) fgetstring(file,"Node Velocity  ");
	if ((c = getc(file)) != '\f') /* NOBINARY */
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < dim; ++i)
	    {
	    	(void) fscanf(file,"%*s%*s");
	    	(void) fscan_float(file,&Node_vel(node)[i]);
	    }
	}
	else
	{
	    (void) getc(file);
	    (void) read_binary_real_array(Node_vel(node),dim,io_type);
	}
	return;
}		/*end gi_user_read_print_node*/

LOCAL	void	gi_set_tangent_function(
	const char       *s,
	TANGENT_FUNCTION *tf,
	INTERFACE        *intfc)
{
	if (strcmp(s,"cc_tangent") == 0)
	    set_cc_tangent_function(tf);
	else if (strcmp(s,"tangent_to_origin") == 0)
	    set_tangent_to_origin(tf);
	else
	    f_set_tangent_function(s,tf,intfc);
}		/*end gi_set_tangent_function*/

LOCAL	void	gi_set_normal_function(
	const char      *s,
	NORMAL_FUNCTION *nf,
	INTERFACE       *intfc)
{
	if (strcmp(s,"cc_normal") == 0)
	    set_cc_normal_function(nf);
	if (strcmp(s,"temp_tnode_normal") == 0)
	    set_temp_tnode_normal(nf);
	if (strcmp(s,"temp_mnode_normal") == 0)
	    set_temp_mnode_normal(nf);
	else if (strcmp(s,"normal_to_origin") == 0)
	    set_normal_to_origin(nf);
	else
	    f_set_normal_function(s,nf,intfc);
}		/*end gi_set_tangent_function*/

LOCAL	int gi_user_read_print_interface(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	INTERFACE     *infc,
	boolean          overlay)
{
	FILE          *file = io_type->file;
	static OUTPUT *oput = NULL;

	if (!f_user_read_print_interface(init,io_type,infc,overlay))
	    return NO;

	oput = save_read_file_variables(file,oput);
	if (next_output_line_containing_string(file,"INTERFACE TYPE") == NULL)
	{
	    /* Old style output file */
	    reset_read_file_variables(oput);
	    interface_type(infc) = PHYSICAL_INTERFACE;
	    g_read_print_Gas_param_list(init,io_type,infc);
	    set_params_list(infc);
	    return YES;
	}
	(void) fgetstring(file,"Interface type = ");
	(void) fscanf(file,"%d",&interface_type(infc));
	switch (interface_type(infc))
	{
	case PHYSICAL_INTERFACE:
	    if (overlay != YES)
	    {
	        if (gas_params_list(infc) == NULL)
	    	    g_read_print_Gas_param_list(init,io_type,infc);
	        /*g_read_print_Dirichlet_bdry_states(init,io_type,infc); */
	    /*
	     *  NOTE: it will be necessary to compile out
	     *  the next line to restart old version files.
	     */
	        g_read_print_RP_DATA_at_nodes(init,io_type,infc);
	        read_print_stratified_state_function(io_type,infc);
	        g_read_print_ContactWallNodeParams(io_type,infc);
	    }
	    break;
	case EOS_INTERFACE:
	    break;
	}
	return YES;
}		/*end gi_user_read_print_interface*/

EXPORT	void	read_print_stratified_state_function(
	const IO_TYPE *io_type,
	INTERFACE     *infc)
{
	FILE            *file = io_type->file;
	char            fname[256];
	if (next_output_line_containing_string(file,"STRATIFIED STATE FUNCTION")
			== NULL)
	    return;
	
	(void) fgetstring(file,"stratified_state = ");
	(void) fscanf(file,"%s",fname);
	if (strcmp(fname,"isentropic_stratified_state") == 0)
	    g_user_interface(infc)._stratified_state =
					isentropic_stratified_state;
	else if (strcmp(fname,"isothermal_stratified_state") == 0)
	    g_user_interface(infc)._stratified_state =
					isothermal_stratified_state;
	else if (strcmp(fname,"constant_density_stratified_state") == 0)
	    g_user_interface(infc)._stratified_state =
					constant_density_stratified_state;
	else if (strcmp(fname,"NULL") == 0)
	    g_user_interface(infc)._stratified_state = NULL;
}		/*end read_print_stratified_state_function*/

/*ARGSUSED*/
LOCAL	void gi_user_read_print_point(
	POINT	      *point,
	const IO_TYPE *io_type,
	boolean          overlay)
{
	FILE	  *file = io_type->file;
	char	  type[120];

	(void) fgetstring(file,"wave_type = ");
	(void) fscanf(file,"%s",type);
	wave_type(point) = read_wave_type_from_string(type,point->interface);
	if (wave_type(point) < FIRST_PHYSICS_WAVE_TYPE)
	    set_is_bdry(point);
	if (fgetstring(file,"curve->layer_index = "))
	{
	    (void) fscanf(file,"%d",&layer_index(Hyper_surf(point)));
	    num_layers(point->interface) = max(num_layers(point->interface),
			                       layer_index(Hyper_surf(point)));
	}
}		/*end gi_user_read_print_point*/

LOCAL	void	g_set_default_front_parameters(
	INIT_DATA	*init,
	Front		*fr)
{
	int		dim = fr->rect_grid->dim;

	d_set_default_front_parameters(init,fr);
	fr->_copy_front =				g_copy_front;
#if defined(USE_OVERTURE)
        fr->_deep_copy_front =                          g_deep_copy_front;
#endif /* if defined(USE_OVERTURE) */
	fr->_copy_into_front =				g_copy_into_front;

	switch (dim)
	{
	case 1:
	    fr->fr_bdry_untangle = g_bdry_untangle1d;
	    fr->untangle_front = g_untangle_front1d;
	    break;
	case 2:
	    fr->snd_node_propagate = g_snd_node_propagate;
	    fr->node_propagate = g_node_propagate;
	    fr->B_node_bifurcation = g_B_node_bifurcation;
	    fr->untangle_front = g_untangle_interior_curves;
	    fr->fr_vec_bdry_untangle = g_vec_bdry_untangle;
	    fr->init_2drproblem = g_init_2drproblem;
	    fr->identify_physical_node = g_identify_physical_node;
	    fr->_is_correspondence_possible = g_is_correspondence_possible;
	    fr->_fprint_header_for_graph_curve_states = 
	    	g_fprint_header_for_graph_curve_states;
	    fr->phys_set_node_types = set_start_end_status;
	    fr->phys_split_bdry_cross = g_phys_split_bdry_cross;
	    fr->twodrproblem = g_2drproblem;	
	    fr->_fgraph_front_states = g_fgraph_front_states;
	    fr->_check_front_state_consistency = 
		    	g_check_front_state_consistency;
	    fr->_fgraph_curve_states = g_fgraph_curve_states;
	    fr->_check_delete_redundant_node = g_check_delete_redundant_node;
	    fr->_find_i_to_prop_dir = g_find_i_to_prop_dir;
	    fr->_untrack_curve = g_untrack_curve;
	    Init_redistribution_function(fr) = g_init_redistribute;
	    break;
	case 3:
	    fr->_principal_tangent = g_principal_tangent;
	    break;
	}
	fr->transform_state = g_transform_state;
	fr->neumann_bdry_state = g_neumann_bdry_state;
	fr->_read_print_front = g_read_print_front;
	if (debugging("verbose"))
	    fr->print_state = g_verbose_print_state;
	else
	    fr->print_state = g_print_state;
	include_wall_normal_velocity(fr) = NO;
	fr->_fprint_front = g_fprint_front;
}		/*end g_set_default_front_parameters*/

LOCAL	void	gi_prompt_for_front_options(
	INIT_DATA *init,
	Front     *front)
{
	char		s[Gets_BUF_SIZE];
	int		dim = i_intfc(init)->dim;

	curvature_factor(front) = 0.0; /*Off by default*/
	f_prompt_for_front_options(init,front);

	if (dim > 1)
	{
	    screen("To use curvature dependent limiting at scalar fronts\n"
	           "\tenter the curvature factor (dflt = %g): ",
	           curvature_factor(front));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        (void) sscan_float(s,&curvature_factor(front));
	    }
	    if (curvature_factor(front) < 0.0)
	    {
	        screen("ERROR in gi_prompt_for_front_options(), "
	               "negative curvature factor "
		       "%g\n",curvature_factor(front));
	        clean_up(ERROR);
	    }
	}
}		/*end gi_prompt_for_front_options*/

LOCAL	void	gi_read_print_front_options(
	INIT_DATA *init,
	Front     *front)
{
	const IO_TYPE *io_type = restart_io_type(init);
	FILE    *file = io_type->file;

	f_read_print_front_options(init,front);

	if (next_output_line_containing_string(file,
	    "Gas Dynamics Front Parameters") != NULL)
	{
	    if (fgetstring(file,"include_wall_normal_velocity = "))
	        include_wall_normal_velocity(front) = fread_boolean(file);
	    curvature_factor(front) =
	        read_print_float("curvature_factor = ",
		                 curvature_factor(front),io_type);
	}

}		/*end gi_read_print_front_options*/

LOCAL	void	g_read_print_front(
	INIT_DATA	*init,
	Front		*front)
{
	f_read_print_front(init,front);
}		/*end g_read_print_front*/

/*
*			g_copy_front():
*
*	Basic default function for copying a front structure.
*	Allocates storage for the new front and copies the
*	argument into the new structure.
*/

LOCAL	Front *g_copy_front(
	Front		*fr)
{
	Front		*newfr;

	scalar(&newfr,sizeof(G_Front));
	copy_into_front(newfr,fr);
	return newfr;
}		/*end g_copy_front*/


/*
*                       g_deep_copy_front():
*
*       Basic default function for copying a front structure.
*       Allocates storage for the new front and copies the
*       argument into the new structure.
*/

#if defined(USE_OVERTURE)
LOCAL   Front *g_deep_copy_front(
        Front           *fr)
{
        Front           *newfr;

        scalar(&newfr,sizeof(G_Front));
        copy_into_front(newfr,fr);
        scalar(&(newfr->rect_grid), sizeof(RECT_GRID));
        *(newfr->rect_grid) = *(fr->rect_grid);
        scalar(&(newfr->pd_flag),sizeof(Patch_bdry_flag));
        return newfr;
}               /*end g_deep_copy_front*/
#endif /* if defined(USE_OVERTURE) */


/*
*			g_copy_into_front():
*
*	Copies fr into newfr.  Assumes newfr is already allocated.
*/

EXPORT	void g_copy_into_front(
	Front		*newfr,
	Front		*fr)
{
	d_copy_into_front(newfr,fr);
	GasFrontExtension(newfr) = GasFrontExtension(fr);
}		/*end g_copy_into_front*/

LOCAL	void	g_set_default_wave_parameters(
	INIT_DATA *init,
	Wave	  *wave)
{
	int		dim = wave->rect_grid->dim;

	d_set_default_wave_parameters(init,wave);
	wave->_copy_wave = g_copy_wave;
	wave->_copy_into_wave = g_copy_into_wave;
	wave->_free_wave = h_free_wave;
	wave->_assign_wave_parameters = g_assign_wave_parameters;
	wave->_bad_state_data = g_bad_state_data;
	wave->_alloc_MaxWaveSpeed = g_alloc_MaxWaveSpeed;
	wave->_detect_and_load_mix_state = g_detect_and_load_mix_state;
	switch (dim)
	{
	case 1:
	    break;
	case 2:
	    wave->el_integral.linear_cell = g_tri_integral;
	    wave->el_integral.bilinear_cell = g_quad_integral;
	    break;
	case 3:
	    wave->el_integral.linear_cell = NULL;
	    wave->el_integral.bilinear_cell = NULL;
	    break;
	}
	if (debugging("verbose"))
	    wave->print_state = g_verbose_print_state;
	else
	    wave->print_state = g_print_state;
	wave->show_wave_states = g_show_wave_states;

	wave->max_wave_speed = max_speed;
	wave->bundle_states = g_bundle_states;
	wave->unbundle_states = g_unbundle_states;

#if defined(USE_OVERTURE)
        wave->overture_assign_wave_params = g_overture_assign_wave_params;
        wave->overture_assign_wave_st_type = g_overture_assign_wave_st_type;
        wave->overture_to_ft_st = g_overture_to_ft_st;
        wave->ft_to_overture_st = g_ft_to_overture_st;

        wave->trans_wv_st_to_overfunc = g_trans_wv_st_to_overfunc;
        wave->fill_root_extr_overture_cell_st_from_wv =
                g_fill_root_extr_overture_cell_st_from_wv;
        wave->overture_init_interpolation_coarse_to_fine =
                g_overture_init_interpolation_coarse_to_fine;
        wave->overture_interpolation_fine_to_coarse =
                g_overture_interpolation_fine_to_coarse;
        wave->overture_undistribute_interpolation_fine_to_coarse =
                g_overture_undistribute_interpolation_fine_to_coarse;
        wave->scatter_patch_states = g_scatter_patch_states;
        wave->scatter_patch_states_in_sweep_dir =
                g_scatter_patch_states_in_sweep_dir;
        wave->overture_fill_patch_amr_buffer =
                g_overture_fill_patch_amr_buffer2d;
        wave->overture_fill_patch_amr_buffer_pt =
                g_overture_fill_patch_amr_buffer_pt2d;
        wave->overture_injection_after_repatch =
                g_overture_injection_after_repatch;
#endif /* defined(USE_OVERTURE) */
}		/*end g_set_default_wave_parameters*/

LOCAL	void	g_assign_wave_parameters(
	Wave		*newwave,
	Wave		*wave)
{
	d_assign_wave_parameters(newwave,wave);
	g_user_wave(newwave) = g_user_wave(wave);
}		/*end g_assign_wave_parameters*/

LOCAL	Wave	*g_copy_wave(
	Wave		*wave)
{
	Wave		*newwave;

	scalar(&newwave,sizeof(G_Wave));
	copy_into_wave(newwave,wave);
	return newwave;
}		/*end g_copy_wave*/

EXPORT	void	g_copy_into_wave(
	Wave	*newwave,
	Wave	*wave)
{
	d_copy_into_wave(newwave,wave);
	GasWaveExtension(newwave) = GasWaveExtension(wave);
}		/*end g_copy_into_wave*/

/*ARGSUSED*/
LOCAL	void	g_restart_clean_up(
	INIT_DATA *init)
{
	g_free_restart_params_list();
}		/*end g_restart_clean_up*/

LOCAL	void	g_init_redistribute(
	INIT_DATA	*init,
	Front		*front)
{
	f_init_redistribute(init,front);
	Delete_small_loops_function(front) = f_delete_small_loops;
	if (reflect_small_loop_shocks(init) == YES)
	    Delete_small_loops_function(front) = g_delete_small_loops;
}		/*end g_init_redistribute*/

LOCAL	MAX_FRONT_SPEED	*g_alloc_MaxFrontSpeed(
	MAX_FRONT_SPEED	*mxsp,
	INTERFACE	*intfc,
	size_t		sizest)
{
	if (mxsp != NULL)
	    return mxsp;

	mxsp = f_alloc_MaxFrontSpeed(mxsp,intfc,sizest);
	mxsp->operators._initialize = g_initialize_max_front_speed;
	mxsp->operators._set = g_set_max_front_speed;

	return mxsp;
}		/*end g_alloc_MaxFrontSpeed*/

LOCAL	h_MaxWaveSpeed	*g_alloc_MaxWaveSpeed(
	h_MaxWaveSpeed	*mxsp,
	Wave		*wave)
{
	mxsp = h_alloc_MaxWaveSpeed(mxsp,wave);

	mxsp->operators._set = g_set_max_wave_speed;
	mxsp->operators._initialize = g_initialize_max_wave_speed;

	return mxsp;
}		/*end g_alloc_MaxWaveSpeed*/

EXPORT	void gi_user_read_print_surface(
	SURFACE	      *surf,
	const IO_TYPE *io_type,
	boolean          overlay)
{
	HYPER_SURF    *hs;
	f_user_read_print_surface(surf,io_type,overlay);

	hs = Hyper_surf(surf);
	surface_tension(hs) = read_print_float(
	                  "surface_tension(hs) = ",0.0,io_type);

	if (wave_type(surf) == NEUMANN_BOUNDARY)
	    if (fgetstring(io_type->file,"no slip = "))
	    {
		no_slip(Hyper_surf(surf)) = fread_boolean(io_type->file);
		if (no_slip(Hyper_surf(surf)) == YES)
		    adherence_coeff(Hyper_surf(surf)) = read_print_float(
			"adherence coefficient = ",1.0,io_type);
	    }
}		/*end gi_user_read_print_surface*/



