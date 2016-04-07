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
*				dinit.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the driver initializing routines, as called in main,
*	and the principle or default initializing subroutines.
*/


#include <driver/ddecs.h>

	/* LOCAL Function Declarations */
LOCAL	int	to_parallel_restart(RECT_GRID*);
LOCAL	void	d_init_pause_criteria(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_hyp_and_top_grid(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_interior_states(Wave*,Front*,INIT_DATA*,
				     void (*)(double*,COMPONENT,
					      Locstate,INTERFACE*,
					      INIT_DATA*));
LOCAL   void    init_max_time_step_mods(INIT_DATA*);
LOCAL	void	init_pp_grid(INIT_DATA*);
LOCAL	void	init_restart_amr(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_restart_interior_states(INIT_PHYSICS*,CHART*,
					     INPUT_SOLN**,INIT_DATA*);
LOCAL	void	init_root_chart(CHART*);
LOCAL	void	init_spatial_grids(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_states(INIT_DATA*,INIT_PHYSICS*,CHART*,INPUT_SOLN**,
			    RESTART_DATA*);
LOCAL	void	init_time_step_limit(INIT_DATA*);
LOCAL	void	no_bc_propagate(Grid*);
LOCAL	void	read_input(int*,char***,INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	set_up_cauchy_data(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	d_set_pp_grid(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	start_up(int*,char***,INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	usage(const char*);

LOCAL	char	s[Gets_BUF_SIZE];	/*Scratch array for line input*/

EXPORT	void	set_driver_hooks(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	D_INIT_DATA	*d_init = d_init_data(init);

	set_hyp_hooks(init);
	ip->root->grid->_print_Grid_structure = d_print_Grid_structure;
	ip->root->grid->_print_remap_values = i_print_remap_values;
	ip->prt->_init_printplot = d_init_printplot;
	ip->init_run = d_init_run;
	ip->_init_remap_and_rect_grid = i_init_remap_and_rect_grid;
	ip->_init_stopping_criteria = d_init_stopping_criteria;
	ip->_init_pause_criteria = d_init_pause_criteria;

	d_init->_prompt_for_printing_options = d_prompt_for_printing_options;
	d_init->_prompt_for_initial_intfc_options =
	    d_prompt_for_initial_intfc_options;
	d_init->_restart_clean_up = d_restart_clean_up;
	d_init_data(init)->_set_interface_hooks = d_set_interface_hooks;
}		/*end set_driver_hooks*/

/*
*				d_init_run():
*
*	This is the main control loop for initialization.
*/


EXPORT	void	d_init_run(
	int		*pargc,
	char		***pargv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	debug_print("init","Entered d_init_run()\n");

	read_input(pargc,pargv,init,ip);
	set_up_cauchy_data(init,ip);

	if (debugging("restart"))
	{
	    (void) printf("Interface after d_init_run()\n");
	    print_interface(ip->root->front->interf);
	}
	if (debugging("wave_states"))
	{
	    Wave *wave = ip->root->wave;
	    (void) printf("Wave states at end of d_init_run()\n");
	    (*wave->show_wave_states)(wave);
	}
	debug_print("init","Left d_init_run()\n");
}		/*end d_init_run*/

/*
*			read_input():
*
*	Reads command line arguments and input file to record user defined
*	settings for a given run.   All actions in this function and its
*	subroutines are minimal,  simply make a record of the user's
*	input specifications.
*/

LOCAL	void	read_input(
	int		*pargc,
	char		***pargv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
		/* Initialize the Grid and Jacobian */
	start_up(pargc,pargv,init,ip);
	init_spatial_grids(init,ip);
	init_time_step_limit(init);
	init_stopping_criteria(init,ip);
	init_pause_criteria(init,ip);
        init_max_time_step_mods(init);
	init_triangulation_options(&Tri_grid_hooks(init),Comp_grid(init).dim);
	init_wave_mv_state_list(init);
	init_pp_grid(init);
	prompt_for_printing_options(init);
	prompt_for_physics_options(init,ip);
	prompt_for_initial_intfc_options(init);
	prompt_for_front_options(init,ip->root->front);
}		/*end read_input*/

/*
*			set_up_cauchy_data():
*
*	Uses the input specifications recorded in the INIT_DATA structure
*	to set up the "first" time step Cauchy data for a run, as well as
*	setting up all diagnostic data structures.
*/

LOCAL	void	set_up_cauchy_data(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	CHART		*root = ip->root;
	Printplot	*prt = ip->prt;
	Front		*front = root->front;
	Grid		*grid = root->grid;
	RECT_GRID	*top_grid = &ip->top_grid;
	RESTART_DATA	*restart = &Restart_data(init);
	Wave		*wave = root->wave;
	static RECT_GRID comp_grid;

#if defined(USE_OVERTURE)
        read_overture_option(init,ip);
#endif /* if defined(USE_OVERTURE) */

	grid->init = init;
	grid->rect_grid = wave->rect_grid = front->rect_grid = &comp_grid;
	init_hyp_and_top_grid(init,ip);
	wave->Tri_grid_hooks = Tri_grid_hooks(init);
	alloc_wave_mv_state_list(init,wave);

		/* Initialize PP_Grid for Domain Decomposition */

	d_set_pp_grid(init,ip);

		/* Initialize Printing and Plotting */

	init_printplot(init,grid,front,prt);

		/*REDO INITIALIZATION BELOW*/

	Clear_redistribution_parameters(front);
	Init_redistribution_function(front) = f_init_redistribute;

	root->bc_propagate = no_bc_propagate;
	root->ellip = NULL;
	root->parab = NULL;

		/* Initialze Physics-Dependent Variables */

	init_physics(init,ip);
	MaxWaveSpeed(wave) = alloc_MaxWaveSpeed(MaxWaveSpeed(wave),wave);

	switch (front->rect_grid->dim)
	{
	case 1:
	    wave->_nearest_crossing = nearest_crossing1d;
	    break;
	case 2:
	    wave->_nearest_crossing = nearest_crossing2d;
	    break;
	case 3:
	    wave->_nearest_crossing = nearest_crossing3d;
	    break;
	}

	if (ip->init_hyp_bdry)
	    (*ip->init_hyp_bdry)(front,ip->init_bdry_state);


		/* Initialize Point Sources */

	if (ip->init_point_sources)
	    (*ip->init_point_sources)(wave);

		/* Make Interface or Input It from File or Screen */

	if (got_intfc_from_file(init) == YES)
	{
	    int		i;
	    OUTPUT_DATA	**mod = prt->main_output_data;

	    /* The next print time(s) may need to be reset because
	     * the current step and time were not known until now. */

	    ++(grid->step);
	    if (!debugging("prt_rst"))
	    {
	        for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
	        {
	            if (output_format_on(i,prt))
	                set_next_print_time(grid,Print_control(mod[i]));
	        }
	    }
#if defined(USE_HDF)
	    set_hdf_append(prt);
#endif /* defined(USE_HDF) */
	}

	if (restart_io_type(init) != NULL)
	{
	    restart_multi_data(init) =
		to_parallel_restart(&Restart_comp_grid(init));
	    if (restart_multi_data(init) == YES)
	    {
	    	zoom_rect_grid(&comp_grid,&Restart_comp_grid(init));
	    	zoom_rect_grid(top_grid,&Restart_top_grid(init));
	    	(void) adjust_top_grid_for_square(top_grid,&comp_grid);
	    	copy_rect_grid(&grid->pp_grid->Zoom_grid,&comp_grid);
	    }
	}

	set_size_of_intfc_state(front->sizest);
	front->interf = (restart_intfc(init) != NULL) ?
	    restart_intfc(init) : make_interface(grid->rect_grid->dim);

	set_topological_grid(front->interf,top_grid);
	set_computational_grid(front->interf,&comp_grid);

		/* Initialize Front Structure */
	init_front(init,front);
	/*DEBUG_TMP printf("#redis_flag = %d\n", front->redis_flag); */

		/* Perform Problem-Dependent Initializations */
	if (ip->init_interface)
	    (*ip->init_interface)(init,ip);
	
	null_sides_are_consistent();
	/*DEBUG_TMP check_print_intfc("After init interface", "init_jet", 'f',  */
	       /*DEBUG_TMP front->interf, 1, -1, NO); */

	/*print_rectangular_grid(front->rect_grid); */

	if (comp_grid.dim == 2)
	    test_for_mono_comp_curves(front->interf);


		/* Now that interface properly made, reset */
		/* got_intfc_.. for state initialization   */
		/* if interface came from stdin.	   */

	if ((got_intfc_from_file(init)==YES) &&
	    (restart_io_type(init)->file==stdin))
	{
	    got_intfc_from_file(init) = NO;
	    restart_io_type(init) = NULL;
	}

		/* Redistribute the initial front */

	initial_front_redistribute(front,restart_io_type(init));
	interpolate_intfc_states(front->interf) = YES;
	
	/*print_rectangular_grid(front->rect_grid); */
	null_sides_are_consistent();
	/*DEBUG_TMP check_print_intfc("After init redist", "init_jet", 'f',  */
	       /*DEBUG_TMP front->interf, 1, -1, YES); */

		/* Initialize Diffusion Effects */

	if (ip->init_parabolic)
	    (*ip->init_parabolic)(front);

		/* Initialize State Variables */

	init_states(init,ip,ip->root,prt->restart_soln,restart);
	
	delete_untracked_hyper_surfaces(front,wave);
        
	/* flexible cauchy data, overwriting states set by init_states() */
	/*if (ip->cauchy_deposition && restart_io_type(init) == NULL) */
	/*if (ip->cauchy_deposition) */
	/*    (*ip->cauchy_deposition)(wave,front); */

	if ((got_intfc_from_file(init) == YES) &&
	    (pause(grid)->_await != NULL))
	{
	    await(grid,front,wave,prt,got_intfc_from_file(init));
	}

		/* Initialize The Root Chart */

	init_root_chart(root);

		/* Initialize Automatic Mesh Refinement */

#if !defined(USE_OVERTURE)
	if (got_intfc_from_file(init) == YES)
	    init_restart_amr(init,ip);
	else 
	    init_amr(root);
#endif /* if !defined(USE_OVERTURE) */

		/* Initialization and Zeroth Step of Elliptic Solver */

	start_clock("ELLIP");
	if (root->ellip)
	{
	    (*root->ellip)(init,root,restart_io_type(init));

	    screen("Request ellip solution each time step? (y(dflt),n): ");
	    (void) Gets(s);
	    if ((s[0] == 'n') || (s[0] == 'N')) root->ellip = NULL;
	}
	stop_clock("ELLIP");

	init_statistics(front,grid,wave,prt,init);

	if (got_intfc_from_file(init) == YES)
	{
	    (void) Fclose(restart_io_type(init)->file);
	    restart_io_type(init) = NULL;
	    grid->dt = nonphysics_timestep_reductions(grid,wave,front,
						      prt,grid->dt);
	}

	screen("\n\n\t\t--- End of Input ---\n\n");
}		/*end set_up_cauchy_data*/


/*
*			no_bc_propagate():
*
*	Used as a default for root->bc_propagate. 
*/

/*ARGSUSED*/
LOCAL void no_bc_propagate(Grid* grid)
{
}		/*end no_bc_propagate*/


LOCAL void init_root_chart(CHART *root)
{
	int		i, j, dim;

	/*
	* Assumes root->front, root->wave, root->hyp_solver,
	* root->bc_propagate, root->parab,  and root->ellip,  are already set.
	*/

	root->parent = NULL;
	root->prev_chart = NULL;
	root->next_chart = NULL;

	scalar(&root->to_root,sizeof(AFLIN));
	scalar(&root->from_root,sizeof(AFLIN));
	dim = root->front->interf->dim;
	for (i = 0; i < dim; ++i)
	{
		root->to_root->b[i] = 0.;	root->from_root->b[i] = 0.;
		for (j = 0; j < dim; ++j)
		{
			root->to_root->a[i][j]   = (i == j) ? 1.0 : 0.0;
			root->from_root->a[i][j] = (i == j) ? 1.0 : 0.0;
		}
	}

	root->level = NULL;
	root->dynamic = NO;
	root->is_old_chart = NO; /* Temporary */

	root->newfront = copy_front(root->front);
	root->newfront->interf = NULL;
	root->newwave = copy_wave(root->wave);
	clear_wave_pointers(root->newwave);
}		/*end init_root_chart*/



LOCAL void init_restart_amr(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	CHART		*base_root = ip->root;
	CHART		*chart, *chartp, *root;
	const IO_TYPE   *io_type = restart_io_type(init);
	Front		*front = base_root->front;
	INTERFACE	*chart_restart_intfc;
	INPUT_SOLN	**restart_soln;
	INPUT_SOLN	**root_rs = ip->prt->restart_soln;
	LEVEL		*level,*levelp;
	RESTART_DATA	*restart = &Restart_data(init);
	RESTART_DATA	restart_at_level;
	RECT_GRID	*rs_fgr = &restart_at_level._comp_grid;
	RECT_GRID	*rs_tgr = &restart_at_level._top_grid;
	Wave		*wave = base_root->wave;
	int		var,last;
	int		i, j, dim = front->interf->dim;
	int		n_restart_vars = ip->prt->n_restart_vars;
	static CHART	Chart;

	debug_print("chart","Entered init_restart_amr()\n"); 

	chart_of_front(front) = chart_of_wave(wave) = ip->root;
	screen("\nType 'y' to have automatic mesh refinement: ");
	(void) Gets(s);
	if ((s[0] != 'y') && (s[0] != 'Y'))
	{
	    ip->root->level = NULL;
	    debug_print("chart","Left init_restart_amr()\n"); 
	    return;
	}
	
	restart_at_level = *restart;
	chart_of_front(front) = chart_of_wave(wave) = &Chart;
	root = &Chart;
	base_root->hyp_solver = hyp_amr;
	root->parent = NULL;
	root->dynamic = NO;
	root->hyp_solver = base_root->hyp_solver;
	root->bc_propagate = base_root->bc_propagate;
	root->parab = base_root->parab;
	root->ellip = base_root->ellip;
	root->front = front;
	root->wave = wave;
	scalar(&root->to_root,sizeof(AFLIN));
	scalar(&root->from_root,sizeof(AFLIN));
	for (i = 0; i < dim; ++i)
	{
	    root->to_root->b[i] = 0.;	root->from_root->b[i] = 0.;
	    for (j = 0; j < dim; ++j)
	    {
	    	root->to_root->a[i][j]   = (i == j) ? 1.0 : 0.0;
	    	root->from_root->a[i][j] = (i == j) ? 1.0 : 0.0;
	    }
	}
	root->newfront = copy_front(front);
	root->newfront->interf = NULL;
	root->newwave = copy_wave(wave);
	clear_wave_pointers(root->newwave);

		/* initialize LEVEL 0 */

	root->level = read_print_LEVEL(io_type);
	root->level->first = root;
	root->level->last = root;
	root->level->next_coarser_level = NULL;
	levelp = root->level;

	/* loop on levels and charts and read them in */

	bi_array(&restart_soln,n_restart_vars,1,sizeof(INPUT_SOLN));
	while ((level = read_print_LEVEL(io_type)))
	{
	    level->next_coarser_level = levelp;
	    chartp = (CHART *) NULL;
	    if (debugging("chart"))
		    print_LEVEL(stdout,level);
	    while ((chart = read_print_CHART(io_type,root,&last)))
	    {
		chart->hyp_solver = root->hyp_solver;
		chart->bc_propagate = root->bc_propagate;
		chart->parab = root->parab;
		chart->ellip = root->ellip;
		chart->level = level;
		chart->prev_chart = chartp;
		if (chartp == NULL)
			level->first = chart;
		else
			chartp->next_chart = chart;

		set_size_of_intfc_state(front->sizest);
		chart_restart_intfc =
		    read_print_intfc_and_grids(init,restart_io_type(init),
	                                       NULL,rs_tgr,rs_fgr);
		chart->front->interf = chart_restart_intfc;
		copy_rect_grid(chart->front->rect_grid,rs_fgr);

		set_topological_grid(chart->front->interf,rs_tgr);
		set_computational_grid(chart->front->interf,
				       chart->front->rect_grid);

		for (var = 0; var < n_restart_vars; ++var)
		{
		    restart_soln[var]->name = root_rs[var]->name;
		    restart_soln[var]->fit = root_rs[var]->fit;
		    restart_soln[var]->smoothness = root_rs[var]->smoothness;
		    restart_soln[var]->set_intfc_states =
			    			root_rs[var]->set_intfc_states;
	        }

		init_states(init,ip,chart,restart_soln,&restart_at_level);

		chart->front->_hyp_solution = amr_hyp_solution;
		chart->front->_hyp_grad_solution = amr_hyp_grad_solution;

		chartp = chart;
		if (debugging("chart"))
		    print_CHART(stdout,chart);
		if (last == YES)
		    break;
	    }
	    level->last = chart;
	    chartp->next_chart = NULL;
	    if (levelp != NULL)
		levelp->next_finer_level = level;
	    levelp = level;
	    if (chart->level->level_num == *chart->level->max_num_levels)  
		break;
	}
	free(restart_soln);
	level->next_finer_level = (LEVEL *) NULL;

	debug_print("chart","Left init_restart_amr()\n"); 
}		/*end init_restart_amr*/

/*ARGSUSED*/
EXPORT	void	d_init_stopping_criteria(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	STOP	   *stop;
	char	   stop_time_mode[80];
	char       mesg[1024];

	static const char *fmt = "%lf %d %lf %s";

	if (initial_stop(init) == NULL)
	{
	    scalar(&stop,sizeof(STOP));
	    stop->_print_stopping_criteria = d_print_stopping_criteria;
	    stop->_stop_run = d_stop_run;
	    initial_stop(init) = stop;
	}
	else
	    stop = initial_stop(init);

	stop->_stop_time_mode = CONSTANT_TIME;
	(void) strcpy(stop_time_mode,"constant");
	stop->_stop_time = HUGE_VAL;
	stop->_stop_step = INT_MAX;
	screen("\t\tGeneral run termination/pause conditions\n\n");
	(void) sprintf(mesg,"Enter limits on real time (max time), mesh time "
			    "(max timesteps), an optional initial time, and an "
			    "optional stop time mode (exact or constant), "
			    "(dflt = %g %d %g %s)",
			    stop->_stop_time,stop->_stop_step,
			    initial_time(init),stop_time_mode);
	screen_print_long_string(mesg);
	screen(": ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    int nr;
	    nr = sscanf(s,fmt,&stop->_stop_time,&stop->_stop_step,
			&initial_time(init),stop_time_mode);
	    if (nr == 4)
	    {
		if (strcasecmp(stop_time_mode,"exact") == 0)
		    stop->_stop_time_mode = EXACT_TIME;
		if (strcasecmp(stop_time_mode,"constant") == 0)
		    stop->_stop_time_mode = CONSTANT_TIME;
	    }
	}
}		/*end d_init_stopping_criteria*/

/*ARGSUSED*/
LOCAL	void	d_init_pause_criteria(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	PAUSE	*pause;
	if (initial_pause(init) == NULL)
	{
	    scalar(&pause,sizeof(PAUSE));
	    pause->_await = d_await;
	    initial_pause(init) = pause;
	}
	else
	    pause = initial_pause(init);
	pause->_pause_mode = CONSTANT_TIME;
	screen("Specify the pause time mode [exact%s, constant%s, mesh%s]: ",
	       (pause->_pause_mode == EXACT_TIME)    ? "(dflt)" : "",
	       (pause->_pause_mode == CONSTANT_TIME) ? "(dflt)" : "",
	       (pause->_pause_mode == MESH_TIME)     ? "(dflt)" : "");
	(void) Gets(s);
	if (s[0] == 'm' || s[0] == 'M')
	    pause->_pause_mode = MESH_TIME;
	else if (s[0] == 'c' || s[0] == 'C')
	    pause->_pause_mode = CONSTANT_TIME;
	else if (s[0] == 'e' || s[0] == 'E')
	    pause->_pause_mode = EXACT_TIME;
	
	if (mesh_time_output(pause->_pause_mode))
	{
	    pause->_pause_step = INT_MAX;
	    screen("Enter the first Pause Time Step (dflt = %d): ",
		   pause->_pause_step);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d \n",&pause->_pause_step);
	}
	else
	{
	    pause->_pause_time = HUGE_VAL;
	    screen("Enter the first Pause Time (dflt = %g): ",
	    	   pause->_pause_time);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscan_float(s,&pause->_pause_time);
	}
}		/*end d_init_pause_criteria*/

/*
*			init_hyp_and_top_grid():
*
*	Initializes all of the paramaters in the driver's grid structure
*	except for the statistical variables, which are initialized in
*	(*ip->_init_physics)().
*/

LOCAL void init_hyp_and_top_grid(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Grid		*grid = ip->root->grid;
	RECT_GRID	*comp_grid = grid->rect_grid;
	RECT_GRID	*top_grid = &ip->top_grid;

	debug_print("init","Entered init_hyp_and_top_grid()\n");

	copy_rect_grid(comp_grid,&Comp_grid(init));
	copy_rect_grid(top_grid,&Top_grid(init));
	dt_lim(grid) = initial_dt_lim(init);
	stop(grid) = initial_stop(init);
	pause(grid) = initial_pause(init);

	if (got_intfc_from_file(init) == YES)
	{
	    grid->step = restart_time_step(init);
	    grid->time = restart_time(init);
	    grid->dt_last = restart_dt_last(init);
	    grid->dt = restart_next_dt(init);
	}
	else
	{
	    grid->time = initial_time(init);
	    grid->dt = ip->initial_dt;		/* Initial time increment */
	    grid->dt_last = grid->dt;
	    grid->step = 0;			 /* Initial time step */
	}

	max_num_time_step_mods(grid) = initial_max_num_time_step_mods(init);
	if (debugging("Grid")) 
	{
	    print_Grid_structure(grid);
	    (void) printf("\nInput Topological Grid\n");
	    print_rectangular_grid(top_grid);
	}

	debug_print("init","Left init_hyp_and_top_grid()\n");
}		/*end init_hyp_and_top_grid*/


/*
*			init_states():
*
*	Initializes state pointers and values in the tracked case.
*/

void    tecplot_interface_in_ball(const char*, INTERFACE*);

LOCAL void init_states(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	CHART		*chart,
	INPUT_SOLN	**restart_soln,
	RESTART_DATA	*restart)
{
	const IO_TYPE   *io_type = restart_io_type(init);
	Front		*front = chart->front;
	RECT_GRID	*rs_fgr = &restart->_comp_grid;
	RECT_GRID	*rs_tgr = &restart->_top_grid;
	RECT_GRID	fgr, tgr;
	Wave		*wave = chart->wave;
	boolean		got_intfc_from_file = got_intfc_from_file(init);
	int		n_restart_vars = ip->prt->n_restart_vars;
	int		var;
	int             *iperm,idir,dim;

	debug_print("init","Entered init_states()\n");

		/* Save copies of rectangular grids */

	copy_rect_grid(&fgr,front->rect_grid);
	copy_rect_grid(&tgr,&topological_grid(front->interf));

	if (ip->init_cauchy_data_pointers != NULL)
	    (*ip->init_cauchy_data_pointers)(ip,got_intfc_from_file);

	dim = front->rect_grid->dim;
	if (got_intfc_from_file == YES)
	{
	    FILE *restart_file = io_type->file;
	    int	i, c, dim = fgr.dim;

	    if (front->sizest > 0)
	    {
	        if (read_state_variables(io_type,n_restart_vars,
			front->interf,restart_soln,dim) == FUNCTION_FAILED)
	        {
		    screen("ERROR in init_states(), "
		           "read_state_variables() failed\n");
		    clean_up(ERROR);
	        }
	    }

	    /* Use restart grids for initializing front and wave */

	    /*DEBUG_TMP printf("#tst res grid\n"); */
	    /*DEBUG_TMP print_rectangular_grid(rs_fgr); */
	    /*DEBUG_TMP print_rectangular_grid(rs_tgr); */
	    /*DEBUG_TMP printf("#tst org grid\n"); */
	    /*DEBUG_TMP print_rectangular_grid(&fgr); */
	    /*DEBUG_TMP print_rectangular_grid(&tgr); */

	    copy_rect_grid(front->rect_grid,rs_fgr);
	    set_topological_grid(front->interf,rs_tgr);
	    set_computational_grid(front->interf,front->rect_grid);

	    init_front_states(front,init,ip->restart_intfc_initializer);
	    if (!restart_multi_data(init))
	    {
		clip_front_to_subdomain(front);
		copy_rect_grid(rs_fgr,front->rect_grid);
		copy_rect_grid(&fgr,front->rect_grid);
		copy_rect_grid(rs_tgr,&topological_grid(front->interf));
		copy_rect_grid(&tgr,&topological_grid(front->interf));
	    }

            if (tracking_algorithm(init) == GRID_BASED_TRACKING  ||
	        tracking_algorithm(init) == THREE_COMP_GRID_BASED_TRACKING)
	        interface_reconstructed(front->interf) = YES;
	

	    init_restart_interior_states(ip,chart,restart_soln,init);

	    initialize_max_wave_speed(wave);
	    if (next_output_line_containing_string(restart_file,"WAVE SPEEDS:"))
	    {
	        if (!read_print_max_wave_speed_info(init,io_type,wave))
	        {
	    	    /*Old style printout */
		    if (fgetstring(restart_file,"wave->maxsp ="))
		    {
			if ((c = getc(restart_file)) != '\f') /* NOBINARY */
			{
			    (void) ungetc(c,restart_file);
			    for (i = 0; i < dim; ++i)
			    	(void) fscan_float(restart_file,Maxsp(wave)+i);
			}
			else
			{		/* BINARY */
			    (void) getc(restart_file);
			    (void) read_binary_real_array(Maxsp(wave),dim,
			                                  io_type);
			}
		    }
		}
	    }
	    restart_clean_up(init);
	}
	else
	{
	    init_front_states(front,init,ip->intfc_initializer);

	    /* parallel part for front */
	
	    clip_front_to_subdomain(front);
	
	    /*DEBUG_TMP check_print_intfc("After clip_front_to_subdomain in init_states",  */
	    		/*DEBUG_TMP "init_st",'g',front->interf,0,-1,NO); */
	   
            /* init the amr grid */
#if defined(USE_OVERTURE)
            if(chart->overparam != NULL)
                overture_init_amr(chart);
#endif /* if defined(USE_OVERTURE) */
	    
	    init_interior_states(wave,front,init,ip->initializer);

            if(ip->stratify_density_for_gravity)
                while((*ip->stratify_density_for_gravity)())
                {
                        init_front_states(front,init,ip->intfc_initializer);
                        init_interior_states(wave,front,init,ip->initializer);
                }
	}

	if (debugging("wave_states"))
	{
	    (void) printf("Wave states in init_states() "
	                  "before scatter_states()\n");
	    (*wave->show_wave_states)(wave);
	}

	iperm = set_iperm(0,dim);
	for (idir = 0; idir < dim; idir++)
        {
	    if (scatter_states(wave,front,iperm,idir) != FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in init_states(), scatter_states() failed\n");
	        clean_up(ERROR);
	    }
	}

	init_point_source_composition(wave,ip->pt_source_initializer);

	if (got_intfc_from_file == YES) 
	{
	    int i, regrid = NO;

	    for (i = 0; i < front->interf->dim; ++i)
	    {
	    	if (rs_fgr->gmax[i] != fgr.gmax[i])
	    	{
	    	    regrid = YES;
		    break;
		}
	    }
	   
	    /*DEBUG_TMP printf("#regrid %d\n", regrid); */

	    if (regrid)
	    {
	    	if (!regrid_front_and_wave(front,wave,&fgr,&tgr,YES))
	    	{
	    	    screen("ERROR in init_states()\n"
		           "regrid_front_and_wave() failed\n");
		    clean_up(ERROR);
		}
	    }
	    else
	    {
	    	copy_rect_grid(front->rect_grid,&fgr);
	    	set_topological_grid(front->interf,&tgr);
	    }

	    for (var = 0; var < n_restart_vars; ++var)
	    	free_input_soln(restart_soln[var]);

	    if (ip->set_Dirichlet_boundary_states)
	    	(*ip->set_Dirichlet_boundary_states)(front);
	}

	init_mv_list_states(wave,front,init,ip->initializer);

	if (debugging("wave_states"))
	{
	    (void) printf("Wave states at end of init_states()\n");
	    (*wave->show_wave_states)(wave);
	}

	debug_print("init","Left init_states()\n");
}		/*end init_states*/


LOCAL	void init_restart_interior_states(
	INIT_PHYSICS	*ip,
	CHART		*chart,
	INPUT_SOLN	**restart_soln,
	INIT_DATA	*init)
{
	COMPONENT	comp;
	Front		*front = chart->front;
	Locstate	state;
	Wave		*wave = chart->wave;
	int		n_vars = ip->prt->n_restart_vars;
	int		icoords[MAXD], Trans[MAXD], abs_icoords[MAXD];
	int		*gmax = wave->rect_grid->gmax;
	int		*Gmax = restart_soln[0]->grid.gmax;
	int		i, dim = front->interf->dim;
	int		status;

	debug_print("init_restart","Entered init_restart_interior_states()\n");

	if (wave->sizest == 0 || ip->restart_initializer == NULL)
	{
	    debug_print("init_restart","Left init_restart_interior_states()\n");
	    return;
	}

	wave->old_wave = NULL;
	
	status = init_hyp_solution_function(wave,front);
	
	if (status != GOOD_STEP)
	{
	    screen("ERROR: init_hyp_solution_function() failed\n");
	    clean_up(ERROR);
	}

	/* define translational constants */
	if (!restart_multi_data(init))
	{
	    double   rel_cnr, abs_cnr, h;

	    for (i = 0; i < dim; ++i)
	    {
	    	rel_cnr = front->rect_grid->L[i];
	    	abs_cnr = front->pp_grid->dom[i][0];
	    	h = front->rect_grid->h[i];
 
	    	Trans[i] = irint((rel_cnr - abs_cnr)/h);
	    }
	}
	else
	{
	    for (i = 0; i < dim; ++i) Trans[i] = 0;
	}

#define set_icoords_index(__i__,__j__,MAX)				\
	{								\
	    double tmp;							\
	    icoords[(__i__)] = (__j__);					\
	    tmp = abs_icoords[(__i__)] = (__j__) + Trans[(__i__)];	\
	    if (tmp < 0 || tmp > (MAX) - 1) continue;			\
	}

	switch (dim)
	{
	case 1:
	{
	    int		ix;
	    int		xmax;
	    int		Xmax;

	    xmax = gmax[0];	Xmax = Gmax[0];
	    for (ix = 0; ix < xmax; ++ix) 
	    {
	        set_icoords_index(0,ix,Xmax)

	        state = Rect_state(icoords,wave);
	        comp = Rect_comp(icoords,wave);
	        (*ip->restart_initializer)(abs_icoords,n_vars,comp,
	    			           restart_soln,state,init);
	    }
	    break;
	}
	case 2:
	{
	    int		ix, iy;
	    int		xmax, ymax;
	    int		Xmax, Ymax;

	    xmax = gmax[0];	Xmax = Gmax[0];
	    ymax = gmax[1];	Ymax = Gmax[1];
	    for (iy = 0; iy < ymax; ++iy)
	    {
	        set_icoords_index(1,iy,Ymax)
	        for (ix = 0; ix < xmax; ++ix) 
	        {
	            set_icoords_index(0,ix,Xmax)
	            state = Rect_state(icoords,wave);
	            comp = Rect_comp(icoords,wave);
	            (*ip->restart_initializer)(abs_icoords,n_vars,comp,
	    				       restart_soln,state,init);
		}
	    }
	    break;
	}
	case 3:
	{
	    int		ix, iy, iz;
	    int		Xmax, Ymax, Zmax;
	    int		imax[3], imin[3];
	    RECT_GRID	*rgr = wave->rect_grid;

	    for(i=0; i<3; i++)
	    {
		imin[i] = 0;
		imax[i] = rgr->gmax[i];
		
		if(rect_boundary_type(front->interf,i,0) == OPEN_BOUNDARY)
		{
		    Trans[i] = rgr->lbuf[i];
		    imin[i] = -rgr->lbuf[i];
		}
		if(rect_boundary_type(front->interf,i,1) == OPEN_BOUNDARY)
		{
		    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
		}
	    }

	    /*DEBUG_TMP printf("#interior st\n"); */
	    /*DEBUG_TMP print_rectangular_grid(rgr); */
	    /*DEBUG_TMP print_rectangular_grid(&(restart_soln[0]->grid)); */
	    /*DEBUG_TMP print_int_vector("Trans", Trans, 3, "\n"); */
	    /*DEBUG_TMP print_int_vector("imin", imin, 3, "\n"); */
	    /*DEBUG_TMP print_int_vector("imax", imax, 3, "\n"); */
	    
	    Xmax = Gmax[0];
	    Ymax = Gmax[1];
	    Zmax = Gmax[2];
	    
	    for (iz = imin[2]; iz < imax[2]; ++iz)
	    {
	        set_icoords_index(2,iz,Zmax)
	        for (iy = imin[1]; iy < imax[1]; ++iy)
		{
		    set_icoords_index(1,iy,Ymax)
		    for (ix = imin[0]; ix < imax[0]; ++ix)
		    {
		        set_icoords_index(0,ix,Xmax)
		        state = Rect_state(icoords,wave);
		        comp = Rect_comp(icoords,wave);
		        (*ip->restart_initializer)(abs_icoords,n_vars,comp,
						   restart_soln,state,init);
		    }
		}
	    }
	    break;
	}
	}

#undef set_icoords_index

	debug_print("init_restart","Left init_restart_interior_states()\n");
}		/*end init_restart_interior_states*/

/*
*			d_restart_clean_up():
*/

/*ARGSUSED*/
EXPORT	void	d_restart_clean_up(
	INIT_DATA *init)
{
}		/*end d_restart_clean_up*/

/*
*		d_prompt_for_initial_intfc_options():
*
*	Prompts for whether or not to input an interface from a file
*	(or the screen).
*/

EXPORT void d_prompt_for_initial_intfc_options(
	INIT_DATA	*init)
{
	FILE    *file;
	char    bname[256], filename[256];

	debug_print("init","Entered init_restart_interior_states()\n");
	/*Initialize restart data structure*/
	zero_scalar(&Restart_data(init),sizeof(RESTART_DATA));
	got_intfc_from_file(init) = NO;
	restart_multi_data(init) = NO;
	restart_time_step(init)	= INT_MIN;
	restart_time(init)	= -HUGE_VAL;
	restart_dt_last(init) = -HUGE_VAL;
	restart_next_dt(init) = -HUGE_VAL;

	screen("\nSpecify initial interface of tracked curves\n");
	screen("Choices are\n");
	screen("\tInput interface by hand (type `screen')\n");
	screen("\tInput interface from a file ");
	screen("(restart option - enter filename)\n");
	screen("\tRequest default option(s) (hit `return')\n");
	screen("Enter choice: ");
	(void) Gets(s);

	if (s[0] == '\0')
	{
	    debug_print("init","Left init_restart_interior_states()\n");
	    return;
	}

	(void) sscanf(s,"%s",bname);

			/* Prompt for Interface */

	if (strcmp(bname,"screen") == 0)
	{
	    restart_intfc(init) = read_interface();
	    debug_print("init","Left init_restart_interior_states()\n");
	    return;
	}

			/* Read Interface from File */

	screen("Enter the time step at which to find the interface data: ");
	(void) Scanf("%d\n",&restart_time_step(init));

	if (strstr(bname,"lastdump") != NULL)
	    set_output_file_name(pp_mynode(),filename,bname,-1,0);
	else
	    set_output_file_name(pp_mynode(),filename,bname,
	    		restart_time_step(init),0);

	if ((file = fopen(filename,"r")) == NULL)
	{
	    /* Old style printout ? */
	    (void) strcpy(filename,bname);
	    if (pp_numnodes() > 1)
	    	(void) sprintf(filename,"%s.%d",filename,pp_mynode());
	    if ((file = fopen(filename,"r")) == NULL)
	    {
	    	screen("ERROR: unable to open file '%s'\n",filename);
	    	clean_up(ERROR);
	    }
	}
	determine_io_type(file,&restart_IO_type_store(init));
	restart_io_type(init) = &restart_IO_type_store(init);

	if (pp_numnodes() > 1)
	{
	    screen("Input data file on processor %d is %s\n",
		   pp_mynode(),filename);
	}

	determine_read_version(file);
	if (!position_file_at_read_time(file,
				       restart_time_step(init),
				       &restart_time(init),
				       &restart_dt_last(init)))
	{
	    screen("ERROR in d_prompt_for_initial_intfc_options(), "
	           "unable to find time step %d\n",restart_time_step(init));
	    clean_up(ERROR);
	}

	read_next_dt(init,restart_io_type(init));

	restart_intfc(init) =
	    read_print_intfc_and_grids(init,restart_io_type(init),
				       &Comp_grid(init),&Restart_top_grid(init),
				       &Restart_comp_grid(init));

	got_intfc_from_file(init) = YES;

	debug_print("init","Left init_restart_interior_states()\n");
}       /*end to d_prompt_for_initial_intfc_options()*/


EXPORT	void	read_next_dt(
	INIT_DATA     *init,
	const IO_TYPE *io_type)
{
	FILE	      *file = io_type->file;
	int	      c;
	static OUTPUT *oput = NULL;

	oput = save_read_file_variables(file,oput);
	if (fgetstring(file,"next_dt:  ") == FUNCTION_FAILED)
	{
	    screen("ERROR in read_next_dt(), next_dt not found\n");
	    clean_up(ERROR);
	}	
	if ((c = getc(file)) != '\f')		/* NOBINARY */
	{
	    (void) ungetc(c,file);
	    (void) fscan_float(file,&restart_next_dt(init));
	}
	else
	{
	    (void) getc(file);
	    (void) read_binary_real_array(&restart_time(init),1,io_type);
	    (void) read_binary_real_array(&restart_next_dt(init),1,io_type);
	}

	if (debugging("prt_init"))
	{
	    screen("next_dt %g  Over-ride next_dt? (y,n(dflt)): ",
	    	   restart_next_dt(init));
	    (void) Gets(s);
	    if ((s[0] == 'y') || (s[0] == 'Y'))
	    {
		screen("Enter new dt: ");
		(void) Scanf("%f\n",&restart_next_dt(init));
	    }
	}
	reset_read_file_variables(oput);
}		/*end read_next_dt*/


EXPORT	boolean	position_file_at_read_time(
	FILE	*file,
	int	step,
	double	*time,
	double	*dt_last)
{
	const char *line;
	double	   jtime, jdt_last;
	int	   current_step;

	if (time == NULL)
	    time = &jtime;
	if (dt_last == NULL)
	    dt_last = &jdt_last;
 	while ((line = next_output_line_containing_string(file,"TIME DATA: ")))
	{
	    /*  Distinguish between blocked and unblocked output files */
 	    if (line[0] == '#')
		(void) sscanf(line,"%*s%*s%*s%*s%*s%lf%*s%*s%d%*s%*s%lf",
			      time,&current_step,dt_last);
	    else
		(void) sscanf(line,"%*s%*s%*s%*s%lf%*s%*s%d%*s%*s%lf",
			      time,&current_step,dt_last);
	    if (current_step == step)
		return FUNCTION_SUCCEEDED;
	}
	return FUNCTION_FAILED;
}		/*end position_file_at_read_time*/


EXPORT	INTERFACE *read_print_intfc_and_grids(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	RECT_GRID     *rgr,
	RECT_GRID     *rsTgr,
	RECT_GRID     *rsCgr)
{
	INTERFACE *intfc;
	FILE	  *file = io_type->file;   
	int	  dim = Comp_grid(init).dim;
	int	  c;
	int	  print_version;
	int	  grSet;

	if (next_output_line_containing_string(file,"FRONT DATA") == NULL)
	{
	    screen("ERROR in read_print_intfc_and_grids(), "
	           "Improper input file\n");
	    clean_up(ERROR);
	}

	/* Determine whether rectangular grids are printed */
	print_version = 0;
	while ((c = getc(file)) == ' ')
	    ++print_version;
	(void) ungetc(c,file);

	switch (print_version)
	{
	case 0:
	    if (debugging("init"))
	    	(void) printf("Rectangular grids not printed\n");
	    if (rgr == NULL)
	    {
	    	screen("ERROR in read_print_intfc_and_grids(), "
	    	       "rect grids not printed and rgr == NULL\n");
	    	clean_up(ERROR);
	    }
	    copy_rect_grid(rsCgr,rgr);
	    copy_rect_grid(rsTgr,rgr);
	    break;

	case 1:
	    if (debugging("init"))
	    	(void) printf("Interface grids not printed\n");
	    if (fgetstring(file,"Front Rectangular Grid:") == FUNCTION_FAILED)
	    {
		screen("ERROR in read_print_intfc_and_grids(), "
		       "comp grid not found\n");
		clean_up(ERROR);
	    }
	    read_rectangular_grid(io_type,rsCgr,YES,remap_info());
	
	    if (fgetstring(file,"Interface Topological Grid:")==FUNCTION_FAILED)
	    {
		screen("ERROR in read_print_intfc_and_grids(), "
		       "top grid not found\n");
		clean_up(ERROR);
	    }
	    read_rectangular_grid(io_type,rsTgr,NO,remap_info());
	    break;

	case 2:
	default:
	    if (debugging("init"))
		(void) printf("Interface grids printed\n");
	    break;
	}

	set_interface_hooks(dim,init);
	set_size_of_intfc_state(StateSize(init));
	if ((intfc = read_print_interface(init,io_type,NO,&grSet)) == NULL)
	{
	    screen("ERROR in read_print_intfc_and_grids(), "
	           "read_print_interface failed\n");
	    clean_up(ERROR);
	}
	if (debugging("init") || debugging("consistency"))
	{
	    (void) printf("Checking consistency of restart interface\n");
	    if (!consistent_interface(intfc))
	    {
	        screen("ERROR in read_print_intfc_and_grids(), "
		       "inconsistent input interface\n");
	        clean_up(ERROR);
	    }
	    else
	    {
		(void) printf("Restart interface is consistent\n");
	    }
	}

	if (grSet == NO)
	{
	    set_topological_grid(intfc,rsTgr);
	    set_computational_grid(intfc,rsCgr);
	}
	else
	{
	    copy_rect_grid(rsTgr,&topological_grid(intfc));
	    copy_rect_grid(rsCgr,computational_grid(intfc));
	}
	return intfc;
}		/*end read_print_intfc_and_grids*/


/*
*				start_up():
*
*	Calls Routines to handle system error messages, and prints
*	Current messages at top of output or on the screen.
*
*	Command line arguments:
*
*	  I/O redirection
*	    -i file : read input from file
*	    -o name : use name as the base for output file names
*	    -e file : write stderr messages to file
*
*	  Parallel Processing control
*	    -p|-par nx | nx ny | nx ny nz :
*		Specify parallel domain decomposition. Depending on the *
*	        spatial dimension of the file x, y, or z defines the number
*		of subdomains in the corresponding coordinate direction.
*	    -b|-buf nx | nx ny | nx ny nz :
*		Specify number of buffer cell zones for the subdomains
*		in the corresponding coordinate direction.
*	    -block|-noblock : use blocking/unblocking parallel messages.
*	    -pp_recv_w n : set a parallel message wait time of n seconds
*	    -pp_recv_t t : set a parallel message timeout time for t seconds
*	    -msg_buf size : set a parallel message buffer size of size size.
*	    -BlockSize size : sets the hyperbolic solutiion message block
*			      size to size size.
*
*	Restart file layout
*	    -reverse_endian : For binary files that did not include information
*			      regarding the endian of their printout,  this
*			      directive asserts that such files have the
*			      opposite endian from the running process.
*	    -read_big_endian|-read_little_endian :
*			      Specifies that binary files that did not
*			      include information regarding their endian
*			      are either in big/little endian format.
*	    -read_float_size N : Specifies that binary files that do not
*			         include the size of the real words printed
*			         in the file have a real (ie double) size
*				 N.  N must be one of sizeof(double) or
*				 sizeof(double).
*	NOTE: FronTier versions with date later than September 24, 2001
*	all print information regarding the endian and floating point size
*	in all output files.  For restarts from such files the above three
*	arguments are ignored and over ridden by the information included
*	in the file.
*
*	Interface storage
*	    -Chunksize size : Sets the interface chunk size.
*
*	Other arguments
*	    -compress : Use gzip compression on outfile files.
*	    -usage : Print a usage message and terminate the run.
*/

#define SHIFT (--argc,++argv)

/*ARGSUSED*/
LOCAL void start_up(
	int		*pargc,
	char		***pargv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	IMPORT boolean       suppress_prompts;
	CHART	          *root = ip->root;
	FILE	          *stream[3];
	Grid	          *grid = root->grid;
	Printplot         *prt = ip->prt;
	char	          file_name[3][1024];
	char	          **argv;
	int	          i;
	int	          argc;
	boolean	          use_parallel = NO;
	const char        *execname;
	static const char *rw[3] = { "r", "w", "w"};
	static const char *stream_name[3] = { "stdin", "stdout", "stderr"};
	static const int DEFAULT_BUF_WIDTH = 3;
        int reduce, numnodes;

	suppress_prompts = NO;
	(void) pp_init(pargc,pargv);
	argc = *pargc;
	argv = *pargv;

	execname = basename(argv[0]);
	for (i = 0; i < 3; ++i)
	{
	    subdomains(init)[i] = 1;
	    buffer_zones(init)[i] = DEFAULT_BUF_WIDTH + MAXD/3;
	}
	pp_grid_set(init) = YES;
	output_filename(init)[0] = '\0';
	for (i = 0; i < 3; ++i)
	{
	    stream[i] = NULL;
	    file_name[i][0] = '\0';
	}
	SHIFT;
	while (argc && argv[0][0]=='-')
	{
	    if ( (strcmp(argv[0],"-p") == 0) ||
	    	 (strncmp(argv[0],"-par",4) == 0) )
	    {
	    	use_parallel = YES;
	    	for (i = 0; i < MAXD; ++i)
	    	{
	            if (argc < 2 || argv[1][0] == '-') break;
	            SHIFT;
	            subdomains(init)[i] = atoi(argv[0]);
	    	}
	    }
	    else if ( (strcmp(argv[0],"-b") == 0) ||
	              (strncmp(argv[0],"-buf",4) == 0) )
	    {
	    	for (i = 0; i < MAXD; ++i)
	    	{
	            if (argc < 2 || argv[1][0] == '-') break;
	            SHIFT;
	            buffer_zones(init)[i] = atoi(argv[0]);
	    	}
	    }
	    else if (strncmp(argv[0]+1,"c",1) == 0)
	    {
	    	compress_output(init) = YES;
	    }
	    else if (strncmp(argv[0]+1,"C",1) == 0)
	    {
	    	if (argc < 2 || argv[1][0] == '-') break;
	    	SHIFT;
	    	SetChunkSize((size_t)atoi(argv[0]));
	    }
	    else if (strncmp(argv[0]+1,"B",1) == 0)
	    {
	    	if (argc < 2 || argv[1][0] == '-') break;
	    	SHIFT;
	    	SetHypPPBlockSize((size_t)atoi(argv[0]));
	    }
	    else if (strncmp(argv[0]+1,"s",1) == 0)
	    {
	    	use_parallel = NO;
	    }
	    else if (strncmp(argv[0]+1,"i",1) == 0)
	    {
	    	SHIFT;
	    	stream[0] = stdin;
		stripcomm(file_name[0],argv[0]);
	        grid->temporary_input_file =
		    (strcmp(file_name[0],argv[0]) != 0) ?
	                strdup(file_name[0]) : NULL;
	    	suppress_prompts = YES;
	    }
	    else if (strncmp(argv[0]+1,"o",1) == 0)
	    {
		char dname[1024];
	    	SHIFT;
	    	stream[1] = stdout;
	    	(void) strcpy(file_name[1],argv[0]);
	    	(void) strcpy(dname,argv[0]);
	    	(void) strcpy(output_filename(init),argv[0]);
		(void) get_dirname(dname);
		if (create_directory(dname,NO) == FUNCTION_FAILED)
		{
		    screen("ERROR in start_up(), can't create directory %s\n",
			   dname);
		    clean_up(ERROR);
		}
	    }
	    else if (strncmp(argv[0]+1,"e",1) == 0)
	    {
	    	SHIFT;
	    	stream[2] = stderr;
	    	strcpy(file_name[2],argv[0]);
	    }
	    else if (strcmp(argv[0],"-pp_recv_w") == 0)
	    {
	    	SHIFT;
	    	set_pp_recv_wait_interval((unsigned)atoi(argv[0]));
	    }
	    else if (strcmp(argv[0],"-pp_recv_n") == 0)
	    {
	    	SHIFT;
	    	set_pp_recv_num_retries((unsigned)atoi(argv[0]));
	    }
	    else if (strcmp(argv[0],"-pp_recv_t") == 0)
	    {
	    	SHIFT;
	    	set_pp_recv_timeout((unsigned)atoi(argv[0]));
	    }
	    else if (strcmp(argv[0],"-msg_buf") == 0)
	    {
	    	size_t	buf_size;
	    	SHIFT;
	    	buf_size = (size_t) atoi(argv[0]);
	    	set_MSG_BUF_SIZE(buf_size);
	    }
#if defined(__MPI__)
	    else if (strcasecmp(argv[0],"-block") == 0)
		use_blocking_pp_comp(init) = YES;
	    else if (strcasecmp(argv[0],"-noblock") == 0)
		use_blocking_pp_comp(init) = NO;
#endif /* defined(__MPI__) */
	    else if (strcasecmp(argv[0],"-reverse_endian") == 0)
	    {
		set_reverse_endian(YES);
	    }
	    else if (strcasecmp(argv[0],"-read_big_endian") == 0)
	    {
		set_read_endian(FT_BIG_ENDIAN);
	    }
	    else if (strcasecmp(argv[0],"-read_little_endian") == 0)
	    {
		set_read_endian(FT_LITTLE_ENDIAN);
	    }
	    else if (strcasecmp(argv[0],"-read_float_size") == 0)
	    {
	        size_t size;
	    	SHIFT;
		if (sscanf(argv[0],"%lu",&size) != 1)
		{
		    screen("ERROR in start_up(), invalid read float size\n");
		    usage(execname);
		    clean_up(ERROR);
		}
		switch (size)
		{
		case sizeof(double):
		case sizeof(TRUEfloat):
		    break;
		default:
		    screen("ERROR in start_up(), "
		           "invalid read float size %lu\n",size);
		    usage(execname);
		    clean_up(ERROR);
		}
		set_read_float_size(size);
	    }
	    else if (strcasecmp(argv[0],"-usage") == 0)
	    {
		usage(execname);
		exit(0);
	    }
	    else
	    {
		screen("Unknown argument %s ignored\n",argv[0]);
		usage(execname);
	    }
	    SHIFT;
	}

	/* This code prevents a race condition bug. If the output directory has
	to be created then proc 0 is slowed down. Hence the other procs reach
	the stdin/stdout/stderr redirection code first. This means they attempt
	to redirect to a file in a non-existant directory, which causes freopen
	to fail and FronTier to exit with an error.*/
        reduce = 1;   
        numnodes = pp_numnodes();
        pp_global_isum(&reduce,1);
        if(reduce != numnodes)
        {
            screen("ERROR in start_up()");
            clean_up(ERROR);
        }

	if (use_parallel == YES)
	{
	    int nn, nd, myid;
	    int nn_world = pp_numnodes();
	    for (nn = 1, i = 0; i < 3; ++i)
	        nn *= subdomains(init)[i];
	    if (nn > nn_world)
	    {
		/*TMP*/
		screen("nn = %d  nn_world = %d\n",nn,nn_world);
	        screen("ERROR in start_up(), too many subdomains\n");
		clean_up(ERROR);
	    }
	    pp_comm_split(nn);
	    myid = pp_mynode();
	    for (nd = 0; nn != 0; nn /=10, ++nd);
	    if (stream[0] != NULL || stream[1] != NULL || stream[2] != NULL)
	    {
	    	for (i = 1; i < 3; ++i)
	    	{
	            if (stream[i] == NULL) continue;
	            if (strncmp(file_name[i],"/dev/null",9) == 0)
	            	continue;
	            (void) sprintf(s,"%s.%s",file_name[i],right_flush(myid,nd));
	            (void) strcpy(file_name[i],s);
	    	}
	    	for (i = 0; i < 3; ++i)
	    	{
	            if (stream[i] == NULL) continue;
	            if (freopen(file_name[i],rw[i],stream[i]) == NULL)
	            {
	                screen("ERROR in start_up(), can't reopen %s to %s\n",
	    	    	       stream_name[i],file_name[i]);
	                clean_up(ERROR);
	            }
	    	}
	    }
	}
	else
	{
	    for (i = 0; i < 3; ++i)
	    {
	    	if (stream[i] == NULL) continue;
	    	if (freopen(file_name[i],rw[i],stream[i]) == NULL)
	    	{
	            screen("ERROR in start_up(), can't reopen "
			   "%s to %s\n",stream_name[i],file_name[i]);
		    clean_up(ERROR);
		}
	    }
	}
#if defined(ultrix)
	(void) fprintf(stderr,
	               "\tUltrix - using default buffering on stdin\n\n");
#else /* defined(ultrix) */
	setbuf(stdin,NULL);
#endif /* defined(ultrix) */

	prt->init = init;
	prt->outfile = (output_filename(init)[0] != '\0') ?
	               strdup(output_filename(init)) : NULL;
	prt->compress = compress_output(init);

	set_error_immediate(stdout);

#if defined(VERSION0)
	set_print_version(0);
#endif /* defined(VERSION0) */
	record_print_version(stdout);
	print_title(stdout,title(init));

	init_prompting_and_debugging(init);

	init_clean_up(d_clean_up,d_clean_up_printout);
}		/*end start_up*/

LOCAL	void	usage(const char *execname)
{
	size_t len;
	int    i;
	static const char *indent = "        ";
	static const char *opts[] = {
	                             "Usage - ",
	                             NULL,
				     "[-i infile]",
				     "[-o outfile]",
				     "[-e errfile]",
#if defined(__MPI__)
				     "[-p [x [y [z]]]|-s]",
				     "[-b [nx [ny [nz]]]]",
	                             "[-block|-noblock]",
	                             "[-pp_recv_w n]",
				     "[-pp_recv_t t]",
	                             "[-msg_buf size]",
				     "[-BlockSize size]",
#endif /* defined(__MPI__) */
				     "[-reverse_endian]",
				     "[-read_big_endian}|-read_little_endian]",
				     "[-read_float_size N]",
	                             "[-Chunksize size]",
	                             "[-compress]",
	                             "[-usage]",
				    NULL};
	opts[1] = execname;
	screen("%s",opts[0]);
	for (len = strlen(opts[0]), i = 1; opts[i] != NULL; ++i)
	{
	    len += strlen(opts[i]);
	    if (len < 70)
	        screen(" %s",opts[i]);
	    else
	    {
	        len = strlen(opts[i]) + strlen(indent);
		screen("\n%s%s",indent,opts[i]);
	    }
	}
	screen("\n");
}		/*end usage*/


LOCAL int to_parallel_restart(
	RECT_GRID	*rst_grid)
{
	double		*L = rst_grid->L, *U = rst_grid->U, *h = rst_grid->h;
	double		*VL = rst_grid->VL, *VU = rst_grid->VU;
	double		tol;
	int		i;

	for (i = 0; i < rst_grid->dim; ++i)
	{
	    tol = 0.0001*h[i];/*TOLERANCE*/
	    if (fabs(L[i] - VL[i]) > tol)
		return YES;
	    if (fabs(U[i] - VU[i]) > tol)
		return YES;
	}
	return NO;
}		/*end to_parallel_restart*/

/*
*			init_interior_states():
*
*	Initializes the states in a wave structure by calling
*
*		(*initializer)(coords,comp,state,intfc,init)
*
*	at the centers of the grid blocks of wave->rect_grid.
*/

LOCAL	void init_interior_states(
	Wave		*wave,
	Front		*front,
	INIT_DATA	*init,
	void		(*initializer)(double*,COMPONENT,Locstate,INTERFACE*,
				       INIT_DATA*))
{
	COMPONENT	comp;
	Locstate	state;
	double		*coords;
	int		icoords[MAXD];
	int		dim = wave->rect_grid->dim;
	int		status;

	debug_print("init","Entered init_interior_states()\n");

	if (wave->sizest == 0 || initializer == NULL)
	{
	    debug_print("init","Left init_interior_states()\n");
	    return;
	}

	wave->old_wave = NULL;
	status = init_hyp_solution_function(wave,front);

	if (status != GOOD_STEP)
	{
	    screen("ERROR in init_interior_states(), "
	           "init_hyp_solution_function() failed\n");
	    print_interface(front->interf);
	    clean_up(ERROR);
	}

	switch (dim)
	{
	case 1:
	{
	    int		ix;
	    int		xmax;

	    xmax = wave->rect_grid->gmax[0];
	    for (ix = 0; ix < xmax; ++ix)
	    {
	    	icoords[0] = ix;
	    	coords = Rect_coords(icoords,wave);
	    	comp = Rect_comp(icoords,wave);
	    	state = Rect_state(icoords,wave);
	    	(*initializer)(coords,comp,state,front->interf,init);
	    }
	    break;
	}
	case 2:
	{
	    int		ix, iy;
	    int		xmax, ymax;

	    xmax = wave->rect_grid->gmax[0];
	    ymax = wave->rect_grid->gmax[1];
	    for (iy = 0; iy < ymax; ++iy)
	    {
	    	icoords[1] = iy;
	    	for (ix = 0; ix < xmax; ++ix)
	    	{
	    	    icoords[0] = ix;
	    	    coords = Rect_coords(icoords,wave);
	    	    comp = Rect_comp(icoords,wave);
	    	    state = Rect_state(icoords,wave);
	    	    (*initializer)(coords,comp,state,front->interf,init);
	    	}
	    }
	    break;
	}
	case 3:
	{
	    int		ix, iy, iz, i;
	    int		imax[3], imin[3];
	    RECT_GRID	*rgr = wave->rect_grid;

	    for(i=0; i<3; i++)
	    {
		imin[i] = 0;
		imax[i] = rgr->gmax[i];
		
		if(rect_boundary_type(front->interf,i,0) == OPEN_BOUNDARY)
		    imin[i] = -rgr->lbuf[i];
		if(rect_boundary_type(front->interf,i,1) == OPEN_BOUNDARY)
		    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	    }
	    
	    for (iz = imin[2]; iz < imax[2]; ++iz)
	    {
	    	icoords[2] = iz;
	    	for (iy = imin[1]; iy < imax[1]; ++iy)
	    	{
	    	    icoords[1] = iy;
	    	    for (ix = imin[0]; ix < imax[0]; ++ix)
	    	    {
	    	    	icoords[0] = ix;
	    	    	coords = Rect_coords(icoords,wave);
	    	    	comp = Rect_comp(icoords,wave);
	    	    	state = Rect_state(icoords,wave);
	    	    	(*initializer)(coords,comp,state,front->interf,init);
	    	    }
	    	}
	    }
	    break;
	}
	}

	debug_print("init","Left init_interior_states()\n");
}		/*end init_interior_states*/


/*
*			init_pp_grid():
*
* 	Initializes PP grid for parallel computing.
*
*	The pp_grid is set in this function through I/O. An alternative 
*	is to insert the struct pp_grid to the present struct 
*	Grid defined in ddecs.h. It appears not to be the best chioce
*	as no-trivial modifications are needed to the calling to function 
*	init_hyp_and_top_grid().
*/


LOCAL	void init_pp_grid(
	INIT_DATA	*init)
{
	int		i, dim = i_intfc(init)->dim;
	int		num_proc = pp_numnodes();
	int     	nn; /* total number of subdomains */

	if (pp_grid_set(init) == NO)
	{
	    screen("\n\t\tParallel Domain Decomposition Control\n");
	    if (num_proc > 1)
	    {
	        static const char *dname[3] = {"x", "y", "z"};
	        for (i = 0; i < dim; ++i)
	        {
	            screen("Enter the number of subdomain blocks "
	                   "in the %s direction (dflt = %d): ",
	    		    dname[i],subdomains(init)[i]);
	            (void) Gets(s);
	            if (s[0] != '\0')
	                (void) sscanf(s,"%d",subdomains(init)+i);
	        }
	    }
	    for (i = 0; i < dim; ++i)
	    {
	        int buf;
	        screen("Enter the subdomain buffer size in the "
	               "%d-%s direction, (dflt = min = %d): ",
		       i,ordinal_suffix(i),buffer_zones(init)[i]);
	        (void) Gets(s);
	        if (s[0] != '\0')
	            (void) sscanf(s,"%d",&buf);
	        if (buf > buffer_zones(init)[i])
	            buffer_zones(init)[i] = buf;
	    }
	}

	for (nn = 1, i = 0; i < dim; ++i)
	    nn *= subdomains(init)[i];

	/* check if mapping is reasonable */

	if (nn > num_proc)
	{
	    screen("ERROR in init_pp_grid(), "
		   "more subdomains than processors\n");
	    (void) printf("(number of subdomains = %d) > "
	                  "(number of processors = %d), "
	                  "unable to proceed!\n",nn,num_proc);
	    clean_up(ERROR);
	}
	else if (nn < num_proc) 
	{
		screen("ERROR in init_pp_grid(), "
		       "more processors than subdomains\n");
		(void) printf("\n(number of subdomains = %d) < "
		              "(number of processors = %d), "
		              "unable to proceed!\n",nn,num_proc);
		clean_up(ERROR);
	}
}		/*end init_pp_grid*/


/*
*			d_set_pp_grid():
*
*	Sets up the pp_grid structure for parallel runs.
*/

LOCAL	void d_set_pp_grid(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	CHART		*root = ip->root;
	Front		*front = root->front;
	Grid		*grid = root->grid;
	Wave		*wave = root->wave;
	RECT_GRID	*comp_glbgr = front->rect_grid;
	RECT_GRID	*top_glbgr = &ip->top_grid;
	double		*h = comp_glbgr->h;
	int		i, dim = comp_glbgr->dim;
	PP_GRID		*pp_grid;

	debug_print("init_pp_grid","Entered d_set_pp_grid():\n");

	pp_grid = set_pp_grid(init,comp_glbgr);
	front->pp_grid = wave->pp_grid = grid->pp_grid = pp_grid;

        if (dim != 1)
        {
            screen("\nEnter yes for adaptive partition (default is no): ");
            (void) Gets(s);
            if (s[0] == 'y' || s[0] == 'Y')
                front->adaptive_partition = YES;
            else
                front->adaptive_partition = NO;
        }
        else
            front->adaptive_partition = NO;
	if (dim == 3) /* TODO Unify 2 and 3 D */
	{
	    double *VU = comp_glbgr->VU, *VL = comp_glbgr->VL;
	    double *GU = comp_glbgr->GU, *GL = comp_glbgr->GL;
	    int tgmax[MAXD];

	    copy_rect_grid(comp_glbgr,&pp_grid->Zoom_grid);
	    for (i = 0; i < dim; ++i)
	    {
	    	double h = top_glbgr->h[i];
	    	tgmax[i] = irint((VU[i] - VL[i])/h);
	    }
	    set_rect_grid(VL,VU,GL,GU,NOBUF,NOBUF,tgmax,dim,&comp_glbgr->Remap,
			  top_glbgr);
	    (void) adjust_top_grid_for_square(top_glbgr,comp_glbgr);
	}
	if (dim == 3) /* TODO Unify 2 and 3 D */
	{
	    screen("Enter yes to re-partition the output at the end of run: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	    	static PP_GRID new_Pp_grid;
	    	grid->repart_at_end_of_run = YES;
	    	grid->new_pp_grid = &new_Pp_grid;
	    	screen("New partition must be multiple of current partition\n");
	    	screen("Current partition is (");
	    	for (i = 0; i < dim; ++i)
	    	    screen("%d ",pp_grid->gmax[i]);
	    	screen(")\n");
	    	screen("Enter the new partition at the end of run: ");
	    	for (i = 0; i < dim; ++i)
	    	    Scanf("%d",&new_Pp_grid.gmax[i]);
	    	screen("\n");
	    	copy_rect_grid(&new_Pp_grid.Global_grid,&pp_grid->Global_grid);

	    	/* Setting up final output partition */
	    	grid->new_pp_grid->nn = 1;
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    int	Gmax, Pmax, k;
	    	    int	basic_slices, extra_slices;

	    	    grid->new_pp_grid->buf[i] = buffer_zones(init)[i];
	    	    Pmax = grid->new_pp_grid->gmax[i];
	    	    grid->new_pp_grid->nn *= Pmax;

	    	    uni_array(&grid->new_pp_grid->dom[i],Pmax + 1,FLOAT);

	    	    grid->new_pp_grid->dom[i][0]    = pp_grid->Global_grid.L[i];
	    	    grid->new_pp_grid->dom[i][Pmax] = pp_grid->Global_grid.U[i];
	    	    Gmax = pp_grid->Global_grid.gmax[i];

	    	    basic_slices = Gmax / Pmax;
	    	    extra_slices = Gmax % Pmax;

	    	    for (k = 1; k < Pmax; ++k)
	    	    {
	    	    	if (k < extra_slices)
	            	    grid->new_pp_grid->dom[i][k] = k*(basic_slices 
				+ 1)*h[i] + grid->new_pp_grid->dom[i][0];
	            	else
	            	    grid->new_pp_grid->dom[i][k] = (k*basic_slices + 
			        extra_slices)*h[i]
	        		+ grid->new_pp_grid->dom[i][0];
		    }
	    	}
	    	(void) Gets(s);
	    }
	}

	if (debugging("d_set_pp_grid"))
	{
	    (void) printf("pp_grid after d_set_pp_grid()\n");
	    (void) print_PP_GRID_structure(pp_grid);
	    if (grid->repart_at_end_of_run)
	    {
	    	(void) printf("Final repartition pp_grid:\n");
	    	print_PP_GRID_structure(grid->new_pp_grid);
	    }
	}

	debug_print("init_pp_grid","Left d_set_pp_grid():\n");
}		/*end d_set_pp_grid*/

LOCAL	void init_spatial_grids(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	RECT_GRID	*r_grid;
	RECT_GRID	*top_grid;
	int		dim;

	debug_print("init","Entered init_spatial_grids()\n");

		/* Initialize the Units */
	init_physical_units(init,ip);

		/* Initialize Rectangular Grid */

			/* Get spatial dimension*/
	screen("\nEnter the spatial dimension of the computation: ");
	(void) Scanf("%d\n",&dim);
	switch(dim)
	{
	case 1:
	    break;
	case 2:
	    break;
	case 3:
	    break;
	default:
	    screen("ERROR in init_spatial_grids():\n");
	    if (dim==1 || dim==2 || dim==3) 
	    {
	    	screen("Computational dimension %d is not available"
	    	       "in this executable.\n");
	    }	
	    else
	    	screen("Invalid computational dimension %d\n",dim);
	    clean_up(ERROR);
	}
	SetDefaultHypPPBlockSize(dim);
	set_interface_hooks(dim,init);
	set_size_of_intfc_state(0);
	i_intfc(init) = make_interface(dim);
	r_grid = &Comp_grid(init);
	top_grid = &Top_grid(init);
	r_grid->dim = dim;
	init_remap_and_rect_grid(r_grid,ip);
	init_topological_grid(top_grid,r_grid);

	debug_print("init","Left init_spatial_grids()\n");
}		/*end init_spatial_grids*/

LOCAL	void	init_time_step_limit(
	INIT_DATA	*init)
{
	initial_dt_lim(init) = HUGE_VAL;
	if (debugging("limit_time_step"))
	{
	    screen("Input upper bound on the time step (dflt = none): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscan_float(s,&initial_dt_lim(init));
	}
}		/*end init_time_step_limit*/

LOCAL   void    init_max_time_step_mods(
        INIT_DATA  *init)
{
	initial_max_num_time_step_mods(init) = 50; /* DEFAULT */
        screen("Enter maximum number of time step modifications allowed\n");
        screen("\tduring a propagation step (default = %d): ",
	       initial_max_num_time_step_mods(init));
        (void) Gets(s);
        if (s[0] != '\0')
            (void) sscanf(s,"%d",&initial_max_num_time_step_mods(init));
        screen("\n");      
}              /*end init_max_time_step_mods*/

EXPORT	void d_set_interface_hooks(
	int		dim,
	INIT_DATA       *init)
{
	h_set_interface_hooks(dim,init);
}		/*end d_set_interface_hooks*/
