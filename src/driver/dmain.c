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
*				dmain.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*		This is the main program for tracked front
*		hydrodynamics, combined version.
*		For compilation details see the makefile.
*/


#include <driver/ddecs.h>


LOCAL	Printplot *prt = NULL;
LOCAL	CHART *root = NULL;


	/* LOCAL Function Declarations */
LOCAL	boolean	d_last_time_step_modification(void);
LOCAL	void	end_of_run(Grid*,Front*,Wave*);
LOCAL	void	main_time_step_loop(INIT_PHYSICS*,CHART*,Printplot*);
LOCAL	void	print_memory_statistics(CHART*);
LOCAL	void	set_modified_time_step(Grid*,Front*,Wave*,Printplot*,double);


EXPORT	int dmain(
	int		argc,
	char**		argv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	perform_initialization(argc,argv,init,ip);

#if defined (USE_OVERTURE)
        if(ip->root->overparam != NULL)
        {
            d_overture_amr_main(init, ip);

            end_of_run(ip->root->grid,ip->root->front,ip->root->wave);
            return 0;
        }
        screen("ERROR In dmain()\n");
        screen("USE_OVERTURE is defined, but NO AMR params are provided\n");
        screen("Require to re-compile the code with the option: use_overture=NO\n");
        clean_up(ERROR);
#endif /* defined(USE_OVERTURE) */
	main_time_step_loop(ip,ip->root,ip->prt);

	end_of_run(ip->root->grid,ip->root->front,ip->root->wave);

	return 0;
}		/*end dmain*/


/*ARGSUSED*/
LOCAL	boolean d_last_time_step_modification(void)
{
	return (root->grid->num_mts < max_num_time_step_mods(root->grid)) ?
		NO : YES;

}		/*end d_last_time_step_modification*/

EXPORT	void d_clean_up_printout(int error)
{
	if (root->grid->initialization_complete) 
	{
	    (void) output();
	    (void) printf("\t\tPRINTOUT DURING CLEANUP:\n");
	    if (prt->printout != NULL)
	    	(*prt->printout)(root,prt,YES,error);
	}
}		/*end d_clean_up_printout*/

EXPORT	void d_clean_up(void)
{
	if ((root->grid->temporary_input_file != NULL) &&
						(is_io_node(pp_mynode())))
	    (void) unlink(root->grid->temporary_input_file);
}		/*end d_clean_up*/


EXPORT	void perform_initialization(
	int		argc,
	char		**argv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Grid		*grid;
	IMPORT boolean	suppress_prompts;

	prt = ip->prt;
	prt->title = title(init);
	root = ip->root;
	grid = root->grid;

	grid->num_mts = 0;
	set_current_chart(root);
	root->front->_last_time_step_modification = 
	    d_last_time_step_modification;

	root->grid->initialization_complete = NO;
	print_storage("before INIT","INIT_storage");
	if (ip->init_run)
	    (*ip->init_run)(&argc,&argv,init,ip);
	if (i_intfc(init) != NULL)
	    (void) delete_interface(i_intfc(init));
	i_intfc(init) = NULL;
	root->grid->initialization_complete = YES;
	suppress_prompts = NO;
	if (debugging("init"))
	{
	    (void) printf("\n\tStructures printout after init_run():\n\n");
	    print_Grid_structure(root->grid);
	    print_Front_structure(root->front);
	    print_Wave_structure(root->wave);
	    print_Printplot_structure(prt);
	}
	print_storage("after INIT","INIT_storage");

	if (prt->print_initial_data)
	    (*prt->print_initial_data)(stdout,root,prt);

	if ((prt->printout != NULL) && debugging("prt_rst"))
	{
	    STOP	Stop;
	    boolean	initial_step;
	    double	dt_sav = root->grid->dt;
	    char	*sav_outfile = NULL;

	    Stop = *stop(grid);
	    root->grid->dt = restart_dt_last(init);
	    if (root->grid->step > 0)
	    {
	    	initial_step = NO;
	    	(void) printf("Restart Timestep Printout\n\n");
	    	root->grid->step--;
	    }
	    else
	    {
	    	initial_step = YES;
	    	(void) printf("Printout before initial timestep\n\n");
	    	initialize_max_front_speed(root->front);
	    	initialize_max_wave_speed(root->wave);
		sav_outfile = prt->outfile;
		prt->outfile = NULL;
	    }
	    (*prt->printout)(root,prt,NO,0);
	    *stop(grid) = Stop;
	    end_of_run(root->grid,root->front,root->wave);
	    if (initial_step == NO)
	    	++root->grid->step;
	    else
		prt->outfile = sav_outfile;
	    root->grid->dt = dt_sav;
	    remove_from_debug("prt_rst");
	}
	if (debugging("wave_init"))
	    (*root->wave->show_wave_states)(root->wave);
}		/*end perform_initialization*/

LOCAL	void main_time_step_loop(
	INIT_PHYSICS    *ip,
	CHART		*root,
	Printplot	*prt)
{
	Front		*fr = root->front;
	Grid		*grid = root->grid;
	Wave		*wave = root->wave;
	double		dt_frac;
	int		status;
	static int	initial_step;
	FILE		*file;
#if defined(TIMING)
	double		start_cpu_time, end_cpu_time;
	double		start_wall_time, end_wall_time;
#endif /* defined(TIMING) */
	int *iperm,i,dim;

				/* Time Loop */
	DEBUG_ENTER(main_time_step_loop)

	debug_print("time_step","Entered main_time_step_loop(), step = %d\n",
			  grid->step);
	initial_step = grid->step;

	start_clock("ALL_TIMESTEPS");

	add_time_clear(0);

	for ( ; ; )
	{
	    if (debugging("memory"))
	    	print_memory_statistics(root);

	    start_clock("TIMESTEP");

#if defined(TIMING)
	    start_cpu_time = cpu_seconds();
	    start_wall_time = real_time();
#endif /* defined(TIMING) */

            if(ip->set_gravity_charts)
                (*ip->set_gravity_charts)(&root,1);
			
	    dim = fr->rect_grid->dim;
            iperm = set_iperm(fr->step+1,dim);
            for (i = 0; i < dim; i++)
            {
                if (scatter_states(wave,fr,iperm,i) != FUNCTION_SUCCEEDED)
                {       
                    screen("ERROR in init_states(), scatter_states() failed\n");
                    clean_up(ERROR);
                }
            } 

	    status = time_step(root,grid->dt,&dt_frac,grid->step,grid->time);
	    
#if defined(TIMING)
	    end_cpu_time = cpu_seconds();
	    end_wall_time = real_time();
	    (void) printf("# Step %d took %g CPU seconds, %g wall seconds.\n",
	    	grid->step,end_cpu_time - start_cpu_time,
	    	end_wall_time - start_wall_time);
#endif /* defined(TIMING) */
	    delete_untracked_hyper_surfaces(fr,wave);

	    switch (status)
	    {
	    case GOOD_STEP:
	    	break;

	    case MODIFY_TIME_STEP:
	    	set_modified_time_step(grid,fr,wave,prt,dt_frac);
	        stop_clock("TIMESTEP");
	    	continue;

	    case REPEAT_TIME_STEP:
	        stop_clock("TIMESTEP");
	    	continue;

	    case ERROR_IN_STEP:
	    default:
	        stop_clock("TIMESTEP");
	    	screen("ERROR: in main(), time_step() failed\n");
	    	clean_up(ERROR);
	    }

	    /* Reset modified time step counter after good time step */
	    fr->num_mts = grid->num_mts = 0;

	    grid->time += grid->dt;
	    fr->time = wave->time = grid->time;

	    if (root->level != NULL)
	    {
	    	if (debugging("time_step"))
	    	{
	            (void) printf("\nBefore amr: interf %llu\n",
	            	          interface_number(root->front->interf));
	            print_interface(root->front->interf);
	    	}
	        status = amr(root,root->level,prt,initial_step,
			     grid->time,grid->dt);
	        switch (status)
	        {
	        case GOOD_STEP:
	    	    break;

	        case MODIFY_TIME_STEP:
	    	    set_modified_time_step(grid,fr,wave,prt,dt_frac);
	    	    continue;

	        case REPEAT_TIME_STEP:
	    	    continue;

	        case ERROR_IN_STEP:
	        default:
	    	    screen("ERROR: in main(), time_step() failed\n");
	    	    clean_up(ERROR);
	        }
	    	interpolate_from_finer_grids_to_level(root->level);
	    }
	    else
	    {
	    	if (grid->step != initial_step)
	    	{
	    	    if (root->ellip)
	    	    {
	                print_storage("before ELLIP","EL_storage");
                        start_clock("ELLIP");	  /* Elliptic Step */
	                if (debugging("time_step"))
	                {
	                   (void) printf("\nBefore ellip: interf %llu\n",
	                                 interface_number(fr->interf));
	                   print_interface(fr->interf);
	                }
	                (*root->ellip)(NULL,root,NULL);
	                if (debugging("time_step"))
	                {
	                    (void) printf("\nAfter ellip: interf %llu\n",
	                                  interface_number(fr->interf));
	                    print_interface(fr->interf);
	                }
	                stop_clock("ELLIP");
	                print_storage("after ELLIP","EL_storage");
	    	    }
	    	}
	    }

	    if (prt->printout != NULL)
	    {
	    	if (stop_run(grid,fr,wave))
	    	{
	            start_clock("PRINTOUT");
	            (*prt->printout)(root,prt,YES,0);
	            stop_clock("PRINTOUT");
	            break;
	    	}
	    	else
	    	{
	            start_clock("PRINTOUT");
	            (*prt->printout)(root,prt,NO,0);
	            stop_clock("PRINTOUT");
	    	}
	    }
	    else if (stop_run(grid,fr,wave))
	    	break;

	            /* Update Time Step */
	    grid->dt_last = grid->dt;
	    grid->dt = find_time_step(grid,wave,fr,prt,YES);
	    ++grid->step;

	    stop_clock("TIMESTEP");

	    print_storage("after TIMESTEP","TIME_storage");
	    if (debugging("num_intfcs"))
	    {
	        struct Table *T;
	        int n;

	        for (n = 0, T = interface_table_list(); T != NULL; T = T->next)
	        {
	    	    ++n;
	    	    (void) printf("Interface of table %llu = %llu\n",
	            	  table_number(T),interface_number(T->interface));
	        }
	        (void) printf("%d interfaces present at end of time step\n",n);
	    }

	}
	stop_clock("ALL_TIMESTEPS");
	debug_print("time_step","Left main_time_step_loop(), step = %d\n",grid->step);

	DEBUG_LEAVE(main_time_step_loop)
}		/*end main_time_step_loop*/

/*ARGSUSED*/
LOCAL	void end_of_run(
	Grid		*grid,
	Front		*front,
	Wave		*wave)
{
	if (grid->step >= stop_step(grid))
	    screen("\nSTOP: time step limit grid stop_step = %d reached\n",
		   stop_step(grid));
	if (grid->time >= stop_time(grid))
	    screen("\nSTOP: time limit grid stop_time = %g reached\n",
		   stop_time(grid));

	if (d_stop_run(grid,front,wave))
	    clean_up(0);
}		/*end end_of_run*/

/*ARGSUSED*/
EXPORT	int d_stop_run(
	Grid		*grid,
	Front		*front,
	Wave		*wave)
{
	return (grid->step >= stop_step(grid) || grid->time >= stop_time(grid));
}		/*end d_stop_run*/


/*                      time_step():
*
*/

EXPORT	int	time_step(
	CHART		*chart,
	double		dt,
	double		*dt_frac,
	int		step,
	double		time)
{
	Front		*front;
	Grid		*grid;
	Wave		*wave;
	int		(*hyp)(double,double*,Wave*,Front*);
	int		(*make_bubbles)(Wave*,Front*);
	int		status;
	void		(*bc_propagate)(Grid*);
	DEBUG_ENTER(time_step)

	set_current_chart(chart);
	grid = chart->grid;
	wave = chart->wave;
	front = chart->front;
	hyp = chart->hyp_solver;
	bc_propagate = chart->bc_propagate;
	make_bubbles = chart->make_bubbles;

	if (debugging("storage"))
	    long_alloc_view(stdout);

	*dt_frac = 1.0;
	if (chart->parent == NULL && bc_propagate)
	{
	    start_clock("BDRY COND");
	    (*bc_propagate)(grid);		/* Boundary Conditions */
	    stop_clock("BDRY COND");
	}

	/* TODO: change to phase transition */
	if(make_bubbles)
	    (*make_bubbles)(wave, front);

	print_storage("before HYPERBOLIC","HYP_storage");
	start_clock("HYPERBOLIC");
	front->step = step;
 	front->time = wave->time = time;
	front->dt = dt;
	front->dt_frac = dt_frac;
        if (chart->parab)
            front->parab = YES;
        else
            front->parab = NO;
#if defined(SUBGRID)
        if (chart->subgrid)
        {
            if (front->subgrid_time <= time)
                (*chart->subgrid)(dt,front,wave);
        }
#endif /* defined SUBGRID */
	status = (*hyp)(dt,dt_frac,wave,front); /* Hyperbolic Step */
	stop_clock("HYPERBOLIC");
	print_storage("after HYPERBOLIC","HYP_storage");
	if (status != GOOD_STEP)
	{
	   print_time_step_status("WARNING in time_step(),  "
				   "hyp step failed, status = ",status,"\n");
	   DEBUG_LEAVE(time_step)
	   return status;
	}

	start_clock("PARAB");
	if (chart->parab)
	    (*chart->parab)(dt,dt_frac,wave,front); /* Parabolic Step */
	stop_clock("PARAB");

	DEBUG_LEAVE(time_step)
	return status;
}		/*end time_step*/



LOCAL	void set_modified_time_step(
	Grid		*grid,
	Front		*fr,
	Wave		*wave,
	Printplot       *prt,
	double		dt_frac)
{
	double		   dt_last = grid->dt_last;
	DEBUG_ENTER(set_modified_time_step)

	if (grid->num_mts++ >= max_num_time_step_mods(grid))
	{
	    screen("\n\nERROR in set_modified_time_step(), TOO MANY (%d) "
		   "ATTEMPTS TO MODIFY THE TIME STEP WITHOUT SUCCESSFUL TIME "
		   "STEP.\n\n",grid->num_mts);
	    print_interface(fr->interf);
	    clean_up(ERROR);
	}
	fr->num_mts = grid->num_mts;
	if (dt_frac < Min_time_step_modification_factor(fr))
	    dt_frac = Min_time_step_modification_factor(fr);
	if (dt_frac > Max_time_step_modification_factor(fr))
	    dt_frac = Max_time_step_modification_factor(fr);
	(void) printf("Changing time step - old dt %g\n",grid->dt);
	if (grid->pp_grid->nn > 1)
	{
	    double dt_max = dt_frac;
	    double dt_min = dt_frac;

	    pp_global_min(&dt_min,1L);
	    pp_global_max(&dt_max,1L);

	    if (dt_max > 1.00001 && dt_min < 0.9999) /*TOLERANCE*/
	    {
	    	(void) printf("WARNING in time_step():  Need to "
	    	              "both reduce and increase dt!\n"
	    	              "Taking larger of the two.\n");
	    	dt_frac = dt_max;
	    }
	    else
	    	dt_frac = (dt_max > 1.00001) /*TOLERANCE*/ ? dt_max : dt_min;

	}
	grid->dt *= dt_frac;
	if (dt_frac > 1.0)
	{
#if defined(TWOD)
	    if (force_tangle() == NO)
#endif /* defined(TWOD) */
	        grid->dt = nonphysics_timestep_reductions(grid,wave,fr,
						          prt,grid->dt);
	}

	if (!debugging("no_cfl"))
	{
	    double cfl_dt, tmp_dt;
	    double coords[MAXD];
	    double max_dt, min_dt;
	    static const double MIN_DT_FACTOR = 1.0e-5;	/*TOLERANCE*/
	    static const double MAX_DT_FACTOR = 10.0;	/*TOLERANCE*/

	    cfl_dt = (dt_last != 0.0) ? dt_last : grid->dt;
	    if (wave->max_hyp_time_step)
	    {
	    	tmp_dt = (*wave->max_hyp_time_step)(wave,coords);
	    	cfl_dt = min(tmp_dt,cfl_dt);
	    }
	    if (fr->max_front_time_step)
	    {
	    	tmp_dt = (*fr->max_front_time_step)(fr,coords);
	    	cfl_dt = min(tmp_dt,cfl_dt);
	    }
	    min_dt = MIN_DT_FACTOR*cfl_dt;
	    max_dt = MAX_DT_FACTOR*cfl_dt;
	    if (grid->dt < min_dt)
	    	(void) printf("WARNING in set_modified_time_step(), "
	    	              "the new dt = %g, is too small (< %g), "
			      "cfl_dt = %g\n",grid->dt,min_dt,cfl_dt);
	    if (grid->dt > max_dt)
	    	(void) printf("WARNING in set_modified_time_step(), "
	    	              "the new dt = %g, is too big (> %g), "
			      "cfl_dt = %g\n",grid->dt,max_dt,cfl_dt);
	}
	(void) printf(" new dt %g\nModification number = %d\n",
	              grid->dt,grid->num_mts);
	DEBUG_LEAVE(set_modified_time_step)
}		/*end set_modified_time_step*/

#if !(defined(__SUNPRO_C) || defined(__SUNPRO_CC))
#include <time.h>
#include <sys/resource.h>
#endif /* !(defined(__SUNPRO_C) || defined(__SUNPRO_CC)) */

LOCAL	void	print_memory_statistics(
	CHART* root)
{
	static	FILE	*vfile = NULL;

	if (!debugging("vm_a")) return;
	if (vfile == NULL)
	{
	    vfile = fopen("ALLOC_VIEW","w");
	    if (vfile == NULL)
	    	return;
	    setbuf(vfile,NULL);
	}
	(void) foutput(vfile);
	(void) fprintf(vfile,"Long alloc view at time step %d\n",
	                      root->grid->step);
	long_alloc_view(vfile);
	(void) foutput(vfile);
	(void) fprintf(vfile,"End long alloc view at time step %d\n",
		             root->grid->step);
	(void) fprintf(vfile,"\n\n");
}		/*end print_memory_statistics*/
