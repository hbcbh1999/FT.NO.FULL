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
*				hdriver.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains drivers for the hyperbolic library.
*/


#include <hyp/hdecs.h>

	/* LOCAL Function Declarations */
LOCAL	int	hyp_split_driver(double,double*,Wave*,Front*,
				 void(*)(int,int*,double,double,
					 Wave*,Wave*,Front*,Front*,COMPONENT));
LOCAL   void adapt_subdomain(Front*,Wave*,double);

EXPORT int hyp_vector_driver(
	double		dt,
	double		*dt_frac,
	Wave		*wave,
	Front		*front)
{
	int		status;
	DEBUG_ENTER(hyp_vector_driver)

	debug_print("hyp","Entered hyp_vector_driver()\n");

	set_pt_source_interior_vectors(wave);

	front->interf->e_comps = NULL;

	status = hyp_split_driver(dt,dt_frac,wave,front,hyp_reg_vec);

	front->interf->e_comps = NULL;

	free_pt_source_interior_vectors(wave);

	debug_print("hyp","Left hyp_vector_driver()\n");

	DEBUG_LEAVE(hyp_vector_driver)
	return status;
}		/*end hyp_vector_driver*/

EXPORT int hyp_scalar_driver(
	double		dt,
	double		*dt_frac,
	Wave		*wave,
	Front		*front)
{
	int		status;
	DEBUG_ENTER(hyp_scalar_driver)

	debug_print("hyp","Entered hyp_scalar_driver()\n");

	set_pt_source_interior_vectors(wave);

	front->interf->e_comps = NULL;

	status = hyp_split_driver(dt,dt_frac,wave,front,hyp_npt);

	front->interf->e_comps = NULL;

	free_pt_source_interior_vectors(wave);

	debug_print("hyp","Left hyp_scalar_driver()\n");

	DEBUG_LEAVE(hyp_scalar_driver)

	return status;
}		/*end hyp_scalar_driver*/

/*
*			hyp_split_driver():
*/

LOCAL	int hyp_split_driver(
	double		dt,
	double		*dt_frac,
	Wave		*wave,
	Front		*front,
	void		(*sweep)(int,int*,double,double,
				 Wave*,Wave*,Front*,Front*,COMPONENT))
{
	Front		*newfront;
	Front		*infront[3], *outfront[3];
	INTERFACE	*current_intfc, *tmp_intfc;
	Wave		*wk_wv1 = NULL, *wk_wv2 = NULL;
	Wave		*tmpwave = NULL;
	Wave		*inwave[3],  *outwave[3];
	double		dh[MAXD];
	int		step = front->step;
	int		status;
	int		i, k, dim = front->interf->dim;
	int		*iperm;	/* permutation of {0,...,dim-1} */
	COMPONENT	max_comp;
	static char	warn[] = "WARNING in hyp_split_driver()";
	static char	err[] = "ERROR in hyp_split_driver()";
	double           cpu_time;
	DEBUG_ENTER(hyp_split_driver)

	debug_print("hyp","Entered hyp_split_driver()\n");
	set_hyp_npt_globals(wave);

	/* Record start time for adaptive partition */

	if (front->adaptive_partition)
	    cpu_time = cpu_seconds();

		/* Advance Front */

	start_clock("advance_front");
	if (front->_pre_advance_front)
	    pre_advance_front(front); 
	status = advance_front(dt,dt_frac,front,&newfront,(POINTER) wave); 
	stop_clock("advance_front");
	
	if (debugging("hyp"))
	{
	    (void) printf("front->interf\n");
	    print_interface(front->interf);
	    (void) printf("newfront->interf\n");
	    print_interface(newfront->interf);
	}

	if (status != GOOD_STEP)
	{
	    (void) printf("%s advance_front() failed, "
	                  "dt_frac = %g\n",warn,*dt_frac);
	    print_time_step_status("time step status = ",status,"\n");
	    DEBUG_LEAVE(hyp_split_driver)
	    return status;
	}
	initialize_max_wave_speed(wave);

			/* hyperbolic solver */
	if( debugging("hyp_states") )
	{
	    (void) printf("States before calling hyp_solver:\n\n");
	    (void) printf("Old front\n");
	    graph_front_states(front);
	    (void) printf("New front\n");
	    graph_front_states(newfront);
	    (void) printf("wave %p\n",(POINTER)wave);
	    (*wave->show_wave_states)(wave);
	}

	iperm = set_iperm(step,dim);
	for (i = 0; i < dim; ++i)
	{
	    dh[i] = wave->rect_grid->h[iperm[i]];
	}
	if (debugging("x_sweep_only"))
	{
	    for (i = 0; i < dim; ++i)
	    {
	    	iperm[i] = i;
	    	dh[i]	 = wave->rect_grid->h[i];
	    }
	    dim = 1;
	}


		/* Initialize Intermediate Storage for States */

	wk_wv1 = copy_wave(wave);
	clear_wave_pointers(wk_wv1);
	switch (dim)
	{
	case 1:
	    inwave[0] = wave;	outwave[0] = wk_wv1;
	    infront[0] = front;	outfront[0] = newfront;
	    break;
	case 2:
	    if (wave->min_storage == YES)
	    {
	    	inwave[0] = wave;	outwave[0] = wk_wv1;
	    	inwave[1] = outwave[0];	outwave[1] = wave;
	    }
	    else
	    {
	    	wk_wv2 = copy_wave(wave);
	    	clear_wave_pointers(wk_wv2);
	    	inwave[0] = wave;	outwave[0] = wk_wv1;
	    	inwave[1] = outwave[0];	outwave[1] = wk_wv2;
	    }
	    infront[0] = front;	outfront[0] = newfront;
	    infront[1] = newfront;	outfront[1] = newfront;
	    tmpwave = wk_wv1;
	    break;
	case 3:
	    if (wave->min_storage == YES)
	    {
	    	inwave[0] = wave;	outwave[0] = wk_wv1;
	    	inwave[1] = outwave[0];	outwave[1] = wave;
	    	inwave[2] = outwave[1];	outwave[2] = wk_wv1;
	    }
	    else
	    {
	    	wk_wv2 = copy_wave(wave);
	    	clear_wave_pointers(wk_wv2);
	    	inwave[0] = wave;	outwave[0] = wk_wv2;
	    	inwave[1] = outwave[0];	outwave[1] = wk_wv1;
	    	inwave[2] = outwave[1];	outwave[2] = wk_wv2;
	    	tmpwave = wk_wv1;
	    }
	    infront[0] = front;	   outfront[0] = newfront;
	    infront[1] = newfront; outfront[1] = newfront;
	    infront[2] = newfront; outfront[2] = newfront;
	    break;
	}

	start_clock("init_hyp_solution");
	assign_wave_parameters(outwave[0],wave);
	outwave[0]->old_wave = wave;
	status = init_hyp_solution_function(outwave[0],newfront);

	/*printf("#init_hyp_solution_func after\n"); */

	status = syncronize_time_step_status(status,front->pp_grid);
	if (status != GOOD_STEP) 
	{
	    *dt_frac = min(*dt_frac,Min_time_step_modification_factor(front));
	    free_front(newfront);
	    free_wave_pointers(outwave[0]);
	    if (wk_wv1 != NULL)
	    	free_wave(wk_wv1);
	    if (wk_wv2 != NULL)
	    	free_wave(wk_wv2);
	    (void) printf("%s, init_hyp_solution_function() failed\n",warn);
	    DEBUG_LEAVE(hyp_split_driver)
	    return status;
	}
	if (dim > 1 && (wave->min_storage != YES))
	{
	    assign_wave_parameters(outwave[1],outwave[0]);
	    if( !copy_hyp_solution_function(outwave[0],outwave[1]) )
	    {
	    	screen("%s, copy_hyp_solution_function() failed\n",err);
	    	free_front(newfront);
	    	free_wave_pointers(outwave[0]);
	    	if (wk_wv1 != NULL)
	    	    free_wave(wk_wv1);
		if (wk_wv2 != NULL)
		    free_wave(wk_wv2);
		DEBUG_LEAVE(hyp_split_driver)
		return ERROR_IN_STEP;
	    }
	}
	stop_clock("init_hyp_solution");

	start_clock("hyp_solver");


	/*
	*	Temporary storage allocated by store() in directional sweeps 
	*	should be allocated from the tmp_intfc table, since
	*	this storage is freed when tmp_intfc is deleted at
	*	the two sweeps.
	*/

	current_intfc = current_interface();
	tmp_intfc = make_interface(wave->rect_grid->dim);
		/* Call sweep functions in cyclic order */
	max_comp = max_component(front->interf);
	for (k = 0; k < dim; ++k)
	{
	    (*sweep)(k,iperm,dh[k],dt,inwave[k],outwave[k],
		     infront[k],outfront[k],max_comp);

	    if( debugging("hyp_states") )
	    {
	    	(void) printf("After %d%s sweep: outwave[%d] %p\n",
			      k,ordinal_suffix(k),k,(POINTER)outwave[k]);
		if (wave->show_tri_soln)
		    (*wave->show_tri_soln)(outfront[k],outwave[k]);
		(*wave->show_wave_states)(outwave[k]);
	    }
	    
	    /* parallel part for interior states */

            if (!scatter_states(outwave[k],outfront[k],iperm,k))
            {
                screen("%s, scatter_states() failed\n",err);
                clean_up(ERROR);
            }

	    if( debugging("special_plot") && wave->plot_hyp_soln)
	    	(*wave->plot_hyp_soln)(outfront[k],outwave[k],step);
	    if (k == dim-1) break;
	    if (k == 0 && (wave->min_storage == YES))
	    {
	    	free_wave_pointers(outwave[1]);

	    	if (! copy_hyp_solution_function(outwave[0],outwave[1]))
	    	{
	    	    free_front(newfront);
	    	    if (wk_wv1 != NULL)
	    		free_wave(wk_wv1);
	    	    if (wk_wv2 != NULL)
	    		free_wave(wk_wv2);
	    	    screen("%s, copy_hyp_solution_function() failed\n",err);
		    DEBUG_LEAVE(hyp_split_driver)
		    return ERROR_IN_STEP;
		}
	    }
	}
	
	if (dim == 1)
	    free_wave_pointers(wave);
	set_current_interface(current_intfc);
	(void) delete_interface(tmp_intfc);
	


		/* Copy updated front, wave */
	if ((dim % 2) || (wave->min_storage == NO))
	    assign_copy_wave_pointers(wave,outwave[dim-1]);

	/* Free temporary storage, update front */

	if (tmpwave != NULL)
	    free_copy_wave_pointers(tmpwave);

	assign_interface_and_free_front(front,newfront);

	if (wk_wv1 != NULL)
	    free_wave(wk_wv1);
	if (wk_wv2 != NULL)
	    free_wave(wk_wv2);

	stop_clock("hyp_solver");

	if( debugging("hyp_states") )
	{
	    int i;
	    for (i = 0; i < dim; ++i)
	    	(void) printf("sweep %d Maxsp(wave)[%d] %g\n",i,
			      i,Maxsp(wave)[i]);
	    (void) printf("Wave %p after calling hyp_solver:\n",(POINTER)wave);
	    if (wave->show_tri_soln)
		(*wave->show_tri_soln)(front,wave);
	    (*wave->show_wave_states)(wave);
	}

	/* Perform adaptive subdomain partition */
	if (front->adaptive_partition)
	{
	    cpu_time = cpu_seconds() - cpu_time;
	    adapt_subdomain(front,wave,cpu_time);
	}

	debug_print("hyp","Left hyp_split_driver()\n");
	DEBUG_LEAVE(hyp_split_driver)
	return status;
}		/*end hyp_split_driver*/

LOCAL	void adapt_subdomain(
	Front *front,
	Wave *wave,
	double cpu_time)
{
	int 		i,j,k,id;
	int		dim = front->rect_grid->dim;
	int 		status;
	int 		icoords[MAXD];
	int  		*iperm;	
	Locstate	s_old,s_new;
	int		lexpand[MAXD],uexpand[MAXD];
	int         	ix, iy, iz;
	int         	xmax, ymax, zmax;
	Wave		*tmp_wv;

	if (!cpu_adapt_front(front,cpu_time,lexpand,uexpand)) 
	    return;
	tmp_wv = copy_wave(wave);
	status = init_hyp_solution_function(tmp_wv,front);
	status = syncronize_time_step_status(status,front->pp_grid);
	if (status != GOOD_STEP)
	{
	    screen("init_hyp_solution_function() failed\n");
	    clean_up(ERROR);
	}
	switch (dim)
	{
	case 2:
	    xmax = tmp_wv->rect_grid->gmax[0];
	    ymax = tmp_wv->rect_grid->gmax[1];
            for (iy = 0; iy < ymax; ++iy)
            {
                for (ix = 0; ix < xmax; ++ix)
                {
                    icoords[0] = ix;
                    icoords[1] = iy;
		    s_new = Rect_state(icoords,tmp_wv);

		    if (lexpand[0] == -1) icoords[0] = ix + 1;
		    else if (lexpand[0] == 1) icoords[0] = ix - 1;
		    if (lexpand[1] == -1) icoords[1] = iy + 1;
		    else if (lexpand[1] == 1) icoords[1] = iy - 1;
		    s_old = Rect_state(icoords,wave);

		    ft_assign(s_new,s_old,front->sizest);
                }
            }
	    break;
	case 3:
            xmax = tmp_wv->rect_grid->gmax[0];
            ymax = tmp_wv->rect_grid->gmax[1];
   	    zmax = tmp_wv->rect_grid->gmax[2];
            for (iz = 0; iz < zmax; ++iz)
	    {
	        for (iy = 0; iy < ymax; ++iy)
                {
                    for (ix = 0; ix < xmax; ++ix)
                    {
                        icoords[0] = ix;
                        icoords[1] = iy;
			icoords[2] = iz;
                        s_new = Rect_state(icoords,tmp_wv);

                        if (lexpand[0] == -1) icoords[0] = ix + 1;
                        else if (lexpand[0] == 1) icoords[0] = ix - 1;
                        if (lexpand[1] == -1) icoords[1] = iy + 1;
                    	else if (lexpand[1] == 1) icoords[1] = iy - 1;
			if (lexpand[2] == -1) icoords[2] = iz + 1;
			else if (lexpand[2] == 1) icoords[2] = iz - 1;
                    
		 	s_old = Rect_state(icoords,wave);
			ft_assign(s_new,s_old,front->sizest);
		    }
                }
            }
	    break;
	}
	iperm = set_iperm(front->step,dim);
	for (i = 0; i < dim; ++i)
	{
            scatter_states(tmp_wv,front,iperm,i);
	}
	free_wave_pointers(wave);
	assign_copy_wave_pointers(wave,tmp_wv);
}	/* end adapt_subdomain */
