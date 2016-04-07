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
*				hoverture_driver.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains drivers for the hyperbolic library.
*/

#define DEBUG_STRING    "hover_driver"
#include <hyp/hdecs.h>

#if defined(USE_OVERTURE)
       /* LOCAL Variable Declarations */  
/* Trans_wv data structure is used to recreate
   Wave after it is sent off to other processors. 
   Trans_wv records the Wave size cut in the FronTier 
   buffer zone. 
*/  
struct _Trans_wv{
        int    wv_level;
        int    wv_id;         /* The wv_id this region is cut from */ 
        int    off_set[MAXD]; /* subdomain's [0,0] index in global domain*/
        int    base[MAXD];  /* cut local grid location index   */
        int    bound[MAXD]; /* cut local grid location index   */
};

typedef struct _Trans_wv Trans_wv;

	/* LOCAL Function Declarations */
LOCAL   Trans_wv  *pp_send_patch_interior(int*,int,int*,int,
                    Overparam*,Wv_on_pc**,int,Front**,Wave**,int,int*); 
LOCAL   Trans_wv  *pp_receive_patch_interior(int*,int,int*,
                     int,Overparam*,Wv_on_pc**,Front**,Wave**,int,int*); 
LOCAL   void    set_send_domain(int*,int*,int*,int,int,RECT_GRID*); 
LOCAL   void    set_receive_domain(int*,int*,int*,int,int,RECT_GRID*); 
LOCAL   void    set_send_patch_domain(int*,int*,int*,int,int,
                     int*,int*,int*,int*,int*,int*,RECT_GRID*); 
LOCAL   void    set_receive_patch_domain(int*,int*,int*,int,int,
                     int*,int*,int*,int*,int*,int*,RECT_GRID*); 
LOCAL   void    create_recv_buffer_zone(Wave***,Front***,
                     Trans_wv*,int,Overparam*,Wave*,Front*); 
LOCAL   int     init_buffer_hyp_solution_function(Wave*,Front*); 
LOCAL   void    pp_send_patch_interior_states(int*,int,int*,
                     int,Front**,Wave**,int,Wv_on_pc**,Trans_wv*,int); 
LOCAL   void    pp_receive_patch_interior_states(int*,int,int*,
                     int,Overparam*,Wv_on_pc**,Front**,Wave**,int,
                     Front**,Wave**,Trans_wv*,int); 
LOCAL   void    perform_recv_scatter_patch_states(Overparam*,Wv_on_pc**,
                 Wave**,Front**,Wave**,Front**,Trans_wv*,int*,int,int); 
LOCAL   void    fill_recv_scatter_patch_obstacle_states(Overparam*,
                 Wv_on_pc**,Wave**,Front**,Wave**,Front**,Trans_wv*,int*,
                 int,int); 
LOCAL   int     match_point_in_zone(int*,int*,int*,int,
                  Trans_wv*,int,int*,int*); 

LOCAL   void    bundle_comp(int*,int*,Wave*,byte*); 
LOCAL   void    unbundle_comp(int*,int*,Wave*,byte*); 
LOCAL   int     make_reinstore_patch_storage(Wave***,Front***,int*,
                  Wave**,Front**,Wv_on_pc**,int,int); 

EXPORT	int hyp_amr_split_driver(
	double			      dt,
	double			      *dt_frac,
	Wave			      **wvs,
	Front			      **frs,
	Front	        	      **nfrs,
        int                           num_patches,
        Wv_on_pc                      **redistr_table,
        int                           max_n_patch,
        CompositeGrid                 *cg_over,
	doubleCompositeGridFunction   *cg_over_function0,
	doubleCompositeGridFunction   *cg_over_function,
        Overparam                     *overparam,
	void		(*sweep)(int,int*,double,double,
			 Wave*,Wave*,Front*,Front*,COMPONENT))
{
        Wave            **wk_wv1;                    
        int             status; 
        int             i, k, dim = frs[0]->interf->dim;
        int             *iperm; /* permutation of {0,...,dim-1} */      
        Front           *infront[3], *outfront[3]; 
        Wave            *inwave[3],  *outwave[3];  
        Wave            *tmpwave = NULL;
        double           dh[MAXD];  
        COMPONENT       max_comp;
        INTERFACE       *current_intfc, *tmp_intfc;
        Wave            **undistr_wvs;
        Front           **undistr_frs;
        int             resid_n; 
        int             numnodes, myid; 

        myid = pp_mynode();
        numnodes = pp_numnodes();

        uni_array(&wk_wv1,num_patches,sizeof(Wave*)); 
        for (i = 0; i < num_patches; i++)
        {
            wk_wv1[i] = copy_wave(wvs[i]);  clear_wave_pointers(wk_wv1[i]); 
            assign_wave_parameters(wk_wv1[i],wvs[i]); 
            wk_wv1[i]->cg_over_function = cg_over_function0;
            wk_wv1[i]->cg_over = cg_over;
            status = init_hyp_solution_function(wk_wv1[i],nfrs[i]); 
            if(status != GOOD_STEP)
            {
                printf("ERROR: hyp_amr_split_driver\n");
                printf("init_hyp_solution_function failed on %d\n", i);
                printf("returned status = %d\n", status);  
                clean_up(ERROR);  
            } 
	    initialize_max_wave_speed(wvs[i]);
        } 
        
        iperm = set_iperm(frs[0]->step,dim); 

        current_intfc = current_interface();
        tmp_intfc = make_interface(wvs[0]->rect_grid->dim);
       
        for(k = 0; k < dim; k++)
        {
            if (debugging("hyp_split_amr"))
                printf("Doing hyp_amr_split_driver swp[%d] in dir[%d]\n",
                      k,iperm[k]);
            for (i = 0; i < num_patches; i++)
            {
	        set_hyp_npt_globals(wvs[i]);
                dh[k] = wvs[i]->rect_grid->h[iperm[k]]; 
                if (debugging("hyp_split_amr"))
                    printf("Doing hyp_amr_split_driver Front[%d] level[%d] swp[%d] in dir[%d]\n",
                          frs[i]->patch_number, frs[i]->patch_level,k,iperm[k]);
                switch (dim)
                {
                case 1:
                case 3:
                        printf("ERROR: hyp_amr_split_driver\n");
                        printf("Implement 1,3 option\n");
                        clean_up(ERROR);  
                break;
                case 2:
                    if (wvs[i]->min_storage)
                    {
                        inwave[0] = wvs[i];     outwave[0] = wk_wv1[i];
                        inwave[1] = outwave[0]; outwave[1] = wvs[i];
                    }
                    else
                    {
                        printf("ERROR: hyp_amr_split_driver\n");
                        printf("Implement min_storage option\n");
                        clean_up(ERROR);  
                    }
                    infront[0] = frs[i];    outfront[0] = nfrs[i];
                    infront[1] = nfrs[i];   outfront[1] = nfrs[i];
                    tmpwave = outwave[0];
                break;
                }  

                max_comp = max_component(frs[i]->interf); 
                outwave[0]->cg_over_function = cg_over_function0;
                outwave[0]->cg_over = cg_over;

                (*sweep)(k,iperm,dh[k],dt,inwave[k],outwave[k],
                     infront[k],outfront[k],max_comp);  
                /*     
                if(i == 0) 
                {
                    printf("Wave[%d] SHOW WAVE STATES after [%d] sweep\n",0,k);
                    (*outwave[k]->show_wave_states)(outwave[k]);
                }
                */  
            } 
            switch (dim)
            {
            case 2:
                if(k == 0) 
                {
                    reinstore_undistribute_patch(&undistr_wvs,&undistr_frs,wk_wv1,nfrs,
                       redistr_table,num_patches,max_n_patch); 
                    h_scatter_patch_states(undistr_wvs, undistr_frs,
                       overparam, redistr_table, max_n_patch, iperm, k); 
                    (*wk_wv1[0]->scatter_patch_states_in_sweep_dir)(overparam,
                         redistr_table, undistr_wvs, undistr_frs,iperm,k);
                    resid_n = 0;
                    for(int tmpi = 0; tmpi < max_n_patch; tmpi++)
                    {
                        if(-1 == redistr_table[pp_mynode()][tmpi].wv_id) continue;
                        if(pp_mynode() == redistr_table[pp_mynode()][tmpi].pc_id)
                            resid_n++;
                    }
                    patch_wave_trans(undistr_wvs,wk_wv1,redistr_table, 
                        myid, resid_n, max_n_patch, numnodes); 
                    free_copy_patches(undistr_wvs,undistr_frs,redistr_table,max_n_patch); 
                    free_these(2,undistr_wvs,undistr_frs); 
                } 
                else
                {
                    reinstore_undistribute_patch(&undistr_wvs,&undistr_frs,wvs,nfrs,
                      redistr_table,num_patches,max_n_patch);
                    /* set the interior states of the coarse new patches.
                     * It is the FronTier interp. function,
                     * which ft_assigns fine state average to coarse. 
                     */
                    /* 051603: add interpolation_fine_to_coarse()
                     * after the last directional sweep, before scatter_states.
                     * This will save communication and provide better symmetry.
                     */ 
                    (*wvs[0]->overture_undistribute_interpolation_fine_to_coarse)(
                        undistr_wvs,undistr_frs); 
                    h_scatter_patch_states(undistr_wvs, undistr_frs,
                       overparam, redistr_table, max_n_patch, iperm, k); 
                    (*wvs[0]->scatter_patch_states_in_sweep_dir)(overparam,redistr_table,
                       undistr_wvs, undistr_frs,iperm,k);
                    resid_n = 0;
                    for(int tmpi = 0; tmpi < max_n_patch; tmpi++)
                    {
                        if(-1 == redistr_table[pp_mynode()][tmpi].wv_id) continue;
                        if(pp_mynode() == redistr_table[pp_mynode()][tmpi].pc_id)
                            resid_n++;
                    }
                    patch_wave_trans(undistr_wvs,wvs,redistr_table, 
                        myid, resid_n, max_n_patch, numnodes); 
                    free_copy_patches(undistr_wvs,undistr_frs,redistr_table,max_n_patch); 
                    free_these(2,undistr_wvs,undistr_frs); 
                } 
            break;
            }  
            if (k == 0 and wvs[0]->min_storage)
            {
                for(i = 0; i < num_patches; i++)
                {
                    free_wave_pointers(wvs[i]);

                    wvs[i]->cg_over_function = cg_over_function;
                    wvs[i]->cg_over = wave_tri_soln(wk_wv1[i])->cg_over;
                    wvs[i]->patch_number = wave_tri_soln(wk_wv1[i])->patch_number;
                    wvs[i]->use_overture_state = wave_tri_soln(wk_wv1[i])->use_overture_state;
                    wvs[i]->overture_init_step = wave_tri_soln(wk_wv1[i])->overture_init_step;
                    wvs[i]->patch_level = wave_tri_soln(wk_wv1[i])->patch_level;
                    wvs[i]->patch_component = wave_tri_soln(wk_wv1[i])->patch_component;
                    wvs[i]->NumberOfLevels = wave_tri_soln(wk_wv1[i])->NumberOfLevels;

                    if (! copy_hyp_solution_function(wk_wv1[i],wvs[i]))
                    {
                        deep_free_front(nfrs[i]);
                        screen("%s, copy_hyp_solution_function() failed for[%d]\n",
                            "ERROR: in hyp_single_patch_split_driver()",i);
                        DEBUG_LEAVE(hyp_amr_split_driver)
                        return ERROR_IN_STEP;
                    }
                }
            }
        } 

        set_current_interface(current_intfc);
        (void) delete_interface(tmp_intfc);  

        if (debugging("hyp_split_amr"))
            printf("Doing hyp_amr_split_driver clear wk_wvs\n");

        for (i = 0; i < num_patches; i++)
        {
            newfront_to_distri_table(nfrs[i],frs[i],
                     num_patches,redistr_table,max_n_patch);
            assign_interface_and_free_front(frs[i],nfrs[i]);
            if (wk_wv1[i] != NULL)
            {
                /* 
                printf("wk_wv1[%d] = %d, wvs[%d] = %d\n", i, 
                   wave_tri_soln(wk_wv1[i]), i, wave_tri_soln(wvs[i])); 
                */  
                free_copy_wave_pointers(wk_wv1[i]); 
                free_wave(wk_wv1[i]);
            } 
        }  
        free(wk_wv1);  
        return GOOD_STEP;  
} 

EXPORT void free_copy_patches(
	Wave        **wvs,
        Front       **frs,
        Wv_on_pc    **redistr_table,
        int         max_n_patch)
{
        int        i, myid; 
        int        patch_id;
        INTERFACE  *current_intfc; 

	DEBUG_ENTER(free_copy_patches)
 
        myid = pp_mynode(); 

        for (i = 0; i < max_n_patch; i++)
        {
            if(-1 == redistr_table[myid][i].wv_id) continue;
            if(myid == redistr_table[myid][i].pc_id) continue;

            patch_id = redistr_table[myid][i].wv_id;
            deep_free_front(frs[patch_id]);
            free_wave_pointers(wvs[patch_id]);
            free_wave(wvs[patch_id]);
        }
	DEBUG_LEAVE(free_copy_patches)
}


/*
*		hyp_patch_split_driver():
*/

EXPORT	int hyp_patch_split_driver(
	double				dt,
	double				*dt_frac,
	Wave				*wave,
	Front				*front,
	Front	        		*newfront,
        CompositeGrid                   *cg_over,
	doubleCompositeGridFunction     *cg_over_function0,
	doubleCompositeGridFunction     *cg_over_function,
	void			(*sweep)(int,int*,double,double,
				 Wave*,Wave*,Front*,Front*,COMPONENT))
{
	Front		*infront[3], *outfront[3];
	INTERFACE	*current_intfc;
	Wave		*wk_wv1 = NULL, *wk_wv2 = NULL;
	Wave		*tmpwave = NULL;
	Wave		*inwave[3],  *outwave[3];
	double		dh[MAXD];
	int		step = front->step;
	int		status;
	int		i, k, dim = front->interf->dim;
	int		*iperm;	/* permutation of {0,...,dim-1} */
	COMPONENT	max_comp;
	static char	warn[] = "WARNING in hyp_patch_split_driver()";
	static char	err[] = "FT_ERROR in hyp_patch_split_driver()";

	DEBUG_ENTER(hyp_patch_split_driver)

	set_hyp_npt_globals(wave);

	initialize_max_wave_speed(wave);

	if( debugging("hyp_patch_split_driver") ) 
	{	
	    printf("In hyp_patch_split_driver(), working in patch %d\n",
			wave->patch_number);
	}      

		/* hyperbolic solver */
	if( debugging("hyp_states") )
	{
		(void) printf("States before calling hyp_solver:\n\n");
		(void) printf("Old front\n");	  graph_front_states(front);
		(void) printf("New front\n");	  graph_front_states(newfront);
		(void) printf("wave %p\n",(POINTER)wave);
		(*wave->show_wave_states)(wave);
	}

	iperm = set_iperm(step,dim);
	for (i = 0; i < dim; i++)
	{
	    dh[i] = wave->rect_grid->h[iperm[i]];
	}
	if (debugging("x_sweep_only"))
	{
	    for (i = 0; i < dim; i++)
	    {
		iperm[i] = i;
		dh[i] = wave->rect_grid->h[i];
	    }
	    dim = 1;
	}

		/* Initialize Intermediate Storage for States */

	wk_wv1 = copy_wave(wave);	clear_wave_pointers(wk_wv1);
	switch (dim)
	{
	case 1:
		inwave[0] = wave;	outwave[0] = wk_wv1;
		infront[0] = front;	outfront[0] = newfront;
		break;
	case 2:
		if (wave->min_storage)
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
		if (wave->min_storage)
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
		infront[0] = front;	outfront[0] = newfront;
		infront[1] = newfront;	outfront[1] = newfront;
		infront[2] = newfront;	outfront[2] = newfront;
		break;
	}


	start_clock("init_hyp_solution");
	assign_wave_parameters(outwave[0],wave);

        outwave[0]->cg_over_function = cg_over_function0; 
	outwave[0]->cg_over = cg_over;

	if(debugging("hyp_patch_split_driver"))
	{
	    printf("In hyp_patch_split_driver(), patch = %d\n",
		                      newfront->patch_number);
	    printf("print interface before init_hyp_solution_function()\n");  
	    print_interface(newfront->interf); 
	}  	

	status = init_hyp_solution_function(outwave[0],newfront);
	if(debugging("hyp_patch_split_driver"))
	{
	    printf("In hyp_patch_split_driver(), patch = %d\n", 
	          newfront->patch_number);		
            int ic[2] = {0,0};
	    printf("IN hyp_patch_split_driver()\n"
	           " graped crds of ic<%d, %d> <%g, %g>\n",ic[0], ic[1],
	        Rect_coords(ic,outwave[0])[0],Rect_coords(ic,outwave[0])[1]);
	    print_rectangular_grid(&wave_tri_soln(outwave[0])->tri_grid->tg_grid);
	    print_rectangular_grid(&wave_tri_soln(outwave[0])->tri_grid->rect_grid);
	    print_rectangular_grid(outwave[0]->rect_grid);
/* 
	    print_components(wave_tri_soln(outwave[0])->tri_grid);  
*/  
	}
/*   
	status = syncronize_time_step_status(status,front->pp_grid);
*/   
	if (status != GOOD_STEP) 
	{
	    *dt_frac = min(*dt_frac,
			       Min_time_step_modification_factor(front));
	    free_front(newfront);
	    free_wave_pointers(outwave[0]);
	    if (wk_wv1 != NULL)
	 	free_wave(wk_wv1);
	    if (wk_wv2 != NULL)
		free_wave(wk_wv2);
	    (void) printf("%s, init_hyp_solution_function() failed\n",
			      warn);
	    DEBUG_LEAVE(hyp_patch_split_driver)
	    return status;
	}
	if (dim > 1 and not wave->min_storage)
	{
            assign_wave_parameters(outwave[1],outwave[0]);
	    if( ! copy_hyp_solution_function(outwave[0],outwave[1]) )
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

	current_intfc = current_interface();

		/* Call sweep functions in cyclic order */
	max_comp = max_component(front->interf);
	for (k = 0; k < dim; k++)
	{
	    (*sweep)(k,iperm,dh[k],dt,inwave[k],outwave[k],
		infront[k],outfront[k],max_comp); 

	    if( debugging("hyp_states") )
	    {
		(void) printf("After %d%s sweep: outwave[%d] %p\n",
			      k,ordinal_suffix(k),k,
			      (POINTER)outwave[k]);
		if (wave->show_tri_soln)
			(*wave->show_tri_soln)(outfront[k],outwave[k]);
		(*wave->show_wave_states)(outwave[k]);
	    }

	    if( debugging("special_plot") and wave->plot_hyp_soln)
		(*wave->plot_hyp_soln)(outfront[k],outwave[k],step);
	    if (k == dim-1) break;
	    if (k == 0 and wave->min_storage)
	    {
		free_wave_pointers(outwave[1]);

	        outwave[1]->cg_over_function = cg_over_function;
                outwave[1]->cg_over = wave_tri_soln(outwave[0])->cg_over;
	        outwave[1]->patch_number = wave_tri_soln(outwave[0])->patch_number;
	        outwave[1]->use_overture_state = wave_tri_soln(outwave[0])->use_overture_state;
	        outwave[1]->overture_init_step = wave_tri_soln(outwave[0])->overture_init_step;
	        outwave[1]->patch_level = wave_tri_soln(outwave[0])->patch_level;
	        outwave[1]->patch_component = wave_tri_soln(outwave[0])->patch_component;
	        outwave[1]->NumberOfLevels = wave_tri_soln(outwave[0])->NumberOfLevels;

		if (! copy_hyp_solution_function(outwave[0],outwave[1]))
		{
		    free_front(newfront);
		    if (wk_wv1 != NULL)
			free_wave(wk_wv1);
		    if (wk_wv2 != NULL)
			free_wave(wk_wv2);
		    screen("%s, copy_hyp_solution_function() failed\n",
			err);
		    DEBUG_LEAVE(hyp_split_driver)
		    return ERROR_IN_STEP;
		}
	    }
	}
	if (dim == 1)
		free_wave_pointers(wave);
	set_current_interface(current_intfc);

		/* Copy updated front, wave */
	if ((dim % 2) or (wave->min_storage == NO))
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

	/* parallel part for interior states */
/*    
	if (not scatter_states(wave,front))
	{
	    screen("%s, scatter_states() failed\n",err);
	    clean_up(FT_ERROR);
	}
*/        

	if( debugging("hyp_states") )
	{
	    int i;
	    for (i = 0; i < dim; i++)
		(void) printf("sweep %d Maxsp(wave)[%d] %g\n",
			      i,i,Maxsp(wave)[i]);
	    (void) printf("Wave %p after calling hyp_solver:\n",
			      (POINTER)wave);
	    if (wave->show_tri_soln) (*wave->show_tri_soln)(front,wave);
	    (*wave->show_wave_states)(wave);
	}

	DEBUG_LEAVE(hyp_patch_split_driver)
	return status;
}		/*end hyp_patch_split_driver*/

EXPORT TRI_SOLN *copy_AMR_tri_soln_storage(
        TRI_SOLN *osoln,
        Wave     *nwave)
{
        TRI_SOLN *nsoln;
        TRI_GRID *grid;

        scalar(&nsoln,sizeof(TRI_SOLN));
        grid = allocate_tri_grid(&osoln->tri_grid->tri_grid_hooks);
        if (grid == NULL)
        {
            free(nsoln);
            return NULL;
        }
        nsoln->cg_over = nwave->cg_over;
        nsoln->patch_number = nwave->patch_number;
        nsoln->use_overture_state = nwave->use_overture_state;
        nsoln->overture_init_step = nwave->overture_init_step;
        nsoln->cg_over_function = nwave->cg_over_function;
        nsoln->patch_component = nwave->patch_component;
        nsoln->patch_level = nwave->patch_level;
        nsoln->NumberOfLevels = nwave->NumberOfLevels;

        grid->cg_over = nwave->cg_over;
        grid->patch_number = nwave->patch_number;
        grid->use_overture_state = nwave->use_overture_state;
        grid->overture_init_step = nwave->overture_init_step;
        grid->cg_over_function = nwave->cg_over_function;
        grid->patch_component = nwave->patch_component;
        grid->patch_level = nwave->patch_level;
        grid->NumberOfLevels = nwave->NumberOfLevels;

        copy_tri_grid(osoln->tri_grid,grid,osoln->sizest);

        set_tri_soln_struct(nsoln,osoln->intfc,grid,osoln->sizest,
                   &osoln->interpolator,&osoln->el_integral,&osoln->unsplit);
        return nsoln;
}        /*end copy_AMR_tri_soln_storage*/

EXPORT int reinstore_mini_undistribute_patch(
        Wave       ***newwvs,
        Front      ***newfrs,  
        int        *resid_n, 
        Wave       **wvs,  
        Front      **frs,
        Wv_on_pc   **redistr_table,
        int        num_patches,  /* number of patches compute in this proc */
        int        max_n_patch,
        int        set_comp)
{
        int        source, numnodes, myid;
        int        dist;
        int        patch_id, i;
        Front      **tmpfrs, *front;
        Front      *basefr;
        Wave       **tmpwvs, *wave;   
        int        total_patch;
        int        nn, bal_n, mm;
        INTERFACE  *current_intfc;
        COMPONENT  dummy = -1;
        RECT_GRID  *cgr, *tgr;
        RECT_GRID  *rgr;        
        size_t      len;
        byte        *storage = NULL, *buf;
        POINTER     info;
        size_t      sizest = wvs[0]->sizest;
        int         ii, jj, L[MAXD], U[MAXD];
        int         *fly_patch_id;
 
        DEBUG_ENTER(reinstore_mini_undistribute_patch) 

        numnodes = pp_numnodes();
        myid = pp_mynode();

        total_patch = make_reinstore_patch_storage(&tmpwvs,&tmpfrs,
            &nn,wvs,frs,redistr_table,num_patches,max_n_patch);

        uni_array(&fly_patch_id, total_patch, sizeof(int));
        for(i = 0; i < total_patch; i++)
            fly_patch_id[i] = -1;         

        *resid_n = bal_n = nn;
        /*  recv all fronts originally in this procs  */
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                patch_id = -100;
                if(-1 == redistr_table[source][i].wv_id) continue;
                if(source == redistr_table[source][i].pc_id) continue;

                if(myid == redistr_table[source][i].pc_id)
                {
                    patch_id = redistr_table[source][i].wv_id;
                    front = redistr_table[source][i].front;

                    pp_send(0,(POINTER)(&patch_id),sizeof(int),source);
                    dummy = wave_of_front(front)->patch_component;
                    send_mini_front_misc(front,&dummy,source);
                    /*    
                    printf("redist[%d][%d] should be back to original: ",source,i);
                    printf("front %p: Proc [%d] Send to proc[%d]\n",
                            front, myid, source);
                    */   
                }
                if(source == myid)
                {
                    pp_recv(0,redistr_table[source][i].pc_id,
                           (POINTER)(&patch_id),sizeof(int));

                    recv_mini_front_misc(tmpfrs[patch_id],
                        &dummy, redistr_table[source][i].pc_id);
                    wave_of_front(tmpfrs[patch_id])->patch_component = dummy;
                    fly_patch_id[patch_id] = patch_id;
                    nn++;
                    /*  
                    printf("redist[%d][%d] received at original: ",source,i);
                    printf("patch[%d]: Proc [%d] recv from proc[%d]\n",
                            patch_id, source, redistr_table[source][i].pc_id);
                    */ 
                }
                pp_gsync();
            }
        }

        /* set recving mini patches, only wave states are transfered */
        for(i = 0; i < total_patch; i++)
        {
            if(fly_patch_id[i] == -1) continue;

            tmpfrs[i]->interf->modified = YES;
            rgr = tmpfrs[i]->rect_grid;
            cgr = computational_grid(tmpfrs[i]->interf);
            copy_rect_grid(cgr,rgr);
            tgr = &topological_grid(tmpfrs[i]->interf);
            tgr->Remap.remap = rgr->Remap.remap;
            set_patch_topo_grid(rgr,tgr); 

            tmpwvs[i]->rect_grid = tmpfrs[i]->rect_grid;
            tmpwvs[i]->pd_flag = tmpfrs[i]->pd_flag;
            tmpwvs[i]->patch_number = tmpfrs[i]->patch_number;
            tmpwvs[i]->patch_level = tmpfrs[i]->patch_level;   
            if (init_buffer_hyp_solution_function(tmpwvs[i],tmpfrs[i]) != GOOD_STEP)
            {
                screen("ERROR: reinstore_undistribute_patch(), "
                       "init_hyp_solution_function() failed for [%d]\n",i);
                clean_up(ERROR);
            } 
        }  

        nn = bal_n;  
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                patch_id = -100;
                if(-1 == redistr_table[source][i].wv_id) continue;
                if(source == redistr_table[source][i].pc_id) continue;

                if(myid == redistr_table[source][i].pc_id)
                {
                    patch_id = redistr_table[source][i].wv_id;
                    wave = NULL;  
                    for(jj = 0; jj < num_patches; jj++) 
                    {
                        if(frs[jj] == redistr_table[source][i].front)
                        {
                            wave = wvs[jj]; 
                            break;  
                        } 
                    } 
                    if(NULL == wave)
                    {
                        printf("ERROR: reinstore_mini_undistribute_patch\n");
                        printf("SEND wave is NULL\n");
                        clean_up(ERROR);  
                    } 
                    pp_send(0,(POINTER)(&patch_id),sizeof(int),source);
                    if(YES == set_comp)
                        len = wave_tri_soln(wave)->tri_grid->n_reg_nodes*
                            (sizest+sizeof(COMPONENT));  
                    else
                        len = wave_tri_soln(wave)->tri_grid->n_reg_nodes*sizest;  
                    scalar(&storage,len);  
                    for(jj = 0; jj < wave->rect_grid->dim; jj++)
                    {
                        L[jj] = -(wave->rect_grid->lbuf[jj]);
                        U[jj] = wave->rect_grid->gmax[jj] 
                              + wave->rect_grid->ubuf[jj];
                    }
                    (*wave->bundle_states)(L,U,wave,storage); 
                    if(YES == set_comp)
                    {
                        buf = storage + 
                           wave_tri_soln(wave)->tri_grid->n_reg_nodes*sizest; 
                        bundle_comp(L,U,wave,buf);
                    }
                    pp_send_large_data(storage,len,source);       
                    free(storage);  
                    /*  
                    printf("redist[%d][%d]  should be back to original: ",source,i);
                    printf("wave %p: Proc [%d] Send to proc[%d]\n",
                             wave, myid, source);
                    */  
                }
                if(source == myid)
                {
                    pp_recv(0,redistr_table[source][i].pc_id,
                             (POINTER)(&patch_id),sizeof(int));
                    wave = tmpwvs[patch_id];
                    if(YES == set_comp)
                        len = wave_tri_soln(wave)->tri_grid->n_reg_nodes
                               *(sizest+sizeof(COMPONENT));  
                    else 
                        len = wave_tri_soln(wave)->tri_grid->n_reg_nodes*sizest;  
                    scalar(&storage,len);
                    pp_receive_large_data(storage,len,redistr_table[source][i].pc_id); 
                    for(jj = 0; jj < wave->rect_grid->dim; jj++)
                    {
                        L[jj] = -(wave->rect_grid->lbuf[jj]);
                        U[jj] = wave->rect_grid->gmax[jj] 
                              + wave->rect_grid->ubuf[jj];
                    }
                    (*wave->unbundle_states)(L,U,wave,storage);
                    if(YES == set_comp)
                    {
                        buf = storage +
                         wave_tri_soln(wave)->tri_grid->n_reg_nodes*sizest;  
                        unbundle_comp(L,U,wave,buf);
                    }  
                    free(storage);
                    nn++;
                    /*                 
                    printf("redist[%d][%d] received at original: ",source,i);
                    printf("patch[%d]: Proc [%d] recv from proc[%d]\n",
                            patch_id, source, redistr_table[source][i].pc_id);
                    */              
                }
                pp_gsync();
            }
        }

        *newwvs = tmpwvs;
        *newfrs = tmpfrs; 
 
        free(fly_patch_id);

        DEBUG_LEAVE(reinstore_mini_undistribute_patch) 
        return GOOD_STEP;  
} 

LOCAL void bundle_comp(
        int             *gmin,
        int             *gmax,
        Wave            *wave,
        byte            *buf)
{
        int       dim = wave->rect_grid->dim;
        int       ic[MAXD];
        POINTER   info;  
      
        DEBUG_ENTER(bundle_comp)

        switch (dim)
        {
        case 1:
        case 3:
            printf("ERROR: bundle_comp\n");
            printf("Implement dim[%d] case\n", dim);
            clean_up(ERROR);  
        break; 
        case 2:
        {
            int        ix, ixmin, ixmax;
            int        iy, iymin, iymax;

            ixmin = gmin[0];    ixmax = gmax[0];
            iymin = gmin[1];    iymax = gmax[1];
            for (iy = iymin; iy < iymax; ++iy)
            {
                ic[1] = iy;
                for (ix = ixmin; ix < ixmax; ++ix)
                {
                    ic[0] = ix;
                    info = (POINTER) buf;
                    ft_assign(info,&Rect_comp(ic,wave),sizeof(COMPONENT));
                    buf += sizeof(COMPONENT);
                }
            }
        } 
        break;
        } /* end of switch */  

        DEBUG_LEAVE(bundle_comp)
}  

LOCAL void unbundle_comp(
        int             *gmin,
        int             *gmax,
        Wave            *wave,
        byte            *buf)
{
        int       dim = wave->rect_grid->dim;
        int       ic[MAXD];
        POINTER   info;  

        DEBUG_ENTER(unbundle_comp)

        switch (dim)
        {
        case 1:
        case 3:
            printf("ERROR: unbundle_comp\n");
            printf("Implement dim[%d] case\n", dim);
            clean_up(ERROR);  
        break; 
        case 2:
        {
            int        ix, ixmin, ixmax;
            int        iy, iymin, iymax;

            ixmin = gmin[0];    ixmax = gmax[0];
            iymin = gmin[1];    iymax = gmax[1];
            for (iy = iymin; iy < iymax; ++iy)
            {
                ic[1] = iy;
                for (ix = ixmin; ix < ixmax; ++ix)
                {
                    ic[0] = ix;
                    info = (POINTER) buf;
                    ft_assign(&Rect_comp(ic,wave),info,sizeof(COMPONENT));
                    buf += sizeof(COMPONENT);
                }
            }
        } 
        break;
        } /* end of switch */  

        DEBUG_LEAVE(unbundle_comp)
}  

LOCAL int make_reinstore_patch_storage(
        Wave       ***newwvs,
        Front      ***newfrs,  
        int        *resid_n, 
        Wave       **wvs,  
        Front      **frs,
        Wv_on_pc   **redistr_table,
        int        num_patches,    /* number of patches compute in this proc */
        int        max_n_patch)
{
        int        i, myid;
        Front      *basefr;
        int        total_patch;
        int        nn;
        int        patch_id;  
        Wave       *wave;  
 
        DEBUG_ENTER(make_reinstore_patch_storage) 

        myid = pp_mynode();

        /* make sure current_intfc is always on patch[0]  */
        if(NULL == frs[0]->interf)
        {
            printf("ERROR: assembly_distribute_patch_fronts()\n");
            printf("frs[0] interface is NULL\n");
            clean_up(ERROR);
        }

        total_patch = set_copy_proc_frs(newfrs,num_patches,
                          redistr_table, max_n_patch, resid_n);

        uni_array(newwvs,total_patch,sizeof(Wave*)); 
        for(i = 0; i < total_patch; i++) (*newwvs)[i] = NULL;

        for(i = 0; i < max_n_patch; i++)
        {
            if(-1 == redistr_table[myid][i].wv_id) continue;
            if(myid == redistr_table[myid][i].pc_id)
            {
                patch_id = redistr_table[myid][i].wv_id;
                (*newwvs)[patch_id] = wvs[i];
                if(patch_id != (*newwvs)[patch_id]->patch_number)
                {
                    printf("ERROR: make_reinstore_patch_storage()\n");
                    printf("resid Wave[%d] NOT same in redistr_table\n",i);
                    clean_up(ERROR);
                }
            }
        }

        /* This is a temp check option, should be removed later. */
        nn = *resid_n;  
        for(i = 0; i < nn; i++)
        {
             if((*newfrs)[i] != frs[i] or (*newwvs)[i] != wvs[i])
             {
                 printf("ERROR: make_reinstore_patch_storage\n");
                 printf("resid. Front[%d] NOT same in redistr_table and frs\n", i);
                 printf("newfrs[%d] = %p, frs[%d] = %p, newwvs[%d] = %p, wvs[%d] = %p\n",
                        i, (*newfrs)[i], i, frs[i], i, (*newwvs)[i], i, wvs[i]);  
                 clean_up(ERROR);
             }
        }

        wave = wvs[0];   
        for(i = 0; i < total_patch; i++)
        {
            if((*newwvs)[i] != NULL)
                continue;
            (*newwvs)[i] = copy_wave(wave);
            clear_wave_pointers((*newwvs)[i]); 
            wave_of_front((*newfrs)[i]) = (*newwvs)[i]; 
        }

        DEBUG_LEAVE(make_reinstore_patch_storage) 
        return total_patch;  
}


EXPORT int reinstore_undistribute_patch(
        Wave       ***newwvs,
        Front      ***newfrs,  
        Wave       **wvs,  
        Front      **frs,
        Wv_on_pc   **redistr_table,
        int        num_patches,    /* number of patches compute in this proc */
        int        max_n_patch)
{
        int        source, numnodes, myid;
        int        dist;
        int        patch_id, i;
        Front      **tmpfrs, *front;
        Front      *basefr;
        Wave       **tmpwvs, *wave;   
        int        total_patch;
        int        nn, bal_n, mm;
        INTERFACE  *current_intfc;
        COMPONENT  dummy = -1;
        RECT_GRID  *cgr, *tgr;
        RECT_GRID  *rgr;        
        size_t      len;
        byte        *storage = NULL;
        POINTER     info;
        size_t      sizest = wvs[0]->sizest;
        int         ii, jj, L[MAXD], U[MAXD];
        int         *fly_patch_id; 
 
        DEBUG_ENTER(reinstore_undistribute_patch) 

        numnodes = pp_numnodes();
        myid = pp_mynode();

        total_patch = make_reinstore_patch_storage(&tmpwvs,&tmpfrs,
            &nn,wvs,frs,redistr_table,num_patches,max_n_patch); 

        uni_array(&fly_patch_id, total_patch, sizeof(int));  
        for(i = 0; i < total_patch; i++)
            fly_patch_id[i] = -1;  

        bal_n = nn;  
        /* Recv all fronts originally in this procs  */ 
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                patch_id = -100;
                if(-1 == redistr_table[source][i].wv_id) continue;
                if(source == redistr_table[source][i].pc_id) continue;

                if(myid == redistr_table[source][i].pc_id)
                {
                    patch_id = redistr_table[source][i].wv_id;
                    front = redistr_table[source][i].front;

                    pp_send(0,(POINTER)(&patch_id),sizeof(int),source);
                    dummy = wave_of_front(front)->patch_component;
                    send_front_misc(front,&dummy,source);
                    /* 
                    printf("redist[%d][%d] should be back to original: ",source,i);
                    printf("front %p: Proc [%d] Send to proc[%d]\n",
                            front, myid, source);
                    */               
                }
                if(source == myid)
                {
                    pp_recv(0,redistr_table[source][i].pc_id,
                       (POINTER)(&patch_id),sizeof(int));
                    
                    fly_patch_id[patch_id] = patch_id;  
                    recv_front_misc(tmpfrs[patch_id],
                        &dummy, redistr_table[source][i].pc_id);
                    wave_of_front(tmpfrs[patch_id])->patch_component = dummy;
                    nn++;
                    /*  
                    printf("redist[%d][%d] received at original: ",source,i);
                    printf("patch[%d]: Proc [%d] recv from proc[%d]\n",
                            patch_id, source, redistr_table[source][i].pc_id);
                    */    
                }
                pp_gsync();
            }
        }

        /* Set recving patches */
        for(i = 0; i < total_patch; i++)
        {
            if(fly_patch_id[i] == -1) continue; 

            tmpfrs[i]->interf->modified = YES;
            rgr = tmpfrs[i]->rect_grid;
            cgr = computational_grid(tmpfrs[i]->interf);
            copy_rect_grid(cgr,rgr);
            tgr = &topological_grid(tmpfrs[i]->interf);
            tgr->Remap.remap = rgr->Remap.remap;
            set_patch_topo_grid(rgr,tgr); 

            tmpwvs[i]->rect_grid = tmpfrs[i]->rect_grid;
            tmpwvs[i]->pd_flag = tmpfrs[i]->pd_flag;
            tmpwvs[i]->patch_number = tmpfrs[i]->patch_number;
            tmpwvs[i]->patch_level = tmpfrs[i]->patch_level;   
            if (init_hyp_solution_function(tmpwvs[i],tmpfrs[i]) != GOOD_STEP)
            {
                screen("ERROR in set_recv_wv_fr(), "
                       "init_hyp_solution_function() failed for [%d]\n",i);
                clean_up(ERROR);
            } 
        }  

        nn = bal_n;  
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                patch_id = -100;
                if(-1 == redistr_table[source][i].wv_id) continue;
                if(source == redistr_table[source][i].pc_id) continue;

                if(myid == redistr_table[source][i].pc_id)
                {
                    patch_id = redistr_table[source][i].wv_id;
                    wave = NULL;
                    for(jj = 0; jj < num_patches; jj++)
                    {
                        if(frs[jj] == redistr_table[source][i].front)
                        {
                            wave = wvs[jj];
                            break;
                        }
                    }
                    pp_send(0,(POINTER)(&patch_id),sizeof(int),source);
                    len = wave_tri_soln(wave)->tri_grid->n_reg_nodes*sizest;
                    scalar(&storage,len);  
                    for(jj = 0; jj < wave->rect_grid->dim; jj++)
                    {
                        L[jj] = -(wave->rect_grid->lbuf[jj]);
                        U[jj] = wave->rect_grid->gmax[jj] 
                              + wave->rect_grid->ubuf[jj];
                    }
                    (*wave->bundle_states)(L,U,wave,storage); 
                    pp_send_large_data(storage,len,source);       
                    free(storage);  
                    /*          
                    printf("redist[%d][%d]  should be back to original: ",source,i);
                    printf("wave %p: Proc [%d] Send to proc[%d]\n",
                             wave, myid, source);
                    */  
                }
                if(source == myid)
                {
                    pp_recv(0,redistr_table[source][i].pc_id,
                        (POINTER)(&patch_id),sizeof(int));
                    wave = tmpwvs[patch_id];  
                    len = wave_tri_soln(wave)->tri_grid->n_reg_nodes*sizest;  
                    scalar(&storage,len);

                    pp_receive_large_data(storage,len,redistr_table[source][i].pc_id); 
                    for(jj = 0; jj < wave->rect_grid->dim; jj++)
                    {
                        L[jj] = -(wave->rect_grid->lbuf[jj]);
                        U[jj] = wave->rect_grid->gmax[jj] 
                              + wave->rect_grid->ubuf[jj];
                    }
                    (*wave->unbundle_states)(L,U,wave,storage);
                    free(storage);
                    nn++;
                    /*            
                    printf("redist[%d][%d] received at original: ",source,i);
                    printf("patch[%d]: Proc [%d] recv from proc[%d]\n",
                            patch_id, source, redistr_table[source][i].pc_id);
                    */                   
                }
                pp_gsync();
            }
        }

        *newwvs = tmpwvs;
        *newfrs = tmpfrs; 

        free(fly_patch_id);  

        DEBUG_LEAVE(reinstore_undistribute_patch) 
        return GOOD_STEP;  
} 

EXPORT void patch_wave_trans(
        Wave        **wvs,
        Wave        **nwvs,
        Wv_on_pc    **redistr_table,
        int         myid,
        int         nn,
        int         max_n_patch,
        int         numnodes)
{
        size_t      len;
        byte        *storage = NULL;
        POINTER     info;
        size_t      sizest = wvs[0]->sizest;
        int         source, i, dist;
        int         patch_id;
        int         ii, jj, L[MAXD], U[MAXD];
        byte        *ps;

        DEBUG_ENTER(patch_wave_trans)

        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                patch_id = -100;
                if(-1 == redistr_table[source][i].wv_id) continue;
                if(source == redistr_table[source][i].pc_id) continue;

                if(myid == source)
                {
                    patch_id = redistr_table[source][i].wv_id;
                    pp_send(0,(POINTER)(&patch_id),sizeof(int),
                            redistr_table[source][i].pc_id);

                    len = wave_tri_soln(wvs[patch_id])->tri_grid->n_reg_nodes*sizest;
                    scalar(&storage,len);
                    dist = redistr_table[source][i].pc_id;

                    for(jj = 0; jj < wvs[patch_id]->rect_grid->dim; jj++)
                    {
                        L[jj] = -(wvs[patch_id]->rect_grid->lbuf[jj]);
                        U[jj] = wvs[patch_id]->rect_grid->gmax[jj]
                              + wvs[patch_id]->rect_grid->ubuf[jj];
                    }

                    (*wvs[patch_id]->bundle_states)(L,U,wvs[patch_id],storage);
                    pp_send_large_data(storage,len,dist);

                    free(storage);
                    /*     
                    printf("Proc[%d] sent patch[%d] state(size %7d) to Proc[%d]\n",
                         myid, patch_id, len, dist);
                    */       
                }
                if(myid == redistr_table[source][i].pc_id)
                {
                    pp_recv(0,source,(POINTER)(&patch_id),sizeof(int)) ;

                    len = wave_tri_soln(nwvs[nn])->tri_grid->n_reg_nodes*sizest;
                    scalar(&storage,len);

                    pp_receive_large_data(storage,len,source);
                    for(jj = 0; jj < nwvs[nn]->rect_grid->dim; jj++)
                    {
                        L[jj] = -(nwvs[nn]->rect_grid->lbuf[jj]);
                        U[jj] = nwvs[nn]->rect_grid->gmax[jj]
                              + nwvs[nn]->rect_grid->ubuf[jj];
                    }
                    (*nwvs[nn]->unbundle_states)(L,U,nwvs[nn],storage);
                    free(storage);
                    /*              
                    printf("Proc[%d] received patch[%d] state (size %7d) from Proc[%d]\n",
                            myid, patch_id, len, source);
                    */                     
                    nn++;
                }
                pp_gsync();
            }
        }
        DEBUG_LEAVE(patch_wave_trans)
}



EXPORT  boolean    h_scatter_patch_states(
        Wave            **wvs,
        Front           **frs,
        Overparam       *overparam, 
        Wv_on_pc        **redistr_table,
        int             max_n_patch, 
        int             *iperm,
        int             swp)
{
        PP_GRID         *pp_grid = wvs[0]->pp_grid;
        RECT_GRID       *gr;  
        int             myid;
        int             me[MAXD];
        int             i,ii,side, dim;
        int             numnodes, source;  
        int             num_patches;
        INTERFACE       *intfc;  

        Trans_wv        *stransed = NULL;  /* FOR send */  
        int             stransed_n = 0; 
        int             refine; 

        Trans_wv        *rtransed = NULL;  /* FOR receive */ 
        int             rtransed_n = 0; 
        Wave            **buff_wvs;
        Front           **buff_frs; 

        DEBUG_ENTER(h_scatter_patch_states) 
        if (wvs[0]->sizest == 0)
        {
            DEBUG_LEAVE(h_scatter_states)
            return FUNCTION_SUCCEEDED;
        }

        dim = wvs[0]->rect_grid->dim;
        numnodes = pp_numnodes();  
        myid = pp_mynode();  
        num_patches = wvs[0]->totalNumberOfPatches;
        for(i = 0; i < dim; i++)
            me[i] = redistr_table[myid][0].pc_ic[i]; 
        /* 
        if(debugging("h_scatter_patch_states"))
        {
            printf("_____________________________\n");
            printf("H_SCATTER_PATCH_STATES, patch is scattered as:\n");
            for(source = 0; source < numnodes; source++)
            {
                printf("Proc[%d]: ", source);
                for(i = 0; i < max_n_patch; i++)
                {
                    if(-1 != redistr_table[source][i].wv_id)
                        printf(" %2d %2d [%d,%d]L[%d], ", 
                            redistr_table[source][i].pc_id,
                            redistr_table[source][i].wv_id,
                            redistr_table[source][i].pc_ic[0],
                            redistr_table[source][i].pc_ic[1],
                            redistr_table[source][i].wv_level);
                    else
                       printf("                 ");
                }
                printf("\n");
            }
        }
        */  
        for (side = 0; side < 2; ++side)
        {
            pp_gsync();  
            /* 
            printf("\nproc[%d] iperm[%d] = %d send_side[%d], recv_side[%d]\n", 
                  myid, swp, iperm[swp], side, (side+1)%2);    
            */         
            stransed = pp_send_patch_interior(me,swp,iperm,side,overparam,
                     redistr_table,max_n_patch,frs,wvs,num_patches,&stransed_n);  
            pp_send_patch_interior_states(me,swp,iperm,side,
              frs,wvs,num_patches, redistr_table,stransed, stransed_n);  

            rtransed = pp_receive_patch_interior(me,swp,iperm,(side+1)%2,
                 overparam,redistr_table,frs,wvs,num_patches, &rtransed_n);

            create_recv_buffer_zone(&buff_wvs, &buff_frs, 
               rtransed, rtransed_n, overparam, wvs[0], frs[0]); 
            pp_receive_patch_interior_states(me,swp,iperm,(side+1)%2,
                 overparam,redistr_table,frs,wvs,num_patches,
                    buff_frs,buff_wvs,rtransed,rtransed_n);
        }
 
        DEBUG_LEAVE(h_scatter_patch_states) 
        return FUNCTION_SUCCEEDED;  
}

LOCAL void perform_recv_scatter_patch_states(
        Overparam    *overparam,
        Wv_on_pc     **redistr_table,
        Wave         **wvs,
        Front        **frs,
        Wave         **buff_wvs,
        Front        **buff_frs,
        Trans_wv     *transed,  
        int          *iperm,
        int          swp,
        int          side)
{
        int          refine, level; 
        int          dim, myn_patch, tran_n; 
        int          my_offset[MAXD], myid;  
        RECT_GRID    *gr;  
        int          ggmax[MAXD], ggmin[MAXD],
                     glbuf[MAXD], gubuf[MAXD];
        int          base[MAXD], bound[MAXD];
        int          L[MAXD], U[MAXD];  
        int          i, j, k; 
        int          found = NO; 

        DEBUG_ENTER(perform_recv_scatter_patch_states) 

        myid = pp_mynode();  
        tran_n = buff_wvs[0]->totalNumberOfPatches; 
        if (tran_n == 0)
        {
            printf("ERROR: perform_recv_scatter_patch_states\n"); 
            printf("tran_n = 0\n");
            clean_up(ERROR);  
        } 
        myn_patch = wvs[0]->totalNumberOfPatches; 
        for(i = 0; i < tran_n; i++)
        {
            my_offset[0] = -1; 
            my_offset[1] = -1; 
            gr = buff_wvs[i]->rect_grid;  
            level = buff_wvs[i]->patch_level;   
            for(j = 0; j < myn_patch; j++)
            {
                if(level == wvs[j]->patch_level) 
                {
                    my_offset[0] = redistr_table[myid][j].off_set[0]; 
                    my_offset[1] = redistr_table[myid][j].off_set[1]; 
                    break; 
                } 
            } 
            if(my_offset[0] == -1 or my_offset[1] == -1)
            {
                printf("In perform_scatter_patch_states\n");
                printf("ERROR: same level not found\n");
                printf("buff_wvs[%d] level[%d], tran_n[%d]\n", 
                   i, level, tran_n);
                printf("Buff_wv[%d] gmax[%d, %d], base[%d,%d], bound[%d,%d]"
                  " trans_offset[%d,%d]\n",
                   i, gr->gmax[0], gr->gmax[1], transed[i].base[0],
                    transed[i].base[1], transed[i].bound[0],
                    transed[i].bound[1], transed[i].off_set[0], 
                    transed[i].off_set[1]);    
                print_rectangular_grid(gr); 
                clean_up(ERROR);  
            } 
            /*        
            printf("Local Buff_wv[%d] gmax[%d, %d], base[%d,%d], bound[%d,%d]"
                  " trans_offset[%d,%d] mynode offset[%d,%d] \n",
                i, gr->gmax[0], gr->gmax[1], transed[i].base[0],
                transed[i].base[1], transed[i].bound[0],
                transed[i].bound[1], transed[i].off_set[0], 
                transed[i].off_set[1],my_offset[0], my_offset[1]);    
            */   
        }   

        gr = wvs[0]->rect_grid; 
        dim = gr->dim; 
        /*  
        set_receive_domain(L,U,iperm,side,swp,gr);   
        printf("me[%d] base[%d, %d], bound[%d,%d] patch[%d] l[%d] RECEIVE:"
               " L[%d, %d], U[%d, %d]\n", myid, 
              redistr_table[myid][0].base[0],
              redistr_table[myid][0].base[1],
              redistr_table[myid][0].bound[0],
              redistr_table[myid][0].bound[1], 0,
              wvs[0]->patch_level,L[0],L[1],U[0],U[1]);
        */    
        for(i = 0; i < myn_patch; i++)
        {
            refine = 1;
            level = wvs[i]->patch_level; 
            for(j = 0; j < level; j++)
                refine *= overparam->refinementRatio; 
            for(j = 0; j < dim; j++)
            {
                ggmin[j] = 0;
                ggmax[j] = refine*gr->gmax[j];
                glbuf[j] = refine*gr->lbuf[j];
                gubuf[j] = refine*gr->ubuf[j];
                base[j] = redistr_table[myid][i].base[j];
                bound[j] = redistr_table[myid][i].bound[j];
            }
            set_receive_patch_domain(L,U,iperm,side,swp,base,
                  bound,ggmax,ggmin,glbuf,gubuf,frs[i]->rect_grid);
            /*       
            printf("me[%d] base[%d, %d], bound[%d,%d] patch[%d] l[%d] RECEIVE:", 
                  myid, base[0], base[1], bound[0], bound[1], i, level);
            */      
            if((L[0]-U[0])*(L[1]-U[1]) != 0)
            {
                int         ix, ixmin, ixmax;
                int         iy, iymin, iymax; 
                int         ic[MAXD], match[MAXD], where;  
                Locstate    s_st, d_st; 
                /*      
                printf(" L[%d, %d], U[%d, %d]\n",L[0],L[1],U[0],U[1]);
                */  
                ixmin = L[0]; ixmax = U[0]; 
                iymin = L[1]; iymax = U[1]; 
                for (iy = iymin; iy < iymax; ++iy)
                {
                    ic[1] = iy;
                    for (ix = ixmin; ix < ixmax; ++ix)
                    {
                        ic[0] = ix;
                        d_st = Rect_state(ic,wvs[i]); 
                        found = match_point_in_zone(ic,base,bound,
                          level, transed, tran_n, match, &where);  
                        /* 
                        if(ic[0] == 46 and ic[1] == 15 and wvs[i]->patch_number == 10)
                        {
                            printf("level[%d] crds[%g, %g] found = %d\n",level, 
                                  Rect_coords(ic,wvs[i])[0],
                                  Rect_coords(ic,wvs[i])[1], found); 
                        }  
                        */         
                        if(found)
                        {
                            s_st = Rect_state(match,buff_wvs[where]); 
                            ft_assign(d_st,s_st,wvs[0]->sizest); 
                            /*                                
                            if(ic[0] == 46 and ic[1] == 15 and 
                               wvs[i]->patch_number == 10)
                            {
                                printf("level[%d] crds[%g, %g] found\n",level, 
                                  Rect_coords(ic,wvs[i])[0],
                                  Rect_coords(ic,wvs[i])[1]); 
                                printf("match[%d, %d] in wave[%d]\n", 
                                  match[0], match[1], where); 

                                for(int tmpi = 0; tmpi < tran_n; tmpi++)
                                {
                                    RECT_GRID *ttgr = buff_wvs[tmpi]->rect_grid; 
                                    printf("Buff_wv[%d]level[%d] gmax[%d, %d] base[%d,%d] bound[%d,%d]"
                                     " trans_offset[%d,%d] mynode offset[%d,%d] \n",
                                     tmpi, buff_wvs[tmpi]->patch_level,
                                     ttgr->gmax[0], ttgr->gmax[1], transed[tmpi].base[0],
                                     transed[tmpi].base[1], transed[tmpi].bound[0],
                                     transed[tmpi].bound[1], transed[tmpi].off_set[0], 
                                     transed[tmpi].off_set[1],my_offset[0], my_offset[1]);    
                                     if(tmpi == where)
                                     {
                                         printf("print buffer[%d], level[%d] states\n", tmpi,
                                           buff_wvs[tmpi]->patch_level);
                                         (*buff_wvs[tmpi]->show_wave_states)(buff_wvs[tmpi]); 
                                     } 
                                }  
                                printf("ft_assigned states:\n"); 
                                (*frs[0]->print_state)(s_st); 
                            }  
                            */                                      
                            /* 
                            printf("ic[%d,%d] matches buffer ic[%d,%d]\n",
                                ic[0],ic[1], match[0], match[1]);  
                            */  
                        }  
                        else
                        {
                        }  
                    }
                } 
            } 
        } 
        /*     
        fill_recv_scatter_patch_obstacle_states(overparam,
          redistr_table,wvs,frs,buff_wvs,buff_frs,transed,
          iperm,swp,side); 
        */                
        DEBUG_LEAVE(perform_scatter_patch_states) 
} 

LOCAL void fill_recv_scatter_patch_obstacle_states(
        Overparam    *overparam,
        Wv_on_pc     **redistr_table,
        Wave         **wvs,
        Front        **frs,
        Wave         **buff_wvs,
        Front        **buff_frs,
        Trans_wv     *transed,  
        int          *iperm,
        int          swp,
        int          side)
{
        int          refine, level; 
        int          dim, myn_patch, tran_n; 
        int          my_offset[MAXD], myid;  
        RECT_GRID    *gr;  
        int          ggmax[MAXD], ggmin[MAXD],
                     glbuf[MAXD], gubuf[MAXD];
        int          base[MAXD], bound[MAXD];
        int          L[MAXD], U[MAXD];  
        int          i, j, k; 
        int          found = NO; 

        DEBUG_ENTER(fill_recv_scatter_patch_obstacle_states) 

        myid = pp_mynode();  
        tran_n = buff_wvs[0]->totalNumberOfPatches; 
        if (tran_n == 0)
        {
            printf("ERROR: perform_recv_scatter_patch_states\n"); 
            printf("tran_n = 0\n");
            clean_up(ERROR);  
        } 
        myn_patch = wvs[0]->totalNumberOfPatches; 

        gr = wvs[0]->rect_grid; 
        dim = gr->dim; 
        for(i = 1; i < myn_patch; i++)
        {
            refine = 1;
            level = wvs[i]->patch_level; 
            for(j = 0; j < level; j++)
                refine *= overparam->refinementRatio; 
            for(j = 0; j < dim; j++)
            {
                ggmin[j] = 0;
                ggmax[j] = refine*gr->gmax[j];
                glbuf[j] = refine*gr->lbuf[j];
                gubuf[j] = refine*gr->ubuf[j];
                base[j] = redistr_table[myid][i].base[j];
                bound[j] = redistr_table[myid][i].bound[j];
            }
            set_receive_patch_domain(L,U,iperm,side,swp,base,
                  bound,ggmax,ggmin,glbuf,gubuf,frs[i]->rect_grid);

            if((L[0]-U[0])*(L[1]-U[1]) != 0)
            {
                int         ix, ixmin, ixmax;
                int         iy, iymin, iymax; 
                int         ic[MAXD], match[MAXD], where;  
                Locstate    s_st, d_st; 

                ixmin = L[0]; ixmax = U[0]; 
                iymin = L[1]; iymax = U[1]; 
                for (iy = iymin; iy < iymax; ++iy)
                {
                    ic[1] = iy;
                    for (ix = ixmin; ix < ixmax; ++ix)
                    {
                        ic[0] = ix;
                        d_st = Rect_state(ic,wvs[i]); 
                        found = match_point_in_zone(ic,base,bound,
                          level, transed, tran_n, match, &where);  
                        if(not found) 
                        {
                            (*wvs[i]->overture_fill_patch_amr_buffer_pt)
                             (ic,wvs[i]->patch_number,wvs,frs); 
                        }  
                    }
                } 
            } 
        } 
        DEBUG_LEAVE(fill_recv_scatter_patch_obstacle_states) 
} 

#define in_between(ic,imin,imax) ((imin)<= (ic) && (ic) < (imax))

LOCAL int match_point_in_zone(
	int        *icrds,
        int        *base,
        int        *bound,
        int        level, 
        Trans_wv   *transed, 
        int        tran_n, 
        int        *match,
        int        *where)  
{
        int        i, ic[MAXD];  
        int        dim = 2;
        RECT_GRID  *gr;  
        int        bufbase[MAXD]; 
        int        bufbound[MAXD];  
        int        found = NO;  

        for(i = 0; i < dim; i++)
            ic[i] = icrds[i]+base[i]; 

        for(i = 0; i < tran_n; i++)
        {
            if(level != transed[i].wv_level) 
                continue;   

            bufbase[0] = transed[i].base[0]; 
            bufbase[1] = transed[i].base[1];
            bufbound[0] = transed[i].bound[0]; 
            bufbound[1] = transed[i].bound[1]; 
            if(in_between(ic[0],bufbase[0],bufbound[0]) and
               in_between(ic[1],bufbase[1],bufbound[1])) 
            {
                found = YES; 
                match[0] = ic[0]-bufbase[0]; 
                match[1] = ic[1]-bufbase[1]; 
                *where = i; 
                break;  
            } 
        }   

        if(not found) return NO;  
        return YES;  
} 

LOCAL void create_recv_buffer_zone(
	Wave	***newbuff_wvs,
        Front   ***newbuff_frs,
        Trans_wv  *transed,
        int       total_n,  
        Overparam *overparam,
        Wave      *wave,
        Front     *front)
{
        RECT_GRID       *gr, *rgs, *cgr, *tgr;  
        int             i,ii,dim;
        int             refine; 
        Wave            **buff_wvs;
        Front           **buff_frs; 
        INTERFACE       *sav_intfc;
        boolean            sav_copy;
      
        gr = front->rect_grid;  
        dim = gr->dim;  

        if(0 == total_n)
            return;  
        
        uni_array(&buff_frs,total_n,sizeof(Front*));
        uni_array(&buff_wvs,total_n,sizeof(Wave*));
     
        for(i = 0; i < total_n; i++)
        {
            buff_frs[i] = deep_copy_front(front);
            buff_wvs[i] = copy_wave(wave);

            rgs = buff_wvs[i]->rect_grid = buff_frs[i]->rect_grid;
            rgs->dim = dim;
            rgs->Remap = wave->rect_grid->Remap;  
            refine = 1;
            for(ii = 0; ii < transed[i].wv_level; ii++)
                refine *= overparam->refinementRatio;
            for(ii = 0; ii < dim; ii++)
                rgs->h[ii] = gr->h[ii]/refine;  
            rgs->GL[0] = rgs->L[0] = gr->L[0]+transed[i].base[0]*rgs->h[0];                
            rgs->GL[1] = rgs->L[1] = gr->L[1]+transed[i].base[1]*rgs->h[1];                
            rgs->GU[0] = rgs->U[0] = gr->L[0]+transed[i].bound[0]*rgs->h[0];                
            rgs->GU[1] = rgs->U[1] = gr->L[1]+transed[i].bound[1]*rgs->h[1];                
            rgs->gmax[0] = transed[i].bound[0]-transed[i].base[0]; 
            rgs->gmax[1] = transed[i].bound[1]-transed[i].base[1]; 
            rgs->lbuf[0] = rgs->lbuf[1] = rgs->ubuf[0] = rgs->ubuf[1] = 0; 
            set_rect_grid(rgs->L,rgs->U,rgs->GL,rgs->GU,
              rgs->lbuf,rgs->ubuf,rgs->gmax,rgs->dim,&rgs->Remap,rgs); 

            buff_frs[i]->patch_level 
                = buff_wvs[i]->patch_level = transed[i].wv_level; 
            buff_frs[i]->patch_number = buff_wvs[i]->patch_number = i;  
            buff_frs[i]->totalNumberOfPatches = 
                 buff_wvs[i]->totalNumberOfPatches = total_n; 
            /*  
            chart_of_front(buff_frs[i]) = NULL;
            */  
            clear_wave_pointers(buff_wvs[i]);  
            wave_of_front(buff_frs[i]) = buff_wvs[i];
            /* 
            printf("I[%d] RECEIVED patch[%d] l[%d]: L[%f,%f] U[%f,%f]\n", 
              pp_mynode(), i, transed[i].wv_level, rgs->L[0], rgs->L[1], 
              rgs->U[0], rgs->U[1]); 
            */  
        }

        sav_intfc = current_interface();
        sav_copy = copy_intfc_states();

        set_size_of_intfc_state(size_of_state(front->interf));
        set_copy_intfc_states(YES);
 
        for(i = 0; i < total_n; i++)
        {
            buff_frs[i]->interf = copy_interface(front->interf);
            /* Do we need interface for the buffer zone patches ? */ 
            delete_patch_all_curves(buff_frs[i]); 
            rgs = buff_frs[i]->rect_grid;  
            cgr = computational_grid(buff_frs[i]->interf);
            copy_rect_grid(cgr,rgs);
            tgr = &topological_grid(buff_frs[i]->interf);
            tgr->Remap.remap = rgs->Remap.remap;
            set_patch_topo_grid(rgs,tgr);   

            if(init_buffer_hyp_solution_function(buff_wvs[i],buff_frs[i]) 
                 != GOOD_STEP)
            {
                screen("ERROR in create_recv_buffer_zone(), "
                       "init_buffer_hyp_solution_function() failed for [%d]\n",i);
                clean_up(ERROR);
            }

        } 
        set_current_interface(sav_intfc);
        set_copy_intfc_states(sav_copy);  

        *newbuff_wvs = buff_wvs;
        *newbuff_frs = buff_frs; 
}

/*      init_buffer_hyp_solution_function() 
*   Mainly alocate state and comp. storage 
*/
LOCAL int init_buffer_hyp_solution_function(
        Wave  *wave,
        Front *front)
{
        RECT_GRID  Dual_grid, *comp_grid = wave->rect_grid;
        RECT_GRID  *expanded_dual_grid;  
        int        i, *gmax, dim;
        int        j,n_reg_nodes; 
        int        status = ERROR_IN_STEP;
        TRI_SOLN   *soln;
        TRI_GRID   *grid; 
        register Locstate *state;
        register byte     *storage;
        size_t     sizest = front->sizest; 
	Table	   *T = table_of_interface(grid->grid_intfc);

        clear_wave_pointers(wave);
        if (wave->sizest == 0)
            return GOOD_STEP;
        set_dual_grid(&Dual_grid,comp_grid);

        scalar(&wave_tri_soln(wave),sizeof(TRI_SOLN));
        if (wave_tri_soln(wave) == NULL)
        {
            (void) printf("WARNING in init_buffer_hyp_solution_function(), "
                          "can't allocate tri_soln\n");
            return ERROR_IN_STEP;
        }
        soln = wave_tri_soln(wave);
/*     
        soln->cg_over = wave->cg_over;
        soln->patch_number = wave->patch_number;
        soln->use_overture_state = wave->use_overture_state;
        soln->overture_init_step = wave->overture_init_step;
        soln->cg_over_function = wave->cg_over_function;
        soln->patch_level = wave->patch_level;
        soln->NumberOfLevels = wave->NumberOfLevels;
*/     
        soln->Tri_grid_hooks = wave->Tri_grid_hooks;

        grid = allocate_tri_grid(&(soln->Tri_grid_hooks));
        if (grid == NULL)
            return ERROR_IN_STEP;
        soln->tri_grid = grid;   

        set_tri_grid_rect_grids(grid,&Dual_grid,front->interf); 
        expanded_dual_grid = &grid->rect_grid; 
        gmax = expanded_dual_grid->gmax;
        dim = expanded_dual_grid->dim;   

        n_reg_nodes = gmax[0] + 1;
        for (i = 1; i < dim; ++i)
            n_reg_nodes *= (gmax[i] + 1);
        grid->n_node_points = grid->n_reg_nodes = n_reg_nodes;
	uni_array(&T->components,n_reg_nodes,sizeof(COMPONENT));
	ntg->cg_comps = T->components + ntg->node_offset;
        alloc_components_array(grid,n_reg_nodes);
        alloc_states_array(grid,n_reg_nodes);
        VECTOR(grid,rect_state_storage,n_reg_nodes,sizest);
        alloc_node_points(grid,grid->n_node_points); 
 
        /* storage for regular points */ 
        state = grid->states;
        storage = grid->rect_state_storage;
        j = 0;
        switch (dim)
        {
#if defined(ONED)
        case 1:
        {
            register double      *xx_grid = expanded_dual_grid->edges[0];
            int                 xmax = gmax[0];
            int                 ix;

            for (ix = 0;  ix <= xmax;  ++ix)
            {
                Coords(grid->node_points+j)[0] = xx_grid[ix];
                ++j;

                *state++ = (Locstate) storage;
                storage += sizest;
            }
        }
                break;
#endif /* defined(ONED) */
#if defined(TWOD)
        case 2:
        {
            register double y;
            register double *xx_grid = expanded_dual_grid->edges[0];
            register double *yy_grid = expanded_dual_grid->edges[1];
            int            xmax = gmax[0],ymax = gmax[1];
            int            ix,iy;

            for (iy = 0;  iy <= ymax;  ++iy)
            {
                y = yy_grid[iy];
                for (ix = 0;  ix <= xmax;  ++ix)
                {
                    Coords(grid->node_points+j)[0] = xx_grid[ix];
                    Coords(grid->node_points+j)[1] = y;
                    ++j;

                    *state++ = (Locstate) storage;
                    storage += sizest;
                }
            }
        }
            break;
#endif /* defined(TWOD) */
#if defined(THREED)
        case 3:
        {
            register double y,z;
            register double *xx_grid = expanded_dual_grid->edges[0];
            register double *yy_grid = expanded_dual_grid->edges[1];
            register double *zz_grid = expanded_dual_grid->edges[2];
            int            xmax = gmax[0],ymax = gmax[1],zmax = gmax[2];
            int            ix,iy,iz;

            for (iz = 0;  iz <= zmax;  ++iz)
            {
                z = zz_grid[iz];
                for (iy = 0;  iy <= ymax;  ++iy)
                {
                    y = yy_grid[iy];
                    for (ix = 0;  ix <= xmax;  ++ix)
                    {
                        Coords(grid->node_points+j)[0] = xx_grid[ix];
                        Coords(grid->node_points+j)[1] = y;
                        Coords(grid->node_points+j)[2] = z;
                        ++j;

                        *state++ = (Locstate) storage;
                        storage += sizest;
                    }
                }
            }
        }
            break;
#endif /* defined(THREED) */
        }
        free_grid_lines(&grid->rect_grid); 
        return GOOD_STEP;  
}  

/*
*                       pp_receive_patch_interior():
*
*       Receives geometric information in a single 
*              buffer region (multiple patches).
*       Set the received wave region index coords to local.  
*
*/

LOCAL Trans_wv *pp_receive_patch_interior(
        int             *me,
        int             swp,
        int             *iperm,
        int             side,
        Overparam       *overparam, 
        Wv_on_pc        **redistr_table, 
        Front           **frs,
        Wave            **wvs,
        int             num_patches,
        int             *transed_n)
{
        INTERFACE     *intfc;
        PP_GRID       *pp_grid;
        RECT_GRID     *gr;
        int           L[MAXD], U[MAXD];
        int           him[MAXD];
        int           myid, src_id;
        int           i, ii, dim;
        int           recv_num;  
        Trans_wv      *transed;  
        byte          *recv_storage = NULL;
        byte          *buf;
        POINTER       info;
        int           myoff_set[MAXD], found;  
  
        DEBUG_ENTER(pp_receive_patch_interior) 
        if (rect_boundary_type(frs[0]->interf,iperm[swp],side) 
             != SUBDOMAIN_BOUNDARY)
        {
            DEBUG_LEAVE(pp_receive_patch_interior)
            *transed_n = 0; 
            return NULL;
        }
        myid = pp_mynode();
        pp_grid = frs[0]->pp_grid; 
        gr = frs[0]->rect_grid; 
        dim = gr->dim; 
        src_id = neighbor_id(him,me,iperm[swp],side,pp_grid);
        if (myid == src_id)
        {
            DEBUG_LEAVE(pp_receive_patch_interior)
            *transed_n = 0;  
            return NULL; /* Already done */
        }
        set_receive_domain(L,U,iperm,side,swp,gr);
        /* 
        printf("me[%d] base receive buffer: L[%d, %d], U[%d, %d]\n",
                myid, L[0], L[1], U[0], U[1]);         
        */  
        pp_recv(0,src_id,(POINTER)(&recv_num),sizeof(int)) ; 
        scalar(&recv_storage, recv_num*sizeof(Trans_wv)); 
        pp_recv(0,src_id,(POINTER)recv_storage,recv_num*sizeof(Trans_wv)); 
        /* 
        printf("me[%d] recv patch regions %d from him[%d]\n",
                 myid, recv_num, src_id);  
        */ 
        uni_array(&transed, recv_num, sizeof(Trans_wv)); 
        buf = recv_storage;
        for(i = 0; i < recv_num; i++)
        {
            info = (POINTER) buf;
            ft_assign(&transed[i],info,sizeof(Trans_wv));
            buf += sizeof(Trans_wv);
        }

        for(i = 0; i < recv_num; i++)
        {
            for(ii = 0; ii < num_patches; ii++)
            {
                if(redistr_table[myid][ii].wv_level == transed[i].wv_level)    
                {
                    myoff_set[0] = redistr_table[myid][ii].off_set[0];  
                    myoff_set[1] = redistr_table[myid][ii].off_set[1];  
                    break; 
                } 
            } 
            transed[i].base[0] = transed[i].off_set[0] 
                               + transed[i].base[0] - myoff_set[0]; 
            transed[i].base[1] = transed[i].off_set[1] 
                               + transed[i].base[1] - myoff_set[1]; 
            transed[i].bound[0] = transed[i].off_set[0] 
                                + transed[i].bound[0] - myoff_set[0]; 
            transed[i].bound[1] = transed[i].off_set[1] 
                                + transed[i].bound[1] - myoff_set[1]; 
            /* 
            printf("me[%d] RECV patch[%d] l[%d]: L[%d, %d] U[%d, %d] from[%d]\n",
                myid, i, transed[i].wv_level, transed[i].base[0], transed[i].base[1],
                transed[i].bound[0],transed[i].bound[1],src_id);  
            */ 
        }

        free(recv_storage);
        *transed_n = recv_num;   
        DEBUG_LEAVE(pp_receive_patch_interior) 
        return transed;   
}

/*
*               set_receive_domain():
*
*       Sets the domain limits for a buffer zone to be received.
*
*/
LOCAL   void    set_receive_domain(
        int             *L,
        int             *U,
        int             *iperm,
        int             side,
        int             swp,
        RECT_GRID       *gr)
{
        int             dim = gr->dim;
        int             *lbuf = gr->lbuf;
        int             *ubuf = gr->ubuf;
        int             *gmax = gr->gmax;
        int             j;

        DEBUG_ENTER(set_receive_domain)
        for (j = 0; j < swp; ++j)
        {
            L[iperm[j]] = -lbuf[iperm[j]];
            U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
        }
        if (side == 0)
        {
            L[iperm[swp]] = -lbuf[iperm[swp]];
            U[iperm[swp]] = 0;
        }
        else
        {
            L[iperm[swp]] = gmax[iperm[swp]];
            U[iperm[swp]] = gmax[iperm[swp]] + ubuf[iperm[swp]];
        }
        for (j = swp+1; j < dim; ++j)
        {
            L[iperm[j]] = -lbuf[iperm[j]];
            U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
        }
        if (DEBUG)
        {
            (void) printf("swp = %d, side = %d, ",swp,side);
            print_int_vector("iperm = ",iperm,dim,"\n");
            print_int_vector("L = ",L,dim,", ");
            print_int_vector("U = ",U,dim,"\n");
        }
        DEBUG_LEAVE(set_receive_domain)
}               /*end set_receive_domain*/

/*
*               set_receive_patch_domain():
*
*       Sets the domain limits for a patch buffer zone to be transmitted.
*       The returned values L[], U[]  are in the patch index coords.
*       The L[], U[] are counted in buffered pattern as patch rect_grid.
*/
LOCAL   void    set_receive_patch_domain(
        int             *L,
        int             *U,
        int             *iperm,
        int             side,
        int             swp,
        int             *base,
        int             *bound,
        int             *ggmax,
        int             *ggmin,
        int             *glbuf,
        int             *gubuf,
        RECT_GRID       *gr)
{
        int             j, dim = gr->dim;
        int             *gmax = gr->gmax;
        int             gL[MAXD], gU[MAXD];
        int             gbase[MAXD], gbound[MAXD];

        DEBUG_ENTER(set_receive_patch_domain) 
        for(j = 0; j < dim; ++j)
        {
            gbase[j] = base[j] - gr->lbuf[j];
            gbound[j] = bound[j] + gr->ubuf[j];
        }

        for (j = 0; j < swp; ++j)
        {
            gL[iperm[j]] = -glbuf[iperm[j]];
            gU[iperm[j]] = ggmax[iperm[j]] + gubuf[iperm[j]];
        }
        if (side == 0)
        {
            gL[iperm[swp]] = -glbuf[iperm[swp]];
            gU[iperm[swp]] = 0;
        }
        else
        {
            gL[iperm[swp]] = ggmax[iperm[swp]];
            gU[iperm[swp]] = ggmax[iperm[swp]] + gubuf[iperm[swp]];
        }
        for (j = swp+1; j < dim; ++j)
        {
            gL[iperm[j]] = -glbuf[iperm[j]];
            gU[iperm[j]] = ggmax[iperm[j]] + gubuf[iperm[j]];
        }
        /* find out intersection region */
        for(j = 0; j < dim; j++)
        {
            L[iperm[j]] = max(gbase[iperm[j]], gL[iperm[j]]);
            U[iperm[j]] = min(gbound[iperm[j]], gU[iperm[j]]);
            
            L[iperm[j]] -= base[iperm[j]]; 
            U[iperm[j]] -= base[iperm[j]]; 
            if(L[iperm[j]] >= U[iperm[j]])
                L[iperm[j]] = U[iperm[j]] = 0;
        }

        DEBUG_LEAVE(set_receive_patch_domain) 
}  



/*
*                       pp_send_patch_interior_states():
*
*    Sends state information in a single buffer region(multiple patches).
*
*/

LOCAL void pp_send_patch_interior_states(
        int           *me,
        int           swp,
        int           *iperm,
        int           side,
        Front         **frs,
        Wave          **wvs,
        int           num_patches,
        Wv_on_pc      **redistr_table, 
        Trans_wv      *transed, 
        int           transed_n)
{
        INTERFACE     *intfc;
        PP_GRID       *pp_grid;  
        RECT_GRID     *gr;  
        int           i, ii, wv_id; 
        int           dim, myid, dst_id;
        int           **L, **U; 
        int           him[MAXD];  
        int           base[MAXD], bound[MAXD];  
        byte          *storage = NULL;
        byte          *buf;
        size_t        len, slen;  
        
        DEBUG_ENTER(pp_send_patch_interior_states) 

        for(i = 0; i < num_patches; i++)
        {
            intfc = frs[i]->interf;
            if (rect_boundary_type(intfc,iperm[swp],side)
                   == REFLECTION_BOUNDARY)
            {
                DEBUG_LEAVE(pp_send_patch_interior_states)
                return; /* already done */  
            }
        }

        dim = frs[0]->rect_grid->dim; 
        myid = pp_mynode();    
        pp_grid = frs[0]->pp_grid;  
        dst_id = neighbor_id(him,me,iperm[swp],side,pp_grid);   

        if (rect_boundary_type(frs[0]->interf,iperm[swp],side) 
            != SUBDOMAIN_BOUNDARY)
        {
            DEBUG_LEAVE(pp_send_patch_interior_states)
            return;
        }

        if (myid == dst_id)
        {
            printf("ERROR: In pp_send_patch_interior_states()\n"); 
            printf("CASE: myid == dst_id, need to implement\n"); 
            clean_up(ERROR);  
        }

        /* set the base front send size, then
           calculate other patches */  
        bi_array(&L,transed_n,MAXD,sizeof(int)); 
        bi_array(&U,transed_n,MAXD,sizeof(int)); 

        L[0][0] = transed[0].base[0]; 
        L[0][1] = transed[0].base[1]; 
        U[0][0] = transed[0].bound[0]; 
        U[0][1] = transed[0].bound[1]; 
        for(i = 1; i < transed_n; i++)
        {
            wv_id = transed[i].wv_id;
            for(ii = 0; ii < dim; ii++)
            {
                base[ii] = redistr_table[myid][wv_id].base[ii];
                bound[ii] = redistr_table[myid][wv_id].bound[ii];
            }
            L[i][0] = transed[i].base[0]; 
            L[i][1] = transed[i].base[1]; 
            U[i][0] = transed[i].bound[0]; 
            U[i][1] = transed[i].bound[1]; 
            L[i][0] -= base[0];
            L[i][1] -= base[1];
            U[i][0] -= base[0];
            U[i][1] -= base[1];
        }  
        len = 0;  
        for(i = 0; i < transed_n; i++)
        {
            slen = 1; 
            for(ii = 0; ii < dim; ii++)
                slen *=  U[i][ii]-L[i][ii];  
            len += slen; 
        }
        scalar(&storage,frs[0]->sizest*len); 

        buf = storage; 
        for(i = 0; i < transed_n; i++)
        {
            wv_id = transed[i].wv_id;
            (*wvs[wv_id]->bundle_states)(L[i],U[i],wvs[wv_id],buf); 
            slen = 1; 
            for(ii = 0; ii < dim; ii++)
                slen *=  U[i][ii]-L[i][ii];  
            slen *= frs[wv_id]->sizest;  
            buf += slen;  
        } 

        pp_send_large_data(storage,frs[0]->sizest*len,dst_id); 

        free_these(2,L,U);  free(storage);  
        free(transed); 
        DEBUG_LEAVE(pp_send_patch_interior_states) 
}

/*
*                       pp_receive_patch_interior_states():
*
*    Receive state information in a single buffer region(multiple patches).
*
*/

LOCAL void pp_receive_patch_interior_states(
        int           *me,
        int           swp,
        int           *iperm,
        int           side,
        Overparam     *overparam, 
        Wv_on_pc      **redistr_table, 
        Front         **frs,
        Wave          **wvs,
        int           num_patches,
        Front         **buff_frs,
        Wave          **buff_wvs,
        Trans_wv      *transed, 
        int           transed_n)
{
        INTERFACE     *intfc;
        PP_GRID       *pp_grid;
        RECT_GRID     *gr;
        int           sL[MAXD], sU[MAXD];
        int           him[MAXD];
        int           myid, src_id;
        int           i, ii, dim;
        byte          *storage = NULL;
        byte          *buf;
        size_t        len, slen;

        DEBUG_ENTER(pp_receive_patch_interior_states)

        if (rect_boundary_type(frs[0]->interf,iperm[swp],side)
             != SUBDOMAIN_BOUNDARY)
        {
            DEBUG_LEAVE(pp_receive_patch_interior_states)
            return;
        }
        myid = pp_mynode();
        pp_grid = frs[0]->pp_grid;
        gr = frs[0]->rect_grid;
        dim = gr->dim;
        src_id = neighbor_id(him,me,iperm[swp],side,pp_grid);
        if (myid == src_id)
        {
            DEBUG_LEAVE(pp_receive_patch_interior_states)
            return; /* Already done */
        }
        len = 0;
        for(i = 0; i < transed_n; i++)
        {
            slen = 1;
            for(ii = 0; ii < dim; ii++)
                slen *= buff_wvs[i]->rect_grid->gmax[ii]; 
            len += slen;
        }
        scalar(&storage,frs[0]->sizest*len);
        pp_receive_large_data(storage,len*frs[0]->sizest,src_id); 
        buf = storage; 
        for(i = 0; i < transed_n; i++)
        {
            sL[0] = sL[1] = 0; 
            sU[0] = buff_wvs[i]->rect_grid->gmax[0]; 
            sU[1] = buff_wvs[i]->rect_grid->gmax[1]; 
            (*wvs[0]->unbundle_states)(sL,sU,buff_wvs[i],buf); 
            slen = 1;
            for(ii = 0; ii < dim; ii++)
                slen *=  buff_wvs[i]->rect_grid->gmax[ii];
            slen *= frs[0]->sizest;
            buf += slen; 
        }
        free(storage); 
        perform_recv_scatter_patch_states(overparam,redistr_table,
           wvs,frs,buff_wvs,buff_frs,transed,iperm,swp,side); 
        free(transed);  
        for(i = 0; i < transed_n; i++)
        {
            deep_free_front(buff_frs[i]); 
            free_wave_pointers(buff_wvs[i]);
            free_wave(buff_wvs[i]);       
        }
        if(0 != transed_n)
        {
            free(buff_frs);  
            free(buff_wvs);  
        } 
        set_current_interface(frs[0]->interf); 
        DEBUG_LEAVE(pp_receive_patch_interior_states)
}


/*
*                       pp_send_patch_interior():
*
*    Sends geometric information in a single buffer region(multiple patches).
*
*/

LOCAL Trans_wv *pp_send_patch_interior(
        int           *me,
        int           swp,
        int           *iperm,
        int           side,
        Overparam     *overparam, 
        Wv_on_pc      **redistr_table, 
        int           max_n_patch,  
        Front         **frs,
        Wave          **wvs,
        int           num_patches,
        int           *transed_n)
{
        INTERFACE     *intfc;
        PP_GRID       *pp_grid;  
        RECT_GRID     *gr;  
        int           i, ii, nn, level; 
        int           dim, myid, dst_id;
        int           Lbase[MAXD], Ubound[MAXD]; 
        int           him[MAXD];  
        int           refine,ggmax[MAXD], ggmin[MAXD], 
                      glbuf[MAXD], gubuf[MAXD];  
        int           base[MAXD], bound[MAXD];  
        int           send_num;  
        Trans_wv      *transed;   
        byte          *send_storage = NULL;
        byte          *buf;
        POINTER       info;
        int           is_reflect = NO;   
        int           send_max_levels = 0;  
        
        DEBUG_ENTER(pp_send_patch_interior) 

        for(i = 0; i < num_patches; i++)
        {
            intfc = frs[i]->interf;
            if (rect_boundary_type(intfc,iperm[swp],side)
                   == REFLECTION_BOUNDARY)
            {
                reflect_states_across_domain(iperm,side,swp,frs[i],wvs[i]);
                is_reflect = YES; 
            }
        }
        if(YES == is_reflect)
        {
            DEBUG_LEAVE(pp_send_patch_interior)
            *transed_n = 0;  
            return NULL;
        } 

        dim = frs[0]->rect_grid->dim; 
        myid = pp_mynode();    
        pp_grid = frs[0]->pp_grid;  
        dst_id = neighbor_id(him,me,iperm[swp],side,pp_grid);   

        if (rect_boundary_type(frs[0]->interf,iperm[swp],side) 
            != SUBDOMAIN_BOUNDARY)
        {
            DEBUG_LEAVE(pp_send_patch_interior)
            *transed_n = 0; 
            return NULL;
        }

        if (myid == dst_id)
        {
            printf("ERROR: In pp_send_patch_interior_states()\n"); 
            printf("CASE: myid == dst_id, need to implement\n"); 
            clean_up(ERROR);  
        }

        for(i = 0; i < max_n_patch; i++)
        {
            if(redistr_table[dst_id][i].wv_id != -1)
            {
                if(redistr_table[dst_id][i].wv_level > send_max_levels) 
                    send_max_levels = redistr_table[dst_id][i].wv_level; 
            }   
            if(redistr_table[dst_id][i].wv_id == -1)
                break; 
        }   
        /* set the base front send size, then
           calculate other patches */  
        gr = frs[0]->rect_grid;    

        set_send_domain(Lbase,Ubound,iperm,side,swp,gr);
        /*  
        printf("me[%d] SEND patch[%d] l[%d]: L[%d, %d], U[%d, %d] to[%d]\n", 
                myid, 0, 0, Lbase[0], Lbase[1], Ubound[0], Ubound[1], dst_id);  
        */  
        send_num = 1;  
        /* Cound how many patches are cut and need to be sent off */  
        for(i = 1; i < num_patches; i++)
        {
            level = wvs[i]->patch_level;  
            if(level > send_max_levels) continue;  

            refine = 1;
            for(ii = 0; ii < level; ii++)
                refine *= overparam->refinementRatio;
            for(ii = 0; ii < dim; ii++)
            {
                ggmin[ii] = 0;
                ggmax[ii] = refine*gr->gmax[ii];
                glbuf[ii] = refine*gr->lbuf[ii];
                gubuf[ii] = refine*gr->ubuf[ii];  
                base[ii] = redistr_table[myid][i].base[ii]; 
                bound[ii] = redistr_table[myid][i].bound[ii]; 
            }
            set_send_patch_domain(Lbase,Ubound,iperm,side,swp,base,
                  bound,ggmax,ggmin,glbuf,gubuf,frs[i]->rect_grid); 
            if((Lbase[0]-Ubound[0])*(Lbase[1]-Ubound[1]) != 0)
            {
                send_num++;  
                  
/*              if(level == 3)   
                {
                    printf("me[%d] SEND patch[%d] l[%d]:"
                      " Lbase[%d, %d], Ubound[%d, %d] to[%d]\n", myid,i, 
                      level,Lbase[0],Lbase[1],Ubound[0],Ubound[1],dst_id); 
                    printf("My base[%d, %d], bound[%d,%d]\n", base[0],
                      base[1], bound[0], bound[1]);  
                 } 
*/                         
            }
        }  

        set_send_domain(Lbase,Ubound,iperm,side,swp,gr);
        uni_array(&transed, send_num, sizeof(Trans_wv)); 
        transed[0].wv_id = 0;
        transed[0].wv_level = 0;
        transed[0].base[0] = Lbase[0]; 
        transed[0].base[1] = Lbase[1]; 
        transed[0].bound[0] = Ubound[0]; 
        transed[0].bound[1] = Ubound[1]; 
        transed[0].off_set[0] = redistr_table[myid][0].off_set[0]; 
        transed[0].off_set[1] = redistr_table[myid][0].off_set[1]; 
        nn = 1;  
        for(i = 1; i < num_patches; i++)
        {
            level = wvs[i]->patch_level;  
            if(level > send_max_levels) continue;  

            refine = 1;
            for(ii = 0; ii < level; ii++)
                refine *= overparam->refinementRatio;
            for(ii = 0; ii < dim; ii++)
            {
                ggmin[ii] = 0;
                ggmax[ii] = refine*gr->gmax[ii];
                glbuf[ii] = refine*gr->lbuf[ii];
                gubuf[ii] = refine*gr->ubuf[ii];  
                base[ii] = redistr_table[myid][i].base[ii]; 
                bound[ii] = redistr_table[myid][i].bound[ii]; 
            }
            set_send_patch_domain(Lbase,Ubound,iperm,side,swp,base,
                  bound,ggmax,ggmin,glbuf,gubuf,frs[i]->rect_grid); 

            if((Lbase[0]-Ubound[0])*(Lbase[1]-Ubound[1]) != 0)
            {
                transed[nn].wv_id = i;
                transed[nn].wv_level = level;
                transed[nn].base[0] = Lbase[0]; 
                transed[nn].base[1] = Lbase[1]; 
                transed[nn].bound[0] = Ubound[0]; 
                transed[nn].bound[1] = Ubound[1]; 
                transed[nn].off_set[0] = redistr_table[myid][i].off_set[0]; 
                transed[nn].off_set[1] = redistr_table[myid][i].off_set[1]; 
                nn++;  
            }
        }  

        if(nn != send_num)  
        {
            printf("ERROR: pp_send_patch_interior\n"); 
            clean_up(ERROR);  
        } 
        scalar(&send_storage, send_num*sizeof(Trans_wv)); 
        buf = send_storage;
        for(i = 0; i < send_num; i++)
        {
            info = (POINTER) buf;
            ft_assign(info,&transed[i],sizeof(Trans_wv));
            buf += sizeof(Trans_wv);
        }
        pp_send(0,(POINTER)(&send_num),sizeof(int),dst_id); 
        pp_send(0,(POINTER)send_storage,send_num*sizeof(Trans_wv),dst_id); 
        /* 
        printf("me[%d] send %d patch regions to him[%d]\n", 
                     myid, send_num, dst_id);  
        */  
        free(send_storage);  
        DEBUG_LEAVE(pp_send_patch_interior) 
        *transed_n = send_num;
        return transed; 
}



/*
*               set_send_patch_domain():
*
*       Sets the domain limits for a patch buffer zone to be transmitted.
*       The returned values Lbase[], Ubound[]  are in the subdomain index.  
*       The Lbase[], Ubound[] are counted in non-buffer pattern. 
*/

LOCAL   void    set_send_patch_domain(
        int             *L,
        int             *U,
        int             *iperm,
        int             side,
        int             swp,
        int             *base,
        int             *bound, 
        int             *ggmax,
        int             *ggmin,
        int             *glbuf,
        int             *gubuf,  
        RECT_GRID       *gr)
{
        int             j, dim = gr->dim;
        int             *gmax = gr->gmax;
        int             gL[MAXD], gU[MAXD]; 
        int             gbase[MAXD], gbound[MAXD]; 
                         /* non-buffered grid of the patch */

        DEBUG_ENTER(set_send_patch_domain)
        for(j = 0; j < dim; ++j)
        {
            gbase[j] = base[j] - gr->lbuf[j]; 
            gbound[j] = bound[j] + gr->ubuf[j]; 
        } 
        for (j = 0; j < swp; ++j)
        {
            gL[iperm[j]] = -glbuf[iperm[j]];
            gU[iperm[j]] = ggmax[iperm[j]] + gubuf[iperm[j]];
        }
        if (side == 0)
        {
            gL[iperm[swp]] =  0;
            gU[iperm[swp]] = glbuf[iperm[swp]];
        }
        else
        {
            gL[iperm[swp]] = ggmax[iperm[swp]] - gubuf[iperm[swp]];
            gU[iperm[swp]] = ggmax[iperm[swp]];
        }
        for (j = swp+1; j < dim; ++j)
        {
            gL[iperm[j]] = -glbuf[iperm[j]];
            gU[iperm[j]] = ggmax[iperm[j]] + gubuf[iperm[j]];
        }
        
        /* find out intersection region */ 
        for(j = 0; j < dim; j++)
        {
            L[iperm[j]] = max(gbase[iperm[j]], gL[iperm[j]]);
            U[iperm[j]] = min(gbound[iperm[j]], gU[iperm[j]]);
            if(L[iperm[j]] >= U[iperm[j]])
                    L[iperm[j]] = U[iperm[j]] = 0;  
        } 
        /* 
        printf("global send buffer: L[%d, %d], U[%d, %d]\n",
               gL[0], gL[1], gU[0], gU[1]); 
        printf("send buffer icrds: L[%d, %d], U[%d, %d]\n",
               L[0], L[1], U[0], U[1]); 
        */  
        DEBUG_LEAVE(set_send_patch_domain)
}  

/*
*               set_send_domain():
*
*       Sets the domain limits for a buffer zone to be transmitted.
*
*/

LOCAL   void    set_send_domain(
        int             *L,
        int             *U,
        int             *iperm,
        int             side,
        int             swp,
        RECT_GRID       *gr)
{
        int             dim = gr->dim;
        int             *lbuf = gr->lbuf;
        int             *ubuf = gr->ubuf;
        int             *gmax = gr->gmax;
        int             j;

        DEBUG_ENTER(set_send_domain)
        for (j = 0; j < swp; ++j)
        {
            L[iperm[j]] = -lbuf[iperm[j]];
            U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
        }
        if (side == 0)
        {
            L[iperm[swp]] =  0;
            U[iperm[swp]] = lbuf[iperm[swp]];
        }
        else
        {
            L[iperm[swp]] = gmax[iperm[swp]] - ubuf[iperm[swp]];
            U[iperm[swp]] = gmax[iperm[swp]];
        }
        for (j = swp+1; j < dim; ++j)
        {
            L[iperm[j]] = -lbuf[iperm[j]];
            U[iperm[j]] = gmax[iperm[j]] + ubuf[iperm[j]];
        }
        DEBUG_LEAVE(set_send_domain)
} 

#endif /* if defined(USE_OVERTURE) */  
