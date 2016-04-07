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
 * *                               doverturepatchmesh.c:
 * *
 * *       Copyright 1999 by The University at Stony Brook, All rights reserved.
 * *
 * *       Contains routines for initializing and maintaining nodes in the
 * *       automatic mesh grid level tree.
 * */

#define DEBUG_STRING    "doverturepatchmesh"
#include <driver/ddecs.h>

LOCAL   char    s[Gets_BUF_SIZE];       /*Scratch array for line input*/ 

#if defined(USE_OVERTURE) 

struct _Wv_info{
        int    wv_level;  
        int    pc_ic[MAXD];  
        int    gmax[MAXD];  /* local ggmax (domain with refinement) at this level */ 
        int    base[MAXD];  /* local grid location index   */
        int    bound[MAXD]; /* local grid location index   */
};

typedef struct _Wv_info Wv_info; 

/* declare local variables here */

/* declare local functions here */

LOCAL   CHART   *copy_chart(CHART*);
LOCAL   void    free_chart(CHART*); 
LOCAL   void    amr_main_time_step_loop(CHART**,Printplot*,int*,Wv_on_pc**,int*);
LOCAL   CHART   **perform_amr_initialization(INIT_DATA*,INIT_PHYSICS*,int*,Wv_on_pc***,int*);
LOCAL   void    set_patch_wv_fr_before_redistr(int,CHART*,Front**,Wave**,int*);  
LOCAL   int     overture_init_mesh_refinement(Overparam*,Wave*,Front*,int*);
LOCAL   int     perform_overture_init_mesh_refinement(CompositeGrid*,
                   doubleCompositeGridFunction*,CompositeGrid**,Overparam*,int); 

LOCAL   int     stop_adding_amr_grids(CompositeGrid,CompositeGrid);
LOCAL   CHART   **perform_amr_restart_initialization(INIT_DATA*,INIT_PHYSICS*,int*,
                    Wv_on_pc***,int*);
LOCAL   void    read_overture_restart_data(CHART*,Printplot*);  

LOCAL   int     find_same_comp_ic_on_base(Wave*,COMPONENT,int*); 
LOCAL   void    set_patch_components(Wave*, CompositeGrid*, int,COMPONENT*); 
LOCAL   void    set_patch_rect_grid(RECT_GRID*,RECT_GRID*,Overparam*,
                              CompositeGrid*,int,int); 
LOCAL   void    set_bdry_patch_flag(Wave*,Wave*,Overparam*); 
LOCAL   void    calculate_computing_load(Wave**,Front**,int,int*);
LOCAL   Wv_on_pc  **redistr_computing_load(Wave**,Front**,int,int*);
LOCAL   Wv_on_pc  **redistr_load_as(int**,int*,int,int); 
LOCAL   void    bcast_balance(Wv_on_pc**,int,int); 
LOCAL   void    allgather_load_info(int*,int*,int*,int,int,int); 
LOCAL   void    allgather_patch_info(Wv_on_pc**,Front**,Overparam*,int,int); 
LOCAL   int     grid_off_set_of_pc(Wv_info**,int*,int,int,int); 
LOCAL   int     perform_redistribute_load(Wv_on_pc**,Wave**,Front**,int,int,
                    Wave***,Front***); 
LOCAL   int     copy_and_alloc_wv_fr(Wave**,Front**,Wave**,Front**,
	            Wv_on_pc**,int,int,int);
LOCAL   void    patch_front_trans(Front**,Front**,Wv_on_pc**,int,int,int,int); 
LOCAL   void    remove_sent_wv_fr_after_redistr(Wv_on_pc**,int,int,Front**,Wave**);  

LOCAL   void    set_recv_wv_fr(Wave**,Front**,int,int);  

LOCAL   int     amr_time_step(CHART**,double,double*,int,double,int,Wv_on_pc**,int); 
LOCAL   int     overture_hyp_amr(double,double*,Overparam*,Wave**,Front**,
                         int,Wv_on_pc**,int);
LOCAL   int     overture_parab_amr_driver(double,double*,Overparam*,
                         Wave**,Front**,int,Wv_on_pc**,int,
                         void(*)(double,Front*,Wave*,Wave*)); 
LOCAL   CHART   **amr_regridding(CHART**,Printplot*,Wv_on_pc**,int*,int*,
                   Wv_on_pc***,int*,int*); 
LOCAL   void    d_overture_injection_in_regridding(Wave**,Front**,
                  Wave**,Front**,Wv_on_pc**,Overparam*,int); 
LOCAL   CHART   **amr_overture_repatch(CHART**,Printplot*,Wv_on_pc**,int*,int*); 
LOCAL   int     overture_perform_regrid(CHART*,Printplot*,CompositeGrid*,
                   doubleCompositeGridFunction*,CompositeGrid**,
                   doubleCompositeGridFunction**,int*);
LOCAL   void    copy_gridfunction_from_uu1(const doubleCompositeGridFunction&,
                  const CompositeGrid&,int,doubleCompositeGridFunction&);   
LOCAL   void    d_overture_injection_after_repatch(CHART**,Wave**,Front**); 
LOCAL   void    update_error(doubleCompositeGridFunction*,const CompositeGrid&,Overparam*);

LOCAL   void    tmp_save_overture_show_files(Printplot*,Wave**,Front**,int,int);
LOCAL   size_t  bundle_send_info(byte*,POINTER,int,size_t); 
LOCAL   size_t  unbundle_recv_info(byte*,POINTER,int,size_t); 

LOCAL   double   find_single_proc_time_step(Grid*,Wave*,Front*,double*); 
LOCAL   int     create_overture_dump_state(CHART**,Printplot*,boolean,
                     int,int,Wv_on_pc**,int); 
LOCAL   void    perform_create_overture_state_function(CHART**,
                  int,Wv_on_pc**,int); 
LOCAL   void    show_grids_info(CompositeGrid&);

EXPORT void read_overture_option(
        INIT_DATA       *init,
        INIT_PHYSICS    *ip)
{
        CHART           *root = ip->root;
        Printplot       *prt = ip->prt;
        Front           *front = root->front;
        Wave            *wave = root->wave;
        Overparam       *overparam;

        debug_print("chart","Entering read_overture_option() \n");

        screen("\nType 'y' to have local mesh refinement: ");
        (void) Gets(s);
        if (s[0] != 'y' and s[0] != 'Y')
        {
                chart_of_front(front) = NULL;
                debug_print("chart","Left init_amr()\n");
                root->overparam = NULL;
                return;
        }

        scalar(&(root->overparam), sizeof(Overparam));

        /* initialize Overparam */
        overparam = root->overparam;
        overparam->baseLevel = 0;
        overparam->numberOfSmooths = 2; 
        overparam->efficiency = 0.5; 
        overparam->errorThreshold = 0.1; 
        overparam->extra_buf = 0; 

        screen("\n-- Begin prompting for AMR mesh refinement --\n");
        screen("\nEnter using Overture Storage (dflt n): ");
        (void) Gets(s);
        if (s[0] == '\0' || (s[0] != 'y' && s[0] != 'Y'))
        {
            root->wave->use_overture_state = 
                   root->use_overture_state = NO;
        }
        else
        {
            root->wave->use_overture_state = 
                    root->use_overture_state = YES;
        } 
        screen("\nEnter time step for regridding: ");
        (void) Scanf("%d\n",&overparam->period);

        screen("\nEnter base level number(deflt %d): ",overparam->baseLevel);
        (void) Gets(s);
        if(s[0] != '\0')
        {
            if(sscanf(s,"%d\n",&overparam->baseLevel) != 1)
            {
                screen("ERROR in read_overture_option(), "
                       "couldn't scan baseLevel\n");
                clean_up(ERROR);
            }  
            /* (void) Scanf("%d\n",&overparam->baseLevel); */
        }

        screen("\nEnter levels of refinement: ");
        (void) Scanf("%d\n",&overparam->numberOfRefinementLevels);

        wave->NumberOfLevels = front->NumberOfLevels = 
                 overparam->numberOfRefinementLevels; 

        screen("\nEnter refinement ratio: ");
        (void) Scanf("%d\n",&overparam->refinementRatio);

        screen("\nEnter numberOfSmooths for error(deflt %d): ",
                overparam->numberOfSmooths);
        (void) Gets(s);
        if(s[0] != '\0')
        {
            if(sscanf(s,"%d\n",&overparam->numberOfSmooths) != 1)
            {
                screen("ERROR in read_overture_option(), "
                       "couldn't scan numberOfSmooths\n");
                clean_up(ERROR);
            } 
            /* (void) Scanf("%f\n",&overparam->numberOfSmooths); */
        }

        screen("\nEnter refinement efficiency(dflt %g): ",
                overparam->efficiency);
        (void) Gets(s);
        if(s[0] != '\0')
        {
            if(sscan_float(s,&overparam->efficiency) != 1)
            {
                screen("ERROR in read_overture_option(), "
                       "couldn't scan efficiency\n");
                clean_up(ERROR);
            } 
            /* (void) Scanf("%f\n",&overparam->efficiency); */
        } 

        screen("\nEnter refinement error tolerance(deflt %g): ",
                  overparam->errorThreshold);
        (void) Gets(s);
        if(s[0] != '\0')
        {
            if(sscan_float(s,&overparam->errorThreshold) != 1)
            {
                screen("ERROR in read_overture_option(), "
                       "couldn't scan errorThreshold\n");
                clean_up(ERROR);
            } 
            /* (void) Scanf("%f\n",&overparam->errorThreshold); */
        }

        screen("\nEnter extra buffer size if needed (deflt %g): ",
                  overparam->extra_buf);
        (void) Gets(s);
        if(s[0] != '\0')
        {
            if(sscanf(s,"%d\n",&overparam->extra_buf) != 1)
            {
                screen("ERROR in read_overture_option(), "
                       "couldn't scan errorThreshold\n");
                clean_up(ERROR);
            } 
        }

        screen("\n-- End prompting for mesh refinement --\n");

        if (got_intfc_from_file(init) == YES)
        {
            root->wave->overture_init_step = NO;
        }
        else
        {
            root->wave->overture_init_step = YES;
        }

        debug_print("chart","Leaving read_overture_option() \n");
}

/*
*                               overture_init_amr():
*
*       Main initializer, called by dinit.c
*/
EXPORT void overture_init_amr(
        CHART           *root)
{
        Front           *front = root->front;
        Wave            *wave = root->wave;
        RECT_GRID       *rg = front->rect_grid;
        CompositeGrid   *cg;
        char            s[Gets_BUF_SIZE];

        debug_print("chart","Entering overture_init_amr() \n");

        static SquareMapping square(rg->L[0], rg->U[0], rg->L[1], rg->U[1]);
        square.setGridDimensions(axis1,rg->gmax[0]+1);
        square.setGridDimensions(axis2,rg->gmax[1]+1);
        static MappedGrid mg(square);

        for(int axis = 0; axis < mg.numberOfDimensions(); axis++)
        {
            for(int side = Start; side <= End; side++)
            {
                mg.setBoundaryCondition(side,axis,1);
                mg.setNumberOfGhostPoints(side,axis,4);
            }
        }
        mg.update();

        cg = new CompositeGrid(2,1);
        (*cg)[0].reference(mg);
        cg->updateReferences();
        cg->update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEmask);
        front->cg_over = wave->cg_over = root->old_cg_over = (POINTER)cg;

        front->patch_number = wave->patch_number = 0;

        debug_print("chart","Left overture_init_amr()\n");
}               /*end init_amr*/

EXPORT int d_overture_amr_main(
        INIT_DATA       *init,  
        INIT_PHYSICS    *ip)	
{
	CHART		**roots;
	int             patch_number; 
        Wv_on_pc        **redistr_table;
        int             max_n_patch; 

        if(got_intfc_from_file(init))
        {
            roots = perform_amr_restart_initialization(init,ip,&patch_number,  
                     &redistr_table,&max_n_patch);
        }
        else
        {
            roots = perform_amr_initialization(init,ip, &patch_number,
                         &redistr_table,&max_n_patch);
        }
        amr_main_time_step_loop(roots,ip->prt,&patch_number,
                     redistr_table,&max_n_patch);

	return YES;  
}  

LOCAL   void amr_main_time_step_loop(
        CHART           **roots,
        Printplot       *prt,
        int             *num_patchs,
        Wv_on_pc        **redistr_table,
        int             *max_n_patch)
{
        TRI_GRID        *ntg;
        double           dt, dt_frac;
        int             i, status;
        static int      initial_step;
        Wv_on_pc        **new_redistr_table; 
        int             new_max_n_patch; 
        int             new_num_patches; 
        int             over_dump; 
#if defined(TIMING)
        double           start_cpu_time, end_cpu_time;
        double           start_wall_time, end_wall_time;
#endif /* defined(TIMING) */

        /* Time Loop */
        DEBUG_ENTER(amr_main_time_step_loop)

        start_clock("ALL_TIMESTEPS");
        for ( ; ; )
        {
            if (debugging("amr_main_time_step_loop"))
            {
                for(i = 0; i < *num_patchs; i++)
                    printf("amr_main_time_step_loop"
                       " patch %d dt = %g\n", i, roots[i]->grid->dt);
            }

            for(i = 0; i < *num_patchs; i++)
            {
                roots[i]->wave->overture_init_step = NO;
                if (debugging("amr_main_time_step_loop"))
                {
                    (void)printf("amr_main_time_step_loop()"
                            " computing in patch = %d\n",i);
                    printf("wvs[%d]->cg_over_function=%d\n",
                                i,roots[i]->wave->cg_over_function);
                    print_rectangular_grid(roots[i]->wave->rect_grid);
                    if (roots[i]->wave->show_tri_soln)
                       (roots[i]->wave->show_tri_soln)(roots[i]->front,roots[i]->wave);
                    (roots[i]->wave->show_wave_states)(roots[i]->wave);
                }
                initial_step =  roots[i]->grid->step;

#if defined(TIMING)
                start_cpu_time = cpu_seconds();
                start_wall_time = real_time();
#endif /* defined(TIMING) */
            }

            start_clock("TIMESTEP");
            status = amr_time_step(roots,roots[0]->grid->dt,&dt_frac,
                        roots[0]->grid->step, roots[0]->grid->time,*num_patchs,
                        redistr_table,*max_n_patch);
#if defined(TIMING)
            end_cpu_time = cpu_seconds();
            end_wall_time = real_time();
            (void) printf("# Step %d took %g CPU seconds, %g wall seconds.\n",
                    grid->step,end_cpu_time - start_cpu_time,
                    end_wall_time - start_wall_time);
#endif /* defined(TIMING) */

            switch (status)
            {
            case GOOD_STEP:
                break;

            case MODIFY_TIME_STEP:
                screen("amr_time_step returned MODIFY_TIME_STEP\n");
                screen("amr_time_step returned MODIFY_TIME_STEP\n");
                screen("debug exit(0)\n");
                exit(0);
                stop_clock("TIMESTEP");
                continue;

            case REPEAT_TIME_STEP:
                screen("amr_time_step returned REPEAT_TIME_STEP\n");
                delete_untracked_hyper_surfaces(roots[0]->front,roots[0]->wave);
                screen("amr_time_step returned REPEAT_TIME_STEP\n");
                screen("debug exit(0)\n");
                exit(0);
                stop_clock("TIMESTEP");
                continue;

            case ERROR_IN_STEP:
                default:
                stop_clock("TIMESTEP");
                screen("FT_ERROR: in amr_main_time_step_loop(), time_step() failed\n");
                clean_up(ERROR);
            }

            /* Reset modified time step counter after good time step */
        

            for(i = 0; i < *num_patchs; i++)
            {
                roots[i]->grid->num_mts = 0;
                roots[i]->grid->time += roots[i]->grid->dt;
                roots[i]->front->time = roots[i]->wave->time = roots[i]->grid->time;
            }

            for(i = 0; i < *num_patchs; i++)
            {
                if (roots[i]->grid->step != initial_step)
                {
                    if (roots[i]->ellip)
                    {
                        print_storage("before ELLIP","EL_storage");
                        start_clock("ELLIP");     /* Elliptic Step */
                        if (debugging("time_step"))
                        {
                           (void) printf("\nBefore ellip: interf %lu\n",
                                         interface_number(roots[i]->front->interf));
                           print_interface(roots[i]->front->interf);
                        }
                        (*roots[i]->ellip)(NULL,roots[i],NULL);
                        if (debugging("time_step"))
                        {
                            (void) printf("\nAfter ellip: interf %lu\n",
                                          interface_number(roots[i]->front->interf));
                            print_interface(roots[i]->front->interf);
                        }
                        stop_clock("ELLIP");
                        print_storage("after ELLIP","EL_storage");
                    }
                }
            }


            if (prt->printout != NULL)
            {
                if (stop_run(roots[0]->grid,roots[0]->front,roots[0]->wave))
                {
                    start_clock("PRINTOUT");
                    if(NO == roots[0]->wave->use_overture_state)
                        over_dump = create_overture_dump_state(roots,prt,YES,0,
                          *num_patchs,redistr_table,*max_n_patch);  

                    (*prt->printout)(roots[0],prt,YES,0);

                    if(NO == roots[0]->wave->use_overture_state)
                    {
                        if(over_dump == YES)
                        {
                            ((doubleCompositeGridFunction*) 
                              roots[0]->wave->cg_over_function)->destroy(); 
                            delete ((doubleCompositeGridFunction*)roots[0]->wave->cg_over_function); 
                            roots[0]->wave->cg_over_function = NULL; 
                        }
                    } 
                    stop_clock("PRINTOUT");
                    break;
                }
                else
                {
                    start_clock("PRINTOUT");
                    if(NO == roots[0]->wave->use_overture_state)
                        over_dump = create_overture_dump_state(roots,prt,NO,0,
                          *num_patchs,redistr_table,*max_n_patch);  
                     
                    (*prt->printout)(roots[0],prt,NO,0);

                    if(YES == over_dump)
                    {
                        ((doubleCompositeGridFunction*) 
                          roots[0]->wave->cg_over_function)->destroy(); 
                        delete ((doubleCompositeGridFunction*)
                          roots[0]->wave->cg_over_function); 
                        roots[0]->wave->cg_over_function = NULL; 
                    } 
                    stop_clock("PRINTOUT");
                }
            }
            else if (stop_run(roots[0]->grid,roots[0]->front,roots[0]->wave))
                break;
            dt = find_amr_time_step(roots,*num_patchs,prt,NO);  

            for(i = 0; i < *num_patchs; i++)
            {
                roots[i]->grid->dt_last = roots[i]->grid->dt;
                roots[i]->grid->step++;
                roots[i]->grid->dt = dt;
                roots[i]->front->step = roots[i]->grid->step;  
            }

            stop_clock("TIMESTEP");
            print_storage("after TIMESTEP","TIME_storage");
            if (debugging("num_intfcs"))
            {
                struct Table *T;
                int n;

                for (n = 0, T = interface_table_list(); T != NULL; T = T->next)
                {
                    n++;
                    (void) printf("Interface of table %lu = %lu\n",
                          table_number(T),interface_number(T->interface));
                }
                (void) printf("%d interfaces present at end of time step\n",n);
            }
            if(roots[0]->overparam->numberOfRefinementLevels != 0)
            {
                if(roots[0]->overparam->period != 0) 
                {
                    if( roots[0]->grid->step % 
                        roots[0]->overparam->period == 0 )
                    {
                        roots = amr_regridding(roots,prt,redistr_table,
                                 max_n_patch,num_patchs, &new_redistr_table,
                                 &new_max_n_patch, &new_num_patches);
                        redistr_table = new_redistr_table;
                        *max_n_patch = new_max_n_patch;
                        *num_patchs = new_num_patches; 
                    } 
                }
                else
                {
                    printf("ERROR: amr_main_time_step_loop\n");
                    printf("numberOfRefinementLevels = %d, regridding-period = %d\n",
                            roots[0]->overparam->numberOfRefinementLevels,
                            roots[0]->overparam->period);
                    clean_up(ERROR);  
                } 
            }
            else
            {
                if(roots[0]->overparam->period != 0) 
                {
                    if( roots[0]->grid->step %
                        roots[0]->overparam->period == 0 )
                    { 
                        Front **frs;
                        uni_array(&frs,*num_patchs,sizeof(Front*));
                        for (int tmpi = 0; tmpi < *num_patchs; tmpi++)
                            frs[tmpi] = roots[tmpi]->front;
                        geomview_amr_fronts_plot2d(basename(prt->outfile),frs,
                          *num_patchs,redistr_table,*max_n_patch);
                        free(frs);  
                    }
                }
            } 
        }
        stop_clock("ALL_TIMESTEPS");

        DEBUG_LEAVE(amr_main_time_step_loop)
}  

LOCAL int create_overture_dump_state(
	CHART           **chart,
        Printplot       *prt,
        boolean            about_to_stop, 
        int             error,
        int             num_patchs,
        Wv_on_pc        **redistr_table,
        int             max_n_patch)
{
        Grid            *grid = chart[0]->grid;
        FILE            *file = stdout;
        char            fname[1024];
        OUTPUT_DATA     **mod = prt->main_output_data;
        boolean            extra_print = NO;
        boolean            print_time = NO;
        boolean            prt_fr_or_sts;
        int             i;

        interactive_printing_control(prt,grid,&extra_print,&about_to_stop);  

        if (extra_print == YES)
        {
            print_time = YES;
            prt_fr_or_sts = YES;
        }
        else if (about_to_stop == YES)
        {
            print_time = YES;
            prt_fr_or_sts = YES;
        }
        else
        {
            for (i = 0; i < NUM_OUTPUT_FORMATS; ++i)
            {
                if (is_ts_for_output_format(i,grid,prt))
                {
                    print_time = YES;
                    break;
                }
            }
        }

        if(YES == print_time)
        {
            perform_create_overture_state_function(chart,
                 num_patchs,redistr_table,max_n_patch);   
            return YES; 
        }  
        return NO; 
}

LOCAL void perform_create_overture_state_function(
	CHART       **chart,
        int         num_patches,
        Wv_on_pc    **redistr_table,
        int         max_n_patch)
{
        doubleCompositeGridFunction     *cg_over_function;
        CompositeGrid                   *cg_over;
        size_t                          sizest = chart[0]->wave->sizest;
        Range                           all;
        int                             i, total; 
        int                             patch_id, myid;
        RECT_GRID                       *gr; 
        int                             *gmax;   
        int                             iy, ix;   
        Index                           I1,I2,I3;  
        int                             smax[MAXD], smin[MAXD];  
        int                             base[MAXD], bound[MAXD]; 
        int                             ic[MAXD], resid_n;  
        Locstate                        st;  
        Wave                            **tmpwvs, **wvs;
        Front                           **tmpfrs, **frs;

        DEBUG_ENTER(perform_create_overture_state_function)

        uni_array(&frs,num_patches,sizeof(Front*));
        uni_array(&wvs,num_patches,sizeof(Wave*));
        for (i = 0; i < num_patches; i++)
        {
            frs[i] = chart[i]->front;
            wvs[i] = chart[i]->wave;
        }

         
        reinstore_mini_undistribute_patch(&tmpwvs,&tmpfrs,
          &resid_n,wvs,frs,redistr_table,num_patches,max_n_patch,NO); 
        free_these(2,wvs,frs);  
        
        total = chart[0]->wave->totalNumberOfPatches;  
        cg_over = (CompositeGrid*) chart[0]->wave->cg_over;
        cg_over_function = new
             doubleCompositeGridFunction(*cg_over,
               sizest/sizeof(double)+1,all,all,all); 
       
        for(i = 0; i < total; i++)
        {
            gr = tmpwvs[i]->rect_grid; 
            gmax = tmpwvs[i]->rect_grid->gmax;  
            for(int dim = 0; dim < gr->dim; dim++)
            {
                smin[dim] = -gr->lbuf[dim]; 
                smax[dim] = gr->gmax[dim]+gr->ubuf[dim];  
            }
            getIndex((*cg_over)[i].indexRange(),I1,I2,I3); 
            base[0]  = I1.getBase();
            base[1] = I2.getBase();
            bound[0] = I1.getBound();
            bound[1] = I2.getBound();
            for(iy = smin[1]; iy < smax[1]; iy++)
            {
                for(ix = smin[0]; ix < smax[0]; ix++)
                {
                    ic[0] = ix; ic[1] = iy; 
                    st = Rect_state(ic,tmpwvs[i]);   
                    ic[0] += base[0]; ic[1] += base[1];   
                    (*tmpwvs[i]->ft_to_overture_st)(st,
                      (POINTER)cg_over_function,i,ic);  
                }  
            } 
        } 

        set_current_interface(chart[0]->front->interf); 

        myid = pp_mynode();
        for (i = 0; i < max_n_patch; i++)
        {
            if(-1 == redistr_table[myid][i].wv_id) continue;
            if(myid == redistr_table[myid][i].pc_id) continue;

            patch_id = redistr_table[myid][i].wv_id;
            deep_free_front(tmpfrs[patch_id]);
            Set_free(wave_tri_soln(tmpwvs[patch_id])->tri_grid,components); 
            Set_free(wave_tri_soln(tmpwvs[patch_id])->tri_grid,states);
            Set_free(wave_tri_soln(tmpwvs[patch_id])->tri_grid,rect_state_storage);
            Set_free(wave_tri_soln(tmpwvs[patch_id])->tri_grid,node_points);  
            free_grid_lines(&(wave_tri_soln(tmpwvs[patch_id])->tri_grid->rect_grid)); 
            free(wave_tri_soln(tmpwvs[patch_id])->tri_grid);
            free(wave_tri_soln(tmpwvs[patch_id])); 
            free_wave(tmpwvs[patch_id]);
        }
        free_these(2,tmpwvs,tmpfrs); 
        chart[0]->wave->cg_over_function = (POINTER)cg_over_function;  

        DEBUG_LEAVE(perform_create_overture_state_function)
}

EXPORT void overture_dump_st_printout(
        CHART           *chart,
        Printplot       *prt)
{
        CompositeGrid *cg_over = (CompositeGrid *)chart->wave->cg_over;
        doubleCompositeGridFunction *cg_over_function; 
        int             i;
        char            outdir[256];
        int             myid, numnodes;
        const char      *nstep, *nd;
        char            outname[256], buf[256];
        char            *dname = basename(prt->outfile);  
        HDF_DataBase    db; 

        DEBUG_ENTER(overture_dump_st_printout) 

        cg_over_function = 
              (doubleCompositeGridFunction *)
                     chart->wave->cg_over_function;

        myid = pp_mynode(); numnodes = pp_numnodes();
        sprintf(outdir,"%s/%s",dname,"over_dump_st");  

        nd = right_flush(pp_mynode(),4);  
        sprintf(buf,"%s/%s_%s.%s",outdir,"over","dump_st",nd);  
        nstep = right_flush(chart->grid->step,7);
        sprintf(outname,"%s.ts%s",buf,nstep);

        if (create_directory(dname,YES) == FUNCTION_FAILED)
        {
            (void) printf("WARNING in overture_dump_st_printout(), directory "
                         "%s doesn't exist and can't be created\n",dname);
            return;
        }  
        if (create_directory(outdir,YES) == FUNCTION_FAILED)
        {
            (void) printf("WARNING in overture_dump_st_printout(), directory "
                         "%s doesn't exist and can't be created\n",outdir);
            return;
        }

        db.mount(outname,"I");
        cg_over->put(db,"My Grid");            /* save a grid */
        cg_over_function->put(db,"My Solution");         /* save a grid function */
        db.unmount();

        DEBUG_LEAVE(overture_dump_st_printout) 
} 

LOCAL   int     amr_time_step(
        CHART           **chart,
        double           dt,
        double           *dt_frac,
        int             step,
        double           time,
        int             num_patches,
        Wv_on_pc        **redistr_table,
        int             max_n_patch)
{
        int             status;
        Front           **frs;
        Wave            **wvs;
        Overparam       *overparam; 
        int             i;

        DEBUG_ENTER(amr_time_step)
        set_current_chart(chart[0]);
        uni_array(&frs,num_patches,sizeof(Front*));
        uni_array(&wvs,num_patches,sizeof(Wave*));

        for (i = 0; i < num_patches; i++)
        {
            frs[i] = chart[i]->front;
            wvs[i] = chart[i]->wave;
        }
        overparam = chart[0]->overparam;   
        start_clock("HYPERBOLIC");
        status = overture_hyp_amr(dt,dt_frac,overparam,
             wvs,frs,num_patches,redistr_table,max_n_patch);
        stop_clock("HYPERBOLIC");
        print_storage("after HYPERBOLIC","HYP_storage");
        if (status != GOOD_STEP)
        {
            free_these(2,frs,wvs);
            print_time_step_status("WARNING in amr_time_step(),  "
                      "hyp step failed, status = ",status,"\n");
            clean_up(ERROR); 
            DEBUG_LEAVE(amr_time_step)
            return status;
        }

        if(chart[0]->parab)
        {
            status = overture_parab_amr_driver(dt,dt_frac,overparam,wvs,frs,
                 num_patches,redistr_table,max_n_patch,chart[0]->parab_npt);
            if (status != GOOD_STEP) 
            {
                free_these(2,frs,wvs);
                print_time_step_status("WARNING in amr_time_step(),  "
                      "parab step failed, status = ",status,"\n");
                clean_up(ERROR); 
                DEBUG_LEAVE(amr_time_step)
                return status;
            }
        } 
       
        free_these(2,frs,wvs);
        DEBUG_LEAVE(amr_time_step)
        return status;
}               /*end amr_time_step*/

LOCAL int overture_parab_amr_driver(
        double           dt,
        double           *dt_frac,
        Overparam       *overparam,
        Wave            **wvs,
        Front           **frs,
        int             num_patches,
        Wv_on_pc        **redistr_table,
        int             max_n_patch,
        void            (*parab_npt)(double,Front*,Wave*,Wave*))
{
        Wave            **newwaves;
        int             i, k, dim = frs[0]->interf->dim;
        int             *iperm; /* permutation of {0,...,dim-1} */
        Wave            **undistr_wvs;
        Front           **undistr_frs;
        int             resid_n;
        int             numnodes, myid;
        static char     warn[] = "WARNING in overture_parab_amr_driver()";
        static char     err[] = "ERROR in overture_parab_amr_driver()";

        DEBUG_ENTER(overture_parab_amr_driver)

        if(debugging("parab_amr_driver"))
            printf("Entered overture_parab_amr_driver\n"); 

        myid = pp_mynode();
        numnodes = pp_numnodes();

        /* parabolic solver */
        iperm = set_iperm(frs[0]->step,dim);

        /* Initialize Intermediate Storage for States */
        uni_array(&newwaves,num_patches,sizeof(Wave*));
        for (i = 0; i < num_patches; i++)
        {
            newwaves[i] = copy_wave(wvs[i]);  
            assign_wave_parameters(newwaves[i],wvs[i]);
            if(!copy_hyp_solution_function(wvs[i],newwaves[i]))
            {
                screen("%s, copy_hyp_solution_function() failed on wave[%d]\n",
                          err, i);
                clean_up(ERROR);  
                DEBUG_LEAVE(overture_parab_amr_driver)
                return ERROR_IN_STEP;
            }
        }

        start_clock("pabar_solver");   
        for (i = 0; i < num_patches; i++)
            (*parab_npt)(dt,frs[i],wvs[i],newwaves[i]);

        start_clock("scatter_states");
        reinstore_undistribute_patch(&undistr_wvs,&undistr_frs,newwaves,
             frs,redistr_table,num_patches,max_n_patch);
        /* 
         * 051603: add interpolation_fine_to_coarse()
         */
        (*newwaves[0]->overture_undistribute_interpolation_fine_to_coarse)(
            undistr_wvs,undistr_frs);
        for(k = 0; k < dim; k++)
        {
            h_scatter_patch_states(undistr_wvs, undistr_frs,
                overparam, redistr_table, max_n_patch, iperm, k);
        }  
        (*newwaves[0]->scatter_patch_states)(overparam,
             redistr_table, undistr_wvs, undistr_frs,iperm);  
        resid_n = 0;
        for(int tmpi = 0; tmpi < max_n_patch; tmpi++)
        {
            if(-1 == redistr_table[pp_mynode()][tmpi].wv_id) continue;
            if(pp_mynode() == redistr_table[pp_mynode()][tmpi].pc_id)
                resid_n++;
        }
        patch_wave_trans(undistr_wvs,newwaves,redistr_table,
            myid, resid_n, max_n_patch, numnodes);
        free_copy_patches(undistr_wvs,undistr_frs,redistr_table,max_n_patch);
        free_these(2,undistr_wvs,undistr_frs);
        stop_clock("scatter_states");   
        
        for (i = 0; i < num_patches; i++)
        {
            assign_copy_wave_pointers(wvs[i],newwaves[i]);
            free_wave(newwaves[i]);
        } 
        free(newwaves); 

        stop_clock("parab_solver");   
        if(debugging("overture_parab_amr_driver"))
            printf("Left overture_parab_amr_driver\n"); 
        DEBUG_LEAVE(overture_parab_amr_driver)
        return GOOD_STEP;  
}


LOCAL int overture_hyp_amr(
        double           dt,
        double           *dt_frac,
        Overparam       *overparam,
        Wave            **wvs,
        Front           **frs,
        int             num_patches,
        Wv_on_pc        **redistr_table,
        int             max_n_patch)
{
        CompositeGrid                   *cg_over;
        doubleCompositeGridFunction     *new_over_function, 
                                        *old_over_function;
        doubleCompositeGridFunction     *swap_over_function;
        int                             timesteps;
        size_t                          sizest = wvs[0]->sizest;
        int                             i, k;
        int                             status = GOOD_STEP;
        Front                           **newfrs;
        INTERFACE                       *tmp_intfc, *sav_intfc;
        int                             *iperm;
        int                             dim = frs[0]->interf->dim;
        int                             step = frs[0]->step;  
        boolean                            sav_copy;
        Wave                            **undistr_wvs;
        Front                           **undistr_frs;  
        int                             resid_n; 
        int                             myid, numnodes;  

        DEBUG_ENTER(overture_hyp_amr)
        myid = pp_mynode(); 
        numnodes = pp_numnodes(); 

        cg_over = (CompositeGrid*)wvs[0]->cg_over;
        if(YES == wvs[0]->use_overture_state)
        {
            /* make a new cg_over_function for wave update */
            int    num_float;  
            Range  all;
            num_float = sizest/sizeof(double)+1;
            new_over_function = new doubleCompositeGridFunction();
            new_over_function->updateToMatchGrid((*cg_over),num_float,all,all,all);
            swap_over_function = new doubleCompositeGridFunction();
            swap_over_function->updateToMatchGrid((*cg_over),num_float,all,all,all);
            old_over_function = (doubleCompositeGridFunction*)wvs[0]->cg_over_function;
        }  

        status = pseudo_advance_amr_fronts2d(&newfrs, dt, dt_frac, overparam,
                      wvs, frs, num_patches, redistr_table, max_n_patch);

        if (status != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
              " failed : pseudo_advance_amr_fronts2d()\n");
            clean_up(ERROR);
        }

        /* Continue for interior state sweeps */
        status = hyp_amr_split_driver(dt,dt_frac,wvs,frs,newfrs,
                  num_patches,redistr_table,max_n_patch,cg_over,
                  swap_over_function, new_over_function,
                  overparam,hyp_reg_vec);
        if (debugging("hyp_amr"))
            printf("Done with hyp_amr_split_driver, status = %d\n",status);

        if (pp_min_status(status) != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
              " failed after hyp_amr_split_driver()\n");
            clean_up(ERROR);
        }

        if(YES == wvs[0]->use_overture_state)
        { 
            old_over_function->destroy();
            delete old_over_function;
            swap_over_function->destroy();
            delete swap_over_function;
        }  

        free(newfrs);

        printf("_________________\n\n");
        printf("Leaving overture_hyp_amr() total patches = %d\n",num_patches);

        DEBUG_LEAVE(overture_hyp_amr)
        return status;
}  

LOCAL void copy_gridfunction_from_uu1(
        const doubleCompositeGridFunction&  uu1,
        const CompositeGrid&                cg,
        int                                 R,
        doubleCompositeGridFunction&        u)
{
        Index I1,I2,I3;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++)
        {
            int base1,base2,bound1,bound2;
            getIndex(cg[grid].dimension(),I1,I2,I3);
            base1 = I1.getBase();
            base2 = I2.getBase();
            bound1 = I1.getBound();
            bound2 = I2.getBound();
            int ix,iy;
            for (ix = base1; ix <= bound1; ix++)
            {
                for (iy = base2; iy <= bound2; iy++)
                {
                    u[grid](ix,iy,0) = uu1[grid](R,ix,iy,0);
                }
            }
        }
}


/* We choose to use density as bases where the
 * error estimate builds on.
 */
/* Overture should not be responsible for the
 * interpolation of states between old set of patches and
 * newly created set of patches.
 * THis is because the overture does not know antything about the front
 * the FrontTier should take care of this task.
 * However, to see how the Overture does the interpolation,
 * see before AMR_src.011903.tar.gz.
 */
LOCAL   int overture_perform_regrid(
        CHART           *chart,
        Printplot       *prt,
        CompositeGrid                  *rect_over_grid,
        doubleCompositeGridFunction    *rect_over_function,
        CompositeGrid                  **new_cg,
        doubleCompositeGridFunction    **new_cg_function,
        int                            *level_array)
{
        int                             i, j;
        Regrid                          regrid;
        CompositeGrid                   cga[2]; /*cga old and new grids(ring buffer) */
        Overparam                       *overparam;
        int                             levels, num_float, refined_levels;
        size_t                          sizest;
        doubleCompositeGridFunction     uu2;
        doubleCompositeGridFunction     Dens, error;
        Ogen                            ogen;
        Index                           I1,I2,I3;

        DEBUG_ENTER(overture_perform_regrid)

        sizest = chart->wave->sizest;
        overparam = chart->overparam;
        levels = overparam->numberOfRefinementLevels;
        real errorThreshold = overparam->errorThreshold; /* error tolerance for regridding. */
        refined_levels = (rect_over_grid)->refinementLevelNumber(
                          rect_over_grid->numberOfComponentGrids()-1);        

        printf(" ------------------------------------------------------------ \n");
        printf("  Adaptive Mesh Refinement regrid for FronTier Solver        \n");
        printf(" ------------------------------------------------------------ \n");
        Range     all;
        CompositeGrid& cgcur = cga[0];
        cgcur = *rect_over_grid;

        /* interpolator for refinements  */
        InterpolateRefinements interp(cgcur.numberOfDimensions());
        interp.setOrderOfInterpolation(3); /* used to be 1 */

        /* used to compute error estimates */
        ErrorEstimator errorEstimator(interp);
        /* set the scale factors for the solution (estimate of solution size) */ 
        RealArray scaleFactor(1);
        scaleFactor=1.;
        errorEstimator.setScaleFactor(scaleFactor);
        /* set the default number of smoothing steps for smoothing the error */
        errorEstimator.setDefaultNumberOfSmooths(overparam->numberOfSmooths);

        /* overlapping grid interpolator */
        Interpolant interpolant;
        interpolant.setImplicitInterpolationMethod(Interpolant::iterateToInterpolate);
        interpolant.setInterpolateRefinements(interp);

        interpolant.updateToMatchGrid(cgcur,1);  /* error estimator needs to interpolate */
                                                 /* orignal was 0, 010903 change to 1 */

        /* compute an error estimate of density. */

        CompositeGridOperators op(cgcur);
        op.updateToMatchGrid(cgcur);
        error.updateToMatchGrid(cgcur);
        error.setOperators(op);

        if( refined_levels == (levels - 1))
        {
            if(YES == chart->wave->use_overture_state)
            {
                Dens.updateToMatchGrid(cgcur);
                Dens.setOperators(op);
                copy_gridfunction_from_uu1(*rect_over_function, cgcur, 0, Dens);
                errorEstimator.computeAndSmoothErrorFunction(Dens, error);
            } 
            else
                errorEstimator.computeAndSmoothErrorFunction(*rect_over_function,
                          error);
            printf("Error: min=%e, max=%e in regrid by Dens\n",min(error),max(error));
       
            if(overparam->extra_buf != 0)
                update_error(&error, cgcur, overparam);

            CompositeGrid & cgNew = cga[1];
            regrid.setEfficiency(overparam->efficiency);
            regrid.setRefinementRatio(overparam->refinementRatio);
            regrid.regrid(cgcur,cgNew, error, errorThreshold,
                overparam->numberOfRefinementLevels-1, overparam->baseLevel);
            regrid.printStatistics(cgNew);
            ogen.updateRefinement(cgNew);
            cgNew.update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEmask);

            if(debugging("regrid"))
            {
                printf("New generated grid position\n");
                show_grids_info(cgNew);
                printf("Old grid position\n");
                show_grids_info(cgcur);

                printf("Generated total new patch numbers =%d\n", cgNew.numberOfComponentGrids());
                printf("numberOfRefinementLevels =%d\n", overparam->numberOfRefinementLevels);
            }
            printf("New generated grid position\n");
            show_grids_info(cgNew);
            *new_cg = new CompositeGrid();
            **new_cg = cgNew;
            /* update the operators and Interpolant */
            op.updateToMatchGrid(**new_cg);
            interpolant.updateToMatchGrid(**new_cg, 1);
        }
        else
        {
            perform_overture_init_mesh_refinement(rect_over_grid,
               rect_over_function, new_cg, overparam, refined_levels);
        } 
        
        if(YES==chart->wave->use_overture_state)
        {
            *new_cg_function = new
              doubleCompositeGridFunction(**new_cg,sizest/sizeof(double)+1,all,all,all);
        } 
        if(debugging("regrid"))
        {
            printf("In overture_regrid()\n");
            printf("Generated new cg_grid = %d\n", *new_cg);
            printf("Generated new cg_over_function = %d\n",*new_cg_function);
            printf("print new_cg->numberOfComponentGrids() = %d\n",
                (*new_cg)->numberOfComponentGrids());
        }

        printf(" ------------------------------------------------------------------- \n");
        printf("  NumberOfComponentGrids = %d, RefinementLevelNumber = %d\n",
                 (*new_cg)->numberOfComponentGrids(), 
                 (*new_cg)->refinementLevelNumber((*new_cg)->numberOfComponentGrids()-1)+1);
        for(j = 0; j < (*new_cg)->numberOfComponentGrids(); j++)
            printf("  refine patch [%d], level = %d\n", j, 
                    (*new_cg)->refinementLevelNumber(j)); 

        printf("---------------------------------------------------------- \n");
        printf(" END Adaptive Mesh Refinement regrid for FronTier Solver  \n");
        printf("---------------------------------------------------------- \n");

        DEBUG_LEAVE(overture_perform_regrid)
        return  (*new_cg)->numberOfComponentGrids();
}  

LOCAL   void update_error(
        doubleCompositeGridFunction  *error,
        const CompositeGrid&      cg,
        Overparam           *overparam)
{
        int                 grid;
        int                 imin[3], imax[3];  
        int                 gimin[3], gimax[3];  
        int                 i, ix, iy;  
        int                 seg_min, seg_max;  
        int                 start, end;
        int                 set_start, set_end;  
        double               LL[3], UU[3];
        Index               I1,I2,I3;
        double               tol = overparam->errorThreshold;
        int                 aug = overparam->extra_buf;  

        if(debugging("regrid"))
        {
            for(grid = 0; grid < (cg).numberOfComponentGrids(); grid++)
            {
                if(grid != 1) continue;  
                getIndex((cg)[grid].dimension(),I1,I2,I3);
                gimin[0] = I1.getBase();
                gimin[1] = I2.getBase();
                gimax[0] = I1.getBound();
                gimax[1] = I2.getBound();
                LL[0] = (cg)[grid].vertex()(gimin[0],gimin[1],0,axis1);
                LL[1] = (cg)[grid].vertex()(gimin[0],gimin[1],0,axis2);
                UU[0] = (cg)[grid].vertex()(gimax[0],gimax[1],0,axis1);
                UU[1] = (cg)[grid].vertex()(gimax[0],gimax[1],0,axis2);
                printf("Old Grid[%d], imin[%d,%d], imax[%d,%d] LL<%g,%g> UU<%g,%g>\n",
                    grid,gimin[0],gimin[1],gimax[0],gimax[1],LL[0],LL[1],UU[0],UU[1]);
                printf("Old Grid[%d] Print ERROR\n", grid);
                (*error)[grid].display("density");
            }   
        }

        for(grid = 1; grid < (cg).numberOfComponentGrids(); grid++)
        {
            getIndex((cg)[grid].dimension(),I1,I2,I3);
            gimin[0] = I1.getBase();
            gimin[1] = I2.getBase();
            gimax[0] = I1.getBound();
            gimax[1] = I2.getBound();
            getIndex((cg)[grid].indexRange(),I1,I2,I3);
            imin[0] = I1.getBase();
            imin[1] = I2.getBase();
            imax[0] = I1.getBound();
            imax[1] = I2.getBound();
            
            /* x-sweep, error augment */   
            for(iy = imin[1]; iy <= imax[1]; iy++)
            {
                seg_min = imin[0];
                /* The end point is not processed,
                 * guess not important now.
                 */  
                while(seg_min < imax[0])
                { 
                    set_start = NO;
                    set_end = NO; 
                    for(seg_max = seg_min; seg_max <= imax[0]; seg_max++)
                    {
                        if( fabs((*error)[grid](seg_max,iy,0)) > tol &&
                            (! set_start))
                        {
                            start = seg_max;  
                            set_start = YES;  
		        }  
                        if(set_start == YES && 
                           fabs((*error)[grid](seg_max,iy,0)) < tol)
                        {
                            end = seg_max; 
                            set_end = YES;  
                            break;  
                        }  
                    }  

                    if(set_start == YES)
                    {
                        if((start - aug) < gimin[0])
                            i = gimin[0];
                        else 
                            i = start - aug; 
                        for(int tmpi = i; tmpi < start; tmpi++)
                            (*error)[grid](tmpi,iy,0) += tol * 2.0;
                    } 
                    if(set_end == YES)
                    {
                        if((end + aug) > gimax[0])
                            i = gimax[0];
                        else 
                            i = end + aug; 
                        for(int tmpi = end+1; tmpi <= i; tmpi++)
                            (*error)[grid](tmpi,iy,0) += tol * 2.0;
                        seg_max += (aug + 1);  
                    } 
                    seg_min = seg_max;   
                }
            } 

            /* y-sweep, error augment */   
            for(ix = imin[0]; ix <= imax[0]; ix++)
            {
                seg_min = imin[1];
                while(seg_min < imax[1])
                { 
                    set_start = NO;
                    set_end = NO; 
                    for(seg_max = seg_min; seg_max <= imax[1]; seg_max++)
                    {
                        if( fabs((*error)[grid](ix,seg_max,0)) > tol &&
                            (! set_start))
                        {
                            start = seg_max;  
                            set_start = YES;  
		        }  
                        if(set_start == YES && 
                           fabs((*error)[grid](ix,seg_max,0)) < tol)
                        {
                            end = seg_max; 
                            set_end = YES;  
                            break;  
                        }  
                    }  

                    if(set_start == YES)
                    {
                        if((start - aug) < gimin[1])
                            i = gimin[1];
                        else 
                            i = start - aug; 
                        for(int tmpi = i; tmpi < start; tmpi++)
                            (*error)[grid](ix,tmpi,0) += tol * 2.0;
                    } 
                    if(set_end == YES)
                    {
                        if((end + aug) > gimax[1])
                            i = gimax[1];
                        else 
                            i = end + aug; 
                        for(int tmpi = end+1; tmpi <= i; tmpi++)
                            (*error)[grid](ix,tmpi,0) += tol * 2.0;
                        seg_max += (aug + 1);  
                    } 
                    seg_min = seg_max;   
                }
            } 
        }   
}


LOCAL   void d_overture_injection_in_regridding(
        Wave      **oldwvs,
        Front     **oldfrs, 
        Wave      **wvs,
        Front     **frs,
        Wv_on_pc  **redistr_table,
        Overparam *overparam,
        int       max_n_patch)
{

        int       i, k;
        int       patch_num = oldwvs[0]->totalNumberOfPatches;
        int       *iperm;
        int       dim = wvs[0]->rect_grid->dim;  
        int       step = frs[0]->step;  

        DEBUG_ENTER(d_overture_injection_in_regridding)

        iperm = set_iperm(step,dim);  
        /*  set the states of the finest new patches */
        (*wvs[0]->overture_injection_after_repatch)(oldwvs,oldfrs,wvs,frs);

        for(i = 0; i < dim; i++)
        {
            h_scatter_patch_states(wvs,frs,overparam,
              redistr_table,max_n_patch,iperm,i);
        } 
        DEBUG_LEAVE(d_overture_injection_in_regridding)
} 

LOCAL   void d_overture_injection_after_repatch(
        CHART     **roots,
        Wave      **wvs,
        Front     **frs)
{
        int       i, k;
        int       patch_num = roots[0]->wave->totalNumberOfPatches;
        int       *iperm;
        int       dim = wvs[0]->rect_grid->dim;  
        int       step = frs[0]->step;  
        Wave      **oldwvs;
        Front     **oldfrs;

        DEBUG_ENTER(d_overture_injection_after_repatch)

        iperm = set_iperm(step,dim);  

        uni_array(&oldfrs,patch_num,sizeof(Front*));
        uni_array(&oldwvs,patch_num,sizeof(Wave*));

        for(i = 0; i < patch_num; i++)
        {
            oldfrs[i] = roots[i]->front;
            oldwvs[i] = roots[i]->wave;
        }

        /*  set the states of the finest new patches */
        (*wvs[0]->overture_injection_after_repatch)(oldwvs, oldfrs, wvs, frs);
        free_these(1, oldfrs); free_these(1, oldwvs);

        /* parallel part for interior states to buffer*/
        for(i = 0; i < wvs[0]->totalNumberOfPatches; i++)
        {
            for(k = 0; k < dim; k++)
            {
                if (not scatter_states(wvs[i],frs[i],iperm,k))
                {
                    printf("ERROR in d_overture_injection_after_repatch, exit\n");
                    screen("scatter_states() failed for patch[%d]\n",i);
                    clean_up(ERROR);
                }
            } 
        }

        (*wvs[0]->overture_fill_patch_amr_buffer)(wvs, frs);

        if(debugging("injection_repatch_st"))
        {
            Index I1,I2,I3;
            doubleCompositeGridFunction Dens, *cg_f; 
            CompositeGrid *cg = (CompositeGrid*)wvs[0]->cg_over;
            cg_f = (doubleCompositeGridFunction*)wvs[0]->cg_over_function; 
 
            Dens.updateToMatchGrid(*cg);

            copy_gridfunction_from_uu1(*cg_f,*cg,0,Dens);
            printf("In d_overture_injection_after_repatch() after set buffer\n");
            for(int grid = 0; grid < wvs[0]->totalNumberOfPatches; grid++)
            {
                getIndex((*cg)[grid].indexRange(),I1,I2,I3);
                printf("patch[%d] x dir range<%d,%d>, y dir range<%d,%d>\n",
                     grid, I1.getBase(), I1.getBound(), I2.getBase(), I2.getBound());
                Dens[grid].display("show dens with buffer");
            }
            screen("exit in d_overture_injection_after_repatch()\n");
            clean_up(ERROR);
        }
        DEBUG_LEAVE(d_overture_injection_after_repatch)
}

LOCAL   CHART **amr_regridding(
        CHART           **chart,
        Printplot       *prt,
        Wv_on_pc        **redistr_table,
        int             *max_n_patch,
        int             *num_patches,
        Wv_on_pc        ***out_redistr_table,
        int             *out_max_n_patch,
        int             *out_num_patches)
{
        CHART           **outchart;
        Wave            **newwvs, **tmpwvs, **wvs, *wave;
        Wave            **outwvs; 
        Front           **outfrs;  
        Front           **newfrs, **tmpfrs, **frs, *front;
        int             total_patch, dim; 
        Overparam       *overparam;
        int             levels, int_level;
        int             *level_array;
        COMPONENT       *patch_comp; 
        RECT_GRID       *rgs;
        int             status = GOOD_STEP;
        boolean            stat = FUNCTION_SUCCEEDED;
        Wv_on_pc        **new_redistr_table; 
        int             i, new_num_patches;
        int             new_max_n_patch, bal_n;
        int             old_bal_n, myid; 
        INTERFACE       *sav_intfc;
        boolean            sav_copy;  

        CompositeGrid   *new_cg, *cg_over;
        doubleCompositeGridFunction  *new_cg_function;
        doubleCompositeGridFunction old_u; 
        CompositeGridOperators op(*(
                (CompositeGrid*)chart[0]->wave->cg_over)); 
      
        DEBUG_ENTER(amr_regridding)

        uni_array(&frs,*num_patches,sizeof(Front*));
        uni_array(&wvs,*num_patches,sizeof(Wave*));

        for (i = 0; i < *num_patches; i++)
        {
            frs[i] = chart[i]->front;
            wvs[i] = chart[i]->wave;
        }
        reinstore_undistribute_patch(&tmpwvs,&tmpfrs,
           wvs,frs,redistr_table,*num_patches,*max_n_patch); 

        /* Free patches belonging to other proces */  
        myid = pp_mynode(); 
        old_bal_n = 0;
        for(i = 0; i < *max_n_patch; i++)
        {
            if(-1 == redistr_table[myid][i].wv_id) continue;
            if(myid == redistr_table[myid][i].pc_id)
                old_bal_n++;
        }
        sav_intfc = current_interface();
        sav_copy = copy_intfc_states();
        for(i = old_bal_n; i < *num_patches; i++)
        {
            free(MaxFrontSpeedState(frs[i]));
            free(MaxFrontSpeedCoords(frs[i]));
            free(MaxFrontSpeed(frs[i]));
            MaxFrontSpeed(frs[i]) = NULL;
            deep_free_front(frs[i]);
            free(MaxWaveSpeedState(wvs[i]));
            free(MaxWaveSpeedCoords(wvs[i]));
            free(MaxWaveSpeed(wvs[i]));
            MaxWaveSpeed(wvs[i]) = NULL;
            free_wave_pointers(wvs[i]);
            free_wave(wvs[i]);
        } 

        set_current_interface(sav_intfc);
        set_copy_intfc_states(sav_copy);

        total_patch = tmpwvs[0]->totalNumberOfPatches;  
        overparam = chart[0]->overparam;
        levels = chart[0]->overparam->numberOfRefinementLevels;
        wave = chart[0]->wave;
        front = chart[0]->front;
        dim = front->rect_grid->dim;
        uni_array(&level_array,levels,sizeof(int));

        if(NO == tmpwvs[0]->use_overture_state)
        {
            cg_over = (CompositeGrid*)wvs[0]->cg_over; 
            old_u.updateToMatchGrid(*cg_over); 
            old_u.setOperators(op);
            for(i = 0; i < total_patch; i++)
            {
                (*tmpwvs[i]->trans_wv_st_to_overfunc)(tmpwvs[i],
                   tmpfrs[i],(POINTER)(&old_u),i,0);
                (*tmpwvs[i]->fill_root_extr_overture_cell_st_from_wv)(
                   tmpwvs[i],tmpfrs[i],(POINTER)(&old_u),i,0); 
            }    
        } 

        new_num_patches = overture_perform_regrid(chart[0],prt,
           cg_over,&old_u,&new_cg,&new_cg_function,level_array);
        if (debugging("regridding"))
        {
            printf("after overture_perform_regrid, numberOfComponentGrids %d\n",
                new_cg->numberOfComponentGrids());
        }

        uni_array(&newfrs,new_num_patches,sizeof(Front*));
        uni_array(&newwvs,new_num_patches,sizeof(Wave*));
        uni_array(&patch_comp,new_num_patches,sizeof(COMPONENT));

        set_patch_components(wave,new_cg,new_num_patches,patch_comp);
        for(i = 0; i < new_num_patches; i++)
        {
            if(0 != i)
                newfrs[i] = deep_copy_front(front);
            else 
                newfrs[i] = copy_front(front);  
            newwvs[i] = copy_wave(wave);
            rgs = newwvs[i]->rect_grid = newfrs[i]->rect_grid;
            rgs->dim = dim;
            rgs->Remap = wave->rect_grid->Remap;   
            newfrs[i]->patch_level = newwvs[i]->patch_level = 0;
            newfrs[i]->NumberOfLevels = newwvs[i]->NumberOfLevels = levels;
            if(0 != i)
            {
                newfrs[i]->patch_level = newwvs[i]->patch_level
                    = new_cg->refinementLevelNumber(i);  
            } 
            if (i == 0)
                copy_rect_grid(rgs,wave->rect_grid);
            else
                set_patch_rect_grid(rgs,wave->rect_grid,
                     overparam,new_cg,newfrs[i]->patch_level,i);
            newwvs[i]->patch_component = patch_comp[i];
            newfrs[i]->totalNumberOfPatches = newwvs[i]->totalNumberOfPatches
                    = new_num_patches;
            newfrs[i]->patch_number = newwvs[i]->patch_number = i;

            newfrs[i]->cg_over = newwvs[i]->cg_over = (POINTER)new_cg;
            if(YES==newwvs[i]->use_overture_state)
                newwvs[i]->cg_over_function = (POINTER)new_cg_function;
            else
                newwvs[i]->cg_over_function = NULL;
            wave_of_front(newfrs[i]) = newwvs[i]; 
        }

        for (i = 0; i < new_num_patches; i++)
        {
            newwvs[i]->pd_flag = newfrs[i]->pd_flag;;
            set_bdry_patch_flag(newwvs[i],newwvs[0],overparam);
        }
 
        for (i = 1; i < new_num_patches; i++)
        {
            set_patch_front(front,newfrs[i],newfrs[i]->rect_grid,i);

            set_amr_intfc_tol(newfrs[i]->interf,
                pow(2.0,-(newfrs[i]->NumberOfLevels-1.0-newfrs[i]->patch_level)));

            if ((stat = clip_patch_front(newfrs[i],NO)) == FUNCTION_FAILED)
            {
                (void) printf("ERROR in amr_regridding(), "
                                "clip_patch_front()[%d] failed \n",i);
                break;
            }
            set_amr_intfc_tol(newfrs[i]->interf,
                pow(2.0,newfrs[i]->NumberOfLevels-1.0-newfrs[i]->patch_level));

            if((status = init_hyp_solution_function(newwvs[i],newfrs[i])) 
                != GOOD_STEP)
            {
                screen("ERROR in amr_regridding(), "
                   "init_hyp_solution_function() failed on patch[%d]\n",i);
                stat = FUNCTION_FAILED;
                break;
            }

            if(debugging("amr_fr_regridding"))
            {
                printf("______________________________\n");
                printf("THE AMR INIT FRONT[%d]\n", i);
                print_interface(frs[i]->interf);
                show_intfc_states(frs[i]->interf);
            }
        }


        if (pp_min_status(stat) == FUNCTION_FAILED)
        {
            printf("ERROR amr_regridding(),"
                       " failed after clip_patch_front\n");
            clean_up(ERROR);
        }

        new_redistr_table = redistr_computing_load(newwvs,newfrs,
            new_num_patches,&new_max_n_patch);
        allgather_patch_info(new_redistr_table,newfrs,overparam,
          new_num_patches,new_max_n_patch);
        geomview_amr_fronts_plot2d(basename(prt->outfile),newfrs,
                new_num_patches,new_redistr_table,new_max_n_patch); 

        d_overture_injection_in_regridding(tmpwvs,tmpfrs,newwvs,
          newfrs,new_redistr_table,overparam,new_max_n_patch);  
        tmp_save_overture_show_files(prt,newwvs,newfrs,
                   new_num_patches,newfrs[0]->step); 

        ((CompositeGrid*) tmpwvs[0]->cg_over)->destroy();
        delete ((CompositeGrid*) tmpwvs[0]->cg_over);

        set_current_interface(newfrs[0]->interf);
        tmpfrs[0]->interf = NULL;
        free_front(tmpfrs[0]);
        clear_wave_pointers(tmpwvs[0]);
        free_wave(tmpwvs[0]);
        for(i = 1; i < total_patch; i++)
        {
            if(i < old_bal_n)
            {
                free(MaxFrontSpeedState(tmpfrs[i]));
                free(MaxFrontSpeedCoords(tmpfrs[i]));
                free(MaxFrontSpeed(tmpfrs[i]));
                MaxFrontSpeed(tmpfrs[i]) = NULL;
                free(MaxWaveSpeedState(tmpwvs[i]));
                free(MaxWaveSpeedCoords(tmpwvs[i]));
                free(MaxWaveSpeed(tmpwvs[i]));
                MaxWaveSpeed(tmpwvs[i]) = NULL;
            }
            free_wave_pointers(tmpwvs[i]);
            free_wave(tmpwvs[i]);
            deep_free_front(tmpfrs[i]);
        }
        set_copy_intfc_states(sav_copy);

        pp_gsync(); 
        bal_n = perform_redistribute_load(new_redistr_table,newwvs,newfrs,
                   new_num_patches,new_max_n_patch,&outwvs,&outfrs);

        remove_sent_wv_fr_after_redistr(new_redistr_table,new_num_patches,
                      new_max_n_patch,newfrs,newwvs);

        uni_array(&outchart,bal_n,sizeof(CHART*));
        for(i = 0; i < bal_n; i++)
        {
            outchart[i] = copy_chart(chart[0]);
            outchart[i]->grid->rect_grid = outfrs[i]->rect_grid;
            outchart[i]->old_cg_over = outfrs[i]->cg_over;
            if(0 != i)
            {
                MaxWaveSpeed(outwvs[i]) = NULL;
                MaxWaveSpeed(outwvs[i]) =
                 alloc_MaxWaveSpeed(MaxWaveSpeed(outwvs[i]),outwvs[i]);
                /* initialize_max_wave_speed(outwvs[i]); */
                MaxFrontSpeed(outfrs[i]) = NULL;
                MaxFrontSpeed(outfrs[i]) =
                 alloc_MaxFrontSpeed(MaxFrontSpeed(outfrs[i]),
                     outfrs[i]->interf, outfrs[i]->sizest);
                /* initialize_max_front_speed(outfrs[i]); */
            }
            outchart[i]->wave = outwvs[i];
            outchart[i]->front = outfrs[i];
            outchart[i]->front->newfront = outfrs[i];
            chart_of_front(outfrs[i]) = outchart[i];
        }

        for(i = 0; i < bal_n; i++)
            if(exists_interface(outfrs[i]->interf) != YES)
                printf("Invalid Out Interface [%d] in"
                  " amr_regridding, "
                   "Probably already deleted\n", i);

        free_these(2,outfrs,outwvs);
        free_these(2,newfrs, newwvs);  
        free_these(2,patch_comp,level_array); 
        free_these(2,wvs,frs); 
        free(redistr_table);  
        free_these(2,tmpwvs,tmpfrs);  
        for(i = 0; i < *num_patches; i++)
        {
            free(chart[i]->grid);
            free(chart[i]);
        }
        free(chart);  
        *out_max_n_patch = new_max_n_patch;
        *out_num_patches = bal_n; 
        *out_redistr_table = new_redistr_table;  

        set_current_interface(outfrs[0]->interf); 
        printf("LEAVING REGRIDING STEP\n");
        pp_gsync(); 
        DEBUG_LEAVE(amr_regridding)
        return outchart;  
}

LOCAL   CHART **amr_overture_repatch(
        CHART           **roots,
        Printplot       *prt,
        Wv_on_pc        **redistr_table,
        int             *max_n_patch,
        int             *patch_num)
{
        CompositeGrid   *new_cg;
        doubleCompositeGridFunction    *new_cg_function;
        Overparam       *overparam;
        int             num_patches, levels, int_level;
        int             *level_array;

        CHART           **rts;
        Grid            *grid;
        Front           *front, **frs;
        Wave            *wave, **wvs;
        RECT_GRID       *rgs;
        COMPONENT       *patch_comp;
        int             i, j, dim;
        int             step;
        double           dt, time;

        DEBUG_ENTER(amr_overture_repatch)

        overparam = roots[0]->overparam;
        levels = roots[0]->overparam->numberOfRefinementLevels;
        wave = roots[0]->wave;
        front = roots[0]->front;
        dim = front->rect_grid->dim;
        uni_array(&level_array,levels,sizeof(int));
        dt = roots[0]->grid->dt;
        time = roots[0]->grid->time;
        step = roots[0]->grid->step;

        if (debugging("repatch"))
        {
            printf("before overture_regrid, numberOfComponentGrids = %d\n",
                ((CompositeGrid*)roots[0]->wave->cg_over)->numberOfComponentGrids());
        }
        screen("Overture regrid generated new number patches = %d\n", num_patches);

        if (debugging("repatch"))
        {
            printf("after overture_regrid, numberOfComponentGrids = %d\n",
                new_cg->numberOfComponentGrids());
        }
        uni_array(&rts,num_patches,sizeof(CHART*));
        uni_array(&frs,num_patches,sizeof(Front*));
        uni_array(&wvs,num_patches,sizeof(Wave*));
        uni_array(&patch_comp,num_patches,sizeof(COMPONENT));

        set_patch_components(wave,new_cg,num_patches,patch_comp);
/* THIS LOOP NEED to be redone */  
        for(i = 0; i < num_patches; i++)
        {
            scalar(&rgs,sizeof(RECT_GRID));
            frs[i] = copy_front(front);
            wvs[i] = copy_wave(wave);
            rts[i] = copy_chart(roots[0]);
            rts[i]->grid->rect_grid = frs[i]->rect_grid =
                           wvs[i]->rect_grid = rgs;
            if (debugging("amr_overture_repatch"))
            {
                printf("set patch waves\n");
                printf("wvs[%d]->cg_over_function=%d\n",
                        i,wvs[i]->cg_over_function);
            }
            rgs->dim = dim;
            if (i == 0)
                copy_rect_grid(rgs,wave->rect_grid);
            else
                set_patch_rect_grid(rgs,wave->rect_grid,
                     overparam,new_cg,-1,i);
            frs[i]->NumberOfLevels = wvs[i]->NumberOfLevels = levels;
            frs[i]->patch_level = wvs[i]->patch_level = 0;
            wvs[i]->patch_component = patch_comp[i];
            frs[i]->totalNumberOfPatches = wvs[i]->totalNumberOfPatches
                    = new_cg->numberOfComponentGrids();
            frs[i]->patch_number = wvs[i]->patch_number = i;

            frs[i]->cg_over = wvs[i]->cg_over = (POINTER)new_cg;
            if(YES==wvs[i]->use_overture_state)
                wvs[i]->cg_over_function = (POINTER)new_cg_function;
            else
                wvs[i]->cg_over_function = NULL;  
            rts[i]->wave = wvs[i];
            rts[i]->front = frs[i];
            chart_of_front(frs[i]) = NULL;
        }
        if (debugging("amr_overture_repatch"))
        {
            printf("Print wvs[0] rect_grid\n");
            print_rectangular_grid(wvs[0]->rect_grid);
        }

        for (i = 1; i < num_patches; i++)
        {
            for (int_level = 0; int_level < levels; int_level++)
            {
                if (i >= level_array[int_level])
                {
                    frs[i]->patch_level = wvs[i]->patch_level = int_level;
                }
            }
        }

        if (debugging("regrid"))
        {
            printf("In amr_overture_repatch() print patch level:\n");
            for(i = 0; i < num_patches; i++)
                printf("patch[%d] is on level[%d]\n",i, wvs[i]->patch_level);
        }

        INTERFACE       *sav_intfc;
        boolean            sav_copy;
        sav_intfc = current_interface();
        sav_copy = copy_intfc_states();
        set_size_of_intfc_state(size_of_state(roots[0]->front->interf));
        set_copy_intfc_states(YES);
        if ((frs[0]->interf = copy_interface(roots[0]->front->interf)) == NULL)
        {
            screen("FT_ERROR in amr_overture_repatch()",
                    "copy_interface() failed\n");
            clean_up(ERROR);
        }
        set_current_interface(sav_intfc);
        set_copy_intfc_states(sav_copy);

        if (scatter_front(frs[0]) == FUNCTION_FAILED)
        {
            (void) printf("ERROR in amr_overture_repatch(), "
                "scatter_front() for base frs[0] failed \n ");
            clean_up(ERROR);
        }
        if( debugging("repatch"))
        {
            printf("In amr_overture repatch,"
              " base interface after copy\n");
            print_interface(frs[0]->interf);
        }

        if (init_hyp_solution_function(wvs[0],frs[0]) != GOOD_STEP)
        {
            screen("ERROR in amr_overture_repatch(), "
                   "init_hyp_solution_function() failed for patch[0]\n");
            clean_up(ERROR);
        }

        for (i = 0; i < num_patches; i++)
        {
            set_bdry_patch_flag(wvs[i],wvs[0],overparam);
            frs[i]->pd_flag = wvs[i]->pd_flag;
        }

        for (i = 1; i < num_patches; i++)
        {
            printf("amr_overture_repatch(), clip patch[%d], level[%d]\n",
              i,wvs[i]->patch_level);
            set_patch_front(frs[0],frs[i],frs[i]->rect_grid,i);
            remove_patch_all_boundary_curves(frs[i]);
            if (clip_patch_front(frs[i],NO) == FUNCTION_FAILED)
            {
                (void) printf("ERROR in amr_overture_repatch(), "
                   "clip_patch_front() failed for patch[%d]\n",i);
                clean_up(ERROR);
            }
            if (debugging("front_repatch"))
            {
                printf("in amr_overture_repatch(), processed patch[%d]\n", i);
                printf("In amr_overture_repatch()"
                    " before delete_subdomain_boundaries()\n");
                print_interface(frs[i]->interf);
            }

            delete_subdomain_boundaries(frs[i]->interf);
            if(init_hyp_solution_function(wvs[i],frs[i]) != GOOD_STEP)
            {
                screen("FT_ERROR in amr_overture_repatch(), "
                    "init_hyp_solution_function() failed for patch[%d]\n",i);
                clean_up(ERROR);
            }
            rts[i]->front = frs[i];
            chart_of_front(frs[i]) = rts[i];
        }

        d_overture_injection_after_repatch(roots,wvs,frs);
       
        if(YES == wave->use_overture_state)
        {
            ((doubleCompositeGridFunction*)wave->cg_over_function)->destroy();
            delete (doubleCompositeGridFunction*)wave->cg_over_function;
        } 
        ((CompositeGrid*)wave->cg_over)->destroy();
        delete (CompositeGrid*)wave->cg_over;

        for(i = 0; i < *patch_num; i++)
        {
          /*        free(roots[i]->wave->pd_flag); */
            free_wave_pointers(roots[i]->wave);
            free_wave(roots[i]->wave);
            if (roots[i]->front->interf != NULL)
                (void) delete_interface(roots[i]->front->interf);
            free_front(roots[i]->front);
            /*      free(roots[i]->grid); */
        }
        free_these(1,roots);
        free_these(1, level_array);
        free_these(1, frs);
        free_these(1, wvs);
        free_these(1, patch_comp);

        *patch_num = num_patches;

        DEBUG_LEAVE(amr_overture_repatch)
        return  rts;
} 

LOCAL void tmp_save_overture_show_files(
        Printplot    *prt,
        Wave         **wvs,
        Front        **frs, 
        int          num_patches, 
        int          step)
{

        const char   *dname = basename(prt->outfile); 
        const char   *nstep, *nd; 
        doubleCompositeGridFunction *uu1;
        CompositeGrid *cg;
        Range        all;
        char         buf[256],outname[256];
        char         outdir[256];  
        Wave         *wave;
        int          use_overture_state, i;  

        debug_print("driver","Entering tmp_save_overture_show_files()\n");

        sprintf(outdir,"%s/%s",dname,"over_show"); 
        nd = right_flush(pp_mynode(),4);
        sprintf(buf,"%s/%s_%s.%s",outdir,"over","den_show",nd);
        nstep = right_flush(step,7); 
        sprintf(outname,"%s.ts%s",buf,nstep);

        if (create_directory(dname,YES) == FUNCTION_FAILED)
        {
            (void) printf("WARNING in save_overture_show_files(), directory "
                         "%s doesn't exist and can't be created\n",dname);
            return;
        }
        if (create_directory(outdir,YES) == FUNCTION_FAILED)
        {
            (void) printf("WARNING in save_overture_show_files(), directory "
                         "%s doesn't exist and can't be created\n",outdir);
            return;
        }
        wave = wvs[0]; 
        use_overture_state = wave->use_overture_state;  
        cg = (CompositeGrid*)wave->cg_over;
        CompositeGridOperators op(*cg);   
        Ogshow show(outname);
        show.saveGeneralComment("The Density solution");
 
        if(use_overture_state)
        {
            uu1 = (doubleCompositeGridFunction*)wave->cg_over_function;
        }
        else
        {
            doubleCompositeGridFunction Dens;
            Dens.updateToMatchGrid(*cg);
            for(i = 0; i < num_patches; i++)
            {
                (*wave->trans_wv_st_to_overfunc)(wvs[i],NULL,(POINTER)(&Dens),i,0);
                (*wave->fill_root_extr_overture_cell_st_from_wv)(
                   wvs[i],frs[i],(POINTER)(&Dens),i,0); 
            }   
            show.startFrame();
            show.saveSolution(Dens); 
            printf("Saved Dens\n");  
        }

        debug_print("driver","Leaving tmp_save_overture_show_files()\n");
}  

LOCAL CHART **perform_amr_initialization(
        INIT_DATA       *init,
        INIT_PHYSICS    *ip,
        int             *patch_num,
        Wv_on_pc        ***redistr_table,
        int             *max_n_patch)
{
        CHART           **rts, *root;
        Grid            *grid;
        Front           *front,**frs,**nfrs;
        Wave            *wave,**wvs,**nwvs;
        RECT_GRID       *rgs;
        Overparam       *overparam;
        int             num_patches;
        int             i;
        int             status = GOOD_STEP;
        boolean            stat = FUNCTION_SUCCEEDED;
        size_t          sizest;
        int             num_float;
        int             use_overture_state;
        int             overture_init_step;
        int             dim;
        int             levels;
        int             *level_array;
        CompositeGrid                   *cg_over;
        doubleCompositeGridFunction     *cg_over_function;
        int             bal_n;  

        DEBUG_ENTER(perform_amr_initialization);

        if (debugging("amr"))
            printf("Entering perform_amr_initialization()\n");

        root = ip->root;
        overparam = root->overparam;

        wave = root->wave;
        front = root->front;
        dim = front->rect_grid->dim;
        levels = ip->root->overparam->numberOfRefinementLevels;
	uni_array(&level_array,levels,sizeof(int));

        num_patches = overture_init_mesh_refinement(overparam,wave,front,level_array);
        cg_over = (CompositeGrid*)wave->cg_over;

        use_overture_state = wave->use_overture_state;
        front->use_overture_state = wave->use_overture_state;
        if(use_overture_state) 
	{	
            cg_over_function = (doubleCompositeGridFunction*)wave->cg_over_function;
	}    
        overture_init_step = wave->overture_init_step = NO;

        uni_array(&frs,num_patches,sizeof(Front*));
        uni_array(&wvs,num_patches,sizeof(Wave*));

        set_patch_wv_fr_before_redistr(num_patches,root,frs,wvs,level_array);  
        for (i = 1; i < num_patches; i++)
        {
            set_patch_front(front,frs[i],frs[i]->rect_grid,i);

            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,-(frs[i]->NumberOfLevels-1.0-frs[i]->patch_level)));

            if ((stat = clip_patch_front(frs[i],NO)) == FUNCTION_FAILED)
            {
                (void) printf("ERROR in perform_amr_initialization(), "
                                "clip_patch_front()[%d] failed \n",i);
                break;
            }
            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,frs[i]->NumberOfLevels-1.0-frs[i]->patch_level)); 

            if((status = init_hyp_solution_function(wvs[i],frs[i])) != GOOD_STEP)
            {
                screen("ERROR in perform_amr_initialization(), "
                   "init_hyp_solution_function() failed on patch[%d]\n",i);
                stat = FUNCTION_FAILED; 
                break; 
            }  
            if(debugging("amr_front_init") )
            {
                printf("______________________________\n");  
                printf("THE AMR INIT FRONT[%d]\n", i); 
                print_rectangular_grid(frs[i]->rect_grid);  
            } 
        }

        if (pp_min_status(stat) == FUNCTION_FAILED)
        {
            printf("ERROR perform_amr_initialization(),"
                       " failed after clip_patch_front\n"); 
            clean_up(ERROR);  
        }  

        stat = (*wave->overture_init_interpolation_coarse_to_fine)(wvs,frs);   

        if (pp_min_status(stat) == FUNCTION_FAILED)
        {
            printf("ERROR perform_amr_initialization(),"
               " failed after overture_init_interpolation_coarse_to_fine\n"); 
            clean_up(ERROR);  
        }  

        *redistr_table = redistr_computing_load(wvs,frs,num_patches,max_n_patch); 
        allgather_patch_info(*redistr_table,frs,overparam,num_patches,*max_n_patch);  

        geomview_amr_fronts_plot2d(basename(ip->prt->outfile),frs,num_patches,
           *redistr_table,*max_n_patch); 
        tmp_save_overture_show_files(ip->prt,wvs,frs,num_patches,root->grid->step);

        bal_n = perform_redistribute_load(*redistr_table,wvs,frs,num_patches,
                      *max_n_patch,&nwvs,&nfrs);  

        if(debugging("perform_redistribute_load"))
        {
            int numnodes = pp_numnodes();
        } 

        remove_sent_wv_fr_after_redistr(*redistr_table,num_patches,
                      *max_n_patch,frs,wvs); 
        free_these(1,frs);
        free_these(1,wvs);
        free_these(1,level_array);
        uni_array(&rts,bal_n,sizeof(CHART*));   
        for(i = 0; i < bal_n; i++)
        {
            rts[i] = copy_chart(root);
            rts[i]->grid->rect_grid = nfrs[i]->rect_grid;
            rts[i]->old_cg_over = nfrs[i]->cg_over;
            if(0 != i)
            {
                MaxWaveSpeed(nwvs[i]) = NULL;
                MaxWaveSpeed(nwvs[i]) =
                 alloc_MaxWaveSpeed(MaxWaveSpeed(nwvs[i]),nwvs[i]);
               initialize_max_wave_speed(nwvs[i]);
                MaxFrontSpeed(nfrs[i]) = NULL;
                MaxFrontSpeed(nfrs[i]) =
                 alloc_MaxFrontSpeed(MaxFrontSpeed(nfrs[i]),
                     nfrs[i]->interf, nfrs[i]->sizest);
               initialize_max_front_speed(nfrs[i]);
            }
            rts[i]->wave = nwvs[i];
            rts[i]->front = nfrs[i];
            rts[i]->front->newfront = nfrs[i];
            chart_of_front(nfrs[i]) = rts[i]; 
        }

        free_these(1,nfrs);
        free_these(1,nwvs);
        
        if (debugging("amr"))
            printf("Leaving perform_amr_initialization()\n");
        DEBUG_LEAVE(perform_amr_initialization)

        /*  *patch_num = num_patches; */
        *patch_num = bal_n; 
        return rts;
}  

LOCAL void remove_sent_wv_fr_after_redistr(
	Wv_on_pc   **redistr_table,
	int        num_patches,
	int        max_n_patch, 
        Front      **frs,
        Wave       **wvs)
{
        int         source, i;
        int         patch_id;  
        INTERFACE  *current_intfc;  
      
        DEBUG_ENTER(remove_sent_wv_fr_after_redistr)
        source = pp_mynode(); 

        for(i = 0; i < max_n_patch; i++)
        {
            patch_id = -100;  
            if(-1 == redistr_table[source][i].wv_id) continue; 
            if(source == redistr_table[source][i].pc_id) continue; 
            patch_id = redistr_table[source][i].wv_id; 

            free_wave_pointers(wvs[i]);
	    free_wave(wvs[i]);
            deep_free_front(frs[i]);
        }
        DEBUG_LEAVE(remove_sent_wv_fr_after_redistr)
} 


LOCAL void set_patch_wv_fr_before_redistr(
	int        num_patches,
        CHART      *root,
        Front      **frs,
        Wave       **wvs, 
        int        *level_array)
{
        Wave       *wave = root->wave;
        Front      *front = root->front;
        int        dim = front->rect_grid->dim;
        Overparam  *overparam = root->overparam; 
        CompositeGrid  *cg_over;
        int        levels = overparam->numberOfRefinementLevels;
        RECT_GRID  *rgs; 
        int        i, int_level;
        COMPONENT  *patch_comp;  

        uni_array(&patch_comp,num_patches,sizeof(COMPONENT));

        cg_over = (CompositeGrid*)wave->cg_over;
        set_patch_components(wave, cg_over,
                wave->totalNumberOfPatches,patch_comp);
        for(i = 0; i < num_patches; i++)
        {
            frs[i] = deep_copy_front(front);
            wvs[i] = copy_wave(wave);
            rgs = wvs[i]->rect_grid = frs[i]->rect_grid;
            rgs->dim = dim;
            rgs->Remap = wave->rect_grid->Remap; 
            frs[i]->patch_level = wvs[i]->patch_level = 0;
            if(0 != i)
            {
                wvs[i]->patch_component = patch_comp[i];
                frs[i]->patch_level = wvs[i]->patch_level = 
			      (cg_over)->refinementLevelNumber(i); 
            }
            if (i == 0)
                copy_rect_grid(rgs, wave->rect_grid);
            else
                set_patch_rect_grid(rgs,wave->rect_grid,
                       overparam,cg_over,frs[i]->patch_level,i);
            frs[i]->NumberOfLevels = wvs[i]->NumberOfLevels = levels;
            frs[i]->totalNumberOfPatches = wvs[i]->totalNumberOfPatches = num_patches;
            frs[i]->patch_number = wvs[i]->patch_number = i;
            frs[i]->cg_over = wave->cg_over;
            chart_of_front(frs[i]) = NULL;
            wave_of_front(frs[i]) = wvs[i]; 
        }
        for (i = 0; i < num_patches; i++)
        {
            wvs[i]->pd_flag = frs[i]->pd_flag;;
            set_bdry_patch_flag(wvs[i],wvs[0],overparam);
        }
        free(patch_comp); 
}

LOCAL int perform_redistribute_load(
	Wv_on_pc    **redistr_table,
        Wave        **wvs,
        Front       **frs,
        int         num_patches,
        int         max_n_patch,
        Wave        ***newwvs,
        Front       ***newfrs)
{
        int         source, i; 
        int         numnodes, myid;
        Wave        **nwvs, *wave;
        Front       **nfrs, *front;
        RECT_GRID   *gr;   
        int         bal_n = 0;  
        int         patch_id;  
        int         nn, sv_nn;  

        DEBUG_ENTER(perform_redistribute_load)

#if defined(__MPI__) 

        myid = pp_mynode();    
        numnodes = pp_numnodes();

        if(debugging("perform_redistribute_load"))
        {
            printf("IN perform_redistribute_load() \n"); 
            printf("_____________________________\n");
            printf("AFTER REDISTRIBUTE, the patch is scattered as:\n");
            for(source = 0; source < numnodes; source++)
            {
                printf("Proc[%d]: ", source);
                for(i = 0; i < max_n_patch; i++)
                {
                    if(-1 != redistr_table[source][i].wv_id)
                        printf(" %2d %2d, ", redistr_table[source][i].pc_id,
                                     redistr_table[source][i].wv_id);
                    else
                       printf("        ");  
                }
                printf("\n");
            }
            printf("_____________________________\n");
        }

        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                if(redistr_table[source][i].pc_id == myid and 
                   redistr_table[source][i].wv_id != -1)
                    bal_n++;   
            }
        }
        if(debugging("perform_redistribute_load"))
        {
            printf("NODE[%d] is ft_assigned with %d patches\n", myid, bal_n);  
        }   

        wave = wvs[0];
        front = frs[0];  
        uni_array(&nfrs,bal_n,sizeof(Front*));
        uni_array(&nwvs,bal_n,sizeof(Wave*));

        nn = 0;   
	nn = copy_and_alloc_wv_fr(wvs,frs,nwvs,nfrs,
		 redistr_table,max_n_patch,myid,numnodes); 
        sv_nn = nn;  

        patch_front_trans(frs, nfrs, redistr_table, myid, nn, max_n_patch, numnodes); 

	set_recv_wv_fr(nwvs, nfrs, sv_nn, bal_n); 

        patch_wave_trans(wvs, nwvs, redistr_table, myid, nn, max_n_patch, numnodes); 
   
        *newwvs = nwvs;
        *newfrs = nfrs;  
#endif /* defined(__MPI__) */  
        DEBUG_LEAVE(perform_redistribute_load)
        return bal_n;  
}


LOCAL void set_recv_wv_fr(
	Wave        **nwvs,
        Front       **nfrs,
        int         st_n,
	int         bal_n)
{
	int         i; 
        RECT_GRID   *cgr, *tgr;
        RECT_GRID   *rgr; 

        DEBUG_ENTER(set_recv_wv_fr)

        for(i = st_n; i < bal_n; i++)
        {
            nfrs[i]->interf->modified = YES;
            nwvs[i]->patch_number = nfrs[i]->patch_number;  
            nwvs[i]->patch_level = nfrs[i]->patch_level;  
            rgr = nfrs[i]->rect_grid;  
            cgr = computational_grid(nfrs[i]->interf);
            copy_rect_grid(cgr,rgr); 
            tgr = &topological_grid(nfrs[i]->interf);
            tgr->Remap.remap = rgr->Remap.remap;
            set_patch_topo_grid(rgr,tgr); 
            if (init_hyp_solution_function(nwvs[i],nfrs[i]) != GOOD_STEP)
            {
                screen("ERROR in set_recv_wv_fr(), "
                       "init_hyp_solution_function() failed for [%d]\n",i);
                clean_up(ERROR);
            }
        }
        DEBUG_LEAVE(set_recv_wv_fr)
} 

LOCAL void patch_front_trans(
        Front       **frs,
        Front       **nfrs,
	Wv_on_pc    **redistr_table,
        int         myid,   
        int         nn,  
        int         max_n_patch,
        int         numnodes)
{
        int         source, i;
        int         patch_id;  
       
        DEBUG_ENTER(patch_front_trans)

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

                    send_front_misc(frs[patch_id], 
                        &(wave_of_front(frs[patch_id])->patch_component), 
                        redistr_table[source][i].pc_id);  
                }   
                if(myid == redistr_table[source][i].pc_id)
                {
                    pp_recv(0,source,(POINTER)(&patch_id),sizeof(int)) ;
                    recv_front_misc(nfrs[nn], 
                        &(wave_of_front(nfrs[nn])->patch_component), source); 
                    redistr_table[source][i].front = nfrs[nn];  
                    nn++;  
                }   
                pp_gsync();   
            }
        }

        DEBUG_LEAVE(patch_front_trans)
} 

LOCAL int copy_and_alloc_wv_fr(
	Wave        **wvs,
        Front       **frs,
	Wave        **nwvs,
        Front       **nfrs,
	Wv_on_pc    **redistr_table,
        int         max_n_patch,
	int         myid, 
        int         numnodes)
{
	int        nn, source;  
	int        i, patch_id;  
        int        bal_n;  
        Front      *front = frs[0];
        Wave       *wave  = wvs[0];   
        RECT_GRID  *gr;  
        INTERFACE  *sav_intfc;
        boolean       sav_copy;

        DEBUG_ENTER(copy_and_alloc_wv_fr)

        sav_intfc = current_interface();
        sav_copy = copy_intfc_states();

        /* assembly fronts and waves which are not subject to move */ 

        nn = bal_n = 0;   
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                if(redistr_table[source][i].pc_id == myid and 
                   redistr_table[source][i].wv_id != -1)
                    bal_n++;   
                if(-1 == redistr_table[source][i].wv_id) continue; 
                if(source == redistr_table[source][i].pc_id)
                { 
                    if(myid == source)
                    {
                        patch_id = redistr_table[source][i].wv_id; 
                        nfrs[nn] = frs[patch_id]; 
                        nwvs[nn] = wvs[patch_id];  
                        redistr_table[source][i].front = nfrs[nn];  
                        nn++;  
                    }  
                }  
            }
        }  

        set_size_of_intfc_state(size_of_state(front->interf));
        set_copy_intfc_states(YES);

        /* set the storage for recv patches from other procs */ 
        for(i = nn; i < bal_n; i++)
        {
            nfrs[i] = deep_copy_front(front);
            nfrs[i]->interf = copy_interface(front->interf);
            delete_patch_all_curves(nfrs[i]);  
            nwvs[i] = copy_wave(wave);
            clear_wave_pointers(nwvs[i]); 
            nwvs[i]->rect_grid = nfrs[i]->rect_grid;  
            nwvs[i]->pd_flag = nfrs[i]->pd_flag; 
            wave_of_front(nfrs[i]) = nwvs[i];  
        }

        set_current_interface(sav_intfc);
        set_copy_intfc_states(sav_copy);

        DEBUG_LEAVE(copy_and_alloc_wv_fr)
	return nn;  
}

LOCAL Wv_on_pc **redistr_computing_load(
	Wave        **wvs,
        Front       **frs,
        int         num_patches,
	int         *max_n_patch)
{
        int         i;  
        int         myid, numnodes, source;  
        int         **total_load, *proc_patch_n, *send_load;
	Wv_on_pc    **redistr_table; 
        int         *tmp_total_load;  
        
        DEBUG_ENTER(redistr_computing_load)

        myid = pp_mynode();   
        numnodes = pp_numnodes();

        *max_n_patch = num_patches;  
	pp_global_imax(max_n_patch, 1);  

	if(debugging("redistr_computing_load"))
            printf("NUMBER of Max patches = %d\n", *max_n_patch);   

        uni_array(&tmp_total_load,numnodes*(*max_n_patch),sizeof(int));  
        zero_scalar(tmp_total_load,sizeof(int)*numnodes*(*max_n_patch));  
        uni_array(&proc_patch_n,numnodes,sizeof(int));  

        uni_array(&send_load, (*max_n_patch), sizeof(int));  
        calculate_computing_load(wvs, frs, num_patches, send_load); 

        proc_patch_n[myid] = num_patches;  
        allgather_load_info(tmp_total_load,send_load, proc_patch_n, 
                     num_patches,numnodes,*max_n_patch); 
        free(send_load);
        /* transfer tmp_total_load(uni_array) to total_load (maxtrix). */
        /* easy to use. :)                                          */  
        uni_array(&total_load,numnodes,sizeof(int*));  
        for(i = 0;   i < numnodes; i++)
            uni_array(&total_load[i], (*max_n_patch), sizeof(int));
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < *max_n_patch; i++)
                total_load[source][i] = 
                          tmp_total_load[source*(*max_n_patch)+i];  
        }
        free(tmp_total_load);  

        redistr_table = redistr_load_as(total_load,proc_patch_n,
                              *max_n_patch,numnodes); 

        for(source = 0; source < numnodes; source++)
            free(total_load[source]);  
        free_these(2,total_load,proc_patch_n);

	if(debugging("redistr_computing_load"))
	{	
            printf("_____________________________\n");
            printf("AFTER REDISTRIBUTE, the patch is scattered as:\n");
            for(source = 0; source < numnodes; source++)
            {
                printf("Proc[%d]: ", source);
                for(i = 0; i < *max_n_patch; i++)
                {
                    if(-1 != redistr_table[source][i].wv_id)
                        printf(" %2d %2d, ", redistr_table[source][i].pc_id, 
				     redistr_table[source][i].wv_id);
                    else
                       printf("        "); 
                }
                printf("\n");
            }
	}    

        DEBUG_LEAVE(redistr_computing_load)
	return redistr_table;   	
} 

LOCAL void allgather_patch_info(
        Wv_on_pc    **redistr_table, 
	Front       **frs,
        Overparam   *overparam, 
        int         num_patches,
	int         max_n_patch)
{
	int        len, send_len, source; 
        int        numnodes, my_ic[MAXD], myid; 
	int        i, j, ii, jj;  
        byte       *send_storage = NULL;
        byte       *recv_storage = NULL;
       	byte       *buf;
        POINTER    info;  
        RECT_GRID  *rgr;
        int        dim = frs[0]->rect_grid->dim;
        int        n_reg_nodes;
        int        mylevel, levels, refine; 
        Wv_info    *send_wvinfo, **recv_wvinfo;  
        CompositeGrid   *cg_over;
        Range      all;
        Index      I1,I2,I3;

        DEBUG_ENTER(allgather_patch_info)
     
        numnodes = pp_numnodes(); 
        myid = pp_mynode();
        find_Cartesian_coordinates(myid, frs[0]->pp_grid, my_ic); 

        cg_over = (CompositeGrid*)frs[0]->cg_over;  
        levels = overparam->numberOfRefinementLevels; 
        uni_array(&send_wvinfo, max_n_patch, sizeof(Wv_info));         
        bi_array(&recv_wvinfo, numnodes, max_n_patch, sizeof(Wv_info));         

        rgr = frs[0]->rect_grid;
        for(i = 0; i < num_patches; i++)
        {
            refine = 1;  
            for(ii = 0; ii < frs[i]->patch_level; ii++)
                refine *= overparam->refinementRatio; 
            for(ii = 0; ii < dim; ii++)
                send_wvinfo[i].gmax[ii] = refine*rgr->gmax[ii]; 

            getIndex((*cg_over)[i].indexRange(),I1,I2,I3);
            send_wvinfo[i].wv_level = frs[i]->patch_level;  
            send_wvinfo[i].pc_ic[0] = my_ic[0]; 
            send_wvinfo[i].pc_ic[1] = my_ic[1]; 
            send_wvinfo[i].base[0] = I1.getBase();
            send_wvinfo[i].base[1] = I2.getBase();
            send_wvinfo[i].bound[0] = I1.getBound();
            send_wvinfo[i].bound[1] = I2.getBound();
        }  

	send_len = max_n_patch*sizeof(Wv_info);  
	scalar(&send_storage, send_len); 
        len = send_len*numnodes; 
	scalar(&recv_storage, len); 

        buf = send_storage;  
        for(i = 0; i < max_n_patch; i++)
        {
            info = (POINTER) buf;       
	    ft_assign(info,&send_wvinfo[i],sizeof(Wv_info));   	    
	    buf += sizeof(Wv_info);  
        }

        pp_all_gather(send_storage,send_len,recv_storage,send_len);  
  
        /* reconstruct base  */ 
        buf = recv_storage;  
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                info = (POINTER) buf;       
	        ft_assign(&recv_wvinfo[source][i],info,sizeof(Wv_info));   	    
	        buf += sizeof(Wv_info);  
            }  
        } 

        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                if(-1 == redistr_table[source][i].wv_id) continue; 
                redistr_table[source][i].wv_level = recv_wvinfo[source][i].wv_level;
                redistr_table[source][i].base[0] = recv_wvinfo[source][i].base[0]; 
                redistr_table[source][i].base[1] = recv_wvinfo[source][i].base[1]; 
                redistr_table[source][i].bound[0] = recv_wvinfo[source][i].bound[0]; 
                redistr_table[source][i].bound[1] = recv_wvinfo[source][i].bound[1]; 
                redistr_table[source][i].pc_ic[0] = recv_wvinfo[source][i].pc_ic[0]; 
                redistr_table[source][i].pc_ic[1] = recv_wvinfo[source][i].pc_ic[1]; 
            }  
        } 
        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                if(-1 == redistr_table[source][i].wv_id) continue; 
                mylevel = redistr_table[source][i].wv_level;  
                redistr_table[source][i].off_set[0] = 
                    grid_off_set_of_pc(recv_wvinfo,redistr_table[source][i].pc_ic,
                       mylevel,0,max_n_patch); 
                redistr_table[source][i].off_set[1] = 
                    grid_off_set_of_pc(recv_wvinfo,redistr_table[source][i].pc_ic,
                       mylevel,1,max_n_patch); 
            }  
        } 
        if(debugging("allgather_patch_info"))
        {
            for(source = 0; source < numnodes; source++)
            {
                for(i = 0; i < max_n_patch; i++)
                {
                    if(-1 == redistr_table[source][i].wv_id) continue; 
                    printf("source[%d] patch[%d] pc_ic[%d, %d] level[%d]"
                            "===> base[%d, %d] bound[%d, %d] off_set[%d, %d]\n",
                       source, i,  
                       redistr_table[source][i].pc_ic[0], 
                       redistr_table[source][i].pc_ic[1], 
                       redistr_table[source][i].wv_level,  
                       redistr_table[source][i].base[0], 
                       redistr_table[source][i].base[1],  
                       redistr_table[source][i].bound[0], 
                       redistr_table[source][i].bound[1],
                       redistr_table[source][i].off_set[0], 
                       redistr_table[source][i].off_set[1]);  
                }  
                printf("\n");  
            }
        } 
        free_these(2, send_storage, recv_storage);  
        free(recv_wvinfo); free(send_wvinfo);  
        DEBUG_LEAVE(allgather_patch_info)
}

LOCAL int grid_off_set_of_pc(
        Wv_info     **recv_wvinfo, 
        int         *my_ic,
        int         my_level,  
        int         dir,   
	int         max_n_patch)
{
	int        source, i; 
        int        numnodes; 
        int        total = 0; 
        int        *gmax; /* gmax in each previous pc_ic */ 
        int        increat = 0;  
     
        numnodes = pp_numnodes(); 

        if(my_ic[dir] == 0)     
            return 0; 
        uni_array(&gmax, my_ic[dir], sizeof(int));          

        for(source = 0; source < numnodes; source++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                if(my_level != recv_wvinfo[source][i].wv_level)  
                    continue; 
                if(recv_wvinfo[source][i].pc_ic[dir] < my_ic[dir])
                {
                    gmax[increat] = recv_wvinfo[source][i].gmax[dir]; 
                    increat++;  
                    break; 
                }  
            }  
            if(increat == my_ic[dir])
                break;  
        } 
        for(i = 0; i < my_ic[dir]; i++)
            total += gmax[i];
        free(gmax); 
        return total; 
}


/* NOT USABLE NOW  */  
/* return the current send size in send_storage */ 
/*  send has to be linear storage               */
LOCAL size_t bundle_send_info(
	byte        *storage, 
        POINTER     send,  
        int         ele_n,
        size_t      ele_size)
{
        int         i;  
       	byte       *buf, *send_head;
        size_t     send_len;  
	POINTER    info;  

        send_head = (byte*)send;  
	buf = storage;  

        for(i = 0; i < ele_n; i++)
        {
            info = (POINTER) buf;       
	    ft_assign(info,send_head,ele_size);   	    
	    buf += ele_size;  
	    send_head += ele_size;  
        }
        return ele_n*ele_size;   
}

/* return the receive size in recv_storage */ 
/*  recv has to be linear storage          */
LOCAL size_t unbundle_recv_info(
	byte        *storage, 
        POINTER     recv,  
        int         ele_n,
        size_t      ele_size)
{
        int         i;  
       	byte       *buf, *recv_head;
	POINTER    info;  

        recv_head = (byte*)recv;  
	buf = storage;  

        for(i = 0; i < ele_n; i++)
        {
            info = (POINTER) buf;       
	    ft_assign(recv_head,info,ele_size);   	    
	    buf += ele_size;  
	    recv_head += ele_size;  
        }
        return ele_size*ele_n;   
}

LOCAL void allgather_load_info(
        int         *total_load,
        int         *send_load,
        int         *proc_patch_n, 
        int         num_patches,  
	int         numnodes,
	int         max_n_patch)
{
	int        len, send_len; 
	int        i;  
        byte       *send_storage = NULL;
        byte       *recv_storage = NULL;
       	byte       *buf;
	POINTER    info;  
        int        myid; 

        myid = pp_mynode();  

	len = numnodes*max_n_patch*sizeof(int); 
	scalar(&recv_storage, len); 
	send_len = max_n_patch*sizeof(int); 
	scalar(&send_storage, send_len); 

	buf = send_storage;  
        for(i = 0; i < max_n_patch; i++)
        {
            info = (POINTER) buf;       
	    ft_assign(info,&send_load[i],sizeof(int));   	    
	    buf += sizeof(int);  
        }
        
        pp_all_gather(send_storage,send_len,recv_storage,send_len);  

        buf = recv_storage;  
        for(i = 0; i < max_n_patch*numnodes; i++)
        {
	    info = (POINTER) buf;       
	    ft_assign(&total_load[i],info,sizeof(int));   	    
	    buf += sizeof(int);  
        }
	free(send_storage);  
	free(recv_storage);  

        /* gather proc patch number */ 
	len = numnodes*sizeof(int); 
	scalar(&recv_storage, len); 
	send_len = sizeof(int); 
	scalar(&send_storage, send_len); 

	buf = send_storage;  
        info = (POINTER) buf;       
        ft_assign(info,&proc_patch_n[myid],sizeof(int));   	    
        
        pp_all_gather(send_storage,send_len,recv_storage,send_len);  

        buf = recv_storage;  
        for(i = 0; i < numnodes; i++)
        {
	    info = (POINTER) buf;       
	    ft_assign(&proc_patch_n[i],info,sizeof(int));   	    
	    buf += sizeof(int);  
        }
	free(send_storage);  
	free(recv_storage);  
}


LOCAL void bcast_balance(
        Wv_on_pc    **redistr_table,
	int         numnodes,
	int         max_n_patch)
{
	int        root = 0, len, myid; 
	int        i, source;  
        byte       *storage = NULL;
       	byte       *buf;
	POINTER    info;  

#if defined(__MPI__)	

	myid = pp_mynode();
	len = numnodes*max_n_patch*sizeof(Wv_on_pc); 
	scalar(&storage, len); 

	if(root == myid)
	{	
	    buf = storage;  
            for(source = 0; source < numnodes; source++)
            {
                for(i = 0; i < max_n_patch; i++)
                {
		    info = (POINTER) buf;       
	            ft_assign(info,&redistr_table[source][i],sizeof(Wv_on_pc));   	    
		    buf += sizeof(Wv_on_pc);  
                }
            }
        }

	pp_bcast(root,(POINTER)storage,len);  

	if(root != myid)
	{	
            buf = storage;  
            for(source = 0; source < numnodes; source++)
            {
                for(i = 0; i < max_n_patch; i++)
                {
		    info = (POINTER) buf;       
	            ft_assign(&redistr_table[source][i],info,sizeof(Wv_on_pc));   	    
		    buf += sizeof(Wv_on_pc);  
                }
            }
	}
	free(storage);  
#endif /* defined(__MPI__) */  	
}


LOCAL Wv_on_pc **redistr_load_as(
	int         **total_load,
        int         *proc_patch_n,
	int         max_n_patch,  
        int         numnodes)
{
        int         id, i, j, nxt;  
        int         sum_load = 0;  
        int         sum_patches = 0; 
        int         load_per_patch, load_per_proc;  /* averaged */
	int         *proc_patches;   /* number of patches in each proc */  
        int         *load_of_proc;   /* the total load on each proc */ 
	Wv_on_pc    **redistr;      /* bi_array, 
					each row (proc) stores the patch id */  
        int         tmp1, tmp2;            
	 

        DEBUG_ENTER(redistr_load_as)

        uni_array(&load_of_proc, numnodes, sizeof(int)); 
        uni_array(&proc_patches, numnodes, sizeof(int)); 

        for(id = 0; id < numnodes; id++)
        {
            load_of_proc[id] = 0;  
            for(i=0; i < proc_patch_n[id]; i++)
            {
		if(debugging("redistr_computing_load"))    
                    printf("patch[%d] from source[%d]'s load = %d\n",
                        i, id, total_load[id][i]);
                load_of_proc[id] += total_load[id][i];  
            } 
            printf("load_of_proc[%d] = %d\n", id, load_of_proc[id]);  

            sum_load += load_of_proc[id];  
            sum_patches += proc_patch_n[id]; 

	    proc_patches[id] = proc_patch_n[id];  
        } 
        load_per_patch = (int)sum_load/sum_patches;
        load_per_proc = (int)sum_load/numnodes;
        
	if(max_n_patch <= 0)
	{
	    printf("ERROR redistr_load_as()\n");
            printf("Max # of patches in proc. is less than 0, max_num = %d\n", 
		    max_n_patch);
            clean_up(ERROR);  	    
	}	
        

        bi_array(&redistr, numnodes, max_n_patch, sizeof(Wv_on_pc));  

        for(id = 0; id < numnodes; id++)
        {
            for(i = 0; i < max_n_patch; i++)
	    {	    
	        redistr[id][i].pc_id = id;  	    
	        if(proc_patch_n[id] > i)     
	            redistr[id][i].wv_id = i;  	    
                else
		    redistr[id][i].wv_id = -1;  	 
	    } 	
	}   	
	
        load_per_patch = 0;          

	if(debugging("redistr_computing_load"))
	{	
	    printf("_____________________________\n");  
	    printf("BEFORE REDISTRIBUTE, the patch is scattered as:\n");

            for(id = 0; id < numnodes; id++)
	    {
	        printf("Proc[%d]: ", id);  	
	        for(i = 0; i < max_n_patch; i++)
	        {	    
	            printf(" %2d %2d, ", redistr[id][i].pc_id, redistr[id][i].wv_id);  	    
	        } 	
	        printf("\n");    
	    }  	
            printf("load_per_patch = %d, load_per_proc = %d\n",
                      load_per_patch, load_per_proc); 
        }   

        for(id = 0; id < numnodes; id++)
        {
            if(load_of_proc[id] < (load_per_proc + load_per_patch)) continue; 

            for(nxt = 0; nxt < numnodes; nxt++)
            {
	        if(nxt == id) continue;  	
	        if(load_of_proc[nxt] > (load_per_proc + load_per_patch)) continue; 
                if(load_of_proc[id] < (load_per_proc + load_per_patch)) break; 

                for(i=proc_patch_n[id]-1; i >=0; i--)
                {
                    tmp1 = load_of_proc[id] - total_load[id][i];
                    tmp2 = load_of_proc[nxt] + total_load[id][i];

                    if(tmp2 < (load_per_proc + load_per_patch))
                    {
                        /* 042603 added */ 
                        if(proc_patches[id] > 1)
                        {
                            load_of_proc[id] -= total_load[id][i]; 
			    redistr[id][proc_patches[id]-1].pc_id = nxt;  

                            load_of_proc[nxt] += total_load[id][i]; 
                            proc_patches[id]--; 
			    proc_patches[nxt]++;
                        }  
                    }  
                    else
                        break;  
                    if(load_of_proc[id] < (load_per_proc+load_per_patch)) break;  
                    if(load_of_proc[nxt] > (load_per_proc+load_per_patch)) break;  
                }  
            }  
        } 
      
	if(debugging("redistr_computing_load"))
	{	
            for(id = 0; id < numnodes; id++)
            {
                printf("Balanced load of proc[%d] = %d, # of patches = %d\n",
                  id, load_of_proc[id], proc_patches[id]);  
            }
	}    

        /* make sure base patch is not transfered to the other procs */ 
        for(id = 0; id < numnodes; id++)
	    redistr[id][0].pc_id = id;  	    

        free(load_of_proc);  
        free(proc_patches);  

        DEBUG_LEAVE(redistr_load_as)
	return redistr;  	
}

LOCAL void calculate_computing_load(
	Wave        **wvs,
        Front       **frs,
        int         num_patches,
        int         *send_load)
{
	RECT_GRID   *rgr;  
        int         dim = wvs[0]->rect_grid->dim;
	int         i, j;
        int         load; 
        int         NumberOfLevels; 

        DEBUG_ENTER(calculate_computing_load)

        NumberOfLevels = frs[0]->NumberOfLevels; 
        for(i = 0; i < num_patches; i++)
        {
            rgr = wvs[i]->rect_grid;      
            load = 1;  
            for (j = 0; j < dim; j++)
            {
                load *= (rgr->gmax[j]+rgr->lbuf[j]+rgr->ubuf[j]);  
            }
            if(frs[i]->interf != NULL && 
               (frs[i]->patch_level == 0 || 
                frs[i]->patch_level == NumberOfLevels-1))
            {
                if(frs[i]->interf->curves != NULL) 
                    load += 15*(frs[i]->interf->num_points); 
            }
	    if(debugging("redistr_computing_load"))
                   printf("Patch[%d]'s load = %d, X[%d], Y[%d]\n",
                    i, load, rgr->gmax[0]+rgr->lbuf[0]+rgr->ubuf[0],
                     rgr->gmax[1]+rgr->lbuf[1]+rgr->ubuf[1]);  
            send_load[i] = load;  
        }   
        DEBUG_LEAVE(calculate_computing_load)
} 

LOCAL   CHART   *copy_chart(
        CHART       *ochart)
{
        CHART   *nchart;
        Grid    *ngrid;

        scalar(&nchart,sizeof(CHART));
        scalar(&ngrid,sizeof(Grid));

        *nchart = *ochart;
        *(ngrid) = *(ochart->grid);
        nchart->grid = ngrid;
        return nchart;
}

LOCAL void free_chart(
        CHART       *chart)
{
        free(chart->grid);
	free(chart);  
}

LOCAL void set_bdry_patch_flag(
        Wave            *wave,
        Wave            *wave0,
        Overparam       *overparam)
{
        Patch_bdry_flag *pd_flag;
        RECT_GRID       *rg0 = wave0->rect_grid;
        int             dim = rg0->dim;
        int             amr_ratio = overparam->refinementRatio;
        int             patch_level = wave->patch_level;
        int             i, j;

        CompositeGrid   *cg1 = (CompositeGrid *)wave->cg_over;
        int             patch_num = wave->patch_number;
        int             imin[3],imax[3];
        int             fine_gmax[3];

        RECT_GRID *tmpgr = wave->rect_grid; 
        pd_flag = wave->pd_flag;  

        for ( i = 0; i < dim; i++)
        {
            fine_gmax[i] = rg0->gmax[i];
        }

        if ( patch_level > 0)
        {
            for ( i = 1; i <= patch_level; i++ )
            {
                for ( j = 0; j < dim; j++)
                {
                    fine_gmax[j] = fine_gmax[j]*amr_ratio;
                }
            }
        }
        pd_flag->patch_bdry = 0;

        Range all;
        Index I1,I2,I3;
        getIndex((*cg1)[patch_num].indexRange(),I1,I2,I3);
        imin[0] = I1.getBase();
        imin[1] = I2.getBase();
        imin[2] = I3.getBase();
        imax[0] = I1.getBound();
        imax[1] = I2.getBound();
        imax[2] = I3.getBound();

        for ( i = 0; i < dim; i++)
        {
            pd_flag->bdry[i] = 0;
            pd_flag->bdry_side[i][0] = 0;
            pd_flag->bdry_side[i][1] = 0;

            if ( imin[i] == 0 )
            {
                 pd_flag->patch_bdry = 1;
                 pd_flag->bdry[i] = 1;
                 pd_flag->bdry_side[i][0] = 1;
            }
            if ( imax[i] == fine_gmax[i] )
            {
                 pd_flag->patch_bdry = 1;
                 pd_flag->bdry[i] = 1;
                 pd_flag->bdry_side[i][1] = 1;
            }
        }
        if(debugging("set_bdry_patch_flag"))
        {
            printf("\nPrinting boundary patch flag of the patch[%d]\n",
                    patch_num);
            printf("patch_level = %d\n",wave->patch_level);
            printf("patch_bdry %d\n", wave->pd_flag->patch_bdry);
            for (i = 0; i < dim; i++)
            {
                printf("%d-dim %d\n",i,wave->pd_flag->bdry[i]);
                printf("%d-dim side %d %d\n",
                  i, wave->pd_flag->bdry_side[i][0],
                  wave->pd_flag->bdry_side[i][1]);
            }
        }
}


LOCAL void set_patch_rect_grid(
        RECT_GRID     *rgr,
        RECT_GRID     *brgr,
        Overparam     *overparam,  
        CompositeGrid *cg_over,
        int           level, 
        int           grid)
{
        double           *L = rgr->L;
        double           *U = rgr->U;
        double           *GL = rgr->GL;
        double           *GU = rgr->GU;
        int             *lbuf = rgr->lbuf;
        int             *ubuf = rgr->ubuf;
        int             *gmax = rgr->gmax;
        int             dim = rgr->dim;
        int             base1,base2,base3,bound1,bound2,bound3;
        int             i, j, refine;  
        int             refinegmax[MAXD];  

        DEBUG_ENTER(set_patch_rect_grid)

        Range all;
        Index I1,I2,I3;
        getIndex((*cg_over)[grid].indexRange(),I1,I2,I3);
        base1 = I1.getBase();
        base2 = I2.getBase();
        base3 = I3.getBase();
        bound1 = I1.getBound();
        bound2 = I2.getBound();
        bound3 = I3.getBound();

        refine = 1;
        for(i = 0; i < level; i++)
            refine *= overparam->refinementRatio;
        for(i = 0; i < dim; i++)
            refinegmax[i] = refine*brgr->gmax[i];
        switch (dim)
        {
        case 1:
            L[0] = (*cg_over)[grid].vertex()(base1,0,0,axis1);
            U[0] = (*cg_over)[grid].vertex()(bound1,0,0,axis1);
            gmax[0] = bound1 - base1;
            GL[0] = L[0];
            GU[0] = U[0];
            lbuf[0] = 4;
            ubuf[0] = 4;
            break;
        case 2:
            L[0] = (*cg_over)[grid].vertex()(base1,base2,0,axis1);
            L[1] = (*cg_over)[grid].vertex()(base1,base2,0,axis2);
            U[0] = (*cg_over)[grid].vertex()(bound1,bound2,0,axis1);
            U[1] = (*cg_over)[grid].vertex()(bound1,bound2,0,axis2);
            gmax[0] = bound1 - base1;
            gmax[1] = bound2 - base2;
            GL[0] = L[0];
            GL[1] = L[1];
            GU[0] = U[0];
            GU[1] = U[1];
            /* 
            printf("cg_over [%d] refine grid[%d], level[%d] gmax[%d, %d]\n", 
                     cg_over, grid, level, gmax[0], gmax[1]);  
            printf("dx = %f, dy = %f\n", (U[0]-L[0])/gmax[0], (U[1]-L[1])/gmax[1]);
            printf("L[%f, %f], U[%f, %f]\n", L[0], L[1], U[0], U[1]);
            */  
            if((base1 < 0) or (base2 < 0)
               or (bound1 > refinegmax[0])
               or (bound2 > refinegmax[1]))
            {
                printf("ERROR: set_patch_rect_grid()\n");
                printf("Grid[%d] level[%d] base is less than 0\n",
                       grid, level);
                clean_up(ERROR);
            }
           if((base1 > 0) and (base1 < 4))
                lbuf[0] = base1;
            else if(0 == base1)
                lbuf[0] = brgr->lbuf[0];
            else
                lbuf[0] = 4;
            if((base2 > 0) and (base2 < 4))
                lbuf[1] = base2;
            else if(0 == base2)
                lbuf[1] = brgr->lbuf[1];
            else
                lbuf[1] = 4;

            if((refinegmax[0] - bound1) < 4 and (refinegmax[0] - bound1) > 0)
                ubuf[0] = refinegmax[0] - bound1;
            else if((refinegmax[0] - bound1) == 0)
                ubuf[0] = brgr->ubuf[0];
            else
                ubuf[0] = 4;
            if((refinegmax[1] - bound2) < 4 and (refinegmax[1] - bound2) > 0)
                ubuf[1] = refinegmax[1] - bound2;
            else if((refinegmax[1] - bound2) == 0)
                ubuf[1] = brgr->ubuf[1];
            else
                ubuf[1] = 4;
            break;
        case 3:
            L[0] = (*cg_over)[grid].vertex()(base1,base2,base3,axis1);
            L[1] = (*cg_over)[grid].vertex()(base1,base2,base3,axis2);
            L[2] = (*cg_over)[grid].vertex()(base1,base2,base3,axis3);
            U[0] = (*cg_over)[grid].vertex()(bound1,bound2,bound3,axis1);
            U[1] = (*cg_over)[grid].vertex()(bound1,bound2,bound3,axis2);
            U[2] = (*cg_over)[grid].vertex()(bound1,bound2,bound3,axis3);
            gmax[0] = bound1 - base1;
            gmax[1] = bound2 - base2;
            gmax[2] = bound3 - base3;
            GL[0] = L[0];
            GL[1] = L[1];
            GL[2] = L[2];
            GU[0] = U[0];
            GU[1] = U[1];
            GU[2] = U[2];
            lbuf[0] = 4;
            lbuf[1] = 4;
            lbuf[2] = 4;
            ubuf[0] = 4;
            ubuf[1] = 4;
            ubuf[2] = 4;
            break;
        }
        set_rect_grid(L,U,GL,GU,lbuf,ubuf,gmax,dim,&rgr->Remap,rgr);

        DEBUG_LEAVE(set_patch_rect_grid)
}

/*
 *
 *       This is a function to create the component in all patches,
 *            which is for those patches without interface
 */

LOCAL void set_patch_components(
        Wave            *wave,
        CompositeGrid   *cg1,
        int             pmax,
        COMPONENT       *patch_comp)
{
        int             patch_num;
        int             dim = wave->rect_grid->dim;
        int             base1,base2,base3,bound1,bound2,bound3;

        DEBUG_ENTER(set_patch_components)

        Range all;
        Index I1,I2,I3;

        for ( patch_num = 0; patch_num < pmax; patch_num++)
        {
            getIndex((*cg1)[patch_num].indexRange(),I1,I2,I3);
            base1 = I1.getBase();
            base2 = I2.getBase();
            base3 = I3.getBase();
            bound1 = I1.getBound();
            bound2 = I2.getBound();
            bound3 = I3.getBound();
            switch (dim)
            {
            case 1:
                break;
            case 2:
                int     pix, piy, ic[MAXD];
                double   crds[MAXD];

                pix = (base1 + bound1)/2;
                piy = (base2 + bound2)/2;
                crds[0] = (*cg1)[patch_num].vertex()(pix,piy,0,axis1);
                crds[1] = (*cg1)[patch_num].vertex()(pix,piy,0,axis2);

                if(rect_in_which(crds,ic,wave->rect_grid) == FUNCTION_FAILED)
                {
                    printf("ERROR in set_patch_components()\n");
                    printf("rect_in_which() failed \n");
                    clean_up(ERROR);
                }
                patch_comp[patch_num] = Rect_comp(ic, wave);
                break;
            case 3:
                break;
            }
        }

        DEBUG_LEAVE(set_patch_components)
}       /* end set_patch_components */


LOCAL int find_same_comp_ic_on_base(
        Wave            *wave,
        COMPONENT       comp,
        int             *icoords)
{
        int             i, dim; 
        int             smax[MAXD], smin[MAXD], ic[MAXD];
        int             iy, ix;  
        RECT_GRID       *gr;  

        gr = wave->rect_grid;
        dim = gr->dim; 
        for(i = 0; i < dim; i++)
        {
            smin[i] = -gr->lbuf[i];
            smax[i] = gr->gmax[i]+gr->ubuf[i];
        }
        for(iy = smin[1]; iy < smax[1]; iy++)
        {
            ic[1] = iy;  
            for(ix = smin[0]; ix < smax[0]; ix++)
            {
                ic[0] = ix;  
                if(Rect_comp(ic,wave) == comp)
                {
                    for(i = 0; i < dim; i++)
                        icoords[i] = ic[i];
                    return YES;  
                }  
            }
        }  
        return NO;  
} 
 
LOCAL CHART **perform_amr_restart_initialization(
        INIT_DATA       *init,
        INIT_PHYSICS    *ip,
        int             *patch_num,
        Wv_on_pc        ***redistr_table,
        int             *max_n_patch)
{
        CHART           **chart, *root;
        Front           *front,**frs, **nfrs;
        Wave            *wave,**wvs, **nwvs;
        Overparam       *overparam;
        int             i, j, dim, levels;
        CompositeGrid   *cg_over; 
        doubleCompositeGridFunction  *cg_over_function; 
        int             *level_array;
        int             num_patches, bal_n;   
        int             status = GOOD_STEP;
        boolean            stat = FUNCTION_SUCCEEDED; 
        RECT_GRID       *gr;
        int             *gmax;
        int             iy, ix;
        Index           I1,I2,I3;
        int             smax[MAXD], smin[MAXD];
        int             base[MAXD], bound[MAXD];
        int             ic[MAXD], baseic[MAXD];
        Locstate        st, basest;  
        COMPONENT       basecomp, comp; 
        size_t          sizest; 

        DEBUG_ENTER(perform_amr_restart_initialization)

        root = ip->root;
        overparam = root->overparam;
        wave = root->wave;
        sizest = wave->sizest;  
        front = root->front;
        dim = front->rect_grid->dim;
        levels = ip->root->overparam->numberOfRefinementLevels;
        wave->overture_init_step = NO;

        read_overture_restart_data(root,ip->prt);    

        cg_over = (CompositeGrid*)wave->cg_over;  
        cg_over_function = (doubleCompositeGridFunction*)wave->cg_over_function;

        num_patches = wave->totalNumberOfPatches = front->totalNumberOfPatches
            = root->totalNumberOfPatches = cg_over->numberOfComponentGrids();

        uni_array(&level_array,levels,sizeof(int));  

        uni_array(&frs,num_patches,sizeof(Front*));
        uni_array(&wvs,num_patches,sizeof(Wave*));

        set_patch_wv_fr_before_redistr(num_patches,root,frs,wvs,level_array); 
        if(debugging("amr_restart_init"))
        {
            printf("______________________________\n");
            printf("THE AMR RESTART INIT FRONT[%d]\n", 0);
            print_interface(frs[0]->interf);
            show_intfc_states(frs[0]->interf);   
        }

        for (i = 1; i < num_patches; i++)
        {
            set_patch_front(front,frs[i],frs[i]->rect_grid,i);

            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,-(frs[i]->NumberOfLevels-1.0-frs[i]->patch_level)));            

            if ((stat = clip_patch_front(frs[i],NO)) == FUNCTION_FAILED)
            {
                (void) printf("ERROR in perform_amr_restart_initialization(), "
                                "clip_patch_front()[%d] failed \n",i);
                break;
            }
            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,frs[i]->NumberOfLevels-1.0-frs[i]->patch_level));

            if((status = init_hyp_solution_function(wvs[i],frs[i])) != GOOD_STEP) 
            {
                screen("ERROR in perform_amr_restart_initialization(), "
                   "init_hyp_solution_function() failed on patch[%d]\n",i);
                stat = FUNCTION_FAILED;
                break;
            }
            if(debugging("amr_restart_init"))
            {
                printf("______________________________\n");
                printf("THE AMR RESTART INIT FRONT[%d]\n", i);
                print_interface(frs[i]->interf);
                show_intfc_states(frs[i]->interf);   
            }
        }

        if (pp_min_status(stat) == FUNCTION_FAILED)
        {
            printf("ERROR perform_amr_initialization(),"
                       " failed after clip_patch_front\n");
            clean_up(ERROR);
        }
        for(i = 1; i < num_patches; i++)
        {
            int  icoords[3] = {0,0,0};  
            gr = wvs[i]->rect_grid;
            gmax = wvs[i]->rect_grid->gmax;
            for(int dim = 0; dim < gr->dim; dim++)
            {
                smin[dim] = -gr->lbuf[dim];
                smax[dim] = gr->gmax[dim]+gr->ubuf[dim];
            }
            getIndex((*cg_over)[i].indexRange(),I1,I2,I3);
            base[0]  = I1.getBase();
            base[1] = I2.getBase();
            bound[0] = I1.getBound();
            bound[1] = I2.getBound();
            for(iy = smin[1]; iy < smax[1]; iy++)
            {
                for(ix = smin[0]; ix < smax[0]; ix++)
                {
                    ic[0] = ix; ic[1] = iy;
                    st = Rect_state(ic,wvs[i]);
                    if(NO == find_same_comp_ic_on_base(wvs[0],
                             Rect_comp(ic,wvs[i]),baseic))
                    {
                        printf("ERROR: in perform_amr_restart_initialization\n");
                        printf("No same component is found on the base wave comp = %d\n",
                           Rect_comp(ic,wvs[i]));
                        print_interface(frs[0]->interf);  
                        clean_up(ERROR);  
                    } 
                    ic[0] += base[0]; ic[1] += base[1];  
                    (*wvs[i]->overture_to_ft_st)(st,
                      (POINTER)cg_over_function,i,ic);  
                    (*wvs[i]->overture_assign_wave_params)(st,
                           Rect_state(baseic,wvs[0]));
                    (*wvs[i]->overture_assign_wave_st_type)(st,
                           Rect_state(baseic,wvs[0]));
                }
            }
        }

        *redistr_table = redistr_computing_load(wvs,frs,num_patches,max_n_patch);
        allgather_patch_info(*redistr_table,frs,overparam,num_patches,*max_n_patch);
        geomview_amr_fronts_plot2d(basename(ip->prt->outfile),frs,num_patches,
           *redistr_table,*max_n_patch); 
        bal_n = perform_redistribute_load(*redistr_table,wvs,frs,num_patches,
                      *max_n_patch,&nwvs,&nfrs); 

        set_current_interface(nfrs[0]->interf);  
        remove_sent_wv_fr_after_redistr(*redistr_table,num_patches,
                      *max_n_patch,frs,wvs);
        free_these(2,frs, wvs);
        free(level_array);
        uni_array(&chart,bal_n,sizeof(CHART*));
        nfrs[0]->interf->modified = YES; 
        for(i = 0; i < bal_n; i++)
        {
            chart[i] = copy_chart(root);
            chart[i]->grid->rect_grid = nfrs[i]->rect_grid;
            chart[i]->old_cg_over = nfrs[i]->cg_over;
            if(0 != i)
            {
                MaxWaveSpeed(nwvs[i]) = NULL; 
                MaxWaveSpeed(nwvs[i]) =
                 alloc_MaxWaveSpeed(MaxWaveSpeed(nwvs[i]),nwvs[i]);
                MaxFrontSpeed(nfrs[i]) = NULL;   
                MaxFrontSpeed(nfrs[i]) =
                 alloc_MaxFrontSpeed(MaxFrontSpeed(nfrs[i]),
                   nfrs[i]->interf,nfrs[i]->sizest);
            }
            chart[i]->wave = nwvs[i];
            chart[i]->front = nfrs[i];
            chart[i]->front->newfront = nfrs[i];
            chart_of_front(nfrs[i]) = chart[i];
        }

        free_these(2,nfrs,nwvs);
        if(NO == chart[0]->wave->use_overture_state)
        {
            cg_over_function->destroy();
            delete cg_over_function;   
            for(i = 0; i < bal_n; i++)
            {
                chart[i]->wave->cg_over_function = NULL; 
            }  
        } 
 
        *patch_num = bal_n; 
        printf("Report from vmalloc_storage_use: %-d M : after restart\n",
          get_vmalloc_storage_use()/1000000); 
        DEBUG_LEAVE(perform_amr_restart_initialization)
	return          chart;   
}  

LOCAL void read_overture_restart_data(
	CHART       *chart,
        Printplot   *prt)
{
        Wave            *wave = chart->wave;  
        CompositeGrid   *cg_over; 
        doubleCompositeGridFunction  *cg_over_function; 
        int         num_float;
        size_t      sizest = wave->sizest;
        char            outdir[256];
        int             myid, numnodes;
        const char      *nstep, *nd;
        char            outname[256], buf[256];
        char            *dname = basename(prt->outfile);
        HDF_DataBase    db;
        Range           all;  

        DEBUG_ENTER(read_overture_restart_data) 
 
        num_float = sizest/sizeof(double)+1;
        myid = pp_mynode(); numnodes = pp_numnodes();
        sprintf(outdir,"%s/%s",dname,"over_dump_st");

        nd = right_flush(pp_mynode(),4);
        sprintf(buf,"%s/%s_%s.%s",outdir,"over","dump_st",nd);
        /* step need - 1, because in set_up_cauchy_data,
         * step is augmented. Consider to move this function
         * into set_up_cauchy_data
         */  
        nstep = right_flush(chart->grid->step-1,7);
        sprintf(outname,"%s.ts%s",buf,nstep);

        db.mount(outname,"R");
        cg_over = new CompositeGrid();
        cg_over->get(db,"My Grid");

        cg_over_function = new doubleCompositeGridFunction();
        cg_over_function->updateToMatchGrid((*cg_over),num_float,all, all, all);
        cg_over_function->get(db,"My Solution");

        wave->cg_over_function = chart->cg_over_function = (POINTER)cg_over_function;
        chart->wave->cg_over = chart->front->cg_over 
             = chart->old_cg_over = (POINTER)cg_over;  
        wave->patch_number = chart->front->patch_number = 0;
        db.unmount();

        DEBUG_LEAVE(read_overture_restart_data) 
} 

LOCAL int stop_adding_amr_grids(
	CompositeGrid  cgc,
        CompositeGrid  cgNew)
{
        int   stop = NO;  
        int   i, start, num;  

        if(cgNew.numberOfComponentGrids() != cgc.numberOfComponentGrids())
        {        
            stop = NO;  
            return NO; 
        } 
        /* If the numbers are the same, should check 
         * for every two grids coming from two CompositeGrids, 
         * whether the levels  are matching.
         */
        for(i = 0; i < cgNew.numberOfComponentGrids(); i++)
        {
            if((cgc).refinementLevelNumber(i) != (cgNew).refinementLevelNumber(i))
            {
                stop = NO;
                return NO;  
            } 
        } 
        stop = YES;  
        return YES;  
}

LOCAL int overture_init_mesh_refinement(
        Overparam  *overparam,
        Wave       *wave,
        Front      *front,
        int        *patch_level)
{
        int     dim = front->rect_grid->dim;
        int     i, j;
        size_t  sizest = wave->sizest; 
        int     numberOfRefinementLevels = overparam->numberOfRefinementLevels;
        CompositeGrid   *rect_over_grid;
        Range           all; 
	Index           I1,I2,I3; 

        DEBUG_ENTER(overture_init_mesh_refinement)

        rect_over_grid = (CompositeGrid*)wave->cg_over;

	if(debugging("over_init_mesh"))
	{
            printf("--------------------------------------\n");  
            printf("Before start init_mesh_refine:\n");  
	    for(int grid = 0; grid<rect_over_grid->numberOfComponentGrids(); grid++)
	    {
	        getIndex((*rect_over_grid)[grid].indexRange(),I1,I2,I3);
                printf("grid = %d\n",grid);
                printf("x dir indexRange base, bond[%d,%d]\n",I1.getBase(), I1.getBound());
                printf("y dir indexRange base, bond[%d,%d]\n",I2.getBase(), I2.getBound());
	        getIndex((*rect_over_grid)[grid].dimension(),I1,I2,I3);
                printf("grid = %d\n",grid);
                printf("x dir dimension base, bond[%d,%d]\n",I1.getBase(), I1.getBound());
                printf("y dir dimension base, bond[%d,%d]\n",I2.getBase(), I2.getBound());
	    } 
            printf("--------------------------------------\n");  
	} 

        printf(" ------------------------------------------------------------------- \n");
        printf("  Init. Adaptive Mesh Refinement for FronTier Solver                 \n");
        printf(" ------------------------------------------------------------------- \n");

	CompositeGrid  *new_cg;
	doubleCompositeGridFunction Dens(*rect_over_grid); 

	(*wave->trans_wv_st_to_overfunc)(wave, front,(POINTER)(&Dens),0,0); 
	(*wave->fill_root_extr_overture_cell_st_from_wv)(wave, front,(POINTER)(&Dens),0,0); 

        perform_overture_init_mesh_refinement(rect_over_grid,
               &Dens,&new_cg,overparam, 0);
  
        front->cg_over = wave->cg_over = (POINTER)new_cg;

        if(YES == wave->use_overture_state)
        {
            doubleCompositeGridFunction *uu1 = new 
                doubleCompositeGridFunction(*new_cg,
                     sizest/sizeof(double)+1,all,all,all);
            wave->cg_over_function = (POINTER)uu1;
        }
        else
            wave->cg_over_function = NULL;  

	if(debugging("over_init_mesh"))
	{
            printf("--------------------------------------\n");  
            printf("AFTER init_mesh_refine:\n");  
	    for(int grid = 0; grid<new_cg->numberOfComponentGrids(); grid++)
	    {
	        getIndex((*new_cg)[grid].indexRange(),I1,I2,I3);
                printf("grid = %d\n",grid);
                printf("x dir indexRange base, bond[%d,%d]\n",I1.getBase(), I1.getBound());
                printf("y dir indexRange base, bond[%d,%d]\n",I2.getBase(), I2.getBound());
	        getIndex((*new_cg)[grid].dimension(),I1,I2,I3);
                printf("grid = %d\n",grid);
                printf("x dir dimension base, bond[%d,%d]\n",I1.getBase(), I1.getBound());
                printf("y dir dimension base, bond[%d,%d]\n",I2.getBase(), I2.getBound());
	    } 
            printf("--------------------------------------\n");  
	} 

        front->totalNumberOfPatches = wave->totalNumberOfPatches
	       	= new_cg->numberOfComponentGrids();

        printf(" ------------------------------------------------------------------- \n");
        printf("  NumberOfComponentGrids = %d, RefinementLevelNumber = %d\n",
                 (new_cg)->numberOfComponentGrids(),
                 (new_cg)->refinementLevelNumber((new_cg)->numberOfComponentGrids()-1)+1);
        for(j = 0; j < (new_cg)->numberOfComponentGrids(); j++)
            printf("  refine patch [%d], level = %d\n", j, (new_cg)->refinementLevelNumber(j));

        printf("  End Init. Adaptive Mesh Refinement for FronTier Solver             \n");
        printf(" ------------------------------------------------------------------- \n");

	DEBUG_LEAVE(overture_init_mesh_refinement)

        return new_cg->numberOfComponentGrids();
}  

LOCAL int perform_overture_init_mesh_refinement(
        CompositeGrid                  *rect_over_grid,
        doubleCompositeGridFunction    *cg_function,
        CompositeGrid                  **new_cg,
        Overparam   *overparam,
        int         start_level)  /* This level is already created */  
{
        int     i, j;
        int     base[MAXD], bound[MAXD]; 
        int     baseLevel = overparam->baseLevel; 
        int     numberOfRefinementLevels = overparam->numberOfRefinementLevels;
        real    efficiency = overparam->efficiency;
        real    errorThreshold = overparam->errorThreshold;   
        int     refinementRatio = overparam->refinementRatio;
        Regrid          regrid; 
	Ogen            ogen; 
        Range           all;
	Index           I1,I2,I3; 

        DEBUG_ENTER(perform_overture_init_mesh_refinement)

	if(debugging("over_init_mesh"))
	{
            printf("--------------------------------------\n");  
            printf("Before start init_mesh_refine:\n");  
	    for(int grid = 0; grid<rect_over_grid->numberOfComponentGrids(); grid++)
	    {
	        getIndex((*rect_over_grid)[grid].indexRange(),I1,I2,I3);
                printf("grid = %d\n",grid);
                printf("x dir indexRange base, bond[%d,%d]\n",I1.getBase(), I1.getBound());
                printf("y dir indexRange base, bond[%d,%d]\n",I2.getBase(), I2.getBound());
	        getIndex((*rect_over_grid)[grid].dimension(),I1,I2,I3);
                printf("grid = %d\n",grid);
                printf("x dir dimension base, bond[%d,%d]\n",I1.getBase(), I1.getBound());
                printf("y dir dimension base, bond[%d,%d]\n",I2.getBase(), I2.getBound());
	    } 
            printf("--------------------------------------\n");  
	} 

        CompositeGrid cga[2];
        cga[0] = *rect_over_grid; 
        cga[0].update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEmask );
        cga[1] = *rect_over_grid; 
        cga[1].update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEmask );

        InterpolateRefinements interp(cga[0].numberOfDimensions()); 
        interp.setOrderOfInterpolation(3);  /* was 1 */ 

        ErrorEstimator errorEstimator(interp);  
        /* set the scale factors for the solution (estimate of solution size) */ 
        RealArray scaleFactor(1);
        scaleFactor=1.;
        errorEstimator.setScaleFactor(scaleFactor);
        /* set the default number of smoothing steps for smoothing the error */
        errorEstimator.setDefaultNumberOfSmooths(overparam->numberOfSmooths);

        Interpolant interpolant;   
        interpolant.setImplicitInterpolationMethod(Interpolant::iterateToInterpolate);
        interpolant.setInterpolateRefinements(interp);

        regrid.setEfficiency(efficiency);
        regrid.setRefinementRatio(refinementRatio);

	CompositeGridOperators op(cga[0]); 

	doubleCompositeGridFunction error; 
        cg_function->updateToMatchGrid(cga[0]); 
	cg_function->setOperators(op); 

        for( int grid = 0; grid < cga[0].numberOfComponentGrids(); grid++)
        {
            /*  
            printf("Original grid[%d] print Dens\n", grid);
            (*cg_function)[grid].display("density");
            */  
        }

        int currentGrid = 0, nextGrid;
        realArray pos(1,3), uval(1,7); 

        for( int level=start_level+1; level<numberOfRefinementLevels; level++ )
        {
            nextGrid = (currentGrid+1) % 2;

            CompositeGrid& cgc = cga[currentGrid];
            CompositeGrid&  cgNew = cga[nextGrid];

            error.updateToMatchGrid(cgc);
            interpolant.updateToMatchGrid(cgc,level-1); 
            op.updateToMatchGrid(cgc);
            error.setOperators(op);
            printf("Error computed = %d ",
             errorEstimator.computeAndSmoothErrorFunction(*cg_function, error));  
            printf("Error: min=%e, max=%e \n",min(error),max(error));

            if(overparam->extra_buf != 0)
                update_error(&error, cgc, overparam);

            regrid.regrid(cgc,cgNew, error, errorThreshold, level, baseLevel);
            ogen.updateRefinement(cgNew);
            cgNew.update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEmask );
            
            printf("level[%d], cgNew # of levels [%d], numberOfComponentGrids [%d]\n",
                   level,cgNew.numberOfMultigridLevels(),cgNew.numberOfComponentGrids());
            printf("cgc # of levels [%d],numberOfComponentGrids [%d]\n",
                   cgc.numberOfMultigridLevels(), cgc.numberOfComponentGrids());            

            if(YES == stop_adding_amr_grids(cgc, cgNew))
	    {	    
                currentGrid = (currentGrid+1) % 2;
                show_grids_info(cgNew);
		break;  
	    } 	

            doubleCompositeGridFunction uu2;
            uu2.updateToMatchGrid(cgNew);
            interp.interpolateRefinements(*cg_function,uu2);

            for( int grid = 1; grid < cgNew.numberOfComponentGrids(); grid++)
            {
                int base1,base2,bound1,bound2;
                int ix,iy;
                double x0,x1,y0,y1,dx,dy;
                getIndex(cgNew[grid].indexRange(),I1,I2,I3);
                base1 = I1.getBase();
                base2 = I2.getBase();
                bound1 = I1.getBound();
                bound2 = I2.getBound();

                for (ix = base1-4; ix <= bound1+4; ix++)
                {
                    for (iy = base2-4; iy < base2; iy++)
                    {
                        pos(0,0) = cgNew[grid].vertex()(ix,iy,0,axis1);
                        pos(0,1) = cgNew[grid].vertex()(ix,iy,0,axis2);
                        pos(0,2) = 0.0;
                        interpolatePoints(pos, *cg_function, uval); 
                        uu2[grid](ix,iy,0) = uval(0,0);
                    }
                    for (iy = bound2; iy <= bound2+4; iy++)
                    {
                        pos(0,0) = cgNew[grid].vertex()(ix,iy,0,axis1);
                        pos(0,1) = cgNew[grid].vertex()(ix,iy,0,axis2);
                        pos(0,2) = 0.0;
                        interpolatePoints(pos, *cg_function, uval); 
                        uu2[grid](ix,iy,0) = uval(0,0);
                    }
                }
                for (iy = base2-4; iy <= bound2+4; iy++)
                {
                    for (ix = base1-4; ix < base1; ix++)
                    {
                        pos(0,0) = cgNew[grid].vertex()(ix,iy,0,axis1);
                        pos(0,1) = cgNew[grid].vertex()(ix,iy,0,axis2);
                        pos(0,2) = 0.0;
                        interpolatePoints(pos, *cg_function, uval); 
                        uu2[grid](ix,iy,0) = uval(0,0);
                    }
                    for (ix = bound1; ix <= bound1+4; ix++)
                    {
                        pos(0,0) = cgNew[grid].vertex()(ix,iy,0,axis1);
                        pos(0,1) = cgNew[grid].vertex()(ix,iy,0,axis2);
                        pos(0,2) = 0.0;
                        interpolatePoints(pos, *cg_function, uval); 
                        uu2[grid](ix,iy,0) = uval(0,0);
                    }
                }
            }   
            (cg_function)->updateToMatchGrid(cgNew);
            (cg_function)->dataCopy(uu2);

            printf("Create Amr Grids\n");  
            show_grids_info(cgNew);

            currentGrid = (currentGrid+1) % 2;
        }

	(*new_cg) = new CompositeGrid();
        (**new_cg) = cga[currentGrid];

	DEBUG_LEAVE(perform_overture_init_mesh_refinement)

        return (*new_cg)->numberOfComponentGrids();
}  

LOCAL void show_grids_info(CompositeGrid& cg)
{
        int        imin[MAXD], imax[MAXD]; 
        double      LL[MAXD], UU[MAXD], dx, dy;    
        Index      I1,I2,I3; 
  
	for( int grid = 0; grid < cg.numberOfComponentGrids(); grid++)
        {
            getIndex(cg[grid].indexRange(),I1,I2,I3);

            imin[0] = I1.getBase();
            imin[1] = I2.getBase();
            imax[0] = I1.getBound();
            imax[1] = I2.getBound();
            LL[0] = (cg)[grid].vertex()(imin[0],imin[1],0,axis1);
            LL[1] = (cg)[grid].vertex()(imin[0],imin[1],0,axis2);
            UU[0] = (cg)[grid].vertex()(imax[0],imax[1],0,axis1);
            UU[1] = (cg)[grid].vertex()(imax[0],imax[1],0,axis2);
            dx = (UU[0] - LL[0]) / ((double)(imax[0] - imin[0]));
            dy = (UU[1] - LL[1]) / ((double)(imax[1] - imin[1]));

            printf("Grid[%d], x range<%d,%d>, y range<%d,%d>, dx,dy[%g,%g]\n",
                    grid,I1.getBase(),I1.getBound(), I2.getBase(), I2.getBound(),
                    dx,dy);
            printf("grid[%d], imin[%d,%d], imax[%d,%d] LL<%g,%g> UU<%g,%g>\n",
                grid,imin[0],imin[1],imax[0],imax[1],LL[0],LL[1],UU[0],UU[1]);
        }

}  

#endif /* defined(USE_OVERTURE) */







