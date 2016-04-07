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
* *                               doverturepatchmesh2.c:
* *
* *       Copyright 1999 by The University at Stony Brook, All rights reserved. * *
* *       Contains routines for initializing and maintaining nodes in the
* *       automatic mesh grid level tree.
* */

#define DEBUG_STRING    "doverturepatchmesh2"
#include <driver/ddecs.h>

#if defined(USE_OVERTURE)

#define normal_advance_front(dt,dt_frac,front,newfront,wave)                    \
        (*(front)->_normal_advance_front)(dt,dt_frac,front,newfront,wave)
#define tangnt_advance_front(dt,dt_frac,front,newfront,wave)                    \
        (*(front)->_tangnt_advance_front)(dt,dt_frac,front,newfront,wave)
#define redist_advance_front(dt,dt_frac,front,newfront,wave)                    \
        (*(front)->_redist_advance_front)(dt,dt_frac,front,newfront,wave)

LOCAL void copy_gridfunction_from_uu1(const doubleCompositeGridFunction&,
                    const CompositeGrid&, int, doubleCompositeGridFunction&);

LOCAL int  pseudo_advance_amr_fronts2d_old(Front***,double,double*,Overparam*,
                 Wave**,Front**,int,Wv_on_pc**,int);

EXPORT int pseudo_advance_amr_fronts2d(
        Front           ***newfronts, 
        double           dt,
        double           *dt_frac,
        Overparam       *overparam,
        Wave            **wvs,
        Front           **frs,
        int             num_patches,
        Wv_on_pc        **redistr_table,
        int             max_n_patch)
{
        int                             timesteps;
        int                             i, k;
        int                             status = GOOD_STEP;
        Front                           **newfrs;
        Front                           *redistr_fr, *n_redistr_fr;  
        RECT_GRID                       *rgr, *cgr, *tgr;  
        INTERFACE                       *tmp_intfc, *sav_intfc;
        int                             *iperm;
        int                             dim = frs[0]->interf->dim;
        int                             step = frs[0]->step;
        boolean                            sav_copy;
        Wave                            **undistr_wvs;
        Front                           **undistr_frs;
        int                             resid_n, refine;
        int                             myid, numnodes;

        DEBUG_ENTER(pseudo_advance_amr_fronts2d)
        myid = pp_mynode();
        numnodes = pp_numnodes();
        
        uni_array(&newfrs,num_patches,sizeof(Front*));
        /* Normal front sweep */
        status = GOOD_STEP;
        for (i = num_patches-1; i >=0; i--)
        {
            if(frs[i]->patch_level != frs[i]->NumberOfLevels-1
               && frs[i]->patch_level != 0)
            {
                sav_intfc = current_interface();
                sav_copy = copy_intfc_states();
                newfrs[i] = copy_front(frs[i]);
                newfrs[i]->interf = copy_interface(frs[i]->interf);
                set_current_interface(sav_intfc);
                set_copy_intfc_states(sav_copy);
                delete_patch_all_curves(newfrs[i]);
                newfront_to_distri_table(frs[i],newfrs[i],
                     num_patches,redistr_table,max_n_patch);
                continue;
            }
            if (debugging("hyp_amr"))
                printf("normal_advance_front in patch[%d], level[%d]",
                                i,frs[i]->patch_level);
            status = normal_advance_front(dt,dt_frac,frs[i],&newfrs[i],
                                (POINTER)wvs[i]);
            if (debugging("hyp_amr"))
                printf(" status = %d\n",status);
            switch(status)
            {
            case GOOD_STEP:
            break;

            case MODIFY_TIME_STEP:
            case REPEAT_TIME_STEP:
                goto normal_status;
            break;
            case ERROR_IN_STEP:
            default:
            screen("FT_ERROR: in hyp_amr(),"
                   "  normal_advance_front() failed\n");
            clean_up(ERROR);
            }
            newfront_to_distri_table(frs[i],newfrs[i],
                 num_patches,redistr_table,max_n_patch);
        }

normal_status:
        status = syncronize_time_step_status(status,frs[0]->pp_grid);
        if (status != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
               " failed after normal_advance_front, status = %d\n",status);
            clean_up(ERROR);
        }

        printf("After Normal: assembly_distribute_patch_fronts ");
        status = assembly_distribute_patch_fronts(newfrs,num_patches,
              redistr_table,max_n_patch, NO);
        printf("status = %d\n", status);          
        if (pp_min_status(status) != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
              " failed after normal: assembly_distribute_patch_fronts()\n");
            clean_up(ERROR);
        }

        /* Tangential front sweep */
        status = GOOD_STEP;
        for (i = num_patches-1; i >= 0; i--)
        {
            if(frs[i]->patch_level != frs[i]->NumberOfLevels-1
               && frs[i]->patch_level != 0)
                continue;

            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,-(frs[i]->NumberOfLevels-1.0-frs[i]->patch_level)));
            set_amr_intfc_tol(newfrs[i]->interf,
                pow(2.0,-(newfrs[i]->NumberOfLevels-1.0-newfrs[i]->patch_level)));

            if (debugging("hyp_amr"))
                printf("tangnt_advance_front in patch[%d],level[%d]",
                       i,frs[i]->patch_level);
            status = tangnt_advance_front(dt,dt_frac,frs[i],&newfrs[i],
                                (POINTER)wvs[i]);

            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,frs[i]->NumberOfLevels-1.0-frs[i]->patch_level));
            set_amr_intfc_tol(newfrs[i]->interf,
                pow(2.0,newfrs[i]->NumberOfLevels-1.0-newfrs[i]->patch_level));

            if (debugging("hyp_amr"))
                printf(" status = %d\n",status);
            switch(status)
            {
            case GOOD_STEP:
            break;

            case MODIFY_TIME_STEP:
            case REPEAT_TIME_STEP:
                goto tangnt_status;
            break;
            case ERROR_IN_STEP:
            default:
            screen("ERROR: overture_hyp_amr(),"
                   " tangnt_advance_front() failed\n");
            clean_up(ERROR);
            }
        }

tangnt_status:
        status = syncronize_time_step_status(status,frs[0]->pp_grid);
        if (status != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
                   " failed after tangnt_advance_front, status = %d\n", status);
            clean_up(ERROR);
        }
        printf("After Tangnt: assembly_distribute_patch_fronts ");
        status = assembly_distribute_patch_fronts(newfrs,num_patches,
                  redistr_table,max_n_patch, NO);
        printf("status = %d\n", status);          
        if (pp_min_status(status) != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
             " failed after tangnt: assembly_distribute_patch_fronts()\n");
            clean_up(ERROR);
        }
        
        status = GOOD_STEP;
        /* Redistribution front sweep */
        for (i = 0; i < num_patches; i++)
        {
            if(frs[i]->patch_level != frs[i]->NumberOfLevels-1
                && frs[i]->patch_level != 0)
            {
                newfrs[i]->step = frs[i]->step + 1;
                newfrs[i]->time = frs[i]->time + dt;
                continue;
            }
            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,-(frs[i]->NumberOfLevels-1.0-frs[i]->patch_level)));
            set_amr_intfc_tol(newfrs[i]->interf,
                pow(2.0,-(newfrs[i]->NumberOfLevels-1.0-newfrs[i]->patch_level)));

            if (debugging("hyp_amr"))
                printf("Doing redistribute advance patch front[%d],",i);
            status = redist_advance_front(dt,dt_frac,frs[i],&newfrs[i],
                                (POINTER)wvs[i]);

            set_amr_intfc_tol(frs[i]->interf,
                pow(2.0,frs[i]->NumberOfLevels-1.0-frs[i]->patch_level));
            set_amr_intfc_tol(newfrs[i]->interf,
                pow(2.0,newfrs[i]->NumberOfLevels-1.0-newfrs[i]->patch_level));

            if (debugging("hyp_amr"))
                printf("status = %d\n",status);
            switch(status)
            {
            case GOOD_STEP:
            break;

            case MODIFY_TIME_STEP:
            case REPEAT_TIME_STEP:
                printf("WARNING: Bad redist_advance_front(), status = %d\n",
                   status);
                goto redis_status;
            break;
            case ERROR_IN_STEP:
            default:
            screen("ERROR: overture_hyp_amr(),"
                   " redist_advance_front() failed\n");
            clean_up(ERROR);
            }
        }

redis_status:
        pp_gsync();
        status = syncronize_time_step_status(status,frs[0]->pp_grid);
        if (status != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
             " failed after redist_advance_front(), status = %d\n", status);
            clean_up(ERROR);
        }

        for (i = 1; i < num_patches; i++)
        { 
            Redistribution_count(frs[i]) =
                  Redistribution_count(frs[0]);
        } 
       
        printf("After Redist: assembly_distribute_patch_fronts ");
        status = assembly_distribute_patch_fronts(newfrs,num_patches,
              redistr_table,max_n_patch, YES);
        printf("status = %d\n", status);          
        if (pp_min_status(status) != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
              " failed after redistr: assembly_distribute_patch_fronts()\n");
            clean_up(ERROR);
        }

        *newfronts = newfrs;
        DEBUG_LEAVE(pseudo_advance_amr_fronts2d)
        return status;
}


LOCAL int pseudo_advance_amr_fronts2d_old(
        Front           ***newfronts, 
        double           dt,
        double           *dt_frac,
        Overparam       *overparam,
        Wave            **wvs,
        Front           **frs,
        int             num_patches,
        Wv_on_pc        **redistr_table,
        int             max_n_patch)
{
        int                             timesteps;
        int                             i, k;
        int                             status = GOOD_STEP;
        Front                           **newfrs;
        Front                           *redistr_fr, *n_redistr_fr;  
        RECT_GRID                       *rgr, *cgr, *tgr;  
        INTERFACE                       *tmp_intfc, *sav_intfc;
        int                             *iperm;
        int                             dim = frs[0]->interf->dim;
        int                             step = frs[0]->step;
        boolean                            sav_copy;
        Wave                            **undistr_wvs;
        Front                           **undistr_frs;
        int                             resid_n, refine;
        int                             myid, numnodes;

        DEBUG_ENTER(pseudo_advance_amr_fronts2d)
        myid = pp_mynode();
        numnodes = pp_numnodes();
        
        uni_array(&newfrs,num_patches,sizeof(Front*));
        /* Normal front sweep */
        status = GOOD_STEP;
        for (i = num_patches-1; i >=0; i--)
        {
            if(frs[i]->patch_level != frs[i]->NumberOfLevels-1
               && frs[i]->patch_level != 0)
            {
                sav_intfc = current_interface();
                sav_copy = copy_intfc_states();
                newfrs[i] = copy_front(frs[i]);
                newfrs[i]->interf = copy_interface(frs[i]->interf);
                set_current_interface(sav_intfc);
                set_copy_intfc_states(sav_copy);
                delete_patch_all_curves(newfrs[i]);
                newfront_to_distri_table(frs[i],newfrs[i],
                     num_patches,redistr_table,max_n_patch);
                continue;
            }
            if (debugging("hyp_amr"))
                printf("normal_advance_front in patch[%d], level[%d]",
                                i,frs[i]->patch_level);
            status = normal_advance_front(dt,dt_frac,frs[i],&newfrs[i],
                                (POINTER)wvs[i]);
            if (debugging("hyp_amr"))
                printf(" status = %d\n",status);
            switch(status)
            {
            case GOOD_STEP:
            break;

            case MODIFY_TIME_STEP:
            case REPEAT_TIME_STEP:
                goto normal_status;
            break;
            case ERROR_IN_STEP:
            default:
            screen("FT_ERROR: in hyp_amr(),"
                   "  normal_advance_front() failed\n");
            clean_up(ERROR);
            }
            newfront_to_distri_table(frs[i],newfrs[i],
                 num_patches,redistr_table,max_n_patch);
        }

normal_status:
        status = syncronize_time_step_status(status,frs[0]->pp_grid);
        if (status != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
               " failed after normal_advance_front, status = %d\n",status);
            clean_up(ERROR);
        }

        for(i = 0; i < num_patches; i++)
        {
            /* 
            printf("After Normal prop, new interface[%d] address [%p]\n",
               i, newfrs[i]->interf);
            for(CURVE **c = newfrs[i]->interf->curves; c && *c; c++)
            {
                printf("curve %p, wave_type [%s]\n", *c, 
                  wave_type_as_string(wave_type(*c),newfrs[i]->interf)); 
            } 
            {
                for(CURVE **c = newfrs[i]->interf->curves; c && *c; c++)
                {
                    if(wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
                    {
                        printf("curve %p, wave_type [%s]\n", *c, 
                          wave_type_as_string(wave_type(*c),newfrs[i]->interf)); 
                        show_curve_states(*c);  
                    }
                } 
            }
            */  
        }

        printf("After Normal: assembly_distribute_patch_fronts ");
        status = assembly_distribute_patch_fronts(newfrs,num_patches,
              redistr_table,max_n_patch, NO);
        printf("status = %d\n", status);          
 
        if (pp_min_status(status) != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
              " failed after normal: assembly_distribute_patch_fronts()\n");
            clean_up(ERROR);
        }

        /* Tangential front sweep */
        status = GOOD_STEP;
        for (i = num_patches-1; i >= 0; i--)
        {
            if(frs[i]->patch_level != frs[i]->NumberOfLevels-1
               && frs[i]->patch_level != 0)
                continue;
            if (debugging("hyp_amr"))
                printf("tangnt_advance_front in patch[%d],level[%d]",
                       i,frs[i]->patch_level);
            status = tangnt_advance_front(dt,dt_frac,frs[i],&newfrs[i],
                                (POINTER)wvs[i]);
            if (debugging("hyp_amr"))
                printf(" status = %d\n",status);
            switch(status)
            {
            case GOOD_STEP:
            break;

            case MODIFY_TIME_STEP:
            case REPEAT_TIME_STEP:
                goto tangnt_status;
            break;
            case ERROR_IN_STEP:
            default:
            screen("ERROR: overture_hyp_amr(),"
                   " tangnt_advance_front() failed\n");
            clean_up(ERROR);
            }
        }

tangnt_status:
        status = syncronize_time_step_status(status,frs[0]->pp_grid);
        if (status != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
                   " failed after tangnt_advance_front, status = %d\n", status);
            clean_up(ERROR);
        }

        printf("After Tangnt: assembly_distribute_patch_fronts ");
        status = assembly_distribute_patch_fronts(newfrs,num_patches,
                  redistr_table,max_n_patch, NO);
        printf("status = %d\n", status);          
        
        if (pp_min_status(status) != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
             " failed after tangnt: assembly_distribute_patch_fronts()\n");
            clean_up(ERROR);
        }

        /*  
        redistr_fr = deep_copy_front(frs[0]); 
        n_redistr_fr = deep_copy_front(newfrs[0]); 
        refine = 1;
        for(i = 0; i < (overparam->numberOfRefinementLevels-1); i++)
            refine *= overparam->refinementRatio;

        rgr = redistr_fr->rect_grid;
        for(i = 0; i < dim; i++)
        {
            rgr->gmax[i] *= refine;
            rgr->lbuf[i] *= refine;
            rgr->ubuf[i] *= refine;
            rgr->h[i] /= (1.0*refine);
        }
        cgr = computational_grid(redistr_fr->interf);
        copy_rect_grid(cgr,rgr);
        tgr = &topological_grid(redistr_fr->interf);
        tgr->Remap.remap = rgr->Remap.remap;
        set_patch_topo_grid(rgr,tgr);        

        rgr = n_redistr_fr->rect_grid;
        for(i = 0; i < dim; i++)
        {
            rgr->gmax[i] *= refine;
            rgr->lbuf[i] *= refine;
            rgr->ubuf[i] *= refine;
            rgr->h[i] /= (1.0*refine);
        }
        cgr = computational_grid(n_redistr_fr->interf);
        copy_rect_grid(cgr,rgr);
        tgr = &topological_grid(n_redistr_fr->interf);
        tgr->Remap.remap = rgr->Remap.remap;
        set_patch_topo_grid(rgr,tgr);        
        */

        status = GOOD_STEP;
        /* Redistribution front sweep */
        /* 061503, redistribution is done on the base front only */

        /* Redistribution front sweep */
        for (i = 0; i < num_patches; i++)
        {
            /* 
            if(frs[i]->patch_level != frs[i]->NumberOfLevels-1
                && frs[i]->patch_level != 0)
                continue;
            */
            if(frs[i]->patch_level != 0) continue;
            if (debugging("hyp_amr"))
                printf("Doing redistribute advance patch front[%d],",i);
            status = redist_advance_front(dt,dt_frac,frs[i],&newfrs[i],
                                (POINTER)wvs[i]);
            if (debugging("hyp_amr"))
                printf("status = %d\n",status);
            switch(status)
            {
            case GOOD_STEP:
            break;

            case MODIFY_TIME_STEP:
            case REPEAT_TIME_STEP:
                printf("WARNING: Bad redist_advance_front(), status = %d\n",
                   status);
                goto redis_status;
            break;
            case ERROR_IN_STEP:
            default:
            screen("ERROR: overture_hyp_amr(),"
                   " redist_advance_front() failed\n");
            clean_up(ERROR);
            }
        }

redis_status:
        pp_gsync();
        status = syncronize_time_step_status(status,frs[0]->pp_grid);
        if (status != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
             " failed after redist_advance_front(), status = %d\n", status);
            clean_up(ERROR);
        }

        /* 
         * For the NOT-redistributed fronts, the
         * step is not updated. So update steps here. 
         */
        /* 
        Redistribution_count(frs[i]) = 
                  Redistribution_count(newfrs[0]); 
        */ 
        for (i = 1; i < num_patches; i++)
        { 
            newfrs[i]->step = newfrs[0]->step;
            newfrs[i]->time = newfrs[0]->time + dt;
        } 

        /* 
        rgr = frs[0]->rect_grid;
        cgr = computational_grid(frs[0]->interf);
        copy_rect_grid(cgr,rgr);
        tgr = &topological_grid(frs[0]->interf);
        tgr->Remap.remap = rgr->Remap.remap;
        set_patch_topo_grid(rgr,tgr);        

        rgr = newfrs[0]->rect_grid;
        cgr = computational_grid(newfrs[0]->interf);
        copy_rect_grid(cgr,rgr);
        tgr = &topological_grid(newfrs[0]->interf);
        tgr->Remap.remap = rgr->Remap.remap;
        set_patch_topo_grid(rgr,tgr);        

        if(n_redistr_fr->interf != newfrs[0]->interf ||
           redistr_fr->interf != frs[0]->interf) 
        {
            printf("ERROR: pseudo_advance_amr_fronts2d\n");
            printf("AFTER redistribution interface\n");
            printf("interface are not the same\n");
            clean_up(ERROR); 
        } 
        n_redistr_fr->interf = NULL;
        redistr_fr->interf = NULL;
        deep_free_front(n_redistr_fr); 
        deep_free_front(redistr_fr); 
        */  
        printf("After Redist: assembly_distribute_patch_fronts ");
        status = assembly_distribute_patch_fronts(newfrs,num_patches,
              redistr_table,max_n_patch, YES);
        printf("status = %d\n", status);          

        if (pp_min_status(status) != GOOD_STEP)
        {
            screen("ERROR overture_hyp_amr(),"
              " failed after redistr: assembly_distribute_patch_fronts()\n");
            clean_up(ERROR);
        }

        *newfronts = newfrs;
        DEBUG_LEAVE(pseudo_advance_amr_fronts2d)
        return status;
}

EXPORT void save_overture_show_files(
        Printplot    *prt,  
        Wave         *wave,
        int          step)
{

        char         *runname;
        doubleCompositeGridFunction *uu1;
        CompositeGrid *cg;
        Range        all;
        char         buf[100];
        int          myid; 
        size_t  sizest = wave->sizest; 
    
        debug_print("driver","Entering save_overture_show_files()\n");

        myid = pp_mynode();  
        runname = basename(prt->outfile);  
        sprintf(buf,"%s/%s_%s.%s.%d.ts%d",runname,"over",runname,"den",myid,step);

        if (create_directory(runname,YES) == FUNCTION_FAILED)
        {
            (void) printf("WARNING in save_overture_show_files(), directory "
                         "%s doesn't exist and can't be created\n",runname);
            return;
        }

        cg = (CompositeGrid*)wave->cg_over; 
        Ogshow show(buf);
        show.saveGeneralComment("The Density solution");

        if(wave->use_overture_state)
        {
            uu1 = (doubleCompositeGridFunction*)wave->cg_over_function;

            CompositeGridOperators op(*cg);

            /* define solution and operators: */
            doubleCompositeGridFunction Dens , Eng;
            doubleCompositeGridFunction Mom0, Mom1, Mom2;
            Dens.updateToMatchGrid(*cg);
            Eng.updateToMatchGrid(*cg);
            Mom0.updateToMatchGrid(*cg);
            Mom1.updateToMatchGrid(*cg);
            Mom2.updateToMatchGrid(*cg);

            copy_gridfunction_from_uu1(*uu1, *cg, 0, Dens);
           /* copy_gridfunction_from_uu1(*uu1, *cg, 1, Eng); */

            show.startFrame();
            /*show.saveSolution(Eng); */
            show.saveSolution(Dens);
        }
        else
        {
            doubleCompositeGridFunction Dens;
            Dens.updateToMatchGrid(*cg);  
        }  

        debug_print("driver","Leaving save_overture_show_files()\n");
}   /*  end save_overture_show_files()  */

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

#endif /* if defined(USE_OVERTURE)  */  
