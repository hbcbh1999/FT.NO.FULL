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
 *                               goverinterp.c:
 *
 *       Copyright 1999 by The University at Stony Brook, All rights reserved.
 *
 *       Contains routines for the interpolation of state data
 */

#define DEBUG_STRING      "goverinterp"

#include <gdecs/gdecs.h>

#if defined(USE_OVERTURE)

enum {
        STRICT_INTERIOR = 1,
        INTERIOR,
        BUFFER_ZONE, 
        EXTERIOR
};

LOCAL  void     copy_base_grid_soln_to_new(Wave*, Wave*);
LOCAL  void     copy_gridfunction_from_uu1(const doubleCompositeGridFunction&,
                  const CompositeGrid&, int, doubleCompositeGridFunction&); 
LOCAL  int      same_point2d(double*,double*,double);  
LOCAL  int      amr_buffer_zone_block(int*,int*,int,RECT_GRID*,INTERFACE*); 

LOCAL  void     fill_extr_overture_upper_cell_st(const CompositeGrid&,
                   doubleCompositeGridFunction&); 
LOCAL  Locstate init_coarse_to_fine(int*,double*,COMPONENT,Wave*,RECT_GRID*,Front*);
LOCAL  double    pt_distance2d(double*, double*);

LOCAL  Locstate fine_to_fine(int*,double*,COMPONENT,Wave*,RECT_GRID*,Front*);  
LOCAL  Locstate fine_to_coarse(int*,double*,COMPONENT,Wave*,RECT_GRID*,Front*,Locstate); 
LOCAL  void     fill_root_extr_overture_cell_st(Wave**,Front**); 
LOCAL  int      pt_in_patch_same_or_coarse_level(int*, double*, int,
                    int, int*, Wave**, Front**);
LOCAL  int      pt_in_patch_same_level_repatch(int*,double*,int,int,
                   int*,Wave**,Front**);
LOCAL  int      pt_in_patch_coarse_level_repatch(int*,double*,int,int,
                   int*,Wave**,Front**);
LOCAL  int      pt_in_others_same_or_coarse_level(int*,double*,int,
                   int,int*,Wave**,Front**); 
LOCAL  int      best_match_pt_in_other_sets(int*,double*,int,int,int*,
                     int*,Wave**,Front**,int);  
LOCAL  int      best_match_pt_in_coarse(int*,double*,int,int*,int*,
                     Wave**,Front**,int); 
LOCAL  int      best_match_pt_in_same_set(int*,double*,
                      int,int,int*,int*,Wave**,Front**);   
LOCAL   void    set_patch_domain_in_sweep_dir(int*,int*,int*,int*,int*,
                    int,int,Overparam*,Wv_on_pc**,Front*);
LOCAL   void    set_receive_domain(int*,int*,int*,int,int,RECT_GRID*);
LOCAL   boolean    g_overture_injection_FT_buf_after_repatch(Wave**,Front**,
                   Wave**,Front**); 

LOCAL double pt_distance2d(
        double   *crds0,
        double   *crds1)
{
        return sqrt(sqr(crds0[0]-crds1[0]) + sqr(crds0[1]+crds1[1]));
}

LOCAL int same_point2d(
        double   *p1,
        double   *p2,
        double   ls)
{
        const double delta = 1.0e-6; /*TOLERANCE*/
        int i;

        for (i = 0; i < 2; i++)
        {
           if (fabs(p1[i] - p2[i]) > delta*ls)
             return NO;
        }
        return YES;
}

EXPORT void g_overture_to_ft_st(
	Locstate    st,
        POINTER     cg_f,
        int         id, 
        int         *ic)
{
	doubleCompositeGridFunction *cg_function
          = (doubleCompositeGridFunction*)cg_f;  
 
        Dens(st) = (*cg_function)[id](0,ic[0],ic[1],0);
        Energy(st)  = (*cg_function)[id](1,ic[0],ic[1],0);
        Mom(st)[0] = (*cg_function)[id](2,ic[0],ic[1],0);
        Mom(st)[1] = (*cg_function)[id](3,ic[0],ic[1],0);
        reset_gamma(st);
}

EXPORT void g_ft_to_overture_st(
	Locstate    st,
        POINTER     cg_f,
        int         id, 
        int         *ic)
{
	doubleCompositeGridFunction *cg_function
          = (doubleCompositeGridFunction*)cg_f;  
 
        if(Params(st) != NULL)
        {
            (*cg_function)[id](0,ic[0],ic[1],0) = Dens(st);
            (*cg_function)[id](1,ic[0],ic[1],0) = Energy(st);
            (*cg_function)[id](2,ic[0],ic[1],0) = Mom(st)[0];
            (*cg_function)[id](3,ic[0],ic[1],0) = Mom(st)[1];
        }
        else
        {
            (*cg_function)[id](0,ic[0],ic[1],0) = 0.0;
            (*cg_function)[id](1,ic[0],ic[1],0) = 0.0;
            (*cg_function)[id](2,ic[0],ic[1],0) = 0.0;
            (*cg_function)[id](3,ic[0],ic[1],0) = 0.0;
        }  
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


/* for a coarse grid pt, obtain 4 surranding pts on the next fine grid
 * based on the assumption that the refinement ratio is 2.
 */
LOCAL  Locstate fine_to_coarse(
        int             icoords0[MAXD], /* coarse pt */
        double           coords0[MAXD],  /* coarse pt */
        COMPONENT       comp0,          /* coarse pt */
        Wave            *wave,
        RECT_GRID       *grid,
        Front           *front,
        Locstate        dfltst)
{
        static Locstate st = NULL;
        size_t          sizest = front->sizest;
        Locstate        st1;
        COMPONENT       comp1;
        double           coords1[MAXD];
        int             icoords1[MAXD], ic[MAXD];
        int             dim = wave->rect_grid->dim;
        int             i, point_num;
        RECT_GRID       *comp_g = &wave_tri_soln(wave)->tri_grid->comp_grid;
        RECT_GRID       *tg = &wave_tri_soln(wave)->tri_grid->tg_grid;
        RECT_GRID       *rect_g = &wave_tri_soln(wave)->tri_grid->rect_grid;


        if(st == NULL)
        {
            alloc_state(front->interf,&st,max((int)sizeof(VGas),(int)sizest));
        }

        switch(dim)
        {
#if defined(ONED)
            case 1:
            {
                printf("ERROR, no implimentation of 1D fine_to_coarse()\n");
                clean_up(ERROR);
                break;
            }
#endif /* defined(ONED) */
#if defined(TWOD)
            case 2:
            {
                /* note : use tg_grid */  
                if(rect_in_which(coords0,icoords1,tg) == FUNCTION_FAILED)
                {
                    printf("ERROR in fine_to_coarse(), "
                                  "rect_in_which() failed\n");
                    printf("Looking for coords<%g, %g> in tri->tg_grid\n",
                                    coords0[0],coords0[1]);
                    print_rectangular_grid(tg);
                    clean_up(ERROR);
                }

                ic[0] = icoords1[0], ic[1] = icoords1[1];
                if(Rect_coords(ic,wave)[0] < coords0[0] or
                   Rect_coords(ic,wave)[1] < coords0[1])
                {
                    printf("ERROR in fine_to_coarse(), not right_upper\n");
                    printf("crds<%g, %g> ic found is <%d, %d><%g,%g>\n",
                            coords0[0],coords0[1], ic[0], ic[1],
                            Rect_coords(ic,wave)[0],Rect_coords(ic,wave)[1]);
                    clean_up(ERROR);
                }
                ic[0] = icoords1[0]-1, ic[1] = icoords1[1];
                if(Rect_coords(ic,wave)[0] > coords0[0] or
                   Rect_coords(ic,wave)[1] < coords0[1])
                {
                    printf("ERROR in fine_to_coarse(), not left_upper\n");
                    printf("crds<%g, %g> ic found is <%d, %d><%g,%g>\n",
                            coords0[0],coords0[1], ic[0], ic[1],
                            Rect_coords(ic,wave)[0],Rect_coords(ic,wave)[1]);
                    clean_up(ERROR);
                }
                ic[0] = icoords1[0]-1, ic[1] = icoords1[1]-1;
                if(Rect_coords(ic,wave)[0] > coords0[0] or
                   Rect_coords(ic,wave)[1] > coords0[1])
                {
                    printf("ERROR in fine_to_coarse(), not left_lower\n");
                    printf("crds<%g, %g> ic found is <%d, %d><%g,%g>\n",
                            coords0[0],coords0[1], ic[0], ic[1],
                            Rect_coords(ic,wave)[0],Rect_coords(ic,wave)[1]);
                    clean_up(ERROR);
                }
                ic[0] = icoords1[0], ic[1] = icoords1[1]-1;
                if(Rect_coords(ic,wave)[0] < coords0[0] or
                   Rect_coords(ic,wave)[1] > coords0[1])
                {
                    printf("ERROR in fine_to_coarse(), not right_lower\n");
                    printf("crds<%g, %g> ic found is <%d, %d><%g,%g>\n",
                            coords0[0],coords0[1], ic[0], ic[1],
                            Rect_coords(ic,wave)[0],Rect_coords(ic,wave)[1]);
                    clean_up(ERROR);
                }


                if(debugging("xxx_fine_to_coarse"))
                {
                     printf("coarse crds<%g, %g> ic<%d,%d> comp[%d]"
                             " found on fine level[%d]\n",
                             coords0[0],coords0[1], icoords0[0],
                             icoords0[1], comp0, wave->patch_level);

                     printf("fine ic0<%d, %d><%g, %g> comp[%d]\n",
                             icoords1[0], icoords1[1],
                       Rect_coords(icoords1,wave)[0],Rect_coords(icoords1,wave)[1],
                       Rect_comp(ic,wave));
                     ic[0] = icoords1[0]-1, ic[1] = icoords1[1];
                     printf("fine ic1<%d, %d><%g, %g> comp[%d]\n",
                        ic[0],ic[1],Rect_coords(ic,wave)[0],
                        Rect_coords(ic,wave)[1],Rect_comp(ic,wave));

                     ic[0] = icoords1[0]-1, ic[1] = icoords1[1]-1;
                     printf("fine ic1<%d, %d><%g, %g> comp[%d]\n",
                        ic[0],ic[1],Rect_coords(ic,wave)[0],
                        Rect_coords(ic,wave)[1],Rect_comp(ic,wave));

                     ic[0] = icoords1[0], ic[1] = icoords1[1]-1;
                     printf("fine ic1<%d, %d><%g, %g> comp[%d]\n",
                        ic[0],ic[1],Rect_coords(ic,wave)[0],
                        Rect_coords(ic,wave)[1],Rect_comp(ic,wave));
                }
                break;
            }
#endif /* defined(TWOD) */
#if defined(THREED)
            case 3:
            {
                printf("ERROR, no implimentation of 3D fine_to_coarse()\n");
                clean_up(ERROR);
                break;
            }
#endif /* defined(THREED) */
        }  /* end of switch */

        point_num = 0;
        Dens(st) = 0.0;
        Energy(st) = 0.0;
        for (i = 0; i < dim; i++)
             Mom(st)[i] = 0.0;

        ic[0] = icoords1[0], ic[1] = icoords1[1];
        comp1 = Rect_comp(ic,wave);
        if (comp1 == comp0)
        {
            point_num++;
            st1 = Rect_state(ic,wave);

            Dens(st) += Dens(st1);
            Energy(st) += Energy(st1);
            for (i = 0; i < dim; i++)
                Mom(st)[i] += Mom(st1)[i];
            reset_gamma(st);
            Set_params(st,st1);
            set_type_of_state(st,state_type(st1));
        }

        ic[0] = icoords1[0]-1, ic[1] = icoords1[1];
        comp1 = Rect_comp(ic,wave);
        if (comp1 == comp0)
        {
            point_num++;
            st1 = Rect_state(ic,wave);

            Dens(st) += Dens(st1);
            Energy(st) += Energy(st1);
            for (i = 0; i < dim; i++)
                Mom(st)[i] += Mom(st1)[i];
            reset_gamma(st);
            Set_params(st,st1);
            set_type_of_state(st,state_type(st1));
        }

        ic[0] = icoords1[0]-1, ic[1] = icoords1[1]-1;
        comp1 = Rect_comp(ic,wave);
        if (comp1 == comp0)
        {
            point_num++;
            st1 = Rect_state(ic,wave);

            Dens(st) += Dens(st1);
            Energy(st) += Energy(st1);
            for (i = 0; i < dim; i++)
                Mom(st)[i] += Mom(st1)[i];
            reset_gamma(st);
            Set_params(st,st1);
            set_type_of_state(st,state_type(st1));
        }

        ic[0] = icoords1[0], ic[1] = icoords1[1]-1;
        comp1 = Rect_comp(ic,wave);
        if (comp1 == comp0)
        {
            point_num++;
            st1 = Rect_state(ic,wave);

            Dens(st) += Dens(st1);
            Energy(st) += Energy(st1);
            for (i = 0; i < dim; i++)
                Mom(st)[i] += Mom(st1)[i];
            reset_gamma(st);
            Set_params(st,st1);
            set_type_of_state(st,state_type(st1));
        }

        if (point_num > 0)
        {
            Dens(st) = Dens(st)/point_num;
            Energy(st) = Energy(st)/point_num;
            for (i = 0; i < dim; i++)
                Mom(st)[i] = Mom(st)[i]/point_num;
            reset_gamma(st);
        }
        else
        {
            if(not is_obstacle_state(dfltst))
                ft_assign(st, dfltst, front->sizest);
            else
            {
                screen("FT_ERROR in fine_to_coarse(), bad interpolation\n");
                printf("dfltst is a obstacle_state\n");
                printf("For crds0<%g, %g> ic0<%d, %d> comp[%d]\n",
                        coords0[0],coords0[1], icoords0[0], icoords0[1], comp0);
                printf("surrounded by fine grid pts:\n");
                coords1[0] = Rect_coords(icoords1,wave)[0];
                coords1[1] = Rect_coords(icoords1,wave)[1];
                comp1 = Rect_comp(icoords1,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d>\n",
                     wave->patch_number, coords1[0],coords1[1],icoords1[0],icoords1[1]);
                screen("comp1=%d\n",comp1);

                ic[0] = icoords1[0]-1, ic[1] = icoords1[1];
                coords1[0] = Rect_coords(ic,wave)[0];
                coords1[1] = Rect_coords(ic,wave)[1];
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d>\n",
                     wave->patch_number, coords1[0],coords1[1],ic[0],ic[1]);
                screen("comp1=%d\n",comp1);

                ic[0] = icoords1[0]-1, ic[1] = icoords1[1]-1;
                coords1[0] = Rect_coords(ic,wave)[0];
                coords1[1] = Rect_coords(ic,wave)[1];
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d>\n",
                     wave->patch_number, coords1[0],coords1[1],ic[0],ic[1]);
                screen("comp1=%d\n",comp1);

                ic[0] = icoords1[0]-1, ic[1] = icoords1[1]-1;
                coords1[0] = Rect_coords(ic,wave)[0];
                coords1[1] = Rect_coords(ic,wave)[1];
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d>\n",
                     wave->patch_number, coords1[0],coords1[1],ic[0],ic[1]);
                screen("comp1=%d\n",comp1);

                ic[0] = icoords1[0], ic[1] = icoords1[1]-1;
                coords1[0] = Rect_coords(ic,wave)[0];
                coords1[1] = Rect_coords(ic,wave)[1];
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d>\n",
                     wave->patch_number, coords1[0],coords1[1],ic[0],ic[1]);
                screen("comp1=%d\n",comp1);
                /*  
                print_components(wave_tri_soln(wave)->tri_grid);
                */  
                clean_up(ERROR);
            }
        }
        return st;
}

/****  g_scatter_patch_states()  ****/
/* 
* This scatter amr patch buffer state function call is used 
* exclusively in the overture_parab_amr_driver(), after
* h_scatter_patch_states() is performed in every direction. 
* The difference between this function call and
* g_scatter_patch_states_in_sweep_dir() is that 
* g_scatter_patch_states() takes care of all the buffer
* zone points of a single patch at one time regardless the 
* sweep direction. Because the overture_parab_amr_driver()
* is not a directional splitting sweep driver.  
* Keeping input param "iperm" is for the purpose of
* convenience, it could be removed.    
*/

EXPORT boolean g_scatter_patch_states(
        Overparam    *overparam,
        Wv_on_pc     **redistr_table,    
	Wave 	     **wvs,
	Front        **frs,
        int          *iperm)
{
        int          num_patches, levels;
        int          i, dim, side, swp;
        int          imin[MAXD], imax[MAXD];
        int          L[MAXD], U[MAXD];
        RECT_GRID    *gr, *gridfrom;
        int          ix, iy, ic[MAXD];
        int          patch_num0;
        Wave         *wave0, *wavefrom;
        Front        *frontfrom;
        Locstate     st;
        int          on_patch, where, where_level; 

        DEBUG_ENTER(g_scatter_patch_states)

        num_patches = wvs[0]->totalNumberOfPatches;
        levels = wvs[0]->NumberOfLevels;
        gr = wvs[0]->rect_grid;
        dim = gr->dim;

        for (i = 0; i < dim; ++i)
        {
            imin[i] = -gr->lbuf[i];
            imax[i] = gr->gmax[i] + gr->ubuf[i];
        }

        if(debugging("parab_amr_driver"))
            printf("Patch[0] imin[%d, %d], imax[%d, %d]\n",
                 imin[0], imin[1], imax[0], imax[1]);
        for(i = 1; i < num_patches; i++)
        {

            if(debugging("parab_amr_driver"))
            {
                printf("patch[%d], level[%d], rect_grid\n", i, wvs[i]->patch_level);
                print_rectangular_grid(wvs[i]->rect_grid);
            }
            for(swp = 0; swp < dim; swp++)
            {
                for(side = 0; side < 2; side++)
                {
                    set_patch_domain_in_sweep_dir(L,U,imin,imax,
                      iperm,swp,side,overparam,redistr_table,frs[i]);  
                    if(rect_boundary_type(frs[i]->interf,iperm[swp],side) !=
                           AMR_SUBDOMAIN_BOUNDARY &&
                       rect_boundary_type(frs[i]->interf,iperm[swp],side) !=
                           SUBDOMAIN_BOUNDARY)
                        continue; 
                    if(debugging("parab_amr_driver"))
                    {
                        printf("Dir[%d] Side[%d]:Bdry_type[%d] amr_buf: L[%d,%d] U[%d,%d]\n",
                         iperm[swp], side, rect_boundary_type(frs[i]->interf,iperm[swp],side),
                         L[0], L[1], U[0], U[1]);
                    }
                    if(rect_boundary_type(frs[i]->interf,iperm[swp],side) ==
                         SUBDOMAIN_BOUNDARY)
                    {
                        for(iy = L[1]; iy < U[1]; iy++)
                        {
                            for(ix = L[0]; ix < U[0]; ix++)
                            {
                                ic[1] = iy; ic[0] = ix;
                                st = Rect_state(ic,wvs[i]);
                                if(NULL == Params(st))
                                {
                                    best_match_pt_in_coarse(ic,Rect_coords(ic,wvs[i]),
                                      wvs[i]->patch_level,&on_patch,&where,wvs,frs,YES);
                                    wavefrom = wvs[on_patch];
                                    frontfrom = frs[on_patch];
                                    hyp_solution(Rect_coords(ic,wvs[i]),
					Rect_comp(ic,wvs[i]),NULL,UNKNOWN_SIDE,
					frontfrom,wavefrom,st,NULL);
                                }
                            } 
                        }
                    }
                    else
                    {
                        for(iy = L[1]; iy < U[1]; iy++)
                        {
                            for(ix = L[0]; ix < U[0]; ix++)
                            {
                                ic[1] = iy; ic[0] = ix;
                                st = Rect_state(ic,wvs[i]);
                                where_level = best_match_pt_in_same_set(
                                   ic, Rect_coords(ic,wvs[i]), Rect_comp(ic,wvs[i]),
                                  wvs[i]->patch_number,&on_patch,&where,wvs,frs);
                                wavefrom = wvs[on_patch];
                                frontfrom = frs[on_patch];
                                gridfrom = wavefrom->rect_grid;
                                if(wvs[i]->patch_level == where_level)
                                {
                                    ft_assign(st, fine_to_fine(ic,Rect_coords(ic,wvs[i]),
                                        Rect_comp(ic,wvs[i]),wavefrom,gridfrom,frontfrom),
                                       frontfrom->sizest);
                                }
                                else
                                {
                                    if(!is_excluded_comp(Rect_comp(ic,wvs[i]),frs[i]->interf))
                                    {
                                        hyp_solution(Rect_coords(ic,wvs[i]),
					    Rect_comp(ic,wvs[i]),NULL,
					    UNKNOWN_SIDE,frontfrom,
					    wavefrom,st,NULL);
                                    }
                                    else
                                    {
                                        g_obstacle_state(st,frs[i]->sizest);
                                    }
                                }
                            }
                        }
                    }  /* end of else */
                }
            } 
        }  

        DEBUG_LEAVE(g_scatter_patch_states)
	return FUNCTION_SUCCEEDED; 
} 

/*
*               g_overture_undistribute_interpolation_fine_to_coarse()
*
*     interpolate back the interior states in coarser patches from finer patches
*       after each time step, after on-the-fly patches are gathered back. 
*/

EXPORT boolean g_overture_undistribute_interpolation_fine_to_coarse(
        Wave            **wvs,
        Front           **frs)
{
        Wave            *wavefrom, *waveto;
        Front           *frontfrom, *frontto;
        RECT_GRID       *gridfrom, *gridto;
        Locstate        st;
        double           coords0[MAXD], coords1[MAXD];
        int             icoords0[MAXD], icoords1[MAXD];
        int             dim = wvs[0]->rect_grid->dim;
        int             num_patches, patch_num0, patch_num1;
        int             levels, patch_level0, patch_level1;
        int             i, int1;
        int             int2;
        COMPONENT       comp0,comp1;
        int             smin[MAXD], smax[MAXD];
        int             lbuf[MAXD], ubuf[MAXD];
        size_t          sizest = wvs[0]->sizest;

        DEBUG_ENTER(g_overture_undistribute_interpolation_fine_to_coarse)

        num_patches = wvs[0]->totalNumberOfPatches;
        levels = wvs[0]->NumberOfLevels;

        /* interpolation of interior states in coarse 
         * patches from finest patches */

        for (int1 = num_patches-1; int1 >= 0; int1--)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;
            if ( waveto->patch_level == levels-1) continue;

            for (i = 0; i < dim; i++)
            {
                smin[i] = 0;
                smax[i] = gridto->gmax[i];
            }
            patch_num0 = waveto->patch_number;
            patch_level0 = waveto->patch_level;

            if(debugging("amr_fine_to_coarse"))
            {
                 printf("In g_overture_interpolation_fine_to_coarse() "
                    " inner-patch states interpolation of the patch[%d],level[%d] \n",
                     patch_num0, patch_level0);
                 printf("patch[%d] state before \n",waveto->patch_number);
            }

            switch(dim)
            {
#if defined(ONED)
            case 1:
                {
                    printf("ERROR: no implimentation of 1D"
                        " g_overture_interpolation_fine_to_coarse()\n");
                    clean_up(ERROR);
                    break;
                }
#endif /* defined(ONED) */
#if defined(TWOD)
            case 2:
                {
                    int      ix, iy;
                    for (iy = smin[1]; iy < smax[1]; iy++)
                    {
                        icoords0[1] = iy;
                        for (ix = smin[0]; ix < smax[0]; ix++)
                        {
                            icoords0[0] = ix;
                            st = Rect_state(icoords0,waveto);
                            coords0[0] = Rect_coords(icoords0,waveto)[0];
                            coords0[1] = Rect_coords(icoords0,waveto)[1];
                            comp0 = Rect_comp(icoords0,waveto);

                            if(!is_excluded_comp(comp0,frontto->interf))
                            {
                                for ( int2 = num_patches - 1; int2 > int1; int2-- )
                                {
                                    patch_level1 = wvs[int2]->patch_level;
                                    if (patch_level1 == patch_level0 + 1)
                                    {
                                        gridfrom = wvs[int2]->rect_grid;
                                        wavefrom = wvs[int2];
                                        frontfrom = frs[int2];
                                        if ( coords0[0] >= gridfrom->L[0] and
                                             coords0[0] <= gridfrom->U[0] and
                                             coords0[1] >= gridfrom->L[1] and
                                             coords0[1] <= gridfrom->U[1])
                                        {
                                            ft_assign(st, fine_to_coarse(icoords0,coords0,comp0,
                                               wavefrom,gridfrom,frontfrom,st),frontfrom->sizest);
                                            break;
                                        }
                                    }
                                } /* End: for (int2 = num_patches - 1;...) */ 
                            }
                        }
                    }
                    break;
                }
#endif /* defined(TWOD) */
#if defined(THREED)
            case 3:
                {
                    printf("ERROR: no implimentation of 3D"
                            " g_overture_interpolation_fine_to_coarse()\n");
                    clean_up(ERROR);
                    break;
                }
#endif /* defined(THREED) */
            }   /* end of switch */  
            if(debugging("amr_fine_to_coarse"))
            {
                Index I1,I2,I3;
                doubleCompositeGridFunction Dens;
                doubleCompositeGridFunction *cg_f;
                CompositeGrid               *cg;  
                if(waveto->use_overture_state)
                {
                    cg = (CompositeGrid*)waveto->cg_over; 
                    cg_f = (doubleCompositeGridFunction*)waveto->cg_over_function; 

                    printf("In g_overture_interpolation_fine_to_coarse(),"
                      " patch[%d] state after \n", waveto->patch_number);
                    getIndex((*cg)[waveto->patch_number].indexRange(),I1,I2,I3);
                    printf("x dir range<%d,%d>, y dir range<%d,%d>\n",
                        I1.getBase(), I1.getBound(), I2.getBase(), I2.getBound());
                    Dens.updateToMatchGrid(*cg);
                    copy_gridfunction_from_uu1(*cg_f, *cg,0,Dens);
                }
            }
        }

        if(debugging("gas_param_interp"))
        {
            for (int1 = num_patches-1; int1 >= 0; int1--)
            {
                int ic[2] = {0,0};
                int      ix, iy;
                waveto = wvs[int1];
                frontto = frs[int1];
                gridto = waveto->rect_grid;
                if ( waveto->patch_level == levels-1) continue;

                for (i = 0; i < dim; i++)
                {
                    smin[i] = 0;
                    smax[i] = gridto->gmax[i];
                    lbuf[i] = gridto->lbuf[i];
                    ubuf[i] = gridto->ubuf[i];
                }
                for (iy = smin[1]; iy < smax[1]; iy++)
                {
                    icoords0[1] = iy;
                    for (ix = smin[0]; ix < smax[0]; ix++)
                    {
                        icoords0[0] = ix;
                        st = Rect_state(icoords0,waveto);
                        coords0[0] = Rect_coords(icoords0,waveto)[0];
                        coords0[1] = Rect_coords(icoords0,waveto)[1];
                        comp0 = Rect_comp(icoords0,waveto);
                        if(Different_params(st,Rect_state(ic,waveto)))
                        {
                            screen("exit In g_overture_interpolation_fine_to_coarse\n");
                            screen("crds<%g,%g> ic<%d,%d> patch[%d] is diff params\n",
                                coords0[0], coords0[1], icoords0[0], icoords0[1], int1);
                            clean_up(ERROR);
                        }
                        if(is_obstacle_state(st))
                        {
                            screen("exit In g_overture_interpolation_fine_to_coarse\n");
                            screen("crds<%g,%g> ic<%d,%d> patch[%d] is obstacle\n",
                            coords0[0], coords0[1], icoords0[0], icoords0[1], int1);
                            clean_up(ERROR);
                        }
                    }
                }
            }
        }

        DEBUG_LEAVE(g_overture_undistribute_interpolation_fine_to_coarse)
        return FUNCTION_SUCCEEDED;  
}


EXPORT boolean g_overture_interpolation_fine_to_coarse(
        Wv_on_pc        **redistr_table,
        Wave            **wvs,
        Front           **frs,
        Wave            ***undistr_wvs,
        Front           ***undistr_frs, 
        int             num_cal_patch,
        int             max_n_patch)
{
        Wave            **tmpwvs;
        Front           **tmpfrs;  
        int             i, num_patches, resid_n; 
        boolean            status;  

        DEBUG_ENTER(g_overture_interpolation_fine_to_coarse)

        /* 
        reinstore_mini_undistribute_patch(&tmpwvs,&tmpfrs,&resid_n,
             wvs,frs,redistr_table,num_cal_patch,max_n_patch, YES);
        */
        reinstore_undistribute_patch(&tmpwvs,&tmpfrs,wvs,frs,
                 redistr_table,num_cal_patch,max_n_patch); 

        status = g_overture_undistribute_interpolation_fine_to_coarse(
            tmpwvs, tmpfrs);

        if (pp_min_status(status) != YES)
        {
            screen("ERROR g_overture_interpolation_fine_to_coarse(),"
                   " failed after local undistribute_interpolate\n");
            clean_up(ERROR);
        }  
        /* 
        num_patches = tmpwvs[0]->totalNumberOfPatches;
        */  
        *undistr_wvs = tmpwvs;
        *undistr_frs = tmpfrs;  

        DEBUG_LEAVE(g_overture_interpolation_fine_to_coarse)
        return FUNCTION_SUCCEEDED;  
}

/* this function call should use boundary info. This is a
   tmp solution.
*/
LOCAL void fill_root_extr_overture_cell_st(
        Wave       **wvs,
        Front      **frs)
{

}  

/*    best_match_pt_in_same_set() */   
/* The point is looked in the same level or the coarer grid level.   
* The current implimentation of this function is for 
*  the buffer zone point interpolation states. */ 
LOCAL int best_match_pt_in_same_set(
        int             *icoords0,
        double           *coords0,
        int             comp0,
        int             patch_num0,
        int             *on_patch,
        int             *where, 
        Wave            **wvs,
        Front           **frs)
{
        Wave            *waveto;
        Front           *frontto;
        RECT_GRID       *gridto;
        double           coords1[MAXD];
        int             icoords1[MAXD];
        int             num_patches, patch_num1;
        int             patch_level0, patch_level1;
        double           VL[MAXD], VU[MAXD];  
        double           dh[MAXD];  
        int             i;
         

        num_patches = wvs[0]->totalNumberOfPatches;
        patch_level0 = wvs[patch_num0]->patch_level;

        /* First check if pt is in interior of other patches. */ 
        for (i =num_patches-1; i >= 0 ; i--)
        {
            waveto = wvs[i];
            frontto = frs[i];
            gridto = waveto->rect_grid;

            if(waveto->patch_level > patch_level0) continue;
            if(patch_num0 == waveto->patch_number) continue;

            if ( coords0[0] >= gridto->L[0] and coords0[0] <= gridto->U[0] and
                 coords0[1] >= gridto->L[1] and coords0[1] <= gridto->U[1])
            {
                *on_patch = waveto->patch_number;
                *where = STRICT_INTERIOR;   
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in best_match_pt_in_same_set(), "
                           "rect_in_which() failed\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(icoords1[0] < 0 or icoords1[1] < 0)
                {
                    printf("ERROR in best_match_pt_in_same_set(), "
                           "pt found in the buffer\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1], icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                return waveto->patch_level;
            }
        }

        /* check coarse grid with buffer zone included */ 
        for (i =num_patches-1; i >= 0 ; i--)
        {
            waveto = wvs[i];
            frontto = frs[i];
            gridto = waveto->rect_grid;

            if(waveto->patch_level >= patch_level0) continue;
            if(patch_num0 == waveto->patch_number) continue;

            for(int dir = 0; dir < gridto->dim; dir++)
                dh[dir] = gridto->h[dir]; 
            for(int dir = 0; dir < gridto->dim; dir++)
            {
                if(gridto->lbuf[dir] > 0)
                    VL[dir] = gridto->VL[dir]+dh[dir]/2.0;
                else
                    VL[dir] = gridto->VL[dir];
                if(gridto->ubuf[dir] > 0)
                    VU[dir] = gridto->VU[dir]-dh[dir]/2.0;
                else
                    VU[dir] = gridto->VU[dir];
            }
            if ( coords0[0] >= VL[0] and
                 coords0[0] <= VU[0] and
                 coords0[1] >= VL[1] and
                 coords0[1] <= VU[1])  
            {
                *on_patch = waveto->patch_number;
                *where = INTERIOR;   
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in best_match_pt_in_same_set(), "
                           "rect_in_which() failed\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                return waveto->patch_level;
            }
        }

        printf("ERROR: best_match_pt_in_same_set\n");  
        printf("print patch[%d] wave rect_grid\n", patch_num0);
        print_rectangular_grid(wvs[patch_num0]->rect_grid);
        clean_up(ERROR);
        return NO;
}

/* assume the patch is properly nested.
 * So, if a point of a patch can not be found on the other patches at
 * the same level, it must be contained in the next coarser level.
 * If no patch is found, an error occur.
 * Note: if YES, this pt should be in the interior of the found patch.
 * If the return NO, pt is in the corser buffer
 */

LOCAL int pt_in_patch_same_or_coarse_level(
        int             *icoords0,
        double           *coords0,
        int             comp0,
        int             patch_num0,
        int             *on_patch,
        Wave            **wvs,
        Front           **frs)
{
        Wave            *waveto;
        Front           *frontto;
        RECT_GRID       *gridto;
        double           coords1[MAXD];
        int             icoords1[MAXD];
        int             num_patches, patch_num1;
        int             patch_level0, patch_level1;
        int             i, int1;
        COMPONENT       comp1;
        CompositeGrid   *cg = (CompositeGrid*)wvs[0]->cg_over;

        num_patches = wvs[0]->totalNumberOfPatches;
        patch_level0 = wvs[patch_num0]->patch_level;

        /* First check if pt is in interior of other patches. */ 
        for (int1 =num_patches-1; int1 >= 0 ; int1--)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;

            if(waveto->patch_level > patch_level0) continue;
            if(patch_level0 - waveto->patch_level > 1) continue;
            if(patch_num0 == waveto->patch_number) continue;

            if ( coords0[0] >= gridto->L[0] and coords0[0] <= gridto->U[0] and
                 coords0[1] >= gridto->L[1] and coords0[1] <= gridto->U[1])
            {
                *on_patch = waveto->patch_number;
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level(), "
                           "rect_in_which() failed\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(icoords1[0] < 0 or icoords1[1] < 0)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level(), "
                           "pt found in the buffer\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1], icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(waveto->patch_level != patch_level0 and
                   patch_level0-waveto->patch_level != 1)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level\n");
                    printf("the levels are not right\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d>\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    printf("pt level[%d], patch level[%d]\n", patch_level0,
                           waveto->patch_level);
                    clean_up(ERROR);
                }
                return 0;
            }
        }
        /* Then check if pt is in FT buffer zone of coarse patches. */  
        for (int1 =num_patches-1; int1 >= 0 ; int1--)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;

            if(waveto->patch_level > patch_level0 or
               waveto->patch_level == patch_level0) continue;
            if(patch_level0 - waveto->patch_level > 1) continue;
            if(patch_num0 == waveto->patch_number) continue;

            if ( coords0[0] >= gridto->VL[0] and coords0[0] <= gridto->VU[0] and
                 coords0[1] >= gridto->VL[1] and coords0[1] <= gridto->VU[1])
            {
                *on_patch = waveto->patch_number;
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level(), "
                           "rect_in_which() failed for next coarse grid\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(patch_level0-waveto->patch_level != 1)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level\n");
                    printf("pt should be in the next coarse(maybe in buffer),"
                          " the levels are not right\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d>\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    printf("pt level[%d], patch level[%d]\n", patch_level0,
                           waveto->patch_level);
                    clean_up(ERROR);
                }
                return 1;
            }
        }
/* Overture buffer zone always contains FT buffer zone.
   check if pt is in oveture buffer zone of coarse patches. 
*/ 
        Index I1, I2, I3;
        int   imin[3], imax[3];
        double LL[3], UU[3];
        for (int1 =num_patches-1; int1 >= 0 ; int1--)
        {
            waveto = wvs[int1];
            if(patch_num0 == waveto->patch_number) continue;
            if(waveto->patch_level > patch_level0 or
               waveto->patch_level == patch_level0) continue;
            if(patch_level0 - waveto->patch_level > 1) continue;

            getIndex((*cg)[int1].dimension(),I1,I2,I3);
            imin[0] = I1.getBase();
            imin[1] = I2.getBase();
            imax[0] = I1.getBound();
            imax[1] = I2.getBound();
            LL[0] = (*cg)[int1].vertex()(imin[0],imin[1],0,axis1);
            LL[1] = (*cg)[int1].vertex()(imin[0],imin[1],0,axis2);
            UU[0] = (*cg)[int1].vertex()(imax[0],imax[1],0,axis1);
            UU[1] = (*cg)[int1].vertex()(imax[0],imax[1],0,axis2);
            if ( coords0[0] >= LL[0] and coords0[0] <= UU[0] and
                 coords0[1] >= LL[1] and coords0[1] <= UU[1])
            {
                *on_patch = waveto->patch_number;
                if(patch_level0-waveto->patch_level != 1)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level\n");
                    printf("pt should be in the next coarse( overture buffer),"
                          " the levels are not right\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d>\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    printf("pt level[%d], patch level[%d]\n", patch_level0,
                           waveto->patch_level);
                    clean_up(ERROR);
                }
                return 2;
            }
        }

        printf("ERROR in pt_in_patch_same_or_coarse_level\n");
        printf("Looking for coords<%g, %g>ic<%d,%d>\n",
                coords0[0], coords0[1],icoords0[0],icoords0[1]);
        printf("looped all the patches, but not found\n");
        for(int axis = 0; axis<cg->numberOfDimensions(); axis++)
        {
            for(int side = Start; side<=End; side++)
            {
                printf("patch[%d], axis<%d>, side<%d> boundary condition %d\n",
                  patch_num0, axis, side,
                (*cg)[patch_num0].boundaryCondition()(side,axis));
            }
        }
        printf("print patch[%d] wave rect_grid\n", patch_num0);
        print_rectangular_grid(wvs[patch_num0]->rect_grid);
        clean_up(ERROR);
        return NO;
}

LOCAL int pt_in_others_same_or_coarse_level(
        int             *icoords0,
        double           *coords0,
        int             comp0,
        int             patch_level0,
        int             *on_patch,
        Wave            **wvs,
        Front           **frs)
{
        Wave            *waveto;
        Front           *frontto;
        RECT_GRID       *gridto;
        double           coords1[MAXD];
        int             icoords1[MAXD];
        int             num_patches, patch_num1;
        int             patch_level1;
        int             i, int1;
        COMPONENT       comp1;
        CompositeGrid   *cg = (CompositeGrid*)wvs[0]->cg_over;

        num_patches = wvs[0]->totalNumberOfPatches;

        for (int1 =num_patches-1; int1 >= 0 ; int1--)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;

            if(waveto->patch_level > patch_level0) continue;
            if(patch_level0 - waveto->patch_level > 1) continue;

            if ( coords0[0] >= gridto->L[0] and coords0[0] <= gridto->U[0] and
                 coords0[1] >= gridto->L[1] and coords0[1] <= gridto->U[1])
            {
                *on_patch = waveto->patch_number;
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in pt_in_others_same_or_coarse_level(), "
                           "rect_in_which() failed for next coarse grid\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                return 0;
            }
        }
        for (int1 =num_patches-1; int1 >= 0 ; int1--)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;

            if(waveto->patch_level > patch_level0) continue;
            if(patch_level0 - waveto->patch_level > 1) continue;

            if ( coords0[0] >= gridto->VL[0] and coords0[0] <= gridto->VU[0] and
                 coords0[1] >= gridto->VL[1] and coords0[1] <= gridto->VU[1])
            {
                *on_patch = waveto->patch_number;
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in pt_in_others_same_or_coarse_level(), "
                           "rect_in_which() failed for next coarse grid\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                return 1;
            }
        }

        printf("ERROR in pt_in_others_same_or_coarse_level\n");
        printf("Looking for coords<%g, %g>ic<%d,%d> on_level[%d]\n",
                coords0[0], coords0[1],icoords0[0],icoords0[1],patch_level0);
        printf("looped all the patches, but not found\n");
        for (int1 =num_patches-1; int1 >= 0 ; int1--)
        {
            printf("Wave[%d] RECT_GRID\n",int1); 
            print_rectangular_grid(wvs[int1]->rect_grid);   
        } 
        clean_up(ERROR);
        return NO;
}

/*       best_match_pt_in_other_sets()   
*  This function call is used for searching a 
*  patch which contains the point.
*  The patches searched and the patch contains the 
*  point are from different sets. So the searching
*  will not find the best match point to be the
*  point itself. 
*
*  The search is performed as the following:
*  1. Search the interior of all patches.
*  2. Search the buffer zone of all patches.   
*
*  The value "on_patch" tells which patch best matches 
*  the point. "where" tells if the best match point is in the
*  patch interior or the patch buffer zone.  The return value
*  tells the level of the best match patch.  
*/  
LOCAL int best_match_pt_in_other_sets(
        int             *icoords0,
        double           *coords0,  
        int             comp0,
        int             patch_level0,  
        int             *on_patch,
        int             *where,  
        Wave            **wvs,
        Front           **frs,
        int             allow_buf)
{
        Wave            *waveto;
        Front           *frontto;
        RECT_GRID       *gridto;
        double           coords1[MAXD];
        double           dh[MAXD];  
        int             icoords1[MAXD];
        int             num_patches;
        double           VL[MAXD], VU[MAXD]; 
        int             i;

        num_patches = wvs[0]->totalNumberOfPatches;
        for (i =num_patches-1; i >= 0 ; i--)
        {
            waveto = wvs[i];
            frontto = frs[i];
            /* 043103 add this level condition */  
            if(patch_level0 < waveto->patch_level) continue; 

            gridto = waveto->rect_grid; 
            for(int dir = 0; dir < gridto->dim; dir++)
                dh[dir] = gridto->h[dir];  
            if ( coords0[0] >= gridto->L[0] and coords0[0] <= gridto->U[0] and
                 coords0[1] >= gridto->L[1] and coords0[1] <= gridto->U[1])
            {
                *on_patch = waveto->patch_number;
                *where = INTERIOR;  
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in best_match_pt_in_other_sets(), "
                           "rect_in_which() failed for grid %d, level[%d]\n",
                          i, waveto->patch_level);
                    printf("Looking for pt<%g, %g>ic<%d,%d> comp[%d]in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1],comp0);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                return waveto->patch_level;
            }
        }

        if(YES == allow_buf)
        {
            for (i =num_patches-1; i >= 0 ; i--)
            {
                waveto = wvs[i];
                frontto = frs[i];
                /* 043103 add this level condition */  
                if(patch_level0 < waveto->patch_level) continue; 

                gridto = waveto->rect_grid;
                for(int dir = 0; dir < gridto->dim; dir++)
                    dh[dir] = gridto->h[dir];  
                for(int dir = 0; dir < gridto->dim; dir++)
                {
                    if(gridto->lbuf[dir] > 0)
                        VL[dir] = gridto->VL[dir]+dh[dir]/2.0;
                    else
                        VL[dir] = gridto->VL[dir];
                    if(gridto->ubuf[dir] > 0)
                        VU[dir] = gridto->VU[dir]-dh[dir]/2.0;
                    else
                        VU[dir] = gridto->VU[dir];
                }

                if ( coords0[0] >= VL[0] and
                     coords0[0] <= VU[0] and
                     coords0[1] >= VL[1] and
                     coords0[1] <= VU[1])
                {
                    *on_patch = waveto->patch_number;
                    *where = BUFFER_ZONE;  
                    if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                    {
                        printf("ERROR in best_match_pt_in_other_sets(), "
                             "rect_in_which() failed for buffer grid %d, level[%d]\n",
                             i, waveto->patch_level);
                        printf("Looking for pt<%g, %g>ic<%d,%d> comp[%d]in rect_grid:\n",
                            coords0[0], coords0[1],icoords0[0],icoords0[1],comp0);
                        print_rectangular_grid(gridto);
                        clean_up(ERROR);
                    }
                    return waveto->patch_level;
                }
            }
        }

        printf("ERROR in best_match_pt_in_other_sets()\n");
        printf("Looking for coords<%g, %g>ic<%d,%d> on_level[%d]\n",
                coords0[0], coords0[1],icoords0[0],icoords0[1],patch_level0);
        printf("looped all the patches, but not found\n");
        for (i =num_patches-1; i >= 0 ; i--)
        {
            printf("Wave[%d] RECT_GRID\n",i); 
            print_rectangular_grid(wvs[i]->rect_grid);   
        } 
        clean_up(ERROR);
        return -1;
}



LOCAL  Locstate fine_to_fine(
        int             icoords0[MAXD],
        double           coords0[MAXD],
        COMPONENT       comp0,
        Wave            *wave,
        RECT_GRID       *grid,
        Front           *front)
{
        static Locstate st = NULL;
        size_t          sizest = front->sizest;
        Locstate        st1;
        COMPONENT       comp1;
        double           coords1[MAXD];
        int             icoords1[MAXD], ic[MAXD];
        int             dim = wave->rect_grid->dim;
        int             i;
        RECT_GRID       *comp_g = &wave_tri_soln(wave)->tri_grid->comp_grid;
        RECT_GRID       *tg = &wave_tri_soln(wave)->tri_grid->tg_grid;
        RECT_GRID       *rect_g = &wave_tri_soln(wave)->tri_grid->rect_grid;


        if(st == NULL)
            alloc_state(front->interf,&st,max((int)sizeof(VGas),(int)sizest));


        switch(dim)
        {
#if defined(ONED)
            case 1:
            {
                int     ix;
                int     xmax;

                xmax = grid->gmax[0];
                for (ix = 0; ix < xmax; ix++)
                {
                    coords1[0] = grid->L[0] + ix*(grid->h[0]);
                    icoords1[0] = ix;

                    if ( coords1[0] >= coords0[0] )
                    {
                        if ( ix == 0)
                        {
                            icoords1[0] =  icoords1[0] + 1;
                            coords1[0] = coords1[0] + grid->h[0];
                        }
                        break;
                    }
                }
                break;
            }
#endif /* defined(ONED) */
#if defined(TWOD)
            case 2:
            {
                if(rect_in_which(coords0,icoords1,comp_g) == FUNCTION_FAILED)
                {
                    printf("ERROR in fine_to_fine(), "
                                  "rect_in_which() failed\n");
                    printf("Looking for coords<%g, %g> in tri->comp_grid:\n",
                            coords0[0], coords0[1]);
                    print_rectangular_grid(comp_g);
                    clean_up(ERROR);
                }

                if(!same_point2d(Rect_coords(icoords1,wave), coords0, comp_g->h[0]))
                {
                    printf("ERROR in in fine_to_fine(), not same pt\n");
                    printf("crds<%g, %g> icrds found from rect_in_which is <%d, %d><%g,%g>\n",
                      coords0[0],coords0[1], icoords1[0], icoords1[1],
                      Rect_coords(icoords1,wave)[0],Rect_coords(icoords1,wave)[1]);
                    print_rectangular_grid(comp_g);
                    rect_in_which(coords0,icoords1, comp_g);
                    printf("crds<%g, %g> icrds found from rect_in_which"
                            " for comp_g is <%d, %d><%g,%g>\n",
                      coords0[0],coords0[1], icoords1[0], icoords1[1],
                      Rect_coords(icoords1,wave)[0],Rect_coords(icoords1,wave)[1]);
                    clean_up(ERROR);
                }
                break;
            }
#endif /* defined(TWOD) */
#if defined(THREED)
            /* 3d has not being tested yet */
            case 3:
            {
                int    ix,  iy,  iz;
                int    xmax, ymax, zmax;

                xmax = grid->gmax[0];
                ymax = grid->gmax[1];
                zmax = grid->gmax[2];

                for (iz = 0; iz < zmax; iz++)
                {
                    coords1[2] = grid->L[2] + iz*(grid->h[2]);
                    icoords1[2] = iz;

                    if ( coords1[2] >= coords0[2] )
                    {
                        if ( iz == 0)
                        {
                            icoords1[2] =  icoords1[2] + 1;
                            coords1[2] = coords1[2] + grid->h[2];
                        }
                        break;
                    }
                }

                for (iy = 0; iy < ymax; iy++)
                {
                    coords1[1] = grid->L[1] + iy*(grid->h[1]);
                    icoords1[1] = iy;

                    if ( coords1[1] >= coords0[1] )
                    {
                        if ( iy == 0)
                        {
                            icoords1[1] =  icoords1[1] + 1;
                            coords1[1] = coords1[1] + grid->h[1];
                        }
                        break;
                    }
                }

                for (ix = 0; ix < xmax; ix++)
                {
                    coords1[0] = grid->L[0] + ix*(grid->h[0]);
                    icoords1[0] = ix;

                    if ( coords1[0] >= coords0[0] )
                    {
                        if ( ix == 0)
                        {
                            icoords1[0] =  icoords1[0] + 1;
                            coords1[0] = coords1[0] + grid->h[0];
                        }
                        break;
                    }
                }
                break;
            }
#endif /* defined(THREED) */
        }  /* end of switch */

        ft_assign(st,Rect_state(icoords1,wave),front->sizest);  
        return st;
}  

LOCAL int amr_buffer_zone_block(
        int          *L,
        int          *U,
        int          pos,
        RECT_GRID    *gr,
        INTERFACE    *intfc) 
{
        int          i, dim = gr->dim; 
        int          lbuf[MAXD], ubuf[MAXD];  

        for (i = 0; i < dim; i++)
        {
            lbuf[i] = gr->lbuf[i];
            ubuf[i] = gr->ubuf[i];
        }
        switch(pos)
        {
        case 0:
            L[0] = -lbuf[0];  L[1] = -lbuf[1]; 
            U[0] = 0;         U[1] = 0; 
            if(rect_boundary_type(intfc,0,0) == AMR_SUBDOMAIN_BOUNDARY
               or rect_boundary_type(intfc,1,0) == AMR_SUBDOMAIN_BOUNDARY) 
                return YES; 
        break; 
        case 1:
            L[0] = 0;         L[1] = -lbuf[1]; 
            U[0] = gr->gmax[0];  U[1] = 0; 
            if(rect_boundary_type(intfc,1,0) == AMR_SUBDOMAIN_BOUNDARY) 
                return YES; 
        break; 
        case 2:
            L[0] = gr->gmax[0];  L[1] = -lbuf[1]; 
            U[0] = gr->gmax[0]+ubuf[0]; U[1] = 0; 
            if(rect_boundary_type(intfc,0,1) == AMR_SUBDOMAIN_BOUNDARY
               or rect_boundary_type(intfc,1,0) == AMR_SUBDOMAIN_BOUNDARY) 
                return YES; 
        break;  
        case 3:
            L[0] = gr->gmax[0];  L[1] = 0;
            U[0] = gr->gmax[0]+ubuf[0]; U[1] = gr->gmax[1];
            if(rect_boundary_type(intfc,0,1) == AMR_SUBDOMAIN_BOUNDARY) 
                return YES; 
        break; 
        case 4: 
            L[0] = gr->gmax[0];  L[1] = gr->gmax[1];
            U[0] = gr->gmax[0]+ubuf[0]; U[1] = gr->gmax[1]+ubuf[1];  
            if(rect_boundary_type(intfc,0,1) == AMR_SUBDOMAIN_BOUNDARY
               or rect_boundary_type(intfc,1,1) == AMR_SUBDOMAIN_BOUNDARY) 
                return YES; 
        break;  
        case 5:
            L[0] = 0;            L[1] = gr->gmax[1];  
            U[0] = gr->gmax[0];  U[1] = gr->gmax[1]+ubuf[1]; 
            if(rect_boundary_type(intfc,1,1) == AMR_SUBDOMAIN_BOUNDARY)  
                return YES; 
        break; 
        case 6:
            L[0] = -lbuf[0];     L[1] = gr->gmax[1];
            U[0] = 0;            U[1] = gr->gmax[1]+ubuf[1];
            if(rect_boundary_type(intfc,0,0) == AMR_SUBDOMAIN_BOUNDARY
               or rect_boundary_type(intfc,1,1) == AMR_SUBDOMAIN_BOUNDARY) 
                return YES; 
        break;  
        case 7:
            L[0] = -lbuf[0];     L[1] = 0;
            U[0] = 0;            U[1] = gr->gmax[1];
            if(rect_boundary_type(intfc,0,0) == AMR_SUBDOMAIN_BOUNDARY)
                return YES; 
        break;  
        default:
            printf("ERROR: amr_buffer_zone_block\n");
            printf("Pos[%d] unknown block position\n", pos); 
            clean_up(ERROR);  
        break; 
        }
        return NO;  
} 

/*
*               g_overture_fill_patch_amr_buffer_pt2d()
*
*     if a fine grid buffer zone is not a physical boundary,
*     interpolate this buffer zone state from coarse grid
*     or from the grid at the same refinement level.
* scatter_states() is called before overture_fill_patch_amr_buffer.  
* All FronTier subdomain, physical boundaries are filled
*  by scatter_states().  
*/  

EXPORT boolean g_overture_fill_patch_amr_buffer_pt2d(
        int             *ic,
        int             wv_id,  
        Wave            **wvs,
        Front           **frs)
{
        Wave            *wavefrom, *waveto;
        Front           *frontfrom, *frontto;
        RECT_GRID       *gridfrom, *gridto;
        Locstate        st;
        double           coords0[MAXD], coords1[MAXD];
        int             icoords0[MAXD], icoords1[MAXD];
        int             dim = wvs[0]->rect_grid->dim;
        int             num_patches, patch_num0, patch_num1;
        int             on_patch;
        int             patch_level0, patch_level1;
        COMPONENT       comp0;
        doubleCompositeGridFunction *uu1;
        CompositeGrid   *cg = (CompositeGrid*)wvs[0]->cg_over;
        size_t          sizest = wvs[0]->sizest;
        int             num_float;
        INTERFACE       *intfc; 
        int             where;
        /* 
        DEBUG_ENTER(g_overture_fill_patch_amr_buffer_pt2d)
        */  
        if(wvs[0]->use_overture_state)
            uu1 = (doubleCompositeGridFunction*)wvs[0]->cg_over_function; 

        num_patches = wvs[0]->totalNumberOfPatches;
        waveto = wvs[wv_id];
        frontto = frs[wv_id];
        gridto = waveto->rect_grid;
        intfc = frs[wv_id]->interf; 
        patch_num0 = waveto->patch_number;
        patch_level0 = waveto->patch_level;
            
        st = Rect_state(ic,waveto);
        coords0[0] = Rect_coords(ic,waveto)[0];
        coords0[1] = Rect_coords(ic,waveto)[1];
        comp0 = Rect_comp(ic,waveto);
        where = pt_in_patch_same_or_coarse_level(ic,coords0,comp0,
                        patch_num0,&on_patch, wvs,frs);
        if(0 == where)
        {
            /* pt interior (same level or next coarse grid) */ 
            wavefrom = wvs[on_patch]; frontfrom = frs[on_patch];
            gridfrom = waveto->rect_grid;
            if(wvs[on_patch]->patch_level == patch_level0)
            {
                /* fine to fine */  
                ft_assign(st, fine_to_fine(icoords0,coords0,comp0,
                      wavefrom,gridfrom,frontfrom),frontfrom->sizest); 
            }
            else
            {
                /* fine from coarse */  
                if(patch_level0 - wvs[on_patch]->patch_level != 1)
                {
                    printf("ERROR in g_overture_fill_patch_non_phy_buffer_pt2d()\n");
                    printf("levels are not right\n");
                    clean_up(ERROR);
                }
                /* 
                ft_assign(st, init_coarse_to_fine(icoords0,coords0,comp0,
                       wavefrom,gridfrom,frontfrom), frontfrom->sizest);
                */  
                hyp_solution(coords0,comp0,NULL,UNKNOWN_SIDE,
                       frontfrom,wavefrom,st,NULL); 
            }
        }
        else if(1 == where)
        {
            /* pt in the next coarse grid, buffer zone */ 
            wavefrom = wvs[on_patch]; frontfrom = frs[on_patch];
            gridfrom = waveto->rect_grid;
            /* 
            ft_assign(st, init_coarse_to_fine(icoords0,coords0,comp0,
                      wavefrom,gridfrom,frontfrom), frontfrom->sizest);
            */ 
            hyp_solution(coords0,comp0,NULL,UNKNOWN_SIDE,
                       frontfrom,wavefrom,st,NULL); 
            /*  This is a very temp solution, need to be fixed */  
            /* 042603 */ 
            /* 
            hyp_solution(coords0,comp0,NULL,UNKNOWN_SIDE,
                       frs[0],wvs[0],st,NULL); 
            */  
        }
        else if(2 == where)
        {
            /* pt in overture buffer zone */  
            int ic[2] = {0,0};
            realArray pos(1,3), uval(1,8);
            Index     I1, I2, I3;
            pos(0,0) = coords0[0]; pos(0,1) = coords0[1];
            pos(0,2) = 0.0;

            printf("ERROR in  g_overture_fill_patch_non_phy_buffer_pt2d()\n");
            printf("Need to consider whether overture_state is used\n");
            clean_up(ERROR);  
        }
        /* 
        DEBUG_LEAVE(g_overture_fill_patch_amr_buffer_pt2d)
        */    
        return FUNCTION_SUCCEEDED;  
}

/*
*               g_overture_fill_patch_amr_buffer2d()
*
*     if a fine grid buffer zone is not a physical boundary,
*     interpolate this buffer zone state from coarse grid
*     or from the grid at the same refinement level.
*
* scatter_states() is called before overture_fill_patch_amr_buffer.  
* All FronTier subdomain, physical boundaries are filled
*  by scatter_states().  
*/  

EXPORT boolean g_overture_fill_patch_amr_buffer2d(
        Wave            **wvs,
        Front           **frs)
{
        Wave            *wavefrom, *waveto;
        Front           *frontfrom, *frontto;
        RECT_GRID       *gridfrom, *gridto;
        Locstate        st;
        double           coords0[MAXD], coords1[MAXD];
        int             icoords0[MAXD], icoords1[MAXD];
        int             dim = wvs[0]->rect_grid->dim;
        int             num_patches, patch_num0, patch_num1;
        int             on_patch;
        int             levels, patch_level0, patch_level1;
        int             i, int1;
        int             int2;
        COMPONENT       comp0,comp1;
        int             lbuf[MAXD], ubuf[MAXD];
        doubleCompositeGridFunction *uu1;
        CompositeGrid   *cg = (CompositeGrid*)wvs[0]->cg_over;
        size_t          sizest = wvs[0]->sizest;
        int             num_float;
        INTERFACE       *intfc; 
        int             L[MAXD], U[MAXD];  
        int             amr_buf; 
        int             swp, side;  
        int             ix, iy;
        int             where;
        int             *iperm; 

        DEBUG_ENTER(g_overture_fill_patch_amr_buffer2d)

        if(dim != 2)
        {
            printf("ERROR: g_overture_fill_patch_amr_buffer\n");
            printf("Dim[%d] case is not implemented\n");
            clean_up(ERROR); 
        } 

        num_patches = wvs[0]->totalNumberOfPatches;
        levels = wvs[0]->NumberOfLevels;

        if(YES == wvs[0]->use_overture_state)
            fill_root_extr_overture_cell_st(wvs,frs);

        iperm = set_iperm(frs[0]->step,dim); 
        for(int1 = 1; int1 < num_patches; int1++) 
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;
            intfc = frs[int1]->interf; 
            patch_num0 = waveto->patch_number;
            patch_level0 = waveto->patch_level;
            
            for(swp = 0; swp < dim; swp++)
            {
                for(side = 0; side < 2; side++)  
                {
                    set_receive_domain(L,U,iperm,side,swp,frs[int1]->rect_grid); 
                    /* 
                    printf("POS[%d], L[%d, %d], U[%d, %d], AMR_BUFFER = %s\n",
                       i, L[0], L[1], U[0], U[1], amr_buf == YES ? "yes": "no"); 
                    */ 
                    for (iy = L[1]; iy < U[1]; iy++)
                    {
                        icoords0[1] = iy;
                        for (ix = L[0]; ix < U[0]; ix++)
                        {
                            icoords0[0] = ix; 
                            st = Rect_state(icoords0,waveto); 
                            if(is_obstacle_state(st))
                                g_overture_fill_patch_amr_buffer_pt2d(icoords0, 
                                   patch_num0, wvs, frs); 
                        }
                    }
                }  
            }  
        }
/*         
        for(int1 = 1; int1 < num_patches; int1++) 
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;
            intfc = frs[int1]->interf; 
            patch_num0 = waveto->patch_number;
            patch_level0 = waveto->patch_level;
            
            for(int pos = 0; pos < 8; pos++)
            {
                amr_buf = amr_buffer_zone_block(L,U,pos,gridto,intfc); 
                if(YES == amr_buf)
                {
                    int      ix, iy;
                    int      where;
                    for (iy = L[1]; iy < U[1]; iy++)
                    {
                        icoords0[1] = iy;
                        for (ix = L[0]; ix < U[0]; ix++)
                        {
                            icoords0[0] = ix;
                            g_overture_fill_patch_amr_buffer_pt2d(icoords0, 
                                patch_num0, wvs, frs); 
                        }
                    }
                }  
            }  
        }
*/       
/*      fill_extr_overture_upper_cell_st(*cg, *uu1); */  
{
        for (int1 = num_patches-1; int1 >= 0; int1--)
        {
            int      ic[2] = {0,0};
            int      smin[MAXD], smax[MAXD];
            int      ix, iy;
            Index    I1, I2, I3;
            int      base1,  base2, base3, bound1, bound2, bound3;

            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;
            getIndex((*cg)[int1].indexRange(),I1,I2,I3);

            base1 = I1.getBase();
            base2 = I2.getBase();
            base3 = I3.getBase();
            bound1 = I1.getBound();
            bound2 = I2.getBound();
            bound3 = I3.getBound();

            for (i = 0; i < dim; i++)
            {
                smin[i] = 0;
                smax[i] = gridto->gmax[i];
                lbuf[i] = gridto->lbuf[i];
                ubuf[i] = gridto->ubuf[i];
                smin[i] -= lbuf[i];
                smax[i] += ubuf[i];
            }
            for (iy = smin[1]; iy < smax[1]; iy++)
            {
                icoords0[1] = iy;
                for (ix = smin[0]; ix < smax[0]; ix++)
                {
                    icoords0[0] = ix;
                    st = Rect_state(icoords0,waveto);
                    coords0[0] = Rect_coords(icoords0,waveto)[0];
                    coords0[1] = Rect_coords(icoords0,waveto)[1];
                    comp0 = Rect_comp(icoords0,waveto);
                    if(Different_params(st,Rect_state(ic,waveto)))
                    {
                        screen("ERROR: g_overture_fill_patch_amr_buffer2d\n");
                        screen("crds<%g,%g> IC<%d,%d> COMP[%d] patch[%d] is diff params\n",
                            coords0[0], coords0[1], icoords0[0], icoords0[1], comp0, int1);
                        screen("crds<%g,%g> subdomain IC<%d,%d>\n",
                            coords0[0], coords0[1], icoords0[0]+base1, icoords0[1]+base2);
                        printf("print <0,0> patch[%d] st\n", waveto->patch_number);
                        for(i = 0; i < dim; i++)
                        {
                            for(int side = 0; side < 2; side++)
                                printf("rect_boundary_type(dir[%d],side[%d]) = %d\n",
                                  i, side, rect_boundary_type(frontto->interf,i,side)); 
                        }  
                        print_rectangular_grid(waveto->rect_grid);
                        g_verbose_print_state(Rect_state(ic,waveto));
                        printf("beforedens<%g>, mom<%g,%g>, eng<%g> sttype<%d>\n",
                           Dens(st), Mom(st)[0], Mom(st)[1], Energy(st), state_type(st));
                        g_verbose_print_state(st);
                        clean_up(ERROR);
                    }
                    if(is_obstacle_state(st))
                    {
                        screen("exit g_overture_fill_patch_amr_buffer2d\n");
                        screen("crds<%g,%g> ic<%d,%d> patch[%d] is obstacle\n",
                            coords0[0], coords0[1], icoords0[0], icoords0[1], int1);
                        print_rectangular_grid(waveto->rect_grid);  
                        for(i = 0; i < dim; i++)
                        {
                            for(int side = 0; side < 2; side++)
                                printf("rect_boundary_type(dir[%d],side[%d]) = %d\n",
                                  i, side, rect_boundary_type(frontto->interf,i,side)); 
                        }  
                        printf("patch[%d][%d]\n",int1, waveto->patch_number);
                        g_verbose_print_state(st);
                        clean_up(ERROR);
                    }
                }
            }
        }
}
        DEBUG_LEAVE(g_overture_fill_patch_amr_buffer2d)
        return FUNCTION_SUCCEEDED;  
}

EXPORT  boolean g_scatter_patch_states_in_sweep_dir(
        Overparam    *overparam,
        Wv_on_pc     **redistr_table,
        Wave         **wvs,
        Front        **frs,
        int          *iperm,
        int          swp)
{
        int          num_patches, levels;
        int          i, dim, side;
        int          imin[MAXD], imax[MAXD];
        int          L[MAXD], U[MAXD];
        RECT_GRID    *gr, *gridfrom;
        int          ix, iy, ic[MAXD];
        int          patch_num0;
        Wave         *wave0, *wavefrom;
        Front        *frontfrom;
        Locstate     st;
        int          on_patch, where, where_level;  

        DEBUG_ENTER(g_scatter_patch_states_in_sweep_dir)

        num_patches = wvs[0]->totalNumberOfPatches;
        levels = wvs[0]->NumberOfLevels;
        gr = wvs[0]->rect_grid;
        dim = gr->dim;

        for (i = 0; i < swp; ++i)
        {
            imin[iperm[i]] = -gr->lbuf[iperm[i]];
            imax[iperm[i]] = gr->gmax[iperm[i]] + gr->ubuf[iperm[i]];
        }
        imin[iperm[swp]] = 0;
        imax[iperm[swp]] = gr->gmax[iperm[swp]];

        for (i = swp+1; i < dim; ++i)
        {
            imin[iperm[i]] = -gr->lbuf[iperm[i]];
            imax[iperm[i]] = gr->gmax[iperm[i]] + gr->ubuf[iperm[i]];
        }

        if(debugging("g_scatter_patch_states_in_sweep_dir"))
            printf("Swp[%d], dir[%d] patch[0] imin[%d, %d], imax[%d, %d]\n",
              swp, iperm[swp], imin[0], imin[1], imax[0], imax[1]);
        for(i = 1; i < num_patches; i++)
        {
            
            if(debugging("g_scatter_patch_states_in_sweep_dir"))
            {
                printf("patch[%d], level[%d], rect_grid\n", 
                       i, wvs[i]->patch_level);
                print_rectangular_grid(wvs[i]->rect_grid);
            }  
            for(side =0; side < 2; side++)
            {
                set_patch_domain_in_sweep_dir(L,U,imin,imax,
                  iperm,swp,side,overparam,redistr_table,frs[i]);
                if(rect_boundary_type(frs[i]->interf,iperm[swp],side) !=
                       AMR_SUBDOMAIN_BOUNDARY &&
                   rect_boundary_type(frs[i]->interf,iperm[swp],side) !=
                       SUBDOMAIN_BOUNDARY)
                    continue;
                if(debugging("g_scatter_patch_states_in_sweep_dir"))
                {
                    printf("Side[%d]:Bdry_type[%d] amr_buf: L[%d,%d], U[%d,%d]\n",
                     side, rect_boundary_type(frs[i]->interf,iperm[swp],side),
                     L[0], L[1], U[0], U[1]);
                }
                if(rect_boundary_type(frs[i]->interf,iperm[swp],side) ==
                     SUBDOMAIN_BOUNDARY)
                {
                    for(iy = L[1]; iy < U[1]; iy++)
                    {
                        for(ix = L[0]; ix < U[0]; ix++)
                        {
                            ic[1] = iy; ic[0] = ix;
                            st = Rect_state(ic,wvs[i]);
                            if(NULL == Params(st))
                            {
                                best_match_pt_in_coarse(ic,Rect_coords(ic,wvs[i]),
                                 wvs[i]->patch_level,&on_patch,&where,wvs,frs,YES); 
                                wavefrom = wvs[on_patch]; 
                                frontfrom = frs[on_patch];
                                hyp_solution(Rect_coords(ic,wvs[i]),
					Rect_comp(ic,wvs[i]),NULL,UNKNOWN_SIDE,
					frontfrom,wavefrom,st,NULL);
                            } 
                        }
                    }  
                }  
                else
                {
                    for(iy = L[1]; iy < U[1]; iy++)
                    {
                        for(ix = L[0]; ix < U[0]; ix++)
                        {
                            ic[1] = iy; ic[0] = ix;
                            st = Rect_state(ic,wvs[i]);
                            where_level = best_match_pt_in_same_set(
                               ic, Rect_coords(ic,wvs[i]), Rect_comp(ic,wvs[i]),
                              wvs[i]->patch_number,&on_patch,&where,wvs,frs); 
                            wavefrom = wvs[on_patch]; 
                            frontfrom = frs[on_patch];
                            gridfrom = wavefrom->rect_grid; 
                            if(wvs[i]->patch_level == where_level)
                            {
                                ft_assign(st, fine_to_fine(ic,Rect_coords(ic,wvs[i]),
                                    Rect_comp(ic,wvs[i]),wavefrom,gridfrom,frontfrom),
                                       frontfrom->sizest);
                            }  
                            else
                            {
                                if(!is_excluded_comp(Rect_comp(ic,wvs[i]),frs[i]->interf))
                                {
                                    hyp_solution(Rect_coords(ic,wvs[i]),
					Rect_comp(ic,wvs[i]),NULL,UNKNOWN_SIDE,
					frontfrom,wavefrom,st,NULL); 
                                }
                                else
                                {
                                    g_obstacle_state(st,frs[i]->sizest);
                                }  
                            }  
                        }
                    }
                }  
            }
        }

        DEBUG_LEAVE(g_scatter_patch_states_in_sweep_dir)
        return FUNCTION_SUCCEEDED;
}  

LOCAL  void set_patch_domain_in_sweep_dir(
        int          *L,
        int          *U,
        int          *imin,
        int          *imax,
        int          *iperm,
        int          swp,
        int          side,
        Overparam    *overparam,
        Wv_on_pc     **redistr_table,
        Front        *front)
{
        int          myid;
        int          i, j, dim;
        RECT_GRID    *gr;
        int          refine;
        int          ggmin[MAXD], ggmax[MAXD];
        int          patch_number, level;
        int          base[MAXD], bound[MAXD];
        int          tmpL[MAXD], tmpU[MAXD];
        INTERFACE    *intfc = front->interf;

        myid = pp_mynode();
        gr = front->rect_grid;
        dim = gr->dim;
        if(rect_boundary_type(intfc,iperm[swp],side) !=
            AMR_SUBDOMAIN_BOUNDARY && 
           rect_boundary_type(intfc,iperm[swp],side) !=
            SUBDOMAIN_BOUNDARY)
        {
            for(j = 0; j < dim; j++)
            {
                L[j] = U[j] = 0;
            }
            return;
        }
        if(rect_boundary_type(intfc,iperm[swp],side) ==
            SUBDOMAIN_BOUNDARY) 
        {
            set_receive_domain(L,U,iperm,side,swp,gr);
            return;  
        } 
     
        refine = 1;
        level = front->patch_level;
        patch_number = front->patch_number;
        for(j = 0; j < level; j++)
            refine *= overparam->refinementRatio;
        set_receive_domain(tmpL,tmpU,iperm,side,swp,gr);

        for(j = 0; j < dim; j++)
        {
            ggmin[j] = refine*imin[j];
            ggmax[j] = refine*imax[j];
            base[j] = redistr_table[myid][patch_number].base[j];
            bound[j] = redistr_table[myid][patch_number].bound[j];
            tmpL[j] += base[j];
            tmpU[j] += base[j];
        }

        for(j = 0; j < dim; j++)
        {
            L[j] = max(tmpL[j],ggmin[j])-base[j];
            U[j] = min(tmpU[j],ggmax[j])-base[j];
        }
} 


EXPORT void g_overture_assign_wave_params(
        Locstate  st1,
        Locstate  st2)
{
       Set_params(st1,st2);
}

EXPORT void g_overture_assign_wave_st_type(
        Locstate        st1,
        Locstate        st2)
{
        set_type_of_state(st1,state_type(st2));
}


LOCAL void fill_extr_overture_upper_cell_st(
       const CompositeGrid    &cg,
       doubleCompositeGridFunction  &u)
{
        Index      I1,I2,I3;
        int        ix, iy;
        int        base1,base2,bound1,bound2;

        DEBUG_ENTER(fill_extr_overture_upper_cell_st)

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++)
        {
            getIndex(cg[grid].dimension(),I1,I2,I3);
            base1 = I1.getBase();
            base2 = I2.getBase();
            bound1 = I1.getBound();
            bound2 = I2.getBound();
            for (ix = base1; ix <= bound1; ix++)
            {
                u[grid](0,ix,bound2,0) = u[grid](0,ix,bound2-1,0);
                u[grid](1,ix,bound2,0) = u[grid](1,ix,bound2-1,0);
                u[grid](2,ix,bound2,0) = u[grid](2,ix,bound2-1,0);
                u[grid](3,ix,bound2,0) = u[grid](3,ix,bound2-1,0);
            }
            for (iy = base2; iy <= bound2; iy++)
            {
                u[grid](0,bound1,iy,0) = u[grid](0,bound1-1,iy,0);
                u[grid](1,bound1,iy,0) = u[grid](1,bound1-1,iy,0);
                u[grid](2,bound1,iy,0) = u[grid](2,bound1-1,iy,0);
                u[grid](3,bound1,iy,0) = u[grid](3,bound1-1,iy,0);
            }
        }

        DEBUG_LEAVE(Leaving fill_extr_overture_upper_cell_st)
}

EXPORT void g_fill_root_extr_overture_cell_st_from_wv(
	Wave       *wv,
	Front      *fr,
	POINTER    cg_func,
	int        grid,  
	int        var)
{
        CompositeGrid *cg = (CompositeGrid*)wv->cg_over;
        /* 
	doubleCompositeGridFunction& uu = *((doubleCompositeGridFunction*)cg_func); 
        */  
	doubleCompositeGridFunction *uu = (doubleCompositeGridFunction*)cg_func; 
        Index      I1,I2,I3;

        INTERFACE  *rt_interf = fr->interf;
        int        ix, iy;
        int        xmin, xmax, ymin, ymax;
        int        base1,base2,bound1,bound2;
        int        base11,base22,bound11,bound22;

	DEBUG_ENTER(fill_extr_overture_cell_st_from_wv) 

        getIndex((*cg)[grid].indexRange(),I1,I2,I3);
        base11 = I1.getBase();
        base22 = I2.getBase();
        bound11 = I1.getBound();
        bound22 = I2.getBound();

        getIndex((*cg)[grid].dimension(),I1,I2,I3);
        base1 = I1.getBase();
        base2 = I2.getBase();
        bound1 = I1.getBound();
        bound2 = I2.getBound();

        /* 
        if(grid == 2)
        {
            printf("In g_fill_root_extr_overture_cell_st_from_wv\n");
            printf("print grid 2 density\n"); 
            uu[grid].display("density");
        }
        */

        for(int axis = 0; axis<cg->numberOfDimensions(); axis++)
        {
            for(int side = Start; side<=End; side++)
            {
                if(rect_boundary_type(rt_interf,axis,side) == SUBDOMAIN_BOUNDARY)
                    continue;
                if(rect_boundary_type(rt_interf,axis,side) == REFLECTION_BOUNDARY)
                    continue;
                if(rect_boundary_type(rt_interf,axis,side) == AMR_SUBDOMAIN_BOUNDARY)
                    continue;
                if(axis == 0 and side == 0)
                {
                    if((*cg)[grid].boundaryCondition()(0,1) == 0)
                        ymin = base22;
                    else
                        ymin = base2;
                    if((*cg)[grid].boundaryCondition()(1,1) == 0)
                        ymax = bound22;
                    else
                        ymax = bound2;
                    xmin=base1; xmax=base11;
                    for (iy = ymin; iy < ymax; iy++)
                    {
                        for (ix = xmin; ix < xmax; ix++)
                        {
                            (*uu)[grid](ix,iy,0) = (*uu)[grid](0,iy,0); 
                        }
                    }  
                }
                if(axis == 0 and side == 1)
                {
                    if((*cg)[grid].boundaryCondition()(0,1) == 0)
                        ymin = base22;
                    else
                        ymin = base2;
                    if((*cg)[grid].boundaryCondition()(1,1) == 0)
                        ymax = bound22;
                    else
                        ymax = bound2;
                    xmin=bound11; xmax=bound1;
                    for (iy = ymin; iy < ymax; iy++)
                    {
                        for (ix = xmin; ix < xmax; ix++)
                        {
                            (*uu)[grid](ix,iy,0) = (*uu)[grid](bound11-1,iy,0);
                        }
                    }
                }
                if(axis == 1 and side == 0)
                {
                    if((*cg)[grid].boundaryCondition()(0,0) == 0)
                        xmin = base11;
                    else
                        xmin = base1;
                    if((*cg)[grid].boundaryCondition()(1,0) == 0)
                        xmax = bound11;
                    else
                        xmax = bound1;
                    ymin=base2; ymax=base22;
                    for (iy = ymin; iy < ymax; iy++)
                    {
                        for (ix = xmin; ix < xmax; ix++)
                        {
                            (*uu)[grid](ix,iy,0) = (*uu)[grid](ix,0,0);
                        }
                    }
                }
                if(axis == 1 and side == 1)
                {
                    if((*cg)[grid].boundaryCondition()(0,0) == 0)
                        xmin = base11;
                    else
                        xmin = base1;
                    if((*cg)[grid].boundaryCondition()(1,0) == 0)
                        xmax = bound11;
                    else
                        xmax = bound1;
                    ymin=bound22; ymax=bound2;

                    for (iy = ymin; iy < ymax; iy++)
                    {
                        for (ix = xmin; ix < xmax; ix++)
                        {
                            (*uu)[grid](ix,iy,0) = (*uu)[grid](ix,bound22-1,0);
                        }
                    }
                }
            }
        }  
        for (ix = base1; ix <= bound1; ix++)
        {
            (*uu)[grid](ix,bound2,0) = (*uu)[grid](ix,bound2-1,0);
        }
        for (iy = base2; iy <= bound2; iy++)
        {
            (*uu)[grid](bound1,iy,0) = (*uu)[grid](bound1-1,iy,0);
        }

	DEBUG_LEAVE(fill_extr_overture_cell_st_from_wv) 
}  

EXPORT void g_trans_wv_st_to_overfunc(
	Wave		*wv,
	Front   	*fr,
	POINTER		cg_func,
	int             grid,  
	int             var)
{
        int 		base[MAXD], bound[MAXD];
        int 		gbase[MAXD], gbound[MAXD];
        int 		ix,iy;
	RECT_GRID       *rg = wv->rect_grid; 
        int             i, dim; 
	int             smin[MAXD], smax[MAXD]; 
        int             lbuf[MAXD], ubuf[MAXD];  
        int             cglbuf[MAXD], cgubuf[MAXD];  
	int             ic[MAXD];  
	Index           I1, I2, I3;  
   
        /* 
	doubleCompositeGridFunction& uu = 
                *((doubleCompositeGridFunction*)cg_func); 
        */
	doubleCompositeGridFunction *uu = 
                (doubleCompositeGridFunction*)cg_func; 
	CompositeGrid *cg = (CompositeGrid*)wv->cg_over;  

	DEBUG_ENTER(g_trans_wv_st_to_overfunc) 

	dim = rg->dim; 

        getIndex((*cg)[grid].dimension(),I1,I2,I3);
        gbase[0] = I1.getBase();
        gbase[1] = I2.getBase();
        gbound[0] = I1.getBound();
        gbound[1] = I2.getBound();
        
        getIndex((*cg)[grid].indexRange(),I1,I2,I3);
        base[0] = I1.getBase();
        base[1] = I2.getBase();
        bound[0] = I1.getBound();
        bound[1] = I2.getBound();

        for(i = 0; i < dim; i++)
        {
            cglbuf[i] = base[i]-gbase[i]; 
            cgubuf[i] = gbound[i]-bound[i]; 
        } 

        for(i = 0; i < dim; i++)
        {
            lbuf[i] = min(rg->lbuf[i], cglbuf[i]); 
            ubuf[i] = min(rg->ubuf[i], cgubuf[i]); 
        }  

	for(i = 0; i < dim; i++)
	{
	    smin[i] = 0;
            smax[i] = rg->gmax[i];	    
            smin[i] -= lbuf[i];  
            smax[i] += ubuf[i];  
	} 
        
	for (iy = smin[1];  iy < smax[1];  iy++)
	{
	    for (ix = smin[0];  ix < smax[0];  ix++) 	
            {
		ic[0] = ix; ic[1] = iy;      
		if(var == 0) 
                {
                    if(Params(Rect_state(ic,wv)) == NULL)
                        (*uu)[grid](ix+base[0],iy+base[1],0) = 0.0;
                    else
	                (*uu)[grid](ix+base[0],iy+base[1],0) = 
                                   Dens(Rect_state(ic,wv));     	    
                }
		if(var == 1) 
                {
                    if(Params(Rect_state(ic,wv)) == NULL)
                        (*uu)[grid](ix+base[0],iy+base[1],0) = 0.0;
                    else
	                (*uu)[grid](ix+base[0],iy+base[1],0) = 
                                   Energy(Rect_state(ic,wv));     	    
                } 
		if(var == 2) 
                {
                    if(Params(Rect_state(ic,wv)) == NULL)
                        (*uu)[grid](ix+base[0],iy+base[1],0) = 0.0;
                    else
	                (*uu)[grid](ix+base[0],iy+base[1],0) = 
                                   Mom(Rect_state(ic,wv))[0];     	    
                }
		if(var == 3) 
                {
                    if(Params(Rect_state(ic,wv)) == NULL)
                        (*uu)[grid](ix+base[0],iy+base[1],0) = 0.0;
                    else
	                (*uu)[grid](ix+base[0],iy+base[1],0) = 
                                   Mom(Rect_state(ic,wv))[1];     	    
                }  
            }    		    
	} 	
	DEBUG_LEAVE(g_trans_wv_st_to_overfunc) 
}  

/*
*               g_overture_init_interpolation_coarse_to_fine()
*
*     initialize the states in patches produced by overture refinement
*/

EXPORT boolean g_overture_init_interpolation_coarse_to_fine(
        Wave            **wvs,
        Front           **frs)
{

        Wave            *waveto, *wavefrom;
        Front           *frontto, *frontfrom;
        RECT_GRID       *gridto, *gridfrom;
        Locstate        st;
        double           coords0[MAXD];
        int             icoords0[MAXD];
        int             dim = wvs[0]->rect_grid->dim;
        int             num_patches, patch_num0;
        int             levels, patch_level0;
        int             where_level, on_patch, where;  
        int             i, int1;
        COMPONENT       comp0;
        int             smin[MAXD], smax[MAXD];
        size_t          sizest = wvs[0]->sizest;
        boolean            stat = FUNCTION_SUCCEEDED;
        
        DEBUG_ENTER(g_overture_init_interpolation_coarse_to_fine)

        num_patches = wvs[0]->totalNumberOfPatches;
        levels = wvs[0]->NumberOfLevels;

        for (int1 = 1; int1 < num_patches ; int1++)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;
            for (i = 0; i < dim; i++)
            {
                smin[i] = 0;
                smax[i] = gridto->gmax[i];
                smin[i] -= gridto->lbuf[i];
                smax[i] += gridto->ubuf[i];
            }
            patch_num0 = waveto->patch_number;
            patch_level0 = waveto->patch_level;

            if(debugging("overture_coarse_to_fine"))
            {
                printf("in g_overture_init_interpolation_coarse_to_fine()"
                       "the working patch[%d], level[%d]\n", patch_num0, patch_level0);
                print_rectangular_grid(waveto->rect_grid);  
            } 
            switch(dim)
            {
#if defined(ONED)
            case 1:
            {
                printf("ERROR: g_overture_init_interpolation_coarse_to_fine()\n");
                printf("Case 1D is not implimented\n");
                clean_up(ERROR);  
                break;
            }
#endif /* defined(ONED) */
#if defined(TWOD)
            case 2:
            {
                int      ix, iy;
                for (iy = smin[1]; iy < smax[1]; iy++)
                {
                    icoords0[1] = iy;
                    for (ix = smin[0]; ix < smax[0]; ix++)
                    {
                        icoords0[0] = ix;
                        st = Rect_state(icoords0,waveto);
                        coords0[0] = Rect_coords(icoords0,waveto)[0];
                        coords0[1] = Rect_coords(icoords0,waveto)[1];
                        comp0 = Rect_comp(icoords0,waveto);
                        where_level = best_match_pt_in_coarse(icoords0,coords0,
                               waveto->patch_level,&on_patch,&where,wvs,frs,YES);
                        wavefrom = wvs[on_patch]; frontfrom = frs[on_patch];
                        gridfrom = wavefrom->rect_grid;  
                        /*
                        ft_assign(st, init_coarse_to_fine(icoords0,coords0,comp0,
                                 wavefrom,gridfrom,frontfrom),frontfrom->sizest);
                        */
                        hyp_solution(coords0,comp0,NULL,UNKNOWN_SIDE,
                                  frontfrom,wavefrom,st,NULL);

                    }
                }
                break;
            }
#endif /* defined(TWOD) */
#if defined(THREED)
            case 3:
                printf("ERROR: g_overture_init_interpolation_coarse_to_fine()\n");
                printf("Case 1D is not implimented\n");
                clean_up(ERROR);  
                break; 
#endif /* defined(THREED) */
            }  /* end of switch */
            if(debugging("overture_coarse_to_fine") and patch_num0 == num_patches-1)
            {
                printf("In g_overture_init_interpolation_coarse_to_fine()\n");
                printf("print wvs[%d] state\n",patch_num0);
                print_rectangular_grid(wvs[patch_num0]->rect_grid);
                (*wvs[0]->show_wave_states)(wvs[patch_num0]);
            }
        }
        DEBUG_LEAVE(g_overture_init_interpolation_coarse_to_fine)
        return stat;  
}

/* Pick a best match point from the coarser grid */ 
LOCAL  int best_match_pt_in_coarse(
	int             *icoords0,
        double           *coords0,
        int             patch_level, 
        int             *on_patch, 
        int             *where,
        Wave            **wvs,
        Front           **frs,
        int             allow_buf)
{
        Wave            *waveto;
        Front           *frontto;
        RECT_GRID       *gridto;
        double           coords1[MAXD];
        double           dh[MAXD];
        int             icoords1[MAXD];
        int             num_patches;
        double           VL[MAXD], VU[MAXD];
        int             i;

        num_patches = wvs[0]->totalNumberOfPatches;
        for (i =num_patches-1; i >= 0 ; i--)
        {
            waveto = wvs[i];
            frontto = frs[i];
            gridto = waveto->rect_grid;
            if(waveto->patch_level >= patch_level) continue;   

            for(int dir = 0; dir < gridto->dim; dir++)
                dh[dir] = gridto->h[dir];
            if ( coords0[0] >= gridto->L[0] and coords0[0] <= gridto->U[0] and
                 coords0[1] >= gridto->L[1] and coords0[1] <= gridto->U[1])
            {
                *on_patch = waveto->patch_number;
                *where = INTERIOR;
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in best_match_pt_in_coarse(), "
                           "rect_in_which() failed for grid %d, level[%d]\n",
                          i, waveto->patch_level);
                    printf("Looking for pt<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                return waveto->patch_level;
            }
            if(YES == allow_buf)
            {
                waveto = wvs[i];
                frontto = frs[i];
                gridto = waveto->rect_grid;

                for(int dir = 0; dir < gridto->dim; dir++)
                    dh[dir] = gridto->h[dir];
                for(int dir = 0; dir < gridto->dim; dir++)
                {
                    if(gridto->lbuf[dir] > 0)
                        VL[dir] = gridto->VL[dir]+dh[dir]/2.0;
                    else
                        VL[dir] = gridto->VL[dir];
                    if(gridto->ubuf[dir] > 0)
                        VU[dir] = gridto->VU[dir]-dh[dir]/2.0;
                    else
                        VU[dir] = gridto->VU[dir];
                }

                if ( coords0[0] >= VL[0] and
                     coords0[0] <= VU[0] and
                     coords0[1] >= VL[1] and
                     coords0[1] <= VU[1])
                {
                    *on_patch = waveto->patch_number;
                    *where = BUFFER_ZONE;
                    if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                    {
                        printf("ERROR in best_match_pt_in_coarse(), "
                             "rect_in_which() failed for buffer grid %d, level[%d]\n",
                             i, waveto->patch_level);
                        printf("Looking for pt<%g, %g>ic<%d,%d> in rect_grid:\n",
                            coords0[0], coords0[1],icoords0[0],icoords0[1]);
                        print_rectangular_grid(gridto);
                        clean_up(ERROR);
                    }
                    return waveto->patch_level;
                }
            }
        } 
        printf("ERROR in best_match_pt_in_coase()\n");
        printf("Looking for coords<%g, %g>ic<%d,%d> on_level[%d]\n",
                coords0[0], coords0[1],icoords0[0],icoords0[1],patch_level);
        printf("looped all the patches, but not found\n");
        for (i =num_patches-1; i >= 0 ; i--)
        {
            printf("Wave[%d] RECT_GRID\n",i);
            print_rectangular_grid(wvs[i]->rect_grid);
        }
        clean_up(ERROR);
        return -1;
} 

LOCAL  Locstate init_coarse_to_fine(
        int             icoords0[MAXD],
        double           coords0[MAXD],
        COMPONENT       comp0,
        Wave            *wave,
        RECT_GRID       *grid,
        Front           *front)
{
        static Locstate st = NULL;
        size_t          sizest = front->sizest;
        Locstate        st1;
        COMPONENT       comp1;
        double           coords1[MAXD];
        int             icoords1[MAXD], ic[MAXD];
        int             dim = wave->rect_grid->dim;
        int             i, point_num;
        RECT_GRID       *comp_g = &wave_tri_soln(wave)->tri_grid->comp_grid;
        RECT_GRID       *tg = &wave_tri_soln(wave)->tri_grid->tg_grid;
        RECT_GRID       *rect_g = wave->rect_grid;
        double           *h;
        double           dist= HUGE_VAL, tmp_dist;

        h = rect_g->h;
        if(st == NULL)
            alloc_state(front->interf,&st,max((int)sizeof(VGas),(int)sizest));

        switch(dim)
        {
#if defined(ONED)
            case 1:
            {
                int     ix;
                int     xmax;

                xmax = grid->gmax[0];
                for (ix = 0; ix < xmax; ix++)
                {
                    coords1[0] = grid->L[0] + ix*(grid->h[0]);
                    icoords1[0] = ix;

                    if ( coords1[0] >= coords0[0] )
                    {
                        if ( ix == 0)
                        {
                            icoords1[0] =  icoords1[0] + 1;
                            coords1[0] = coords1[0] + grid->h[0];
                        }
                        break;
                    }
                }
                break;
            }
#endif /* defined(ONED) */
#if defined(TWOD)
            case 2:
            {
                if(rect_in_which(coords0,icoords1,rect_g) == FUNCTION_FAILED)
                {
                    printf("ERROR in init_coarse_to_fine(), "
                                  "rect_in_which() failed\n");
                    clean_up(ERROR);
                }

#if !defined(TOL)
#define _USE_TOL_
  #define TOL        0.000001        /*TOLERANCE*/
#endif /* if !defined(TOL) */ 

                if(fabs(Rect_coords(icoords1,wave)[0]-coords0[0]) >= h[0]*(0.5+TOL) or
                   fabs(Rect_coords(icoords1,wave)[1]-coords0[1]) >= h[1]*(0.5+TOL))
                {
                    printf("ERROR in init_coarse_to_fine(), "
                                  "rect_in_which() failed\n");
                    printf("crds<%g, %g> icrds<%d,%d> found from rect_in_which\n",
                           coords0[0],coords0[1], icoords0[0], icoords0[1]);
                    printf("greped crds of ic<%d, %d> <%g, %g> on tg_grid\n",icoords1[0], icoords1[1],
                          Rect_coords(icoords1,wave)[0],Rect_coords(icoords1,wave)[1]);
                    print_rectangular_grid(rect_g);
                    clean_up(ERROR);
                }
#if defined(_USE_TOL_)
#undef TOL
#undef _USE_TOL_
#endif
                break;
            }
#endif /* defined(TWOD) */
#if defined(THREED)
            case 3:
            {
                break;
            }
#endif /* defined(THREED) */
        }  /* end of switch */

        point_num = 0;
        ic[0] = icoords1[0], ic[1] = icoords1[1];
        comp1 = Rect_comp(icoords1,wave);

        Dens(st) = 0.0;
        Energy(st) = 0.0;
        for (i = 0; i < dim; i++)
             Mom(st)[i] = 0.0;
        if ( comp1 == comp0 )
        {
            point_num++;
            st1 = Rect_state(ic,wave);
            Dens(st) = Dens(st1);
            Energy(st) = Energy(st1);
            for (i = 0; i < dim; i++)
                Mom(st)[i] = Mom(st1)[i];
            reset_gamma(st);
            Set_params(st,st1);
            set_type_of_state(st,GAS_STATE);
        }
        else
        {
            ic[0] = icoords1[0]-1, ic[1] = icoords1[1];
            comp1 = Rect_comp(ic,wave);
            if(comp1 == comp0)
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
                    reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }
            ic[0] = icoords1[0]-1, ic[1] = icoords1[1]-1;
            comp1 = Rect_comp(ic,wave);
            if ( comp1 == comp0 )
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
	            reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }
            ic[0] = icoords1[0], ic[1] = icoords1[1]-1;
            comp1 = Rect_comp(ic,wave);
            if ( comp1 == comp0 )
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
                    reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }

            ic[0] = icoords1[0]-1, ic[1] = icoords1[1]+1;
            comp1 = Rect_comp(ic,wave);
            if ( comp1 == comp0 )
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
                    reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }
            ic[0] = icoords1[0], ic[1] = icoords1[1]+1;
            comp1 = Rect_comp(ic,wave);
            if ( comp1 == comp0 )
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
                    reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }
            ic[0] = icoords1[0]+1, ic[1] = icoords1[1]+1;
            comp1 = Rect_comp(ic,wave);
            if ( comp1 == comp0 )
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
                    reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }
            ic[0] = icoords1[0]+1, ic[1] = icoords1[1];
            comp1 = Rect_comp(ic,wave);
            if ( comp1 == comp0 )
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
                    reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }
            ic[0] = icoords1[0]+1, ic[1] = icoords1[1]-1;
            comp1 = Rect_comp(ic,wave);
            if ( comp1 == comp0 )
            {
                tmp_dist = pt_distance2d(Rect_coords(ic,wave),coords0);
                if(tmp_dist < dist)
                {
                    point_num++;
                    st1 = Rect_state(ic,wave);
                    Dens(st) = Dens(st1);
                    Energy(st) = Energy(st1);
                    for (i = 0; i < dim; i++)
                        Mom(st)[i] = Mom(st1)[i];
                    reset_gamma(st);
                    Set_params(st,st1);
                    set_type_of_state(st,GAS_STATE);
                    dist = tmp_dist;
                }
            }
            if(point_num == 0)
            {
                screen("FT_ERROR in init_coarse_to_fine(), bad interpolation\n");
                screen("interpolate to,  cds<%f, %f>, ic<%d, %d> comp0 = %d on rect_grid\n",
                     coords0[0],coords0[1],icoords0[0],icoords0[1],comp0);
                print_rectangular_grid(rect_g);
                ic[0] = icoords1[0], ic[1] = icoords1[1];
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0]-1, ic[1] = icoords1[1];
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0]-1, ic[1] = icoords1[1]-1;
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0], ic[1] = icoords1[1]-1;
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0]-1, ic[1] = icoords1[1]+1;
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0], ic[1] = icoords1[1]+1;
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0]+1, ic[1] = icoords1[1]+1;
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0]+1, ic[1] = icoords1[1];
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                ic[0] = icoords1[0]+1, ic[1] = icoords1[1]-1;
                comp1 = Rect_comp(ic,wave);
                screen("patch_num=%d, cds<%f, %f>, ic<%d, %d> comp1[%d]\n",
                     wave->patch_number, Rect_coords(ic,wave)[0],
                     Rect_coords(ic,wave)[1],ic[0],ic[1],comp1);
                clean_up(ERROR);
            }
        }
        return st;
}

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
} 

/*  g_overture_injection_FT_buf_after_repatch()   */  
/*                                                */  
/* injection old waves values to (new) wvs, only take care of the interior */
/* scatter_states and fill non-physical amr buffer will take care of buffer zone */
/* It is a fine to fine injection (now extend to the interp. of all */
/* the patches at the same level, if not found, look into the next coarse grid */

LOCAL boolean g_overture_injection_FT_buf_after_repatch(
        Wave            **oldwvs,
        Front           **oldfrs,
        Wave            **wvs,
        Front           **frs)
{
        int       i, k, swp, side;
        int       patch_num = oldwvs[0]->totalNumberOfPatches;
        int       *iperm;
        int       dim = wvs[0]->rect_grid->dim;
        int       step = frs[0]->step;
        int       L[MAXD], U[MAXD];   
        int       where, on_patch, where_level; 
        int       ix, iy; 
        int       ic[MAXD]; 
        Locstate  st;  
        double     crds[MAXD]; 
        COMPONENT comp; 
        Wave            *wavefrom, *waveto;
        Front           *frontfrom, *frontto;
        RECT_GRID       *gridfrom, *gridto;

        DEBUG_ENTER(g_overture_injection_FT_buf_after_repatch)

        iperm = set_iperm(step,dim);
 
        for(i = 1; i < wvs[0]->totalNumberOfPatches; i++)
        {
            for(swp = 0; swp < dim; swp++)
            {
                for(side = 0; side < 2; side++)
                {
                    set_receive_domain(L,U,iperm,side,swp,frs[i]->rect_grid); 
                    for (iy = L[1]; iy < U[1]; iy++)
                    {
                        ic[1] = iy;
                        for (ix = L[0]; ix < U[0]; ix++)
                        {
                            ic[0] = ix;
                            st = Rect_state(ic,wvs[i]);
                            crds[0] = Rect_coords(ic,wvs[i])[0];
                            crds[1] = Rect_coords(ic,wvs[i])[1];
                            comp = Rect_comp(ic,wvs[i]);
                            where_level = best_match_pt_in_other_sets(ic,crds,comp,
                               wvs[i]->patch_level,&on_patch,&where,oldwvs,oldfrs,YES);  
                            if(wvs[i]->patch_level == where_level)
                            {
                               /* pt same level, fine_to_fine */
                                wavefrom = oldwvs[on_patch]; frontfrom = oldfrs[on_patch];
                                gridfrom = wavefrom->rect_grid;
                                ft_assign(st, fine_to_fine(ic,crds,comp,
                                    wavefrom,gridfrom,frontfrom),frontfrom->sizest);
                            }
                            else
                            {
                                /* pt in the coarse grid */
                                wavefrom = oldwvs[on_patch]; frontfrom = oldfrs[on_patch];
                                gridfrom = wavefrom->rect_grid;
                                /* 
                                ft_assign(st, init_coarse_to_fine(ic,crds,comp,
                                  wavefrom,gridfrom,frontfrom), frontfrom->sizest);
                                */ 
                                hyp_solution(crds,comp,NULL,UNKNOWN_SIDE,
                                      frontfrom,wavefrom,st,NULL);
                            }
                        }
                    }
                }  
            }
        } 

        DEBUG_LEAVE(g_overture_injection_FT_buf_after_repatch)
        return YES;  
}

/* injection old waves values to (new) wvs, only take care of the interior */
/* scatter_states and fill non-physical amr buffer will take care of buffer zone */
/* It is a fine to fine injection (now extend to the interp. of all */
/* the patches at the same level, if not found, look into the next coarse grid */

EXPORT boolean g_overture_injection_after_repatch(
        Wave            **oldwvs,
        Front           **oldfrs,
        Wave            **wvs,
        Front           **frs)
{
        Wave            *wavefrom, *waveto;
        Front           *frontfrom, *frontto;
        RECT_GRID       *gridfrom, *gridto;
        Locstate        st;
        double           coords0[MAXD], coords1[MAXD];
        int             icoords0[MAXD], icoords1[MAXD];
        int             dim = wvs[0]->rect_grid->dim;
        int             num_patches, old_num_patches;
        int             patch_num0;
        int             where, on_patch, where_level; 
        int             levels, patch_level0, patch_level1;
        int             i, int1;
        int             int2;
        COMPONENT       comp0,comp1;
        int             smin[MAXD], smax[MAXD];
        int             lbuf[MAXD], ubuf[MAXD];

        size_t          sizest = wvs[0]->sizest;
        int             num_float;
        int             level_diff;

        DEBUG_ENTER(g_overture_injection_after_repatch)

        num_patches = wvs[0]->totalNumberOfPatches;
        old_num_patches = oldwvs[0]->totalNumberOfPatches;
        levels = wvs[0]->NumberOfLevels;

        copy_base_grid_soln_to_new(wvs[0],oldwvs[0]);

        for (int1 = 1; int1 < num_patches; int1++)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;
            patch_level0 = waveto->patch_level;
            patch_num0 = waveto->patch_number;

            if(debugging("injection_repatch_st"))
            {
                printf("g_overture_injection_after_repatch(),"
                        "working on patch[%d],level[%d]\n",
                        waveto->patch_number, waveto->patch_level);
                print_rectangular_grid(waveto->rect_grid);  
            }

            for (i = 0; i < dim; i++)
            {
                smin[i] = 0;
                smax[i] = gridto->gmax[i];
                lbuf[i] = gridto->lbuf[i];
                ubuf[i] = gridto->ubuf[i];
                smin[i] -= lbuf[i];
                smax[i] += ubuf[i];  
            }

            switch(dim)
            {
#if defined(TWOD)
            case 2:
            {
                int      ix, iy;

                for (iy = smin[1]; iy < smax[1]; iy++)
                {
                    icoords0[1] = iy;
                    for (ix = smin[0]; ix < smax[0]; ix++)
                    {
                        icoords0[0] = ix;
                        st = Rect_state(icoords0,waveto);
                        coords0[0] = Rect_coords(icoords0,waveto)[0];
                        coords0[1] = Rect_coords(icoords0,waveto)[1];
                        comp0 = Rect_comp(icoords0,waveto);

                        if(debugging("injection_repatch_st"))
                        {
                            if(icoords0[0] == 13 && icoords0[1] == 17)
                            {
                                printf("ic<%d,%d>, crds<%g,%g>, comp0[%d] on patch[%d] level[%d]\n",
                                   icoords0[0],icoords0[1], coords0[0], coords0[1],
                                   comp0, waveto->patch_number, waveto->patch_level);
                            }  
                        }
                        where_level = best_match_pt_in_other_sets(icoords0,coords0,comp0,
                               waveto->patch_level,&on_patch,&where,oldwvs,oldfrs,YES); 
                       
                        if(debugging("injection_repatch_st"))
                        {
                            if(icoords0[0] == 13 && icoords0[1] == 17)
                            {
                                printf("FOUND on the old wave:\n");
                                print_rectangular_grid(oldwvs[on_patch]->rect_grid);
                            }  
                        }
                        if(waveto->patch_level == where_level)
                        {
                            /* point is found on the old patches at the same level */
                            gridfrom = oldwvs[on_patch]->rect_grid;
                            wavefrom = oldwvs[on_patch];
                            frontfrom = oldfrs[on_patch];
                            ft_assign(st, fine_to_fine(icoords0,coords0,comp0,
                                 wavefrom,gridfrom,frontfrom),frontfrom->sizest);
                        }
                        else
                        {
                            /* pt in the coarse grid */
                            wavefrom = oldwvs[on_patch]; frontfrom = oldfrs[on_patch];
                            gridfrom = wavefrom->rect_grid;
                            if(!is_excluded_comp(comp0,frontto->interf))
                            {
                                hyp_solution(coords0,comp0,NULL,UNKNOWN_SIDE,
                                      frontfrom,wavefrom,st,NULL);
                            }
                            else
                            {
                                g_obstacle_state(st,frontto->sizest);
                            } 
                        }
                    }
                }
                break;
            }
#endif /* defined(TWOD) */
            } /* end of switch */
        }  /* end of for loop */
        /* 
        g_overture_injection_FT_buf_after_repatch(oldwvs,oldfrs,wvs,frs);
        */  

        DEBUG_LEAVE(g_overture_injection_after_repatch)
        return FUNCTION_SUCCEEDED;  
}

/* after repatch, pt looks into the interior of
 * same set of patches (repatched set)
 * at the next coarse level. called by the injection  */

LOCAL int pt_in_patch_coarse_level_repatch(
        int             *icoords0,
        double           *coords0,
        int             comp0,
        int             patch_num0,
        int             *on_patch,
        Wave            **wvs,
        Front           **frs)
{
        Wave            *waveto;
        Front           *frontto;
        RECT_GRID       *gridto;
        double           coords1[MAXD];
        int             icoords1[MAXD];
        int             num_patches, patch_num1;
        int             patch_level0, patch_level1;
        int             i, int1;
        COMPONENT       comp1;
        CompositeGrid   *cg = (CompositeGrid*)wvs[0]->cg_over;

        num_patches = wvs[0]->totalNumberOfPatches;
        patch_level0 = wvs[patch_num0]->patch_level;

        for (int1 = num_patches-1; int1 >= 0 ; int1--)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;

            if(waveto->patch_level >= patch_level0) continue;
            if(patch_level0 - waveto->patch_level > 1) continue;
            if(patch_num0 == waveto->patch_number) continue;


            if ( coords0[0] >= gridto->L[0] and coords0[0] <= gridto->U[0] and
                 coords0[1] >= gridto->L[1] and coords0[1] <= gridto->U[1])
            {
                *on_patch = waveto->patch_number;
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level(), "
                           "rect_in_which() failed\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(icoords1[0] < 0 or icoords1[1] < 0)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level(), "
                           "pt found in the buffer\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1], icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(patch_level0 - waveto->patch_level != 1)

                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level\n");
                    printf("the levels are not right\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d>\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    printf("pt level[%d], patch level[%d]\n", patch_level0,
                           waveto->patch_level);
                    clean_up(ERROR);
                }
                return YES;
            }
        }
        return NO;
}


/* after repatch, pt looks into the interior of
 * other set of patches (old patch set)
 * at the same level. called by teh injection  */
LOCAL int pt_in_patch_same_level_repatch(
        int             *icoords0,
        double           *coords0,
        int             comp0,
        int             patch_level0,
        int             *on_patch,
        Wave            **wvs,
        Front           **frs)
{
        Wave            *waveto;
        Front           *frontto;
        RECT_GRID       *gridto;
        double           coords1[MAXD];
        int             icoords1[MAXD];
        int             num_patches, patch_num1;
        int             patch_level1;
        int             i, int1;
        COMPONENT       comp1;
        CompositeGrid   *cg = (CompositeGrid*)wvs[0]->cg_over;

        num_patches = wvs[0]->totalNumberOfPatches;

        for (int1 = num_patches-1; int1 >= 0 ; int1--)
        {
            waveto = wvs[int1];
            frontto = frs[int1];
            gridto = waveto->rect_grid;

            if(waveto->patch_level != patch_level0) continue;

            if ( coords0[0] >= gridto->L[0] and coords0[0] <= gridto->U[0] and
                 coords0[1] >= gridto->L[1] and coords0[1] <= gridto->U[1])
            {
                *on_patch = waveto->patch_number;
                if(rect_in_which(coords0,icoords1,gridto) == FUNCTION_FAILED)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level(), "
                           "rect_in_which() failed\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(icoords1[0] < 0 or icoords1[1] < 0)
                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level(), "
                           "pt found in the buffer\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d> in rect_grid:\n",
                          coords0[0], coords0[1], icoords0[0],icoords0[1]);
                    print_rectangular_grid(gridto);
                    clean_up(ERROR);
                }
                if(waveto->patch_level != patch_level0)

                {
                    printf("ERROR in pt_in_patch_same_or_coarse_level\n");
                    printf("the levels are not right\n");
                    printf("Looking for coords<%g, %g>ic<%d,%d>\n",
                          coords0[0], coords0[1],icoords0[0],icoords0[1]);
                    printf("pt level[%d], patch level[%d]\n", patch_level0,
                           waveto->patch_level);
                    clean_up(ERROR);
                }
                return YES;
            }
        }
        return NO;
}

LOCAL void copy_base_grid_soln_to_new(
        Wave           *wv,
        Wave           *oldwv)
{
        RECT_GRID       *grid = wv->rect_grid;
        int             i;
        int             dim = wv->rect_grid->dim;
        int             smin[MAXD], smax[MAXD];
        int             lbuf[MAXD], ubuf[MAXD];
        size_t          sizest = wv->sizest;
        int             ix, iy;
        int             ic[2];

        for (i = 0; i < dim; i++)
        {
            smin[i] = 0;
            smax[i] = grid->gmax[i];
            lbuf[i] = grid->lbuf[i];
            ubuf[i] = grid->ubuf[i];
            smin[i] -= lbuf[i];
            smax[i] += ubuf[i];
        }

        for (iy = smin[1]; iy < smax[1]; iy++)
        {
            ic[1] = iy;
            for (ix = smin[0]; ix < smax[0]; ix++)
            {
                ic[0] = ix;
                ft_assign(Rect_state(ic,wv),Rect_state(ic,oldwv),sizest);
            }
        }
}

#endif /* defined(USE_OVERTURE) */

