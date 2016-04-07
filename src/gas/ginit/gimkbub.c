/*
*
*                               gimkbub.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       This file contains functions for the construction of non-trivial
*       geometric curves needed for the initialization of an interface.
*
*/

#if defined(TWOD) || defined(THREED)

#include <ginit/ginit.h>

typedef struct {
        boolean    bubble;
        boolean    eos;
        double   R; /*radius of the bubble region */
        double   r; /*bubble radius*/
        double   frac;/*void fraction*/
        double   L; /*length of the bubble region*/
        double   l; /*length between centers of bubbles*/
        int     size;  /* bubble diameter in mesh size */
        int     N_layers;
        int     N_bubbles;
        double   rand_coef;
        Locstate st; /* initial state for bubbles */
        double   rho; /*density inside bubbles*/
        double   pr; /*pressure inside bubbles*/
        int     eos_n; /*number eos to use*/
        COMPONENT  bubblecomp;
        double   pr_crit; /*critical pressure for making bubbles*/
        double   GAMMA; /*bubbles (gamma - 1)*/
}  BUBBLES;

typedef struct {
        int    icoords[MAXD];  /* icrds of checked cells */
        double  center[MAXD];   /* bubble center */
        int    bmin[MAXD];     /* lower corner indx of blks the bubble ocupy */
        int    bmax[MAXD];     /* upper corner indx of blks the bubble ocupy */
        int    insert;
        double  rad_distr;      /* Mass ditribution radius */
        int    num_cells;      /* number of cells changed from liquid to vapor */
#if defined(USE_OVERTURE)
        int    on_pc_ic[2];
        int    on_patch;
#endif /* if defined(USE_OVERTURE) */
        double  press;       /* pressure in the cell */
        double  vapor_press;    /* vapor pressure from averaged states */
                               /* average states of locations.
                                * Density, momentum and total energy.
                                */
} INSERT_SITE;

         /* LOCAL variables for bubble insertion */
LOCAL   int     loop_bubble = 0;
LOCAL   BUBBLES bubbles;
LOCAL   double   pi = 3.14159265358979323846;
LOCAL   double   bubble_spacing = 1.0;
LOCAL   int     bubble_diameter = 4;

LOCAL   int     make_bubbles(Wave*,Front*,COMPONENT,COMPONENT);
LOCAL   int     make_bubbles_2d(Wave*,Front*,COMPONENT,COMPONENT);
LOCAL   int     make_bubbles_3d(Wave*,Front*,COMPONENT,COMPONENT);

LOCAL 	boolean    check_crossings_in_region(int*,int,Wave*,Front*);
LOCAL   boolean    check_crossings_in_region_ver2(int*,int,Wave*,Front*);
LOCAL	boolean    check_near_by(int*,COMPONENT,int,int,Wave*,Front*,Locstate);
LOCAL   int     press_ascend(const void*,const void*);
LOCAL   int     check_in_bubble_diameter_2d(INSERT_SITE*,Wave*,Front*,
		                     COMPONENT,int*,int*,int*,int*);
LOCAL   int     check_in_bubble_diameter_3d(INSERT_SITE*,Wave*,Front*,
                     COMPONENT,int*,int*,int*,int*);
LOCAL   void    blks_bubble_ocupy(INSERT_SITE,RECT_GRID*,int*,int*,double);
LOCAL   int     check_crx_in_bubble_space_2d(INSERT_SITE,Wave*,Front*,int*,int*,double);
LOCAL   int     check_crx_in_bubble_space_3d(INSERT_SITE,Wave*,Front*,int*,int*,double);
LOCAL   int     intfc_on_rect_domain_2d(INTERFACE*,double*,double*);
LOCAL   int     intfc_on_rect_domain_3d(INTERFACE*,double*,double*);
LOCAL   int     intfc_cross_rect(INTERFACE*,double*,double*);
LOCAL   int     intfc_cross_box(INTERFACE*,double*,double*);
LOCAL   int     tri_intersect_surf(TRI*,int,double,double,double,double,double);
LOCAL   void    avg_states_in_reg_2d(INSERT_SITE*,int*,int*,Wave*,Front*,Locstate,int*);
LOCAL   void    avg_states_in_reg_3d(INSERT_SITE*,int*,int*,Wave*,Front*,Locstate,int*);
LOCAL   int     intfc_intersect_segment(INTERFACE*,int,double,double,double);
LOCAL   int     intfc_intersect_surf(INTERFACE*,int,double,double,double,double,double);
LOCAL   void    set_state_on_phase_boundary(Locstate,int);
LOCAL   void    set_state_on_phase_boundary_ver2(Locstate,Locstate,double);
LOCAL   void    set_state_on_phase_boundary_ver3(Locstate,Locstate,double);
LOCAL   void    create_bubbles(int*,COMPONENT,int,int,Wave*,Front*,Locstate);
LOCAL   void    create_bubbles_ver2(INSERT_SITE*,COMPONENT,Wave*,Front*,Locstate,int*,int*);
LOCAL   void    create_bubbles_3d(INSERT_SITE*,COMPONENT,Wave*,Front*,Locstate,int*,int*);
LOCAL   void    front_init_2d(POINT*,HYPER_SURF*,Locstate,Locstate,Front*,Wave*,Locstate);
LOCAL   void    front_init_3d(POINT*,HYPER_SURF*,Locstate,Locstate,Front*,Wave*,Locstate);
LOCAL   void    reset_front_and_states_ver2(Wave*,Front*,int);

#if defined(USE_OVERTURE)
LOCAL   int     amr_make_bubbles(Wave**,Front**,int,Wv_on_pc**,int,
                                 Overparam*,COMPONENT,COMPONENT);
LOCAL   int     amr_make_bubbles_new_comm(Wave**,Front**,int,Wv_on_pc**,int,
                                 Overparam*,COMPONENT,COMPONENT);
LOCAL  void     amr_reset_front_and_states(Wave**,Front**,int,
                   Wv_on_pc**,int,Overparam*,int*);
LOCAL  void     amr_reset_front_and_states_new_comm(Wave**,Front**,int,
                   Wv_on_pc**,int,Overparam*,int*);
LOCAL  int      conserv_make_bubbles(Wave**,Front**,int,Wv_on_pc**,
                   int,Overparam*,COMPONENT,COMPONENT);
LOCAL   int     count_bubble_locations(Wave*,Front*,COMPONENT,int*,int*,INSERT_SITE**);
LOCAL   int     locate_on_redistr_table(int*,int,Wv_on_pc**,int,int*);
LOCAL   int     make_bubble_conserv_mass(INSERT_SITE*,int,Wave**,Front**,int,Wv_on_pc**,
                     int,double,COMPONENT,int***,int,Locstate,int,INSERT_SITE*);
LOCAL   void    distr_mass_to_neighbors(INSERT_SITE,COMPONENT,int,int*,int*,
                     Wave*,Front*,int**);
LOCAL   int     conserv_make_bubbles_orig_ver2(Wave*,Front*,COMPONENT,COMPONENT);
LOCAL   void    set_insert_sites(INSERT_SITE*,INSERT_SITE*,int,int,int*);
#endif /* if defined(USE_OVERTURE) */

EXPORT  void    g_prompt_for_dynamic_bubble_insertion(
        CHART  *chart)
{
        char      s[121];

#if defined(USE_OVERTURE)
            chart->amr_make_bubbles = NULL;
#else /* if defined (USE_OVERTURE) */
            chart->phase_solver = NULL;
#endif /* if defined (USE_OVERTURE) */

        screen("Type 'y' to dynamically insert vapor bubbles in the"
               " liquid (y, n(dflt)): ");

        (void) Gets(s);
        if (strncmp(s,"y",1) == 0)
        {
#if defined(USE_OVERTURE)
            chart->amr_make_bubbles = amr_make_bubbles;
#else /* if defined (USE_OVERTURE) */
            chart->phase_solver = make_bubbles;
#endif /* if defined (USE_OVERTURE) */
            screen("Enter bubble diameter(dflt 4): ");
            (void) Scanf("%d\n",&bubble_diameter);
            screen("Enter bubble spacing(dflt 1.0): ");
            (void) Scanf("%f\n",&bubble_spacing);
        }
}

#if defined(USE_OVERTURE)
LOCAL  int    amr_make_bubbles(
        Wave                 **wvs,
        Front                **frs,
        int                  num_patches,
        Wv_on_pc             **redistr_table,
        int                  max_n_patch,
        Overparam            *overparam,
        COMPONENT            l_comp,
        COMPONENT            bubblecomp)
{
        int   num_bubble = 0;
        int   i;
        Wave  *wave;
        Front *front;
        int   cur_bubble;

	/*
        num_bubble = conserv_make_bubbles(wvs,frs,num_patches,
           redistr_table,max_n_patch,overparam,l_comp,bubblecomp);
	*/   

        for ( i = 0; i < num_patches; i++)
        {
            wave = wvs[i];
            front = frs[i];
            if (wave->patch_level == wave->NumberOfLevels - 1)
            {
                cur_bubble = make_bubbles_orig_ver2(wave,front,l_comp,bubblecomp);
                /* cur_bubble = make_bubbles_orig(wave,front,l_comp,bubblecomp); */
                if ( cur_bubble > 0 )
                {
                     front->interf->modified = YES;
                     num_bubble += cur_bubble;
                     if (debugging("bubbles_mpi"))
                         printf("cur_bubble = %d, num_bubble = %d\n",
                            cur_bubble, num_bubble);
                }
            }
        }

        amr_reset_front_and_states(wvs,frs,num_patches,redistr_table,
                                   max_n_patch,overparam,&num_bubble);

        return YES; 
}

LOCAL  int    amr_make_bubbles_new_comm(
        Wave                 **wvs,
        Front                **frs,
        int                  num_patches,
        Wv_on_pc             **redistr_table,
        int                  max_n_patch,
        Overparam            *overparam,
        COMPONENT            l_comp,
        COMPONENT            bubblecomp)
{
        int   num_bubble = 0;
        int   i;
        Wave  *wave;
        Front *front;
        int   cur_bubble;

	/*
        start_clock("make_bubble");
        num_bubble = conserv_make_bubbles(wvs,frs,num_patches,
           redistr_table,max_n_patch,overparam,l_comp,bubblecomp);
        stop_clock("make_bubble");
	*/

        start_clock("make_bubble");
        for ( i = 0; i < num_patches; i++)
        {
            wave = wvs[i];
            front = frs[i];
            if ( wave->patch_level == wave->NumberOfLevels - 1 )
            {
                cur_bubble = make_bubbles_orig_ver2(wave,front,l_comp,bubblecomp);
                /* cur_bubble = make_bubbles_orig(wave,front,l_comp,bubblecomp); */
                if ( cur_bubble > 0 )
                {
                     front->interf->modified = YES;
                     num_bubble += cur_bubble;
                     if (debugging("bubbles_mpi"))
                         printf("cur_bubble = %d, num_bubble = %d\n",
                            cur_bubble, num_bubble);
                }
            }
        }
        stop_clock("make_bubble");

        start_clock("bubble_comm");
	/* amr_reset_front_and_states */
        amr_reset_front_and_states_new_comm(wvs,frs,num_patches,redistr_table,
                                   max_n_patch,overparam,&num_bubble);
        stop_clock("bubble_comm");
        if(0 == num_bubble)
            return NO;

        return YES; 
}
#endif /* #if defined(USE_OVERTURE) */

LOCAL int make_bubbles(
        Wave                 *wave,
        Front                *front,
        COMPONENT            l_comp,
        COMPONENT            bubblecomp)
{
        int      num_bubbles;

        switch(wave->rect_grid->dim)
        {
	case 2:
	    num_bubbles =  make_bubbles_2d(wave,front,l_comp,bubblecomp);
	    break;
        case 3:
            /*#bjet2 */
	    num_bubbles = make_bubbles3d(wave,front);
            /*num_bubbles = make_bubbles_3d(wave,front,l_comp,bubblecomp); */
	    /* liuxt12 */
	    interface_reconstructed(front->interf) = YES;
            break;
        default:
            printf("ERROR: make_bubbles(), implement for dim %d\n",
                   wave->rect_grid->dim);
            clean_up(ERROR);
        }
        reset_front_and_states_ver2(wave,front,num_bubbles);

        return YES;
}

/* Feb 18 2004: Yarema: insert bubble */
LOCAL   int     make_bubbles_orig(
        Wave                    *wave,
        Front                   *front,
        COMPONENT               l_comp,
        COMPONENT               bubblecomp)
{
        RECT_GRID               *gr = front->rect_grid;
        int                     dim = gr->dim;
        int                     icoords[MAXD],imin[3],imax[3];
        int                     ix,iy, i;
        double                   *GL = front->rect_grid->GL;
        double                   *GU = front->rect_grid->GU;
        double                   *L = front->rect_grid->L;
        double                   *U = front->rect_grid->U;
        double                   coords[MAXD];
        INTERFACE               *intfc = front->interf;
        int                     size, spacing;
        int                     num_bubble = 0;
        double                   period = 5.e-6;
        Locstate                state;
        double                   T_sat, P_sat;
        static Locstate         ref_st = NULL;
        static Locstate         rst = NULL;

        debug_print("bubbles_mpi","Entering make_bubbles()\n");

        bubbles.bubblecomp = bubblecomp;

        if (bubbles.bubble != YES && loop_bubble == 0)
        {
            if (ref_st == NULL)
            {
                alloc_state(intfc,&ref_st,front->sizest);
                alloc_state(intfc,&rst,front->sizest);
                Params(ref_st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                Params(rst) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_ambient_comp_type(comp_type(bubbles.bubblecomp),front);
                set_type_of_state(ref_st,GAS_STATE);
                set_type_of_state(rst,GAS_STATE);
                zero_state_velocity(ref_st,Params(ref_st)->dim);

                /*Origin: from S2PH
                Dens(ref_st) = 0.0005;
                Energy(ref_st) = 0.004963755209666677/(1.04117 - 1.0);
                */
                /* values are taken at T = 294.26 K = 70 F
                 * P = 0.049368 bar
                 * V = 4829.82024036586 cm^3/g
                 * from n-heptane values
                 *
                Dens(ref_st) = 1./4829.82024036586;
                Energy(ref_st) =  0.049368/(1.05 - 1.0);
                */
                /* values are taken at T = 305.37 K = 90 F
                 * P = 0.085272 bar
                 * V = 2897.88714998276 cm^3/g
                 * from n-heptane values
                Dens(ref_st) = 1./2897.88714998276;
                Energy(ref_st) =  0.085272/(1.0516 - 1.0);
                 */
                /* values are taken at T = 288 K
                 * P = 1 bar
                 * rho = 0.00419 g/cm^3
                 * from n-heptane values in gas state
                Dens(ref_st) = 0.00419;
                Energy(ref_st) =  1./(1.05 - 1.0);
                 */
                /* values are taken at T =  298.05329710794149 K
                 * P = 1.013 bar
                 * rho = 0.0040948445783132534 g/cm^3
                 * from n-heptane values in gas state
                 */
                Dens(ref_st) = 0.0040948445783132534;
                Energy(ref_st) =  1.013/(1.05 - 1.0);
            }
            if ( bubbles.st == NULL )
            {
                alloc_state(intfc,&bubbles.st,front->sizest);
                Params(bubbles.st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_type_of_state(bubbles.st,GAS_STATE);
                zero_state_velocity(bubbles.st,Params(bubbles.st)->dim);
            }

            size = 4;
	    if (debugging("big_bubble"))
		size = 8;

            spacing = 1;
            if (debugging("big_spacing"))
                spacing = 5;

            if (debugging("bubbles_mpi"))
            {
                printf("No bubbles, so try to insert it with size = %d and spacing = %d\n",
                       size, spacing);
            }

            if (debugging("big_bubble")) size *= 2;
            bubbles.bubble = YES;
            bubbles.pr_crit = -10.0;
            bubbles.N_bubbles = 1;
            bubbles.N_layers = 1;
            bubbles.L = 2 * gr->h[0];
            bubbles.l = 3 * gr->h[0];
            bubbles.r = size/2 * min ( gr->h[0], gr->h[1])*0.95;
            /*bubbles.r = size/2 * min ( gr->h[0], gr->h[1]); */
            bubbles.eos_n = 2;
            bubbles.eos = YES;

            /*Origin: from S2PH
            bubbles.pr = ( debugging("one_atm") ) ? 1.013: 0.004963755209666677;
            */
            /*from n-heptane
            bubbles.pr =  0.085272;
            */
            /*from n-heptane in gas state
            bubbles.pr =  1.;
            */
            /*from n-heptane in gas state*/
            bubbles.pr =  1.013;

            state_on_adiabat_with_pr(ref_st, bubbles.pr, rst, GAS_STATE);
            bubbles.rho = Dens(rst);

            /*Origin: from S2PH
            bubbles.GAMMA = 1.04117 - 1.0;
            bubbles.R = 0.0317172856847711;
            bubbles.R = 0.156722; wrong number from old setting 
            */
            /*from n-heptane
            bubbles.GAMMA = 1.0516 - 1.0;
            bubbles.R = 0.7803;
            */
            /*from n-heptane in gas state*/
            bubbles.GAMMA = 1.05 - 1.0;
            bubbles.R = 0.83;

            if (debugging("bubbles_mpi"))
            {
                printf("Bubble:\n");
                printf("\tr = %g, rho = %g, pr = %g, GAMMA = %g\n",
                       bubbles.r, bubbles.rho, bubbles.pr, bubbles.GAMMA);
            }
        }

        if  (bubbles.bubble == YES)
        {
            /** mnkim's code **/
            /**
            imin[0] = size*spacing - gr->lbuf[0];
            imin[1] = size*spacing - gr->lbuf[1];
            imax[0] = wave->rect_grid->gmax[0] - (spacing+1)* size + 1;
            imax[1] = wave->rect_grid->gmax[1] - (spacing+1)* size + 1;
            if(rect_boundary_type(intfc,0,0) == REFLECTION_BOUNDARY)
                imin[0] += gr->lbuf[0];
            if(rect_boundary_type(intfc,1,0) == REFLECTION_BOUNDARY)
                imin[1] += gr->lbuf[1];
            **/
            
            /** My code **/
            /* dir ,side */
            /***
            imin[0] = size*spacing;
            imin[1] = size*spacing;
            imax[0] = wave->rect_grid->gmax[0] - (spacing+1)* size + 1;
            imax[1] = wave->rect_grid->gmax[1] - (spacing+1)* size + 1;
            if(rect_boundary_type(intfc,0,0) == REFLECTION_BOUNDARY)
                imin[0] -= gr->lbuf[0];
            if(rect_boundary_type(intfc,1,0) == REFLECTION_BOUNDARY)
                imin[1] -= gr->lbuf[1];
            ***/
 
            for ( ix = 0; ix < dim; ix ++ )
            { 
                imin[ix] = 0; 
                if ( size*spacing > gr->lbuf[ix] )
                    imin[ix] = size*spacing - gr->lbuf[ix]; 
                imax[ix] = wave->rect_grid->gmax[ix] - size + 1; 
                if ( size*spacing > gr->ubuf[ix] )
                    imax[ix] -= size*spacing - gr->ubuf[ix]; 
            }

            if ( debugging("bubble_in_nozzle") )
            {
                if(  gr->L[1] >= 0.05 || gr->U[1] <= 0.05 )
                {
                    imax[0] = 0;
                    imax[1] = 0;
                }
                bubbles.pr_crit = 0;
                if ( !debugging("insert_once") )
                    add_to_debug("insert_once");
            }
            if (debugging("imax"))
            {
                printf("\nmake_bubbles: main loop starts\n");
                printf("ix = (%d, %d), iy = (%d, %d), gmax = (%d, %d)\n",
                       imin[0], imax[0], imin[1], imax[1],
                       wave->rect_grid->gmax[0], wave->rect_grid->gmax[1]);
                printf("L[0],L[1],U[0],U[1],GL[0],GL[1],GU[0],GU[1] = %f %f %f %f %f %f %f %f\n",
                       gr->L[0],gr->L[0],gr->U[0],gr->U[1],gr->GL[0],gr->GL[1],gr->GU[0],gr->GU[1]);
                printf("Buffer : lbuf = (%d, %d), ubuf = (%d, %d)\n", gr->lbuf[0], gr->lbuf[1],
                       gr->ubuf[0],gr->ubuf[1]);
            }

            for (iy = imin[1]; iy < imax[1]; iy++)
            {
                icoords[1] = iy;
                if ( debugging("insert_once") && num_bubble > 0 ) break;
                for (ix = imin[0]; ix < imax[0]; ix++)
                {
                    icoords[0] = ix;
                    for(i = 0; i < dim; i++)
                        coords[i] = Rect_coords(icoords,wave)[i];

                    if (Rect_comp(icoords,wave) != l_comp) continue;

                    state = Rect_state(icoords,wave);
                    T_sat = temperature(state);

                    /*TMP*/
                    if ( T_sat <= 278.983333333333341 )
                        P_sat = pressure(state);
                    else
                        P_sat = saturated_pres(state,T_sat);

                    if ( debugging("PT_curve") )
                    {
                        printf("PT curve: p = %"FFMT", T = %"FFMT"\n",
                               P_sat, T_sat);
                        clean_up(ERROR);
                    }
                    if ( P_sat > pressure(state) )
                    {
                        if (debugging("bubbles_mpi"))
                        {
                            printf("coords = (%g,%g), pressure = %"FFMT", P_sat = %"FFMT", "
                                   "T_sat = %"FFMT"\n",
                                   coords[0],coords[1],pressure(state), P_sat, T_sat);
                        }

                        /* call create_bubble */
                        copy_state(ref_st,state);
                        if (debugging("bubbles_mpi"))
			    printf("Before check_near_by gmax[%d %d]\n",front->rect_grid->gmax[0],
					    front->rect_grid->gmax[1]); 	
                        if (check_near_by(icoords,l_comp,size,spacing,wave,front,ref_st))
                        {
                            /* Jul 7 2004: Zhiliang: check crossings */
                            /* another function for more strict crossing check */
                            if (debugging("bubbles_mpi"))
			        printf("Pass check_near_by\n"); 	
                            if(NO == check_crossings_in_region_ver2(icoords,size,wave,front))
                            {
                                /***
                                printf(" Bubbles are inserted at [%d %d], (%g %g) size[%d]"
                                 " on wave[%d] level[%d] step[%d]\n",
                                    icoords[0], icoords[1],  Rect_coords(icoords,wave)[0], 
                                     Rect_coords(icoords,wave)[1], size, 
                                     wave->patch_number, wave->patch_level, front->step);
                                printf("ix = (%d, %d), iy = (%d, %d) gmax = (%d, %d)\n",
                                    imin[0], imax[0], imin[1], imax[1],
                                    wave->rect_grid->gmax[0], wave->rect_grid->gmax[1]);
                                print_rectangular_grid(wave->rect_grid); 
                                ***/

                                set_state_on_phase_boundary(ref_st,spacing);
                                create_bubbles(icoords,l_comp,size,spacing,wave,front,ref_st);
                                num_bubble ++;
                                if ( debugging("insert_once") ) break;
                            }
                        }
                        else if (debugging("bubbles_mpi"))
                            printf("Less than %d x %d region\n", size, size);
                    }
                }
            }
        }

        bubbles.bubble = NO;
        debug_print("bubbles_mpi","Left make_bubbles()\n");
        return num_bubble;
}       /* end make_bubbles_orig */

LOCAL void reset_front_and_states_ver2(
        Wave                      *wave,
        Front                     *front,
        int                       num_bubble)
{
        static int                num_bubbles = 0;
        int                       *iperm;
        INTERFACE                 *intfc = front->interf;
        int                       ix;
        /* Local statistics */
        if ( num_bubble > 0 )
        {
            intfc->modified = YES;
            printf("%d bubbles are added at time step %d\n",
                   num_bubble, front->step);
        }

        pp_global_isum(&num_bubble,1);

        /* Communication */
        if ( num_bubble > 0 )
        {
            if (!scatter_front(front))
            {
                printf("ERROR in make_bubbles, scatter_front() failed: "
                       "number of bubbles = %d\n", num_bubble);
                clean_up(ERROR);
            }

            reinit_hyp_solution_function(wave,front);

            iperm = set_iperm(front->step,front->rect_grid->dim);
            for ( ix = 0; ix < front->rect_grid->dim; ix++)
            {
                if (!scatter_states(wave,front,iperm,ix))
                {
                    screen("ERROR in make_bubbles, scatter_states() "
                           "failed in %d direction, number of bubbles = "
                           "%d\n",ix,num_bubble);
                    clean_up(ERROR);
                }
            }

            /* Global statistics */
            num_bubbles += num_bubble;
            if (debugging("insert_once")) loop_bubble ++;
            if (debugging("more_bubbles") && num_bubbles <= 1) loop_bubble --;
            if ( is_io_node(pp_mynode()) )
                printf("Total %d, new %d at time step %d\n",
                       num_bubbles, num_bubble, front->step);
        }
        pp_gsync();
}

LOCAL void front_init_2d(
        POINT      *p,
        HYPER_SURF *hs,
        Locstate   lstate,
        Locstate   rstate,
        Front      *front,
        Wave       *wave,
        Locstate   rst)
{

        COMPONENT         pcomp = NO_COMP, ncomp = NO_COMP;
        INTERFACE         *intfc = hs->interface;
        int               dim = intfc->dim;

        int               stype;
        int               ic[MAXD];
        int               i, j;
        double             coordp[MAXD];
        Locstate          state;
        double             dens;
        double             v[MAXD];
        double             pl, tl;

        if (debugging("bubbles_mpi"))
        {
          printf("\nfront_init1\n");
        }

        ncomp = negative_component(hs);
        pcomp = positive_component(hs);

        if (debugging("bubbles_mpi"))
        {
            printf("nct->type = %d, EXTERIOR = %d, \n",comp_type(ncomp)->type, EXTERIOR);
            printf("In default, ncomp = %d, bubbles.bubblecomp = %d\n"
                   ,ncomp, bubbles.bubblecomp);
        }

        if (ncomp != bubbles.bubblecomp)
        {
            if ( debugging("from_nbhd") )
            {
                rect_in_which(Coords(p),ic,wave->rect_grid);
                for(i=0; i < dim; i++)
                    coordp[i] = Rect_coords(ic,wave)[i];
                if ( Rect_comp(ic,wave) != 2 )
                {
                    printf("ERROR in front_init(), all points on the curve of "
                           "the bubble must be liquid, not comp = %d\n", Rect_comp(ic,wave));
                    clean_up(ERROR);
                }
                state = Rect_state(ic,wave);
                stype = state_type(state);
                copy_state(lstate,state);
                /*
                set_state(lstate,FGAS_STATE,lstate);
                /*verbose_print_state("lstate1",lstate); */
                tl = Temperature(lstate);
                Params(lstate) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                /*verbose_print_state("lstate2",lstate); */
                Dens(lstate) = bubbles.rho;
                /*verbose_print_state("lstate3",lstate); */
                pl = pressure(lstate);
                Params(lstate) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[0];
                /*verbose_print_state("lstate4",lstate); */
                Dens(lstate) = (pl+8935.62992886953)/18.291/1.792700566721809/tl;
                /*Press(lstate) = pl; */
                /*verbose_print_state("lstate5",lstate); */
                /*clean_up(ERROR); */
                /* */
                  set_state(lstate,TGAS_STATE,lstate);
                  for ( i = 0; i < dim; i++ )
                  v[i] = Vel(state)[i];
                  state_on_adiabat_with_pr(state,bubbles.pr,lstate,TGAS_STATE);
                  for ( i = 0; i < dim; i++ )
                  Vel(lstate)[i] = v[i];
                /* */
                set_state(lstate,stype,lstate);
                */
            }
            else
                copy_state(lstate,rst);

            if ( debugging("bubbles_mpi") )
            {
                verbose_print_state("front_state(), lstate",lstate);
            }
        }

        if (debugging("bubbles_mpi"))
        {
            printf("pct->type = %d, BUBBLE = %d, \n",comp_type(pcomp)->type, BUBBLE);
            printf("In default, pcomp = %d, bubbles.bubblecomp = %d\n"
                   ,pcomp, bubbles.bubblecomp);
        }

        if (pcomp == bubbles.bubblecomp)
        {
            /*
            copy_state(rstate,lstate);
            Params(rstate) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
            set_state(rstate,FGAS_STATE,rstate);
            Temperature(rstate) = temperature(lstate);
            Dens(rstate) = bubbles.rho;
            /*Press(rstate) = bubbles.pr; */
            set_state(rstate,stype,rstate);
            */

            copy_state(rstate,bubbles.st);

            if (debugging("bubbles_mpi"))
            {
                verbose_print_state("front_state(), rstate",rstate);
            }
        }


        if (debugging("bubbles_mpi"))
        {
            printf("Bubbles:\n");
            print_general_vector("negative side state at ",Coords(p),dim,"");
            (void) printf(" on hs %llu\n",hypersurface_number(hs));
            print_gas_state(lstate);
            print_general_vector("positive side state at ",Coords(p),dim,"");
            (void) printf(" on hs %llu\n",hypersurface_number(hs));
            print_gas_state(rstate);
            printf("\nLeft front_init\n");
        }
        debug_print("bubbles_mpi","Left front_init()\n");
}       /* end of front_init_2d */

/* Create Bubbles */
LOCAL   void create_bubbles(
        int                     *icoords,
        COMPONENT               l_comp,
        int                     size,
        int                     spacing,
        Wave                    *wave,
        Front                   *front,
        Locstate                rst)
{
        double                   r,theta,d_theta,step;
        int                     n_theta;
        RECT_GRID               *gr = front->rect_grid;
        int                     dim = gr->dim;
        int                     idirs[3],imin[3],imax[3];
        int                     i, j, k;
        int                     ix,iy,iz;

        double                   *GL = front->rect_grid->GL;
        double                   *GU = front->rect_grid->GU;
        double                   *L = front->rect_grid->L;
        double                   *U = front->rect_grid->U;

        double                   CG[MAXD];
        double                   CL[MAXD];

        double                   c[MAXD];
        double                   coords[MAXD];

        int                     ic[MAXD];
        Locstate                state;

        NODE                    *ns, *ne;
        INTERFACE               *intfc = front->interf;
        CURVE                   *cv;
        BOND                    *b;
        INTERFACE               *sv_intfc = current_interface();


        set_current_interface(intfc);

        for(i = 0; i < dim; i++)
            CG[i] = Rect_coords(icoords,wave)[i];

        ic[0] = icoords[0] + size - 1;
        ic[1] = icoords[1] + size - 1;

        for(i = 0; i < dim; i++)
            CL[i] = Rect_coords(ic,wave)[i];

        debug_print("bubbles_mpi","Entering create_bubbles()\n");

        c[0] = (CG[0] + CL[0])/2.0;
        c[1] = (CG[1] + CL[1])/2.0;
        r = bubbles.r ;

        if (debugging("bubbles_mpi"))
        {
            printf("center = (%g, %g), ",c[0], c[1]);
            printf("CG = (%g, %g), ",CG[0], CG[1]);
            printf("r = %g, 0.95*min(dx,dy) = %g\n",
                   r, 0.95*min(gr->h[0],gr->h[1]));
        }

        /* make node and curve */
        coords[0] = c[0] + r;
        coords[1] = c[1];
        ns = make_node(Point(coords));
        ne = ns;
        node_type(ns) = CLOSED_NODE;
        cv = make_curve(l_comp,bubbles.bubblecomp,ns,ne);

        n_theta = 4*(irint( pi*r/1.5/min(gr->h[0], gr->h[1])) + 1);
        d_theta = 2*pi/n_theta;

        if (debugging("bubbles_mpi"))
            printf("after make_curve: n_theta = %d, n_theta/4 = %d, d_theta = %g\n",
                   n_theta,irint(n_theta/4), d_theta);

        /* insert points into the curve */
        theta = -d_theta;
        for (k = 1; k < n_theta; k++)
        {
            coords[0] = c[0] + r*cos(theta);
            coords[1] = c[1] + r*sin(theta);
            Check_return(insert_point_in_bond(Point(coords),cv->last,cv),
                         create_bubbles);
            if (debugging("bubbles_mpi"))
            {
                printf("start = (%g, %g), point = (%g, %g),\n",
                       Coords(cv->last->start)[0],
                       Coords(cv->last->start)[1],
                       coords[0], coords[1]);
            }
            theta -= d_theta;
        }

        wave_type(cv) = CONTACT;
        start_status(cv) = end_status(cv) = INCIDENT;
	create_time(cv) = front->time;

        /*init front_states*/

        if (debugging("bubbles_mpi"))
        {
            (void) printf("Initializing states on curve %llu\n",
                          curve_number(cv));
            print_curve(cv);
        }
        front_init_2d(cv->first->start,
                   Hyper_surf(cv),
                   left_start_state(cv),right_start_state(cv),
                   front,wave,rst);
        front_init_2d(cv->first->start,
                   Hyper_surf(cv),
                   left_state(cv->first->start),right_state(cv->first->start),
                   front,wave,rst);
        for (b = cv->first; b != cv->last; b = b->next)
            front_init_2d(b->end,
                       Hyper_surf(cv),
                       left_state(b->end),right_state(b->end),
                       front,wave,rst);
        front_init_2d(cv->last->end,
                   Hyper_surf(cv),left_state(cv->last->end),
                   right_state(cv->last->end),front,wave,rst);
        front_init_2d(cv->last->end,
                   Hyper_surf(cv),left_end_state(cv),
                   right_end_state(cv),front,wave,rst);
        if (debugging("bubbles_mpi"))
        {
            (void) printf("After initializing states on curve %llu\n",
                          curve_number(cv));
            print_curve(cv);
        }

        /*init interior states for bubble*/
        k = 0;
        for ( iy = icoords[1]; iy < icoords[1] + size; iy++)
        {
            ic[1] = iy;
            for ( ix = icoords[0]; ix < icoords[0] + size; ix++)
            {
                double   dist = 0.0;
                int     stype;

                ic[0] = ix;
                coords[0] = Rect_coords(ic,wave)[0];
                coords[1] = Rect_coords(ic,wave)[1];

                for ( i = 0; i < dim; i++)
                    dist += (coords[i] - c[i])*(coords[i] - c[i]);
                dist = sqrt(dist);
                if ( dist > bubbles.r ) continue;

                state = Rect_state(ic,wave);
                if (debugging("bubbles_mpi"))
                {
                    (void) printf("Original, comp = %d, ",Rect_comp(ic,wave));
                    printf("coords = (%g, %g), ",coords[0],coords[1]);
                    printf("distance = %g <= %g\n", dist, bubbles.r);
                    verbose_print_state("state",state);
                }
                stype = state_type(state);
                Rect_comp(ic,wave) = bubbles.bubblecomp;
                copy_state(state,bubbles.st);
                set_state(state,stype,state);

                k ++;

                if (debugging("bubbles_mpi"))
                {
                    (void) printf("New state, comp = %d\n", Rect_comp(ic,wave));
                    verbose_print_state("state",state);
                }
            }
        }
        if (debugging("bubbles_mpi"))
        {
            printf("End of initilization of bubbles in %d cells\n", k);
        }

        if ( debugging("skip_nbhd") )
        {
            if (debugging("bubbles_mpi"))
            {
                printf("\nInterface with bubbles after create_bubbles():\n");
                /*print_interface(intfc); */
                /*show_intfc_states(front->interf); */
            }
            debug_print("bubbles_mpi","Left create_bubbles()\n");

            set_current_interface(sv_intfc);
            return;
        }

        k = 0;
        /*init interior states for liquid surrounding bubble*/
        for ( iy = icoords[1] - spacing*size; iy < icoords[1] + (spacing+1)*size; iy++)
        {
            ic[1] = iy;
            for ( ix = icoords[0] - spacing*size ; ix < icoords[0] + (spacing+1)*size; ix++)
            {
                double   dist = 0.0;
                int     stype;

                ic[0] = ix;
                coords[0] = Rect_coords(ic,wave)[0];
                coords[1] = Rect_coords(ic,wave)[1];

                for ( i = 0; i < dim; i++)
                    dist += (coords[i] - c[i])*(coords[i] - c[i]);
                dist = sqrt(dist);
                if ( dist > bubbles.r*(2*spacing+1) || dist <= bubbles.r ) continue;

                state = Rect_state(ic,wave);
                if (debugging("bubbles_mpi"))
                {
                    (void) printf("Original, comp = %d, ",Rect_comp(ic,wave));
                    printf("coords = (%g, %g), ",coords[0],coords[1]);
                    printf("distance = %g <= %g\n", dist, bubbles.r*(spacing+1) );
                    verbose_print_state("state",state);
                }

                stype = state_type(state);
                copy_state(state,rst);
                set_state(state,stype,state);

                k++;

                if (debugging("bubbles_mpi"))
                {
                    (void) printf("New state, comp = %d\n", Rect_comp(ic,wave));
                    verbose_print_state("state",state);
                }
            }
        }

        if (debugging("bubbles_mpi"))
        {
            printf("End of initilization of liquid in %d cells\n", k);
            printf("\nInterface with bubbles after create_bubbles():\n");
            /*print_interface(intfc); */
            /*show_intfc_states(front->interf); */
        }

        set_current_interface(sv_intfc);
        debug_print("bubbles_mpi","Left create_bubbles()\n");
}

LOCAL   void     set_state_on_phase_boundary(
        Locstate lst,
        int      spacing)
{
        int   stype;
        double ratio = 2*spacing+1.;
        double alpha = ratio * ratio * ratio;
        double rho_i, e_i;
        double rho_v, e_v, e_l, beta;
        double v[2], v2;
        int   i;

        for ( i = 0; i < 2; i++)
        {
            v[i] = vel(i,lst);
            v2  += sqr(v[i]);
        }
        stype = state_type(bubbles.st);
        set_state(bubbles.st,FGAS_STATE,bubbles.st);
        Temperature(bubbles.st) = temperature(lst);
        /*TMP*/
        if ( Temperature(bubbles.st) <= 278.983333333333341 )
        {
            printf("ERROR in set_state_on_phase_boundary\n");
            verbose_print_state("lst", lst);
            clean_up(ERROR);
        }

        Dens(bubbles.st) = saturated_pres(lst,Temperature(bubbles.st))
                          /bubbles.R/Temperature(bubbles.st);

        for( i=0; i<2; i++) Vel(bubbles.st)[i] = v[i];

        set_state(bubbles.st,stype,bubbles.st);

        if ( debugging("bubbles_mpi") )
        {
            verbose_print_state("bubbles.st in regular form",bubbles.st);
            printf("v2 = %"FFMT"\n",v2);
        }

        if ( debugging("no_change_lst") )
        {
            if ( debugging("bubbles_mpi") )
            {
                verbose_print_state("no change in lst",lst);
            }
            return;
        }
        if ( debugging("saturated_liquid") )
        {
            stype = state_type(lst);
            set_state(lst,TGAS_STATE,lst);
            Dens(lst) = density_via_PT(lst, pressure(bubbles.st),
                                       temperature(bubbles.st));
            Press(lst) = pressure(bubbles.st);
            set_state(lst,stype,lst);
            if ( debugging("bubbles_mpi") )
            {
                verbose_print_state("lst in saturation point",lst);
            }
            return;
        }

        stype = state_type(lst);
        set_state(lst,GAS_STATE,lst);
        rho_i = Dens(lst);
        e_i = energy(lst);
        rho_v = Dens(bubbles.st);
        e_v = energy(bubbles.st);

        Dens(lst) = ( alpha * rho_i - rho_v )/(alpha - 1);
        e_l = ( alpha * e_i - e_v - heat_of_vaporization_via_T(lst,temperature(lst)) * rho_v )/(alpha - 1);

        Energy(lst) = e_l;
        for ( i = 0; i< 2; i++ ) Mom(lst)[i] = v[i] * Dens(lst);

        if ( debugging("bubbles_mpi") )
        {
            printf("internal energy = %"FFMT"\n", internal_energy(lst));
            printf("pressure = %"FFMT"\n", pressure(lst));
            printf("temperature = %"FFMT"\n", temperature(lst));
        }

        set_state(lst,stype,lst);
        if ( debugging("bubbles_mpi") )
        {
            printf("rho_i = %"FFMT", e_i = %"FFMT", rho_v = %"FFMT", e_v = %"FFMT", e_l = %"FFMT"\n",
                   rho_i, e_i, rho_v, e_v, e_l);
            printf("beta = %"FFMT"\n", beta);
            verbose_print_state("lst in regular form",lst);
            printf("Saturated pressure = %"FFMT" at liquid temperature = %"FFMT"\n",
                   saturated_pres(lst,temperature(lst)), temperature(lst));
            printf("Saturated pressure = %"FFMT" at vapor temperature = %"FFMT"\n",
                   saturated_pres(lst,temperature(bubbles.st)), temperature(bubbles.st));
            /*clean_up(ERROR); */
        }
}

LOCAL boolean check_near_by(
        int                     *icoords,
        COMPONENT               l_comp,
        int                     size,
        int                     spacing,
        Wave                    *wave,
        Front                   *front,
        Locstate                rst)
{
        int      ix,iy, i, j;
        int      ic[MAXD];

        int      stype;
        Locstate state;
        double    P_sat, T_sat;
        int      num_cells;
        double    *h = front->rect_grid->h;
        double    CG[MAXD], CL[MAXD], c[MAXD], coords[MAXD];

        if ( icoords[0] == wave->rect_grid->gmax[0] ||
             icoords[1] == wave->rect_grid->gmax[1] ) return NO;

        CG[0] = Rect_coords(icoords,wave)[0];
        CG[1] = Rect_coords(icoords,wave)[1];
        ic[0] = icoords[0] + size - 1;
        ic[1] = icoords[1] + size - 1;

        CL[0] = Rect_coords(ic,wave)[0];
        CL[1] = Rect_coords(ic,wave)[1];

        for ( i = 0; i < 2; i++) c[i] = (CG[i] + CL[i])/2.0;

        /* Only look at (size)x(size) region */
        for ( iy = icoords[1];  iy < icoords[1] + size; iy++)
        {
            ic[1] = iy;
            for ( ix = icoords[0]; ix < icoords[0] + size; ix++)
            {
                double   dist[5];
                int     pivot[2][5] ={ {0, 1, -1, -1, 1}, {0, 1, 1, -1, -1}};
                boolean    pass = YES;

                ic[0] = ix;
                coords[0] = Rect_coords(ic,wave)[0];
                coords[1] = Rect_coords(ic,wave)[1];

                for ( j = 0; j < 5; j++)
                {
                    dist[j] = 0.0;
                    for ( i = 0; i < 2; i++)
                      dist[j] += sqr(coords[i] + pivot[i][j]*h[i]/2. - c[i]);
                    dist[j] = sqrt(dist[j]);

                    if ( dist[j] <= bubbles.r )
                    {
                        pass = NO;
                        break;
                    }
                }

                if ( pass ) continue;
                if ( Rect_comp(ic,wave) == l_comp )
                {
                    state = Rect_state(ic,wave);
                    T_sat = temperature(state);

                    /*TMP*/
                    if ( T_sat <= 278.983333333333341 )
                    {
                        /*
                        printf("Warning low temperature = %"FFMT" in "
                               "check_near_by at (%"FFMT", %"FFMT")\n",
                               T_sat,coords[0], coords[1]);
                        verbose_print_state("liquid state", state);
                        */
                        return NO;
                    }

                    P_sat = saturated_pres(state,T_sat);
                    if ( P_sat <= pressure(state) )
                        return NO;
                }
                else
                    return NO;
            }
        }

        if( debugging("bubbles_mpi") )
	    printf("check_near_by, pass check 4 by 4\n");	
        num_cells = 0;
        stype = state_type(rst);
        set_type_of_state(rst,EGAS_STATE);

        Dens(rst) = 0.0;
        Energy(rst) = 0.0;
        for ( i=0; i<2; i++) Mom(rst)[i] = 0.0;

        if( debugging("bubbles_mpi") )
	    printf("check_near_by, check near-by cells to check crossing of front\n");	
        /* check near-by cells to check crossing of front */
        for ( iy = icoords[1] - spacing*size; iy < icoords[1] + (spacing+1)*size; iy++)
        {
            ic[1] = iy;
            for ( ix = icoords[0] - spacing*size; ix < icoords[0] + (spacing+1)*size; ix++)
            {
                double   dist[5];
                int     pivot[2][5] ={ {0, 1, -1, -1, 1}, {0, 1, 1, -1, -1}};
                boolean    pass = YES;

                ic[0] = ix;
                coords[0] = Rect_coords(ic,wave)[0];
                coords[1] = Rect_coords(ic,wave)[1];

                for ( j = 0; j < 5; j++)
                {
                    dist[j] = 0.0;
                    for ( i = 0; i < 2; i++)
                      dist[j] += sqr(coords[i] + pivot[i][j]*h[i]/2. - c[i]);
                    dist[j] = sqrt(dist[j]);

                    if ( dist[j] <= bubbles.r*(2*spacing+1) )
                    {
                        pass = NO;
                        break;
                    }
                }
                if ( pass ) continue;

                if ( Rect_comp(ic,wave) == l_comp )
                {
                    state = Rect_state(ic,wave);

                    if ( debugging("high_pres") )
                    {
                        T_sat = temperature(state);
                        P_sat = saturated_pres(state,T_sat);
                        if ( P_sat*10 <= pressure(state) )
                            return NO;
                    }
                    num_cells ++;
                    Dens(rst) += Dens(state);
                    Energy(rst) += specific_internal_energy(state);

                    if (debugging("same_vel"))
                        for ( i=0; i<2; i++) Mom(rst)[i] += mom(i,state);

                    if( debugging("bubbles_mpi") )
		    {
		        printf("wave %p, ic_comp %d, num_cells %d\n", wave, Rect_comp(ic,wave), num_cells);	
                        verbose_print_state("state",Rect_state(ic,wave));
		    }
                }
                else
                    return NO;
            }
        }

        for ( i=0; i<2; i++) Mom(rst)[i] /= num_cells;
        Dens(rst) /= num_cells;
        Energy(rst) /= num_cells;

        if( debugging("bubbles_mpi") )
        {
            printf("Number of cells = %d\nAverage density = %"FFMT"\nAverage energy = %"FFMT"\n"
                   "Average temperature = %"FFMT"\nAverage pressure = %"FFMT"\n"
                   "Average velocity = (%"FFMT", %"FFMT")\n",
                   num_cells, Dens(rst), Energy(rst), temperature(rst),
                   pressure(rst), Vel(rst)[0], Vel(rst)[1]);
            printf("Saturated pressure = %"FFMT"\n",
                   saturated_pres(rst,temperature(rst)));
        }

        set_state(rst,stype,rst);

        if( debugging("bubbles_mpi") )
        {
            verbose_print_state("rst",rst);
            /*clean_up(ERROR); */
        }

        return YES;
}

/* Jul 7 2004: Zhiliang: check crossings */
/* If the crossings are found in the region, return YES */
LOCAL   boolean     check_crossings_in_region(
        int         *icoords,
        int         size,
        Wave        *wave,
        Front       *front)
{
        int         i, j, k, dim = wave->rect_grid->dim;
        RECT_GRID   *gr = wave->rect_grid;
        double       *CG, *CL;
        double       c[MAXD];
        double       coords[MAXD];
        int         ic[MAXD], icrds[MAXD];
        int             nc[4];   /* number of CRXING */
        static int      first = YES;
        static CRXING   **crx[4];
        TRI_GRID        *ntg = wave_tri_soln(wave)->tri_grid;

        if(YES == first)
        {
            first = NO;
            for(i = 0; i < 4; i++)
                uni_array(&crx[i],4,sizeof(CRXING *));
        }

        CG = Rect_coords(icoords,wave);
        for(i = 0; i < dim; i++)
            ic[i] = icoords[i] + size - 1;
        CL = Rect_coords(ic,wave);
        for(i = 0; i < dim; i++)
            c[i] = (CG[i] + CL[i])/2.0;

        if(debugging("check_crossings_in_region"))
        {
            printf("Create bubble, icrds[%d %d], ic[%d %d]\n",
                icoords[0], icoords[1], icoords[0] + size - 1,
               icoords[1] + size - 1);
            printf("center = (%g, %g), \n",c[0], c[1]);
            printf("CL = (%g, %g), ",CL[0], CL[1]);
            printf("CG = (%g, %g) \n",CG[0], CG[1]);
            printf("r = %g, 0.95*min(dx,dy) = %g\n",
                 bubbles.r, 0.95*min(gr->h[0],gr->h[1]));
            printf("BLOCK region L(%g %g) R(%g %g)\n",
                cell_edge(icoords[0], 0, gr),  
                cell_edge(icoords[1], 1, gr),  
                cell_edge(ic[0], 0, gr) + gr->h[0], 
                cell_edge(ic[1], 1, gr) + gr->h[1]);  
        }

        for(i = icoords[0]; i < ic[0]; i++)
        {
            for(j = icoords[1]; j < ic[1]; j++)
            {
                icrds[0] = i; icrds[1] = j;
                nc[0] = crossings_in_direction(crx[0],icrds,EAST,ntg);
                nc[2] = crossings_in_direction(crx[2],icrds,NORTH,ntg);
                if(nc[0] > 0 || nc[2] > 0)
                    return YES;
            }
        }

        /* Left side */
        for(j = icoords[1]; j <= ic[1]; j++)
        {
            icrds[0] = icoords[0]; icrds[1] = j;
            nc[1] = crossings_in_direction(crx[1],icrds,WEST,ntg);
            if(nc[1] > 0)
            {
                c[0] = CG[0] - 0.5*gr->h[0];
                for(k = 0; k < nc[1]; k++)
                {
                    if(Coords(crx[1][k]->pt)[0] > c[0])
                        return YES;
                }
            }
        }

        /* Right side */
        for(j = icoords[1]; j <= ic[1]; j++)
        {
            icrds[0] = ic[0]; icrds[1] = j;
            nc[0] = crossings_in_direction(crx[0],icrds,EAST,ntg);
            if(nc[0] > 0)
            {
                c[0] = CL[0] + 0.5*gr->h[0];
                for(k = 0; k < nc[0]; k++)
                {
                    if(Coords(crx[0][k]->pt)[0] < c[0])
                        return YES;
                }
            }
        }

        /* Lower side */
        for(i = icoords[0]; i <= ic[0]; i++)
        {
            icrds[0] = i, icrds[1] = icoords[1];
            nc[3] = crossings_in_direction(crx[3],icrds,SOUTH,ntg);
            if(nc[3] > 0)
            {
                c[1] = CG[1] - 0.5*gr->h[1];
                for(k = 0; k < nc[3]; k++)
                {
                    if(Coords(crx[3][k]->pt)[1] > c[1])
                        return YES;
                }
            }
        }
        /* Upper side */
        for(i = icoords[0]; i <= ic[0]; i++)
        {
            icrds[0] = i, icrds[1] = ic[1];
            nc[2] = crossings_in_direction(crx[2],icrds,NORTH,ntg);

            if(nc[2] > 0)
            {
                c[1] = CL[1] + 0.5*gr->h[1];
                for(k = 0; k < nc[2]; k++)
                {
                    if(Coords(crx[2][k]->pt)[1] < c[1])
                        return YES;
                }
            }
        }

        /* Upper side */
        for(i = icoords[0]; i < ic[0]; i++)
        {
            icrds[0] = i, icrds[1] = ic[1];
            nc[2] = crossings_in_direction(crx[2],icrds,EAST,ntg);
            if(nc[2] > 0)
                return YES;
        }

        /* Right side */
        for(j = icoords[1]; j < ic[1]; j++)
        {
            icrds[0] = ic[0]; icrds[1] = j;
            nc[0] = crossings_in_direction(crx[0],icrds,NORTH,ntg);
            if(nc[0] > 0)
                return YES;
        }

        return NO;
}

/* If the crossings are found in the region, return YES */
LOCAL   boolean     check_crossings_in_region_ver2(
        int         *icoords,
        int         size,
        Wave        *wave,
        Front       *front)
{
        int         i, j, k, dim = wave->rect_grid->dim;
        RECT_GRID   *gr = wave->rect_grid;
        double       *CG, *CL;
        double       center[MAXD];
        double       coords[MAXD];
        double       cl[MAXD], cu[MAXD]; 
        int         ic[MAXD], icrds[MAXD];
        INTERFACE   *intfc = front->interf;

        CG = Rect_coords(icoords,wave);
        for(i = 0; i < dim; i++)
            ic[i] = icoords[i] + size - 1;
        CL = Rect_coords(ic,wave);
        for(i = 0; i < dim; i++)
            center[i] = (CG[i] + CL[i])/2.0;

        if(debugging("check_crossings_in_region"))
        {
            printf("Create bubble, icrds[%d %d], ic[%d %d]\n",
                icoords[0], icoords[1], icoords[0] + size - 1,
               icoords[1] + size - 1);
            printf("center = (%g, %g), \n",center[0], center[1]);
            printf("CL = (%g, %g), ",CL[0], CL[1]);
            printf("CG = (%g, %g) \n",CG[0], CG[1]);
            printf("r = %g, 0.95*min(dx,dy) = %g\n",
                 bubbles.r, 0.95*min(gr->h[0],gr->h[1]));
        }

        for(i = 0; i < dim; i++)
        {
            cl[i] = cell_edge(icoords[i], i, gr); 
            cu[i] = cell_edge(ic[i], i, gr) + gr->h[i]; 
        }
        if(debugging("check_crossings_in_region"))
        {
            printf("BLOCK region L(%g %g) U(%g %g)\n",
                 cl[0], cl[1], cu[0], cu[1]); 
        }

        /* Interface cross with left rectangle boundary */
        if(YES == intfc_intersect_segment(intfc, 0, cl[0], cl[1], cu[1]))
            return YES; 

        /* Interface cross with right rectangle boundary */
        if(YES == intfc_intersect_segment(intfc, 0, cu[0], cl[1], cu[1]))
            return YES; 

        /* Interface cross with lower rectangle boundary */
        if(YES == intfc_intersect_segment(intfc, 1, cl[1], cl[0], cu[0]))
            return YES; 

        /* Interface cross with upper rectangle boundary */
        if(YES == intfc_intersect_segment(intfc, 1, cu[1], cl[0], cu[0]))
            return YES; 

        if(YES == intfc_cross_rect(intfc,cl, cu))
        {
            printf("WARNING: intfc_cross_rect detected crossings\n"); 
            return YES; 
        }
        return NO; 
}

LOCAL  int    intfc_cross_rect(
        INTERFACE       *intfc,        /* interface to be cut */
        double           *L,             /* lower coordinate of cut line in other direction */
        double           *U)             /* upper coordinate of cut line in other direction */
{
        CURVE           **cc;
        BOND            *b;
        POINT           *p;

        for (cc = intfc->curves; cc && *cc; ++cc)
        {
            p = (*cc)->first->start;  
            if(Coords(p)[0] >= L[0] && Coords(p)[0] <= U[0] &&
               Coords(p)[1] >= L[1] && Coords(p)[1] <= U[1])
            { 
                /*
                printf("WARNING: point in rect[%g %g] [%g %g]\n",
                      L[0], L[1], U[0], U[1]); 
                print_curve(*cc); 
                */
                return YES; 
            }
            for (b = (*cc)->first; b != NULL; b = b->next)
            {
                p = b->end;  
                if(Coords(p)[0] >= L[0] && Coords(p)[0] <= U[0] &&
                   Coords(p)[1] >= L[1] && Coords(p)[1] <= U[1])
                {
                    /*
                    printf("WARNING: point in rect[%g %g] [%g %g]\n",
                          L[0], L[1], U[0], U[1]); 
                    print_curve(*cc); 
                    */
                    return YES; 
                }
            }
        }
        return NO; 
}	/* end of intfc_cross_rect */

LOCAL  int    intfc_cross_box(
        INTERFACE       *intfc,        /* interface to be cut */
        double           *L,             /* lower coordinate of cut line in other direction */
        double           *U)             /* upper coordinate of cut line in other direction */
{
        POINT             *p;
        SURFACE           **s;
        TRI               *tri;
        int               num_surfs, num_tris, k;

        for (num_tris = 0, s = intfc->surfaces; s && *s; ++s)
        {
            num_tris += (*s)->num_tri;
            for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
            {
                for (k = 0; k < 3; ++k)
                    p = Point_of_tri(tri)[k]; 
                if(Coords(p)[0] >= L[0] && Coords(p)[0] <= U[0] &&
                   Coords(p)[1] >= L[1] && Coords(p)[1] <= U[1] && 
                   Coords(p)[2] >= L[2] && Coords(p)[2] <= U[2])
                {
                    return YES;
                }
            }
        }
        return NO;
}	/* intfc_cross_box */

#if defined(TWOD)
LOCAL  int    intfc_intersect_segment(
        INTERFACE       *intfc,        /* interface to be cut */
        int             dir,           /* direction of cut line normal */
        double           cut,           /* coordinate of cut line */
        double           L,             /* lower coordinate of cut line in other direction */
        double           U)             /* upper coordinate of cut line in other direction */
{
        CURVE           **cc;
        BOND            *b;
        POINT           *p;
        double           cut_crds[MAXD]; 

        for (cc = intfc->curves; cc && *cc; ++cc)
        {
            for (b = (*cc)->first; b != NULL; b = b->next)
            {
                if ((p = bond_crosses_cut_segment(b,
                          dir, 0, cut, L, U, cut_crds)) == NULL)
                    continue;
                else
                    return YES; 
            }
        }
     
        return NO; 
}
#endif /* defined(TWOD) */

#if defined(USE_OVERTURE)
LOCAL void amr_reset_front_and_states(
        Wave                      **wvs,
        Front                     **frs,
        int                       num_patches,
        Wv_on_pc                  **redistr_table,
        int                       max_n_patch,
        Overparam                 *overparam,
        int                       *num_bubble)
{
        static int                num_bubbles = 0;
        int                       *iperm;
        int                       i;
        Front                     *front;
        Wave                      *wave;
        /* PP_cps                    **recv_cps; */
        /* Amr_intrp_set             *getaintrp, *sndaintrp; */
        int                       get_blk_n = 0, snd_blk_n = 0;

        debug_print("gimkcur", "Entered amr_reset_front_and_states\n"); 

        /* Local statistics */
        if (*num_bubble > 0)
        {
            printf("amr_reset_front_and_states(), myid = %d, num_bubble[%d] at step[%d]\n",
                   pp_mynode(), *num_bubble, frs[0]->step);
        }
        pp_global_isum(num_bubble,1);
	if(*num_bubble > 0) 
            printf("after global, num_bubble = %d at step %d\n",
	          *num_bubble,frs[0]->step);

        /* Communication */
        if (*num_bubble > 0)
        {
	    set_adjust_bdry_crx_flag(NO);
            if (!scatter_patch_fronts(frs,num_patches,redistr_table,
                                      max_n_patch,overparam,NO))
            {
                printf("ERROR in make_bubbles, "
                       "scatter_patch_fronts() failed: "
                       "number of bubbles = %d\n", *num_bubble);
                clean_up(ERROR);
            }
	    set_adjust_bdry_crx_flag(YES);

            /***
            if (!assembly_distribute_patch_fronts(frs,num_patches,
                                   redistr_table,max_n_patch, YES))
            {
                screen("ERROR in make_bubbles, "
                       "assembly_distribute_patch_fronts() "
                       "failed, number of bubbles = "
                       "%d\n",*num_bubble);
                clean_up(ERROR);
            }
            ***/

            if (!assembly_coarse_patch_fronts(frs,num_patches,
                                   redistr_table,max_n_patch))
            {
                screen("ERROR in make_bubbles, "
                       "assembly_distribute_patch_fronts() "
                       "failed, number of bubbles = "
                       "%d\n",*num_bubble);
                clean_up(ERROR);
            }

            for (i = 0; i < num_patches; i++)
            {
                wave = wvs[i];
                front = frs[i];
                reinit_hyp_solution_function(wave,front);
            }

            start_clock("amr_bubble_scat");
            average_fine_to_coarse_st(wvs,frs,num_patches,
                    redistr_table,max_n_patch,overparam);
            wave_intrp_buf_blks(wvs,frs,num_patches,redistr_table,
                    max_n_patch,overparam);
            stop_clock("amr_bubble_scat");

            /**** Do not have to use these functions 
            iperm = set_iperm(frs[0]->step,frs[0]->interf->dim);

            for ( i = 0; i < frs[0]->interf->dim; i++)
            {

#if !defined(USE_OLD_CODE_DEBUG)
                if(i == 0)
                {
                    start_clock("amr_intrp1_reset");
                    wave_set_recv_intrp_blks_in_dir(frs,wvs,num_patches,
                     overparam,redistr_table,max_n_patch,iperm,i,&getaintrp,
                     &get_blk_n,&sndaintrp,&snd_blk_n);
                    recv_cps = wave_post_recv_intrp_blks(frs,wvs,
                     num_patches,getaintrp,get_blk_n,overparam,redistr_table,
                     max_n_patch);
                    stop_clock("amr_intrp1_reset");
                    start_clock("amr_restl_reset");
                    wave_intrp_buf_blks_in_dir(frs,wvs,num_patches,
                      overparam,redistr_table,max_n_patch,iperm,i,
                      recv_cps,getaintrp,get_blk_n,sndaintrp,snd_blk_n);
                    stop_clock("amr_restl_reset");
                }
                else
                {
                    start_clock("amr_restl2_reset");
                    average_fine_to_coarse_st(wvs,frs,num_patches,
                      redistr_table,max_n_patch,overparam);
                    wave_set_recv_intrp_blks_in_dir(frs,wvs,num_patches,
                      overparam,redistr_table,max_n_patch,iperm,i,&getaintrp,
                      &get_blk_n,&sndaintrp,&snd_blk_n);
                    recv_cps = wave_post_recv_intrp_blks(frs,wvs,num_patches,
                      getaintrp,get_blk_n,overparam,redistr_table,max_n_patch);
                    wave_intrp_buf_blks_in_dir(frs,wvs,num_patches,
                      overparam,redistr_table,max_n_patch,iperm,i,
                      recv_cps,getaintrp,get_blk_n,sndaintrp,snd_blk_n);
                    stop_clock("amr_restl2_reset");
                }
#endif /* if !defined(USE_OLD_CODE_DEBUG)  */
            }
            ****/

            /* Global statistics */
            num_bubbles += *num_bubble;
            if (debugging("insert_once")) loop_bubble ++;
            if (debugging("more_bubbles") && num_bubbles <= 1) loop_bubble --;
            if ( is_io_node(pp_mynode()) )
                printf("Total %d, new %d at time step %d\n",
                       num_bubbles, *num_bubble, front->step);
        }

        debug_print("gimkcur", "Left amr_reset_front_and_states\n"); 
}

LOCAL void amr_reset_front_and_states_new_comm(
        Wave                      **wvs,
        Front                     **frs,
        int                       num_patches,
        Wv_on_pc                  **redistr_table,
        int                       max_n_patch,
        Overparam                 *overparam,
        int                       *num_bubble)
{
        static int                num_bubbles = 0;
        int                       *iperm;
        int                       i;
        Front                     *front;
        Wave                      *wave;
        int                       get_blk_n = 0, snd_blk_n = 0;
	int                       next_step_regrid = NO;

        debug_print("gimkcur", "Entered amr_reset_front_and_states\n"); 

        /* Local statistics */
        if (*num_bubble > 0)
        {
            printf("amr_reset_front_and_states(), myid = %d, num_bubble[%d] at step[%d]\n",
                   pp_mynode(), *num_bubble, frs[0]->step);
        }
        pp_global_isum(num_bubble,1);
	if(*num_bubble > 0) 
            printf("after global, num_bubble = %d at step %d\n",
	          *num_bubble,frs[0]->step);

        /* Communication */
        if (*num_bubble > 0)
        {
	    if(overparam->period != 0 &&
	       (frs[0]->step+1) % overparam->period == 0)
                next_step_regrid = YES;

            set_scatter_base_fr_flag(NO);
            if (!scatter_patch_fronts(frs,num_patches,redistr_table,
                                      max_n_patch,overparam,NO))
            {
                printf("ERROR in make_bubbles, "
                       "scatter_patch_fronts() failed: "
                       "number of bubbles = %d\n", *num_bubble);
                clean_up(ERROR);
            }

	    if(next_step_regrid == YES)
	    {
                if (!assembly_coarse_patch_fronts(frs,num_patches,
                                   redistr_table,max_n_patch))
                {
                    screen("ERROR in make_bubbles, "
                       "assembly_distribute_patch_fronts() "
                       "failed, number of bubbles = "
                       "%d\n",*num_bubble);
                    clean_up(ERROR);
                }
	    }
            /*
            /* It seems that this function changes fine grid interface too much */
            if( !assembly_distribute_patch_fronts(frs,num_patches,
                      redistr_table,max_n_patch, YES))
            {
                screen("ERROR in make_bubbles, "
                       "assembly_distribute_patch_fronts() "
                       "failed, number of bubbles = "
                       "%d\n",*num_bubble);
                clean_up(ERROR);
            }
            */

            /* Only modify the finest grids */
            for (i = 0; i < num_patches; i++)
            {
                wave = wvs[i];
                front = frs[i];
                if(wave->patch_level == wave->NumberOfLevels-1)
                    reinit_hyp_solution_function(wave,front);
		else if(next_step_regrid == YES)
                    reinit_hyp_solution_function(wave,front);
            }

            /*
            start_clock("amr_bubble_scat");
            average_fine_to_coarse_st(wvs,frs,num_patches,
                    redistr_table,max_n_patch,overparam);
            wave_intrp_buf_blks(wvs,frs,num_patches,redistr_table,
                    max_n_patch,overparam);
            stop_clock("amr_bubble_scat");
            */

            /* Global statistics */
            num_bubbles += *num_bubble;
            if (debugging("insert_once")) loop_bubble ++;
            if (debugging("more_bubbles") && num_bubbles <= 1) loop_bubble --;
            if ( is_io_node(pp_mynode()) )
                printf("Total %d, new %d at time step %d\n",
                       num_bubbles, *num_bubble, front->step);
        }

        debug_print("gimkcur", "Left amr_reset_front_and_states\n"); 
}
#endif /* if defined(USE_OVERTURE) */

/* make_bubbles_2d() will search for the
 * locations of the negative pressure first. Then
 * try to insert bubbles from the location that has
 * the lowest pressure to the location has higher pressure.
 */
LOCAL   int     make_bubbles_2d(
        Wave                    *wave,
        Front                   *front,
        COMPONENT               l_comp,
        COMPONENT               bubblecomp)
{
        RECT_GRID               *gr = front->rect_grid;
        int                     dim = gr->dim;
        int                     icoords[MAXD],imin[MAXD],imax[MAXD];
        int                     bmin[MAXD], bmax[MAXD], Rbmin[MAXD], Rbmax[MAXD];
        int                     ix,iy, i, j;

        INTERFACE               *intfc = front->interf;
        int                     size; 
        double                   spacing;
        int                     num_bubble = 0;
        double                   period = 5.e-6;
        Locstate                state;
        double                   T_sat, P_sat, press, ref_temp;
        static Locstate         ref_st = NULL, rst = NULL, distr_st;
        static int              First = YES;
        INSERT_SITE             *insert_site = NULL;        
        int                     n_site = 0, alloc_site = 0;
        int                     num_cells, **distr_site;
        double                   alpha = 1.0;
	static	double		templ,temph;

        debug_print("bubbles_mpi","Entering make_bubbles_2d()\n"); 
        if(YES == First)
        {
            bubbles.st = NULL;
            bubbles.bubble = NO;
            First = NO;
	    phase_temp_range(&templ,&temph);
        }

        bubbles.bubblecomp = bubblecomp;
        if ( bubbles.bubble != YES && loop_bubble == 0)
        {
            if (ref_st == NULL)
            {
                alloc_state(intfc,&ref_st,front->sizest);
                alloc_state(intfc,&rst,front->sizest);
                alloc_state(intfc,&distr_st,front->sizest);
                Params(ref_st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                Params(rst) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_ambient_comp_type(comp_type(bubbles.bubblecomp),front);
                set_type_of_state(ref_st,GAS_STATE);
                set_type_of_state(rst,GAS_STATE);
                zero_state_velocity(ref_st,Params(ref_st)->dim);

                /* values are taken at T =  298.05329710794149 K
                 * P = 1.013 bar
                 * rho = 0.0040948445783132534 g/cm^3
                 * from n-heptane values in gas state
                 */
                Dens(ref_st) = 0.0040948445783132534;
                Energy(ref_st) =  1.013/(1.05 - 1.0);
            }
            if (bubbles.st == NULL)
            {
                alloc_state(intfc,&bubbles.st,front->sizest);
                Params(bubbles.st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_type_of_state(bubbles.st,GAS_STATE);
                zero_state_velocity(bubbles.st,Params(bubbles.st)->dim);
            }

            bubbles.bubble = YES;
            /* bubble diameter in mesh size */
            /* bubbles.size = size = 4; */
	    bubbles.size = size = bubble_diameter;
            /* mini. distance between two bubble boundary in bubble diameter */
            /* spacing = 1.0; */
	    spacing = bubble_spacing;

            bubbles.r = 0.98*size/2.0*min(gr->h[0], gr->h[1]);
            bubbles.N_bubbles = 1;
            bubbles.N_layers = 1;
            bubbles.L = 2 * gr->h[0];
            bubbles.l = 3 * gr->h[0];
            bubbles.eos_n = 2;
            bubbles.eos = YES;

            /* Problem specific init. */

            /* Critical pressure to allow bubble creation */
            /* Or can choose to use Clausis-Claperon in case of cavitation */
            bubbles.pr_crit = -10.0;

            /* Again, Clausis-Claperon can be used to set bubble state */
            /*from n-heptane in gas state*/
            bubbles.pr =  1.013;

            state_on_adiabat_with_pr(ref_st, bubbles.pr, rst, GAS_STATE);
            bubbles.rho = Dens(rst);

            /* EOS parameter should not appear. must remove */
            /*from n-heptane in gas state*/
            bubbles.GAMMA = 1.05 - 1.0;
            bubbles.R = 0.83;
        }        

        if  (bubbles.bubble == YES)
        {
            /* domain which allows bubble to be created */
            /* The current domain setting might give insufficient
             * neighboring space if a bubble is created near the boundary, when
             * the spacing is greater than 1 bubble diameter 
             * and bubble diameter is greater than buffer size (4).
             * However, to include these in the domain setting might never
             * create bubbles near the boundary. 
             */
            for(i = 0; i < dim; i++)
            {
                imin[i] = 0;
                imax[i] = gr->gmax[i]; 
                /** For conservative insertion, do not allow across reflection line **/
		/*
                if(rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
                    imin[i] -= gr->lbuf[i];
                if(rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
                    imax[i] += gr->ubuf[i];
		*/
            }

            alloc_site = max(imax[0],imax[1]); 
            uni_array(&insert_site, alloc_site, sizeof(INSERT_SITE)); 
            bi_array(&distr_site, gr->gmax[0], gr->gmax[1], sizeof(int));
            for(iy = 0; iy < gr->gmax[1]; iy++)
                for(ix = 0; ix < gr->gmax[0]; ix++)
                    distr_site[ix][iy] = NO;

            /* Search pressure below specified critical value */
            for (iy = imin[1]; iy < imax[1]; iy++)
            {
                icoords[1] = iy;
                for (ix = imin[0]; ix < imax[0]; ix++)
                {
                    icoords[0] = ix;
                    if (Rect_comp(icoords,wave) != l_comp) continue;

                    state = Rect_state(icoords,wave);
                    T_sat = temperature(state);
                    P_sat = saturated_pres(state,T_sat);
                    press = pressure(state);

                    if (press+5 < P_sat) /* for simulation sp_mp_jet40_new3.r */
                    /* if (press+0.18 < P_sat)  for simulation sp_mp_jet40_new3.r */
                    /* from homo model */
                    /* if (press < -5.0) for simulation sp_mp_jet40_new2.r */
                    /* Clausis-Clapeyron */
                    /* if (P_sat > press + 0.18)  */
                    /* from \delta P = 2*surface_tension/radius. */
                    /* Here radius=2*1.1125 microns, s_tension=1.96e-2 kg/s^2 */
                    {
                        ft_assign(insert_site[n_site].icoords,icoords,INT*MAXD); 
                        insert_site[n_site].press = press; 
                        insert_site[n_site].insert = NO;
                        n_site++;
                        if(n_site == alloc_site)
                        {
                            INSERT_SITE *tmp_site; 
                            alloc_site += max(imax[0],imax[1]);
                            uni_array(&tmp_site, alloc_site, sizeof(INSERT_SITE));
                            for(i = 0; i < n_site; i++)
                                ft_assign(&tmp_site[i], &insert_site[i], sizeof(INSERT_SITE));
                            free(insert_site);
                            insert_site = tmp_site;
                        }
                    }
                }
            }
            
            /* put insert_site in order according to pressure */
            if(n_site != 0)
                qsort((POINTER)insert_site, n_site,sizeof(INSERT_SITE),press_ascend);

            if (debugging("bubbles_mpi"))
            {
#if defined(USE_OVERTURE)
                (void) printf("On patch[%d], level[%d]\n", 
                          wave->patch_number, wave->patch_level);
#endif /* if defined(USE_OVERTURE) */
                (void) printf("domain size imin[%d %d], imax[%d %d]\n",
                             imin[0], imin[1], imax[0], imax[1]);
                (void) printf("Find %d possible insertion locations\n", n_site);
                for(i = 0; i < n_site; i++)
                    (void)printf("At coords[%g, %g] icoords[%d, %d], pressure %g\n",
                        Rect_coords(insert_site[i].icoords,wave)[0],
                        Rect_coords(insert_site[i].icoords,wave)[1],
                        insert_site[i].icoords[0],insert_site[i].icoords[1],
                        insert_site[i].press);
            }

            /* Create bubbles */
            for(i = 0; i < n_site; i++)
            {
                if(YES != check_in_bubble_diameter_2d(&(insert_site[i]), 
                    wave,front,l_comp,imin,imax,bmin,bmax))
                    continue;

                /* Geometric specific */
                /* Do not create bubbles close to nozzle wall */
		if(debugging("apart_from_wall"))
		{
                    if((fabs(insert_site[i].center[0]-0.00979) < (bubbles.r + 3.0*gr->h[0])) &&
                        insert_site[i].center[1] < 0.15  && 
                        insert_site[i].center[1] > 0.05)
                        continue;
		}
		if(debugging("awayfrmwll"))
		{
                    if((fabs(insert_site[i].center[0]-0.010235) < (bubbles.r + 4.0*gr->h[0])) &&
                        insert_site[i].center[1] < 0.1502  && 
                        insert_site[i].center[1] > 0.045)
                        continue;
		}
		/* liuxt12 not insert in chamber */
		if (debugging("out_of_chamber") && insert_site[i].center[1] < 0.05)
		    continue;

                blks_bubble_ocupy(insert_site[i],gr,Rbmin,Rbmax,spacing);

                /*** Old alg. bubbles are always created in spacing distance to curves ***/
                if(YES == check_crx_in_bubble_space_2d(insert_site[i], 
                        wave,front,Rbmin,Rbmax,spacing))
                    continue;

                /*** New alg.: Only bubbles created at the same time must 
                 *** have the spacing. Bubbles can be created right next previous 
                 *** time step bubbles.
                 ***/
		/**
                if(YES == new_bubbles_in_domain(insert_site,i-1,Rbmin,Rbmax))
                    continue;               
                if(YES == check_crx_in_bubble_space_2d(insert_site[i], 
                        wave,front,Rbmin,Rbmax,0.0))
                    continue;
		**/    
                /*** End of New alg. ***/

                avg_states_in_reg_2d(&(insert_site[i]),bmin,bmax,wave,front,ref_st,&num_cells);

                /** Compute total quantity needs to be distributed **/
                copy_state(distr_st, ref_st);

		if (is_rotational_symmetry())
		{
                    if (gr->Remap.remap == CYLINDRICAL_REMAP)
                        alpha = insert_site[i].center[0];
		    alpha = 1.0;
		}
                Dens(distr_st) *= alpha*pi*sqr(bubbles.r);
                Mom(distr_st)[0] *= alpha*pi*sqr(bubbles.r);
                Mom(distr_st)[1] *= alpha*pi*sqr(bubbles.r);
                Energy(distr_st) *= alpha*pi*sqr(bubbles.r);

                /* The reference temperature will be used to determine  */
                /* vapor bubble states and init. phase bondary state. */

                ref_temp = temperature(ref_st);
		/* Temperature range test */
		if(ref_temp < templ || ref_temp > temph)
		{
	            printf("ERROR: in make_bubbles_2d,"
		         " ref_st temp %g is outside range[294.261 533.15]\n", ref_temp);
                    printf("Center[%g %g][%d %d] ,%d cells: center to ax %g\n",
                        insert_site[i].center[0], insert_site[i].center[1], 
			insert_site[i].icoords[0],
			insert_site[i].icoords[1], num_cells,
                        insert_site[i].center[0]/gr->h[0]);
		    verbose_print_state("ref_st", ref_st);
		    /* clean_up(ERROR); */
		}
                set_state_on_phase_boundary_ver2(ref_st, rst, ref_temp); 

		/* liuxt12 */
		printf("bubble will be inserted at:\n");
                print_general_vector("center",insert_site[i].center,3,"\n");
                print_int_vector("icoords",insert_site[i].icoords,3,"\n");

                create_bubbles_ver2(&(insert_site[i]),l_comp,wave,front,rst,bmin,bmax);

                insert_site[i].insert = YES;
                for(j = 0; j < dim; j++)
                {
                    insert_site[i].bmin[j] = bmin[j]; 
                    insert_site[i].bmax[j] = bmax[j]; 
                }
                num_bubble++;
            }
            free(insert_site);
            free(distr_site); 
        }

	/* liuxt12 */
	if (num_bubble > 0)
	    printf("%d bubbles inserted\n",num_bubble);
        bubbles.bubble = NO;
        debug_print("bubbles_mpi","Left make_bubbles_2d()\n"); 
        return num_bubble;
}	/* end of make_bubbles_2d */

#if defined(USE_OVERTURE)
/* conserv_make_bubbles_orig_ver2() will search for the
 * locations of the negative pressure first. Then
 * try to insert bubbles from the location that has
 * the lowest pressure to the location has higher pressure.
 */
LOCAL   int     conserv_make_bubbles(
        Wave                 **wvs,
        Front                **frs,
        int                  num_patches,
        Wv_on_pc             **redistr_table,
        int                  max_n_patch,
        Overparam            *overparam,
        COMPONENT            l_comp,
        COMPONENT            bubblecomp)
{
        int           bubble_loc = 0;
        int           i, j, k, on_patch, on_process;
        Wave          *wave;
        Front         *front;
        INTERFACE     *intfc;
        RECT_GRID     *gr;
        int           n_site = 0, alloc_site = 20;
        INSERT_SITE   *insert_site, *patch_site;
        int           myid, numnodes, *n_site_nodes, max_site;
        static int    First = YES;
        double         refine_ratio, ref_temp;
        static double  spacing;
        int           dim = wvs[0]->rect_grid->dim;
        int           Rbmin[MAXD], Rbmax[MAXD], bmin[MAXD], bmax[MAXD];
        int           imin[MAXD], imax[MAXD], make_bubble;
        static Locstate     ref_st = NULL;
        int           ***distr_site, num_bubbles;

        if(YES == First)
        {
            First = NO;
            bubbles.st = NULL;
            bubbles.bubble = NO;
            bubbles.bubblecomp = bubblecomp;

            front = frs[0];
            intfc = front->interf;
            gr = front->rect_grid;
            if (ref_st == NULL)
            {
                alloc_state(intfc,&ref_st,front->sizest);
                Params(ref_st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_ambient_comp_type(comp_type(bubbles.bubblecomp),front);
                set_type_of_state(ref_st,GAS_STATE);
                zero_state_velocity(ref_st,Params(ref_st)->dim);

                /* values are taken at T =  298.05329710794149 K
                 * P = 1.013 bar
                 * rho = 0.0040948445783132534 g/cm^3
                 * from n-heptane values in gas state
                 */
                Dens(ref_st) = 0.0040948445783132534;
                Energy(ref_st) =  1.013/(1.05 - 1.0);
            }
            if (bubbles.st == NULL)
            {
                alloc_state(intfc,&bubbles.st,front->sizest);
                Params(bubbles.st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_type_of_state(bubbles.st,GAS_STATE);
                zero_state_velocity(bubbles.st,Params(bubbles.st)->dim);
            }
        }
 
        if ( bubbles.bubble != YES && loop_bubble == 0)
        {
            gr = frs[0]->rect_grid;
            bubbles.bubblecomp = bubblecomp;
            set_ambient_comp_type(comp_type(bubbles.bubblecomp),frs[0]);
            bubbles.bubble = YES;
            /* bubble diameter in mesh size */
            bubbles.size = 4;
            /* mini. distance between two bubble boundary in bubble diameter */
            spacing = 50.0;

            refine_ratio = pow(overparam->refinementRatio,
                       overparam->numberOfRefinementLevels-1);
            bubbles.r = 0.98*bubbles.size/2.0*min(gr->h[0]/refine_ratio,
                       gr->h[1]/refine_ratio);
            bubbles.L = 2 * gr->h[0]/refine_ratio;
            bubbles.l = 3 * gr->h[0]/refine_ratio;
            bubbles.N_bubbles = 1;
            bubbles.N_layers = 1;
            bubbles.eos_n = 2;
            bubbles.eos = YES;
            /* Problem specific init. */

            /* Critical pressure to allow bubble creation */
            /* Or can choose to use Clausis-Claperon in case of cavitation */
            bubbles.pr_crit = -10.0;

            /* Again, Clausis-Claperon can be used to set bubble state */
            /*from n-heptane in gas state*/
            bubbles.pr =  1.013;

            /* EOS parameter should not appear. must remove */
            /*from n-heptane in gas state*/
            bubbles.GAMMA = 1.05 - 1.0;
            bubbles.R = 0.83;
        }

        uni_array(&insert_site, alloc_site, sizeof(INSERT_SITE));
        numnodes = pp_numnodes();
        myid = pp_mynode();

        /* Count possible bubble insertion locations */
        for ( i = 0; i < num_patches; i++)
        {
            wave = wvs[i];
            front = frs[i];
            if (wave->patch_level == wave->NumberOfLevels-1)
            {
                bubble_loc = count_bubble_locations(wave,front,l_comp,
                     &alloc_site,&n_site,&insert_site);

                if (debugging("conserv_make_bubbles") && bubble_loc > 0)
                {
                     printf("bubbleLoc on patch[%d][%d] = %d,total_location %d\n",
                        front->patch_number,front->patch_level,bubble_loc,n_site);
                }
            }
        }

        /* put insert_site in order according to pressure */
        uni_array(&n_site_nodes,numnodes,sizeof(int));
        if(numnodes != 1)
        {
#if defined(__MPI__)
            pp_all_gather((POINTER)&n_site,sizeof(int),(POINTER)n_site_nodes,sizeof(int));    
#endif /* if defined(__MPI__) */
        }
        else
            n_site_nodes[0] = n_site;

        max_site = bubble_loc = 0; 
        for(i = 0; i < numnodes; i++)
        {
            /*
            if (debugging("conserv_make_bubbles"))
                printf("node[%d] has insertion locations %d\n", i, n_site_nodes[i]);
            */
            if(n_site_nodes[i] > max_site)
                max_site = n_site_nodes[i];
            bubble_loc += n_site_nodes[i];
        }

        if (debugging("conserv_make_bubbles"))
            printf("Total locations %d, max_site %d\n", bubble_loc, max_site);

        if(max_site == 0)
        {
            free_these(2,n_site_nodes,insert_site);
            return max_site;
        }

        uni_array(&patch_site,bubble_loc,sizeof(INSERT_SITE));
	set_insert_sites(patch_site,insert_site,n_site,max_site,n_site_nodes);
        qsort((POINTER)patch_site, bubble_loc,sizeof(INSERT_SITE),press_ascend);

        /*
        if (debugging("conserv_make_bubbles"))
        {
            for(i = 0; i < bubble_loc; i++)
            {
                printf("loc[%d], pressure %g, ic[%d, %d] on patch [%d]_pic[%d,%d]\n",
                      i, patch_site[i].pressure, patch_site[i].icoords[0],
                      patch_site[i].icoords[1],patch_site[i].on_patch,
                      patch_site[i].on_pc_ic[0],patch_site[i].on_pc_ic[1]);
            }
        }
        */

        /* insert bubbles */
        uni_array(&distr_site,num_patches,sizeof(int**));
        for(j = 0; j < num_patches; j++)
            distr_site[j] = NULL;

	start_clock("make_single_bubble");

        num_bubbles = 0;
        for(i = 0; i < bubble_loc; i++)
        {
	    if(patch_site[i].insert == YES)
	        continue;

            make_bubble = on_process = NO;
            on_patch = -10;
            for(j = 0; j < num_patches; j++)
            {
                front = frs[j];
                wave = wvs[j];
                if (front->patch_level == front->NumberOfLevels-1 &&
                    front->pc_ic[0] == patch_site[i].on_pc_ic[0] &&
                    front->pc_ic[1] == patch_site[i].on_pc_ic[1] &&
                    front->patch_number == patch_site[i].on_patch)
                {
                    on_process = YES;
                    on_patch = j;
                    break;
                }
            } 

            if(on_process == YES)
            {
                gr = front->rect_grid;
                for(k = 0; k < dim; k++)
                {
                    imin[k] = 0;
                    imax[k] = gr->gmax[k];
                }
                if(YES == check_in_bubble_diameter(&(patch_site[i]),
                    wave,front,l_comp,imin,imax,patch_site[i].bmin,
                    patch_site[i].bmax))
                    make_bubble = YES;
                if(YES == make_bubble)
                {
		    start_clock("check_crx_in_bubble_space");	
                    if(YES == check_crx_in_bubble_space(patch_site[i],
                        wave,front,Rbmin,Rbmax,0.0))
                        make_bubble = NO;
		    stop_clock("check_crx_in_bubble_space");	
                }
                if(YES == make_bubble)
                {
                    avg_states_in_reg(&(patch_site[i]),patch_site[i].bmin,
                         patch_site[i].bmax,wave,front,ref_st,
                         &(patch_site[i].num_cells));
                    ref_temp = temperature(ref_st);

                    /* THIS IS NOT A GOOF FIX */
                    if(ref_temp < 294.261)
                        ref_temp = 294.262;

                    patch_site[i].vapor_press = saturated_pres(ref_st,ref_temp);
                }
            }
#if defined(__MPI__)
            pp_global_imax(&make_bubble,1);
#endif /* if defined(__MPI__) */ 
            if(make_bubble == YES)
            {
                if(on_process == YES && debugging("conserv_make_bubbles"))
                {
                    printf("make_bubble[%d], ic[%d, %d] on patch [%d]_pic[%d,%d]\n",
                      i, patch_site[i].icoords[0],
                      patch_site[i].icoords[1],patch_site[i].on_patch,
                      patch_site[i].on_pc_ic[0],patch_site[i].on_pc_ic[1]);
                    /*
                    printf("Average state density %g, mom(%g, %g),"
                           " energy %g, pressure %g, vapor_pressure %g\n",
                      Dens(ref_st),Mom(ref_st)[0],Mom(ref_st)[1],
                      Energy(ref_st),patch_site[i].press, patch_site[i].vapor_press);
                    */
                }
                if(YES == make_bubble_conserv_mass(patch_site,i,wvs,frs,
                      num_patches,redistr_table,max_n_patch,
                      spacing,l_comp,distr_site,on_patch,ref_st,
		      n_site,insert_site))
                {
                    if(on_process == YES)
                        num_bubbles++;
	            set_insert_sites(patch_site,insert_site,n_site,max_site,n_site_nodes);
                    qsort((POINTER)patch_site, bubble_loc,sizeof(INSERT_SITE),press_ascend);
                }
            }
        }

	stop_clock("make_single_bubble");

        for(j = 0; j < num_patches; j++)
        {
            if(distr_site[j] != NULL)
                free(distr_site[j]);
        }
        free(distr_site);

        free(patch_site);
        free(n_site_nodes);
        free(insert_site);

        return num_bubbles;
}

LOCAL int make_bubble_conserv_mass(
        INSERT_SITE      *patch_site,
        int              loc,
        Wave             **wvs,
        Front            **frs,
        int              num_patches,
        Wv_on_pc         **redistr_table,
        int              max_n_patch,
        double            spacing,
        COMPONENT        l_comp,
        int              ***in_distr_site,
        int              on_patch,
        Locstate         ref_st,
	int              n_site,
	INSERT_SITE      *insert_site)
{
        RECT_GRID    *gr = wvs[0]->rect_grid;
        int          dim = wvs[0]->rect_grid->dim, j, k;
        int          gic[MAXD], item[2], Rbmin[MAXD], Rbmax[MAXD];
        int          gRbmin[MAXD], gRbmax[MAXD], it[2];
        int          grid, gmin[MAXD], gmax[MAXD], iL[MAXD], iU[MAXD];
        int          lRbmin[MAXD], lRbmax[MAXD], ix, iy;
        Front        *front;
        Wave         *wave;
        int          **distr_site, def_cells, all_def_cells;
        double        D_add, ref_temp, all_mass_add, alpha = 1.0;
        static Locstate rst = NULL, distr_st;
        Locstate     st;
        int          create_bubble = NO;
        static byte  *send_bub_st;
        int          all_num_cell;
	double        *crds;

        if(rst == NULL)
        {
            front = frs[0];
            alloc_state(front->interf,&rst,front->sizest);
            alloc_state(front->interf,&distr_st,front->sizest);
            set_type_of_state(rst,GAS_STATE);
            set_type_of_state(distr_st,GAS_STATE);
            scalar(&send_bub_st,wvs[0]->sizest);
        }

        locate_on_redistr_table(patch_site[loc].on_pc_ic,
              patch_site[loc].on_patch,redistr_table,max_n_patch,item);
        blks_bubble_ocupy(patch_site[loc],gr,Rbmin,Rbmax,spacing);
        
        if(on_patch >= 0)
        {
            create_bubble = YES;
            front = frs[on_patch];
            wave = wvs[on_patch];
        }

        /* Create vapor bubble here */
        if(create_bubble == YES)
        {
            ref_temp = temperature(ref_st);

            /* THIS IS NOT A GOOF FIX */
            if(ref_temp < 294.261)
                ref_temp = 294.262;

            set_state_on_phase_boundary_ver2(ref_st,rst,ref_temp);
            (*wave->bundle_single_st)(bubbles.st,wave,send_bub_st);
        }

#if defined(__MPI__)
        pp_bcast(redistr_table[item[0]][item[1]].pc_id,send_bub_st,wvs[0]->sizest);

        if(pp_mynode() != redistr_table[item[0]][item[1]].pc_id)
        {
            (*wvs[0]->unbundle_single_st)(bubbles.st,wvs[0],send_bub_st);
            patch_site[loc].vapor_press = pressure(bubbles.st);
        }
#endif /* if defined(__MPI__) */
 
        /* if(debugging("conserv_make_bubbles")) */
        {
            copy_state(distr_st, ref_st);
            /* verbose_print_state("ref_st",ref_st); */
            if(create_bubble == YES)
                (*wvs[0]->bundle_single_st)(distr_st,wvs[0],send_bub_st);  
#if defined(__MPI__)
            pp_bcast(redistr_table[item[0]][item[1]].pc_id,send_bub_st,wvs[0]->sizest);
            if(pp_mynode() != redistr_table[item[0]][item[1]].pc_id)
                (*wvs[0]->unbundle_single_st)(distr_st,wvs[0],send_bub_st);
#endif /* if defined(__MPI__) */
            if (gr->Remap.remap == CYLINDRICAL_REMAP && is_rotational_symmetry())
                alpha = patch_site[loc].center[0];
            /* verbose_print_state("distr_st",distr_st); */
            Dens(distr_st) *= alpha*pi*sqr(bubbles.r);
            Mom(distr_st)[0] *= alpha*pi*sqr(bubbles.r);
            Mom(distr_st)[1] *= alpha*pi*sqr(bubbles.r);
            Energy(distr_st) *= alpha*pi*sqr(bubbles.r);
            /* printf("alpha of center %g, area %g\n", alpha, sqr(bubbles.r)); */
        }

        /* Create vapor bubble here */
        if(create_bubble == YES)
        {
            create_bubbles_ver2(&(patch_site[loc]),l_comp,wave,front,rst,
                  patch_site[loc].bmin,patch_site[loc].bmax);
            patch_site[loc].insert = YES;
	    for(k = 0; k < n_site; k++)
	    {
	        if(insert_site[k].on_pc_ic[0] == front->pc_ic[0] &&
	           insert_site[k].on_pc_ic[1] == front->pc_ic[1] && 
		   insert_site[k].on_patch == front->patch_number &&
		   insert_site[k].icoords[0] == patch_site[loc].icoords[0] &&
		   insert_site[k].icoords[1] == patch_site[loc].icoords[1])
		    insert_site[k].insert = YES;	
                if(insert_site[k].on_pc_ic[0] == front->pc_ic[0] &&
		   insert_site[k].on_pc_ic[1] == front->pc_ic[1] &&
		   insert_site[k].on_patch == front->patch_number &&
		   insert_site[k].insert == NO)
		{
                    for(iy = patch_site[loc].bmin[1]; iy <= patch_site[loc].bmax[1]; iy++) 
		    {
		        for(ix = patch_site[loc].bmin[0]; ix <= patch_site[loc].bmax[0]; ix++)
		        {
		   	    gic[0] = ix; gic[1] = iy;    
		            crds = Rect_coords(gic,wave);	     
			    if((sqr(patch_site[loc].center[0]-crds[0]) + 
			        sqr(patch_site[loc].center[1]-crds[1])) >
	                            sqr(bubbles.r))
                                continue;
			    if(insert_site[k].icoords[0] == gic[0] &&
			       insert_site[k].icoords[1] == gic[1])
			    {
		                insert_site[k].insert = YES;		 
			        break;
			    }
		        }
			if(insert_site[k].insert == YES)
			    break;	
		    }
		} 
	    }
        }

reduce_spacing_conservation:
        for(k = 0; k < dim; k++)
        {
            gic[k] = patch_site[loc].icoords[k] +
                     redistr_table[item[0]][item[1]].base[k] +
                     redistr_table[item[0]][item[1]].off_set[k];
            gRbmin[k] = Rbmin[k] + redistr_table[item[0]][item[1]].base[k] +
                        redistr_table[item[0]][item[1]].off_set[k];
            gRbmax[k] = Rbmax[k] + redistr_table[item[0]][item[1]].base[k] +
                        redistr_table[item[0]][item[1]].off_set[k];
        }
        
        if(debugging("conserv_make_bubbles") && create_bubble == YES)
        {
            printf("step[%d] create_bubble_at[%d %d][%g %g], mass_conserv_in"
                   " X[%d %d],Y[%d %d] size[%d, %d]\n",
                 frs[0]->step,gic[0],gic[1],patch_site[loc].center[0],
                 patch_site[loc].center[1],
                  gRbmin[0],gRbmax[0],gRbmin[1],
                 gRbmax[1],gRbmax[0]-gRbmin[0]+1,gRbmax[1]-gRbmin[1]+1);
            verbose_print_state("bubbles.st",bubbles.st);
            /* verbose_print_state("ref_st",ref_st); */
            printf("Average state density %g, mom(%g, %g),"
                   " energy %g, pressure %g, vapor_pressure %g\n",
            Dens(ref_st),Mom(ref_st)[0],Mom(ref_st)[1],
            Energy(ref_st),patch_site[loc].press, patch_site[loc].vapor_press);
        }

        all_mass_add = 0.0; all_def_cells = 0;
        for(grid = 0; grid < num_patches; grid++)
        {
            if(frs[grid]->patch_level != frs[grid]->NumberOfLevels-1)
                continue;
            front = frs[grid];
            wave = wvs[grid];
            locate_on_redistr_table(front->pc_ic,
              front->patch_number,redistr_table,max_n_patch,it);

            for(k = 0; k < dim; k++)
            {
                gmin[k] = redistr_table[it[0]][it[1]].base[k] +
                            redistr_table[it[0]][it[1]].off_set[k];
                gmax[k] = redistr_table[it[0]][it[1]].bound[k] +
                            redistr_table[it[0]][it[1]].off_set[k];
            }

            /*
            printf("front[%d][%d] on pc[%d %d] item[%d %d] with size[%d %d][%d %d],"
                    " base[%d %d] off_set[%d %d]\n",
                 front->patch_number, front->patch_level, front->pc_ic[0], 
                 front->pc_ic[1],
                 it[0], it[1], gmin[0], gmin[1], gmax[0], gmax[1],
                   redistr_table[it[0]][it[1]].base[0],
                   redistr_table[it[0]][it[1]].base[1],
                   redistr_table[it[0]][it[1]].off_set[0],
                   redistr_table[it[0]][it[1]].off_set[1]);
            */

            if(YES != box_intersect(gRbmin,gRbmax,gmin,gmax,dim,iL,iU))
                continue;

            /*
            printf("front[%d][%d] on pc[%d %d] item[%d %d] with base[%d %d]"
                   " sizeX-[%d %d]Y-[%d %d]," " intersect with box\n",
                 front->patch_number, front->patch_level, front->pc_ic[0], 
                 front->pc_ic[1],
                 it[0], it[1], redistr_table[it[0]][it[1]].base[0], 
                 redistr_table[it[0]][it[1]].base[1],
                 gmin[0], gmin[1], gmax[0], gmax[1]);
            */

            for(k = 0; k < dim; k++)
            {
                lRbmin[k] = gRbmin[k] - (redistr_table[it[0]][it[1]].base[k] +
                           redistr_table[it[0]][it[1]].off_set[k]);
                lRbmax[k] = gRbmax[k] - (redistr_table[it[0]][it[1]].base[k] +
                           redistr_table[it[0]][it[1]].off_set[k]);
            }
            if(in_distr_site[grid] == NULL)
            {
                bi_array(&distr_site,front->rect_grid->gmax[0],
                       front->rect_grid->gmax[1],sizeof(int));
                in_distr_site[grid] = distr_site;
                for(iy = 0; iy < front->rect_grid->gmax[1]; iy++)
                    for(ix = 0; ix < front->rect_grid->gmax[0]; ix++)
                        distr_site[ix][iy] = NO;
            }
            else
                distr_site = in_distr_site[grid];

            def_cells = neighbor_mass_deficiency(patch_site[loc],l_comp,
                  lRbmin,lRbmax,wave,front,distr_site,&D_add);

            all_mass_add += D_add;
            all_def_cells += def_cells;
            /*
            printf("front[%d][%d] on pc[%d %d] with size[%d %d] has %d deficient cells, D_add %g\n",
                 front->patch_number, front->patch_level, front->pc_ic[0], front->pc_ic[1],
                 gmax[0]-gmin[0], gmax[1]-gmin[1], def_cells, D_add);
            */
        }

#if defined(__MPI__)
        pp_global_sum(&all_mass_add,1);
        pp_global_isum(&all_def_cells,1);
#endif /* if defined(__MPI__) */
        all_num_cell = (Rbmax[0]-Rbmin[0]+1)*(Rbmax[1]-Rbmin[1]+1);

        if(all_mass_add/Dens(distr_st) > 1.0)
        {
            spacing *= 0.9;
            blks_bubble_ocupy(patch_site[loc],gr,Rbmin,Rbmax,spacing);
            goto reduce_spacing_conservation;
        }
        else
        {
            for(grid = 0; grid < num_patches; grid++)
            {
                if(frs[grid]->patch_level != frs[grid]->NumberOfLevels-1)
                    continue;
                front = frs[grid];
                wave = wvs[grid];
                locate_on_redistr_table(front->pc_ic,
                  front->patch_number,redistr_table,max_n_patch,it);

                for(k = 0; k < dim; k++)
                {
                    gmin[k] = redistr_table[it[0]][it[1]].base[k] +
                            redistr_table[it[0]][it[1]].off_set[k];
                    gmax[k] = redistr_table[it[0]][it[1]].bound[k] +
                            redistr_table[it[0]][it[1]].off_set[k];
                }

                if(YES != box_intersect(gRbmin,gRbmax,gmin,gmax,dim,iL,iU))
                    continue;

                for(k = 0; k < dim; k++)
                {
                    lRbmin[k] = gRbmin[k] - (redistr_table[it[0]][it[1]].base[k] +
                           redistr_table[it[0]][it[1]].off_set[k]);
                    lRbmax[k] = gRbmax[k] - (redistr_table[it[0]][it[1]].base[k] +
                           redistr_table[it[0]][it[1]].off_set[k]);
                }
                distr_site = in_distr_site[grid];

                distr_mass_to_neighbors(patch_site[loc],l_comp,
                    patch_site[loc].num_cells,lRbmin,lRbmax,wave,front,
		    distr_site);
		for(k = 0; k < n_site; k++)
		{
		    if(insert_site[k].on_pc_ic[0] == front->pc_ic[0] &&
		       insert_site[k].on_pc_ic[1] == front->pc_ic[1] && 
		       insert_site[k].on_patch == front->patch_number)
		    {
		        if(distr_site[insert_site[k].icoords[0]][insert_site[k].icoords[1]] 
					== YES)
			    insert_site[k].insert = YES;	
		    }
		}
            }
        }
 
        /* if(debugging("conserv_make_bubbles")) */
        {
            /* if(create_bubble == YES) */
            printf("step[%d]: %g mass conserved[%g%] out of %g(avg %g);"
                       " To %d cells(%g% Total)\n", 
                     frs[0]->step,all_mass_add, 
                     100.0*all_mass_add/Dens(distr_st), Dens(distr_st),
                     Dens(distr_st)/(alpha*pi*sqr(bubbles.r)),
                     all_def_cells, 100.0*all_def_cells/all_num_cell);
        }
        return YES;
}

LOCAL void distr_mass_to_neighbors(
      INSERT_SITE  insert_site,
      COMPONENT    l_comp,
      int          n_cells,
      int          *Rbmin,
      int          *Rbmax,
      Wave         *wave,
      Front        *front,
      int          **distr_site)
{
      int          x_size, y_size, i, j, k, side_y, side_x;
      int          ic[MAXD], icent[MAXD];
      Locstate     st;
      double        D, T, T_sat, P_sat, pres, d_sat, c_sqr, vapor_press, p_jump;
      double        *crds, area, alpha = 1.0, d_dif, s_alpha = 0.0;
      double        avail, d_incr, m0_incr, m1_incr, e_incr;
      double        in_vel0, in_vel1, n_vel0, n_vel1, vel0, vel1, D_add;
      int          num_cell = 0, count = 0, all_cell, more_count = 0;
      double        total_add_mass = 0.0;

      if (wave->rect_grid->Remap.remap == CYLINDRICAL_REMAP && is_rotational_symmetry())
          alpha = insert_site.center[0];

      vapor_press = insert_site.vapor_press;

      x_size = ((Rbmax[0] - Rbmin[0]) - (Rbmax[0] - Rbmin[0])%2)/2 + 1;
      y_size = ((Rbmax[1] - Rbmin[1]) - (Rbmax[1] - Rbmin[1])%2)/2 + 1;
      icent[0] = ((Rbmax[0] + Rbmin[0]) + (Rbmax[0] + Rbmin[0])%2)/2;
      icent[1] = ((Rbmax[1] + Rbmin[1]) + (Rbmax[1] + Rbmin[1])%2)/2;
      area = wave->rect_grid->h[0]*wave->rect_grid->h[1];

      for(side_y = -1; side_y <=1; side_y+=2)
      { 
          int starty;
          if(side_y == -1) starty = 0; 
          else
              starty = 1;   
          for(j = starty; j <= y_size; j++)
          {
              int startx;
              ic[1] = icent[1] + side_y*j;  
              for(side_x = -1; side_x <=1; side_x+=2)
              {
                  if(side_x == -1) startx = 0; 
                  else
                     startx = 1;
                  for(i = startx; i <= x_size; i++)
                  {
                      ic[0] = icent[0] + side_x*i;  
                      if(ic[0] < Rbmin[0] || ic[0] > Rbmax[0] ||
                         ic[1] < Rbmin[1] || ic[1] > Rbmax[1])
                          continue;
                      if(ic[0] >= wave->rect_grid->gmax[0] ||
                         ic[1] >= wave->rect_grid->gmax[1] ||
                         ic[0] < 0 || ic[1] < 0)
                          continue;
                      if(Rect_comp(ic,wave) != l_comp)
                          continue;
                      /* Each cell is distributed once */
                      if(distr_site[ic[0]][ic[1]] == YES)
                          continue;

                      crds = Rect_coords(ic,wave);
                      st = Rect_state(ic,wave);

                      T = temperature(st);
                      pres = pressure(st);
                      P_sat = saturated_pres(st,T);

#if defined(USE_ADIABAT_PRES)
                      if(pres > vapor_press)
                          continue;
#else /* if defined(USE_ADIABAT_PRES) */
                      if(pres > P_sat)
                          continue;
#endif /* if defined(USE_ADIABAT_PRES) */

                      if (wave->rect_grid->Remap.remap == CYLINDRICAL_REMAP &&
				is_rotational_symmetry())
                          alpha = crds[0];

                      /* if(debugging("conserv_bubble") && */
		      /*    fabs(0.00144625-crds[0]) < 0.0001 && */
	              /*  fabs(0.168544-crds[1]) < 0.0001) */
                      if(debugging("conserv_bubble"))
                      {
                          /***
                          printf("cell [%d %d] comp[%d] before add states\n", 
                              ic[0], ic[1], Rect_comp(ic,wave));
                          verbose_print_state("before_add", st);
                          ***/
                      }
                      D = Dens(st);
                      T = temperature(st);
                      c_sqr = sound_speed_squared(st);
                      vel0 = vel(0,st);
                      vel1 = vel(1,st);

                      distr_site[ic[0]][ic[1]] = YES;
                      
#if defined(USE_ADIABAT_PRES)
                      set_state(st,TGAS_STATE,st);
                      pres = Press(st);
                      Press(st) = vapor_press;
                      p_jump = vapor_press - pres;

                      /* adiabatically add density */
                      D_add = p_jump/c_sqr;
                      Dens(st) += D_add;

                      /* Mom increment */
                      /*
                      n_vel0 = (vel0*D + D_add*in_vel0)/Dens(st);
                      n_vel1 = (vel1*D + D_add*in_vel1)/Dens(st);
                      Vel(st)[0] = n_vel0;
                      Vel(st)[1] = n_vel1;
                      */
                      set_state(st,GAS_STATE,st);
                      total_add_mass +=  D_add*alpha*area;
                      num_cell++;
#endif /* if defined(USE_ADIABAT_PRES) */

                      /* if(debugging("conserv_bubble") && */
		      /*    fabs(0.00144625-crds[0]) < 0.0001 && */
	              /* fabs(0.168544-crds[1]) < 0.0001) */
                      if(debugging("conserv_bubble"))
                      { 
                          /* printf("cell [%d %d] after add states\n", ic[0], ic[1]); */
                          /* verbose_print_state("after_add", st); */
                      }
                  }
              }
          }
      }
}
 

LOCAL   int locate_on_redistr_table(
        int     *on_pc_ic,
        int      on_patch,
        Wv_on_pc   **redistr_table,
        int        max_n_patch,
        int        *item)
{
        int          numnodes;
        int          i, src;
        numnodes = pp_numnodes();
                                                                                                     
        for(src = 0; src < numnodes; src++)
        {
            for(i = 0; i < max_n_patch; i++)
            {
                if(redistr_table[src][i].wv_id == -1 ||
                   redistr_table[src][i].wv_id != on_patch)
                    continue;
                if(redistr_table[src][i].pc_ic[0] == on_pc_ic[0] &&
                   redistr_table[src][i].pc_ic[1] == on_pc_ic[1])
                {
                    item[0] = src;
                    item[1] = i;
                    return YES;
                }
            }
        }
        return NO;
}

LOCAL   int     count_bubble_locations(
        Wave                *wave,
        Front               *front,
        COMPONENT           l_comp,
        int                 *alloc_site,
        int                 *n_site,
        INSERT_SITE         **insert_site)
{
        RECT_GRID           *gr = front->rect_grid;
        int                 dim = gr->dim;
        int                 icoords[MAXD],imin[MAXD],imax[MAXD];
        int                 bmin[MAXD], bmax[MAXD];
        int                 ix,iy, i, j;
        Locstate            state;
        double               T_sat, P_sat, press, ref_temp;
        int                 n_loc = 0;
        double               dt, c_speed;
        int                 make_bubble, num_cells;
        static Locstate     ref_st = NULL, rst;
	static	double	    templ,temph;

        dt = front->dt;
        for(i = 0; i < dim; i++)
        {
            imin[i] = 0;
            imax[i] = gr->gmax[i];
        }

        if (ref_st == NULL)
        {
            alloc_state(front->interf,&ref_st,front->sizest);
            alloc_state(front->interf,&rst,front->sizest);
            Params(ref_st) = ((G_INTERFACE *) front->interf)->g_user_intfc.params_list[2];
            Params(rst) = ((G_INTERFACE *) front->interf)->g_user_intfc.params_list[2];
            set_type_of_state(ref_st,GAS_STATE);
            set_type_of_state(rst,GAS_STATE);
            zero_state_velocity(ref_st,Params(ref_st)->dim);
                                                                                                           
            /* values are taken at T =  298.05329710794149 K
             * P = 1.013 bar
             * rho = 0.0040948445783132534 g/cm^3
             * from n-heptane values in gas state
             */
            Dens(ref_st) = 0.0040948445783132534;
            Energy(ref_st) =  1.013/(1.05 - 1.0);
	    phase_temp_range(&templ,&temph);
        }

        /* Search pressure below specified critical value */
        for (iy = imin[1]; iy < imax[1]; iy++)
        {
            icoords[1] = iy;
            for (ix = imin[0]; ix < imax[0]; ix++)
            {
                icoords[0] = ix;
                if (Rect_comp(icoords,wave) != l_comp) continue;
                make_bubble = NO;
                state = Rect_state(icoords,wave);
                T_sat = temperature(state);
                /* Temperature range test */
                if(T_sat < templ || T_sat > temph)
                {
                    printf("WARNING: in count_bubble_locations(),"
                         " T-sat temp is outside range[294.261 533.15]\n");
                    printf("Center[%g %g][%d %d] , center to ax %g comp[%d]\n",
                        Rect_coords(icoords,wave)[0],
                        Rect_coords(icoords,wave)[1],
                        icoords[0], icoords[1],
                        Rect_coords(icoords,wave)[0]/gr->h[0],
			Rect_comp(icoords,wave));
#if defined(USE_OVERTURE)
		    printf("On patch[%d][%d]\n", wave->patch_number,wave->patch_level);
#endif /* if defined(USE_OVERTURE) */
                    verbose_print_state("state",state);
		    /*
		    print_rectangular_grid(front->rect_grid);
		    geomview_intfc_plot2d("gview_plot",front->step,
		         front->interf,front->rect_grid);
                    clean_up(ERROR);
		    */
		    /* THIS IS NOT A GOOF FIX */
		    if(T_sat < 294.261)
		        T_sat = 294.262;
                }
                P_sat = saturated_pres(state,T_sat);
                press = pressure(state);

                /* from homo model */
                /* if (press < -5.0) */
                /* Clausis-Clapeyron */
                if (P_sat > press + 0.18)  
                /* from \delta P = 2*surface_tension/radius. */
                /* Here radius=2*1.1125 microns, s_tension=1.96e-2 kg/s^2 */
                {
                    c_speed = sound_speed(state);
                    (*insert_site)[*n_site].rad_distr = c_speed*dt;
                    ft_assign((*insert_site)[*n_site].icoords,icoords,INT*MAXD);
                    (*insert_site)[*n_site].press = press;
                    (*insert_site)[*n_site].insert = NO;
                    (*insert_site)[*n_site].on_patch = front->patch_number;
                    for(i = 0; i < dim; i++)
                        (*insert_site)[*n_site].on_pc_ic[i] = front->pc_ic[i];

                    if(YES == check_in_bubble_diameter(&((*insert_site)[*n_site]),
                        wave,front,l_comp,imin,imax,bmin,bmax))
                        make_bubble = YES;
                    if(make_bubble == YES)
                    {
                        /* Geometric specific */
                        /* Do not create bubbles close to nozzle wall */
                        if((fabs((*insert_site)[*n_site].center[0]-0.0089) < 
                                 (bubbles.r + 2.0*gr->h[0])) &&
                            (*insert_site)[*n_site].center[1] < 0.15  &&
                            (*insert_site)[*n_site].center[1] > 0.05) 
                            make_bubble = NO;
                    }
                    if(make_bubble == YES)
                    {
                        (*n_site)++;
                        n_loc++;
                    }
                    if(*n_site == *alloc_site)
                    {
                        INSERT_SITE *tmp_site;
                        (*alloc_site) += max(imax[0],imax[1]);
                        uni_array(&tmp_site, (*alloc_site), sizeof(INSERT_SITE));
                        for(i = 0; i < (*n_site); i++)
                            ft_assign(&tmp_site[i], &((*insert_site)[i]), sizeof(INSERT_SITE));
                        free(*insert_site);
                        *insert_site = tmp_site;
                    }
                }
            }
        }
        return n_loc;
}

LOCAL   int     conserv_make_bubbles_orig_ver2(
        Wave                    *wave,
        Front                   *front,
        COMPONENT               l_comp,
        COMPONENT               bubblecomp)
{
        RECT_GRID               *gr = front->rect_grid;
        int                     dim = gr->dim;
        int                     icoords[MAXD],imin[MAXD],imax[MAXD];
        int                     bmin[MAXD], bmax[MAXD], Rbmin[MAXD], Rbmax[MAXD];
        int                     ix,iy, i, j;

        INTERFACE               *intfc = front->interf;
        int                     size; 
        double                   spacing;
        int                     num_bubble = 0;
        double                   period = 5.e-6;
        Locstate                state;
        double                   T_sat, P_sat, press, ref_temp;
        static Locstate         ref_st = NULL, rst = NULL, distr_st;
        static int              First = YES;
        INSERT_SITE             *insert_site = NULL;        
        int                     n_site = 0, alloc_site = 0;
        int                     num_cells, **distr_site;
        double                   alpha = 1.0,D_add;
	stataic	double		templ,temph;

        debug_print("bubbles_mpi","Entering conserv_make_bubbles_orig_ver2()\n"); 
        if(YES == First)
        {
            bubbles.st = NULL;
            bubbles.bubble = NO;
            First = NO;
	    phase_temp_range(&templ,&temph);
        }

 
        bubbles.bubblecomp = bubblecomp;
        if ( bubbles.bubble != YES && loop_bubble == 0)
        {
            if (ref_st == NULL)
            {
                alloc_state(intfc,&ref_st,front->sizest);
                alloc_state(intfc,&rst,front->sizest);
                alloc_state(intfc,&distr_st,front->sizest);
                Params(ref_st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                Params(rst) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_ambient_comp_type(comp_type(bubbles.bubblecomp),front);
                set_type_of_state(ref_st,GAS_STATE);
                set_type_of_state(rst,GAS_STATE);
                zero_state_velocity(ref_st,Params(ref_st)->dim);

                /* values are taken at T =  298.05329710794149 K
                 * P = 1.013 bar
                 * rho = 0.0040948445783132534 g/cm^3
                 * from n-heptane values in gas state
                 */
                Dens(ref_st) = 0.0040948445783132534;
                Energy(ref_st) =  1.013/(1.05 - 1.0);
            }
            if (bubbles.st == NULL)
            {
                alloc_state(intfc,&bubbles.st,front->sizest);
                Params(bubbles.st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_type_of_state(bubbles.st,GAS_STATE);
                zero_state_velocity(bubbles.st,Params(bubbles.st)->dim);
            }

            bubbles.bubble = YES;
            /* bubble diameter in mesh size */
            bubbles.size = size = 3;
            /* mini. distance between two bubble boundary in bubble diameter */
            spacing = 40.0;

            bubbles.r = 0.98*size/2.0*min(gr->h[0], gr->h[1]);
            bubbles.N_bubbles = 1;
            bubbles.N_layers = 1;
            bubbles.L = 2 * gr->h[0];
            bubbles.l = 3 * gr->h[0];
            bubbles.eos_n = 2;
            bubbles.eos = YES;

            /* Problem specific init. */

            /* Critical pressure to allow bubble creation */
            /* Or can choose to use Clausis-Claperon in case of cavitation */
            bubbles.pr_crit = -10.0;

            /* Again, Clausis-Claperon can be used to set bubble state */
            /*from n-heptane in gas state*/
            bubbles.pr =  1.013;

            state_on_adiabat_with_pr(ref_st, bubbles.pr, rst, GAS_STATE);
            bubbles.rho = Dens(rst);

            /* EOS parameter should not appear. must remove */
            /*from n-heptane in gas state*/
            bubbles.GAMMA = 1.05 - 1.0;
            bubbles.R = 0.83;
        }        

        if  (bubbles.bubble == YES)
        {
            /* domain which allows bubble to be created */
            /* The current domain setting might give insufficient
             * neighboring space if a bubble is created near the boundary, when
             * the spacing is greater than 1 bubble diameter 
             * and bubble diameter is greater than buffer size (4).
             * However, to include these in the domain setting might never
             * create bubbles near the boundary. 
             */
            for(i = 0; i < dim; i++)
            {
                imin[i] = 0;
                imax[i] = gr->gmax[i]; 
                /** For conservative insertion
                if(rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
                    imin[i] -= gr->lbuf[i];
                if(rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
                    imax[i] += gr->ubuf[i];
                **/
            }

            alloc_site = max(imax[0],imax[1]); 
            uni_array(&insert_site, alloc_site, sizeof(INSERT_SITE)); 
            bi_array(&distr_site, gr->gmax[0], gr->gmax[1], sizeof(int));
            for(iy = 0; iy < gr->gmax[1]; iy++)
                for(ix = 0; ix < gr->gmax[0]; ix++)
                    distr_site[ix][iy] = NO;

            /* Search pressure below specified critical value */
            for (iy = imin[1]; iy < imax[1]; iy++)
            {
                icoords[1] = iy;
                for (ix = imin[0]; ix < imax[0]; ix++)
                {
                    icoords[0] = ix;
                    if (Rect_comp(icoords,wave) != l_comp) continue;

                    state = Rect_state(icoords,wave);
                    T_sat = temperature(state);
                    /* Temperature range test */
                    if(T_sat < templ || T_sat > temph)
                    {
                        printf("WARNING: in make_bubbles_orig_ver2,"
                             " T-sat temp is outside range[294.261 533.15]\n");
                        printf("Center[%g %g][%d %d] , center to ax %g\n",
                            Rect_coords(icoords,wave)[0], 
                            Rect_coords(icoords,wave)[1],
                            icoords[0], icoords[1], 
                            Rect_coords(icoords,wave)[0]/gr->h[0]);
                        verbose_print_state("state",state);
                        /* clean_up(ERROR); */
                    }
                    P_sat = saturated_pres(state,T_sat);
                    press = pressure(state);

                    /* from homo model */
                    /* if (press < -10.0) */
                    /* Clausis-Clapeyron */
                    if (P_sat > press + 0.1)
                    {
                        ft_assign(insert_site[n_site].icoords,icoords,INT*MAXD); 
                        insert_site[n_site].press = press; 
                        insert_site[n_site].insert = NO;
                        n_site++;
                        if(n_site == alloc_site)
                        {
                            INSERT_SITE *tmp_site; 
                            alloc_site += max(imax[0],imax[1]);
                            uni_array(&tmp_site, alloc_site, sizeof(INSERT_SITE));
                            for(i = 0; i < n_site; i++)
                                ft_assign(&tmp_site[i], &insert_site[i], sizeof(INSERT_SITE));
                            free(insert_site);
                            insert_site = tmp_site;
                        }
                    }
                }
            }
            
            /* put insert_site in order according to pressure */
            if(n_site != 0)
                qsort((POINTER)insert_site, n_site,sizeof(INSERT_SITE),press_ascend);

            if (debugging("bubbles_mpi"))
            {
#if defined(USE_OVERTURE)
                (void) printf("On patch[%d], level[%d]\n", 
                          wave->patch_number, wave->patch_level);
#endif /* if defined(USE_OVERTURE) */
                (void) printf("domain size imin[%d %d], imax[%d %d]\n",
                             imin[0], imin[1], imax[0], imax[1]);
                (void) printf("Find %d possible insertion locations\n", n_site);
                for(i = 0; i < n_site; i++)
                    (void)printf("At coords[%g, %g] icoords[%d, %d], pressure %g\n",
                        Rect_coords(insert_site[i].icoords,wave)[0],
                        Rect_coords(insert_site[i].icoords,wave)[1],
                        insert_site[i].icoords[0],insert_site[i].icoords[1],
                        insert_site[i].press);
            }

            /* Create bubbles */
            for(i = 0; i < n_site; i++)
            {
                if(YES != check_in_bubble_diameter(&(insert_site[i]), 
                    wave,front,l_comp,imin,imax,bmin,bmax))
                    continue;

                /* Geometric specific */
                /* Do not create bubbles close to nozzle wall */
                if((fabs(insert_site[i].center[0]-0.0089) < (bubbles.r + gr->h[0])) &&
                    insert_site[i].center[1] < 0.15  && 
                    insert_site[i].center[1] > 0.5)
                    continue;

                blks_bubble_ocupy(insert_site[i],gr,Rbmin,Rbmax,spacing);

                /*** Old alg. bubbles are always created in spacing distance to curves
                if(YES == check_crx_in_bubble_space(insert_site[i], 
                        wave,front,Rbmin,Rbmax,spacing))
                    continue;
                ***/

                /*** New alg.: Only bubbles created at the same time must 
                 *** have the spacing. Bubbles can be created right next previous 
                 *** time step bubbles.
                 ***/
                /* The conservative insertion does not need restrict spacing at 
                 * the same time step. new_bubbles_in_domain() is not needed.
                 */
                /* if(YES == new_bubbles_in_domain(insert_site,i-1,Rbmin,Rbmax)) */
                /*     continue;                */
                if(YES == check_crx_in_bubble_space(insert_site[i], 
                        wave,front,Rbmin,Rbmax,0.0))
                    continue;
                /*** End of New alg. ***/

                avg_states_in_reg(&(insert_site[i]),bmin,bmax,wave,front,ref_st,&num_cells);

                /** Compute total quantity needs to be distributed **/
                copy_state(distr_st, ref_st);

                if (gr->Remap.remap == CYLINDRICAL_REMAP &&
			is_rotational_symmetry())
                    alpha = insert_site[i].center[0];
                Dens(distr_st) *= alpha*pi*sqr(bubbles.r);
                Mom(distr_st)[0] *= alpha*pi*sqr(bubbles.r);
                Mom(distr_st)[1] *= alpha*pi*sqr(bubbles.r);
                Energy(distr_st) *= alpha*pi*sqr(bubbles.r);

                /* The reference temperature will be used to determine  */
                /* vapor bubble states and init. phase bondary state. */

                ref_temp = temperature(ref_st);
		/* temperature range test */
		if(ref_temp < templ || ref_temp > temph)
		{
	            printf("ERROR: in make_bubbles_orig_ver2,"
		         " ref_st temp %g is outside range[294.261 533.15]\n", ref_temp);
                    printf("Center[%g %g][%d %d] ,%d cells: center to ax %g\n",
                        insert_site[i].center[0], insert_site[i].center[1], 
			insert_site[i].icoords[0],
			insert_site[i].icoords[1], num_cells,
                        insert_site[i].center[0]/gr->h[0]);
		    verbose_print_state("ref_st", ref_st);
		    /* clean_up(ERROR); */
		}
                set_state_on_phase_boundary_ver2(ref_st, rst, ref_temp); 

                if(debugging("conserv_bubble"))
                {
                    printf("Center[%g %g],%d cells: center to ax %g\n",
                        insert_site[i].center[0], insert_site[i].center[1], num_cells,
                        insert_site[i].center[0]/gr->h[0]);

                    verbose_print_state("Vapor bubble state", bubbles.st);

                    printf("Removed: Den %g, Mom[0] %g, Mom[1] %g, Eng %g\n",
                        Dens(distr_st), Mom(distr_st)[0], 
                        Mom(distr_st)[1], Energy(distr_st));
                    printf("Energy consumed to create bubble %g g.cm/msec^2, distr_engy %g\n", 
                        1.96e-5*2.0*pi*bubbles.r, Energy(distr_st)/alpha);
                    printf("Total heat used %g\n", 
                         Dens(bubbles.st)*pi*sqr(bubbles.r)*3619); 
                    verbose_print_state("avg_state", ref_st);
                }
                else if(debugging("bubble_insert"))
                {
                    printf("Center[%g %g],%d cells: center to ax %g\n",
                        insert_site[i].center[0], insert_site[i].center[1], num_cells,
                        insert_site[i].center[0]/gr->h[0]);
                    verbose_print_state("Vapor bubble state", bubbles.st);
                }
            
                create_bubbles_ver2(&(insert_site[i]),l_comp,wave,front,rst,bmin,bmax);

                insert_site[i].insert = YES;
                for(j = 0; j < dim; j++)
                {
                    insert_site[i].bmin[j] = bmin[j]; 
                    insert_site[i].bmax[j] = bmax[j]; 
                }
                num_bubble++;

                /* Distribute quantity to neighbors */
                neighbor_mass_deficiency(insert_site[i],l_comp,
                      Rbmin,Rbmax,wave,front,distr_site,&D_add);
                neighbor_mass_deficiency_on_curve(insert_site[i],l_comp,distr_st,
                      num_cells,Rbmin,Rbmax,wave,front,distr_site);
                distr_to_neighbors(insert_site[i],l_comp,distr_st,
                      num_cells,Rbmin,Rbmax,wave,front,distr_site);
            }
            free(insert_site);
            free(distr_site); 
        }

        bubbles.bubble = NO;
        debug_print("bubbles_mpi","Left conserv_make_bubbles_orig_ver2()\n"); 
        return num_bubble;
}
#endif /* if defined(USE_OVERTURE) */

LOCAL int neighbor_mass_deficiency_on_curve(
      INSERT_SITE  insert_site,
      COMPONENT    l_comp,
      Locstate     distr_st,
      int          n_cells,
      int          *Rbmin,
      int          *Rbmax,
      Wave         *wave,
      Front        *front,
      int          **distr_site)
{
      int          i, j;
      int          ic[MAXD], imin[MAXD], imax[MAXD];
      double        fmin[MAXD], fmax[MAXD];
      double        *crds;
      Locstate     st;
      double        temp, pres, P_sat, vapor_press, p_jump, D_add, c_sqr;
      double        vel0, vel1, in_vel0, in_vel1, n_vel0, n_vel1;
      CURVE        **cc;
      BOND         *b;
      POINT        *pt;
      RECT_GRID    *gr = wave->rect_grid;
      int          stype;
      double        eps;

      eps = MACH_EPS*100;
      for(i = 0; i < gr->dim; i++)
      {
          imin[i] = max(0,Rbmin[i]);
          imax[i] = min(gr->gmax[i], Rbmax[i]);
      }
      for(i = 0; i < gr->dim; i++)
      {
          fmin[i] = cell_center(imin[i],i,gr) - 0.5*gr->h[i];
          fmax[i] = cell_center(imax[i],i,gr) + 0.5*gr->h[i];
      }

      vapor_press = pressure(bubbles.st);
      in_vel0 = vel(0,bubbles.st);
      in_vel1 = vel(1,bubbles.st);

      for(cc = front->interf->curves; cc && *cc; cc++)
      {
          if(wave_type(*cc) != CONTACT)
              continue;
          b = (*cc)->first;
          while(b)
          {
              pt = b->start;
              crds = Coords(pt);
              if(crds[0] >= fmin[0]-eps && crds[0] <= fmax[0]+eps &&
                 crds[1] >= fmin[1]-eps && crds[1] <= fmax[1]+eps)
              {
                  if(negative_component(*cc) == l_comp)
                  {
                      if (pt == (*cc)->start->posn && b == (*cc)->first)
                          st = left_start_state(*cc);
                      else
                          st = left_state(pt);
                  }
                  if(positive_component(*cc) == l_comp)
                  {
                      if (pt == (*cc)->start->posn && b == (*cc)->first)
                          st = right_start_state(*cc);
                      else
                          st = right_state(pt);
                  }
                  pres = pressure(st);
                  c_sqr = sound_speed_squared(st);
#if defined(USE_ADIABAT_PRES)
                  if(pres < vapor_press)
                  {
                      double dens_incr; 
                      if(debugging("conserv_bubble"))
                      {
                          /* print_general_vector("intfc point",crds,2,"\n"); */
                          /* verbose_print_state("before_chg_intfc",st); */
                      }
                      stype = state_type(st);
                      set_state(st,TGAS_STATE,st);
                      Press(st) = vapor_press;
                      p_jump = vapor_press - pres;

                      /* adiabatically add density */
                      D_add = p_jump/c_sqr;
                      dens_incr = D_add/Dens(st)*100;
                      Dens(st) += D_add;
                      set_state(st,stype,st);
                 
                      if(debugging("conserv_bubble"))
                      {
                          /* verbose_print_state("after_chg_intfc",st); */
                          /* printf("density incre %g\n", dens_incr); */
                      }
                  }
#endif /* if defined(USE_ADIABAT_PRES) */
              }
              b = b->next;
          };
          pt = (*cc)->end->posn;
          crds = Coords(pt);
          if(crds[0] >= fmin[0]-eps && crds[0] <= fmax[0]+eps &&
             crds[1] >= fmin[1]-eps && crds[1] <= fmax[1]+eps)
          {
              if(negative_component(*cc) == l_comp)
                  st = left_end_state(*cc);
              if(positive_component(*cc) == l_comp)
                  st = right_end_state(*cc);
              pres = pressure(st);
              c_sqr = sound_speed_squared(st);
#if defined(USE_ADIABAT_PRES)
              /* if(pres < vapor_press) */
              {
                  if(debugging("conserv_bubble"))
                  {
                  /*     print_general_vector("intfc point",crds,2,"\n"); */
                  /*     verbose_print_state("before_chg_intfc",st); */
                  }
                  stype = state_type(st);
                  set_state(st,TGAS_STATE,st);
                  Press(st) = vapor_press;
                  p_jump = vapor_press - pres;
                                                                                                                  
                  /* adiabatically add density */
                  D_add = p_jump/c_sqr;
                  Dens(st) += D_add;
                  set_state(st,stype,st);
                                                                                                               
                  /* if(debugging("conserv_bubble")) */
                  /*     verbose_print_state("after_chg_intfc",st); */
              }
#endif /* if defined(USE_ADIABAT_PRES) */
          }
      }
      return YES;
}


LOCAL int neighbor_mass_deficiency(
      INSERT_SITE  insert_site,
      COMPONENT    l_comp,
      int          *Rbmin,
      int          *Rbmax,
      Wave         *wave,
      Front        *front,
      int          **distr_site,
      double        *D_add)
{
      int          x_size, y_size, i, j, side_y, side_x;
      int          ic[MAXD], icent[MAXD], count = 0;
      double        *crds, area, alpha = 1.0, d_dif, s_alpha = 0.0;
      Locstate     st;
      double        temp, pres, P_sat, vapor_press, p_jump, c_sqr;

      x_size = ((Rbmax[0] - Rbmin[0]) - (Rbmax[0] - Rbmin[0])%2)/2 + 1;
      y_size = ((Rbmax[1] - Rbmin[1]) - (Rbmax[1] - Rbmin[1])%2)/2 + 1;
      icent[0] = ((Rbmax[0] + Rbmin[0]) + (Rbmax[0] + Rbmin[0])%2)/2;
      icent[1] = ((Rbmax[1] + Rbmin[1]) + (Rbmax[1] + Rbmin[1])%2)/2;
      area = wave->rect_grid->h[0]*wave->rect_grid->h[1];
      vapor_press = insert_site.vapor_press;
      *D_add = 0.0;

      /*
      if(debugging("conserv_bubble"))
      {
          printf("insert_site[%d %d], Rbmin[%d %d], Rbmax[%d %d] on Domain [%d %d]X[%d %d]\n",
            insert_site.icoords[0], insert_site.icoords[1],
            Rbmin[0], Rbmin[1], Rbmax[0], Rbmax[1],
            -(wave->rect_grid->lbuf[0]), wave->rect_grid->gmax[0]+wave->rect_grid->ubuf[0],
            -(wave->rect_grid->lbuf[1]), wave->rect_grid->gmax[1]+wave->rect_grid->ubuf[1]);
          printf("half diameter: x_size %d, y_size %d. Icenter[%d %d]\n",
             x_size, y_size, icent[0], icent[1]);
      }
      */

      for(side_y = -1; side_y <=1; side_y+=2)
      {
          int starty;
          if(side_y == -1) starty = 0;
          else
              starty = 1;
          for(j = starty; j <= y_size; j++)
          {
              int startx;
              ic[1] = icent[1] + side_y*j;
              for(side_x = -1; side_x <=1; side_x+=2)
              {
                  if(side_x == -1) startx = 0;
                  else
                     startx = 1;
                  for(i = startx; i <= x_size; i++)
                  {
                      ic[0] = icent[0] + side_x*i;
                      if(ic[0] < Rbmin[0] || ic[0] > Rbmax[0] ||
                         ic[1] < Rbmin[1] || ic[1] > Rbmax[1])
                          continue;
                      if(ic[0] >= wave->rect_grid->gmax[0] ||
                         ic[1] >= wave->rect_grid->gmax[1] ||
                         ic[0] < 0 || ic[1] < 0)
                          continue;
                      if(Rect_comp(ic,wave) != l_comp)
                          continue;
                      if(distr_site[ic[0]][ic[1]] == YES)
                          continue;

                      crds = Rect_coords(ic,wave);
                      if(crds[0] < 0.0)
                          continue;
                      st = Rect_state(ic,wave);
                      temp = temperature(st);
                      pres = pressure(st);
                      c_sqr = sound_speed_squared(st);
                      P_sat = saturated_pres(st,temp);

                      if (wave->rect_grid->Remap.remap == CYLINDRICAL_REMAP &&
				is_rotational_symmetry())
                          alpha = crds[0];

#if defined(USE_ADIABAT_PRES)                      
                      if(pres < vapor_press)
                      {
                          p_jump = vapor_press - pres;
                          (*D_add) += p_jump/c_sqr*alpha*area;
                          count++;
                      }
#else /* if defined(USE_ADIABAT_PRES) */
                      if(pres < P_sat)
                          count++;
#endif /* if defined(USE_ADIABAT_PRES) */
                  }
              }
          }
      }

      /*
      if(debugging("conserv_bubble"))
      {
          printf("%d cells have mass deficiency (reg[%d])\n",
                count, (Rbmax[0]-Rbmin[0]+1)*(Rbmax[1]-Rbmin[1]+1));
      }
      */
      return count;
}

LOCAL void distr_to_neighbors(
      INSERT_SITE  insert_site,
      COMPONENT    l_comp,
      Locstate     distr_st,
      int          n_cells,
      int          *Rbmin,
      int          *Rbmax,
      Wave         *wave,
      Front        *front,
      int          **distr_site)
{
      int          x_size, y_size, i, j, side_y, side_x;
      int          ic[MAXD], icent[MAXD];
      Locstate     st;
      double        D, T, T_sat, P_sat, pres, d_sat, c_sqr, vapor_press, p_jump;
      double        *crds, area, alpha = 1.0, d_dif, s_alpha = 0.0;
      double        avail, d_incr, m0_incr, m1_incr, e_incr;
      double        in_vel0, in_vel1, n_vel0, n_vel1, vel0, vel1, D_add;
      int          num_cell = 0, count = 0, all_cell, more_count = 0;
      static Locstate temp_st = NULL, sum_st, in_bub_st, v_st;
      double        total_add_mass = 0.0;

      if(temp_st == NULL)
      {
          alloc_state(front->interf,&temp_st,front->sizest);
          alloc_state(front->interf,&sum_st,front->sizest);
          alloc_state(front->interf,&in_bub_st,front->sizest);
          alloc_state(front->interf,&v_st,front->sizest);
      }

      if (wave->rect_grid->Remap.remap == CYLINDRICAL_REMAP && is_rotational_symmetry())
          alpha = insert_site.center[0];

      copy_state(v_st,bubbles.st);
      copy_state(in_bub_st,distr_st);
      area = alpha*pi*sqr(bubbles.r);
      vapor_press = pressure(bubbles.st);

      Dens(in_bub_st) -= Dens(v_st)*area;
      Mom(in_bub_st)[0] -= Mom(v_st)[0]*area;
      Mom(in_bub_st)[1] -= Mom(v_st)[1]*area;
      Energy(in_bub_st) -= Energy(v_st)*area;

      in_vel0 = Mom(in_bub_st)[0]/Dens(in_bub_st);
      in_vel1 = Mom(in_bub_st)[1]/Dens(in_bub_st);

      if(debugging("conserv_bubble"))
      {
          printf("Before distribution " 
               "Taken-out: Den %g, Mom[0] %g, Mom[1] %g, Eng %g\n",
                Dens(in_bub_st), Mom(in_bub_st)[0],
               Mom(in_bub_st)[1], Energy(in_bub_st));
      }

      Dens(sum_st) = 0.0;
      Mom(sum_st)[0] = 0.0;
      Mom(sum_st)[1] = 0.0;
      Energy(sum_st) = 0.0;
      set_type_of_state(sum_st,GAS_STATE);
      Set_params(sum_st,distr_st);

      x_size = ((Rbmax[0] - Rbmin[0]) - (Rbmax[0] - Rbmin[0])%2)/2 + 1;
      y_size = ((Rbmax[1] - Rbmin[1]) - (Rbmax[1] - Rbmin[1])%2)/2 + 1;
      icent[0] = ((Rbmax[0] + Rbmin[0]) + (Rbmax[0] + Rbmin[0])%2)/2;
      icent[1] = ((Rbmax[1] + Rbmin[1]) + (Rbmax[1] + Rbmin[1])%2)/2;
      area = wave->rect_grid->h[0]*wave->rect_grid->h[1];

      for(side_y = -1; side_y <=1; side_y+=2)
      { 
          int starty;
          if(side_y == -1) starty = 0; 
          else
              starty = 1;   
          for(j = starty; j <= y_size; j++)
          {
              int startx;
              ic[1] = icent[1] + side_y*j;  
              for(side_x = -1; side_x <=1; side_x+=2)
              {
                  if(side_x == -1) startx = 0; 
                  else
                     startx = 1;
                  for(i = startx; i <= x_size; i++)
                  {
                      ic[0] = icent[0] + side_x*i;  
                      if(ic[0] < Rbmin[0] || ic[0] > Rbmax[0] ||
                         ic[1] < Rbmin[1] || ic[1] > Rbmax[1])
                          continue;
                      if(ic[0] >= wave->rect_grid->gmax[0] ||
                         ic[1] >= wave->rect_grid->gmax[1] ||
                         ic[0] < 0 || ic[1] < 0)
                          continue;
                      if(Rect_comp(ic,wave) != l_comp)
                          continue;

                      crds = Rect_coords(ic,wave);

                      if(crds[0] < 0.0)
                          continue;
                      more_count++;
                     
                      st = Rect_state(ic,wave);
                      T = temperature(st);
                      pres = pressure(st);
                      P_sat = saturated_pres(st,T);
#if defined(USE_ADIABAT_PRES)
                      if(pres > vapor_press)
                          continue;
#else /* if defined(USE_ADIABAT_PRES) */
                      if(pres > P_sat)
                          continue;
#endif /* if defined(USE_ADIABAT_PRES) */

                      if (wave->rect_grid->Remap.remap == CYLINDRICAL_REMAP &&
				is_rotational_symmetry())
                      {
                          alpha = crds[0];
                          s_alpha += alpha;
                      }
                      
                      Dens(sum_st) += alpha*area*Dens(st);
                      Mom(sum_st)[0] += alpha*area*Mom(st)[0];
                      Mom(sum_st)[1] += alpha*area*Mom(st)[1];
                      Energy(sum_st) += alpha*area*Energy(st);
                      count++;
                  }
              }
          }
      }
 
      if(count == 0)
          return;
 
      if (wave->rect_grid->Remap.remap == CYLINDRICAL_REMAP &&
		is_rotational_symmetry())
          s_alpha /= count;
      else
          s_alpha = 1.0;

      if(debugging("conserv_bubble"))
      {
          printf("\nBefore add liquid states, all quantity in region:\n"
            "Den %g, Mom[0] %g, Mom[1] %g, Eng %g\n",
                Dens(sum_st), Mom(sum_st)[0],
               Mom(sum_st)[1], Energy(sum_st));

          printf("Before add liquid states, averaged states in region:");  
          Dens(sum_st) /= (s_alpha*count*area);
          Mom(sum_st)[0] /= (s_alpha*count*area);
          Mom(sum_st)[1] /= (s_alpha*count*area);
          Energy(sum_st) /= (s_alpha*count*area);
          verbose_print_state("before_add_liquid, avg st",sum_st);
          Dens(sum_st) *= s_alpha*count*area;
          Mom(sum_st)[0] *= s_alpha*count*area;
          Mom(sum_st)[1] *= s_alpha*count*area;
          Energy(sum_st) *= s_alpha*count*area;
      }

      all_cell = (Rbmax[0]-Rbmin[0]+1)*(Rbmax[1]-Rbmin[1]+1)-n_cells;
      avail = (1.0*count)/all_cell;
      d_incr = Dens(in_bub_st)*avail/Dens(sum_st);
      e_incr = Energy(in_bub_st)*avail/Energy(sum_st);
      m0_incr = Mom(in_bub_st)[0]*avail/Mom(sum_st)[0];
      m1_incr = Mom(in_bub_st)[1]*avail/Mom(sum_st)[1];

      Dens(sum_st) += Dens(in_bub_st)*avail;
      Mom(sum_st)[0] += Mom(in_bub_st)[0]*avail;
      Mom(sum_st)[1] += Mom(in_bub_st)[1]*avail;
      Energy(sum_st) += Energy(in_bub_st)*avail;

      Dens(sum_st) /= (s_alpha*count*area);
      Mom(sum_st)[0] /= (s_alpha*count*area);
      Mom(sum_st)[1] /= (s_alpha*count*area);
      Energy(sum_st) /= (s_alpha*count*area);

      /*
      Dens(in_bub_st) *= (1.0-avail);
      Mom(in_bub_st)[0] *= (1.0-avail);
      Mom(in_bub_st)[1] *= (1.0-avail);
      Energy(in_bub_st) *= (1.0-avail);
      */

      if(debugging("conserv_bubble"))
      {
            printf("Distribution to %d cells (more_cells count %d) (reg[%d]) Density increament %g%\n",
                count, more_count, (Rbmax[0]-Rbmin[0]+1)*(Rbmax[1]-Rbmin[1]+1), d_incr*100);
          printf("After add liquid states, averaged states in [%g%] region:", avail*100);  
          verbose_print_state("after_add_liquid, avg st",sum_st);
      }
      
      for(side_y = -1; side_y <=1; side_y+=2)
      { 
          int starty;
          if(side_y == -1) starty = 0; 
          else
              starty = 1;   
          for(j = starty; j <= y_size; j++)
          {
              int startx;
              ic[1] = icent[1] + side_y*j;  
              for(side_x = -1; side_x <=1; side_x+=2)
              {
                  if(side_x == -1) startx = 0; 
                  else
                     startx = 1;
                  for(i = startx; i <= x_size; i++)
                  {
                      ic[0] = icent[0] + side_x*i;  
                      if(ic[0] < Rbmin[0] || ic[0] > Rbmax[0] ||
                         ic[1] < Rbmin[1] || ic[1] > Rbmax[1])
                          continue;
                      if(ic[0] >= wave->rect_grid->gmax[0] ||
                         ic[1] >= wave->rect_grid->gmax[1] ||
                         ic[0] < 0 || ic[1] < 0)
                          continue;
                      if(Rect_comp(ic,wave) != l_comp)
                          continue;
                      /* Each cell is distributed once */
                      if(distr_site[ic[0]][ic[1]] == YES)
                          continue;

                      crds = Rect_coords(ic,wave);
                      st = Rect_state(ic,wave);

                      T = temperature(st);
                      pres = pressure(st);
                      P_sat = saturated_pres(st,T);

#if defined(USE_ADIABAT_PRES)
                      if(pres > vapor_press)
                          continue;
#else /* if defined(USE_ADIABAT_PRES) */
                      if(pres > P_sat)
                          continue;
#endif /* if defined(USE_ADIABAT_PRES) */

                      copy_state(temp_st,st);
                      if (wave->rect_grid->Remap.remap == CYLINDRICAL_REMAP &&
				is_rotational_symmetry())
                          alpha = crds[0];

                      /* if(debugging("conserv_bubble") && */
		      /*    fabs(0.00144625-crds[0]) < 0.0001 && */
	              /*  fabs(0.168544-crds[1]) < 0.0001) */
                      if(debugging("conserv_bubble"))
                      {
                          /***
                          printf("cell [%d %d] comp[%d] before add states\n", 
                              ic[0], ic[1], Rect_comp(ic,wave));
                          verbose_print_state("before_add", st);
                          ***/
                      }
                      D = Dens(st);
                      T = temperature(st);
                      c_sqr = sound_speed_squared(st);
                      vel0 = vel(0,st);
                      vel1 = vel(1,st);

                      distr_site[ic[0]][ic[1]] = YES;

#if defined(USE_ADIABAT_PRES)
                      set_state(st,TGAS_STATE,st);
                      pres = Press(st);
                      Press(st) = vapor_press;
                      p_jump = vapor_press - pres;

                      /* adiabatically add density */
                      D_add = p_jump/c_sqr;
                      Dens(st) += D_add;

                      /* Mom increment */
                      /*
                      n_vel0 = (vel0*D + D_add*in_vel0)/Dens(st);
                      n_vel1 = (vel1*D + D_add*in_vel1)/Dens(st);
                      Vel(st)[0] = n_vel0;
                      Vel(st)[1] = n_vel1;
                      */
                      set_state(st,GAS_STATE,st);
                      total_add_mass +=  D_add*alpha*area;
                      num_cell++;
#endif /* if defined(USE_ADIABAT_PRES) */

                      /* if(debugging("conserv_bubble") && */
		      /*    fabs(0.00144625-crds[0]) < 0.0001 && */
	              /* fabs(0.168544-crds[1]) < 0.0001) */
                      if(debugging("conserv_bubble"))
                      { 
                          int k = 0, MAX_NUM_ITERATIONS = 60;
                          double dd;
                          /* printf("cell [%d %d] after add states\n", ic[0], ic[1]); */
                          /* verbose_print_state("after_add", st); */

                          /**
                          printf("Add by adiabatic:\n");
                          dd =  Dens(temp_st)*d_incr/MAX_NUM_ITERATIONS;
                          set_state(temp_st,TGAS_STATE,temp_st);
                          do{
                              Press(temp_st) += c_sqr*dd;
                              state_on_adiabat_with_pr(temp_st,Press(temp_st),temp_st,TGAS_STATE);
                              c_sqr = sound_speed_squared(temp_st);
                              Dens(temp_st) += dd;
                              k++;
                          }while(k < MAX_NUM_ITERATIONS);

                          /* Press(temp_st) += c_sqr*d_incr*Dens(temp_st); */
                          /* Dens(temp_st) *= (1.0+d_incr); */
                          set_state(temp_st,GAS_STATE,temp_st);
                          verbose_print_state("evenly distribute", temp_st);
                          **/
                      }
                  }
              }
          }
      }

      if(debugging("conserv_bubble"))
      {
          /*
          printf("After distribution to %d cells (count %d, reg[%d]),\n"
                "Remaining: Den %g(%g%), Mom[0] %g, Mom[1] %g, Eng %g\n",
                num_cell, count, (Rbmax[0]-Rbmin[0]+1)*(Rbmax[1]-Rbmin[1]+1),
                 Dens(in_bub_st)-total_add_mass, total_add_mass/Dens(in_bub_st)*100,
                 Mom(in_bub_st)[0],
                Mom(in_bub_st)[1], Energy(in_bub_st));
          */
          printf("After distribution to %d cells (count %d, reg[%d]),\n"
              "Remaining: Den %g; Distributed: %g%\n",
             num_cell, count, (Rbmax[0]-Rbmin[0]+1)*(Rbmax[1]-Rbmin[1]+1),
              Dens(in_bub_st)-total_add_mass, total_add_mass/Dens(in_bub_st)*100);
      }
}

LOCAL int intfc_on_rect_domain_2d(
      INTERFACE    *interf, 
      double        *cl,
      double        *cu)
{
      /* Interface cross with left rectangle boundary */
      if(YES == intfc_intersect_segment(interf, 0, cl[0], cl[1], cu[1]))
          return YES;
      /* Interface cross with right rectangle boundary */
      if(YES == intfc_intersect_segment(interf, 0, cu[0], cl[1], cu[1]))
          return YES;
      /* Interface cross with lower rectangle boundary */
      if(YES == intfc_intersect_segment(interf, 1, cl[1], cl[0], cu[0]))
          return YES;
      /* Interface cross with upper rectangle boundary */
      if(YES == intfc_intersect_segment(interf, 1, cu[1], cl[0], cu[0]))
          return YES;
      if(YES == intfc_cross_rect(interf,cl, cu))
          return YES;
      return NO;
}	/* end of intfc_on_rect_domain_2d */

LOCAL void blks_bubble_ocupy(
      INSERT_SITE  insert_site,
      RECT_GRID    *gr,
      int          *Rbmin,
      int          *Rbmax,
      double        spacing)
{
      int          i, dim = gr->dim;
      int          size, is_size_odd;

      size = bubbles.size;
      is_size_odd = (size%2 == 0 ? NO : YES);

      if(is_size_odd)
      {
          for(i = 0; i < dim; i++)
          {
              Rbmin[i] = insert_site.icoords[i] - ((size-1)/2 + (int)(spacing*size));
              Rbmax[i] = insert_site.icoords[i] + ((size-1)/2 + (int)(spacing*size));
          }
      }
      else
      {
          for(i = 0; i < dim; i++)
          {
              Rbmin[i] = insert_site.icoords[i] - ((size/2-1) + (int)(spacing*size));
              Rbmax[i] = insert_site.icoords[i] + (size/2 + (int)(spacing*size));
          }
      }
}

LOCAL void create_bubbles_ver2(
      INSERT_SITE  *insert_site,
      COMPONENT    l_comp,
      Wave         *wave,
      Front        *front,
      Locstate     rst,
      int          *imin,
      int          *imax)
{
      double        r,theta,d_theta,step;
      int          n_theta, ic[MAXD], ix, iy, k;
      double        *c, coords[MAXD];
      Locstate     state;
      RECT_GRID    *gr = wave->rect_grid; 
      NODE         *ns, *ne;
      INTERFACE    *intfc = front->interf;
      CURVE        *cv;
      BOND         *b;
      INTERFACE    *sv_intfc;

      debug_print("bubbles_mpi","Entered create_bubbles()\n");

      sv_intfc = current_interface();
      set_current_interface(intfc);
      c = insert_site->center;
      r = bubbles.r;

      /* make node and curve */
      coords[0] = c[0] + r;
      coords[1] = c[1];
      ns = make_node(Point(coords));
      ne = ns;
      node_type(ns) = CLOSED_NODE;
      cv = make_curve(l_comp,bubbles.bubblecomp,ns,ne);

      n_theta = 4*(irint(pi*r/1.5/min(gr->h[0], gr->h[1])) + 1);
      d_theta = 2*pi/n_theta;

      if (debugging("bubbles_mpi"))
          printf("after make_curve: n_theta = %d, n_theta/4 = %d, d_theta = %g\n",
                n_theta,irint(n_theta/4), d_theta);
      /* insert points into the curve */
      theta = -d_theta;
      for (k = 1; k < n_theta; k++)
      {
          coords[0] = c[0] + r*cos(theta);
          coords[1] = c[1] + r*sin(theta);
          Check_return(insert_point_in_bond(Point(coords),cv->last,cv),
                       create_bubbles);
          if (debugging("bubbles_mpi"))
          {
              printf("start = (%g, %g), point = (%g, %g),\n",
                     Coords(cv->last->start)[0],
                     Coords(cv->last->start)[1],
                     coords[0], coords[1]);
          }
          theta -= d_theta;
      }

      wave_type(cv) = CONTACT;
      start_status(cv) = end_status(cv) = INCIDENT;
      create_time(cv) = front->time;

      /*init front_states*/

      if (debugging("bubbles_mpi"))
      {
          (void) printf("Initializing states on curve %llu\n",
                        curve_number(cv));
          print_curve(cv);
      }

      front_init_2d(cv->first->start,
                   Hyper_surf(cv),
                   left_start_state(cv),right_start_state(cv),
                   front,wave,rst);
      front_init_2d(cv->first->start,
                   Hyper_surf(cv),
                   left_state(cv->first->start),right_state(cv->first->start),
                   front,wave,rst);
      for (b = cv->first; b != cv->last; b = b->next)
          front_init_2d(b->end,
                       Hyper_surf(cv),
                       left_state(b->end),right_state(b->end),
                       front,wave,rst);
      front_init_2d(cv->last->end,
                   Hyper_surf(cv),left_state(cv->last->end),
                   right_state(cv->last->end),front,wave,rst);
      front_init_2d(cv->last->end,
                   Hyper_surf(cv),left_end_state(cv),
                   right_end_state(cv),front,wave,rst);

      if (debugging("bubbles_mpi"))
      {
          BOND *b;
          (void) printf("After initializing states on curve %llu\n",
                        curve_number(cv));
          print_curve(cv);
          /**
          b = cv->first;
          printf("%-11g %-11g  ",
                      Coords(b->start)[0],Coords(b->start)[1]);
          verbose_print_state("left", left_state(b->start));
          verbose_print_state("right", left_state(b->start));
          while(b!=NULL)
          {
              printf("%-11g %-11g  ",
                      Coords(b->end)[0],Coords(b->end)[1]);
              verbose_print_state("left", left_state(b->end));
              verbose_print_state("right", left_state(b->end));
              b = b->next;
          }
          **/
          /* show_curve_states(cv); */
      }

      /*
      printf("create_bubble at cent[%g %g][%d %d] on patch[%d][%d] ref_T %g, step %d, time %12.10f\n",
		insert_site->center[0], insert_site->center[1],
	        insert_site->icoords[0], insert_site->icoords[1], 
	         front->patch_number, front->patch_level, temperature(rst),
		 front->step,create_time(cv));
      verbose_print_state("bubble init st",bubbles.st);
      printf("Initializing states on curve %llu, front->interf %p, cv->interf %p\n",
                        curve_number(cv), front->interf, cv->interface);
      print_curve(cv);
      */

      /*init interior states for bubble*/
      k = 0;
      for ( iy = imin[1]; iy <= imax[1]; iy++)
      {
          ic[1] = iy;
          for ( ix = imin[0]; ix <= imax[0]; ix++)
          {
              int     stype;
              double   dist;

              ic[0] = ix;
              coords[0] = Rect_coords(ic,wave)[0];
              coords[1] = Rect_coords(ic,wave)[1];
              if( (dist = (sqr(c[0]-coords[0]) + sqr(c[1]-coords[1]))) > sqr(r))
                  continue;
              state = Rect_state(ic,wave);
              if (debugging("bubbles_mpi"))
              {
                  (void) printf("Original, comp = %d, bub_comp %d ",
                          Rect_comp(ic,wave),bubbles.bubblecomp);
                  printf("coords = (%g, %g), ",coords[0],coords[1]);
                  printf("distance = %g <= %g\n", sqrt(dist), bubbles.r);
                  verbose_print_state("state",state);
              }
              stype = state_type(state);
              Rect_comp(ic,wave) = bubbles.bubblecomp;
              copy_state(state,bubbles.st);
              set_state(state,stype,state);
              k++;
          }
      }
      if (debugging("bubbles_mpi"))
          printf("End of initilization of bubbles in %d cells\n", k);

      set_current_interface(sv_intfc);
      debug_print("bubbles_mpi","Left create_bubbles()\n");
}	/* end of create_bubbles_ver2 */

LOCAL void set_state_on_phase_boundary_ver2(
      Locstate   lst, 
      Locstate   CC_lst,
      double      ref_temp)
{
      double      v[MAXD], press;
      int        i;
      int        stype;

      for(i = 0; i < 2; i++)
          v[i] = vel(i,lst);

      stype = state_type(bubbles.st);
      set_state(bubbles.st,TGAS_STATE,bubbles.st);

      /* Set vapor bubble state according to Clausius-Clapeyron */
      /* Explicit EOS model is used, should be removed */
      Dens(bubbles.st) = saturated_pres(lst,ref_temp)
                          /bubbles.R/ref_temp;
      press = Press(bubbles.st) = saturated_pres(lst,ref_temp);
      for( i=0; i<2; i++) 
          Vel(bubbles.st)[i] = v[i];
      set_state(bubbles.st,stype,bubbles.st);

      /* Set liquid state according to Clausius-Clapeyron */
      set_state(CC_lst,TGAS_STATE,lst);
      Dens(CC_lst) = density_via_PT(lst, press, ref_temp);
      Press(CC_lst) = press;
      set_state(CC_lst,GAS_STATE,CC_lst);
}	/* end of set_state_on_phase_boundary_ver2 */

LOCAL void set_state_on_phase_boundary_ver3(
      Locstate   lst,
      Locstate   CC_lst,
      double      ref_temp)
{
      double      v[MAXD], press;
      int        i;
      int        stype;

      for(i = 0; i < 3; i++)
          v[i] = vel(i,lst);

      stype = state_type(bubbles.st);
      set_state(bubbles.st,TGAS_STATE,bubbles.st);

      /* Set vapor bubble state according to Clausius-Clapeyron */
      /* Explicit EOS model is used, should be removed */
      Dens(bubbles.st) = saturated_pres(lst,ref_temp)
                          /bubbles.R/ref_temp;
      press = Press(bubbles.st) = saturated_pres(lst,ref_temp);
      for( i=0; i<3; i++)
          Vel(bubbles.st)[i] = v[i];
      set_state(bubbles.st,stype,bubbles.st);

      /* Set liquid state according to Clausius-Clapeyron */
      set_state(CC_lst,TGAS_STATE,lst);
      Dens(CC_lst) = density_via_PT(lst, press, ref_temp);
      Press(CC_lst) = press;
      set_state(CC_lst,GAS_STATE,CC_lst);
}

LOCAL void avg_states_in_reg_2d(
      INSERT_SITE  *insert_site,
      int    *imin,
      int    *imax,
      Wave   *wave,
      Front  *front,
      Locstate  ref_st,
      int       *num_cells)
{
      int    iy, ix, ic[MAXD];
      double  *crds, *center; 
      Locstate  state;
      double     rho, m[MAXD], E;
      int    first = YES;
      RECT_GRID  *gr = wave->rect_grid;
      double      alpha = 1.0; 

      Dens(ref_st) = 0.0;  
      Mom(ref_st)[0] = 0.0;
      Mom(ref_st)[1] = 0.0;
      Energy(ref_st) = 0.0;
 
      center = insert_site->center;
      rho = m[0] = m[1] = E = 0;
      *num_cells = 0;
      for(iy = imin[1]; iy <= imax[1]; iy++)
      {
          for(ix = imin[0]; ix<= imax[0]; ix++)
          {
              ic[0] = ix; ic[1] = iy;
              crds = Rect_coords(ic, wave);
              if((sqr(center[0]-crds[0]) + sqr(center[1]-crds[1])) >
                 sqr(bubbles.r))
                  continue;

              state = Rect_state(ic,wave);
              if(YES == first)
              {
                  Set_params(ref_st,state);
                  first = NO;
              }
              if(Different_params(ref_st,state))
              {
                  printf("ERROR: avg_states_in_reg, diff params\n");
                  clean_up(ERROR);
              }
	      if (is_rotational_symmetry())
	      {
                  if (gr->Remap.remap == CYLINDRICAL_REMAP)
                      alpha = crds[0];
	          alpha = 1.0;
	      }
              Dens(ref_st) += alpha*Dens(state);
              Mom(ref_st)[0] += alpha*Mom(state)[0];
              Mom(ref_st)[1] += alpha*Mom(state)[1];
              Energy(ref_st) += alpha*Energy(state);
              (*num_cells)++;
          }
      }

      if (is_rotational_symmetry())
      {
          if (gr->Remap.remap == CYLINDRICAL_REMAP)
              alpha = center[0];
          alpha = 1.0;
      }

      Dens(ref_st) /= ((*num_cells)*alpha);  
      Mom(ref_st)[0] /= ((*num_cells)*alpha);
      Mom(ref_st)[1] /= ((*num_cells)*alpha);
      Energy(ref_st) /= ((*num_cells)*alpha);
}

LOCAL int press_ascend(
        const void      *c1,
        const void      *c2)
{
        if ( ((INSERT_SITE*)c1)->press >
             ((INSERT_SITE*)c2)->press ) return 1;
        else return -1;
}               /*end press_ascend*/

LOCAL int  check_in_bubble_diameter_2d(
      INSERT_SITE  *insert_site, 
      Wave         *wave, 
      Front        *front,
      COMPONENT    l_comp,
      int          *dimin, 
      int          *dimax,
      int          *inmin,
      int          *inmax)
{
      int          i, ix, iy, size, is_size_odd; 
      int          ic[MAXD];
      int          dim = wave->rect_grid->dim;
      double        center[MAXD], *crds;
      Locstate     state;
      double        T_sat, P_sat, press;

      size = bubbles.size;
      is_size_odd = (size%2 == 0 ? NO : YES); 
        
      /* Find the cells the bubble will occupy. */
      /* If the bubble diameter is odd number of mesh sizes,
       * bubble center is at the center of cell "insert_site.icoords". Otherwise,
       * It is centered at the upper corner of cell "insert_site.icoords". 
       */
      if(is_size_odd)
      {
          for(i = 0; i < dim; i++)
          {
              inmin[i] = insert_site->icoords[i] - (size-1)/2; 
              inmax[i] = insert_site->icoords[i] + (size-1)/2; 
              center[i] = Rect_coords(insert_site->icoords, wave)[i];
          }
      }
      else
      {
          for(i = 0; i < dim; i++)
          {
              inmin[i] = insert_site->icoords[i] - (size/2 - 1); 
              inmax[i] = insert_site->icoords[i] + size/2; 
              center[i] = Rect_coords(insert_site->icoords, wave)[i] + 
                              0.5*wave->rect_grid->h[i];
          }
      }    
    
      for(iy = inmin[1]; iy <= inmax[1]; iy++)
      {
          for(ix = inmin[0]; ix<= inmax[0]; ix++)
          {
              ic[0] = ix; ic[1] = iy;
              if(ic[0] < dimin[0] || ic[0] >= dimax[0] ||
                 ic[1] < dimin[1] || ic[1] >= dimax[1])
                  return NO; 
              crds = Rect_coords(ic, wave);
              if((sqr(center[0]-crds[0]) + sqr(center[1]-crds[1])) > 
                 sqr(bubbles.r))
                  continue;

              if (Rect_comp(ic,wave) != l_comp)
                  return NO;

              /* Physics checking */
              state = Rect_state(ic,wave);
              T_sat = temperature(state);

              /* THIS IS NOT A GOOF FIX */
              if(T_sat < 294.261)
                  T_sat = 294.262;

              P_sat = saturated_pres(state,T_sat);
              press = pressure(state);
              
              /* From homo. model */
              /* if(press > -5.0) for simulation sp_mp_jet40_new2.r */
              /*     return NO; */
              /* Clausis-Clapeyron */
              if (press + 5 > P_sat)
                  return NO;
          }
      }
      insert_site->center[0] = center[0];
      insert_site->center[1] = center[1];
      return YES;
}	/* check_in_bubble_diameter_2d */

LOCAL int  check_in_bubble_diameter_3d(
      INSERT_SITE  *insert_site,
      Wave         *wave,
      Front        *front,
      COMPONENT    l_comp,
      int          *dimin,
      int          *dimax,
      int          *inmin,
      int          *inmax)
{
      int          i, ix, iy, iz, size, is_size_odd;
      int          ic[MAXD];
      int          dim = wave->rect_grid->dim;
      double        center[MAXD], *crds;
      Locstate     state;
      double        T_sat, P_sat, press;

      size = bubbles.size;
      is_size_odd = (size%2 == 0 ? NO : YES);

      /* Find the cells the bubble will occupy. */
      /* If the bubble diameter is odd number of mesh sizes,
       * bubble center is at the center of cell "insert_site.icoords". Otherwise,
       * It is centered at the upper corner of cell "insert_site.icoords".
       */
      if(is_size_odd)
      {
          for(i = 0; i < dim; i++)
          {
              inmin[i] = insert_site->icoords[i] - (size-1)/2;
              inmax[i] = insert_site->icoords[i] + (size-1)/2;
              center[i] = Rect_coords(insert_site->icoords, wave)[i];
          }
      }
      else
      {
          for(i = 0; i < dim; i++)
          {
              inmin[i] = insert_site->icoords[i] - (size/2 - 1);
              inmax[i] = insert_site->icoords[i] + size/2;
              center[i] = Rect_coords(insert_site->icoords, wave)[i] +
                              0.5*wave->rect_grid->h[i];
          }
      }

      for(iz = inmin[2]; iz <= inmax[2]; iz++)
      {
          for(iy = inmin[1]; iy <= inmax[1]; iy++)
          {
              for(ix = inmin[0]; ix<= inmax[0]; ix++)
              {
                  ic[0] = ix; ic[1] = iy; ic[2] = iz;
                  if(ic[0] < dimin[0] || ic[0] >= dimax[0] ||
                     ic[1] < dimin[1] || ic[1] >= dimax[1] ||
                     ic[2] < dimin[2] || ic[2] >= dimax[2])
                      return NO;
                  crds = Rect_coords(ic, wave);
                  if((sqr(center[0]-crds[0]) + sqr(center[1]-crds[1]) +
                      sqr(center[2]-crds[2])) > sqr(bubbles.r))
                      continue;

                  if (Rect_comp(ic,wave) != l_comp)
                      return NO;

                  /* Physics checking */
                  state = Rect_state(ic,wave);
                  T_sat = temperature(state);

                  /* THIS IS NOT A GOOF FIX */
                  if(T_sat < 294.261)
                      T_sat = 294.262;

                  /* P_sat = saturated_pres(state,T_sat); */
                  P_sat = 0.5;
                  press = pressure(state);

                  if (press + 0.2 > P_sat)
                      return NO;
              }
          }
      }

      insert_site->center[0] = center[0];
      insert_site->center[1] = center[1];
      insert_site->center[2] = center[2];
      return YES;
}	/* end of check_in_bubble_diameter_3d */


LOCAL int check_crx_in_bubble_space_2d(
      INSERT_SITE  insert_site,
      Wave         *wave,
      Front        *front,
      int          *Rbmin,
      int          *Rbmax,
      double        spacing)
{

      double        cl[MAXD], cu[MAXD];
      RECT_GRID    *gr = wave->rect_grid;
      int          i;

      /* make the domain slightly bigger than bubble radius */
      for(i = 0; i < gr->dim; i++)
      {
          cl[i] = insert_site.center[i]-(spacing*2.0+1.0)*
                 bubbles.size/2.0*min(gr->h[0],gr->h[1]);
          cu[i] = insert_site.center[i]+(spacing*2.0+1.0)*
                 bubbles.size/2.0*min(gr->h[0],gr->h[1]);
      }

      return intfc_on_rect_domain_2d(front->interf, cl, cu);
      /* return intfc_on_rect_domain_ver2(front->interf, cl, cu); */
}	/* end of check_crx_in_bubble_space_2d */

LOCAL void set_insert_sites(
	INSERT_SITE 	 *patch_site,
        INSERT_SITE      *insert_site,
        int              n_site,
        int              max_site,
	int              *n_site_nodes)
{	
        byte          *send_storage, *recv_storage, *buf;	
	POINTER       info;
	int           numnodes, i, j, loc_n;

#if defined(__MPI__)
	numnodes = pp_numnodes();
        scalar(&send_storage,max_site*sizeof(INSERT_SITE));
        scalar(&recv_storage,numnodes*max_site*sizeof(INSERT_SITE));
        buf = send_storage;
        for(i = 0; i < n_site; i++) 
        {
            info = (POINTER) buf;
            ft_assign(info,&insert_site[i],sizeof(INSERT_SITE));
            buf += sizeof(INSERT_SITE);
        }
        pp_all_gather((POINTER)send_storage,max_site*sizeof(INSERT_SITE),
                  (POINTER)recv_storage,max_site*sizeof(INSERT_SITE));
        loc_n = 0;
        buf = recv_storage;
        for(i = 0; i < numnodes; i++)
        {
            for(j = 0; j < n_site_nodes[i]; j++)
            {
                info = (POINTER) buf;
                ft_assign(&patch_site[loc_n],info,sizeof(INSERT_SITE));
                buf += sizeof(INSERT_SITE);
                loc_n++;
            }
            for(j = n_site_nodes[i]; j < max_site; j++)
                buf += sizeof(INSERT_SITE);
        }
        free_these(2,send_storage,recv_storage);
#else /* if defined(__MPI__) */
        for(i = 0; i < max_site; i++)
            ft_assign(&patch_site[i],&insert_site[i],sizeof(INSERT_SITE));
#endif /* if defined(__MPI__) */
}

LOCAL int check_crx_in_bubble_space_3d(
      INSERT_SITE  insert_site,
      Wave         *wave,
      Front        *front,
      int          *Rbmin,
      int          *Rbmax,
      double        spacing)
{

      double        cl[MAXD], cu[MAXD];
      RECT_GRID    *gr = wave->rect_grid;
      int          i;

      /* make the domain slightly bigger than bubble radius */
      for(i = 0; i < gr->dim; i++)
      {
          cl[i] = insert_site.center[i]-(spacing*2.0+1.0)*
                 bubbles.size/2.0*min(gr->h[0],gr->h[1]);
          cu[i] = insert_site.center[i]+(spacing*2.0+1.0)*
                 bubbles.size/2.0*min(gr->h[0],gr->h[1]);
      }

      return intfc_on_rect_domain_3d(front->interf, cl, cu);
      /* return intfc_on_rect_domain_ver2(front->interf, cl, cu); */
}

/* Order xyzxyz */
LOCAL int intfc_on_rect_domain_3d(
      INTERFACE    *interf,
      double        *cl,
      double        *cu)
{
      /* Interface cross with left rectangle boundary */
      if(YES == intfc_intersect_surf(interf, 0, cl[0], cl[1], cu[1], cl[2], cu[2]))
          return YES;      
      /* Interface cross with right rectangle boundary */
      if(YES == intfc_intersect_surf(interf, 0, cu[0], cl[1], cu[1], cl[2], cu[2]))
          return YES;      
      /* Interface cross with front rectangle boundary */
      if(YES == intfc_intersect_surf(interf, 1, cl[1], cl[2], cu[2], cl[0], cu[0]))
          return YES;      
      /* Interface cross with back rectangle boundary */
      if(YES == intfc_intersect_surf(interf, 1, cu[1], cl[2], cu[2], cl[0], cu[0]))
          return YES;      
      /* Interface cross with top rectangle boundary */
      if(YES == intfc_intersect_surf(interf, 2, cl[2], cl[0], cu[0], cl[1], cu[1]))
          return YES;      
      /* Interface cross with bottom rectangle boundary */
      if(YES == intfc_intersect_surf(interf, 2, cu[2], cl[0], cu[0], cl[1], cu[1]))
          return YES;      
      if(YES == intfc_cross_box(interf,cl, cu))
          return YES;
      return NO;
}

LOCAL int make_bubbles_3d(
        Wave                 *wave,
        Front                *front,
        COMPONENT            l_comp,
        COMPONENT            bubblecomp)
{
        RECT_GRID               *gr = front->rect_grid;
        int                     dim = gr->dim;
        int                     icoords[MAXD],imin[MAXD],imax[MAXD];
        int                     bmin[MAXD], bmax[MAXD], Rbmin[MAXD], Rbmax[MAXD];
        int                     ix, iy, iz, i, j;

        INTERFACE               *intfc = front->interf;
        int                     size;
        double                   spacing;
        int                     num_bubble = 0;
        double                   period = 5.e-6;
        Locstate                state;
        double                   T_sat, P_sat, press, ref_temp;
        static Locstate         ref_st = NULL, rst = NULL, distr_st;
        static int              First = YES;
        INSERT_SITE             *insert_site = NULL;
        int                     n_site = 0, alloc_site = 0;
        int                     num_cells, ***distr_site;
        double                   alpha = 1.0;

        debug_print("bubbles_mpi","Entering make_bubbles_3d()\n");
        if(YES == First)
        {
            bubbles.st = NULL;
            bubbles.bubble = NO;
            First = NO;
        }

        bubbles.bubblecomp = bubblecomp;
        if ( bubbles.bubble != YES && loop_bubble == 0)
        {
            if (ref_st == NULL)
            {
                alloc_state(intfc,&ref_st,front->sizest);
                alloc_state(intfc,&rst,front->sizest);
                alloc_state(intfc,&distr_st,front->sizest);
                Params(ref_st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                Params(rst) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_ambient_comp_type(comp_type(bubbles.bubblecomp),front);
                set_type_of_state(ref_st,GAS_STATE);
                set_type_of_state(rst,GAS_STATE);
                zero_state_velocity(ref_st,Params(ref_st)->dim);

                /* values are taken at T =  298.05329710794149 K
                 * P = 1.013 bar
                 * rho = 0.0040948445783132534 g/cm^3
                 * from n-heptane values in gas state
                 */
                Dens(ref_st) = 0.0040948445783132534;
                Energy(ref_st) =  1.013/(1.05 - 1.0);
            }
            if (bubbles.st == NULL)
            {
                alloc_state(intfc,&bubbles.st,front->sizest);
                Params(bubbles.st) = ((G_INTERFACE *) intfc)->g_user_intfc.params_list[2];
                set_type_of_state(bubbles.st,GAS_STATE);
                zero_state_velocity(bubbles.st,Params(bubbles.st)->dim);
            }

            bubbles.bubble = YES;
            /* bubble diameter in mesh size */
            /* bubbles.size = size = 4; */
            bubbles.size = size = bubble_diameter;
            /* mini. distance between two bubble boundary in bubble diameter */
            /* spacing = 1.0; */
            spacing = bubble_spacing;

            bubbles.r = 0.98*size/2.0*min(gr->h[0], gr->h[1]);
            bubbles.N_bubbles = 1;
            bubbles.N_layers = 1;
            bubbles.L = 2 * gr->h[0];
            bubbles.l = 3 * gr->h[0];
            bubbles.eos_n = 2;
            bubbles.eos = YES;

            /* Problem specific init. */

            /* Critical pressure to allow bubble creation */
            /* Or can choose to use Clausis-Claperon in case of cavitation */
            bubbles.pr_crit = -10.0;

            /* Again, Clausis-Claperon can be used to set bubble state */
            /*from n-heptane in gas state*/
            bubbles.pr =  1.013;

            state_on_adiabat_with_pr(ref_st, bubbles.pr, rst, GAS_STATE);
            bubbles.rho = Dens(rst);

            /* EOS parameter should not appear. must remove */
            /*from n-heptane in gas state*/
            bubbles.GAMMA = 1.05 - 1.0;
            bubbles.R = 0.83;
        }

        if  (bubbles.bubble == YES)
        {
            /* domain which allows bubble to be created */
            /* The current domain setting might give insufficient
             * neighboring space if a bubble is created near the boundary, when
             * the spacing is greater than 1 bubble diameter
             * and bubble diameter is greater than buffer size (4).
             * However, to include these in the domain setting might never
             * create bubbles near the boundary.
             */
            for(i = 0; i < dim; i++)
            {
                imin[i] = 0;
                imax[i] = gr->gmax[i];
                /** For conservative insertion, do not allow across reflection line **/
                /*
                if(rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
                    imin[i] -= gr->lbuf[i];
                if(rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
                    imax[i] += gr->ubuf[i];
                */
            }

            alloc_site = max(imax[0],max(imax[1],imax[2]));
            uni_array(&insert_site, alloc_site, sizeof(INSERT_SITE));
            tri_array(&distr_site, gr->gmax[0], gr->gmax[1], gr->gmax[2],sizeof(int));
            for(iz = 0; iz < gr->gmax[2]; iz++)
                for(iy = 0; iy < gr->gmax[1]; iy++)
                    for(ix = 0; ix < gr->gmax[0]; ix++)
                        distr_site[ix][iy][iz] = NO;

            /* Search pressure below specified critical value */
            for (iz = imin[2]; iz < imax[2]; iz++)
            {
                icoords[2] = iz;
                for (iy = imin[1]; iy < imax[1]; iy++)
                {
                    icoords[1] = iy;
                    for (ix = imin[0]; ix < imax[0]; ix++)
                    {
                        icoords[0] = ix;
                        if (Rect_comp(icoords,wave) != l_comp) continue;

                        state = Rect_state(icoords,wave);
                        T_sat = temperature(state);

                        /* P_sat = saturated_pres(state,T_sat); */
                        /* For test purpose */
                        P_sat = 0.5;

                        press = pressure(state);

                        if (press+0.2 < P_sat) /* for simulation sp_mp_jet40_new3.r */
                        {
                            ft_assign(insert_site[n_site].icoords,icoords,INT*MAXD);
                            insert_site[n_site].press = press;
                            insert_site[n_site].insert = NO;
                            n_site++;
                            if    (n_site == alloc_site)
                            {
                                INSERT_SITE *tmp_site;
                                alloc_site += max(imax[0],imax[1]);
                                uni_array(&tmp_site, alloc_site, sizeof(INSERT_SITE));
                                for(i = 0; i < n_site; i++)
                                    ft_assign(&tmp_site[i], &insert_site[i], sizeof(INSERT_SITE));
                                free(insert_site);
                                insert_site = tmp_site;
                            }
                        }
                    }
                }
            }

            /* put insert_site in order according to pressure */
            if(n_site != 0)
                qsort((POINTER)insert_site, n_site,sizeof(INSERT_SITE),press_ascend);

            /* Create bubbles */
            for(i = 0; i < n_site; i++)
            {
                if(YES != check_in_bubble_diameter_3d(&(insert_site[i]),
                    wave,front,l_comp,imin,imax,bmin,bmax))
                    continue;
              
                blks_bubble_ocupy(insert_site[i],gr,Rbmin,Rbmax,spacing);

                /*** Old alg. bubbles are always created in spacing distance to surfaces ***/
                if(YES == check_crx_in_bubble_space_3d(insert_site[i],
                        wave,front,Rbmin,Rbmax,spacing))
                    continue;
 
                avg_states_in_reg_3d(&(insert_site[i]),bmin,bmax,wave,front,ref_st,&num_cells);
 
                /** Compute total quantity needs to be distributed **/
                copy_state(distr_st, ref_st);

                ref_temp = temperature(ref_st);

		/*TMP liuxt12 */
		printf("manually set temperature, should be removed\n",i);
		ref_temp = 273;
		printf("bubble will be inserted at:\n");
		print_general_vector("center",insert_site[i].center,3,"\n");
		print_int_vector("icoords",insert_site[i].icoords,3,"\n");
                set_state_on_phase_boundary_ver3(ref_st, rst, ref_temp);

                create_bubbles_3d(&(insert_site[i]),l_comp,wave,front,rst,bmin,bmax);
		/*TMP liuxt12 */
		printf("bubble %d is created\n",i);
		fflush(stdout);

                insert_site[i].insert = YES;
                for(j = 0; j < dim; j++)
                {
                    insert_site[i].bmin[j] = bmin[j];
                    insert_site[i].bmax[j] = bmax[j];
                }
                num_bubble++;
            }
            
            free(insert_site);
            free(distr_site);
        }

        bubbles.bubble = NO;
        debug_print("bubbles_mpi","Left make_bubbles_3d()\n");
        return num_bubble;

        printf("ERROR: Entered make_bubbles_3d()\n");
	printf("call clean_up in the end for debugging\n");
        clean_up(ERROR);
}	/* end of make_bubbles_3d */

LOCAL  int    intfc_intersect_surf(
        INTERFACE       *intfc,        /* interface to be cut */
        int             dir,           /* direction of cut line normal */
        double           cut,           /* coordinate of cut line */
        double           L1,            /* lower coordinate of cut line in other direction */
        double           U1,
        double           L2,            /* upper coordinate of cut line in other direction */
        double           U2)            
{
        POINT             *p;
        SURFACE           **s;
        TRI               *tri;
        int               num_surfs, num_tris, k;

        for (num_tris = 0, s = intfc->surfaces; s && *s; ++s)
        {
            num_tris += (*s)->num_tri;
            for (tri=first_tri(*s); !at_end_of_tri_list(tri,*s); tri=tri->next)
            {
                if(YES == tri_intersect_surf(tri,dir,cut,L1,U1,L2,U2))
                    return YES;
            }
        }
        return NO;
}

/* The normals of the plane are
 * (1,0,0) (0,1,0), (0,0,1) respectively.
 * This is a simplied version of general 
 * line segment intersection with facet.
 * plane function n_x X + n_y Y + n_z Z + d = 0;
 */
LOCAL  int    tri_intersect_surf(
        TRI        *tri,        
        int        dir,          
        double      cut,           
        double      L1,            
        double      U1,
        double      L2,            
        double      U2)            
{
        int        i;
        double      *p1, *p2, p[MAXD];
        double      nx, ny, nz, d;
        double      denom, mu;

        nx = ny = nz = 0;
        switch(dir)
        {
        case 0:
           nx = 1.0;
           break;
        case 1:
           ny = 1.0;
           break;
        case 2:
           nz = 1.0;
           break;
        default:
              printf("ERROR: tri_intersect_surf, unknown case, dir = %d\n", dir);
              clean_up(ERROR);
        }
        d = -cut;
        
        for (i = 0; i < 3; ++i)
        {
            p1 = Coords(Point_of_tri(tri)[i]); 
            p2 = Coords(Point_of_tri(tri)[Next_m3(i)]);
        }
 
        denom = nx*(p2[0]-p1[0]) + ny*(p2[1]-p1[1]) + nz*(p2[2]-p1[2]);
        if(fabs(denom) < 10.0*MACH_EPS) /* line and plane do not intersect */
            return NO;

        mu = -(d + nx*p1[0] + ny*p1[1] + nz*p1[2])/denom;
        if(mu < 0 || mu > 1) /* intersection not along line segment */        
            return NO;

        for(i = 0; i < 3; i++)
            p[i] = p1[i] + mu*(p2[i]-p1[i]);

        switch(dir)
        {
        case 0:
            if(p[1] <= U1 && p[1] >= L1 && 
               p[2] <= U2 && p[2] >= L2)
                return YES;
            break;
        case 1:
            if(p[2] <= U1 && p[2] >= L1 && 
               p[0] <= U2 && p[0] >= L2)
                return YES;
            break;
        case 2:
            if(p[0] <= U1 && p[0] >= L1 && 
               p[1] <= U2 && p[1] >= L2)
                return YES;
            break;
            
        }
        return NO;
}

LOCAL void avg_states_in_reg_3d(
      INSERT_SITE  *insert_site,
      int    *imin,
      int    *imax,
      Wave   *wave,
      Front  *front,
      Locstate  ref_st,
      int       *num_cells)
{
      int    iy, ix, iz, ic[MAXD];
      double  *crds, *center;
      Locstate  state;
      double     rho, m[MAXD], E;
      int    first = YES;
      RECT_GRID  *gr = wave->rect_grid;
      double      alpha = 1.0;

      Dens(ref_st) = 0.0;
      Mom(ref_st)[0] = 0.0;
      Mom(ref_st)[1] = 0.0;
      Energy(ref_st) = 0.0;

      center = insert_site->center;
      rho = m[0] = m[1] = E = 0;
      *num_cells = 0;

      for(iz = imin[2]; iz <= imax[2]; iz++)
      {
          for(iy = imin[1]; iy <= imax[1]; iy++)
          {
              for(ix = imin[0]; ix<= imax[0]; ix++)
              {
                  ic[0] = ix; ic[1] = iy; ic[2] = iz;
                  crds = Rect_coords(ic, wave);
                  if((sqr(center[0]-crds[0]) + sqr(center[1]-crds[1]) + 
                      sqr(center[2]-crds[2])) >
                     sqr(bubbles.r))
                      continue;

                  state = Rect_state(ic,wave);
                  if(YES == first)
                  {
                      Set_params(ref_st,state);
                      first = NO;
                  }
                  if(Different_params(ref_st,state))
                  {
                      printf("ERROR: avg_states_in_reg, diff params\n");
                      clean_up(ERROR);
                  }
		  if (is_rotational_symmetry())
		  {
                      if (gr->Remap.remap == CYLINDRICAL_REMAP)
                          alpha = crds[0];
                      alpha = 1.0;
		  }
                  Dens(ref_st) += alpha*Dens(state);
                  Mom(ref_st)[0] += alpha*Mom(state)[0];
                  Mom(ref_st)[1] += alpha*Mom(state)[1];
                  Mom(ref_st)[2] += alpha*Mom(state)[2];
                  Energy(ref_st) += alpha*Energy(state);
                  (*num_cells)++;
              }
          }
      }

      if (is_rotational_symmetry())
      {
          if (gr->Remap.remap == CYLINDRICAL_REMAP)
              alpha = center[0];
          alpha = 1.0;
      }

      Dens(ref_st) /= ((*num_cells)*alpha);
      Mom(ref_st)[0] /= ((*num_cells)*alpha);
      Mom(ref_st)[1] /= ((*num_cells)*alpha);
      Mom(ref_st)[2] /= ((*num_cells)*alpha);
      Energy(ref_st) /= ((*num_cells)*alpha);
}

LOCAL void create_bubbles_3d(
      INSERT_SITE  *insert_site,
      COMPONENT    l_comp,
      Wave         *wave,
      Front        *front,
      Locstate     rst,
      int          *imin,
      int          *imax)
{
      double        r,theta,d_theta,step;
      int          n_theta, ic[MAXD], ix, iy, iz, k, i;
      double        *c, coords[MAXD];
      Locstate     state;
      RECT_GRID    *gr = wave->rect_grid;
      INTERFACE    *intfc = front->interf;
      INTERFACE    *sv_intfc;
      ELLIPSOID    ellip;
      HYPER_SURF   *hs;
      TRI          *tri;
      POINT        *pt;

      debug_print("bubbles_mpi","Entered create_bubbles()\n");

      sv_intfc = current_interface();
      set_current_interface(intfc);
      c = insert_site->center;
      r = bubbles.r;

      ellip.compin = bubbles.bubblecomp;
      ellip.compout = l_comp;
      ellip.fpoly = NULL;
      for(k = 0; k < gr->dim; k++)
      {
          ellip.cen[k] = c[k];
          ellip.rad[k] = r;
      }

      hs = g_make_ellipsoid(&ellip, ellip.compin, ellip.compout, front);

      wave_type(hs) = CONTACT;
      create_time(hs) = front->time;

      for (tri=first_tri(Surface_of_hs(hs)); !at_end_of_tri_list(tri,Surface_of_hs(hs)); 
           tri=tri->next)
      {
          for(i = 0; i < 3; i++)
          {
              pt = Point_of_tri(tri)[i];
              front_init_3d(pt, hs, left_state(pt),right_state(pt), front,wave,rst);
          }
      }

      /*init interior states for bubble*/
      k = 0;
      for ( iz = imin[2]; iz <= imax[2]; iz++)
      {
          ic[2] = iz;
          for ( iy = imin[1]; iy <= imax[1]; iy++)
          {
              ic[1] = iy;
              for ( ix = imin[0]; ix <= imax[0]; ix++)
              {
                  int     stype;
                  double   dist;

                  ic[0] = ix;
                  coords[0] = Rect_coords(ic,wave)[0];
                  coords[1] = Rect_coords(ic,wave)[1];
                  coords[2] = Rect_coords(ic,wave)[2];
                  if( (dist = (sqr(c[0]-coords[0]) + sqr(c[1]-coords[1]) + 
                               sqr(c[2]-coords[2]))) > sqr(r))
                      continue;
                  state = Rect_state(ic,wave);
                  if (debugging("bubbles_mpi"))
                  {
                      (void) printf("Original, comp = %d, bub_comp %d ",
                          Rect_comp(ic,wave),bubbles.bubblecomp);
                      printf("coords = (%g, %g), ",coords[0],coords[1]);
                      printf("distance = %g <= %g\n", sqrt(dist), bubbles.r);
                      verbose_print_state("state",state);
                  }
                  stype = state_type(state);
                  Rect_comp(ic,wave) = bubbles.bubblecomp;
                  copy_state(state,bubbles.st);
                  set_state(state,stype,state);
                  k++;
              }
          }
      }
      if (debugging("bubbles_mpi"))
          printf("End of initilization of bubbles in %d cells\n", k);

      set_current_interface(sv_intfc);
      debug_print("bubbles_mpi","Left create_bubbles_3d()\n");
}

LOCAL void front_init_3d(
        POINT      *p,
        HYPER_SURF *hs,
        Locstate   lstate,
        Locstate   rstate,
        Front      *front,
        Wave       *wave,
        Locstate   rst)
{

        COMPONENT         pcomp = NO_COMP, ncomp = NO_COMP;
        INTERFACE         *intfc = hs->interface;
        int               dim = intfc->dim;

        int               stype;
        int               ic[MAXD];
        int               i, j;
        double             coordp[MAXD];
        Locstate          state;
        double             dens;
        double             v[MAXD];
        double             pl, tl;

        if (debugging("bubbles_mpi"))
        {
          printf("\nfront_init1\n");
        }

        ncomp = negative_component(hs);
        pcomp = positive_component(hs);

        if (debugging("bubbles_mpi"))
        {
            printf("nct->type = %d, EXTERIOR = %d, \n",comp_type(ncomp)->type, EXTERIOR);
            printf("In default, ncomp = %d, bubbles.bubblecomp = %d\n"
                   ,ncomp, bubbles.bubblecomp);
        }

        if (ncomp != bubbles.bubblecomp)
        {
            if ( debugging("from_nbhd") )
            {
                rect_in_which(Coords(p),ic,wave->rect_grid);
                for(i=0; i < dim; i++)
                    coordp[i] = Rect_coords(ic,wave)[i];
                if ( Rect_comp(ic,wave) != 2 )
                {
                    printf("ERROR in front_init(), all points on the curve of "
                           "the bubble must be liquid, not comp = %d\n", Rect_comp(ic,wave));
                    clean_up(ERROR);
                }
                state = Rect_state(ic,wave);
                stype = state_type(state);
                copy_state(lstate,state);
            }
            else
                copy_state(lstate,rst);

            if ( debugging("bubbles_mpi") )
            {
                verbose_print_state("front_state(), lstate",lstate);
            }
        }
	else
        {
            copy_state(lstate,bubbles.st);

            if (debugging("bubbles_mpi"))
            {
                verbose_print_state("front_state(), lstate",lstate);
            }
        }

        if (debugging("bubbles_mpi"))
        {
            printf("pct->type = %d, BUBBLE = %d, \n",comp_type(pcomp)->type, BUBBLE);
            printf("In default, pcomp = %d, bubbles.bubblecomp = %d\n"
                   ,pcomp, bubbles.bubblecomp);
        }

        if (pcomp == bubbles.bubblecomp)
        {
            copy_state(rstate,bubbles.st);

            if (debugging("bubbles_mpi"))
            {
                verbose_print_state("front_state(), rstate",rstate);
            }
        }
	else
        {
            if ( debugging("from_nbhd") )
            {
                rect_in_which(Coords(p),ic,wave->rect_grid);
                for(i=0; i < dim; i++)
                    coordp[i] = Rect_coords(ic,wave)[i];
                if ( Rect_comp(ic,wave) != 2 )
                {
                    printf("ERROR in front_init(), all points on the curve of "
                           "the bubble must be liquid, not comp = %d\n", Rect_comp(ic,wave));
                    clean_up(ERROR);
                }
                state = Rect_state(ic,wave);
                stype = state_type(state);
                copy_state(rstate,state);
            }
            else
                copy_state(rstate,rst);

            if ( debugging("bubbles_mpi") )
            {
                verbose_print_state("front_state(), rstate",rstate);
            }
        }

        if (debugging("bubbles_mpi"))
        {
            printf("Bubbles:\n");
            print_general_vector("negative side state at ",Coords(p),dim,"");
            (void) printf(" on hs %llu\n",hypersurface_number(hs));
            print_gas_state(lstate);
            print_general_vector("positive side state at ",Coords(p),dim,"");
            (void) printf(" on hs %llu\n",hypersurface_number(hs));
            print_gas_state(rstate);
            printf("\nLeft front_init\n");
        }
        debug_print("bubbles_mpi","Left front_init()\n");
}       /* end of front_init_3d */


#endif /* if defined(TWOD) || defined(THREED) */

