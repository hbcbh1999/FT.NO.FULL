/*
*                               gvisc.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       The difference between Euler and Navier-Stokes equations are the
*       viscosity and heat conduction source terms. 
*
*       The viscosity terms involve velocity Laplace for incompressible
*       fluid and velocity cross derivative terms for compressible fluid.
*       The heat conduction follows Fourier Law; the heat flux is proportional 
*       to temperaure gradient.  We assume the molecular viscosity coefficient
*       mu and the second viscosity coefficient lambda have the relation      
*       lambda = -2/3 mu by Stokes' hypothesis to simplify the equations.
*       
*       This file provides the finite difference loop to calculate the Laplace
*       and cross derivative terms of velocity by using the state data from
*       Euler equation solver.  Also the heat conduction term of temperature
*       gradient.
*                                        
*/

#include <ghyp/ghyp.h>
#include <ghyp/ghypprotos.h>
#include <gdecs/gdecs.h>
#include <gdecs/geos.h>
#include <gdecs/gstate.h>

#define  twothird  0.666666666666666666666667     /* value of (2.0/3.0) */
#define  MAX_NCOMPS 2

typedef struct {
        int             npts;  /* n point stencil, not used now */
        float           dt;
#if defined ONED
        int             **icoords1d;
        Locstate        *sts1d;
#endif /* ONED */
#if defined TWOD
        int             ***icoords2d;
        Locstate        **sts2d;
#endif /* TWOD */
#if defined THREED
        int             ****icoords3d;
        Locstate        ***sts3d;
#endif /* THREED */
        Locstate        ans;
        Front           *fr, *newfr;
        struct _Wave    *wave, *newwave;
} Pstencil;    /* stencil for parabolic solvers */

	/* LOCAL Function Declarations */
LOCAL bool	g_compute_NS_terms(Locstate);
LOCAL float	Fourier_heat_conduction(Pstencil*);
LOCAL float	diffusive_mass_flux(float,float,float);
LOCAL void	g_ns_soln(Pstencil*);
LOCAL void	fill_1d_3pt_Pstencil(int,Pstencil*);
LOCAL void	fill_2d_9pt_Pstencil(int,int,Pstencil*);
LOCAL void	fill_3d_27pt_Pstencil(int,int,int,Pstencil*);
LOCAL int     	g_neumann_bdry_state_beta(float*,COMPONENT,POINT*,
			HYPER_SURF*,Front*,POINTER,Locstate);
LOCAL void    	find_ghost_state(int*,int*,Locstate,Front*,Wave*);
LOCAL bool 	cross_line_segments(float*,float*,float*,float*,float*);
LOCAL bool 	cross_line_segment_triangle(float*,float*,float*,float*,
			float*,float*);
LOCAL void      parab_front2d(float,Wave*,Wave*,Front*,Front**);
LOCAL void      parab_front3d(float,Wave*,Wave*,Front*,Front**);
LOCAL void      parab_front2d_update(Front*,Wave*,Wave*,CURVE*,CURVE*,float);;
LOCAL void      parab_front3d_update(float,Wave*,Wave*,Front*,Front*);
LOCAL void      parab_front_point(Front*,Wave*,POINT*,POINT*,
                                  HYPER_SURF_ELEMENT*,HYPER_SURF*,float);
LOCAL void      parab_point_update(Front*,Locstate);
LOCAL void      one_side_parab_tan_solver2d(float,float,Tan_stencil*,Locstate,Front*);
LOCAL void      one_side_parab_tan_solver3d(float,float,Tan_stencil*,Locstate,Front*);
LOCAL   void    set_propagation_bounds(Front*,float*,float*);
LOCAL   bool    out_of_bound(POINT*,float*,float*,int);

#if defined(SUBGRID)
LOCAL void      SGS2d(float, Front*, Wave*);
LOCAL void      SGS3d(float, Front*, Wave*);
LOCAL void      mesh_avr(Wave*,float***,float***,int*, int*, int*);
LOCAL void      mesh_avr_mc(Wave*,float****,float****,int*, int*, int*);
#endif /* defined SUBGRID */
LOCAL void      cal_volume_fraction(float, Front*, Wave*);

/*
*			parab_driver():
*/

EXPORT	int parab_driver(
	float		dt,
	float		*dt_frac,
	Wave		*wave,
	Front		*front)
{
	Wave		*newwave;
        Front           *newfront;
	int		i, dim = front->interf->dim;
	int             *iperm; /* permutation of {0,...,dim-1} */
        INTERFACE       *intfc = front->interf;
        Gas_param       **params;
        params = gas_params_list(intfc);
	static char	warn[] = "WARNING in parab_driver()";
	static char	err[] = "ERROR in parab_driver()";
	DEBUG_ENTER(parab_driver)

	debug_print("parab","Entered parab_driver()\n");

	debug_front("front","in parabolic solver",front);
	if( debugging("parab_states") )
	{
	    (void) printf("States before calling parab_solver:\n\n");
	    (void) printf("Front\n");
	    graph_front_states(front);
	    (void) printf("wave %p\n",(POINTER)wave);
	    (*wave->show_wave_states)(wave);
	}

			/* parabolic solver */

	iperm = set_iperm(front->step,dim);

		/* Initialize Intermediate Storage for States */
        fprintf(stdout, "LOCAL AVE SGS\n");
        SGS(dt,front,wave);
        int recon_step = 10;
        if ( (front->step)%recon_step == 1)
            cal_volume_fraction(dt,front,wave);

	newwave = copy_wave(wave);
	clear_wave_pointers(newwave);
	assign_wave_parameters(newwave,wave);
	if( !copy_hyp_solution_function(wave,newwave) )
	{
	    screen("%s, copy_hyp_solution_function() failed\n",err);
	    free_wave(newwave);
	    DEBUG_LEAVE(parab_driver)
	    return ERROR_IN_STEP;
	}

	start_clock("parab_solver");

        parab_npt(dt,front,wave,newwave);

	start_clock("scatter_states");

	iperm = set_iperm(front->step,dim);
	for (i = 0; i < dim; i++)
        {
            if (!scatter_states(newwave,front,iperm,i))
            {
                screen("scatter_states() failed in direction\n",i);
                clean_up(ERROR);
            }
        }

	stop_clock("scatter_states");

        /*
        start_clock("front_parab_solver");
        switch (dim)
        {
#if defined ONED
        case 1:
        {
            //parab_front1d(dt,wave,newwave,front,&newfront);
        }
        break;

#endif // defined ONED 
#if defined TWOD
        case 2:
        {
            parab_front2d(dt,wave,newwave,front,&newfront);
        }
        break;

#endif // defined TWOD 
#if defined THREED
        case 3:
        {
            parab_front3d(dt,wave,newwave,front,&newfront);
        }
        break;

#endif // defined THREED 

        }
        reinit_hyp_solution_function(newwave,newfront);
        assign_interface_and_free_front(front,newfront);
        stop_clock("front_parab_solver");
        */
                /* Copy updated front and wave */

	assign_copy_wave_pointers(wave,newwave);
	free_wave(newwave);

	stop_clock("parab_solver");
	if( debugging("parab_states") )
	{
	    printf("States after calling parab_solver:\n");
	    if (wave->show_tri_soln)
		(*wave->show_tri_soln)(front,wave);
	    (*wave->show_wave_states)(wave);
	}

	debug_print("parab","Left parab_driver()\n");
	DEBUG_LEAVE(parab_driver)
	return GOOD_STEP;
}		/*end parab_driver*/

        /* for parabolic solver */

EXPORT	void    parab_npt(
        float   dt,
        Front   *front,
        Wave    *wv,
        Wave    *newwv)
{
        RECT_GRID       *rgr = wv->rect_grid;
        int             dim = rgr->dim;
        int             *gmax = rgr->gmax;
        int             i,j,k;
        int             icoords[MAXD];
        static Pstencil *nsten;
        int             npts,mid;
        Locstate        ans;
        int             imin[MAXD],imax[MAXD];
	int		*iperm;
	size_t		sizest = front->sizest;

	debug_print("parab","Entered parab_npt()\n");
        if (nsten == NULL)
        {
            stat_scalar(&nsten,sizeof(Pstencil));
            switch (dim)
            {
            case 1:
                bi_array(&nsten->icoords1d,5,1,INT);
                uni_array(&nsten->sts1d,5,sizeof(Locstate));
		for (i = 0; i < 5; i++)
		    alloc_state(front->interf,&nsten->sts1d[i],front->sizest);
                break;
            case 2:
                tri_array(&nsten->icoords2d,3,3,2,INT);
                bi_array(&nsten->sts2d,3,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                {
                    for (j = 0; j < 3; j++)
                    {
                        alloc_state(front->interf,&nsten->sts2d[i][j],
                                front->sizest);
                    }
                }
                break;
            case 3:
                quad_array(&nsten->icoords3d,3,3,3,3,INT);
                tri_array(&nsten->sts3d,3,3,3,sizeof(Locstate));
		for (i = 0; i < 3; i++)
		    for (j = 0; j < 3; j++)
		        for (k = 0; k < 3; k++)
			    alloc_state(front->interf,&nsten->sts3d[i][j][k],
			            front->sizest);
                break;
            }
        }
        nsten->fr = front;
        nsten->wave = wv;
        nsten->newwave = newwv;
        nsten->npts = npts = wv->npts_sten;
        mid = (int) (npts/2);

        for (i = 0; i < dim; i++)
        {
            imin[i] = 0;        imax[i] = gmax[i];
        }
        switch (dim)
        {
#if defined ONED
        case 1:
        {
            int i0;

            for (i = imin[0]; i < imax[0]; i++)
            {
                icoords[0] = i;
		ft_assign(Rect_state(icoords,newwv),
		       Rect_state(icoords,wv),sizest);
		
		if(g_compute_NS_terms(Rect_state(icoords,wv)))
		{
		   fill_1d_3pt_Pstencil(i,nsten);
		   nsten->ans = Rect_state(icoords,newwv);
		   g_ns_soln(nsten);
		}
		
		/*
                for (i0 = 0; i0 < 3; i0++)
                {
                    nsten->icoords1d[i0][0] = i + i0 - 1;
                }
                for (i0 = 0; i0 < 3; i0++)
                {
                    nsten->sts1d[i0] =
                        Rect_state(nsten->icoords1d[i0],wv);
                }
                icoords[0] = i;
                nsten->ans = Rect_state(icoords,newwv);
		*/
            }
        }
        break;

#endif /* defined ONED */
#if defined TWOD
        case 2:
        {
	    int i0, i1, ixl, ixr, iy, icl[2], icr[2];

	    float *h, *L, *U, coordsl[2], coordsr[2], *coords_l, *coords_r; 

            for (j = imin[1]; j < imax[1]; j++)
            {
                icoords[1] = j;
                for (i = imin[0]; i < imax[0]; i++)
                {
                    icoords[0] = i;
		    ft_assign(Rect_state(icoords,newwv),
			   Rect_state(icoords,wv),sizest);
		    if(g_compute_NS_terms(Rect_state(icoords,wv)))
		    {
                    	fill_2d_9pt_Pstencil(i,j,nsten);
                    	nsten->ans = Rect_state(icoords,newwv);
		    	g_ns_soln(nsten);
		    }
                }
            }
        }
        break;

#endif /* defined TWOD */
#if defined THREED
        case 3:
        {
            for (k = imin[2]; k < imax[2]; k++)
            {
              icoords[2] = k;
              for (j = imin[1]; j < imax[1]; j++)
              {
                 icoords[1] = j;
                 for (i = imin[0]; i < imax[0]; i++)
                 {
                    icoords[0] = i;
		    ft_assign(Rect_state(icoords,newwv),
		       	   Rect_state(icoords,wv),sizest);
		    if(g_compute_NS_terms(Rect_state(icoords,wv)))
		    {
                        //fprintf(stdout, "g_ns %d %d %d \n",i,j,k);
			fill_3d_27pt_Pstencil(i,j,k,nsten);
			nsten->ans = Rect_state(icoords,newwv);
			g_ns_soln(nsten);
		    }
                 }
              }
            }
        }
        break;

#endif /* defined THREED */
        }
	debug_print("parab","Left parab_npt()\n");
}	/* end parab_npt */

LOCAL void fill_1d_3pt_Pstencil(
		int             i,
		Pstencil        *nsten)
{
        RECT_GRID   *rgr = nsten->wave->rect_grid;
	int         j, i0, dim = nsten->wave->rect_grid->dim;
	int         icoords[MAXD], ic[MAXD];
	Wave        *wv = nsten->wave;
	Front       *front = nsten->fr;
	COMPONENT   cc,cck;
	size_t      sizest = nsten->fr->sizest;
	HYPER_SURF  *hs = NULL;
	float       coords[MAXD], coords_on[MAXD], t[MAXD];
	HYPER_SURF_ELEMENT *hse = NULL;
	Locstate        state;

	icoords[0] = i;
	cc = Rect_comp(icoords,wv);
	for (i0 = 0; i0 < 3; i0++)
	{
	    nsten->icoords1d[i0][0] = i + i0 - 1;
	    ic[0] = i + i0 - 1;
	    state = nsten->sts1d[i0];
	    cck =  Rect_comp(nsten->icoords1d[i0],wv);
	    if (is_excluded_comp(cck,front->interf))
	    {
		coords[0] = cell_center(nsten->icoords1d[i0][0],0,rgr);
		if (!nearest_interface_point(coords,cc,front->interf,
					NO_SUBDOMAIN,NULL,coords_on,
					t,&hse,&hs))
		{
		    printf("ERROR in fill_1d_3pt_Pstencil(), "
				    "can't find nearest interface point\n");
		    printf("ic[0] = %d, crds[0] = %g\n", ic[0], coords[0]);
		    print_interface(front->interf);
		    clean_up(ERROR);
		}
		printf("third test in the loop\n");
		switch (wave_type(hs))
		{
		case CONTACT:
		    hyp_solution(coords,cc,hs,UNKNOWN_SIDE,front,
				     wv,state,
			hs->pos_comp == cc ? right_state(Point_of_hs(hs)):
			   left_state(Point_of_hs(hs)));
		    break;
		case NEUMANN_BOUNDARY:
		    g_neumann_bdry_state_beta(coords,cc,0,hs,front,
				    (POINTER)wv,state);
		    break;
		case DIRICHLET_BOUNDARY:    /* far field conditions */
		    evaluate_dirichlet_boundary_state(coords,hs,
				    front,wv,state);
		    break;
		case PASSIVE_BOUNDARY:
		case SUBDOMAIN_BOUNDARY:
		default:
		    screen("ERROR in fill_1d_3pt_Pstencil(), "
				    "unknown boundary type\n");
		    screen("\twave type = %s\n"
				   ,wave_type_as_string(wave_type(hs)
				,hs->interface));
		    screen("\n\tcoord = (%g)\n",coords[0]);
		    clean_up(ERROR);
		    break;
		}
	    }
	    else
            {
		ft_assign(state, Rect_state(ic,wv), sizest);
	    }
	}
}  /* end Left fill_1d_3pt_Pstencil */

LOCAL void fill_2d_9pt_Pstencil(
                int             i,
                int             j,
                Pstencil        *nsten)
{
        RECT_GRID   *rgr = nsten->wave->rect_grid;
        int         ii, jj, dim = nsten->wave->rect_grid->dim;
        int         gmax[MAXD], gmin[MAXD];
        int         icoords[MAXD], ic[MAXD];
        Wave        *wv = nsten->wave;
        Front       *front = nsten->fr;
        COMPONENT   cc,cck;
        size_t      sizest = nsten->fr->sizest;
        HYPER_SURF  *hs = NULL;
        float       coords[MAXD], coords_on[MAXD], t[MAXD];
        HYPER_SURF_ELEMENT *hse = NULL;
        Locstate        state;

        icoords[0] = i;
        icoords[1] = j;
        cc = Rect_comp(icoords,wv);
        for(jj = 0; jj < 3; jj++)
        {
            for (ii = 0; ii < 3; ii++)
            {
                nsten->icoords2d[ii][jj][0] = i + ii - 1;
                nsten->icoords2d[ii][jj][1] = j + jj - 1;
                ic[0] = i + ii - 1;
                ic[1] = j + jj - 1;
                state = nsten->sts2d[ii][jj];
                cck =  Rect_comp(nsten->icoords2d[ii][jj],wv);
                if (is_excluded_comp(cck,front->interf))
                {
                    coords[0] = cell_center(nsten->icoords2d[ii][jj][0],0,rgr);
                    coords[1] = cell_center(nsten->icoords2d[ii][jj][1],1,rgr);
                    if (!nearest_interface_point(coords,cc,front->interf,
                                            NO_SUBDOMAIN,NULL,coords_on,
                                            t,&hse,&hs))
                    {
                        screen("ERROR in fill_2d_9pt_Pstencil(), "
                               "can't find nearest interface point\n");
                        clear_state(front->interf,state,front->sizest);
                        clean_up(ERROR);
                    }
                    switch (wave_type(hs))
                    {
                    case CONTACT:
                        hyp_solution(coords,cc,hs,UNKNOWN_SIDE,front,
                                    wv,state,NULL);
                        break;
                    case NEUMANN_BOUNDARY:
                        g_neumann_bdry_state_beta(coords,cc,0,hs,front,
                                        (POINTER)wv,state);
                        break;
                    case DIRICHLET_BOUNDARY:    /* far field conditions */
                        evaluate_dirichlet_boundary_state(coords,hs,
                                        front,wv,state);
                        break;
                    case PASSIVE_BOUNDARY:
                    case SUBDOMAIN_BOUNDARY:
                    default:
                        screen("ERROR in fill_2d_9pt_stencil_new(), "
                               "unknown boundary type\n");
                        screen("\twave type = %s\n"
                               ,wave_type_as_string(wave_type(hs)
                                                    ,hs->interface));
                        screen("\n\tcoord = (%g,%g)\n",coords[0],coords[1]);
                        clean_up(ERROR);
                        break;
                    }
                }
                else
                {
                    ft_assign(state, Rect_state(ic,wv), sizest);
                }
            }
        }
}       /* end fill_2d_9pt_stencil_new */

LOCAL void fill_3d_27pt_Pstencil(
                int             i,
                int             j,
                int             k,
                Pstencil        *nsten)
{
        RECT_GRID   *rgr = nsten->wave->rect_grid;
        int         *gmax = rgr->gmax;
        int         ii, jj, kk, dim = nsten->wave->rect_grid->dim;
        //int         gmax[MAXD], gmin[MAXD];
        int         icoords[MAXD], ic[MAXD];
        Wave        *wv = nsten->wave;
        Front       *front = nsten->fr;
        COMPONENT   cc,cck;
        size_t      sizest = nsten->fr->sizest;
        HYPER_SURF  *hs = NULL;
        float       coords[MAXD], coords_on[MAXD], t[MAXD];
        HYPER_SURF_ELEMENT *hse = NULL;
        Locstate        state;

        const int  nn = pp_numnodes();
        int        myid = pp_mynode();
        int   *ppgmax = front->pp_grid->gmax;
        int   ppx = myid % ppgmax[0];
        int   ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int   ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);

        icoords[0] = i;
        icoords[1] = j;
        icoords[2] = k;
        cc = Rect_comp(icoords,wv);
        for(kk = 0; kk < 3; kk++)
        {
            for(jj = 0; jj < 3; jj++)
            {
                for (ii = 0; ii < 3; ii++)
                {
                    nsten->icoords3d[ii][jj][kk][0] = i + ii - 1;
                    nsten->icoords3d[ii][jj][kk][1] = j + jj - 1;
                    nsten->icoords3d[ii][jj][kk][2] = k + kk - 1;
                    ic[0] = i + ii - 1;
                    ic[1] = j + jj - 1;
                    ic[2] = k + kk - 1;
                    state = nsten->sts3d[ii][jj][kk];
                    cck =  Rect_comp(nsten->icoords3d[ii][jj][kk],wv);
                    if (is_excluded_comp(cck,front->interf))
                    {
                        coords[0] = cell_center(nsten->icoords3d[ii][jj][kk][0],0,rgr);
                        coords[1] = cell_center(nsten->icoords3d[ii][jj][kk][1],1,rgr);
                        coords[2] = cell_center(nsten->icoords3d[ii][jj][kk][2],2,rgr);
                        if (!nearest_interface_point(coords,cc,front->interf,
                                                NO_SUBDOMAIN,NULL,coords_on,
                                                t,&hse,&hs))
                        {
                            screen("ERROR in fill_3d_27pt_Pstencil(), "
                               "can't find nearest interface point\n");
                            clear_state(front->interf,state,front->sizest);
                            clean_up(ERROR);
                        }
                        switch (wave_type(hs))
                        {
                        case CONTACT:
                            hyp_solution(coords,cc,hs,UNKNOWN_SIDE,front,
                                        wv,state,NULL);
                            break;
                        case NEUMANN_BOUNDARY:
                            g_neumann_bdry_state_beta(coords,cc,0,hs,front,
                                        (POINTER)wv,state);
                            break;
                        case DIRICHLET_BOUNDARY:    /* far field conditions */
                            evaluate_dirichlet_boundary_state(coords,hs,
                                            front,wv,state);
                            break;
                        case PASSIVE_BOUNDARY:
                        case SUBDOMAIN_BOUNDARY:
                        default:
                            screen("ERROR in fill_3d_27pt_stencil_new(), "
                               "unknown boundary type\n");
                                screen("\twave type = %s\n"
                               ,wave_type_as_string(wave_type(hs)
                                                    ,hs->interface));
                            screen("\n\tcoord = (%g,%g,%g)\n",coords[0],coords[1],coords[2]);
                            clean_up(ERROR);
                            break;
                        }
                    }
                    else
                    {
                        ft_assign(state, Rect_state(ic,wv), sizest);
                    }
                }
            }
	}
}       /* end fill_3d_27pt_stencil_new */

LOCAL  void  g_ns_soln(
	Pstencil        *nsten)
{
	RECT_GRID       *gr = nsten->wave->rect_grid;
	int             dim = gr->dim;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF      *hs;
	float		dt = nsten->fr->dt;
	float           diff_rho = 0.0, d_new, d_old;
	Locstate	ans = nsten->ans;
	COMPONENT       cc0,cc1,cc2;
	float           time =  nsten->fr->time, min_dist, dist;
	int             step = nsten->fr->step;
	POINT           *closet_point,*old_point;
	INTERFACE       *intfc = nsten->fr->interf;
	float           *coords,dens_old,x,d_jump;
	int             i,icoords[MAXD];
	bool            FD = NO;
	bool		viscosity, mass_diffusion, thermal_conduction;
	bool            use_stokes_vis;
	bool            subgrid_vis, subgrid_md, subgrid_con;
	float		mu, d, kappa;
	Gas_param       **params;
	bool            subgrid;
        float           VIS_FORCE;

/*
 *      The state array sts3d[3][3][3] is a moving cube in the general three
 *      dimensional space, but also could be applied in 2D and 1D space.
 *      For the first, second, and cross derivatives of velocity and Laplacian,
 *      three point stencil will produce second order accuracy computation.
 *      s[1][1][1] is the center of the cube, which is moving through all
 *      computational grids to get all velocity derivatives of all grids.
 */

        //Get NS Viscosity and mass diffusion coefficients
	params =    gas_params_list(intfc); 

	mu = params[0]->avisc.viscosity_coef;
	d = params[0]->avisc.diffusivity_coef;
        kappa = params[0]->avisc.conductivity_coef;
	use_stokes_vis = params[0]->avisc.use_stokes_vis;
        subgrid_vis = params[0]->avisc.subgrid_vis;
        subgrid_md = params[0]->avisc.subgrid_md;
        subgrid_con = params[0]->avisc.subgrid_con;

	if( fabs(mu) < 1e-10)
	    viscosity = NO;
	else
	    viscosity = YES;

	if( fabs(d) < 1e-10)
	    mass_diffusion = NO;
	else
	    mass_diffusion = YES;

        if(fabs(kappa) <= 1e-16)
            thermal_conduction = NO;
        else
            thermal_conduction = YES;

	if (is_obstacle_state(ans))
	  return;

	//mu = shear_viscosity(ans);
	//kappa = heat_coeff(ans);

	for (i = 0; i < dim; i++)
	{
	    if (dim == 3)
	        icoords[i] = nsten->icoords3d[1][1][1][i];
	    if (dim == 2)
	        icoords[i] = nsten->icoords2d[1][1][i];
	    if (dim == 1)
		icoords[i] = nsten->icoords1d[1][i];
	}
	coords = Rect_coords(icoords,nsten->wave);
	//set_max_viscosity(mu/density(ans),ans,coords,nsten->wave);
	
	switch (dim)
        {
#if defined ONED
        case 1:
        {
            float           ux;    /* first derivative of vel u */
            float           uxx;   /* second derivative of vel u */
	    float	    dh[1];
	    float           rhoxx,distance=HUGE;
	    int             sign;
	    float           *tmp_coords,point_coords[3],rho_r,rho_l;
	    int             icoords0[MAXD],icoords2[MAXD];
	    

	   
	    sign = 0;
	    rho_r = 1.0;
	    rho_l = 3.0;
	    //rho_r = Dens(nsten->sts1d[2]); 
	    //rho_l = Dens(nsten->sts1d[0]);
	    icoords0[0] = nsten->icoords1d[0][0];
	    icoords2[0] = nsten->icoords1d[2][0];
	    cc0 = Rect_comp(icoords0,nsten->wave);
	    cc1 = Rect_comp(icoords,nsten->wave);
	    cc2 = Rect_comp(icoords2,nsten->wave);
            tmp_coords = Rect_coords(icoords0,nsten->wave);
	    point_coords[0] = tmp_coords[0];
	    tmp_coords = Rect_coords(icoords2,nsten->wave);
	    point_coords[2] = tmp_coords[0];
	    point_coords[1] = coords[0];
	    
	    
	    dh[0] = nsten->wave->rect_grid->h[0];


            if (PI*sqrt(d*time) > dh[0])
	    {
	        //printf("switch to FD\n");
		//printf("time = %g\n",time);
		subgrid = NO;
		FD = YES;
	    }
	    if (cc2 != cc1)
	    {
		(void) (next_point(intfc,NULL,NULL,NULL));
		while (next_point(intfc,&old_point,&hse,&hs))
		{
		    if (wave_type(hs) != CONTACT)
		        continue;
		    if (Coords(old_point)[0] < point_coords[2] &&
		        Coords(old_point)[0] > point_coords[1])
		    {
		        sign = -1;
			//rho_l = Dens(left_state(old_point));
			//rho_r = Dens(right_state(old_point));
			distance =  point_coords[1] - Coords(old_point)[0];
			break;
		    }
		}

	    }
	    else if (cc1 != cc0)
	    {
		(void) (next_point(intfc,NULL,NULL,NULL));
		while (next_point(intfc,&old_point,&hse,&hs))
		{
		    if (wave_type(hs) != CONTACT)
		        continue;
		    if (Coords(old_point)[0] < point_coords[1] &&
		        Coords(old_point)[0] > point_coords[0])
		    {
		         sign = 1;
			 //rho_l = Dens(left_state(old_point));
			 //rho_r = Dens(right_state(old_point));
			 distance =  point_coords[1] - Coords(old_point)[0];
			 break;
		    }
		}
	    }
	    ux = (vel(0,nsten->sts1d[2]) - vel(0,nsten->sts1d[0]))/(2.0*dh[0]);

            uxx = (vel(0,nsten->sts1d[2]) - 2.0*vel(0,nsten->sts1d[1])
                   + vel(0,nsten->sts1d[0]))/(dh[0]*dh[0]);

            rhoxx = (Dens(nsten->sts1d[2]) - 2.0*Dens(nsten->sts1d[1])
		   + Dens(nsten->sts1d[0]))/(dh[0]*dh[0]);

	    if (subgrid)
	    {
	        d_new = PI*sqrt(d*time); 
		if (time > 0)
		    d_old = PI*sqrt(d*(time-dt));
		else
	            d_old = 0.0;
	        d_jump = rho_l-rho_r;

	        if (sign != 0)
		    printf("point_coords = %g\n",point_coords[1]);
		//if ((time > 0) && (fabs(distance) < PI*sqrt(d*time)))
		    diff_rho = sign*d_jump*(d_new-d_old)/(4.0*dh[0]);
		//else
		    //diff_rho = 0;
		/*
		if (time > 0 && fabs(distance) < PI*sqrt(d*time))
	            diff_rho = sign*rhoxx*d*dt;
		else
		    diff_rho = 0;
                */ 
	    }
	    if (FD)
	    {
	        diff_rho = rhoxx*d*dt;
	    }
	    //printf("density =%g, rhoxx = %g\n",Dens(ans),rhoxx);
	    Dens(ans) += diff_rho;
	    Mom(ans)[0] += diff_rho*vel(0,nsten->sts1d[1]);
	    //printf("distance= %g,density = %g, coods = %g\n",distance,
			    //Dens(ans),coords[0]);
			  
        }
	break;
	 
#endif /* defined ONED */
#if defined TWOD
	case 2:
	{
	    float           ux, uy;   /* first derivative of vel u */
            float           vx, vy;   /* first derivative of vel v */
            float           uxy;      /* cross derivative of vel u */
            float           vxy;      /* cross derivative of vel v */
            float           uxx, uyy; /* second derivative of vel u */
            float           vxx, vyy; /* second derivative of vel v */
 	    float           dh[2];
 	    float           invs_dh[2];
 	    float           invsq_dh[2];
 	    float           cross_dh;
      
            dh[0] = nsten->wave->rect_grid->h[0];

	    invs_dh[0] = 1.0/(2.0*dh[0]);     /* inverse of 2.0*dh[0] */
    	
	    invsq_dh[0] = 1.0/(sqr(dh[0])); 

            dh[1] = nsten->wave->rect_grid->h[1];

	    invs_dh[1] = 1.0/(2.0*dh[1]);     /* inverse of 2.0*dh[1] */
	
	    invsq_dh[1] = 1.0/(sqr(dh[1]));

	    cross_dh = 1.0/(4.0*dh[0]*dh[1]); 


            if(debugging("ns_soln"))
	    {
	        (void) printf("sts2d[0][0] = %d\n",nsten->sts2d[0][0]);
                verbose_print_state("sts2d[0][0]",nsten->sts2d[0][0]);
	        (void) printf("sts2d[1][0] = %d\n",nsten->sts2d[1][0]);
                verbose_print_state("sts2d[1][0]",nsten->sts2d[1][0]);
	        (void) printf("sts2d[2][0] = %d\n",nsten->sts2d[2][0]);
	        (void) printf("icoords2d[2][0]=(%d,%d)\n",
			      nsten->icoords2d[2][0][0],
			      nsten->icoords2d[2][0][1]);
                verbose_print_state("sts2d[2][0]",nsten->sts2d[2][0]);
	        (void) printf("sts2d[0][1] = %d\n",nsten->sts2d[0][1]);
                verbose_print_state("sts2d[0][1]",nsten->sts2d[0][1]);
	        (void) printf("sts2d[1][1] = %d\n",nsten->sts2d[1][1]);
                verbose_print_state("sts2d[1][1]",nsten->sts2d[1][1]);
	        (void) printf("sts2d[2][1] = %d\n",nsten->sts2d[2][1]);
                verbose_print_state("sts2d[2][1]",nsten->sts2d[2][1]);
	        (void) printf("sts2d[0][2] = %d\n",nsten->sts2d[0][2]);
                verbose_print_state("sts2d[0][2]",nsten->sts2d[0][2]);
	        (void) printf("sts2d[1][2] = %d\n",nsten->sts2d[1][2]);
                verbose_print_state("sts2d[1][2]",nsten->sts2d[1][2]);
	        (void) printf("sts2d[2][2] = %d\n",nsten->sts2d[2][2]);
                verbose_print_state("sts2d[2][2]",nsten->sts2d[2][2]);
	    }
            /*** velocity_first_derivative ***/

            ux = (vel(0,nsten->sts2d[2][1]) - vel(0,nsten->sts2d[0][1]))*invs_dh[0];

            uy = (vel(0,nsten->sts2d[1][2]) - vel(0,nsten->sts2d[1][0]))*invs_dh[1];

            vx = (vel(1,nsten->sts2d[2][1]) - vel(1,nsten->sts2d[0][1]))*invs_dh[0];

            vy = (vel(1,nsten->sts2d[1][2]) - vel(1,nsten->sts2d[1][0]))*invs_dh[1];

	    if(debugging("ns_soln"))
                (void) printf("ux = %lf, uy = %lf, vx = %lf, vy = %lf\n",ux,uy,vx,vy); 
    
	    /*** velocity_cross_derivative ***/

            uxy = (vel(0,nsten->sts2d[2][2]) - vel(0,nsten->sts2d[2][0])
                  -vel(0,nsten->sts2d[0][2]) + vel(0,nsten->sts2d[0][0]))
                  *cross_dh;

            vxy = (vel(1,nsten->sts2d[2][2]) - vel(1,nsten->sts2d[2][0])
                  -vel(1,nsten->sts2d[0][2]) + vel(1,nsten->sts2d[0][0]))
                  *cross_dh;
	
	    if(debugging("ns_soln"))
                (void) printf("uxy = %lf, vxy = %lf\n",uxy,vxy); 
    
            /* uxy = vxy = 0.0;  Temporary setting */

	    /*** velocity_second_derivative ***/ 
            uxx = (vel(0,nsten->sts2d[2][1]) - 2.0*vel(0,nsten->sts2d[1][1])
                   + vel(0,nsten->sts2d[0][1]))*invsq_dh[0];

            uyy = (vel(0,nsten->sts2d[1][2]) - 2.0*vel(0,nsten->sts2d[1][1])
                   + vel(0,nsten->sts2d[1][0]))*invsq_dh[1];

            vxx = (vel(1,nsten->sts2d[2][1]) - 2.0*vel(1,nsten->sts2d[1][1])
                   + vel(1,nsten->sts2d[0][1]))*invsq_dh[0];

            vyy = (vel(1,nsten->sts2d[1][2]) - 2.0*vel(1,nsten->sts2d[1][1])
                   + vel(1,nsten->sts2d[1][0]))*invsq_dh[1];
	
	    if (debugging("ns_soln"))
	    {
                (void) printf("uxx = %lf, uyy = %lf, vxx = %lf, vyy = %lf\n",
				    uxx,uyy,vxx,vyy); 
	    }

            if (viscosity == YES)
            {
                if (!use_stokes_vis)
                {
                    Mom(ans)[0] += dt*mu*(uxx + uyy);
                    Mom(ans)[1] += dt*mu*(vxx + vyy);
                    Energy(ans) += dt*mu*(vel(0,nsten->sts2d[1][1])*(uxx + uyy)
                                    + vel(1,nsten->sts2d[1][1])*(vxx + vyy));
                }
                else
                {
                    Mom(ans)[0] += dt*mu*(uxx + uyy);
                    Mom(ans)[1] += dt*mu*(vxx + vyy);
                    Energy(ans) += dt*mu*((vel(0,nsten->sts2d[1][1])*(uxx + uyy))
                                    + (vel(1,nsten->sts2d[1][1])*(vxx + vyy))
                                    + sqr(vx + uy)
                                    + sqr(ux - vy) );
                }
            }

             if (thermal_conduction == YES)
                  Energy(ans) += dt*(kappa*Fourier_heat_conduction(nsten));

	    if(debugging("ns_soln"))
                (void) printf("Mom(ans)[0] adds %lf\n",
			    dt*mu*(2/3*(2*uxx-vxy) + (uyy+vxy)));  

	    if(debugging("ns_soln"))
                (void) printf("Mom(ans)[1] adds %lf\n",
			      dt*mu*(2/3*(2*vyy-uxy) + (vxx+uxy)));  

            if (mass_diffusion == YES)
            {
                float     hx,hx1,hx3,hy1,hy3;
                int       N = Params(ans)->n_comps;
                float     diffu_term_x[N], diffu_term_y[N];
                float     xf1[N], xf2[N], xf3[N], yf1[N], yf2[N], yf3[N];       
                float     dx1, dx2[N], dx3[N], dy1, dy2[N],dy3[N];
        
                float     x1 = Dens(nsten->sts2d[2][1]);
                float     x2 = Dens(nsten->sts2d[1][1]);
                float     x3 = Dens(nsten->sts2d[0][1]);
                float     y1 = Dens(nsten->sts2d[1][2]);
                float     y2 = Dens(nsten->sts2d[1][1]);
                float     y3 = Dens(nsten->sts2d[1][0]);

                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    int  i;
                    if(Params(ans)->n_comps != 1)
                    {
                         hx = specific_internal_energy_h_species(nsten->sts2d[1][1]) 
                              - specific_internal_energy_l_species(nsten->sts2d[1][1]);
                         hx1 = specific_internal_energy_h_species(nsten->sts2d[2][1])
                              - specific_internal_energy_l_species(nsten->sts2d[2][1]);
                         hx3 = specific_internal_energy_h_species(nsten->sts2d[0][1])
                              - specific_internal_energy_l_species(nsten->sts2d[0][1]);
                         hy1 = specific_internal_energy_h_species(nsten->sts2d[1][2])
                              - specific_internal_energy_l_species(nsten->sts2d[1][2]);
                         hy3 = specific_internal_energy_h_species(nsten->sts2d[1][0])
                              - specific_internal_energy_l_species(nsten->sts2d[1][0]);
                         
                         for(i = 0; i < Params(ans)->n_comps; i++)
                         {
                            xf1[i] = pdens(nsten->sts2d[2][1])[i]/x1;
                            xf2[i] = pdens(nsten->sts2d[1][1])[i]/x2;
                            xf3[i] = pdens(nsten->sts2d[0][1])[i]/x3;

                            yf1[i] = pdens(nsten->sts2d[1][2])[i]/y1;
                            yf2[i] = pdens(nsten->sts2d[1][1])[i]/y2;
                            yf3[i] = pdens(nsten->sts2d[1][0])[i]/y3;

                            dx1 = (x1 - x3)/(2*dh[0]);
                            dx2[i] = ((xf1[i] - xf3[i])/(2*dh[0]));
                            dx3[i] = ((xf1[i] - 2*xf2[i] + xf3[i])/(dh[0]*dh[0]));

                            dy1 = (y1 - y3)/(2*dh[1]);
                            dy2[i] = ((yf1[i] - yf3[i])/(2*dh[1]));
                            dy3[i] = ((yf1[i] - 2*yf2[i] + yf3[i])/(dh[1]*dh[1]));

                            diffu_term_x[i] = (dx1*dx2[i] + Dens(nsten->sts2d[1][1])*dx3[i]);
                            diffu_term_y[i] = (dy1*dy2[i] + Dens(nsten->sts2d[1][1])*dy3[i]);
                            pdens(ans)[i] += dt*d*(diffu_term_x[i] + diffu_term_y[i]);

                            if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                               pdens(ans)[i] = 0.0;
                            else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                            {
                               printf("ERROR in g_ns()\n");
                               printf("partial density %20.18g < 0.0\n",pdens(ans)[i]);
                               clean_up(ERROR);
                            }
                         }

                        float  sum = 0.0;
                        for(i = 0; i < Params(ans)->n_comps; i++)
                            sum += pdens(ans)[i];
                        for(i = 0; i < Params(ans)->n_comps; i++)
                            pdens(ans)[i] *= Dens(ans)/sum;

                        float  hdx1,hdx2,hdx3,hdy1,hdy2,hdy3;
                        float  hdiffu_term_x,hdiffu_term_y;

                        hdx1 = ((hx1*x1) - (x3*hx3))/(2*dh[0]);
                        hdx2 = ((xf1[1] - xf3[1])/(2*dh[0]));
                        hdx3 = ((xf1[1] - 2*xf2[1] + xf3[1])/(dh[0]*dh[0]));

                        hdy1 = ((hy1*y1) - (hy3*y3))/(2*dh[1]);
                        hdy2 = ((yf1[1] - yf3[1])/(2*dh[1]));
                        hdy3 = ((yf1[1] - 2*yf2[1] + yf3[1])/(dh[1]*dh[1]));

                        hdiffu_term_x = (hdx1*hdx2 + (hx)*Dens(nsten->sts3d[1][1][1])*hdx3);
                        hdiffu_term_y = (hdy1*hdy2 + (hx)*Dens(nsten->sts3d[1][1][1])*hdy3);
                        Energy(ans) += dt*d*(hdiffu_term_x + hdiffu_term_y);
                    }
                }
            }

#if defined(SUBGRID)
             float taux,tauy,qt,qp[2];
             float t1,t2,t3,t4,t5,t6,t7,t8,qt1,qt2,qt3,qt4;
             float qp1,qp2,qp3,qp4;
             float qp5,qp6,qp7,qp8;

             t1 = Tau(nsten->sts2d[0][1])[0][0];
             t2 = Tau(nsten->sts2d[2][1])[0][0];
             t3 = Tau(nsten->sts2d[1][0])[0][1];
             t4 = Tau(nsten->sts2d[1][2])[0][1];
             t5 = Tau(nsten->sts2d[0][1])[1][0];
             t6 = Tau(nsten->sts2d[2][1])[1][0];
             t7 = Tau(nsten->sts2d[1][0])[1][1];
             t8 = Tau(nsten->sts2d[1][2])[1][1];
             qt1 = Qh(nsten->sts2d[0][1])[0];
             qt2 = Qh(nsten->sts2d[2][1])[0];
             qt3 = Qh(nsten->sts2d[1][0])[1];
             qt4 = Qh(nsten->sts2d[1][2])[1];
             qp1 = Qc(nsten->sts2d[0][1])[0];
             qp2 = Qc(nsten->sts2d[2][1])[0];
             qp3 = Qc(nsten->sts2d[1][0])[1];
             qp4 = Qc(nsten->sts2d[1][2])[1];
             qp5 = Qc(nsten->sts2d[0][1])[2];
             qp6 = Qc(nsten->sts2d[2][1])[2];
             qp7 = Qc(nsten->sts2d[1][0])[3];
             qp8 = Qc(nsten->sts2d[1][2])[3];

             taux = ((t2 - t1)/(2*dh[0]))+((t4 - t3)/(2*dh[1]));
             tauy = ((t6 - t5)/(2*dh[0]))+((t8 - t7)/(2*dh[1]));
             qt = ((qt2 - qt1)/(2*dh[0]))+((qt4 - qt3)/(2*dh[1]));
             qp[1] = ((qp2 - qp1)/(2*dh[0]))+((qp4 - qp3)/(2*dh[1]));
             qp[0] = ((qp6 - qp5)/(2*dh[0]))+((qp8 - qp7)/(2*dh[1]));

             float taukkr, taukkl;
             float taukkx, taukky;
             taukkr = t2 + t8;
             taukkl = t1 + t7;
             taukkx = ((taukkr*vel(0,nsten->sts2d[2][1]))
                      - (taukkl*vel(0,nsten->sts2d[0][1]))) / (2*dh[0]);
             taukky = ((taukkr*vel(1,nsten->sts2d[1][2]))
                      - (taukkl*vel(1,nsten->sts2d[1][0]))) / (2*dh[1]);

            if(subgrid_vis == YES)
            {
                 Mom(ans)[0] -= dt*taux;
                 Mom(ans)[1] -= dt*tauy;
                 Energy(ans) += dt*0.5*(taukkx + taukky);
            }

            if(subgrid_con == YES)
                 Energy(ans) -= dt*qt;

             if(subgrid_md == YES)
             {
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    int  i;
                    if(Params(ans)->n_comps != 1)
                    {
                        for(i = 0; i < Params(ans)->n_comps; i++)
                        {
                               pdens(ans)[i] -= dt*qp[i];
                               if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                                  pdens(ans)[i] = 0.0;
                               else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                               {
                                  printf("ERROR in g_ns() after subgrid\n");
                                  printf("partial density %20.18g < 0.0\n",pdens(ans)[i]);
                                  clean_up(ERROR);
                               }
                        }

                        float sum = 0.0;
                        for(i = 0; i < Params(ans)->n_comps; i++)
                           sum += pdens(ans)[i];
                        for(i = 0; i < Params(ans)->n_comps; i++)
                           pdens(ans)[i] *= Dens(ans)/sum;

                    }
                }
            }
#endif /* defined SUBGRID */

  	    if(debugging("ns_soln"))
	    {
        	(void) printf("txx + tyy = %lf\n",
			      Fourier_heat_conduction(nsten));

        	(void) printf("Energy(ans) adds %lf\n",
		         dt * (kappa*Fourier_heat_conduction(nsten)
                         + mu*(2/3*(ux*(2*ux - vy) + vy*(2*vy - ux))
                         + sqr(uy+vx)
                         + vel(0,nsten->sts2d[1][1])*(2/3*(2*uxx - vxy)
                         + (uyy + vxy))
                         + vel(1,nsten->sts2d[1][1])*(2/3*(2*vyy - uxy)
                         + (vxx + uxy))
                         ))); 
                verbose_print_state("sts2d[0][0]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[1][0]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[2][0]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[0][1]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[1][1]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[2][1]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[0][2]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[1][2]",nsten->sts2d[1][1]); 
                verbose_print_state("sts2d[2][2]",nsten->sts2d[1][1]); 
                verbose_print_state("NS_ans",ans); 
	    }
	    /* Add 1/r terms for cylindrical geometry */
	    /* : d/d(phi) = 0, all phi-component = 0  */
	    if (gr->Remap.remap == CYLINDRICAL_REMAP)
	    {
	      	float r = coords[0];
	      	float u = vel(0,nsten->sts2d[1][1]);
              	float v = vel(1,nsten->sts2d[1][1]);
	      	Mom(ans)[0] += dt*mu*2*twothird/r*(ux - u/r);
	      	Mom(ans)[1] += dt*mu/r*(vx + uy/3);
              	Energy(ans) += dt*mu*(4*(u*ux/r  -  u*u/(r*r))/3
                                    + v*uy/(3*r) + v*vx/r);
	    }
	}
	break;
	
#endif /* defined TWOD */
#if defined THREED
	case 3:
	{
	    float           rhoxx, rhoyy, rhozz;
	    float           ux,uxx,uyy,uzz,uxy,uxz;
	    float           vy,vxx,vyy,vzz,vyx,vyz;
	    float           wz,wxx,wyy,wzz,wzx,wzy;
            float           dux[3], duy[3], duz[3], Fx, Fy, Fz;
	    float           dh[3];
	    int             i;
	    
	    dh[0] = nsten->wave->rect_grid->h[0];
	    dh[1] = nsten->wave->rect_grid->h[1];
	    dh[2] = nsten->wave->rect_grid->h[2];
	    
	    rhoxx = (Dens(nsten->sts3d[2][1][1]) - 2.0*Dens(nsten->sts3d[1][1][1])
	               + Dens(nsten->sts3d[0][1][1]))/(dh[0]*dh[0]);

            rhoyy = (Dens(nsten->sts3d[1][2][1]) - 2.0*Dens(nsten->sts3d[1][1][1])
		       + Dens(nsten->sts3d[1][0][1]))/(dh[1]*dh[1]);

	    rhozz = (Dens(nsten->sts3d[1][1][2]) - 2.0*Dens(nsten->sts3d[1][1][1])
		       + Dens(nsten->sts3d[1][1][0]))/(dh[2]*dh[2]);
           
	     for(i = 0; i < 3; i++)
             {
                 dux[i] = (vel(i,nsten->sts3d[2][1][1]) 
		          -  vel(i,nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                 duy[i] = (vel(i,nsten->sts3d[1][2][1]) 
		          -  vel(i,nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                 duz[i] = (vel(i,nsten->sts3d[1][1][2]) 
		          -  vel(i,nsten->sts3d[1][1][0])) / (2.0*dh[2]);
             }

	    uxx = (vel(0,nsten->sts3d[2][1][1]) - 2.0*vel(0,nsten->sts3d[1][1][1])
			      + vel(0,nsten->sts3d[0][1][1]))/(dh[0]*dh[0]);

	    uyy = (vel(0,nsten->sts3d[1][2][1]) - 2.0*vel(0,nsten->sts3d[1][1][1])
			       + vel(0,nsten->sts3d[1][0][1]))/(dh[1]*dh[1]);

	    uzz = (vel(0,nsten->sts3d[1][1][2]) - 2.0*vel(0,nsten->sts3d[1][1][1])
			      + vel(0,nsten->sts3d[1][1][0]))/(dh[2]*dh[2]);
	    
	    vxx = (vel(1,nsten->sts3d[2][1][1]) - 2.0*vel(1,nsten->sts3d[1][1][1])
			           + vel(1,nsten->sts3d[0][1][1]))/(dh[0]*dh[0]);

	    vyy = (vel(1,nsten->sts3d[1][2][1]) - 2.0*vel(1,nsten->sts3d[1][1][1])
			           + vel(1,nsten->sts3d[1][0][1]))/(dh[1]*dh[1]);

	    vzz = (vel(1,nsten->sts3d[1][1][2]) - 2.0*vel(1,nsten->sts3d[1][1][1])
			        + vel(1,nsten->sts3d[1][1][0]))/(dh[2]*dh[2]);

	    wxx = (vel(2,nsten->sts3d[2][1][1]) - 2.0*vel(2,nsten->sts3d[1][1][1])
			          + vel(2,nsten->sts3d[0][1][1]))/(dh[0]*dh[0]);

	    wyy = (vel(2,nsten->sts3d[1][2][1]) - 2.0*vel(2,nsten->sts3d[1][1][1])
			          + vel(2,nsten->sts3d[1][0][1]))/(dh[1]*dh[1]);

	    wzz = (vel(2,nsten->sts3d[1][1][2]) - 2.0*vel(2,nsten->sts3d[1][1][1])
			          + vel(2,nsten->sts3d[1][1][0]))/(dh[2]*dh[2]);
	    
	    //Stokes
	    uxy = ((vel(0,nsten->sts3d[2][2][1]) - vel(0,nsten->sts3d[0][2][1])) 
	         -(vel(0,nsten->sts3d[2][0][1]) - vel(0,nsten->sts3d[0][0][1])))/(4.0*dh[0]*dh[1]);
 	    
	    uxz = ((vel(0,nsten->sts3d[2][1][2]) - vel(0,nsten->sts3d[0][1][2])) 
	         -(vel(0,nsten->sts3d[2][1][0]) - vel(0,nsten->sts3d[0][1][0])))/(4.0*dh[0]*dh[2]);
 	    
 	    vyx = ((vel(1,nsten->sts3d[2][2][1]) - vel(1,nsten->sts3d[0][2][1])) 
	         -(vel(1,nsten->sts3d[2][0][1]) - vel(1,nsten->sts3d[0][0][1])))/(4.0*dh[0]*dh[1]);
    
	    vyz = ((vel(1,nsten->sts3d[1][2][2]) - vel(1,nsten->sts3d[1][0][2])) 
	         -(vel(1,nsten->sts3d[1][2][0]) - vel(1,nsten->sts3d[1][0][0])))/(4.0*dh[1]*dh[2]);
    	    
	    wzy = ((vel(2,nsten->sts3d[1][2][2]) - vel(2,nsten->sts3d[1][0][2])) 
	         -(vel(2,nsten->sts3d[1][2][0]) - vel(2,nsten->sts3d[1][0][0])))/(4.0*dh[1]*dh[2]);
 	    
	    wzx = ((vel(2,nsten->sts3d[2][1][2]) - vel(2,nsten->sts3d[0][1][2])) 
	         -(vel(2,nsten->sts3d[2][1][0]) - vel(2,nsten->sts3d[0][1][0])))/(4.0*dh[0]*dh[2]);

            if (viscosity == YES)
            {
                if (!use_stokes_vis)
                {
                    //Incompressible viscosity term
                    Mom(ans)[0] += dt*mu*(uxx + uyy + uzz);
                    Mom(ans)[1] += dt*mu*(vxx + vyy + vzz);
                    Mom(ans)[2] += dt*mu*(wxx + wyy + wzz);

                    Energy(ans) += dt*mu*(vel(0,nsten->sts3d[1][1][1])*(uxx + uyy + uzz)
                                 + vel(1,nsten->sts3d[1][1][1])*(vxx + vyy + vzz)
                                 + vel(2,nsten->sts3d[1][1][1])*(wxx + wyy + wzz));
                }
                else
                {
                    Fx = (twothird*(2.0*uxx-vyx-wzx) + (uyy+vyx) + (uzz+wzx));
                    Fy = (twothird*(2.0*vyy-uxy-wzy) + (uxy+vxx) + (vzz+wzy));
                    Fz = (twothird*(2.0*wzz-uxz-vyz) + (uxz+wxx) + (vyz+wyy));
                    Mom(ans)[0] += dt*mu*Fx;
                    Mom(ans)[1] += dt*mu*Fy;
                    Mom(ans)[2] += dt*mu*Fz;
                    //Energy = Kinetic terms + dissipation terms
                    VIS_FORCE = dt*(mu*(
                            vel(0,nsten->sts3d[1][1][1])*Fx
                            + vel(1,nsten->sts3d[1][1][1])*Fy
                            + vel(2,nsten->sts3d[1][1][1])*Fz
                       + 2.0*( sqr(dux[0]) + sqr(duy[1]) + sqr(duz[2]) )
                       + sqr(dux[1] + duy[0])
                       + sqr(duy[2] + duz[1])
                       + sqr(dux[2] + duz[0])
                       - twothird*sqr(dux[0] + duy[1] + duz[2])));
                    Energy(ans) += dt*(mu*(
                            vel(0,nsten->sts3d[1][1][1])*Fx
                            + vel(1,nsten->sts3d[1][1][1])*Fy
                            + vel(2,nsten->sts3d[1][1][1])*Fz
                       + 2.0*( sqr(dux[0]) + sqr(duy[1]) + sqr(duz[2]) )
                       + sqr(dux[1] + duy[0])
                       + sqr(duy[2] + duz[1])
                       + sqr(dux[2] + duz[0])
                       - twothird*sqr(dux[0] + duy[1] + duz[2])));
                     //fprintf(stdout, "viscosity %e %e %e %e %e %e \n",dt*mu*Fx,dt*mu*Fy,dt*mu*Fz,dt,mu,VIS_FORCE);        
                }
            }
            if (thermal_conduction == YES)
            {
                Energy(ans) += dt*(kappa*Fourier_heat_conduction(nsten));
            }

            if (mass_diffusion == YES)
	    {
                float     hx,hx1,hx3,hy1,hy3,hz1,hz3;
                int       N = Params(ans)->n_comps;
                float     diffu_term_x[N], diffu_term_y[N], diffu_term_z[N];
                float  xf1[N], xf2[N], xf3[N], yf1[N], yf2[N], yf3[N], zf1[N], zf2[N], zf3[N];
                float  dx1, dx2[N], dx3[N], dy1, dy2[N],dy3[N], dz1, dz2[N],dz3[N];

                float     x1 = Dens(nsten->sts3d[2][1][1]);
                float     x2 = Dens(nsten->sts3d[1][1][1]);
                float     x3 = Dens(nsten->sts3d[0][1][1]);
                float     y1 = Dens(nsten->sts3d[1][2][1]);
                float     y2 = Dens(nsten->sts3d[1][1][1]);
                float     y3 = Dens(nsten->sts3d[1][0][1]);
                float     z1 = Dens(nsten->sts3d[1][1][2]);
                float     z2 = Dens(nsten->sts3d[1][1][1]);
                float     z3 = Dens(nsten->sts3d[1][1][0]);


                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    int  i;
                    if(Params(ans)->n_comps != 1)
                    {
                         hx = specific_internal_energy_h_species(nsten->sts3d[1][1][1])                              - specific_internal_energy_l_species(nsten->sts3d[1][1][1]);

                         hx1 = specific_internal_energy_h_species(nsten->sts3d[2][1][1])
                               - specific_internal_energy_l_species(nsten->sts3d[2][1][1]);

                         hx3 = specific_internal_energy_h_species(nsten->sts3d[0][1][1]) 
                               - specific_internal_energy_l_species(nsten->sts3d[0][1][1]);

                         hy1 = specific_internal_energy_h_species(nsten->sts3d[1][2][1])
                               - specific_internal_energy_l_species(nsten->sts3d[1][2][1]);

                         hy3 = specific_internal_energy_h_species(nsten->sts3d[1][0][1])
                               - specific_internal_energy_l_species(nsten->sts3d[1][0][1]);

                         hz1 = specific_internal_energy_h_species(nsten->sts3d[1][1][2])
                               - specific_internal_energy_l_species(nsten->sts3d[1][1][2]);

                         hz3 = specific_internal_energy_h_species(nsten->sts3d[1][1][0])
                               - specific_internal_energy_l_species(nsten->sts3d[1][1][0]);

                        for(i = 0; i < Params(ans)->n_comps; i++)
                        {
                            xf1[i] = pdens(nsten->sts3d[2][1][1])[i]/x1;
                            xf2[i] = pdens(nsten->sts3d[1][1][1])[i]/x2;
                            xf3[i] = pdens(nsten->sts3d[0][1][1])[i]/x3;

                            yf1[i] = pdens(nsten->sts3d[1][2][1])[i]/y1;
                            yf2[i] = pdens(nsten->sts3d[1][1][1])[i]/y2;
                            yf3[i] = pdens(nsten->sts3d[1][0][1])[i]/y3;

                            zf1[i] = pdens(nsten->sts3d[1][1][2])[i]/z1;
                            zf2[i] = pdens(nsten->sts3d[1][1][1])[i]/z2;
                            zf3[i] = pdens(nsten->sts3d[1][1][0])[i]/z3;

                            dx1 = (x1 - x3)/(2*dh[0]);
                            dx2[i] = ((xf1[i] - xf3[i])/(2*dh[0]));
                            dx3[i] = ((xf1[i] - 2*xf2[i] + xf3[i])/(dh[0]*dh[0]));

                            dy1 = (y1 - y3)/(2*dh[1]);
                            dy2[i] = ((yf1[i] - yf3[i])/(2*dh[1]));
                            dy3[i] = ((yf1[i] - 2*yf2[i] + yf3[i])/(dh[1]*dh[1]));

                            dz1 = (z1 - z3)/(2*dh[2]);
                            dz2[i] = ((zf1[i] - zf3[i])/(2*dh[2]));
                            dz3[i] = ((zf1[i] - 2*zf2[i] + zf3[i])/(dh[2]*dh[2]));

                            /*
                            diffu_term_x[i] = (dx1*dx2[i] + Dens(nsten->sts3d[1][1][1])*dx3[i]);
                            diffu_term_y[i] = (dy1*dy2[i] + Dens(nsten->sts3d[1][1][1])*dy3[i]);
                            diffu_term_z[i] = (dz1*dz2[i] + Dens(nsten->sts3d[1][1][1])*dz3[i]);
                            */

                            //stable alg for diffusion calc
                            diffu_term_x[i] = 0.5*(x2+x3)/x3*pdens(nsten->sts3d[0][1][1])[i]
                                             +0.5*(x1+x2)/x1*pdens(nsten->sts3d[2][1][1])[i]
                                             -2.0*pdens(nsten->sts3d[1][1][1])[i];
                            diffu_term_x[i] /= dh[0]*dh[0];

                            diffu_term_y[i] = 0.5*(y2+y3)/y3*pdens(nsten->sts3d[1][0][1])[i]
                                             +0.5*(y1+y2)/y1*pdens(nsten->sts3d[1][2][1])[i]
                                             -2.0*pdens(nsten->sts3d[1][1][1])[i];
                            diffu_term_y[i] /= dh[1]*dh[1];

                            diffu_term_z[i] = 0.5*(z2+z3)/z3*pdens(nsten->sts3d[1][1][0])[i]
                                             +0.5*(z1+z2)/z1*pdens(nsten->sts3d[1][1][2])[i]
                                             -2.0*pdens(nsten->sts3d[1][1][1])[i];
                            diffu_term_z[i] /= dh[2]*dh[2];

                            pdens(ans)[i] += dt*d*(diffu_term_x[i] + diffu_term_y[i] + diffu_term_z[i]);
                            //fprintf(stdout, "partial density mass diffusion %e %e %e \n",dt*d*(diffu_term_x[i] + diffu_term_y[i] + diffu_term_z[i]),dt,d);
                            if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                               pdens(ans)[i] = 0.0;
                            else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                            {
                               printf("ERROR in g_ns()\n");
                               printf("partial density %20.18g < 0.0\n",pdens(ans)[i]);
                               clean_up(ERROR);
                            }
                        }

                        float  sum = 0.0;
                        for(i = 0; i < Params(ans)->n_comps; i++)
                            sum += pdens(ans)[i];
                        for(i = 0; i < Params(ans)->n_comps; i++)
                            pdens(ans)[i] *= Dens(ans)/sum;

                        float  hdx1,hdx2,hdx3,hdy1,hdy2,hdy3,hdz1,hdz2,hdz3;
                        float  hdiffu_term_x,hdiffu_term_y,hdiffu_term_z;

                        hdx1 = ((hx1*x1) - (x3*hx3))/(2*dh[0]);
                        hdx2 = ((xf1[1] - xf3[1])/(2*dh[0]));
                        hdx3 = ((xf1[1] - 2*xf2[1] + xf3[1])/(dh[0]*dh[0]));

                        hdy1 = ((hy1*y1) - (hy3*y3))/(2*dh[1]);
                        hdy2 = ((yf1[1] - yf3[1])/(2*dh[1]));
                        hdy3 = ((yf1[1] - 2*yf2[1] + yf3[1])/(dh[1]*dh[1]));

                        hdz1 = ((hz1*z1) - (hz3*z3))/(2*dh[2]);
                        hdz2 = ((zf1[1] - zf3[1])/(2*dh[2]));
                        hdz3 = ((zf1[1] - 2*zf2[1] + zf3[1])/(dh[2]*dh[2]));

                        hdiffu_term_x = (hdx1*hdx2 + (hx)*Dens(nsten->sts3d[1][1][1])*hdx3);
                        hdiffu_term_y = (hdy1*hdy2 + (hx)*Dens(nsten->sts3d[1][1][1])*hdy3);
                        hdiffu_term_z = (hdz1*hdz2 + (hx)*Dens(nsten->sts3d[1][1][1])*hdz3);
                        Energy(ans) += dt*d*(hdiffu_term_x + hdiffu_term_y + hdiffu_term_z);
                        //fprintf(stdout, "Energy mass diffusion %e %e %e \n",dt*d*(hdiffu_term_x + hdiffu_term_y + hdiffu_term_z),dt,d);
                    }
                }
            }

#if defined(SUBGRID)
            float taux,tauy,tauz,qt,qp[2];
            float t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18;
            float qt1,qt2,qt3,qt4,qt5,qt6;
            float qp1,qp2,qp3,qp4,qp5,qp6,qp7,qp8,qp9,qp10,qp11,qp12;

             t1 = Tau(nsten->sts3d[0][1][1])[0][0];
             t2 = Tau(nsten->sts3d[2][1][1])[0][0];
             t3 = Tau(nsten->sts3d[1][0][1])[0][1];
             t4 = Tau(nsten->sts3d[1][2][1])[0][1];
             t5 = Tau(nsten->sts3d[1][1][0])[0][2];
             t6 = Tau(nsten->sts3d[1][1][2])[0][2];

             t7 = Tau(nsten->sts3d[0][1][1])[1][0];
             t8 = Tau(nsten->sts3d[2][1][1])[1][0];
             t9 = Tau(nsten->sts3d[1][0][1])[1][1];
             t10 = Tau(nsten->sts3d[1][2][1])[1][1];
             t11 = Tau(nsten->sts3d[1][1][0])[1][2];
             t12 = Tau(nsten->sts3d[1][1][2])[1][2];

             t13 = Tau(nsten->sts3d[0][1][1])[2][0];
             t14 = Tau(nsten->sts3d[2][1][1])[2][0];
             t15 = Tau(nsten->sts3d[1][0][1])[2][1];
             t16 = Tau(nsten->sts3d[1][2][1])[2][1];
             t17 = Tau(nsten->sts3d[1][1][0])[2][2];
             t18 = Tau(nsten->sts3d[1][1][2])[2][2];

             qt1 = Qh(nsten->sts3d[0][1][1])[0];
             qt2 = Qh(nsten->sts3d[2][1][1])[0];
             qt3 = Qh(nsten->sts3d[1][0][1])[1];
             qt4 = Qh(nsten->sts3d[1][2][1])[1];
             qt5 = Qh(nsten->sts3d[1][1][0])[2];
             qt6 = Qh(nsten->sts3d[1][1][2])[2];

             qp1 = Qc(nsten->sts3d[0][1][1])[0];
             qp2 = Qc(nsten->sts3d[2][1][1])[0];
             qp3 = Qc(nsten->sts3d[1][0][1])[1];
             qp4 = Qc(nsten->sts3d[1][2][1])[1];
             qp5 = Qc(nsten->sts3d[1][1][0])[2];
             qp6 = Qc(nsten->sts3d[1][1][2])[2];

             qp7 = Qc(nsten->sts3d[0][1][1])[3];
             qp8 = Qc(nsten->sts3d[2][1][1])[3];
             qp9 = Qc(nsten->sts3d[1][0][1])[4];
             qp10 = Qc(nsten->sts3d[1][2][1])[4];
             qp11 = Qc(nsten->sts3d[1][1][0])[5];
             qp12 = Qc(nsten->sts3d[1][1][2])[5];

             taux = ((t2 - t1)/(2*dh[0]))+((t4 - t3)/(2*dh[1]))+((t6 - t5)/(2*dh[2]));
             tauy = ((t8 - t7)/(2*dh[0]))+((t10 - t9)/(2*dh[1]))+((t12 - t11)/(2*dh[2]));
             tauz = ((t14 - t13)/(2*dh[0]))+((t16 - t15)/(2*dh[1]))+((t18 - t17)/(2*dh[2]));
             qt = ((qt2 - qt1)/(2*dh[0]))+((qt4 - qt3)/(2*dh[1]))+((qt6 - qt5)/(2*dh[2]));
             qp[1] = ((qp2 - qp1)/(2*dh[0]))+((qp4 - qp3)/(2*dh[1]))+((qp6 - qp5)/(2*dh[2]));
             qp[0] = ((qp8 - qp7)/(2*dh[0]))+((qp10 - qp9)/(2*dh[1]))+((qp12 - qp11)/(2*dh[2]));

             float taukkr, taukkl;
             float taukkx, taukky, taukkz;
             taukkr = t2 + t10 + t18;
             taukkl = t1 + t9 + t17;
             taukkx = ((taukkr*vel(0,nsten->sts3d[2][1][1]))
                      - (taukkl*vel(0,nsten->sts3d[0][1][1]))) / (2*dh[0]);
             taukky = ((taukkr*vel(1,nsten->sts3d[1][2][1]))
                      - (taukkl*vel(1,nsten->sts3d[1][0][1]))) / (2*dh[1]);
             taukkz = ((taukkr*vel(2,nsten->sts3d[1][1][2]))
                      - (taukkl*vel(2,nsten->sts3d[1][1][0]))) / (2*dh[2]);

            if(subgrid_vis == YES)
            {
                 Mom(ans)[0] -= dt*taux;
                 Mom(ans)[1] -= dt*tauy;
                 Mom(ans)[2] -= dt*tauz;
                 Energy(ans) += dt*0.5*(taukkx + taukky + taukkz);
                 //fprintf(stdout, "subgrid_vis %e %e %e %e %e \n",-dt*taux,-dt*tauy,-dt*tauz,dt*0.5*(taukkx + taukky + taukkz),dt);

            }

            if(subgrid_con == YES)
            { 
                 Energy(ans) -= dt*qt;
            }

             if(subgrid_md == YES)
             {           
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    int  i;
                    if(Params(ans)->n_comps != 1)
                    {
                        for(i = 0; i < Params(ans)->n_comps; i++)
                        {
                               pdens(ans)[i] -= dt*qp[i];
                               //fprintf(stdout, "subgrid partial density mass diffusion %e %e \n",dt*qp[i],dt);
                               if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                                  pdens(ans)[i] = 0.0;
                               else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                               {
                                  printf("ERROR in g_ns() after subgrid\n");
                                  printf("partial density %20.18g < 0.0\n",pdens(ans)[i]);
                                  clean_up(ERROR);
                               }
                        }

                        float sum = 0.0;
                        for(i = 0; i < Params(ans)->n_comps; i++)
                           sum += pdens(ans)[i];
                        for(i = 0; i < Params(ans)->n_comps; i++)
                           pdens(ans)[i] *= Dens(ans)/sum;
                            
                    }
                }
            }
#endif /* defined SUBGRID */
	}
	break;
	 
#endif /* defined THREED */
	}
        /*printf("Leaving ns_soln()\n");*/
}


/*
 *      Fourier heat conduction is used in the energy equation in Navier-
 *	Stokes equations, and in fact is a Laplacian of temperature T. 
 *
*/

LOCAL float Fourier_heat_conduction(
	    Pstencil         *nsten)
{
	RECT_GRID       *gr = nsten->wave->rect_grid;
	int             dim = gr->dim;
	
        switch (dim)
        {
#if defined ONED
        case 1:
        { 
            float           txx;   /* second derivative of temp t */
	    float           dh[1];
	
	dh[0] = nsten->wave->rect_grid->h[0];

        txx = (temperature(nsten->sts1d[2]) - 2.0*temperature(nsten->sts1d[1])
               + temperature(nsten->sts1d[0]))/(dh[0]*dh[0]);

        return txx;
	}
	break;
         
#endif /* defined ONED */
#if defined TWOD
        case 2:
	{
            float           txx, tyy;   /* second derivative of temp t */
	    float           dh[2];

	dh[0] = nsten->wave->rect_grid->h[0];
	dh[1] = nsten->wave->rect_grid->h[1];

        txx = (temperature(nsten->sts2d[2][1]) 
               - 2.0*temperature(nsten->sts2d[1][1])
               + temperature(nsten->sts2d[0][1]))/(dh[0]*dh[0]);
	
        tyy = (temperature(nsten->sts2d[1][2]) 
               - 2.0*temperature(nsten->sts2d[1][1])
               + temperature(nsten->sts2d[1][0]))/(dh[1]*dh[1]);

        return (txx + tyy);
	}
	break;
         
#endif /* defined TWOD */
#if defined THREED
        case 3:
	{
	    float           txx, tyy, tzz;   /* second derivative of temp t */
	    float           dh[3];

	dh[0] = nsten->wave->rect_grid->h[0];
	dh[1] = nsten->wave->rect_grid->h[1];
	dh[2] = nsten->wave->rect_grid->h[2];

        txx = (temperature(nsten->sts3d[1][1][2]) 
               - 2.0*temperature(nsten->sts3d[1][1][1]) 
	       + temperature(nsten->sts3d[1][1][0]))/(dh[0]*dh[0]);

	tyy = (temperature(nsten->sts3d[1][2][1])
               - 2.0*temperature(nsten->sts3d[1][1][1]) 
	       + temperature(nsten->sts3d[1][0][1]))/(dh[1]*dh[1]);

	tzz = (temperature(nsten->sts3d[2][1][1])
               - 2.0*temperature(nsten->sts3d[1][1][1]) 
	       + temperature(nsten->sts3d[0][1][1]))/(dh[2]*dh[2]);

 	return (txx + tyy + tzz);
	}
	break;
         
#endif /* defined THREED */
        }
}

LOCAL  bool  g_compute_NS_terms(
	      Locstate        state) 
{
  	if ( (Params(state) == NULL) ||
       	    (Params(state)->eos->_compute_ns_terms == NO))
            return NO;
  	else if (Params(state)->eos->_compute_ns_terms == YES)
       	    return YES;
  	return NO;
}

        /* Compute the diffused mass flux */
LOCAL   float diffusive_mass_flux(
	float           density_jump,
	float           mix_len_new,
	float           mix_len_old)
{
        float           flux;
        flux = density_jump*(mix_len_new - mix_len_old)/2.0;
	return flux;

	
}
/* End diffusive_mass_flux */

	/* Use only in parabolic step */
LOCAL	int g_neumann_bdry_state_beta(
	float		*coords,
	COMPONENT	int_comp,
	POINT		*pt,
	HYPER_SURF	*Nbdry,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	Wave		*wave = (Wave*)p2wave;
	float		nor[MAXD];
	float		coords_ref[MAXD];
	float		coords_tmp[MAXD];

	int		i, dim = front->rect_grid->dim;

	float           tmp_vel[MAXD], vn;

	/* coords_ref is the location of the reflection of the point coords
	   across the Neumann boundary. */

	if (debugging("no_nbdry") 
	    || !reflect_pt_about_Nbdry(coords,coords_ref,nor,
					  int_comp,Nbdry,front))
	    return NO;

	//print_general_vector("#coord_ref=", coords_ref, 3, "\n");
	/* Get the state at coords_ref. */
	hyp_solution(coords_ref,int_comp,Nbdry,UNKNOWN_SIDE,
		     front,wave,state,NULL);

	set_state(state,TGAS_STATE,state);
	vn = 0.0;
	for (i = 0; i < dim; i ++)
	{
	    tmp_vel[i] = Vel(state)[i];
	    vn += tmp_vel[i] * nor[i];
	}

	if (no_slip(Nbdry))
        {
	    float alpha = adherence_coeff(Nbdry);
            for (i = 0; i < dim; i ++)
	    	Vel(state)[i] = (1-2*alpha)*Vel(state)[i];
	}

	zero_normal_velocity(state,nor,dim);
	for (i = 0; i < dim; i ++)
	    Vel(state)[i] += - vn * nor[i];
	set_state(state,GAS_STATE,state);
	return YES;
}		/*end g_neumann_bdry_state_beta*/


LOCAL	void find_ghost_state(
	int *ic,
	int *icn,
	Locstate state,
	Front *front,
	Wave *wv)
{
	GRID_DIRECTION dir_l[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION dir_u[3] = {EAST,NORTH,UPPER};
	GRID_DIRECTION dir[MAXD];
	int i,j,k,nc,num_dic;
    	int ii,ic_tmp[MAXD];
	int dim = wv->rect_grid->dim;
	CRXING *crx[6];
	COMPONENT comp = Rect_comp(ic,wv);
	size_t sizest = wv->sizest;
	Locstate s1,s2,s3;
	static TRI *stat_tri;
	float p1[MAXD],p2[MAXD],p[MAXD];
	float ps[MAXD],pe[MAXD];
	float *pt1,*pt2,*pt3;
	float alpha;

	num_dic = 0;
	for (i = 0; i < dim; ++i)
	{
	    if (ic[i] != icn[i])
	    {
	    	num_dic++;
		dir[i] = (icn[i] - ic[i] == -1) ? dir_l[i] : dir_u[i];
	    }
	}
	switch (num_dic)
	{
	case 1:		/* Apply to 1D, 2D, 3D */
	    for (i = 0; i < dim; ++i)
	    	if (ic[i] != icn[i])
	    	    crx[0] = Rect_crossing(ic,dir[i],wv);
	    if (crx[0] == NULL)
	    {
	    	screen("ERROR: (%d %d %d) and (%d %d %d) ",
			ic[0],ic[1],ic[2],icn[0],icn[1],icn[2]);
		screen("are connected with different comps\n");
		clean_up(ERROR);
	    }
	    s1 = state_with_comp(crx[0]->pt,crx[0]->hs,comp);
	    ft_assign(state,s1,sizest);
	    break;
	case 2:		/* Apply to 2D, 3D */
	    nc = 0;
	    for (i = 0; i < dim; ++i)
	    {
	    	if (ic[i] != icn[i])
		{
		    crx[nc] = Rect_crossing(ic,dir[i],wv);
		    if (crx[nc] != NULL) ++nc;
		    else
		    {
		    	for (k = 0; k < dim; ++k) ic_tmp[k] = ic[k];
			ic_tmp[i] = icn[i];
			for (ii = 0; ii < dim; ++ii)
			{
			    if (ic_tmp[ii] != icn[ii])
			    {
		    		crx[nc] = Rect_crossing(ic_tmp,dir[ii],wv);
			    	if (crx[nc] != NULL) ++nc;
			    	else
			    	{
			    	    screen("ERROR: (%d %d %d) and (%d %d %d) ",
					ic[0],ic[1],ic[2],icn[0],icn[1],icn[2]);
				    screen("are connected but "
				    	"with different comps\n");
				    clean_up(ERROR);
			    	}
			    }
			}
		    }
		}
	    }

	    s1 = state_with_comp(crx[0]->pt,crx[0]->hs,comp);
	    s2 = state_with_comp(crx[1]->pt,crx[1]->hs,comp);
	    /* Map it to two dimensional plane */
	    if (dim == 2) ii = 0;
	    else
	    {
	    	for (ii = 0; ii < dim; ++ii) 
		    if (ic[ii] == icn[ii]) break;
	    	++ii;
	    }
	    for (i = 0; i < 2; ++i)
	    {
	    	p1[i] = Coords(crx[0]->pt)[(i+ii)%dim];
	    	p2[i] = Coords(crx[1]->pt)[(i+ii)%dim];
	    	ps[i] = Rect_coords(ic,wv)[(i+ii)%dim];
	    	pe[i] = Rect_coords(icn,wv)[(i+ii)%dim];
	    }
	    if (!cross_line_segments(ps,pe,p1,p2,p))
	    {
	    	screen("ERROR: no crossing found\n");
		clean_up(ERROR);
	    }
	    if (fabs(p1[0] - p2[0]) > fabs(p1[1] - p2[1]))
	    	alpha = fabs(p[0] - p2[0])/fabs(p1[0] - p2[0]);
	    else
	    	alpha = fabs(p[1] - p2[1])/fabs(p1[1] - p2[1]);
	    bi_interpolate_intfc_states(front->interf,alpha,1.0-alpha,
	           Coords(crx[0]->pt),s1,Coords(crx[1]->pt),s2,state);
	    break;
	case 3:		/* Apply to 3D */
	    nc = 0;
	    for (i = 0; i < dim; ++i)
	    {
	    	if (ic[i] != icn[i])
		{
		    crx[nc] = Rect_crossing(ic,dir[i],wv);
		    if (crx[nc] != NULL)
		    	++nc;
		    else
		    {
		    	for (k = 0; k < dim; ++k) ic_tmp[k] = ic[k];
			ic_tmp[i] = icn[i];
		    	crx[nc] = Rect_crossing(ic_tmp,dir[(i+1)%dim],wv);
			if (crx[nc] != NULL) ++nc;
			else
			{
			    ic_tmp[(i+1)%dim] = icn[(i+1)%dim];
		    	    crx[nc] = Rect_crossing(ic_tmp,dir[(i+2)%dim],wv);
			    if (crx[nc] != NULL) ++nc;
			    else
			    {
			    	screen("ERROR: (%d %d %d) and (%d %d %d) ",
					ic[0],ic[1],ic[2],icn[0],icn[1],icn[2]);
				screen("are connected with different comps\n");
				clean_up(ERROR);
			    }
			}
		    	for (k = 0; k < dim; ++k) ic_tmp[k] = ic[k];
			ic_tmp[i] = icn[i];
		    	crx[nc] = Rect_crossing(ic_tmp,dir[(i+2)%dim],wv);
			if (crx[nc] != NULL) ++nc;
			else
			{
			    ic_tmp[(i+2)%dim] = icn[(i+2)%dim];
		    	    crx[nc] = Rect_crossing(ic_tmp,dir[(i+1)%dim],wv);
			    if (crx[nc] != NULL) ++nc;
			    else
			    {
			    	screen("ERROR: (%d %d %d) and (%d %d %d) ",
					ic[0],ic[1],ic[2],icn[0],icn[1],icn[2]);
				screen("are connected with different comps\n");
				clean_up(ERROR);
			    }
			}
		    }
		}
	    }
	    for (i = 0; i < dim; ++i)
	    {
	    	ps[i] = Rect_coords(ic,wv)[i];
	    	pe[i] = Rect_coords(icn,wv)[i];
	    }
	    /* Order and eliminate possibly identical crossings */
	    float d1,d2;
	    for (i = 0; i < nc-1; ++i)
	    {
		d1 = distance_between_positions(ps,Coords(crx[i]->pt),3);
	    	for (ii = i+1; ii < nc; ++ii)
		{
		    if (crx[i] == crx[ii])
		    {
		    	for (k = ii+1; k < nc; ++k)
			    crx[k-1] = crx[k];
		    	--nc;
		    }
		    d2 = distance_between_positions(ps,Coords(crx[ii]->pt),3);
		    if (d1 < d2)
		    {
		    	float dtmp;
			CRXING *crx_tmp;
			dtmp = d1;
			d1 = d2;
			d2 = dtmp;
			crx_tmp = crx[i];
			crx[i] = crx[ii];
			crx[ii] = crx_tmp;

		    }
		}
	    }
	    for (i = 0; i < nc; i++)
	    {
		pt1 = Coords(crx[i]->pt);
	    	for (j = i+1; j < nc; ++j)
		{
		    pt2 = Coords(crx[j]->pt);
		    for (k = j+1; k < nc; ++k)
		    {
		    	pt3 = Coords(crx[k]->pt);
			if (cross_line_segment_triangle(pt1,pt2,pt3,ps,pe,p))
			{
	    		    float f[3];
	    		    s1 = state_with_comp(crx[i]->pt,crx[i]->hs,comp);
	    		    s2 = state_with_comp(crx[j]->pt,crx[j]->hs,comp);
	    		    s3 = state_with_comp(crx[k]->pt,crx[k]->hs,comp);
	    		    linear_interp_coefs_three_pts(f,p,pt1,pt2,pt3);
	    		    tri_interpolate_intfc_states(front->interf,f[0],
			    		f[1],f[2],pt1,s1,pt2,s2,pt3,s3,state);
			    return;
			}
		    }
		}
	    }
	    screen("ERROR: in num_dic == 3 case, cannot find ghost state\n");
	    clean_up(ERROR);
	}
}	/* end find_ghost_state */


LOCAL bool cross_line_segments(
	float		*p1s,
	float		*p1e,
	float		*p2s,
	float		*p2e,
	float		*p)
{
	double		x1=(double)p1s[0];
	double		y1=(double)p1s[1];
	double		x2=(double)p1e[0];
	double		y2=(double)p1e[1];
	double		x3=(double)p2s[0];
	double		y3=(double)p2s[1];
	double		x4=(double)p2e[0];
	double		y4=(double)p2e[1];
	double		nor_dist_t,nor_dist_s;
	double		sinth;		/* sin of angle between bonds */
	double		xcross,ycross;	/* coord of intersection
					   of lines 12 34 */
	double		x00,y00;	/* beginning of long bond */
	double		x0,y0;		/* beginning of short bond after
					   coord translation */
	double		x,y;		/* end of long bond after
					   coord translation */
	double		dx,dy;		/* short bond end - start */
	double		t;		/* fractional distance on short bond */
	double		s;		/* fractional distance on long bond */
	double		len12;
	double		parallel = PARALLEL(current_interface());

	double d1 = distance_between_positions(p1s,p1e,2);
	double d2 = distance_between_positions(p2s,p2e,2);
	if (d1 > d2) 
	{
	    x00 = x1;		y00 = y1;
	    x0 = x3 - x1;	y0 = y3 - y1;
	    x  = x2 - x1;	y  = y2 - y1;
	    dx = x4 - x3;	dy = y4 - y3;
	}
	else 
	{
	    x00 = x3;		y00 = y3;
	    x0 = x1 - x3;	y0 = y1 - y3;
	    x  = x4 - x3;	y  = y4 - y3;
	    dx = x2 - x1;	dy = y2 - y1;
	}
	sinth = dx*y - dy*x;
	nor_dist_t = x0*y - y0*x;
	nor_dist_s = dx*y0 - dy*x0;
	len12 = d1*d2;

	if (fabs(sinth) <= parallel * len12) 
	{
	    /* Case of parallel lines */
	    if (fabs(nor_dist_t) <= parallel * len12) 
	    {
	    	/* Lines coincide */
	    	if (Between(x0,0.0,x) && Between(y0,0.0,y)) 
	    	{
	    	    /* Cross at x0,y0 */
	    	    p[0] = (float)(x0 + x00);
	    	    p[1] = (float)(y0 + y00);
	    	    return YES;
	    	}
		if (Between(x0+dx,0.0,x) && Between(y0+dy,0.0,y)) 
		{
		    /* Cross at x0+dx,y0+dy */
		    p[0] = (float)(x0 + dx + x00);
		    p[1] = (float)(y0 + dy + y00);
		    return YES;
		}
		return NO; /* No cross; line segments don't overlap */
	    }
	    return NO; /* No cross; lines distinct although parallel */
	}

		/* Now lines are not parallel */

	t = - nor_dist_t / sinth;
	s = nor_dist_s / sinth;
	if (t <= 0.0 || t > 1.0 || s <= 0.0 || s > 1.0)
	    return NO;
	xcross = 0.5*(x0 + t*dx + s*x);
	ycross = 0.5*(y0 + t*dy + s*y);
	p[0] = (float)(xcross + x00);
	p[1] = (float)(ycross + y00);
	return YES;
}		/*end cross_line_segments */

LOCAL bool cross_line_segment_triangle(
        float           *p1t1,
        float           *p1t2,
        float           *p1t3,
        float           *p2s,
        float           *p2e,
        float           *p)
{
        int             i, max_index;
        bool            intersect_with_plane;
        double          v1[3],v2[3],v3[3],v4[2],v5[2],v6[2];
        double          det1,det2;
        double          tmp1[3],tmp2[3],tmp3[3];
        double          p1t1_tmp[2],p1t2_tmp[2],p1t3_tmp[2],p_tmp[2];
        double          n1,n2,n3,n4,t,nv1[3],nv2[3],norm[3];


        for(i = 0; i < 3; i++)
        {
            v1[i] = p2s[i] - p1t1[i];
            v2[i] = p1t2[i] - p1t1[i];
            v3[i] = p1t3[i] - p1t1[i];
        }

        det1 = Det3d(v1,v2,v3);

        for(i = 0; i < 3; i++)
            v1[i] = p2e[i] - p1t1[i];

        det2 = Det3d(v1,v2,v3);

        if(det1 == 0 && det2 == 0)
            intersect_with_plane = NO;
        else if(det1*det2 > 0)
            intersect_with_plane = NO;
        else
            intersect_with_plane = YES;

        if(intersect_with_plane != YES)
            return NO;

        tmp1[0] = p1t1[1];  tmp1[1] = p1t1[2];  tmp1[2] = 1.0;
        tmp2[0] = p1t2[1];  tmp2[1] = p1t2[2];  tmp2[2] = 1.0;
        tmp3[0] = p1t3[1];  tmp3[1] = p1t3[2];  tmp3[2] = 1.0;
        n1 = Det3d(tmp1,tmp2,tmp3);

        tmp1[0] = p1t1[0];  tmp1[1] = p1t1[2];  tmp1[2] = 1.0;
        tmp2[0] = p1t2[0];  tmp2[1] = p1t2[2];  tmp2[2] = 1.0;
        tmp3[0] = p1t3[0];  tmp3[1] = p1t3[2];  tmp3[2] = 1.0;
        n2 = Det3d(tmp1,tmp2,tmp3);

        tmp1[0] = p1t1[0];  tmp1[1] = p1t1[1];  tmp1[2] = 1.0;
        tmp2[0] = p1t2[0];  tmp2[1] = p1t2[1];  tmp2[2] = 1.0;
        tmp3[0] = p1t3[0];  tmp3[1] = p1t3[1];  tmp3[2] = 1.0;
        n3 = Det3d(tmp1,tmp2,tmp3);

        n4 = Det3d(p1t1,p1t2,p1t3);

        t = (p2s[1]*n2-p2s[0]*n1-p2s[2]*n3+n4)
            /((p2e[0]-p2s[0])*n1-(p2e[1]-p2s[1])*n2+(p2e[2]-p2s[2])*n3);

        p[0] = p2s[0] + (p2e[0] - p2s[0])*t;
        p[1] = p2s[1] + (p2e[1] - p2s[1])*t;
        p[2] = p2s[2] + (p2e[2] - p2s[2])*t;

        for(i = 0; i < 3; i++)
        {
            nv1[i] = p1t2[i] - p1t1[i];
            nv2[i] = p1t3[i] - p1t1[i];
        }

        Cross3d(nv1,nv2,norm);

        if(within_tri(p,p1t1,p1t2,p1t3,norm,1.0e-14) != YES)
            return NO;

        return YES;
}

#if defined(SUBGRID)
void   SGS(
        float   dt,
        Front   *front,
        Wave    *wave)
{
        int             i, dim = front->interf->dim;
        int             *iperm;
        iperm = set_iperm(front->step,dim);

        start_clock("SGS");
        switch (dim)
        {
#if defined ONED
        case 1:
        {
            //SGS1d(dt,front,wave);
        }
        break;

#endif /* defined ONED */
#if defined TWOD
        case 2:
        {
            SGS2d(dt,front,wave);
        }
        break;

#endif /* defined TWOD */
#if defined THREED
        case 3:
        {
            SGS3d(dt,front,wave);
        }
        break;

#endif /* defined THREED */
        }
        for (i = 0; i < dim; i++)
        {
            if (!scatter_states(wave,front,iperm,i))
            {
                screen("scatter_states() failed in direction\n",i);
                clean_up(ERROR);
            }
        }
        stop_clock("SGS");
}

LOCAL  void   SGS2d(
        float           dt,
        Front           *front,
        Wave            *wv)
{
        RECT_GRID       *rgr = wv->rect_grid;
        float           *L = rgr->L, *U = rgr->U;
        Gas_param       **params;
        INTERFACE       *intfc = front->interf;
        int             dim = rgr->dim;
        int             *gmax = rgr->gmax;
        int             i,j,k;
        int             icoords[MAXD];
        static          Pstencil *nsten;
        int             npts,mid;
        int             imin[MAXD],imax[MAXD];
        size_t          sizest = front->sizest;
        bool            viscosity, mass_diffusion;
        float           ux,uy,uz,vx,vy,vz,wx,wy,wz;
        float           tx,ty,cx,cy,cx0,cy0,S11,S12,S22,SS;
        float           hat_ux,hat_uy,hat_vx,hat_vy,hat_tx,hat_ty,hat_cx,hat_cy,hat_cx0,hat_cy0;
        float           hat_S11,hat_S12,hat_S22,hat_SS;
        float           s11, s12, s21, s22, dtemx, dtemy, pdenx, pdeny, pdenx0, pdeny0,S;
        const int       nn = pp_numnodes();
        int             myid = pp_mynode();
        int             *ppgmax = front->pp_grid->gmax;
        int             ppx = myid % ppgmax[0];
        int             ppy = (myid-ppx)/ppgmax[0];
        int             mink, maxk;

        for (i = 0; i < dim; i++)
        {
            imin[i] = 0;        imax[i] = gmax[i];
        }

        if(ppy == ppgmax[1]-1)
        {
            maxk = imax[1];
            if(ppy == 0)
                mink = imin[1];
            else
                mink = imin[1]-2;
        }
        else if(ppy == 0)
        {
            mink = imin[1];
            if(ppy == ppgmax[1]-1)
                maxk = imax[1];
            else
                maxk = imax[1]+2;
        }
        else
        {
            mink = imin[1]-2;
            maxk = imax[1]+2;
        }

        if (nsten == NULL)
        {
            stat_scalar(&nsten,sizeof(Pstencil));
            switch (dim)
            {
            case 1:
                bi_array(&nsten->icoords1d,3,1,INT);
                uni_array(&nsten->sts1d,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                    alloc_state(front->interf,&nsten->sts1d[i],front->sizest);
                break;
            case 2:
                tri_array(&nsten->icoords2d,3,3,2,INT);
                bi_array(&nsten->sts2d,3,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                {
                    for (j = 0; j < 3; j++)
                    {
                        alloc_state(front->interf,&nsten->sts2d[i][j],
                                front->sizest);
                    }
                }
                break;
            case 3:
                quad_array(&nsten->icoords3d,3,3,3,3,INT);
                tri_array(&nsten->sts3d,3,3,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        for (k = 0; k < 3; k++)
                            alloc_state(front->interf,&nsten->sts3d[i][j][k],
                                     front->sizest);
                break;
            }
        }
        nsten->fr = front;
        nsten->wave = wv;
        nsten->npts = npts = wv->npts_sten;
        mid = (int) (npts/2);

        Locstate        state;
        float           **C, **co_C, **CI, **co_CI, **Prt, **co_Prt, **Sct, **co_Sct,**Sct0, **co_Sct0;
        float           **co_coords_x, **co_coords_y;
        float           **co_rho, **co_vel_u, **co_vel_v, **co_rho_vel_u, **co_rho_vel_v,**co_rho_vel_uv,**co_rho_vel_uu,**co_rho_vel_vv;
        float           **co_rho_cp, **co_cp, **co_rho_t, **co_rho_c, **co_rho_c0,**co_t, **co_con, **co_con0;
        float           **co_rho_cu, **co_rho_cv, **co_rho_cu0, **co_rho_cv0,**co_rho_cp_tu, **co_rho_cp_tv;
        float           **co_temp11, **co_temp12, **co_temp21, **co_temp22, **co_temp, **co_temp_t1, **co_temp_t2, **co_temp_c1, **co_temp_c2, **co_temp_c10, **co_temp_c20;
        float           **MA11, **MA12, **MA21, **MA22, **MI, **LA11, **LA12, **LA22, **LA1, **LA2, **LI, **MH1, **MH2, **LH1, **LH2, **MC1, **MC2, **LC1, **LC2, **MC10, **MC20, **LC10, **LC20;
        float           **r, *cs_new, *ci_new,*prt_new,*sct_new,*sct0_new, *cs_new_new, *ci_new_new,*prt_new_new,*sct_new_new,*sct0_new_new;
        int             **r_color, **r_color_LA,**r_color_LH,**r_color_LC,**r_color_LC0,**r_color_LI;
        float           **re_coords_x, **re_coords_y;
        float           **re_rho, **re_vel_u, **re_vel_v, **re_rho_vel_u, **re_rho_vel_v,**re_rho_vel_uv,**re_rho_vel_uu,**re_rho_vel_vv;
        float           **re_rho_cp, **re_cp, **re_rho_t, **re_rho_c, **re_rho_c0,**re_t, **re_c, **re_c0;
        float           **re_rho_cu, **re_rho_cv, **re_rho_cu0, **re_rho_cv0, **re_rho_cp_tu, **re_rho_cp_tv;
        float           **re_temp11, **re_temp12, **re_temp21, **re_temp22, **re_temp, **re_temp_t1, **re_temp_t2, **re_temp_c1, **re_temp_c2, **re_temp_c10, **re_temp_c20;
        bool            **re_NS_terms, **co_NS_terms;
        float           *deno, *nume;
        float           *ci_deno, *ci_nume;
        float           *prt_deno, *prt_nume;
        float           *sct_deno, *sct_nume;
        float           *sct_deno0, *sct_nume0;
        int             ii, jj, iii, jjj;
        float           dh[2];
        dh[0] = nsten->wave->rect_grid->h[0];
        dh[1] = nsten->wave->rect_grid->h[1];
        float           delta = sqrt(dh[0]*dh[1]);
        float           co_delta = sqrt((2*dh[0])*(2*dh[1]));
        float           delta_r = 2*dh[0];
        int             num_r = (int)((U[0]-L[0])/delta_r)+1;

        switch (dim)
        {
        case 2:
            bi_array(&co_coords_x,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_coords_y,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_cp,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_cp,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_t,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_c,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
	    bi_array(&co_rho_c0,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_t,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_con,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_con0,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_cu,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_cv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_cu0,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_cv0,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));

            bi_array(&co_rho_cp_tu,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_cp_tv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));

            bi_array(&co_vel_u,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_vel_v,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_vel_u,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_vel_v,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_vel_uu,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_vel_vv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_rho_vel_uv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp11,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp12,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp21,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp22,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp_t1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp_t2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp_c1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp_c2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp_c10,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&co_temp_c20,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));

            bi_array(&MA11,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MA12,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MA21,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MA22,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MI,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LA11,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LA12,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LA22,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LA1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LA2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LI,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MH1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MH2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LH1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LH2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MC1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MC2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LC1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LC2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MC10,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&MC20,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LC10,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&LC20,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&r,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
	    bi_array(&r_color,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(int));
            bi_array(&r_color_LA,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(int));
            bi_array(&r_color_LH,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(int));
            bi_array(&r_color_LC,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(int));
            bi_array(&r_color_LC0,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(int));
            bi_array(&r_color_LI,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(int));
            bi_array(&co_NS_terms,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(bool));
            bi_array(&re_coords_x,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_coords_y,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_cp,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_cp,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_t,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_c,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_c0,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_t,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_c,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_c0,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_cu,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_cv,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_cu0,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_cv0,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_cp_tu,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_cp_tv,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_vel_u,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_vel_v,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_vel_u,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_vel_v,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_vel_uu,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_vel_vv,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_rho_vel_uv,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp11,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp12,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp21,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp22,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp_t1,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp_t2,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp_c1,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp_c2,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp_c10,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_temp_c20,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(float));
            bi_array(&re_NS_terms,imax[0]-imin[0]+4,imax[1]-imin[1]+4,sizeof(bool));
            bi_array(&C,imax[0]-imin[0],imax[1]-imin[1],sizeof(float));
            bi_array(&co_C,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&CI,imax[0]-imin[0],imax[1]-imin[1],sizeof(float));
            bi_array(&co_CI,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&Prt,imax[0]-imin[0],imax[1]-imin[1],sizeof(float));
            bi_array(&co_Prt,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&Sct,imax[0]-imin[0],imax[1]-imin[1],sizeof(float));
            bi_array(&co_Sct,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
            bi_array(&Sct0,imax[0]-imin[0],imax[1]-imin[1],sizeof(float));
            bi_array(&co_Sct0,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,sizeof(float));
        }

        for (j = mink; j < maxk; j++)
        {
            icoords[1] = j;
            for (i = imin[0]-2; i < imax[0]+2; i++)
            {
                icoords[0] = i;
                state =  Rect_state(icoords,wv);
                if(g_compute_NS_terms(Rect_state(icoords,wv)))
                {
                    fill_2d_9pt_Pstencil(i,j,nsten);

                    ux = (vel(0,nsten->sts2d[2][1]) -
                                vel(0,nsten->sts2d[0][1])) / (2.0*dh[0]);
                    uy = (vel(0,nsten->sts2d[1][2]) -
                                vel(0,nsten->sts2d[1][0])) / (2.0*dh[1]);
                    vx = (vel(1,nsten->sts2d[2][1]) -
                                vel(1,nsten->sts2d[0][1])) / (2.0*dh[0]);
                    vy = (vel(1,nsten->sts2d[1][2]) -
                                vel(1,nsten->sts2d[1][0])) / (2.0*dh[1]);
                    tx = (temperature(nsten->sts2d[2][1]) -
                          temperature(nsten->sts2d[0][1])) / (2.0*dh[0]);
                    ty = (temperature(nsten->sts2d[1][2]) -
                          temperature(nsten->sts2d[1][0])) / (2.0*dh[1]);
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    cx = ((pdens(nsten->sts2d[2][1])[1]/Dens(nsten->sts2d[2][1])) -
                          (pdens(nsten->sts2d[0][1])[1]/Dens(nsten->sts2d[0][1]))) / (2.0*dh[0]);
                    cy = ((pdens(nsten->sts2d[1][2])[1]/Dens(nsten->sts2d[1][2])) -
                          (pdens(nsten->sts2d[1][0])[1]/Dens(nsten->sts2d[1][0]))) / (2.0*dh[1]);
                    cx0 = ((pdens(nsten->sts2d[2][1])[0]/Dens(nsten->sts2d[2][1])) -
                           (pdens(nsten->sts2d[0][1])[0]/Dens(nsten->sts2d[0][1]))) / (2.0*dh[0]);
                    cy0 = ((pdens(nsten->sts2d[1][2])[0]/Dens(nsten->sts2d[1][2])) -
                           (pdens(nsten->sts2d[1][0])[0]/Dens(nsten->sts2d[1][0]))) / (2.0*dh[1]);
                    }
                    re_coords_x[i+2][j+2] = Rect_coords(icoords,wv)[0];
                    re_coords_y[i+2][j+2] = Rect_coords(icoords,wv)[1];
                    re_rho[i+2][j+2] = Dens(state);
                    re_cp[i+2][j+2] = C_P(state);
                    re_t[i+2][j+2] = temperature(state);
                    re_vel_u[i+2][j+2] = vel(0,state);
                    re_vel_v[i+2][j+2] = vel(1,state);
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    re_c[i+2][j+2] = pdens(state)[1]/Dens(state);
                    re_c0[i+2][j+2] = pdens(state)[0]/Dens(state);
                    re_rho_cu[i+2][j+2] = Dens(state)*(pdens(state)[1]/Dens(state))*vel(0,state);
                    re_rho_cv[i+2][j+2] = Dens(state)*(pdens(state)[1]/Dens(state))*vel(1,state);
                    re_rho_cu0[i+2][j+2] = Dens(state)*(pdens(state)[0]/Dens(state))*vel(0,state);
                    re_rho_cv0[i+2][j+2] = Dens(state)*(pdens(state)[0]/Dens(state))*vel(1,state);
                    re_rho_c[i+2][j+2] = Dens(state)*(pdens(state)[1]/Dens(state));
                    re_rho_c0[i+2][j+2] = Dens(state)*(pdens(state)[0]/Dens(state));

                    }
                    re_rho_cp_tu[i+2][j+2] = Dens(state)*C_P(state)*temperature(state)*vel(0,state);
                    re_rho_cp_tv[i+2][j+2] = Dens(state)*C_P(state)*temperature(state)*vel(1,state);
                    re_rho_cp[i+2][j+2] = Dens(state)*C_P(state);
                    re_rho_t[i+2][j+2] = Dens(state)*temperature(state);
                    re_rho_vel_u[i+2][j+2] = Dens(state)*vel(0,state);
                    re_rho_vel_v[i+2][j+2] = Dens(state)*vel(1,state);
                    re_rho_vel_uu[i+2][j+2] = Dens(state)*vel(0,state)*vel(0,state);
                    re_rho_vel_vv[i+2][j+2] = Dens(state)*vel(1,state)*vel(1,state);
                    re_rho_vel_uv[i+2][j+2] = Dens(state)*vel(0,state)*vel(1,state);
                    S11 = ux;
                    S12 = (uy + vx)/2.0;
                    S22 = vy;
                    SS = sqrt(2.0*( (S11*S11) + 2.0*(S12*S12) + (S22*S22)));

                    re_temp11[i+2][j+2] = Dens(state)*SS*((S11/2.0)-(S22/2.0));
                    re_temp12[i+2][j+2] = Dens(state)*SS*S12;
                    re_temp21[i+2][j+2] = Dens(state)*SS*S12;
                    re_temp22[i+2][j+2] = Dens(state)*SS*((S22/2.0)-(S11/2.0));
                    re_temp[i+2][j+2] = Dens(state)*SS*SS;
                    re_temp_t1[i+2][j+2] = Dens(state)*C_P(state)*SS*tx;
                    re_temp_t2[i+2][j+2] = Dens(state)*C_P(state)*SS*ty;
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    re_temp_c1[i+2][j+2] = Dens(state)*SS*cx;
                    re_temp_c2[i+2][j+2] = Dens(state)*SS*cy;
                    re_temp_c10[i+2][j+2] = Dens(state)*SS*cx0;
                    re_temp_c20[i+2][j+2] = Dens(state)*SS*cy0;
                    }
                    re_NS_terms[i+2][j+2] = YES;
                }
                else
                {
                    re_NS_terms[i+2][j+2] = NO;
                }
            }
        }

        for (j = imin[1]; j < (imax[1]/2)+2; j++)
        {
            jj = 2*j;
            for (i = imin[0]; i < (imax[0]/2)+2; i++)
            {
                ii = 2*i;
                if(re_NS_terms[ii][jj] == YES &&
                   re_NS_terms[ii+1][jj] == YES &&
                   re_NS_terms[ii][jj+1] == YES &&
                   re_NS_terms[ii+1][jj+1] == YES)
                {
                    co_coords_x[i][j] = re_coords_x[ii][jj] + (dh[0]/2);
                    co_coords_y[i][j] = re_coords_y[ii][jj] + (dh[1]/2);
                    co_rho[i][j] = (re_rho[ii][jj] + re_rho[ii+1][jj]
                                   + re_rho[ii][jj+1] + re_rho[ii+1][jj+1])/4.0;
                    co_vel_u[i][j] = (re_vel_u[ii][jj] + re_vel_u[ii+1][jj]
                                   + re_vel_u[ii][jj+1] + re_vel_u[ii+1][jj+1])/4.0;
                    co_vel_v[i][j] = (re_vel_v[ii][jj] + re_vel_v[ii+1][jj]
                                   + re_vel_v[ii][jj+1] + re_vel_v[ii+1][jj+1])/4.0;

                    co_rho_vel_u[i][j] = (re_rho_vel_u[ii][jj] + re_rho_vel_u[ii+1][jj]
                                          + re_rho_vel_u[ii][jj+1] + re_rho_vel_u[ii+1][jj+1])/4.0;
                    co_rho_vel_v[i][j] = (re_rho_vel_v[ii][jj] + re_rho_vel_v[ii+1][jj]
                                          + re_rho_vel_v[ii][jj+1] + re_rho_vel_v[ii+1][jj+1])/4.0;
                    co_rho_vel_uu[i][j] = (re_rho_vel_uu[ii][jj] + re_rho_vel_uu[ii+1][jj]
                                          + re_rho_vel_uu[ii][jj+1] + re_rho_vel_uu[ii+1][jj+1])/4.0;
                    co_rho_vel_vv[i][j] = (re_rho_vel_vv[ii][jj] + re_rho_vel_vv[ii+1][jj]
                                          + re_rho_vel_vv[ii][jj+1] + re_rho_vel_vv[ii+1][jj+1])/4.0;
                    co_rho_vel_uv[i][j] = (re_rho_vel_uv[ii][jj] + re_rho_vel_uv[ii+1][jj]
                                          + re_rho_vel_uv[ii][jj+1] + re_rho_vel_uv[ii+1][jj+1])/4.0;
                    co_rho_cp[i][j] = (re_rho_cp[ii][jj] + re_rho_cp[ii+1][jj]
                                          + re_rho_cp[ii][jj+1] + re_rho_cp[ii+1][jj+1])/4.0;
                    co_cp[i][j] = (re_cp[ii][jj] + re_cp[ii+1][jj]
                                          + re_cp[ii][jj+1] + re_cp[ii+1][jj+1])/4.0;
                    co_rho_t[i][j] = (re_rho_t[ii][jj] + re_rho_t[ii+1][jj]
                                          + re_rho_t[ii][jj+1] + re_rho_t[ii+1][jj+1])/4.0;
                    co_rho_c[i][j] = (re_rho_c[ii][jj] + re_rho_c[ii+1][jj]
                                          + re_rho_c[ii][jj+1] + re_rho_c[ii+1][jj+1])/4.0;
                    co_rho_c0[i][j] = (re_rho_c0[ii][jj] + re_rho_c0[ii+1][jj]
                                          + re_rho_c0[ii][jj+1] + re_rho_c0[ii+1][jj+1])/4.0;

                    co_t[i][j] = (re_t[ii][jj] + re_t[ii+1][jj]
                                          + re_t[ii][jj+1] + re_t[ii+1][jj+1])/4.0;
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    co_con[i][j] = (re_c[ii][jj] + re_c[ii+1][jj]
                                          + re_c[ii][jj+1] + re_c[ii+1][jj+1])/4.0;
                    co_con0[i][j] = (re_c0[ii][jj] + re_c0[ii+1][jj]
                                          + re_c0[ii][jj+1] + re_c0[ii+1][jj+1])/4.0;
                    co_rho_cu[i][j] = (re_rho_cu[ii][jj] + re_rho_cu[ii+1][jj]
                                          + re_rho_cu[ii][jj+1] + re_rho_cu[ii+1][jj+1])/4.0;
                    co_rho_cv[i][j] = (re_rho_cv[ii][jj] + re_rho_cv[ii+1][jj]
                                          + re_rho_cv[ii][jj+1] + re_rho_cv[ii+1][jj+1])/4.0;
                    co_rho_cu0[i][j] = (re_rho_cu0[ii][jj] + re_rho_cu0[ii+1][jj]
                                          + re_rho_cu0[ii][jj+1] + re_rho_cu0[ii+1][jj+1])/4.0;
                    co_rho_cv0[i][j] = (re_rho_cv0[ii][jj] + re_rho_cv0[ii+1][jj]
                                          + re_rho_cv0[ii][jj+1] + re_rho_cv0[ii+1][jj+1])/4.0;
                    }
                    co_rho_cp_tu[i][j] = (re_rho_cp_tu[ii][jj] + re_rho_cp_tu[ii+1][jj]
                                          + re_rho_cp_tu[ii][jj+1] + re_rho_cp_tu[ii+1][jj+1])/4.0;
                    co_rho_cp_tv[i][j] = (re_rho_cp_tv[ii][jj] + re_rho_cp_tv[ii+1][jj]
                                          + re_rho_cp_tv[ii][jj+1] + re_rho_cp_tv[ii+1][jj+1])/4.0;

                    co_temp11[i][j] = (re_temp11[ii][jj] + re_temp11[ii+1][jj]
                                       + re_temp11[ii][jj+1] + re_temp11[ii+1][jj+1])/4.0;
                    co_temp12[i][j] = (re_temp12[ii][jj] + re_temp12[ii+1][jj]
                                       + re_temp12[ii][jj+1] + re_temp12[ii+1][jj+1])/4.0;
                    co_temp21[i][j] = (re_temp21[ii][jj] + re_temp21[ii+1][jj]
                                       + re_temp21[ii][jj+1] + re_temp21[ii+1][jj+1])/4.0;
                    co_temp22[i][j] = (re_temp22[ii][jj] + re_temp22[ii+1][jj]
                                       + re_temp22[ii][jj+1] + re_temp22[ii+1][jj+1])/4.0;

                    co_temp[i][j] = (re_temp[ii][jj] + re_temp[ii+1][jj]
                                       + re_temp[ii][jj+1] + re_temp[ii+1][jj+1])/4.0;

                    co_temp_t1[i][j] = (re_temp_t1[ii][jj] + re_temp_t1[ii+1][jj]
                                       + re_temp_t1[ii][jj+1] + re_temp_t1[ii+1][jj+1])/4.0;
                    co_temp_t2[i][j] = (re_temp_t2[ii][jj] + re_temp_t2[ii+1][jj]
                                       + re_temp_t2[ii][jj+1] + re_temp_t2[ii+1][jj+1])/4.0;
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    co_temp_c1[i][j] = (re_temp_c1[ii][jj] + re_temp_c1[ii+1][jj]
                                       + re_temp_c1[ii][jj+1] + re_temp_c1[ii+1][jj+1])/4.0;
                    co_temp_c2[i][j] = (re_temp_c2[ii][jj] + re_temp_c2[ii+1][jj]
                                       + re_temp_c2[ii][jj+1] + re_temp_c2[ii+1][jj+1])/4.0;
                    co_temp_c10[i][j] = (re_temp_c10[ii][jj] + re_temp_c10[ii+1][jj]
                                       + re_temp_c10[ii][jj+1] + re_temp_c10[ii+1][jj+1])/4.0;
                    co_temp_c20[i][j] = (re_temp_c20[ii][jj] + re_temp_c20[ii+1][jj]
                                       + re_temp_c20[ii][jj+1] + re_temp_c20[ii+1][jj+1])/4.0;
                    }
                    co_NS_terms[i][j] = YES;
                }
                else
                    co_NS_terms[i][j] = NO;
            }
        }

        for (j = imin[1]+1; j < (imax[1]/2)+1; j++)
        {
            for (i = imin[0]+1; i < (imax[0]/2)+1; i++)
            {
                if(co_NS_terms[i][j] == YES)
                {
                    hat_ux = ((co_rho_vel_u[i+1][j]/co_rho[i+1][j])
                                  - (co_rho_vel_u[i-1][j]/co_rho[i-1][j])) / (4.0*dh[0]);
                    hat_uy = ((co_rho_vel_u[i][j+1]/co_rho[i][j+1])
                                  - (co_rho_vel_u[i][j-1]/co_rho[i][j-1])) / (4.0*dh[1]);
                    hat_vx = ((co_rho_vel_v[i+1][j]/co_rho[i+1][j])
                                  - (co_rho_vel_v[i-1][j]/co_rho[i-1][j])) / (4.0*dh[0]);
                    hat_vy = ((co_rho_vel_v[i][j+1]/co_rho[i][j+1])
                                  - (co_rho_vel_v[i][j-1]/co_rho[i][j-1])) / (4.0*dh[1]);

                    hat_tx  = (co_t[i+1][j] - co_t[i-1][j]) / (4.0*dh[0]);
                    hat_ty  = (co_t[i][j+1] - co_t[i][j-1]) / (4.0*dh[1]);
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    hat_cx  = (co_con[i+1][j] - co_con[i-1][j]) / (4.0*dh[0]);
                    hat_cy  = (co_con[i][j+1] - co_con[i][j-1]) / (4.0*dh[1]);
                    hat_cx0  = (co_con0[i+1][j] - co_con0[i-1][j]) / (4.0*dh[0]);
                    hat_cy0  = (co_con0[i][j+1] - co_con0[i][j-1]) / (4.0*dh[1]);
                    }
                    hat_S11 = hat_ux;
                    hat_S12 = (hat_uy + hat_vx)/2.0;
                    hat_S22 = hat_vy;
                    hat_SS = sqrt(2.0*( (hat_S11*hat_S11) + 2.0*(hat_S12*hat_S12)
                                            + (hat_S22*hat_S22)));

                    MA11[i][j] = 2.0*delta*delta*co_temp11[i][j]
                                 - (2*co_delta*co_delta*co_rho[i][j]*hat_SS
                                 *((hat_S11/2.0)-(hat_S22/2.0)));
                    MA12[i][j] = 2.0*delta*delta*co_temp12[i][j]
                                 - (2.0*co_delta*co_delta*co_rho[i][j]*hat_SS
                                   *hat_S12);

                    LA11[i][j] = co_rho_vel_uu[i][j]
                                 - (co_rho_vel_u[i][j]*co_rho_vel_u[i][j]/co_rho[i][j]);
                    LA12[i][j] = co_rho_vel_uv[i][j]
                                 - (co_rho_vel_u[i][j]*co_rho_vel_v[i][j]/co_rho[i][j]);
                    LA22[i][j] = co_rho_vel_vv[i][j]
                                 - (co_rho_vel_v[i][j]*co_rho_vel_v[i][j]/co_rho[i][j]);
		    LA1[i][j] = (LA11[i][j]/2.0) - (LA22[i][j]/2.0);
		    LA2[i][j] = LA12[i][j];

		    LI[i][j] = (LA11[i][j]+LA22[i][j])/2.0;
                    MI[i][j] = -(delta*delta*co_temp[i][j])
		               + (co_delta*co_delta*co_rho[i][j]*hat_SS*hat_SS);

                    MH1[i][j] = delta*delta*co_temp_t1[i][j]
                                - (co_delta*co_delta*co_rho[i][j]*co_cp[i][j]*hat_SS*hat_tx);
                    MH2[i][j] = delta*delta*co_temp_t2[i][j]
                                - (co_delta*co_delta*co_rho[i][j]*co_cp[i][j]*hat_SS*hat_ty);
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    LH1[i][j] = co_rho_cp_tu[i][j]
                                - ((co_rho_cp[i][j]*co_rho_t[i][j]*co_rho_vel_u[i][j])/(co_rho[i][j]*co_rho[i][j]));
                    LH2[i][j] = co_rho_cp_tv[i][j]
                                - ((co_rho_cp[i][j]*co_rho_t[i][j]*co_rho_vel_v[i][j])/(co_rho[i][j]*co_rho[i][j]));
                    MC1[i][j] = delta*delta*co_temp_c1[i][j]
                                - (co_delta*co_delta*co_rho[i][j]*hat_SS*hat_cx);
                    MC2[i][j] = delta*delta*co_temp_c2[i][j]
                                - (co_delta*co_delta*co_rho[i][j]*hat_SS*hat_cy);
                    MC10[i][j] = delta*delta*co_temp_c10[i][j]
                                - (co_delta*co_delta*co_rho[i][j]*hat_SS*hat_cx0);
                    MC20[i][j] = delta*delta*co_temp_c20[i][j]
                                - (co_delta*co_delta*co_rho[i][j]*hat_SS*hat_cy0);
                    LC1[i][j] = co_rho_cu[i][j]
                                - ((co_rho_vel_u[i][j]*co_rho_c[i][j])/co_rho[i][j]);
                    LC2[i][j] = co_rho_cv[i][j]
                                - ((co_rho_vel_v[i][j]*co_rho_c[i][j])/co_rho[i][j]);
                    LC10[i][j] = co_rho_cu0[i][j]
                                - ((co_rho_vel_u[i][j]*co_rho_c0[i][j])/co_rho[i][j]);
                    LC20[i][j] = co_rho_cv0[i][j]
                                - ((co_rho_vel_v[i][j]*co_rho_c0[i][j])/co_rho[i][j]);
                    }
                    //spherical average
                    r[i][j] = sqrt( (co_coords_x[i][j]*co_coords_x[i][j])
                                     + (co_coords_y[i][j]*co_coords_y[i][j]) );

                    //planar average
                    r[i][j] = co_coords_y[i][j];
	            r_color[i][j] = (int)(r[i][j]/delta_r);
                    r_color_LA[i][j] = (int)(r[i][j]/delta_r);
                    r_color_LI[i][j] = (int)(r[i][j]/delta_r);
                    r_color_LH[i][j] = (int)(r[i][j]/delta_r);
                    r_color_LC[i][j] = (int)(r[i][j]/delta_r);
                    r_color_LC0[i][j] = (int)(r[i][j]/delta_r);

                    if(((LA1[i][j]*MA11[i][j]) + (LA2[i][j]*MA12[i][j])) < 0.0)
                    {
                        r_color_LA[i][j] = -2;
                    }
                    if(LI[i][j] < 0.0 || MI[i][j] < 0.0)
                    {
                        r_color_LI[i][j] = -2;
                    }
                    if((LH1[i][j]*MH1[i][j]) + (LH2[i][j]*MH2[i][j]) < 0.0)
                    {
                        r_color_LH[i][j] = -2;
                    }
                    if((LC1[i][j]*MC1[i][j]) + (LC2[i][j]*MC2[i][j]) < 0.0)
                    {
                        r_color_LC[i][j] = -2;
                    }
                    if((LC10[i][j]*MC10[i][j]) + (LC20[i][j]*MC20[i][j]) < 0.0)
                    {
                        r_color_LC0[i][j] = -2;
                    }
               }
               else
               {
                    r_color[i][j] = -3;
                    r_color_LA[i][j] = -3;
                    r_color_LI[i][j] = -3;
                    r_color_LH[i][j] = -3;
                    r_color_LC[i][j] = -3;
                    r_color_LC0[i][j] = -3;
               }
            }
        }

        for (k = 0; k < num_r; k++)
        {
              deno[k] = 0.0;
              nume[k] = 0.0;
              ci_deno[k] = 0.0;
              ci_nume[k] = 0.0;
              prt_deno[k] = 0.0;
              prt_nume[k] = 0.0;
              sct_deno[k] = 0.0;
              sct_nume[k] = 0.0;
              sct_deno0[k] = 0.0;
              sct_nume0[k] = 0.0;
              cs_new_new[k] = 0.0;
              ci_new_new[k] = 0.0;
              prt_new_new[k] = 0.0;
              sct_new_new[k] = 0.0;
              sct0_new_new[k] = 0.0;
	      cs_new[k] = 0.0;
              ci_new[k] = 0.0;
              prt_new[k] = 0.0;
              sct_new[k] = 0.0;
              sct0_new[k] = 0.0;
        }

        for (k = 0; k < num_r; k++)
        {
           for (j = imin[1]+1; j < (imax[1]/2)+1; j++)
           {
               for (i = imin[0]+1; i < (imax[0]/2)+1; i++)
               {
                    if(co_NS_terms[i][j] == YES) 
                    {
                       if(k == r_color_LA[i][j])
                       {
                          deno[k] += ((MA11[i][j]*MA11[i][j]) + (MA12[i][j]*MA12[i][j]));
                          nume[k] += ((LA1[i][j]*MA11[i][j]) + (LA2[i][j]*MA12[i][j]));
                       }
                       if(k == r_color_LI[i][j])
                       {
                          ci_deno[k] += MI[i][j];
                          ci_nume[k] += LI[i][j];
                       }
                       if(k == r_color_LH[i][j])
                       {
                          prt_nume[k] += (LH1[i][j]*MH1[i][j]) + (LH2[i][j]*MH2[i][j]);
                          prt_deno[k] += (MH1[i][j]*MH1[i][j]) + (MH2[i][j]*MH2[i][j]);
                       }
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       {
                       if(k == r_color_LC[i][j])
                       {
                          sct_nume[k] += (LC1[i][j]*MC1[i][j]) + (LC2[i][j]*MC2[i][j]);
                          sct_deno[k] += (MC1[i][j]*MC1[i][j]) + (MC2[i][j]*MC2[i][j]);
                       }
                       if(k == r_color_LC0[i][j])
                       {
                          sct_nume0[k] += (LC10[i][j]*MC10[i][j]) + (LC20[i][j]*MC20[i][j]);
                          sct_deno0[k] += (MC10[i][j]*MC10[i][j]) + (MC20[i][j]*MC20[i][j]);
                       }
                       }
                    }
                }
            }
        }
       
        pp_gsync();

        if (nn > 1)
        {
           for (k = 0; k < num_r; k++)
           {
              pp_global_sum(&deno[k],1L);
              pp_global_sum(&nume[k],1L);
              pp_global_sum(&ci_deno[k],1L);
              pp_global_sum(&ci_nume[k],1L);
              pp_global_sum(&prt_deno[k],1L);
              pp_global_sum(&prt_nume[k],1L);
              pp_global_sum(&sct_deno[k],1L);
              pp_global_sum(&sct_nume[k],1L);
              pp_global_sum(&sct_deno0[k],1L);
              pp_global_sum(&sct_nume0[k],1L);
           }
        }

        for (k = 0; k < num_r; k++)
        {
            if(deno[k] < 10e-16)
            cs_new[k] = 0.0;
            else
            cs_new[k] = nume[k]/deno[k];
            if(ci_deno[k] < 10e-16)
            ci_new[k] = 0.0;
            else
            ci_new[k] = ci_nume[k]/ci_deno[k];
            if(prt_deno[k] < 10e-16)
            prt_new[k] = 0.0;
            else
            prt_new[k] = prt_nume[k]/prt_deno[k];
            if(sct_deno[k] < 10e-16)
            sct_new[k] = 0.0;
            else
            sct_new[k] = sct_nume[k]/sct_deno[k];
            if(sct_deno0[k] < 10e-16)
            sct0_new[k] = 0.0;
            else
            sct0_new[k] = sct_nume0[k]/sct_deno0[k];

            if(cs_new[k]<10e-16)
            {
                cs_new[k] = 0.0;
                ci_new[k] = 0.0;
                prt_new[k] = 0.0;
                sct_new[k] = 0.0;
                sct0_new[k] = 0.0;
            }
        }

        for (k = 0; k < num_r; k++)
        {
            for(i = -12; i < 13; i++)
            {
                cs_new_new[k] += cs_new[k+i];
                ci_new_new[k] += ci_new[k+i];
                prt_new_new[k] += prt_new[k+i];
                sct_new_new[k] += sct_new[k+i];
                sct0_new_new[k] += sct0_new[k+i];
            }

            cs_new_new[k] = cs_new_new[k]/25.0;
            ci_new_new[k] = ci_new_new[k]/25.0;
            prt_new_new[k] = prt_new_new[k]/25.0;
            sct_new_new[k] = sct_new_new[k]/25.0;
            sct0_new_new[k] = sct0_new_new[k]/25.0;
        }

        for (j = imin[1]; j < (imax[1]/2); j++)
        {
             for (i = imin[0]; i < (imax[0]/2); i++)
             {
                  if(j == imin[1] || j == (imax[1]/2)-1)
                  { 
                       co_C[i][j] = 0.0;
                       co_CI[i][j] = 0.0;
                       co_Prt[i][j] = 0.0;
                       co_Sct[i][j] = 0.0;
                       co_Sct0[i][j] = 0.0;
                  }
		  else
		  {
                       co_C[i][j] = cs_new_new[r_color[i+1][j+1]];
                       co_CI[i][j] = ci_new_new[r_color[i+1][j+1]];
                       co_Prt[i][j] = prt_new_new[r_color[i+1][j+1]];
                       co_Sct[i][j] = sct_new_new[r_color[i+1][j+1]];
                       co_Sct0[i][j] = sct0_new_new[r_color[i+1][j+1]];
		  }
             }
        }

        for (j = imin[1]; j < imax[1]; j++)
        {
            icoords[1] = j;
            for (i = imin[0]; i < imax[0]; i++)
            {
                icoords[0] = i;
                state = Rect_state(icoords,wv);
                if(g_compute_NS_terms(Rect_state(icoords,wv)))
                {
                    fill_2d_9pt_Pstencil(i,j,nsten);

                    ux = (vel(0,nsten->sts2d[2][1]) -
                          vel(0,nsten->sts2d[0][1])) / (2.0*dh[0]);
                    uy = (vel(0,nsten->sts2d[1][2]) -
                          vel(0,nsten->sts2d[1][0])) / (2.0*dh[1]);
                    vx = (vel(1,nsten->sts2d[2][1]) -
                          vel(1,nsten->sts2d[0][1])) / (2.0*dh[0]);
                    vy = (vel(1,nsten->sts2d[1][2]) -
                          vel(1,nsten->sts2d[1][0])) / (2.0*dh[1]);

                    s11 = ux;
                    s12 = (uy + vx)/2;
                    s21 = (uy + vx)/2;
                    s22 = vy;

                    dtemx = (temperature(nsten->sts2d[2][1])
                            - temperature(nsten->sts2d[0][1]))/(2.0*dh[0]);
                    dtemy = (temperature(nsten->sts2d[1][2])
                            - temperature(nsten->sts2d[1][0]))/(2.0*dh[1]);
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    pdenx = (pdens(nsten->sts2d[2][1])[1]/Dens(nsten->sts2d[2][1])
                             - pdens(nsten->sts2d[0][1])[1]/Dens(nsten->sts2d[0][1]))/(2.0*dh[0]);
                    pdeny = (pdens(nsten->sts2d[1][2])[1]/Dens(nsten->sts2d[1][2])
                             - pdens(nsten->sts2d[1][0])[1]/Dens(nsten->sts2d[1][0]))/(2.0*dh[1]);
                    pdenx0 = (pdens(nsten->sts2d[2][1])[0]/Dens(nsten->sts2d[2][1])
                             - pdens(nsten->sts2d[0][1])[0]/Dens(nsten->sts2d[0][1]))/(2.0*dh[0]);
                    pdeny0 = (pdens(nsten->sts2d[1][2])[0]/Dens(nsten->sts2d[1][2])
                             -pdens(nsten->sts2d[1][0])[0]/Dens(nsten->sts2d[1][0]))/(2.0*dh[1]);
                    }
                    S = sqrt(2*( (s11*s11) + (s12*s12)
                             + (s21*s21) + (s22*s22)));

                    Tau(state)[0][0] = ((2/2)*CI[i][j]*Dens(nsten->sts2d[1][1])*delta*delta*S*S)
                                        - 2*Dens(nsten->sts2d[1][1])*C[i][j]*delta*delta*S*((s11/2)-(s22/2));
                    Tau(state)[0][1] = - 2*Dens(nsten->sts2d[1][1])*C[i][j]*delta*delta*S*(s12);
                    Tau(state)[1][0] = - 2*Dens(nsten->sts2d[1][1])*C[i][j]*delta*delta*S*(s21);
                    Tau(state)[1][1] = ((2/2)*CI[i][j]*Dens(nsten->sts2d[1][1])*delta*delta*S*S)
                                        - 2*Dens(nsten->sts2d[1][1])*C[i][j]*delta*delta*S*((s22/2)-(s11/2));

                    Qh(state)[0] = - (Dens(nsten->sts2d[1][1])*C_P(nsten->sts2d[1][1])*Prt[i][j]*delta*delta*S)*dtemx;
                    Qh(state)[1] = - (Dens(nsten->sts2d[1][1])*C_P(nsten->sts2d[1][1])*Prt[i][j]*delta*delta*S)*dtemy;
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    Qc(state)[0] = - (Dens(nsten->sts2d[1][1])*Sct[i][j]*delta*delta*S)*pdenx;
                    Qc(state)[0] = - (Dens(nsten->sts2d[1][1])*Sct[i][j]*delta*delta*S)*pdeny;
                    Qc(state)[0] = - (Dens(nsten->sts2d[1][1])*Sct0[i][j]*delta*delta*S)*pdenx0;
                    Qc(state)[0] = - (Dens(nsten->sts2d[1][1])*Sct0[i][j]*delta*delta*S)*pdeny0;
                    }
                    //C(state) = C[i][j];
                    //CI(state) = CI[i][j];   
                    //SCT(state) = Sct[i][j];
                    //SCT0(state) = Sct0[i][j];    
                    //PRT(state) = Prt[i][j];
                    //STR(state) = S;
                }
            }
        }

        free(re_NS_terms);
        free(co_NS_terms);    
        free(co_coords_x);
        free(co_coords_y);
        free(co_rho);
        free(co_rho_cp);
        free(co_cp);
        free(co_rho_t);
        free(co_rho_c);
        free(co_rho_c0);
        free(co_t);
        free(co_con);
        free(co_con0);
        free(co_rho_cu);
        free(co_rho_cv);
        free(co_rho_cu0);
        free(co_rho_cv0);
        free(co_rho_cp_tu);
        free(co_rho_cp_tv);
        free(co_vel_u);
        free(co_vel_v);
        free(co_rho_vel_u);
        free(co_rho_vel_v);
        free(co_rho_vel_uv);
        free(co_rho_vel_uu);
        free(co_rho_vel_vv);
        free(co_temp11);
        free(co_temp12);
        free(co_temp21);
        free(co_temp22);
        free(co_temp);
        free(co_temp_t1);
        free(co_temp_t2);
        free(co_temp_c1);
        free(co_temp_c2);
        free(co_temp_c10);
        free(co_temp_c20);
        free(re_temp11);
        free(re_temp12);
        free(re_temp21);
        free(re_temp22);
        free(re_temp_t1);
        free(re_temp_t2);
        free(re_temp_c1);
        free(re_temp_c2);
        free(re_temp_c10);
        free(re_temp_c20);
        free(MA11);
        free(MA12);
        free(MA21);
        free(MA22);
        free(MI);
        free(LA11);
        free(LA12);
        free(LA22);
        free(LA1);
        free(LA2);
        free(LI);
        free(MH1);
        free(MH2);
        free(LH1);
        free(LH2);
        free(MC1);
        free(MC2);
        free(LC1);
        free(LC2);
        free(MC10);
        free(MC20);
        free(LC10);
        free(LC20);
        free(r);
        free(r_color);
        free(cs_new);
        free(ci_new);
        free(prt_new);
        free(sct_new);
        free(sct0_new);
        free(cs_new_new);
        free(ci_new_new);
        free(prt_new_new);
        free(sct_new_new);
        free(sct0_new_new);
        free(re_coords_x);
        free(re_coords_y);
        free(re_rho);
        free(re_rho_cp);
        free(re_cp);
        free(re_rho_t);
        free(re_rho_c);
        free(re_rho_c0);
        free(re_t);
        free(re_c);
        free(re_c0);
        free(re_rho_cu);
        free(re_rho_cv);
        free(re_rho_cu0);
        free(re_rho_cv0);
        free(re_rho_cp_tu);
        free(re_rho_cp_tv);
        free(re_vel_u);
        free(re_vel_v);
        free(re_rho_vel_u);
        free(re_rho_vel_v);
        free(re_rho_vel_uv);
        free(re_rho_vel_uu);
        free(re_rho_vel_vv);
        free(re_temp);
        free(C);
        free(co_C);
        free(CI);
        free(co_CI);
        free(Prt);
        free(co_Prt);
        free(Sct);
        free(co_Sct);
        free(Sct0);
        free(co_Sct0);
        free(r_color_LA);
        free(r_color_LH);
        free(r_color_LC);
        free(r_color_LC0);
        free(r_color_LI);
        free(deno);
        free(nume);
        free(ci_deno);
        free(ci_nume);
        free(prt_deno);
        free(prt_nume);
        free(sct_deno);
        free(sct_nume);
        free(sct_deno0);
        free(sct_nume0);
}

char *  get_directory();
LOCAL   void   SGS3d(
        float   dt,
        Front   *front,
        Wave    *wv)
{
        RECT_GRID       *rgr = wv->rect_grid;
        float           *GL = rgr->GL, *GU = rgr->GU;
        float           *L = rgr->L, *U = rgr->U;
        int             dim = rgr->dim;
        int             *gmax = rgr->gmax;
        int             i,j,k;
        int             icoords[MAXD];
        static          Pstencil *nsten;
        static          Pstencil *nsten2;
        int             npts,mid;
        int             imin[MAXD],imax[MAXD];
        size_t          sizest = front->sizest;
        bool            viscosity, mass_diffusion;
        float           ux,uy,uz,vx,vy,vz,wx,wy,wz;
        float           tx,ty,tz,cx[MAX_NCOMPS],cy[MAX_NCOMPS],cz[MAX_NCOMPS],S11,S12,S13,S21,S22,S23,S31,S32,S33,SS;
        float           hat_ux,hat_uy,hat_uz,hat_vx,hat_vy,hat_vz,hat_wx,hat_wy,hat_wz,hat_tx,hat_ty,hat_tz,hat_cx[MAX_NCOMPS],hat_cy[MAX_NCOMPS],hat_cz[MAX_NCOMPS];
        float           hat_S11,hat_S12,hat_S13,hat_S21,hat_S22,hat_S23,hat_S31,hat_S32,hat_S33,hat_SS;
        float           hat_S11_a,hat_S12_a,hat_S13_a,hat_S21_a,hat_S22_a,hat_S23_a,hat_S31_a,hat_S32_a,hat_S33_a;
        float           s11, s12, s13, s21, s22, s23, s31, s32, s33, dtemx, dtemy, dtemz, pdenx[MAX_NCOMPS], pdeny[MAX_NCOMPS], pdenz[MAX_NCOMPS], S;
        float           S11_a, S12_a, S13_a, S21_a, S22_a, S23_a, S31_a, S32_a, S33_a;
        const int       nn = pp_numnodes();
        int             myid = pp_mynode();
        int             *ppgmax = front->pp_grid->gmax;
        int             ppx = myid % ppgmax[0];
        int             ppy = ((myid-ppx)/ppgmax[0]) % ppgmax[1];
        int             ppz = (myid-ppx-(ppgmax[0]*ppy))/(ppgmax[0]*ppgmax[1]);
        int             mink, maxk;

        int             av[3] = {3,3,3};  /* Should go to inputfile, transmitted to both SGS3d() and mesh_avr() */

        int             avmax[3], avmin[3] = {0,0,0};

        for (i = 0; i < dim; i++)
        {
            imin[i] = 0;        imax[i] = gmax[i];
        }
/*
        if(ppz == ppgmax[2]-1)
        {
            maxk = imax[2];
            if(ppz == 0)
                mink = imin[2];
            else
                mink = imin[2]-1;
        }
        else if(ppz == 0)
        {
            mink = imin[2];
            if(ppz == ppgmax[2]-1)
                maxk = imax[2];
            else
                maxk = imax[2]+1;
        }
        else
        {
            mink = imin[2]-1;
            maxk = imax[2]+1;
        }
*/
        if (nsten == NULL)
        {
            stat_scalar(&nsten,sizeof(Pstencil));
            stat_scalar(&nsten2,sizeof(Pstencil));
            switch (dim)
            {
            case 1:
                bi_array(&nsten->icoords1d,3,1,INT);
                uni_array(&nsten->sts1d,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                    alloc_state(front->interf,&nsten->sts1d[i],front->sizest);
                break;
            case 2:
                tri_array(&nsten->icoords2d,3,3,2,INT);
                bi_array(&nsten->sts2d,3,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                {
                    for (j = 0; j < 3; j++)
                    {
                        alloc_state(front->interf,&nsten->sts2d[i][j],
                                front->sizest);
                    }
                }
                break;
            case 3:
                quad_array(&nsten->icoords3d,3,3,3,3,INT);
                quad_array(&nsten2->icoords3d,3,3,3,3,INT);
                tri_array(&nsten->sts3d,3,3,3,sizeof(Locstate));
                tri_array(&nsten2->sts3d,3,3,3,sizeof(Locstate));
                for (i = 0; i < 3; i++){
                    for (j = 0; j < 3; j++){
                        for (k = 0; k < 3; k++){
                            alloc_state(front->interf,&nsten->sts3d[i][j][k],
                                     front->sizest);
                            alloc_state(front->interf,&nsten2->sts3d[i][j][k],
                                     front->sizest);
			}
			}
			}
                break;
            }
        }
        nsten->fr = front;
        nsten->wave = wv;
        nsten->npts = npts = wv->npts_sten;
        mid = (int) (npts/2);

        Locstate state;
        float    ***C, ***co_C, ***CI, ***co_CI, ***Prt, ***co_Prt;
        float    ****Sct, ****co_Sct;
        float    ***co_coords_x, ***co_coords_y, ***co_coords_z;
        float    ***co_rho, ***co_vel_u, ***co_vel_v, ***co_vel_w, ***co_rho_vel_u, ***co_rho_vel_v, ***co_rho_vel_w, ***co_rho_vel_uv,***co_rho_vel_uw,***co_rho_vel_vw,***co_rho_vel_uu,***co_rho_vel_vv,***co_rho_vel_ww;
        float    ***co_rho_cp, ***co_cp, ***co_rho_t, ****co_rho_c,***co_t, ****co_con;
        float    ****co_rho_cu, ****co_rho_cv, ****co_rho_cw, ***co_rho_cp_tu, ***co_rho_cp_tv, ***co_rho_cp_tw;
        float    ***co_temp11, ***co_temp12, ***co_temp13, ***co_temp21, ***co_temp22, ***co_temp23, ***co_temp31, ***co_temp32, ***co_temp33, ***co_temp, ***co_temp_t1, ***co_temp_t2, ***co_temp_t3, ****co_temp_c1, ****co_temp_c2, ****co_temp_c3;
        float    ***MA11, ***MA12, ***MA13, ***MA21, ***MA22, ***MA23, ***MA31, ***MA32, ***MA33, ***MI, ***LA11, ***LA12, ***LA13, ***LA21, ***LA22, ***LA23, ***LA31, ***LA32, ***LA33, ***LA1, ***LA2, ***LA3, ***LA4, ***LA5, ***LI, ***MH1, ***MH2, ***MH3, ***LH1, ***LH2, ***LH3, ****MC1, ****MC2, ****MC3, ****LC1, ****LC2, ****LC3;
        float    ***r, ***cs_new, ***ci_new, ***prt_new, ****sct_new, *CS_LA1, *CS_LA2, *CS_MA1, *CS_MA2, *CI_LI, *CI_MI,*PRT_LH1, *PRT_LH2, *PRT_MH1, *PRT_MH2, **SCT_LC1, **SCT_LC2, **SCT_MC1, **SCT_MC2;
	float    ***MMA, ***LMA, ***LMH, ***MMH, ****LMC, ****MMC;

 /*       int      ***r_color, ***r_color_LA,***r_color_LH,***r_color_LC,***r_color_LC0,***r_color_LI;  */
        float    ***re_coords_x, ***re_coords_y, ***re_coords_z;
        float    ***re_rho, ***re_vel_u, ***re_vel_v, ***re_vel_w, ***re_rho_vel_u, ***re_rho_vel_v, ***re_rho_vel_w, ***re_rho_vel_uv,***re_rho_vel_uw,***re_rho_vel_vw,***re_rho_vel_uu,***re_rho_vel_vv,***re_rho_vel_ww;
        float    ***re_rho_cp, ***re_cp, ***re_rho_t, ****re_rho_c,***re_t, ****re_con;
        float    ****re_rho_cu, ****re_rho_cv, ****re_rho_cw, ***re_rho_cp_tu, ***re_rho_cp_tv, ***re_rho_cp_tw;
        float    ***re_temp11, ***re_temp12, ***re_temp13, ***re_temp21, ***re_temp22, ***re_temp23, ***re_temp31, ***re_temp32, ***re_temp33, ***re_temp, ***re_temp_t1, ***re_temp_t2, ***re_temp_t3, ****re_temp_c1, ****re_temp_c2, ****re_temp_c3;

         bool     ***re_NS_terms, ***co_NS_terms;

         float    ***deno, ***nume;
         float    ***ci_deno, ***ci_nume;
         float    ***prt_deno, ***prt_nume;
         float    ****sct_deno, ****sct_nume;

        int ii, jj, kk, iii, jjj, kkk,ic;
        float  dh[3];
        dh[0] = nsten->wave->rect_grid->h[0];
        dh[1] = nsten->wave->rect_grid->h[1];
        dh[2] = nsten->wave->rect_grid->h[2];

        float delta = pow((dh[0]*dh[1]*dh[2]),(1.0/3.0));
        float co_delta = pow(((2*dh[0])*(2*dh[1])*(2*dh[2])),(1.0/3.0));
        float delta_r = 2*dh[2];
        int num_r = gmax[2]*ppgmax[2]/2;

        float      sum_co_rho=0.0,sum_co_vel_u=0.0,sum_co_vel_v=0.0,sum_co_vel_w=0.0;
        float      sum_co_rho_vel_u=0.0,sum_co_rho_vel_v=0.0,sum_co_rho_vel_w=0.0;
        float      sum_co_rho_vel_uu=0.0,sum_co_rho_vel_vv=0.0,sum_co_rho_vel_ww=0.0;
        float      sum_co_rho_vel_uv=0.0,sum_co_rho_vel_vw=0.0,sum_co_rho_vel_uw=0.0;
        float      sum_co_rho_cp=0.0,sum_co_cp=0.0,sum_co_rho_t=0.0,sum_co_t=0.0;
        float      sum_co_rho_cp_tu=0.0,sum_co_rho_cp_tv=0.0,sum_co_rho_cp_tw=0.0;
        float      sum_co_temp11=0.0,sum_co_temp12=0.0,sum_co_temp13=0.0;
        float      sum_co_temp21=0.0,sum_co_temp22=0.0,sum_co_temp23=0.0;
        float      sum_co_temp31=0.0,sum_co_temp32=0.0,sum_co_temp33=0.0;
        float      sum_co_temp=0.0,sum_co_temp_t1=0.0,sum_co_temp_t2=0.0,sum_co_temp_t3=0.0;
        float      sum_co_rho_c[MAX_NCOMPS],sum_co_con[MAX_NCOMPS];
        float      sum_co_rho_cu[MAX_NCOMPS],sum_co_rho_cv[MAX_NCOMPS],sum_co_rho_cw[MAX_NCOMPS];
        float      sum_co_temp_c1[MAX_NCOMPS],sum_co_temp_c2[MAX_NCOMPS],sum_co_temp_c3[MAX_NCOMPS];
        for (ic = 0; ic < MAX_NCOMPS; ic++)
        {
            sum_co_rho_c[ic] = 0.0;
            sum_co_con[ic] = 0.0;
            sum_co_rho_cu[ic] = 0.0;
            sum_co_rho_cv[ic] = 0.0;
            sum_co_rho_cw[ic] = 0.0;
            sum_co_temp_c1[ic] = 0.0;
            sum_co_temp_c2[ic] = 0.0;
            sum_co_temp_c3[ic] = 0.0;
        }

        int        iimin[3] = {imin[0]+1, imin[1]+1, imin[2]+1}, iimax[3] = {(imax[0]-imin[0])/2+1, (imax[1]-imin[1])/2+1, (imax[2]-imin[2])/2+1};

        for (i = 0; i < 3; i++)
        {
            avmax[i] = (imax[i] - imin[i])/2 - av[i] + 1;
        }



        switch (dim)
        {
        case 3:
	    tri_array(&co_NS_terms,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_coords_x,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_coords_y,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_coords_z,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_rho,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_rho_cp,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_cp,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_rho_t,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            quad_array(&co_rho_c,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
	    tri_array(&co_t,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    quad_array(&co_con,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
	    quad_array(&co_rho_cu,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
	    quad_array(&co_rho_cv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&co_rho_cw,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            tri_array(&co_rho_cp_tu,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_rho_cp_tv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_rho_cp_tw,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_vel_u,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_vel_v,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_vel_w,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_rho_vel_u,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_rho_vel_v,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_rho_vel_w,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_rho_vel_uu,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_rho_vel_vv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_rho_vel_ww,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_rho_vel_uv,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_rho_vel_vw,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_rho_vel_uw,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_temp11,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_temp12,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_temp13,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_temp21,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_temp22,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_temp23,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_temp31,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_temp32,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_temp33,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_temp,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_temp_t1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&co_temp_t2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&co_temp_t3,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    quad_array(&co_temp_c1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
	    quad_array(&co_temp_c2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&co_temp_c3,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));

            tri_array(&MA11,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&MA12,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&MA13,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&MA22,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&MA23,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA11,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA12,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA13,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA22,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA23,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA33,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA3,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA4,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LA5,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&LI,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&MI,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&MH1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&MH2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
	    tri_array(&MH3,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&LH1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&LH2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&LH3,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            quad_array(&MC1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&MC2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&MC3,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&LC1,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&LC2,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&LC3,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));

            tri_array(&MMA,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&LMA,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&LMH,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            tri_array(&MMH,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
            quad_array(&LMC,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));
            quad_array(&MMC,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,MAX_NCOMPS,sizeof(float));

            tri_array(&r,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(float));
 /*           tri_array(&r_color,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(int));
            tri_array(&r_color_LA,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(int));
            tri_array(&r_color_LI,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(int));
            tri_array(&r_color_LH,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(int));
            tri_array(&r_color_LC,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(int));
            tri_array(&r_color_LC0,((imax[0]-imin[0])/2)+2,((imax[1]-imin[1])/2)+2,((imax[2]-imin[2])/2)+2,sizeof(int));  

	    uni_array(&cs_new,num_r,sizeof(float));
            uni_array(&ci_new,num_r,sizeof(float));
            uni_array(&prt_new,num_r,sizeof(float));
            uni_array(&sct_new,num_r,sizeof(float));
            uni_array(&sct0_new,num_r,sizeof(float));
            uni_array(&deno,num_r,sizeof(float));
            uni_array(&nume,num_r,sizeof(float));
            uni_array(&ci_deno,num_r,sizeof(float));
            uni_array(&ci_nume,num_r,sizeof(float));
            uni_array(&prt_deno,num_r,sizeof(float));
            uni_array(&prt_nume,num_r,sizeof(float));
            uni_array(&sct_deno,num_r,sizeof(float));
            uni_array(&sct_nume,num_r,sizeof(float));
            uni_array(&sct_deno0,num_r,sizeof(float));
            uni_array(&sct_nume0,num_r,sizeof(float));  */

	    tri_array(&deno,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
	    tri_array(&nume,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
	    tri_array(&ci_deno,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
	    tri_array(&ci_nume,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
	    tri_array(&prt_nume,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
	    tri_array(&prt_deno,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
	    quad_array(&sct_nume,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],MAX_NCOMPS,sizeof(float));
	    quad_array(&sct_deno,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],MAX_NCOMPS,sizeof(float));

            tri_array(&cs_new,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
            tri_array(&ci_new,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
            tri_array(&prt_new,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],sizeof(float));
            quad_array(&sct_new,avmax[0]-avmin[0],avmax[1]-avmin[1],avmax[2]-avmin[2],MAX_NCOMPS,sizeof(float));

            tri_array(&re_NS_terms,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_coords_x,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_coords_y,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_coords_z,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_rho,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_rho_cp,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_cp,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_rho_t,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            quad_array(&re_rho_c,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));
	    tri_array(&re_t,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    quad_array(&re_con,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));
	    quad_array(&re_rho_cu,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));
	    quad_array(&re_rho_cv,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));
            quad_array(&re_rho_cw,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));
            tri_array(&re_rho_cp_tu,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_rho_cp_tv,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_rho_cp_tw,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_vel_u,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_vel_v,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_vel_w,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_rho_vel_u,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_rho_vel_v,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_rho_vel_w,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_rho_vel_uu,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_rho_vel_vv,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_rho_vel_ww,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_rho_vel_uv,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_rho_vel_vw,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_rho_vel_uw,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_temp11,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_temp12,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_temp13,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_temp21,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_temp22,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_temp23,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_temp31,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_temp32,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_temp33,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_temp,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_temp_t1,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    tri_array(&re_temp_t2,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
            tri_array(&re_temp_t3,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,sizeof(float));
	    quad_array(&re_temp_c1,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));
	    quad_array(&re_temp_c2,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));
            quad_array(&re_temp_c3,(imax[0]-imin[0])+4,(imax[1]-imin[1])+4,(imax[2]-imin[2])+4,MAX_NCOMPS,sizeof(float));

            tri_array(&C,(imax[0]-imin[0]),(imax[1]-imin[1]),(imax[2]-imin[2]),sizeof(float));
	    tri_array(&co_C,((imax[0]-imin[0])/2),((imax[1]-imin[1])/2),((imax[2]-imin[2])/2),sizeof(float));
            tri_array(&CI,(imax[0]-imin[0]),(imax[1]-imin[1]),(imax[2]-imin[2]),sizeof(float));
	    tri_array(&co_CI,((imax[0]-imin[0])/2),((imax[1]-imin[1])/2),((imax[2]-imin[2])/2),sizeof(float));
            tri_array(&Prt,(imax[0]-imin[0]),(imax[1]-imin[1]),(imax[2]-imin[2]),sizeof(float));
	    tri_array(&co_Prt,((imax[0]-imin[0])/2),((imax[1]-imin[1])/2),((imax[2]-imin[2])/2),sizeof(float));
            quad_array(&Sct,(imax[0]-imin[0]),(imax[1]-imin[1]),(imax[2]-imin[2]),MAX_NCOMPS,sizeof(float));
	    quad_array(&co_Sct,((imax[0]-imin[0])/2),((imax[1]-imin[1])/2),((imax[2]-imin[2])/2),MAX_NCOMPS,sizeof(float));
        }
	
	float cp, T, v0, v1, v2;

//        for (k = mink; k < maxk; k++)
        for (k = imin[2]-2; k < imax[2]+2; k++)
        {
        icoords[2] = k;
        for (j = imin[1]-2; j < imax[1]+2; j++)
        {
            icoords[1] = j;
            for (i = imin[0]-2; i < imax[0]+2; i++)
            {
                icoords[0] = i;
                state =  Rect_state(icoords,wv);

                if((ppx == 0 && i == imin[0]-2) || (ppy == 0 && j == imin[1]-2)  || (ppz == 0 && k == imin[2]-2)  
                  || (ppx == ppgmax[0]-1 && i == imax[0]+1) || (ppy == ppgmax[1]-1 && j == imax[1]+1)
		  || (ppz == ppgmax[2]-1 && k == imax[2]+1))
                {
                    re_NS_terms[i+2][j+2][k+2] = NO;
                    continue;
                }

                if(g_compute_NS_terms(Rect_state(icoords,wv)))
                {
                    fill_3d_27pt_Pstencil(i,j,k,nsten);

                    ux = (vel(0,nsten->sts3d[2][1][1]) -
                                vel(0,nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                    uy = (vel(0,nsten->sts3d[1][2][1]) -
                                vel(0,nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                    uz = (vel(0,nsten->sts3d[1][1][2]) -
                                vel(0,nsten->sts3d[1][1][0])) / (2.0*dh[2]);

                    vx = (vel(1,nsten->sts3d[2][1][1]) -
                                vel(1,nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                    vy = (vel(1,nsten->sts3d[1][2][1]) -
                                vel(1,nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                    vz = (vel(1,nsten->sts3d[1][1][2]) -
                                vel(1,nsten->sts3d[1][1][0])) / (2.0*dh[2]);

                    wx = (vel(2,nsten->sts3d[2][1][1]) -
                                vel(2,nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                    wy = (vel(2,nsten->sts3d[1][2][1]) -
                                vel(2,nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                    wz = (vel(2,nsten->sts3d[1][1][2]) -
                                vel(2,nsten->sts3d[1][1][0])) / (2.0*dh[2]);

                    tx = (temperature(nsten->sts3d[2][1][1]) -
                          temperature(nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                    ty = (temperature(nsten->sts3d[1][2][1]) -
                          temperature(nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                    tz = (temperature(nsten->sts3d[1][1][2]) -
                          temperature(nsten->sts3d[1][1][0])) / (2.0*dh[2]);

		    v0 = vel(0,state);
		    v1 = vel(1,state);
		    v2 = vel(2,state);
		    cp = C_P(state);
		    T = temperature(state);

                   if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
	            {
                    for (ic = 0; ic < MAX_NCOMPS; ic++)
                    {
                    cx[ic] = ((pdens(nsten->sts3d[2][1][1])[ic]/Dens(nsten->sts3d[2][1][1])) -
                          (pdens(nsten->sts3d[0][1][1])[ic]/Dens(nsten->sts3d[0][1][1]))) / (2.0*dh[0]);
                    cy[ic] = ((pdens(nsten->sts3d[1][2][1])[ic]/Dens(nsten->sts3d[1][2][1])) -
                          (pdens(nsten->sts3d[1][0][1])[ic]/Dens(nsten->sts3d[1][0][1]))) / (2.0*dh[1]);
                    cz[ic] = ((pdens(nsten->sts3d[1][1][2])[ic]/Dens(nsten->sts3d[1][1][2])) -
                          (pdens(nsten->sts3d[1][1][0])[ic]/Dens(nsten->sts3d[1][1][0]))) / (2.0*dh[2]);
                    re_con[i+2][j+2][k+2][ic] = pdens(state)[ic]/Dens(state);
                    
		    re_rho_cu[i+2][j+2][k+2][ic] = Dens(state)*(pdens(state)[ic]/Dens(state))*v0;
                    re_rho_cv[i+2][j+2][k+2][ic] = Dens(state)*(pdens(state)[ic]/Dens(state))*v1;
                    re_rho_cw[i+2][j+2][k+2][ic] = Dens(state)*(pdens(state)[ic]/Dens(state))*v2;
		    
		    /*re_rho_cu[i+2][j+2][k+2][ic] = Dens(state)*(pdens(state)[ic]/Dens(state))*vel(0,state);
                    re_rho_cv[i+2][j+2][k+2][ic] = Dens(state)*(pdens(state)[ic]/Dens(state))*vel(1,state);
                    re_rho_cw[i+2][j+2][k+2][ic] = Dens(state)*(pdens(state)[ic]/Dens(state))*vel(2,state);*/
                    
		    re_rho_c[i+2][j+2][k+2][ic] = Dens(state)*(pdens(state)[ic]/Dens(state));
                    }
                    }

                    re_coords_x[i+2][j+2][k+2] = cell_center(i,0,rgr) - GL[0];
                    re_coords_y[i+2][j+2][k+2] = cell_center(j,1,rgr) - GL[1];
                    re_coords_z[i+2][j+2][k+2] = cell_center(k,2,rgr) - GL[2];
	
                    re_rho[i+2][j+2][k+2] = Dens(state);
                    
		    re_cp[i+2][j+2][k+2] = cp;
                    re_t[i+2][j+2][k+2] = T;
                    re_vel_u[i+2][j+2][k+2] = v0;
                    re_vel_v[i+2][j+2][k+2] = v1;
                    re_vel_w[i+2][j+2][k+2] = v2;
                    re_rho_cp_tu[i+2][j+2][k+2] = Dens(state)*cp*T*v0;
                    re_rho_cp_tv[i+2][j+2][k+2] = Dens(state)*cp*T*v1;
                    re_rho_cp_tw[i+2][j+2][k+2] = Dens(state)*cp*T*v2;
                    re_rho_cp[i+2][j+2][k+2] = Dens(state)*cp;
                    re_rho_t[i+2][j+2][k+2] = Dens(state)*T;
                    re_rho_vel_u[i+2][j+2][k+2] = Dens(state)*v0;
                    re_rho_vel_v[i+2][j+2][k+2] = Dens(state)*v1;
                    re_rho_vel_w[i+2][j+2][k+2] = Dens(state)*v2;
                    re_rho_vel_uu[i+2][j+2][k+2] = Dens(state)*v0*v0;
                    re_rho_vel_vv[i+2][j+2][k+2] = Dens(state)*v1*v1;
                    re_rho_vel_ww[i+2][j+2][k+2] = Dens(state)*v2*v2;
                    re_rho_vel_uv[i+2][j+2][k+2] = Dens(state)*v0*v1;
                    re_rho_vel_vw[i+2][j+2][k+2] = Dens(state)*v1*v2;
                    re_rho_vel_uw[i+2][j+2][k+2] = Dens(state)*v0*v2;

		    /*re_cp[i+2][j+2][k+2] = C_P(state);
                    re_t[i+2][j+2][k+2] = temperature(state);
                    re_vel_u[i+2][j+2][k+2] = vel(0,state);
                    re_vel_v[i+2][j+2][k+2] = vel(1,state);
                    re_vel_w[i+2][j+2][k+2] = vel(2,state);
                    re_rho_cp_tu[i+2][j+2][k+2] = Dens(state)*C_P(state)*temperature(state)*vel(0,state);
                    re_rho_cp_tv[i+2][j+2][k+2] = Dens(state)*C_P(state)*temperature(state)*vel(1,state);
                    re_rho_cp_tw[i+2][j+2][k+2] = Dens(state)*C_P(state)*temperature(state)*vel(2,state);
                    re_rho_cp[i+2][j+2][k+2] = Dens(state)*C_P(state);
                    re_rho_t[i+2][j+2][k+2] = Dens(state)*temperature(state);
                    re_rho_vel_u[i+2][j+2][k+2] = Dens(state)*vel(0,state);
                    re_rho_vel_v[i+2][j+2][k+2] = Dens(state)*vel(1,state);
                    re_rho_vel_w[i+2][j+2][k+2] = Dens(state)*vel(2,state);
                    re_rho_vel_uu[i+2][j+2][k+2] = Dens(state)*vel(0,state)*vel(0,state);
                    re_rho_vel_vv[i+2][j+2][k+2] = Dens(state)*vel(1,state)*vel(1,state);
                    re_rho_vel_ww[i+2][j+2][k+2] = Dens(state)*vel(2,state)*vel(2,state);
                    re_rho_vel_uv[i+2][j+2][k+2] = Dens(state)*vel(0,state)*vel(1,state);
                    re_rho_vel_vw[i+2][j+2][k+2] = Dens(state)*vel(1,state)*vel(2,state);
                    re_rho_vel_uw[i+2][j+2][k+2] = Dens(state)*vel(0,state)*vel(2,state);*/

                    S11 = ux;
                    S12 = (uy + vx)/2.0;
                    S13 = (uz + wx)/2.0;
                    S21 = (vx + uy)/2.0;
                    S22 = vy;
                    S23 = (vz + wy)/2.0;
                    S31 = (wx + uz)/2.0;
                    S32 = (wy + vz)/2.0;
                    S33 = wz;

                    SS = sqrt(2.0*( (S11*S11) + (S12*S12)
                                    + (S13*S13) + (S21*S21)
                                    + (S22*S22) + (S23*S23)
                                    + (S31*S31) + (S32*S32)
                                    + (S33*S33)));

                    S11_a = S11 - ((S11+S22+S33)/3.0);
                    S12_a = S12;
                    S13_a = S13;
                    S21_a = S21;
                    S22_a = S22 - ((S11+S22+S33)/3.0);
                    S23_a = S23;
                    S31_a = S31;
                    S32_a = S32;
                    S33_a = S33 - ((S11+S22+S33)/3.0);

                    re_temp11[i+2][j+2][k+2] = Dens(state)*SS*(S11_a);
                    re_temp12[i+2][j+2][k+2] = Dens(state)*SS*S12_a;
                    re_temp13[i+2][j+2][k+2] = Dens(state)*SS*S13_a;
                    re_temp21[i+2][j+2][k+2] = Dens(state)*SS*S21_a;
                    re_temp22[i+2][j+2][k+2] = Dens(state)*SS*(S22_a);
                    re_temp23[i+2][j+2][k+2] = Dens(state)*SS*(S23_a);
                    re_temp31[i+2][j+2][k+2] = Dens(state)*SS*S31_a;
                    re_temp32[i+2][j+2][k+2] = Dens(state)*SS*(S32_a);
                    re_temp33[i+2][j+2][k+2] = Dens(state)*SS*(S33_a);
                    re_temp[i+2][j+2][k+2] = Dens(state)*SS*SS;
                    re_temp_t1[i+2][j+2][k+2] = Dens(state)*cp*SS*tx;
                    re_temp_t2[i+2][j+2][k+2] = Dens(state)*cp*SS*ty;
                    re_temp_t3[i+2][j+2][k+2] = Dens(state)*cp*SS*tz;

                    /*re_temp_t1[i+2][j+2][k+2] = Dens(state)*C_P(state)*SS*tx;
                    re_temp_t2[i+2][j+2][k+2] = Dens(state)*C_P(state)*SS*ty;
                    re_temp_t3[i+2][j+2][k+2] = Dens(state)*C_P(state)*SS*tz;*/

                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    for (ic = 0; ic < MAX_NCOMPS; ic++)
                    {
                    re_temp_c1[i+2][j+2][k+2][ic] = Dens(state)*SS*cx[ic];
                    re_temp_c2[i+2][j+2][k+2][ic] = Dens(state)*SS*cy[ic];
                    re_temp_c3[i+2][j+2][k+2][ic] = Dens(state)*SS*cz[ic];
                    }
                    }

                    re_NS_terms[i+2][j+2][k+2] = YES;
                }
                else
                {
                    re_NS_terms[i+2][j+2][k+2] = NO;
                }
            }
        }
        }

        for (k = imin[2]; k < (imax[2]/2)+2; k++)
        {
        kk = 2*k;
        for (j = imin[1]; j < (imax[1]/2)+2; j++)
        {
            jj = 2*j;
            for (i = imin[0]; i < (imax[0]/2)+2; i++)
            {
                ii = 2*i;
                if(re_NS_terms[ii][jj][kk] == YES ||
                   re_NS_terms[ii+1][jj][kk] == YES ||
                   re_NS_terms[ii][jj+1][kk] == YES ||
                   re_NS_terms[ii+1][jj+1][kk] == YES ||
                   re_NS_terms[ii][jj][kk+1] == YES ||
                   re_NS_terms[ii+1][jj][kk+1] == YES ||
                   re_NS_terms[ii][jj+1][kk+1] == YES ||
                   re_NS_terms[ii+1][jj+1][kk+1] == YES)
                {                
                    co_NS_terms[i][j][k] = YES;

                    co_coords_x[i][j][k] = re_coords_x[ii][jj][kk] + (dh[0]/2);
                    co_coords_y[i][j][k] = re_coords_y[ii][jj][kk] + (dh[1]/2);
                    co_coords_z[i][j][k] = re_coords_z[ii][jj][kk] + (dh[2]/2);

                    /*
                    sum_co_rho=0.0,sum_co_vel_u=0.0,sum_co_vel_v=0.0,sum_co_vel_w=0.0;
                    sum_co_rho_vel_u=0.0,sum_co_rho_vel_v=0.0,sum_co_rho_vel_w=0.0;
                    sum_co_rho_vel_uu=0.0,sum_co_rho_vel_vv=0.0,sum_co_rho_vel_ww=0.0;
                    sum_co_rho_vel_uv=0.0,sum_co_rho_vel_vw=0.0,sum_co_rho_vel_uw=0.0;
                    sum_co_rho_cp=0.0,sum_co_cp=0.0,sum_co_rho_t=0.0,sum_co_rho_c=0.0;
                    sum_co_rho_c0=0.0,sum_co_t=0.0,sum_co_con=0.0,sum_co_con0=0.0;
                    sum_co_rho_cu=0.0,sum_co_rho_cv=0.0,sum_co_rho_cw=0.0;
                    sum_co_rho_cu0=0.0,sum_co_rho_cv0=0.0,sum_co_rho_cw0=0.0;
                    sum_co_rho_cp_tu=0.0,sum_co_rho_cp_tv=0.0,sum_co_rho_cp_tw=0.0;
                    sum_co_temp11=0.0,sum_co_temp12=0.0,sum_co_temp13=0.0;
                    sum_co_temp21=0.0,sum_co_temp22=0.0,sum_co_temp23=0.0;
                    sum_co_temp31=0.0,sum_co_temp32=0.0,sum_co_temp33=0.0;
                    sum_co_temp=0.0,sum_co_temp_t1=0.0,sum_co_temp_t2=0.0;
                    sum_co_temp_t3=0.0,sum_co_temp_c1=0.0,sum_co_temp_c2=0.0;
                    sum_co_temp_c3=0.0,sum_co_temp_c10=0.0,sum_co_temp_c20=0.0;
                    sum_co_temp_c30=0.0;
                    */

                    sum_co_rho=0.0,sum_co_vel_u=0.0,sum_co_vel_v=0.0,sum_co_vel_w=0.0;
                    sum_co_rho_vel_u=0.0,sum_co_rho_vel_v=0.0,sum_co_rho_vel_w=0.0;
                    sum_co_rho_vel_uu=0.0,sum_co_rho_vel_vv=0.0,sum_co_rho_vel_ww=0.0;
                    sum_co_rho_vel_uv=0.0,sum_co_rho_vel_vw=0.0,sum_co_rho_vel_uw=0.0;
                    sum_co_rho_cp=0.0,sum_co_cp=0.0,sum_co_rho_t=0.0,sum_co_t=0.0;
                    sum_co_rho_cp_tu=0.0,sum_co_rho_cp_tv=0.0,sum_co_rho_cp_tw=0.0;
                    sum_co_temp11=0.0,sum_co_temp12=0.0,sum_co_temp13=0.0;
                    sum_co_temp21=0.0,sum_co_temp22=0.0,sum_co_temp23=0.0;
                    sum_co_temp31=0.0,sum_co_temp32=0.0,sum_co_temp33=0.0;
                    sum_co_temp=0.0,sum_co_temp_t1=0.0,sum_co_temp_t2=0.0,sum_co_temp_t3=0.0;
                    for (ic = 0; ic < MAX_NCOMPS; ic++)
                    {
                    sum_co_rho_c[ic] = 0.0;
                    sum_co_con[ic] = 0.0;
                    sum_co_rho_cu[ic] = 0.0;
                    sum_co_rho_cv[ic] = 0.0;
                    sum_co_rho_cw[ic] = 0.0;
                    sum_co_temp_c1[ic] = 0.0;
                    sum_co_temp_c2[ic] = 0.0;
                    sum_co_temp_c3[ic] = 0.0;
                    }

                    for(kkk = kk; kkk < kk+2; kkk++)
                    {
                        for(jjj = jj; jjj < jj+2; jjj++)
                        {
                           for(iii = ii; iii < ii+2; iii++) 
                           {  
                               sum_co_rho += re_rho[iii][jjj][kkk];
                               sum_co_vel_u += re_vel_u[iii][jjj][kkk];
                               sum_co_vel_v += re_vel_v[iii][jjj][kkk];
                               sum_co_vel_w += re_vel_w[iii][jjj][kkk];
                               sum_co_rho_vel_u += re_rho_vel_u[iii][jjj][kkk];
                               sum_co_rho_vel_v += re_rho_vel_v[iii][jjj][kkk];
                               sum_co_rho_vel_w += re_rho_vel_w[iii][jjj][kkk];
                               sum_co_rho_vel_uu += re_rho_vel_uu[iii][jjj][kkk];
                               sum_co_rho_vel_vv += re_rho_vel_vv[iii][jjj][kkk];
                               sum_co_rho_vel_ww += re_rho_vel_ww[iii][jjj][kkk];
                               sum_co_rho_vel_uv += re_rho_vel_uv[iii][jjj][kkk];
                               sum_co_rho_vel_vw += re_rho_vel_vw[iii][jjj][kkk];
                               sum_co_rho_vel_uw += re_rho_vel_uw[iii][jjj][kkk];
                               sum_co_rho_cp += re_rho_cp[iii][jjj][kkk];
                               sum_co_cp += re_cp[iii][jjj][kkk];
                               sum_co_rho_t += re_rho_t[iii][jjj][kkk];
                               sum_co_t += re_t[iii][jjj][kkk];
                               sum_co_rho_cp_tu += re_rho_cp_tu[iii][jjj][kkk];
                               sum_co_rho_cp_tv += re_rho_cp_tv[iii][jjj][kkk];
                               sum_co_rho_cp_tw += re_rho_cp_tw[iii][jjj][kkk];
                               sum_co_temp11 += re_temp11[iii][jjj][kkk];
                               sum_co_temp12 += re_temp12[iii][jjj][kkk];
                               sum_co_temp13 += re_temp13[iii][jjj][kkk];
                               sum_co_temp21 += re_temp21[iii][jjj][kkk];
                               sum_co_temp22 += re_temp22[iii][jjj][kkk];
                               sum_co_temp23 += re_temp23[iii][jjj][kkk];
                               sum_co_temp31 += re_temp31[iii][jjj][kkk];
                               sum_co_temp32 += re_temp32[iii][jjj][kkk];
                               sum_co_temp33 += re_temp33[iii][jjj][kkk];
                               sum_co_temp += re_temp[iii][jjj][kkk];
                               sum_co_temp_t1 += re_temp_t1[iii][jjj][kkk];
                               sum_co_temp_t2 += re_temp_t2[iii][jjj][kkk];
                               sum_co_temp_t3 += re_temp_t3[iii][jjj][kkk];
                               if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                               {
                                   for (ic = 0; ic < MAX_NCOMPS; ic++)
                                   {
                                       sum_co_rho_c[ic] += re_rho_c[iii][jjj][kkk][ic];
                                       sum_co_con[ic] += re_con[iii][jjj][kkk][ic];
                                       sum_co_rho_cu[ic] += re_rho_cu[iii][jjj][kkk][ic];
                                       sum_co_rho_cv[ic] += re_rho_cv[iii][jjj][kkk][ic];
                                       sum_co_rho_cw[ic] += re_rho_cw[iii][jjj][kkk][ic];
                                       sum_co_temp_c1[ic] += re_temp_c1[iii][jjj][kkk][ic];
                                       sum_co_temp_c2[ic] += re_temp_c2[iii][jjj][kkk][ic];
                                       sum_co_temp_c3[ic] += re_temp_c3[iii][jjj][kkk][ic];
                                   }
                               }
                            }
                         }
                     }

                     co_rho[i][j][k] =  sum_co_rho/8.0;
                     co_vel_u[i][j][k] = sum_co_vel_u/8.0;
                     co_vel_v[i][j][k] = sum_co_vel_v/8.0;
                     co_vel_w[i][j][k] = sum_co_vel_w/8.0;
                     co_rho_vel_u[i][j][k] = sum_co_rho_vel_u/8.0;
                     co_rho_vel_v[i][j][k] = sum_co_rho_vel_v/8.0;
                     co_rho_vel_w[i][j][k] = sum_co_rho_vel_w/8.0;
                     co_rho_vel_uu[i][j][k] = sum_co_rho_vel_uu/8.0;
                     co_rho_vel_vv[i][j][k] = sum_co_rho_vel_vv/8.0;
                     co_rho_vel_ww[i][j][k] = sum_co_rho_vel_ww/8.0;
                     co_rho_vel_uv[i][j][k] = sum_co_rho_vel_uv/8.0;
                     co_rho_vel_vw[i][j][k] = sum_co_rho_vel_vw/8.0;
                     co_rho_vel_uw[i][j][k] = sum_co_rho_vel_uw/8.0;
                     co_rho_cp[i][j][k] = sum_co_rho_cp/8.0;
                     co_cp[i][j][k] = sum_co_cp/8.0;
                     co_rho_t[i][j][k] = sum_co_rho_t/8.0;
                     co_t[i][j][k] = sum_co_t/8.0; 
                     co_rho_cp_tu[i][j][k] = sum_co_rho_cp_tu/8.0;
                     co_rho_cp_tv[i][j][k] = sum_co_rho_cp_tv/8.0;
                     co_rho_cp_tw[i][j][k] = sum_co_rho_cp_tw/8.0;
                     co_temp11[i][j][k] = sum_co_temp11/8.0;
                     co_temp12[i][j][k] = sum_co_temp12/8.0;
                     co_temp13[i][j][k] = sum_co_temp13/8.0;
                     co_temp21[i][j][k] = sum_co_temp21/8.0;
                     co_temp22[i][j][k] = sum_co_temp22/8.0;
                     co_temp23[i][j][k] = sum_co_temp23/8.0;
                     co_temp31[i][j][k] = sum_co_temp31/8.0;
                     co_temp32[i][j][k] = sum_co_temp32/8.0;
                     co_temp33[i][j][k] = sum_co_temp33/8.0;
                     co_temp[i][j][k] = sum_co_temp/8.0;
                     co_temp_t1[i][j][k] = sum_co_temp_t1/8.0;
                     co_temp_t2[i][j][k] = sum_co_temp_t2/8.0;
                     co_temp_t3[i][j][k] = sum_co_temp_t3/8.0;
                     if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                     {
                     for (ic = 0; ic < MAX_NCOMPS; ic++)
                     {
                     co_rho_c[i][j][k][ic] = sum_co_rho_c[ic]/8.0;
                     co_con[i][j][k][ic] = sum_co_con[ic]/8.0;
                     co_rho_cu[i][j][k][ic] = sum_co_rho_cu[ic]/8.0;
                     co_rho_cv[i][j][k][ic] = sum_co_rho_cv[ic]/8.0;
                     co_rho_cw[i][j][k][ic] = sum_co_rho_cw[ic]/8.0;
                     co_temp_c1[i][j][k][ic] = sum_co_temp_c1[ic]/8.0;
                     co_temp_c2[i][j][k][ic] = sum_co_temp_c2[ic]/8.0;
                     co_temp_c3[i][j][k][ic] = sum_co_temp_c3[ic]/8.0;
                     }
                     }
                } 
                else
                {
                     co_NS_terms[i][j][k] = NO;
                     co_rho[i][j][k] = 12.41e-05;
                }
            }
        }
        }
	
 /*       for (k = imin[2]; k < (imax[2]/2)+2; k++)
        {
        for (j = imin[1]; j < (imax[1]/2)+2; j++)
        {
            for (i = imin[0]; i < (imax[0]/2)+2; i++)
            {
                if(co_NS_terms[i][j][k] == YES)
                {
                    //spherical average
                    //r[i][j][k] = sqrt( (co_coords_x[i][j][k]*co_coords_x[i][j][k])
                    //                 + (co_coords_y[i][j][k]*co_coords_y[i][j][k]) 
                    //               + (co_coords_z[i][j][k]*co_coords_z[i][j][k]) );

                    //planar average
                    r[i][j][k] = co_coords_z[i][j][k];

                    r_color[i][j][k] = (int)(r[i][j][k]/delta_r);
                    r_color_LA[i][j][k] = (int)(r[i][j][k]/delta_r);
                    r_color_LI[i][j][k] = (int)(r[i][j][k]/delta_r);
                    r_color_LH[i][j][k] = (int)(r[i][j][k]/delta_r);
                    r_color_LC[i][j][k] = (int)(r[i][j][k]/delta_r);
                    r_color_LC0[i][j][k] = (int)(r[i][j][k]/delta_r);
               }
               else
               {
                    r_color[i][j][k] = -1;
                    r_color_LA[i][j][k] = -1;
                    r_color_LI[i][j][k] = -1;
                    r_color_LH[i][j][k] = -1;
                    r_color_LC[i][j][k] = -1;
                    r_color_LC0[i][j][k] = -1;
               }
            }
        }
        }           */

        for (k = imin[2]+1; k < (imax[2]/2)+1; k++)
        {
        for (j = imin[1]+1; j < (imax[1]/2)+1; j++)
        {
            for (i = imin[0]+1; i < (imax[0]/2)+1; i++)
            {
                if(co_NS_terms[i][j][k] == YES)
                {
                    if(ppz == ppgmax[2]-1 && k == (imax[2]/2))
                         co_NS_terms[i][j][k]  = NO;
                    else if(ppz == 0 && k == imin[2]+1)
                         co_NS_terms[i][j][k]  = NO;
                    else
                    {

                    hat_ux = ((co_rho_vel_u[i+1][j][k]/co_rho[i+1][j][k])
                                  - (co_rho_vel_u[i-1][j][k]/co_rho[i-1][j][k])) / (4.0*dh[0]);
                    hat_uy = ((co_rho_vel_u[i][j+1][k]/co_rho[i][j+1][k])
                                  - (co_rho_vel_u[i][j-1][k]/co_rho[i][j-1][k])) / (4.0*dh[1]);
		    hat_uz = ((co_rho_vel_u[i][j][k+1]/co_rho[i][j][k+1])
                                  - (co_rho_vel_u[i][j][k-1]/co_rho[i][j][k-1])) / (4.0*dh[2]);

                    hat_vx = ((co_rho_vel_v[i+1][j][k]/co_rho[i+1][j][k])
                                  - (co_rho_vel_v[i-1][j][k]/co_rho[i-1][j][k])) / (4.0*dh[0]);
                    hat_vy = ((co_rho_vel_v[i][j+1][k]/co_rho[i][j+1][k])
                                  - (co_rho_vel_v[i][j-1][k]/co_rho[i][j-1][k])) / (4.0*dh[1]);
                    hat_vz = ((co_rho_vel_v[i][j][k+1]/co_rho[i][j][k+1])
                                  - (co_rho_vel_v[i][j][k-1]/co_rho[i][j][k-1])) / (4.0*dh[2]);

                    hat_wx = ((co_rho_vel_w[i+1][j][k]/co_rho[i+1][j][k])
                                  - (co_rho_vel_w[i-1][j][k]/co_rho[i-1][j][k])) / (4.0*dh[0]);
                    hat_wy = ((co_rho_vel_w[i][j+1][k]/co_rho[i][j+1][k])
                                  - (co_rho_vel_w[i][j-1][k]/co_rho[i][j-1][k])) / (4.0*dh[1]);
                    hat_wz = ((co_rho_vel_w[i][j][k+1]/co_rho[i][j][k+1])
                                  - (co_rho_vel_w[i][j][k-1]/co_rho[i][j][k-1])) / (4.0*dh[2]);


                    hat_tx  = (co_t[i+1][j][k] - co_t[i-1][j][k]) / (4.0*dh[0]);
                    hat_ty  = (co_t[i][j+1][k] - co_t[i][j-1][k]) / (4.0*dh[1]);
                    hat_tz  = (co_t[i][j][k+1] - co_t[i][j][k-1]) / (4.0*dh[2]);

                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    for (ic = 0; ic < MAX_NCOMPS; ic++)
                    {
                    hat_cx[ic]  = (co_con[i+1][j][k][ic] - co_con[i-1][j][k][ic]) / (4.0*dh[0]);
                    hat_cy[ic]  = (co_con[i][j+1][k][ic] - co_con[i][j-1][k][ic]) / (4.0*dh[1]);
                    hat_cz[ic]  = (co_con[i][j][k+1][ic] - co_con[i][j][k-1][ic]) / (4.0*dh[2]);
                    }
                    }

                    hat_S11 = hat_ux;
                    hat_S12 = (hat_uy + hat_vx)/2.0;
                    hat_S13 = (hat_uz + hat_wx)/2.0;
                    hat_S21 = (hat_vx + hat_uy)/2.0;
                    hat_S22 = hat_vy;
                    hat_S23 = (hat_vz + hat_wy)/2.0;
                    hat_S31 = (hat_wx + hat_uz)/2.0;
                    hat_S32 = (hat_wy + hat_vz)/2.0;
                    hat_S33 = hat_wz;

                    hat_SS = sqrt(2.0*( (hat_S11*hat_S11) + (hat_S12*hat_S12)
                                    + (hat_S13*hat_S13) + (hat_S21*hat_S21)
                                    + (hat_S22*hat_S22) + (hat_S23*hat_S23)
                                    + (hat_S31*hat_S31) + (hat_S32*hat_S32)
                                    + (hat_S33*hat_S33)));

                    hat_S11_a = hat_S11 - ((hat_S11+hat_S22+hat_S33)/3.0);
                    hat_S12_a = hat_S12;
                    hat_S13_a = hat_S13;
                    hat_S21_a = hat_S21;
                    hat_S22_a = hat_S22 - ((hat_S11+hat_S22+hat_S33)/3.0);
                    hat_S23_a = hat_S23;
                    hat_S31_a = hat_S31;
                    hat_S32_a = hat_S32;
                    hat_S33_a = hat_S33 - ((hat_S11+hat_S22+hat_S33)/3.0);


                    MA11[i][j][k] = 2.0*delta*delta*co_temp11[i][j][k]
                                 - (2.0*co_delta*co_delta*co_rho[i][j][k]*hat_SS
                                 *(hat_S11_a));
                    MA12[i][j][k] = 2.0*delta*delta*co_temp12[i][j][k]
                                 - (2.0*co_delta*co_delta*co_rho[i][j][k]*hat_SS
                                   *hat_S12_a);
                    MA13[i][j][k] = 2.0*delta*delta*co_temp13[i][j][k]
                                 - (2.0*co_delta*co_delta*co_rho[i][j][k]*hat_SS
                                   *hat_S13_a);
                    MA22[i][j][k] = 2.0*delta*delta*co_temp22[i][j][k]
                                 - (2.0*co_delta*co_delta*co_rho[i][j][k]*hat_SS
                                   *hat_S22_a);
                    MA23[i][j][k] = 2.0*delta*delta*co_temp23[i][j][k]
                                 - (2.0*co_delta*co_delta*co_rho[i][j][k]*hat_SS
                                   *hat_S23_a);


                    LA11[i][j][k] = co_rho_vel_uu[i][j][k]
                                 - (co_rho_vel_u[i][j][k]*co_rho_vel_u[i][j][k]/co_rho[i][j][k]);
                    LA12[i][j][k] = co_rho_vel_uv[i][j][k]
                                 - (co_rho_vel_u[i][j][k]*co_rho_vel_v[i][j][k]/co_rho[i][j][k]);
                    LA13[i][j][k] = co_rho_vel_uw[i][j][k]
                                 - (co_rho_vel_u[i][j][k]*co_rho_vel_w[i][j][k]/co_rho[i][j][k]);
                    LA22[i][j][k] = co_rho_vel_vv[i][j][k]
                                 - (co_rho_vel_v[i][j][k]*co_rho_vel_v[i][j][k]/co_rho[i][j][k]);
                    LA23[i][j][k] = co_rho_vel_vw[i][j][k]
                                 - (co_rho_vel_v[i][j][k]*co_rho_vel_w[i][j][k]/co_rho[i][j][k]);
                    LA33[i][j][k] = co_rho_vel_ww[i][j][k]
                                 - (co_rho_vel_w[i][j][k]*co_rho_vel_w[i][j][k]/co_rho[i][j][k]);

		    LA1[i][j][k] = LA11[i][j][k] - ((LA11[i][j][k]+LA22[i][j][k]+LA33[i][j][k])/3.0);
		    LA2[i][j][k] = LA12[i][j][k];
                    LA3[i][j][k] = LA13[i][j][k];
                    LA4[i][j][k] = LA22[i][j][k] - ((LA11[i][j][k]+LA22[i][j][k]+LA33[i][j][k])/3.0);
                    LA5[i][j][k] = LA23[i][j][k];

		    LI[i][j][k] = (LA11[i][j][k]+LA22[i][j][k]+LA33[i][j][k])/3.0;

                    MI[i][j][k] = -(delta*delta*co_temp[i][j][k])
		               + (co_delta*co_delta*co_rho[i][j][k]*hat_SS*hat_SS);
             
                    MH1[i][j][k] = delta*delta*co_temp_t1[i][j][k]
                                - (co_delta*co_delta*co_rho[i][j][k]*co_cp[i][j][k]*hat_SS*hat_tx);
                    MH2[i][j][k] = delta*delta*co_temp_t2[i][j][k]
                                - (co_delta*co_delta*co_rho[i][j][k]*co_cp[i][j][k]*hat_SS*hat_ty);
                    MH3[i][j][k] = delta*delta*co_temp_t3[i][j][k]
                                - (co_delta*co_delta*co_rho[i][j][k]*co_cp[i][j][k]*hat_SS*hat_tz);
                    MMA[i][j][k] = (MA11[i][j][k]*MA11[i][j][k]) + (MA12[i][j][k]*MA12[i][j][k]) + (MA13[i][j][k]*MA13[i][j][k]) + (MA22[i][j][k]*MA22[i][j][k]) + (MA23[i][j][k]*MA23[i][j][k]);

                    LMA[i][j][k] = (LA1[i][j][k]*MA11[i][j][k]) + (LA2[i][j][k]*MA12[i][j][k]) + (LA3[i][j][k]*MA13[i][j][k]) + (LA4[i][j][k]*MA22[i][j][k]) + (LA5[i][j][k]*MA23[i][j][k]);

                    LH1[i][j][k] = co_rho_cp_tu[i][j][k]
                                - ((co_rho_cp[i][j][k]*co_rho_t[i][j][k]*co_rho_vel_u[i][j][k])/(co_rho[i][j][k]*co_rho[i][j][k]));
                    LH2[i][j][k] = co_rho_cp_tv[i][j][k]
                                - ((co_rho_cp[i][j][k]*co_rho_t[i][j][k]*co_rho_vel_v[i][j][k])/(co_rho[i][j][k]*co_rho[i][j][k]));
                    LH3[i][j][k] = co_rho_cp_tw[i][j][k]
                                - ((co_rho_cp[i][j][k]*co_rho_t[i][j][k]*co_rho_vel_w[i][j][k])/(co_rho[i][j][k]*co_rho[i][j][k]));
                    LMH[i][j][k] = (LH1[i][j][k]*MH1[i][j][k]) + (LH2[i][j][k]*MH2[i][j][k]) + (LH3[i][j][k]*MH3[i][j][k]);

                    MMH[i][j][k] = (MH1[i][j][k]*MH1[i][j][k]) + (MH2[i][j][k]*MH2[i][j][k]) + (MH3[i][j][k]*MH3[i][j][k]);

                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                    for (ic = 0; ic < MAX_NCOMPS; ic++)
                    {
                    MC1[i][j][k][ic] = delta*delta*co_temp_c1[i][j][k][ic]
                                - (co_delta*co_delta*co_rho[i][j][k]*hat_SS*hat_cx[ic]);
                    MC2[i][j][k][ic] = delta*delta*co_temp_c2[i][j][k][ic]
                                - (co_delta*co_delta*co_rho[i][j][k]*hat_SS*hat_cy[ic]);
                    MC3[i][j][k][ic] = delta*delta*co_temp_c3[i][j][k][ic]
                                - (co_delta*co_delta*co_rho[i][j][k]*hat_SS*hat_cz[ic]);
                    LC1[i][j][k][ic] = co_rho_cu[i][j][k][ic]
                                - ((co_rho_vel_u[i][j][k]*co_rho_c[i][j][k][ic])/co_rho[i][j][k]);
                    LC2[i][j][k][ic] = co_rho_cv[i][j][k][ic]
                                - ((co_rho_vel_v[i][j][k]*co_rho_c[i][j][k][ic])/co_rho[i][j][k]);
                    LC3[i][j][k][ic] = co_rho_cw[i][j][k][ic]
                                - ((co_rho_vel_w[i][j][k]*co_rho_c[i][j][k][ic])/co_rho[i][j][k]);
                    LMC[i][j][k][ic] = (LC1[i][j][k][ic]*MC1[i][j][k][ic]) + (LC2[i][j][k][ic]*MC2[i][j][k][ic]) + (LC3[i][j][k][ic]*MC3[i][j][k][ic]);

                    MMC[i][j][k][ic] = (MC1[i][j][k][ic]*MC1[i][j][k][ic]) + (MC2[i][j][k][ic]*MC2[i][j][k][ic]) + (MC3[i][j][k][ic]*MC3[i][j][k][ic]);
                    }
                    }

  /*                      if(((LA1[i][j][k]*MA11[i][j][k]) + (LA2[i][j][k]*MA12[i][j][k])
                            + (LA3[i][j][k]*MA13[i][j][k]) + (LA4[i][j][k]*MA22[i][j][k])
                            + (LA5[i][j][k]*MA23[i][j][k])) < 0.0)
                        {
                            r_color_LA[i][j][k] = -2;
                        }
                        if((LI[i][j][k] < 0.0 && MI[i][j][k] > 0.0) || (LI[i][j][k] > 0.0 && MI[i][j][k] < 0.0))
                        {
                            r_color_LI[i][j][k] = -2;
                        }
                        if((LH1[i][j][k]*MH1[i][j][k]) + (LH2[i][j][k]*MH2[i][j][k]) + (LH3[i][j][k]*MH3[i][j][k]) < 0.0)
                        {
                            r_color_LH[i][j][k] = -2;
                        }
                        if((LC1[i][j][k]*MC1[i][j][k]) + (LC2[i][j][k]*MC2[i][j][k]) + (LC3[i][j][k]*MC3[i][j][k]) < 0.0)
                        {
                            r_color_LC[i][j][k] = -2;
                        }
                        if((LC10[i][j][k]*MC10[i][j][k]) + (LC20[i][j][k]*MC20[i][j][k]) + (LC30[i][j][k]*MC30[i][j][k]) < 0.0)
                        {
                            r_color_LC0[i][j][k] = -2;
                        }         */
                    }
               }
            }
        }
        }

/*        for (k = 0; k < num_r; k++)
        {
              deno[k] = 0.0;
              nume[k] = 0.0;
              ci_deno[k] = 0.0;
              ci_nume[k] = 0.0;
              prt_deno[k] = 0.0;
              prt_nume[k] = 0.0;
              sct_deno[k] = 0.0;
              sct_nume[k] = 0.0;
              sct_deno0[k] = 0.0;
              sct_nume0[k] = 0.0;

	      cs_new[k] = 0.0;
              ci_new[k] = 0.0;
              prt_new[k] = 0.0;
              sct_new[k] = 0.0;
              sct0_new[k] = 0.0;
        }      

        for (kkk = 0; kkk < num_r; kkk++)
        {
        for (k = imin[2]+1; k < (imax[2]/2)+1; k++)
        {
           for (j = imin[1]+1; j < (imax[1]/2)+1; j++)
           {
               for (i = imin[0]+1; i < (imax[0]/2)+1; i++)
               {
                    if(co_NS_terms[i][j][k] == YES)
                    {
                       if(kkk == r_color_LA[i][j][k])
                       {
                          deno[kkk] += ((MA11[i][j][k]*MA11[i][j][k])
                                     + (MA12[i][j][k]*MA12[i][j][k])
                                     + (MA13[i][j][k]*MA13[i][j][k]) 
                                     + (MA22[i][j][k]*MA22[i][j][k])
                                     + (MA23[i][j][k]*MA23[i][j][k]));
                          nume[kkk] += ((LA1[i][j][k]*MA11[i][j][k]) 
                                     + (LA2[i][j][k]*MA12[i][j][k])
                                     + (LA3[i][j][k]*MA13[i][j][k]) 
                                     + (LA4[i][j][k]*MA22[i][j][k])
                                     + (LA5[i][j][k]*MA23[i][j][k]));
                       }
                       if(kkk == r_color_LI[i][j][k])
                       {
                          ci_deno[kkk] += MI[i][j][k];
                          ci_nume[kkk] += LI[i][j][k];
                       }
                       if(kkk == r_color_LH[i][j][k])
                       {
                          prt_nume[kkk] += (LH1[i][j][k]*MH1[i][j][k]) + (LH2[i][j][k]*MH2[i][j][k]) + (LH3[i][j][k]*MH3[i][j][k]);
                          prt_deno[kkk] += (MH1[i][j][k]*MH1[i][j][k]) + (MH2[i][j][k]*MH2[i][j][k]) + (MH3[i][j][k]*MH3[i][j][k]);
                       }

                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                       if(kkk == r_color_LC[i][j][k])
                       {
                          sct_nume[kkk] += (LC1[i][j][k]*MC1[i][j][k]) + (LC2[i][j][k]*MC2[i][j][k]) + (LC3[i][j][k]*MC3[i][j][k]);
                          sct_deno[kkk] += (MC1[i][j][k]*MC1[i][j][k]) + (MC2[i][j][k]*MC2[i][j][k]) + (MC3[i][j][k]*MC3[i][j][k]);
                       }
                       if(kkk == r_color_LC0[i][j][k])
                       {
                          sct_nume0[kkk] += (LC10[i][j][k]*MC10[i][j][k]) + (LC20[i][j][k]*MC20[i][j][k]) + (LC30[i][j][k]*MC30[i][j][k]);
                          sct_deno0[kkk] += (MC10[i][j][k]*MC10[i][j][k]) + (MC20[i][j][k]*MC20[i][j][k]) + (MC30[i][j][k]*MC30[i][j][k]);
                       }
                    }
                    }
                  }
             }
        }
        }   Why specify if(co_NS_terms[i][j][k] == YES)?   */

        //if(co_NS_terms[i][j][k] == YES)
 	mesh_avr(wv,MMA, deno, iimin, iimax, av);
       	mesh_avr(wv,LMA, nume, iimin, iimax, av);   
       	mesh_avr(wv,MI, ci_deno, iimin, iimax, av);
       	mesh_avr(wv,LI, ci_nume, iimin, iimax, av);
       	mesh_avr(wv,LMH, prt_nume, iimin, iimax, av);
       	mesh_avr(wv,MMH, prt_deno, iimin, iimax, av);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
             mesh_avr_mc(wv,LMC, sct_nume, iimin, iimax, av);
             mesh_avr_mc(wv,MMC, sct_deno, iimin, iimax, av);
        }
      	
/*        if (nn > 1)
        {
            for (k = 0; k < num_r; k++)
            {
                 pp_global_sum(&deno[k],1L);
                 pp_global_sum(&nume[k],1L);
                 pp_global_sum(&ci_deno[k],1L);
                 pp_global_sum(&ci_nume[k],1L);
                 pp_global_sum(&prt_deno[k],1L);
                 pp_global_sum(&prt_nume[k],1L);
                 pp_global_sum(&sct_deno[k],1L);
                 pp_global_sum(&sct_nume[k],1L);
                 pp_global_sum(&sct_deno0[k],1L);
                 pp_global_sum(&sct_nume0[k],1L);
            }
        }     

        for (k = 0; k < num_r; k++)
        {
                if(deno[k] < 10e-16)
                cs_new[k] = 0.0;
                else
                cs_new[k] = nume[k]/deno[k];
                if(ci_deno[k] < 10e-16)
                ci_new[k] = 0.0;
                else
                ci_new[k] = ci_nume[k]/ci_deno[k];
                if(prt_deno[k] < 10e-16)
                prt_new[k] = 0.0;
                else
                prt_new[k] = prt_nume[k]/prt_deno[k];

                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                if(sct_deno[k] < 10e-16)
                sct_new[k] = 0.0;
                else
                sct_new[k] = sct_nume[k]/sct_deno[k];
                if(sct_deno0[k] < 10e-16)
                sct0_new[k] = 0.0;
                else
                sct0_new[k] = sct_nume0[k]/sct_deno0[k];
                }

                if(fabs(cs_new[k])<10e-10 || cs_new[k]<0.0)
                {
                   cs_new[k] = 0.0;
                   ci_new[k] = 0.0;
                   prt_new[k] = 0.0;
                   sct_new[k] = 0.0;
                   sct0_new[k] = 0.0;
                }

                if(fabs(ci_new[k])<10e-10 || ci_new[k]<0.0)
                    ci_new[k] = 0.0;
                if(fabs(prt_new[k])<10e-10 || prt_new[k]<0.0)
                    prt_new[k] = 0.0;
                if(fabs(sct_new[k])<10e-10 || sct_new[k]<0.0)
                    sct_new[k] = 0.0;
                if(fabs(sct0_new[k])<10e-10 || sct0_new[k]<0.0)
                    sct0_new[k] = 0.0;
        }  */

        for (k = 0; k < avmax[2] - avmin[2]; k++)
        {
            for (j = 0; j < avmax[1] - avmin[1]; j++)
            {
                for (i = 0; i < avmax[0] - avmin[0]; i++)
                {
                    if(deno[i][j][k] < 10e-16)
                    cs_new[i][j][k] = 0.0;
                    else
                    cs_new[i][j][k] = nume[i][j][k]/deno[i][j][k];
                    if(ci_deno[i][j][k] < 10e-16)
                    ci_new[i][j][k] = 0.0;
                    else
                    ci_new[i][j][k] = ci_nume[i][j][k]/ci_deno[i][j][k];
                    if(prt_deno[i][j][k] < 10e-16)
                    prt_new[i][j][k] = 0.0;
                    else
                    prt_new[i][j][k] = prt_nume[i][j][k]/prt_deno[i][j][k];
            
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                        for (ic = 0; ic < MAX_NCOMPS; ic++)
                        {
                        if(sct_deno[i][j][k][ic] < 10e-16)
                        sct_new[i][j][k][ic] = 0.0;
                        else
                        sct_new[i][j][k][ic] = sct_nume[i][j][k][ic]/sct_deno[i][j][k][ic];
                        }
                    }
        
		   if(fabs(cs_new[i][j][k]) < 10e-10 || cs_new[i][j][k] < 0.0)
                        cs_new[i][j][k] = 0.0;
                  if(fabs(ci_new[i][j][k]) < 10e-10 || ci_new[i][j][k] < 0.0)
                  	ci_new[i][j][k] = 0.0;
                  if(fabs(prt_new[i][j][k]) < 10e-10 || prt_new[i][j][k] < 0.0)
                  	prt_new[i][j][k] = 0.0;
                  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                  for (ic = 0; ic < MAX_NCOMPS; ic++)
                      if(fabs(sct_new[i][j][k][ic]) < 10e-10 || sct_new[i][j][k][ic] < 0.0)
                          sct_new[i][j][k][ic] = 0.0;		  
                }
            }
        }       
	
        for (i = imin[0]; i < imin[0] + av[0]/2 + 1; i++)
        {
        for (j = imin[1]; j < ((imax[1] - imin[1])/2); j++)
        {
	    if (j >= imin[1] && j <= imin[1] + av[1]/2)
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                 if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {
                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
                              co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[0][0][0];
                       co_CI[i][j][k] = ci_new[0][0][0];
                       co_Prt[i][j][k] = prt_new[0][0][0];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[0][0][0][ic];
		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              { 
                       co_C[i][j][k] = cs_new[0][0][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[0][0][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[0][0][avmax[2]-avmin[2]-1];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                            co_Sct[i][j][k][ic] = sct_new[0][0][avmax[2]-avmin[2]-1][ic];
		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[0][0][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[0][0][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[0][0][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[0][0][k-imin[2]-av[2]/2][ic];
		      }
                  }
             } 
	    }
            else if (j >= ((imax[1]-imin[1])/2) - 1 - av[1]/2 && j <= ((imax[1]-imin[1])/2) - 1)
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[0][avmax[1]-avmin[1]-1][0];
                       co_CI[i][j][k] = ci_new[0][avmax[1]-avmin[1]-1][0];
                       co_Prt[i][j][k] = prt_new[0][avmax[1]-avmin[1]-1][0];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[0][avmax[1]-avmin[1]-1][0][ic];
		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
			  
                       co_C[i][j][k] = cs_new[0][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[0][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[0][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                        co_Sct[i][j][k][ic] = sct_new[0][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1][ic];
		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[0][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[0][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[0][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[0][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2][ic];
		      }

                  }
             }

	    }
	    else
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[0][j-imin[1]-av[1]/2][0];
                       co_CI[i][j][k] = ci_new[0][j-imin[1]-av[1]/2][0];
                       co_Prt[i][j][k] = prt_new[0][j-imin[1]-av[1]/2][0];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[0][j-imin[1]-av[1]/2][0][ic];

		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
                       co_C[i][j][k] = cs_new[0][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[0][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[0][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                            co_Sct[i][j][k][ic] = sct_new[0][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1][ic];
		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[0][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[0][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[0][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[0][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2][ic];

		      }

                  }
              }
	    }
        }
        }








        for (i = imin[0] + av[0]/2 + 1; i < ((imax[0]-imin[0])/2) - 1 - av[0]/2; i++)
        {
        for (j = imin[1]; j < ((imax[1]-imin[1])/2); j++)
        {
	    if (j >= imin[1] && j <= imin[1] + av[1]/2)
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][0][0];
                       co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][0][0];
                       co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][0][0];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][0][0][ic];
		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
                       co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][0][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][0][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][0][avmax[2]-avmin[2]-1];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][0][avmax[2]-avmin[2]-1][ic];
		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][0][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][0][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][0][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][0][k-imin[2]-av[2]/2][ic];

		      }

                  }
             } 
	    }
            else if (j >= ((imax[1]-imin[1])/2) - 1 - av[1]/2 && j <= ((imax[1]-imin[1])/2) - 1)
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][0];
                       co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][0];
                       co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][0];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][0][ic];
		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
                       co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                            co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1][ic];

		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2][ic];

		      }

                  }
             }

	    }
	    else
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)                    
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][0];
                       co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][0];
                       co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][0];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][0][ic];

		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
                       co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                            co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1][ic];
		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[i-imin[0]-av[0]/2][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2][ic];

		      }

                  }
              }
	    }
        }
        }

        for (i = (imax[0]-imin[0])/2 - 1 - av[0]/2; i < (imax[0]-imin[0])/2; i++)
        {
        for (j = imin[1]; j < ((imax[1]-imin[1])/2); j++)
        {
	    if (j >= imin[1] && j <= imin[1] + av[1]/2)
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][0][0];
                       co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][0][0];
                       co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][0][0];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][0][0][ic];
		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
			  
                       co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][0][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][0][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][0][avmax[2]-avmin[2]-1];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][0][avmax[2]-avmin[2]-1][ic];
		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][0][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][0][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][0][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][0][k-imin[2]-av[2]/2][ic];

		      }

                  }
             } 
	    }
            else if (j >= ((imax[1]-imin[1])/2) - 1 - av[1]/2 && j <= ((imax[1]-imin[1])/2) - 1)
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][0];
                       co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][0];
                       co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][0];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][0][ic];
		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
			  
                       co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1];
                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][avmax[2]-avmin[2]-1][ic];
		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2];

			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][avmax[1]-avmin[1]-1][k-imin[2]-av[2]/2][ic];

		      }
		  }

             }

	    }
	    else
	    {
             for (k = imin[2]; k < ((imax[2]-imin[2])/2); k++)
             {
                  if(k == imin[2] || k == ((imax[2]-imin[2])/2)-1)
                  {

                      co_C[i][j][k] = 0.0;
                      co_CI[i][j][k] = 0.0;
                      co_Prt[i][j][k] = 0.0;
                      if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                      for (ic = 0; ic < MAX_NCOMPS; ic++)
                          co_Sct[i][j][k][ic] = 0.0;
                  }
                  else
                  {
                      if (k >= imin[2] + 1 && k <= imin[2] + av[2]/2)
	              {
                       co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][0];
                       co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][0];
                       co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][0];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][0][ic];

		      }
		      else if (k >= ((imax[2]-imin[2])/2) - 1 - av[2]/2 && k <= ((imax[2]-imin[2])/2) - 2)
	              {
                       co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];
                       co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];
                       co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1];

                       if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                       for (ic = 0; ic < MAX_NCOMPS; ic++)
                           co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][avmax[2]-avmin[2]-1][ic];

		      }
		      else
		      {
			  co_C[i][j][k] = cs_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];
			  co_CI[i][j][k] = ci_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];
			  co_Prt[i][j][k] = prt_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2];
			  if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                          for (ic = 0; ic < MAX_NCOMPS; ic++)
			      co_Sct[i][j][k][ic] = sct_new[avmax[0]-avmin[0]-1][j-imin[1]-av[1]/2][k-imin[2]-av[2]/2][ic];

		      }

                  }
              }
	    }
        }
	}

	//xiaoxue
	int iix,iiy,iiz;
	iix = imax[0]/2;
	iiy = imax[1]/2;
	iiz = imax[2]/2;
	if( 2*iix != imax[0])
	{
	    for (k = imin[2]; k < imax[2]/2; k++)
	    {
		for (j = imin[1]; j < imax[1]/2; j++)
		{
		    C[imax[0]-1][2*j][2*k]  = co_C[iix-1][j][k];
		    C[imax[0]-1][(2*j)+1][2*k]  = co_C[iix-1][j][k];
		    C[imax[0]-1][2*j][(2*k)+1]  = co_C[iix-1][j][k];
		    C[imax[0]-1][(2*j)+1][(2*k)+1]  = co_C[iix-1][j][k];
		    CI[imax[0]-1][2*j][2*k]  = co_CI[iix-1][j][k];
		    CI[imax[0]-1][(2*j)+1][2*k]  = co_CI[iix-1][j][k];
		    CI[imax[0]-1][2*j][(2*k)+1]  = co_CI[iix-1][j][k];
		    CI[imax[0]-1][(2*j)+1][(2*k)+1]  = co_CI[iix-1][j][k];
		    Prt[imax[0]-1][2*j][2*k]  = co_Prt[iix-1][j][k];
		    Prt[imax[0]-1][(2*j)+1][2*k]  = co_Prt[iix-1][j][k];
		    Prt[imax[0]-1][2*j][(2*k)+1]  = co_Prt[iix-1][j][k];
		    Prt[imax[0]-1][(2*j)+1][(2*k)+1]  = co_Prt[iix-1][j][k];
		    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
		    {
			for (ic = 0; ic < MAX_NCOMPS; ic++)
			{
			    Sct[imax[0]-1][2*j][2*k][ic]  = co_Sct[iix-1][j][k][ic];
			    Sct[imax[0]-1][(2*j)+1][2*k][ic]  = co_Sct[iix-1][j][k][ic];
			    Sct[imax[0]-1][2*j][(2*k)+1][ic]  = co_Sct[iix-1][j][k][ic];
			    Sct[imax[0]-1][(2*j)+1][(2*k)+1][ic]  = co_Sct[iix-1][j][k][ic];
			}
		    }
		}
	    }
	}
	//xiaoxue

        for (k = imin[2]; k < imax[2]/2; k++)
        {
        for (j = imin[1]; j < imax[1]/2; j++)
        {
            for (i = imin[0]; i < imax[0]/2; i++)
            {
                C[2*i][2*j][2*k]  = co_C[i][j][k];
                C[(2*i)+1][2*j][2*k]  = co_C[i][j][k];
                C[2*i][(2*j)+1][2*k]  = co_C[i][j][k];
                C[(2*i)+1][(2*j)+1][2*k]  = co_C[i][j][k];
                C[2*i][2*j][(2*k)+1]  = co_C[i][j][k];
                C[(2*i)+1][2*j][(2*k)+1]  = co_C[i][j][k];
                C[2*i][(2*j)+1][(2*k)+1]  = co_C[i][j][k];
                C[(2*i)+1][(2*j)+1][(2*k)+1]  = co_C[i][j][k];

                CI[2*i][2*j][2*k]  = co_CI[i][j][k];
                CI[(2*i)+1][2*j][2*k]  = co_CI[i][j][k];
                CI[2*i][(2*j)+1][2*k]  = co_CI[i][j][k];
                CI[(2*i)+1][(2*j)+1][2*k]  = co_CI[i][j][k];
                CI[2*i][2*j][(2*k)+1]  = co_CI[i][j][k];
                CI[(2*i)+1][2*j][(2*k)+1]  = co_CI[i][j][k];
                CI[2*i][(2*j)+1][(2*k)+1]  = co_CI[i][j][k];
                CI[(2*i)+1][(2*j)+1][(2*k)+1]  = co_CI[i][j][k];


                Prt[2*i][2*j][2*k]  = co_Prt[i][j][k];
                Prt[(2*i)+1][2*j][2*k]  = co_Prt[i][j][k];
                Prt[2*i][(2*j)+1][2*k]  = co_Prt[i][j][k];
                Prt[(2*i)+1][(2*j)+1][2*k]  = co_Prt[i][j][k];
                Prt[2*i][2*j][(2*k)+1]  = co_Prt[i][j][k];
                Prt[(2*i)+1][2*j][(2*k)+1]  = co_Prt[i][j][k];
                Prt[2*i][(2*j)+1][(2*k)+1]  = co_Prt[i][j][k];
                Prt[(2*i)+1][(2*j)+1][(2*k)+1]  = co_Prt[i][j][k];

                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                for (ic = 0; ic < MAX_NCOMPS; ic++)
                {
                    Sct[2*i][2*j][2*k][ic]  = co_Sct[i][j][k][ic];
                    Sct[(2*i)+1][2*j][2*k][ic]  = co_Sct[i][j][k][ic];
                    Sct[2*i][(2*j)+1][2*k][ic]  = co_Sct[i][j][k][ic];
                    Sct[(2*i)+1][(2*j)+1][2*k][ic]  = co_Sct[i][j][k][ic];
                    Sct[2*i][2*j][(2*k)+1][ic]  = co_Sct[i][j][k][ic];
                    Sct[(2*i)+1][2*j][(2*k)+1][ic]  = co_Sct[i][j][k][ic];
                    Sct[2*i][(2*j)+1][(2*k)+1][ic]  = co_Sct[i][j][k][ic];
                    Sct[(2*i)+1][(2*j)+1][(2*k)+1][ic]  = co_Sct[i][j][k][ic];
                }
                }
            }
        }
        }

        for (k = imin[2]; k < imax[2]; k++)            
        {
            icoords[2] = k;
            for (j = imin[1]; j < imax[1]; j++)
            {
                icoords[1] = j;
                for (i = imin[0]; i < imax[0]; i++)
                {
                    icoords[0] = i;
                    state = Rect_state(icoords,wv);
                    if(g_compute_NS_terms(Rect_state(icoords,wv)))
                    {

                            fill_3d_27pt_Pstencil(i,j,k,nsten);

                            ux = (vel(0,nsten->sts3d[2][1][1]) -
                                vel(0,nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                            uy = (vel(0,nsten->sts3d[1][2][1]) -
                                vel(0,nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                            uz = (vel(0,nsten->sts3d[1][1][2]) -
                                vel(0,nsten->sts3d[1][1][0])) / (2.0*dh[2]);

                            vx = (vel(1,nsten->sts3d[2][1][1]) -
                                vel(1,nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                            vy = (vel(1,nsten->sts3d[1][2][1]) -
                                vel(1,nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                            vz = (vel(1,nsten->sts3d[1][1][2]) -
                                vel(1,nsten->sts3d[1][1][0])) / (2.0*dh[2]);

                            wx = (vel(2,nsten->sts3d[2][1][1]) -
                                vel(2,nsten->sts3d[0][1][1])) / (2.0*dh[0]);
                            wy = (vel(2,nsten->sts3d[1][2][1]) -
                                vel(2,nsten->sts3d[1][0][1])) / (2.0*dh[1]);
                            wz = (vel(2,nsten->sts3d[1][1][2]) -
                                vel(2,nsten->sts3d[1][1][0])) / (2.0*dh[2]);


                            S11 = ux;
                            S12 = (uy + vx)/2.0;
                            S13 = (uz + wx)/2.0;
                            S21 = (vx + uy)/2.0;
                            S22 = vy;
                            S23 = (vz + wy)/2.0;
                            S31 = (wx + uz)/2.0;
                            S32 = (wy + vz)/2.0;
                            S33 = wz;

                            dtemx = (temperature(nsten->sts3d[2][1][1])
                                          -temperature(nsten->sts3d[0][1][1]))/(2.0*dh[0]);
                            dtemy = (temperature(nsten->sts3d[1][2][1])
                                          -temperature(nsten->sts3d[1][0][1]))/(2.0*dh[1]);
                            dtemz = (temperature(nsten->sts3d[1][1][2])
                                          -temperature(nsten->sts3d[1][1][0]))/(2.0*dh[2]);
  
                            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                            {
                            for (ic = 0; ic < MAX_NCOMPS; ic++)
                            {
                                pdenx[ic] = (pdens(nsten->sts3d[2][1][1])[ic]/Dens(nsten->sts3d[2][1][1])
                                          - pdens(nsten->sts3d[0][1][1])[ic]/Dens(nsten->sts3d[0][1][1]))/(2.0*dh[0]);
                                pdeny[ic] = (pdens(nsten->sts3d[1][2][1])[ic]/Dens(nsten->sts3d[1][2][1])
                                          -pdens(nsten->sts3d[1][0][1])[ic]/Dens(nsten->sts3d[1][0][1]))/(2.0*dh[1]);
                                pdenz[ic] = (pdens(nsten->sts3d[1][1][2])[ic]/Dens(nsten->sts3d[1][1][2])
                                          -pdens(nsten->sts3d[1][1][0])[ic]/Dens(nsten->sts3d[1][1][0]))/(2.0*dh[2]);
                            }
                            }

                            S = sqrt(2.0*( (S11*S11) + (S12*S12)
                                    + (S13*S13) + (S21*S21)
                                    + (S22*S22) + (S23*S23)
                                    + (S31*S31) + (S32*S32)
                                    + (S33*S33)));

                            S11_a = S11 - ((S11+S22+S33)/3.0);
                            S12_a = S12;
                            S13_a = S13;
                            S21_a = S21;
                            S22_a = S22 - ((S11+S22+S33)/3.0);
                            S23_a = S23;
                            S31_a = S31;
                            S32_a = S32;
                            S33_a = S33 - ((S11+S22+S33)/3.0);
/*
                            if(fabs(Prt[i][j][k]) > 0.5)
                            Prt[i][j][k] = 0.0018;
                            if(fabs(C[i][j][k]) > 0.5)
                            C[i][j][k] = 0.033;
                            if(fabs(CI[i][j][k]) > 1.0)
                            CI[i][j][k] = 0.1;
                            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                            for (ic = 0; ic < MAX_NCOMPS; ic++)
                                if(fabs(Sct[i][j][k][ic]) > 0.5)
                                    Sct[i][j][k][ic] = 0.0014;
*/
			    //add by xiaoxue
/*
			    Para(state)[0] =  C[i][j][k];
			    Para(state)[1] =  CI[i][j][k];
			    Para(state)[2] =  Prt[i][j][k];
*/
			    //Para(state)[0] =  C[i][j][k]*delta*delta*S;
			    //Para(state)[1] =  Dens(nsten->sts3d[1][1][1])*C_P(nsten->sts3d[1][1][1])*Prt[i][j][k]*delta*delta*S;
			    //Para(state)[2] =  Prt[i][j][k]*delta*delta*S;

			    //for( ic = 0; ic < 10; ic++)
			    //{
			    //	Para(state)[3 + ic] = Sct[i][j][k][3 + ic]*delta*delta*S;
			    //}
			    //int ii;
			    //for( ii = 0; ii < MAX_NCOMPS; ii++)
				//Sct(state)[ii] = Sct[i][j][k][ii]*delta*delta*S;

                            Tau(state)[0][0] = ((2/3)*CI[i][j][k]*Dens(nsten->sts3d[1][1][1])*delta*delta*S*S)
                                      - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S11_a);
                            Tau(state)[0][1] = - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S12_a);
                            Tau(state)[0][2] = - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S13_a);
                            Tau(state)[1][0] = - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S21_a);

                            Tau(state)[1][1] = ((2/3)*CI[i][j][k]*Dens(nsten->sts3d[1][1][1])*delta*delta*S*S)
                                      - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S22_a);
                            Tau(state)[1][2] = - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S23_a);
                            Tau(state)[2][0] = - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S31_a);
                            Tau(state)[2][1] = - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S32_a);
                            Tau(state)[2][2] = ((2/3)*CI[i][j][k]*Dens(nsten->sts3d[1][1][1])*delta*delta*S*S)
                                      - 2*Dens(nsten->sts3d[1][1][1])*C[i][j][k]*delta*delta*S*(S33_a);

                            Qh(state)[0] = - (Dens(nsten->sts3d[1][1][1])*C_P(nsten->sts3d[1][1][1])*Prt[i][j][k]*delta*delta*S)*dtemx;
                            Qh(state)[1] = - (Dens(nsten->sts3d[1][1][1])*C_P(nsten->sts3d[1][1][1])*Prt[i][j][k]*delta*delta*S)*dtemy;
                            Qh(state)[2] = - (Dens(nsten->sts3d[1][1][1])*C_P(nsten->sts3d[1][1][1])*Prt[i][j][k]*delta*delta*S)*dtemz;

                            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                            {
                            /*
                            for (ic = 0; ic < MAX_NCOMPS; ic++)
                            {
                                Qc(state)[0][ic] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][ic]*delta*delta*S)*pdenx[ic];
                                Qc(state)[1][ic] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][ic]*delta*delta*S)*pdeny[ic];
                                Qc(state)[2][ic] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][ic]*delta*delta*S)*pdenz[ic];
                            } 
                            */

                            Qc(state)[0] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][1]*delta*delta*S)*pdenx[1];
                            Qc(state)[1] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][1]*delta*delta*S)*pdeny[1];
                            Qc(state)[2] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][1]*delta*delta*S)*pdenz[1];
                            Qc(state)[3] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][0]*delta*delta*S)*pdenx[0];
                            Qc(state)[4] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][0]*delta*delta*S)*pdeny[0];
                            Qc(state)[5] = - (Dens(nsten->sts3d[1][1][1])*Sct[i][j][k][0]*delta*delta*S)*pdenz[0];
                           
                            }
                      }

                 }
            }
        }
	
	    free(co_NS_terms);    
            free(co_coords_x);
            free(co_coords_y);
	    free(co_coords_z);
            free(co_rho);
            free(co_rho_cp);
            free(co_cp);
            free(co_rho_t);
            free(co_rho_c);
            free(co_t);
            free(co_con);
            free(co_rho_cu);
            free(co_rho_cv);
	    free(co_rho_cw);
            free(co_rho_cp_tu);
            free(co_rho_cp_tv);
	    free(co_rho_cp_tw);
            free(co_vel_u);
            free(co_vel_v);
	    free(co_vel_w);
            free(co_rho_vel_u);
            free(co_rho_vel_v);
	    free(co_rho_vel_w); 
            free(co_rho_vel_uu);
            free(co_rho_vel_vv);
	    free(co_rho_vel_ww);
	    free(co_rho_vel_uv);
	    free(co_rho_vel_vw);
	    free(co_rho_vel_uw);
            free(co_temp11);
            free(co_temp12);
	    free(co_temp13);
            free(co_temp21);
            free(co_temp22);
	    free(co_temp23);
            free(co_temp31);
            free(co_temp32);
            free(co_temp33);
            free(co_temp);
            free(co_temp_t1);
            free(co_temp_t2);
            free(co_temp_t3);
            free(co_temp_c1);
            free(co_temp_c2);
            free(co_temp_c3);
            free(MA11);
            free(MA12);
            free(MA13);
            free(MA22);
            free(MA23);
            free(LA11);
            free(LA12);
            free(LA13);
            free(LA22);
            free(LA23);
            free(LA33);
            free(LA1);
            free(LA2);
            free(LA3);
            free(LA4);
            free(LA5);
            free(LI);
            free(MI);
            free(MH1);
            free(MH2);
            free(MH3);
            free(LH1);
            free(LH2);
            free(LH3);
            free(MC1);
            free(MC2);
            free(MC3);
            free(LC1);
            free(LC2);
            free(LC3);
            free(MMA);
            free(LMA);
            free(LMH);
            free(MMH);
            free(LMC);
            free(MMC);
            free(r);
/*          free(r_color);
            free(r_color_LA);
            free(r_color_LI);
            free(r_color_LH);
            free(r_color_LC);
            free(r_color_LC0);     */
	    free(cs_new);
            free(ci_new);
            free(prt_new);
            free(sct_new);
            free(deno);
            free(nume);
            free(ci_deno);
            free(ci_nume);
            free(prt_deno);
            free(prt_nume);
            free(sct_deno);
            free(sct_nume);
            free(re_NS_terms);
            free(re_coords_x);
            free(re_coords_y);
            free(re_coords_z);
            free(re_rho);
            free(re_rho_cp);
            free(re_cp);
            free(re_rho_t);
            free(re_rho_c);
            free(re_t);
            free(re_con);
            free(re_rho_cu);
            free(re_rho_cv);
            free(re_rho_cw);
            free(re_rho_cp_tu);
            free(re_rho_cp_tv);
            free(re_rho_cp_tw);
            free(re_vel_u);
            free(re_vel_v);
            free(re_vel_w);
            free(re_rho_vel_u);
            free(re_rho_vel_v);
            free(re_rho_vel_w);
            free(re_rho_vel_uu);
            free(re_rho_vel_vv);
            free(re_rho_vel_ww);
            free(re_rho_vel_uv);
            free(re_rho_vel_vw);
            free(re_rho_vel_uw);
            free(re_temp11);
            free(re_temp12);
            free(re_temp13);
            free(re_temp21);
            free(re_temp22);
            free(re_temp23);
            free(re_temp31);
            free(re_temp32);
            free(re_temp33);
            free(re_temp);
            free(re_temp_t1);
            free(re_temp_t2);
            free(re_temp_t3);
            free(re_temp_c1);
            free(re_temp_c2);
            free(re_temp_c3);
            free(C);
            free(co_C);
            free(CI);
            free(co_CI);
            free(Prt);
            free(co_Prt);
            free(Sct);
            free(co_Sct);
}

/*  The following uses a semi-Lagrangian scheme to compute the 
    Lagrangian derivative of the entropy.  If not equal to zero, 
    a shock is present, and the SGS terms should be turned off.*/

LOCAL void mesh_avr (Wave *wv,float ***var, float ***av_var, int *min, int *max, int *av)
{
RECT_GRID       *rgr = wv->rect_grid;
 /* Define averaged variable */
 float sum = 0;
 int avmin[3] = {0,0,0}, avmax[3];
 int i, j, k, ii, jj, kk;

 /* Compute number of grid blocks for the new (averaged) variable */
 for (i = 0; i < 3; i++)
     avmax[i] = max[i] - min[i] - av[i] + 1;

 /* For every average block */
 for (i = avmin[0]; i < (avmax[0] - avmin[0]); i++)
 {
     for (j = avmin[1]; j < (avmax[1] - avmin[1]); j++)
     {
          for (k = avmin[2]; k < (avmax[2] - avmin[2]); k++)
           /* Inside every average block, compute the sum, 
	      and divide it by the number of meshblocks */
          {
              //if (cell_center(k,2,rgr) < 3.0)
              {
              for (ii = i + min[0]; ii < i + min[0] + av[0]; ii++)
              {
                  for (jj = j + min[1]; jj < j + min[1] + av[1]; jj++)
                  {
                      for (kk = k + min[2]; kk < k + min[2] + av[2]; kk++)
                      {
              		   if (var[ii][jj][kk] >= 0)
                           	sum += var[ii][jj][kk];
                      }
                  }
              }
              av_var[i][j][k] = sum / (av[0] * av[1] * av[2]);
              sum = 0;      
              }
              //else
              //    av_var[i][j][k] = 0.0;
          }
     }
 }
}

LOCAL void mesh_avr_mc (Wave *wv,float ****var, float ****av_var, int *min, int *max, int *av)
{
RECT_GRID       *rgr = wv->rect_grid;
 /* Define averaged variable */
 float sum = 0;
 int avmin[3] = {0,0,0}, avmax[3];
 int i, j, k, ii, jj, kk, ic;

 /* Compute number of grid blocks for the new (averaged) variable */
 for (i = 0; i < 3; i++)
     avmax[i] = max[i] - min[i] - av[i] + 1;

 /* For every average block */
for (ic = 0; ic < MAX_NCOMPS; ic++)
{
 for (i = avmin[0]; i < (avmax[0] - avmin[0]); i++)
 {
     for (j = avmin[1]; j < (avmax[1] - avmin[1]); j++)
     {
          for (k = avmin[2]; k < (avmax[2] - avmin[2]); k++)
           /* Inside every average block, compute the sum, 
              and divide it by the number of meshblocks */
          {
              //if (cell_center(k,2,rgr) < 3.0)
              {
              for (ii = i + min[0]; ii < i + min[0] + av[0]; ii++)
              {
                  for (jj = j + min[1]; jj < j + min[1] + av[1]; jj++)
                  {
                      for (kk = k + min[2]; kk < k + min[2] + av[2]; kk++)
                      {
              		   if (var[ii][jj][kk][ic] >= 0)
                           	sum += var[ii][jj][kk][ic];
                      }
                  }
              }
              av_var[i][j][k][ic] = sum / (av[0] * av[1] * av[2]);
              sum = 0.0;
              }
              //else
              //    av_var[i][j][k] = 0.0;
          }
     }
 }
}
}


#endif /* defined SUBGRID */

LOCAL   void parab_front2d(
        float    dt,
        Wave     *wave,
        Wave     *newwave,
        Front    *front,
        Front    **newfront)
{
        CURVE      **c;
        CURVE      *oldc;
        CURVE      *newc;
        int        status;
        INTERFACE  *tempintfc;
        tempintfc = NULL;

        *newfront = copy_front(front);
        set_size_of_intfc_state(size_of_state(front->interf));
        set_copy_intfc_states(YES);
        (*newfront)->interf = pp_copy_interface(front->interf);

        for (c = front->interf->curves; c && *c; ++c)
        {
            oldc = *c;
            if (((newc = correspond_curve(oldc)) != NULL) &&
                 (correspond_curve(newc) != NULL))
            {
                parab_front2d_update(front,wave,newwave,oldc,newc,dt);
            }
        }

        if (!scatter_front(*newfront))
        {
            (void) printf("ERROR scatter_front in parab_driver\n");
           clean_up(ERROR);
        }
}

LOCAL   void parab_front3d(
        float    dt,
        Wave     *wave,
        Wave     *newwave,
        Front    *front,
        Front    **newfront)
{
        *newfront = copy_front(front);
        set_size_of_intfc_state(size_of_state(front->interf));
        set_copy_intfc_states(YES);
        (*newfront)->interf = pp_copy_interface(front->interf);
        parab_front3d_update(dt,wave,newwave,front,*newfront);
        if (!scatter_front(*newfront))
        {
            (void) printf("ERROR scatter_front in parab_driver\n");
           clean_up(ERROR);
        }
}

LOCAL   void parab_front2d_update(
        Front    *fr,
        Wave     *wave,
        Wave     *newwave,
        CURVE    *oldc,
        CURVE    *newc,
        float    dt)
{
        BOND            *oldb = oldc->first;
        BOND            *newb = newc->first;
        float           V[MAXD];
        int             dim = fr->interf->dim;
        float           L[MAXD],U[MAXD];        /* propagation boundary */

        debug_print("parab_front2d_update","Entered parab_front2d_update\n");

        if ((fr->_point_propagate == NULL)   ||
            (oldc == NULL)                   ||
            (newc == NULL)                   ||
            (correspond_curve(oldc) != newc) ||
            (correspond_curve(newc) != oldc))
            return;

        set_propagation_bounds(fr,L,U);
        while (oldb)
        {
            if ((oldb != oldc->last) && (!n_pt_propagated(newb->end)))
            {
                n_pt_propagated(newb->end) = YES;
                if (out_of_bound(oldb->end,L,U,dim))
                {
                    Locstate newsl,newsr;
                    Locstate oldsl,oldsr;
                    slsr(newb->end,Hyper_surf_element(newb),Hyper_surf(newc),
                                &newsl,&newsr);
                    slsr(oldb->end,Hyper_surf_element(oldb),Hyper_surf(oldc),
                                &oldsl,&oldsr);
                    ft_assign(newsl,oldsl,fr->sizest);
                    ft_assign(newsr,oldsr,fr->sizest);
                    continue;
                }
                parab_front_point(fr,wave,oldb->end,newb->end,Hyper_surf_element(oldb->next),Hyper_surf(oldc),dt);
            }
            if (oldb == oldc->last)
                break;
            oldb = oldb->next;
            newb = newb->next;
        }
        debug_print("parab_front2d_update","Leaving parab_front2d_update\n");
}

LOCAL   void parab_front3d_update(
        float    dt,
        Wave     *wave,
        Wave     *newwave,
        Front    *front,
        Front    *newfront)
{
        INTERFACE               *intfc_old = front->interf;
        INTERFACE               *intfc_new = newfront->interf;
        HYPER_SURF              *oldhs, *newhs;
        HYPER_SURF_ELEMENT      *oldhse, *newhse;
        POINT                   *oldp, *newp;
        DEBUG_ENTER(parab_propagate_surface_points)

        start_clock("parab_surface_propagate");

        (void) next_point(intfc_old,NULL,NULL,NULL);
        (void) next_point(intfc_new,NULL,NULL,NULL);
        while (next_point(intfc_old,&oldp,&oldhse,&oldhs) &&
             next_point(intfc_new,&newp,&newhse,&newhs))
        {
            if(Boundary_point(newp))
                continue;
            parab_front_point(front,wave,oldp,newp,oldhse,oldhs,dt);
        }
        stop_clock("parab_surface_propagate");
        DEBUG_LEAVE(parab_propagate_surface_points)
}               //end propagate_surface_points

LOCAL   void parab_front_point(
        Front              *fr,
        Wave               *wave,
        POINT              *oldp,
        POINT              *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        float              dt)
{
        Locstate sl, sr;
       
        slsr(oldp,oldhse,oldhs,&sl,&sr);
        ft_assign(left_state(newp),sl,fr->sizest);
        ft_assign(right_state(newp),sr,fr->sizest);
        parab_point_update(fr,left_state(newp));
        parab_point_update(fr,right_state(newp));
}

LOCAL   void parab_point_update(
        Front           *front,
        Locstate        ans)
{
        INTERFACE       *intfc = front->interf;
        Gas_param       **params;
        bool            viscosity, mass_diffusion, thermal_conduction;
        bool            subgrid_vis, subgrid_md, subgrid_con;
        float           mu, d, kappa;
        float           dt = front->dt;
        int             i,dim;

        dim = front->rect_grid->dim;
        params =    gas_params_list(intfc);
        mu = params[0]->avisc.viscosity_coef;
        d = params[0]->avisc.diffusivity_coef;
        kappa = params[0]->avisc.conductivity_coef;
        subgrid_vis = params[0]->avisc.subgrid_vis;
        subgrid_md = params[0]->avisc.subgrid_md;
        subgrid_con = params[0]->avisc.subgrid_con;
        subgrid_vis = NO;
        subgrid_md = NO;
        if( fabs(mu) < 1e-10)
            viscosity = NO;
        else
            viscosity = YES;

        if( fabs(d) < 1e-10)
            mass_diffusion = NO;
        else
            mass_diffusion = YES;

        if(fabs(kappa) <= 1e-16)
            thermal_conduction = NO;
        else
            thermal_conduction = YES;

        if (is_obstacle_state(ans))
          return;

        switch (dim)
        {
#if defined ONED
        case 1:
        {
            //need to update
        }
        break;

#endif /* defined ONED */
#if defined TWOD
        case 2:
        {
        if (viscosity == YES)
        {
            Mom(ans)[0] += dt*mu*Mvisx(ans);
            Mom(ans)[1] += dt*mu*Mvisy(ans);
            Energy(ans) += dt*mu*Evis(ans);
        }
        if (thermal_conduction == YES)
            Energy(ans) += dt*kappa*Ethermal(ans);
        if (mass_diffusion == YES)
        {
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                if(Params(ans)->n_comps != 1)
                {
                    for(i = 0; i < Params(ans)->n_comps; i++)
                    {
                        pdens(ans)[i] += dt*d*Cmd(ans)[i];
#if defined(SUBGRID)
                        //if(subgrid_md == YES)
                        //    pdens(ans)[i] -= dt*Scmd(ans)[i];
#endif /* defined SUBGRID */
                    }
                }
                float sum = 0.0;
                for(i = 0; i < Params(ans)->n_comps; i++)
                sum += pdens(ans)[i];
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] *= Dens(ans)/sum;

                Energy(ans) += dt*d*Emd(ans);
            }
        }
#if defined(SUBGRID)
        /*
        if(subgrid_con == YES)
        {
             Energy(ans) -= dt*Sethermal(ans);
        }
        if(subgrid_vis == YES)
        {
             Mom(ans)[0] -= dt*Svisx(ans);
             Mom(ans)[1] -= dt*Svisy(ans);
        }
        */
#endif /* defined SUBGRID */
        }
        break;

#endif /* defined TWOD */
#if defined THREED
        case 3:
        {
        if (viscosity == YES)
        {
            Mom(ans)[0] += dt*mu*Mvisx(ans);
            Mom(ans)[1] += dt*mu*Mvisy(ans);
            Mom(ans)[2] += dt*mu*Mvisz(ans);
            Energy(ans) += dt*mu*Evis(ans);
        }
        if (thermal_conduction == YES)
        {
            Energy(ans) += dt*kappa*Ethermal(ans);
        }
        if (mass_diffusion == YES)
        {
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                if(Params(ans)->n_comps != 1)
                {
                    for(i = 0; i < Params(ans)->n_comps; i++)
                    {
                        pdens(ans)[i] += dt*d*Cmd(ans)[i];
#if defined(SUBGRID)
                        //if(subgrid_md == YES)
                        //    pdens(ans)[i] -= dt*Scmd(ans)[i];
#endif /* defined SUBGRID */
                    }
                }
                float sum = 0.0;
                for(i = 0; i < Params(ans)->n_comps; i++)
                sum += pdens(ans)[i];
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] *= Dens(ans)/sum;

                Energy(ans) += dt*d*Emd(ans);
            }
        }
#if defined(SUBGRID)
        //if(subgrid_con == YES)
        //{
        //     Energy(ans) -= dt*Sethermal(ans);
        //}
        //if(subgrid_vis == YES)
        //{
        //     Mom(ans)[0] -= dt*Svisx(ans);
        //     Mom(ans)[1] -= dt*Svisy(ans);
        //     Mom(ans)[2] -= dt*Svisz(ans);
        //}
#endif /* defined SUBGRID */
        }
        break;

#endif /* defined THREED */
        }
}

EXPORT  void parab_tan_solver2d(
        float           ds,
        float           dt,
        Tan_stencil     *sten,
        Locstate        ansl,
        Locstate        ansr,
        Front           *fr)
{
        if (sten->hs[0] == NULL)
            return;

        if(wave_type(sten->hs[0]) == CONTACT)
        {
            sten->states = sten->leftst;
            one_side_parab_tan_solver2d(ds,dt,sten,ansl,fr);
            sten->states = sten->rightst;
            one_side_parab_tan_solver2d(ds,dt,sten,ansr,fr);
        }
}

EXPORT  void  parab_tan_solver3d(
        Front         *fr,
        const Tparams *tp,
        const Tparams *tpp,
        POINT         *newp)
{
        Locstate            ansl, ansr;
        int                 il, ir;
        Locstate            sl,sr;
        static  Tan_stencil *sten = NULL;
        float               ds, dt;
 
        if (sten == NULL)
            sten = alloc_tan_stencil(fr,fr->npts_tan_sten/2);

        ir =  sten->npts/2 + 1;
        il = -sten->npts/2 - 1;

        if (tp->tnl.tri != NULL)
            sten->hse[0] = Hyper_surf_element(tp->tnl.tri);
        else if (tp->tnr.tri != NULL)
            sten->hse[0] = Hyper_surf_element(tp->tnr.tri);
        else
        {
            (void) printf("WARNING in parab_tan_solver(), "
                          "can't set sten->hse[0]\n");
            return;
        }
        sten->hs[0] = tp->hs;
        sten->p[0] = tp->p;

        slsr(tp->p,Hyper_surf_element(tp->tnl.tri),tp->hs,&sl,&sr);
        ft_assign(sten->leftst[0],sl,fr->sizest);
        ft_assign(sten->rightst[0],sr,fr->sizest);

        fill_stencil_in_direction(fr,sten,-1,il,tp);
        fill_stencil_in_direction(fr,sten,1,ir,tp);

        ansl = left_state(newp);
        ansr = right_state(newp);
        sten->newhs = tp->hs;
        sten->dir = tp->tan;
        sten->dir1 = tpp->tan;
        sten->nor = tp->nor;
        ds = tp->ds;
        dt = tp->dt;  

        if (tp->hs == NULL)
            return;

        if(wave_type(tp->hs) == CONTACT)
        {
            sten->states = sten->leftst;
            one_side_parab_tan_solver3d(ds,dt,sten,ansl,fr);
            sten->states = sten->rightst;
            one_side_parab_tan_solver3d(ds,dt,sten,ansr,fr);
        }
}

LOCAL void one_side_parab_tan_solver2d(
        float           ds,
        float           dt,
        Tan_stencil     *sten,
        Locstate        ans,
        Front           *fr)
{
        RECT_GRID      *gr = fr->rect_grid;
        int            i, j, dim = gr->dim;
        Locstate       *state;
        Locstate       *states = sten->states;
        Locstate         sll;
        Locstate          sr;
        Locstate         srr;

        if (is_obstacle_state(ans))
          return;

        sll  =  states[-1];
        sr   =  states[0];
        srr  =  states[1];

        float  nor[dim];
        float  norsll, norsr, norsrr;
        float  tan1sll, tan1sr, tan1srr;
        float  tan1[dim];
        norsll = norsr =  norsrr = 0.0;
        tan1sll = tan1sr =  tan1srr = 0.0;
        normal(sten->p[0],sten->hse[0],sten->hs[0],nor,fr);
        for ( i = 0; i < dim; ++i)
        {
           tan1[i] = sten->dir[i];
        }
        for ( i = 0; i < dim; ++i)
        {
            norsll += nor[i] * vel(i,sll);
            norsr += nor[i] * vel(i,sr);
            norsrr += nor[i] * vel(i,srr);
        }

        for ( i = 0; i < dim; ++i)
        {
            tan1sll += tan1[i] * vel(i,sll);
            tan1sr += tan1[i] * vel(i,sr);
            tan1srr += tan1[i] * vel(i,srr);
        }
	
        int    N = Params(sr)->n_comps;
        float  enthsll, enthsl, enthsr, enthsrr;
        float  pdenssll[N], pdenssl[N], pdenssr[N], pdenssrr[N];
        float  dvnt1t1, dvt1t1, dvnt1, dvt1t1t1;
        float  dcst1, cst1, dcit1;        
        float  dtt1t1, dtst1t1, dmdt1t1[N], dmdst1t1[N], energymdt1t1;

        dvnt1t1 = (norsrr - (2*norsr) + norsll)/(ds*ds);
        dvt1t1 = (tan1srr - tan1sll)/(2*ds);
        dvnt1 = (norsrr - norsll)/(2*ds);
        dvt1t1t1 = (tan1srr - (2*tan1sr) + tan1sll)/(ds*ds);
        dtt1t1 = (temperature(srr) - (2*temperature(sr)) + temperature(sll))/(ds*ds);
#if defined(SUBGRID)
/*
        C(ans) = C(sr);
        CI(ans) = CI(sr);
        PRT(ans) = PRT(sr);
        SCT(ans) = SCT(sr);
        SCT0(ans) = SCT0(sr);
        STR(ans) = STR(sr);
        cst1 = -2*C(sr)*ds*ds*Dens(sr)*STR(sr);
        dcst1 = ( (-2*C(srr)*ds*ds*Dens(srr)*STR(srr))
                - (-2*C(sll)*ds*ds*Dens(sll)*STR(sll)) )/(2*ds);
        dcit1 = ( ((2.0/2.0)*CI(srr)*ds*ds*Dens(srr)*STR(srr)*STR(srr))
                - ((2.0/2.0)*CI(sll)*ds*ds*Dens(sll)*STR(sll)*STR(sll)) )/(2*ds);
        dtst1t1 = ((((Dens(srr)*C_P(srr)*PRT(srr)*ds*ds*STR(srr))
                  - (Dens(sll)*C_P(sll)*PRT(sll)*ds*ds*STR(sll)))/(2*ds))
                  *((temperature(srr) - temperature(sll))/(2*ds)))
                  + ((Dens(sr)*C_P(sr)*PRT(sr)*ds*ds*STR(sr))*(dtt1t1));
*/
#endif /* defined SUBGRID */

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            enthsll = specific_internal_energy_h_species(sll) 
                      - specific_internal_energy_l_species(sll);
            enthsr = specific_internal_energy_h_species(sr)
                     - specific_internal_energy_l_species(sr);
            enthsrr = specific_internal_energy_h_species(srr) 
                      - specific_internal_energy_l_species(srr);

            if(Params(sr)->n_comps != 1)
            {
                for(i = 0; i < Params(sr)->n_comps; i++)
                {
                    pdenssrr[i] = pdens(srr)[i]/Dens(srr);
                    pdenssr[i] = pdens(sr)[i]/Dens(sr);
                    pdenssll[i] = pdens(sll)[i]/Dens(sll);

                    dmdt1t1[i] = (Dens(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(ds*ds))
                                  + (((Dens(srr)) - (Dens(sll)))/(2*ds))*((pdenssrr[i] - pdenssll[i])/(2*ds));
 
#if defined(SUBGRID)
/*
                    dmdst1t1[i] = ((((Dens(srr)*SCT(srr)*ds*ds*STR(srr))
                                   - (Dens(sll)*SCT(sll)*ds*ds*STR(sll)))/(2*ds))
                                   *((pdenssrr[i] - pdenssll[i])/(2*ds)))
                                   + ((Dens(sr)*SCT(sr)*ds*ds*STR(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(ds*ds)));
*/
#endif /* defined SUBGRID */
                }
            }

            energymdt1t1 = (enthsr*Dens(sr))*((pdenssrr[1] - (2*pdenssr[1]) + pdenssll[1])/(ds*ds))
                            + (((enthsrr*Dens(srr)) - (enthsll*Dens(sll)))/(2*ds))*((pdenssrr[1] - pdenssll[1])/(2*ds));
        }

        Mvisx(ans) = Dvnnn(ans) + dvnt1t1;
        Mvisy(ans) = Dvt1nn(ans) + dvt1t1t1;
        Evis(ans) = ( (( Dvnnn(ans) + dvnt1t1 )*norsr)
                    + (( Dvt1nn(ans) + dvt1t1t1 )*tan1sr))
                    + (( Dvnn(ans) - dvt1t1 )*(Dvnn(ans)))
                    + (( Dvt1n(ans) + dvnt1 )*(dvnt1))
                    + (( dvnt1 + Dvt1n(ans) )*(Dvt1n(ans)))
                    + (( dvt1t1 - Dvnn(ans) )*(dvt1t1));
        Ethermal(ans) = Ethermal(ans) + dtt1t1;

#if defined(SUBGRID)
/*
        Svisx(ans) = (DCSn(ans)*(((1.0/2.0)*Dvnn(ans)) - ((1.0/2.0)*dvt1t1)))
                     + ((CSn(ans)*(((1.0/2.0)*Dvnnn(ans))))) + DCIn(ans) 
                     + (dcst1*(((1.0/2.0)*Dvt1n(ans)) + ((1.0/2.0)*dvnt1)))
                     + (cst1*(((1.0/2.0)*dvnt1t1)))
                     + dcit1 ;
        Svisy(ans) = (DCSn(ans)*(((1.0/2.0)*dvnt1) + ((1.0/2.0)*Dvt1n(ans))))
                     + (CSn(ans)*((1.0/2.0)*Dvt1nn(ans))) + DCIn(ans)
                     + (dcst1*(((1.0/2.0)*dvt1t1) - ((1.0/2.0)*Dvnn(ans))))
                     + (cst1*((1.0/2.0)*Dvnnn(ans)))
                     + dcit1 ;
        Sethermal(ans) = Sethermal(ans) + dtst1t1;
*/
#endif /* defined SUBGRID */

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(sr)->n_comps != 1)
            {
                for(i = 0; i < Params(sr)->n_comps; i++)
                {
                    Cmd(ans)[i] =  Cmd(ans)[i] +dmdt1t1[i];
#if defined(SUBGRID)
                    //Scmd(ans)[i] = Scmd(ans)[i] + dmdst1t1[i];
#endif /* defined SUBGRID */
                }
             }
             Emd(ans) = Emd(ans) + energymdt1t1;
        }

        float mvisn, mvist1;
        float svisn, svist1;
        mvisn = Mvisx(ans);
        mvist1 = Mvisy(ans);
        Mvisx(ans) = (mvisn*nor[0]) + (mvist1*tan1[0]);
        Mvisy(ans) = (mvisn*nor[1]) + (mvist1*tan1[1]);
#if defined(SUBGRID)
/*
        svisn = Svisx(ans);
        svist1 = Svisy(ans);
        Svisx(ans) = (svisn*nor[0]) + (svist1*tan1[0]);
        Svisy(ans) = (svisn*nor[1]) + (svist1*tan1[1]);
*/
#endif /* defined SUBGRID */
}

LOCAL void one_side_parab_tan_solver3d(
        float           ds,
        float           dt,
        Tan_stencil     *sten,
        Locstate        ans,
        Front           *fr)
{
        RECT_GRID      *gr = fr->rect_grid;
        int            i, j, dim = gr->dim;
        Locstate       *state;
        Locstate       *states = sten->states;
        Locstate         sll;
        Locstate          sr;
        Locstate         srr;

        if (is_obstacle_state(ans))
          return;

        sll  =  states[-1];
        sr   =  states[0];
        srr  =  states[1];

        float  nor[3];
        float  norsll, norsr, norsrr;
        float  tan1sll, tan1sr, tan1srr;
        float  tan2sll, tan2sr, tan2srr;
        float  tan1[3], tan2[3];
        norsll = norsr =  norsrr = 0.0;
        tan1sll = tan1sr =  tan1srr = 0.0;
        tan2sll = tan2sr =  tan2srr = 0.0;
        normal(sten->p[0],sten->hse[0],sten->hs[0],nor,fr);
        for ( i = 0; i < dim; ++i)
        {
           tan1[i] = sten->dir[i];
           tan2[i] = sten->dir1[i];
        }
        for ( i = 0; i < dim; ++i)
        {
            norsll += nor[i] * vel(i,sll);
            norsr += nor[i] * vel(i,sr);
            norsrr += nor[i] * vel(i,srr);
        }

        for ( i = 0; i < dim; ++i)
        {
            tan1sll += tan1[i] * vel(i,sll);
            tan1sr += tan1[i] * vel(i,sr);
            tan1srr += tan1[i] * vel(i,srr);
            tan2sll += tan2[i] * vel(i,sll);
            tan2sr += tan2[i] * vel(i,sr);
            tan2srr += tan2[i] * vel(i,srr);
        }

        if(fr->tan_sec == NO)
        {
            int    N = Params(sr)->n_comps;
            float  enthsll, enthsl, enthsr, enthsrr;
            float  pdenssll[N], pdenssl[N], pdenssr[N], pdenssrr[N];
            float  dvt1nt1, dvnnt1, dvnt1t1, dvt1t1t1, dvt2t1t1, dvnt1;
            float  cit1, dcit1, dtt1t1, dtst1t1, dmdt1t1[N], dmdst1t1[N], energymdt1t1;

            dvt1nt1 = (Dvt1n(srr) - Dvt1n(sll))/(2*ds);
            dvnnt1 = (Dvnn(srr) - Dvnn(sll))/(2*ds);
            dvnt1t1 = (norsrr - (2*norsr) + norsll)/(ds*ds);
            dvt1t1t1 = (tan1srr - (2*tan1sr) + tan1sll)/(ds*ds);
            dvt2t1t1 = (tan2srr - (2*tan2sr) + tan2sll)/(ds*ds);
            dvnt1 = (norsrr - norsll)/(2*ds);
            Dvt1t1(ans) = (tan1srr - tan1sll)/(2*ds);   
            Dvt2t1(ans) = (tan2srr - tan2sll)/(2*ds);
            dtt1t1 = (temperature(srr) - (2*temperature(sr)) + temperature(sll))/(ds*ds);
#if defined(SUBGRID)
/*
            C(ans) = C(sr);
            CI(ans) = CI(sr);
            PRT(ans) = PRT(sr);
            SCT(ans) = SCT(sr);
            SCT0(ans) = SCT0(sr);
            STR(ans) = STR(sr);
            CSt1(ans) = -2*C(sr)*ds*ds*Dens(sr)*STR(sr);
            cit1 = (2.0/3.0)*CI(sr)*ds*ds*Dens(sr)*STR(sr)*STR(sr);
            DCSt1(ans) = ( (-2*C(srr)*ds*ds*Dens(srr)*STR(srr))
                         - (-2*C(sll)*ds*ds*Dens(sll)*STR(sll)) )/(2*ds);
            dcit1 = ( ((2.0/3.0)*CI(srr)*ds*ds*Dens(srr)*STR(srr)*STR(srr))
                         - ((2.0/3.0)*CI(sll)*ds*ds*Dens(sll)*STR(sll)*STR(sll)) )/(2*ds);

            dtst1t1 = ((((Dens(srr)*C_P(srr)*PRT(srr)*ds*ds*STR(srr))
                          - (Dens(sll)*C_P(sll)*PRT(sll)*ds*ds*STR(sll)))/(2*ds))
                          *((temperature(srr) - temperature(sll))/(2*ds)))
                          + ((Dens(sr)*C_P(sr)*PRT(sr)*ds*ds*STR(sr))*(dtt1t1));
*/
#endif /* defined SUBGRID */

            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                enthsll = specific_internal_energy_h_species(sll)
                          - specific_internal_energy_l_species(sll);
                enthsr = specific_internal_energy_h_species(sr)
                         - specific_internal_energy_l_species(sr);
                enthsrr = specific_internal_energy_h_species(srr)
                          - specific_internal_energy_l_species(srr);

                if(Params(sr)->n_comps != 1)
                {
                    for(i = 0; i < Params(sr)->n_comps; i++)
                    {
                        pdenssrr[i] = pdens(srr)[i]/Dens(srr);
                        pdenssr[i] = pdens(sr)[i]/Dens(sr);
                        pdenssll[i] = pdens(sll)[i]/Dens(sll);

                        dmdt1t1[i] = (Dens(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(ds*ds))
                                  + (((Dens(srr)) - (Dens(sll)))/(2*ds))*((pdenssrr[i] - pdenssll[i])/(2*ds));
#if defined(SUBGRID)
/*
                       dmdst1t1[i] = ((((Dens(srr)*SCT(srr)*ds*ds*STR(srr))
                                      - (Dens(sll)*SCT(sll)*ds*ds*STR(sll)))/(2*ds))
                                      *((pdenssrr[i] - pdenssll[i])/(2*ds)))
                                      + ((Dens(sr)*SCT(sr)*ds*ds*STR(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(ds*ds)));
*/
#endif /* defined SUBGRID */
                    }
                 }

                 energymdt1t1 = (enthsr*Dens(sr))*((pdenssrr[1] - (2*pdenssr[1]) + pdenssll[1])/(ds*ds))
                                  + (((enthsrr*Dens(srr)) - (enthsll*Dens(sll)))/(2*ds))*((pdenssrr[1] - pdenssll[1])/(2*ds));
            }

            Mvisx(ans) = ((1.0/3.0)*Dvnnn(ans)) - ((2.0/3.0)*dvt1nt1)         
                         + dvt1nt1 + dvnt1t1;                           
            Mvisy(ans) = ((1.0/3.0)*dvt1t1t1) - ((2.0/3.0)*dvnnt1)
                         + dvnnt1 + Dvt1nn(ans);
            Mvisz(ans) = Dvt2nn(ans) + dvt2t1t1;
            Evis(ans) = ((((1.0/3.0)*Dvnnn(ans)) - ((2.0/3.0)*dvt1nt1)   
                        + dvt1nt1 + dvnt1t1)*norsr) + ((((1.0/3.0)*dvt1t1t1) 
                        - ((2.0/3.0)*dvnnt1) + dvnnt1 + Dvt1nn(ans) )*tan1sr) 
                        + ((Dvt2nn(ans) + dvt2t1t1)*tan2sr)   
                        + ((((1.0/3.0)*Dvnn(ans)) 
                        - ((2.0/3.0)*Dvt1t1(ans)))*(Dvnn(ans))) 
                        + ((((1.0/3.0)*Dvt1t1(ans)) 
                        - ((2.0/3.0)*Dvnn(ans)))*(Dvt1t1(ans)))   
                        + ((Dvt1n(ans)+dvnt1)*dvnt1) 
                        + ((dvnt1+Dvt1n(ans))*Dvt1n(ans)); 
            Ethermal(ans) = Ethermal(ans) + dtt1t1;
#if defined(SUBGRID)
/*
            Svisx(ans) = DCSn(ans)*(((2.0/3.0)*Dvnn(ans)) - ((1.0/3.0)*Dvt1t1(ans)))                                                                 
                         + (CSn(ans)*(((2.0/3.0)*Dvnnn(ans)) - ((1.0/3.0)*dvt1nt1))) 
                         + DCIn(ans)
                         + (DCSt1(ans)*(((1.0/2.0)*Dvt1n(ans)) + ((1.0/2.0)*dvnt1)))
                         + (CSt1(ans)*(((1.0/2.0)*dvt1nt1) + ((1.0/2.0)*dvnt1t1)))
                         + dcit1;
            Svisy(ans) = DCSt1(ans)*(((2.0/3.0)*Dvt1t1(ans)) - ((1.0/3.0)*Dvnn(ans))) 
                         + (CSt1(ans)*(((2.0/3.0)*dvt1t1t1) - ((1.0/3.0)*dvnnt1))) 
                         + dcit1
                         + (DCSn(ans)*(((1.0/2.0)*dvnt1) + ((1.0/2.0)*Dvt1n(ans))))
                         + (CSn(ans)*(((1.0/2.0)*dvnnt1) + ((1.0/2.0)*Dvt1nn(ans))))
                         + DCIn(ans);
            Svisz(ans) = (DCSt1(ans)*((1.0/2.0)*Dvt2t1(ans)))  
                         + (CSt1(ans)*((1.0/2.0)*dvt2t1t1))   
                         + dcit1
                         + (DCSn(ans)*((1.0/2.0)*Dvt2n(ans)))  
                         + (CSn(ans)*((1.0/2.0)*Dvt2nn(ans)))  
                         + DCIn(ans) ;
            Sethermal(ans) = Sethermal(ans) + dtst1t1;    
*/
#endif /* defined SUBGRID */

            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                if(Params(sr)->n_comps != 1)
                {
                    for(i = 0; i < Params(sr)->n_comps; i++)
                    {
                        Cmd(ans)[i] =  Cmd(ans)[i] + dmdt1t1[i]; 
#if defined(SUBGRID)
                        //Scmd(ans)[i] = Scmd(ans)[i] + dmdst1t1[i]; 
#endif /* defined SUBGRID */
                    }
                 }
                 Emd(ans) = Emd(ans) + energymdt1t1; 
            }
        }  
        else 
        {
            int    N = Params(sr)->n_comps;
            float  enthsll, enthsl, enthsr, enthsrr;
            float  pdenssll[N], pdenssl[N], pdenssr[N], pdenssrr[N];
            float  dvt2nt2, dvnnt2, dvnt2t2, dvt2t1t2, dvt1t1t2, dvt1t2t2;
            float  dvt2t2t2, dvnt2, dvt2t2, dvt1t2;
            float  cst2, cit2, dcst2, dcit2, dtt2t2, dtst2t2, dmdt2t2[N], dmdst2t2[N], energymdt2t2;

            dvt2nt2 = (Dvt2n(srr) - Dvt2n(sll))/(2*ds);
            dvnnt2 = (Dvnn(srr) - Dvnn(sll))/(2*ds);
            dvnt2t2 = (norsrr - (2*norsr) + norsll)/(ds*ds);
            dvt2t1t2 = (Dvt2t1(srr) - Dvt2t1(sll))/(2*ds);
            dvt1t1t2 = (Dvt1t1(srr) - Dvt1t1(sll))/(2*ds);
            dvt1t2t2 = (tan1srr - (2*tan1sr) + tan1sll)/(ds*ds);
            dvt2t2t2 = (tan2srr - (2*tan2sr) + tan2sll)/(ds*ds);
            dvnt2 = (norsrr - norsll)/(2*ds);
            dvt2t2 = (tan2srr - tan2sll)/(2*ds);
            dvt1t2 = (tan1srr - tan1sll)/(2*ds);
            dtt2t2 = (temperature(srr) - (2*temperature(sr)) + temperature(sll))/(ds*ds);
#if defined(SUBGRID)
/*
            C(ans) = C(sr);
            CI(ans) = CI(sr);
            PRT(ans) = PRT(sr);
            SCT(ans) = SCT(sr);
            SCT0(ans) = SCT0(sr);
            STR(ans) = STR(sr);
            cst2 = -2*C(sr)*ds*ds*Dens(sr)*STR(sr);
            cit2 = (2.0/3.0)*CI(sr)*ds*ds*Dens(sr)*STR(sr)*STR(sr);
            dcst2 = ( (-2*C(srr)*ds*ds*Dens(srr)*STR(srr))
                         - (-2*C(sll)*ds*ds*Dens(sll)*STR(sll)) )/(2*ds);
            dcit2 = ( ((2.0/3.0)*CI(srr)*ds*ds*Dens(srr)*STR(srr)*STR(srr))
                         - ((2.0/3.0)*CI(sll)*ds*ds*Dens(sll)*STR(sll)*STR(sll)) )/(2*ds);
            dtst2t2 = ((((Dens(srr)*C_P(srr)*PRT(srr)*ds*ds*STR(srr))
                          - (Dens(sll)*C_P(sll)*PRT(sll)*ds*ds*STR(sll)))/(2*ds))
                          *((temperature(srr) - temperature(sll))/(2*ds)))
                          + ((Dens(sr)*C_P(sr)*PRT(sr)*ds*ds*STR(sr))*(dtt2t2));
*/
#endif /* defined SUBGRID */            

            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                enthsll = specific_internal_energy_h_species(sll)
                          - specific_internal_energy_l_species(sll);
                enthsr = specific_internal_energy_h_species(sr)
                         - specific_internal_energy_l_species(sr);
                enthsrr = specific_internal_energy_h_species(srr)
                          - specific_internal_energy_l_species(srr);

                if(Params(sr)->n_comps != 1)
                {
                    for(i = 0; i < Params(sr)->n_comps; i++)
                    {
                        pdenssrr[i] = pdens(srr)[i]/Dens(srr);
                        pdenssr[i] = pdens(sr)[i]/Dens(sr);
                        pdenssll[i] = pdens(sll)[i]/Dens(sll);

                        dmdt2t2[i] = (Dens(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(ds*ds))
                                  + (((Dens(srr)) - (Dens(sll)))/(2*ds))*((pdenssrr[i] - pdenssll[i])/(2*ds));
#if defined(SUBGRID)
/*
                        dmdst2t2[i] = ((((Dens(srr)*SCT(srr)*ds*ds*STR(srr))
                                      - (Dens(sll)*SCT(sll)*ds*ds*STR(sll)))/(2*ds))
                                      *((pdenssrr[i] - pdenssll[i])/(2*ds)))
                                      + ((Dens(sr)*SCT(sr)*ds*ds*STR(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(ds*ds)));
*/
#endif /* defined SUBGRID */                    
                     }
                 }

                 energymdt2t2 = (enthsr*Dens(sr))*((pdenssrr[1] - (2*pdenssr[1]) + pdenssll[1])/(ds*ds))
                                  + (((enthsrr*Dens(srr)) - (enthsll*Dens(sll)))/(2*ds))*((pdenssrr[1] - pdenssll[1])/(2*ds));
            }

            Mvisx(ans) = Mvisx(ans) - ((2.0/3.0)*dvt2nt2) + dvt2nt2 + dvnt2t2;
            Mvisy(ans) = Mvisy(ans) - ((2.0/3.0)*dvt2t1t2) + dvt2t1t2 + dvt1t2t2;
            Mvisz(ans) = Mvisz(ans)
                         + ((1.0/3.0)*dvt2t2t2) - ((2.0/3.0)*dvnnt2) - ((2.0/3.0)*dvt1t1t2);
            Evis(ans) = Evis(ans)
                         + ( ((- ((2.0/3.0)*dvt2nt2) + dvt2nt2 + dvnt2t2)*norsr)
                         + ((- ((2.0/3.0)*dvt2t1t2) + dvt2t1t2 + dvt1t2t2)*tan1sr)
                         + ((((1.0/3.0)*dvt2t2t2) - ((2.0/3.0)*dvnnt2)
                         - ((2.0/3.0)*dvt1t1t2) + dvnnt2 + dvt1t1t2)*tan2sr) )
                         - (((2.0/3.0)*dvt2t2)*(Dvnn(ans))) - (((2.0/3.0)*dvt2t2)*Dvt1t1(ans))
                         + ((((1.0/3.0)*dvt2t2) - ((2.0/3.0)*Dvt1t1(ans))
                         - ((2.0/3.0)*Dvnn(ans)))*(dvt2t2))
                         + ((Dvt2n(ans)+dvnt2)*dvnt2) + ((Dvt2t1(ans)+dvt1t2)*dvt1t2)
                         + ((dvnt2+Dvt2n(ans))*Dvt2n(ans)) + ((dvt1t2+Dvt2t1(ans))*Dvt2t1(ans));

            Ethermal(ans) = Ethermal(ans) + dtt2t2;

#if defined(SUBGRID)
/*
            Svisx(ans) = Svisx(ans) + (DCSn(ans)*(((1.0/3.0)*dvt2t2))) 
                         + ((CSn(ans)*(((1.0/3.0)*dvt2nt2))))
                         + (dcst2*(((1.0/2.0)*Dvt2n(ans)) + ((1.0/2.0)*dvnt2)))
                         + (cst2*(((1.0/2.0)*dvt2nt2) + ((1.0/2.0)*dvnt2t2)))
                         + dcit2 ;
            Svisy(ans) = Svisy(ans) + (DCSt1(ans)*((1.0/3.0)*dvt2t2))
                         + (CSt1(ans)*((1.0/3.0)*dvt2t1t2))
                         + (dcst2*(((1.0/2.0)*Dvt2t1(ans)) + ((1.0/2.0)*dvt1t2)))
                         + (cst2*(((1.0/2.0)*dvt2t1t2) + ((1.0/2.0)*dvt1t2t2)))
                         + dcit2 ;
            Svisz(ans) = Svisz(ans) + (DCSt1(ans)*((1.0/2.0)*dvt1t2))
                         + (CSt1(ans)*((1.0/2.0)*dvt1t1t2)) 
                         + (DCSn(ans)*((1.0/2.0)*dvnt2))
                         + (CSn(ans)*((1.0/2.0)*dvnnt2)) 
                         + dcit2 + dcst2*(((2.0/3.0)*dvt2t2) - (((1.0/3.0)*Dvt1t1(ans)) - (((1.0/3.0)*Dvnn(ans)))))
                         + (cst2*(((2.0/3.0)*dvt2t2t2) - (((1.0/3.0)*dvnnt2) - (((1.0/3.0)*dvt1t1t2)))));
             Sethermal(ans) = Sethermal(ans) + dtst2t2;
*/
#endif /* defined SUBGRID */

             if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
             {
                 if(Params(sr)->n_comps != 1)
                 {
                     for(i = 0; i < Params(sr)->n_comps; i++)
                     {
                         Cmd(ans)[i] =  Cmd(ans)[i] +dmdt2t2[i];
#if defined(SUBGRID)
                         //Scmd(ans)[i] = Scmd(ans)[i] + dmdst2t2[i];
#endif /* defined SUBGRID */
                     }
                 }
                 Emd(ans) = Emd(ans) + energymdt2t2;
             }

             float mvisn, mvist1, mvist2;
             float svisn, svist1, svist2;
             mvisn = Mvisx(ans);
             mvist1 = Mvisy(ans);
             mvist2 = Mvisz(ans);
             Mvisx(ans) = (mvisn*nor[0]) + (mvist1*tan1[0]) + (mvist2*tan2[0]);
             Mvisy(ans) = (mvisn*nor[1]) + (mvist1*tan1[1]) + (mvist2*tan2[1]);
             Mvisz(ans) = (mvisn*nor[2]) + (mvist1*tan1[2]) + (mvist2*tan2[2]);
#if defined(SUBGRID) 
/*
             svisn = Svisx(ans);
             svist1 = Svisy(ans);
             svist2 = Svisz(ans);
             Svisx(ans) = (svisn*nor[0]) + (svist1*tan1[0]) + (svist2*tan2[0]);
             Svisy(ans) = (svisn*nor[1]) + (svist1*tan1[1]) + (svist2*tan2[1]);
             Svisz(ans) = (svisn*nor[2]) + (svist1*tan1[2]) + (svist2*tan2[2]);
*/
#endif /* defined SUBGRID */        
         } 
}

EXPORT  void    parab_nor_solver(
        WSSten          *sten,
        Locstate        ansl,
        Locstate        ansr)
{
        Locstate         sll;
        Locstate          sl;
        Locstate          sr;
        Locstate         srr;
        sll  =  sten->sl[1];
        sl   =  sten->sl[0];
        sr   =  sten->sr[0];
        srr  =  sten->sr[1];
        int             w_type = sten->w_type;
        float           dn = sten->dn;
        float           dt = sten->dt;
        Front           *front = sten->front;
        int             i,j, dim = front->interf->dim;
        float           nor[dim];
        float  norsll, norsl, norsr, norsrr;
        float  tan1sll, tan1sl, tan1sr, tan1srr;
        float  tan2sll, tan2sl, tan2sr, tan2srr;
        float  tan1[3], tan2[3];
        static Tparams tp[2];
        float  tngt[MAXD];

        if (is_obstacle_state(ansl))
          return;
        if (is_obstacle_state(ansr))
          return;

        if(dim == 2)
        {
        tangent(sten->p,Bond_of_hse(sten->hse),Curve_of_hs(sten->hs),tngt,sten->front);
        normal(sten->p,sten->hse,sten->hs,nor,sten->front);
        for ( i = 0; i < dim; ++i)
        {
           tan1[i] = tngt[i];
        }
        norsll = norsl = norsr =  norsrr = 0.0;
        tan1sll = tan1sl = tan1sr =  tan1srr = 0.0;
        for ( i = 0; i < dim; ++i)
        {
            norsll += nor[i] * vel(i,sll);
            norsl += nor[i] * vel(i,sl);
            norsr += nor[i] * vel(i,sr);
            norsrr += nor[i] * vel(i,srr);
        }

        for ( i = 0; i < dim; ++i)
        {
            tan1sll += tan1[i] * vel(i,sll);
            tan1sl += tan1[i] * vel(i,sl);
            tan1sr += tan1[i] * vel(i,sr);
            tan1srr += tan1[i] * vel(i,srr);
        }
       
        if (w_type == CONTACT)
        {
            int    N = Params(sr)->n_comps;
            float  enthsll, enthsl, enthsr, enthsrr;
            float  pdenssll[N], pdenssl[N], pdenssr[N], pdenssrr[N];

            Dvnnn(ansl) = (norsrr - (2*norsl) + norsll)/(dn*dn);
            Dvnnn(ansr) = (norsrr - (2*norsr) + norsll)/(dn*dn);
            Dvnn(ansl) = (norsrr - norsll)/(2*dn);
            Dvnn(ansr) = Dvnn(ansl);
            Dvt1n(ansl) = (tan1srr - tan1sll)/(2*dn);
            Dvt1n(ansr) = (tan1srr - tan1sll)/(2*dn);
            Dvt1nn(ansl) = (tan1srr - (2*tan1sl) + tan1sll)/(dn*dn);
            Dvt1nn(ansr) = (tan1srr - (2*tan1sr) + tan1sll)/(dn*dn);
            Ethermal(ansl) = (temperature(srr) - (2*temperature(sl)) + temperature(sll))/(dn
*dn);
            Ethermal(ansr) = (temperature(srr) - (2*temperature(sr)) + temperature(sll))/(dn
*dn);
#if defined(SUBGRID)
/*
            C(ansl) = C(sl);
            CI(ansl) = CI(sl);
            PRT(ansl) = PRT(sl);
            SCT(ansl) = SCT(sl);
            SCT0(ansl) = SCT0(sl);
            STR(ansl) = STR(sl);
            C(ansr) = C(sr);
            CI(ansr) = CI(sr);
            PRT(ansr) = PRT(sr);
            SCT(ansr) = SCT(sr);
            SCT0(ansr) = SCT0(sr);
            STR(ansr) = STR(sr);
            CSn(ansl) = -2*C(sl)*dn*dn*Dens(sl)*STR(sl);
            CSn(ansr) = -2*C(sr)*dn*dn*Dens(sr)*STR(sr);
            CIn(ansl) = (2.0/2.0)*CI(sl)*dn*dn*Dens(sl)*STR(sl)*STR(sl);
            CIn(ansr) = (2.0/2.0)*CI(sr)*dn*dn*Dens(sr)*STR(sr)*STR(sr);
            DCSn(ansl) = ( (-2*C(srr)*dn*dn*Dens(srr)*STR(srr))
                         - (-2*C(sll)*dn*dn*Dens(sll)*STR(sll)) )/(2*dn);
            DCSn(ansr) = DCSn(ansl);
            DCIn(ansl) = ( ((2.0/2.0)*CI(srr)*dn*dn*Dens(srr)*STR(srr)*STR(srr))
                         - ((2.0/2.0)*CI(sll)*dn*dn*Dens(sll)*STR(sll)*STR(sll)) )/(2*dn);
            DCIn(ansr) = DCIn(ansl);
            Sethermal(ansl) = ((((Dens(srr)*C_P(srr)*PRT(srr)*dn*dn*STR(srr))
                          - (Dens(sll)*C_P(sll)*PRT(sll)*dn*dn*STR(sll)))/(2*dn))
                          *((temperature(srr) - temperature(sll))/(2*dn)))
                          + ((Dens(sl)*C_P(sl)*PRT(sl)*dn*dn*STR(sl))*(Ethermal(ansl)));
            Sethermal(ansr) = ((((Dens(srr)*C_P(srr)*PRT(srr)*dn*dn*STR(srr))
                          - (Dens(sll)*C_P(sll)*PRT(sll)*dn*dn*STR(sll)))/(2*dn))
                          *((temperature(srr) - temperature(sll))/(2*dn)))
                          + ((Dens(sr)*C_P(sr)*PRT(sr)*dn*dn*STR(sr))*(Ethermal(ansr)));
*/
#endif /* defined SUBGRID */

            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                enthsll = specific_internal_energy_h_species(sll)
                          - specific_internal_energy_l_species(sll);
                enthsl = specific_internal_energy_h_species(sl)
                          - specific_internal_energy_l_species(sl);
                enthsr = specific_internal_energy_h_species(sr)
                         - specific_internal_energy_l_species(sr);
                enthsrr = specific_internal_energy_h_species(srr)
                          - specific_internal_energy_l_species(srr);

                if(Params(sr)->n_comps != 1)
                {
                    for(i = 0; i < Params(sr)->n_comps; i++)
                    {
                        pdenssrr[i] = pdens(srr)[i]/Dens(srr);
                        pdenssr[i] = pdens(sr)[i]/Dens(sr);
                        pdenssl[i] = pdens(sl)[i]/Dens(sl);
                        pdenssll[i] = pdens(sll)[i]/Dens(sll);

                        Cmd(ansl)[i] = (Dens(sl))*((pdenssrr[i] - (2*pdenssl[i]) + pdenssll[i])/(dn*dn))
                                  + (((Dens(srr)) - (Dens(sll)))/(2*dn))*((pdenssrr[i] - pdenssll[i])/(2*dn));
                        Cmd(ansr)[i] = (Dens(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(dn*dn))
                                  + (((Dens(srr)) - (Dens(sll)))/(2*dn))*((pdenssrr[i] - pdenssll[i])/(2*dn));
#if defined(SUBGRID) 
/*
                        Scmd(ansl)[i] = ((((Dens(srr)*SCT(srr)*dn*dn*STR(srr))
                                      - (Dens(sll)*SCT(sll)*dn*dn*STR(sll)))/(2*dn))
                                      *((pdenssrr[i] - pdenssll[i])/(2*dn)))
                                      + ((Dens(sl)*SCT(sl)*dn*dn*STR(sl))*((pdenssrr[i] - (2*pdenssl[i]) + pdenssll[i])/(dn*dn)));
                        Scmd(ansr)[i] = ((((Dens(srr)*SCT(srr)*dn*dn*STR(srr))
                                      - (Dens(sll)*SCT(sll)*dn*dn*STR(sll)))/(2*dn))
                                      *((pdenssrr[i] - pdenssll[i])/(2*dn)))
                                      + ((Dens(sr)*SCT(sr)*dn*dn*STR(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(dn*dn)));
*/
#endif /* defined SUBGRID */                    
                     }
                 }

                 Emd(ansl) = (enthsl*Dens(sl))*((pdenssrr[1] - (2*pdenssl[1]) + pdenssll[1])/(dn*dn))
                                  + (((enthsrr*Dens(srr)) - (enthsll*Dens(sll)))/(2*dn))*((pdenssrr[1] - pdenssll[1])/(2*dn));
                 Emd(ansr) = (enthsr*Dens(sr))*((pdenssrr[1] - (2*pdenssr[1]) + pdenssll[1])/(dn*dn))
                                  + (((enthsrr*Dens(srr)) - (enthsll*Dens(sll)))/(2*dn))*((pdenssrr[1] - pdenssll[1])/(2*dn));
            }
        }
        }

        if(dim == 3)
        {
        set_up_tangent_params(sten->front,sten->p,sten->hse,sten->hs,tp);
        normal(sten->p,sten->hse,sten->hs,nor,sten->front);
        for ( i = 0; i < dim; ++i)
        {
           tan1[i] = tp[0].tan[i];
           tan2[i] = tp[1].tan[i];
        }
        norsll = norsl = norsr =  norsrr = 0.0;
        tan1sll = tan1sl = tan1sr =  tan1srr = 0.0;
        tan2sll = tan2sl = tan2sr =  tan2srr = 0.0;
        for ( i = 0; i < dim; ++i)
        {
            norsll += nor[i] * vel(i,sll);
            norsl += nor[i] * vel(i,sl);
            norsr += nor[i] * vel(i,sr);
            norsrr += nor[i] * vel(i,srr);
        }

        for ( i = 0; i < dim; ++i)
        {
            tan1sll += tan1[i] * vel(i,sll);
            tan1sl += tan1[i] * vel(i,sl);
            tan1sr += tan1[i] * vel(i,sr);
            tan1srr += tan1[i] * vel(i,srr);
            tan2sll += tan2[i] * vel(i,sll);
            tan2sl += tan2[i] * vel(i,sl);
            tan2sr += tan2[i] * vel(i,sr);
            tan2srr += tan2[i] * vel(i,srr);
        }

        if (w_type == CONTACT)
        {
            int    N = Params(sr)->n_comps;
            float  enthsll, enthsl, enthsr, enthsrr;
            float  pdenssll[N], pdenssl[N], pdenssr[N], pdenssrr[N];

            Dvnnn(ansl) = (norsrr - (2*norsl) + norsll)/(dn*dn);
            Dvnnn(ansr) = (norsrr - (2*norsr) + norsll)/(dn*dn);
            Dvnn(ansl) = (norsrr - norsll)/(2*dn);
            Dvnn(ansr) = Dvnn(ansl);
            Dvt1n(ansl) = (tan1srr - tan1sll)/(2*dn);
            Dvt1n(ansr) = (tan1srr - tan1sll)/(2*dn);
            Dvt2n(ansl) = (tan2srr - tan2sll)/(2*dn);
            Dvt2n(ansr) = (tan2srr - tan2sll)/(2*dn);
            Dvt1nn(ansl) = (tan1srr - (2*tan1sl) + tan1sll)/(dn*dn);
            Dvt1nn(ansr) = (tan1srr - (2*tan1sr) + tan1sll)/(dn*dn);
            Dvt2nn(ansl) = (tan2srr - (2*tan2sl) + tan2sll)/(dn*dn);
            Dvt2nn(ansr) = (tan2srr - (2*tan2sr) + tan2sll)/(dn*dn);
            Ethermal(ansl) = (temperature(srr) - (2*temperature(sl)) + temperature(sll))/(dn
*dn);
            Ethermal(ansr) = (temperature(srr) - (2*temperature(sr)) + temperature(sll))/(dn
*dn);
#if defined(SUBGRID)
/*
            C(ansl) = C(sl);
            CI(ansl) = CI(sl);
            PRT(ansl) = PRT(sl);
            SCT(ansl) = SCT(sl);
            SCT0(ansl) = SCT0(sl);
            STR(ansl) = STR(sl);
            C(ansr) = C(sr);
            CI(ansr) = CI(sr);
            PRT(ansr) = PRT(sr);
            SCT(ansr) = SCT(sr);
            SCT0(ansr) = SCT0(sr);
            STR(ansr) = STR(sr);
            CSn(ansl) = -2*C(sl)*dn*dn*Dens(sl)*STR(sl);
            CSn(ansr) = -2*C(sr)*dn*dn*Dens(sr)*STR(sr);
            CIn(ansl) = (2.0/3.0)*CI(sl)*dn*dn*Dens(sl)*STR(sl)*STR(sl); 
            CIn(ansr) = (2.0/3.0)*CI(sr)*dn*dn*Dens(sr)*STR(sr)*STR(sr);
            DCSn(ansl) = ( (-2*C(srr)*dn*dn*Dens(srr)*STR(srr))
                         - (-2*C(sll)*dn*dn*Dens(sll)*STR(sll)) )/(2*dn);
            DCSn(ansr) = DCSn(ansl);
            DCIn(ansl) = ( ((2.0/3.0)*CI(srr)*dn*dn*Dens(srr)*STR(srr)*STR(srr))
                         - ((2.0/3.0)*CI(sll)*dn*dn*Dens(sll)*STR(sll)*STR(sll)) )/(2*dn);
            DCIn(ansr) = DCIn(ansl);
            Sethermal(ansl) = ((((Dens(srr)*C_P(srr)*PRT(srr)*dn*dn*STR(srr)) 
                          - (Dens(sll)*C_P(sll)*PRT(sll)*dn*dn*STR(sll)))/(2*dn))
                          *((temperature(srr) - temperature(sll))/(2*dn)))
                          + ((Dens(sl)*C_P(sl)*PRT(sl)*dn*dn*STR(sl))*(Ethermal(ansl)));
            Sethermal(ansr) = ((((Dens(srr)*C_P(srr)*PRT(srr)*dn*dn*STR(srr))
                          - (Dens(sll)*C_P(sll)*PRT(sll)*dn*dn*STR(sll)))/(2*dn))
                          *((temperature(srr) - temperature(sll))/(2*dn)))
                          + ((Dens(sr)*C_P(sr)*PRT(sr)*dn*dn*STR(sr))*(Ethermal(ansr)));
*/
#endif /* defined SUBGRID */

            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                enthsll = specific_internal_energy_h_species(sll)
                          - specific_internal_energy_l_species(sll);
                enthsl = specific_internal_energy_h_species(sl)
                          - specific_internal_energy_l_species(sl);
                enthsr = specific_internal_energy_h_species(sr)
                         - specific_internal_energy_l_species(sr);
                enthsrr = specific_internal_energy_h_species(srr)
                          - specific_internal_energy_l_species(srr);

                if(Params(sr)->n_comps != 1)
                {
                    for(i = 0; i < Params(sr)->n_comps; i++)
                    {
                        pdenssrr[i] = pdens(srr)[i]/Dens(srr);
                        pdenssr[i] = pdens(sr)[i]/Dens(sr);
                        pdenssl[i] = pdens(sl)[i]/Dens(sl);
                        pdenssll[i] = pdens(sll)[i]/Dens(sll);

                        Cmd(ansl)[i] = (Dens(sl))*((pdenssrr[i] - (2*pdenssl[i]) + pdenssll[i])/(dn*dn))
                                  + (((Dens(srr)) - (Dens(sll)))/(2*dn))*((pdenssrr[i] - pdenssll[i])/(2*dn));
                        Cmd(ansr)[i] = (Dens(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(dn*dn))
                                  + (((Dens(srr)) - (Dens(sll)))/(2*dn))*((pdenssrr[i] - pdenssll[i])/(2*dn));
#if defined(SUBGRID)
/*
                        Scmd(ansl)[i] = ((((Dens(srr)*SCT(srr)*dn*dn*STR(srr))
                                      - (Dens(sll)*SCT(sll)*dn*dn*STR(sll)))/(2*dn))
                                      *((pdenssrr[i] - pdenssll[i])/(2*dn)))
                                      + ((Dens(sl)*SCT(sl)*dn*dn*STR(sl))*((pdenssrr[i] - (2*pdenssl[i]) + pdenssll[i])/(dn*dn)));
                        Scmd(ansr)[i] = ((((Dens(srr)*SCT(srr)*dn*dn*STR(srr))
                                      - (Dens(sll)*SCT(sll)*dn*dn*STR(sll)))/(2*dn))
                                      *((pdenssrr[i] - pdenssll[i])/(2*dn)))
                                      + ((Dens(sr)*SCT(sr)*dn*dn*STR(sr))*((pdenssrr[i] - (2*pdenssr[i]) + pdenssll[i])/(dn*dn)));
*/
#endif /* defined SUBGRID */
                    }
                 }

                 Emd(ansl) = (enthsl*Dens(sl))*((pdenssrr[1] - (2*pdenssl[1]) + pdenssll[1])/(dn*dn)) 
                                  + (((enthsrr*Dens(srr)) - (enthsll*Dens(sll)))/(2*dn))*((pdenssrr[1] - pdenssll[1])/(2*dn));
                 Emd(ansr) = (enthsr*Dens(sr))*((pdenssrr[1] - (2*pdenssr[1]) + pdenssll[1])/(dn*dn))
                                  + (((enthsrr*Dens(srr)) - (enthsll*Dens(sll)))/(2*dn))*((pdenssrr[1] - pdenssll[1])/(2*dn));
            }
        }
        }
}

LOCAL   void    set_propagation_bounds(
        Front *fr,
        float *L,
        float *U)
{
        RECT_GRID       *gr = fr->rect_grid;
        int             i,dim = gr->dim;

        /* Set propagation boundaries */
        for (i = 0; i < dim; ++i)
        {
            if (gr->lbuf[i] > 0)
                L[i] = gr->L[i] - 2.0*gr->h[i];
            else
                L[i] = -HUGE_VAL;;
            if (gr->ubuf[i] > 0)
                U[i] = gr->U[i] + 2.0*gr->h[i];
            else
                U[i] = HUGE_VAL;;
        }
}       /* end set_propagation_bounds */

LOCAL   bool    out_of_bound(
        POINT *p,
        float *L,
        float *U,
        int dim)
{
        int i;
        for (i = 0; i < dim; ++i)
            if (Coords(p)[i] < L[i] || Coords(p)[i] > U[i])
                return YES;
        return NO;
}       /* end out_of_bound */

LOCAL   void cal_volume_fraction(
        float           dt,
        Front           *front,
        Wave            *wv)
{
        const int       nn = pp_numnodes();
        RECT_GRID       *rgr = wv->rect_grid;
        float           *GL = rgr->GL, *GU = rgr->GU;
        float           *L = rgr->L, *U = rgr->U;
        int             dim = rgr->dim;
        int             *gmax = rgr->gmax;
        int             i,j,k;
        int             icoords[MAXD];
        static          Pstencil *nsten;
        int             npts,mid;
        int             imin[MAXD],imax[MAXD];

/*      
        if (nsten == NULL)
        {
            stat_scalar(&nsten,sizeof(Pstencil));
            switch (dim)
            {
            case 1:
                bi_array(&nsten->icoords1d,3,1,INT);
                uni_array(&nsten->sts1d,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                    alloc_state(front->interf,&nsten->sts1d[i],front->sizest);
                break;
            case 2:
                tri_array(&nsten->icoords2d,3,3,2,INT);
                bi_array(&nsten->sts2d,3,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                {
                    for (j = 0; j < 3; j++)
                    {
                        alloc_state(front->interf,&nsten->sts2d[i][j],
                                front->sizest);
                    }
                }
                break;
            case 3:
                quad_array(&nsten->icoords3d,3,3,3,3,INT);
                tri_array(&nsten->sts3d,3,3,3,sizeof(Locstate));
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        for (k = 0; k < 3; k++)
                            alloc_state(front->interf,&nsten->sts3d[i][j][k],
                                     front->sizest);
                break;
            }
        }
        nsten->fr = front;
        nsten->wave = wv;
        nsten->npts = npts = wv->npts_sten;
        mid = (int) (npts/2);
        */

        FILE *outfile;
                                           
        Locstate state;

        for (i = 0; i < dim; i++)
        {
            imin[i] = 0;        imax[i] = gmax[i];
        }

        float  dh[3];
        dh[0] = wv->rect_grid->h[0];
        dh[1] = wv->rect_grid->h[1];
        dh[2] = wv->rect_grid->h[2];

        float vol_frac_l, vol_frac_h, fl;
        float new_vol_frac_l, new_vol_frac_h;
        float x_ccoords,y_ccoords,z_ccoords,max_dens,min_dens;
        float new_max_dens,new_min_dens;
        float count, sum_vf1,sum_vf2,sum_vf2old,sum_vf1vf2,sum_vf1vf2old,sum_nvf1,sum_nvf2,sum_nvf2old,sum_nvf1nvf2,sum_nvf1nvf2old,sum_f1,sum_f2,sum_f1f2,mean_nvf1,mean_nvf2,mean_nvf2old,mean_nvf1nvf2,mean_nvf1nvf2old,mean_vf1,mean_vf2,mean_vf2old,mean_vf1vf2,mean_vf1vf2old,mean_f1,mean_f2,mean_f1f2;
        float vf1,vf2,vf2old,vf1vf2,vf1vf2old,nvf1,nvf2,nvf2old,nvf1nvf2,nvf1nvf2old,f1,f2,f1f2,Y1,Y2;

        float unew_max_dens,unew_min_dens;
        float ucount, usum_vf1,usum_vf2,usum_vf2old,usum_vf1vf2,usum_vf1vf2old,usum_nvf1,usum_nvf2,usum_nvf2old,usum_nvf1nvf2,usum_nvf1nvf2old,usum_f1,usum_f2,usum_f1f2,umean_nvf1,umean_nvf2,umean_nvf2old,umean_nvf1nvf2,umean_nvf1nvf2old,umean_vf1,umean_vf2,umean_vf2old,umean_vf1vf2,umean_vf1vf2old,umean_f1,umean_f2,umean_f1f2;
        float uvf1,uvf2,uvf2old,uvf1vf2,uvf1vf2old,unvf1,unvf2,unvf2old,unvf1nvf2,unvf1nvf2old,uf1,uf2,uf1f2,uY1,uY2;

        max_dens = 0.0;
        min_dens = 1000.0;

        for (k = imin[2]; k < imax[2]; k++)
        {
            icoords[2] = k;
            for (j = imin[1]; j < imax[1]; j++)
            {
                icoords[1] = j;
                for (i = imin[0]; i < imax[0]; i++)
                {
                    icoords[0] = i;
                    state = Rect_state(icoords,wv);
                    if(max_dens < Dens(state))
                        max_dens = Dens(state);
                    if(min_dens > Dens(state))
                        min_dens = Dens(state);
                }
            }
        }

        if (nn > 1)
        {
              pp_global_min(&min_dens,1L);
              pp_global_max(&max_dens,1L);
        }

        count = 0.0;
        sum_vf1 = 0.0;
        sum_vf2 = 0.0;
        sum_vf1vf2 = 0.0;
        sum_f1 = 0.0;
        sum_f2 = 0.0;
        sum_f1f2 = 0.0;

        ucount = 0.0;
        usum_vf1 = 0.0;
        usum_vf2 = 0.0;
        usum_vf1vf2 = 0.0;
        usum_f1 = 0.0;
        usum_f2 = 0.0;
        usum_f1f2 = 0.0;


        for (k = imin[2]; k < imax[2]; k++)
        {
            icoords[2] = k;
            for (j = imin[1]; j < imax[1]; j++)
            {
                icoords[1] = j;
                for (i = imin[0]; i < imax[0]; i++)
                {
                    icoords[0] = i;
                    state = Rect_state(icoords,wv);
                    int comp = Rect_comp(icoords,wv); 
                    z_ccoords = Rect_coords(icoords,wv)[2];
                    //fprintf(stdout,"PRINT %d %e %e %e \n", comp,pdens(state)[0],pdens(state)[1],Dens(state));

                    //if(z_ccoords > 12.0-dh[2] && z_ccoords < 12.0+dh[2])
                    if(z_ccoords > 12.0-dh[2] && z_ccoords < 12.0)
                    {
                    count++;    

                    vf1 = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vf2 = (Dens(state)-min_dens)/(max_dens-min_dens);
                    vf1vf2 = vf1*vf2;
                    Y1 = pdens(state)[0]/Dens(state);
                    Y2 = pdens(state)[1]/Dens(state);
                    //f1 = (Y1/76.45)/((Y1/76.45)+(Y2/18.02));
                    //f2 = (Y2/18.02)/((Y1/76.45)+(Y2/18.02));
                    f1 = (Y1/18.02)/((Y1/18.02)+(Y2/18.02));
                    f2 = (Y2/18.02)/((Y1/18.02)+(Y2/18.02));
                    f1f2 = f1*f2;
                    //fprintf(stdout,"PRINT-1 %d %e %e %e %e  %e %e %e %e\n", comp,pdens(state)[0],pdens(state)[1],Dens(state),f1,f2,Y1/76.45,Y2/18.02,(Y1/76.45)+(Y2/18.02));

                    sum_vf1 += vf1;
                    sum_vf2 += vf2;
                    sum_vf1vf2 += vf1vf2;
                    sum_f1 += f1;
                    sum_f2 += f2;
                    sum_f1f2 += f1f2;
                    }

                    if(z_ccoords > 12.0 && z_ccoords < 12.0+dh[2])
                    {
                    ucount++;

                    uvf1 = (Dens(state)-max_dens)/(min_dens-max_dens);
                    uvf2 = (Dens(state)-min_dens)/(max_dens-min_dens);
                    uvf1vf2 = uvf1*uvf2;
                    uY1 = pdens(state)[0]/Dens(state);
                    uY2 = pdens(state)[1]/Dens(state);
                    //uf1 = (uY1/76.45)/((uY1/76.45)+(uY2/18.02));
                    //uf2 = (uY2/18.02)/((uY1/76.45)+(uY2/18.02));
                    uf1 = (uY1/18.02)/((uY1/18.02)+(uY2/18.02));
                    uf2 = (uY2/18.02)/((uY1/18.02)+(uY2/18.02));
                    uf1f2 = uf1*uf2;
                    //fprintf(stdout,"PRINT-1 %d %e %e %e %e  %e %e %e %e\n", comp,pdens(state)[0],pdens(state)[1],Dens(state),f1,f2,Y1/76.45,Y2/18.02,(Y1/76.45)+(Y2/18.02));

                    usum_vf1 += uvf1;
                    usum_vf2 += uvf2;
                    usum_vf1vf2 += uvf1vf2;
                    usum_f1 += uf1;
                    usum_f2 += uf2;
                    usum_f1f2 += uf1f2;
                    }
                }
            }
        }

        if (nn > 1)
        {
              pp_global_sum(&count,1L);
              pp_global_sum(&sum_vf1,1L);
              pp_global_sum(&sum_vf2,1L);
              pp_global_sum(&sum_vf1vf2,1L);
              pp_global_sum(&sum_f1,1L);
              pp_global_sum(&sum_f2,1L);
              pp_global_sum(&sum_f1f2,1L);

              pp_global_sum(&ucount,1L);
              pp_global_sum(&usum_vf1,1L);
              pp_global_sum(&usum_vf2,1L);
              pp_global_sum(&usum_vf1vf2,1L);
              pp_global_sum(&usum_f1,1L);
              pp_global_sum(&usum_f2,1L);
              pp_global_sum(&usum_f1f2,1L);
        }
        mean_vf1 = sum_vf1/count;
        mean_vf2 = sum_vf2/count;
        mean_vf1vf2 = sum_vf1vf2/count;
        mean_f1 = sum_f1/count;
        mean_f2 = sum_f2/count;
        mean_f1f2 = sum_f1f2/count;

        umean_vf1 = usum_vf1/ucount;
        umean_vf2 = usum_vf2/ucount;
        umean_vf1vf2 = usum_vf1vf2/ucount;
        umean_f1 = usum_f1/ucount;
        umean_f2 = usum_f2/ucount;
        umean_f1f2 = usum_f1f2/ucount;

if(pp_mynode() == 0)
{
    if (front->step == 0)
    {
        outfile = fopen("/nfs/user01/hyulim/RT_simulation/FronTier/src/gas/out_sim_hc/sim_hc_theta.txt","w");
        fprintf(outfile,"    time        theta \n");
    }
    else
        outfile = fopen("/nfs/user01/hyulim/RT_simulation/FronTier/src/gas/out_sim_hc/sim_hc_theta.txt","a");

    float time32,time24,A,g,H;
    A = (2.0-1.0)/(2.0+1.0); 
    g = 0.0981;
    H = 32.0;
    time32 = front->time*sqrt(A*g/H); 
    H = 24.0;
    time24 = front->time*sqrt(A*g/H);

    //fprintf(outfile,"%e %e %e %e %e %e %e %e\n", front->time,time32,time24,mean_nvf1nvf2old/(mean_nvf1*mean_nvf2old),mean_vf1vf2old/(mean_vf1*mean_vf2old),mean_nvf1nvf2/(mean_nvf1*mean_nvf2),mean_vf1vf2/(mean_vf1*mean_vf2),mean_f1f2/(mean_f1*mean_f2));
    fprintf(outfile,"%e %e %e %e %e %e %e %e %e\n", front->time,time32,time24,mean_vf1vf2/(mean_vf1*mean_vf2),umean_vf1vf2/(umean_vf1*umean_vf2),((mean_vf1vf2/(mean_vf1*mean_vf2))+(umean_vf1vf2/(umean_vf1*umean_vf2)))/2.0,mean_f1f2/(mean_f1*mean_f2),umean_f1f2/(umean_f1*umean_f2),((mean_f1f2/(mean_f1*mean_f2))+(umean_f1f2/(umean_f1*umean_f2)))/2.0);

    fclose(outfile);
}
        //fprintf(stdout,"theta %e %e %e\n", mean_nvf1nvf2/(mean_nvf1*mean_nvf2),mean_vf1vf2/(mean_vf1*mean_vf2),mean_f1f2/(mean_f1*mean_f2));


/*
        int recon_step = 50;
        if ( (front->step)%recon_step == 1)
        {
        //volume_fraction
        char filename[200];
        char filename1[200],filename2[200],filename3[200],filename4[200],filename5[200],filename6[200],filename7[200],filename8[200],filename9[200],filename10[200];
        FILE *outfile;
        FILE *outfile1,*outfile2,*outfile3,*outfile4,*outfile5,*outfile6,*outfile7,*outfile8,*outfile9,*outfile10;
        sprintf(filename,"%svolume_fraction-ts%s",get_directory(),right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif
        //sprintf(filename,"%s-liquid.vtk",filename);
        sprintf(filename1,"%s-z1",filename);
        sprintf(filename2,"%s-z2",filename);
        sprintf(filename3,"%s-z3",filename);
        sprintf(filename4,"%s-z4",filename);
        sprintf(filename5,"%s-z5",filename);
        sprintf(filename6,"%s-z6",filename);
        sprintf(filename7,"%s-z7",filename);
        sprintf(filename8,"%s-z8",filename);
        sprintf(filename9,"%s-z9",filename);
        sprintf(filename10,"%s-z10",filename);

        outfile1 = fopen(filename1,"w");
        outfile2 = fopen(filename2,"w");
        outfile3 = fopen(filename3,"w");
        outfile4 = fopen(filename4,"w");
        outfile5 = fopen(filename5,"w");
        outfile6 = fopen(filename6,"w");
        outfile7 = fopen(filename7,"w");
        outfile8 = fopen(filename8,"w");
        outfile9 = fopen(filename9,"w");
        outfile10 = fopen(filename10,"w");

        float vol_frac_l, vol_frac_h, fl;
        float new_vol_frac_l, new_vol_frac_h;
        float x_ccoords,y_ccoords,z_ccoords,max_dens,min_dens;
        float new_max_dens,new_min_dens;

        max_dens = 0.0;
        min_dens = 1000.0;
        new_max_dens = 0.0;
        new_min_dens = 1000.0;
        int compo;

        for (k = imin[2]; k < imax[2]; k++)
        {
            icoords[2] = k;
            for (j = imin[1]; j < imax[1]; j++)
            {
                icoords[1] = j;
                for (i = imin[0]; i < imax[0]; i++)
                {
                    icoords[0] = i;
                    state = Rect_state(icoords,wv);
                    if(max_dens < Dens(state))
                        max_dens = Dens(state);
                    if(min_dens > Dens(state))
                        min_dens = Dens(state);

                    z_ccoords = cell_center(k,2,rgr);
                    if(z_ccoords > 12.0-dh[2] && z_ccoords < 12.0+dh[2])
                    {
                        if(new_max_dens < Dens(state))
                            new_max_dens = Dens(state);
                        if(new_min_dens > Dens(state))
                            new_min_dens = Dens(state);
                    }

                }
            }
        }

        pp_gsync();

        if (nn > 1)
        {
              pp_global_min(&min_dens,1L);
              pp_global_max(&max_dens,1L);
              pp_global_min(&new_min_dens,1L);
              pp_global_max(&new_max_dens,1L);
        }
        for (k = imin[2]; k < imax[2]; k++)
        {
            icoords[2] = k;
            for (j = imin[1]; j < imax[1]; j++)
            {
                icoords[1] = j;
                for (i = imin[0]; i < imax[0]; i++)
                {
                    icoords[0] = i;
                    state = Rect_state(icoords,wv);
                    compo = Rect_comp(icoords,wv);
                    x_ccoords = cell_center(i,0,rgr);
                    y_ccoords = cell_center(j,1,rgr);
                    z_ccoords = cell_center(k,2,rgr);

                    if(z_ccoords > 12.0-(5*dh[2]) && z_ccoords < 12.0-(4*dh[2]))
                    {
                    //fprintf(stdout, "%d %d %d %e %e %e %e %e \n",i,j,k,C[i][j][k],CI[i][j][k],Prt[i][j][k],Sct[i][j][k],Sct0[i][j][k]);
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    //fprintf(stdout, "compo %d C_V %e Gamma   %e %e %e\n",compo,C_V(state),Dens(state),pdens(state)[0],pdens(state)[1]);
                    fprintf(outfile1,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    //fprintf(outfile,"%20.16f %20.16f %20.16f %20.16f %20.16f\n", vol_frac_l,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0-(4*dh[2]) && z_ccoords < 12.0-(3*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile2,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0-(3*dh[2]) && z_ccoords < 12.0-(2*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile3,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0-(2*dh[2]) && z_ccoords < 12.0-(1*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile4,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0-(1*dh[2]) && z_ccoords < 12.0)
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile5,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0 && z_ccoords < 12.0+(1*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile6,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0+(1*dh[2]) && z_ccoords < 12.0+(2*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile7,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0+(2*dh[2]) && z_ccoords < 12.0+(3*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile8,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0+(3*dh[2]) && z_ccoords < 12.0+(4*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile9,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                    if(z_ccoords > 12.0+(4*dh[2]) && z_ccoords < 12.0+(5*dh[2]))
                    {
                    new_vol_frac_l = (Dens(state)-new_max_dens)/(new_min_dens-max_dens);
                    new_vol_frac_h = 1.0 - new_vol_frac_l;
                    vol_frac_l = (Dens(state)-max_dens)/(min_dens-max_dens);
                    vol_frac_h = 1.0 - vol_frac_l;
                    fprintf(outfile10,"%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n", x_ccoords,y_ccoords,z_ccoords,vol_frac_l,vol_frac_h,vol_frac_l*vol_frac_h,new_vol_frac_l,new_vol_frac_h,pdens(state)[0]/Dens(state),pdens(state)[1]/Dens(state),pdens(state)[0],pdens(state)[1], pdens(state)[0]*pdens(state)[1]);
                    }

                }
            }
        }
        fclose(outfile1);
        fclose(outfile2);
        fclose(outfile3);
        fclose(outfile4);
        fclose(outfile5);
        fclose(outfile6);
        fclose(outfile7);
        fclose(outfile8);
        fclose(outfile9);
        fclose(outfile10);

        fprintf(stdout, "max_dens %e min_dens %e new_max_dens %e new_min_dens %e\n",max_dens,min_dens,new_max_dens,new_min_dens);
        }
*/       
}
