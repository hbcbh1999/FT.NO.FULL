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
*			hsoln.c:
*
*    Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    A triangulation of the computational region is constructed
*    from a rectangular grid and an interface.  The corners of
*    the triangles and quadrangles in the triangulation are either
*    intersections between grid lines of the rectangular grid or
*    intersections of grid lines with the interface.  States are
*    associated with these corners, so this triangulation may be
*    used to interpolate states at arbitary points.
*/



#include <hyp/hdecs.h>

	/* LOCAL Function Declarations */
LOCAL	TRI_SOLN *copy_tri_soln_storage(TRI_SOLN*);
LOCAL	boolean  coords_outside_subdomain(double*,register TRI_SOLN*);
LOCAL	boolean  tg_state_from_interface(double*,COMPONENT,register TRI_SOLN*,
					 Locstate,Locstate);
LOCAL	boolean  tg_solution(double*,COMPONENT,register TRI_SOLN*,Locstate,
	                     Locstate);
LOCAL	void  h_hyp_grad_solution(double*,COMPONENT,HYPER_SURF*,SIDE,Front*,
	                            POINTER,Locstate*);    
LOCAL	void  linear_interp_coefs_for_element(double*,double*,
	                                         LINEAR_ELEMENT*,int);
#if defined(TWOD) || defined(THREED)
LOCAL	double integral_over_cell(int*,COMPONENT,TRI_SOLN*,Locstate,POINTER);
LOCAL	void  least_square_coefs_for_cluster(LEAST_SQR_CLUSTER*);
#endif /* defined(TWOD) || defined(THREED) */
#if defined(ONED)
LOCAL	boolean  tg_solution1d(double*,COMPONENT,register TRI_SOLN*,Locstate,
	                       Locstate);
#endif /* defined(ONED) */
#if defined(TWOD)
LOCAL	boolean  fd_grad_solution2d(double*,COMPONENT,register TRI_SOLN*,
	                            Locstate*);
LOCAL	void     flux2d(int*,COMPONENT,TRI_SOLN*,Locstate,Locstate);
#endif /* defined(TWOD) */

LOCAL	double coords_on[3]; /*Coordinates of last hyp_solution evaluation*/
EXPORT	const double *hs_coords_on = coords_on;

/*
*			h_hyp_solution():
*
*    Calculates the state at an arbitrary position coords, as considered
*    to belong to a given component comp.  The state is obtained by
*    interpolation between interior states and states on the front
*    and is loaded in the storage pointed to by state.  If comp is
*    exterior_component(front->interf) then boundary conditions are applied.
*/

/* ARGSUSED */
EXPORT void h_hyp_solution(
	double      *coords,
	COMPONENT  comp,
	HYPER_SURF *hs,
	SIDE        side,
	Front      *front,
	POINTER    wv,
	Locstate   state,
	Locstate   dflt_state)
{
	HYPER_SURF_ELEMENT *hse;
	INTERFACE          *intfc = front->interf;
	double              t[MAXD];
	double              newcrds[MAXD];
	int                i, dim = intfc->dim;
	int                wave_t;
	Wave		   *wave = (Wave*)wv;

	for (i = 0; i < dim; ++i)
	    coords_on[i] = coords[i];
	
	remove_from_debug("line_pj");
	if(debugging("h_hyp_solution"))
	{
	    /*tg_solution */
	    /*    tg_build */
	    double   *p, tol = 1.0e-6;
	    double   p1[3] = {0.00342148271731775,   -0.01892920439794145,     0.1433219887668478 };
	    double   p2[3] = {0.00342148271731775,   -0.01892920439794145,     0.1433219887668478 };

	    if(pp_mynode() == 0)
		p = p1;
	    else
		p = p2;
	    
	    p = p1;

	    if(fabs(coords[0]-p[0]) < tol  &&
	       fabs(coords[1]-p[1]) < tol  &&
	       fabs(coords[2]-p[2]) < tol)
	    {
	        add_to_debug("line_pj");
	        printf("#line_pj begin comp = %d \n", component(coords, intfc));
	    }
	}

	if (!is_excluded_comp(comp,intfc))
	{
	    if (!tri_solution(coords,comp,wave_tri_soln(wave),state,dflt_state))
	    {
	        screen("ERROR in h_hyp_solution(), tri_solution() failed "
	               "in interior region\n");
	        clean_up(ERROR);
	    }
	    return;
	}

	for (i = 0; i < dim; ++i)
	    newcrds[i] = coords[i];

	if (hs == NULL)
	{
	    /*
	    *  Call to nearest_interface_point() is okay
	    *  as interface topology is now fixed
	    */
	    if (!nearest_interface_point(coords,comp,intfc,NO_SUBDOMAIN,
					 NULL,coords_on,t,&hse,&hs))
	    {
	        screen("ERROR in h_hyp_solution(), "
	               "can't find nearest interface point\n");
	        clear_state(front->interf,state,front->sizest);
	        clean_up(ERROR);
	        return;
	    }
	}
	wave_t = wave_type(hs);
	
	switch (wave_t)
	{
	case PASSIVE_BOUNDARY:
	case NEUMANN_BOUNDARY:
	    obstacle_state(front->interf,state,front->sizest);
	    return;

	case DIRICHLET_BOUNDARY:    /* far field conditions */
	    evaluate_dirichlet_boundary_state(coords,hs,front,wave,state);
	    return;

	case SUBDOMAIN_BOUNDARY:

	    (void) printf("WARNING in h_hyp_solution(), attempting to "
	                  "evaluate state on a subdomain boundary\n");
	    comp = (negative_component(hs) == comp) ? positive_component(hs) :
	                                  negative_component(hs);
	    if (!tri_solution(newcrds,comp,wave_tri_soln(wave),state,dflt_state))
	    {
	        screen("ERROR in h_hyp_solution(), "
	               "tri_solution)() failed at SUBDOMAIN_BOUNDARY\n");
	        clean_up(ERROR);
	    }

	    return;
	default:
	    screen("ERROR in h_hyp_solution(), unknown boundary type\n");
	    break;
	}
}            /*end h_hyp_solution*/

EXPORT	void	evaluate_dirichlet_boundary_state(
	double	   *coords,
	HYPER_SURF *hs,
	Front      *front,
	Wave       *wave,
	Locstate   state)
{
	if (boundary_state_function(hs) != NULL)
	{
	    (*boundary_state_function(hs))(coords,hs,front,(POINTER)wave,state);
	}
	else if (boundary_state(hs) != NULL)
	{
	    ft_assign(state,boundary_state(hs),front->sizest);
	}
	else
	{
	    screen("ERROR in evaluate_dirichlet_boundary_state(), "
		   "NULL boundary state\n");
	    clean_up(ERROR);
	}
}		/*end evaluate_dirichlet_boundary_state*/


/*
*			h_hyp_grad_solution():
*
*    Calculates the state gradient at an arbitrary point x,y, as considered
*    to belong to a given component comp.  The gradient is obtained by
*    interpolation between interior states and states on the front
*    and is loaded in the storage pointed to by grad_state.
*    If comp is exterior_component(front->interf) then boundary
*    conditions are applied.
*/

/* ARGSUSED */
LOCAL void h_hyp_grad_solution(
	double      *coords,
	COMPONENT  comp,
	HYPER_SURF *hs,
	SIDE        side,
	Front      *front,
	POINTER    wv,
	Locstate   *grad_state)
{
	HYPER_SURF_ELEMENT *hse;
	INTERFACE          *intfc = front->interf;
	double              t[MAXD];
	double              newcoords[MAXD];
	int                i, dim = intfc->dim;
	int                wave_t;
	Wave		   *wave = (Wave*)wv;

	if (!is_exterior_comp(comp,intfc))
	{
	    if (!grad_tri_solution(coords,comp,wave_tri_soln(wave),grad_state))
	    {
	        screen("ERROR in h_hyp_grad_solution(), "
	               "grad_tri_solution)() failed in interior region\n");
	        clean_up(ERROR);
	    }
	    return;
	}

	for (i = 0; i < dim; ++i)
	    newcoords[i] = coords[i];
	if (hs == NULL)
	{
	    if (nearest_interface_point(coords,comp,intfc,INCLUDE_BOUNDARIES,
					NULL,coords_on,t,&hse,&hs) != YES)
	    {
	        screen("ERROR in h_hyp_grad_solution(), "
	               "nearest_interface_point() failed\n");
	        clean_up(ERROR);
	        return;
	    }
	}
	wave_t = wave_type(hs);
	switch (wave_t)
	{
	case PASSIVE_BOUNDARY:
	case NEUMANN_BOUNDARY:
	case DIRICHLET_BOUNDARY:    /* far field conditions */
	        /* TODO: code needed */
	    screen("ERROR in h_hyp_grad_solution(), "
	           "Unable to return state gradient "
	           "for boundary component at (%g,%g)\n",
	            coords[0],coords[1]);
	    clean_up(ERROR);
	    return;
	case SUBDOMAIN_BOUNDARY:

	    comp = (negative_component(hs) == comp) ?
	            positive_component(hs) :
	            negative_component(hs);
	    if (!grad_tri_solution(newcoords,comp,
			           wave_tri_soln(wave),grad_state))
	    {
	        screen("ERROR in h_hyp_grad_solution(), "
	               "grad_tri_solution)(), "
	               "failed at SUBDOMAIN_BOUNDARY\n");
	        clean_up(ERROR);
	    }
	    return;
	default:
	    screen("ERROR in h_hyp_grad_solution(), unknown boundary type\n");
	    clean_up(ERROR);
	}

}            /*end h_hyp_grad_solution*/


/*
*		init_hyp_solution_function():
*
*    Performs the necessary initialization to enable use of
*    h_hyp_solution_function():  constructs the triangulated
*    grid and allocates the storage for the interior states.
*    Also initializes wave_areas(wave) for use in Rect_area().
*
*    Note that the interior states need to be loaded following
*    the call to this routine in order for h_hyp_solution() to work.
*/


EXPORT int init_hyp_solution_function(
	Wave  *wave,
	Front *front)
{
	RECT_GRID  Dual_grid, *comp_grid = wave->rect_grid;
	int        status = ERROR_IN_STEP;

	clear_wave_pointers(wave);
	if (wave->sizest == 0)
	    return GOOD_STEP;

	/*
	*    hsoln.c locates interior states at crossings
	*    of grid lines; thus it should not be given
	*    wave->rect_grid but rather its dual.
	*/

	set_dual_grid(&Dual_grid,comp_grid);

	scalar(&wave_tri_soln(wave),sizeof(TRI_SOLN));
	if (wave_tri_soln(wave) == NULL)
	{
	    (void) printf("WARNING in init_hyp_solution_function(), "
	                  "can't allocate tri_soln\n");
	    return ERROR_IN_STEP;
	}
#if defined(USE_OVERTURE)
        wave_tri_soln(wave)->cg_over = wave->cg_over;
        wave_tri_soln(wave)->patch_number = wave->patch_number;
        wave_tri_soln(wave)->use_overture_state = wave->use_overture_state;
        wave_tri_soln(wave)->overture_init_step = wave->overture_init_step;
        wave_tri_soln(wave)->cg_over_function = wave->cg_over_function;
        wave_tri_soln(wave)->patch_level = wave->patch_level;
        wave_tri_soln(wave)->NumberOfLevels = wave->NumberOfLevels;
        if (wave->patch_number > 0)
        {
            wave_tri_soln(wave)->patch_component = wave->patch_component;
        }
#endif /* if defined(USE_OVERTURE) */

	wave_tri_soln(wave)->Tri_grid_hooks = wave->Tri_grid_hooks;

	status = hyp_tri_grid_driver(front,wave,&Dual_grid);

	if (status != GOOD_STEP)
	{
	    (void) printf("WARNING in init_hyp_solution_function(), "
	                  "hyp_tri_grid_driver() failed\n");
	    free_wave_pointers(wave);
	    return status;
	}

#if defined(USE_OVERTURE)
        wave->cg_over_function =
                wave_tri_soln(wave)->tri_grid->cg_over_function;
        wave_tri_soln(wave)->cg_over_function =
                wave_tri_soln(wave)->tri_grid->cg_over_function;
#endif /* if defined(USE_OVERTURE) */

	front->_hyp_solution = h_hyp_solution;
	front->_hyp_grad_solution = h_hyp_grad_solution;

	return GOOD_STEP;
}            /*end init_hyp_solution_function*/

EXPORT	void	reinit_hyp_solution_function(
	Wave *wave,
	Front *front)
{
	Wave		*tempwave;

	tempwave = copy_wave(wave);
	clear_wave_pointers(tempwave);
	start_clock("init_hyp_solution");
	if (init_hyp_solution_function(tempwave,front) != GOOD_STEP)
	{
	    screen("ERROR in reinit_hyp_solution_function(), "
	           "init_hyp_solution_function() failed\n");
	    clean_up(ERROR);
	}
	stop_clock("init_hyp_solution");


	ft_assign(wave_tri_soln(tempwave)->tri_grid->rect_state_storage,
		wave_tri_soln(wave)->tri_grid->rect_state_storage,
		wave_tri_soln(wave)->tri_grid->n_reg_nodes*wave->sizest);

	free_wave_pointers(wave);
	assign_wave_pointers(wave,tempwave);
	clear_wave_pointers(tempwave);
	free_wave(tempwave);
}		/*end reinit_hyp_solution_function*/



/*
*			copy_hyp_solution():
*
*    Makes an empty copy of tri_soln/grid storage in nwave
*    from owave.
*/

EXPORT	int copy_hyp_solution_function(
	Wave *owave,
	Wave *nwave)
{
	clear_wave_pointers(nwave);
	if (nwave->sizest == 0)
	    return YES;
#if defined(USE_OVERTURE)
        wave_tri_soln(nwave) = copy_AMR_tri_soln_storage(wave_tri_soln(owave),
                                 nwave);
#else /* if defined(USE_OVERTURE) */	
	wave_tri_soln(nwave) = copy_tri_soln_storage(wave_tri_soln(owave));
#endif /* defined(USE_OVERTURE) */
	if (wave_tri_soln(nwave) == NULL)
	    return NO;

	return YES;
}            /*end copy_hyp_solution_function*/


/*
*			copy_tri_soln_storage():
*
*    This routine makes a copy of the tri_grid storage in soln->tri_grid
*    and loads it in nsoln->tri_grid. It is intended for use by the
*    split stencil schemes for temporary state storage after the first
*    stencil sweep.
*    Only the off-front storage is copied. The front storage is shared
*    between the copy and the original.
*/

LOCAL	TRI_SOLN *copy_tri_soln_storage(
	TRI_SOLN *soln)
{
	TRI_SOLN *nsoln;
	TRI_GRID *grid;

	scalar(&nsoln,sizeof(TRI_SOLN));
	grid = allocate_tri_grid(&soln->tri_grid->tri_grid_hooks);
	if (grid == NULL)
	{
	    free(nsoln);
	    return NULL;
	}

	copy_tri_grid(soln->tri_grid,grid,soln->sizest);

	set_tri_soln_struct(nsoln,soln->intfc,grid,soln->sizest,
	           &soln->interpolator,&soln->el_integral,&soln->unsplit);
	return nsoln;
}        /*end copy_tri_soln_storage*/

/*rect_grid is the dual of computational grid. */
EXPORT	int hyp_tri_grid_driver(
	Front     *front,
	Wave      *wave,
	RECT_GRID *rect_grid)
{
	TRI_GRID           *grid;
	int                status;
	TRI_SOLN           *soln = wave_tri_soln(wave);
	INTERFACE          *intfc = front->interf;
	size_t             sizest = front->sizest;
	INTERPOLATORS      *interpolator = &wave->interpolator;
	EL_INTEGRALS       *el_integral = &wave->el_integral;
	UNSPLIT            *unsplit = &wave->unsplit;

	debug_print("tri_grid","Entered hyp_tri_grid_driver()\n");

#if defined(DEBUG_TRI_GRID)
	if (debugging("tri_grid"))
	{
	    (void) printf("Interface input to tri grid %llu\n",
	              interface_number(intfc));
	    print_interface(intfc);
	    (void) printf("Trigrid rectangular grid\n");
	    print_rectangular_grid(rect_grid);
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_TRI_GRID) */

	grid = allocate_tri_grid(&soln->Tri_grid_hooks);
	if (grid == NULL)
	    return ERROR_IN_STEP;

#if defined(USE_OVERTURE)
        grid->cg_over = soln->cg_over;
        grid->cg_over_function = soln->cg_over_function;
        grid->patch_number = soln->patch_number;
        grid->use_overture_state = soln->use_overture_state;
        grid->overture_init_step = soln->overture_init_step;
        grid->patch_component = soln->patch_component;
        grid->patch_level = soln->patch_level;
        grid->NumberOfLevels = soln->NumberOfLevels;
#endif /* if defined(USE_OVERTURE) */

	if (wave->old_wave != NULL)
	    grid->old_tri_grid = wave_tri_soln(wave->old_wave)->tri_grid;
	else
	    grid->old_tri_grid = NULL;
	status = construct_tri_grid(grid,rect_grid,front);

	if (status != GOOD_STEP)
	{
	    (void) printf("WARNING in hyp_tri_grid_driver(), "
	                  "construct_tri_grid() failed\n");
	}

	set_tri_soln_struct(soln,intfc,grid,sizest,interpolator,
	                    el_integral,unsplit);

	debug_print("tri_grid","Left hyp_tri_grid_driver()\n");
	return status;
}        /*end hyp_tri_grid_driver*/




EXPORT	void set_tri_soln_struct(
	TRI_SOLN      *soln,
	INTERFACE     *intfc,
	TRI_GRID      *grid,
	size_t        sizest,
	INTERPOLATORS *interpolator,
	EL_INTEGRALS  *el_integral,
	UNSPLIT       *unsplit)
{
	soln->intfc    = intfc;
	soln->tri_grid = grid;
	soln->sizest   = sizest;

	soln->interpolator = *interpolator;
	soln->el_integral = *el_integral;
	soln->unsplit = *unsplit;

	switch (intfc->dim)
	{
#if defined(ONED)
	case 1:
	    soln->_tri_solution      = tg_solution1d;
	    soln->_grad_tri_solution = NULL;
	    soln->flux               = NULL;
	    soln->integral           = NULL;
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    soln->_tri_solution      = tg_solution;
	    soln->_grad_tri_solution = fd_grad_solution2d;
	    soln->flux               = flux2d;
	    soln->integral           = integral_over_cell;
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:    /* TODO */
	    soln->_tri_solution      = tg_solution;
	    soln->_grad_tri_solution = NULL;
	    soln->flux               = NULL;
	    soln->integral           = integral_over_cell;
	    break;
#endif /* defined(THREED) */
	}
}        /*end set_tri_soln_struct*/


EXPORT	double    measure_of_linear_element(
	LINEAR_ELEMENT *et,
	TRI_SOLN       *soln)
{
	RECT_GRID *gr = &soln->tri_grid->comp_grid;
	TG_PT     **p = et->p;
	double     measure;
	double     *h = gr->h;
	int       dim = gr->dim;

	switch (dim)
	{
#if defined(ONED)
	case 1:
	    measure = fabs(Coords(p[1])[0] - Coords(p[0])[0])/h[0];
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    measure = 0.5*((Coords(p[1])[0] - Coords(p[0])[0])*
	                   (Coords(p[2])[1] - Coords(p[0])[1]) -
	               (Coords(p[2])[0] - Coords(p[0])[0])*
	                   (Coords(p[1])[1] - Coords(p[0])[1])) / 
	                   (h[0] * h[1]);
	    
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    {
	        double v1[MAXD], v2[MAXD], v3[MAXD];
	        int i;
	        for (i = 0; i < dim; ++i)
	        {
	            v1[i] = Coords(p[1])[i] - Coords(p[0])[i];
	            v2[i] = Coords(p[2])[i] - Coords(p[0])[i];
	            v3[i] = Coords(p[3])[i] - Coords(p[0])[i];
	        }
	        measure = 0.5 * triple_product(v1,v2,v3,dim) /
	                (h[0] * h[1] * h[2]);
	    }
	    break;
#endif /* defined(THREED) */
	}
	return measure;
}        /*end measure_of_linear_element*/


/*
*			tg_solution():
*
*    Returns a value for the Locstate "answer" at the point coords
*    corresponding to component comp.
*
*    Uses either    - interpolation in a rectangular element
*            - interpolation in a triangular element
*            - a call to nearest intfc_state if locate() fails
*              to find a quadrangle or triangle corresponding
*              to coords,comp. (Can happen due to roundoff errors.)
*/

LOCAL	boolean tg_solution(
	double             *coords,
	COMPONENT         comp,
	register TRI_SOLN *soln,
	Locstate          answer,
	Locstate          dflt_answer)
{
	TRI_GRID         *grid = soln->tri_grid;
	BLK_EL0          *blk_el0;
	LINEAR_ELEMENT   *lin;
	BILINEAR_ELEMENT *bilin;
	LEAST_SQR_CLUSTER *lsq;
	int              i, dim = grid->comp_grid.dim;

	if (!Locate_on_trigrid(coords,comp,grid,&bilin,&lin,&lsq))
	{
	    if (dflt_answer != NULL)
	    {
	        ft_assign(answer,dflt_answer,soln->sizest);
	        return FUNCTION_SUCCEEDED;
	    }
	    else if (!tg_state_from_interface(coords,comp,soln,answer,
	    			dflt_answer))
	    {
	        screen("ERROR in tg_solution(), "
	               "tg_state_from_interface() (1) failed\n");
		print_general_vector("coords =",coords,dim,"\n");
		(void) printf("comp = %d, dflt_answer = 0x%p\n",comp,
		        dflt_answer);
		(void) printf("Current interface\n");
		print_interface(current_interface());
	        clean_up(ERROR);
	    }
	    return FUNCTION_SUCCEEDED;
	}

	if(debugging("line_pj"))
	    printf("#line_pj  %d %d %d\n", bilin, lin,  lsq);

	if (bilin != NULL)    /* interpolation in quadrangle */
	{
	    register BILINEAR_ELEMENT *q = bilin;
	    static    int umax[3] = {1, 2, 7};
	    double    f[MAXD];
	    double    *l, *u;

	    l = Coords(q->p[0]); u = Coords(q->p[umax[dim-1]]);
	    for (i = 0; i < dim; ++i)
	    {
	        f[i] = (coords[i] - l[i])/(u[i] - l[i]);
		f[i] = max(f[i],0.0);
		f[i] = min(f[i],1.0);
	    }

	    Bilinear_cell_interpolate(f,q,soln,answer);
	}
	else if (lsq != NULL)
	{
	    least_square_coefs_for_cluster(lsq);
	    if (lsq->nr < dim+1)
	    {
	        if (!tg_state_from_interface(coords,comp,soln,answer,
	                                     dflt_answer))
	        {
	            screen("ERROR in tg_solution(), "
	                   "tg_state_from_interface() (3) failed\n");
	            clean_up(ERROR);
	        }
	        return FUNCTION_SUCCEEDED;
	    }
	    Least_square_interpolate(lsq,soln,coords,answer);
	}
	else 		/* interpolation in triangle */
	{
	    register LINEAR_ELEMENT *t = lin;
	    double    f[MAXD+1];

	    linear_interp_coefs_for_element(f,coords,t,dim);

	    for (i = 0; i <= dim; ++i)
	    {
	        if ((f[i] <= -0.01) || (f[i] >= 1.01))/*TOLERANCE*/
	        {
	            if (!tg_state_from_interface(coords,comp,soln,
				                 answer,dflt_answer))
	            {
	                screen("ERROR in tg_solution(), "
	                       "tg_state_from_interface() (2) failed\n");
	                clean_up(ERROR);
	            }
	            return FUNCTION_SUCCEEDED;
	        }
	    }

	    if (!Linear_cell_interpolate(f,t,soln,answer))
	    {
	        LINEAR_ELEMENT   *et;
	        int      num_els, k;

	        print_general_vector("\nTRI_INTRP fails at ",coords,dim,"");
	        (void) printf("comp %d in element\n\t",comp);
	        for (i = 0; i <= dim; ++i)
	            print_general_vector("",(double *)t->p[i],dim,"");

	        blk_el0 = blk_el0_for_coords(coords,grid);
	        num_els = num_lin_els_in_blk(blk_el0);
	        (void) printf("\n\t\t%d Linear elements in block \n\n",num_els);
	        (void) printf("\n");
	        for (k = 0;  k < num_els;  ++k)
	        {
	            et = blk_el0_linear_els(blk_el0) + k;
	            (void) printf("  element %d:\n",k);
	            print_LINEAR_ELEMENT(et,grid);
	            (void) printf("\n");
	        }

	        if (!tg_state_from_interface(coords,comp,soln,answer,
	                                     dflt_answer))
	        {
	            screen("ERROR in tg_solution(), "
	                   "tg_state_from_interface() (3) failed\n");
	            clean_up(ERROR);
	        }
	        return FUNCTION_SUCCEEDED;
	    }
	}

	return FUNCTION_SUCCEEDED;
}        /*end tg_solution*/

EXPORT	void    states_on_bilinear_element(
	Locstate         *sts,
	BILINEAR_ELEMENT *el,
	TRI_GRID         *tg)
{
	register Locstate *states = tg->states;
	register TG_PT    *nodes  = tg->node_points;
	int               i, nc = (1 << tg->rect_grid.dim);

	for (i = 0; i < nc; ++i)
	    sts[i] = states[el->p[i] - nodes];
}        /*end states_on_bilinear_element*/


/*
*			tg_state_from_interface():
*
*    Attempts to locate a state on the interface when the tri-grid cannot
*    be used, for example when coords is exterior to the domain.
*/

LOCAL	boolean    tg_state_from_interface(
	double             *coords,
	COMPONENT         comp,
	register TRI_SOLN *soln,
	Locstate          answer,
	Locstate          dflt_answer)
{
	HYPER_SURF **hs;
	RECT_GRID  *gr;
	double      *h;
	int        i, dim;

	if ((dflt_answer != NULL) &&
	    (coords_outside_subdomain(coords,soln) == YES))
	{
	    ft_assign(answer,dflt_answer,soln->sizest);
	    return FUNCTION_SUCCEEDED;
	}

	if (nearest_intfc_state(coords,comp,soln->intfc,answer,coords_on,NULL))
	    return FUNCTION_SUCCEEDED;

	for (hs = soln->intfc->hss; hs && *hs; ++hs)
	    if (!is_subdomain_boundary(*hs))
	        return FUNCTION_FAILED;

	gr = computational_grid(soln->intfc);
	h = gr->h;
	dim = gr->dim;

	        /*Only subdomain boundaries*/
	for (i = 0; i < dim; ++i)
	{
	    if (coords[i] < (gr->VL[i]+h[i]))
	        coords_on[i] = gr->VL[i]+h[i];
	    else if (coords[i] > (gr->VU[i]-h[i]))
	        coords_on[i] = gr->VU[i]-h[i];
	    else
	        coords_on[i] = coords[i];
	}
	if (!tg_solution(coords_on,comp,soln,answer,dflt_answer))
	{
	    (void) printf("WARNING in empty_interface_state(), "
	                  "tg_solution() failed\n");
	    return FUNCTION_FAILED;
	}
	return FUNCTION_SUCCEEDED;
}	/*end tg_state_from_interface*/

LOCAL	boolean coords_outside_subdomain(
	double	*coords,
	register TRI_SOLN *soln)
{
	TRI_GRID  *grid = soln->tri_grid;
	RECT_GRID *c_gr = &grid->comp_grid;
	RECT_GRID *gr = &grid->rect_grid;
	int       i, dim = c_gr->dim;

	for (i = 0; i < dim; ++i)
	{
	    if ((coords[i] < gr->L[i]) && (c_gr->lbuf[i] > 0))
		return YES;
	    if ((coords[i] > gr->U[i]) && (c_gr->ubuf[i] > 0))
		return YES;
	}

	return NO;
}	/*end coords_outside_subdomain*/

#if defined(ONED)
/*
*			tg_solution1d():
*
*    Returns a value for the Locstate "answer" at the point coords
*    corresponding to component comp.
*/

LOCAL	boolean tg_solution1d(
	double             *coords,
	COMPONENT         comp,
	register TRI_SOLN *soln,
	Locstate          answer,
	Locstate          dflt_answer)
{
	POINT     **p;
	CRXING    *crx, *crxl, *crxr;
	COMPONENT comp_ic, l_comp, r_comp;
	TRI_GRID  *grid = soln->tri_grid;
	INTERFACE *intfc = grid->grid_intfc;
	Table	  *T = table_of_interface(intfc);
	RECT_GRID *gr = &grid->comp_grid;
	int       ic[MAXD], icp1[MAXD], icm1[MAXD], icgr[MAXD];
	int       i;
	int       *list, nc;
	double     x_ic, xl, x, xr;
	double     alpha, beta;
	Locstate  stl, str;

	if ((!rect_in_which(coords,ic,gr)) ||
	    (!rect_in_which(coords,icgr,&grid->rect_grid)))
	{
	    if (dflt_answer != NULL)
	    {
		ft_assign(answer,dflt_answer,soln->sizest);
	        return FUNCTION_SUCCEEDED;
	    }
	    else
	    {
	        screen("ERROR in tg_solution1d(), rect_in_which() failed\n");
		print_general_vector("coords = ",coords,intfc->dim,"\n");
		print_interface(intfc);
	        clean_up(ERROR);
	        return FUNCTION_FAILED;
	    }
	}
	comp_ic = Regular_grid_comp(ic,grid);
	nc = T->seg_crx_count[icgr[0]];
	x = coords[0];
	x_ic = Coords(Regular_grid_node(ic,grid))[0];
	icm1[0] = ic[0] - 1;
	icp1[0] = ic[0] + 1;
	if ((comp_ic == comp) && (nc == 0))
	{
	    stl = str = NULL;
	    if (x < x_ic)
	    {
		if (Regular_grid_comp(icm1,grid) == comp)
		{
	            xl = Coords(Regular_grid_node(icm1,grid))[0];
	            stl = Regular_grid_state(icm1,grid);
		}
	        xr = x_ic;
	        str = Regular_grid_state(ic,grid);
	    }
	    else
	    {
	        xl = x_ic;
	        stl = Regular_grid_state(ic,grid);
		if (Regular_grid_comp(icp1,grid) == comp)
		{
	            xr = Coords(Regular_grid_node(icp1,grid))[0];
	            str = Regular_grid_state(icp1,grid);
		}
	    }
	    if ((stl != NULL) && (str != NULL))
	    {
	        alpha = (xr - x)/(xr - xl);
	        beta = 1.0 - alpha;
	        bi_interpolate_intfc_states(intfc,alpha,beta,
	                        &xl,stl,&xr,str,answer);
	        return FUNCTION_SUCCEEDED;
	    }
	}

	if (x < x_ic)
	{
	    l_comp = Regular_grid_comp(icm1,grid);
	    r_comp = comp_ic;
	}
	else
	{
	    l_comp = comp_ic;
	    r_comp = Regular_grid_comp(icp1,grid);
	}
	/*
	*  Locate points on the left and right of coords with the
	*  correct component
	*/
	list = T->seg_crx_lists[icgr[0]];
	crxl = NULL;
	for (i = nc-1; i >= 0; i--)
	{
	    crx = &(T->crx_store[list[i]]);
	    if ((coords[0] >= Coords(crx->pt)[0]) &&
	        (comp == positive_component(crx->pt)))
	    {
	        crxl = crx;
	        break;
	    }
	}
	crxr = NULL;
	for (i = 0; i < nc; ++i)
	{
	    crx = &(T->crx_store[list[i]]);
	    if ((coords[0] <= Coords(crx->pt)[0]) &&
	        (comp == negative_component(crx->pt)))
	    {
	        crxr = crx;
	        break;
	    }
	}
	if ((crxl == NULL) && (comp == l_comp) && (crxr != NULL))
	{
	    if (x < x_ic)
	    {
	        xl = Coords(Regular_grid_node(icm1,grid))[0];
	        stl = Regular_grid_state(icm1,grid);
	    }
	    else
	    {
	        xl = Coords(Regular_grid_node(ic,grid))[0];
	        stl = Regular_grid_state(ic,grid);
	    }
	    xr = Coords(crxr->pt)[0];
	    str = left_state(crxr->pt);
	    alpha = (xr - x)/(xr - xl);
	    beta = 1.0 - alpha;
	    bi_interpolate_intfc_states(intfc,alpha,beta,
	                    &xl,stl,&xr,str,answer);
	    return FUNCTION_SUCCEEDED;
	}
	if ((crxr == NULL) && (comp == r_comp) && (crxl != NULL))
	{
	    xl = Coords(crxl->pt)[0];
	    stl = right_state(crxl->pt);
	    if (x < x_ic)
	    {
	        xr = Coords(Regular_grid_node(ic,grid))[0];
	        str = Regular_grid_state(ic,grid);
	    }
	    else
	    {
	        xr = Coords(Regular_grid_node(icp1,grid))[0];
	        str = Regular_grid_state(icp1,grid);
	    }
	    alpha = (xr - x)/(xr - xl);
	    beta = 1.0 - alpha;
	    bi_interpolate_intfc_states(intfc,alpha,beta,
	                    &xl,stl,&xr,str,answer);
	    return FUNCTION_SUCCEEDED;
	}
	if ((crxl != NULL) && (crxr != NULL))
	{
	    xl = Coords(crxl->pt)[0];
	    stl = right_state(crxl->pt);
	    xr = Coords(crxr->pt)[0];
	    str = left_state(crxr->pt);
	    alpha = (xr - x)/(xr - xl);
	    beta = 1.0 - alpha;
	    bi_interpolate_intfc_states(intfc,alpha,beta,
	                    &xl,stl,&xr,str,answer);
	    return FUNCTION_SUCCEEDED;
	}
	p = intfc->points;
	for (i = 0; i < intfc->num_points-1; ++i)
	{
	    if ((Coords(p[i])[0] <= x) && (x <= Coords(p[i+1])[0]) &&
	        (positive_component(p[i]) == comp) &&
		(negative_component(p[i+1]) == comp))
	    {
	        if (x < x_ic)
	        {
	            xl = Coords(Regular_grid_node(icm1,grid))[0];
	            if ((l_comp == comp) && (Coords(p[i])[0] < xl))
	                stl = Regular_grid_state(icm1,grid);
	            else
	            {
	                xl = Coords(p[i])[0];
	                stl = right_state(p[i]);
	            }
	            xr = Coords(Regular_grid_node(ic,grid))[0];
	            if ((r_comp == comp) && (xr < Coords(p[i])[0]))
	                str = Regular_grid_state(ic,grid);
	            else
	            {
	                xr = Coords(p[i+1])[0];
	                str = left_state(p[i+1]);
	            }
	        }
	        else
	        {
	            xl = Coords(Regular_grid_node(ic,grid))[0];
	            if ((l_comp == comp) && (Coords(p[i])[0] < xl))
	                stl = Regular_grid_state(ic,grid);
	            else
	            {
	                xl = Coords(p[i])[0];
	                stl = right_state(p[i]);
	            }
	            xr = Coords(Regular_grid_node(icp1,grid))[0];
	            if ((r_comp == comp) && (xr < Coords(p[i])[0]))
	                str = Regular_grid_state(ic,grid);
	            else
	            {
	                xr = Coords(p[i+1])[0];
	                str = left_state(p[i+1]);
	            }
	        }
	        alpha = (xr - x)/(xr - xl);
	        beta = 1.0 - alpha;
	        bi_interpolate_intfc_states(intfc,alpha,beta,
	                                    &xl,stl,&xr,str,answer);
	        return FUNCTION_SUCCEEDED;
	    }
	}
	if (!nearest_intfc_state(coords,comp,soln->intfc,answer,coords_on,NULL))
	{
	    screen("ERROR in tg_solution1d(), "
	           "nearest_intfc_state() (1) failed\n");
	    clean_up(ERROR);
	    return FUNCTION_FAILED;
	}
	return FUNCTION_SUCCEEDED;
}        /*end tg_solution1d*/
#endif /* defined(ONED) */


LOCAL	void    least_square_coefs_for_cluster(
	LEAST_SQR_CLUSTER *cluster)
{
	TG_PT **pts = cluster->p;
	int p,n;
	int i,j;	
	static	double *x;

	if (x == NULL)
	{
	    stat_vector(&x,4*MAXD*MAX_LSQ_PTS,FLOAT);
	}
	n = cluster->nr;
	switch (cluster->dim)
	{
	case 2:
	    cluster->i_order = 1;
	    cluster->nc = p = 3;
	    /* For future test
	    if (n < 6) 
	    {
	    	cluster->i_order = 1;
		cluster->nc = p = 3;
	    }
	    else 
	    {
	    	cluster->i_order = 2;
		cluster->nc = p = 6;
	    }
	    */
	    break;
	case 3:
	    cluster->i_order = 1;
	    cluster->nc = p = 4;
	    /* For future test
	    if (n < 11) 
	    {
	    	cluster->i_order = 1;
		cluster->nc = p = 4;
	    }
	    else 
	    {
	    	cluster->i_order = 2;
		cluster->nc = p = 11;
	    }
	    */
	}

	for (i = 0; i < n; ++i)
	{
	    for (j = 0; j < p-1; ++j)
	    {
                x[i*p+j] = Coords(pts[i])[j];
            }
            x[i*p+j] = 1.0;
	}
	if (debugging("least_square"))
	{
	    printf("n = %d  p = %d\n",n,p);
	    for (i = 0; i < n; ++i)
	    {
	        for (j = 0; j < p-1; ++j)
		    printf("x[%d][%d] = %f ",i,j,x[j*n+i]);
	        printf("x[%d][%d] = %f\n",i,j,x[j*n+i]);
	    }
	}
	cluster->x = x;
}	/* end least_square_coefs_for_cluster */


LOCAL	void    linear_interp_coefs_for_element(
	double          *f,
	double          *crds,
	LINEAR_ELEMENT *t,
	int            dim)
{
	/*
	*    The below computation can be speeded up
	*    by storing the computation of the linear
	*    element center, side lengths (areas in 3d)
	*    and area (volume in 3d).
	*/

	switch (dim)
	{
#if defined(ONED)
	case 1:
	    f[0] = (crds[0] - Coords(t->p[0])[0]) /
	        (Coords(t->p[1])[0] - Coords(t->p[0])[0]);
	    f[0] = max(0.0,f[0]);	f[0] = min(1.0,f[0]);
	    return;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	{
	    double x0, y0, x1, y1, x2, y2;
	    double xx, yy;
	    double den;

	    x0 = Coords(t->p[0])[0];    y0 = Coords(t->p[0])[1];
	    x1 = Coords(t->p[1])[0] - x0;    y1 = Coords(t->p[1])[1] - y0;
	    x2 = Coords(t->p[2])[0] - x0;    y2 = Coords(t->p[2])[1] - y0;
	    xx = crds[0] - x0;        yy = crds[1] - y0;
	    den = x1*y2 - y1*x2;
	    f[1] = (xx*y2 - yy*x2) / den;
	    f[2] = (x1*yy - y1*xx) / den;
	    f[0] = 1.0 - f[1] - f[2];
	    f[0] = max(0.0,f[0]);	f[0] = min(1.0,f[0]);
	    f[1] = max(0.0,f[1]);	f[1] = min(1.0,f[1]);
	    f[2] = max(0.0,f[2]);	f[2] = min(1.0,f[2]);
	}
	    return;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    /*
	    *     This algorithm assumes the vertices on
	    *    the tetrahedra are numbered such that the
	    *    triple product of (p1-p0), (p2-p0), (p3-p0)
	    *    is positive.
	    *
	    *               p1
	    *               /|\
	    *              / | \
	    *             /  |  \
	    *            /   |   \
	    *               /    |    p3
	    *              /     |   /
	    *             /      |  /
	    *            /       | /
	    *           p0-------p2
	    */

	{
	    double v10, v11, v12;
	    double v20, v21, v22;
	    double v30, v31, v32;
	    double q0, q1, q2;
	    double *p0, *p1, *p2, *p3;
	    double den;

	    p0 = Coords(t->p[0]);    p2 = Coords(t->p[2]);
	    p1 = Coords(t->p[1]);    p3 = Coords(t->p[3]);
	    q0 = crds[0] - p0[0]; q1 = crds[1] - p0[1]; q2 = crds[2] - p0[2];
	    v10 = p1[0] - p0[0]; v11 = p1[1] - p0[1]; v12 = p1[2] - p0[2];
	    v20 = p2[0] - p0[0]; v21 = p2[1] - p0[1]; v22 = p2[2] - p0[2];
	    v30 = p3[0] - p0[0]; v31 = p3[1] - p0[1]; v32 = p3[2] - p0[2];
	    den = QDet3d(v1,v2,v3);
	    if (fabs(den) < MACH_EPS)
	    {
	        f[0] = 0.0;
	        f[1] = 0.0;
	        f[2] = 0.0;
	        f[3] = 1.0;
	        return;
	    }

	    f[1] = QDet3d(q,v2,v3)/den;
	    f[2] = QDet3d(v1,q,v3)/den;
	    f[3] = QDet3d(v1,v2,q)/den;
	    f[0] = 1.0 - f[1] - f[2] - f[3];
	    f[0] = max(0.0,f[0]);	f[0] = min(1.0,f[0]);
	    f[1] = max(0.0,f[1]);	f[1] = min(1.0,f[1]);
	    f[2] = max(0.0,f[2]);	f[2] = min(1.0,f[2]);
	    f[3] = max(0.0,f[3]);	f[3] = min(1.0,f[3]);
	}
	    return;
#endif /* defined(THREED) */
	}
}        /*end linear_interp_coefs_for_element*/


#if defined(TWOD) || defined(THREED)
/*
*			integral_over_cell():
*
*    Returns the integral of the Locstate over the portion of the
*    mesh block icoords corresponding to component "comp", normalized
*    by the area of the mesh block.
*    Also returns the area of the mesh block corresponding to "comp".
*
*    Uses the formula for the integral of a linear function
*    over a linear element (triangle or tetrahedron)
*
*        I = (area or volume) * average
*
*    and for a bilinear function over a rectangular parallel piped
*
*        I = (area or volume) * average
*    
*/

LOCAL	double integral_over_cell(
	int       *icoords,
	COMPONENT comp,
	TRI_SOLN  *soln,
	Locstate  answer,
	POINTER   user)
{
	BLK_EL0  *blk_el0;
	TRI_GRID *grid = soln->tri_grid;
	double    area;

	clear_state(grid->grid_intfc,answer,soln->sizest);

	blk_el0 = &Regular_blk_el0(icoords,grid);

	if (blk_el0_is_bilinear(blk_el0))
	{
	    if (blk_el0_bilinear_el(blk_el0)->comp != comp)
	        area = 0.0;
	    else
	    {
	        register BILINEAR_ELEMENT   *eq;

	        eq = blk_el0_bilinear_el(blk_el0);

	        area = Bilinear_cell_integrate(eq,soln,answer,user);
	    }
	}
	else
	{
	    register LINEAR_ELEMENT   *et;
	    int        num_els, k;

	    num_els = num_lin_els_in_blk(blk_el0);
	    area = 0.0;
	    for (k = 0;  k < num_els;  ++k)
	    {
	        et = blk_el0_linear_els(blk_el0) + k;

	        if (et->comp != comp)
		    continue;

	        area += Linear_cell_integrate(et,soln,answer,user);
	    }
	}
	return area;
}        /*end integral_over_cell*/
#endif /* defined(TWOD) || defined(THREED) */


#if defined(TWOD)

/*
*			tri_grad_solution():
*
*    Returns x and y gradients of Locstate corresponding to x,y,comp.
*
*    Uses either    - interpolation in a quadrangle
*            - interpolation in a triangle
*            - a call to nearest intfc_state if locate() fails
*              to find a quadrangle or triangle corresponding
*              to x,y,comp. (Can happen due to roundoff errors.)
*/


#if defined(UNUSED_FUNCTION)
LOCAL	boolean    tri_grad_solution(double*,COMPONENT,register TRI_SOLN*,
	              Locstate*);
LOCAL	int	nearest_tri(double*,COMPONENT,TRI_GRID*,LINEAR_ELEMENT**);
LOCAL	void	shortest_dist(double*,double*,double*,double*,double*,double*);
LOCAL	void	test_block_for_tri(int*,double*,COMPONENT,TRI_GRID*,
				   LINEAR_ELEMENT**,double*);

LOCAL	boolean tri_grad_solution(
	double             *coords,
	COMPONENT         comp,
	register TRI_SOLN *soln,
	Locstate          *grad_answer)
{
	BILINEAR_ELEMENT *bilin;
	LINEAR_ELEMENT   *lin;
	LINEAR_ELEMENT   Tmp;
	TRI_GRID         *grid = soln->tri_grid;
	double            coords[MAXD];
	int              deb_tgs = debugging("tgs");/* debugging switch */
	static double     **df = NULL;

	if (df == NULL)
	    bi_array(&df,MAXD,MAXD,FLOAT);

	if (deb_tgs)
	    (void) printf("tri_grad_soln for %g %g comp %d\n",
	              coords[0],coords[1],comp);
	if (!Locate_on_trigrid(coords,comp,grid,&bilin,&lin))
	{
	    if (deb_tgs)
	        (void) printf("locate failed\n");
	    if (!nearest_tri(coords,comp,grid,&lin))
	    {
	        HYPER_SURF         *hs;
	        HYPER_SURF_ELEMENT *hse;
	        BOND               *b;
	        CURVE              *c;
	        Locstate           s1, s2;
	        double              tt[MAXD];

	        if (deb_tgs)
	            (void) printf("nearest_tri failed\n");
	        if (nearest_interface_point(coords,comp,soln->intfc,
					    INCLUDE_BOUNDARIES,NULL,
					    coords_on,tt,&hse,&hs) != YES)
	        {
	            int dim = grid->rect_grid.dim;
	            (void) printf("WARNING in tri_grad_solution(), "
	                          "nearest_interface_point() failed\n");
	            print_general_vector("locate failed at ",coords,dim,"");
	            (void) printf("comp %d\n",comp);
#if defined(DEBUG_TRI_GRID)
	            if (debugging("tri_grid"))
			print_blk_els0(grid);
#endif /* defined(DEBUG_TRI_GRID) */

	            (void) printf("Returning zero gradiant\n");
	            zero_scalar(grad_answer[0],size_of_state(soln->intfc));
	            zero_scalar(grad_answer[1],size_of_state(soln->intfc));
	            return FUNCTION_FAILED;
	        }
	        b = Bond_of_hse(hse);
	        c = Curve_of_hs(hs);

	        if (comp == negative_component(c))
	        {
	            s1 = left_state_at_point_on_curve(b->start,b,c);
	            s2 = left_state_at_point_on_curve(b->end,b,c);
	        }
	        else
	        {
	            s1 = right_state_at_point_on_curve(b->start,b,c);
	            s2 = right_state_at_point_on_curve(b->end,b,c);
	        }
	        if (soln->interpolator.grad_bond)
	        {
	            Grad_bond_interpolate(soln,tt,b,s1,s2,grad_answer);
	        }
	        else
	        {
	             /* When all else fails return grad = 0 */

	            Tmp.p[0] = (TG_PT *) b->start;
	            Tmp.p[1] = (TG_PT *) b->end;
	            Tmp.p[2] = (TG_PT *) b->end;
	            Tmp.s[0] = s1;    Tmp.s[1] = s1;    Tmp.s[2] = s1;
	            Tmp.side[0] = F_SIDE;
	            Tmp.side[1] = F_SIDE;
	            Tmp.side[2] = F_SIDE;
	            Tmp.comp = comp;
	            df[0][0] = 0.0    df[0][1] = 0.0;
	            df[1][0] = 0.0    df[1][1] = 0.0;
	            Grad_linear_cell_interpolate(df,&Tmp,soln,grad_answer);
	        }
	        return FUNCTION_FAILED;
	    }
	}

	if (bilin != NULL)         /* interpolation in quadrangle */
	{
	    register BILINEAR_ELEMENT *q = bilin;
	    double     x1, y1, f[MAXD], d[MAXD];

	    x1 = Coords(q->p[0])[0];    y1 = Coords(q->p[0])[1];
	    d[0] = Coords(q->p[2])[0] - x1;    d[0] = Coords(q->p[3])[1] - y1;
	    f[0] = (coords[0] - x1) / d[0];    f[1] = (coords[1] - y1) / d[1];

	    if (deb_tgs)
	        (void) printf("found in rectangular cell\n");

	    Grad_bilinear_cell_interpolate(f,d,q,soln,grad_answer);
	}
	else                /* interpolation in triangle */
	{
	    double    x0, y0, x1, y1, x2, y2, det;
	    double    xx, yy, f[MAXD+1];
	    register LINEAR_ELEMENT *t = lin;

	    x0 = Coords(t->p[0])[0];        y0 = Coords(t->p[0])[1];
	    x1 = Coords(t->p[1])[0] - x0;   y1 = Coords(t->p[1])[1] - y0;
	    x2 = Coords(t->p[2])[0] - x0;   y2 = Coords(t->p[2])[1] - y0;
	    
	    xx = coords[0] - x0;            yy = coords[1] - y0;

	    det = x1*y2 - y1*x2;
	    f[2] = (xx*y2 - yy*x2) / det;    f[3] = (x1*yy - y1*xx) / det;

	    df[0][0] =   y2 / det;        df[0][1] = - x2 / det;
	    df[1][0] = - y1 / det;        df[1][1] =   x1 / det;

	    f[0] = 1.0 - f[2] - f[3];

	    if (deb_tgs)
	        (void) printf("found in lin\n");

	    if ((f[0] <= -0.01) || (f[2] <= -0.005) || (f[3] <= -0.005) ||
	        (f[2] >= 1.01)  || (f[3] >=  1.01))
	    {
	        HYPER_SURF         *hs;
	        HYPER_SURF_ELEMENT *hse;
	        BOND               *b;
	        CURVE              *c;
	        Locstate           s1, s2;
	        double              tt[MAXD];

	        if (deb_tgs)
	            (void) printf("require nearest intfc pt\n");
	        if (nearest_interface_point(coords,comp, soln->intfc,
					    INCLUDE_BOUNDARIES,NULL,
					    coords_on,tt,&hse,&hs) != YES)
	        {
	            (void) printf("WARNING in tri_grad_solution(), "
	                          "nearest_interface_point failed\n"
	                          "Returning zero gradiant\n");
	            zero_scalar(grad_answer[0],size_of_state(soln->intfc));
	            zero_scalar(grad_answer[1],size_of_state(soln->intfc));
	            return FUNCTION_FAILED;
	        }
	        b = Bond_of_hse(hse);
	        c = Curve_of_hs(hs);

	        if (comp == negative_component(c))
	        {
	            s1 = left_state_at_point_on_curve(b->start,b,c);
	            s2 = left_state_at_point_on_curve(b->end,b,c);
	        }
	        else
	        {
	            s1 = right_state_at_point_on_curve(b->start,b,c);
	            s2 = right_state_at_point_on_curve(b->end,b,c);
	        }
	        if (soln->interpolator.grad_bond)
	        {
	            Grad_bond_interpolate(soln,tt,b,s1,s2,grad_answer);
	        }
	        else
	        {
	            /* When all else fails return grad = 0 */
	            Tmp.p[0] = (TG_PT *) b->start;
	            Tmp.p[1] = (TG_PT *) b->end;
	            Tmp.p[2] = (TG_PT *) b->end;
	            Tmp.s[0] = s1;    Tmp.s[1] = s1;    Tmp.s[2] = s1;
	            Tmp.side[0] = F_SIDE;
	            Tmp.side[1] = F_SIDE;
	            Tmp.side[2] = F_SIDE;
	            Tmp.comp = comp;
	            df[0][0] = 0.0    df[0][1] = 0.0;
	            df[1][0] = 0.0    df[1][1] = 0.0;
	            Grad_linear_cell_interpolate(df,&Tmp,soln,grad_answer);
	        }
	        return FUNCTION_FAILED;
	    }

	    if (deb_tgs)
	    {
	        (void) printf("calling grad tri interpolator f %g %g %g\n",
	                        f[0],f[2],f[3]);
	    }
	    if (Grad_linear_cell_interpolate(df,t,soln,grad_answer) == NO)
	    {
	        LINEAR_ELEMENT *et;
	        BLK_EL0        *blk_el0;
	        int            num_els, k;
	        BOND           *b;
	        CURVE          *c;
	        Locstate       s1, s2;
	        double          tt[MAXD];


	        (void) printf("\nGRAD_TRI_INTRP fails at %g %g comp %d ",
	              x,y,comp);
	        (void) printf("in lin\n\t%g %g   %g %g   %g %g\n",
	              Coords(t->p[0])[0],Coords(t->p[0])[1],
	              Coords(t->p[1])[0],Coords(t->p[1])[1],
	              Coords(t->p[2])[0],Coords(t->p[2])[1]);

	        blk_el0 = blk_el0_for_coords(coords,grid);
	        num_els = num_lin_els_in_blk(blk_el0);
	        (void) printf("\n\t\t%d Triangles in block \n\n",num_els);
	        for (k = 0;  k < num_els;  ++k)
	        {
	            et = blk_el0_linear_els(blk_el0) + k;
	            (void) printf("  triangle %d:\n",k);
	            print_LINEAR_ELEMENT(et,grid);
	            (void) printf("\n");
	        }

	        if (nearest_interface_point(coords,comp,soln->intfc,
					    INCLUDE_BOUNDARIES,NULL,
					    coords_on,tt,&b,&c) != YES)
	        {
	            (void) printf("WARNING in tri_grad_solution(), "
	                          "nearest_interface_point failed\n"
	                          "Returning zero gradiant\n");
	            zero_scalar(grad_answer[0],size_of_state(soln->intfc));
	            zero_scalar(grad_answer[1],size_of_state(soln->intfc));
	            return FUNCTION_FAILED;
	        }

	        if (comp == negative_component(c))
	        {
	            s1 = left_state_at_point_on_curve(b->start,b,c);
	            s2 = left_state_at_point_on_curve(b->end,b,c);
	        }
	        else
	        {
	            s1 = right_state_at_point_on_curve(b->start,b,c);
	            s2 = right_state_at_point_on_curve(b->end,b,c);
	        }
	        if (soln->interpolator.grad_bond)
	        {
	            Grad_bond_interpolate(soln,tt,b,s1,s2,grad_answer);
	        }
	        else
	        {
	            /* When all else fails return grad = 0 */
	            Tmp.p[0] = (TG_PT *) b->start;
	            Tmp.p[1] = (TG_PT *) b->end;
	            Tmp.p[2] = (TG_PT *) b->end;
	            Tmp.s[0] = s1;    Tmp.s[1] = s1;    Tmp.s[2] = s1;
	            Tmp.side[0] = F_SIDE;
	            Tmp.side[1] = F_SIDE;
	            Tmp.side[2] = F_SIDE;
	            Tmp.comp = comp;
	            df[0][0] = 0.0    df[0][1] = 0.0;
	            df[1][0] = 0.0    df[1][1] = 0.0;
	            Grad_linear_cell_interpolate(df,&Tmp,soln,
	                         grad_answer);
	        }
	        return FUNCTION_FAILED;
	    }
	}

	return FUNCTION_SUCCEEDED;
}        /*end tri_grad_solution*/


/*
*			nearest_tri():
*
*	Given x, y, comp, this routine returns the nearest triangle
*	lying in one-of-a-possible-5 dual lattice mesh blocks that have
*	component value  comp .  The order of searching mesh blocks is
*				ix,iy
*		ix-1,iy    ix+1,iy     ix  ,iy-1  ix  ,iy+1
*		ix-1,iy-1  ix+1,iy-1   ix-1,iy+1  ix+1,iy+1
*	Since a curve can enter and leave a mesh block by the same side
*	it is possible that a desired triangle is in one of the last four
*	possibilties. It only searches the last four if no triangle is found
*	in ix,iy.
*	This routine is analogous to tg_locate, but does not require x,y
*	to lie in the triangle.
*/

LOCAL	int nearest_tri(
	double		*coords,
	COMPONENT	comp,
	TRI_GRID	*grid,
	LINEAR_ELEMENT	**tri)
{
	INTERFACE	*intfc = grid->grid_intfc;
	RECT_GRID	*gr = &grid->rect_grid;
	LINEAR_ELEMENT  *btri;
	double		min_dist, sav_min;
	int		ix, iy, icoords[MAXD];

#if defined(DEBUG_TRI_LOC)
	if (debugging("near_tri"))
	    (void) printf("searching for nearest tri to %g %g comp %d\n",
		          coords[0],coords[1],comp);
#endif /* defined(DEBUG_TRI_LOC) */

	if (is_exterior_comp(comp,intfc) || (comp < min_component(intfc)) ||
	    				    (comp > max_component(intfc)))
	{
#if defined(DEBUG_TRI_LOC)
	    if (debugging("near_tri"))
	    	(void) printf("comp exterior or < min or > max\n");
#endif /* defined(DEBUG_TRI_LOC) */
	    return NO;
	}

	if (!rect_in_which(coords,icoords,gr))
	{
#if defined(DEBUG_TRI_LOC)
	    if (debugging("near_tri"))
	    	(void) printf("%g %g not in any mesh block\n",
	    		      coords[0],coords[1]);
#endif /* defined(DEBUG_TRI_LOC) */
	    return NO;
	}

	ix = icoords[0];	iy = icoords[1];
	test_block_for_tri(icoords,coords,comp,grid,tri,&min_dist);
	if (*tri == NULL)
	{
	    sav_min = 999999.0;
	    if (ix-1 >= 0)
	    {
		icoords[0] = ix-1;	icoords[1] = iy;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	    if (ix+1 < gr->gmax[0])
	    {
		icoords[0] = ix+1;	icoords[1] = iy;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	    if (iy-1 >= 0)
	    {
		icoords[0] = ix;	icoords[1] = iy-1;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	    if (iy+1 < gr->gmax[1])
	    {
		icoords[0] = ix;	icoords[1] = iy+1;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	    if ((ix-1 >= 0) && (iy-1 >= 0))
	    {
		icoords[0] = ix-1;	icoords[1] = iy-1;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	    if ((ix+1 < gr->gmax[0]) && (iy-1 >= 0))
	    {
		icoords[0] = ix+1;	icoords[1] = iy-1;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	    if ((ix-1 >= 0) && (iy+1 < gr->gmax[1]))
	    {
		icoords[0] = ix-1;	icoords[1] = iy+1;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	    if ((ix+1 < gr->gmax[0]) && (iy+1 < gr->gmax[1]))
	    {
		icoords[0] = ix+1;	icoords[1] = iy+1;
		test_block_for_tri(icoords,coords,comp,grid,&btri,&min_dist);
		if (btri && (min_dist < sav_min))
		{
		    sav_min = min_dist;
		    *tri = btri;
		}
	    }
	}

#if defined(DEBUG_TRI_LOC)
	if (debug_tri_loc)
	{
	    (void) printf("coords %g %g comp %d in mesh block %d %d\n",
	    	          coords[0],coords[1],comp,ix,iy);
	    (void) printf("nearest tri:  ");
	    print_LINEAR_ELEMENT(*tri,grid);
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_TRI_LOC) */
	return ( *tri == NULL) ? NO : YES;
}		/*end nearest_tri*/

#if defined(DEBUG_TRI_LOC)
#define	debug_test_nearby_tri(coords,p0,p1)				\
	if (debugging("near_tri"))					\
	{								\
	    (void) printf("test tri edge %g %g -> %g %g\n",		\
			  p0[0],p0[1],p1[0],p1[1]);			\
	    (void) printf("dist %g norm_dist %g t %g pres_min_dist %g\n",\
			  dist,norm_dist,t,*min_dist);			\
	    (void) printf("new min_dist %g tri %d\n",*min_dist,*tri);	\
	}
#else /* defined(DEBUG_TRI_LOC) */
#define	debug_test_nearby_tri(coords,p0,p1)
#endif /* defined(DEBUG_TRI_LOC) */

#define Test_nearby_tri(coords,p0,p1)					\
	shortest_dist(coords,p0,p1,&dist,&norm_dist,&t);		\
	if ((0.0 < t) && (t < 1.0))					\
	{								\
	    if (norm_dist < *min_dist)					\
		{ *min_dist = norm_dist; *tri = et; }			\
	}								\
	else								\
	{								\
	    if (dist < *min_dist)					\
		{ *min_dist = dist; 	*tri = et; }			\
	}								\
	debug_test_nearby_tri(coords,p0,p1)


LOCAL	void test_block_for_tri(
	int		*icoords,
	double		*coords,
	COMPONENT	comp,
	TRI_GRID	*grid,
	LINEAR_ELEMENT	**tri,
	double		*min_dist)
{
	BLK_EL0		*blk_el0;
	LINEAR_ELEMENT  *et;
	double		*p0, *p1, *p2;
	double		dist, norm_dist, t;
	int		num_els, k;

	*tri = NULL;	*min_dist = 1000000.0;

	blk_el0 = &Blk_el0(icoords,grid);

#if defined(DEBUG_TRI_LOC)
	if (debugging("near_tri"))
	    (void) printf("\nTesting %g %g in block %d %d num_els %d\n",
			  coords[0],coords[1],icoords[0],icoords[1],
	    	          num_lin_els_in_blk(num_els));
#endif /* defined(DEBUG_TRI_LOC) */

	if (blk_el0_is_bilinear(blk_el0))
	    return;

	num_els = num_lin_els_in_blk(blk_el0);
	for (k = 0;  k < num_els;  ++k)
	{
	    et = blk_el0_linear_els(blk_el0) + k;
	    if (et->comp != comp)
	        continue;
	    p0 = Coords(et->p[0]);
	    p1 = Coords(et->p[1]);
	    p2 = Coords(et->p[2]);
	    Test_nearby_tri(coords,p0,p1);
	    Test_nearby_tri(coords,p1,p2);
	    Test_nearby_tri(coords,p2,p0);
	}
}		/*end test_block_for_tri*/

/*
*			shortest_dist():
*
*	Computes the shortest distance (squared) from the point p0
*	to the line segment determined by ps -> pe.
*	In case the line_segment has zero length, distance is set to 1000000
*
*	Also computes the normal distance (squared) to the infinite
*	line through the bond.
*
*	The paramater  t  below is used to find the parametric 
*	position of the closest point on the line segment to x,y.
*/

#define  dist(t)   sqr((t)*x10 - x_0) + sqr((t)*y10 - y_0)
#define  EPS3 	   (EPSILON*EPSILON*EPSILON)		/* from int.h */


LOCAL	void shortest_dist(
	double		*p0,
	double		*ps,
	double		*pe,
	double		*distance,
	double		*norm_dist,
	double		*t)
{
	double		x   = (double)p0[0],	y = (double)p0[1];
	double		x0  = (double)ps[0],	y0  = (double)ps[1];
	double		x1  = (double)pe[0],	y1  = (double)pe[1];
	double		y10 = y1 - y0,		x10 = x1 - x0;
	double		x_0 = x - x0,		y_0 = y - y0;
	double		scalar_prod;
	double		l;		/* squared length of bond */

	scalar_prod = x10 * x_0 + y10 * y_0;
	l = x10 * x10 + y10 * y10;

	if (l <= EPS3)
	{
	    if (l == 0.0)
	    {
	    	*distance = 1000000.0;
	    	return;
	    }
	    else if (l <= EPS3*scalar_prod)
	    {
	    	*distance = 100000.0;
	    	return;
	    }
	}

	*t =  scalar_prod / l;

	*distance = *norm_dist = (double) dist(*t);
	if (*t >= 1.0)
	    *distance = (double) (sqr(x - x1) + sqr(y - y1));
	else if (*t <= 0.0)
	    *distance = (double) (sqr(x_0) + sqr(y_0));
}		/*end shortest_dist*/


#endif /* defined(UNUSED_FUNCTION) */


/*
*			fd_grad_solution2d():
*
*    Returns x and y gradients of Locstate corresponding to x,y,comp.
*    by centered finite differences based on a regular stencil centered
*    on the point x,y.
*
*                     u[1][1]
*
*        u[0][-1]    u[0][0] = u[1][0]    u[0][1]
*
*                     u[1][-1]
*
*    NOTE the assumption of a regular stencil, ignoring any irregularity
*    that may actually be the stencil as returned by the calls to
*    tg_solution().
*/

LOCAL	boolean fd_grad_solution2d(
	double             *coords,
	COMPONENT         comp,
	register TRI_SOLN *soln,
	Locstate          *grad_ans)
{
	static boolean     first = YES;
	static int      dim;
	static double    h[MAXD], hi[MAXD];
	static double    **ncrds[MAXD];
	static Locstate *u[MAXD];
	INTERFACE       *gintfc = soln->tri_grid->grid_intfc;
	int           i, j;

	if (first)
	{
	    size_t   sizest = soln->sizest;
	    Locstate *tmpst_array;
	    Locstate tmpst;
	    double    **tmp;

	    first = NO;
	    dim = gintfc->dim;

	    alloc_state(gintfc,&tmpst,sizest);
	    for (i = 0; i < dim; ++i)
	    {
	        h[i] = soln->tri_grid->comp_grid.h[0];
	        hi[i] = 0.5 / h[i];
	        uni_array(&tmpst_array,3,sizeof(Locstate *));
	        u[i] = tmpst_array + 1;
	        alloc_state(gintfc,&u[i][-1],sizest);
	        u[i][0] = tmpst;
	        alloc_state(gintfc,&u[i][1],sizest);
	        uni_array(&tmp,3,sizeof(double **));
	        ncrds[i] = tmp + 1;
	        uni_array(&ncrds[i][-1],dim,FLOAT);
	        uni_array(&ncrds[i][1],dim,FLOAT);
	    }
	}

	for (i = 0; i < dim; ++i)
	{
	    ncrds[i][0] = coords;
	    ncrds[i][-1][i] = coords[i] - h[i];
	    ncrds[i][1][i] = coords[i] + h[i];
	    for (j = 0; j < dim; ++j)
	    {
	        if (j == i)
		    continue;
	        ncrds[i][-1][j] = coords[j];
	        ncrds[i][1][j] = coords[j];
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    for (j = -1; j < 2; j += 2)
	    {
	        if (!tg_solution(ncrds[i][j],comp,soln,u[i][j],NULL))
	        {
	            screen("ERROR in fd_grad_solution2d(), "
	                   "tg_solution() failed\n");
	            return FUNCTION_FAILED;
	        }
	    }
	}

	for (i = 0; i < dim; ++i)
	    bi_interpolate_intfc_states(gintfc,-hi[i],hi[i],
	                    ncrds[i][-1],u[i][-1],
	                    ncrds[i][1],u[i][1],grad_ans[i]);

	return FUNCTION_SUCCEEDED;
}        /*end fd_grad_solution2d*/



/*
*			flux2d():
*
*    Returns the line integral of the flux through that portion of
*    the mesh block ix,iy = icoords[0] + 1, icoords[1] + 1
*    associated with component comp. Loads the
*    integral into the Locstate "flx_ans".
*    Does the loading by calling the appropriate user supplied routines
*    for each boundary segment of the component area and summing
*    the results.
*
*    Vertical lines:
*        flux_across_grid_segment(tri_soln,0,state1,state2,dy,flx_ans)
*    Horizontal lines:
*        flux_across_grid_segment(tri_soln,1,state1,state2,dx,flx_ans)
*    Oblique lines:
*        flux_across_line_segment(tri_soln,state1,state2,ndir,flx_ans)
*
*    For each function call, state1 and state2 correspond to the
*    two end points of the boundary segment, ndir = dx, dy, dz is the signed
*    area (length) vector of the boundary segment, directed outwards.
*    The line integral is evaluated in the counter-clockwise direction.
*
*    A "typical state" is evaluated along with the flux. This
*    will either correspond to an interface state if the interface
*    cuts through the mesh block, otherwise to the state at the center
*    of the block.
*
*    If comp is not found in the mesh block, flx_ans and typ_state
*    are returned filled with zeros.
*/

LOCAL	void flux2d(
	int       *icoords,
	COMPONENT comp,
	TRI_SOLN  *soln,
	Locstate  flx_ans,
	Locstate  typ_state)
{
	BLK_EL0      *blk_el0;
	Locstate     s[4];
	TRI_GRID     *ntg = soln->tri_grid;
	RECT_GRID    *gr  = &ntg->comp_grid;
	double        ndir[MAXD];
	double        *h = gr->h;
	static double f[3] = {0.5, 0.5, 0.5};

	clear_state(ntg->grid_intfc,flx_ans,soln->sizest);
	clear_state(ntg->grid_intfc,typ_state,soln->sizest);

	blk_el0 = &Regular_blk_el0(icoords,ntg);

	if (blk_el0_is_bilinear(blk_el0))
	{
	    register BILINEAR_ELEMENT *eq;

	    eq = blk_el0_bilinear_el(blk_el0);
	    if (comp != eq->comp)
		return;

	    states_on_bilinear_element(s,eq,ntg);

	    flux_across_grid_segment(soln,1,s[0],s[1],h[0],flx_ans);
	    flux_across_grid_segment(soln,0,s[1],s[2],h[1],flx_ans);
	    flux_across_grid_segment(soln,1,s[2],s[3],-h[0],flx_ans);
	    flux_across_grid_segment(soln,0,s[3],s[0],-h[1],flx_ans);
	    Bilinear_cell_interpolate(f,eq,soln,typ_state);
	}
	else
	{
	    register LINEAR_ELEMENT    *et;
	    int             num_els, k;

	    num_els = num_lin_els_in_blk(blk_el0);
	    for (k = 0;  k < num_els;  ++k)
	    {
	        et = blk_el0_linear_els(blk_el0) + k;
	        if (et->comp == comp)
	        {
	            if (et->side[0] != I_SIDE)
	            {
	                ndir[0] =   Coords(et->p[1])[1] - Coords(et->p[0])[1];
	                ndir[1] = - Coords(et->p[1])[0] + Coords(et->p[0])[0];
	                flux_across_line_segment(soln,et->s[0],
	                                         et->s[1],ndir,flx_ans);
	            }
	            if (et->side[0] == F_SIDE)
	            {
	                bi_interpolate_intfc_states(
	                    soln->intfc,0.5,0.5,
	                    Coords(et->p[0]),et->s[0],
	                    Coords(et->p[1]),et->s[1],
	                    typ_state);
	            }

	            if (et->side[1] != I_SIDE)
	            {
	                ndir[0] =   Coords(et->p[2])[1] - Coords(et->p[1])[1];
	                ndir[1] = - Coords(et->p[2])[0] + Coords(et->p[1])[0];
	                flux_across_line_segment(soln,et->s[1],
	                                         et->s[2],ndir,flx_ans);
	            }
	            if (et->side[1] == F_SIDE)
	            {
	                bi_interpolate_intfc_states(soln->intfc,0.5,0.5,
	                                            Coords(et->p[1]),et->s[1],
	                                            Coords(et->p[2]),et->s[2],
	                                            typ_state);
	            }

	            if (et->side[2] != I_SIDE)
	            {
	                ndir[0] =   Coords(et->p[0])[1] - Coords(et->p[2])[1];
	                ndir[1] = - Coords(et->p[0])[0] + Coords(et->p[2])[0];
	                flux_across_line_segment(soln,et->s[2],
	                                         et->s[0],ndir,flx_ans);
	            }
	            if (et->side[2] == F_SIDE)
	            {
	                bi_interpolate_intfc_states(soln->intfc,0.5,0.5,
	                                            Coords(et->p[2]),et->s[2],
	                                            Coords(et->p[0]),et->s[0],
	                                            typ_state);
	            }
	        }
	    }
	}
}        /*end flux2d*/
#endif /* defined(TWOD) */
