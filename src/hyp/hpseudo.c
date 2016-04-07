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
*				hpseudo.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Computes a finite difference solution using locally generated stencils
*	about cell centers near fronts.
*/


#include <hyp/hdecs.h>

LOCAL	UnsplitStencil *AllocUnsplitStencil(Front*,Wave*);
LOCAL	boolean  regular_point(int*,COMPONENT,UnsplitStencil*);
LOCAL	boolean  set_state_from_adjacent_crosses(Locstate,COMPONENT,int*,double*,
	                                      UnsplitStencil*);
LOCAL	double SS2(const double*,const double*,const double*,int);
LOCAL	int   SetUnsplitStencilDataForStep(UnsplitStencil*,double,
					  Front*,Front*,Wave*);
LOCAL	void  FreeUnsplitStencilFrontAndWave(int,UnsplitStencil*);
LOCAL	void  LoadUnsplitStencil(int*,UnsplitStencil*);
LOCAL	void  update_state(int*,UnsplitStencil*);


LOCAL	UnsplitStencilOptions USopts = { NO };

IMPORT	const double *hs_coords_on;

EXPORT	void	h_set_unsplit_options(
	UnsplitStencilOptions *usopts)
{
	USopts = *usopts;
}		/*end h_set_unsplit_options*/

EXPORT int pseudo_unsplit_driver(
	double		dt,
	double		*dt_frac,
	Wave		*wave,
	Front		*front)
{
	COMPONENT              new_comp;
	Front	               *newfront;
	Front	               **infront, **outfront;
	INTERFACE              *current_intfc, *tmp_intfc;
	RECT_GRID              *gr = front->rect_grid;
	Wave	               **inwave,  **outwave;
	double	               *dh;
	double                  *dir;
	int	               status;
	int	               k, dim = gr->dim;
	int                    *iperm; /* permutation of {0,...,dim-1} */
	int                    imin[3], imax[3];
	int                    vsize;
	int                    *idirs;
	int                    icoords[3];
	int                    i0min, i1min, i2min;
	int                    i0max, i1max, i2max;
	int                    i0, i1, i2;
	int                    nrad = vsten_radius(wave);
	static UnsplitStencil  *ussten = NULL;
	static const char *warn = "WARNING in pseudo_unsplit_driver()";
	static const char *err = "ERROR in pseudo_unsplit_driver()";

	if (ussten == NULL)
	    ussten = AllocUnsplitStencil(front,wave);
	set_pt_source_interior_vectors(wave);

	front->interf->e_comps = NULL;

		/* Advance Front */

	if ((front->interf->modified) &&
	    (!make_interface_topology_lists(front->interf)))
	{
	    screen("ERROR in pseudo_unsplit_driver(), "
		   "can't make interface topology lists\n");
	    clean_up(ERROR);
	}
	start_clock("advance_front");
	status = advance_front(dt,dt_frac,front,&newfront,(POINTER) wave); 
	stop_clock("advance_front");
	if (status != GOOD_STEP)
	{
	    (void) printf("%s advance_front() failed, "
	                  "dt_frac = %g\n",warn,*dt_frac);
	    print_time_step_status("time step status = ",status,"\n");
	    return status;
	}
	initialize_max_wave_speed(wave);

	status = SetUnsplitStencilDataForStep(ussten,dt,front,newfront,wave);
	if (status != GOOD_STEP) 
	{
	    *dt_frac = min(*dt_frac,Min_time_step_modification_factor(front));
	    FreeUnsplitStencilFrontAndWave(status,ussten);
	    (void) printf("%s, SetUnsplitStencilDataForStep() failed\n",warn);
	    return status;
	}

	dh = ussten->dh;
	iperm = ussten->iperm;
	inwave = ussten->inwave;   outwave = ussten->outwave;
	infront = ussten->infront; outfront = ussten->outfront;


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

#if defined(TIME_HVEC)
	start_clock("Finite difference step");
#endif /* defined(TIME_HVEC) */

	for (k = 0; k < dim; ++k)
	{
	    idirs = ussten->idirs[k];
	    dir = ussten->dirs[idirs[0]];
	    set_sweep_limits(inwave[k],k,idirs,imin,imax);
	    i0min = imin[0];
	    i0max = imax[0];
	    if (gr->lbuf[idirs[0]] > 0)
		i0min -= nrad;
	    if (gr->ubuf[idirs[0]] > 0)
	        i0max += nrad;
	    vsize = i0max - i0min;
	    alloc_phys_vecs(inwave[k],vsize);
	    icoords[idirs[0]] = 0;

	    switch (dim)
	    {
	    case 1:
	        vec_solver(k,iperm,dir,inwave[k],outwave[k],
			   infront[k],outfront[k],icoords,i0min,i0max,dt,dh[k]);
	        break;
#if defined(TWOD)
	    case 2:
	        i1min = imin[1];
		i1max = imax[1];
	        for (i1 = i1min; i1 < i1max; ++i1)
	        {
	            icoords[idirs[1]] = i1;
	            vec_solver(k,iperm,dir,inwave[k],outwave[k],infront[k],
			       outfront[k],icoords,i0min,i0max,dt,dh[k]);		
	        }
	        break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
	        i1min = imin[1];
		i1max = imax[1];
	        i2min = imin[2];
		i2max = imax[2];
	        for (i2 = i2min; i2 < i2max; ++i2)
	        {
	            icoords[idirs[2]] = i2;
	   	    for (i1 = i1min; i1 < i1max; ++i1)
	   	    {
	                icoords[idirs[1]] = i1;
	                vec_solver(k,iperm,dir,inwave[k],outwave[k],infront[k],
				   outfront[k],icoords,i0min,i0max,dt,dh[k]);
	   	    }
	        }
	        break;
#endif /* defined(THREED) */
	    }
	    free_phys_vecs(inwave[k]);
	}

#if defined(TIME_HVEC)
	stop_clock("Finite difference step");
#endif /* defined(TIME_HVEC) */

	/*  UPDATE THE IRREGULAR GRID */

#if defined(TIME_HVEC)
	start_clock("sweep on irregular grid");
#endif /* defined(TIME_HVEC) */
	i0min = 0; i0max = gr->gmax[0];
	i1min = 0; i1max = (dim > 1) ? gr->gmax[1] : 1;
	i2min = 0; i2max = (dim > 2) ? gr->gmax[2] : 1;
	for (i2 = i2min; i2 < i2max; ++i2)
	{
	    icoords[2] = i2;
	    for (i1 = i1min; i1 < i1max; ++i1)
	    {
	        icoords[1] = i1;
	        for (i0 = i0min; i0 < i0max; ++i0)
	        {
	            icoords[0] = i0;
		    new_comp = Rect_comp(icoords,outwave[dim-1]);
		    if (!regular_point(icoords,new_comp,ussten))
		    {
		        ussten->new_comp = new_comp;
		        update_state(icoords,ussten);
		    }
	        }
	    }
	}

#if defined(TIME_HVEC)
	stop_clock("sweep on irregular grid");
#endif /* defined(TIME_HVEC) */

	if (dim == 1)
	    free_wave_pointers(wave);
	set_current_interface(current_intfc);
	if (!delete_interface(tmp_intfc))
	{
	    screen("%s, copy_hyp_solution_function() failed\n",err);
	    FreeUnsplitStencilFrontAndWave(ERROR_IN_STEP,ussten);
	    return ERROR_IN_STEP;
	}

	FreeUnsplitStencilFrontAndWave(GOOD_STEP,ussten);

	stop_clock("hyp_solver");

	/* parallel part for interior states */

	for (k = 0; k < dim; ++k)
        {
	    if (scatter_states(wave,front,iperm,k) != FUNCTION_SUCCEEDED)
	    {
	        screen("%s, scatter_states() failed\n",err);
	        clean_up(ERROR);
	    }
	}

	front->interf->e_comps = NULL;

	free_pt_source_interior_vectors(wave);

	return status;
}		/*end pseudo_unsplit_driver*/

LOCAL	int	SetUnsplitStencilDataForStep(
	UnsplitStencil *ussten,
	double          dt,
	Front          *front,
	Front          *newfront,
	Wave           *wave)
{
	Front	          **infront, **outfront;
	Wave	          **inwave,  **outwave;
	int               *iperm; /* permutation of {0,...,dim-1} */
	int               i, k, dim = front->rect_grid->dim;
	int	          step = front->step;
	int               status;
	static const char *warn = "WARNING in SetUnsplitStencilDataForStep()";

	ussten->max_comp = max_component(front->interf);
	ussten->dt = dt;
	iperm = ussten->iperm = set_iperm(step,dim);

	for (k = 0; k < dim; ++k)
	{
	    ussten->dh[k] = wave->rect_grid->h[iperm[k]];
	    for (i = 0; i < dim; ++i)
		ussten->idirs[k][i] = iperm[(i+k)%dim];
	}

	ussten->front = front;
	ussten->newfront = newfront;
	ussten->wave  = wave;

		/* Initialize Intermediate Storage for States */

	inwave = ussten->inwave;   outwave = ussten->outwave;
	infront = ussten->infront; outfront = ussten->outfront;
	for (i = 0; i < 3; ++i)
	{
	     inwave[i] = NULL;   outwave[i] = NULL;
	    infront[i] = NULL;  outfront[i] = NULL;
	}
	ussten->wk_wv1 = ussten->wk_wv2 = ussten->tmpwave = NULL;

	ussten->wk_wv1 = copy_wave(wave);
	clear_wave_pointers(ussten->wk_wv1);

	switch (dim)
	{
	case 1:
	    inwave[0]  = wave;  outwave[0]  = ussten->wk_wv1;
	    infront[0] = front;	outfront[0] = newfront;
	    break;
	case 2:
	    ussten->wk_wv2 = copy_wave(wave);
	    clear_wave_pointers(ussten->wk_wv2);
	    inwave[0]  = wave;	     outwave[0]  = ussten->wk_wv1;
	    inwave[1]  = outwave[0]; outwave[1]  = ussten->wk_wv2;
	    infront[0] = front;	     outfront[0] = newfront;
	    infront[1] = newfront;   outfront[1] = newfront;
	    ussten->tmpwave = ussten->wk_wv1;
	    break;
	case 3:
	    ussten->wk_wv2 = copy_wave(wave);
	    clear_wave_pointers(ussten->wk_wv2);
	    inwave[0] = wave;	    outwave[0] = ussten->wk_wv2;
	    inwave[1] = outwave[0]; outwave[1] = ussten->wk_wv1;
	    inwave[2] = outwave[1]; outwave[2] = ussten->wk_wv2;
	    ussten->tmpwave = ussten->wk_wv1;
	    infront[0] = front;	    outfront[0] = newfront;
	    infront[1] = newfront;  outfront[1] = newfront;
	    infront[2] = newfront;  outfront[2] = newfront;
	    break;
	}

	start_clock("init_hyp_solution");
	assign_wave_parameters(outwave[0],wave);
	status = init_hyp_solution_function(outwave[0],newfront);
	status = syncronize_time_step_status(status,front->pp_grid);
	if (status != GOOD_STEP) 
	{
	    (void) printf("%s, init_hyp_solution_function() failed\n",warn);
	    return status;
	}
	if (dim > 1)
	{
	    assign_wave_parameters(outwave[1],outwave[0]);
	    if (!copy_hyp_solution_function(outwave[0],outwave[1]))
	    {
	        (void) printf("%s, copy_hyp_solution_function() failed\n",warn);
		return ERROR_IN_STEP;
	    }
	}
	stop_clock("init_hyp_solution");
	return status;
}		/*end SetUnsplitStencilDataForStep*/

LOCAL	void	FreeUnsplitStencilFrontAndWave(
	int            status,
	UnsplitStencil *ussten)
{
	int i, dim = ussten->front->rect_grid->dim;
	if (status == GOOD_STEP)
	{
		/* Copy updated front, wave */
	    assign_wave_pointers(ussten->wave,ussten->outwave[dim-1]);

	    /* Free temporary storage, update front */

	    if (ussten->tmpwave != NULL)
	        free_copy_wave_pointers(ussten->tmpwave);

	    assign_interface_and_free_front(ussten->front,ussten->outfront[0]);
	}
	else
	{
	    free_front(ussten->outfront[0]);
	    free_wave_pointers(ussten->outwave[0]);
	}
	if (ussten->wk_wv1 != NULL)
	    free_wave(ussten->wk_wv1);
	if (ussten->wk_wv2 != NULL)
	    free_wave(ussten->wk_wv2);
	for (i = 0; i < 3; ++i)
	{
	     ussten->inwave[i] = NULL;   ussten->outwave[i] = NULL;
	    ussten->infront[i] = NULL;  ussten->outfront[i] = NULL;
	}
	ussten->wk_wv1 = ussten->wk_wv2 = ussten->tmpwave = NULL;
}		/*end FreeUnsplitStencilFrontAndWave*/

LOCAL	UnsplitStencil	*AllocUnsplitStencil(
	Front *front,
	Wave  *wave)
{
	ALIGN          *buf;
	UnsplitStencil *ussten;
	int            i, j, len;
	int            dim = front->rect_grid->dim;
	int            indx;
	int            npts, nrad0, nrad1, nrad2;
	size_t         sizest = wave->sizest;
	size_t         size;
	size_t         NaUsS;
	size_t         NaSts, NaStStore;
	size_t         NaCoords, NaCoordsStore;
	size_t         NaIc, NaIcStore;
	size_t         NaComp;
	size_t         offset;
	size_t         NaStateSet;

	npts = wave->npts_sten;
	for (len = 1, i = 0; i < dim; ++i)
	    len *= npts;

	NaUsS = num_aligns(sizeof(UnsplitStencil));
	NaSts = num_aligns(sizeof(Locstate)*len);
	NaCoords = num_aligns(sizeof(double*)*len);
	NaStStore = num_aligns(sizest*len);
	NaComp = num_aligns(sizeof(COMPONENT)*len);
	NaCoordsStore = num_aligns(3*sizeof(double)*len);
	NaIc = num_aligns(sizeof(int*)*len);
	NaIcStore = num_aligns(3*sizeof(int)*len);
	NaStateSet = num_aligns(sizeof(boolean)*len);
	size = sizeof(ALIGN)*(NaUsS + NaSts + NaCoords + NaStStore +
			      NaComp + NaCoordsStore +
			      NaIc + NaIcStore + NaStateSet);
	scalar(&buf,size);
	ussten = (UnsplitStencil*)buf;
	ussten->npts = npts;
	ussten->nrad[0] = nrad0 = ussten->rad = stencil_radius(wave);
	ussten->nrad[1] = nrad1 = (dim > 1) ? ussten->rad : 0;
	ussten->nrad[2] = nrad2 = (dim > 2) ? ussten->rad : 0;
	offset = NaUsS;
	ussten->state_list = (Locstate*)(buf+offset);
	offset += NaSts;
	ussten->coords_list = (double**)(buf+offset);
	offset += NaCoords;
	ussten->state_store = (byte*)(buf+offset);
	offset += NaStStore;
	ussten->coords_store = (double*)(buf+offset);
	offset += NaCoordsStore;
	ussten->comp_list = (COMPONENT*)(buf+offset);
	offset += NaComp;
	ussten->ic_list = (int**)(buf+offset);
	offset += NaIc;
	ussten->ic_store = (int*)(buf+offset);
	offset += NaIcStore;
	ussten->state_set_list = (boolean*)(buf+offset);
	for (i = 0; i < len; ++i)
	{
	    ussten->state_list[i] = (Locstate)(ussten->state_store + i*sizest);
	    ussten->coords_list[i] = ussten->coords_store + 3*i;
	    ussten->ic_list[i] = ussten->ic_store + 3*i;
	}
	indx = nrad0 + npts*(nrad1 + npts*nrad2);

	/* Set arrays to offset from stencil center point */

	ussten->comp = ussten->comp_list + indx;
	ussten->state = ussten->state_list + indx;
	ussten->coords = ussten->coords_list + indx;
	ussten->ic = ussten->ic_list + indx;
	ussten->state_set = ussten->state_set_list + indx;

	for (i = 0; i < 3; ++i)
	{
	    ussten->idirs[i] = ussten->idir_store+3*i;
	    ussten->dirs[i] = ussten->dir_store+3*i;
	    for (j = 0; j < 3; ++j)
		ussten->dirs[i][j] = 0.0;
	    ussten->dirs[i][i] = 1.0;
	}

	ussten->sten = alloc_stencil(npts,front);
	return ussten;
}		/*end AllocUnsplitStencil*/

/*
*			regular_point():
*
*	Determines whether front points lie within the domain of dependency
*	of a given coordinate on the computation grid.   Returns YES if no
*	front points influence the updated state at icoords,  NO otherwise.
*/

LOCAL	boolean regular_point(
	int             *icoords,
	COMPONENT       nc,
	UnsplitStencil  *ussten)
{
	Wave             *wv = ussten->wave;
	INTERFACE        *intfc = ussten->newfront->interf;
	RECT_GRID        *gr = wv->rect_grid;
	int              *lbuf = gr->lbuf, *ubuf = gr->ubuf;
	int              *gmax = gr->gmax;
	int              i, j, k, rad = stencil_radius(wv);
	int              imin, imax, jmin, jmax, kmin, kmax;
	int              dim = gr->dim;
	int              ic[3];

	/* Check for near boundary */

	for (i = 0; i < dim; ++i)
	{
	    if ((lbuf[i] == 0) && (icoords[i] < rad))
		return NO;
	    if ((ubuf[i] == 0) && (gmax[i] <= (icoords[i] + rad)))
		return NO;
	}
	imin = icoords[0] - rad;
	imax = icoords[0] + rad;
	if (dim > 1)
	{
	    jmin = icoords[1] - rad;
	    jmax = icoords[1] + rad;
	}
	else
	{
	    jmin = 0;
	    jmax = 0;
	}
	if (dim > 2)
	{
	    kmin = icoords[2] - rad;
	    kmax = icoords[2] + rad;
	}
	else
	{
	    kmin = 0;
	    kmax = 0;
	}
	for (k = kmin; k <= kmax; ++k)
	{
	    ic[2] = k;
	    for (j = jmin; j <= jmax; ++j)
	    {
		ic[1] = j;
	        for (i = imin; i <= imax; ++i)
		{
		    ic[0] = i;
		    if (!equivalent_comps(Rect_comp(ic,wv),nc,intfc))
			return NO;
		    if ((k < kmax) && (Rect_crossing(ic,UPPER,wv) != NULL))
			return NO;
		    if ((j < jmax) && (Rect_crossing(ic,NORTH,wv) != NULL))
			return NO;
		    if ((i < imax) && (Rect_crossing(ic,EAST,wv) != NULL))
			return NO;
		}
	    }
	}
	return YES;
}		/*end regular_point*/

/*
*			update_state():
*
*	Updates the state at an irregular (ie near front) grid location.
*
*	Input:
*	       Fields that vary with position
*		icoords
*		ussten->new_comp
*              Fields that are constant across a time step
*		ussten->max_comp
*		ussten->dt
*		ussten->iperm
*		ussten->dh
*		ussten->idirs
*		ussten->front
*		ussten->wave
*		ussten->inwave
*		ussten->outwave
*		ussten->infront
*		ussten->outfront
*
*	Output:
*		ussten->state array
*		ussten->coords array
*		ussten->ic array
*/
LOCAL	void	update_state(
	int            *icoords,
	UnsplitStencil *ussten)
{
	COMPONENT  new_comp = ussten->new_comp, max_comp = ussten->max_comp;
	Front	   **infront = ussten->infront, **outfront = ussten->outfront;
	Front      *front = ussten->front;
	HYPER_SURF *hs;
	Locstate   *st, *state;
	POINT      **p;
	RECT_GRID  *gr = front->rect_grid;
	Stencil    *sten = ussten->sten;
	Wave	   **inwave = ussten->inwave,  **outwave = ussten->outwave;
	Wave       *wave = ussten->wave;
	double      **coords, *crds;
	double      *dh = ussten->dh;
	double      *dir;
	double      crds_grid[3], coords_on[3];
	double      dt = ussten->dt;
	int        **icrds;
	int        *iperm = ussten->iperm;
	int        dim = gr->dim;
	int        i, i0, i1, i2;
	int        ic[3];
	int        indx, index;
	int        npts = ussten->npts;
	int        rad, nrad0, nrad1, nrad2;
	size_t     sizest = wave->sizest;

	if (new_comp > max_comp)
	{
	    double ds;
	    for (i = 0; i < dim ;++i)
		crds_grid[i] = cell_center(icoords[i],i,gr);
	    nearest_intfc_state_and_pt(crds_grid,new_comp,outfront[dim-1],front,
				       Rect_state(icoords,outwave[dim-1]),
				       coords_on,&hs);
	    if ((ds=SS2(crds_grid,coords_on,gr->h,dim)) > dim)
	    {
		screen("ERROR in update_state(), answer taken from remote "
		       "interface point\n");
		print_int_vector("icoords = ",icoords,dim,", ");
		print_general_vector("crds_grid = ",crds_grid,dim,", ");
		print_general_vector("coords_on = ",coords_on,dim,", ");
		(void) printf("ds2 = %g, ds = %g\n",ds,sqrt(ds));
		clean_up(ERROR);
	    }
	    return;
	}

	LoadUnsplitStencil(icoords,ussten);

	rad = ussten->rad;
	nrad2 = (dim > 2) ? ussten->nrad[iperm[2]] : 0;
	nrad1 = (dim > 1) ? ussten->nrad[iperm[1]] : 0;
	nrad0 = ussten->nrad[iperm[0]];
	state = ussten->state;
	coords = ussten->coords;

	sten->reg_stencil = sten->prev_reg_stencil = NO;
	sten->npts = npts;
	icrds = sten->icoords;
	sten->newcomp = new_comp;
	for (i = -rad; i <= rad; ++i)
	{
	    sten->comp[i] = new_comp;
	    sten->nc[i] = 0;
	    sten->crx[i][0] = NULL;
	    sten->hs[i] = NULL;
	    sten->st[i] = sten->worksp[i];
	}
	st = sten->st;
	p = sten->p;

	/* First sweep */
	sten->fr = infront[0];
	sten->newfr = outfront[0];
	sten->wave = inwave[0];
	sten->newwave = outwave[0];
	dir = &ussten->dirs[iperm[0]][0];

	for (i2 = -nrad2; i2 <= nrad2; ++i2)
	{
	    ic[iperm[2]] = i2;
	    if (dim > 2)
	    {
	        for (i0 = -nrad0; i0 <= nrad0; ++i0)
	            icrds[i0][iperm[2]] = icoords[iperm[2]] + i2;
	    }
	    for (i1 = -nrad1; i1 <= nrad1; ++i1)
	    {
	        ic[iperm[1]] = i1;
		if (dim > 1)
		{
	            for (i0 = -nrad0; i0 <= nrad0; ++i0)
	                icrds[i0][iperm[1]] = icoords[iperm[1]] + i1;
		}
	        for (i0 = -nrad0; i0 <= nrad0; ++i0)
		{
	            ic[iperm[0]] = i0;
	            icrds[i0][iperm[0]] = icoords[iperm[0]] + i0;
		    indx = ic[0] + npts*(ic[1] + npts*ic[2]);
		    ft_assign(st[i0],state[indx],sizest);
		    crds = coords[indx];
		    for (i = 0; i < dim; ++i)
			Coords(p[i0])[i] = crds[i];
		    ft_assign(left_state(sten->p[i0]),st[i0],sizest);
		    ft_assign(right_state(sten->p[i0]),st[i0],sizest);
		}
		index = is_source_block(wave,outfront[0]->interf,new_comp,
					icrds[0]);
		ic[iperm[0]] = 0;
		indx = ic[0] + npts*(ic[1] + npts*ic[2]);
		npt_solver(dh[0],dt,state[indx],dir,0,iperm,&index,sten,wave);
	    }
	}
	if (dim == 1)
	{
	    ft_assign(Rect_state(icoords,outwave[0]),state[0],sizest);
	    return;
	}

	/* Second sweep */
	dir = &ussten->dirs[iperm[1]][0];
	ic[iperm[0]] = 0;
	for (i1 = -nrad1; i1 <= nrad1; ++i1)
	    icrds[i1][iperm[0]] = icoords[iperm[0]];
	for (i2 = -nrad2; i2 <= nrad2; ++i2)
	{
	    ic[iperm[2]] = i2;
	    if (dim > 2)
	    {
	        for (i1 = -nrad0; i1 <= nrad1; ++i1)
	            icrds[i1][iperm[2]] = icoords[iperm[2]] + i2;
	    }
	    for (i1 = -nrad1; i1 <= nrad1; ++i1)
	    {
	        ic[iperm[1]] = i1;
		indx = ic[0] + npts*(ic[1] + npts*ic[2]);
		ft_assign(st[i1],state[indx],sizest);
		crds = coords[indx];
		for (i = 0; i < dim; ++i)
		    Coords(p[i1])[i] = crds[i];
		ft_assign(left_state(sten->p[i1]),st[i1],sizest);
		ft_assign(right_state(sten->p[i1]),st[i1],sizest);
	    }
	    index = is_source_block(wave,outfront[1]->interf,new_comp,icrds[0]);
	    ic[iperm[1]] = 0;
	    indx = ic[0] + npts*(ic[1] + npts*ic[2]);
	    npt_solver(dh[1],dt,state[indx],dir,1,iperm,&index,sten,wave);
	}

	if (dim == 2)
	{
	    ft_assign(Rect_state(icoords,outwave[1]),state[0],sizest);
	    return;
	}

	/* Third sweep */
	dir = &ussten->dirs[iperm[2]][0];
	ic[iperm[1]] = 0;
	for (i2 = -nrad2; i2 <= nrad2; ++i2)
	{
	    icrds[i2][iperm[1]] = icoords[iperm[1]];
	    icrds[i2][iperm[2]] = icoords[iperm[2]] + i2;
	    ic[iperm[2]] = i2;
	    indx = ic[0] + npts*(ic[1] + npts*ic[2]);
	    ft_assign(st[i2],state[indx],sizest);
	    crds = coords[indx];
	    for (i = 0; i < dim; ++i)
		Coords(p[i2])[i] = crds[i];
	    ft_assign(left_state(sten->p[i2]),st[i2],sizest);
	    ft_assign(right_state(sten->p[i2]),st[i2],sizest);
	}

	index = is_source_block(wave,outfront[2]->interf,new_comp,icrds[0]);
	
	npt_solver(dh[2],dt,state[0],dir,2,iperm,&index,sten,wave);
	ft_assign(Rect_state(icoords,outwave[2]),state[0],sizest);
}		/*end update_state*/

/*
*			LoadUnsplitStencil():
*
*	Loads state, component,  and position data in the UnsplitStencil
*	structure ussten.
*
*	Input:
*	       Fields that vary with position
*		icoords
*		ussten->new_comp
*              Fields that are constant across a time step
*		ussten->max_comp
*		ussten->dt
*		ussten->iperm
*		ussten->dh
*		ussten->idirs
*		ussten->front
*		ussten->wave
*		ussten->inwave
*		ussten->outwave
*		ussten->infront
*		ussten->outfront
*
*	Output:
*		ussten->state array
*		ussten->coords array
*		ussten->ic array
*/

LOCAL	void	LoadUnsplitStencil(
	int            *icoords,
	UnsplitStencil *ussten)
{
	COMPONENT  new_comp = ussten->new_comp;
	COMPONENT  *comp;
	COMPONENT  ecomp;
	Front      *front = ussten->front;
	HYPER_SURF *hs;
	INTERFACE  *intfc = ussten->front->interf;
	INTERFACE  *new_intfc = ussten->newfront->interf;
	RECT_GRID  *gr = front->rect_grid;
	Locstate   *state;
	Wave       *wave = ussten->wave;
	boolean       *state_set;
	double      *h = gr->h;
	double      **coords;
	double      coords_on[3];
	double      y, z;
	int        *iperm = ussten->iperm;
	int        dim = intfc->dim;
	int        i, j, k, isgn, jsgn, ksgn;
	int        ic[3], icp[3], **ussten_ic;
	int        imin, imax, jmin, jmax, kmin, kmax;
	int        indx, indxp;
	int        npts = ussten->npts;
	int        rad, nrad0, nrad1, nrad2;
	int        xmin, xmax, ymin, ymax, zmin, zmax;
	size_t     sizest = wave->sizest;

	rad = ussten->rad;
	nrad0 = ussten->nrad[0];
	nrad1 = ussten->nrad[1];
	nrad2 = ussten->nrad[2];
	imin  = icoords[0] - nrad0;
	imax  = icoords[0] + nrad0;
	xmin = -gr->lbuf[0];
	xmax =  gr->gmax[0] + gr->ubuf[0];
	if (dim > 1)
	{
	    jmin = icoords[1] - nrad1;
	    jmax = icoords[1] + nrad1;
	    ymin = -gr->lbuf[1];
	    ymax =  gr->gmax[1] + gr->ubuf[1];
	}
	else
	{
	    jmin = 0;
	    jmax = 0;
	    ymin = 0;
	    ymax = 1;
	}
	if (dim > 2)
	{
	    kmin = icoords[2] - nrad2;
	    kmax = icoords[2] + nrad2;
	    zmin = -gr->lbuf[2];
	    zmax =  gr->gmax[2] + gr->ubuf[2];
	}
	else
	{
	    kmin = 0;
	    kmax = 0;
	    zmin = 0;
	    zmax = 1;
	}

	/*
	 *  Initialize coords, comps,  and icoords.  Set states whose component
	 *  agrees with new_comp or that have an adjacent crossing with the
	 *  component new_comp.
	 */
	ecomp = exterior_component(intfc);
	state = ussten->state_list;
	coords = ussten->coords_list;
	comp = ussten->comp_list;
	ussten_ic = ussten->ic_list;
	state_set = ussten->state_set_list;
	for (k = kmin; k <= kmax; ++k)
	{
	    z = (dim > 2) ? cell_center(k,2,gr) : 0.0;
	    for (j = jmin; j <= jmax; ++j)
	    {
	        y = (dim > 1) ? cell_center(j,1,gr) : 0.0;
	        for (i = imin; i <= imax; ++i)
		{
		    (*ussten_ic)[0] = i;
		    (*ussten_ic)[1] = j;
	            (*ussten_ic)[2] = k;
	            (*coords)[0] = cell_center(i,0,gr);
	            (*coords)[1] = y;
	            (*coords)[2] = z;
		    if ((k < zmin) || (k >= zmax))
		        *comp = ecomp;
		    else if ((j < ymin) || (j >= ymax))
		        *comp = ecomp;
		    else if ((i < xmin) || (i >= xmax))
			*comp = ecomp;
		    else
			*comp = Rect_comp(*ussten_ic,wave);
		    if (equivalent_comps(*comp,new_comp,new_intfc))
		    {
			ft_assign(*state,Rect_state(*ussten_ic,wave),sizest);
			*state_set = YES;
		    }
		    else if (USopts.use_hyp_solution)
		    {
			double ds;
			hyp_solution(*coords,new_comp,NULL,UNKNOWN_SIDE,
				     front,wave,*state,NULL);
			*state_set = YES;
	                if ((ds=SS2(*coords,hs_coords_on,h,dim)) > dim)
	                {
		            screen("ERROR in LoadUnsplitStencil(), hyp_solution "
				    "state taken from remote location\n");
		            print_int_vector("icoords = ",*ussten_ic,dim,", ");
		            print_general_vector("coords = ",*coords,dim,", ");
		            print_general_vector("hs_coords_on = ",hs_coords_on,
				                 dim,", ");
		            (void) printf("ds2 = %g, ds = %g\n",ds,sqrt(ds));
		            clean_up(ERROR);
	                }
		    }
		    else
		    {
			*state_set = set_state_from_adjacent_crosses(*state,
								     new_comp,
								     *ussten_ic,
								     *coords,
								     ussten);
		    }
		    ++state;
		    ++coords;
		    ++comp;
		    ++ussten_ic;
		    ++state_set;
		}
	    }
	}

	/* Set mid state if not already set */
	state = ussten->state;
	state_set = ussten->state_set;
	coords = ussten->coords;
	if (state_set[0] == NO)
	{
	    nearest_intfc_state_and_pt(coords[0],new_comp,front,NULL,
				       state[0],coords_on,&hs);
	    state_set[0] = YES;
	}

	/*
	 *  Set states along the axes of the stencil center that are still
	 *  unset by copying the inner adjacent state.
	 */

	for (i = 0; i < dim; ++i)
	{
	    ic[0] = ic[1] = ic[2] = 0;
	    icp[0] = icp[1] = icp[2] = 0;
	    for (j = 1; j <= rad; ++j)
	    {
		for (jsgn = -1; jsgn < 2; jsgn += 2)
		{
		    ic[i] = jsgn*j;
	            indx = ic[0] + npts*(ic[1] + npts*ic[2]);
		    if (state_set[indx] == NO)
		    {
		        icp[i] = jsgn*(j-1);
	                indxp = icp[0] + npts*(icp[1] + npts*icp[2]);
		        ft_assign(state[indx],state[indxp],sizest);
		        state_set[indx] = YES;
		    }
		}
	    }
	}

	if (dim > 1)
	{
	    ic[0] = ic[1] = ic[2] = 0;
	    icp[0] = icp[1] = icp[2] = 0;
	    for (i = 1; i <= rad; ++i)
	    {
		for (isgn = -1; isgn < 2; isgn += 2)
		{
		    icp[iperm[dim-1]]  = ic[iperm[dim-1]] = isgn*i;
		    for (j = 1; j <= rad; ++j)
		    {
			for (jsgn = -1; jsgn < 2; jsgn += 2)
			{
		            ic[iperm[dim-2]] = jsgn*j;
	                    indx = ic[0] + npts*(ic[1] + npts*ic[2]);
		            if (state_set[indx] == NO)
		            {
		                icp[iperm[dim-2]] = jsgn*(j-1);
	                        indxp = icp[0] + npts*(icp[1] + npts*icp[2]);
		                ft_assign(state[indx],state[indxp],sizest);
		                state_set[indx] = YES;
		            }
			}
		    }
		}
	    }
	}

	if (dim > 2)
	{
	    ic[0] = ic[1] = ic[2] = 0;
	    icp[0] = icp[1] = icp[2] = 0;
	    for (i = -rad; i <= rad; ++i)
	    {
		icp[iperm[2]]  = ic[iperm[2]] = i;
		for (j = -rad; j <= rad; ++j)
		{
		    icp[iperm[1]]  = ic[iperm[1]] = j;
		    for (k = 1; k <= rad; ++k)
		    {
		        for (ksgn = -1; ksgn < 2; ksgn += 2)
		        {
		            ic[iperm[0]] = ksgn*k;
	                    indx = ic[0] + npts*(ic[1] + npts*ic[2]);
		            if (state_set[indx] == NO)
		            {
		                icp[iperm[0]] = ksgn*(k-1);
	                        indxp = icp[0] + npts*(icp[1] + npts*icp[2]);
		                ft_assign(state[indx],state[indxp],sizest);
		                state_set[indx] = YES;
		            }
	                }
		    }
		}
	    }
	}
}		/*end LoadUnsplitStencil*/

LOCAL	boolean  set_state_from_adjacent_crosses(
	Locstate  state,
	COMPONENT new_comp,
	int       *icoords,
	double     *coords,
	UnsplitStencil *ussten)
{
	Front            *front = ussten->front;
	Front            *newfront = ussten->newfront;
	Wave             *wave = ussten->wave;
	INTERFACE        *intfc = front->interf;
	INTERFACE        *new_intfc = newfront->interf;
	CRXING           *closest_crx, *crx[MAX_NUM_CRX];
	HYPER_SURF       *hs, *chs;
	TRI_GRID         *tri_grid = wave_tri_soln(wave)->tri_grid;
	double            dist, min_dist;
	Locstate         st;
	int              i, j, k, dim = front->rect_grid->dim;
	int              ncross;
	static const GRID_DIRECTION crx_dir[3][2] = { {  EAST,  WEST },
					              { NORTH, SOUTH },
					              { UPPER, LOWER } };

	min_dist = HUGE_VAL;
	closest_crx = NULL;
	for (i = 0; i < dim; ++i)
	{
	    for (j = 0; j < 2; ++j)
	    {
		ncross = crossings_in_direction(crx,icoords,crx_dir[i][j],
						tri_grid);
		for (k = 0; k < ncross; ++k)
		{
		    if (equivalent_comps(positive_component(crx[k]->hs),
					  new_comp,new_intfc) ||
		        equivalent_comps(negative_component(crx[k]->hs),
					  new_comp,new_intfc))
		    {
			dist = distance_between_positions(Coords(crx[k]->pt),
							  coords,dim);
			if (dist < min_dist)
			{
			    min_dist = dist;
			    closest_crx = crx[k];
			}
		    }
		}
	    }
	}
	if (closest_crx == NULL)
	    return NO;
	hs = closest_crx->hs;
	if ((wave_type(hs) == NEUMANN_BOUNDARY) &&
	    (front->neumann_bdry_state))
	{
	    chs = find_correspond_hyper_surface(hs,NULL,NULL,front,intfc);
	    if (chs == NULL)
		chs = hs;
	    if (!(*front->neumann_bdry_state)(coords,new_comp,
					      closest_crx->pt,
				              chs,front,(POINTER)wave,state))
		return NO;
	}
	else if (wave_type(hs) == DIRICHLET_BOUNDARY)
	{
	    chs = find_correspond_hyper_surface(hs,NULL,NULL,front,intfc);
	    if (chs == NULL)
		chs = hs;
	    evaluate_dirichlet_boundary_state(coords,chs,front,wave,state);
	}
	else
	{
	    st = state_with_comp(closest_crx->pt,hs,new_comp);
	    if (st != NULL)
	        ft_assign(state,st,front->sizest);
	    else
		return NO;
	}
	return YES;
}		/*end set_state_from_adjacent_crosses*/

LOCAL	double SS2(
	const double *p,
	const double *q,
	const double *h,
	int dim)
{
    	double ds2;
    	int i;
	for (ds2 = 0.0, i = 0; i < dim; ++i)
	    ds2 += ((p[i]-q[i])/h[i])*((p[i]-q[i])/h[i]);
	return ds2;
}		/*end SS2*/

