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
*				hsub.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines that interface to the grid construction
*	and solution function algorthm(s).  Also contains some
*	miscellaneous subroutines for the hyperbolic library.
*/


#include <hyp/hdecs.h>

	/* LOCAL Function Prototypes */
LOCAL	Front	*h_copy_front(Front*);
#if defined(USE_OVERTURE)
LOCAL   Front   *h_deep_copy_front(Front*);
#endif /* if defined(USE_OVERTURE) */


/*ARGSUSED*/
EXPORT	void	set_sweep_limits(
	Wave		*wv,
	int		swp_num,
	int		*idirs,
	int		*imin,
	int		*imax)
{
	RECT_GRID	*rgr = wv->rect_grid;
	int		*gmax = rgr->gmax;
	int		dim = rgr->dim;
	int		*lbuf = rgr->lbuf;
	int		*ubuf = rgr->ubuf;
	int		i;
	int		nrad = dim - 1 - swp_num;

	for (i = 0; i < dim; ++i)
	{
	    imin[i] = 0;
	    imax[i] = gmax[idirs[i]];
	}
	for (i = dim; i < 3; ++i)
	{
	    imin[i] = 0;
	    imax[i] = 1;
	}
	for (i = 1; i < dim; ++i)
	{
	    if (lbuf[idirs[i]] > 0)
		imin[i] -= lbuf[idirs[i]];
	    if (ubuf[idirs[i]] > 0)
		imax[i] += ubuf[idirs[i]];
	}
        /* 042603 added */
#if defined(USE_OVERTURE)
        if (lbuf[idirs[0]] == 4)
            imin[0] -= nrad;
        if (ubuf[idirs[0]] == 4)
            imax[0] += nrad;
#else
        if (lbuf[idirs[0]] > 0)
            imin[0] -= nrad;
        if (ubuf[idirs[0]] > 0)
            imax[0] += nrad;
#endif /* if defined(USE_OVERTURE) */


	if (debugging("sweep_limits"))
	{
	    static const char *var[] = {"ix","iy","iz"};
	    int k;

	    (void) printf("\nSweep limits\n");
	    (void) printf("\tdir = %d, ",idirs[0]);
	    print_int_vector("idirs = ",idirs,dim," ");
	    (void) printf(", swp_num = %d\n",swp_num);
	    for (k = 0; k < dim; ++k)
	    	(void) printf("\t\t%2d <= %s < %d\n",
	    		imin[k],var[idirs[k]],imax[k]);
	    (void) printf("\n");
	}
}		/*end set_sweep_limits*/

EXPORT	void	set_limits_for_open_bdry(
	Wave		*wv,
	Front		*fr,
	int		swp_num,
	int		*idirs,
	int		*imin,
	int		*imax)
{
	INTERFACE	*intfc;
	RECT_GRID	*rgr = wv->rect_grid;
	int		dir;

	if(rgr->dim != 3)
	    return;

	dir = idirs[0];
	intfc = fr->interf;
	if(rect_boundary_type(intfc,dir,0) == OPEN_BOUNDARY)
	    imin[0] = -rgr->lbuf[dir];
	if(rect_boundary_type(intfc,dir,1) == OPEN_BOUNDARY)
	    imax[0] = rgr->gmax[dir] + rgr->ubuf[dir];
}

EXPORT	int	*set_iperm(
	int step,
	int dim)
{
	int	    i, j;
	static	int iperm[3];
	static int perm3[6][3] =
	    {
	        {0, 1, 2},
		{2, 1, 0},
		{1, 2, 0},
		{0, 2, 1},
		{2, 0, 1},
		{1, 0, 2}
	    };
	static int perm6[48][6] = 
	    {
	        {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 5, 4}, {0, 1, 3, 2, 4, 5},
		{0, 1, 3, 2, 5, 4}, {0, 1, 4, 5, 2, 3}, {0, 1, 4, 5, 3, 2},
		{0, 1, 5, 4, 2, 3}, {0, 1, 5, 4, 3, 2}, {1, 0, 2, 3, 4, 5},
		{1, 0, 2, 3, 5, 4}, {1, 0, 3, 2, 4, 5}, {1, 0, 3, 2, 5, 4},
		{1, 0, 4, 5, 2, 3}, {1, 0, 4, 5, 3, 2}, {1, 0, 5, 4, 2, 3},
		{1, 0, 5, 4, 3, 2}, {2, 3, 0, 1, 4, 5}, {2, 3, 0, 1, 5, 4},
	        {2, 3, 1, 0, 4, 5}, {2, 3, 1, 0, 5, 4}, {2, 3, 4, 5, 0, 1},
		{2, 3, 4, 5, 1, 0}, {2, 3, 5, 4, 0, 1}, {2, 3, 5, 4, 1, 0},
	        {3, 2, 0, 1, 4, 5}, {3, 2, 0, 1, 5, 4}, {3, 2, 1, 0, 4, 5},
		{3, 2, 1, 0, 5, 4}, {3, 2, 4, 5, 0, 1}, {3, 2, 4, 5, 1, 0},
		{3, 2, 5, 4, 0, 1}, {3, 2, 5, 4, 1, 0}, {4, 5, 0, 1, 2, 3},
		{4, 5, 0, 1, 3, 2}, {4, 5, 1, 0, 2, 3}, {4, 5, 1, 0, 3, 2},
		{4, 5, 2, 3, 0, 1}, {4, 5, 2, 3, 1, 0}, {4, 5, 3, 2, 0, 1},
		{4, 5, 3, 2, 1, 0}, {5, 4, 0, 1, 2, 3}, {5, 4, 0, 1, 3, 2},
	        {5, 4, 1, 0, 2, 3}, {5, 4, 1, 0, 3, 2}, {5, 4, 2, 3, 0, 1},
		{5, 4, 2, 3, 1, 0}, {5, 4, 3, 2, 0, 1}, {5, 4, 3, 2, 1, 0}
	    };

	switch (dim)
	{
	case 1:
	    iperm[0] = 0;
	    iperm[1] = 1;
	    iperm[2] = 2;
	    break;
	case 2:
	    if (step%2)
	    {
	        iperm[0] = 1;
		iperm[1] = 0;
	    }
	    else
	    {
	        iperm[0] = 0;
		iperm[1] = 1;
	    }
	    iperm[2] = 2;
	    break;
	case 3:
	    i = (step/6)%48;
	    j = step%6;
            iperm[0] = perm3[perm6[i][j]][0];
            iperm[1] = perm3[perm6[i][j]][1];
            iperm[2] = perm3[perm6[i][j]][2];
	    break;
	}
	return iperm;
}		/*end set_iperm*/


/*
*			h_max_hyp_time_step():
*
*	Sets max_dt to the maximum time step allowed by the
*	Courant-Friedrichs-Levy condition for the appropriate
*	hyperbolic method.
*/

EXPORT double h_max_hyp_time_step(
	Wave		*wave,
	double		*coords)
{
	RECT_GRID	*gr = wave->rect_grid;
	double		max_dt;
	double		dt[MAXD];
	double		*h = gr->h;
	int		i, j, dim = gr->dim;

	max_dt = HUGE_VAL;
	if (wave->sizest == 0)
	    return max_dt;

	if  (strcmp(wave->method,"ADVANCE_FRONTS_ONLY") == 0)
	    return max_dt;

	for (i = 0; i < dim; ++i)
	{
	    if (Maxsp(wave)[i] > 0.0)
	    {
	    	dt[i] = h[i] / Maxsp(wave)[i];
	    	if (max_dt > dt[i])
	    	{
	    	    max_dt = dt[i];
	    	    for (j = 0; j < dim; ++j)
	    	        coords[j] = MaxWaveSpeedCoords(wave)[i][j];
	    	}
	    }
	    else
		dt[i] = HUGE_VAL;
	}
	if (debugging("time_step"))
	{
	    (void) printf("In h_max_hyp_time_step()\n");
	    for (i = 0; i < dim; ++i)
	    {
	    	(void) printf("wave: dt(%d) %g wavesp[%d] %g dx(%d) %g\n",
			      i,dt[i],i,Maxsp(wave)[i],i,h[i]);
		(void) printf("hyp: max_dt = %g\n",max_dt);
	    }
	    print_general_vector("coords = ",coords,dim,"\n");
	}

	/* Apr 10 2003: Myoung-Nyoun: must move to find_time_step() */
	{
	    double pcrds[MAXD];
	    double pdt;
	    pdt = h_max_parab_time_step(wave,pcrds);
	    if (max_dt > pdt)
	    {
		max_dt = pdt;
		for (j = 0; j < dim; j++)
		    coords[j] = pcrds[j];
	    }
	}

	return max_dt;
}			/*end h_max_hyp_time_step*/

/*
 * Apr 10 2003: Myoung-Nyoun: moved from for-loop of h_max_hyp_time_step
 *
 *     dt = dx * dx / ( 2 * visc ) : visc = kinetic viscosity
 *
 */
EXPORT double h_max_parab_time_step(
	Wave		*wave,
	double		*coords)
{
	RECT_GRID	*gr = wave->rect_grid;
	double		max_dt;
	double		dt[MAXD];
	double		*h = gr->h;
	int		i, j, dim = gr->dim;

	max_dt = HUGE_VAL;
	if (wave->sizest == 0)
	    return max_dt;

	if  (strcmp(wave->method,"ADVANCE_FRONTS_ONLY") == 0)
	    return max_dt;

	if(MaxViscosity(wave) && Maxvisc(wave) > 0.0)
	{
	    for (i = 0; i < dim; ++i)
	    {
		dt[i] = h[i]*h[i] / (2*Maxvisc(wave)); /*parabolic check*/ 
		if (max_dt > dt[i])
		{
		    max_dt = dt[i];
		    for (j = 0; j < dim; j++)
		        coords[j] = MaxViscosityCoords(wave)[j];
		}
	    }
	}
	else
	    for (i = 0; i < dim; ++i)
	        dt[i] = HUGE_VAL;

	if (debugging("parab_time_step"))
	{
	    (void) printf("In h_max_parab_time_step()\n");
	    for (i = 0; i < dim; ++i)
	    {
	    	(void) printf("wave: dt(%d)=%g, Maxvisc=%g, dx(%d)=%g\n",
			      i,dt[i],Maxvisc(wave),i,h[i]);
		(void) printf("parab: max_dt = %g\n",max_dt);
	    }
	    print_general_vector("coords = ",coords,dim,"\n");
	}
	return max_dt;
}			/*end h_max_parab_time_step*/

EXPORT	int reflect_pt_about_Nbdry(
	double		*coords,
	double		*coordsref,
	double		*nor,
	COMPONENT	int_comp,
	HYPER_SURF	*Nbdry,
	Front		*front)
{
	HYPER_SURF	   *hsbdry;
	HYPER_SURF_ELEMENT *hsebdry;
	double		   coordsbdry[MAXD];
	double		   ns[MAXD], ne[MAXD];
	double		   t[MAXD];
	int		   i, dim = front->rect_grid->dim;
	int		   icoords[MAXD];
	double 		   coords_tmp[MAXD];
	RECT_GRID 	   *gr = &topological_grid(front->interf);

	if (wave_type(Nbdry) != NEUMANN_BOUNDARY)
	    return NO;

	if (dim != 1)
	{
	    if (Nbdry->interface != front->interf)
	    {
	    	Nbdry = find_correspond_hyper_surface(Nbdry,NULL,NULL,
	    			                      front,front->interf);
	    	if (Nbdry == NULL)
		    return NO;
	    }
	    if (!Nbdry || Nbdry->interface != front->interf)
		return NO;
	    if (!rect_in_which(coords,icoords,gr))
	    {
		for (i = 0; i < dim; ++i)
		{
		    coords_tmp[i] = coords[i];
		    if (coords_tmp[i] < gr->L[i])
			coords_tmp[i] = gr->L[i];
		    else if (coords_tmp[i] > gr->U[i])
			coords_tmp[i] = gr->U[i];
		}
	    	if (!nearest_interface_point(coords_tmp,int_comp,
				front->interf,INCLUDE_BOUNDARIES,Nbdry,
				coordsbdry,t,&hsebdry,&hsbdry))
		    return NO;
	    }
	    else if (!nearest_interface_point(coords,int_comp,front->interf,
			                INCLUDE_BOUNDARIES,Nbdry,coordsbdry,t,
					&hsebdry,&hsbdry))
		return NO;
	}

	switch (dim)
	{
	case 1:
	    coordsbdry[0] = Coords(Point_of_hs(Nbdry))[0];
	    nor[0] = (coords[0] <
			0.5*(front->rect_grid->L[0]+front->rect_grid->U[0])) ?
			1.0 : -1.0;
	    break;
	case 2:
	    normal(Bond_of_hse(hsebdry)->start,hsebdry,hsbdry,ns,front);
	    normal(Bond_of_hse(hsebdry)->end,hsebdry,hsbdry,ne,front);
	    for (i = 0; i < dim; ++i)
	    	nor[i] = (1.0 - t[0])*ns[i] + t[0]*ne[i];
	    break;
	case 3:
	{
	    const double *tnor = Tri_normal(Tri_of_hse(hsebdry));
	    for (i = 0; i < dim; ++i)
	    	nor[i] = tnor[i];
	}
	    break;
	}

	for (i = 0; i < dim; ++i)
	    coordsref[i] = 2.0*coordsbdry[i] - coords[i];
	return YES;
}		/*end reflect_pt_about_Nbdry*/


EXPORT	Stencil	*alloc_stencil(
	int		npts,
	Front		*fr)
{
	Stencil		*sten;
	int		i, nrad = npts/2;

	scalar(&sten,sizeof(Stencil));
	bi_array(&sten->icoords_store,npts,3,INT);
	bi_array(&sten->crxstore,npts,MAX_NUM_CRX,sizeof(CRXING **));
	uni_array(&sten->ncstore,npts,INT);
	uni_array(&sten->worksp_store,npts,sizeof(Locstate));
	uni_array(&sten->worksp_st_store,npts,fr->sizest);
	uni_array(&sten->compstore,npts,sizeof(COMPONENT));
	uni_array(&sten->hsstore,npts,sizeof(HYPER_SURF *)); /*TODO REMOVE*/
	uni_array(&sten->pstore,npts,sizeof(POINT *));
	uni_array(&sten->ststore,npts,sizeof(Locstate));
	for (i = 0; i < npts; ++i)
	    sten->pstore[i] = Static_point(fr->interf);
	sten->npts = npts;
	sten->icoords = sten->icoords_store + nrad;
	sten->nc = sten->ncstore + nrad;
	sten->worksp = sten->worksp_store + nrad;
	sten->comp = sten->compstore + nrad;
	sten->crx = sten->crxstore + nrad;
	sten->hs = sten->hsstore + nrad; /*TODO REMOVE*/
	sten->p = sten->pstore + nrad;
	sten->st = sten->ststore + nrad;
	for (i = 0; i < npts; ++i)
	{
	    sten->worksp_store[i] = sten->worksp_st_store + i*fr->sizest;
	    sten->hsstore[i] = NULL;/*TODO REMOVE*/
	    sten->ncstore[i] = 0;
	    sten->crxstore[i][0] = NULL;
	}
	return sten;
}			/*end alloc_stencil*/

EXPORT	void	h_set_default_front_parameters(
	INIT_DATA	*init,
	Front		*fr)
{
	f_set_default_front_parameters(init,fr);
	fr->_copy_front =			h_copy_front;
	fr->_copy_into_front =			h_copy_into_front;
	fr->_print_Front_structure =		h_print_H_Front_structure;
#if defined(USE_OVERTURE)
        fr->_deep_copy_front =                  h_deep_copy_front;
#endif /* if defined(USE_OVERTURE) */
}		/*end h_set_default_front_parameters*/


/*
*                       h_deep_copy_front():
*
*       Basic default function for copying a front structure.
*       Allocates storage for the new front and copies the
*       argument into the new structure.
*/
#if defined(USE_OVERTURE)
LOCAL   Front *h_deep_copy_front(
        Front           *fr)
{
        Front           *newfr;

        scalar(&newfr,sizeof(H_Front));
        copy_into_front(newfr,fr);
        scalar(&(newfr->rect_grid), sizeof(RECT_GRID));
        *(newfr->rect_grid) = *(fr->rect_grid);
        scalar(&(newfr->pd_flag),sizeof(Patch_bdry_flag));
        return newfr;
}               /*end h_deep_copy_front*/
#endif /* if defined(USE_OVERTURE) */


/*
*			h_copy_front():
*
*	Basic default function for copying a front structure.
*	Allocates storage for the new front and copies the
*	argument into the new structure.
*/

LOCAL	Front *h_copy_front(
	Front		*fr)
{
	Front		*newfr;

	scalar(&newfr,sizeof(H_Front));
	copy_into_front(newfr,fr);
	return newfr;
}		/*end h_copy_front*/

/*
*			h_copy_into_front():
*
*	Copies fr into newfr.  Assumes newfr is already allocated.
*/

EXPORT	void h_copy_into_front(
	Front		*newfr,
	Front		*fr)
{
	f_copy_into_front(newfr,fr);
	HypFrontExtension(newfr) = HypFrontExtension(fr);
}		/*end h_copy_into_front*/

