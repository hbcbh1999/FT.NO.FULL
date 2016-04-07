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
*
*				gimksurf.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains functions for the construction of non-trivial
*	geometric curves needed for the initialization of an interface.
*	
*	The following exported functions are contained in this file.
*
*		make_random_surface()
*		g_make_ellipsoid()
*/

#if defined(THREED)

#include <ginit/ginit.h>

LOCAL	double   adjust_for_z_grid_spacing(double,double,double);
LOCAL	double	pert_height(double*,RECT_GRID*,SINE_PERT*,int);
LOCAL	double   spherical_harmonics(double,double,FOURIER_POLY*);
LOCAL 	double 	random_pert_func(POINTER,double*); 
LOCAL 	double 	bilin_interpolate(double,double,double,double,double,
				  double,double,double,double,double,
				  double*,double*,int*); 
LOCAL 	double 	file_amplitudes_pert_func(POINTER,double*);
LOCAL	void    perturb_ellipsoid(SURFACE*,ELLIPSOID*);

/*
*		make_random_surface():
*
*	Make_random_surface constructs the initial contact wave
*	for the "random surface" (statistics of fingers) problem.
*
*/

struct _RANDOM_PARAMS
{
	RECT_GRID *gr;
	SINE_PERT *pert;
};
typedef struct _RANDOM_PARAMS RANDOM_PARAMS;



EXPORT void make_random_surface(
	int w_type,
	Front *front,
	SINE_PERT *RS_b,
	COMPONENT compb,
	COMPONENT compa,
	double surf_ten,
	int dim)
{
	INTERFACE	   *intfc = front->interf;
	RECT_GRID          *gr = computational_grid(intfc);
	SINE_PERT	   *pert = Sine_pert(comp_type(compb));
	SURFACE		   *surf;
	RANDOM_PARAMS rand_params;

	rand_params.gr = gr;
	rand_params.pert = pert;

#if !defined(USE_OVERTURE)
	gr = computational_grid(intfc);
#else  /* if !defined(USE_OVERTURE) */
	{
	    RECT_GRID fine_grid, *tmpgr;
	    int i,refine;
            refine = 1;
	    tmpgr = computational_grid(intfc);
	    for(i = 1; i < front->NumberOfLevels; i++)
	    	refine *= 2;  /* should be adjustable */
	    for(i = 0; i < tmpgr->dim; i++)
	    {
	    	fine_grid.gmax[i] = refine*tmpgr->gmax[i];
		fine_grid.lbuf[i] = refine*tmpgr->lbuf[i];
		fine_grid.ubuf[i] = refine*tmpgr->ubuf[i];
	    }
	    set_rect_grid(tmpgr->L, tmpgr->U, tmpgr->GL, tmpgr->GU,
	    	fine_grid.lbuf, fine_grid.ubuf, fine_grid.gmax,
		tmpgr->dim, &tmpgr->Remap, &fine_grid);
	    gr = &fine_grid;
	}
#endif /* if !defined(USE_OVERTURE) */

	if (rand_params.pert->read_amplitudes_from_file == 0)
	    make_level_surface(gr,intfc,compb,compa,
	     file_amplitudes_pert_func,(POINTER)&rand_params,&surf);
	else
	    make_level_surface(gr,intfc,compb,compa,random_pert_func,
			(POINTER)&rand_params,&surf);

	wave_type(surf) = w_type;

	reset_normal_on_intfc(intfc);
	surface_tension(surf) = surf_ten;
	install_subdomain_bdry_curves(intfc);
	reset_intfc_num_points(intfc);
	
	(void)printf("\tEND make_random_surface\n");
}		/*end make_random_surface*/




LOCAL double random_pert_func(
	POINTER func_params,
	double *coords)
{
	RANDOM_PARAMS *rand_params = (RANDOM_PARAMS*)func_params;
	RECT_GRID *gr = rand_params->gr;
	SINE_PERT *pert = rand_params->pert;
	int dim = gr->dim;
	double dist;

	dist = coords[dim-1] - pert_height(coords,gr,pert,dim);
	return dist;
}	/* end random_pert_func */

/*
*		file_amplitudes_pert_func():
*
*	file_amplitudes_pert_func takes initial amplitudes read 
*       in from a file and uses bilinear interpolation to create 
* 	a function which passes through those data points. The 
*       initial amplitude data along with the corresponding 
*	x-y coordinate location of that data are stored in 
*	the SINE_PERT structure. This function is modeled after 
*	random_pert_function() so that it may be used in 
*	make_level_surface() and all relevant function 
*	calls nested within make_level_surface().
*/
LOCAL double file_amplitudes_pert_func(
	POINTER        func_params,
	double          *coords)
{
	RANDOM_PARAMS  *rand_params = (RANDOM_PARAMS*)func_params;
	RECT_GRID      *gr = rand_params->gr;
	SINE_PERT      *pert = rand_params->pert;
	int            dim = gr->dim;
	double          *GL = gr->GL;
	double          *GU = gr->GU;
	double          dist;
	double          crds[2];
	double          tolerance = 0.00000000001*min(gr->h[0],gr->h[1]);
	int            i,j,i0,j0;
	int            bdry_type[2] = {0,0};  


	for (i = 0; i < dim-1; i++) 
	{
	    crds[i] = coords[i];
	    if (pert->pert_bdry_type[i] == PERIODIC)
	    {
		while (crds[i] < GL[i])
		    crds[i] += (GU[i] - GL[i]);
		while (crds[i] > GU[i])
		    crds[i] -= (GU[i] - GL[i]);
	    }
	    else if (pert->pert_bdry_type[i] == SYMMETRIC)
	    {
		if (crds[i] < GL[i])
		    crds[i] = 2.0*GL[i] - coords[i];
		if (crds[i] > GU[i])
		    crds[i] = 2.0*GU[i] - coords[i];

		bdry_type[i] = 1;  
	    }
	}
	

	dist = crds[dim-1] - pert->z_intfc; /* default value for dist */

	/* Firstly, deal with the cases in which the function value 
	   (i.e. amplitude or z value) must be determined at an 
	   x-y coordinate location which is not interior to four 
	   x-y nodal locations for which initial amplitudes  
	   were read in from the initial interface amplitudes 
	   file - meaning "normal" bilinear interpolation cannot 
	   be used. So instead a bilinear interpolation based on 
	   extra "ghost" points determined from the extended 
	   domain (i.e., the x and y boundary conditions) is 
	   used. This case occupies the next two "main" if 
	   statements.

	   As of 11 March 2004, the following two main if loops deal 
	   only with the Periodic x and y boundary condition case. 
	   TODO - handle the symmetric boundary conditions case.
	*/
	if (crds[0] < pert->x_coord[0] ||
	    crds[0] > pert->x_coord[pert->Nx - 1])
	{
	    if((crds[1] >= pert->y_coord[0]) && 
		   (crds[1] <= pert->y_coord[pert->Ny - 1]))
	    {
	        for (j = 0; j < pert->Ny - 1; j++)
		{
         	    if ((crds[1] >= pert->y_coord[j]) && 
			(crds[1] <= pert->y_coord[j+1]))
		        break;
		}
		
		dist = coords[dim-1] - bilin_interpolate
		     (pert->x_coord[pert->Nx - 1], pert->x_coord[0], 
		     pert->y_coord[j], pert->y_coord[j+1],
		     (pert->amplitudes[pert->Nx - 1][j] + pert->z_intfc),
		     (pert->amplitudes[0][j] + pert->z_intfc),
		     (pert->amplitudes[0][j+1] + pert->z_intfc),
		     (pert->amplitudes[pert->Nx - 1][j+1] + pert->z_intfc),
		     crds[0],crds[1], GU,GL,bdry_type);
	    	return dist;
	    }
	    else if (crds[1] < pert->y_coord[0])
	    {
	        dist = coords[dim-1] - bilin_interpolate
		     (pert->x_coord[pert->Nx - 1], pert->x_coord[0], 
		     pert->y_coord[pert->Ny - 1], pert->y_coord[0],
		     (pert->amplitudes[pert->Nx - 1][pert->Ny - 1] 
		     + pert->z_intfc), (pert->amplitudes[0][pert->Ny - 1] 
		     + pert->z_intfc), (pert->amplitudes[0][0] + pert->z_intfc),
		     (pert->amplitudes[pert->Nx - 1][0] + pert->z_intfc),
		     crds[0],crds[1], GU,GL,bdry_type);
	    	return dist;
	    }
	    else if (crds[1] > pert->y_coord[pert->Ny - 1])
	    {
	        dist = coords[dim-1] - bilin_interpolate
		     (pert->x_coord[pert->Nx - 1], pert->x_coord[0], 
		     pert->y_coord[pert->Ny - 1], pert->y_coord[0],
		     (pert->amplitudes[pert->Nx - 1][pert->Ny - 1] 
		     + pert->z_intfc), (pert->amplitudes[0][pert->Ny - 1] 
		     + pert->z_intfc), (pert->amplitudes[0][0] + pert->z_intfc),
		     (pert->amplitudes[pert->Nx - 1][0] + pert->z_intfc),
		     crds[0],crds[1], GU,GL,bdry_type);
	    	return dist;
	    }
	}
       


	if ((crds[1] > pert->y_coord[pert->Nx - 1]) ||
	    (crds[1] < pert->y_coord[0]))
	{
	    if((crds[0] >= pert->x_coord[0]) && 
	       (crds[0] <= pert->x_coord[pert->Nx - 1]))
	    {
	        for (i = 0; i < pert->Nx - 1; i++)
		{
		    if ((crds[0] >= pert->x_coord[i]) && 
		        (crds[0] <= pert->x_coord[i+1]))
		      break;
		}

		dist = coords[dim-1] - bilin_interpolate
		     (pert->x_coord[i], pert->x_coord[i+1], 
		     pert->y_coord[pert->Ny - 1], pert->y_coord[0],
		     (pert->amplitudes[i][pert->Ny - 1] + pert->z_intfc),
		     (pert->amplitudes[i+1][pert->Ny - 1] + pert->z_intfc),
		     (pert->amplitudes[i+1][0] + pert->z_intfc),
		     (pert->amplitudes[i][0] + pert->z_intfc),
		     crds[0],crds[1], GU,GL,bdry_type);
	    	return dist;
	    }
	}



	/* Secondly, deal with the cases in which the x-y location 
	   at which the function value is to be determined 
	   corresponds to an x-y location for which data was 
	   read in from a file - so that no interpolation is 
	   done; the data value (amplitude) read in for that 
	   location is used.
	*/
	i0 = j0 = -1;
	for (i = 0; i < pert->Nx; i++)
	{
	    if (crds[0] - tolerance <= pert->x_coord[i] && 
		crds[0] + tolerance >= pert->x_coord[i]) 
	    {
	    	i0 = i;
		break;
	    }
	}
	for (j = 0; j < pert->Ny; j++)
	{
	    if (crds[1] - tolerance <= pert->y_coord[j] && 
		crds[1] + tolerance >= pert->y_coord[j]) 
	    {
	    	j0 = j;
		break;
	    }
	}
	if (i0 != -1 && j0 != -1)
	{
	    dist = coords[dim-1] - (pert->amplitudes[i][j] 
		                    + pert->z_intfc);
	    return dist;
	}

	/* Thirdly and finally, deal with the cases in which 
	   the function value must be determined at an x-y
	   coordinate location which is interior to four x-y
	   nodal locations for which initial amplitudes were
	   read in from the initial interface amplitudes 
	   file, BUT does not correspond to one of those 
	   nodal locations for which initial amplitudes exist 
	   - meaning "normal" bilinear interpolation is used.

	   For efficiency, the following loop could be 
	   incorporated into the preceding for loop. 1 Mar. 2004.
	*/
	for (i = 0; i < pert->Nx; i++)
	{
	    if (crds[0] > pert->x_coord[i] && 
	        crds[0] < pert->x_coord[i+1])
	    {
	    	i0 = i;
		break;
	    }
	}
	for(j = 0; j < pert->Ny; j++)
	{
	    if (crds[1] > pert->y_coord[j] && 
	        crds[1] < pert->y_coord[j+1])
	    {
	    	j0 = j;
		break;
	    }
	}
        dist = coords[dim-1] - bilin_interpolate
		     (pert->x_coord[i0], pert->x_coord[i0+1], 
		     pert->y_coord[j0], pert->y_coord[j0+1],
		     (pert->amplitudes[i0][j0] + pert->z_intfc),
		     (pert->amplitudes[i0+1][j0] + pert->z_intfc),
		     (pert->amplitudes[i0+1][j0+1] + pert->z_intfc),
		     (pert->amplitudes[i0][j0+1] + pert->z_intfc),
		     crds[0],crds[1], GU,GL,bdry_type);
	return dist;
}	/* end file_amplitudes_pert_func */

/* 
*                      bilin_interpolate():
*   Following bilinear interpolation for the z coordinate of the point in the 
*   grid square with corner points, counterclockwise from lower left, (x1,y1), 
*   (x2,y1), (x2,y2), and (x1,y2). Taken from "Numerical Recipes in C (2nd ed.)" 
*   pg. 123. 1 March 2004.
*/
LOCAL 	double 	bilin_interpolate(
		double x1,
		double x2,
		double y1,
		double y2,
		double f1,
		double f2,
		double f3,
		double f4,
		double x,
		double y,
		double *GU,
		double *GL,
		int *bdry_type) 
{
        double t,u, ans;
	double tol = 0.00000000001;



	if ((x1 < x2) && (y1 < y2)) /* usual case */
	{
	    t = (x - x1)/(x2 - x1);
	    u = (y - y1)/(y2 - y1);
	}
	/* Otherwise, if the points are from the edges of the 
	   computational domain such that some of the points 
	   are determined using the x or y boundary conditions
	   ("ghost points"), ...

	   TODO - symmetric boundary condition case AND
	   double check I've considered all possibilities in 
	   periodic boundary condition cases.
	*/
	else if ((x1 > x2) && (y1 < y2))
	{
	    if (bdry_type[0] == 0) /* periodic boundary */
	    {
	        if ((x + tol) < x1)
		    t = (x + (GU[0] - GL[0]) - x1)/
		        (x2 + (GU[0] - GL[0]) - x1);
		else /* if  ((x + tol) >= x1) */
		    t = (x - x1)/(x2 + (GU[0] - GL[0]) - x1);
		u = (y - y1)/(y2 - y1);
	    }
	}
	else if ((x1 < x2) && (y1 > y2))
	{
	    if (bdry_type[1] == 0) /* periodic boundary */
	    {
	        if (y + tol < y1)
		    u = (y + (GU[1] - GL[1]) - y1)/
		        (y2 + (GU[1] - GL[1]) - y1);
		else /* if  ((y + tol) >= y1) */
		    u = (y - y1)/(y2 + (GU[1] - GL[1]) - y1);
		t = (x - x1)/(x2 - x1);
	    }
	}
	else /* if ((x1 > x2) && (y1 > y2)) */
	{
	    if (bdry_type[0] == 0) /* periodic x boundary */
	    {
	        if ((x + tol) < x1)
		    t = (x + (GU[0] - GL[0]) - x1)/
		        (x2 + (GU[0] - GL[0]) - x1);
		else /* if  ((x + tol) >= x1) */
		    t = (x - x1)/(x2 + (GU[0] - GL[0]) - x1);
	    }
	    if (bdry_type[1] == 0) /* periodic y boundary */
	    {
	        if (y + tol < y1)
		    u = (y + (GU[1] - GL[1]) - y1)/
		        (y2 + (GU[1] - GL[1]) - y1);
		else /* if  ((y + tol) >= y1) */
		    u = (y - y1)/(y2 + (GU[1] - GL[1]) - y1);
		u = (y - y1)/(y2 + (GU[1] - GL[1]) - y1);
	    }
	}

	ans = (1.0-t)*(1.0-u)*f1 + t*(1.0-u)*f2 +t*u*f3 
	                         + (1.0-t)*u*f4;
	return ans;
}	/* end bilin_interpolate */

EXPORT HYPER_SURF       *g_make_ellipsoid(
	ELLIPSOID       *ellip,
	COMPONENT       compin,
	COMPONENT       compout,
	Front           *front)
{
	SURFACE         *s;
        RECT_GRID       *rgr = front->rect_grid;
	ELLIP_PARAMS	ep;
	int i;

	for (i = 0; i < rgr->dim; ++i)
	{
	    ep.cen[i] = ellip->cen[i];	
	    ep.rad[i] = ellip->rad[i];
	}
	make_level_surface(rgr,front->interf,compin,compout,ellipsoid_func,
			(POINTER)&ep,&s);
        ellip->hs = Hyper_surf(s);

        perturb_ellipsoid(s,ellip);
        reset_normal_on_intfc(front->interf);
        wave_type(s) = ellip->wv_type;
        surface_tension(s) = ellip->surf_tension;
        layer_index(ellip->hs) = ellip->layer_index;

        install_subdomain_bdry_curves(front->interf);
        reset_intfc_num_points(front->interf);
	untracked_hyper_surf(Hyper_surf(s)) = ellip->untracked;
        return Hyper_surf(s);
}		/* end g_make_ellipsoid */

EXPORT	void	g_make_ellip_region_boundaries3d(
	ELLIPSOID *ellip,
	Front     *front)
{
	boolean make_bdry;
	int   i, dim = front->rect_grid->dim;

	if (dim != 3)
	{
	    screen("ERROR in g_make_ellip_region_boundaries3d(), "
		   "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	}

	/* Check for the need to make region boundaries */
	make_bdry = NO;
	for (i = 0; i < 4; i++)
	{
	    if (ellip->btype[i] != UNKNOWN_BOUNDARY_TYPE)
	        make_bdry = YES;
	    if (ellip->obtype[i] != UNKNOWN_BOUNDARY_TYPE)
	        make_bdry = YES;
	}
	if (make_bdry == NO)
	    return;

	screen("ERROR in g_make_ellip_region_boundaries3d(), "
	       "function not implemented\n");
	clean_up(ERROR);
}		/*end g_make_ellip_region_boundaries3d*/

LOCAL	double	pert_height(
	double     *coords,
	RECT_GRID *gr,
	SINE_PERT *pert,
	int       dim)
{
	double *GL = gr->GL;
	double *GU = gr->GU;
	double h = gr->h[dim-1];
	double L = gr->VL[dim-1] + 0.5*h;
	int   i;
	double height;
	double crds[3];

	for (i = 0; i < dim-1; i++)
	{
	    crds[i] = coords[i];
	    if (pert->pert_bdry_type[i] == PERIODIC)
	    {
		while (crds[i] < GL[i])
		    crds[i] += (GU[i] - GL[i]);
		while (crds[i] > GU[i])
		    crds[i] -= (GU[i] - GL[i]);
	    }
	    else if (pert->pert_bdry_type[i] == SYMMETRIC)
	    {
		if (crds[i] < GL[i])
		    crds[i] = 2.0*GL[i] - coords[i];
		if (crds[i] > GU[i])
		    crds[i] = 2.0*GU[i] - coords[i];
	    }
	}
	height = pert_interface(pert,crds,0.0,dim);
	return adjust_for_z_grid_spacing(height,L,h);
}		/*end pert_height */

LOCAL	double adjust_for_z_grid_spacing(
	double height,
	double L,
	double h)
{
	double              k, delta;
	static const double htol = 0.004;/*TOLERANCE*/

	k = floor((height - L)/h);
	delta = (height - (L + k*h))/h;
	if (delta < htol)
	    height = L + (k + htol)*h;
	else if (1.0 - delta < htol)
	    height = L + (k + 1 - htol)*h;
	return height;
}	/* end adjust_for_z_grid_spacing */


LOCAL	void perturb_ellipsoid(
	SURFACE *s,
	ELLIPSOID *ellip)
{
	TRI          *tri;
	POINT        *p;
	int          i, j, k, l;
	double        phi, theta;
	double        x[3], xp[3], rr;
	double	     *cen = ellip->cen, **Q = ellip->Q;
	FOURIER_POLY *fpoly = ellip->fpoly;
	double	     Dr[3];
	double	     er[3];

	if (fpoly == NULL)
	    return;

	for (tri = first_tri(s); !at_end_of_tri_list(tri,s);
			tri = tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		sorted(p) = NO;
	    }
	}
	for (tri = first_tri(s); !at_end_of_tri_list(tri,s);
			tri = tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (sorted(p))
		    continue;
		for (k = 0; k < 3; k++)
		    x[k] = Coords(p)[k] - cen[k];
		for (k = 0; k < 3; k++)
		    for(xp[k] = 0, l = 0; l < 3; l++)
			xp[k] += Q[l][k]*x[l];
		rr = mag_vector(xp,3);
		for (k = 0; k < 3; k++)
		    er[k] = xp[k]/rr;
		theta = atan2(er[1],er[0]);
		phi = acos(er[2]);
		rr += spherical_harmonics(phi,theta,fpoly);
		for (i = 0; i < 3; i++)
		    Dr[i] = rr*er[i];
		for (i = 0; i < 3; i++)
		{
		    Coords(p)[i] = cen[i];
		    for (j = 0; j < 3; j++)
			Coords(p)[i] += Q[i][j]*Dr[j];
		}
		sorted(p) = YES;
	    }
	}
}	/* end perturb_ellipsoid */

LOCAL	double spherical_harmonics(
	double       phi,
	double       theta,
	FOURIER_POLY *fpoly)
{
	int 	k,l,m;
	int 	N = fpoly->num_modes;
	double   *phase = fpoly->phase;
	double   **nu = fpoly->nu;
	double   *A = fpoly->A;
	double 	z = 0.0;

	for (k = 0; k < N; k++)
	{
	    l = irint(nu[k][0]);
	    m = irint(nu[k][1]);
	    z += A[k]*SphericalHarmonic_s(l,m,phi,theta,phase[k]);
	}
	return z;
}	/* end spherical_harmonics */
#endif /* defined(THREED) */
