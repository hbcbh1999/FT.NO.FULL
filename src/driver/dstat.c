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
*				dstat.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains statistics computation driver routines.
*/


#include <driver/ddecs.h>

	/* LOCAL function prototypes*/
LOCAL	void	accumulate_sources(Grid_Stats_data*,Front*,Grid*,Wave*);
LOCAL	void	d_reg_grid_cyl_statistics1d(Grid_Stats_data*,
	                                    Grid*,Wave*,Front*,int);
LOCAL	void	d_reg_grid_ident_statistics1d(Grid_Stats_data*,
	                                      Grid*,Wave*,Front*,int);
LOCAL	void	d_reg_grid_ident_statistics3d(Grid_Stats_data*,
	                                      Grid*,Wave*,Front*,int);
LOCAL	void	end_of_statistics(Grid*);
LOCAL	void	initialize_for_statistics(Grid_Stats_data*,Grid*,int);

/*
*			d_init_grid_statistics():
*
*	Initializes the physics independent controls for the grid statistics.
*	Should be called from a physics function, which sets all physics
*	dependent data except for the print_function().
*
*	TODO:	Use should be made of states on the front,
*		especially for the incremented variables.
*		The notion of "outside" boundary has to be
*		generalized, since (e.g.) there is a pressure
*		contribution to the line integral along a
*		wedge in gas dynamics.  Perhaps prev/next_bdry
*		should be used.
*/

EXPORT	void d_init_grid_statistics(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Wave		*wave,
	Printplot	*prt,
	void		(*print_function)(Grid*,Wave*,Front*,Printplot*,
	                                  OUTPUT_DATA*,boolean))
{
	Grid_Stats_data *gs_data = prt->gs_data;
	RECT_GRID	*gr = grid->rect_grid;
	char		s[Gets_BUF_SIZE];
	int		var, source;
	int		nfloats = gs_data->nfloats;
	int		num_point_sources = wave->num_point_sources;
	int		dim = grid->rect_grid->dim;
	double		*initial_present, *present, *inhom_source_present;
	double		*incremented, **source_incremented;

	if (nfloats == 0) return; /* no conserved variables for statistics */

	if (prt->grid_statistics == NULL)
	{
	    switch (gr->Remap.remap)
	    {
	    case IDENTITY_REMAP:
	        switch (dim)
	        {
	        case 1:
	            prt->grid_statistics = d_reg_grid_ident_statistics1d;
	            break;
	        case 2:
	            prt->grid_statistics = d_reg_grid_ident_statistics2d;
	            break;
	        case 3:
	            prt->grid_statistics = d_reg_grid_ident_statistics3d;
	            break;
	        }
	        break;
	    case CYLINDRICAL_REMAP:
	        switch (dim)
	        {
	        case 1:
	            prt->grid_statistics = d_reg_grid_cyl_statistics1d;
	            break;
	        case 2:
	            prt->grid_statistics = d_reg_grid_cyl_statistics2d;
	            break;
	        }
	        break;
	    case SPHERICAL_REMAP:
	        prt->grid_statistics = NULL;/* TODO */
	        return;
	    }
	}

	init_output_data(init,&gs_data->data,grid,prt,
			 "the grid stats data",YES,YES,NO);
	add_user_output_function(print_function,&gs_data->data,prt);

	screen("Request printout of conserved variable statistical data\n");
	screen("\tin a separate output file in columnar form (y,n(dflt): ");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	{
	    gs_data->col_file = open_data_file(front,"the columnar data",
	                                       YES,NO,NULL,&gs_data->col_dname,
	                                       NULL,&gs_data->col_filename);
	}
	else
	{
	    gs_data->col_dname = NULL;
	    gs_data->col_filename = NULL;
	    gs_data->col_file = NULL;
	}

	gs_data->num_point_sources = num_point_sources;

	gs_data->initial_present = NULL;
	gs_data->present = NULL;
	gs_data->incremented = NULL;
	gs_data->inhom_source_present = NULL;
	gs_data->source_incremented = NULL;

	uni_array(&gs_data->initial_present,nfloats,FLOAT);
	uni_array(&gs_data->present,        nfloats,FLOAT);
	uni_array(&gs_data->incremented,    nfloats,FLOAT);

	if (gs_data->stat_inhom_source != NULL)
	{
	    uni_array(&gs_data->inhom_source_present,nfloats,FLOAT);
	    inhom_source_present = gs_data->inhom_source_present;
	    for (var = 0;  var < nfloats;  var++)
	    {
	        inhom_source_present[var] = 0.0;
	    }
	}

	if (num_point_sources > 0)
	{
	    bi_array(&gs_data->source_incremented,
	           num_point_sources,nfloats,FLOAT);
	    source_incremented = gs_data->source_incremented;
	}

	initial_present = gs_data->initial_present;
	present = gs_data->present;
	incremented = gs_data->incremented;

	for (var = 0;  var < nfloats;  var++)
	{
	    present[var] = 0.0;
	    incremented[var] = 0.0;
	    initial_present[var] = 0.0;
	    for (source = 0;  source < num_point_sources; source++)
	        source_incremented[source][var] = 0.0;
	}
	screen("\n");
	return;
}		/*end d_init_grid_statistics*/

/*ARGSUSED*/
LOCAL	void d_reg_grid_ident_statistics1d(
	Grid_Stats_data	*gs_data,
	Grid		*grid,
	Wave		*wave,
	Front		*fr,
	int		acc_all)
{
#if defined(ONED)

	static const char *fname = "d_reg_grid_ident_statistics1d";

	RECT_GRID	*rect_grid = grid->rect_grid;
	double	(*Area)(const double*,const RECT_GRID*) = rect_grid->Remap.Area;
	double	*present = gs_data->present;
	double	*incremented = gs_data->incremented;
	double	*initial_present = gs_data->initial_present;
	double	*inhom_source_present = gs_data->inhom_source_present;
	double	(*stat_var)(int,int*,Wave*,Front*,double) = gs_data->stat_var;
	double	(*flux)(int,int,int*,Wave*) = gs_data->stat_flux;
	double	(*inhom_source)(int,int*,Wave*,Front*,double) =
	                        gs_data->stat_inhom_source;
	int	nfloats = gs_data->nfloats;
	double	*xx_grid;
	int	*gmax = rect_grid->gmax;
	int	ix, ixmax = gmax[0];
	double	coords[MAXD];
	double	unit_area;
	double	dt = grid->dt;
	int	var;
	int	icoords[MAXD];
	static	boolean	first = YES;

	debug_print("statistics","Entered %s\n",fname);
	if (debugging("statistics"))
	{
	    (*wave->show_wave_states)(wave);
	    graph_front_states(fr);
	}

	initialize_for_statistics(gs_data,grid,acc_all);
	xx_grid = rect_grid->edges[0];

	if (inhom_source)	/* Coded thusly for speed */
	{
	    for (ix = 0; ix < ixmax; ix++)
	    {
	        coords[0] = xx_grid[ix];
	        icoords[0] = ix;
	        unit_area = (*Area)(coords,rect_grid);
	        for (var = 0; var < nfloats; var++)
	        {
	            if (first || acc_all)
	            {
	                present[var] += (*stat_var)(var,icoords,wave,
	                                            fr,unit_area);
	            }
	            inhom_source_present[var] += (*inhom_source)(var,icoords,
	                                                         wave,fr,
	                                                         unit_area);
	        }
	    }
	}
	else
	{
	    if (first || acc_all)
	    {
	        for (ix = 0; ix < ixmax; ix++)
	        {
	            coords[0] = xx_grid[ix];
	            icoords[0] = ix;
	            unit_area = (*Area)(coords,rect_grid);
	            for (var = 0; var < nfloats; var++)
	            {
	                present[var] += (*stat_var)(var,icoords,
	                                            wave,fr,unit_area);
	            }
	        }
	    }
	}

	if (first) /* Initialize incremented quantities */
	{
	    first = NO;
	    for (var = 0; var < nfloats; var++)
	    {
	        incremented[var] = present[var];
	        initial_present[var] = present[var];
	    }
	}

	            /* Update Incremented Quantities */

	icoords[0] = 0;
	for (var = 0; var < nfloats; var++)
	    incremented[var] += dt * (*flux)(0,var,icoords,wave);
	icoords[0] = ixmax-1;
	for (var = 0; var < nfloats; var++)
	    incremented[var] -= dt * (*flux)(0,var,icoords,wave);

	    /* include inhomogeneous source in incremented quantities */

	accumulate_sources(gs_data, fr, grid, wave);

	end_of_statistics(grid);

	debug_print("statistics","Left %s\n",fname);

#endif /* defined(ONED) */
}		/*end d_reg_grid_ident_statistics1d*/

/*ARGSUSED*/
EXPORT	void d_reg_grid_ident_statistics2d(
	Grid_Stats_data	*gs_data,
	Grid		*grid,
	Wave		*wave,
	Front		*fr,
	int		acc_all)
{
#if defined(TWOD)

	static const char *fname = "d_reg_grid_ident_statistics2d";

	RECT_GRID	*rect_grid = grid->rect_grid;
	double	(*Area)(const double*,const RECT_GRID*) = rect_grid->Remap.Area;
	double	*present = gs_data->present;
	double	*incremented = gs_data->incremented;
	double	*initial_present = gs_data->initial_present;
	double	*inhom_source_present = gs_data->inhom_source_present;
	double	(*stat_var)(int,int*,Wave*,Front*,double) = gs_data->stat_var;
	double	(*flux)(int,int,int*,Wave*) = gs_data->stat_flux;
	double	(*inhom_source)(int,int*,Wave*,Front*,double) =
	                    gs_data->stat_inhom_source;
	int	nfloats = gs_data->nfloats;
	double	*h = rect_grid->h;
	int	*gmax = rect_grid->gmax;
	double	dt = grid->dt;
	double	dxdt = h[0]*dt, *xx_grid;
	double	dydt = h[1]*dt, *yy_grid;
	int	ix, ixmax = gmax[0];
	int	iy, iymax = gmax[1];
	double	coords[MAXD];
	double	unit_area;
	int	var;
	int	icoords[MAXD];
	static	boolean	first = YES;

	debug_print("statistics","Entered %s\n",fname);
	if (debugging("statistics"))
	{
	    (*wave->show_wave_states)(wave);
	    graph_front_states(fr);
	}

	initialize_for_statistics(gs_data,grid,acc_all);
	xx_grid = rect_grid->edges[0];
	yy_grid = rect_grid->edges[1];

	if (inhom_source)	/* Coded thusly for speed */
	{
	    for (iy = 0; iy < iymax; iy++)
	    {
	        coords[1] = yy_grid[iy];
	        icoords[1] = iy;
	        for (ix = 0; ix < ixmax; ix++)
	        {
	            coords[0] = xx_grid[ix];
	            icoords[0] = ix;
	            unit_area = (*Area)(coords,rect_grid);
	            for (var = 0; var < nfloats; var++)
	            {
	                if (first || acc_all)
	                {
	                    present[var] += (*stat_var)(var,icoords,wave,
							fr,unit_area);
	                }
	                inhom_source_present[var] +=
	                        (*inhom_source)(var,icoords,wave,fr,unit_area);
	            }
	        }
	    }
	}
	else
	{
	    if (first || acc_all)
	    {
	        for (iy = 0; iy < iymax; iy++)
	        {
	            coords[1] = yy_grid[iy];
	            icoords[1] = iy;
	            for (ix = 0; ix < ixmax; ix++)
	            {
	                coords[0] = xx_grid[ix];
	                icoords[0] = ix;
	                unit_area = (*Area)(coords,rect_grid);
	                for (var = 0; var < nfloats; var++)
	                {
	                    present[var] += (*stat_var)(var,icoords,
	                                                wave,fr,unit_area);
	                }
	            }
	        }
	    }
	}

	if (first) /* Initialize incremented quantities */
	{
	    first = NO;
	    for (var = 0; var < nfloats; var++)
	    {
	        incremented[var] = present[var];
	        initial_present[var] = present[var];
	    }
	}

	        /* Update Incremented Quantities */

	for (ix = 0; ix < ixmax; ix++)
	{ 
	    icoords[0] = ix;
	    icoords[1] = 0;
	    for (var = 0; var < nfloats; var++)
	        incremented[var] += dxdt * (*flux)(1,var,icoords,wave);
	    icoords[1] = iymax-1;
	    for (var = 0; var < nfloats; var++)
	        incremented[var] -= dxdt * (*flux)(1,var,icoords,wave);
	}
	for (iy = 0; iy < iymax; iy++)
	{ 
	    icoords[1] = iy;
	    icoords[0] = 0;
	    for (var = 0; var < nfloats; var++)
	        incremented[var] += dydt * (*flux)(0,var,icoords,wave);
	    icoords[0] = ixmax-1;
	    for (var = 0; var < nfloats; var++)
	        incremented[var] -= dydt * (*flux)(0,var,icoords,wave);
	}

	    /* include inhomogeneous source in incremented quantities */

	accumulate_sources(gs_data, fr, grid, wave);

	end_of_statistics(grid);

	debug_print("statistics","Left %s\n",fname);

#endif /* defined(TWOD) */
}		/*end d_reg_grid_ident_statistics2d*/


/*ARGSUSED*/
LOCAL	void d_reg_grid_ident_statistics3d(
	Grid_Stats_data	*gs_data,
	Grid		*grid,
	Wave		*wave,
	Front		*fr,
	int		acc_all)
{
#if defined(THREED)

	static const char *fname = "d_reg_grid_ident_statistics3d";

	RECT_GRID	*rect_grid = grid->rect_grid;
	double	(*Area)(const double*,const RECT_GRID*) = rect_grid->Remap.Area;
	double	*present = gs_data->present;
	double	*incremented = gs_data->incremented;
	double	*initial_present = gs_data->initial_present;
	double	*inhom_source_present = gs_data->inhom_source_present;
	double	(*stat_var)(int,int*,Wave*,Front*,double) = gs_data->stat_var;
	double	(*flux)(int,int,int*,Wave*) = gs_data->stat_flux;
	double	(*inhom_source)(int,int*,Wave*,Front*,double) =
	                        gs_data->stat_inhom_source;
	int	nfloats = gs_data->nfloats;
	static	boolean	first = YES;

	double	*h = rect_grid->h;
	double	dt = grid->dt;
	int	*gmax = rect_grid->gmax;
	int	ix, ixmax = gmax[0];
	int	iy, iymax = gmax[1];
	int	iz, izmax = gmax[2];
	double	*xx_grid, *yy_grid, *zz_grid;
	double	dxdydt = h[0]*h[1]*dt;
	double	dxdzdt = h[0]*h[2]*dt;
	double	dydzdt = h[1]*h[2]*dt;
	double	coords[MAXD];
	double	unit_area;
	int	var;
	int	icoords[MAXD];

	debug_print("statistics","Entered %s\n",fname);
	if (debugging("statistics"))
	{
	    (*wave->show_wave_states)(wave);
	    graph_front_states(fr);
	}

	initialize_for_statistics(gs_data,grid,acc_all);
	xx_grid = rect_grid->edges[0];
	yy_grid = rect_grid->edges[1];
	zz_grid = rect_grid->edges[2];

	if (inhom_source)	/* Coded thusly for speed */
	{
	    for (iz = 0; iz < izmax; iz++)
	    {
	        coords[2] = zz_grid[iz];
	        icoords[2] = iz;
	        for (iy = 0; iy < iymax; iy++)
	        {
	            coords[1] = yy_grid[iy];
	            icoords[1] = iy;
	            for (ix = 0; ix < ixmax; ix++)
	            {
	                coords[0] = xx_grid[ix];
	                icoords[0] = ix;
	                unit_area = (*Area)(coords,rect_grid);
	                for (var = 0; var < nfloats; var++)
	                {
	                    if (first || acc_all)
	                    {
	                        present[var] += (*stat_var)(var,icoords,
	                                                    wave,fr,unit_area);
	                    }
	                    inhom_source_present[var] +=
	                                (*inhom_source)(var,icoords,
	                                                wave,fr,unit_area);
	                }
	            }
	        }
	    }
	}
	else
	{
	    for (iz = 0; iz < izmax; iz++)
	    {
	        coords[2] = zz_grid[iz];
	        icoords[2] = iz;
	        if (first || acc_all)
	        {
	            for (iy = 0; iy < iymax; iy++)
	            {
	                coords[1] = yy_grid[iy];
	                icoords[1] = iy;
	                for (ix = 0; ix < ixmax; ix++)
	                {
	                    coords[0] = xx_grid[ix];
	                    icoords[0] = ix;
	                    unit_area = (*Area)(coords,rect_grid);
	                    for (var = 0; var < nfloats; var++)
	                    {
	                        present[var] += (*stat_var)(var,icoords,
	                                                    wave,fr,unit_area);
	                    }
	                }
	            }
	        }
	    }
	}

	if (first) /* Initialize incremented quantities */
	{
	    first = NO;
	    for (var = 0; var < nfloats; var++)
	    {
	        incremented[var] = present[var];
	        initial_present[var] = present[var];
	    }
	}

	        /* Update Incremented Quantities */

	for (iz = 0; iz < izmax; iz++)
	{
	    icoords[2] = iz;
	    for (ix = 0; ix < ixmax; ix++)
	    { 
	        icoords[0] = ix;
	        icoords[1] = 0;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] += dxdzdt * (*flux)(1,var,icoords,wave);
	        icoords[1] = iymax-1;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] -= dxdzdt * (*flux)(1,var,icoords,wave);
	    }
	}
	for (iz = 0; iz < izmax; iz++)
	{
	    icoords[2] = iz;
	    for (iy = 0; iy < iymax; iy++)
	    { 
	        icoords[1] = iy;
	        icoords[0] = 0;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] += dydzdt * (*flux)(0,var,icoords,wave);
	        icoords[0] = ixmax-1;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] -= dydzdt * (*flux)(0,var,icoords,wave);
	    }
	}
	for (iy = 0; iy < iymax; iy++)
	{
	    icoords[1] = iy;
	    for (ix = 0; ix < ixmax; ix++)
	    { 
	        icoords[0] = ix;
	        icoords[2] = 0;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] += dxdydt * (*flux)(2,var,icoords,wave);
	        icoords[0] = izmax-1;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] -= dxdydt * (*flux)(2,var,icoords,wave);
	    }
	}

	    /* include inhomogeneous source in incremented quantities */

	accumulate_sources(gs_data, fr, grid, wave);

	end_of_statistics(grid);

	debug_print("statistics","Left %s\n",fname);

#endif /* defined(THREED) */
}		/*end d_reg_grid_ident_statistics3d*/


/*ARGSUSED*/
LOCAL	void d_reg_grid_cyl_statistics1d(
	Grid_Stats_data	*gs_data,
	Grid		*grid,
	Wave		*wave,
	Front		*fr,
	int		acc_all)
{
#if defined(ONED)

	static const char *fname = "d_reg_grid_cyl_statistics1d";

	RECT_GRID	*rect_grid = grid->rect_grid;
	double	(*Area)(const double*,const RECT_GRID*) = rect_grid->Remap.Area;
	double	*present = gs_data->present;
	double	*incremented = gs_data->incremented;
	double	*initial_present = gs_data->initial_present;
	double	*inhom_source_present = gs_data->inhom_source_present;
	double	(*stat_var)(int,int*,Wave*,Front*,double) = gs_data->stat_var;
	double	(*flux)(int,int,int*,Wave*) = gs_data->stat_flux;
	double	(*inhom_source)(int,int*,Wave*,Front*,double) =
	                                gs_data->stat_inhom_source;
	int	nfloats = gs_data->nfloats;
	int	*gmax = rect_grid->gmax;
	double	*L = rect_grid->L, *U = rect_grid->U;
	int	ix, ixmax = gmax[0];
	double	*xx_grid;
	double	coords[MAXD];
	double	unit_area;
	double	dt = grid->dt;
	int	var;
	int	icoords[MAXD];
	int	icoords1[MAXD];
	static	boolean	first = YES;

	debug_print("statistics","Entered %s\n",fname);
	if (debugging("statistics"))
	{
	    (*wave->show_wave_states)(wave);
	    graph_front_states(fr);
	}

	initialize_for_statistics(gs_data,grid,acc_all);
	xx_grid = rect_grid->edges[0];

	if (inhom_source)	/* Coded thusly for speed */
	{
	    for (ix = 0; ix < ixmax; ix++)
	    {
	        coords[0] = xx_grid[ix];
	        icoords[0] = ix;
	        unit_area = (*Area)(coords,rect_grid);
	        for (var = 0; var < nfloats; var++)
	        {
	            if (first || acc_all)
	            {
	                present[var] += (*stat_var)(var,icoords,
	                                            wave,fr,unit_area);
	            }
	            inhom_source_present[var] += (*inhom_source)(var,icoords,
	                                                         wave,fr,
	                                                         unit_area);
	        }
	    }
	}
	else
	{
	    if (first || acc_all)
	    {
	        for (ix = 0; ix < ixmax; ix++)
	        {
	            coords[0] = xx_grid[ix];
	            icoords[0] = ix;
	            unit_area = (*Area)(coords,rect_grid);
	            for (var = 0; var < nfloats; var++)
	            {
	                present[var] += (*stat_var)(var,icoords,
	                                            wave,fr,unit_area);
	            }
	        }
	    }
	}

	if (first) /* Initialize incremented quantities */
	{
	    first = NO;
	    for (var = 0; var < nfloats; var++)
	    {
	        incremented[var] = present[var];
	        initial_present[var] = present[var];
	    }
	}

	        /* Update Incremented Quantities */

	icoords[0] = 0;	icoords1[0] = ixmax-1;
	for (var = 0; var < nfloats; var++)
	    incremented[var] += dt * L[0] * (*flux)(0,var,icoords,wave);
	for (var = 0; var < nfloats; var++)
	    incremented[var] -= dt * U[0] * (*flux)(0,var,icoords1,wave);

	    /* include inhomogeneous source in incremented quantities */

	accumulate_sources(gs_data, fr, grid, wave);

	end_of_statistics(grid);

	debug_print("statistics","Left %s\n",fname);

#endif /* defined(ONED) */
}		/*end d_reg_grid_cyl_statistics1d*/

/*ARGSUSED*/
EXPORT	void d_reg_grid_cyl_statistics2d(
	Grid_Stats_data	*gs_data,
	Grid		*grid,
	Wave		*wave,
	Front		*fr,
	int		acc_all)
{
#if defined(TWOD)

	static const char *fname = "d_reg_grid_cyl_statistics2d";

	RECT_GRID	*rect_grid = grid->rect_grid;
	double	(*Area)(const double*,const RECT_GRID*) = rect_grid->Remap.Area;
	double	*present = gs_data->present;
	double	*incremented = gs_data->incremented;
	double	*initial_present = gs_data->initial_present;
	double	*inhom_source_present = gs_data->inhom_source_present;
	double	(*stat_var)(int,int*,Wave*,Front*,double) = gs_data->stat_var;
	double	(*flux)(int,int,int*,Wave*) = gs_data->stat_flux;
	double	(*inhom_source)(int,int*,Wave*,Front*,double) =
	                                gs_data->stat_inhom_source;
	int	nfloats = gs_data->nfloats;
	double	*h = rect_grid->h, *L = rect_grid->L, *U = rect_grid->U;
	double	dt = grid->dt;
	int	*gmax = rect_grid->gmax;
	int	ix, ixmax = gmax[0];
	int	iy, iymax = gmax[1];
	double	x, dxdt	= h[0]*dt, *xx_grid;
	double	   dydt = h[1]*dt, *yy_grid;
	double	coords[MAXD];
	double	unit_area;
	int	var;
	int	icoords[MAXD];
	int	icoords1[MAXD];
	static	boolean	first = YES;

	debug_print("statistics","Entered %s\n",fname);
	if (debugging("statistics"))
	{
	    (*wave->show_wave_states)(wave);
	    graph_front_states(fr);
	}

	initialize_for_statistics(gs_data,grid,acc_all);
	xx_grid = rect_grid->edges[0];
	yy_grid = rect_grid->edges[1];

	if (inhom_source)	/* Coded thusly for speed */
	{
	    for (iy = 0; iy < iymax; iy++)
	    {
	        coords[1] = yy_grid[iy];
	        icoords[1] = iy;
	        for (ix = 0; ix < ixmax; ix++)
	        {
	            coords[0] = xx_grid[ix];
	            icoords[0] = ix;
	            unit_area = (*Area)(coords,rect_grid);
	            for (var = 0; var < nfloats; var++)
	            {
	                if (first || acc_all)
	                {
	                    present[var] += (*stat_var)(var,icoords,
	                                                wave,fr,unit_area);
	                }
	                inhom_source_present[var] += (*inhom_source)(var,
	                                                             icoords,
	                                                             wave,fr,
	                                                             unit_area);
	            }
	        }
	    }
	}
	else
	{
	    if (first || acc_all)
	    {
	        for (iy = 0; iy < iymax; iy++)
	        {
	            coords[1] = yy_grid[iy];
	            icoords[1] = iy;
	            for (ix = 0; ix < ixmax; ix++)
	            {
	                coords[0] = xx_grid[ix];
	                icoords[0] = ix;
	                unit_area = (*Area)(coords,rect_grid);
	                for (var = 0; var < nfloats; var++)
	                {
	                    present[var] += (*stat_var)(var,icoords,
	                                                wave,fr,unit_area);
	                }
	            }
	        }
	    }
	}

	if (first) /* Initialize incremented quantities */
	{
	    first = NO;
	    for (var = 0; var < nfloats; var++)
	    {
	        incremented[var] = present[var];
	        initial_present[var] = present[var];
	    }
	}

	        /* Update Incremented Quantities */

	icoords[0] = 0;	icoords1[0] = ixmax-1;
	for (iy = 0; iy < iymax; iy++)
	{ 
	    icoords[1] = iy;
	    icoords1[1] = iy;
	    for (var = 0; var < nfloats; var++)
	        incremented[var] += dydt * L[0] * (*flux)(0,var,icoords,wave);
	    for (var = 0; var < nfloats; var++)
	        incremented[var] -= dydt * U[0] * (*flux)(0,var,icoords1,wave);
	}


	/* TODO: Check validity of this step for other versions */
	icoords[1] = 0;	icoords1[1] = iymax-1;
	for (ix = 0; ix < ixmax; ix++)
	{ 

	    x = xx_grid[ix];
	    icoords[0] = ix;
	    icoords1[0] = ix;
	    for (var = 0; var < nfloats; var++)
	    {
	        incremented[var] += x * dxdt * (   (*flux)(1,var,icoords,wave)
	                               - (*flux)(1,var,icoords1,wave)  );
	    }
	}

	    /* include inhomogeneous source in incremented quantities */

	accumulate_sources(gs_data, fr, grid, wave);

	end_of_statistics(grid);

	debug_print("statistics","Left %s\n",fname);

#endif /* defined(TWOD) */
}		/*end d_reg_grid_cyl_statistics2d*/


/*ARGSUSED*/
EXPORT	void d_dual_grid_ident_statistics3d(
	Grid_Stats_data	*gs_data,
	Grid		*grid,
	Wave		*wave,
	Front		*fr,
	int		acc_all)
{
#if defined(THREED)

	static const char *fname = "d_dual_grid_ident_statistics3d";

	RECT_GRID	*rect_grid = grid->rect_grid;
	double	(*Area)(const double*,const RECT_GRID*) = rect_grid->Remap.Area;
	double	*present = gs_data->present;
	double	*incremented = gs_data->incremented;
	double	*initial_present = gs_data->initial_present;
	double	*inhom_source_present = gs_data->inhom_source_present;
	double	(*stat_var)(int,int*,Wave*,Front*,double) = gs_data->stat_var;
	double	(*flux)(int,int,int*,Wave*) = gs_data->stat_flux;
	double	(*inhom_source)(int,int*,Wave*,Front*,double) =
	                                gs_data->stat_inhom_source;
	int	nfloats = gs_data->nfloats;
	double	*h = rect_grid->h;
	double	dt = grid->dt;
	int	*gmax = rect_grid->gmax;
	int	ix, ixmax = gmax[0];
	int	iy, iymax = gmax[1];
	int	iz, izmax = gmax[2];
	double	*dxx_grid, *dyy_grid, *dzz_grid;
	double	dxdydt = h[0]*h[1]*dt;
	double	dxdzdt = h[0]*h[2]*dt;
	double	dydzdt = h[1]*h[2]*dt;
	double	coords[MAXD];
	double	unit_area;
	int	var;
	int	icoords[MAXD];
	static	boolean	first = YES;

	debug_print("statistics","Entered %s\n",fname);
	if (debugging("statistics"))
	{
	    (*wave->show_wave_states)(wave);
	    graph_front_states(fr);
	}

	initialize_for_statistics(gs_data,grid,acc_all);
	dxx_grid = rect_grid->centers[0];
	dyy_grid = rect_grid->centers[1];
	dzz_grid = rect_grid->centers[2];

	if (inhom_source)	/* Coded thusly for speed */
	{
	    for (iz = 0; iz <= izmax; iz++)
	    {
	        coords[2] = dzz_grid[iz];
	        icoords[2] = iz;
	        for (iy = 0; iy <= iymax; iy++)
	        {
	            coords[1] = dyy_grid[iy];
	            icoords[1] = iy;
	            for (ix = 0; ix <= ixmax; ix++)
	            {
	                coords[0] = dxx_grid[ix];
	                icoords[0] = ix;
	                unit_area = (*Area)(coords,rect_grid);
	                for (var = 0; var < nfloats; var++)
	                {
	                    if (first || acc_all)
	                    {
	                        present[var] += (*stat_var)(var,icoords,
	                                                    wave,fr,unit_area);
	                    }
	                    inhom_source_present[var] +=
	                        (*inhom_source)(var,icoords,wave,fr,unit_area);
	                }
	            }
	        }
	    }
	}
	else
	{
	    if (first || acc_all)
	    {
	        for (iz = 0; iz <= izmax; iz++)
	        {
	            coords[2] = dzz_grid[iz];
	            icoords[2] = iz;
	            for (iy = 0; iy <= iymax; iy++)
	            {
	                coords[1] = dyy_grid[iy];
	                icoords[1] = iy;
	                for (ix = 0; ix <= ixmax; ix++)
	                {
	                    coords[0] = dxx_grid[ix];
	                    icoords[0] = ix;
	                    unit_area = (*Area)(coords,rect_grid);
	                    for (var = 0; var < nfloats; var++)
	                    {
	                        present[var] += (*stat_var)(var,icoords,wave,
	                                                    fr,unit_area);
	                    }
	                }
	            }
	        }
	    }
	}

	if (first) /* Initialize incremented quantities */
	{
	    first = NO;
	    for (var = 0; var < nfloats; var++)
	    {
	        initial_present[var] = present[var];
	    }
	}

	        /* Update Incremented Quantities */

	for (iz = 0; iz < izmax; iz++)
	{
	    icoords[2] = iz;
	    for (ix = 0; ix < ixmax; ix++)
	    { 
	        icoords[0] = ix;
	        icoords[1] = 0;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] += dxdzdt * (*flux)(1,var,icoords,wave);
	        icoords[1] = iymax-1;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] -= dxdzdt * (*flux)(1,var,icoords,wave);
	    }
	}
	for (iz = 0; iz < izmax; iz++)
	{
	    icoords[2] = iz;
	    for (iy = 0; iy < iymax; iy++)
	    { 
	        icoords[1] = iy;
	        icoords[0] = 0;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] += dydzdt * (*flux)(0,var,icoords,wave);
	        icoords[0] = ixmax-1;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] -= dydzdt * (*flux)(0,var,icoords,wave);
	    }
	}
	for (iy = 0; iy < iymax; iy++)
	{
	    icoords[1] = iy;
	    for (ix = 0; ix < ixmax; ix++)
	    { 
	        icoords[0] = ix;
	        icoords[2] = 0;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] += dxdydt * (*flux)(2,var,icoords,wave);
	        icoords[0] = izmax-1;
	        for (var = 0; var < nfloats; var++)
	            incremented[var] -= dxdydt * (*flux)(2,var,icoords,wave);
	    }
	}

	    /* include inhomogeneous source in incremented quantities */

	accumulate_sources(gs_data, fr, grid, wave);

	end_of_statistics(grid);

	debug_print("statistics","Left %s\n",fname);

#endif /* defined(THREED) */
}		/*end d_dual_grid_ident_statistics3d*/

LOCAL	void initialize_for_statistics(
	Grid_Stats_data	*gs_data,
	Grid		*grid,
	int		acc_all)
{
	int	var;
	int	nfloats = gs_data->nfloats;
	double	*present = gs_data->present;
	double	*inhom_source_present = gs_data->inhom_source_present;

	/* Compute conserved quantities and inhomogeneous source present */

	if (gs_data->stat_inhom_source)	/* Coded as if-else for speed */
	{
	    for (var = 0; var < nfloats; var++)
	    {
	        if (acc_all) present[var] = 0.0;
	        inhom_source_present[var] = 0.0;
	    }
	}
	else if (acc_all)
	{
	    for (var = 0; var < nfloats; var++) present[var] = 0.0;
	}

	    /* The calculation of present[] now uses the hyp_solution
	    * function integral() which operates on the dual grid    
	    * from L[i] - 0.5 * cell_width(0,i,grid->rect_grid) ->
	    *      U[i] + 0.5 * cell_width(gmax[i]-1,i,grid->rect_grid)
	    * and is gmax[i]+1 grid block long in the ith direction. */ 

	if (!set_grid_lines(grid->rect_grid))
	{
	    screen("ERROR in initialize_for_statistics(), "
	           "set_grid_lines() failed\n");
	    clean_up(ERROR);
	}
}		/*end initialize_for_statistics*/


LOCAL	void accumulate_sources(
	Grid_Stats_data	*gs_data,
	Front		*fr,
	Grid		*grid,
	Wave		*wave)
{
	RECT_GRID	*rect_grid = grid->rect_grid;
	double		*incremented = gs_data->incremented;
	double		*inhom_source_present = gs_data->inhom_source_present;
	double		**source_incremented = gs_data->source_incremented;
	double		(*pf)(int,int,Wave*) = gs_data->stat_point_flux;
	double		*h = rect_grid->h, *L = rect_grid->L, *U = rect_grid->U;
	double		dt = grid->dt;
	double		geom_factor;
	double		temp;
	int		j, var, source;
	int		nfloats = gs_data->nfloats;
	int		dim = rect_grid->dim;

	if (gs_data->stat_inhom_source)
	{
	    for (var = 0; var < nfloats; var++)
	        incremented[var] += dt * inhom_source_present[var];
	}

	    /* include point source contributions */

	    /* The numerical constant in geom_factor accounts for the   */
	    /* geometrical fraction of point source output flowing into */
	    /* the computational square. If a well is on the boundary   */
	    /* IT IS ASSUMED that half of the flow is into the region   */
	    /* of computation. For a corner well this means 1/4 of its  */
	    /* flow.						    */

	for (source = 0; source < wave->num_point_sources; source++)
	{
	    static const double EDGE_TOL = 0.0001;	/*TOLERANCE*/
	    geom_factor = dt;
	    for (j = 0; j < dim; j++)
	    {
	        if ((wave->srcpt[j][source] <= L[j] + EDGE_TOL * h[j])
	            || (wave->srcpt[j][source] >= U[j] + EDGE_TOL * h[j]))
	            geom_factor *= 0.5;
	    }
	    if (debugging("statistics"))
	    {
	        (void) printf("well: ");
	        for (j = 0; j < dim; j++)
	            (void) printf("%g ",wave->srcpt[j][source]);
	        (void) printf("dt %g geom_factor %g\n",dt,geom_factor);
	    }

	    for (var = 0; var < nfloats; var++)
	    {
	        temp = geom_factor * (*pf)(var,source,wave);
	        incremented[var] += temp;
	        source_incremented[source][var] += temp;
	    }
	    (*wave->stat_prod_rate)(source,wave,fr,wave->prod_rate);
	}
}		/*end accumulate_sources*/

LOCAL	void end_of_statistics(Grid* grid)
{
	free_grid_lines(grid->rect_grid);
	if (debugging("Grid"))
	{
	    (void) printf("Grid structure after statistics computation\n");
	    print_Grid_structure(grid);
	}
}		/*end end_of_statistics*/
