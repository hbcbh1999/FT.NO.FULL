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
*				gdriverstat.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	The following routines are accessed through function pointers set in
*	init_physics():
*
*	g_stat_var(), g_stat_flux() and
*	g_stat_inhom_source() are called by statistics() in dsub.c.
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <gdecs/gdecs.h>

		/* Possible types for return data from integrals */
		/* Stored in Item_Ctrl->type */

	/* SGas for computing integrals of state data */

typedef	struct {
	Locstate state; /* average state over cell (conserved quantities) */

	POINTER	 extra;  /* User supplied extra data */
	POINTER	 *avg;	/* Array of user supplied average structures */

	int	 navgs;

	/* Arrays of function pointers to find avg */
	void (**tri_comp_avg)(LINEAR_ELEMENT*,TRI_SOLN*,double,
			      Locstate,POINTER,int);
	void (**quad_comp_avg)(BILINEAR_ELEMENT*,TRI_SOLN*,double,
			       Locstate,POINTER,int);
	void (**line_tri_comp_avg)(LINEAR_ELEMENT*,TRI_SOLN*,double,Locstate,
				   POINTER,int,double*,double*,double,double*);
	void (**line_quad_comp_avg)(BILINEAR_ELEMENT*,TRI_SOLN*,double,Locstate,
				    POINTER,int,double*,double*,double,int*);
} SGas;

		/* Macros for accessing SGas data */

#define	Stat_num_avgs(ptr)		((SGas *) (ptr))->navgs

#define	compute_tri_average(ptr,i)				\
	(*(((SGas *) (ptr))->tri_comp_avg)[i])
#define	compute_line_tri_average(ptr,i)				\
	(*(((SGas *) (ptr))->line_tri_comp_avg)[i])
#define	compute_quad_average(ptr,i)				\
	(*(((SGas *) (ptr))->quad_comp_avg)[i])
#define	compute_line_quad_average(ptr,i)			\
	(*(((SGas *) (ptr))->line_quad_comp_avg)[i])


	/* LOCAL Function Declarations */
LOCAL	double	g_stat_flux(int,int,int*,Wave*);
LOCAL	double	g_stat_inhom_source(int,int*,Wave*,Front*,double);
LOCAL	double	g_stat_var(int,int*,Wave*,Front*,double);
LOCAL	void	avg_state(int*,Locstate,Front*,Wave*);
LOCAL	void	g_print_statistics(Grid*,Wave*,Front*,Printplot*,
				   OUTPUT_DATA*,boolean);


/*
*			DRIVER STATISTICS ROUTINES
*
*		g_stat_var():
*		g_stat_flux():
*		g_stat_inhom_source():
*
*	Routines for statistics.  These routines are called by statistics()
*	in a loop for "var" from 0 to grid->n_stat_var - 1.  n_stat_var is
*	set appropriately in init_physics.  The correspondence between 
*	var and the statistics quantity is:
*		0		density
*		1		energy
*		2 to 1+dim	momenta
*		2+dim		reaction product density	
*		3+dim		first constituent density
*/

/*
*			g_init_grid_statistics():
*
*	Initializes the physics dependent controls for grid statistics, then
*	calls d_init_grid_statistics() for the rest.
*/

EXPORT void g_init_grid_statistics(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Wave		*wave,
	Printplot	*prt)
{
	char		s[Gets_BUF_SIZE];
	int		nfloats = g_nfloats();
	static Grid_Stats_data *gs_data = NULL;

	if (nfloats == 0) return; /* no conserved variables for statistics */

	screen("Type 'y' to request grid statistics for conserved variables: ");
	(void) Gets(s);
	if ((s[0] != 'y') && (s[0] != 'Y')) 
	{
		prt->grid_statistics = NULL;
		prt->gs_data = NULL;
		return;
	}

	if (gs_data == NULL)
		scalar(&gs_data,sizeof(Grid_Stats_data));

	prt->gs_data = gs_data;

	gs_data->nfloats = nfloats;
	gs_data->stat_var = g_stat_var;
	gs_data->stat_flux = g_stat_flux;
	gs_data->stat_inhom_source = g_stat_inhom_source;
	gs_data->stat_point_flux = NULL;

	/* All other gs_data fields are set in the following. */

	d_init_grid_statistics(init,front,grid,wave,prt,g_print_statistics);

}		/*end g_init_grid_statistics*/


/*
*			g_print_statistics():
*
*	Provides a printout of the statistical quantities.
*	See the comments above the statistics routines
*	g_stat_... for the correspondence between the indices of
*	*present, *incremented, and *inhom_source_present and the
*	"conserved" quantities.
*	
*	ASSUMES gs_data->nfloats >= 4
*/

/*ARGSUSED*/
LOCAL void g_print_statistics(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*odata,
	boolean		about_to_stop)
{
	Grid_Stats_data *gs_data = GS_data(odata);
	FILE		*file = Output_file(odata);
	double		mass,enrgy,total_momentum,momentum;
	int		var;
	int		i, dim = grid->rect_grid->dim;
	long		nfloats = gs_data->nfloats;
	static double	*pres = NULL, *incr = NULL, *inhom_source = NULL;


	if (grid->pp_grid->nn > 1)
	{
	    if (pres == NULL)
	    {
	    	uni_array(&pres,        nfloats,FLOAT);
	    	uni_array(&incr,        nfloats,FLOAT);
	    	uni_array(&inhom_source,nfloats,FLOAT);
	    }
		
	    for (i = 0; i < nfloats; i++)
	    {
	    	pres[i]         = gs_data->present[i];
	    	incr[i]         = gs_data->incremented[i];
	    	inhom_source[i] = gs_data->inhom_source_present[i];
	    }
	    pp_global_sum(pres,        nfloats);
	    pp_global_sum(incr,        nfloats);
	    pp_global_sum(inhom_source,nfloats);
	}
	else
	{
	    pres         = gs_data->present;
	    incr         = gs_data->incremented;
	    inhom_source = gs_data->inhom_source_present;
	}

	(void) fprintf(file,"Grid Statistics Data for t = %g, step = %d:\n\n",
		       grid->time,grid->step);
	(void) fprintf(file,"\t\tpresent \tincremented\tsource  \t%% error\n");

	var = 0;
	mass = .5 * (pres[var] + incr[var]);
	(void) fprintf(file,"  Mass  \t%-10g\t%-10g\t%-10g\t %g\n",
		pres[var],incr[var],inhom_source[var],
		100. * (pres[var] - incr[var]) / mass);

	var = 1;
	enrgy = .5 * (pres[var] + incr[var]);
	(void) fprintf(file,"  Energy\t%-10g\t%-10g\t%-10g\t %g\n",
		pres[var],incr[var],inhom_source[var],
		100. * (pres[var] - incr[var]) / enrgy);

	total_momentum = sqrt(2. * mass * enrgy);

	for (i = 0; i < dim; i++)
	{
		var = i + 2;
		momentum = .5 * (fabs(pres[var]) + fabs(incr[var]));
		momentum = (momentum > .01 * total_momentum) ? 
			momentum : total_momentum;
		(void) fprintf(file,
			       "Momentum[%d]\t%-10g\t%-10g\t%-10g\t %g\n",
			       i,pres[var],incr[var],
		inhom_source[var],
		100. * (pres[var] - incr[var]) / momentum);
	}

	if (nfloats >= 3+dim)
	{
		var = 2+dim;
		mass = .5 * (pres[var] + incr[var]);
		(void) fprintf(file,"  Burned\t%-10g\t%-10g\t%-10g\t %g\n",
			pres[var],incr[var],inhom_source[var],
			100. * (pres[var] - incr[var]) / mass);
	}
	if (nfloats >= 4+dim)
	{
		var = 3+dim;
		mass = .5 * (pres[var] + incr[var]);
		(void) fprintf(file,"   Mass1\t%-10g\t%-10g\t%-10g\t %g\n",
			pres[var],incr[var],inhom_source[var],
			100. * (pres[var] - incr[var]) / mass);
	}

	if (gs_data->col_file != NULL)
	{
		(void) fprintf(gs_data->col_file,"%g ",grid->time);
		for (var = 0; var < nfloats; var++)
		{
			(void) fprintf(gs_data->col_file,"%g ",
				       gs_data->present[var]);
		}
		(void) fprintf(gs_data->col_file,"\n");
	}
}		/*end g_print_statistics*/


LOCAL double g_stat_var(
	int		var,
	int		*icoords,
	Wave		*wave,
	Front		*fr,double unit_area)
{
	static Locstate	state = NULL;
	int		dim = fr->interf->dim;
	
	if (state == NULL)
	    alloc_state(fr->interf,&state,fr->sizest);

		/* initialize state on first call with given icoords */
	if (var == 0)
	    avg_state(icoords,state,fr,wave);

	if (is_obstacle_state(state))
	    return 0.;

	switch (var)
	{
	case 0:
	    return Dens(state) * unit_area;
	case 1:
	    return Energy(state) * unit_area;
	}
	if (var <= dim + 1)
	    return Mom(state)[var-2] * unit_area;
#if defined(COMBUSTION_CODE)
	if (var == dim + 2)	
	{
	    if (Composition_type(state) == PTFLAME ||
	       (Composition_type(state) == THINFLAME))	
	    	return Burned(state) ? Dens(state) * unit_area : 0.;
	    else
	    	return Prod(state) * unit_area;
	}
	if (var == dim + 3)
	    return Dens1(state) * unit_area;
#endif /* defined(COMBUSTION_CODE) */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(state)->n_comps != 1)
            {
                return pdens(state)[var-(dim+2)]*unit_area;
            }
            return ERROR_FLOAT;
        }
	screen("ERROR: unknown variable var in g_stat_var()\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end g_stat_var*/


LOCAL double g_stat_flux(
	int		dir,
	int		var,
	int		*icoords,
	Wave		*wave)
{
	Locstate	state;
	int		dim = wave->rect_grid->dim;

	state = Rect_state(icoords,wave);

	if (is_obstacle_state(state)) return 0.;

	switch(var)	
	{
	case 0:
		return Mom(state)[dir];
	case 1:
		return vel(dir,state) * (Energy(state) + pressure(state));
	}
	if (var < dim + 2)
		return Mom(state)[dir] * Mom(state)[var-2] / Dens(state)
			+ (var - 2 == dir) ? pressure(state) : 0.;
#if defined(COMBUSTION_CODE)
	if (var == dim + 2)
	{
		if ((Composition_type(state) == PTFLAME) ||
		    (Composition_type(state) == THINFLAME))	
			return Burned(state) ? Mom(state)[dir] : 0.;
		else
			return Prod(state) * vel(dir,state);
	}
	if (var == dim + 3)
		return Dens1(state) * vel(dir,state);
#endif /* defined(COMBUSTION_CODE) */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(state)->n_comps != 1)
            {
                return pdens(state)[var-(dim+2)]*vel(dir,state);
            }
            return ERROR_FLOAT;
        }
	screen("ERROR: unknown variable var in g_stat_flux_x()\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end g_stat_flux*/


/*
*		g_stat_inhom_source():
*
*	Calculates the contribution of the inhomogeneous source terms 
*	introduced by gravity and 3-D cylindrical symmetry. 
*/

LOCAL double g_stat_inhom_source(
	int		var,
	int		*icoords,
	Wave		*wave,
	Front		*fr,
	double		unit_area)
{
	static Locstate	state = NULL;
	const double  *g;
	double  *coords;
	static double	alpha;
	double r, rmin;
	double		source;
	int		dim = fr->interf->dim;
	
	if (state == NULL)
	{
	    alloc_state(fr->interf,&state,fr->sizest);
	    if (is_rotational_symmetry())
		alpha = rotational_symmetry();
	}

	if (var == 0)
	{
	    avg_state(icoords,state,fr,wave);
	}
	if (is_obstacle_state(state))
	    return 0.;

	coords = Rect_coords(icoords,wave);
	g = gravity(coords,fr->time);
	if (is_rotational_symmetry())
	{
	    r = pos_radius(coords[0],fr->rect_grid);
	    rmin = pos_radius(0.0,fr->rect_grid);
	}

	/* TODO: add code for gravity (for CYLINDRICAL codes) */

	switch (var)
	{
	case 0:
	    source = 0.0;
	    if (is_rotational_symmetry() &&
		(alpha > 0.0) && (fabs(r) > fabs(rmin)))
		source = -alpha*unit_area*Mom(state)[0]/r;
	    return source;
	case 1:
	    source = (is_gravity() == YES) ?
		scalar_product(Mom(state),g,dim)*unit_area : 0.0;
	    if (is_rotational_symmetry() &&
		alpha > 0.0)
		source -= alpha*unit_area*(Energy(state)+pressure(state))*
				           Mom(state)[0]/(r*Dens(state));
	    return source;
	}
	if (var < dim + 2)
	{
	    source = 0.;
	    if (is_gravity() == YES)
		source = Dens(state)*g[var-2]*unit_area;
	    if (is_rotational_symmetry() &&
		(alpha > 0.0) && (fabs(r) > fabs(rmin)))
	    {
	    	source -= alpha*unit_area*Mom(state)[0]*Mom(state)[var-2] /
			      (r*Dens(state));
	    }
	    return source;
	}
	if (var < dim + 4)
	    return 0.0;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            return 0.0;
        }
	screen("ERROR: unknown variable var in g_stat_inhom_source()\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end g_stat_inhom_source*/



/*
*			avg_state():
*
*	Computes the state for the (dual grid) square icoords based on the 
*	fraction of the square in each connected component defined by the front.
*	Used in statistics routines above.
*
*	TODO: Currently just returns the state corresponding to the
*	      (regular grid) block icoords.  Should use integral().
*	      When implemented correctly, changes will be needed
*	      in the loop indices for the present[] and inhom_source[]
*	      calculations in statistics().  Changes may also be needed
*	      in g_stat_var() and g_stat_inhom_source().
*/

LOCAL void avg_state(
	int		*icoords,	/* dual mesh indices */
	Locstate	answer,	       	/* average state in mesh block */
	Front		*fr,
	Wave		*wave)
{
	int		i, dim = fr->interf->dim, test = 0;

		/* TODO: CODE NEEDED */

	for (i = 0; i < dim; i++)
		if (icoords[i] < fr->rect_grid->gmax[i]) test++;
	if (test == dim)
	{
		ft_assign(answer,Rect_state(icoords,wave),fr->sizest);
	}
	else
	{
		obstacle_state(fr->interf,answer,fr->sizest);
	}
}		/*end avg_state*/

/*
*			g_tri_integral():
*
*	Calculates the integral over a triangle.
*	and adds the quantity to the state in ans.
*/

EXPORT double g_tri_integral(
	LINEAR_ELEMENT	*et,
	TRI_SOLN	*soln,
	Locstate	ans,
	POINTER		user)
{
	double		area;
	double		tri_area;
	Locstate	*s = et->s;
	RECT_GRID	*gr = &soln->tri_grid->comp_grid;
	int		dim = gr->dim;
	int		nsts = dim+1;
	int		i, j;

	if (coord_system() != RECTANGULAR)
	{
	    screen("ERROR in g_tri_integral(), function not implemented\n");
	    clean_up(ERROR);
	}

	area = measure_of_linear_element(et,soln);
	tri_area = area/nsts;

	Set_params(ans,s[0]);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	for (j = 0; j < nsts; j++)
	{
		Dens(ans) += tri_area*Dens(s[j]);
		Energy(ans) += tri_area*Energy(s[j]);
		for (i = 0; i < dim; i++) Mom(ans)[i] += tri_area*Mom(s[j])[i];
#if defined(COMBUSTION_CODE)
		if (Composition_type(ans) == ZND)
			Prod(ans) += tri_area*Prod(s[j]);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	}

	if (user != NULL)
	{
		int i;
		int n = Stat_num_avgs(user);

		for (i = 0; i < n; i++)
		    compute_tri_average(user,i)(et,soln,tri_area,ans,user,i);
	}
	return area;
}		/*end g_tri_integral*/


/*
*			g_quad_integral():
*
*	Calculates the integral over a quadrangle.  It assumes area = 1.
*/

EXPORT double g_quad_integral(
	BILINEAR_ELEMENT *eq,
	TRI_SOLN	*soln,
	Locstate	ans,
	POINTER		user)
{
	Locstate	s[(1<<MAXD)];
	double		area;
	double		nfac;
	double		*fans, *fs;
	int		dim = soln->tri_grid->comp_grid.dim;
	int		nsts = 1 << dim;
	int		i, j;
	size_t		sizest = size_of_state(soln->tri_grid->grid_intfc);
	int		nfloats = g_nfloats();

	if (coord_system() != RECTANGULAR)
	{
	    screen("ERROR in g_quad_integral(), function not implemented\n");
	    clean_up(ERROR);
	}

	states_on_bilinear_element(s,eq,soln->tri_grid);
	area = 1.0;

	clear_state(soln->tri_grid->grid_intfc,ans,sizest);
	Set_params(ans,s[0]);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	fans = (double *)ans;
	for (j = 0; j < nsts; j++)
	{
	    fs = (double *)s[j];
	    for (i = 0; i < nfloats; i++)
	    	fans[i] += fs[i];
	}
	nfac = 1.0/nsts;
	for (i = 0; i < nfloats; i++)
		fans[i] *= nfac;

	if (user != NULL)
	{
	    int i;
	    int n = Stat_num_avgs(user);

	    for (i = 0; i < n; i++)
	    {
	        compute_quad_average(user,i)(eq,soln,area,ans,user,i);
	    }
	}
	return area;
}		/*end g_quad_integral*/
