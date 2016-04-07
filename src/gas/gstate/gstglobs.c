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
*				gstglobs.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains global variables for the use in the gas dynamics code.
*	All globals in this file are local to the file, but can be
*	evaluated or set by calls to the functions in this file.
*
*
*	TODO:  Introduce a physics structure to pass these parameters
*
*	Contains
*
*		g_sizest()
*		set_coord_sys()
*		g_set_sizeof_state()
*		set_composition_type()
*		g_obstacle_state()
*		return_obst_state()
*/

#include <gdecs/gdecs.h>
#if defined(HYPRE)
#include <SN_ELLIP.h>
#include <VectorMatrix.h>
#endif /* defined(HYPRE) */

typedef struct {
	AVISC	Avisc;	       /* Global defaults for artifical parameters */
	double	alpha;	       /* geometry factor */
	GRAVITY	*grav_data;	  /* gravitational force,  NULL if no gravity */
	int	composition_type; /* generic type of materials in flow */
	size_t	sizest;	       /* size of flow state structure */
	int	nfloats;       /* number of floating point variables in state*/
	int	n_comps;       /* number of possible components in a state*/
	GEOMETRY geometry;     /* geometry flag for coordinate system */
	int	dim;	       /* number of space dimensions in the flow */
	Gas_param **params_list;/* List of EOS params used in simulation */
#if defined(ONED)
	int	state_type;
#endif /* defined(ONED) */
} FLOW_PARAMS;

	/* LOCAL function declarations */
LOCAL const double *astrophysical_gravity(const double*,GRAVITY*);
LOCAL const double *radial_gravity(const double*,GRAVITY*);
LOCAL const double *time_dependent_gravity(const double,GRAVITY*);
LOCAL const double *generalized_astrophysical_gravity(const double*,GRAVITY*);
LOCAL boolean eval_partial_mass1D(double*,int,int);
LOCAL boolean eval_partial_mass2D(double*,int,int);
LOCAL boolean eval_partial_mass3D(double*,int,int);
#if defined(HYPRE)
LOCAL void solve_gravity();
LOCAL void interpolate_gravity(const double*,double*,double*,int);
#endif /* defined(HYPRE) */


LOCAL	FLOW_PARAMS Fparams = {
	{
	    NULL,  /* hyp_method */
	    NO,    /* use_lin_av */
	    NO,	   /* use_lapidus_av */
	    NO,	   /* use_upwind_av */
	    YES,   /* use_msf */
	    0.0,   /* linear_visc_coef */
	    0.0,   /* SEE NOTE!!! lapidus_visc_coef */
	    0.0,   /* upwind_visc_coef */
	    2.0,   /* msf_ieta */
	    0.25,  /* min_shock_jump */
	    1e-06, /* min_sp_vol_jump */
	    0.0,   /* heat_cond */
	    1.0,   /* SEE NOTE!!! sp_coef */
	    0.0,   /* char_speed_cutoff */
	    0.0    /* dynamic surface tension */
	}, /* Avisc */
	0.0,			/* alpha */
	NULL,			/* grav_data */
	PURE_NON_REACTIVE,	/* composition_type */
	sizeof(Gas),		/* sizest */
	MAXD+2,			/* nfloats */
	1,			/* n_comps */
	RECTANGULAR,		/* geometry */
	MAXD,			/* dim */
	NULL,			/*params_list */
#if defined(ONED)
	GAS_STATE,		/* state_type */
#endif /* defined(ONED) */
};

EXPORT size_t g_sizest(void)
{
	return Fparams.sizest;
}		/*end g_sizest*/

EXPORT int g_nfloats(void)
{
	return Fparams.nfloats;
}		/*end g_nfloats*/

EXPORT	void g_set_sizeof_state(
	CHART		*chart,
	size_t		sizest,
	int		nfloats)
{
	Fparams.sizest = sizest;
	Fparams.nfloats = nfloats;
	if (chart != NULL)
	{
	    chart->wave->sizest = sizest;
	    chart->wave->nfloats = nfloats;
	    chart->front->sizest = sizest;
	    chart->front->nfloats = nfloats;
	}
}		/*end g_set_sizeof_state*/

EXPORT void set_composition_type(
	int		type)
{
	Fparams.composition_type = type;
}		/*end set_composition_type*/

EXPORT int g_composition_type(void)
{
	return Fparams.composition_type;
}		/*end g_composition_type*/


EXPORT void g_obstacle_state(
	Locstate	state,
	size_t		sizest)
{
	if (sizest != 0)
	{
	    zero_scalar(state,sizest);
	    Dens(state) = MACH_EPS;/*Allow divide by rho*/
	    Params(state) = NULL;
	    set_type_of_state(state,OBSTACLE_STATE);
	}
}		/*end g_obstacle_state*/

EXPORT Locstate return_obst_state(void)
{
	static	Locstate obst_state = NULL;

	if (obst_state == NULL)
	{
	    scalar(&obst_state,Fparams.sizest);
	    g_obstacle_state(obst_state,Fparams.sizest);
	}

	return obst_state;
}		/*end return_obst_state*/


EXPORT void set_gravity(
	GRAVITY		*grav_data)
{
	F_USER_INTERFACE *fuh;
	struct Table    *T;
	int		i;

	Fparams.grav_data = grav_data;
	if (grav_data == NULL || grav_data->type == NO_GRAVITY)
	    return;
	
	for (i = 0; i < MAXD; i++)
	{
	    fuh = f_user_hook(i+1);
#if !defined(COMBUSTION_CODE)
	    fuh->_reflect_state = g_reflect_and_stratify_state;
#endif /* defined(COMBUSTION_CODE) */
	}
	for (T = interface_table_list(); T != NULL; T = T->next)
	{
	    if (interface_type(T->interface) != PHYSICAL_INTERFACE)
	    	continue;
#if !defined(COMBUSTION_CODE)
	    f_user_interface(T->interface)._reflect_state =
	    	g_reflect_and_stratify_state;
#endif /* defined(COMBUSTION_CODE) */
	}
}		/*end set_gravity*/

EXPORT	boolean is_gravity(void)
{
	if (Fparams.grav_data == NULL)
	    return NO;
	else if (Fparams.grav_data->type == NO_GRAVITY)
	    return NO;
	else
	    return YES;
}		/*end is_gravity*/

EXPORT	double astrophys_grav_constant()
{
	return Fparams.grav_data->G;
}	/* end astrophys_grav_constant */

EXPORT const double *gravity(
	const double *coords,
	const double time)
{
	static double no_g[3] = {0.0, 0.0, 0.0};
	if ((Fparams.grav_data==NULL) || (Fparams.grav_data->type==NO_GRAVITY))
	{
	    return (const double*) no_g;
	}
	switch (Fparams.grav_data->type)
	{
	case NO_GRAVITY:
	case CONSTANT_GRAVITY:
	    return (const double*) Fparams.grav_data->g;
	case TIME_DEPENDENT_GRAVITY:
	    return time_dependent_gravity(time,Fparams.grav_data);
	case ASTROPHYSICAL_GRAVITY:
	    return astrophysical_gravity(coords,Fparams.grav_data);
	case RADIAL_GRAVITY:
	    return radial_gravity(coords,Fparams.grav_data);
	case GENERALIZED_ASTROPHYSICAL_GRAVITY:
	    return generalized_astrophysical_gravity(coords,Fparams.grav_data);
	case USER_DEFINED:
	    return user_defined_gravity(coords,time,Fparams.grav_data);
	default:
	    screen("ERROR in gravity(), invalid gravity type %d\n",
		   Fparams.grav_data->type);
	    clean_up(ERROR);
	}
	return (const double*) no_g;
}		/*end gravity*/

EXPORT	void	eval_gravity(
	const double *coords,
	const double time,
	double	    *grav)
{
	const double *g = gravity(coords,time);
	int   i, dim = Fparams.dim;
	for (i = 0; i < dim; i++)
	    grav[i] = g[i];
}		/*end eval_gravity*/

LOCAL const double *time_dependent_gravity(
	const double time,
	GRAVITY	    *grav_data)
{
	double	*g = grav_data->g;
	double	**g_of_t;
	double	tau;
	int	i, dim, num_time_points;
	static	double last_time;
	static	int last_n;
	static	boolean first = YES;

	if (first == YES)
	{
	    last_time = -HUGE_VAL;
	    last_n = 0;
	}
	if (time == last_time)
	    return (const double*) g;
	last_time = time;
	dim = grav_data->dim;
	num_time_points = grav_data->num_time_points;
	g_of_t = grav_data->g_of_t;

	if (time < g_of_t[0][0])
	{
	    last_n = 0;
	    for (i = 0; i < dim; i++)
		g[i] = g_of_t[0][i+1];
	    return (const double*) g;
	}
	if (time > g_of_t[num_time_points-1][0])
	{
	    last_n = num_time_points-1;
	    for (i = 0; i < dim; i++)
		g[i] = g_of_t[num_time_points-1][i+1];
	    return (const double*) g;
	}
	for (i = last_n; i < (num_time_points-1); i++)
	    if ((g_of_t[i][0] <= time) && (time < g_of_t[i+1][0]))
		break;
	if (i == (num_time_points-1))
	{
	    for (i = 0; i < (num_time_points-1); i++)
	        if ((g_of_t[i][0] <= time) && (time < g_of_t[i+1][0]))
		break;
	}
	if (i == (num_time_points-1))
	{
	    screen("ERROR in time_dependent_gravity(), can't find time %g\n",
		   time);
	}
	last_n = i;
	tau = (time - g_of_t[last_n][0]) /
	      (g_of_t[last_n+1][0] - g_of_t[last_n][0]);
	for (i = 0; i < dim; i++)
	{
	    g[i] = (1.0 - tau)*g_of_t[last_n][i+1] + tau*g_of_t[last_n+1][i+1];
	}
	return (const double*) g;
}		/*end time_dependent_gravity*/

LOCAL const double *astrophysical_gravity(
	const double *coords,
	GRAVITY     *grav_data)
{
	double	*g = grav_data->g;
	double	*cen = grav_data->center;
	double	G = grav_data->G;
	double	M = grav_data->M;
	int	i, dim = grav_data->dim;
	double	r;

	for (i = 0; i < dim; i++)
	    g[i] = cen[i] - coords[i];
	r = mag_vector(g,dim);
	for (i = 0; i < dim; i++)
	    g[i] = G*M*g[i]/(r*r*r);
	return (const double*)g;
}		/*end astrophysical_gravity*/

LOCAL const double *generalized_astrophysical_gravity(
	const double *coords,
	GRAVITY     *grav_data)
{
	double	*g = grav_data->g;
	double	*cen = grav_data->center;
	double	G = grav_data->G;
	RECT_GRID *gr = grav_data->roots[0]->wave->rect_grid;
	int	i, dim = grav_data->dim;
	double	r, M, R, dR;
	double   *h = gr->h;
	GEOMETRY_REMAP remap = gr->Remap.remap;
        double   sign[MAXD], tmpcoords[MAXD];

#if defined(HYPRE)
	if(debugging("2d_gravity"))
	{

        for (i = 0; i < dim; i++)
        {
            tmpcoords[i] = coords[i];
            if (coords[i] < 0)
                sign[i] = -1;
            else
                sign[i] = 1;
        }
        for (i = 0; i < dim; i++) tmpcoords[i] *= sign[i];
        interpolate_gravity((const double *)tmpcoords,g,h,dim);
        for (i = 0; i < dim; i++) g[i] *= -sign[i];
        return (const double*)g;

	}
#endif /* defined(HYPRE) */

	if(grav_data->Mp==NULL)
	{
	    for (i = 0; i < dim; i++)
		g[i] = 0;
	    return (const double*)g;
	}

	R = grav_data->R;
	dR = grav_data->dR;

	for (i = 0; i < dim; i++)
	{
	    g[i] = cen[i] - coords[i];
	}

	r = mag_vector(g,dim);
	i = (int)(r/dR);
	if(i > grav_data->N-1)
	{
	    M = grav_data->Mp[grav_data->N-1];
	}
	else if (i == 0)
	{
            M = grav_data->Mp[0]*r/dR;
	}
	else
	{  
	    double  shell_mass;
	    shell_mass = grav_data->Mp[i]-grav_data->Mp[i-1];
	    M = grav_data->Mp[i-1] + (r - i*dR)/dR*shell_mass;
	}

	switch (dim)
	{
	case 1:
	    switch (remap)
	    {
	    case IDENTITY_REMAP:
		g[0] = G*M;
	    	break;
	    case CYLINDRICAL_REMAP:
		g[0] = G*M/r;
	    	break;
	    case SPHERICAL_REMAP:
		g[0] = G*M/(r*r);
	    	break;
	    }
	    break;
	case 2:
	    switch (remap)
	    {
	    case IDENTITY_REMAP:
	    	for (i = 0; i < dim; i++)
		    g[i] = G*M*g[i]/(r*r);
	    	break;
	    case CYLINDRICAL_REMAP:
	    	for (i = 0; i < dim; i++)
		    g[i] = G*M*g[i]/(r*r*r);
	    	break;
	    }
	    break;
	case 3:
	    for (i = 0; i < dim; i++)
		g[i] = G*M*g[i]/(r*r*r);
	    break;
	}
	if (debugging("gravity"))
	{
	    printf("In generalized_astrophysical_gravity()\n");
	    printf("coords = %f %f\n",coords[0],coords[1]);
	    printf("G = %f  M = %f  r = %f  g = %f %f\n",
			    G,M,r,g[0],g[1]);
	}

        return (const double*)g;
}	/*end generalized_astrophysical_gravity*/

LOCAL const double *radial_gravity(
	const double *coords,
	GRAVITY     *grav_data)
{
	double	*g = grav_data->g;
	double	*cen = grav_data->center;
	double	G = grav_data->G;
	int	i, dim = grav_data->dim;
	double	r;

	for (i = 0; i < dim; i++)
	    g[i] = coords[i] - cen[i];
	r = mag_vector(g,dim);
	for (i = 0; i < dim; i++)
	    g[i] = G*g[i]/r;
	return (const double*)g;
}		/*end radial_gravity*/

EXPORT	boolean source_terms_exist(void)
{
	if (is_gravity() == YES)
	    return YES;
	if (rotational_symmetry() > 0.0)
	    return YES;
	return NO;
}		/*end source_terms_exist*/

EXPORT void set_default_artificial_viscosity(
	AVISC		*avisc)
{
	if ((avisc != NULL) && (avisc != &Fparams.Avisc))
		Fparams.Avisc = *avisc;
	Fparams.Avisc.sp_coef =
		lapidus_stability_factor(Fparams.Avisc.lapidus_visc_coef);
}		/*end set_default_artificial_viscosity*/


EXPORT void default_artificial_viscosity(
	AVISC		*avisc)
{
	if (avisc == NULL)
	    return;
	if (avisc != &Fparams.Avisc)
	    *avisc = Fparams.Avisc;
	avisc->sp_coef = lapidus_stability_factor(avisc->lapidus_visc_coef);
}		/*end default_artificial_viscosity*/


EXPORT	void use_artificial_dissipation(
	AVISC		*avisc)
{
	Gas_param	**params_list;
	int		i, nprms;

	default_artificial_viscosity(avisc);
	nprms = return_params_list(&params_list);
	for (i = 0; i < nprms; i++)
	{
		if (use_lapidus_artificial_viscosity(params_list[i]->avisc))
			use_lapidus_artificial_viscosity(*avisc) = YES;
		if (use_linear_artificial_viscosity(params_list[i]->avisc))
			use_linear_artificial_viscosity(*avisc) = YES;
		if (use_upwind_artificial_viscosity(params_list[i]->avisc))
			use_upwind_artificial_viscosity(*avisc) = YES;
		if (use_muscl_slope_flattening(params_list[i]->avisc))
			use_muscl_slope_flattening(*avisc) = YES;
	}
	return;
}		/*end use_artificial_dissipation*/


/*
*                       set_coord_sys():
*
*	The coordinate remap is visible in the gas code through the
*	coord_system() function.  This function initializes the remap
*	identifier (Fparams.geometry) returned by coord_system().
*	It must be called before coord_system() is used. 
*
*	Input:  remap		- coordinate remap
*	Output: none		- Fparams.geometry is set
*
*/

EXPORT void set_coord_sys(
	int		remap,
	int		dim)
{
	Fparams.dim = dim;
	switch (remap)
	{
	case IDENTITY_REMAP:
	    Fparams.geometry = RECTANGULAR;
	    break;
	case CYLINDRICAL_REMAP:
	    Fparams.geometry = CYLINDRICAL;
	    break;
	case SPHERICAL_REMAP:
	    Fparams.geometry = SPHERICAL;
	    break;
	default:
	    screen("ERROR in set_coord_sys(), Illegal remap %d "
		   "in set_coord_sys\n",remap);
	    clean_up(ERROR);
	    break;
	}
	switch (Fparams.geometry)
	{
	case SPHERICAL:
	    Fparams.alpha = 2.0;
	    break;
	case CYLINDRICAL:
	    Fparams.alpha = 1.0;
	    break;
	case RECTANGULAR:
	default:
	    Fparams.alpha = 0.0;
	}
}		/*end set_coord_sys*/

/*
*		coord_system():
*
*	This function is used to make the coordinate remap visible
*	throughout the gas code.  It can only be called after
*	set_coord_sys() (which is called in init_physics).
*/

EXPORT GEOMETRY coord_system(void)
{
 	if (debugging("coord_system"))
	    screen("coord_system called\n");
	return Fparams.geometry;
}		/*end coord_system*/


EXPORT double rotational_symmetry(void)
{
	return Fparams.alpha;
}		/*end rotational_symmetry*/


EXPORT GEOMETRY Geometry(
	double		*aa)
{
	if (aa != NULL)
	    *aa = Fparams.alpha;
	return Fparams.geometry;
}		/*end Geometry*/

#if defined(HYPRE)
int		gbuffer=8;/*number of buffer grid blocks */
double		ratio = 2.0;/*use coarser grid to calculate gravity */

LOCAL void solve_gravity()
{
	GRAVITY		*gd = Fparams.grav_data;
	double		***g_of_grid = gd->g_of_grid;
	static	double	**den;
	SN_ELLIP        ellip;
	Wave		*wave = gd->roots[0]->wave;
	Front		*front = gd->roots[0]->front;
	RECT_GRID	*gr = wave->rect_grid;
	double		*h = gr->h;
	double		*L = gr->L;
	int		i,j,myid;
	int		icoords[MAXD],*ppgmax=gd->roots[0]->front->pp_grid->gmax;
	int		NRank=-1,WRank=-1,SRank=-1,ERank=-1;
	static Locstate 	st=NULL;
	COMPONENT 	comp;
	int		dim = gr->dim;
	int		ggmax[MAXD];/*number of gravity grid blocks */
	double		coords[MAXD];
	
	
	DEBUG_ENTER(solve_gravity);
#if defined (USE_OVERTURE)
	int		levels = gd->roots[0]->overparam->numberOfRefinementLevels;
	for (i = 0; i < levels-1; i++)
	    ratio /= 2.0;
	printf("amr case\n");
#endif /* defined(USE_OVERTURE) */
	for (i = 0; i < dim; i++)
	    ggmax[i] = (int)rint(gr->gmax[i]/ratio);
	printf("ratio = %f ggmax %d %d gbuffer %d\n", ratio,ggmax[0],ggmax[1],gbuffer);
	printf("\n entering solve_gravity\n");
	if (gd->g_of_grid == NULL)
	{
	    tri_array(&gd->g_of_grid,MAXD,ggmax[0]+2*gbuffer,ggmax[1]+2*gbuffer,FLOAT);
	    g_of_grid = gd->g_of_grid;
	}
	if (den == NULL)
	    bi_array(&den,ggmax[0],ggmax[1],FLOAT);

	for (i = 0; i < ggmax[0]; i++)
	{
	    icoords[0] = i;
	    coords[0] = L[0]+(i+0.5)*h[0]*ratio;
	    for (j = 0; j < ggmax[1]; j++)
	    {
	        icoords[1] = j;
		coords[1] = L[1]+(j+0.5)*h[1]*ratio;
		/*st = Rect_state(icoords,wave); */
		comp = component(coords, front->interf);
		if(st == NULL)
		    alloc_state(front->interf,&st,front->sizest);
		hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,front,
				wave,st,NULL);
		den[i][j] = density(st);
		/*if (debugging("test")) */
		    /*printf("den[%d,%d] = %f coords %f %f\n",i,j,den[i][j],coords[0],coords[1]); */
	    }
	}
	/*set neighbour ranks */
	myid = pp_mynode();
	NRank = myid + ppgmax[0];
	if (NRank > ppgmax[0] * ppgmax[1] -1)
	    NRank = -1;
	SRank = myid - ppgmax[0];
	if (SRank < 0)
	    SRank = -1;
	if (myid%ppgmax[0] != 0)
	    WRank = myid - 1;
	if ((myid+1)%ppgmax[0] != 0)
	    ERank = myid + 1;
	
	ellip.setParameters(gr->GL,gr->GU,gr->L,gr->U, ggmax[0],ggmax[1],gd->G,den);
        ellip.setParaParameters(MPI_COMM_WORLD, ppgmax[0], ppgmax[1], NRank, WRank, SRank, ERank);
	printf("ppgmax %d %d \n",ppgmax[0], ppgmax[1]);
        ellip.solve();
	ellip.setFluxBuffer(gbuffer);
	ellip.getFlux(g_of_grid[0],g_of_grid[1]);
	char filename[100];
        sprintf(filename, "solution.plt");
        ellip.print(filename);
			 
	if (debugging("gravity"))
	{
	    /*print_gravity(NULL); */
	    for (i = gbuffer; i < ggmax[0]+gbuffer; i++)
            {
	        for (j = gbuffer; j < ggmax[1]+gbuffer; j++)
                {
		    if(i == j)
		        printf("g_of_grid[%d,%d] = %f %f\n",i,j,g_of_grid[0][i][j],g_of_grid[1][i][j]);
	        }
	    }
	}
	printf("\n Leaving solve_gravity()\n");
	DEBUG_LEAVE(solve_gravity)
	return;
}	/*end solve_gravity*/

LOCAL void interpolate_gravity(
	const 	double	*coords,
		double*	g,
		double*	h,
		int	dim)
{
        GRAVITY         *gd = Fparams.grav_data;
	RECT_GRID	*gr = gd->roots[0]->wave->rect_grid;
	double		*L = gd->roots[0]->wave->rect_grid->L;
        double           ***g_of_grid = gd->g_of_grid;
	double		**fx = g_of_grid[0],**fy = g_of_grid[1];		
	double 		p01,p23;
	double 		pcoords[4][MAXD];
	double 		pg[4][MAXD];
	double 		x[4*3],b[3],y[4];
	int 		i,j,n=4,p=3;
	int		ix,iy;/*icoords of p0 */
	int		debug_gravity = NO;
	int		ggmax[MAXD];/*number of gravity grid blocks */
	int		myid;
	int		*ppgmax=gd->roots[0]->front->pp_grid->gmax;
	
	for (i = 0; i < dim; i++)
	    ggmax[i] = (int)rint(gr->gmax[i]/ratio);

	/*if (fabs(coords[0]-0.030000)<0.1 && fabs(coords[1]-0.030000)<0.1) */
	  /*  debug_gravity = YES; */
	if (g_of_grid == NULL)
	{
	    for (i = 0; i < dim; i++)
	        g[i] = 0;
	    return;
	}
	/*quad_interpolate from four points */
	/*	p1 ------ p2 */
	/*	   |.p	| */
	/*	p0 ------ p3 */
	
	ix = (int)floor((coords[0]-L[0])/(ratio*h[0])-0.5)+gbuffer;
	iy = (int)floor((coords[1]-L[1])/(ratio*h[1])-0.5)+gbuffer;
	myid = pp_mynode();
	/*make sure four point not exceed upper and right boundaries; */
        if (ix >= ggmax[0]-1+gbuffer && (myid+1)%ppgmax[0] == 0)
	{
	    /*printf("ix exceed, coords %f %f ix %d iy %d myid %d\n", coords[0],coords[1],ix,iy,myid); */
	    /*printf("ratio %f L %f %f h %f %f gbuffer %d\n", ratio,L[0],L[1],h[0],h[1],gbuffer); */
	    ix = ggmax[0]-2+gbuffer;
	}
	if (iy >= ggmax[1]-1+gbuffer && (myid+1)&ppgmax[1] == 0)
	{
	    /*printf("iy exceed, coords %f %f ix %d iy %d myid %d\n", coords[0],coords[1],ix,iy,myid); */
	    iy = ggmax[1]-2+gbuffer;
	}
	if (ix >= ggmax[0]+2*gbuffer)
	{
	    printf("exceed in x direction, coords %f %f ix %d iy %d myid %d\n", coords[0],coords[1],ix,iy,myid);
	    printf("ratio %f L %f %f h %f %f gbuffer %d\n", ratio,L[0],L[1],h[0],h[1],gbuffer);
	}
	
	if (ix < 0 || iy < 0)
	{
	    printf("error, coords %f %f ix %d iy %d index shouldn't be < 0\n", coords[0],coords[1],ix,iy);
	    printf("ratio %f L %f %f h %f %f gbuffer %d\n", ratio,L[0],L[1],h[0],h[1],gbuffer);
	    clean_up(0);
	}
	else /*ix >= 0 && iy >= 0 */
	{
	    /*p0 */
	    pcoords[0][0] = L[0]+(ix-gbuffer+0.5)*ratio*h[0];
	    pcoords[0][1] = L[1]+(iy-gbuffer+0.5)*ratio*h[1];
	    pg[0][0] = fx[ix][iy];
	    pg[0][1] = fy[ix][iy];
	    /*p1 */
	    iy++;
	    pcoords[1][0] = L[0]+(ix-gbuffer+0.5)*ratio*h[0];
	    pcoords[1][1] = L[1]+(iy-gbuffer+0.5)*ratio*h[1];
	    pg[1][0] = fx[ix][iy];
	    pg[1][1] = fy[ix][iy];
	    /*p2 */
	    ix++;
	    pcoords[2][0] = L[0]+(ix-gbuffer+0.5)*ratio*h[0];
	    pcoords[2][1] = L[1]+(iy-gbuffer+0.5)*ratio*h[1];
	    pg[2][0] = fx[ix][iy];
	    pg[2][1] = fy[ix][iy];
	    /*p3 */
	    iy--;
	    pcoords[3][0] = L[0]+(ix-gbuffer+0.5)*ratio*h[0];
	    pcoords[3][1] = L[1]+(iy-gbuffer+0.5)*ratio*h[1];
	    pg[3][0] = fx[ix][iy];
	    pg[3][1] = fy[ix][iy];
	}
	if (debug_gravity)
	{
	    for (i = 0; i < 4; i++)
	        printf("i=%d pcoords %f %f pg %f %f\n",i,pcoords[i][0],pcoords[i][1],pg[i][0],pg[i][1]);
	}
	/*four points interpolate, first y direction, then x direction */
	for (i = 0; i < dim; i++)
	{
	    p01 = (coords[1]-pcoords[1][1])/(ratio*h[1])*pg[0][i] + (pcoords[0][1]-coords[1])/(ratio*h[1])*pg[1][i];
	    p23 = (coords[1]-pcoords[2][1])/(ratio*h[1])*pg[3][i] + (pcoords[3][1]-coords[1])/(ratio*h[1])*pg[2][i];
 	    g[i] = (coords[0]-pcoords[2][0])/(ratio*h[0])*p01 + (pcoords[1][0]-coords[0])/(ratio*h[0])*p23;
	}

	return;	
}	/*end interpolate_gravity*/
#endif /* defined(HYPRE) */

EXPORT void g_set_gravity_charts(CHART** roots, int nroots)
{
        GRAVITY *gd = Fparams.grav_data;
	int     step = roots[0]->front->step;

        if (gd == NULL)
            return;

        gd->roots = roots;
        gd->nroots = nroots;
        if(gd->type==GENERALIZED_ASTROPHYSICAL_GRAVITY)
	{
#if defined(HYPRE)
	    if (debugging("2d_gravity"))
		solve_gravity();
	    else
#endif /* defined(HYPRE) */
                eval_partial_mass(NULL, 0, 0);
	}
}	/*end g_set_gravity_charts*/

EXPORT boolean eval_partial_mass(double* Mp,
			      int N,
			      int mass_type)
{
	boolean 		ret;
	const int       myid = pp_mynode();
	GRAVITY*	gd = Fparams.grav_data;

	/*if(gd->type!=GENERALIZED_ASTROPHYSICAL_GRAVITY)
		return NO;*/
	
	if(Mp==NULL)
	{
		if(gd->Mp==NULL)
			uni_array(&gd->Mp, Fparams.grav_data->N, FLOAT);
		Mp = gd->Mp;
		N  = gd->N;
	}

	switch(Fparams.grav_data->roots[0]->wave->rect_grid->dim)
	{
		case 1:  ret = eval_partial_mass1D(Mp,N,mass_type); break;
		case 2:  ret = eval_partial_mass2D(Mp,N,mass_type); break;
		case 3:  ret = eval_partial_mass3D(Mp,N,mass_type); break;
		default: ret = NO;
	}
	return ret;
}	/*end eval_partial_mass*/

/* flame */

LOCAL boolean eval_partial_mass1D(double* Mp,
			       int N,
			       int mass_type)
{
	int 		i, j;
	int 		icoords[] = {0,0,0};
	double		*coords;
	double 		r, R;
	Locstate        state;
	double 		dens;
	int 		imin, k;
	double 		dx;
	const int       myid = pp_mynode();
	const int       nn = pp_numnodes();
	static double	*pRecvBuff = NULL;
	Wave		*wave = Fparams.grav_data->roots[0]->wave;
	int		*gmax = Fparams.grav_data->gmax;
	int		*size = (nn>1)?wave->pp_grid->Zoom_grid.gmax:gmax;
	double		*cen = Fparams.grav_data->center;
	double		*GL = Fparams.grav_data->GL;
	double		*GU = Fparams.grav_data->GU;

	memset(Mp, 0, N*FLOAT);
	if(wave->wave_pointers.tri_soln==NULL)
		return NO;

	R = (GU[0]-GL[0]);
	dx = R/gmax[0];

	for(i=0; i<size[0]; i++)
	{
		icoords[0] = i;
		state = Rect_state(icoords,wave);
		if(!is_obstacle_state(state))
		{
#if defined(COMBUSTION_CODE)
			if((Burned(state)==YES)&&(mass_type==UNBURNED))
				dens = 0;
			else if((Burned(state)==NO)&&(mass_type==BURNED))
				dens = 0;
			else
#endif /* defined(COMBUSTION_CODE) */
				dens = density(state);
		}
		else
			dens = 0;

		coords = Rect_coords(icoords, wave);
		r = sqrt((coords[0]-cen[0])*(coords[0]-cen[0])+
			 (coords[1]-cen[1])*(coords[1]-cen[1])+
			 (coords[2]-cen[2])*(coords[2]-cen[2]));
		imin = (int)(r*(N-1)/R);
		for(k=imin+1; k<N; k++)
				Mp[k] += dens*4*PI*dx*r*r;
	}

	if(nn>1)
	{
		pp_gsync();
		if(is_io_node(myid))
		{
			if(pRecvBuff==NULL)
				uni_array(&pRecvBuff, N, FLOAT);
			memset(pRecvBuff, 0, N*FLOAT);
			for(i=0;i<nn;i++)
				if(i!=myid)
				{
					pp_recv(i,i,pRecvBuff,N*FLOAT);
					for(j=0;j<N;j++)
						Mp[j] += pRecvBuff[j];
				}
			for(i=0;i<nn;i++)
				if(i!=myid)
					pp_send(i,Mp,N*FLOAT,i);
		}
		else
		{
			pp_send(myid,Mp,N*FLOAT,IO_NODE_ID);
			pp_recv(myid,IO_NODE_ID,Mp,N*FLOAT);
		}
		pp_gsync();
	}

	return YES;
} /* end eval_partial_mass1D */


/* flame */
LOCAL boolean eval_partial_mass2D(double* Mp,
			       int N,
			       int mass_type)
{
	int 		i, j, k, n, nsec;
	int 		icoords[] = {0,0,0};
	double		*coords;
	static	double	*ave_dens;
	static	int	*num_cells;
	Locstate        state;
	Front 		*front;
	Wave 		*wave;
	RECT_GRID	*rgrid;
	const int 	*gmax;
	COMPONENT 	comp;
	double 		R, r, dR;
	GRAVITY		*gd = Fparams.grav_data;
	double		*cen = gd->center;
	const double	*GL, *GU;
	const double	*L, *U;
#if defined(__MPI__)
	const int       myid = pp_mynode();
	const int       nn = pp_numnodes();
	static double	*pRecvBuff = NULL;
	static int	*nRecvBuff = NULL;
#endif /* if defined(__MPI__) */

	debug_print("gravity","Entering eval_partial_mass2D()\n");

	if (num_cells == NULL)
	{
	    uni_array(&num_cells,N,INT);
	    uni_array(&ave_dens,N,FLOAT);
#if defined(__MPI__)
	    uni_array(&pRecvBuff,N,FLOAT);
	    uni_array(&nRecvBuff,N,INT);
#endif /* if defined(__MPI__) */
	}
	memset(Mp,0,N*FLOAT);
	memset(num_cells,0,N*INT);
	memset(ave_dens,0,N*FLOAT);

	GL = gd->GL;
	GU = gd->GU;
	R = 0.0;
	R = max(R,sqrt(sqr(GL[0]-cen[0]) + sqr(GL[1]-cen[1])));
	R = max(R,sqrt(sqr(GL[0]-cen[0]) + sqr(GU[1]-cen[1])));
	R = max(R,sqrt(sqr(GU[0]-cen[0]) + sqr(GU[1]-cen[1])));
	R = max(R,sqrt(sqr(GU[0]-cen[0]) + sqr(GL[1]-cen[1])));
	gd->R = R;
	gd->dR = dR = R/(N-1);

	for(n=gd->nroots-1; n>=0; n--)
	{
	    front = gd->roots[n]->front;
	    wave = gd->roots[n]->wave;
	    rgrid = wave->rect_grid;
	    gmax = rgrid->gmax;

	    if(wave->wave_pointers.tri_soln==NULL)
		return NO;

	    for(i = 0; i < gmax[0]; i++)
	    {
		icoords[0] = i;
		for(j = 0; j < gmax[1]; j++)
		{
		    icoords[1] = j;
		    coords = Rect_coords(icoords,wave);
		    comp = Rect_comp(icoords,wave);
		    r = sqrt(sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1]));
		    nsec = (int)(r/dR);
		    state = Rect_state(icoords,wave);

		    if(!is_obstacle_state(state))
		    {
			    ave_dens[nsec] += Dens(state);
			    num_cells[nsec] += 1;
		    }
		}
	    }
	    L = rgrid->L;
	    U = rgrid->U;

	    /* If the gravity center is in (sub)domain, but no
	     * regular cell center is near the gravity center */

	    if (L[0] <= cen[0] && cen[0] <= U[0] &&
		L[1] <= cen[1] && cen[1] <= U[1] && num_cells[0] == 0)
	    {
		COMPONENT comp = component(cen,front->interf);
		hyp_solution(cen,comp,NULL,UNKNOWN_SIDE,front,
				wave,state,NULL);
		ave_dens[0] = Dens(state);
		num_cells[0] = 1;
	    }
	}


#if defined(__MPI__)
	if(nn > 1)
	{
	    pp_gsync();
	    pp_send_all(myid,ave_dens,N*FLOAT);
	    pp_gsync();
	    for(i = 0; i < nn; i++)
	    {
		if(i != myid)
		{
		    pp_recv(i,i,pRecvBuff,N*FLOAT);
		    for(j = 0; j < N; j++)
			ave_dens[j] += pRecvBuff[j];
		}
	    }
	    pp_gsync();
	    pp_send_all(myid,num_cells,N*INT);
	    pp_gsync();
	    for(i = 0; i < nn; i++)
	    {
		if(i != myid)
		{
		    pp_recv(i,i,nRecvBuff,N*INT);
		    for(j = 0; j < N; j++)
			num_cells[j] += nRecvBuff[j];
		}
	    }
	    pp_gsync();
	}
#endif /* if defined(__MPI__) */
	
	if (debugging("gravity"))
	{
	    printf("In eval_partial_mass2D(), number of shells = %d\n",N);
	    printf("Index_of_shell   num_of_cells   accum_dens\n");
	    for (i = 0; i < N; ++i)
	    {
	    	printf("        %2d              %2d    %12.8f\n",
				i,num_cells[i],ave_dens[i]);
	    }
	}

	for(k = 0; k < N; k++)
	{
	    if (num_cells[k] != 0)
		ave_dens[k] /= (double)num_cells[k];
	    else
		ave_dens[k] = ave_dens[k-1];
	}

	switch (rgrid->Remap.remap)
	{
	case CYLINDRICAL_REMAP:
	    Mp[0] = 1.333333*PI*dR*dR*dR*ave_dens[0];
	    for(k = 1; k < N; k++)
	    {
	    	r = k*dR;
	    	Mp[k] = Mp[k-1] + 4.0*PI*r*r*dR*ave_dens[k];
	    }
	    break;
	case IDENTITY_REMAP:
	    Mp[0] = PI*dR*dR*ave_dens[0];
	    for(k = 1; k < N; k++)
	    {
	    	r = k*dR;
	    	Mp[k] = Mp[k-1] + 2.0*PI*r*dR*ave_dens[k];
	    }
	    break;
	}
	if (debugging("gravity"))
	{
	    printf("In eval_partial_mass2D()\n");
	    printf("Index_of_shell   ave_dens   mass_within\n");
	    for (i = 0; i < N; ++i)
	    {
	    	printf("        %2d      %12.8f  %12.8f\n",
				i,ave_dens[i],Mp[i]);
	    }
	}

	debug_print("gravity","Leaving eval_partial_mass2D()\n");
	return YES;
} /* end eval_partial_mass2D */

LOCAL boolean eval_partial_mass3D(double* Mp,
                               int N,
                               int mass_type)
{
        printf("ERROR: 3D not implemented\n");
        clean_up(ERROR);
        return NO;
} /* end eval_partial_mass3D */

LOCAL void print_circular_mass()
{
	GRAVITY*        gd = Fparams.grav_data;	
	double *		Mp = gd->Mp;
	int		N = gd->N, i;

	printf("begin print circular mass:\n");
	for (i = 0; i < N; i++)
	{
		printf("Mp[%d] = %- "FFMT"\n", i, Mp[i]);
	}
	printf("end print circular mass.\n");
	return;
}

LOCAL void print_grid_gravity(INTERFACE* intfc)
{
	RECT_GRID	*gr;
	int 		i,j,k,icoords[MAXD],dim;
	double		*grav, coords[MAXD];
	
	gr = computational_grid(intfc);
	dim = gr->dim;
	for (i = 0; i < 1; i++)/*x */
	{
	    icoords[0] = i;
	    for (j = 0; j < gr->gmax[1]; j++)/*y */
	    {
		icoords[1] = j;
		for (k = 0; k < dim; k++)
		    coords[k] = cell_center(icoords[k], k, gr);
		grav = (double*)gravity(coords,0);
		printf("icoord: %d %d %f %f %f gravity = %- "FFMT" %- "FFMT" %- "FFMT" \n",
                        icoords[0],icoords[1],coords[0],coords[1],coords[2],grav[0], grav[1], grav[2]);
	    }
	}
}

EXPORT void print_gravity(INTERFACE * intfc)
{

	CURVE  **cur;
	int i,k;
	double pp[4][3]={{0.01,0.01,0.},
	   	    {2.98,2.98,0.},
		    {0.0299,2.98,0.},
		    {2.98,0.0299,0.}};
	double *grav;
	double ggg[2];

	printf("first some specific points\n");
	for (i = 0; i < 4; i++)
	{
 	    grav = (double*)gravity(pp[i],0);
	    printf("%- "FFMT" %- "FFMT" %- "FFMT" gravity = %- "FFMT" %- "FFMT" %- "FFMT" \n",
			    pp[i][0],pp[i][1],pp[i][2],
			    grav[0], grav[1], grav[2]);
	}
	/*for (i = 0; i < 49; i++)
	{
		ggg[0] = ggg[1] = 0.03001+0.06*i;
		grav = (double*)gravity(ggg,0);
		printf(" coords %f %f gravity %f %f\n", ggg[0],ggg[1],grav[0],grav[1]);
	}*/
				
	/*TMP for print mass distribution in circuler sectors */
	print_circular_mass();

	/*TMP for print gravity in Interior grid */
	/*print_grid_gravity(intfc); */
	
	/*TMP for print gravity on Interior curve */
	
	/*printf(" Print Interior Curve Gravity\n");
	if ((cur = intfc->curves) != NULL)
            while (*cur) 
	    {
		if(!is_bdry(*cur))
		    print_curve_gravity(*cur);
		cur++;
	    }
	printf(" End Print Gravity\n");*/

	return;
}

LOCAL void print_curve_gravity(
		CURVE          *curve)
{
	
	BOND            *bond;
	POINT           *p;
	double		*grav;
	int		i, k, dim;
               
	dim = curve->interface->dim;
	
	for (i = 1,bond = curve->first; bond != NULL; bond = bond->next,++i)
        {
            p = bond->start;
            for (k = 0; k < dim; ++k)
                (void) printf("%- "FFMT" ",Coords(p)[k]);
	    /*eval_gravity(Coords(p),0,grav); */
	    grav = (double*)gravity(Coords(p),0);
	    printf("gravity = ");
	    for (k = 0; k < dim; ++k)
                (void) printf("%- "FFMT" ",grav[k]);
	    (void) printf("\n");
        }	
                                                                                                               
    	if ((curve->last != NULL) && (curve->last->end != NULL))
	{
            p = curve->last->end;
	    for (k = 0; k < dim; ++k)
                (void) printf("%- "FFMT" ",Coords(p)[k]);
	    /*eval_gravity(Coords(p),0,grav); */
	    grav = (double*)gravity(Coords(p),0);
            printf("gravity = ");
            for (k = 0; k < dim; ++k)
                (void) printf("%- "FFMT" ",grav[k]);
            (void) printf("\n");
	}
	
	return;
}
