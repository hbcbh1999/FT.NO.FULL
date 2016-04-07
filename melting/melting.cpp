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
*				crystal.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This is example of three circles all moving the a normal velocity.
*	Bifurcation occurs when they meet each other. FronTier solves
*	the bifurcation automatically.
*
*/

#include "melting.h"

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/
#define		MAX_NUM_VERTEX_IN_CELL		20

typedef struct {
        double center[3];
        double radius;
} TEST_SPHERE_PARAMS;

	/*  Local Application Function Declarations */

static void 	temperature_main_driver(Front*,TEST_SPHERE_PARAMS,CARTESIAN&);
static void	readPhaseParams(char*,PARAMS*);
static void	temperature_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	record_1d_growth(Front*,double***,int*);
static void 	plot_growth_data(char*,double**,int);
static double 	temperature_func(double*,COMPONENT,double);
static double 	crystal_curve(POINTER,double*);
static double 	sine_curve(POINTER,double*);
static double 	sphere_surf(POINTER,double*);
static boolean 	fractal_dimension(Front*,TEST_SPHERE_PARAMS,double*,double*);
static void	initPhaseIntfc(int,char*,PARAMS*);
static void	initPhaseIntfc1d(char*,PARAMS*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;
boolean binary = YES;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TEST_SPHERE_PARAMS s_params;
	PARAMS eqn_params;

	CARTESIAN       cartesian(front);

	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.size_of_intfc_state = sizeof(STATE);

	in_name      		= f_basic.in_name;
	restart_state_name      = f_basic.restart_state_name;
        out_name     		= f_basic.out_name;
        restart_name 		= f_basic.restart_name;
        RestartRun   		= f_basic.RestartRun;
        ReadFromInput   	= f_basic.ReadFromInput;
	RestartStep 		= f_basic.RestartStep;

	sprintf(restart_state_name,"%s/state.ts%s",restart_name,
			right_flush(RestartStep,7));
	sprintf(restart_name,"%s/intfc-ts%s",restart_name,
			right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
        sprintf(restart_state_name,"%s-nd%s",restart_state_name,
			right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

        FT_ReadSpaceDomain(in_name,&f_basic);
        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	if (debugging("trace")) printf("Passed PetscInitialize()\n");

	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);
	readPhaseParams(in_name,&eqn_params);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    initPhaseIntfc(f_basic.dim,in_name,&eqn_params);
	    if (debugging("trace")) printf("Passed initPhaseIntfc()\n");

	    FT_InitIntfc(&front,&level_func_pack);
	}
	if (debugging("trace")) printf("Passed FT_InitIntfc()\n");

        Time_step_factor(&front);
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	front.extra2 = (POINTER)&eqn_params;

	eqn_params.dim = f_basic.dim;

	velo_func_pack.func_params = (POINTER)&eqn_params;
	velo_func_pack.func = NULL;

        cartesian.initMesh();
	if (RestartRun)
	{
	    FT_ParallelExchIntfcBuffer(&front);
	    cartesian.readFrontInteriorState(restart_state_name);
	}
	else
	{
	    cartesian.setInitialCondition();
	}
	if (debugging("trace"))
	    printf("Passed state initialization\n");


	FT_InitVeloFunc(&front,&velo_func_pack);

	if (debugging("trace")) printf("Passed FT_InitVeloFunc()\n");

        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */

        front._point_propagate = temperature_point_propagate;

	/* Propagate the front */

	temperature_main_driver(&front,s_params,cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void temperature_main_driver(
        Front *front,
	TEST_SPHERE_PARAMS s_params,
	CARTESIAN &cartesian)
{
        int ip,im,status;
        Front *newfront;
        double dt,dt_frac,CFL;
        boolean is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        double fcrds[MAXD];
        int  i,dim = front->rect_grid->dim;
	INTERFACE *grid_intfc;
	PARAMS *eqn_params;
	double temperature_time_limit; 
	double frac_dim,radius;
	FILE *Radius_file,*FracDim_file;
	boolean bdry_reached = NO;
	double **growth_data = NULL;
	int count = 0;

	if (debugging("trace"))
	    printf("Entering temperature_main_driver()\n");
	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

	eqn_params = (PARAMS*)front->vparams;

        if (!RestartRun)
        {
	    FT_ResetTime(front);
            FT_SetOutputCounter(front);
            // Front standard output
            FT_Save(front,out_name);
	    if (dim != 1)
            	FT_AddMovieFrame(front,out_name,binary);
	    else
	    	cartesian.oneDimPlot(out_name);

	    /*
	    if (dim == 1)
	    {
	    	bdry_reached = NO;
		record_1d_growth(front,&growth_data,&count);
	    }
	    else
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
						&radius);
	    */
	    if (pp_mynode() == 0 && dim != 1)
	    {
		char fname[200];
		sprintf(fname,"%s/radius",out_name);
	    	Radius_file = fopen(fname,"w");
		sprintf(fname,"%s/FracDim",out_name);
	    	FracDim_file = fopen(fname,"w");

		fprintf(Radius_file,"\"Crystal radius\"\n\n");
		fprintf(Radius_file,"%f  %f\n",front->time,radius);
		fprintf(FracDim_file,"\"Fractal dimension\"\n\n");
		fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
		fflush(Radius_file);
		fflush(FracDim_file);
	    }

	    // Problem specific output

	    FT_Propagate(front);
	    FT_SetTimeStep(front);
	    cartesian.setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*cartesian.m_dt);
        }
        else
	    FT_SetOutputCounter(front);
	if (debugging("trace"))
	    printf("Passed second restart check()\n");

	FT_TimeControlFilter(front);

        for (;;)
        {
	    if (debugging("trace"))
		printf("Before FrontAdvance()\n");
	    FT_Propagate(front);
	    if (debugging("trace"))
		printf("After FrontAdvance(), before solve()\n");
	    cartesian.solve(front->dt);
	    if (debugging("trace"))
		printf("After solve()\n");

	    FT_AddTimeStepToCounter(front);

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*cartesian.m_dt);

            printf("\ntime = %f   step = %7d   dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (FT_IsSaveTime(front))
	    {
		// Front standard output
                FT_Save(front,out_name);
		// Problem specific output
		cartesian.printFrontInteriorState(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
	    {
		// Front standard output
		if (dim != 1)
                    FT_AddMovieFrame(front,out_name,binary);
	    	else
	    	    cartesian.oneDimPlot(out_name);
	    }
	    if (dim == 1)
	    {
	    	bdry_reached = NO;
		record_1d_growth(front,&growth_data,&count);
	    }
	    else
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
					&radius);

	    if (bdry_reached) front->time_limit_reached = YES;
	    if (pp_mynode() == 0 && dim != 1)
	    {
		fprintf(Radius_file,"%f  %f\n",front->time,radius);
		fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
		fflush(Radius_file);
		fflush(FracDim_file);
	    }

            if (FT_TimeLimitReached(front))
	    {
		if (dim == 1)
		    plot_growth_data(out_name,growth_data,count);
                break;
	    }
	    /* Output section, next dt may be modified */

	    FT_TimeControlFilter(front);
        }
	if (pp_mynode() == 0 && dim != 1)
	{
	    fclose(Radius_file);
	    fclose(FracDim_file);
	}
}       /* end temperature_main_driver */

LOCAL double crystal_curve(
        POINTER func_params,
        double *coords)
{

        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
        double dist, theta;
	double *cen = s_params->center;
	double radius = s_params->radius;

        dist =   sqrt(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]));
        theta = asin(fabs(coords[1]-0.5)/dist);
	if (coords[0]-0.5 < 0 && coords[1]-0.5 > 0)
	    theta = PI - theta;
	else if (coords[0]-0.5 < 0 && coords[1]-0.5 < 0)
	    theta = PI + theta;
	else if (coords[0]-0.5 > 0 && coords[1]-0.5 < 0)
	    theta = 2*PI - theta;
        return dist - radius + .003*sin(6.0*theta);
}       /* end crystal_curve */

extern  double temperature_of_state(
        POINTER state)
{
        STATE *temperature_state = (STATE*)state;
        return (double)(*temperature_state);
}       /* end temperature_of_state */

static	double temperature_func(
	double *coords,
	COMPONENT comp,
	double T_l)
{
	if (comp != LIQUID_COMP) return 0.0;
	return T_l;
}	/* end temperature_init_func */

static  void temperature_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double nor_l[MAXD],nor_s[MAXD];
        double p1[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
	double grad_s,grad_l;
        double s0,s1;
        STATE *sl,*sr,*state;
	PARAMS *eqn_params = (PARAMS*)front->vparams;
        double *temperature = eqn_params->field->temperature;
        double T_si = eqn_params->Ti[0];
        double T_li = eqn_params->Ti[1];
        double rho_s = eqn_params->rho[0];
        double rho_l = eqn_params->rho[1];
        double k_s = eqn_params->k[0];
        double k_l = eqn_params->k[1];
        double Cp_s = eqn_params->Cp[0];
        double Cp_l = eqn_params->Cp[1];
        double L = eqn_params->L[0];
	double H,kappa;

        if (wave_type(oldhs) < FIRST_PHYSICS_WAVE_TYPE)
        {
            for (i = 0; i < dim; ++i)
                Coords(newp)[i] = Coords(oldp)[i];
            FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
            state = (STATE*)left_state(newp);
            *state = *sl;
            state = (STATE*)right_state(newp);
            *state = *sr;
            return;
        }

        GetFrontNormal(oldp,oldhse,oldhs,nor_l,front);
        if (negative_component(oldhs) == LIQUID_COMP)
            for (i = 0; i < dim; ++i)
                nor_l[i] *= -1.0;
        dn = grid_size_in_direction(nor_l,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor_l[i]*dn;
        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == LIQUID_COMP) ? sl : sr;
        s0 = (double)(*state);
	if (!FT_IntrpStateVarAtCoords(front,LIQUID_COMP,p1,temperature,
				temperature_of_state,&s1))
            s1 = (double)(*state);
        grad_l = (s1 - s0)/dn;

        GetFrontNormal(oldp,oldhse,oldhs,nor_s,front);
        if (negative_component(oldhs) == SOLID_COMP)
            for (i = 0; i < dim; ++i)
                nor_s[i] *= -1.0;
        dn = grid_size_in_direction(nor_s,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor_s[i]*dn;
        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == SOLID_COMP) ? sl : sr;
        s0 = (double)(*state);
	if (!FT_IntrpStateVarAtCoords(front,SOLID_COMP,p1,temperature,
				temperature_of_state,&s1))
            s1 = (double)(*state);

        grad_s = (s1 - s0)/dn;
	H = Cp_l*k_l*grad_l + Cp_s*k_s*grad_s;
        for (i = 0; i < dim; ++i)
        {
	    if (H > 0)
		vel[i] = H/L/rho_s*nor_s[i];
	    else
		vel[i] = H/L/rho_l*nor_s[i];
	}
	printf("grad_l = %f  grad_s = %f  vel = %f\n",grad_l,grad_s,vel[0]);

        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }

	/* Update the state of the new interface point */
	state = (STATE*)left_state(newp);
	if (negative_component(oldhs) == LIQUID_COMP)
            *state = (STATE)eqn_params->Ti[0];
	else if (negative_component(oldhs) == SOLID_COMP)
            *state = (STATE)eqn_params->Ti[0];
	else
            *state = 0.0;

	state = (STATE*)right_state(newp);
	if (positive_component(oldhs) == LIQUID_COMP)
            *state = (STATE)eqn_params->Ti[0];
	else if (positive_component(oldhs) == SOLID_COMP)
            *state = (STATE)eqn_params->Ti[0];
	else
            *state = 0.0;
}       /* temperature_point_propagate */


LOCAL double sine_curve(
        POINTER func_params,
        double *coords)
{

        double dist, theta;

	double y_intfc = 0.5 + 0.2*sin(6.0*PI*coords[0]);

        dist =   y_intfc - coords[1];
        return dist;
}       /* end sine_curve */


static double sphere_surf(
        POINTER func_params,
        double *coords)
{
        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
	double x0,y0,z0,R;
	double distance;

        x0 = s_params->center[0];
        y0 = s_params->center[1];
        z0 = s_params->center[2];
	R = s_params->radius;

	distance = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0) +
			sqr(coords[2] - z0)) - R;

        return distance;

}       /* end sphere_surf */



static	void	readPhaseParams(
	char *in_name,
	PARAMS *eqn_params)
{
	FILE *infile;
	char scheme[200];
	int i,num_phases;
	char string[200];

	infile = fopen(in_name,"r");

	CursorAfterString(infile,"Enter number of phases:");
	fscanf(infile,"%d",&num_phases);
	eqn_params->num_phases = num_phases;
	(void) printf(" %d\n",num_phases);
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->T0,num_phases,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->rho,num_phases,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->Cp,num_phases,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->k,num_phases,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->Ti,num_phases-1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&eqn_params->L,num_phases-1,sizeof(double));
	for (i = 0; i < num_phases; ++i)
	{
	    sprintf(string,"Enter ambient temperature of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->T0[i]);
	    (void) printf("%f\n",eqn_params->T0[i]);
	    if (i > 0 && same_sign(
			eqn_params->T0[i]-eqn_params->Ti[i-1],
			eqn_params->T0[i-1]-eqn_params->Ti[i-1]))
	    {
		(void) printf("ERROR: same adjacent phases!\n");
		clean_up(ERROR);
	    }
	    sprintf(string,"Enter density of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->rho[i]);
	    (void) printf("%f\n",eqn_params->rho[i]);
	    sprintf(string,"Enter specific heat of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->Cp[i]);
	    (void) printf("%f\n",eqn_params->Cp[i]);
	    sprintf(string,"Enter thermal conductivity of phase %d:",i+1);
	    CursorAfterString(infile,string);
	    fscanf(infile,"%lf",&eqn_params->k[i]);
	    (void) printf("%f\n",eqn_params->k[i]);
	    if (i != num_phases-1)
	    {
	        sprintf(string,"Enter melting temperature of interface %d:",
						i+1);
	        CursorAfterString(infile,string);
	    	fscanf(infile,"%lf",&eqn_params->Ti[i]);
	        (void) printf("%f\n",eqn_params->Ti[i]);
	        sprintf(string,"Enter latent heat at interface %d:",i+1);
	        CursorAfterString(infile,string);
	    	fscanf(infile,"%lf",&eqn_params->L[i]);
	        (void) printf("%f\n",eqn_params->L[i]);
	    }
	}

	CursorAfterString(infile,"Choose numerical scheme");
	(void) printf("\n");
	CursorAfterString(infile,"Enter scheme:");
	fscanf(infile,"%s",scheme);
	(void) printf("%s\n",scheme);
	if ((scheme[0] == 'E' || scheme[0] == 'e') &&
	    (scheme[1] == 'X' || scheme[1] == 'x')) 
	    eqn_params->num_scheme = UNSPLIT_EXPLICIT;
	else if ((scheme[0] == 'I' || scheme[0] == 'i') &&
	    (scheme[1] == 'M' || scheme[1] == 'm')) 
	    eqn_params->num_scheme = UNSPLIT_IMPLICIT;
	else if ((scheme[0] == 'C' || scheme[0] == 'c') &&
	    (scheme[1] == 'N' || scheme[1] == 'n')) 
	    eqn_params->num_scheme = CRANK_NICOLSON;
	else
	{
	    printf("Unknown numerical scheme!\n");
	    clean_up(ERROR);
	}
	fclose(infile);
}

static boolean fractal_dimension(
	Front *front,
	TEST_SPHERE_PARAMS s_params,
	double *frac_dim,
	double *radius)
{
	double coords[MAXD],*center = s_params.center;
	double dist,r_sqr,r_max,r_min = s_params.radius;
	INTERFACE *grid_intfc,*intfc = front->interf;
	int i,j,k,*gmax,dim = intfc->dim;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	struct Table *T;
	RECT_GRID *grid;
	double *L,*U,*h;
	COMPONENT *gr_comp,comp;
	int N,Nc;
	boolean grid_intfc_made = NO;
	double ratio;
	boolean bdry_reached;
	double margin[MAXD];

	/* Find maximum radius of crystal growth */
	r_max = 0.0;
	grid = computational_grid(intfc);
	L = grid->GL;
	U = grid->GU;
	for (i = 0; i < dim; ++i)
	    margin[i] = 0.01*(U[i] - L[i]);
	bdry_reached = NO;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE)
	    	continue;
	    r_sqr = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
		if (Coords(p)[i] >= U[i] - margin[i] ||
		    Coords(p)[i] <= L[i] + margin[i])
		    bdry_reached = YES;
	    	r_sqr += sqr(Coords(p)[i] - center[i]);
	    }
	    if (r_max < r_sqr) r_max = r_sqr;
	}
	r_max = sqrt(r_max);
	*radius = r_max;
#if defined (__MPI__)
	pp_global_max(radius,1);
#endif /* defined (__MPI__)

	/* Preparation for counting */

        if (front->grid_intfc == NULL)
        {
            FT_MakeGridIntfc(front);
            grid_intfc = front->grid_intfc;
            grid_intfc_made = YES;
        }
        else
            grid_intfc = front->grid_intfc;
        grid = &topological_grid(grid_intfc);
        gmax = grid->gmax;
        L = grid->L;
        h = grid->h;
	T = table_of_interface(grid_intfc);
        gr_comp = T->components;

	/* Start counting */
	N = Nc = 0;
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= gmax[0]; ++i)
            for (j = 0; j <= gmax[1]; ++j)
            {
	    	comp = gr_comp[d_index2d(i,j,gmax)];
		coords[0] = L[0] + i*h[0];
                coords[1] = L[1] + j*h[1];
		dist = sqrt(sqr(coords[0] - center[0]) +
			    sqr(coords[1] - center[1]));
	    	if (dist > r_min && dist < r_max)
		{
		    ++N;
		    if (comp == SOLID_COMP)
		    	++Nc;
		}
	    }
	    break;
	case 3:
	    for (i = 0; i <= gmax[0]; ++i)
            for (j = 0; j <= gmax[1]; ++j)
	    for (k = 0; k <= gmax[2]; ++k)
            {
	    	comp = gr_comp[d_index3d(i,j,k,gmax)];
		coords[0] = L[0] + i*h[0];
                coords[1] = L[1] + j*h[1];
		coords[2] = L[2] + k*h[2];
		dist = sqrt(sqr(coords[0] - center[0]) +
			    sqr(coords[1] - center[1]) +
			    sqr(coords[2] - center[2]));
	    	if (dist > r_min && dist < r_max)
		{
		    ++N;
		    if (comp == SOLID_COMP)
		    	++Nc;
		}
	    }
	}
#if defined (__MPI__)
	pp_global_isum(&N,1);
	pp_global_isum(&Nc,1);
#endif /* defined (__MPI__) */
	if (grid_intfc_made)
	    FT_FreeGridIntfc(front);
	ratio = ((double)N)/((double)Nc);
	*frac_dim = (double)dim + log(ratio)/log(h[0]);
	return bdry_reached;
}	/* end fractal_dimension */
	
static void	record_1d_growth(
	Front *front,
	double ***growth_data,
	int *count)
{
	INTERFACE *intfc = front->interf;
	int i,j;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	PARAMS *eqn_params = (PARAMS*)front->vparams;
	static int total_length = 0;
	STATE *sl,*sr;
	
	if (*count >= total_length)
	{
	    static double **tmp_data;
	    total_length += 1000;
	    FT_MatrixMemoryAlloc((POINTER*)&tmp_data,total_length,3,sizeof(double));
	    for (i = 0; i < *count; ++i)
	    for (j = 0; j < 3; ++j)
	    	tmp_data[i][j] = (*growth_data)[i][j];
	    FT_FreeThese(1,*growth_data);
	    *growth_data = tmp_data;
	}

	(*growth_data)[*count][0] = front->time;

	/* Initialize states at the interface */
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (is_bdry(p)) continue;
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (negative_component(hs) == LIQUID_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = (double)(*sl);
	    }
	    else if (positive_component(hs) == LIQUID_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = (double)(*sr);
	    }
	}
	*count += 1;
}	/* end record_1d_growth */

static void plot_growth_data(
	char out_name[100],
	double **growth_data,
	int count)
{
	char fname[100];
	FILE *ofile;
	int i;

	sprintf(fname,"%s/posn.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface position vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][1]);
	fclose(ofile);

	sprintf(fname,"%s/solt.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface temperature vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][2]);
	fclose(ofile);
}	/* end plot_growth_data */

static void initPhaseIntfc(
	int dim,
	char *in_name,
	PARAMS *eqn_params)
{
	switch (dim)
	{
	case 1:
	    return initPhaseIntfc1d(in_name,eqn_params);
	case 2:
	case 3:
	    (void) printf("2D and 3D initialization not implemented\n");
	    clean_up(ERROR);
	}
}	/* end initPhaseIntfc */

static void initPhaseIntfc1d(
	char *in_name,
	PARAMS *eqn_params)
{
	char s[200];
	FILE *infile;
	double **points;
	int i,num_phases = eqn_params->num_phases;
	COMPONENT neg_comp,pos_comp;
	POINT *p;
	char string[200];
	
	infile = fopen(in_name,"r");

	FT_MatrixMemoryAlloc((POINTER*)&points,num_phases-1,MAXD,sizeof(double));
	for (i = 0; i < num_phases-1; ++i)
	{
	    sprintf(string,"Enter position of interface %d:",i+1);
	    CursorAfterString(infile,string);
            fscanf(infile,"%lf",&points[i][0]);
            (void) printf("%f\n",points[i][0]);
	    neg_comp = (eqn_params->T0[i] < eqn_params->Ti[i]) ?
                                SOLID_COMP : LIQUID_COMP;
            pos_comp = (eqn_params->T0[i] > eqn_params->Ti[i]) ?
                                SOLID_COMP : LIQUID_COMP;
            p = make_point(points[i],neg_comp,pos_comp);
            wave_type(Hyper_surf(p)) = FIRST_PHYSICS_WAVE_TYPE;
	}
}	/* end initPhaseIntfc1d */
