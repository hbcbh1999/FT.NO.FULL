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
*				example.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static double tdumbbell_func(POINTER,double*);
static int test_curvature_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/
//Dummbell shape interface, two spheres of radius R centered at (x0,y,z)
//and (x1,y,z) connected by a cylinder of radius rr along x0->x1. Assume
//x0<x1;
typedef struct {
        double x0;
	double x1;
        double y;
	double z;
	double R;
        double rr;
} TDUMBBELL_PARAMS;


typedef struct {
	int dim;
	double coeff;
	double epsilon;
} TEST_CURV_PARAMS;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

static	void test_grid_intfc(INTERFACE*);
static 	void print_volume_fraction(INTERFACE*);

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	TDUMBBELL_PARAMS d_params;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TEST_CURV_PARAMS curv_params; /* velocity function parameters */
	Locstate  sl;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 0.5; 	f_basic.U[2] = 0.5;
	f_basic.gmax[0] = 50;	f_basic.gmax[1] = 25; f_basic.gmax[2] = 25;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

	FT_StartUp(&front,&f_basic);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */
	    d_params.x0 = 0.25;
	    d_params.x1 = 0.75;
            d_params.y = 0.25;
            d_params.z = 0.25;
            d_params.R = 0.15;
            d_params.rr = 0.075;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&d_params;
	    level_func_pack.func = tdumbbell_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;

	    FT_InitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	curv_params.dim = 3;
	curv_params.coeff = 0.1;
	curv_params.epsilon = 0.0001;

	velo_func_pack.func_params = (POINTER)&curv_params;
	velo_func_pack.func = test_curvature_vel;

	FT_InitVeloFunc(&front,&velo_func_pack);

	/* For geometry-dependent velocity, use first 
	* order point propagation function, higher order
	* propagation requires surface propagate, currently
	* in writing, not yet in use. The following override
	* the assigned fourth_order_point_propagate.
	*/

	front._point_propagate = first_order_point_propagate;

	/* Propagate the front */

	add_to_debug("make_grid_intfc");
	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        double CFL;

	front->max_time = 0.003; 
	front->max_step = 1000;
	front->print_time_interval = 0.001;
	front->movie_frame_interval = 0.001;

        CFL = Time_step_factor(front) = 0.1;

	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
            FT_RedistMesh(front);
	    FT_ResetTime(front);

	    // Always output the initial interface.
	    FT_Save(front,out_name);
            FT_AddMovieFrame(front,out_name,binary);

	    // This is a virtual propagation to get maximum front 
	    // speed to determine the first time step.

	    FT_Propagate(front);
            FT_SetTimeStep(front);
	    FT_SetOutputCounter(front);
	}
	else
	{
            FT_SetOutputCounter(front);
	}

	FT_TimeControlFilter(front);

	for (;;)
        {
	    /* Propagating interface for time step dt */

	    FT_Propagate(front);
	    FT_AddTimeStepToCounter(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            FT_SetTimeStep(front);

	    /* Output section */

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (FT_IsSaveTime(front))
		FT_Save(front,out_name);
            if (FT_IsMovieFrameTime(front))
                FT_AddMovieFrame(front,out_name,binary);

            if (FT_TimeLimitReached(front))
                    break;

	    FT_TimeControlFilter(front);
	}
	set_binary_output(NO);
	test_grid_intfc(front->interf);
        (void) delete_interface(front->interf);
}       /* end test_propagate */

/********************************************************************
 *	Sample (dummbell 3D) level function for the initial interface    *
 ********************************************************************/

static double tdumbbell_func(
        POINTER func_params,
        double *coords)
{
        TDUMBBELL_PARAMS *d_params = (TDUMBBELL_PARAMS*)func_params;
	double x0,x1,y,z,f0,f1,R,rr;
	double arg;

        x0 = d_params->x0;
        y = d_params->y;
	z = d_params->z;
        x1  = d_params->x1;
	R = d_params->R;
        rr  = d_params->rr;

	f0 = x0+sqrt(R*R-rr*rr);
	f1 = x1-sqrt(R*R-rr*rr);

	if (coords[0]<f0)
	{
	    arg = sqrt(sqr(coords[0]-x0)+sqr(coords[1]-y)+sqr(coords[2]-z))-R;
        }
        else if(coords[0] > f1)
        {
            arg = sqrt(sqr(coords[0]-x1)+sqr(coords[1]-y)+sqr(coords[2]-z))-R;
        }
        else
        {
            arg = sqrt(sqr(coords[1]-y)+sqr(coords[2]-z))-rr;
        }
        return -arg;

}       /* end tdumbbell_func */

static int test_curvature_vel(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	TEST_CURV_PARAMS *curv_params = (TEST_CURV_PARAMS*) params;
        int i;
        double coeff,epsilon,eps;
        double kappa;
        double nor[MAXD];

        coeff = curv_params->coeff;
        epsilon = curv_params->epsilon;

        normal(p,hse,hs,nor,front);
        kappa = mean_curvature_at_point(p,hse,hs,front);

        for (i = 0; i < curv_params->dim; ++i)
        {
            vel[i] = nor[i]*(coeff - epsilon*kappa);
        }
}

static	void test_grid_intfc(
	INTERFACE *intfc)
{
	INTERFACE *grid_intfc;
	VOLUME_FRAC volume_frac;

	printf("Test make_grid_intfc on COMP_GRID:\n");
	grid_intfc = make_grid_intfc(intfc,COMP_GRID,NULL);
	free_grid_intfc(grid_intfc);

	printf("Test make_grid_intfc on DUAL_GRID:\n");
	grid_intfc = make_grid_intfc(intfc,DUAL_GRID,NULL);
	free_grid_intfc(grid_intfc);

	printf("Test make_grid_intfc on EXPANDED_DUAL_GRID:\n");
	grid_intfc = make_grid_intfc(intfc,EXPANDED_DUAL_GRID,NULL);
	free_grid_intfc(grid_intfc);
	
	add_to_debug("vol_frac");

	volume_frac.comp_vfrac = 1;
	printf("Doing for component 1\n");

	printf("Test make_grid_intfc on COMP_GRID with volume fraction:\n");
	grid_intfc = make_grid_intfc(intfc,COMP_GRID,&volume_frac);
	print_volume_fraction(grid_intfc);
	free_grid_intfc(grid_intfc);

	printf("Test make_grid_intfc on DUAL_GRID with volume fraction:\n");
	grid_intfc = make_grid_intfc(intfc,DUAL_GRID,&volume_frac);
	print_volume_fraction(grid_intfc);
	free_grid_intfc(grid_intfc);

	printf("Test make_grid_intfc on EXPANDED_DUAL_GRID with volume fraction:\n");
	grid_intfc = make_grid_intfc(intfc,EXPANDED_DUAL_GRID,&volume_frac);
	print_volume_fraction(grid_intfc);
	free_grid_intfc(grid_intfc);

	volume_frac.comp_vfrac = 2;
	printf("Doing for component 2\n");

	printf("Test make_grid_intfc on COMP_GRID with volume fraction:\n");
	grid_intfc = make_grid_intfc(intfc,COMP_GRID,&volume_frac);
	print_volume_fraction(grid_intfc);
	free_grid_intfc(grid_intfc);

	printf("Test make_grid_intfc on DUAL_GRID with volume fraction:\n");
	grid_intfc = make_grid_intfc(intfc,DUAL_GRID,&volume_frac);
	print_volume_fraction(grid_intfc);
	free_grid_intfc(grid_intfc);

	printf("Test make_grid_intfc on EXPANDED_DUAL_GRID with volume fraction:\n");
	grid_intfc = make_grid_intfc(intfc,EXPANDED_DUAL_GRID,&volume_frac);
	print_volume_fraction(grid_intfc);
	free_grid_intfc(grid_intfc);

}	/* test_grid_intfc */


static void print_volume_fraction(
	INTERFACE *grid_intfc)
{
	double ***area,***vol_frac;
	struct Table *T;
	int i,j,k,*gmax;
	RECT_GRID *grid = &topological_grid(grid_intfc);
	double area_norm;

	T = table_of_interface(grid_intfc);
	area = T->area;
	vol_frac = T->vol_frac;
	gmax = grid->gmax;
	area_norm = grid->h[0]*grid->h[1];

	printf("\n\nNormalized Block Area:\n\n");
	for (i = 0; i < gmax[0]; ++i)
	{
	    printf("ix = %d\n\n",i);
	    for (j = 0; j < gmax[1]; ++j)
	    {
	    	for (k = 0; k < gmax[2]; ++k)
		{
		    printf("%4.2f  ",area[i][j][k]/area_norm);
		}
		printf("\n");
	    }
	    printf("\n\n");
	}
	printf("\n\nBlock Volume Fraction:\n\n");
	for (i = 0; i < gmax[0]; ++i)
	{
	    printf("ix = %d\n\n",i);
	    for (j = 0; j < gmax[1]; ++j)
	    {
	    	for (k = 0; k < gmax[2]; ++k)
		{
		    printf("%4.2f  ",vol_frac[i][j][k]);
		}
		printf("\n");
	    }
	    printf("\n\n");
	}
}	/* end print_volume_fraction */
