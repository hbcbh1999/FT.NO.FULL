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
*				test_map_3d.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This is example of three circles all moving the a normal velocity.
*	Bifurcation occurs when they meet each other. FronTier solves
*	the bifurcation automatically.
*
*	This example shows how to convert the interface into arrays
*	of floats. Functions introduced in this example include:
*
*	Dimension(intfc);
*
* int     NumOfPoints(SURFACE *s);
* int     NumOfPoints(INTERFACE *intfc);
* int     NumOfNodes(INTERFACE *intfc);
* int     NumOfBonds(INTERFACE *intfc);
* int     NumOfCurves(INTERFACE *intfc);
* int     NumOfSurfaces(INTERFACE *intfc);
* int     NumOfTris(SURFACE *s);
* int     NumOfTris(INTERFACE *intfc);

* void    ArrayOfPoints(INTERFACE *intfc, double *coords);
* void    ArrayOfTri(SURFACE *surface, TRI **tris);
* void    ArrayOfTri(SURFACE *surface, double *coords, int *vertices_index);
* void    ArrayOfTri(INTERFACE *intfc, TRI **tris);
* void    ArrayOfTri(INTERFACE *intfc, double *coords, int *vertex_indices);
*
*	See their use in function map_output_interface().
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static int norm_vel_func(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);
static void map_output_interface(Front*,char*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
        int num_cir;
        double **cen;
        double *rad;
} TMC_PARAMS;


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

struct _TNORV_PARAMS
{
        int dim;
        double coeff;
};
typedef struct _TNORV_PARAMS TNORV_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TNORV_PARAMS norm_params; /* velocity function parameters */
	TMC_PARAMS mc_params;
	Locstate  sl;

	f_basic.dim = 3;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 40;	f_basic.gmax[1] = 40;	f_basic.gmax[2] = 40;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
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

	    mc_params.num_cir = 1;
            FT_VectorMemoryAlloc((POINTER*)&mc_params.rad,mc_params.num_cir,FLOAT);
            FT_MatrixMemoryAlloc((POINTER*)&mc_params.cen,mc_params.num_cir,3,FLOAT);
	    mc_params.cen[0][0] = 0.4;
	    mc_params.cen[0][1] = 0.4;
	    mc_params.cen[0][2] = 0.4;
	    mc_params.rad[0] = 0.3;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&mc_params;
	    level_func_pack.func = multi_circle_func;
	    level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;
	
	    FT_InitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	norm_params.dim = 3;
	norm_params.coeff = 0.1;

	velo_func_pack.func_params = (POINTER)&norm_params;
	velo_func_pack.func = norm_vel_func;
	
	FT_InitVeloFunc(&front,&velo_func_pack);

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        double CFL;

	front->max_time = 0.20;
	front->max_step = 100;
	front->print_time_interval = 0.10;
	front->movie_frame_interval = 0.02;

        CFL = Time_step_factor(front);

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
        (void) delete_interface(front->interf);
}       /* end test_propagate */

static int norm_vel_func(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        TNORV_PARAMS *norv_params;
        int i;
        double coeff;
        double curvature;
        double nor[MAXD];
                                                                                
        norv_params = (TNORV_PARAMS*)params;
        coeff = norv_params->coeff;
                                                                                
        GetFrontNormal(p,hse,hs,nor,front);
        for (i = 0; i < norv_params->dim; ++i)
        {
            vel[i] = nor[i]*coeff;
	    //vel[i] = (i==0?1:0) * coeff;
        }
}       /* end normal_vel_func */

static void SaveAsTecplot(
	char *filename, 
	int nP, 
	double *coords, 
	int nTri, 
	int *vIndices)
{
	printf("open file: %s\n",filename);
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
	    printf("\n SaveAsTecplot: can't open %s for writing.",filename);
	    return;
	}
	fprintf(hfile, "TITLE = %s \n", filename);
	fprintf(hfile, "VARIABLES = \"X\" \"Y\" \"Z\" \n");
	fprintf(hfile, "ZONE N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON \n", 
					nP, nTri);
	for(int i=0; i<nP; i++)
	    fprintf(hfile,"%f %f %f \n",coords[i*3+0],coords[i*3+1],
	    				coords[i*3+2]);
	for(int i=0; i<nTri; i++)
	    fprintf(hfile, "%d %d %d %d", vIndices[i*3+0], vIndices[i*3+1],
					vIndices[i*3+2], vIndices[i*3+0]);
	fclose(hfile);
		
}


static void map_output_interface(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
        char filename[100];
	int i,step;
	step = front->step;
			
	int num_points = NumOfIntfcPoints(intfc);
	int num_tris = NumOfIntfcTris(intfc);
	int dim = 3;
	double *coords = (double*)malloc(num_points*dim*sizeof(double));
	int *vertex_indices = (int*)malloc(num_tris*3*sizeof(int));
	
	// the interface
	printf("Number of interface points = %d\n",num_points);
	printf("Number of interface triangles = %d\n",num_tris);
	sprintf(filename,"%s-%d.intfc.plt",out_name,step);
	ArrayOfIntfcTris(intfc,coords,vertex_indices);	
	SaveAsTecplot(filename,num_points,coords,num_tris,vertex_indices);
	
	int num_surfaces = NumOfSurfaces(intfc);
	printf("Number of surfaces = %d\n",num_surfaces);
	SURFACE **surfaces;
	surfaces = (SURFACE**)malloc(num_surfaces*sizeof(SURFACE*));
	
	ArrayOfSurfaces(intfc, surfaces);
	// print the tris of each surface
	for(i = 0; i < num_surfaces; i++)
	{
	    sprintf(filename,"%s-%d.surface-%d.plt",out_name,step,i);
	    num_points = NumOfSurfPoints(surfaces[i]);
	    num_tris   = NumOfSurfTris(surfaces[i]);
	    printf("Number of points on surface %d = %d\n",i+1,num_points);
	    printf("Number of tris on surface %d = %d\n",i+1,num_tris);
	    ArrayOfSurfTris(surfaces[i],coords,vertex_indices);
	    SaveAsTecplot(filename,num_points,coords,num_tris,vertex_indices);
	}
}	/* end map_output */
