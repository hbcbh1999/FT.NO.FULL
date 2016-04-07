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

#include <iFluid.h>
#include <crystal.h>
#include "subsurf.h"

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/
#define		MAX_NUM_VERTEX_IN_CELL		20

	/*  Local Application Function Declarations */

static void 	subsurf_main_driver(Front*,SEED_PARAMS,C_CARTESIAN&,
			Incompress_Solver_Smooth_Basis*);
static void	read_crystal_params(char*,CRT_PARAMS*);
static void	solute_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	record_1d_growth(Front*,double***,int*);
static void 	plot_growth_data(char*,double**,int);
static void 	read_seed_params(int,char*,SEED_PARAMS*);
static double 	solute_func(double*,COMPONENT,double);
static boolean 	fractal_dimension(Front*,SEED_PARAMS,double*,double*);

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
	SEED_PARAMS s_params;
	IF_PARAMS iFparams;
	CRT_PARAMS cRparams;
	PROB_TYPE prob_type;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

	/* Construct the c_cartesian and l_cartesian class for crystal and
	 * incompressible fluid computation.
        */
	
	//Initialize Petsc
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	//if (debugging("trace")) printf("Passed PetscInitialize()\n");


	C_CARTESIAN  c_cartesian(front);
	Incompress_Solver_Smooth_Basis *l_cartesian;
	if(f_basic.dim == 2)
	    l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	else if(f_basic.dim == 3)
	    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

	/* Initialize basic computational data */

	in_name      		= f_basic.in_name;
	restart_state_name      = f_basic.restart_state_name;
        out_name     		= f_basic.out_name;
        restart_name 		= f_basic.restart_name;
        RestartRun   		= f_basic.RestartRun;
        ReadFromInput   	= f_basic.ReadFromInput;
	RestartStep 		= f_basic.RestartStep;

	sprintf(restart_state_name,"%s/state.ts%s",restart_name,
			right_flush(RestartStep,7));
	sprintf(restart_name,"%s/intfc-ts%s",
			restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-nd%s",restart_name,
				right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
				right_flush(pp_mynode(),4));
	}
#endif /* defined(__MPI__) */
	if (!ReadFromInput)
	{
	    (void) printf("ERROR: Input file needed!\n");
	    clean_up(ERROR);
	}
	FT_ReadSpaceDomain(in_name,&f_basic);

	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

	fflush(stdout);

	iFparams.dim = f_basic.dim;
	read_iF_movie_options(in_name,&iFparams);
	read_crt_movie_options(in_name,&cRparams);
        if (debugging("trace")) printf("Passed read_ifluid_params()\n");

	if (!RestartRun)
	{
	    if (f_basic.dim == 1)
	    {
		(void) printf("Subsurface problem has no 1D version!\n");
		clean_up(ERROR);
	    }

	    /* Initialize interface through level function */

	    read_seed_params(f_basic.dim,in_name,&s_params);
	    level_func_pack.neg_component = CRYSTAL_COMP;
	    level_func_pack.pos_component = SOLUTE_COMP;
	    level_func_pack.wave_type = GROWING_BODY_BOUNDARY;
	    level_func_pack.func_params = (POINTER)&s_params;
	    level_func_pack.func = crystal_seed_curve;
	    if (f_basic.dim == 3)
                    level_func_pack.set_3d_bdry = YES;

	    if (debugging("trace")) printf("Passed setting intfc params()\n");
	    FT_InitIntfc(&front,&level_func_pack);
	    if (debugging("trace")) printf("Passed FT_InitIntfc()\n");
		read_ss_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (debugging("trace")) 
	    	printf("Passed read_ss_dirichlet_bdry_data()\n");
	    if (f_basic.dim != 3)
	    	FT_ClipIntfcToSubdomain(&front);
	    if (debugging("trace")) 
	    	printf("Passed FT_ClipIntfcToSubdomain()\n");
	}

	read_crystal_params(in_name,&cRparams);
	if (debugging("trace")) printf("Passed read_crystal_params()\n");
	read_fluid_params(in_name,&iFparams);
	if (debugging("trace")) printf("Passed read_fluid_params()\n");
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	front.extra1 = (POINTER)&iFparams;
	front.extra2 = (POINTER)&cRparams;

	iFparams.dim = f_basic.dim;
	cRparams.dim = f_basic.dim;

	velo_func_pack.func_params = (POINTER)&iFparams;
	velo_func_pack.func = NULL;

	FT_InitVeloFunc(&front,&velo_func_pack);
	if (debugging("trace")) printf("Passed FT_InitVeloFunc()\n");

	if (debugging("trace")) printf("Before initMesh()\n");
        c_cartesian.initMesh();
	if (debugging("trace")) printf("Passed Crystal initMesh()\n");
        l_cartesian->initMesh();
	if (debugging("trace")) printf("Passed iFluid initMesh()\n");
	if (debugging("sample_velocity"))
            l_cartesian->initSampleVelocity(in_name);

        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */

        front._point_propagate = solute_point_propagate;

	if (RestartRun)
	{
	    c_cartesian.readFrontInteriorStates(restart_state_name);
	    FT_FreeGridIntfc(&front);
	    l_cartesian->readFrontInteriorStates(restart_state_name);
	}
	else
	{
	    initFrontStates(&front);
	    c_cartesian.setInitialCondition(s_params);
	    if (debugging("trace")) 
		printf("Passed Crystal setInitialCondition()\n");
	    FT_FreeGridIntfc(&front);
	    init_fluid_state_func(l_cartesian,prob_type);
	    l_cartesian->setInitialCondition();
	    if (debugging("trace")) 
		printf("Passed iFluid setInitialCondition()\n");
	}
	cRparams.field->vel = iFparams.field->vel;

	/* Propagate the front */

	subsurf_main_driver(&front,s_params,c_cartesian,l_cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void subsurf_main_driver(
        Front *front,
	SEED_PARAMS s_params,
	C_CARTESIAN &c_cartesian,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        int status;
        Front *newfront;
        double dt,dt_frac,CFL;
        char s[10];
        double fcrds[MAXD];
        int  i,dim = front->rect_grid->dim;
	INTERFACE *grid_intfc;
	IF_PARAMS *iFparams;
	CRT_PARAMS *cRparams;
	double solute_time_limit; 
	double frac_dim,radius;
	FILE *Radius_file,*FracDim_file;
	boolean bdry_reached = NO;
	double **growth_data = NULL;
	int count = 0;

	if (debugging("trace"))
	    printf("Entering subsurf_main_driver()\n");

	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

	iFparams = (IF_PARAMS*)front->extra1;
	cRparams = (CRT_PARAMS*)front->extra2;

	front->hdf_movie_var = NULL;

	if (dim == 1)
	    c_cartesian.oneDimPlot(out_name);

	if (dim == 1)
	{
	    bdry_reached = NO;
	    record_1d_growth(front,&growth_data,&count);
	}
	else
	    bdry_reached = fractal_dimension(front,s_params,&frac_dim,
						&radius);
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

		/*TMP*/
		// Front standard output
                FT_Save(front,out_name);
		// Problem specific output
		c_cartesian.printFrontInteriorStates(out_name);
	    	if (debugging("trace")) printf("Passed Crystal "
				"printFrontInteriorStates()\n");
		l_cartesian->printFrontInteriorStates(out_name);
	    	if (debugging("trace")) printf("Passed iFluid "
				"printFrontInteriorStates()\n");
	if (!RestartRun)
	{
	    FT_ResetTime(front);
	    if (debugging("trace"))
		printf("Before FT_Propagate() front->dt = %f\n",front->dt);
	    FT_Propagate(front);
	    if (debugging("trace")) printf("Calling Cystal solve()\n");
	    c_cartesian.solve(front->dt);
	    if (debugging("trace")) printf("Calling iFluid solve()\n");
	    l_cartesian->solve(front->dt);
	    FT_SetTimeStep(front);
	    c_cartesian.setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);
	    l_cartesian->setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
        }
	if (debugging("trace"))
	    printf("Passed second restart check()\n");

	FT_SetOutputCounter(front);
	FT_TimeControlFilter(front);

	if (debugging("step_size"))
		printf("Time step from start: %f\n",front->dt);
        for (;;)
        {
	    if (debugging("trace"))
		printf("Before FT_Propagate()\n");

	    FT_Propagate(front);
            if (debugging("trace")) printf("Passed FT_Propagate()\n");

	    if (debugging("trace")) printf("Calling Cystal solve()\n");
	    c_cartesian.solve(front->dt);
	    if (debugging("trace")) printf("Passed Cystal solve()\n");
	    if (debugging("trace")) printf("Calling iFluid solve()\n");
	    l_cartesian->solve(front->dt);
	    if (debugging("trace")) printf("Passed iFluid solve()\n");
	    if (debugging("trace"))
            {
                printf("After solve()\n");
                print_storage("at end of time step","trace");
            }

	    FT_AddTimeStepToCounter(front);
	    FT_SetTimeStep(front);

	    if (debugging("step_size"))
		printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);
	    if (debugging("step_size"))
		printf("Time step from c_cartesian.max_dt(): %f\n",front->dt);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
	    if (debugging("step_size"))
		printf("Time step from l_cartesian->max_dt(): %f\n",front->dt);

            printf("\ntime = %f   step = %7d   dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (FT_IsSaveTime(front))
	    {
		// Front standard output
                FT_Save(front,out_name);
		// Problem specific output
		c_cartesian.printFrontInteriorStates(out_name);
	    	if (debugging("trace")) printf("Passed Crystal "
				"printFrontInteriorStates()\n");
		l_cartesian->printFrontInteriorStates(out_name);
	    	if (debugging("trace")) printf("Passed iFluid "
				"printFrontInteriorStates()\n");
	    }
	    if (debugging("trace"))
		printf("After print output()\n");
            if (FT_IsMovieFrameTime(front))
	    {
		// Front standard output
		if (dim != 1)
		{
		    if (debugging("trace"))
                    	printf("Calling c_cartesian.initMovieVariable()\n");
                    c_cartesian.initMovieVariables();
		    if (debugging("trace"))
                    	printf("Calling l_cartesian->augmentMovieVariable()\n");
                    l_cartesian->augmentMovieVariables();
		    if (debugging("trace"))
                    	printf("Calling FT_AddMovieFrame()\n");
                    FT_AddMovieFrame(front,out_name,binary);
		}
	    	else
	    	    c_cartesian.oneDimPlot(out_name);
	    }
	    if (debugging("trace"))
		printf("After make movie frame()\n");
	    if (dim == 1)
	    {
	    	bdry_reached = NO;
		record_1d_growth(front,&growth_data,&count);
	    }
	    else
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
					&radius);
	    if (debugging("trace"))
		printf("After auxilary output()\n");

	    if (bdry_reached) front->time_limit_reached = YES;
	    if (pp_mynode() == 0 && dim != 1 && !s_params.grow_from_floor)
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
	    if (debugging("trace"))
		printf("After time output()\n");
	    /* Output section, next dt may be modified */

	    FT_TimeControlFilter(front);

	    if (debugging("step_size"))
		printf("Time step from FrontOutputTimeControl(): %f\n",
					front->dt);
        }
	if (pp_mynode() == 0 && dim != 1 && !s_params.grow_from_floor)
	{
	    fclose(Radius_file);
	    fclose(FracDim_file);
	}
}       /* end subsurf_main_driver */

enum _PROXIMITY {
	FLOOR = 1,
	CEILING,
	SPACE
};
typedef enum _PROXIMITY PROXIMITY;

extern  double getStateSolute(
        POINTER state)
{
        STATE *solute_state = (STATE*)state;
        return solute_state->solute;
}       /* end getStateSolute */

static	double solute_func(
	double *coords,
	COMPONENT comp,
	double C_0)
{
	if (comp != SOLUTE_COMP) return 0.0;
	return C_0;
}	/* end solute_init_func */

static  void solute_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	ifluid_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
	crystal_point_propagate(front,wave,oldp,newp,oldhse,oldhs,dt,V);
}	/* end solute_point_propagate */

static	void	read_crystal_params(
	char *in_name,
	CRT_PARAMS *cRparams)
{
	FILE *infile;
	char scheme[200];

	infile = fopen(in_name,"r");
	CursorAfterString(infile,"Diffusion coefficient:");
	fscanf(infile,"%lf",&cRparams->D);
	(void) printf("%f\n",cRparams->D);
	CursorAfterString(infile,"Growth rate:");
	fscanf(infile,"%lf",&cRparams->k);
	(void) printf("%f\n",cRparams->k);
	CursorAfterString(infile,"Equilibrium concentration:");
	fscanf(infile,"%lf",&cRparams->C_eq);
	(void) printf("%f\n",cRparams->C_eq);
	CursorAfterString(infile,"Ambient concentration:");
	fscanf(infile,"%lf",&cRparams->C_0);
	(void) printf("%f\n",cRparams->C_0);
	CursorAfterString(infile,"Crystal density:");
	fscanf(infile,"%lf",&cRparams->rho_s);
	(void) printf("%f\n",cRparams->rho_s);
	CursorAfterString(infile,"Choose numerical scheme");
	(void) printf("\n");
	CursorAfterString(infile,"Enter scheme:");
	fscanf(infile,"%s",scheme);
	(void) printf("%s\n",scheme);
	if ((scheme[0] == 'E' || scheme[0] == 'e') &&
	    (scheme[1] == 'X' || scheme[1] == 'x')) 
	    cRparams->num_scheme = UNSPLIT_EXPLICIT;
	else if ((scheme[0] == 'I' || scheme[0] == 'i') &&
	    (scheme[1] == 'M' || scheme[1] == 'm')) 
	    cRparams->num_scheme = UNSPLIT_IMPLICIT;
	else if ((scheme[0] == 'C' || scheme[0] == 'c') &&
	    (scheme[1] == 'N' || scheme[1] == 'n')) 
	    cRparams->num_scheme = CRANK_NICOLSON;
	else
	{
	    printf("Unknown numerical scheme!\n");
	    clean_up(ERROR);
	}
	CursorAfterString(infile,"Choose point propagation scheme");
	(void) printf("\n");
        CursorAfterString(infile,"Enter scheme:");
        fscanf(infile,"%s",scheme);
	(void) printf("%s\n",scheme);
        if ((scheme[0] == 'E' || scheme[0] == 'e') &&
            (scheme[1] == 'X' || scheme[1] == 'x'))
            cRparams->point_prop_scheme = EXPLICIT_EULER;
        else if ((scheme[0] == 'I' || scheme[0] == 'i') &&
            (scheme[1] == 'M' || scheme[1] == 'm'))
            cRparams->point_prop_scheme = IMPLICIT_EULER;
        else if ((scheme[0] == 'M' || scheme[0] == 'm') &&
            (scheme[1] == 'P' || scheme[1] == 'p'))
            cRparams->point_prop_scheme = MIDDLE_POINT;
        else
        {
            printf("Unknown point propagation scheme!\n");
            clean_up(ERROR);
        }
	CursorAfterString(infile,"Enter yes to add curvature effect:");
        fscanf(infile,"%s",scheme);
        (void) printf("%s\n",scheme);
        if (scheme[0] == 'Y' || scheme[0] == 'y')
            cRparams->add_curvature = YES;
        else
            cRparams->add_curvature = NO;
	fclose(infile);
}

static boolean fractal_dimension(
	Front *front,
	SEED_PARAMS s_params,
	double *frac_dim,
	double *radius)
{
	double coords[MAXD],*center = s_params.space_center[0];
	double dist,r_sqr,r_max,r_min = s_params.seed_radius;
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
	boolean bdry_reached = NO;
	double margin[MAXD];

	if (s_params.grow_from_floor || s_params.grow_from_ceiling)
	    return bdry_reached;
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
		    if (comp == CRYSTAL_COMP)
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
		    if (comp == CRYSTAL_COMP)
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
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
	static int total_length = 0;
	STATE *sl,*sr;
	
	if (*count >= total_length)
	{
	    static double **tmp_data;
	    total_length += 1000;
	    FT_MatrixMemoryAlloc((POINTER*)&tmp_data,total_length,3,sizeof(double));
	    printf("Memory allocated\n");
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
	    if (negative_component(hs) == SOLUTE_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = sl->solute;
	    }
	    else if (positive_component(hs) == SOLUTE_COMP)
	    {
		(*growth_data)[*count][1] = Coords(p)[0];
		(*growth_data)[*count][2] = sr->solute;
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

	sprintf(fname,"%s-posn.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface position vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][1]);
	fclose(ofile);

	sprintf(fname,"%s-solt.xg",out_name);
	ofile = fopen(fname,"w");
	fprintf(ofile,"\"Interface solute vs. time\"\n");
	for (i = 0; i < count; ++i)
	    fprintf(ofile,"%f  %f\n",growth_data[i][0],growth_data[i][2]);
	fclose(ofile);
}	/* end plot_growth_data */

static	void read_seed_params(
	int dim,
	char *in_name,
	SEED_PARAMS *s_params)
{
	FILE *infile;
	char s[100];
	int i,j;

	infile = fopen(in_name,"r");	
	s_params->dim = dim;

	CursorAfterString(infile,"Are there seeds grow from floor:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    s_params->grow_from_floor = YES;
	else
	    s_params->grow_from_floor = NO;
	if (s_params->grow_from_floor == YES)
	{
	    CursorAfterString(infile,"Enter floor level:");
	    fscanf(infile,"%lf",&s_params->floor_level);
	    (void) printf("%f\n",s_params->floor_level);
	    CursorAfterString(infile,"Enter number of floor seeds:");
	    fscanf(infile,"%d",&s_params->num_floor_seeds);
	    (void) printf("%d\n",s_params->num_floor_seeds);
	    FT_MatrixMemoryAlloc((POINTER*)&s_params->floor_center,
			s_params->num_floor_seeds,MAXD,sizeof(double));
	    for (i = 0; i < s_params->num_floor_seeds; ++i)
	    {
	    	sprintf(s,"Enter center coordinates of floor seed %d:",i+1);
	    	CursorAfterString(infile,s);
		for (j = 0; j < dim-1; ++j)
		{
		    fscanf(infile,"%lf ",&s_params->floor_center[i][j]);
	    	    (void) printf("%f ",s_params->floor_center[i][j]);
		}
		(void) printf("\n");
		s_params->floor_center[i][dim-1] = s_params->floor_level;
	    }
	}

	CursorAfterString(infile,"Are there seeds grow from ceiling:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    s_params->grow_from_ceiling = YES;
	else
	    s_params->grow_from_ceiling = NO;
	if (s_params->grow_from_ceiling == YES)
	{
	    CursorAfterString(infile,"Enter ceiling level:");
	    fscanf(infile,"%lf",&s_params->ceiling_level);
	    (void) printf("%f\n",s_params->ceiling_level);
	    CursorAfterString(infile,"Enter number of ceiling seeds:");
	    fscanf(infile,"%d",&s_params->num_ceiling_seeds);
	    (void) printf("%d\n",s_params->num_ceiling_seeds);
	    FT_MatrixMemoryAlloc((POINTER*)&s_params->ceiling_center,
			s_params->num_ceiling_seeds,MAXD,sizeof(double));
	    for (i = 0; i < s_params->num_ceiling_seeds; ++i)
	    {
	    	sprintf(s,"Enter center coordinates of ceiling seed %d:",i+1);
	    	CursorAfterString(infile,s);
		for (j = 0; j < dim-1; ++j)
		{
		    fscanf(infile,"%lf ",&s_params->ceiling_center[i][j]);
		    (void) printf("%f ",s_params->ceiling_center[i][j]);
		}
		(void) printf("\n");
		s_params->ceiling_center[i][dim-1] = s_params->ceiling_level;
	    }
	}

	CursorAfterString(infile,"Are there seeds grow from space:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	    s_params->grow_from_space = YES;
	else
	    s_params->grow_from_space = NO;
	if (s_params->grow_from_space == YES)
	{
	    CursorAfterString(infile,"Enter number of space seeds:");
	    fscanf(infile,"%d",&s_params->num_space_seeds);
	    (void) printf("%d\n",s_params->num_space_seeds);
	    FT_MatrixMemoryAlloc((POINTER*)&s_params->space_center,
			s_params->num_space_seeds,MAXD,sizeof(double));
	    for (i = 0; i < s_params->num_space_seeds; ++i)
	    {
	    	sprintf(s,"Enter center coordinates of space seed %d:",i+1);
	    	CursorAfterString(infile,s);
		for (j = 0; j < dim; ++j)
		{
		    fscanf(infile,"%lf ",&s_params->space_center[i][j]);
	    	    (void) printf("%f ",s_params->space_center[i][j]);
		}
		(void) printf("\n");
	    }
	}

	CursorAfterString(infile,"Enter radius of seeds:");
	fscanf(infile,"%lf",&s_params->seed_radius);
	(void) printf("%f\n",s_params->seed_radius);
}	/* end read_seed_params */
