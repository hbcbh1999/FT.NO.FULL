
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

#include "crystal.h"
#include "crystal_basic.h"

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/
#define		MAX_NUM_VERTEX_IN_CELL		20

	/*  Local Application Function Declarations */

static void 	solute_main_driver(Front*,SEED_PARAMS,C_CARTESIAN&);
static void	read_crystal_params(char*,CRT_PARAMS*);
static void	crystal_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	neumann_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	reaction_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void	record_1d_growth(Front*,double***,int*);
static void	record_1d_growth(Front*,double***,int*);
static void	record_1d_growth(Front*,double***,int*);
static void	record_1d_growth(Front*,double***,int*);
static void 	plot_growth_data(char*,double**,int);
static void 	read_seed_params(int,char*,SEED_PARAMS*);
static void 	read_crt_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);
static void 	initFrontStates(Front*);
static void 	readFrontStates(Front*,char*);
static void 	euler_forward_scheme(CRT_PARAMS*,double,double,double,double,
				double,double*);
static void 	euler_backward_scheme(CRT_PARAMS*,double,double,double,double,
				double,double*);
static void 	middle_point_scheme(CRT_PARAMS*,double,double,double,double,
				double,double*);
static void 	constant_state(CRT_PARAMS*,double,double,double,double,
				double,double*);
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
	CRT_PARAMS cRparams;

	C_CARTESIAN       c_cartesian(front);

	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.size_of_intfc_state = sizeof(STATE);
	
	//Initialize Petsc
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	//if (debugging("trace")) printf("Passed PetscInitialize()\n");

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


	if (f_basic.dim >= 2)
	    read_seed_params(f_basic.dim,in_name,&s_params);
	cRparams.dim = f_basic.dim;
	read_crt_movie_options(in_name,&cRparams);
        if (debugging("trace")) printf("Passed read_iFparams()\n");

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    level_func_pack.neg_component = CRYSTAL_COMP;
	    level_func_pack.pos_component = SOLUTE_COMP;
	    level_func_pack.wave_type = GROWING_BODY_BOUNDARY;
	    if (f_basic.dim == 1)
	    {
	    	double **point;
                FT_MatrixMemoryAlloc((POINTER*)&point,1,1,FLOAT);
                s_params.point = point[0][0] = 0.03;
                level_func_pack.num_points = 1;
                level_func_pack.point_array = point;
	    }
	    else if (f_basic.dim >= 2)
	    {
		if (s_params.num_floor_seeds == 0 &&
		    s_params.num_ceiling_seeds == 0 &&
		    s_params.num_space_seeds == 0)
		{
		    level_func_pack.func_params = NULL;
		    level_func_pack.func = NULL;
		}
		else
		{
		    level_func_pack.func_params = (POINTER)&s_params;
		    level_func_pack.func = crystal_seed_curve;
		}
		if (f_basic.dim == 3)
		    level_func_pack.set_3d_bdry = YES;
	    }

	    if (debugging("trace")) printf("Passed setting intfc params()\n");
	    FT_InitIntfc(&front,&level_func_pack);
	}
	if (debugging("trace")) printf("Passed FT_InitIntfc()\n");
	read_crt_dirichlet_bdry_data(in_name,&front,f_basic);
	if (debugging("trace")) 
	    printf("Passed read_crt_dirichlet_bdry_data()\n");
	if (f_basic.dim == 2)
	    FT_ClipIntfcToSubdomain(&front);
	if (debugging("trace"))
	{
	    if (f_basic.dim == 2)
	    {
		char xg_name[100];
		sprintf(xg_name,"init_intfc-%d",pp_mynode());
		xgraph_2d_intfc(xg_name,front.interf);
	    }
	    else if (f_basic.dim == 3)
	    {
		gview_plot_color_interface("init_intfc",front.interf,YES);
	    }
	}

	read_crystal_params(in_name,&cRparams);
        if (debugging("trace")) printf("Passed read_crystal_params()\n");
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	front.extra2 = (POINTER)&cRparams;

	cRparams.dim = f_basic.dim;

	velo_func_pack.func_params = (POINTER)&cRparams;
	velo_func_pack.func = NULL;

	FT_InitVeloFunc(&front,&velo_func_pack);
	if (debugging("trace")) printf("Passed FT_InitVeloFunc()\n");

	if (debugging("trace")) printf("Before initMesh()\n");
        c_cartesian.initMesh();
	if (debugging("trace")) printf("Passed Crystal initMesh()\n");

        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */

        front._point_propagate = crystal_point_propagate;

	if (RestartRun)
	{
	    c_cartesian.readFrontInteriorStates(restart_state_name);
	}
	else
	{
	    initFrontStates(&front);
	    c_cartesian.setInitialCondition(s_params);
	}
	cRparams.field->vel = NULL;

	/* Propagate the front */

	solute_main_driver(&front,s_params,c_cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void solute_main_driver(
        Front *front,
	SEED_PARAMS s_params,
	C_CARTESIAN &c_cartesian)
{
        int status;
        double dt,CFL;
        char s[10];
        double fcrds[MAXD];
        int  i,dim = front->rect_grid->dim;
	INTERFACE *grid_intfc;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
	double h_min;
	double solute_time_limit; 
	double frac_dim,radius;
	FILE *Radius_file,*FracDim_file,*FracDim_vs_Damkl;
	boolean bdry_reached = NO;
	double **growth_data = NULL;
	int count = 0;
	double D = cRparams->D;
	double k = cRparams->k;

	if (debugging("trace"))
	    printf("Entering solute_main_driver()\n");
	Curve_redistribution_function(front) = expansion_redistribute;
        CFL = Time_step_factor(front);

	h_min = front->rect_grid->h[0];
	for (i = 1; i < dim; ++i)
	    h_min = std::min(h_min,front->rect_grid->h[i]);

	front->hdf_movie_var = NULL;
        c_cartesian.initMovieVariables();
	if (debugging("trace"))
	    printf("Passed initMovieVariables()\n");

	if (dim != 1)
            FT_AddMovieFrame(front,out_name,binary);
        else
            c_cartesian.oneDimPlot(out_name);
	if (debugging("trace"))
	    printf("Passed FT_AddMovieFrame()\n");

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
	    printf("Radius_file = %d\n",Radius_file);
	    sprintf(fname,"%s/FracDim",out_name);
	    FracDim_file = fopen(fname,"w");
	    sprintf(fname,"%s/DimDam",out_name);
	    FracDim_vs_Damkl = fopen(fname,"w");

	    fprintf(Radius_file,"\"Crystal radius\"\n\n");
	    fprintf(Radius_file,"%f  %f\n",front->time,radius);
	    fprintf(FracDim_file,"\"Fractal dimension\"\n\n");
	    fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
	    fprintf(FracDim_vs_Damkl,"\"Fractal dimension vs Damkhl\"\n\n");
	    fprintf(FracDim_vs_Damkl,"%f  %f\n",k*radius/D,frac_dim);
	    fflush(Radius_file);
	    fflush(FracDim_file);
	    fflush(FracDim_vs_Damkl);
	}

        if (!RestartRun)
        {
	    if (debugging("sample"))
	    {
		cRparams->max_solute = -HUGE;	
		cRparams->min_solute =  HUGE;	
	    }
	    FT_ResetTime(front);
	    FT_Propagate(front);
	    c_cartesian.solve(front->dt);
	    if (debugging("sample"))
	    {
		printf("Front: max_solute = %f  min_solute = %f\n",
			cRparams->max_solute,cRparams->min_solute);
	    }
	    FT_SetTimeStep(front);
	    c_cartesian.setAdvectionDt();
	    front->dt = std::min(front->dt,CFL*c_cartesian.max_dt);

	    if (dim == 1)
	    	c_cartesian.oneDimPlot(out_name);
	    FT_SetOutputCounter(front);

        }
        else
        {
	    FT_SetOutputCounter(front);
        }
	if (debugging("trace"))
	    printf("Passed second restart check()\n");

	FT_TimeControlFilter(front);

        for (;;)
        {
	    if (debugging("trace"))
		printf("Before FT_Propagate()\n");
	    if (debugging("sample"))
	    {
		cRparams->max_solute = -HUGE;	
		cRparams->min_solute =  HUGE;	
	    }

	    xgraph_2d_intfc("test4.xg",front->interf);
	    FT_Propagate(front);
	    xgraph_2d_intfc("test5.xg",front->interf);
	    if (debugging("trace")) printf("Passed FT_Propagate()\n");

	    if (debugging("sample"))
	    {
		printf("Front: max_solute = %f  min_solute = %f\n",
			cRparams->max_solute,cRparams->min_solute);
	    }
	    if (debugging("trace")) printf("Calling Cystal solve()\n");
	    c_cartesian.solve(front->dt);
	    if (debugging("trace")) printf("Passed Cystal solve()\n");

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

            (void) printf("\ntime = %f   step = %7d   dt = %f\n",
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
	    }
	    if (debugging("trace"))
		printf("After print output()\n");
            if (FT_IsMovieFrameTime(front))
	    {
		// Front standard output
		if (dim != 1)
		{
		    c_cartesian.initMovieVariables();
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
		fprintf(FracDim_vs_Damkl,"%f  %f\n",k*radius/D,frac_dim);
		fflush(Radius_file);
		fflush(FracDim_file);
		fflush(FracDim_vs_Damkl);
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
        }
	if (pp_mynode() == 0 && dim != 1 && !s_params.grow_from_floor)
	{
	    fclose(Radius_file);
	    fclose(FracDim_file);
	    fclose(FracDim_vs_Damkl);
	}
}       /* end solute_main_driver */

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

static  void crystal_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	switch(wave_type(oldhs))
	{
	case SUBDOMAIN_BOUNDARY:
            return;
	case NEUMANN_BOUNDARY:
	    neumann_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	case DIRICHLET_BOUNDARY:
	    dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	case GROWING_BODY_BOUNDARY:
	    reaction_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
	    return;
	}
}	/* end crystal_point_propagate */

static  void neumann_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        int i, dim = front->rect_grid->dim;
	double nor[MAXD];
        double p1[MAXD],p2[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
        double s0,s1,s2;
        STATE *sl,*sr,*state;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double *solute = cRparams->field->solute;

        for (i = 0; i < dim; ++i)
            Coords(newp)[i] = Coords(oldp)[i];
       	FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
       	state = (negative_component(oldhs) == SOLUTE_COMP) ? sl : sr;
	s0 = state->solute;
       	FT_NormalAtPoint(oldp,front,nor,SOLUTE_COMP);
       	dn = grid_size_in_direction(nor,h,dim);
       	for (i = 0; i < dim; ++i)
       	    p1[i] = p0[i] + nor[i]*dn;
	if (!FT_IntrpStateVarAtCoords(front,SOLUTE_COMP,p1,solute,
			getStateSolute,&s1))
       	    s1 = s0;
       	for (i = 0; i < dim; ++i)
       	    p2[i] = p1[i] + nor[i]*dn;
	if (!FT_IntrpStateVarAtCoords(front,SOLUTE_COMP,p2,solute,
			getStateSolute,&s2))
       	    s2 = s1;
	state = (negative_component(oldhs) == SOLUTE_COMP) ? 
		(STATE*)left_state(newp) : (STATE*)right_state(newp);
	state->solute = (4.0*s1 - s2)/3.0;
        return;
}       /* neumann_point_propagate */

static  void dirichlet_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        int i, dim = front->rect_grid->dim;
        STATE *sl,*sr,*state,*bstate;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;

        for (i = 0; i < dim; ++i)
            Coords(newp)[i] = Coords(oldp)[i];
	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
       	    FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	    state =  (STATE*)left_state(newp);
	    state->solute = bstate->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	    state =  (STATE*)right_state(newp);
	    state->solute = bstate->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	}
	else
	{
       	    FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
	    state =  (STATE*)left_state(newp);
	    state->solute = sl->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	    state =  (STATE*)right_state(newp);
	    state->solute = sr->solute;
	    if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	    if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	}
        return;
}       /* dirichlet_point_propagate */


static  void reaction_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double nor_speed,vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double nor[MAXD];
        double p1[MAXD],p2[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
        double s0,s1,s2,grad_s,ans;
        STATE *sl,*sr,*state;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double *solute = cRparams->field->solute;
        double D = cRparams->D;
        double k = cRparams->k;
        double C_0 = cRparams->C_0;
        double C_eq = cRparams->C_eq;
        double rho_s = cRparams->rho_s;
	double kappa;
	POINT_PROP_SCHEME point_prop_scheme = cRparams->point_prop_scheme;
	static void (*reaction_scheme)(CRT_PARAMS*,double,double,double,
				double,double,double*);
	static boolean first = YES;
	static double max_nor_speed = 0.0;

	if (first)
	{
	    first = NO;
	    switch (point_prop_scheme)
	    {
	    case EXPLICIT_EULER:
		reaction_scheme = euler_forward_scheme;
	    	break;
	    case IMPLICIT_EULER:
		reaction_scheme = euler_backward_scheme;
	    	break;
	    case MIDDLE_POINT:
		reaction_scheme = middle_point_scheme;
	    	break;
	    case CONSTANT_STATE:
		reaction_scheme = constant_state;
	    	break;
	    default:
		(void) printf("ERROR: Unknow reaction scheme!\n");
		clean_up(ERROR);
	    } 
	}
        FT_NormalAtPoint(oldp,front,nor,SOLUTE_COMP);
	/*Not yet working properly*/
	if (cRparams->add_curvature)
            FT_CurvatureAtPoint(oldp,front,&kappa);
	else
	    kappa = 0.0;

        dn = grid_size_in_direction(nor,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor[i]*dn;

        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == SOLUTE_COMP) ? sl : sr;
        s0 = state->solute;

	if (!FT_IntrpStateVarAtCoords(front,SOLUTE_COMP,p1,solute,
				getStateSolute,&s1))
        {
            s1 = state->solute;
        }

        grad_s = (s1 - s0)/dn;
	if (point_prop_scheme == CONSTANT_STATE)
	    nor_speed = std::max(0.0,D*grad_s/rho_s);
	else
	    nor_speed = std::max(0.0,k*(s0 - C_eq)/rho_s);
        for (i = 0; i < dim; ++i)
        {
            vel[i] = nor[i]*nor_speed;
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
        }
	if (max_nor_speed < nor_speed)
	{
	    max_nor_speed = nor_speed;
	    if (debugging("step_size"))
	    	(void) printf("Scaled max_prop_spacing = %f\n",
				max_nor_speed/dn*dt);
	}
	reaction_scheme(cRparams,dt,dn,kappa,s0,s1,&ans);

	/* Update the state of the new interface point */

	state = (negative_component(oldhs) == CRYSTAL_COMP) ? 
			(STATE*)left_state(newp) : (STATE*)right_state(newp);
        state->solute = rho_s;

	state = (negative_component(oldhs) == SOLUTE_COMP) ? 
			(STATE*)left_state(newp) : (STATE*)right_state(newp);
        s0 = state->solute = ans;

	if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
	nor_speed = 0.0;
        for (i = 0; i < dim; ++i)
        {
            vel[i] = nor[i]*k*(s0 - C_eq)/rho_s;
	    nor_speed += sqr(vel[i]);
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }
	FT_RecordMaxFrontSpeed(dim,sqrt(nor_speed),NULL,Coords(newp),front);
}       /* reaction_point_propagate */

static 	void euler_forward_scheme(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double D = cRparams->D;
        double k = cRparams->k;
        double C_eq = cRparams->C_eq;

	*ans = s0 + 2.0*dt/dn*(D*(s1 - s0)/dn - k*(s0 - C_eq));
}	/* end euler_forward_scheme */

static 	void euler_backward_scheme(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double D = cRparams->D;
        double k = cRparams->k;
        double C_eq = cRparams->C_eq;
	double c1,c2;

	c1 = 0.5*D*dt*(1.0/dn + kappa)/dn;
	c2 = 0.5*k*dt/dn;
	*ans = (s0 + c1*(s1 - s0) - c2*(s0 - C_eq) + c1*s1 + c2*C_eq)
			/(1.0 + c1 + c2);
}	/* end euler_backward_scheme */

static 	void middle_point_scheme(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double D = cRparams->D;
        double k = cRparams->k;
        double C_eq = cRparams->C_eq;
	double c1,c2;

	c1 = D*dt*(2.0/dn + kappa)/dn;
	c2 = 2.0*k*dt/dn;
	*ans = (s0 + c1*s1 + c2*C_eq)/(1.0 + c1 + c2);
}	/* end middle_point_scheme */

static 	void constant_state(
	CRT_PARAMS *cRparams,
	double dt,
	double dn,
	double kappa,
	double s0,
	double s1,
	double *ans)
{
        double C_eq = cRparams->C_eq;

	*ans = C_eq;
}	/* end middle_point_scheme */

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
        else if ((scheme[0] == 'C' || scheme[0] == 'c') &&
            (scheme[1] == 'S' || scheme[1] == 's'))
            cRparams->point_prop_scheme = CONSTANT_STATE;
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
	boolean crystal_exist = NO;
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
	    crystal_exist = YES;
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
	crystal_exist = pp_max_status(crystal_exist);
#endif /* defined (__MPI__) */
	if (!crystal_exist)
	    return NO;

	/* Preparation for counting */

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
	ratio = ((double)N)/((double)Nc);
	*frac_dim = (double)dim + log(ratio)/log(h[0]);
	return pp_max_status(bdry_reached);
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

	s_params->num_floor_seeds = 0;
	s_params->num_ceiling_seeds = 0;
	s_params->num_space_seeds = 0;
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
	    s_params->add_space_seed_pert = NO;		// default
	    CursorAfterString(infile,"Enter yes to add perturbation:");
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
		s_params->add_space_seed_pert = YES;
	    	CursorAfterString(infile,"Enter number of period:");
		fscanf(infile,"%lf ",&s_params->nu);
	    	(void) printf("%f\n",s_params->nu);
	    	CursorAfterString(infile,"Enter amplitude:");
		fscanf(infile,"%lf ",&s_params->amp);
	    	(void) printf("%f\n",s_params->amp);
	    	CursorAfterString(infile,"Enter phase shift:");
		fscanf(infile,"%lf ",&s_params->phase);
	    	(void) printf("%f\n",s_params->phase);
	    }
	}

	CursorAfterString(infile,"Enter radius of seeds:");
	fscanf(infile,"%lf",&s_params->seed_radius);
	(void) printf("%f\n",s_params->seed_radius);
}	/* end read_seed_params */

static void read_crt_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,k,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	STATE state;
	HYPER_SURF *hs;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		if (rect_boundary_type(front->interf,i,0) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
						i,0);
		sprintf(msg,"For lower boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_SetDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	    if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
		if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
						i,1);
		sprintf(msg,"For upper boundary in %d-th dimension",i);
		CursorAfterString(infile,msg);
		(void) printf("\n");
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		(void) printf("%s\n",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_SetDirichletBoundary(front,NULL,NULL,NULL,
					(POINTER)&state,hs);
		    break;
		default: 
		    printf("ERROR: Dirichlet type %s not implemented\n",s);
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_crt_dirichlet_bdry_data */


static void initFrontStates(
	Front *front)
{
	INTERFACE *intfc = front->interf;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *sl,*sr;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double rho_s = cRparams->rho_s;

	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);

            if (positive_component(hs) == SOLUTE_COMP)
                sr->solute = cRparams->C_eq;
            else if (positive_component(hs) == CRYSTAL_COMP)
                sr->solute = rho_s;
            else
                sr->solute = 0.0;
            if (negative_component(hs) == SOLUTE_COMP)
                sl->solute = cRparams->C_eq;
            else if (negative_component(hs) == CRYSTAL_COMP)
                sl->solute = rho_s;
            else
                sl->solute = 0.0;
        }
}	/* end initFrontStates */

void solute_read_front_states(
	FILE *infile,
	Front *front)
{
        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        double x;


        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface solute states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            sl->solute = x;
            fscanf(infile,"%lf",&x);
            sr->solute = x;
        }
}	/* end solute_read_front_states */

void solute_print_front_states(
	FILE *outfile,
	Front *front)
{
	int i,j,k,l,index;
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;

        /* Initialize states at the interface */
        fprintf(outfile,"Interface solute states:\n");
        int count = 0;
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            count++;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",sl->solute,sr->solute);
        }

}	/* end solute_print_front_states */

