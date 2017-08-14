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
*				ifluid.c
* This program is modified from example0.c for solving incompressible flow.
*
* The solver is define in lcartsn.h/c.
*
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include "iFluid.h"
#include "ifluid_basic.h"
#include "iFluid_debug.h"

	/*  Function Declarations */
static void init_io( int,char**);
static void ifluid_driver(Front*,Incompress_Solver_Smooth_Basis*,F_BASIC_DATA*);
static void ifluid_driver_debug(Front*,Incompress_Solver_Smooth_Basis*);
//Accuracy testing ifluid_driver

static int l_cartesian_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);
static int l_cartesian_vel_MAC_vd(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);

/************************************************************
 * 	init the whole domain as a single component.
 * This is used as testing.
 */
void initIntfc_oneComponent(Front *front,LEVEL_FUNC_PACK *level_func_pack);
double level_oneComponent_2D(POINTER func_params,double *coords);
double level_oneComponent_3D(POINTER func_params,double *coords);
double level_oneComponent_3D_Cylindrical(POINTER func_params,double *coords);
/******************************************************/

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean RegridRun;
boolean RegridRestart;
boolean ReadFromInput;
int RestartStep;
const boolean binary = NO;
const boolean isTesting = NO; //The switch for the testing mode
const boolean isVd = YES; //The switch for variable density simulation

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

int main(int argc, char **argv)
{
	static Front front;
	// static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	IF_PARAMS iFparams;
	IF_PROB_TYPE prob_type;

	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);//Read parameters from command line
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc before the FrontStartUp
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	//printf("Passed PetscInitialize()\n");

	/*Construct Incompress Solver l_cartesian*/

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;


	Incompress_Solver_Smooth_Basis *l_cartesian;//Construct Incompress Solver

	if(f_basic.dim == 2)
	{
	    if (isTesting)
            {
                printf("\nUsing the debugging class!\n");
		l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian_Debug(front);
            }
	    else
		l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	}
	else if(f_basic.dim == 3)
	{
	    if (f_basic.coord_system == CYLINDRICAL_REMAP)
	    /*
            //Cylindrical coordinate in 3D
	    {
		if (isTesting)//Testing mode using exact solution
		{
		    printf("\nUsing the debugging class!\n");
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cylindrical_Debug(front);
		}
		else
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cylindrical(front);
	    }
            */
            {
            }
	    else if(f_basic.coord_system == SPHERICAL_REMAP)
	    //Sperical coordinate in 3D, not implemented yet
	    {
	    }
	    else
	    //Default: Rectangular coordinate in 3D
	    {
		if (isTesting)
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian_Debug(front);
		else
		    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);
	    }
	}

        sprintf(restart_state_name,"%s/restart/state.t%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s/restart/intfc-t%s",restart_name,
			right_flush(RestartStep,7));

	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-p%s",restart_name,
			right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-p%s",restart_state_name,
                        right_flush(pp_mynode(),4));
	}

	FT_ReadSpaceDomain(in_name,&f_basic);
        printf("Passed FT_ReadSpaceDomain()\n");
        //if RestartRun, read from restart/intfc files during FT_StartUp()
        printf("Entered FT_StartUp()\n");
	FT_StartUp(&front,&f_basic);
        printf("Passed FT_StartUp()\n");

/*
            {
                printf("\n\nStorage after the first FT_StartUp()\n");
                int totalNum = -1;
                const char *blockName;
                blockName = "table";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                blockName = "chunk";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                long_alloc_view(stdout);
            }
            fflush(stdout);
*/

        FT_InitDebug(in_name);

        RegridRun = f_basic.RegridRun;
        RegridRestart = f_basic.RegridRestart;

	iFparams.dim = f_basic.dim;
	front.extra1 = (POINTER)&iFparams;

	if (!isTesting)
	    read_iF_prob_type(in_name,&prob_type);

	read_iFparams(in_name,&iFparams);
	read_iF_movie_options(in_name,&iFparams);
	if (debugging("trace")) printf("Passed read_iFparams()\n");

	/* Initialize interface through level function */
	if (isTesting)
	    initIntfc_oneComponent(&front, &level_func_pack);
	else
	    setInitialIntfc(&front,&level_func_pack,in_name,prob_type);
	if (debugging("trace")) printf("Passed setInitialIntfc()\n");

	if (!RestartRun)
	{
	    if (f_basic.dim == 3)
                level_func_pack.set_3d_bdry = YES;
	    FT_InitIntfc(&front,&level_func_pack);
	    if (debugging("trace"))
	    {
		char test_name[100];

		printf("Passed FT_InitIntfc()\n");
        /*
		switch (f_basic.dim)
		{
		case 2:
		    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
		    //xgraph_2d_intfc(test_name,front.interf);
		    break;
		case 3:
		    sprintf(test_name,"init_intfc-%d.xg",pp_mynode());
		    //gview_plot_interface("gv-init",front.interf);
		    break;
		}
        */
	    }
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
            if (f_basic.dim < 3)
                FT_ClipIntfcToSubdomain(&front);
            if (debugging("trace"))
                printf("Passed read_iF_dirichlet_bdry_data()\n");
	}

	/* Initialize velocity field function */
	velo_func_pack.func_params = (POINTER)l_cartesian;
        if (!isVd)
        {
            velo_func_pack.func = l_cartesian_vel;
            velo_func_pack.point_propagate = ifluid_point_propagate;
        }
        else
        {
            if (debugging("trace"))
                printf("velo_func_pack.func = l_cartesian_vel_MAC_vd\n");
            velo_func_pack.func = l_cartesian_vel_MAC_vd;
            velo_func_pack.point_propagate = ifluid_point_propagate_vd;
        }
	FT_InitVeloFunc(&front,&velo_func_pack);
	if (debugging("trace"))
	    printf("Passed FT_InitVeloFunc()\n");

        if (!isVd)
	    l_cartesian->initMesh();
        else
            l_cartesian->initMesh_vd();
        if (debugging("trace"))
            printf("Passed l_cartesian.initMesh()\n");

            if (debugging("storage"))
            {
                printf("\n\nStorage after initMesh_vd()\n");
                int totalNum = -1;
                const char *blockName;
                blockName = "table";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                blockName = "chunk";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                long_alloc_view(stdout);
            }


	if (debugging("sample_velocity"))
	    l_cartesian->initSampleVelocity(in_name);

	init_fluid_state_func(l_cartesian,prob_type); //initialize "getInitialState"

        //read from restart/state files
	if (RestartRun)
        {
            if (debugging("storage"))
            {
                printf("\n\nStorage before the first readFrontInteriorStates_vd()\n");
                int totalNum = -1;
                const char *blockName;
                blockName = "table";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                blockName = "chunk";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                long_alloc_view(stdout);
            }

            if (!isVd)
	        l_cartesian->readFrontInteriorStates(restart_state_name,binary,RegridRestart);
            else {
                l_cartesian->readFrontInteriorStates_vd(restart_state_name,binary,RegridRestart);
                printf("Passed readFrontInteriorStates_vd()\n");
            }

            if (debugging("storage"))
            {
                printf("\n\nStorage after the first readFrontInteriorStates_vd()\n");
                int totalNum = -1;
                const char *blockName;
                blockName = "table";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                blockName = "chunk";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                long_alloc_view(stdout);
            }



        }
	else
        {
            if (!isVd)
            {
	        l_cartesian->setInitialCondition();
                if (debugging("trace"))
                    printf("Passed setInitialCondition()\n");
            }
            else
            {
                if (level_func_pack.is_RS_RV) //for RS_RV case
                {
                    l_cartesian->setInitialCondition_RSRV_vd(&level_func_pack);
                    if (debugging("trace"))
                        printf("Passed setInitialCondition_RSRV_vd()\n");
                }
                else if (level_func_pack.is_RS_SY) //for Smeeton-Youngs experiment
                {
                    l_cartesian->setInitialCondition_RSSY_vd(&level_func_pack);
                    if (debugging("trace"))
                        printf("Passed setInitialCondition_RSSY_vd()\n");
                }
                else
                {
                    l_cartesian->setInitialCondition_vd(&level_func_pack);
                    if (debugging("trace"))
                        printf("Passed setInitialCondition_vd()\n");
                }
            }
        }

        /* Enter the iFluid Driver */
	if (isTesting)
	{
	    printf("\nEntering ifluid driver for testing problem!\n");
	    ifluid_driver_debug(&front, l_cartesian);
	}
	else
	    ifluid_driver(&front, l_cartesian, &f_basic);

	PetscFinalize();
	clean_up(0);
} /* end main() */


static  void ifluid_driver(
        Front *front,
	Incompress_Solver_Smooth_Basis *l_cartesian,
        F_BASIC_DATA *f_basic)
{
        printf("\nnode = %d, PID = %d\n",pp_mynode(), getpid());
        fflush(stdout);

        double CFL;
	int dim = front->rect_grid->dim;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

        /**** The zero time step ****/

        if (!RestartRun || RegridRestart)
	{
            if (debugging("trace"))
                (void) printf("\nEnter FT_RedistMesh(front) in Zeroth step\n");
	    FT_RedistMesh(front);
            if (debugging("trace"))
                (void) printf("Leave FT_RedistMesh(front) in Zeroth step\n");
	}

        if (!RestartRun)
        {
	    FT_ResetTime(front);
            (void)printf("\n");
	    if (debugging("trace"))
                printf("Zeroth step: Before FT_Propagate() front->dt = %f,"
                        "l_cartesian->max_dt = %f\n",front->dt, l_cartesian->max_dt);
        // TODO && FIXME: Install Reflection Boundary Condition Constraint
        //((Incompress_Solver_Smooth_3D_Cartesian*)l_cartesian)->enforceReflectionState();
         FT_Propagate(front);

	    if (debugging("trace"))
                (void) printf("Zeroth step: Calling ifluid solve()\n");
            if (!isVd)
                l_cartesian->solve(front->dt);
            else
                l_cartesian->solve_vd(front->dt);
	    if (debugging("trace")) printf("Zeroth step: Passed ifluid solve()\n");

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
            FT_SetOutputCounter(front);
        }
        else
        {
	    FT_SetOutputCounter(front);
        }

	FT_TimeControlFilter(front); // Adjust front->dt when front->time is very close to print/movie interval


        /**** The zero time step output  ****/

        if (!RestartRun)
        {
            //print P-node/state files for making plot
            l_cartesian->printInteriorVelocity(out_name);
            l_cartesian->printExpandedMesh(out_name, binary);

            if (debugging("trace"))
                (void) printf("Calling initMovieVariables()\n");
            l_cartesian->initMovieVariables();
            if (debugging("trace"))
                (void) printf("Calling FT_AddMovieFrame()\n");
            //print P-node/intfc files for making plot
            FT_AddMovieFrame(front,out_name,binary);
        }
        //print restart files for RestartStep, which is necessary for supporting automatic restart
/*
        else
        {
            if (!RegridRun && !RegridRestart && !(front->repartition)) {
                //print restart/intfc files
                FT_Save(front,out_name);
                //print restart/state files
                if (!isVd)
                    l_cartesian->printFrontInteriorStates(out_name, binary);
                else
                    l_cartesian->printFrontInteriorStates_vd(out_name, binary);
	    }
        }
*/
        /******************************************************************/

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("step_size"))
                printf("Time step from start: %f\n",front->dt);

        front->restart_step = RestartStep;
        front->last_restart_step = RestartStep;
        front->start_small_dt = NO;
        front->end_small_dt = NO;
        front->bad_step = NO;
        double sum_num_reduction,sum_num_points;

	// spatial-average and temporal-average velocity variance
	const int nStep_velo_var = 1; // CFL = 0.45 => nStep_velo_var = 2
        // only calculate velocity variance, and then exit
        const boolean onlyCalcVeloVar = NO;

        for (;;)
        {
            /* compute and print terminal velocity, top and bottom of interior surface, alpha, and theta for RT simulation*/

            if (isVd && dim==3 && (front->step - front->restart_step - 1)%nStep_velo_var==0)
	    {
                l_cartesian->computeRTParameters(front->dt,out_name, nStep_velo_var);

		if (onlyCalcVeloVar)
                    clean_up(0);
	    }

            /* Propagating interface for time step dt */

	    if (debugging("trace"))
                printf("\nBefore FT_Propagate()\n");

            front->num_reduction = 0;
            front->num_points = 0;

            start_clock("FT_Propagate");
            // TODO && FIXME: Reflection Boundary Condition
            //((Incompress_Solver_Smooth_3D_Cartesian*)l_cartesian)->enforceReflectionState();
            FT_Propagate(front);
            stop_clock("FT_Propagate");

            sum_num_reduction = front->num_reduction;
            if (pp_numnodes() > 1)
                pp_global_sum(&sum_num_reduction,1);
            front->num_reduction = sum_num_reduction;
            sum_num_points = front->num_points;
            if (pp_numnodes() > 1)
                pp_global_sum(&sum_num_points,1);
            front->num_points = sum_num_points;


            if (debugging("storage"))
            {
                printf("\n\nStorage after FT_Propagate() in time step %d\n", front->step);
                int totalNum = -1;
                const char *blockName;
                blockName = "table";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                blockName = "chunk";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                long_alloc_view(stdout);
            }



            if (front->restart_small_dt == YES)
            {
                fprintf(stdout, "use of small CFL number from the time step: %d\n",front->last_restart_step);
                front->start_small_dt = YES;
                front->end_small_dt = NO;

                restart_state_name      = f_basic->restart_state_name;
                restart_name            = f_basic->restart_name;
                sprintf(restart_state_name,"%s/restart/state.t%s",f_basic->out_name,
                                right_flush(front->last_restart_step,7));
                sprintf(restart_name,"%s/restart/intfc-t%s",f_basic->out_name,
                                right_flush(front->last_restart_step,7));
                if (pp_numnodes() > 1)
                {
                    sprintf(restart_name,"%s-p%s",restart_name,
                                right_flush(pp_mynode(),4));
                    sprintf(restart_state_name,"%s-p%s",restart_state_name,
                                right_flush(pp_mynode(),4));
                }

                printf("f_basic->restart_state_name = %s\n", f_basic->restart_state_name);
                printf("f_basic->restart_name = %s\n", f_basic->restart_name);

                FT_ReadSpaceDomain(in_name,f_basic);

                //if (debugging("storage"))
                {
                    printf("\n\nStorage before clean front in automatic rs %d\n", front->step);
                    int totalNum = -1;
                    const char *blockName;
                    blockName = "table";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    blockName = "chunk";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    long_alloc_view(stdout);
                }

                //delete intfc and grid_intfc
                if (front != NULL)
                {
                    if (prev_interface(front->interf) != NULL)
                        prev_interface(front->interf) = NULL;
                    if (front->interf != NULL) {
                        (void) delete_interface(front->interf);
                        front->interf = NULL;
                    }
                    if (front->grid_intfc != NULL)
                        FT_FreeGridIntfc(front);
                }

		//delete front
		//TODO: might have bug if we call free_front(front)???
//                if (front)	free_front(front);

                //if (debugging("storage"))
                {
                    printf("\n\nStorage after clean front in automatic rs %d\n", front->step);
                    int totalNum = -1;
                    const char *blockName;
                    blockName = "table";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    blockName = "chunk";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    long_alloc_view(stdout);
                }

                FT_StartUp(front, f_basic);

                fflush(stdout);
                //if (debugging("storage"))
                {
                    printf("\n\nStorage after FT_StartUp() in automatic rs %d\n", front->step);
                    int totalNum = -1;
                    const char *blockName;
                    blockName = "table";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    blockName = "chunk";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    long_alloc_view(stdout);
                }


                FT_InitDebug(in_name);
                if (!isVd)
                    l_cartesian->readFrontInteriorStates(restart_state_name, binary, RegridRestart);
                else
                    l_cartesian->readFrontInteriorStates_vd(restart_state_name, binary, RegridRestart);



                if (debugging("storage"))
                {
                    printf("\n\nStorage after readFrontInteriorStates_vd() in automatic rs %d\n", front->step);
                    int totalNum = -1;
                    const char *blockName;
                    blockName = "table";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    blockName = "chunk";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    long_alloc_view(stdout);
                }



                front->restart_small_dt = NO;
                start_clock("FT_Propagate");
                // TODO && FIXME: Reflection Boundary Condition
                //((Incompress_Solver_Smooth_3D_Cartesian*)l_cartesian)->enforceReflectionState();
                FT_Propagate(front);
                stop_clock("FT_Propagate");


                if (debugging("storage"))
                {
                    printf("\n\nStorage after FT_Propagate() using small_dt in automatic rs %d\n", front->step);
                    int totalNum = -1;
                    const char *blockName;
                    blockName = "table";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    blockName = "chunk";
                    totalNum = get_number_of_blockName(stdout, blockName);
                    printf("number_of_%s = %d\n", blockName, totalNum);
                    long_alloc_view(stdout);
                }



            } // if (front->restart_small_dt == YES)

	    if (debugging("trace")) printf("Passed FT_Propagate()\n");

            if (debugging("trace")) printf("Calling ifluid solve()\n");

            if (!isVd)
	        l_cartesian->solve(front->dt);
            else
                l_cartesian->solve_vd(front->dt);

	    if (debugging("trace"))
                printf("Passed ifluid solve_vd()\n");
	    if (debugging("storage"))
            {
                printf("\n\nStorage after solve_vd() at end of time step\n");
                int totalNum = -1;
                const char *blockName;
                blockName = "table";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                blockName = "chunk";
                totalNum = get_number_of_blockName(stdout, blockName);
                printf("number_of_%s = %d\n", blockName, totalNum);
                long_alloc_view(stdout);
            }

	    FT_AddTimeStepToCounter(front);

            if(front->step == 0)
                front->time = 0;

            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
	    //if (debugging("step_size"))
            //    (void) printf("Time step from FT_SetTimeStep(): %20.14f\n",
	    //			front->dt);
            front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
            if(front->start_small_dt == YES && front->end_small_dt == NO)
                front->dt = front->dt/9.0;

	    if (debugging("step_size"))
                (void) printf("Time step from l_cartesian->max_dt(): %20.14f\n",
					front->dt);


            /* Output section */

            (void) printf("\ntime = %20.14f   step = %5d   next dt = %20.14f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (RegridRun || front->repartition == YES)
            {
                //print restart/intfc files
                FT_Save(front,out_name);
                //print restart/state files
                if (!isVd)
                    l_cartesian->printFrontInteriorStatesRegridRep(out_name, binary, RegridRun);
                else
                    l_cartesian->printFrontInteriorStatesRegridRep_vd(out_name, binary, RegridRun);
            }
            front->regrid_restart = NO;
            if (RegridRun || front->repartition == YES)
                clean_up(0);


            printf("RestartStep = %d\n", RestartStep);
            printf("last_restart_step = %d\n", front->last_restart_step);

            if (FT_IsSaveTime(front)) //|| (front->step)%40 == 0)
		 // || front->step == 892 //|| front->step == 981)
                 // || front->step == RestartStep + 9 || front->step == RestartStep + 10)
	    {
                front->end_small_dt = YES;
                front->last_restart_step = front->step;
                fprintf(stdout, "last_restart_step %d\n",front->last_restart_step);
                //print restart/intfc files
            	FT_Save(front,out_name);
                //print restart/state files
                if (!isVd)
		    l_cartesian->printFrontInteriorStates(out_name, binary);
                else
                    l_cartesian->printFrontInteriorStates_vd(out_name, binary);
	    }
	    if (debugging("trace"))
                (void) printf("After print output()\n");
            if (FT_IsMovieFrameTime(front))
	    {
                //print P-node/state files for making plot
//                l_cartesian->printInteriorVelocity(out_name);

	    	if (debugging("trace"))
		    (void) printf("Calling initMovieVariables()\n");
	        l_cartesian->initMovieVariables();
	    	if (debugging("trace"))
		    (void) printf("Calling FT_AddMovieFrame()\n");
                //print P-node/intfc files for making plot
            	FT_AddMovieFrame(front,out_name,binary);
	    }

            if (FT_TimeLimitReached(front))
                break;

	    FT_TimeControlFilter(front);
	    if (debugging("step_size"))
                (void) printf("Time step from FT_TimeControlFilter(): %f\n",
                                        front->dt);
            fflush(stdout);
        }
	if (debugging("trace")) printf("After time loop\n");
} /* end ifluid_driver */


// The testing mode using exact solution
static  void ifluid_driver_debug(
        Front *front,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        double CFL;
	int dim = front->rect_grid->dim;

	Curve_redistribution_function(front) = full_redistribute;

	FT_ReadTimeControl(in_name,front);
	CFL = Time_step_factor(front);

	if (RestartRun)
	{
	    FT_ParallelExchIntfcBuffer(front);
	}
	else
	{
	    FT_RedistMesh(front);
	}

        if (!RestartRun)
        {
	    FT_ResetTime(front);
            FT_SetOutputCounter(front);
	    if (debugging("trace"))
		printf("Before FrontProp() front->dt = %f\n",front->dt);
            FT_Propagate(front);

	    if (debugging("trace")) printf("Calling ifluid solve()\n");
            if (!isVd)
                l_cartesian->solve(front->dt);
            else
                l_cartesian->solve_vd(front->dt);
	    if (debugging("trace")) printf("Passed ifluid solve()\n");

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
        }
        else
        {
	    FT_SetOutputCounter(front);
        }

	FT_TimeControlFilter(front);

	if (debugging("trace"))
	{
	    printf("CFL = %f\n",CFL);
	    printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));
	}

	if (debugging("step_size"))
                printf("Time step from start: %f\n",front->dt);

	double dh;
	dh = std::min(front->rect_grid->h[0], front->rect_grid->h[1]);
	if (dim == 3)
	    dh = std::min(dh, front->rect_grid->h[2]);

	double dt = 0.1 * CFL * dh; // time step for 3D manufacturing solution only
	int nStep = int (front->max_time/dt) + 1;
	dt = front->max_time/ nStep;
	front->dt = dt;

	if(dim==2)
	{
            if (!isVd)
            {
	        ((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,true);
	        //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
	        //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
            }
            else
            {
                ((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot_vd(out_name,front->time,true);
                //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot_vd(out_name,front->time,false);
                //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
            }
	}
	else if(dim==3)
	{
	    if (front->coordinate == 'R' || front->coordinate == 'r')
	    {
                if (!isVd)
                {
              	    ((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,true);
	    	    //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
	    	    //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
                }
                else
                {
                    ((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot_vd(out_name,front->time,true);
                    //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot_vd(out_name,front->time,false);
                    //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
                }
	    }
	    else if (front->coordinate == 'C' || front->coordinate == 'c')
	    {
	    	((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,false);
	    	//((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,true);
	    	//((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,true);
	    }
	}

	for (int i=1; i<=nStep; i++)
	{
	    front->step = i;
	    front->dt = dt;

            if (!isVd)
	        l_cartesian->solve(front->dt);
            else
                l_cartesian->solve_vd(front->dt);

	    front->time += front->dt;
	    printf("\ntime = %f   step = %5d   next dt = %f\n",
		    front->time,front->step,front->dt);
	    if(dim==2)
	    {
                if (!isVd)
                {
		    ((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,true);
		    //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
		    //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
                }
                else
                {
                    ((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot_vd(out_name,front->time,true);
                    //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot_vd(out_name,front->time,false);
                    //((Incompress_Solver_Smooth_2D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
                }
	    }
	    else if(dim==3)
	    {
		if (front->coordinate == 'R' || front->coordinate == 'r')
		{
                    if (!isVd)
                    {
                        ((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,true);
                        //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,false);
                        //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
                    }
                    else
                    {
                        ((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveStates_Tecplot_vd(out_name,front->time,true);
                        //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveParameters_Tecplot_vd(out_name,front->time,false);
                        //((Incompress_Solver_Smooth_3D_Cartesian_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,false);
                    }
		}
		else
		{
	    	    ((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveStates_Tecplot(out_name,front->time,true);
                    //((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveParameters_Tecplot(out_name,front->time,true);
                    //((Incompress_Solver_Smooth_3D_Cylindrical_Debug*)l_cartesian)->saveDivUPhi_Tecplot(out_name,front->time,true);
		}
	    }
	}

	if (debugging("trace")) printf("After time loop\n");
} /* end ifluid_driver_debug */

static int l_cartesian_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	double *coords = Coords(p);
	((Incompress_Solver_Smooth_Basis*)params)->getVelocity(coords, vel);
}	/* end l_cartesian_vel */

//for MAC grid
static int l_cartesian_vel_MAC_vd(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        //printf("Enter l_cartesian_vel_MAC_vd()\n");

        double *coords = Coords(p);
        ((Incompress_Solver_Smooth_Basis*)params)->getVelocity_MAC_vd(coords, vel);
}       /* end l_cartesian_vel_MAC_vd */

void initIntfc_oneComponent(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack)
{
    IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
    iFparams->m_comp1 = LIQUID_COMP1;
    iFparams->m_comp2 = LIQUID_COMP2; //added by wenlin hu

    level_func_pack->neg_component = LIQUID_COMP1;
    level_func_pack->pos_component = LIQUID_COMP2;
    level_func_pack->is_RS_RV = NO;
    level_func_pack->is_RS_SY = NO;

    if(front->rect_grid->dim==2)
	level_func_pack->func = level_oneComponent_2D;
    else if(front->rect_grid->dim==3)
    {
	if (front->coordinate == 'R' || front->coordinate == 'r')
	    level_func_pack->func = level_oneComponent_3D;
	else if(front->coordinate == 'C' || front->coordinate == 'c')
	{
	    printf("\nSetting the level_func_pack to cylindrical case!!\n");
	    level_func_pack->func = level_oneComponent_3D_Cylindrical;
	}
    }
    level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
}	/* end initCirclePlaneIntfc */

double level_oneComponent_2D(
        POINTER func_params,
        double *coords)
{
    double dist = sqrt(
	    sqr(coords[0]-0.5) +
	    sqr(coords[1]-0.5)) - 0.2;
    return dist;
}       /* end level_circle_func */

double level_oneComponent_3D(
        POINTER func_params,
        double *coords)
{
/*
    double dist = sqrt(
	    sqr(coords[0]-0.5) +
	    sqr(coords[1]-0.5) +
	    sqr(coords[2]-0.5)) - 0.2;
*/

    double dist = coords[2] - 0.1875; //mid-plane
    return dist;
}       /* end level_circle_func */

double level_oneComponent_3D_Cylindrical(
        POINTER func_params,
        double *coords)
{
    double dist = sqrt(
	    sqr(coords[0]-3) +
	    sqr(coords[1]-0) +
	    sqr(coords[2]-1.5)  ) - 0.3;
    return dist;
//    return 1;
}       /* end level_circle_func */
