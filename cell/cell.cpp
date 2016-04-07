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
*				cell.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static double sphere_func(POINTER,double*);
static int cell_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
			HYPER_SURF*,double*);

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
	int dim;
        double center[3];
	double radius;
} TEST_SPHERE_PARAMS;


typedef struct {
	int dim;
        double coeff;
	double epsilon;
        unsigned short int rand_seeds[3];
} RNORV_PARAMS;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	TEST_SPHERE_PARAMS s_params;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	RNORV_PARAMS norv_params; /* velocity function parameters */
	Locstate  sl;

	f_basic.dim = 2;	
	FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	//f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; 	//f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 64;	f_basic.gmax[1] = 64; 	//f_basic.gmax[2] = 64;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
	//f_basic.boundary[2][0] = f_basic.boundary[2][1] = DIRICHLET_BOUNDARY;
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
	    s_params.dim = 2;
	    s_params.center[0] = 0.25;
	    s_params.center[1] = 0.25;
	    //s_params.center[2] = 0.25;
	    s_params.radius = 0.15;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&s_params;
	    level_func_pack.func = sphere_func;

	    FT_InitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	norv_params.dim = 2;
	norv_params.coeff = 0.1;
	norv_params.epsilon = 0.0001;
	norv_params.rand_seeds[0] = 5635;
	norv_params.rand_seeds[1] = 5523;
	norv_params.rand_seeds[2] = 1552;

	velo_func_pack.func_params = (POINTER)&norv_params;
	velo_func_pack.func = cell_vel;

	FT_InitVeloFunc(&front,&velo_func_pack);

	/* For geometry-dependent velocity, use first 
	* order point propagation function, higher order
	* propagation requires surface propagate, currently
	* in writing, not yet in use. The following override
	* the assigned fourth_order_point_propagate.
	*/

	front._point_propagate = first_order_point_propagate;

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        int ip,im,status,count;
        Front *newfront;
        double dt,dt_frac,CFL;
        boolean is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        double fcrds[MAXD];
        int  dim = front->rect_grid->dim;

	front->max_time = 1000; 
	front->max_step = 1000000;
	front->print_time_interval = 100;
	front->movie_frame_interval = 10;

        CFL = Time_step_factor(front) = 0.5;

        printf("dim = %d\n", dim);
	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
            redistribute(front,YES,NO);

            front->time = 0.0;
            front->dt = 0.0;
	    front->step = 0;

	    // Always output the initial interface.
	    FT_Save(front,out_name);
            FT_AddMovieFrame(front,out_name,binary);
            ip = im = 1;

	    // This is a virtual propagation to get maximum front 
	    // speed to determine the first time step.

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            front->dt = CFL*FrontHypTimeStep(front); 
	}
	else
	{
	    ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
	}

	front->dt = FrontOutputTimeControl(front,
			&is_movie_time,&is_print_time,
			&time_limit_reached,&im,&ip);

        for (;;)
        {
	    /* Propagating interface for time step dt */

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            assign_interface_and_free_front(front,newfront);

            ++front->step;
            front->time += front->dt;

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            front->dt = CFL*FrontHypTimeStep(front); 

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
                print_front_output(front,out_name);
            if (is_movie_time || time_limit_reached)
                FT_AddMovieFrame(front,out_name,binary);

            if (time_limit_reached)
                    break;

	    /* Output section, next dt may be modified */

	    front->dt = FrontOutputTimeControl(front,
			&is_movie_time,
			&is_print_time,
			&time_limit_reached,
			&im,&ip);
        }
        (void) delete_interface(front->interf);
}       /* end test_propagate */

/********************************************************************
 *	Sample (dummbell 3D) level function for the initial interface    *
 ********************************************************************/

static double sphere_func(
        POINTER func_params,
        double *coords)
{
        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
	double *center,R;
	double distance;
	int i,dim;

	dim = s_params->dim;
	center = s_params->center;
	R = s_params->radius;

	distance = 0.0;
	for (i = 0; i < dim; ++i)
	    distance += sqr(coords[i] - center[i]);
	distance = sqrt(distance) - R;

        return distance;

}       /* end sphere_func */

static int cell_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	RNORV_PARAMS *norv_params;
        RECT_GRID *gr = front->rect_grid;
	BOND  *b = Bond_of_hse(hse);
	CURVE *c = Curve_of_hs(hs);
	double *L = gr->GL;
	double *U = gr->GU;
	double *h = gr->h;
	int i,j;
	static int nnn=0;
	double coeff;
	double curvature,pt[3];
	double nor[MAXD];
	double rn,seta;
	static boolean is_wall = NO,is_upper_wall =NO;
	static boolean is_left_wall = NO,is_buttom_wall = NO;
	static int wall = 0, b_wall = 0;
	static int l_wall = 0;
	unsigned short int *rand_seeds;
	int dim;

	norv_params = (RNORV_PARAMS*)params;
	coeff = norv_params->coeff;
	rand_seeds = norv_params->rand_seeds;
	dim = norv_params->dim;
	
	normal(p,hse,hs,nor,front);
	/* xcjia April 2005*/
	/* velocity components from inside stochastic forces */
	for (i = 0; i < norv_params->dim; ++i)
	{
	    double theta = 0.15;
	    coeff *= 0.05;
	    rn = 2.0*erand48(rand_seeds) - 1.0;
	    p->vel[0] = rn*coeff*cos(theta);
	    p->vel[1] = rn*coeff*sin(theta);
	}
	
	/* velocity components from two closest points */
	/* the: 'p' is the starting point of the bond */

	if( p == b->start)
	{
	    if ((b->prev != NULL)&&(b->prev->prev != NULL))
	    {
	        for (i = 0; i < norv_params->dim; ++i)
	        {
		    p->vel[i] = 0.5*p->vel[i] + 0.25*b->prev->start->vel[i]
		  	      + 0.25*b->prev->prev->start->vel[i]; 
	            if (debugging("cell_k"))
	            {
		        printf("b->prex->vel[%d] = %f\t",i,
					b->prev->start->vel[i]);
		        printf("p->vel[%d] = %f\t",i,p->vel[i]);
	 	        printf("b->next->vel[%d] = %f\n",i,
					b->prev->prev->start->vel[i]);
	            }
		}
	    }
	}
	else
	{
	    if ((b->prev != NULL)&&(b->next != NULL))
	    {
	        for (i = 0; i < norv_params->dim; ++i)
	        {
		    p->vel[i] = 0.4*p->vel[i] + 0.2*b->start->vel[i]
		 	      + 0.2*b->prev->start->vel[i];
	        }
	    }
	}
	/* velocity components from outside forces */
	p->vel[0] += 0.0002;
	p->vel[1] += 0.0002;
	for (i = 0; i < dim; ++i)
	    vel[i] = p->vel[i];
	printf("vel = %24.18g  %24.18g\n",vel[0],vel[1]);
}	/* end normal_vel */
