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
*				finance.c:
*
*	Application of FronTier to Black-Scholes Equation.
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include "finance.h"

	/*  Function Declarations */
static void finance_driver(Front*,CARTESIAN&);
static void excercise_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void read_finance_params(char*,PARAMS*);
static double exact_put_solution(Front*);
static void setInitialIntfc(Front*,LEVEL_FUNC_PACK*);
static void timeDependBoundaryState(double*,HYPER_SURF*,Front*,POINTER,POINTER);
static void get_time_dependent_params(int,FILE*,POINTER*);

struct _FUNC_PARAMS {
	double S1,S2,S3;
	double V1,V2,V3;
	double r,D,E;
	double a,b,c;
	double sigma,t;
	F_TYPE f_type;
};
typedef struct _FUNC_PARAMS FUNC_PARAMS;

static boolean extract_func(double,double*,POINTER);
static boolean extract_func_root(double*,POINTER);
static boolean erf_func(double,double*,POINTER);
static void record_excer_point(Front*,char*);
static void read_dirichlet_bdry_data(char*,Front*,F_BASIC_DATA);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
boolean ReadFromInput;
int RestartStep;
boolean binary = NO;


/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

	/* It is very strange that after using pbs instead of qsub -I 
	 * to login a node, I can make and after I make the file using 
	 * starzero, test_incompressible doesn't converge. 
	 */

int main(int argc, char **argv)
{
	static Front front;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	static PARAMS eqn_params;

	CARTESIAN	cartesian(front);

	/* Initialize basic computational data */

	FT_Init(argc,argv,&f_basic);

	f_basic.size_of_intfc_state = sizeof(STATE);

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_state_name,"%s/state.ts%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s/intfc-ts%s",restart_name,
			right_flush(RestartStep,7));
	if (pp_numnodes() > 1)
	{
            sprintf(restart_name,"%s-nd%s",restart_name,
			right_flush(pp_mynode(),4));
            sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                        right_flush(pp_mynode(),4));
	}

	FT_ReadSpaceDomain(in_name,&f_basic);
	FT_StartUp(&front,&f_basic);
	FT_InitDebug(in_name);

        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	if (debugging("trace")) printf("Passed PetscInitialize()\n");

	/* Initialize interface through level function */

	read_finance_params(in_name,&eqn_params);
	front.extra1 = (POINTER)&eqn_params;

    	level_func_pack.neg_component = EXCERCISE_COMP;
    	level_func_pack.pos_component = BLACK_SCHOLES_COMP; // default comp

	if (!RestartRun)
	{
	    setInitialIntfc(&front,&level_func_pack);
	    FT_InitIntfc(&front,&level_func_pack);
	    read_dirichlet_bdry_data(in_name,&front,f_basic);
	}
	if (debugging("trace")) printf("Passed FT_InitIntfc()\n");

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)&cartesian;
	velo_func_pack.point_propagate = excercise_point_propagate;

	eqn_params.option_price = NULL;
	cartesian.initMesh();

	FT_InitVeloFunc(&front,&velo_func_pack);

	if (debugging("trace")) printf("Passed FT_InitVeloFunc()\n");

	/* Propagate the front */

	finance_driver(&front, cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void finance_driver(
        Front *front,
	CARTESIAN &cartesian)
{
        int ip,im;
        Front *newfront;
        double dt_frac,CFL,tmp;
	boolean is_print_time,is_movie_time,time_limit_reached;
	int i,dim = front->rect_grid->dim;

	Curve_redistribution_function(front) = expansion_redistribute;

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

	if (RestartRun)
	    cartesian.readOptionPrice(restart_state_name);
	else
	    cartesian.setInitialCondition();

        if (!RestartRun)
        {
            // Front standard output
	    FT_ResetTime(front);
            FT_SetOutputCounter(front);
            FT_Save(front,out_name);
            FT_AddMovieFrame(front,out_name,binary);
	    cartesian.printOptionPrice(out_name);
	    if (dim == 1)
	    	cartesian.oneDimPlot(out_name);

	    if (debugging("trace"))
                printf("Before FT_Propagate() front->dt = %f\n",front->dt);
	    FT_Propagate(front);
	    if (debugging("trace")) printf("Calling finance solve()\n");
	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*cartesian.m_dt);
        }
        else
	    FT_SetOutputCounter(front);

	FT_TimeControlFilter(front);
        printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
        fflush(stdout);

	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));

        for (;;)
        {
            /* Propagating interface for time step dt */

	    record_excer_point(front,out_name);
	    cartesian.solve(front->dt);
	    FT_Propagate(front);

	    FT_AddTimeStepToCounter(front);
				
            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*cartesian.m_dt);
	
            /* Output section */

            if (FT_IsSaveTime(front))
	    {
            	FT_Save(front,out_name);
		cartesian.printOptionPrice(out_name);
		if (dim == 1)
	    	    cartesian.oneDimPlot(out_name);
	    }
            if (FT_IsMovieFrameTime(front))
            	FT_AddMovieFrame(front,out_name,binary);

            if (FT_TimeLimitReached(front))
                    break;

	    FT_TimeControlFilter(front);

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

        }
}       /* end finance_driver */

static	void	read_finance_params(
	char *in_name,
	PARAMS *eqn_params)
{
	FILE *infile;
	char scheme[200],opt_type[100];

	infile = fopen(in_name,"r");

	CursorAfterString(infile,"Enter option type:");
	fscanf(infile,"%s",opt_type);
	(void) printf("%s\n",opt_type);
	if (opt_type[0] == 'P' || opt_type[0] == 'p')
	{
	    if (opt_type[4] == 'E' || opt_type[4] == 'E')
	    	eqn_params->f_type = EURO_PUT_OPTION;
	    if (opt_type[4] == 'A' || opt_type[4] == 'a')
	    	eqn_params->f_type = AMRI_PUT_OPTION;
	}
	else if (opt_type[0] == 'C' || opt_type[0] == 'c')
	{
	    if (opt_type[5] == 'E' || opt_type[5] == 'E')
	    	eqn_params->f_type = EURO_CALL_OPTION;
	    if (opt_type[5] == 'A' || opt_type[5] == 'a')
	    	eqn_params->f_type = AMRI_CALL_OPTION;
	}

	CursorAfterString(infile,"Enter strike price:");
	fscanf(infile,"%lf",&eqn_params->E);
	(void) printf("%f\n",eqn_params->E);
	CursorAfterString(infile,"Enter volatility:");
	fscanf(infile,"%lf",&eqn_params->sigma[0]);
	(void) printf("%f\n",eqn_params->sigma[0]);
	CursorAfterString(infile,"Enter interest rate:");
	fscanf(infile,"%lf",&eqn_params->r);
	(void) printf("%f\n",eqn_params->r);
	CursorAfterString(infile,"Enter dividend rate:");
	fscanf(infile,"%lf",&eqn_params->D);
	(void) printf("%f\n",eqn_params->D);

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
}	/* end read_finance_params */

extern  double price_of_state(
        POINTER state)
{
        STATE *price_state = (STATE*)state;
        return (double)(*price_state);
}       /* end price_of_state */


static  void excercise_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	RECT_GRID *top_grid = &topological_grid(front->grid_intfc);
	PARAMS *eqn_params = (PARAMS*)front->extra1;
        int i, id, dim = top_grid->dim;
	int ip0[MAXD],ip1[MAXD],ip2[MAXD],ip3[MAXD];
        double dn,*h = top_grid->h,*L = top_grid->L;
        double S;
        double *p0 = Coords(oldp);
	double *option_price = eqn_params->option_price;
	double epsilon, delta;
	const double EPS = 0.0000001;
	F_TYPE f_type = eqn_params->f_type;
        STATE *sl,*sr,*new_state,*old_state;
	FUNC_PARAMS eparams;
	double r,D,E,x;

	r = eqn_params->r;
	D = eqn_params->D;
	E = eqn_params->E;

	if (debugging("excercise_point"))
	    (void) printf("Entering excercise_point_propagate()\n");

        FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        if (wave_type(oldhs) < FIRST_PHYSICS_WAVE_TYPE)
        {
            for (i = 0; i < dim; ++i)
                Coords(newp)[i] = Coords(oldp)[i];
	    if (negative_component(oldhs) == BLACK_SCHOLES_COMP)
	    {
	    	new_state = (STATE*)left_state(newp);
            	*new_state = (STATE)(*sl);
	    }
	    else
	    {
	    	new_state = (STATE*)right_state(newp);
            	*new_state = (STATE)(*sr);
	    }
	    if (debugging("excercise_point"))
	    	(void) printf("Leaving excercise_point_propagate()\n");
            return;
        }
	if (front->time == 0.0)
	{
	    new_state = (STATE*)left_state(newp);
	    *new_state = (STATE)(E - Coords(newp)[0]);
	    new_state = (STATE*)right_state(newp);
	    *new_state = (STATE)(E - Coords(newp)[0]);
	    return;
	}

	rect_in_which(p0,ip0,top_grid);
	if (debugging("excercise_point"))
	{
	    printf("p0 = %f\n",p0[0]);
	    printf("ip0 = %d\n",ip0[0]);
	}
	ip1[0] = ip0[0] - 1;
	ip3[0] = ip0[0] + 2;
	eparams.S1 = L[0] + ip1[0]*h[0];
        eparams.S3 = L[0] + ip3[0]*h[0];
	eparams.E = E;
	if (f_type == AMRI_CALL_OPTION)
	{
	    if (debugging("excercise_point"))
		printf("American call option\n");
	    ip2[0] = ip0[0];
	    eparams.S2 = L[0] + ip2[0]*h[0];
	    id = d_index1d(ip1[0],top_gmax);
	    eparams.V1 = option_price[id];
	    id = d_index1d(ip2[0],top_gmax);
	    eparams.V2 = option_price[id];
	    id = d_index1d(ip3[0],top_gmax);
	    eparams.V3 = option_price[id]*(1.0 - r*dt) - dt*(r-D)*eparams.S3;
	    eparams.f_type = AMRI_CALL_OPTION;
	}
	else if (f_type == AMRI_PUT_OPTION)
	{
	    if (debugging("excercise_point"))
		printf("American put option\n");
	    ip2[0] = ip0[0] + 1;
	    eparams.S2 = L[0] + ip2[0]*h[0];
	    id = d_index1d(ip1[0],top_gmax);
	    eparams.V1 = option_price[id]*(1.0 - r*dt) - dt*(r-D)*eparams.S1;
	    id = d_index1d(ip2[0],top_gmax);
	    eparams.V2 = option_price[id];
	    id = d_index1d(ip3[0],top_gmax);
	    eparams.V3 = option_price[id];
	    eparams.f_type = AMRI_PUT_OPTION;
	}

	epsilon = delta = h[0]*EPS;

	if (!extract_func_root(&x,(POINTER)&eparams))
	{
	    char efname[100];
	    sprintf(efname,"extract_func%d.xg",front->step);
	    FILE *exfile = fopen(efname,"w");
	    double dS,f;
	    double a,b,c;
	    printf("S1 = %20.14f  V1 = %20.14f\n",eparams.S1,eparams.V1);
	    printf("S2 = %20.14f  V2 = %20.14f\n",eparams.S2,eparams.V2);
	    printf("S3 = %20.14f  V3 = %20.14f\n",eparams.S3,eparams.V3);
	    printf(" E = %20.14f\n",eparams.E);
	    a = eparams.a;
	    b = eparams.b;
	    c = eparams.c;
	    fprintf(exfile,"\"extract function\"\n");
	    dS = (eparams.S3 - eparams.S1)/200.0;
	    for (i = 0; i <= 200; ++i)
	    {
		S = eparams.S1 + i*dS;
		f = a*S*S + b*S + c;
	    	fprintf(exfile,"%f  %f\n",S,f);
	    }
	    fprintf(exfile,"\n\n");
	    fclose(exfile);
	    int count = 0;
	    while (count < 5)
	    {
		if (f_type == AMRI_PUT_OPTION)
		{
		    eparams.S3 += h[0];
		    ip3[0] += 1;
		    id = d_index1d(ip3[0],top_gmax);
		    eparams.V3 = option_price[id];
		    printf("Trying: S3 = %20.14f  V3 = %20.14f\n",
				eparams.S3,eparams.V3);
		}
		else if (f_type == AMRI_CALL_OPTION)
		{
		    eparams.S1 -= h[0];
		    ip1[0] -= 1;
		    id = d_index1d(ip1[0],top_gmax);
		    eparams.V1 = option_price[id];
		    printf("Trying: S1 = %20.14f  V1 = %20.14f\n",
				eparams.S1,eparams.V1);
		}
		if (!extract_func_root(&x,(POINTER)&eparams))
		{
		    printf("Solution find at:\n");
		    printf("S1 = %20.14f  V1 = %20.14f\n",
				eparams.S1,eparams.V1);
		    printf("S2 = %20.14f  V2 = %20.14f\n",
				eparams.S2,eparams.V2);
		    printf("S3 = %20.14f  V3 = %20.14f\n",
				eparams.S3,eparams.V3);
		    printf("x = %20.14f\n",x);
		    break;
		}
		count++;
	    }
	    if (count == 5)
	    	clean_up(ERROR);
	}
	Coords(newp)[0] = x;

	new_state = (STATE*)left_state(newp);
	if (f_type == AMRI_PUT_OPTION)
	    *new_state = (STATE)(E - Coords(newp)[0]);
	else if (f_type == AMRI_CALL_OPTION)
	    *new_state = (STATE)(Coords(newp)[0] - E);

	new_state = (STATE*)right_state(newp);
	if (f_type == AMRI_PUT_OPTION)
	    *new_state = (STATE)(E - Coords(newp)[0]);
	else if (f_type == AMRI_CALL_OPTION)
	    *new_state = (STATE)(Coords(newp)[0] - E);
}       /* solute_point_propagate */

extern double extend_from_put_exc(
	double P0,
	double *S,
	double *Sf,
	double *h,
	double E)
{
	return E - *Sf;
}	/* end extend_from_put_exc */

extern double extend_from_call_exc(
	double P0,
	double *S,
	double *Sf,
	double *h,
	double E)
{
	return *Sf - E;
}	/* end extend_from_call_exc */

static boolean extract_func(
	double S,
	double *fvalue,
	POINTER params)
{
	FUNC_PARAMS *func_params = (FUNC_PARAMS*)params;
	double S1 = func_params->S1;
	double S2 = func_params->S2;
	double S3 = func_params->S3;
	double V1 = func_params->V1;
	double V2 = func_params->V2;
	double V3 = func_params->V3;
	double D, Na, Nb, a, b, c;

	D  = (sqr(S2) - sqr(S1))*(S3 - S1) - (sqr(S3) - sqr(S1))*(S2 - S1);
	if (D == 0.0) return NO;
	Na = (V2 - V1)*(S3 - S1) - (V3 - V1)*(S2 - S1);
	Nb = (sqr(S2) - sqr(S1))*(V3 - V1) - (sqr(S3) - sqr(S1))*(V2 - V1);
	a = Na/D;		b = Nb/D;
	c = V2 - a*sqr(S2) - b*S2;
	if (func_params->f_type == AMRI_PUT_OPTION)
	    *fvalue = a*S*S + b*S + c - func_params->E + S;
	else if (func_params->f_type == AMRI_CALL_OPTION)
	    *fvalue = a*S*S + b*S + c - S + func_params->E;
	func_params->a = a;
	func_params->b = b;
	func_params->c = c;
	return YES;
}	/* end extract_func */	

static boolean extract_func_root(
	double *S,
	POINTER params)
{
	FUNC_PARAMS *func_params = (FUNC_PARAMS*)params;
	double S1 = func_params->S1;
	double S2 = func_params->S2;
	double S3 = func_params->S3;
	double V1 = func_params->V1;
	double V2 = func_params->V2;
	double V3 = func_params->V3;
	double D, Na, Nb, a, b, c;
	double arg;

	D  = (sqr(S2) - sqr(S1))*(S3 - S1) - (sqr(S3) - sqr(S1))*(S2 - S1);
	if (D == 0.0) return NO;
	Na = (V2 - V1)*(S3 - S1) - (V3 - V1)*(S2 - S1);
	Nb = (sqr(S2) - sqr(S1))*(V3 - V1) - (sqr(S3) - sqr(S1))*(V2 - V1);
	a = Na/D;		b = Nb/D;
	c = V2 - a*sqr(S2) - b*S2;
	if (func_params->f_type == AMRI_PUT_OPTION)
	{
	    b += 1.0;
	    c -= func_params->E;
	}
	else if (func_params->f_type == AMRI_CALL_OPTION)
	{
	    b -= 1.0;
	    c += func_params->E;
	}
	arg = b*b - 4.0*a*c;
	if (arg < 0.0) 
	{
	    func_params->a = a;
	    func_params->b = b;
	    func_params->c = c;
	    return NO;
	}
	S1 = (-b + sqrt(arg))/2.0/a;
	S2 = (-b - sqrt(arg))/2.0/a;
	if (S1 >= func_params->S1 && S1 <= func_params->S3) 
	    *S = S1;
	else if (S2 >= func_params->S1 && S2 <= func_params->S3) 
	    *S = S2;
	return YES;
}	/* end extract_func_root */	

static boolean erf_func(
	double x,
	double *y,
	POINTER params)
{
  
	FUNC_PARAMS *erf_params = (FUNC_PARAMS*)params;
	double r,D,E,t,sigma,temp;

	r = erf_params->r;
	D = erf_params->D;
	E = erf_params->E;
	t = erf_params->t;
	sigma = erf_params->sigma;
	temp = r - D + sqr(sigma)/2;
	
	double b1,b2,s_f;  
	b1 = x/(2*sqrt(t)) - temp/(2*sigma);
	b2 = sigma + b1;
	s_f = std::min(E,r*E/D)*exp(-temp*t-sigma*sqrt(t)*x);
	
         
        *y = -exp(-D*t)*(1+erf(x/sqrt(2)))/2 
	  + erf(sqrt((2*D+sqr(b1))*t/2))/sqrt(2*D+sqr(b1)) * (D/sigma - b1/2)
	  - r*E*erf(sqrt((2*r+sqr(b2))*t/2))/sigma/s_f/sqrt(2*r+sqr(b2))
	  - (1+exp(-D*t)*(-2+(1-erf(sqrt(t)*b1/sqrt(2)))))/2;

	return YES;
}	/* end erf_func */	

static void record_excer_point(
	Front *front,
	char *outname)
{
	INTERFACE *intfc = front->interf;
	POINT **p,*point;
	static FILE *ifile;
	static FILE *efile;
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    char fname[100];
	    sprintf(fname,"%s/intfc_numer",outname);
	    ifile = fopen(fname,"w");
	    sprintf(fname,"%s/intfc_exact",outname);
	    efile = fopen(fname,"w");
	}
	for (p = intfc->points; p && *p; ++p)
	{
	    point = *p;
	    if (wave_type(point) == FIRST_PHYSICS_WAVE_TYPE)
		break;
	}
	fprintf(ifile,"%f %f\n",front->time,Coords(point)[0]);
	fprintf(efile,"%f %f\n",front->time,exact_put_solution(front));
}	/* end record_excer_point */

static double exact_put_solution(
	Front *front)
{
	PARAMS *eqn_params = (PARAMS*)front->extra1;
	double time = front->time;
        int step = front->step;
	double epsilon, delta;
	const double EPS = 1e-6;
	static FUNC_PARAMS erf_params;
	double r, D, E, sigma,t;
	double y = -1.0;
	double s_f,temp,Sl,Su;
	static double x;

	epsilon = delta = EPS;
	r = erf_params.r = eqn_params->r;  //interest rate
	D = erf_params.D = eqn_params->D;  //dividend
	E = erf_params.E = eqn_params->E;  // strike price
	sigma = erf_params.sigma = eqn_params->sigma[0];  //volatility
	t = erf_params.t = time; // current time (time to expiry)

	if (x == 0.0)
	{
	    Sl = 2.0; 	Su = 3.0;
	}
	else
	{
	    Sl = x/2.0; Su = x*2.0;
	}
	    
	if (!find_root(erf_func,(POINTER)&erf_params,y,&x,Sl,Su,epsilon,delta))
	    x = 0.0;

	temp = r - D + sqr(sigma)/2;
	s_f = std::min(E,r*E/D)*exp(-temp*t-sigma*sqrt(t)*x);
	return s_f;
}       /* exact_put_solution */

static void setInitialIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack)
{
	int dim = front->rect_grid->dim;
	static double **point;
	PARAMS *eqn_params = (PARAMS*)front->extra1;
	double r,D,E;
	r = eqn_params->r;
	D = eqn_params->D;
	E = eqn_params->E;

	switch (dim)
	{
	case 1:
            FT_MatrixMemoryAlloc((POINTER*)&point,1,1,sizeof(double));
	    if (eqn_params->f_type == AMRI_CALL_OPTION)
            {
                point[0][0] = E;
                level_func_pack->num_points = 1;
                level_func_pack->point_array = point;
                level_func_pack->pos_component = EXCERCISE_COMP;
                level_func_pack->neg_component = BLACK_SCHOLES_COMP;
            }
            if (eqn_params->f_type == AMRI_PUT_OPTION)
            {
                point[0][0] = std::min(E,r*E/D);
                level_func_pack->num_points = 1;
                level_func_pack->point_array = point;
                level_func_pack->neg_component = EXCERCISE_COMP;
                level_func_pack->pos_component = BLACK_SCHOLES_COMP;
            }
	    break;
	case 2:
	    break;
	case 3:
	    break;
	}
	level_func_pack->wave_type = FIRST_PHYSICS_WAVE_TYPE;
}	/* setInitialIntfc */


static void read_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,k,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	double state;
	POINTER func_params;
	HYPER_SURF *hs;

	for (i = 0; i < dim; ++i)
	{
	    if (f_basic.boundary[i][0] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
	        if (rect_boundary_type(front->interf,i,0) == DIRICHLET_BOUNDARY)
		    hs = FT_RectBoundaryHypSurf(front->interf,
					DIRICHLET_BOUNDARY,i,0);
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
		    CursorAfterString(infile,"Enter option value:");
		    fscanf(infile,"%lf",&state);
		    (void) printf("%f\n",state);
		    FT_SetDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)&state,hs);
		    break;
		case 't':			// Time dependent state
		case 'T':
		    get_time_dependent_params(dim,infile,&func_params);
		    FT_SetDirichletBoundary(front,timeDependBoundaryState,
					"timeDependBoundaryState",
					func_params,NULL,hs);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
            if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
                if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)
                    hs = FT_RectBoundaryHypSurf(front->interf,
					DIRICHLET_BOUNDARY,i,1);
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
		    CursorAfterString(infile,"Enter option value:");
		    fscanf(infile,"%lf",&state);
		    (void) printf("%f\n",state);
		    FT_SetDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)&state,hs);
		    break;
		case 't':			// Time dependent state
		case 'T':
		    get_time_dependent_params(dim,infile,&func_params);
		    FT_SetDirichletBoundary(front,timeDependBoundaryState,
					"timeDependBoundaryState",
					func_params,NULL,hs);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_dirichlet_bdry_data */


static void timeDependBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	TIME_DEPENDENT_PARAMS *td_params = (TIME_DEPENDENT_PARAMS*)params;
	double *op_value = (double*)state;
	int i,dim = front->rect_grid->dim;
	double time = front->time;
	static int step = 0;
	printf("timeDependBoundaryState() needed!\n");
	clean_up(ERROR);

}	/* end timeDependBoundaryState */


static void get_time_dependent_params(
	int dim,
	FILE *infile,
	POINTER *params)
{
	static TIME_DEPENDENT_PARAMS *td_params;
	char string[100];
	int i;

	FT_ScalarMemoryAlloc((POINTER*)&td_params,
			sizeof(TIME_DEPENDENT_PARAMS));
	CursorAfterString(infile,"Enter type of time-dependent function:");
	fscanf(infile,"%s",string);
	(void) printf(" %s\n",string);
	printf("get_time_dependent_params() needed!\n");
	clean_up(ERROR);

	*params = (POINTER)td_params;	
}	/* end get_time_dependent_params */

