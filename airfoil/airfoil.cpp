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
*				example16.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a multiple vortex field. It 
*	demonstrates the reversibility of the front tracking method.
*
*/

#include <iFluid.h>
#include <airfoil.h>

	/*  Function Declarations */
static void airfoil_driver(Front*,Incompress_Solver_Smooth_Basis*);
static void airfoil_point_propagate(Front*,POINTER,POINT*,POINT*,
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void elastic_point_propagate(Front*,POINTER,POINT*,POINT*,
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void InitLocalIO(int,char**);
static void print_total_length(Front*,char*);
static int test_vortex_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                HYPER_SURF*,double*);
static int test_double_vortex_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                        HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
boolean RestartRun;
int RestartStep;
boolean binary = YES;
int constrained_propagate;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

static void test_tan_curve_propagate(Front*,Front*,INTERFACE*,
				CURVE*,CURVE*,double);
static void test_tan_curve_propagate1(Front*,Front*,INTERFACE*,
				CURVE*,CURVE*,double);
static void set_intfc_bond_length0(INTERFACE*);
static void print_compare_bond_length(INTERFACE*);
static void xgraph_front(Front*,char*);
static void FrontCurveSegLengthConstr(CURVE*,BOND*,BOND*,int,double,
			REDISTRIBUTION_DIRECTION);
static void forward_curve_seg_len_constr(CURVE*,BOND*,BOND*,int,double);
static void backward_curve_seg_len_constr(CURVE*,BOND*,BOND*,int,double);

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	IF_PARAMS 	iFparams;
	AF_PARAMS 	af_params;
	AF_PROB_TYPE 	prob_type;

	FT_Init(argc,argv,&f_basic);
	f_basic.size_of_intfc_state = sizeof(STATE);

	//Initialize Petsc before FrontStartUP
        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
        //if (debugging("trace")) printf("Passed PetscInitialize()\n");

	/*Construct Incompress Solver l_cartesian*/

	Incompress_Solver_Smooth_Basis *l_cartesian;
	if(f_basic.dim == 2)
	    l_cartesian = new Incompress_Solver_Smooth_2D_Cartesian(front);
	else if(f_basic.dim == 3)
	    l_cartesian = new Incompress_Solver_Smooth_3D_Cartesian(front);

	/* Initialize basic computational data */

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

	if (debugging("trace")) printf("Passed FT_StartUp()\n");

        iFparams.dim = f_basic.dim;
        front.extra1 = (POINTER)&iFparams;
        front.extra2 = (POINTER)&af_params;
	read_af_prob_type(in_name,&prob_type);
        read_iFparams(in_name,&iFparams);
        if (debugging("trace")) printf("Passed read_iFparams()\n");
        read_iF_movie_options(in_name,&iFparams);
        if (debugging("trace")) printf("Passed read_iF_movie_options()\n");

	setInitialIntfc(&front,&level_func_pack,in_name,prob_type);

	if (!RestartRun)
	{
	    if (f_basic.dim == 3)
                level_func_pack.set_3d_bdry = YES;
	    FT_InitIntfc(&front,&level_func_pack);
	    read_iF_dirichlet_bdry_data(in_name,&front,f_basic);
	    if (f_basic.dim < 3)
            	FT_ClipIntfcToSubdomain(&front);
	}

	/* Time control */
	FT_ReadTimeControl(in_name,&front);

	/* Initialize velocity field function */

	initVelocityFunc(&front,prob_type);

	front.tan_curve_propagate = test_tan_curve_propagate1;
	front._point_propagate = airfoil_point_propagate;

	l_cartesian->initMesh();
	if (debugging("sample_velocity"))
            l_cartesian->initSampleVelocity(in_name);
        init_fluid_state_func(l_cartesian,prob_type);
        if (debugging("trace"))
            printf("Passed l_cartesian->initMesh()\n");
        if (RestartRun)
            l_cartesian->readFrontInteriorStates(restart_state_name);
        else
            l_cartesian->setInitialCondition();
        if (debugging("trace"))
            printf("Passed state initialization()\n");

	/* Propagate the front */

	set_intfc_bond_length0(front.interf);
	airfoil_driver(&front,l_cartesian);

	clean_up(0);
}

static  void airfoil_driver(
        Front *front,
	Incompress_Solver_Smooth_Basis *l_cartesian)
{
        int status,count;
        Front *newfront;
        double dt,dt_frac,CFL;
        char s[10];
        double fcrds[MAXD];
        int  dim = front->rect_grid->dim;

        CFL = Time_step_factor(front);
	Frequency_of_redistribution(front,GENERAL_WAVE) = 100;

	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
	    FT_ResetTime(front);

	    // Always output the initial interface.
	    if (dim == 2)
	    {
	    	xgraph_front(front,out_name);
	    	print_total_length(front,out_name);
	    }

	    FT_Propagate(front);
            if (debugging("trace")) printf("Calling ifluid solve()\n");

            FT_SetOutputCounter(front);
	    FT_SetTimeStep(front);
	    front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
	}
	else
	{
	    FT_SetOutputCounter(front);
	}

	FT_TimeControlFilter(front);

	if (debugging("step_size"))
                printf("Time step from start: %f\n",front->dt);
        for (;;)
        {
	    /* Propagating interface for time step dt */

	    if (debugging("trace"))
                printf("Before FT_Propagate()\n");
            FT_Propagate(front);
            if (debugging("trace")) printf("Passed FT_Propagate()\n");

            if (debugging("trace")) printf("Calling ifluid solve()\n");
            l_cartesian->solve(front->dt);
            if (debugging("trace")) printf("Passed ifluid solve()\n");
	    if (debugging("trace"))
            {
                printf("After solve()\n");
                print_storage("at end of time step","trace");
            }

	    FT_AddTimeStepToCounter(front);

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

	    FT_SetTimeStep(front);
            if (debugging("step_size"))
                printf("Time step from FrontHypTimeStep(): %f\n",front->dt);
            front->dt = std::min(front->dt,CFL*l_cartesian->max_dt);
            if (debugging("step_size"))
                printf("Time step from l_cartesian->max_dt(): %f\n",front->dt);

	    /* Output section */

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);
	    if (dim == 2)
	    	print_total_length(front,out_name);

            if (FT_IsSaveTime(front))
	    {
		FT_Save(front,out_name);
                l_cartesian->printFrontInteriorStates(out_name);
	    }
	    if (debugging("trace"))
                printf("After print output()\n");
            if (FT_IsMovieFrameTime(front))
	    {
		if (debugging("trace"))
                    printf("Calling initMovieVariable()\n");
                l_cartesian->initMovieVariables();
                if (debugging("trace"))
                    printf("Calling FT_AddMovieFrame()\n");
                FT_AddMovieFrame(front,out_name,binary);
	    }

            if (FT_TimeLimitReached(front))
                    break;

	    /* Time and step control section */

	    FT_TimeControlFilter(front);
        }
        (void) delete_interface(front->interf);
}       /* end airfoil_driver */


static void test_tan_curve_propagate1(
	Front           *fr,
        Front           *newfr,
        INTERFACE       *intfc,
        CURVE           *oldc,
        CURVE           *newc,
        double           dt)
{
	BOND *b,*bs,*bs2,*be;
	int nb,n;
	double seg_len,total_length,total_length0;
	static boolean first = YES;

	if (fr->rect_grid->dim != 2) return;

	if (wave_type(newc) != ELASTIC_BOUNDARY)
	    return;

	if (first)
	{
	    first = NO;
	    nb = 0;
	    total_length  = 0.0;
	    for (b = newc->first; b != NULL; b = b->next)
	    {
	    	total_length += bond_length(b);
		nb++;
	    }
	    for (b = newc->first; b != NULL; b = b->next)
		b->length0 = total_length/(double)nb;
	}
	if (debugging("airfoil"))
	{
	    printf("Entering test_tan_curve_propagate1()\n");
	    nb = 0;
	    total_length  = 0.0;
	    total_length0 = 0.0;
	    for (b = newc->first; b != NULL; b = b->next)
	    {
	    	total_length += bond_length(b);
	    	total_length0 += b->length0;
		nb++;
	    }
	    printf("Entering: total_length  = %20.14f\n", total_length);
	    printf("Entering: total_length0 = %20.14f\n", total_length0);
	}

	nb = 0;
	for (b = newc->first; b != NULL; b = b->next) nb++;
	if (nb%2 == 0)	n = nb/2;
	else n = nb/2 + 1;

	seg_len = 0.0;
	bs = newc->first;
	for (nb = 0, b = bs; nb < n; b = b->next) 
	{
	    nb++;
	    seg_len += b->length0;
	    be = b;
	}
	bs2 = be->next;
	FrontCurveSegLengthConstr(newc,bs,be,nb,seg_len,
				BACKWARD_REDISTRIBUTION);

	seg_len = 0.0;
	bs = bs2;
	for (nb = 0, b = bs; b != NULL; b = b->next) 
	{
	    nb++;
	    seg_len += b->length0;
	    be = b;
	}
	FrontCurveSegLengthConstr(newc,bs,be,nb,seg_len,
				FORWARD_REDISTRIBUTION);
	
	total_length = 0.0;
	for (b = newc->first; b != NULL; b = b->next)
	{
	    total_length += bond_length(b);
	    b->length0 = seg_len/(double)nb;
	}
	if (debugging("airfoil"))
	{
	    printf("Exiting : %20.14f\t", total_length);
	    printf("Leaving test_tan_curve_propagate()\n");
	}
}	/* end test_tan_curve_propagate1 */

EXPORT void FrontCurveSegLengthConstr(
	CURVE *c,
	BOND *bs,
	BOND *be,
	int nbds,
	double seg_length,
	REDISTRIBUTION_DIRECTION dir)
{
	switch (dir)
	{
	case FORWARD_REDISTRIBUTION:
	    return forward_curve_seg_len_constr(c,bs,be,nbds,seg_length);
	case BACKWARD_REDISTRIBUTION:
	    return backward_curve_seg_len_constr(c,bs,be,nbds,seg_length);
	}
}	/* end FrontCurveSegLengthConstr */

static void test_tan_curve_propagate(
	Front           *fr,
        Front           *newfr,
        INTERFACE       *intfc,
        CURVE           *oldc,
        CURVE           *newc,
        double           dt)
{
	char xname[100];
	FILE *xfile;
	BOND *b,*oldb;
	double *dist_l,*dist_r;
	double **tan_l,**tan_r;
	POINT *p,*p_nb,**pts;
	int i,n,np,k;
	int dim = Dimension(intfc);
	double lambda;
	double **dp;
	double total_length,fixed_length;
	
	if (wave_type(newc) != ELASTIC_BOUNDARY)
	    return;

	if (debugging("airfoil"))
	    printf("Entering test_tan_curve_propagate()\n");
	lambda = 0.0;

	np = NumOfCurvePoints(newc);
	FT_VectorMemoryAlloc((POINTER*)&dist_l,np,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&dist_r,np,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&tan_l,np,MAXD,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&tan_r,np,MAXD,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&dp,np,MAXD,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&pts,np,sizeof(POINT*));

	//calculate the total length of the newc 
	total_length = 0.0;
	for (b = newc->first; b != NULL; b = b->next)
	    total_length += bond_length(b);
	printf("Entering: %24.18g\t", total_length);

	b = newc->first;
	p = pts[0] = b->start;
	p_nb = b->end;
	dist_l[0] = 0.0;
	dist_r[0] = bond_length(b) - b->length0;
	for (i = 0; i < dim; ++i)
	{
	    tan_r[0][i] = Coords(p_nb)[i] - Coords(p)[i];
	    dp[0][i] = lambda*dist_r[0]*tan_r[0][i];
	}

	n = 1;
	for (b = newc->first; b != newc->last; b = b->next)
	{
	    p = pts[n] = b->end; 	
	    p_nb = b->start;
	    dist_l[n] = bond_length(b) - b->length0;
	    for (i = 0; i < dim; ++i)
		tan_l[n][i] = Coords(p_nb)[i] - Coords(p)[i];
	    p_nb = b->next->end;
	    dist_r[n] = bond_length(b->next) - b->next->length0;
	    for (i = 0; i < dim; ++i)
	    {
		tan_r[n][i] = Coords(p_nb)[i] - Coords(p)[i];
	    	dp[n][i] = lambda*(dist_l[n]*tan_l[n][i] +
				   dist_r[n]*tan_r[n][i]);
	    }
	    ++n;
	}
	b = newc->last;
	p = pts[n] = b->end;
	p_nb = b->start;
	dist_l[n] = bond_length(b) - b->length0;
	dist_r[n] = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    tan_l[n][i] = Coords(p_nb)[i] - Coords(p)[i];
	    dp[n][i] = lambda*dist_l[n]*tan_l[n][i];
	}

	k = np / 2;
	b = newc->first;
	for (n = 1; n < k; ++n)
	    b = b->next;
	for (n = 1; n <= k; ++n)
	{
	    for (i = 0; i < dim; ++i)
		lambda += tan_l[k-n+1][i] * tan_l[k-n+1][i];
	    lambda = sqrt(lambda);
	    for (i = 0; i < dim; ++i)
		Coords(pts[k-n])[i] = Coords(pts[k-n+1])[i] + 
				b->length0*tan_l[k-n+1][i]/lambda;
	    lambda = 0.0;
	    b = b->prev;
	}

	b = newc->first;
	for (n = 0; n < k; ++n)
	    b = b->next;
	for (n = 1; n <= k; ++n)
	{
	    for (i = 0; i < dim; ++i)
		lambda += tan_r[k+n-1][i] * tan_r[k+n-1][i];
	    lambda = sqrt(lambda);
	    for (i = 0; i < dim; ++i)
	    	Coords(pts[k+n])[i] = Coords(pts[k+n-1])[i] + 
			b->length0*tan_r[k+n-1][i]/lambda;
	    lambda = 0.0;
	    b = b->next;
	}

	fixed_length = total_length = 0.0;
	for (b = newc->first, oldb = oldc->first; b != NULL; 
		b = b->next, oldb = oldb->next)
	{
	    set_bond_length(b,dim);
	    total_length += bond_length(b);
	    fixed_length += b->length0;
	}

	printf("Exiting: %24.18g\t", total_length);

	if (fr->step%50 == 0)
	{
	    sprintf(xname,"curves-%d",fr->step);
	    xfile = fopen(xname,"w");
	    xgraph_curve(xfile,newc,XY_PLANE);
	    xgraph_curve(xfile,newc,XY_PLANE);
	    fclose(xfile);
	}
	if (debugging("airfoil"))
	{
	    printf("total_length = %f  fixed_length = %f\n",
			total_length,fixed_length);
	    printf("Leaving test_tan_curve_propagate()\n");
	}
	FT_FreeThese(6,dist_l,dist_r,tan_l,tan_r,dp,pts);
}	/* end test_tan_curve_propagate */

static void set_intfc_bond_length0(
	INTERFACE *intfc)
{
	CURVE **c,*curve;
	BOND *b;

	if (intfc->dim == 3) return;

	for (c = intfc->curves; c && *c; ++c)
	{
	    curve = *c;
	    for (b = curve->first; b != NULL; b = b->next)
		b->length0 = bond_length(b);
	}
}

void airfoil_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	if (wave_type(oldhs) == ELASTIC_BOUNDARY)
	    return elastic_point_propagate(front,wave,oldp,newp,oldhse,oldhs,
					dt,V);
	else
	    return ifluid_point_propagate(front,wave,oldp,newp,oldhse,oldhs,
					dt,V);
}       /* airfoil_point_propagate */

void elastic_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	STATE *oldst,*newst;
	POINTER sl,sr;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        IF_FIELD *field = iFparams->field;
	int i, dim = front->rect_grid->dim;
	double *m_vor = field->vort;

	fourth_order_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	newst = (STATE*)left_state(newp);
	for (i = 0; i < dim; ++i)
	    newst->vel[i] = V[i];
	newst = (STATE*)right_state(newp);
	for (i = 0; i < dim; ++i)
	    newst->vel[i] = V[i];
	if (dim == 2)
        {
            if (!FT_IntrpStateVarAtCoords(front,-1,Coords(oldp),m_vor,
				getStateVort,&newst->vort))
                newst->vort = oldst->vort;
        }
}       /* elastic_point_propagate */

void print_compare_bond_length(
	INTERFACE *intfc)
{
	CURVE **c,*curve;
	BOND *b;

	for (c = intfc->curves; c && *c; ++c)
	{
	    curve = *c;
	    if (wave_type(curve) != ELASTIC_BOUNDARY)
		continue;
	    for (b = curve->first; b != NULL; b = b->next)
	    	printf("old length = %f  new length = %f\n",
				b->length0,bond_length(b));
	}
}

static void xgraph_front(
	Front *front,
	char *outname)
{
	char fname[100];
	sprintf(fname,"%s/intfc-%s",outname,right_flush(front->step,4));
	xgraph_2d_intfc(fname,front->interf);
}	/* end xgraph_front */

static void print_total_length(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
	CURVE **c,*curve;
	double total_length;
	BOND *b;
	char fname[100];
	static FILE *lfile;

	if (lfile == NULL)
	{
	    sprintf(fname,"%s/curve_length",out_name);
	    lfile = fopen(fname,"w");
	    fprintf(lfile,"\"Curve_length vs. time\"\n");
	}
	curve = NULL;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
	    {
		curve = *c;
		break;
	    }	
	}
	if (curve == NULL)
	    return;
	total_length = 0.0;
	for (b = curve->first; b != NULL; b = b->next)
	    total_length += bond_length(b);
	fprintf(lfile,"%f %f\n",front->time,total_length);
	fflush(lfile);
}	/* end print_total_length */


static	void forward_curve_seg_len_constr(
	CURVE		*c,
	BOND		*bs,
	BOND		*be,
	int		nbds,
	double		seg_len)
{
	BOND		*b, *bstart, *bend;
	double		b_len, sc_len, offset, total_offset, s, oms;
	double		coords[MAXD];
	double		space_tol;
	boolean		reset_bend;

	b_len = seg_len/(double)nbds;
	space_tol = b_len*1.0e-8;


	offset = b_len;		total_offset = seg_len;
	bstart = bs;		bend = be;
	while (bstart != bend)
	{
	    b = bstart;
	    while ((sc_len = bond_length(b)) < offset - space_tol)
	    {
	    	offset -= sc_len;
		if (b == bend) break;
	    	b = b->next;
	    }
	    if (sc_len > offset + space_tol)
	    {
	    	s = offset/sc_len;	oms = 1.0 - s;
	    	coords[0] = oms*Coords(b->start)[0] + s*Coords(b->end)[0];
	    	coords[1] = oms*Coords(b->start)[1] + s*Coords(b->end)[1];
		if (b == bend) reset_bend = YES;
		else reset_bend = NO;
	    	if (insert_point_in_bond(Point(coords),b,c) !=
		    FUNCTION_SUCCEEDED)
	    	{
	    	    screen("ERROR in forward_curve_seg_len_constr(), "
	    	           "insert_point_in_bond failed\n");
	    	    clean_up(ERROR);
	    	}
		if (reset_bend) bend = b->next;
	    }
	    replace_curve_seg_by_bond(c,bstart,b);
	    total_offset -= b_len;
	    bstart = bstart->next;
	    offset = b_len;
	    if (total_offset < b_len - space_tol) break;
	}
	if (bstart == bend)
	{
	    if (total_offset < space_tol)
	    {
		Coords(bend->end)[0] = Coords(bstart->start)[0];
		Coords(bend->end)[1] = Coords(bstart->start)[1];
		bstart = bstart->prev;
	    	replace_curve_seg_by_bond(c,bstart,bend);
	    }
	    else
	    {
	    	sc_len = bond_length(bend);
	    	s = total_offset/sc_len;      oms = 1.0 - s;
	    	Coords(bend->end)[0] = Coords(bend->start)[0] + s*
			(Coords(bend->end)[0] - Coords(bend->start)[0]);
	    	Coords(bend->end)[1] = Coords(bend->start)[1] + s*
			(Coords(bend->end)[1] - Coords(bend->start)[1]);
	    	set_bond_length(bend,2);
		b = bend;
		while (total_offset/b_len > 2.0 - space_tol)
		{
	    	    s = b_len/total_offset;	oms = 1.0 - s;
	    	    coords[0] = oms*Coords(b->start)[0] + s*Coords(b->end)[0];
	    	    coords[1] = oms*Coords(b->start)[1] + s*Coords(b->end)[1];
	    	    if (insert_point_in_bond(Point(coords),b,c) !=
		    	FUNCTION_SUCCEEDED)
	    	    {
	    	    	screen("ERROR in forward_curve_seg_len_constr(), "
	    	           	"insert_point_in_bond failed\n");
	    	    	clean_up(ERROR);
	    	    }
		    total_offset -= b_len;
		    b = b->next;
		}
	    }
	}
	else
	{
	    if (total_offset < space_tol && bstart != NULL)
	    {
		Coords(bend->end)[0] = Coords(bstart->start)[0];
		Coords(bend->end)[1] = Coords(bstart->start)[1];
		bstart = bstart->prev;
	    	replace_curve_seg_by_bond(c,bstart,bend);
	    }
	}
	return;
}		/*end forward_curve_seg_len_constr*/

static	void backward_curve_seg_len_constr(
	CURVE		*c,
	BOND		*bs,
	BOND		*be,
	int		nbds,
	double		seg_len)
{
	BOND		*b, *bstart, *bend;
	double		b_len, sc_len, offset, total_offset, s, oms;
	double		coords[MAXD];
	double		space_tol;
	boolean		reset_bend;

	b_len = seg_len/(double)nbds;
	space_tol = b_len*1.0e-8;


	offset = b_len;		total_offset = seg_len;
	bstart = bs;		bend = be;
	while (bstart != bend)
	{
	    b = bend;
	    while ((sc_len = bond_length(b)) < offset - space_tol)
	    {
	    	offset -= sc_len;
		if (b == bstart) break;
	    	b = b->prev;
	    }
	    if (sc_len > offset + space_tol)
	    {
	    	s = offset/sc_len;	oms = 1.0 - s;
	    	coords[0] = oms*Coords(b->end)[0] + s*Coords(b->start)[0];
	    	coords[1] = oms*Coords(b->end)[1] + s*Coords(b->start)[1];
		if (b == bend) reset_bend = YES;
		else reset_bend = NO;
	    	if (insert_point_in_bond(Point(coords),b,c) !=
		    FUNCTION_SUCCEEDED)
	    	{
	    	    screen("ERROR in backward_curve_seg_len_constr(), "
	    	           "insert_point_in_bond failed\n");
	    	    clean_up(ERROR);
	    	}
		b = b->next;
		if (reset_bend) bend = b;
	    }
	    replace_curve_seg_by_bond(c,b,bend);
	    total_offset -= b_len;
	    bend = b->prev;
	    offset = b_len;
	    if (total_offset < b_len - space_tol) break;
	}
	if (bstart == bend)
	{
	    if (total_offset < space_tol)
	    {
		Coords(bstart->start)[0] = Coords(bstart->end)[0];
		Coords(bstart->start)[1] = Coords(bstart->end)[1];
	    	replace_curve_seg_by_bond(c,bstart,bstart->next);
	    }
	    else
	    {
	    	sc_len = bond_length(bstart);
	    	s = total_offset/sc_len;      oms = 1.0 - s;
	    	Coords(bstart->start)[0] = Coords(bend->end)[0] + s*
			(Coords(bend->start)[0] - Coords(bend->end)[0]);
	    	Coords(bstart->start)[1] = Coords(bend->end)[1] + s*
			(Coords(bend->start)[1] - Coords(bend->end)[1]);
	    	set_bond_length(bstart,2);
		b = bstart;
		while (total_offset/b_len > 2.0 - space_tol)
		{
	    	    s = b_len/total_offset;	oms = 1.0 - s;
	    	    coords[0] = oms*Coords(b->start)[0] + s*Coords(b->end)[0];
	    	    coords[1] = oms*Coords(b->start)[1] + s*Coords(b->end)[1];
	    	    if (insert_point_in_bond(Point(coords),b,c) !=
		    	FUNCTION_SUCCEEDED)
	    	    {
	    	    	screen("ERROR in forward_curve_seg_len_constr(), "
	    	           	"insert_point_in_bond failed\n");
	    	    	clean_up(ERROR);
	    	    }
		    total_offset -= b_len;
		}
	    }
	}
	else
	{
	    if (total_offset < space_tol && bend != NULL)
	    {
		Coords(bstart->start)[0] = Coords(bend->end)[0];
		Coords(bstart->start)[1] = Coords(bend->end)[1];
		bend = bend->next;
	    	replace_curve_seg_by_bond(c,bstart,bend);
	    }
	}
	return;
}		/*end backward_curve_seg_len_constr*/

