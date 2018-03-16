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


#include "iFluid.h"
#include "ifluid_basic.h"

	/*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
//for vd
static void neumann_point_propagate_vd(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate_vd(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);


static void rgbody_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static double (*getStateVel[MAXD])(POINTER) = {getStateXvel,getStateYvel,
					getStateZvel};
static double (*getStateVort3d[MAXD])(POINTER) = {getStateXvort,getStateYvort,
					getStateZvort};
static void iF_flowThroughBoundaryState2d(double*,HYPER_SURF*,Front*,POINTER,
        		POINTER);
static void iF_flowThroughBoundaryState3d(double*,HYPER_SURF*,Front*,POINTER,
        		POINTER);
static void get_time_dependent_params(int,FILE*,POINTER*);
static void addToEnergyFlux(RECT_GRID*,HYPER_SURF*,double*,double*,int,int,
			boolean);


extern double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

extern double getStateVort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort;
}	/* end getStateVort */

extern double getStateXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[0];
}	/* end getStateXvel */

extern double getStateYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[1];
}	/* end getStateYvel */

extern double getStateZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[2];
}	/* end getStateZvel */

extern double getStateXvort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort3d[0];
}	/* end getStateXvort */

extern double getStateYvort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort3d[1];
}	/* end getStateYvort */

extern double getStateZvort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort3d[2];
}	/* end getStateZvort */

// for vd
extern double getStateDens(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->dens;
}       /* end getStateDens */

// for vd
extern double getStateDensOld(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->dens_old;
}       /* end getStateDens */

// for vd
extern double getStateConc(POINTER state)
{
        STATE *fstate = (STATE*)state;
        return fstate->conc;
}       /* end getStateConc */

extern void read_iF_dirichlet_bdry_data(
	char *inname,
	Front *front,
	F_BASIC_DATA f_basic)
{
	char msg[100],s[100];
	int i,k,dim = front->rect_grid->dim;
	FILE *infile = fopen(inname,"r");
	static STATE *state;
	POINTER func_params;
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
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf",&state->vel[k]);
			(void) printf("%f ",state->vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state->pres);
		    (void) printf("%f\n",state->pres);
		    FT_SetDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_SetDirichletBoundary(front,iF_flowThroughBoundaryState,
					"iF_flowThroughBoundaryState",NULL,
					NULL,hs);
		    break;
		case 't':			// Flow through state
		case 'T':
		    get_time_dependent_params(dim,infile,&func_params);
		    FT_SetDirichletBoundary(front,iF_timeDependBoundaryState,
					"iF_timeDependBoundaryState",
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
		    FT_ScalarMemoryAlloc((POINTER*)&state,sizeof(STATE));
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf ",&state->vel[k]);
			(void) printf("%f ",state->vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state->pres);
		    (void) printf("%f\n",state->pres);
		    FT_SetDirichletBoundary(front,NULL,NULL,
					NULL,(POINTER)state,hs);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_SetDirichletBoundary(front,iF_flowThroughBoundaryState,
					"iF_flowThroughBoundaryState",NULL,
					NULL,hs);
		    break;
		case 't':			// Flow through state
		case 'T':
		    get_time_dependent_params(dim,infile,&func_params);
		    FT_SetDirichletBoundary(front,iF_timeDependBoundaryState,
					"iF_timeDependBoundaryState",
					func_params,NULL,hs);
		    break;
		default:
		    (void) printf("Unknown Dirichlet boundary!\n");
		    clean_up(ERROR);
		}
	    }
	}
	fclose(infile);
}	/* end read_iF_dirichlet_bdry_data */

extern void restart_set_dirichlet_bdry_function(Front *front)
{
	INTERFACE *intfc = front->interf;
	int i;
	BOUNDARY_STATE  *bstate;
	const char *s;
	for (i = 0; i < num_bstates(intfc); ++i)
	{
	    bstate = bstate_list(intfc)[i];
	    if (bstate == NULL) continue;
	    if (bstate->_boundary_state != NULL)
	    {
		STATE *bst = (STATE*)bstate->_boundary_state;
	    }
	    s = bstate->_boundary_state_function_name;
	    if (s == NULL) continue;
	    if (strcmp(s,"iF_flowThroughBoundaryState") == 0)
            	bstate->_boundary_state_function = iF_flowThroughBoundaryState;
	}
}	/* end restart_set_dirichlet_bdry_function */

extern void iF_timeDependBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	TIME_DEPENDENT_PARAMS *td_params = (TIME_DEPENDENT_PARAMS*)params;
	STATE *iF_state = (STATE*)state;
	int i,dim = front->rect_grid->dim;
	double time = front->time;
	double *T = td_params->T;
	double omega = td_params->omega;
	double phase = td_params->phase;
	static int step = 0;

	switch (td_params->td_type)
	{
	case CONSTANT:
	    for (i = 0; i < dim; ++i)
	    	iF_state->vel[i] = td_params->v_base[i];
	    iF_state->pres = td_params->p_base;
	    break;
	case PULSE_FUNC:
	    if (time <= T[0])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_base[i];
	    	iF_state->pres = td_params->p_base;
	    }
	    else if (time > T[0] && time <= T[1])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_base[i] +
				(time - T[0])/(T[1] - T[0])*
				(td_params->v_peak[i] - td_params->v_base[i]);
	    	iF_state->pres = td_params->p_base +
				(time - T[0])/(T[1] - T[0])*
				(td_params->p_peak - td_params->p_base);
	    }
	    else if (time > T[1] && time <= T[2])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_peak[i];
	    	iF_state->pres = td_params->p_peak;
	    }
	    else if (time > T[2] && time <= T[3])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_peak[i] +
				(time - T[2])/(T[3] - T[2])*
				(td_params->v_tail[i] - td_params->v_peak[i]);
	    	iF_state->pres = td_params->p_peak +
				(time - T[2])/(T[3] - T[2])*
				(td_params->p_tail - td_params->p_peak);
	    }
	    if (time > T[3])
	    {
	    	for (i = 0; i < dim; ++i)
	    	    iF_state->vel[i] = td_params->v_tail[i];
	    	iF_state->pres = td_params->p_tail;
	    }
	    break;
	case SINE_FUNC:
	    for (i = 0; i < dim; ++i)
	    	iF_state->vel[i] = td_params->v_base[i] +
			td_params->v_amp[i]*sin(omega*time + phase);
	    iF_state->pres = td_params->p_base + td_params->p_amp*
				sin(omega*time + phase);
	    break;
	default:
	    (void) printf("In iF_timeDependBoundaryState(), unknown type!\n");
	    clean_up(ERROR);
	}
	if (debugging("time_depend_bdry"))
	{
	    if (step != front->step)
	    {
	    	printf("time = %f  vel = %f %f  p = %f\n",time,iF_state->vel[0],
				front->time,iF_state->vel[1],iF_state->pres);
	    	step = front->step;
	    }
	}
}	/* end iF_timeDependBoundaryState */

extern void iF_flowThroughBoundaryState(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    return iF_flowThroughBoundaryState2d(p0,hs,front,params,state);
	case 3:
	    return iF_flowThroughBoundaryState3d(p0,hs,front,params,state);
	}
}	/* end iF_flowThroughBoundaryState */

static void iF_flowThroughBoundaryState3d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double vort3d[3][MAXD];		/* vorticity stencil */
	double pres[3];		/* pressure stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_vort3d[MAXD];		/* vort flux */
	double f_pres;		/* pressure flux */
	double dn,dt = front->dt;
	STATE *oldst,*newst = (STATE*)state;
	STATE  **sts;
	POINTER sl,sr;
	int i,j,k,dim = front->rect_grid->dim;
	int nrad = 2;

	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState3d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
	oldst = NULL;
	if (comp == negative_component(hs))
	    oldst = (STATE*)sl;
	else
	    oldst = (STATE*)sr;

	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);

	u[1] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    pres[j] = oldst->pres;
	    for (i = 0; i < dim; ++i)
	    	vort3d[j][i] = oldst->vort3d[i];
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			field->vel[i],getStateVel[i],&vtmp))
            {
        	vtmp = oldst->vel[i];
            }
	    u[1] += vtmp*dir[i];
	    newst->vel[i] = vtmp;
	    if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			field->vort3d[i],getStateVort3d[i],&vort3d[2][i]))
	    	vort3d[2][i] = oldst->vort3d[i];
	}
	if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2]))
	    pres[2] = oldst->pres;

	for (i = 0; i < dim; ++i)
	{
	    f_vort3d[i] = linear_flux(u[1],vort3d[i][0],vort3d[i][1],
				vort3d[i][2]);
	    newst->vort3d[i] = oldst->vort3d[i] - dt/dn*f_vort3d[i];
	}
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	newst->pres = oldst->pres - dt/dn*f_pres;
	if (debugging("flow_through"))
	{
	    (void) print_Nor_stencil(front,nsten);
	    (void) printf("new velocity after normal prop: %f %f %f\n",
			newst->vel[0],newst->vel[1],newst->vel[2]);
	}

	tsten = FrontGetTanStencils(front,oldp,nrad);
	if (tsten == NULL) return;

	for (k = 0; k < dim-1; ++k)
	{
	    for (i = 0; i < dim; ++i)
	    	dir[i] = tsten[k]->dir[i];
	    dn = FT_GridSizeInDir(dir,front);

	    if (comp == negative_component(hs))
	    	sts = (STATE**)tsten[k]->leftst;
	    else
	    	sts = (STATE**)tsten[k]->rightst;

	    if (debugging("flow_through"))
	    {
	    	(void) printf("Ambient component: %d\n",comp);
	    	(void) printf("Tangential grid size = %f\n",dn);
		(void) printf("For direction %d\n",k);
	    	(void) print_Tan_stencil(front,tsten[k]);
	    }

	    for (j = 0; j < 3; ++j)
	    	u[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
	    	pres[j] = sts[j-1]->pres;
	    	for (i = 0; i < dim; ++i)
		    u[j] += sts[j-1]->vel[i]*dir[i];
	    	for (i = 0; i < dim; ++i)
	    	{
		    v[j][i] = sts[j-1]->vel[i] - u[j]*dir[i];
	    	    vort3d[j][i] = sts[j-1]->vort3d[i];
	    	}
	    }

	    f_u = burger_flux(u[0],u[1],u[2]);
	    for (i = 0; i < dim; ++i)
	    {
	    	f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	        f_vort3d[i] = linear_flux(u[1],vort3d[0][i],vort3d[1][i],
					vort3d[2][i]);
	    }
	    f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	    for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
	    	newst->vort3d[i] += - dt/dn*f_vort3d[i];
	    }
	    newst->pres += - dt/dn*f_pres;
	}
	if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) print_general_vector("Vorticity: ",newst->vort3d,dim,"\n");
	    (void) printf("Pressure: %f\n",newst->pres);
	}
}	/* end iF_flowThroughBoundaryState3d */

static void iF_flowThroughBoundaryState2d(
        double          *p0,
        HYPER_SURF      *hs,
        Front           *front,
        POINTER         params,
        POINTER         state)
{
	Tan_stencil **tsten;
	Nor_stencil *nsten;
	FLOW_THROUGH_PARAMS *ft_params = (FLOW_THROUGH_PARAMS*)params;
	POINT *oldp = ft_params->oldp;
	COMPONENT comp = ft_params->comp;
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double vort[3];		/* vorticity stencil */
	double pres[3];		/* pressure stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_vort;		/* vort flux */
	double f_pres;		/* pressure flux */
	double dn,dt = front->dt;
	STATE *oldst,*newst = (STATE*)state;
	STATE  **sts;
	POINTER sl,sr;
	int i,j,dim = front->rect_grid->dim;
	int nrad = 2;

	if (debugging("flow_through"))
	    printf("Entering iF_flowThroughBoundaryState2d()\n");

	FT_GetStatesAtPoint(oldp,oldp->hse,oldp->hs,&sl,&sr);
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	dn = FT_GridSizeInDir(nsten->nor,front);
	if (debugging("flow_through"))
	{
	    (void) printf("Normal grid size = %f\n",dn);
	    (void) print_Nor_stencil(front,nsten);
	}

	if (comp == negative_component(hs))
	    oldst = (STATE*)sl;
	else
	    oldst = (STATE*)sr;

	u[1] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    vort[j] = oldst->vort;
	    pres[j] = oldst->pres;
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			field->vel[i],getStateVel[i],&vtmp))
            {
        	vtmp = oldst->vel[i];
            }
	    u[1] += vtmp*dir[i];
	    newst->vel[i] = vtmp;
	}
	for (i = 0; i < dim; ++i)
	    v[2][i] -= u[2]*dir[i];
	if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->vort,
                            getStateVort,&vort[2]))
	    vort[2] = oldst->vort;
	if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],field->pres,
                            getStatePres,&pres[2]))
	    pres[2] = oldst->pres;

	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	newst->vort = oldst->vort - dt/dn*f_vort;
	newst->pres = oldst->pres - dt/dn*f_pres;

	tsten = FrontGetTanStencils(front,oldp,nrad);

	if (debugging("flow_through"))
	{
	    (void) printf("Ambient component: %d\n",comp);
	    (void) printf("Tangential grid size = %f\n",dn);
	    (void) print_Tan_stencil(front,tsten[0]);
	}

	if (comp == negative_component(hs))
	    sts = (STATE**)tsten[0]->leftst;
	else
	    sts = (STATE**)tsten[0]->rightst;

	for (i = 0; i < dim; ++i)
	    dir[i] = tsten[0]->dir[i];
	dn = FT_GridSizeInDir(dir,front);

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 3; ++j)
	{
	    vort[j] = sts[j-1]->vort;
	    pres[j] = sts[j-1]->pres;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += sts[j-1]->vel[i]*dir[i];
	    }
	    for (i = 0; i < dim; ++i)
	    {
		v[j][i] = sts[j-1]->vel[i] - dir[i]*u[j];
	    }
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;

	if (debugging("flow_through"))
	{
	    (void) printf("State after tangential sweep:\n");
	    (void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    (void) printf("Vorticity: %f\n",newst->vort);
	    (void) printf("Pressure: %f\n",newst->pres);
	}
}       /* end iF_flowThroughBoundaryState2d */

extern void ifluid_point_propagate(
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
	case MOVABLE_BODY_BOUNDARY:
	    return rgbody_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case NEUMANN_BOUNDARY:
	case GROWING_BODY_BOUNDARY:
	    return neumann_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	case DIRICHLET_BOUNDARY:
	    return dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	default:
	    return contact_point_propagate(front,wave,oldp,newp,oldhse,
					oldhs,dt,V);
	}
}       /* ifluid_point_propagate */

extern void ifluid_point_propagate_vd(
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
        case MOVABLE_BODY_BOUNDARY:
            return rgbody_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
        case NEUMANN_BOUNDARY:
        case GROWING_BODY_BOUNDARY:
            return neumann_point_propagate_vd(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
        case DIRICHLET_BOUNDARY:
            return dirichlet_point_propagate(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
        default:
            return contact_point_propagate_vd(front,wave,oldp,newp,oldhse,
                                        oldhs,dt,V);
        }
}       /* ifluid_point_propagate_vd */


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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        int i, dim = front->rect_grid->dim;
	double *m_pre = field->pres;
	double *m_vor = field->vort;
	double **m_vor3d = field->vort3d;
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	double dn,*h = front->rect_grid->h;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (ifluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}
	FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;

	for (i = 0; i < dim; ++i)
	{
            Coords(newp)[i] = Coords(oldp)[i];
	    newst->vel[i] = 0.0;
            FT_RecordMaxFrontSpeed(i,0.0,NULL,Coords(newp),front);
	}
	if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres))
            newst->pres = oldst->pres;
	if (dim == 2)
	{
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort))
		newst->vort = oldst->vort;
	}
	else
	{
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[0],
			getStateXvort,&newst->vort3d[0]))
		newst->vort3d[0] = oldst->vort3d[0];
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[1],
			getStateYvort,&newst->vort3d[1]))
		newst->vort3d[1] = oldst->vort3d[1];
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[2],
			getStateZvort,&newst->vort3d[2]))
		newst->vort3d[2] = oldst->vort3d[2];
	}
	FT_RecordMaxFrontSpeed(dim,0.0,NULL,Coords(newp),front);
        return;
}	/* end neumann_point_propagate */

static  void neumann_point_propagate_vd(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
//       IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
//        IF_FIELD *field = iFparams->field;
        int i, dim = front->rect_grid->dim;
//        double *m_pre = field->pres;
//        double *m_vor = field->vort;
//        double **m_vor3d = field->vort3d;
        STATE *oldst,*newst;
        POINTER sl,sr;
        COMPONENT comp;
        double nor[MAXD],p1[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;

        FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
        if (ifluid_comp(negative_component(oldhs)))
        {
            comp = negative_component(oldhs);
            oldst = (STATE*)sl;
            newst = (STATE*)left_state(newp);
        }
        else if (ifluid_comp(positive_component(oldhs)))
        {
            comp = positive_component(oldhs);
            oldst = (STATE*)sr;
            newst = (STATE*)right_state(newp);
        }
        FT_NormalAtPoint(oldp,front,nor,comp);

        dn = grid_size_in_direction(nor,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor[i]*dn;

        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i];
            newst->vel[i] = 0.0;
            FT_RecordMaxFrontSpeed(i,0.0,NULL,Coords(newp),front);
        }
/*
        if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
                        getStatePres,&newst->pres))
            newst->pres = oldst->pres;
        if (dim == 2)
        {
            if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
                        getStateVort,&newst->vort))
                newst->vort = oldst->vort;
        }
        else
        {
            if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[0],
                        getStateXvort,&newst->vort3d[0]))
                newst->vort3d[0] = oldst->vort3d[0];
            if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[1],
                        getStateYvort,&newst->vort3d[1]))
                newst->vort3d[1] = oldst->vort3d[1];
            if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[2],
                        getStateZvort,&newst->vort3d[2]))
                newst->vort3d[2] = oldst->vort3d[2];
        }
*/
        FT_RecordMaxFrontSpeed(dim,0.0,NULL,Coords(newp),front);
        return;
}       /* end neumann_point_propagate_vd */

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
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        double speed;
        int i, dim = front->rect_grid->dim;
	STATE *newst = NULL;
	STATE *bstate;
	COMPONENT comp;

	if (debugging("dirichlet_bdry"))
	{
	    (void) printf("Entering dirichlet_point_propagate()\n");
	    (void) print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
	}

	if (ifluid_comp(negative_component(oldhs)))
	{
	    newst = (STATE*)left_state(newp);
	    comp = negative_component(oldhs);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    newst = (STATE*)right_state(newp);
	    comp = positive_component(oldhs);
	}
	if (newst == NULL) return;	// node point

	if (boundary_state(oldhs) != NULL)
	{
	    bstate = (STATE*)boundary_state(oldhs);
            for (i = 0; i < dim; ++i)
	    {
	    	newst->vel[i] = bstate->vel[i];
	    	newst->vort3d[i] = 0.0;
		FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
					front);
	    }
	    speed = mag_vector(newst->vel,dim);
	    FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
            newst->pres = bstate->pres;
            newst->vort = 0.0;

	    if (debugging("dirichlet_bdry"))
	    {
		(void) printf("Preset boundary state:\n");
		(void) print_general_vector("Velocity: ",newst->vel,dim,"\n");
		(void) printf("Pressure: %f\n",newst->pres);
		(void) printf("Vorticity: %f\n",newst->vort);
	    }
	}
	else if (boundary_state_function(oldhs))
	{
	    if (strcmp(boundary_state_function_name(oldhs),
		       "iF_flowThroughBoundaryState") == 0)
	    {
		FLOW_THROUGH_PARAMS ft_params;
		oldp->hse = oldhse;
		oldp->hs = oldhs;
	    	ft_params.oldp = oldp;
	    	ft_params.comp = comp;
	    	(*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)&ft_params,(POINTER)newst);
	    }
	    else if (strcmp(boundary_state_function_name(oldhs),
		       "iF_timeDependBoundaryState") == 0)
	    {
		TIME_DEPENDENT_PARAMS *td_params = (TIME_DEPENDENT_PARAMS*)
				boundary_state_function_params(oldhs);
	    	(*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)td_params,(POINTER)newst);
	    }
            for (i = 0; i < dim; ++i)
		FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
					front);
	    speed = mag_vector(newst->vel,dim);
	    FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
	}
	if (debugging("dirichlet_bdry"))
	    (void) printf("Leaving dirichlet_point_propagate()\n");
        return;
}	/* end dirichlet_point_propagate */

static  void contact_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        double vel[MAXD],vt[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double *m_pre = field->pres;
	double *m_vor = field->vort;
	double **m_vel = field->vel;
	double **m_vor3d = field->vort3d;
	double vort3d[3];
	double *p0;
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;
	double pres,vort;

        (*front->vfunc)(front->vparams,front,oldp,oldhse,oldhs,vel);
	if (front->coordinate == 'R' || front->coordinate == 'r')
	{
	    for (i = 0; i < dim; ++i)
	    {
		Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
		set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newp),front);
	    }
	}
	else if (front->coordinate == 'C' || front->coordinate == 'c')
	{
            double coord_x,coord_y,coord_z;
            double vel_x,vel_y,vel_z;
            double new_coord_x,new_coord_y,new_coord_z;

            coord_x = Coords(oldp)[2]*cos(Coords(oldp)[0]);
            coord_y = Coords(oldp)[2]*sin(Coords(oldp)[0]);
            coord_z = Coords(oldp)[1];
            vel_x = -vel[0]*Coords(oldp)[2]*sin(Coords(oldp)[0]) + vel[2]*cos(Coords(oldp)[0]);
            vel_y = vel[0]*Coords(oldp)[2]*cos(Coords(oldp)[0]) + vel[2]*sin(Coords(oldp)[0]);
            vel_z = vel[1];

            new_coord_x = coord_x + dt*vel_x;
            new_coord_y = coord_y + dt*vel_y;
            new_coord_z = coord_z + dt*vel_z;

            double r = sqrt(sqr(new_coord_x)+sqr(new_coord_y));
            for (i = 0; i < dim; ++i)
            {
                if(i == 0)
                    Coords(newp)[i] = Coords(oldp)[i] + (dt*vel[i]/Coords(oldp)[2]);
                else
                    Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
                set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newp),front);
            }
	}

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	oldst = (STATE*)sl;
	p0 = Coords(newp);
	if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_pre,getStatePres,&pres))
	    pres = oldst->pres;
	if (dim == 2)
	{
	    if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor,getStateVort,&vort))
	    	vort = oldst->vort;
	}
	else if (dim == 3)
	{
            if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor3d[0],
                        getStateXvort,&vort3d[0]))
	    	vort3d[0] = oldst->vort3d[0];
            if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor3d[1],
                        getStateYvort,&vort3d[1]))
	    	vort3d[1] = oldst->vort3d[1];
            if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor3d[2],
                        getStateZvort,&vort3d[2]))
	    	vort3d[2] = oldst->vort3d[2];
	}

	newst = (STATE*)left_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	    newst->vort3d[i] = vort3d[i];
	}
	newst = (STATE*)right_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	{
	    newst->vel[i] = vel[i];
	    newst->vort3d[i] = vort3d[i];
	}

	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
}	/* end contact_point_propagate */

static  void contact_point_propagate_vd(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        IF_FIELD *field = iFparams->field;
        double vel[MAXD],vt[MAXD],s;
        int i, dim = front->rect_grid->dim;
//        double *m_pre = field->pres;
//        double *m_vor = field->vort;
//        double **m_vel = field->vel;
//        double **m_vor3d = field->vort3d;
//        double vort3d[3];

//        double *m_den = field->dens;
//        double dens;

        double *p0;
        STATE *oldst,*newst;
        POINTER sl,sr;
        COMPONENT comp;
//        double pres,vort;

        //printf("Enter contact_point_propagate_vd\n");
        //printf("oldp = [%lf, %lf, %lf]\n", Coords(oldp)[0], Coords(oldp)[1], Coords(oldp)[2]);

        //l_cartesian_vel_MAC_vd called
        (*front->vfunc)(front->vparams,front,oldp,oldhse,oldhs,vel);
        //printf("Leave front_>vfunc\n");

        //For Hao's simulation only.    Dan
        /*
        double L[3], U[3];
        double tol = 1e-4;
        for (i = 0; i < dim; ++i)
        {
            L[i] = front->rect_grid->L[i];
            U[i] = front->rect_grid->U[i];
            if (fabs(Coords(oldp)[i]-L[i]) < tol ||
                fabs(Coords(oldp)[i]-U[i]) < tol)
                if (fabs(vel[i]) > 1e-12)
                    vel[i] = 0.0;
        }
        vel[1] = 0;
        */

        if (front->coordinate == 'R' || front->coordinate == 'r')
        {
            for (i = 0; i < dim; ++i)
            {
                Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
                set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newp),front);
            }
        }
        else if (front->coordinate == 'C' || front->coordinate == 'c')
        {
            double coord_x,coord_y,coord_z;
            double vel_x,vel_y,vel_z;
            double new_coord_x,new_coord_y,new_coord_z;

            coord_x = Coords(oldp)[2]*cos(Coords(oldp)[0]);
            coord_y = Coords(oldp)[2]*sin(Coords(oldp)[0]);
            coord_z = Coords(oldp)[1];
            vel_x = -vel[0]*Coords(oldp)[2]*sin(Coords(oldp)[0]) + vel[2]*cos(Coords(oldp)[0]);
            vel_y = vel[0]*Coords(oldp)[2]*cos(Coords(oldp)[0]) + vel[2]*sin(Coords(oldp)[0]);
            vel_z = vel[1];

            new_coord_x = coord_x + dt*vel_x;
            new_coord_y = coord_y + dt*vel_y;
            new_coord_z = coord_z + dt*vel_z;

            double r = sqrt(sqr(new_coord_x)+sqr(new_coord_y));
            for (i = 0; i < dim; ++i)
            {
                if(i == 0)
                    Coords(newp)[i] = Coords(oldp)[i] + (dt*vel[i]/Coords(oldp)[2]);
                else
                    Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
                set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newp),front);
            }
        }

/*
        if (fabs(vel[0]) + fabs(vel[1]) + fabs(vel[2]) == 0) {
            printf("oldp = [%lf, %lf, %lf]\n", Coords(oldp)[0], Coords(oldp)[1], Coords(oldp)[2]);
            printf("newp = [%lf, %lf, %lf]\n", Coords(newp)[0], Coords(newp)[1], Coords(newp)[2]);
            printf("vel = [%lf, %lf, %lf]\n", vel[0], vel[1], vel[2]);
        }
*/

        FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);


        //printf("Leave FT_GetStatesAtPoint\n");


        oldst = (STATE*)sl;
        p0 = Coords(newp);

/*
        if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_pre,getStatePres,&pres))
            pres = oldst->pres;
        if (dim == 2)
        {
            if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor,getStateVort,&vort))
                vort = oldst->vort;
        }
        else if (dim == 3)
        {
            if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor3d[0],
                        getStateXvort,&vort3d[0]))
                vort3d[0] = oldst->vort3d[0];
            if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor3d[1],
                        getStateYvort,&vort3d[1]))
                vort3d[1] = oldst->vort3d[1];
            if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor3d[2],
                        getStateZvort,&vort3d[2]))
                vort3d[2] = oldst->vort3d[2];
        }
*/

//        if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_den,getStateDens,&dens))
//            dens = oldst->dens;

        newst = (STATE*)left_state(newp);
//        newst->vort = vort;
//        newst->pres = pres;
//        newst->dens = dens;

        for (i = 0; i < dim; ++i)
        {
            newst->vel[i] = vel[i];
//            newst->vort3d[i] = vort3d[i];
        }

        newst = (STATE*)right_state(newp);
//        newst->vort = vort;
//        newst->pres = pres;
//        newst->dens = dens;

        for (i = 0; i < dim; ++i)
        {
            newst->vel[i] = vel[i];
//            newst->vort3d[i] = vort3d[i];
        }

        s = mag_vector(vel,dim);
        FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);

        //printf("Leave contact_point_propagate_vd\n");

}       /* end contact_point_propagate_vd */

static  void rgbody_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	IF_FIELD *field = iFparams->field;
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pre = field->pres;
	double *m_vor = field->vort;
	double **m_vor3d = field->vort3d;
	double nor[MAXD],p1[MAXD];
	double *p0 = Coords(oldp);
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	if (ifluid_comp(negative_component(oldhs)))
	{
	    comp = negative_component(oldhs);
	    oldst = (STATE*)sl;
	    newst = (STATE*)left_state(newp);
	}
	else if (ifluid_comp(positive_component(oldhs)))
	{
	    comp = positive_component(oldhs);
	    oldst = (STATE*)sr;
	    newst = (STATE*)right_state(newp);
	}
	FT_NormalAtPoint(oldp,front,nor,comp);

	dn = grid_size_in_direction(nor,h,dim);
	for (i = 0; i < dim; ++i)
	    p1[i] = p0[i] + nor[i]*dn;

	fourth_order_point_propagate(front,NULL,oldp,newp,oldhse,
				oldhs,dt,vel);
	for (i = 0; i < dim; ++i) newst->vel[i] = vel[i];
	if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_pre,
			getStatePres,&newst->pres))
            newst->pres = oldst->pres;
	if (dim == 2)
	{
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor,
			getStateVort,&newst->vort))
		newst->vort = oldst->vort;
	}
	else
	{
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[0],
			getStateXvort,&newst->vort3d[0]))
		newst->vort3d[0] = oldst->vort3d[0];
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[1],
			getStateYvort,&newst->vort3d[1]))
		newst->vort3d[1] = oldst->vort3d[1];
	    if (!FT_IntrpStateVarAtCoords(front,comp,p1,m_vor3d[2],
			getStateZvort,&newst->vort3d[2]))
		newst->vort3d[2] = oldst->vort3d[2];
	}
        return;
}	/* end rgbody_point_propagate */

extern void fluid_print_front_states(
	FILE *outfile,
	Front *front,
        bool binary)
{
        char str[100];
        double val[3];

	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int dim = intfc->dim;

        if (binary == YES)
        {
            if (hardware_is_little_endian())
            {
                sprintf(str,"Interface ifluid states:\n");
                fwrite(str, sizeof(char), 25, outfile);
                next_point(intfc,NULL,NULL,NULL);
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    val[0] = endian_double_swap(getStatePres(sl));
                    val[1] = endian_double_swap(getStatePres(sr));
                    fwrite(val, sizeof(double), 2, outfile);

                    if (dim == 2)
                    {
                        val[0] = endian_double_swap(getStateXvel(sl));
                        val[1] = endian_double_swap(getStateXvel(sr));
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = endian_double_swap(getStateYvel(sl));
                        val[1] = endian_double_swap(getStateYvel(sr));
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = endian_double_swap(getStateVort(sl));
                        val[1] = endian_double_swap(getStateVort(sr));
                        fwrite(val, sizeof(double), 2 , outfile);
                    }
                    else if (dim == 3)
                    {
                        val[0] = endian_double_swap(getStateXvel(sl));
                        val[1] = endian_double_swap(getStateXvel(sr));
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = endian_double_swap(getStateYvel(sl));
                        val[1] = endian_double_swap(getStateYvel(sr));
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = endian_double_swap(getStateZvel(sl));
                        val[1] = endian_double_swap(getStateZvel(sr));
                        fwrite(val, sizeof(double), 2, outfile);
                    }
                }
            }
            else //big endian system
            {
                sprintf(str,"Interface ifluid states:\n");
                fwrite(str, sizeof(char), 25, outfile);
                next_point(intfc,NULL,NULL,NULL);
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    val[0] = getStatePres(sl);
                    val[1] = getStatePres(sr);
                    fwrite(val, sizeof(double), 2, outfile);

                    if (dim == 2)
                    {
                        val[0] = getStateXvel(sl);
                        val[1] = getStateXvel(sr);
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = getStateYvel(sl);
                        val[1] = getStateYvel(sr);
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = getStateVort(sl);
                        val[1] = getStateVort(sr);
                        fwrite(val, sizeof(double), 2 , outfile);
                    }
                    else if (dim == 3)
                    {
                        val[0] = getStateXvel(sl);
                        val[1] = getStateXvel(sr);
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = getStateYvel(sl);
                        val[1] = getStateYvel(sr);
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = getStateZvel(sl);
                        val[1] = getStateZvel(sr);
                        fwrite(val, sizeof(double), 2, outfile);
                    }
                }
            }
        }
        else //ASCII output
        {
            fprintf(outfile,"Interface ifluid states:\n");
            next_point(intfc,NULL,NULL,NULL);
            while (next_point(intfc,&p,&hse,&hs))
            {
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                fprintf(outfile,"%24.18g %24.18g\n",getStatePres(sl),
                        getStatePres(sr));
                if (dim == 2)
                {
                    fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                            getStateXvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                            getStateYvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateVort(sl),
                            getStateVort(sr));
                }
                if (dim == 3)
                {
                    fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                            getStateXvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                            getStateYvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateZvel(sl),
                            getStateZvel(sr));
                }
            }
        }

}       /* end fluid_print_front_states */


extern void fluid_print_front_states_vd(
        FILE *outfile,
        Front *front,
        bool binary)
{
        char str[100];
        double val[3];

        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        int dim = intfc->dim;

        if (binary == YES)
        {
            if (hardware_is_little_endian())
            {
                sprintf(str,"Interface ifluid states:\n");
                fwrite(str, sizeof(char), 25, outfile);
                next_point(intfc,NULL,NULL,NULL);
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    val[0] = endian_double_swap(getStatePres(sl));
                    val[1] = endian_double_swap(getStatePres(sr));
                    fwrite(val, sizeof(double), 2, outfile);

                    if (dim == 2)
                    {
                        val[0] = endian_double_swap(getStateXvel(sl));
                        val[1] = endian_double_swap(getStateXvel(sr));
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = endian_double_swap(getStateYvel(sl));
                        val[1] = endian_double_swap(getStateYvel(sr));
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = endian_double_swap(getStateVort(sl));
                        val[1] = endian_double_swap(getStateVort(sr));
                        fwrite(val, sizeof(double), 2 , outfile);
                    }
                    else if (dim == 3)
                    {
                        val[0] = endian_double_swap(getStateXvel(sl));
                        val[1] = endian_double_swap(getStateXvel(sr));
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = endian_double_swap(getStateYvel(sl));
                        val[1] = endian_double_swap(getStateYvel(sr));
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = endian_double_swap(getStateZvel(sl));
                        val[1] = endian_double_swap(getStateZvel(sr));
                        fwrite(val, sizeof(double), 2, outfile);
                    }
                    // for vd
                    val[0] = endian_double_swap(getStateDens(sl));
                    val[1] = endian_double_swap(getStateDens(sr));
                    fwrite(val, sizeof(double), 2, outfile);
                    val[0] = endian_double_swap(getStateDensOld(sl));
                    val[1] = endian_double_swap(getStateDensOld(sr));
                    fwrite(val, sizeof(double), 2, outfile);
                    val[0] = endian_double_swap(getStateConc(sl));
                    val[1] = endian_double_swap(getStateConc(sr));
                    fwrite(val, sizeof(double), 2, outfile);
                }
            }
            else //big endian system
            {
                sprintf(str,"Interface ifluid states:\n");
                fwrite(str, sizeof(char), 25, outfile);
                next_point(intfc,NULL,NULL,NULL);
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    val[0] = getStatePres(sl);
                    val[1] = getStatePres(sr);
                    fwrite(val, sizeof(double), 2, outfile);

                    if (dim == 2)
                    {
                        val[0] = getStateXvel(sl);
                        val[1] = getStateXvel(sr);
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = getStateYvel(sl);
                        val[1] = getStateYvel(sr);
                        fwrite(val, sizeof(double), 2 , outfile);
                        val[0] = getStateVort(sl);
                        val[1] = getStateVort(sr);
                        fwrite(val, sizeof(double), 2 , outfile);
                    }
                    else if (dim == 3)
                    {
                        val[0] = getStateXvel(sl);
                        val[1] = getStateXvel(sr);
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = getStateYvel(sl);
                        val[1] = getStateYvel(sr);
                        fwrite(val, sizeof(double), 2, outfile);
                        val[0] = getStateZvel(sl);
                        val[1] = getStateZvel(sr);
                        fwrite(val, sizeof(double), 2, outfile);
                    }
                    // for vd
                    val[0] = getStateDens(sl);
                    val[1] = getStateDens(sr);
                    fwrite(val, sizeof(double), 2, outfile);
                    val[0] = getStateDensOld(sl);
                    val[1] = getStateDensOld(sr);
                    fwrite(val, sizeof(double), 2, outfile);
                    val[0] = getStateConc(sl);
                    val[1] = getStateConc(sr);
                    fwrite(val, sizeof(double), 2, outfile);
                }
            }
        }
        else //ASCII output
        {
            fprintf(outfile,"Interface ifluid states:\n");
            next_point(intfc,NULL,NULL,NULL);
            while (next_point(intfc,&p,&hse,&hs))
            {
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                fprintf(outfile,"%24.18g %24.18g\n",getStatePres(sl),
                        getStatePres(sr));
                if (dim == 2)
                {
                    fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                            getStateXvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                            getStateYvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateVort(sl),
                            getStateVort(sr));
                }
                if (dim == 3)
                {
                    fprintf(outfile,"%24.18g %24.18g\n",getStateXvel(sl),
                            getStateXvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateYvel(sl),
                            getStateYvel(sr));
                    fprintf(outfile,"%24.18g %24.18g\n",getStateZvel(sl),
                            getStateZvel(sr));
                }
                //for vd
                fprintf(outfile,"%24.18g %24.18g\n",getStateDens(sl),
                        getStateDens(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateDensOld(sl),
                        getStateDensOld(sr));
                fprintf(outfile,"%24.18g %24.18g\n",getStateConc(sl),
                        getStateConc(sr));
            }
        }

} /* end fluid_print_front_states_vd */


extern void fluid_read_front_states(
	FILE *infile,
	Front *front,
        bool binary)
{
        double val[3];

	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *lstate,*rstate;
	int dim = intfc->dim;

        if (binary == YES)
        {
            next_output_line_containing_string(infile,"Interface ifluid states:");
            next_point(intfc,NULL,NULL,NULL);
            if (hardware_is_little_endian())
            {
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    lstate = (STATE*)sl;        rstate = (STATE*)sr;
                    fread(val, sizeof(double), 2, infile);
                    lstate->pres = endian_double_swap(val[0]);
                    rstate->pres = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[0] = endian_double_swap(val[0]);
                    rstate->vel[0] = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[1] = endian_double_swap(val[0]);
                    rstate->vel[1] = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    if (dim == 2)
                    {
                        lstate->vort = endian_double_swap(val[0]);
                        rstate->vort = endian_double_swap(val[1]);
                    }
                    else if (dim == 3)
                    {
                        lstate->vel[2] = endian_double_swap(val[0]);
                        rstate->vel[2] = endian_double_swap(val[1]);
                    }
                }
            }
            else //big endian system
            {
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    lstate = (STATE*)sl;        rstate = (STATE*)sr;
                    fread(val, sizeof(double), 2, infile);
                    lstate->pres = val[0];
                    rstate->pres = val[1];
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[0] = val[0];
                    rstate->vel[0] = val[1];
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[1] = val[0];
                    rstate->vel[1] = val[1];
                    fread(val, sizeof(double), 2, infile);
                    if (dim == 2)
                    {
                        lstate->vort = val[0];
                        rstate->vort = val[1];
                    }
                    else if (dim == 3)
                    {
                        lstate->vel[2] = val[0];
                        rstate->vel[2] = val[1];
                    }
                }
            }
        }
        else //read ASCII restart file
        {
            next_output_line_containing_string(infile,"Interface ifluid states:");
            next_point(intfc,NULL,NULL,NULL);
            while (next_point(intfc,&p,&hse,&hs))
            {
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                lstate = (STATE*)sl;        rstate = (STATE*)sr;
                fscanf(infile,"%lf %lf",&lstate->pres,&rstate->pres);
                fscanf(infile,"%lf %lf",&lstate->vel[0],&rstate->vel[0]);
                fscanf(infile,"%lf %lf",&lstate->vel[1],&rstate->vel[1]);
                if (dim == 2)
                    fscanf(infile,"%lf %lf",&lstate->vort,&rstate->vort);
                if (dim == 3)
                    fscanf(infile,"%lf %lf",&lstate->vel[2],&rstate->vel[2]);
            }
        }
}	/* end fluid_read_front_states */

extern void fluid_read_front_states_vd(
        FILE *infile,
        Front *front,
        bool binary)
{
        double val[3];

        INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *lstate,*rstate;
        int dim = intfc->dim;

        if (binary == YES)
        {
            next_output_line_containing_string(infile,"Interface ifluid states:");
            next_point(intfc,NULL,NULL,NULL);
            if (hardware_is_little_endian())
            {
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    lstate = (STATE*)sl;        rstate = (STATE*)sr;
                    fread(val, sizeof(double), 2, infile);
                    lstate->pres = endian_double_swap(val[0]);
                    rstate->pres = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[0] = endian_double_swap(val[0]);
                    rstate->vel[0] = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[1] = endian_double_swap(val[0]);
                    rstate->vel[1] = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    if (dim == 2)
                    {
                        lstate->vort = endian_double_swap(val[0]);
                        rstate->vort = endian_double_swap(val[1]);
                    }
                    else if (dim == 3)
                    {
                        lstate->vel[2] = endian_double_swap(val[0]);
                        rstate->vel[2] = endian_double_swap(val[1]);
                    }
                    //for vd
                    fread(val, sizeof(double), 2, infile);
                    lstate->dens = endian_double_swap(val[0]);
                    rstate->dens = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    lstate->dens_old = endian_double_swap(val[0]);
                    rstate->dens_old = endian_double_swap(val[1]);
                    fread(val, sizeof(double), 2, infile);
                    lstate->conc = endian_double_swap(val[0]);
                    rstate->conc = endian_double_swap(val[1]);
                }
            }
            else //big endian system
            {
                while (next_point(intfc,&p,&hse,&hs))
                {
                    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                    lstate = (STATE*)sl;        rstate = (STATE*)sr;
                    fread(val, sizeof(double), 2, infile);
                    lstate->pres = val[0];
                    rstate->pres = val[1];
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[0] = val[0];
                    rstate->vel[0] = val[1];
                    fread(val, sizeof(double), 2, infile);
                    lstate->vel[1] = val[0];
                    rstate->vel[1] = val[1];
                    fread(val, sizeof(double), 2, infile);
                    if (dim == 2)
                    {
                        lstate->vort = val[0];
                        rstate->vort = val[1];
                    }
                    else if (dim == 3)
                    {
                        lstate->vel[2] = val[0];
                        rstate->vel[2] = val[1];
                    }
                    // for vd
                    fread(val, sizeof(double), 2, infile);
                    lstate->dens = val[0];
                    rstate->dens = val[1];
                    fread(val, sizeof(double), 2, infile);
                    lstate->dens_old = val[0];
                    rstate->dens_old = val[1];
                    fread(val, sizeof(double), 2, infile);
                    lstate->conc = val[0];
                    rstate->conc = val[1];
                }
            }
        }
        else //read ASCII restart file
        {
            next_output_line_containing_string(infile,"Interface ifluid states:");
            next_point(intfc,NULL,NULL,NULL);
            while (next_point(intfc,&p,&hse,&hs))
            {
                FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
                lstate = (STATE*)sl;        rstate = (STATE*)sr;
                fscanf(infile,"%lf %lf",&lstate->pres,&rstate->pres);
                fscanf(infile,"%lf %lf",&lstate->vel[0],&rstate->vel[0]);
                fscanf(infile,"%lf %lf",&lstate->vel[1],&rstate->vel[1]);
                if (dim == 2)
                    fscanf(infile,"%lf %lf",&lstate->vort,&rstate->vort);
                if (dim == 3)
                    fscanf(infile,"%lf %lf",&lstate->vel[2],&rstate->vel[2]);
                // for vd
                fscanf(infile,"%lf %lf",&lstate->dens,&rstate->dens);
                fscanf(infile,"%lf %lf",&lstate->dens_old,&rstate->dens_old);
                fscanf(infile,"%lf %lf",&lstate->conc,&rstate->conc);
            }
        }
} /* end fluid_read_front_states_vd */

extern void read_iFparams(
	char *inname,
	IF_PARAMS *iFparams)
{
	char string[100];
	int i,dim = iFparams->dim;
	FILE *infile = fopen(inname,"r");

	iFparams->num_scheme = ERROR_SCHEME;
	CursorAfterString(infile,"Enter projection type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'S':
	case 's':
	    iFparams->num_scheme = SIMPLE;
	    break;
	case 'B':
	case 'b':
	    iFparams->num_scheme = BELL_COLELLA;
	    break;
	case 'K':
	case 'k':
	    iFparams->num_scheme = KIM_MOIN;
	    break;
	case 'P':
	case 'p':
	    iFparams->num_scheme = PEROT_BOTELLA;
	}
	assert(iFparams->num_scheme != ERROR_SCHEME);

	fclose(infile);
}	/* end read_iFparams */

extern void read_iF_movie_options(
	char *inname,
	IF_PARAMS *iFparams)
{
	static IF_MOVIE_OPTION *movie_option;
	FILE *infile = fopen(inname,"r");
	char string[100];

	FT_ScalarMemoryAlloc((POINTER*)&movie_option,sizeof(IF_MOVIE_OPTION));
	iFparams->movie_option = movie_option;
	CursorAfterString(infile,"Type y to make movie of pressure:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	    movie_option->plot_pres = YES;
	CursorAfterString(infile,"Type y to make movie of vorticity:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	    movie_option->plot_vort = YES;
	CursorAfterString(infile,"Type y to make movie of velocity:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	    movie_option->plot_velo = YES;
        CursorAfterString(infile,"Type y to make movie of density:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
            movie_option->plot_dens = YES;
        CursorAfterString(infile,"Type y to make movie of concentration:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
            movie_option->plot_conc = YES;
        CursorAfterString(infile,"Type y to make movie of component:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'Y' || string[0] == 'y')
            movie_option->plot_comp = YES;

	if (iFparams->dim == 3)
	{
	    CursorAfterString(infile,"Type y to make yz cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[0] = YES;
	    CursorAfterString(infile,"Type y to make xz cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[1] = YES;
	    CursorAfterString(infile,"Type y to make xy cross section movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_cross_section[2] = YES;
	    CursorAfterString(infile,
			"Type y to make cross sectional component movie:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'Y' || string[0] == 'y')
		movie_option->plot_comp = YES;
	}
	CursorAfterString(infile,"Type y to make vector velocity field movie:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	if (string[0] == 'Y' || string[0] == 'y')
	    movie_option->plot_vel_vector = YES;
	fclose(infile);
}	/* end read_iF_movie_options */

extern boolean isDirichletPresetBdry(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp)
{
	HYPER_SURF *hs;
	POINTER intfc_state;
	double crx_coords[MAXD];

	if (!FT_StateStructAtGridCrossing(front,icoords,dir,
                                comp,&intfc_state,&hs,crx_coords))
	    return NO;
	if (wave_type(hs) != DIRICHLET_BOUNDARY)
	    return NO;
	if (boundary_state(hs) != NULL)
	    return NO;
	return YES;
}	/* end isDirichletPresetBdry */

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
	switch (string[0])
	{
	case 'C':
	case 'c':
	    td_params->td_type = CONSTANT;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    break;
	case 'P':
	case 'p':
	    td_params->td_type = PULSE_FUNC;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    CursorAfterString(infile,"Enter peak velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_peak[i]);
		(void) printf("%f ",td_params->v_peak[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter peak pressure:");
	    fscanf(infile,"%lf ",&td_params->p_peak);
	    (void) printf("%f\n",td_params->p_peak);
	    CursorAfterString(infile,"Enter tail velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_tail[i]);
		(void) printf("%f ",td_params->v_tail[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter tail pressure:");
	    fscanf(infile,"%lf ",&td_params->p_tail);
	    (void) printf("%f\n",td_params->p_tail);
	    CursorAfterString(infile,"Enter time to rise:");
	    fscanf(infile,"%lf ",&td_params->T[0]);
	    (void) printf("%f\n",td_params->T[0]);
	    CursorAfterString(infile,"Enter time to reach peak:");
	    fscanf(infile,"%lf ",&td_params->T[1]);
	    (void) printf("%f\n",td_params->T[1]);
	    CursorAfterString(infile,"Enter time to fall:");
	    fscanf(infile,"%lf ",&td_params->T[2]);
	    (void) printf("%f\n",td_params->T[2]);
	    CursorAfterString(infile,"Enter time to reach tail:");
	    fscanf(infile,"%lf ",&td_params->T[3]);
	    (void) printf("%f\n",td_params->T[3]);
	    break;
	case 'S':
	case 's':
	    td_params->td_type = SINE_FUNC;
	    CursorAfterString(infile,"Enter base velocity:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_base[i]);
		(void) printf("%f ",td_params->v_base[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter base pressure:");
	    fscanf(infile,"%lf ",&td_params->p_base);
	    (void) printf("%f\n",td_params->p_base);
	    CursorAfterString(infile,"Enter velocity amplitude:");
	    for (i = 0; i < dim; ++i)
	    {
		fscanf(infile,"%lf ",&td_params->v_amp[i]);
		(void) printf("%f ",td_params->v_amp[i]);
	    }
	    (void) printf("\n");
	    CursorAfterString(infile,"Enter pressure amplitude:");
	    fscanf(infile,"%lf ",&td_params->p_amp);
	    (void) printf("%f\n",td_params->p_amp);
	    CursorAfterString(infile,"Enter oscillation period:");
	    fscanf(infile,"%lf ",&td_params->omega);
	    (void) printf("%f\n",td_params->omega);
	    td_params->omega = 2.0*PI/td_params->omega;
	    CursorAfterString(infile,"Enter initial phase:");
	    fscanf(infile,"%lf ",&td_params->phase);
	    (void) printf("%f\n",td_params->phase);
	    td_params->phase *= PI/180.0;
	    break;
	default:
	    (void) printf("Unknown type of time-dependent function!\n");
	    clean_up(ERROR);
	}

	*params = (POINTER)td_params;
}	/* end get_time_dependent_params */

extern void recordBdryEnergyFlux(
	Front *front,
	char *out_name)
{
	RECT_GRID *rgr = front->rect_grid;
	INTERFACE *intfc = front->interf;
	HYPER_SURF *hs;
	double Energy_in,Energy_out;
	int dir,side;
	int dim = rgr->dim;
	static FILE *ein_file,*eout_file;
	char file_name[100];

	if (ein_file == NULL && pp_mynode() == 0)
	{
	    sprintf(file_name,"%s/in_energy.xg",out_name);
	    ein_file = fopen(file_name,"w");
	    fprintf(ein_file,"\"Energy influx\" vs. time\n");
	    sprintf(file_name,"%s/out_energy.xg",out_name);
	    eout_file = fopen(file_name,"w");
	    fprintf(eout_file,"\"Energy outflux\" vs. time\n");
	}
	Energy_in = Energy_out = 0.0;
	for (dir = 0; dir < dim; ++dir)
	for (side = 0; side < 2; ++side)
	{
	    hs = FT_RectBoundaryHypSurf(intfc,DIRICHLET_BOUNDARY,dir,side);
	    if (hs == NULL) continue;
	    if (boundary_state(Hyper_surf(hs)))
	    {
		addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,YES);
	    }
	    else if (boundary_state_function(Hyper_surf(hs)))
	    {
		if (strcmp(boundary_state_function_name(hs),
                       "iF_timeDependBoundaryState") == 0)
		{
		    addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,YES);
		}
		else if (strcmp(boundary_state_function_name(hs),
                       "iF_flowThroughBoundaryState") == 0)
		{
		    addToEnergyFlux(rgr,hs,&Energy_in,&Energy_out,dir,side,NO);
		}
	    }
	}
	pp_global_sum(&Energy_in,1);
	pp_global_sum(&Energy_out,1);

	if (pp_mynode() == 0)
	{
	    fprintf(ein_file,"%f %f\n",front->time,Energy_in);
	    fprintf(eout_file,"%f %f\n",front->time,Energy_out);
	}
}	/* end recordBdryEnergyFlux */


static void addToEnergyFlux(
	RECT_GRID *rgr,
	HYPER_SURF *hs,
	double *Energy_in,
	double *Energy_out,
	int dir,
	int side,
	boolean is_influx)
{
	int i,dim = rgr->dim;
	double *L = rgr->L;
	double *U = rgr->U;
	CURVE *c;
	SURFACE *s;
	BOND *b;
	TRI *t;
	double ave_coord,engy_flux,vel;
	boolean is_outside_hse,is_left_side;
	STATE *sl,*sr,*state;
	POINT *pts[MAXD];

	is_left_side = (ifluid_comp(negative_component(hs))) ? YES : NO;
	switch (dim)
	{
	case 2:
	    c = Curve_of_hs(hs);
	    for (b = c->first; b != NULL; b = b->next)
	    {
		is_outside_hse = NO;
		pts[0] = b->start;
		pts[1] = b->end;
		for (i = 1; i < dim; ++i)
		{
		    ave_coord = (Coords(pts[0])[(dir+i)%dim] +
				 Coords(pts[1])[(dir+i)%dim])/2.0;
		    if (ave_coord < L[(dir+i)%dim] ||
			ave_coord > U[(dir+i)%dim])
			is_outside_hse = YES;
		}
		if (is_outside_hse) continue;
		engy_flux = 0.0;
		for (i = 0; i < 2; ++i)
		{
		    FT_GetStatesAtPoint(pts[i],Hyper_surf_element(b),hs,
				(POINTER*)&sl,(POINTER*)&sr);
		    state = (is_left_side) ? sl : sr;
		    if (is_influx)
		    	vel = (side == 0) ? state->vel[dir] : -state->vel[dir];
		    else
		    	vel = (side == 1) ? state->vel[dir] : -state->vel[dir];
		    //engy_flux += 0.5*state->dens*(sqr(state->vel[0]) +
		    engy_flux += 0.5*(sqr(state->vel[0]) +
					sqr(state->vel[1]))*vel;
		}
		engy_flux *= bond_length(b)/2.0;
		if (is_influx)
		    *Energy_in += engy_flux;
		else
		    *Energy_out += engy_flux;
	    }
	    break;
	case 3:
	    s = Surface_of_hs(hs);
	    for (t = first_tri(s); !at_end_of_tri_list(t,s); t = t->next)
	    {
		is_outside_hse = NO;
		for (i = 0; i < 3; ++i)
		    pts[i] = Point_of_tri(t)[i];
		for (i = 1; i < dim; ++i)
		{
		    ave_coord = (Coords(pts[0])[(dir+i)%dim] +
				 Coords(pts[1])[(dir+i)%dim] +
				 Coords(pts[2])[(dir+i)%dim])/3.0;
		    if (ave_coord < L[(dir+i)%dim] ||
			ave_coord > U[(dir+i)%dim])
			is_outside_hse = YES;
		}
		if (is_outside_hse) continue;
		engy_flux = 0.0;
		for (i = 0; i < 3; ++i)
		{
		    FT_GetStatesAtPoint(pts[i],Hyper_surf_element(t),hs,
				(POINTER*)&sl,(POINTER*)&sr);
		    state = (is_left_side) ? sl : sr;
		    if (is_influx)
		    	vel = (side == 0) ? state->vel[dir] : -state->vel[dir];
		    else
		    	vel = (side == 1) ? state->vel[dir] : -state->vel[dir];
		    //engy_flux += 0.5*state->dens*(sqr(state->vel[0]) +
		    engy_flux += 0.5*(sqr(state->vel[0]) +
					sqr(state->vel[1]) +
					sqr(state->vel[2]))*vel;
		}
		engy_flux *= bond_length(b)/3.0;
		if (is_influx)
		    *Energy_in += engy_flux;
		else
		    *Energy_out += engy_flux;
	    }
	    break;
	}
}	/* end addToEnergyFlux */
