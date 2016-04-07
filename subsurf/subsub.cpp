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


#include <iFluid.h>
#include <crystal.h>
#include "subsurf.h"

	/*  Function Declarations */
static void neumann_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void dirichlet_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static void contact_point_propagate(Front*,POINTER,POINT*,POINT*,
                        HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static double zero_state(COMPONENT,double*,L_STATE&,int,IF_PARAMS*);

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
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double vel[MAXD],vt[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double dn,*h = front->rect_grid->h;
	double *m_pre = iFparams->field->pres;
	double *m_vor = iFparams->field->vort;
	double **m_vor3d = iFparams->field->vort3d;
	double nor[MAXD],tan[MAXD],p1[MAXD];
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
	tan[0] = -nor[1]; 	tan[1] = nor[0];

	if (wave_type(oldhs) == MOVABLE_BODY_BOUNDARY)
	{
            double omega_dt,crds_com[MAXD];
            omega_dt = angular_velo(oldhs)*dt;
            for (i = 0; i < dim; ++i)
            {
                vel[i] = center_of_mass_velo(oldhs)[i];
                crds_com[i] = Coords(oldp)[i] +dt*vel[i] - 
			center_of_mass(oldhs)[i];
            }
            vel[0] += -angular_velo(oldhs)*crds_com[1]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[0]*sin(omega_dt);
            vel[1] +=  angular_velo(oldhs)*crds_com[0]*cos(omega_dt) -
                     angular_velo(oldhs)*crds_com[1]*sin(omega_dt);
	}
	else
	{
            for (i = 0; i < dim; ++i)
	    	vel[i] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
            Coords(newp)[i] = Coords(newp)[i] + dt*vel[i];
	    newst->vel[i] = vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
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
	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
        return;
}	/* end neumann_point_propagate */

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
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
	double speed;
        int i, dim = front->rect_grid->dim;
	STATE *newst;
	STATE *bstate;
	FLOW_THROUGH_PARAMS ft_params;
	COMPONENT comp;

	if (debugging("dirichlet_bdry"))
	{
	    printf("Entering dirichlet_point_propagate()\n");
	    print_general_vector("oldp:  ",Coords(oldp),dim,"\n");
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
		printf("Preset boundary state:\n");
		print_general_vector("Velocity: ",newst->vel,dim,"\n");
		printf("Pressure: %f\n",newst->pres);
		printf("Vorticity: %f\n",newst->vort);
	    }
	}
	else if (boundary_state_function(oldhs))
	{
	    ft_params.oldp = oldp;
	    ft_params.comp = comp;
	    (*boundary_state_function(oldhs))(Coords(oldp),oldhs,front,
			(POINTER)&ft_params,(POINTER)newst);	
	    for (i = 0; i < dim; ++i)
                FT_RecordMaxFrontSpeed(i,fabs(newst->vel[i]),NULL,Coords(newp),
                                        front);
            speed = mag_vector(newst->vel,dim);
            FT_RecordMaxFrontSpeed(dim,speed,NULL,Coords(newp),front);
	}
	if (debugging("dirichlet_bdry"))
	    printf("Leaving dirichlet_point_propagate()\n");
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
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double vel[MAXD],vt[MAXD],s;
        int i, dim = front->rect_grid->dim;
	double *m_pre = iFparams->field->pres;
	double *m_vor = iFparams->field->vort;
	double **m_vel = iFparams->field->vel;
	double *p0;
	STATE *oldst,*newst;
	POINTER sl,sr;
	COMPONENT comp;
	double pres,vort;

        (*front->vfunc)(front->vparams,front,oldp,oldhse,oldhs,vel);
        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(newp)[i] + dt*vel[i];
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }

	FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);
	oldst = (STATE*)sl;
	p0 = Coords(newp);
	if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_pre,getStatePres,&pres))
	    pres = oldst->pres;
	if (!FT_IntrpStateVarAtCoords(front,-1,p0,m_vor,getStateVort,&vort))
	    vort = oldst->vort;

	newst = (STATE*)left_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	    newst->vel[i] = vel[i];
	newst = (STATE*)right_state(newp);
	newst->vort = vort;
	newst->pres = pres;
	for (i = 0; i < dim; ++i)
	    newst->vel[i] = vel[i];

	s = mag_vector(vel,dim);
	FT_RecordMaxFrontSpeed(dim,s,NULL,Coords(newp),front);
}	/* end contact_point_propagate */

void crystal_point_propagate(
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
        double s0,s1,s2,grad_s;
        STATE *sl,*sr,*state,*bstate;
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
        double *solute = cRparams->field->solute;
        double D = cRparams->D;
        double k = cRparams->k;
        double C_0 = cRparams->C_0;
        double C_eq = cRparams->C_eq;
        double rho_s = cRparams->rho_s;
	double kappa,c1,c2;
	POINT_PROP_SCHEME point_prop_scheme = cRparams->point_prop_scheme;

        if (wave_type(oldhs) != GROWING_BODY_BOUNDARY)
        {
            for (i = 0; i < dim; ++i)
                Coords(newp)[i] = Coords(oldp)[i];
	    if (wave_type(oldhs) == NEUMANN_BOUNDARY)
	    {
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
	    }
	    else if (wave_type(oldhs) == DIRICHLET_BOUNDARY)
	    {
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
	    }
	    else
	    {
            	FT_GetStatesAtPoint(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
            	state = (negative_component(oldhs) == SOLUTE_COMP) ?
                        	(STATE*)left_state(newp) :
                        	(STATE*)right_state(newp);
                state->solute = (negative_component(oldhs) == SOLUTE_COMP) ? 
				sl->solute : sr->solute;
		if (cRparams->max_solute < state->solute)
			cRparams->max_solute = state->solute;
		if (cRparams->min_solute > state->solute)
			cRparams->min_solute = state->solute;
	    }
            return;
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

	/* Update the state of the new interface point */
	state = (negative_component(oldhs) == SOLUTE_COMP) ? 
			(STATE*)left_state(newp) : (STATE*)right_state(newp);
	switch (point_prop_scheme)
	{
	case EXPLICIT_EULER:
	    s0 = s0 + 2.0*dt/dn*(D*(s1 - s0)/dn - k*(s0 - C_eq));
	    break;
	case IMPLICIT_EULER:
	    c1 = 0.5*D*dt*(1.0/dn + kappa)/dn;
	    c2 = 0.5*k*dt/dn;
	    s0 = (s0 + c1*(s1 - s0) - c2*(s0 - C_eq) + c1*s1 + c2*C_eq)
			/(1.0 + c1 + c2);
	    break;
	case MIDDLE_POINT:
	    c1 = D*dt*(2.0/dn + kappa)/dn;
	    c2 = 2.0*k*dt/dn;
	    s0 = (s0 + c1*s1 + c2*C_eq)
			/(1.0 + c1 + c2);
	    break;
	case CONSTANT_STATE:
	    s0 = C_eq;
	    break;
	}
        state->solute = s0;
	if (cRparams->max_solute < state->solute)
		cRparams->max_solute = state->solute;
	if (cRparams->min_solute > state->solute)
		cRparams->min_solute = state->solute;
        for (i = 0; i < dim; ++i)
        {
            vel[i] = nor[i]*k*(s0 - C_eq)/rho_s;
            FT_RecordMaxFrontSpeed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }
}       /* crystal_point_propagate */

void read_dirichlet_bdry_data(
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
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
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
		CursorAfterString(infile,"Enter type of Dirichlet boundary:");
		fscanf(infile,"%s",s);
		switch (s[0])
		{
		case 'c':			// Constant state
		case 'C':
		    CursorAfterString(infile,"Enter solute concentration:");
		    fscanf(infile,"%lf",&state.solute);
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
}	/* end read_dirichlet_bdry_data */

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
	FT_MakeGridIntfc(front);
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

void read_ss_dirichlet_bdry_data(
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
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf",&state.vel[k]);
			(void) printf("%f ",state.vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state.pres);
		    (void) printf("%f\n",state.pres);
		    CursorAfterString(infile,"Enter solute concentration:");
                    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_SetDirichletBoundary(front,NULL,NULL,NULL,
				(POINTER)&state,hs);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_SetDirichletBoundary(front,ss_flowThroughBoundaryState,
				"ss_flowThroughBoundaryState",NULL,NULL,hs);
		    break;
		}
	    }
            if (f_basic.boundary[i][1] == DIRICHLET_BOUNDARY)
	    {
		hs = NULL;
                if (rect_boundary_type(front->interf,i,1) == DIRICHLET_BOUNDARY)                    hs = FT_RectBoundaryHypSurf(front->interf,DIRICHLET_BOUNDARY,
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
		    CursorAfterString(infile,"Enter velocity:");
		    for (k = 0; k < dim; ++k)
		    {
			fscanf(infile,"%lf ",&state.vel[k]);
			(void) printf("%f ",state.vel[k]);
		    }
		    (void) printf("\n");
		    CursorAfterString(infile,"Enter pressure:");
		    fscanf(infile,"%lf",&state.pres);
		    (void) printf("%f\n",state.pres);
		    CursorAfterString(infile,"Enter solute concentration:");
                    fscanf(infile,"%lf",&state.solute);
		    (void) printf("%f\n",state.solute);
		    FT_SetDirichletBoundary(front,NULL,NULL,NULL,
				(POINTER)&state,hs);
		    break;
		case 'f':			// Flow through state
		case 'F':
		    FT_SetDirichletBoundary(front,ss_flowThroughBoundaryState,
				"ss_flowThroughBoundaryState",NULL,NULL,hs);
		    break;
		}
	    }
	}
	fclose(infile);
}	/* end read_ss_dirichlet_bdry_data */

void init_fluid_state_func(
        Incompress_Solver_Smooth_Basis *l_cartesian,
        PROB_TYPE prob_type)
{
	l_cartesian->getInitialState = zero_state;
}	/* end init_fluid_state_func */

static double zero_state(
        COMPONENT comp,
        double *coords,
        L_STATE& state,
        int dim,
        IF_PARAMS *iFparams)
{
        int i;
        for (i = 0; i < dim; ++i)
            state.m_U[i] = 0.0;
}       /* end zero_state */


static double (*getStateVel[MAXD])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};

extern void ss_flowThroughBoundaryState(
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
	CRT_PARAMS *cRparams = (CRT_PARAMS*)front->extra2;
	IF_FIELD *iF_field = iFparams->field;
	CRT_FIELD *cR_field = cRparams->field;
	double dir[MAXD];
	double u[3];		/* velocity in the sweeping direction */
	double v[3][MAXD];	/* velocity in the orthogonal direction */
	double vort[3];		/* vorticity stencil */
	double pres[3];		/* pressure stencil */
	double solu[3];		/* solute stencil */
	double f_u;		/* u flux in the sweeping direction */
	double f_v[MAXD];	/* v flux in the orthogonal direction */
	double f_vort;		/* vort flux */
	double f_pres;		/* pressure flux */
	double f_solu;		/* solute flux */
	double dn,dt = front->dt;
	STATE *newst = (STATE*)state;
	STATE  **sts;
	int i,j,dim = front->rect_grid->dim;
	int nrad = 2;
	
	if (debugging("flow_through"))
	    printf("Entering ss_flowThroughBoundaryState()\n");

	tsten = FrontGetTanStencils(front,oldp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = tsten[0]->dir[i];
	dn = FT_GridSizeInDir(dir,front);

	if (comp == negative_component(hs))  
	    sts = (STATE**)tsten[0]->leftst;
	else 
	    sts = (STATE**)tsten[0]->rightst;

	if (debugging("flow_through"))
	{
	    printf("Ambient component: %d\n",comp);
	    printf("hs = %d  oldp->hs = %d\n",hs,oldp->hs);
	    printf("Time step = %f  Tangential grid size = %f\n",dt,dn);
	    printf("Tangential direction: ");
	    for (j = 0; j < dim; ++j)
		printf("%f ",tsten[0]->dir[j]);
	    printf("\n");
	    printf("Tan_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    printf("Left points:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",Coords(tsten[0]->p[-i])[j]);
		printf("\n");
	    }
	    printf("Right points:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",Coords(tsten[0]->p[i])[j]);
		printf("\n");
	    }
	}

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 3; ++j)
	{
	    vort[j] = sts[j-1]->vort;
	    pres[j] = sts[j-1]->pres;
	    solu[j] = sts[j-1]->solute;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += sts[j-1]->vel[i]*dir[i];
		v[j][i] = sts[j-1]->vel[i]*(1.0 - dir[i]);
	    }
	}

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_solu = linear_flux(u[1],solu[0],solu[1],solu[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] = sts[0]->vel[i] - dt/dn*(
		f_u*dir[i] + f_v[i]) ;
	newst->vort = sts[0]->vort - dt/dn*f_vort;
	newst->pres = sts[0]->pres - dt/dn*f_pres;
	newst->solute = sts[0]->solute - dt/dn*f_solu;
	
	nsten = FT_CreateNormalStencil(front,oldp,comp,nrad);
	for (i = 0; i < dim; ++i)
	    dir[i] = nsten->nor[i];
	dn = FT_GridSizeInDir(dir,front);

	if (debugging("flow_through"))
	{
	    printf("Time step = %f  Normal grid size = %f\n",dt,dn);
	    printf("Normal direction: ");
	    for (j = 0; j < dim; ++j)
		printf("%f ",nsten->nor[j]);
	    printf("\n");
	    printf("Nor_stencil at point p(%f %f)\n",Coords(oldp)[0],
				Coords(oldp)[1]);
	    printf("Nor_stencil:\n");
	    for (i = 0; i < nrad; ++i)
	    {
		for (j = 0; j < dim; ++j)
	    	    printf("%f ",nsten->pts[i][j]);
		printf("\n");
	    }
	}

	for (j = 0; j < 3; ++j)
	    u[j] = 0.0;
	for (j = 0; j < 2; ++j)
	{
	    vort[j] = sts[0]->vort;
	    pres[j] = sts[0]->pres;
	    solu[j] = sts[0]->solute;
	    for (i = 0; i < dim; ++i)
	    {
		u[j] += sts[0]->vel[i]*dir[i];
		v[j][i] = sts[0]->vel[i]*(1.0 - dir[i]);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    double vtmp;
	    if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			iF_field->vel[i],getStateVel[i],&vtmp))
            {
        	vtmp = sts[0]->vel[i];
            }
	    u[2] += vtmp*dir[i];
	    v[2][i] = sts[1]->vel[i] - vtmp*dir[i];
	}
	if (dim == 2)
	{
	    if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],
			iF_field->vort,getStateVort,&vort[2]))
	    	vort[2] = sts[1]->vort;
	}
	if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],iF_field->pres,
                            getStatePres,&pres[2]))
	    pres[2] = sts[1]->pres;
	if (!FT_IntrpStateVarAtCoords(front,comp,nsten->pts[1],cR_field->solute,
                            getStateSolute,&solu[2]))
	    solu[2] = sts[1]->solute;

	f_u = burger_flux(u[0],u[1],u[2]);
	for (i = 0; i < dim; ++i)
	    f_v[i] = linear_flux(u[1],v[0][i],v[1][i],v[2][i]);
	f_vort = linear_flux(u[1],vort[0],vort[1],vort[2]);
	f_pres = linear_flux(u[1],pres[0],pres[1],pres[2]);
	f_solu = linear_flux(u[1],solu[0],solu[1],solu[2]);

	for (i = 0; i < dim; ++i)
	    newst->vel[i] += - dt/dn*(f_u*dir[i] + f_v[i]) ;
	newst->vort += - dt/dn*f_vort;
	newst->pres += - dt/dn*f_pres;
	newst->solute += - dt/dn*f_solu;
	if (debugging("dirichlet_bdry"))
	{
	    printf("flow through boundary state:\n");
	    print_general_vector("Velocity: ",newst->vel,dim,"\n");
	    printf("Pressure: %f\n",newst->pres);
	    printf("Vorticity: %f\n",newst->vort);
	    printf("Solute: %f\n",newst->solute);
	}
}       /* end ss_flowThroughBoundaryState */

extern void read_fluid_params(
	char *inname,
	IF_PARAMS *iFparams)
{
	char string[100];
	FILE *infile = fopen(inname,"r");

	CursorAfterString(infile,"Enter density and viscosity of the fluid:");
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	iFparams->m_comp1 = CRYSTAL_COMP;
	iFparams->m_comp2 = LIQUID_COMP2;
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);
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
}	/* end read_fluid_properties */

void initFrontStates(
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
	    {
		if (wave_type(hs) == GROWING_BODY_BOUNDARY)
                    sr->solute = cRparams->C_eq;
		else
                    sr->solute = cRparams->C_0;
	    }
            else if (positive_component(hs) == CRYSTAL_COMP)
                sr->solute = rho_s;
            else
                sr->solute = 0.0;
            if (negative_component(hs) == SOLUTE_COMP)
	    {
		if (wave_type(hs) == GROWING_BODY_BOUNDARY)
                    sr->solute = cRparams->C_eq;
		else
                    sr->solute = cRparams->C_0;
	    }
            else if (negative_component(hs) == CRYSTAL_COMP)
                sl->solute = rho_s;
            else
                sl->solute = 0.0;
        }
}       /* end initFrontStates */


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

extern void fluid_print_front_states(
	FILE *outfile,
	Front *front)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int dim = intfc->dim;

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
}	/* end fluid_print_front_states */

extern void fluid_read_front_states(
	FILE *infile,
	Front *front)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *lstate,*rstate;
	int dim = intfc->dim;

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
}	/* end fluid_read_front_states */

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
	case NEUMANN_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
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
}       /* end isDirichletPresetBdry */
