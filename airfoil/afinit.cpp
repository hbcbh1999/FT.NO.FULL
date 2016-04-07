#include <iFluid.h>
#include <airfoil.h>

typedef struct {
        double N[MAXD];         /* normal of the plane */
        double P[MAXD];         /* a point on the plane */
        double cen[MAXD];       /* center of the hole */
        double radius;          /* radius of the hole */
} CONSTR_PARAMS;

static double zero_state(COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
static void setInitialIntfc2d(Front*,LEVEL_FUNC_PACK*,char*,AF_PROB_TYPE);
static void setInitialIntfc3d(Front*,LEVEL_FUNC_PACK*,char*,AF_PROB_TYPE);
static boolean parachute_constr_func(POINTER,double*);
static int zero_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				double*);
static int unfolding_velo(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				double*);

extern void read_af_prob_type(
        char *inname,
        AF_PROB_TYPE *prob_type)
{
        char string[100];
        FILE *infile = fopen(inname,"r");

        *prob_type = ERROR_TYPE;
        CursorAfterString(infile,"Enter problem type:");
        fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
        if (string[0] == 'D' || string[0] == 'd')
        {
            *prob_type = DOUBLE_VORTEX_VELO;
        }
        else if (string[0] == 'F' || string[0] == 'f')
        {
            *prob_type = FLIP_FLAP;
        }
        else if (string[0] == 'L' || string[0] == 'l')
        {
            *prob_type = LEAF_FALL;
        }
        else if (string[0] == 'P' || string[0] == 'p')
        {
            *prob_type = PARACHUTE;
        }
        else if (string[0] == 'R' || string[0] == 'r')
        {
            *prob_type = REVERSAL_VELO;
        }
        else if (string[0] == 'Z' || string[0] == 'z')
        {
            *prob_type = ZERO_VELO;
        }

        assert(*prob_type != ERROR_TYPE);
        fclose(infile);
}       /* end read_af_prob_type */

extern void setInitialIntfc(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        AF_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;

	level_func_pack->wave_type = ELASTIC_BOUNDARY;
	iFparams->m_comp1 = SOLID_COMP;
        iFparams->m_comp2 = LIQUID_COMP2;

	switch (front->rect_grid->dim)
	{
	case 2:
	    return setInitialIntfc2d(front,level_func_pack,inname,prob_type);
	case 3:
	    return setInitialIntfc3d(front,level_func_pack,inname,prob_type);
	}
}	/* end setInitialIntfc */

static void setInitialIntfc3d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        AF_PROB_TYPE prob_type)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int i,dim = iFparams->dim;
	static ELLIP_PARAMS ellip_params;
	double *cen,*rad;
	static CONSTR_PARAMS constr_params;

	level_func_pack->set_3d_bdry = YES;
	level_func_pack->neg_component = LIQUID_COMP2;
        level_func_pack->pos_component = LIQUID_COMP2;	
	level_func_pack->func_params = NULL;
        level_func_pack->func = NULL;
	CursorAfterString(infile,"Enter initial canopy surface type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'E':
	case 'e':
	    cen = ellip_params.cen;
	    rad = ellip_params.rad;
	    af_params->init_surf_type = ELLIPTIC;
	    CursorAfterString(infile,
			"Enter center coordinate of the ellipsoid:");
	    fscanf(infile,"%lf %lf %lf",&cen[0],&cen[1],&cen[2]);
	    (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
	    CursorAfterString(infile,
			"Enter three radii of the ellipsoid:");
	    fscanf(infile,"%lf %lf %lf",&rad[0],&rad[1],&rad[2]);
	    (void) printf("%f %f %f\n",rad[0],rad[1],rad[2]);
	    level_func_pack->is_mono_hs = YES;
	    level_func_pack->func = ellipsoid_func;
	    level_func_pack->func_params = (POINTER)&ellip_params;
	    level_func_pack->func = ellipsoid_func;
	    level_func_pack->constr_func = parachute_constr_func;
	    level_func_pack->constr_params = (POINTER)&constr_params;
	    constr_params.N[0] = constr_params.N[1] = 0.0;
	    constr_params.N[2] = 1.0;
	    constr_params.P[0] = cen[0];
	    constr_params.P[1] = cen[1];
	    CursorAfterString(infile,
			"Enter the height of canopy boundary:");
	    fscanf(infile,"%lf",&constr_params.P[2]);
	    (void) printf("%lf\n",constr_params.P[2]);
	    constr_params.radius = 0.0;
	    CursorAfterString(infile,"Enter yes to cut a hole on canopy:");
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	    {
	    	constr_params.cen[0] = cen[0];
	    	constr_params.cen[1] = cen[1];
	    	CursorAfterString(infile,"Enter radius of the hole:");
	    	fscanf(infile,"%lf",&constr_params.radius);
	    	(void) printf("%lf\n",constr_params.radius);
	    }
	    break;
	case 'V':
	case 'v':
	    af_params->init_surf_type = READ_FROM_VTK;
	    break;
	case 'N':
	case 'n':
	    break;
	}

	CursorAfterString(infile,
                        "Enter density and viscosity of the fluid:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	    (void) printf("%f ",iFparams->gravity[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);
	fclose(infile);
}	/* end setInitialIntfc3d */

static void setInitialIntfc2d(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        AF_PROB_TYPE prob_type)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int dim = iFparams->dim;

	/*TMP*/
        int i,np;
        double x,y,phi;

        level_func_pack->num_points = np = 251;         //myex num points
	/*END TMP*/

	FT_MatrixMemoryAlloc((POINTER*)&level_func_pack->point_array,np,2,sizeof(double));
        level_func_pack->is_closed_curve = NO;
	level_func_pack->neg_component = LIQUID_COMP2;
        level_func_pack->pos_component = LIQUID_COMP2;	
	level_func_pack->func_params = NULL;
        level_func_pack->func = NULL;
	CursorAfterString(infile,"Enter initial canopy surface type:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 'E':
	case 'e':
	    af_params->init_surf_type = ELLIPTIC;
	    for (i = 0; i < np; ++i)
            {
                phi = i*PI/(double)(np-1);
                level_func_pack->point_array[i][0] = 0.5 + 0.20*cos(phi);
                level_func_pack->point_array[i][1] = 0.55 + 0.40*sin(phi);
            }
	    break;
	case 'S':
	case 's':
	    af_params->init_surf_type = SINE_WAVE;
	    for (i = 0; i < np; ++i)
            {
                x = 0.4 + i*0.2/(double)(np-1);
                phi = 9.0*PI*(x - 0.4)/0.2;
                y = 0.8 + 0.04*sin(phi);
                level_func_pack->point_array[i][0] = x;
                level_func_pack->point_array[i][1] = y;
            }
	    break;
	}
	CursorAfterString(infile,
                        "Enter density and viscosity of the fluid:");
        fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	CursorAfterString(infile,"Enter gravity:");
        for (i = 0; i < dim; ++i)
	{
            fscanf(infile,"%lf",&iFparams->gravity[i]);
	    (void) printf("%f ",iFparams->gravity[i]);
	}
	(void) printf("\n");
        CursorAfterString(infile,"Enter surface tension:");
        fscanf(infile,"%lf",&iFparams->surf_tension);
	(void) printf("%f\n",iFparams->surf_tension);
        CursorAfterString(infile,"Enter factor of smoothing radius:");
        fscanf(infile,"%lf",&iFparams->smoothing_radius);
	(void) printf("%f\n",iFparams->smoothing_radius);
	fclose(infile);
}	/* end setInitialIntfc2d */

extern void init_fluid_state_func(
        Incompress_Solver_Smooth_Basis *l_cartesian,
        AF_PROB_TYPE prob_type)
{
        l_cartesian->getInitialState = zero_state;
}       /* end init_fluid_state_func */

static double zero_state(
        COMPONENT comp,
        double *coords,
        L_STATE& state,
        int dim,
        IF_PARAMS *af_params)
{
        int i;
        for (i = 0; i < dim; ++i)
            state.m_U[i] = 0.0;
        state.m_U[1] = 0.0;
}       /* end zero_state */

static boolean parachute_constr_func(
        POINTER params,
        double *coords)
{
        CONSTR_PARAMS *constr_params = (CONSTR_PARAMS*)params;
        double *N,v[MAXD],dist;
        int i;
	double *cen = constr_params->cen;

        N = constr_params->N;
        for (i = 0; i < 3; ++i)
            v[i] = coords[i] - constr_params->P[i];
        dist = Dot3d(N,v);
        if (dist < 0.0) return NO;
	dist = sqrt(sqr(coords[0] - cen[0]) + sqr(coords[1] - cen[1]));
	return (dist > constr_params->radius) ? YES : NO;
}       /* end parachute_constr_func */

extern void initVelocityFunc(
	Front *front,
	AF_PROB_TYPE prob_type)
{
	VELO_FUNC_PACK velo_func_pack;
	static VORTEX_PARAMS *vortex_params; /* velocity function parameters */
        static BIPOLAR_PARAMS *dv_params;

	switch (prob_type)
        {
        case REVERSAL_VELO:
	    FT_ScalarMemoryAlloc((POINTER*)&vortex_params,sizeof(VORTEX_PARAMS));
            front->max_time = 0.4;
            front->movie_frame_interval = 0.02;
            vortex_params->dim = 2;
            vortex_params->type[0] = 'M';
            vortex_params->cos_time = 0;
            vortex_params->cen[0] = 0.5;
            vortex_params->cen[1] = 0.25;
            vortex_params->rad = 0.15;
            vortex_params->time = 0.5*front->max_time;
            velo_func_pack.func_params = (POINTER)vortex_params;
            velo_func_pack.func = vortex_vel;
            break;
        case DOUBLE_VORTEX_VELO:
	    FT_ScalarMemoryAlloc((POINTER*)&dv_params,sizeof(BIPOLAR_PARAMS));
            dv_params->cen1[0] = 0.25;
            dv_params->cen1[1] = 0.25;
            dv_params->cen2[0] = 0.75;
            dv_params->cen2[1] = 0.25;
            dv_params->i1 = -0.5;
            dv_params->i2 =  0.5;
            velo_func_pack.func_params = (POINTER)dv_params;
            velo_func_pack.func = double_vortex_vel;
            break;
        case FLIP_FLAP:
            velo_func_pack.func_params = NULL;
            velo_func_pack.func = unfolding_velo;
            break;
        case ZERO_VELO:
            velo_func_pack.func_params = NULL;
            velo_func_pack.func = zero_velo;
            break;
        }	
	FT_InitVeloFunc(front,&velo_func_pack);
}	/* end initVelocityFunc */

static int zero_velo(
	POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
	vel[0] = vel[1] = vel[2] = 0.0;
}	/* end zero_velo */

static int unfolding_velo(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        double *vel)
{
        double *coords = Coords(p);
	double dist,v_vert;
	int i,dim = front->rect_grid->dim;

	dist = sqrt(sqr(coords[0] - 0.5) + sqr(coords[1] - 0.5));
	v_vert = 0.5*cos(dist/0.3*2.0*PI)*sin(front->time/0.15*PI);
	if (front->time > 0.3)
	{
	    for (i = 0; i < dim; ++i)
		vel[i] = 0.0;
	}
	else
	{
	    for (i = 0; i < dim-1; ++i)
		vel[i] = 0.0;
	    vel[dim-1] = v_vert;
	}
}       /* end unfolding_velo */

