#include <iFluid.h>
#include <rgbody.h>

static void make_rigid_body(Front*,char*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);
static void make_rigid_body2d(Front*,char*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);
static void make_rigid_body3d(Front*,char*,POINTER,
		double(*func)(POINTER,double*),COMPONENT,COMPONENT,int);
static void prompt_for_rigid_body_params(int,char*,RG_PARAMS*);
static void set_rgbody_params(RG_PARAMS,HYPER_SURF*);
static POINTER get_propeller_params(FILE *);
static double zero_state(COMPONENT,double*,L_STATE&,int,IF_PARAMS*);
static void initRotorIntfc(Front*,LEVEL_FUNC_PACK*,char*,RG_PROB_TYPE);
static void initWindMillIntfc3d(Front*,LEVEL_FUNC_PACK*,char*,RG_PROB_TYPE);
/*TMP*/
static	void insert_surface_tris_into_another(SURFACE*,SURFACE*);

static void initRotorIntfc(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	RG_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	INTERFACE *intfc = front->interf;
        RECT_GRID *gr = front->rect_grid;
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        int i,j,num_moving_body,dim = gr->dim;
        char msg[100],s[100];
        double  (*func)(POINTER,double*);
        POINTER func_params;
	static CIRCLE_PARAMS circle_params;
	CURVE **wall,**contact;
	int num_segs;

	circle_params.dim = 2;
	circle_params.add_plan_surf = NO;
	CursorAfterString(infile,"Enter center of the cylinder:");
	fscanf(infile,"%lf  %lf",&circle_params.cen[0],
				&circle_params.cen[1]);
	(void) printf("%f %f\n",circle_params.cen[0],circle_params.cen[1]);
	CursorAfterString(infile,"Enter radius of the cylinder:");
	fscanf(infile,"%lf",&circle_params.R);
	(void) printf("%f\n",circle_params.R);
	func_params = (POINTER)&circle_params;
        func = level_circle_func;
	neg_comp = (prob_type == ROTOR_ONE_FLUID) ?
				LIQUID_COMP2 : LIQUID_COMP1;
	pos_comp = SOLID_COMP;
	wall = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,front->interf,
			neg_comp,pos_comp,func,func_params,NEUMANN_BOUNDARY,
			&num_segs);
	if (prob_type == ROTOR_TWO_FLUID)
	{
	    CursorAfterString(infile,
			"Enter radius of the two fluid interface:");
	    fscanf(infile,"%lf",&circle_params.R);
	    (void) printf("%f\n",circle_params.R);
	    neg_comp = LIQUID_COMP2;
	    pos_comp = LIQUID_COMP1;
	    contact = (CURVE**)FT_CreateLevelHyperSurfs(front->rect_grid,
			front->interf,neg_comp,pos_comp,func,func_params,
			FIRST_PHYSICS_WAVE_TYPE,&num_segs);
	    sprintf(msg,"Enter density and viscosity of fluid 1:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf %lf",&iFparams->rho1,&iFparams->mu1);
	    (void) printf("%f %f\n",iFparams->rho1,iFparams->mu1);
	    sprintf(msg,"Enter density and viscosity of fluid 2:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	    (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	    neg_comp = LIQUID_COMP2;
	    pos_comp = SOLID_COMP;
	}
	else
	{
	    sprintf(msg,"Enter density and viscosity of the fluid:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	    (void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);
	}

	sprintf(msg,"Enter the shape of the rotor:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	switch(s[0])
	{
	case 'p':
        case 'P':
	    func_params = get_propeller_params(infile);
            func = propeller_func;
	    make_rigid_body(front,inname,func_params,func,
				pos_comp,neg_comp,0);
	    break;
	}
	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initRotorIntfc */

static void initWindMillIntfc3d(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	RG_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	INTERFACE *intfc = front->interf;
        RECT_GRID *gr = front->rect_grid;
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        int i,j,num_moving_body,dim = gr->dim;
        char msg[100],s[100];
	SURFACE *wing,*rotor,*support;
	char vtk_name[100],string[100];
	FILE *vtk_file;
	RG_PARAMS rgb_params;

	sprintf(string,"Enter the vtk file name for wing:");
	CursorAfterString(infile,string);
	fscanf(infile,"%s",vtk_name);
	(void) printf("%s\n",vtk_name);
	neg_comp = SOLID_COMP;
	pos_comp = LIQUID_COMP2;
	read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name,&wing);
	wave_type(wing) = MOVABLE_BODY_BOUNDARY;
	if (debugging("windmill"))
	{
	    (void) printf("Wing surface read from vtk file %s\n",vtk_name);
	    (void) gview_plot_interface("vtk_wing",intfc);
	}

	sprintf(string,"Enter the vtk file name for rotor:");
        CursorAfterString(infile,string);
        fscanf(infile,"%s",vtk_name);
        (void) printf("%s\n",vtk_name);
        neg_comp = SOLID_COMP;
        pos_comp = LIQUID_COMP2;
        read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name,&rotor);
        wave_type(rotor) = MOVABLE_BODY_BOUNDARY;
        if (debugging("windmill"))
        {
            (void) printf("Rotor surface read from vtk file %s\n",vtk_name);
            (void) gview_plot_interface("vtk_rotor",intfc);
        }
	insert_surface_tris_into_another(wing,rotor);
	delete_surface(rotor);

	prompt_for_rigid_body_params(dim,inname,&rgb_params);
	rgb_params.dim = 3;
	set_rgbody_params(rgb_params,Hyper_surf(wing));

	if (debugging("windmill"))
            (void) printf("Passed prompt_for_rigid_body()\n");

	sprintf(string,"Enter the vtk file name for support structure:");
	CursorAfterString(infile,string);
	fscanf(infile,"%s",vtk_name);
	(void) printf("%s\n",vtk_name);
	neg_comp = SOLID_COMP;
	pos_comp = LIQUID_COMP2;
	read_vtk_surface(intfc,neg_comp,pos_comp,vtk_name,&support);
	wave_type(support) = NEUMANN_BOUNDARY;
	if (debugging("windmill"))
	{
	    (void) printf("Support surface read from vtk file %s\n",vtk_name);
	    (void) gview_plot_interface("vtk_full",intfc);
	    if (consistent_interface(intfc))
	    {
		(void) printf("Interface read from vtk files is consistent\n");
		(void) print_interface(intfc);
	    }
	    else
	    {
		(void) printf("Interface read from vtk files not consistent\n");
		clean_up(ERROR);
	    }
	}
	sprintf(msg,"Enter density and viscosity of the fluid:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initWindMillIntfc3d */

static void initFluidRgb(
	Front *front,
	LEVEL_FUNC_PACK *level_func_pack,
	char *inname,
	RG_PROB_TYPE prob_type)
{
	IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
	FILE *infile = fopen(inname,"r");
	INTERFACE *intfc = front->interf;
        RECT_GRID *gr = front->rect_grid;
        COMPONENT neg_comp;
        COMPONENT pos_comp;
        int i,j,num_moving_body,dim = gr->dim;
        char msg[100],s[100];
        double  (*func)(POINTER,double*);
        POINTER func_params;
	PROPELLER_PARAMS propeller_params;

	sprintf(msg,"Enter density and viscosity of the fluid:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%lf %lf",&iFparams->rho2,&iFparams->mu2);
	(void) printf("%f %f\n",iFparams->rho2,iFparams->mu2);

	neg_comp = LIQUID_COMP2;
	pos_comp = SOLID_COMP;

	sprintf(msg,"Enter the shape of the rigid body:");
        CursorAfterString(infile,msg);
        fscanf(infile,"%s",s);
        (void) printf("%s\n",s);
        switch(s[0])
        {
        case 'p':
        case 'P':
	    func_params = get_propeller_params(infile);
            func = propeller_func;
	    make_rigid_body(front,inname,func_params,func,
				pos_comp,neg_comp,0);
	    break;
	}
	fclose(infile);
	FT_ParallelExchIntfcBuffer(front);
}	/* end initFluidRgb by JDKim */

extern void init_moving_bodies(
        Front *front,
        LEVEL_FUNC_PACK *level_func_pack,
        char *inname,
        RG_PROB_TYPE prob_type)
{
        IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
        iFparams->m_comp2 = LIQUID_COMP2;
	switch (prob_type)
        {
	case ROTOR_ONE_FLUID:
            iFparams->m_comp1 = SOLID_COMP;
	case ROTOR_TWO_FLUID:
	    initRotorIntfc(front,level_func_pack,inname,prob_type);
	    break;
	case FLUID_RIGID_BODY:
            iFparams->m_comp1 = SOLID_COMP;
	    initFluidRgb(front,level_func_pack,inname,prob_type);
	    break;
	case WINDMILL_3D:
            iFparams->m_comp1 = SOLID_COMP;
	    initWindMillIntfc3d(front,level_func_pack,inname,prob_type);
	    break;
	}
}	/* init_moving_bodies */

extern void init_fluid_state_func(
	Incompress_Solver_Smooth_Basis *cartesian, 
	RG_PROB_TYPE prob_type)
{
	cartesian->getInitialState = zero_state;
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

static void make_rigid_body(
	Front *front,
	char *inname,
	POINTER func_params,
	double  (*func)(POINTER,double*),
	COMPONENT neg_comp,
	COMPONENT pos_comp,
	int id)
{
	int dim = front->rect_grid->dim;
	switch (dim)
	{
	case 2:
	    make_rigid_body2d(front,inname,func_params,func,neg_comp,
				pos_comp,id);
	    return;
	case 3:
	    make_rigid_body3d(front,inname,func_params,func,neg_comp,
				pos_comp,id);
	    return;
	}
}	/* end make_rigid_body */

static void make_rigid_body2d(
	Front *front,
	char *inname,
	POINTER func_params,
	double  (*func)(POINTER,double*),
	COMPONENT neg_comp,
	COMPONENT pos_comp,
	int id)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = front->rect_grid;
	int i,j,dim = intfc->dim;
	int num_segs;
	char msg[100],s[100];
	CURVE **curves;
	RG_PARAMS rg_params;

	curves = (CURVE**)FT_CreateLevelHyperSurfs(gr,intfc,neg_comp,pos_comp,
			func,func_params,MOVABLE_BODY_BOUNDARY,&num_segs);
	if (curves == NULL || num_segs == 0) return;
	prompt_for_rigid_body_params(dim,inname,&rg_params);

	rg_params.dim = 2;
	for (j = 0; j < num_segs; ++j)
	{
	    body_index(curves[j]) = id;
	    set_rgbody_params(rg_params,Hyper_surf(curves[j]));
	}
}	/* end make_rigid_body2d */

static void make_rigid_body3d(
	Front *front,
	char *inname,
	POINTER func_params,
	double  (*func)(POINTER,double*),
	COMPONENT neg_comp,
	COMPONENT pos_comp,
	int id)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = front->rect_grid;
	int i,j,dim = intfc->dim;
	int num_segs;
	char msg[100],s[100];
	SURFACE *surf,**ps;
	RG_PARAMS rg_params;

	ps = (SURFACE**)FT_CreateLevelHyperSurfs(gr,intfc,neg_comp,pos_comp,
			func,func_params,MOVABLE_BODY_BOUNDARY,&num_segs);
	if (ps == NULL || num_segs == 0) return;
	/*
	if (!make_level_surface(gr,intfc,neg_comp,pos_comp,func,func_params,
				&surf))
	    return;
        wave_type(surf) = MOVABLE_BODY_BOUNDARY;
	*/
	prompt_for_rigid_body_params(dim,inname,&rg_params);
	body_index(surf) = id;
	rg_params.dim = 3;
	set_rgbody_params(rg_params,Hyper_surf(surf));
}	/* end make_rigid_body3d */

static POINTER get_propeller_params(FILE *infile)
{
        static PROPELLER_PARAMS pparams;

        CursorAfterString(infile,"Enter the center:");
        fscanf(infile,"%lf %lf",&pparams.x[0],&pparams.y[0]);
	(void) printf("%f %f\n",pparams.x[0],pparams.y[0]);
        CursorAfterString(infile,"Enter the number of wings:");
        fscanf(infile,"%d",&pparams.NofW);
	(void) printf("%d\n",pparams.NofW);
        CursorAfterString(infile,"Enter the inner radius:");
        fscanf(infile,"%lf",&pparams.r[0]);
	(void) printf("%f\n",pparams.r[0]);
        CursorAfterString(infile,"Enter the outer radius:");
        fscanf(infile,"%lf",&pparams.r[1]);
	(void) printf("%f\n",pparams.r[1]);
	if (debugging("init_rigid_body"))
	{
	    printf("Propeller shaped rigid body\n");
	    printf("Center: %f %f\n",pparams.x[0],pparams.y[0]);
	    printf("Number of wings: %d\n",pparams.NofW);
	    printf("Inner radius: %f\n",pparams.r[0]);
	    printf("Outer radius: %f\n",pparams.r[1]);
	}
        return (POINTER)&pparams;
}	/* end get_propeller_params */

static	void prompt_for_rigid_body_params(
	int dim,
	char *inname,
	RG_PARAMS *rgb_params)
{
	int i,j;
	char msg[100],s[100];
	FILE *infile = fopen(inname,"r");

	if (debugging("rgbody")) 
	    (void) printf("Enter prompt_for_rigid_body_params()\n");

	sprintf(msg,"Enter the total mass for rigid body:");
	CursorAfterString(infile,msg);
	fscanf(infile,"%lf",&rgb_params->total_mass);
	(void) printf("%f\n",rgb_params->total_mass);
	sprintf(msg,"Enter the center of mass for rigid body:");
	CursorAfterString(infile,msg);
	for (i = 0; i < dim; ++i)
	{
	    fscanf(infile,"%lf",&rgb_params->center_of_mass[i]);
	    (void) printf("%f ",rgb_params->center_of_mass[i]);
	}
	(void) printf("\n");
	CursorAfterString(infile,
		"Type yes if rigid body will only rotate about an axis:");
	fscanf(infile,"%s",s);
	(void) printf("%s\n",s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    if (dim == 3)
	    {
		double mag_dir = 0.0;
	    	sprintf(msg,"Enter direction of the axis:");
		CursorAfterString(infile,msg);
		for (i = 0; i < dim; ++i)
		{
	    	    fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
	    	    (void) printf("%f ",rgb_params->rotation_dir[i]);
		    mag_dir += sqr(rgb_params->rotation_dir[i]);
		}
		mag_dir = sqrt(mag_dir);
		for (i = 0; i < dim; ++i)
		    rgb_params->rotation_dir[i] /= mag_dir;
		(void) printf("\n");
	    }

	    sprintf(msg,"Enter center of the axis:");
	    CursorAfterString(infile,msg);
	    for (i = 0; i < dim; ++i)
	    {
	    	fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
	    	(void) printf("%f ",rgb_params->rotation_cen[i]);
	    }
	    (void) printf("\n");

	    sprintf(msg,"Enter the moment of inertial about the axis:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf",&rgb_params->moment_of_inertial);
	    (void) printf("%f\n",rgb_params->moment_of_inertial);

	    CursorAfterString(infile,
			"Type yes if angular velocity is preset: ");
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	    	rgb_params->motion_type = PRESET_MOTION;
	        CursorAfterString(infile,"Enter preset angular velocity: ");
	    }
	    else
	    {
	    	rgb_params->motion_type = ROTATION;
	        CursorAfterString(infile,"Enter initial angular velocity: ");
	    }
	    fscanf(infile,"%lf",&rgb_params->angular_velo);
	    (void) printf("%f\n",rgb_params->angular_velo);
	}
        else
        {
	    sprintf(msg,"Enter the moment of inertial about center of mass:");
	    CursorAfterString(infile,msg);
	    fscanf(infile,"%lf",&rgb_params->moment_of_inertial);
	    (void) printf("%f\n",rgb_params->moment_of_inertial);

            CursorAfterString(infile,
			"Type yes if you want vertical motion only?: ");
	    fscanf(infile,"%s",s);
	    (void) printf("%s\n",s);
            if (s[0] == 'y' || s[0] == 'Y')
	    	rgb_params->motion_type = VERTICAL_MOTION;
            else
	    	rgb_params->motion_type = FREE_MOTION;

	    sprintf(msg,"Enter the initial center of mass velocity:");
	    CursorAfterString(infile,msg);
	    for (i = 0; i < dim; ++i)
	    {
	    	fscanf(infile,"%lf",&rgb_params->cen_of_mass_velo[i]);
	    	(void) printf("%f ",rgb_params->cen_of_mass_velo[i]);
	    }

        }
	    
	if (debugging("rgbody")) 
	    (void) printf("Leaving prompt_for_rigid_body_params()\n");
}	/* end prompt_for_rigid_body_params */

static	void insert_surface_tris_into_another(
	SURFACE *s1,
	SURFACE *s2)
{
	int i,n,num_tri = s2->num_tri;
	TRI *tri,**tris;
	FT_VectorMemoryAlloc((POINTER*)&tris,num_tri,sizeof(TRI));
	n = 0;
	for (tri = first_tri(s2); !at_end_of_tri_list(tri,s2); tri = tri->next)
	    tris[n++] = tri;
	for (i = 0; i < n; ++i)
	    insert_tri_at_tail_of_list(tris[i],s1);
	FT_FreeThese(1,tris);
	
}	/* end insert_surface_tris_into_another */

static void set_rgbody_params(
	RG_PARAMS rg_params,
	HYPER_SURF *hs)
{
	int i,dim = rg_params.dim;
	total_mass(hs) = rg_params.total_mass;
	mom_inertial(hs) = rg_params.moment_of_inertial;
	angular_velo(hs) = rg_params.angular_velo;
	motion_type(hs) = rg_params.motion_type;
        surface_tension(hs) = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    center_of_mass(hs)[i] = rg_params.center_of_mass[i];
	    center_of_mass_velo(hs)[i] = 
				rg_params.cen_of_mass_velo[i];
	    rotation_center(hs)[i] = 
				rg_params.rotation_cen[i];
	    if (dim == 3)
	        rotation_direction(hs)[i] = 
				rg_params.rotation_dir[i];
	}
}	/* end set_rgbody_params */
