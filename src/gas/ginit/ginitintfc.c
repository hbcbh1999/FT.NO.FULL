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
*				ginitintfc.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains initialization routines for gas dynamics.
*
*	The main routines in this file that are accessed through function
*	pointers are:
*
*	g_init_interface(), called from init() in dinit.c.
*
*	g_init_cauchy_data_pointers(), g_initializer(), g_intfc_initializer(),
*	g_restart_initializer(), and g_restart_intfc_initializer(), called
*	from init_states() in dinit.c.
*/

#include <ginit/ginit.h>

struct _DENSITY_STEP_WALL_DATA {
	OUTPUT_DATA odata;

	COMPONENT   *pcomp;
	Locstate    *exact_wstate;
	Locstate    *pstate;
	Locstate    *wave_head, *wave_tail;
	boolean        inside_wave;
	boolean        isforward;
	double       c, s;
	double       wave_tol;
	double       *exact_incomingWs, *exact_outgoingWs;
	double       *wave_head_time, *wave_tail_time, *wave_mid_time;
	double       Plast;
	double       *px;
	double       *t;
	double       DX;
	int         num_exact_wall_states;
	int         NumAllocWallEvents;
	int         num_events;
	int         num_samples;
	size_t      sizest;
};
typedef struct _DENSITY_STEP_WALL_DATA DENSITY_STEP_WALL_DATA;

	/* LOCAL Function Declarations */
LOCAL   int     prompt_for_time_dependent_boundary_state(int,Gas_param*,
                                                         INTERFACE*);
LOCAL	void	free_ambient_comp_type(COMP_TYPE*);
LOCAL	void	g_bdry_state_initializer(int,int,HYPER_SURF*,INIT_DATA*,
					 INIT_PHYSICS*);
LOCAL	void	g_initializer(double*,COMPONENT,Locstate,INTERFACE*,INIT_DATA*);
LOCAL	void	g_intfc_initializer(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				    Locstate,Locstate,INIT_DATA*);
LOCAL	void	get_state_exterior(double*,Locstate,COMP_TYPE*,
				   HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	get_state_obstacle(double*,Locstate,COMP_TYPE*,
				   HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	init_ambient(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void 	init_2d_riemann(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_tri_grid_test(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	print_comp_types(INTERFACE*);
LOCAL	void	print_component_type(const char*,COMP_TYPE*,const char*);
LOCAL	void	prompt_for_boundary_flags(INIT_DATA*,INIT_PHYSICS*,
					  const Prompt_type*);
LOCAL	void	prompt_for_contact_wall_params(const IO_TYPE*,Front*);
LOCAL	void	prompt_for_rect_boundary_flags(INIT_DATA*,INIT_PHYSICS*,
					       const Prompt_type*);
LOCAL	void	set_exterior_comp_type(COMP_TYPE*);
LOCAL	void	set_passive_curve_flags(Front*);

LOCAL	DENSITY_STEP_WALL_DATA *alloc_DENSITY_STEP_WALL_DATA(int,int,Front*);
LOCAL	void    expand_num_wall_events(DENSITY_STEP_WALL_DATA*,Front*);
LOCAL	void	g_1d_bdry_state_initializer(int,int,POINT*,INIT_DATA*,
					    INIT_PHYSICS*);
LOCAL	void	init_oned_density_step(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	prompt_for_bdry_flags1d(INIT_DATA*,INIT_PHYSICS*,
					const Prompt_type*);
LOCAL	void	print_density_step_wall_data(Grid*,Wave*,Front*,Printplot*,
                                             OUTPUT_DATA*,boolean);

LOCAL	void	g_2d_bdry_state_initializer(int,int,CURVE*,INIT_DATA*,
					    INIT_PHYSICS*);
LOCAL	void	init_bow_shock(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_expanding_shock(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	prompt_for_bdry_flags2d(INIT_DATA*,INIT_PHYSICS*,
					const Prompt_type*);
LOCAL	void	prompt_for_contact_wall_node_params(const IO_TYPE*,Front*);

LOCAL	void	g_3d_bdry_state_initializer(int,int,SURFACE*,INIT_DATA*,
					    INIT_PHYSICS*);
LOCAL	void	prompt_for_bdry_flags3d(INIT_DATA*,INIT_PHYSICS*,
					const Prompt_type*);

LOCAL	void	prompt_for_cylindrical_pencil(INIT_DATA*,INIT_PHYSICS*);

LOCAL 	void    set_time_dep_pres(BOUNDARY_STATE*,double,double,double,
			double,double,Locstate);


/*
*			g_init_interface():
*
*	Initializes the interface points, constraints, and the components.
*	Any choice of  a component which is a union of connected components
*	of domain\interface is allowed; the the union of the components must 
*	equal domain\interface.
*
*	The function pointer ip->init_interface is set to be g_init_interface
*	in set_basic_phys_parameters().
*/

/*ARGSUSED*/
EXPORT void g_init_interface(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front	  *front = ip->root->front;
	CROSS	  *cross;
	INTERFACE *intfc = front->interf;
	int	  max_n_comps = max_num_comps();
	int       ptOffset = g_init_data(init)->_promptForBoundaryFlagsOffset;
	const Prompt_type *ptypes = bdry_wave_type_prompt_type(init);
	boolean	  status;

	debug_print("init","Entered g_init_interface()\n");

	exclude_comp(COMPOBST,intfc);

	init_problem_type(ip);
	if (restart_io_type(init) == NULL)
	{
	    if (coord_system() == CYLINDRICAL) 
	        prompt_for_cylindrical_pencil(init,ip);

	    /* Set initial interface */
	    prompt_for_problem_specific_data(init,ip);

	    if (coord_system() == CYLINDRICAL) 
	    	make_vertical_axis_curve(init,ip);

	    if (debugging("init")) 
	    {
	        (void) output();
	        (void) printf("\t\t\tINTERFACE BEFORE set_boundary()\n\n");
	        print_interface(intfc);
	    }

	    if(problem_type(ip) == INJECTION_INLET && front->rect_grid->dim == 3)
	    /*if(Tracking_algorithm(front) == THREE_COMP_GRID_BASED_TRACKING) */
	    {
	        printf("#3comp track\n");
	    	clip_front_to_rect_boundary_type(front);
	    }
	    else
	    {
	        /* Set Boundary Information */
	        prompt_for_rect_boundary_flags(init,ip,ptypes);

	        /* TODO Unify 2 and 3 D */
	        if (front->rect_grid->dim == 3)
	    	    clip_front_to_rect_boundary_type(front);

		null_sides_are_consistent();
		/*DEBUG_TMP check_print_intfc("After clip_front", "init_rt", 'f',  */
	        	/*DEBUG_TMP front->interf, 1, -1, NO); */

	        if (set_boundary(intfc,front->rect_grid,FIRST_DYNAMIC_COMPONENT,
			     grid_tolerance(front->rect_grid)) !=
			     FUNCTION_SUCCEEDED)
	        {
	    	    screen("ERROR in g_init_interface(), set_boundary() failed\n");
	    	    print_interface(intfc);
	    	    clean_up(ERROR);
	        }

	        set_passive_curve_flags(front);
	        /*
	        * Second pass over boundary curves.
	        * Only wave types with index >= ptOffset are prompted for.
	        * See bdry_wave_type_prompt_type().
	        */
	        prompt_for_boundary_flags(init,ip,ptypes+ptOffset);
	    }
	    
	    if (debugging("reverse"))
	    {
	        CURVE **c;
	        for (c = front->interf->curves; c && *c; c++)
	            if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
			    invert_curve(*c);
	    }
	    
	    if (debugging("init")) 
	    {
	        (void) output();
	        (void) printf("\t\t\tINTERFACE after set_boundary()\n\n");
	        print_interface(intfc);
            }
	}
	
#if defined __MPI__
	/* for 3d restart */
	if (front->rect_grid->dim == 3 && pp_numnodes() > 1)
	    communicate_default_comp(front);
#endif /* defined __MPI__ */

	prompt_for_contact_wall_params(restart_io_type(init),front);

	if (max_component(intfc) >= max_n_comps) 
	{
	    screen("ERROR in g_init_interface(), max_component(intfc) = %d "
	           ">= max_n_comps = %d "
	           "in g_init_interface()\n",max_component(intfc),
		   max_n_comps);
	    clean_up(ERROR);
	}
	set_exterior_comp_type(comp_type(exterior_component(intfc)));

	    /* Check for intersections */

	if (debugging("plot3d") && intfc->dim == 3)
	    gview_plot_interface("plot-intfc",intfc);

	if(Tracking_algorithm(front) == GRID_FREE_TRACKING)
	{
	    status = intersections(intfc,&cross,NO);
	    status = pp_min_status(status);
	    if (status == FUNCTION_FAILED)
	    {
	        screen("ERROR in g_init_interface(), "
	           "intersections() failed\n\n\n");
	        (void) print_interface(intfc);
	        clean_up(ERROR);
	    }
	    if (cross != NULL)
	    {
	    	screen("ERROR in g_init_interface(), interface tangled\n");
	    	(void) print_number_of_tangles("initial ",intfc,cross);
	    	(void) print_intersections(cross,intfc);
	    	clean_up(ERROR);
	    }
	}

	set_params_list(intfc);

	order_interface(intfc);

	if (debugging("init")) 
	{
	    (void) output();
	    (void) printf("\t\t\tINITIAL INTERFACE:\n\n");
	    print_interface(intfc);
	}
	debug_print("init","Left g_init_interface()\n");
}		/*end g_init_interface*/

/*ARGSUSED*/
EXPORT int g_prompt_for_wave_type(
	const char	*mesg,
	INTERFACE	*intfc,
	INIT_PHYSICS	*ip)
{
	int   w_type;
	char  s[Gets_BUF_SIZE];

	if (mesg == NULL)
	    mesg = "";

	screen("Possible choices for wave types are\n"
	       "\t Forward shock wave (f), \n"
	       "\t Forward sound wave leading edge (fl), \n"
	       "\t Forward sound wave trailing edge (ft), \n"
	       "\t Backward shock wave (b), \n"
	       "\t Backward sound wave leading edge (bl), \n"
	       "\t Backward sound wave trailing edge (bt), \n"
	       "\t Contact (c), \n"
	       "\t Thin flame (tf), \n"
	       "\t Neumann boundary (n), \n"
	       "\t Movable Neumann boundary (mn), \n"
	       "\t Dirichlet boundary (d), \n"
	       "\t Periodic boundary (p), \n"
	       "\t Passive boundary (pa), \n");
	if (supports_riemann_problem_waves(ip))
	    screen("\t Riemann problem wave (rpw), \n");
	screen("\t Unspecified wave type (u, default), \n"
	       "Enter the wave type%s%s: ",(strlen(mesg) != 0) ? " " : "",mesg);
	(void) Gets(s);
	w_type = read_wave_type_from_string(s,intfc);
	(void) printf("wave type = %s\n",wave_type_as_string(w_type,intfc));
	return w_type;
}		/*end g_prompt_for_wave_type*/


EXPORT	void	g_init_problem_type(
	INIT_PHYSICS	*ip)
{
	int       dim = ip->root->front->rect_grid->dim;
	Prob_type *probtype;
	static Prob_type Probtype1d[] =
	{
	    {{"an ambient state test","AM",2,{AMBIENT_STATE_TEST}},
	        init_ambient},
	    {"a oned test", "ONED", 1, {ONED_TEST}, init_multi_layer},
	    {"a density step", "DS", 2, {DENSITY_STEP}, init_oned_density_step},
	    {{NULL, NULL, 0, {UNKNOWN_PROBLEM_TYPE}}, NULL}
	};
	static Prob_type Probtype2d[] =
	{
	    {{"an ambient state test","AM",2,{AMBIENT_STATE_TEST}},
	        init_ambient},
	    {{"a trigrid test", "TRI", 3, {TRI_GRID_TEST}}, init_tri_grid_test},
	    {{"a plane front", "P", 1, {PLANE_FRONT}}, init_multi_layer},
	    {{"a bowshock", "BO", 2, {BOWSHOCK}}, init_bow_shock},
	    {{"a 2D Riemann Problem", "R2", 2, {RIEMANN2D}}, init_2d_riemann},
#if defined(FULL_PHYSICS)
	    {{"a Meshkov instability problem","M",1,{MESHKOV}},init_meshkov},
	    {{"a shock diffraction", "D", 1,
		{SHOCK_DIFFRACTION}}, init_shock_diffraction},
	    {{"a shock transmission", "T", 1,
		{SHOCK_TRANSMISSION}}, init_shock_transmission},
	    {{"a ramp reflection problem", "RR", 2,
		{RAMP_REFLECTION}}, init_ramp_reflection},
	    {{"a contact-contact interaction", "CC", 2,
		{CC_NODE_TEST}}, init_CC_interaction},
	    {{"a Richtmyer linear theory", "RL", 2,
		{RICHTMYER_LINEAR_THEORY}}, init_multi_layer},
#endif /* defined(FULL_PHYSICS) */
	    {{"an astrophysical jet",  "AJ", 2,
		{ASTROPHYSICAL_JET}}, init_injection_inlet_jet},
	    {{"an injection inlet jet",  "IJ", 2,
		{INJECTION_INLET}}, init_injection_inlet_jet},
	    {{"a gas injection jet",  "FJ", 2,
		{INJECTION_INLET}}, init_fuel_injection_jet},
	    {{"a neutrino booster colapse",  "NB", 2,
                {NE_BOOSTER}}, init_neutrino_booster_detector},
	    {{"a Supernova simulation", "SN", 2,
		{SUPERNOVA}}, init_el_riem_prob},
	    {{"an imploding elliptical shock", "IMP", 3,
		{IMPLOSION}}, init_el_riem_prob},
	    {{"a shock running over an expanding ramp", "X", 1,
		{EXPANDING_SHOCK}}, init_expanding_shock},
	    {{"a random surface instability problem", "RS", 2,
		{RANDOM_SURFACE}}, init_random_surface},
	    {{"a shocked thermal layer", "STL", 2,
		{SHOCKED_THERMAL_LAYER}}, init_multi_layer},
	    {{"a Richtmyer-Meshkov instability problem", "RM", 2,
		{RICHTMYER_MESHKOV}}, init_multi_layer},
	    {{"a Rayleigh-Taylor instability problem", "RT", 2,
		{RAYLEIGH_TAYLOR}}, init_multi_layer},
	    {{"a bubbles and drops problem", "BD", 2,
		{BUBBLES_DROPS}}, init_multi_layer},
	    {{"a fluid rigid body interaction problem", "RGB", 3,
		{FLUID_RIGID_BODY}}, init_fluid_rigid_body},
	    {{"an expanding shells", "ES", 2,
		{EXPANDING_SHELLS}}, init_multi_layer},
	    {{"shock jet interaction", "SJ", 2,
		{SHOCK_JET}}, init_multi_layer},
	    {{"a Radial Rayleigh-Taylor instability problem",
		"Radial Rayleigh Taylor", 22,
		{RADIAL_RAYLEIGH_TAYLOR}}, init_el_riem_prob},
	    {{"a Kelvin-Helmholtz instability problem", "KH", 2,
		{KELVIN_HELMHOLTZ}}, init_kelvin_helmholtz},
	    {{NULL, NULL, 0, {UNKNOWN_PROBLEM_TYPE}}, NULL}
	};
	static	Prob_type Probtype3d[] =
	{
	    {{"an ambient state test","AM",2,{AMBIENT_STATE_TEST}},
	        init_ambient},
	    {{"a random surface instability problem", "RS", 2,
		{RANDOM_SURFACE}}, init_random_surface},
	    {{"a shocked thermal layer", "STL", 2,
		{SHOCKED_THERMAL_LAYER}}, init_multi_layer},
	    {{"a Richtmyer-Meshkov instability problem", "RM", 2,
		{RICHTMYER_MESHKOV}}, init_multi_layer},
	    {{"a Rayleigh-Taylor instability problem", "RT", 2,
		{RAYLEIGH_TAYLOR}}, init_multi_layer},
	    {{"a bubbles and drops problem", "BD", 2,
		{BUBBLES_DROPS}}, init_multi_layer},
	    {{"an expanding shells", "ES", 2,
		{EXPANDING_SHELLS}}, init_multi_layer},
	    {{"shock jet interaction", "SJ", 2,
		{SHOCK_JET}}, init_multi_layer},
	    {{"a Radial Rayleigh-Taylor instability problem",
		"Radial Rayleigh Taylor", 22,
		{RADIAL_RAYLEIGH_TAYLOR}}, init_el_riem_prob},
	    {{"a Kelvin-Helmholtz instability problem", "KH", 2,
		{KELVIN_HELMHOLTZ}}, init_kelvin_helmholtz},
	    {{"an imploding elliptical shock", "IMP", 3,
		{IMPLOSION}}, init_el_riem_prob},
	    {{"an injection inlet jet",  "IJ", 2,
		{INJECTION_INLET}}, init_injection_inlet_jet},
	    {{"a gas injection jet",  "FJ", 2,
		{INJECTION_INLET}}, init_fuel_injection_jet},
	    {{NULL, NULL, 0, {UNKNOWN_PROBLEM_TYPE}}, NULL}
	};

	switch (dim)
	{
	case 1:
	    probtype = Probtype1d;
	    break;
	case 2:
	    probtype = Probtype2d;
	    break;
	case 3:
	    probtype = Probtype3d;
	    break;
	default:
	    probtype = NULL;
	    break;
	}

	prompt_for_problem_type(ip,probtype);
}		/* end g_init_problem_type */

EXPORT void	prompt_for_problem_type(
	INIT_PHYSICS    *ip,
	Prob_type 	*Probtype)
{
        char            s[Gets_BUF_SIZE];
        int             i;

	screen("\nRequest problem type.  Current choices are\n");
	for (i = 0; Probtype[i].ptype.prompt != NULL; i++)
	{
	    screen("\t\t");
	    if (Probtype[i+1].ptype.prompt == NULL)
		screen("or ");
	    screen("%s (%s)",Probtype[i].ptype.prompt,Probtype[i].ptype.select);
	    screen("%s\n",(Probtype[i+1].ptype.prompt == NULL) ? "." : ",");
	}
	screen("\tEnter choice here: ");
	(void) Gets(s);

	problem_type(ip) = UNKNOWN_PROBLEM_TYPE;
     	for (i = 0; Probtype[i].ptype.prompt != NULL; i++)
	{
	    if (strncasecmp(s,Probtype[i].ptype.select,
			    Probtype[i].ptype.ncmp)==0)
	    {
	    	problem_type(ip) = prt_problem_type(ip->prt) =
	    		Probtype[i].ptype.type.itype;
	    	g_iphys(ip)->_prompt_for_problem_specific_data =
				Probtype[i].ppsd;
		break;
	    }
	}

	if (problem_type(ip) == UNKNOWN_PROBLEM_TYPE)
	{
	    screen("ERROR in prompt_for_problem_type(), "
	           "unrecognized problem type\n");
	    clean_up(ERROR);
	}
}		/*end prompt_for_problem_type*/

LOCAL void prompt_for_rect_boundary_flags(
	INIT_DATA	  *init,
	INIT_PHYSICS	  *ip,
	const Prompt_type *ptypes)
{
	PP_GRID    *pp_grid = ip->root->front->pp_grid;
	INTERFACE  *intfc = ip->root->front->interf;
	char	   mesg[256];
	int	   i, j, dim = intfc->dim;
	int	   b_type;
	int        me[3], *G;
	static const char *direction[3] = { "x", "y", "z"};
	static const char *side[3][2] = { {"left", "right"},
				          {"lower", "upper"},
				          {"bottom", "top"}
				        };

	G = pp_grid->gmax;
	find_Cartesian_coordinates(pp_mynode(),pp_grid,me);
	for (i = 0; i < dim; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
	        if (rect_boundary_type(intfc,i,j) != UNKNOWN_BOUNDARY_TYPE)
		    continue;
	        (void) sprintf(mesg,"for the %s boundary in the %s direction",
	    		       side[i][j],direction[i]);
	        b_type = prompt_for_bdry_wave_type(init,mesg,ptypes);
		rect_boundary_type(intfc,i,j) = b_type;
	        g_bdry_state_initializer(i,j,NULL,init,ip);
		if (((me[i]>0) && (j==0)) || ((me[i]<(G[i]-1)) && (j==1)))
		    rect_boundary_type(intfc,i,j) = SUBDOMAIN_BOUNDARY;
	        if ((b_type == SUBDOMAIN_BOUNDARY) && (j == 0))
	        {
	            j++;
	            (void) printf("Boundary type for the %s boundary ",
	        	   side[i][1]);
	            (void) printf("in the %s direction:",direction[i]);
		    if (me[i] < (G[i]-1))
	                (void) printf("SUBDOMAIN_BOUNDARY\n");
		    else
	                (void) printf("PERIODIC_BOUNDARY\n");
	            rect_boundary_type(intfc,i,j) = SUBDOMAIN_BOUNDARY;
	            g_bdry_state_initializer(i,j,NULL,init,ip);
	        }
	    }
	}
}		/*end prompt_for_rect_boundary_flags*/

LOCAL void prompt_for_boundary_flags(
	INIT_DATA	  *init,
	INIT_PHYSICS	  *ip,
	const Prompt_type *ptypes)
{
	INTERFACE	*intfc = ip->root->front->interf;

	switch (intfc->dim)
	{
	case 1:
	    prompt_for_bdry_flags1d(init,ip,ptypes);
	    break;
	case 2:
	    prompt_for_bdry_flags2d(init,ip,ptypes);
	    break;
	case 3:
	    prompt_for_bdry_flags3d(init,ip,ptypes);
	    break;/* TODO */
	}
}		/*end prompt_for_boundary_flags*/



/*
*
*			prompt_for_bdry_flags1d():
*
*	Loops over boundary points for which the wave type has not
*	been defined and prompts for the wave type.  Also calls
*
*		(*bdry_state_initializer)(idir,iside,Hyper_surf(p),init,ip);
*
*	which is expected to set boundary_state(point) or
*	boundary_state_function(point) appropriately.
*/

LOCAL void prompt_for_bdry_flags1d(
	INIT_DATA	  *init,
	INIT_PHYSICS	  *ip,
	const Prompt_type *ptypes)
{
	char		mesg[256];
	int		iside, idir;
	int		w_type;
	INTERFACE	*intfc = ip->root->front->interf;
	POINT		*p,**pp;
	double 		h = ip->root->front->rect_grid->h[0];
	static const double BOUNDARY_TOL = 0.001;  /* TOLERANCE */

#define is_left_boundary(p)						\
	(is_bdry(p) &&							\
	 fabs((p)->interface->table->rect_grid.L[0] - Coords(p)[0]) 	\
	 						< BOUNDARY_TOL*h)
#define is_right_boundary(p)						\
	(is_bdry(p) &&							\
	 fabs((p)->interface->table->rect_grid.U[0] - Coords(p)[0])	\
	 						< BOUNDARY_TOL*h)

	if (ptypes == NULL) return;

	idir = 0;
	for (pp = intfc->points; pp && (*pp); pp++)
	{
	    p = *pp;
	    if (!is_bdry(p)) continue;
	    /*
	    if (wave_type(p) != ERROR)
	    	continue;
	    */
	    if (is_left_boundary(p))
		iside = 0;
	    else if (is_right_boundary(p))
		iside = 1;
	    else
	    {
		iside = -1;
	    	screen("ERROR in prompt_for_bdry_flags1d(), "
	    	       "invalid boundary\n");
	    	clean_up(ERROR);
	    }
	    w_type = rect_boundary_type(intfc,idir,iside);
	    if ((w_type != UNKNOWN_BOUNDARY_TYPE) &&
		(w_type != MIXED_TYPE_BOUNDARY))
	    {
		switch (w_type)
		{
		case REFLECTION_BOUNDARY:
	            wave_type(p) = SUBDOMAIN_BOUNDARY;
		    break;
		case NO_SLIP_NEUMANN_BOUNDARY:
	            wave_type(p) = NEUMANN_BOUNDARY;
		    no_slip(Hyper_surf(p)) = YES;
		    adherence_coeff(Hyper_surf(p)) = 1.0;
		    break;
		default:
	            wave_type(p) = w_type;
		    break;
		}
	        if (w_type == DIRICHLET_BOUNDARY)
	    	    bstate_index(p) = iside;
	    }
	    else
	    {
	    	(void) sprintf(mesg, "for the boundary point (%g)",
	    		       Coords(p)[0]);
	    	wave_type(p) = prompt_for_bdry_wave_type(init,mesg,ptypes);
	    	if (wave_type(p) == UNKNOWN_BOUNDARY_TYPE)
	    	{
	    	    screen("ERROR in prompt_for_bdry_flags1d(), "
	    	           "Unknown wave type\n");
	    	    clean_up(ERROR);
	    	}
	    }
	    g_bdry_state_initializer(idir,iside,Hyper_surf(p),init,ip);
	}
#undef is_left_boundary
#undef is_right_boundary
}		/*end prompt_for_bdry_flags1d*/




/*
*
*			prompt_for_bdry_flags2d():
*
*	Loops over boundary curve for which the wave type has not
*	been defined and prompts for the wave type.  Also calls
*	g_bdry_state_initializer() which is expected to set the boundary state
*	data for the curve appropriately.  Finally, the node types of the
*	boundary nodes are set.
*/

LOCAL void prompt_for_bdry_flags2d(
	INIT_DATA	  *init,
	INIT_PHYSICS	  *ip,
	const Prompt_type *ptypes)
{
	Front		*front = ip->root->front;
	char		mesg[256];
	int		passive, fixed, bdry_type1, bdry_type2;
	int		iside, idir;
	int		w_type;
	INTERFACE	*intfc = front->interf;
	RECT_GRID	*tgr = &topological_grid(intfc);
	CURVE		*c,**pc;
	NODE		**n;

	if (ptypes == NULL)
	    return;
		
	for ((void) next_curve(intfc,NULL); next_curve(intfc,&c); )
	{
	    if (!is_bdry(c))
		continue;
	    /* By adding f_set_boundary2d(), this needs to be deleted 
	    if (wave_type(c) != ERROR)
	    	continue;
	    */
	    (void) rect_bdry_side_for_curve(&idir,&iside,c,tgr);
	    w_type = rect_boundary_type(intfc,idir,iside);
	    if ((w_type != UNKNOWN_BOUNDARY_TYPE) &&
		(w_type != MIXED_TYPE_BOUNDARY))
	    {
		switch (w_type)
		{
		case REFLECTION_BOUNDARY:
	            wave_type(c) = SUBDOMAIN_BOUNDARY;
		    break;
		case NO_SLIP_NEUMANN_BOUNDARY:
	            wave_type(c) = NEUMANN_BOUNDARY;
		    no_slip(Hyper_surf(c)) = YES;
		    adherence_coeff(Hyper_surf(c)) = 1.0;
		    break;
		default:
	            wave_type(c) = w_type;
		    break;
		}
	        if (w_type == DIRICHLET_BOUNDARY)
	        {
	    	    bstate_index(c) = 2*idir + iside;
	        }
	    }
	    else if (is_excluded_comp(negative_component(c),intfc) &&
	    		is_excluded_comp(positive_component(c),intfc))
	    {
	        wave_type(c) = PASSIVE_BOUNDARY;
	    }
	    else
	    {
	    	(void) sprintf(mesg,
			       "for the boundary from (%g, %g) to (%g, %g)",
			       Coords(c->start->posn)[0],
			       Coords(c->start->posn)[1],
			       Coords(c->end->posn)[0],
			       Coords(c->end->posn)[1]);
		wave_type(c) = prompt_for_bdry_wave_type(init,mesg,ptypes);
		if (wave_type(c) == NO_SLIP_NEUMANN_BOUNDARY)
		{
		    wave_type(c) = NEUMANN_BOUNDARY;
		    no_slip(Hyper_surf(c)) = YES;
		    adherence_coeff(Hyper_surf(c)) = 1.0;
		}
		if (wave_type(c) == UNKNOWN_BOUNDARY_TYPE)
		{
		    screen("ERROR in prompt_for_bdry_flags2d(), "
		           "Unknown wave type\n");
		    clean_up(ERROR);
		}
	    }
	    g_bdry_state_initializer(idir,iside,Hyper_surf(c),init,ip);
	}

	/* The node_type of a boundary node is set to PASSIVE_NODE if only
	 * passive curves meet there. If there is a physical curve at the
	 * node, the node type will be set to one of DIRICHLET_NODE,
	 * NEUMANN_NODE or SUBDOMAIN_NODE depending on the wave type of the
	 * boundary at that node.  If there is no physical curve at the node,
	 * it is a FIXED_NODE.
	 */

	for (n = intfc->nodes; *n; n++)
	{
	    if (node_type(*n) != ERROR)
		continue;
	    bdry_type1 = bdry_type2 = ERROR;
	    passive = YES;
	    fixed = YES;
	    if ((*n)->in_curves != NULL)
	    {
	    	for (pc = (*n)->in_curves; *pc; pc++)
	    	{
	    	    if (wave_type(*pc) != PASSIVE_BOUNDARY)
	    	    	passive = NO;
	    	    if (wave_type(*pc) >= FIRST_PHYSICS_WAVE_TYPE)
	    	    	fixed = NO;
		    else if (is_bdry(*pc))
		    {
	    	        if (bdry_type1 == ERROR)
	    	    	    bdry_type1 = wave_type(*pc);
	    	        else if (bdry_type2 == ERROR)
	    	    	    bdry_type2 = wave_type(*pc);
		    }
	    	}
	    }
	    if ((*n)->out_curves != NULL)
	    {
	    	for (pc = (*n)->out_curves; *pc; pc++)
	    	{
	    	    if (wave_type(*pc) != PASSIVE_BOUNDARY)
	    	    	passive = NO;
	    	    if (wave_type(*pc) >= FIRST_PHYSICS_WAVE_TYPE)
	    	    	fixed = NO;
		    else if (is_bdry(*pc))
		    {
	    	        if (bdry_type1 == ERROR)
	    	    	    bdry_type1 = wave_type(*pc);
	    	        else if (bdry_type2 == ERROR)
	    	    	    bdry_type2 = wave_type(*pc);
		    }
	    	}
	    }

	    if (is_bdry(*n) && (bdry_type1 != bdry_type2))
	    	fixed = YES;
			
	    if ((bdry_type1 == SUBDOMAIN_BOUNDARY) ||
		(bdry_type2 == SUBDOMAIN_BOUNDARY))
	    	node_type(*n) = SUBDOMAIN_NODE;
	    else if (passive)
	    	node_type(*n) = PASSIVE_NODE;
	    else if (fixed)
	    	node_type(*n) = FIXED_NODE;
	    else if (bdry_type1 == DIRICHLET_BOUNDARY)
	    	node_type(*n) = DIRICHLET_NODE;
	    else if (bdry_type1 == NEUMANN_BOUNDARY)
	    	node_type(*n) = NEUMANN_NODE;
	    else if (bdry_type1 == SUBDOMAIN_BOUNDARY)
	    	node_type(*n) = SUBDOMAIN_NODE;
	}
}		/*end prompt_for_bdry_flags2d*/


/*
*
*			prompt_for_bdry_flags3d():
*
*	Loops over boundary surfaces for which the wave type has not
*	been defined and prompts for the wave type.  Also calls
*	g_bdry_state_initializer() which is expected to set the boundary state
*	data for the surface appropriately.  Finally, the hsbdry types of the
*	hypersurface boundaries are set.
*/

LOCAL void prompt_for_bdry_flags3d(
	INIT_DATA	  *init,
	INIT_PHYSICS	  *ip,
	const Prompt_type *ptypes)
{
	Front		*front = ip->root->front;
	char		mesg[256];
	double		pbar[3];
	int		passive, fixed, bdry_type1, bdry_type2;
	int		iside, idir;
	int		w_type;
	INTERFACE	*intfc = front->interf;
	RECT_GRID	*tgr = &topological_grid(intfc);
	SURFACE		**s;
	CURVE		**c;

	if (ptypes == NULL)
	    return;
		
	for (s = intfc->surfaces; s && *s; s++)
	{
	    if (!is_bdry(*s))
		continue;
	    /*
	    if (wave_type(*s) != ERROR)
	    	continue;
	    */
	    rect_bdry_side_for_hyper_surf(&idir,&iside,Hyper_surf(*s),tgr);
	    w_type = rect_boundary_type(intfc,idir,iside);
	    if ((w_type != UNKNOWN_BOUNDARY_TYPE) &&
		(w_type != MIXED_TYPE_BOUNDARY))
	    {
	        wave_type(*s) = (w_type == REFLECTION_BOUNDARY) ?
				 SUBDOMAIN_BOUNDARY : w_type;
		switch (w_type)
		{
		case REFLECTION_BOUNDARY:
	            wave_type(*s) = SUBDOMAIN_BOUNDARY;
		    break;
		case NO_SLIP_NEUMANN_BOUNDARY:
	            wave_type(*s) = NEUMANN_BOUNDARY;
		    no_slip(Hyper_surf(*s)) = YES;
		    adherence_coeff(Hyper_surf(*s)) = 1.0;
		    break;
		default:
	            wave_type(*s) = w_type;
		    break;
		}
	        if (w_type == DIRICHLET_BOUNDARY)
	        {
	    	    bstate_index(*s) = 2*idir + iside;
	        }
	    }
	    else if (is_excluded_comp(negative_component(*s),intfc) &&
	    		is_excluded_comp(positive_component(*s),intfc))
	    {
	        wave_type(*s) = PASSIVE_BOUNDARY;
	    }
	    else
	    {
		int   i;
		double *h = front->rect_grid->h;
		average_position_of_surface(pbar,*s);
		for (i = 0; i < 3; i++)
		{
		    if (fabs(pbar[i]) < 0.001*h[i]) /*TOLERANCE*/
			pbar[i] = 0.0;
		}
	    	(void) sprintf(mesg,
			       "for the boundary with average position\n\t"
			       "(%g, %g, %g)",pbar[0],pbar[1],pbar[1]);
		wave_type(*s) = prompt_for_bdry_wave_type(init,mesg,ptypes);
		if (wave_type(*s) == NO_SLIP_NEUMANN_BOUNDARY)
		{
	            wave_type(*s) = NEUMANN_BOUNDARY;
		    no_slip(Hyper_surf(*s)) = YES;
		    adherence_coeff(Hyper_surf(*s)) = 1.0;
		}
		if (wave_type(*s) == UNKNOWN_BOUNDARY_TYPE)
		{
		    screen("ERROR in prompt_for_bdry_flags3d(), "
		           "Unknown wave type\n");
		    clean_up(ERROR);
		}
	    }
	    g_bdry_state_initializer(idir,iside,Hyper_surf(*s),init,ip);
	}

	/* The hsbdry_type of a boundary curve is set to PASSIVE_HSBDRY if only
	 * passive surfaces meet there. If there is a physical surface at the
	 * curve, the hsbdry_type type will be set to one of DIRICHLET_HSBDRY,
	 * NEUMANN_HSBDRY or SUBDOMAIN_HSBDRY depending on the wave type of the
	 * boundary at that curve.  If there is no physical surface at the
	 * curve, it is a FIXED_HSBDRY.
	 */

	for (c = intfc->curves; c && *c; c++)
	{
	    if (hsbdry_type(*c) != ERROR)
		continue;
	    bdry_type1 = bdry_type2 = ERROR;
	    passive = YES;
	    fixed = YES;
	    for (s = (*c)->neg_surfaces; s && *s; s++)
	    {
	    	if (wave_type(*s) != PASSIVE_BOUNDARY)
	    	    passive = NO;
	    	if (wave_type(*s) >= FIRST_PHYSICS_WAVE_TYPE)
	    	    fixed = NO;
	    	else if (bdry_type1 == ERROR)
	    	    bdry_type1 = wave_type(*s);
	    	else if (bdry_type2 == ERROR)
	    	    bdry_type2 = wave_type(*s);
	    }
	    for (s = (*c)->neg_surfaces; s && *s; s++)
	    {
	    	if (wave_type(*s) != PASSIVE_BOUNDARY)
	    	    passive = NO;
	    	if (wave_type(*s) >= FIRST_PHYSICS_WAVE_TYPE)
	    	    fixed = NO;
	    	else if (bdry_type1 == ERROR)
	    	    bdry_type1 = wave_type(*s);
	    	else if (bdry_type2 == ERROR)
	    	    bdry_type2 = wave_type(*s);
	    }

	    if (is_bdry(*c) && (bdry_type1 != bdry_type2))
	    	fixed = YES;
			
	    if (passive)
	    	hsbdry_type(*c) = PASSIVE_HSBDRY;
	    else if (fixed)
	    	hsbdry_type(*c) = FIXED_HSBDRY;
	    else if (bdry_type1 == DIRICHLET_HSBDRY)
	    	hsbdry_type(*c) = DIRICHLET_HSBDRY;
	    else if (bdry_type1 == NEUMANN_HSBDRY)
	    	hsbdry_type(*c) = NEUMANN_HSBDRY;
	    else if (bdry_type1 == SUBDOMAIN_HSBDRY)
	    	hsbdry_type(*c) = SUBDOMAIN_HSBDRY;
	}
}		/*end prompt_for_bdry_flags3d*/

/*ARGSUSED*/
EXPORT	int g_prompt_for_bdry_wave_type(
	INIT_DATA         *init,
	const char        *mesg,
	const Prompt_type *ptypes)
{
	int		  i, n;
	char		  s[Gets_BUF_SIZE];
	const Prompt_type *ptype;

	if (ptypes == NULL)
	    return UNKNOWN_BOUNDARY_TYPE;

	if (mesg == NULL)
	    mesg = "";
	for (n = 0, ptype = ptypes; ptype->prompt != NULL; ptype++, n++);
	n--;
	screen("\nEnter the boundary type -- ");
	for (i = 0; i < n; i++)
	{
	    if (i > 0 && i%3 == 0)
	    	screen("\n                           ");
	    screen("%s, ",ptypes[i].prompt);
	}
	if (i > 0 && i%3 == 0)
	    screen("\n                           ");
	screen("or %s --\n",ptypes[n++].prompt);
	screen("\t%s: ",mesg);
	(void) Gets(s);

	for (i = 0; i < n; i++)
	{
	    if (strncasecmp(s,ptypes[i].select,ptypes[i].ncmp) == 0 ||
	    	strncasecmp(s,ptypes[i].prompt,ptypes[i].ncmp) == 0)
	    	return ptypes[i].type.itype;
	}
	return UNKNOWN_BOUNDARY_TYPE;
}		/*end g_prompt_for_bdry_wave_type*/


/*ARGSUSED*/
LOCAL	void prompt_for_contact_wall_params(
	const IO_TYPE *io_type,
	Front         *front)
{
	switch (front->rect_grid->dim)
	{
	case 2:
	    prompt_for_contact_wall_node_params(io_type,front);
	    break;
	case 3:
	    break; /* TODO */
	}
}		/*end prompt_for_contact_wall_params*/


/*ARGSUSED*/
LOCAL	void prompt_for_contact_wall_node_params(
	const IO_TYPE *io_type,
	Front         *front)
{
	INTERFACE *intfc = front->interf;
	char	  s[Gets_BUF_SIZE];
	CWNP	  *cwnp;
	CURVE	  *c;
	int	  dim = intfc->dim;
	boolean	  is_neumann;
	int	  i, j;

	is_neumann = NO;
	for (i = 0; i < dim; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
	    	if (rect_boundary_type(intfc,i,j) == NEUMANN_BOUNDARY)
		    is_neumann = YES;
	    }
	}
	if (is_neumann == NO)
	{
	    (void) next_curve(intfc,NULL);
	    while (next_curve(intfc,&c))
	    {
	    	if (wave_type(c) == NEUMANN_BOUNDARY)
	    	{
	    	    is_neumann = YES;
	    	    break;
	    	}
	    }
	}
	is_neumann = pp_max_status(is_neumann);
	if (is_neumann == NO)
	    return;

	cwnp = contact_wall_node_params(intfc);

/**
	if (io_type != NULL)
	{
	    set_use_normal_D_extend(cwnp->adjust);
	    return;
	}
**/

	screen("Insert wall normal at wall contact nodes (dflt = %s): ",
		(cwnp->adjust == YES) ? "yes" : "no");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if (s[0] == 'y' || s[0] == 'Y')
	    	cwnp->adjust = YES;
	    else if (s[0] == 'n' || s[0] == 'N')
	    	cwnp->adjust = NO;
	}
	set_use_normal_D_extend(cwnp->adjust);
	if (cwnp->adjust == NO)
	    return;

	screen("Enter first time step to begin adjustment (dflt = %d): ",
		cwnp->first_adjust_step);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&cwnp->first_adjust_step);
	screen("Enter the first real time to begin adjustment (dflt = %g): ",
		cwnp->first_adjust_time);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&cwnp->first_adjust_time);
	}
	screen("Enter the scaled wall bond length (dflt = %g): ",
		cwnp->wall_bond_len);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&cwnp->wall_bond_len);
	}
}		/*end prompt_for_contact_wall_node_params*/


LOCAL	void init_bow_shock(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front     *front = ip->root->front;
	double     angle1, angle2, thickness, dist_from_left;
	double     rad[MAXD];
	double     surf_ten;
	double     *L, *U;
	double     center[MAXD], theta, delta_theta;
	double     coords[MAXD];
	int	  type;
	int       *gmax;
	int       i, num_points;
	NODE      *ns, *ne;
	CURVE     *cur, *ramp;
	RECT_GRID *rect_grid = front->rect_grid;
	int       dim = rect_grid->dim;
	Gas_param *params;

	if (dim != 2)
	{
	    screen("ERROR in init_bow_shock(), dim = %d != 2 not supported\n",
		   dim);
	    clean_up(ERROR);
	}
	ramp = NULL;

	screen("Enter the ramp angles (degrees),\n");
	screen("\tthickness, and distance to inlet: ");
	(void) Scanf("%f %f %f %f\n",&angle1,&angle2,&thickness,
		     &dist_from_left);
	screen("Enter height of bowshock at exit: ");
	(void) Scanf("%f\n",&rad[1]);

	angle1 = radians(angle1);	angle2 = radians(angle2);

	type = prompt_for_wave_type("",front->interf,ip);
	params = init_eos_params(init,ip,"",YES);
	prompt_for_ambient_state(comp_type(COMPA),params," ahead",front,init);
	prompt_for_ambient_state(comp_type(COMPB),params," behind",front,init);

	surf_ten = prompt_for_surface_tension(type,"of the bow wave ");

	set_obstacle_comp_type(comp_type(COMPOBST),front);

	gmax = rect_grid->gmax;
	L = rect_grid->L;
	U = rect_grid->U;
	make_ramp(NORMAL,angle1,angle2,thickness,dist_from_left,
	    	  COMPOBST,COMPB,rect_grid,&ramp);
	center[0] = U[0];
	center[1] = L[1];
	rad[0] = U[0] - Coords(ramp->end->posn)[0];
	ns = ramp->end;
	node_type(ns) = ATTACHED_B_NODE;
	coords[0] = U[0];	coords[1] = rad[1] + L[1];
	ne = make_node(Point(coords));
	node_type(ne) = DIRICHLET_NODE;
	cur = make_curve(COMPA,COMPB,ns,ne);
	wave_type(cur) = type;
	start_status(cur) = INCIDENT;
	end_status(cur) = INCIDENT;
	surface_tension(cur) = surf_ten;
	num_points = 2*max(gmax[0],gmax[1]);
	delta_theta = .5*PI / num_points;
	for (i = 1; i <= num_points; i++)
	{
	    theta = i * delta_theta;
	    coords[0] = center[0] - rad[0] * cos(theta);
	    coords[1] = center[1] - rad[1] * sin(theta);
	    if (insert_point_in_bond(Point(coords),cur->last,cur) !=
		FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in init_bow_shock(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	}
}		/*end init_bow_shock*/

LOCAL	void init_expanding_shock(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	INTERFACE	*intfc = front->interf;
	double		*L,*U;
	double		x;
	double		angle1, angle2;
	double		thickness, dist_from_left;
	double		M;
	double		rho, pr;
	double		coords0[MAXD], coords1[MAXD];
	double		s_n[MAXD];
	char		s[Gets_BUF_SIZE];
	POINT		*p, Pcross;
	BOND		*b1, *b2;
	NODE		*ns, *ne;
	CURVE		*ramp, *incident;
	CURVE		**split_curves;
	RECT_GRID	*rect_grid = front->rect_grid;
	Gas_param	*params;
	boolean		sav_scss = interpolate_states_at_split_curve_node();
	int		constant_comp;

	ramp = NULL;
	L = rect_grid->L;
	U = rect_grid->U;
	s_n[0] = -1.;	s_n[1] = 0;

	screen("Enter the ramp angles (degrees),\n");
	screen("\tthickness, and distance to outlet: ");
	(void) Scanf("%f %f %f %f\n",&angle1,&angle2,&thickness,
		     &dist_from_left);
	screen("Enter the distance of the shock from the inlet: ");
	(void) Scanf("%f\n",&x);

	x = U[0] - x;
	angle1 = radians(angle1);	angle2 = radians(angle2);

	params = init_eos_params(init,ip,"",YES);
	screen("Enter the ahead state -- rho,pr: ");
	(void) Scanf("%f %f\n",&rho,&pr);
	screen("Enter the incident Mach number: ");
	(void) Scanf("%f\n",&M);

	set_ambient_comp_type(comp_type(COMPA),front);
	screen("Is the flow ahead of the wave constant (dflt = no): ");
	(void) Gets(s);
	constant_comp = (s[0] == 'y' || s[0] == 'Y') ? YES : NO;
	set_ambient_comp_type(comp_type(COMPB),front);
	set_obstacle_comp_type(comp_type(COMPOBST),front);
	init_shock_states(rho,pr,M,s_n,params,Ambient(comp_type(COMPA)),
			  Ambient(comp_type(COMPB)));
	if (constant_comp == YES)
	    (void)SetConstantFlowRegion(COMPA,Ambient(comp_type(COMPA)),intfc);
	make_ramp(NORMAL,angle1,angle2,thickness,dist_from_left,COMPOBST,COMPB,
		  rect_grid,&ramp);

	coords0[0] = x;	coords0[1] = L[1];
	coords1[0] = x;	coords1[1] = U[1];
	b1 = Bond(Point(coords0),Point(coords1));
	b2 = ramp->first;
	if (!cross_bonds(b1,b2,&Pcross))
	{
	    screen("ERROR in init_expanding_shock(), "
		   "the incident shock is not on ramp\n");
	    clean_up(ERROR);
	}
	p = Point(Coords(&Pcross));
	set_interpolate_states_at_split_curve_node(NO);
	split_curves = split_curve(p,ramp->first,ramp,COMPOBST,COMPB,COMPOBST,
				   COMPA);
	set_interpolate_states_at_split_curve_node(sav_scss);
	wave_type(split_curves[0]) = NEUMANN_BOUNDARY;
	start_status(split_curves[0]) = FIXED;
	end_status(split_curves[0]) = FIXED;
	wave_type(split_curves[1]) = NEUMANN_BOUNDARY;
	start_status(split_curves[1]) = FIXED;
	end_status(split_curves[1]) = FIXED;
	ns = split_curves[0]->end;
	node_type(ns) = NEUMANN_NODE;
	ne = make_node(Point(coords1));
	node_type(ne) = NEUMANN_NODE;
	incident = make_curve(COMPA,COMPB,ns,ne);
	wave_type(incident) = BACKWARD_SHOCK_WAVE;
	start_status(incident) = INCIDENT;
	end_status(incident) = INCIDENT;

}		/*end init_expanding_shock*/

/*
*			set_passive_curve_flags():
*
*	The wave_type of a boundary curve corresponds to the
*	boundary data type unless it borders on an OBSTACLE
*	region, in which case it is PASSIVE_BOUNDARY.  Here we
*	set the flags for the passive curves.
*/

LOCAL void set_passive_curve_flags(
	Front		*front)
{
	INTERFACE	*intfc = front->interf;

	switch (intfc->dim)
	{
	case 2:
	{
	    CURVE *c;
	    for ((void) next_curve(intfc,NULL); next_curve(intfc,&c); )
	    {
	    	if (!is_bdry(c))
	    	    continue;
	    	if (comp_type(positive_component(c))->type != OBSTACLE)
	    	    continue;
	    	wave_type(c) = PASSIVE_BOUNDARY;
	    	start_status(c) = PASSIVE;
	    	end_status(c) = PASSIVE;
	    	continue;
	    }
	}
	    break;
	case 3:
	    /* TODO: 3d */
	    break;
	}
}		/*end set_passive_curve_flags*/


LOCAL void g_bdry_state_initializer(
	int		idir,
	int		iside,
	HYPER_SURF	*hs,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	INTERFACE	*intfc = ip->root->front->interf;
	RECT_GRID	*gr = computational_grid(intfc);
	COMPONENT	comp = NO_COMP;
	double		coords[MAXD];
	int		i, dim = gr->dim;
	double		*L = gr->L, *U = gr->U;

	if (hs == NULL)
	{
	    if (rect_boundary_type(intfc,idir,iside) == DIRICHLET_BOUNDARY)
	    {
	    	for (i = 0; i < dim; i++)
	    	    coords[i] = 0.5*(L[i] + U[i]);
	    	coords[idir] = (iside == 0) ? L[idir] : U[idir];
		comp = nearest_interior_comp(YES,FIRST_DYNAMIC_COMPONENT,
	    				         coords,intfc);
	    }
	    (void) prompt_for_boundary_state(
	    		rect_boundary_type(intfc,idir,iside),NULL,
	    		coords,comp,2*idir+iside,hs,init,ip);
	}
	else
	{
	    switch (dim)
	    {
	    case 1:
	    	g_1d_bdry_state_initializer(idir,iside,Point_of_hs(hs),init,ip);
		break;
	    case 2:
	    	g_2d_bdry_state_initializer(idir,iside,Curve_of_hs(hs),init,ip);
		break;
	    case 3:
	    	g_3d_bdry_state_initializer(idir,iside,Surface_of_hs(hs),
					    init,ip);
	    	break;
	    }
	}
}		/*end g_bdry_state_initializer*/


EXPORT	int g_prompt_for_boundary_state(
	int		w_type,
	const char	*name,
	double		*coords,
	COMPONENT	comp,
	int		index,
	HYPER_SURF	*hs,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	INTERFACE	*intfc = ip->root->front->interf;
	BOUNDARY_STATE	Bstate;
	static	Locstate tmpst = NULL;
	char		s[Gets_BUF_SIZE];
	int		i, dim = intfc->dim;
	size_t		sizest = size_of_state(intfc);
	COMP_TYPE	*ctype;

	if (tmpst == NULL)
	    alloc_state(intfc,&tmpst,sizest);

	switch (w_type) 
	{
	case NEUMANN_BOUNDARY:
	case NO_SLIP_NEUMANN_BOUNDARY:
	case PASSIVE_BOUNDARY:
	case REFLECTION_BOUNDARY:
	case SUBDOMAIN_BOUNDARY:
	case MIXED_TYPE_BOUNDARY:
	    return -1;
	case DIRICHLET_BOUNDARY:
	    screen("Specify the %s boundary state"
	           "\n\ttime-independent boundary state set by "
	           "ambient state (A, default),"
	           "\n\tflow-through boundary conditions (FT),"
	           "\n\tflow-through boundary with constant pressure (FP), or"
	           "\n\tflow-through boundary with time dependent pressure (FD), or"
	           "\n\trandom velocity inlet (R),"
	           "\n\ttime dependent state read from file (TDBS)."
	           "\n\tpreset boundary state (P), or"
	           "\n\ttime-independent state specified by the user (U)."
	           "\nEnter choice here: ",(name == NULL) ? "Dirichlet" : name);
	    (void) Gets(s);
	    if (s[0] == '\0')
		(void) strcpy(s,"A");
	    if (strncasecmp(s,"FT",2) == 0)
	    {
		Bstate._boundary_state = NULL;
		Bstate._boundary_state_function = flow_through_boundary_state;
		Bstate._boundary_state_function_name =
				strdup("flow_through_boundary_state");
		Bstate._boundary_state_data = NULL;
		Bstate._fprint_boundary_state_data = 
				f_fprint_boundary_state_data;
		index = add_bstate_to_list(&Bstate,intfc,index);
	    }
	    else if (strncasecmp(s,"FP",2) == 0)
	    {
		Bstate._boundary_state = tmpst;
		Bstate._boundary_state_function =
		    constant_pressure_flow_through_boundary_state;
		Bstate._boundary_state_function_name =
			strdup("constant_pressure_flow_through_boundary_state");
		Bstate._boundary_state_data = NULL;
		Bstate._fprint_boundary_state_data = 
				f_fprint_boundary_state_data;
		set_ambient_boundary_state(tmpst,coords,comp,intfc,init);
		set_state(tmpst,TGAS_STATE,tmpst);
		screen("Enter rho, pr: ");
		(void) Scanf("%f %f\n",&Dens(tmpst),&Press(tmpst));
		Dens(tmpst)=density(tmpst);
		screen("Density = %"FFMT" for pressure = %"FFMT"\n"
		       ,Dens(tmpst),Press(tmpst));

		set_state(tmpst,GAS_STATE,tmpst);
		index = add_bstate_to_list(&Bstate,intfc,index);
	    }
	    else if (strncasecmp(s,"FD",2) == 0) /*reformat*/
	    {
/*  #bjet	     */
		double tw,tp,tc,prb,prp;
                /* Feb 19 2004: Myoung-Nyoun: Fix for general eos */
                double ref_dens, ref_pres;
                static Locstate ref_st = NULL;
                if (ref_st == NULL)
                    alloc_state(intfc,&ref_st,sizest);
		
		/** New code from Mnkim **/
		set_ambient_boundary_state(ref_st,coords,comp,intfc,init);
		zero_state_velocity(ref_st,dim);
                set_type_of_state(ref_st,TGAS_STATE);
                Bstate._boundary_state_function =
                time_dep_pressure_flow_through_boundary_state;
                Bstate._boundary_state_function_name =
                strdup("time_dep_pressure_flow_through_boundary_state");
                Bstate._boundary_state_data = NULL;
                Bstate._fprint_boundary_state_data =
                             g_fprint_tdp_boundary_state_data;
                set_ambient_boundary_state(tmpst,coords,comp,intfc,init);
                /* Feb 19 2004: Myoung-Nyoun: Exclude kinetic energy */
                zero_state_velocity(tmpst,dim);
		
                set_type_of_state(tmpst,TGAS_STATE);
		screen("Enter warming-up time, peak time, cooling-down time,"
                      "\n\tbase pressure, peak pressure, reference density and pressure: ");
                (void) Scanf("%f %f %f %f %f %f %f",&tw,&tp,&tc,&prb,&prp,
			     &ref_dens,&ref_pres);
                (void) getc(stdin);/*get trailing newline*/
                screen("\n");

                Dens(ref_st) = ref_dens;
                Press(ref_st) = ref_pres;
                (void) set_time_dep_pres(&Bstate,tw,tp,tc,prb,prp,ref_st);
		 
                /* Feb 19 2004: Myoung-Nyoun: Fix for general eos */
                state_on_adiabat_with_pr(ref_st,prb,tmpst,TGAS_STATE);
                /* 070103 : add : check density-pressure pairs */
                screen("Density = %"FFMT" for pressure = %"FFMT", comp = %d\n",
			  Dens(tmpst),Press(tmpst), comp);

	        /* Feb 19 2004: Myoung-Nyoun: Fix for general eos */
	        state_on_adiabat_with_pr(ref_st,prp,tmpst,TGAS_STATE);
	        screen("Density = %"FFMT" for pressure = %"FFMT", comp = %d\n\n",
		          Dens(tmpst),Press(tmpst), comp);

	 	/* Feb 19 2004: Myoung-Nyoun: Fix for general eos */
                state_on_adiabat_with_pr(ref_st,prb,tmpst,TGAS_STATE);
                set_state(tmpst,GAS_STATE,tmpst);
                Bstate._boundary_state = tmpst;
                index = add_bstate_to_list(&Bstate,intfc,index);
	    }
	    else if (strncasecmp(s,"U",1) == 0)
	    {
	    	(void) strcpy(s,"n");
	    	ctype = comp_type(comp);
	    	if (ctype->type == AMBIENT) 
	    	{
	    	    screen("Use adjacent ambient state params (y,n(deflt)?: ");
	    	    (void) Gets(s);
	    	}
	    	if (s[0] == 'y' || s[0] == 'Y')
		    Set_params(tmpst,Ambient(ctype));
		else
		{
		    Init_params(tmpst,init_eos_params(init,ip,"",YES));
#if defined(COMBUSTION_CODE)
		    prompt_for_burning(&Params(tmpst),"");
#endif /* defined(COMBUSTION_CODE) */
		}

		set_type_of_state(tmpst,TGAS_STATE);

		screen("Enter rho, pr, ");
#if defined(COMBUSTION_CODE)
		if (Composition_type(tmpst) == ZND)
		{
		    screen("react, ");
		}
#endif /* defined(COMBUSTION_CODE) */
		for (i = 0; i < dim-1; i++) screen("v[%d], ",i);
		    screen("v[%d]: ",i);

		(void) Scanf("%f %f",&Dens(tmpst),&Press(tmpst));
#if defined(COMBUSTION_CODE)
		if (Composition_type(tmpst) == ZND)
		{
		    (void) Scanf("%f",&React(tmpst));
		}
#endif /* defined(COMBUSTION_CODE) */
		for (i = 0; i < dim; i++)
		    (void) Scanf("%f",&Vel(tmpst)[i]);
		(void) getc(stdin);/*get trailing newline*/
		screen("\n");
		Dens(tmpst)=density(tmpst);
		screen("Density = %"FFMT" for pressure = %"FFMT"\n"
		       ,Dens(tmpst),Press(tmpst));

		set_state(tmpst,GAS_STATE,tmpst);
		Bstate._boundary_state = tmpst;
		Bstate._boundary_state_function = g_fixed_boundary_state;
		Bstate._boundary_state_function_name = 
		    strdup("g_fixed_boundary_state");
		Bstate._boundary_state_data = NULL;
		Bstate._fprint_boundary_state_data = 
				f_fprint_boundary_state_data;
		index = add_bstate_to_list(&Bstate,intfc,index);
		break;
	    }
	    else if (strncasecmp(s,"R",1) == 0)
	    {
	    	index = prompt_for_random_flow_inlet(coords,comp,hs,
						     intfc,index,init);
	    }
	    else if (strncasecmp(s,"TDBS",4) == 0)
            {
                COMP_TYPE *ctype = comp_type(comp);
	    	index = prompt_for_time_dependent_boundary_state(index,
                                 ctype->params,intfc);
            }
	    else if (strncasecmp(s,"A",1) == 0)
	    {
	        set_ambient_boundary_state(tmpst,coords,comp,intfc,init);
	    	Bstate._boundary_state_function = g_fixed_boundary_state;
		Bstate._boundary_state_function_name =
		    strdup("g_fixed_boundary_state");
		Bstate._boundary_state = tmpst;
		Bstate._boundary_state_data = NULL;
		Bstate._fprint_boundary_state_data = 
				f_fprint_boundary_state_data;
		index = add_bstate_to_list(&Bstate,intfc,index);
	    }
	    return index;
	case UNKNOWN_BOUNDARY_TYPE:
	    return index;
	default:
	    screen("ERROR in g_prompt_for_boundary_state(), "
	           "non-bdry wave type - %d\n",w_type);
	    clean_up(ERROR);
	    break;
	}
	return index;
}		/*end g_prompt_for_boundary_state*/

EXPORT	void	set_ambient_boundary_state(
	Locstate	bdry_state,
	double		*coords,
	COMPONENT	comp,
	INTERFACE	*intfc,
	INIT_DATA	*init)
{
	COMP_TYPE	*ctype;

	ctype = comp_type(comp);

	Get_state(coords,bdry_state,ctype,NULL,intfc,init,GAS_STATE);
}		/*end set_ambient_boundary_state*/


LOCAL void g_1d_bdry_state_initializer(
	int		idir,
	int		iside,
	POINT		*p,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	INTERFACE *intfc = ip->root->front->interf;
	double	  coords[MAXD];
	int       w_type;

	w_type = rect_boundary_type(intfc,idir,iside);
	if ((w_type != UNKNOWN_BOUNDARY_TYPE) &&
	    (w_type != MIXED_TYPE_BOUNDARY))
	    return;
	coords[0] = Coords(p)[0];
	bstate_index(p) = prompt_for_boundary_state(wave_type(p),NULL,coords,
						    positive_component(p),-1,
						    Hyper_surf(p),init,ip);
}		/*end g_1d_bdry_state_initializer*/


LOCAL void g_2d_bdry_state_initializer(
	int		idir,
	int		iside,
	CURVE		*c,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	INTERFACE	*intfc = ip->root->front->interf;
	double		coords[MAXD];
	int		i, dim = intfc->dim;
	int       w_type;

	start_status(c) = FIXED;	end_status(c) = FIXED;
	w_type = rect_boundary_type(intfc,idir,iside);
	if ((w_type!=UNKNOWN_BOUNDARY_TYPE) && (w_type!=MIXED_TYPE_BOUNDARY))
	    return;
	for (i = 0; i < dim; i++)
	{
	    coords[i] = 0.5*(Coords(c->start->posn)[i]+Coords(c->end->posn)[i]);
	}
	bstate_index(c) = prompt_for_boundary_state(wave_type(c),NULL,coords,
						    positive_component(c),-1,
						    Hyper_surf(c),init,ip);
}		/*end g_2d_bdry_state_initializer*/

LOCAL void g_3d_bdry_state_initializer(
	int		idir,
	int		iside,
	SURFACE		*s,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front     *front = ip->root->front;
	INTERFACE *intfc = front->interf;
	double	  coords[3];
	int       w_type;
	int       i;
	double     *h = front->rect_grid->h;

	w_type = rect_boundary_type(intfc,idir,iside);
	if ((w_type != UNKNOWN_BOUNDARY_TYPE) &&
	    (w_type != MIXED_TYPE_BOUNDARY))
	    return;
	average_position_of_surface(coords,s);
	for (i = 0; i < 3; i++)
	{
	    if (fabs(coords[i]) < 0.001*h[i]) /*TOLERANCE*/
	        coords[i] = 0.0;
	}
	bstate_index(s) = prompt_for_boundary_state(wave_type(s),NULL,coords,
						    positive_component(s),-1,
						    Hyper_surf(s),init,ip);
}		/*end g_3d_bdry_state_initializer*/

/*
*			g_init_cauchy_data_pointers():
*
*	Initializes function-pointer array handling the cauchy data.
*
*/

EXPORT void g_init_cauchy_data_pointers(
	INIT_PHYSICS	*ip,
	boolean		got_intfc_from_file)
{
	debug_print("init","Entered init_cauchy_data()\n");

	if (got_intfc_from_file == YES)
	{
	    ip->restart_initializer = g_restart_initializer;
	    ip->restart_intfc_initializer = g_restart_intfc_initializer;
	    ip->restart_pt_source_initializer = NULL;
	    ip->initializer = NULL;
	    ip->intfc_initializer = NULL;
	    ip->pt_source_initializer = NULL;
	}
	else
	{
	    ip->initializer = g_initializer;
	    ip->intfc_initializer = g_intfc_initializer;
	    ip->pt_source_initializer = NULL;
	    ip->restart_initializer = NULL;
	    ip->restart_intfc_initializer = NULL;
	    ip->pt_source_initializer = NULL;
	}
	debug_print("init","Left init_cauchy_data()\n");
}		/*end g_init_cauchy_data_pointers*/


/* ARGSUSED */
LOCAL void g_initializer(
	double		*coords,
	COMPONENT	comp,
	Locstate	state,
	INTERFACE	*intfc,
	INIT_DATA	*init)
{
	int i, dim = intfc->dim;
	int stype = type_of_state(init);

	debug_print("init_states","Entered g_initializer()\n");
	if (debugging("init_states"))
	{
  	    (void) printf("comp = %d, ",comp);
	    print_component_type("type = ",comp_type(comp),", ");
  	    for (i = 0; i < dim; i++)
  	    	(void) printf(" coords[%d] = %g ",i,coords[i]);
	    (void) printf("state - ");
	}

	Get_state(coords,state,comp_type(comp),NULL,intfc,init,stype);

	if (debugging("init_states"))
	{
	    print_gas_state(state);
	}
	debug_print("init_states","Left g_initializer()\n");
}		/*end g_initializer*/


/* ARGSUSED */
LOCAL void g_intfc_initializer(
	POINT		   *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	   *hs,
	Locstate	   lstate,
	Locstate	   rstate,
	INIT_DATA	   *init)
{
	COMPONENT      pcomp = NO_COMP, ncomp = NO_COMP;
	COMP_TYPE      *nct, *pct;
	INTERFACE      *intfc = hs->interface;
	int            stype = type_of_state(init);
	int	       typel, typer;
	int	       dim = intfc->dim;
	int		debug_flag = NO;

	debug_print("init_states","Entered g_intfc_initializer()\n");

	if (debugging("init_states"))
	{
	    (void) printf("hs = %llu, hse = %llu ",
	    	          hypersurface_number(hs),
			  hypersurface_element_number(hse,intfc));
	    print_general_vector("Coords(p) = ",Coords(p),dim,"\n");
	}

	ncomp = negative_component(hs);
	nct = comp_type(ncomp);
	pcomp = positive_component(hs);
	pct = comp_type(pcomp);
	switch(wave_type(hs))
	{
	case FORWARD_SOUND_WAVE_LE:
	case BACKWARD_SOUND_WAVE_TE:
	    /* Ambient side == right */
	    Get_state(Coords(p),rstate,pct,hs,intfc,init,stype);
	    copy_state(lstate,rstate);
	    goto leave;

	case FORWARD_SOUND_WAVE_TE:
	case BACKWARD_SOUND_WAVE_LE:
	    /* Ambient side == left */
	    Get_state(Coords(p),lstate,nct,hs,intfc,init,stype);
	    copy_state(rstate,lstate);
	    goto leave;

	case REFLECTION_BOUNDARY:
	case SUBDOMAIN_BOUNDARY:
	case PASSIVE_BOUNDARY:
	    obstacle_state(intfc,lstate,size_of_state(intfc));
	    obstacle_state(intfc,rstate,size_of_state(intfc));
	    goto leave;

	default:
	    break;
	}

	typel = nct->type;
	if (typel != EXTERIOR)
	{
	    switch(typel)
	    {
	    case PRANDTL_MEYER_WAVE:
	    case UNTRACKED_SHOCK:
	    	if (dim == 2)
	    	{
	    	    CURVE *c = Curve_of_hs(hs);
	            BOND *b = Bond_of_hse(hse);
		    if (p == c->start->posn)
			Get_state(Coords(b->end),lstate,nct,hs,
				  intfc,init,stype);
		    else if (p == c->end->posn)
			Get_state(Coords(b->start),lstate,nct,hs,
				  intfc,init,stype);
		    else
			Get_state(Coords(p),lstate,nct,hs,intfc,init,stype);
		}
		else
		    Get_state(Coords(p),lstate,nct,hs,
			      intfc,init,stype);
		break;
	    case BUBBLE:
	    	debug_print("init_states","BUBBLE\n");
	    	bubble_state(lstate,p,Curve_of_hs(hs),
	    		     (Bubble *)nct->extra,ncomp);
		break;
	    default:
		Get_state(Coords(p),lstate,nct,hs,intfc,init,stype);
	    }
	}

	typer = pct->type;
	if (typer != EXTERIOR)
	{
	    switch(typer)
	    {
	    case PRANDTL_MEYER_WAVE:
	    case UNTRACKED_SHOCK:
	    	if (dim == 2)
	    	{
	    	    CURVE *c = Curve_of_hs(hs);
		    BOND *b = Bond_of_hse(hse);
		    if (p == c->start->posn)
		        Get_state(Coords(b->end),rstate,pct,hs,
				  intfc,init,stype);
		    else if (p == c->end->posn)
			Get_state(Coords(b->start),rstate,pct,hs,
				  intfc,init,stype);
		    else
			Get_state(Coords(p),rstate,pct,hs,intfc,init,stype);
		}
		else
		    Get_state(Coords(p),rstate,pct,hs,intfc,init,stype);
		break;
	    case BUBBLE:
	    	debug_print("init_states","BUBBLE\n");
	    	bubble_state(rstate,p,Curve_of_hs(hs),
	    		     (Bubble *)pct->extra,pcomp);
		break;
	    default:
		Get_state(Coords(p),rstate,pct,hs,intfc,init,stype);
	    }
	}

	if ((typel == EXTERIOR) || (typel == OBSTACLE))
	{
	    switch (wave_type(hs)) 
	    {
	    case DIRICHLET_BOUNDARY:
	    case NEUMANN_BOUNDARY:
	    case MOVABLE_BODY_BOUNDARY:
	    case PASSIVE_BOUNDARY:
	    case REFLECTION_BOUNDARY:
	    case SUBDOMAIN_BOUNDARY:
	        obstacle_state(intfc,lstate,size_of_state(intfc));
	        break;

	    default:
	        screen("ERROR in g_intfc_initializer(), "
	               "unknown bdry type\n");
	        print_wave_type("typel == EXTERIOR, wave_type = ",
	        	wave_type(hs),"\n",intfc);
	        print_comp_types(intfc);
	        print_hypersurface(hs);
	        print_interface(intfc);
	        clean_up(ERROR);
	    }
	}

	if ((typer == EXTERIOR) || (typer == OBSTACLE))
	{
	    switch (wave_type(hs)) 
	    {
	    case DIRICHLET_BOUNDARY:
	    case NEUMANN_BOUNDARY:
	    case MOVABLE_BODY_BOUNDARY:
	    case PASSIVE_BOUNDARY:
	    case REFLECTION_BOUNDARY:
	    case SUBDOMAIN_BOUNDARY:
	    	obstacle_state(intfc,rstate,size_of_state(intfc));
	    	break;

	    default:
	    	screen("ERROR in g_intfc_initializer(), "
	    	       "unknown bdry type\n");
	    	print_wave_type("typer == EXTERIOR, wave_type = ",
	    		wave_type(hs),"\n",intfc);
	    	print_comp_types(intfc);
	    	print_hypersurface(hs);
	    	print_interface(intfc);
	    	clean_up(ERROR);
	    }
	}

leave:
	if (debugging("init_states"))
	{
	    print_general_vector("negative side state at ",Coords(p),
				 intfc->dim,"");
	    (void) printf(" on hs %llu\n",hypersurface_number(hs));
	    print_gas_state(lstate);
	    print_general_vector("positive side state at ",Coords(p),
			         intfc->dim,"");
	    (void) printf(" on hs %llu\n",hypersurface_number(hs));
	    print_gas_state(rstate);
	}
	debug_print("init_states","Left g_intfc_initializer()\n");
}		/*end g_intfc_initializer*/


LOCAL	void	print_comp_types(
	INTERFACE	*intfc)
{
	HYPER_SURF	**hs;
	COMPONENT	*comps;
	COMPONENT	comp;
	int		len = max_component(intfc) - min_component(intfc) + 1;
	int		i, num_comps;

	uni_array(&comps,len,sizeof(COMPONENT));
	for (num_comps = 0, hs = intfc->hss; hs && *hs; hs++)
	{
	    comp = positive_component(*hs);
	    for (i = 0; i < num_comps; i++)
	    	if (comp == comps[i]) break;
	    if (i == num_comps)
	    	comps[num_comps++] = comp;
	    comp = negative_component(*hs);
	    for (i = 0; i < num_comps; i++)
	    	if (comp == comps[i]) break;
	    if (i == num_comps)
	    	comps[num_comps++] = comp;
	}
	for (i = 0; i < num_comps; i++)
	{
	    (void) printf("Type of component %d = ",comps[i]);
	    print_component_type("",comp_type(comps[i]),"\n");
	}
	free(comps);
}		/*end print_comp_types*/

EXPORT	const char *comp_type_name(
	COMP_TYPE_TYPE type)
{
	static char s[512];
	switch (type)
	{
	case EXTERIOR:
	    return "EXTERIOR";
	case OBSTACLE:
	    return "OBSTACLE";
	case AMBIENT:
	    return "AMBIENT";
	case ELLIPTICAL:
	    return "ELLIPTICAL";
	case KH_SINE_PERTURBED:
	    return "KH_SINE_PERTURBED";
	case NWAVE:
	    return "NWAVE";
	case RANDOM_SURFACE_PERTURBED:
	    return "RANDOM_SURFACE_PERTURBED";
	case BUBBLE:
	    return "BUBBLE";
	case PRANDTL_MEYER_WAVE:
	    return "PRANDTL_MEYER_WAVE";
	case UNTRACKED_SHOCK:
	    return "UNTRACKED_SHOCK";
	case RESTART:
	    return "RESTART";
	case TAYLOR_WAVE:
	    return "TAYLOR_WAVE";
	case RT_KH:
	    return "RT_KH";
	case RT_PERTURBED:
	    return "RT_PERTURBED";
	case KH_PERTURBED:
	    return "KH_PERTURBED";
	case RAREFACTION_WAVE_1D:
	    return "RAREFACTION_WAVE_1D";
	case TRANS_LAYER:
	    return "TRANS_LAYER";
	case STRETCHING:
	    return "STRETCHING";
	case ONE_DIMENSIONAL_OVERLAY:
	    return "ONE_DIMENSIONAL_OVERLAY";
	default:
	    (void) sprintf(s,"UNKNOWN COMPONENT TYPE %d",type);
	    return s;
	}
}		/*end comp_type_name*/

LOCAL	void	print_component_type(
	const char *mesg1,
	COMP_TYPE  *comp_type,
	const char *mesg2)
{
	if (comp_type == NULL)
	    return;
	(void) printf("%s%s%s",mesg1,comp_type_name(comp_type->type),mesg2);
}		/*end print_component_type*/

LOCAL   void init_ambient(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Gas_param	*params;

	params = init_eos_params(init,ip,"",YES);

	prompt_for_ambient_state(comp_type(FIRST_DYNAMIC_COMPONENT),params,
			         " interior",ip->root->front,init);
}		/*end init_ambient*/

LOCAL	void init_2d_riemann(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Gas_param	*params;


	params = init_eos_params(init,ip,"",YES);

	set_2d_riemann_comp_type(comp_type(FIRST_DYNAMIC_COMPONENT),params,
				init);

}	/* end init_2d_riemann */

/*ARGSUSED*/
LOCAL	void	init_tri_grid_test(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	COMPONENT	comp;
	Front		*front = ip->root->front;
	INTERFACE	*intfc = front->interf;

	for (comp = min_component(intfc); comp <= max_component(intfc); comp++)
	{
	    set_obstacle_comp_type(comp_type(comp),front);
	}
}		/*end init_tri_grid_test*/

/*ARGSUSED*/
EXPORT	void	g_prompt_for_ref_state(
	const char *message,
	Locstate   state,
	int	   st_type,
	Gas_param  *params,
	INIT_DATA  *init)
{
	int	   i, dim = params->dim;
	static const char vname[3][3] = {"vx","vy","vz"};

	Init_params(state,params);
	set_type_of_state(state,TGAS_STATE);
#if defined(COMBUSTION_CODE)
	if (params->composition_type == ZND) 
	{
	    screen("Enter the gas state%s\n",message);
	    screen("\t\trho, pr, reaction progress");
	    for (i = 0; i < dim; i++)
		screen(", %s",vname[i]);
	    screen(": ");
	    (void) Scanf("%f %f %f",&Dens(state),&Press(state),
	    	     &React(state));
	    for (i = 0; i < dim; i++)
		(void) Scanf("%f",&Vel(state)[i]);
	    (void) Scanf("\n");
	}
	else
	{
	    screen("Enter the gas state%s\n\t\trho, pr",message);
	    for (i = 0; i < dim; i++)
		screen(", %s",vname[i]);
	    screen(": ");
	    (void) Scanf("%lf %lf ",&Dens(state),&Press(state));
	    for (i = 0; i < dim; i++)
		(void) Scanf("%lf",&Vel(state)[i]);
	    (void) Scanf("\n");
	    prompt_for_burning(&Params(state),message);
	}
#else /* defined(COMBUSTION_CODE) */
	screen("Enter the gas state%s\n\t\trho, pr",message);
	for (i = 0; i < dim; i++)
	    screen(", %s",vname[i]);
	screen(": ");
	(void) Scanf("%lf %lf ",&Dens(state),&Press(state));
	for (i = 0; i < dim; i++)
	    (void) Scanf("%lf ",&Vel(state)[i]);
	(void) Scanf("\n");
#endif /* defined(COMBUSTION_CODE) */

	Dens(state)=density(state); /* set correct density if fuel */
	screen("Density = %"FFMT" for pressure = %"FFMT"\n"
	       ,Dens(state),Press(state));

	if (st_type != state_type(state))
	    set_state(state,st_type,state);
}		/*end g_prompt_for_ref_state*/


EXPORT void prompt_for_ambient_state(
	COMP_TYPE	*comp_type,
	Gas_param	*params,
	const char	*message,
	Front		*front,
	INIT_DATA	*init)
{
 	size_t		sizest = params->sizest;
	Locstate	state,tstate;


	(*params->_alloc_state)(&tstate,sizest);
	prompt_for_ref_state(message,tstate,TGAS_STATE,params,init);

	set_ambient_comp_type(comp_type,front);
	state = Ambient(comp_type);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            Gas_param *params = Params(tstate);
            if(params->n_comps != 1)
            {
                if(debugging("injector_jet"))
                {
                    /* SF6 is the first comp, and is not mixed at init. */
                    pdens(tstate)[0] = Dens(tstate);
                    pdens(tstate)[1] = 0.0;
                    set_state(state,GAS_STATE,tstate);
                }
                else
                    set_state(state,TGAS_STATE,tstate);
            }
            else
                set_state(state,GAS_STATE,tstate);
        }
        else
	    set_state(state,GAS_STATE,tstate);

	if (debugging("init")) 
	    verbose_print_state("ambient",state);
	free(tstate);
}		/*end prompt_for_ambient_state*/

EXPORT	void	set_ambient_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	if (comp_type->type == AMBIENT)	/*ALREADY DONE*/
	    return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = AMBIENT;
	if ((front != NULL) && (comp_type->extra == NULL))
	{
	    Locstate	state;
	    if (front->sizest != 0)
	    	alloc_state(front->interf,&state,front->sizest);
	    comp_type->extra = (POINTER)state;
	}
	comp_type->_get_state = get_state_ambient;
	comp_type->free_comp_type_extra = free_ambient_comp_type;
}		/*end set_ambient_comp_type*/

LOCAL	void	free_ambient_comp_type(
	COMP_TYPE	*comp_type)
{
	if (comp_type->type != AMBIENT)
	    return;
	if (comp_type->extra != NULL)
	    free(comp_type->extra);
	comp_type->extra = NULL;
}		/*end free_ambient_comp_type*/

LOCAL	void	set_exterior_comp_type(
	COMP_TYPE	*comp_type)
{
	if (comp_type->type == EXTERIOR)
	    return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = EXTERIOR;
	comp_type->extra = NULL;
	comp_type->_get_state = get_state_exterior;
	comp_type->free_comp_type_extra = NULL;
}		/*end set_exterior_comp_type*/

/*ARGSUSED*/
LOCAL	void	get_state_exterior(
	double           *coords,
	Locstate        s,
	COMP_TYPE       *ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	debug_print("init_states","Entered get_state_exterior()\n");
	screen("ERROR in get_state_exterior(),  EXTERIOR comp found\n");
	print_general_vector("coords = ",coords,intfc->dim,"\n");
	(void) printf("comp = %d\n",ct->comp);
	print_interface(intfc);
	clean_up(ERROR);
}		/*end get_state_exterior*/

EXPORT	void	set_obstacle_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	if (comp_type->type == OBSTACLE)
	    return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = OBSTACLE;
	comp_type->extra = NULL;
	comp_type->_get_state = get_state_obstacle;
	comp_type->free_comp_type_extra = NULL;
	if (front != NULL)
	{
	    COMPONENT	comp = comp_type->comp;
	    exclude_comp(comp,front->interf);
	    (void)SetConstantFlowRegion(comp,return_obst_state(),front->interf);
	}
}		/*end set_obstacle_comp_type*/

/*ARGSUSED*/
LOCAL	void	get_state_obstacle(
	double           *coords,
	Locstate        s,
	COMP_TYPE       *ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int             stype)
{
	debug_print("init_states","Entered get_state_obstacle()\n");
	g_obstacle_state(s,g_sizest());
}		/*end get_state_obstacle*/


/* 
*			prompt_for_cylindrical_pencil():
*
*	In cylindrical geometry, a "pencil" near the axis is used.
*	To prevent curves from being constructed in the pencil,
*	RL_eff is set to the outer edge of the pencil.
*	TODO: 
*	A better solution would be to create curves in the whole region
*	and have make_vertical_axis_curve() delete all curves in the 
*	pencil.
*/

LOCAL	void	prompt_for_cylindrical_pencil(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	RECT_GRID *gr = ip->root->front->rect_grid;
	double     sfactor;
	char	  s[Gets_BUF_SIZE];

	RL_eff(init) = gr->GL[0];
	if (coord_system() != CYLINDRICAL)
	    return;

	if (RL_eff(init) != 0.0) /*User requests positive lower boundary*/
	    return;

        screen("\nRequest Neumann Wall to cut off zero radius  [y,n(dflt)]: ");
        (void) Gets(s);
        if(s[0] != 'y' && s[0] != 'Y')
            return;

        sfactor = 1.002;/*TOLERANCE*/
        screen("Enter a nonzero, positive, non-integer "
               "amount for the number\n");
        screen("of grid-blocks the wall is away "
               "from the center (default=%g): ",sfactor);

        (void) Gets(s);
        if (s[0] != '\0')
        {
	    double rtol = 0.000001;/*TOLERANCE*/
            (void) sscan_float(s, &sfactor);
            if(sfactor <= 0.0)
            {
                screen("ERROR in prompt_for_cylindrical_pencil(), "
                       "The factor chosen %g is not positive\n",sfactor);
                clean_up(ERROR);
            }
            else if (sfactor-floor(sfactor) < rtol)/*TOLERANCE*/
            {
                (void) printf("WARNING in prompt_for_cylindrical_pencil(), "
                              "The factor chosen %g is almost integer\n",
                              sfactor);
	        sfactor += rtol;
            }
        }

        RL_eff(init) = sfactor*cell_width(0,0,gr);
	RL_eff(init) = pos_radius(RL_eff(init),gr);

        if (RL_eff(init) < gr->GL[0])
	    RL_eff(init) = gr->GL[0];
}		/*end prompt_for_cylindrical_pencil*/



LOCAL   void init_oned_density_step(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	DENSITY_STEP_WALL_DATA   *wdata;
	Grid                     *grid = ip->root->grid;
	Front                    *front = ip->root->front;
	Printplot                *prt = ip->prt;
	INTERFACE                *intfc = front->interf;
	RECT_GRID                *gr = front->rect_grid;
	COMP_TYPE                *ctype[3];
	POINT                    *pc, *ps;
	Locstate                 initial_state[3];
	Locstate                 state[7];
	Locstate                 left, right;
	boolean                     isforward;
	char                     buf[Gets_BUF_SIZE];
	double                    pml, pmr, uml, umr, ml, mr;
	double                    c, s;
	double                    L, U, sL, sU;
	double                    W[9];
	double                    DX, dx;
	double                    center[3], r[3];
	double                    sgn;
	double                    dWall;
	int                      rside;
	int                      n[3];
	int                      N;
	int                      i;
	int                      dim = gr->dim;
	RIEMANN_SOLVER_WAVE_TYPE rswtype[9];
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	int                      w_type;
	size_t                   sizest = front->sizest;
	static double             nor[] = {1.0, 0.0, 0.0};

	if (dim != 1)
	{
	    screen("ERROR in init_oned_density_step(), "
	           "dimension %d not supported\n",dim);
	    clean_up(ERROR);
	}

	ctype[0] = comp_type(FIRST_DYNAMIC_COMPONENT);
	ctype[1] = comp_type(FIRST_DYNAMIC_COMPONENT+1);
	ctype[2] = comp_type(FIRST_DYNAMIC_COMPONENT+2);
	for (i = 0; i < 3; ++i)
	{
	    set_ambient_comp_type(ctype[i],front);
	    initial_state[i] = Ambient(ctype[i]);
	}

	(void) prompt_for_eos_params(init,ip,YES,"");
	ctype[0]->params = prompt_for_eos_params(init,ip,YES,
	                                         "of the target");
	isforward = NO;
	screen("Enter propagation direction of the incident shock "
	       "(dflt = %s): ",(isforward)?"left to right":"right to left");
	(void) Gets(buf);
	if ((buf[0] == 'R') || (buf[0] == 'r'))
	    isforward = NO;
	else if ((buf[0] == 'L') || (buf[0] == 'l'))
	    isforward = YES;
	if (isforward)
	{
	    sgn = 1.0;
	    rside = 1;
	}
	else
	{
	    sgn = -1.0;
	    rside = 0;
	}
	screen("Enter the density and pressure of the target: ");
	Init_params(initial_state[0],ctype[0]->params);
	set_type_of_state(initial_state[0],TGAS_STATE);
	(void) Scanf("%f %f\n",
	             &Dens(initial_state[0]),&Press(initial_state[0]));
	zero_state_velocity(initial_state[0],dim);

	ctype[1]->params = prompt_for_eos_params(init,ip,YES,
	                                         "ahead of the incident shock");
	prompt_for_behind_contact_state(initial_state[0],initial_state[1],
	                                ctype[1]->params,NO,nor,GAS_STATE,init);

	ctype[2]->params = ctype[1]->params;
	prompt_for_behind_shock_state(initial_state[1],initial_state[2],
				      NO,nor,GAS_STATE,isforward,init);
	W[0] = 0.0;
	W[1] = sgn*mass_flux(pressure(initial_state[2]),initial_state[1])/
	           Dens(initial_state[1]);

	L = gr->L[0];
	U = gr->U[0];
	screen("Enter the initial contact position c (%g < c < %g): ",L,U);
	(void) Scanf("%f\n",&c);
	if ((c <= L) || (U <= c))
	{
	    screen("ERROR in init_oned_density_step(), "
	           "initial contact position %g out of range\n",c);
	    clean_up(ERROR);
	}
	if (isforward)
	{
	    sL = L;
	    sU = c;
	    dWall = c - U;
	}
	else
	{
	    sL = c;
	    sU = U;
	    dWall = c - L;
	}
	screen("Enter the initial shock position s (%g < s < %g): ",sL,sU);
	(void) Scanf("%f\n",&s);
	if ((s <= sL) || (sU <= s))
	{
	    screen("ERROR in init_oned_density_step(), "
	           "initial shock position %g out of range\n",s);
	    clean_up(ERROR);
	}

	rect_boundary_type(intfc,0,rside) = REFLECTION_BOUNDARY;
	screen("Use Neumman (N) or reflecting (R) "
	       "boundaries on the %s (dflt = %s): ",(isforward)?"right":"left",
	       (rect_boundary_type(intfc,0,1) == REFLECTION_BOUNDARY) ?
	       "R" : "N");
	(void) Gets(buf);
	if ((buf[0] == 'R') || (buf[0] == 'r'))
	    rect_boundary_type(intfc,0,rside) = REFLECTION_BOUNDARY;
	else if ((buf[0] == 'N') || (buf[0] == 'n'))
	    rect_boundary_type(intfc,0,rside) = NEUMANN_BOUNDARY;

	pc = (isforward) ? make_point(&c,ctype[1]->comp,ctype[0]->comp) :
	                   make_point(&c,ctype[0]->comp,ctype[1]->comp);
	wave_type(pc) = CONTACT;
	if (ctype[0]->params == ctype[1]->params)
	{
	    screen("Type y to track the contact: ");
	    (void) Gets(buf);
	    if ((buf[0] != 'y') && (buf[0] != 'Y'))
	        untracked_hyper_surf(pc) = YES;
	}
	ps = (isforward) ? make_point(&s,ctype[2]->comp,ctype[1]->comp) :
	                   make_point(&s,ctype[1]->comp,ctype[2]->comp);
	wave_type(ps) = (isforward) ? FORWARD_SHOCK_WAVE : BACKWARD_SHOCK_WAVE;
	screen("Type y to track the incident shock: ");
	(void) Gets(buf);
	if ((buf[0] != 'y') && (buf[0] != 'Y'))
	    untracked_hyper_surf(ps) = YES;

	DX = 3.0*gr->h[0];
	N = 5;
	screen("Enter an averaging length for the "
	       "pressure near the wall (dflt = %g): ",DX);
	(void) Gets(buf);
	if (buf[0] != '\0')
	{
	    if (sscan_float(buf,&DX) != 1)
	    {
	        screen("ERROR in init_oned_density_step(), improper input "
	               "of pressure averaging length\n");
	        clean_up(ERROR);
	    }
	}
	screen("Enter the number of pressure samples per averaging length "
	       "(dflt = %d): ",N);
	(void) Gets(buf);
	if (buf[0] != '\0')
	{
	    if (sscanf(buf,"%d",&N) != 1)
	    {
	        screen("ERROR in init_oned_density_step(), improper input "
	               "of number of pressure samples per averaging length\n");
	        clean_up(ERROR);
	    }
	}
	for (i = 0; i < 7; ++i)
	    alloc_state(intfc,state+i,max(sizest,sizeof(VGas)));
	wdata = alloc_DENSITY_STEP_WALL_DATA(2,N,front);
	wdata->isforward = isforward;
	set_state(state[0],VGAS_STATE,initial_state[0]);
	set_state(state[1],VGAS_STATE,initial_state[1]);
	set_state(state[2],VGAS_STATE,initial_state[2]);
	set_state(wdata->exact_wstate[0],VGAS_STATE,state[0]);
	wdata->Plast = pressure(wdata->exact_wstate[0]);
	wdata->c = c;
	wdata->s = s;
	wdata->t[0] = grid->time;
	wdata->DX = DX;
	dx = DX/N;
	if (isforward)
	{
	    for (i = 0; i < N; ++i)
	        wdata->px[i] = U - (i+0.5)*dx;
	}
	else
	{
	    for (i = 0; i < N; ++i)
	        wdata->px[i] = L + (i+0.5)*dx;
	}

	alloc_state(intfc,&left,max(sizest,sizeof(VGas)));
	alloc_state(intfc,&right,max(sizest,sizeof(VGas)));

	/* Incident shock collides with the contact */
	if (isforward)
	{
	    if (!find_mid_state(state[2],state[0],0.0,&pml,&pmr,&uml,
	                        &umr,&ml,&mr,rswtype+4,rswtype+2))
	    {
	        screen("ERROR in init_oned_density_step(), find_mid_state() "
	                      "did not converge\n");
	        clean_up(ERROR);
	    }
	    state_behind_sound_wave(state[0],state[3],NULL,W+2,0.0,mr,umr,
	                            pmr,VGAS_STATE,FORWARD_SHOCK_WAVE,
				    rswtype[2],RIGHT_FAMILY);
	    w_type = (rswtype[4] == SHOCK) ? BACKWARD_SHOCK_WAVE :
	                                     BACKWARD_SOUND_WAVE_TE;
	    state_behind_sound_wave(state[2],state[4],NULL,W+4,0.0,ml,uml,
	                            pml,VGAS_STATE,w_type,rswtype[4],
				    LEFT_FAMILY);
	}
	else
	{
	    set_state_for_find_mid_state(state[2],state[2]);
	    if (!find_mid_state(state[0],state[2],0.0,&pml,&pmr,&uml,&umr,
	                        &ml,&mr,rswtype+2,rswtype+4))
	    {
	        screen("ERROR in init_oned_density_step(), find_mid_state() "
	                      "did not converge\n");
	        clean_up(ERROR);
	    }
	    state_behind_sound_wave(state[0],state[3],NULL,W+2,0.0,ml,uml,
	                            pml,VGAS_STATE,BACKWARD_SHOCK_WAVE,
				    rswtype[2],LEFT_FAMILY);
	    w_type = (rswtype[4] == SHOCK) ? FORWARD_SHOCK_WAVE :
	                                     FORWARD_SOUND_WAVE_TE;
	    state_behind_sound_wave(state[2],state[4],NULL,W+4,0.0,mr,umr,
	                            pmr,VGAS_STATE,w_type,rswtype[4],
				    RIGHT_FAMILY);
	}
	if (rswtype[4] == RAREFACTION)
	    W[4] = vel(0,state[2]) - sgn*sound_speed(state[2]);
	wdata->t[1] = wdata->t[0] + (c - s)/W[1] - dWall/W[2];
	wdata->exact_incomingWs[0] = W[2];
	W[3] = 0.5*(uml+umr);
	wdata->wave_tol = ((pressure(state[3])-pressure(state[0]))/
	                        (pressure(state[3])+pressure(state[0])))*
			   fabs(W[2])/
			   (sound_speed(state[3])+sound_speed(state[0]));
	screen("Enter an error tolerance to identify the pressure wave "
	        "profile tails (dflt = %g): ",wdata->wave_tol);
	    (void) Gets(buf);
	if (buf[0] != '\0')
	{
	    if (sscan_float(buf,&wdata->wave_tol) != 1)
	    {
	        screen("ERROR in init_oned_density_step(), improper input "
	               "of pressure rise tolerance\n");
	        clean_up(ERROR);
	    }
	}

	r[0] = 0.5*DX;         r[1] = 0.5;      r[2] = 0.5;
	center[0] = (isforward) ? center[0] = U - r[0] : L + r[0];
	center[1] = 0.0; center[2] = 0.0;
	n[0] = max((N-1)/2,0); n[1] = 0;        n[2] = 0;
	(void) init_probe(center,r,n,&wdata->wave_tol,init,ip);


	/* Transmitted shock reaches the wall */
	if (isforward)
	{
	    set_state_for_find_mid_state(right,state[3]);
	    Vel(right)[0] *= -1.0;
	    if (!find_mid_state(state[3],right,0.0,&pml,&pmr,&uml,&umr,&ml,&mr,
	                        rswtype+5,&r_wave))
	    {
	        screen("ERROR in init_oned_density_step(), find_mid_state() "
	                      "did not converge\n");
	        clean_up(ERROR);
	    }
	    state_behind_sound_wave(state[3],wdata->exact_wstate[1],NULL,
	                            W+5,0.0,ml,uml,pml,VGAS_STATE,
				    BACKWARD_SHOCK_WAVE,rswtype[5],LEFT_FAMILY);
	    wdata->exact_outgoingWs[0] = W[5];

	    /* Reflected shock from wall collides with contact */
	    set_state_for_find_mid_state(left,state[4]);
	    set_state_for_find_mid_state(right,wdata->exact_wstate[1]);
	    if (!find_mid_state(left,right,0.0,&pml,&pmr,&uml,&umr,&ml,&mr,
	                        rswtype+8,rswtype+6))
	    {
	        screen("ERROR in init_oned_density_step(), find_mid_state() "
	                      "did not converge\n");
	        clean_up(ERROR);
	    }
	    w_type = (rswtype[6]==SHOCK) ? FORWARD_SHOCK_WAVE :
	                                   FORWARD_SOUND_WAVE_LE;
	    state_behind_sound_wave(right,state[5],NULL,W+6,0.0,mr,umr,
	                            pmr,VGAS_STATE,w_type,rswtype[6],
				    RIGHT_FAMILY);
	    state_behind_sound_wave(left,state[6],NULL,W+8,0.0,ml,uml,
	                            pml,VGAS_STATE,w_type,rswtype[8],
				    LEFT_FAMILY);
	}
	else
	{
	    set_state_for_find_mid_state(left,state[3]);
	    Vel(left)[0] *= -1.0;
	    if (!find_mid_state(left,state[3],0.0,&pml,&pmr,&uml,&umr,&ml,&mr,
	                        &l_wave,rswtype+5))
	    {
	        screen("ERROR in init_oned_density_step(), find_mid_state() "
	                      "did not converge\n");
	        clean_up(ERROR);
	    }
	    state_behind_sound_wave(state[3],wdata->exact_wstate[1],NULL,
	                            W+5,0.0,mr,umr,pmr,VGAS_STATE,
				    FORWARD_SHOCK_WAVE,rswtype[5],RIGHT_FAMILY);
	    wdata->exact_outgoingWs[0] = W[5];

	    /* Reflected shock from wall collides with contact */
	    set_state_for_find_mid_state(left,wdata->exact_wstate[1]);
	    set_state_for_find_mid_state(right,state[4]);
	    if (!find_mid_state(left,right,0.0,&pml,&pmr,&uml,&umr,&ml,&mr,
	                           rswtype+6,rswtype+8))
	    {
	        screen("ERROR in init_oned_density_step(), find_mid_state() "
	                      "did not converge\n");
	        clean_up(ERROR);
	    }
	    w_type = (rswtype[6]==SHOCK) ? BACKWARD_SHOCK_WAVE :
	                                   BACKWARD_SOUND_WAVE_LE;
	    state_behind_sound_wave(left,state[5],NULL,W+6,0.0,ml,uml,
	                            pml,VGAS_STATE,w_type,rswtype[6],
				    LEFT_FAMILY);
	    state_behind_sound_wave(right,state[6],NULL,W+8,0.0,mr,umr,
	                            pmr,VGAS_STATE,w_type,rswtype[8],
				    RIGHT_FAMILY);
	}
	W[7] = 0.5*(uml+umr);
	wdata->t[2] = wdata->t[1] +
	             dWall*(1.0-W[3]/W[2])*(1.0-W[5]/W[6])/(W[5]-W[3]);

	/*stop_time(grid) = wdata->t[2];*/

	init_output_data(init,&wdata->odata,grid,prt,"the wall data",YES,NO,NO);
	add_user_output_function(print_density_step_wall_data,
				 &wdata->odata,prt);

	(void) printf("\nIncident shock collides with contact at time %g\n",
	              (c-s)/W[1]);
	(void) printf("Transmitted wave = %s, velocity = %g\n",
	              rsoln_wave_name(rswtype[2]),W[2]);
	verbose_print_state("transmitted mid state",state[3]);
	(void) printf("Contact velocity = %g\n",W[3]);
	(void) printf("Reflected wave = %s, velocity = %g\n",
	              rsoln_wave_name(rswtype[4]),W[4]);
	verbose_print_state("reflected mid state",state[4]);

	(void) printf("First wall collision at time %g\n",wdata->t[1]);
	(void) printf("Reflected wave at wall = %s, velocity = %g\n",
	              rsoln_wave_name(rswtype[5]),W[5]);
	verbose_print_state("state on wall after collision",
	                    wdata->exact_wstate[1]);

	(void) printf("Second shock contact interaction at time %g\n",
	              wdata->t[1]+dWall*(1.0-W[3]/W[2])/(W[5]-W[3]));
	(void) printf("Reflected wave = %s, velocity = %g\n",
	              rsoln_wave_name(rswtype[6]),W[6]);
	verbose_print_state("reflected mid state",state[5]);
	(void) printf("Contact velocity = %g\n",W[7]);
	(void) printf("Transmitted wave = %s, velocity = %g\n",
	              rsoln_wave_name(rswtype[8]),W[8]);
	verbose_print_state("transmitted mid state",state[6]);
	(void) printf("Second wall collision at time %g\n",wdata->t[2]);

	(void) printf("Initial States\n");
	verbose_print_state("Target state",state[0]);
	verbose_print_state("Ahead shock state",state[1]);
	verbose_print_state("Behind shock state",state[2]);
	(void) printf("Incident shock velocity = %g\n",W[1]);
	(void) printf("Incident shock Mach number = %g\n",
	              fabs(W[1])/sound_speed(state[1]));

	verbose_print_state("State behind transmitted shock",state[3]);

	(void) printf("Wall states\n");
	for (i = 0; i < wdata->num_exact_wall_states; ++i)
	{
	    (void) sprintf(buf,"Wall state for %g < t < %g",
	                   wdata->t[i],wdata->t[i+1]);
	    verbose_print_state(buf,wdata->exact_wstate[i]);
	    (void) printf("\n");
	}
	free_these(2,left,right);
	for (i = 0; i < 7; ++i)
	    free(state[i]);
}		/*end init_oned_density_step*/

LOCAL	DENSITY_STEP_WALL_DATA *alloc_DENSITY_STEP_WALL_DATA(
	int    num_exact_wall_states,
	int    num_samples,
	Front *front)
{
	DENSITY_STEP_WALL_DATA   *wdata;
	size_t   sizest = front->sizest;
	INTERFACE *intfc = front->interf;
	int i, n, N, M;

	scalar(&wdata,sizeof(DENSITY_STEP_WALL_DATA));
	n = wdata->num_exact_wall_states = num_exact_wall_states;
	uni_array(&wdata->exact_wstate,n,sizeof(Locstate));
	uni_array(&wdata->t,n+1,FLOAT);
	for (i = 0; i < n; ++i)
	    alloc_state(intfc,wdata->exact_wstate+i,max(sizest,sizeof(VGas)));
	uni_array(&wdata->exact_incomingWs,n,FLOAT);
	uni_array(&wdata->exact_outgoingWs,n,FLOAT);
	N = wdata->num_samples = num_samples;
	uni_array(&wdata->pstate,N,sizeof(Locstate));
	uni_array(&wdata->pcomp,N,sizeof(COMPONENT));
	uni_array(&wdata->px,N,FLOAT);
	for (i = 0; i < N; ++i)
	    alloc_state(intfc,wdata->pstate+i,sizest);
	wdata->inside_wave = NO;
	M = wdata->NumAllocWallEvents = 10;
	uni_array(&wdata->wave_head_time,M,FLOAT);
	uni_array(&wdata->wave_tail_time,M,FLOAT);
	uni_array(&wdata->wave_mid_time,M,FLOAT);
	uni_array(&wdata->wave_head,M,sizeof(Locstate));
	uni_array(&wdata->wave_tail,M,sizeof(Locstate));
	for (i = 0; i < M; ++i)
	{
	    alloc_state(intfc,wdata->wave_head+i,max(sizest,sizeof(VGas)));
	    alloc_state(intfc,wdata->wave_tail+i,max(sizest,sizeof(VGas)));
	    wdata->wave_head_time[i] = HUGE_VAL;
	    wdata->wave_mid_time[i] = HUGE_VAL;
	    wdata->wave_tail_time[i] = HUGE_VAL;
	}
	wdata->num_events = 0;
	return wdata;
}		/*end alloc_DENSITY_STEP_WALL_DATA*/

LOCAL	void expand_num_wall_events(
	DENSITY_STEP_WALL_DATA *wdata,
	Front                  *front)
{
	INTERFACE *intfc = front->interf;
	size_t    sizest = front->sizest;
	int       i, M;
	double *pht = wdata->wave_head_time;
	double *ptt = wdata->wave_tail_time;
	double *mtt = wdata->wave_mid_time;
	Locstate *pwh = wdata->wave_head;
	Locstate *pwt = wdata->wave_tail;
	M = 2*wdata->NumAllocWallEvents;
	uni_array(&wdata->wave_head_time,M,FLOAT);
	uni_array(&wdata->wave_tail_time,M,FLOAT);
	uni_array(&wdata->wave_mid_time,M,FLOAT);
	uni_array(&wdata->wave_head,M,sizeof(Locstate));
	uni_array(&wdata->wave_tail,M,sizeof(Locstate));
	for (i = wdata->NumAllocWallEvents; i < M; ++i)
	{
	    alloc_state(intfc,wdata->wave_head+i,max(sizest,sizeof(VGas)));
	    alloc_state(intfc,wdata->wave_tail+i,max(sizest,sizeof(VGas)));
	    wdata->wave_head_time[i] = HUGE_VAL;
	    wdata->wave_mid_time[i] = HUGE_VAL;
	    wdata->wave_tail_time[i] = HUGE_VAL;
	}
	for (i = 0; i < wdata->NumAllocWallEvents; ++i)
	{
	    wdata->wave_head_time[i] = pht[i];
	    wdata->wave_tail_time[i] = ptt[i];
	    wdata->wave_mid_time[i] = mtt[i];
	    wdata->wave_head[i] = pwh[i];
	    wdata->wave_tail[i] = pwt[i];
	}
	wdata->NumAllocWallEvents = M;
}		/*end expand_num_wall_events*/

/*ARGSUSED*/
LOCAL	void	print_density_step_wall_data(
	Grid        *grid,
	Wave        *wave,
	Front       *front,
	Printplot   *prt,
	OUTPUT_DATA *data,
	boolean        about_to_stop)
{
	DENSITY_STEP_WALL_DATA *wdata = (DENSITY_STEP_WALL_DATA*)data;
	INTERFACE              *intfc = front->interf;
	FILE                   *file;
	double                  time;
	double                  Pe, Pn, dP, adP, rdP, ardP;
	double                  rhon;
	double                  Plast;
	int                    i, n;

	n = wdata->num_exact_wall_states;
	if (Output_file(data) == NULL)
	{
	    file = Output_file(data) = fopen(Output_filename(data),"w");
	    print_machine_parameters(file);
	    (void) foutput(file);
	    (void) fprintf(file,"%-18s %-18s %-18s\n",
	                   "Time","Pressure","EXACT-SOLUTION");
	    for (i = 0; i < n; ++i)
	    {
	        time = wdata->t[i];
	        Pe = pressure(wdata->exact_wstate[i]);
	        if (Output_in_binary(data))
	        {
	            (void) fprintf(file,"\f%c",4);
		    (void) fwrite((const void *)&time,FLOAT,1,file);
		    (void) fwrite((const void *)&Pe,FLOAT,1,file);
	        }
	        else
	            (void) fprintf(file,"%-"FFMT" %-"FFMT"\n",time,Pe);
	        time = wdata->t[i+1];
	        if (Output_in_binary(data))
	        {
	            (void) fprintf(file,"\f%c",4);
		    (void) fwrite((const void *)&time,FLOAT,1,file);
		    (void) fwrite((const void *)&Pe,FLOAT,1,file);
	        }
	        else
	            (void) fprintf(file,"%-"FFMT" %-"FFMT"\n",time,Pe);
	    }
	    (void) fprintf(file,"\n");

	    (void) foutput(file);
	    (void) fprintf(file,"%-18s %-18s %-18s %-18s %-18s "
	                        "%-18s %-18s %-18s\n",
	                   "Time","P-numerical","P-exact","P-error","aP-error",
			   "%rel-P-error","%arel-P-error",
			   "NUMERICAL-SOLUTION");
	}
	else
	    file = Output_file(data);

	time = grid->time;
	Plast = wdata->Plast;
	Pn = 0.0;
	rhon = 0.0;
	for (i = 0; i < wdata->num_samples; ++i)
	{
	    wdata->pcomp[i] = component(wdata->px+i,intfc);
	    hyp_solution(wdata->px+i,wdata->pcomp[i],NULL,UNKNOWN_SIDE,
		         front,wave,wdata->pstate[i],NULL);
	    Pn += pressure(wdata->pstate[i]);
	    rhon += Dens(wdata->pstate[i]);
	}
	Pn /= wdata->num_samples;
	rhon /= wdata->num_samples;
	dP = 2.0*fabs(Pn - Plast)/(Pn + Plast);
	if (wdata->num_events >= wdata->NumAllocWallEvents)
	    expand_num_wall_events(wdata,front);
	if (wdata->inside_wave == NO)
	{
	    if (dP > wdata->wave_tol)
	    {
	        int n = wdata->num_events;
	        wdata->inside_wave = YES;
	        wdata->wave_head_time[n] = wdata->wave_mid_time[n] =
	            wdata->wave_tail_time[n] = time;
	        set_type_of_state(wdata->wave_head[n],TGAS_STATE);
	        Dens(wdata->wave_head[n]) = rhon;
	        Press(wdata->wave_head[n]) = Pn;
	        zero_state_velocity(wdata->wave_head[n],front->rect_grid->dim);
	        Set_params(wdata->wave_head[n],wdata->pstate[0]);
	        set_state(wdata->wave_head[n],VGAS_STATE,wdata->wave_head[n]);
	        set_state(wdata->wave_tail[n],VGAS_STATE,wdata->wave_head[n]);
	    }
	}
	else
	{
	    if (dP < wdata->wave_tol)
	    {
	        int n = wdata->num_events;
	        wdata->wave_tail_time[n] = time;
		wdata->wave_mid_time[n] =
		    0.5*(wdata->wave_head_time[n]+wdata->wave_tail_time[n]);
	        set_type_of_state(wdata->wave_tail[n],TGAS_STATE);
	        Dens(wdata->wave_tail[n]) = rhon;
	        Press(wdata->wave_tail[n]) = Pn;
	        zero_state_velocity(wdata->wave_tail[n],front->rect_grid->dim);
	        Set_params(wdata->wave_tail[n],wdata->pstate[0]);
	        set_state(wdata->wave_tail[n],VGAS_STATE,wdata->wave_tail[n]);
		++wdata->num_events;
	        wdata->inside_wave = NO;
	    }
	}
	wdata->Plast = Pn;

	if (time < wdata->t[n])
	{
	    for (i = 0; i < n; ++i)
	        if (wdata->t[i] <= time && time <= wdata->t[i+1])
		    break;
	    Pe = pressure(wdata->exact_wstate[i]);
	    dP = Pn-Pe;
	    rdP = (Pn-Pe)/Pe;
	    adP = fabs(dP);
	    ardP = fabs(rdP);
	    if (Output_in_binary(data))
	    {
	        double prdP = 100.0*rdP, pardP = 100.0*ardP;
	        (void) fprintf(file,"\f%c",4);
		(void) fwrite((const void *)&time,FLOAT,1,file);
		(void) fwrite((const void *)&Pn,FLOAT,1,file);
		(void) fwrite((const void *)&Pe,FLOAT,1,file);
		(void) fwrite((const void *)&dP,FLOAT,1,file);
		(void) fwrite((const void *)&adP,FLOAT,1,file);
		(void) fwrite((const void *)&prdP,FLOAT,1,file);
		(void) fwrite((const void *)&pardP,FLOAT,1,file);
	    }
	    else
	        (void) fprintf(file,"%-"FFMT" %-"FFMT" %-"FFMT" %-"FFMT" "
		               "%-"FFMT" %-"FFMT" %-"FFMT"\n",
	                       time,Pn,Pe,dP,adP,100.0*rdP,100.0*ardP);
	}
	else if (wdata->num_events > 0)
	{
	    static boolean first = YES;
	    if (first)
	    {

	        double dat;
	        double Wexact;
	        double Wmid, Ws, We, Werr, dW;
		double tail_time, head_time, mid_time;

	        first = NO;
	        Wexact = wdata->exact_incomingWs[0];
	        (void) fprintf(file,"\n");
	        (void) fprintf(file,"Exact arrival time at wall "
	                            "of transmitted shock = %g\n",wdata->t[1]);
	        (void) fprintf(file,"Exact transmitted wave speed = %g\n",
	                       Wexact);

		tail_time = wdata->wave_tail_time[0];
		mid_time = wdata->wave_mid_time[0];
		head_time = wdata->wave_head_time[0];
	        (void) fprintf(file,"\nEstimates for inital "
		                    "transmitted shock\n");
	        (void) fprintf(file,"start of arrival = %g\n",head_time);
	        (void) fprintf(file,"middle of arrival = %g\n",mid_time);
	        (void) fprintf(file,"end of arrival = %g\n",tail_time);
	        dat = mid_time - wdata->t[1];
	        (void) fprintf(file,"Error in middle of arrival = %g\n",dat);
	        (void) fprintf(file,"Absolute relative error in "
	                            "middle of arrival = %g%%\n",
	                            100.0*fabs(dat)/wdata->t[1]);
	        (void) fprintf(file,"Width of arrival interval = %g\n",
	                       tail_time-head_time);
	        (void) fprintf(file,"Relative width of "
		                    "arrival interval = %g%%\n",
	                       100.0*(tail_time - head_time)/wdata->t[1]);

	        Ws = -(wdata->c - front->rect_grid->L[0])/
	              (head_time - wdata->t[0]);
	        (void) fprintf(file,"Estimated transmitted wave "
	                       "head velocity = %g\n",Ws);
	        Werr = Ws - Wexact;
	        (void) fprintf(file,"Error in estimated transmitted wave "
	                       "head velocity = %g\n",Werr);
	        (void) fprintf(file,"Relative error in "
	                            "transmitted wave head velocity = %g%%\n",
	                            100.0*fabs(Werr/Wexact));
    
	        Wmid = -(wdata->c-front->rect_grid->L[0])/
	                (mid_time - wdata->t[0]);
	        (void) fprintf(file,"Estimated transmitted wave midpoint "
	                       "velocity = %g\n",Wmid);
	        Werr = Wmid - Wexact;
	        (void) fprintf(file,"Error in estimated transmitted wave "
	                       "velocity = %g\n",Werr);
	        (void) fprintf(file,"Absolute relative error in "
	                            "transmitted wave velocity = %g%%\n",
	                            100.0*fabs(Werr/Wexact));
	                       
	        We = -(wdata->c-front->rect_grid->L[0])/
	              (tail_time - wdata->t[0]);
	        (void) fprintf(file,"Estimated transmitted wave "
	                       "tail velocity = %g\n",We);
	        Werr = We - Wexact;
	        (void) fprintf(file,"Error in estimated transmitted wave "
	                       "tail velocity = %g\n",Werr);
	        (void) fprintf(file,"Relative error in "
	                            "transmitted wave tail velocity = %g%%\n",
	                       100.0*fabs(Werr/Wexact));
	        dW = 100.0*0.5*fabs((We-Ws)/Wmid);
	        (void) fprintf(file,"Transmitted wave velocity = %g +- %g%%\n",
	                       Wmid,dW);
	    }
	}
	if (about_to_stop)
	{
	    char name[256];
	    (void) fprintf(file,"\nDetected Wall Events\n");
	    for (i = 0; i < wdata->num_events; ++i)
	    {
	        (void) fprintf(file,"Time of wave %d-%s head = %g\n",i+1,
		              ordinal_suffix(i+1),
		              wdata->wave_head_time[i]);
	        (void) fprintf(file,"Time of %d-%s wave middle = %g\n",i+1,
		              ordinal_suffix(i+1),
		              wdata->wave_mid_time[i]);
	        (void) fprintf(file,"Time of wave %d-%s tail = %g\n",i+1,
		              ordinal_suffix(i+1),
		              wdata->wave_tail_time[i]);
		(void) fprintf(file,"pressure jump across wave = %g\n",
		               pressure(wdata->wave_tail[i])-
			       pressure(wdata->wave_head[i]));
		(void) sprintf(name,"State at wave head");
		verbose_fprint_state(file,name,wdata->wave_head[i]);
		(void) sprintf(name,"State at wave tail");
		verbose_fprint_state(file,name,wdata->wave_tail[i]);
	    }
	}
}		/*end print_density_step_wall_data*/

/*
 * storage for time dependent pressure input
 * time_dep_pres = { warming-up time,
 *                   peak time,
 *                   cooling-down time,
 *                   base pressure,
 *                   peak pressure }
 *
 *              set_time_dep_pres();
 *      Set warming-up, peak, cooling-down times, base and peak pressures
 *      into the array time_dep_pres
 */
LOCAL 	void    set_time_dep_pres(
	BOUNDARY_STATE *bstate,
	double   tr,
	double   tp,
	double   ts,
	double   pr_a,
	double   pr_p,
        Locstate state)
{
	FD_DATA *fd_data;
	stat_scalar(&bstate->_boundary_state_data,sizeof(FD_DATA));
	fd_data = (FD_DATA*)bstate->_boundary_state_data;
        fd_data->tr = tr;
        fd_data->tp = tp;
        fd_data->ts = ts;
        fd_data->pr_a = pr_a;
        fd_data->pr_p = pr_p;
	
	/* Feb 19 2004: Myoung-Nyoun: Fix for general eos */
	if (state == NULL)
	{
            printf("ERROR in set_time_dep_pres(), state == NULL\n");
            clean_up(ERROR);
        }
        fd_data->state = state;
}               /* end set_time_dep_pres */

EXPORT  void      g_init_parabolic(
        Front     *front)
{
        INTERFACE *intfc = front->interf;
	Gas_param **prms_list; /* = gas_params_list(intfc);*/
	int       nprms;
	int       i, nz, *idx=NULL;
        char      s[121];
	char      fmt[1024],fmt2[1024];
	double    temp = 0.0;
        boolean      temp_boolean;

	debug_print("init","Entered g_init_parabolic()\n");

	screen("\n\n\t\tSpecify parabolic steps\n\n");
        screen("Type 'y' to have the Navier-Stokes terms computed for"
	       "\n\tseveral eos models, and this will turn on"
	       "\n\tparabolic driver parab_driver (y, n(dflt)): ");
        (void) Gets(s);
	if (strncmp(s,"y",1) == 0)
	{        
            nprms = return_params_list(&prms_list);
            screen("Enter NS viscosity coefficient(dflt = 0): ");
            (void) Gets(s);
            if (s[0] != '\0')
            {
                if (sscan_float(s,&temp) != 1)
                {
                    screen("ERROR in g_init_parabolic(), incorrect input "
                       "of viscosity coefficient\n");
                    clean_up(ERROR);
                }
            }
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.viscosity_coef=temp;
            temp_boolean = YES;
            screen("Use Stokes hypothesis for viscosity term(y(dflt), n): ");
            (void) Gets(s);
            if (strncmp(s,"n",1) == 0)
                temp_boolean = NO;
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.use_stokes_vis=temp_boolean;
            temp = 0.0;
            screen("Enter NS mass diffusion coefficient(dflt = 0): ");
            (void) Gets(s);
            if (s[0] != '\0')
            {
                if (sscan_float(s,&temp) != 1)
                {
                    screen("ERROR in g_init_parabolic(), incorrect input "
                                      "of mass diffusion coefficient\n");
                    clean_up(ERROR);
                }
            }
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.diffusivity_coef=temp;
            temp = 0.0;
            screen("Enter NS thermal conductivity coefficient(dflt = 0): ");
            (void) Gets(s);
            if (s[0] != '\0')
            {
                if (sscan_float(s,&temp) != 1)
                {
                    screen("ERROR in g_init_parabolic(), incorrect input "
                                      "of thermal conductivity coefficient\n");
                    clean_up(ERROR);
                }
            }
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.conductivity_coef=temp;
#if defined(SUBGRID)
            temp_boolean = NO;
            screen("Type 'y' to have the turbulence simulation for dynamic model (y, n(dflt)): ");
            (void) Gets(s);
            if (strncmp(s,"y",1) == 0)
                temp_boolean = YES;
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.turbulence=temp_boolean;
            if (strncmp(s,"y",1) == 0)
            {
               chart_of_front(front)->subgrid = SGS;

               temp = 0.0;
               screen("Enter the time to start turbulent simulation: ");
               (void) Gets(s);
               if (s[0] != '\0')
               {
                   if (sscan_float(s,&temp) != 1)
                   {
                       screen("ERROR in g_init_parabolic(), incorrect input "
                          "of the time to start turbulent simulation\n");
                       clean_up(ERROR);
                   }
               }
               front->subgrid_time = temp;

            temp_boolean = YES;
            screen("Use planar average for subgrid model (y(dflt), n): ");
            (void) Gets(s);
            if (strncmp(s,"y",1) == 0)
                temp_boolean = YES;
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.aver_planar=temp_boolean;            
            temp_boolean = NO;
            screen("Use Subgrid model for viscosity (y, n(dflt)): ");
            (void) Gets(s);
            if (strncmp(s,"y",1) == 0)
                temp_boolean = YES;
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.subgrid_vis=temp_boolean;
            temp_boolean = NO;
            screen("Use Subgrid model for mass diffusion (y, n(dflt)): ");
            (void) Gets(s);
            if (strncmp(s,"y",1) == 0)
                temp_boolean = YES;
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.subgrid_md=temp_boolean; 
            temp_boolean = NO;
            screen("Use Subgrid model for thermal conductivity (y, n(dflt)): ");
            (void) Gets(s);
            if (strncmp(s,"y",1) == 0)
                temp_boolean = YES;
            for (i = 0; i < nprms; i++)
                prms_list[i]->avisc.subgrid_con=temp_boolean;
            }
#endif /* defined SUBGRID */

	    /* Apr 15 2003: Myoung-Nyoun:
	     * moved from g_set_basic_phys_parameters */
	    chart_of_front(front)->parab = parab_driver;
#if defined(USE_OVERTURE)
            chart_of_front(front)->parab_npt = parab_npt;
#endif /* if defined(USE_OVERTURE) */
	    screen("\nCurrent gas param list\n");
	    nprms = return_params_list(&prms_list);
	    screen("Number of params = %d\n\n",nprms);
	    for (i=0; i<nprms; i++)
	    {
		IMPORT boolean suppress_prompts;
		screen ("Param[%d]\n",i);
		fprint_Gas_param(stdout,prms_list[i]);
		if (suppress_prompts == NO)
		    fprint_Gas_param(stderr,prms_list[i]);
		screen("\n");
	    }

	    uni_array(&idx,nprms,sizeof(i));

	    screen("Enter indices of Param (less than %d) to use "
		   "Navier-Stokes terms: ",nprms);
	    (void) Gets(s);

	    strcpy(fmt,"%d");
	    strcpy(fmt2,"%*d ");

	    for (i=0; i<nprms; i++)
	    {
		 nz = sscanf(s,fmt,&idx[i]);
		 if (nz <= 0) break;
		 strcat(fmt2,fmt);
		 strcpy(fmt,fmt2);
		 strcpy(fmt2,"%*d ");
	    }
	    nz = i;

	    for (i=0; i<nz; i++)
	    {
		if (idx[i] >= nprms || idx[i] < 0)
		{
		    (void) printf("ERROR in g_init_parabolic, index is out"
				  " of range, index = %d, nprms = %d\n"
				  ,idx[i],nprms);
		    clean_up(ERROR);
		}
	    }

	    for (i=0; i<nprms; i++)
		prms_list[i]->eos->_compute_ns_terms = NO;
	    for (i=0; i<nz; i++)
		prms_list[idx[i]]->eos->_compute_ns_terms = YES;
	}

	if (debugging("parabolic_step"))
	{
	    (void) printf("\nParabolic step (EOS indices) =");
	    for (i=0; i<nz; i++)
	        (void) printf(" %d",idx[i]);
	    (void) printf("\nNumber of indices to turn on flag = %d\n",nz);
	}
	screen("\n\n");
	if (idx) free(idx);
	debug_print("init","Left g_init_parabolic()\n");
}		/*end g_init_parabolic*/


LOCAL int prompt_for_time_dependent_boundary_state(
         int index,
         Gas_param *params,
         INTERFACE *intfc)
{
	int n;
	char  s[Gets_BUF_SIZE];
	BOUNDARY_STATE	Bstate;
        TD_BSTATE *tds;

	Bstate._boundary_state = NULL;
	Bstate._boundary_state_function = g_time_dependent_boundary_state;
	Bstate._boundary_state_function_name =
				strdup("g_time_dependent_boundary_state");
        stat_scalar(&tds,sizeof(TD_BSTATE));
	Bstate._boundary_state_data = (POINTER)tds;
	Bstate._fprint_boundary_state_data = 
				g_print_time_dependent_boundary_state;
        /* prompt for the file name with the data*/
	screen("\nFile name with boundary data: ");
        (void) Gets(s);
        if (s[0] == '\0')
        {
             screen("ERROR in prompt_for_time_dependent_boundary_state() "
                    "no file name specified\n");
             clean_up(ERROR);
        }
        strcpy(tds->fn,s);
        tds->params = params;
        read_time_dependent_data(tds);


	return add_bstate_to_list(&Bstate,intfc,index);
}                /*end prompt_time_dependent_boundary_state*/

EXPORT void read_time_dependent_data(
    TD_BSTATE *tds)
{
	int i, c, numlines;
	char  s[Gets_BUF_SIZE];
	FILE *fp;
	/* open the file */
	fp = fopen(tds->fn, "r");
        if (fp == NULL)
        {
            screen("ERROR in read_time_dependent_data() "
                   "can't open %s for reading\n",tds->fn);
            clean_up(ERROR);
        }

	/* count the number of lines in the file (fgets gets one line) */
	for (numlines = -1; fgets(s,Gets_BUF_SIZE,fp)!=NULL; ++numlines);
        	tds->n = numlines;

	/*stat_vectors*/
        stat_vector(&tds->p,tds->n,FLOAT);
        stat_vector(&tds->rho,tds->n,FLOAT);
        stat_vector(&tds->v,tds->n,FLOAT);
        stat_vector(&tds->t,tds->n,FLOAT);

        /* rewind the file */
	rewind(fp);

        /* skip the first line (header)*/
	fgets(s,Gets_BUF_SIZE,fp);

        /* read values from file */
        for (i = 0; i < tds->n; ++i)
        {
            fscan_float(fp,tds->t+i);
            fscan_float(fp,tds->p+i);
            fscan_float(fp,tds->rho+i);
            fscan_float(fp,tds->v+i);
        }

        /* close the file */
	fclose(fp);
}/*end read_time_dependent_data*/

