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
*				dpatchmesh.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for initializing and maintaining nodes in the
*	automatic mesh grid level tree.
*/


#include <driver/ddecs.h>

	/* LOCAL Function Declarations */
LOCAL	AFLIN	*read_print_AFLIN(FILE*,int);
LOCAL	COMPONENT	comp_of_bdry(HYPER_SURF*,int,CHART*,CHART*);
LOCAL	int	integrate_at_level(CHART*,LEVEL*,double*);
LOCAL	int	interior_coord_in_best_chart(double*,CHART*,double*,CHART**,int);
LOCAL	int	is_strict_interior_point(double*,CHART*,double*,CHART*,double*);
LOCAL	int	prompt_for_child_chart(CHART*);
LOCAL	void	find_next_level(CHART*,LEVEL**);
LOCAL	void	init_amr_hyp_solution(double*,COMPONENT,HYPER_SURF*,SIDE,
				      Front*,Wave*,Locstate);
LOCAL	void	init_level(CHART*,LEVEL*);
LOCAL	void	intfc_point_and_states_from_finer_levels(POINT*,BOND*,
							 CURVE*,CHART*);
#if defined(THREED)
LOCAL	void	prompt_for_3d_rotation(double**);
#endif /* defined(THREED) */
LOCAL	void	regrid_at_level(LEVEL*);
LOCAL	void	set_chart(CHART*,double*,double*,double**);
LOCAL	void	set_coord_transform(CHART*,double*,double**,int);
LOCAL	void	set_wave_node_types(CHART*);
LOCAL	void	set_zoomed_front_states(CHART*);
LOCAL	void	set_zoomed_grid(RECT_GRID*,double**,RECT_GRID*,AFLIN*,int);
LOCAL	void	set_zoomed_interior_states(CHART*);
LOCAL	void	solution_from_finer_levels(double*,COMPONENT,CHART*,Locstate);
LOCAL	void	transform_coord(double*,CHART*,double*,CHART*);

enum {
	FROM_COARSER = 0,
	FROM_BEST
};

LOCAL	int regrid_time_step;
LOCAL	int states_from_coarser_or_best_chart = FROM_COARSER;

#if !defined(__INTEL_COMPILER)
#pragma	noinline	read_print_AFLIN
#pragma	noinline	comp_of_bdry
#pragma	noinline	integrate_at_level
#pragma	noinline	interior_coord_in_best_chart
#pragma	noinline	is_strict_interior_point
#pragma	noinline	prompt_for_child_chart
#pragma	noinline	find_next_level
#pragma	noinline	init_amr_hyp_solution
#pragma	noinline	init_level
#pragma	noinline	intfc_point_and_states_from_finer_levels
#if defined(THREED)
#pragma	noinline	prompt_for_3d_rotation
#endif /* defined(THREED) */
#pragma	noinline	regrid_at_level
#pragma	noinline	set_chart
#pragma	noinline	set_coord_transform
#pragma	noinline	set_wave_node_types
#pragma	noinline	set_zoomed_front_states
#pragma	noinline	set_zoomed_grid
#pragma	noinline	set_zoomed_interior_states
#pragma	noinline	solution_from_finer_levels
#pragma	noinline	transform_coord
#endif /*!defined(__INTEL_COMPILER)*/

/*
*				init_amr():
*
*	Main initializer, called by dinit.c
*/

/*ARGSUSED*/
EXPORT void init_amr(
	CHART		*root)
{
	int		step = root->grid->step;
	double		t = root->grid->time, dt = root->grid->dt;
	Front		*front = root->front;
	Wave		*wave = root->wave;
	LEVEL		*level_str;
	char		s[Gets_BUF_SIZE];

	debug_print("chart","Entered init_amr(), front wave = %d %d\n",front,wave);
	chart_of_front(front) = chart_of_wave(wave) = root;
	screen("\nType 'y' to have local mesh refinement: ");
	(void) Gets(s);
	if (s[0] != 'y' && s[0] != 'Y')
	{
	    root->level = NULL;
	    debug_print("chart","Left init_amr(), front wave = %p %p\n",
	          (POINTER)front,(POINTER)wave);
	    return;
	}
	screen("\n-- Begin prompting for mesh refinement --\n");
	screen("\nEnter time step for regridding\n");
	screen("(greater than max number of time steps in no regridding): ");
	(void) Scanf("%d\n",&regrid_time_step);
	(void) printf("regrid_time_step is = %d\n",regrid_time_step);

		/* initialize LEVEL */

	init_level(root,(LEVEL *) NULL);
	level_str = root->level;
	screen("Enter steps between synchronization and between regrid: ");
	(void) Scanf("%d %d\n",&level_str->steps_between_synchronize,
		&level_str->steps_between_regrid);
	level_str->steps_before_synchronize =
		level_str->steps_between_synchronize;
	level_str->steps_before_regrid =
		level_str->steps_between_regrid;
	screen("Enter (positive integer) refinement ratio: ");
	(void) Scanf("%d\n",&level_str->refinement_ratio);
	level_str->steps_before_level_change = level_str->refinement_ratio;
	level_str->step = step;
	level_str->t = t;
	level_str->dt = dt;

	if (debugging("chart"))
	{
		print_CHART(stdout,root);
		print_LEVEL(stdout,level_str);
	}

		/* Prompt for child charts */

	while (prompt_for_child_chart(root))
	        ;
	screen("\n-- End prompting for mesh refinement --\n");
	if (debugging("chart"))
	{
		CHART *ch;
		LEVEL *le;
		for (le = root->level; le; le = le->next_finer_level)
			for (ch = le->first;ch;ch=ch->next_chart)
				graph_front_states(ch->front);
	}
	debug_print("chart","Left init_amr(), front wave = %p %p\n",
	      (POINTER)front,(POINTER)wave);
}		/*end init_amr*/



/*
*				init_level():
*
*	Initializes the LEVEL structure.
*/

LOCAL void init_level(
	CHART		*chart,
	LEVEL		*parent_level)
{
	LEVEL		*level_str;
	static int	max_num_levels;

	scalar(&level_str,sizeof(LEVEL));
	chart->level = level_str;
	level_str->first = chart;
	level_str->last = chart;
	level_str->next_finer_level = NULL;
	level_str->next_coarser_level = parent_level;
	if (parent_level)
	{
		parent_level->next_finer_level = level_str;
		level_str->level_num = parent_level->level_num + 1;
		level_str->max_num_levels = parent_level->max_num_levels;
		*level_str->max_num_levels += 1;
		level_str->steps_before_synchronize = 
			parent_level->steps_between_synchronize;
		level_str->steps_between_synchronize = 
			parent_level->steps_between_synchronize;
		level_str->steps_before_regrid = 
			parent_level->steps_between_regrid;
		level_str->steps_between_regrid = 
			parent_level->steps_between_regrid;
		level_str->steps_before_level_change =
			parent_level->refinement_ratio;
		level_str->refinement_ratio = parent_level->refinement_ratio;
		level_str->step = parent_level->step;
		level_str->t = parent_level->t;
		level_str->dt = parent_level->dt / level_str->refinement_ratio;
	}
	else
	{
		level_str->max_num_levels = &max_num_levels;
		*level_str->max_num_levels = 0;
		level_str->level_num = 0;
	}
	chart->next_chart = NULL;
	chart->prev_chart = NULL;
}		/*end init_level*/

/*
*				prompt_for_child_chart():
*
*	Initializes a single child CHART, with prompts for rectangle position.
*/

LOCAL int prompt_for_child_chart(
	CHART		*chart)
{
	char		s[Gets_BUF_SIZE];
	CHART		*child;
	double		pc[MAXD], pu[MAXD];
	int		i, dim = chart->front->rect_grid->dim;
	static boolean	first = YES;
	static double	**rot_matrix = NULL;

	debug_print("chart","Entered prompt_for_child_chart()\n");
	if (first)
	{
	    first = NO;
	    bi_array(&rot_matrix,dim,dim,FLOAT);
	    s[0] = 'Y';
	}
	else
	{
	    screen("Request (another) level %d mesh refined grid\n",
	    	   chart->level->level_num + 1);
	    screen("\tEnter yes or (default) no: ");
	    (void) Gets(s);
	}
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    scalar(&child,sizeof(CHART));

	    child->parent = chart;

	    /* Initialize child in level linked list */

	    if (chart->level->level_num == *chart->level->max_num_levels)
	    	init_level(child,chart->level);
	    else
	    {
	    	child->level = chart->level->next_finer_level;
	    	child->prev_chart = child->level->last;
	    	child->level->last = child;
	    }
	    screen("Is new grid static (YES or NO): ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y') 
	    	child->dynamic = NO;
	    else
	    	child->dynamic = YES;
	    child->is_old_chart = NO;
	    screen("Enter lower corner of child domain "
	           "(in root coordinate system): ");
	    for (i = 0; i < dim; i++)  
	        (void) Scanf("%f",pc+i);
	    (void) Scanf("\n");
	    screen("Enter root coord. of upper corner of domain: ");
	    for (i = 0; i < dim; i++) 
	            (void) Scanf("%f",pu+i);
	    (void) Scanf("\n");
	    switch (dim)
	    {
#if defined(ONED)
	    case 1:
	    	rot_matrix[0][0] = 1.0;
	    	break;
#endif /* defined(ONED) */
#if defined(TWOD)
	    case 2:
	    /*  allow rotations later ...  
	    {
                double angle;
	    	screen("Enter angle (degrees) root x axis to ");
		screen("child x axis: ");
		(void) Scanf("%f\n",&angle);
		angle = radians(angle);
		rot_matrix[0][0] = rot_matrix[1][1] = cos(angle);
		rot_matrix[1][0] = sin(angle);
		rot_matrix[0][1] = -rot_matrix[1][0];
	    }
	    */
	        rot_matrix[0][0] = rot_matrix[1][1] = 1.0;
	        rot_matrix[0][1] = rot_matrix[1][0] = 0.0;
	        break;
#endif /* defined(TWOD) */
#if defined(THREED)
	    case 3:
		prompt_for_3d_rotation(rot_matrix);
		break;
#endif /* defined(THREED) */
	    }
	    set_chart(child,pc,pu,rot_matrix);

/* The hyperbolic solution method should be set here.  Temporarily
it is simply set equal to that of the parent chart.    */

	    child->hyp_solver = child->parent->hyp_solver;
	    child->bc_propagate = child->parent->bc_propagate;
	    child->parab = child->parent->parab;
	    child->ellip = child->parent->ellip;
	    child->grid = child->parent->grid;
	    debug_print("chart","Leaving prompt_for_child_chart\n");
	    while (prompt_for_child_chart(child))
			;
	    return YES;
	}
	screen("Returning to level %d parent\n",chart->level->level_num);
	debug_print("chart","Leaving prompt_for_child_chart\n");
	return NO;
}		/*end prompt_for_child_chart*/

/*
*				set_chart():
*
*	Sets the data in a single CHART, using given rectangle position and
*	the front and wave of the parent CHART.
*/


LOCAL void set_chart(
	CHART		*chart,
	double		*pc,	/* rectangle opposite corners parent coord. */
	double		*pu,
	double		**rot_matrix)/* root x axis to child x axis rotation */
{
	INTERFACE	*oldintfc;
	Front		*fr;
	RECT_GRID	*gr,*parent_gr;
	int		dim, i;

	oldintfc = current_interface();
	dim = oldintfc->dim;

		/* Set coordinate transformations */

	debug_print("chart","Entered set_chart\n");
	scalar(&chart->to_root,sizeof(AFLIN));
	scalar(&chart->from_root,sizeof(AFLIN));
	set_coord_transform(chart,pc,rot_matrix,dim);
	if (debugging("chart"))
	{
		print_AFLIN(stdout,chart->to_root,dim);
		print_AFLIN(stdout,chart->from_root,dim);
	}

		/* Initialize the front */

	chart->front = fr = copy_front(chart->parent->front);
	fr->interf = NULL;
	chart_of_front(fr) = chart;
	if (debugging("chart"))
	{
		print_CHART(stdout,chart);
		print_LEVEL(stdout,chart->level);
		(void) printf("chart_of_front(chart->front) is %p, ",
			      (POINTER)chart_of_front(chart->front));
		(void) printf("and chart_of_front(fr) is %p\n",
			      (POINTER)chart_of_front(fr));
	}

		/* Initialize the grid to parent grid */

	scalar(&gr,sizeof(RECT_GRID));
	fr->rect_grid = gr;
	parent_gr = chart->parent->front->rect_grid;
	copy_rect_grid(gr,parent_gr);
	chart->newfront = copy_front(chart->parent->front);
	chart->newfront->interf = NULL;

		/* Zoom the interface and front */

	set_copy_intfc_states(NO);
	set_size_of_intfc_state(size_of_state(chart->parent->front->interf));
	fr->interf = zoom_interface(chart->parent->front->interf,gr,
			pc,pu,rot_matrix);
	set_copy_intfc_states(YES);
	if (debugging("chart"))
	{
	        (void) printf("\n\nGrids before zoom :\n");
		(void) printf("\ntopological_grid(fr->interf)   address   %p",
			      (POINTER)&topological_grid(fr->interf));
		print_RECT_GRID_structure(&topological_grid(fr->interf));
		(void) printf("\ncomputational_grid(fr->interf)   address   %p",
			      (POINTER)computational_grid(fr->interf));
		print_RECT_GRID_structure(computational_grid(fr->interf));
		(void) printf("\ngr = fr->rect_grid  address   %p",(POINTER)gr);
		print_RECT_GRID_structure(gr);
	}
	for (i = 0; i < dim; i++)
	{
		gr->L[i] = topological_grid(fr->interf).VL[i];
		gr->U[i] = topological_grid(fr->interf).VU[i];
		gr->VL[i] = topological_grid(fr->interf).L[i];
		gr->VU[i] = topological_grid(fr->interf).U[i];
	}
	set_computational_grid(fr->interf,gr);
	set_zoomed_grid(gr,rot_matrix,parent_gr,chart->parent->to_root,
		chart->level->refinement_ratio);
	set_topological_grid(fr->interf,gr);
	set_zoomed_grid(&topological_grid(fr->interf),rot_matrix,
		&topological_grid(chart->parent->front->interf),
		chart->parent->to_root,chart->level->refinement_ratio);
	if (debugging("chart"))
	{
	        (void) printf("\n\nGrids after zoom :\n");
		(void) printf("\ntopological_grid(fr->interf)   address   %p",
			      (POINTER)&topological_grid(fr->interf));
		print_RECT_GRID_structure(&topological_grid(fr->interf));
		(void) printf("\ncomputational_grid(fr->interf)   address   %p",
			      (POINTER)computational_grid(fr->interf));
		print_RECT_GRID_structure(computational_grid(fr->interf));
		(void) printf("\ngr = fr->rect_grid  address   %p",(POINTER)gr);
		print_RECT_GRID_structure(gr);
	}
	set_wave_node_types(chart);
#if defined(TWOD) || defined(THREED)
	Redistribution_count(fr) = 0;
	if (redistribute(fr,YES,NO) != GOOD_REDISTRIBUTION)
	{
		screen("Redistribution failed in set_chart()\n");
		clean_up(ERROR);
	}
#endif /* defined(TWOD) || defined(THREED) */
	if (debugging("chart"))
	{
	    (void) output();
	    (void) printf("\n\nHere is the final zoomed and redist. intfc:\n");
	    print_interface(fr->interf);
	}
	interpolate_intfc_states(fr->interf) = YES;

		/* Zoom on wave structure; set hyp */

	chart->hyp_solver = chart->parent->hyp_solver;
	chart->bc_propagate = chart->parent->bc_propagate;
	chart->parab = chart->parent->parab;
	chart->ellip = chart->parent->ellip;
	chart->grid = chart->parent->grid;
	chart->wave = copy_wave(chart->parent->wave);
	clear_wave_pointers(chart->wave);
	chart->wave->rect_grid = gr;
	chart->newwave = copy_wave(chart->parent->wave);
	set_zoomed_front_states(chart);
	set_zoomed_interior_states(chart);

		/* TODO: set bdry data for wave and front */
		/* TODO: set source data for res only */
	if (debugging("chart"))
	{
		print_CHART(stdout,chart);
		print_Front_structure(chart->front);
		print_Wave_structure(chart->wave);
	}
	set_current_interface(oldintfc);
	debug_print("chart","Leaving set_chart\n");
}		/*end set_chart*/

LOCAL void set_zoomed_grid(
	RECT_GRID	*gr,
	double		**rot_matrix,
	RECT_GRID	*parent_gr,
	AFLIN		*ptor,	/* parent to root transformation */
	int		refinement_ratio)
{
	double		*h = gr->h, *parent_h = parent_gr->h;
	int		i, j, dim = parent_gr->dim;
	int		*gmax = gr->gmax;
	static boolean	first = YES;
	static double	**t_ctop = NULL, **M = NULL;

	debug_print("chart","Entering set_zoomed_grid()\n");

	if (first)
	{
		first = NO;
		bi_array(&t_ctop,MAXD,MAXD,FLOAT);
		bi_array(&M,MAXD,MAXD,FLOAT);
	}

	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
		{
			M[i][j] = ptor->a[j][i];
		}
	}

	rotate_matrix(t_ctop,M,rot_matrix,dim);

	for (i = 0; i < dim; i++)
	{
		h[i] = grid_size_in_direction(t_ctop[i],parent_h,dim);
		gmax[i] = (int) ((gr->U[i] - gr->L[i]) / h[i]);
		if (gmax[i] == 0) gmax[i] = 1;
		gmax[i] *= refinement_ratio;
		h[i] = (gr->U[i] - gr->L[i]) / gmax[i];
	}
	debug_print("chart","Leaving set_zoomed_grid()\n");
}		/*end set_zoomed_grid*/

/*ARGSUSED*/
LOCAL void set_wave_node_types(
	CHART		*chart)
{
#if defined(TWOD)
	COMPONENT	comp_par;
	CURVE		**c;
	CURVE		*c_par,*tempc_par;
	HYPER_SURF	*hs_par;
	HYPER_SURF_ELEMENT *hse_par;
	INTERFACE	*intfc = chart->front->interf;
	INTERFACE	*intfc_par = chart->parent->front->interf;
	NODE		**n;
	double		coords[MAXD], coords_par[MAXD];
	double		t[MAXD];
	ORIENTATION	tempc_orient;
	int		fixed_node;
	int		i, dim = intfc->dim;

	debug_print("chart","Entering set_wave_node_types\n");
	for (c = intfc->curves; *c; c++)
	{
	    if (wave_type(*c) == ERROR)
	    {

	    	/* New curve at current level; DIRICHLET or PASSIVE */

	    	wave_type(*c) = DIRICHLET_BOUNDARY; /* Default setting */
	    	for (i = 0; i < dim; i++)
	    	    coords[i] = 0.5*( Coords((*c)->first->start)[i] +
						Coords((*c)->first->end)[i] );
		transform_coord(coords,chart,coords_par,chart->parent);
		comp_par = component(coords_par,intfc_par);
		if (nearest_interface_point(coords_par,comp_par,intfc_par,INCLUDE_BOUNDARIES,
					    NULL,coords_par,t,
					    &hse_par,&hs_par) != YES)
		{
		    screen("ERROR in set_wave_node_types(), "
			   "nearest_interface_point failed\n");
		    clean_up(ERROR);
		}
		c_par = Curve_of_hs(hs_par);
		if (wave_type(c_par) == PASSIVE_BOUNDARY)
		{
		    wave_type(*c) = PASSIVE_BOUNDARY;
		    continue;
		}

		/* Check all parent curves around comp for PASSIVE */

		if (comp_par == positive_component(c_par))
		{
		    tempc_par = adjacent_curve(c_par,NEGATIVE_ORIENTATION,
					       COUNTER_CLOCK,&tempc_orient);
		}
		else
		{
			tempc_par = adjacent_curve(c_par,POSITIVE_ORIENTATION,
					           COUNTER_CLOCK,&tempc_orient);
		}
		if (wave_type(tempc_par) == PASSIVE_BOUNDARY)
		{
		    wave_type(*c) = PASSIVE_BOUNDARY;
		    continue;
		}
		while (tempc_par != c_par)
		{
		    if (comp_par == positive_component(tempc_par))
		    {
		    	tempc_par = adjacent_curve(tempc_par,
						   NEGATIVE_ORIENTATION,
						   COUNTER_CLOCK,
					           &tempc_orient);
		    }
		    else
		    {
		    	tempc_par = adjacent_curve(tempc_par,
				                   POSITIVE_ORIENTATION,
						   COUNTER_CLOCK,
					           &tempc_orient);
		    }
		    if (wave_type(tempc_par) == PASSIVE_BOUNDARY)
		    {
		        wave_type(*c) = PASSIVE_BOUNDARY;
		        break;
		    }
		}
	    }
	}
	for (n = intfc->nodes; *n; n++)
	{
	    if (node_type(*n) == ERROR)
	    {

	    	/* New node at current level; DIRICHLET, NEUMANN, *
	    	*  PASSIVE or FIXED */

	    	node_type(*n) = PASSIVE_NODE;
	    	fixed_node = YES;
	    	for (c = (*n)->in_curves; c && *c; c++)
	    	{
	    	    if (wave_type(*c) == DIRICHLET_BOUNDARY)
	    	    	node_type(*n) = DIRICHLET_NODE;
	    	    if (wave_type(*c) == NEUMANN_BOUNDARY)
	    	    	node_type(*n) = NEUMANN_NODE;
	    	    if (is_subdomain_boundary(Hyper_surf(*c)) )
	    	    	node_type(*n) = SUBDOMAIN_NODE;
	    	    if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
	    	    	fixed_node = NO;
	    	}
	    	for (c = (*n)->out_curves; c && *c; c++)
	    	{
	    	    if (wave_type(*c) == DIRICHLET_BOUNDARY)
	    	    	node_type(*n) = DIRICHLET_NODE;
	    	    if (wave_type(*c) == NEUMANN_BOUNDARY)
	    	    	node_type(*n) = NEUMANN_NODE;
	    	    if (is_subdomain_boundary(Hyper_surf(*c)) )
	    	    	node_type(*n) = SUBDOMAIN_NODE;
	    	    if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
	    	    	fixed_node = NO;
	    	}
	    	if (fixed_node)
	    	    node_type(*n) = FIXED_NODE;
	    }
	}
	if (chart->front->phys_set_node_types != NULL)
		(*chart->front->phys_set_node_types)(chart->front);
	debug_print("chart","Leaving set_wave_node_types\n");
#endif /* defined(TWOD) */
}		/*end set_wave_node_types*/

/* ARGSUSED */
LOCAL void set_zoomed_front_states(
	CHART		*chart)
{
#if defined(TWOD)

	BOND		*b;
	CURVE		**c;
	COMPONENT	lcomp,rcomp;
	Front		*front = chart->front;
	Wave		*wave = chart->wave;

	debug_print("chart","Entering set_zoomed_front_states()\n");

	if (front->sizest == 0)
	{
	    debug_print("chart","Leaving set_zoomed_front_states()\n");
	    return;
	}

	for (c = front->interf->curves;*c;c++)
	{

	    	/* Find components in chart */

	    rcomp = positive_component((*c));
	    lcomp = negative_component((*c));

	    	/* Find states in best of coarser or all charts */

	    	/* take care of start of curve */

	    init_amr_hyp_solution(Coords((*c)->start->posn),lcomp,
			          Hyper_surf(*c),NEGATIVE_SIDE,front,wave,
			          left_start_state(*c));
	    init_amr_hyp_solution(Coords((*c)->start->posn),rcomp,
			          Hyper_surf(*c),POSITIVE_SIDE,
				  front,wave,right_start_state(*c));

	    for (b = (*c)->first; b != (*c)->last; b = b->next)
	    {

		init_amr_hyp_solution(Coords(b->end),lcomp,Hyper_surf(*c),
				      NEGATIVE_SIDE,front,wave,
				      left_state(b->end));
		init_amr_hyp_solution(Coords(b->end),rcomp,Hyper_surf(*c),
				      POSITIVE_SIDE,front,wave,
				      right_state(b->end));

	    }
		 
		 /* repeat at curve end */
		
	    init_amr_hyp_solution(Coords((*c)->end->posn),lcomp,
				  Hyper_surf(*c),NEGATIVE_SIDE,front,
				  wave,left_end_state(*c));
	    init_amr_hyp_solution(Coords((*c)->end->posn),rcomp,
				  Hyper_surf(*c),POSITIVE_SIDE,front,
				  wave,right_end_state(*c));

	}
	if (debugging("s_z_front_states"))
	{
	    (void) printf("States on Zoomed Front\n");
	    graph_front_states(front);
	}
	debug_print("chart","Leaving set_zoomed_front_states()\n");
#endif /* defined(TWOD) */
}		/*end set_zoomed_front_states*/



LOCAL void set_zoomed_interior_states(
	CHART		*chart)
{
	COMPONENT	comp;
	Locstate	state;
	Front		*front = chart->front;
	Wave		*wave = chart->wave;
	double		*coords;
	int		ix,iy;
	int		icoords[MAXD];
	int		xmax = wave->rect_grid->gmax[0];
	int		ymax = wave->rect_grid->gmax[1];
	int		status;

	debug_print("chart","Entering set_zoomed_interior_states()\n");

	if (wave->sizest == 0)
	{
	    debug_print("chart","Leaving set_zoomed_interior_states()\n");
	    return;
	}
	status = init_hyp_solution_function(wave,front);
	if (status != GOOD_STEP)
	{
	    screen("ERROR in set_zoomed_interior_states(), "
		   "init_hyp_solution() failed\n");
	    clean_up(ERROR);
	}

	for (iy = 0; iy < ymax; iy++)
	{
	    icoords[1] = iy;
	    for (ix = 0; ix < xmax; ix++)
	    {
	    	icoords[0] = ix;
	    	coords = Rect_coords(icoords,wave);
	    	comp = Rect_comp(icoords,wave);
	    	state = Rect_state(icoords,wave);
	    	init_amr_hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
	    			      chart->front,chart->wave,state);
	    }
	}
	debug_print("init","Leaving set_zoomed_interior_states()\n");
}		/*end set_zoomed_interior_states*/

EXPORT	void reinit_amr_hyp_solution_functions(
	LEVEL	*level)
{
	CHART	*chart;
	
	if (level == NULL)
	    return;
	for (chart=level->first; chart != NULL; chart = chart->next_chart)
	{
	    reinit_hyp_solution_function(chart->wave,chart->front);
	    reinit_amr_hyp_solution_functions(level->next_finer_level);
	}
}		/*end reinit_amr_hyp_solution_functions*/

EXPORT	int	amr(
	CHART		*root,
	LEVEL		*level,
	Printplot	*prt,
	int		initial_step,
	double		Time,
	double		Dt)
{
	CHART		*chart;
	int		status;
	double		dt_frac;

	if (level == NULL)
	    return GOOD_STEP;
	level = level->next_finer_level;
	level->dt = Dt / level->refinement_ratio;
	if(level->step == initial_step)
	    level->dt = Dt;

			/* ----- put regrid here -------*/

	for ( ; ; ) 			/* time loop at level */
	{
	    chart=level->first;
	    for ( ; ; )		/* loop over charts at level */
	    {
	        status = time_step(chart,level->dt,&dt_frac,
		    		   level->step,level->t);
	        if (status != GOOD_STEP)
		    return status;

	        if (chart == level->last)
		    break;
	        chart = chart->next_chart;
	    }

	    level->t += level->dt;

	    if (level->next_finer_level != NULL)
	    {
	    	status = amr(root,level,prt,initial_step,level->t,level->dt);
	        if (status != GOOD_STEP)
		    return status;

	    	interpolate_from_finer_grids_to_level(level);
	    }

	    if (level->next_finer_level == NULL)
	    {
	    	if (level->step != initial_step)
	    	{
	    	    print_storage("before ELLIP","EL_storage");
		    start_clock("ELLIP");	/* Elliptic Step */
		    (*root->ellip)(NULL,root,NULL);
		    stop_clock("ELLIP");
		    print_storage("after ELLIP","EL_storage");

		}
	    }
	    level->step++;
	    if (level->t >= Time)
		break;
	}
	set_current_chart(root);
	return GOOD_STEP;
}		/*end amr*/



/*
*				hyp_amr():
*
*	This is the control flow for a single coarse grid time step in
*	automatic mesh refinement.
*/

/*ARGSUSED*/
EXPORT int hyp_amr(
	double		dt,
	double		*dt_frac,
	Wave		*wave,
	Front		*front)
{
	LEVEL		*lev;
	CHART		*root;
	int		num_coarse_steps;
	int		status;

		/* initialize dt and wave speeds */

	debug_print("chart","Entered hyp_amr()\n");
	root = chart_of_front(front);
	for (lev = root->level; ; )
	{
	    if (lev->level_num != 0)
	    	lev->dt = lev->next_coarser_level->dt / lev->refinement_ratio;
	    else
	    	lev->dt = dt;
	    if (lev->level_num == *lev->max_num_levels)
		break;
	    lev = lev->next_finer_level;
	}
	initialize_max_front_speed(root->level->first->front);

		/* Integrate one step on coarse grid */

	for (num_coarse_steps = 0; num_coarse_steps <= 2;)
	{
	    find_next_level(root,&lev);
	    if (lev->level_num == 0)
	    	num_coarse_steps++;
	    if (num_coarse_steps == 2)
	    	break;
	    start_clock("integrate_at_level");
	    status = integrate_at_level(root,lev,dt_frac);
	    if (status != GOOD_STEP)
	    	return status;
	    stop_clock("integrate_at_level");
	    if (!lev->steps_before_regrid)
	    	regrid_at_level(lev);
	    if (!lev->steps_before_synchronize)
		interpolate_from_finer_grids_to_level(lev);
	}
	debug_print("chart","Leaving hyp_amr()\n");
	return status;
}		/*end hyp_amr*/


/*
*				find_next_level():
*
*	Used in the control flow for a single coarse grid time step. This
*	routine finds the next grid level to be integrated. The algorithm
*	is to choose the coarsest level with the property that all finer
*	levels are at the same time.  If steps_before_refinement ==
*	refinement_ratio, then a level and its parent are at the same time.
*	Thus a fine grid will never get ahead in time of some coarser grid.
*/


LOCAL void find_next_level(
	CHART		*root,
	LEVEL		**level)
{
	LEVEL		*lev = root->level;
	int		num = *root->level->max_num_levels;
	
		/* Find finest level */

	while (num--)
		lev = lev->next_finer_level;

		/* Are finest level and parent at same time? */

	if (lev->steps_before_level_change != lev->refinement_ratio)
	{
		*level = lev;
		return;
	}

		/* Find coarsest common time level */

	while (lev->next_coarser_level && lev->step == 
		lev->refinement_ratio * lev->next_coarser_level->step)
		lev = lev->next_coarser_level;
	*level = lev;
}		/*end find_next_level*/

/*
*				integrate_at_level():
*/

LOCAL int integrate_at_level(
	CHART		*root,
	LEVEL		*level,
	double		*dt_frac)
{
	CHART		*chart;
	LEVEL		*lev;
	boolean		integrate = YES;
	int		status;

		/* Dynamics on given grid level */

	debug_print("chart","Entered integrate_at_level()\n");
	if (debugging("chart"))
		print_LEVEL(stdout,level);
	for (chart = level->first; chart; chart = chart->next_chart)
	{
		if (debugging("chart"))
			print_CHART(stdout,chart);
		chart->front->step = level->step;
		chart->front->time = level->t;

		status = (*chart->hyp_solver)(level->dt,dt_frac,
				      chart->wave,chart->front);

		stop_clock("hyp_solver");
		if (status != GOOD_STEP)
		{
		       (void) printf("WARNING - hyp step failed, status = %d\n",
				     status);
		       return status;
		}

		if (level->level_num != 0)
		{
			chart->front->_hyp_solution = amr_hyp_solution;
			chart->front->_hyp_grad_solution = 
						amr_hyp_grad_solution;
		}
		start_clock("hyp_solver");
	}

		/* Reset step counting flags and times */

	level->t += level->dt;
	level->step++;
	level->steps_before_synchronize--;
	level->steps_before_regrid--;
	level->steps_before_level_change--;


		/* assign newwave and newfront? */

	if (level->level_num != *level->max_num_levels) return status;
	lev = level;
	while (integrate == YES)
	{
	    for (chart = lev->first; chart; chart = chart->next_chart)
	    {
	    	assign_front_interface(chart->front,chart->newfront);
	    	assign_wave_pointers(chart->wave,chart->newwave);
	    	include_max_front_speedInfo(MaxFrontSpeed(chart->newfront),
					     root->front);
		include_max_wave_speed_info(MaxWaveSpeed(chart->newwave),
					    root->wave);
	    }
	    if (lev->level_num == 0)
	    {
	    	lev->steps_before_level_change = 1;
	    	break;
	    }
	    if (lev->steps_before_level_change != 0) break;
	    lev->steps_before_level_change = lev->refinement_ratio;

	    	/* Set level, first to next coarser level */

	    lev = lev->next_coarser_level;
	}
	debug_print("chart","Leaving integrate_at_level after assign fr, wave\n");
	return status;
}		/*end integrate_at_level*/



/*
*				regrid_at_level():
*
*	Calls routines to flag and cluster bad points, to determine new chart
*	rectangle locations.
*/

LOCAL void regrid_at_level(
	LEVEL		*level)
{
		level->steps_before_regrid = level->steps_between_regrid;
	/* Not needed for static charts */
}		/*end regrid_at_level*/

EXPORT void interpolate_from_finer_grids_to_level(
	LEVEL		*level)
{
	BOND		*b;
	CHART		*chart;
	COMPONENT	comp;
	CURVE		**c;
	INTERFACE	*intfc;
	Locstate	state;
	Wave		*wave;
	double		*coords;
	int		ix,iy,xmax,ymax;
	int		icoords[MAXD];

		/* Interpolate states */

	debug_print("chart","Entered interpolate_from_finer_grids_to_level\n");
	for (chart = level->first; chart; chart = chart->next_chart)
	{ 
		intfc = chart->front->interf;
		for (c = intfc->curves; *c; c++)
		{

			/* Start and end of curve */
			/* Is this good enough for nodes???? */

			intfc_point_and_states_from_finer_levels(
				(*c)->start->posn,(*c)->first,*c,chart);
			intfc_point_and_states_from_finer_levels(
				(*c)->end->posn,(*c)->last,*c,chart);

			/* Interior points of curve */

			for (b = (*c)->first; b->next; b = b->next)
				intfc_point_and_states_from_finer_levels(
					b->end,b,*c,chart);
		}
		/* Interpolate states at regular grid points */

		wave = chart->wave;
		if (wave->sizest == 0) break;
		xmax = wave->rect_grid->gmax[0];
		ymax = wave->rect_grid->gmax[1];
		for (iy = 0; iy < ymax; iy++)
		{
			icoords[1] = iy;
			for (ix = 0; ix < xmax; ix++)
			{
				icoords[0] = ix;
				coords = Rect_coords(icoords,wave);
				comp = Rect_comp(icoords,wave);
				state = Rect_state(icoords,wave);
				solution_from_finer_levels(coords,comp,
					chart,state);
			}
		}
	}
	level->steps_before_synchronize = level->steps_between_synchronize;
	debug_print("chart","Leaving interpolate_from_finer_grids_to_level\n");
}		/*end interpolate_from_finer_grids_to_level*/

enum {
	SEARCH_FINER_LEVELS = 1,
	SEARCH_COARSER_LEVELS,
	SEARCH_ALL_LEVELS
};

enum {
	STRICT_INTERIOR = 1,
	INTERIOR,
	EXTERIOR
};

/*
*				init_amr_hyp_solution():
*
*	This is the amr version of hyp_solution used only for INITIALIZATION of
*       new grid. We search over charts,
*	to find the best chart for evaluation of the ordinary solution
*	function, based on the chart found. The search is over the 
*	coarser levels or all levels (except the one initialized).
*	The search for all grids is used for regridding at time step != 0.
*	A finest chart with the given point one mesh 
*	block from the border is selected as the best.  In case there is no
*	such chart, the point is either very close to the border of all
*	charts, including the root chart, or is outside the rectangle in the
*	root chart.  In the first of these two cases, the best chart is the
*	one with the smallest distance from the given point to the rectangle
*	one mesh block in from the border.  In the second case, hyp_solution
*	on the root level is called, which, depending on the boundary 
*	conditions at an
*	obstacle state, or stored ambient (Dirichlet) data.  
*/
/*ARGSUSED*/
LOCAL void init_amr_hyp_solution(
	double		*coords,
	COMPONENT	comp,
	HYPER_SURF	*hs,
	SIDE		side,
	Front		*front,
	Wave		*wave,
	Locstate	ans)
{
	CHART		*chart,	*best_chart;
	COMPONENT	best_comp;
	double		best_coords[MAXD];
	int		status;
	static Locstate	ans1 = NULL;
	
	if (ans1 == NULL)
	{
	    alloc_state(front->interf,&ans1,front->sizest);
	}
	if ((chart = chart_of_front(front)) == NULL)
	{
	    screen("ERROR in init_amr_hyp_solution(), NULL chart \n");
	    clean_up(ERROR);
	}
	if (chart->level->level_num == 0)
	{
	    hyp_solution(coords,comp,hs,side,chart->front,
			 chart->wave,ans,NULL);
	    return;
	}
		
	if (states_from_coarser_or_best_chart == FROM_COARSER)
	    status = interior_coord_in_best_chart(coords,chart,best_coords,
			                          &best_chart,
						  SEARCH_COARSER_LEVELS);
	else
	    status = interior_coord_in_best_chart(coords,chart,best_coords,
			                          &best_chart,
						  SEARCH_ALL_LEVELS);
	if (chart == best_chart || is_interior_comp(comp,front->interf)) 
	    best_comp = comp;
	else if (hs == NULL)
	{
	    if (status == EXTERIOR)
	    	best_comp = exterior_component(best_chart->front->interf);
	    else
	    	best_comp = component(best_coords,best_chart->front->interf);
	}
	else  
	    best_comp = comp_of_bdry(hs,side,chart,best_chart);

	if (chart == best_chart)
	    hyp_solution(best_coords,best_comp,hs,side,
			 best_chart->front,best_chart->wave,ans,NULL);
	else
	    hyp_solution(best_coords,best_comp,NULL,UNKNOWN_SIDE,
			 best_chart->front,best_chart->wave,ans,NULL);

		/* Transform data back to chart */

	(*best_chart->front->transform_state)(ans,best_chart->to_root);
	(*chart->front->transform_state)(ans,chart->from_root);
}		/*end init_amr_hyp_solution*/

/*
*				amr_hyp_solution():
*
*	This is the amr version of hyp_solution, with a search over charts,
*	to find the best chart for evaluation of the ordinary solution
*	function, based on the chart found. The search is over the given 
*	and coarser levels. A finest chart with the given point one mesh 
*	block from the border is selected as the best.  In case there is no
*	such chart, the point is either very close to the border of all
*	charts, including the root chart, or is outside the rectangle in the
*	root chart.  In the first of these two cases, the best chart is the
*	one with the smallest distance from the given point to the rectangle
*	one mesh block in from the border.  In the second case, hyp_solution
*	on the root level is called, which, depending on the boundary 
*	conditions, returns data at an
*	obstacle state, or stored ambient (Dirichlet) data.  A similar 
*	routine, with the search limited to finer levels, is used for
*	resynchronize.
*/

/*ARGSUSED*/

EXPORT void amr_hyp_solution(
	double		*coords,
	COMPONENT	comp,
	HYPER_SURF	*hs,
	SIDE		side,
	Front		*front,
	POINTER		wv,
	Locstate	ans,
	Locstate	dflt_ans)
{
	CHART		*chart,*best_chart;
	COMPONENT	best_comp;
	CURVE		*curve = Curve_of_hs(hs);
	double		best_coords[MAXD];
	double		a,b,c;
	int		status;
	static Locstate ans1 = NULL;
	Wave		*wave = (Wave*)wv;
	
	if (ans1 == NULL)
	{
	    alloc_state(front->interf,&ans1,front->sizest);
	}
	if ((chart = chart_of_front(front)) == NULL)
	{
	    screen("ERROR in amr_hyp_solution(), NULL chart \n");
	    clean_up(ERROR);
	}
	if (states_from_coarser_or_best_chart == FROM_COARSER)
	    status = interior_coord_in_best_chart(coords,chart,best_coords,
			                          &best_chart,
						  SEARCH_COARSER_LEVELS);
	else
	    status = interior_coord_in_best_chart(coords,chart,best_coords,
			                          &best_chart,
						  SEARCH_ALL_LEVELS);
	if (chart == best_chart || is_interior_comp(comp,front->interf))
	    best_comp = comp;
	else if (curve == NULL)
	{
	    if (status == EXTERIOR)
	    	best_comp = exterior_component(best_chart->front->interf);
	    else
		best_comp = component(best_coords,best_chart->front->interf);
	}
	else 
	    best_comp = comp_of_bdry(hs,side,chart,best_chart);
	if (chart == best_chart)
	    hyp_solution(best_coords,best_comp,hs,side,best_chart->front,
			 best_chart->wave,ans,dflt_ans);
	else
	    hyp_solution(best_coords,best_comp,NULL,UNKNOWN_SIDE,
			 best_chart->front,best_chart->wave,ans,dflt_ans);
	if (chart->level->level_num != best_chart->level->level_num)
	{

		/* Interpolate in time */

	    hyp_solution(best_coords,best_comp,NULL,UNKNOWN_SIDE,
			 best_chart->newfront,best_chart->newwave,
			 ans1,dflt_ans);
	    a = (chart->level->t - best_chart->level->t);
	    b = best_chart->level->dt;
	    if (b > 1.e-12 && b > .01 * a) /*TOLERANCE*/
	    	c = a/b;
	    else c = 0.;
		interpolate_states(best_chart->front,1.-c,c,
				   best_coords,ans,best_coords,ans1,ans);
	}

		/* Transform data back to chart */

	(*best_chart->front->transform_state)(ans,best_chart->to_root);
	(*chart->front->transform_state)(ans,chart->from_root);
}		/*end amr_hyp_solution*/

/*
*		amr_hyp_grad_solution():
*/

/*ARGSUSED*/
EXPORT void amr_hyp_grad_solution(
	double		*coords,
	COMPONENT	comp,
	HYPER_SURF	*hs,
	SIDE		side,
	Front		*front,
	POINTER		wave,
	Locstate	*grad_state)
{
	screen("ERROR: Code Needed in amr_hyp_grad_solution()\n");
	clean_up(ERROR);
}		/*end amr_hyp_grad_solution*/

LOCAL void solution_from_finer_levels(
	double		*coords,
	COMPONENT	comp,
	CHART		*chart,
	Locstate	ans)
{
	CHART		*best_chart;
	double		best_coords[MAXD];
	int		status;
	
	status = interior_coord_in_best_chart(coords,chart,best_coords,
		                              &best_chart,SEARCH_FINER_LEVELS);
	if (status == EXTERIOR)
	    return;
	if (is_exterior_comp(comp,chart->front->interf))
	    comp = component(best_coords,best_chart->front->interf);
	hyp_solution(best_coords,comp,NULL,UNKNOWN_SIDE,
		     best_chart->front,best_chart->wave,ans,NULL);

		/* Transform data back to chart */

	(*best_chart->front->transform_state)(ans,best_chart->to_root);
	(*chart->front->transform_state)(ans,chart->from_root);
}		/*end solution_from_finer_levels*/

/* ARGSUSED */
LOCAL void intfc_point_and_states_from_finer_levels(
	POINT		*p,
	BOND		*b,
	CURVE		*c,
	CHART		*chart)
{
#if defined(TWOD)
	CHART		*best_chart;
	COMPONENT	lcomp,rcomp,comp;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF	*hs;
	double		best_coords[MAXD];
	double		t[MAXD];
	int		status;
	Locstate	lstate,rstate;

	if ( is_subdomain_boundary(Hyper_surf(c)) )
	    return;

	status = interior_coord_in_best_chart(Coords(p),chart,best_coords,
		                              &best_chart,SEARCH_FINER_LEVELS);
	if (status == EXTERIOR || status == INTERIOR)
	    return;

		/* Now status == STRICT_INTERIOR */

	lcomp = negative_component(c);
	rcomp = positive_component(c);
	comp = (is_exterior_comp(lcomp,chart->front->interf)) ? rcomp : lcomp;

		/* Set coord of p from coord on best_chart */

	if (nearest_interface_point(best_coords,comp,chart->front->interf,NO_BOUNDARIES,
				    NULL,best_coords,t,&hse,&hs) != YES)
	{
	    screen("ERROR in intfc_point_and_states_from_finer_levels(), "
		   "nearest_interface_point failed\n");
	    clean_up(ERROR);
	}
	transform_coord(best_coords,best_chart,Coords(p),chart);
	if (chart->front->sizest == 0)
	    return;

	lstate = left_state_at_point_on_curve(p,b,c);
	rstate = right_state_at_point_on_curve(p,b,c);
	if (is_interior_comp(rcomp,chart->front->interf))
	{
	    hyp_solution(best_coords,rcomp,NULL,UNKNOWN_SIDE,
	    	         best_chart->front,best_chart->wave,rstate,NULL);
	    (*best_chart->front->transform_state)(rstate,best_chart->to_root);
	    (*chart->front->transform_state)(rstate,chart->from_root);
	    if (wave_type(c)==PASSIVE_BOUNDARY ||
		wave_type(c)==NEUMANN_BOUNDARY)
		return;
	    else if (wave_type(c) == DIRICHLET_BOUNDARY) 
	    {
	    	obstacle_state(chart->front->interf,lstate,
			       chart->front->sizest);
	    	return;
	    }
	}
	if (is_interior_comp(lcomp,chart->front->interf))
	{
	    hyp_solution(best_coords,lcomp,NULL,UNKNOWN_SIDE,
			 best_chart->front,best_chart->wave,lstate,NULL);
	    (*best_chart->front->transform_state)(lstate,best_chart->to_root);
	    (*chart->front->transform_state)(lstate,chart->from_root);
	    if (wave_type(c)==PASSIVE_BOUNDARY || 
		wave_type(c)==NEUMANN_BOUNDARY)
		return;
	    else if (wave_type(c) == DIRICHLET_BOUNDARY)
	    {
	    	obstacle_state(chart->front->interf,rstate,
			       chart->front->sizest);
	    	return;
	    }
	}
#endif /* defined(TWOD) */
}		/*end intfc_point_and_states_from_finer_levels*/




/*
*			interior_coord_in_best_chart():
*
*	Returns STRICT_INTERIOR if strictly interior coordinates are found;
*	if so they are	returned from the finest possible chart. Here strictly 
*	interior means at least one mesh block from the border.  If the point is
*	interior but one mesh block from the border, then the chart which 
*	minimizes the *	distance to this region at least one mesh block from the
*	border is chosen, and INTERIOR is returned also.  If the point is not 
*	interior, then EXTERIOR is returned. The chart and chart coordinates for
*	the best chart are returned in all cases. In the purely exterior case 
*	the best chart is the one minimizing the above distance to the strict 
*	interior region. If the flag finer_or_coarser == SEARCH_FINER_LEVELS,
*	then the search is over all levels at the level of chart1 or finer.
*	Otherwise it is over all levels at the level of chart1 or coarser.
*/


LOCAL int interior_coord_in_best_chart(
	double		*coords,
	CHART		*chart1,
	double		*coords2,
	CHART		**chart2,
	int		finer_or_coarser_or_all)
{
	CHART		*ch;
	CHART		*best_int_ch = NULL;
	CHART		*best_ext_ch = NULL;
	LEVEL 		*lev = chart1->level;	/* default for coarse
						   level search */
	double		dist;			/* dist to strict interior
						   in current chart */
	double		best_int_dist = 0.;	/* dist to strict interior
						   in best int chart */
	double		best_ext_dist = 0.;	/* dist to strict interior
						   in best ext chart */
	double		best_int_x,best_int_y;
	double		best_ext_x,best_ext_y;
	int		coarsest_level_num = 0;	/* default for coarse or
						   all level search */
	int		status = ERROR;

	if (finer_or_coarser_or_all == SEARCH_FINER_LEVELS)
	{

		/* set limits for fine level search */

		coarsest_level_num = lev->level_num;
		for ( ; lev->level_num < *lev->max_num_levels; )
			lev = lev->next_finer_level;
	}

	if (finer_or_coarser_or_all == SEARCH_ALL_LEVELS)
	{

		/* set limits for all level search */

		for ( ; lev->level_num < *lev->max_num_levels; )
			lev = lev->next_finer_level;
	}
	while (lev != NULL)
	{
		for (ch = lev->first; ch; ch = ch->next_chart)
		{
			if (ch == chart1)   continue; 
		   /** Is it correct to exclude chart1 for:	 
			finer_or_coarser_or_all == SEARCH_FINER_LEVELS ? **/     
			if ( (status = is_strict_interior_point(coords,
				chart1,coords2,ch,&dist) ) == STRICT_INTERIOR )
			{
				*chart2 = ch;
				return status;
			}
			else if (status == INTERIOR &&
					(best_int_ch == NULL ||
					dist < best_int_dist)         )
			{
				best_int_ch = ch;
				best_int_x = coords2[0];
				best_int_y = coords2[1];
				best_int_dist = dist;
			}

				/* Exterior */

			else if (best_ext_ch == NULL || dist < best_ext_dist)
			{	
				best_ext_ch = ch;
				best_ext_x = coords2[0];
				best_ext_y = coords2[1];
				best_ext_dist = dist;
				
			}
		}


		if (lev->level_num == coarsest_level_num)
		{

			/* If interior at any level */

			if (best_int_ch)
			{
				*chart2 = best_int_ch;
				coords2[0] = best_int_x;
				coords2[1] = best_int_y;
				return INTERIOR;
			}

			/* Exterior at all levels searched */

			else 
			{
				*chart2 = best_ext_ch;
				coords2[0] = best_ext_x;
				coords2[1] = best_ext_y;
				return EXTERIOR;
			}
		}
		lev = lev->next_coarser_level;
	}
	return status;
}		/*end interior_coord_in_best_chart*/

/*
*			is_strict_interior_point():
*
*	STRICT_INTERIOR means at least one mesh block from the boundary.  
*	INTERIOR and EXTERIOR have the obvious meanings.  This
*	routine returns the appropriate one of these values, and finds the 
*	coodinates in the
*	new chart in any case.  The distance is the distance to the strict
*	interior region, measured and the sum of the x and y coordinate 
*	distances.
*/

LOCAL int is_strict_interior_point(
	double		*coords1,
	CHART		*chart1,
	double		*coords2,
	CHART		*chart2,
	double		*dist)
{
	RECT_GRID	*gr2 = chart2->front->rect_grid;
	double		lower_remainder;
	double		upper_remainder;
	int		ans = STRICT_INTERIOR;
	int		dim = chart2->front->interf->dim,i;

	*dist = 0.;
	transform_coord(coords1,chart1,coords2,chart2);
	for (i = 0; i < dim; i++)
	{
		lower_remainder = gr2->L[i] + cell_width(0,i,gr2) - coords2[i];
		upper_remainder = coords2[i] - gr2->U[i] +
					cell_width(gr2->gmax[i]-1,i,gr2);
		if (lower_remainder > 0.)
		{
			if (lower_remainder > gr2->h[i])
				ans = EXTERIOR;
			else
				ans = INTERIOR;
			*dist += lower_remainder;
		}
		else if (upper_remainder > 0.)
		{
			if (upper_remainder > gr2->h[i])
				ans = EXTERIOR;
			else if (ans != EXTERIOR)
				ans = INTERIOR;
			*dist += upper_remainder;
		}
	}
	return ans;
}		/*end is_strict_interior_point*/

/*
*			comp_of_bdry():
*
*	Finds the component in a newchart corresponding to a boundary curve
*	and side in another chart.
*/

/* ARGSUSED */


LOCAL COMPONENT comp_of_bdry(
	HYPER_SURF	*hs,
	int		side,
	CHART		*chart,
	CHART		*newchart)
{
	BOND	    *b = Curve_of_hs(hs)->first;
	RECT_GRID   *newgr = newchart->front->rect_grid;
	double	    coords[MAXD];
	double	    new_coords[MAXD];
	int	    i, dim = chart->front->rect_grid->dim;
	static const double MLEN = 0.001;	/* TOLERANCE */

	for (i = 0; i < dim; i++)
		coords[i] = 0.5*(Coords(b->end)[i] + Coords(b->start)[i]);
	transform_coord(coords,chart,new_coords,newchart);
	if (new_coords[0] < newgr->L[0] + MLEN * (newgr->U[0] - newgr->L[0]) ||
	    new_coords[0] > newgr->U[0] - MLEN * (newgr->U[0] - newgr->L[0]) ||
	    new_coords[1] < newgr->L[1] + MLEN * (newgr->U[1] - newgr->L[1]) ||
	    new_coords[1] > newgr->U[1] - MLEN * (newgr->U[1] - newgr->L[1])  )
			return exterior_component(newchart->front->interf);
	else
		return component(new_coords,newchart->front->interf);
}		/*end comp_of_bdry*/


LOCAL void transform_coord(
	double		*coords,
	CHART		*from,
	double		*coords_ans,
	CHART		*to)
{
	AFLIN		*to_root = from->to_root;
	AFLIN		*from_root = to->from_root;
	int		dim = from->front->rect_grid->dim;
	int		i, j;
	double		root[MAXD];

	for (i = 0; i < dim; i++)
	{
		root[i] = to_root->b[i];
		for (j = 0; j < dim; j++)
			root[i] += to_root->a[i][j] * coords[j];
	}
	for (i = 0; i < dim; i++)
	{
		coords_ans[i] = from_root->b[i];
		for (j = 0; j < dim; j++)
			coords_ans[i] += from_root->a[i][j] * root[j];
	}
		
}		/*end transform_coord*/


/*
*			set_coord_transform():
*
*	The coordinate transformations act on n-tuples, i.e. coordinate 
*	functions. Xo is the the lower limit, ie the SW corner of the
*	rectangle in both the root and chart coordinate system.
*	xu,yu are the root coordinates of the NE corner. `A' represents a
*	rotation about the point Xo. In 2-D it is the angle between the 
*	coordinates, measured from the positive x root axis to the 
*	positive x chart axis.  Thus the chart to root transformation is
*	             Y = A(X - Xo) + Xo = AX + b,
*	where b = (I - A)Xo.  We assume that A is an orthogonal bi_array,
*       so inv(A) = transpose(A) and det(A)=1.
*/

LOCAL void set_coord_transform(
	CHART		*chart,
	double		*Xo,
	double		**A,
	int             dim)
{
	int		i, j;

	for (i = 0; i < dim; i++)
	{
		chart->to_root->b[i] = Xo[i];
		chart->from_root->b[i] = Xo[i];
		for (j = 0; j < dim; j++)
		{
			chart->to_root->a[i][j] = A[i][j];
			chart->from_root->a[i][j] = A[j][i];
			chart->to_root->b[i] -= A[i][j]*Xo[j];
			chart->from_root->b[i] -= A[j][i]*Xo[j];
		}
	}

	chart->to_root->det = 1.0;
	chart->from_root->det = 1.0;
}		/*end set_coords_transform*/


LOCAL AFLIN *read_print_AFLIN(
	FILE		*file,
	int		dim)
{
	int		i, j;
	char		s[20];
	AFLIN		*aflin;
	int		old_aflin;

	scalar(&aflin,sizeof(AFLIN));
	(void) fgetstring(file,"AFLIN");
	(void) fscanf(file,"%d",&old_aflin);
	for (i = 0; i < dim; i++)
	{
	    for (j = 0; j < dim; j++)
	    {
	    	(void) sprintf(s,"a[%d][%d]",i,j);
	    	(void) fgetstring(file,s);
	    	(void) fscan_float(file,&aflin->a[i][j]);
	    }
	}
	for (i = 0; i < dim; i++)
	{
	    (void) sprintf(s,"b[%d]",i);
	    (void) fgetstring(file,s);
	    (void) fscan_float(file,&aflin->b[i]);
	}
	(void) fgetstring(file,"det");
	(void) fscan_float(file,&aflin->det);
	if (debugging("chart_restart"))  print_AFLIN(stdout,aflin,dim);
	return aflin;
}		/*end read_print_AFLIN*/

EXPORT void print_CHART(
	FILE		*file,
	CHART		*chart)
{
	int		dim = chart->front->rect_grid->dim;

	if (file == NULL)
	    return;

	(void) foutput(file);
	(void) fprintf(file,"Chart\n");
	(void) fprintf(file,"\tCHART %p\n",(POINTER)chart);
	fprint_rectangular_grid(file,chart->front->rect_grid);
	(void) fprintf(file,"\tTO_ROOT %p\n",(POINTER)chart->to_root);
	print_AFLIN(file,chart->to_root,dim);
	(void) fprintf(file,"\tFROM_ROOT %p\n",(POINTER)chart->from_root);
	print_AFLIN(file,chart->from_root,dim);
	(void) fprintf(file,"\tdynamic = %d\n",chart->dynamic);
	if (chart == chart->level->last)
		(void) fprintf(file,"\tLAST_CHART_AT_LEVEL\n");
	(void) fprintf(file,"\n");
	if (debugging("chart")) 
		(void) printf("prev_chart,next_chart = %p %p\n",
			      (POINTER)chart->prev_chart,
			      (POINTER)chart->next_chart);
}		/*end print_CHART*/

EXPORT CHART *read_print_CHART(
	const IO_TYPE *io_type,
	CHART	      *root,
	int	      *last)
{
	FILE	  *file = io_type->file;
	int	  dim = root->newfront->interf->dim;
	CHART	  *chart;
	RECT_GRID *gr;
	int	  cmp;
	int	  old_chart_p,old_to_root_p,old_from_root_p;
	char	  line[2048];

	*last = NO;
	scalar(&chart,sizeof(CHART));
	chart->front = copy_front(root->front);
	chart->front->interf = NULL;
	chart->wave = copy_wave(root->wave);
	clear_wave_pointers(chart->wave);
	scalar(&gr,sizeof(RECT_GRID));
	chart->front->rect_grid = gr;
	chart->wave->rect_grid = gr;
	chart->newfront = copy_front(root->front);
	chart->newfront->interf = NULL;
	chart->newwave = copy_wave(chart->wave);
	if (next_output_line_containing_string(file,"Chart") == NULL)
		return NULL;
	else
	{
	    (void) fgetstring(file,"CHART");
	    (void) fscanf(file,"%d",&old_chart_p);
	    read_rectangular_grid(io_type,gr,YES,
	                          &root->front->rect_grid->Remap);
	    chart->wave->rect_grid = chart->front->rect_grid ;
	    (void) fgetstring(file,"TO_ROOT");
	    (void) fscanf(file,"%d",&old_to_root_p);
	    chart->to_root = read_print_AFLIN(file,dim);
	    (void) fgetstring(file,"FROM_ROOT");
	    (void) fscanf(file,"%d",&old_from_root_p);
	    chart->from_root = read_print_AFLIN(file,dim);
	    (void) fscanf(file,"%*s %*s %d\n",&chart->dynamic);

	    if (fgets(line,2046,file) == NULL) 
	    	return NULL;

	    cmp = strcmp("LAST_CHART_AT_LEVEL\n",line);
	    if ( !cmp )
	    {
	    	(void) printf("chart %p\n",(POINTER)chart);
	    	*last = YES;
	    }
	    return chart;
	}
}		/*end read_print_CHART*/

EXPORT void print_LEVEL(
	FILE		*file,
	LEVEL		*level)
{
	if (file == NULL)
	    return;
	(void) foutput(file);
	(void) fprintf(file,"Level\n");
	(void) fprintf(file,"\tLEVEL %p\n",(POINTER)level);
	(void) fprintf(file,"max_num_levels  %p %d level_num  %d\n",
		       (POINTER)level->max_num_levels,
		       *level->max_num_levels,level->level_num);

	(void) fprintf(file,"steps_before_synch steps_between_synch = %d %d\n",
	level->steps_before_synchronize,level->steps_between_synchronize);
	(void) fprintf(file,
		       "steps_before_regrid steps_between_regrid = %d %d\n",
		       level->steps_before_regrid,level->steps_between_regrid);
	(void) fprintf(file,
	      "steps_before_level_change refinement_ratio = %d %d step = %d\n",
	      level->steps_before_level_change,level->refinement_ratio,
	      level->step);
	(void) fprintf(file,"t = %g dt = %g\n",level->t,level->dt);
}		/*end print_LEVEL*/


EXPORT LEVEL *read_print_LEVEL(
	const IO_TYPE	*io_type)
{
	FILE       *file = io_type->file;
	LEVEL	   *level;
	int	   old_level_p;
	static int max_num_levels;

	scalar(&level,sizeof(LEVEL));
	if (next_output_line_containing_string(file,"Level") == NULL)
		return NULL;
	(void) fgetstring(file,"LEVEL");
	(void) fscanf(file,"%d",&old_level_p);
	(void) fgetstring(file,"max_num_levels");  
	(void) fscanf(file,"%*d %d %*s %d",&max_num_levels,&level->level_num);
	level->max_num_levels = &max_num_levels;
	*level->max_num_levels = max_num_levels;
	(void) fgetstring(file,"steps_before_synch");
	(void) fscanf(file,"%*s %*s %d %d",&level->steps_before_synchronize,
		      &level->steps_between_synchronize);
	(void) fgetstring(file,"steps_before_regrid");
	(void) fscanf(file,"%*s %*s %d %d",&level->steps_before_regrid,
		      &level->steps_between_regrid);
	(void) fgetstring(file,"steps_before_level_change");
	(void) fscanf(file,"%*s %*s %d %d %*s %*s %d",
		      &level->steps_before_level_change,
		      &level->refinement_ratio,&level->step);
	(void) fgetstring(file,"t =");
	(void) fscanf(file,"%lf %*s %*s %lf",&level->t,&level->dt);

	return level;
}		/*end read_print_LEVEL*/

LOCAL	CHART	*cur_chart = NULL;

EXPORT	CHART	*current_chart(void)
{
	return cur_chart;
}		/*end current_chart*/

EXPORT	void	set_current_chart(CHART *chart)
{
	cur_chart = chart;
}		/*end set_current_chart*/

#if defined(THREED)
LOCAL	void prompt_for_3d_rotation(
	double		**R)
{
	double		st0, ct0, sp0, cp0, r[MAXD];
	double		U[MAXD][MAXD], Q[MAXD][MAXD], A[MAXD][MAXD];
	double		angle, mr;
	int		i, j, k, dim = 3;
	const double eps = MACH_EPS;

	screen("Enter the coordinates of the axis of rotation: ");
	for (i = 0; i < dim; i++) (void) Scanf("%f",r+i);
	(void) Scanf("\n");
	mr = mag_vector(r,dim);
	for (i = 0; i < dim; i++) r[i] /= mr;
	if (fabs(r[0]) < eps && fabs(r[1]) < eps)
	{
		for (i = 0; i < dim; i++)
			for (j = 0; j < dim; j++)
				U[i][j] = 0.0;
		for (i = 0; i < dim; i++) U[i][i] = 1.0;
	}
	else
	{
		cp0 = r[2];
		sp0 = sqrt(1.0 - sqr(cp0));
		ct0 = r[0]/sp0;
		st0 = r[1]/sp0;
		U[0][0] = ct0*cp0; U[0][1] = -st0; U[0][2] = ct0*sp0;
		U[1][0] = st0*cp0; U[1][1] = ct0;  U[1][2] = st0*sp0;
		U[2][0] = -sp0;    U[2][1] = 0.0;  U[2][2] = cp0;
	}
	screen("Enter the rotation angle (degrees) about ");
	screen("the axis of rotation: ");
	(void) Scanf("%f\n",&angle);
	Q[0][0] = Q[1][1] = cos(angle);
	Q[1][0] = sin(angle);
	Q[0][1] = -Q[1][0];
	for (i = 0; i < 2; i++) Q[2][i] = Q[i][2] = 0.0;
	Q[2][2] = 1.0;
	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
		{
			A[i][j] = 0.0;
			for (k = 0; k < dim; k++)
				A[i][j] += Q[i][k]*U[j][k];
		}
	}
	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
		{
			R[i][j] = 0.0;
			for (k = 0; k < dim; k++)
				R[i][j] += U[i][k]*A[k][j];
		}
	}

}		/*end prompt_for_3d_rotation*/
#endif /* defined(THREED) */



#if DONT_COMPILE

                /* UNUSED function declarations */
LOCAL	CHART	*chart_of(Front*);
LOCAL	void	regrid_amr(CHART*);
LOCAL	void	delete_chart(CHART*);
LOCAL	void	delete_level(LEVEL*);

#pragma	noinline	chart_of
#pragma	noinline	regrid_amr
#pragma	noinline	delete_chart
#pragma	noinline	delete_level


/* ARGSUSED */
LOCAL void regrid_amr(
	CHART		*root)
{
	CHART		*chart;
	LEVEL		*level;

	debug_print("regrid","Entered regrid_amr()\n");
	screen("regridding time step %d reached\n",regrid_time_step);
	if (chart_of_front(root->front) == NULL)
	{
		(void) printf("chart_of_front(root->front) = %d\n",
			      chart_of_front(root->front));
		init_amr(root);
		return;
	}


	/* The following loop is temporary. It enables to remove 
	 * only old dynamic grids after using them to initialize new grids */
	for (level = root->level; level; level = level->next_finer_level)
	{
		(void) printf("level = %d next_finer_level = %d\n",level,
			      level->next_finer_level);
		print_LEVEL(stdout,level);
		for (chart = level->first; chart; chart =chart->next_chart)
		{
			(void) printf("chart = %d next_chart = %d\n",chart,
			chart->next_chart);
			print_CHART(stdout,chart);
			chart->is_old_chart = YES;
		}
	}

		/* add new grids */

	states_from_coarser_or_best_chart = FROM_BEST;

	(void) printf("FROM BEST\n");
	while (prompt_for_child_chart(root))
                ;
	if (debugging("regrid"))
	{
		CHART *ch;
		LEVEL *le;
		for (le = root->level; le; le = le->next_finer_level)
			for (ch = le->first;ch;ch=ch->next_chart)
				graph_front_states(ch->front);
	}

	for (level = root->level; level; level = level->next_finer_level)
	{
		for (chart = level->first; chart; chart->next_chart)
		{
			if ( chart->is_old_chart == YES && 
						chart->dynamic == YES )
			{
				delete_chart(chart);	
			}
		}
	}

	states_from_coarser_or_best_chart = FROM_COARSER;

	debug_print("chart","Leaving regrid_amr\n");
}		/*end regrid_amr*/

LOCAL void delete_chart(
	CHART		*chart)
{
	debug_print("regrid","Entering delete_chart()\n");
	if (chart->prev_chart == NULL)
	{
		/*first at level */
		if (chart->next_chart == NULL) /* also last at level */
			delete_level(chart->level);
		else
			chart->level->first = chart->next_chart;
	}
	else if (chart->next_chart == NULL) 
		/*last at level */
		chart->level->last = chart->prev_chart;
	else
	{		
		chart->prev_chart->next_chart = chart->next_chart;
		chart->next_chart->prev_chart = chart->prev_chart;
	}	
	debug_print("regrid","Leaving delete_chart()\n");
	return;
}		/*end delete_chart*/

LOCAL void delete_level(
	LEVEL		*level)
{
	/* we never delete level 0 (root level) */
	debug_print("regrid","Entering delete_level()\n");
	 
	level->next_coarser_level->next_finer_level = level->next_finer_level;

	if (level->next_finer_level != NULL) 
		level->next_finer_level->next_coarser_level =
			level->next_coarser_level;

	*level->max_num_levels -= 1;
	debug_print("regrid","Leaving delete_level()\n");
	return;
}		/*end delete_level*/

/*
*				chart_of():
*
*	Returns the chart of front. This is a temporary routine (I hope!).
*/

LOCAL CHART *Last_chart = NULL;

LOCAL CHART *chart_of(
	Front		*front)
{
	CHART		*chart;
	LEVEL		*lev;

	if (Last_chart && Last_chart->front == front) return Last_chart;
	lev = Root->level;
	while (YES)
	{
		for (chart = lev->first; chart; chart = chart->next_chart)
			if (chart->front == front)
			{
				Last_chart = chart;
				return chart;
			}
		if (lev->level_num == *lev->max_num_levels)
		{
			screen("ERROR: in chart_of()\n");
			return NULL;
		}
		lev = lev->next_finer_level;
	}
}		/*end chart_of*/

#endif /* DONT_COMPILE */
