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
*
*				guntan2d.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Resolves untangles of interior curves.
*
*/

#if defined(TWOD)

#define DEBUG_STRING    "untangle"
#include <gdecs/gdecs.h>


	/* LOCAL Function Declarations */
LOCAL	int	consistent_scalar_vector_unravel(CURVE**);
LOCAL	int	scalar_vector_unravel(Front*,CROSS*,int);
LOCAL	int	unitary_scalar_vector_unravel(Front*,CROSS*);
LOCAL	int	untangle_satisfies_CFL(NODE**,double**,Front*,int);
LOCAL	int	vector_vector_unravel(Front*,CROSS*);
LOCAL	int	wave_end_vector_vector_unravel(Front*,CROSS*);
LOCAL	void	adjust_node_velocities_for_symmetry(double**,int);
LOCAL	void	insert_3bond_wave(NODE**,CURVE***,ORIENTATION**,
				  double**,int,int);



EXPORT int g_untangle_interior_curves(
	Front		*fr,
	CROSS		**cross,
	int		flag)
{
	CROSS *fst_cr = *cross, *cr;
	NODE  *node, **n;
	CURVE **c, **delete_curves;
	int   status;

	DEBUG_ENTER(g_untangle_interior_curves)
	if (*cross == NULL)
	{
	    DEBUG_LEAVE(g_untangle_interior_curves)
	    return CURVES_UNTANGLED;
	}

	start_clock("g_untangle_interior");

	/* Untangle all vector tangles and remove from cross list */

	for (cr = *cross; cr; cr = cr->next)
	{
	    if ((is_vector_vector_cross(cr) == NO) &&
		(is_scalar_vector_cross(cr) == NO))
		continue;

#if !defined(FULL_PHYSICS)

	    screen("ERROR in g_untangle_interior_curves(), "
	           "uni_array wave interaction detected\n"
	           "Unable to handle with partial physics\n"
	           "Recompile with -DFULL_PHYSICS\n");
	    clean_up(ERROR);

#else /* !defined(FULL_PHYSICS) */

		/*
		 *  HACK: regard self intersections as scalar
		 *  TODO: remove this assumption
		 */

	    if (cr->c1 == cr->c2)
		continue;

		/* Untangle vector crosses */

	    if (DEBUG)
	    	(void) printf("Testing vector uni_array case\n");
	    if (is_vector_vector_cross(cr) == YES)
	    {
	    	status = vector_vector_unravel(fr,cr);
	    	if (status != CURVES_UNTANGLED)
	    	{
	    	    stop_clock("g_untangle_interior");
	            DEBUG_LEAVE(g_untangle_interior_curves)
	    	    return status;
	    	}
	    }
	    else
	    {
	    	status = scalar_vector_unravel(fr,cr,flag);
	    	if (status != CURVES_UNTANGLED)
	    	{
	    	    stop_clock("g_untangle_interior");
	            DEBUG_LEAVE(g_untangle_interior_curves)
	    	    return status;
	    	}
	    }

	    	/* Delete cr from cross list */

	    if (DEBUG)
	    	(void) printf("Deleting cr from cross list\n");
	    fst_cr = cr->prev ? cr->prev : cr->next;
	    delete_from_cross_list(cr);
	    while (fst_cr && fst_cr->prev)
		fst_cr = fst_cr->prev;

#endif /* !defined(FULL_PHYSICS) */

	}
	*cross = fst_cr;

	debug_front("untangle","after vector untangles",fr);

	/* Purely scalar crosses (and vector self intersections) */

	if (DEBUG)
	    (void) printf("Pure scalar unravel *cross = %p\n",*cross);
	status = scalar_unravel(fr,cross,flag);

	debug_front("untangle","after scalar untangles",fr);

		/* Delete short vector loops */
	
	delete_curves = NULL;
	for (c = fr->interf->curves; c && *c; ++c)
	{
	    if ((wave_type(*c) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE) &&
	        ((*c)->num_points <= 3) &&	is_closed_curve(*c))
	    {
	        if (!add_to_pointers(*c,&delete_curves))
		{
	            screen("ERROR in delete_subdomain_curves(), "
		           "add_to_pointers() failed\n");
		    clean_up(ERROR); }
		}
	}
	for (c = delete_curves; c && *c; ++c)
	{
	    node = (*c)->start;
	    (void) delete_curve(*c);
	    (void) delete_node(node);
	}

		/* Remove redundant CLOSED_NODE */

	for (n = fr->interf->nodes; n && *n; ++n)
	{
	    if (node_type(*n) == CLOSED_NODE)
	    {
	        if (delete_redundant_node(*n,NULL,NULL,fr))
		    n = fr->interf->nodes - 1;
	    }
	}
	stop_clock("g_untangle_interior");
	debug_front("untangle","after delete short vector loop",fr);
	DEBUG_LEAVE(g_untangle_interior_curves)
	return status;
}		/*end g_untangle_interior_curves*/


EXPORT	boolean	g_replace_unphys_loop(
	NNLIST	*nl,
	NNLIST	**new_node_list,
	CURVE	**newc,
	Front	*fr,
	int	i,
	double	min_area,
	int	flag)
{
	CURVE	**cc;
	NODE	*ns, *ne;
	boolean	status;

	debug_print("hypunr","Entered g_replace_unphys_loop()\n");

	status = f_replace_unphys_loop(nl,new_node_list,newc,
				       fr,i,min_area,flag);

	if (*newc == NULL)
	{
	    if (debugging("hypunr"))
	    {
		(void) printf("newc returned by f_replace_unphys_loop "
			      "is NULL\n");
	    }
	    debug_print("hypunr","Entered g_replace_unphys_loop()\n");
	    return status;
	}

	if ((Same_params(left_start_state(*newc),right_start_state(*newc))) &&
	     (positive_component(*newc) == negative_component(*newc)))
	{
	    if (debugging("hypunr"))
	    {
		(void) printf("deleteing newc returned by "
			      "f_replace_unphys_loop\n");
	    }
	    (void) delete_curve(*newc);
	    *newc = NULL;
	    debug_print("hypunr","Entered g_replace_unphys_loop()\n");
	    return status;
	}

	ns = (*newc)->start;		ne = (*newc)->end;
	node_type(ns) = node_type(ne) = CC_NODE;
	for (cc = ns->in_curves; cc && *cc; ++cc)
		end_status(*cc) = INCIDENT;
	for (cc = ns->out_curves; cc && *cc; ++cc)
		start_status(*cc) = INCIDENT;
	for (cc = ne->in_curves; cc && *cc; ++cc)
		end_status(*cc) = INCIDENT;
	for (cc = ne->out_curves; cc && *cc; ++cc)
		start_status(*cc) = INCIDENT;
	start_status(*newc) = end_status(*newc) = SLIP;

	if (debugging("hypunr"))
	{
	    (void) printf("newc\n");
	    print_curve(*newc);
	    (void) printf("Interface after creation of newc\n");
	    print_interface((*newc)->interface);
	}
	debug_print("hypunr","Entered g_replace_unphys_loop()\n");
	return status;
}		/*end g_replace_unphys_loop*/


EXPORT	void	g_phys_split_bdry_cross(
	CURVE		**physcurves,
	CURVE		**bdrycurves)
{
	set_incident_status(physcurves[0],NEGATIVE_ORIENTATION);
	set_incident_status(physcurves[1],POSITIVE_ORIENTATION);
	set_fixed_status(bdrycurves[0],NEGATIVE_ORIENTATION);
	set_fixed_status(bdrycurves[1],POSITIVE_ORIENTATION);
}		/*end g_phys_split_bdry_cross*/

#if defined(FULL_PHYSICS)

/*
*			vector_vector_unravel():
*/


LOCAL int vector_vector_unravel(
	Front		*fr,
	CROSS		*cross)
{
	BOND		*b1 = cross->b1,*b2 = cross->b2;
	CURVE		*c1 = cross->c1,*c2 = cross->c2;
	CROSS		*cross1;
	ORIENTATION	c1_orient, c2_orient, c11_orient, c12_orient;

	DEBUG_ENTER(vector_vector_unravel)
	if (	(b1 == c1->first && node_type(c1->start) == WAVE_END_NODE) ||
		(b1 == c1->last && node_type(c1->end) == WAVE_END_NODE))
	{
		DEBUG_LEAVE(vector_vector_unravel)
		return wave_end_vector_vector_unravel(fr,cross);
	}
	if (	(b2 == c2->first && node_type(c2->start) == WAVE_END_NODE) ||
		(b2 == c2->last && node_type(c2->end) == WAVE_END_NODE))
	{
		DEBUG_LEAVE(vector_vector_unravel)
		return wave_end_vector_vector_unravel(fr,cross);
	}

		/* Is a cross or overtake node (hopefully) 
		        Brian Wetton 30.06.89              */

	if (DEBUG)
		(void) printf("Cross or overtake node case\n");
	if (!find_companion_cross(cross, &cross1, &c1_orient, &c2_orient,
				     &c11_orient, &c12_orient))
	{
		(void) printf("WARNING in vector_vector_unravel(), ");
		(void) printf("Unable to find companion cross\n");
		DEBUG_LEAVE(vector_vector_unravel)
		return ERROR_IN_UNTANGLE;
	}

	if (DEBUG)
	{
		(void) printf("Cross: \n");		print_cross(cross);
		(void) printf("cross->c1\n");		print_curve(cross->c1);
		(void) printf("cross->c2\n");		print_curve(cross->c2);
		(void) printf("\n Cross1: \n");	print_cross(cross1);
		(void) printf("cross1->c1\n");		print_curve(cross1->c1);
		(void) printf("cross1->c2\n");		print_curve(cross1->c2);
		print_orientation("c1_orient: ", c1_orient,"\n");
		print_orientation("c2_orient: ", c2_orient,"\n");
		(void) printf("\n \n");
	}

		/* Determine if the new node is a cross node or an
		   overtake node and call appropriate routine
		   in gvecuntan.c */

	if ( (is_forward_wave(wave_type(cross->c1)) ==
		is_forward_wave(wave_type(cross->c2))) 
			==
	     (c1_orient == c2_orient) )
	{
		DEBUG_LEAVE(vector_vector_unravel)
		return vector_overtake_unravel (fr, cross, cross1,
			c1_orient, c2_orient);
	}
	else
	{
		DEBUG_LEAVE(vector_vector_unravel)
		return vector_cross_unravel (fr, cross, cross1,
			c1_orient, c2_orient);
	}
}		/*end vector_vector_unravel*/

/*
*			wave_end_vector_vector_unravel():
*/

LOCAL int wave_end_vector_vector_unravel(
	Front		*fr,
	CROSS		*cross)
{
	BOND		*b1 = cross->b1,*b2 = cross->b2;
	CURVE		*c1 = cross->c1,*c2 = cross->c2;
	CURVE		*c;
	CURVE		*presplit_c;
	CURVE		*no_split_c;
	CURVE		**curves;
	BOND		*presplit_b;
	NODE		*n;
	double		c1_dist,c2_dist;
	ORIENTATION	c1_orient,c2_orient,no_split_orient;
	boolean		sav_scss = interpolate_states_at_split_curve_node();

		/* Check for overtake node */

	DEBUG_ENTER(wave_end_vector_vector_unravel)
	if (	(b1 == c1->first && node_type(c1->start) == WAVE_END_NODE))
	    c1_orient = POSITIVE_ORIENTATION;
	else if ((b1 == c1->last && node_type(c1->end) == WAVE_END_NODE))
	    c1_orient = NEGATIVE_ORIENTATION;
	else
	    goto overtake1;

		/* default settings for overtake2 */

	presplit_c = c2;
	no_split_c = c1;
	presplit_b = b2;
	no_split_orient = c1_orient;
	if (	(b2 == c2->first && node_type(c2->start) == WAVE_END_NODE))
	    c2_orient = POSITIVE_ORIENTATION;
	else if ((b2 == c2->last && node_type(c2->end) == WAVE_END_NODE))
	    c2_orient = NEGATIVE_ORIENTATION;
	else
	    goto overtake2;
	if ((wave_type(c1) == wave_type(c2) && c1_orient == c2_orient) ||
	    (wave_type(c1) != wave_type(c2) && c1_orient != c2_orient))
	    goto overtake3;

		/* Join curves */

	n = Node_of(c1,c1_orient);
	rcl_after_delete_bond_fragment_at_node(cross,cross->p,c1,c1_orient);
	Coords(n->posn)[0] = Coords(cross->p)[0];
	Coords(n->posn)[1] = Coords(cross->p)[1];
	rcl_after_delete_bond_fragment_at_node(cross,cross->p,c1,c1_orient);
	if (c2_orient == POSITIVE_ORIENTATION)
	    c2->first->start = n->posn;
	else
	    c2->last->end = n->posn;
	if (wave_type(c1) != wave_type(c2))
	    invert_curve(c2);
	if ((negative_component(c1) != negative_component(c2)) ||
	    (positive_component(c1) != positive_component(c2)))
	{
	    screen("ERROR in wave_end_vector_vector_unravel(), "
	           "inconsistent components\n");
	    clean_up(ERROR);
	}
	if ((c = join_curves(c1,c2,negative_component(c1),
			     positive_component(c1),NULL)) == NULL)
	{
	    (void) printf("WARNING in wave_end_vector_vector_unravel(), "
	                  "join_curves() failed\n");
	    return ERROR_IN_UNTANGLE;
	}
	(void) delete_node(Node_of(c1,c1_orient));
	(void) delete_node(Node_of(c2,c2_orient));
	rcl_after_join(cross,c,c1,c2);
	if (DEBUG)
	{
	    (void) printf("Joined curves %llu %llu at (%g, %g) ",
			  curve_number(c1),curve_number(c2),
			  Coords(n->posn)[0],Coords(n->posn)[1]);
	    (void) printf("to form curve %llu\n",curve_number(c));
	    print_curve(c);
	}
	return CURVES_UNTANGLED;

		/* Insert overtake node */

overtake1:

		/* cross is in the middle of c1 */
	if ((b2 == c2->first && node_type(c2->start) == WAVE_END_NODE))
	    c2_orient = POSITIVE_ORIENTATION;
	else if ((b2 == c2->last && node_type(c2->end) == WAVE_END_NODE))
	    c2_orient = NEGATIVE_ORIENTATION;
	else
	    goto overtake4;

		/* cross is in the middle of c1 only */

	presplit_c = c1;
	no_split_c = c2;
	presplit_b = b1;
	no_split_orient = c2_orient;
overtake2:
	if (DEBUG)
	{
	    (void) printf("overtake case 1-2\n");
	    (void) printf("presplit_c\n");
	    print_curve(presplit_c);
	    (void) printf("no_split_c\n");
	    print_curve(no_split_c);
	}

		/* cross is in the middle of c2 only */


	if (insert_point_in_bond(cross->p,presplit_b,presplit_c) !=
	    FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in wave_end_vector_vector_unravel(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	rcl_after_insert_point(cross,cross->p,presplit_b);
	interpolate_intfc_states(fr->interf) = NO;
	set_interpolate_states_at_split_curve_node(NO);
	curves = split_curve(cross->p,presplit_b,presplit_c,
			     negative_component(presplit_c),
			     positive_component(presplit_c),
			     negative_component(presplit_c),
			     positive_component(presplit_c));
	set_interpolate_states_at_split_curve_node(sav_scss);
	rcl_after_split(cross,cross->p,cross->b1,cross->c1,curves);
	interpolate_intfc_states(fr->interf) = YES;
	if (DEBUG)
	{
	    (void) printf("presplit_c\n");
	    print_curve(presplit_c);
	    (void) printf("curves[0]\n");
	    print_curve(curves[0]);
	    (void) printf("curves[1]\n");
	    print_curve(curves[1]);
	    (void) printf("no_split_c\n");
	    print_curve(no_split_c);
	}
	interpolate_intfc_states(fr->interf) = YES;
	node_type(curves[0]->end) = OVERTAKE_NODE;
	n = Node_of(no_split_c,no_split_orient);
	change_node_of_curve(no_split_c,no_split_orient,curves[0]->end);
	if (no_split_orient == POSITIVE_ORIENTATION)
	{
	    start_status(no_split_c) = INCIDENT;
	    if (wave_type(no_split_c) == wave_type(presplit_c))
	    {
	    	end_status(curves[0]) = TRANSMITTED;
	    	start_status(curves[1]) = OVERTOOK;
	    }
	    else
	    {
	    	end_status(curves[0]) = OVERTOOK;
	    	start_status(curves[1]) = TRANSMITTED;
	    }
	}
	else
	{
	    end_status(no_split_c) = INCIDENT;
	    if (wave_type(no_split_c) == wave_type(presplit_c))
	    {
	    	end_status(curves[0]) = OVERTOOK;
	    	start_status(curves[1]) = TRANSMITTED;
	    }
	    else
	    {
	    	end_status(curves[0]) = TRANSMITTED;
	    	start_status(curves[1]) = OVERTOOK;
	    }
	}
	(void) delete_node(n);
	if (DEBUG)
	{
	    (void) printf("curves[0]\n");
	    print_curve(curves[0]);
	    (void) printf("curves[1]\n");
	    print_curve(curves[1]);
	    (void) printf("no_split_c\n");
	    print_curve(no_split_c);
	    print_interface(fr->interf);
	    (void) printf("return status = CURVES_UNTANGLED\n");
	}
	DEBUG_LEAVE(wave_end_vector_vector_unravel)
	return CURVES_UNTANGLED;

overtake3:
		/* c1 and c2 oppositely oriented w/r intersecting bonds */

	if (DEBUG)
	    (void) printf("overtake case 3\n");
	c1_dist = 
	    separation(cross->p,Node_of(c1,c1_orient)->posn,fr->interf->dim);
	c2_dist = 
	    separation(cross->p,Node_of(c2,c2_orient)->posn,fr->interf->dim);
	if (c1_dist > c2_dist)
	{
		
		/* split c1 in this case; else c2 */

	    presplit_c = c1;
	    no_split_c = c2;
	    presplit_b = b1;
	    no_split_orient = c2_orient;
	}
	goto overtake2;
	
overtake4:

		/* cross is in the middle of c1 and c2 */

	(void) printf("WARNING in wave_end_vector_vector_unravel(), "
	              "Code needed vector uni_array cross case overtake4\n");
	print_curve(c1);
	print_curve(c2);
	print_cross(cross);
	return ERROR_IN_UNTANGLE;
}		/*end wave_end_vector_vector_unravel*/

/*
*			scalar_vector_unravel():
*/

LOCAL int scalar_vector_unravel(
	Front	*front,
	CROSS	*cross,
	int	flag)
{
	CROSS		*cross1;
	CURVE		*vcurve,*scurve;
	CURVE		*curves1[2], *curves2[2];
	CURVE		**c;
	BOND		*bb;
	NODE		*diffn[2];
	NODE		*n, *ns, *ne;
	double		len;
	double		coords[MAXD];
	ANGLE_DIRECTION	i_to_f_dir[2];
	boolean		is_reflected_shock[2];
	int		i, j;
	boolean		sav_intrp;
	boolean		track[7];
	int		w_type[2][7];
	int		numcr;
	int		node_vel_set[2];
	int		dim = front->interf->dim;
	NODE_FLAG	ndflag;
	COMPONENT	comp[7];
	Locstate	ls, rs;
	RP_DATA		*RP[2];
	static	CURVE	    ***dncur = NULL;
	static	double	    ***t = NULL;
	static	double	    **nod_v = NULL;
	static	double	    **tgnt = NULL;
	static  Locstate    **l_st = NULL, **r_st = NULL;
	static  COMPONENT   **l_comp = NULL, **r_comp = NULL;
	static	ORIENTATION **c_orient = NULL;
	static	int	status[7] = { INCIDENT, REFLECTED, REFLECTED,
				      REFLECTED, SLIP, TRANSMITTED,
					      CONTACT_TARGET};

		/* Allocate storage */

	DEBUG_ENTER(scalar_vector_unravel)
	set_to_next_node_only(ndflag);
	if (DEBUG)
	{
	    (void) printf("Tangled interface\n");
	    print_interface(front->interf);
	}
	if (nod_v == NULL)
	{
	    bi_array(&tgnt,2,dim,FLOAT);
	    bi_array(&dncur,2,7,sizeof(CURVE *));
	    bi_array(&c_orient,2,7,sizeof(ORIENTATION));
	    tri_array(&t,2,2,2,FLOAT);
	    bi_array(&nod_v,2,SMAXD,FLOAT);
	    bi_array(&l_comp,2,7,sizeof(COMPONENT));
	    bi_array(&r_comp,2,7,sizeof(COMPONENT));
	    bi_array(&l_st,2,7,sizeof(Locstate));
	    bi_array(&r_st,2,7,sizeof(Locstate));
	}
	sav_intrp = interpolate_intfc_states(front->interf);
	interpolate_intfc_states(front->interf) = YES;
	for (i = 0; i < 2; ++i)
	for (j = 0; j < 7; ++j)
	    dncur[i][j] = NULL;

	if (wave_type(cross->c1) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE) 
	{
	    vcurve = cross->c1;
	    scurve = cross->c2;
	}
	else 
	{
	    vcurve = cross->c2;
	    scurve = cross->c1;
	    bb = cross->b1;
	    cross->c1 = vcurve;	cross->b1 = cross->b2;
	    cross->c2 = scurve;	cross->b2 = bb;
	}
	w_type[0][0] = w_type[1][0] = wave_type(vcurve);
	w_type[0][6] = w_type[1][6] = wave_type(scurve);

	    /* Delete short WAVE_END curves */

	if (node_type(vcurve->start) == WAVE_END_NODE &&
	    node_type(vcurve->end) == WAVE_END_NODE &&
	    vcurve->num_points <= 3) 
	{
	    (void) delete_curve(vcurve);
	    (void) delete_node(vcurve->start);
	    (void) delete_node(vcurve->end);
	    interpolate_intfc_states(front->interf) = sav_intrp;
	    DEBUG_LEAVE(scalar_vector_unravel)
	    return CURVES_UNTANGLED;
	}
	if (DEBUG)
	    (void) printf("Not a short WAVE_END cross\n");

	    /* find companion cross */

	if (!find_companion_cross(cross,&cross1,&c_orient[0][0],
				     &c_orient[0][6], &c_orient[1][0],
				     &c_orient[1][6]))
	{
	    int stat;
	    stat = unitary_scalar_vector_unravel(front,cross);
	    DEBUG_LEAVE(scalar_vector_unravel)
	    return stat;
	}
	if ((point_in_buffer(Coords(cross->p),front->rect_grid) == YES) && 
	    (point_in_buffer(Coords(cross1->p),front->rect_grid) == YES))
	{
	    /* This tangle will be handled by PP communication */
	    front->interf->modified = YES;
	    delete_from_cross_list(cross1);
	    DEBUG_LEAVE(scalar_vector_unravel)
	    return CURVES_UNTANGLED;
	}

	if (DEBUG) 
	{
	    (void) printf("Crosses cross and cross1:\n");
	    print_cross(cross);
	    (void) printf("cross->c1\n");
	    print_curve(cross->c1);
	    (void) printf("cross->c2\n");
	    print_curve(cross->c2);
	    print_cross(cross1);
	    (void) printf("cross1->c1\n");
	    print_curve(cross1->c1);
	    (void) printf("cross1->c2\n");
	    print_curve(cross1->c2);
	    print_orientation("c_orient[0][0] = ",c_orient[0][0],"\n");
	    print_orientation("c_orient[0][6] = ",c_orient[0][6],"\n");
	    print_orientation("c_orient[1][0] = ",c_orient[1][0],"\n");
	    print_orientation("c_orient[1][6] = ",c_orient[1][6],"\n");
	}

	    /* split curves */
	split_curves_at_cross(cross,front,&diffn[0],curves1,NULL,NULL,curves2,
			      NULL,NULL,MIN_SC_SEP(front->interf),NULL);
	split_curves_at_cross(cross1,front,&diffn[1],curves1,NULL,NULL,curves2,
			      NULL,NULL,MIN_SC_SEP(front->interf),NULL);
	delete_from_cross_list(cross1);
	if (DEBUG)
	{
	    (void) printf("after delete_from_cross_list\n");
	    for (i = 0; i < 2; ++i)
	    {
	    	if (diffn[0] == NULL)
	    	    continue;
	    	(void) printf("diffn[%d]\n",i);
    	        print_node(diffn[i]);
		(void) printf("In curves diffn[%d]\n",i);
		for (c = diffn[i]->in_curves; c && *c; ++c)
		    print_curve(*c);
		(void) printf("Out curves diffn[%d]\n",i);
		for (c = diffn[i]->out_curves; c && *c; ++c)
		    print_curve(*c);
	    }
	}

	    /* label curves, orientation and status */

	numcr = (cross1 != NULL) ? 2 : 1;
	for (i = 0; i < numcr; ++i)
	{
	    if (DEBUG)
	        (void) printf("start loop i = %d\n",i);
	    node_type(diffn[i]) = DIFFRACTION_NODE;
	    RP[i] = Rp_data(diffn[i]);

	    set_orients_and_w_type_for_shock_diffraction(w_type[i][0],
			                                 c_orient[i][0],
							 w_type[i][6],
							 c_orient[i][6],
			                                 c_orient[i],w_type[i]);

	    i_to_f_dir[i] = incident_shock_orientation(w_type[i][0],
					               c_orient[i][0]);

	    RP[i]->ang_dir = Opposite_ang_dir(i_to_f_dir[i]);
	    for (c = diffn[i]->in_curves; c && *c; ++c)
	    {
	    	if (is_scalar_wave(wave_type(*c)))
	    	{
	    	    if (c_orient[i][6] == NEGATIVE_ORIENTATION)
		    	dncur[i][6] = *c;
		    else
		    	dncur[i][4] = *c;
		}
		else
		{
		    if (c_orient[i][0] == NEGATIVE_ORIENTATION)
			dncur[i][0] = *c;
		    else
			dncur[i][5] = *c;
		}
	    }
	    for (c = diffn[i]->out_curves; c && *c; ++c)
	    {
	    	if (is_scalar_wave(wave_type(*c)))
	    	{
	    	    if (c_orient[i][6] == POSITIVE_ORIENTATION)
	    	    	dncur[i][6] = *c;
	            else
			dncur[i][4] = *c;
		}
		else
		{
		    if (c_orient[i][0] == POSITIVE_ORIENTATION)
			dncur[i][0] = *c;
		    else
			dncur[i][5] = *c;
		}
	    }


	    /* Find states around the diffraction node */

	    assign_ahead_states_at_shock_diffraction(
			(c_orient[i][0]==POSITIVE_ORIENTATION) ? 0.0 : 1.0,
			NULL,dncur[i][0],c_orient[i][0],
			(c_orient[i][6]==POSITIVE_ORIENTATION) ? 0.0 : 1.0,
			NULL,dncur[i][6],c_orient[i][6],
			(i==0) ? comp : NULL,RP[i]);


	    node_vel_set[i] = estimate_node_vel(dncur[i][6],c_orient[i][6],
					        t[i][1],dncur[i][0],
						c_orient[i][0],t[i][0],
						i_to_f_dir[i],
						node_type(diffn[i]),
						nod_v[i],front);
	}

	if (numcr == 2 && node_vel_set[0] && node_vel_set[1])
	{
	    adjust_node_velocities_for_symmetry(nod_v,dim);
	    if (!untangle_satisfies_CFL(diffn,nod_v,front,flag))
	    {
	        (void) printf("WARNING in scalar_vector_unravel(), "
	                      "possible CFL violation\n");
	        interpolate_intfc_states(front->interf) = sav_intrp;
	        DEBUG_LEAVE(scalar_vector_unravel)
	        return MODIFY_TIME_STEP_TO_UNTANGLE_BUT_FORCE_TANGLE;
	    }
	}

	for (i = 0; i < numcr; ++i)
	{
	    if (node_vel_set[i])
	    {
		if (DEBUG)
		    (void) printf("Calling is_regular_diffraction_node()\n");

		switch(is_regular_diffraction_node(Coords(diffn[i]->posn),
						   nod_v[i],(double *)NULL,
	                                           t[i],RP[i],(RP_DATA *)NULL,
						   &is_reflected_shock[i],
						   front,DIFFRACTION_NODE,
						   ndflag))
		{
		case REGULAR_DIFFRACTION:
		    Node_vel(diffn[i])[0] = nod_v[i][0];
		    Node_vel(diffn[i])[1] = nod_v[i][1];
		    break;
	        case ANOMALOUS_REFLECTION:
	        case REGULAR_TO_MACH_DIFFRACTION:
	        case ERROR_DIFFRACTION:
	        default:
		    (void) printf("WARNING in scalar_vector_unravel(), "
		                  "is_regular_diffraction_node() failed,"
		                  " POSSIBLE BIFURCATION, CODE NEEDED\n");
	            interpolate_intfc_states(front->interf) = sav_intrp;
		    DEBUG_LEAVE(scalar_vector_unravel)
		    return ERROR_IN_UNTANGLE;
	        }
	    }
	    else
	    {
		NODE_FLAG flag;
		clear_node_flag(flag);
		(void) printf("WARNING in scalar_vector_unravel(), "
		              "Can't estimate node velocity\n");
	        if (!is_small_inc_ang_reg_diff_node(t[i],RP[i],
	                                               &is_reflected_shock[i],
						       flag))
		{
		    (void) printf("WARNING in scalar_vector_unravel(), "
		                  "is_small_inc_ang_reg_diff_node() failed,"
		                  " POSSIBLE BIFURCATION, CODE NEEDED\n");
		    interpolate_intfc_states(front->interf) = sav_intrp;
		    DEBUG_LEAVE(scalar_vector_unravel)
		    return ERROR_IN_UNTANGLE;
	        }
	    }

	    if (i == 0)
	    {
		/* NOTE: there is a potential problem here since the
		 * tracking is determined by the arbitrary "first"
		 * node.  It would be better to use the weaker of
		 * the two nodes. */

		set_track_and_newcomp_list_for_shock_diffraction(track,comp,
							w_type[i],RP[i],
							is_reflected_shock[i],
							front);
	    }
	    set_states_and_comps_about_shock_diffraction(comp,i_to_f_dir[i],
	                                                 c_orient[i],l_comp[i],
	                                                 r_comp[i],l_st[i],
	                                                 r_st[i],RP[i]);
	}

	if ((cross1!=NULL) && (is_reflected_shock[0]!=is_reflected_shock[1]))
	{
	    (void) printf("WARNING in scalar_vector_unravel(), "
	                  "inconsistent companion "
	                  "diffraction configurations\n");
	    interpolate_intfc_states(front->interf) = sav_intrp;
	    DEBUG_LEAVE(scalar_vector_unravel)
	    return ERROR_IN_UNTANGLE;
	}

	if (cross1)
	{
	    if (!track[5])
	    {
	    	(void) delete_curve(dncur[0][5]);
	    	dncur[0][5] = dncur[1][5] = NULL;
	    }
	    for (j = 1; j < 6; ++j)
	    	delete_interior_points_of_curve(front,dncur[0][j]);
	}
	else
	{
	    track[1] = track[2] = track[3] = NO;
	    n = Node_of(dncur[0][5],Opposite_orient(c_orient[0][5]));
	    if (node_type(n) != WAVE_END_NODE)
	    {
	    	(void) printf("WARNING in scalar_vector_unravel(), Unexpected "
			      "case, single cross without wave end node\n");
	    	interpolate_intfc_states(front->interf) = sav_intrp;
	    	DEBUG_LEAVE(scalar_vector_unravel)
	    	return ERROR_IN_UNTANGLE;
	    }

	    if (track[5])
	    	delete_interior_points_of_curve(front,dncur[0][5]);
	    else
	    {
	    	(void) delete_curve(dncur[0][5]);
	    	(void) delete_node(n);
	    	dncur[0][5] = dncur[1][5] = NULL;
	    }
	}
	if (DEBUG)
	{
	    for (i = 0; i < 2; ++i)
	    {
	    	if (i == 1 && !cross1) continue;
	    	for (j = 0; j < 7; ++j) 
	    	{
	    	    (void) printf("curve dncur[%d][%d]\n",i,j);
		    print_curve(dncur[i][j]);
		}
	    }
	}
	if (DEBUG)
	    (void) printf("After delete interior points of joint curves\n");

	/* 
	*	Assign components and start-end states of reflected, transmitted
	*	and deflected waves.
	*/


	if (DEBUG)
	{
	    print_angle_direction("i_to_f_dir[0] = ",i_to_f_dir[0],"\n");
	    print_angle_direction("i_to_f_dir[1] = ",i_to_f_dir[1],"\n");
	}

	ns = (c_orient[0][1] == POSITIVE_ORIENTATION) ? diffn[0] : diffn[1];
	ne = (c_orient[0][1] == POSITIVE_ORIENTATION) ? diffn[1] : diffn[0];
	for (j = 0; j < 7; ++j)
	{
	    if (!track[j])
	        continue;
	    if (dncur[0][j] == NULL)
	    {
	    	dncur[0][j] = dncur[1][j] = make_curve(l_comp[0][j],
							r_comp[0][j],ns,ne);
	    }
	    else
	    {
	    	/* Reset components */
	    	negative_component(dncur[0][j]) = l_comp[0][j];
	    	positive_component(dncur[0][j]) = r_comp[0][j];
	    }
	    wave_type(dncur[0][j]) = w_type[0][j]; 

	    /* Assign states at nodes */

	    for (i = 0; i < 2; ++i)
	    {
	    	if (i == 0 || cross1)
	    	{
	    	    set_status_at_node(dncur[i][j],c_orient[i][j],status[j]);
		    ls = Left_state_at_node(dncur[i][j],c_orient[i][j]);
		    rs = Right_state_at_node(dncur[i][j],c_orient[i][j]);
		    ft_assign(ls,l_st[i][j],front->sizest);
		    ft_assign(rs,r_st[i][j],front->sizest);
		}
		else if (j == 5)
		{
		    ls = Left_state_at_node(dncur[0][5],
					    Opposite_orient(c_orient[0][5]));
		    rs = Right_state_at_node(dncur[0][5],
					     Opposite_orient(c_orient[0][5]));
		    ft_assign(ls,l_st[0][5],front->sizest);
		    ft_assign(rs,r_st[0][5],front->sizest);
		}
	    }
	}
	for (j = 1; j < 6; ++j)
	{
	    if (dncur[0][j] != NULL)
	    {
	        /* Ensure redistribution won't interfere with untangle */
	        if (DEBUG)
	        {
	            (void) printf("Disabling redistribute of dncur[0][%d] "
	                          "%llu\n",j,curve_number(dncur[0][j]));
	        }
	        do_not_redistribute(dncur[0][j]) = YES;
	    }
	}

	    /* Insert points in intermediate curves to satisfy angles */

	if (!cross1)
	{
	    double *h = front->rect_grid->h;
	    double dir[MAXD], alpha;

	    if (track[5])
	    {
	    	n = Node_of(dncur[0][5],Opposite_orient(c_orient[0][5]));
		len = bond_length(Bond_at_node(dncur[0][5],c_orient[0][5]));
		Coords(n->posn)[0] =
			    Coords(diffn[0]->posn)[0] + len*cos(RP[0]->ang[5]);
		Coords(n->posn)[1] =
			    Coords(diffn[0]->posn)[1] + len*sin(RP[0]->ang[5]);
	    }

	    len = Front_spacing(front,GENERAL_WAVE);
	    bb = Bond_at_node(dncur[0][4],c_orient[0][4]);
	    while (scaled_bond_length(bb,h,dim) < len)
	    {
	    	if (!Following_bond(bb,c_orient[0][4])) break;
	    	(void) delete_point_adjacent_to_node(front,dncur[0][4],
				                     c_orient[0][4]);
		bb = Bond_at_node(dncur[0][4],c_orient[0][4]);
	    }
	    if (scaled_bond_length(bb,h,dim) > len)
	    {
	    	dir[0] = cos(RP[0]->ang[4]);
	    	dir[1] = sin(RP[0]->ang[4]);
	    	alpha = len/scaled_hypot(dir,h,dim);
	    	coords[0] = Coords(diffn[0]->posn)[0] + alpha*dir[0];
	    	coords[1] = Coords(diffn[0]->posn)[1] + alpha*dir[1];
	    	insert_point_adjacent_to_node(Point(coords),
					      dncur[0][4],c_orient[0][4]);
	    }
	    interpolate_intfc_states(front->interf) = sav_intrp;
	    DEBUG_LEAVE(scalar_vector_unravel)
	    return CURVES_UNTANGLED;
	}

	for (j = 1; j < 6; ++j)
	{
	    if (!track[j]) continue;
	    tgnt[0][0] = cos(RP[0]->ang[j]);tgnt[0][1] = sin(RP[0]->ang[j]);
	    tgnt[1][0] = cos(RP[1]->ang[j]);tgnt[1][1] = sin(RP[1]->ang[j]);
	    if (!intersect_ray_with_sector(diffn[0]->posn,diffn[1]->posn,
					      tgnt[0],tgnt+1,coords,dim))
	    {
		if (j == 4)		/* contacts don't meet */
		{
		    insert_3bond_wave(diffn,dncur,c_orient,tgnt,j,dim);
		    continue;
		}
		(void) printf("WARNING in scalar_vector_unravel(), "
		              "Inconsistent orientations at opposite nodes\n");
		print_curve(dncur[0][j]);
		(void) printf("Can't join rays, %d\n",j);
		(void) printf("diffn[0]\n"); print_node(diffn[0]);
		(void) printf("diffn[1]\n"); print_node(diffn[1]);
		(void) printf("RP_DATA[0]\n");
		print_RP_DATA(RP[0],Node_vel(diffn[0]));
		(void) printf("RP_DATA[1]\n");
		print_RP_DATA(RP[1],Node_vel(diffn[1]));
		interpolate_intfc_states(front->interf) = sav_intrp;
		DEBUG_LEAVE(scalar_vector_unravel)
		return ERROR_IN_UNTANGLE;
	    }
	    insert_point_adjacent_to_node(Point(coords),dncur[0][j],
					  c_orient[0][j]);
	    if (DEBUG) print_curve(dncur[0][j]);
	}
	if (!consistent_scalar_vector_unravel(dncur[0]))
	{
	    (void) printf("WARNING in scalar_vector_unravel(), "
		          "Inconsistent orientations at opposite nodes\n");
	    interpolate_intfc_states(front->interf) = sav_intrp;
	    DEBUG_LEAVE(scalar_vector_unravel)
	    return ERROR_IN_UNTANGLE;
	}
	if (DEBUG && debugging("untan_front"))
	{
	    (void) printf("Interface after scalar_vector_unravel()\n");
	    print_interface(front->interf);
	    if (debugging("untan_states"))
	    {
		(void) printf("Front states after scalar_vector_unravel()\n");
		show_intfc_states(front->interf);
	    }
	}
	interpolate_intfc_states(front->interf) = sav_intrp;
	DEBUG_LEAVE(scalar_vector_unravel)
	return CURVES_UNTANGLED;
}		/*end scalar_vector_unravel*/

LOCAL int unitary_scalar_vector_unravel(
	Front		*front,
	CROSS		*cross)
{
	COMPONENT	newcomp;
	CURVE		*c1, *c2;
	O_CURVE		Oc;
	UNTRACK_FLAG	tflg;

	if (!is_rarefaction_wave(wave_type(cross->c1)))
	{
	    (void) printf("WARNING in unitary_scalar_vector_unravel(), "
	                  "code needed for shock interaction\n");
	    return ERROR_IN_UNTANGLE;
	}
	c1 = cross->c1;	c2 = cross->c2;
	newcomp = ERROR;
	if (negative_component(c1) == negative_component(c2))
		newcomp = negative_component(c1);
	else if (negative_component(c1) == positive_component(c2))
		newcomp = negative_component(c1);
	else if (positive_component(c1) == negative_component(c2))
		newcomp = positive_component(c1);
	else if (positive_component(c1) == positive_component(c2))
		newcomp = positive_component(c1);
	if (newcomp == ERROR)
	{
		(void) printf("WARNING in unitary_scalar_vector_unravel(), ");
		(void) printf("inconsistent components\n");
		return ERROR_IN_UNTANGLE;
	}
	Oc.curve = c1;	Oc.orient = POSITIVE_ORIENTATION;
	set_untrack_flag(tflg,Oc.orient,YES,YES,YES,YES,NO);
	(void) untrack_curve(&Oc,NULL,newcomp,front->dt,front,NULL,NULL,tflg);
	return CURVES_UNTANGLED;
}		/*end unitary_scalar_vector_unravel*/

/*
*			insert_3bond_wave():
*
*	This function is used for installing opposite nodes when 
*	corresponding waves point in different directions.  We construct
*	a perp bisector of the segment between the two nodes, and insert
*	two new points.  The tangents are the rays out of the nodes with
*	the correct angles.  We locate the two new points on these rays,
*	halfway between the nodes and the intersections of the bisector
*	and the tangents.
*/

LOCAL void insert_3bond_wave(
	NODE		**diffn,
	CURVE		***dncur,
	ORIENTATION	**orient,
	double		**tgnt,
	int		index,
	int		dim)
{
	double		*p0, *p1;
	double		disp[MAXD], sep, r;
	double		coords[MAXD];
	int		j;

	p0 = Coords(diffn[0]->posn);
	p1 = Coords(diffn[1]->posn);

	for (j = 0; j < dim; ++j)
	    disp[j] = p1[j] - p0[j];
	sep = mag_vector(disp,dim);

	r = 0.25 * sqr(sep) / scalar_product(disp,tgnt[0],dim);
	for (j = 0; j < dim; ++j)
	    coords[j] = p0[j] + r * tgnt[0][j];
	insert_point_adjacent_to_node(Point(coords),
					dncur[0][index],orient[0][index]);

	r = - 0.25 * sqr(sep) / scalar_product(disp,tgnt[1],dim);
	for (j = 0; j < dim; ++j)
	    coords[j] = p1[j] + r * tgnt[1][j];
	insert_point_adjacent_to_node(Point(coords),
					dncur[1][index],orient[1][index]);
}		/*end insert_3bond_wave*/

LOCAL	int consistent_scalar_vector_unravel(
	CURVE		**dncur)
{
	int		i, j;
	BOND		*b1, *b2;
	POINT		P;

	for (i = 1; i < 6; ++i)
	{
	    if (dncur[i] == NULL) continue;
	    for (j = i+1; j < 6; ++j)
	    {
	    	if (dncur[j] == NULL) continue;
	    	for (b1 = dncur[i]->first; b1 != NULL; b1 = b1->next)
	    	{
	    	    for (b2 = dncur[j]->first; b2 != NULL; b2 = b2->next)
	    	    {
	    	    	if (b1->start == b2->start) continue;
	    	    	if (b1->end == b2->start) continue;
	    	    	if (b1->start == b2->end) continue;
	    	    	if (b1->end == b2->end) continue;
	    	    	if (cross_bonds(b1,b2,&P))
	    	    	{
	    	    	    return NO;
	    	    	}
	    	    }
	    	}
	    }
	}
	return YES;
}		/*end consistent_scalar_vector_unravel*/


LOCAL	void adjust_node_velocities_for_symmetry(
	double		**nod_v,
	int		dim)
{
	double		mag[2], alpha[2], tmpx[2], tmpy[2];
	double		nx, ny, tx, ty, len;
	static const double 	SYMM_TOL = 0.5;	/* TOLERANCE */
	int		i;

	DEBUG_ENTER(adjust_node_velocites_for_symmetry)
	if (DEBUG)
	{
	    (void) printf("Input node velocities\n");
	    for (i = 0; i < 2; ++i)
	    {
	    	(void) printf("\tNode velocity %d = <%g %g>\n",
			      i,nod_v[i][0],nod_v[i][1]);
	    }
	}

	/* Check if node velocities are nearly symmetric about their bisector */
	for (i = 0; i < 2; ++i)
	    mag[i] = mag_vector(nod_v[i],dim);
	if (fabs(mag[0] - mag[1]) < SYMM_TOL * max(mag[0],mag[1]))
	{
	    tx = nod_v[0][0] - nod_v[1][0];	ty = nod_v[0][1] - nod_v[1][1];
	    len = hypot(tx,ty);
	    tx /= len;			ty /= len;
	    nx = ty;			ny = -tx;

	    for (i = 0; i <2; ++i)
	    	alpha[i] = nod_v[i][0]*nx + nod_v[i][1]*ny;

	    tmpx[0] = 0.5*(nod_v[0][0] - nod_v[1][0]) + alpha[1]*nx;
	    tmpy[0] = 0.5*(nod_v[0][1] - nod_v[1][1]) + alpha[1]*ny;
	    tmpx[1] = 0.5*(nod_v[1][0] - nod_v[0][0]) + alpha[0]*nx;
	    tmpy[1] = 0.5*(nod_v[1][1] - nod_v[0][1]) + alpha[0]*ny;
	    nod_v[0][0] = tmpx[0]; nod_v[0][1] = tmpy[0];
	    nod_v[1][0] = tmpx[1]; nod_v[1][1] = tmpy[1];
	    if (DEBUG)
	    {
	    	(void) printf("After reset node velocities\n");
	    	for (i = 0; i < 2; ++i)
	    	{
	    	    (void) printf("\tNode velocity %d = <%g %g>\n",
		    		  i,nod_v[i][0],nod_v[i][1]);
		}
	    }
	}
	DEBUG_LEAVE(adjust_node_velocites_for_symmetry)
}		/*end adjust_node_velocities_for_symmetry*/

LOCAL	int untangle_satisfies_CFL(
	NODE		**node,
	double		**nod_v,
	Front		*front,
	int		flag)
{
	double		sep;
	double		*h = front->rect_grid->h;
	double		dir[MAXD], d;
	double		*posn0, *posn1;
	double		rel_speed;
	double		*dt_frac;
	double		dt;
	double		frac;
	double		max_sep = Max_new_node_separation(front);
	int		i, dim = front->rect_grid->dim;

	if (Apply_CFL_at_nodes(front) == NO)
	    return YES;

	if ((flag == LAST_ATTEMPT_TO_UNTANGLE) ||
				(last_time_step_modification(front) == YES))
	    return YES;

	if (node[1] == NULL) /* TODO */
	    return YES;

	posn0 = Coords(node[0]->posn);
	posn1 = Coords(node[1]->posn);

	sep = _scaled_separation(posn0,posn1,h,dim);
	if (sep < max_sep)
	    return YES;

	for (i = 0; i < dim; ++i)
	    dir[i] = posn1[i] - posn0[i];
	d = mag_vector(dir,dim);
	for (i = 0; i < dim; ++i)
	    dir[i] /= d;
	rel_speed = scalar_product(nod_v[1],dir,dim) -
		    scalar_product(nod_v[0],dir,dim);

	dt_frac = front->dt_frac;
	dt = front->dt;

	frac = 1.0 - (d/(rel_speed*dt))*(1.0 - max_sep/sep);
	frac = max(0.0,frac);

	/* To be conservative expand frac a little */
	if (frac < 1.0/Time_step_increase_factor(front))
	    frac *= Time_step_increase_factor(front);

	frac = max(Min_time_step_modification_factor(front),frac);
	frac = min(Max_time_step_modification_factor(front),frac);
	*dt_frac = min(frac,*dt_frac);
	(void) printf("WARNING in untangle_satisfies_CFL(), "
	              "Bifurcation violates CFL\n"
	              "scaled separation = %g, max = %g, frac = %g\n",
		      sep,max_sep,frac);
	print_general_vector("dir = ",dir,dim,"\n");
	(void) printf("rel_speed = %g\n",rel_speed);
	(void) printf("node[0]\n"); print_node(node[0]);
	(void) printf("node[1]\n"); print_node(node[1]);
	return NO;
}		/*end untangle_satisfies_CFL*/
#endif /* defined(FULL_PHYSICS) */
#endif /* defined(TWOD) */
