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
*				g2dprint.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains printing routines for gas dynamics that are specific
*	to dynamic simulations.
*	TODO: This file should be dim independent and renamed!
*
*	g_fgraph_front_states() is set to a function pointer
*	in the Front structure.
*
*/

#if defined(TWOD)

#include <gdecs/gdecs.h>
#if defined(sparc)
#include <sys/param.h>
#endif /* defined(sparc) */

struct _CLIST {
	CURVE *c;
	double pbar[MAXD];
	struct _CLIST *next, *prev;
};
typedef struct _CLIST CLIST;

	/* LOCAL Function Declarations */
LOCAL	boolean	curve_is_in_clist(CURVE*,CLIST*);
LOCAL	boolean	select_contacts(CURVE*);
LOCAL	boolean	select_contacts_and_dirichlet_bdrys(CURVE*);
LOCAL	boolean	select_le(CURVE*);
LOCAL	boolean	select_te(CURVE*);
LOCAL	void	expand_clist(CLIST*,CLIST**,CURVE***,
			     boolean (*)(CURVE*),ORIENTATION);
LOCAL	void	find_average_curve_position(double*,CURVE*);
LOCAL	void	print_Mach_node_stats(NODE*,Front*,Grid*,OUTPUT_DATA*);
LOCAL	void	print_reg_refl_node_stats(NODE*,Front*,Grid*,OUTPUT_DATA*);
LOCAL	void	show_selected_curve_states(FILE*,Front*,boolean (*)(CURVE*),
					   const char*,boolean);

/*
*			g_fgraph_front_states():
*
*	Prints the states on front in a columnar format that
*	can be used by "graphs".
*/

EXPORT void g_fgraph_front_states(
	FILE		*file,
	Front		*front)
{
	CURVE		**c;
	char		title[80];
	double		length;
	int		num;

	if ( (front->_fprint_header_for_graph_curve_states == NULL) ||
				(front->_fgraph_curve_states == NULL))
		return;

	(void) foutput(file);
	(void) fprintf(file,"\t\tSTATES ON THE FRONT\n\n");
	num = 0;
	for (c = front->interf->curves; c && *c; ++c)
	{
	    switch (wave_type(*c))
	    {
	    case SUBDOMAIN_BOUNDARY:
	    	break;
		
	    default:
	    	(void) sprintf(title,"along CURVE_%llu (%s)",
			       curve_number(*c),
			       wave_type_as_string(wave_type(*c),front->interf));
	    	fprint_header_for_graph_curve_states(file,front,title);
	    	length = 0.;
	    	fgraph_curve_states(file,*c,front,&length);
	    	(void) fprintf(file,"\n\nEnd %s\n",title);
	    	++num;
	    }
	}
	(void) fprintf(file,"\n\n");
	(void) fprintf(file,"\t\tEND OF STATES ON THE FRONT\n\n");
}		/*end g_fgraph_front_states*/

LOCAL	boolean select_contacts(CURVE *c)
{
	return (is_scalar_wave(wave_type(c)) || 
		is_thinflame_wave(wave_type(c))) ? YES : NO;
}		/*end select_contacts*/

LOCAL   boolean select_contacts_and_dirichlet_bdrys(CURVE *c)
{
	int     w_type = wave_type(c);
	return (is_scalar_wave(w_type) || (w_type == DIRICHLET_BOUNDARY)) ?
	       YES : NO;
}		/*end select_contacts_and_dirichlet_bdrys*/

LOCAL	boolean select_le(CURVE *c)
{
	return (is_rarefaction_leading_edge(wave_type(c))) ? YES : NO;
}		/*end select_le*/

LOCAL	boolean select_te(CURVE *c)
{
	return (is_rarefaction_trailing_edge(wave_type(c))) ? YES : NO;
}		/*end select_te*/

LOCAL	void show_selected_curve_states(
	FILE		*file,
	Front		*front,
	boolean		(*select)(CURVE*),
	const char	*mesg,
	boolean         index)
{
	CLIST		*head, *cl, *clnext;
	CURVE		*c, **carray;
	double		length;

	if ((front->_fprint_header_for_graph_curve_states == NULL) ||
	    (front->_fgraph_curve_states == NULL))
	    return;

	carray = NULL;
	(void) next_curve(front->interf,NULL);
	while (next_curve(front->interf,&c))
	{
	    if ((*select)(c))
	    {
	    	if (!add_to_pointers(c,&carray))
	    	{
	    	    screen("ERROR in show_selected_curve_states(), "
	    	           "add_to_pointers() failed\n");
	    	    clean_up(ERROR);
	    	}
	    }
	}

	(void) foutput(file);
	(void) fprintf(file,"\t\tSTATES ON THE FRONT\n\n");
	if (index == YES)
	{
	    char *title;
	    int  icurve;
	    size_t  ncurves = size_of_pointers(carray);
	    size_t  ndigits;

	    for (ndigits = 0; ncurves != 0; ++ndigits, ncurves /= 10);
	    uni_array(&title,strlen(mesg)+ndigits+2,CHAR);
            icurve = -1;
            while (carray != NULL)
            {
                ++icurve;
                c = carray[0];
                if (!delete_from_pointers(c,&carray))
                {
                    screen("ERROR in show_selected_curve_states(), "
                           "delete_from_pointers() failed\n");
                    clean_up(ERROR);
                }
                (void) sprintf(title,"%s_%d",mesg,icurve);
                fprint_header_for_graph_curve_states(file,front,title);
                length = 0.0;
                fgraph_curve_states(file,c,front,&length);
                (void) fprintf(file,"\n\nEnd %s\n\n",title);
            }
	    free(title);
	}
	else
	{
	    while (carray != NULL)
	    {
	        c = carray[0];
	        if (!delete_from_pointers(c,&carray))
	        {
	    	    screen("ERROR in show_selected_curve_states(), "
	    	           "delete_from_pointers() failed\n");
	    	    clean_up(ERROR);
	        }
	        scalar(&head,sizeof(CLIST));
	        head->c = c;
	        head->prev = head->next = NULL;
	        find_average_curve_position(head->pbar,c);
	        expand_clist(head,&head,&carray,select,POSITIVE_ORIENTATION);
	        expand_clist(head,&head,&carray,select,NEGATIVE_ORIENTATION);
	        (void) fprintf(file,"\n\n");
	        fprint_header_for_graph_curve_states(file,front,mesg);
	        length = 0.0;
	        for (cl = head; cl != NULL; cl = cl->next)
	    	    fgraph_curve_states(file,cl->c,front,&length);
	        (void) fprintf(file,"\n\nEnd %s\n\n",mesg);
	        for (cl = head, clnext = NULL; cl != NULL; cl = clnext)
	        {
	    	    clnext = cl->next;
	    	    free(cl);
	        }
	    }
	}
	(void) fprintf(file,"\n\n\t\tEND OF STATES ON THE FRONT\n\n");
}		/*end show_selected_curve_states*/

LOCAL	boolean curve_is_in_clist(
	CURVE		*c,
	CLIST		*head)
{
	CLIST		*cl;

	for (cl = head; cl != NULL; cl = cl->next)
		if (cl->c == c) return YES;
	return NO;
}		/*end curve_is_in_clist*/


LOCAL	void find_average_curve_position(
	double		*p,
	CURVE		*c)
{
	BOND		*b;
	int		dim = c->interface->dim;
	int		i, num;

	for (i = 0; i < dim; ++i)
	    p[i] = Coords(c->first->start)[i];
	for (num = 1, b = c->first; b != NULL; b = b->next, ++num)
	{
	    for (i = 0; i < dim; ++i)
		p[i] += Coords(b->end)[i];
	}
	for (i = 0; i < dim; ++i)
	    p[i] /= num;
}		/*end find_average_curve_position*/

LOCAL	void expand_clist(
	CLIST		*cl,
	CLIST		**head,
	CURVE		***carray,
	boolean		(*select)(CURVE*),
	ORIENTATION	orient)
{
	CLIST		*newcl;
	CURVE		**c;
	NODE		*n = Node_of(cl->c,Opposite_orient(orient));

	for (c = (orient==POSITIVE_ORIENTATION) ? n->out_curves : n->in_curves;
		c && *c; ++c)
	{
	    if ((*select)(*c) == NO)
	        continue;
	    if (curve_is_in_clist(*c,*head))
		continue;
	    if (!delete_from_pointers(*c,carray))
	    {
	    	screen("ERROR in expand_clist(), "
	    	       "delete_from_pointers() failed\n");
	    	clean_up(ERROR);
	    }
	    scalar(&newcl,sizeof(CLIST));
	    newcl->c = *c;
	    find_average_curve_position(newcl->pbar,*c);
	    if (orient == POSITIVE_ORIENTATION)
	    {
	    	newcl->next = NULL;
	    	newcl->prev = cl;
	    	cl->next = newcl;
	    }
	    else
	    {
	    	newcl->prev = NULL;
	    	newcl->next = cl;
	    	cl->prev = newcl;
	    	*head = newcl;
	    }
	    expand_clist(newcl,head,carray,select,orient);
	}
}		/*end expand_clist*/

/*
*		show_front_states_for_expanding_shocks():
*
*	Prints the states on selected curves of the front in a
*	columnar format that can be used by "graphs".
*
*	For use in EXPANDING_SHOCK problems.
*
*	TODO: SHOULD PRINT MAXES AND MINS FOR SCALING THE GRAPHS
*/

/*ARGSUSED*/
EXPORT void show_front_states_for_expanding_shocks(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	CURVE		*shock,*first_wall,*second_wall;
	FILE		*logfile = Output_file(data);
	char		title[80];
	double		length;

	if ((front->_fprint_header_for_graph_curve_states == NULL) ||
	    (front->_fgraph_curve_states == NULL))
	    return;

	{
	    int	num;
	    CURVE	**c;

	    for (num = 0,c = front->interf->curves; *c != NULL; ++num,++c)
	    {
	    	switch (num) 
	    	{
	    	case 0:
	    	    second_wall = *c;
	    	    continue;
	    	case 1:		
	    	    first_wall = *c;
	    	    continue;
	    	case 2:
	    	    shock = *c;
	    	    continue;
	    	}
	    }
	}

	(void) foutput(logfile);
	(void) fprintf(logfile,"\t\tSTATES ON THE FRONT\n\n");

	(void) sprintf(title,"along SHOCK (%llu)",curve_number(shock));
	fprint_header_for_graph_curve_states(logfile,front,title);
	length = 0.;
	fgraph_curve_states(logfile,shock,front,&length);
	(void) fprintf(logfile,"\n\nEnd %s\n",title);

	(void) sprintf(title,"along WALLS (%llu) and (%llu)",
		curve_number(first_wall),curve_number(second_wall));
	fprint_header_for_graph_curve_states(logfile,front,title);
	length = 0.;
	fgraph_curve_states(logfile,first_wall,front,&length);
	fgraph_curve_states(logfile,second_wall,front,&length);
	(void) fprintf(logfile,"\n\nEnd %s\n",title);

	(void) fprintf(logfile,"\n\n");
	(void) fprintf(logfile,"\t\tEND OF STATES ON THE FRONT\n\n");
}		/*end show_front_states_for_expanding_shocks*/


/*
*		show_front_states_for_bowshocks():
*
*	Prints the states on selected curves of the front in a
*	columnar format that can be used by "graphs".
*
*	For use in BOWSHOCK problems.
*
*	TODO: SHOULD PRINT MAXES AND MINS FOR SCALING THE GRAPHS
*/

/*ARGSUSED*/
EXPORT void show_front_states_for_bowshocks(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	CURVE		*bow_shock,*wedge,*xit;
	FILE		*logfile = Output_file(data);
	char		title[80];
	double		length;

	if ( (front->_fprint_header_for_graph_curve_states == NULL) ||
				(front->_fgraph_curve_states == NULL))
		return;
	/* TODO: METHOD FOR LOCATING INTERESTING CURVES TOO AD HOC */
	{
		int num;
		CURVE **c;

		for (num = 0,c = front->interf->curves; *c != NULL; ++num,++c)
		{
			switch (num)
			{
			case 0:
				bow_shock = *c;
				continue;
			case 2:
				wedge = *c;
				continue;
			case 6:
				xit = *c;
				continue;
			}
		}
	}

	(void) fprintf(logfile,"\n\n");
	(void) foutput(logfile);
	(void) fprintf(logfile,"\t\tSTATES ON THE FRONT\n\n");

	(void) sprintf(title,"along BOW SHOCK (%llu)",curve_number(bow_shock));
	fprint_header_for_graph_curve_states(logfile,front,title);
	length = 0.;
	fgraph_curve_states(logfile,bow_shock,front,&length);
	(void) fprintf(logfile,"\n\n");
	(void) fprintf(logfile,"End %s\n",title);

	(void) sprintf(title,"along WEDGE (%llu)",curve_number(wedge));
	fprint_header_for_graph_curve_states(logfile,front,title);
	length = 0.;
	fgraph_curve_states(logfile,wedge,front,&length);
	(void) fprintf(logfile,"\n\n");
	(void) fprintf(logfile,"End %s\n",title);

	(void) sprintf(title,"along EXIT (%llu)",curve_number(xit));
	fprint_header_for_graph_curve_states(logfile,front,title);
	length = 0.;
	fgraph_curve_states(logfile,xit,front,&length);
	(void) fprintf(logfile,"\n\n");
	(void) fprintf(logfile,"End %s\n",title);

	(void) fprintf(logfile,"\n\n");
	(void) fprintf(logfile,"\t\tEND OF STATES ON THE FRONT\n\n");
}		/*end show_front_for_bow_shocks*/

LOCAL double time_elapsed,initial_time_elapsed;

/*
*		show_front_states_for_ramp_reflections():
*
*	Prints the states on selected curves of the front in a
*	columnar format that can be used by "graphs".
*
*	For use in the RAMP_REFLECTION problem.  We assume the physical
*	curves are all oriented positively wrt the reflection point.
*	We also assume the refl wall is oriented from the top, through
*	the base of the mach stem, and ending at the corner.
*
*	TODO: SHOULD PRINT MAXES AND MINS FOR SCALING THE GRAPHS
*/

/*ARGSUSED*/
EXPORT void show_front_states_for_ramp_reflections(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	COMPONENT	comp;			/* bow for reg, else mach */
	CURVE		*bow_shock;
	CURVE		*bow_wall = NULL;
	CURVE		*refl_wall = NULL;	/* only for regular refl */
	CURVE		*mach_stem = NULL;	/* next 4 only for mach refl */
	CURVE		*contact = NULL;
	CURVE		*mach_wall = NULL;
	CURVE		*contact_wall = NULL;
	CURVE		**c;
	FILE		*logfile = Output_file(data);
	NODE		**n;
	char		title[80];
	double		length;
	double		current_time = grid->time;
	int		i, dim = front->interf->dim;
	ORIENTATION	orient;
	ORIENTATION	cw_orient;
	static double	ss_origin[MAXD];
	
	if (front->pp_grid->nn > 1)
		return;         /*this function won't work across processors*/

	ramp_reflection_corner_posn(ss_origin,NO,dim);

	for (n =  front->interf->nodes; *n != NULL; ++n)
		if ((node_type(*n) == MACH_NODE) || 
				(node_type(*n) == B_REFLECT_NODE))
			break;

	if (*n == NULL)	return;		/* shock has not passed corner
					   or has been untracked */


	if (node_type(*n) == MACH_NODE)
	{
		find_curve_with_status(*n,&contact,&orient,SLIP);
		find_curve_with_status(*n,&mach_stem,&orient,MACH_STEM);
			/* orient now mach stem wrt reflection point */
		if (orient == POSITIVE_ORIENTATION)
		{
			comp = positive_component(mach_stem);
			c = mach_stem->end->out_curves;
		}
		else
		{
			comp = negative_component(mach_stem);
			c = mach_stem->start->out_curves;
		}
		for (; *c != NULL; ++c)
		{
			if (positive_component(*c) == comp)
			{
				mach_wall = *c;
				orient = POSITIVE_ORIENTATION;
				break;
			}
		}
		if (mach_wall == NULL)
		{
			for (c = (orient == POSITIVE_ORIENTATION) ?
				mach_stem->end->in_curves :
				mach_stem->start->in_curves; *c != NULL; ++c)
			{
				if (negative_component(*c) == comp)
				{
					mach_wall = *c;
					orient = NEGATIVE_ORIENTATION;
					break;
				}
			}
		}

			/* orient now mach wall wrt base of mach stem */
		Check_return(
		    next_boundary(mach_wall,Opposite_orient(orient),
				  &contact_wall,&cw_orient),
		    show_front_states_for_ramp_reflections)
	}

	find_curve_with_status(*n,&bow_shock,&orient,REFLECTED);
			/* orient now bow shock wrt reflection point */
	if (node_type(Node_of(bow_shock,Opposite_orient(orient))) ==
							ATTACHED_B_NODE)
	{
		bow_wall = NULL;
		(void) fprintf(logfile,
			       "Bow shock is attached, bow_wall == NULL\n");
	}
	else if (node_type(*n) == MACH_NODE)
	{
		if (orient == POSITIVE_ORIENTATION)
		{
			comp = negative_component(bow_shock);
			c = bow_shock->end->out_curves;
		}
		else
		{
			comp = positive_component(bow_shock);
			c = bow_shock->start->out_curves;
		}
		for (; *c != NULL; ++c)
		{
			if (negative_component(*c) == comp)
			{
				bow_wall = *c;
				orient = POSITIVE_ORIENTATION;
				break;
			}
		}
		if (bow_wall == NULL)
		{
			for (c = (orient == POSITIVE_ORIENTATION) ?
		     	bow_shock->end->in_curves :
		     	bow_shock->start->in_curves; *c != NULL; ++c)
			{
				if (positive_component(*c) == comp)
				{
					bow_wall = *c;
					orient = NEGATIVE_ORIENTATION;
					break;
				}
			}
		}
	}
	else	/* (node_type(n) == B_REFLECT_NODE) */
	{
		comp = (orient == POSITIVE_ORIENTATION) ?
				negative_component(bow_shock) :
				positive_component(bow_shock);

		for (c = (*n)->out_curves; *c != NULL; ++c)
		{
			if ((positive_component(*c) == comp) &&
			    (wave_type(*c) == NEUMANN_BOUNDARY))
			{
				refl_wall = *c;
				orient = POSITIVE_ORIENTATION;
				break;
			}
		}
		if (refl_wall == NULL)
		{
			for (c = (*n)->in_curves; *c != NULL; ++c)
			{
				if ((negative_component(*c) == comp) &&
				    (wave_type(*c) == NEUMANN_BOUNDARY))
				{
					refl_wall = *c;
					orient = NEGATIVE_ORIENTATION;
					break;
				}
			}
		}

		/* orient is now refl wall wrt reflection point */
		Check_return(next_boundary(refl_wall,orient,&bow_wall,&orient),
			     show_front_states_for_ramp_reflections)
	}

	time_elapsed = current_time + initial_time_elapsed;
	if (debugging("self_similar"))
	{
		(void) printf("ss_origin: ");
		for (i = 0; i < dim; ++i)
			(void) printf("%g ",ss_origin[i]);
		(void) printf("\n");
		(void) printf("initial_time_elapsed: %g\n",
			      initial_time_elapsed);
		(void) printf("current_time: %g time_elapsed: %g\n",
			current_time,time_elapsed);
	}

	(void) fprintf(logfile,"\n\n");
	(void) foutput(logfile);
	(void) fprintf(logfile,"\t\tSTATES ON THE FRONT\n\n");

	(void) sprintf(title,
		       "along REFLECTED SHOCK (%llu)",curve_number(bow_shock));
	print_header_for_self_similar_states_along_curve(logfile,title,dim);
	length = 0.;
	print_self_similar_front_states_along_curve(logfile,bow_shock,front,
		ss_origin,time_elapsed,&length);
	(void) fprintf(logfile,"\n");

	length = 0.;
	if (mach_stem == NULL)
	{
	    if (bow_wall == NULL)
		    (void) sprintf(title,"along REFLECTION WALLS (%llu) ",
			    curve_number(refl_wall));
	    else
		(void) sprintf(title,"along REFLECTION WALLS (%llu) and (%llu)",
			curve_number(refl_wall),curve_number(bow_wall));
	    print_header_for_self_similar_states_along_curve(logfile,title,dim);
	    print_self_similar_front_states_along_curve(logfile,refl_wall,
					front,ss_origin,time_elapsed,&length);
	}
	else
	{
	    if (bow_wall == NULL)
		(void) sprintf(title,"along REFLECTION WALLS (%llu) and (%llu)",
			curve_number(mach_wall),curve_number(contact_wall));
	    else
		(void) sprintf(title,
			"along REFLECTION WALLS (%llu), (%llu), and (%llu)",
			curve_number(mach_wall),curve_number(contact_wall),
			curve_number(bow_wall));
	    print_header_for_self_similar_states_along_curve(logfile,title,dim);
	    print_self_similar_front_states_along_curve(logfile,mach_wall,
				front,ss_origin,time_elapsed,&length);
	    while (node_type(Node_of(contact_wall,
				 Opposite_orient(cw_orient))) == NEUMANN_NODE)
	    {
		print_self_similar_front_states_along_curve(logfile,
			   contact_wall,front,ss_origin,time_elapsed,&length);
		Check_return(
		    next_boundary(contact_wall,Opposite_orient(cw_orient),
				  &contact_wall,&cw_orient),
		    show_front_states_for_ramp_reflections)
	    }
	    print_self_similar_front_states_along_curve(logfile,contact_wall,
			front,ss_origin,time_elapsed,&length);
	}
	if (bow_wall != NULL)
		print_self_similar_front_states_along_curve(logfile,bow_wall,
				    front,ss_origin,time_elapsed,&length);
	(void) fprintf(logfile,"\n");

	if (mach_stem != NULL)
	{
		(void) sprintf(title,
			       "along MACH STEM (%llu)",
			       curve_number(mach_stem));
		print_header_for_self_similar_states_along_curve(logfile,
								 title,dim);
		length = 0.;
		print_self_similar_front_states_along_curve(logfile,mach_stem,
			front,ss_origin,time_elapsed,&length);
		(void) fprintf(logfile,"\n");

		(void) sprintf(title,
			       "along CONTACT (%llu)",curve_number(contact));
		print_header_for_self_similar_states_along_curve(logfile,
								 title,dim);
		length = 0.;
		print_self_similar_front_states_along_curve(logfile,contact,
			front,ss_origin,time_elapsed,&length);
		(void) fprintf(logfile,"\n");
	}

	(void) fprintf(logfile,"\n\n");
	(void) fprintf(logfile,"\t\tEND OF STATES ON THE FRONT\n\n");
}		/*end show_front_states_for_ramp_reflection*/


/*
*			show_front_states_for_rm_problem():
*
*	This function is simply a driver for show_selected_curve_states()
*	for the Richtmyer-Meshkov problem.  It accepts an argument list
*	compatible with the user output functions in g_User_printplot.
*/

/*ARGSUSED*/
EXPORT	void	show_front_states_for_rm_problem(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	FILE		*logfile = Output_file(data);

	show_selected_curve_states(logfile,front,select_le,"along LE",NO);
	show_selected_curve_states(logfile,front,select_te,"along TE",NO);
	show_selected_curve_states(logfile,front,select_contacts,
				   "along CONTACT",NO);
	return;
}		/*end show_front_states_for_rm_problem*/


/*
*			show_front_states_along_lower_wall():
*
*	Picks out the boundary curves along the lower wall, and prints their
*	state information.  This is intended for use in analysing blast
*	problems (ideal and non-ideal).
*/

/*ARGSUSED*/
EXPORT void show_front_states_along_lower_wall(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	CURVE		**c, *curr_c;
	FILE		*logfile = Output_file(data);
	char		title[32];
	double		*L = front->rect_grid->L;
	double		*U = front->rect_grid->U;
	double		length = 0.0;
	ORIENTATION	curr_or;

	curr_c = NULL;
	for (c = front->interf->curves; c && *c; ++c)
	{
		if ((Coords((*c)->start->posn)[1] == L[1]) &&
		    (Coords((*c)->end->posn)[1] == L[1]))
		{
			if (Coords((*c)->start->posn)[0] == U[0])
			{
				curr_c = *c;
				curr_or = POSITIVE_ORIENTATION;
				break;
			}
			else if (Coords((*c)->end->posn)[0] == U[0])
			{
				curr_c = *c;
				curr_or = NEGATIVE_ORIENTATION;
				break;
			}
		}
	}

	if (curr_c == NULL)
	{
		screen("WARNING in show_front_states_along_lower_wall(), ");
		screen("unable to find starting curve\n");
		return;
	}

	(void) sprintf(title,"along LOWER WALLS");
	fprint_header_for_graph_curve_states(logfile,front,title);
	while (curr_c)
	{
		fgraph_curve_states(logfile,curr_c,front,&length);
		(void) next_boundary(curr_c,Opposite_orient(curr_or),
				     &curr_c,&curr_or);

		if (Coords(Node_of(curr_c,Opposite_orient(curr_or))->posn)[1]
		     					!= L[1])
			curr_c = NULL;
	}
	return;
}		/*end show_front_states_along_lower_wall*/


/*
*			show_contact_states():
*
*	This function is simply a shell for showing only contact states in
*	a selected problem.  It accepts an argument list compatible with that
*	of the user functions in g_User_printplot.
*/

/*ARGSUSED*/
EXPORT	void	show_contact_states(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	FILE		*logfile = Output_file(data);

	debug_print("show_contact","Entered show_contact_states()\n");
	show_selected_curve_states(logfile,front,select_contacts,
				   "along CONTACT",YES);
	debug_print("show_contact","Left show_contact_states()\n");
	return;
}		/*end show_contact_states*/


/*ARGSUSED*/
EXPORT	void	show_contact_and_dirichlet_bdry_states(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	FILE		*logfile = Output_file(data);

	debug_print("show_contact",
	      "Entered show_contact_and_dirichlet_bdry_states()\n");
	show_selected_curve_states(
		logfile,front,select_contacts_and_dirichlet_bdrys,
		"along CURVE",YES);
	debug_print("show_contact","Left show_contact_and_dirichlet_bdry_states()\n");
	return;
}		/*end show_contact_and_dirichlet_bdry_states*/


EXPORT void set_initial_time_elapsed(
	double value)
{
	initial_time_elapsed = value;
}		/*end set_initial_time_elapsed*/


EXPORT	void g_fprint_RP_DATA_at_nodes(
	FILE		*file,
	INTERFACE	*intfc)
{
	NODE		**n;

	if (intfc->dim != 2)
	    return;
	(void) fprintf(file,"\n\n");
	(void) foutput(file);
	(void) fprintf(file,
		       "RP_DATA information at nodes of interface %llu\n\n",
		       interface_number(intfc));
	for (n = intfc->nodes; n && *n; ++n)
	{
	    (void) fprintf(file,
	    	           "Node %llu, Rp_data(node) = %llu\n",node_number(*n),
				ptr2ull(Rp_data(*n)));
	    fprint_RP_DATA(file,Rp_data(*n),Node_vel(*n));
		
	}
	(void) fprintf(file,"\n");
	(void) foutput(file);
	(void) fprintf(file,
		"End of RP_DATA information at nodes of interface %llu\n",
		interface_number(intfc));
	(void) fprintf(file,"\n\n");
}		/*end g_fprint_RP_DATA_at_nodes*/

EXPORT	void fprint_RP_DATA(
	FILE		*file,
	RP_DATA		*RP,
	double		*v)
{
	int		i, dim;

	if (RP == NULL)
	    return;

	(void) fprintf(file,"\t\tPrintout of RP_DATA structure %llu\n\n",
		       ptr2ull(RP));
	dim = current_interface()->dim;
	(void) fprint_general_vector(file,"node velocity = ",v,dim,"\n");
	(void) fprintf(file,"ang_dir = %d",RP->ang_dir);
	fprint_angle_direction(file,"",RP->ang_dir,"\n");
	(void) fprintf(file,"Angles:");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",MAX_N_CURVES);
	    (void) fwrite((const void *)RP->ang,sizeof(double),
			  MAX_N_CURVES,file);
	    (void) fprintf(file,"\n");
	}
	else
	{
	    (void) fprintf(file,"\n");
	    for (i = 0; i < MAX_N_CURVES; ++i)
	    	(void) fprintf(file,"\tang[%d] = %"FFMT" (%g degrees)\n",
			       i,RP->ang[i],degrees(RP->ang[i]));
	}
	(void) fprintf(file,"RP State Data:\n");
	for (i = 0; i < MAX_N_CURVES; ++i)
	{
	    (void) fprintf(file,"state[%d]:\n",i);
	    (void) fprintf(file,"Mach number[%d] = ",i);
	    if (is_binary_output() == YES)
	    {
	    	(void) fprintf(file,"\f%c",1);
	    	(void) fwrite((const void *)(RP->M+i),sizeof(double),1,file);
	    }
	    else
	    	(void) fprintf(file,"%"FFMT,mach_number(RP->state[i],v));
	    (void) fprintf(file,"\n");
	    fprint_gas_data(file,RP->state[i]);
	}
	(void) fprintf(file,"Turning angles:");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",MAX_N_CURVES);
	    (void) fwrite((const void *)RP->theta,sizeof(double),
			  MAX_N_CURVES,file);
	    (void) fprintf(file,"\n");
	}
	else
	{
	    (void) fprintf(file,"\n");
	    for (i = 0; i < MAX_N_CURVES; ++i)
	    	(void) fprintf(file,"\ttheta[%d] = %"FFMT" (%g degrees)\n",
			       i,RP->theta[i],degrees(RP->theta[i]));
	}
	(void) fprintf(file,"End RP State Data:\n");
	(void) fprintf(file,"\t\tEnd printout of RP_DATA structure %llu\n\n",
		       ptr2ull(RP));

}		/*end fprint_RP_DATA*/


EXPORT	void g_fprint_ContactWallNodeParams(
	FILE	  *file,
	INTERFACE *intfc)
{
	CWNP	*cwnp;
	boolean bio;
	if (intfc->dim != 2)
	    return;

	(void) fprintf(file,"\n\n");
	(void) foutput(file);
	(void) fprintf(file,
		       "CONTACT WALL NODE PARAMS for interface %llu\n\n",
		       interface_number(intfc));

	cwnp = contact_wall_node_params(intfc);
	bio = is_binary_output();
	(void) fprintf(file,"wall_bond_len = ");
	if (bio == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *)&cwnp->wall_bond_len,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,cwnp->wall_bond_len);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"first_adjust_time = ");
	if (bio == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *)&cwnp->first_adjust_time,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,cwnp->first_adjust_time);
	(void) fprintf(file,"\n");
	(void) fprintf(file,"first_adjust_step = %d\n",cwnp->first_adjust_step);
	(void) fprintf(file,"adjust = %s\n",y_or_n(cwnp->adjust));

	(void) fprintf(file,"\n");
	(void) fprintf(file,
		"End of CONTACT WALL NODE PARAMS for interface %llu\n\n",
		interface_number(intfc));
	(void) fprintf(file,"\n\n");
}		/*end g_fprint_ContactWallNodeParams*/

#if defined(FULL_PHYSICS)
EXPORT	void print_cross_node(
	NODE		*node,
	O_CURVE		*c0,
	O_CURVE		*c1,
	O_CURVE		*c2,
	O_CURVE		*c3,
	O_CURVE		*c4)
{
	if (node_type(node) != CROSS_NODE) return;

	(void) printf("PRINTOUT OF CURVES AT CROSS NODE %llu\n\n",
		      node_number(node));
	print_node(node);
	if (c0 && c0->curve)
	{
		(void) printf("\t\tINCIDENT CURVE 0:\n");
		if (debugging("states"))
			show_curve_states(c0->curve);
		else
			print_o_curve(c0);
	}
	else
		(void) printf("\t\tINCIDENT CURVE 0 IS UNTRACKED:\n\n");
	if (c1 && c1->curve) 
	{
		(void) printf("\t\tREFLECTED CURVE 1:\n");
		if (debugging("states"))
			show_curve_states(c1->curve);
		else
			print_o_curve(c1);
	}
	else
		(void) printf("\t\tREFLECTED CURVE 1 IS UNTRACKED:\n\n");
	if (c2 && c2->curve) 
	{
		(void) printf("\t\tCONTACT CURVE 2:\n");
		if (debugging("states"))
			show_curve_states(c2->curve);
		else
			print_o_curve(c2);
	}
	else
		(void) printf("\t\tCONTACT CURVE 2 IS UNTRACKED:\n\n");
	if (c3 && c3->curve) 
	{
		(void) printf("\t\tREFLECTED CURVE 3:\n");
		if (debugging("states"))
			show_curve_states(c3->curve);
		else
			print_o_curve(c3);
	}
	else
		(void) printf("\t\tREFLECTED CURVE 3 IS UNTRACKED:\n\n");
	if (c4 && c4->curve)
	{
		(void) printf("\t\tINCIDENT CURVE 4:\n");
		if (debugging("states"))
			show_curve_states(c4->curve);
		else
			print_o_curve(c4);
	}
	else
		(void) printf("\t\tINCIDENT CURVE 4 IS UNTRACKED:\n\n");
	(void) printf("\nEND OF PRINTOUT OF CURVES AT CROSS NODE %llu\n",
	       node_number(node));
}		/*end print_cross_node*/
#endif /* defined(FULL_PHYSICS) */


/*ARGSUSED*/
EXPORT	void record_jet_velocity(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	COMPONENT	comp,original_comp;
	HYPER_SURF	*hs_on;
	FILE		*logfile = Output_file(data);
	INTERFACE	*intfc = front->interf;
	RECT_GRID	*gr = front->rect_grid;
	double		coords[MAXD], coords_on[MAXD];
	double		time = grid->time;
	int		max_step = stop_step(grid)+1;
	int		ymax = front->rect_grid->gmax[1];
	int		i;
	static Locstate s = NULL;
	static double	*record_time = NULL, *jet_velocity = NULL;
	static int	vel_indx = 0;

	if (s == NULL) 
	{
	    alloc_state(front->interf,&s,front->sizest);
	    vel_indx = 0;
	    uni_array(&record_time,max_step,FLOAT);
	    uni_array(&jet_velocity,max_step,FLOAT);
	}

	coords[0] = grid_center_coord(0,gr);
	coords[1] = cell_edge(ymax-1,1,gr);
	original_comp = component(coords,intfc);
	for (i = 2; i <= ymax; ++i) 
	{
	    coords[1] = cell_edge(ymax-i,1,gr);
	    comp = component(coords,intfc);
	    if (original_comp == comp) continue;
	    coords[1] += cell_width(ymax-i,1,gr);

	    nearest_intfc_state_and_pt(coords,original_comp,front,
	    	                       front,s,coords_on,&hs_on);

	    jet_velocity[vel_indx] = vel(1,s);
	    break;
	}
	record_time[vel_indx] = time;
	++vel_indx;
	if (about_to_stop == YES)
	{
	    (void) foutput(logfile);
	    (void) fprintf(logfile,"%-9s %-10s\n","TIME","JET V");
	    for (i = 0; i < vel_indx; ++i) 
	    {
	    	fprint_line_of_floats(logfile,2,record_time[i],jet_velocity[i]);
	    }
	    (void) fprintf(logfile,"\n");
	}
}		/*end record_jet_velocity*/

/*ARGSUSED*/
EXPORT	void multi_bubble_velocities(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	BOND		*b;
	CURVE		*c;
	HYPER_SURF	*hs,*hs_tmp;
	HYPER_SURF_ELEMENT *hse,*hse_tmp;
	FILE		*logfile = Output_file(data);
	INTERFACE	*intfc = front->interf;
	Locstate	sminl,sminr;
	POINT		*p_prev,*p_next,*p_prev_prev,*p_next_next;
	POINT		*p,*p_tmp;
	double		YU = front->rect_grid->U[1];
	double		YL = front->rect_grid->L[1];
	double		XL = front->rect_grid->L[0];
	double		XU = front->rect_grid->U[0];
	double		time = grid->time;
	double		xmin,ymin,vmin;
	int		i,j,num_bubble;
	static int     	num_alloc_bub = 0;
	static Locstate s = NULL;
	static HYPER_SURF **hsmin = NULL;
	static HYPER_SURF_ELEMENT **hsemin = NULL;
	static POINT	**pmin = NULL;

	debug_print("record_rs","Entered multi_bubble_velocities()\n");
	if (pmin == NULL)
	{
	    num_alloc_bub = 100;
	    if (logfile != stdout)
	    {
	    	(void) foutput(logfile);
	    	(void) fprintf(logfile,"%-9s %-10s %-10s %-10s\n",
				       "time","x_bubbles","y_bubbles",
				       "v_bubble");
	    }
	    alloc_state(front->interf,&s,front->sizest);
	    uni_array(&pmin,num_alloc_bub,sizeof(POINT *));
	    uni_array(&hsemin,num_alloc_bub,sizeof(HYPER_SURF_ELEMENT *));
	    uni_array(&hsmin,num_alloc_bub,sizeof(HYPER_SURF *));
	}

	i = 0;
	num_bubble = 0;
	(void) next_point(intfc,NULL,NULL,NULL);
	p_prev = NULL;
	p_next = NULL;
	p_prev_prev = NULL;
	p_next_next = NULL;
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE)
		continue;
	    c = Curve_of_hs(hs);
	    b = Bond_of_hse(hse);
	    if (c->num_points < 12)
		continue;
	    if (p == c->start->posn)
		continue;
	    if (p == c->end->posn)
		continue;
	    if (b->prev->start == c->start->posn)
		continue;
	    if (b->end == c->end->posn)
		continue;
	    if (Coords(p)[1] > 0.51*(YL+YU))
		continue;
	    p_prev = b->prev->start;
	    p_next = b->end;
	    p_prev_prev = b->prev->prev->start;
	    p_next_next = b->next->end;
	    if ((Coords(p)[1] <= Coords(p_prev)[1]) &&
	        (Coords(p)[1] <= Coords(p_next)[1]) &&
	        (Coords(p_prev)[1] < Coords(p_prev_prev)[1]) &&
	        (Coords(p_next)[1] < Coords(p_next_next)[1]))
	    {
		pmin[i] = p;
		hsemin[i] = hse;
		hsmin[i] = hs;
		++i;
		if(i >= num_alloc_bub)
		{
	            HYPER_SURF **old_hsmin = hsmin;
	            HYPER_SURF_ELEMENT **old_hsemin = hsemin;
	            POINT	**old_pmin = pmin;
		    int	old_num_alloc_bub = num_alloc_bub;
		    int j;

		    num_alloc_bub *= 2;
	            uni_array(&pmin,num_alloc_bub,sizeof(POINT*));
	            uni_array(&hsemin,num_alloc_bub,sizeof(HYPER_SURF_ELEMENT*));
	            uni_array(&hsmin,num_alloc_bub,sizeof(HYPER_SURF*));
		    for (j = 0; j < old_num_alloc_bub; ++j)
		    {
			pmin[j] = old_pmin[j];
			hsemin[j] = old_hsemin[j];
			hsmin[j] = old_hsmin[j];
		    }
		    free_these(3,old_pmin,old_hsemin,old_hsmin);
		}
		else
		    ++num_bubble;
	    }
	}
	for (i = 0; i < num_bubble-1; ++i)
	{
	    for (j = i+1; j < num_bubble; ++j) 
	    {
	        if (Coords(pmin[i])[0] > Coords(pmin[j])[0])
		{
		    p_tmp = pmin[i];
		    pmin[i] = pmin[j];
		    pmin[j] = p_tmp;
		    hse_tmp = hsemin[i];
		    hsemin[i] = hsemin[j];
		    hsemin[j] = hse_tmp;
		    hs_tmp = hsmin[i];
		    hsmin[i] = hsmin[j];
		    hsmin[j] = hs_tmp;
	        }
	    }
        }

	if (num_bubble == 0)
	{
	    debug_print("record_rs","Left multi_bubble_velocities()\n");
	    return;
	}

        if (Coords(pmin[0])[0] < 0.1*(XU-XL))
	{
	    p_tmp = pmin[0];
	    hse_tmp = hsemin[0];
	    hs_tmp = hsmin[0];
	    for (i = 1; i < num_bubble; ++i)
	    {
	    	pmin[i-1] = pmin[i];
	    	hsemin[i-1] = hsemin[i];
	    	hsmin[i-1] = hsmin[i];
	    }
	    pmin[num_bubble-1] = p_tmp;
	    hsemin[num_bubble-1] = hse_tmp;
	    hsmin[num_bubble-1] = hs_tmp;
	}

	(void) fprintf(logfile,"%-10g",time);

	for (i = 0; i < num_bubble; ++i)
	{
	    slsr(pmin[i],hsemin[i],hsmin[i],&sminl,&sminr);
	    ymin = Coords(pmin[i])[1];
	    xmin = Coords(pmin[i])[0];
	    vmin = ( vel(1,sminl) + vel(1,sminr) )*0.5;

	    (void) fprintf(logfile," %-10g %-10g %-10g",xmin,ymin,vmin);
	    if (i == num_bubble-1)
	    {
	        (void) fprintf(logfile,"\n");
	    }
	}

	(void) fflush(logfile);
	debug_print("record_rs","Left multi_bubble_velocities()\n");
}		/*end multi_bubble_velocities*/

/*
*			print_ramp_refl_stats():
*
*	This function is simply a driver for the output for ramp reflection
*	runs.  We just fork here for Mach or regular reflection.
*/

/*ARGSUSED*/
EXPORT void print_ramp_refl_stats(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	NODE		**n = front->interf->nodes;

	while (n && *n)
	{
		if (node_type(*n) == B_REFLECT_NODE)
		{
			print_reg_refl_node_stats(*n,front,grid,data);
			return;
		}
		else if (node_type(*n) == MACH_NODE)
		{
			print_Mach_node_stats(*n,front,grid,data);
			return;
		}
		++n;
	}
	return;		/* can't find appropriate node */

}		/*end print_ramp_refl_stats*/


/*
*		print_Mach_node_stats():
*
*	This function prints out diagnostics for analysis of Mach node
*	runs.  
*
*	TODO:  Add the output for secondary triple point and w1pp.  This is
*	already partially implemented.
*/

LOCAL void print_Mach_node_stats(
	NODE		*mnode,
	Front		*front,
	Grid		*grid,
	OUTPUT_DATA	*data)
{
	BOND		*b, *b_last;
	FILE		*logfile = Output_file(data);
	O_CURVE		Cmach, *newcmach = &Cmach;
	O_CURVE		Crefl, *newcrefl = &Crefl;
	RP_DATA		*RP = Rp_data(mnode);
	double		*base_posn;
	double		corner_posn[MAXD];
	double		*tp1_posn;		/* posn primary triple point */
	double		*tp2_posn_next;		/* current secondary posn */
	double		*scrds, *ecrds;
	double		refl_len;		/* base of mach to corner */
	double		bow_len;		/* bow to corner */
	double		nor[MAXD];
	double		node_ang;		/* angle of node velocity */
	double		wall_ang;		/* angle of ramp */
	double		chi;			/* node trajectory wrt ramp */
	double		w1p;			/* refl ang at node wrt node */
	double		base_pr;		/* pressure behind base mach */
	double		Mn;			/* Mach number at base mach */
	double		len = 0, shock_speed;
	double		next_ang, prev_ang;
	int		i, dim = front->rect_grid->dim;
	int		first_time_through_loop = YES;
	static Locstate state = NULL;
	static double	tp2_posn_prev[MAXD] = {0.0,0.0};
	static int	initial_call_step = 0;

	if (state == NULL)
	{
	    initial_call_step = front->step;
	    alloc_state(front->interf,&state,front->sizest);

	    (void) fprintf(logfile,"%10s %10s %10s %10s %10s %10s %10s",
			   "Time","Time_Step","Dt","TP_Posn[0]","TP_Posn[1]",
			   "Node_Vel[0]","Node_Vel[1]");

	    (void) fprintf(logfile,"%10s %10s %10s %10s %10s %10s %10s ",
			   "Refl_Len","Bow_Len","Inc_ang","Refl_ang",
			   "Cd_ang","Ms_ang","Mn");

	    (void) fprintf(logfile,"%10s %10s %10s %10s %10s %10s ",
			   "Chi","omega_1p","omega_0","Omega_1","beta_1",
			   "beta_s");

	    (void) fprintf(logfile,"%10s %10s %10s %10s %10s %10s %10s\n",
			   "Dens2","Dens3","Pressure3","Vel2[0]","Vel2[1]",
			   "Vel3[0]","Vel3[1]");
	}

	find_curve_with_status(mnode,&newcmach->curve,
			       &newcmach->orient,MACH_STEM);
	find_curve_with_status(mnode,&newcrefl->curve,
			       &newcrefl->orient,REFLECTED);

	ramp_reflection_corner_posn(corner_posn,NO,dim);
	if (corner_posn[0] == ERROR_FLOAT) return;

	tp1_posn = Coords(mnode->posn);
	base_posn = Coords(Opp_node_of_o_curve(newcmach)->posn);

	node_ang = angle(tp1_posn[0]-corner_posn[0],
			 tp1_posn[1]-corner_posn[1]);
	wall_ang = angle(base_posn[0]-corner_posn[0],
			 base_posn[1]-corner_posn[1]);

	chi = node_ang - wall_ang;

	nor[0] = cos(wall_ang);
	nor[1] = sin(wall_ang);
	if (RP->ang_dir == CLOCKWISE)
	{
	    if (newcmach->orient == POSITIVE_ORIENTATION)
	    	base_pr = pressure(left_end_state(newcmach->curve));
	    else
	    	base_pr = pressure(right_start_state(newcmach->curve));
	}
	else
	{
	    if (newcmach->orient == POSITIVE_ORIENTATION)
	    	base_pr = pressure(right_end_state(newcmach->curve));
	    else
	    	base_pr = pressure(left_start_state(newcmach->curve));
	}
	(void) s_polar_4(BEHIND_PRESSURE,base_pr,&shock_speed,nor,RP->state[0],
			 state,GAS_STATE);
	Mn = shock_speed/sound_speed(RP->state[0]);

	refl_len =
		distance_between_positions(tp1_posn,corner_posn,dim)*cos(chi);

	bow_len = distance_between_positions(corner_posn,
			     Coords(Opp_node_of_o_curve(newcrefl)->posn),dim);

	if (newcrefl->orient == POSITIVE_ORIENTATION)
	{
	    b = newcrefl->curve->first;
	    b_last = newcrefl->curve->last;
	    scrds = Coords(b->start);
	}
	else
	{
	    b = newcrefl->curve->last;
	    b_last = newcrefl->curve->first;
	    scrds = Coords(b->end);
	}

	tp2_posn_next = tp2_posn_prev;

	while (b != NULL)
	{
	    len += bond_length(b);

	    ecrds = (newcrefl->orient == POSITIVE_ORIENTATION) ?
	    	        Coords(b->end) : Coords(b->start);

	    next_ang = angle(ecrds[0]-scrds[0],ecrds[1]-scrds[1]);

	    if (first_time_through_loop)
	    {
	    	first_time_through_loop = NO;
	    	w1p = (node_ang + PI) - next_ang;
	    	prev_ang = next_ang;
	    }

	    /* The timestep test is to give the reflected wave a chance
	     * to smooth out after initialization.  The bond test is to
	     * avoid a kink that appears in the last bond of the reflected
	     * wave at the wall.  The final test is the tolerance used
	     * to identify a "valid" kink in the reflected shock.
	     */
	    if ((front->step - initial_call_step >= 15) && (b != b_last)
	        && (fabs(prev_ang - next_ang) > PI/36.0))
	    {
	    	for (i = 0; i < dim; ++i)
	    	    tp2_posn_prev[i] = tp2_posn_next[i];
		tp2_posn_next = scrds;
	    }

	    scrds = ecrds;
	    b = (newcrefl->orient == POSITIVE_ORIENTATION) ? b->next : b->prev;

	    prev_ang = next_ang;
	}

	len = distance_between_positions(tp2_posn_prev,tp2_posn_next,dim);
#if DONT_COMPILE
	if (len != 0.0) /*TODO see comment at top of function*/
	{
	    double		v[MAXD];
	    for (i = 0; i < dim; ++i)
	    {
	    	v[i] = (tp2_posn_prev[i] - tp2_posn_next[i]) / front->dt;
		tp2_posn_prev[i] = tp2_posn_next[i];
	    }
	    (void) fprintf(logfile,"Mtp2: ");
	    fprint_line_of_floats(logfile,4,tp2_posn_next[0],
				  tp2_posn_next[1],v[0],v[1]);
	    (void) fprintf(logfile,"\n");
	}
#endif /* DONT_COMPILE */
	fprint_line_of_floats(logfile,27,grid->time,(double)front->step,
			      front->dt,Coords(mnode->posn)[0],
			      Coords(mnode->posn)[1],Node_vel(mnode)[0],
			      Node_vel(mnode)[1],refl_len,bow_len,
			      degrees(RP->ang[0]),degrees(RP->ang[1]),
			      degrees(RP->ang[2]),degrees(RP->ang[3]),Mn,
			      degrees(chi),degrees(w1p),
			      degrees(RP->ang[0] - node_ang),
			      degrees(RP->ang[1] - RP->ang[0]),
			      degrees(RP->ang[2] - RP->ang[1]),
			      degrees(RP->ang[3] - RP->ang[2]),
			      Dens(RP->state[2]),Dens(RP->state[3]),
			      pressure(RP->state[3]),vel(0,RP->state[2]),
			      vel(1,RP->state[2]),
			      vel(0,RP->state[3]),vel(1,RP->state[3]));
}		/*end print_Mach_node_stats*/


/*
*			print_reg_refl_node_stats():
*
*	This function provides the diagnostic output for regular reflection
*	runs.
*/

LOCAL void print_reg_refl_node_stats(
	NODE		*rnode,
	Front		*front,
	Grid		*grid,
	OUTPUT_DATA	*data)
{
	FILE		*logfile = Output_file(data);
	O_CURVE		Crefl, *newcrefl = &Crefl;
	RP_DATA		*RP = Rp_data(rnode);
	double		corner_posn[MAXD];
	double		refl_len, bow_len;
	int		dim = front->rect_grid->dim;
	static boolean	first = YES;

	if (first == YES)
	{
	    first = NO;

	    (void) fprintf(logfile,"%10s %10s %10s %10s %10s %10s %10s ",
	    	           "Time","Time_Step","Dt","Refl_posn[0]",
	    	           "Refl_posn[1]","Node_Vel[0]","Node_Vel[1]");

	    (void) fprintf(logfile,"%10s %10s %10s %10s\n",
	    	           "Refl_Len","Bow_Len","Inc_ang","Refl_ang");
	}

	find_curve_with_status(rnode,&newcrefl->curve,
			       &newcrefl->orient,REFLECTED);

	ramp_reflection_corner_posn(corner_posn,NO,dim);

	refl_len = distance_between_positions(Coords(rnode->posn),
					      corner_posn,dim);

	bow_len = distance_between_positions(corner_posn,
			     Coords(Opp_node_of_o_curve(newcrefl)->posn),dim);

	fprint_line_of_floats(logfile,11,
			      grid->time,(double)front->step,front->dt,
			      Coords(rnode->posn)[0],Coords(rnode->posn)[1],
			      Node_vel(rnode)[0],Node_vel(rnode)[1],
			      refl_len,bow_len,
			      degrees(RP->ang[0]),degrees(RP->ang[1]));

}		/*end print_reg_refl_node_stats*/
#endif /* defined(TWOD) */
