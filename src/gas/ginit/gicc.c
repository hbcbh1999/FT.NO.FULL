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
*				gicc.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains initialization routines for contact-contact interactions.
*
*/

#if defined(FULL_PHYSICS) && defined(TWOD)
#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	void	init_CC_node(Gas_param**,RECT_GRID*,size_t,int,
			     Front*,INIT_DATA*);

enum { MAX_CC_CURVES = 5 };


EXPORT	void init_CC_interaction(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	Gas_param	*params[MAX_CC_CURVES];
	int		i,n_curves;

	(void) printf("Enter number of curves around the CC node: ");
	(void) Scanf("%d\n",&n_curves);

	/* Input gas parameters */

	for( i = 0; i< n_curves; i++)
	{
	    params[i] = init_eos_params(init,ip," with contact discontinuity",
			                YES);
	}

	init_CC_node(params,front->rect_grid,front->sizest,n_curves,front,init);

}		/*end init_CC_interaction*/

/* 
*			init_CC_node():
*
*	This routine initializes the states about the CC node.
*/	

LOCAL void init_CC_node(
	Gas_param	**params,
	RECT_GRID	*rect_grid,
	size_t		sizest,
	int		n_curves,
	Front		*front,
	INIT_DATA	*init)
{
	Locstate	st[MAX_CC_CURVES]; 
	CURVE		*cur;
	NODE		*ns, *ne;
	POINT		*p;
	double		angle[MAX_CC_CURVES];
	double 		*L = rect_grid->L, *U = rect_grid->U;
	double		n[MAXD];
	double		coords[MAXD];
	int		i;
	char		s[Gets_BUF_SIZE];

	debug_print("init","Entering init_CC_node()\n");

		/* Allocate storage for states around node */

	for(i = 0; i < n_curves; i++)
	{
	    (void) sprintf(s," for region %d",i);
	    prompt_for_ambient_state(comp_type(COMPA+i),params[i],s,front,init);
	    st[i] = Ambient(comp_type(COMPA+i));
	    screen("Is the flow in region %d constant (dflt = no): ",i);
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    	(void)SetConstantFlowRegion(COMPA+i,Ambient(comp_type(COMPA+i)),
					    front->interf);
	}

	/* Input the angles of the three contact discontinuity curves  */
	/* respect to the positive x axis                              */

	screen("\nEnter %d angles for the %d contact curves.\n",
							n_curves,n_curves);
	screen("The angles are measured respect to the positive x axis\n");
	for(i = 0; i < n_curves; i++)
	{
	    screen("\tEnter angle[%d] here: ",i);
	    (void) Scanf("%f\n",&angle[i]);
	}
	
	screen("\nTo specify an initial position for the \n");
	screen("CC node, enter the x and y positions\n");
	screen("in relative coordinates.\n");
	{
	    double posn_x = 0.25,posn_y = 0.5;
	    screen("\tEnter the x coordinate (default = 0.25): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
		(void) sscan_float(s,&posn_x);
	    screen("\tEnter the y coordinate (default = 0.5): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
		(void) sscan_float(s,&posn_y);
	    coords[0] = L[0]+posn_x*(U[0]-L[0]);
	    coords[0] = L[1]+posn_y*(U[1]-L[1]);
	    p = Point(coords);
	}

	ns = make_node(p);
	node_type(ns) = CC_NODE;

		/* Make three contact discontinuity curves */

	for(i=0; i<n_curves; i++)
	{
	    angle[i] = radians(angle[i]);
	    n[0] = cos(angle[i]);
	    n[1] = sin(angle[i]);
	    Check_return(
		intersect_ray_with_boundary(Coords(p),n,L,U,coords,2),
		init_CC_node)

	    ne = make_node(Point(coords));
	    node_type(ne) = DIRICHLET_NODE;
	    cur = make_curve(((COMPA+i) % n_curves) +1,COMPA+i,ns,ne);
	    wave_type(cur) = CONTACT;
	    start_status(cur) = end_status(cur) = INCIDENT;
	    ft_assign(left_start_state(cur),st[(i+1)%n_curves],sizest);
	    ft_assign(left_end_state(cur),st[(i+1)%n_curves],sizest);
	    ft_assign(right_start_state(cur),st[i],sizest);
	    ft_assign(right_end_state(cur),st[i],sizest);
	}

	debug_print("init","Leaving init_CC_node()\n");
}		/*end init_CC_node*/
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */
