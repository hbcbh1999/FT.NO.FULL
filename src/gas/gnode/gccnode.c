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
*				gccnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Contains the routine used in the first advance of nodes
*	corresponding to contact-contact interactions.
*/

#if defined(TWOD)

#define DEBUG_STRING    "ccnode"
#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	void	cc_set_uswssten0(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			         HYPER_SURF*,Front*,Wave*,double);
LOCAL	void	cc_set_uswssten1(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			         HYPER_SURF*,Front*,Wave*,double);
LOCAL	void	cc_set_uswssten2(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			         HYPER_SURF*,Front*,Wave*,double);
LOCAL	boolean	cc_identify_curves_at_node(NODE*,NODE*,O_CURVE**,O_CURVE**,
					   Front*);
LOCAL	void	cc_compute_states_at_crossing_bonds(NODE*,O_CURVE**,O_CURVE**,
						    BOND**,double*,Front*,
						    Wave*,double);
LOCAL	void	cc_normal(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,double*,
			  Front*);
LOCAL	void	cc_tangent(POINT*,BOND*,CURVE*,double*,Front*);
LOCAL	void	cc_node_point_propagate0(Front*,Wave*,POINT*,POINT*,
					 HYPER_SURF_ELEMENT*,HYPER_SURF*,
					 double,double*);
LOCAL	void	cc_node_point_propagate1(Front*,Wave*,POINT*,POINT*,
					 HYPER_SURF_ELEMENT*,HYPER_SURF*,
					 double,double*);
LOCAL	void	cc_node_point_propagate2(Front*,Wave*,POINT*,POINT*,
					 HYPER_SURF_ELEMENT*,HYPER_SURF*,
					 double,double*);
LOCAL	void	cc_compute_net_normal_at_node(double*,O_CURVE**,
					      O_CURVE**,Front*);
LOCAL	void	cc_propagate_old_node_position(NODE*,O_CURVE**,O_CURVE**,Front*,
					       Wave*,double);
LOCAL	void	cc_set_new_node_position(NODE*,O_CURVE**,O_CURVE**,
					 BOND**,POINT*,
				      double*,Front*,Wave*,double);
LOCAL	void	cc_reindex_to_minimum_angle_sector(O_CURVE**,O_CURVE**,BOND**);

#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void	fprint_cc_set_uswssten(USWSSten2d*);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

LOCAL	double	cc_nor[3];
LOCAL	double	cc_tan[3];
LOCAL	O_CURVE	**oldc = NULL;
LOCAL	O_CURVE	**newc = NULL;
LOCAL	TANGENT_FUNCTION _fr_tangent;
LOCAL	NORMAL_FUNCTION _fr_normal;

/*
*			cc_node_propagate():
*
*	                    curve 2
*				|
*				|
*				|
*				|
*	              state 2   |    state 1
*				|
*				|
*		       \       / \      /
*		       _\/   /     \  |/_
*			   /         \
*		         /   state 0   \
*		       /	         \
*	            curve 0            curve 1
*
*	A CC (Contact Collision) node is a node at which three contact
*	discontiuity curves  meet. Such nodes are created by the collision of
*	two blobs of material embedded in an external medium.  Two curves
*	bounding the colliding blobs have status at the node equal to
*	INCIDENT, the third curve has status SLIP.  Initially, the INCIDENT
*	curves correspond to the exterior sides of the colliding materials,
*	while the SLIP curve corresponds to the collision surface between the
*	two blobs.  Note: this algorithm assumes that the angle included
*	between curves 0 and 1 is smaller (generally much smaller) than the
*	angles of the 1,2 or 2,0 sectors.  As the node geometry
*	changes, the curves will be reindexed so that 0,1 sector is always
*	the smallest; see cc_reindex_to_minimum_angle_sector() for details.
*/

/*ARGSUSED*/
EXPORT int cc_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,	/* initialized to oldn in advance_front2d() */
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
        BOND		*newb[2];
	NODE		*oppn[2];
	POINT		*pc;
	double		tcr[2], newsep, oppsep;

	DEBUG_ENTER(cc_node_propagate)
#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG)
	{
	    (void) printf("\n\tOLD NODE:\n");
	    print_node(oldn);
	    (void) printf("Interface into cc_node_propagate\n");
	    print_interface(oldn->interface);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (oldc == NULL)
	{
	    bi_array(&oldc,3,1,sizeof(O_CURVE));
	    bi_array(&newc,3,1,sizeof(O_CURVE));
	}

	if (cc_identify_curves_at_node(oldn,newn,oldc,newc,fr) ==
	    FUNCTION_FAILED)
	{
	    (void) printf("WARNING in cc_node_propagate(), "
	                  "can't identify curves\n");
	    print_node(oldn);
	    (void) printf("returning ERROR_NODE");
	    DEBUG_LEAVE(cc_node_propagate)
	    return ERROR_NODE;
	}

	cc_propagate_old_node_position(oldn,oldc,newc,fr,wave,dt);

	set_ordinary_cross_only(flag);
	while (intersection_of_two_o_curves(newc[0],NULL,newc[1],NULL,
					    newb,newb+1,&pc,tcr,tcr+1,fr,
					    (POINTER)wave,dt,flag) == YES)
	{

	    /*
	    *  newb, newb+1 are the crossing bonds for newc[0],
	    *  newc[1] pc is the crossing point tcr, tcr+1 are fractional
	    *  distances along the crossing bonds (in the positive
	    *  direction ?)  to the crossing point.
	    *
	     *  Use `while' rather than `if' because there might be
	    *  several closely spaced crosses which should be removed
	    *  at the same time.
	     *
	     * TODO: failure to intersect may indicate the interfaces
	    * are pulling apart.  This case needs to be handled.
	    */

	    oppn[0] = (newc[0]->orient == POSITIVE_ORIENTATION) ?
	    		newc[0]->curve->end : newc[0]->curve->start ;
	    oppn[1] = (newc[1]->orient == POSITIVE_ORIENTATION) ?
	    		newc[1]->curve->end : newc[1]->curve->start ;
	    newsep = separation(pc,newn->posn,2);
	    oppsep = separation(pc,oppn[0]->posn,2);

	    if (oppn[0] != oppn[1] || newsep <= oppsep)
	    {
	    	cc_set_new_node_position(oldn,oldc,newc,newb,pc,tcr,fr,wave,dt);
	       	if (DEBUG)
	        {
	    	    (void) printf("\n\tNEW NODE:\n");
	    	    print_node(newn);
	    	    (void) printf("New interface after ");
	    	    (void) printf("cc_set_new_node_position()\n");
	    	    print_interface(newn->interface);
	       	}
	    }
	    else
	    {
	    	/* These two O_CURVES form a closed loop but the
	    	 * cross is closer to the other node.  The curves
	    	 * should be untangled working from there. */
	    	break;
	    }
	}

	propagation_status(newn) = PROPAGATED_NODE;

#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG) 
	{
	    int i;
	    for(i = 0; i< 3; ++i)
	    {
	    	(void) printf("\t\tnewc[%d]\n",i);
	    	print_o_curve(newc[i]);
	    	if (debugging("ccstates"))
	    	    show_curve_states(newc[i]->curve);
	    }
	    (void) printf("New interface after cc_node_propagate()\n");
	    print_interface(newn->interface);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	DEBUG_LEAVE(cc_node_propagate)
	return GOOD_NODE;
}		/*end cc_node_propagate*/


LOCAL	void	cc_propagate_old_node_position(
	NODE            *oldn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	Front		*fr,
	Wave		*wave,
	double		dt)
{
	INTERFACE       *old_intfc = oldn->interface;
	NODE		*newn = Node_of_o_curve(newc[2]);
	Locstate	s0, s1;
	POINT		*p;
	RP_DATA		*RP = Rp_data(newn);
	boolean		sav_interpolate;
	double		n[MAXD];
	double		W[MAXD];
	int		i, dim = fr->rect_grid->dim;
	static	void	(*cc_node_point_propagate[])(Front*,Wave*,POINT*,
				POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				double,double*) = {
					cc_node_point_propagate0,
					cc_node_point_propagate1,
					cc_node_point_propagate2};

	DEBUG_ENTER(cc_propagate_old_node_position)

	cc_compute_net_normal_at_node(n,oldc,newc,fr);

	/* Propagate node using oldc[2] and net tangent, normal uni_arrays */
	_fr_normal = interface_normal_function(old_intfc);
	set_cc_normal_function(&interface_normal_function(old_intfc));
	_fr_tangent = interface_tangent_function(old_intfc);
	set_cc_tangent_function(&interface_tangent_function(old_intfc));

	sav_interpolate = interpolate_intfc_states(newn->interface);
	interpolate_intfc_states(newn->interface) = NO;

	/*
	 * Propagate the node once for each of the three curves using the net
	 * normal/tangent uni_arrays obtained in cc_compute_net_normal_at_node().
	 * For curves 0,1 (status INCIDENT prior to reindexing) only state
	 * data is updated; the updated node position is determined solely by
	 * the node propagation for curve 2 (status  SLIP prior to reindexing).
	 */

	for (i = 0; i < 3; ++i)
	{
	    p = (i < 2) ? Point(NULL) : newn->posn;
	    (*cc_node_point_propagate[i])(fr,wave,
	    	Node_of_o_curve(oldc[i])->posn,p,
	    	Hyper_surf_element(Bond_at_node_of_o_curve(oldc[i])),
	    	Hyper_surf(oldc[i]->curve),dt,W);
	    copy_state(Left_state_at_node_of_o_curve(newc[i]),left_state(p));
	    copy_state(Right_state_at_node_of_o_curve(newc[i]),right_state(p));
	}
	/* Restore prior values of intfc/front functions */
	interpolate_intfc_states(newn->interface) = sav_interpolate;
	interface_normal_function(old_intfc) = _fr_normal;
	interface_tangent_function(old_intfc) = _fr_tangent;
	/* On exit above loop, p = newn->posn */
	for (i = 0; i < dim; ++i)
		Node_vel(newn)[i] = W[i];
	/* Set RP states of new node */
	RP->ang_dir = COUNTER_CLOCK;
	RP->stype = state_type(left_state(p));
	if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[2]->orient))
	{
	    set_state(RP->state[0],RP->stype,right_state(p));
	    set_state(RP->state[2],RP->stype,left_state(p));
	}
	else
	{
	    set_state(RP->state[0],RP->stype,left_state(p));
	    set_state(RP->state[2],RP->stype,right_state(p));
	}
	/* set new node RP state for sector between curves 0,1 to */
	/* AVERAGE of these two solutions. */
	s0 = (curve_ang_oriented_l_to_r(RP->ang_dir,newc[0]->orient)) ?
		Right_state_at_node_of_o_curve(newc[0]) :
		Left_state_at_node_of_o_curve(newc[0]);
	s1 = (curve_ang_oriented_l_to_r(RP->ang_dir,newc[1]->orient)) ?
		Left_state_at_node_of_o_curve(newc[1]) :
		Right_state_at_node_of_o_curve(newc[1]);
	interpolate_states(fr,0.5,0.5,Coords(p),s0,Coords(p),s1,RP->state[1]);
	/* Finally, update curve 0,1 states */
	/* to agree with RP data at new node. */
	for (i = 0; i < 2; ++i)
	{
		if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[i]->orient))
		{
			set_state(Left_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[i]);
			set_state(Right_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[i+1]);
		}
		else
		{
			set_state(Left_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[i+1]);
			set_state(Right_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[i]);
		}
	}

	DEBUG_LEAVE(cc_propagate_old_node_position)
}	/*end cc_propagate_old_node_position*/


/*
*			cc_compute_net_normal_at_node():
*
*	Computes a net unit normal vector at the outgoing interface of the 
*	CC node triple point.  The normal is computed by rotating
*	a net tangent vector formed as indicated in the diagram below.
*
*	                                  O
*	                               * /
*                                    *  /
*	                           *   /
*	  direction of tngt[1]   *    /bdir[1]
*	                       *     /
*	                     *      / oldc[1]->curve
*	                   *       /
*                        *bdir[2] /
*	oldc[2]->curve O---------*------------
*                        *        \ 
*	                   *       \ 
*	                     *      \  oldc[0]->curve
*	                       *     \ 
*	                         *    \bdir[0] 
*          direction of tngt[0]    *   \
*	                             *  \
*	                               * \
*	                                 *O
*
*	In the diagram,  bdir[i] is the bond tangent to the respective
*	curves,  oldc[i]->curve, as returned by the function
*	bond_tangent_to_curve.  The net tangents, tngt[i], are obtained
*	by joining the opposite (with respect to the node) endpoints of
*	the bonds bdir[i] consistent with the orientation of oldc[2]->curve.
*	In effect,  tngt[i] is the tangent that would be computed at
*	the position of the node were the curves oldc[2]->curve and
*	oldc[i]->curve joined.  The net tangent to oldc[2]->curve is the
*	bisector of the the two tangents tngt[0] and tngt[1].
*	
*	NOTE: Ordering of curves: (COUNTER_CLOCK: 0,1,2)
*	      Prior to reindexing, curve stati are: INCIDENT, INCIDENT, SLIP.
*
*	NOTE: By construction, the PROJECTION of tgnt[0], tgnt[1] onto oldc[2]
*	at the node is POSITIVE when taking into account the orientation of
*	the latter.
*/

LOCAL	void cc_compute_net_normal_at_node(
	double	*n,
	O_CURVE	**oldc,
	O_CURVE	**newc,
	Front	*fr)
{
	POINT	*p = Node_of_o_curve(oldc[2])->posn;
	double	*h = fr->rect_grid->h;
	double	magv, blen;
	int	i, j, dim = fr->rect_grid->dim;
	static	double	**tngt = NULL;
	static	BOND	*b[3] = {NULL, NULL, NULL};

	DEBUG_ENTER(cc_compute_net_normal_at_node)
	if (tngt == NULL)
	{
		static BOND *B = NULL;
		bi_array(&tngt,2,MAXD,FLOAT);
		uni_array(&B,3,sizeof(BOND));
		for (i = 0; i < 3; ++i)
		{
			b[i] = B+i;
			b[i]->start = Static_point(fr->interf);
			b[i]->end = Static_point(fr->interf);
		}
	}

	for (blen = 0, i = 0; i < 3; ++i)
	{
		double slen;

		slen = scaled_bond_length(Bond_at_node_of_o_curve(oldc[i]),
					  h,dim);
		blen = max(blen,slen);
	}

#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG)
		(void) printf("blen = %g\n",blen);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	for (i = 0; i < 3; ++i)
	{
		bond_secant_to_curve(p,Bond_at_node_of_o_curve(oldc[i]),
			oldc[i]->curve,oldc[i]->orient,b[i],fr,blen);
	}

	/*
	* Determine the sector of minimum angle and reindex the
	* curves so that curves 0 and 1 bound this sector.
	*/

	cc_reindex_to_minimum_angle_sector(oldc,newc,b);


	for (i = 0; i < 2; ++i)
	{
	    if (oldc[2]->orient == POSITIVE_ORIENTATION)
	    {
	        if (oldc[i]->orient == POSITIVE_ORIENTATION)
		{
		    for (j = 0; j < dim; ++j)
			tngt[i][j] = Coords(b[2]->end)[j] -
				     Coords(b[i]->end)[j];
		}
	        else
		{
		    for (j = 0; j < dim; ++j)
	    	        tngt[i][j] = Coords(b[2]->end)[j] -
				     Coords(b[i]->start)[j];
		}
	    }
	    else
	    {
	        if (oldc[i]->orient == POSITIVE_ORIENTATION)
		{
	    	    for (j = 0; j < dim; ++j)
	    	        tngt[i][j] = Coords(b[i]->end)[j] -
				     Coords(b[2]->start)[j];
		}
	        else
		{
		    for (j = 0; j < dim; ++j)
	    	        tngt[i][j] = Coords(b[i]->start)[j] -
				     Coords(b[2]->start)[j];
		}
	    }
	    magv = mag_vector(tngt[i],dim);
#if defined(DEBUG_NODE_PROPAGATE)
	    if (DEBUG)
	    {
	       (void) printf("Before normalization, ");
	       (void) printf("tngt[%d] = <%g, %g>, magv = %g\n",
			     i,tngt[i][0],tngt[i][1],magv);
			
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    for (j = 0; j < dim; ++j)
	        tngt[i][j] /= magv;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (DEBUG)
	    {
		(void) printf("After normalization, tngt[%d] = <%g, %g>\n",
			i,tngt[i][0],tngt[i][1]);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	}

	for (i = 0; i < dim; ++i)
		cc_tan[i] = 0.5*(tngt[0][i]+tngt[1][i]);
	magv = mag_vector(cc_tan,dim);
	for (i = 0; i < dim; ++i)
		cc_tan[i] /= magv;
	/* TWO DIMENSIONS ONLY*/
	n[0] = cc_nor[0] =  cc_tan[1];
	n[1] = cc_nor[1] = -cc_tan[0];
	DEBUG_LEAVE(cc_compute_net_normal_at_node)
}	/*end cc_compute_net_normal_at_node*/

/*
*			cc_reindex_to_minimum_angle_sector():
*
*	Reindex the curve lists so that curve 2 is opposite the sector with
*	smallest included angle (the 0--1 sector), subject to the proviso below.
*/

LOCAL	void	cc_reindex_to_minimum_angle_sector(
	O_CURVE	**oldc,
	O_CURVE	**newc,
	BOND	**b)
{
	double	ang[3], ang_min;
	double	t[3][2];
	int	i, imin;

	DEBUG_ENTER(cc_reindex_to_minimum_angle_sector)
	for (i = 0; i < 3; ++i)
	{
		double mag;
		if (oldc[i]->orient == POSITIVE_ORIENTATION)
		{
		       t[i][0] = Coords(b[i]->end)[0] - Coords(b[i]->start)[0];
		       t[i][1] = Coords(b[i]->end)[1] - Coords(b[i]->start)[1];
		}
		else
		{
		       t[i][0] = Coords(b[i]->start)[0] - Coords(b[i]->end)[0];
		       t[i][1] = Coords(b[i]->start)[1] - Coords(b[i]->end)[1];
		}
		mag = hypot(t[i][0],t[i][1]);
		t[i][0] /= mag;
		t[i][1] /= mag;
	}

	ang_min = HUGE_VAL;
	imin = ERROR;
	for (i = 0; i < 3; ++i)
	{
		double sp, vp;
		sp = t[i][0]*t[(i+1)%3][0] + t[i][1]*t[(i+1)%3][1];
		vp = t[i][0]*t[(i+1)%3][1] - t[i][1]*t[(i+1)%3][0];
		ang[i] = angle(sp,vp);
#if defined(DEBUG_NODE_PROPAGATE)
		if (DEBUG)
		{
			(void) printf("t[%d] = <%g, %g>, ang[%d] = %g\n",
				i,t[i][0],t[i][1],i,ang[i]);
		}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		if (ang[i] < ang_min)
		{
			ang_min = ang[i];
			imin = i;
		}
	}
#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG)
		(void) printf("imin = %d\n",imin);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	/*
	 * The generic situation initially following the collision of two
	 * blobs with locally smooth and convex boundaries (near the initial
	 * contact point) is two cc_nodes each with a fairly narrow 0,1
	 * (INCIDENT / INCIDENT) sector.  As penetration proceeds the 0,1
	 * angle narrows until the curves cross.  The node is then moved to
	 * the location of the cross which also opens the 0,1 angle and the
	 * process repeats.  Eventually, backflow may begin to develop in
	 * either or both blobs.  When this happens the 0,1 angle may widen
	 * to PI or even more.  At some point the algorithm must switch to
	 * the (current) minimum angle sector.  However, to avoid
	 * "chattering" between different sectors and erratic node motion,
	 * the original sector is retained until the backflow patern is
	 * clearly established per the condition below.  The choice of PI/3
	 * is somewhat arbitrary but suggested by observation of rod
	 * penetration runs.
	 */

	if ((imin != 0) && (ang_min < PI/3.0)) /*TOLERANCE*/
	{
		O_CURVE	Oldc_tmp[3], Newc_tmp[3];
		BOND	*b_tmp[3];

		/* Reindex segments */
		for (i = 0; i < 3; ++i)
		{
			Oldc_tmp[i] = *oldc[i];
			Newc_tmp[i] = *newc[i];
			b_tmp[i] = b[i];
		}
		for (i = 0; i < 3; ++i)
		{
			*oldc[i] = Oldc_tmp[(i+imin)%3];
			*newc[i] = Newc_tmp[(i+imin)%3];
			b[i] = b_tmp[(i+imin)%3];
		}
#if defined(DEBUG_NODE_PROPAGATE)
		if (DEBUG) 
		{
			(void) printf("NOTICE - reindexing curves\n");
			(void) printf("Curves after index shift\n");
			for(i = 0; i < 3; ++i)
			{
				(void) printf("\n\t\toldc[%d]\n",i);
				print_o_curve(oldc[i]);
			}
			for(i = 0; i < 3; ++i)
			{
				(void) printf("\n\t\tnewc[%d]\n",i);
				print_o_curve(newc[i]);
			}
		}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	}
	DEBUG_LEAVE(cc_reindex_to_minimum_angle_sector)
}	/*end cc_reindex_to_minimum_angle_sector*/


LOCAL	boolean	cc_identify_curves_at_node(
	NODE	*oldn,
	NODE	*newn,
	O_CURVE	**oldc,
	O_CURVE	**newc,
	Front	*fr)
{
	int i;

	DEBUG_ENTER(cc_identify_curves_at_node)

	/*
	 * Get the first (should be the only) oriented curve (O_CURVE)
	 * attached to the node which has status SLIP.  Both outgoing
	 * (POSITIVE_ORIENTATION = 1) and incoming (NEGATIVE_ORIENTATION = 2)
	 * lists of curves at node are searched.  Store the curve and
	 * orientation as the last elements (i.e., index = 2) in the curve
	 * lists for both new/old nodes.
	 */

	find_curve_with_status(oldn,&oldc[2]->curve,&oldc[2]->orient,SLIP);
	find_curve_with_status(newn,&newc[2]->curve,&newc[2]->orient,SLIP);
	if (!oldc[2]->curve)
	{
		(void) printf("WARNING in cc_identify_curves_at_node(), ");
		(void) printf("no slip curve\n");
		DEBUG_LEAVE(cc_identify_curves_at_node)
		return FUNCTION_FAILED;
	}

	/* Identify the other two (INCIDENT) curves around the old node */
	/* and order them (COUNTER_CLOCK: 0,1,2):  INCIDENT, INCIDENT, SLIP */

	oldc[0]->curve = adjacent_curve(oldc[2]->curve,oldc[2]->orient,
					COUNTER_CLOCK,&oldc[0]->orient);
	oldc[1]->curve = adjacent_curve(oldc[2]->curve,oldc[2]->orient,
					CLOCKWISE,&oldc[1]->orient);

	if ( (status_at_node(oldc[0]->curve,oldc[0]->orient) != INCIDENT)
	      ||
	     (status_at_node(oldc[1]->curve,oldc[1]->orient) != INCIDENT)
	)
	{
		(void) printf("WARNING in cc_identify_curves_at_node(), ");
		(void) printf("can't find incident curves\n");
		DEBUG_LEAVE(cc_identify_curves_at_node)
		return FUNCTION_FAILED;
	}

	/* Load newc->curve with curve corresponding to oldc->curve */
	/* and set newc->orient to match oldc->orient. */

	for (i = 0; i < 2; ++i)
	{
		if (!find_correspond_of_oriented_curve(oldc[i],
				newc[i],newn,fr,newn->interface))
		{
			(void) printf("WARNING in cc_identify_curves_at_node(), ");
			(void) printf("can't find correspond curve\n");
			DEBUG_LEAVE(cc_identify_curves_at_node)
			return FUNCTION_FAILED;
		}
	}

#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG) 
	{
		(void) printf("Curves before index shift\n");
		for(i = 0; i < 3; ++i)
		{
			(void) printf("\n\t\toldc[%d]\n",i);
			print_o_curve(oldc[i]);
			if (debugging("ccstates"))
				show_curve_states(oldc[i]->curve);
		}
		for(i = 0; i < 3; ++i)
		{
			(void) printf("\n\t\tnewc[%d]\n",i);
			print_o_curve(newc[i]);
		}
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	DEBUG_LEAVE(cc_identify_curves_at_node)
	return FUNCTION_SUCCEEDED;
}	/*end cc_identify_curves_at_node*/


LOCAL	void	cc_set_new_node_position(
	NODE	*oldn,
	O_CURVE	**oldc,
	O_CURVE	**newc,
	BOND	**newb,
	POINT	*pc,
	double	*tcr,
	Front	*fr,
	Wave	*wave,
	double	dt)
{
	INTERFACE	*intfc = newc[2]->curve->interface;
			/* start or end node per curve orientation */ 
	NODE		*newn = Node_of_o_curve(newc[2]);
	POINT		*p;
	RP_DATA		*RP = Rp_data(newn); /* init RP data structure */
	boolean		sav_interpolate = interpolate_intfc_states(intfc);
	boolean 	sav_copy = copy_intfc_states();
	Locstate	s0, s1;
	int		i, dim = fr->rect_grid->dim;

	DEBUG_ENTER(cc_set_new_node_position)
	/* set p to newn->posn w/o copying states */
	set_copy_intfc_states(NO);
	p = copy_point(newn->posn);
	set_copy_intfc_states(sav_copy);

	interpolate_intfc_states(intfc) = NO;
	/*
	 * Insert a new POINT on the SLIP curve,
	 * coincident with the position of the new node.
	 * Ref: i_insert_point_in_bond(), f_user_2d_insert_point_in_bond()
	 */
	insert_point_adjacent_to_node(p,newc[2]->curve,newc[2]->orient);
	interpolate_intfc_states(intfc) = sav_interpolate;
	/* copy curve state data to p */
	set_state(left_state(p),RP->stype,
		Left_state_at_node_of_o_curve(newc[2]));
	set_state(right_state(p),RP->stype,
		Right_state_at_node_of_o_curve(newc[2]));

	cc_compute_states_at_crossing_bonds(oldn,oldc,newc,newb,tcr,fr,wave,dt);

	for (i = 0; i < dim; ++i)
	{
		Coords(newn->posn)[i] = Coords(pc)[i];
		if (dt > 0.0)
			Node_vel(newn)[i] =
			    (Coords(newn->posn)[i] - Coords(oldn->posn)[i])/dt;
	}
	for (i = 0; i < 3; ++i)
		set_bond_length(Bond_at_node_of_o_curve(newc[i]),dim);

	if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[0]->orient))
	{
		s0 = Right_state_at_node_of_o_curve(newc[0]);
		set_state(RP->state[0],RP->stype,
			Left_state_at_node_of_o_curve(newc[0]));
	}
	else
	{
		s0 = Left_state_at_node_of_o_curve(newc[0]);
		set_state(RP->state[0],RP->stype,
			Right_state_at_node_of_o_curve(newc[0]));
	}
	if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[1]->orient))
	{
		s1 = Left_state_at_node_of_o_curve(newc[1]);
		set_state(RP->state[2],RP->stype,
			Right_state_at_node_of_o_curve(newc[1]));
	}
	else
	{
		s1 = Right_state_at_node_of_o_curve(newc[1]);
		set_state(RP->state[2],RP->stype,
			Left_state_at_node_of_o_curve(newc[1]));
	}
	interpolate_states(fr,0.5,0.5,Coords(newn->posn),s0,
		Coords(newn->posn),s1,RP->state[1]);
	for (i = 0; i < 3; ++i)
	{
		if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[i]->orient))
		{
			set_state(Left_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[i]);
			set_state(Right_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[(i+1)%3]);
		}
		else
		{
			set_state(Left_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[(i+1)%3]);
			set_state(Right_state_at_node_of_o_curve(newc[i]),
				RP->stype,RP->state[i]);
		}
	}
	DEBUG_LEAVE(cc_set_new_node_position)
}	/*end cc_set_new_node_position*/


LOCAL	void	cc_compute_states_at_crossing_bonds(
	NODE    *oldn,
	O_CURVE	**oldc,
	O_CURVE	**newc,
	BOND	**newb, /* The crossing bonds for newc[0,1] */
	double	*tcr,
	Front	*fr,
	Wave    *wave,
	double	dt)
{
	INTERFACE *old_intfc = oldn->interface;
	BOND	*oldb, *b;
	double	W[MAXD];
	int	i;
	static	POINT	*ps = NULL, *pe = NULL;
	static	void	(*cc_node_point_propagate[])(Front*,Wave*,POINT*,
				POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				double,double*) = {
					cc_node_point_propagate0,
					cc_node_point_propagate1};

	DEBUG_ENTER(cc_compute_states_at_crossing_bonds)
	if (fr->_point_propagate == UnsplitPointPropagate(fr))
	{
	    for (i = 0; i < 2; ++i)
	    {
	        left_state_along_bond(tcr[i],newb[i],newc[i]->curve,
	            Left_state_at_node_of_o_curve(newc[i]));
	        right_state_along_bond(tcr[i],newb[i],newc[i]->curve,
	            Right_state_at_node_of_o_curve(newc[i]));
	        while (Point_adjacent_to_node(newc[i]->curve,
					      newc[i]->orient) !=
		       Point_of_bond(newb[i],Opposite_orient(newc[i]->orient)))
			(void) delete_point_adjacent_to_node(fr,newc[i]->curve,
							     newc[i]->orient);
	    }
	    DEBUG_LEAVE(cc_compute_states_at_crossing_bonds)
	    return;
	}

	if (ps == NULL)
	{
		ps = Static_point(fr->interf);
		pe = Static_point(fr->interf);
	}

	/*
	* When using the operator split point propagate,  it is
	* necessary to recompute the states on the crossing bonds
	* using the unsplit point propagate operator.  This is to
	* ensure consistency between the unsplit propagation method
	* used at the old node,  and the operator splitting used
	* to propagate interior points on the curves.
	*/

	_fr_normal = interface_normal_function(old_intfc);
	set_cc_normal_function(&interface_normal_function(old_intfc));
	_fr_tangent = interface_tangent_function(old_intfc);
	set_cc_tangent_function(&interface_tangent_function(old_intfc));
	for (i = 0; i < 2; ++i)
	{
		for (oldb = Bond_at_node_of_o_curve(oldc[i]),
				b = Bond_at_node_of_o_curve(newc[i]);
				oldb != NULL && b != NULL;
				oldb = Following_bond(oldb,oldc[i]->orient),
				b = Following_bond(b,newc[i]->orient))
			if (b == newb[i]) break;
		if (b == NULL || oldb == NULL)
		{
		    screen("ERROR in cc_compute_states_at_crossing_bonds(), ");
		    screen("inconsistent old and new curves\n");
		    clean_up(ERROR);
		}
		(*cc_node_point_propagate[i])(fr,wave,oldb->start,ps,
					      Hyper_surf_element(oldb),
					      Hyper_surf(oldc[i]->curve),dt,W);
		(*cc_node_point_propagate[i])(fr,wave,oldb->end,pe,
					      Hyper_surf_element(oldb),
					      Hyper_surf(oldc[i]->curve),dt,W);
	        interpolate_states(fr,1.0-tcr[i],tcr[i],Coords(ps),
		                   left_state(ps),Coords(pe),left_state(pe),
		                   Left_state_at_node_of_o_curve(newc[i]));
	        interpolate_states(fr,1.0-tcr[i],tcr[i],Coords(ps),
		                   right_state(ps),Coords(pe),right_state(pe),
		                   Right_state_at_node_of_o_curve(newc[i]));
		/*
		 * Loop outward from the node deleting points (and
		 * bonds) until the bond at the node IS the crossing
		 * bond.
		 */
	        while (Point_adjacent_to_node(newc[i]->curve,
					      newc[i]->orient) !=
		       Point_of_bond(newb[i],Opposite_orient(newc[i]->orient)))
			(void) delete_point_adjacent_to_node(fr,newc[i]->curve,
							     newc[i]->orient);
	}
	interface_normal_function(old_intfc) = _fr_normal;
	interface_tangent_function(old_intfc) = _fr_tangent;
	DEBUG_LEAVE(cc_compute_states_at_crossing_bonds)
}	/*end cc_compute_states_at_crossing_bonds*/


LOCAL	void cc_node_point_propagate0(
	Front			*fr,
	Wave			*wave,
	POINT			*oldp,
	POINT			*newp,
	HYPER_SURF_ELEMENT	*oldhse,
	HYPER_SURF		*oldhs,
	double			dt,
	double			*V)
{
	void (*sav_set_uswssten)(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			         HYPER_SURF*,Front*,Wave*,double);

	DEBUG_ENTER(cc_node_point_propagate0)
	sav_set_uswssten = SetUSWSStencil(fr);
	SetUSWSStencil(fr) = cc_set_uswssten0;
	unsplit_point_propagate(fr,(POINTER)wave,oldp,newp,oldhse,oldhs,dt,V);
	SetUSWSStencil(fr) = sav_set_uswssten;
	DEBUG_LEAVE(cc_node_point_propagate0)
}	/*end cc_node_point_propagate0*/


LOCAL	void cc_node_point_propagate1(
	Front			*fr,
	Wave			*wave,
	POINT			*oldp,
	POINT			*newp,
	HYPER_SURF_ELEMENT	*oldhse,
	HYPER_SURF		*oldhs,
	double			dt,
	double			*V)
{
	void (*sav_set_uswssten)(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			         HYPER_SURF*,Front*,Wave*,double);

	DEBUG_ENTER(cc_node_point_propagate1)
	sav_set_uswssten = SetUSWSStencil(fr);
	SetUSWSStencil(fr) = cc_set_uswssten1;
	unsplit_point_propagate(fr,(POINTER)wave,oldp,newp,oldhse,oldhs,dt,V);
	SetUSWSStencil(fr) = sav_set_uswssten;
	DEBUG_LEAVE(cc_node_point_propagate1)
}	/*end cc_node_point_propagate1*/


LOCAL	void cc_node_point_propagate2(
	Front			*fr,
	Wave			*wave,
	POINT			*oldp,
	POINT			*newp,
	HYPER_SURF_ELEMENT	*oldhse,
	HYPER_SURF		*oldhs,
	double			dt,
	double			*V)
{
	void (*sav_set_uswssten)(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			         HYPER_SURF*,Front*,Wave*,double);

	DEBUG_ENTER(cc_node_point_propagate2)
	sav_set_uswssten = SetUSWSStencil(fr);
	SetUSWSStencil(fr) = cc_set_uswssten2;
	unsplit_point_propagate(fr,(POINTER)wave,oldp,newp,oldhse,oldhs,dt,V);
	SetUSWSStencil(fr) = sav_set_uswssten;
	DEBUG_LEAVE(cc_node_point_propagate2)
}	/*end cc_node_point_propagate2*/

/*
*			cc_set_uswssten():
*
*	Function of cc_set_uswssten? depends crucially on the
*	COUNTER_CLOCK ordering of the curves established in
*	cc_identify_curves_at_node().
*
*	Explanation: Location of points in a USWSSten2d depend on curve
*	orientation.  NEGATIVE (tangential) index means toward the start of the
*	curve RELATIVE TO ITS ORIENTATION; POSITIVE means the opposite.  This
*	must be accounted for when mixing points from stencils for different
*	curves as is done here.
*
*	The state data from the sector OPPOSITE the PRINCIPAL curve (ie, the
*	curve being propagated) is ignored.  We assume that this sector either
*	contains void or is sufficiently narrow so that state changes across
*	it may be ignored for purposes of node propagation.
*/

LOCAL	Tan_stencil	*tmp_sten = NULL;

LOCAL	void	cc_set_uswssten0(
	USWSSten2d		*uswssten,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	Front			*fr,
	Wave			*wave,
	double			dt)
{
	CURVE		*c = Curve_of_hs(hs);/*ASSUMES TWOD*/
	BOND		*b = Bond_of_hse(hse);/*ASSUMES TWOD*/
	Tan_stencil	*on_fr = tan_fr(uswssten)[0];
	int		i, j, dim = fr->rect_grid->dim;
	int		tan_rad = uswssten->tan_rad;
	size_t		sizest = fr->sizest;

	DEBUG_ENTER(cc_set_uswssten0)
	if (tmp_sten == NULL)
	    tmp_sten = alloc_tan_stencil(fr,tan_rad);

	uswssten->fr = fr;
	uswssten->wave = wave;
	uswssten->hs = hs;
	set_uswssten_geometry_and_center(uswssten,p,hse,hs,fr,dt);

	/* Use single normal ... uswssten is FLAT. */
	for (i = -tan_rad; i <= tan_rad; ++i)
	    for (j = 0; j < dim; ++j)
		nor_vec(uswssten)[i][j] = nor_vec(uswssten)[0][j];

	/* Curve 0 is the PRINCIPAL curve (the one being propagated). */

	switch (oldc[0]->orient)
	{
	case NEGATIVE_ORIENTATION:
	    states_at_distance_along_curve(p,b,c,NEGATIVE_ORIENTATION,
					   uswssten->ds,tan_rad,
					   on_fr->leftst-1,on_fr->rightst-1,
			                   on_fr->hs-1,on_fr->hse-1,
					   on_fr->t-1,on_fr->p-1,fr);
	    if (oldc[2]->orient == POSITIVE_ORIENTATION)
	    {
	    	states_at_distance_along_curve(p,
			       Bond_at_node_of_o_curve(oldc[2]),oldc[2]->curve,
			       POSITIVE_ORIENTATION,uswssten->ds,tan_rad,
			       on_fr->leftst+1,tmp_sten->rightst+1,
			       on_fr->hs+1,on_fr->hse+1,on_fr->t+1,
			       on_fr->p+1,fr);
	    }
	    else
	    {
	    	states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[2]),
					       oldc[2]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
			                       tmp_sten->leftst-1,
					       tmp_sten->rightst-1,
			                       tmp_sten->hs-1,
					       tmp_sten->hse-1,tmp_sten->t-1,
			                       tmp_sten->p-1,fr);
	    	for (i = 1; i <= tan_rad; ++i)
	    	{
	    	    ft_assign(on_fr->leftst[i],tmp_sten->rightst[-i],sizest);
		    on_fr->hs[i] = tmp_sten->hs[-i];
		    on_fr->hse[i] = tmp_sten->hse[-i];
		    on_fr->t[i] = tmp_sten->t[-i];
		    for (j = 0; j < dim; ++j)
		    	Coords(on_fr->p[i])[j] = Coords(tmp_sten->p[-i])[j];
		}
	    }
	    states_at_distance_along_curve(p,b,c,POSITIVE_ORIENTATION,
					   uswssten->ds,tan_rad,
			                   tmp_sten->leftst+1,on_fr->rightst+1,
			                   tmp_sten->hs+1,tmp_sten->hse+1,
					   tmp_sten->t+1,tmp_sten->p+1,fr);
	    break;
	case POSITIVE_ORIENTATION:
	    states_at_distance_along_curve(p,b,c,POSITIVE_ORIENTATION,
					   uswssten->ds,tan_rad,
					   on_fr->leftst+1,on_fr->rightst+1,
			                   on_fr->hs+1,on_fr->hse+1,
					   on_fr->t+1,on_fr->p+1,fr);
	    if (oldc[2]->orient == NEGATIVE_ORIENTATION)
	    {
	        states_at_distance_along_curve(p,
	    	                               Bond_at_node_of_o_curve(oldc[2]),
					       oldc[2]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
	    	                               tmp_sten->leftst-1,
					       on_fr->rightst-1,
	    	                               on_fr->hs-1,on_fr->hse-1,
					       on_fr->t-1,on_fr->p-1,fr);
	    }
	    else
	    {
	    	states_at_distance_along_curve(p,
	    	                               Bond_at_node_of_o_curve(oldc[2]),
					       oldc[2]->curve,
					       POSITIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst+1,
					       tmp_sten->rightst+1,
					       tmp_sten->hs+1,tmp_sten->hse+1,
					       tmp_sten->t+1,tmp_sten->p+1,fr);
	    	for (i = 1; i <= tan_rad; ++i)
	    	{
	    	    ft_assign(on_fr->rightst[-i],tmp_sten->leftst[i],sizest);
		    on_fr->hs[-i] = tmp_sten->hs[i];
		    on_fr->hse[-i] = tmp_sten->hse[i];
		    on_fr->t[-i] = tmp_sten->t[i];
		    for (j = 0; j < dim; ++j)
			Coords(on_fr->p[-i])[j] = Coords(tmp_sten->p[i])[j];
		}
	    }
	    states_at_distance_along_curve(p,b,c,NEGATIVE_ORIENTATION,
					   uswssten->ds,tan_rad,on_fr->leftst-1,
					   tmp_sten->rightst-1,tmp_sten->hs-1,
					   tmp_sten->hse-1,tmp_sten->t-1,
			                   tmp_sten->p-1,fr);
	    break;
	default:
	    screen("ERROR in cc_set_uswssten0(), invalid orientation\n");
	    clean_up(ERROR);
	}
	for (i = -tan_rad; i <= tan_rad; ++i)
	{
	    WSSten	*nsten = nor_fr(uswssten)[i];
	    Locstate	sl = on_fr->leftst[i];
	    Locstate	sr = on_fr->rightst[i];
	    nsten->hs = hs;
	    for (j = 0; j < dim; ++j)
	    {
		nsten->coords[j] = Coords(on_fr->p[i])[j];
		nsten->lcrds[0][j] = nsten->rcrds[0][j] = nsten->coords[j];
		nsten->nor[j] = nor_vec(uswssten)[i][j];
	    }
	    nsten->dn = uswssten->dn;
	    nsten->dt = dt;
	    nsten->pjump = 0.0;
	    nsten->front = fr;
	    nsten->wave = wave;
	    copy_state(nsten->sl[0],sl);
	    copy_state(nsten->sr[0],sr);
	}

	set_USWSSten2d_off_front_states(uswssten,hs,fr,wave);
#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG == YES)
	    fprint_cc_set_uswssten(uswssten);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	DEBUG_LEAVE(cc_set_uswssten0)
}	/*end cc_set_uswssten0*/

LOCAL	void	cc_set_uswssten1(
	USWSSten2d		*uswssten,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	Front			*fr,
	Wave			*wave,
	double			dt)
{
	CURVE		*c = Curve_of_hs(hs);/*ASSUMES TWOD*/
	BOND		*b = Bond_of_hse(hse);/*ASSUMES TWOD*/
	Tan_stencil	*on_fr = tan_fr(uswssten)[0];
	int		i, j, dim = fr->rect_grid->dim;
	int		tan_rad = uswssten->tan_rad;
	size_t		sizest = fr->sizest;

	DEBUG_ENTER(cc_set_uswssten1)
	if (tmp_sten == NULL)
	    tmp_sten = alloc_tan_stencil(fr,tan_rad);

	uswssten->fr = fr;
	uswssten->wave = wave;
	uswssten->hs = hs;
	set_uswssten_geometry_and_center(uswssten,p,hse,hs,fr,dt);

	if (oldc[1]->orient == NEGATIVE_ORIENTATION)
	{
	    states_at_distance_along_curve(p,b,c,NEGATIVE_ORIENTATION,
					   uswssten->ds,tan_rad,
	    	                           on_fr->leftst-1,on_fr->rightst-1,
	    	                           on_fr->hs-1,on_fr->hse-1,
					   on_fr->t-1,on_fr->p-1,fr);
	    if (oldc[2]->orient == POSITIVE_ORIENTATION)
	    {
	    	states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[2]),
					       oldc[2]->curve,
					       POSITIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst+1,
					       on_fr->rightst+1,on_fr->hs+1,
					       on_fr->hse+1,on_fr->t+1,
					       on_fr->p+1,fr);
	    }
	    else
	    {
	    	states_at_distance_along_curve(p,
	    	                               Bond_at_node_of_o_curve(oldc[2]),
					       oldc[2]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst-1,
					       tmp_sten->rightst-1,
					       tmp_sten->hs-1,tmp_sten->hse-1,
					       tmp_sten->t-1, tmp_sten->p-1,fr);
	    	for (i = 1; i <= tan_rad; ++i)
	    	{
	    	    ft_assign(on_fr->rightst[i],tmp_sten->leftst[-i],sizest);
	    	    on_fr->hs[i] = tmp_sten->hs[-i];
	    	    on_fr->hse[i] = tmp_sten->hse[-i];
	    	    on_fr->t[i] = tmp_sten->t[-i];
	    	    for (j = 0; j < dim; ++j)
	    	    	Coords(on_fr->p[i])[j] = Coords(tmp_sten->p[-i])[j];
	    	}
	    }
	    states_at_distance_along_curve(p,b,c,POSITIVE_ORIENTATION,
					   uswssten->ds,tan_rad,on_fr->leftst+1,
					   tmp_sten->rightst+1,tmp_sten->hs+1,
					   tmp_sten->hse+1,tmp_sten->t+1,
					   tmp_sten->p+1,fr);
	}
	else
	{
	    states_at_distance_along_curve(p,b,c,POSITIVE_ORIENTATION,
					   uswssten->ds,tan_rad,
					   on_fr->leftst+1,on_fr->rightst+1,
					   on_fr->hs+1,on_fr->hse+1,
					   on_fr->t+1,on_fr->p+1,fr);
	    if (oldc[2]->orient == NEGATIVE_ORIENTATION)
	    {
	    	states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[2]),
					       oldc[2]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       on_fr->leftst-1,
					       tmp_sten->rightst-1,
					       on_fr->hs-1,on_fr->hse-1,
					       on_fr->t-1,on_fr->p-1,fr);
	    }
	    else
	    {
	    	states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[2]),
					       oldc[2]->curve,
					       POSITIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst+1,
					       tmp_sten->rightst+1,
					       tmp_sten->hs+1,tmp_sten->hse+1,
					       tmp_sten->t+1,tmp_sten->p+1,fr);
	    	for (i = 1; i <= tan_rad; ++i)
	    	{
	    	    ft_assign(on_fr->leftst[-i],tmp_sten->rightst[i],sizest);
	    	    on_fr->hs[-i] = tmp_sten->hs[i];
	    	    on_fr->hse[-i] = tmp_sten->hse[i];
	    	    on_fr->t[-i] = tmp_sten->t[i];
	    	    for (j = 0; j < dim; ++j)
	    	    	Coords(on_fr->p[-i])[j] = Coords(tmp_sten->p[i])[j];
	    	}
	    }
	    states_at_distance_along_curve(p,b,c,NEGATIVE_ORIENTATION,
					   uswssten->ds,tan_rad,
					   tmp_sten->leftst-1,on_fr->rightst-1,
					   tmp_sten->hs-1,tmp_sten->hse-1,
					   tmp_sten->t-1, tmp_sten->p-1,fr);
	}
	for (i = -tan_rad; i <= tan_rad; ++i)
	{
	    WSSten	*nsten = nor_fr(uswssten)[i];
	    Locstate	sl = on_fr->leftst[i];
	    Locstate	sr = on_fr->rightst[i];
	    nsten->hs = hs;
	    for (j = 0; j < dim; ++j)
	    {
		nsten->coords[j] = Coords(on_fr->p[i])[j];
		nsten->lcrds[0][j] = nsten->rcrds[0][j] = nsten->coords[j];
		nsten->nor[j] = nor_vec(uswssten)[i][j];
	    }
	    nsten->dn = uswssten->dn;
	    nsten->dt = dt;
	    nsten->pjump = 0.0;
	    nsten->front = fr;
	    nsten->wave = wave;
	    copy_state(nsten->sl[0],sl);
	    copy_state(nsten->sr[0],sr);
	}

	set_USWSSten2d_off_front_states(uswssten,hs,fr,wave);
#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG == YES)
	    fprint_cc_set_uswssten(uswssten);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	DEBUG_LEAVE(cc_set_uswssten1)
}	/*end cc_set_uswssten1*/

LOCAL	void	cc_set_uswssten2(
	USWSSten2d		*uswssten,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	Front			*fr,
	Wave			*wave,
	double			dt)
{
	CURVE       *c = Curve_of_hs(hs);/*ASSUMES TWOD*/
	BOND        *b = Bond_of_hse(hse);/*ASSUMES TWOD*/
	Tan_stencil *on_fr = tan_fr(uswssten)[0];
	int         i, j, dim = fr->rect_grid->dim;
	int	    tan_rad = uswssten->tan_rad;
	size_t      sizest = fr->sizest;

	DEBUG_ENTER(cc_set_uswssten2)

	 /* tmp_sten used for temporary holder and as data sink. */
	if (tmp_sten == NULL)
	    tmp_sten = alloc_tan_stencil(fr,tan_rad);

	/*
	 * NOTE: Stencil geometry determined by nor_vec(uswssten)[0] and
	 * corresponding tangent which in this case have been set to
	 * cc_normal, cc_tangent and return different uni_arrays than obtain from
	 * usual adjacent-points secant algorithm.  Also sets states at
	 * stencil center.
	 */
	uswssten->fr = fr;
	uswssten->wave = wave;
	uswssten->hs = hs;
	set_uswssten_geometry_and_center(uswssten,p,hse,hs,fr,dt);

	/* Curve 2 is the PRINCIPAL curve ... the one being propagated. */

	if (oldc[2]->orient == NEGATIVE_ORIENTATION)
	{
	    /*
	     * On front data for PRINCIPAL curve.
	     *
	     * NOTE: states_at_distance_along_curve() returns, as
	     * penultimate arg, POINT from which states are retrieved.
	     */
	    states_at_distance_along_curve(p,b,c,NEGATIVE_ORIENTATION,
					   uswssten->ds,tan_rad,on_fr->leftst-1,
					   on_fr->rightst-1,on_fr->hs-1,
					   on_fr->hse-1,on_fr->t-1,on_fr->p-1,
					   fr);
	    if (oldc[0]->orient == POSITIVE_ORIENTATION)
	    {
	        /*
	         * Absolute orient for curve 0 matches that for
	         * PRINCIPAL curve so transfer on_fr data directly.
	         */
	        states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[0]),
					       oldc[0]->curve,
					       POSITIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst+1,
					       on_fr->rightst+1,
					       on_fr->hs+1,on_fr->hse+1,
					       on_fr->t+1,on_fr->p+1,fr);
	    }
	    else
	    {
	        /*
	         * Absolute orient for curve 0 opposite that for
	         * PRINCIPAL curve so transfer on_fr data to tmp_sten,
	         * then reverse indices when copying to on_fr.
	         */
	        states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[0]),
					       oldc[0]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst-1,
					       tmp_sten->rightst-1,
					       tmp_sten->hs-1,tmp_sten->hse-1,
					       tmp_sten->t-1,tmp_sten->p-1,fr);
	        for (i = 1; i <= tan_rad; ++i)
	        {
	            ft_assign(on_fr->rightst[i],tmp_sten->leftst[-i],sizest);
	            on_fr->hs[i] = tmp_sten->hs[-i];
	            on_fr->hse[i] = tmp_sten->hse[-i];
	            on_fr->t[i] = tmp_sten->t[-i];
	            for (j = 0; j < dim; ++j)
	                Coords(on_fr->p[i])[j] =
	                    Coords(tmp_sten->p[-i])[j];
	        }
	    }
	    /* Same for curve 1 */
	    if (oldc[1]->orient == POSITIVE_ORIENTATION)
	    {
	        states_at_distance_along_curve(p,
	                                       Bond_at_node_of_o_curve(oldc[1]),
					       oldc[1]->curve,
					       POSITIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       on_fr->leftst+1,
					       tmp_sten->rightst+1,
					       tmp_sten->hs+1,tmp_sten->hse+1,
					       tmp_sten->t+1,tmp_sten->p+1,fr);
	    }
	    else
	    {
	        states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[1]),
					       oldc[1]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst-1,
					       tmp_sten->rightst-1,
					       tmp_sten->hs-1,tmp_sten->hse-1,
					       tmp_sten->t-1,tmp_sten->p-1,fr);
	        for (i = 1; i <= tan_rad; ++i)
	        {
	            ft_assign(on_fr->leftst[i],tmp_sten->rightst[-i],sizest);
	            tmp_sten->hs[i] = tmp_sten->hs[-i];
	            tmp_sten->hse[i] = tmp_sten->hse[-i];
	            tmp_sten->t[i] = tmp_sten->t[-i];
	            for (j = 0; j < dim; ++j)
		    {
			/* Save coords from curve 0 */
			Coords(tmp_sten->p[i])[j] = Coords(on_fr->p[-i])[j];
			/* Copy coords from curve 1 */
			Coords(on_fr->p[i])[j] = Coords(tmp_sten->p[-i])[j];
		    }
	        }
	    }
	    /*
	     * Finish setting section of on_fr stencil that lies in the
	     * "ignored" sector between curves 0, 1.
	     */
	    for (i = 1; i <= tan_rad; ++i)
	    {
	        on_fr->hs[i] = Hyper_surf(c);
	        on_fr->hse[i] =
		    Hyper_surf_element(Bond_at_node_of_o_curve(oldc[2]));
	        on_fr->t[i] = 1.0;
	        /*
	         * Set on_fr->p[1] to average of positions from which
	         * left/right states were extracted.
	         */
	        for (j = 0; j < dim; ++j)
	        {
	            Coords(on_fr->p[i])[j] = 0.5*(Coords(on_fr->p[i])[j] +
	                                          Coords(tmp_sten->p[i])[j]);
	        }
	    }
	    /*
	     * TODO: write debug code to print angles between
	     * uni_arrays returned by cc_tangent() and
	     * (on_fr->p[0] - on_fr->p[-1]), (on_fr->p[1] - on_fr->p[0]).
	     */
	}
	else
	{
	    /* EXACTLY same as prev. case modulo PRINCIPAL curve orient. */
	    states_at_distance_along_curve(p,b,c,POSITIVE_ORIENTATION,
					   uswssten->ds,tan_rad,
					   on_fr->leftst+1,on_fr->rightst+1,
					   on_fr->hs+1,on_fr->hse+1,
					   on_fr->t+1,on_fr->p+1,fr);
	    if (oldc[0]->orient == NEGATIVE_ORIENTATION)
	    {
	        states_at_distance_along_curve(p,
	                                       Bond_at_node_of_o_curve(oldc[0]),
					       oldc[0]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       on_fr->leftst-1,
					       tmp_sten->rightst-1,
					       on_fr->hs-1,on_fr->hse-1,
					       on_fr->t-1,on_fr->p-1,fr);
	    }
	    else
	    {
	        states_at_distance_along_curve(p,
	                                       Bond_at_node_of_o_curve(oldc[0]),
					       oldc[0]->curve,
					       POSITIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst+1,
					       tmp_sten->rightst+1,
					       tmp_sten->hs+1,tmp_sten->hse+1,
					       tmp_sten->t+1,tmp_sten->p+1,fr);
	        for (i = 1; i <= tan_rad; ++i)
	        {
	            ft_assign(on_fr->leftst[-i],tmp_sten->rightst[i],sizest);
	            on_fr->hs[-i] = tmp_sten->hs[i];
	            on_fr->hse[-i] = tmp_sten->hse[i];
	            on_fr->t[-i] = tmp_sten->t[i];
	            for (j = 0; j < dim; ++j)
	                Coords(on_fr->p[-i])[j] = Coords(tmp_sten->p[i])[j];
	        }
	    }
	    if (oldc[1]->orient == NEGATIVE_ORIENTATION)
	    {
	        states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[1]),
					       oldc[1]->curve,
					       NEGATIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst-1,
					       on_fr->rightst-1,
					       tmp_sten->hs-1,tmp_sten->hse-1,
					       tmp_sten->t-1,tmp_sten->p-1,fr);
	    }
	    else
	    {
	        states_at_distance_along_curve(p,
					       Bond_at_node_of_o_curve(oldc[1]),
					       oldc[1]->curve,
					       POSITIVE_ORIENTATION,
					       uswssten->ds,tan_rad,
					       tmp_sten->leftst+1,
					       tmp_sten->rightst+1,
					       tmp_sten->hs+1,tmp_sten->hse+1,
					       tmp_sten->t+1,tmp_sten->p+1,fr);
	        for (i = 1; i <= tan_rad; ++i)
	        {
	            ft_assign(on_fr->rightst[-i],tmp_sten->leftst[i],sizest);
	            tmp_sten->hs[-i] = tmp_sten->hs[i];
	            tmp_sten->hse[-i] = tmp_sten->hse[i];
	            tmp_sten->t[-i] = tmp_sten->t[i];
	            for (j = 0; j < dim; ++j)
		    {
			/* Save coords from curve 0 */
			Coords(tmp_sten->p[-i])[j] = Coords(on_fr->p[-i])[j];
			/* Copy coords from curve 1 */
	                Coords(on_fr->p[-i])[j] = Coords(tmp_sten->p[i])[j];
		    }
	        }
	    }
	    for (i = 1; i <= tan_rad; ++i)
	    {
	        on_fr->hs[-i] = Hyper_surf(c);
	        on_fr->hse[-i] =
		    Hyper_surf_element(Bond_at_node_of_o_curve(oldc[2]));
	        on_fr->t[-i] = 0.0;
	        for (j = 0; j < dim; ++j)
	        {
	            Coords(on_fr->p[-i])[j] = 0.5*(Coords(on_fr->p[-i])[j] +
	                                           Coords(tmp_sten->p[-i])[j]);
	        }
	    }
	}
	for (i = -tan_rad; i <= tan_rad; ++i)
	{
	    WSSten	*nsten = nor_fr(uswssten)[i];
	    Locstate	sl = on_fr->leftst[i];
	    Locstate	sr = on_fr->rightst[i];
	    nsten->hs = hs;
	    for (j = 0; j < dim; ++j)
	    {
		nsten->coords[j] = Coords(on_fr->p[i])[j];
		nsten->lcrds[0][j] = nsten->rcrds[0][j] = nsten->coords[j];
		nsten->nor[j] = nor_vec(uswssten)[i][j];
	    }
	    nsten->dn = uswssten->dn;
	    nsten->dt = dt;
	    nsten->pjump = 0.0;
	    nsten->front = fr;
	    nsten->wave = wave;
	    copy_state(nsten->sl[0],sl);
	    copy_state(nsten->sr[0],sr);
	}

	set_USWSSten2d_off_front_states(uswssten,hs,fr,wave);
#if defined(DEBUG_NODE_PROPAGATE)
	if (DEBUG == YES)
	    fprint_cc_set_uswssten(uswssten);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	DEBUG_LEAVE(cc_set_uswssten2)
}    /*end cc_set_uswssten2*/

/*ARGSUSED*/
LOCAL	void	cc_normal(
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	double			*nor,
	Front			*fr)
{
	INTERFACE *intfc = hs->interface;
	int i, dim = fr->rect_grid->dim;
	NORMAL_FUNCTION Sav_nor_func;
	TANGENT_FUNCTION Sav_tan_func;

	DEBUG_ENTER(cc_normal)

	Sav_nor_func = interface_normal_function(intfc);
	interface_normal_function(intfc) = _fr_normal;
	Sav_tan_func = interface_tangent_function(intfc);
	interface_tangent_function(intfc) = _fr_tangent;
	normal(p,hse,hs,nor,fr);
	interface_normal_function(intfc) = Sav_nor_func;
	interface_tangent_function(intfc) = Sav_tan_func;

	if (scalar_product(cc_nor,nor,dim) < 0.0)
	{
	    for (i = 0; i < dim; ++i)
	    	nor[i] = -cc_nor[i];
	}
	else
	{
	    for (i = 0; i < dim; ++i)
	    	nor[i] = cc_nor[i];
	}
	
	DEBUG_LEAVE(cc_normal)
}		/*end cc_normal*/

EXPORT	void	set_cc_tangent_function(
	TANGENT_FUNCTION *tf)
{
	static const char *tname = "cc_tangent";
	tf->_tangent = cc_tangent;
	tf->_tangent_name = tname;
}		/*end set_cc_tangent_function*/

EXPORT	void	set_cc_normal_function(
	NORMAL_FUNCTION *nf)
{
	static const char *tname = "cc_normal";
	nf->_normal = cc_normal;
	nf->_normal_name = tname;
}		/*end set_cc_normal_function*/

/*ARGSUSED*/
LOCAL	void	cc_tangent(
	POINT		*p,
	BOND		*b,
	CURVE		*c,
	double		*t,
	Front		*fr)
{
	INTERFACE *intfc = c->interface;
	int i, dim = fr->rect_grid->dim;
	NORMAL_FUNCTION Sav_nor_func;
	TANGENT_FUNCTION Sav_tan_func;

	DEBUG_ENTER(cc_tangent)

	Sav_nor_func = interface_normal_function(intfc);
	interface_normal_function(intfc) = _fr_normal;
	Sav_tan_func = interface_tangent_function(intfc);
	interface_tangent_function(intfc) = _fr_tangent;
	tangent(p,b,c,t,fr);
	interface_normal_function(intfc) = Sav_nor_func;
	interface_tangent_function(intfc) = Sav_tan_func;

	if (scalar_product(cc_tan,t,dim) < 0.0)
	{
	    for (i = 0; i < dim; ++i)
	    	t[i] = -cc_tan[i];
	}
	else
	{
	    for (i = 0; i < dim; ++i)
	    	t[i] = cc_tan[i];
	}
	DEBUG_LEAVE(cc_tangent)
}	/*end cc_tangent*/

#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void fprint_cc_set_uswssten(
	USWSSten2d	*uswssten)
{
	int	i,j;
	int	tan_rad = uswssten->tan_rad, dim = uswssten->fr->interf->dim;

	static const char *fmt[3] = {"%-22s%s","%-22s%s","%-33s%s"};
	static const char *vname[3] = {"X", "Y", "Z"};

	(void) printf("\tCoords from Tan_stencils.\n");
	(void) printf("%-9s","INDEX");
	(void) printf(fmt[dim-1],"ON FRONT POINT","	   ");
	(void) printf(fmt[dim-1],"LEFT OFF FRONT POINT","	 ");
	(void) printf(fmt[dim-1],"RIGHT OFF FRONT POINT","\n");
	(void) printf("%-9s","I");
	for (j = 0; j < 3; ++j)
	{
	    for (i = 0; i < dim; ++i)
		(void) printf("%-11s",vname[i]);
	    (void) printf("%s",(j < 2) ? "	  " : "\n");
	}
	for (j = -tan_rad; j <= tan_rad ; ++j)
	{
	    (void) printf("%-9d",j);
	    for (i = 0; i < dim; ++i)
		(void) printf("%- 11g",Coords(tan_fr(uswssten)[0]->p[j])[i]);
	    (void) printf("	   ");
	    for (i = 0; i < dim; ++i)
		(void) printf("%- 11g",Coords(tan_fr(uswssten)[-1]->p[j])[i]);
	    (void) printf("	   ");
	    for (i = 0; i < dim; ++i)
		(void) printf("%- 11g",Coords(tan_fr(uswssten)[1]->p[j])[i]);
	    (void) printf("\n");
	}
	g_PrintUSWSSten2d(uswssten);
}	/*end fprint_cc_set_uswssten*/
#endif /* defined(DEBUG_NODE_PROPAGATE) */
#endif /* defined(TWOD) */
