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
*				gnodesub.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       Contains physics dependent subroutines used in the propagation
*	of nodes.
*
*	The subroutines are:
*
*		find_curve_with_status()
*		identify_curves_with_status()
*		find_adjacent_curves()
*		g_find_i_to_prop_dir()
*		modify_regular_reflection_node()
*		modify_curves_at_node()
*		modify_transmission_node()
*		modify_CC_node()
*		corner_node_propagate()
*		propagate_curve_near_node()
*
*	The topology of the new interface might change during the propagation.
*	This may necessitate adding new bonds, curves, and nodes on the new
*	interface.
*/


#if defined(TWOD)

#include <gdecs/gdecs.h>


	/* LOCAL Function Declarations */
LOCAL	POINT	*eliminate_short_bonds_at_node(O_CURVE*,O_CURVE*,double,
					       Front*,Wave*,double*);
LOCAL	boolean	check_for_consistency_at_node(O_CURVE*,double*,double*,
					      POINT*,Front*);
LOCAL	double	sonic_radius_at_node(NODE*,Locstate,Locstate,double,int);
LOCAL	void	coords_along_ray(double*,double,double*,double*,double*,double,int);
LOCAL	void	ensure_states_set_at_opp_end_of_bond(O_CURVE*,O_CURVE*,POINT*,
						     Front*,Wave*,double);
LOCAL	void	extend_states_out_from_node(O_CURVE*,O_CURVE*,Front*,Wave*,
					    double);


EXPORT void find_curve_with_status(
	NODE		*n,
	CURVE		**c,
	ORIENTATION	*orient,
	int		status)
{
	CURVE		**pc;

	*c = NULL;
	if (n == NULL)
		return;

	if (n->in_curves != NULL)
		for (pc = n->in_curves; *pc; pc++)
			if (end_status(*pc) == status)
			{
				*c = *pc;
				*orient = NEGATIVE_ORIENTATION;
				return;
			}

	if (n->out_curves != NULL)
		for (pc = n->out_curves; *pc; pc++)
			if (start_status(*pc) == status)
			{
				*c = *pc;
				*orient = POSITIVE_ORIENTATION;
				return;
			}
}		/*end find_curve_with_status*/


/*
*		identify_curves_with_status():
*
*	Finds up to two curves at the given node with status at node
*	status.
*/

EXPORT	void identify_curves_with_status(
	NODE		*node,
	O_CURVE		*c0,
	O_CURVE		*c1,
	int		status)
{
	CURVE		**c;

	c0->curve = c1->curve = NULL;
	for (c = node->in_curves; c && *c; c++) 
	{
		if (end_status(*c) == status) 
		{
			if (c0->curve) 
			{
				c1->curve = *c;
				c1->orient = NEGATIVE_ORIENTATION;
			}
			else 
			{
				c0->curve = *c;
				c0->orient = NEGATIVE_ORIENTATION;
			}
		}
		if (c0->curve && c1->curve) break;
	}
	if (!c0->curve || !c1->curve) 
	{
		for (c = node->out_curves; c && *c; c++) 
		{
			if (start_status(*c) == status) 
			{
				if (c0->curve) 
				{
					c1->curve = *c;
					c1->orient = POSITIVE_ORIENTATION;
				}
				else 
				{
					c0->curve = *c;
					c0->orient = POSITIVE_ORIENTATION;
				}
			}
			if (c0->curve && c1->curve) break;
		}
	}

	if (!c0->curve) return;

}		/*end identify_curves_with_status*/



/*
*		find_adjacent_curves():
*/

EXPORT void find_adjacent_curves(
	O_CURVE		*cinc,
	ANGLE_DIRECTION	*angle_dir,
	O_CURVE		*ca,
	O_CURVE		*cb,
	COMPONENT	*ahead_comp)
{
	if (is_forward_wave(wave_type(cinc->curve)))
	{
		*ahead_comp = positive_component(cinc->curve);
		*angle_dir = (cinc->orient == POSITIVE_ORIENTATION) ?
			CLOCKWISE : COUNTER_CLOCK;
	}
	else
	{
		*ahead_comp = negative_component(cinc->curve);
		*angle_dir = (cinc->orient == POSITIVE_ORIENTATION) ?
			COUNTER_CLOCK : CLOCKWISE;
	}

	ca->curve = adjacent_curve(cinc->curve,cinc->orient,
		*angle_dir,&ca->orient);
	cb->curve = adjacent_curve(cinc->curve,cinc->orient,
		Opposite_ang_dir(*angle_dir),&cb->orient);
}		/*end find_adjacent_curves*/


/*
*			g_find_i_to_prop_dir():
*
*	This function is a physics dependent driver for f_find_i_to_prop_dir()
*	which finds the angle direction from the incident toward the
*	directon of propagation of the curve.  For dt == 0, the
*	point_propagate() used to find the node displacement is zero.
*	In this case we need to find a small virtual dt for the
*	propagation so we get a valid displacement uni_array.
*/

EXPORT ANGLE_DIRECTION g_find_i_to_prop_dir(
	Front		*fr,
	POINTER		p2wave,
	NODE		*oldn,
	CURVE		*oldc,
	ORIENTATION	c_orient,
	double		dt,
	COMPONENT	*ahead_comp,
	POINT		*newp,
	double		*V)
{
	Wave		*wave = (Wave *)p2wave;
	double		virtual_dt, max_dt;
	double		coords[MAXD];
	double		CFL;

	CFL = Time_step_factor(fr);
	max_dt = (*wave->max_hyp_time_step)(wave,coords);
#if !defined(_AIX)
	if (finite(max_dt) && finite(-max_dt))
#endif /* !defined(_AIX) */
	    max_dt *= 0.1*CFL;/*TOLERANCE*/
	virtual_dt = max(dt,max_dt);

	/* this is needed for initial time step && restart */
	if (virtual_dt > 100000.0) /*TOLERANCE*/
	    virtual_dt = dt;

	if (debugging("i_to_prop"))
	    (void) printf("dt = %g, virtual_dt = %g\n",dt,virtual_dt);

	return f_find_i_to_prop_dir(fr,(POINTER)wave,oldn,oldc,
				    c_orient,virtual_dt,ahead_comp,newp,V);
}		/*end g_find_i_to_prop_dir*/


/*
*		       	modify_curves_at_node():
*
*/

EXPORT	int modify_curves_at_node(
	POINT		*pc,
	BOND		**newb,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	int		max_n_cur,
	boolean		*correct_angle_at_node,
	Front		*fr,
	Wave		*wave,
	double		dt,
	RP_DATA		*RP,
	NODE_FLAG	flag)
{
	Locstate	st_l, st_r;
	int		i, dim = fr->rect_grid->dim;

	debug_print("modify","Entered modify_curves_at_node()\n");

	for (i = 0; i < dim; i++)
	    Coords(newn->posn)[i] = Coords(pc)[i];

	for (i = 0; i < max_n_cur; i++) 
	{
	    if (!newc[i] || !newc[i]->curve) continue;
	    if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[i]->orient))
	    {
	    	st_l = RP->state[i];
	    	st_r = RP->state[(i+1) % max_n_cur];
	    }
	    else 
	    {
	    	st_l = RP->state[(i+1)% max_n_cur];
	    	st_r = RP->state[i];
	    }
	    if (!propagate_curve_near_node(Node_of_o_curve(oldc[i]),
			                      Node_of_o_curve(newc[i]),
					      oldc[i],newc[i],newb[i],
			                      st_l,st_r,
					      correct_angle_at_node[i],
			                      RP->ang[i],fr,wave,dt,flag)) 
	    {
	    	(void) printf("WARNING in modify_curves_at_node(), "
	    	              "propagate_curve_near_node() failed\n");
	    	debug_print("modify","Left modify_curves_at_node()\n");
	    	return NO;
	    }
	}
	debug_print("modify","Left modify_curves_at_node()\n");
	return YES;
}		/*end modify_curves_at_node*/


/*
*
*			propagate_curve_near_node():
*
*	This routine is called for the curves around the node. It  propagates
*	the point adjacent to the old node and then inserts a new bond into
*	the new curve to satisfy a prescribed angle.
*/

EXPORT int propagate_curve_near_node(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		*oldc,
	O_CURVE		*newc,
	BOND		*newb,
	Locstate	st_l,
	Locstate	st_r,
	boolean		correct_angle_at_node,
	double		ang,
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE_FLAG	flag)
{
	INTERFACE	*intfc;
	NODE		*oppn = Opp_node_of_o_curve(newc);
	POINT		*p0, *p00, *pn1;
	double		coords[MAXD];
	double		d[MAXD];
	double		d0[MAXD];
	double		t[MAXD];   /* tangent to old bond at oldn */
	double		dir[MAXD];
	double		len,rel_len;
	double		*h = fr->rect_grid->h;
	double		v[MAXD];
	double		adjust_len;
	double		sonic_rad;
	int		status;
	boolean		sav_intrp;
	int		i, dim = fr->rect_grid->dim;

	debug_print("propagate","Entered propagate_curve_near_node()\n");

	if (debugging("propagate")) 
	{
		(void) printf("\t\tOLD CURVE:\n");
		print_o_curve(oldc);
		(void) printf("\t\tNEW CURVE before cut paste:\n");
		print_o_curve(newc);
		(void) printf("\t\tOLD NODE:\n");
		print_node(oldn);
		(void) printf("\t\tNEW NODE:\n");
		print_node(newn);
	}
	interpolate_intfc_states(newn->interface) = YES;
	p00 = oldn->posn;
	p0 = newn->posn;
	dir[0] = cos(ang);
	dir[1] = sin(ang);

	/* Cut new curve at new bond and assign states at new node */

	cut_curve(p0,newb,newc->curve,newc->orient,fr,st_l,st_r);

	if (debugging("propagate")) 
	{
		(void) printf("Bond_at_node_of_o_curve(newc)\n");
		print_bond(Bond_at_node_of_o_curve(newc));
	}

	/* Find direction of node propagation */

	for (i = 0; i < dim; i++)
		d[i] = Coords(p0)[i] - Coords(p00)[i];
	len = mag_vector(d,dim);
	rel_len = scaled_hypot(d,h,dim);

	/* Find direction vector for old curve */

	find_tangent_to_curve(p00,Bond_at_node_of_o_curve(oldc),oldc->curve,
		oldc->orient,t,fr);

	if (debugging("propagate") || debugging("adjust"))
	{
	   print_boolean("dont_correct_angles_at_node(flag) = ",
			 dont_correct_angles_at_node(flag),"\n");
	   print_boolean("correct_angle_at_node = ",correct_angle_at_node,"\n");
	   print_angle("ang = ",ang,"\n");
	   (void) printf("\n");
	   print_general_vector("dir = ",dir,dim,"");
	   print_general_vector("\tt = ",t,dim,"\n");
	   (void) printf("scalar_product(t,dir,dim) = %g\n",
	                 scalar_product(t,dir,dim));
	}

	/* Propagate old node and insert in new curve if necessary */

	if (rel_len*scalar_product(d,t,dim) < -2.0*len)/*TOLERANCE*/
	{
	    coords[0] = 0.0;	coords[1] = 1.0;
	    pn1 = Point(coords);
	    point_propagate(fr,(POINTER)wave,oldn->posn,pn1,
	    	Bond_at_node_of_o_curve(oldc),oldc->curve,dt,v);
	    if (oldc->orient != newc->orient)
	    	reverse_states_at_point(pn1,fr);

	    for (i = 0; i < dim; i++)
	    	d0[i] = Coords(pn1)[i] - Coords(p0)[i];
	    if ((scalar_product(d0,dir,dim) > 0.0) &&
	        (scaled_hypot(d0,h,dim) > 0.001) /*TOLERANCE*/)
	    {
	        if ((newc->curve->first == newc->curve->last) &&
	            (propagation_status(oppn) == UNPROPAGATED_NODE))
	        {
	            (void) printf("WARNING in propagate_curve_near_node(), ");
	            (void) printf("direction reversal, ");
	            (void) printf("short curve with oppn unpropagated\n");
	            return NO;
	        }
	        intfc = newc->curve->interface;
	        sav_intrp = interpolate_intfc_states(intfc);
	        interpolate_intfc_states(intfc) = NO;
	        insert_point_adjacent_to_node(pn1,newc->curve,newc->orient);
		interpolate_intfc_states(intfc) = sav_intrp;
	        if (debugging("propagate"))
	        {
	            BOND	*b = Bond_at_node_of_o_curve(newc);
	            CURVE	*c = newc->curve;

	            (void) printf("Inserted point from node ");
	            (void) printf("propagation\n");
	            (void) printf("pn1 = (%g, %g)\n",
	                   Coords(pn1)[0],Coords(pn1)[1]);
	            verbose_print_state("Left state pn1",
	                  left_state_at_point_on_curve(pn1,b,c));
	            verbose_print_state("Right state pn1",
	                  right_state_at_point_on_curve(pn1,b,c));
	        }
	    }
	}

	status = YES;
	if ((dont_correct_angles_at_node(flag) != YES) &&
	    (correct_angle_at_node == YES) &&
	    (scalar_product(t,dir,dim) >= 0.0))
	{
	    if (adjust_angle_len(newn,oldc,newc,st_l,st_r,dir,dt,&adjust_len,
				 &sonic_rad,fr,wave) == ADJUST_ANGLE)
	    {
	        if (debugging("adjust"))
	        {
	    	    (void) printf("len = %g, adjust_len = %g, sonic_rad = %g\n",
				  len,adjust_len,sonic_rad);
	        }
	        if (adjust_len*scaled_hypot(dir,h,dim) > MIN_SC_SEP(fr->interf))
	        {				/*TOLERANCE*/
	    	    status = adjust_angle_at_node(newn,oldc,newc,st_l,st_r,
	    			   		  dir,dt,adjust_len,fr,wave);
	        }
	        else if (debugging("propagate") || debugging("adjust"))
	        {
	    	    (void) printf("WARNING in propagate_curve_near_node(), "
	    	                  "did not adjust angle because "
				  "adjust_len is too short\n");
	        }
	    }
	}
	else if (debugging("propagate") || debugging("adjust"))
	{
	    (void) printf("WARNING in propagate_curve_near_node(), ");
	    (void) printf("did not adjust angle\n");
	    if (debugging("propagate") || debugging("adjust"))
	    {
	        print_boolean("dont_correct_angles_at_node(flag) = ",
			      dont_correct_angles_at_node(flag),"\n");
	        print_boolean("correct_angle_at_node = ",
			      correct_angle_at_node,"\n");
	        (void) printf("scalar_product(t,dir,dim) = %g\n",
			      scalar_product(t,dir,dim));
	    }
	}

	if (debugging("propagate"))
	{
	    (void) printf("final curve after propagating:\n");
	    print_o_curve(newc);
	}

	debug_print("propagate","Left propagate_curve_near_node\n");
	return status;
}		/*end propagate_curve_near_node*/

/*ARGSUSED*/
EXPORT ADJUST_ANGLE_VALUE g_adjust_angle_len(
	NODE		*newn,
	O_CURVE		*oldc,
	O_CURVE		*newc,
	Locstate	st_l,
	Locstate	st_r,
	double		*dir,
	double		dt,
	double		*adjust_len,
	double		*sonic_rad,
	Front		*fr,
	Wave		*wave)
{
	ADJUST_ANGLE_VALUE adjust_flag = ADJUST_ANGLE;
	POINT		*p0, *p00;
	double		len;
	double		*h = fr->rect_grid->h;
	int		dim = fr->rect_grid->dim;

	debug_print("adjust","Entered g_adjust_angle_len()\n");
	*sonic_rad = sonic_radius_at_node(newn,st_l,st_r,dt,dim);
	if (adjust_len(newn) > 0.0)
	{
	    if (debugging("adjust"))
		(void) printf("adjust length set from node value\n");
	    *adjust_len = adjust_len(newn);
	}
	else
	{
	    p0 = Node_of_o_curve(oldc)->posn;
	    p00 = Node_of_o_curve(newc)->posn;
	    len = separation(p0,p00,dim);
	    if (len*scaled_hypot(dir,h,dim) > MIN_SC_SEP(fr->interf)) /*TOLERANCE*/
	    {
	        if (debugging("adjust"))
		{
		    (void) printf("Angle will not be adjusted since "
				  "old and new node positions are too close\n");
		}
		adjust_flag = DONT_ADJUST_ANGLE;
	    }
	    *adjust_len = dt*max(sound_speed(st_l),sound_speed(st_r));
	    if (*sonic_rad < *adjust_len)
		*adjust_len = *sonic_rad;
	}

	debug_print("adjust","Left g_adjust_angle_len()\n");
	return adjust_flag;
}		/*end g_adjust_angle_len*/


/*
*			adjust_angle_at_node():
*
*	This function adjusts the point(s) near the node so that the
*	curve satisfies a prescribed angle.  The user must supply the
*	direction dir, and the right and left states for a new point
*	adjacent to the node.  len is compared with the sonic radius at
*	the node, and the smaller is used.
*
*	Note: failure of this routine is not always correctly handled.  It
*	can lead to an infinite loop because the resulting failed node
*	propagation is always given status PSEUDOCROSS_NODE_NODE.
*/

#define LEAVE_ADJUST_ANGLE_AT_NODE					\
{									\
	if (tmp_point != NULL)						\
	{								\
	    	(void) delete_point_adjacent_to_node(fr,newc->curve,	\
				Opposite_orient(newc->orient));		\
	}								\
	fuh->_bi_interpolate_intfc_states = sav_interpolator;		\
	fuh->_tri_interpolate_intfc_states = sav_tri_interpolator;	\
	interpolate_intfc_states(intfc) = sav_intrp;			\
	debug_print("adjust","Left adjust_angle_at_node()\n");		\
}

EXPORT int adjust_angle_at_node(
	NODE		*newn,
	O_CURVE		*oldc,
	O_CURVE		*newc,
	Locstate	st_l,
	Locstate	st_r,
	double		*dir,
	double		dt,
	double		adjust_len,
	Front		*fr,
	Wave		*wave)
{
	BOND		*b;
	INTERFACE 	*intfc = newc->curve->interface;
	F_USER_INTERFACE *fuh = &f_user_interface(intfc);
	POINT		*tmp_point = NULL;
	POINT		*p0 = newn->posn, *pnew;
	double		new_coords[MAXD];
	double		*h = fr->rect_grid->h;
	double		dist;
	double		d[MAXD];
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	int		dim = fr->rect_grid->dim;
	boolean		(*sav_tri_interpolator)(double,double,double,double*,
						Locstate,double*,Locstate,
						double*,Locstate,RECT_GRID*,
						Locstate);
	void		(*sav_interpolator)(double,double,double*,Locstate,
					    double*,Locstate,RECT_GRID*,
					    Locstate);

	debug_print("adjust","Entered adjust_angle_at_node()\n");

	if (debugging("adjust"))
	    (void) printf("adjust_len = %g\n",adjust_len);

	sav_tri_interpolator = fuh->_tri_interpolate_intfc_states;
	sav_interpolator = fuh->_bi_interpolate_intfc_states;
	dist = adjust_len +
	    2.0*MIN_SC_SEP(fr->interf)/scaled_hypot(dir,h,dim); /*TOLERANCE*/
	tmp_point = eliminate_short_bonds_at_node(oldc,newc,dt,fr,wave,&dist);

	if (check_for_consistency_at_node(newc,d,dir,tmp_point,fr) ==
							FUNCTION_FAILED)
	{
	    double ang = angle(dir[0],dir[1]);

	    (void) printf("WARNING in adjust_angle_at_node(), "
	                  "Incompatible angles\n");
	    print_angle("ang =",ang,"\n");
	    print_general_vector("dir = ",dir,dim,"\n");
	    (void) printf("Old Curve\n");
	    print_o_curve(oldc);
	    (void) printf("New Curve\n");
	    print_o_curve(newc);
	    LEAVE_ADJUST_ANGLE_AT_NODE
	    return NO;
	}
	b = Bond_at_node_of_o_curve(newc);
	if ((bond_length(b) < dist) ||
	    (scaled_bond_length(b,h,dim) < 0.25)) /*TOLERANCE*/
	{
	    if (debugging("adjust"))
	    {
	        (void) printf("WARNING in adjust_angle_at_node(), "
	                      "unable to adjust due to bond length tests\n");
	        (void) printf("\tscaled_bond_length(b) = %g\tdist = %g\n",
	    	              scaled_bond_length(b,h,dim),dist);
		print_bond(b);
		print_o_curve(oldc);
	    }
	    LEAVE_ADJUST_ANGLE_AT_NODE
	    return YES;
	}

		/* Compute coords of new point */

	coords_along_ray(new_coords,dist,d,Coords(p0),dir,adjust_len,dim);

	pnew = Point(new_coords);

	ensure_states_set_at_opp_end_of_bond(oldc,newc,tmp_point,fr,wave,dt);

	if (debugging("adjust")) 
	{
	    (void) printf("Inserting point newp = (%g, %g)\n",
	                  Coords(pnew)[0],Coords(pnew)[1]);
	    (void) printf("scaled_separation(pnew,p0,h,dim) = %g\n",
		          scaled_separation(pnew,p0,h,dim));
	    verbose_print_state("st_l",st_l);
	    verbose_print_state("st_r",st_r);
	}

	        /* Insert new point to satisfy prescribed angle */

	interpolate_intfc_states(intfc) = YES;
	if (is_scalar_wave(wave_type(newc->curve)))
	{
	    fuh->_bi_interpolate_intfc_states = gt_lin_comb_states;
	    fuh->_tri_interpolate_intfc_states = gt_tri_lin_comb_states;
	}
	insert_point_adjacent_to_node(pnew,newc->curve,newc->orient);

	if (debugging("experimental"))
	    extend_states_out_from_node(oldc,newc,fr,wave,dt);

	LEAVE_ADJUST_ANGLE_AT_NODE
	return YES;
#undef LEAVE_ADJUST_ANGLE_AT_NODE
}		/*end adjust_angle_at_node*/


LOCAL	double sonic_radius_at_node(
	NODE		*newn,
	Locstate	st_l,
	Locstate	st_r,
	double		dt,
	int		dim)
{
	double		distl, distr;

	distl = sonic_radius(st_l,Node_vel(newn),dt,dim);
	distr = sonic_radius(st_r,Node_vel(newn),dt,dim);
	return min(distl,distr);
}		/*end sonic_radius_at_node*/

EXPORT	double	sonic_radius(
	Locstate	state,
	double		*v,
	double		dt,
	int		dim)
{
	double	dist;
	int	i;

	for (dist = 0.0, i = 0; i < dim; i++)
		dist += sqr(vel(i,state) - v[i]);
	dist = sqrt(dist) + sound_speed(state);
	return	dt*dist;
}		/*end sonic_radius*/



LOCAL	POINT *eliminate_short_bonds_at_node(
	O_CURVE		*oldc,
	O_CURVE		*newc,
	double		dt,
	Front		*fr,
	Wave		*wave,
	double		*pdist)
{
	BOND		*b;
	INTERFACE       *intfc = newc->curve->interface;
	POINT		*p = NULL;
	double   	v[MAXD];
	double		dist = *pdist;
	double		*h = fr->rect_grid->h;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	int		dim = fr->rect_grid->dim;

	if (newc->curve->num_points < 5) /*TOLERANCE*/
		return NULL;

	b = Bond_at_node_of_o_curve(newc);
	while ((bond_length(b) < dist) ||
	       (scaled_bond_length(b,h,dim) < 0.1*MIN_SC_SEP(fr->interf))) /*TOLERANCE*/
	{
		(void) delete_point_adjacent_to_node(fr,newc->curve,
						     newc->orient);
		b = Bond_at_node_of_o_curve(newc);
		if (newc->curve->num_points == 2)
		{
			NODE	*oppn = Opp_node_of_o_curve(newc);

			if (oldc && oldc->curve &&
			    (propagation_status(oppn) != PROPAGATED_NODE))
			{
				p = Point(NULL);
	    			point_propagate(fr,(POINTER)wave,
	        			Opp_node_of_o_curve(oldc)->posn,p,
	        			Bond_at_opp_node_of_o_curve(oldc),
					oldc->curve,dt,v);
	    			if (oldc->orient != newc->orient)
					reverse_states_at_point(p,fr);
				interpolate_intfc_states(intfc) = NO;
				insert_point_adjacent_to_node(p,
					      newc->curve,newc->orient);
				interpolate_intfc_states(intfc) = sav_intrp;
				b = Bond_at_node_of_o_curve(newc);
				*pdist = min(dist,0.5*bond_length(b));
			}
			return p;
		}
	}
	return p;
}		/*end eliminate_short_bonds_at_node*/


/*
*			check_for_consistency_at_node():
*
*	This function checks for problems during the modification of a node.
*	The first test is for large angles in the curve, ie kinks.  These
*	are deleted if they exist.  
*
*	Next a test is done to make sure the components are consistent at
*	the node.  It is possible to generate tangles as the curves are
*	modified at the node and computed angles are enforced.  We are mainly
*	checking to see that adjustment of a curve ahead of the current one
*	hasn't clipped off a point on newc, leaving it in the wrong component.
*
*	Note: the current implementation of this second check may be too
*	simplistic in that it assumes that the point adjacent to a node
*	cannot be the first point on the curve to form a valid tangle.  One
*	possible upgrade would be to consider all comps around the node.  If
*	the component of the adjacent point is not one of these, then the
*	tangle is more fundamental and should probably be left alone.  This
*	still may be insufficient, but it seems unreasonable for curves
*	incident on the same node to tangle with each other near the node.
*/

LOCAL	boolean check_for_consistency_at_node(
	O_CURVE		*newc,
	double		*d,
	double		*dir,
	POINT		*pt,
	Front		*fr)
{
	CURVE		*c = newc->curve;
	POINT		*p0, *p1;
	int		i, dim = c->interface->dim;
	ORIENTATION	orient = newc->orient;

	/* First check for large angles at node */

	p0 = Node_of_o_curve(newc)->posn;
	p1 = Point_adjacent_to_node(c,orient);
	for (i = 0; i < dim; i++) d[i] = Coords(p1)[i] - Coords(p0)[i];
	while (scalar_product(dir,d,dim) < 0.0)
	{
		if (!Following_bond(Bond_at_node_of_o_curve(newc),orient))
			break;
		if (pt != NULL)
			return FUNCTION_FAILED;
		(void) delete_point_adjacent_to_node(fr,c,orient);
		p1 = Point_adjacent_to_node(c,orient);
		for (i = 0; i < dim; i++)
			d[i] = Coords(p1)[i] - Coords(p0)[i];
	}

	return FUNCTION_SUCCEEDED;
}		/*end check_for_consistency_at_node*/


LOCAL	void coords_along_ray(
	double		*coords,	/* answer */
	double		max_dist,	/* max dist along ray */
	double		*d,		/* reference dir (current, at node) */
	double		*p0,		/* start point of ray */
	double		*dir,		/* direction of ray */
	double		len,		/* target dist along ray */
	int		dim)
{
	double		mag_d, para;
	int		i;

	mag_d = mag_vector(d,dim);
	if (mag_d > 0.0)
	{
		para = 0.5*len*scalar_product(dir,d,dim)/mag_d;
		para = min(para,max_dist);
	}
	else
		para = max_dist;
	for (i = 0; i < dim; i++) coords[i] = p0[i] + para*dir[i];

	if (debugging("adjust")) 
	{
	    double r_angle = angle(d[0],d[1]);

	    print_general_vector("new point to be inserted, at ",
				 coords,dim,"\n");
	    (void) printf("para = %g, max_dist = %g, ",para,max_dist);
	    print_angle("r_angle =",r_angle,"\n");
	    if (mag_d > 0.0)
	    {
	        (void) printf("0.5*len*scalar_product(dir,d,dim)/mag_d = %g\n",
		    0.5*len*scalar_product(dir,d,dim)/mag_d);
	    }
	    else
	        (void) printf("mag_d = 0.0\n");
	}
}		/*end coords_along_ray*/

LOCAL	void ensure_states_set_at_opp_end_of_bond(
	O_CURVE		*oldc,
	O_CURVE		*newc,
	POINT		*tmp_point,
	Front		*fr,
	Wave		*wave,
	double		dt)
{
	NODE		*oppn = Opp_node_of_o_curve(newc);
	ORIENTATION	opp_orient;
	double   	v[MAXD];
	static POINT	*p = NULL;
	static boolean	first = YES;

	if ((oldc == NULL) || (oldc->curve == NULL) ||
	    (newc->curve->first != newc->curve->last) ||
	    (propagation_status(oppn) == PROPAGATED_NODE))
		return;

	opp_orient = Opposite_orient(newc->orient);
	if (tmp_point != NULL)
	{
		ft_assign(Left_state_at_node(newc->curve,opp_orient),
			left_state(tmp_point),fr->sizest);
		ft_assign(Right_state_at_node(newc->curve,opp_orient),
			right_state(tmp_point),fr->sizest);
		return;
	}

	if (first) 
	{
		first = NO;
		p = Static_point(fr->interf);
	}

	point_propagate(fr,(POINTER)wave,Opp_node_of_o_curve(oldc)->posn,p,
		Bond_at_opp_node_of_o_curve(oldc),oldc->curve,dt,v);
	if (oldc->orient != newc->orient) reverse_states_at_point(p,fr);
	
	ft_assign(Left_state_at_node(newc->curve,opp_orient),
		left_state(p),fr->sizest);
	ft_assign(Right_state_at_node(newc->curve,opp_orient),
		right_state(p),fr->sizest);
}		/*end ensure_states_set_at_opp_end_of_bond*/


LOCAL	void extend_states_out_from_node(
	O_CURVE		*oldc,
	O_CURVE		*newc,
	Front		*fr,
	Wave		*wave,
	double		dt)
{
	BOND		*bb;
	POINT		*pp = Point_adjacent_to_node(newc->curve,newc->orient);
	Locstate	sl, sr;
	double		t[MAXD], nor[MAXD], W[MAXD];
	double		pjump;

	bb = Bond_at_node_of_o_curve(newc);
	find_tangent_to_propagated_curve(pp,bb,oldc,newc,t,fr,(POINTER)wave,dt);
	if (newc->orient == POSITIVE_ORIENTATION)
	{
		nor[0] = t[1];
		nor[1] = -t[0];
	}
	else
	{
		nor[0] = -t[1];
		nor[1] = t[0];
	}
	/*
	*  The curvature is not well defined on partially
	*  propagated curves, so compute pressure jump using oldc
	*
	*  POSSIBLE UPGRADE
	*  Create curvature_at_point_on_propagated_curve().
	*/
	pjump = set_pjump_at_wave(Node_of_o_curve(oldc)->posn,
			Hyper_surf_element(Bond_at_node_of_o_curve(oldc)),
			Hyper_surf(oldc->curve),fr,nor);

	sl = Left_state_at_node_of_o_curve(newc);
	sr = Right_state_at_node_of_o_curve(newc);
	w_speed(Coords(pp),sl,sr,left_state(pp),right_state(pp),
		W,pjump,nor,wave_type(newc->curve),fr);

	n_pt_propagated(pp) = YES;
	t_pt_propagated(pp) = YES;
}		/*end extend_states_out_from_node*/

/*
*		g_check_delete_redundant_node():
*
*	Routine to check if we should actually delete
*	a node with only two incoming curves.
*/

EXPORT boolean g_check_delete_redundant_node(
	NODE		*n,
	CURVE		*c1,
	CURVE		*c2)
{
	/*Physics independent part*/
	if (!f_check_delete_redundant_node(n,c1,c2))
	    return NO;

	/*Physics dependent part*/

	/* For overtake node, we have two FORWARD/BACKWARD shocks.  To
	 * delete the node, we must invert one of them, turning it into
	 * a shock of the opposite family, which causes problems.
	 */

	return (node_type(n) == OVERTAKE_NODE) ? NO : YES;
}		/*end g_check_delete_redundant_node*/	
#endif /* defined(TWOD) */
