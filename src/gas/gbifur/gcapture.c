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
*				gcapture.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for the capture, or initialization of tracking
*	for untracked shocks. 
*	Also the node type WAVE_END_NODE which is a free end of a tracked 
*	curve is propagated here. This node can be thought of either
*	as the point where a tracked wave ends in a point of zero
*	strength or the point along a curve of nonzero strength where the
*	tracking stops and the untracked portion of the wave begins.
*
*/

#if defined(FULL_PHYSICS) && defined(TWOD)

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	void	extend_curve(Front*,Wave*,double*,double,double,
			     O_CURVE*,O_CURVE*);
LOCAL	void	find_extended_curve_tangent(O_CURVE*,O_CURVE*,double*,
					    double*,double*,Front*,Wave*);
LOCAL	void	find_extended_nt_stencil(int,Front*,CURVE*,ORIENTATION,
					 POINT**,POINT**,POINT**);
LOCAL	void	grad(POINT*,POINT*,POINT*,double,double,double,double*,double*);
LOCAL	void	init_curve_segment(Front*,Wave*,int*,double*);


/*
*			wave_end_propagate():
*/

/*ARGSUSED*/
EXPORT int wave_end_propagate(
	Front		*fr,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	Wave		*wave = wave_of_front(fr);
	O_CURVE		Cinc, New_cinc;
	int		status = GOOD_NODE;
	int		dim = fr->rect_grid->dim;
	double		t[MAXD]; 	/* tan to new_cinc at node,
					   points away from new_cinc */
	double	 	nor[MAXD];	/* new_cinc nor, oriented towards 
					   right side new_cinc */
	double		len;
	double		node_strength;
	double		incipient_strength;
	Locstate	ls,rs;

	debug_print("wave_end_prop","Entered wave_end_propagate\n");
	find_curve_with_status(oldn,&Cinc.curve,&Cinc.orient,INCIDENT);
	if (debugging("wave_end_prop"))
	{
		(void) printf("Old Node and Old Curve\n");
		print_node(oldn);
		print_o_curve(&Cinc);
	}
	find_curve_with_status(newn,&New_cinc.curve,&New_cinc.orient,INCIDENT);
	ls = Left_state_at_node_of_o_curve(&New_cinc);
	rs = Right_state_at_node_of_o_curve(&New_cinc);

		/* Preliminary propagation of node and next point */

	point_propagate(fr,(POINTER)wave,oldn->posn,newn->posn,
		Bond_at_node_of_o_curve(&Cinc),Cinc.curve,dt,Node_vel(newn));
	ft_assign(ls,left_state(newn->posn),fr->sizest);
	ft_assign(rs,right_state(newn->posn),fr->sizest);

		/* Final propagation of node */

	find_extended_curve_tangent(&Cinc,&New_cinc,t,
		&incipient_strength,&node_strength,fr,wave);
	normal(Node_of_o_curve(&Cinc)->posn,
	       Hyper_surf_element(Bond_at_node_of_o_curve(&Cinc)),
	       Hyper_surf(Cinc.curve),nor,fr);
	if (Cinc.orient == NEGATIVE_ORIENTATION)
	{
		t[0] *= -1.;
		t[1] *= -1.;
	}
	debug_print("wave_end_prop","nx ny = %g %g t[0] t[1] = %g %g\n",
		nor[0],nor[1],t[0],t[1]);
	if (node_strength > .05)
	{
		nor[0] = node_strength * nor[0] - incipient_strength * t[1];
		nor[1] = node_strength * nor[1] + incipient_strength * t[0];
		len = mag_vector(nor,dim);
		nor[0] /= len;	nor[1] /= len;
		debug_print("wave_end_prop",
			"orient %d final nx ny = %g %g len/(n+i) = %g\n",
			Cinc.orient,nor[0],nor[1],
			len/(node_strength+incipient_strength));
	}
	oblique_propagate_at_node(fr,wave,newn->posn,&Cinc,&New_cinc,
		nor,dt);
	ft_assign(ls,left_state(newn->posn),fr->sizest);
	ft_assign(rs,right_state(newn->posn),fr->sizest);

		/* Extend or contract curve at node */

	t[0] = nor[1];	t[1] = -nor[0];
	if (Cinc.orient == NEGATIVE_ORIENTATION)
	{
		t[0] *= -1.;
		t[1] *= -1.;
	}
	extend_curve(fr,wave,t,incipient_strength,
		node_strength,&Cinc,&New_cinc);
	propagation_status(newn) = PROPAGATED_NODE;
	if (debugging("wave_end_prop"))
	{
		(void) printf("New Node and New Curve\n");
		print_node(newn);
		print_o_curve(&New_cinc);
	}
	debug_print("wave_end_prop","Leaving wave_end_propagate\n");
	return status;
}		/*end wave_end_propagate*/


/*
*			snd_wave_end_propagate():
*/

/*ARGSUSED*/
EXPORT boolean snd_wave_end_propagate(
	Front		*fr,
	Front		*newfr,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	Wave		*wave = wave_of_front(fr);
	O_CURVE		Oldc, Tempc, Newc;
	static POINT	**pr = NULL, **pl = NULL;
	int		i, isgn, indx;
	int		nrad = fr->npts_tan_sten/2, dim = fr->rect_grid->dim;
	Locstate	ansr,ansl;
	double		nor[MAXD]; 	/* nor at node, directed
					   toward right side of curve */
	double		tgnt[MAXD];	/* tan at node, directed
					   away from curve */
	double		len;
	double		ds;
	double		node_strength;
	double		incipient_strength;
	static boolean	first = YES;
	static Tan_stencil *sten = NULL;

	debug_print("wave_end_propagate","Entered snd_wave_end_propagate\n");

		/* Allocate storage */

	if (first)
	{
	    first = NO;

	    sten = alloc_tan_stencil(fr,nrad);
	    uni_array(&pl,nrad,sizeof(POINT *));
	    uni_array(&pr,nrad,sizeof(POINT *));
	    for (i = 0; i < nrad; ++i)
	    {
	    	pr[i] = Static_point(fr->interf);
	    	pl[i] = Static_point(fr->interf);
	    }
	}
	sten->curvature = 0.0;

		/* Identify curves */

	find_curve_with_status(oldn,&Oldc.curve,&Oldc.orient,INCIDENT);
	find_curve_with_status(tempn,&Tempc.curve,&Tempc.orient,INCIDENT);
	find_curve_with_status(newn,&Newc.curve,&Newc.orient,INCIDENT);

		/* find tangent to curve at node */

	normal(Node_of_o_curve(&Oldc)->posn,
	       Hyper_surf_element(Bond_at_node_of_o_curve(&Oldc)),
	       Hyper_surf(Oldc.curve),nor,fr);
	find_extended_curve_tangent(&Oldc,&Tempc,tgnt,
		&incipient_strength,&node_strength,fr,wave);
	if (Oldc.orient == NEGATIVE_ORIENTATION) 
	{
		tgnt[0] *= -1.;
		tgnt[1] *= -1.;
	}
	if (node_strength > .05)
	{
		nor[0] = node_strength * nor[0] - incipient_strength * tgnt[1];
		nor[1] = node_strength * nor[1] + incipient_strength * tgnt[0];
		len = mag_vector(nor,dim);
		nor[0] /= len;	nor[1] /= len;
	}
	tgnt[0] = nor[1];	tgnt[1] = -nor[0];
	if (Oldc.orient == NEGATIVE_ORIENTATION) 
	{
		tgnt[0] *= -1.;
		tgnt[1] *= -1.;
	}
	debug_print("wave_end_propagate","nx ny %g %g tx ty %g %g len %g\n",
		nor[0],nor[1],tgnt[0],tgnt[1],len);

		/* Find three stencil states */

	ansr = Right_state_at_node_of_o_curve(&Newc);
	ansl = Left_state_at_node_of_o_curve(&Newc);

		/* Find mid state */

	sten->p[0] = tempn->posn;
	sten->hse[0] = Hyper_surf_element(Bond_at_node_of_o_curve(&Tempc));
	sten->hs[0] = Hyper_surf(Tempc.curve);
	sten->t[0] = (Tempc.orient == POSITIVE_ORIENTATION) ? 0.0 : 1.0;

		
	ds = grid_size_in_direction(tgnt,fr->rect_grid->h,dim);
	isgn = (Tempc.orient == NEGATIVE_ORIENTATION) ? -1 : 1;

		/* Find next state */
	find_extended_nt_stencil(nrad,fr,Oldc.curve,Oldc.orient,
		pl,sten->p+isgn,pr);
	for (i = 0; i < nrad; ++i)
	{
		indx = isgn*(i+1);
		hyp_solution(Coords(pr[i]),positive_component(Oldc.curve),
			Hyper_surf(Oldc.curve),POSITIVE_SIDE,fr,wave,
			right_state(sten->p[indx]),NULL);
		hyp_solution(Coords(pl[i]),negative_component(Oldc.curve),
			Hyper_surf(Oldc.curve),NEGATIVE_SIDE,fr,
			wave,left_state(sten->p[indx]),NULL);
		sten->hse[indx] = NULL;
		sten->hs[indx] = NULL;
		sten->t[indx] = 0.0;
	}

		
		/* Find prev state */
	states_at_distance_along_curve(tempn->posn,
		Bond_at_node_of_o_curve(&Tempc),Tempc.curve,
		Tempc.orient,ds,nrad,
		sten->leftst+isgn,sten->rightst+isgn,
		sten->hs+isgn,sten->hse+isgn,sten->t+isgn,
		sten->p+isgn,fr);

		/* Solve tang. difference eqn */

	sten->newhs = Hyper_surf(Newc.curve);
	sten->dir = tgnt;
	npt_tang_solver(ds,dt,sten,ansl,ansr,fr);

	debug_print("wave_end_prop","Leaving snd_wave_end_propagate\n");
	return YES;
}		/*end snd_wave_end_propagate*/

#define is_regular_quad(icoords,wave) 					\
	(blk_el0_is_bilinear(&Rect_blk_el0(icoords,wave)))

/*
*			g_shock_capture():
*
*	This is routine provides the main control flow for the initialization
*	of tracked shocks.
*/

EXPORT void g_shock_capture(
	Front		*fr)
{
	Wave		*wave = wave_of_front(fr);
	int		ix,iy,icoords[MAXD],icoords00[MAXD];
	int		icoords01[MAXD],icoords10[MAXD],icoords11[MAXD];
	int		xmax = fr->rect_grid->gmax[0];
	int		ymax = fr->rect_grid->gmax[1];
	double		n[MAXD];
	double		hx = fr->rect_grid->h[0];
	double		hy = fr->rect_grid->h[1];
	double		p00,p01,p10,p11;
	double		grad_x_p,grad_y_p;
	double		len;
	Locstate	s00,s01,s10,s11;


	debug_print("capture","Entered g_shock_capture\n");

	for (ix = 1; ix < xmax; ++ix)
	{
	    icoords00[0] = icoords01[0] = ix-1;
	    icoords10[0] = icoords11[0] = icoords[0] = ix;
	    for (iy = 1; iy < ymax; ++iy)
	    {
	    	icoords00[1] = icoords10[1] = iy-1;
	    	icoords01[1] = icoords11[1] = icoords[1] = iy;

	    	if (!is_regular_quad(icoords00,wave))
	    		continue;
	    	if (!is_regular_quad(icoords01,wave))
	    		continue;
	    	if (!is_regular_quad(icoords10,wave))
	    		continue;
	    	if (!is_regular_quad(icoords,wave))
	    		continue;

		/*  Compute grad p and normal nx,ny in direction - grad p */

		s00 = Rect_state(icoords00,wave);
		s01 = Rect_state(icoords01,wave);
		s10 = Rect_state(icoords10,wave);
		s11 = Rect_state(icoords11,wave);
		p00 = pressure(s00);
		p01 = pressure(s01);
		p10 = pressure(s10);
		p11 = pressure(s11);
		grad_x_p = (p10 + p11 - p00 - p01) / (2. * hx);
		grad_y_p = (p01 + p11 - p00 - p10) / (2. * hy);
		len = hypot(grad_x_p,grad_y_p);
		if (len <= .001) continue;
		n[0] = - grad_x_p / len;
		n[1] = - grad_y_p / len;
		init_curve_segment(fr,wave,icoords,n);
	    }
	}
	debug_print("capture","Leaving g_shock_capture\n");
}		/*end g_shock_capture*/

/*
*			extend_curve():
*
*	This routine decides whether to extend or contract a curve which 
*	terminates at a WAVE_END_NODE. The solution in a neighborhood of 
*	the curve end is examined to determine whether there is an 
*	incipient wave of some minimum strength at that location. 
*	The prefered normal and tangent directions of this incipient 
*	wave are also determined and if the extension is made, it is 
*	made in the direction determined by these directions, for a 
*	length of one front spacing bond.
*	If the wave is very weak at the node and also at the point one
*	bond inward, then the node is moved inward by one bond.
*/

#define delete_from_linked_list(n)					\
	if (prev_node(n)) next_node(prev_node(n)) = next_node(n);	\
	if (next_node(n)) prev_node(next_node(n)) = prev_node(n);


LOCAL void extend_curve(
	Front		*fr,
	Wave		*wave,
	double		*t,	/* tan to c at node, directed away from c */
	double		incipient_strength,
	double		node_strength,
	O_CURVE		*oldc,
	O_CURVE		*c)
{
	double		adj_shock_strength;/* shock pressure ratio adj to node*/
	double		pr_l,pr_r;	   /* pressure left and right of node */
	double		nor[MAXD];	   /* normal to c, directed to right */
	double		normal_offset;	   /* location left, right states */
	double		new_bond_length;   /* length of extra bond */
	double		coords[MAXD];
	double		spacing;
	double		*h = fr->rect_grid->h;
	NODE		*n;
	POINT		*adj_point;
	Locstate	l_state,r_state;
	int		i, dim = fr->rect_grid->dim;

	debug_print("capture","Entered extend_curve\n");
	if (wave_type(c->curve) < FIRST_PHYSICS_WAVE_TYPE) return;
	if (is_scalar_wave(wave_type(c->curve)))
	{
		spacing = Front_spacing(fr,GENERAL_WAVE);
		goto contact;
	}
	else
		spacing = Front_spacing(fr,VECTOR_WAVE);
	n = Node_of_o_curve(c);
	if (node_strength > expand_threshold(fr) || 
		incipient_strength > expand_threshold(fr))
	{

			/* Add new bond */

		debug_print("capture","Insert new point to extend curve\n");
		interpolate_intfc_states(c->curve->interface) = NO;
		insert_point_adjacent_to_node(Point(Coords(n->posn)),
			c->curve,c->orient);
		new_bond_length = .5 / scaled_hypot(t,h,dim);
		for (i = 0; i < dim; ++i)
			Coords(n->posn)[i] += new_bond_length * t[i];
		adj_point = Point_adjacent_to_node(c->curve,c->orient);
		l_state = left_state_at_point_on_curve(adj_point,
				Bond_at_node_of_o_curve(c),c->curve);
		r_state = right_state_at_point_on_curve(adj_point,
				Bond_at_node_of_o_curve(c),c->curve);
		ft_assign(l_state,Left_state_at_node_of_o_curve(c),fr->sizest);
		ft_assign(r_state,Right_state_at_node_of_o_curve(c),fr->sizest);
		nor[0] = -t[1];	nor[1] = t[0];
		if (c->orient == NEGATIVE_ORIENTATION)
		{
			nor[0] *= -1.; nor[1] *= -1.;
		}
		normal_offset = .25 / scaled_hypot(nor,h,dim);
		for (i = 0; i < dim; ++i)
			coords[i] = Coords(n->posn)[i]+normal_offset*nor[i];
		hyp_solution(coords,positive_component(c->curve),
			Hyper_surf(oldc->curve),POSITIVE_SIDE,fr,wave,
			Right_state_at_node_of_o_curve(c),NULL);
		for (i = 0; i < dim; ++i)
			coords[i] = Coords(n->posn)[i]-normal_offset*nor[i];
		hyp_solution(coords,negative_component(c->curve),
			Hyper_surf(oldc->curve),NEGATIVE_SIDE,fr,wave,
			Left_state_at_node_of_o_curve(c),NULL);
		if (debugging("capture")) print_o_curve(c);
		debug_print("capture","Leaving extend_curve after extension\n");
		return;
	}

			/* Don't delete single bond if ... */

	if (	c->curve->first == c->curve->last &&
		(node_type(c->curve->start) != WAVE_END_NODE ||
		node_type(c->curve->end) != WAVE_END_NODE )	) 
		return;
	if (	c->curve->first == c->curve->last &&
		(propagation_status(Opp_node_of_o_curve(c)) != PROPAGATED_NODE))
		return;

			/* Test for wave strength */

	adj_point = Point_adjacent_to_node(c->curve,c->orient);
	l_state = left_state_at_point_on_curve(adj_point,
			Bond_at_node_of_o_curve(c),c->curve);
	r_state = right_state_at_point_on_curve(adj_point,
			Bond_at_node_of_o_curve(c),c->curve);
	pr_l = pressure(l_state);
	pr_r = pressure(r_state);
	adj_shock_strength = fabs(2. * (pr_l - pr_r) / (pr_l + pr_r));
	if (debugging("capture"))
	{
		verbose_print_state("left",l_state);
		verbose_print_state("right",r_state);
		(void) printf("pr_l,pr_r,adj str,node str = %g %g %g %g\n",pr_l,
			pr_r,adj_shock_strength,node_strength);
	}
	if (adj_shock_strength < contract_threshold(fr) &&
		node_strength < contract_threshold(fr))
	{
			/* Delete full curve if short */

		if (c->curve->first == c->curve->last)
		{
			delete_from_linked_list(c->curve->start);
			delete_from_linked_list(c->curve->end);
/**** THIS MUST BE WRONG, OLD INTERFACE SHOULD NOT BE MODIFIED
			delete_from_linked_list(oldc->curve->start);
			delete_from_linked_list(oldc->curve->end);
***/
			(void) delete_curve(c->curve);
			(void) delete_node(c->curve->start);
			(void) delete_node(c->curve->end);
			c->curve = NULL;
			if (debugging("capture"))
			{
				(void) printf("Interface after delete curve\n");
				print_interface(c->curve->interface);
			}
			return;
		}

			/*  Delete bond at node */

		for (i = 0; i < dim; ++i)
			Coords(n->posn)[i] = Coords(adj_point)[i];
		(void) delete_point_adjacent_to_node(fr,c->curve,c->orient);
		ft_assign(Left_state_at_node_of_o_curve(c),l_state,fr->sizest);
		ft_assign(Right_state_at_node_of_o_curve(c),r_state,fr->sizest);
		debug_print("capture","Leaving extend_curve after contraction\n");
		return;
	}

		/* Delete short bond at node */
	
	if (scaled_bond_length(Bond_at_node_of_o_curve(c),h,dim) < .2 * spacing)
		(void) delete_point_adjacent_to_node(fr,c->curve,c->orient);
	debug_print("capture","Leaving extend_curve no change\n");
	return;
contact:
		/* TODO: extend in contact case */
	debug_print("capture","Leaving extend_curve no change\n");
	return;
}		/*end exten_curve*/

/*
*			find_extended_curve_tangent():
*
*	Given two triangular elements, the pressure gradient determines a
*	direction in each element.  The magnitude
*	of the variation in the pressure gradient determines the
*	incipient curve strength in each triangle, taken to be zero in case
*	the variation in the relevant family is that of a rarefaction wave.
*	From these curve strengths, a new direction is chosen by interpolation
*	between the directions defined by the pressure gradients,
*	weighted by the pair of incipient curve strengths. Moreover
*	a composite value for the incipient curve strength is determined.
*
*	TODO: In case the curve is a contact, use density or vorticity in place
*	of pressure.
*/

LOCAL void find_extended_curve_tangent(
	O_CURVE		*oldc,
	O_CURVE		*c,
	double		*t,
	double		*incipient_strength,
	double		*node_strength,
	Front		*fr,
	Wave		*wave)
{
	static POINT	*pl = NULL,*pm = NULL,
	                *pr = NULL; /* On left, mid right of extended wave */
	static boolean	first = YES;
	double		press_l; /* pressures left mid right of extended wave */
	double		press_m;
	double		press_r;
	double		press_el;    /* pressures left and right of curve end */
	double		press_er;
	double		grad_press_lx;
	double		grad_press_ly;
	double		grad_press_rx;
	double		grad_press_ry;
	double		nx,ny;	      /* curve normal directed to low pr side */
	double		len,len_l,len_r;	/* lengths of uni_arrays */

	if (first)
	{
		first = NO;
		pl = Static_point(fr->interf);
		pm = Static_point(fr->interf);
		pr = Static_point(fr->interf);
	}
	find_extended_nt_stencil(1,fr,c->curve,c->orient,&pl,&pm,&pr);

	hyp_solution(Coords(pr),positive_component(c->curve),
		Hyper_surf(oldc->curve),POSITIVE_SIDE,fr,
		wave,right_state(pr),NULL);
	hyp_solution(Coords(pm),positive_component(c->curve),
		Hyper_surf(oldc->curve),POSITIVE_SIDE,fr,
		wave,right_state(pm),NULL);
	hyp_solution(Coords(pl),positive_component(c->curve),
		Hyper_surf(oldc->curve),NEGATIVE_SIDE,fr,
		wave,right_state(pl),NULL);

	press_r =  pressure(right_state(pr));
	press_m =  pressure(right_state(pm));
	press_l =  pressure(right_state(pl));
	press_el = pressure(Left_state_at_node_of_o_curve(c));
	press_er = pressure(Right_state_at_node_of_o_curve(c));
	*node_strength = 2. * (press_el - press_er)/(press_el + press_er);
	*node_strength = (wave_type(c->curve) == FORWARD_SHOCK_WAVE) ?
		max(*node_strength,0.) : max(-*node_strength,0.);
	if (debugging("capture"))
	{
		verbose_print_state("left",Left_state_at_node_of_o_curve(c));
		verbose_print_state("right",Right_state_at_node_of_o_curve(c));
		(void) printf("node str %g\n",*node_strength);
	}

	grad(Node_of_o_curve(c)->posn,pm,pr,press_er,press_m,press_r,
		&grad_press_rx,&grad_press_ry);
	grad(Node_of_o_curve(c)->posn,pm,pl,press_el,press_m,press_l,
		&grad_press_lx,&grad_press_ly);

		/* Find normal oriented toward low pressure (ahead) side */

	find_tangent_to_curve(Node_of_o_curve(oldc)->posn,
		Bond_at_node_of_o_curve(oldc),oldc->curve,oldc->orient,
		t,fr);
	t[0] *= -1.;	t[1] *= -1.;
	nx = -t[1];	ny =  t[0];
	if (	(wave_type(c->curve) == FORWARD_SHOCK_WAVE && 
		c->orient == NEGATIVE_ORIENTATION) ||
		(wave_type(c->curve) == BACKWARD_SHOCK_WAVE && 
		c->orient == POSITIVE_ORIENTATION)	)
	{
		nx *= -1; ny *= -1;
	}
	if (debugging("capture"))
	{
		(void) printf("press_r,m,l,er,el %g %g %g %g %g\n",
			      press_r,press_m,press_l,press_er,press_el);
		(void) printf("nx,ny = %g %g\n",nx,ny);
		print_o_curve(c);
	}

		/* Select shock portion of pressure gradient only */

	if (	(nx * grad_press_lx + ny * grad_press_ly) >= 0. &&
		(nx * grad_press_rx + ny * grad_press_ry) >= 0.)
	{
		debug_print("capture","normal from neither tri\n");
		*incipient_strength = 0.;
		return;
	}
	len_r = hypot(grad_press_rx,grad_press_ry);
	len_l = hypot(grad_press_lx,grad_press_ly);
	if (	(nx * grad_press_lx + ny * grad_press_ly) >= 0.)
	{
		debug_print("capture","normal from right tri only\n");
		*incipient_strength = 
			2. * (press_m - press_r) / (press_r + press_m);
		nx = -grad_press_rx / len_r;
		ny = -grad_press_ry / len_r;
	}
	else if ((nx * grad_press_rx + ny * grad_press_ry) >= 0.)
	{
		debug_print("capture","normal from left tri only\n");
		*incipient_strength = 
			2. * (press_l - press_m) / (press_l + press_m);
		nx = -grad_press_lx / len_l;
		ny = -grad_press_ly / len_l;
	}
	else
	{
		debug_print("capture","normal from both tri's\n");
		*incipient_strength = 
			2. * (press_l - press_r) / (press_l + press_r);
		nx = -grad_press_lx - grad_press_rx;
		ny = -grad_press_ly - grad_press_ry;
		len = hypot(nx,ny);
		nx /= len;
		ny /= len;
	}
	*incipient_strength = (wave_type(c->curve) == FORWARD_SHOCK_WAVE) ?
		max(*incipient_strength,0.) : max(-*incipient_strength,0.);
	t[0] =  ny;
	t[1] = -nx;
	if (	(wave_type(c->curve) == FORWARD_SHOCK_WAVE && 
		c->orient == NEGATIVE_ORIENTATION) ||
		(wave_type(c->curve) == BACKWARD_SHOCK_WAVE && 
		c->orient == POSITIVE_ORIENTATION)	)
	{
		t[0] *= -1; t[1] *= -1;
	}
	if (debugging("capture"))
	{
		(void) printf("grad_press_lx,ly,rx,ry = %g %g  %g %g\n",
			grad_press_lx,grad_press_ly,
			grad_press_rx,grad_press_ry);
		(void) printf("final values: nx ny t[0] t[1] inc_str ");
		(void) printf("%g %g %g %g %g\n",
			      nx,ny,t[0],t[1],*incipient_strength);
	}
}		/*end find_extended_curve_tangent*/
	
/*
*			find_extended_nt_stencil():
*
*	At a WAVE_END_NODE, two triangular elements are constructed with
*	normal and tangential orientation. The curve is extended by a bond
*	length to give the line which is the common border of the two
*	triangles, and at the end of this bond, a normal of one bond length
*	in each direction is constructed, thereby forming the sides of two
*	right triangles. The three points have their positions set to 
*	the vertices of these two triangles.
*/

LOCAL void find_extended_nt_stencil(
	int		npts,
	Front		*fr,
	CURVE		*c,
	ORIENTATION	c_orient,
	POINT		**pl,
	POINT		**pm,
	POINT		**pr)
{
	double		t[MAXD];       /* tang to wave at node, with c_orient */
	double		n[MAXD];       /* normal, towards right side*/
	double		ds, dn;	       /* fr spacing in nor or tan direction */
	double		*h = fr->rect_grid->h;
	int		i, j, dim = fr->rect_grid->dim;
	int		indx, isgn;
	NODE		*node;

	isgn = (c_orient == NEGATIVE_ORIENTATION) ? -1 : 1;

	node = Node_of(c,c_orient);
	find_tangent_to_curve(node->posn,Bond_at_node(c,c_orient),c,
		              c_orient,t,fr);
	normal(node->posn,Hyper_surf_element(Bond_at_node(c,c_orient)),
	       Hyper_surf(c),n,fr);
	ds = grid_size_in_direction(t,h,dim);
	dn = grid_size_in_direction(n,h,dim);

	for (j = 0; j < npts; ++j)
	{
		indx = j*isgn;
		for (i = 0; i < dim; ++i)
		{
			Coords(pm[indx])[i] =
				Coords(node->posn)[i] - (0.5+j)*ds * t[i];
			Coords(pl[indx])[i] = Coords(pm[indx])[i] - dn * n[i];
			Coords(pr[indx])[i] = Coords(pm[indx])[i] + dn * n[i];
		}
	}
}		/*end find_extended_nt_stencil*/

/*
*			grad():
*
*	This routine determines a gradient, given the values of some variable
*	at any three points, assumed not to be colinear
*/

LOCAL void grad(
	POINT		*p1,
	POINT		*p2,
	POINT		*p3,
	double		f1,
	double		f2,
	double		f3,
	double		*grad_x,
	double		*grad_y)
{
	double		x1 = Coords(p1)[0];
	double		x2 = Coords(p2)[0];
	double		x3 = Coords(p3)[0];
	double		y1 = Coords(p1)[1];
	double		y2 = Coords(p2)[1];
	double		y3 = Coords(p3)[1];
	double		xstar,ystar,fstar,slope,inv_slope;

	if (debugging("grad"))
	{
		(void) printf("x,y1 x,y2 x,y3 %g %g  %g %g  %g %g\n",
			      x1,y1,x2,y2,x3,y3);
		(void) printf("f1,2,3 = %g %g %g\n",f1,f2,f3);
	}
	
	if (fabs(y3 -y2) > .1 * fabs(x3 - x2) )
	{
		inv_slope = (x3 - x2) / (y3 - y2);
		xstar = x2 + (y1 -y2) * inv_slope;
		fstar = f2 + (f3 - f2) * (y1 - y2) / (y3 - y2);
		*grad_x = (f1 - fstar) / (x1 - xstar);
		debug_print("grad","inv_s xstar fstar = %g %g %g 1st option\n",
			inv_slope,xstar,fstar);
	}
	else if (fabs(y1 -y3) > .1 * fabs(x1 - x3) )
	{
		inv_slope = (x1 - x3) / (y1 - y3);
		xstar = x3 + (y2 -y3) * inv_slope;
		fstar = f3 + (f1 - f3) * (y2 - y3) / (y1 - y3);
		*grad_x = (f2 - fstar) / (x2 - xstar);
		debug_print("grad","inv_s xstar fstar = %g %g %g 2nd option\n",
			inv_slope,xstar,fstar);
	}
	else
	{
		inv_slope = (x2 - x1) / (y2 - y1);
		xstar = x1 + (y3 -y1) * inv_slope;
		fstar = f1 + (f2 - f1) * (y3 - y1) / (y2 - y1);
		*grad_x = (f2 - fstar) / (x2 - xstar);
		debug_print("grad","inv_s xstar fstar = %g %g %g 3rd option\n",
			inv_slope,xstar,fstar);
	}
	if (fabs(x3 -x2) > .1 * fabs(y3 - y2) )
	{
		slope = (y3 - y2) / (x3 - x2);
		ystar = y2 + (x1 -x2) * slope;
		fstar = f2 + (f3 - f2) * (x1 - x2) / (x3 - x2);
		*grad_y = (f1 - fstar) / (y1 - ystar);
		debug_print("grad","inv_slope ystar fstar = %g %g %g 1st option\n",
			slope,ystar,fstar);
	}
	else if (fabs(x1 -x3) > .1 * fabs(y1 - y3) )
	{
		slope = (y1 - y3) / (x1 - x3);
		ystar = y3 + (x2 -x3) * slope;
		fstar = f3 + (f1 - f3) * (x2 - x3) / (x1 - x3);
		*grad_y = (f2 - fstar) / (y2 - ystar);
		debug_print("grad","slope ystar fstar = %g %g %g 2nd option\n",
			slope,ystar,fstar);
	}
	else
	{
		slope = (y2 - y1) / (x2 - x1);
		ystar = y1 + (x3 -x1) * slope;
		fstar = f1 + (f2 - f1) * (x3 - x1) / (x2 - x1);
		*grad_y = (f2 - fstar) / (y2 - ystar);
		debug_print("grad","slope ystar fstar = %g %g %g 3rd option\n",
			slope,ystar,fstar);
	}
	debug_print("grad","*grad_x,y = %g %g\n",*grad_x,*grad_y);
}		/*end grad*/

/*
*			init_curve_segment():
*
*	This routine tests for the size of an incipient shock in the
*	given mesh block and with given orientation. If there is none
*	such or if the strength is too small (relative to some threshold)
*	then nothing is done. Otherwise the routine adds a (short) 
*	curve in the mesh block ix,iy, ie in the block 
*	cell_center(ix-1,i,gr) < coords[i] < cell_center(ix,i,gr),
*						gr = front->rect_grid.
*	or dual lattice mesh block, which
*	is one bond long and whose node types are each WAVE_END_NODEs. The curve
*	is assumed to be a shock, and its wave_type is set accordingly.
*/


LOCAL void init_curve_segment(
	Front		*fr,
	Wave		*wave,
	int		*icoords,
	double		*n)	  	/* curve normal, oriented to low
					   press (ahead) side */
{
	double		length;
	double		t[MAXD];	       /* unit tangent to curve */
	double		s[MAXD],e[MAXD];       /* location of curve start end */
	double		coordsa[MAXD];	       /* location ahead behind curve */
	double		coordsb[MAXD];
	double		vs[MAXD],ve[MAXD];     /* node velocities */
	double		pl,pr;		       /* left,right pressures */
	double		s_strength,e_strength; /* start, end shock strengths */
	double		coords[MAXD];
	double		*h = fr->rect_grid->h;
	NODE		*ns,*ne;
	CURVE		*c;
	COMPONENT	comp;
	int	i,dim = fr->interf->dim;
	static Locstate left = NULL, right = NULL, lss = NULL, rss = NULL,
	                les = NULL, res = NULL;
	static boolean	first = YES;

	if (first)
	{
		first = NO;
		alloc_state(fr->interf,&left,fr->sizest);
		alloc_state(fr->interf,&right,fr->sizest);
		alloc_state(fr->interf,&lss,fr->sizest);
		alloc_state(fr->interf,&rss,fr->sizest);
		alloc_state(fr->interf,&les,fr->sizest);
		alloc_state(fr->interf,&res,fr->sizest);
	}

		/* Find start and end of curve */

	t[0] = -n[1];
	t[1] =  n[0];
	length = .25 / scaled_hypot(t,h,dim);
	for (i = 0; i < dim; ++i)
	{
		coords[i] = cell_edge(icoords[i],i,fr->rect_grid);
		s[i] = coords[i] - length * t[i];
		e[i] = coords[i] + length * t[i];
	}

	length = .25 / scaled_hypot(n,h,dim);
	for (i = 0; i < dim; ++i)
	{
		coordsa[i] = s[i] + n[i] * length;
		coordsb[i] = s[i] - n[i] * length;
	}
	comp = Rect_comp(icoords,wave);
	hyp_solution(coordsa,comp,NULL,UNKNOWN_SIDE,fr,wave,right,NULL);
	hyp_solution(coordsb,comp,NULL,UNKNOWN_SIDE,fr,wave,left,NULL);
	w_speed(s,left,right,lss,rss,vs,0.0,n,FORWARD_SHOCK_WAVE,fr);
	pl = pressure(left);
	pr = pressure(right);
	s_strength = 4. * (pl - pr) / (pl + pr);
	for (i = 0; i < dim; ++i)
	{
		coordsa[i] = e[i] + n[i] * length;
		coordsb[i] = e[i] - n[i] * length;
	}
	hyp_solution(coordsa,comp,NULL,UNKNOWN_SIDE,fr,wave,right,NULL);
	hyp_solution(coordsb,comp,NULL,UNKNOWN_SIDE,fr,wave,left,NULL);
	w_speed(e,left,right,les,res,ve,0.0,n,FORWARD_SHOCK_WAVE,fr);
	pl = pressure(left);
	pr = pressure(right);
	e_strength = 4. * (pl - pr) / (pl + pr);
	if (s_strength + e_strength <= 2. * init_threshold(fr))
		return;

		/* Test for strength of incipient shock */
	
	pl = pressure(lss);
	pr = pressure(rss);
	s_strength = fabs(4. * (pl - pr) / (pl + pr));
	pl = pressure(les);
	pr = pressure(res);
	e_strength = fabs(4. * (pl - pr) / (pl + pr));
	if (s_strength + e_strength <= 2. * init_threshold(fr))
		return;

		/* Make nodes */

	debug_print("capture","Start constr new curve icoords = %d %d\n",
		icoords[0],icoords[1]);
	if (debugging("capture"))
	{
		(void) printf("n = %g %g\n",n[0],n[1]);
		(void) printf("end state pressures, l,r for forward wave %g %g\n",
			pl,pr);
		(void) printf("states before Riemann sol.\n");
		verbose_print_state("left",left);
		verbose_print_state("right",right);
	}
	debug_print("capture","s & e strength = %g %g\n",s_strength,e_strength);
	for (i = 0; i < dim; ++i) coords[i] = s[i];
	ns = make_node(Point(coords));
	for (i = 0; i < dim; ++i) coords[i] = e[i];
	ne = make_node(Point(coords));
	node_type(ns) = node_type(ne) = WAVE_END_NODE;
	propagation_status(ns) = propagation_status(ne) = PROPAGATED_NODE;
	for (i = 0; i < dim; ++i)
	{
		Node_vel(ns)[i] = vs[i];
		Node_vel(ne)[i] = ve[i];
	}

		/* Make curve */

	c = make_curve(comp,comp,ns,ne);
	wave_type(c) = FORWARD_SHOCK_WAVE;
	start_status(c) = end_status(c) = INCIDENT;
	set_no_tan_propagate(c);
	correspond_hyper_surf(c) = NULL;
	if (debugging("capture"))
	{
		(void) printf("Newly captured curve:");
		print_curve(c);
	}

		/* Set start and end states */

	ft_assign(left_start_state(c),lss,fr->sizest);
	ft_assign(right_start_state(c),rss,fr->sizest);
	ft_assign(left_end_state(c),les,fr->sizest);
	ft_assign(right_end_state(c),res,fr->sizest);
	if (debugging("capture"))
	{
		(void) printf("Newly captured curve states:");
		verbose_print_state("left start",left_start_state(c));
		verbose_print_state("right start",right_start_state(c));
		verbose_print_state("left end",left_end_state(c));
		verbose_print_state("right end",right_end_state(c));
	}
}		/*end init_curve_segment*/
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */
