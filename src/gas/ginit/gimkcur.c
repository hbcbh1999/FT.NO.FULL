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
*				gimkcur.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains functions for the construction of non-trivial
*	geometric curves needed for the initialization of an interface.
*	
*	The following functions are contained in this file.
*
*		make_rotated_curve()
*		g_make_fourier_curve()
*		make_random_curve()
*		g_make_ellipse()
*		make_ramp()
*		make_vertical_axis_curve()
*		make_neumann_curve()
*/

#if defined(TWOD)

#include <ginit/ginit.h>


	/* LOCAL Function Declarations */
LOCAL	boolean	increment_theta(double*,double,double,double,ELLIPSOID*,
				double*,double);
LOCAL	void	normal_to_origin(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				 double*,Front*);
LOCAL	void	tangent_to_origin(POINT*,BOND*,CURVE*,double*,Front*);


/*	
*			make_rotated_curve():
*
*	Given a function f(x) defined for x >= 0, this function constructs
*	the curve which corresponds to that portion of the graph of f(x)
*	rotated by the angle theta0 in the counter clockwise direction
*	which lies in the rectangle XL <= x <= XU, YL <= y <= YU.
*	Points are inserted with approximately equal separation with
*	respect to the arc length along f(x) according to the metric
*	defined by the ellipse with axes of length 2*hx and 2*hy.
*/



EXPORT	void make_rotated_curve(
	double		*h,
	double		theta0,
	double		*L,
	double		*U,
	double		(*f)(double,POINTER),
	double		(*fprime)(double,POINTER),
	POINTER		params,
	NODE		*node,
	COMPONENT	lcomp,
	COMPONENT	rcomp,
	CURVE		**curve,
	ORIENTATION	orient,
	int		w_type,
	int		status,
	int		opp_status,
	double		surf_ten,
	Front		*fr)
{
	double		s, ds, p[MAXD], p0[MAXD];
	double		slow, smid, shigh;
	double		cos_theta0, sin_theta0, temp;
	double		mid[MAXD], len[MAXD], high[MAXD];
	double		d[MAXD];
	double		n[MAXD];
	POINT		*pt;
	BOND		*b;
	NODE		*opp_node, *ns, *ne;
	int		i, dim = node->interface->dim;

	if (dim != 2)
	{
	    screen("ERROR in make_rotated_curve(), "
		   "dim = %d != 2 not supported\n",dim);
	    clean_up(ERROR);
	}

	cos_theta0 = cos(theta0);
	sin_theta0 = sin(theta0);
	n[0] = cos_theta0;
	n[1] = sin_theta0;
	if (!intersect_ray_with_boundary(Coords(node->posn),n,L,U,p,dim)) 
	{
	    screen("ERROR in make_rotated_curve, Point out of domain\n");
	    print_general_vector("Coords(node->posn) = ",Coords(node->posn),
			         dim,"\n");
	    print_general_vector("L = ",L,dim,"\n");
	    print_general_vector("U = ",U,dim,"\n");
	    clean_up(ERROR);
	}
	opp_node = make_node(Point(p));
	if (orient == POSITIVE_ORIENTATION)
	{
	    ns = node;
	    ne = opp_node;
	}
	else
	{
	    ne = node;
	    ns = opp_node;
	}
	*curve = make_curve(lcomp,rcomp,ns,ne);
	wave_type(*curve) = w_type;
	surface_tension(*curve) = surf_ten;
	if (orient == POSITIVE_ORIENTATION)
	{
	    start_status(*curve) = status;
	    end_status(*curve) = opp_status;
	}
	else
	{
	    end_status(*curve) = status;
	    start_status(*curve) = opp_status;
	}

	if (f == NULL)
	    return;

	b = (*curve)->first;
	for (i = 0; i < dim; i++)
	    p0[i] = Coords(node->posn)[i];
	temp = (*fprime)(0.,params);
	d[0] = cos_theta0 - temp*sin_theta0;
	d[1] = sin_theta0 + temp*cos_theta0;
	s = 0.5/scaled_hypot(d,h,dim);
	temp = (*f)(s,params);
	p[0] = p0[0] + s*cos_theta0 - temp*sin_theta0;
	p[1] = p0[1] + s*sin_theta0 + temp*cos_theta0;
	ds = 0.0;
	while (L[0] <= p[0] && p[0] <= U[0] && L[1] <= p[1] && p[1] <= U[1]) 
	{
	    if (insert_point_in_bond(Point(p),b,*curve) != FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in make_rotated_curve(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	    b = b->next;
	    temp = (*fprime)(s,params);
	    d[0] = cos_theta0 - temp*sin_theta0;
	    d[1] = sin_theta0 + temp*cos_theta0;
	    ds = 0.5/scaled_hypot(d,h,dim);
	    s += ds;
	    temp = (*f)(s,params);
	    p[0] = p0[0] + s*cos_theta0 - 0.5*temp*sin_theta0;
	    p[1] = p0[1] + s*sin_theta0 + 0.5*temp*cos_theta0;
	}
	slow = s - ds;
	shigh = s;
	smid = 0.5*(slow + shigh);
	for (i = 0; i < dim; i++)
	{
	    high[i] = p[i];
	    mid[i] = 0.5*(U[i] + L[i]);
	    len[i] = 0.5*(U[i] - L[i]);
	}
	while ((fabs(high[0] - mid[0]) > len[0] + 0.001*h[0]) || 
		(fabs(high[1] - mid[1]) > len[1] + 0.001*h[1])) 
	{
	    temp = (*f)(s,params);
	    p[0] = p0[0] + smid*cos_theta0 - 0.5*temp*sin_theta0;
	    p[1] = p0[1] + smid*sin_theta0 + 0.5*temp*cos_theta0;
	    if ((L[0]<=p[0]) && (p[0]<=U[0]) && (L[1]<=p[1]) && (p[1]<=U[1]))
		slow = smid;
	    else 
	    {
	    	shigh = smid;
	    	for (i = 0; i < dim; i++)
	    	    high[i] = p[i];
	    }
	    smid = 0.5*(slow + shigh);
	}
	for (i = 0; i < dim; i++)
	{
	    p[i] = max(high[i],L[i]);
	    p[i] = min(p[i],U[i]);
	    Coords(ne->posn)[i] = p[i];
	}
	pt = (*curve)->last->start;
	if (_scaled_separation(p,Coords(pt),h,dim) < 0.01)
	    (void) delete_point_adjacent_to_node(fr,*curve,
						 NEGATIVE_ORIENTATION);
}		/*end make_rotated_curve*/

/*
*			g_make_ellipse():
*
*	make_ellipsoid constructs the ellipsoid with center
*	ellip->cen and radii = ellip->rad.
*	The ellipsoid may be rotated by the orthogonal transformation
*	ellip->Q.
*	
*	The input variables are:
*
*	ELLIPSOID	*ellip		Describes the geometry of the ellipsoid
*	ORIENTATION	nor_or 		Normal vector direction
*					POSITIVE_ORIENTATION = outward normal
*					NEGATIVE_ORIENTATION = inward normal
*	COMPONENT	compin,		Component numbers for inside and
*			compout		outside of ellipsoid
*	Front		*front		Front data structure
*	int		closed		Full ellisoid if closed = YES
*/

EXPORT HYPER_SURF	*g_make_ellipse(
	ELLIPSOID	*ellip,
	COMPONENT	compin,
	COMPONENT	compout,
	Front		*front)
{
	RECT_GRID   *gr = front->rect_grid;
	int	    i, dim = gr->dim;
	boolean	    done;
	ORIENTATION nor_or = ellip->nor_orient;
	int	    closed = ellip->closed;
	double	    theta[MAXD-1];
	double	    *h = gr->h, *L = gr->L, *U = gr->U;
	double	    coords[MAXD];
	double	    sgn;
	double	    ThetaS = ellip->ThetaS[0], ThetaE = ellip->ThetaE[0];
	CURVE	    *cur;
	COMPONENT    l_comp, r_comp;
	NODE	    *ns, *ne;
	static const double space = 0.01;/*TOLERANCE*/

	debug_print("ellip","Entered g_make_ellipse()\n");
	if (dim != 2)
	{
	    screen("ERROR in g_make_ellipse(), dim = %d != 2 not supported\n");
	    clean_up(ERROR);
	}
	ellip->hs = NULL;

	if (nor_or == NEGATIVE_ORIENTATION)
	{
	    sgn = -1.0;
	    l_comp = compout;
	    r_comp = compin;
	}
	else
	{
	    sgn = 1.0;
	    l_comp = compin;
	    r_comp = compout;
	}
	if (closed == YES)
	{
	    if (nor_or == NEGATIVE_ORIENTATION)
	    {
		ThetaS = normalized_angle(ThetaS);
		ThetaE = ThetaS - 2.0*PI;
	    }
	    else
	    {
		ThetaS = normalized_angle(ThetaS);
		ThetaE = ThetaS + 2.0*PI;
	    }
	    coords_on_pert_ellipsoid(ellip->ThetaS,coords,ellip);
	    ns = make_node(Point(coords));
	    node_type(ns) = CLOSED_NODE;
	    ne = ns;
	}
	else
	{
	    ThetaS = normalized_angle(ThetaS);
	    ThetaE = normalized_angle(ThetaE);
	    if ((nor_or==NEGATIVE_ORIENTATION) && (ThetaE>ThetaS))
		ThetaS += 2.0*PI;
	    if ((nor_or==POSITIVE_ORIENTATION) && (ThetaE<ThetaS))
		ThetaE += 2.0*PI;
	    coords_on_pert_ellipsoid(ellip->ThetaS,coords,ellip);
	    if (ellip->rbdry[0])
	        nearest_boundary_point(coords,coords,gr);
	    ns = make_node(Point(coords));
	    if (ellip->rbdry[0])
	    {
	        set_is_bdry(ns);
		if (ellip->wv_type < FIRST_PHYSICS_WAVE_TYPE)
		{
		    /* Project to corner if necessary */
		    for (i = 0; i < dim; ++i)
		    {
		        if (fabs(Coords(ns->posn)[i] - L[i]) < space*h[i])
			    Coords(ns->posn)[i] = L[i];
		        if (fabs(Coords(ns->posn)[i] - U[i]) < space*h[i])
			    Coords(ns->posn)[i] = U[i];
		    }
		}
	    }
	    coords_on_pert_ellipsoid(ellip->ThetaE,coords,ellip);
	    if (ellip->rbdry[1] == YES)
	        nearest_boundary_point(coords,coords,gr);
	    ne = make_node(Point(coords));
	    if (ellip->rbdry[1])
	    {
	        set_is_bdry(ne);
		if (ellip->wv_type < FIRST_PHYSICS_WAVE_TYPE)
		{
		    /* Project to corner if necessary */
		    for (i = 0; i < dim; ++i)
		    {
		        if (fabs(Coords(ne->posn)[i] - L[i]) < space*h[i])
			    Coords(ne->posn)[i] = L[i];
		        if (fabs(Coords(ne->posn)[i] - U[i]) < space*h[i])
			    Coords(ne->posn)[i] = U[i];
		    }
		}
	    }
	}
	cur = make_curve(l_comp,r_comp,ns,ne);
	ellip->hs = Hyper_surf(cur);
	if (ellip->untracked == YES)
	    untracked_hyper_surf(cur) = YES;

	if (ellip->n2o_enforced == YES)
	{
	    set_normal_to_origin(&hypersurface_normal_function(cur));
	    set_tangent_to_origin(&curve_tangent_function(cur));
	}

	done = NO;
	theta[0] = ThetaS;
	while (done == NO)
	{
	    done = increment_theta(theta,ThetaS,ThetaE,sgn,ellip,h,space);
	    if (done == YES)
		break;
	    coords_on_pert_ellipsoid(theta,coords,ellip);
	    if (insert_point_in_bond(Point(coords),cur->last,cur) !=
		FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in g_make_ellipse(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	}
	start_status(cur) = end_status(cur) = INCIDENT;
	wave_type(cur) = ellip->wv_type;
	surface_tension(cur) = ellip->surf_tension;
	total_mass(cur) = ellip->rgb_params.total_mass;
	mom_inertial(cur) = ellip->rgb_params.mom_of_inertial;
	angular_velo(cur) = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    center_of_mass(cur)[i] = ellip->rgb_params.center_of_mass[i];
	    center_of_mass_velo(cur)[i] = 0.0;
	}
	layer_index(ellip->hs) = ellip->layer_index;
	return ellip->hs;
}		/*end g_make_ellipse*/

EXPORT	void	g_make_ellip_region_boundaries2d(
	ELLIPSOID *ellip,
	Front     *front)
{
	CURVE     *c, *cin, *cout;
	COMPONENT ncomp_in[2], pcomp_in[2];
	COMPONENT ncomp, pcomp;
	POINT     *p;
	ELLIPSOID *inner;
	NODE      *ns, *ne;
	RECT_GRID *gr = front->rect_grid;
	double     v0[3], v1[3], vp[3];
	double     crds[3];
	int       i, k,  dim = gr->dim;
	static const ORIENTATION orient[2] = {POSITIVE_ORIENTATION,
					      NEGATIVE_ORIENTATION};

	if (dim != 2)
	{
	    screen("ERROR in g_make_ellip_region_boundaries2d(), "
		   "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	}

	if ((ellip->btype[0] == UNKNOWN_BOUNDARY_TYPE) &&
	    (ellip->btype[1] == UNKNOWN_BOUNDARY_TYPE))
	    return;

	cout = Curve_of_hs(ellip->hs);
	inner = inner_ellipsoid(ellip);
	cin = (inner != NULL) ? Curve_of_hs(inner->hs) : NULL;
	ns = ne = NULL;
	for (k = 0; k < 2; k++)
	{
	    pcomp_in[k] = ncomp_in[k] = NO_COMP;
	    if (ellip->btype[k] == UNKNOWN_BOUNDARY_TYPE)
		continue;
	    if (cin != NULL)
	        ns = Node_of(cin,orient[k]);
	    else if (ns == NULL)
	    {
		ns = make_node(Point(ellip->cen));
		node_type(ns) = FIXED_NODE;
	    }
	    ne = Node_of(cout,orient[k]);
	    p = Point_adjacent_to_node(cout,orient[k]);
	    for (i = 0; i < dim; i++)
	    {
		v0[i] = Coords(ne->posn)[i] - Coords(ns->posn)[i];
		v1[i] = Coords(p)[i] - Coords(ns->posn)[i];
	    }
	    (void) vector_product(v0,v1,vp,dim);
	    if (vp[0] < 0.0)
	    {
		pcomp_in[k] = ellip->compin;
		ncomp_in[k] = ellip->bcomp[k];
	    }
	    else
	    {
		ncomp_in[k] = ellip->compin;
		pcomp_in[k] = ellip->bcomp[k];
	    }
	    c = make_curve(ncomp_in[k],pcomp_in[k],ns,ne);
	    wave_type(c) = ellip->btype[k];
	    (*front->identify_physical_node)(ne);
	    (*front->identify_physical_node)(ns);
	}
	for (k = 0; k < 2; k++)
	{
	    if (ellip->obtype[k] == UNKNOWN_BOUNDARY_TYPE)
		continue;
	    ns = Node_of(cout,orient[k]);
	    intersect_ray_with_boundary(Coords(ns->posn),ellip->odir[k],
					gr->L,gr->U,crds,dim);
	    ne = make_node(Point(crds));
	    p = Point_adjacent_to_node(cout,orient[k]);
	    for (i = 0; i < dim; i++)
	    {
		v0[i] = Coords(ns->posn)[i] - crds[i];
		v1[i] = Coords(p)[i] - crds[i];
	    }
	    (void) vector_product(v0,v1,vp,dim);
	    if (vp[0] < 0.0)
	    {
		ncomp = ellip->compout;
		pcomp = ellip->obcomp[k];
	    }
	    else
	    {
		pcomp = ellip->compout;
		ncomp = ellip->obcomp[k];
	    }
	    c = make_curve(ncomp,pcomp,ns,ne);
	    wave_type(c) = ellip->obtype[k];
	    (*front->identify_physical_node)(ne);
	    (*front->identify_physical_node)(ns);
	}
}		/*end g_make_ellip_region_boundaries2d*/

/*
*			increment_theta():
*
*	Computes the incremented theta in the polar coordinated
*	of the z = constant cross section of an ellipsoid.
*	The difference new_theta - theta is determined so that
*	the displacement of the two consecutive points on the
*	ellipsoid are approximately separated by the distance
*	space in the scaled metric with scale vector h.
*	This function assumes that if sgn and theta_e - theta_s
*	have the same algebraic sign.
*/

LOCAL	boolean increment_theta(
	double		*theta,
	double		theta_s,
	double		theta_e,
	double		sgn,
	ELLIPSOID	*ellip,
	double		*h,
	double		space)
{
	double		d_theta, dl;
	double		new_theta;
	double		rad[MAXD], r;
	double		er[2];
	boolean		done;
	
	rad[0] = ellip->rad[0];	rad[1] = ellip->rad[1];
	dl = space*hypot(h[0],h[1]);
	er[0] = cos(*theta);	er[1] = sin(*theta);
	r = 1.0/hypot(er[0]/rad[0],er[1]/rad[1]);
	d_theta = dl/(r*r*r*hypot(er[0]/sqr(rad[0]),er[1]/sqr(rad[1])));

	new_theta = theta[0] + sgn*d_theta;
	done = NO;
	if (!Between(new_theta,theta_s,theta_e))
	{
	    new_theta = theta_e;
	    done = YES;
	}

	*theta = new_theta;
	return done;
}		/*end increment_theta*/

EXPORT	void coords_on_pert_ellipsoid(
	double		*theta,
	double		*coords,
	ELLIPSOID	*ellip)
{
	double		*cen = ellip->cen, *rad = ellip->rad, **Q = ellip->Q;
	FOURIER_POLY	*fpoly = ellip->fpoly;
	LEGENDRE_POLY	*lpoly = ellip->lpoly;
	int		dim = ellip->dim;
	double		Dr[MAXD];
	double		er[3];
	double		r;
	double           scale = ellip->scale;
	int		i, j;

	switch (dim)
	{
	case 2:
	    er[0] = cos(theta[0]);
	    er[1] = sin(theta[0]);
	    r = 1.0/hypot(er[0]/rad[0],er[1]/rad[1]);
	    break;
#if defined(THREED)
	case 3:
	    er[0] = cos(theta[0])*sin(theta[1]);
	    er[1] = sin(theta[0])*sin(theta[1]);
	    er[2] = cos(theta[1]);
	    r = 1.0/sqrt(sqr(er[0]/rad[0])+sqr(er[1]/rad[1])+sqr(er[2]/rad[2]));
	    break;
#endif /* defined(THREED) */
	default:
	    r = ERROR_FLOAT;
	    screen("ERROR in coords_on_pert_ellipsoid(), "
	           "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	}
	if (fpoly != NULL)
	    r += fourier_poly(theta,fpoly);
	if ((dim == 2) && (lpoly != NULL))
	    r += legendre_poly(er[1],lpoly);
	for (i = 0; i < dim; i++)
	    Dr[i] = scale*r*er[i];
	if (Q != NULL)
	{
	    for (i = 0; i < dim; i++)
	    {
	        coords[i] = cen[i];
	        for (j = 0; j < dim; j++)
	    	    coords[i] += Q[i][j]*Dr[j];
	    }
	}
	else
	{
	    for (i = 0; i < dim; i++)
		coords[i] = cen[i] + Dr[i];
	}
}		/*end coords_on_pert_ellipsoid*/

/*
*			g_make_fourier_curve():
*
*	g_make_fourier_curve() constructs a Fourier polynomial curve with 
*	Fourier data given in the structure fpoly.
*/


EXPORT void g_make_fourier_curve(
	int		w_type,
	int		num_points,
	double		x0,
	double		x1,
	FOURIER_POLY	*fpoly,
	COMPONENT	l_comp,
	COMPONENT	r_comp,
	double		surf_ten)
{
	CURVE		*cur;

	cur = i_make_fourier_curve(num_points,x0,x1,fpoly,l_comp,r_comp);

	wave_type(cur) = w_type;
	start_status(cur) = end_status(cur) = INCIDENT;
	surface_tension(cur) = surf_ten;
}		/*end g_make_fourier_curve*/



/*
*				make_ramp():
*
*	In the NORMAL geometry case, creates a ramp starting on the
*	lower wall, ramping at the given angle1 up to the given thickness,
*	and then ramping at the given angle2 to the right/top/bottom wall.
*	We also allow creation of a ramp with no corner by taking angle1
*	and angle2 the same.  A thickness of zero is suggested for this
*	case.
*	In the DOUBLE_REFL case, we are creating a double ramp.  This
*	geometry is basically the same as NORMAL with a positive angle1,
*	but angle2 measures how far below the x-axis the lower half of
*	the ramp reaches.  The corner here is at dist_from_left in x, and
*	the y coord is centered in the computational domain.
*
*	WARNING:  since we allow the ramp to attach to any of three walls,
*	care must be taken when declaring boundary types.
*/

EXPORT void make_ramp(
	int		geometry,
	double		angle1,
	double		angle2,
	double		thickness,
	double		dist_from_left,
	COMPONENT	l_comp,
	COMPONENT	r_comp,
	RECT_GRID	*rect_grid,
	CURVE		**ramp)
{
	CURVE		**curves;
	NODE		*ns, *ne;
	double		*L = rect_grid->L;
	double		*U = rect_grid->U;
	double		*h = rect_grid->h;
	double		start[MAXD];
	double		corner[MAXD];
	double		end[MAXD];
	boolean		sav_scss = interpolate_states_at_split_curve_node();

	if (geometry == NORMAL)
	{
	    end[0] = L[0] + dist_from_left;
	    end[1] = L[1];

	    corner[0] = end[0] + thickness/tan(angle1);
	    corner[1] = end[1] + thickness;

	    if (corner[0] >= U[0])
	    {
	        screen("ERROR in make_ramp(), "
	               "ramp too long - corner past right boundary\n");
	        clean_up(ERROR);
	    }

	    if (corner[1] >= U[1])
	    {
	    	screen("ERROR in make_ramp(), "
	    	       "ramp too thick - corner above top boundary\n");
	    	clean_up(ERROR);
	    }

	    if (angle2 >= 0.0)
	    {
	    	start[0] = U[0];
	    	start[1] = corner[1] + (start[0] - corner[0])*tan(angle2);

		if (start[1] > U[1] - .5*h[1])
		{
		    start[1] =  U[1];
		    start[0] = corner[0] + (start[1]-corner[1])/tan(angle2);
		}
	    }
	    else
	    {
	    	start[0] = corner[0] + thickness / tan(angle2);
	    	if (start[0] > U[0])
	    	{
	    	    start[0] = U[0];
		    start[1] = corner[1] + tan(angle2)*(start[0]-U[0]);
		}
		else
		    start[1] = L[1];
	    }
	}
	else if (geometry == DOUBLE_REFL)
	{
	    corner[0] = L[0] + dist_from_left;
	    corner[1] = thickness;

	    start[0] = ((U[1] - corner[1])/tan(angle1)) + corner[0];
	    if (start[0] > U[0])
	    {
	    	start[0] = U[0];
	    	start[1] = corner[1] + ((U[0] - corner[0])*tan(angle1));
	    }
	    else
	    	start[1] = U[1];

	    end[0] = (L[1] + corner[1])/tan(angle2) + corner[0];
	    if (end[0] > U[0])
	    {
	    	end[0] = U[0];
	    	end[1] = corner[1] - ((U[0] - corner[0])*tan(angle2));
	    }
	    else
	    	end[1] = L[1];
	}

	ns = make_node(Point(start));
	node_type(ns) = FIXED_NODE;
	ne = make_node(Point(end));
	node_type(ne) = FIXED_NODE;
	*ramp = make_curve(l_comp,r_comp,ns,ne);
	wave_type(*ramp) = NEUMANN_BOUNDARY;
	start_status(*ramp) = FIXED;
	end_status(*ramp) = FIXED;

	if ((angle2 != angle1) || (geometry == DOUBLE_REFL))
	{
	    	/* Insert node in ramp corner */

	    if ((U[0] - corner[0] > .0001*h[0]))		/*TOLERANCE*/
	    {
	    	if (insert_point_in_bond(Point(corner),(*ramp)->last,*ramp) !=
		    FUNCTION_SUCCEEDED)
	        {
	            screen("ERROR in insert_point_adjacent_to_node(), "
		           "insert_point_in_bond() failed\n");
	            clean_up(ERROR);
	        }
	    }
	    set_interpolate_states_at_split_curve_node(NO);
	    curves = split_curve((*ramp)->first->end,(*ramp)->first,*ramp,
				 l_comp,r_comp,l_comp,r_comp);
	    set_interpolate_states_at_split_curve_node(sav_scss);
	    node_type(curves[0]->end) = FIXED_NODE;
	    start_status(curves[1]) = end_status(curves[0]) = FIXED;
	    *ramp = curves[1];
	}
}		/*end make_ramp*/


/*
*			sin_sqr_pert():
*			sin_sqr_pert_prime():
*
*	Functions used for a purely geometric sinusoidal perturbation
*	of a curve.
*/

EXPORT	double	sin_sqr_pert(
	double		x,
	POINTER		params)
{
	double		nu = ((SIN_SQR_PERT_PARAMS *) params)->nu;
	double		epsilon = ((SIN_SQR_PERT_PARAMS *) params)->epsilon;

	return 0.5*epsilon*(1. - cos(2.*PI*nu*x));
}		/*end sin_sqr_pert*/

EXPORT	double	sin_sqr_pert_prime(
	double		x,
	POINTER		params)
{
	double		nu = ((SIN_SQR_PERT_PARAMS *) params)->nu;
	double		epsilon = ((SIN_SQR_PERT_PARAMS *) params)->epsilon;

	return PI*epsilon*nu*sin(2.*PI*nu*x);
}		/*end sin_sqr_pert_prime*/



/*
*			make_neumann_curve():
*
*	Makes a NEUMANN curve from node_s to node_e.  If either or both
*	are NULL it uses the respective points xs, ys, and/or xe, ye to
*	create a FIXED node of is_bdry_* boundary type.
*	The start_ and end_status of the curve must also be input.
*/

EXPORT void make_neumann_curve(
	INTERFACE	*intfc,
	NODE		*node_s,
	NODE		*node_e,
	double		xs,
	double		ys,
	double		xe,
	double		ye,
	int		is_bdry_ns,
	int		is_bdry_ne,
	COMPONENT	lcomp,
	COMPONENT	rcomp,
	CURVE		**crve)
{
	CURVE		*neumann_curve;
	NODE		*ns, *ne;
	double		coords[MAXD];

	if (exists_interface(intfc) != YES)
	{
	    screen("ERROR in make_neumann_curve(), no interface\n");
	    clean_up(ERROR);
	}

	if( node_s != NULL )
	    ns = node_s;
	else
	{
	    coords[0] = xs;	coords[1] = ys;
	    ns = make_node(Point(coords));
	    node_type(ns) = FIXED_NODE;
	}
	if( node_e != NULL )
	    ne = node_e;
	else
	{
	    coords[0] = xe;	coords[1] = ye;
	    ne = make_node(Point(coords));
	    node_type(ne) = FIXED_NODE;
	}
	if (is_bdry_ns)
	    set_is_bdry(ns);
	else
	    set_not_bdry(ns);

	if (is_bdry_ne)
	    set_is_bdry(ne);
	else
	    set_not_bdry(ne);
		
	neumann_curve = make_curve(lcomp,rcomp,ns,ne);
	*crve = neumann_curve;
	wave_type(neumann_curve) = NEUMANN_BOUNDARY;
	start_status(neumann_curve) = FIXED;
	end_status(neumann_curve)   = FIXED;
}		/*end make_neumann_curve*/


/*
*		make_vertical_axis_curve():
*
*	This function creates the unused "pencil" needed near the axis in 
*	the cylindrical	geometry.  More specifically, it creates FIXED NODEs 
*	at (rl_eff,YL) and (rl_eff,YU) and then creates NEUMANN curves
*	connecting all the nodes along the line r = rl_eff.  The left component
*	is compobst and the right component is obtained from the interior.
*
*	NOTE: If curves already run along the axis, duplicates will be created.
*	      If nodes already exist at (rl_eff,YL) or (rl_eff,YU), duplicates
*	      will be created.  Curves crossing the line x = rl_eff are
*	      truncated.  
*/

EXPORT void make_vertical_axis_curve(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	COMPONENT	comp_int;
	CROSS		*cross;
	CURVE		*axis_curve, *curves1[2], *curves2[2];
	CURVE		*phys_curve = NULL, **tempc;
	INTERFACE	*intfc = front->interf;
	NODE		*ns,*ne,*newn,*tempn;
	RECT_GRID	*gr = front->rect_grid;
	double		eps;		/* tolerance */
        double		coords[MAXD],min_sc_sep;
	double		rl_eff, YL, YU;
	ORIENTATION	orient = ORIENTATION_NOT_SET;

	debug_print("make_axis","Entered make_vertical_axis_curve()\n");
	if (exists_interface(intfc) != YES) 
	{
	    screen("ERROR in make_vertical_axis_curve(), - no interface\n");
	    clean_up(ERROR);
	}

	if (debugging("make_axis"))
	{
	    (void) printf("\n\nINTERFACE at start:\n");
	    print_interface(intfc);
	}

	rl_eff = RL_eff(init);
        if (rl_eff <= gr->GL[0])
	    return;
	
	set_obstacle_comp_type(comp_type(COMPOBST),front);
	YL = gr->GL[1];	YU = gr->GU[1];

	eps = 0.01*min(gr->h[0],gr->h[1]);/*TOLERANCE*/

	/* Make nodes for axis curve at top and bottom of domain */
	coords[0] = rl_eff;	coords[1] = YL;
	ns = make_node(Point(coords));
	coords[0] = rl_eff;	coords[1] = YU;
	ne = make_node(Point(coords));
	node_type(ns) = node_type(ne) = FIXED_NODE;
	    
	/* Make axis curve using a guess at the correct component */
	coords[0] = rl_eff + eps;
	coords[1] = YL + eps;
	comp_int = long_component(coords,intfc);

	rect_boundary_type(intfc,0,0) = PASSIVE_BOUNDARY;
	axis_curve = make_curve(COMPOBST,comp_int,ns,ne);
	wave_type(axis_curve) = NEUMANN_BOUNDARY;
	start_status(axis_curve) = end_status(axis_curve) = FIXED;

	/* Check for intersections with existing curves */

	if (intersections(intfc,&cross,NO) == FUNCTION_FAILED)
	{
	    screen("\nERROR in make_vertical_axis_curve(), "
	           "intersections() failed\n");
	}

	/* 
	*  Loop over the crosses (if there are any) splitting axis
	*  and crossed curves.  
	*/

	min_sc_sep = rl_eff+eps; /* TOLERANCE */

	while (cross)
	{
	    split_curves_at_cross(cross,front,&newn,curves1,NULL,NULL,curves2,
	    	                  NULL,NULL,min_sc_sep,NULL);

	    /* Need to get rid of curves in pencil.  Also need to
	    *  find the curve coming in from the interior since
	    *  we need to know its component numbers
	    */

	    for (tempc = newn->in_curves; tempc && *tempc; tempc++)
	    {
	    	tempn = (*tempc)->start;
	    	if (Coords(tempn->posn)[0] < rl_eff - eps)
	    	{
	    	    /* Inside pencil.  Delete curve and node */
	    	    (void) delete_curve(*tempc);
	    	    if (delete_node(tempn) == FUNCTION_FAILED)
	    	    {
	    	        /*
	    	         * Must be other curves connected to
	    	         * this node.  Can't deal with it.
	    	         */
	    	        screen("ERROR in make_vertical_axis_curve(), "
	    	               "delete_node() failed");
	    	        clean_up(ERROR);
	    	    }
	    	}
	    	else if (Coords(tempn->posn)[0] > rl_eff + eps)
	    	{
	    	    /* A curve from the interior */

	    	    phys_curve = *tempc;
	    	    orient = NEGATIVE_ORIENTATION;
	    	}
	    }

	    for (tempc = newn->out_curves; tempc && *tempc; tempc++)
            {
	    	tempn = (*tempc)->end;
                if (Coords(tempn->posn)[0] < rl_eff - eps)
                { 
                    /* Inside pencil.  Delete curve and node */ 
                    (void) delete_curve(*tempc); 
                                 
                    if (delete_node(tempn) == FUNCTION_FAILED)
                    { 
	    		/*
			 * Must be other curves connected 
	    		 * to this node.  Can't deal with it.
			 */ 
	    		screen("ERROR in make_vertical_axis_curve()"
	    		      "delete_node() failed");
	    		clean_up(ERROR); 
	    	    }
                }                                    
                else if (Coords(tempn->posn)[0] > rl_eff + eps)  
                { 
                    /* A curve from the interior */ 
                 
                    phys_curve = *tempc;                 
	    	    orient = POSITIVE_ORIENTATION;
                }                                 
            }        

	    if (phys_curve == NULL)
	    {
	    	screen("ERROR in make_vertical_axis_curve()"
	               "phys_curve is NULL\n");
	    	clean_up(ERROR); 
	    }

	    node_type(newn) = (wave_type(phys_curve)<FIRST_PHYSICS_WAVE_TYPE) ?
	    	              FIXED_NODE : NEUMANN_NODE;

	    /* 
	    *  Now we need to reset the component numbers along
	    *  the axis curve since our first guess was probably
	    *  wrong.  The components of the incoming axis curve
	    *  should be OK.  It's just the outward one we need
	    *  to worry about.
	    */

	    if (orient == NEGATIVE_ORIENTATION)
	    {
	    	/* 
	    	*  Since the physical curve is inward pointing,
	    	*  the axis curve is the only outward pointing
	    	*  curve.
	    	*/

	    	positive_component(newn->out_curves[0]) =
	    		positive_component(phys_curve);
	    }
	    else if (orient == POSITIVE_ORIENTATION)
	    {
	    	/*
	    	*  Both curves are outward pointing.  Which is the
	    	*  axis curve?
	    	*/


	    	if (newn->out_curves[0] == phys_curve)
	    	{
	    	    positive_component(newn->out_curves[1]) = 
	    	    	negative_component(phys_curve);
	    	}
	    	else
	    	{
	    	    positive_component(newn->out_curves[0]) =
	    	    	negative_component(phys_curve);
	    	}
	    }
	    else
	    {
                screen("ERROR in make_vertical_axis_curve(), "
                       "invalid orientation\n");
                clean_up(ERROR);
	    }

	    cross = cross->next;
	}

	if (debugging("make_axis"))
	{
	    (void) printf("\n\nINTERFACE at end:\n");
	    print_interface(intfc);
	}
	debug_print("make_axis","Left make_vertical_axis_curve()\n");

}		/*end make_vertical_axis_curve*/


EXPORT	void	set_normal_to_origin(
	NORMAL_FUNCTION *nf)
{
	nf->_normal = normal_to_origin;
	nf->_normal_name = strdup("normal_to_origin");
}		/*end set_normal_to_origin*/

LOCAL	void	normal_to_origin(
	POINT              *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF         *hs,
	double              *nor,
	Front              *front)
{
	int   i, dim = hs->interface->dim;
	double mag;

	if (dim == 1)
	{
	    nor[0] = 1.0;
	    return;
	}
	mag = mag_vector(Coords(p),dim);
	if (mag > 0.0)
	{
	    for (i = 0; i < dim; i++)
		nor[i] = Coords(p)[i]/mag;
#if defined(TWOD)
	    if (dim == 2)
	    {
		BOND *b = Bond_of_hse(hse);
		double n[2];

		n[0] =  (Coords(b->end)[1] - Coords(b->start)[1]);
		n[1] = -(Coords(b->end)[0] - Coords(b->start)[0]);
		if ((n[0]*nor[0] + n[1]*nor[1]) < 0.0)
		{
		    nor[0] *= -1.0;
		    nor[1] *= -1.0;
		}
	    }
#endif /* defined(TWOD) */
#if defined(THREED)
	    if (dim == 3)
	    {
		TRI *tri = Tri_of_hse(hse);
		const double *n = Tri_normal(tri);
		if (Dot3d(n,nor) < 0.0)
		{
		    nor[0] *= -1.0;
		    nor[1] *= -1.0;
		    nor[2] *= -1.0;
		}
	    }
#endif /* defined(THREED) */
	}
	else
	    (*interface_normal(hs->interface))(p,hse,hs,nor,front);
}		/*end normal_to_origin*/

EXPORT	void	set_tangent_to_origin(
	TANGENT_FUNCTION *tf)
{
	tf->_tangent = tangent_to_origin;
	tf->_tangent_name = "tangent_to_origin";
}		/*end set_tangent_to_origin*/

LOCAL	void	tangent_to_origin(
	POINT *p,
	BOND  *b,
	CURVE *c,
	double *tgnt,
	Front *front)
{
	INTERFACE *intfc = c->interface;
	double     sp;
	int       i, dim = intfc->dim;
	double     nor[3], mag;

	(*interface_tangent(intfc))(p,b,c,tgnt,front);
	mag = mag_vector(Coords(p),dim);
	if (mag > 0.0)
	{
	    for (i = 0; i < dim; i++)
		nor[i] = Coords(p)[i]/mag;
	    sp = scalar_product(tgnt,nor,dim);
	    for (i = 0; i < dim; i++)
		tgnt[i] -= sp*nor[i];
	    mag = mag_vector(tgnt,dim);
	    if (mag > 0.0)
	    {
	        for (i = 0; i < dim; i++)
		    tgnt[i] /= mag;
	    }
	    else
	        (*interface_tangent(intfc))(p,b,c,tgnt,front);
	}
}		/*end tangent_to_origin*/

#endif /* defined(TWOD) */
