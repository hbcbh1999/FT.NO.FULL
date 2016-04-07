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
*				gipolar.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains routines which determine the intersections of
*	shock polars describing shocks and rarefaction waves in two dimensions.
*
*	Basic routines:
*
* 		i_polar_2()
*
*	See the documentation at the top of gpolar.c.
*
*	Because there is an intersection of two or more waves under
*	consideration, each shock defines a ray or half line, starting from
*	the intersection point and extending to infinity.  The shock angle is
*	then defined to be the angle from the positive x axis to this ray.
*/

#include <gdecs/gdecs.h>

		/* structures for function parameters in gipolar.c */

typedef struct	{
	double		inc_shock_speed,
			inc_ang,
			Mach_ang,
			*refl_ang,
			*contact_ang;
	ANGLE_DIRECTION	i_to_a_dir;
	Locstate	ahead, behind,
			refl_bow, refl_mach;
} IP2_PARAMS;

typedef struct {
	Locstate	st0, st1;
	double		M0sq, M1sq;
	double		sign0, sign1;
} DT_ANG_PARAMS;

	/* LOCAL Function Declarations */
LOCAL	boolean	diff_turn_ang(double,double*,POINTER);
LOCAL	boolean	ip2(double,double*,POINTER);
LOCAL	int	i_polar_1(Locstate,Locstate,Locstate,Locstate,ANGLE_DIRECTION,
			  double*,double*,double*,double*);
LOCAL	int	spolar_intersect_bounds(double,double,double,Locstate,Locstate,
					double*,double*,double*,double*,
					int,int,RIEMANN_SOLVER_WAVE_TYPE*,
					RIEMANN_SOLVER_WAVE_TYPE*);


/*
*			i_polar_1():
*	This routines determines the intersection point of two shock polars
*	SP_0 and SP_1 in the pressure turning angle plane and solves for the 
*	states at the intersection point. There are two states, refl_bow and 
*	refl_mach at the intersection point, and these states are separated
*	by a contact. The bottom of the polar SP_1 is assumed to
*	lie on SP_0, and there may be zero, one or two other intersection
*	points. This configuration of the polars is that of a Mach triple 
*	point, so this language is used to describe the states and the shocks.
*	Input:
*		ahead, behind
*		abs_v, 		(defines the change to the rest frame)
*		i_to_a_dir	(CLOCKWISE or not inc shock to ahead dir)
*	Output:
*		refl_bow, refl_mach	(these differ only by a contact)
*		refl_ang	(from pos x axis to indicated shock)
*		contact_ang	(		"		   )
*		mach_ang	(		"		   )
*/


LOCAL int i_polar_1(
	Locstate	ahead,
	Locstate	behind,
	Locstate	refl_bow,
	Locstate	refl_mach,
	ANGLE_DIRECTION	i_to_a_dir,
	double		*abs_v,
	double		*refl_ang,
	double		*cont_ang,
	double		*mach_ang)
{
	double		p3;
	double		theta12;		/*turn angle across refl*/
	double		theta03;		/*turn angle across mach*/
	double		vel0_ang;		/*in steady frame*/
	double		vel1_ang;		/* "    "     "  */
	double		pmax, pmin;
	double		M0_sq, M1_sq;
	double		rel_v[MAXD];
	boolean		Cplus_w0;		/*wave 0 is mach stem*/
	boolean		Cplus_w1;		/*wave 1 is refl wave*/
	RIEMANN_SOLVER_WAVE_TYPE wtype0;
	RIEMANN_SOLVER_WAVE_TYPE wtype1;
	int		i, dim = Params(ahead)->dim;
	static Locstate st_0 = NULL, st_1 = NULL;
	static Locstate	st_2 = NULL, st_3 = NULL;
	static boolean	first = YES;

	debug_print("ipolar1","Entered i_polar_1()\n");

	if (first)
	{
		size_t	sizest = Params(ahead)->sizest;

		first = NO;
		(*Params(ahead)->_alloc_state)(&st_0,sizest);
		(*Params(ahead)->_alloc_state)(&st_1,sizest);
		(*Params(ahead)->_alloc_state)(&st_2,sizest);
		(*Params(ahead)->_alloc_state)(&st_3,sizest);
	}

	set_state(st_0,TGAS_STATE,ahead);
	set_state(st_1,TGAS_STATE,behind);

	if (debugging("ipolar1"))
	{
		verbose_print_state("ahead",ahead);
		verbose_print_state("behind",behind);
		(void) printf("\tabs_v[0] = (%g, %g)\n",abs_v[0],abs_v[1]);
		print_angle_direction("\ti_to_a_dir =",i_to_a_dir,"\n");
	}

	M0_sq = mach_number_squared(st_0,abs_v,rel_v);
	vel0_ang = angle(rel_v[0],rel_v[1]);
	if (debugging("ipolar1"))
	{
		print_angle("\tvel0 angle =",vel0_ang,"\n");
		(void) printf("\tahead Mach number = %g\n",
			      mag_vector(rel_v,dim)/sound_speed(st_0));
	}
	if (mag_vector(rel_v,dim) < sound_speed(st_0))
	{
		if (debugging("ipolar1"))
		{
			(void) printf("WARNING in i_polar_1(), ");
			(void) printf("ahead is subsonic in the rest frame.\n");
		}
		debug_print("ipolar1","Left i_polar_1()\n");
		return NO;
	}

	M1_sq = mach_number_squared(st_1,abs_v,rel_v);
	if (debugging("ipolar1"))
	{
		vel1_ang = angle(rel_v[0],rel_v[1]);
		print_angle("\tvel1 angle =",vel1_ang,"\n");
		(void) printf("\tst1 Mach number = %g\n",sqrt(M1_sq));
	}
	if (mag_vector(rel_v,dim) < sound_speed(st_1))
	{
		if (debugging("ipolar1"))
		{
		       (void) printf("WARNING in i_polar_1(), ");
		       (void) printf("behind is subsonic in the rest frame.\n");
		}
		debug_print("ipolar1","Left i_polar_1()\n");
		return NO;
	}

	if (i_to_a_dir == CLOCKWISE)
	{
		Cplus_w0 = YES;
		Cplus_w1 = NO;
	}
	else
	{
		Cplus_w0 = NO;
		Cplus_w1 = YES;
	}
	
	if (pr_at_max_turn_angle(&pmin,M0_sq,st_0) == FUNCTION_FAILED)
	{
		if (debugging("ipolar1"))
		{
			(void) printf("WARNING in i_polar_1(), ");
			(void) printf("pr_at_max_turn_angle() failed\n");
		}
		return NO;
	}
	pmax = max_behind_shock_pr(M1_sq,st_1);
	if (!intersection_of_two_shock_polars(st_0,st_1,abs_v,&p3,
	                                         &pmin,&pmax,Cplus_w0,
						 Cplus_w1,&wtype0,&wtype1))
	{
		if (debugging("ipolar1"))
		{
			(void) printf("WARNING in i_polar_1(), ");
			(void) printf("unable to compute refl_mach pressure\n");
		}
		debug_print("ipolar1","Left i_polar_1()\n");
		return NO;
	}

	if ((wtype0 != SHOCK) || (wtype1 != SHOCK))
	{
		if (debugging("ipolar1"))
		{
			(void) printf("WARNING in i_polar_1() -- ");
			(void) printf("non-shock wave types.\n");
		}
		debug_print("ipolar1","Left i_polar_1()\n");
		return NO;
	}

	if (!s_polar_3(st_1,YES,p3,Cplus_w1,YES,abs_v,st_2,
			  refl_ang,&theta12))
	{
		if (debugging("ipolar1"))
		{
			(void) printf("WARNING in i_polar_1() -- ");
			(void) printf("unable to compute refl_bow.\n");
		}
		debug_print("ipolar1","Left i_polar_1()\n");
		return NO;
	}

	if (!s_polar_3(st_0,YES,p3,Cplus_w0,YES,abs_v,st_3,
				mach_ang,&theta03))
	{
		if (debugging("ipolar1"))
		{
			(void) printf("WARNING in i_polar_1() -- ");
			(void) printf("unable to compute refl_mach.\n");
		}
		debug_print("ipolar1","Left i_polar_1()\n");
		return NO;
	}

	*cont_ang = normalized_angle(vel0_ang + theta03);

	set_state(refl_bow,GAS_STATE,st_2);
	set_state(refl_mach,GAS_STATE,st_3);

	if (debugging("ipolar1"))
	{
		double ang;
		(void) printf("\n");
		print_angle("\trefl angle = %g d\n",*refl_ang,"\n");
		print_angle("\tcont angle = %g d\n",*cont_ang,"\n");
		print_angle("\tMach angle = %g d\n",*mach_ang,"\n");
		print_angle("\ttheta12 = %g d\n",theta12,"\n");
		print_angle("\ttheta03 = %g d\n",theta03,"\n");
		for (i = 0; i < dim; i++)
			rel_v[i] = vel(i,st_2) - abs_v[i];
		ang = angle(rel_v[0],rel_v[1]);
		print_angle("\tvel2_ang =",ang,"\n");
		for (i = 0; i < dim; i++)
			rel_v[i] = vel(i,st_3) - abs_v[i];
		ang = angle(rel_v[0],rel_v[1]);
		print_angle("\tvel3_ang =",ang,"\n");
		ang = (vel1_ang+theta12) - (vel0_ang+theta03);
		print_angle("(vel1_ang+theta12) - (vel0_ang+theta03) =",
			    ang,"\n");
		verbose_print_state("refl_bow",refl_bow);
		verbose_print_state("refl_mach",refl_mach);
	}

	debug_print("ipolar1","Left i_polar_1()\n");
	return YES;
}		/*end i_polar_1*/

	
/*     
*			i_polar_2():
*	This routine determines the same intersection point of shock polars
*	as does i_polar_1, but on the basis of different input data.  In
*	particular the transformation to the steady frame is not given here,
*	while the Mach angle is given.
*	Input:
*		ahead, behind
*		inc_ang		(pos x axis to incident shock)
*		Mach_ang	(pos x axis to Mach stem)
*		i_to_a_dir	(CLOCKWISE or not for inc shock to ahead dir)
*	Output:
*		refl_bow, refl_mach	(these differ only by a contact)
*		abs_v,		(defines the change to the rest frame)
*		refl_ang	(from pos x axis to indicated shock)
*		contact_ang	(		"		  )
*
*	The rest frame angle of the node trajectory must lie in the
*	90 degree sector bounded by the incident shock and the
*	forward facing  normal to the incident shock. This information 
*	is needed by the root solving subroutine.
*	TODO: Find a bound which separates the node trajectory from the
*	incident shock itself.
*/


EXPORT int i_polar_2(
	Locstate	ahead,
	Locstate	behind,
	Locstate	refl_bow,
	Locstate	refl_mach,
	ANGLE_DIRECTION	i_to_a_dir,
	double		*abs_v,
	double		aw_ang,
	double		inc_ang,
	double		*refl_ang,
	double		*contact_ang,
	double		Mach_ang)
{
	IP2_PARAMS	ip2_params;
	double		inc_shock_speed;
	double		node_speed;
	double		node_ang, ans;
	double		min_ang, max_ang;
	double		i_n[SMAXD];	/* inc shock normal points to ahead */
	double		delta, epsilon;
	const double     meps = MACH_EPS;

	debug_print("ipolar2","Entered i_polar_2()\n");
#if defined(DEBUG_GIPOLAR)
	if (debugging("ipolar2"))
	{
		verbose_print_state("ahead",ahead);
		verbose_print_state("behind",behind);
		print_angle("inc_ang =", inc_ang,"\n");
		print_angle("Mach_ang =",Mach_ang,"\n");
		print_angle_direction("\ti_to_a_dir =",i_to_a_dir,"\n");
	}
#endif /* defined(DEBUG_GIPOLAR) */

		/* Find incident shock speed */

	if (i_to_a_dir == COUNTER_CLOCK)
	{
		i_n[0] = cos(inc_ang + PI / 2.);
		i_n[1] = sin(inc_ang + PI / 2.);
		min_ang = inc_ang;
		max_ang = aw_ang;
	}
	else
	{
		i_n[0] = cos(inc_ang - PI / 2.);
		i_n[1] = sin(inc_ang - PI / 2.);
		min_ang = aw_ang;
		max_ang = inc_ang;
	}
	(void) s_polar_4(BEHIND_PRESSURE,pressure(behind),
			 &inc_shock_speed,i_n,ahead,behind,GAS_STATE);

#if defined(DEBUG_GIPOLAR)
	if (debugging("ipolar2"))
	{
		print_angle("min_ang =", min_ang,"\n");
		print_angle("max_ang =",max_ang,"\n");
		(void) printf("inc_shock_speed = %g\n",inc_shock_speed);
	}
#endif /* defined(DEBUG_GIPOLAR) */

		/* Initialize ip2_params */

	ip2_params.ahead = ahead;
	ip2_params.behind = behind;
	ip2_params.refl_bow = refl_bow;
	ip2_params.refl_mach = refl_mach;
	ip2_params.inc_shock_speed = inc_shock_speed;
	ip2_params.inc_ang = inc_ang;
	ip2_params.refl_ang = refl_ang;
	ip2_params.contact_ang = contact_ang;
	ip2_params.Mach_ang = Mach_ang;
	ip2_params.i_to_a_dir = i_to_a_dir;

		/* Solve for node_angle */

	if (min_ang > max_ang) min_ang -= 2. * PI;
	delta = EPS*fabs(max_ang-min_ang);
	delta = max(delta, meps);
	(void) ip2(0.5*(max_ang+min_ang), &epsilon, (POINTER) &ip2_params);
	epsilon = fabs(epsilon)*EPS;
	epsilon = max(epsilon, meps);
	if (debugging("ip2_vals"))
	{
		print_function_values(ip2,(POINTER) &ip2_params,
				      0.0,min_ang,max_ang,100,"ip2",stdout);
	}
	if (find_root(ip2,(POINTER)&ip2_params,0.0,&node_ang,
			min_ang,max_ang,epsilon,delta) == FUNCTION_FAILED)
		return NO;

	/* Since ip2() is a discontinuous function, the root finder may find
	   an invalid root at the point of discontinuity.  Check for this.
	*/
	(void) ip2(node_ang,&ans,(POINTER) &ip2_params);
	if (fabs(ans) > epsilon)
	{
		if (debugging("ipolar2"))
		{
			screen("WARNING in i_polar_2(), ");
			screen("root_finder returned invalid root\n");
		}
		return NO;
	}

		/* Determine node velocity */

	node_speed = inc_shock_speed/fabs(sin(inc_ang - node_ang));
	abs_v[0] = node_speed * cos(node_ang);
	abs_v[1] = node_speed * sin(node_ang);

	debug_print("ipolar2","Left i_polar_2()\n");
	return YES;
}		/*end i_polar_2*/


/*
*		intersection_of_two_shock_polars():
*
*	This function finds the pressure at the point of intersection
*	of two shock polars in the pressure turning angle phase space.
*	The shock polar for st1 includes that part of the rarefaction
*	curve which might give a generic intersection with the st0
*	shock polar.
*
*	The input is:
*			st0, st1		* Given states *
*			Cplus_w0, Cplus_w1	* Wave Families * 
*			abs_v, 	* Transformation to the steady frame *
*			pmin, pmax		* Optional limits on 
*
*	The output is:
*			p		* pressure at intersection *
*			wtype0, wtype1	* Wave types (SHOCK or RAREFACTION) *
*
*	If the pointer pmin or pmax is non NULL it is assumed that the
*	intersection pressure *p satisfies *pmin <= *p or *p <= *pmax,
*	respectively.
*
*	Returns YES if intersection is found, NO otherwise.
*	The values of wtype0 and wtype1 are set even if no intersection
*	is found.  These correspond to the only possible wave types
*	at the intersection point and may be used for bifurcation
*	analysis.  If either wtype0 or wtype1 is UNKNOWN_WAVE_TYPE,  then
*	NO intersection was found and the only segment of possible
*	intersection occured on the subsonic portion of shock polar
*	one or two respectively.
*/

#define set_parameters(parameters,Cplus_w0,Cplus_w1,st0,st1,M0sq,M1sq)	\
	(parameters).sign0 = (Cplus_w0) ? 1.0 : -1.0;			\
	(parameters).sign1 = (Cplus_w1) ? 1.0 : -1.0;			\
	(parameters).st0 = (st0);					\
	(parameters).st1 = (st1);					\
	(parameters).M0sq = (M0sq);					\
	(parameters).M1sq = (M1sq);

EXPORT int intersection_of_two_shock_polars(
	Locstate	         st0,
	Locstate	         st1,
	double		         *abs_v,
	double		         *p,
	double		         *pmin,
	double		         *pmax,
	boolean		         Cplus_w0,
	boolean		         Cplus_w1,
	RIEMANN_SOLVER_WAVE_TYPE *wtype0,
	RIEMANN_SOLVER_WAVE_TYPE *wtype1)
{
	double		v0[MAXD];	/* Velocities in the steady frame */
	double		v1[MAXD];
	double		M0sq, M1sq;	/* Mach #'s squared in steady frame */
	double		theta01;	/* turn angle vector v0 to v1 */
	double		v0_x_v1;
	double		v0_d_v1;
	double		p0, p1;
	double		p_upper, p_lower;
	const double     meps = 10.0*MACH_EPS;/*TOLERANCE*/
	double		eps_theta;
	double		eps_pressure;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double		min_pressure;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int		dim;
	DT_ANG_PARAMS	parameters;
	static char	yesstatus[] =
		"Left intersection_of_two_shock_polars() status = YES\n";
	static	char	nostatus[]  =
		"Left intersection_of_two_shock_polars() status = NO\n";
	
	debug_print("shock_polar","Entered intersection_of_two_shock_polars()\n");

	dim = Params(st0)->dim;
	p0 = pressure(st0);
	p1 = pressure(st1);
	M0sq = mach_number_squared(st0,abs_v,v0);
	M1sq = mach_number_squared(st1,abs_v,v1);

	(void) vector_product(v0,v1,&v0_x_v1,dim);
	v0_d_v1 = scalar_product(v0,v1,dim);
	theta01 = atan2(v0_x_v1,v0_d_v1);
#if defined(DEBUG_GIPOLAR)
	if (debugging("shock_polar")) 
	{
		(void) printf("abs_v = (%g, %g)\n",abs_v[0],abs_v[1]);
		(void) printf("Cplus_w0 = %s, Cplus_w1 = %s\n",
			      (Cplus_w0) ? "YES" : "NO",
			      (Cplus_w1) ? "YES" : "NO");
		verbose_print_state("st0",st0);
		(void) printf("v0 = (%g, %g), q0 = %g, M0 = %g\n",v0[0],v0[1],
			      mag_vector(v0,dim),sqrt(M0sq));
		verbose_print_state("st1",st1);
		(void) printf("v1 = (%g, %g), q1 = %g, M1 = %g\n\n",v1[0],v1[1],
			      mag_vector(v1,dim),sqrt(M1sq));
		print_angle("theta01 =",theta01,"\n");
	}
#endif /* defined(DEBUG_GIPOLAR) */

	*wtype0 = *wtype1 = UNSET_RIEMANN_SOLVER_WAVE_TYPE;
	if (!spolar_intersect_bounds(theta01,M0sq,M1sq,st0,
					st1,pmin,pmax,&p_upper,&p_lower,
					Cplus_w0,Cplus_w1,wtype0,wtype1))
	{
	    if (debugging("shock_polar"))
	    {
		(void) printf("WARNING in intersection_of_two_shock_polars(), "
		              "spolar_intersect_bounds() failed.\n");
		debug_print("shock_polar",nostatus);
	    }
	    return NO;
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	min_pressure = max(Min_pressure(st0),Min_pressure(st1));
	if (p_upper < min_pressure)
	{
	    *p = min_pressure;
	    *wtype0 = *wtype1 = RAREFACTION;
	    debug_print("shock_polar",yesstatus);
	    return YES;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	if (p_upper <= p0)	*wtype0 = RAREFACTION;
	if (p_lower >= p0)	*wtype0 = SHOCK;
	if (p_upper <= p1)	*wtype1 = RAREFACTION;
	if (p_lower >= p1)	*wtype1 = SHOCK;
	eps_theta = fabs(theta01)*EPS;
	eps_theta = max(meps,eps_theta);
	eps_pressure = max(p0,p1)*EPS;
	eps_pressure = max(meps,eps_pressure);
	set_parameters(parameters,Cplus_w0,Cplus_w1,st0,st1,M0sq,M1sq);
	if (find_root(diff_turn_ang,(POINTER)&parameters,theta01,p,p_lower,
		      p_upper,eps_theta,eps_pressure) == FUNCTION_FAILED)
	{
	    if (debugging("shock_polar"))
	    {
		(void) printf("WARNING in intersection_of_two_shock_polars()");
		(void) printf(", No intersection of shock polars.\n");
	    }
	    debug_print("shock_polar",nostatus);
	    return NO;
	}
	*wtype0 = (*p >= p0) ? SHOCK : RAREFACTION;
	*wtype1 = (*p >= p1) ? SHOCK : RAREFACTION;
	debug_print("shock_polar",yesstatus);
	return YES;
}		/*end intersection_of_two_shock_polars*/

#define sib_failed(mesg)						\
{									\
	if (debugging("spibnds"))					\
	{								\
	    (void) printf("WARNING in spolar_intersect_bounds(), ");	\
	    (void) printf("%s() failed.\n",mesg);			\
	}								\
	debug_print("spibnds",nostatus);					\
}

#define sswc(p,M0sq,theta,st)						\
if (steady_state_wave_curve((p),(M0sq),(theta),(st)) == FUNCTION_FAILED)\
{									\
	sib_failed("steady_state_wave_curve");				\
	return NO;							\
}

#if defined(DEBUG_GIPOLAR)
#define print_point_on_sswc(vname,p,t0,t1,t01)				\
{									\
	(void) printf("%s = %g, ",(vname),(p));				\
	print_angle("theta0 =",(t0),", ");				\
	print_angle("theta1 =",(t1),"\n");				\
	print_angle("\ttheta0 - (theta01 - theta1) =",			\
		    (t0) - ((t01) - (t1)),"\n");			\
}
#endif /* defined(DEBUG_GIPOLAR) */

LOCAL int spolar_intersect_bounds(
	double		          theta01,
	double		          M0sq,
	double		          M1sq,
	Locstate	          st0,
	Locstate	          st1,
	double		          *pmin,
	double		          *pmax,
	double		          *p_upper,
	double		          *p_lower,
	int		          Cplus_w0,
	int		          Cplus_w1,
	RIEMANN_SOLVER_WAVE_TYPE *wtype0,
	RIEMANN_SOLVER_WAVE_TYPE *wtype1)
{
	double		pr_max;
	double		p0, p1;
	double		p0_max, p1_max;
	double		pmid, pclosest;
	double		dtheta_closest;
	double		theta0, theta1;
	double		press[8];
	double		ptmp[8];
	double		x[3], y[3];
	double		adiff, adiff_last;
	const double     meps = 10.0*MACH_EPS;/*TOLERANCE*/
	double		pr_epsilon;
	int		i, j, jmin, n;
	DT_ANG_PARAMS	parameters;
	int		status;
	static const char *yesstatus =
				"Left spolar_intersect_bounds() status = YES";
	static const char *nostatus  =
				"Left spolar_intersect_bounds() status = NO";
	
	debug_print("spibnds","Entered spolar_intersect_bounds()\n");

	if (Cplus_w0 == Cplus_w1)
	{
	    screen("ERROR in spolar_intersect_bounds(), "
	           "nongeneric ssrp, Cplus_w0 == Cplus_w1\n");
	    clean_up(ERROR);
	}

	if (!Cplus_w0)
	{
	    status = spolar_intersect_bounds(-theta01,M1sq,M0sq,st1,
				             st0,pmin,pmax,p_upper,p_lower,
					     Cplus_w1,Cplus_w0,wtype1,wtype0);
	    debug_print("spibnds","%s\n",(status) ? yesstatus : nostatus);
	    return status;
	}


	p0 = pressure(st0);	p1 = pressure(st1);
	if (p1 < p0)
	{
	    status = spolar_intersect_bounds(theta01,M1sq,M0sq,st1,
				             st0,pmin,pmax,p_upper,p_lower,
				             !Cplus_w1,!Cplus_w0,
					     wtype1,wtype0);
	    debug_print("spibnds","%s\n",(status) ? yesstatus : nostatus);
	    return status;
	}

	sswc(0.0,M0sq,&theta0,st0);
	sswc(0.0,M1sq,&theta1,st1);
	adiff_last = theta0 - (theta01 - theta1);
	if (adiff_last > 0.0)
	{
	    /* VACUUM AT INTERSECTION, RAREFACTION - RAREFACTION */
	    *p_upper = *p_lower = 0.0;

	    if (debugging("spibnds"))
	    {
	    	(void) printf("Vacuum at intersection, "
	    	              "rarefaction - rarefaction\n");
	    	(void) printf("\tp_lower = %g, p_upper = %g\n",
	    		      *p_lower,*p_upper);
	    }
		debug_print("spibnds","%s\n",yesstatus);
		return YES;
	}

	/* Check for roots between base pressures, sonic points and
		max turn angle pressures */

	p0_max = max_behind_shock_pr(M0sq,st0);
	p1_max = max_behind_shock_pr(M1sq,st1);
	press[0] = 0.0;
	press[1] = p0;
	x[0] = pressure_at_sonic_point(M0sq,st0);
	if (pr_at_max_turn_angle(&x[1],M0sq,st0) == FUNCTION_FAILED)
	{
	    sib_failed("pr_at_max_turn_angle");
	    return NO;
	}
	y[0] = p1;
	y[1] = pressure_at_sonic_point(M1sq,st1);
	if (pr_at_max_turn_angle(&y[2],M1sq,st1) == FUNCTION_FAILED)
	{
	    sib_failed("pr_at_max_turn_angle");
	    return NO;
	}
	pr_max = min(p0_max,p1_max);

	/* Sort above pressures into an increasing array */

	n = 2;
	jmin = 0;
	for (i = 0; i < 3; i++)
	{
	    for (j = jmin; j < 2; j++)
	    {
	    	if (x[j] > y[i])
		    break;
	    	press[n++] = x[j];
	    	jmin = j+1;
	    }
	    press[n++] = y[i];
	}
	for (j = jmin; j < 2; j++)
	{
	    press[n++] = x[j];
	}
	for (n = 2; n < 7; n++)
	{
	    if (press[n] > pr_max)
	    {
	    	n--;
	    	break;
	    }
	}
#if defined(DEBUG_GIPOLAR)
	if (debugging("spibnds")) 
	{
	    (void) printf("State 0 data\n");
	    (void) printf("\tpressure state 0 =\t\t\t%g\n",p0);
	    (void) printf("\tsonic point pressure state 0 =\t\t%g\n",x[0]);
	    (void) printf("\tpressure at max turn angle state 0 =\t%g\n",x[1]);
	    (void) printf("\tmax behind shock pressure state 0 =\t%g\n",p0_max);
	    (void) printf("State 1 data\n");
	    (void) printf("\tpressure state 1 =\t\t\t%g\n",p1);
	    (void) printf("\tsonic point pressure state 1 =\t\t%g\n",y[1]);
	    (void) printf("\tpressure at max turn angle state 1 =\t%g\n",y[2]);
	    (void) printf("\tmax behind shock pressure state 1 =\t%g\n",p1_max);
	    (void) printf("n = %d\n",n);
	    for (i = 0; i < 7; i++)
	        (void) printf("\tpress[%d] = %g\n",i,press[i]);
	    if (pmax != NULL)
	        (void) printf("pmax = %g\n",*pmax);
	    else
	        (void) printf("pmax = NULL\n");
	    if (pmin != NULL)
	        (void) printf("pmin = %g\n",*pmin);
	    else
		(void) printf("pmin = NULL\n");
	}
#endif /* defined(DEBUG_GIPOLAR) */

	/* Sort pmax and pmin (if given) into the pressure array */

	if (pmax != NULL)
	    pr_max = min(pr_max,*pmax);
	if ((pmin != NULL && *pmin > press[n-1]))
	{
	    if (*pmin > pr_max)
	    {
	        if (debugging("spibnds"))
	        {
	    	    (void) printf("WARNING in spolar_intersect_bounds(), "
		                  "no possible root above pmin = %g\n",*pmin);
		}
		return NO;
	    }
	    set_parameters(parameters,Cplus_w0,Cplus_w1,st0,st1,M0sq,M1sq);

	    pr_epsilon = max(meps, EPS*pr_max);

	    if (find_separation_point(diff_turn_ang,
	    			      (POINTER)&parameters,theta01,&pmid,
				      *pmin,pr_max,&pclosest,&dtheta_closest,
				      pr_epsilon) == FUNCTION_FAILED)
	    {
		sib_failed("find_separation_point");
#if defined(DEBUG_GIPOLAR)
		if (debugging("spibnds")) 
		{
		    sswc(pclosest,M0sq,&theta0,st0);
		    sswc(pclosest,M1sq,&theta1,st1);
		    print_point_on_sswc("pclosest",pclosest,
					theta0,theta1,theta01);
		    sswc(pr_max,M0sq,&theta0,st0);
		    sswc(pr_max,M1sq,&theta1,st1);
		    print_point_on_sswc("pr_max",pr_max,
					theta0,theta1,theta01);
		}
#endif /* defined(DEBUG_GIPOLAR) */
		return NO;
	    }
	    *p_upper = pmid;
	    *p_lower = *pmin;
#if defined(DEBUG_GIPOLAR)
	    if (debugging("spibnds"))
	    {
	    	sswc(pmid,M0sq,&theta0,st0);
	    	sswc(pmid,M1sq,&theta1,st1);
	    	print_point_on_sswc("pmid",pmid,theta0,theta1,theta01);
		(void) printf("p_lower = %g, p_upper = %g\n",*p_lower,*p_upper);
	    }
#endif /* defined(DEBUG_GIPOLAR) */
	    debug_print("spibnds","%s\n",yesstatus);
	    return YES;
	}
	else if (pmin != NULL)
	{
	    for (i = 0; i < n; i++)
	    	if (press[i] > *pmin)
		    break;
	    ptmp[0] = *pmin;
	    for (j = 0; j < n-i; j++)
	    	ptmp[j+1] = press[i+j];
	    n -= i - 1;
	    for (i = 0; i < n; i++)
	    	press[i] = ptmp[i];
	    sswc(press[0],M0sq,&theta0,st0);
	    sswc(press[0],M1sq,&theta1,st1);
	    adiff_last = theta0 - (theta01 - theta1);
	}
	if (pmax != NULL && *pmax < press[n-1])
	{
	    for (i = n - 1; i >= 0; i--)
	    	if (press[i] < *pmax)
		    break;
	    n = i + 2;
	    press[i+1] = *pmax;
	}

	/* Check for interval with root */

	for (i = 1; i < n; i++)
	{
	    sswc(press[i],M0sq,&theta0,st0);
	    sswc(press[i],M1sq,&theta1,st1);
#if defined(DEBUG_GIPOLAR)
	    if (debugging("spibnds")) 
	    {
	    	char vname[20];
	    	(void) sprintf(vname,"press[%d]",i);
	    	print_point_on_sswc(vname,press[i],theta0,theta1,theta01);
	    }
#endif /* defined(DEBUG_GIPOLAR) */
	    adiff = theta0 - (theta01 - theta1);
	    if (adiff*adiff_last <= 0.0)
	    {
	    	*p_upper = press[i];
	    	*p_lower = press[i-1];
#if defined(DEBUG_GIPOLAR)
	    	if (debugging("spibnds"))
	    	{
	    	    (void) printf("p_lower = %g, p_upper = %g\n",
	    			  *p_lower,*p_upper);
	    	}
#endif /* defined(DEBUG_GIPOLAR) */
	    	debug_print("spibnds","%s\n",yesstatus);
	    	return YES;
	    }
	    adiff_last = adiff;
	}
	set_parameters(parameters,Cplus_w0,Cplus_w1,st0,st1,M0sq,M1sq);

	pr_epsilon = max(meps, pr_max*EPS);

	if (find_separation_point(diff_turn_ang,
				  (POINTER)&parameters,theta01,&pmid,
				  press[n-1],pr_max,&pclosest,&dtheta_closest,
				  pr_epsilon) == FUNCTION_FAILED)
	{
	    sib_failed("find_separation_point");
#if defined(DEBUG_GIPOLAR)
	    if (debugging("spibnds")) 
	    {
	    	sswc(pclosest,M0sq,&theta0,st0);
	    	sswc(pclosest,M1sq,&theta1,st1);
	    	print_point_on_sswc("pclosest",pclosest,theta0,theta1,theta01);
		sswc(pr_max,M0sq,&theta0,st0);
		sswc(pr_max,M1sq,&theta1,st1);
		print_point_on_sswc("pr_max",pr_max,theta0,theta1,theta01);
	    }
#endif /* defined(DEBUG_GIPOLAR) */
	    if (press[n-1] >= p0)
		*wtype0 = SHOCK;
	    if (p1 > p0_max)
		*wtype1 = RAREFACTION;
	    return NO;
	}
	*p_upper = pmid;
	*p_lower = press[n-1];
#if defined(DEBUG_GIPOLAR)
	if (debugging("spibnds"))
	{
	    sswc(pmid,M0sq,&theta0,st0);
	    sswc(pmid,M1sq,&theta1,st1);
	    print_point_on_sswc("pmid",pmid,theta0,theta1,theta01);
	    (void) printf("p_lower = %g, p_upper = %g\n",*p_lower,*p_upper);
	}
#endif /* defined(DEBUG_GIPOLAR) */
	debug_print("spibnds","%s\n",yesstatus);
	return YES;
}		/*end spolar_intersect_bounds*/


/*
*				ip2():
*	All data in this routine is in the computational frame, ie the
*	frame of the data of i_polar_2.
*/



LOCAL boolean   ip2(
	double u,		/* trial node angle in computational frame */
	double *diff,
	POINTER parameters) 
{
	static const double SMALL_ANG         = PI/30.0; /*TOLERANCE*/
	static const double twopi_m_small_ang = 2.0*PI - PI/30.0;/*TOLERANCE*/
	static const double twopi             = 2.0*PI;
	static const double halfpi            = 0.5*PI;

	double	node_speed;
	double	trial_Mach_ang;
	double	nod_ang_to_inc_ang;
	double	abs_v[SMAXD];		/*node velocity*/
	double	Msq;			/*steady Mach # sq - temp*/
	double	pinc_max;		/*max pr behind inc shock*/
	double	prefl_max;		/*max pr behind refl shock*/
	double	inc_shock_speed	  = ((IP2_PARAMS *)parameters)->inc_shock_speed;
	double	inc_ang		  = ((IP2_PARAMS *)parameters)->inc_ang;
	double	Mach_ang 	  = ((IP2_PARAMS *)parameters)->Mach_ang;
	double	*refl_ang	  = ((IP2_PARAMS *)parameters)->refl_ang;
	double	*contact_ang	  = ((IP2_PARAMS *)parameters)->contact_ang;
	ANGLE_DIRECTION i_to_a_dir= ((IP2_PARAMS *)parameters)->i_to_a_dir;
	Locstate ahead		  = ((IP2_PARAMS *)parameters)->ahead;
	Locstate behind		  = ((IP2_PARAMS *)parameters)->behind;
	Locstate refl_bow	  = ((IP2_PARAMS *)parameters)->refl_bow;
	Locstate refl_mach	  = ((IP2_PARAMS *)parameters)->refl_mach;
	static boolean	first = YES;
	static double	sin_small_ang;

	debug_print("ip2","Entered ip2()\n");
	if (first == YES)
	{
	    first = NO;
	    sin_small_ang = sin(SMALL_ANG);
	}

	nod_ang_to_inc_ang = normalized_angle(inc_ang - u);
	if (fabs(nod_ang_to_inc_ang) < SMALL_ANG ||
	    fabs(nod_ang_to_inc_ang) > twopi_m_small_ang)
	{ /* A GUESS */
	    node_speed = inc_shock_speed / sin_small_ang;
	    abs_v[0] = node_speed * cos(u);
	    abs_v[1] = node_speed * sin(u);
	}
	else 
	{
	    node_speed = inc_shock_speed / fabs(sin(inc_ang - u));
	    abs_v[0] = node_speed * cos(u);
	    abs_v[1] = node_speed * sin(u);
	}

#if defined(DEBUG_GIPOLAR)
	if (debugging("ip2"))
	{
	    print_angle("\ttrial_node_ang =",u,"\n");
	    print_angle("\tinc_ang =",inc_ang,"\n");
	    (void) printf("\tinc_shock_speed = %g node_speed = %g\n",
	    	          inc_shock_speed,node_speed);
	    (void) printf("\tabs_v = (%g, %g)\n",abs_v[0],abs_v[1]);
	}
#endif /* defined(DEBUG_GIPOLAR) */

	if (!i_polar_1(ahead,behind,refl_bow,refl_mach,i_to_a_dir,
		abs_v,refl_ang,contact_ang,&trial_Mach_ang))
	{
	    Msq = mach_number_squared(ahead,abs_v,(double *)NULL);
	    pinc_max = max_behind_shock_pr(Msq,ahead);
	    Msq = mach_number_squared(behind,abs_v,(double *)NULL);
	    prefl_max = max_behind_shock_pr(Msq,behind);
	    if (((prefl_max > pinc_max) && (i_to_a_dir == CLOCKWISE))
	    	 ||
	        ((prefl_max < pinc_max) && (i_to_a_dir == COUNTER_CLOCK)))
	    {
	    	*diff = -halfpi;
	    }
	    else
	    {
	    	*diff = halfpi;
	    }
	}
	else 
	{
	    *diff = Mach_ang - trial_Mach_ang;
#if defined(DEBUG_GIPOLAR)
	    if (debugging("ip2"))
	    {
	    	print_angle("\ttrial_Mach_ang =",trial_Mach_ang,"\n");
	    	print_angle("\tdiff =",*diff,"\n");
	    }
#endif /* defined(DEBUG_GIPOLAR) */
	}
	*diff = (*diff < PI) ? *diff : *diff - twopi;
	*diff = (*diff > -PI) ? *diff : *diff + twopi;

	debug_print("ip2","Left ip2()\n");
	return FUNCTION_SUCCEEDED;
}		/*end ip2*/


LOCAL boolean diff_turn_ang(
	double		p,
	double		*ans,
	POINTER		parameters)
{
	Locstate	st0 = ((DT_ANG_PARAMS *) parameters)->st0;
	Locstate	st1 = ((DT_ANG_PARAMS *) parameters)->st1;
	double		M0sq = 	((DT_ANG_PARAMS *) parameters)->M0sq;
	double		M1sq = 	((DT_ANG_PARAMS *) parameters)->M1sq;
	double		sign0 = ((DT_ANG_PARAMS *) parameters)->sign0;
	double		sign1 = ((DT_ANG_PARAMS *) parameters)->sign1;
	double		theta0, theta1;

	debug_print("dt_ang","Entered diff_turn_ang()\n");
	debug_print("dt_ang","\tp = %g\n",p);
	
	if (steady_state_wave_curve(p,M0sq,&theta0,st0) == FUNCTION_FAILED)
	{
		screen("WARNING in diff_turn_ang(), ");
		screen("steady_state_wave_curve() failed\n");
		return FUNCTION_FAILED;
	}
	if (steady_state_wave_curve(p,M1sq,&theta1,st1) == FUNCTION_FAILED)
	{
		screen("WARNING in diff_turn_ang(), ");
		screen("steady_state_wave_curve() failed\n");
		return FUNCTION_FAILED;
	}
	
	*ans = sign0*theta0 - sign1*theta1;

	if (debugging("dt_ang"))
	{
		print_angle("\ttheta0 =",theta0,"\n");
		print_angle("\ttheta1 =",theta1,"\n");
		print_angle("\tans =",*ans,"\n");
	}
	debug_print("dt_ang","Left diff_turn_ang()\n");
	return FUNCTION_SUCCEEDED;
}		/*end diff_turn_ang*/
