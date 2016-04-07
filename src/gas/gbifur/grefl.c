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
*				girefl.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for initializing and analyzing shock
*	reflection or refraction problems.
*/


#include <gdecs/gdecs.h>
#if defined(TWOD)

	/* LOCAL Function Declarations */
LOCAL	int	find_bow_state_and_oblateness(int,double,double*,double,
					      Bubble*,Front*);
LOCAL	int	init_mach_corner_states(Bubble*,Front*);

/*
*			init_bubble():
*
*	This routine initializes the states for a regular reflection,
*	assumed to arise as a result of a normal shock running into a
*	wedge corner. It is assumed that the corner is stationary in
*	the computaional frame. This routine also sets a few parameters 
*	to fix the approximate shape of the reflection bubble, for 
*	initialization purposes.
*	Input:
*		ahead state
*		behind state
*		reflection state
*		reflected shock angle
*		node velocity
*		corner velocity
*		incident shock tangent, oriented away from refl point
*		ahead wall tangent, oriented away from corner
*		behind wall tangent, oriented away from corner
*		scale_param, one of:
*			time elapsed
*			bow wall length
*			reflection wall length
*	Output:
*		states at corner
*		state at base of bow shock
*		time elapsed
*		bow wall length
*		reflection wall length
*		is attached (for bow shock), if YES, also:
*			angle pos x axis to bow shock
*/

EXPORT void init_bubble(
	double		*node_v,
	int		which_scale_param,
	double		scale_param,
	Bubble		*bubble,
	Front		*front)
{
	RP_DATA		*RP = bubble->RP;
	Gas_param	*params = Params(RP->state[0]);
	Locstate	ahead = RP->state[0];
	Locstate	behind = RP->state[1];
	Locstate	refl = RP->state[2];
	Locstate	bow = bubble->bow;
	Locstate	refl_corner = bubble->refl_corner;
	Locstate	bow_corner = bubble->bow_corner;
	double		node_speed;
	double		*inc_t = bubble->inc_t;
	double		*aw_t = bubble->aw_t;
	double		*bw_t = bubble->bw_t;
	double		*cor_v = bubble->cor_v;
	double		refl_ang = RP->ang[2];
	int		dim = params->dim;
	size_t		sizest = params->sizest;

	debug_print("init_bubble","Entered init_bubble()\n");
	if (debugging("init_bubble"))
	{
		double	ang;

		print_general_vector("Corner velocity = ",cor_v,dim,"\n");
		print_general_vector("Ahead wall tangent = ",aw_t,dim,"\n");
		print_general_vector("Behind wall tangent = ",bw_t,dim,"\n");
		print_general_vector("Node velocity = ",node_v,dim,"\n");
		verbose_print_state("ahead",ahead);
		(void) printf("Mach number ahead = %g\n",
			      mach_number(ahead,node_v));
		verbose_print_state("behind",behind);
		(void) printf("Mach number behind = %g\n",
			      mach_number(behind,node_v));
		ang = angle(aw_t[0],aw_t[1]);
		print_angle("ahead wall angle =",ang,"\n");
	}

		/* Set oblateness and bow data */

	debug_print("init_bubble","Bow shock is detached\n");
	if (!find_bow_state_and_oblateness(which_scale_param,
		scale_param,node_v,0.0,bubble,front))
	{
		Error(ERROR,"ERROR: in init_bubble()\n");
		(void) printf("find_bow_state_and_oblateness() failed\n");
		clean_up(ERROR);
	}

		/* No corner states are set in is_attached case */
		/* Set corner state:  density from refl, entropy from bow */

	if (!bubble->is_attached)
	{
		node_speed = hypot(node_v[0] - cor_v[0],node_v[1] - cor_v[1]);
		state_on_adiabat_with_dens(bow,Dens(refl),
					   refl_corner,GAS_STATE);
		ft_assign(bow_corner,refl_corner,sizest);
	}

	if (debugging("init_bubble"))
	{
	    double ang;

	    (void) printf("init_bubble --\n");
	    (void) printf("node_speed = %g\n",node_speed);
	    ang = normalized_angle(angle(inc_t[1],inc_t[0]) -
				   angle(aw_t[0],aw_t[1]));
	    print_angle("incident shock angle relative to ahead wall =",
			ang,"\n");
	    ang = normalized_angle(refl_ang - angle(aw_t[0],aw_t[1]));
	    print_angle("reflected shock angle relative to ahead wall =",
			ang,"\n");
	    print_angle("refl_ang =",refl_ang,"\n");

	    verbose_print_state("ahead",ahead);
	    verbose_print_state("behind",behind);
	    verbose_print_state("reflection",refl);
	    verbose_print_state("bow",bow);
	    verbose_print_state("corner",refl_corner);
	}
	debug_print("init_bubble","Left init_bubble()\n");
}		/*end init_bubble*/

/*
*			init_mach_bubble():
*
*	This routine initializes the states for a Mach reflection
*	assumed to result either from a normal shock or a regular
*	reflection where the incident angle has shrunk/grown resp.
*	so that bifurcation to a Mach configuration takes place.
*	It is assumed that the corner is stationary in the frame
*	of the computation. This routine also determines a few
*	parameters to fix approximately the shape of the reflection
*	bubble, for initialization purposes.  The shape is used only
*	in the case of transition from a normal shock to a Mach
*	reflection.
*	Input:
*		ahead state
*		behind state
*		shock speed
*		shock tangent, oriented away from Mach triple point
*		ahead wall tangent, oriented away from corner
*		behind wall tangent, oriented away from corner
*		scale_param, one of:
*			time elapsed
*			bow wall length
*			reflection wall length
*	Output:
*		state at corner
*		state at base of bow shock
*		states base of Mach stem
*		reflected state - behind refl shock at triple point
*		Mach state - state behind Mach stem at triple point
*		bow wall length
*		reflection wall length
*		mach height
*		refl angle - pos x axis to reflected shock
*		contact angle - pos x axis to contact
*		Mach angle - pos x axis to Mach stem
*		is attached (for bow shock), if YES, also:
*			angle pos x axis to bow shock
*/

EXPORT int init_mach_bubble(
	int		which_scale_param,
	double		scale_param,
	Bubble		*bubble,
	int		flag,
	Front		*front)
{
	RP_DATA		*RP = bubble->RP;
	Gas_param	*params = Params(RP->state[0]);
	Locstate	corner;
	Locstate	ahead = RP->state[0];
	Locstate	behind = RP->state[1];
	Locstate	refl_bow = RP->state[2];
	Locstate	refl_mach = RP->state[3];
	Locstate	bow = bubble->bow;
	Locstate	mach = bubble->mach;
	Locstate	contact_mach = bubble->contact_mach;
	Locstate	contact_bow = bubble->contact_bow;
	double		node_speed;
	double		wall_mach_pr;
	double		aw_ang;			/* pos x axis to ahead wall */
	double		wall_to_node_ang;
	double		node_v[MAXD];
	double		*aw_t = bubble->aw_t;
	double		*inc_t = bubble->inc_t;
	double		inc_ang = RP->ang[0];
	double		*refl_ang = &RP->ang[1];
	double		*cont_ang = &RP->ang[2];
	double		*mach_ang = &RP->ang[3];
	double		*cor_v = bubble->cor_v;
	int		i, dim = params->dim;

	debug_print("init_mach","Entered init_mach_bubble()\n");

	aw_ang = angle(aw_t[0],aw_t[1]);
	*mach_ang = (RP->ang_dir == CLOCKWISE) ?
		normalized_angle(aw_ang + PI/2.0) :
			normalized_angle(aw_ang - PI/2.0);

	if (debugging("init_mach"))
	{
		double	ang;

		print_general_vector("Shock tangent, inc_t = ",inc_t,dim,"\n");
		print_general_vector("Ahead wall tangent = ",aw_t,dim,"\n");
		print_general_vector("Corner velocity = ",cor_v,dim,"\n");
		print_angle("inc_ang =",inc_ang,"\n");
		print_angle("aw_ang =",aw_ang,"\n");

		ang = acos(scalar_product(inc_t,aw_t,dim));
		print_angle("i_to_w_ang =",ang,"\n");
		print_angle("mach_ang =",*mach_ang,"\n");
		print_angle_direction("RP->ang_dir =",RP->ang_dir,"\n");
	}

	if (!i_polar_2(ahead,behind,refl_bow,refl_mach,
		Opposite_ang_dir(RP->ang_dir),node_v,aw_ang,inc_ang,
		refl_ang,cont_ang,*mach_ang))
	{
		(void) printf("WARNING in init_mach_bubble(), ");
		(void) printf("i_polar_2(), failed\n");
		return NO;
	}

	if (RP->ang_dir == CLOCKWISE)
		wall_to_node_ang =
			normalized_angle(aw_ang - angle(node_v[0],node_v[1]));
	else
		wall_to_node_ang =
			normalized_angle(angle(node_v[0],node_v[1]) - aw_ang);

	node_speed = hypot(node_v[0] - cor_v[0],node_v[1] - cor_v[1]);
	if (debugging("init_mach"))
	{
		print_angle("wall_to_node_ang =",wall_to_node_ang,"\n");
	}

	if (flag == NORMAL_TO_MACH_REFLECTION)
	{
		if (!find_bow_state_and_oblateness(which_scale_param,
			scale_param,node_v,wall_to_node_ang,bubble,front))
		{
		    screen("ERROR in init_mach_bubble(), ");
		    screen("find_bow_state_and_oblateness() failed\n");
		    clean_up(ERROR);
		}
	}
	else if (flag == REGULAR_TO_MACH_REFLECTION)
	{
		double	initial_time_elapsed;

		bubble->is_attached = NO;
		if (which_scale_param == TIME_ELAPSED)
		{
			initial_time_elapsed = scale_param;
			bubble->refl_length = node_speed *
				initial_time_elapsed * cos(wall_to_node_ang);
			bubble->mach_height = node_speed *
				initial_time_elapsed *
					fabs(sin(wall_to_node_ang));
		}
		else
		{
			screen("ERROR in init_mach_bubble(), ");
			screen("scale_param not equal to TIME_ELAPSED ");
			screen("for regular to mach reflection\n");
			return NO;
		}
		set_initial_time_elapsed(initial_time_elapsed);
	}

		/* Set mach state */

	bubble->mach_speed = node_speed * cos(wall_to_node_ang);
	if (!s_polar_4(SHOCK_SPEED,bubble->mach_speed,&bubble->mach_speed,
			  aw_t,ahead,mach,GAS_STATE))
	{
		(void) printf("WARNING in init_mach_bubble(), ");
		(void) printf("s_polar_4() failed\n");
		return NO;
	}

	/* Set contact states. We assume that the density and the entropy
	 * are nearly constant in the region between the mach stem and the
	 * contact.  Thus the pressure at the contact is pressure(mach).
	 * Furthermore we assume that the entropy has a jump at some point
	 * approximately half way between the contact and the corner;
	 * between this point and the bow the entropy is constant with the
	 * value at the bow, while between this point and the contact the
	 * entropy is constant with the value given at the Mach node between
	 * the reflected shock and the contact.
	 */

	if (RP->ang_dir == CLOCKWISE)
		bubble->contact_speed = bubble->mach_speed *
		    (1.0 - tan(wall_to_node_ang)*tan(*cont_ang - *mach_ang));
	else
		bubble->contact_speed = bubble->mach_speed *
		    (1.0 - tan(wall_to_node_ang)*tan(*mach_ang - *cont_ang));

	wall_mach_pr = pressure(mach);

	set_state(contact_mach,TGAS_STATE,mach);
	Press(contact_mach) = wall_mach_pr;
	for (i = 0; i < dim; i++)
		Vel(contact_mach)[i] =  bubble->contact_speed * aw_t[i];
	set_state(contact_mach,GAS_STATE,contact_mach);

	state_on_adiabat_with_pr(refl_bow,wall_mach_pr,contact_bow,GAS_STATE);
	for (i = 0; i < dim; i++)
		Mom(contact_bow)[i] =
			Dens(contact_bow) * bubble->contact_speed * aw_t[i];
	Energy(contact_bow) +=
		0.5 * Dens(contact_bow) * sqr(bubble->contact_speed);

	/* No corner states are set in is_attached case */

	if ((flag == NORMAL_TO_MACH_REFLECTION) && !bubble->is_attached)
	{
		/* density from contact, entropy from bow */

		corner = bubble->refl_corner;
		state_on_adiabat_with_dens(bow,Dens(contact_bow),
					   corner,GAS_STATE);
		ft_assign(bubble->bow_corner,corner,params->sizest);
	}
	else if (flag == REGULAR_TO_MACH_REFLECTION)
	{
	    if (bubble->is_node_at_corner)
	    {
		if (!init_mach_corner_states(bubble,front))
		{
		    if (debugging("init_mach"))
		    {
			screen("ERROR in init_mach_bubble, ");
			screen("unable to set corner states in ");
			screen("reg to mach refl bifurcation\n");
		    }
		    return NO;
		}
	    }
	}

	if (debugging("init_mach"))
	{
		double ang;

		(void) printf("init_mach_bubble --\n");
		ang = angle(inc_t[0],inc_t[1]) - aw_ang;
		print_angle("incident shock angle wrt to ahead wall =",
			    ang,"\n");
		ang = *refl_ang - aw_ang;
		print_angle("reflected shock angle wrt to ahead wall =",
			    ang,"\n");
		ang = *cont_ang - aw_ang;		
		print_angle("slip line angle wrt to ahead wall =",ang,"\n");
		ang = *mach_ang - aw_ang;
		print_angle("mach stem angle wrt to ahead wall =",ang,"\n");

		verbose_print_state("ahead",ahead);
		verbose_print_state("behind",behind);
		verbose_print_state("mach",mach);
		if (flag == NORMAL_TO_MACH_REFLECTION)
		{
			(void) printf("In NORMAL_TO_MACH_REFLECTION case.\n");
			verbose_print_state("bow",bow);
			verbose_print_state("corner",bubble->refl_corner);
		}
		else
		{	
			(void) printf("In REGULAR_TO_MACH_REFLECTION case.\n");
		}
		if (bubble->is_node_at_corner)
		{
			(void) printf("There is a node at the corner.\n");
			verbose_print_state("refl_corner",bubble->refl_corner);
			verbose_print_state("bow_corner",bubble->bow_corner);
		}
		else
		{
			(void) printf("There is no node at the corner.");
		}
		verbose_print_state("refl_mach",refl_mach);
		verbose_print_state("refl_bow",refl_bow);
		verbose_print_state("contact_mach",contact_mach);
		verbose_print_state("contact_bow",contact_bow);
	}
	debug_print("init_mach","Left init_mach_bubble()\n");
	return YES;
}		/*end init_mach_bubble*/


/*
*			init_mach_corner_states():
*
*	Support function for init_mach_bubble() in the case where we
*	are bifurcating from regular to mach reflection and there is
*	a sharp corner in the wall.  The bubble is presumably the bow
*	bubble.  The state refl_corner is on the same wall curve as the
*	base of the contact, and bow_corner is around the corner.
*	We copy the state at the base of the contact into refl_corner,
*	and then use w_speed to compute the state around the corner
*	using the interior bisector of the corner as normal.
*/

LOCAL int init_mach_corner_states(
	Bubble		*bubble,
	Front		*front)
{
	Locstate	refl_corner = bubble->refl_corner;
	Locstate	bow_corner = bubble->bow_corner;
	Locstate	contact = bubble->contact_bow;
	Locstate	sl, sr, ansl, ansr;
	double		a_ang, b_ang;
	double		*cor_v = bubble->cor_v;
	double		nor_ang;
	double		nor[MAXD];
	double		W[MAXD];
	int		i, dim;
	size_t		sizest = Params(bubble->RP->state[0])->sizest;
	static boolean	first = YES;

	if (first == YES)
	{
		first = NO;
		alloc_state(front->interf,&sr,sizest);
		alloc_state(front->interf,&ansr,sizest);
		alloc_state(front->interf,&sl,sizest);
		alloc_state(front->interf,&ansl,sizest);
	}

	a_ang = angle(bubble->aw_t[0],bubble->aw_t[1]);
	b_ang = angle(bubble->bw_t[0],bubble->bw_t[1]);
	ft_assign(refl_corner,contact,sizest);
	if (bubble->RP->ang_dir == CLOCKWISE)
	{
		nor_ang = normalized_angle(a_ang - 
			0.5 * normalized_angle(a_ang - b_ang));
		sl = refl_corner;
		ansl = bow_corner;
		sr = return_obst_state();
	}
	else
	{
		nor_ang = normalized_angle(a_ang + 
			0.5 * normalized_angle(b_ang - a_ang));
		sr = refl_corner;
		ansr = bow_corner;
		sl = return_obst_state();
	}
	nor[0] = cos(nor_ang);
	nor[1] = sin(nor_ang);
	dim = Params(refl_corner)->dim;
	for (i = 0; i < dim; i++) cor_v[i] *= -1.0;
	add_velocity_to_state(sl,cor_v);
	add_velocity_to_state(sr,cor_v);

	w_speed(bubble->corner_posn,sl,sr,ansl,ansr,W,0.0,nor,
		NEUMANN_BOUNDARY,front);

	for (i = 0; i < dim; i++) cor_v[i] *= -1.0;
	add_velocity_to_state(sl,cor_v);
	add_velocity_to_state(sr,cor_v);
	add_velocity_to_state(ansl,cor_v);
	add_velocity_to_state(ansr,cor_v);

	return YES;
}		/*end init_mach_corner_states*/


/*
*			find_bow_state_and_oblateness():
*
*	This functions calculates the the bow pressure and the oblateness
*	of the reflected shock bubble for regular and Mach reflection.
*	For regular reflection, wall_to_node_ang = 0.
*	The input data are:
*	
*	The output are:
*		bow state - state behind bow at wall
*		oblateness - bow_length / refl_length
*		refl_wall_length - corner to reflection point (reg)
*				   corner to base of mach stem (mach)
*		bow_wall_length - corner to base of bow shock
*
*	Returns YES if sucesssful.
*/

typedef struct {
	double u1, v1, p2, r1, theta1, r0, k, eps;
	double mf_eps;
	double k0, k1;
	Locstate state1;
	int st_type;
} BUBBLE_PARAMS;


/*
*		find_bow_state_and_oblateness();
*
*	This function computes the state behind the bow shock, and the
*	lengths of the refl and bow walls, the height of the mach stem (for
*	mach reflection), and the time elapsed since bifurcation.
*	One estimate for the bow length comes from using a straight bow
*	shock.  Another comes from solving a Riemann problem for
*	normal ahead wall.  The answer that should be correct comes from
*	solving a Riemann problem between the behind and bow states at
*	the reflection (RP->state[1] and RP->state[2]).
*	Note: this function may not be correct for CCW angle dir.
*/

LOCAL int find_bow_state_and_oblateness(
	int		which_scale_param,
	double		scale_param,
	double		*node_v,
	double		wall_to_node_ang,
	Bubble		*bubble,
	Front		*front)
{
	Locstate	behind = bubble->RP->state[1];
	Locstate	bow = bubble->bow;
	Locstate	refl_st = bubble->RP->state[2];
	double		*aw_t = bubble->aw_t;
	double		*bw_t = bubble->bw_t;
	double		rs_t[MAXD];
	double		refl_ang;
	double		W[MAXD];
	double		num, den;
	double		node_dir[MAXD];
	double		node_speed;
	double		aw_ang;
	double		*bow_length = &bubble->bow_length;
	double		*refl_length = &bubble->refl_length;
	double		*mach_height = &bubble->mach_height;
	double		oblateness;	/* bow_length / refl_length */
	double		initial_time_elapsed;
	double		bow_pr;
	double		fdummy;
	int		i, dim = Params(bubble->RP->state[0])->dim;
	int		w_type = FORWARD_SHOCK_WAVE;
	static Locstate sl = NULL, sr = NULL, ansl = NULL, ansr = NULL;
	static boolean	first = YES;
	static char	fname[] = "find_bow_state_and_oblateness()";

	debug_print("find_bow_st","Entered find_bow_state_and_oblateness()\n");
	start_clock("find_oblateness");

	if (first)
	{
		size_t	sizest = Params(bubble->RP->state[0])->sizest;

		first = NO;
		alloc_state(front->interf,&sl,sizest);
		alloc_state(front->interf,&sr,sizest);
		alloc_state(front->interf,&ansl,sizest);
		alloc_state(front->interf,&ansr,sizest);
	}

	aw_ang = angle(aw_t[0],aw_t[1]);
	bubble->is_attached = is_bow_shock_attached(behind,aw_ang,
					    bow,&bubble->bow_base_ang);

	if (bubble->is_attached)
	{
		oblateness = 0.0;
	}
	else
	{
		/*Note: angles labelled differently for reg vs Mach refl */
	    refl_ang = (wall_to_node_ang == 0.0) ?
			bubble->RP->ang[2] : bubble->RP->ang[1];
	    rs_t[0] = cos(refl_ang);
	    rs_t[1] = sin(refl_ang);

	    if (debugging("find_bow_st"))
	    {
		print_general_vector("Bow wall tangent = ",bw_t,dim,"\n");
		print_general_vector("Ahead wall tangent = ",aw_t,dim,"\n");
		print_general_vector("Reflected shock tangent = ",
				     rs_t,dim,"\n");
		print_general_vector("Node Velocity = ",node_v,dim,"\n");
		print_general_vector("Corner Velocity = ",
				     bubble->cor_v,dim,"\n");

		verbose_print_state("behind",behind);
		verbose_print_state("refl_st",refl_st);
	    }

	    for (i = 0; i < dim; i++)
		    node_dir[i] = node_v[i] - bubble->cor_v[i];
	    node_speed = mag_vector(node_dir,dim);
	    *refl_length = node_speed * cos(wall_to_node_ang);
	    for (i = 0; i < dim; i++)
		    node_dir[i] /= node_speed;

	    if (debugging("find_bow_st"))
	    {
		    (void) printf("Reflected wall length = %g\n",*refl_length);
	    }
	    if (*refl_length < 0.0)
	    {
		    (void) printf("WARNING in %s, ",fname);
		    (void) printf("Invalid orientation\n");
		    debug_print("find_bow_st",
			 "Leaving find_bow_state_and_oblateness(), ans = NO\n");
		    return NO;
	    }

		/* use law of sines to compute bow length for straight bow */
	    (void) vector_product(bw_t,rs_t,&den,dim);
	    (void) vector_product(node_dir,rs_t,&num,dim);
	    *bow_length = *refl_length * num / fabs(den);
	    if (debugging("find_bow_st"))
		    (void) printf("Bow len computed from straight shock = %g\n",
				  *bow_length);

		/* Compute wave speed for normal ahead wall */

	    set_state(sr,GAS_STATE,behind);
	    Mom(sr)[0] = -Mom(behind)[0];
	    Mom(sr)[1] = -Mom(behind)[1];
	    w_speed(bubble->corner_posn,behind,sr,ansl,ansr,W,0.0,
		    bw_t,w_type,front);
	    *bow_length = min(hypot(W[0],W[1]),*bow_length);
	    bow_pr = pressure(ansr);
	    if (debugging("find_bow_st"))
	    {
		(void) printf("Answer from behind - ahead Riemann problem, ");
		(void) printf("bow_length = %g, bow_pr = %g\n",W[0],bow_pr);
	    }

		/* Compute wave speed using refl_st behind */

	    w_speed(bubble->corner_posn,behind,refl_st,ansl,ansr,W,0.0,
		    bw_t,w_type,front);
	    if (debugging("find_bow_st"))
	    {
		(void) printf("Answer from behind - refl_st Riemann problem, ");
		(void) printf("bow_length = %g, bow_pr = %g\n",
			      W[0],pressure(ansr));
	    }
	    if (hypot(W[0],W[1]) < *bow_length)
	    {
		    *bow_length = hypot(W[0],W[1]);
		    bow_pr = pressure(ansr);
	    }
	    oblateness = *bow_length / *refl_length;
	    stop_clock("find_oblateness");

	    if (oblateness < 0.0)
	    {
		    screen("ERROR in %s, ",fname);
		    screen("corner angle less than 90 degrees");
		    clean_up(ERROR);
	    }
	    if (debugging("find_bow_st"))
	    {
		    (void) printf("Oblateness = %g\n",oblateness);
		    (void) printf("Bow pressure = %g\n",bow_pr);
	    }

	    (void) s_polar_4(BEHIND_PRESSURE,bow_pr,&fdummy,bw_t,
			     behind,bow,GAS_STATE);
        }

		/* Set scale length and time elapsed since bifurcation */

	switch(which_scale_param)
	{
	case BOW_LENGTH:
		debug_print("find_bow_st","which_scale_param = BOW_LENGTH = %g\n",
			scale_param);
		*bow_length = scale_param;
		*refl_length = *bow_length / oblateness;
		initial_time_elapsed = *refl_length /
			cos(wall_to_node_ang) / node_speed;
		break;
	case REFL_WALL_LENGTH:
		debug_print("find_bow_st",
			"which_scale_param = REFL_WALL_LENGTH = %g\n",
				scale_param);
		*refl_length = scale_param;
		*bow_length = *refl_length * oblateness;
		initial_time_elapsed = *refl_length /
			cos(wall_to_node_ang) / node_speed;
		break;
	case TIME_ELAPSED:
		debug_print("find_bow_st","which_scale_param = TIME_ELAPSED = %g\n",
			scale_param);
		/*Note: initial_time_elapsed is the time since bifurcation*/
		initial_time_elapsed = scale_param;
		*refl_length = node_speed * initial_time_elapsed *
			cos(wall_to_node_ang);
		*bow_length = *refl_length * oblateness;
		break;
	default:
		Error(ERROR,"ERROR in find_bow_state_and_oblateness()");
		screen("unknown scale param type (%d) %s\n",
			which_scale_param,"in find_bow_state_and_oblateness()");
		return NO;
	}
	*mach_height = fabs(*refl_length * tan(wall_to_node_ang));
	set_initial_time_elapsed(initial_time_elapsed);

	debug_print("find_bow_st","Left find_bow_state_and_oblateness()\n");
	return YES;
}		/*end find_bow_state_and_oblateness*/

#define pressure_at_angle_on_bubble(theta,p2,k)				\
		((p2)*(1.0 + (k)*(theta1)))

#define k_from_pressure_theta(theta,p,p2)				\
		((p)/(p2) - 1.0)/(theta)


/*
*			find_bubble_state():
*
*	find_bubble_state() computes the state at the point (x,y) at the
*	interior of a bubble as a linear combination of states at nodes
*	and points on the bubble interface .
*/

#define pt_frc(x)	(1.0/x)

EXPORT void find_bubble_state(
	Locstate	state,
	double		*coords,
	COMPONENT	comp,
	INTERFACE	*intfc,
	int		stype)
{
	HYPER_SURF	*hs;
	HYPER_SURF_ELEMENT *hse;
	BOND		*bs,*b,*bond;
	CURVE		*curve;
	RECT_GRID	*gr = computational_grid(intfc);
	Locstate	s1, s2, st,st1; 
	Locstate	comb_of_states;
	O_CURVE		*oc;
	O_CURVE_FAMILY	*loop;
	POINT		P;
	double		coords_on[MAXD];
	double		t, alpha, beta;
	double		tot_pt_wt, pt_wt, d1, crp;
	double		d[MAXD];
	double		*h = gr->h;
	double		*crds, *crds1;
	int		i, dim = gr->dim;
	ORIENTATION	c_orient;
	ANGLE_DIRECTION	ang_dir;
	SIDE		r_or_l_state;
	size_t		sizest = size_of_state(intfc);

	debug_print("find_bubble_state","Entered find_bubble_state()\n");

	for (i = 0; i < dim; i++)
	    Coords(&P)[i] = coords[i];
	/* nearest_interface_point() will make sure curve is not a
	   subdomain boundary with the NO_SUBDOMAIN flag set*/
	if (nearest_interface_point(Coords(&P),comp,intfc,NO_SUBDOMAIN,
		                    NULL,coords_on,&t,&hse,&hs) != YES)
	{
	    screen("ERROR in find_bubble_state(), "
	           "Nearest interface point failed\n");
	    clean_up(ERROR);
	}
	bond = Bond_of_hse(hse);
	curve = Curve_of_hs(hs);

	(void) vector_product_on_points(Coords(bond->start),Coords(&P),
					Coords(bond->end),2,&crp);
	if (crp > 0.0 )
	    r_or_l_state = POSITIVE_SIDE;
	else if(crp < 0.0 )
	    r_or_l_state = NEGATIVE_SIDE;
	else if (comp == positive_component(curve))
	    r_or_l_state = POSITIVE_SIDE;
	else if (comp == negative_component(curve))
	    r_or_l_state = NEGATIVE_SIDE;
	else
	{
	    screen("ERROR in find_bubble_state(), "
	           "Can't identify interior side of bubble\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < dim; i++)
	    d[i] = Coords(&P)[i] - coords_on[i];
	if (scaled_hypot(d,h,dim) < MIN_SC_SEP(intfc))
	{
	    if (r_or_l_state == POSITIVE_SIDE)
	    {
	        s1 = right_state_at_point_on_curve(bond->start,bond,curve);
	        s2 = right_state_at_point_on_curve(bond->end,bond,curve);
	    }
	    else
	    {
	        s1 = left_state_at_point_on_curve(bond->start,bond,curve);
	        s2 = left_state_at_point_on_curve(bond->end,bond,curve);
	    }
	    gt_lin_comb_states(1.0-t,t,Coords(bond->start),s1,
	    	               Coords(bond->end),s2,gr,state);
	    set_state(state,stype,state);
	    return;
	}

	alloc_state(intfc,&comb_of_states,sizest);

	tot_pt_wt = 0.0;
	pt_wt = 0.0;
	ang_dir = CLOCKWISE;
	c_orient = (r_or_l_state == NEGATIVE_SIDE) ? POSITIVE_ORIENTATION :
			                             NEGATIVE_ORIENTATION;

	loop = find_loop(curve,c_orient,ang_dir);

	if (is_subdomain_boundary(Hyper_surf(loop->last->curve)))
	{
	    st = (r_or_l_state == POSITIVE_SIDE) ?
	    	Right_state_at_node_of_o_curve(loop->first) :
	    	Left_state_at_node_of_o_curve(loop->first);
	    ft_assign(comb_of_states,st,sizest);
	}
	else
	{
	    crds = Coords(Node_of_o_curve(loop->first)->posn);
	    crds1 = Coords(Node_of_o_curve(loop->last)->posn);
	    if (r_or_l_state == POSITIVE_SIDE)
	    {
	    	st = Right_state_at_node_of_o_curve(loop->first);
	    	if(loop->first->orient != loop->last->orient)
	    	    st1 = Left_state_at_opp_node_of_o_curve(loop->last);
	    	else
	    	    st1 = Right_state_at_opp_node_of_o_curve(loop->last);
	    }
	    else
	    {
	    	st = Left_state_at_node_of_o_curve(loop->first);
	    	if(loop->first->orient != loop->last->orient)
	    	    st1 = Right_state_at_opp_node_of_o_curve(loop->last);
	    	else
	    	    st1 = Left_state_at_opp_node_of_o_curve(loop->last);
	    }
	    gt_lin_comb_states(0.5,0.5,crds,st,crds1,st1,gr,comb_of_states);
	}

	tot_pt_wt = pt_frc(separation(Node_of_o_curve(loop->first)->posn,
			   &P,dim));

	for (oc = loop->first; oc; oc = oc->next)
	{
	    if (is_subdomain_boundary(Hyper_surf(oc->curve)))
		continue;

	    if (oc->prev != NULL && oc->prev->orient != oc->orient)
	        r_or_l_state = Opposite_side(r_or_l_state);

	    bs = Bond_at_node_of_o_curve(oc);


	    for (b = Following_bond(bs,oc->orient); b != NULL;
	              	b = Following_bond(b,oc->orient)) 
	    {
                if (r_or_l_state == POSITIVE_SIDE)
	            st = right_state_at_point_on_curve(
	                Point_of_bond(b,oc->orient),b,oc->curve);	
	        else
	            st = left_state_at_point_on_curve(
	                Point_of_bond(b,oc->orient),b,oc->curve);


	        /*  Updating the state at the interior point  */
	        /*  by adding the weighted value of state at */
	        /*       the   point on a curve.           */
		
                pt_wt =
	        pt_frc(separation(Point_of_bond(b,oc->orient),&P,dim));
	        d1 = tot_pt_wt + pt_wt;
	        alpha = tot_pt_wt / d1;
	        beta = pt_wt / d1;
	        gt_lin_comb_states(alpha,beta,coords,comb_of_states,
		                   Coords(Point_of_bond(b,oc->orient)),st,
				   gr,comb_of_states);
	        tot_pt_wt = d1;
	    }

	    if (oc->next == NULL) break;
    	    crds = Coords(Opp_node_of_o_curve(oc)->posn);
	    if (is_subdomain_boundary(Hyper_surf(oc->next->curve)))
	    {
	        if (r_or_l_state == POSITIVE_SIDE)
		    st = Right_state_at_opp_node_of_o_curve(oc);
	        else
		    st = Left_state_at_opp_node_of_o_curve(oc); 
	    }
	    else
	    {
	        crds1 = Coords(Node_of_o_curve(oc->next)->posn);
	        if (r_or_l_state == POSITIVE_SIDE)
	        {
		    st = Right_state_at_opp_node_of_o_curve(oc);
		    if(oc->orient != oc->next->orient)
		        st1 = Left_state_at_node_of_o_curve(oc->next);
		    else
		        st1 = Right_state_at_node_of_o_curve(oc->next);
		}
		else
		{
		    st = Left_state_at_opp_node_of_o_curve(oc); 
		    if(oc->orient != oc->next->orient)
		        st1 = Right_state_at_node_of_o_curve(oc->next);
		    else                         
		        st1 = Left_state_at_node_of_o_curve(oc->next); 
	        }
		gt_lin_comb_states(0.5,0.5,crds,st,crds1,st1,gr,st);
	    }
	    pt_wt = pt_frc(separation(
			   Opp_node_of_o_curve(oc)->posn,&P,dim));
	    d1 = tot_pt_wt + pt_wt;
	    alpha = tot_pt_wt / d1;
	    beta = pt_wt / d1;
	    gt_lin_comb_states(alpha,beta,coords,comb_of_states,
			       crds,st,gr,comb_of_states);
	    tot_pt_wt = d1;
	}

	ft_assign(state,comb_of_states,sizest);
	set_state(state,stype,state);
	free(comb_of_states);
	free_o_curve_family(loop); 
	debug_print("find_bubble_state","Left find_bubble_state()\n");
}		/*end find_bubble_state*/


/*
*			bubble_state():
*/

EXPORT void bubble_state(
	Locstate	state,
	POINT		*p,
	CURVE		*c,
	Bubble		*bubble,
	COMPONENT	comp)
{
	int		put_in_jump, found;
	double		distance, total_length, alpha;
	double		*crds1, *crds2;
	BOND		*bond;
	CURVE		*cps, *cpe;
	Locstate	s1, s2;
	RECT_GRID	*gr = computational_grid(c->interface);
	ORIENTATION	cps_or, cpe_or;

	found = NO;
	distance = 0.0;
	total_length = 0.0;
	for (bond = c->first; bond != NULL; bond = bond->next)
	{
	    if (!found)
	    {
	    	if (bond->start == p)
		    found = YES;
	    	else
	    	{
	    	    if (bond->end == p)
			found = YES;
	    	    distance += bond_length(bond);
	    	}
	    }
	    total_length += bond_length(bond);
	}
	alpha = distance/total_length;

	put_in_jump = NO;
	crds1 = Coords(c->start->posn);
	crds2 = Coords(c->end->posn);
	if (start_status(c) == REFLECTED)
	{		/* refl shock */
	    s1 = bubble->RP->state[2];
	    s2 = bubble->bow;
	}
	else if (end_status(c) == REFLECTED)
	{		/* refl shock */
	    s1 = bubble->bow;
	    s2 = bubble->RP->state[2];
	}
	else if (start_status(c) == MACH_STEM)
	{		/* mach stem */
	    s1 = bubble->RP->state[3];
	    s2 = bubble->mach;
	}
	else if (end_status(c) == MACH_STEM)
	{		/* mach stem */
	    s1 = bubble->mach;
	    s2 = bubble->RP->state[3];
	}
	else if (start_status(c) == SLIP)
	{		/* contact */
	    if (comp == negative_component(c))	/* on mach side */
	    {
	    	s1 = bubble->RP->state[3];
	    	s2 = bubble->contact_mach;
	    }
	    else	/* on bow side */
	    {
	    	s1 = bubble->RP->state[2];
	    	s2 = bubble->contact_bow;
	    }
	}
	else if (end_status(c) == SLIP)
	{		/* contact */
	    if (comp == positive_component(c))	/* on mach side */
	    {
	    	s1 = bubble->contact_mach;
	    	s2 = bubble->RP->state[3];
	    }
	    else	/* on bow side */
	    {
	    	s1 = bubble->contact_bow;
	    	s2 = bubble->RP->state[2];
	    }
	}
	else if (node_type(c->start) == ATTACHED_B_NODE)
	{		/* boundary with attached bow shock */
	    s1 = bubble->bow;
	    cpe = find_physical_curve_at_node(c->end,&cpe_or);
	    cpe_or = Opposite_orient(cpe_or);
	    if (status_at_node(cpe,cpe_or) == SLIP)
	    {		/* mach reflection */
	    	s2 = bubble->contact_bow;
	    }
	    else
	    {		/* regular reflection */
	    	s2 = bubble->RP->state[2];
	    }
	}
	else if (node_type(c->end) == ATTACHED_B_NODE)
	{
	    cps = find_physical_curve_at_node(c->start,&cps_or);
	    cps_or = Opposite_orient(cps_or);
	    if (status_at_node(cps,cps_or) == SLIP)
	    {		/* mach reflection */
	    	s1 = bubble->contact_bow;
	    }
	    else
	    {		/* regular reflection */
	    	s1 = bubble->RP->state[2];
	    }
	    s2 = bubble->bow;
	}
	else if (node_type(c->start) == B_REFLECT_NODE)
	{		/* refl wall of reg reflection */
	    s1 = bubble->RP->state[2];
	    s2 = bubble->refl_corner;
	}
	else if (node_type(c->end) == B_REFLECT_NODE)
	{		/* refl wall of reg reflection */
	    s1 = bubble->refl_corner;
	    s2 = bubble->RP->state[2];
	}
	else if ((node_type(c->start) == NEUMANN_NODE) &&
		(node_type(c->end) == NEUMANN_NODE))
	{
	    cps = find_physical_curve_at_node(c->start,&cps_or);
	    cps_or = Opposite_orient(cps_or);
	    cpe = find_physical_curve_at_node(c->end,&cpe_or);
	    cpe_or = Opposite_orient(cpe_or);
	    if (status_at_node(cps,cps_or) == MACH_STEM)
	    {	/* wall -- mach stem to contact, ahead of corner */
	    	s1 = bubble->mach;
	    	s2 = bubble->contact_mach;
	    }
	    else if (status_at_node(cpe,cpe_or) == MACH_STEM)
	    {	/* wall -- contact to mach stem, ahead of corner */
	    	s1 = bubble->contact_mach;
	    	s2 = bubble->mach;
	    }
	    else if (status_at_node(cps,cps_or) == REFLECTED)
	    {	/* wall -- bow shock to contact, behind corner */
	    	s1 = bubble->bow;
	    	s2 = bubble->contact_bow;
	    }
	    else if (status_at_node(cpe,cpe_or) == REFLECTED)
	    {	/* wall -- contact to bow shock, behind corner */
	    	s1 = bubble->contact_bow;
	    	s2 = bubble->bow;
	    }
	    else
	    {
	    	screen("ERROR in bubble_state(), "
	    	       "Inconsistent geometry on wall\n");
		clean_up(ERROR);
	    }
	}
	else if(node_type(c->start) == NEUMANN_NODE)
	{
	    cps = find_physical_curve_at_node(c->start,&cps_or);
	    cps_or = Opposite_orient(cps_or);
	    switch(status_at_node(cps,cps_or))
	    {
	    case MACH_STEM:		/* wall -- mach stem to corner */
	    	s1 = bubble->mach;
	    	s2 = bubble->refl_corner;
	    	break;
	    case REFLECTED:		/* wall -- bow shock to corner */
	    	s1 = bubble->bow;
	    	s2 = bubble->bow_corner;
	    	break;
	    case SLIP:		/* wall -- contact to corner */
	    	s1 = bubble->contact_bow;  /* ahead or behind corner */
	    	s2 = bubble->refl_corner;
	    	put_in_jump = YES;
	    	break;
	    default:
	    	break;
	    }
	}
	else if(node_type(c->end) == NEUMANN_NODE) 
	{
	    cpe = find_physical_curve_at_node(c->end,&cpe_or);
	    cpe_or = Opposite_orient(cpe_or);
	    switch(status_at_node(cpe,cpe_or))
	    {
	    case MACH_STEM:		/* wall -- corner to mach stem */
	    	s1 = bubble->refl_corner;
	    	s2 = bubble->mach;
	    	break;
	    case REFLECTED:		/* wall -- corner to bow shock */
	    	s1 = bubble->bow_corner;
	    	s2 = bubble->bow;
	    	break;
	    case SLIP:		/* wall -- corner to contact */
	    	s1 = bubble->refl_corner;  /* ahead or behind corner */
	    	s2 = bubble->contact_bow;
	    	put_in_jump = YES;
	    	break;
	    default:
	    	break;
	    }
	}
	else 
	{
	    screen("ERROR in bubble_state(), "
	           "Inconsistent geometry of Mach bubble\n");
	    clean_up(ERROR);
	}

	if (!put_in_jump) 
	{
	    gt_lin_comb_states(1.0-alpha,alpha,crds1,s1,crds2,s2,gr,state);
	}
	else 
	{
	    double m[MAXD];
	    int i, dim = c->interface->dim;
	    g_lin_comb_states(1.0-alpha,alpha,crds1,s1,crds2,s2,gr,state);
	    for (i = 0; i < dim; i++)
		m[i] = Mom(state)[i];
	    if (alpha < 0.5)
	    {
	    	state_on_adiabat_with_dens(s1,Dens(state),state,
					   state_type(state));
	    }
	    else
	    {
	    	state_on_adiabat_with_dens(s2,Dens(state),state,
					   state_type(state));
	    }
	    for (i = 0; i < dim; i++)
		Mom(state)[i] = m[i];
	    Energy(state) += kinetic_energy(state);
            reset_gamma(state);
	}
}		/*end bubble_state*/


/*
*			make_half_bow_curve():
*
*	A half bow curve consists of two parts:
*
*		1. a straight curve from node nrefl to the point pmid
*		2. an ellipse with center on the line through nwall and
*			the corner and tangent to the straight line from
*			nrefl to mid_posn.
*
*	To obtain the ellipse part, we first translate and rotate all points 
*	so that the corner is the new origin and nend is on the new positive 
*	x-axis.
*/

EXPORT CURVE *make_half_bow_curve(
	COMPONENT	l_comp,
	COMPONENT	r_comp,
	int		num_points,
	ORIENTATION	bow_or,
	NODE		*nrefl,
	NODE		*nbow,
	double		*mid_posn,
	RECT_GRID	*rect_grid,
	Bubble		*bubble)
{
	CURVE		*bow_curve;
	NODE		*nstart, *nend;
	POINT		*p;
	double		*corner_posn = bubble->corner_posn;
	double		refl[MAXD];	/*translated refl posn*/
	double		bow[MAXD];	/*translated bow posn*/
	double		mid[MAXD];	/*translated mid posn*/
	double		q[MAXD];
	double		rq[MAXD];
	double		rrefl[MAXD];	/*trans and rotated refl posn*/
	double		rbow[MAXD];	/*trans and rotated bow posn*/
	double		rmid[MAXD];	/*trans and rotated mid posn*/
	double		rcen[MAXD];
	double		minlength;
	double		rslope;
	double		coords[MAXD];
	double		cos_b, sin_b;	/*directions of behind wall*/
	double		a, b, bsqr;
	double		delta_theta;
	int		l, i, dim = rect_grid->dim;
	static int	first = YES;
	static double	**Q = NULL;

	debug_print("make_half_bow_curve","Entered make_half_bow_curve()\n");
	if (first)
	{
		first = NO;
		bi_array(&Q,MAXD,MAXD,FLOAT);
	}
	
	if (bow_or == POSITIVE_ORIENTATION)
	{
		nstart = nrefl;
		nend = nbow;
	}
	else
	{
		nstart = nbow;
		nend = nrefl;
	}
	minlength = min(rect_grid->h[0],rect_grid->h[1]);
	bow_curve = make_curve(l_comp,r_comp,nstart,nend);

	if (bubble->is_attached)
	{
	    double len;
	    /* the ellipse creation doesn't work for this case */

	    if (insert_point_in_bond(Point(mid_posn),bow_curve->first,
				     bow_curve) != FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in make_half_bow_curve(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }

	    len = separation(nrefl->posn,bow_curve->first->end,dim);
	    cos_b = cos(bubble->bow_base_ang + 0.5*PI);
	    sin_b = sin(bubble->bow_base_ang + 0.5*PI);
	    coords[0] = corner_posn[0] + len*cos(bubble->bow_base_ang);
	    coords[1] = corner_posn[1] + len*sin(bubble->bow_base_ang);

	    if (insert_point_in_bond(Point(coords),bow_curve->last,
				     bow_curve) != FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in make_half_bow_curve(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	    return bow_curve;
	}
	else
	{
		cos_b = bubble->bw_t[0];
		sin_b = bubble->bw_t[1];
	}
	for (i = 0; i < dim; i++)
	{
		refl[i] = Coords(nrefl->posn)[i] - corner_posn[i];
		bow[i] = Coords(nbow->posn)[i] - corner_posn[i];
		mid[i] = mid_posn[i] - corner_posn[i];
	}
	if (debugging("make_half_bow_curve"))
	{
		print_general_vector("\tnrefl = ",Coords(nrefl->posn),dim,"\n");
		print_general_vector("\tnbow = ",Coords(nbow->posn),dim,"\n");
		print_general_vector("\tmid_posn = ",mid_posn,dim,"\n");
		print_general_vector("\tcorner_posn = ",corner_posn,dim,"\n");
		print_general_vector("\trefl = ",refl,dim,"\n");
		print_general_vector("\tbow = ",bow,dim,"\n");
		print_general_vector("\tmid = ",mid,dim,"\n");
	}
	Q[0][0] =  cos_b;	Q[0][1] = sin_b;
	Q[1][0] = -sin_b;	Q[1][1] = cos_b;
	rotate_vector(rrefl,Q,refl,dim);
	rotate_vector(rbow,Q,bow,dim);
	rotate_vector(rmid,Q,mid,dim);
	if (debugging("make_half_bow_curve"))
	{
		print_general_vector("\trrefl = ",rrefl,dim,"\n");
		print_general_vector("\trbow = ",rbow,dim,"\n");
		print_general_vector("\trmid = ",rmid,dim,"\n");
	}
	if (separation(nrefl->posn,nbow->posn,dim) < minlength)
	{
	    if (insert_point_in_bond(Point(mid_posn),bow_curve->first,
				     bow_curve) != FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in make_half_bow_curve(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	    debug_print("make_half_bow_curve","Left make_half_bow_curve()\n");
	    return bow_curve;
	}
	if (rmid[0] > rbow[0])
	{
		screen("ERROR: mid point is over bow of bow curve\n");
		clean_up(ERROR);
	}

	/* Make an ellipse across (rbow[0],rbow[1]) and (rmid[0],rmid[1]) */

	rcen[1] = 0.0;
	if (fabs(rrefl[0] - rmid[0]) < EPS)
	{
		debug_print("make_half_bow_curve","Degenerate ellipse\n");
		debug_print("make_half_bow_curve","Left make_half_bow_curve()\n");
		return bow_curve;
	}
	else if (fabs(rrefl[1] - rmid[1]) < EPS)
	{
		rcen[0] = rmid[0];
		a = rbow[0] - rmid[0];
		b = rmid[1];
	}
	else
	{
		rslope = (rrefl[1] - rmid[1])/(rrefl[0] - rmid[0]);
		debug_print("make_half_bow_curve","rslope = %g\n",rslope);
		if (fabs(rslope) < EPS)
		{
			rcen[0] = rmid[0];
			a = rbow[0] - rmid[0];
			b = rmid[1];
		}
		else
		{
			rcen[0] = (sqr(rbow[0]) - sqr(rmid[0]) +
					rmid[0]*rmid[1]/rslope)   /
				(rmid[1]/rslope - 2.*rmid[0] + 2.*rbow[0]);
			a = rbow[0] - rcen[0];
			bsqr = (rmid[1]*rslope*a*a)/(rcen[0]-rmid[0]);
			if (bsqr < 0.0)
			{
			/*
			*	In this case, the ellipse is not tangent to the
			*	line ( nrefl, mid_posn ).
			*/
				rcen[0] = rmid[0];
				a = rbow[0] - rmid[0];
				b = fabs(rmid[1]);
				debug_print("make_half_bow_curve",
					"Cannot make a tangent bow curve\n");
			}
			else b = sqrt(bsqr);
		}
	}
	if (debugging("make_half_bow_curve"))
	{
		(void) printf("a = %g  b = %g\n",a,b);
		print_general_vector("\trcen = ",rcen,dim,"\n");
	}
	delta_theta = mag_vector(bow,dim)/(a + b);
	if (rmid[1] < 0.0) delta_theta = - delta_theta;
	Q[0][0] = cos_b;	Q[0][1] = -sin_b;
	Q[1][0] = sin_b;	Q[1][1] =  cos_b;
	for (l = 1; l < num_points; l++)
	{
		rq[0] = rcen[0] + a*cos(l*delta_theta);
		rq[1] = rcen[1] + b*sin(l*delta_theta);
		if (debugging("make_half_bow_curve"))
		    print_general_vector("\trq = ",rq,dim,"\n");
		if (rq[0] <= rmid[0]) break;
		else
		{
		    rotate_vector(q,Q,rq,dim);
		    for (i = 0; i < dim; i++)
			coords[i] = corner_posn[i] + q[i];
		    p = Point(coords);
		    insert_point_adjacent_to_node(p,bow_curve,bow_or);
		    if (debugging("make_half_bow_curve"))
		    {
			print_general_vector("\tq = ",q,dim,"\n");
			print_general_vector("\tp = ",Coords(p),dim,"\n");
		    }
		}
	}
	debug_print("make_half_bow_curve","Left make_half_bow_curve()\n");
	return bow_curve;
}		/*end make_half_bow_curve*/
#endif /* defined(TWOD) */
