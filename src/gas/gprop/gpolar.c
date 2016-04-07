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
*				gpolar.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains routines for the analysis of a single shock or
*	rarefaction wave in two dimensions by the method of shock polars.
*
*	Basic routines:
*
*		s_polar_2()
*		s_polar_3()
*		s_polar_4()
*		s_polar_5()
*		prandtl_meyer_wave()
*
*	Auxiliary frame transformation routines:
*		steady_to_unsteady()
*		unsteady_to_steady()
*
*	Shock polar routines calculate the properties of a two dimensional
*      	oblique shock.  A single oblique shock can be described by the 
*      	following parameters:
*   
*          ahead_states/back_states: density, pressure, vx, vy
*          shock_angle, 
*          shock_speed,
*          shock_strength, defined to be, for instance, ahead_pr/back_pr.
*          turn_angle,  this is the turning angle from the ahead velocity
*                       to the velocity behind the shock, in some prescribed
*			rest frame for the shock.
*         
*      	These parameters are certainly not independant.  The following 
*	routines s_polar_j() assume various sets of input parameters
*	and obtain some of the remaining parameters as output.  There
*	are also some discrete parameters, or flags, to guide the calculation:
*
*	   is_ahead_given,
*	   is_ans_weak,
*	   is_interior,
*	   is_turn_angle_pos,
*	   which_parameter.
*
*      	In these routines, the given state could be either the ahead state
*	or the behind state; which one is given is determined by whether 
*	the input argument is_given_ahead is YES or NO.
*
*	All angles are measured from the positive x axis in the
*	counterclockwise direction to some other given half line or ray.
*	When the shock angle is an input variable, it defines a ray or
*	half line for the shock. Similarly when the shock angle is an
*	output variable, some additional discrete flag is needed to select
*	the half line or ray of the shock which the angle refers to.
*	The streamline through the shock will bend at the shock.  The side
*	of the streamline which forms an angle less than PI is called the
*	interior side of the streamline; also the portion of the shock on
*	the interior side of the streamline is a ray, called the interior
*	ray of the shock.
*	If is_interior is set to YES, then the shock angle is defined to be
*	the angle from the positive x axis to the interior ray of the shock.
*	Otherwise the opposite ray is used to define the shock angle.
*
*	There are several frames used within this file.
*
*	The frame of the data is called the computational frame. All input 
*	and output data (with the exception of the turning angle) are given 
*	in this frame.
*
*	The steady frame is defined by given data, abs_vx, abs_vy and is the
*	same for all shocks and states. 
*
*	The rotated steady frame is given by a rotation of the steady frame 
*	so that the steady frame velocity ahead of the shock is directed 
*	along the positive x axis.  This frame is depends on the state ahead, 
*	and is in general distinct for each shock.
*
*	A normal frame of a shock is defined to be a frame in which the
*	flow is normal to the shock and in the direction of the positive
*	x axis.
*
*	An ambient frame is a frame in which the state ahead of a shock is
*	at rest.
*
*	All input and output is given in terms of conserved quantities 
*	(momentum, density and energy), while the computations internal to 
*	these routines generally employ thermodynamic variables (pressure,
*	velocity and density).
*
*      	Reference: Courant & Friedrichs, Supersonic Flow and Shock Waves.
*/


#include <gdecs/gdecs.h>

typedef	struct {
	Locstate	st0;
	double		p0, V0;
	double		rq2;
	double		tan_theta; /* (+-) tan(theta) */
	double		numerator; /* (+-) (rho0/p0) * q0^2 * tan(theta) */
} TA_PRMS;


	/* LOCAL Function Declarations */
LOCAL	boolean	f_turn_ang(double,double*,POINTER);
LOCAL	boolean	ta_f(double,double*,TA_PRMS*);
LOCAL	boolean	pr_given_turn_ang_Hugoniot(double*,double,double,int,int,
					   Locstate);
LOCAL	boolean	velocity_behind_oblique_shock(double,Locstate,int,int,
					      double,double,double,double**,
					      double,double*,double*,double*);
LOCAL	boolean	unsteady_to_steady(double*,double*,double*,double**,double*,int );
LOCAL	void	steady_to_unsteady(double*,double*,double*,double**,
				   double,double*,int );
LOCAL	void	zero_strength_shock(Locstate,Locstate,double,double,
				    int,int,double*,double*);

/*
*			s_polar_2():
*
*       If the ahead state, the transformation to the steady frame and 
*	the turning angle are given, there are two possible
*       states behind an oblique shock.  Similarly, the behind state and
*	turning angle determine two possible states ahead of an oblique shock.
*	One forms a weak shock, while the other is the strong shock. 
*       Assigning is_ans_weak to be YES gives the weak shock and NO gives
*       the strong one. Similarly, is_given_ahead specifies whether the ahead
*	or behind state is given.
*
*       Input:
*		given_st1
*		is_given_ahead
*		is_ans_weak
*		turn_angle	(the difference of the angle : answer velocity -
*			  given_st velocity)
*		abs_v		(defines the change to the steady frame)
*	Output
*		answer1
*		pshock_angle	(the angle from the positive x axis to the
*				 shock line)
*
*	The shock_angle is normalized so that the vector with direction 
*	shock_angle is on the clockwise side of the vector (abs_vx,abs_vy).
*/

EXPORT int s_polar_2(
	Locstate	given_st1,
	int		is_given_ahead,
	int		is_ans_weak,
	double		turn_angle,
	double		*abs_v,
	Locstate	answer1,
	double		*pshock_angle)
{
	double		vel_ang;  /* rest frame velocity angle of given state */
	double		r0,r1;    /* given_st ,answer_st density */
	double		p0,p1;    /* given_st, answer_st pressure */
	double		u1,v1;    /* rotated rest frame velocity of answer_st */
	double		q0;       /* rest frame speed of given state */
	int		dim;
	static Locstate given_st = NULL,answer = NULL;	/* TGas versions of */
	static double	**Q = NULL;			/* rotation matrix */
	static boolean	first = YES;		/* given_st1 and answer1 */

	debug_print("spolar2","Entering s_polar_2()\n");

	if ((fabs(fabs(turn_angle) - PI/2.)) < EPS)
	{
		(void) printf("WARNING in s_polar_2(), ");
		(void) printf("Turning angle is too close to 90 degrees.\n");
		return 0;
	}

		/* Allocate, translate and rotate states */

	if (first)
	{
		size_t		sizest = Params(given_st1)->sizest;

		first = NO;
		(*Params(given_st1)->_alloc_state)(&given_st,sizest);
		(*Params(given_st1)->_alloc_state)(&answer,sizest);
		bi_array(&Q,SMAXD,SMAXD,FLOAT);
	}
	dim = Params(given_st1)->dim;
	set_state(given_st,TGAS_STATE,given_st1);

	if (unsteady_to_steady(Vel(given_st),abs_v,&q0,Q,&vel_ang,dim) ==
							FUNCTION_FAILED)
		return 0;

#if defined(DEBUG_GPOLAR)
	if (debugging("spolar2"))
	{
		double c0 = sound_speed(given_st);
		
		(void) printf("State and turn angle into s_polar_2()\n");
		verbose_print_state("given_st1",given_st1);
		print_angle("Turning angle =",turn_angle,"\n");
		(void) printf("abs_v = %g, %g\n",abs_v[0],abs_v[1]);
		(void) printf("q0 = %g\n",q0);
		(void) printf("Steady frame Mach number = %g\n",q0/c0);
		(void) printf("is_given_ahead = %d,\tis_ans_weak = %d\n",
			      is_given_ahead,is_ans_weak);
	}
#endif /* defined(DEBUG_GPOLAR) */
		/* Identify opposite state */

	r0 = Dens(given_st);
	p0 = Press(given_st);
	if (pr_given_turn_ang_Hugoniot(&p1,turn_angle,q0,is_ans_weak,
				   is_given_ahead,given_st) == FUNCTION_FAILED)
		return 0;
	r1 = dens_Hugoniot(p1,given_st);
	u1 = q0 - (p1 - p0)/(q0*r0);
	v1 = u1*tan(turn_angle);
	*pshock_angle = angle(v1,(q0-u1));
#if defined(DEBUG_GPOLAR)
	if (debugging("spolar2"))
	{
		double	cta;
		double	M0sq = sqr(q0) / sound_speed_squared(given_st);

		if (steady_state_wave_curve(p1,M0sq,&cta,given_st1) ==
							FUNCTION_FAILED)
		{
			(void) printf("WARNING in s_polar_2(), ");
			(void) printf("steady_state_wave_curve() failed.\n");
			return NO;
		}
		if (turn_angle < 0.0) cta = -cta;

		(void) printf("Inside of s_polar_2, checking results of p1");
		(void) printf(" computation\n");
		(void) printf("turn_angle = %g\n",turn_angle);
		(void) printf("computed turn_angle = %g\n",cta);
		(void) printf("absolute value of difference = %g\n",
			      fabs(turn_angle - cta));
	}
	debug_print("spolar2",
	      "shock angle (relative to the given flow direction) = %g\n",
	      *pshock_angle);
#endif /* defined(DEBUG_GPOLAR) */
	Set_params(answer,given_st);
	set_type_of_state(answer,TGAS_STATE);
	Dens(answer) = r1;
	Press(answer) = p1;
	Vel(answer)[0] = u1;
	Vel(answer)[1] = v1;
	reset_gamma(answer);
#if defined(DEBUG_GPOLAR)
	debug_print("spolar2","r1 = %g, p1 = %g, u1 = %g, v1 = %g\n",r1,p1,u1,v1);
#endif /* defined(DEBUG_GPOLAR) */
	steady_to_unsteady(Vel(answer),abs_v,Vel(answer),Q,
			vel_ang,pshock_angle,dim);
	set_state(answer1,GAS_STATE,answer);
	if (cos(*pshock_angle)*abs_v[0] + sin(*pshock_angle)*abs_v[1] > 0. )
		*pshock_angle= normalized_angle( PI + *pshock_angle);
#if defined(DEBUG_GPOLAR)
	if (debugging("spolar2"))
	{
	    (void) printf("Answer state and shock angle from s_polar_2()\n");
	    verbose_print_state("answer1",answer1);
	    (void) printf("shock_angle = %g\n",*pshock_angle);
	}
#endif /* defined(DEBUG_GPOLAR) */

	debug_print("spolar2","Leaving s_polar_2()\n");
	return 1;
}		/*end s_polar_2*/


/*
*			s_polar_3():
*
*       Given a state on one side of a shock and the pressure on the 
*	other side, this routine solves for the oblique shock.  Roughly 
*	speaking, this data is the same as giving the shock strength.
*	There are two solutions; one has a positive turn_angle, and the
*	other has a negative turn_angle.
*       Input:
*		given_state
*		is_ahead_given
*		p1		(pressure of the other side)
*		is_turn_ang_positive
*		is_interior
*		abs_v		(defines the change to the steady frame)
*       Output:
*		answer
*		pshock_angle	(the angle from the positive x axis to the
*				 shock line)
*      		pturn_angle	(the difference of the angles : answer velocity
*					- given_st velocity)
*
*	The streamline through the shock will bend at the shock.  The side
*	of the streamline which forms an angle less than PI is called the
*	interior side of the streamline; also the portion of the shock on
*	the interior side of the streamline is a ray, called the interior
*	ray of the shock. In case the shock is part of a 2 dimensional
*	elementary wave configuration, then there is a natural definition
*	of a ray or half line given by the part of the shock realized in
*	physical space.  For incoming shocks, this ray is the exterior portion
*	and for outgoing shocks, the ray is the incoming portion.
*	If is_interior is set to YES, then *pshock_angle gives the angle 
*	from the positive x axis to the interior ray of the shock.  Otherwise
*	the opposite ray is used to define the shock angle. To further clarify
*	these concepts, we note that for a shock which is part of a two
*	dimensional elementary wave, there is automatically a ray or half line
*	selected. For an incoming such shock, the shock is exterior while the
*	outgoing case is interior.
*/

EXPORT int s_polar_3(
	Locstate	given_st,
	int		is_given_ahead,
	double		p1,		/* pressure on the other side  */
	int		is_turn_ang_positive,
	int		is_interior,
	double		*abs_v,		/* absolute velocity */
	Locstate	answer,
	double		*pshock_angle,
	double		*pturn_angle)
{
	int		status;
	double  		vel_ang;  /* rest frame vel angle of the given state */
	double	 	r0;	  /* given_st ,answer_st density */
	double		p0;	  /* given_st pressure */
	double		q0;	  /* rest frame speed of given state */
	double		c0;	  /* given_st sound speed */
	double		M0;	  /* M0 = q0/c0 Mach number of given state */
	int		dim;
	static Locstate vstate = NULL;
	static double	**Q = NULL;

	debug_print("spolar3","Entering s_polar_3()\n");

	if (Q == NULL)
	{
		(*Params(given_st)->_alloc_state)(&vstate,sizeof(VGas));
		bi_array(&Q,SMAXD,SMAXD,FLOAT);
	}
	dim = Params(given_st)->dim;

	set_state(vstate,VGAS_STATE,given_st);

	if (unsteady_to_steady(Vel(vstate),abs_v,&q0,Q,&vel_ang,dim) ==
							FUNCTION_FAILED)
	{
		(void) printf("WARNING - Failure in s_polar_3(), ");
		(void) printf("unsteady_to_steady() failed\n");
		debug_print("spolar3","Leaving s_polar_3()\n");
		return NO;
	}

	r0 = Dens(vstate);
	p0 = pressure(vstate);
        c0 = sound_speed(vstate);
        M0 = q0/c0;
#if defined(DEBUG_GPOLAR)
	if (debugging("spolar3"))
	{
	    (void) printf("abs_v = %g, %g\n",abs_v[0],abs_v[1]);
	    (void) printf("Behind pressure = %g\n",p1);
	    (void) printf("Pressure difference = %g\n",
			  p1 - pressure(vstate));
	    (void) printf("is_given_ahead = %s, is_turn_ang_positive = %s, ",
			  (is_given_ahead) ? "YES" : "NO",
			  (is_turn_ang_positive) ? "YES" : "NO");
	    (void) printf("is_interior = %s\n",
			  (is_interior) ? "YES" : "NO");
	    (void) printf("q0 = %g\n",q0);
	    print_state_type("Input state type = ",state_type(given_st));
	    verbose_print_state("Input state",vstate);
	    (void) printf("Given state steady frame mach number = %g\n\n",M0);
	}
#endif /* defined(DEBUG_GPOLAR) */
	if ((is_given_ahead == YES) && (M0 < SONIC_MINUS))
	{
		(void) printf("WARNING - Failure in s_polar_3(), ");
		(void) printf( "The ahead state is subsonic, M0 = %g\n",M0);
		debug_print("spolar3","Leaving s_polar_3()\n");
		return NO;
	}
	M0 = max(1.0,M0);

	if (fabs(p1 - p0) < EPS*r0*q0)
	{
		zero_strength_shock(given_st,answer,M0,vel_ang,
			is_turn_ang_positive,is_interior,pturn_angle,
			pshock_angle);
		debug_print("spolar3","Leaving s_polar_3()\n");
		return YES;
	}

	state_w_pr_on_Hugoniot(vstate,p1,answer,state_type(given_st));

#if defined(DEBUG_GPOLAR)
	if (debugging("spolar3")) 
	{
		(void) printf("Mass flux squared = %g\n",
			      mass_flux_squared(p1,given_st));
		verbose_print_state("Thermodynamic answer state",answer);
	}
#endif /* defined(DEBUG_GPOLAR) */

	status = velocity_behind_oblique_shock(p1,answer,
			is_turn_ang_positive,is_interior,q0,p0,r0,
			Q,vel_ang,abs_v,pturn_angle,pshock_angle);

	debug_print("spolar3","Leaving s_polar_3()\n");
	return status;
}		/*end s_polar_3*/



/* 
*			s_polar_4():
*
*	A one-dimensional shock (and a shock in two dimensions
*	for which the direction normal to the shock is given) is
*	determined by the ahead state and one additional parameter.
*	Input:
*		ahead state
*		shock normal (oriented from high to low pressure sides)
*		any one of the following as an additional parameter:
*
*			1. the pressure behind
*			2. the density behind
*			3. the flow velocity behind
*			4. the shock speed
*			5. the Mach number of the shock speed, relative to the
*			sound speed in the ahead state.
*		which_parameter
*	Output:
*		behind state
*		shock speed
*	We note that the shock speed and Mach number refer to the computational
*	frame, i. e. the frame of the input and output data.
*/

EXPORT int s_polar_4(
	int		which_parameter,
	double		parameter,
	double		*shock_speed,
	double		*shock_nor,
	Locstate	ahead_state1,
	Locstate	behind_state1,
	int		behind_stype)
{
	double		r0, p0, u0, vtan[SMAXD];/* ahead state */
	double		r1, p1, u1;		/* behind state */
	double		U;			/* shock speed */
	double		M0n;			/* shock mack number,
						   relative to ahead flow */
	double		m;			/* mass flux */
	double		M0nsq;			/* steady normal ahead Mach
						   number squared */
	int		i, dim;
	static Locstate ahead_state = NULL, behind_state = NULL;

	debug_print("spolar4","Entering s_polar_4() (parameter=%d):\n",
		which_parameter);

	if (behind_state == NULL)
	{
	    size_t sizest = Params(ahead_state1)->sizest;

	    (*Params(ahead_state1)->_alloc_state)(&ahead_state,sizest);
	    (*Params(ahead_state1)->_alloc_state)(&behind_state,sizest);
	}
	dim = Params(ahead_state1)->dim;
	set_state(ahead_state,TGAS_STATE,ahead_state1);

	r0 = Dens(ahead_state);
	p0 = pressure(ahead_state);
	u0 = scalar_product(Vel(ahead_state),shock_nor,dim);
	for (i = 0; i < dim; i++)
	    vtan[i] = Vel(ahead_state)[i] - u0*shock_nor[i];
	switch(which_parameter)
	{
	case BEHIND_PRESSURE:
#if defined(DEBUG_GPOLAR)
	    debug_print("spolar4","parameter = BEHIND_PRESSURE = %g\n",parameter);
#endif /* defined(DEBUG_GPOLAR) */
	    p1 = parameter;
	    r1 = dens_Hugoniot(p1,ahead_state);
	    m = mass_flux(p1,ahead_state);
	    *shock_speed = u0 + m/r0;
	    u1 = *shock_speed - m/r1;
	    break;
	case BEHIND_VELOCITY:
#if defined(DEBUG_GPOLAR)
	    debug_print("spolar4","parameter = BEHIND_VELOCITY = %g\n",parameter);
#endif /* defined(DEBUG_GPOLAR) */
	    u1 = parameter;
	    p1 = pr_normal_vel_wave_curve(fabs(u1-u0),ahead_state);
	    r1 = dens_Hugoniot(p1,ahead_state);
	    *shock_speed = u0 + r1*(u1 - u0)/(r1 - r0);
	    break;
	case SHOCK_SPEED:
#if defined(DEBUG_GPOLAR)
	    debug_print("spolar4","parameter = SHOCK_SPEED = %g\n",parameter);
#endif /* defined(DEBUG_GPOLAR) */
	    *shock_speed = U = parameter;	/* shock_velocity */
	    M0nsq = sqr(u0 - U) / sound_speed_squared(ahead_state);
	    p1 = max_behind_shock_pr(M0nsq,ahead_state);
	    u1 =  u0 + (p0 - p1) / (r0*(u0 - U)); 
	    r1 = r0*((u0 - U)/(u1 - U));
	    break;
	case SHOCK_MACH_NUMBER:
#if defined(DEBUG_GPOLAR)
	    debug_print("spolar4","parameter = SHOCK_MACH_NUMBER = %g\n",parameter);
#endif /* defined(DEBUG_GPOLAR) */
	    M0n = parameter;
	    *shock_speed = U = u0 + M0n * sound_speed(ahead_state);
	    M0nsq = sqr(M0n);
	    p1 = max_behind_shock_pr(M0nsq,ahead_state);
	    u1 =  u0 + (p0 - p1) / (r0*(u0 - U)); 
	    r1 = r0*((u0 - U)/(u1 - U));
	    break;
	default:
	    screen("ERROR in s_polar_4(), "
	           "unknown parameter %d\n",which_parameter);
	    clean_up(ERROR);
		return NO;
	}	
	Dens(behind_state) = r1;
	Press(behind_state) = p1;
	for (i = 0; i < dim; i++)
	    Vel(behind_state)[i] = u1 * shock_nor[i] + vtan[i];
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* preserve ahead state mass fraction */
            Gas_param *params = Params(ahead_state);
            if(params->n_comps != 1)
            {
                for(i = 0; i < params->n_comps; i++)
                    pdens(behind_state)[i] = pdens(ahead_state)[i]/Dens(ahead_state)
                                             *Dens(behind_state);
            }
        }
	reset_gamma(behind_state);
	Set_params(behind_state,ahead_state);
	set_type_of_state(behind_state,TGAS_STATE);
	set_state(behind_state1,behind_stype,behind_state);
#if defined(DEBUG_GPOLAR)
	if (debugging("spolar4"))
	{
	    print_general_vector("shock normal = ",shock_nor,dim,"\n");
	    verbose_print_state("Ahead state",ahead_state1);
	    verbose_print_state("Behind state",behind_state1);
	    (void) printf("shock_speed = %g\n",*shock_speed);
	}
#endif /* defined(DEBUG_GPOLAR) */

	debug_print("spolar4","Leaving s_polar_4()\n");
	return YES;
}		/*end s_polar_4*/

/*
*			s_polar_5():
*
*       Given a state on one side of a shock and the mass flux across 
*	the shock, this routine solves for the oblique shock.  Roughly 
*	speaking, this data is the same as giving the shock strength.
*	There are two solutions; one has a positive turn_angle, and the
*	other has a negative turn_angle.
*       Input:
*		given_st
*		is_ahead_given
*		mfsqr		(square of the mass flux across the shock)
*		is_turn_ang_positive
*		is_interior
*		abs_v		(defines the change to the steady frame)
*       Output:
*		answer state
*		pshock_angle	(the angle from the positive x axis to the
*				 shock line)
*      		pturn_angle	(the difference of the angles : answer velocity
*					- given_st velocity)
*
*	The streamline through the shock will bend at the shock.  The side
*	of the streamline which forms an angle less than PI is called the
*	interior side of the streamline; also the portion of the shock on
*	the interior side of the streamline is a ray, called the interior
*	ray of the shock. In case the shock is part of a 2 dimensional
*	elementary wave configuration, then there is a natural definition
*	of a ray or half line given by the part of the shock realized in
*	physical space.  For incoming shocks, this ray is the exterior portion
*	and for outgoing shocks, the ray is the incoming portion.
*	If is_interior is set to YES, then *pshock_angle gives the angle 
*	from the positive x axis to the interior ray of the shock.  Otherwise
*	the opposite ray is used to define the shock angle. To further clarify
*	these concepts, we note that for a shock which is part of a two
*	dimensional elementary wave, there is automatically a ray or half line
*	selected. For an incoming such shock, the shock is exterior while the
*	outgoing case is interior.
*/

EXPORT int s_polar_5(
	Locstate	given_st,
	int		is_given_ahead,
	double		mfsqr,	 /* square of the mass flux across the shock */
	int		is_turn_ang_positive,
	int		is_interior,
	double		*abs_v,		/* absolute velocity */
	Locstate	answer,
	double		*pshock_angle,
	double		*pturn_angle)
{
	int		status;
	double   	vel_ang;  /* rest frame vel angle of the given state */
	double   	r0;	  /* given_st ,answer_st density */
	double   	p0;	  /* given_st pressure */
	double   	q0;	  /* rest frame speed of given state */
	double   	c0;	  /* given_st sound speed */
	double		aisqr;	  /* given_st acoustic impedance squared */
	double		M0;	  /* M0 = q0/c0 Mach number of given state */
	int		dim;
	static	Locstate vstate = NULL;
	static	double	**Q = NULL;

	debug_print("spolar5","Entering s_polar_5()\n");

	if (vstate == NULL)
	{
		(*Params(given_st)->_alloc_state)(&vstate,sizeof(VGas));
		bi_array(&Q,SMAXD,SMAXD,FLOAT);
	}
	dim = Params(given_st)->dim;

	set_state(vstate,VGAS_STATE,given_st);

	if (unsteady_to_steady(Vel(vstate),abs_v,&q0,Q,&vel_ang,dim) ==
							FUNCTION_FAILED)
	{
		(void) printf("WARNING - Failure in s_polar_5(), ");
		(void) printf("unsteady_to_steady() failed\n");
		debug_print("spolar5","Leaving s_polar_5()\n");
		return NO;
	}

	r0 = Dens(vstate);
	p0 = pressure(vstate);
        c0 = sound_speed(vstate);
        M0 = q0/c0;

#if defined(DEBUG_GPOLAR)
	if (debugging("spolar5"))
	{
	    (void) printf("abs_v = %g, %g\n",abs_v[0],abs_v[1]);
	    (void) printf("Mass flux squared = %g\n",mfsqr);
	    (void) printf("is_given_ahead = %s, is_turn_ang_positive = %s, ",
			  (is_given_ahead) ? "YES" : "NO",
			  (is_turn_ang_positive) ? "YES" : "NO");
	    (void) printf("is_interior = %s\n",
			  (is_interior) ? "YES" : "NO");
	    (void) printf("q0 = %g\n",q0);
	    print_state_type("Input state type = ",state_type(given_st));
	    verbose_print_state("Input state",vstate);
	    (void) printf("Given state steady frame mach number = %g\n\n",M0);
	}
#endif /* defined(DEBUG_GPOLAR) */

	if ((is_given_ahead == YES) && (M0 < SONIC_MINUS))
	{
		(void) printf("WARNING - Failure in s_polar_5(), ");
		(void) printf( "The ahead state is subsonic\n");
		debug_print("spolar5","Leaving s_polar_5()\n");
		return NO;
	}
	M0 = max(1.0,M0);

	aisqr = sqr(r0*c0);
	if (fabs(mfsqr - aisqr) < EPS*aisqr)
	{
		zero_strength_shock(given_st,answer,M0,vel_ang,
			is_turn_ang_positive,is_interior,pturn_angle,
			pshock_angle);

		debug_print("spolar5","Leaving s_polar_5()\n");
		return YES;
	}

	if (state_w_mf_sqr_on_Hugoniot(vstate,mfsqr,answer,
				       state_type(given_st)) == FUNCTION_FAILED)
	{
		(void) printf("WARNING - Failure in s_polar_5(), ");
		(void) printf("state_w_mf_sqr_on_Hugoniot() failed\n");
		debug_print("spolar5","Leaving s_polar_5()\n");
		return NO;
	}

#if defined(DEBUG_GPOLAR)
	if (debugging("spolar5")) 
	{
		verbose_print_state("Thermodynamic answer state",answer);
	}
#endif /* defined(DEBUG_GPOLAR) */

	status = velocity_behind_oblique_shock(pressure(answer),
			answer,is_turn_ang_positive,
			is_interior,q0,p0,r0,Q,vel_ang,abs_v,
			pturn_angle,pshock_angle);

	debug_print("spolar5","Leaving s_polar_5()\n");
	return status;
}		/*end s_polar_5*/


/*
*	              	prandtl_meyer_wave():
*
*       Given a state on one side of a Prandtl-Meyer rarefaction wave and the 
*	pressure on the other side, this routine finds the state on the other
*	side and the angles of the characteristic bounding the simple wave.
*       There are two solutions , one corresponding to a
*	C+ simple wave and one corresponding to a C- simple wave.  The simple
*	wave extends over the sector bounded by the rays with angles 
*	angle0 and angle1.
*	C+ rarefactions turn the flow in the clockwise direction
*	and C- rarefactions turn the flow in the counter clockwise direction.
*       Input:
*		given_state
*		pressure of the other side
*		is_C_minus_wave
*		abs_vx,abs_vy (defines the change to the steady frame)
*       Output:
*		answer state
*		angle0
*		angle1
*/

EXPORT int prandtl_meyer_wave(
	Locstate	given_st,
	double		p1,
	boolean		is_C_plus_wave,
	double		*abs_v,
	Locstate	answer,
	double		*angle0,
	double		*angle1,
	double		*pturn_angle)
{
	double   	vel_ang;	/* velocity angle of the given state */
	double        	q0, q1;		/* speed of given/answer state*/
	double       	q0sqr, q1sqr;	/* speed squared of given/answer state*/
	double		A0, A1;		/* Mach angles */
	double 		H1, H0;
	double		M0sqr, M1sqr;
	int 		dim;
	static Locstate vstate = NULL;		/* VGas versions of given_st */
	static double	**Q = NULL;

	debug_print("prandtl","Entering prandtl_meyer_wave()\n");

	if (Q == NULL) 
	{
		(*Params(given_st)->_alloc_state)(&vstate,sizeof(VGas));
		bi_array(&Q,SMAXD,SMAXD,FLOAT);
	}
	dim = Params(given_st)->dim;
	set_state(vstate,VGAS_STATE,given_st);

#if defined(DEBUG_GPOLAR)
	if (debugging("prandtl"))
	{
		(void) printf("abs_v = (%g, %g)\n",abs_v[0],abs_v[1]);
		print_state_type("Input state type = ",state_type(given_st));
		verbose_print_state("Input state",vstate);
		(void) printf("Mach number = %g\n",mach_number(vstate,abs_v));
		(void) printf("is_C_plus_wave = %s, ",
			      (is_C_plus_wave) ? "YES" : "NO");
		(void) printf("Behind pressure = %g\n",p1);
	}
#endif /* defined(DEBUG_GPOLAR) */

	if (unsteady_to_steady(Vel(vstate),abs_v,&q0,Q,&vel_ang,dim) ==
							FUNCTION_FAILED)
		return NO;
	q0sqr = sqr(q0);
	M0sqr = q0sqr/sound_speed_squared(vstate);
	if (M0sqr < SONIC_MINUS_SQR)
		return NO;
	else if (M0sqr < 1.0)
	{
		M0sqr = 1.0;
		A0 = 0.5*PI;
	}
	else
		A0 = asin(1.0/sqrt(M0sqr));
	H0 = specific_enthalpy(vstate);

	state_on_adiabat_with_pr(vstate,p1,answer,EGAS_STATE);
	H1 = specific_enthalpy(answer);


#if defined(DEBUG_GPOLAR)
	if (debugging("prandtl"))
	{
		(void) printf("Given state steady frame mach number = %g\n\n",
			      sqrt(M0sqr));
		(void) printf("pressure answer = %g\n",pressure(answer));
	}
#endif /* defined(DEBUG_GPOLAR) */

	q1sqr = q0sqr + 2.0*(H0 - H1);
	M1sqr = q1sqr/sound_speed_squared(answer);
	if (M1sqr < SONIC_MINUS_SQR)
		return NO;
	else if (M1sqr < 1.0)
	{
		M1sqr = 1.0;
		A1 = 0.5*PI;
	}
	else
		A1 = asin(1.0/sqrt(M1sqr));
	q1 = sqrt(q1sqr);

	if (steady_state_wave_curve(p1,M0sqr,pturn_angle,vstate) ==
							FUNCTION_FAILED)
	{
		return NO;
	}
	if (is_C_plus_wave)
	{
		*pturn_angle = -fabs(*pturn_angle);
		*angle0 = A0;
		*angle1 = A1 + *pturn_angle;
	}
	else
	{
		*pturn_angle = fabs(*pturn_angle);
		*angle0 = -A0;
		*angle1 = *pturn_angle - A1;
	}

	Vel(answer)[0] = q1*cos(*pturn_angle);
	Vel(answer)[1] = q1*sin(*pturn_angle);

#if defined(DEBUG_GPOLAR)
	if (debugging("prandtl"))
	{
		(void) printf("Before rotation\n");
		print_angle("angle0 =",*angle0,"\n");
		print_angle("angle1 =",*angle1,"\n");
		print_angle("*pturn_angle =",*pturn_angle,"\n");
	}
#endif /* defined(DEBUG_GPOLAR) */

		/* Rotate and translate back */

	steady_to_unsteady(Vel(answer),abs_v,Vel(answer),
		Q,vel_ang,angle0,dim);
	*angle1 = normalized_angle(*angle1 + vel_ang);

	set_state(answer,state_type(given_st),answer);
#if defined(DEBUG_GPOLAR)
	if (debugging("prandtl")) 
	{
		(void) printf("After rotations\n");
		(void) printf("angle0 = %g, angle1 = %g\n",*angle0,*angle1);
		verbose_print_state("answer state",answer);
		(void) printf("answer state steady frame mach number = %g\n",
			      sqrt(M1sqr));
	}
#endif /* defined(DEBUG_GPOLAR) */

	debug_print("prandtl","Leaving prandtl_meyer_wave()\n");
	return YES;
}		/*end prandtl_meyer_wave*/


/*
*			unsteady_to_steady():
*
*	This routine translates to the frame moving with the velocity
*	abs_vx,abs_vy.  In this new frame the velocity cosines, speed and
*	the velocity angle (angle from the x axis to the velocity direction)
*	are returned, but the state itself is not modified.
*/

LOCAL boolean unsteady_to_steady(
	double		*lab_v,		/* Lab frame velocities	*/
	double		*abs_v,		/* absolute velocity */
	double		*pspeed,	/* speed of the relative state */
	double		**Q,		/* rotation matrix */
	double		*prel_ang,	/* relative velocity angle */
	int		dim)
{
	double		rel_v[SMAXD];	    	/* relative velocity */
	int		i;

	debug_print("unsteady_to_steady","Entering unsteady_to_steady()\n");
	for (i = 0; i < dim; i++)
		rel_v[i] = lab_v[i] - abs_v[i];
	*pspeed = mag_vector(rel_v,dim);
	if (*pspeed < EPS) 
	{
		(void) printf("WARNING in unsteady_to_steady(), ");
		(void) printf("relative speed is too small\n");
		debug_print("unsteady_to_steady","Leaving unsteady_to_steady()\n\n");
		return FUNCTION_FAILED;
	}
	for (i = 0; i < dim; i++) Q[i][0] = rel_v[i] / *pspeed;
	switch (dim)
	{
	case 1:
		break;
	case 2:
		Q[0][1] = -Q[1][0];
		Q[1][1] =  Q[0][0];
		break;
	case 3:
		if (Q[0][0] == 0.0 && Q[1][0] == 0.0)
		{
			Q[0][1] = 0.0;	Q[1][1] = 1.0;	Q[2][1] = 0.0;
			Q[0][1] = 1.0;	Q[1][1] = 0.0;	Q[2][1] = 0.0;
		}
		else
		{
			Q[0][1] = -Q[1][0]; Q[1][1] =  Q[0][0]; Q[2][1] = 0.0;
			(void) vector_product(Q[0],Q[1],Q[2],dim);
		}
		break;
	}
	*prel_ang = angle(rel_v[0],rel_v[1]);
#if defined(DEBUG_GPOLAR)
	if (debugging("unsteady_to_steady"))
	{
		print_general_vector("abs_v = ",abs_v,dim,"\n");
		print_general_vector("lab_v = ",lab_v,dim,"");
		print_general_vector(", rel_v = ",rel_v,dim,"\n");
		(void) printf("speed = %g, ",*pspeed);
		print_angle("rel_ang =",*prel_ang,"\n");
	}
#endif /* defined(DEBUG_GPOLAR) */
	debug_print("unsteady_to_steady","Leaving unsteady_to_steady()\n\n");
        return FUNCTION_SUCCEEDED;
}		/*end unsteady_to_steady*/


/*
*               	steady_to_unsteady():
*
*	The state and shock angle are rotated and translated from the
*	rotated steady frame back to the computational frame. The
*	transformation is defined by the other arguments, abs_vx, etc.
*/

LOCAL   void steady_to_unsteady(
	double		*rel_v,		/* Steady frame velocities */
	double		*abs_v,		/* absolute velocity */
	double		*lab_v,		/* Lab frame velocity */
	double		**Q,		/* rotation matrix */
	double		rel_ang,	/* relative velocity angle */
	double		*pshock_ang,	/* absolute shock angle */
	int		dim)
{
	double		v[SMAXD];
	int		i, j;

	debug_print("steady_to_unsteady","Entering steady_to_unsteady()\n");
	for (i = 0; i < dim; i++)
	{
		v[i] = 0.0;
		for (j = 0; j < dim; j++)
			v[i] += Q[i][j]*rel_v[j];
	}
	for (i = 0; i < dim; i++)
		lab_v[i] = v[i] + abs_v[i];
	*pshock_ang = normalized_angle(*pshock_ang + rel_ang);
#if defined(DEBUG_GPOLAR)
	if (debugging("steady_to_unsteady"))
	{
		print_general_vector("rel_v = ",rel_v,dim,"");
		print_general_vector("abs_v = ",abs_v,dim,"\n");
		print_angle("rel_ang =",rel_ang,"\n");
		print_general_vector("lab_v = ",lab_v,dim,"\n");
		print_angle("shock_ang =",*pshock_ang,"\n");
	}
#endif /* defined(DEBUG_GPOLAR) */
	debug_print("steady_to_unsteady","Leaving steady_to_unsteady()\n\n");
}		/*end steady_to_unsteady*/


/*
*			max_behind_shock_pr():
*
*	This function computes the maximum allowable pressure
*	behind an oblique shock with given ahead state state0
*	and given the Mach number of the ahead state in the frame
*	of the shock.  This is just the pressure behind a
*	normal shock with ahead state0.
*/

EXPORT double max_behind_shock_pr(
	double		M0sq,
	Locstate	state0)
{
	double		prmx, m2;
	static Locstate state1 = NULL;

	if (state1 == NULL)
		(*Params(state0)->_alloc_state)(&state1,Params(state0)->sizest);

	m2 = M0sq*acoustic_impedance_squared(state0);
	if (state_w_mf_sqr_on_Hugoniot(state0,m2,state1,TGAS_STATE) ==
							FUNCTION_FAILED)
	{
		screen("ERROR in max_behind_shock_pr(), ");
		(void) printf("state_w_mf_sqr_on_Hugoniot() failed\n");
		clean_up(ERROR);
		return ERROR_FLOAT;
	}
	prmx = Press(state1);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	prmx = max(prmx,Min_pressure(state0));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	return prmx;
}		/*end max_behind_shock_pr*/

/*
*			pr_given_turn_ang_Hugoniot():
*
*	This functions returns the pressure on the other side of an
*	oblique shock given the state information on one side
*	and the angle through which the flow is turned in the steady
*	frame of the shock.
*	There are two possiblities for the given state.  If is_given_ahead
*	is YES it is assumed that the given state is the ahead state of 
*	the shock, that is the state on the low pressure side of the
*	shock.  If is_given_ahead is NO then the given state is assumed
*	to be the state behind the shock.
*	If the ahead state is given there are in general two possible answers
*	one corresponding to a weak shock and the other corresponding to
*	a strong shock.  If is_ans_weak is YES the weak shock is
*	calculated otherwise the strong shock is calculated.
*
*	The solution is found by solving the equation
*
*	sqr(tan(theta)) = sqr(p1 - p0)/(sqr(rho0*sqr(q0) - p1 + p0)) *
*				(sqr(rho0*q0)*(tau0 - tau1)/(p1 - p0) - 1)
*
*	for the pressure p1 on the other side of the shock.
*	Here p0,rho0,tau0 are the pressure, density and specific volume
*	of the given state, q0 is the velocity of the given state in
*	the frame of the shock and theta is the turning angle.
*
*	Note: this function assumes that the state is given in terms
*	of the conserved quantities.
*	Reference: Courant and Friedrichs page 306 ff. and page 347 ff.
*/


/*
*				f_turn_ang():
*
*	Given the shock strength x of an oblique shock defined by the
*	relative change in pressure x = (p1 - p0)/p0, this fuction 
*	returns the square of the tangent of the turning angle.
*
*	Reference Courant and Friedrichs page 347.
*
*/

LOCAL boolean  f_turn_ang(
	double		x,
	double		*f,
	POINTER		parameters)
{

	if (ta_f(x,f,(TA_PRMS *)parameters) == FUNCTION_FAILED)
	{
		if (debugging("f_turn_ang"))
		{
			(void) printf("WARNING in f_turn_ang(), ");
			(void) printf("ta_f() failed\n");
		}
		return FUNCTION_FAILED;
	}
	*f = x - *f;
	return FUNCTION_SUCCEEDED;
}		/*end f_turn_ang*/

LOCAL	boolean	ta_f(
	double		x,
	double		*f,
	TA_PRMS		*prms)
{
	Locstate	st0 = prms->st0;
	double		cot_beta;
	double		csc_beta2;
	double		p0 = prms->p0, p1 = p0*(1.0 + x);
	double		V0 = prms->V0, V1;
	double		msq;
	double		i0;


	V1 = 1.0/dens_Hugoniot(p1,st0);
	i0 = acoustic_impedance_squared(st0);
	msq = (fabs(1.0 - V1/V0) < EPS || fabs(x) < EPS) ?
			i0 : (p1 - p0) / (V0 - V1);

	csc_beta2 = prms->rq2/msq;
	if (csc_beta2 < 1.0)
	{
		if (debugging("ta_f"))
		{
			(void) printf("WARNING in ta_f(), x too large\n");
			(void) printf("x = %g\n",x);
		}
		return FUNCTION_FAILED;
	}
	cot_beta = sqrt(csc_beta2 - 1.0);
	*f = prms->numerator/(cot_beta + prms->tan_theta);
	return FUNCTION_SUCCEEDED;
}		/*end ta_f*/

LOCAL	boolean pr_given_turn_ang_Hugoniot(
	double		*p1,
	double		turn_ang,
	double		q0,
	int		is_ans_weak,
	int		is_given_ahead,
	Locstate	state0)
{
	const double     meps = 10.0*MACH_EPS;/*TOLERANCE*/
	double		epsilon, delta; /* Convergence criteria */
	double		p0, M0, M0sq, c0, rho0;
	double		tmp;
	double		x0;
	double		xlow, xhigh;
	TA_PRMS		Prms;
#if !defined(COMBUSTION_CODE)
	VGas		VST;
#endif /* !defined(COMBUSTION_CODE) */
	Locstate	st0;


#if !defined(COMBUSTION_CODE)
	set_state((Locstate) &VST,VGAS_STATE,state0);
	st0 = (Locstate) &VST;
#else /* !defined(COMBUSTION_CODE) */
	st0 = state0;
#endif /* !defined(COMBUSTION_CODE) */
	Prms.p0 = p0 = pressure(st0);
	rho0 = Dens(st0);
	Prms.V0 = 1.0/rho0;
	c0 = sound_speed(st0);
	tmp = (is_given_ahead) ? fabs(tan(turn_ang)) : -fabs(tan(turn_ang));
	M0 = q0/c0;	M0sq = sqr(M0);
	Prms.st0 = st0;
	Prms.rq2 = sqr(rho0*q0);
	Prms.tan_theta = tmp;
	Prms.numerator = (q0*q0*rho0/p0)*Prms.tan_theta;

	if (!is_given_ahead) 
	{
		xlow = -1.0;
		xhigh = 0.0;
	}
	else 
	{
		if (M0 < 1.0) 
		{
			return FUNCTION_FAILED;
		}
		if (is_ans_weak) 
		{
			xlow = 0.0;
			if (pr_at_max_turn_angle(&xhigh,M0sq,st0) ==
							FUNCTION_FAILED)
			{
			    return FUNCTION_FAILED;
			}
			xhigh = xhigh / p0 - 1.0;
		}
		else 
		{
			if (pr_at_max_turn_angle(&xlow,M0sq,st0) ==
							FUNCTION_FAILED)
			{
			    return FUNCTION_FAILED;
			}
			xhigh = max_behind_shock_pr(M0sq,st0);
			xlow = xlow / p0 - 1.0;
			xhigh = xhigh / p0 - 1.0;
		}
	}
	if (is_ans_weak)
	{
	    double f0, x1, f1;

	    if (ta_f(0.0,&x1,&Prms) == FUNCTION_FAILED)
	    {
	        return FUNCTION_FAILED;
	    }
	    f0 = -x1;
	    x1 = min(x1,xhigh);
	    if (f_turn_ang(x1,&f1,(POINTER) &Prms) == FUNCTION_FAILED)
	    {
	        return FUNCTION_FAILED;
	    }
	    while (f0*f1 > 0.0)
	    {
	        x1 *= 2.0;
	        if (x1 > xhigh) break;
	        if (f_turn_ang(x1,&f1,(POINTER) &Prms) == FUNCTION_FAILED)
	        {
	            return FUNCTION_FAILED;
	        }
	    }
	    xhigh = min(x1,xhigh);
	}

	delta = 0.5*(xlow + xhigh)*EPS;  /* Scale the convergence criteria */
	delta = max(delta, meps);

	(void) f_turn_ang(xhigh, &epsilon, (POINTER) &Prms);

	epsilon *= EPS;
	epsilon = max(epsilon, meps);

	if (find_root(f_turn_ang,(POINTER) &Prms,0.0,&x0,xlow,xhigh,
		      epsilon, delta) == FUNCTION_FAILED)
	{
		return FUNCTION_FAILED;
	}
	if (ta_f(x0,&x0,&Prms) == FUNCTION_FAILED)
	{
		return FUNCTION_FAILED;
	}
	*p1 =  p0*(x0 + 1.0);
	return FUNCTION_SUCCEEDED;
}		/*end pr_given_turn_ang_Hugoniot*/

LOCAL	void zero_strength_shock(
	Locstate	given_st,
	Locstate	answer,
	double		M0,
	double		vel_ang,
	int		is_turn_ang_positive,
	int		is_interior,
	double		*pturn_angle,
	double		*pshock_angle)
{
	double		sinA, cosA, mach_ang;
	double		sign = (is_turn_ang_positive) ? 1.0 : -1.0;

#if defined(DEBUG_GPOLAR)
	if (debugging("spolar"))
	{
		(void) printf("Zero strength shock, M0 = %g\n",M0);
	}
#endif /* defined(DEBUG_GPOLAR) */

	set_state(answer,state_type(given_st),given_st);
	*pturn_angle = 0.0;
	sinA = 1.0/M0;		cosA = sqrt(1.0 - sqr(sinA));
	mach_ang = (is_interior) ? atan2(sign*sinA,cosA) :
					atan2(-sign*sinA,-cosA);
	*pshock_angle = normalized_angle(vel_ang + mach_ang);

#if defined(DEBUG_GPOLAR)
	if (debugging("spolar"))
	{
		print_angle("*pshock_angle(before rotation) = ",mach_ang,"\n");
		print_angle("*pshock_angle(after rotation) = ",
			    *pshock_angle,"\n");
	}
#endif /* defined(DEBUG_GPOLAR) */
}		/*end zero_strength_shock*/

LOCAL	boolean velocity_behind_oblique_shock(
	double		p1,
	Locstate	answer,
	int		is_turn_ang_positive,
	int		is_interior,
	double		q0,
	double		p0,
	double		r0,
	double		**Q,
	double		vel_ang,
	double		*abs_v,
	double		*pturn_angle,
	double		*pshock_angle)
{
	double		r1;		/* given_st ,answer_st density */
	double		u[SMAXD];	/* velocity of answer_st */
	double		v1sqr;		/* v1sqr = sqr(u[1]) */
	double		sign = (is_turn_ang_positive) ? 1.0 : -1.0;
	int		dim;

	debug_print("sp_vel","Entering velocity_behind_oblique_shock()\n");

	dim = Params(answer)->dim;
	r1 = Dens(answer);

	if (fabs(r0*q0) <= EPS*EPS) 
	{
		(void) printf("WARNING in velocity_behind_oblique_shock(), "
		              "rho = 0 or the given speed is too small.\n");
		debug_print("sp_vel","Leaving velocity_behind_oblique_shock()\n");
		return FUNCTION_FAILED;
	}

	u[0] = q0 - (p1-p0)/(r0*q0);
	v1sqr = (q0 - u[0])*(u[0] - r0*q0/r1);

#if defined(DEBUG_GPOLAR)
	if (debugging("sp_vel"))
		(void) printf("u[0] = %g, v1sqr = %g\n",u[0],v1sqr);
#endif /* defined(DEBUG_GPOLAR) */

	if (v1sqr  < -EPS*sqr(u[0])) 
	{
		(void) printf("WARNING in velocity_behind_oblique_shock(), ");
		(void) printf("The given pressure is not compatible with the ");
		(void) printf("given state across an oblique shock\n");
		debug_print("sp_vel","Leaving velocity_behind_oblique_shock()\n");
		return FUNCTION_FAILED;
	}

	u[1] = sign * sqrt(max(0.0,v1sqr));

	*pturn_angle = atan2(u[1],u[0]);
	*pshock_angle = (is_interior) ? atan2(sign*(q0 - u[0]),fabs(u[1])) :
					atan2(-sign*(q0 - u[0]),-fabs(u[1]));

#if defined(DEBUG_GPOLAR)
	if (debugging("sp_vel"))
	{
	    (void) printf("Steady frame u = %g, %g, q1 = %g\n",
			  u[0],u[1],mag_vector(u,dim));
	    print_angle("*pshock_angle(before rotation) =",
			*pshock_angle,"\n");
	    print_angle("*pturn_angle =",*pturn_angle,"\n");
	    (void) printf("Mass flux squared computed by upstream data = %g\n",
			  r0*r0*q0*q0*sqr(sin(*pshock_angle)));
	}
#endif /* defined(DEBUG_GPOLAR) */

		/* Rotate and translate back */

	steady_to_unsteady(u,abs_v,u,Q,vel_ang,pshock_angle,dim);

	switch(state_type(answer))
	{
	case GAS_STATE:
		Mom(answer)[0] = r1 * u[0];
		Mom(answer)[1] = r1 * u[1];
		Energy(answer) += kinetic_energy(answer);
	        reset_gamma(answer);
			
		break;
	default:
		Vel(answer)[0] = u[0];
		Vel(answer)[1] = u[1];
	}
	
#if defined(DEBUG_GPOLAR)
	if (debugging("sp_vel")) 
	{
	    double M1 = hypot(u[0]-abs_v[0],u[1]-abs_v[1])/
			sound_speed(answer);

	    verbose_print_state("answer state",answer);
	    (void) printf("answer state steady frame mach number = %g\n",M1);
	    print_angle("Turn angle =",*pturn_angle,"\n");
	    print_angle("shock angle =",*pshock_angle,"\n");
	}
#endif /* defined(DEBUG_GPOLAR) */


	debug_print("sp_vel","Leaving velocity_behind_oblique_shock()\n");
	return FUNCTION_SUCCEEDED;
}		/*end velocity_behind_oblique_shock*/
