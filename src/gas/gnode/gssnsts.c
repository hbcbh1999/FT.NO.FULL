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
*				gssnsts.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains high level routines for two dimensional shock polar analysis
*	of shock wave - shock wave interactions.
*
*
*		find_Mach_node_states()
*		find_cross_node_states()
*		find_overtake_node_states()
*/


#if defined(TWOD)
#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	double	turn_angle(double,double,double,double);
LOCAL	double	compute_behind_state_pressure(O_CURVE*,O_CURVE*,
					      RP_DATA*,Front*,Wave*);
LOCAL	int	degenerate_mach_node(O_CURVE*,O_CURVE*,O_CURVE*,POINT*,POINT*,
				     ANGLE_DIRECTION,BOND**,BOND**,POINT*,
				     Wave*,Front*,RPROBLEM**,double,double*,
				     RP_DATA*,NODE_FLAG);
LOCAL	void	temp_mnode_normal(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				  double*,Front*);
LOCAL	boolean	turn_angle_difference(double,double*,POINTER);

/*
*
*      			find_Mach_node_states():
*
*		This routine calculates the states around a Mach node. Notice
*	that we do not take into consideration any possible wall.
*		
*
*         \         state 1         /                                 
*           \ refl                /  inc                                
*             \                 /                                        
*               \             /                                           
*     state 2     \         /                                              
*                   \     /
*                     \ /
*            ---------- \            state 0
*           slip          \
*                           \
*                             \ Mach
*               state 3
*
*
*
*	We solve the states and angles of curves around the Mach node by
*	using state0 and the pressures of state1 (p1) and state3 (p3).
*	If state0 is in rest, going to a steady frame which sits at the
*	node, state0's velocity is -node velocity. Using (118.12) in
*	Courant-Friedrichs which give the relation between q0_sq and 
*	q1_sq (the speed of states 0 and 1 in the steady frame), we can
*	write the turn angles (theta01, theta12, and theta03) 
*	as functions of M0_sq (Mach speed of state0 in the steady frame)
*	using the shock polars in the theta-p plane (see 136.01).
*	The turn angles are relative to the angle of state0 which is
*	unknown. M0_sq is now the solution of the equation:
*		theta03 = theta01 - theta12.
*	After obtaining M0_sq we find the angle of the Mach node velocity
*	by intersecting the normally propagated incident shock with the
*	circle of radius q0*dt for which center is located at the old node 
*	position.
*	This defines finally the velocity of the node and the velocity of 
*	state0 in the steady frame. The states around the node and the angles
*	of the reflected shock, the contact and the Mach stem are then 
*	calculated.
*	Note that the function turn_angle() returns the absolute value of
*	the turn angle.  The proper sign must then be added according to
*	the i_to_a_dir.
*/


EXPORT int find_Mach_node_states(
	NODE		*oldn,
	O_CURVE		*oldcinc,
	O_CURVE		*newcinc,
	O_CURVE		*oldcmach,
	O_CURVE		*newcmach,
	O_CURVE		*oldcslip,
	BOND		**crossbinc,
	BOND		**crossbmach,
	POINT		*pc,
	Wave		*wave,
	Front		*fr,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	RP_DATA		*RP,
	NODE_FLAG	flag)
{
	NODE		*newn = Node_of_o_curve(newcinc);
	POINT		*oldp;
	POINT		Pcenter;
	double		*nod_v = Node_vel(newn);/*vel transf to rest frame*/
	double		r;		/* circle radius */
	double		tcr; 		/* fract dist on bond to cross */
	double		q0;		/* rest frame node speed */
	double		inc_ang;
	double		Mach_ang;
	double		refl_ang;
	double		cont_ang;
	double		theta0;		/* rest frame angle of vel0 */
	double		theta01;	/* rest frame flow turn angles */
	double		theta12;
	double		theta03;
	double		p1, p3;		/* pressures */
	double		V[MAXD];
	int		status;
	ANGLE_DIRECTION	i_to_a_dir = Opposite_ang_dir(RP->ang_dir);
	int		i, dim = fr->rect_grid->dim;
	int		iang_pos, rang_pos, mang_pos;
	const double MACH_TOL = 0.1; /* TOLERANCE */
	const double MIN_REL_RADIUS = 0.1; /* TOLERANCE */
	static BOND	*btmp = NULL;
	static POINT	*newp = NULL;
	static boolean	first = YES;

	debug_print("Mach_node","Entered find_Mach_node_states()\n");

	if (first)
	{
		first = NO;
		scalar(&btmp,sizeof(BOND));
		btmp->start = Static_point(fr->interf);
		btmp->end = Static_point(fr->interf);
		newp = Static_point(fr->interf);
	}

	/*	Find states on both sides of the new incident curve */

	oldp = oldn->posn;
	point_propagate(fr,(POINTER)wave,oldp,newp,
		Bond_at_node_of_o_curve(oldcinc),oldcinc->curve,dt,V);
	if (curve_ang_oriented_l_to_r(i_to_a_dir,oldcinc->orient))
	{
		ft_assign(RP->state[1],left_state(newp),fr->sizest);
		ft_assign(RP->state[0],right_state(newp),fr->sizest);
	}
        else 
	{
		ft_assign(RP->state[0],left_state(newp),fr->sizest);
		ft_assign(RP->state[1],right_state(newp),fr->sizest);
	}

	p1 = pressure(RP->state[1]);
	p3 = compute_behind_state_pressure(oldcslip,oldcmach,RP,fr,wave);

	if (debugging("Mach_node")) 
	{
		if (curve_ang_oriented_l_to_r(i_to_a_dir,oldcinc->orient))
		{
			verbose_print_state("Old state0",
				Right_state_at_node_of_o_curve(oldcinc));
			verbose_print_state("Old state1",
				Left_state_at_node_of_o_curve(oldcinc));
		}
		else
		{
			verbose_print_state("Old state0",
				Left_state_at_node_of_o_curve(oldcinc));
			verbose_print_state("Old state1",
				Right_state_at_node_of_o_curve(oldcinc));
		}
		verbose_print_state("New state0",RP->state[0]);
		verbose_print_state("New state1",RP->state[1]);
		(void) printf("New behind state pressure = %g\n",p3);
	}

	if (p3 < p1)
	{
		(void) printf("WARNING in find_Mach_node_states(), "
			      "-- p3 < p1\n");
		return ERROR_NODE;
	}

		/* Special code for the degenerate case */

	if ((p3 - p1) < MACH_TOL*p1) 
	{
		status = degenerate_mach_node(oldcinc,newcinc,newcmach,
			oldp,&Pcenter,i_to_a_dir,crossbinc,crossbmach,
			pc,wave,fr,rp,dt,dt_frac,RP,flag);
		if (status != GOOD_NODE)
		{
			(void) printf("WARNING in find_Mach_node_states(), ");
			(void) printf(" degenerate_mach_node() failed\n");
		}
		debug_print("Mach_node","Left find_Mach_node_states()\n");
		return status;
	}

	if (!find_steady_ahead_speed(i_to_a_dir,RP->state[0],RP->state[1],
			RP->state[3],p3,&q0,&theta01,&theta12,&theta03))
	{
		if (debugging("Mach_node"))
		{
			(void) printf("WARNING in find_Mach_node_states(), ");
			(void) printf("find_steady_ahead_speed() failed\n");
		}
		debug_print("Mach_node","Left find_Mach_node_states()\n");
		return ERROR_NODE;
	}
	if (debugging("Mach_node"))
	{
		(void) printf("Turn angles from find_steady_ahead_speed()\n");
		print_angle("\ttheta01 =",theta01,"\n");
		print_angle("\ttheta12 =",theta12,"\n");
		print_angle("\ttheta03 =",theta03,"\n");
	}

	/* find node position and node velocity by intersecting the circle
	   with radius = node_speed*dt and the incident shock */

	r = q0 * dt;
	for (i = 0; i < dim; i++)
		Coords(&Pcenter)[i] = Coords(oldp)[i] + vel(i,RP->state[0])*dt;
	*crossbmach = Bond_at_node_of_o_curve(newcmach);

	/*
	*	Calculate the angle of the velocity of the node.
	*	For dt~0 we use the angle obtained using a very small dt.
	*/

	if (r < (MIN_REL_RADIUS * 
			bond_length(Bond_at_node_of_o_curve(oldcinc))))
	{
		double 	dt0;	/* temporary time step for dt~0 */

		dt0 = MIN_REL_RADIUS * 
			bond_length(Bond_at_node_of_o_curve(oldcinc))/q0;

		propagated_tangent_bond_at_node(btmp,oldcinc->curve,
				oldcinc->orient,fr,(POINTER)wave,dt0);

		if (!robust_cross_bond_circle(btmp,&Pcenter,sqr(q0*dt0),
				&tcr,pc))
		{
			(void) printf("WARNING in find_Mach_node_states(), ");
			(void) printf("robust_cross_bond_circle() failed\n");
			debug_print("Mach_node","Left find_Mach_node_states()\n");
			return ERROR_NODE;
		}

		/* find incident bond for calculating the incident angle */

		for (i = 0; i < dim; i++)
		{
			nod_v[i] = (Coords(pc)[i] - Coords(oldp)[i]) / dt0;
			Coords(pc)[i] = Coords(oldp)[i] + nod_v[i]*dt;
		}

		/* Need to actually set crossbinc for modifying Mach node. */

		*crossbinc = Bond_at_node_of_o_curve(newcinc);
	}
	else
	{
		status = crossing_of_a_propagated_curve_and_circle(oldcinc,
			newcinc,r,&Pcenter,pc,crossbinc,&tcr,fr,(POINTER)wave,
			rp,dt,dt_frac,flag);

		if (status != GOOD_NODE)
		{
			(void) printf("WARNING in find_Mach_node_states(), ");
			(void) printf("crossing_of_a_propagated_curve_and_circle()");
			(void) printf(" failed\n");
			debug_print("Mach_node","Left find_Mach_node_states()\n");
			return MODIFY_TIME_STEP_NODE;
		}

		for (i = 0; i < dim; i++)
			nod_v[i] = (Coords(pc)[i] - Coords(oldp)[i]) / dt;
	}

	status = velocity_satisfies_CFL(newn,dt,dt_frac,fr);
	if (status != GOOD_NODE) return status;

	if (i_to_a_dir == CLOCKWISE)
	{
		iang_pos = YES; rang_pos = NO; mang_pos = YES;
	}
	else
	{
		iang_pos = NO; rang_pos = YES; mang_pos = NO;
	}

	Check_return(
	    s_polar_3(RP->state[0],YES,p1,iang_pos,NO,nod_v,RP->state[1],
		      &inc_ang,&theta01),
	    find_Mach_node_states)

	Check_return(
	    s_polar_3(RP->state[1],YES,p3,rang_pos,YES,nod_v,RP->state[2],
		      &refl_ang,&theta12),
	    find_Mach_node_states)

	Check_return(
	    s_polar_3(RP->state[0],YES,p3,mang_pos,YES,nod_v,RP->state[3],
		      &Mach_ang,&theta03),
	    find_Mach_node_states)

	for (i = 0; i < dim; i++)
		V[i] = vel(i,RP->state[0]) - nod_v[i];
	theta0 = angle(V[0],V[1]);
	cont_ang = avg_angle_and_normalize(theta0 + theta01 + theta12,
			theta0 + theta03);

	RP->ang[0] = inc_ang;
	RP->ang[1] = refl_ang;
	RP->ang[2] = cont_ang;
	RP->ang[3] = Mach_ang;
	RP->ang_dir = Opposite_ang_dir(i_to_a_dir);

	if (debugging("Mach_node"))
	{
		double chi;		/* node velocity angle */
		double ang;

		chi = angle(nod_v[0],nod_v[1]);
		(void) printf("q0 = %g\n",q0);
		(void) printf("*pc = (%g, %g)\n",Coords(pc)[0],Coords(pc)[1]);
		(void) printf("nod_v = (%g, %g), ",nod_v[0],nod_v[1]);
		print_angle("node_velocity_angle =",chi,"\n");
		print_angle("\tinc_ang =",inc_ang,"\n");
		print_angle("\trefl_ang =",refl_ang,"\n");
		print_angle("\tcont_ang =",cont_ang,"\n");
		print_angle("\tMach_ang =",Mach_ang,"\n");
		(void) printf("thetas from s_polar_3()\n");
		print_angle("\ttheta01 =",theta01,"\n");
		print_angle("\ttheta12 =",theta12,"\n");
		print_angle("\ttheta03 =",theta03,"\n");
		ang = fabs(theta01 + theta12 - theta03);
		print_angle("\tfabs(theta01 + theta12 - theta03) =",ang,"\n");
		print_RP_node_states("States after find_Mach_node_states()",
			nod_v,RP,MACH_NODE);
	}

		/* Check for bifurcations */

	if ((fabs(inc_ang - refl_ang) < EPS) ||
			(fabs(inc_ang - refl_ang - PI) < EPS))
	{
		screen("ERROR in find_Mach_node_states(), ");
		screen("refl shock coincides with incident shock\n");
		return ERROR_NODE;
	}

	debug_print("Mach_node","Left find_Mach_node_states()\n");
	return GOOD_NODE;
}		/*end find_Mach_node_states*/


/*	Support functions for find_Mach_node_states() */

/*
*			compute_behind_state_pressure():
*
*	This function uses the behind states to compute the new behind
*	pressure for the Mach node, which is the only behind information
*	used in the propagation.  After trying numerous algorithms, this
*	one worked the best, and was coincidentally, the simplest.  We
*	simply solve two non-local riemann problems (via point_propagate())
*	across the Mach stem and across the slip line, and average the
*	two pressures of the states thus obtained.
*/

LOCAL double compute_behind_state_pressure(
	O_CURVE		*cslip,
	O_CURVE		*cmach,
	RP_DATA		*RP,
	Front		*fr,
	Wave		*wave)
{
	NORMAL_FUNCTION sav_normal;
	INTERFACE       *intfc = cmach->curve->interface;
	POINT		*oldp;
	double		V[MAXD];
	double		tmp1, tmp2;
	double		dt = fr->dt;
	ANGLE_DIRECTION	i_to_a_dir = Opposite_ang_dir(RP->ang_dir);
	static POINT	*newp = NULL;

	debug_print("cbsp","Entering compute_behind_state_pressure()\n");
	if (newp == NULL)
	    newp = Static_point(fr->interf);

	oldp = Node_of_o_curve(cmach)->posn;

	sav_normal = interface_normal_function(intfc);
	set_temp_mnode_normal(&interface_normal_function(intfc));
	point_propagate(fr,(POINTER)wave,oldp,newp,
	                Bond_at_node_of_o_curve(cmach),cmach->curve,dt,V);
	tmp1 = (curve_ang_oriented_l_to_r(i_to_a_dir,cmach->orient)) ?
		pressure(right_state(newp)) :
		pressure(left_state(newp));
	interface_normal_function(intfc) = sav_normal;

	if (cslip->curve != NULL)
	{
		point_propagate(fr,(POINTER)wave,oldp,newp,
				Bond_at_node_of_o_curve(cslip),
				cslip->curve,dt,V);
		if (curve_ang_oriented_l_to_r(i_to_a_dir,cslip->orient))
			tmp2 = pressure(right_state(newp));
		else
			tmp2 = pressure(left_state(newp));

		debug_print("cbsp","Leaving compute_behind_state_pressure()\n");
		return 0.5*(tmp1 + tmp2);
	}
	else
	{
		debug_print("cbsp","Leaving compute_behind_state_pressure()\n");
		return tmp1;
	}
}		/*end compute_behind_state_pressure*/

EXPORT	void	set_temp_mnode_normal(
	NORMAL_FUNCTION *nf)
{
	static const char *nname = "temp_mnode_normal";
	nf->_normal = temp_mnode_normal;
	nf->_normal_name = nname;
}		/*end set_temp_mnode_normal*/

/*ARGSUSED*/
LOCAL void temp_mnode_normal(
	POINT		*p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	*hs,
	double		*nor,
	Front		*front)
{
	CURVE		*c = Curve_of_hs(hs);
	double		ang;
	ORIENTATION	orient;

	orient = (separation(p,c->start->posn,front->rect_grid->dim) == 0.0) ?
		POSITIVE_ORIENTATION : NEGATIVE_ORIENTATION;

	ang = (is_scalar_wave(wave_type(c))) ?
		Rp_data(Node_of(c,orient))->ang[2] :
		Rp_data(Node_of(c,orient))->ang[3];

	if (orient == POSITIVE_ORIENTATION)
	{
		nor[0] = sin(ang);
		nor[1] = -cos(ang);
	}
	else
	{
		nor[0] = -sin(ang);
		nor[1] = cos(ang);
	}
	
}		/*end temp_mnode_normal*/

typedef struct	{
	double	min_M01_sq, min_M12_sq, min_M03_sq;
	double	alpha01, alpha12, alpha03;
	double	beta;
} TA_DIFF_PARAMS;


LOCAL double turn_angle(
	double		M_sq,
	double		alpha,
	double		min_M_sq,
	double		beta)
{
	double		tan_theta;

	tan_theta = alpha / (M_sq + beta - alpha) *
					sqrt((M_sq + beta)/min_M_sq - 1.);
	
	return atan(tan_theta);
}		/*end turn_angle*/

LOCAL boolean turn_angle_difference(
	double		M0sq,
	double		*fans,
	POINTER		parameters)
{
	double	min_M01_sq = ((TA_DIFF_PARAMS *) parameters)->min_M01_sq;
	double	min_M12_sq = ((TA_DIFF_PARAMS *) parameters)->min_M12_sq;
	double	min_M03_sq = ((TA_DIFF_PARAMS *) parameters)->min_M03_sq;
	double	alpha01 = ((TA_DIFF_PARAMS *) parameters)->alpha01;
	double	alpha12 = ((TA_DIFF_PARAMS *) parameters)->alpha12;
	double	alpha03 = ((TA_DIFF_PARAMS *) parameters)->alpha03;
	double	beta = ((TA_DIFF_PARAMS *) parameters)->beta;
	double	theta01, theta12, theta03;

	theta01 = turn_angle(M0sq,alpha01,min_M01_sq,0.0);
	theta12 = turn_angle(M0sq,alpha12,min_M12_sq,beta);
	theta03 = turn_angle(M0sq,alpha03,min_M03_sq,0.0);

	*fans = theta01 - theta12 - theta03;
	return FUNCTION_SUCCEEDED;
}		/*end turn_angle_difference*/


/*
*		find_steady_ahead_speed():
*
*	Drives the iteration to compute the magnitude of the ahead velocity
*	in the steady frame given the ahead state and behind state pressures.
*	If the ahead state is at rest, this will simply be the magnitude of
*	the node velocity.
*
*	Made exportable for use in girefl.c for computing various diagnostic
*	quantities.
*/

EXPORT int find_steady_ahead_speed(
	ANGLE_DIRECTION	i_to_a_dir,
	Locstate	state0,
	Locstate	state1,
	Locstate	state3,
	double		p3,
	double		*psteady_ahead_speed,
	double		*theta01,
	double		*theta12,
	double		*theta03)
{
	TA_DIFF_PARAMS	ta_diff_params;
	double		p0, p1, p2;	/* given data (p2 = p3) */
	double		rho0, c0_sq;	/* given data */
	double		rho1;
	double		rho3, c3_sq;
	double		alpha01;	/* turning angle eqn parameters */
	double		alpha03;
	double		min_M01_sq;
	double		min_M12_sq;
	double		min_M03_sq;
	double		beta;
	double		Ml_sq, Mr_sq;	/* interval limits */
	double		M0_sq;		/* ahead Mach num in steady frame */
	double		sign;		/* (+) for i_to_a_dir clockwise */

	debug_print("Mach_node","Entered find_steady_ahead_speed()\n");

	p0 = pressure(state0);
	p1 = pressure(state1);
	p2 = p3;

	state_w_pr_on_Hugoniot(state0,p1,state1,GAS_STATE);
	state_w_pr_on_Hugoniot(state0,p3,state3,GAS_STATE);

	rho0 = Dens(state0);
	rho1 = Dens(state1);
	rho3 = Dens(state3);

	c0_sq = sound_speed_squared(state0);
	c3_sq = sound_speed_squared(state3);

	if (debugging("Mach_node"))
	{
		(void) printf("p0 = %g\tp1 = %g\tp3 = %g\n",p0,p1,p3);
		(void) printf("rho0 = %g\trho1 = %g\tc0_sq = %g\n",rho0,rho1,c0_sq);
	}

	ta_diff_params.alpha01 = alpha01 = (p1 - p0) / (rho0 * c0_sq);
	ta_diff_params.alpha12 = (p2 - p1) / (rho1 * c0_sq);
	ta_diff_params.alpha03 = alpha03 = (p3 - p0) / (rho0 * c0_sq);
	ta_diff_params.beta= beta = (p0 - p1) * (1./rho0 + 1./rho1) / c0_sq;
	ta_diff_params.min_M01_sq = min_M01_sq =
		mass_flux_squared(p1,state0) / (sqr(rho0) * c0_sq);
	ta_diff_params.min_M12_sq =  min_M12_sq =
		mass_flux_squared(p2,state1) / (sqr(rho1) * c0_sq);
	ta_diff_params.min_M03_sq =  min_M03_sq =
		mass_flux_squared(p3,state0) / (sqr(rho0) * c0_sq);

	/* Find bounds for the interval */

	Ml_sq = max(min_M03_sq,(min_M12_sq - beta)) + EPS;
	Mr_sq = (c3_sq - ((p0 - p3) * (1./rho0 + 1./rho3))) / c0_sq - EPS;

	if (debugging("ta_diff"))
	{
		print_function_values(turn_angle_difference,
				      (POINTER) &ta_diff_params,
				      0.0,Ml_sq,Mr_sq,50,
				      "turn_angle_difference",stdout);
	}

	if (debugging("Mach_node"))
		(void) printf("Ml_sq = %g\tMr_sq = %g\n",Ml_sq,Mr_sq);

	if (find_root(turn_angle_difference,(POINTER) &ta_diff_params,
		      0.0,&M0_sq,Ml_sq,Mr_sq,EPS,EPS) == FUNCTION_FAILED)
	{
		if (debugging("ta_diff_vals"))
		{
			print_function_values(turn_angle_difference,
					      (POINTER) &ta_diff_params,
					      0.0,Ml_sq,Mr_sq,100,
					      "ta_diff",stdout);
		}

		if (debugging("Mach_node"))
		{
			(void) printf("WARNING in find_steady_ahead_speed(), ");
			(void) printf("no root for turn_angle_difference()\n");
		}
		debug_print("Mach_node","Left find_steady_ahead_speed()\n");
		return NO;
	}

	if (debugging("Mach_node"))
		(void) printf("M0_sq = %g\n",M0_sq);

	sign = (i_to_a_dir == CLOCKWISE) ? 1.0 : -1.0;
	*theta03 = sign * turn_angle(M0_sq,alpha03,min_M03_sq,0.0);
	*theta01 = sign * turn_angle(M0_sq,alpha01,min_M01_sq,0.0);
	*theta12 = *theta03 - *theta01;
	*psteady_ahead_speed = sqrt(M0_sq * c0_sq);

	debug_print("Mach_node","Left find_steady_ahead_speed()\n");
	return YES;
}		/*end find_steady_ahead_speed*/


/*
*			degenerate_mach_node():
*
*	This is the function to compute the states and angles for a 
*	degenerate mach node (one in which the reflected shock is very
*	weak).  We take the behind state to be sonic in the steady frame,
*	and use this with (118.12) to compute the tangential component of
*	the node velocity.  The inc and mach shocks are parallel, and
*	the contact is parallel to the outgoing streamline.  The refl
*	wave is taken to be at cont_ang +- A.
*/

LOCAL int degenerate_mach_node(
	O_CURVE		*oldcinc,
	O_CURVE		*newcinc,
	O_CURVE		*newcmach,
	POINT		*oldp,
	POINT		*pcenter,
	ANGLE_DIRECTION	i_to_a_dir,
	BOND		**crossbinc,
	BOND		**crossbmach,
	POINT		*pc,
	Wave		*wave,
	Front		*fr,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	RP_DATA		*RP,
	NODE_FLAG	flag)
{
	NODE		*newn = Node_of_o_curve(newcinc);
	double		*nod_v = Node_vel(newn);
	double 		sign;
	double		p1;
	double		inc_ang;
	double		refl_ang;
	double		Mach_ang;
	double		cont_ang;
	double		vel1[MAXD];	/*unsteady behind vel in shock coords*/
	double		i_t[MAXD];	/*incident tangent from node*/
	double		i_n[MAXD];	/*incident normal - high p to low*/
	double		q1[MAXD];	/*steady behind velocity*/
	double		c1sq;		/*sound speed squared of state1*/
	double		alpha;		/*tan comp of nod_v in shock coords*/
	double		A;		/*mach angle = asin(1/M)*/
	double		s;		/*shock speed*/
	double		r;
	double		tcr;
	int		i, dim = fr->rect_grid->dim;
	int		status;

	debug_print("Mach_node","Entered degenerate_mach_node()\n");

	sign = (i_to_a_dir == CLOCKWISE) ? 1.0 : -1.0;
	p1 = pressure(RP->state[1]);

	find_tangent_to_curve(Node_of_o_curve(oldcinc)->posn,
		Bond_at_node_of_o_curve(oldcinc),oldcinc->curve,
		oldcinc->orient,i_t,fr);
	i_n[0] = sign*i_t[1];	i_n[1] = -sign*i_t[0];
	if (!s_polar_4(BEHIND_PRESSURE,p1,&s,i_n,
			RP->state[0],RP->state[1],RP->stype))
	{
		if (debugging("Mach_node"))
		{
			(void) printf("WARNING in degenerate_mach_node() -- ");
			(void) printf("s_polar_4() failed.\n");
		}
		return ERROR_NODE;
	}
	for (i = 0; i < dim; i++)
	{
		q1[i] = vel(i,RP->state[1]);	/*temp storage*/
	}
	vel1[0] = scalar_product(i_n,q1,dim);
	vel1[1] = scalar_product(i_t,q1,dim);

	set_state(RP->state[2],RP->stype,RP->state[1]);
	set_state(RP->state[3],RP->stype,RP->state[1]);
	
	c1sq = sound_speed_squared(RP->state[1]);
	alpha = vel1[1] + sqrt(c1sq - sqr(vel1[0] - s));
	for (i = 0; i < dim; i++)
		nod_v[i] = s*i_n[i] + alpha*i_t[i];

	status = velocity_satisfies_CFL(newn,dt,dt_frac,fr);
	if (status != GOOD_NODE) return status;

	inc_ang = angle(i_t[0],i_t[1]);
	Mach_ang = normalized_angle(inc_ang + PI);

	for (i = 0; i < dim; i++)
		q1[i] -= nod_v[i];

	cont_ang = angle(q1[0],q1[1]);
	A = asin(sqrt(c1sq) / mag_vector(q1,dim));
	refl_ang = cont_ang - sign*A;

	RP->ang[0] = inc_ang;
	RP->ang[1] = refl_ang;
	RP->ang[2] = cont_ang;
	RP->ang[3] = Mach_ang;
	RP->ang_dir = Opposite_ang_dir(i_to_a_dir);

	r = mag_vector(nod_v,dim) * dt;
	for (i = 0; i < dim; i++)
	{
		Coords(pc)[i] = Coords(oldp)[i] + nod_v[i] * dt;
		Coords(pcenter)[i] = Coords(oldp)[i] +
				vel(i,RP->state[0])*dt;
	}

	status = crossing_of_a_propagated_curve_and_circle(oldcinc,
		newcinc,r,pcenter,pc,crossbinc,&tcr,fr,(POINTER)wave,
		rp,dt,dt_frac,flag);
	*crossbmach = Bond_at_node_of_o_curve(newcmach);

	debug_print("Mach_node","Left degenerate_mach_node()\n");
	return GOOD_NODE;
}		/*end degenerate_mach_node*/


/*	
*			find_cross_node_states():
*
*	This function performs the shock polar analysis needed to update
*	the states at a cross node, which consist of two incident shocks
*	colliding producing reflected shocks and a contact discontinuity.
*
*/


EXPORT int find_cross_node_states(
	double		*nod_v,
	RP_DATA		*RP,
	boolean		*is_refl_raref,
	boolean		is_theta1_pos)
{
	double		 theta0;
	double		 theta1, theta2, theta3, theta4; /* Flow turn angles */
	double		 p1, p2, p4; /* pressures */
	int 		         dim = Params(RP->state[1])->dim;
	int		         is_theta4_pos;
	boolean		         Cplus_w3, Cplus_w2;
	RIEMANN_SOLVER_WAVE_TYPE wtype1, wtype4;

	debug_print("cross_node","Entered find_cross_node_states()\n");
	if (debugging("cross_node"))
	{
		(void) printf("node velocity = <%g, %g>, magnitude = %g\n",
			nod_v[0],nod_v[1],mag_vector(nod_v,dim));
		(void) printf("is_theta1_pos = %s\n",
			      (is_theta1_pos) ? "YES" : "NO");
	}
	p1 = pressure(RP->state[1]);
	p4 = pressure(RP->state[4]);
	is_theta4_pos = Cplus_w2 = (is_theta1_pos) ? NO : YES;
	Cplus_w3 = is_theta1_pos;

	/* Find state behind incident shocks */

	if (!s_polar_3(RP->state[0],YES,p1,is_theta1_pos,NO,
		nod_v,RP->state[1],&RP->ang[0],&theta1)
				 ||
	    !s_polar_3(RP->state[0],YES,p4,is_theta4_pos,NO,
		nod_v,RP->state[4],&RP->ang[4],&theta4)) 
	{

		(void) printf("WARNING in find_cross_node_states(), ");
		(void) printf("s_polar_3() failed\n");
		debug_print("cross_node",
			"Left find_cross_node_states() ans = NO\n");
		return NO;
	}


	/* Find the pressure at the intersection of the two shock polars */

	if (p1 > p4) 
	{
	    if (!intersection_of_two_shock_polars(RP->state[4],
	    					     RP->state[1],nod_v,&p2,
	    					     NULL,NULL,Cplus_w3,
	    					     Cplus_w2,
	    					     &wtype4,&wtype1)) 
	    {
	        (void) printf("WARNING in find_cross_node_states(), ");
	        (void) printf("intersection_of_two_shock_polars() failed\n");
	        debug_print("cross_node","Left find_cross_node_states() ans = NO\n");
	        return NO;
	    }	
	}
	else 
	{
	    if (!intersection_of_two_shock_polars(RP->state[1],
			                             RP->state[4],nod_v,&p2,
			                             NULL,NULL,
						     Cplus_w2,Cplus_w3,
			                             &wtype1,&wtype4)) 
	    {
		(void) printf("WARNING in find_cross_node_states(), ");
		(void) printf("intersection_of_two_shock_polars() failed\n");
		debug_print("cross_node","Left find_cross_node_states() ans = NO\n");
		return NO;
	    }
	}

	*is_refl_raref = (p2 < p1 || p2 < p4) ? YES : NO;

	/* Find states behind reflected waves */

	if (p2 >= p1) 
	{
		if (!s_polar_3(RP->state[1],YES,p2,Cplus_w2,
			YES,nod_v,RP->state[2],&RP->ang[1],
			&theta2)) 
		{
			(void) printf("WARNING in find_cross_node_states(), ");
			(void) printf("s_polar_3() failed\n");
			debug_print("cross_node",
				"Left find_cross_node_states() ans = NO\n");
			return NO;
		}
	}
	else 
	{
		if (!prandtl_meyer_wave(RP->state[1],p2,
					   Cplus_w2,nod_v,RP->state[2],
					   &RP->ang[5],&RP->ang[6],&theta2)) 
		{

			(void) printf("WARNING in find_cross_node_states(), ");
			(void) printf("prandtl_meyer_wave() failed\n");
			debug_print("cross_node",
				"Left find_cross_node_states() ans = NO\n");
			return NO;
		}
	}

	if (p2 >= p4) 
	{
		if (!s_polar_3(RP->state[4],YES,p2,Cplus_w3,
			YES,nod_v,RP->state[3],&RP->ang[3],
			&theta3)) 
		{
			(void) printf("WARNING in find_cross_node_states(), ");
			(void) printf("s_polar_3() failed\n");
			debug_print("cross_node",
				"Left find_cross_node_states() ans = NO\n");
			return NO;
		}
	}
	else 
	{
		if (!prandtl_meyer_wave(RP->state[4],p2,
					   Cplus_w3,nod_v,RP->state[3],
					   &RP->ang[5],&RP->ang[6],&theta3)) 
		{

			(void) printf("WARNING in find_cross_node_states(), ");
			(void) printf("prandtl_meyer_wave() failed\n");
			debug_print("cross_node",
				"Left find_cross_node_states() ans = NO\n");
			return NO;
		}
	}

	debug_print("cross_node",
		"theta1 = %g, theta2 = %g, theta3 = %g, theta4 = %g\n",
		theta1,theta2,theta3,theta4);
	debug_print("cross_node","theta1 + theta2 - theta3 - theta4= %g\n",
		theta1 + theta2 - theta3 - theta4);


	/* Find contact turning angle */

	theta0 = angle(vel(0,RP->state[0]) - nod_v[0],
				vel(1,RP->state[0]) - nod_v[1]);
	RP->ang[2] = avg_angle_and_normalize(theta0 + theta1 + theta2,
				theta0 + theta3 + theta4);

	if (debugging("cross_node")) 
	{
		(void) printf("is_refl_raref = %s\n",
			(*is_refl_raref) ? "YES" : "NO");
		print_RP_node_states("States after find_cross_node_states()",
			nod_v,RP,CROSS_NODE);
	}
	debug_print("cross_node","Left find_cross_node_states(), ans = YES\n");
	return YES;
}		/*end find_cross_node_states*/


/*	
*			find_overtake_node_states():
*
*	This function performs the shock polar analysis needed to update
*	the states at a overtake node, which occurs when one shock overtakes
*	another from behind.
*
*/


EXPORT int find_overtake_node_states(
	double		*nod_v,
	RP_DATA		*RP,
	boolean		*is_refl_raref,
	boolean		is_theta1_pos)
{
	double		theta0;
	double		theta1, theta2, theta3, theta4; /* Flow turn angles */
	double		p1, p2, p6; /* pressures */
	int		is_theta2_pos;
	boolean		refl_Cplus_wave, transm_Cplus_wave;
	RIEMANN_SOLVER_WAVE_TYPE wtype0, wtype2;

	debug_print("overtake","Entered find_overtake_node_states()\n");
	if (debugging("overtake"))
	{
	    (void) printf("Node velocity = <%g, %g>\n",nod_v[0],nod_v[1]);
	    (void) printf("is_theta1_pos = %s\n",
	    	      (is_theta1_pos) ? "YES" : "NO");
	    (void) printf("Input states\n");
	    verbose_print_state("RP->state[0]",RP->state[0]);
	    (void) printf("M0 = %g\n",mach_number(RP->state[0],nod_v));
	    verbose_print_state("RP->state[1]",RP->state[1]);
	    (void) printf("M1 = %g\n",mach_number(RP->state[1],nod_v));
	    verbose_print_state("RP->state[2]",RP->state[2]);
	    (void) printf("M2 = %g\n",mach_number(RP->state[2],nod_v));
	}

	p1 = pressure(RP->state[1]);
	p2 = pressure(RP->state[2]);
	is_theta2_pos = transm_Cplus_wave = is_theta1_pos;
	refl_Cplus_wave = (is_theta1_pos) ? NO : YES;

	/* Find state behing incident shocks */

	if (!s_polar_3(RP->state[0],YES,p1,is_theta1_pos,NO,
	    nod_v,RP->state[1],&RP->ang[0],&theta1))
	{
	    (void) printf("WARNING in find_overtake_node_states(), ");
	    (void) printf("s_polar_3() failed\n");
	    debug_print("overtake","Left find_overtake_node_states()\n");
	    return NO;
	}

	if (!s_polar_3(RP->state[1],YES,p2,is_theta2_pos,NO,
	    nod_v,RP->state[2],&RP->ang[1],&theta2)) 
	{
	    (void) printf("WARNING in find_overtake_node_states(), ");
	    (void) printf("s_polar_3() failed\n");
	    debug_print("overtake","Left find_overtake_node_states()\n");
	    return NO;
	}
	copy_state(RP->state[3],RP->state[2]);

	if (debugging("overtake"))
	{
	    (void) printf("After s_polar_3()\n");
	    verbose_print_state("RP->state[1]",RP->state[1]);
	    (void) printf("M1 = %g\n",mach_number(RP->state[1],nod_v));
	    verbose_print_state("RP->state[2]",RP->state[2]);
	    (void) printf("M2 = %g\n",mach_number(RP->state[2],nod_v));
	}


	/* Find the pressure at the intersection of the two shock polars */

	if (fabs(p1 - p2) < EPS*max(p1,p2))
	{
	    p6 = 0.5*(p1 + p2);
	    wtype0 = SHOCK;
	    wtype2 = (p6 > p2) ? SHOCK : RAREFACTION;
	}
	else if (!intersection_of_two_shock_polars(RP->state[0],
	                                              RP->state[2],nod_v,&p6,
	                                              NULL,NULL,
	    				              transm_Cplus_wave,
	    				              refl_Cplus_wave,
	    				              &wtype0,&wtype2)) 
	{
	    (void) printf("WARNING in find_overtake_node_states(), "
	                  "intersection_of_two_shock_polars() failed\n");
	    debug_print("overtake","Left find_overtake_node_states()\n");
	    return NO;
	}	

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	*is_refl_raref = (p6 < p2*Raref_press(RP->state[2])) ? YES : NO;
#else /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	*is_refl_raref = (p6 < p2) ? YES : NO;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	/* Find states behind reflected waves */

	if (*is_refl_raref)
	{
	    if (!prandtl_meyer_wave(RP->state[2],p6,
	    			       refl_Cplus_wave,nod_v,RP->state[4],
	    			       &RP->ang[2], &RP->ang[4],&theta3)) 
	    {
	        (void) printf("WARNING in find_overtake_node_states(), ");
	        (void) printf("prandtl_meyer_wave() failed\n");
	        debug_print("overtake","Left find_overtake_node_states()\n");
	        return NO;
	    }
	    RP->ang[3] = 0.5*(RP->ang[2] + RP->ang[4]);
	}
	else
	{
	    if (!s_polar_3(RP->state[2],YES,p6,refl_Cplus_wave,YES,nod_v,
			      RP->state[4],&RP->ang[2],&theta3)) 
	    {
	    	(void) printf("WARNING in find_overtake_node_states(), ");
	    	(void) printf("s_polar_3() failed\n");
	    	debug_print("overtake","Left find_overtake_node_states()\n");
	    	return NO;
	    }
	    RP->ang[3] = RP->ang[4] = RP->ang[2];
	}
	copy_state(RP->state[5],RP->state[4]);

	/* Find state behind transmitted wave */

	if (!s_polar_3(RP->state[0],YES,p6,transm_Cplus_wave,
	    	          YES,nod_v,RP->state[6],&RP->ang[6],&theta4)) 
	{
	    (void) printf("WARNING in find_overtake_node_states(), ");
	    (void) printf("s_polar_3() failed\n");
	    debug_print("overtake","Left find_overtake_node_states()\n");
	    return NO;
	}

	if (debugging("overtake"))
	{
	    (void) printf("theta1 = %g, theta2 = %g,\n",theta1,theta2);
	    (void) printf("\ttheta3 = %g, theta4 = %g\n",theta3,theta4);
	    (void) printf("theta1 + theta2 + theta3 - theta4 = %g\n",
	    	theta1 + theta2 + theta3 - theta4);
	}


	/* Find contact turning angle */

	theta0 = angle(vel(0,RP->state[0]) - nod_v[0],
	    	vel(1,RP->state[0]) - nod_v[1]);
	RP->ang[5] = avg_angle_and_normalize(theta0 + theta1 + theta2 + theta3,
	    				     theta0 + theta4);

	if (debugging("overtake")) 
	{
	    print_RP_node_states("States after find_overtake_node_states()",
	    			 nod_v,RP,OVERTAKE_NODE);
	}
	debug_print("overtake","Left find_overtake_node_states()\n");
	return YES;
}		/*end find_overtake_node_states*/
#endif /* defined(TWOD) */
