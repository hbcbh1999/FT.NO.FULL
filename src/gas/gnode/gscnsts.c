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
*				gscnsts.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains high level routines for two dimensional shock polar analysis
*	of shock wave contact discontinuity interactions.
*
*
*		is_regular_diffraction_node()
*		find_transmission_node_states()
*/

#if defined(FULL_PHYSICS) && defined(TWOD)

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	boolean	find_state_behind_incident_shock(double,double*,RP_DATA*,
						 double*,int*,int,NODE_FLAG);
LOCAL	boolean	is_diffraction_at_vacuum(double,RP_DATA*);
LOCAL	boolean	transm_node_speed(double,double*,POINTER);
LOCAL	int	diffraction_to_vacuum(double,double*,RP_DATA*,boolean*,
                                      boolean,NODE_FLAG);
LOCAL	int	identify_irregular_diffraction(double**,RP_DATA*,RP_DATA*,
	                                       double*,double*,boolean,boolean,
					       boolean*,
					       RIEMANN_SOLVER_WAVE_TYPE*,
					       RIEMANN_SOLVER_WAVE_TYPE*,
					       NODE_FLAG);
LOCAL	int	subsonic_state_behind_incident_shock(double,double,RP_DATA*,
						     RP_DATA*,double*,int,
						     NODE_FLAG);
LOCAL	void	root_limit(RP_DATA*,double,double,double,double*,double*,int);
#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void	debug_reg_dfrctn(double*,double**,RP_DATA*);
#endif /* defined(DEBUG_NODE_PROPAGATE) */


/*	
*			is_regular_diffraction_node():
*
*	This function performs the shock polar analysis needed to project
*	a set of five states onto states which satisfy the constraint of
*	being states at a diffraction node, i.e. which satisfy a shock hitting
*	a contact discontinuity, producing a transmitted shock and
*	either a reflected shock or rarefaction wave. The projection is defined
*	as follows: The inflow states (0 and 6) are projected onto states with
*	rest frame velocities parallel to a tangent to the incoming contact.
*	(This is presently done before this routine is called.)  A shock polar
*	analysis (s_polar_3) computes the entire state 1 behind the incident
*	shock, using the pressure of this state 1 as input data. Then states
*	1 through 5 and all angles are computed from this information by shock
*	polar analysis also.
*
*/

EXPORT int is_regular_diffraction_node(
	double		*posn,		/*Position of the node            */
	double		*nod_v,         /*Estimated new node velocity     */
	double		*old_nod_v,     /*old node velocity               */
	double		**t,            /*t[0] = incident shock tangent   */
	                                /*t[1] = incidnet contact tangent */
	RP_DATA		*RP,		/*RP data at node                 */
	RP_DATA		*oldRP,		/*RP data at old node             */
	boolean		*is_reflected_shock,
	Front		*front,
	int		n_type,
	NODE_FLAG	flag)
{
	double	theta1, theta2, theta3; /* Flow turn angles */
	double	v0x, v0y, v6x, v6y;
	double	p1, pstar; /* pressures */
	double   sgn;
	int	status;
	boolean	refl_Cplus_wave, transm_Cplus_wave;
	int	s1_subsonic;
	int	i, dim = Params(RP->state[0])->dim;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("is_reg_dfrctn","Entered is_regular_diffraction_node()\n");
	if (debugging("is_reg_dfrctn"))
	{
	    (void) printf("Input data\n");
	    debug_reg_dfrctn(nod_v,t,RP);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	switch (RP->ang_dir)
	{
	case COUNTER_CLOCK:
	    transm_Cplus_wave = YES;
	    refl_Cplus_wave = NO;
	    break;
	case CLOCKWISE:
	    transm_Cplus_wave = NO;
	    refl_Cplus_wave = YES;
	    break;
	case ANGLE_DIRECTION_NOT_SET:
	default:
	    screen("ERROR in is_regular_diffraction_node(), "
	           "RP->ang_dir = %s not set\n",
		   angle_direction_name(RP->ang_dir));
	    clean_up(ERROR);
	    return ERROR_DIFFRACTION;
	}

	p1 = pressure(RP->state[1]);
	RP->M[0] = mach_number(RP->state[0],nod_v);
	if (RP->M[0] < 1.0)
	{
	    node_warning("is_regular_diffraction_node",
		         "state ahead of incident shock is subsonic","\n",flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    if (node_warnings_off(flag) != YES)
	        debug_reg_dfrctn(nod_v,t,RP);
	    debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return ERROR_DIFFRACTION;
	}
	RP->M[6] = mach_number(RP->state[6],nod_v);
	if ((RP->M[6] < 1.0) && (n_type != TOT_INT_REFL_NODE))
	{
	    node_warning("is_regular_diffraction_node",
	                  "state ahead of transmitted shock is subsonic","\n",
			  flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    if (node_warnings_off(flag) != YES)
	        debug_reg_dfrctn(nod_v,t,RP);
	    debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return ERROR_DIFFRACTION;
	}

	if (!find_state_behind_incident_shock(p1,nod_v,RP,&theta1,
						 &s1_subsonic,n_type,flag))
	{
	    node_warning("is_regular_diffraction_node",
	                  "find_state_behind_incident_shock failed","\n",flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    if (node_warnings_off(flag) != YES)
	        debug_reg_dfrctn(nod_v,t,RP);
	    debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return ERROR_DIFFRACTION;
	}

	/* Check for bifurcation due to transonic incident shock */

	if (s1_subsonic)
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    	(void) printf("Transsonic incident shock\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = subsonic_state_behind_incident_shock(p1,theta1,RP,oldRP,
							  old_nod_v,n_type,
							  flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    	debug_reg_dfrctn(nod_v,t,RP);
	    debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;

	}
	if (n_type == TOT_INT_REFL_NODE)
	{
	    pstar = pressure(RP->state[6]);
	    *is_reflected_shock = NO;
	    set_state(RP->state[5],RP->stype,RP->state[6]);
	    RP->M[5] = RP->M[6];
	    theta3 = 0.0;
	    RP->ang[5] = ERROR_FLOAT;/*Undefined*/
	}
	else
	{
	    RIEMANN_SOLVER_WAVE_TYPE	wtype6, wtype1;
	
	/* Check for possible vacuum in state6 and */
	/* Find the pressure at the intersection of the two shock polars */

	    if (is_diffraction_at_vacuum(p1,RP))
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	    (void) printf("diffraction at vacuum\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    	status = diffraction_to_vacuum(theta1,nod_v,RP,
	    				       is_reflected_shock,
					       refl_Cplus_wave,flag);

#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	    debug_reg_dfrctn(nod_v,t,RP);
	    	debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	return status;
	    }


	    if (!intersection_of_two_shock_polars(RP->state[6],
	    					     RP->state[1],nod_v,&pstar,
	    					     NULL,NULL,
						     transm_Cplus_wave,
	    					     refl_Cplus_wave,
						     &wtype6,&wtype1)) 
	    {
	    	status = identify_irregular_diffraction(t,RP,oldRP,
							nod_v,&pstar,
							transm_Cplus_wave,
							refl_Cplus_wave,
							is_reflected_shock,
	    						&wtype6,&wtype1,
							flag);
	        if (status != REGULAR_DIFFRACTION)
		{
#if defined(DEBUG_NODE_PROPAGATE)
	    	    if (debugging("is_reg_dfrctn"))
	    	    {
	    	        (void) printf("No intersection of shock polars\n");
	    	        debug_reg_dfrctn(nod_v,t,RP);
	    	    }
	    	    debug_print("is_reg_dfrctn",
			  "Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	    return status;
		}
	    }


	    *is_reflected_shock = (pstar >= p1) ? YES : NO;

	    /* Find state behind transmitted shock */

	    if (!s_polar_3(RP->state[6],YES,pstar,transm_Cplus_wave,YES,
			      nod_v,RP->state[5],&RP->ang[5],&theta3)) 
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	        if (debugging("is_reg_dfrctn"))
	        {
		    node_warning("is_regular_diffraction_node",
	    	                  "s_polar_3() failed for transmitted shock",
				  "\n",flag);
	            if (node_warnings_off(flag) != YES)
		    {
	                (void) printf("pstar = %g\n",pstar);
	                debug_reg_dfrctn(nod_v,t,RP);
		    }
	        }
	        debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	        return ERROR_DIFFRACTION;
	    }

	    /* Check for sonic transition at transmitted shock */

	    RP->M[5] = mach_number(RP->state[5],nod_v);
	    if ((use_subsonic_state(flag) != YES) && (RP->M[5] < SONIC_MINUS))
	    {
	    	if (oldRP != NULL)
	    	{
	    	    /* Numerical overshoot can cause an incorrect
	    	     * value of *is_reflected_shock.  We use the
	    	     * oldRP to find the correct value.
	    	     */
	    	    if (pressure(oldRP->state[1]) > pressure(oldRP->state[4]))
	    		*is_reflected_shock = NO;
	    	}
	    	status = (*is_reflected_shock) ?
	    		    PRECURSOR_WITH_REFLECTED_SHOCK :
	    		    PRECURSOR_WITH_REFLECTED_RAREFACTION;
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	{
		    char s[80];
	    	    (void) sprintf(s,"state[5] is subsonic, M = %g",RP->M[5]);
		    node_warning("is_regular_diffraction_node",s,"\n",flag);
	            if (node_warnings_off(flag) != YES)
	    	        debug_reg_dfrctn(nod_v,t,RP);
	    	}
	    	debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	return status;
	    }
	}

	/* Find state behind reflected wave */

	if (*is_reflected_shock) 
	{	 /* reflected wave is a shock */

	    if (!s_polar_3(RP->state[1],YES,pstar,refl_Cplus_wave,YES,
			      nod_v,RP->state[4],&RP->ang[2],&theta2)) 
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	{
		    node_warning("is_regular_diffraction_node",
	    	                 "s_polar_3() failed for reflected shock","\n",
				 flag);
	            if (node_warnings_off(flag) != YES)
	    	        debug_reg_dfrctn(nod_v,t,RP);
	    	}
	    	debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	return ERROR_DIFFRACTION;
	    }
	    RP->ang[1] = RP->ang[3] = RP->ang[2];

	    /* Check for sonic transition at reflected shock */

	    RP->M[4] = mach_number(RP->state[4],nod_v);
	    if (RP->M[4] < SONIC_MINUS)
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	{
	    	    (void) printf("regular to mach diffraction\n");
	    	    debug_reg_dfrctn(nod_v,t,RP);
	    	}
	    	debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	return REGULAR_TO_MACH_DIFFRACTION;
	    }

	}
	else 
	{	/* Reflected wave is a rarefaction */

	    if (!prandtl_meyer_wave(RP->state[1],pstar,refl_Cplus_wave,
	    			       nod_v,RP->state[4],&RP->ang[1],
				       &RP->ang[3],&theta2)) 
	    {

#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	{
		    node_warning("is_regular_diffraction_node",
	    	                 "prandtl_meyer_wave() failed "
				 "for reflected rarefaction","\n",flag);
	            if (node_warnings_off(flag) != YES)
	    	        debug_reg_dfrctn(nod_v,t,RP);
	    	}
	    	debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	return ERROR_DIFFRACTION;
	    }
	    RP->M[4] = mach_number(RP->state[4],nod_v);
	    RP->ang[2] = RP->ang[1];
	}
	set_state(RP->state[3],RP->stype,RP->state[4]);
	RP->M[3] = RP->M[4];
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn"))
	{
	    (void) printf("is_reflected_shock = %s\n",
	    			(*is_reflected_shock) ? "YES" : "NO");
	    (void) printf("theta1 = %g, theta2 = %g, theta3 = %g\n",
	    		theta1,theta2,theta3);
	    (void) printf("theta1 + theta2 - theta3 = %g\n",
	    		theta1 + theta2 - theta3);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	/* Find contact angles */

	v0x = vel(0,RP->state[0]) - nod_v[0];
	v0y = vel(1,RP->state[0]) - nod_v[1];
	v6x = vel(0,RP->state[6]) - nod_v[0];
	v6y = vel(1,RP->state[6]) - nod_v[1];
	RP->ang[6] = avg_angle_and_normalize(angle(-v0x,-v0y),angle(-v6x,-v6y));

	RP->ang[4] = avg_angle_and_normalize(angle(v0x,v0y) + theta1 + theta2,
	    	                             angle(v6x,v6y) + theta3);

	RP->theta[0] = theta1;
	RP->theta[1] = 0.0;
	RP->theta[2] = theta2;
	RP->theta[3] = 0.0;
	RP->theta[4] = 0.0;
	RP->theta[5] = theta3;
	RP->theta[6] = 0.0;

	/* Check for monotonicity of RP->ang */

	sgn = (RP->ang_dir == COUNTER_CLOCK) ? 1.0 : -1.0;
	for (i = 0; i < 7; ++i)
	{
	    if (sgn*sin(RP->ang[(i+1)%7] - RP->ang[i]) < 0.0)
	    {
	        (void) printf("WARNING in is_regular_diffraction_node(), "
		              "inconsistent angle directions\n");
		(void) printf("RP->ang[%d] = %g (%g degrees)\n",i,
		               RP->ang[i],degrees(RP->ang[i]));
		(void) printf("RP->ang[%d] = %g (%g degrees)\n",(i+1)%7,
		               RP->ang[(i+1)%7],degrees(RP->ang[(i+1)%7]));
	    }
	}



#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn")) 
	{
	    print_RP_node_states(
	    	"States after is_regular_diffraction_node()",
	    	nod_v,RP,n_type);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	if ((front != NULL) && (Apply_CFL_at_nodes(front) == YES))
	{
	    double	*h = front->rect_grid->h;

	    for (i = 0; i < dim; i++)
	    {
	    	set_max_front_speed(i,fabs(nod_v[i]),return_obst_state(),
	    			    posn,front);
	    }
	    set_max_front_speed(dim,scaled_hypot(nod_v,h,dim),
	    		        return_obst_state(),posn,front);
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
	    	(void) printf("Node velocity added to ");
	    	(void) printf("max front wave speed\n");
	    	print_general_vector("nod_v = ",nod_v,dim,"\n");
	    	print_general_vector("Spfr(front) = ",
	    		Spfr(front),dim+1,"\n");
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	}
#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("is_reg_dfrctn","Left is_regular_diffraction_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return REGULAR_DIFFRACTION;
}		/*end is_regular_diffraction_node*/


LOCAL	boolean find_state_behind_incident_shock(
	double	  p1,
	double	  *nod_v,
	RP_DATA	  *RP,
	double	  *theta1,
	int	  *s1_subsonic,
	int	  n_type,
	NODE_FLAG flag)
{
	double		M0sq, M1, psonic;
	int		is_theta1_pos;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("is_reg_dfrctn","Entered find_state_behind_incident_shock()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	M0sq = sqr(RP->M[0]);
	is_theta1_pos = (RP->ang_dir == COUNTER_CLOCK) ? YES : NO;
	if (n_type == TOT_INT_REFL_NODE)
	{
	    psonic = pressure_at_sonic_point(M0sq,RP->state[0]);
	    if (p1 > psonic)
		p1 = psonic;
	}
	if (p1 >= max_behind_shock_pr(M0sq,RP->state[0]))
	{
	    *s1_subsonic = YES;
	    *theta1 = 0.0;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    	(void) printf("p1 > maximum behind steady shock pressure\n");
	    debug_print("is_reg_dfrctn","Left find_state_behind_incident_shock()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		return YES;
	}

	if (!s_polar_3(RP->state[0],YES,p1,is_theta1_pos,NO,nod_v,
			  RP->state[1],&RP->ang[0],theta1)) 
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		node_warning("find_state_behind_incident_shock",
	    	             "s_polar_3() failed","\n",flag);
	    }
	    debug_print("is_reg_dfrctn","Left find_state_behind_incident_shock()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return NO;
	}
	set_state(RP->state[2],RP->stype,RP->state[1]);
	M1 = RP->M[1] = RP->M[2] = mach_number(RP->state[1],nod_v);
	*s1_subsonic =  (M1 < SONIC_MINUS) ? YES : NO;
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn"))
	{
	    verbose_print_state("state0 ",RP->state[0]);
	    verbose_print_state("state1 ",RP->state[1]);
	    (void) printf("M1 = %g\n",M1);
	}
	debug_print("is_reg_dfrctn","Left find_state_behind_incident_shock()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return YES;
}		/*end find_state_behind_incident_shock*/

/*ARGSUSED*/
LOCAL	int subsonic_state_behind_incident_shock(
	double	  p1,
	double	  theta1,
	RP_DATA	  *RP,
	RP_DATA	  *oldRP,
	double	  *old_nod_v,
	int	  n_type,
	NODE_FLAG flag)
{
	double		M6sq, pmax, pmax_t_ang, max_t_ang, t_ang;
	int		status;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("is_reg_dfrctn",
	      "Entered subsonic_state_behind_incident_shock()\n");
	if (debugging("is_reg_dfrctn"))
	{
	    if (oldRP != NULL)
	    {
	    	print_RP_node_states("Old RP Data",old_nod_v,oldRP,n_type);
	    }
	    else
	    	(void) printf("oldRP = NULL\n");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	
	/* Subsonic state behind incident shock 	*/
	/* check for possible anomalous reflection	*/

	if (oldRP != NULL)
	{
	    if (pressure(oldRP->state[1]) < pressure(oldRP->state[4]))
	    {
	    	/* Previous reflected wave was a shock */
	    	status = REGULAR_TO_MACH_DIFFRACTION;
	    }
	    else
	    {
	    	/* Previous reflected wave was a rarefaction */
	    	status = ANOMALOUS_REFLECTION;
	    }
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("is_reg_dfrctn",
		  "Left subsonic_state_behind_incident_shock()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	}

	/*
	*  oldRP data not available,  try to determine bifurcation
	*  from current timestep data.
	*/

	M6sq = sqr(RP->M[6]);
	pmax = max_behind_shock_pr(M6sq,RP->state[6]);

	if (p1 >= pmax)
	    status = ANOMALOUS_REFLECTION;
	else if (steady_state_wave_curve(p1,M6sq,&t_ang,RP->state[6]) ==
		 FUNCTION_FAILED)
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		char s[80];
	        (void) sprintf(s,"steady_state_wave_curve() failed "
	                         "at pressure p1 = %g, M6sq = %g",p1,M6sq);
		node_warning("subsonic_state_behind_incident_shock",s,
			     "\n",flag);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = ERROR_DIFFRACTION;
	}
	else if (fabs(theta1) <= fabs(t_ang))
	{
	    status = ANOMALOUS_REFLECTION;
	}
	else if (pr_at_max_turn_angle(&pmax_t_ang,M6sq,RP->state[6]) ==
		 FUNCTION_FAILED)
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		char s[80];
	        (void) sprintf(s,"pr_at_max_turn_angle() failed, "
	                      "p1 = %g, pr_max_theta0 = %g M6sq = %g",
			      p1,pmax_t_ang,M6sq);
		node_warning("subsonic_state_behind_incident_shock",s,
			     "\n",flag);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = ERROR_DIFFRACTION;
	}
	else if (p1 <= pmax_t_ang)
	{
	    status = REGULAR_TO_MACH_REFLECTION;
	}

	else if (steady_state_wave_curve(pmax_t_ang,M6sq,&max_t_ang,
					 RP->state[6]) == FUNCTION_FAILED)
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		char s[80];
	        (void) sprintf(s,"steady_state_wave_curve() failed at "
	                         "pmax_t_ang, p1 = %g, pr_max_theta0 = %g "
				 "M6sq = %g",p1,pmax_t_ang,M6sq);
		node_warning("subsonic_state_behind_incident_shock",s,
			     "\n",flag);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = ERROR_DIFFRACTION;
	}
	else if (fabs(theta1) > fabs(max_t_ang))
	{
	    status = REGULAR_TO_MACH_REFLECTION;
	}
	else
	{

	    /*
	    *  At this point it is impossible to determine whether
	    *  the diffraction is a transition from regular to anomalous
	    *  diffraction or to Mach diffraction without recourse
	    *  to previous time step data.  Furthermore this is a very
	    *  unlikely case since transitions occur at sonic points
	    *  and should usually put you in one of the cases above.
	    *  Thus at this point I will return an error message
	    *  and try to solve this problem at a later date.
	    */


#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		char s[80];
	    	(void) sprintf(s,"unable to identify bifurcation, "
	    	                 "p1 = %g, pr_max_theta0 = %g M6sq = %g",
	    		         p1,pmax_t_ang,M6sq);
		node_warning("subsonic_state_behind_incident_shock",s,
			     "\n",flag);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    status = ERROR_DIFFRACTION;
	}
#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("is_reg_dfrctn",
		"Left subsonic_state_behind_incident_shock()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return status;
}		/*end subsonic_state_behind_incident_shock*/

/*ARGSUSED*/
LOCAL	int diffraction_to_vacuum(
	double	  theta1,
	double	  *nod_v,
	RP_DATA	  *RP,
	boolean	  *is_reflected_shock,
	boolean	  refl_Cplus_wave,
	NODE_FLAG flag)
{
	double		theta2, pstar;

	pstar = pressure(RP->state[6]);
	*is_reflected_shock = NO;
	set_state(RP->state[5],RP->stype,RP->state[6]);
	RP->M[5] = RP->M[6];

	/* Reflected wave is a rarefaction */

	if (!prandtl_meyer_wave(RP->state[1],pstar,refl_Cplus_wave,nod_v,
				   RP->state[4],&RP->ang[1],&RP->ang[3],
				   &theta2)) 
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		node_warning("diffraction_to_vacuum",
	    	             "prandtl_meyer_wave() failed","\n",flag);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return ERROR_DIFFRACTION;
	}
	RP->ang[2] = RP->ang[1];
	set_state(RP->state[3],RP->stype,RP->state[4]);
	RP->M[4] = RP->M[3] = mach_number(RP->state[4],nod_v);
	RP->ang[4] = RP->ang[5] = RP->ang[2];

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn"))
	{
	    double theta3;

	    theta3 = theta1 + theta2;
	    (void) printf("After diffraction_to_vacuum(), ");
	    (void) printf("is_reflected_shock = %s\n",
	    		  (*is_reflected_shock) ? "YES" : "NO");
	    (void) printf("theta1 = %g, theta2 = %g, theta3 = %g\n",
	    	          theta1,theta2,theta3);
	    (void) printf("theta1 + theta2 - theta3 = %g\n",
	    	          theta1 + theta2 - theta3);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	/* Find contact angles */

	RP->ang[6] = angle(nod_v[0] - vel(0,RP->state[0]),
			nod_v[1] - vel(1,RP->state[0]));

	RP->theta[0] = theta1;
	RP->theta[1] = 0.0;
	RP->theta[2] = theta2;
	RP->theta[3] = 0.0;
	RP->theta[4] = 0.0;
	RP->theta[5] = theta1 + theta2;
	RP->theta[6] = 0.0;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn")) 
	{
	    print_RP_node_states("States after diffraction_to_vacuum()",
			         nod_v,RP,DIFFRACTION_NODE);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return REGULAR_DIFFRACTION;
}		/*end diffraction_to_vacuum*/

/*ARGSUSED*/
LOCAL	int identify_irregular_diffraction(
	double	                 **t,
	RP_DATA	                 *RP,
	RP_DATA	                 *oldRP,
	double                    *nod_v,
	double                    *pstar,
	boolean                     transm_Cplus_wave,
	boolean                     refl_Cplus_wave,
	boolean	                 *is_reflected_shock,
	RIEMANN_SOLVER_WAVE_TYPE *wtype6,
	RIEMANN_SOLVER_WAVE_TYPE *wtype1,
	NODE_FLAG                flag)
{
	double                    M1sq, M6sq;
	double                    p6mta, p1mta;
	double                    p, nv[3];
	RIEMANN_SOLVER_WAVE_TYPE wt1, wt6;
	int                      status = ERROR_DIFFRACTION;
#if defined(DEBUG_NODE_PROPAGATE)
	static const char *mesg = "identify_irregular_diffraction() returns ";
#endif /* defined(DEBUG_NODE_PROPAGATE) */

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("is_reg_dfrctn","Entered identify_irregular_diffraction()\n");
	if (debugging("is_reg_dfrctn"))
	{
	    print_rsoln_wave("wtype1 = ",*wtype1,"\n");
	    print_rsoln_wave("wtype6 = ",*wtype6,"\n");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (is_small_inc_ang_reg_diff_node(t,RP,is_reflected_shock,flag))
	{
	    status = REGULAR_DIFFRACTION;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		node_warning("identify_irregular_diffraction",
	                      "Answer computed by small incident angle "
	                      "approximation","\n",flag);
	        print_diffraction_status(mesg,status);
	    }
	    debug_print("is_reg_dfrctn","Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	}

	/* Test for robustness with respect to increase in node velocity */
	nv[0] = 1.1*nod_v[0];/*TOLERANCE*/
	nv[1] = 1.1*nod_v[1];/*TOLERANCE*/
	if (intersection_of_two_shock_polars(RP->state[6],RP->state[1],
					     nv,&p,NULL,NULL,
				             transm_Cplus_wave,refl_Cplus_wave,
					     &wt6,&wt1)) 
	{
	    *wtype6 = wt6;
	    *wtype1 = wt1;
	    *pstar = p;
	    nod_v[0] = nv[0];
	    nod_v[1] = nv[1];
	    status = REGULAR_DIFFRACTION;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		node_warning("identify_irregular_diffraction",
	                     "Answer computed by increasing node velocity",
			     "\n",flag);
	        if (node_warnings_off(flag) != YES)
	            print_diffraction_status(mesg,status);
	    }
	    debug_print("is_reg_dfrctn","Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	}

	if (*wtype1 == RAREFACTION)
	{
	    *is_reflected_shock = NO;
	    status = PRECURSOR_WITH_REFLECTED_RAREFACTION;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    	print_diffraction_status(mesg,status);
	    debug_print("is_reg_dfrctn","Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	}

	M1sq = sqr(RP->M[1]);
	if (pr_at_max_turn_angle(&p1mta,M1sq,RP->state[1]) == FUNCTION_FAILED)
	{
	    status = ERROR_DIFFRACTION;
	    *is_reflected_shock = NO;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		node_warning("identify_irregular_diffraction",
	                      "pr_at_max_turn_angle() failed at state 1","\n",
			      flag);
	        if (node_warnings_off(flag) != YES)
	            print_diffraction_status(mesg,status);
	    }
	    debug_print("is_reg_dfrctn","Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	}
	M6sq = sqr(RP->M[6]);
	if (pr_at_max_turn_angle(&p6mta,M6sq,RP->state[6]) == FUNCTION_FAILED)
	{
	    status = ERROR_DIFFRACTION;
	    *is_reflected_shock = NO;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
	    	node_warning("identify_irregular_diffraction",
	    	              "pr_at_max_turn_angle() failed at state 6","\n",
			      flag);
	        if (node_warnings_off(flag) != YES)
	    	    print_diffraction_status(mesg,status);
	    }
	    debug_print("is_reg_dfrctn","Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	}

	if (*wtype1 == SHOCK)
	{
	    *is_reflected_shock = YES;
	    status = (p6mta < p1mta) ? PRECURSOR_WITH_REFLECTED_SHOCK :
	    			       REGULAR_TO_MACH_DIFFRACTION;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    	print_diffraction_status(mesg,status);
	    debug_print("is_reg_dfrctn","Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	}

	if (oldRP != NULL)
	{
	    if (pressure(oldRP->state[1]) < pressure(oldRP->state[4]))
	    {
	    	/* Previous reflected wave was a shock */

	    	status = PRECURSOR_WITH_REFLECTED_SHOCK;
	    	*is_reflected_shock = YES;
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	    print_diffraction_status(mesg,status);
	        debug_print("is_reg_dfrctn",
		      "Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		return status;
	    }
	    else
	    {
	    	/* Previous reflected wave was a rarefaction */
	    	status = PRECURSOR_WITH_REFLECTED_RAREFACTION;
	    	*is_reflected_shock = NO;
#if defined(DEBUG_NODE_PROPAGATE)
	    	if (debugging("is_reg_dfrctn"))
	    	    print_diffraction_status(mesg,status);
		debug_print("is_reg_dfrctn",
		      "Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		return status;
	    }
	}

	status = ERROR_DIFFRACTION;
	*is_reflected_shock = NO;
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn"))
	{
	    node_warning("identify_irregular_diffraction",
	                 "unable to identify bifurcation,  CODE NEEDED","\n",
			 flag);
	    if (node_warnings_off(flag) != YES)
	        print_diffraction_status(mesg,status);
	}
	debug_print("is_reg_dfrctn","Left identify_irregular_diffraction()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return status;
}		/*end identify_irregular_diffraction*/


#if defined(DEBUG_NODE_PROPAGATE)
LOCAL void debug_reg_dfrctn(
	double		*nod_v,
	double		**t,
	RP_DATA		*RP)
{
	double		v0[MAXD], v6[MAXD];
	double		vp[MAXD], sp;
	double		q0, q1, q12, q6;
	double		i0, i1;
	double		inc_ang;
	int		i, dim = Params(RP->state[0])->dim;

	(void) printf("DATA for is_regular_diffraction_node()\n");

	print_general_vector("node velocity = ",nod_v,dim,"\n");
	(void) printf("is_theta1_pos = %s\n",
		      (RP->ang_dir == COUNTER_CLOCK) ? "YES" : "NO");
	(void) printf("States\n");
	verbose_print_state("state0 ",RP->state[0]);
	verbose_print_state("state1 ",RP->state[1]);
	verbose_print_state("state6 ",RP->state[6]);
	for (i = 0; i < dim; i++)
	{
	    v0[i] = vel(i,RP->state[0]) - nod_v[i];
	    v6[i] = vel(i,RP->state[6]) - nod_v[i];
	}
	print_general_vector("rest frame vel0 ",v0,dim,"");
	print_general_vector(", vel6 ",v6,dim,"");
	(void) vector_product(v0,v6,vp,dim);
	(void) printf("cross prod %g\n",vp[0]);

	q0 = mag_vector(v0,dim);
	q6 = mag_vector(v6,dim);
	i0 = specific_enthalpy(RP->state[0]);
	i1 = specific_enthalpy(RP->state[1]);

	q12 = q0*q0 + 2.0*(i0 - i1);
	q1 = (q12 > 0.0) ? sqrt(q12) : -sqrt(fabs(q12));
	(void) printf("Steady state speeds ahead q0 = %g, q1 = %g%s, q6 = %g\n",
		      q0,fabs(q1),(q1 < 0.0) ? "i" : "",q6);
	(void) printf("Mach numbers - M0 = %g, M1 = %g%s, M6 = %g\n",
		      q0/sound_speed(RP->state[0]),
		      fabs(q1)/sound_speed(RP->state[1]),(q1 < 0.0) ? "i" : "",
		      q6/sound_speed(RP->state[6]));
	print_general_vector("Wave directions t0 = ",t[0],dim,"");
	print_general_vector(", t1 = ",t[1],dim,"\n");
	sp = scalar_product(t[0],t[1],dim);
	vp[0] = vector_product(t[0],t[1],vp,dim);
	inc_ang = angle(sp,vp[0]);
	if (inc_ang > PI)
	    inc_ang -= 2.0*PI;
	print_angle("Incident angle = ",inc_ang,"\n");
}		/*end debug_reg_dfrctn*/
#endif /* defined(DEBUG_NODE_PROPAGATE) */

LOCAL	boolean is_diffraction_at_vacuum(
	double		p1,
	RP_DATA		*RP)
{
	double M6sq, prmx;
	const double eps = MACH_EPS;

	M6sq = sqr(RP->M[6]);
	prmx = max_behind_shock_pr(M6sq,RP->state[6]);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (prmx < Min_pressure(RP->state[6]))
	    return YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	if (prmx < eps*p1) /*TOLERANCE*/
	    return YES;
	return NO;
}		/*end is_diffraction_at_vacuum*/

/*
*			is_small_inc_ang_reg_diff_node():
*
*	Sets up the states and angles for a regular diffraction node
*	with very small incident angle (and hence large node velocity).
*	RP->state[0] and RP->state[4] are assumed to be given.  RP->state[1]
*	is given as input and then modified to be that state with
*	the same pressure and connected to RP->state[0] by an oblique
*	shock with normal <nx[1], ny[1]>.  States 2 and 3 are then found
*	by solving a one dimensional Riemann problem between states 1 and 4
*	in the direction given by <nx[0], ny[0]>.
*/

EXPORT int is_small_inc_ang_reg_diff_node(
	double	  **t,
	RP_DATA	  *RP,
	boolean	  *is_reflected_shock,
	NODE_FLAG flag)
{
	double		         p0, p1;
	double		         pstarl, pstarr;
	double		         ustarl, ustarr;
	double		         ml, mr;
	double		         shock_speed, inc_ang;
	double		         q65sqr, q01sqr, q14sqr;
	double		         csqr1, csqr4, tau0, tau1;
	double		         m, a0, i1, i4;
	double		         tan_ia, tan_iasqr, tmp1, tmp2;
	double		         M1sqr, M4sqr;
	double		         A1, A4, turn_ang;
	double		         nor[MAXD];
	double		         vtan;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	int		         i;
	static Locstate	         Tsl = NULL, Tsr = NULL;
	static double	         **n = NULL;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("is_reg_dfrctn","Entered is_small_inc_ang_reg_diff_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	if (Tsr == NULL)
	{
	    Gas_param *params = Params(RP->state[1]);
	    (*params->_alloc_state)(&Tsl,sizeof(VGas));
	    (*params->_alloc_state)(&Tsr,sizeof(VGas));
	    bi_array(&n,2,MAXD,FLOAT);
	}

		/* Calculate incident angle */

	for (i = 0; i < 2; i++)
	{
	    RP->ang[6*i] = angle(t[i][0],t[i][1]);
	    if (RP->ang_dir == COUNTER_CLOCK)
	    {
	    	n[i][0] = -t[i][1];
	    	n[i][1] = t[i][0];
	    }
	    else
	    {
	    	n[i][0] = t[i][1];
	    	n[i][1] = -t[i][0];
	    }
	}
	
	inc_ang = angle(t[0][0]*t[1][0] + t[0][1]*t[1][1],
				t[1][0]*t[0][1] - t[1][1]*t[0][0]);
	p1 = pressure(RP->state[1]);
	m = mass_flux(p1,RP->state[0]);
	a0 = acoustic_impedance(RP->state[0]);
	if (inc_ang > PI) inc_ang -= 2.*PI;
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn"))
	    (void) printf("Incident angle = %g\n",inc_ang);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	if (fabs(sin(inc_ang)) > m/a0)
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("is_reg_dfrctn"))
	    {
		node_warning("is_small_inc_ang_reg_diff_node",
	                     "Ahead state in not supersonic","\n",flag);
	    }
	    debug_print("is_reg_dfrctn","Left is_small_inc_ang_reg_diff_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return NO;
	}

		/* Find state behind incident shock */

	(void) s_polar_4(BEHIND_PRESSURE,p1,&shock_speed,n[0],RP->state[0],
		         RP->state[1],RP->stype);


		/* Find pressure behind transmitted wave */

	set_state_for_find_mid_state(Tsl,RP->state[6]);
	set_state_for_find_mid_state(Tsr,RP->state[1]);
	Vel(Tsl)[0] = n[1][0]*Vel(Tsl)[0] + n[1][1]*Vel(Tsl)[1];
	Vel(Tsr)[0] = n[1][0]*Vel(Tsr)[0] + n[1][1]*Vel(Tsr)[1];
	(void) find_mid_state(Tsl,Tsr,0.0,&pstarl,&pstarr,&ustarl,&ustarr,
			      &ml,&mr,&l_wave,&r_wave);

	if (l_wave != SHOCK)
	{
	    node_warning("is_small_inc_ang_reg_diff_node",
	                 "Unexpected case - no transmitted shock","\n",flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("is_reg_dfrctn","Left is_small_inc_ang_reg_diff_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return NO;
	}

		/* Find state behind transmitted wave */

	q01sqr = sqr(m) / sqr(Dens(RP->state[0]));
	q65sqr = mass_flux_squared(pstarl,Tsl)/sqr(Dens(Tsl));
	RP->ang[5] = normalized_angle(RP->ang[6] + PI + 
			              inc_ang*sqrt(q65sqr/q01sqr));
	if (RP->ang_dir == COUNTER_CLOCK)
	{
	    nor[0] = sin(RP->ang[5]);
	    nor[1] = -cos(RP->ang[5]);
	}
	else
	{
	    nor[0] = -sin(RP->ang[5]);
	    nor[1] = cos(RP->ang[5]);
	}
	(void) s_polar_4(BEHIND_PRESSURE,pstarl,&shock_speed,nor,
			 RP->state[6],RP->state[5],RP->stype);

		/* Find state behind reflected wave */

	if (r_wave == SHOCK)
	{
	    *is_reflected_shock = YES;
	    q14sqr = mass_flux_squared(pstarr,Tsr)/sqr(Dens(Tsr));
	    RP->ang[2] = normalized_angle(RP->ang[6] + PI -
			inc_ang*(sqrt(q14sqr/q01sqr) -
				 1.0 + Dens(RP->state[0])/Dens(RP->state[1])));
	    RP->ang[1] = RP->ang[3] = RP->ang[2];
	    if (RP->ang_dir == COUNTER_CLOCK)
	    {
	    	nor[0] = sin(RP->ang[2]);
	    	nor[1] = -cos(RP->ang[2]);
	    }
	    else
	    {
	    	nor[0] = -sin(RP->ang[2]);
	    	nor[1] = cos(RP->ang[2]);
	    }
	    (void) s_polar_4(BEHIND_PRESSURE,pstarr,&shock_speed,nor,
			     RP->state[1],RP->state[4],RP->stype);
	}
	else
	{
	    *is_reflected_shock = NO;
	    state_on_adiabat_with_pr(Tsr,pstarr,RP->state[4],TGAS_STATE);
	    vtan = t[0][0]*Vel(Tsr)[0] + t[0][1]*Vel(Tsr)[1];
	    Vel(RP->state[4])[0] = t[0][0]*vtan + n[0][0]*ustarr;
	    Vel(RP->state[4])[1] = t[0][1]*vtan + n[0][1]*ustarr;
	    csqr4 = sound_speed_squared(RP->state[4]);
	    i4 = specific_enthalpy(RP->state[4]);
	    set_state(RP->state[4],RP->stype,RP->state[4]);
	    csqr1 = sound_speed_squared(RP->state[1]);
	    i1 = specific_enthalpy(RP->state[1]);
	    p0 = pressure(RP->state[0]);
	    tau0 = 1.0/Dens(RP->state[0]);
	    tau1 = 1.0/Dens(RP->state[1]);
	    tan_ia = tan(inc_ang);
	    tan_iasqr = sqr(tan_ia);
	    tmp1 = (tau0 - tau1)*tan_iasqr;
	    tmp2 = (p1 - p0)*(sqr(tau1)*tan_iasqr + sqr(tau0));
	    M1sqr = csqr1*tmp1/tmp2;
	    M4sqr = csqr4*tmp1/(tmp2 + 2.0*(i4 - i1)*tmp1);
	    if ((M1sqr < 1.0) || (M4sqr < 1.0))
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	debug_print("is_reg_dfrctn",
		      "Left is_small_inc_ang_reg_diff_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		return NO;
	    }
	    A1 = asin(1.0/sqrt(M1sqr));
	    A4 = asin(1.0/sqrt(M4sqr));

	    if (steady_state_wave_curve(pstarr,M1sqr,&turn_ang,RP->state[1]) ==
		FUNCTION_FAILED)
	    {
#if defined(DEBUG_NODE_PROPAGATE)
	    	debug_print("is_reg_dfrctn",
		      "Left is_small_inc_ang_reg_diff_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		return NO;
	    }
	    if (RP->ang_dir == COUNTER_CLOCK) 
	    {
	    	RP->ang[1] = normalized_angle(-A1 + RP->ang[6]);
	    	RP->ang[3] = normalized_angle(turn_ang - A4 + RP->ang[6]);
	    }
	    else
	    {
	    	RP->ang[1] = normalized_angle(A1 + RP->ang[6]);
	    	RP->ang[3] = normalized_angle(A4 - turn_ang + RP->ang[6]);
	    }
	    RP->ang[2] = RP->ang[1];
	}
	RP->ang[4] = normalized_angle(PI + RP->ang[0]);
	set_state(RP->state[2],RP->stype,RP->state[1]);
	set_state(RP->state[3],RP->stype,RP->state[4]);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("is_reg_dfrctn")) 
	{
	    int	i;
	    char	mesg[20];

	    (void) printf("States after is_small_inc_ang_reg_diff_node()\n");
	    for (i = 0; i < 7; i++)
	    {
	    	(void) sprintf(mesg,"RP->state[%d]",i);
	    	verbose_print_state(mesg,RP->state[i]);
	    }

	    (void) printf("\nAngles:\n");
	    print_angle("Incident angle = ",RP->ang[0],"\n");
	    if (*is_reflected_shock) 
	    {
	    	print_angle("Reflected angle = ", RP->ang[2],"\n");
	    }
	    else 
	    {
	    	print_angle("Rarefaction angle0 = ", RP->ang[1],"\n");
	    	print_angle("Rarefaction angle1 = ", RP->ang[3],"\n");
	    }
	    print_angle("Back contact angle = ",RP->ang[4],"\n");
	    print_angle("Transmitted angle = ",RP->ang[5],"\n");
	    print_angle("Front contact angle = %g\n",RP->ang[6],"\n");
	    print_angle("Transmitted angle minus Contact angle = ",
			RP->ang[5] - RP->ang[4],"\n");
	}
	debug_print("is_reg_dfrctn","Left is_small_inc_ang_reg_diff_node()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return YES;
}		/*end is_small_inc_ang_reg_diff_node*/


/*
*			find_transmission_node_states():
*
* 		This routine calculates the states around a transmission 
*	node, where the incident shock impinges on a contact, and causes a
* 	transmitted shock but no reflected shock.
*
*
*		        state 1 = state 2         /                         
*		                              inc/                        
*		                                /
*		                               /
*		                              /     state 0
*		                             / 
*		                            / 
*		           ----------------/----------        
*		  contact(behind)        /       contact(ahead) 
*		                       /    
*		                     /       state 4
*		        state 3    /    
*		                 / transmitted
*
*		
*		We assume that regions 0 and 4 lie in front of the shock,
*	i.e. gas particles cross the shocks from 0 and 4 into 1 and 2.  In
*	order to move to the steady frame, we need the node velocity.  This is
*	determined in the following way: We take eqn. (136.01) across the
*	incident and the transmitted shock, and since the turning angles are
*	equal, we may set sqr(tan(theta01)) = sqr(tan(theta43)).  See
*	transm_node_speed() below -- this is exactly the turn angle relation
*	just mentioned with some manipulation.  To derive it, simply take
*	the equation below, plug in the expressions in tn_prms, and it is
*	relatively simple to recover the turn angle relation above.
*	As input parameters we use rho0, rho4, p0, p1 and either the slip of 
*	the velocity across the contact discontinuity ahead or the incident 
*	angle, which is the angle the incident shock makes with the contact
*	ahead. In the USE_SLIP case, we use the velocities ahead to compute
*	the slip, giving a relation between q0 and q4, and a sixth order
*	polynomial for q0.  For USE_INCIDENT_ANGLE, we get an explicit
*	expression for q0 and a quadratic polynomial for q4 squared.
*	These polynomials will have in general two positive roots.  However
*	experimental evidence suggests that it is always the weaker of these
*	two roots which is observed.  Since the flow ahead in the steady
*	frame is parallel to the contact ahead, this gives the steady frame
*	velocities ahead.  The node velocity is then found by the difference
*	between these steady velocities and the lab frame velocities which
*	are assumed to be given.  The behind states are then found by using
*	the standard oblique shock conditions.
*/

typedef struct {
	double	dp0;	/* dp0 = rho0*(pb - pa)/mi^2 */
	double	dp4;	/* dp4 = rho4*(pb - pa)/mt^2 */
	double	A;	/* A = rho4*(q0 - q4)/mt */
	double	B;	/* B = rho4*mi/(rho0*mt) */
	double	C;	/* C = sqr(rho4*mi^2/(rho0*mt^2)) */
} TN_PARAMS;

LOCAL boolean transm_node_speed(
	double		x,
	double		*speed,
	POINTER		parameters)
{
	double		A = ((TN_PARAMS *) parameters)->A;
	double		B = ((TN_PARAMS *) parameters)->B;
	double		C = ((TN_PARAMS *) parameters)->C;
	double		dp0 = ((TN_PARAMS *) parameters)->dp0;
	double		dp4 = ((TN_PARAMS *) parameters)->dp4;
	double		y;

	y = B*x - A;
	*speed = (x*x - 1.0)*sqr(y*y - dp4) - C*(y*y - 1.0)*sqr(x*x - dp0);
	return FUNCTION_SUCCEEDED;
}		/*end transm_node_speed*/

LOCAL	void	root_limit(
	RP_DATA		*RP,
	double		a,
	double		mi,
	double		mt,
	double		*xmin,
	double		*xmax,
	int		transonic_transmitted_shock)
{
	double		V0 = 1.0/Dens(RP->state[0]);
	double		V1 = 1.0/Dens(RP->state[1]);
	double		V3 = 1.0/Dens(RP->state[3]);
	double		V4 = 1.0/Dens(RP->state[4]);
	double		c1s = sound_speed_squared(RP->state[1]);
	double		c3s = sound_speed_squared(RP->state[3]);
	double		xm, xM;
	double		si_sonic, st_sonic;

	si_sonic = sqrt(1.0 + (c1s - sqr(mi*V1))/sqr(V0*mi));
	st_sonic = (a + sqrt(sqr(V4*mt) + (c3s - sqr(mt*V3))))/(V0*mi);
	if (transonic_transmitted_shock)
	{
		xm = a/(V0*mi);
		xm = max(1.0,xm);
		xM = min(si_sonic,st_sonic);
	}
	else
	{
		xm = max(1.0,st_sonic);
		xM = si_sonic;
	}
	if (xm > xM)
	{
		*xmin = 1.0;
		*xmax = si_sonic;
	}
	else
	{
		*xmin = xm;
		*xmax = xM;
	}
	return;
}		/*end root_limit*/

EXPORT int find_transmission_node_states(
	double	  *q,
	double	  **t,
	RP_DATA	  *RP,
	int	  trans_node_param,
	int	  root_type,
	NODE_FLAG flag)
{
	TN_PARAMS	tn_prms;
	double		rho0;	  	/* rho ahead of incident shock */
	double		rho1;	  	/* rho behind incident shock */
	double		rho3;	  	/* rho behind transmitted shock */
	double		rho4;	  	/* rho ahead of transmitted shock */
	double		pa, pb;    	/* pressure ahead, behind */
	double		dp;		/* pb - pa */
	double		q0_mins;	/* min squared for ahead vel */
	double		q4_mins;	/* min squared for behind vel */
	double		mi, mt;		/* mass fluxes across inc/transm */
	double		x;		/* q0 = x*mi/rho0 */
	double		theta01;	/* turn angle across incident */
	double		theta43;	/* turn angle across transmitted */
	double      	cost, sint;	/* cos/sin of average of turn angles */
	double      	tants;		/* tan(theta01) squared*/
	double		r_ratio_a;	/* rho1/rho0 */
	double		r_ratio_b;	/* rho3/rho4 */
	double		sinbeta1s;	/*squares of incident sine and cosine*/
	double		cosbeta1s;
	double		sinbeta3s;	/*squares of transm sine and cosine*/
	double		cosbeta3s;
	double		incident_angle;	/* these names are misleading */	
	double		contact_angle;	/* look at the usage below */
	double		transm_angle;
	double		a, slip[MAXD];	/* the slip */
	double		node_v[MAXD];	/* node velocity*/
	double		v0[MAXD];	/* velocities in steady frame */	
	double		v1[MAXD];
	double		v3[MAXD];
	double		v4[MAXD];
	double		discrm;
	double		eps, delta;
	const double     meps = MACH_EPS;
	double	     	xmin, xmax;
	double	      	length;
	double		sgn;
	double		sin_ang, cos_ang;
	int		i, dim;
	
#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("transmission","Entered find_transmission_node_states()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	dim = Params(RP->state[0])->dim;
	length = mag_vector(t[1],dim);
	for (i = 0; i < dim; i++) t[1][i] /= length;
	RP->ang[0] = angle(t[1][0],t[1][1]);
	sgn = (RP->ang_dir == CLOCKWISE) ? -1.0 : 1.0;

		/*Find state information*/

	pb = pressure(RP->state[1]);
	pa = pressure(RP->state[0]);
	rho0 = Dens(RP->state[0]);
	rho4 = Dens(RP->state[4]);
	dp = pb - pa;
	state_w_pr_on_Hugoniot(RP->state[0],pb,RP->state[1],RP->stype);
	mi = mass_flux(pb,RP->state[0]);
	rho1 = Dens(RP->state[1]);
	state_w_pr_on_Hugoniot(RP->state[4],pb,RP->state[3],RP->stype);
	mt = mass_flux(pb,RP->state[4]);
	rho3 = Dens(RP->state[3]);
	r_ratio_a = rho1/rho0;
	r_ratio_b = rho3/rho4;

		/* find minimum Mach number */

	q0_mins = sqr(mi)/sqr(rho0);
	q4_mins = sqr(mt)/sqr(rho4);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("transmission"))
	{
	    verbose_print_state("State0:",RP->state[0]);
	    verbose_print_state("State4:",RP->state[4]);
	    (void) printf("\tpa = %g\n",pa);
	    (void) printf("\tpb = %g\n",pb);
	    (void) printf("Normal velocity(0) squared = %g\n",q0_mins);
	    (void) printf("Normal velocity(4) squared = %g\n",q4_mins);
	    if (trans_node_param == USE_INCIDENT_ANGLE)
	    {
	    	(void) printf("trans_node_param == USE_INCIDENT_ANGLE\n");
	    	print_general_vector("Direction of incident shock:  ",
	    			     t[0],dim,"\n");
	    }
	    else
	    {
	    	(void) printf("trans_node_param == USE_SLIP\n");
	    	(void) printf("Direction of incident shock not given\n");
	    }
	    print_general_vector("Direction of front contact:  ",
							t[1],dim,"\n");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	switch (trans_node_param) 
	{
	case USE_SLIP:
	{
     	    /*
	    * We define the slip vector as the difference between the
	    * velocities in states 0 and 4.  This vector is parallel to
	    * the direction of the front contact in the steady frame.
	    * The slip a is the dot product of the slip vector and the
	    * unit vector in the direction of the contact ahead. 
            */

	    for (i = 0; i < dim; i++)
	    {
		    slip[i] = vel(i,RP->state[0]) - vel(i,RP->state[4]);
	    }

	    a = -scalar_product(t[1],slip,dim);
	    tn_prms.A = rho4*a/mt;
	    tn_prms.B = rho4*mi/(rho0*mt);
	    tn_prms.C = rho4*mi*mi/(rho0*mt*mt);
	    tn_prms.C = sqr(tn_prms.C);
	    tn_prms.dp0 = rho0*dp/(mi*mi);
	    tn_prms.dp4 = rho4*dp/(mt*mt);
		
	    root_limit(RP,a,mi,mt,&xmin,&xmax,NO);

	    delta = (xmax - xmin)*EPS; /* Scale convergence criteria */
	    (void) transm_node_speed(0.5*(xmin+xmax),&eps,(POINTER)&tn_prms);
	    eps = fabs(eps)*EPS;
	    eps = max(eps, meps);

	    if (debugging("transm_node_speed"))
	    {
		    print_function_values(transm_node_speed,(POINTER)&tn_prms,
				0.0,xmin,xmax,100,
				"transm_node_speed",stdout);
	    }

	    if (find_root(transm_node_speed,(POINTER) &tn_prms,
			  0.0,&x,xmin,xmax,eps,delta) == FUNCTION_FAILED)
	    {
		node_warning("find_transmission_node_states",
		             "can't find solution with supersonic "
			     "transmitted shock","\n",flag);
		root_limit(RP,a,mi,mt,&xmin,&xmax,YES);
		delta = (xmax - xmin)*EPS;
		(void) transm_node_speed(0.5*(xmin + xmax),&eps,
					 (POINTER) &tn_prms);
		eps = fabs(eps)*EPS;
		eps = max(eps, meps);
		if (find_root(transm_node_speed,(POINTER) &tn_prms,0.0,&x,xmin,
			      xmax,eps,delta) == FUNCTION_FAILED)
		{
		    if (debugging("transmission"))
		    {
			node_warning("find_transmission_node_states",
			             "can't find steady flow speed","\n",flag);
		    }
		    debug_print("transmission",
		    	  "Left find_transmission_node_states()\n");
		        return BIFURCATION_TRANSMISSION;
		}
	    }
	    q[0] = x*mi/rho0;

		/* Compute steady speed in region 4 */

	    q[4] = q[0] - a;

		/* Compute new incident angle and new transmitted angle */

	    sinbeta1s = q0_mins/sqr(q[0]);
	    cosbeta1s = 1.0 - sinbeta1s;
	    sinbeta3s = q4_mins/sqr(q[4]);
	    cosbeta3s = 1.0 - sinbeta3s;

#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("transmission"))
	    {
		(void) printf("USE_SLIP in find_transmission_node_states()\n");
		for (i = 0; i < dim; i++)
			v0[i] = vel(i,RP->state[0]);
		print_general_vector("Velocity(0):  ",v0,dim,"\n");
		for (i = 0; i < dim; i++)
			v4[i] = vel(i,RP->state[4]);
		print_general_vector("Velocity(4):  ",v4,dim,"\n");
		print_general_vector("Slip uni_array:  ",slip,dim,"\n");
		(void) printf("The slip a = %g\n",a);
		(void) printf("xmin = %g\n",xmin);
		(void) printf("xmax = %g\n",xmax);
		(void) printf("From root finder, q0 = %g\n",q[0]);
		(void) printf("q4 = q0 - a = %g\n",q[4]);
		(void) printf("\tsinbeta1s = %g\n",sinbeta1s);
		(void) printf("\tcosbeta1s = %g\n",cosbeta1s);
		(void) printf("\tsinbeta3s = %g\n",sinbeta3s);
		(void) printf("\tcosbeta3s = %g\n",cosbeta3s);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    if((cosbeta1s < 0.0) || (cosbeta3s < 0.0))
	    {
		if (debugging("transmission"))
		{
		    node_warning("find_transmission_node_states",
		                 "Pressure increase incompatible with shock",
				 "\n",flag);
		}
		debug_print("transmission",
			  "Left find_transmission_node_states()\n");
		return ERROR_TRANSMISSION;
	    }

	    incident_angle = sgn*angle(sqrt(cosbeta1s),sqrt(sinbeta1s));
	    transm_angle = sgn*angle(sqrt(cosbeta3s),sqrt(sinbeta3s));
	    break;
	}
	case USE_INCIDENT_ANGLE: 
	{
	    (void) vector_product(t[1],t[0],&sin_ang,dim);
	    cos_ang = scalar_product(t[1],t[0],dim);
	    length = hypot(sin_ang,cos_ang);
	    sin_ang /= length;
	    cos_ang /= length;

	    sinbeta1s = sqr(sin_ang);
	    cosbeta1s = sqr(cos_ang);
	    incident_angle = atan2(sin_ang,cos_ang);

	    /* Compute steady velocity above right*/

	    q[0] = sqrt(q0_mins/sinbeta1s);

	    /* Compute steady velocity below right */
		
	    tants = sinbeta1s*cosbeta1s/sqr(cosbeta1s + 1./(r_ratio_a-1.));
	    discrm = 1. - 4.*tants*r_ratio_b/sqr(r_ratio_b - 1.);
	    if (discrm < 0.0) 
	    {
	        if (debugging("transmission"))
		{
		    node_warning("find_transmission_node_states",
		                 "Transmission node doesn't exist","\n",flag);
		}
		debug_print("transmission","Left find_transmission_node_states()\n");
		return ERROR_TRANSMISSION;
	    }
	    discrm = sqrt(discrm);
	    if (root_type == WEAK) 
	    {
	        q[4] = sqrt(q4_mins*(1.0 + (1.0-discrm)/
			                   (r_ratio_b*(1.0+discrm))));
	    }
	    else if (root_type == STRONG) 
	    {
	        q[4] = sqrt(q4_mins*(1.0 + (1.0+discrm)/
					   (r_ratio_b*(1.0-discrm))));
	    }
	    else 
	    {
	        if (debugging("transmission"))
	        {
		    node_warning("find_transmission_node_states",
	                         "Unknown root type.","\n",flag);
	        }
	        debug_print("transmission","Left find_transmission_node_states()\n");
	        return ERROR_TRANSMISSION;
	    }

	    /* Compute new transmitted angle */

	    sinbeta3s = q4_mins/sqr(q[4]);
	    cosbeta3s = 1. - sinbeta3s;

	    if (cosbeta3s < 0.0)
	    {
	        if (debugging("transmission"))
	        {
		    node_warning("find_transmission_node_states",
	    	                 "Pressure increase incompatible with shock.",
				 "\n",flag);
	        }
	        debug_print("transmission","Left find_transmission_node_states()\n");
		return ERROR_TRANSMISSION;
	    }

	    transm_angle = sgn*angle(sqrt(cosbeta3s),sqrt(sinbeta3s));
	    break;
	}
	default:
	    screen("ERROR in find_transmission_node_states(), "
	           "Unknown parameter.\n");
	    clean_up(ERROR);
	}

		/* Steady velocities in states 0 and 4 are parallel 
		*  to the contact ahead.
		*/
	
	for (i = 0; i < dim; i++)
	{
	    v0[i] = -q[0]*t[1][i];
	    v4[i] = -q[4]*t[1][i];
	}

		/* Compute the steady speeds behind */

	q[1] = sqrt(q[0]*q[0] - ((1./rho1)+(1./rho0))*dp);
	q[2] = q[1];
	q[3] = sqrt(q[4]*q[4] - ((1./rho3)+(1./rho4))*dp);

		/* Compute turning angles */

	theta01 = sgn*atan((dp/(rho0*q[0]*q[0]-dp))*
			sqrt((q[0]*q[0]/q0_mins)-1.));
	theta43 = sgn*atan((dp/(rho4*q[4]*q[4]-dp))*
			sqrt((q[4]*q[4]/q4_mins)-1.));

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("transmission"))
	{
	    print_angle("\ttheta01 =",theta01,"\n");
	    print_angle("\ttheta43 =",theta43,"\n");
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

		/* If all is well we should have theta01 = theta43. */

 	if (fabs(theta01-theta43) > 5.0*(fabs(theta01)+fabs(theta43))*EPS) /*TOLERANCE*/
	{
	    if (debugging("transmission"))
	    {
		node_warning("find_transmission_node_states",
	                     "Turning angles for T node are unequal","\n",flag);
	    }
	    debug_print("transmission","Left find_transmission_node_states()\n");
	    return ERROR_TRANSMISSION;
	}

		/* Compute steady velocities behind */

	cost = 0.5*(cos(theta01) + cos(theta43));
	sint = 0.5*(sin(theta01) + sin(theta43));
	v1[0] = q[1]*(-t[1][0]*cost + t[1][1]*sint);
	v1[1] = -q[1]*(t[1][0]*sint + t[1][1]*cost);
	v3[0] = q[3]*(-t[1][0]*cost + t[1][1]*sint);
	v3[1] = -q[3]*(t[1][0]*sint + t[1][1]*cost);

		/* Translate velocities from the steady frame to lab frame */

	switch (trans_node_param) 
	{
	case USE_SLIP:
	    for (i = 0; i < dim; i++)
	    {
	    	node_v[i] = 0.5*(vel(i,RP->state[0]) +
				 vel(i,RP->state[4]) - v0[i] - v4[i]);
	    }
	    break;

	case USE_INCIDENT_ANGLE: 
	    Dens(RP->state[4]) = rho4;
	    Press(RP->state[4]) = pa;
	    for (i = 0; i < dim; i++)
	    {
	    	node_v[i] = vel(i,RP->state[0]) - v0[i];
	    	Vel(RP->state[4])[i] = (v4[i] + node_v[i]);
	    }
	    set_type_of_state(RP->state[4],TGAS_STATE);
	    set_state(RP->state[4],RP->stype,RP->state[4]);
	    break;

	default:
	    screen("ERROR in find_transmission_node_states(), "
	           "Unknown parameter.\n");
	    clean_up(ERROR);
	}

	Dens(RP->state[1]) = rho1;
	Press(RP->state[1]) = pb;
	for (i = 0; i < dim; i++)
	    Vel(RP->state[1])[i] = (v1[i] + node_v[i]);
	set_type_of_state(RP->state[1],TGAS_STATE);
	set_state(RP->state[1],RP->stype,RP->state[1]);
	set_state(RP->state[2],RP->stype,RP->state[1]);

	Dens(RP->state[3]) = rho3;
	Press(RP->state[3]) = pb;
	for (i = 0; i < dim; i++)
	    Vel(RP->state[3])[i] = (v3[i] + node_v[i]);
	set_type_of_state(RP->state[3],TGAS_STATE);
	set_state(RP->state[3],RP->stype,RP->state[3]);

		/* Compute the contact angle. */

	contact_angle = .5*(theta01 + theta43);

	RP->ang[1] = RP->ang[0] + incident_angle;
	RP->ang[4] = PI + transm_angle + RP->ang[0];
	RP->ang[3] = PI + contact_angle + RP->ang[0];
	RP->ang[2] = RP->ang[1];
	for (i = 0; i < 5; i++)
	    RP->ang[i] = normalized_angle(RP->ang[i]);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("transmission")) 
	{
	    print_RP_node_states("States after find_transmission_node_states()",
			         node_v,RP,TRANSMISSION_NODE);
	}

	debug_print("transmission","Left find_transmission_node_states()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return REGULAR_TRANSMISSION;
}		/*end find_transmission_node_states*/

#endif /* defined(FULL_PHYSICS) && defined(TWOD) */
