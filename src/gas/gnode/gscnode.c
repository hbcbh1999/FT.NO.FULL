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


#if defined(FULL_PHYSICS) && defined(TWOD)
/*
*
*				gscnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Contains the routine used in the first advance of nodes
*	corresponding to shock wave-contact discontinuity interactions.
*
*/


#include <gdecs/gdecs.h>

typedef struct {
	Front *fr;
	POINT *p0, *pc;
	BOND *b, *bc;
	CURVE *c;
	Locstate sa, sa_opp;
	double dt;
	double abs_v0[MAXD];
	double abs_v[MAXD];
	double tcr;
	ORIENTATION c_orient;
	SIDE inc_side;
} OVERSHOOT;

	/* LOCAL Function Declarations */
LOCAL	boolean	amn(double,double*,POINTER);
LOCAL	int	correct_new_node_posn_for_overshoot(NODE*,double*,double*,
						    BOND**,BOND**,POINT*,
						    O_CURVE*,O_CURVE*,O_CURVE*,
						    O_CURVE*,RP_DATA*,Front*,
						    Wave*,double,NODE_FLAG);
LOCAL	int	identify_curves_at_shock_diffraction(NODE*,NODE*,O_CURVE**,
						     O_CURVE**,
						     ANGLE_DIRECTION*,boolean*,
						     Front*,NODE*,NODE_FLAG);
LOCAL	int	init_transmission_node_deprecursion(NODE*,NODE*,O_CURVE*,
						    O_CURVE*,O_CURVE*,O_CURVE*,
						    double,double,Front*,
						    Wave*,RPROBLEM**);
LOCAL	int	interacting_nodes(NODE*,NODE*,O_CURVE**,O_CURVE**,double,
				  RPROBLEM**,Front*,Wave*);
LOCAL	int	is_short_and_opp_unprop(O_CURVE*,O_CURVE*,RECT_GRID*);
LOCAL	int	modify_transmission_node(POINT*,BOND*,BOND*,NODE*,NODE*,
					 O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
					 O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
					 ANGLE_DIRECTION,Front*,Wave*,
					 double,RP_DATA*,NODE_FLAG);
LOCAL	int	return_dnp(int,NODE*);
LOCAL	int	set_limits_for_amn_find_root(OVERSHOOT*,double*,double*,
					     double*,double*,NODE_FLAG);
LOCAL	void	check_for_incomplete_tir_deprecursor(int,NODE*,NODE*,
						     O_CURVE**,O_CURVE**,
						     double,double,Front*,
						     Wave*,RPROBLEM**);
LOCAL	void	temp_tnode_normal(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				  double*,Front*);
#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void	print_diffraction_debugging(NODE*,NODE*,O_CURVE**,RP_DATA*,
					    Front*,double**);
#endif /* defined(DEBUG_NODE_PROPAGATE) */


LOCAL NORMAL_FUNCTION sav_normal;

/*
*			  diffraction_node_propagate():
*
*
*
*	reflected wave
* (curve 2)  reflected shock or
* (curves 1 and 3) leading and 
* trailing edges of reflected 
* rarefaction
*			     \               / 
*                             \             /incident shock (curve 0)
*			       \ state 1   /
*				\	  /
*		state 4		 \	 /
*                                 \     /     state 0
*                                  \   /
*                                   \ /
*           -------------------------/-----------------------------
* (curve 4) contact(back)          /       contact(front) (curve 6)
*                                / 
*                               /       state 6
*                  state 5    /
*                           / transmitted shock (curve 5)
*
*
*	A diffraction node is a node at which an incident shock meets a
*	contact discontinuity producing a transmitted shock and a reflected
*	wave.  The forward facing (low pressure) side of the incident shock
*	defines the front side of the wave and hence an incident to front
*	angular orientation of the node.  The propagation direction (ahead)
*	of the node also provides an angular orientation, the incident to
*	ahead direction.
*
*	Front and back refer to the low and high pressure sides of the 
*	incident shock, while ahead and behind refer to the direction
*	of node propagation.
*
*	Note: when compiling for defined(DEBUG_NODE_PROPAGATE), the files
*	gprt/gprstate.c and gprt/gprcur.c must also be included so that the
*	functions print_diffraction_status() and verbose_print_bond_states()
*	respectively are defined.
*/


/* ARGSUSED */
EXPORT int diffraction_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag,
	POINTER		user)
{
	BOND		*newb[7];
	RECT_GRID	*gr = fr->rect_grid;
	RP_DATA		*RP;
	double		alpha;
	double		tcr[2];		 /* fractional distances to cross */
	double		va[MAXD], qa, qasqr;
	int		status, i;
	ANGLE_DIRECTION	i_to_f_dir; /* dir of angle to contact faced by */
				    /*  low pressure side of cinc  */
	boolean		is_plus_or, is_reflected_shock;
	boolean		c_ext[2];
	int		dim = fr->rect_grid->dim;
	static POINT	*pc = NULL;			/* cross point */
	static O_CURVE	**oldc = NULL, **newc = NULL;
	static double	**t = NULL;

#if defined(DEBUG_NODE_PROPAGATE)
static	char header[7][100] = {
		"INCIDENT SHOCK CURVE",
		"REFLECTED RAREFACTION LEADING EDGE",
		"REFLECTED SHOCK CURVE",
		"REFLECTED RAREFACTION TRAILING EDGE",
		"BACK CONTACT CURVE",
		"TRANSMITTED SHOCK CURVE",
		"FRONT CONTACT CURVE",
	};
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	debug_print("diffraction","Entered diffraction_node_propagate()\n");
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction")) 
	{
	    (void) printf("dt = %g\n",dt);
	    (void) printf("\n\tOLD NODE:\n");
	    print_node(oldn);
	    (void) printf("\n\tNEW NODE:\n");
	    print_node(newn);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */


		/* Allocate storage */

	if (oldc == NULL) 
	{
	    uni_array(&oldc,7,sizeof(O_CURVE *));
	    uni_array(&newc,7,sizeof(O_CURVE *));
	    bi_array(&t,2,MAXD,FLOAT);
	    for (i = 0; i < 7; i++)
	    {
	    	scalar(&oldc[i],sizeof(O_CURVE));
	    	scalar(&newc[i],sizeof(O_CURVE));
	    }
	    pc = Static_point(fr->interf);
	}
	for (i = 0; i < 7; i++)
	{
	    oldc[i]->curve = newc[i]->curve = NULL;
	}


		/* Identify curves */
	if (!identify_curves_at_shock_diffraction(oldn,newn,oldc,newc,
			                             &i_to_f_dir,&is_plus_or,
						     fr,(NODE*)user,flag))
	{
	    status = ERROR_NODE;
	    node_warning("diffraction_node_propagate",
	                 "unable to identify curves","\n",flag);
	    return return_dnp(status,newn);
	}


	/* Check incident shock and contact */
	for (i = 0; i < 7; i += 6)
	{
	    if (is_short_and_opp_unprop(oldc[i],newc[i],gr))
	    {
	    	set_prop_status_for_pseudo_cross_node(oldc[0],newc[0],
				                      oldc[6],newc[6],fr,
						      (POINTER)wave,dt,flag);
		status = PSEUDOCROSS_NODE_NODE;
		return return_dnp(status,newn);
	    }
	}

	/*
	*   If only the incident shock and contact are tracked allow
	*   subsonic states in the shock polar solution.
	*/
	for (i = 1; i < 6; i++)
	    if ((i != 4) && (newc[i]->curve != NULL))
		break;
	if (i == 6)
	    use_subsonic_state(flag) = YES;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction")) 
	{
	    if (oldn != NULL)
	    {
		double *nv = Node_vel(oldn);
	        double vstate[3], vsteady[3];
	        int   j;

	    	(void) printf("RP DATA at old node\n");
	    	print_RP_DATA(Rp_data(oldn),nv);
	    	(void) printf("Verbose states at node\n");
	    	for (i = 0; i < 7; i++)
	    	{
	    	    char s[120];

	    	    (void) sprintf(s,"Rp_data(oldn)->state[%d]",i);
	    	    verbose_print_state(s,Rp_data(oldn)->state[i]);
		    (void) VelocityVector(Rp_data(oldn)->state[i],vstate);
		    for (j = 0; j < dim; j++)
			vsteady[i] = vstate[i] - nv[i];
	    	    (void) sprintf(s,"Steady velocity of "
				     "Rp_data(oldn)->state[%d] = ",i);
		    print_general_vector(s,vsteady,dim,"\n");
	    	}
	    	(void) printf("End verbose states at node\n");
	    }
	    for (i = 0; i < 7; i++)
	    {
	    	if (oldc[i]->curve == NULL)
		    continue;
	    	(void) printf("\t\tOLD %s:\n",header[i]);
	    	if (debugging("states"))
	    	    verbose_print_curve_states(oldc[i]->curve);
	    	else
	    	    print_o_curve(oldc[i]);
	    }
	    for (i = 0; i < 7; i++)
	    {
	    	if (newc[i]->curve == NULL)
		    continue;
	    	(void) printf("\t\tNEW %s:\n",header[i]);
	    	if (debugging("states"))
	    	    verbose_print_curve_states(newc[i]->curve);
	    	else
	    	    print_o_curve(newc[i]);
	    }
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

		/* Identify new position of node */

	status = cross_or_extend_to_cross_two_propagated_curves(oldc[0],
			newc[0],oldc[6],newc[6],&pc,&newb[0],&newb[6],
			&tcr[0],&tcr[1],fr,(POINTER)wave,rp,dt,dt_frac,
			flag,c_ext);
	
	if (status != GOOD_NODE)
	{
	    /*
	     * If this is a deprecursion, make sure all relevant 
	     * nodes are included in Riemann problem
	     */

	    check_for_incomplete_tir_deprecursor(status,oldn,newn,oldc,newc,dt,
	    				         *dt_frac,fr,wave,rp);

	    return return_dnp(status,newn);
	}

	    /* Read off states from the propagated incident curve */


	RP = Rp_data(newn);
	RP->ang_dir = Opposite_ang_dir(i_to_f_dir);

	if (node_vel_by_angle(flag) == YES)
	{
	    if (debugging("diffraction"))
	    {
	    	(void) printf("Using angle to set node velocity\n");
	    	(void) printf("Node velocity before reset = <%g, %g>\n",
			      Node_vel(newn)[0],Node_vel(newn)[1]);
	    }
	    find_node_vel_at_rp(pc,tcr[0],newb[0],oldc[0],newc[0],t[0],
	    	                tcr[1],newb[6],oldc[6],newc[6],t[1],
	    		        RP->ang_dir,node_type(newn),
	    		        fr,wave,dt,Node_vel(newn));
	    if (debugging("diffraction"))
	    {
	    	(void) printf("Node velocity after reset = <%g, %g>\n",
	    		      Node_vel(newn)[0],Node_vel(newn)[1]);
	    }
	    status = velocity_satisfies_CFL(newn,dt,dt_frac,fr);
	    if (status != GOOD_NODE)
	    {
	    	return return_dnp(status,newn);
	    }
	}


#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
	    Locstate st;

	    print_angle_direction("i_to_f_dir =",i_to_f_dir,"\n");
	    (void) printf("is_plus_or = %s\n",(is_plus_or) ? "YES" : "NO");
	    (void) printf("Node velocity = <%g, %g>\n",
	    		  Node_vel(newn)[0],Node_vel(newn)[1]);
	    (void) printf("States on crossing bonds\n");
	    print_orientation("newc[0]->orient =",newc[0]->orient,"\n");
	    (void) printf("tcr[0] = %g\n",tcr[0]);
	    (void) printf("c_ext[0] = %s\n",(c_ext[0]) ? "YES" : "NO");
	    verbose_print_bond_states("newb[0]",newb[0],newc[0]->curve);
	    st = left_state_at_point_on_curve(newb[0]->start,newb[0],
	    	                              newc[0]->curve);
	    (void) printf("Mach number left state newb[0]->start = %g\n",
	    	          mach_number(st,Node_vel(newn)));
	    st = right_state_at_point_on_curve(newb[0]->start,newb[0],
	    	                               newc[0]->curve);
	    (void) printf("Mach number right state newb[0]->start = %g\n",
	    	          mach_number(st,Node_vel(newn)));
	    st = left_state_at_point_on_curve(newb[0]->end,newb[0],
	    	                              newc[0]->curve);
	    (void) printf("Mach number left state newb[0]->end = %g\n",
	    	          mach_number(st,Node_vel(newn)));
	    st = right_state_at_point_on_curve(newb[0]->end,newb[0],
	    	                               newc[0]->curve);
	    (void) printf("Mach number right state newb[0]->end = %g\n",
	    	          mach_number(st,Node_vel(newn)));

	    print_orientation("newc[6]->orient =",newc[6]->orient,"\n");
	    (void) printf("tcr[1] = %g\n",tcr[1]);
	    (void) printf("c_ext[1] = %s\n",(c_ext[1]) ? "YES" : "NO");
	    verbose_print_bond_states("newb[6]",newb[6],newc[6]->curve);
	    st = left_state_at_point_on_curve(newb[6]->start,newb[6],
	    	                              newc[6]->curve);
	    (void) printf("Mach number left state newb[6]->start = %g\n",
	    	          mach_number(st,Node_vel(newn)));
	    st = right_state_at_point_on_curve(newb[6]->start,newb[6],
	    	                               newc[6]->curve);
	    (void) printf("Mach number right state newb[6]->start = %g\n",
	    	          mach_number(st,Node_vel(newn)));
	    st = left_state_at_point_on_curve(newb[6]->end,newb[6],
	    	                              newc[6]->curve);
	    (void) printf("Mach number left state newb[6]->end = %g\n",
	    	          mach_number(st,Node_vel(newn)));
	    st = right_state_at_point_on_curve(newb[6]->end,newb[6],
	    	                               newc[6]->curve);
	    (void) printf("Mach number right state newb[6]->end = %g\n",
	    	          mach_number(st,Node_vel(newn)));
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	/* Compute new states by shock polars */

	assign_ahead_states_at_shock_diffraction(tcr[0],(c_ext[0])?NULL:newb[0],
		                                 newc[0]->curve,
                                                 newc[0]->orient,tcr[1],
		                                 (c_ext[1])?NULL:newb[6],
		                                 newc[6]->curve,
                                                 newc[6]->orient,
		                                 NULL,RP);

	/* Check for transition to sonic incident shock */

	if (mach_number(RP->state[0],Node_vel(newn)) < SONIC_PLUS)
	{
	    /* Don't allow bifurcations in buffer zone */
	    if ((oldn != NULL) &&
	        (point_in_buffer(Coords(oldn->posn),
	              computational_grid(oldn->interface)) == YES))
	        return set_node_states_and_continue(oldn,newn,fr);

	    status = sonic_incident_shock_at_diffraction(newn,oldc,newc,newb,
							 pc,tcr,RP,fr,wave,
							 rp,dt,dt_frac,flag);

	    if (status != GOOD_NODE) 
	    {
	        return return_dnp(status,newn);
	    }
	}
	else
	{
	    find_tangent_to_propagated_curve(pc,newb[0],oldc[0],newc[0],t[0],
					     fr,(POINTER)wave,dt);

	    for (i = 0; i < dim; i++)
	        va[i] = vel(i,RP->state[0]) - Node_vel(newn)[i];
	    qasqr = scalar_product(va,va,dim);
	    qa = sqrt(qasqr);
	    for (i = 0; i < dim; i++)
		t[1][i] = -va[i]/qa;

#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("diffraction"))
	    {
	        print_diffraction_debugging(oldn,newn,oldc,RP,fr,t);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    status = is_regular_diffraction_node(Coords(pc),Node_vel(newn),
	                                         (oldn != NULL) ?
						     Node_vel(oldn) : NULL,
	                                         t,RP,
	                                         (oldn != NULL) ?
						     Rp_data(oldn) : NULL,
	                                         &is_reflected_shock,fr,
	                                         node_type(newn),flag); 
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("diffraction"))
	    {
	        print_diffraction_status(
	            "is_regular_diffraction_node() returns ",status);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    /* Don't allow bifurcations in buffer zone */
	    if ((status != REGULAR_DIFFRACTION) && (oldn != NULL) &&
	        (point_in_buffer(Coords(oldn->posn),
	             computational_grid(oldn->interface)) == YES))
	        return set_node_states_and_continue(oldn,newn,fr);

	    switch (status)
	    {
	    case REGULAR_DIFFRACTION:
	         status = modify_diffraction_node(pc,oldn,newn,oldc,newc,
						  newb,fr,wave,dt,RP,flag);
	        if (status == GOOD_NODE)
	            propagation_status(newn) = PROPAGATED_NODE;
	        return return_dnp(status,newn);

	    case ANOMALOUS_REFLECTION:
	        if (fr->step == 0 && oldn != NULL)
	        {
	            for (i = 0; i < 7; i++)
	            {
	                ft_assign(RP->state[i],Rp_data(oldn)->state[i],
	                       fr->sizest);
	                RP->ang[i] = Rp_data(oldn)->ang[i];
	            }
	            is_reflected_shock = NO;
	            status = GOOD_NODE;
	        }
	        else
	        {
	            status = anomalous_reflection_propagate(fr,wave,oldn,newn,
							    rp,dt,dt_frac,flag,
							    oldc,newc,
	                                                    &is_reflected_shock,
	                                                    RP,pc,tcr,newb);
	        }
	        break;

	    case REGULAR_TO_MACH_DIFFRACTION:
	        status = reg_to_mach_diff_reconfigure(newc,oldc,oldn);
	        break;

	    case PRECURSOR_WITH_REFLECTED_RAREFACTION:
	        status = precursor_shock_rr_propagate(fr,wave,oldn,newn,
						      rp,dt,dt_frac,flag,
	                                              oldc,newc,pc,RP,newb);
		if (propagation_status(newn) == DELETED_NODE)
	            return return_dnp(status,newn);
	        break;

	    case PRECURSOR_WITH_REFLECTED_SHOCK:
	        status = ERROR_NODE;
	        screen("ERROR in diffraction_node_propagate(), "
	               "PRECURSOR_WITH_REFLECTED_SHOCK CODE NEEDED\n");
	        clean_up(ERROR);
	        break;

	    case ERROR_DIFFRACTION:
	    default:
	        status = ERROR_NODE;
	        break;
	    }
	    switch (status)
	    {
	    case GOOD_NODE:
	        break;

	    case PSEUDOCROSS_NODE_NODE:
	    case MODIFY_TIME_STEP_NODE:
	    case REPEAT_TIME_STEP_NODE:
	        return return_dnp(status,newn);

	    default:
	        if (interacting_nodes(oldn,newn,oldc,newc,dt,rp,fr,wave))
	            status = CROSS_NODE_NODE;
	        else
	        {
		    node_warning("diffraction_node_propagate",
	                          "is_regular_diffraction_node() failed\n"
	                          "unable to compute node configuration, "
	                          "possible bifurcation,  CODE NEEDED",
				  "\n",flag);
	            status = ERROR_NODE;
	        }
	        return return_dnp(status,newn);
	    }

	    if ((oldc[2]->curve != NULL && !is_reflected_shock) ||
	        (is_reflected_shock && (oldc[1]->curve || oldc[3]->curve)))
	    {
	        if (interacting_nodes(oldn,newn,oldc,newc,dt,rp,fr,wave))
	            status = CROSS_NODE_NODE;
	        else
	        {
		    char s[80];
		    (void) sprintf(s,"The reflected wave has undergone a "
				     "bifurcation.  CODE NEEDED\n"
	                             "is_reflected_shock = %s\n",
	                             y_or_n(is_reflected_shock));
		    node_warning("diffraction_node_propagate",s,"\n",flag);
	            status = ERROR_NODE;
	        }
	        return return_dnp(status,newn);
	    }
	}
	if ((Geometry(&alpha) != RECTANGULAR))
	{
	    static double nor[3] = {1.0, 0.0, 0.0};
	    double	 *pt = Coords(oldn->posn);
	    double        vn0;
	    Locstate	 sl, sr;
	    double	 *W = Node_vel(newn);
	    boolean	 do_ahead_states;

	    do_ahead_states = 
                (fr->_point_propagate == UnsplitPointPropagate(fr)) ? NO : YES;

	    if (do_ahead_states == YES)
	    {
	        sl = Left_state_at_node_of_o_curve(newc[0]);
		vn0 = scalar_product(VelocityVector(sl,NULL),t[0],dim);
	        include_source(pt,sl,vn0,dt,t[0],W,NULL,gr,0.0,
			       wave_type(newc[0]->curve));

	        sr = Right_state_at_node_of_o_curve(newc[0]);
		vn0 = scalar_product(VelocityVector(sr,NULL),t[0],dim);
	        include_source(pt,sr,vn0,dt,t[0],W,NULL,gr,0.0,
			       wave_type(newc[0]->curve));

	        sl = Left_state_at_node_of_o_curve(newc[6]);
		vn0 = scalar_product(VelocityVector(sl,NULL),t[6],dim);
	        include_source(pt,sl,vn0,dt,t[6],W,NULL,gr,0.0,
			       wave_type(newc[6]->curve));

	        sr = Right_state_at_node_of_o_curve(newc[6]);
		vn0 = scalar_product(VelocityVector(sr,NULL),t[6],dim);
	        include_source(pt,sr,vn0,dt,t[6],W,NULL,gr,0.0,
			       wave_type(newc[6]->curve));
	    }

	    if ((newc[1] != NULL) && (do_ahead_states == YES))
	    {
	        if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[1]->orient))
		{
		    sl = Left_state_at_node_of_o_curve(newc[1]);
		    vn0 = scalar_product(VelocityVector(sl,NULL),t[1],dim);
	            include_source(pt,sl,vn0,dt,t[1],W,NULL,gr,0.0,
			           wave_type(newc[1]->curve));

	            sr = Right_state_at_node_of_o_curve(newc[1]);
		    set_state(sr,state_type(sl),sl);
		}
		else
	        {
		    sr = Right_state_at_node_of_o_curve(newc[1]);
		    vn0 = scalar_product(VelocityVector(sr,NULL),t[1],dim);
	            include_source(pt,sr,vn0,dt,t[1],W,NULL,gr,0.0,
			           wave_type(newc[1]->curve));
		    sl = Left_state_at_node_of_o_curve(newc[1]);
		    set_state(sl,state_type(sr),sr);
	        }
	    }
	    if (newc[2] != NULL)
	    {
	        if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[2]->orient))
		{
	            if (do_ahead_states == YES)
	            {
		        sl = Left_state_at_node_of_o_curve(newc[2]);
		        vn0 = scalar_product(VelocityVector(sl,NULL),t[2],dim);
	                include_source(pt,sl,vn0,dt,t[2],W,NULL,gr,0.0,
				       wave_type(newc[2]->curve));
	            }
	            sr = Right_state_at_node_of_o_curve(newc[2]);
		    vn0 = scalar_product(VelocityVector(sr,NULL),nor,dim);
	            include_source(pt,sr,vn0,dt,nor,W,NULL,gr,0.0,
			           wave_type(newc[2]->curve));
		}
		else
	        {
	            if (do_ahead_states == YES)
	            {
		        sr = Right_state_at_node_of_o_curve(newc[2]);
		        vn0 = scalar_product(VelocityVector(sr,NULL),t[2],dim);
	                include_source(pt,sr,vn0,dt,t[2],W,NULL,gr,0.0,
				       wave_type(newc[2]->curve));
	            }
		    sl = Left_state_at_node_of_o_curve(newc[2]);
		    vn0 = scalar_product(VelocityVector(sl,NULL),nor,dim);
	            include_source(pt,sl,vn0,dt,nor,W,NULL,gr,0.0,
			           wave_type(newc[2]->curve));
	        }
	    }
	    if (newc[5] != NULL)
	    {
	        if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[5]->orient))
		{
	            if (do_ahead_states == YES)
	            {
		        sr = Right_state_at_node_of_o_curve(newc[5]);
		        vn0 = scalar_product(VelocityVector(sr,NULL),t[5],dim);
	                include_source(pt,sr,vn0,dt,t[5],W,NULL,gr,0.0,
				       wave_type(newc[5]->curve));
	            }
		    sl = Left_state_at_node_of_o_curve(newc[5]);
		    vn0 = scalar_product(VelocityVector(sl,NULL),nor,dim);
	            include_source(pt,sl,vn0,dt,nor,W,NULL,gr,0.0,
			           wave_type(newc[5]->curve));
		}
		else
	        {
	            if (do_ahead_states == YES)
	            {
		        sl = Left_state_at_node_of_o_curve(newc[5]);
		        vn0 = scalar_product(VelocityVector(sl,NULL),t[5],dim);
	                include_source(pt,sl,vn0,dt,t[5],W,NULL,gr,0.0,
				       wave_type(newc[5]->curve));
	            }
	            sr = Right_state_at_node_of_o_curve(newc[5]);
		    vn0 = scalar_product(VelocityVector(sr,NULL),nor,dim);
	            include_source(pt,sr,vn0,dt,nor,W,NULL,gr,0.0,
			           wave_type(newc[5]->curve));
	        }
	    }
	    for (i = 3; i < 5; i++)
	    {
	        if (newc[i] != NULL)
		{
		    sl = Left_state_at_node_of_o_curve(newc[i]);
		    vn0 = scalar_product(VelocityVector(sl,NULL),nor,dim);
	            include_source(pt,sl,vn0,dt,nor,W,NULL,gr,0.0,
			           wave_type(newc[i]->curve));
	            sr = Right_state_at_node_of_o_curve(newc[i]);
		    vn0 = scalar_product(VelocityVector(sr,NULL),nor,dim);
	            include_source(pt,sr,vn0,dt,nor,W,NULL,gr,0.0,
			           wave_type(newc[i]->curve));
		}
	    }
	}
	return return_dnp(status,newn);
}		/*end diffraction_node_propagate*/


/*
*			modify_diffraction_node():
*
*	This function simply drives the calls to modify_curves_at_node()
*	for a diffraction node.  A call to this function is needed to
*	complete the propagation of a diffraction node.  In the regular
*	diffraction case, the call is done from diffraction_node_propagate().
*	For all bifurcation codes, it is assumed that they call this function
*	themselves.  This behaviour was needed because certain of these
*	bifurcations required the processing in this function before they
*	could complete installation of the new configuration.
*/

EXPORT	int	modify_diffraction_node(
	POINT		*pc,
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	Front		*fr,
	Wave		*wave,
	double		dt,
	RP_DATA		*RP,
	NODE_FLAG	flag)
{
	int		i;
	static boolean	correct_angle[7] = {NO,YES,YES,YES,YES,YES,NO};
	double		min_adjust_len, adjust_len[7];
#if defined(DEBUG_NODE_PROPAGATE)
	static char	header[7][100] = {
				"INCIDENT SHOCK CURVE",
				"REFLECTED RAREFACTION LEADING EDGE",
				"REFLECTED SHOCK CURVE",
				"REFLECTED RAREFACTION TRAILING EDGE",
				"BACK CONTACT CURVE",
				"TRANSMITTED SHOCK CURVE",
				"FRONT CONTACT CURVE",
			};
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	for (i = 1; i < 6; i++)
	{
	    newb[i] = (newc[i] && newc[i]->curve) ? 
				Bond_at_node_of_o_curve(newc[i]) : NULL;
	}
	adjust_len[0] = adjust_len[6] = 0.0;
	min_adjust_len = HUGE_VAL;
	for (i = 1; i < 6; i++)
	{
	    adjust_len[i] = sonic_radius(RP->state[i],
					 Node_vel(newn),dt,fr->rect_grid->dim);
	    min_adjust_len = min(min_adjust_len,adjust_len[i]);
	}
	for (i = 1; i < 6; i++)
	    adjust_len[i] = min_adjust_len;
	adjust_len(newn) = min_adjust_len;
	if (!modify_curves_at_node(pc,newb,newn,oldc,newc,7,
				      correct_angle,fr,wave,dt,RP,flag))
	{
	    node_warning("modify_diffraction_node",
	                 "modify_curves_at_node() failed","\n",flag);
	    if (oldn != NULL && next_node(oldn)) 
	    {
	    	set_prop_status_for_pseudo_cross_node(oldc[0],newc[0],
				                      oldc[6],newc[6],
						      fr,(POINTER)wave,
						      dt,flag);
	    	return PSEUDOCROSS_NODE_NODE;
	    }
	    else 
	    {
	    	return ERROR_NODE;
	    }
	}
	adjust_len(newn) = 0.0;
	for (i = 0; i < 7; i++)
	    if ((newc[i]->curve != NULL) && (correct_angle[i] == YES))
		redistribute_hyper_surface(newc[i]->curve) = YES;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction")) 
	{
	    Locstate st;
	    double vl,vn,alpha,c,s;

	    (void) printf("NEW DIFFRACTION NODE\n");	print_node(newn);
	    for (i = 0; i < 7; i++)
	    {
	    	if (newc[i]->curve == NULL) continue;
	    	(void) printf("\t\tNEW %s:\n",header[i]);
	    	if (debugging("states"))
	    		verbose_print_curve_states(newc[i]->curve);
	    	else
	    		print_o_curve(newc[i]);
	    }
	    alpha = RP->ang[4];
	    c = cos(alpha);
	    s = sin(alpha);
	    print_angle("Behind contact direction =",alpha,"\n");
	    st = Left_state_at_node_of_o_curve(newc[4]);
	    vl = c*vel(0,st) + s*vel(1,st);
	    vn = -s*vel(0,st) + c*vel(1,st);
	    (void) printf("Left tan. and nor. velocities: %g %g\n",vl,vn);
	    st = Right_state_at_node_of_o_curve(newc[4]);
	    vl = c*vel(0,st) + s*vel(1,st);
	    vn = -s*vel(0,st) + c*vel(1,st);
	    (void) printf("Right tan. and nor. velocities: %g %g\n",vl,vn);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return GOOD_NODE;
}		/*end modify_diffraction_node*/

LOCAL	int return_dnp(
	int	status,
	NODE	*newn)
{
	if (debugging("dintfc"))
	{
	    (void) printf("Final interface from "
			  "diffraction_node_propagate()\n");
	    print_interface(newn->interface);
	}
	debug_print("diffraction","Left diffraction_node_propagate(), ");
	if (debugging("diffraction"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end return_dnp*/


/*
*		check_for_incomplete_tir_deprecursor();
*
*	This function checks to make sure all relevant nodes are included
*	in the rp for a deprecursion from precursor with reflected 
*	rarefaction back to a single diffraction node.  We allow for the
*	case of untracked overtake node.  The tir in the name means
*	the deprecursion is being triggered by a total internal reflection
*	node.
*/

LOCAL	void	check_for_incomplete_tir_deprecursor(
	int		status,
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	double		dt,
	double		dt_frac,
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp)
{
	NODE		*interact_nodes[9];
	RECT_GRID	*gr;
	RP_NODE		*rpn;
	boolean		tir_node, c_node, t_node, o_node;

	if (status != CROSS_NODE_NODE)
	    return; 
	if (rp == NULL || *rp == NULL)
	    return;
	tir_node = c_node = t_node = o_node = NO;
	for (rpn = (*rp)->first_rp_node; rpn != NULL; rpn = rpn->next)
	{
	    if (node_type(rpn->node) == TOT_INT_REFL_NODE)
	    	tir_node = YES;
	    if (node_type(rpn->node) == CROSS_NODE)
	    	c_node = YES;
	    if (node_type(rpn->node) == TRANSMISSION_NODE)
	    	t_node = YES;
	    if (node_type(rpn->node) == OVERTAKE_NODE)
	    	o_node = YES;
	}
	if ((tir_node == YES) && (c_node == YES) &&
	    (t_node == YES) && (o_node == YES))
	    return;

	gr = fr->rect_grid;
	if (tir_node == NO)
	    return;

	c_node = (
		  (c_node == YES) ||
		  (
		    (node_type(Opp_node_of_o_curve(oldc[0])) == CROSS_NODE)
		 		 &&
		    is_short_curve(oldc[0]->curve,oldc[0]->orient,gr,1.0)
		  )
		 ) ? YES : NO;

	t_node = (
		  (t_node == YES) ||
		  (
		    (node_type(Opp_node_of_o_curve(oldc[6]))==TRANSMISSION_NODE)
		 		 &&
		    is_short_curve(oldc[6]->curve,oldc[6]->orient,gr,1.0)
		  )
		 ) ? YES : NO;

	o_node = (
		  (o_node == YES) ||
		  (
		    oldc[1]->curve	    /* allow for untracked rarefaction*/
		 	 &&
		    (node_type(Opp_node_of_o_curve(oldc[1])) == OVERTAKE_NODE)
		 	 &&
		    is_short_curve(oldc[1]->curve,oldc[1]->orient,gr,1.0)
		  )
		 ) ? YES : NO;

	if ((c_node == YES) && (t_node == YES))
	{
	    interact_nodes[0] = newn;
	    interact_nodes[1] = oldn;
	    interact_nodes[2] = Opp_node_of_o_curve(newc[0]);
	    interact_nodes[3] = Opp_node_of_o_curve(oldc[0]);
	    interact_nodes[4] = Opp_node_of_o_curve(newc[6]);
	    interact_nodes[5] = Opp_node_of_o_curve(oldc[6]);
	    if (o_node)
	    {
	    	interact_nodes[6] = Opp_node_of_o_curve(newc[1]);
	    	interact_nodes[7] = Opp_node_of_o_curve(oldc[1]);
	    	interact_nodes[8] = NULL;
	    }
	    else
	    	interact_nodes[6] = NULL;

	    augment_rproblem_list(rp,interact_nodes,dt,dt_frac,oldn->interface,
				  newn->interface,fr,(POINTER)wave);
	}
}		/*end check_for_incomplete_tir_deprecursor*/


/*
*		identify_curves_at_shock_diffraction():
*
*	First identifies the curves at the old node by calling 
*	curves_at_shock_diffraction().  The corresponding curves at the
*	new node are then found via repeated calls to 
*	find_correspond_of_oriented_curve().  It is an error if the
*	incident or contact target do not exist at the new node.  All other
*	curves may be null, irregardless of whether they existed at the
*	previous time step.  This allows for dynamic changes in topology
*	during the node propagate loop.
*/

LOCAL	int identify_curves_at_shock_diffraction(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	ANGLE_DIRECTION	*i_to_f_dir,
	boolean		*is_plus_or,
	Front		*fr,
	NODE		*user,
	NODE_FLAG       flag)
{
	int	i;

	if (oldn != NULL)
	{
	    if (curves_at_shock_diffraction(oldn,oldc,i_to_f_dir,
					    node_warnings_off(flag)) != YES)
	    {
		node_warning("identify_curves_at_shock_diffraction",
	    	             "Unable to find curves","\n",flag);
	    	return NO;
	    }
	    if (!find_correspond_of_oriented_curve(oldc[0],newc[0],newn,fr,
						      newn->interface))
	    {
	    	find_curve_with_status(newn,&newc[0]->curve,&newc[0]->orient,
				       INCIDENT);
		Check_return(
		    find_correspond_of_oriented_curve(newc[0],oldc[0],
						      user,fr,oldn->interface),
		    identify_curves_at_shock_diffraction)
	    }
	    for (i = 1; i < 6; i++)
	    {
	    	(void) find_correspond_of_oriented_curve(oldc[i],newc[i],newn,
							 fr,newn->interface);
	    }
	    Check_return(
		find_correspond_of_oriented_curve(oldc[6],newc[6],newn,fr,
						  newn->interface),
		identify_curves_at_shock_diffraction)
	}
	else
	{
	    if (curves_at_shock_diffraction(newn,newc,i_to_f_dir,
					    node_warnings_off(flag)) != YES)
	    {
		node_warning("identify_curves_at_shock_diffraction",
	    	              "Unable to find curves","\n",flag);
	    	return NO;
	    }
	    for (i = 0; i < 7; i++)
	    {
	    	Check_return(
		    find_correspond_of_oriented_curve(newc[i],oldc[i],oldn,fr,
						      fr->interf),
		    identify_curves_at_shock_diffraction)
	    }
	}
	*is_plus_or = (*i_to_f_dir == CLOCKWISE) ? YES : NO;
	return YES;
}		/*end identify_curves_at_shock_diffraction*/


/*
*			curves_at_shock_diffraction():
*
*	Finds the curves at a diffraction node and loads them in an
*	O_CURVE array.  If the incident or target contact are not found,
*	this is considered an error, and NO is returned.  Any or all other
*	curves may not be found, with no error signal.  It is assumed
*	these waves are untracked.
*/	

EXPORT	boolean curves_at_shock_diffraction(
	NODE		*node,
	O_CURVE		**oc,
	ANGLE_DIRECTION	*i_to_f_dir,
	boolean         warnings_off)
{
	O_CURVE		Ctmp0, Ctmp1;

	find_curve_with_status(node,&oc[0]->curve,&oc[0]->orient,INCIDENT);
	if (oc[0]->curve == NULL)
	{
	    if (warnings_off != YES)
	    {
	        (void) printf("WARNING in curves_at_shock_diffraction(), "
	                      "Unable to find incident curve\n");
	    }
	    return NO;
	}

	*i_to_f_dir =
	    incident_shock_orientation(wave_type(oc[0]->curve),oc[0]->orient);

	find_curve_with_status(node,&oc[6]->curve,&oc[6]->orient,
		               CONTACT_TARGET);
	if (oc[0]->curve == NULL)
	{
	    if (warnings_off != YES)
	    {
	        (void) printf("WARNING in curves_at_shock_diffraction(), "
	                      "Unable to find contact target\n");
	    }
	    return NO;
	}

	find_curve_with_status(node,&oc[4]->curve,&oc[4]->orient,SLIP);

	find_curve_with_status(node,&oc[5]->curve,&oc[5]->orient,TRANSMITTED);

	identify_curves_with_status(node,&Ctmp0,&Ctmp1,REFLECTED);
	if (Ctmp0.curve == NULL)
	{
	    /* Untracked reflected wave */
	    oc[1]->curve = oc[2]->curve = oc[3]->curve = NULL;
	}
	else if (is_shock_wave(wave_type(Ctmp0.curve)))
	{
	    /* Tracked reflected shock */
	    *oc[2] = Ctmp0;
	    oc[1]->curve = oc[3]->curve = NULL;
	}
	else 
	{
	    /* Tracked reflected rarefaction */
	    oc[2]->curve = NULL;
	    if (is_rarefaction_leading_edge(wave_type(Ctmp0.curve)))
	    {
	    	*oc[1] = Ctmp0;
	    	*oc[3] = Ctmp1;
	    }
	    else
	    {
	    	*oc[1] = Ctmp1;
	    	*oc[3] = Ctmp0;
	    }
	}
	return YES;
}		/*end curves_at_shock_diffraction*/

LOCAL	int interacting_nodes(
	NODE	  *oldn,
	NODE	  *newn,
	O_CURVE	  **oldc,
	O_CURVE	  **newc,
	double	  dt,
	RPROBLEM  **rp,
	Front	  *fr,
	Wave	  *wave)
{
	NODE		   *interact_nodes[7];
	RECT_GRID	   *gr = fr->rect_grid;
	static const double SC_TOL = 1.5; /*TOLERANCE*/

	if (oldn == NULL)
	    return NO;
	interact_nodes[0] = NULL;
	if (is_short_curve(oldc[0]->curve,oldc[0]->orient,gr,SC_TOL))
	{
	    interact_nodes[0] = newn;
	    interact_nodes[1] = oldn;
	    interact_nodes[2] = Opp_node_of_o_curve(newc[0]);
	    interact_nodes[3] = Opp_node_of_o_curve(oldc[0]);
	    interact_nodes[4] = NULL;
	}
	if (is_short_curve(oldc[6]->curve,oldc[6]->orient,gr,SC_TOL))
	{
	    if (interact_nodes[0] == NULL)
	    {
	    	interact_nodes[0] = newn;
	    	interact_nodes[1] = oldn;
	    	interact_nodes[2] = Opp_node_of_o_curve(newc[6]);
	    	interact_nodes[3] = Opp_node_of_o_curve(oldc[6]);
	    	interact_nodes[4] = NULL;
	    }
	    else
	    {
	    	interact_nodes[4] = Opp_node_of_o_curve(newc[6]);
	    	interact_nodes[5] = Opp_node_of_o_curve(oldc[6]);
	    	interact_nodes[6] = NULL;
	    }
	}
	if (interact_nodes[0] != NULL)
	{
	    propagation_status(newn) = VEL_COMPUTED_NODE;
	    augment_rproblem_list(rp,interact_nodes,
			          dt,1.0,oldc[0]->curve->interface,
			          newc[0]->curve->interface,fr,(POINTER)wave);
	    return YES;
	}
	return NO;
}		/*end interacting_nodes*/

EXPORT	int sonic_incident_shock_at_diffraction(
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		*pc,
	double		*tcr,
	RP_DATA		*RP,
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	double		theta;
	int		i, status;
	int		is_plus_or = (RP->ang_dir == COUNTER_CLOCK) ? YES : NO;

	debug_print("diffraction","Entered sonic_incident_shock_at_diffraction()\n");

		/* Subsonic state ahead of incident shock.	 */
		/* Probably caused by overshoot in extrapolation */
		/* of new node position.			 */

	if (!correct_new_node_posn_for_overshoot(newn,&tcr[0],&tcr[1],
			&newb[0],&newb[6],pc,oldc[0],newc[0],oldc[6],newc[6],RP,
			fr,wave,dt,flag))
	{
	    node_warning("sonic_incident_shock_at_diffraction",
	                 "correct_new_node_posn_for_overshoot() failed\n"
	                 "possible bifurcation,  CODE NEEDED\n","\n",flag);
	    status = ERROR_NODE;
	    debug_print("diffraction","Left sonic_incident_shock_at_diffraction(), ");
	    if (debugging("diffraction"))
		print_node_status("status = ",status,"\n");
	    return status;
	}

		/* disconnect reflected waves from node */

	for (i = 1; i < 4; i++)
	{
	    status = refl_curve_overtakes_incident_shock(oldc,newc,newb,pc,fr,
							 wave,rp,dt,dt_frac,i,
							 RP->ang_dir,flag);


	    if (status != GOOD_NODE) 
	    {
#if defined(DEBUG_NODE_PROPAGATE)
		if (debugging("diffraction"))
		{
		    (void) printf("Unable to propagate reflected "
		                  "curve up incident shock");
		}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
		debug_print("diffraction",
		      "Left sonic_incident_shock_at_diffraction(), ");
		if (debugging("diffraction")) 
		    print_node_status("status = ",status,"\n");
		return status;
			
	    }
	}

	if (!s_polar_3(RP->state[0],YES,pressure(RP->state[0]),is_plus_or,
			  NO,Node_vel(newn),RP->state[1],&RP->ang[0],&theta))
	{
	    node_warning("sonic_incident_shock_at_diffraction",
	                 "s_polar_3() failed","\n",flag);
	    status = ERROR_NODE;
	    debug_print("diffraction","Left sonic_incident_shock_at_diffraction(), ");
	    if (debugging("diffraction"))
		print_node_status("status = ",status,"\n");
	    return status;
	}

	for (i = 2; i < 5; i++)
	{
	    ft_assign(RP->state[i],RP->state[1],fr->sizest);
	    RP->ang[i-1] = RP->ang[0];
	}

	if (!s_polar_3(RP->state[6],YES,pressure(RP->state[6]),is_plus_or,
			  YES,Node_vel(newn),RP->state[5],&RP->ang[5],&theta))
	{
	    node_warning("sonic_incident_shock_at_diffraction",
			 "s_polar_3() failed","\n",flag);
	    status = ERROR_NODE;
	    debug_print("diffraction","Left sonic_incident_shock_at_diffraction(), ");
	    if (debugging("diffraction"))
		print_node_status("status = ",status,"\n");
	    return status;
	}

	RP->ang[6] = avg_angle_and_normalize(
		angle(Node_vel(newn)[0] - vel(0,RP->state[0]),
			Node_vel(newn)[1] - vel(1,RP->state[0])),
		angle(Node_vel(newn)[0] - vel(0,RP->state[6]),
			Node_vel(newn)[1] - vel(1,RP->state[6])));
	
	RP->ang[4] = normalized_angle(RP->ang[6] + PI);

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("diffraction"))
	{
		print_RP_node_states(
			"States after sonic_incident_shock_at_diffraction()",
			Node_vel(newn),RP,node_type(newn));
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	debug_print("diffraction",
			"Left sonic_incident_shock_at_diffraction(), ");
	if (debugging("diffraction"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end sonic_incident_shock_at_diffraction*/

#if defined(DEBUG_NODE_PROPAGATE)
LOCAL	void print_diffraction_debugging(
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		**oldc,
	RP_DATA		*RP,
	Front		*fr,
	double		**t)
{
	Locstate	sa, sb;
	double		sin_ang, cos_ang, ang;
	double		va[MAXD], qasqr, qa, mf_sqr;
	double		*abs_v = Node_vel(newn);
	double		ot[2];
	int		dim = fr->rect_grid->dim;
	int		i;

	(void) printf("Tangent uni_arrays\n");

	find_tangent_to_curve(Node_of_o_curve(oldc[0])->posn,
			Bond_at_node_of_o_curve(oldc[0]),oldc[0]->curve,
			oldc[0]->orient,ot,fr);
	(void) printf("Old incident shock tangent = <%g, %g>,\n",ot[0],ot[1]);
	ang = atan2(ot[1],ot[0]);
	(void) printf("\t\t\tmagnitude = %g, ",mag_vector(ot,dim));
	print_angle("angle =",ang,"\n");

	(void) printf("New incident shock tangent = <%g, %g>,\n",
		      t[0][0],t[0][1]);
	ang = atan2(t[0][1],t[0][0]);
	(void) printf("\t\t\tmagnitude = %g, ",mag_vector(t[0],dim));
	print_angle("angle =",ang,"\n");

	if (oldn != NULL)
	{
		sa = (is_forward_wave(wave_type(oldc[0]->curve))) ?
			Right_state_at_node_of_o_curve(oldc[0]) :
			Left_state_at_node_of_o_curve(oldc[0]);
		for (i = 0; i < dim; i++)
			va[i] = vel(i,sa) - Node_vel(oldn)[i];
		qasqr = scalar_product(va,va,dim);	qa = sqrt(qasqr);
		for (i = 0; i < dim; i++)
			ot[i] = -va[i]/qa;
		(void) printf("Old incident contact tangent = <%g, %g>,\n",
		       ot[0],ot[1]);
		ang = atan2(ot[1],ot[0]);
		(void) printf("\t\t\tmagnitude = %g, ",mag_vector(ot,dim));
		print_angle("angle =",ang,"\n");

		(void) printf("New incident contact tangent = <%g, %g>,\n",
			t[1][0],t[1][1]);
		ang = atan2(t[1][1],t[1][0]);
		(void) printf("\t\t\tmagnitude = %g, ",mag_vector(t[1],dim));
		print_angle("angle =",ang,"\n");
	}
		
	(void) vector_product(t[0],t[1],&sin_ang,dim);
	cos_ang = scalar_product(t[0],t[1],dim);
	(void) printf("Incident angle (shock to contact)\n");
	ang = atan2(sin_ang,cos_ang);
	(void) printf("sin(ang) = %g, cos(ang) = %g, ",sin_ang,cos_ang);
	print_angle("ang =",ang,"\n");

	(void) printf("Node velocity = <%g, %g>,\n",abs_v[0],abs_v[1]);
	ang = atan2(abs_v[1],abs_v[0]);
	(void) printf("\t\t\tmagnitude = %g, ",mag_vector(abs_v,dim));
	print_angle("angle =",ang,"\n");

	sa = RP->state[0], sb = RP->state[1];
	for (i = 0; i < dim; i++)
		va[i] = vel(i,sa) - abs_v[i]; 
	qasqr = scalar_product(va,va,dim);
	mf_sqr = sqr(Dens(sa)) * qasqr * sqr(sin_ang);
	(void) printf("Mass flux squared computed by angle = %g\n",mf_sqr);
	(void) printf("Mass flux squared computed by pressure = %g\n",
		mass_flux_squared(pressure(sb),sa));
}		/*end print_diffraction_debugging*/
#endif /* defined(DEBUG_NODE_PROPAGATE) */


LOCAL	int correct_new_node_posn_for_overshoot(
	NODE		*newn,
	double		*tcr0,
	double		*tcr6,
	BOND		**newb0,
	BOND		**newb6,
	POINT		*pc,
	O_CURVE		*oldc0,
	O_CURVE		*newc0,
	O_CURVE		*oldc6,
	O_CURVE		*newc6,
	RP_DATA		*RP,
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE_FLAG	flag)
{
	HYPER_SURF	*hs;
	HYPER_SURF_ELEMENT *hse;
	INTERFACE	*intfc = newn->interface;
	OVERSHOOT	Ovrsht;
	CURVE		*c0, *c6;
	POINT		*p0_sav, *p6_sav;
	NODE		*oppn0, *oppn6;
	BOND		*oppb0, *oppb6;
	BOND		B0virtual, B6virtual;
	BOND		*btmp, *bc6;
	double		coords[MAXD], V[MAXD];
	double		dir[MAXD], ds;
	double		epsilon, delta;
	double		ds_min, ds_max;
	ORIENTATION	c0_orient, c6_orient;
	int		status = NO;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	int		dim = intfc->dim;
	int		i;
	static	POINT	*p0 = NULL, *p0_opp = NULL,
	                *p6 = NULL, *p6_opp = NULL, *pbase = NULL;
	static	boolean	first = YES;

	debug_print("diffraction","Entered correct_new_node_posn_for_overshoot()\n");

	if (first)
	{
		first = NO;
		p6 = Static_point(fr->interf);
		p6_opp = Static_point(fr->interf);
		p0 = Static_point(fr->interf);
		p0_opp = Static_point(fr->interf);
		pbase = Static_point(fr->interf);
	}

	/* Check that ahead state has gone sonic */

	if (mach_number(RP->state[0],Node_vel(newn)) > SONIC_MINUS)
	{
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("diffraction"))
	        (void) printf("Ahead state is sonic\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    debug_print("diffraction","Left correct_new_node_posn_for_overshoot()\n");
	    return YES;
	}

	Ovrsht.fr = fr;			Ovrsht.dt = dt;
	for (i = 0; i < dim; i++) Ovrsht.abs_v0[i] = Node_vel(newn)[i];
	Ovrsht.c = c6 = newc6->curve;	Ovrsht.b = bc6 = *newb6;
	Ovrsht.p0 = pbase;			Ovrsht.pc = pc;
	Ovrsht.c_orient = c6_orient = newc6->orient;
	Ovrsht.inc_side = (curve_ang_oriented_l_to_r(RP->ang_dir,c6_orient)) ?
				POSITIVE_SIDE : NEGATIVE_SIDE;
	Ovrsht.sa = RP->state[0];	Ovrsht.sa_opp = RP->state[6];

	p6_sav = Point_of_bond(Bond_at_node(c6,c6_orient),c6_orient);
	init_curve_for_crossing(p6,p6_opp,&B6virtual,oldc6,newc6,
		&oppn6,&oppb6,fr,(POINTER)wave,dt,V,flag);
	set_point_of_bond(p6,Bond_at_node(c6,c6_orient),c6_orient,dim);
	find_tangent_to_curve(p6,Bond_at_node(c6,c6_orient),c6,c6_orient,
		dir,fr);
	ds = grid_size_in_direction(dir,fr->rect_grid->h,dim);
	for (i = 0; i < dim; i++) Coords(p6)[i] -= 2.0*ds*dir[i];
	interpolate_intfc_states(intfc) = YES;
	for (i = 0; i < dim; i++) Coords(pbase)[i] = Coords(pc)[i];
	if(insert_point_in_bond(pbase,bc6,c6) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in correct_new_node_posn_for_overshoot(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	btmp = bc6->next;

	if (!set_limits_for_amn_find_root(&Ovrsht,&ds_min,&ds_max,
		                             &epsilon,&delta,flag))
	{
	    status = NO;
	    goto leave;
	}

	if (find_root(amn,(POINTER)&Ovrsht,1.0,&ds,ds_min,
		      ds_max,epsilon,delta) == FUNCTION_FAILED)
	{
	    status = NO;
	    goto leave;
	}
	
	status = YES;
	*newb6 = (Ovrsht.bc == btmp) ? *newb6 : Ovrsht.bc;
	for (i = 0; i < dim; i++) Node_vel(newn)[i] = Ovrsht.abs_v[i];
	*tcr6 = Ovrsht.tcr;

		/* Find nearest point and bond to new pc on c0 */

	c0 = newc0->curve;		c0_orient = newc0->orient;
	p0_sav = Point_of_bond(Bond_at_node(c0,c0_orient),c0_orient);
	init_curve_for_crossing(p0,p0_opp,&B0virtual,oldc0,newc0,
		&oppn0,&oppb0,fr,(POINTER)wave,dt,V,flag);
	set_point_of_bond(p0,Bond_at_node(c0,c0_orient),c0_orient,dim);

	if (long_nearest_interface_point(Coords(pc),negative_component(c0),
		intfc,NO_BOUNDARIES,Hyper_surf(c0),coords,tcr0,&hse,&hs) == YES)
	{
	    *newb0 = Bond_of_hse(hse);
	    if (is_forward_wave(wave_type(c0)))
	    {
	    	left_state_along_bond(*tcr0,*newb0,c0,RP->state[1]);
	    }
	    else
	    {
	    	right_state_along_bond(*tcr0,*newb0,c0,RP->state[1]);
	    }
	}
	else
	    status = NO;

	set_point_of_bond(p0_sav,Bond_at_node(c0,c0_orient),c0_orient,dim);
	set_point_of_bond(oppn0->posn,oppb0,Opposite_orient(c0_orient),dim);


leave:
	(void) delete_start_of_bond(btmp,c6);
	set_point_of_bond(p6_sav,Bond_at_node(c6,c6_orient),c6_orient,dim);
	set_point_of_bond(oppn6->posn,oppb6,Opposite_orient(c6_orient),dim);
	interpolate_intfc_states(intfc) = sav_intrp;
	debug_print("diffraction","Left correct_new_node_posn_for_overshoot()\n");
	return status;
}		/*end correct_new_node_posn_for_overshoot*/


LOCAL	int set_limits_for_amn_find_root(
	OVERSHOOT  *ovrsht,
	double	  *pds_min,
	double	  *pds_max,
	double	  *epsilon,
	double	  *delta,
	NODE_FLAG flag)
{
	double		M0, Mds_mid, Mds_min, Mds_max;
	double		ds_mid, ds_min, ds_max, dslen;
	int		Mds_max_set, Mds_min_set;
	int		i;
	static const int   MAX_ITER = 20;
	static const double INITIAL_DSDT = 0.01;	/*TOLERANCE*/

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("diffraction","Entered set_limits_for_amn_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (amn(0.0,&M0,(POINTER)ovrsht) == FUNCTION_FAILED)
	{
	    node_warning("set_limits_for_amn_find_root",
	                 "Unable to evaluate amn()","\n",flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("diffraction","Left set_limits_for_amn_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return NO;
	}


	dslen = INITIAL_DSDT * ovrsht->dt;
	ds_min = -dslen;	Mds_min_set = NO;
	ds_max = dslen;		Mds_max_set = NO;
	ds_mid = 0.0;		Mds_mid = M0;
	for (i = 0; i < MAX_ITER; i++)
	{
	    if (!Mds_max_set)
	    {
	    	if (amn(ds_max,&Mds_max,(POINTER)ovrsht) == FUNCTION_FAILED)
	    	{
	            node_warning("set_limits_for_amn_find_root",
	                         "Unable to evaluate amn()","\n",flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    	    debug_print("diffraction",
	    		  "Left set_limits_for_amn_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	    return NO;
	    	}
	    }
	    Mds_max_set = YES;

	    if (!Mds_min_set)
	    {
	    	if (amn(ds_min,&Mds_min,(POINTER)ovrsht) == FUNCTION_FAILED)
	    	{
	            node_warning("set_limits_for_amn_find_root",
	                         "Unable to evaluate amn()","\n",flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    	    debug_print("diffraction",
	    		  "Left set_limits_for_amn_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    	    return NO;
	    	}
	    }
	    Mds_min_set = YES;

	    if (Mds_mid < 1.0)
	    {
	    	if (Mds_min > 1.0)
	    	{
	    	    ds_max = ds_mid;
	    	    break;
	    	}
	    	if (Mds_max > 1.0)
	    	{
	    	    ds_min = ds_mid;
	    	    break;
	    	}
	    	if (Mds_mid < Mds_min && Mds_mid > Mds_max)
	    	{
	    	    ds_max = ds_mid;	Mds_max = Mds_mid;
	    	    ds_mid = ds_min;	Mds_mid = Mds_min;
	    	    ds_min = 2.0*ds_mid - ds_max;
	    	    Mds_min_set = NO;
	    	}
	    	else if (Mds_mid > Mds_min && Mds_mid < Mds_max)
	    	{
	    	    ds_min = ds_mid;	Mds_min = Mds_mid;
	    	    ds_mid = ds_max;	Mds_mid = Mds_max;
	    	    ds_max = 2.0*ds_mid - ds_min;
	    	    Mds_max_set = NO;
	    	}
	    	else
	    	{
	    	    dslen *= 0.5;
	    	    ds_min = -dslen;	Mds_min_set = NO;
	    	    ds_max = dslen;		Mds_max_set = NO;
	    	    ds_mid = 0.0;		Mds_mid = M0;
	    	}
	    }
	    else if (Mds_min < 1.0)
	    {
	    	ds_max = ds_mid;
	    	break;
	    }
	    else if (Mds_max < 1.0)
	    {
	    	ds_min = ds_mid;
	    	break;
	    }
	    else
	    {
	    	dslen *= 0.5;
	    	ds_min = -dslen;	Mds_min_set = NO;
	    	ds_max = dslen;		Mds_max_set = NO;
	    	ds_mid = 0.0;		Mds_mid = M0;
	    }
	}
	if (i >= MAX_ITER)
	{
	    char s[80];
	    (void) sprintf(s,"Uable to find limits after %d iterations\n",i);
	    node_warning("set_limits_for_amn_find_root",s,"\n",flag);
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("diffraction","Left set_limits_for_amn_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return NO;
	}

	*pds_min = ds_min;		*pds_max = ds_max;
	*epsilon = 0.5 * SONIC_TOL;	*delta = 2.0*EPS*(ds_max - ds_min);
#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("diffraction","Left set_limits_for_amn_find_root()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return YES;
}		/*end set_limits_for_amn_find_root*/

LOCAL	boolean amn(
	double		ds,
	double		*M,
	POINTER		povrsht)
{
	OVERSHOOT	*ovrsht = (OVERSHOOT *) povrsht;
	Front		*fr = ovrsht->fr;
	POINT		*p0 = ovrsht->p0;
	POINT		*pc = ovrsht->pc;
	BOND		*b = ovrsht->b;
	CURVE		*c = ovrsht->c;
	Locstate	sa = ovrsht->sa, sa_opp = ovrsht->sa_opp;
	double		*abs_v0 = ovrsht->abs_v0;
	double		dt = ovrsht->dt;
	ORIENTATION	c_orient = ovrsht->c_orient;
	SIDE		inc_side = ovrsht->inc_side;
	BOND		*curr_b;
	CURVE		*curr_c;
	POINT		*ptmp;
	double		va[MAXD];
	double		d[MAXD], qa;
	double		abs_v[MAXD];
	double		para;
	double           *h = fr->rect_grid->h;
	int		dim = fr->rect_grid->dim;
	int		i;
	static	boolean	first = YES;
	static	BOND	*bdir = NULL;

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("diffraction","\nEntered amn(), ds = %g\n",ds);
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	if (first)
	{
		first = NO;
		scalar(&bdir,sizeof(BOND));
		bdir->start = Static_point(fr->interf);
		bdir->end = Static_point(fr->interf);
	}

	bond_secant_to_curve(p0,b,c,c_orient,bdir,fr,ds);

	curr_c = c;
	ptmp = (c_orient == POSITIVE_ORIENTATION) ? bdir->end : bdir->start;
	if (((c_orient == POSITIVE_ORIENTATION) && (ds >= 0.0)) ||
	    ((c_orient == NEGATIVE_ORIENTATION) && (ds < 0.0)))
	{
		curr_b = bdir->next;
		for (i = 0; i < dim; i++)
			Coords(pc)[i] = Coords(bdir->end)[i];
	}
	else
	{
		curr_b = bdir->prev;
		for (i = 0; i < dim; i++)
			Coords(pc)[i] = Coords(bdir->start)[i];
	}
	ovrsht->bc = (curr_c == c) ? curr_b :
				Bond_at_node(c,Opposite_orient(c_orient));


	ovrsht->tcr = para =
		(scaled_bond_length(curr_b,h,dim) < 0.001) ? /*TOLERANCE*/
		0.5 : separation(ptmp,curr_b->start,fr->interf->dim) /
			bond_length(curr_b);

	if (inc_side == NEGATIVE_SIDE)
	{
		left_state_along_bond(para,curr_b,curr_c,sa);
		right_state_along_bond(para,curr_b,curr_c,sa_opp);
	}
	else
	{
		right_state_along_bond(para,curr_b,curr_c,sa);
		left_state_along_bond(para,curr_b,curr_c,sa_opp);
	}

	for (i = 0; i < dim; i++)
	{
		d[i] = Coords(pc)[i] - Coords(p0)[i];
		abs_v[i] = abs_v0[i] + d[i]/dt;
		ovrsht->abs_v[i] = abs_v[i];
		va[i] = vel(i,sa) - abs_v[i];
	}

	qa = mag_vector(va,dim);
	*M = qa/sound_speed(sa);

#if defined(DEBUG_NODE_PROPAGATE)
	debug_print("diffraction","Left amn(), Mach number = %g\n\n",*M);
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	return FUNCTION_SUCCEEDED;
}		/*end amn*/

LOCAL int is_short_and_opp_unprop(
	O_CURVE		*oldc,
	O_CURVE		*newc,
	RECT_GRID	*gr)
{
	NODE		*n, *oppn, *nn;

	if (!oldc || !oldc->curve)
	    return NO;
	if (!newc || !newc->curve)
	    return NO;
	if (!is_short_curve(oldc->curve,oldc->orient,gr,1.0))
	    return NO;
	
	n = Node_of_o_curve(newc);
	oppn = Opp_node_of_o_curve(newc);
	if (propagation_status(oppn) == UNPROPAGATED_NODE)
	    return YES;
	if (propagation_status(oppn) == VEL_COMPUTED_NODE)
	{
	    for (nn = next_node(n); nn; nn = next_node(nn))
	    {
	    	if (propagation_status(nn) == UNPROPAGATED_NODE)
	    	    return YES;
	    }
	}
	
	return NO;
}		/*end is_short_and_opp_unprop*/

EXPORT	ANGLE_DIRECTION incident_shock_orientation(
	int		w_type,
	ORIENTATION	orient)
{
	if ((is_forward_wave(w_type) && orient == POSITIVE_ORIENTATION)
				 ||
	(is_backward_wave(w_type) && orient == NEGATIVE_ORIENTATION))
	{
		return CLOCKWISE;
	}
	else 
	{
		return COUNTER_CLOCK;
	}
}		/*end incident_shock_orientation*/

/*
*			  transmission_node_propagate():
*
*
*					    curve 1
*                   state 1 = state 2          / incident shock
*                                             /                        
* 					     /
*					    /
*					   /
*                                         /     state 0
*                                        / 
*                                      / 
*            curve 3 ----------------/---------- curve 0
*           contact(behind)        /       contact(ahead) 
*                                /    
*                               /       state 4
*                  state 3    /    
*                           / transmitted shock
*			curve 4
*
*	A transmission node is a node at which an incident shock meets a
*	contact discontinuity producing a transmitted shock but no reflected
*	waves.  The forward facing (low pressure) side of the incident shock
*	defines the front side of the wave and hence an incident to front
*	angular orientation of the node.  The propagation direction (ahead)
*	of the node also provides an angular orientation, the incident to
*	ahead direction.
*
*	Front and back refer to the low and high pressure sides of the 
*	incident shock, while ahead and behind refer to the direction
*	of node propagation.
*	The curve indices refer to the way the angles are stored in the RP,
*	and this is not consistent with the way other nodes are indexed.
*	TODO:  The incident should be curve 0.
*/


/* ARGSUSED */
EXPORT int transmission_node_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	BOND		*crossbinc;	/* intersecting two bonds */
	BOND		*crossbfront;	/* on newcinc, newcfront */
	COMPONENT	front_comp;
	O_CURVE		Oldcinc, Newcinc;	/* the incident curve */
	O_CURVE		Oldctrans, Newctrans;	/* the transmitted curve */
	O_CURVE		Oldcfront, Newcfront;	/* the front contact */
	O_CURVE		Oldcback, Newcback;	/* the back contact */
	RP_DATA		*RP;
	double		tcr_inc;  	/* fractional distances to cross */
	double		tcr_front;
	double		q[5];		/* steady flow speeds */
	ANGLE_DIRECTION	i_to_f_dir;	/* dir of angle to contact faced by */
					/* low pressure side of cinc  */
	boolean		c_ext[2];
	int		status;
	static POINT	*pc = NULL;		/* cross point */
	static double	**t = NULL;

	debug_print("transmission","Entered transmission_node_propagate()\n");
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("transmission")) 
	{
	    (void) printf("\n\tOLD NODE:\n");
	    print_node(oldn);
	    (void) printf("RP DATA at old node\n");
	    print_RP_DATA(Rp_data(oldn),Node_vel(oldn));
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */


		/* Allocate storage */

	if (pc == NULL) 
	{
	    pc = Static_point(fr->interf);
	    bi_array(&t,2,MAXD,FLOAT);
	}


		/* Identify curves */

	find_curve_with_status(oldn,&Oldcinc.curve,&Oldcinc.orient,INCIDENT);
	find_curve_with_status(oldn,&Oldctrans.curve,&Oldctrans.orient,
		TRANSMITTED);
	if (Oldcinc.curve == NULL || Oldctrans.curve == NULL)
	{
	    node_warning("transmission_node_propagate",
	                 "unable to locate incident or transmitted curve","\n",
			 flag);
	    debug_print("transmission","Left transmission_node_propagate(), ");
	    if (debugging("transmission"))
		print_node_status("status = ",ERROR_NODE,"\n");
	    return ERROR_NODE;
	}
	find_adjacent_curves(&Oldcinc,&i_to_f_dir,&Oldcfront,&Oldcback,
			     &front_comp);

	Check_return(
	    find_correspond_of_oriented_curve(&Oldcinc,&Newcinc,
					      newn,fr,newn->interface),
	    transmission_node_propagate)
	Check_return(
	    find_correspond_of_oriented_curve(&Oldcfront,&Newcfront,
					      newn,fr,newn->interface),
	    transmission_node_propagate)
	Check_return(
	    find_correspond_of_oriented_curve(&Oldcback,&Newcback,
					      newn,fr,newn->interface),
	    transmission_node_propagate)
	Check_return(
	    find_correspond_of_oriented_curve(&Oldctrans,&Newctrans,
					      newn,fr,newn->interface),
	    transmission_node_propagate)

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("transmission")) 
	{
		(void) printf("\t\tOLD INCIDENT CURVE:\n");
		print_o_curve(&Oldcinc);
		(void) printf("\t\tOLD TRANSMITTTED CURVE:\n");
		print_o_curve(&Oldctrans);
		(void) printf("\t\tOLD FRONT CONTACT:\n");
		print_o_curve(&Oldcfront);
		(void) printf("\t\tOLD BACK CONTACT:\n");
		print_o_curve(&Oldcback);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

		/* Identify new position of node */

	sav_normal = interface_normal_function(oldn->interface);
	set_temp_tnode_normal(&interface_normal_function(oldn->interface));
	status = cross_or_extend_to_cross_two_propagated_curves(&Oldcinc,
			&Newcinc,&Oldcfront,&Newcfront,&pc,
			&crossbinc,&crossbfront,
			&tcr_inc,&tcr_front,fr,(POINTER)wave,rp,
			dt,dt_frac,flag,c_ext);
	interface_normal_function(oldn->interface) = sav_normal;
	
	if (status != GOOD_NODE)
	{
	    debug_print("transmission","Left transmission_node_propagate(), ");
	    if (debugging("transmission"))
	        print_node_status("status = ",status,"\n");
	    return status;
	}

	RP = Rp_data(newn);
	RP->ang_dir = Opposite_ang_dir(i_to_f_dir);

		/* Find a direction for the new cfront */

	find_tangent_to_propagated_curve(pc,crossbfront,&Oldcfront,&Newcfront,
		t[1],fr,(POINTER)wave,dt);

	if (trans_node_parameter() == USE_INCIDENT_ANGLE) 
	{
			/* Find direction of new cinc */
		find_tangent_to_propagated_curve(pc,crossbinc,&Oldcinc,
						 &Newcinc,t[0],
						 fr,(POINTER)wave,dt);
	}

		/* Read off states from the propagated incident curve */
	if (curve_ang_oriented_l_to_r(i_to_f_dir,Newcinc.orient))
	{
		right_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,
			RP->state[0]);
		left_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,
			RP->state[1]);
	}
	else 
	{
		left_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,
			RP->state[0]);
		right_state_along_bond(tcr_inc,crossbinc,Newcinc.curve,
			RP->state[1]);
	}
	ft_assign(RP->state[2],RP->state[1],fr->sizest);


		/* Read off state at node on back contact */

	if (curve_ang_oriented_l_to_r(i_to_f_dir,Newcback.orient))
	{
		ft_assign(RP->state[3],Left_state_at_node_of_o_curve(&Oldcback),
			fr->sizest);
	}
	else 
	{
		ft_assign(RP->state[3],Right_state_at_node_of_o_curve(&Oldcback),
			fr->sizest);
	}

		/* Read off state from propagated front contact */

	if (curve_ang_oriented_l_to_r(i_to_f_dir,Newcfront.orient))
	{
	    right_state_along_bond(tcr_front,crossbfront,Newcfront.curve,
	    	                   RP->state[4]);
	}
	else 
	{
	    left_state_along_bond(tcr_front,crossbfront,Newcfront.curve,
	    	                  RP->state[4]);
	}
	

		/* Modify the interface and assign the new states */

#if defined(DEBUG_NODE_PROPAGATE)
	if(debugging("transmission"))
	{
	    (void) printf("States into find_transmission_node_states()\n");
	    verbose_print_state("state0",RP->state[0]);
	    verbose_print_state("state1",RP->state[1]);
	    verbose_print_state("state2",RP->state[2]);
	    verbose_print_state("state3",RP->state[3]);
	    verbose_print_state("state4",RP->state[4]);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	status = find_transmission_node_states(q,t,RP,trans_node_parameter(),
					       WEAK,flag);

	switch (status)
	{
	case REGULAR_TRANSMISSION:
	    break;

	case BIFURCATION_TRANSMISSION:
	    status = init_transmission_node_deprecursion(oldn,newn,&Oldctrans,
							 &Newctrans,&Oldcback,
							 &Newcback,dt,*dt_frac,
							 fr,wave,rp);
	    propagation_status(newn) = VEL_COMPUTED_NODE;
#if defined(DEBUG_NODE_PROPAGATE)
	    debug_print("transmission","Left transmission_node_propagate()\n");
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    return status;
	default:
	    status = ERROR_NODE;
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("transmission"))
	    {
	        node_warning("transmission_node_propagate",
		             "find_transmission_node_states() failed","\n",
		    	 flag);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    debug_print("transmission","Left transmission_node_propagate()\n");
	    return status;
	}

	if (!modify_transmission_node(pc,crossbinc,crossbfront,oldn,newn,
		                         &Oldcinc,&Newcinc,&Oldcfront,
					 &Newcfront,&Oldctrans,&Newctrans,
		                         &Oldcback,&Newcback,i_to_f_dir,
					 fr,wave,dt,RP,flag)) 
	{
	    node_warning("transmission_node_propagate",
		         "modify_transmission_node() failed","\n",flag);
	    status = ERROR_NODE;
	    if (debugging("transmission"))
	        print_node_status("status = ",status,"\n");
	    debug_print("transmission","Left transmission_node_propagate()\n");
	    return status;
	}

	status = GOOD_NODE;

	propagation_status(newn) = PROPAGATED_NODE;

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("transmission")) 
	{
		(void) printf("\t\tNEW INCIDENT CURVE:\n");
		print_o_curve(&Newcinc);
		(void) printf("\t\tNEW TRANSMITTTED CURVE:\n");
		print_o_curve(&Newctrans);
		(void) printf("\t\tNEW FRONT CONTACT:\n");
		print_o_curve(&Newcfront);
		(void) printf("\t\tNEW BACK CONTACT:\n");
		print_o_curve(&Newcback);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	debug_print("transmission","Left transmission_node_propagate(), ");
	if (debugging("transmission"))
	    print_node_status("status = ",status,"\n");
	return status;
}		/*end transmission_node_propagate*/

EXPORT	void	set_temp_tnode_normal(
	NORMAL_FUNCTION *nf)
{
	static const char *nname = "temp_tnode_normal";
	nf->_normal = temp_tnode_normal;
	nf->_normal_name = nname;
}		/*end set_temp_tnode_normal*/

/*
*			temp_tnode_normal():
*
*	This function is used to help enforce the angle of the incident
*	curve at a transmission node.  The transmission node is unusual in
*	that the angle of the incident is set via the transmission node
*	computation.  It is unfeasible to insert a point far enough away
*	from the node that the tangent functions will compute the correct
*	angle.  Thus for this one curve (incident) at a transmission node,
*	we read the angle from the old RP structure.  Otherwise we use the
*	normal function saved above in the global (to this file) variable
*	sav_normal.
*/

LOCAL void temp_tnode_normal(
	POINT		   *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	   *hs,
	double		   *nor,
	Front		   *front)
{
	CURVE		*c = Curve_of_hs(hs);
	double		inc_ang;
	ORIENTATION	orient;

	orient = (separation(p,c->start->posn,front->rect_grid->dim) == 0.0) ?
			POSITIVE_ORIENTATION : NEGATIVE_ORIENTATION;

	if ((status_at_node(c,orient) == INCIDENT) &&
	    (node_type(Node_of(c,orient)) == TRANSMISSION_NODE))
	{
	    inc_ang = Rp_data(Node_of(c,orient))->ang[1];
	    if (orient == POSITIVE_ORIENTATION)
	    {
	    	nor[0] = sin(inc_ang);
	    	nor[1] = -cos(inc_ang);
	    }
	    else
	    {
	    	nor[0] = -sin(inc_ang);
	    	nor[1] = cos(inc_ang);
	    }
	}
	else
	    (*sav_normal._normal)(p,hse,hs,nor,front);
	
}		/*end temp_tnode_normal*/

LOCAL int tran_node_parameter_choice = USE_SLIP;

EXPORT int trans_node_parameter(void)
{
	return tran_node_parameter_choice;
}		/*end trans_node_parameter*/

EXPORT	int set_trans_node_parameter(
	int		value)
{
	tran_node_parameter_choice = value;
	return value;
}		/*end set_trans_node_parameter*/

/*
*			init_transmission_node_deprecursion():
*
*	Initializes the rp structure for deprecursion from precursor with
*	reflected rarefaction back to a single diffraction node.  We allow
*	for the case of an untracked overtake node.
*/

LOCAL int init_transmission_node_deprecursion(
	NODE		*old_tnode,
	NODE		*new_tnode,
	O_CURVE		*oldctrans,
	O_CURVE		*newctrans,
	O_CURVE		*oldcback,
	O_CURVE		*newcback,
	double		dt,
	double		dt_frac,
	Front		*fr,
	Wave		*wave,
	RPROBLEM 	**rp)
{
	NODE		*new_cnode;		/* CROSS_NODE */
	NODE		*old_cnode;
	NODE		*new_irnode;		/* TOT_INT_REFL_NODE */
	NODE		*old_irnode;
	NODE		*new_onode;		/* OVERTAKE_NODE */
	NODE		*old_onode;
	NODE		*interact_nodes[9];
	O_CURVE		Ctmp0, Ctmp1;		/* raref edges at irnode */
	RECT_GRID	*gr = fr->rect_grid;

	debug_print("transmission",
	      "Entered init_transmission_node_deprecursion()\n");

		/* Identify cross and total internal reflection nodes */


	new_cnode = Opp_node_of_o_curve(newctrans);
	old_cnode = Opp_node_of_o_curve(oldctrans);
	new_irnode = Opp_node_of_o_curve(newcback);
	old_irnode = Opp_node_of_o_curve(oldcback);

	if ((node_type(new_cnode) != CROSS_NODE) ||
	    (node_type(old_cnode) != CROSS_NODE) ||
	    !is_short_curve(oldctrans->curve,oldctrans->orient,gr,1.0))
	{
	    if (debugging("transmission"))
	    {
	        (void) printf("Unable to find cross node\n");
	        print_node_status("status = ",ERROR_NODE,"\n");
	    }
	    debug_print("transmission",
		  "Left init_transmission_node_deprecursion()\n");
	    return ERROR_NODE;
	}

	if ((node_type(new_irnode) != TOT_INT_REFL_NODE) ||
	    (node_type(old_irnode) != TOT_INT_REFL_NODE) ||
	    !is_short_curve(oldcback->curve,oldcback->orient,gr,1.0))
	{
	    if (debugging("transmission"))
	    {
	    	(void) printf("Unable to find total internal refl node\n");
	    	print_node_status("status = ",ERROR_NODE,"\n");
	    }
	    debug_print("transmission",
		  "Left init_transmission_node_deprecursion()\n");
	    return ERROR_NODE;
	}

		/* Identify new overtake node */

	new_onode = NULL;
	identify_curves_with_status(new_irnode,&Ctmp0,&Ctmp1,REFLECTED);
	if (is_rarefaction_leading_edge(wave_type(Ctmp0.curve)))
	{
		new_onode = Opp_node_of_o_curve(&Ctmp0);
	}
	else if (is_rarefaction_leading_edge(wave_type(Ctmp1.curve)))
	{
		new_onode = Opp_node_of_o_curve(&Ctmp1);
		if ((node_type(new_onode) != OVERTAKE_NODE) ||
		    !is_short_curve(Ctmp1.curve,Ctmp1.orient,gr,1.0))
		{
			new_onode = NULL;
		}
	}

		/* Identify new overtake node */

	identify_curves_with_status(old_irnode,&Ctmp0,&Ctmp1,REFLECTED);
	if (is_rarefaction_leading_edge(wave_type(Ctmp0.curve)))
	{
		old_onode = Opp_node_of_o_curve(&Ctmp0);
		if (node_type(old_onode) != OVERTAKE_NODE)
		{
			old_onode = NULL;
		}
	}
	else if (is_rarefaction_leading_edge(wave_type(Ctmp1.curve)))
	{
		old_onode = Opp_node_of_o_curve(&Ctmp1);
		if (node_type(old_onode) != OVERTAKE_NODE)
		{
			old_onode = NULL;
		}
	}
	if (((old_onode != NULL) && (new_onode == NULL)) ||
	    ((old_onode == NULL) && (new_onode != NULL)))
	{
	    if (debugging("transmission"))
	    {
	    	(void) printf("Old and new overtake nodes inconsistent\n");
	    	print_node_status("status = ",ERROR_NODE,"\n");
	    	print_node(old_onode);
	    	print_node(new_onode);
	    }
	    debug_print("transmission",
	          "Left init_transmission_node_deprecursion()\n");
	    return ERROR_NODE;
	}

		/* Initialize rp */

	interact_nodes[0] = new_irnode;
	interact_nodes[1] = old_irnode;
	interact_nodes[2] = new_cnode;
	interact_nodes[3] = old_cnode;
	interact_nodes[4] = new_tnode;
	interact_nodes[5] = old_tnode;
	if (old_onode != NULL)
	{
		interact_nodes[6] = new_onode;
		interact_nodes[7] = old_onode;
		interact_nodes[8] = NULL;
	}
	else
	{
		interact_nodes[6] = NULL;
	}
	augment_rproblem_list(rp,interact_nodes,dt,dt_frac,
		              old_tnode->interface,new_tnode->interface,
		              fr,(POINTER)wave);

	if (debugging("transmission"))
	    print_node_status("status = ",CROSS_NODE_NODE,"\n");
	debug_print("transmission","Left init_transmission_node_deprecursion()\n");
	return CROSS_NODE_NODE;
}		/*end init_transmission_node_deprecursion*/


LOCAL int modify_transmission_node(
	POINT		*pc,
	BOND		*newbinc,
	BOND		*newbfront,
	NODE		*oldn,
	NODE		*newn,
	O_CURVE		*oldcinc,
	O_CURVE		*newcinc,
	O_CURVE		*oldcfront,
	O_CURVE		*newcfront,
	O_CURVE		*oldctrans,
	O_CURVE		*newctrans,
	O_CURVE		*oldcback,
	O_CURVE		*newcback,
	ANGLE_DIRECTION	i_to_f_dir,
	Front		*fr,
	Wave		*wave,
	double		dt,
	RP_DATA		*RP,
	NODE_FLAG	flag)
{
	Locstate	st_l, st_r;

	debug_print("transmission","Entered modify_transmission_node()\n");

		/* Propagate old front contact */

	Coords(newn->posn)[0] = Coords(pc)[0];
	Coords(newn->posn)[1] = Coords(pc)[1];
	if (curve_ang_oriented_l_to_r(i_to_f_dir,newcfront->orient))
	{
	    st_l = RP->state[0];
	    st_r = RP->state[4];
	}
	else 
	{
	    st_l = RP->state[4];
	    st_r = RP->state[0];
	}
	if (!propagate_curve_near_node(oldn,newn,oldcfront,newcfront,
			                  newbfront,st_l,st_r,NO,
					  RP->ang[0],fr,wave,dt,flag)) 
	{
	    node_warning("modify_transmission_node",
	                 "propagate_curve_near_node() failed","\n",flag);
	    debug_print("transmission","Left modify_transmission_node()\n");
	    return NO;
	}

		/* Propagate old incident curve */

	if ((newcinc->orient==POSITIVE_ORIENTATION && i_to_f_dir==CLOCKWISE)
			 || 
	(newcinc->orient==NEGATIVE_ORIENTATION && i_to_f_dir==COUNTER_CLOCK))
	{
	    st_l = RP->state[1];
	    st_r = RP->state[0];
	}
	else 
	{
	    st_l = RP->state[0];
	    st_r = RP->state[1];
	}
	if (!propagate_curve_near_node(oldn,newn,oldcinc,newcinc,
			                  newbinc,st_l,st_r,YES,RP->ang[1],
					  fr,wave,dt,flag)) 
	{
	    node_warning("modify_transmission_node",
	                 "propagate_curve_near_node() failed","\n",flag);
	    debug_print("transmission","Left modify_transmission_node()\n");
	    return NO;
	}


		/* Propagate old transmitted curve */

	if ((newctrans->orient==NEGATIVE_ORIENTATION && i_to_f_dir==CLOCKWISE)
			 || 
	(newctrans->orient==POSITIVE_ORIENTATION && i_to_f_dir==COUNTER_CLOCK))
	{
	    st_l = RP->state[3];
	    st_r = RP->state[4];
	}
	else 
	{
	    st_l = RP->state[4];
	    st_r = RP->state[3];
	}
	if (!propagate_curve_near_node(oldn,newn,oldctrans,newctrans,
		                          Bond_at_node_of_o_curve(newctrans),
		                          st_l,st_r,YES,RP->ang[4],
					  fr,wave,dt,flag)) 
	{
	    node_warning("modify_transmission_node",
	                 "propagate_curve_near_node() failed","\n",flag);
	    debug_print("transmission","Left modify_transmission_node()\n");
	    return NO;
	}

		/* Propagate old slipline and insert a new bond */

	if (curve_ang_oriented_l_to_r(i_to_f_dir,newcback->orient))
	{
	    st_l = RP->state[3];
	    st_r = RP->state[2];
	}
	else
	{
	    st_l = RP->state[2];
	    st_r = RP->state[3];
	}
	if (!propagate_curve_near_node(oldn,newn,oldcback,newcback,
					  Bond_at_node_of_o_curve(newcback),
					  st_l,st_r,YES,RP->ang[3],
					  fr,wave,dt,flag)) 
	{
	    node_warning("modify_transmission_node",
	                 "propagate_curve_near_node() failed","\n",flag);
	    debug_print("transmission","Left modify_transmission_node()\n");
	    return NO;
	}
	debug_print("transmission","Left modify_transmission_node()\n");
	return YES;
}		/*end modify_transmission_node*/
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */
