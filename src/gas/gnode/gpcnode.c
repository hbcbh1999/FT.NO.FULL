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
*				gpcnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Contains the routine used in the first advance of nodes
*	corresponding to shock wave-contact discontinuity interactions.
*
*/

#if defined(TWOD)

#include <gdecs/gdecs.h>

/* Values returned by transmit_precursor() */

enum _TRANSMIT_PRECURSOR {
	PRECURSOR_TRANSMITTED	     = 1,
	CROSS_POSITION_NOT_FOUND,
	CROSS_POSITION_BUT_NO_STATES,
	SHOCK_DIFFRACTION_NOT_RESET,
	CROSS_NODE_NOT_INSTALLED
};
typedef enum _TRANSMIT_PRECURSOR TRANSMIT_PRECURSOR;

	/* LOCAL Function Declarations */
LOCAL	Locstate	find_state_behind_precursor(Front*,Wave*,double,NODE*,
					            O_CURVE**,O_CURVE**);
LOCAL	TRANSMIT_PRECURSOR transmit_precursor(Front*,Wave*,NODE*,double,
					      NODE_FLAG,O_CURVE**,O_CURVE**,
					      BOND**,RP_DATA*,RP_DATA*,POINT*);
LOCAL	NODE	*install_new_transmission_node(Front*,Wave*,double,NODE*,
					       POINT*,BOND*,O_CURVE**,
					       O_CURVE**,BOND**,RP_DATA*);
LOCAL	NODE	*new_transmission_node(NODE*,O_CURVE**,O_CURVE**,BOND**,
				       POINT*,BOND*,double**,RP_DATA*,RP_DATA*,
				       Front*,Wave*,double);
LOCAL	int	detach_precursor(Front*,Wave*,NODE*,RPROBLEM**,double,double*,
				 NODE_FLAG,O_CURVE**,O_CURVE**,BOND**,RP_DATA*,
				 RP_DATA**,POINT**,NODE**);
LOCAL	int	find_new_cross_node_position(Front*,Wave*,double,NODE_FLAG,
					     O_CURVE**,O_CURVE**,RP_DATA*,
					     double*,BOND**,POINT*,POINT**);
LOCAL	int	find_overtake_node_position(BOND**,CURVE*,NODE*,POINT*,
					    ORIENTATION,Front*);
LOCAL	int	new_cross_node(POINT*,BOND*,O_CURVE**,BOND**,RP_DATA*,
			       RP_DATA*,double*,Front*);
LOCAL	int	new_transmission_node_states_and_position(NODE*,O_CURVE**,
							  O_CURVE**,POINT*,
							  BOND**,Front*,
							  Wave*,double**,
							  RP_DATA*,RP_DATA*,
							  RPROBLEM**,double,
							  double*,NODE_FLAG);
LOCAL	int	precursor_shock_rr_reconfigure(Front*,Wave*,RPROBLEM**,
					       O_CURVE**,O_CURVE**,NODE*,
					       RP_DATA*,ANGLE_DIRECTION,double);
LOCAL	int	reset_shock_diffraction(NODE*,O_CURVE**,POINT*,POINT*,
					RP_DATA*,RP_DATA*,RP_DATA*,Front*);
LOCAL	void	final_shock_at_precursor_rr(NODE*,O_CURVE**,O_CURVE**,Front*);

/*
*			precursor_shock_rr_propagate():
*	
*	This function solves for the slow-fast anomalous transmission.
*	The Mach number behind the transmitted shock is subsonic causing
*	the transmitted shock  to detach from the original node.  
*
*	This bifurcation causes a change in topology.  The old transmitted
*	wave moves ahead on the incident contact, creating a new transmission
*	node.  This in turn scatters off a new transmitted wave which interacts
*	with the old incident, creating a cross node.  This slip line at
*	this node is not tracked.  One reflected wave at this cross node
*	is directed back toward the old behind contact, creating a total
*	internal reflection node with the old reflected rarefaction.  It is
*	this node that is returned in newn, with corresponding changes in
*	RP and newc.  The other reflected wave at the cross node is directed
*	toward the reflected rarefaction, creating an overtake node on the
*	the leading edge.
*/

EXPORT	int precursor_shock_rr_propagate(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*newn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	POINT		*pc,
	RP_DATA		*RP,
	BOND		**newb)
{
	RP_DATA		*tRP;		/* RP for new transmission node */
	RP_DATA		*cRP;		/* RP for new cross node */
	POINT		*ptr;		/* location of transmission node */
	NODE		*tnode;		/* new transmission node */
	int		status = ERROR_NODE;
	int		i, dim = fr->rect_grid->dim;

	debug_print("precursor","Entered precursor_shock_rr_propagate()\n");

	if (untrack_scattered_waves_at_refraction_precursor(flag) == YES)
	{
	    /*
	    * Turn off tracking at all waves but the incident shock
	    * and the contact
	    */

	    (void) printf("WARNING in precursor_shock_rr_propagate(), "
			  "untracking scattered curves\n");

	    if (oldc[1]->curve != NULL)
	        untracked_hyper_surf(oldc[1]->curve) = YES;
	    if (oldc[2]->curve != NULL)
	        untracked_hyper_surf(oldc[2]->curve) = YES;
	    if (oldc[3]->curve != NULL)
	        untracked_hyper_surf(oldc[3]->curve) = YES;
	    if (oldc[5]->curve != NULL)
	        untracked_hyper_surf(oldc[5]->curve) = YES;
	    node_type(oldn) = TOT_INT_REFL_NODE;
	    return REPEAT_TIME_STEP_NODE;
	}

	if (debugging("precursor"))
	{
	    (void) printf("WARNING in diffraction_node_propagate(), "
	                  "\tbifurcation from regular diffraction to "
	                  "precursor with reflected rarefaction.\n");
	    print_angle_direction("RP->ang_dir =",RP->ang_dir,"\n");
	}
	if (oldn == NULL)
	{
	    (void) printf("WARNING in precursor_shock_rr_propagate(), "
	                  "double bifurcation in rproblem\n");
	    debug_print("precursor","Left precursor_shock_rr_propagate()\n");
	    return status;
	}

	status = detach_precursor(fr,wave,oldn,rp,dt,dt_frac,
				  flag,oldc,newc,newb,RP,&tRP,&ptr,&tnode);

	if (status != GOOD_NODE)
	{
	    (void) printf("WARNING in precursor_shock_rr_propagate(), "
	                  "detach_precursor() failed\n");
	    debug_print("precursor","Left precursor_shock_rr_propagate()\n");
	    return status;
	}

	switch (transmit_precursor(fr,wave,newn,dt,flag,oldc,newc,newb,
				   RP,tRP,ptr))
	{
	case PRECURSOR_TRANSMITTED:
	    status = GOOD_NODE;
	    break;

	case CROSS_POSITION_BUT_NO_STATES:
	    (void) printf("WARNING in precursor_shock_rr_propagate(), "
	                  "new cross node doesn't exist\n");
	    status = precursor_shock_rr_reconfigure(fr,wave,rp,newc,oldc,
					            tnode,tRP,RP->ang_dir,dt);
	    debug_print("precursor","Left precursor_shock_rr_propagate()\n");
	    return status;

	case CROSS_POSITION_NOT_FOUND:
	case SHOCK_DIFFRACTION_NOT_RESET:
	case CROSS_NODE_NOT_INSTALLED:
	default:
	    status = ERROR_NODE;
	    (void) printf("WARNING in precursor_shock_rr_propagate(), "
	                  "can't install transmitted precursor\n");
	    return status;
	}

	for (i = 0; i < dim; i++)
	    Coords(pc)[i] = Coords(newn->posn)[i];

	newc[5]->curve = NULL;

	status = modify_diffraction_node(pc,oldn,newn,oldc,newc,newb,
					 fr,wave,dt,RP,flag);
	if (status != GOOD_NODE)
	{
	    debug_print("precursor","Left precursor_shock_rr_propagate()\n");
	    return status;
	}

	propagation_status(newn) = PROPAGATED_NODE;
	
	cRP = Rp_data(Opp_node_of_o_curve(newc[0]));
	if (track_scattered_wave(CROSS_NODE,SHOCK_WAVE,REFLECTED,
				 cRP->state[1],cRP->state[2],fr))
	    final_shock_at_precursor_rr(newn,oldc,newc,fr);

	debug_print("precursor","Left precursor_shock_rr_propagate()\n");
	return status;
}		/*end precursor_shock_rr_propagate*/

/*
*			detach_precursor():
*
*	This function installs the new transmission node in the
*	bifurcation from regular diffraction to precursor with reflected
*	rarefaction.  The states and angles are found (in tRP), the node 
*	position, and the transmission node is installed as far as possible.
*	The new transmitted wave and behind contact are not handled.
*/

LOCAL	int detach_precursor(
	Front		*fr,
	Wave		*wave,
	NODE		*oldn,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	RP_DATA		*RP,
	RP_DATA		**tRP,
	POINT		**ptr,
	NODE		**tnode)
{
	BOND		*btr;		/* ptr lies here on old ahead cont */
	int		status;
	static RP_DATA	*sRP = NULL;
	static double	**t = NULL;	/* t[1] ahead contact tangent */
				        /* t[0] new incident tangent */

	debug_print("precursor","Entered detach_precursor()\n");

	if (sRP == NULL)
	{
	    sRP = allocate_RP_DATA_structure(fr->sizest,NO,GAS_STATE);
	    bi_array(&t,2,MAXD,FLOAT);
	}
	*tRP = (newc[5]->curve != NULL) ?
		allocate_RP_DATA_structure(fr->sizest,YES,GAS_STATE) : sRP;

	*ptr = Point(Coords(oldn->posn));
	status = new_transmission_node_states_and_position(oldn,oldc,newc,*ptr,
							   &btr,fr,wave,t,RP,
							   *tRP,rp,dt,dt_frac,
							   flag);

	if (status != GOOD_NODE)
	    return status;

	*tnode = new_transmission_node(oldn,oldc,newc,newb,
				       *ptr,btr,t,RP,*tRP,fr,wave,dt);

	debug_print("precursor","Left detach_precursor()\n");
	return status;

}		/*end detach_precursor*/


/*
*		new_transmission_node_states_and_position():
*
*	This function computes the position and states around the new
*	transmission node.  The new node position is found by computing
*	the steady via (118.12) assuming the behind state is sonic.  We
*	then intersect the old ahead contact with a circle of this radius,
*	giving the transmission node position and the angles of the new
*	ahead contact.  A point is inserted in the old ahead contact at this
*	intersection point.  States are then computed using a call to 
*	find_transmission_node_states().
*/


LOCAL int new_transmission_node_states_and_position(
	NODE		*oldn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	POINT		*ptr,
	BOND		**btr,
	Front		*fr,
	Wave		*wave,
	double		**t,
	RP_DATA		*RP,
	RP_DATA		*tRP,
	RPROBLEM	**rp,
	double		dt,
	double		*dt_frac,
	NODE_FLAG	flag)
{
	CURVE		*c;
	POINT		*oldp;
	INTERFACE	*intfc = newc[0]->curve->interface;
	Locstate	spa, spb;
	double		m, qa, cbsqr;
	double		pa, pb, Va, Vb;
	double		radius;
	double		beta;
	double		tcr;
	double		q[5];
	ORIENTATION	orient;
	int		i, dim = fr->rect_grid->dim;
	int		status = ERROR_NODE;
	boolean		sav_intrp;
	static const	int	MAX_NUM_ITER = 1;
	static POINT	*pcen = NULL;

	debug_print("precursor",
		"Entered  new_transmission_node_states_and_position()\n");
	if (pcen == NULL)
	    pcen = Static_point(fr->interf);

	sav_intrp = interpolate_intfc_states(intfc);
	spa = RP->state[6];
	tRP->ang_dir = Opposite_ang_dir(RP->ang_dir);
	spb = find_state_behind_precursor(fr,wave,dt,oldn,oldc,newc);
	ft_assign(tRP->state[0],RP->state[6],fr->sizest);
	ft_assign(tRP->state[1],spb,fr->sizest);
	ft_assign(tRP->state[2],tRP->state[1],fr->sizest);
	ft_assign(tRP->state[4],RP->state[0],fr->sizest);

	cbsqr = sound_speed_squared(spb);
	pa = pressure(spa);	pb = pressure(spb);
	Va = 1.0/Dens(spa);	Vb = 1.0/Dens(spb);
	m = mass_flux(pb,spa);
	qa = sqrt(cbsqr + (pb - pa)*(Vb + Va));
	oldp = oldn->posn;
	for (i = 0; i < dim; i++)
	    Coords(pcen)[i] = Coords(oldp)[i] + vel(i,spa)*dt;
	if (debugging("precursor"))
	{
	    (void) printf("pcen = <%g, %g>\n",Coords(pcen)[0],Coords(pcen)[1]);
	    (void) printf("newc[6]\n");
	    print_o_curve(newc[6]);
	}

	*btr = NULL;
	c = newc[6]->curve;
	orient = newc[6]->orient;
	for (i = 0; i < MAX_NUM_ITER; i++)
	{
	    beta = asin(Va*m/qa);
	    if (debugging("precursor"))
	    {
	        (void) printf("Iteration %d\n",i);
	        (void) printf("beta = %g, (%g degrees)\n",
	    		  beta,degrees(beta));
	    }
	    radius = dt * qa * 1.05;/*TOLERANCE*/
	    if (debugging("precursor"))
	    	(void) printf("qa = %g, radius = %g\n",qa,radius);
	    status = crossing_of_a_propagated_curve_and_circle(oldc[6],
				newc[6],radius,pcen,ptr,btr,&tcr,fr,
				(POINTER)wave,rp,dt,dt_frac,flag);
	    if (status != GOOD_NODE)
		return status;
	    if (debugging("precursor"))
	    {
	    	(void) printf("ptr = <%g, %g>\n",Coords(ptr)[0],Coords(ptr)[1]);
	    	(void) printf("btr\n");	print_bond(*btr);
	    }

	    interpolate_intfc_states(intfc) = NO;
	    if (*btr != NULL)
	        (void) delete_start_of_bond((*btr)->next,c);
 	    if (insert_point_in_bond(ptr,*btr,c) != FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in new_transmission_node_states_and_position(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	    find_tangent_to_curve(ptr,*btr,c,orient,t[1],fr);
	    /* Note: computing t[0] is not needed for USE_SLIP */
	    if (RP->ang_dir == COUNTER_CLOCK)
	    {
	    	t[0][0] = cos(beta)*t[1][0] + sin(beta)*t[1][1];
	    	t[0][1] = cos(beta)*t[1][1] - sin(beta)*t[1][0];
	    }
	    else
	    {
	    	t[0][0] = cos(beta)*t[1][0] - sin(beta)*t[1][1];
	    	t[0][1] = cos(beta)*t[1][1] + sin(beta)*t[1][0];
	    }

	    (void) find_transmission_node_states(q,t,tRP,USE_SLIP,WEAK,flag);
	    qa = q[0];
	}


	interpolate_intfc_states(intfc) = sav_intrp;
	debug_print("precursor",
		"Left  new_transmission_node_states_and_position()\n");
	return status;
}		/*end new_transmission_node_states_and_position*/


/*
*			find_state_behind_precursor():
*
*	This function finds the states behind the new incident (old
*	transmitted) wave at the new transmission node.  If the old
*	transmitted wave was tracked, this is simply a call to 
*	point_propagate().  Otherwise, a mock point_propagate is
*	performed using the computed angle of the old transmitted wave,
*	finding states along that normal, and calling npt_w_speed().
*/

LOCAL	Locstate find_state_behind_precursor(
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE		*oldn,
	O_CURVE		**oldc,
	O_CURVE		**newc)
{
	COMPONENT	comp;
	Locstate	sl, sr;
	POINT		*oldp;
	RP_DATA		*RP;
	double		V[MAXD];
	double		*nor;
	ORIENTATION	orient;
	static Locstate	sll = NULL, srr = NULL;
	static POINT	*newp = NULL;
	static WSSten	*sten = NULL;

	if (srr == NULL)
	{
	    newp = Static_point(fr->interf);
	    alloc_state(fr->interf,&sll,fr->sizest);
	    alloc_state(fr->interf,&srr,fr->sizest);
	    sten = AllocDefaultWSSten(fr);
	}	
	else
	    ClearWSStenData(sten);
	sten->front = fr;
	sten->wave = wave;
	sten->w_type = FORWARD_SHOCK_WAVE;
	sten->pjump = 0.0;
	sten->hs = NULL;
	sten->dt = dt;
	nor = sten->nor;

	oldp = oldn->posn;
	RP = Rp_data(oldn);
	if (oldc[5]->curve != NULL)
	{
	    point_propagate(fr,(POINTER)wave,oldp,newp,
			    Bond_at_node_of_o_curve(oldc[5]),oldc[5]->curve,
			    dt,V);
	    orient = oldc[5]->orient;
	}
	else
	{
	    comp = (curve_ang_oriented_l_to_r(RP->ang_dir,newc[6]->orient)) ?
			negative_component(newc[6]->curve) :
			positive_component(newc[6]->curve);
	    if (RP->ang_dir == COUNTER_CLOCK)
	    {
	    	orient = NEGATIVE_ORIENTATION;
	    	nor[0] = -sin(RP->ang[5]);
	    	nor[1] =  cos(RP->ang[5]);
	    	sl = RP->state[6];
	    	sr = RP->state[5];
	    }
	    else
	    {
	    	orient = POSITIVE_ORIENTATION;
	    	nor[0] =  sin(RP->ang[5]);
	    	nor[1] = -cos(RP->ang[5]);
	    	sl = RP->state[5];
	    	sr = RP->state[6];
	    }
	    states_near_location(sten,Coords(oldn->posn),nor,comp,comp,sl,sr);
	    npt_w_speed(sten,left_state(newp),right_state(newp),V);
	}
	return (curve_ang_oriented_l_to_r(RP->ang_dir,orient)) ?
		left_state(newp) : right_state(newp);
}		/*end find_state_behind_precursor*/


/*
*			new_transmission_node():
*
*	This function creates the new transmission node via 
*	install_new_transmission_node(), and then adjusts the angle of
*	the new incident.  This is necessary since USE_SLIP was used
*	to compute the node configuration (see above).
*/

LOCAL NODE *new_transmission_node(
	NODE		*oldn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	POINT		*p,
	BOND		*btr,
	double		**t,
	RP_DATA		*RP,
	RP_DATA		*tRP,
	Front		*fr,
	Wave		*wave,
	double		dt)
{
	CURVE		*c;
	INTERFACE	*intfc = newc[0]->curve->interface;
	Locstate	sl, sr;
	NODE		*tnode;
	double		p1, r0, q0;
	double		dir[MAXD];
	double		len;
	double		dv[MAXD];
	ORIENTATION	orient;
	int		i, dim = fr->rect_grid->dim;
	boolean		sav_intrp;

	debug_print("precursor","Entered new_transmission_node()\n");
	if (debugging("precursor"))
	{
		(void) printf("New transmission node at <%g, %g>\n",
			Coords(p)[0],Coords(p)[1]);
	}
	sav_intrp = interpolate_intfc_states(intfc);
	if (newc[5]->curve != NULL)
	{
		tnode = install_new_transmission_node(fr,wave,dt,oldn,
				p,btr,oldc,newc,newb,tRP);
	}
	else
		tnode = NULL;

	if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[6]->orient))
		assign_interacting_states(p,newc[6]->curve,newc[6]->orient,
			fr,tRP->state[1],tRP->state[3]);
	else
		assign_interacting_states(p,newc[6]->curve,newc[6]->orient,
			fr,tRP->state[3],tRP->state[1]);
	if (newc[5]->curve != NULL)
	{
		c = newc[5]->curve;
		orient = newc[5]->orient;
		sl = Left_state_at_node(c,orient);
		sr = Right_state_at_node(c,orient);
		if (curve_ang_oriented_l_to_r(tRP->ang_dir,orient))
		{
			ft_assign(sl,tRP->state[0],fr->sizest);
			ft_assign(sr,tRP->state[1],fr->sizest);
		}
		else
		{
			ft_assign(sl,tRP->state[1],fr->sizest);
			ft_assign(sr,tRP->state[0],fr->sizest);
		}
		p1 = pressure(tRP->state[1]);
		r0 = Dens(tRP->state[0]);
		q0 = mass_flux(p1,tRP->state[0])/
			(r0*fabs(sin(tRP->ang[0] - tRP->ang[1])));
		for (i = 0; i < dim; i++)
		{
			Node_vel(tnode)[i] = q0*t[1][i];
			dv[i] = Node_vel(oldn)[i] - Node_vel(tnode)[i];
		}
		len = 1.1*mag_vector(dv,dim)*dt;		/*TOLERANCE*/

		/*len = 0.1*Front_spacing(fr,GENERAL_WAVE);*/	/*TOLERANCE*/
		
		dir[0] = cos(tRP->ang[1]);
		dir[1] = sin(tRP->ang[1]);
		interpolate_intfc_states(intfc) = YES;
		(void) adjust_angle_at_node(tnode,oldc[5],
				newc[5],sl,sr,dir,dt,len,fr,wave);
	}
	else
		(void) delete_start_of_bond(btr->next,newc[6]->curve);

	interpolate_intfc_states(intfc) = sav_intrp;
	debug_print("precursor","Left new_transmission_node()\n");
	return tnode;
}		/*end new_transmission_node*/

/*
*			install_new_transmission_node():
*
*	This function installs the new transmission node by splitting the
*	old ahead contact, setting states and RP, and moving the new
*	incident to the split point (transmission node position).
*	Note that newc[6] becomes the behind contact at the new
*	transmission node, and has a new component number created for the
*	new region created by the bifurcation.
*/

LOCAL	NODE *install_new_transmission_node(
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE		*oldn,
	POINT		*p,
	BOND		*btr,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	RP_DATA		*tRP)
{
	COMPONENT	left0, right0, left1, right1;
	CURVE		**curves;
	INTERFACE	*intfc = newc[0]->curve->interface;
	NODE		*tnode;
	POINT		*oldp, *newp;
	Locstate	stl, str;
	double		V[MAXD];
	boolean		sav_intrp;
	boolean		sav_scss = interpolate_states_at_split_curve_node();

	sav_intrp = interpolate_intfc_states(intfc);
	if (newc[6]->orient == POSITIVE_ORIENTATION)
	{
		left1 = negative_component(newc[6]->curve);
		right1 = positive_component(newc[6]->curve);
		if (tRP->ang_dir == CLOCKWISE)
		{
			str = tRP->state[0];
			stl = tRP->state[4];
			left0 = new_component(NEW_COMP);
			right0 = (newc[5]->orient == POSITIVE_ORIENTATION) ?
				positive_component(newc[5]->curve) :
				negative_component(newc[5]->curve);

		}
		else
		{
			stl = tRP->state[0];
			str = tRP->state[4];
			right0 = new_component(NEW_COMP);
			left0 = (newc[5]->orient == POSITIVE_ORIENTATION) ?
				negative_component(newc[5]->curve) :
				positive_component(newc[5]->curve);
		}
	}
	else
	{
		left0 = negative_component(newc[6]->curve);
		right0 = positive_component(newc[6]->curve);
		if (tRP->ang_dir == CLOCKWISE)
		{
			stl = tRP->state[0];
			str = tRP->state[4];
			right1 = new_component(NEW_COMP);
			left1 = (newc[5]->orient == POSITIVE_ORIENTATION) ?
				positive_component(newc[5]->curve) :
				negative_component(newc[5]->curve);
		}
		else
		{
			str = tRP->state[0];
			stl = tRP->state[4];
			left1 = new_component(NEW_COMP);
			right1 = (newc[5]->orient == POSITIVE_ORIENTATION) ?
				negative_component(newc[5]->curve) :
				positive_component(newc[5]->curve);
		}
	}
	set_copy_intfc_states(YES);
	interpolate_intfc_states(intfc) = NO;
	set_interpolate_states_at_split_curve_node(NO);
	curves = split_curve(p,btr,newc[6]->curve,left0,right0,
			left1,right1);
	set_interpolate_states_at_split_curve_node(sav_scss);
	if (newc[6]->orient == POSITIVE_ORIENTATION)
	{
		newc[6]->curve = curves[0];
		end_status(curves[0]) = SLIP;
		start_status(curves[1]) = CONTACT_TARGET;
		ft_assign(left_start_state(curves[1]),stl,fr->sizest);
		ft_assign(right_start_state(curves[1]),str,fr->sizest);
	}
	else
	{
		newc[6]->curve = curves[1];
		start_status(curves[1]) = SLIP;
		end_status(curves[0]) = CONTACT_TARGET;
		ft_assign(left_end_state(curves[0]),stl,fr->sizest);
		ft_assign(right_end_state(curves[0]),str,fr->sizest);
	}
	newb[6] = NULL;
	delete_interior_points_of_curve(fr,newc[6]->curve);
	tnode = curves[1]->start;
	node_type(tnode) = TRANSMISSION_NODE;
	Rp_data(tnode) = tRP;
	change_node_of_curve(newc[5]->curve,newc[5]->orient,tnode);
	set_status_at_node(newc[5]->curve,newc[5]->orient,INCIDENT);
	if (oldc[5]->curve != NULL)
	{
		newp = Point(NULL);
		oldp = oldn->posn;
		point_propagate(fr,(POINTER)wave,oldp,
			newp,Bond_at_node_of_o_curve(oldc[5]),
			oldc[5]->curve,dt,V);
		insert_point_adjacent_to_node(newp,newc[5]->curve,
			newc[5]->orient);
	}
	interpolate_intfc_states(intfc) = sav_intrp;
	return tnode;
}		/*end install_new_transmission_node*/

/*
*			transmit_precursor():
*
*	This function is basically a driver for the creation of the cross
*	node and related structures.  The cross position is found
*	geometrically. Then the cross node RP is created and states are 
*	found (if possible).  Then the total internal reflection is set
*	up.  Finally, the cross node is installed.  Various flags are
*	returned depending on how many of these steps are completed.
*/

LOCAL	TRANSMIT_PRECURSOR transmit_precursor(
	Front		*fr,
	Wave		*wave,
	NODE		*newn,
	double		dt,
	NODE_FLAG	flag,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	BOND		**newb,
	RP_DATA		*RP,
	RP_DATA		*tRP,
	POINT		*ptr)
{
	BOND		*bcr;		/* first bond on old inc at cnode */
	POINT		*pcr;		/* location of new cross node */
	RP_DATA		*cRP;		/* RP for new cross node */
	double		node_v[MAXD];	/* velocity of cross node */
	double		t1[MAXD];	/* tangent to old incident */
	double		t2[MAXD];	/* tangent to trans at tnode */
	boolean		is_refl_rarefaction;
	int		dim = fr->interf->dim;

	pcr = Point(Coords(ptr));
	if (!find_new_cross_node_position(fr,wave,dt,flag,oldc,newc,
				             tRP,t1,&bcr,ptr,&pcr))
		return CROSS_POSITION_NOT_FOUND;

	cRP = allocate_RP_DATA_structure(fr->sizest,YES,GAS_STATE);
	cRP->ang_dir = RP->ang_dir;
	t2[0] = -cos(tRP->ang[4]);	t2[1] = -sin(tRP->ang[4]);
	ft_assign(cRP->state[0],RP->state[0],fr->sizest);
	ft_assign(cRP->state[1],RP->state[1],fr->sizest);
	ft_assign(cRP->state[4],tRP->state[3],fr->sizest);
	if (!compute_node_velocity(wave_type(newc[0]->curve),
			wave_type(newc[0]->curve),cRP->state[0],
			cRP->state[1],cRP->state[4],
			t1,t2,node_v,dim,CROSS_NODE,cRP->ang_dir)
	    			 ||
	    !find_cross_node_states(node_v,cRP,&is_refl_rarefaction,
			       (RP->ang_dir == COUNTER_CLOCK) ? YES : NO))
		return CROSS_POSITION_BUT_NO_STATES;

	if (!reset_shock_diffraction(newn,newc,pcr,ptr,RP,tRP,cRP,fr))
		return SHOCK_DIFFRACTION_NOT_RESET;

	if (!new_cross_node(pcr,bcr,newc,newb,RP,cRP,node_v,fr))
		return CROSS_NODE_NOT_INSTALLED;

	return PRECURSOR_TRANSMITTED;
}		/*end transmit_precursor*/


/*
*			find_new_cross_node_position():
*
*	This function finds the location of the new cross node by intersecting
*	a ray at the angle of the transmitted wave from the transmission
*	node with the propagated old incident shock.  A point is inserted
*	in the old incident at this location.
*/

LOCAL	int find_new_cross_node_position(
	Front		*fr,
	Wave		*wave,
	double		dt,
	NODE_FLAG	flag,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	RP_DATA		*tRP,
	double		*t1,
	BOND		**bcr,
	POINT		*ptr,
	POINT		**pcr)
{
	BOND		*oppb0;			/* opp bond on old incident */
	BOND		B0, *b0virtual = &B0;	/* virtual propagated bond
						   at start of old incident */
	INTERFACE	*intfc = newc[0]->curve->interface;
	NODE		*oppn0;			/* opp node of old incident */
	double		v0[MAXD];	/* tangent to trans wave at tnode */
	int		dim = fr->interf->dim;
	boolean		sav_intrp;
	static POINT	*p0 = NULL, *p0_opp = NULL;

	debug_print("precursor","Entered find_new_cross_node_position()\n");
	if (p0 == NULL)
	{
	    p0 = Static_point(fr->interf);
	    p0_opp = Static_point(fr->interf);
	}

	init_curve_for_crossing(p0,p0_opp,b0virtual,oldc[0],newc[0],
		&oppn0,&oppb0,fr,(POINTER)wave,dt,v0,flag);
	v0[0] = cos(tRP->ang[4]);	v0[1] = sin(tRP->ang[4]);
	if (debugging("precursor"))
	{
		(void) printf("ptr = <%g, %g>, v0 = <%g, %g>\n",Coords(ptr)[0],
			Coords(ptr)[1],v0[0],v0[1]);
		(void) printf("b0virtual, ");	print_bond(b0virtual);
		(void) printf("oppb0, ");	print_bond(oppb0);
		(void) printf("newc[0], ");	print_o_curve(newc[0]);
	}
	if (!intersect_ray_with_curve(ptr,v0,b0virtual,oppb0,
			newc[0]->curve,newc[0]->orient,bcr,*pcr))
	{
		(void) printf("WARNING in find_new_cross_node_position(), ");
		(void) printf("new cross node position not found\n");
		debug_print("precursor","Left find_new_cross_node_position()\n");
		return NO;
	}
	if (debugging("precursor"))
	{
		(void) printf("New cross node at <%g, %g>\n",
			Coords(*pcr)[0],Coords(*pcr)[1]);
		(void) printf("Bcr\n");
		print_bond(*bcr);
	}
	if (*bcr == b0virtual) *bcr = Bond_at_node_of_o_curve(newc[0]);


	sav_intrp = interpolate_intfc_states(intfc);
	interpolate_intfc_states(intfc) = NO;
	if (insert_point_in_bond(*pcr,*bcr,newc[0]->curve) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in reset_shock_diffraction(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	while (*pcr != Point_adjacent_to_node(newc[0]->curve,newc[0]->orient))
	    (void) delete_point_adjacent_to_node(fr,newc[0]->curve,
						 newc[0]->orient);
	interpolate_intfc_states(intfc) = sav_intrp;
	*bcr = Bond_at_node_of_o_curve(newc[0]);
	find_tangent_to_curve(*pcr,*bcr,newc[0]->curve,newc[0]->orient,t1,fr);
	if (newc[0]->orient == POSITIVE_ORIENTATION)
	    oppb0->end = oppn0->posn;
	else
	    oppb0->start = oppn0->posn;
	set_bond_length(oppb0,dim);
	debug_print("precursor","Left find_new_cross_node_position()\n");
	return YES;
}		/*end find_new_cross_node_position*/


/*
*			reset_shock_diffraction():
*
*	This function computes the location and states at the total internal
*	reflection node.  The newn corresponding to the old diffraction node
*	becomes the new total internal reflection node.  After this function
*	the RP structure will reflect this, but not newc.  The old incident
*	has not yet been split at the cross node, so the incident at the
*	total internal reflection doesn't really exist at this point.
*/

LOCAL	int reset_shock_diffraction(
	NODE		*newn,
	O_CURVE		**newc,
	POINT		*pcr,
	POINT		*ptr,
	RP_DATA		*RP,
	RP_DATA		*tRP,
	RP_DATA		*cRP,
	Front		*fr)
{
	Locstate	st;
	double		t0[MAXD], *t1;
	double		theta;
	int		dim = fr->interf->dim;
	int		i;

	debug_print("precursor","Entered reset_shock_diffraction()\n");
	uni_array(&t1,dim,FLOAT);
	t0[0] = cos(cRP->ang[3]);	t0[1] = sin(cRP->ang[3]);
	t1[0] = cos(tRP->ang[3]);	t1[1] = sin(tRP->ang[3]);
	if (!intersect_ray_with_sector(pcr,ptr,t0,&t1,
						Coords(newn->posn),dim))
	{
		free(t1);
		(void) printf("WARNING in reset_shock_diffraction(), ");
		(void) printf("no intersection of rays\n");
		return NO;
	}
	for (i = 0; i < dim; i++)
	{
		t0[i] = -t0[i];	t1[i] = -t1[i];
	}
	if (debugging("precursor"))
	{
		(void) printf("Cross node position <%g, %g>\n",
			Coords(pcr)[0],Coords(pcr)[1]);
		(void) printf("t0 = <%g, %g>\n",t0[0],t0[1]);
		(void) printf("Transmission node position <%g, %g>\n",
			Coords(ptr)[0],Coords(ptr)[1]);
		(void) printf("t1 = <%g, %g>\n",t1[0],t1[1]);
		(void) printf("Total internal reflection node position <%g, %g>\n",
			Coords(newn->posn)[0],Coords(newn->posn)[1]);
	}
	RP->ang[0] = angle(t0[0],t0[1]);
	RP->ang[6] = angle(t1[0],t1[1]);
	ft_assign(RP->state[0],tRP->state[3],fr->sizest);
	ft_assign(RP->state[1],cRP->state[3],fr->sizest);
	ft_assign(RP->state[6],tRP->state[1],fr->sizest);
	if (!compute_node_velocity(wave_type(newc[0]->curve),
			wave_type(newc[6]->curve),RP->state[0],
			RP->state[1],RP->state[6],
			t0,t1,Node_vel(newn),dim,node_type(newn),RP->ang_dir))
	{
		free(t1);
		(void) printf("WARNING in reset_shock_diffraction(), ");
		(void) printf("can't compute node velocity\n");
		return NO;
	}
	if (debugging("precursor"))
	{
		(void) printf("Estimated node velocity = <%g, %g>\n",
			Node_vel(newn)[0],Node_vel(newn)[1]);
	}
	set_state(RP->state[6],TGAS_STATE,RP->state[6]);
	for (i = 0; i < dim; i++)
		Vel(RP->state[6])[i] = Node_vel(newn)[i];
	set_state(RP->state[6],GAS_STATE,RP->state[6]);
	ft_assign(RP->state[5],RP->state[6],fr->sizest);
	if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[6]->orient))
		st = Left_state_at_node_of_o_curve(newc[6]);
	else
		st = Right_state_at_node_of_o_curve(newc[6]);

	ft_assign(st,RP->state[6],fr->sizest);
	if (!prandtl_meyer_wave(RP->state[1],
			pressure(RP->state[6]),
			(RP->ang_dir == COUNTER_CLOCK) ? NO : YES,
			Node_vel(newn),RP->state[4],
			&RP->ang[1],&RP->ang[3],&theta))
	{
		free(t1);
		(void) printf("WARNING in reset_shock_diffraction(), ");
		(void) printf("prandtl_meyer_wave failed\n");
		return NO;
	}
	RP->ang[2] = RP->ang[1];
	ft_assign(RP->state[2],RP->state[1],fr->sizest);
	ft_assign(RP->state[3],RP->state[4],fr->sizest);
	RP->ang[4] = angle(vel(0,RP->state[4]) - Node_vel(newn)[0],
			vel(1,RP->state[4]) - Node_vel(newn)[1]);
	RP->ang[5] = RP->ang[4];

	node_type(newn) = TOT_INT_REFL_NODE;
	free(t1);
	if (debugging("precursor"))
	{
		print_RP_node_states(
			"States after reset_shock_diffraction()",
			Node_vel(newn),RP,node_type(newn));
	}
	debug_print("precursor","Left reset_shock_diffraction()\n");
	return YES;
}		/*end reset_shock_diffraction*/


/*
*			new_cross_node():
*
*	This function does the actual work for the changes in topology
*	needed for the insertion of the cross node.  The old incident
*	is split at the cross position, and the start/end curve states
*	are set.  Then a new curve is created between the cross and
*	transmission nodes.  newc[0] is also reset to be the incident
*	at the total internal reflection, completing the transfer from
*	diffraction node to total internal reflection.
*/

LOCAL int new_cross_node(
	POINT		*pcr,
	BOND		*bcr,
	O_CURVE		**newc,
	BOND		**newb,
	RP_DATA		*RP,
	RP_DATA		*cRP,
	double		*node_v,
	Front		*fr)
{
	COMPONENT	left0, right0, left1, right1;
	COMPONENT	newcomp, aheadcomp;
	CURVE		*c, **curves;
	INTERFACE	*intfc = newc[0]->curve->interface;
	Locstate	lst0, rst0, lst1, rst1;
	NODE		*cnode;
	NODE		*ns, *ne;
	int		w_type;
	int		sstatus, estatus;
	int		i, dim = fr->rect_grid->dim;
	boolean		sav_intrp;
	boolean		sav_scss = interpolate_states_at_split_curve_node();

	if (newc[5]->curve == NULL) return YES;

	debug_print("precursor","Entered new_cross_node()\n");
	sav_intrp = interpolate_intfc_states(intfc);
	newcomp = (curve_ang_oriented_l_to_r(RP->ang_dir,newc[6]->orient)) ?
			positive_component(newc[6]->curve) :
			negative_component(newc[6]->curve);
	if (newc[0]->orient == POSITIVE_ORIENTATION)
	{
		left1 = negative_component(newc[0]->curve);
		right1 = positive_component(newc[0]->curve);
		if (RP->ang_dir == COUNTER_CLOCK)
		{
			aheadcomp = right1;
			SetActiveFlowComponent(left1,fr);
			left0 = negative_component(newc[0]->curve);
			right0 = newcomp;
			lst0 = cRP->state[3];
			rst0 = cRP->state[4];
			lst1 = cRP->state[1];
			rst1 = cRP->state[0];
		}
		else
		{
			aheadcomp = left1;
			SetActiveFlowComponent(right1,fr);
			left0 = newcomp;
			right0 = positive_component(newc[0]->curve);
			lst0 = cRP->state[4];
			rst0 = cRP->state[3];
			lst1 = cRP->state[0];
			rst1 = cRP->state[1];
		}
	}
	else
	{
		left0 = negative_component(newc[0]->curve);
		right0 = positive_component(newc[0]->curve);
		if (RP->ang_dir == CLOCKWISE)
		{
			aheadcomp = right0;
			SetActiveFlowComponent(left0,fr);
			left1 = negative_component(newc[0]->curve);
			right1 = newcomp;
			lst0 = cRP->state[1];
			rst0 = cRP->state[0];
			lst1 = cRP->state[3];
			rst1 = cRP->state[4];
		}
		else
		{
			aheadcomp = left0;
			SetActiveFlowComponent(right0,fr);
			left1 = newcomp;
			right1 = positive_component(newc[0]->curve);
			lst0 = cRP->state[0];
			rst0 = cRP->state[1];
			lst1 = cRP->state[4];
			rst1 = cRP->state[3];
		}
	}
	set_copy_intfc_states(YES);
	interpolate_intfc_states(intfc) = NO;
	set_interpolate_states_at_split_curve_node(NO);
	curves = split_curve(pcr,bcr,newc[0]->curve,left0,right0,
			left1,right1);
	set_interpolate_states_at_split_curve_node(sav_scss);
	if (newc[0]->orient == POSITIVE_ORIENTATION)
	{
		start_status(curves[1]) = INCIDENT;
		end_status(curves[0]) = REFLECTED;
		newc[0]->curve = curves[0];
		newb[0] = newc[0]->curve->first;
	}
	else
	{
		start_status(curves[1]) = REFLECTED;
		end_status(curves[0]) = INCIDENT;
		newc[0]->curve = curves[1];
		newb[0] = newc[0]->curve->last;
	}
	ft_assign(left_end_state(curves[0]),lst0,fr->sizest);
	ft_assign(right_end_state(curves[0]),rst0,fr->sizest);
	ft_assign(left_start_state(curves[1]),lst1,fr->sizest);
	ft_assign(right_start_state(curves[1]),rst1,fr->sizest);

	cnode = curves[0]->end;
	node_type(cnode) = CROSS_NODE;
	Rp_data(cnode) = cRP;
	for (i = 0; i < dim; i++) Node_vel(cnode)[i] = node_v[i];

	if (newc[5]->orient == POSITIVE_ORIENTATION)
	{
		/* CROSS_NODE           TRANSMISSION_NODE */
		ns = cnode;		ne = Opp_node_of_o_curve(newc[6]);
		sstatus = INCIDENT;	estatus = TRANSMITTED;
	}
	else
	{
		/* CROSS_NODE           TRANSMISSION_NODE */
		ne = cnode;		ns = Opp_node_of_o_curve(newc[6]);
		estatus = INCIDENT;	sstatus = TRANSMITTED;
	}
	if (curve_ang_oriented_l_to_r(cRP->ang_dir,newc[5]->orient))
	{
		left1 = newcomp;	right1 = aheadcomp;
		lst1 = cRP->state[4];	rst1 = cRP->state[0];
		w_type = FORWARD_SHOCK_WAVE;
	}
	else
	{
		right1 = newcomp;	left1 = aheadcomp;
		rst1 = cRP->state[4];	lst1 = cRP->state[0];
		w_type = BACKWARD_SHOCK_WAVE;
	}
	c = make_curve(left1,right1,ns,ne);
	start_status(c) = sstatus;	end_status(c) = estatus;
	wave_type(c) = w_type;
	ft_assign(left_start_state(c),lst1,fr->sizest);
	ft_assign(left_end_state(c),lst1,fr->sizest);
	ft_assign(right_start_state(c),rst1,fr->sizest);
	ft_assign(right_end_state(c),rst1,fr->sizest);
	interpolate_intfc_states(intfc) = sav_intrp;
	if (debugging("precursor"))
	{
		(void) printf("New cross node at <%g, %g>\n",
			Coords(pcr)[0],Coords(pcr)[1]);
	}
	debug_print("precursor","Left new_cross_node()\n");
	return YES;
}		/*end new_cross_node*/

/*
*			precursor_shock_rr_reconfigure():
*
*	This function is called when the cross-node position has been found,
*	but the states cannot be computed (for example if the angle between
*	the inc shock and the transmitted precursor is too large).  The
*	transmission node has already been installed, and states around it
*	set except for the (new) transmitted shock.  The point pcr has
*	already been inserted in find_new_cross_node_position().  Then the end
*	of the inc shock is redirected.  It now runs from pcr to the
*	transmission node (instead of to the diff node).  Then we need to
*	set the states for this curve at the transmission node, and
*	interpolate states for pcr using the adjacent curve states.  Then
*	any reflected waves are deleted (untracked) and components reset if
*	necessary.  The old front and back contacts are then joined.  During
*	this process the old diff node is propagated to a new position, and
*	becomes a regular point.  Finally, the angle for the joined curve
*	must be adjusted due to the deletion of the old diff node.
*	Notice that the inc shock for the diff node is now the transmitted
*	shock for the transmission node.
*/

LOCAL int precursor_shock_rr_reconfigure(
	Front		*fr,
	Wave		*wave,
	RPROBLEM	**rp,
	O_CURVE		**newc,
	O_CURVE		**oldc,
	NODE		*tnode,
	RP_DATA		*tRP,
	ANGLE_DIRECTION	ang_dir,
	double		dt)
{
	COMPONENT	comp1;		/*comp in state 1 of diff node*/
	BOND		*b0, *b1;
	CURVE		*inc_shock = newc[0]->curve;
	CURVE		*front_cont = newc[6]->curve;
					/*between tnode & old diff node*/
	NODE	        *oppn;
	POINT		*p0, *p1, *pcr;
	UNTRACK_FLAG	untrack_flag;
	ORIENTATION	inc_or = newc[0]->orient;
	ORIENTATION	fcont_or = newc[6]->orient;
	boolean         oppn_ss;
	int		ref_wave = NO;
	double		t0, t1, alpha, beta;
	double		len, dir[MAXD];

	debug_print("precursor","Entered precursor_shock_rr_reconfigure()\n");
	if (debugging("precursor"))
	{
		(void) printf("WARNING in precursor_shock_rr_propagate(), ");
		(void) printf("CROSS_NODE location found, but\n");
		(void) printf("\tunable to compute states.  Reconfiguring.\n");
	}
	change_node_of_curve(inc_shock,inc_or,tnode);

	if (inc_or == POSITIVE_ORIENTATION)
		start_status(inc_shock) = TRANSMITTED;
	else
		end_status(inc_shock) = TRANSMITTED;

	if (curve_ang_oriented_l_to_r(ang_dir,inc_or))
	{
		ft_assign(Left_state_at_node(inc_shock,inc_or),
				tRP->state[4],fr->sizest);
		ft_assign(Right_state_at_node(inc_shock,inc_or),
				tRP->state[3],fr->sizest);
	}
	else
	{
		ft_assign(Left_state_at_node(inc_shock,inc_or),
				tRP->state[3],fr->sizest);
		ft_assign(Right_state_at_node(inc_shock,inc_or),
				tRP->state[4],fr->sizest);
	}

	/*interpolate states for pcr*/
	pcr = Point_adjacent_to_node(inc_shock,inc_or);
	b0 = Bond_at_node(inc_shock,inc_or);
	b1 = Following_bond(b0,inc_or);
	p0 = tnode->posn;
	p1 = Point_of_bond(b1,Opposite_orient(inc_or));
	t0 = bond_length(b0);			t1 = bond_length(b1);
	alpha = t1/(t1+t0);			beta = t0/(t0+t1);
	
	interpolate_states(fr,alpha,beta,Coords(p0),
			Right_state_at_node(inc_shock,inc_or),Coords(p1),
			right_state_at_point_on_curve(p1,b1,inc_shock),
			right_state_at_point_on_curve(pcr,b0,inc_shock));
	interpolate_states(fr,alpha,beta,Coords(p0),
			Left_state_at_node(inc_shock,inc_or),Coords(p1),
			left_state_at_point_on_curve(p1,b1,inc_shock),
			left_state_at_point_on_curve(pcr,b0,inc_shock));

	comp1 = (curve_ang_oriented_l_to_r(ang_dir,inc_or)) ?
			positive_component(inc_shock) :
			negative_component(inc_shock);

	/*comps must be set consistently for untrack_curve()*/
	if (curve_ang_oriented_l_to_r(ang_dir,fcont_or))
		positive_component(front_cont) = comp1;
	else
		negative_component(front_cont) = comp1;

	/*set states using back contact -- this maintains consistency
		between ref_wave and no ref_wave cases*/
	init_redundant_node_for_deletion(Node_of_o_curve(newc[4]),
		                         Node_of_o_curve(oldc[4]),fr,
					 (POINTER)wave,dt);

	if (newc[1]->curve != NULL)	/*leading edge of rarefaction*/
	{
	    oppn = Opp_node_of_o_curve(newc[1]);
	    oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	    ref_wave = YES;
	    set_untrack_flag(untrack_flag,newc[1]->orient,NO,oppn_ss,NO,YES,NO);
	    (void) untrack_curve(newc[1],oldc[1],comp1,dt,fr,(POINTER)wave,
			         *rp,untrack_flag);
	    newc[1]->curve = NULL;
	}
	if (newc[3]->curve != NULL)	/*trailing edge of rarefaction*/
	{
	    oppn = Opp_node_of_o_curve(newc[3]);
	    oppn_ss = (propagation_status(oppn)==PROPAGATED_NODE) ? YES : NO;
	    ref_wave = YES;
	    set_untrack_flag(untrack_flag,newc[3]->orient,NO,oppn_ss,NO,YES,NO);
	    (void) untrack_curve(newc[3],oldc[3],comp1,dt,fr,(POINTER)wave,
			         *rp,untrack_flag);
	    newc[3]->curve = NULL;
	}

	if (!ref_wave)
	{
		(void) delete_redundant_node(Node_of_o_curve(newc[4]),
					     (CROSS *)NULL,*rp,fr);
	}

	/*adjust angle of back cont -- newc[6] used as temp storage*/
	find_curve_with_status(tnode,&newc[6]->curve,&newc[6]->orient,SLIP);
	(void) delete_point_adjacent_to_node(fr,newc[6]->curve,
					     newc[6]->orient);
	len = bond_length(Bond_at_node_of_o_curve(newc[6]));
	dir[0] = cos(tRP->ang[3]);	dir[1] = sin(tRP->ang[3]);
	if (!adjust_angle_at_node(tnode,(O_CURVE *)NULL,newc[6],
				     tRP->state[4],tRP->state[0],
				     dir,dt,len,fr,wave))
	{
		(void) printf("WARNING in precursor_shock_rr_reconfigure():");
		(void) printf("unable to adjust_angle_at_node for back contact.\n");
		return ERROR_NODE;
	}

	if (debugging("precursor"))
		(void) printf("Reconfigure successful.\n");
	debug_print("precursor","Left precursor_shock_rr_reconfigure()\n");
	return GOOD_NODE;
}		/*end precursor_shock_rr_reconfigure*/

/*
*			final_shock_at_precursor_rr():
*
*	This function installs the tracking for the reflected shock at the
*	cross node which is directed toward the reflected rarefaction.  This
*	creates an overtake node at the point of intersection of the reflected
*	shock and the leading edge of the rarefaction, with creation of a
*	new component.
*/

LOCAL	void final_shock_at_precursor_rr(
	NODE		*newn,
	O_CURVE		**oldc,
	O_CURVE		**newc,
	Front		*fr)
{
	BOND		*bint;			/* bond for onode */
	COMPONENT	left_end_comp0;		/* at split, curve 0 */
	COMPONENT	right_end_comp0;
	COMPONENT	left_start_comp1;	/* at split, curve 1 */
	COMPONENT	right_start_comp1;
	COMPONENT	newcomp;
	COMPONENT	final_shock_lcomp;
	COMPONENT	final_shock_rcomp;
	CURVE		**curves;		/* from split */
	CURVE		*final_shock;
	INTERFACE	*intfc = newn->interface;
	Locstate	left_end_st0;		/* at split, curve 0 */
	Locstate	right_end_st0;
	Locstate	left_start_st1;		/* at split, curve 1 */
	Locstate	right_start_st1;
	NODE		*ns, *ne;
	POINT		*pc;			/* location of onode */
	RP_DATA		*cRP, *RP;
	ORIENTATION	final_shock_or;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	boolean		sav_scss = interpolate_states_at_split_curve_node();

	debug_print("precursor","Entered final_shock_at_precursor_rr()\n");

	if (newc[1] == NULL || newc[1]->curve == NULL ||
		oldc[1] == NULL || oldc[1]->curve == NULL)
	{
		if (debugging("precursor"))
			(void) printf("Rarefaction leading edge untracked\n");
		debug_print("precursor","Left final_shock_at_precursor_rr()\n");
		return;
	}

	if (debugging("precursor"))
	{
		(void) printf("Rarefaction leading edge\n");
		print_o_curve(newc[1]);
	}

	cRP = Rp_data(Opp_node_of_o_curve(newc[0]));
	RP = Rp_data(newn);

	newcomp = new_component(NEW_COMP);
	if (newc[1]->orient == POSITIVE_ORIENTATION)
	{
		left_start_comp1 = negative_component(newc[1]->curve);
		right_start_comp1 = positive_component(newc[1]->curve);
		left_end_st0 = cRP->state[2];
		right_end_st0 = cRP->state[2];
		if (RP->ang_dir == COUNTER_CLOCK)
		{
			left_end_comp0 = negative_component(newc[1]->curve);
			right_end_comp0 = newcomp;

			final_shock_lcomp = newcomp;
			final_shock_rcomp = right_start_comp1;

			left_start_st1 = cRP->state[1];
			right_start_st1 = cRP->state[2];
		}
		else
		{
			left_end_comp0 = newcomp;
			right_end_comp0 = positive_component(newc[1]->curve);

			final_shock_lcomp = left_start_comp1;
			final_shock_rcomp = newcomp;

			left_start_st1 = cRP->state[2];
			right_start_st1 = cRP->state[1];
		}
	}
	else
	{
		left_end_comp0 = negative_component(newc[1]->curve);
		right_end_comp0 = positive_component(newc[1]->curve);
		left_start_st1 = cRP->state[2];
		right_start_st1 = cRP->state[2];
		if (RP->ang_dir == CLOCKWISE)
		{
			left_start_comp1 = negative_component(newc[1]->curve);
			right_start_comp1 = newcomp;

			final_shock_lcomp = newcomp;
			final_shock_rcomp = right_end_comp0;

			left_end_st0 = cRP->state[1];
			right_end_st0 = cRP->state[2];
		}
		else
		{
			left_start_comp1 = newcomp;
			right_start_comp1 = positive_component(newc[1]->curve);

			final_shock_lcomp = left_end_comp0;
			final_shock_rcomp = newcomp;

			left_end_st0 = cRP->state[2];
			right_end_st0 = cRP->state[1];
		}
	}
	pc = Point(Coords(newn->posn));
	set_copy_intfc_states(YES);
	interpolate_intfc_states(intfc) = NO;
	if (!find_overtake_node_position(&bint,newc[1]->curve,
		    Opp_node_of_o_curve(newc[0]),pc,newc[1]->orient,fr))
	{
	    if (debugging("precursor"))
	    {
	    	(void) printf("WARNING in final_shock_at_precursor_rr(), "
	    	              "unable to find overtake node position\n");
	    }
	    debug_print("precursor","Leaving final_shock_at_precursor_rr()\n");
	    return;
	}
	if (insert_point_in_bond(pc,bint,newc[1]->curve) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in final_shock_at_precursor_rr(), "
		   "insert_point_in_bond() failed\n");
	    clean_up(ERROR);
	}
	set_interpolate_states_at_split_curve_node(NO);
	curves = split_curve(pc,bint,newc[1]->curve,left_end_comp0,
		     right_end_comp0,left_start_comp1,right_start_comp1);
	set_interpolate_states_at_split_curve_node(sav_scss);
	ft_assign(left_end_state(curves[0]),left_end_st0,fr->sizest);
	ft_assign(right_end_state(curves[0]),right_end_st0,fr->sizest);
	ft_assign(left_start_state(curves[1]),left_start_st1,fr->sizest);
	ft_assign(right_start_state(curves[1]),right_start_st1,fr->sizest);
	node_type(curves[0]->end) = OVERTAKE_NODE;
	if (newc[1]->orient == POSITIVE_ORIENTATION)
	{
		delete_interior_points_of_curve(fr,curves[0]);
		end_status(curves[0]) = INCIDENT;
		start_status(curves[1]) = TRANSMITTED;
		newc[1]->curve = curves[0];
		ns = Opp_node_of_o_curve(newc[0]);
		ne = curves[0]->end;
	}
	else
	{
		delete_interior_points_of_curve(fr,curves[1]);
		end_status(curves[0]) = TRANSMITTED;
		start_status(curves[1]) = INCIDENT;
		newc[1]->curve = curves[1];
		ns = curves[0]->end;
		ne = Opp_node_of_o_curve(newc[0]);
	}
	if (debugging("precursor"))
	{
		(void) printf("curves from split\n");
		(void) printf("curves[0]");	print_curve(curves[0]);
		(void) printf("curves[1]");	print_curve(curves[1]);
	}
	final_shock = make_curve(final_shock_lcomp,final_shock_rcomp,ns,ne);
	final_shock_or = newc[1]->orient;
	if (curve_ang_oriented_l_to_r(cRP->ang_dir,final_shock_or))
	{
		wave_type(final_shock) = BACKWARD_SHOCK_WAVE;
		ft_assign(left_start_state(final_shock),cRP->state[1],fr->sizest);
		ft_assign(right_start_state(final_shock),cRP->state[2],fr->sizest);
		ft_assign(left_end_state(final_shock),cRP->state[1],fr->sizest);
		ft_assign(right_end_state(final_shock),cRP->state[2],fr->sizest);
	}
	else
	{
		wave_type(final_shock) = FORWARD_SHOCK_WAVE;
		ft_assign(left_start_state(final_shock),cRP->state[2],fr->sizest);
		ft_assign(right_start_state(final_shock),cRP->state[1],fr->sizest);
		ft_assign(left_end_state(final_shock),cRP->state[2],fr->sizest);
		ft_assign(right_end_state(final_shock),cRP->state[1],fr->sizest);
	}
	if (newc[1]->orient == POSITIVE_ORIENTATION)
	{
		start_status(final_shock) = REFLECTED;
		end_status(final_shock) = OVERTOOK;
	}
	else
	{
		start_status(final_shock) = OVERTOOK;
		end_status(final_shock) = REFLECTED;
	}
	if (curve_ang_oriented_l_to_r(RP->ang_dir,newc[0]->orient))
		positive_component(newc[0]->curve) = newcomp;
	else
		negative_component(newc[0]->curve) = newcomp;

	interpolate_intfc_states(intfc) = sav_intrp;
	debug_print("precursor","Left final_shock_at_precursor_rr()\n");
}		/*end final_shock_at_precursor_rr*/


/*
*			find_overtake_node_position();
*
*	This function finds the position of the overtake node, and the
*	bond to insert it in in the leading edge of the rarefaction.  The
*	position is found by intersecting two rays.  One starts at the
*	cross node, in the direction of cRP->ang[1].  The other starts
*	at the total internal reflection in the direction of the Mach angle
*	as computed using cRP->state[2] (the state ahead of the rarefaction
*	edge at the overtake) in the frame of reference of the cross node.
*	The position is found by writing the two rays as p1 + a1*t1 and
*	p2 + a2*t2 where the p's are the starting points, the t's tangent
*	uni_arrays, and the a's fractional distances along the ray.  Equating
*	these and taking cross products yields the a's, and thus the position.
*	We then find a bond in which to insert this point by stepping along
*	the rarefaction edge and finding the first point farther from the
*	refl node than pc.
*/

LOCAL int find_overtake_node_position(
	BOND		**bint,
	CURVE		*raref_curve,
	NODE		*cnode,
	POINT		*pc,
	ORIENTATION	raref_orient,
	Front		*fr)
{
	POINT		*ref_point;
	RP_DATA		*cRP = Rp_data(cnode);
	const double	eps = MACH_EPS;/*TOLERANCE*/
	double		A, M;		/* Mach angle and Mach number */
	double		inc_ang;	/* overtaking angle + PI/2 */
	double		t1[MAXD];	/* overtaking tangent + PI/2 */
	double		t2[MAXD];	/* overtaken tangent + PI/2 */
	double		p[MAXD];
	double		den, num1, num2;
	double		a1;		/* fract dist, rnode to onode */
	int		i, dim = fr->rect_grid->dim;

	M = mach_number(cRP->state[2],Node_vel(cnode));

	if (M < SONIC_MINUS)
		return NO;	/* state ahead of mach line subsonic */
	else if (M < 1.0)
	{
		M = 1.0;
		A = 0.5*PI;
	}
	else
		A = fabs(asin(1.0 / M));

	inc_ang = (cRP->ang_dir == CLOCKWISE) ?
			cRP->ang[1] + A : cRP->ang[1] - A;

	t1[0] = cos(inc_ang);
	t1[1] = sin(inc_ang);

	t2[0] = cos(cRP->ang[1]);
	t2[1] = sin(cRP->ang[1]);

	(void) vector_product(t1,t2,&den,dim);
	if (fabs(den) < eps)
		return NO;		/* tangents are parallel */

	ref_point = Node_of(raref_curve,raref_orient)->posn;
	for (i = 0; i < dim; i++)
		p[i] = Coords(cnode->posn)[i] - Coords(ref_point)[i];
	(void) vector_product(p,t2,&num1,dim);
	(void) vector_product(p,t1,&num2,dim);

	a1 = num1 / den;		/* a2 is not needed */
	if ((a1 <= 0.0) || (num2*den <= 0.0))
		return NO;		/* intersection is invalid */

	for (i = 0; i < dim; i++)
		Coords(pc)[i] = Coords(ref_point)[i] + a1*t1[i];

	*bint = Bond_at_node(raref_curve,raref_orient);
	while (separation(Point_of_bond(*bint,Opposite_orient(raref_orient)),
		     ref_point,dim)		< a1)
	{
		if (!(*bint = Following_bond(*bint,raref_orient)))
			return NO;
	}
	return YES;
}		/*end find_overtake_node_position*/
#endif /* defined(TWOD) */
