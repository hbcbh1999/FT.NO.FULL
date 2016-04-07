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
*				gsndnode.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Contains the routine used for the second propagation of the nodes:
*
*	g_snd_node_propagate() performs the tangential sweep
*		update of states near the node.
*
*	as well as several of its support subroutines.
*	All access to this function is though function pointers
*	in the Front structure that are set by init_physics() in ginit.c.
*
*/

#if defined(TWOD)

#include <gdecs/gdecs.h>


	/* LOCAL Function Declarations */
LOCAL	boolean	comp_is_interior_on_curve(COMPONENT,CURVE*);
LOCAL	boolean	snd_B_node_propagate(Front*,Front*,Wave*,
				     NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_B_reflect_node_propagate(Front*,Front*,
					     NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_CC_node_propagate(Front*,Front*,NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_Mach_node_propagate(Front*,Front*,NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_SC_node_propagate(Front*,Front*,NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_SS_node_propagate(Front*,Front*,NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_attached_b_node_propagate(Front*,Front*,Wave*,
					      NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_contact_wall_node_propagate(Front*,Front*,Wave*,
						NODE*,NODE*,NODE*,double);
LOCAL	boolean	snd_cross_node_propagate(Front*,Front*,NODE*,NODE*,NODE*,
					 double);
LOCAL	boolean	snd_fixed_node_propagate(Front*,Front*,Wave*,
					 NODE*,NODE*,NODE*,double);
LOCAL	int	find_active_curves(NODE*,O_CURVE*);
LOCAL	void	assign_states_around_node(O_CURVE*,O_CURVE*,Locstate,size_t);
LOCAL	void	find_states_for_comp(COMPONENT,CURVE*,CURVE*,NODE*,
				    Locstate*,int,Wave*,Front*);
LOCAL	void	corner_node_propagate(NODE*,NODE*,O_CURVE*,O_CURVE*,Front*,
				      Wave*,Locstate,double);
LOCAL	void	impose_no_slip_bc(NODE*,Front*);
LOCAL	void	update_points_crossed_by_node(CURVE*,CURVE*,ORIENTATION,
					      NODE*,Front*);
LOCAL	void	update_wall_states(Front*,NODE*,NODE*,NODE*);



/*
*			g_snd_node_propagate():
*
*	Performs the tangential sweep for the states near a node.
*/

/* ARGSUSED */
EXPORT int g_snd_node_propagate(
	Front		*fr,
	Front		*newfr,
	POINTER		p2wave,
	INTERFACE	*tempintfc,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	Wave  *wave = (Wave*)p2wave;
	boolean  ok;
	NODE  *oldn;

	debug_print("node_propagate","Entering g_snd_node_propagate()\n");
	if (debugging("node_propagate"))
	{
	    (void) printf("dt = %g\n",dt);
	    (void) printf("tempn\n");
	    print_node(tempn);
	    (void) printf("newn\n");
	    print_node(newn);
	}

	if ((tempn->in_curves == NULL) && (tempn->out_curves == NULL))
	{
	    ok = YES;
	    if (debugging("node_propagate"))
	    	(void) printf("No curves incident on node\n");
	}
	else
	{

#if defined(CHECK_FOR_BAD_STATES)
	    if (debugging("bad_state") &&
		    is_bad_state_at_node("g_snd_node_propagate",tempn))
	    {
	        screen("ERROR in g_snd_node_propagate(), bad state at tempn\n");
	        print_node(tempn);
	        print_interface(tempn->interface);
	        clean_up(ERROR);
	    }
#endif /* defined(CHECK_FOR_BAD_STATES) */

		/* Return if bifurcation has occurred at node */

	    oldn = node_corresponding_to(tempn,fr);
	    if (!oldn)
	    {
	        ok = YES;
	        if (debugging("node_propagate"))
	    	    (void) printf("No corresponding old node\n");
	    }
	    else
	    {


	        /* Propagate node according to its type */

	        switch (node_type(tempn))
	        {
	        case PASSIVE_NODE:
	        case SUBDOMAIN_NODE:
	        case CLOSED_NODE:
	            ok = YES;
	            break;
        
	        case FIXED_NODE:
	            if (!is_virtual_fixed_node(tempn))
	                ok = snd_fixed_node_propagate(fr,newfr,wave,oldn,
					              tempn,newn,dt);
	            else
		        ok = YES;
	            break;

	        case NEUMANN_NODE:
	        {
	            void (*save_impose_bc)(POINT*,BOND*,CURVE*,double*,Front*,
				           boolean,boolean);

	            save_impose_bc = fr->impose_bc;
	            fr->impose_bc = f_impose_bc;
	            {
		        CURVE *cphys;
		        ORIENTATION cphys_orient;

		        cphys = find_physical_curve_at_node(tempn,&cphys_orient);
		        if (is_shock_wave(wave_type(cphys)))
		            fr->impose_bc = f_impose_bc;
	            }

	            ok = snd_B_node_propagate(fr,newfr,wave,oldn,tempn,newn,dt);
	            update_wall_states(fr,oldn,tempn,newn);
	            fr->impose_bc = save_impose_bc;
	            break;
	        }
	        case DIRICHLET_NODE:
	            ok = snd_B_node_propagate(fr,newfr,wave,oldn,tempn,newn,dt);
	            break;

	        case B_REFLECT_NODE:
	            ok = snd_B_reflect_node_propagate(fr,newfr,oldn,
			                              tempn,newn,dt);
	            break;

	        case MACH_NODE:
	            ok = snd_Mach_node_propagate(fr,newfr,oldn,tempn,newn,dt);
	            break;

	        case ATTACHED_B_NODE:
	            ok = snd_attached_b_node_propagate(fr,newfr,wave,oldn,
					               tempn,newn,dt);
	            break;

#if defined(FULL_PHYSICS)
	        case WAVE_END_NODE:
	            ok = snd_wave_end_propagate(fr,newfr,oldn,tempn,newn,dt);
	            break;

	        case CROSS_NODE:
	            ok = snd_cross_node_propagate(fr,newfr,oldn,tempn,newn,dt);
	            break;

	        case OVERTAKE_NODE:
	            ok = snd_SS_node_propagate(fr,newfr,oldn,tempn,newn,dt);
	            break;

	        case DIFFRACTION_NODE:
	        case TOT_INT_REFL_NODE:
	        case TRANSMISSION_NODE:
	            ok = snd_SC_node_propagate(fr,newfr,oldn,tempn,newn,dt);
	            break;

	        case CC_NODE:
	            ok = snd_CC_node_propagate(fr,newfr,oldn,tempn,newn,dt);
	            break;

#endif /* defined(FULL_PHYSICS) */

	        default:
	            screen("ERROR in g_snd_node_propagate(), "
			   "unknown node type %d\n",node_type(tempn));
	            clean_up(ERROR);
	            break;
	        }
	    }
	}

	if (ok == NO)
	{
	    screen("ERROR in g_snd_node_propagate(), node_propagate failed\n");
	    clean_up(ERROR);
	}

	impose_no_slip_bc(newn,fr);

#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") &&
		is_bad_state_at_node("g_snd_node_propagate",newn))
	{
	    screen("ERROR in g_snd_node_propagate(), bad state at newn\n");
	    print_node(newn);
	    print_interface(newn->interface);
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
	debug_print("node_propagate","Left g_snd_node_propagate()\n");
	return YES;
}		/*end g_snd_node_propagate*/

LOCAL	void impose_no_slip_bc(
	NODE 	*n,
	Front	*fr)
{
    	CURVE           **c;
	CURVE           *adj_c;
	ORIENTATION     orient;
	ORIENTATION     adj_c_or;
	Locstate        sl, sr;
	int             dim = n->interface->dim;
	double		alpha;
	double           nor[MAXD];

	orient = NEGATIVE_ORIENTATION;
	for (c = n->in_curves; c && *c; ++c)
	{
	    if ((wave_type(*c) == NEUMANN_BOUNDARY) && (no_slip(Hyper_surf(*c))))
	    {
		sl = Left_state_at_node(*c,orient);
		sr = Right_state_at_node(*c,orient);
		alpha = 1.0 - adherence_coeff(Hyper_surf(*c));
		normal(n->posn,Hyper_surf_element(Bond_at_node(*c,
                               orient)),Hyper_surf(*c),nor,fr);
		if (is_obstacle_state(sl))
		{
		    alpha_state_velocity(alpha,sr,dim);
		    zero_normal_velocity(sr,nor,dim);
		    adj_c = adjacent_curve(*c,orient,COUNTER_CLOCK,&adj_c_or);
		    if (adj_c &&
                        wave_type(adj_c) == NEUMANN_BOUNDARY &&
                        no_slip(Hyper_surf(adj_c)))
		    {
		        normal(Node_of(adj_c,adj_c_or)->posn,
			       Hyper_surf_element(Bond_at_node(adj_c,adj_c_or)),
			       Hyper_surf(adj_c),nor,fr);

			if (adj_c_or == POSITIVE_ORIENTATION)
			{
			    sr = Right_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sr,dim);
			    zero_normal_velocity(sr,nor,dim);
			}
			else if (adj_c_or == NEGATIVE_ORIENTATION)
			{
			    sl = Left_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sl,dim);
			    zero_normal_velocity(sl,nor,dim);
			}
		    }
		}
		if (is_obstacle_state(sr))
		{
		    alpha_state_velocity(alpha,sl,dim);
		    zero_normal_velocity(sl,nor,dim);
		    adj_c = adjacent_curve(*c,orient,CLOCKWISE,&adj_c_or);
		    if (adj_c &&
			wave_type(adj_c) == NEUMANN_BOUNDARY &&
			no_slip(Hyper_surf(adj_c)))
		    {
			normal(Node_of(adj_c,adj_c_or)->posn,
				Hyper_surf_element(Bond_at_node(adj_c,
				adj_c_or)),Hyper_surf(adj_c),nor,fr);
			if (adj_c_or == NEGATIVE_ORIENTATION)
			{
			    sr = Right_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sr,dim);
			    zero_normal_velocity(sr,nor,dim);
			}
			else if (adj_c_or == POSITIVE_ORIENTATION)
			{
			    sl = Left_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sl,dim);
			    zero_normal_velocity(sl,nor,dim);
			}
		    }
		}
	    }
	}
	orient = POSITIVE_ORIENTATION;
	for (c = n->out_curves; c && *c; ++c)
	{
	    if ((wave_type(*c) == NEUMANN_BOUNDARY) && (no_slip(Hyper_surf(*c))))
	    {
		alpha = 1.0 - adherence_coeff(Hyper_surf(*c));
		sl = Left_state_at_node(*c,orient);
		sr = Right_state_at_node(*c,orient);
		normal(n->posn,Hyper_surf_element(Bond_at_node(*c,
				orient)),Hyper_surf(*c),nor,fr);
		if (is_obstacle_state(sl))
		{
		    alpha_state_velocity(alpha,sr,dim);
		    zero_normal_velocity(sr,nor,dim);
		    adj_c = adjacent_curve(*c,orient,CLOCKWISE,&adj_c_or);
		    if (adj_c &&
		    	wave_type(adj_c) == NEUMANN_BOUNDARY &&
                        no_slip(Hyper_surf(adj_c)))
		    {
		    	normal(Node_of(adj_c,adj_c_or)->posn,
				Hyper_surf_element(Bond_at_node(adj_c,
				adj_c_or)),Hyper_surf(adj_c),nor,fr);
			if (adj_c_or == POSITIVE_ORIENTATION)
			{
			    sl = Left_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sl,dim);
			    zero_normal_velocity(sl,nor,dim);
			}
			else if (adj_c_or == NEGATIVE_ORIENTATION)
			{
			    sr = Right_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sr,dim);
			    zero_normal_velocity(sr,nor,dim);
			}
		    }
		}
		if (is_obstacle_state(sr))
		{
		    alpha_state_velocity(alpha,sl,dim);
		    zero_normal_velocity(sl,nor,dim);
		    adj_c = adjacent_curve(*c,orient,COUNTER_CLOCK,&adj_c_or);
		    if (adj_c &&
                        wave_type(adj_c) == NEUMANN_BOUNDARY &&
                        no_slip(Hyper_surf(adj_c)))
		    {
			normal(Node_of(adj_c,adj_c_or)->posn,
				Hyper_surf_element(Bond_at_node(adj_c,
				adj_c_or)),Hyper_surf(adj_c),nor,fr);
			if (adj_c_or == NEGATIVE_ORIENTATION)
			{
			    sl = Left_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sl,dim);
			    zero_normal_velocity(sl,nor,dim);
			}
			else if (adj_c_or == POSITIVE_ORIENTATION)
			{
			    sr = Right_state_at_node(adj_c,adj_c_or);
		            alpha_state_velocity(alpha,sr,dim);
			    zero_normal_velocity(sr,nor,dim);
			}
		    }
		}
	    }
	}
}		/*end impose_no_slip_bc*/


/*ARGSUSED*/
LOCAL boolean snd_fixed_node_propagate(
	Front		*fr,
	Front		*newfr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	BDRY_SIDE       bside[2];
	RECT_GRID       *gr = fr->rect_grid;
	int		i, j, n_active_curves, dim = gr->dim;
	int		idir[2];
	int		isgn;
	double		tgnt[MAXD], ds;
	double		V[MAXD];
	double		dn, *nor;
	O_CURVE		Oac[10], Tac[10];
	O_CURVE		Wall_c1, Wall_c2;
	O_CURVE		Newwall1, Newwall2;
	O_CURVE		Newcurt, Newcurn;
	O_CURVE		Tempcurt, Tempcurn;
					/* tempcurt propagated tangentially */
					/* by propagating tempcurn normally */
	BOND		*tempbt;	/* bond at temp node */
	BOND		*tempbn;
	POINT		*tnp = tempn->posn;
	Locstate	ans = NULL;
	size_t		sizest = fr->sizest;
	static Locstate ans_out = NULL;
	static Locstate ansl = NULL,ansr = NULL;
	static Locstate *sa = NULL;
	static POINT	**p = NULL;
	static WSSten	*wssten = NULL;
	static int nsts = 0;

	debug_print("fixed_node","Entering snd_fixed_node_propagate()\n");
	if (debugging("fixed_node"))
	{
	    (void) printf("oldn\n");
	    print_node(oldn);
	    (void) printf("tempn\n");
	    print_node(tempn);
	    (void) printf("newn\n");
	    print_node(newn);
	}

		/* Allocate storage */

	if (wssten == NULL)
	{
	    wssten = AllocDefaultWSSten(fr);
	    nsts = wssten->nsts;
	    alloc_state(fr->interf,&ans_out,sizest);
	    alloc_state(fr->interf,&ansl,sizest);
	    alloc_state(fr->interf,&ansr,sizest);
	    uni_array(&sa,2*nsts-1,sizeof(Locstate));
	    sa += nsts-1;
	    uni_array(&p,2*nsts-1,sizeof(POINT*));
	    p += nsts-1;
	    for (i = -nsts+1; i < nsts; ++i)
		p[i] = Static_point(fr->interf);
	}
	else
	{
	    ClearWSStenData(wssten);
	    for (i = -nsts+1; i < nsts; ++i)
		sa[i] = NULL;
	}

		/* take care of passive curves */

	assign_states_on_passive_curves_at_node(newn);

		/* find active curves */

	n_active_curves = find_active_curves(oldn,Oac);
	if (n_active_curves == 0)
	    return YES;

	/* Exit if subdomain curves.  We allow more than two active
	 * curves here to allow for a physical curve(s) passing through
	 * a subdomain corner, a configuration which is really not
	 * incorrect, and will be thrown out at the next scatter anyway.
	 */

	for (i = 0; i < n_active_curves; ++i)
	    if (is_subdomain_boundary(Hyper_surf(Oac[i].curve)))
		return YES;

	if (n_active_curves != 2)
	{
	   screen("ERROR in snd_fixed_node_propagate(), "
	          "wrong number (%d) of active curves\n",n_active_curves);
	   (void) printf("oldn\n");	print_node(oldn);
	   (void) printf("tempn\n");	print_node(tempn);
	   (void) printf("newn\n");	print_node(newn);
	   (void) printf("oldn->interface\n");
	   print_interface(oldn->interface);
	   (void) printf("tempn->interface\n");
	   print_interface(tempn->interface);
	   (void) printf("newn->interface\n");
	   print_interface(newn->interface);
	   clean_up(ERROR);
	}

		/* set temp_active_curves[] */

	Check_return(
	    find_correspond_of_oriented_curve(&Oac[0],&Tac[0],tempn,fr,
					      tempn->interface),
	    snd_fixed_node_propagate)
	Check_return(
	    find_correspond_of_oriented_curve(&Oac[1],&Tac[1],tempn,fr,
					      tempn->interface),
	    snd_fixed_node_propagate)

	if ((wave_type(Tac[0].curve) == NEUMANN_BOUNDARY ||
	     wave_type(Tac[0].curve) == MOVABLE_BODY_BOUNDARY) &&
	    (wave_type(Tac[1].curve) == NEUMANN_BOUNDARY ||
	     wave_type(Tac[1].curve) == MOVABLE_BODY_BOUNDARY))
	{
	    /* Propagate the node using a 2d scheme obtained by calculating
	     * the fluxes through the interior boundary of the triangle. */

	    copy_o_curve(&Wall_c1,&Oac[0]);
	    copy_o_curve(&Wall_c2,&Oac[1]);
	    Check_return(
		find_correspond_of_oriented_curve(&Wall_c1,&Newwall1,newn,fr,
						  newn->interface),
		snd_fixed_node_propagate)
	    Check_return(
		find_correspond_of_oriented_curve(&Wall_c2,&Newwall2,newn,fr,
						  newn->interface),
		snd_fixed_node_propagate)

	    corner_node_propagate(oldn,tempn,&Wall_c1,&Wall_c2,fr,wave,
				  ans_out,dt);

	    assign_states_around_node(&Newwall1,&Wall_c1,ans_out,sizest);
	    assign_states_around_node(&Newwall2,&Wall_c2,ans_out,sizest);
	    return YES;
	}

	/* set tempcurt to the DIRICHLET curve if there is one */
	/* set tempcurn to the other curve */

	bside[0] = rect_bdry_side_for_curve(idir,NULL,Tac[0].curve,gr);
	bside[1] = rect_bdry_side_for_curve(idir+1,NULL,Tac[1].curve,gr);
	if ((wave_type(Tac[0].curve) == DIRICHLET_BOUNDARY) &&
	    (wave_type(Tac[1].curve) != DIRICHLET_BOUNDARY))
	{
	    copy_o_curve(&Tempcurt,&Tac[0]);
	    copy_o_curve(&Tempcurn,&Tac[1]);
	}
	else if ((wave_type(Tac[0].curve) != DIRICHLET_BOUNDARY) &&
	         (wave_type(Tac[1].curve) == DIRICHLET_BOUNDARY))
	{
	    copy_o_curve(&Tempcurt,&Tac[1]);
	    copy_o_curve(&Tempcurn,&Tac[0]);
	}
	else /*Both are Dirichlet boundaries */
	{
	    if ((bside[0] != NOT_A_BDRY) && (bside[1] != NOT_A_BDRY))
	    {
	        if (idir[fr->step % 2] == 0)
	        {
	            copy_o_curve(&Tempcurt,&Tac[0]);
	            copy_o_curve(&Tempcurn,&Tac[1]);
	        }
	        else
	        {
	            copy_o_curve(&Tempcurt,&Tac[1]);
	            copy_o_curve(&Tempcurn,&Tac[0]);
	        }
	    }
	    else
	    {
	        if (fr->step % 2)
	        {
	            copy_o_curve(&Tempcurt,&Tac[0]);
	            copy_o_curve(&Tempcurn,&Tac[1]);
	        }
	        else
	        {
	            copy_o_curve(&Tempcurt,&Tac[1]);
	            copy_o_curve(&Tempcurn,&Tac[0]);
	        }
	    }
	}

	    /* set arguments for npt_w_speed */

	wssten->p = tnp;
	wssten->coords = Coords(tnp);
	wssten->lcrds[0] = wssten->rcrds[0] = wssten->coords;
	wssten->pcomp = positive_component(Tempcurn.curve);
	wssten->ncomp = negative_component(Tempcurn.curve);
	wssten->front = fr;
	wssten->wave = wave;
	wssten->pjump = 0.0;
	wssten->w_type = wave_type(Tempcurn.curve);
	wssten->dt = dt;
	wssten->hs = Hyper_surf(Tempcurn.curve);
	tempbn = Bond_at_node_of_o_curve(&Tempcurn);
	wssten->hse = Hyper_surf_element(tempbn);
	nor = wssten->nor;
	normal(tnp,Hyper_surf_element(tempbn),
	       Hyper_surf(Tempcurn.curve),nor,fr);
	wssten->dn = dn = grid_size_in_direction(nor,gr->h,dim);

	tempbt = Bond_at_node_of_o_curve(&Tempcurt);
	find_tangent_to_curve(tnp,tempbt,Tempcurt.curve,
			      Tempcurt.orient,tgnt,fr);
	ds = grid_size_in_direction(tgnt,gr->h,dim);

	    /* updated states are obtained by propagating the */
	    /* states along tempcurt together with appropriate */
	    /* states with respect to tempcurn */

	isgn = (Tempcurt.orient == POSITIVE_ORIENTATION) ? 1 : -1;
	if (comp_is_interior_on_curve(negative_component(Tempcurt.curve),
				      Tempcurn.curve) == YES)
	{
	    /* update using left tempcurt (= interior) states */

	    if (Tempcurt.orient == Tempcurn.orient)
	    {
		find_states_for_comp(negative_component(Tempcurt.curve),
				     Tempcurt.curve,Tempcurn.curve,tempn,
				     wssten->sl,nsts,wave,fr);
		for (i = 1; i < nsts; ++i)
		{
		    for (j = 0; j < dim; ++j)
			wssten->lcrds[i][j] = wssten->coords[j] - i*dn*nor[j];
		}
	        wssten->sr[0] = Left_state_at_node_of_o_curve(&Tempcurt);
		for (i = 1; i < nsts; ++i)
		{
	            sa[isgn*(i-1)] = wssten->sr[i];
		    wssten->rcrds[i] = Coords(p[isgn*(i-1)]);
		}
	        states_at_distance_along_curve(tnp,tempbt,Tempcurt.curve,
					       Tempcurt.orient,ds,nsts-1,
					       sa,NULL,NULL,NULL,NULL,p,fr);
		ans = ansr;
	    }
	    else
	    {
		find_states_for_comp(negative_component(Tempcurt.curve),
				     Tempcurt.curve,Tempcurn.curve,tempn,
				     wssten->sr,nsts,wave,fr);
		for (i = 1; i < nsts; ++i)
		{
		    for (j = 0; j < dim; ++j)
			wssten->rcrds[i][j] = wssten->coords[j] + i*dn*nor[j];
		}
	        wssten->sl[0] = Left_state_at_node_of_o_curve(&Tempcurt);
		for (i = 1; i < nsts; ++i)
		{
		    sa[isgn*(i-1)] = wssten->sl[i];
		    wssten->lcrds[i] = Coords(p[isgn*(i-1)]);
		}
	        states_at_distance_along_curve(tnp,tempbt,Tempcurt.curve,
					       Tempcurt.orient,ds,nsts-1,
					       sa,NULL,NULL,NULL,NULL,p,fr);
		ans = ansl;
	    }
	}

	    	/* update using right tempcurt (= interior) states */

	else if (comp_is_interior_on_curve(positive_component(Tempcurt.curve),
				           Tempcurn.curve) == YES)
	{
	    if (Tempcurt.orient == Tempcurn.orient)
	    {
		find_states_for_comp(positive_component(Tempcurt.curve),
				     Tempcurt.curve,Tempcurn.curve,tempn,
				     wssten->sr,nsts,wave,fr);
		for (i = 1; i < nsts; ++i)
		{
		    for (j = 0; j < dim; ++j)
			wssten->rcrds[i][j] = wssten->coords[j] + i*dn*nor[j];
		}
	        wssten->sl[0] = Right_state_at_node_of_o_curve(&Tempcurt);
		for (i = 1; i < nsts; ++i)
		{
		    sa[isgn*(i-1)] = wssten->sl[i];
		    wssten->lcrds[i] = Coords(p[isgn*(i-1)]);
		}
	        states_at_distance_along_curve(tnp,tempbt,Tempcurt.curve,
					       Tempcurt.orient,ds,nsts-1,
					       NULL,sa,NULL,NULL,NULL,p,fr);
		ans = ansl;
	    }
	    else
	    {
		find_states_for_comp(positive_component(Tempcurt.curve),
				     Tempcurt.curve,Tempcurn.curve,tempn,
				     wssten->sl,nsts,wave,fr);
		for (i = 1; i < nsts; ++i)
		{
		    for (j = 0; j < dim; ++j)
			wssten->lcrds[i][j] = wssten->coords[j] - i*dn*nor[j];
		}
	        wssten->sr[0] = Right_state_at_node_of_o_curve(&Tempcurt);
		for (i = 1; i < nsts; ++i)
		{
		    sa[isgn*(i-1)] = wssten->sr[i];
		    wssten->rcrds[i] = Coords(p[isgn*(i-1)]);
		}
	        states_at_distance_along_curve(tnp,tempbt,Tempcurt.curve,
					       Tempcurt.orient,ds,nsts-1,
					       NULL,sa,NULL,NULL,NULL,p,fr);
		for (i = 1; i < nsts; ++i)
		    wssten->rcrds[i] = Coords(p[i-1]);
		ans = ansr;
	    }
	}
	else
	{
	    screen("ERROR in snd_fixed_node_propagate, "
	           "(comp_is_interior_on_curve("
		   "negative_component(Tempcurt.curve),Tempcurn.curve) != YES) "
		   "and "
	           "(comp_is_interior_on_curve("
		   "positive_component(Tempcurt.curve),Tempcurn.curve) != YES)"
		   "\n");

	    (void) printf("exterior_component(Tempcurn.curve->interface) = "
	                  "%d\n",exterior_component(Tempcurn.curve->interface));
	    (void) printf("Tempcurt.curve\n");
	    print_curve(Tempcurt.curve);
	    (void) printf("Tempcurn.curve\n");
	    print_curve(Tempcurn.curve);
	    (void) printf("oldn\n");
	    print_node(oldn);
	    (void) printf("tempn\n");
	    print_node(tempn);
	    (void) printf("newn\n");
	    print_node(newn);
	    (void) printf("oldn->interface\n");
	    print_interface(oldn->interface);
	    (void) printf("tempn->interface\n");
	    print_interface(tempn->interface);
	    (void) printf("newn->interface\n");
	    print_interface(newn->interface);
	    clean_up(ERROR);
	}
	npt_w_speed(wssten,ansl,ansr,V);

	Check_return(
	    find_correspond_of_oriented_curve(&Tempcurt,&Newcurt,newn,
					      fr,newn->interface),
	    snd_fixed_node_propagate)
	Check_return(
	    find_correspond_of_oriented_curve(&Tempcurn,&Newcurn,newn,
					      fr,newn->interface),
	    snd_fixed_node_propagate)
	assign_states_around_node(&Newcurt,&Tempcurt,ans,sizest);
	assign_states_around_node(&Newcurn,&Tempcurn,ans,sizest);

	if (debugging("fixed_node"))
	{
	    (void) printf("left/right states on newcurt/newcurn:\n");
	    (*fr->print_state)(Left_state_at_node_of_o_curve(&Newcurt));
	    (*fr->print_state)(Right_state_at_node_of_o_curve(&Newcurt));
	    (*fr->print_state)(Left_state_at_node_of_o_curve(&Newcurn));
	    (*fr->print_state)(Right_state_at_node_of_o_curve(&Newcurn));
	}
	debug_print("fixed_node","Leaving snd_fixed_node_propagate()\n");
	return YES;
}		/*end snd_fixed_node_propagate*/

/* choice of version: */

#undef QUADRANGLE
#undef FLAT_CORNER
#define ENTROPY

#if defined(QUADRANGLE)

/*
*			corner_node_propagate():
*
*	Calculates the states at a corner of two walls meeting at
*	an arbitrary angle. The function calculates the states
*	by using a 2d scheme obtained by flux calculation in a quadrangle.
*/

/* ARGSUSED */
LOCAL void corner_node_propagate(
	NODE		*node,
	NODE		*tempnode,
	O_CURVE		*c1,
	O_CURVE		*c2,
	Front		*fr,
	Wave		*wave,
	Locstate	ans,
	double		dt)
{
	COMPONENT	comp;
	double		n1[MAXD], n2[MAXD];
	double		t1[MAXD], t2[MAXD];
	double		coords[MAXD];
	double		nor1[MAXD],nor2[MAXD];
	double		ds1,ds2;
	double		area,factor;
	Locstate	s0,s1,s2;
	int		dim = node->interface->dim;
	static Locstate s3 = NULL,s1l = NULL,s1r = NULL,s2l = NULL,s2r = NULL;

	/* Allocate storage */

	if (s2r == NULL)
	{
	    alloc_state(fr->interf,&s3,fr->sizest);
	    alloc_state(fr->interf,&s1l,fr->sizest);
	    alloc_state(fr->interf,&s1r,fr->sizest);
	    alloc_state(fr->interf,&s2l,fr->sizest);
	    alloc_state(fr->interf,&s2r,fr->sizest);
	}


	/* define state at node */

	if (is_obstacle_state(Right_state_at_node(cur1,cur1_orient)))
	{
	    s0 = Left_state_at_node(cur1,cur1_orient);
	    comp = negative_component(cur1);
	}
	else
	{
	    s0 = Right_state_at_node(cur1,cur1_orient);
	    comp = positive_component(cur1);
	}

	/* calculate normals at node pointing into region of flow */
	/* calculate tangents at node pointing out from corner */

	normal(node->posn,Hyper_surf_element(Bond_at_node(cur1,cur1_orient)),
	       Hyper_surf(cur1),n1,fr);
	find_tangent_to_curve(node->posn,Bond_at_node(cur1,cur1_orient),cur1,
		              cur1_orient,t1,fr);
	if (negative_component(cur1) == comp)
	{
	    n1[0] *= -1.;
	    n1[1] *= -1.;
	}

	normal(node->posn,Hyper_surf_element(Bond_at_node(cur2,cur2_orient)),
	       Hyper_surf(cur2),n2,fr);
	find_tangent_to_curve(node->posn,Bond_at_node(cur2,cur2_orient),cur2,
		              cur2_orient,t2,fr);
	if (negative_component(cur2) == comp)
	{
	    n2[0] *= -1.;
	    n2[1] *= -1.;
	}


	/* define states on walls */

	ds1 = grid_size_in_direction(t1,fr->rect_grid->h,dim);
	states_at_distance_along_curve(node->posn,
		                       Bond_at_node(cur1,cur1_orient),
				       cur1,cur1_orient,ds1,1,
				       &s1l,&s1r,NULL,NULL,NULL,NULL,fr);
	s1 = (negative_component(cur1) == comp)	? s1l : s1r;

	ds2 = grid_size_in_direction(t2,fr->rect_grid->h,dim);
	states_at_distance_along_curve(node->posn,
		                       Bond_at_node(cur2,cur2_orient),
				       cur2,cur2_orient,ds2,1,
				       &s2l,&s2r,NULL,NULL,NULL,NULL,fr);
	s2 = (negative_component(cur2) == comp)	? s2l : s2r;

	/* define state at interior corner of quadrangle */

	coords[0] = Coords(node->posn)[0] + ds1 * t1x + ds2 * t2x;
	coords[1] = Coords(node->posn)[1] + ds1 * t1y + ds2 * t2y;
	hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,fr,wave,s3,NULL);

	/* increment by flux */

	area = ds1 * ds2 * fabs(t1x*t2y - t1y*t2x);
	factor = - dt / area;
	ft_assign(ans,s0,fr->sizest);
	nor1[1] = factor*ds2*n2[1];	nor1[0] = factor*ds2*n2[0];
	nor2[1] = factor*ds1*n1[1];	nor1[0] = factor*ds1*n1[0];
	flux_across_line_segment(wave_tri_soln(wave),s1,s3,nor1,ans);
	flux_across_line_segment(wave_tri_soln(wave),s3,s2,nor2,ans);
	if (debugging("fixed_node"))
	{
	    static Locstate test_ans = NULL;

	    if (test_ans == NULL)
	        alloc_state(fr->interf,&test_ans,fr->sizest);

	    ft_assign(test_ans,ans,fr->sizest);
	    nor1[1] = -factor*ds2*n2[1];	nor1[0] = -factor*ds2*n2[0];
	    nor2[1] = -factor*ds1*n1[1];	nor1[0] = -factor*ds1*n1[0];
	    flux_across_line_segment(wave_tri_soln(wave),s2,s0,nor1,test_ans);
	    flux_across_line_segment(wave_tri_soln(wave),s0,s1,nor2,test_ans);
	    (void) printf("test answer at corner:\n");
	    (*fr->print_state)(test_ans);
	}
	zero_state_velocity(ans,dim);

	if (debugging("fixed_node"))
	{
	    (void) printf("ds1= %g ds2= %g\n",ds1,ds2);
	    (void) printf("x0= %g y0= %g \n",
			  Coords(node->posn)[0],Coords(node->posn)[1]);
	    (void) printf("x1= %g y1= %g \n",
			  Coords(node->posn)[0] + ds1 * t1x,
			  Coords(node->posn)[1] + ds1 * t1y);
	    (void) printf("x2= %g y2= %g \n",
			  Coords(node->posn)[0] + ds2 * t2x,
			  Coords(node->posn)[1] + ds2 * t2y);
	    (void) printf("x3= %g y3= %g \n",
			  Coords(node->posn)[0] + ds1 * t1x + ds2 * t2x,
			  Coords(node->posn)[1] + ds1 * t1y + ds2 * t2y);
	    (void) printf("corresponding stencil states:\n");
	    (*fr->print_state)(s0);
	    (*fr->print_state)(s1);
	    (*fr->print_state)(s2);
	    (*fr->print_state)(s3);
	    (void) printf("area= %g dt= %g factor= %g\n",area,dt,factor);
	    (void) printf("answer at corner:\n");
	    (*fr->print_state)(ans);
	}
	debug_print("fixed_node","Leaving walls_corner_propagate()\n");
}		/*end corner_node_propagate*/
#endif /* defined(QUADRANGLE) */


#if defined(FLAT_CORNER)

/*
*			corner_node_propagate():
*
*	Calculates the states at a corner of two walls meeting at
*	an arbitrary angle. The function calculates the states
*	by a normal-tangential sweep where the normal is in the
*	direction of the bisector. This can (?) be used for an
*	angle close to 180 deg.
*/

LOCAL void corner_node_propagate(
	NODE		*node,
	NODE		*tempnode,
	O_CURVE		*c1,
	O_CURVE		*c2,
	Front		*fr,
	Wave		*wave,
	Locstate	ans,
	double		dt)
{
	COMPONENT	comp;
	double		n1[MAXD], n2[MAXD], magn;
	double		t1[MAXD], t2[MAXD], tgnt[MAXD],magt,ds;
	double		v[MAXD];
	int		inverted;
	int		i, j, dim = fr->interf->dim;
	O_CURVE		Tempc1, Tempc2;
	BOND		*b1,*b2;		/* bonds at temp node */
	BOND		*tempb1,*tempb2;	/* bonds at temp node */
	POINT		pprev,pnext;
	Locstate	s0, ansl, ansr;
	size_t		sizest = fr->sizest;
	static WSSten	*wssten = NULL;
	static int	nrad = 0;
	static Tan_stencil *sten = NULL;
	static Locstate s1 = NULL;
	static Locstate sprevin = NULL,sprevex = NULL;	/* prev states for tan sweep */
	static Locstate snextin = NULL,snextex = NULL;	/* next states for tan sweep */
	static Locstate snormex = NULL,snormin = NULL;	/* answer of normal sweep */
	Locstate	obst_state = return_obst_state();

	/* Allocate storage */

	if (snormex == NULL)
	{
	    nrad = fr->npts_tan_sten/2;
	    wssten = AllocDefaultWSSten(fr);
	    alloc_state(fr->interf,&s1,sizest);
	    sten = alloc_tan_stencil(fr,nrad);
	    set_default_tan_stencil(sten);
	    alloc_state(fr->interf,&sprevin,sizest);
	    alloc_state(fr->interf,&sprevex,sizest);
	    alloc_state(fr->interf,&snextin,sizest);
	    alloc_state(fr->interf,&snextex,sizest);
	    alloc_state(fr->interf,&snormin,sizest);
	    alloc_state(fr->interf,&snormex,sizest);
	}
	else
	    ClearWSStenData(wssten);

	debug_print("fixed_node","Entering corner_node_propagate()\n");

	/* define state at node */

	if (is_obstacle_state(Right_state_at_node_of_o_curve(c1)))
	{
	    s0 = Left_state_at_node_of_o_curve(c1);
	    comp = negative_component(c1->curve);
	}
	else
	{
	    s0 = Right_state_at_node_of_o_curve(c1);
	    comp = positive_component(c1->curve);
	}

	normal(oldn->posn,Hyper_surf_element(Bond_at_node_of_o_curve(c1)),
	       Hyper_surf(c1->curve),n1,fr);
	find_tangent_to_curve(oldn->posn,Bond_at_node_of_o_curve(c1),c1->curve,
		              c1->orient,t1,fr);
	if (negative_component(c1->curve) == comp)
	    for (i = 0; i < dim; ++i)
		n1[i] *= -1.;

	normal(oldn->posn,Hyper_surf_element(Bond_at_node_of_o_curve(c2)),
	       Hyper_surf(c2->curve),n2,fr);
	find_tangent_to_curve(oldn->posn,Bond_at_node_of_o_curve(c2),c2->curve,
		              c2->orient,t2,fr);
	if (negative_component(c2->curve) == comp)
	    for (i = 0; i < dim; ++i)
		n2[i] *= -1.;

	for (i = 0; i < dim; ++i)
	{
	    wssten->nor[i] = n1[i] + n2[i];
	    tgnt[i] = t2[i] - t1[i];
	}
	magn = mag_vector(n,dim);
	magt = mag_vector(t,dim);
	for (i = 0; i < dim; ++i)
	{
	    wssten->nor[i] /= magn;	tgnt[i] /= magt;
	    Coords(sten->p[0]) = Coords(tempn->posn)[i];
	}

	/* define state at distance dn from corner */

	wssten->dn = grid_size_in_direction(wssten->nor,fr->rect_grid->h,dim);
	wssten->coords = Coords(tempn->posn);
	wssten->lcrds[0] = wssten->rcrds[0] = wssten->coords;
	for (i = 1; i < wssten->nsts; ++i)
	{
	    for (j = 0; j < dim; ++j)
	    {
	        wssten->rcrds[i][j] =
		    Coords(oldn->posn)[j] + i*wssten->dn*wssten->nor[i];
	        wssten->lcrds[i][j] =
		    Coords(oldn->posn)[j] - i*wssten->dn*wssten->nor[i];
	    }
	}

	wssten->dt = dt;
	wssten->w_type = wave_type(c1->curve);
	wssten->pjump = 0.0;
	wssten->ncomp = negative_component(c1->curve);
	wssten->pcomp = positive_component(c1->curve);
	wssten->hs = Hyper_surf(c1->curve);
	wssten->front = fr;
	wssten->wave = wave;
	if (is_obstacle_state(Left_state_at_node_of_o_curve(c1)))
	{
	    ansl = sten->leftst[0];
	    ansr = sten->rightst[0];
	    for (i = 0; i < wssten->nsts; ++i)
		wssten->sl[i] = obst_state;
	    wssten->sr[0] = s0;
	    for (i = 1; i < wssten->nsts; ++i)
	    {
		wssten->sr[i] = (Locstate) (wssten->sr_store + i*sizest);
	        hyp_solution(wssten->rcrds[i],comp,NULL,UNKNOWN_SIDE,
			     fr,wave,wssten->sr[i],s0);
	    }
	}
	else
	{
	    ansl = sten->rightst[0];
	    ansr = sten->leftst[0];
	    wssten->sl[0] = s0;
	    for (i = 1; i < wssten->nsts; ++i)
	    {
		wssten->sl[i] = (Locstate) (wssten->sl_store + i*sizest);
	        hyp_solution(wssten->lcrds[i],comp,NULL,UNKNOWN_SIDE,
			     fr,wave,wssten->sl[i],s0);
	    }
	    for (i = 0; i < wssten->nsts; ++i)
		wssten->sr[i] = obst_state;
	}
	npt_w_speed(wssten,ansl,ansr,v);

	/*
	 * find prev and next states on tempfront;
	 * call one_side_npt_tang_solver()
	 */

	find_correspond_of_oriented_curve(c1,&Tempc1,tempn,fr,
					  tempn->interface);
	find_correspond_of_oriented_curve(c2,&Tempc2,tempn,fr,
					  tempn->interface);
	tempb1 = Bond_at_node_of_o_curve(&Tempc1);
	tempb2 = Bond_at_node_of_o_curve(&Tempc2);

	ds = grid_size_in_direction(t1,fr->rect_grid->h,dim);
	inverted = NO;
	if (c1->orient == POSITIVE_ORIENTATION)
	{
	    inverted = YES;
	    invert_curve(Tempc1.curve);
	    Tempc1.orient = Opposite_orient(Tempc1.orient);
	}
	if (positive_component(c1->curve) == comp)
	    states_at_distance_along_curve(tempn->posn,tempb1,Tempc1.curve,
					   Tempc1.orient,ds,nrad,
					   sten->leftst-1,sten->rightst-1,
			                   sten->hs-1,sten->hse-1,sten->t-1,
			                   sten->p-1,fr);
	else
	    states_at_distance_along_curve(tempn->posn,tempb1,Tempc1.curve,
					   Tempc1.orient,ds,nrad,
			                   sten->rightst-1,sten->leftst-1,
			                   sten->hs-1,sten->hse-1,sten->t-1,
			                   sten->p-1,fr);
	if (inverted == YES)
	{
	    invert_curve(Tempc1.curve);
	    Tempc1.orient = Opposite_orient(Tempc1.orient);
	}

	ds = grid_size_in_direction(t2,fr->rect_grid->h,dim);
	inverted = NO;
	if (Tempc2.orient == POSITIVE_ORIENTATION)
	{
	    inverted = YES;
	    invert_curve(Tempc2.curve);
	    Tempc2.orient = Opposite_orient(Tempc2.orient);
	}
	if (positive_component(c2->curve) == comp)
	    states_at_distance_along_curve(tempn->posn,tempb2,Tempc2.curve,
					   Tempc2.orient,ds,nrad,
			                   sten->leftst+1,sten->rightst+1,
			                   sten->hs+1,sten->hse+1,sten->t+1,
			                   sten->p+1,fr);
	else
	    states_at_distance_along_curve(tempn->posn,tempb2,Tempc2.curve,
					   Tempc2.orient,ds,nrad,
			                   sten->rightst+1,sten->leftst+1,
			                   sten->hs+1,sten->hse+1,sten->t+1,
			                   sten->p+1,fr);
	if (inverted == YES)
	{
	    invert_curve(Tempc2.curve);
	    Tempc2.orient = Opposite_orient(Tempc2.orient);
	}


	sten->newhs = NULL;
	sten->comp = comp;
	sten->states = sten->rightst;
	sten->dir = tgnt;
	one_side_npt_tang_solver(ds,dt,sten,ans,fr);
	if (debugging("fixed_node"))
	{
	    (void) printf("ds= %g\n",ds);
	    (void) printf("x0 = %g y0 = %g\n",
		          Coords(oldn->posn)[0],Coords(oldn->posn)[1]);
	    (void) printf("stencil states for normal sweep, s0 and s1:\n");
	    (*fr->print_state)(s0);
	    (*fr->print_state)(s1);
	    (void) printf("stencil states for tangential sweep:\n");
	    (*fr->print_state)(sprevin);
	    (*fr->print_state)(snormin);
	    (*fr->print_state)(snextin);
	    (void) printf("answer at corner:\n");
	    (*fr->print_state)(ans);
	}

	debug_print("fixed_node","Leaving corner_node_propagate()\n");
}		/*end corner_node_propagate*/
#endif /* defined(FLAT_CORNER) */

#if defined(ENTROPY)

/* ARGSUSED */
LOCAL void corner_node_propagate(
	NODE		*node,
	NODE		*tempnode,
	O_CURVE		*c1,
	O_CURVE		*c2,
	Front		*fr,
	Wave		*wave,
	Locstate	ans,
	double		dt)
{
	COMPONENT	comp;
	double		charact_speed1,charact_speed2,pr;
	double		n1[MAXD], n2[MAXD];
	double		t1[MAXD], t2[MAXD];
	double		ds1,ds2;
	Locstate	s0,s1,s2, st_r0;
	int		i, dim = node->interface->dim;
	static Locstate s3 = NULL,s3l = NULL,s3r = NULL,
	                s1l = NULL,s1r = NULL,s2l = NULL,s2r = NULL;

	/* Allocate storage */

	if (s3r == NULL)
	{
	    alloc_state(fr->interf,&s1l,fr->sizest);
	    alloc_state(fr->interf,&s1r,fr->sizest);
	    alloc_state(fr->interf,&s2l,fr->sizest);
	    alloc_state(fr->interf,&s2r,fr->sizest);
	    alloc_state(fr->interf,&s3,fr->sizest);
	    alloc_state(fr->interf,&s3l,fr->sizest);
	    alloc_state(fr->interf,&s3r,fr->sizest);
	}

	/* define state at node */

	if (is_obstacle_state(Right_state_at_node_of_o_curve(c1)))
	{
	    s0 = Left_state_at_node_of_o_curve(c1);
	    comp = negative_component(c1->curve);
	}
	else
	{
	    s0 = Right_state_at_node_of_o_curve(c1);
	    comp = positive_component(c1->curve);
	}

	/* calculate normals at node pointing into region of flow */
	/* calculate tangents at node pointing out from corner */

	normal(node->posn,Hyper_surf_element(Bond_at_node_of_o_curve(c1)),
	       Hyper_surf(c1->curve),n1,fr);
	find_tangent_to_curve(node->posn,Bond_at_node_of_o_curve(c1),c1->curve,
	                      c1->orient,t1,fr);
	if (negative_component(c1->curve) == comp)
	    for (i = 0; i < dim; ++i)
		n1[i] *= -1.;

	normal(node->posn,Hyper_surf_element(Bond_at_node_of_o_curve(c2)),
	       Hyper_surf(c2->curve),n2,fr);
	find_tangent_to_curve(node->posn,Bond_at_node_of_o_curve(c2),c2->curve,
		              c2->orient,t2,fr);
	if (negative_component(c2->curve) == comp)
	    for (i = 0; i < dim; ++i)
		n2[i] *= -1.;

	/* define states on walls */

	ds1 = grid_size_in_direction(t1,fr->rect_grid->h,dim);
	states_at_distance_along_curve(node->posn,
		                       Bond_at_node_of_o_curve(c1),c1->curve,
				       c1->orient,ds1,1,
				       &s1l,&s1r,NULL,NULL,NULL,NULL,fr);
	s1 = (negative_component(c1->curve) == comp) ? s1l : s1r;

	ds2 = grid_size_in_direction(t2,fr->rect_grid->h,dim);
	states_at_distance_along_curve(node->posn,
		                       Bond_at_node_of_o_curve(c2),c2->curve,
				       c2->orient,ds2,1,&s2l,&s2r,
				       NULL,NULL,NULL,NULL,fr);
	s2 = (negative_component(c2->curve) == comp) ? s2l : s2r;


	charact_speed1 = -0.5*scalar_product(Mom(s1),t1,dim)/Dens(s1);
	charact_speed2 = -0.5*scalar_product(Mom(s2),t2,dim)/Dens(s2);

	if (charact_speed1 < 0. && charact_speed2 < 0.)
	{
	    st_r0 = s0;
	}
	else if (charact_speed1 >= 0. && charact_speed1 >= charact_speed2)
	{
	    states_at_distance_along_curve(node->posn,
			                   Bond_at_node_of_o_curve(c1),
					   c1->curve,c1->orient,
			                   charact_speed1*dt,1,&s3l,&s3r,
			                   NULL,NULL,NULL,NULL,fr);
	    s3 = (negative_component(c1->curve) == comp) ? s3l : s3r;
	    st_r0 = s3;
	}
	else /* charact_speed2 >= 0. && charact_speed2 >= charact_speed1 */
	{
	    states_at_distance_along_curve(node->posn,
			                   Bond_at_node_of_o_curve(c2),
					   c2->curve,c2->orient,
					   charact_speed2*dt,1,&s3l,&s3r,
					   NULL,NULL,NULL,NULL,fr);
	    s3 = (negative_component(c2->curve) == comp) ? s3l : s3r;
	    st_r0 = s3;
	}

	Set_params(ans,s0);
	zero_state_velocity(ans,dim);
	pr = 0.5 * (pressure(s1) + pressure(s2));
	state_on_adiabat_with_pr(st_r0,pr,ans,GAS_STATE);

	if (debugging("fixed_node"))
	{
	    double ftmp[MAXD];

	    (void) printf("ds1= %g ds2= %g\n",ds1,ds2);
	    print_general_vector("node->posn = ",
	    	                 Coords(node->posn),dim,"\n");
	    for (i = 0; i < dim; ++i)
	    	ftmp[i] = Coords(node->posn)[i] + ds1 * t1[i];
	    print_general_vector("p1 = ",ftmp,dim,"\n");
	    for (i = 0; i < dim; ++i)
	    	ftmp[i] = Coords(node->posn)[i] + ds2 * t2[i];
	    print_general_vector("p2 = ",ftmp,dim,"\n");
	    (void) printf("corresponding stencil states:\n");
	    verbose_print_state("old_corner",s0);
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    (void) printf("answer at corner:\n");
	    verbose_print_state("new_corner",ans);
	}
	debug_print("fixed_node","Leaving corner_node_propagate()\n");
}		/*end corner_node_propagate*/

#endif /* defined(ENTROPY) */


LOCAL int find_active_curves(
	NODE		*node,
	O_CURVE		*activec)
{
	int		n_act;
	CURVE		**c;

	n_act = 0;

	if (node->in_curves != NULL)
	{
	    for (c = node->in_curves; *c; ++c)
	    {
	    	if (wave_type(*c) != PASSIVE_BOUNDARY) 
	    	{
	    	    activec[n_act].curve = *c;
	    	    activec[n_act].orient = NEGATIVE_ORIENTATION;
	    	    ++n_act;
	    	}
	    }
	}

	if (node->out_curves != NULL)
	{
	    for (c = node->out_curves; *c; ++c)
	    {
	    	if (wave_type(*c) != PASSIVE_BOUNDARY) 
	    	{
	    	    activec[n_act].curve = *c;
	            activec[n_act].orient = POSITIVE_ORIENTATION;
		    ++n_act;
		}
	    }
	}
	return n_act;
}		/*end find_active_curves*/

LOCAL void assign_states_around_node(
	O_CURVE		*newc,
	O_CURVE		*oldc,
	Locstate	ans,
	size_t		sizest)
{
	if (is_obstacle_state(left_start_state(oldc->curve)))
	{
	    obstacle_state(oldc->curve->interface,
			   Left_state_at_node_of_o_curve(newc),sizest);
	}
	else
	{
	    ft_assign(Left_state_at_node_of_o_curve(newc),ans,sizest);
	}

	if (is_obstacle_state(right_start_state(oldc->curve)))
	{
	    obstacle_state(oldc->curve->interface,
			   Right_state_at_node_of_o_curve(newc),sizest);
	}
	else
	{
	    ft_assign(Right_state_at_node_of_o_curve(newc),ans,sizest);
	}
}		/*end assign_states_around_node*/


LOCAL	boolean	comp_is_interior_on_curve(
	COMPONENT	comp,
	CURVE		*c)
{
	INTERFACE	*intfc = c->interface;

	if (is_excluded_comp(comp,intfc))
	    return NO;
	return ((comp==negative_component(c)) ||
		(comp==positive_component(c))) ? YES : NO;
}		/*end comp_is_interior_on_curve*/

LOCAL void find_states_for_comp(
	COMPONENT	comp,
	CURVE		*tc,
	CURVE		*nc,
	NODE		*node,
	Locstate	*state,
	int		nsts,
	Wave		*wave,
	Front		*fr)
{
	int	i;

	if (wave_type(nc) == DIRICHLET_BOUNDARY)
	{

/* TODO: get rid of the use of these */
#define left_state_at_node(n,c)						\
	(((n) == (c)->start) ? left_start_state(c) : left_end_state(c))
#define right_state_at_node(n,c)					\
	(((n) == (c)->start) ? right_start_state(c) : right_end_state(c))
#define state_at_node_with_comp(n,c,comp) (				\
	((comp)==negative_component(c)) ? left_state_at_node(n,c) :	\
	    (((comp)==positive_component(c)) ? right_state_at_node(n,c) : NULL))

	    if (comp == negative_component(tc))
		state[0] = left_state_at_node(node,tc);
	    else if (comp == positive_component(tc))
		state[0] = right_state_at_node(node,tc);
	    else
		state[0] = NULL;
	    if (boundary_state_function(nc) != NULL)
	    	(*boundary_state_function(nc))(Coords(node->posn),
					       Hyper_surf(nc),fr,
					       (POINTER)wave,state[1]);
	    else
	    {
	    	set_state(state[1],state_type(boundary_state(nc)),
			  boundary_state(nc));
	    }

#undef left_state_at_node
#undef right_state_at_node
	    for (i = 2; i < nsts; ++i)
	        set_state(state[i],state_type(state[1]),state[1]);
	}
	else		/* defined(NEUMANN_BOUNDARY) */
	{
	    for (i = 0; i < nsts; ++i)
	        obstacle_state(fr->interf,state[i],fr->sizest);
	}
}		/*end find_states_for_comp*/




/* ARGSUSED */
LOCAL boolean snd_B_node_propagate(
	Front		*fr,
	Front		*newfr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	INTERFACE       *intfc = fr->interf;
	ORIENTATION     ci_or, cl_or, cr_or; /* curve orientations */
	int             i, j, dim = fr->rect_grid->dim;
	int	        isgn;
	double           t[MAXD];      /* ci tan points to interior */
	double           l_n[MAXD];    /* ci tan oriented by cl nor */
	double           r_n[MAXD];    /* ci tan oriented by cr nor */
	double           V[MAXD];
	CURVE           *tempci,*newci; /* incident curve */
	CURVE           *tempcl,*newcl; /* wall left of inc curve */
	CURVE           *tempcr,*newcr; /* wall right of inc curve */
	Locstate        ansl,ansr;
	size_t	        sizest = fr->sizest;
	static WSSten	*wssten = NULL;
	static Locstate *sleft = NULL, *sright = NULL;       /* left/right states of ci */
	static Locstate stemp1 = NULL, stemp2 = NULL;     /* workspace */
	static POINT	**p = NULL;
	static int	nsts = 0;

	debug_print("B_node","Entering snd_B_node_propagate()\n");
#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("B_node")) 
	{
	    (void) printf("\n\tOLD NODE (interface = %llu):\n",
	           interface_number(oldn->interface));
	    print_node(oldn);
	    (void) printf("\n\tTEMP NODE (interface = %llu):\n",
	           interface_number(tempn->interface));
	    print_node(tempn);
	    (void) printf("\n\tNEW NODE (interface = %llu):\n",
	           interface_number(newn->interface));
	    print_node(newn);
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    /* Allocate storage */

	if (wssten == NULL) 
	{
	    wssten = AllocDefaultWSSten(fr);
	    nsts = wssten->nsts;
	    uni_array(&sleft,2*nsts-1,sizeof(Locstate*));
	    sleft += nsts-1;
	    uni_array(&sright,2*nsts-1,sizeof(Locstate*));
	    sright += nsts-1;
	    uni_array(&p,2*nsts-1,sizeof(POINT*));
	    p += nsts-1;
	    for (i = -nsts+1; i < nsts; ++i)
	    {
		p[i] = Static_point(fr->interf);
		sleft[i] = left_state(p[i]);
		sright[i] = right_state(p[i]);
	    }
	    alloc_state(fr->interf,&stemp1,sizest);
	    alloc_state(fr->interf,&stemp2,sizest);
	}
	else
	{
	    ClearWSStenData(wssten);
	    for (i = -nsts+1; i < nsts; ++i)
	    {
	        obstacle_state(fr->interf,left_state(p[i]),sizest);
	        obstacle_state(fr->interf,right_state(p[i]),sizest);
	    }
	}


	    /* Identify curves */

	find_curve_with_status(tempn,&tempci,&ci_or,INCIDENT);
	find_curve_with_status(newn,&newci,&ci_or,INCIDENT);
	if (ci_or == NEGATIVE_ORIENTATION) 
	{
	    tempcl = adjacent_curve(tempci,ci_or,CLOCKWISE,&cl_or);
	    newcl = adjacent_curve(newci,ci_or,CLOCKWISE,&cl_or);
	    tempcr = adjacent_curve(tempci,ci_or,COUNTER_CLOCK,&cr_or);
	    newcr = adjacent_curve(newci,ci_or,COUNTER_CLOCK,&cr_or);
	}
	else 
	{
	    tempcl = adjacent_curve(tempci,ci_or,COUNTER_CLOCK,&cl_or);
	    newcl = adjacent_curve(newci,ci_or,COUNTER_CLOCK,&cl_or);
	    tempcr = adjacent_curve(tempci,ci_or,CLOCKWISE,&cr_or);
	    newcr = adjacent_curve(newci,ci_or,CLOCKWISE,&cr_or);
	}

	if (
	    (is_scalar_wave(wave_type(tempci))) &&
	    (    (wave_type(tempcl) ==NEUMANN_BOUNDARY) ||
	         (wave_type(tempcr)==NEUMANN_BOUNDARY)    )
	)
	{
	    boolean status;
	    status = snd_contact_wall_node_propagate(fr,newfr,wave,
						     oldn,tempn,newn,dt);
	    debug_print("B_node","Leaving snd_B_node_propagate()\n");
	    return status;
	}

	    /* Tangent to ci is oriented out from node */

	find_tangent_to_curve(Node_of(tempci,ci_or)->posn,
			      Bond_at_node(tempci,ci_or),tempci,ci_or,t,fr);

	    /* Find l/r states on tempci */

	wssten->dn = grid_size_in_direction(t,fr->rect_grid->h,dim);
	isgn = (ci_or == POSITIVE_ORIENTATION) ? 1 : -1;
	states_at_distance_along_curve(tempn->posn,Bond_at_node(tempci,ci_or),
				       tempci,ci_or,wssten->dn,nsts-1,
				       sleft,sright,NULL,NULL,NULL,p,fr);
	
	wssten->front = fr;
	wssten->wave = wave;
	wssten->dt = dt;
	wssten->pjump = 0.0;
	wssten->coords = Coords(tempn->posn);
	wssten->lcrds[0] = wssten->rcrds[0] = wssten->coords;
	if (node_type(tempn) == NEUMANN_NODE) 
	{

	    /* Wall interaction */

	    normal(tempn->posn,Hyper_surf_element(Bond_at_node(tempcl,cl_or)),
		   Hyper_surf(tempcl),wssten->nor,fr);
	    if (ci_or != cl_or) 
	    {
	        wssten->sl[0] = Left_state_at_node(tempci,ci_or);
		for (i = 1; i < nsts; ++i)
		{
		    wssten->lcrds[i] = Coords(p[isgn*(i-1)]);
		    wssten->sl[i] = sleft[isgn*(i-1)];
		}
		wssten->sr[0] = return_obst_state();
		for (i = 1; i < nsts; ++i)
		{
		    wssten->sr[i] = wssten->sr[0];
		    for (j = 0; j < dim; ++j)
		    {
			wssten->rcrds[i][j] =
			    wssten->coords[j] + i*wssten->dn*wssten->nor[j];
		    }
		}
	        ansl = Left_state_at_node(newcl,cl_or);
	    }
	    else 
	    {
		wssten->sl[0] = return_obst_state();
		for (i = 1; i < nsts; ++i)
		{
		    wssten->sl[i] = wssten->sl[0];
		    for (j = 0; j < dim; ++j)
		    {
			wssten->lcrds[i][j] =
			    wssten->coords[j] - i*wssten->dn*wssten->nor[j];
		    }
		}
		wssten->sr[0] = Left_state_at_node(tempci,ci_or);
		for (i = 1; i < nsts; ++i)
		{
		    wssten->rcrds[i] = Coords(p[isgn*(i-1)]);
		    wssten->sr[i] = sleft[isgn*(i-1)];
		}
	        ansl = Right_state_at_node(newcl,cl_or);
	    }

#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("B_node"))
	    {
		PrintWSSten(wssten);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    wssten->ncomp = negative_component(newcl);
	    wssten->pcomp = positive_component(newcl);
	    wssten->w_type = wave_type(tempcl);
	    wssten->hs = Hyper_surf(tempcl);
	    npt_w_speed(wssten,Left_state_at_node(newcl,cl_or),
	                Right_state_at_node(newcl,cl_or),V);
	    ft_assign(Left_state_at_node(newci,ci_or),ansl,sizest);

	    normal(tempn->posn,Hyper_surf_element(Bond_at_node(tempcr,cr_or)),
		   Hyper_surf(tempcr),wssten->nor,fr);
	    if (ci_or != cr_or) 
	    {
		wssten->sl[0] = return_obst_state();
		for (i = 1; i < nsts; ++i)
		{
		    wssten->sl[i] = wssten->sl[0];
		    for (j = 0; j < dim; ++j)
		    {
			wssten->lcrds[i][j] =
			    wssten->coords[j] - i*wssten->dn*wssten->nor[j];
		    }
		}
		wssten->sr[0] = Right_state_at_node(tempci,ci_or);
		for (i = 1; i < nsts; ++i)
		{
		    wssten->rcrds[i] = Coords(p[isgn*(i-1)]);
		    wssten->sr[i] = sright[isgn*(i-1)];
		}
	        ansr = Right_state_at_node(newcr,cr_or);
	    }
	    else 
	    {
	        wssten->sl[0] = Right_state_at_node(tempci,ci_or);
		for (i = 1; i < nsts; ++i)
		{
		    wssten->lcrds[i] = Coords(p[isgn*(i-1)]);
		    wssten->sl[i] = sright[isgn*(i-1)];
		}
		wssten->sr[0] = return_obst_state();
		for (i = 1; i < nsts; ++i)
		{
		    wssten->sr[i] = wssten->sr[0];
		    for (j = 0; j < dim; ++j)
		    {
			wssten->rcrds[i][j] =
			    wssten->coords[j] + i*wssten->dn*wssten->nor[j];
		    }
		}
	        ansr = Left_state_at_node(newcr,cr_or);
	    }

#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("B_node"))
	    {
		PrintWSSten(wssten);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    wssten->ncomp = negative_component(newcr);
	    wssten->pcomp = positive_component(newcr);
	    wssten->w_type = wave_type(tempcr);
	    wssten->hs = Hyper_surf(tempcr);
	    npt_w_speed(wssten,Left_state_at_node(newcr,cr_or),
	                Right_state_at_node(newcr,cr_or),V);
	    ft_assign(Right_state_at_node(newci,ci_or),ansr,sizest);
	}
	else 
	{    /* DIRICHLET_NODE */

	    /* Tangent to ci is oriented parallel to cl/right nor */

	    l_n[0] = r_n[0] = t[0];
	    l_n[1] = r_n[1] = t[1];

	    wssten->nor = l_n;
	    if (ci_or != cl_or) 
	    {
	        l_n[0] *= -1.;
	        l_n[1] *= -1.;
	        wssten->sl[0] = Left_state_at_node(tempci,ci_or);
	        wssten->sr[0] = Right_state_at_node(tempcl,cl_or);
		for (i = 1; i < nsts; ++i)
		{
		    wssten->sr[i] = (Locstate) (wssten->sr_store + i*sizest);
		    for (j = 0; j < dim; ++j)
			wssten->rcrds[i][j] =
			    wssten->coords[j] + i*wssten->dn*wssten->nor[j];
		}
	        switch (wave_type(tempcl))
	        {
	        case DIRICHLET_BOUNDARY:
	            copy_state(stemp2,wssten->sr[0]);
	            copy_state(wssten->sr[0],wssten->sl[0]);
	            Set_params(wssten->sr[0],stemp2);
		    for (i = 1; i < nsts; ++i)
	                hyp_solution(wssten->rcrds[i],
	                             exterior_component(tempci->interface),
	                             Hyper_surf(tempcl),POSITIVE_SIDE,fr,
				     wave,wssten->sr[i],wssten->sr[0]);
	            break;
	        case NEUMANN_BOUNDARY:
		    for (i = 0; i < nsts; ++i)
	                obstacle_state(fr->interf,wssten->sr[i],sizest);
	            break;
	        default:
	            screen("ERROR in snd_B_node_propagate(), "
	                   "Invalid boundary type\n");
	            clean_up(ERROR);
	            break;
	        }
		for (i = 1; i < nsts; ++i)
		{
		    wssten->lcrds[i] = Coords(p[isgn*(i-1)]);
		    wssten->sl[i] = sleft[isgn*(i-1)];
		}
	        ansl = Left_state_at_node(newci,cl_or);
	    }
	    else 
	    {
	        wssten->sr[0] = Left_state_at_node(tempci,ci_or);
	        wssten->sl[0] = Left_state_at_node(tempcl,cl_or);
		for (i = 1; i < nsts; ++i)
		{
		    wssten->sl[i] = (Locstate) (wssten->sl_store + i*sizest);
		    for (j = 0; j < dim; ++j)
			wssten->lcrds[i][j] =
			    wssten->coords[j] - i*wssten->dn*wssten->nor[j];
		}
	        switch (wave_type(tempcl))
	        {
	        case DIRICHLET_BOUNDARY:
	            copy_state(stemp2,wssten->sl[0]);
	            copy_state(wssten->sl[0],wssten->sr[0]);
	            Set_params(wssten->sl[0],stemp2);
		    for (i = 1; i < nsts; ++i)
		    {
	                hyp_solution(wssten->lcrds[i],
	                             exterior_component(tempci->interface),
	                             Hyper_surf(tempcl),NEGATIVE_SIDE,fr,
				     wave,wssten->sl[i],wssten->sl[0]);
		    }
	            break;
	        case NEUMANN_BOUNDARY:
		    for (i = 0; i < nsts; ++i)
	                obstacle_state(fr->interf,wssten->sl[i],sizest);
	            break;
	        default:
	            screen("ERROR in snd_B_node_propagate(), "
	                   "Invalid boundary type\n");
	            clean_up(ERROR);
	            break;
	        }
		for (i = 1; i < nsts; ++i)
		{
		    wssten->rcrds[i] = Coords(p[isgn*(i-1)]);
		    wssten->sr[i] = sleft[isgn*(i-1)];
		}
	        ansl = Right_state_at_node(newcl,cl_or);
	    }

	    wssten->ncomp = negative_component(tempcl);
	    wssten->pcomp = positive_component(tempcl);
	    wssten->w_type = wave_type(tempcl);
	    wssten->hs = Hyper_surf(tempcl);
#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("B_node"))
	    {
		PrintWSSten(wssten);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	    npt_w_speed(wssten,Left_state_at_node(newcl,cl_or),
	                Right_state_at_node(newcl,cl_or),V);

	    wssten->nor = r_n;
	    if (ci_or != cr_or) 
	    {
	        wssten->sr[0] = Right_state_at_node(tempci,ci_or);
	        wssten->sl[0] = Left_state_at_node(tempcr,cr_or);
		for (i = 1; i < nsts; ++i)
		{
	            wssten->sl[i] = (Locstate) (wssten->sl_store + i*sizest);
		    for (j = 0; j < dim; ++j)
			wssten->lcrds[i][j] =
			    wssten->coords[j] - i*wssten->dn*wssten->nor[j];
		}
	        switch (wave_type(tempcr))
	        {
	        case DIRICHLET_BOUNDARY:
	            copy_state(stemp2,wssten->sl[0]);
	            copy_state(wssten->sl[0],wssten->sr[0]);
	            Set_params(wssten->sl[0],stemp2);
		    for (i = 1; i < nsts; ++i)
		    {
	                hyp_solution(wssten->lcrds[i],
	                             exterior_component(tempci->interface),
	                             Hyper_surf(tempcr),NEGATIVE_SIDE,fr,
				     wave,wssten->sl[i],wssten->sl[0]);
		    }
	            break;
	        case NEUMANN_BOUNDARY:
		    for (i = 0; i < nsts; ++i)
	                obstacle_state(fr->interf,wssten->sl[i],sizest);
	            break;
	        default:
	            screen("ERROR in snd_B_node_propagate(), "
	                   "Invalid boundary type\n");
	            clean_up(ERROR);
	            break;
	        }
		for (i = 1; i < nsts; ++i)
		{
		    wssten->rcrds[i] = Coords(p[isgn*(i-1)]);
		    wssten->sr[i] = sright[isgn*(i-1)];
		}
	        ansr = Right_state_at_node(newcr,cr_or);
	    }
	    else 
	    {
	        r_n[0] *= -1.;
	        r_n[1] *= -1.;
	        wssten->sl[0] = Right_state_at_node(tempci,ci_or);
	        wssten->sr[0] = Right_state_at_node(tempcr,cr_or);
		for (i = 1; i < nsts; ++i)
		{
	            wssten->sr[i] = (Locstate) (wssten->sr_store + i*sizest);
		    for (j = 0; j < dim; ++j)
			wssten->rcrds[i][j] =
			    wssten->coords[j] + i*wssten->dn*wssten->nor[j];
		}
	        switch (wave_type(tempcr))
	        {
	        case DIRICHLET_BOUNDARY:
	            copy_state(stemp2,wssten->sr[0]);
	            copy_state(wssten->sr[0],wssten->sl[0]);
	            Set_params(wssten->sr[0],stemp2);
		    for (i = 1; i < nsts; ++i)
		    {
	                hyp_solution(wssten->rcrds[i],
	                             exterior_component(tempci->interface),
	                             Hyper_surf(tempcr),POSITIVE_SIDE,fr,
				     wave,wssten->sr[i],wssten->sr[0]);
		    }
	            break;
	        case NEUMANN_BOUNDARY:
		    for (i = 0; i < nsts; ++i)
	                obstacle_state(fr->interf,wssten->sr[i],sizest);
	            break;
	        default:
	            screen("ERROR in snd_B_node_propagate(), "
	                   "Invalid boundary type\n");
	            clean_up(ERROR);
	            break;
	        }
		for (i = 1; i < nsts; ++i)
		{
	            wssten->sl[i] = sright[isgn*(i-1)];
		    wssten->lcrds[i] = Coords(p[isgn*(i-1)]);
		}
	        ansr = Left_state_at_node(newcr,cr_or);
	    }

#if defined(DEBUG_NODE_PROPAGATE)
	    if (debugging("B_node"))
	    {
		PrintWSSten(wssten);
	    }
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	    wssten->ncomp = negative_component(tempcr);
	    wssten->pcomp = positive_component(tempcr);
	    wssten->w_type = wave_type(tempcr);
	    wssten->hs = Hyper_surf(tempcr);
	    npt_w_speed(wssten,Left_state_at_node(newcr,cr_or),
	                Right_state_at_node(newcr,cr_or),V);

	    if (is_excluded_comp(positive_component(newcr),intfc) == YES)
	        obstacle_state(fr->interf,Right_state_at_node(newcr,cr_or),
			       sizest);
	    if (is_excluded_comp(negative_component(newcr),intfc) == YES)
	        obstacle_state(fr->interf,Left_state_at_node(newcr,cr_or),
			       sizest);
	    if (is_excluded_comp(positive_component(newcl),intfc) == YES)
	        obstacle_state(fr->interf,Right_state_at_node(newcl,cl_or),
			       sizest);
	    if (is_excluded_comp(negative_component(newcl),intfc) == YES)
	        obstacle_state(fr->interf,Left_state_at_node(newcl,cl_or),
			       sizest);
	}

#if defined(DEBUG_NODE_PROPAGATE)
	if (debugging("B_node"))
	{
	    verbose_print_state("Left_state_at_node(newci,ci_or)",
	                Left_state_at_node(newci,ci_or));
	    verbose_print_state("Right_state_at_node(newci,ci_or)",
	                Right_state_at_node(newci,ci_or));
	}
#endif /* defined(DEBUG_NODE_PROPAGATE) */

	debug_print("B_node","Leaving snd_B_node_propagate()\n");
	return YES;
}        /*end snd_B_node_propagate*/

/*
*		snd_contact_wall_node_propagate():
*
*
*	Computes the states at a contact wall node,  including the
*	case of a jet boundary at an inlet.  The method is to use operator
*	splitting based on the stencil of states:
*
*
*                 left side      contact            right side
*                  of contact       ^                 of contact
*                                   |
*	<--wssten_n[2*nsts-1]-->    |
*	                            |
*	         .                  |
*	         .                  |
*	         .                  |
*	                            |
*	                            |
*	<--wssten_n[nsts+1]-->      |
*	                            |
*	                            |
*	                            |
*	<--wssten_n[nsts]-->        |                               
*	------wall/inlet-------------------------wall/inlet---------
*	<--wssten_n[nsts-1]-->      |                               
*	                            |
*	                            |
*	                            |
*	                            |
*	<--wssten_n[nsts-2]-->      |
*	                            |
*	         .                  |
*	         .                  |
*	         .                  |
*	                            |
*	<--wssten_n[0]-->           |
*		               inlet boundary
*			         or artifical
*			         boundary
*
*                 left side      contact            right side
*                  of contact       ^                 of contact
*                                   |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*   / \                  / \        | / \                  / \
*    |                    |         |  |                    |
*  wssten_t[0]  ... wssten_t[nsts-1]|wssten_t[nsts]  ... wssten_t[2*nsts-1]
*    |                    |         |  |                    |
*   \ /                  \ /        | \ /                  \ /
*	------wall/inlet-------------------------wall/inlet---------
*	                            |                               
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*	                            |
*		               inlet boundary
*			         or artifical
*			         boundary
*
*	On alternating time steps,  the solutions of three point wave
*	speed computations are computed on either the horizontal or
*	vertical rows of the state array,  the solution of which are
*	used as data for a final sweep at the node.
*/


/* ARGSUSED */
LOCAL boolean snd_contact_wall_node_propagate(
	Front		*fr,
	Front		*newfr,
	Wave		*wave,
	NODE		*on,
	NODE		*tempn,
	NODE		*nn,
	double		dt)
{
	INTERFACE	*intfc = fr->interf;
	HYPER_SURF_ELEMENT	*ohse;
	COMPONENT	pcompi, ncompi;
	COMPONENT	pcompl, ncompl;
	COMPONENT	pcompr, ncompr;
	ORIENTATION	ci_or;		/* curve orientation */
	ORIENTATION	cl_or, cr_or;
	boolean		sav_include_wall_normal_velocity;
	int		i, j, k, dim = intfc->dim;
	int		isgn;
	double		*ocrds = Coords(on->posn);
	double		*h = fr->rect_grid->h;
	double		length;
	double		l_t[MAXD];		/* tangent to left wall */
	double		r_t[MAXD];		/* tangent to right wall */
	double		mstmp;
	double		dsw, dsi;	 /* tangential mesh on wall and ci */
	double		maxsp, hmin;
	double		t[MAXD];		/* ci tangent */
	double		pjump;
	HYPER_SURF	*ohi;
	HYPER_SURF	*ohl;
	HYPER_SURF	*ohr;
	CURVE		*oci,*nci;	/* incident curve */
	CURVE		*ocl,*ncl;	/* wall left of inc curve */
	CURVE		*ocr,*ncr;	/* wall right of inc curve */
	Locstate	sl, sr;
	Locstate	ansl, ansr;
	int		db_wsp = NO;
	size_t		sizest = fr->sizest;
	static int	nsts = 0;
	static WSSten	**wssten_n = NULL;	/* Stencil normal to contact*/
	static WSSten	**wssten_t = NULL;	/* Stencils tangent to contact*/
	static WSSten	*ansl_t = NULL, *ansr_t = NULL, *ans_n = NULL;
	static double	*nor = NULL;		/* ci normal */
	static double 	*w_t = NULL;		/* average wall tangent */
	static double	*l_n = NULL;		/* cl normal */
	static double	*r_n = NULL;		/* cr normal */
	static POINT	**p = NULL;
	static Locstate	*sleft = NULL, *sright = NULL;
	static Locstate	sts[4][4];
	static Locstate	temp_sts[4][4];
	static Locstate stemp = NULL;

	debug_print("snd_contact","Entered snd_contact_wall_node_propagate()\n");

		/* Allocate storage */

	if (wssten_n == NULL)
	{
	    nsts = nsts_WSSten(fr);
	    ansl_t = AllocDefaultWSSten(fr);
	    ansr_t = AllocDefaultWSSten(fr);
	    ans_n = AllocDefaultWSSten(fr);
	    uni_array(&wssten_n,2*nsts,sizeof(WSSten*));
	    uni_array(&wssten_t,2*nsts,sizeof(WSSten*));
	    for (i = 0; i < 2*nsts; ++i)
	    {
	        wssten_n[i] = AllocDefaultWSSten(fr);
		wssten_t[i] = AllocDefaultWSSten(fr);
	    }

	    for (i = 0; i < nsts; ++i)
	    for (j = 0; j < nsts; ++j)
	    {
		wssten_t[j]->sl[nsts-1-i]  = wssten_n[i]->sl[nsts-1-j];
		wssten_t[j]->lcrds[nsts-1-i] = wssten_n[i]->lcrds[nsts-1-j];
	    }
	    for (i = 0; i < nsts; ++i)
	    for (j = nsts; j < 2*nsts; ++j)
	    {
		wssten_t[j]->sl[nsts-1-i]  = wssten_n[i]->sr[j-nsts];
		wssten_t[j]->lcrds[nsts-1-i] = wssten_n[i]->rcrds[j-nsts];
	    }
	    for (i = nsts; i < 2*nsts; ++i)
	    for (j = 0; j < nsts; ++j)
	    {
		wssten_t[j]->sr[i-nsts]  = wssten_n[i]->sl[nsts-1-j];
		wssten_t[j]->rcrds[i-nsts] = wssten_n[i]->lcrds[nsts-1-j];
	    }
	    for (i = nsts; i < 2*nsts; ++i)
	    for (j = nsts; j < 2*nsts; ++j)
	    {
		wssten_t[j]->sr[i-nsts]  = wssten_n[i]->sr[j-nsts];
		wssten_t[j]->rcrds[i-nsts] = wssten_n[i]->rcrds[j-nsts];
	    }
	    w_t = wssten_n[nsts]->nor;
	    for (i = 0; i < nsts; ++i)
		wssten_n[i]->nor = w_t;
	    nor = wssten_n[nsts+1]->nor;
	    for (i = nsts+2; i < 2*nsts; ++i)
		wssten_n[i]->nor = nor;
	    l_n = wssten_t[0]->nor;
	    for (i = 1; i < nsts; ++i)
		wssten_t[i]->nor = l_n;
	    r_n = wssten_t[nsts]->nor;
	    for (i = nsts+1; i < 2*nsts; ++i)
		wssten_t[i]->nor = r_n;
	    for (i = 0; i < 2*nsts; ++i)
	    {
		wssten_t[i]->coords = wssten_t[i]->rcrds[0];
		wssten_n[i]->coords = wssten_n[i]->rcrds[0];
	    }

	    uni_array(&sleft,2*nsts-1,sizeof(Locstate));
	    sleft += nsts-1;
	    uni_array(&sright,2*nsts-1,sizeof(Locstate));
	    sright += nsts-1;

	    uni_array(&p,2*nsts-1,sizeof(POINT*));
	    p += nsts-1;
	    for (i = -nsts+1; i < nsts; ++i)
		p[i] = Static_point(fr->interf);

	    ansl_t->nor = l_n;
	    ansr_t->nor = r_n;
	    for (i = 0; i < nsts; ++i)
	    {
		ansl_t->lcrds[i] = wssten_t[nsts-1]->lcrds[i];
		ansl_t->rcrds[i] = wssten_t[nsts-1]->rcrds[i];
		ansr_t->lcrds[i] = wssten_t[nsts]->lcrds[i];
		ansr_t->rcrds[i] = wssten_t[nsts]->rcrds[i];
	    }
	    ans_n->nor = w_t;
	    for (i = 0; i < nsts; ++i)
	    {
		ans_n->lcrds[i] = wssten_n[nsts]->lcrds[i];
		ans_n->rcrds[i] = wssten_n[nsts]->rcrds[i];
	    }
	}
	ansl_t->coords = Coords(nn->posn);
	ansr_t->coords = Coords(nn->posn);
	ans_n->coords = Coords(nn->posn);

	for (i = 0; i < 2*nsts; ++i)
	{
	    wssten_n[i]->front = fr;
	    wssten_n[i]->wave = wave;
	    wssten_t[i]->front = fr;
	    wssten_t[i]->wave = wave;
	}
	if (stemp == NULL)
	{
	    for (i = 0; i < 4; ++i)
	    for (j = 0; j < 4; ++j)
	    {
	    	alloc_state(intfc,&sts[i][j],sizest);
	    	alloc_state(intfc,&temp_sts[i][j],sizest);
	    }
	    alloc_state(intfc,&stemp,sizest);
	}
	for (i = 0; i < 4; ++i)
	for (j = 0; j < 4; ++j)
	{
	    clear_state(intfc,sts[i][j],sizest);
	    clear_state(intfc,temp_sts[i][j],sizest);
	}

		/* Identify curves */

	find_curve_with_status(on,&oci,&ci_or,INCIDENT);
	ohi = Hyper_surf(oci);
	pcompi = positive_component(ohi);
	ncompi = negative_component(ohi);
	for (i = nsts; i < 2*nsts; ++i)
	{
	    wssten_n[i]->ncomp = ncompi;
	    wssten_n[i]->pcomp = pcompi;
	    wssten_n[i]->w_type = wave_type(ohi);
	    wssten_n[i]->hs = ohi;
	}
	find_curve_with_status(nn,&nci,&ci_or,INCIDENT);
	if (ci_or == NEGATIVE_ORIENTATION)
	{
	    ocl = adjacent_curve(oci,ci_or,CLOCKWISE,&cl_or);
	    ncl = adjacent_curve(nci,ci_or,CLOCKWISE,&cl_or);
	    ocr = adjacent_curve(oci,ci_or,COUNTER_CLOCK,&cr_or);
	    ncr = adjacent_curve(nci,ci_or,COUNTER_CLOCK,&cr_or);
	}
	else
	{
	    ocl = adjacent_curve(oci,ci_or,COUNTER_CLOCK,&cl_or);
	    ncl = adjacent_curve(nci,ci_or,COUNTER_CLOCK,&cl_or);
	    ocr = adjacent_curve(oci,ci_or,CLOCKWISE,&cr_or);
	    ncr = adjacent_curve(nci,ci_or,CLOCKWISE,&cr_or);
	}
	ohr = Hyper_surf(ocr);
	ohl = Hyper_surf(ocl);

	normal(on->posn,Hyper_surf_element(Bond_at_node(oci,ci_or)),
	       Hyper_surf(oci),nor,fr);

	if (debugging("snd_contact"))
	{
	    (void) printf("dt = %g\n",dt);
	    print_general_vector("nor = ",nor,dim,"\n");
	    print_orientation("ci_or = ",ci_or,"\n");
	    print_orientation("cl_or = ",cl_or,"\n");
	    print_orientation("cr_or = ",cr_or,"\n");
	    (void) printf("oci\n");    print_curve(oci);
	    (void) printf("nci\n");    print_curve(nci);
	    (void) printf("ocl\n");    print_curve(ocl);
	    (void) printf("ncl\n");    print_curve(ncl);
	    (void) printf("ocr\n");    print_curve(ocr);
	    (void) printf("ncr\n");    print_curve(ncr);
	}

	    /* orient tangent to incident curve out from node */

	if (ci_or == POSITIVE_ORIENTATION)
	{
	    t[0]  = -nor[1];
	    t[1]  =  nor[0];
	}
	else
	{
	    t[0]  =  nor[1];
	    t[1]  = -nor[0];
	}
	dsi = grid_size_in_direction(t,h,dim);
	hmin = min(h[0],h[1]);

	/* find average tangent to wall, oriented as is (nor[0],nor[1]) */

	normal(on->posn,Hyper_surf_element(Bond_at_node(ocl,cl_or)),
	       Hyper_surf(ocl),l_n,fr);
	if (scalar_product(t,l_n,dim) < 0.0)
	{
	    for (i = 0; i < dim; ++i)
	        l_n[i] = -l_n[i];
	    ncompl = positive_component(ohl);
	    pcompl = negative_component(ohl);
	}
	else
	{
	    pcompl = positive_component(ohl);
	    ncompl = negative_component(ohl);
	}
	for (i = 0; i < nsts; ++i)
	{
	    wssten_t[i]->pcomp = pcompl;
	    wssten_t[i]->ncomp = ncompl;
	    wssten_t[i]->w_type = wave_type(ohl);
	    wssten_t[i]->hs = ohl;
	}
	normal(on->posn,Hyper_surf_element(Bond_at_node(ocr,cr_or)),
	       Hyper_surf(ocr),r_n,fr);
	if (scalar_product(t,r_n,dim) < 0.0)
	{
	    for (i = 0; i < dim; ++i)
	        r_n[i] = -r_n[i];
	    ncompr = positive_component(ohr);
	    pcompr = negative_component(ohr);
	}
	else
	{
	    pcompr = positive_component(ohr);
	    ncompr = negative_component(ohr);
	}
	for (i = nsts; i < 2*nsts; ++i)
	{
	    wssten_t[i]->pcomp = pcompr;
	    wssten_t[i]->ncomp = ncompr;
	    wssten_t[i]->w_type = wave_type(ohr);
	    wssten_t[i]->hs = ohr;
	}
	if (debugging("snd_contact"))
	{
	    print_general_vector("l_n = ",l_n,dim,", ");
	    print_general_vector("r_n = ",r_n,dim,"\n");
	}

	r_t[0] = -r_n[1];
	r_t[1] =  r_n[0];
	if (scalar_product(r_t,nor,dim) < 0.0)
	{
	    for (i = 0; i < dim; ++i)
	        r_t[i] = -r_t[i];
	}
	l_t[0] = -l_n[1];
	l_t[1] =  l_n[0];
	if (scalar_product(l_t,nor,dim) < 0.0)
	{
	    for (i = 0; i < dim; ++i)
	        l_t[i] = -l_t[i];
	}
	w_t[0] = r_t[0] + l_t[0];
	w_t[1] = r_t[1] + l_t[1];
	length = mag_vector(w_t,dim);
	w_t[0] /= length; w_t[1] /= length;
	dsw = grid_size_in_direction(w_t,h,dim);

	    /* divide time interval by sin(contact_wall_angle) */

	/*dt /= scalar_product(w_t,nor,dim);  IS THIS CORRECT? */

	if (debugging("snd_contact"))
	{
	    (void) printf("new dt = %g\n",dt);
	    print_general_vector("t = ",t,dim,", ");
	    print_general_vector("w_t = ",w_t,dim,"\n");
	    (void) printf("sin(contact_wall_angle) = %g\n",
	    	          scalar_product(w_t,nor,dim));
	}

	    /* find states at node */

	ohse = Hyper_surf_element(Bond_at_node(oci,ci_or));
	slsr(on->posn,ohse,ohi,&sl,&sr);
	set_state(wssten_n[nsts]->sl[0],GAS_STATE,sl);
	set_state(wssten_n[nsts]->sr[0],GAS_STATE,sr);
	for (i = 0; i < dim; ++i)
	{
	    wssten_n[nsts]->lcrds[0][i] = ocrds[i];
	    wssten_n[nsts]->rcrds[0][i] = ocrds[i];
	}

	    /* find states on oci */

	isgn = (ci_or == POSITIVE_ORIENTATION) ? 1 : -1;
	for (i = 1; i < nsts; ++i)
	{
	    sleft[isgn*(i-1)] = wssten_t[nsts-1]->sr[i];
	    sright[isgn*(i-1)] = wssten_t[nsts]->sr[i];
	}
	states_at_distance_along_curve(on->posn,Bond_at_node(oci,ci_or),
				       oci,ci_or,dsi,nsts-1,sleft,sright,
				       NULL,NULL,NULL,p,fr);
	for (i = 1; i < nsts; ++i)
	{
	    for (k = 0; k < dim; ++k)
	    {
	        wssten_t[nsts-1]->rcrds[i][k] = Coords(p[isgn*(i-1)])[k];
	        wssten_t[nsts]->rcrds[i][k] = Coords(p[isgn*(i-1)])[k];
	    }
	}

	    /* find wall states */

	isgn = (cl_or == POSITIVE_ORIENTATION) ? 1 : -1;
	if (cl_or == ci_or)
	{
	    set_state(wssten_n[nsts-1]->sl[0],GAS_STATE,
		      Left_state_at_node(ocl,cl_or));
	    for (i = 1; i < nsts; ++i)
	    {
	        sleft[isgn*(i-1)] = wssten_n[nsts-1]->sl[i];
	        sright[isgn*(i-1)] = wssten_n[nsts]->sl[i];
	    }
	    states_at_distance_along_curve(on->posn,Bond_at_node(ocl,cl_or),
				           ocl,cl_or,dsw,nsts-1,sleft,sright,
				           NULL,NULL,NULL,p,fr);
	}
	else
	{
	    set_state(wssten_n[nsts-1]->sl[0],GAS_STATE,
		      Right_state_at_node(ocl,cl_or));
	    for (i = 1; i < nsts; ++i)
	    {
	        sleft[isgn*(i-1)] = wssten_n[nsts]->sl[i];
	        sright[isgn*(i-1)] = wssten_n[nsts-1]->sl[i];
	    }
	    states_at_distance_along_curve(on->posn,Bond_at_node(ocl,cl_or),
				           ocl,cl_or,dsw,nsts-1,sleft,sright,
				           NULL,NULL,NULL,p,fr);
	}
	for (i = 1; i < nsts; ++i)
	{
	    for (k = 0; k < dim; ++k)
	    {
	        wssten_n[nsts-1]->lcrds[i][k] = Coords(p[isgn*(i-1)])[k];
	        wssten_n[nsts]->lcrds[i][k] = Coords(p[isgn*(i-1)])[k];
	    }
	}

	isgn = (cr_or == POSITIVE_ORIENTATION) ? 1 : -1;
	if (cr_or == ci_or)
	{
	    set_state(wssten_n[nsts-1]->sr[0],GAS_STATE,
		      Right_state_at_node(ocr,cr_or));
	    for (i = 1; i < nsts; ++i)
	    {
	        sleft[isgn*(i-1)] = wssten_n[nsts]->sr[i];
	        sright[isgn*(i-1)] = wssten_n[nsts-1]->sr[i];
	    }
	    states_at_distance_along_curve(on->posn,Bond_at_node(ocr,cr_or),
				           ocr,cr_or,dsw,nsts-1,sleft,sright,
					   NULL,NULL,NULL,p,fr);
	}
	else
	{
	    set_state(wssten_n[nsts-1]->sr[0],GAS_STATE,
		      Left_state_at_node(ocr,cr_or));
	    for (i = 1; i < nsts; ++i)
	    {
	        sleft[isgn*(i-1)] = wssten_n[nsts-1]->sr[i];
	        sright[isgn*(i-1)] = wssten_n[nsts]->sr[i];
	    }
	    states_at_distance_along_curve(on->posn,Bond_at_node(ocr,cr_or),
				           ocr,cr_or,dsw,nsts-1,sleft,sright,
					   NULL,NULL,NULL,p,fr);
	}
	for (i = 1; i < nsts; ++i)
	{
	    for (k = 0; k < dim; ++k)
	    {
	        wssten_n[nsts-1]->rcrds[i][k] = Coords(p[isgn*(i-1)])[k];
	        wssten_n[nsts]->rcrds[i][k] = Coords(p[isgn*(i-1)])[k];
	    }
	}

	    /* find states at interior points */

	for (i = nsts+1; i < 2*nsts; ++i)
	for (j = 1; j < nsts; ++j)
	{
	    for (k = 0; k < dim; ++k)
	    {
	        wssten_n[i]->lcrds[j][k] = ocrds[k] +
					   (i-nsts)*dsi*t[k] - j*dsw*w_t[k];
	        wssten_n[i]->rcrds[j][k] = ocrds[k] +
					   (i-nsts)*dsi*t[k] + j*dsw*w_t[k];
	    }
	    hyp_solution(wssten_n[i]->lcrds[j],wssten_n[i]->ncomp,
			 wssten_n[i]->hs,NEGATIVE_SIDE,fr,wave,
			 wssten_n[i]->sl[j],wssten_n[nsts]->sl[0]);
	    hyp_solution(wssten_n[i]->rcrds[j],wssten_n[i]->pcomp,
			 wssten_n[i]->hs,POSITIVE_SIDE,fr,wave,
			 wssten_n[i]->sr[j],wssten_n[nsts]->sr[0]);
	}

	/* Find exterior off front states */

	for (i = 0; i < nsts; ++i)
	for (j = 0; j < nsts; ++j)
	{
	    for (k = 0; k < dim; ++k)
	    {
	        wssten_n[i]->lcrds[j][k] = ocrds[k] -
					   (nsts-1-i)*dsi*l_n[k] - j*dsw*w_t[k];
	        wssten_n[i]->rcrds[j][k] = ocrds[k] -
					   (nsts-1-i)*dsi*r_n[k] + j*dsw*w_t[k];
	    }
	}

	if (wave_type(ocl) == NEUMANN_BOUNDARY)
	{
	    for (i = 0; i < nsts-1; ++i)
	    for (j = 0; j < nsts; ++j)
	    {
	        obstacle_state(intfc,wssten_n[i]->sl[j],sizest);
	    }
	}
	else
	{
	    SIDE side;

	    side = (cl_or == ci_or) ? NEGATIVE_SIDE : POSITIVE_SIDE;
	    for (i = 1; i < nsts; ++i)
	    for (j = 0; j < nsts; ++j)
	    {
	        hyp_solution(wssten_t[j]->lcrds[i],wssten_t[j]->ncomp,
			     wssten_t[j]->hs,side,fr,wave,
			     wssten_t[j]->sl[i],wssten_n[nsts-1]->sl[0]);
	    }
	}

	if (wave_type(ocr) == NEUMANN_BOUNDARY)
	{
	    for (i = 0; i < nsts-1; ++i)
	    for (j = 0; j < nsts; ++j)
	    {
	        obstacle_state(intfc,wssten_n[i]->sr[j],sizest);
	    }
	}
	else
	{
	    SIDE       side;

	    side = (cr_or == ci_or) ? POSITIVE_SIDE : NEGATIVE_SIDE;
	    for (i = 1; i < nsts; ++i)
	    for (j = nsts; j < 2*nsts; ++j)
	    {
	        hyp_solution(wssten_t[j]->lcrds[i],wssten_t[j]->ncomp,
			     wssten_t[j]->hs,side,fr,wave,
			     wssten_t[j]->sl[i],wssten_n[nsts-1]->sr[0]);
	    }
	}

	maxsp = -HUGE_VAL;
	for (i = 0; i < 2*nsts; ++i)
	for (j = 0; j < nsts; ++j)
	{
	    mstmp = maximum_wave_speed(wave,wssten_n[i]->sl[j]);
	    maxsp = max(mstmp,maxsp);
	    mstmp = maximum_wave_speed(wave,wssten_n[i]->sr[j]);
	    maxsp = max(mstmp,maxsp);
	}

	if (debugging("snd_contact"))
	{
	    (void) printf("maxsp = %g\n",maxsp);
	    for (i = 0; i < 2*nsts; ++i)
	    {
		(void) printf("wssten_n[%d]\n",i);
		PrintWSSten(wssten_n[i]);
	    }
	    for (j = 0; j < 2*nsts; ++j)
	    {
		(void) printf("wssten_t[%d]\n",j);
		PrintWSSten(wssten_t[j]);
	    }
	    db_wsp = debugging("w_speed");
	    if (!db_wsp)
		add_to_debug("w_speed");
	}

	if (dt*maxsp > Time_step_factor(fr) * hmin)
	{
	    /* Violation of CFL,  TODO Resolve this issue */
	    (void) printf("WARNING in snd_contact_wall_node_propagate(), ");
	    (void) printf("CFL violation\n");
	    (void) printf("dt = %g, maxsp = %g, ",dt,maxsp);
	    (void) printf("Time_step_factor(fr) = %g, hmin = %g\n",
	    	          Time_step_factor(fr),hmin);
	    (void) printf("dt*maxsp/(Time_step_factor(fr)*hmin) = %g\n",
			  dt*maxsp/(Time_step_factor(fr)*hmin));
	    (void) printf("dt*maxsp/hmin = %g\n",dt*maxsp/hmin);
	    if (dt*maxsp >  hmin)
	    {
	    	(void) printf("second node propagate aborted\n");
	    	return YES;
	    }
	}

	for (i = 0; i < 2*nsts; ++i)
	{
	    wssten_n[i]->dn = dsw;
	    wssten_t[i]->dn = dsi;
	}
	for (i = 0; i < nsts; ++i)
	{
	    wssten_n[i]->hs = NULL;
	    wssten_n[i]->w_type = NEUMANN_BOUNDARY;
	}
	for (i = 0; i < 2*nsts; ++i)
	{
	    wssten_n[i]->dt = wssten_t[i]->dt = dt;
	}

	    /* alternate the order of the sweeps */


	sav_include_wall_normal_velocity = include_wall_normal_velocity(fr);
	if (fr->step%2) /* normal to contact, then normal to walls */
	{
	    if (debugging("snd_contact"))
	    {
		(void) printf("Odd numbered time step, sweeping normal to "
			      "contact, then normal to walls\n");
	    }
	    
	    ansl_t->dt = ansr_t->dt = dt;
	    ansl_t->pjump = ansr_t->pjump = 0.0;
	    ansl_t->front = ansr_t->front = fr;
	    ansl_t->wave = ansr_t->wave = wave;
	    ansl_t->dn = ansr_t->dn = dsi;

	    ansl_t->w_type = wssten_t[nsts-1]->w_type;
	    ansl_t->pcomp = wssten_t[nsts-1]->pcomp;
	    ansl_t->ncomp = wssten_t[nsts-1]->ncomp;
	    ansl_t->hs = wssten_t[nsts-1]->hs;

	    ansr_t->w_type = wssten_t[nsts]->w_type;
	    ansr_t->pcomp = wssten_t[nsts]->pcomp;
	    ansr_t->ncomp = wssten_t[nsts]->ncomp;
	    ansr_t->hs = wssten_t[nsts]->hs;

	    pjump = set_pjump_at_wave(on->posn,ohse,ohi,fr,nor);
	    /* sweep to get temp states at node */
	    for (i = 0; i < nsts; ++i)
	    {
	        if (debugging("snd_contact"))
		    (void) printf("Calling npt_w_speed on wssten_n[%d]\n",i);
		if (is_obstacle_state(wssten_n[i]->sl[0]) &&
		    is_obstacle_state(wssten_n[i]->sr[0]))
		{
		    obstacle_state(intfc,ansl_t->sl[nsts-1-i],sizest);
		    obstacle_state(intfc,ansr_t->sl[nsts-1-i],sizest);
		}
		else
		{
		    wssten_n[i]->pjump = 0.0;
		    npt_w_speed(wssten_n[i],ansl_t->sl[nsts-1-i],
			        ansr_t->sl[nsts-1-i],NULL);
		}
	    }
	    for (i = nsts; i < 2*nsts; ++i)
	    {
	        if (debugging("snd_contact"))
		    (void) printf("Calling npt_w_speed on wssten_n[%d]\n",i);
		if (is_obstacle_state(wssten_n[i]->sl[0]) &&
		    is_obstacle_state(wssten_n[i]->sr[0]))
		{
		    obstacle_state(intfc,ansl_t->sr[i-nsts],sizest);
		    obstacle_state(intfc,ansr_t->sr[i-nsts],sizest);
		}
		else
		{
		    wssten_n[i]->pjump = pjump;
		    npt_w_speed(wssten_n[i],ansl_t->sr[i-nsts],
			        ansr_t->sr[i-nsts],NULL);
		}
	    }

	    /* final sweep to get states at node */

	    if (wave_type(ocl) == NEUMANN_BOUNDARY)
		include_wall_normal_velocity(fr) = YES;
	    if (debugging("snd_contact"))
	    {
		(void) printf("Calling npt_w_speed on ansl_t\n");
		PrintWSSten(ansl_t);
	    }
	    npt_w_speed(ansl_t,stemp,Left_state_at_node(nci,ci_or),NULL);
	    include_wall_normal_velocity(fr) = sav_include_wall_normal_velocity;

	    if (wave_type(ocr) == NEUMANN_BOUNDARY)
	        include_wall_normal_velocity(fr) = YES;
	    if (debugging("snd_contact"))
	    {
		(void) printf("Calling npt_w_speed on ansr_t\n");
		PrintWSSten(ansr_t);
	    }
	    npt_w_speed(ansr_t,stemp,Right_state_at_node(nci,ci_or),NULL);
	    include_wall_normal_velocity(fr) = sav_include_wall_normal_velocity;

	    /* sweep to get final states at node */

	    ansl =  Left_state_at_node(nci,ci_or);
	    ansr = Right_state_at_node(nci,ci_or);
	    if (debugging("snd_contact"))
		(void) printf("Calling w_speed for final states\n");
	    w_speed(Coords(tempn->posn),ansl,ansr,ansl,ansr,NULL,0.0,
	    	    w_t,wave_type(oci),fr);
	}
	else	/* normal to walls, then normal to contact */
	{
	    if (debugging("snd_contact"))
	    {
		(void) printf("Even numbered time step, sweeping normal to "
			      "walls, then normal to contact\n");
	    }
	    
	    ans_n->dt = dt;
	    ans_n->pjump = 0.0;
	    ans_n->front = fr;
	    ans_n->wave = wave;
	    ans_n->dn = dsw;

	    ans_n->w_type = wssten_n[nsts]->w_type;
	    ans_n->pcomp = wssten_n[nsts]->pcomp;
	    ans_n->ncomp = wssten_n[nsts]->ncomp;
	    ans_n->hs = wssten_n[nsts]->hs;
	    include_wall_normal_velocity(fr) = YES;
	    for (i = 0; i < nsts; ++i)
	    {
	        if (debugging("snd_contact"))
		    (void) printf("Calling npt_w_speed on wssten_t[%d]\n",i);
		wssten_t[i]->pjump = 0.0;
		npt_w_speed(wssten_t[i],NULL,ans_n->sl[nsts-1-i],NULL);
	    }
	    for (i = nsts; i < 2*nsts; ++i)
	    {
	        if (debugging("snd_contact"))
		    (void) printf("Calling npt_w_speed on wssten_t[%d]\n",i);
		wssten_t[i]->pjump = 0.0;
		npt_w_speed(wssten_t[i],NULL,ans_n->sr[i-nsts],NULL);
	    }
	    include_wall_normal_velocity(fr) = sav_include_wall_normal_velocity;

	    /* sweep to get final states at node */

	    if (debugging("snd_contact"))
	    {
		(void) printf("Calling npt_w_speed on ans_n\n");
		PrintWSSten(ans_n);
	    }
	    npt_w_speed(ans_n,Left_state_at_node(nci,ci_or),
			Right_state_at_node(nci,ci_or),NULL);
	}
	include_wall_normal_velocity(fr) = sav_include_wall_normal_velocity;
	if (debugging("snd_contact"))
	    if (!db_wsp) remove_from_debug("w_speed");

	    /* average left and right so that */
	    /* Rankine-Hugoniot and boundary conditions hold */

	ansl = Left_state_at_node(nci,ci_or);
	if (wave_type(ocl) == NEUMANN_BOUNDARY)
	{
	    if (ci_or != cl_or)
	    	w_speed(Coords(tempn->posn),ansl,return_obst_state(),
	    		ansl,stemp,NULL,0.0,l_n,wave_type(ocl),fr);
	    else
	    	w_speed(Coords(tempn->posn),return_obst_state(),ansl,
	    		stemp,ansl,NULL,0.0,l_n,wave_type(ocl),fr);
	}
	ansr = Right_state_at_node(nci,ci_or);
	if (wave_type(ocr) == NEUMANN_BOUNDARY)
	{
	    if (ci_or != cr_or)
	    	w_speed(Coords(tempn->posn),return_obst_state(),ansr,
	    		stemp,ansr,NULL,0.0,r_n,wave_type(ocr),fr);
	    else
	    	w_speed(Coords(tempn->posn),ansr,return_obst_state(),
	    		ansr,stemp,NULL,0.0,r_n,wave_type(ocr),fr);
	}
	w_speed(Coords(tempn->posn),ansl,ansr,ansl,ansr,NULL,0.0,
		w_t,wave_type(oci),fr);
	if ((wave_type(ocl) == NEUMANN_BOUNDARY) && no_slip(Hyper_surf(ocl)))
	{
	    double alpha = 1.0 - adherence_coeff(Hyper_surf(ocl));
	    normal(Node_of(ocl,cl_or)->posn,
	    		Hyper_surf_element(Bond_at_node(ocl,cl_or)),
			Hyper_surf(ocl),l_n,fr);
	    alpha_state_velocity(alpha,ansl,dim);
	    zero_normal_velocity(ansl,l_n,dim);
	}
	if ((wave_type(ocr) == NEUMANN_BOUNDARY) && no_slip(Hyper_surf(ocr)))
	{
	    double alpha = 1.0 - adherence_coeff(Hyper_surf(ocr));
	    normal(Node_of(ocr,cr_or)->posn,
	    		Hyper_surf_element(Bond_at_node(ocr,cr_or)),
			Hyper_surf(ocr),r_n,fr);
	    alpha_state_velocity(alpha,ansl,dim);
	    zero_normal_velocity(ansr,r_n,dim);
	}

	if (debugging("snd_contact"))
	{
	    (void) printf("After averaging:\n");
	    verbose_print_state("ansl",ansl);
	    verbose_print_state("ansr",ansr);
	}

	    /* copy answer into node states on walls */

	ansl = (ci_or != cl_or) ? Left_state_at_node(ncl,cl_or) :
	    			  Right_state_at_node(ncl,cl_or);

	ansr = (ci_or != cr_or) ? Right_state_at_node(ncr,cr_or) :
	    			  Left_state_at_node(ncr,cr_or);

	ft_assign(Left_state_at_node(nci,ci_or),ansl,sizest);
	ft_assign(Right_state_at_node(nci,ci_or),ansr,sizest);
	if (wave_type(ocl) == NEUMANN_BOUNDARY)
	{
	    if (ci_or != cl_or)
	    	obstacle_state(intfc,Right_state_at_node(ncl,cl_or),sizest);
	    else
	    	obstacle_state(intfc,Left_state_at_node(ncl,cl_or),sizest);
	}
	else
	{
	    if (ci_or != cl_or)
	    	ft_assign(Right_state_at_node(ncl,cl_or),
		       Left_state_at_node(ncl,cl_or),sizest);
	    else
	    	ft_assign(Left_state_at_node(ncl,cl_or),
	    	       Right_state_at_node(ncl,cl_or),sizest);
	}
	if (wave_type(ocr) == NEUMANN_BOUNDARY)
	{
	    if (ci_or != cr_or)
	    	obstacle_state(intfc,Left_state_at_node(ncr,cr_or),sizest);
	    else
	    	obstacle_state(intfc,Right_state_at_node(ncr,cr_or),sizest);
	}
	else
	{
	    if (ci_or != cr_or)
	    	ft_assign(Left_state_at_node(ncr,cr_or),
	    	       Right_state_at_node(ncr,cr_or),sizest);
	    else
	    	ft_assign(Right_state_at_node(ncr,cr_or),
	    	       Left_state_at_node(ncr,cr_or),sizest);
	}

	debug_print("snd_contact","Left snd_contact_wall_node_propagate()\n");
	return YES;
}		/*end snd_contact_wall_node_propagate*/


/*
*			update_wall_states():
*
*	Copies states from a wall node to those wall points crossed by the
*	node during the time step.  Resets the t_pt_propagated flats so
*	that the tangential solver doesn't update these states.
*
*/

LOCAL void update_wall_states(
	Front		*front,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn)
{
	ORIENTATION	c1_orient, c2_orient, cinc_orient;
	CURVE		*newc1, *tempc1;
	CURVE		*newc2, *tempc2;
	CURVE		*newcinc, *tempcinc;

	find_curve_with_status(tempn,&tempcinc,&cinc_orient,INCIDENT);
	find_curve_with_status(newn,&newcinc,&cinc_orient,INCIDENT);

	tempc1 = adjacent_curve(tempcinc,cinc_orient,CLOCKWISE,&c1_orient);
	newc1 = adjacent_curve(newcinc,cinc_orient,CLOCKWISE,&c1_orient);
	tempc2 = adjacent_curve(tempcinc,cinc_orient,COUNTER_CLOCK,&c2_orient);
	newc2 = adjacent_curve(newcinc,cinc_orient,COUNTER_CLOCK,&c2_orient);

	update_points_crossed_by_node(tempc1, newc1, c1_orient, oldn, front);
	update_points_crossed_by_node(tempc2, newc2, c2_orient, oldn, front);

	return;
}		/*end update_wall_states*/

LOCAL void update_points_crossed_by_node(
	CURVE		*tempc,
	CURVE		*newc,
	ORIENTATION	orient,
	NODE		*oldn,
	Front		*front)
{
	double		newdiff[MAXD], olddiff[MAXD];
	size_t		sizest;
	int		dim,i;
	BOND		*tempb, *newb;
	Locstate	r_node_state, l_node_state;
	NODE		*newn;
	POINT		*p;
	
	sizest = front->sizest;
	dim = front->interf->dim;
	newn = Node_of(newc, orient);
	tempb = Bond_at_node(tempc, orient);
	newb = Bond_at_node(newc, orient);

	r_node_state = Right_state_at_node(newc, orient);
	l_node_state = Left_state_at_node(newc, orient);

	for ( ; tempb!=NULL; tempb = Following_bond(tempb,orient),
				newb = Following_bond(newb,orient))
	{
		p = Point_of_bond(newb, Opposite_orient(orient));

		for (i = 0; i < dim; ++i)
		{
			newdiff[i] = Coords(newn->posn)[i] - Coords(p)[i];
			olddiff[i] = Coords(oldn->posn)[i] - Coords(p)[i];
		}

		if (scalar_product(newdiff,olddiff,dim)>0.0)
				break;
			
		ft_assign(left_state(p), l_node_state, sizest);
		ft_assign(right_state(p), r_node_state, sizest);
		t_pt_propagated(p) = YES;
	}

	return;
}		/*end update_states_crossed_by_node*/



/* ARGSUSED */
LOCAL boolean snd_B_reflect_node_propagate(
	Front		*fr,
	Front		*newfr,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	/* TODO: CODE NEEDED */
	return YES;
}		/*end snd_B_reflect_node_propagate*/


/* ARGSUSED */
LOCAL boolean snd_Mach_node_propagate(
	Front		*fr,
	Front		*newfr,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	/* TODO: CODE NEEDED */
	return YES;
}		/*end snd_Mach_node_propagate*/


/* ARGSUSED */
LOCAL boolean snd_attached_b_node_propagate(
	Front		*fr,
	Front		*newfr,
	Wave		*wave,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	ORIENTATION	cinc_orient;
	CURVE		*tempcinc;
	boolean		status;

	debug_print("attached_b_node","Entered snd_attached_b_node_propagate()\n");

	/* TODO: CODE NEEDED FOR SHOCK */

	find_curve_with_status(tempn,&tempcinc,&cinc_orient,INCIDENT);
	if (is_scalar_wave(wave_type(tempcinc)) &&
	    (node_type(tempn) == ATTACHED_B_NODE))
	    status = snd_B_node_propagate(fr,newfr,wave,oldn,tempn,newn,dt);
	else
	    status = YES;

	debug_print("attached_b_node","Left snd_attached_b_node_propagate()\n");
	return status;
}		/*end snd_attached_b_node_propagate*/

#if defined(FULL_PHYSICS)

/*
*			snd_cross_node_propagate():
*
*	Handles the second node propagate for the cross node.
*	Currently only the case of a degenerate cross node is
*	handled.  This is the case of no reflected waves (either
*	tracked or untracked).  In this case the node can be 
*	essentially regarded as an interior point on the
*	conceptual curve formed by joining the two incident curves.
*/

/*ARGSUSED*/
LOCAL boolean snd_cross_node_propagate(
	Front		*fr,
	Front		*newfr,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	O_CURVE		Tempc0, Tempc4, Newc0, Newc4;
	CURVE		**c;
	POINT		*p0, *p4;
	BOND		*b0, *b4;
	double		t0[MAXD], t4[MAXD];
	double		cos04, sin04, ang04;
	double		u[MAXD], cr;
	double		ds, tgnt[MAXD];
	int		isgn;
	ANGLE_DIRECTION	i0_to_i4_dir;
	int		num_curves;
	int		i, dim = fr->rect_grid->dim;
	boolean		is_plus_orientation, is_refl_rarefaction;
	Locstate	st01, st04, ansl0, ansr0, ansl4, ansr4;
	static RP_DATA	*RP = NULL;
	static int	nrad = 0;
	static Tan_stencil *sten = NULL;
	static const double ERR_ANG = 0.45; /*TOLERANCE*/

	debug_print("snd_cross_node","Entered snd_cross_node_propagate()\n");
	if (RP == NULL) 
	{
	    nrad = fr->npts_tan_sten/2;
	    sten = alloc_tan_stencil(fr,nrad);
	    RP = allocate_RP_DATA_structure(fr->sizest,NO,GAS_STATE);
	}

	num_curves = 0;
	for (c = tempn->in_curves; c && *c; ++c)
	    ++num_curves;
	for (c = tempn->out_curves; c && *c; ++c)
	    ++num_curves;

	if (num_curves != 2)
	{
	    debug_print("snd_cross_node","Leaving snd_cross_node_propagate()\n");
	    return YES;
	}

	identify_curves_with_status(tempn,&Tempc0,&Tempc4,INCIDENT);

	if (debugging("snd_cross_node")) 
	{
	    print_cross_node(newn,&Tempc0,NULL,NULL,NULL,&Tempc4);
	}
	Check_return(
	    find_correspond_of_oriented_curve(&Tempc0,&Newc0,newn,fr,
					      newn->interface),
	    snd_cross_node_propagate)
	Check_return(
	    find_correspond_of_oriented_curve(&Tempc4,&Newc4,newn,fr,
					      newn->interface),
	    snd_cross_node_propagate)
	if ((Tempc0.orient == POSITIVE_ORIENTATION && 
	         wave_type(Tempc0.curve) == FORWARD_SHOCK_WAVE)
	    ||
	    (Tempc0.orient == NEGATIVE_ORIENTATION &&
		wave_type(Tempc0.curve) == BACKWARD_SHOCK_WAVE))
	    i0_to_i4_dir = CLOCKWISE;
	else
	    i0_to_i4_dir = COUNTER_CLOCK;

	p0 = Node_of_o_curve(&Tempc0)->posn;
	b0 = Bond_at_node_of_o_curve(&Tempc0);
	p4 = Node_of_o_curve(&Tempc4)->posn;
	b4 = Bond_at_node_of_o_curve(&Tempc4);
	find_tangent_to_curve(p0,b0,Tempc0.curve,Tempc0.orient,t0,fr);
	find_tangent_to_curve(p4,b4,Tempc4.curve,Tempc4.orient,t4,fr);
	cos04 = scalar_product(t0,t4,dim);
	(void) vector_product(t0,t4,&sin04,dim);
	ang04 = angle(cos04,sin04);

	if (curve_ang_oriented_l_to_r(i0_to_i4_dir,Tempc0.orient))
	{
	    st01 = Right_state_at_node(Tempc0.curve,Tempc0.orient);
	    ft_assign(RP->state[1],Left_state_at_node(Tempc0.curve,
	    	   Tempc0.orient),fr->sizest);
	}
	else 
	{
	    st01 = Left_state_at_node(Tempc0.curve,Tempc0.orient);
	    ft_assign(RP->state[1],Right_state_at_node(Tempc0.curve,
		   Tempc0.orient),fr->sizest);
	}
	if (curve_ang_oriented_l_to_r(i0_to_i4_dir,Tempc4.orient))
	{
	    st04 = Left_state_at_node(Tempc4.curve,Tempc4.orient);
	    ft_assign(RP->state[4],Right_state_at_node(Tempc4.curve,
	    	   Tempc4.orient),fr->sizest);
	}
	else 
	{
	    st04 = Right_state_at_node(Tempc4.curve,Tempc4.orient);
	    ft_assign(RP->state[4],Left_state_at_node(Tempc4.curve,
		   Tempc4.orient),fr->sizest);
	}
	interpolate_states(fr,0.5,0.5,Coords(tempn->posn),st01,
		           Coords(tempn->posn),st04,RP->state[0]);

	if ((i0_to_i4_dir == COUNTER_CLOCK && ang04 < PI - ERR_ANG) ||
	    (i0_to_i4_dir == CLOCKWISE && ang04 > PI + ERR_ANG))
	{
	    for (i = 0; i < dim; ++i)
	        u[i] = vel(i,RP->state[0]) - Node_vel(tempn)[i];
	    (void) vector_product(u,t0,&cr,dim);
	    is_plus_orientation = (cr > 0.) ? NO : YES;

	    if (find_cross_node_states(Node_vel(tempn),RP,
			               &is_refl_rarefaction,
				       is_plus_orientation)) 
	    {
	    	debug_print("snd_cross_node","Leaving snd_cross_node_propagate()\n");
		if (debugging("snd_cross_node")) 
		{
		    print_cross_node(newn,&Newc0,NULL,NULL,NULL,&Newc4);
		}
		return YES;
	    }
	}

	/* Set mid state */

	for (i = 0; i < dim; ++i)
	    Coords(sten->p[0])[i] = Coords(p0)[i];
	sten->hse[0] = NULL;
	sten->hs[0] = NULL;
	sten->t[0] = ERROR_FLOAT;
	if (curve_ang_oriented_l_to_r(i0_to_i4_dir,Tempc0.orient))
	{
	    ft_assign(sten->rightst[0],RP->state[0],fr->sizest);
	    interpolate_states(fr,0.5,0.5,Coords(tempn->posn),
	    	               RP->state[1],Coords(tempn->posn),RP->state[4],
	    	               sten->leftst[0]);
	}
	else 
	{
	    ft_assign(sten->leftst[0],RP->state[0],fr->sizest);
	    interpolate_states(fr,0.5,0.5,Coords(tempn->posn),RP->state[1],
			       Coords(tempn->posn),RP->state[4],
			       sten->rightst[0]);
	}

	tangent_at_degenerate_node(Tempc0.curve,Tempc0.orient,
				   Tempc4.curve,Tempc4.orient,tgnt,fr);
	ds = grid_size_in_direction(tgnt,newfr->rect_grid->h,dim);

	isgn = (Tempc0.orient == POSITIVE_ORIENTATION) ? 1 : -1;

	/* Set next state */

	states_at_distance_along_curve(p0,b0,Tempc0.curve,Tempc0.orient,
		                       ds,nrad,sten->leftst+isgn,
				       sten->rightst+isgn,
		                       sten->hs+isgn,sten->hse+isgn,
				       sten->t+isgn,sten->p+isgn,newfr);
		
	/* Set prev state */
	if (Tempc0.orient == Tempc4.orient)
	{
	    states_at_distance_along_curve(p4,b4,Tempc4.curve,
			                   Tempc4.orient,ds,nrad,
					   sten->rightst-isgn,
			                   sten->leftst-isgn,sten->hs-isgn,
					   sten->hse-isgn,
			                   sten->t-isgn,sten->p-isgn,newfr);
	}
	else
	{
	    states_at_distance_along_curve(p4,b4,Tempc4.curve,
			                   Tempc4.orient,ds,nrad,
					   sten->leftst-isgn,
				           sten->rightst-isgn,sten->hs-isgn,
					   sten->hse-isgn,
					   sten->t-isgn,sten->p-isgn,newfr);
	}

	ansl0 = Left_state_at_node(Newc0.curve,Newc0.orient);
	ansr0 = Right_state_at_node(Newc0.curve,Newc0.orient);
	ansl4 = Left_state_at_node(Newc4.curve,Newc4.orient);
	ansr4 = Right_state_at_node(Newc4.curve,Newc4.orient);

	sten->newhs = Hyper_surf(Newc0.curve);
	sten->dir = tgnt;
	npt_tang_solver(ds,dt,sten,ansl0,ansr0,fr);

	if (Newc0.orient == Newc4.orient)
	{
	    ft_assign(ansl4,ansr0,fr->sizest);
	    ft_assign(ansr4,ansl0,fr->sizest);
	}
	else
	{
	    ft_assign(ansl4,ansl0,fr->sizest);
	    ft_assign(ansr4,ansr0,fr->sizest);
	}
	if (debugging("snd_cross_node")) 
	{
	    (void) printf("\t\tNEW INCIDENT CURVE0:\n");
	    if (debugging("states"))
	    	show_curve_states(Newc0.curve);
	    else
	    	print_curve(Newc0.curve);
	    (void) printf("\t\tNEW INCIDENT CURVE4:\n");
	    if (debugging("states"))
	    	show_curve_states(Newc4.curve);
	    else
	    	print_curve(Newc4.curve);
	}
	debug_print("snd_cross_node","Leaving snd_cross_node_propagate()\n");
	return YES;
}		/*end snd_cross_node_propagate*/

/* ARGSUSED */
LOCAL boolean snd_SS_node_propagate(
	Front		*fr,
	Front		*newfr,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	/* TODO: CODE NEEDED */
	return YES;
}		/*end snd_SS_node_propagate*/



/* ARGSUSED */
LOCAL boolean snd_SC_node_propagate(
	Front		*fr,
	Front		*newfr,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	/* TODO: CODE NEEDED */
	return YES;
}		/*end snd_SC_node_propagate*/



/* ARGSUSED */
LOCAL boolean snd_CC_node_propagate(
	Front		*fr,
	Front		*newfr,
	NODE		*oldn,
	NODE		*tempn,
	NODE		*newn,
	double		dt)
{
	return YES;
}		/*end snd_CC_node_propagate*/
#endif /* defined(FULL_PHYSICS) */

#endif /* defined(TWOD) */
