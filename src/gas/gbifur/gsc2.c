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
*				gsc2.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Solves two dimensional shock contact Riemann problems.
*
*/

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	void	no_intersection_in_install(POINT*,double,double,double,double,
					   double,ANGLE_DIRECTION,CURVE*);


EXPORT	int install_dfrctn_to_bdry(
	CURVE		**dncur,
	ORIENTATION	*c_orient,
	ANGLE_DIRECTION	i_to_f_dir,
	int		bdry_type,
	COMPONENT	*newcomp,
	RPROBLEM	*rp,
	RP_DATA		*RP,
	Front		*front)
{
	NODE		*diffn = Node_of(dncur[0],c_orient[0]);
	NODE		*oppn = Node_of(dncur[4],Opposite_orient(c_orient[4]));
	NODE		*newn;
	CURVE		*t_bdry, *r_bdry, *ca, *cb;
	CURVE		**curves;
	BOND		*bintsc;
	POINT		*pnew;
	POINT		Newp;
	Locstate	l_st_b[7], r_st_b[7];
	Locstate	l_st_a[7], r_st_a[7];
	COMPONENT	l_c_b[7], r_c_b[7];
	COMPONENT	l_c_a[7], r_c_a[7];
	double		cp, tob[MAXD], tnb[MAXD];
	ORIENTATION	t_bdry_orient, r_bdry_orient, ca_orient, cb_orient;
	int		newnode;
	int		dim = front->rect_grid->dim;
	boolean		sav_scss = interpolate_states_at_split_curve_node();
	int		i;

	debug_print("shock_contact","Entered install_dfrctn_to_bdry()\n");

	/* Identify ahead boundary curves and set state and component lists */

	identify_bdry_curves_at_shock_diffraction(&t_bdry,&t_bdry_orient,
		                                  &r_bdry,&r_bdry_orient,
						  i_to_f_dir,rp);

#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("Initial boundary curves\n");
	    print_orientation("t_bdry, t_bdry_orient = ",t_bdry_orient,"\n");
	    print_curve(t_bdry);
	    show_curve_states(t_bdry);
	    print_orientation("r_bdry, r_bdry_orient = ",r_bdry_orient,"\n");
	    print_curve(r_bdry);
	    show_curve_states(r_bdry);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	for (i = 0; i < 7; i++)
	{
	    l_st_b[i] = r_st_b[i] = NULL;
	    l_st_a[i] = r_st_a[i] = NULL;
	    l_c_b[i] = r_c_b[i] = NO_COMP;
	    l_c_a[i] = r_c_a[i] = NO_COMP;
	}

	for (i = 1; i < 6; i++)
	{
	    if (curve_ang_oriented_l_to_r(i_to_f_dir,t_bdry_orient))
	    {
	    	r_c_b[i] = r_c_a[i] = positive_component(t_bdry);
	    	l_c_b[i] = newcomp[i];
		l_c_a[i] = newcomp[i+1];
	    	l_st_b[i] = RP->state[i];
	    	l_st_a[i] = RP->state[i+1];
	    	if (bdry_type != NEUMANN_BOUNDARY)
	    	{
	    	    r_st_b[i] = l_st_b[i];
	    	    r_st_a[i] = l_st_a[i];
	    	}
	    	else
		{
		    r_st_b[i] = r_st_a[i] = return_obst_state();
		}
	    }
	    else
	    {
	    	l_c_b[i] = l_c_a[i] = negative_component(t_bdry);
	    	r_c_b[i] = newcomp[i];
		r_c_a[i] = newcomp[i+1];
	    	r_st_b[i] = RP->state[i];
	    	r_st_a[i] = RP->state[i+1];
	    	if (bdry_type != NEUMANN_BOUNDARY)
	    	{
	    	    l_st_b[i] = r_st_b[i];
	    	    l_st_a[i] = r_st_a[i];
	    	}
	    	else
	    	{
	    	    l_st_b[i] = l_st_a[i] = return_obst_state();
	    	}
	    }
	}



	tnb[0] = cos(RP->ang[1]);		tnb[1] = sin(RP->ang[1]);
	for (i = 0; i < dim; i++)
	    tob[i] = Coords(oppn->posn)[i] - Coords(diffn->posn)[i];
	(void) vector_product(tob,tnb,&cp,dim);
	if ((i_to_f_dir == CLOCKWISE && cp > 0.) ||
	    (i_to_f_dir == COUNTER_CLOCK && cp < 0.)) 
	{
	    if (!intersect_ray_with_curve(diffn->posn,tnb,NULL,NULL,t_bdry,
					     t_bdry_orient,&bintsc,&Newp)) 
	    {
	    	no_intersection_in_install(diffn->posn,tnb[0],tnb[1],
				           tob[0],tob[1],cp,i_to_f_dir,t_bdry);
		debug_print("shock_contact",
		      "Left install_dfrctn_to_bdry(), ans = NO\n");
		return NO;
	    }
	    shift_node(&Newp,bintsc,t_bdry,t_bdry_orient,r_bdry,
		       r_bdry_orient,oppn,front,l_st_a[5],r_st_a[5],
		       l_st_b[1],r_st_b[1]);
	}
	else 
	{
	    if (!intersect_ray_with_curve(diffn->posn,tnb,NULL,NULL,r_bdry,
					     r_bdry_orient,&bintsc,&Newp)) 
	    {
		no_intersection_in_install(diffn->posn,tnb[0],tnb[1],
					   tob[0],tob[1],cp,i_to_f_dir,r_bdry);
		debug_print("shock_contact",
		      "Left install_dfrctn_to_bdry(), ans = NO\n");
		return NO;
	    }
	    shift_node(&Newp,bintsc,r_bdry,r_bdry_orient,t_bdry,t_bdry_orient,
		       oppn,front,l_st_b[1],r_st_b[1],l_st_a[5],r_st_a[5]);
	}


#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    (void) printf("After first position correction\n");
	    print_orientation("t_bdry, t_bdry_orient = ",t_bdry_orient,"\n");
	    print_curve(t_bdry);
	    show_curve_states(t_bdry);
	    print_orientation("r_bdry, r_bdry_orient = ",r_bdry_orient,"\n");
	    print_curve(r_bdry);
	    show_curve_states(r_bdry);
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	newnode = NO;
	for (i = 1; i < 6; i++)
	{
	    if (dncur[i] == NULL)
		continue;

#if defined(DEBUG_SHOCK_CONTACT)
	    if (debuging("shock_contact")
	        (void) printf"Correcting position of boundary "
			     "node of curve %d\n",i);
#endif /* defined(DEBUG_SHOCK_CONTACT) */

		/* Find ahead boundary curve, comp and states */

	    if (i_to_f_dir == CLOCKWISE)
	    {
		ca = rp->bdry_curves->first->curve;
		ca_orient = rp->bdry_curves->first->orient;
	    }
	    else
	    {
	    	ca = rp->bdry_curves->last->curve;
	    	ca_orient = rp->bdry_curves->last->orient;
	    }

	    if (newnode)
	    {
#if defined(DEBUG_SHOCK_CONTACT)
		if (debugging("shock_contact"))
		{
		    (void) printf("Ahead curve before split\n");
		    print_orientation("ca_orient = ",ca_orient,"\n");
		    print_curve(ca);
		}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

		pnew = Point(Coords(Node_of(ca,ca_orient)->posn));
		interpolate_intfc_states(rp->new_intfc) = NO;
		set_interpolate_states_at_split_curve_node(NO);
		insert_point_adjacent_to_node(pnew,ca,ca_orient);
		if (ca_orient == POSITIVE_ORIENTATION)
		{
	            curves = split_curve(pnew,Bond_at_node(ca,ca_orient),ca,
					 l_c_b[i],r_c_b[i],l_c_a[i],r_c_a[i]);
		    cb = curves[0];
		    cb_orient = NEGATIVE_ORIENTATION;
		}
		else
		{
		    curves = split_curve(pnew,Bond_at_node(ca,ca_orient),ca,
					 l_c_a[i],r_c_a[i],l_c_b[i],r_c_b[i]);
		    cb = curves[1];
		    cb_orient = POSITIVE_ORIENTATION;
		}
		set_interpolate_states_at_split_curve_node(sav_scss);
		interpolate_intfc_states(rp->new_intfc) = YES;
		roclists_after_split(rp,ca,curves,YES);
		ca = (ca_orient==POSITIVE_ORIENTATION) ? curves[1] : curves[0];
#if defined(DEBUG_SHOCK_CONTACT)
		if (debugging("shock_contact"))
		{
		     (void) printf("Ahead curve after split\n");
		     print_orientation("ca_orient = ",ca_orient,"\n");
		     print_curve(ca);
		}
#endif /* defined(DEBUG_SHOCK_CONTACT) */

	        newn = curves[0]->end;
		copy_state(Left_state_at_node(cb,Opposite_orient(cb_orient)),
			   l_st_a[i-1]);
		copy_state(Right_state_at_node(cb,Opposite_orient(cb_orient)),
			   r_st_a[i-1]);
	    }
	    else
	    {
	        if (i_to_f_dir == CLOCKWISE)
	        {
	            cb = rp->bdry_curves->last->curve;
	            cb_orient = rp->bdry_curves->last->orient;
	        }
	        else
	        {
	            cb = rp->bdry_curves->first->curve;
	            cb_orient = rp->bdry_curves->first->orient;
	        }
	        newn = Node_of(ca,ca_orient);
	        newnode = YES;
	    }

	    set_status_at_node(ca,ca_orient,FIXED);
	    set_status_at_node(cb,cb_orient,FIXED);
	    node_type(newn) = bdry_node_type(bdry_type);

	    tnb[0] = cos(RP->ang[i]);
	    tnb[1] = sin(RP->ang[i]);
	    if (!intersect_ray_with_curve(diffn->posn,tnb,NULL,NULL,
					     ca,ca_orient,&bintsc,&Newp)) 
	    {
	        no_intersection_in_install(diffn->posn,tnb[0],tnb[1],
				           tob[0],tob[1],cp,i_to_f_dir,ca);
		debug_print("shock_contact",
		      "Left install_dfrctn_to_bdry(), ans = NO\n");
		return NO;
	    }
	    shift_node(&Newp,bintsc,ca,ca_orient,cb,cb_orient,newn,front,
		       l_st_a[i],r_st_a[i],l_st_b[i],r_st_b[i]);
		
	    /* Install physical curve at new node */

	    change_node_of_curve(dncur[i],Opposite_orient(c_orient[i]),
	    	                 newn);
	}

#if defined(DEBUG_SHOCK_CONTACT)
	if (debugging("shock_contact"))
	{
	    int		j;
	    double	x, y, ang;

	    (void) printf("Angles and states of installed curves\n");
	    for (j = 1; j < 6; j++)
	    {
		if (!dncur[j])
		    continue;
		x = Coords(Node_of(dncur[j],
				   Opposite_orient(c_orient[j]))->posn)[0] -
					   Coords(diffn->posn)[0];
		y = Coords(Node_of(dncur[j],
				   Opposite_orient(c_orient[j]))->posn)[1] - 
					   Coords(diffn->posn)[1];
		(void) printf("\n");
		switch (j)
		{
		case 1:
	    	    (void) printf("Reflected leading edge rarefaction\n");
	    	    print_angle("Given angle =",RP->ang[1],"\n");
	    	    break;
	        case 2:
	    	    (void) printf("Reflected shock\n");
	    	    print_angle("Given angle =",RP->ang[2],"\n");
	    	    break;
	        case 3:
	    	    (void) printf("Reflected trailing edge rarefaction\n");
	    	    print_angle("Given angle =",RP->ang[3],"\n");
		    break;
		case 4:
		    (void) printf("Deflected contact\n");
		    print_angle("Given angle =",RP->ang[4],"\n");
		    break;
		case 5:
		    (void) printf("Transmitted shock\n");
		    print_angle("Given angle =",RP->ang[5],"\n");
		    break;
		}
		ang = angle(x,y);
		print_angle("Installed angle =",ang,"\n");
		(void) printf("States on new curve\n");
		verbose_print_curve_states(dncur[j]);
	    }
	}
#endif /* defined(DEBUG_SHOCK_CONTACT) */
	debug_print("shock_contact",
	      "Left install_dfrctn_to_bdry(), ans = YES\n");
	return YES;
}		/*end install_dfrctn_to_bdry*/


EXPORT	void identify_bdry_curves_at_shock_diffraction(
	CURVE		**t_bdry,
	ORIENTATION	*t_bdry_orient,
	CURVE		**r_bdry,
	ORIENTATION	*r_bdry_orient,
	ANGLE_DIRECTION	i_to_f_dir,
	RPROBLEM	*rp)
{
	static int	num_invert = 0;

	if (i_to_f_dir == CLOCKWISE)
	{
	    *t_bdry = rp->bdry_curves->first->curve;
	    *t_bdry_orient = rp->bdry_curves->first->orient;
	    *r_bdry = rp->bdry_curves->last->curve;
	    *r_bdry_orient = rp->bdry_curves->last->orient;
	}
	else
	{
	    *r_bdry = rp->bdry_curves->first->curve;
	    *r_bdry_orient = rp->bdry_curves->first->orient;
	    *t_bdry = rp->bdry_curves->last->curve;
	    *t_bdry_orient = rp->bdry_curves->last->orient;
	}
	if (*t_bdry_orient == Opposite_orient(*r_bdry_orient))
	    return;

	(void) printf("WARNING in identify_bdry_curves_at_shock_diffraction(), "
	              "Unexpected case\n"
	              "Inconsistent boundary orientations\n");

	if (is_exterior_comp(positive_component((*t_bdry)),rp->new_intfc))
	{
	    invert_curve(*t_bdry);
	    roclists_after_invert(rp,*t_bdry,(O_CURVE *)NULL);
	    *t_bdry_orient = Opposite_orient(*t_bdry_orient);
	}
	else if (is_exterior_comp(positive_component((*r_bdry)),rp->new_intfc))
	{
	    invert_curve(*r_bdry);
	    roclists_after_invert(rp,*r_bdry,(O_CURVE *)NULL);
	    *r_bdry_orient = Opposite_orient(*r_bdry_orient);
	}
	else if ((num_invert++)%2)
	{
	    invert_curve(*t_bdry);
	    roclists_after_invert(rp,*t_bdry,(O_CURVE *)NULL);
	    *t_bdry_orient = Opposite_orient(*t_bdry_orient);
	}
	else
	{
	    invert_curve(*r_bdry);
	    roclists_after_invert(rp,*r_bdry,(O_CURVE *)NULL);
	    *r_bdry_orient = Opposite_orient(*r_bdry_orient);
	}

}		/*end identify_bdry_curves_at_shock_diffraction*/

LOCAL	void no_intersection_in_install(
	POINT		*p0,
	double		tnbx,
	double		tnby,
	double		tobx,
	double		toby,
	double		cp,
	ANGLE_DIRECTION	i_to_f_dir,
	CURVE		*c_bdry)
{
	(void) printf("WARNING in install_dfrctn_to_bdry(), "
	              "NO CROSS of ray with boundary\n");
	(void) printf("p = (%g, %g), tn =  <%g, %g>,\n\t",
		      Coords(p0)[0],Coords(p0)[1],tnbx,tnby);
	(void) printf("to = <%g, %g>\n",tobx,toby);
	(void) printf("cp = %g\n",cp);
	print_angle_direction("i_to_f_dir = ",i_to_f_dir,"\n");
	print_curve(c_bdry);
	if (debugging("shock_contact"))
	    print_interface(c_bdry->interface);
}		/*end no_intersection_in_install*/
#endif /* defined(FULL_PHYSICS) && defined(TWOD) */
