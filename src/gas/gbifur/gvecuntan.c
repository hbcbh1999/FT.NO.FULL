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
*			gvecuntan.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Resolves untangles of interior curves.
*
*/


#if defined(TWOD)
#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	int	find_curves_at_vec_bdry_untangle(CURVE*,CURVE*,CURVE**,
						 ORIENTATION,ANGLE_DIRECTION,
						 int,NODE**,CURVE**,
						 CURVE**,CURVE**,CURVE**,
						 ORIENTATION*,ORIENTATION*,
						 ANGLE_DIRECTION*);
LOCAL	int	set_states_at_new_B_refl_node(NODE*,CURVE*,CURVE*,CURVE*,
					      CURVE*,ANGLE_DIRECTION,
					      ORIENTATION,ORIENTATION,Front*);
LOCAL	void	set_new_bdry_states(int,NODE*,ANGLE_DIRECTION,
				    ORIENTATION,CURVE*,Front*);
LOCAL	void	set_orients_and_w_type_for_cross(ORIENTATION*,int*);
LOCAL	void	set_orients_and_w_type_for_overtake(ORIENTATION*,int*);
LOCAL	void	set_states_and_comps_about_cross(COMPONENT*,ANGLE_DIRECTION,
						 ORIENTATION*,
						 COMPONENT*,COMPONENT*,
						 Locstate*,Locstate*,RP_DATA*);
LOCAL	void	set_track_and_comp_list_for_cross(boolean*,COMPONENT*,
						  ANGLE_DIRECTION,ORIENTATION*,
						  Front*,CURVE*,CURVE*,
						  int*,RP_DATA*);
LOCAL	void	set_track_and_comp_list_for_overtake(boolean*,COMPONENT*,
						     int,ANGLE_DIRECTION,
						     ORIENTATION*,
						     Front*,CURVE*,CURVE*,
						     RP_DATA*);


/*
*			vector_overtake_unravel():
*
*	Attempts to unravel vector overtake nodes.  
*	Called by vector_vector_unravel in guntan.c
*
*/

EXPORT int vector_overtake_unravel(
	Front		*front,
	CROSS		*cross,
	CROSS		*cross1,
	ORIENTATION	c1_orient,
	ORIENTATION	c2_orient)
{
	BOND		*tbond;
	COMPONENT	comp[7];
	CROSS		*cr[2];
	CURVE		*curves1[2], *curves2[2], *tcurve, **c;
	CURVE		*ctmp = NULL;
	Locstate	ls, rs;
	NODE		*otake_node[2], *ns, *ne;
	RP_DATA		*RP[2];
	double		nod_v[SMAXD], coords[MAXD], cp;
	double		t0[MAXD];
	int		w_type[2][7];
	boolean		sav_intrp;
	int		c1_incident;
	ORIENTATION	tmp_orient;
	boolean		is_plus_orientation[2];
	boolean		is_refl_rarefaction[2];
	ANGLE_DIRECTION	i_to_o_dir[2];
	boolean		track[7];
	int		reset_comp;
	int		i, j, dim = front->interf->dim;

	static CURVE	   ***oncur = NULL;
	static COMPONENT   **l_comp = NULL, **r_comp = NULL;
	static Locstate    **l_st = NULL, **r_st = NULL;
	static double	   ***t = NULL;
	static double	   *t1 = NULL;
	static ORIENTATION **c_orient = NULL;
	static int	status[7] = { OVERTOOK, INCIDENT, REFLECTED, REFLECTED,
					REFLECTED, SLIP, TRANSMITTED };



	debug_print("untangle","Entered vector_overtake_unravel\n");
	if (t == NULL)
	{
		uni_array(&t1,dim,FLOAT);
		bi_array(&oncur, 2, 7, sizeof(CURVE *));
		bi_array(&c_orient, 2, 7, sizeof(ORIENTATION));
		tri_array(&t,2,2,MAXD,FLOAT);
		bi_array(&l_comp,2,7,sizeof(COMPONENT));
		bi_array(&r_comp,2,7,sizeof(COMPONENT));
		bi_array(&l_st,2,7,sizeof(Locstate));
		bi_array(&r_st,2,7,sizeof(Locstate));
	}
	sav_intrp = interpolate_intfc_states(front->interf);
	interpolate_intfc_states(front->interf) = YES;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 7; j++)
			oncur[i][j] = NULL;

		/* Rearrange curves at the nodes so that 
		   curve1 overtakes curve2. */

	c1_incident = YES;
	vector_product_on_bonds(cross->b1,cross->b2,2,&cp);
	if (cp < 0.) c1_incident = !c1_incident;
	if (c2_orient == POSITIVE_ORIENTATION)	    c1_incident = !c1_incident;
	if (is_backward_wave(wave_type(cross->c1))) c1_incident = !c1_incident;

	if (!c1_incident) 
	{
		debug_print("untangle","Re-orienting the crossed curves\n");
		tcurve = cross->c1;
		cross->c1 = cross->c2;
		cross1->c1 = cross1->c2;
		cross->c2 = tcurve;
		cross1->c2 = tcurve;

		tbond = cross->b1;
		cross->b1 = cross->b2;
		cross->b2 = tbond;
		tbond = cross1->b1;
		cross1->b1 = cross1->b2;
		cross1->b2 = tbond;

		tmp_orient = c1_orient;
		c1_orient = c2_orient;
		c2_orient = tmp_orient;
	}
	w_type[0][1] = w_type[1][1] = wave_type(cross->c1);
	w_type[0][0] = w_type[1][0] = wave_type(cross->c2);

		/* Mark the overtaken wave so that it can be 
		   identified after the curves are split */

	wave_type(cross->c2) = MARKED_CURVE;

	c_orient[0][1] = c1_orient;
	c_orient[0][0] = c2_orient;
	c_orient[1][1] = Opposite_orient(c1_orient);
	c_orient[1][0] = Opposite_orient(c2_orient);

	cr[0] = cross;
	cr[1] = cross1;
	for (i = 0; i < 2; i++)
	{
		split_curves_at_cross(cr[i], front, &otake_node[i], curves1,
			              NULL,NULL,curves2,NULL,NULL,0.0,NULL);
	}

	delete_from_cross_list(cross1);

	for (i = 0; i < 2; i++)
	{
		if (debugging("untangle"))
			(void) printf("start loop %d\n", i);
		node_type(otake_node[i]) = OVERTAKE_NODE;
		RP[i] = Rp_data(otake_node[i]);

		set_orients_and_w_type_for_overtake(c_orient[i], w_type[i]);
		debug_print("untangle","After set orients\n");

		if ( ((c_orient[i][1]== NEGATIVE_ORIENTATION) &&
			is_forward_wave(w_type[i][1]))
			 ||
		     ((c_orient[i][1]== POSITIVE_ORIENTATION) &&
			is_backward_wave(w_type[i][1])) )
		{
			i_to_o_dir[i] = COUNTER_CLOCK;
			is_plus_orientation[i] = NO;
		}
		else
		{
			i_to_o_dir[i] = CLOCKWISE;
			is_plus_orientation[i] = YES;
		}

		RP[i]->ang_dir = Opposite_ang_dir(i_to_o_dir[i]); 

		/* Find the incident and overtaken curves 
		   oncur[i][j]    j = 	0 overtaken curve
					1 incident curve
					5 slip line
					6 transmitted shock */

		if(debugging("untangle"))
		{
			(void) printf("Overtake node %d \n", i);
			print_node(otake_node[i]);
		}

		for (c = otake_node[i]->in_curves; c && *c; c++)
		{
			if (wave_type(*c) == MARKED_CURVE)
			{
				if (c_orient[i][0] == NEGATIVE_ORIENTATION)
					oncur[i][0] = *c;
				else oncur[i][5] = *c;
			}
			else
			{
				if (c_orient[i][1] == NEGATIVE_ORIENTATION)
					oncur[i][1] = *c;
				else oncur[i][6] = *c;
			}
		}
				
		for (c = otake_node[i]->out_curves; c && *c; c++)
		{
			if (wave_type(*c) == MARKED_CURVE)
			{
				if (c_orient[i][0] == POSITIVE_ORIENTATION)
					oncur[i][0] = *c;
				else oncur[i][5] = *c;
			}
			else
			{
				if (c_orient[i][1] == POSITIVE_ORIENTATION)
					oncur[i][1] = *c;
				else oncur[i][6] = *c;
			}
		}

		debug_print("untangle","Curves identified \n");
	}

		/* Need to break the loop up since the wave types
		   for the common curves cannot be changed until 
		   both are identified */

	for(i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
		{
			wave_type(oncur[i][j]) = w_type[i][j];
			wave_type(oncur[i][5+j]) = w_type[i][j];
		}

		/* Find the states needed to solve RP */

		if (curve_ang_oriented_l_to_r(i_to_o_dir[i], c_orient[i][0]))
		{
			ft_assign(RP[i]->state[0],
				Right_state_at_node(oncur[i][0],c_orient[i][0]),
				front->sizest);
		}
		else
		{
			ft_assign(RP[i]->state[0],
				Left_state_at_node(oncur[i][0],c_orient[i][0]),
				front->sizest);
		}

		if (curve_ang_oriented_l_to_r(i_to_o_dir[i], c_orient[i][1]))
		{
			ft_assign(RP[i]->state[1],
				Right_state_at_node(oncur[i][1],c_orient[i][1]),
				front->sizest);
			ft_assign(RP[i]->state[2],
				Left_state_at_node(oncur[i][1],c_orient[i][1]),
				front->sizest);
		}
		else
		{
			ft_assign(RP[i]->state[1],
				Left_state_at_node(oncur[i][1],c_orient[i][1]),
				front->sizest);
			ft_assign(RP[i]->state[2],
				Right_state_at_node(oncur[i][1],c_orient[i][1]),
				front->sizest);
		}

		debug_print("untangle", "Input states for Riemann problem found\n");


		if (!estimate_node_vel(oncur[i][0],c_orient[i][0],
			t[i][0],oncur[i][1],c_orient[i][1], t[i][1],
			i_to_o_dir[i],node_type(otake_node[i]),nod_v,front))
		{
		    (void) printf("WARNING in vector_overtake_unravel(), "
		                  "Cannot estimate node velocity.\n");
		    goto fail;
		}
		if (debugging("untangle"))
		{
		    (void) printf("Estimated velocity (node %d): %g %g \n",
				  i,nod_v[0],nod_v[1]);
		}

		if (!find_overtake_node_states(nod_v, RP[i], 
			&is_refl_rarefaction[i], is_plus_orientation[i]))
		{
		    (void) printf("WARNING in vector_overtake_unravel(), "
		                  "Cannot find overtake node states\n");
		    goto fail;
		} 
		debug_print("untangle", "After overtake node states found\n");
	}

		/* Assign Components */
	
	set_track_and_comp_list_for_overtake(track, comp,
		is_refl_rarefaction[0], i_to_o_dir[0], c_orient[0],
		front, oncur[0][1], oncur[0][0], RP[0]);

	for (i = 0; i < 2; i++)
	{

		/* This routine should be renamed since it has a 
		   more general use than just the shock_diffraction case */
		   
		set_states_and_comps_about_shock_diffraction(
			comp, i_to_o_dir[i], c_orient[i], l_comp[i],
			r_comp[i], l_st[i], r_st[i], RP[i]);
	}
	debug_print("untangle","Components and states set\n");

	if (!track[6])
	{
		(void) delete_curve(oncur[0][6]);
		oncur[0][6] = oncur[1][6] = NULL;
	}
	else delete_interior_points_of_curve(front,oncur[0][6]);
	if (!track[5])
	{
		(void) delete_curve(oncur[0][5]);
		oncur[0][5] = oncur[1][5] = NULL;
	}
	else delete_interior_points_of_curve(front,oncur[0][5]);

	reset_comp = YES;
	for (j = 2; j < 7; j++) if (track[j])
	{
		reset_comp = NO;
		break;
	}

	if (reset_comp)
	{
		if (i_to_o_dir[0] == COUNTER_CLOCK)
		{
			ns = otake_node[0];
			ne = otake_node[1];
		}
		else
		{
			ns = otake_node[1];
			ne = otake_node[0];
		}
		ctmp = make_curve(comp[2],comp[0],ns,ne);
	}


	for (j = 0; j < 7; j++)
	{
		if (!track[j]) continue;
		if (oncur[0][j] == NULL)
		{
			if (c_orient[0][j] == POSITIVE_ORIENTATION)
			{
				ns = otake_node[0];
				ne = otake_node[1];
			}
			else
			{
				ns = otake_node[1];
				ne = otake_node[0];
			}
			oncur[0][j] = oncur[1][j] =
				make_curve(l_comp[0][j], r_comp[0][j],ns,ne);
		}
		else
		{
			negative_component(oncur[0][j]) = l_comp[0][j];
			positive_component(oncur[0][j]) = r_comp[0][j];
		}
		wave_type(oncur[0][j]) = w_type[0][j];

		/* Assign states at nodes */

		for (i = 0; i < 2; i++)
		{
			set_status_at_node(oncur[i][j], c_orient[i][j],
				status[j]);
			ls = Left_state_at_node(oncur[i][j], c_orient[i][j]);
			rs = Right_state_at_node(oncur[i][j], c_orient[i][j]);
			ft_assign(ls, l_st[i][j], front->sizest);
			ft_assign(rs, r_st[i][j], front->sizest);
		}
	}

	/*
	*  If none of the connecting curves from otake_node[0] and 
	*  otake_node[1] are tracked,  then reset the component of
	*  the now connected regions.
	*/

	if (reset_comp)
	{
		reset_component_of_loop(oncur[0][1],c_orient[0][1],
			i_to_o_dir[0],comp[0],front);
		(void) delete_curve(ctmp);
	}

	/* Insert points in curves to satisfy angles */

	for (j = 2; j < 7; j++)
	{
		if (!track[j]) continue;
		if (debugging("untangle"))
		{
			char mesg[20];
			int i;

			(void) printf("Inserting new point in curve %d \n",j);
			for (i = 0; i < 2; i++)
			{
				(void) printf("Position [%d][%d] = <%g, %g>, ",i,j,
					Coords(otake_node[i]->posn)[0],
					Coords(otake_node[i]->posn)[1]);
				(void) sprintf(mesg,"Angle [%d][%d] = ",i,j);
				print_angle(mesg,RP[i]->ang[j],"\n");
			}
		} 
		t0[0] = cos(RP[0]->ang[j]);	t0[1] = sin(RP[0]->ang[j]);
		t1[0] = cos(RP[1]->ang[j]);	t1[1] = sin(RP[1]->ang[j]);
		if (!intersect_ray_with_sector(otake_node[0]->posn,
			otake_node[1]->posn, t0, &t1,coords,dim))
		{
		    (void) printf("WARNING in vector_overtake_unravel(), ");
		    (void) printf("\tInconsistent orientations ");
		    (void) printf("at oppst nodes, j = %d\n",j);
		    goto fail;
		}
		if (debugging("untangle"))
		{
		    (void) printf("Intersection found at point <%g, %g>\n",
				  coords[0],coords[1]);
		}
		insert_point_adjacent_to_node(Point(coords), oncur[0][j],
			c_orient[0][j]);
	}
	if(debugging("untangle"))
	{
		(void) printf("New curve orientations are established\n");
		for (i = 0; i < 2; i++)
		{
			(void) printf("New overtake node\n");
			print_node(otake_node[i]);
			for (j = 0; j < 7; j++)
			{
			    if( !track[j]) continue;
			    (void) printf("\n oncur[%d][%d] : \n",i,j);
			    print_orientation("Orient =",c_orient[i][j],"\n");
			    print_curve(oncur[i][j]);
			}
		}
	}

		/* Success !!! */

	interpolate_intfc_states(front->interf) = sav_intrp;
	debug_print("untangle", 
		"Successful completion of vector_overtake_unravel\n");
	return CURVES_UNTANGLED; 


	/* Reset old interpolate setting report failure */

fail:   interpolate_intfc_states(front->interf) = sav_intrp;
	debug_print("untangle", 
		"Unsuccessful completion of vector_overtake_unravel\n");
	return ERROR_IN_UNTANGLE;
}		/*end vector_overtake_unravel*/

LOCAL	void set_orients_and_w_type_for_overtake(
	ORIENTATION		*c_orient,
	int			*w_type)
{
	c_orient[2] = Opposite_orient(c_orient[1]);
	c_orient[3] = Opposite_orient(c_orient[1]);
	c_orient[4] = Opposite_orient(c_orient[1]);
	c_orient[5] = Opposite_orient(c_orient[0]);
	c_orient[6] = Opposite_orient(c_orient[1]);

	w_type[5] = CONTACT;
	if (is_forward_wave(w_type[1]))
	{
		w_type[2] = BACKWARD_SOUND_WAVE_LE;
		w_type[3] = BACKWARD_SHOCK_WAVE;
		w_type[4] = BACKWARD_SOUND_WAVE_TE;
		w_type[6] = FORWARD_SHOCK_WAVE;
	}
	else
	{
		w_type[2] = FORWARD_SOUND_WAVE_LE;
		w_type[3] = FORWARD_SHOCK_WAVE;
		w_type[4] = FORWARD_SOUND_WAVE_TE;
		w_type[6] = BACKWARD_SHOCK_WAVE;
	}
}		/*end set_orients_and_w_type_for_overtake*/

LOCAL	void set_track_and_comp_list_for_overtake(
	boolean		*track,
	COMPONENT	*comp,
	int		is_refl_rarefaction,
	ANGLE_DIRECTION	i_to_o_dir,
	ORIENTATION	*c_orient,
	Front		*fr,
	CURVE		*incident,
	CURVE		*overtaken,
	RP_DATA		*RP)
{
	int		i, j;

	/* Set tracking of incident waves */
	track[0] = track[1] = YES;
	if (is_rarefaction_wave(wave_type(incident)) || 
				is_rarefaction_wave(wave_type(overtaken)))
	{
	    track[2] = track[3] = track[4] = track[5] = NO;
	}
	else
	{
	    /* Set tracking of reflected waves */
	    if (is_refl_rarefaction)
	    {
	        track[3] = NO;
	        track[2] = track[4] =
		    track_scattered_wave(OVERTAKE_NODE,RAREF_LEADING_EDGE,
					  REFLECTED,RP->state[2],RP->state[5],
					  fr);
	    }
	    else
	    {
	        track[2] = track[4] = NO;
	        track[3] = track_scattered_wave(OVERTAKE_NODE,SHOCK_WAVE,
	    				         REFLECTED,RP->state[2],
	    			                 RP->state[5],fr);
	    }

	    /* Set tracking of slip line */
	    track[5] = track_scattered_wave(OVERTAKE_NODE,CONTACT_WAVE,SLIP,
					     RP->state[5],RP->state[6],fr);
	}

		/* Set tracking of transmitted wave */
	
	track[6] = track_scattered_wave(OVERTAKE_NODE,SHOCK_WAVE,TRANSMITTED,
					 RP->state[0],RP->state[6],fr);

		/* Set components */

	if (curve_ang_oriented_l_to_r(i_to_o_dir, c_orient[0]))
	{
	    comp[0] = positive_component(overtaken);
	    comp[1] = negative_component(overtaken);
	}
	else
	{
	    comp[1] = positive_component(overtaken);
	    comp[0] = negative_component(overtaken);
	}
	comp[2] = (curve_ang_oriented_l_to_r(i_to_o_dir, c_orient[1])) ?
				             negative_component(incident) :
				             positive_component(incident);

	for (j = 3; j < 7; j++)
	{
	    if (track[j-1] == NO)
	    	comp[j] = comp[j-1];
	    else
	    {
	    	comp[j] = comp[0];
	    	for (i = j; i < 7; i++)
	    	{
	    	    if (track[i])
	    	    {
	    		comp[j] = new_component(NEW_COMP);
	    		break;
	    	    }
	    	}
	    }
	}
}		/*end set_track_and_comp_lists_for_overtake*/

/*
*			vector_cross_unravel():
*
*	Attempts to unravel vector cross nodes.
*	Called by vector_vector_unravel() in guntan.c.
*
*/

EXPORT int vector_cross_unravel(
	Front		*front,
	CROSS		*cross,
	CROSS		*cross1,
	ORIENTATION	c1_orient,
	ORIENTATION	c2_orient)
{
	CURVE           *curves1[2], *curves2[2], **c;
	RP_DATA		*RP[2];
	NODE		*cross_node[2], *ns, *ne;
	CROSS		*cr[2];
	COMPONENT	comp[5];
	Locstate	ls, rs;
	int		w_type[2][5], join;
	boolean		sav_intrp;
	int		i, j;
	double		nod_v[SMAXD], coords[MAXD];
	double		t0[MAXD];
	boolean		is_plus_orientation[2];
	boolean		is_refl_rarefaction[2];
	ANGLE_DIRECTION	i0_to_i4_dir[2];
	boolean		track[5];
	int		dim = front->interf->dim;

	static COMPONENT   **l_comp = NULL, **r_comp = NULL;
	static Locstate    **l_st = NULL, **r_st = NULL;
	static	CURVE	    ***cncur = NULL;
	static	double	    ***t = NULL;
	static	double	    *t1 = NULL;
	static	ORIENTATION **c_orient = NULL;
	static	int	status[5] = { INCIDENT,
				      REFLECTED,
				      SLIP,
				      REFLECTED,
				      INCIDENT}; 


	debug_print("untangle","Entered vector_cross_unravel\n");
	if (t1 == NULL)
	{
	    uni_array(&t1,dim,FLOAT);
	    bi_array(&cncur, 2, 5, sizeof(CURVE *));
	    bi_array(&c_orient, 2, 5, sizeof(ORIENTATION));
	    tri_array(&t,2,2,MAXD,FLOAT);
	    bi_array(&l_comp,2,7,sizeof(COMPONENT));
	    bi_array(&r_comp,2,7,sizeof(COMPONENT));
	    bi_array(&l_st,2,7,sizeof(Locstate));
	    bi_array(&r_st,2,7,sizeof(Locstate));
	}
	sav_intrp = interpolate_intfc_states(front->interf);
	interpolate_intfc_states(front->interf) = YES;
	for (i = 0; i < 2; i++)
	    for (j = 0; j < 5; j++) cncur[i][j] = NULL;

	w_type[0][0] = w_type[1][0] = wave_type(cross->c1);
	w_type[0][4] = w_type[1][4] = wave_type(cross->c2);

	    /* Mark incident 0 curve so that it can be identified
	       after the split */

	wave_type(cross->c1) = MARKED_CURVE;

	c_orient[0][0] = c1_orient;
	c_orient[0][4] = c2_orient;
	c_orient[1][0] = Opposite_orient(c1_orient);
	c_orient[1][4] = Opposite_orient(c2_orient);

	cr[0] = cross;
	cr[1] = cross1;
	for (i = 0; i < 2; i++)
	{
	    split_curves_at_cross(cr[i],front,&cross_node[i],curves1,
			          NULL,NULL,curves2,NULL,NULL,0.0,NULL);
	}

	delete_from_cross_list(cross1);

	for (i = 0; i < 2; i++)
	{
	    if (debugging("untangle"))
	    	(void) printf("start loop %d\n", i);
	    node_type(cross_node[i]) = CROSS_NODE;
	    RP[i] = Rp_data(cross_node[i]);

	    set_orients_and_w_type_for_cross(c_orient[i], w_type[i]);

	    debug_print("untangle","After set orients\n");

	    if ( ((c_orient[i][0]== NEGATIVE_ORIENTATION) &&
	    	is_forward_wave(w_type[i][0]))
	    	 ||
	         ((c_orient[i][0]== POSITIVE_ORIENTATION) &&
	    	is_backward_wave(w_type[i][0])) )
	    {
	    	i0_to_i4_dir[i] = COUNTER_CLOCK;
	    	is_plus_orientation[i] = NO;
	    }
	    else
	    {
	    	i0_to_i4_dir[i] = CLOCKWISE;
	    	is_plus_orientation[i] = YES;
	    }

	    RP[i]->ang_dir = Opposite_ang_dir(i0_to_i4_dir[i]); 

	    /* Find the incident and reflected curves 
	       cncur[i][j]    j = 0 incident curve	
	    			  4 incident curve
	    		 	  1 reflected (4)	
	    			  3 reflected (0) */

	    if(debugging("untangle"))
	    {
	    	(void) printf("Overtake node %d \n", i);
	    	print_node(cross_node[i]);
	    }

	    for (c = cross_node[i]->in_curves; c && *c; c++)
	    {
	    	if (wave_type(*c) == MARKED_CURVE)
	    	{
	    	    if (c_orient[i][0] == NEGATIVE_ORIENTATION)
	    	    	cncur[i][0] = *c;
	    	    else cncur[i][3] = *c;
	    	}
	    	else
	    	{
	    	    if (c_orient[i][4] == NEGATIVE_ORIENTATION)
	    	    	cncur[i][4] = *c;
	    	    else cncur[i][1] = *c;
	    	}
	    }
				
	    for (c = cross_node[i]->out_curves; c && *c; c++)
	    {
	    	if (wave_type(*c) == MARKED_CURVE)
	    	{
	    	    if (c_orient[i][0] == POSITIVE_ORIENTATION)
	            	cncur[i][0] = *c;
	            else cncur[i][3] = *c;
	    	}
	    	else
	    	{
	    	    if (c_orient[i][4] == POSITIVE_ORIENTATION)
	    	    	cncur[i][4] = *c;
	            else cncur[i][1] = *c;
		}
	    }

	    debug_print("untangle","Curves identified \n");
	}

	    /* Need to break the loop up since the wave types
	       for the common curves cannot be changed until 
	       both are identified */

	for(i = 0; i < 2; i++)
	{
	    for (j = 0; j < 2; j++)
	    	wave_type(cncur[i][j*4]) = w_type[i][j*4];

	    /* Find the states needed to solve RP */

	    if (curve_ang_oriented_l_to_r(i0_to_i4_dir[i], c_orient[i][0]))
	    {
	    	ft_assign(RP[i]->state[0],
	    	       Right_state_at_node(cncur[i][0], c_orient[i][0]),
		       front->sizest);
		ft_assign(RP[i]->state[1],
		       Left_state_at_node(cncur[i][0], c_orient[i][0]),
		       front->sizest);
	    }
	    else
	    {
	    	ft_assign(RP[i]->state[0],
	    	       Left_state_at_node(cncur[i][0], c_orient[i][0]),
		       front->sizest);
		ft_assign(RP[i]->state[1],
		       Right_state_at_node(cncur[i][0], c_orient[i][0]),
		       front->sizest);
	    }

	    if (curve_ang_oriented_l_to_r(i0_to_i4_dir[i], c_orient[i][4]))
	    {
	    	ft_assign(RP[i]->state[4],
	    	       Left_state_at_node(cncur[i][4], c_orient[i][4]),
	    	       front->sizest);
	    }
	    else
	    {
	    	ft_assign(RP[i]->state[4],
	    	       Right_state_at_node(cncur[i][4], c_orient[i][4]),
	    	       front->sizest);
	    }

	    debug_print("untangle", "Input states for Riemann problem found\n");

	    /* Note: Should do an average for state 1 */

	    if (!estimate_node_vel(cncur[i][4],c_orient[i][4],
	    	                      t[i][1],cncur[i][0],
				      c_orient[i][0],t[i][0],
			              i0_to_i4_dir[i],
				      node_type(cross_node[i]),nod_v,front))
	    {
	    	(void) printf("WARNING in vector_cross_unravel(), "
	    	              "Cannot estimate node velocity.\n");
	    	interpolate_intfc_states(front->interf) = sav_intrp;
	    	return ERROR_IN_UNTANGLE;
	    }
	    if (debugging("untangle"))
	       (void) printf("Estimated velocity: %g %g \n",nod_v[0],nod_v[1]);

	    if (!find_cross_node_states(nod_v,RP[i],&is_refl_rarefaction[i],
					   is_plus_orientation[i]))
	    {
	    	(void) printf("WARNING in vector_cross_unravel(), "
	    	              "Cannot find cross node states\n");
		interpolate_intfc_states(front->interf) = sav_intrp;
		return ERROR_IN_UNTANGLE;
	    } 
	    if (is_refl_rarefaction[i])
	    {
	        (void) printf("WARNING in vector_cross_unravel(), "
		              "find_cross_node_states has rarefaction\n");
	        (void) printf("No code to deal with this case\n");
		interpolate_intfc_states(front->interf) = sav_intrp;
		return ERROR_IN_UNTANGLE;
	    }

	    debug_print("untangle", "After cross node states found\n");
	}

	    /* Assign Components */
	
	set_track_and_comp_list_for_cross(track,comp,i0_to_i4_dir[0],
					  c_orient[0],front,cncur[0][0],
					  cncur[0][4],&join,RP[0]);

	for (i = 0; i < 2; i++)
	{
	    set_states_and_comps_about_cross(comp,i0_to_i4_dir[i],c_orient[i],
					     l_comp[i],r_comp[i],l_st[i],
					     r_st[i],RP[i]);
	}
	debug_print("untangle","Components and states set\n");

	if (!track[1])
	{
	    (void) delete_curve(cncur[0][1]);
	    cncur[0][1] = cncur[1][1] = NULL;
	}
	else
	    delete_interior_points_of_curve(front,cncur[0][1]);

	if (!track[3])
	{
	    (void) delete_curve(cncur[0][3]);
	    cncur[0][3] = cncur[1][3] = NULL;
	}
	else
	    delete_interior_points_of_curve(front,cncur[0][3]);

	for (j = 0; j < 5; j++)
	{
	    if (!track[j])
		continue;
	    if (cncur[0][j] == NULL)
	    {
	    	if (c_orient[0][j] == POSITIVE_ORIENTATION)
	    	{
	    	    ns = cross_node[0];
	    	    ne = cross_node[1];
	    	}
	    	else
	    	{
	    	    ns = cross_node[1];
	    	    ne = cross_node[0];
	    	}
	    	cncur[0][j] = cncur[1][j] = make_curve(l_comp[0][j],
						       r_comp[0][j],ns,ne);
	    }
	    else
	    {
	    	negative_component(cncur[0][j]) = l_comp[0][j];
	    	positive_component(cncur[0][j]) = r_comp[0][j];
	    }
	    wave_type(cncur[0][j]) = w_type[0][j];

	    /* Assign states at nodes */

	    for (i = 0; i < 2; i++)
	    {
	    	set_status_at_node(cncur[i][j], c_orient[i][j],status[j]);
		ls = Left_state_at_node(cncur[i][j], c_orient[i][j]);
		rs = Right_state_at_node(cncur[i][j], c_orient[i][j]);
		ft_assign(ls, l_st[i][j], front->sizest);
		ft_assign(rs, r_st[i][j], front->sizest);
	     }
	}
	debug_print("untangle", "After new curves and states are found\n");

	/* Join components together if there are no tracked reflected or
	   transmitted waves. */

	if (join)
	{
	     reset_component_of_loop(cncur[0][0],
				     Opposite_orient(c_orient[0][0]),
			             Opposite_ang_dir(i0_to_i4_dir[0]), 
			             comp[1],front);
	}

	/* Insert points in curves to satisfy angles */

	for (j = 1; j < 4; j++)
	{
	    if (!track[j])
		continue;
	    t0[0] = cos(RP[0]->ang[j]);	t0[1] = sin(RP[0]->ang[j]);
	    t1[0] = cos(RP[1]->ang[j]);	t1[1] = sin(RP[1]->ang[j]);
	    if (!intersect_ray_with_sector(cross_node[0]->posn,
			                      cross_node[1]->posn,
					      t0,&t1,coords,dim))
	    {
	    	(void) printf("WARNING in vector_cross_unravel(), "
	    	              "Inconsistent orientations at oppst nodes\n");
	    	interpolate_intfc_states(front->interf) = sav_intrp;
	    	return ERROR_IN_UNTANGLE;
	    }
	    insert_point_adjacent_to_node(Point(coords), cncur[0][j],
			                  c_orient[0][j]);
	}
	if(debugging("untangle"))
	{
	    (void) printf("New curve orientations are established\n");
	    for (j = 0; j < 5; j++)
	    {
	    	if( !track[j])
		    continue;
	    	(void) printf("\n cncur[0][ %d ] : \n", j);
	    	print_curve(cncur[0][j]);
	    }
	}

		/* Success !!! */

	interpolate_intfc_states(front->interf) = sav_intrp;
	debug_print("untangle", 
		"Successful completion of vector_cross_unravel\n");
	return CURVES_UNTANGLED; 
}		/*end vector_cross_unravel*/

EXPORT	int g_vec_bdry_untangle(
	CURVE		*cinterior,
	CURVE		*cexterior,
	CURVE		**bdrycurves,
	ORIENTATION	cphys_orient,
	ANGLE_DIRECTION	cb_to_cp_dir,
	int		is_irreg_bdry_crx,
	Front		*front)
{
	COMPONENT	comp;
	CURVE		*cint[2], *cext[2], *cbint[2], *cbext[2];
	NODE		*n[2];
	POINT		*p;
	RECT_GRID	*rgr = front->rect_grid;
	double		t0[MAXD], t1store[MAXD], *t1 = &t1store[0];
	double		coords[MAXD], W[MAXD];
	double		nor[MAXD];
	boolean		sav_intrp;
	int		trk[2];
	ORIENTATION	cp_orient[2], cb_orient[2];
	ANGLE_DIRECTION	i_to_f_dir[2];
	int		n_nodes, i;
	int		dim = rgr->dim;
	double		*dt_frac = front->dt_frac;

	debug_print("vec_bdry","Entered g_vec_bdry_untangle()\n");

	if (debugging("vec_bdry"))
	{
		(void) printf("Data into g_vec_bdry_untangle()\n");
		(void) printf("cinterior\n");		print_curve(cinterior);
		(void) printf("cexterior\n");		print_curve(cexterior);
		(void) printf("bdrycurves[0]\n");	print_curve(bdrycurves[0]);
		(void) printf("bdrycurves[1]\n");	print_curve(bdrycurves[1]);
		print_orientation("cphys_orient =",cphys_orient,"\n");
		print_angle_direction("cb_to_cp_dir =",cb_to_cp_dir,"\n");
		(void) printf("is_irreg_bdry_crx = %s\n",
			(is_irreg_bdry_crx) ? "YES" : "NO");
		(void) printf("Left_state_at_node(cexterior,cphys_orient)\n");
		(*front->print_state)(Left_state_at_node(cexterior,
							 cphys_orient));
		(void) printf("Right_state_at_node(cexterior,cphys_orient)\n");
		(*front->print_state)(Right_state_at_node(cexterior,
							  cphys_orient));
		(void) printf("Left_state_at_node(cinterior,");
		(void) printf("Opposite_orient(cphys_orient))\n");
		(*front->print_state)(Left_state_at_node(cinterior,
					 Opposite_orient(cphys_orient)));
		(void) printf("Right_state_at_node(cinterior,");
		(void) printf("Opposite_orient(cphys_orient))\n");
		(*front->print_state)(Right_state_at_node(cinterior,
					  Opposite_orient(cphys_orient)));
	}

	n_nodes = find_curves_at_vec_bdry_untangle(cinterior,cexterior,
			bdrycurves,cphys_orient,cb_to_cp_dir,is_irreg_bdry_crx,
			n,cint,cext,cbint,cbext,cp_orient,cb_orient,i_to_f_dir);

	if (!is_short_curve(cext[0],Opposite_orient(cp_orient[0]),rgr,1.5))
	{
		if (debugging("vec_bdry"))
		{
			(void) printf("Exterior curve cext[0] not short\n");
			(void) printf("cext[0]");	print_curve(cext[0]);
		}
		*dt_frac = min(Min_time_step_modification_factor(front),
			       *dt_frac);
		debug_print("vec_bdry","Left g_vec_bdry_untangle()\n");
		return MODIFY_TIME_STEP_TO_UNTANGLE;
	}
	if (n_nodes != 1 &&
	(!is_short_curve(cbext[0],Opposite_orient(cb_orient[0]),rgr,1.5)))
	{
		if (debugging("vec_bdry"))
		{
			(void) printf("Exterior curve cbext[0] not short\n");
			(void) printf("cbext[0]");	print_curve(cbext[0]);
		}
		*dt_frac = min(Min_time_step_modification_factor(front),
			       *dt_frac);
		debug_print("vec_bdry","Left g_vec_bdry_untangle()\n");
		return MODIFY_TIME_STEP_TO_UNTANGLE;
	}

	delete_interior_points_of_curve(front,cext[0]);
	for (i = 0; i < n_nodes; i++)
	{
		if (!set_states_at_new_B_refl_node(n[i],cint[i],cext[i],
				cbint[i],cbext[i],i_to_f_dir[i],cp_orient[i],
				cb_orient[i],front))
		{
			(void) printf("WARNING in g_vec_bdry_untangle(), ");
			(void) printf("unable to set states at new node\n");
			debug_print("vec_bdry","Left g_vec_bdry_untangle()\n");
			return ERROR_IN_UNTANGLE;
		}
		trk[i] =  (node_type(n[i]) == B_REFLECT_NODE) ?
				track_scattered_wave(B_REFLECT_NODE,SHOCK_WAVE,
					REFLECTED,Rp_data(n[i])->state[1],
					Rp_data(n[i])->state[2],front) :
				NO;
	}

	set_new_bdry_states(n_nodes,n[0],i_to_f_dir[0],
		Opposite_orient(cb_orient[0]),cbext[0],front);


	if (n_nodes == 1 || (!trk[0]) || (!trk[1]))
	{
		comp = (curve_ang_oriented_l_to_r(i_to_f_dir[0],cp_orient[0]))
			? negative_component(cint[0]) :
			  positive_component(cint[0]);
		if (curve_ang_oriented_l_to_r(i_to_f_dir[0],cb_orient[0]))
			negative_component(cbext[0]) = comp;
		else
			positive_component(cbext[0]) = comp;
		if (debugging("vec_bdry"))
		{
			(void) printf("deleting reflected curve, trk = (%s, %s), ",
				(trk[0]) ? "YES" : "NO",
				(trk[1]) ? "YES" : "NO");
			(void) printf("n_nodes = %d\n",n_nodes);
		}
		(void) delete_curve(cext[0]);
		/*
		*  TODO:  allow for single node transition 
		*  this requires the debifurcation from REGULAR
		*  to MACH or NEUMANN node.
		*/
		if (n_nodes == 1) node_type(n[0]) = NEUMANN_NODE;
		debug_print("vec_bdry","Left g_vec_bdry_untangle()\n");
		return CURVES_UNTANGLED;
	}
	t0[0] = cos(Rp_data(n[0])->ang[2]); t0[1] = sin(Rp_data(n[0])->ang[2]);
	t1[0] = cos(Rp_data(n[1])->ang[2]); t1[1] = sin(Rp_data(n[1])->ang[2]);
	if (!intersect_ray_with_sector(n[0]->posn,n[1]->posn,
		t0,&t1,coords,dim))
	{
		(void) printf("WARNING in g_vec_bdry_untangle(), ");
		(void) printf("No intersection of rays\n");
		debug_print("vec_bdry","Left g_vec_bdry_untangle()\n");
		return ERROR_IN_UNTANGLE;
	}
	comp = new_component(NEW_COMP);
	if (curve_ang_oriented_l_to_r(i_to_f_dir[0],cp_orient[0]))
		positive_component(cexterior) = comp;
	else
		negative_component(cexterior) = comp;
	if (curve_ang_oriented_l_to_r(i_to_f_dir[0],cb_orient[0]))
		negative_component(cbext[0]) = comp;
	else
		positive_component(cbext[0]) = comp;

	wave_type(cexterior) = opposite_wave_type(wave_type(cexterior));
	sav_intrp = interpolate_intfc_states(cexterior->interface);
	interpolate_intfc_states(cexterior->interface) = YES;
	p = Point(coords);
	insert_point_adjacent_to_node(p,cexterior,POSITIVE_ORIENTATION);
	interpolate_intfc_states(cexterior->interface) = sav_intrp;
	normal(p,Hyper_surf_element(cexterior->first),
	       Hyper_surf(cexterior),nor,front);
	w_speed(Coords(p),left_state(p),right_state(p),left_state(p),
		right_state(p),W,0.0,nor,wave_type(cexterior),front);

	if (debugging("vec_bdry"))
	{
		(void) printf("Interface after g_vec_bdry_untangle()\n");
		print_interface(cinterior->interface);
	}
	debug_print("vec_bdry","Left g_vec_bdry_untangle()\n");
	return CURVES_UNTANGLED;
}		/*end g_vec_bdry_untangle*/

LOCAL	void set_new_bdry_states(
	int		n_nodes,
	NODE		*n,
	ANGLE_DIRECTION	i_to_f_dir,
	ORIENTATION	cb_orient,
	CURVE		*cbext,
	Front		*front)
{
	NODE		*oppn;
	BOND		*b;
	double		t, nor[MAXD];
	double		dist;
	Locstate	st0, st1;
	int		dim = front->rect_grid->dim;

	if (node_type(n) == NEUMANN_NODE) return;
	oppn = Node_of(cbext,Opposite_orient(cb_orient));
	if (n_nodes == 1)
	{
		if (curve_ang_oriented_l_to_r(i_to_f_dir,cb_orient))
		{
			assign_interacting_states(oppn->posn,cbext,cb_orient,
				front,Left_state_at_node(cbext,cb_orient),
				Rp_data(n)->state[2]);
		}
		else
		{
			assign_interacting_states(oppn->posn,cbext,cb_orient,
				front,Rp_data(n)->state[2],
				Right_state_at_node(cbext,cb_orient));
		}
		return;
	}
	dist = separation(n->posn,oppn->posn,dim);
	if (curve_ang_oriented_l_to_r(i_to_f_dir,cb_orient))
	{
	    /* Physical side = right */

	    st0 = right_start_state(cbext);
	    st1 = right_end_state(cbext);
	    for (b = cbext->first; b != cbext->last; b = b->next)
	    {
	    	t = separation(b->end,cbext->start->posn,dim)/dist;
	    	interpolate_states(front,1.0-t,t,Coords(cbext->start->posn),st0,
				   Coords(cbext->end->posn),st1,
				   right_state(b->end));
		normal(b->end,Hyper_surf_element(b),Hyper_surf(cbext),
		       nor,front);
		zero_normal_velocity(right_state(b->end),nor,dim);
	    }
	}
	else
	{
	    /* Physical side = left */

	    st0 = left_start_state(cbext);
	    st1 = left_end_state(cbext);
	    for (b = cbext->first; b != cbext->last; b = b->next)
	    {
	    	t = separation(b->end,cbext->start->posn,dim)/dist;
	    	interpolate_states(front,1.0-t,t,Coords(cbext->start->posn),st0,
	    		           Coords(cbext->end->posn),st1,
	    		           left_state(b->end));
	    	normal(b->end,Hyper_surf_element(b),Hyper_surf(cbext),
		       nor,front);
	    	zero_normal_velocity(left_state(b->end),nor,dim);
	    }
	}
}		/*end set_new_bdry_states*/

LOCAL	int		set_states_at_new_B_refl_node(
	NODE		*n,
	CURVE		*cint,
	CURVE		*cext,
	CURVE		*cbint,
	CURVE		*cbext,
	ANGLE_DIRECTION	i_to_f_dir,
	ORIENTATION	cp_orient,
	ORIENTATION	cb_orient,
	Front		*front)
{
	RP_DATA		*RP;
	double		tcp[MAXD], tcab[MAXD], tcbb[MAXD];

	if (!is_shock_wave(wave_type(cint)))
	{
		node_type(n) = NEUMANN_NODE;
		return YES;
	}

	RP = Rp_data(n);
	RP->ang_dir = Opposite_ang_dir(i_to_f_dir);
	if (!estimate_node_vel(cbint,cb_orient,tcab,cint,cp_orient,tcp,
		i_to_f_dir,B_REFLECT_NODE,Node_vel(n),front))
	{
		(void) printf("WARNING in set_states_at_new_B_refl_node(), ");
		(void) printf("can't find node velocity\n");
		node_type(n) = NEUMANN_NODE;
		return YES;
	}
	RP->ang[0] = angle(tcab[0],tcab[1]);
	RP->ang[1] = angle(tcp[0],tcp[1]);
	find_tangent_to_curve(n->posn,
		Bond_at_node(cbext,Opposite_orient(cb_orient)),
		cbext,Opposite_orient(cb_orient),tcbb,front);
	RP->ang[3] = angle(tcbb[0],tcbb[1]);
	if (curve_ang_oriented_l_to_r(i_to_f_dir,cp_orient))
	{
		ft_assign(RP->state[0],Right_state_at_node(cint,cp_orient),
			front->sizest);
		ft_assign(RP->state[1],Left_state_at_node(cint,cp_orient),
			front->sizest);
	}
	else
	{
		ft_assign(RP->state[0],Left_state_at_node(cint,cp_orient),
			front->sizest);
		ft_assign(RP->state[1],Right_state_at_node(cint,cp_orient),
			front->sizest);
	}
	if (!is_regular_reflection(Node_vel(n),front,RP))
	{
		/* NOTE: consider the case of an expanding shock tangling
		 * with the boundary.  This routine is called in the course
		 * of the untangle.  For weak shocks, it is possible that
		 * the node velocities satisfy the CFL condition, BUT the
		 * incident angle is too large for a regular configuration.
		 * The correct answer would be not to change the node type
		 * below and to return NO, causing a time step reduction.
		 * This messes up the rp code, though, when a shock completes
		 * its diffraction through a contact.  In that case, the shock
		 * SHOULD make a very large angle with the wall, so that
		 * the following code is the reasonable thing to do.  The
		 * expanding shock problem occurs in a relatively small
		 * region of parameter space, but this issue needs to be
		 * resolved somehow.
		 */

		(void) printf("WARNING in set_states_at_new_B_refl_node(), ");
		(void) printf("irregular boundary reflection\n");
		(void) printf("CODE NEEDED\n");
		node_type(n) = NEUMANN_NODE;
		return YES;
	}
	node_type(n) = B_REFLECT_NODE;
	set_status_at_node(cext,Opposite_orient(cp_orient),REFLECTED);
	if (curve_ang_oriented_l_to_r(i_to_f_dir,cp_orient))
	{
		ft_assign(Left_state_at_node(cext,Opposite_orient(cp_orient)),
			RP->state[1],front->sizest);
		ft_assign(Right_state_at_node(cext,Opposite_orient(cp_orient)),
			RP->state[2],front->sizest);
	}
	else
	{
		ft_assign(Right_state_at_node(cext,Opposite_orient(cp_orient)),
			RP->state[1],front->sizest);
		ft_assign(Left_state_at_node(cext,Opposite_orient(cp_orient)),
			RP->state[2],front->sizest);
	}
	if (curve_ang_oriented_l_to_r(RP->ang_dir,cb_orient))
	{
		ft_assign(Right_state_at_node(cbext,Opposite_orient(cb_orient)),
			RP->state[2],front->sizest);
	}
	else
	{
		ft_assign(Left_state_at_node(cbext,Opposite_orient(cb_orient)),
			RP->state[2],front->sizest);
	}
	return YES;
}		/*end set_states_at_new_B_refl_node*/

LOCAL	int find_curves_at_vec_bdry_untangle(
	CURVE		*cinterior,
	CURVE		*cexterior,
	CURVE		**bdrycurves,
	ORIENTATION	cphys_orient,
	ANGLE_DIRECTION	cb_to_cp_dir,
	int		is_irreg_bdry_crx,
	NODE		**n,
	CURVE		**cint,
	CURVE		**cext,
	CURVE		**cbint,
	CURVE		**cbext,
	ORIENTATION	*cp_orient,
	ORIENTATION	*cb_orient,
	ANGLE_DIRECTION	*i_to_f_dir)
{
	CURVE		**c;
	int		w_type;

	debug_print("vec_bdry","Entered find_curves_at_vec_bdry_untangle()\n");
	n[0] = Node_of(cexterior,cphys_orient);
	cint[0] = cinterior;	cext[0] = cexterior;
	cp_orient[0] = Opposite_orient(cphys_orient);
	w_type = wave_type(cinterior);
	if (((cb_to_cp_dir == CLOCKWISE) && is_forward_wave(w_type)) ||
		(cb_to_cp_dir == COUNTER_CLOCK && is_backward_wave(w_type)))
	{
		cbint[0] = bdrycurves[0];
		cb_orient[0] = NEGATIVE_ORIENTATION;
		cbext[0] = bdrycurves[1];
		i_to_f_dir[0] =
			(curve_ang_oriented_l_to_r(cb_to_cp_dir,cp_orient[0])) ?
				CLOCKWISE : COUNTER_CLOCK;
	}
	else
	{
		cbint[0] = bdrycurves[1];
		cb_orient[0] = POSITIVE_ORIENTATION;
		cbext[0] = bdrycurves[0];
		i_to_f_dir[0] =
			(curve_ang_oriented_l_to_r(cb_to_cp_dir,cp_orient[0])) ?
				COUNTER_CLOCK : CLOCKWISE;
	}

	if (debugging("vec_bdry"))
	{
		(void) printf("n[0]\n");	print_node(n[0]);
		print_orientation("cp_orient[0] =",cp_orient[0],"\n");
		print_orientation("cb_orient[0] =",cb_orient[0],"\n");
		print_angle_direction("i_to_f_dir[0] =",i_to_f_dir[0],"\n");
		(void) printf("cint[0]\n");	print_curve(cint[0]);
		(void) printf("cext[0]\n");	print_curve(cext[0]);
		(void) printf("cbint[0]\n");	print_curve(cbint[0]);
		(void) printf("cbext[0]\n");	print_curve(cbext[0]);
	}

	if (is_irreg_bdry_crx)
	{
		debug_print("vec_bdry","Left find_curves_at_vec_bdry_untangle()\n");
		return 1;
	}

	n[1] = Node_of(cexterior,Opposite_orient(cphys_orient));
	cbext[1] = cbext[0];
	cext[1] = cext[0];
	cp_orient[1] = Opposite_orient(cp_orient[0]);
	cb_orient[1] = Opposite_orient(cb_orient[0]);
	i_to_f_dir[1] = Opposite_ang_dir(i_to_f_dir[0]);
	if (cp_orient[1] == POSITIVE_ORIENTATION)
	{
		for (c = n[1]->out_curves; c && *c; c++)
		{
			if (wave_type(*c) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
			{
				cint[1] = *c;
				break;
			}
		}
	}
	else
	{
		for (c = n[1]->in_curves; c && *c; c++)
		{
			if (wave_type(*c) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
			{
				cint[1] = *c;
				break;
			}
		}
	}
	if (cb_orient[1] == POSITIVE_ORIENTATION)
	{
		for (c = n[1]->out_curves; c && *c; c++)
		{
			if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE)
			{
				cbint[1] = *c;
				break;
			}
		}
	}
	else
	{
		for (c = n[1]->in_curves; c && *c; c++)
		{
			if (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE)
			{
				cbint[1] = *c;
				break;
			}
		}
	}
	if (debugging("vec_bdry"))
	{
		(void) printf("n[1]\n");	print_node(n[1]);
		print_orientation("cp_orient[1] =",cp_orient[1],"\n");
		print_orientation("cb_orient[1] =",cb_orient[1],"\n");
		print_angle_direction("i_to_f_dir[1] =",i_to_f_dir[1],"\n");
		(void) printf("cint[1]\n");	print_curve(cint[1]);
		(void) printf("cext[1]\n");	print_curve(cext[1]);
		(void) printf("cbint[1]\n");	print_curve(cbint[1]);
		(void) printf("cbext[1]\n");	print_curve(cbext[1]);
	}
	debug_print("vec_bdry","Left find_curves_at_vec_bdry_untangle()\n");
	return 2;
}		/*end find_curves_at_vec_bdry_untangle*/

LOCAL	void set_orients_and_w_type_for_cross(
	ORIENTATION	*c_orient,
	int		*w_type)
{
	c_orient[1] = Opposite_orient(c_orient[4]);
	c_orient[3] = Opposite_orient(c_orient[0]);
	c_orient[2] = Opposite_orient(c_orient[0]);

	w_type[2] = CONTACT;
	w_type[1] = w_type[4];
	w_type[3] = w_type[0];
}		/*end set_orients_and_w_type_for_cross*/

/*
*		set_track_and_comp_list_for_cross():
*
*	Join is set TRUE if there are no tracked reflected or 
*	transmitted waves so that the two components must be 
*	later meshed together
*/

LOCAL	void set_track_and_comp_list_for_cross(
	boolean		*track,
	COMPONENT	*comp,
	ANGLE_DIRECTION	i0_to_i4_dir,
	ORIENTATION	*c_orient,
	Front		*fr,
	CURVE		*incident0,
	CURVE		*incident4,
	int		*join,
	RP_DATA		*RP)
{
	int		i, j;

	track[0] = track[4] = YES;
	track[1] = track_scattered_wave(CROSS_NODE,SHOCK_WAVE,REFLECTED,
				 RP->state[1],RP->state[2],fr);

	track[2] = track_scattered_wave(CROSS_NODE,CONTACT_WAVE,
				 SLIP,RP->state[2],RP->state[3],fr);

	track[3] = track_scattered_wave(CROSS_NODE,SHOCK_WAVE,REFLECTED,
				 RP->state[4],RP->state[3],fr);

	if (!track[1] && !track[2] && !track[3])
	{
		*join = YES;
		for (i = 2; i < 5; i++) comp[i] = comp[1];
		return;
	}
	*join = NO;

	if (curve_ang_oriented_l_to_r(i0_to_i4_dir,c_orient[0]))
	{
		comp[0] = positive_component(incident0);
		comp[1] = negative_component(incident0);
	}
	else
	{
		comp[1] = positive_component(incident0);
		comp[0] = negative_component(incident0);
	}
	comp[4] = (curve_ang_oriented_l_to_r(i0_to_i4_dir , c_orient[4])) ?
		positive_component(incident4) : negative_component(incident4);
	for (j = 2; j < 4; j++)
	{
	    if (track[j-1] == NO) comp[j] = comp[j-1];
	    else
	    {
		comp[j] = comp[4];
		for (i = j; i < 4; i++)
		{
		    if (track[i])
		    {
			comp[j] = new_component(NEW_COMP);
			break;
		    }
	        }
	    }
	}
}		/*end set_track_and_comp_lists_for_cross*/


LOCAL	void set_states_and_comps_about_cross(
	COMPONENT	*newcomp,
	ANGLE_DIRECTION	i0_to_i4_dir,
	ORIENTATION	*c_orient,
	COMPONENT	*l_comp,
	COMPONENT	*r_comp,
	Locstate	*l_st,
	Locstate	*r_st,
	RP_DATA		*RP)
{
	/* VERY similar to set_states_and_comps_about_shock_diffraction
	   Should be combined with this */

	int j;
	for (j = 0; j < 5; j++)
	{
	    if (curve_ang_oriented_l_to_r(i0_to_i4_dir, c_orient[j]))
	    {
	        l_st[j] = RP->state[(j+1)%5];
	        l_comp[j] = newcomp[(j+1)%5];
	        r_st[j] = RP->state[j];
	        r_comp[j] = newcomp[j];
	    }
	    else
	    {
	        l_st[j] = RP->state[j];
	        l_comp[j] = newcomp[j];
	        r_st[j] = RP->state[(j+1)%5];
	        r_comp[j] = newcomp[(j+1)%5];
	    }
	}
}		/*end set_states_and_comps_about_cross*/
#endif /* defined(TWOD) */
