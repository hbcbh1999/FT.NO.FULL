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
*				gbifurprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Function prototypes for the gbifur libarary
*/

#if !defined(_GBIFURPROTOS_H)
#define _GBIFURPROTOS_H

	/* gbifur EXPORTED Function Declarations */
/*	gbifur/gbifur.c */
#if defined(ONED) || defined(TWOD)
IMPORT	boolean	track_scattered_wave(int,int,int,Locstate,Locstate,Front*);
#endif /* defined(ONED) || defined(TWOD) */
#if defined(TWOD)
IMPORT	Bubble	*allocate_bubble(boolean,Front*,RP_DATA*,int);
IMPORT	boolean	g_replace_unphys_loop(NNLIST*,NNLIST**,CURVE**,Front*,
				      int,double,int);
IMPORT	int	g_B_node_bifurcation(Front*,POINTER,O_CURVE*,O_CURVE*,O_CURVE*,
				     O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
				     O_CURVE*,O_CURVE*,O_CURVE*,POINT*,
				     Locstate,Locstate,ANGLE_DIRECTION,
				     RPROBLEM**,double,double*,NODE_FLAG);
IMPORT	void	free_bubble(Bubble*,boolean);
IMPORT	void	g_identify_physical_node(NODE*);
IMPORT	void	g_reflect_node_bifurcation(Front*,Wave*,O_CURVE*,O_CURVE*,
					   O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
					   O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
					   O_CURVE*,O_CURVE*,POINT*,POINT*,
					   POINT*,BOND*,BOND*,RP_DATA*,
					   RPROBLEM*,double,double*,NODE_FLAG);
IMPORT	void	set_corner_for_bifurcation(O_CURVE*,O_CURVE*,O_CURVE*,O_CURVE*,
					   O_CURVE*,O_CURVE*,POINT*,POINT*,
					   Front*,Wave*,POINT**,BOND**,int,
					   double*,RP_DATA*,double,double*);
#endif /* defined(TWOD) */

/*	gbifur/gcapture.c */
#if defined(TWOD)
#if defined(FULL_PHYSICS)
IMPORT	boolean	snd_wave_end_propagate(Front*,Front*,NODE*,NODE*,NODE*,double);
IMPORT	int	wave_end_propagate(Front*,NODE*,NODE*,RPROBLEM**,
				   double,double*,NODE_FLAG);
IMPORT	void	g_shock_capture(Front*);
#endif /* defined(FULL_PHYSICS) */
#endif /* defined(TWOD) */

/*	gbifur/grefl.c */
#if defined(TWOD)
IMPORT	CURVE	*make_half_bow_curve(COMPONENT,COMPONENT,int,ORIENTATION,
				     NODE*,NODE*,double*,RECT_GRID*,Bubble*);
IMPORT	int	init_mach_bubble(int,double,Bubble*,int,Front*);
IMPORT	void	bubble_state(Locstate,POINT*,CURVE*,Bubble*,COMPONENT);
IMPORT	void	find_bubble_state(Locstate,double*,COMPONENT,INTERFACE*,int);
IMPORT	void	init_bubble(double*,int,double,Bubble*,Front*);
#endif /* defined(TWOD) */

/*	gbifur/grp.c */
#if defined(TWOD)
IMPORT	boolean	g_untrack_curve(O_CURVE*,O_CURVE*,COMPONENT,double,Front*,
				POINTER,RPROBLEM*,UNTRACK_FLAG);
IMPORT 	double	find_position_and_dt_of_intersection(CURVE*,ORIENTATION,
						     CURVE*,ORIENTATION,
						     double*,Front*,double);
IMPORT	int	g_2drproblem(Front*,Front*,POINTER,RPROBLEM**);
IMPORT	void	g_init_2drproblem(RPROBLEM*,Front*);
IMPORT	void	set_g_rproblem_hooks(void);
#endif /* defined(TWOD) */

/*	gbifur/gsc1.c */
#if defined(TWOD)
IMPORT	int	realign_phys_curves_at_nbdry(CURVE**,ORIENTATION*,
					     CURVE*,ORIENTATION,
					     CURVE*,ORIENTATION,
					     ANGLE_DIRECTION,POINT**,Front*);
IMPORT	int	shock_contact_rp(Front*,Wave*,RPROBLEM*);
IMPORT	void	assign_ahead_states_at_shock_diffraction(double,BOND*,CURVE*,
							 ORIENTATION,double,
							 BOND*,CURVE*,
							 ORIENTATION,COMPONENT*,
							 RP_DATA*);
IMPORT	void	attach_phys_curve_to_n_bdry(CURVE*,ORIENTATION,
					    int,POINT**,Locstate,
					    Locstate,double*,double*,double*,
					    Front*,double,double,Front*);
IMPORT	void	find_node_vel_at_rp(POINT*,double,BOND*,O_CURVE*,O_CURVE*,
				    double*,double,BOND*,O_CURVE*,O_CURVE*,
				    double*,ANGLE_DIRECTION,int,Front*,
				    Wave*,double,double*);
IMPORT	void	interpolate_state_next_to_node(CURVE*,ORIENTATION,Front*);
IMPORT	void	set_orients_and_w_type_for_shock_diffraction(int,ORIENTATION,
							     int,
							     ORIENTATION,
							     ORIENTATION*,
							     int*);
IMPORT	void	set_states_and_comps_about_shock_diffraction(COMPONENT*,
							     ANGLE_DIRECTION,
							     ORIENTATION*,
							     COMPONENT*,
							     COMPONENT*,
							     Locstate*,
							     Locstate*,
							     RP_DATA*);
IMPORT	void	set_track_and_newcomp_list_for_shock_diffraction(boolean*,
								 COMPONENT*,
								 int*,RP_DATA*,
								 boolean,Front*);
#endif /* defined(TWOD) */

/*	gbifur/gsc2.c */
#if defined(TWOD)
IMPORT	int	install_dfrctn_to_bdry(CURVE**,ORIENTATION*,
				       ANGLE_DIRECTION,int,COMPONENT*newcomp,
				       RPROBLEM*,RP_DATA*,Front*);
IMPORT	void	identify_bdry_curves_at_shock_diffraction(CURVE**,ORIENTATION*,
							  CURVE**,ORIENTATION*,
							  ANGLE_DIRECTION,
							  RPROBLEM*);
#endif /* defined(TWOD) */


/*	gbifur/gsc3.c */
#if defined(TWOD)
IMPORT	int	compute_node_velocity(int,int,Locstate,Locstate,Locstate,
				      double*,double*,double*,int,int,
				      ANGLE_DIRECTION);
IMPORT	int	diffracted_shock_exits_contact(Front*,Front*,Wave*,RPROBLEM*);
IMPORT	int	estimate_node_vel(CURVE*,ORIENTATION,double*,
				  CURVE*,ORIENTATION,double*,
				  ANGLE_DIRECTION,int,double*,Front*);
IMPORT	int	n_diffracted_shock_exits_contact(Front*,Front*,Wave*,RPROBLEM*);
IMPORT	int	overtake_diffraction_collision(Front*,Front*,Wave*,RPROBLEM*);
IMPORT	void	find_states_at_node_of_interaction(ORIENTATION,ORIENTATION,
						   ANGLE_DIRECTION,Locstate,
						   Locstate,Locstate,Locstate,
						   Locstate*,Locstate*,
						   Locstate*,Locstate*);
IMPORT  void    install_bdry_cross(CROSS*,O_CURVE*,RP_NODE*,RPROBLEM*,
				   Front*,SIDE*);
#endif /* defined(TWOD) */

#if defined(ONED)
/* guntan1d.c */
IMPORT int g_untangle_front1d(Front*,CROSS**,int);
IMPORT int g_bdry_untangle1d(Front*,CROSS**,RPROBLEM*,NODE*,int);
#endif /* defined(ONED) */

/*	gbifur/guntan2d.c */
#if defined(TWOD)
IMPORT	int	g_untangle_interior_curves(Front*,CROSS**cross,int);
IMPORT	void	g_phys_split_bdry_cross(CURVE**,CURVE**);

/*	gbifur/gvecuntan.c */
IMPORT	int	g_vec_bdry_untangle(CURVE*,CURVE*,CURVE**,
				    ORIENTATION,ANGLE_DIRECTION,int,Front*);
IMPORT	int	vector_cross_unravel(Front*,CROSS*,CROSS*,
				     ORIENTATION,ORIENTATION);
IMPORT	int	vector_overtake_unravel(Front*,CROSS*,CROSS*,
					ORIENTATION,ORIENTATION);
#endif /* defined(TWOD) */
#endif /* !defined(_GBIFURPROTOS_H) */
