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
*				gnodeprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GNODEPROTOS_H)
#define _GNODEPROTOS_H

	/* gnode EXPORTED Function Declarations */

#if defined(TWOD)

/*	gnode/gRP.c */
IMPORT	RP_DATA	*allocate_RP_DATA_structure(size_t,boolean,int);
IMPORT	void	copy_RP_DATA_structure(RP_DATA*,RP_DATA*);
IMPORT	void	free_RP_DATA_structure(RP_DATA*);

/*	gnode/gbnode.c */
IMPORT	int	B_reflect_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
					 double,double*,NODE_FLAG);
IMPORT	int	attached_b_node_propagate(Front*,Wave*,NODE*,NODE*,
					  RPROBLEM**,double,double*,NODE_FLAG);
IMPORT	int	is_bow_shock_attached(Locstate,double,Locstate,double*);
IMPORT	int	is_regular_reflection(double*,Front*,RP_DATA*);

/*	gnode/gccnode.c */
IMPORT	int	cc_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
				  double,double*,NODE_FLAG);
IMPORT	void	set_cc_normal_function(NORMAL_FUNCTION*);
IMPORT	void	set_cc_tangent_function(TANGENT_FUNCTION*);

/*      gnode/gmdnode.c */
IMPORT	int	reg_to_mach_diff_reconfigure(O_CURVE**,O_CURVE**,NODE*);

/*	gnode/gnode.c */
IMPORT	boolean	is_bad_state_at_node(const char*,NODE*);
IMPORT	int	g_node_propagate(Front*,POINTER,NODE*,NODE*,RPROBLEM**,
				 double,double*,NODE_FLAG,POINTER);
IMPORT	void	node_warning(const char*,const char*,const char*,NODE_FLAG);

/*	gnode/gnodesub.c */
IMPORT ADJUST_ANGLE_VALUE g_adjust_angle_len(NODE*,O_CURVE*,O_CURVE*,Locstate,
					     Locstate,double*,double,double*,
					     double*,Front*,Wave*);
IMPORT	double	sonic_radius(Locstate,double*,double,int);
IMPORT	boolean	g_check_delete_redundant_node(NODE*,CURVE*,CURVE*);
IMPORT	int	adjust_angle_at_node(NODE*,O_CURVE*,O_CURVE*,Locstate,Locstate,
				     double*,double,double,Front*,Wave*);
IMPORT	ANGLE_DIRECTION	g_find_i_to_prop_dir(Front*,POINTER,NODE*,CURVE*,
				     	     ORIENTATION,double,
				     	     COMPONENT*,POINT*,double*);
IMPORT	int	modify_curves_at_node(POINT*,BOND**,NODE*,O_CURVE**,O_CURVE**,
				      int,boolean*,Front*,Wave*,double,
				      RP_DATA*,NODE_FLAG);
IMPORT	int	propagate_curve_near_node(NODE*,NODE*,O_CURVE*,O_CURVE*,BOND*,
					  Locstate,Locstate,boolean,double,
					  Front*,Wave*,double,NODE_FLAG);
IMPORT	void	find_adjacent_curves(O_CURVE*,ANGLE_DIRECTION*,O_CURVE*,
				     O_CURVE*,COMPONENT*);
IMPORT	void	find_curve_with_status(NODE*,CURVE**,ORIENTATION*,int);
IMPORT	void	identify_curves_with_status(NODE*,O_CURVE*,O_CURVE*,int);

/*	gnode/gpcnode.c */
IMPORT	int	precursor_shock_rr_propagate(Front*,Wave*,NODE*,NODE*,
					     RPROBLEM**,double,double*,NODE_FLAG,
					     O_CURVE**,O_CURVE**,POINT*,
					     RP_DATA*,BOND**);

/*	gnode/gsndnode.c */
IMPORT	int	g_snd_node_propagate(Front*,Front*,POINTER,INTERFACE*,
				     NODE*,NODE*,double);
/*	gnode/gssnode.c */
IMPORT	int	Mach_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
				    double,double*,NODE_FLAG);
IMPORT	int	cross_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
				     double,double*,NODE_FLAG);
IMPORT	int	overtake_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
					double,double*,NODE_FLAG);
IMPORT	void	ramp_reflection_corner_posn(double*,int,int);
IMPORT	void	tangent_at_degenerate_node(CURVE*,ORIENTATION,
					   CURVE*,ORIENTATION,
					   double*,Front*);

/*	gnode/gssnsts.c */
IMPORT	int	find_Mach_node_states(NODE*,O_CURVE*,O_CURVE*,O_CURVE*,
				      O_CURVE*,O_CURVE*,BOND**,BOND**,POINT*,
				      Wave*,Front*,RPROBLEM**,double,double*,
				      RP_DATA*,NODE_FLAG);
IMPORT	int	find_cross_node_states(double*,RP_DATA*,boolean*,boolean);
IMPORT	int	find_overtake_node_states(double*,RP_DATA*,boolean*,boolean);
IMPORT	int	find_steady_ahead_speed(ANGLE_DIRECTION,Locstate,Locstate,
					Locstate,double,double*,double*,
					double*,double*);
IMPORT	void	set_temp_mnode_normal(NORMAL_FUNCTION*);
#endif /* defined(TWOD) */

#if defined(FULL_PHYSICS) && defined(TWOD)
/*	gnode/ganom.c */
IMPORT	int	anomalous_reflection_propagate(Front*,Wave*,NODE*,NODE*,
					       RPROBLEM**,double,double*,
					       NODE_FLAG,O_CURVE**,O_CURVE**,
					       boolean*,RP_DATA*,POINT*,double*,
					       BOND**);
IMPORT	int	refl_curve_overtakes_incident_shock(O_CURVE**,O_CURVE**,BOND**,
						    POINT*,Front*,Wave*,
						    RPROBLEM**,double,double*,
						    int,ANGLE_DIRECTION,
						    NODE_FLAG);

/*	gnode/gscnode.c */
IMPORT	ANGLE_DIRECTION	incident_shock_orientation(int,ORIENTATION);
IMPORT	boolean	curves_at_shock_diffraction(NODE*,O_CURVE**,ANGLE_DIRECTION*,
					    boolean);
IMPORT	int	diffraction_node_propagate(Front*,Wave*,NODE*,NODE*,
					   RPROBLEM**,double,double*,NODE_FLAG,
					   POINTER);
IMPORT	int	modify_diffraction_node(POINT*,NODE*,NODE*,O_CURVE**,O_CURVE**,
					BOND**,Front*,Wave*,double,
					RP_DATA*,NODE_FLAG);
IMPORT	int	trans_node_parameter(void);
IMPORT	int	set_trans_node_parameter(int);
IMPORT	int	sonic_incident_shock_at_diffraction(NODE*,O_CURVE**,O_CURVE**,
						    BOND**,POINT*,double*,
						    RP_DATA*,Front*,Wave*,
						    RPROBLEM**,double,double*,
						    NODE_FLAG);
IMPORT	int	transmission_node_propagate(Front*,Wave*,NODE*,NODE*,
					    RPROBLEM**,double,double*,NODE_FLAG);
IMPORT	void	set_temp_tnode_normal(NORMAL_FUNCTION*);

/*	gnode/gscnsts.c */
IMPORT	int	find_transmission_node_states(double*,double**,RP_DATA*,int,
					      int,NODE_FLAG);
IMPORT	int	is_regular_diffraction_node(double*,double*,double*,double**,
					    RP_DATA*,RP_DATA*,boolean*,Front*,
					    int,NODE_FLAG);
IMPORT	int	is_small_inc_ang_reg_diff_node(double**,RP_DATA*,boolean*,NODE_FLAG);

#endif /* defined(FULL_PHYSICS) && defined(TWOD) */

/*	gnode/gcurve.c */
#if defined(THREED)
IMPORT	void	g_curve_propagate_3d(Front*,POINTER,CURVE*,CURVE*,double);
#endif /* defined(THREED) */

#endif /* !defined(_GNODEPROTOS_H) */
