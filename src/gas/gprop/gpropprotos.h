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
*				gpropprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GPROPPROTOS_H)
#define _GPROPPROTOS_H

	/* gprop EXPORTED Function Declarations */
/*	gprop/gipolar.c */
IMPORT	int	i_polar_2(Locstate,Locstate,Locstate,Locstate,ANGLE_DIRECTION,
			  double*,double,double,double*,double*,double);
IMPORT	int	intersection_of_two_shock_polars(Locstate,Locstate,double*,
						 double*,double*,double*,
						 boolean,boolean,
						 RIEMANN_SOLVER_WAVE_TYPE*,
						 RIEMANN_SOLVER_WAVE_TYPE*);

/*	gprop/gpolar.c */
IMPORT	double	max_behind_shock_pr(double,Locstate);
IMPORT	int	prandtl_meyer_wave(Locstate,double,boolean,double*,Locstate,
				   double*,double*,double*);
IMPORT	int	s_polar_2(Locstate,int,int,double,double*,Locstate,double*);
IMPORT	int	s_polar_3(Locstate,int,double,int,int,double*,Locstate,
			  double*,double*);
IMPORT	int	s_polar_4(int,double,double*,double*,Locstate,Locstate,int);
IMPORT	int	s_polar_5(Locstate,int,double,int,int,double*,Locstate,
			  double*,double*);

/*	gprop/gprop.c */
IMPORT	USWSSten2d	*AllocUSWSSten2d(Front*,int,int);
IMPORT	WSSten	*g_AllocWSSten(int,int,Front*);
IMPORT	boolean	passive_point_propagate(int,POINT*,POINT*,double*,INTERFACE*);
IMPORT	double	set_pjump_at_wave(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				  Front*,double*);
IMPORT	void	g_pt_prop_by_w_speed(Front*,POINT*,POINT*,
				     HYPER_SURF_ELEMENT*,HYPER_SURF*,
				     double,double*);
IMPORT	void	include_source(double*,Locstate,double,double,double*,double*,double*,
			       RECT_GRID*,double,int);
IMPORT  void    include_source_1(double*,Locstate,double,double,double*,double*,double*,
                               RECT_GRID*,double*,int);
IMPORT	void	print_point_propagate_data(POINT*,HYPER_SURF_ELEMENT*,
					   HYPER_SURF*,int);
IMPORT	void	reflect_wssten(WSSten*,SIDE,Front*);
IMPORT	void	set_point_propagate(Front*,boolean);
IMPORT	void	states_near_location(WSSten*,double*,double*,COMPONENT,COMPONENT,
				     Locstate,Locstate);
/*#bjet2 */
IMPORT	void	g_point_propagate_along_wall(Front*,POINTER,POINT*,BOND*,CURVE*,
                                  HYPER_SURF_ELEMENT*,HYPER_SURF*,POINT*,double,double*);
IMPORT	void	pseudo_g_point_propagate_along_wall(Front*,POINTER,POINT*,BOND*,CURVE*,
                                  HYPER_SURF_ELEMENT*,HYPER_SURF*,POINT*,double,double*);
IMPORT	void	g_compute_force_and_torque(Front*,CURVE*,double,double*,double*);


IMPORT	void	g_unsplit_w_speed2d(USWSSten2d*,Locstate,Locstate,double*);
IMPORT	void	set_USWSSten2d_off_front_states(USWSSten2d*,HYPER_SURF*,
						Front*,Wave*);
IMPORT	void	set_uswssten_geometry_and_center(USWSSten2d*,POINT*,
						 HYPER_SURF_ELEMENT*,
						 HYPER_SURF*,Front*,double);
IMPORT	void	g_curve_propagate2d(Front*,POINTER,CURVE*,CURVE*,double);

/*	gprop/griemann.c */
IMPORT	boolean	oned_state_in_rarefaction_fan(double,double,Locstate,Locstate,
					      Locstate,int,double*,int);
IMPORT	boolean	find_mid_state(Locstate,Locstate,double,double*,double*,double*,
			       double*,double*,double*,RIEMANN_SOLVER_WAVE_TYPE*,
			       RIEMANN_SOLVER_WAVE_TYPE*);
IMPORT	int	onedrsoln(double,Locstate,Locstate,Locstate,double*,int);
IMPORT	void	riemann_solution(double,const double*,Locstate,
                                 Locstate,Locstate,int);
#if defined(DEBUG_GRIEMANN)
IMPORT	void	set_debug_find_mid_state(boolean);
IMPORT	void	set_debug_riem_sol(boolean);
#endif /* defined(DEBUG_GRIEMANN) */
#if defined(COMBUSTION_CODE)

/*	griemcombst.c */
IMPORT	double	CJ_det(Locstate,int,Locstate,int);
IMPORT	boolean	combust_find_mid_state(Locstate,Locstate,double*,double*,double*,
				       double*,double*,double*,
				       RIEMANN_SOLVER_WAVE_TYPE*,
				       RIEMANN_SOLVER_WAVE_TYPE*);
IMPORT	int	combust_onedrsoln(double,Locstate,Locstate,Locstate,double*,int);
#endif /* defined(COMBUSTION_CODE) */

/*	gprop/gwspeed.c */
#if defined(DEBUG_W_SPEED)
IMPORT	void	g_set_debug_w_speed(boolean);
#endif /* defined(DEBUG_W_SPEED) */

IMPORT	void	g_npt_w_speed(WSSten*,Locstate,Locstate,double*);
IMPORT	void	g_w_speed(double*,Locstate,Locstate,Locstate,Locstate,
			  double*,double,double*,int,Front*);
IMPORT	void	print_g_npt_w_speed_opts(NptWSpeedOpts*);
IMPORT	void	set_npt_wspeed_options(NptWSpeedOpts*);
IMPORT	void	state_behind_sound_wave(Locstate,Locstate,double*,double*,double,
					double,double,double,int,int,
					RIEMANN_SOLVER_WAVE_TYPE,int);
IMPORT  void    parab_nor_solver(WSSten*,Locstate,Locstate);
#endif /* !defined(_GPROPPROTOS_H) */
