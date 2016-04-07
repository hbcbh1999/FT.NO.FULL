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
*				gintfcprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the structures for gas dynamics states.
*/

#if !defined(_GINTFCPROTOS_H)

	/* gintfc EXPORTED Function Declarations */

/*	gintfcpert.c */
#if defined(TWOD) || defined(THREED)
IMPORT	SINE_PERT *alloc_and_copy_sine_pert_structure(SINE_PERT*,Front*);
IMPORT	SINE_PERT *alloc_sine_pert_structure(int,Front*);
IMPORT	double	  pert_interface(SINE_PERT*,double*,double,int);
#endif /* defined(TWOD) || defined(THREED) */

/*	gintfc/guserintfc.c */
IMPORT	G_USER_INTERFACE	*g_user_hook(int);
IMPORT	boolean	g_form_subintfc_via_communication2d(Front*);
#if defined(USE_OVERTURE)
IMPORT  void    g_assembly_fine_patch_fronts_to_one(Front**,Front*);
IMPORT  boolean    g_form_patch_subintfc_via_cut2d(Front*);
IMPORT  boolean    g_form_patch_subintfc_2d(Front*,COMPONENT);
#endif /* if defined(USE_OVERTURE) */

IMPORT	double	prompt_for_surface_tension(const int,const char*);
IMPORT	void	constant_density_stratified_state(Locstate,double,
						  double,Locstate);
IMPORT	void	g_preserve_user_hooks(int,PRESERVE_USER_HOOKS);
IMPORT	void	g_set_interface_hooks(int,INIT_DATA*);
IMPORT	void	g_set_stratified_state(void(*)(Locstate,double,double,Locstate),
				       const char*);
IMPORT  void    g_user_fprint_node(FILE*,NODE*);
IMPORT	void	isentropic_stratified_state(Locstate,double,double,Locstate);
IMPORT	void	isothermal_stratified_state(Locstate,double,double,Locstate);
IMPORT	void	stratified_state(INTERFACE*,Locstate,double,double,Locstate);

/*	gintfc/gtop.c */
IMPORT	void	g_reflect_point(POINT*,double*,double*,INTERFACE*);
#if defined(TWOD)
IMPORT	void	g_delete_small_loops(Front*);
IMPORT	void	g_reflect_node2d(NODE*,double*,double*);
IMPORT  boolean    g_redist_on_cc_node(Front*);
IMPORT  boolean    g_parallel_redist_on_cc_node(Front*);
IMPORT  void    merge_slip_with_other_in_parallel(Front*,CURVE*);
#endif /* defined(TWOD) */

/*	gintfc/gtypes.c */
IMPORT	boolean	g_is_correspondence_possible(HYPER_SURF*,HYPER_SURF*,
					  HYPER_SURF_BDRY**,HYPER_SURF_BDRY**);
IMPORT	int	opposite_wave_type(int);
IMPORT	void	set_fixed_status(CURVE*,ORIENTATION);
IMPORT	void	set_incident_status(CURVE*,ORIENTATION);
IMPORT	void	set_start_end_status(Front*);
IMPORT	void	set_status_at_node(CURVE*,ORIENTATION,int);

/*	gintfc/guserhooks.c */
IMPORT	void	w_speed(double*,Locstate,Locstate,Locstate,Locstate,double*,
			double,double*,int,Front*);
IMPORT	void	unsplit_w_speed2d(USWSSten2d*,Locstate,Locstate,double*);
IMPORT	void	npt_w_speed(WSSten*,Locstate,Locstate,double*);

#if defined(THREED)
/*	gintfc/gcheck3d.c */
IMPORT	boolean	g_consistent_interface(INTERFACE*);
#endif /* defined(THREED) */

#endif /* !defined(_GINTFCPROTOS_H) */
