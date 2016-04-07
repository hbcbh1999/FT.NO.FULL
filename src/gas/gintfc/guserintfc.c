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
*				guserintfc.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*			User Supplied Operations
*	for gas dynamical specific interface operations.
*/

#include <gdecs/gdecs.h>

	/* LOCAL function declarations */
LOCAL	HYPER_SURF      *g_make_hypersurface(COMPONENT,COMPONENT);
LOCAL	INTERFACE	*g_copy_interface(INTERFACE*);
LOCAL	INTERFACE	*g_receive_interface(int);
LOCAL	POINT	*g_Point(double*);
LOCAL	void	g_reconstruct_interface_pointers(INTERFACE*,struct Table*,
						 POINTER*,POINTER*);
LOCAL	void	g_reconstruct_point_pointers(POINT*,INTERFACE*,INTERFACE*,
					     POINTER*,POINTER*,int);
LOCAL	void 	g_send_interface(INTERFACE*,int);
LOCAL	void	g_fset_hyper_surf_color(FILE*,HYPER_SURF*);
LOCAL	void	g_user_copy_hyper_surf(HYPER_SURF*,HYPER_SURF*);
LOCAL	void	g_user_fprint_interface(FILE*,INTERFACE*);
LOCAL	void	g_user_make_interface(INTERFACE*);
LOCAL	void	reconstruct_gas_params(Locstate,INTERFACE*,INTERFACE*,POINTER*,
				       POINTER*,int);
LOCAL	void	set_params_list_for_interface(INTERFACE*);

LOCAL	CURVE	*g_copy_curve(CURVE*,NODE*,NODE*);
LOCAL	CURVE	*g_make_curve(COMPONENT,COMPONENT,NODE*,NODE*);
LOCAL	NODE	*g_copy_node(NODE*);
LOCAL	NODE	*g_make_node(POINT*);
LOCAL	boolean	g_user_join_curves(CURVE*,CURVE*,CURVE*);
LOCAL	boolean	g_user_split_curve(int,POINT*,BOND*,CURVE*,CURVE**);
LOCAL	void	g_invert_curve(CURVE*);
LOCAL	void	g_reconstruct_bond_pointers(BOND*,INTERFACE*,INTERFACE*,
					    POINTER*,POINTER*,int);
LOCAL	void	g_reconstruct_curve_pointers(CURVE*,INTERFACE*,INTERFACE*,
					     POINTER*,POINTER*,int);
LOCAL	void	g_reconstruct_node_pointers(NODE*,INTERFACE*,INTERFACE*,
					    POINTER*,POINTER*,int);
LOCAL	void	g_reverse_curve(CURVE*);
LOCAL	void	g_user_fprint_curve(FILE*,CURVE*);
LOCAL	void	g_user_fprint_point(FILE*,POINT*);
LOCAL	CURVE	*g_attach_curve_to_node(CURVE*,POINT*,BOND*,NODE*);
LOCAL	void	coords_on_axis(int,double,BOND*,double*);
LOCAL	void	states_on_axis(int,int,POINT*,BOND*,CURVE*,Front*);
LOCAL	SURFACE	*g_copy_surface(SURFACE*,CURVE**,CURVE**,boolean);
LOCAL	SURFACE *g_join_surfaces(CURVE*);
LOCAL	void	g_user_fprint_surface(FILE*,SURFACE*);
LOCAL   SURFACE *g_detach_one_surface(SURFACE *);
LOCAL   void  	g_check_print_intfc(const char*,const char*,char,INTERFACE*,
				int,int,boolean);
LOCAL   void  g_print_wall_crx(const char*,int*,int,int,CRXING*);
LOCAL   void  g_print_wall_curve_crx(const char*,int*,int,int,CRXING*);
LOCAL   void  g_print_wall_curve_crx0(const char*,POINT *, int,CRXING*);



EXPORT	void g_set_interface_hooks(
	int		dim,
	INIT_DATA       *init)
{
	I_USER_INTERFACE	*iuh = i_user_hook(dim);
	F_USER_INTERFACE	*fuh = f_user_hook(dim);

	d_set_interface_hooks(dim,init);
	iuh->size_interface = sizeof(G_INTERFACE);
	iuh->size_hyper_surf = sizeof(G_HYPER_SURF);
	iuh->size_curve = sizeof(G_CURVE);
	iuh->size_node = sizeof(G_NODE);
	iuh->_user_make_interface = g_user_make_interface;
	iuh->_copy_interface = g_copy_interface;
	iuh->_receive_interface = g_receive_interface;
	iuh->_user_fprint_interface = g_user_fprint_interface;
	iuh->_Point = g_Point;
	iuh->_send_interface = g_send_interface;
	iuh->_reconstruct_interface_pointers = g_reconstruct_interface_pointers;
	iuh->_reconstruct_point_pointers = g_reconstruct_point_pointers;
	fuh->_wave_type_as_string = g_wave_type_as_string;
	fuh->_fprint_hsbdry_type = g_fprint_hsbdry_type;
	fuh->_fprint_state_data = g_fprint_state_data;
	iuh->_fset_hyper_surf_color = g_fset_hyper_surf_color;
	iuh->_make_hypersurface = g_make_hypersurface;
	iuh->_user_copy_hyper_surf = g_user_copy_hyper_surf;
	switch (dim)
	{
	case 1:
	    iuh->_user_fprint_point = g_user_fprint_point;
	    iuh->_reflect_point = g_reflect_point;
	    break;
	case 2:
	    iuh->_make_node = g_make_node;
	    iuh->_copy_node = g_copy_node;
	    iuh->_user_fprint_node = g_user_fprint_node;
	    iuh->_make_curve = g_make_curve;
	    iuh->_copy_curve = g_copy_curve;
	    iuh->_user_fprint_curve = g_user_fprint_curve;
	    iuh->_user_split_curve = g_user_split_curve;
	    iuh->_user_join_curves = g_user_join_curves;
	    iuh->_reconstruct_node_pointers = g_reconstruct_node_pointers;
	    iuh->_reconstruct_bond_pointers = g_reconstruct_bond_pointers;
	    iuh->_reconstruct_curve_pointers = g_reconstruct_curve_pointers;
	    iuh->_invert_curve = g_invert_curve;
	    iuh->_reverse_curve = g_reverse_curve;
	    iuh->_attach_curve_to_node = g_attach_curve_to_node;
	    iuh->_reflect_node = g_reflect_node2d;
	    break;
	case 3:
	    iuh->_make_node = g_make_node;
	    iuh->_copy_node = g_copy_node;
	    iuh->_user_fprint_node = g_user_fprint_node;
	    iuh->_make_curve = g_make_curve;
	    iuh->_copy_curve = g_copy_curve;
	    iuh->_user_fprint_curve = g_user_fprint_curve;
	    iuh->_user_split_curve = g_user_split_curve;
	    iuh->_user_join_curves = g_user_join_curves;
	    iuh->_reconstruct_node_pointers = g_reconstruct_node_pointers;
	    iuh->_reconstruct_bond_pointers = g_reconstruct_bond_pointers;
	    iuh->_reconstruct_curve_pointers = g_reconstruct_curve_pointers;
	    iuh->_invert_curve = g_invert_curve;
	    iuh->_reverse_curve = g_reverse_curve;
	    iuh->_copy_surface = g_copy_surface;
	    iuh->_join_surfaces = g_join_surfaces;
	    iuh->_consistent_interface = g_consistent_interface;
	    iuh->_user_fprint_surface = g_user_fprint_surface;
	    iuh->_detach_one_surface = g_detach_one_surface;
	    iuh->_check_print_intfc = g_check_print_intfc;
	    iuh->_print_wall_crx = g_print_wall_crx;
	    iuh->_print_wall_curve_crx = g_print_wall_curve_crx;
	    iuh->_print_wall_curve_crx0 = g_print_wall_curve_crx0;
	    break;
	}
}		/*end g_set_interface_hooks*/

EXPORT	G_USER_INTERFACE *g_user_hook(
	int		dim)
{
	static G_USER_INTERFACE Guser_hooks[3];
	static boolean first = YES;

	if (first == YES)
	{
	    /* Set default values for Guser_hooks fields*/
	    first = NO;

	    Guser_hooks[0]._intfc_type = PHYSICAL_INTERFACE;
	    Guser_hooks[0].num_params = 0;
	    Guser_hooks[0].params_list = NULL;
	    Guser_hooks[0]._stratified_state = isothermal_stratified_state;
	    Guser_hooks[0].stratified_state_name =
		"isothermal_stratified_state";
	    Guser_hooks[0]._w_speed = g_w_speed;
	    Guser_hooks[0]._npt_w_speed = g_npt_w_speed;
	    Guser_hooks[0]._unsplit_w_speed2d = NULL;
	    Guser_hooks[0]._ContactWallNodeParams.wall_bond_len = 0.0;
	    Guser_hooks[0]._ContactWallNodeParams.first_adjust_time = 0.0;
	    Guser_hooks[0]._ContactWallNodeParams.first_adjust_step = 0;
	    Guser_hooks[0]._ContactWallNodeParams.adjust = NO;

	    Guser_hooks[1]._intfc_type = PHYSICAL_INTERFACE;
	    Guser_hooks[1].num_params = 0;
	    Guser_hooks[1].params_list = NULL;
	    Guser_hooks[1]._stratified_state = isothermal_stratified_state;
	    Guser_hooks[1].stratified_state_name =
		"isothermal_stratified_state";
	    Guser_hooks[1]._w_speed = g_w_speed;
	    Guser_hooks[1]._npt_w_speed = g_npt_w_speed;
	    Guser_hooks[1]._unsplit_w_speed2d = g_unsplit_w_speed2d;
	    Guser_hooks[1]._ContactWallNodeParams.wall_bond_len = 0.5;
	    Guser_hooks[1]._ContactWallNodeParams.first_adjust_time = -1.0;
	    Guser_hooks[1]._ContactWallNodeParams.first_adjust_step = -1;
	    Guser_hooks[1]._ContactWallNodeParams.adjust = YES;

	    Guser_hooks[2]._intfc_type = PHYSICAL_INTERFACE;
	    Guser_hooks[2].num_params = 0;
	    Guser_hooks[2].params_list = NULL;
	    Guser_hooks[2]._stratified_state = isothermal_stratified_state;
	    Guser_hooks[2].stratified_state_name =
		"isothermal_stratified_state";
	    Guser_hooks[2]._w_speed = g_w_speed;
	    Guser_hooks[2]._npt_w_speed = g_npt_w_speed;
	    Guser_hooks[2]._unsplit_w_speed2d = NULL;
	    Guser_hooks[2]._ContactWallNodeParams.wall_bond_len = 0.0;
	    Guser_hooks[2]._ContactWallNodeParams.first_adjust_time = 0.0;
	    Guser_hooks[2]._ContactWallNodeParams.first_adjust_step = 0;
	    Guser_hooks[1]._ContactWallNodeParams.adjust = NO;
	}
	if (dim < 1 || dim > 3)
	{
	    screen("ERROR in g_user_hook(), invalid dim %d\n",dim);
	    clean_up(ERROR);
	    return NULL;
	}
	else
	    return Guser_hooks + dim - 1;
}		/*end g_user_hook*/

EXPORT	void	g_preserve_user_hooks(
	int                     dim,
	PRESERVE_USER_HOOKS	flag)
{
	G_USER_INTERFACE        *guh;
	static G_USER_INTERFACE Sav_guh;
	static G_USER_INTERFACE *sav_guh = NULL;

	f_preserve_user_hooks(dim,flag);
	switch (flag)
	{
	case SAVE_HOOKS:
	    if (sav_guh != NULL)
	    {
		screen("ERROR in g_preserve_user_hooks(), "
		       "attempt to save without prior restore\n");
		clean_up(ERROR);
	    }
	    guh = g_user_hook(dim);
	    sav_guh = &Sav_guh;
	    *sav_guh = *guh;
	    break;
	case RESTORE_HOOKS:
	    if (sav_guh == NULL)
	    {
		screen("ERROR in g_preserve_user_hooks(), "
		       "attempt to restore without prior save\n");
		clean_up(ERROR);
	    }
	    guh = g_user_hook(dim);
	    *guh = *sav_guh;
	    sav_guh = NULL;
	    break;
	}
}		/*end g_preserve_user_hooks*/


EXPORT	double prompt_for_surface_tension(
	const int  w_type,
	const char *mesg)
{
	double		st;
	char		s[Gets_BUF_SIZE];

	if (!is_scalar_wave(w_type) && w_type != RIEMANN_PROBLEM_WAVE)
	    return 0.0;
	st = 0.0;
	screen("Enter the surface tension %s(dflt = 0): ",
		(mesg != NULL) ? mesg : "\0");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&st);
	}
	if (st < 0.0)
	{
	    screen("\nERROR in prompt_for_surface_tension(), "
	           "Surface tension must be positive!!!\n");
	    clean_up(ERROR);
	}
	return st;
}		/*end prompt_for_surface_tension*/

EXPORT	void	stratified_state(
	INTERFACE	*intfc,
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	if (intfc == NULL || g_user_interface(intfc)._stratified_state == NULL)
	{
	    zero_state_velocity(ans,intfc->dim);
	    return;
	}
	(*g_user_interface(intfc)._stratified_state)(ans,dz,gz,ref_st);
}		/*end stratified_state*/

EXPORT	void	g_set_stratified_state(
	void       (*fpointer)(Locstate,double,double,Locstate),
	const char *fname)
{
	G_USER_INTERFACE *guh;
	struct Table	*T;
	int	i;

	for (i = 0; i < MAXD; ++i)
	{
	    guh = g_user_hook(i+1);
	    guh->_stratified_state = fpointer;
	    guh->stratified_state_name = fname;
	}

	for (T = interface_table_list(); T != NULL; T = T->next)
	{
	    if (interface_type(T->interface) != PHYSICAL_INTERFACE)
	    	continue;
	    g_user_interface(T->interface)._stratified_state = fpointer;
	    g_user_interface(T->interface).stratified_state_name = fname;
	}
}		/*end g_set_stratified_state*/


/*
*		isothermal_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isothermal one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_ref, rho(0) = rho_ref, T = T_ref.
*/

EXPORT	void	isothermal_stratified_state(
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	compute_isothermal_stratified_state(ans,dz,gz,ref_st);
}		/*end isothermal_stratified_state*/

/*
*		isentropic_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isentropic one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_ref, rho(0) = rho_ref, S = S_ref.
*/

EXPORT	void	isentropic_stratified_state(
	Locstate	ans,
	double		dz,	/* distance from reference position */
	double		gz,	/* gravity */
	Locstate	ref_st)
{
	compute_isentropic_stratified_state(ans,dz,gz,ref_st);
}		/*end isentropic_stratified_state*/


EXPORT  void    constant_density_stratified_state(
        Locstate        ans,
        double           dz,     /* distance from reference position */
        double           gz,     /* gravity */
        Locstate        ref_st)
{
	compute_constant_density_stratified_state(ans,dz,gz,ref_st);
}               /*end constant_density_stratified_state*/


LOCAL	void g_user_make_interface(
	INTERFACE	*intfc)
{
	G_USER_INTERFACE *guh = g_user_hook(intfc->dim);

	f_user_make_interface(intfc);

	g_user_interface(intfc) = *guh;
	num_layers(intfc) = 0;
}		/*end g_user_make_interface*/

LOCAL	HYPER_SURF *g_make_hypersurface(
	COMPONENT neg_comp,
	COMPONENT pos_comp)
{
	HYPER_SURF *hs;
	hs = f_make_hypersurface(neg_comp,pos_comp);
	no_slip(hs) = NO;
	adherence_coeff(hs) = 0.0;
	return hs;
}		/*end g_make_hypersurface*/

LOCAL	void g_user_copy_hyper_surf(
	HYPER_SURF	*new_hs,
	HYPER_SURF	*old_hs)
{
	f_user_copy_hyper_surf(new_hs,old_hs);
	no_slip(new_hs) = no_slip(old_hs);
	adherence_coeff(new_hs) = adherence_coeff(old_hs);
}		/*end g_user_copy_hyper_surf*/


LOCAL	INTERFACE *g_receive_interface(
	int		src_id)
{
	Gas_param	**new_list;
	int		num_params, i;

	set_size_of_intfc_state(g_sizest());
	
	num_params = return_params_list(&new_list);
	return i_receive_interface(src_id);
}		/*end g_receive_interface*/

LOCAL	INTERFACE *g_copy_interface(
	INTERFACE	*intfc)
{
	INTERFACE	*new_intfc = f_copy_interface(intfc);

	if (new_intfc == NULL) return new_intfc;

	g_user_interface(new_intfc) = g_user_interface(intfc);
	num_layers(new_intfc) = num_layers(intfc);
	return new_intfc;
}		/*end g_copy_interface*/


LOCAL	POINT *g_Point(
	double		*coords)
{
	POINT		*p;

	if ((p = f_Point(coords)) == NULL) return NULL;

	return p;
}		/*end g_Point*/

LOCAL	void g_user_fprint_interface(
	FILE		*file,
	INTERFACE	*infc)
{
	f_user_fprint_interface(file,infc);
	(void) foutput(file);
	(void) fprintf(file,"INTERFACE TYPE FOR INTERFACE %llu\n",
		       interface_number(infc));
	switch (interface_type(infc))
	{
	case PHYSICAL_INTERFACE:
	    (void) fprintf(file,"Interface type = %d PHYSICAL_INTERFACE\n",
	    	           interface_type(infc));
	    if (infc->hss != NULL || infc->dim == 3)
	    {
	    	g_fprint_Gas_param_list(file,infc);
	    	/*g_fprint_Dirichlet_bdry_states(file,infc); */
	    }
	    g_fprint_RP_DATA_at_nodes(file,infc);
	    (void) foutput(file);
	    (void) fprintf(file,"STRATIFIED STATE FUNCTION\n");
	    (void) fprintf(file,"stratified_state = %s\n",
	    	           g_user_interface(infc).stratified_state_name);
	    g_fprint_ContactWallNodeParams(file,infc);
	    break;
	case EOS_INTERFACE:
	    (void) fprintf(file,"Interface type = %d EOS_INTERFACE\n",
	    	           interface_type(infc));
	    break;
	default:
	    (void) fprintf(file,"Interface type = %d -- ** UNKNOWN **\n",
	    	           interface_type(infc));
	}
}		/*end g_user_fprint_interface*/


LOCAL	void	g_send_interface(
	INTERFACE       *intfc,
	int             dst_id)
{
	Gas_param **sav_list = gas_params_list(intfc);

	set_params_list_for_interface(intfc);
	
	i_send_interface(intfc,dst_id);
	gas_params_list(intfc) = sav_list;
}


/*
*			set_params_list_for_interface():
*
*	Copy params list into g_user_interface
*/

LOCAL	void	set_params_list_for_interface(
	INTERFACE	*intfc)
{
	G_USER_INTERFACE	*guh = g_user_hook(intfc->dim);
	int			i;
	int			num_params = guh->num_params;
	Gas_param		**params_list = guh->params_list;
	INTERFACE		*sav_intfc = current_interface();

	num_gas_params(intfc) = num_params;
	if (num_params == 0)
	{
	    gas_params_list(intfc) = NULL;
	    return;
	}

	set_current_interface(intfc);
	gas_params_list(intfc) =
	    (Gas_param**)Store(num_params*sizeof(Gas_param*));
	for (i = 0; i < num_params; ++i)
	    gas_params_list(intfc)[i] = params_list[i];
	set_current_interface(sav_intfc);

}	/*end set_params_list_for_interface*/


void	check_curve_states(INTERFACE *);

void	check_curve_states(
	INTERFACE	*intfc)
{
	CURVE		   **c;
	BOND		   *b;
	BOND_TRI	   **btris;
	Locstate	   sl, sr;

	/*printf("recv ed  %d \n", pp_mynode()); */

	for(c = intfc->curves; c && *c; c++)
	{
	    for(b=(*c)->first; b; b=b->next)
	    {
		printf("bond %llu  |  ", bond_number(b,intfc));
	        for(btris=Btris(b); btris && *btris; btris++)
		{
		    sl = left_end_btri_state(*btris);
		    sr = right_end_btri_state(*btris);
		    printf("%d  %d | ", sl, sr);
		    /*printf("%3d %3d  | ", gas_param_number(Params(sl)), gas_param_number(Params(sr))); */
		    /*printf("%d %d  | ", Params(sl), Params(sr)); */
		}
		printf("\n");
		
		printf("bond %llu  |  ", bond_number(b,intfc));
		for(btris=Btris(b); btris && *btris; btris++)
		{
		    sl = left_end_btri_state(*btris);
		    sr = right_end_btri_state(*btris);
		    printf("%d  %d | ", Params(sl), Params(sr));
		}
		printf("\n");
	    }
	}
}

LOCAL	void g_reconstruct_interface_pointers(
	INTERFACE	*nintfc,
	struct Table	*otbl,
	POINTER		*ocad,
	POINTER		*ncad)
{
	f_reconstruct_interface_pointers(nintfc,otbl,ocad,ncad);

}		/*end g_reconstruct_interface_pointers*/

LOCAL	void g_reconstruct_point_pointers(
	POINT		*p,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	f_reconstruct_point_pointers(p,nintfc,ointfc,ocad,ncad,nchks);
	/* see comments on  */
	/* i_send_interface */
	if ((nintfc->dim == 3) && (sorted(p) == YES)) /*Already done*/
	    return;
	if (size_of_state(nintfc) != 0)
	{
	    reconstruct_gas_params(left_state(p),nintfc,ointfc,
				   ocad,ncad,nchks);
	    reconstruct_gas_params(right_state(p),nintfc,ointfc,
				   ocad,ncad,nchks);
	    if (nintfc->dim == 3)
		sorted(p) = YES;
	}
}		/*end g_reconstruct_point_pointers*/



/*ARGSUSED*/
/*get the address of the Gas params in the current processor. In the initialization step,  */
/* proc_0	    proc_1 */
/*old_list[0],    new_list[0] */
/*old_list[1],    new_list[0] */
/*.		  . */
/*.		  . */
/*old_list[num_params-1],    new_list[num_params-1] */
/*represent same Gas_param */
/*This function uses this properties to convert the address of Gas_param. */

LOCAL	void	reconstruct_gas_params(
	Locstate	state,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	Gas_param	**new_list, **old_list;
	int		i, num_params;

	debug_print("reconstruct","Entered reconstruct_gas_params()\n");

	if (state == NULL || Params(state) == NULL)
	{
	    debug_print("reconstruct","Left reconstruct_gas_params()\n");
	    return;
	}
	num_params = return_params_list(&new_list);
	old_list = (Gas_param **) new_address(nintfc,gas_params_list(ointfc),
					      ocad,ncad,nchks);

	for (i = 0; i < num_params; ++i)
	{
	    if (Params(state) == old_list[i])
	    {
	        Params(state) = new_list[i];
	        break;
	    }
	}
	debug_print("reconstruct","Left reconstruct_gas_params()\n");
}		/*end reconstruct_gas_params*/

LOCAL	CURVE* g_make_curve(
	COMPONENT	left_c,
	COMPONENT	right_c,
	NODE		*start,
	NODE		*end)
{
	CURVE		*curve = f_make_curve(left_c,right_c,start,end);

	if (curve == NULL)
	    return curve;

	if (curve->interface->dim == 2)
	{
	    INTERFACE *intfc = curve->interface;
	    zero_scalar(left_start_state(curve),size_of_state(intfc));
	    zero_scalar(right_start_state(curve),size_of_state(intfc));
	    zero_scalar(left_end_state(curve),size_of_state(intfc));
	    zero_scalar(right_end_state(curve),size_of_state(intfc));
	    start_status(curve) = ERROR;
	    end_status(curve) = ERROR;
	    layer_index(Hyper_surf(curve)) = 0;
	    surface_tension(curve) = 0.0;
	}
	return curve;
}		/*end g_make_curve*/

LOCAL	CURVE	*g_copy_curve(
	CURVE		*curve,
	NODE		*start,
	NODE		*end)
{
	CURVE		*new_curve = f_copy_curve(curve, start, end);
	int i;

	if (new_curve == NULL)
	    return new_curve;
	if (curve->interface->dim == 2)
	{
	    start_status(new_curve) = start_status(curve);
	    end_status(new_curve) = end_status(curve);
	    surface_tension(new_curve) = surface_tension(curve);
	    layer_index(Hyper_surf(new_curve)) = layer_index(Hyper_surf(curve));
	    mom_inertial(new_curve) = mom_inertial(curve);
	    total_mass(new_curve) = total_mass(curve);
	    angular_velo(new_curve) = angular_velo(curve);
	    motion_type(new_curve) = motion_type(curve);
	    for (i = 0; i < curve->interface->dim; ++i)
	    {
	    	center_of_mass(new_curve)[i] = center_of_mass(curve)[i];
	    	center_of_mass_velo(new_curve)[i] = 
				center_of_mass_velo(curve)[i];
	    }
	}
	return new_curve;
}		/*end g_copy_curve*/


/* ARGSUSED */
LOCAL	void g_user_fprint_curve(
	FILE		*file,
	CURVE		*curve)
{
	f_user_fprint_curve(file,curve);
	if (curve->interface->dim == 2)
	{
	    fprint_curve_status(file,"\tcurve->start_status = ",
				start_status(curve));
	    fprint_curve_status(file,"\tcurve->end_status = ",
			        end_status(curve));
	    (void) fprintf(file,"\tcurve->surface_tension = ");
	    if (is_binary_output() == YES)
	    {
	    	(void) fprintf(file,"\f%c",1);
	    	(void) fwrite((const void *) &surface_tension(curve),
			      FLOAT,1,file);
		(void) fprintf(file,"\n");
	    }
	    else
	    	(void) fprintf(file,"%"FFMT"\n",surface_tension(curve));
	    (void) fprintf(file,"\tcurve->layer_index = %d\n",
	    		   layer_index(Hyper_surf(curve)));
	    if (wave_type(curve) == NEUMANN_BOUNDARY)
	    {
		(void) fprintf(file,"\tno slip = %s\n",
			       y_or_n(no_slip(Hyper_surf(curve))));
	    	if (no_slip(Hyper_surf(curve)) == YES)
		{
		    (void) fprintf(file,"\tadherence coefficient = ");
		    if (is_binary_output() == YES)
		    {
		    	(void) fprintf(file,"\f%c",1);
			(void) fwrite((const void *) &adherence_coeff(
				Hyper_surf(curve)),FLOAT,1,file);
		    	(void) fprintf(file,"\n");
		    }
		    else
		    	(void) fprintf(file,"%"FFMT"\n",adherence_coeff(
				Hyper_surf(curve)));
		}
	    }
	    (void) fprintf(file,"\n");
	}
}		/*end g_user_fprint_curve*/

LOCAL	boolean	g_user_split_curve(
	int		is_a_node,
	POINT		*p,
	BOND		*bond,
	CURVE		*curve,
	CURVE		**curves)
{
	if (f_user_split_curve(is_a_node,p,bond,curve,curves) != YES)
	    return NO;

	if (curve->interface->dim == 2)
	{
	    HYPER_SURF	*hs = Hyper_surf(curve);
	    INTERFACE	*intfc = curve->interface;
	    start_status(curves[0]) = start_status(curve);
	    end_status(curves[1]) = end_status(curve);
	    /* Default status at new node */
	    if (wave_type(curve) < FIRST_PHYSICS_WAVE_TYPE)
	    {
	    	end_status(curves[0])   = FIXED;
	    	start_status(curves[1]) = FIXED;
	    }
	    else
	    {
	    	end_status(curves[0])   = INCIDENT;
	    	start_status(curves[1]) = INCIDENT;
	    }
	    surface_tension(curves[0]) = surface_tension(curve);
	    surface_tension(curves[1]) = surface_tension(curve);
	    if (is_subdomain_boundary(hs))
	    {
		size_t sizest = size_of_state(intfc);
	    	obstacle_state(intfc,left_end_state(curves[0]),sizest);
	    	obstacle_state(intfc,right_end_state(curves[0]),sizest);
	    	obstacle_state(intfc,left_start_state(curves[1]),sizest);
	    	obstacle_state(intfc,right_start_state(curves[1]),sizest);
	    }
	    layer_index(Hyper_surf(curves[0])) = layer_index(hs);
	    layer_index(Hyper_surf(curves[1])) = layer_index(hs);
	    return YES;
	}
	return YES;
}		/*end g_user_split_curve*/

LOCAL	boolean g_user_join_curves(
	CURVE		*curve,
	CURVE		*curve1,
	CURVE		*curve2)
{
	if (f_user_join_curves(curve,curve1,curve2) != YES)
	    return NO;
	if (curve->interface->dim == 2)
	{
	    HYPER_SURF *hs = Hyper_surf(curve);
	    HYPER_SURF *hs1 = Hyper_surf(curve1);
	    HYPER_SURF *hs2 = Hyper_surf(curve2);
	    start_status(curve) = start_status(curve1);
	    end_status(curve) = end_status(curve2);
	    surface_tension(curve) = 0.5*(surface_tension(curve1) +
				          surface_tension(curve2));
	    if (layer_index(hs1) == layer_index(hs2))
	        layer_index(hs) = layer_index(hs1);
	}
	return YES;
}		/*end g_user_join_curves*/

LOCAL	NODE *g_make_node(
	POINT		*posn)
{
	NODE		*newnod = f_make_node(posn);

	if (newnod->interface->dim == 2)
	{
	    static ADJUST_ANGLE_DATA dflt_adjust_ang = {
							 0.0,
							 g_adjust_angle_len
						       };

	    size_t sizest = size_of_state(newnod->interface);
	    Rp_data(newnod) = allocate_RP_DATA_structure(sizest,YES,GAS_STATE);
	    adjust_angle(newnod) = dflt_adjust_ang;
	}

	return newnod;
}		/*end g_make_node*/

LOCAL	NODE *g_copy_node(
	NODE		*old_node)
{
	NODE		*new_node = f_copy_node(old_node);

	if (new_node == NULL)
	    return new_node;

	if ((new_node->interface->dim == 2) &&
	    (interface_type(new_node->interface) == PHYSICAL_INTERFACE) &&
	    (copy_intfc_states() == YES))
		copy_RP_DATA_structure(Rp_data(new_node),Rp_data(old_node));
	adjust_angle(new_node) = adjust_angle(old_node);
	adjust_len(new_node) = 0.0;
	return new_node;
}		/*end g_copy_node*/


EXPORT	void g_user_fprint_node(
	FILE		*file,
	NODE		*node)
{
	static char	vname[3][3] = {"vx","vy","vz"};
	int		i, dim = node->interface->dim;

	f_user_fprint_node(file,node);
	if (node->interface &&
		interface_type(node->interface) == PHYSICAL_INTERFACE)
	{
	    (void) fprintf(file,"\t\tNode Velocity  ");
	    if (is_binary_output() == YES)
	    {
	    	(void) fprintf(file,"\f%c",dim);
	    	(void) fwrite((const void *) Node_vel(node),FLOAT,dim,file);
		(void) fprintf(file,"\n");
	    }
	    else
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    (void) fprintf(file,"%s = %"FFMT"%s",
				   vname[i],Node_vel(node)[i],
	    			   (i==(dim-1))?"\n":"   ");       
	    	}
	    }
	}       
}		/*end g_user_fprint_node*/

LOCAL	void g_reconstruct_bond_pointers(
	BOND		*b,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	BOND_TRI		**btris;

	/*printf("revbond %llu \n", bond_number(b,nintfc)); */
	f_reconstruct_bond_pointers(b,nintfc,ointfc,ocad,ncad,nchks);

	if (size_of_state(nintfc) == 0)
	    return;
	
	/*for (btris = Btris(b); btris && *btris; ++btris) */
	/*{ */
	/*	sl = left_start_btri_state(*btris); */
	/*	sr = right_start_btri_state(*btris); */
	/*	printf("%d  %d | ", Params(sl), Params(sr)); */
	/*} */
	/*printf("\n"); */

	/*see i_send_interface for consistence. */
	for (btris = Btris(b); btris && *btris; ++btris)
	{
	    reconstruct_gas_params(left_start_btri_state(*btris),
			       nintfc,ointfc,ocad,ncad,nchks);
	    reconstruct_gas_params(right_start_btri_state(*btris),
			       nintfc,ointfc,ocad,ncad,nchks);
	
	    /*see i_send_interface */
	    /*    f_reorder_curve_link_list will make the link consistent. */
	    /*if the last bond and not closed curve, the end state should be reconstructed. */
	    if(b->next == NULL && !is_closed_curve((*btris)->curve))
	    {
	        reconstruct_gas_params(left_end_btri_state(*btris),
			       nintfc,ointfc,ocad,ncad,nchks);
	        reconstruct_gas_params(right_end_btri_state(*btris),
			       nintfc,ointfc,ocad,ncad,nchks);
	    }
	}

}		/*end g_reconstruct_bond_pointers*/

LOCAL	void g_reconstruct_node_pointers(
	NODE		*n,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	debug_print("reconstruct","Entered g_user_reconstruct_node_pointers()\n");

	f_reconstruct_node_pointers(n,nintfc,ointfc,ocad,ncad,nchks);
	if (Rp_data(n) != NULL)
	{
		int i;
		Rp_data(n) = (RP_DATA *) new_address(nintfc,Rp_data(n),ocad,
						     ncad,nchks);

		for (i = 0; i < MAX_N_CURVES; ++i)
		{
			Rp_data(n)->state[i] =
			    (Locstate) new_address(nintfc,Rp_data(n)->state[i],
						   ocad,ncad,nchks);
			reconstruct_gas_params(Rp_data(n)->state[i],
					       nintfc,ointfc,ocad,ncad,nchks);
		}
	}
	debug_print("reconstruct","Left g_user_reconstruct_node_pointers()\n");
}		/*end g_reconstruct_node_pointers*/

LOCAL	void g_reconstruct_curve_pointers(
	CURVE		*c,
	INTERFACE	*nintfc,
	INTERFACE	*ointfc,
	POINTER		*ocad,
	POINTER		*ncad,
	int		nchks)
{
	debug_print("reconstruct","Entered g_reconstruct_curve_pointers()\n");

	f_reconstruct_curve_pointers(c,nintfc,ointfc,ocad,ncad,nchks);
	if ((nintfc->dim == 2) && (size_of_state(nintfc) != 0))
	{
		reconstruct_gas_params(left_start_state(c),
				       nintfc,ointfc,ocad,ncad,nchks);
		reconstruct_gas_params(right_start_state(c),
				       nintfc,ointfc,ocad,ncad,nchks);
		reconstruct_gas_params(left_end_state(c),
				       nintfc,ointfc,ocad,ncad,nchks);
		reconstruct_gas_params(right_end_state(c),
				       nintfc,ointfc,ocad,ncad,nchks);
	}
	debug_print("reconstruct","Left g_reconstruct_curve_pointers()\n");
}		/*end g_reconstruct_curve_pointers*/

LOCAL	void	g_invert_curve(
	CURVE		*c)
{
	f_invert_curve(c);
	if (c->interface->dim == 2)
	{
	    wave_type(c) = opposite_wave_type(wave_type(c));
	}
}		/*end g_invert_curve*/

LOCAL	void	g_reverse_curve(
	CURVE		*c)
{
	f_reverse_curve(c);
	if (c->interface->dim == 2)
	{
	    int status;
	    status = start_status(c);
	    start_status(c) = end_status(c);
	    end_status(c) = status;
	}
}		/*end g_reverse_curve*/


LOCAL	CURVE *g_attach_curve_to_node(
	CURVE		*c1,
	POINT		*p,
	BOND		*b,
	NODE		*n)
{
	CURVE		*c2;

	c2 = f_attach_curve_to_node(c1,p,b,n);
	node_type(n) = ERROR;
	return c2;
}		/*end g_attach_curve_to_node*/

LOCAL	void g_user_fprint_point(
	FILE		*file,
	POINT		*point)
{
	if (point->interface != NULL)
	{
	    fprint_wave_type(file,"\t\twave_type = ",wave_type(point),"\n",
	    	                  Hyper_surf(point)->interface);
	    (void) fprintf(file,"\tlayer_index = %d\n",
				layer_index(Hyper_surf(point)));
	}
	(void) fprintf(file,"\t\tLeft state:\n");
	g_fprint_state(file,left_state(point));
	(void) fprintf(file,"\t\tRight state:\n");
	g_fprint_state(file,right_state(point));
}		/*end g_fprint_point*/


/* ASSUME the moving curve is NEUMANN_CURVE */
LOCAL   SURFACE *g_detach_one_surface(
	SURFACE *s)
{
	INTERFACE	*intfc = current_interface();
	SURFACE 	*surf; 
	CURVE		**c;
	boolean		found = NO;

	/*curves on the wall */
	for(c=s->pos_curves; c && *c; c++)
	{
	    if(curve_type(*c) == NEUMANN_CURVE)
	    {
	        curve_type(*c) = NEUMANN_CURVE_W;
	        found = YES;
	    }
	}
	for(c=s->neg_curves; c && *c; c++)
	{
	    if(curve_type(*c) == NEUMANN_CURVE)
	    {
	        curve_type(*c) = NEUMANN_CURVE_W;
	        found = YES;
	    }
	}

	/*printf("#found detach = %d  \n", found); */
	if(!found)
	    return s;

	surf = f_detach_one_surface(s);

	/*curves on the physical surface */
	for(c=surf->pos_curves; c && *c; c++)
	{
	    if(curve_type(*c) == NEUMANN_CURVE_W)
	        curve_type(*c) = NEUMANN_CURVE_P;
	}
	for(c=surf->neg_curves; c && *c; c++)
	{
	    if(curve_type(*c) == NEUMANN_CURVE_W)
	        curve_type(*c) = NEUMANN_CURVE_P;
	}
	
	return surf;
}

LOCAL	SURFACE	*g_copy_surface(
	SURFACE	*surf,
	CURVE	**pos,
	CURVE	**neg,
	boolean copy_tris)
{
	SURFACE	*new_surf = f_copy_surface(surf,pos,neg,copy_tris);

	if (new_surf == NULL)
	    return new_surf;

	surface_tension(new_surf) = surface_tension(surf);
	layer_index(Hyper_surf(new_surf)) = layer_index(Hyper_surf(surf));
	return new_surf;
}		/*end g_copy_surface*/

LOCAL	SURFACE *g_join_surfaces(
	CURVE *c)
{
	SURFACE *news, *sn, *sp;
	double   stn, stp;
	int     lin, lip;

	if (!find_surfaces_to_join_at_curve(c,&sn,&sp))
	    return NULL;
	stn = surface_tension(sn);
	stp = surface_tension(sp);
	lin = layer_index(Hyper_surf(sn));
	lip = layer_index(Hyper_surf(sp));
	news = f_join_surfaces(c);
	if (sn != sp)
	{
	    surface_tension(news) = 0.5*(stn+stp);
	    if (lin == lip)
	        layer_index(Hyper_surf(news)) = lin;
	}
	return news;
}		/*end g_join_surfaces*/

LOCAL	void	g_user_fprint_surface(
	FILE		*file,
	SURFACE		*s)
{
	HYPER_SURF         *hs;

	f_user_fprint_surface(file,s);

	hs = Hyper_surf(s);
	(void) fprintf(file,"\tsurface_tension(hs) = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &surface_tension(hs),
	                   FLOAT,1,file);
	    (void) fprintf(file,"\n");
	}
	else
	    (void) fprintf(file,"%"FFMT"\n",surface_tension(hs));

	if (wave_type(s) == NEUMANN_BOUNDARY)
	{
	    (void) fprintf(file,"\tno slip = %s\n",
			   y_or_n(no_slip(Hyper_surf(s))));
	    if (no_slip(Hyper_surf(s)) == YES)
	    {
		(void) fprintf(file,"\tadherence coefficient = ");
		if (is_binary_output() == YES)
		{
		    (void) fprintf(file,"\f%c",1);
		    (void) fwrite((const void *) &adherence_coeff(
				Hyper_surf(s)),FLOAT,1,file);
		    (void) fprintf(file,"\n");
		}
		else
		    (void) fprintf(file,"%"FFMT"\n",adherence_coeff(
				Hyper_surf(s)));
	    }
	}
}		/*end g_user_fprint_surface*/


EXPORT	boolean	g_form_subintfc_via_communication2d(
	Front		*fr)
{
	CURVE	  **c;
	INTERFACE *intfc;
	RECT_GRID *gr = fr->rect_grid;
	boolean	  status;
	boolean      sav_intrp;
	int       i, d;

	debug_print("fscatter","Entered g_form_subintfc_via_communication2d()\n");
	status = f_intfc_communication2d(fr);
	intfc = fr->interf;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
	    {
	    	if (start_status(*c) == ERROR)
	    	    start_status(*c) = INCIDENT;
	    	if (end_status(*c) == ERROR)
	    	    end_status(*c) = INCIDENT;;
	    }
	    else
	    {
	    	switch (wave_type(*c))
	    	{
	    	case SUBDOMAIN_BOUNDARY:
	    	    if (start_status(*c) == ERROR)
	    	    	start_status(*c) = VIRTUAL;
	    	    if (end_status(*c) == ERROR)
	    	    	end_status(*c) = VIRTUAL;;
	    	    break;
	    	case PASSIVE_BOUNDARY:
	    	    if (start_status(*c) == ERROR)
	    	    	start_status(*c) = PASSIVE;
	    	    if (end_status(*c) == ERROR)
	    	    	end_status(*c) = PASSIVE;;
	    	    break;
	    	default:
	    	    if (start_status(*c) == ERROR)
	    	    	start_status(*c) = FIXED;
	    	    if (end_status(*c) == ERROR)
	    	    	end_status(*c) = FIXED;;
	    	    break;
	    	}
	    }
	}

	/*
	* The following code is designed to enforce reflection symmetry
	* along a curve crossing a reflection boundary.  First points
	* within a tolerance of the boundary are moved to the boundary. Next
	* if necessary new points exactly at the boundary are inserted.
	* Finally the states at the new points are modified so that the boundary
	* normal component of velocity vanishes at these points.
	*/

	sav_intrp = interpolate_intfc_states(intfc);
	interpolate_intfc_states(intfc) = YES;
	for (i = 0; i < 2; ++i)
	{
	    for (d = 0; d < 2; ++d)
	    {
	        if ((rect_boundary_type(intfc,i,d) == REFLECTION_BOUNDARY))
		{
	            BOND  *b;
		    double tol = MIN_SC_SEP(intfc);/*TOLERANCE*/
		    double *h = gr->h;
	            double crds[3];
		    double s, e, R;
		    int   dim = gr->dim;
    
		    R = (d == 0) ? gr->L[i] : gr->U[i];
	            for (c = intfc->curves; c && *c; ++c)
		    {
	                if (omit_redistribution(*c) ||
		            (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE))
		            continue;
			for (b = (*c)->first; b != NULL; b = b->next)
			{
			    s = Coords(b->start)[i] - R;
			    e = Coords(b->end)[i] - R;
		            if (s*e <= 0.0)
			    {
			        coords_on_axis(i,R,b,crds);
				s = _scaled_separation(Coords(b->start),
				                       crds,h,dim);
				e = _scaled_separation(Coords(b->end),
				                       crds,h,dim);
		                while ((s < tol) && (e < tol))
				{
				  if (b->prev)
				  {
				    delete_end_of_bond(b->prev,*c);
				    s = _scaled_separation(Coords(b->start),
				                           crds,h,dim);
				  }
				  else if (b->next)
				  {
				    delete_start_of_bond(b->next,*c);
				    e = _scaled_separation(Coords(b->end),crds,
				                           h,dim);
				  }
				  else
				  {
				    delete_curve(*c);
				    c = intfc->curves - 1;
				    b = NULL;
				    break;
				  }
				}
				if (b == NULL)
				  break;
		                if (s < tol)
		                {
			          Coords(b->start)[0] = crds[0];
			          Coords(b->start)[1] = crds[1];
		                  states_on_axis(i,d,b->start,b,*c,fr);
				  if ((b->prev) &&
				    (_scaled_separation(Coords(b->prev->start),
				                        crds,h,dim) < tol))
				  {
				      if (b->prev->prev)
				          delete_start_of_bond(b->prev,*c);
				      else
				      {
			                Coords(b->prev->start)[0] = crds[0];
			                Coords(b->prev->start)[1] = crds[1];
		                        states_on_axis(i,d,b->prev->start,
				                       b->prev,*c,fr);
					delete_end_of_bond(b->prev,*c);
				      }
				  }
		                }
		                else if (e < tol)
		                {
			          Coords(b->end)[0] = crds[0];
			          Coords(b->end)[1] = crds[1];
			          states_on_axis(i,d,b->end,b,*c,fr);
				  if (b->next)
				  {
			            b = b->next;
				    if (_scaled_separation(Coords(b->end),crds,
					                  h,dim) < tol)
				    {
				      if (b->next)
				        delete_start_of_bond(b->next,*c);
				      else
				      {
			                Coords(b->end)[0] = crds[0];
			                Coords(b->end)[1] = crds[1];
			                states_on_axis(i,d,b->end,b,*c,fr);
					delete_start_of_bond(b,*c);
					break;
	                              }
				    }
				  }
				  else
				    break;
		                }
	                        else
		                {
		                    POINT *p;
		                    p = Point(crds);
			            if (insert_point_in_bond(p,b,*c))
			            {
			                states_on_axis(i,d,p,b,*c,fr);
			                b = b->next;
			            }
			            else
			                status = NO;
		                }
			    }
			}
		    }
		}
	    }
	}
	interpolate_intfc_states(intfc) = sav_intrp;
	status = pp_min_status(status);

	if (debugging("fscatter"))
	{
	    (void) printf("Final interface at end of ");
	    (void) printf("g_form_subintfc_via_communication2d()\n");
	    print_interface(fr->interf);
	}
	debug_print("fscatter","Left g_form_subintfc_via_communication2d()\n");
	return status;
}		/*end g_form_subintfc_via_communication2d*/


#if defined(USE_OVERTURE)
EXPORT  void    g_assembly_fine_patch_fronts_to_one(
        Front    **oldfrs,
        Front    *newfr)
{
        CURVE     **c;
        INTERFACE *intfc;
        RECT_GRID *gr;
        boolean      status;
        boolean      sav_intrp;
        int       i, d;

        debug_print("fscatter","Entered g_assembly_fine_patch_fronts_to_one()\n");
        /* 051503, test new assembly patch interface alg.     */
        /* The physical boundary curves are reassembled also. */
        /*
        assembly_fine_patch_fronts_to_one2d(oldfrs,newfr);
        */
        assembly_fine_patch_fronts_to_one2d_ver2(oldfrs,newfr);
        gr = newfr->rect_grid;
        intfc = newfr->interf;
        for (c = intfc->curves; c && *c; ++c)
        {
            if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
            {
                if (start_status(*c) == ERROR)
                    start_status(*c) = INCIDENT;
                if (end_status(*c) == ERROR)
                    end_status(*c) = INCIDENT;;
            }
            else
            {
                switch (wave_type(*c))
                {
                case SUBDOMAIN_BOUNDARY:
                /* 052703 closed
                case AMR_SUBDOMAIN_BOUNDARY:
                */
                    if (start_status(*c) == ERROR)
                        start_status(*c) = VIRTUAL;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = VIRTUAL;;
                    break;
                case PASSIVE_BOUNDARY:
                    if (start_status(*c) == ERROR)
                        start_status(*c) = PASSIVE;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = PASSIVE;;
                    break;
                case AMR_SUBDOMAIN_BOUNDARY:
                    printf("ERROR: g_assembly_fine_patch_fronts_to_one\n");
                    printf("CURVE wave type can not be AMR_SUBDOMAIN_BOUNDARY\n");
                    print_curve(*c);
                    print_interface(intfc);
                    clean_up(ERROR);
                    break;
                default:
                    if (start_status(*c) == ERROR)
                        start_status(*c) = FIXED;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = FIXED;;
                    break;
                }
            }
        }
        debug_print("fscatter","Left g_assembly_fine_patch_fronts_to_one()\n");
}


EXPORT  boolean    g_form_patch_subintfc_via_cut2d(
        Front           *fr)
{
        CURVE     **c;
        INTERFACE *intfc;
        RECT_GRID *gr = fr->rect_grid;
        boolean      status;
        boolean      sav_intrp;
        int       i, d;

        debug_print("fscatter","Entered g_form_patch_subintfc_via_cut2d()\n");
        status = f_form_patch_subintfc_via_cut2d(fr);
        intfc = fr->interf;
        for (c = intfc->curves; c && *c; ++c)
        {
            if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
            {
                if (start_status(*c) == ERROR)
                    start_status(*c) = INCIDENT;
                if (end_status(*c) == ERROR)
                    end_status(*c) = INCIDENT;;
            }
            else
            {
                switch (wave_type(*c))
                {
                case SUBDOMAIN_BOUNDARY:
                /* closed on 052703
                case AMR_SUBDOMAIN_BOUNDARY:
                */
                    if (start_status(*c) == ERROR)
                        start_status(*c) = VIRTUAL;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = VIRTUAL;;
                    break;
                case PASSIVE_BOUNDARY:
                    if (start_status(*c) == ERROR)
                        start_status(*c) = PASSIVE;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = PASSIVE;;
                    break;
                case AMR_SUBDOMAIN_BOUNDARY:
                    printf("ERROR: g_form_patch_subintfc_via_cut2d\n");
                    printf("CURVE wave type can not be AMR_SUBDOMAIN_BOUNDARY\n");
                    print_curve(*c);
                    print_interface(intfc);
                    clean_up(ERROR);
                    break;
                default:
                    if (start_status(*c) == ERROR)
                        start_status(*c) = FIXED;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = FIXED;;
                    break;
                }
            }
        }


        /*
        * The following code is designed to enforce reflection symmetry
        * along a curve crossing a reflection boundary.  First points
        * within a tolerance of the boundary are moved to the boundary. Next
        * if necessary new points exactly at the boundary are inserted.
        * Finally the states at the new points are modified so that the boundary
        * normal component of velocity vanishes at these points.
        */

        /* Do not need this piece of code, it's done on the base front */
        /*  042103 added */
        /*
        sav_intrp = interpolate_intfc_states(intfc);
        interpolate_intfc_states(intfc) = YES;
        for (i = 0; i < 2; ++i)
        {
            for (d = 0; d < 2; ++d)
            {
                if ((rect_boundary_type(intfc,i,d) == REFLECTION_BOUNDARY))
                {
                    BOND  *b;
                    double tol = MIN_SC_SEP(intfc);/*TOLERANCE */
                    double *h = gr->h;
                    double crds[3];
                    double s, e, R;
                    int   dim = gr->dim;

                    R = (d == 0) ? gr->L[i] : gr->U[i];
                    for (c = intfc->curves; c && *c; ++c)
                    {
                        if (omit_redistribution(*c) ||
                            (wave_type(*c) < FIRST_PHYSICS_WAVE_TYPE))
                            continue;
                        for (b = (*c)->first; b != NULL; b = b->next)
                        {
                            s = Coords(b->start)[i] - R;
                            e = Coords(b->end)[i] - R;
                            if (s*e <= 0.0)
                            {
                                coords_on_axis(i,R,b,crds);
                                s = _scaled_separation(Coords(b->start),
                                                       crds,h,dim);
                                e = _scaled_separation(Coords(b->end),
                                                       crds,h,dim);
                                while ((s < tol) && (e < tol))
                                {
                                  if (b->prev)
                                  {
                                    delete_end_of_bond(b->prev,*c);
                                    s = _scaled_separation(Coords(b->start),
                                                           crds,h,dim);
                                  }
                                  else if (b->next)
                                  {
                                    delete_start_of_bond(b->next,*c);
                                    e = _scaled_separation(Coords(b->end),crds,
                                                           h,dim);
                                  }
                                  else
                                  {
                                    delete_curve(*c);
                                    c = intfc->curves - 1;
                                    b = NULL;
                                    break;
                                  }
                                }
                                if (b == NULL)
                                  break;
                                if (s < tol)
                                {
                                  Coords(b->start)[0] = crds[0];
                                  Coords(b->start)[1] = crds[1];
                                  states_on_axis(i,d,b->start,b,*c,fr);
                                  if ((b->prev) &&
                                    (_scaled_separation(Coords(b->prev->start),
                                                        crds,h,dim) < tol))
                                  {
                                      if (b->prev->prev)
                                          delete_start_of_bond(b->prev,*c);
                                      else
                                      {
                                        Coords(b->prev->start)[0] = crds[0];
                                        Coords(b->prev->start)[1] = crds[1];
                                        states_on_axis(i,d,b->prev->start,
                                                       b->prev,*c,fr);
                                        delete_end_of_bond(b->prev,*c);
                                      }
                                  }
                                }
                                else if (e < tol)
                                {
                                  Coords(b->end)[0] = crds[0];
                                  Coords(b->end)[1] = crds[1];
                                  states_on_axis(i,d,b->end,b,*c,fr);
                                  if (b->next)
                                  {
                                    b = b->next;
                                    if (_scaled_separation(Coords(b->end),crds,
                                                          h,dim) < tol)
                                    {
                                      if (b->next)
                                        delete_start_of_bond(b->next,*c);
                                      else
                                      {
                                        Coords(b->end)[0] = crds[0];
                                        Coords(b->end)[1] = crds[1];
                                        states_on_axis(i,d,b->end,b,*c,fr);
                                        delete_start_of_bond(b,*c);
                                        break;
                                      }
                                    }
                                  }
                                  else
                                    break;
                                }
                                else
                                {
                                    POINT *p;
                                    p = Point(crds);
                                    if (insert_point_in_bond(p,b,*c))
                                    {
                                        states_on_axis(i,d,p,b,*c,fr);
                                        b = b->next;
                                    }
                                    else
                                        status = NO;
                                }
                            }
                        }
                    }
                }
            }
        }
        interpolate_intfc_states(intfc) = sav_intrp;
        status = pp_min_status(status);
        */
        if (debugging("fscatter"))
        {
            (void) printf("Final interface at end of ");
            (void) printf("g_form_patch_subintfc_via_cut2d()\n");
            print_interface(fr->interf);
        }
        debug_print("fscatter","Left g_form_patch_subintfc_via_cut2d()\n");
        return status;
}               /*end g_form_patch_subintfc_via_cut2d*/


EXPORT  boolean    g_form_patch_subintfc_2d(
        Front           *fr,
        COMPONENT       i_comp)
{
        CURVE     **c;
        INTERFACE *intfc;
        RECT_GRID *gr = fr->rect_grid;
        boolean      status;
        boolean      sav_intrp;
        int       i, d;

        debug_print("fscatter","Entered g_form_patch_subintfc_2d()\n");
        status = f_form_patch_subintfc_2d(fr,i_comp);
        intfc = fr->interf;
        for (c = intfc->curves; c && *c; ++c)
        {
            if (wave_type(*c) >= FIRST_PHYSICS_WAVE_TYPE)
            {
                if (start_status(*c) == ERROR)
                    start_status(*c) = INCIDENT;
                if (end_status(*c) == ERROR)
                    end_status(*c) = INCIDENT;;
            }
            else
            {
                switch (wave_type(*c))
                {
                case SUBDOMAIN_BOUNDARY:
                    if (start_status(*c) == ERROR)
                        start_status(*c) = VIRTUAL;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = VIRTUAL;;
                    break;
                case PASSIVE_BOUNDARY:
                    if (start_status(*c) == ERROR)
                        start_status(*c) = PASSIVE;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = PASSIVE;;
                    break;
                case AMR_SUBDOMAIN_BOUNDARY:
                    printf("ERROR: g_form_patch_subintfc_via_cut2d\n");
                    printf("CURVE wave type can not be AMR_SUBDOMAIN_BOUNDARY\n");
                    print_curve(*c);
                    print_interface(intfc);
                    clean_up(ERROR);
                    break;
                default:
                    if (start_status(*c) == ERROR)
                        start_status(*c) = FIXED;
                    if (end_status(*c) == ERROR)
                        end_status(*c) = FIXED;;
                    break;
                }
            }
        }

        if (debugging("fscatter"))
        {
            (void) printf("Final interface at end of ");
            (void) printf("g_form_patch_subintfc_2d()\n");
            print_interface(fr->interf);
        }
        debug_print("fscatter","Left g_form_patch_subintfc_2d()\n");
        return status;
}               /*end g_form_patch_subintfc_via_cut2d*/
#endif /* if defined(USE_OVERTURE) */



LOCAL	void	coords_on_axis(
	int   i,
	double R,
	BOND  *b,
	double *crds)
{
	int   j = (i+1)%2;
	double a;
	
	if (Coords(b->end)[i]-Coords(b->start)[i] == 0.0)
	{
	    crds[i] = R;
	    crds[j] = 0.5*(Coords(b->start)[j] + Coords(b->end)[j]);
	    crds[2] = 0.0;
	    return;
	}
	else
	{
	    a = (R - Coords(b->start)[i])/
	          (Coords(b->end)[i]-Coords(b->start)[i]);

	    crds[i] = R;
	    crds[j] = (1.0-a)*Coords(b->start)[j] + a*Coords(b->end)[j];
	    crds[2] = 0.0;
	    return;
	}
}

LOCAL	void states_on_axis(
	int   i,
	int   d,
	POINT *p,
	BOND  *b,
	CURVE *c,
	Front *fr)
{
	Locstate        s;
	double           dir[3];
	static Locstate sl = NULL, sr = NULL;

	dir[0] = dir[1] = dir[2] = 0.0;
	dir[i] = 1.0;
	if (sl == NULL)
	{
	    alloc_state(fr->interf,&sl,fr->sizest);
	    alloc_state(fr->interf,&sr,fr->sizest);
	}
	s = left_state_at_point_on_curve(p,b,c);

	copy_state(sr,s);
	copy_state(sl,s);
	if (d == 0)
	    Mom(sl)[i] *= -1.0;
	else
	    Mom(sr)[i] *= -1.0;
	riemann_solution(0.0,dir,sl,sr,s,state_type(s));

	s = right_state_at_point_on_curve(p,b,c);
	copy_state(sr,s);
	copy_state(sl,s);
	if (d == 0)
	    Mom(sl)[i] *= -1.0;
	else
	    Mom(sr)[i] *= -1.0;
	riemann_solution(0.0,dir,sl,sr,s,state_type(s));
}		/*end states_on_axis*/


#include <plotdecs.h>

LOCAL	void	g_fset_hyper_surf_color(
	FILE            *file,
	HYPER_SURF	*hs)
{
	int		w_type;

	w_type = wave_type(hs);
	if (w_type >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
	{
	    if (is_shock_wave(w_type))
	    	fset_color_from_table(file,5);
	    else if (is_rarefaction_leading_edge(w_type))
	    	fset_color_from_table(file,3);
	    else if (is_rarefaction_trailing_edge(w_type))
	    	fset_color_from_table(file,4);
	}
	else
	    f_fset_hyper_surf_color(file,hs);
}		/*end g_fset_hyper_surf_color*/


void    tecplot_interface_states(const char*, INTERFACE	*);

LOCAL	void  g_check_print_intfc(
	const char  *msg,
	const char  *fname,
	char   	    ch,
	INTERFACE   *intfc, 
	int	    step,
	int         fstep,
	boolean	    final)
{
	char    s[50], sn[50];
	boolean 	fcon;

	if(intfc->dim != 3 || !debugging("check_print_intfc"))
	    return;

	/*check restart debug code. */
	if(NO)
	{
	    SURFACE	**s;
	    POINT	*p;
	    TRI		*tri;
	    int		i;

	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		print_wave_type("wtype = ", wave_type(*s), "\n", intfc);

		tri=first_tri(*s);
		for(i=0; i<3; i++)
		{
		    p = Point_of_tri(tri)[i];
		    print_general_vector("nor", normal_at_point(p), 3, "\n");
		}
		printf("\n");
	    }
	}

	if( fname != NULL && step == fstep )
	{
	    sprintf(sn, "%s", right_flush(pp_mynode(),PP_NODE_FIELD_WIDTH));
	    sprintf(s, "%s%s_%s.plt", fname, 
	    		right_flush(step, TSTEP_FIELD_WIDTH), sn);
	    (void) printf("#check_print_intfc %d %s\n", step, s);
	    if(ch != 's')
	        tecplot_interface(s, intfc);
	    else
	        tecplot_interface_states(s, intfc);
	}

	if(ch == 's')
	    ch = 'g';

	switch(ch)
	{
	case 'i':
	    fcon = i_consistent_interface(intfc);
	    break;
	case 'f':
	    fcon = f_consistent_interface(intfc);
	    break;
	default: 
	    fcon = consistent_interface(intfc);
	}

	if(fcon)
	{
	    (void) printf("%s %d, interface is consistent.\n", msg, step);
	}
	else
	{
	    (void) printf("%s %d, interface is inconsistent.\n", msg, step);
	    clean_up(0);
	}


	if(step == fstep && final)
	{
	    long_alloc_view(stdout);
	    clean_up(0);
	}
}

/*debug functions for  reconstruct_intfc3d_in_box */

/*wall_fix_crx */
/*init: set_wall_crossings3d */
/*use:  fill_block_crx */
LOCAL   void  g_print_wall_crx(
	const char *msg,
	int  	*ip, 
	int 	dir, 
	int 	k, 
	CRXING  *crx)
{
	POINT	  *p = crx->pt;
	Locstate  sl = left_state(p), sr = right_state(p);

	printf("#%s  %5d %5d %5d  %2d | %9d |", msg, ip[0], ip[1], ip[2], dir, k);
	printf(" %d | %2d %2d | %d | %d  %d\n", crx, crx->lcomp, crx->ucomp, 
	        Surface_of_hs(crx->hs), Params(sl), Params(sr));
	print_general_vector("crds= ", Coords(p), 3, "\n");
}

/*init:  insert_curve_face_crossings */
/*use:   fill_block_curve_crx */
LOCAL   void  g_print_wall_curve_crx(
	const char *msg,
	int  	*ip, 
	int 	dir, 
	int 	k, 
	CRXING  *crx)
{
	POINT	  *p = crx->pt;
	Locstate  sl = left_state(p), sr = right_state(p);

	/*print the crx-point pair for check */
	printf("#%s %5d %5d %5d  %2d | %9d |", msg, ip[0], ip[1], ip[2], dir, k);
	printf("[%d %d]| %2d %2d | %d | %d  %d\n", crx, p, crx->lcomp, crx->ucomp, 
	        crx->hsb != NULL ? Curve_of_hsb(crx->hsb) : NULL, 
		Params(sl), Params(sr));
	print_general_vector("crds= ", Coords(p), 3, "\n\n");
	/*verbose_print_state("sl", sl); */
	/*verbose_print_state("sr", sr); */
}

/*init:  insert_curve_face_crossings */
/*use:   install_curve_points_state */
/*check the crx p pair */
LOCAL   void  g_print_wall_curve_crx0(
	const char *msg,
	POINT	*p0,
	int 	k, 
	CRXING  *crx)
{
	POINT	  *p = crx->pt;
	Locstate  sl = left_state(p), sr = right_state(p);

	/*print the crx-point pair for check */
	printf("#%s | %9d |", msg, k);
	printf("[%d %d]| %2d %2d | %d | %d  %d\n", crx, p, crx->lcomp, crx->ucomp, 
	        crx->hsb != NULL ? Curve_of_hsb(crx->hsb) : NULL, 
	        Params(sl), Params(sr));
	/*print_general_vector("crds= ", Coords(p), 3, "\n"); */
	if(p0 != p)
	{
	    printf("ERROR print_wall_curve_crx0, curve crx point inconsistent\n");
	    printf("p0= %d  p= %d \n", p0, p);
	    clean_up(ERROR);
	}
}


