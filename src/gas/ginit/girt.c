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
*				girt.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <ginit/ginit.h>

 
	/* LOCAL Function Declarations */
#if defined(TWOD) || defined(THREED)
LOCAL	void	get_state_rt_kh(double*,Locstate,COMP_TYPE*,
				HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	kh_add_to_state(double*,Locstate,COMP_TYPE*);
LOCAL	void	rt_add_to_state(double*,Locstate,COMP_TYPE*,INIT_DATA*);
#endif /* defined(TWOD) || defined(THREED) */


EXPORT	STRATIFICATION_TYPE prompt_for_stratification(
	const char *name)
{
	char s[Gets_BUF_SIZE];
	if (is_gravity() == YES)
	{
	    screen("\nInput the gravity stratified property %s\n",name);
	    screen("\tConstant state (CO) \n");
	    screen("\tConstant density (CD) \n");
	    screen("\tHydro static (HS) \n");
	    screen("\tAdiabatic (AD) \n");
	    screen("\tIsothermal (IS, default) \n");
	    screen("Enter choice: ");
	    (void) Gets(s);
	    if ((strncmp(s,"CO",2) == 0) || (strncmp(s,"co",2) == 0))
	    {
	    	set_stratified_state(NULL);
	    	return CONSTANT;
	    }
	    else if ((strncmp(s,"CD",2) == 0) || (strncmp(s,"cd",2) == 0))
	    {
	       set_stratified_state(constant_density_stratified_state);
	       return CONSTANT_DENSITY;
	    }
	    else if ((strncmp(s,"AD",2) == 0) || (strncmp(s,"ad",2) == 0))
	    {
	    	set_stratified_state(isentropic_stratified_state);
	    	return ADIABATIC;
	    }
	    else if ((strncmp(s,"HS",2) == 0) || (strncmp(s,"hs",2) == 0))
	    {
	    	set_stratified_state(NULL);
	    	return HYDRO_STATIC;
	    }
	    else
	    {
	    	set_stratified_state(isothermal_stratified_state);
	    	return ISOTHERMAL;
	    }
	}
	else
	    return CONSTANT;
}		/*end prompt_for_stratification*/

EXPORT	void get_state_in_stratified_region(
	STRATIFICATION_TYPE s_type,
	Locstate	    state,
	double		    dz,
	double		    g_z,
	Locstate	    ref_st)
{
	debug_print("init_states","Entered get_state_in_stratified_region()\n");
	switch (s_type)
	{
	case CONSTANT:
	    set_state(state,GAS_STATE,ref_st);
	    return;
	case CONSTANT_DENSITY:
	case HYDRO_STATIC:
	    constant_density_stratified_state(state,dz,g_z,ref_st);
	    break;
	case ISOTHERMAL:
	    isothermal_stratified_state(state,dz,g_z,ref_st);
	    break;
	case ADIABATIC:
	    isentropic_stratified_state(state,dz,g_z,ref_st);
	    break;
	default:
	    screen("ERROR in get_state_in_stratified_region(), "
	           "Unrecoganized stratification type.");
	    clean_up(ERROR);
	}
}		/*end get_state_in_stratified_region*/

EXPORT	void	get_state_rt_kh_perturbed(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int             stype)
{
	int		prob_type = ct->type;

	debug_print("init_states","Entered get_state_rt_kh_perturbed()\n");
	if ((prob_type != RT_PERTURBED) && (prob_type != KH_PERTURBED))
	{
	    screen("ERROR in get_state_rt_kh_perturbed(), "
	           "inconsistent comp_type->type\n");
	    clean_up(ERROR);
	}

	if (prob_type == RT_PERTURBED)
	{
	    if (Rt_perturbed(ct)->rt_kh->stratification_type == HYDRO_STATIC)
	    { 
		LAYER_SURF	*lower_surf = Rt_perturbed(ct)->lower_surf;
	        LAYER_SURF	*upper_surf = Rt_perturbed(ct)->upper_surf;
	        int		dim = ct->params->dim;
	        const double	*grav;
		double           g_z;
		double           P, U;
		double           z = coords[dim-1];
		Locstate        ref_state = Rt_perturbed(ct)->rt_kh->ref_state;
	        grav = gravity(coords,initial_time(init));
	        g_z = grav[dim-1];
		U = computational_grid(intfc)->GU[dim-1];
	        set_type_of_state(s,TGAS_STATE);
		zero_state_velocity(s,dim);
		Dens(s) = Dens(ref_state);
		Params(s) = Params(ref_state);
		if (lower_surf != NULL && lower_surf->wv_type == CONTACT) 
					/* Above the contact */
		{
		    double rho;
		    P = pressure(ref_state);
		    rho = Dens(ref_state);
		    Press(s) = P + rho*g_z*(z-U);
		}
		else if (upper_surf != NULL && upper_surf->wv_type == CONTACT) 
					/* Below the contact */
		{
		    Locstate  sa;
		    COMP_TYPE *ca;
		    FOURIER_POLY *fpoly = upper_surf->fpoly;
		    double rhoa, rhob, Pa;
		    double zi;
		    if (ct->comp == upper_surf->l_comp)
			ca = comp_type(upper_surf->r_comp);
		    else if (ct->comp == upper_surf->r_comp)
			ca = comp_type(upper_surf->l_comp);
		    else
		    {
			screen("ERROR in get_state_rt_kh_perturbed() - "
				"inconsistent components\n");
			clean_up(ERROR);
		    }
		    if (ca->type != prob_type)
		    {
			screen("ERROR in get_state_rt_kh_perturbed() - "
				"unsupported case of different component "
				"types\n");
			clean_up(ERROR);
		    }
		    rhob = Dens(ref_state);
		    sa = Rt_perturbed(ca)->rt_kh->ref_state;
		    rhoa = Dens(sa);
		    Pa = pressure(sa);
		    if (fpoly != NULL)
			zi = fourier_poly(coords,fpoly);
		    else
			zi = upper_surf->pbar[dim-1];
		    Press(s) = Pa + rhoa*g_z*(zi - U) + rhob*g_z*(z-zi);
		}
		else
		{
		    screen("ERROR in get_state_rt_kh_perturbed() - "
				"unsupported option for HYDRO_STATIC "
				"initialization\n");
		    clean_up(ERROR);
		}
	    }
	    else
	    {
	        POINTER		extra = ct->extra;
	        ct->extra = (POINTER) Rt_perturbed(ct)->rt_kh;
	        get_state_rt_kh(coords,s,ct,hs,intfc,init,TGAS_STATE);
	        ct->extra = extra;
	        rt_add_to_state(coords,s,ct,init);
	    }
	}
	else if (prob_type == KH_PERTURBED)
	    kh_add_to_state(coords, s, ct);

	if (Dens(s) < 0)
	{
	    int	dim = ct->params->dim;

	    screen("\nERROR in get_state_rt_kh_perturbed(), ");
	    print_general_vector("At point ",coords,dim,"");
	    screen(", the density of the state is less than 0.\n");
	    clean_up(ERROR);
	}
	set_state(s,stype,s);
}		/*end get_state_rt_kh_perturbed*/

#if defined(TWOD) || defined(THREED)

/*ARGSUSED*/
LOCAL	void	get_state_rt_kh(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int	        stype)
{
	int		    dim = ct->params->dim;
	STRATIFICATION_TYPE s_type = Rt_kh(ct)->stratification_type;
	double		    *ref_c = Rt_kh(ct)->ref_coords;
	Locstate	    ref_st = Rt_kh(ct)->ref_state;
	double		    z, z_r, g_z;
	const double	    *grav;

	debug_print("init_states","Entered get_state_rt_kh()\n");
	z = coords[dim-1]; 
	z_r = ref_c[dim-1];
	grav = gravity(coords,initial_time(init));
	g_z = grav[dim-1];
	get_state_in_stratified_region(s_type,state,z-z_r,g_z,ref_st);
	if (Dens(state) < 0)
	{
	    screen("ERROR in get_state_rt_kh().\n");
	    fprint_general_vector(stderr,"At point ",coords,dim,"");
	    print_general_vector("At point ",coords,dim,"");
	    screen(", the density of the state is less than 0.\n");
	    clean_up(ERROR);
	}
	if (Press(state) < 0.0)
	{
	    (void) printf("WARNING in get_state_rt_kh().\n");
	    print_general_vector("At point ",coords,dim,"");
	    (void) printf(", the pressure of the state is less than 0.\n");
	}
	set_state(state,stype,state);
}		/*end get_state_rt_kh*/


LOCAL	void	rt_add_to_state(
	double		*coords,
	Locstate	state, /* assumed to be TGAS_STATE */
	COMP_TYPE	*ct,
	INIT_DATA	*init)
{
	static const double ACC = 1.0e-14;		/*TOLERANCE*/

	int		i, j, k0, k1, dim, layer_label, num_modes, nstep;
	double		rho, csq, rho_prime, a_z, b_z, k_dot_r, sig, phase;
	double		tmp, ub, lb, z0, z1, a0, a1, b0, b1, z, h_z, g_z;
	const double	*grav;
	_RT_PERTURBED	*rtp = Rt_perturbed(ct);
	NORMAL_MODE	*n_m;

	dim = ct->params->dim;
	layer_label = rtp->layer_label;
	nstep = rtp->lin_pert_intvl;
	num_modes = rtp->num_modes;

	if (num_modes <= 0)	return;

	z0 = lb = get_surf_height(coords, rtp->lower_surf);
	z1 = ub = get_surf_height(coords, rtp->upper_surf);
	h_z = (ub-lb)/nstep;
	z = coords[dim-1];
	grav = gravity(coords,initial_time(init));
	g_z = grav[dim-1];

	/* find out which interval z belongs */
	if (ub < lb)
	{
	    screen("\nERROR in rt_add_to_state(), ");
	    screen("Interfaces are already tangled.\n");
	    (void) printf("ub = %g, lb = %g\n",ub,lb);
	    print_interface(current_interface());
	    clean_up(ERROR);
	}
	if (z >= ub)
	    k0 = k1 = nstep;
	else if (z <= lb)
	    k0 = k1 = 0;
	else
	{
	    k0 = 0,		k1 = nstep;
	    for (;;)
	    {
	    	tmp = (z0+z1)/2.0;
	    	if (z >= tmp)
	    	    k0 = (k0+k1)/2, 	z0 = tmp;
	    	else
	    	    k1 = (k0+k1)/2, 	z1 = tmp;
	    	if (k1-k0 == 1)	break;
	    }
	}

	z0 = lb+k0*h_z;	
	z1 = lb+k1*h_z;	
	rho = Dens(state);
	csq = sound_speed_squared(state);
	rho_prime = get_rho_prime(state,ct,g_z);
	/* sum over all the normal modes */
	for (i = 0; i < num_modes; ++i)
	{
	    n_m = rtp->normal_mode[i];
	    if (k0 == k1)
	    {
	    	a_z = n_m->a[layer_label][k1];
	    	b_z = n_m->b[layer_label][k1];
	    }
	    else
	    {
	    	a0 = n_m->a[layer_label][k0];
	    	a1 = n_m->a[layer_label][k1];
	    	a_z = ((a1-a0)*z+(a0*z1-a1*z0))/h_z;
	    	b0 = n_m->b[layer_label][k0];
	    	b1 = n_m->b[layer_label][k1];
	    	b_z = ((b1-b0)*z+(b0*z1-b1*z0))/h_z;
	    }
	    if ((fabs(a_z) < ACC) && (fabs(b_z) < ACC))	continue;

	    sig = n_m->sigma_r;
	    phase = n_m->phase;
	    k_dot_r = 0.0;
	    for (j = 0; j <= dim-2; ++j)
	    	k_dot_r += n_m->wv_num[j]*coords[j];
	    phase += k_dot_r;
	    Dens(state) += ((g_z*rho/csq-rho_prime)*a_z+b_z/csq)/sig*sin(phase);
	    Press(state) += b_z/sig*sin(phase);
	    Vel(state)[dim-1] += a_z*sin(phase);
	    set_type_of_state(state,TGAS_STATE);
	    phase += 3.0*PI/2;
	    for (j = 0; j <= dim-2; ++j)
	        Vel(state)[j] += sin(phase)*n_m->wv_num[j]*b_z/(rho*sig*sig);
            reset_gamma(state);
	}
}		/*end rt_add_to_state*/


/*ARGSUSED*/
LOCAL	void	kh_add_to_state(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*ct)
{
	screen("\nERROR in kh_add_to_state(), "
	       "This function has not been written yet.\n");
	clean_up(ERROR);
}		/*end kh_add_to_state*/

EXPORT	void	set_rt_kh_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	_RT_KH		*extra;

	if (comp_type->type == RT_KH) /*ALREADY SET*/
	    return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = RT_KH;
	comp_type->_get_state = get_state_rt_kh;
	scalar(&extra,sizeof(_RT_KH));
	extra->stratification_type = CONSTANT;
	alloc_state(front->interf,&extra->ref_state,front->sizest);
	comp_type->extra = (POINTER)extra;
	comp_type->free_comp_type_extra = free_rt_kh_comp_type;
}		/*end set_rt_kh_comp_type*/

EXPORT	void	free_rt_kh_comp_type(
	COMP_TYPE	*comp_type)
{
	_RT_KH		*extra;

	if (comp_type->type != RT_KH)
		return;
	extra = Rt_kh(comp_type);
	if (extra == NULL)
		return;
	if (extra->ref_state != NULL)
		free(extra->ref_state);
	free(extra);
	comp_type->extra = NULL;
}		/*end free_rt_kh_comp_type*/
#endif /* defined(TWOD) || defined(THREED) */
