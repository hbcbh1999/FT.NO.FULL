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
*				giniteos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains equation of state dependent functions needed for 
*	initialization and restart.
*
*	Functions for initialization
*
*		prompt_for_polytropic
*
*/

#include <geos/geosdecs.h>

	/* LOCAL Function Declarations */
LOCAL	void	add_params_to_list(Gas_param*);
LOCAL	void	set_default_params_tolerances(Gas_param*);

LOCAL int num_params = 0;
LOCAL Gas_param **params_list = NULL;

LOCAL	void add_params_to_list(
	Gas_param	*params)
{
	G_USER_INTERFACE *guh;
	int		i;
	Gas_param	**tmp_params;
	static const int	SMALL_BLOCK = 5;/*TOLERANCE*/
	static int	max_num_params;

	if (index_of_Gas_param(params) != -1) /* Already in list */
	    return;
	if (params_list == NULL)
	{
	    max_num_params = SMALL_BLOCK;
	    uni_array(&params_list,max_num_params,sizeof(Gas_param *));
	    params_list[0] = NULL;
	}
	if (num_params >= (max_num_params-1))
	{
	    max_num_params += SMALL_BLOCK;
	    uni_array(&tmp_params,max_num_params,sizeof(Gas_param *));
	    for (i = 0; i < num_params; ++i)
	    	tmp_params[i] = params_list[i];
	    free(params_list);
	    params_list = tmp_params;
	    tmp_params = NULL;
	}
	params_list[num_params++] = params;
	params_list[num_params] = NULL;
	for (i = 1; i <= MAXD; ++i)
	{
	    guh = g_user_hook(i);
	    guh->num_params = num_params;
	    guh->params_list = params_list;
	}
}		/*end add_params_to_list*/

EXPORT	int return_params_list(
	Gas_param	***pparams_list)
{
	if (pparams_list != NULL)
	    *pparams_list = params_list;
	return num_params;
}		/*end return_params_list*/

EXPORT	boolean	valid_params(
	Locstate	st)
{
	Gas_param	*params = Params(st);
	int i;

	if (params == NULL) return YES;

	for (i = 0; i < num_params; ++i)
	    if (params == params_list[i])
		return YES;

	return NO;
}		/*end valid_params*/

EXPORT	uint64_t gas_param_number(
	Gas_param       *params)
{
	if (debugging("addresses"))
	    return ptr2ull(params);
	else if (params == NULL)
	    return 0; /*Obstacle state is address 0*/
	else
	{
	    Gas_param	**prms;
	    uint64_t	i;
	    (void) return_params_list(&prms);
	    for (i = 1; prms && *prms; ++i, ++prms)
	    	if (params == *prms)
		    return i;
	}
	screen("ERROR in gas_param_number(), can not find params index for "
	       "params %llu\n",ptr2ull(params));
	clean_up(ERROR);
	return (uint64_t)(-1);
}		/*end gas_param_number*/

EXPORT	int index_of_Gas_param(
	Gas_param	*params)
{
	Gas_param	**prms_lst;
	int		i, nprms;

	nprms = return_params_list(&prms_lst);

	for (i = 0; i < nprms; ++i)
		if (params == prms_lst[i]) return i;

	return -1;
}		/*end index_of_Gas_param*/

EXPORT	Gas_param	*Params_of_index(
	int	params_index)
{
	Gas_param	**prms_lst;

	if ((params_index<0) || (params_index>=return_params_list(&prms_lst)))
	    return NULL;
	return prms_lst[params_index];
}		/*end Params_of_index*/

EXPORT	void set_params_list(
	INTERFACE	*intfc)
{
	num_gas_params(intfc) = num_params;
	gas_params_list(intfc) = params_list;
}		/*end set_params_list*/



LOCAL	void set_default_params_tolerances(
	Gas_param	*params)
{
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	const double eps = MACH_EPS; /*TOLERANCE Machine epsilon*/

	params->min_energy = eps;	/*TOLERANCE*/;
	params->min_pressure = eps; /*TOLERANCE*/
				 /* Minimum pressure appearing in */
				 /* non-combustion Riemann find_mid_state */
				 /* and in certain initializations */
	params->vacuum_dens = eps; /*TOLERANCE vacuum density cutoff*/
	params->raref_press = 1.0 - eps; /*TOLERANCE*/

#if defined(COMBUSTION_CODE)
	params->tol_alpha = eps;	/* TOLERANCE how close to zero  */
					/* pressure ratio may get	*/
	params->tol_press = eps;
#endif /* defined(COMBUSTION_CODE) */

#else /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	params->min_energy = -HUGE_VAL;
	params->min_pressure = -HUGE_VAL;
	params->vacuum_dens = -HUGE_VAL;
	params->raref_press = 1.0;

#if defined(COMBUSTION_CODE)
	params->tol_alpha = 0.0;
	params->tol_press = 0.0;
#endif /* defined(COMBUSTION_CODE) */

#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

}		/*end set_default_params_tolerances*/

EXPORT	Gas_param *alloc_Gas_param(
	size_t	size)
{
	Gas_param	*params;

	size = max(size,sizeof(Gas_param));
	scalar(&params,size);
	params->prt_label = "Gas_param";
	params->n_comps = 1;/*DEFAULT*/
	params->_fprint_state = g_fprint_state;
	params->_verbose_fprint_state = g_verbose_fprint_state;
	params->_fprint_Gas_param = g_fprint_Gas_param;
	params->_set_state = g_set_state;
	params->_set_params = g_set_params;
	params->_alloc_state = g_alloc_state;
	params->_alloc_intfc_state = g_alloc_intfc_state;
	params->_clear_state = g_clear_state;
	params->_invalid_state = g_invalid_state;
	set_default_params_tolerances(params);
	add_params_to_list(params);

	return params;
}		/*end alloc_Gas_param*/

EXPORT  void	g_set_params(
	Locstate	st1,
	Locstate	st2)
{
	Params(st1) = Params(st2);
}               /* g_set_params */



