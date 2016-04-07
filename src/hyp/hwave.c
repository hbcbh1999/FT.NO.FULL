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
*				hwave.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines that interface to the grid construction
*	and solution function algorthm(s).  Also contains some
*	miscellaneous subroutines for the hyperbolic library.
*/


#include <hyp/hdecs.h>

	/* LOCAL Function prototypes */
LOCAL	Wave	*h_copy_wave(Wave*);
LOCAL	h_MaxWaveSpeed	*h_copy_MaxWaveSpeed(Wave*);
LOCAL	void	h_assign_copy_wave_pointers(Wave*,Wave*);
LOCAL	void	h_assign_wave_pointers(Wave*,Wave*);
LOCAL	void	h_clear_wave_pointers(Wave*);
LOCAL	void	h_destroy_MaxWaveSpeed(Wave*);
LOCAL	void	h_free_copy_wave_pointers(Wave*);
LOCAL	void	h_free_wave_pointers(Wave*);
LOCAL	void	h_include_max_wave_speed_info(h_MaxWaveSpeed*,Wave*);

/* For parabolic step */
LOCAL	h_MaxViscosity	*h_copy_MaxViscosity(Wave*);
LOCAL	void	h_destroy_MaxViscosity(Wave*);
LOCAL	void	h_include_max_viscosity_info(h_MaxViscosity*,Wave*);
LOCAL	void	h_initialize_max_viscosity(Wave*);

/*ARGSUSED*/
EXPORT	void	h_set_default_wave_parameters(
	INIT_DATA *init,
	Wave	  *wave)
{
	int	dim = wave->rect_grid->dim;
	F_USER_INTERFACE *fuh = f_user_hook(dim);

	wave->_scatter_states = h_scatter_states;
#if defined(__MPI__)
	if (use_blocking_pp_comp(init) == NO)
	    wave->_scatter_states = h_iscatter_states;
#endif /* defined(__MPI__) */
	wave->_copy_wave = h_copy_wave;
	wave->_copy_into_wave = h_copy_into_wave;
	wave->_free_wave = h_free_wave;
	wave->_assign_wave_parameters = h_assign_wave_parameters;
	wave->_assign_wave_pointers = h_assign_wave_pointers;
	wave->_clear_wave_pointers = h_clear_wave_pointers;
	wave->_assign_copy_wave_pointers = h_assign_copy_wave_pointers;
	wave->_free_wave_pointers = h_free_wave_pointers;
	wave->_free_copy_wave_pointers = h_free_copy_wave_pointers;
	wave->min_storage = YES;
	wave->show_tri_soln = NULL;
	wave->plot_hyp_soln = NULL;
	wave->_npt_solver = NULL;
	wave->_vec_solver = NULL;
	wave->_alloc_phys_vecs = NULL;
	wave->_free_phys_vecs = NULL;
	wave->max_hyp_time_step = h_max_hyp_time_step;
	wave->stat_prod_rate = NULL;
	wave->_alloc_state = fuh->_alloc_state;
	wave->_clear_state = fuh->_clear_state;
	wave->_obstacle_state = fuh->_obstacle_state;
	wave->_alloc_MaxWaveSpeed = h_alloc_MaxWaveSpeed;

	/* For parabolic step */
	MaxViscosity(wave) = h_alloc_MaxViscosity(MaxViscosity(wave),wave);
        /* new add 052103, allocate space only when viscosity is used */
        wave->_alloc_MaxViscosity = h_alloc_MaxViscosity;
        /* end 052103: copied from Mn-kim, for parabolic driver */
}		/*end h_set_default_wave_parameters*/

/*
*			h_assign_wave_parameters():
*
*	Copies the "parameters" at the top of a wave structure
*	from wave to newwave, and then initializes the "pointers"
*	at the top of newwave to NULL.
*/

EXPORT void h_assign_wave_parameters(
	Wave		*newwave,
	Wave		*wave)
{
	*newwave = *wave;
	clear_wave_pointers(newwave);
}		/*end h_assign_wave_parameters*/

LOCAL	void	h_clear_wave_pointers(
	Wave		*wave)
{
	zero_scalar(&wave_pointers(wave),sizeof(WAVE_POINTERS));
	wave_min_comp(wave) = NO_COMP;
}		/*end h_clear_wave_pointers*/
			

/*
*			h_assign_wave_pointers():
*
*	Frees the wave_pointers of the left wave structure
*	and ft_assigns to them the ones on the right.
*/

LOCAL void h_assign_wave_pointers(
	Wave		*newwave,
	Wave		*wave)
{
	free_wave_pointers(newwave);

	wave_pointers(newwave) = wave_pointers(wave);
}		/*end h_assign_wave_pointers*/

/*
*                       h_assign_copy_wave_pointers():
*
*       Frees the "pointers" at the bottom of the left wave structure
*       and ft_assigns to them the ones on the right.
*/

LOCAL void h_assign_copy_wave_pointers(
	Wave		*newwave,
	Wave		*wave)
{
	free_copy_wave_pointers(newwave);
 
	wave_pointers(newwave) = wave_pointers(wave);
}		/*end h_assign_copy_wave_pointers*/



/*
*			h_free_wave_pointers():
*
*	Frees the storage associated with the "pointer" variables
*	at the bottom of the wave structure.
*/

LOCAL void h_free_wave_pointers(
	Wave		*wave)
{
	if (wave_tri_soln(wave) != NULL)
	{
	    free_hyp_tri_grid(&wave_tri_soln(wave)->tri_grid);
	    free(wave_tri_soln(wave));
	}
	if (wave_areas(wave) != NULL)
	{
	    free(wave_areas(wave));
	}
	clear_wave_pointers(wave);
}		/*end h_free_wave_pointers*/

LOCAL void h_free_copy_wave_pointers(
	Wave		*wave)
{
	if (wave_tri_soln(wave) != NULL)
	{
		free_copy_tri_grid(&wave_tri_soln(wave)->tri_grid);
		free(wave_tri_soln(wave));
		wave_tri_soln(wave) = NULL;
	}
	if (wave_areas(wave) != NULL)
	{
		free(wave_areas(wave));
		wave_areas(wave) = NULL;
	}
}		/*end h_free_copy_wave_pointers*/

LOCAL	Wave	*h_copy_wave(
	Wave		*wave)
{
	Wave		*newwave;

	scalar(&newwave,sizeof(Wave));
	copy_into_wave(newwave,wave);
	return newwave;
}		/*end h_copy_wave*/

EXPORT	void	h_copy_into_wave(
	Wave	*newwave,
	Wave	*wave)
{
	*newwave = *wave;
}		/*end h_copy_into_wave*/

EXPORT	void	h_free_wave(
	Wave		*wave)
{
	free(wave);
}		/*end h_free_wave*/

EXPORT	h_MaxWaveSpeed	*h_alloc_MaxWaveSpeed(
	h_MaxWaveSpeed	*mxsp,
	Wave		*wave)
{
	byte **buf;
	if (mxsp != NULL)
	    return mxsp;

	scalar(&mxsp,sizeof(h_MaxWaveSpeed));
	mxsp->_sizest = wave->sizest;
	bi_array(&buf,MAXD,1,wave->sizest);
	mxsp->_mxspst = (Locstate*)buf;
	bi_array(&mxsp->_coords,MAXD,MAXD,FLOAT);

	mxsp->operators._set        = h_set_max_wave_speed;
	mxsp->operators._include    = h_include_max_wave_speed_info;
	mxsp->operators._initialize = h_initialize_max_wave_speed;
	mxsp->operators._print      = h_fprint_max_wave_speed_info;
	mxsp->operators._read_print = h_read_print_max_wave_speed_info;
	mxsp->operators._copy       = h_copy_MaxWaveSpeed;
	mxsp->operators._destroy    = h_destroy_MaxWaveSpeed;

	return mxsp;
}		/*end h_alloc_MaxWaveSpeed*/

EXPORT	void	h_set_max_wave_speed(
	int		i,
	double		spd,
	Locstate	state,
	double		*coords,
	Wave		*wave)
{
	if (fabs(spd) > Maxsp(wave)[i])
	{
	    int	j, dim = wave->rect_grid->dim;
	    Maxsp(wave)[i] = fabs(spd);
	    if (coords != NULL)
	        for (j = 0; j < dim; j++)
	            MaxWaveSpeedCoords(wave)[i][j] = coords[j];
	    if (state != NULL)
	        ft_assign(MaxWaveSpeedState(wave)[i],state,wave->sizest);
	}
}		/*end h_set_max_wave_speed*/

LOCAL	void	h_include_max_wave_speed_info(
	h_MaxWaveSpeed	*mxsp,
	Wave		*wave)
{
	int i, dim = wave->rect_grid->dim;

	for (i = 0; i < dim; i++)
	{
	    if (fabs(mxsp->_maxsp[i]) > Maxsp(wave)[i])
	    {
	    	int	j;

	    	Maxsp(wave)[i] = fabs(mxsp->_maxsp[i]);
	    	for (j = 0; j < dim; j++)
	    	    MaxWaveSpeedCoords(wave)[i][j] = mxsp->_coords[i][j];
	    	ft_assign(MaxWaveSpeedState(wave)[i],mxsp->_mxspst[i],
		       wave->sizest);
	    }
	}
}		/*end h_include_max_wave_speed_info*/

EXPORT	void	h_initialize_max_wave_speed(
	Wave	*wave)
{
	int i, j;
	for (i = 0; i < MAXD; i++)
	{
	    Maxsp(wave)[i] = 0.0;
	    for (j = 0; j < MAXD; j++)
	    	MaxWaveSpeedCoords(wave)[i][j] = HUGE_VAL;
	    (*wave->_clear_state)(MaxWaveSpeedState(wave)[i],wave->sizest);
	}
}		/*end h_initialize_max_wave_speed*/

LOCAL	h_MaxWaveSpeed	*h_copy_MaxWaveSpeed(
	Wave	*wave)
{
	h_MaxWaveSpeed	*mxsp;
	int	i, j, dim = wave->rect_grid->dim;
	
	mxsp = h_alloc_MaxWaveSpeed(NULL,wave);
	mxsp->operators = MaxWaveSpeedOperators(wave);
	mxsp->_sizest = MaxWaveSpeed(wave)->_sizest;
	for (i = 0; i < dim; i++)
	{
	    mxsp->_maxsp[i] = Maxsp(wave)[i];
	    ft_assign(mxsp->_mxspst[i],MaxWaveSpeedState(wave)[i],mxsp->_sizest);
	    for (j = 0; j < dim; j++)
	    	mxsp->_coords[i][j] = MaxWaveSpeedCoords(wave)[i][j];
	}
	return mxsp;
}		/*end h_copy_MaxWaveSpeed*/

LOCAL	void	h_destroy_MaxWaveSpeed(
	Wave	*wave)
{
	free(MaxWaveSpeedState(wave));
	free(MaxWaveSpeedCoords(wave));
	free(MaxWaveSpeed(wave));
	MaxWaveSpeed(wave) = NULL;
}		/*end h_destroy_MaxWaveSpeed*/

/* For parabolic step */
EXPORT	h_MaxViscosity	*h_alloc_MaxViscosity(
	h_MaxViscosity	*mxvisc,
	Wave		*wave)
{
	byte **buf;
	if (mxvisc != NULL)
	    return mxvisc;

	scalar(&mxvisc,sizeof(h_MaxViscosity));
	mxvisc->_sizest = wave->sizest;
	bi_array(&buf,MAXD,1,wave->sizest);
	mxvisc->_mxviscst = (Locstate*)buf;
	uni_array(&mxvisc->_coords,MAXD,FLOAT);

	mxvisc->operators._set        = h_set_max_viscosity;
	mxvisc->operators._include    = h_include_max_viscosity_info;
	mxvisc->operators._initialize = h_initialize_max_viscosity;
	mxvisc->operators._print      = h_fprint_max_viscosity_info;
	mxvisc->operators._read_print = h_read_print_max_viscosity_info;
	mxvisc->operators._copy       = h_copy_MaxViscosity;
	mxvisc->operators._destroy    = h_destroy_MaxViscosity;

	return mxvisc;
}		/*end h_alloc_MaxViscosity*/

EXPORT	void	h_set_max_viscosity(
	double		mu,
	Locstate	state,
	double		*coords,
	Wave		*wave)
{
	if (fabs(mu) > Maxvisc(wave))
	{
	    int	j, dim = wave->rect_grid->dim;
	    Maxvisc(wave) = fabs(mu);
	    if (coords != NULL)
	        for (j = 0; j < dim; j++)
	            MaxViscosityCoords(wave)[j] = coords[j];
	    if (state != NULL)
	        ft_assign(MaxViscosityState(wave),state,wave->sizest);
	}
}		/*end h_set_max_viscosity*/

LOCAL	void	h_include_max_viscosity_info(
	h_MaxViscosity	*mxvisc,
	Wave		*wave)
{
	int i, dim = wave->rect_grid->dim;

	if (fabs(mxvisc->_maxvisc) > Maxvisc(wave))
	{
	    int	j;
	    
	    Maxvisc(wave) = fabs(mxvisc->_maxvisc);
	    for (j = 0; j < dim; j++)
	        MaxViscosityCoords(wave)[j] = mxvisc->_coords[j];
	    ft_assign(MaxViscosityState(wave),mxvisc->_mxviscst,
		   wave->sizest);
	    
	}
}		/*end h_include_max_viscosity_info*/

LOCAL	void	h_initialize_max_viscosity(
	Wave	*wave)
{
	int j;
	Maxvisc(wave) = -HUGE_VAL;
	for (j = 0; j < MAXD; j++)
	  MaxViscosityCoords(wave)[j] = HUGE_VAL;
	(*wave->_clear_state)(MaxViscosityState(wave),wave->sizest);
}		/*end h_initialize_max_viscosity*/

LOCAL	h_MaxViscosity	*h_copy_MaxViscosity(
	Wave	*wave)
{
	h_MaxViscosity	*mxvisc;
	int	j, dim = wave->rect_grid->dim;
	
	mxvisc = h_alloc_MaxViscosity(NULL,wave);
	mxvisc->operators = MaxViscosityOperators(wave);
	mxvisc->_sizest = MaxViscosity(wave)->_sizest;
	mxvisc->_maxvisc = Maxvisc(wave);
	ft_assign(mxvisc->_mxviscst,MaxViscosityState(wave),
	       mxvisc->_sizest);
	for (j = 0; j < dim; j++)
	    mxvisc->_coords[j] = MaxViscosityCoords(wave)[j];
	return mxvisc;
}		/*end h_copy_MaxViscosity*/

LOCAL	void	h_destroy_MaxViscosity(
	Wave	*wave)
{
	free(MaxViscosityState(wave));
	free(MaxViscosityCoords(wave));
	free(MaxViscosity(wave));
	MaxViscosity(wave) = NULL;
}		/*end h_destroy_MaxViscosity*/
