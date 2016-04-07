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
*				gintfcpert.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if defined(TWOD) || defined(THREED)

#include <gdecs/gdecs.h>

/*
*			pert_interface():
*/

EXPORT double pert_interface(
	SINE_PERT	*pert,
	double		*coords,
	double		t,
	int		dim)
{
	int		i, j;
	int		number_modes = pert->number_modes;
	double		arg, sigma;
	double		z;

	z = pert->z_intfc;
	for (i = 0; i < number_modes; i++) 
	{
	    arg = 0.0;
	    for (j = 0; j < dim-1; j++)
	    	arg += pert->mode[i].wave_number[j]*coords[j];
	    arg -= pert->mode[i].phase;
	    sigma = pert->mode[i].growth_rate;
	    z += pert->mode[i].amplitude * exp(sigma*t)*sin(arg);
	}
	return z;
}		/*end pert_interface*/


EXPORT	SINE_PERT	*alloc_sine_pert_structure(
	int	number_modes,
	Front	*front)
{
	SINE_PERT	*rsp;
	ALIGN		*rspalgn;
	size_t		size;
	size_t		nasp, naas, nam;

	nasp = num_aligns(sizeof(SINE_PERT));
	naas = num_aligns(front->sizest);
	nam = num_aligns(sizeof(MODE));
	size = sizeof(ALIGN)*(nasp+naas+number_modes*nam);
	scalar(&rspalgn,size);
	rsp = (SINE_PERT*)rspalgn;
	rspalgn += nasp;
	rsp->mode = (MODE*)rspalgn;
	rsp->number_modes = number_modes;
	rspalgn += number_modes*nam;
	rsp->amb_st = (Locstate)rspalgn;
	clear_state(front->interf,rsp->amb_st,front->sizest);
	return	rsp;
}		/*end alloc_sine_pert_structure*/

EXPORT	SINE_PERT	*alloc_and_copy_sine_pert_structure(
	SINE_PERT	*pert,
	Front		*front)
{
	SINE_PERT       *rsp;
	Locstate	amb_st;
	MODE		*mode;
	int		i;

	rsp = alloc_sine_pert_structure(pert->number_modes,front);
	amb_st = rsp->amb_st;
	mode = rsp->mode;
	*rsp = *pert;
	rsp->amb_st = amb_st;
	rsp->mode = mode;
	copy_state(rsp->amb_st,pert->amb_st);
	for (i = 0; i < pert->number_modes; i++)
		rsp->mode[i] = pert->mode[i];
	return rsp;
}		/*end alloc_and_copy_sine_pert_structure*/
#endif /* defined(TWOD) || defined(THREED) */
