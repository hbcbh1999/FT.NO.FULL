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
*				geosvec.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains generic equation of state independent functions to load
*	equation of state information into Vec_Gas arrays.
*/

#include <geos/geosdecs.h>


/* LOCAL Function Prototypes */
LOCAL	Locstate	load_tmpst(double,double,Locstate);

EXPORT	void	set_params_jumps(
	Vec_Gas		*vst,
	int		imin,
	int		imax)
{
	Locstate	*state = vst->state;
	int		*prms_jmp = vst->prms_jmp;
	int		i, j, nprms;

	if (Vec_Gas_field_set(vst,prms_jmp))
	    return;

	nprms = 0;
	j = imin;
	while (j < imax)
	{
	    prms_jmp[nprms++] = j;
	    for (i = j+1; i < imax; ++i)
	    	if (Different_params(state[i],state[j]))
		    break;
	    j = i;
	}
	prms_jmp[nprms] = imax;
	Vec_Gas_field_set(vst,prms_jmp) = YES;

	vst->nprms = nprms;
}		/*end set_params_jumps*/

/*
*			load_pressure_and_gammas():
*
*	Loads the pressure, adiabatic exponent, and Gruneisen
*	coefficient uni_arrays of the Vec_Gas state vst.
*	This function assumes that either the specific internal energy
*	uni_array vst->e or vst->re is already loaded. Upon return the
*	fields vst->p, vst->c2, and vst->GAM are set.  If in addition
*	vst->e or vst->re is allocated these will also be set upon
*	return;
*/

EXPORT void load_pressure_and_gammas(
	Vec_Gas		*vst,
	int		offset,
	int		vsize)
{
	Locstate *state = vst->state + offset;
	boolean     e_tmp_field = NO;
	double    *p = vst->p + offset;
	double    *c2 = vst->c2 + offset;
	double    *GAM = vst->GAM + offset;
	double    *rho = vst->rho + offset;
	double    *e = vst->e + offset;
	double    *FD;
	double    c2tmp;
	int	 *prms_jmp = vst->prms_jmp + offset;
	int	 nprms;
	int	 i, j, k;

	FD = (vst->FD != NULL) ? vst->FD + offset : NULL;

	if (Vec_Gas_field_set(vst,p) &&
	    Vec_Gas_field_set(vst,c2) &&
	    Vec_Gas_field_set(vst,GAM) &&
	    ((vst->FD == NULL) || Vec_Gas_field_set(vst,FD)))
	    return;

#define	single_eos_load_pressure_and_gammas(vst,offset,vsize)		\
	(*Params(vst->state[offset])->eos->_single_eos_load_pressure_and_gammas)(vst,offset,vsize)

	if (!Vec_Gas_field_set(vst,e))
	{
	    if (vst->e == NULL)
	    {			/* Use pressure field for temporary */
	        e_tmp_field = YES;
	        vst->e = vst->p;
	        Vec_Gas_field_allocated(vst,e) = YES;
	        Vec_Gas_field_set(vst,e) = NO;
	    }
	    load_specific_internal_energy(vst,offset,vsize);
	}

	nprms = 0;
	for (i = 0, j = 0; j < vsize; j = i)
	{
	    prms_jmp[nprms++] = j+offset;
	    for (i = j+1; i < vsize; ++i)
	    	if (Different_params(state[i],state[j])) break;
	    if (is_obstacle_state(state[j]))
	    {
	    	for (k = 0; k < i-j; ++k)
	    	{
	    	    p[j+k] = 0.0;
	            c2[j+k] = 0.0;
		    GAM[j+k] = 0.0;
		}
		if (FD != NULL)
		{
	    	    for (k = 0; k < i-j; ++k)
		        FD[j+k] = 0.0;
		}
	    }
	    else
	    {
		single_eos_load_pressure_and_gammas(vst,j+offset,i-j);
		if ((c2[j]<0.0) && (j>0) && !is_obstacle_state(state[j-1]))
		{
		    Locstate tmpst = load_tmpst(rho[j],e[j],state[j-1]);
		    if ((c2tmp = sound_speed_squared(tmpst)) > 0.0)
		    {
		        c2[j] = c2tmp;
			p[j] = pressure(tmpst);
			GAM[j] = gruneisen_gamma(tmpst);
			if (FD != NULL)
			    FD[j] = fundamental_derivative(tmpst);
		    }
		}
		if ((c2[i-1]<0.0) && (i<vsize) &&
		    !is_obstacle_state(state[i]))
		{
		    Locstate tmpst = load_tmpst(rho[i-1],e[i-1],state[i]);
		    if ((c2tmp = sound_speed_squared(tmpst)) > 0.0)
		    {
		        c2[i-1] = c2tmp;
			p[i-1] = pressure(tmpst);
			GAM[i-1] = gruneisen_gamma(tmpst);
			if (FD != NULL)
			    FD[i-1] = fundamental_derivative(tmpst);
		    }
		}
	    	for (k = 0; k < i-j; ++k)
		    if (c2[j+k] < 0.0)
		        c2[j+k] = MACH_EPS;
	    }
	}
	Vec_Gas_field_set(vst,p) = YES;
	Vec_Gas_field_set(vst,c2) = YES;
	Vec_Gas_field_set(vst,GAM) = YES;
	if (FD != NULL)
	    Vec_Gas_field_set(vst,FD) = YES;
	prms_jmp[nprms] = vsize+offset;
	vst->nprms = nprms;
	Vec_Gas_field_set(vst,prms_jmp) = YES;
	if (e_tmp_field)
	{
	    vst->e = NULL;
	    Vec_Gas_field_allocated(vst,e) = NO;
	    Vec_Gas_field_set(vst,e) = NO;
	}

#undef	single_eos_load_pressure_and_gammas

}                       /*end load_pressure_and_gammas*/

EXPORT void load_pressure_and_sound_speed(
	Vec_Gas		*vst,
	int		offset,
	int		vsize)
{
	boolean     e_tmp_field = NO;
	Locstate *state = vst->state + offset;
	double    *p = vst->p + offset;
	double    *rho = vst->rho + offset;
	double    *e = vst->e + offset;
	double    *c2 = vst->c2 + offset, *c = vst->c + offset;
	double    c2tmp;
	int	 *prms_jmp = vst->prms_jmp + offset;
	int	 nprms;
	int	 i, j, k;

	if (Vec_Gas_field_set(vst,p) && Vec_Gas_field_set(vst,c))
	    return;

	for (i = 0; i < vsize; ++i)
	    if (Params(state[i]) != NULL)
	        break;

	Vec_Gas_field_set(vst,p) = YES;
	Vec_Gas_field_set(vst,c) = YES;
	Vec_Gas_field_set(vst,c2) = YES;
	if (i == vsize)
	{
	    prms_jmp[0] = offset;
	    prms_jmp[1] = vsize + offset;
	    vst->nprms = 1;
	    for (i = 0; i < vsize; ++i)
	    {
	    	p[i] = 0.0;
	    	c[i] = c2[i] = 0.0;
	    }
	    return;
	}

	if (!Vec_Gas_field_set(vst,e))
	{
	    if (vst->e == NULL)
	    {			/* Use pressure field for temporary */
	        e_tmp_field = YES;
	        vst->e = vst->p;
	        Vec_Gas_field_allocated(vst,e) = YES;
	        Vec_Gas_field_set(vst,e) = NO;
	    }
	    load_specific_internal_energy(vst,offset,vsize);
	}

#define	single_eos_load_pressure_and_sound_speed2(vst,offset,vsize)	\
	(*Params(vst->state[offset])->eos->_single_eos_load_pressure_and_sound_speed2)(vst,offset,vsize)

	nprms = 0;
	j = 0;
	for (j = 0; j < vsize; j = i)
	{
	    prms_jmp[nprms++] = j+offset;
	    for (i = j+1; i < vsize; ++i)
	        if (Different_params(state[i],state[j])) break;
	    if (is_obstacle_state(state[j]))
	    {
	        for (k = 0; k < i-j; ++k)
	        {
	            p[j+k] = 0.0;
	            c[j+k] = c2[j+k] = 0.0;
	        }
	    }
	    else
	    {
	        single_eos_load_pressure_and_sound_speed2(vst,j+offset,i-j);
		if ((c2[j]<0.0) && (j>0) && !is_obstacle_state(state[j-1]))
		{
		    Locstate tmpst = load_tmpst(rho[j],e[j],state[j-1]);
		    if ((c2tmp = sound_speed_squared(tmpst)) > 0.0)
		    {
		        c2[j] = c2tmp;
			p[j] = pressure(tmpst);
		    }
		}
		if ((c2[i-1]<0.0) && (i<vsize) &&
		    !is_obstacle_state(state[i]))
		{
		    Locstate tmpst = load_tmpst(rho[i-1],e[i-1],state[i]);
		    if ((c2tmp = sound_speed_squared(tmpst)) > 0.0)
		    {
		        c2[i-1] = c2tmp;
			p[i-1] = pressure(tmpst);
		    }
		}
	        for (k = 0; k < i-j; ++k)
		{
		    if (c2[j+k] > 0.0)
		        c[j+k] = sqrt(c2[j+k]);
		    else
	                c[j+k] = c2[j+k] = MACH_EPS;
		}
	    }
	}

	if (e_tmp_field) /* Use pressure field for temporary */
	{
	    vst->e = NULL;
	    Vec_Gas_field_allocated(vst,e) = NO;
	    Vec_Gas_field_set(vst,e) = NO;
	}

	prms_jmp[nprms] = vsize+offset;
	vst->nprms = nprms;
	Vec_Gas_field_set(vst,prms_jmp) = YES;

#undef	single_eos_load_pressure_and_sound_speed2
}			/*end load_pressure_and_sound_speed*/


EXPORT void load_pressure(
	Vec_Gas		*vst,
	int		offset,
	int		vsize)
{
	boolean     e_tmp_field = NO;
	Locstate *state = vst->state + offset;
	double    *p = vst->p + offset;
	int	 i, j, k;
	int	 *prms_jmp = vst->prms_jmp + offset;
	int	 nprms;

	if (Vec_Gas_field_set(vst,p))
	    return;

#define	single_eos_load_pressure(vst,offset,vsize)		\
	(*Params(vst->state[offset])->eos->_single_eos_load_pressure)(vst,offset,vsize)

	if (!Vec_Gas_field_set(vst,e))
	{
	    if (vst->e == NULL)
	    {			/* Use pressure field for temporary */
	        e_tmp_field = YES;
	        vst->e = vst->p;
	        Vec_Gas_field_allocated(vst,e) = YES;
	        Vec_Gas_field_set(vst,e) = NO;
	    }
	    load_specific_internal_energy(vst,offset,vsize);
	}

	nprms = 0;
	j = 0;
	for (i = 0, j = 0; j < vsize; j = i)
	{
	    prms_jmp[nprms++] = j+offset;
	    for (i = j+1; i < vsize; ++i)
	    	if (Different_params(state[i],state[j]))
		    break;
	    if (is_obstacle_state(state[j]))
	    {
	    	for (k = 0; k < i-j; ++k)
	    	    p[j+k] = 0.0;
	    }
	    else
	    	single_eos_load_pressure(vst,j+offset,i-j);
	}
	Vec_Gas_field_set(vst,p) = YES;
	prms_jmp[nprms] = vsize+offset;
	vst->nprms = nprms;

#undef	single_eos_load_pressure

	if (e_tmp_field) /* Use pressure field for temporary */
	{
	    vst->e = NULL;
	    Vec_Gas_field_allocated(vst,e) = NO;
	    Vec_Gas_field_set(vst,e) = NO;
	}
}                       /*end load_pressure*/

EXPORT void load_sound_speed(
	Vec_Gas		*vst,
	int		offset,
	int		vsize)
{
	boolean     e_tmp_field = NO;
	Locstate *state = vst->state + offset;
	double    *c = vst->c + offset, *c2 = vst->c2 + offset;
	double    *rho = vst->rho + offset, *e = vst->e + offset;
	double    c2tmp;
	int	 i, j, k;
	int	 *prms_jmp = vst->prms_jmp + offset;
	int	 nprms;

	if (Vec_Gas_field_set(vst,c))
	    return;

	if (!Vec_Gas_field_set(vst,e))
	{
	    if (vst->e == NULL)
	    {			/* Use pressure field for temporary */
	        e_tmp_field = YES;
	        vst->e = vst->p;
	        Vec_Gas_field_allocated(vst,e) = YES;
	        Vec_Gas_field_set(vst,e) = NO;
	    }
	    load_specific_internal_energy(vst,offset,vsize);
	}

#define	single_eos_load_sound_speed2(vst,offset,vsize)			\
	(*Params(vst->state[offset])->eos->_single_eos_load_sound_speed2)(vst,offset,vsize)

	nprms = 0;
	j = 0;
	for (i = 0, j = 0; j < vsize; j = i)
	{
	    prms_jmp[nprms++] = j+offset;
	    for (i = j+1; i < vsize; ++i)
	    	if (Different_params(state[i],state[j])) break;
	    if (is_obstacle_state(state[j]))
	    {
	    	for (k = 0; k < i-j; ++k)
	    	    c[j+k] = c2[j+k] = 0.0;
	    }
	    else
	    {
	    	single_eos_load_sound_speed2(vst,j+offset,i-j);
		if ((c2[j]<0.0) && (j>0) && !is_obstacle_state(state[j-1]))
		{
		    Locstate tmpst = load_tmpst(rho[j],e[j],state[j-1]);
		    if ((c2tmp = sound_speed_squared(tmpst)) > 0.0)
		        c2[j] = c2tmp;
		}
		if ((c2[i-1]<0.0) && (i<vsize) &&
		    !is_obstacle_state(state[i]))
		{
		    Locstate tmpst = load_tmpst(rho[i-1],e[i-1],state[i]);
		    if ((c2tmp = sound_speed_squared(tmpst)) > 0.0)
		        c2[i-1] = c2tmp;
		}

	    	for (k = 0; k < i-j; ++k)
		{
		    if (c2[j+k] > 0.0)
		        c[j+k] = sqrt(c2[j+k]);
		    else
	                c[j+k] = c2[j+k] = MACH_EPS;
		}
	    }
	}
	Vec_Gas_field_set(vst,c) = YES;
	Vec_Gas_field_set(vst,c2) = YES;
	prms_jmp[nprms] = vsize+offset;
	vst->nprms = nprms;
	Vec_Gas_field_set(vst,prms_jmp) = YES;

#undef	single_eos_load_sound_speed2

	if (e_tmp_field) /* Use pressure field for temporary */
	{
	    vst->e = NULL;
	    Vec_Gas_field_allocated(vst,e) = NO;
	    Vec_Gas_field_set(vst,e) = NO;
	}

}                       /*end load_sound_speed*/

#if !defined(UNRESTRICTED_THERMODYNAMICS)
LIB_LOCAL	void	limit_pressure(
	double	*p,
	double	*min_pressure,
	int     vsize)
{
	int k;

	for (k = 0; k < vsize; ++k)
	    if (p[k] < min_pressure[k])
	        p[k] = min_pressure[k];
}		/*end limit_pressure*/
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

/*
*			load_specific_internal_energy():
*
*	Loads the internal energy per unit volume of the Vec_Gas state vst
*	into the vector re.
*/

EXPORT	void	load_specific_internal_energy(
	Vec_Gas *vst,
	int     offset,
	int     vsize)
{
	double     *rho, *e, *re;
	int       k;

	if (Vec_Gas_field_set(vst,e))
	    return;

	rho = vst->rho + offset;
	e = vst->e + offset;
	if (!Vec_Gas_field_set(vst,re))
	    load_internal_energy_density(vst,offset,vsize);

	re = vst->re + offset;
	for (k = 0; k < vsize; ++k)
	    e[k] = re[k]/rho[k];
	Vec_Gas_field_set(vst,e) = YES;
}		/*end load_specific_internal_energy*/

/*
*			load_internal_energy_density():
*
*	Loads the internal energy per unit volume of the Vec_Gas state vst
*	into the vector re.
*/

EXPORT	void	load_internal_energy_density(
	Vec_Gas *vst,
	int     offset,
	int     vsize)
{
	Gas_param *params;
	Locstate *state;
	double     *rho, *en_den, *mx, *my, *mz, *e, *re;
	int       k;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double      *min_energy = vst->min_energy + offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	if (Vec_Gas_field_set(vst,re))
	    return;

	params = Params(vst->state[offset]);
	rho = vst->rho + offset;

	re = vst->re + offset;
	if (Vec_Gas_field_set(vst,e))
	{
	    e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	        re[k] = rho[k]*e[k];
	    Vec_Gas_field_set(vst,re) = YES;
	    return;
	}

	en_den = vst->en_den + offset;
	state = vst->state + offset;

	params = NULL;
	for (k = 0; k < vsize; ++k)
	{
	    if (!is_obstacle_state(state[k]))
	    {
	        params = Params(vst->state[offset]);
		break;
	    }
	}
	if (params == NULL)
	{
	    for (k = 0; k < vsize; ++k)
	    	re[k] = 0.0;
	}
	else
	{
	    switch (params->dim)
	    {
	    case 1:
	        mx = vst->m[0] + offset;
	        for (k = 0; k < vsize; ++k)
	    	    re[k] = en_den[k] - 0.5*sqr(mx[k])/rho[k];
	        break;
	    case 2:
	        mx = vst->m[0] + offset;
	        my = vst->m[1] + offset;
	        for (k = 0; k < vsize; ++k)
	            re[k] = en_den[k] - 0.5*(sqr(mx[k])+sqr(my[k]))/rho[k];
	        break;
	    case 3:
	        mx = vst->m[0] + offset;
	        my = vst->m[1] + offset;
	        mz = vst->m[2] + offset;
	        for (k = 0; k < vsize; ++k)
	            re[k] = en_den[k] -
		            0.5*(sqr(mx[k])+sqr(my[k])+sqr(mz[k]))/rho[k];
	        break;
	    default:
	        screen("ERROR in load_internal_energy_density(), "
		       "invalid dimension %d\n",params->dim);
		clean_up(ERROR);
	    }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    for (k = 0; k < vsize; ++k)
	        if (re[k] < min_energy[k])
		    re[k] = min_energy[k];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	}
	Vec_Gas_field_set(vst,re) = YES;
}		/*end load_specific_internal_energy*/

LOCAL	Locstate	load_tmpst(
	double    rho,
	double    e,
	Locstate state)
{
	static  Locstate tmpst = NULL;
	if (tmpst == NULL)
	    g_alloc_state(&tmpst,Params(state)->sizest);
	Set_params(tmpst,state);
	Dens(tmpst) = rho;
	Energy(tmpst) = e;
	set_type_of_state(tmpst,EGAS_STATE);
	return tmpst;
}		/*load_tmpst*/
