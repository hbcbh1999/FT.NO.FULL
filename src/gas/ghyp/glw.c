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
*				glw.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for the Lax-Wendroff method.
*
*	LW(), and LWoblique() perform Lax-Wendroff
*	steps for a three-point stencil of states in specified directions.
*
*	All access to these functions is though function pointers
*	in the Wave and Front structures that are set by init_physics()
*	in ginit.c.
*/

#include <ghyp/ghyp.h>

	/* LOCAL Function Declarations */
LOCAL	double	compute_art_visc_corrections(Locstate,double,const double*,double*,
					     Locstate,double,const double*,double*,
					     Locstate,double,const double*,double*,
					     double**,int,Stencil*,Front*);
LOCAL	void	vec_Lax_Wendroff(int,double,double,int,Vec_Gas*,Vec_Gas*,
				 Vec_Gas*,double**,Wave*);

enum { TSWEEP = SMAXD+1 };


/*
*				LW():
*
*	Lax Wendroff finite difference scheme in one space variable
*	applied, in operator splitting to the time, iperm coordinates,
*	part of the conservation form of the equations of compressible flow.
*
*	IMPORTANT NOTE:
*	If both the stencil and the previous stencil are regular,
*	then s1 equates the previous stencil's s2 and s2 equates the
*	previous stencil's s3.  This implies that the mid state
*	sa equates the previously computed sb.  Since both these
*	structures are stored in static storage, it is possible
*	to obtain the current sa by coping the previous sb.
*	This should result in considerable computational savings
*	in regions away from the tracked front.
*/

EXPORT void LW(
	double		dh,
	double		dt,
	Locstate	ans,
	const double	*dir,
	int		swp_num,
	int		*iperm,
	int		*index,
	Stencil		*sten)
{
	double		lm, lmpg, mid_max_sp;
	COMPONENT	comp = sten->newcomp;
	Locstate	s1 = sten->st[-1], s2 = sten->st[0], s3 = sten->st[1];
#if defined(COMBUSTION_CODE)
	double		rt2, rta, rtb;
#endif /* defined(COMBUSTION_CODE) */
	double		k1, k2, k3;
	RECT_GRID	*gr = sten->fr->rect_grid;
	int		i, dim = gr->dim;
	int		idir = iperm[swp_num];
	double		time = sten->fr->time + 0.5*dt;
	static size_t	sizest = 0;
	static double	g2, g3, ga, gb;
	static double	crdsa[3], crdsb[3];
	static double	v1, pr1, v2, pr2, v3, pr3, va, pra, vb, prb;
	static boolean	is_grav;
	static double	**visc_corr = NULL;
	static Locstate sa = NULL,sb = NULL;
	double           ka, kb;
	double		rmin = pos_radius(0.0,gr);
	static double    Ra, Rb;
	static double    rad1, rad2, rad3;
	static GEOMETRY geom;
	static double    a;

	if (is_obstacle_state(s2))
	{
	    g_obstacle_state(ans,sten->fr->sizest);
	    sten->reg_stencil = NO;
	    return;
	}
	if (sb == NULL)
	{
	    sizest = Params(s2)->sizest;
	    (*Params(s2)->_alloc_state)(&sa,sizest);
	    (*Params(s2)->_alloc_state)(&sb,sizest);
	    set_type_of_state(sa,GAS_STATE);
	    set_type_of_state(sb,GAS_STATE);
	    bi_array(&visc_corr,2,3,FLOAT);
	    is_grav = is_gravity();
	    if (is_rotational_symmetry())
		geom = Geometry(&a);
	}

	set_type_of_state(ans,GAS_STATE);

#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND)
	    rt2 = 0.5 * reaction_rate(s2);
#endif /* defined(COMBUSTION_CODE) */

	if (sten->reg_stencil && sten->prev_reg_stencil)
	{
	    v1 = v2;	pr1 = pr2;
	    v2 = v3;	pr2 = pr3;	g2 = g3;
	    if (is_rotational_symmetry())
	    {
	    	rad1 = rad2;
	    	rad2 = rad3;
	    }
	}
	else
	{
	    v1 = vel(idir,s1);
	    pr1 = pressure(s1);
	    v2 = vel(idir,s2);
	    pr2 = pressure(s2);
	    g2 = 0.5*gravity(Coords(sten->p[0]),time)[idir];
	    if (is_rotational_symmetry())
	    {
	    	rad1 = pos_radius(Coords(sten->p[-1])[0],gr);
	    	rad2 = pos_radius(Coords(sten->p[0])[0],gr);
	    }
	}
	v3 = vel(idir,s3);
	pr3 = pressure(s3);
	if (is_rotational_symmetry())
	    rad3 = pos_radius(Coords(sten->p[1])[0],gr);
	g3 = 0.5*gravity(Coords(sten->p[1]),time)[idir];
	mid_max_sp = compute_art_visc_corrections(s1,v1,dir,Coords(sten->p[-1]),
			                          s2,v2,dir,Coords(sten->p[0]),
			                          s3,v3,dir,Coords(sten->p[1]),
			                          visc_corr,idir,sten,
						  sten->fr);
	set_max_wave_speed(idir,mid_max_sp,s2,Coords(sten->p[0]),sten->wave);

	lm = dt/dh;

	/* Compute state sa by applying Lax-Friedrichs to states s1 and s2 with
	   modifications due to artificial viscosity. */

	if (sten->reg_stencil && sten->prev_reg_stencil)
	{
	    ft_assign(sa,sb,sizest);
	    va = vb;
	    pra = prb;
	    if (is_rotational_symmetry())
	    	Ra = Rb;
	    ga = gb;
	    for (i = 0; i < dim; ++i)
		crdsa[i] = crdsb[i];
	}
	else
	{
	    if (is_rotational_symmetry())
	    {
	    	Ra = (geom == SPHERICAL) ?
	       	    2.0*(rad1*rad1+rad1*rad2+rad2*rad2)/(3.0*(rad1+rad2)) :
	       	    0.5*(rad1 + rad2);
	    }
	    for (i = 0; i < dim; ++i)
		crdsa[i] = 0.5*(Coords(sten->p[-1])[i] + Coords(sten->p[0])[i]);
	    ga = 0.5*gravity(crdsa,time)[idir];
	    lmpg = 0.5*(lm + visc_corr[0][2]);
	    k2 = 0.5*(1.0 - visc_corr[0][1]) - lmpg*v2;
	    k1 = 0.5*(1.0 + visc_corr[0][1]) + lmpg*v1;
	    if (is_rotational_symmetry() && 
		(a > 0.0) && (idir == 0) && (fabs(Ra) > fabs(rmin)))
	    {
		k2 -= 0.25*a*dt*v2/Ra;
		k1 -= 0.25*a*dt*v1/Ra;
	    }
	    Dens(sa) = k2*Dens(s2) + k1*Dens(s1);
	    Energy(sa) = k2*Energy(s2) + k1*Energy(s1) - lmpg*(v2*pr2 - v1*pr1);
	    if (is_rotational_symmetry() &&
		(a > 0.0) && (idir == 0) && (fabs(Ra) > fabs(rmin)))
		Energy(sa) -= 0.25*a*dt*(v2*pr2 + v1*pr1)/Ra;
	    for (i = 0; i < dim; ++i)
	    	Mom(sa)[i] = k2*Mom(s2)[i] + k1*Mom(s1)[i];
	    Mom(sa)[idir] -= lmpg*(pr2 - pr1);

	    reset_gamma(sa);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(s2) == ZND)
	    {
	    	rta = rt2 + 0.5 * reaction_rate(s1);
	    	Prod(sa) = k2*Prod(s2) + k1*Prod(s1) + 0.5*dt*rta;
	    	Prod(sa) = max(Prod(sa),0.0);
	    	Prod(sa) = min(Prod(sa),Dens(sa));
	    }
#endif /* defined(COMBUSTION_CODE) */
	    Set_params(sa,s2);
	    if (is_grav == YES)
	    {
	        Mom(sa)[idir] += 0.5*dt*ga*
	           (Dens(sa) + 0.5*(Dens(s2)+Dens(s1)));
	        Energy(sa) += 0.5*dt*ga*
		    (Mom(sa)[idir] + 0.5*(Mom(s2)[idir]+Mom(s1)[idir]));
                reset_gamma(sa);
	    }
	    va = vel(idir,sa);	pra = pressure(sa);
	}

	/* Compute state sb by applying Lax-Friedrichs to states s2 and s3 with
	   modifications due to artificial viscosity. */

	if (is_rotational_symmetry())
	    Rb = (geom == SPHERICAL) ?
	    	2.0*(rad2*rad2+rad2*rad3+rad3*rad3)/(3.0*(rad2+rad3)) :
	    	0.5*(rad2 + rad3);
	for (i = 0; i < dim; ++i)
	    crdsb[i] = 0.5*(Coords(sten->p[0])[i] + Coords(sten->p[1])[i]);
	gb = 0.5*gravity(crdsb,time)[idir];
	lmpg = 0.5*(lm + visc_corr[1][2]);
	k3 = 0.5*(1.0 - visc_corr[1][1]) - lmpg*v3;
	k2 = 0.5*(1.0 + visc_corr[1][1]) + lmpg*v2;
	if (is_rotational_symmetry() &&
	    (a > 0.0) && (idir == 0) && (fabs(Rb) > fabs(rmin)))
	{
	    k3 -= 0.25*a*dt*v3/Rb;
	    k2 -= 0.25*a*dt*v2/Rb;
	}
	Dens(sb) = k3*Dens(s3) + k2*Dens(s2);
	Energy(sb) = k3*Energy(s3) + k2*Energy(s2) - lmpg*(v3*pr3 - v2*pr2);
	if (is_rotational_symmetry() &&
	    (a > 0.0) && (idir == 0) && (fabs(Rb) > fabs(rmin)))
	    Energy(sb) -= 0.25*a*dt*(v2*pr2 + v3*pr3)/Rb;
	for (i = 0; i < dim; ++i)
	    Mom(sb)[i] = k3*Mom(s3)[i] + k2*Mom(s2)[i];
	Mom(sb)[idir] -= lmpg*(pr3 - pr2);

	reset_gamma(sb);
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND)
	{
	    rtb = rt2 + 0.5 * reaction_rate(s3);
	    Prod(sb) = k2*Prod(s2) + k3*Prod(s3) + 0.5*dt*rtb;
	    Prod(sb) = max(Prod(sb),0.0);
	    Prod(sa) = min(Prod(sb),Dens(sb));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(sb,s2);

	if (is_grav == YES)
	{
	    Mom(sb)[idir] += 0.5*dt*gb*
		(Dens(sb) + 0.5*(Dens(s3)+Dens(s2)));
	    Energy(sb) += 0.5*dt*gb*
	        (Mom(sb)[idir] + 0.5*(Mom(s3)[idir]+Mom(s2)[idir]));
            reset_gamma(sb);
	}
	vb = vel(idir,sb);	prb = pressure(sb);

	/* Compute final answer */

	k1 = 0.5*lm*visc_corr[0][0];
	k3 = 0.5*lm*visc_corr[1][0];
	k2 = 1.0 - k1 - k3;
	Dens(ans) = k1*Dens(s1) + k2*Dens(s2) + k3*Dens(s3) -
	    	lm*(Mom(sb)[idir] - Mom(sa)[idir]);
	Energy(ans) = k1*Energy(s1) + k2*Energy(s2) + k3*Energy(s3) -
	    	lm*(vb*(Energy(sb) + prb) - va*(Energy(sa) + pra));
	for (i = 0; i < dim; ++i)
	{
	    Mom(ans)[i] = k1*Mom(s1)[i] + k2*Mom(s2)[i] + k3*Mom(s3)[i] -
			  lm*(vb*Mom(sb)[i] - va*Mom(sa)[i]);
	}
	Mom(ans)[idir] -= lm*(prb - pra);
	if (is_rotational_symmetry() &&
	    (a > 0.0) && (idir == 0))
	{
	    ka = (fabs(Ra) > fabs(rmin)) ? -0.5*a*dt*va/Ra : 0.0;
	    kb = (fabs(Rb) > fabs(rmin)) ? -0.5*a*dt*vb/Rb : 0.0;
	    Dens(ans) += ka*Dens(sa) + kb*Dens(sb);
	    Energy(ans) += ka*(Energy(sa) + pra) + kb*(Energy(sb) + prb);
	    for(i = 0; i < dim; ++i)
		Mom(ans)[i] += ka*Mom(sa)[i] + kb*Mom(sb)[i];
	}
	reset_gamma(ans);
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND)
	{
	    rta = 0.5 * reaction_rate(sa);
	    rtb = 0.5 * reaction_rate(sb);
	    Prod(ans) = k1*Prod(s1) + k2*Prod(s2) + k3*Prod(s3) -
			lm*(vb*Prod(sb) - va*Prod(sa)) + dt*(rtb + rta);
	    Prod(ans) = max(Prod(ans),0.0);
	    Prod(ans) = min(Prod(ans),Dens(ans));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(ans,s2);

	if (is_grav == YES)
	{
	    Mom(ans)[idir] += dt*g2*(Dens(s2) + Dens(ans));
	    Energy(ans) += dt*g2*(Mom(s2)[idir] + Mom(ans)[idir]);
            reset_gamma(ans);
	}

	if (use_linear_artificial_viscosity(Params(s2)->avisc))
	{
	    double coef;
	    double visc = Params(s2)->avisc.linear_visc_coef;

	    coef = visc * lm * mid_max_sp;
	    Dens(ans) += coef*(Dens(s3) - 2.0*Dens(s2) + Dens(s1));
	    for (i = 0; i < dim; ++i)
	        Mom(ans)[i] += coef*(Mom(s3)[i] - 2.0*Mom(s2)[i] + Mom(s1)[i]);
	    Energy(ans) += coef*(Energy(s3) - 2.0*Energy(s2) + Energy(s1));
	    reset_gamma(ans);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(s2) == ZND)
	    	Prod(ans) += coef*(Prod(s3) - 2.0*Prod(s2) + Prod(s1));
#endif /* defined(COMBUSTION_CODE) */
	}

#if defined(CHECK_FOR_BAD_STATES)
	{
	    static const char *fname[3] = {"LWx()","LWy()","LWz()"};
	    if (!check_ans(fname[idir],dh,dt,ans,comp,sten,NO))
	        godunov(dh,dt,ans,dir,swp_num,iperm,index,sten);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
}		/*end LW*/



/*
*				LWoblique():
*
*	Lax Wendroff finite difference scheme in one space variable
*	applied to an arbitrary direction of the conservation form of the
*	equations of compressible flow.
*/

EXPORT void LWoblique(
	double		ds,
	double		dt,
	Tan_stencil	*sten,
	Locstate	ans,
	Front		*fr)
{
	COMPONENT	comp = sten->comp;
	Locstate        *sts = sten->states;
	const double	*dir = sten->dir;
	Locstate	s1 = sts[-1], s2 = sts[0], s3 = sts[1];
	double		lms, lmspg, k1, k2, k3, mid_max_sp;
	double		dir1[MAXD], dir3[MAXD];
	double		v1, pr1, v2, pr2, v3, pr3, va, pra, vb, prb;
	double		g2, ga, gb;
	double		time = fr->time + 0.5*dt;
	RECT_GRID	*gr = fr->rect_grid;
	int		i, dim = gr->dim;
#if defined(COMBUSTION_CODE)
	double		coef, rt2, rta, rtb;
	boolean		reaction;
#endif /* defined(COMBUSTION_CODE) */
	static size_t	sizest;
	static boolean	is_grav;
	static Locstate sa = NULL,sb = NULL;
	static double	**visc_corr = NULL;
	POINT           *p1 = sten->p[-1], *p2 = sten->p[0], *p3 = sten->p[1];
	static double    a;
	static GEOMETRY geom;
	double           rad1, rad2, rad3, Ra, Rb, rmin = pos_radius(0.0,gr);
	double           ka, kb;

	if (is_obstacle_state(s2))
	{
	    g_obstacle_state(ans,fr->sizest);
	    return;
	}
	if (sb == NULL)
	{
	    sizest = Params(s2)->sizest;
	    (*Params(s2)->_alloc_state)(&sa,sizest);
	    (*Params(s2)->_alloc_state)(&sb,sizest);
	    set_type_of_state(sa,GAS_STATE);
	    set_type_of_state(sb,GAS_STATE);
	    bi_array(&visc_corr,2,3,FLOAT);
	    is_grav = is_gravity();
	    if (is_rotational_symmetry())
		geom = Geometry(&a);
	}

	if (RegionIsFlowSpecified(ans,s2,Coords(sten->p[0]),comp,comp,fr))
	    return;
	set_type_of_state(ans,GAS_STATE);


#if defined(COMBUSTION_CODE)
	if (sten->newhs != NULL && Composition_type(s2) == ZND)
	{
	    reaction = (wave_type(sten->newhs) < FIRST_PHYSICS_WAVE_TYPE)
	        ? YES : NO;
	    coef = 0.03;
	    rt2 = reaction ? coef * reaction_rate(s2) : 0.0;
	}
#endif /* defined(COMBUSTION_CODE) */

#define OLD

#if defined(OLD)
	for (i = 0; i < dim; ++i)
	{
	    dir1[i] = dir3[i] = dir[i];
	}
#else /* defined(OLD) */
	/* if !defined(OLD) */
/*******
*	The code in #if defined was put it because of our analysis
*	of the correct way to split in "Front Tracking for Gas Dynamics".
*	The corrections would make the splitting second order if we could
*	get second order states.
******/


	if (sten->newhs == NULL) /* points p1,p2,p3 aren't supplied */
	{
	    for (i = 0; i < dim; ++i)
	    {
	    	dir1[i] = dir3[i] = dir[i];
	    }
	}
	else
	{
	    POINT *p1 = sten->p[-1], *p2 = sten->p[0], *p3 = sten->p[1];
	    double len;

	    len = 0.0;
	    for ( i = 0; i < dim; ++i)
	    {
	    	dir1[i] = Coords(p2)[i] - Coords(p1)[i];
	    	len += sqr(dir1[i]);
	    }
	    len = sqrt(len);
	    if (len < 0.01*ds)
	    {
	    	for (i = 0; i < dim; ++i)
		    dir1[i] = dir[i];
	    }
	    else
	    {
	    	for ( i = 0; i < dim; ++i)
		    dir1[i] /= len;
	    }

	    len = 0.0;
	    for ( i = 0; i < dim; ++i)
	    {
	    	dir3[i] = Coords(p3)[i] - Coords(p2)[i];
	    	len += sqr(dir3[i]);
	    }
	    len = sqrt(len);
	    if (len < 0.01*ds)
	    {
	    	for (i = 0; i < dim; ++i)
		    dir3[i] = dir[i];
	    }
	    else
	    {
	    	for ( i = 0; i < dim; ++i)
		    dir3[i] /= len;
	    }
	}
#endif /* defined(OLD) */

	v1 = v2 = v3 = 0.;
	for (i = 0; i < dim; ++i)
	{
	    v1 += vel(i,s1)*dir1[i];
	    v2 += vel(i,s2)*dir[i];
	    v3 += vel(i,s3)*dir3[i];
	}
	pr1 = pressure(s1);
	pr2 = pressure(s2);
	pr3 = pressure(s3);
	if (is_rotational_symmetry())
	{
	    rad1 = pos_radius(Coords(p1)[0],gr);
	    rad2 = pos_radius(Coords(p2)[0],gr);
	    rad3 = pos_radius(Coords(p3)[0],gr);
	}
	mid_max_sp = compute_art_visc_corrections(s1,v1,dir1,
	                                          Coords(sten->p[-1]),
	                                          s2,v2,dir,Coords(sten->p[0]),
	                                          s3,v3,dir3,Coords(sten->p[1]),
	                                          visc_corr,TSWEEP,NULL,fr);
	for (i = 0; i < dim; ++i)
	{
	    set_max_front_speed(i,fabs(dir[i]*mid_max_sp),s2,
				Coords(sten->p[0]),fr);
	}
	set_max_front_speed(dim,mid_max_sp/ds,s2,Coords(sten->p[0]),fr);

	if (is_rotational_symmetry())
	    Ra = (geom == SPHERICAL) ?
	    	2.0*(rad1*rad1+rad1*rad2+rad2*rad2)/(3.0*(rad1+rad2)) :
	    	0.5*(rad1 + rad2);
	lms = dt/ds;
	lmspg = 0.5*(lms + visc_corr[0][2]);
	k2 = 0.5*(1.0 - visc_corr[0][1]) - lmspg*v2;
	k1 = 0.5*(1.0 + visc_corr[0][1]) + lmspg*v1;
	if (is_rotational_symmetry() && (a > 0.0) && (Ra > rmin))
	{
	    k2 -= 0.25*a*dir[0]*dt*v2/Ra;
	    k1 -= 0.25*a*dir[0]*dt*v1/Ra;
	}
	Dens(sa) = k2*Dens(s2) + k1*Dens(s1);
	Energy(sa) = k2*Energy(s2) + k1*Energy(s1) - lmspg*(v2*pr2 - v1*pr1);
	if (is_rotational_symmetry() && (a > 0.0) && (Ra > rmin))
	    Energy(sa) -= 0.25*a*dir[0]*dt*(v2*pr2 + v1*pr1)/Ra;
	for (i = 0; i < dim; ++i)
	    Mom(sa)[i] = k2*Mom(s2)[i]+k1*Mom(s1)[i]-lmspg*dir[i]*(pr2-pr1);

	reset_gamma(sa);
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND)
	{
	    rta = reaction ? (rt2 + coef * reaction_rate(s1)) : 0.0;
	    Prod(sa) = k2*Prod(s2) + k1*Prod(s1) + 0.5*dt*rta;
	    Prod(sa) = max(Prod(sa),0.0);
	    Prod(sa) = min(Prod(sa),Dens(sa));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(sa,s2);

	if (is_rotational_symmetry())
	    Rb = (geom == SPHERICAL) ?
	    	2.0*(rad2*rad2+rad2*rad3+rad3*rad3)/(3.0*(rad2+rad3)) :
	    	0.5*(rad2 + rad3);
	lmspg = 0.5*(lms + visc_corr[1][2]);
	k3 = 0.5*(1.0 - visc_corr[1][1]) - lmspg*v3;
	k2 = 0.5*(1.0 + visc_corr[1][1]) + lmspg*v2;
	if (is_rotational_symmetry() && (a > 0.0) && (Rb > rmin))
	{
	    k3 -= 0.25*a*dir[0]*dt*v3/Rb;
	    k2 -= 0.25*a*dir[0]*dt*v2/Rb;
	}
	Dens(sb) = k3*Dens(s3) + k2*Dens(s2);
	Energy(sb) = k3*Energy(s3) + k2*Energy(s2) - lmspg*(v3*pr3 - v2*pr2);
	if (is_rotational_symmetry() && (a > 0.0) && (Rb > rmin))
	    Energy(sb) -= 0.25*a*dir[0]*dt*(v2*pr2 + v3*pr3)/Rb;
	for (i = 0; i < dim; ++i)
	    Mom(sb)[i] = k3*Mom(s3)[i]+k2*Mom(s2)[i]-lmspg*dir[i]*(pr3-pr2);

	reset_gamma(sb);
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND)
	{
	    rtb = reaction ? (rt2 + coef * reaction_rate(s3)) : 0.0;
	    Prod(sb) = k3*Prod(s3) + k2*Prod(s2) + 0.5*dt*rtb;
	    Prod(sb) = max(Prod(sb),0.0);
	    Prod(sa) = min(Prod(sb),Dens(sb));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(sb,s2);

	if (is_grav == YES)
	{
	    double		crdsa[3], crdsb[3];
	    for (i = 0; i < dim; ++i)
	    {
		crdsa[i] = 0.5*(Coords(sten->p[-1])[i] + Coords(sten->p[0])[i]);
		crdsb[i] = 0.5*(Coords(sten->p[0])[i]  + Coords(sten->p[1])[i]);
	    }
	    ga = 0.5*scalar_product(gravity(crdsa,time),dir,dim);
	    gb = 0.5*scalar_product(gravity(crdsb,time),dir,dim);

	    for (i = 0; i < dim; ++i)
	    {
	        Mom(sa)[i] += 0.5*dt*ga*dir[i]*
		              (Dens(sa) + 0.5*(Dens(s2)+Dens(s1)));
	        Mom(sb)[i] += 0.5*dt*gb*dir[i]*
		              (Dens(sb) + 0.5*(Dens(s3)+Dens(s2)));
	    }
	}
	else
	{
	    ga = gb = 0.0;
	}

	va = vb = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    va += vel(i,sa)*dir[i];
	    vb += vel(i,sb)*dir[i];
	}

	if (is_grav == YES)
	{
	    Energy(sa) += 0.5*dt*ga*(va*Dens(sa)+0.5*(v2*Dens(s2)+v1*Dens(s1)));
	    Energy(sb) += 0.5*dt*gb*(vb*Dens(sb)+0.5*(v3*Dens(s3)+v2*Dens(s2)));
            reset_gamma(sa);
            reset_gamma(sb);
	}

	pra = pressure(sa);
	prb = pressure(sb);

	k1 = 0.5*lms*visc_corr[0][0];
	k3 = 0.5*lms*visc_corr[1][0];
	k2 = 1.0 - k1 - k3;
	Dens(ans) = k1*Dens(s1) + k2*Dens(s2) + k3*Dens(s3) -
		    lms*(vb*Dens(sb) - va*Dens(sa));
	Energy(ans) = k1*Energy(s1) + k2*Energy(s2) + k3*Energy(s3) -
		      lms*(vb*(Energy(sb) + prb) - va*(Energy(sa) + pra));
	for (i = 0; i < dim; ++i)
	{
	    Mom(ans)[i] = k1*Mom(s1)[i] + k2*Mom(s2)[i] + k3*Mom(s3)[i] -
	                  lms*(vb*Mom(sb)[i]-va*Mom(sa)[i]+dir[i]*(prb-pra));
	}
	if (is_rotational_symmetry() && a > 0.0)
	{
	    ka = (Ra > rmin) ? -0.5*a*dir[0]*dt*va/Ra : 0.0;
	    kb = (Rb > rmin) ? -0.5*a*dir[0]*dt*vb/Rb : 0.0;
	    Dens(ans) += ka*Dens(sa) + kb*Dens(sb);
	    Energy(ans) += ka*(Energy(sa) + pra) + kb*(Energy(sb) + prb);
	    for (i = 0; i < dim; ++i)
		Mom(ans)[i] += ka*Mom(sa)[i] + kb*Mom(sb)[i];
	}

	reset_gamma(ans);
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND)
	{
	    rta = reaction ? coef * reaction_rate(sa) : 0.0;
	    rtb = reaction ? coef * reaction_rate(sb) : 0.0;
	    Prod(ans) = k1*Prod(s1) + k2*Prod(s2) + k3*Prod(s3) -
	    	        lms*(vb*Prod(sb) - va*Prod(sa)) + dt*(rtb + rta);
	    Prod(ans) = max(Prod(ans),0.0);
	    Prod(ans) = min(Prod(ans),Dens(ans));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(ans,s2);

	if (is_grav == YES)
	{
	    double vans;

	    g2 = 0.5*scalar_product(gravity(Coords(sten->p[0]),time),dir,dim);
	    for (i = 0; i < dim; ++i)
	        Mom(ans)[i] += dt*g2*dir[i]*(Dens(s2) + Dens(ans));
	    vans = 0.0;
	    for (i = 0; i < dim; ++i)
		vans += vel(i,ans)*dir[i];
	    Energy(ans) += dt*g2*(v2*Dens(s2) + vans*Dens(ans));
            reset_gamma(ans);
	}

	if (use_linear_artificial_viscosity(Params(s2)->avisc))
	{
	    double coef;
	    double visc = Params(s2)->avisc.linear_visc_coef;

	    coef = visc * lms * mid_max_sp;
	    Dens(ans) += coef*(Dens(s3) - 2.0*Dens(s2) + Dens(s1));
	    for (i = 0; i < dim; ++i)
	    	Mom(ans)[i] += coef*(Mom(s3)[i] - 2.0*Mom(s2)[i] + Mom(s1)[i]);
	    Energy(ans) += coef*(Energy(s3) - 2.0*Energy(s2) + Energy(s1));
	    reset_gamma(ans);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(s2) == ZND)
	    	Prod(ans) += coef*(Prod(s3) - 2.0*Prod(s2) + Prod(s1));
#endif /* defined(COMBUSTION_CODE) */
	}

#if defined(CHECK_FOR_BAD_STATES)
	if (!check_gas("LWoblique",sts,ans,sten,NO,fr))
	    godunovobl(ds,dt,sten,ans,fr);
#endif /* defined(CHECK_FOR_BAD_STATES) */
}		/*end LWoblique*/


/*
*			compute_art_visc_corrections():
*
*	Computes the non-linear artificial viscosity terms for
*	the Lax-Wendroff method.  See Richtmyer and Morton
*	"Difference Methods for Initial-Value Problems" page 336
*	for details.  This function also returns the quantity
*	sp_coef*(fabs(v2) + sound_speed(s2)), where
*	sp_coef = 1.0/(sqrt(1.0 + 0.25*sqr(visc)) - 0.5*visc), and
*	visc is the coefficient of artificial viscosity.
*	This quantity is used in computing the next time step dt.
*
*	IMPORTANT NOTE:
*	If both the stencil and its previous stencil are regular,
*	then s0 equates the previous stencil's s1 and s1 equates
*	the previous stencils s2. Thus the coefficents corresponding
*	to the prev portion (s0, s1) of the current stencil are equal
*	to those computed for the next portion (s1, s2) of the previous
*	stencil.  This fact allows for a substantial computational savings
*	in regions away from any tracked fronts.
*/

/*ARGSUSED*/
LOCAL	double compute_art_visc_corrections(
	Locstate    s0,
	double       v0,
	const double *dir0,
	double       *crds0,
	Locstate    s1, double v1,
	const double *dir1,
	double *crds1,
	Locstate    s2,
	double       v2,
	const double *dir2,
	double       *crds2,
	double	    **g,
	int	    sweep,
	Stencil	    *sten,
	Front	    *fr)
{
	int	        dim = fr->rect_grid->dim;
	double		uh[2], chs[2];
	double		u[3];
	double		b[3];
	double		chi[3][3];
	double		lambda0, ch;
	double		visc, sp_coef;
	int		i, j;
	static const double VAC_FAC = 0.00025; /* TOLERANCE */
	static double	c[3];
	static Locstate smh = NULL, sph = NULL;

	if (sph == NULL)
	{
	    size_t	sizest = Params(s2)->sizest;
	    (*Params(s2)->_alloc_state)(&smh,sizest);
	    (*Params(s2)->_alloc_state)(&sph,sizest);
	}

	if (!use_lapidus_artificial_viscosity(Params(s1)->avisc))
	{
	    for (i = 0; i < 3; ++i)
	    	g[0][i] = g[1][i] = 0.0;
	    return fabs(v1) + sound_speed(s1);
	}

	visc = Params(s2)->avisc.lapidus_visc_coef;
	sp_coef = Params(s2)->avisc.sp_coef;
	if (sten && sten->reg_stencil && sten->prev_reg_stencil)
	{
	    for (i = 0; i < 3; ++i) g[0][i] = g[1][i];

	/* Compute only the s1-s2 artificial viscosity coefficients */

	    interpolate_states(fr,0.5,0.5,crds1,s1,crds2,s2,sph);
	    u[1] = v1;	 u[2] = v2;
	    if (sweep == TSWEEP)
	    {
	    	uh[1] = 0.0;
	    	for (i = 0; i < dim; ++i)
	    	    uh[1] += dir2[i]*vel(i,sph);
	    }
	    else
	    	uh[1] = vel(sweep,sph);
	    c[1] = c[2];
	    c[2] = sound_speed(s2);
	    for (i = 1; i < 3; ++i)
	    {
	    	chi[i][0] = u[i] - c[i];
	    	chi[i][1] = u[i];
	    	chi[i][2] = u[i] + c[i];
	    }
	/* TODO  Sometimes , chs[i] ~ 0 causes trouble */
	    chs[1] = sound_speed_squared(sph);
	    ch = sqrt(chs[1]);
	    if (ch <= VAC_FAC*(fabs(u[1])+c[1]+fabs(u[2])+c[2]))
	    {
	    	g[1][0] = visc*fabs(u[2]-u[1]);
	    	g[1][1] = 0.0;
	    	g[1][2] = 0.0;
	    }
	    else
	    {
	    	for (j = 0; j < 3; ++j)
	    	    b[j] = visc*fabs(chi[2][j] - chi[1][j]);
	        g[1][2] = (b[0] - 2.0*b[1] + b[2])/(2.0*chs[1]);
		lambda0 = uh[1] - ch;
	        g[1][1] = (b[1] - b[0])/ch - (uh[1] + lambda0)*g[1][2];
	        g[1][0] = b[0] - g[1][1]*lambda0 - g[1][2]*sqr(lambda0);
	    }
	}
	else
	{
	    interpolate_states(fr,0.5,0.5,crds0,s0,crds1,s1,smh);
	    interpolate_states(fr,0.5,0.5,crds1,s1,crds2,s2,sph);
	    u[0] = v0;
	    u[1] = v1;
	    u[2] = v2;
	    if (sweep == TSWEEP)
	    {
	    	uh[0] = 0.0;	uh[1] = 0.0;
	    	for (i = 0; i < dim; ++i)
	    	{
	            uh[0] += dir0[i]*vel(i,smh);
	            uh[1] += dir2[i]*vel(i,sph);
	        }
	    }
	    else
	    {
	    	uh[0] = vel(sweep,smh);
	    	uh[1] = vel(sweep,sph);
	    }
	    c[0] = sound_speed(s0);
	    c[1] = sound_speed(s1);
	    c[2] = sound_speed(s2);
	    for (i = 0; i < 3; ++i)
	    {
	    	chi[i][0] = u[i] - c[i];
	    	chi[i][1] = u[i];
	    	chi[i][2] = u[i] + c[i];
	    }
	    chs[0] = sound_speed_squared(smh);
	    chs[1] = sound_speed_squared(sph);
	    for (i = 0; i < 2; ++i)
	    {
	    	ch = sqrt(chs[i]);
	    	if (ch <= VAC_FAC*(fabs(u[i])+c[i]+fabs(u[i+1])+c[i+1]))
	    	{
	            g[i][0] = visc*fabs(u[i+1]-u[i]);
	            g[i][1] = 0.0;
	            g[i][2] = 0.0;
	        }
	        else
	        {
	            for (j = 0; j < 3; ++j)
	                b[j] = visc*fabs(chi[i+1][j] - chi[i][j]);
	            g[i][2] = (b[0] - 2.0*b[1] + b[2])/(2.0*chs[i]);
	            lambda0 = uh[i] - ch;
	            g[i][1] = (b[1] - b[0])/ch - (uh[i] + lambda0)*g[i][2];
	            g[i][0] = b[0] - g[i][1]*lambda0 - g[i][2]*sqr(lambda0);
	        }
	    }
	}
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("VIS"))
	{
	    int bad_coef = NO;
	    for (i = 0; i < 2; ++i)
	    {
	        for (j = 0; j < 3; ++j)
	        {
	            if (fabs(g[i][j]) >= 1000.0)
	            {
	                (void) printf("WARNING: g[%d][%d] = %g\n",i,j,g[i][j]);
	                bad_coef = YES;
	            }
	        }
	    }
	    if (bad_coef)
	    {
	    	(void) printf("u[0] = %g, u[1] = %g, u[2] = %g\n",
	    	              u[0],u[1],u[2]);
	    	(void) printf("uh[0] = %g, uh[1] = %g\n",uh[0],uh[1]);
	    	(void) printf("c[0] = %g, c[1] = %g, c[2] = %g\n",
	    	              c[0],c[1],c[2]);
	    	(void) printf("chs[0] = %g, chs[1] = %g\n",chs[0],chs[1]);
	    	verbose_print_state("s0",s0);
	    	verbose_print_state("s1",s1);
	    	verbose_print_state("s2",s2);
	    }
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */

	return sp_coef*(fabs(v1) + c[1]);
}		/*end compute_art_visc_corrections*/



LOCAL	double *visc[3] = {NULL, NULL, NULL};

EXPORT void LW_alloc_phys_vecs(
	Wave		*wave,
	int		vctr_size)
{
	Vec_Gas		*vst;
	int		i;
	int		dim = wave->rect_grid->dim;

	scalar(&vst,3*sizeof(Vec_Gas));
	g_wave_vgas(wave) = vst;
	g_wave_vsrc(wave) = NULL;
	for (i = 0; i < 3; ++i)
	{
	    (void) g_alloc_vgas(vst+i,vctr_size,dim);
	    vst[i].alloc.vst = NO;
	}
}		/*end LW_alloc_phys_vecs*/

/*ARGSUSED*/
EXPORT void LW_free_phys_vecs(
	Wave		*wave)
{
	Vec_Gas		*vst;
	int		i;

	vst = g_wave_vgas(wave);
	g_wave_vgas(wave) = NULL;
	if (vst != NULL)
	{
	    for (i = 0; i < 3; ++i)
	        g_free_vgas(vst+i);
	    free(vst);
	}
}		/*end LW_free_phys_vecs*/


/*
*			  LW_vec()
*/

/*ARGSUSED*/
EXPORT	void LW_vec(
	int		swp_num,
	int		*iperm,
	double		*dir,
	Wave		*wv,
	Wave		*newwv,
	Front		*fr,
	Front		*newfr,
	int		*icoords,
	int		imin,
	int		imax,
	double		dt,
	double		dh)
{
	int		vsize;
	Vec_Gas		*vst = g_wave_vgas(wv);
	int		nrad = vsten_radius(wv);

#if defined(TIME_HYPVEC)
	start_clock("LW_vec");
#endif /* defined(TIME_HYPVEC) */
	debug_print("hyp","Entered LW_vec()\n");
	vsize = imax - imin;
	if (load_state_vectors(swp_num,iperm,&vst[0],0,vsize,wv,newwv,icoords,
			       imin) == CONSTANT_IN_TIME)
	{
#if defined(TIME_HYPVEC)
	    stop_clock("LW_vec");
#endif /* defined(TIME_HYPVEC) */
	    return;
	}
	vec_Lax_Wendroff(iperm[swp_num],dh,dt,vsize,&vst[0],&vst[1],&vst[2],
			 visc,wv);
	assign_wave_state_vectors(swp_num,iperm,wv,newwv,&vst[2],nrad,
				  vsize-nrad,icoords,imin);
	debug_print("hyp","Left LW_vec()\n");
#if defined(TIME_HYPVEC)
	stop_clock("LW_vec");
#endif /* defined(TIME_HYPVEC) */
}		/*end LW_vec*/


LOCAL void vec_Lax_Wendroff(
	int		idir,
	double		dh,
	double		dt,
	int		vs,
	Vec_Gas		*vst0,
	Vec_Gas		*vst1,
	Vec_Gas		*ans,
	double		**visc,
	Wave		*wave)
{
	int		vsm1 = vs - 1, vsm2 = vs - 2, i, j;
	int		dim = wave->rect_grid->dim;
	double		lmx = dt/dh;
	double		*lmxpvsc2d2, *k1, *k2, *k3;
	double		*g_e_corr, *g_mx_corr;
	double		*rho0 = vst0->rho, *rho1 = vst1->rho;
	double		*rhoans = ans->rho + 1;
	double		*e0 = vst0->en_den, *e1 = vst1->en_den;
	double		*eans = ans->en_den + 1;
	double		*mx0 = vst0->m[0], *mx1 = vst1->m[0];
	double		*mxans = ans->m[0] + 1;
	double		*my0 = vst0->m[1], *my1 = vst1->m[1];
	double		*myans = ans->m[1] + 1;
	double		*v0 = vst0->v[0], *v1 = vst1->v[0];
	double		*p0 = vst0->p, *p1 = vst1->p;
	double		*visc_coef, *sp_coef;
	double		**coords = vst0->coords;
	double		time = wave->time + 0.5*dt;
	static	double	*g_dt = NULL;
	static  int	g_dt_size;
	static	boolean	first = YES;
	static	boolean	is_grav;

	if (first == YES)
	{
	    first = NO;
	    is_grav = is_gravity();
	}


	if (g_dt == NULL)
	{
	    g_dt_size = vs;
	    uni_array(&g_dt,vs,FLOAT);
	}
	else if (g_dt_size < vs)
	{
	    free(g_dt);
	    g_dt_size = vs;
	    uni_array(&g_dt,vs,FLOAT);
	}


	/* Load states from old interface */

	for (i = 0; i < vs; ++i)
	    v0[i] = mx0[i]/rho0[i];
	load_pressure_and_sound_speed(vst0,0,vs);

	/* Compute artificial viscosity coefficients */

	for (i = 0; i < vsm1; ++i)
	{
	    rho1[i] = 0.5*(rho0[i] + rho0[i+1]);
	    e1[i] = 0.5*(e0[i] + e0[i+1]);
	    mx1[i] = 0.5*(mx0[i] + mx0[i+1]);
	    my1[i] = 0.5*(my0[i] + my0[i+1]);
	}

	if (is_grav == YES)
	{
	    double    crds[3];
	    for (i = 0; i < vsm1; ++i)
	    {
	        for (j = 0; j < dim; ++j)
		    crds[j] = 0.5*(coords[i][j] + coords[i+1][j]);
	        g_dt[i] = 0.5*gravity(crds,time)[idir]*dt;
	    }
	}

	for (i = 0; i < vsm1; ++i)
	    v1[i] = mx1[i]/rho1[i];
	visc_coef = ans->p; /*Work space*/
	sp_coef = ans->c;/*Work space*/
	for (i = 0; i < vsm1; ++i)
	{
	    vst1->state[i] = vst0->state[i];
	    visc_coef[i] = Params(vst0->state[i])->avisc.lapidus_visc_coef;
	    sp_coef[i] = Params(vst0->state[i])->avisc.sp_coef;
	}
	sp_coef[vsm1] = (Params(vst0->state[vsm1]))->avisc.sp_coef;
	load_pressure_and_sound_speed(vst1,0,vsm1);
	for (i = 0; i < vs; ++i)
	    set_max_wave_speed(idir,fabs(v0[i])+vst0->c[i],
	    		       vst0->state[i],vst0->coords[i],wave);
	FORTRAN_NAME(artvsc)(v0,vst0->c,v1,vst1->c,visc[0],visc[1],visc[2],
		             ans->rho,ans->m[0],ans->m[1],ans->en_den,
		             visc_coef,sp_coef,&vs);

	/* Compute states at time t + 0.5*dt */

			/* Use ans as working storage */

	lmxpvsc2d2 = ans->rho;
	k1 = ans->m[0];
	k2 = ans->m[1];
	g_e_corr = ans->v[0];
	g_mx_corr = ans->p;

	FORTRAN_NAME(lwvec1)(rho0,e0,mx0,my0,p0,v0,rho1,e1,mx1,my1,v1,
		             g_e_corr,g_mx_corr,k1,k2,visc[1],visc[2],
			     lmxpvsc2d2,&lmx,g_dt,&is_grav,&vsm1);

	Vec_Gas_field_set(vst1,rho) = YES;
	Vec_Gas_field_set(vst1,en_den) = YES;
	Vec_Gas_field_set(vst1,m) = YES;
	load_pressure(vst1,0,vsm1);


	/* Compute final answer */

		/* Use ans->v[0], ans->p, and ans->c as work space*/

	k1 = ans->v[0]; k2 = ans->p; k3 = ans->c;

	if (is_grav == YES)
	{
	    for (i = 0; i < vsm2; ++i)
	    {
	        g_dt[i] = 0.5*gravity(coords[i+1],time)[idir]*dt;
	    }
	}

	FORTRAN_NAME(lwvec2)(rho0,e0,mx0,my0,rho1,e1,mx1,my1,
		             p1,v1,rhoans,eans,mxans,myans,k1,k2,k3,visc[0],
		             &vsm2,&lmx,g_dt,&is_grav);
	Vec_Gas_field_set(ans,rho) = YES;
	Vec_Gas_field_set(ans,en_den) = YES;
	Vec_Gas_field_set(ans,m) = YES;
	for (i = 1; i < vsm1; ++i)
	    ans->state[i] = vst0->state[i];
}		/*end vec_Lax_Wendroff*/


