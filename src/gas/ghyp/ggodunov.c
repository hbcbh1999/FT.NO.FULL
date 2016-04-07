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
*				ggodunov.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains Godunov routines used in the hyperbolic step.
*
*/

#include <ghyp/ghyp.h>

#define HYP_DBG

/*
*				godunov():
*
*      Godunov finite difference scheme in one space variable
*      applied to the t-x part of the conservation form of the
*      equations of compressible flow.
*
*/

/*ARGSUSED*/
EXPORT void godunov(
	double		dh,
	double		dt,
	Locstate	ans,
	const double	*dir,
	int             swp_num,
	int             *iperm,
	int		*index,
	Stencil		*sten)
{
	COMPONENT	comp = sten->newcomp;
	Locstate	s1 = sten->st[-1], s2 = sten->st[0], s3 = sten->st[1];
	int		i, dim = sten->fr->interf->dim;
	int             idir = iperm[swp_num];
	double		lm, speed;
	double		va, pra, vb, prb;
	double		g2, ga, gb;
	double		time = sten->fr->time + 0.5*dt;
#if defined(COMBUSTION_CODE)
	double		rta, rtb;
#endif /* defined(COMBUSTION_CODE) */
#if defined(HYP_DBG)
	static	char	fname[3][12] = {"godunovx()","godunovy()","godunovy()"};
#endif /* defined(HYP_DBG) */
	static size_t	sizest;
	static boolean	is_grav;
	static Locstate sa = NULL, sb = NULL;
	static double    a;
	static GEOMETRY geom;
	double           rad1, rad2, rad3, Ra, Rb, rmin;
	double           ka, kb;
	double           v1, pr1, v2, pr2, v3, pr3;

	debug_print("godunov","Entered godunov(), idir = %d\n",idir);

	if (is_obstacle_state(s2)) 
	{
	    g_obstacle_state(ans,sten->fr->sizest);
	    debug_print("godunov","Left godunov()\n");
	    return;
	}
	if (sa == NULL) 
	{
	    sizest = Params(s2)->sizest;
	    (*Params(s2)->_alloc_state)(&sa,sizest);
	    (*Params(s2)->_alloc_state)(&sb,sizest);
	    set_type_of_state(sa,GAS_STATE);
	    set_type_of_state(sb,GAS_STATE);
	    is_grav = is_gravity();
	    if (is_rotational_symmetry())
		geom = Geometry(&a);
	}
	if (debugging("godunov"))
	{
	    (void) printf("Input data into godunovobl()\n");
	    (void) printf("dh = %g, dt = %g, comp = %d\n",dh,dt,comp);
	    for (i = 0; i < dim; ++i)
	    {
	        (void) printf("dir[%d] = %g%s",i,dir[i],
	                      (i == dim-1) ? "\n" : ", ");
	    }
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	}
	set_type_of_state(ans,GAS_STATE);

	if (is_grav == YES)
	{
	    double    crdsa[3], crdsb[3];
	    for (i = 0; i < dim; ++i)
	    {
		crdsa[i] = 0.5*(Coords(sten->p[-1])[i] + Coords(sten->p[0])[i]);
		crdsb[i] = 0.5*(Coords(sten->p[0])[i] + Coords(sten->p[1])[i]);
	    }
	    ga = 0.5*gravity(crdsa,time)[idir];
	    gb = 0.5*gravity(crdsb,time)[idir];
	    g2 = 0.5*gravity(Coords(sten->p[0]),time)[idir];
	}
	else
	{
	    ga = gb = g2 = 0.0;
	}
	if (is_rotational_symmetry() && (a > 0.0) && (idir == 0))
	{
	    RECT_GRID *gr = sten->fr->rect_grid;
	    v1 = vel(idir,s1);
	    pr1 = pressure(s1);
	    v2 = vel(idir,s2);
	    pr2 = pressure(s2);
	    v3 = vel(idir,s3);
	    pr3 = pressure(s3);
	    rad1 = pos_radius(Coords(sten->p[-1])[0],gr);
	    rad2 = pos_radius(Coords(sten->p[0])[0],gr);
	    rad3 = pos_radius(Coords(sten->p[1])[0],gr);
	    rmin = fabs(pos_radius(0.0,gr));
	}

	lm = dt/dh;

	/* Find sa */

	riemann_solution(0.0,dir,s1,s2,sa,GAS_STATE);

	if (is_rotational_symmetry() && (a > 0.0) && (idir == 0))
	{
	    Ra = (geom == SPHERICAL) ?
		 2.0*(rad1*rad1+rad1*rad2+rad2*rad2)/(3.0*(rad1+rad2)) :
		 0.5*(rad1 + rad2);
	    if (fabs(Ra) > rmin)
	    {
		Dens(sa) += - 0.25*a*dt*(Dens(s1)*v1 + Dens(s2)*v2)/Ra;
		Energy(sa) += - 0.25*a*dt*(Energy(s1)*v1 + Energy(s2)*v2)/Ra
			      - 0.25*a*dt*(pr1*v1 + pr2*v2)/Ra;
		for(i = 0; i < dim; ++i)
		    Mom(sa)[i] += -0.25*a*dt*(Mom(s1)[i]*v1 + Mom(s2)[i]*v2)/Ra;
		reset_gamma(sa);
	    }
	}

	if (is_grav == YES) 
	{
	    Mom(sa)[idir] += 0.5*dt*ga*
	                     (Dens(sa) + 0.5*(Dens(s2) + Dens(s1)));
	    Energy(sa) += 0.5*dt*ga*
	                  (Mom(sa)[idir] + 0.5*(Mom(s2)[idir] + Mom(s1)[idir]));
	    reset_gamma(sa);
	}
#if defined(HYP_DBG)
	if (!check_ans(fname[idir],dh,dt,sa,comp,sten,YES))
	{
	    (void) printf("WARNING in godunov(), bad state detected for sa\n");
	    print_general_vector("dir = ",dir,dim," ");
	    (void) printf("idir = %d\n",idir);
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	    (void) printf("sa ");
	    fprint_raw_gas_data(stdout,sa,dim);
	    LF(dh,dt,ans,dir,swp_num,iperm,index,sten);
	    return;
	}
#endif /* defined(HYP_DBG) */

	/* Compute sb */

	riemann_solution(0.0,dir,s2,s3,sb,GAS_STATE);

	if (is_rotational_symmetry() && (a > 0.0) && (idir == 0))
	{
	    Rb = (geom == SPHERICAL) ?
			2.0*(rad2*rad2+rad2*rad3+rad3*rad3)/(3.0*(rad2+rad3)) :
			0.5*(rad2 + rad3);
	    if (fabs(Rb) > rmin)
	    {
		Dens(sb) += - 0.25*a*dt*(Dens(s2)*v2 + Dens(s3)*v3)/Rb;
		Energy(sb) += - 0.25*a*dt*(Energy(s2)*v2 + Energy(s3)*v3)/Rb
			      - 0.25*a*dt*(pr2*v2 + pr3*v3)/Rb;
		for(i = 0; i < dim; ++i)
		    Mom(sb)[i] += -0.25*a*dt*(Mom(s2)[i]*v2 + Mom(s3)[i]*v3)/Rb;
		reset_gamma(sb);
	    }
	}

	if (is_grav == YES) 
	{
	    Mom(sb)[idir] += 0.5*dt*gb*
	                     (Dens(sb) + 0.5*(Dens(s3) + Dens(s2)));
	    Energy(sb) += 0.5*dt*gb*
	                  (Mom(sb)[idir] + 0.5*(Mom(s3)[idir] + Mom(s2)[idir]));
	    reset_gamma(sb);
	}
#if defined(HYP_DBG)
	if (!check_ans(fname[idir],dh,dt,sb,comp,sten,YES))
	{
	    (void) printf("WARNING in godunov(), bad state detected for sb\n");
	    print_general_vector("dir = ",dir,dim," ");
	    (void) printf("idir = %d\n",idir);
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	    (void) printf("sb ");
	    fprint_raw_gas_data(stdout,sb,dim);
	    LF(dh,dt,ans,dir,swp_num,iperm,index,sten);
	    return;
	}
#endif /* defined(HYP_DBG) */

	/* Compute final answer */

	va = vel(idir,sa);	pra = pressure(sa);
	vb = vel(idir,sb);	prb = pressure(sb);

	Dens(ans) = Dens(s2) - lm*(Mom(sb)[idir] - Mom(sa)[idir]);
	Energy(ans) = Energy(s2) -
	              lm*(vb*(Energy(sb) + prb) - va*(Energy(sa) + pra));
	for (i = 0; i < dim; ++i)
	{
	    Mom(ans)[i] = Mom(s2)[i] - lm*(vb*Mom(sb)[i] - va*Mom(sa)[i]);
	}
	Mom(ans)[idir] -= lm*(prb - pra);
	reset_gamma(ans);

	if (is_rotational_symmetry() && (a > 0.0) && (idir == 0))
	{
	    ka = (fabs(Ra) > rmin) ? -0.5*a*dt*va/Ra : 0.0;
	    kb = (fabs(Rb) > rmin) ? -0.5*a*dt*vb/Rb : 0.0;
	    Dens(ans) += ka*Dens(sa) + kb*Dens(sb);
	    Energy(ans) += ka*(Energy(sa) + pra) + kb*(Energy(sb) + prb);
	    for(i = 0; i < dim; ++i)
		Mom(ans)[i] += ka*Mom(sa)[i] + kb*Mom(sb)[i];
	   reset_gamma(ans);
	}

#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND) 
	{
	    rta = 0.5 * reaction_rate(sa);
	    rtb = 0.5 * reaction_rate(sb);
	    Prod(ans) = Prod(s2)-lm*(vb*Prod(sb)-va*Prod(sa)) + dt*(rtb+rta);
	    Prod(ans) = max(Prod(ans),0.0);
	    Prod(ans) = min(Prod(ans),Dens(ans));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(ans,s2);
	if (is_grav == YES) 
	{
	    Mom(ans)[idir] += dt*(Dens(s2) + Dens(ans))*g2;
	    Energy(ans) += dt*(Mom(s2)[idir] + Mom(ans)[idir])*g2;
            reset_gamma(ans);
	}
#if defined(HYP_DBG)
	if (!check_ans(fname[idir],dh,dt,ans,comp,sten,YES))
	{
	    (void) printf("WARNING in godunov(), bad state detected for ans\n");
	    print_general_vector("dir = ",dir,dim," ");
	    (void) printf("idir = %d\n",idir);
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	    (void) printf("ans ");
	    fprint_raw_gas_data(stdout,ans,dim);
	    LF(dh,dt,ans,dir,swp_num,iperm,index,sten);
	}
#endif /* defined(HYP_DBG) */

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Dens(ans) = max(Dens(ans),Vacuum_dens(ans));
	if (internal_energy(ans) < Min_energy(ans))
	    Energy(ans) = kinetic_energy(ans)+Min_energy(ans);
	reset_gamma(ans);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	speed = fabs(vel(idir,ans)) + sound_speed(ans);
	set_max_wave_speed(idir,speed,ans,Coords(sten->p[0]),sten->wave);
	if (debugging("godunov"))
	{
	    verbose_print_state("ans",ans);
	}
	debug_print("godunov","Left godunov()\n");
}		/*end godunov*/


/*
*				godunovobl():
*
*      Godunov finite difference scheme in one space variable
*      applied to an arbitrary direction of the conservation form of the
*      equations of compressible flow.
*/

EXPORT void godunovobl(
	double		ds,
	double		dt,
	Tan_stencil	*sten,
	Locstate	ans,
	Front		*fr)
{
	COMPONENT	comp = sten->comp;
	Locstate        *sts = sten->states;
	const double	*dir = sten->dir;
	POINT		*p1 = sten->p[-1], *p2 = sten->p[0], *p3 = sten->p[1];
	Locstate	s1 = sts[-1], s2 = sts[0], s3 = sts[1];
	double		lms,speed;
	double		dir1[SMAXD], dir3[SMAXD];
	double		va, pra, vb, prb;
	double		ga, gb, g2;
	double		time = fr->time + 0.5*dt;
	int		i, dim = fr->interf->dim;
	double		v1, v2, v3;
#if defined(COMBUSTION_CODE)
	double		coef, rt2, rta, rtb;
	int		reaction;
#endif /* defined(COMBUSTION_CODE) */
	static size_t	sizest = 0;
	static boolean	is_grav = NO;
	static Locstate left = NULL,right = NULL,sa = NULL,sb = NULL;
	static double    a;
	static GEOMETRY geom;
	double           pr1, pr2, pr3;
	double           rad1, rad2, rad3, Ra, Rb, rmin;
	double           ka, kb;

	debug_print("godunov","Entered godunovobl()\n");
	if (is_obstacle_state(s2)) 
	{
	    g_obstacle_state(ans,fr->sizest);
	    debug_print("godunov","obstacle state, Left godunovobl()\n");
	    return;
	}
	if (sa == NULL) 
	{
	    sizest = Params(s2)->sizest;
	    (*Params(s2)->_alloc_state)(&left,sizest);
	    (*Params(s2)->_alloc_state)(&right,sizest);
	    (*Params(s2)->_alloc_state)(&sa,sizest);
	    (*Params(s2)->_alloc_state)(&sb,sizest);
	    set_type_of_state(sa,GAS_STATE);
	    set_type_of_state(sb,GAS_STATE);
	    is_grav = is_gravity();
	    if (is_rotational_symmetry())
		geom = Geometry(&a);
	}
	if (debugging("godunov"))
	{
	    (void) printf("Input data into godunovobl()\n");
	    (void) printf("ds = %g, dt = %g, comp = %d\n",ds,dt,comp);
	    for (i = 0; i < dim; ++i)
	    {
	        (void) printf("dir[%d] = %g%s",i,dir[i],
	                      (i == dim-1) ? "\n" : ", ");
	    }

	    (void) printf("hypersurface: hypersurface %llu\n",
	                  hypersurface_number(sten->newhs));
	    if (p1 != NULL)
	    {
	        for (i = 0; i < dim; ++i)
	        {
	            (void) printf("p1[%d] = %g%s",i,Coords(p1)[i],
	                          (i == dim-1) ? "\n" : ", ");
	        }
	    }
	    (void) printf("\n");
	    verbose_print_state("s1",s1);

	    if (p2 != NULL)
	    {
	        for (i = 0; i < dim; ++i)
	        {
	            (void) printf("p2[%d] = %g%s",i,Coords(p2)[i],
	                          (i == dim-1) ? "\n" : ", ");
	        }
	    }
	    (void) printf("\n");
	    verbose_print_state("s2",s2);

	    if (p3 != NULL)
	    {
	        for (i = 0; i < dim; ++i)
	        {
	            (void) printf("p3[%d] = %g%s",i,Coords(p3)[i],
	                          (i == dim-1) ? "\n" : ", ");
	        }
	    }
	    (void) printf("\n");
	    verbose_print_state("s3",s3);
	}

	if (RegionIsFlowSpecified(ans,s2,Coords(p2),comp,comp,fr))
	{
	    debug_print("godunov","skip comp, Left godunovobl()\n");
	    return;
	}
	set_type_of_state(ans,GAS_STATE);


#if defined(COMBUSTION_CODE)
	if (sten->newhs != NULL && Composition_type(s2) == ZND) 
	{
	    reaction = (wave_type(sten->newhs) < FIRST_PHYSICS_WAVE_TYPE) ?
	               YES : NO;
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
	        for (i = 0; i < dim; ++i) dir1[i] = dir[i];
	    }
	    else
	    {
	        for ( i = 0; i < dim; ++i) dir1[i] /= len;
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
	        for (i = 0; i < dim; ++i) dir3[i] = dir[i];
	    }
	    else
	    {
	        for ( i = 0; i < dim; ++i) dir3[i] /= len;
	    }
	}
#endif /* defined(OLD) */

#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND) 
	{
	    v1 = scalar_product(VelocityVector(s1,NULL),dir1,dim);
	    v2 = scalar_product(VelocityVector(s2,NULL),dir,dim);
	    v3 = scalar_product(VelocityVector(s3,NULL),dir3,dim);
	}
#endif /* defined(COMBUSTION_CODE) */

	if (is_rotational_symmetry() && a > 0)
	{
	    RECT_GRID *gr = fr->rect_grid;
	    rmin = fabs(pos_radius(0.0,gr));
	    v1 = scalar_product(VelocityVector(s1,NULL),dir1,dim);
	    v2 = scalar_product(VelocityVector(s2,NULL),dir,dim);
	    v3 = scalar_product(VelocityVector(s3,NULL),dir3,dim);
	    pr1 = pressure(s1);
	    pr2 = pressure(s2);
	    pr3 = pressure(s3);
	    rad1 = pos_radius(Coords(p1)[0],gr);
	    rad2 = pos_radius(Coords(p2)[0],gr);
	    rad3 = pos_radius(Coords(p3)[0],gr);
	}

	lms = dt/ds;
	ft_assign(left,s1,fr->sizest);
	ft_assign(right,s2,fr->sizest);

	riemann_solution(0.0,dir1,left,right,sa,GAS_STATE);
	if (debugging("godunov"))
	{
	    verbose_print_state("sa",sa);
	}
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND) 
	{
	    rta = reaction ? (rt2 + coef*reaction_rate(s1)) : 0.0;
	    Prod(sa) = 0.5*(Prod(s2) + Prod(s1) -
	               lms*(v2*Prod(s2) - v1*Prod(s1)) + dt*rta);
	    Prod(sa) = max(Prod(sa),0.0);
	    Prod(sa) = min(Prod(sa),Dens(sa));
	}
#endif /* defined(COMBUSTION_CODE) */

	ft_assign(left,right,fr->sizest);
	ft_assign(right,s3,fr->sizest);

	riemann_solution(0.0,dir3,left,right,sb,GAS_STATE);
	if (debugging("godunov"))
	{
	    verbose_print_state("sb",sb);
	}
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND) 
	{ 
	    rtb = reaction ? (rt2 + coef*reaction_rate(s3)) : 0.0;
	    Prod(sb) = 0.5*(Prod(s3) + Prod(s2) -
	               lms*(v3*Prod(s3) - v2*Prod(s2)) + dt*rtb);
	    Prod(sb) = max(Prod(sb),0.0);
	    Prod(sa) = min(Prod(sb),Dens(sb));
	}
#endif /* defined(COMBUSTION_CODE) */

	if (is_rotational_symmetry() && a > 0)
	{
	    Ra = (geom == SPHERICAL) ?
			2.0*(rad1*rad1+rad1*rad2+rad2*rad2)/(3.0*(rad1+rad2)) :
			0.5*(rad1 + rad2);
	    Rb = (geom == SPHERICAL) ?
			2.0*(rad2*rad2+rad2*rad3+rad3*rad3)/(3.0*(rad2+rad3)) :
			0.5*(rad2 + rad3);
	    if (fabs(Ra) > rmin)
	    {
		Dens(sa) += - 0.25*a*dir[0]*dt*(Dens(s1)*v1 + Dens(s2)*v2)/Ra;
		Energy(sa) += - 0.25*a*dir[0]*dt*
		              ((Energy(s1)+pr1)*v1 + (Energy(s2)+pr2)*v2)/Ra;
	        for(i = 0; i < dim; ++i)
	            Mom(sa)[i] += -0.25*a*dir[0]*dt*
				   (Mom(s1)[i]*v1 + Mom(s2)[i]*v2)/Ra;
		reset_gamma(sa);
	    }
	    if (fabs(Rb) > rmin)
	    {
		Dens(sb) += - 0.25*a*dir[0]*dt*(Dens(s2)*v2 + Dens(s3)*v3)/Rb;
		Energy(sb) += - 0.25*a*dir[0]*dt*(
		                (Energy(s2)+pr2)*v2 + (Energy(s3)+pr3)*v3)/Rb;
		for(i = 0; i < dim; ++i)
		    Mom(sb)[i] += -0.25*a*dir[0]*dt*
				   (Mom(s2)[i]*v2 + Mom(s3)[i]*v3)/Rb;
		reset_gamma(sb);
	    }
	}

	if (is_grav == YES) 
	{
	    double    crdsa[3], crdsb[3];
	    if (debugging("godunov"))
	        (void) printf("Gravity added\n");
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
	                      (Dens(sa) + 0.5*(Dens(s2) + Dens(s1)));
	        Energy(sa) += 0.5*dt*ga*dir[i]*
	                      (Mom(sa)[i] + 0.5*(Mom(s2)[i]+Mom(s1)[i]));
	        Mom(sb)[i] += 0.5*dt*gb*dir[i]*
	                      (Dens(sb) + 0.5*(Dens(s3) + Dens(s2)));
	        Energy(sb) += 0.5*dt*gb*dir[i]*
	                      (Mom(sb)[i] + 0.5*(Mom(s3)[i]+Mom(s2)[i]));
	    }
            reset_gamma(sa);
            reset_gamma(sb);
	}

	va = vb = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    va += vel(i,sa)*dir[i];
	    vb += vel(i,sb)*dir[i];
	}
	pra = pressure(sa);
	prb = pressure(sb);

	Dens(ans) = Dens(s2) - lms*(vb*Dens(sb) - va*Dens(sa));
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(s2)->n_comps != 1)
            {
                for(i = 0; i < Params(s2)->n_comps; i++)
                {
                    pdens(ans)[i] = pdens(s2)[i] - lms*(vb*pdens(sb)[i] - va*pdens(sa)[i]);
                    /* New 051005 */
                    if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                        pdens(ans)[i] = 0.0;
                    else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                    {
                        printf("ERROR in godunovobl()\n");
                        printf("partial density < 0.0, case 0\n");
                        verbose_print_state("s1",s1);
                        verbose_print_state("s2",s2);
                        verbose_print_state("s3",s3);
                        verbose_print_state("ans",ans);
                        printf("point %14.12f, %14.12f\n", Coords(p2)[0], Coords(p2)[1]);
                        printf("LF point1 %14.12f, %14.12f\n", Coords(p1)[0], Coords(p1)[1]);
                        printf("LF point3 %14.12f, %14.12f\n", Coords(p3)[0], Coords(p3)[1]);
                        geomview_intfc_plot2d("gview_plot",
                              fr->interf,fr->rect_grid);
                        clean_up(ERROR);
                    }
                    /* End of New 051005 */
                }
            }
        }
	Energy(ans) = Energy(s2) -
	    lms*(vb*(Energy(sb) + prb) - va*(Energy(sa) + pra));
	for (i = 0; i < dim; ++i)
	{
	    Mom(ans)[i] = Mom(s2)[i] -
	                  lms*(vb*Mom(sb)[i]-va*Mom(sa)[i]+dir[i]*(prb-pra));
	}
	reset_gamma(ans);

#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND) 
	{
	    rta = reaction ? coef * reaction_rate(sa) : 0.0;
	    rtb = reaction ? coef * reaction_rate(sb) : 0.0;
	    Prod(ans) = Prod(s2) - lms*(vb*Prod(sb)-va*Prod(sa)) + dt*(rtb+rta);
	    Prod(ans) = max(Prod(ans),0.0);
	    Prod(ans) = min(Prod(ans),Dens(ans));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(ans,s2);
	if (is_rotational_symmetry() && a > 0.0)
	{
	    ka = (fabs(Ra) > rmin) ? -0.5*a*dir[0]*dt*va/Ra : 0.0;
	    kb = (fabs(Rb) > rmin) ? -0.5*a*dir[0]*dt*vb/Rb : 0.0;
	    Dens(ans) += ka*Dens(sa) + kb*Dens(sb);
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                if(Params(s2)->n_comps != 1)
                {
                    for(i = 0; i < Params(s2)->n_comps; i++)
                    {
                        pdens(ans)[i] += ka*pdens(sa)[i] + kb*pdens(sb)[i];
                        /* New 051005 */
                        if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                            pdens(ans)[i] = 0.0;
                        else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                        {
                            printf("ERROR in godunovobl()\n");
                            printf("partial density < 0.0\n");
                            clean_up(ERROR);
                        }
                        /* End of New 051005 */
                    }
                }
            }
	    Energy(ans) += ka*(Energy(sa) + pra) + kb*(Energy(sb) + prb);
	    for (i = 0; i < dim; ++i)
	        Mom(ans)[i] += ka*Mom(sa)[i] + kb*Mom(sb)[i];
	    reset_gamma(ans);
	}

	if (is_grav == YES) 
	{
	    g2 = 0.5*scalar_product(gravity(Coords(sten->p[0]),time),dir,dim);
	    for (i = 0; i < dim; ++i)
	    {
	        Mom(ans)[i] += dt*g2*dir[i]*(Dens(s2) + Dens(ans));
	        Energy(ans) += dt*g2*dir[i]*(Mom(s2)[i] + Mom(ans)[i]);
	    }
            reset_gamma(ans);
	}
#if defined(HYP_DBG)
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Scaling partial density */
            if(Params(s2)->n_comps != 1)
            {
                double sum = 0.0;
                for(i = 0; i < Params(s2)->n_comps; i++)
                    sum += pdens(ans)[i];
                for(i = 0; i < Params(s2)->n_comps; i++)
                    pdens(ans)[i] = pdens(ans)[i]/sum*Dens(ans);
            }
        }
	if (!check_gas("godunovobl",sts,ans,sten,YES,fr))
	{
	    (void) printf("WARNING in godunovobl(), bad state detected\n");
	    print_general_vector("dir = ",dir,dim,"\n");
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	    verbose_print_state("ans",ans);
	    LFoblique(ds,dt,sten,ans,fr);
	}
#endif /* defined(HYP_DBG) */
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Dens(ans) = max(Dens(ans),Vacuum_dens(ans));
	if (internal_energy(ans) < Min_energy(ans))
	    Energy(ans) = kinetic_energy(ans)+Min_energy(ans);
	reset_gamma(ans);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	
	speed = 0.0;
	for (i = 0; i < dim; ++i)
	    speed += dir[i]*vel(i,ans);
	speed = fabs(speed) + sound_speed(ans);
	for (i = 0; i < dim; ++i)
	{
	    set_max_front_speed(i,fabs(dir[i]*speed),ans,Coords(p2),fr);
	}
	set_max_front_speed(dim,speed/ds,ans,Coords(p2),fr);

	if (debugging("godunov")) verbose_print_state("ans",ans);

	debug_print("godunov","Left godunovobl()\n");

}		/*end godunovobl*/
