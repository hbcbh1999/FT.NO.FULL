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
*				glf.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Implementation of the Lax-Friedrichs method.
*
*/

#include <ghyp/ghyp.h>

#define HYP_DBG

/*
*				LF():
*
*      Lax_Friedrichs finite difference scheme in one space variable
*      applied to the t-x part of the conservation form of the
*      equations of compressible flow.
*
*/

/*ARGSUSED*/
EXPORT void LF(
	double		dh,
	double		dt,
	Locstate	ans,
	const double	*dir,
	int		swp_num,
	int		*iperm,
	int		*index,
	Stencil		*sten)
{
	COMPONENT	comp = sten->newcomp;
	Locstate	s1 = sten->st[-1], s2 = sten->st[0], s3 = sten->st[1];
	int		i, dim = sten->fr->interf->dim;
	int             idir = iperm[swp_num];
	double		lm, speed;
	double		g2;
	double		time = sten->fr->time + 0.5*dt;
#if defined(COMBUSTION_CODE)
	double		rt1, rt3;
#endif /* defined(COMBUSTION_CODE) */
#if defined(HYP_DBG)
	static	char	fname[3][12] = {"LFx()","LFy()","LFz()"};
#endif /* defined(HYP_DBG) */
	static boolean	is_grav;
	double           v1, pr1, v3, pr3, v2;
	static double    a;
	double           rad1, rad3, rmin;
	double           k1, k3;
	static boolean     first = YES;

	if (first)
	{
	    first = NO;
	    is_grav = is_gravity();
	    if (is_rotational_symmetry())
		(void) Geometry(&a);
	}

	if (is_obstacle_state(s2)) 
	{
	    g_obstacle_state(ans,sten->fr->sizest);
	    return;
	}
	set_type_of_state(ans,GAS_STATE);

	if (is_grav == YES)
	    g2 = 0.5*gravity(Coords(sten->p[0]),time)[idir];
	else
	    g2 = 0.0;
	v1 = vel(idir,s1);
	pr1 = pressure(s1);
	v3 = vel(idir,s3);
	pr3 = pressure(s3);
	v2 = vel(idir,s2);
	if (is_rotational_symmetry() && (a > 0.0) && (idir == 0))
	{
	    RECT_GRID *gr = sten->fr->rect_grid;
	    rad1 = pos_radius(Coords(sten->p[-1])[0],gr);
	    rad3 = pos_radius(Coords(sten->p[1])[0],gr);
	    rmin = fabs(pos_radius(0.0,gr));
	}

	lm = 0.5*dt/dh;

	/* Compute final answer */

	Dens(ans) = 0.5*(Dens(s1)+Dens(s3))-lm*(Mom(s3)[idir] - Mom(s1)[idir]);
	Energy(ans) = 0.5*(Energy(s1)+Energy(s3)) -
	              lm*(v3*(Energy(s3) + pr3) - v1*(Energy(s1) + pr1));
	for (i = 0; i < dim; ++i)
	{
	    Mom(ans)[i] = 0.5*(Mom(s1)[i]+Mom(s3)[i]) -
	                  lm*(v3*Mom(s3)[i] - v1*Mom(s1)[i]);
	}
	Mom(ans)[idir] -= lm*(pr3 - pr1);

	if (is_rotational_symmetry() && (a > 0.0) && (idir == 0))
	{
	    k1 = (fabs(rad1) > rmin) ? -0.5*a*dt*v1/rad1 : 0.0;
	    k3 = (fabs(rad3) > rmin) ? -0.5*a*dt*v3/rad3 : 0.0;
	    Dens(ans) += k1*Dens(s1) + k3*Dens(s3);
	    Energy(ans) += k1*(Energy(s1) + pr1) + k3*(Energy(s3) + pr3);
	    for(i = 0; i < dim; ++i)
		Mom(ans)[i] += k1*Mom(s1)[i] + k3*Mom(s3)[i];
	}

	reset_gamma(ans);
#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND) 
	{
	    rt1 = 0.5 * reaction_rate(s1);
	    rt3 = 0.5 * reaction_rate(s3);
	    Prod(ans) = 0.5*(Prod(s1)+Prod(s3))-lm*(v3*Prod(s3)-v1*Prod(s1)) +
	                dt*(rt1+rt3);
	    Prod(ans) = max(Prod(ans),0.0);
	    Prod(ans) = min(Prod(ans),Dens(ans));
	}
#endif /* defined(COMBUSTION_CODE) */
	Set_params(ans,s2);
	if (is_grav == YES) 
	{
	    Mom(ans)[idir] += dt*(Dens(s2) + Dens(ans))*g2;
	    Energy(ans) += dt*(Mom(s2)[idir] + Mom(ans)[idir])*g2;
	}
#if defined(HYP_DBG)
	if (!check_ans(fname[idir],dh,dt,ans,comp,sten,YES))
	{
	    (void) printf("WARNING in LF(), bad state detected\n");
	    print_general_vector("dir = ",dir,dim," ");
	    (void) printf("idir = %d\n",idir);
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	    verbose_print_state("ans",ans);
	}
#endif /* defined(HYP_DBG) */

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Dens(ans) = max(Dens(ans),Vacuum_dens(ans));
	if (internal_energy(ans) < Min_energy(ans))
	    Energy(ans) = kinetic_energy(ans)+Min_energy(ans);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	speed = fabs(vel(idir,ans)) + sound_speed(ans);
	reset_gamma(ans);
	set_max_wave_speed(idir,speed,ans,Coords(sten->p[0]),sten->wave);
}		/*end LF*/


/*
*				LFoblique():
*
*      Lax-Friedrichs finite difference scheme in one space variable
*      applied to an arbitrary direction of the conservation form of the
*      equations of compressible flow.
*/

EXPORT void LFoblique(
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
	double		g2;
	double		time = fr->time + 0.5*dt;
	int		i, dim = fr->interf->dim;
	double		v1, v3;
	double           pr1, pr3;
#if defined(COMBUSTION_CODE)
	double		coef, rt2, rt1, rt3;
	int		reaction;
#endif /* defined(COMBUSTION_CODE) */
	static boolean	is_grav = NO;
	static double    a;
	double           rad1, rad3, rmin;
	double           k1, k3;
	static boolean     first = YES;

	debug_print("LF","Entered LFoblique()\n");
	if (first)
	{
	    first = NO;
	    is_grav = is_gravity();
	    if (is_rotational_symmetry())
		(void) Geometry(&a);
	}

	if (is_obstacle_state(s2)) 
	{
	    g_obstacle_state(ans,fr->sizest);
	    debug_print("LF","obstacle state, Left LFoblique()\n");
	    return;
	}
	if (debugging("LF"))
	{
	    (void) printf("Input data into LFoblique()\n");
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
	    debug_print("LF","skip comp, Left LFoblique()\n");
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

	v1 = scalar_product(VelocityVector(s1,NULL),dir1,dim);
	v3 = scalar_product(VelocityVector(s3,NULL),dir3,dim);
	pr1 = pressure(s1);
	pr3 = pressure(s3);

	if (is_rotational_symmetry() && a > 0)
	{
	    RECT_GRID *gr = fr->rect_grid;
	    rmin = fabs(pos_radius(0.0,gr));
	    rad1 = pos_radius(Coords(p1)[0],gr);
	    rad3 = pos_radius(Coords(p3)[0],gr);
	}

	lms = 0.5*dt/ds;

	Dens(ans) = 0.5*(Dens(s1)+Dens(s3)) - lms*(v3*Dens(s3) - v1*Dens(s1));
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(s2)->n_comps != 1)
            {
                for(i = 0; i < Params(s2)->n_comps; i++)
                {
	            pdens(ans)[i] = 0.5*(pdens(s1)[i]+pdens(s3)[i]) 
                                   - lms*(v3*pdens(s3)[i] - v1*pdens(s1)[i]);
                    /* New 051005 */
                    if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                        pdens(ans)[i] = 0.0;
                    else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                    {
                        printf("ERROR in LFoblique()\n");
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
	Energy(ans) = 0.5*(Energy(s1)+Energy(s3)) -
	    lms*(v3*(Energy(s3) + pr3) - v1*(Energy(s1) + pr1));
	for (i = 0; i < dim; ++i)
	{
	    Mom(ans)[i] = 0.5*(Mom(s1)[i]+Mom(s3)[i]) -
	                  lms*(v3*Mom(s3)[i]-v1*Mom(s1)[i]+dir[i]*(pr3-pr1));
	}

#if defined(COMBUSTION_CODE)
	if (Composition_type(s2) == ZND) 
	{
	    rt1 = reaction ? coef * reaction_rate(s1) : 0.0;
	    rt3 = reaction ? coef * reaction_rate(s3) : 0.0;
	    Prod(ans) = 0.5*(Prod(s1)+Prod(s3)) -
	         lms*(v3*Prod(s3)-v1*Prod(s1)) + dt*(rt1+rt3);
	    Prod(ans) = max(Prod(ans),0.0);
	    Prod(ans) = min(Prod(ans),Dens(ans));
	}
#endif /* defined(COMBUSTION_CODE) */
        reset_gamma(ans);
	Set_params(ans,s2);
	if (is_rotational_symmetry() && a > 0.0)
	{
	    k1 = (fabs(rad1) > rmin) ? -0.5*a*dir[0]*dt*v1/rad1 : 0.0;
	    k3 = (fabs(rad3) > rmin) ? -0.5*a*dir[0]*dt*v3/rad3 : 0.0;
	    Dens(ans) += k1*Dens(s1) + k3*Dens(s3);
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    if(Params(s2)->n_comps != 1)
                    {
                        for(i = 0; i < Params(s2)->n_comps; i++)
                        {
                            pdens(ans)[i] += k1*pdens(s1)[i] + k3*pdens(s3)[i];
                            /* New 051005 */
                            if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                                pdens(ans)[i] = 0.0;
                            else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                            {
                                printf("ERROR in LFoblique()\n");
                                printf("partial density < 0.0\n");
                                clean_up(ERROR);
                            }
                            /* End of New 051005 */
                        }
                    }
                }
	    Energy(ans) += k1*(Energy(s1) + pr1) + k3*(Energy(s3) + pr3);
	    for (i = 0; i < dim; ++i)
	        Mom(ans)[i] += k1*Mom(s1)[i] + k3*Mom(s3)[i];
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
	if (!check_gas("LFoblique",sts,ans,sten,YES,fr))
	{
	    (void) printf("WARNING in LFoblique(), bad state detected\n");
	    print_general_vector("dir = ",dir,dim,"\n");
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	    verbose_print_state("ans",ans);
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

	if (debugging("LF")) verbose_print_state("ans",ans);

	debug_print("LF","Left LFoblique()\n");

}		/*end LFoblique*/
