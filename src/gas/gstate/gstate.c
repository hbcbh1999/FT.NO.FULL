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
*				gstate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routine for the evaluation of states.
*
*	The following routines are accessed through function pointers
*	set in init_physics():
*
*	g_solution() and g_intfc_solution() are called from dinout.c.
*/


#include <gdecs/gdecs.h>


	/*LOCAL Function Prototypes*/
LOCAL	double	compute_max(double,double*,double,double*,int);
LOCAL	double	compute_min(double,double*,double,double*,int);
LOCAL	void	initialize_extreme_values(EXTREME_VALUES*);
LOCAL	void	update_extreme_values(EXTREME_VALUES*,Locstate,double*,int);
LOCAL	void 	g_check_front_state_consistency2d(Front*);
#if defined(COMBUSTION_CODE)
LOCAL	double	prod(Locstate);
#endif /* defined(COMBUSTION_CODE) */

EXPORT	void	g_alloc_state(
	Locstate	*g_alloced_state,
	size_t		sizest)
{
	scalar(g_alloced_state,sizest);
	if (*g_alloced_state != NULL)
	    set_type_of_state(*g_alloced_state,UNKNOWN_STATE);
}		/*end g_alloc_state*/

EXPORT	Locstate	g_alloc_intfc_state(
	size_t          sizest)
{
	Locstate	newst;

	newst = (Locstate) store(sizest);
	if (newst != NULL)
	    set_type_of_state(newst,UNKNOWN_STATE);
	return newst;
}               /*end g_alloc_intfc_state*/

EXPORT	void	g_clear_state(
	Locstate	state,
	size_t		sizest)
{
	zero_scalar(state,sizest);
	set_type_of_state(state,UNKNOWN_STATE);
}		/*end g_clear_state*/


/*
*			g_reflect_state():
*
*	Reflects the given state about the plane of symmetry with
*	normal direction nor.
*      
*/

/*ARGSUSED*/
EXPORT void g_reflect_state(
	Locstate	state,	/* state being reflected */
	INTERFACE	*intfc, /* interface being reflected */
	double		*pt,	/* position of state being reflected */
	double		*p,	/* point on reflection plane */
	double		*nor)	/* unit normal to reflection plane */
{
	double	v[MAXD], vdn;
	int	i;
	int	dim = intfc->dim;
	int 	debug_flag = NO;

	if (size_of_state(intfc) == 0)
	    return;

	if (is_obstacle_state(state))
	    return;

	for (i = 0 ; i < dim ; ++i)
	{		
	    v[i] = vel(i,state);
	    Vel(state)[i] = 0.0;
	}
	vdn = 2.0*scalar_product(v,nor,dim);

	for (i = 0; i < dim; ++i)
	    v[i] -= vdn*nor[i];

	if (state_type(state) == GAS_STATE)
	{
	  for (i = 0 ; i < dim ; ++i)
	      Mom(state)[i] = Dens(state)*v[i];
	}
	else
	{
	  for (i = 0 ; i < dim ; ++i)
	      Vel(state)[i] = v[i];
	}

	reset_gamma(state);
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") &&
		(is_bad_state(state,YES,"g_reflect_state")))
	{
	    screen("ERROR in g_reflect_state(), bad state detected\n");
	    print_gas_state(state);
	    verbose_print_state("bad state",state);
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
}		/*end g_reflect_state*/

/*
*			g_reflect_and_stratify_state():
*
*	Reflects the given state about the plane of symmetry with
*	normal direction nor.  It then stratifies the thermodynamics
*	of the state to account for gravity.
*/

/*ARGSUSED*/
EXPORT void g_reflect_and_stratify_state(
	Locstate	state,	/* state being reflected */
	INTERFACE	*intfc, /* interface being reflected */
	double		*pt,	/* position of state being reflected */
	double		*p,	/* arbitrary point on reflection plane */
	double		*nor)	/* unit normal to reflection plane */
{
	CHART	*chart = current_chart();
	int	i, dim = intfc->dim;
	double	gn;		/* normal component of gravity */
	double	ds;		/* distance from pt to reflection plane */
	double	dp[MAXD];
	double	pr[MAXD];	/* Coordinates of reflected point */
	double	v[MAXD];
	double   gpt[MAXD], gpr[MAXD], gbar[MAXD];
	double	time = chart->grid->time;

	if (is_obstacle_state(state))
	    return;
	g_reflect_state(state,intfc,pt,p,nor);

	for (i = 0; i < dim; ++i)
	{
	    dp[i] = pt[i] - p[i];
	    v[i] = vel(i,state);
	}

	ds = 2.0*scalar_product(dp,nor,dim);
	for (i = 0; i < dim; ++i)
	    pr[i] = pt[i] - ds*nor[i];
	eval_gravity(pt,time,gpt);
	eval_gravity(pr,time,gpr);
	for (i = 0; i < dim; ++i)
	    gbar[i] = 0.5*(gpt[i]+gpr[i]);
	gn = scalar_product(gbar,nor,dim);
	stratified_state(intfc,state,ds,gn,state);
	add_velocity_to_state(state,v);
	set_state(state,GAS_STATE,state);
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") &&
		(is_bad_state(state,YES,"g_reflect_and_stratify_state")))
	{
	    screen("ERROR in g_reflect_and_stratify_state(), "
		   "bad state detected\n");
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
}		/*end g_reflect_and_stratify_state*/

EXPORT double energy(
	Locstate	state)
{
	switch (state_type(state))
	{
	case GAS_STATE:
	    return	Energy(state);

	case EGAS_STATE:
	    return Dens(state)*Energy(state) + kinetic_energy(state);

	case TGAS_STATE:
	case FGAS_STATE:
	    return (internal_energy(state) + kinetic_energy(state));

	case VGAS_STATE:
	    return Dens(state) * Int_en(state) + kinetic_energy(state);

	case UNKNOWN_STATE:
	    screen("ERROR in energy(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in energy(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR in energy(), no such state type %d\n",
	           state_type(state));
	    clean_up(ERROR);
	    break;
	}
	return ERROR_FLOAT; /* For lint */
}		/*end energy*/


/*
*	The functions mom() and vel() return the momentum or 
*	velocity of a state in a specified coordinate direction.
*/

EXPORT double mom(
	int		idir,
	Locstate	state)
{
	switch (state_type(state))
	{
	case	GAS_STATE:
	    return	Mom(state)[idir];

	case	TGAS_STATE:
	case	FGAS_STATE:
	case	EGAS_STATE:
	case	VGAS_STATE:
	    return	Dens(state)*Vel(state)[idir];

	case UNKNOWN_STATE:
	    screen("ERROR in mom(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in mom(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR: in mom(), no such state type %d\n",
	           state_type(state));
	    clean_up(ERROR);
	    break;
	}
	return ERROR_FLOAT; /* For lint */
}		/*end mom*/

EXPORT	double	*MomentumVector(
	Locstate	state,
	double		*mstate)
{
	static	double	m[MAXD];
	int	i, dim = Params(state)->dim;

	if (mstate == NULL)
	    mstate = m;

	switch (state_type(state))
	{
	case	GAS_STATE:
	    for (i = 0; i < dim; ++i)
	    	mstate[i] = Mom(state)[i];
	    break;

	case	TGAS_STATE:
	case	FGAS_STATE:
	case	EGAS_STATE:
	case	VGAS_STATE:
	    for (i = 0; i < dim; ++i)
	    	mstate[i] = Dens(state)*Vel(state)[i];
	    break;

	case UNKNOWN_STATE:
	    screen("ERROR in MomentumVector(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in MomentumVector(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR in MomentumVector(), no such state type %d\n",
	    	   state_type(state));
	    clean_up(ERROR);
	    break;
	}
	return mstate;
}		/*end MomentumVector*/


EXPORT double vel(
	int		idir,
	Locstate	state)
{
	if (is_obstacle_state(state))
	    return 0.0;
	switch (state_type(state))
	{
	case	GAS_STATE:
	    return	Mom(state)[idir]/Dens(state);

	case	TGAS_STATE:
	case	FGAS_STATE:
	case	EGAS_STATE:
	case	VGAS_STATE:
	    return	Vel(state)[idir];

	case UNKNOWN_STATE:
	    screen("ERROR in vel(), unknown state type\n");
	    g_print_state(state);
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in vel(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR: in vel(), no such state type %d\n",
	           state_type(state));
	    clean_up(ERROR);
	}
	return ERROR_FLOAT; /* For lint */
}		/*end vel*/

EXPORT	double	*VelocityVector(
	Locstate	state,
	double		*vstate)
{
	static	double	v[MAXD];
	int	i, dim;

	if (vstate == NULL)
	    vstate = v;

	if (is_obstacle_state(state))
	{
	    dim = current_interface()->dim;
	    for (i = 0; i < dim; ++i)
	    	vstate[i] = 0.0;
	    return vstate;
	}

	dim = Params(state)->dim;
	switch (state_type(state))
	{
	case	GAS_STATE:
	    for (i = 0; i < dim; ++i)
	    	vstate[i] = Mom(state)[i]/Dens(state);
	    break;

	case	TGAS_STATE:
	case	FGAS_STATE:
	case	EGAS_STATE:
	case	VGAS_STATE:
	    for (i = 0; i < dim; ++i)
	    	vstate[i] = Vel(state)[i];
	    break;

	case UNKNOWN_STATE:
	    screen("ERROR in VelocityVector(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in VelocityVector(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR in VelocityVector(), no such state type %d\n",
	    	   state_type(state));
	    clean_up(ERROR);
	    break;
	}
	return vstate;
}		/*end VelocityVector*/

EXPORT	boolean	is_bad_state(
	Locstate   state,
	boolean       print_warning,
	const char *function)
{
	boolean is_bad;

	if (state == NULL)
	    is_bad = YES;
	else if (is_obstacle_state(state))
	    is_bad = NO;
	else
	    is_bad = invalid_state(function,state,print_warning);

	if (is_bad && debugging("fatalbad"))
	{
	    screen("ERROR in is_bad_state(), "
		   "bad state found in %s\n",function);
	    fprint_raw_gas_data(stdout,state,current_interface()->dim);
	    clean_up(ERROR);
	}
	return is_bad;
}		/*end is_bad_state*/

EXPORT	boolean	g_invalid_state(
	const char *function,
	Locstate   state,
	boolean       print_warning)
{
    	double      *m;
	double      E, minE, P, minP;
	int	   i, dim;
	const char *warn = "WARNING in g_invalid_state()";

	switch (state_type(state))
	{
	case	GAS_STATE:
	case	TGAS_STATE:
	case	FGAS_STATE:
	case	EGAS_STATE:
	case	VGAS_STATE:
	    break;

	case UNKNOWN_STATE:
	    if (print_warning)
	        (void) printf("%s, UNKNOWN_STATE state type detected in %s\n",
			      warn,function);
	    return YES;

	case OBSTACLE_STATE:
	default:
	    if (print_warning)
	        (void) printf("%s, OBSTACLE_STATE state type detected in %s\n",
			      warn,function);
	    return YES;
	}

	if (isnan(Dens(state)))
	{
	    (void) printf("%s, Dens(state) = %g is a NaN in %s\n",
		          warn,Dens(state),function);
	    return YES;
	}
	if (isnan(Energy(state)))
	{
	    (void) printf("%s, Energy(state) = %g is a NaN in %s\n",
		          warn,Energy(state),function);
	    return YES;
	}

	dim = Params(state)->dim;
	m = Mom(state);
	for (i = 0; i < dim; ++i)
	{
	    if (isnan(m[i]))
	    {
	        (void) printf("%s, Mom(state)[%d] = %g is a NaN in %s\n",
			      warn,i,m[i],function);
	        return YES;
	    }
	}

	if (Dens(state) < 0.0)
	{
	    if (print_warning)
	        (void) printf("%s, Dens(state) = %24.20g is negative in %s\n",
			      warn,Dens(state),function);
	    return YES;
	}
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(state);
            double pfrac_sum; 
            if(params->n_comps != 1)
            {
                pfrac_sum = 0.0;
                for(i = 0; i < params->n_comps; i++)
                {
                    pfrac_sum += pdens(state)[i]/Dens(state);
                    /* New 051005 */
                    if(fabs(pdens(state)[i]) > 10.0*MACH_EPS && pdens(state)[i] < 0.0)
                    {
	                if (print_warning)
                            (void) printf("%s, partial Dens[%d] = %24.20g < 0.0 %s\n",
                                     warn,i,pdens(state)[i],function);
                        return YES;
                    }
                    /* End of New 051005 */
                }
                if(fabs(1.0-pfrac_sum) > MACH_EPS*100)
                {
	            if (print_warning)
                        (void) printf("%s, Dens(state) = %24.20g not equal partial dens %s\n",
                                     warn,Dens(state),function);
                    printf("pfrac_sum = %24.20g\n",pfrac_sum);
                    return YES;
                }
            }
        }
	switch (state_type(state))
	{
	case	GAS_STATE:
	case	EGAS_STATE:
	case	FGAS_STATE:
	    E = internal_energy(state);	minE = Min_energy(state);
	    if (E < minE)
	    {
		if ((minE - E) > MACH_EPS)
		{
	            if (print_warning)
	                (void) printf("%s, internal_energy(state) = %24.20g < "
			              "Min_energy(state) = %24.20g in %s\n"
			              "Min_energy(state) - "
				      "internal_energy(state) = "
			              "%24.20g, MACH_EPS = %24.20g\n",
			              warn,E,minE,function,minE-E,MACH_EPS);
	            return YES;
		}
	    }
	    break;
	case	TGAS_STATE:
	case	VGAS_STATE:
	    P = Press(state);	minP = Min_pressure(state);
	    if (P < minP)
	    {
		if ((minP - P) > MACH_EPS)
		{
	            if (print_warning)
	                (void) printf("%s, pressure(state) = %24.20g < "
			              "Min_pressure(state) = %24.20g in %s\n"
			              "Min_pressure(state) - "
				      "pressure(state) = "
			              "%24.20g, eps = %24.20g\n",
			              warn,P,minP,function,minP-P,MACH_EPS);
	            return YES;
		}
	    }
	    break;
	}

#if defined(COMBUSTION_CODE)
	switch(Composition_type(state))
	{
	case ZND:
	    if (isnan(Prod(state)))
	    {
	        (void) printf("%s, Prod(state) = %g is a NaN in %s\n",
			      warn,Prod(state),function);
	        return YES;
	    }
	    if (prod(state) < 0.0)
	    {
		if (print_warning)
	            (void) printf("%s, prod(state) = %24.20g is negative in %s\n",
				  warn,prod(state),function);
	        return YES;
	    }
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    if (isnan(Prod(state)))
	    {
	        (void) printf("%s, Prod(state) = %g is a NaN in %s\n",
			      warn,Prod(state),function);
	        return YES;
	    }
	    if (isnan(Dens1(state)))
	    {
	        (void) printf("%s, Dens1(state) = %g is a NaN in %s\n",
			      warn,Dens1(state),function);
	        return YES;
	    }
	    if (prod(state) < 0.0)
	    {
		if (print_warning)
	            (void) printf("%s, prod(state) = %24.20g is negative in %s\n",
				  warn,prod(state),function);
	        return YES;
	    }
	    if (Dens1(state) < 0.0)
	    {
	        if (print_warning)
	            (void) printf("%s, Dens1(state) = %24.20g is negative in %s\n",
				  warn,Dens1(state),function);
	        return YES;
	    }
	    break;
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	default:
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */

	return NO;
}		/*end g_invalid_state*/

EXPORT	boolean	g_bad_state_data(
	const char *function,
	Front      *front,
	Wave       *wave)
{
	INTERFACE          *intfc = front->interf;
	POINT              *p;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF         *hs;
	Locstate           s, sl, sr;
	RECT_GRID	   *gr = front->rect_grid;
	boolean            bad = NO;
	int                dim = front->rect_grid->dim;
	int                icoords[3];
	int                i, j, k, imin, jmin, kmin, imax, jmax, kmax;

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    slsr(p,hse,hs,&sl,&sr);
	    if (is_bad_state(sl,YES,function))
	    {
		bad = YES;
		print_general_vector("Left state at point ",Coords(p),dim,
				     "is bad\n");
		fprint_raw_gas_data(stdout,sl,dim);
	    }
	    if (is_bad_state(sr,YES,function))
	    {
		bad = YES;
		print_general_vector("Right state at point ",Coords(p),dim,
				     "is bad\n");
		fprint_raw_gas_data(stdout,sr,dim);
	    }
	}
	imin = (gr->lbuf[0] > 0) ? -gr->lbuf[0] : 0;
	imax = (gr->ubuf[0] > 0) ?  gr->gmax[0] + gr->ubuf[0] : gr->gmax[0];
	if (dim > 1)
	{
	    jmin = (gr->lbuf[1] > 0) ? -gr->lbuf[1] : 0;
	    jmax = (gr->ubuf[1] > 0) ?  gr->gmax[1] + gr->ubuf[1] : gr->gmax[1];
	}
	else
	{
	    jmin = 0;
	    jmax = 1;
	}
	if (dim > 2)
	{
	    kmin = (gr->lbuf[2] > 0) ? -gr->lbuf[2] : 0;
	    kmax = (gr->ubuf[2] > 0) ?  gr->gmax[2] + gr->ubuf[2] : gr->gmax[2];
	}
	else
	{
	    kmin = 0;
	    kmax = 1;
	}
	for (k = kmin; k < kmax; ++k)
	{
	    icoords[2] = k;
	    for (j = jmin; j < jmax; ++j)
	    {
		icoords[1] = j;
		for (i = imin; i < imax; ++i)
		{
		    icoords[0] = i;
		    s = Rect_state(icoords,wave);
	            if (is_bad_state(s,YES,function))
	            {
		        bad = YES;
		        print_int_vector("State at icoords ",icoords,dim,
				             "is bad\n");
		        fprint_raw_gas_data(stdout,s,dim);
	            }
		}
	    }
	}
	return bad;
}		/*end g_bad_state_data*/

/*
*			g_principal_tangent():
*
*	Returns an appropriate tangent vector to a curve or surface.
*	For contacts and Neumann boundaries with non-zero shear this will be
*	the direction of velocity shear across the interface (note by
*	definition, the velocity behind a wall is zero.).  In all other
*	cases the tangent returned is given by the default value as
*	determined by the function f_principal_tangent().
*/

EXPORT	void	g_principal_tangent(
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	double			*nor,
	double			*vdir)
{
	Locstate	sl, sr;
	double		len, vl[3], vr[3];
	int		w_type = wave_type(hs);
	int		i, dim = hs->interface->dim;

	if ((!is_scalar_wave(w_type)) && (w_type != NEUMANN_BOUNDARY))
	{
	    f_principal_tangent(p,hse,hs,nor,vdir);
	    return;
	}

	slsr(p,hse,hs,&sl,&sr);
	(void) VelocityVector(sl,vl);
	(void) VelocityVector(sr,vr);
	for (i = 0; i < dim; ++i)
	    vdir[i] = vl[i] - vr[i];
	
	len = mag_vector(nor,dim);
	if (len > 0.0)
	{
	    double sp = scalar_product(vdir,nor,dim);
	    for (i = 0; i < dim; ++i)
	    	vdir[i] -= sp*nor[i]/(len*len);
	}
	len = mag_vector(vdir,dim);
	if (len == 0.0)
	{
	    f_principal_tangent(p,hse,hs,nor,vdir);
	    return;
	}
	for (i = 0; i < dim; ++i)
	    vdir[i] /= len;
}		/*end g_principal_tangent*/

EXPORT	double	RadialComponentOfVelocity(
	Locstate	state,
	double		*coords,
	int		dim)
{
	double	len;
	double	*v;

	if (is_obstacle_state(state))
	    return 0.0;

	len = mag_vector(coords,dim);
	v = VelocityVector(state,NULL);
	return (len == 0.0) ? 0.0 : scalar_product(coords,v,dim)/len;
}		/*end RadialComponentOfVelocity*/

/*
*			TangCmptOfXYVelocity():
*
*	Returns the tangential compenent (wrt coordinate origin) of the
*	projection of the veloctiy vector in the x-y plane.
*/
EXPORT	double	TangCmptOfXYVelocity(
	Locstate	state,
	double		*coords,
	int		dim)
{
	double	len;
	double	*v, tngt[MAXD];

	if (is_obstacle_state(state) || (dim != 2))
	    return 0.0;

	len = mag_vector(coords,dim);
	tngt[0] = -coords[1];	tngt[1] =  coords[0];

	v = VelocityVector(state,NULL);
	return (len == 0.0) ? 0.0 : scalar_product(tngt,v,dim)/len;
}		/*end TangCmptOfXYVelocity*/


EXPORT	void zero_state_velocity(
	Locstate	state,
	int		dim)
{
	double  E;
	int    i;

	if (is_obstacle_state(state))
	{
	    for (i = 0; i < dim; ++i)
	        Vel(state)[i] = 0.0;
	    return;
	}
	switch (state_type(state))
	{
	case GAS_STATE:
#if defined(COMBUSTION_CODE)
	case ZGAS_STATE:
	case CGAS_STATE:
#endif /*defined(COMBUSTION_CODE)*/
	    E = internal_energy(state);
	    for (i = 0; i < dim; ++i)
	        Mom(state)[i] = 0.0;
	    Energy(state) = E;
            reset_gamma(state);
	    break;
	case TGAS_STATE:
	case EGAS_STATE:
	case FGAS_STATE:
	    for (i = 0; i < dim; ++i)
	        Vel(state)[i] = 0.0;
	    break;
	case VGAS_STATE:
	    for (i = 0; i < dim; ++i)
	        Vel(state)[i] = 0.0;
	    set_type_of_state(state,TGAS_STATE);
	    set_state(state,VGAS_STATE,state);
	    break;
	case OBSTACLE_STATE:
	    break;
	case UNKNOWN_STATE:
	default:
	    for (i = 0; i < dim; ++i)
	        Vel(state)[i] = 0.0;
	    break;
	}
}		/*end zero_state_velocity*/

EXPORT	void	zero_normal_velocity(
	Locstate	state,
	double		*nor,
	int		dim)
{
	int		i;
	double		v[MAXD];
	double		s;

	for (i = 0; i < dim; ++i)
	    v[i] = vel(i,state);
	s = scalar_product(v,nor,dim);
	for (i = 0; i < dim; ++i)
	    v[i] = -s*nor[i];
	add_velocity_to_state(state,v);
}		/*end zero_normal_velocity*/

EXPORT	void alpha_state_velocity(
	double           alpha,
	Locstate	state,
	int		dim)
{
	int    i;
	double		v[MAXD];

	if(debugging("andrea_no_slip"))
	{
	    zero_state_velocity(state,dim);
	    return;
	}

	if (is_obstacle_state(state))
	{
	    for (i = 0; i < dim; ++i)
	      Vel(state)[i] = 0;
	    return;
	}

	for (i = 0; i < dim; ++i)
	    v[i] = (alpha-1)*vel(i,state);
	add_velocity_to_state(state,v);
}		/*end alpha_state_velocity*/

/*
*			max_speed():
*
*      Computes fabs(velocity) + sound speed.
*/

EXPORT double max_speed(
	Locstate	state)
{
	double		ans = 0.0;
	static Locstate state_therm = NULL;

	if (is_obstacle_state(state)) return 0.0;
	if (state_therm == NULL)
	{
	    (*Params(state)->_alloc_state)(&state_therm,Params(state)->sizest);
	}
	set_state(state_therm,TGAS_STATE,state);
	ans = mag_vector(Vel(state_therm),Params(state)->dim);
	ans += sound_speed(state_therm);
	return ans;
}		/*end max_speed*/

EXPORT	double mach_number_squared(
	Locstate	state,
	double		*abs_v,
	double		*rel_v)
{
	int		i, dim;
	double		qsq, csq;
	double		vtmp[MAXD], *rv;
	static double	zero_v[MAXD];

	if (is_obstacle_state(state)) return 0.0;
	csq = sound_speed_squared(state);

	rv = (rel_v != NULL) ? rel_v : vtmp;
	if (abs_v == NULL) abs_v = zero_v;
	dim = Params(state)->dim;
	for (qsq = 0.0, i = 0; i < dim; ++i)
	{
	    rv[i] = vel(i,state) - abs_v[i];
	    qsq += sqr(rv[i]);
	}
	return qsq/csq;
}		/*end mach_number_squared*/

EXPORT	double mach_number(
	Locstate	 state,
	double		*abs_v)
{
	return sqrt(mach_number_squared(state,abs_v,(double *)NULL));
}		/*end mach_number*/


#if defined(COMBUSTION_CODE)
EXPORT	double prod(
	Locstate	state)
{
	double		prd;

	if (Composition_type(state) == PURE_NON_REACTIVE)
	    prd = 0.0;
	else if ((Composition_type(state) == PTFLAME) ||
		 (Composition_type(state) == THINFLAME))
	    prd = (Burned(state)) ? Dens(state) : 0.0;
	else if (Composition_type(state) == ZND)
	{
	    switch (state_type(state))
	    {
	    case GAS_STATE:
	    	prd = Prod(state);
	    	break;

	    case TGAS_STATE:
	    case	FGAS_STATE:
	    case EGAS_STATE:
	    case VGAS_STATE:
	    	prd = React(state) * Dens(state);
	    	break;

	    case UNKNOWN_STATE:
	        screen("ERROR in prod(), unknown state type\n");
	        clean_up(ERROR);
	        break;

	    case OBSTACLE_STATE:
	        screen("ERROR in prod(), obstacle state type\n");
	        clean_up(ERROR);
	        break;

	    default:
	    	screen("ERROR in prod(), no such state type %d\n",
	               state_type(state));
	    	prd = ERROR_FLOAT;
	    	clean_up(ERROR);
	    }
	}
	return prd;
}		/*end prod*/

EXPORT	double react(
	Locstate	state)
{
	double		rct;

	if (Composition_type(state) == PURE_NON_REACTIVE)
	    rct = 0.0;
	else if ((Composition_type(state) == PTFLAME) ||
		 (Composition_type(state) == THINFLAME))
	    rct = (Burned(state)) ? 1.0 : 0.0;
	else if (Composition_type(state) == ZND)
	{
	    switch (state_type(state))
	    {
	    case GAS_STATE:
	    	rct = min(Prod(state) / Dens(state),1.0);
	    	break;

	    case TGAS_STATE:
	    case FGAS_STATE:
	    case EGAS_STATE:
	    case VGAS_STATE:
	    	rct = React(state);
	    	break;

	    case UNKNOWN_STATE:
	        screen("ERROR in react(), unknown state type\n");
	    	rct = ERROR_FLOAT;
	        clean_up(ERROR);
	        break;

	    case OBSTACLE_STATE:
	        screen("ERROR in react(), obstacle state type\n");
	        clean_up(ERROR);
	        break;

	    default:
	    	screen("ERROR in react(), no such state type %d\n",
	               state_type(state));
	    	rct = ERROR_FLOAT;
	    	clean_up(ERROR);
	    }
	}
	return rct;
}		/*end react*/
#endif /* defined(COMBUSTION_CODE) */

EXPORT	void g_check_front_state_consistency(Front *front)
{
	switch (front->rect_grid->dim)
	{
	case 2: 
	    g_check_front_state_consistency2d(front);
	    break;
	case 3: 
	    /*g_check_front_state_consistency3d(front); */
	    break;
	}
}	/* end g_check_front_state_consistency */

LOCAL	void g_check_front_state_consistency2d(Front *front)
{
	INTERFACE *intfc = front->interf;
	CURVE **c,*curve;
	BOND *b;
	POINT *p;
	Locstate sl,sr;
	Gas_param *paramsl,*paramsr;
	printf("Begin checking front state consistency\n");
	for (c = intfc->curves; c && *c; ++c)
	{
	    curve = *c;
	    sl = left_start_state(curve);
	    sr = right_start_state(curve);
	    paramsl = Params(sl);
	    paramsr = Params(sr);
	    for (b = curve->first; b != curve->last; b = b->next)
	    {
		p = b->end;
		sl = left_state(p);
		sr = right_state(p);
		if (Params(sl) != paramsl || Params(sr) != paramsr)
		{
		    screen("Inconsistent curve state found\n");
		    print_curve(curve);
		    show_curve_states(curve);
		    clean_up(ERROR);
		}
	    }
	    sl = left_end_state(curve);
	    sr = right_end_state(curve);
	    if (Params(sl) != paramsl || Params(sr) != paramsr)
	    {
		screen("Inconsistent curve state found\n");
		print_curve(curve);
		show_curve_states(curve);
		clean_up(ERROR);
	    }
	}
	printf("Front states are consistent\n");
}	/* end g_check_front_state_consistency2d */


/*
*			g_set_state():
*
*	This routine takes as input state2 and the storage location
*	for state1.  It writes into the state1 storage location the state which
*	is equivalent to state2, but with the state_type st1_type.
*/

EXPORT	void g_set_state(
	Locstate	st1,
	int		st1_type,
	Locstate	st2)
{
	int i, dim;
	
#if defined(COMBUSTION_CODE)
	if (st1 == st2)
	    Local_gamma_set(st2) = NO;
#endif /* defined(COMBUSTION_CODE) */

	if ((is_obstacle_state(st2)) || (st1_type == OBSTACLE_STATE))
	{
	    g_obstacle_state(st1,g_sizest());
	    return;
	}
	else if (st1_type == state_type(st2))
	{
	    if (st1_type == VGAS_STATE)
	    	ft_assign(st1,st2,sizeof(VGas));
	    else
	    	ft_assign(st1,st2,Params(st2)->sizest);
	    return;
	}

	dim = Params(st2)->dim;
	Set_params(st1,st2);
	switch (st1_type)
	{
	case GAS_STATE:
	    Dens(st1) = Dens(st2);
	    Energy(st1) = energy(st2);
	    for (i = 0; i < dim; ++i)
	    	Mom(st1)[i] = mom(i,st2);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(st2) == ZND)
	    	Prod(st1) = prod(st2);
#endif /* defined(COMBUSTION_CODE) */
	    if (st1 != st2)
	        reset_gamma(st1);
	    break;

	case TGAS_STATE:
	    Press(st1) = pressure(st2);
	    for (i = 0; i < dim; ++i)
	    	Vel(st1)[i] = vel(i,st2);
	    Dens(st1) = Dens(st2);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(st1) == ZND)
	    	React(st1) = react(st2);
#endif /* defined(COMBUSTION_CODE) */
	    if (st1 != st2)
	        reset_gamma(st1);
	    break;

	case FGAS_STATE:
	    Temperature(st1) = temperature(st2);
	    for (i = 0; i < dim; ++i)
	    	Vel(st1)[i] = vel(i,st2);
	    Dens(st1) = Dens(st2);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(st1) == ZND)
	    	React(st1) = react(st2);
#endif /* defined(COMBUSTION_CODE) */
	    if (st1 != st2)
	        reset_gamma(st1);
	    break;

	case EGAS_STATE:
	    Energy(st1) = specific_internal_energy(st2);
	    for (i = 0; i < dim; ++i)
	    	Vel(st1)[i] = vel(i,st2);
	    Dens(st1) = Dens(st2);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(st1) == ZND)
	    	React(st1) = react(st2);
#endif /* defined(COMBUSTION_CODE) */
	    if (st1 != st2)
	        reset_gamma(st1);
	    break;

	case VGAS_STATE:
	    g_set_state(st1,TGAS_STATE,st2);
	    Int_en(st1) = specific_internal_energy(st1);
	    Entropy(st1) = entropy(st1);
	    Sound_speed(st1) = sound_speed(st1);
#if defined(VERBOSE_GAS_PLUS)
	    Enthalpy(st1) = specific_enthalpy(st1);
	    Temp(st1) = temperature(st1);
#endif /* defined(VERBOSE_GAS_PLUS) */
#if defined(PHASE_CODE)
	    Wave_curve(st1) = NULL;
#endif /* defined(PHASE_CODE) */
	    break;

	case UNKNOWN_STATE:
	    screen("ERROR in g_set_state(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in g_set_state(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR in g_set_state(), no such state type %d\n",st1_type);
	    clean_up(ERROR);
	    break;
	}
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
             Gas_param *params = Params(st2);
             int    num_comps;
             if((num_comps = params->n_comps) != 1)
             {
                 for(i = 0; i < num_comps; i++)
                     pdens(st1)[i] = pdens(st2)[i];
             }
        }
	set_type_of_state(st1,st1_type);
}		/*end g_set_state*/

EXPORT	void add_velocity_to_state(
	Locstate	state,
	double		*v)
{
	int		i;
	int		dim;
	double		en;

	if (is_obstacle_state(state))
	    return;

	dim = Params(state)->dim;
	switch (state_type(state))
	{
	case GAS_STATE:
	    en = Energy(state) - kinetic_energy(state);
	    for (i = 0; i < dim; ++i)
	    	Mom(state)[i] += Dens(state)*v[i];
	    Energy(state) = en + kinetic_energy(state);
            reset_gamma(state);
	    break;

	case TGAS_STATE:
	case FGAS_STATE:
	case EGAS_STATE:
	case VGAS_STATE:
	    for (i = 0; i < dim; ++i)
	    	Vel(state)[i] += v[i];
	    break;

	case UNKNOWN_STATE:
	    screen("ERROR in add_velocity_to_state(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in add_velocity_to_state(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR in add_velocity_to_state(), no such state type %d\n",
	           state_type(state));
	    clean_up(ERROR);
	}
}		/*end add_velocity_to_state*/


EXPORT	double kinetic_energy(
	Locstate	state)
{
	double		ke;
	int		i;
	int		dim;

	if (is_obstacle_state(state))
	    return 0.0;
	dim = Params(state)->dim;
	switch (state_type(state))
	{
	case GAS_STATE:
	    for (ke = 0.0, i = 0; i < dim; ++i)
	    	ke += sqr(Mom(state)[i]);
	    ke /= 2.0*Dens(state);
	    break;

	case EGAS_STATE:
	case FGAS_STATE:
	case TGAS_STATE:
	case VGAS_STATE:
	    for (ke = 0.0, i = 0; i < dim; ++i)
	    	ke += sqr(Vel(state)[i]);
	    ke *= 0.5*Dens(state);
	    break;

	case UNKNOWN_STATE:
	    screen("ERROR in kinetic_energy(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in kinetic_energy(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR in kinetic_energy(), no such state type %d\n",
	           state_type(state));
	    ke = ERROR_FLOAT;
	    clean_up(ERROR);
	    break;
	}
#if defined(SUBGRID)
        if(Params(state)->avisc.turbulence == YES)
        {
            for (i = 0; i < dim; ++i)
            {
                ke += 0.5*(Tau(state)[i][i]);
            }
        }
#endif /* defined SUBGRID */
	return ke;
}		/*end kinetic_energy*/

/*
*			g_transform_state():
*
*	Transforms the vector degrees of freedom, i.e. the momenta, of a gas 
*	state under a given affine linear coordinate transformation.  We assume
*	the Jacobian is 1, so that densities are not transformed.
*/

EXPORT void g_transform_state(
	Locstate	state,
	AFLIN		*aflin)
{
	double		m[SMAXD]; 
	int		i, j; 
	int		dim = Params(state)->dim;
 
	for (i = 0; i < dim; ++i)
	    m[i] = Mom(state)[i];
	for (i = 0; i < dim; ++i) 
	{ 
	    Mom(state)[i] = 0.0; 
	    for (j = 0; j < dim; ++j) 
	    	Mom(state)[i] += aflin->a[i][j] * m[j];
	}
}		/*end g_transform_state*/


#if defined(COMBUSTION_CODE)

/*
*				reaction_rate():
*
*	This function is a reaction rate function.
*/

EXPORT double reaction_rate(
	Locstate	state)
{
	double		temp, temp_inv;

	if (is_obstacle_state(state)) return 0.0;

	temp = temperature(state);
	if (temp <= EPS)
	    return 0.0;
	temp_inv = 1.0/temp;

	if (temp_inv * Params(state)->critical_temperature > 1.)
	    return 0.0;

	switch (state_type(state))
	{
	case GAS_STATE:
	    if (Dens(state) < Prod(state))
	        return 0.0;
	    return (Params(state)->rate_mult * (Dens(state) - Prod(state)) *
	    	exp(-Params(state)->activ_en * temp_inv));

	case EGAS_STATE:
	case TGAS_STATE:
	case FGAS_STATE:
	case VGAS_STATE:
	    return (Params(state)->rate_mult*Dens(state)* 
	    	(1. - React( state))*exp(-Params(state)->activ_en*temp_inv));

	case UNKNOWN_STATE:
	    screen("ERROR in arrienus(), unknown state type\n");
	    clean_up(ERROR);
	    break;

	case OBSTACLE_STATE:
	    screen("ERROR in arrienus(), obstacle state type\n");
	    clean_up(ERROR);
	    break;

	default:
	    screen("ERROR in arrienus(), no such state type %d\n",
	           state_type(state));
	    clean_up(ERROR);
	    break;
	}
	return 0.0; /* For lint */
}		/*end reaction_rate*/

/*
*		flame_velocity():
*
*	Returns the flame_velocity relative to the Vel(st1) or Vel(st2),
*	whichever is unburnt (if both unburnt, returns rel. to Vel(st1)).
*/

EXPORT double flame_velocity(
	Locstate st1,
	Locstate st2)
{
	if(is_obstacle_state(st1)||is_obstacle_state(st2))
	{
	    screen("ERROR in flame_velocity(), st%c is obstacle state\n",
			is_obstacle_state(st1)?'1':'2');
	    clean_up(ERROR);
	}

	if(Unburned(st1))
	{
	    if(Params(st1)->_flame_velocity!=NULL)
		return (*Params(st1)->_flame_velocity)(st1);
	    else
		return Params(st1)->flame_coeff[0];
	}
	else if(Burned(st1)&&Unburned(st2))
	{
	    if(Params(st2)->_flame_velocity!=NULL)
		return (*Params(st2)->_flame_velocity)(st2);
	    else
		return Params(st2)->flame_coeff[0];
	}
	else
	{
	    screen("ERROR in flame_velocity(), both states are Burned\n");
	    clean_up(ERROR);
	}

	return 0;
}		/* end flame_velocity */
#endif /* defined(COMBUSTION_CODE) */



/*
*		g_solution(), g_intfc_solution():
*
*	Routines for printout.	FOR USE BY dinout.c ONLY.
*	The correspondence between the "var" and the state
*	variable is set in init_printing.
*/

/*ARGSUSED*/
EXPORT OUTPUT_VALUE *g_solution(
	OUTPUT_SOLN	*os,
	double		*coords,
	int		*icoords)
{
	static OUTPUT_VALUE Sol;
	Wave		    *wave = (Wave *) ((POINTER*) os->extra)[1];
	Locstate	    state;
	int		    dim = wave->rect_grid->dim;
	int		    var = os->var;
#if defined(COMBUSTION_CODE)
	COMPOSITION_TYPE    ctype;
#endif /* defined(COMBUSTION_CODE) */
	
	Sol.utype = Float;
	state = Rect_state(icoords,wave);

#if defined(COMBUSTION_CODE)
	ctype = (COMPOSITION_TYPE)((is_obstacle_state(state)) ? 
		PURE_NON_REACTIVE : Composition_type(state));
#endif /* defined(COMBUSTION_CODE) */


	if (var == 0)
	    Sol.uval.fval = Dens(state);
	else if (var == 1)
	    Sol.uval.fval = Energy(state);
	else if ((2 <= var) && (var < (2+dim)))
	    Sol.uval.fval = Mom(state)[var-2];
#if defined(COMBUSTION_CODE)
	else if ((strcmp(os->name,"REACTION_PROGRESS") == 0) &&
	                                 (ctype == PTFLAME))
	{
	    if (is_obstacle_state(state))
	        Sol.uval.fval = 0.0;
	    else
	        Sol.uval.fval = (Burned(state)) ? 1. : 0.;
	}
	else if ((strcmp(os->name,"PRODUCT_DENSITY") == 0) && (ctype == ZND))
	    Sol.uval.fval = Prod(state);
	else if ((strcmp(os->name,"DENSITY_1") == 0) &&
	                                 (ctype == TWO_CONSTITUENT_REACTIVE))
	    Sol.uval.fval = Dens1(state);
#endif /* defined(COMBUSTION_CODE) */
	else if (strcmp(os->name,"X-VELOCITY") == 0)
	    Sol.uval.fval = vel(0,state);
	else if (strcmp(os->name,"Y-VELOCITY") == 0)
	    Sol.uval.fval = vel(1,state);
	else if (strcmp(os->name,"Z-VELOCITY") == 0)
	    Sol.uval.fval = vel(2,state);
	else if (strcmp(os->name,"PRESSURE") == 0)
	    Sol.uval.fval = (is_obstacle_state(state)) ? 0.0 : pressure(state);
	else if (strcmp(os->name,"SOUND_SPEED") == 0)
	    Sol.uval.fval = (is_obstacle_state(state)) ? 0.0:sound_speed(state);
	else if (strcmp(os->name,"TEMPERATURE") == 0)
	    Sol.uval.fval = (is_obstacle_state(state)) ? 0.0:temperature(state);
	else if (strcmp(os->name,"SPECIFIC_ENTROPY") == 0)
	    Sol.uval.fval = (is_obstacle_state(state)) ? 0.0 : entropy(state);
	else if (strcmp(os->name,"RADIAL_COMPONENT_OF_VELOCITY") == 0)
	    Sol.uval.fval = RadialComponentOfVelocity(state,coords,dim);
	else if (strcmp(os->name,"TANGENTIAL_COMPONENT_OF_XY_VELOCITY") == 0)
	    Sol.uval.fval = TangCmptOfXYVelocity(state,coords,dim);
	else if (strcmp(os->name,"EOS-PARAMS") == 0)
	{
	    Sol.utype = ULong;
	    Sol.uval.ulval = gas_param_number(Params(state));
	}
        else if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* screen("os->name = %s : g_solution() os->var %d\n", os->name, os->var); */
            if((Params(state) != NULL) && (Params(state)->n_comps != 1))
                Sol.uval.fval = pdens(state)[os->var-(dim+2)];
            else
                Sol.uval.fval = ERROR_FLOAT;
        }
	else
	{
	    screen("ERROR in g_solution(), unknown value of var %d\n",var);
	    Sol.uval.fval = ERROR_FLOAT;
	    clean_up(ERROR);
	}
	return &Sol;
}		/*end g_solution*/


/* ARGSUSED */
EXPORT void g_intfc_solution(
	OUTPUT_SOLN	   *os,
	POINT		   *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	   *hs,
	OUTPUT_VALUE	   *left,
	OUTPUT_VALUE	   *right)
{
	Locstate	 stl, str;
	int		 dim = hs->interface->dim;
	int		 var = os->var;
#if defined(COMBUSTION_CODE)
	COMPOSITION_TYPE ctype;
#endif /* defined(COMBUSTION_CODE) */

	slsr(p,hse,hs,&stl,&str);
	left->utype = right->utype = Float;
#if defined(COMBUSTION_CODE)
	if (!is_obstacle_state(stl)) 
	    ctype = (COMPOSITION_TYPE)Composition_type(stl);
	else if (!is_obstacle_state(str)) 
	    ctype = (COMPOSITION_TYPE)Composition_type(str);
	else
	    ctype = PURE_NON_REACTIVE;
#endif /* defined(COMBUSTION_CODE) */

	if (var == 0)
	{
	    left->uval.fval = Dens(stl);
	    right->uval.fval = Dens(str);
	}
	else if (var == 1)
	{
	    left->uval.fval = Energy(stl);
	    right->uval.fval = Energy(str);
	}
	else if ((2 <= var) && (var < (2+dim)))
	{
	    left->uval.fval = Mom(stl)[var-2];
	    right->uval.fval = Mom(str)[var-2];
	}
#if defined(COMBUSTION_CODE)
	else if ((strcmp(os->name,"REACTION_PROGRESS")==0) && (ctype==PTFLAME))
	{
	    if (is_obstacle_state(stl)) 
	    	left->uval.fval = 0.;
	    else
	    	left->uval.fval = (Burned(stl)) ? 1. : 0.;
	    if (is_obstacle_state(str)) 
	    	right->uval.fval = 0.;
	    else
	    	right->uval.fval = (Burned(str)) ? 1. : 0.;
	}
	else if ((strcmp(os->name,"PRODUCT_DENSITY") == 0) && (ctype == ZND))
	{
	    left->uval.fval = Prod(stl);
	    right->uval.fval = Prod(str);
	}
	else if ((strcmp(os->name,"DENSITY_1") == 0) &&
	         (ctype == TWO_CONSTITUENT_REACTIVE))
	{
	    left->uval.fval = Dens1(stl);
	    right->uval.fval = Dens1(str);
	}
#endif /* defined(COMBUSTION_CODE) */
	else if (strcmp(os->name,"X-VELOCITY") == 0)
	{
	    left->uval.fval = vel(0,stl);
	    right->uval.fval = vel(0,str);
	}
	else if (strcmp(os->name,"Y-VELOCITY") == 0)
	{
	    left->uval.fval = vel(1,stl);
	    right->uval.fval = vel(1,str);
	}
	else if (strcmp(os->name,"Z-VELOCITY") == 0)
	{
	    left->uval.fval = vel(2,stl);
	    right->uval.fval = vel(2,str);
	}
	else if (strcmp(os->name,"PRESSURE") == 0)
	{
	    left->uval.fval  = (is_obstacle_state(stl)) ? 0.0 : pressure(stl);
	    right->uval.fval = (is_obstacle_state(str)) ? 0.0 : pressure(str);
	}
	else if (strcmp(os->name,"SOUND_SPEED") == 0)
	{
	    left->uval.fval  = (is_obstacle_state(stl)) ? 0.0 : sound_speed(stl);
	    right->uval.fval = (is_obstacle_state(str)) ? 0.0 : sound_speed(str);
	}
	else if (strcmp(os->name,"TEMPERATURE") == 0)
	{
	    left->uval.fval  = (is_obstacle_state(stl)) ? 0.0 : temperature(stl);
	    right->uval.fval = (is_obstacle_state(str)) ? 0.0 : temperature(str);
	}
	else if (strcmp(os->name,"SPECIFIC_ENTROPY") == 0)
	{
	    left->uval.fval  = (is_obstacle_state(stl)) ? 0.0 : entropy(stl);
	    right->uval.fval = (is_obstacle_state(str)) ? 0.0 : entropy(str);
	}
	else if (strcmp(os->name,"RADIAL_COMPONENT_OF_VELOCITY") == 0)
	{
	    left->uval.fval  = RadialComponentOfVelocity(stl,Coords(p),dim);
	    right->uval.fval = RadialComponentOfVelocity(str,Coords(p),dim);
	}
	else if (strcmp(os->name,"TANGENTIAL_COMPONENT_OF_XY_VELOCITY") == 0)
	{
	    left->uval.fval  = TangCmptOfXYVelocity(stl,Coords(p),dim);
	    right->uval.fval = TangCmptOfXYVelocity(str,Coords(p),dim);
	}
	else if (strcmp(os->name,"EOS-PARAMS") == 0)
	{
	    left->utype = ULong;
	    right->utype = ULong;
	    left->uval.ulval = gas_param_number(Params(stl));
	    right->uval.ulval = gas_param_number(Params(str));
	}
        else if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(! is_obstacle_state(stl))
            {
                if(Params(stl)->n_comps != 1)
                    left->uval.fval = pdens(stl)[os->var-(dim+2)];
                else
                    left->uval.fval = ERROR_FLOAT;
            }
            if(! is_obstacle_state(str))
            {
                if(Params(str)->n_comps != 1)
                    right->uval.fval = pdens(str)[os->var-(dim+2)];
                else
                    right->uval.fval = ERROR_FLOAT;
            }
        }
	else
	{
	    screen("ERROR in g_intfc_solution(), "
	           "unknown value of var %d\n",var);
	    clean_up(ERROR);
	}
}		/*end g_intfc_solution*/

EXPORT	void	g_bundle_states(
	int		*gmin,
	int		*gmax,
	Wave		*wave,
	byte		*buf)
{
	int       num_params;
	size_t    i;
	Gas_param **prmslst;
	Locstate  state;
	int       dim = wave->rect_grid->dim;
	size_t    sizest = wave->sizest;
	int       ic[MAXD];

	debug_print("bundle_states","Entered g_bundle_states()\n");

	if (debugging("bundle_states"))
	{
	    print_int_vector("gmin = ",gmin,dim,", ");
	    print_int_vector("gmax = ",gmax,dim,"\n");
	}

	num_params = return_params_list(&prmslst);

	switch (dim)
	{
	case 1:
	{
	    int        ix, ixmin, ixmax;

	    ixmin = gmin[0];    ixmax = gmax[0];
	    for (ix = ixmin; ix < ixmax; ++ix)
	    {
	        ic[0] = ix;
	        state = (Locstate) buf;
	        ft_assign(state,Rect_state(ic,wave),sizest);
	        for (i = 0; i < num_params; ++i)
	            if (Params(state) == prmslst[i])
	                break;
	        Params(state) = (Gas_param *) i;
	        buf += sizest;
	    }
	    break;
	}
	case 2:
	{
	    int        ix, ixmin, ixmax;
	    int        iy, iymin, iymax;

	    ixmin = gmin[0];    ixmax = gmax[0];
	    iymin = gmin[1];    iymax = gmax[1];
	    for (iy = iymin; iy < iymax; ++iy)
	    {
	        ic[1] = iy;
	        for (ix = ixmin; ix < ixmax; ++ix)
	        {
	            ic[0] = ix;
	            state = (Locstate) buf;
	            ft_assign(state,Rect_state(ic,wave),sizest);
	            for (i = 0; i < num_params; ++i)
	                if (Params(state) == prmslst[i])
	                    break;
	            Params(state) = (Gas_param *) i;
	            buf += sizest;
	        }
	    }
	    break;
	}
	case 3:
	{
	    int        ix, ixmin, ixmax;
	    int        iy, iymin, iymax;
	    int        iz, izmin, izmax;

	    ixmin = gmin[0];    ixmax = gmax[0];
	    iymin = gmin[1];    iymax = gmax[1];
	    izmin = gmin[2];    izmax = gmax[2];
	    for (iz = izmin; iz < izmax; ++iz)
	    {
	        ic[2] = iz;
	        for (iy = iymin; iy < iymax; ++iy)
	        {
	            ic[1] = iy;
	            for (ix = ixmin; ix < ixmax; ++ix)
	            {
	                ic[0] = ix;
	                state = (Locstate) buf;
	                ft_assign(state,Rect_state(ic,wave),sizest);
	                for (i = 0; i < num_params; ++i)
	                {
	                    if (Params(state) == prmslst[i])
	                        break;
	                }
	                Params(state) = (Gas_param *) i;
	                buf += sizest;
	            }
	        }
	    }
	    break;
	}
	}

	debug_print("bundle_states","Left g_bundle_states()\n");
}        /*end g_bundle_states*/

LOCAL	int	return_params_for_comp(
	Gas_param	***prms,
	INTERFACE	*intfc)
{	
	static Gas_param   **params = NULL;
	static int	   maxlen = 0;
	int		   comp, ncomps, mincomp, maxcomp;
	
	if (intfc->dim != 3)
	    return -1;

	mincomp = min_component(intfc);
	maxcomp = max_component(intfc);
	ncomps = maxcomp - mincomp + 1;
	
	if (ncomps > maxlen)
	{
	    if (params != NULL)
		free(params);
	    maxlen = ncomps;
	    uni_array(&params,maxlen,sizeof(Gas_param*));
	}
	
	for (comp = mincomp; comp <= maxcomp; comp++)
	    params[comp-mincomp] = gas_params_for_comp(comp,intfc);

	*prms = params;
	return	mincomp;
}

/*TMP Cannot do it this way */
void    g_check_params_comp(POINTER, INTERFACE*);
void    g_check_params_comp(
	POINTER 	pwave,
	INTERFACE	*intfc)
{
	int		ix, ixmin, ixmax;
	int		iy, iymin, iymax;
	int		iz, izmin, izmax;
	Locstate	st;
	Gas_param	**prms;
	int		comp, mincomp, ic[3];
	Wave		*wave = (Wave*)pwave;
	RECT_GRID	*gr = wave->rect_grid;
	int		*gmax = gr->gmax, *lbuf = gr->lbuf, *ubuf = gr->ubuf;
	double		*L = gr->L, *h = gr->h, coords[3];

	mincomp = return_params_for_comp(&prms, intfc);
	if(mincomp == -1)
	    return;

	ixmin = -lbuf[0];	ixmax = gmax[0] + ubuf[0];
	iymin = -lbuf[1];	iymax = gmax[1] + ubuf[1];
	izmin = -lbuf[2];	izmax = gmax[2] + ubuf[2];
	
	for (iz = izmin; iz < izmax; ++iz)
	{
	    ic[2] = iz;
	    for (iy = iymin; iy < iymax; ++iy)
	    {
		ic[1] = iy;
		for (ix = ixmin; ix < ixmax; ++ix)
		{
		    ic[0] = ix;
		    st = Rect_state(ic,wave);
			
		    /*#check params */
		    comp = Rect_comp(ic,wave);
		    if(comp >=mincomp && Params(st) != prms[comp-mincomp])
		    {
			printf("WARNING g_check_params_comp, param and comp are incosistent, "
			       "prms(%d)=%d %d, Params(st)=%d %d\n", comp, prms[comp-mincomp], 
			       gas_param_number(prms[comp-mincomp]), Params(st), gas_param_number(Params(st)));
			printf("position = %d %d %d\n", ix, iy, iz);
			
			/*cell center coords */
			coords[0] = L[0] + ix*h[0] + 0.5*h[0];
			coords[1] = L[1] + iy*h[1] + 0.5*h[1];
			coords[2] = L[2] + iz*h[2] + 0.5*h[2];
		        print_general_vector("coords=", coords, 3, "\n");
			
			print_int_vector("gmax=", gmax, 3, "\n");
			print_int_vector("lbuf=", lbuf, 3, "\n");
			print_int_vector("ubuf=", ubuf, 3, "\n");
			printf("\n");
		    }
	    	}
	    }
	}
}

EXPORT	void	g_unbundle_states(
	int		*gmin,
	int		*gmax,
	Wave		*wave,
	byte		*buf)
{
	int		num_params;
	Gas_param	**prmslst;
	Locstate	r_st;
	size_t		sizest = wave->sizest;
	int		dim = wave->rect_grid->dim;
	int		ic[MAXD];
	
	debug_print("bundle_states","Entered g_unbundle_states()\n");
	
	num_params = return_params_list(&prmslst);

	if (debugging("bundle_states"))
	{
	    print_int_vector("gmin = ",gmin,dim,", ");
	    print_int_vector("gmax = ",gmax,dim,"\n");
	}

	switch (dim)
	{
	case 1:
	{
	    int		ix, ixmin, ixmax;

	    ixmin = gmin[0];	ixmax = gmax[0];
	    for (ix = ixmin; ix < ixmax; ++ix)
	    {
	        ic[0] = ix;
	        r_st = Rect_state(ic,wave);
	        ft_assign(r_st,buf,sizest);
	        Params(r_st) = prmslst[(size_t) Params(r_st)];
	        buf += sizest;
	    }
	    break;
	}
	case 2:
	{
	    int		ix, ixmin, ixmax;
	    int		iy, iymin, iymax;

	    ixmin = gmin[0];	ixmax = gmax[0];
	    iymin = gmin[1];	iymax = gmax[1];
	    for (iy = iymin; iy < iymax; ++iy)
	    {
	        ic[1] = iy;
	        for (ix = ixmin; ix < ixmax; ++ix)
	        {
	            ic[0] = ix;
	            r_st = Rect_state(ic,wave);
	            ft_assign(r_st,buf,sizest);
	            Params(r_st) = prmslst[(size_t) Params(r_st)];
	            buf += sizest;
	        }
	    }
	    break;
	}
	case 3:
	{
	    int		ix, ixmin, ixmax;
	    int		iy, iymin, iymax;
	    int		iz, izmin, izmax;

	    ixmin = gmin[0];	ixmax = gmax[0];
	    iymin = gmin[1];	iymax = gmax[1];
	    izmin = gmin[2];	izmax = gmax[2];
	    for (iz = izmin; iz < izmax; ++iz)
	    {
	        ic[2] = iz;
	        for (iy = iymin; iy < iymax; ++iy)
	        {
	            ic[1] = iy;
	            for (ix = ixmin; ix < ixmax; ++ix)
	            {
	                ic[0] = ix;
	                r_st = Rect_state(ic,wave);
			ft_assign(r_st,buf,sizest);
			Params(r_st) = prmslst[(size_t) Params(r_st)];
			buf += sizest;
	            }
	    	}
	    }
	    break;
	}
	}

	debug_print("bundle_states","Left g_unbundle_states()\n");
}		/*end g_unbundle_states*/

EXPORT	void	g_initialize_max_front_speed(
	Front	*fr)
{
	f_initialize_max_front_speed(fr);
	initialize_extreme_values(&extreme_front_vals(fr));
}		/*end g_initialize_max_front_speed*/

EXPORT	void	g_set_max_front_speed(
	int		i,
	double		spd,
	Locstate	state,
	double		*coords,
	Front		*fr)
{
	CHART *chart;
	EXTREME_VALUES	*val = &extreme_front_vals(fr);
	Grid  *grid;
	int   dim;

	f_set_max_front_speed(i,spd,state,coords,fr);
	chart = chart_of_front(fr);
	grid = chart->grid;
#if DONT_COMPILE /*TODO REMOVE*/
	static int last_step = INT_MIN;
	if (grid->step > last_step)
	{
	    initialize_extreme_values(val);
	    last_step = grid->step;
	}
#endif /*DONT_COMPILE*/
	dim = grid->rect_grid->dim;
	update_extreme_values(val,state,coords,dim);
}		/*end g_set_max_front_speed*/

EXPORT	void	g_initialize_max_wave_speed(
	Wave	*wave)
{
	h_initialize_max_wave_speed(wave);
	initialize_extreme_values(&extreme_wave_vals(wave));
}		/*end g_initialize_max_wave_speed*/

EXPORT	void	g_set_max_wave_speed(
	int		i,
	double		spd,
	Locstate	state,
	double		*coords,
	Wave		*wave)
{
	CHART *chart;
	EXTREME_VALUES	*val = &extreme_wave_vals(wave);
	Grid  *grid;
	int    dim;

	h_set_max_wave_speed(i,spd,state,coords,wave);
	chart = chart_of_wave(wave);
	grid = chart->grid;
#if DONT_COMPILE /*TODO REMOVE*/
	static int last_step = INT_MIN;
	if (grid->step > last_step)
	{
	    initialize_extreme_values(val);
	    last_step = grid->step;
	}
#endif /*DONT_COMPILE*/
	dim = grid->rect_grid->dim;
	update_extreme_values(val,state,coords,dim);
}		/*end g_set_max_wave_speed*/

LOCAL	void	update_extreme_values(
	EXTREME_VALUES	*val,
	Locstate        state,
	double           *coords,
	int             dim)
{
	double rho, e, p, ke, E;
	double *v, *m;
	int   i;

	/*#bjet2 */
	if (state==NULL || is_obstacle_state(state))
	    return;

	rho = Dens(state);
	val->rho_max = compute_max(rho,coords,
	                           val->rho_max,val->coords_rho_max,dim);
	val->rho_min = compute_min(rho,coords,
	                           val->rho_min,val->coords_rho_min,dim);

	e = specific_internal_energy(state);
	val->e_max = compute_max(e,coords,val->e_max,val->coords_e_max,dim);
	val->e_min = compute_min(e,coords,val->e_min,val->coords_e_min,dim);

	p = pressure(state);
	
	if( NO && fabs(p+169.145) < 1.0e-2 )
	{
	    printf("#extreme \n");
	    print_general_vector("#newp", coords, 3, "\n");
	    verbose_print_state("state ", state);
	    
	    clean_up(ERROR);
	}

	val->p_max = compute_max(p,coords,val->p_max,val->coords_p_max,dim);
	val->p_min = compute_min(p,coords,val->p_min,val->coords_p_min,dim);

	ke = kinetic_energy(state);
	val->ke_max = compute_max(ke,coords,val->ke_max,val->coords_ke_max,dim);
	val->ke_min = compute_min(ke,coords,val->ke_min,val->coords_ke_min,dim);

	E = energy(state);
	val->E_max = compute_max(E,coords,val->E_max,val->coords_E_max,dim);
	val->E_min = compute_min(E,coords,val->E_min,val->coords_E_min,dim);

	v = VelocityVector(state,NULL);
	m = MomentumVector(state,NULL);
	for (i = 0; i < dim; ++i)
	{
	    val->v_max[i] = compute_max(v[i],coords,val->v_max[i],
	                                val->coords_v_max[i],dim);
	    val->v_min[i] = compute_min(v[i],coords,val->v_min[i],
	                                val->coords_v_min[i],dim);
	    val->m_max[i] = compute_max(m[i],coords,val->m_max[i],
	                                val->coords_m_max[i],dim);
	    val->m_min[i] = compute_min(m[i],coords,val->m_min[i],
	                                val->coords_m_min[i],dim);
	}
}		/*end update_extreme_values*/


LOCAL	double	compute_max(
	double val,
	double *coords,
	double current_max_val,
	double *current_max_val_coords,
	int   dim)
{
	int i;
	if (val > current_max_val)
	{
	    current_max_val = val;
	    for (i = 0; i < dim; ++i)
	        current_max_val_coords[i] = coords[i];
	}
	return current_max_val;
}		/*end compute_max*/

LOCAL	double	compute_min(
	double val,
	double *coords,
	double current_min_val,
	double *current_min_val_coords,
	int   dim)
{
	int i;
	if (val < current_min_val)
	{
	    current_min_val = val;
	    for (i = 0; i < dim; ++i)
	        current_min_val_coords[i] = coords[i];
	}
	return current_min_val;
}		/*end compute_min*/

LOCAL	void	initialize_extreme_values(
	EXTREME_VALUES	*val)
{
	int k, i;
	val->rho_max = -HUGE_VAL; val->rho_min = HUGE_VAL;
	val->p_max   = -HUGE_VAL;   val->p_min = HUGE_VAL;
	val->e_max   = -HUGE_VAL;   val->e_min = HUGE_VAL;
	val->ke_max  = -HUGE_VAL;  val->ke_min = HUGE_VAL;
	val->E_max   = -HUGE_VAL;   val->E_min = HUGE_VAL;
	for (k = 0; k < 3; ++k)
	{
	    val->m_max[k] = -HUGE_VAL; val->m_min[k] = HUGE_VAL;
	    val->v_max[k] = -HUGE_VAL; val->v_min[k] = HUGE_VAL;
	}
	for (i = 0; i < 3; ++i)
	{
	    val->coords_rho_max[i] = -HUGE_VAL;
	    val->coords_rho_min[i] =  HUGE_VAL;
	    val->coords_p_max[i]   = -HUGE_VAL;
	    val->coords_p_min[i]   =  HUGE_VAL;
	    val->coords_e_max[i]   = -HUGE_VAL;
	    val->coords_e_min[i]   =  HUGE_VAL;
	    val->coords_ke_max[i]  = -HUGE_VAL;
	    val->coords_ke_min[i]  =  HUGE_VAL;
	    val->coords_E_max[i]   = -HUGE_VAL;
	    val->coords_E_min[i]   =  HUGE_VAL;
	    for (k = 0; k < 3; ++k)
	    {
	        val->coords_m_max[k][i] = -HUGE_VAL;
		val->coords_m_min[k][i] =  HUGE_VAL;
	        val->coords_v_max[k][i] = -HUGE_VAL;
		val->coords_v_min[k][i] =  HUGE_VAL;
	    }
	}
}		/*end initialize_extreme_values*/

EXPORT void check_for_consistent_tri_states(
	INTERFACE *intfc)
{
	TRI                *tri;
	SURFACE            *s;
	POINT              *p[3];
	HYPER_SURF         *hs;
	HYPER_SURF_ELEMENT *hse;
	Locstate           sl[3], sr[3];
	int                i;

	(void) next_tri(intfc,NULL,NULL);
	while (next_tri(intfc,&tri,&s))
	{
	    hs = Hyper_surf(s);
	    hse = Hyper_surf_element(tri);
	    for (i = 0; i < 3; ++i)
	    {
		p[i] = Point_of_tri(tri)[i];
	        slsr(p[i],hse,hs,sl+i,sr+i);
	    }
	    switch (consistent_params_in_tri_lin_comb(sl[0],sl[1],sl[2]))
	    {
	    case PARAMS_ALL_OBSTACLE_STATES:
	    case PARAMS_CONSISTENT:
		break;
	    case PARAMS_INCONSISTENT:
	    case OBSTACLE_STATE_FOUND:
		screen("ERROR in check_for_consistent_tri_states(), "
		       "inconsistent params on negative side of triangle\n");
		print_tri(tri,s->interface);
		print_tri_states(tri,hs);
		print_hypersurface(hs);
		print_interface(hs->interface);
		gview_plot_interface("LEFT_SIDE_INCONSISTENT_PARAMS",
				     hs->interface);
		clean_up(ERROR);
		break;
	    }
	    switch (consistent_params_in_tri_lin_comb(sr[0],sr[1],sr[2]))
	    {
	    case PARAMS_ALL_OBSTACLE_STATES:
	    case PARAMS_CONSISTENT:
		break;
	    case PARAMS_INCONSISTENT:
	    case OBSTACLE_STATE_FOUND:
		screen("ERROR in check_for_consistent_tri_states(), "
		       "inconsistent params on positive side of triangle\n");
		print_tri(tri,s->interface);
		print_tri_states(tri,hs);
		print_hypersurface(hs);
		print_interface(hs->interface);
		gview_plot_interface("RIGHT_SIDE_INCONSISTENT_PARAMS",
				     hs->interface);
		clean_up(ERROR);
		break;
	    }
	}
}		/*end check_for_consistent_tri_states*/
