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
*				griemann.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the Riemann solvers for gas dynamics.
*
*	riemann_solution() calculates the solution of a Riemann
*	problem along a given ray in the space-time plane. It is
*	called by g_riemann() in ghyp/oned/g1dhyp.c and by npt_w_speed().
*
*	find_mid_state(), with its subroutine mass_flux(), calculates
*	the middle state in the solution of a Riemann problem.
*/

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	boolean	Godunov_iter_conv(double,double,double,double,double,double,double);
LOCAL	boolean	secant_press_conv(double,double,double,double);
LOCAL	boolean    secant_vel_conv(double,double,double,double,double,double);
LOCAL	boolean	bernoulli_random_variable(void);
LOCAL	boolean    bisection_find_mid_state(Locstate,Locstate,double,double*,double*,
	                                 double,double*,double*,double*,double*,
					 RIEMANN_SOLVER_WAVE_TYPE*,
					 RIEMANN_SOLVER_WAVE_TYPE*,
					 double,double);
LOCAL	boolean	mid_state_is_vacuum(Locstate,Locstate,double,Locstate,
	                            int*,double*,int);
LOCAL	boolean	godunov_find_mid_state(Locstate,Locstate,double,double*,double*,
	                               double,double*,double*,double*,double*,
	                               RIEMANN_SOLVER_WAVE_TYPE*,
				       RIEMANN_SOLVER_WAVE_TYPE*,
				       double,double);
LOCAL	int	sample_rsoln(double,Locstate,Locstate,Locstate,Locstate,Locstate,
	                     Locstate,double*,double,double,double,double,double,
	                     double,RIEMANN_SOLVER_WAVE_TYPE,
			     RIEMANN_SOLVER_WAVE_TYPE,int,size_t);
LOCAL	boolean	secant_find_mid_state(Locstate,Locstate,double,double*,double*,
	                              double*,double*,double*,double*,
				      RIEMANN_SOLVER_WAVE_TYPE*,
				      RIEMANN_SOLVER_WAVE_TYPE*,
	                              double,double,double,double,double,double);
LOCAL	void	print_Riemann_wave_curves(Locstate,Locstate,double,double,double);

#if defined(PHASE_CODE) && defined(TWOD)
LOCAL	void	find_phase_mid_state(Locstate,Locstate,double*,double*,double*,
	                             double*,double*,double*,double*,double*,
	                             double*,int,double*,double*,double*,
	                             int,RIEMANN_SOLVER_WAVE_TYPE*,
				     RIEMANN_SOLVER_WAVE_TYPE*);
LOCAL	void	get_sample_phase_state(double,Locstate,Locstate,Locstate,int,
	                               Locstate,size_t,double,double,double,
	                               double*,int,RIEMANN_SOLVER_WAVE_TYPE);
LOCAL	void	sample_comp_ending_in_raref(double,Locstate,Locstate,Locstate,
	                                    Locstate,Locstate,Locstate,double*,
	                                    int,double,double,double,size_t,int);
LOCAL	void	sample_comp_ending_in_shock(double,Locstate,Locstate,Locstate,
	                                    Locstate,Locstate,double,double*,int,
	                                    double,double,double,double,size_t,int);
LOCAL	void	sample_raref_state(double,Locstate,Locstate,Locstate,size_t,
	                           int,int,double,double,double*);
LOCAL	void	sample_shock_state(Locstate,Locstate,int,Locstate,double,double,
	                           double,double*,double,size_t,int);
LOCAL	void	sample_split_raref_state(double,Locstate,Locstate,Locstate,
	                                 Locstate,size_t,int,int,
	                                 double,double,double*);
LOCAL	void	sample_split_shock(double,Locstate,Locstate,Locstate,Locstate,
	                           double*,double,double,double,double,size_t,int);
LOCAL	void	set_vel_and_pr_across_raref(Locstate,int,double);
LOCAL	void	set_vel_and_pr_across_shock(Locstate,int,Locstate,
	                                    double,double,double);
#endif /* defined(PHASE_CODE) && defined(TWOD) */

#if defined(DEBUG_GRIEMANN)
LOCAL	void	left_find_mid_state(const char*,double,double,double,double,double,
	                            double,RIEMANN_SOLVER_WAVE_TYPE,
				    RIEMANN_SOLVER_WAVE_TYPE);

LOCAL	boolean debug_riem_sol = NO;
LOCAL	boolean	debug_find_mid_state = NO;

EXPORT void set_debug_riem_sol(
	boolean	y_or_n)
{
	debug_riem_sol = y_or_n;
}		/*end set_debug_riem_sol*/


EXPORT void set_debug_find_mid_state(
	boolean	y_or_n)
{
	debug_find_mid_state = y_or_n;
}		/*end set_debug_find_mid_state*/
#endif /* defined(DEBUG_GRIEMANN) */


/* Aug 20 2002 : Myoung-Nyoun : tolerance for find_mid_state : Start */
LOCAL  double tolerance=1.e-6;
/* Aug 20 2002 : Myoung-Nyoun : tolerance for find_mid_state : End */


/*
*			riemann_solution():
*
*	Determines the solution of a Riemann problem along a given
*	ray, specified by the slope "sample_speed" in the space-time
*	plane.  The normal nx,ny defines the orientation of the
*	distinguished direction along which the Riemann problem is
*	posed; it points from left to right.
*
*/

EXPORT void riemann_solution(
	double		sample_speed,
	const double	*dir,
	Locstate	sl,
	Locstate	sr,
	Locstate	ans,
	int		state_type)
{
	double		vl[SMAXD] ,vr[SMAXD];	/* velocities */
	double		vans;			/* normal velocity */
	double		spdans;			/* maximum wave speed of ans */
	int		l_or_r;
	int		i, dim;
	static Locstate Tsl = NULL,Tsr = NULL;

#if defined(DEBUG_GRIEMANN)
	debug_print("riem_sol","Entered riemann_solution()\n");
	if (debugging("riem_sol"))
	    debug_riem_sol = YES;
	else if (debug_riem_sol == YES)
	    (void) printf("Entered riemann_solution()\n");

	if (debug_riem_sol == YES)
	{
	    verbose_print_state("sl",sl);
	    verbose_print_state("sr",sr);
	    print_state_type("state_type = ",state_type);
	}
#endif /* defined(DEBUG_GRIEMANN) */

	if (is_obstacle_state(sl) && is_obstacle_state(sr))
	{
	    g_obstacle_state(ans,g_sizest());
#if defined(DEBUG_GRIEMANN)
	    debug_print("riem_sol","Left riemann_solution()\n");
	    if ((debug_riem_sol == YES) && (!debugging("riem_sol")))
	        (void) printf("Left riemann_solution()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	    return;
	}
	if (Tsl == NULL) 
	{
	    if (!is_obstacle_state(sl))
	    {
	        (*Params(sl)->_alloc_state)(&Tsl,Params(sl)->sizest);
	        (*Params(sl)->_alloc_state)(&Tsr,Params(sl)->sizest);
	    }
	    else
	    {
	        (*Params(sr)->_alloc_state)(&Tsl,Params(sr)->sizest);
	        (*Params(sr)->_alloc_state)(&Tsr,Params(sr)->sizest);
	    }
	}
	set_state(Tsl,EGAS_STATE,sl);
	set_state(Tsr,EGAS_STATE,sr);
	dim = Params(sl)->dim;

	for (i = 0; i < dim; ++i)
	{
	    vr[i] = Vel(Tsr)[i];        vl[i] = Vel(Tsl)[i];
	    Vel(Tsr)[i] = 0.0;	        Vel(Tsl)[i] = 0.0;
	    Vel(ans)[i] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
	    Vel(Tsr)[0] += dir[i] * vr[i];
	    Vel(Tsl)[0] += dir[i] * vl[i];
	}

	l_or_r = onedrsoln(sample_speed,Tsl,Tsr,ans,&spdans,EGAS_STATE);

	vans = Vel(ans)[0];
	switch(l_or_r)
	{
	case LEFT_FAMILY:
	    for (i = 0; i < dim; ++i)
	        Vel(ans)[i] = (vans - Vel(Tsl)[0])*dir[i] + vl[i];
	    break;

	case RIGHT_FAMILY:
	    for (i = 0; i < dim; ++i)
	        Vel(ans)[i] = (vans - Vel(Tsr)[0])*dir[i] + vr[i];
	    break;
	}

#if defined(DEBUG_GRIEMANN)
	if (debug_riem_sol == YES)
	{
	    (void) printf("EGAS_STATE form of answer in riemann_solution()\n");
	    verbose_print_state("EGas answer state",ans);
	}
#endif /* defined(DEBUG_GRIEMANN) */

	set_state(ans,state_type,ans);

#if defined(DEBUG_GRIEMANN)
	if (debug_riem_sol == YES)
	    verbose_print_state("answer state",ans);
	debug_print("riem_sol","Left riemann_solution()\n");
	if ((debug_riem_sol == YES) && (!debugging("riem_sol")))
	    (void) printf("Left riemann_solution()\n");
#endif /* defined(DEBUG_GRIEMANN) */
}		/*end riemann_solution*/


EXPORT int onedrsoln(
	double		sample_speed,
	Locstate	sl,
	Locstate	sr,
	Locstate	ans,
	double		*spdnew,
	int		stype)
{
	RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
	int             l_or_r = ERROR;
	double           pml, pmr, uml, umr, mr, ml; /* midstate quantities */
	static Locstate Tsl = NULL, Tsr = NULL, mid = NULL;
	static size_t    sizest = 0;

#if defined(DEBUG_GRIEMANN)
	debug_print("riem_sol","Entered onedrsoln\n");
	if (debugging("riem_sol"))
	    debug_riem_sol = YES;
	else if (debug_riem_sol == YES)
	    (void) printf("Entered onedrsoln()\n");

	if (debug_riem_sol == YES)
	{
	    verbose_print_state("sl",sl);
	    verbose_print_state("sr",sr);
	    print_state_type("state_type = ",stype);
	}
#endif /* defined(DEBUG_GRIEMANN) */

#if defined(COMBUSTION_CODE)
	if (Composition_type(sl) != PURE_NON_REACTIVE ||
	    Composition_type(sr) != PURE_NON_REACTIVE)
	{
 	    if (Composition_type(sl) != THINFLAME ||
	        Composition_type(sr) != THINFLAME)
	    {
	        l_or_r = combust_onedrsoln(sample_speed,sl,sr,ans,spdnew,stype);
#if defined(DEBUG_GRIEMANN)
	        debug_print("riem_sol","Left onedrsoln\n");
	        if ((debug_riem_sol == YES) && (!debugging("riem_sol")))
	            (void) printf("Left onedrsoln()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        return l_or_r;
	    }
	}
#endif /* defined(COMBUSTION_CODE) */

	if ((is_obstacle_state(sl)) && (is_obstacle_state(sr)))
	{
	    g_obstacle_state(ans,g_sizest());
	    *spdnew = 0.0;
	    return l_or_r;
	}

	if (mid == NULL) 
	{
	    Gas_param *params = (is_obstacle_state(sl)) ? Params(sr):Params(sl);

	    sizest = params->sizest;
	    (*params->_alloc_state)(&Tsl,sizeof(VGas));
	    (*params->_alloc_state)(&Tsr,sizeof(VGas));
	    (*params->_alloc_state)(&mid,sizest);
	}


	set_state_for_find_mid_state(Tsl,sl);
	set_state_for_find_mid_state(Tsr,sr);

#if defined(DEBUG_GRIEMANN)
	if (debug_riem_sol == YES)
	{
	    verbose_print_state("Tsl",Tsl);
	    verbose_print_state("Tsr",Tsr);
	}
#endif /* defined(DEBUG_GRIEMANN) */

	/*
	 * Check if the mid_state is vacuum. If it is TRUE then we don't need
	 * to use the godunov iteration. 
	 */

	if (mid_state_is_vacuum(Tsl,Tsr,sample_speed,ans,&l_or_r,
	                        spdnew,stype) == YES)
	    return l_or_r;

	    /*   if mid_state is not vacuum  */

	if (find_mid_state(Tsl,Tsr,0.0,&pml,&pmr,&uml,&umr,
	                       &ml,&mr,&l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in onedrsoln(), find_mid_state() failed\n");
	    clean_up(ERROR);
	}

#if defined(PHASE_CODE) && defined(TWOD)
	if (multiphase_eos(Eos(Tsl)) && multiphase_eos(Eos(Tsr)))
	{
	        /*intersection points for phase transtion problems */
	    double rpi[2], lpi[2], rui[2], lui[2], rri[2], lri[2];
	    int nri, nli;

#if defined(DEBUG_GRIEMANN)
	    if (debug_riem_sol == YES)
	        (void) printf("Multiphase case\n");
#endif /* defined(DEBUG_GRIEMANN) */

	    if (!intrsct_wv_crv_wth_phs_bdry(Tsl,Tsr,rpi,rui,rri,&nri,lpi,
						lui,lri,&nli,l_wave,r_wave))
	    {
	        Wave_curve(Tsl) = Wave_curve(Tsr) = NULL;
	        l_or_r = sample_rsoln(sample_speed,sl,Tsl,sr,Tsr,mid,ans,
	                              spdnew,pml,pmr,uml,umr,ml,mr,
	                              l_wave,r_wave,stype,sizest);
	    }
	    else
	    {
	        if ((state_type(Tsl) != VGAS_STATE) ||
	            (state_type(Tsr) != VGAS_STATE))
	        {
	            screen("ERROR in onedrsoln(), Phase code ");
	            screen("called without VGAS_STATES\n");
	            clean_up(ERROR);
	        }
	        Wave_curve(Tsl) = Wave_curve(Tsr) = NULL;
	        if (debugging("ph_riem"))
	        {
	            (void) printf("The initial guess for a solution is\n");
	            (void) printf("pml = %g, pmr = %g, uml = %g, umr = %g, ",
		                  pml,pmr,uml,umr);
		    print_rsoln_wave("l_wave = ",l_wave,", ");
		    print_rsoln_wave("r_wave = ",r_wave,"\n");
	        }
	        find_phase_mid_state(Tsl,Tsr,&pml,&pmr,&uml,&umr,&ml,&mr,
	                             rpi,rui,rri,nri,lpi,lui,lri,nli,
	                             &l_wave,&r_wave);
	        
	        l_or_r = sample_rsoln(sample_speed,sl,Tsl,sr,Tsr,mid,ans,
	                              spdnew,pml,pmr,uml,umr,ml,mr,
	                              l_wave,r_wave,stype,sizest);
	    }
	}
	else
#endif /* defined(PHASE_CODE) && defined(TWOD) */
	l_or_r = sample_rsoln(sample_speed,sl,Tsl,sr,Tsr,mid,ans,spdnew,
	                      pml,pmr,uml,umr,ml,mr,l_wave,r_wave,
	                      stype,sizest);

#if defined(DEBUG_GRIEMANN)
	debug_print("riem_sol","Left onedrsoln\n");
	if ((debug_riem_sol == YES) && (!debugging("riem_sol")))
	    (void) printf("Left onedrsoln()\n");
#endif /* defined(DEBUG_GRIEMANN) */

	return l_or_r;
}	    /*end onedrsoln*/


/* 	
*				sample_rsoln()
*
*	This function does what its name suggests: it samples 
*	(i.e., evaluates) a Riemann solution, which is determined 
*	by the input data
*
*        pml		left mid state pressure
*        pmr		right mid state pressure
*        uml 		left mid state velocity
*        umr 		right mid state velocity
*        ml 		left wave mass flux
*        mr 		right wave mass flux
*        l_wave 	left wave type
*        r_wave 	right wave type
*        sl, Tsl 	left states (Tsl in TGAS form)
*        sr, Tsr   	right states (Tsr in TGAS form)
*
*        sample_speed  	slope of ray on which to sample the solution
*/

LOCAL	int sample_rsoln(
	double		         sample_speed,
	Locstate	         sl,
	Locstate	         Tsl,
	Locstate	         sr,
	Locstate	         Tsr,
	Locstate	         mid,
	Locstate	         ans,
	double		         *spdnew,
	double		         pml,
	double		         pmr,
	double		         uml,
	double		         umr,
	double		         ml,
	double		         mr,
	RIEMANN_SOLVER_WAVE_TYPE l_wave,
	RIEMANN_SOLVER_WAVE_TYPE r_wave,
	int		         state_type,
	size_t		         sizest)
{
	double		rr,rl;		/* density */
	double 		ur,ul;		/* normal velocities */
	double		cl,cr, cm;	/* sound speeds */
	double		u;
	int		l_or_r = ERROR;

	/* Aug 20 2002 : Myoung-Nyoun : Start */
	int             i;
	int             loop=1;
	/* Aug 20 2002 : Myoung-Nyoun : End */

#if defined(DEBUG_GRIEMANN)
	debug_print("riem_sol","Entered sample_rsoln\n");
	if (debugging("riem_sol"))
	    debug_riem_sol = YES;
	else if (debug_riem_sol == YES)
	    (void) printf("Entered sample_rsoln()\n");
#endif /* defined(DEBUG_GRIEMANN) */

	rr = Dens(Tsr);
	ur = vel(0,Tsr);

	rl = Dens(Tsl);
	ul = vel(0,Tsl);
#if defined(DEBUG_GRIEMANN)
	if (debug_riem_sol == YES)
	    (void) printf("sample_speed = %g, uml = %g, umr = %g\n",
			  sample_speed,uml,umr);
#endif /* defined(DEBUG_GRIEMANN) */

/* Aug 20 2002 : Myoung-Nyoun : loop header : Start */
      for(i=0;loop;i++)
      {
        loop=0;
/* Aug 20 2002 : Myoung-Nyoun : loop header : End */
	  
	if (sample_speed <= uml)		/* On left side of contact */
	{
#if defined(DEBUG_GRIEMANN)
	    if (debug_riem_sol == YES)
	        (void) printf("Left side of contact\n");
#endif /* defined(DEBUG_GRIEMANN) */
	    l_or_r = LEFT_FAMILY;
#if defined(PHASE_CODE) && defined(TWOD)
	    if ((state_type(Tsl) == VGAS_STATE) &&
	        Wave_curve(Tsl) != NULL && Wave_curve(Tsl)->special)
	    {
	        get_sample_phase_state(sample_speed,Tsl,sl,ans,state_type,
	                               mid,sizest,ml,uml,pml,spdnew,
	                               l_or_r,l_wave);
	        return l_or_r;
	    }
#endif /* defined(PHASE_CODE) && defined(TWOD) */
	    if (l_wave == SHOCK) /* Left shock */
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debug_riem_sol == YES)
	            (void) printf("l_wave == SHOCK\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        u = ul - ml/rl;
	        if (sample_speed <= u) /* Left of left shock */
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debug_riem_sol == YES) 
	                (void) printf("left of left shock\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            ft_assign(ans,sl,sizest);
	            *spdnew = fabs(ul) + sound_speed(Tsl);
	        }
	        else    /* Right of left shock */
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debug_riem_sol == YES) 
	                (void) printf("right of left shock\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            Dens(ans) = ml/(uml - u); 
	            Set_params(ans,Tsl);
	            set_type_of_state(ans,state_type);
                    /* Set partial density for solution */
                    /* Mass fraction is constant except across contact */
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                        int   k;
                        if(Params(Tsl)->n_comps != 1)
                        {
                            for(k = 0; k < Params(Tsl)->n_comps; k++)
                                pdens(ans)[k] = (pdens(Tsl)[k]/Dens(Tsl))*Dens(ans);
                        }
                    }
	            switch(state_type)
	            {
	            case TGAS_STATE:
	                Vel(ans)[0] = uml;
	                Press(ans) = pml;
	                break;
	            case GAS_STATE:
	                Mom(ans)[0] = Dens(ans)*uml;
	                Energy(ans) = Dens(ans)*(0.5*sqr(uml) +
	                              specific_internal_energy(Tsl) + 
	                              0.5*(pressure(Tsl) + pml)*
	                              (1.0/rl - 1.0/Dens(ans)));
	                break;
	            case EGAS_STATE:
	                Vel(ans)[0] = uml;
	                Energy(ans) = specific_internal_energy(Tsl) +
	                              0.5*(pressure(Tsl) + pml)*
	                              (1.0/rl - 1.0/Dens(ans));
	                break;
	            }
		    reset_gamma(ans);
	            *spdnew = fabs(uml)+sound_speed(ans);
	        }
	    }
	    else        /* Left rarefaction */
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debug_riem_sol == YES)
	            (void) printf("l_wave == RAREFACTION\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        cl = sound_speed(Tsl);
	        if (sample_speed <= ul - cl) /* Left of left fan */
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debug_riem_sol == YES)
	                (void) printf("left of left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            ft_assign(ans,sl,sizest);
	            *spdnew = fabs(ul) + cl;
	        }    
	        else
	        {
	            state_on_adiabat_with_pr(Tsl,pml,mid,state_type);
	            cm = sound_speed(mid);
	            if (sample_speed >= uml - cm)	/* Right of left fan */
	            {	
#if defined(DEBUG_GRIEMANN)
	                if (debug_riem_sol == YES) 
	                    (void) printf("right of left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	                ft_assign(ans,mid,sizest);
	                switch (state_type)
	                {
	                case TGAS_STATE:
	                case EGAS_STATE:
	                    Vel(ans)[0] = uml;
	                    break;
	                case GAS_STATE:
	                    Mom(ans)[0] = Dens(ans)*uml;
	                    Energy(ans) += 0.5*uml*Mom(ans)[0];
		            reset_gamma(ans);
				    
	                    break;
	                }
	                *spdnew = fabs(uml) + cm;
	            }
	            else                        /* In left fan */
	            {
#if defined(DEBUG_GRIEMANN)
	                if (debug_riem_sol == YES) 
	                    (void) printf("in left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	                (void) oned_state_in_rarefaction_fan(sample_speed,
	                                                     ul,Tsl,
	                                                     mid,ans,state_type,
	                                                     spdnew,
	                                                     LEFT_FAMILY);
	            }
	        }
	    }
	}
	else if (sample_speed >= umr) /* On right side of contact */
	{
#if defined(DEBUG_GRIEMANN)
	    if (debug_riem_sol == YES)
	        (void) printf("Right side of contact\n");
#endif /* defined(DEBUG_GRIEMANN) */
	    l_or_r = RIGHT_FAMILY;
#if defined(PHASE_CODE) && defined(TWOD)
	    if ((state_type(Tsr) == VGAS_STATE) &&
	        Wave_curve(Tsr) != NULL && Wave_curve(Tsr)->special)
	    {
	        get_sample_phase_state(sample_speed,Tsr,sr,ans,state_type,
	                               mid,sizest,mr,umr,pmr,
	                               spdnew,l_or_r,r_wave);
	        return l_or_r;
	    }
#endif /* defined(PHASE_CODE) && defined(TWOD) */
	    if (r_wave == SHOCK)		/* Right shock wave */
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debug_riem_sol == YES)
	            (void) printf("r_wave == SHOCK\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        u = ur + mr/rr;
	        if (sample_speed >= u) /* Right of right shock */
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debug_riem_sol == YES) 
	                (void) printf("right of right shock\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            ft_assign(ans,sr,sizest);
	            *spdnew = fabs(ur) + sound_speed(Tsr);
	        }
	        else    /* Left of right shock */
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debug_riem_sol == YES) 
	                (void) printf("left of right shock\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            Dens(ans) = -mr/(umr - u);
	            Set_params(ans,Tsr);
	            set_type_of_state(ans,state_type);
                    /* Set partial density for solution */
                    /* Mass fraction is constant except across contact */
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                        int   k;
                        if(Params(Tsr)->n_comps != 1)
                        {
                            for(k = 0; k < Params(Tsr)->n_comps; k++)
                                pdens(ans)[k] = (pdens(Tsr)[k]/Dens(Tsr))*Dens(ans);
                        }
                    }
	            switch(state_type)
	            {
	            case TGAS_STATE:
	                Vel(ans)[0] =  umr;
	                Press(ans) = pmr;
	                break;
	            case GAS_STATE:
	                Mom(ans)[0] =  Dens(ans)*umr;
	                Energy(ans) = Dens(ans)*(0.5*sqr(umr) +
	                              specific_internal_energy(Tsr) + 
	                              0.5*(pressure(Tsr) + pmr)*
	                              (1.0/rr - 1.0/Dens(ans)));
	                break;
	            case EGAS_STATE:
	                Vel(ans)[0] =  umr;
	                Energy(ans) = specific_internal_energy(Tsr) +
	                              0.5*(pressure(Tsr) + pmr)*
	                              (1.0/rr - 1.0/Dens(ans));
	                break;
	            }
		    reset_gamma(ans);
	            *spdnew = fabs(umr)+sound_speed(ans);
	        }
	    }
	    else            /* Right rarefaction */
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debug_riem_sol == YES)
	            (void) printf("r_wave == RAREFACTION\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        cr = sound_speed(Tsr);
	        if (sample_speed >= ur + cr) /* Right of right fan */
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debug_riem_sol == YES) 
	                (void) printf("right of right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            ft_assign(ans,sr,sizest);
	            *spdnew = fabs(ur) + cr;
	        }
	        else
	        {
	            state_on_adiabat_with_pr(Tsr,pmr,mid,state_type);
	            cm = sound_speed(mid);
	            if (sample_speed <= umr + cm)	/* Left of right fan */
	            {
#if defined(DEBUG_GRIEMANN)
	                if (debug_riem_sol == YES) 
	                    (void) printf("left of right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	                ft_assign(ans,mid,sizest);
	                switch(state_type)
	                {
	                case TGAS_STATE:
	                case EGAS_STATE:
	                    Vel(ans)[0] = umr;
	                    break;
	                case GAS_STATE:
	                    Mom(ans)[0] = Dens(ans)*umr;
	                    Energy(ans) += 0.5*umr*Mom(ans)[0];
	                    break;
	                }
			reset_gamma(ans);
	                *spdnew = fabs(umr) + cm;
	            }
	            else                        /* In right fan */
	            {
#if defined(DEBUG_GRIEMANN)
	                if (debug_riem_sol == YES) 
	                    (void) printf("in right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	                (void) oned_state_in_rarefaction_fan(sample_speed,
	                                                     ur,Tsr,
	                                                     mid,ans,state_type,
	                                                     spdnew,
	                                                     RIGHT_FAMILY);
	            }
	        }
	    }
	}
	else	/*Inside a vacuum*/
	{
	     if (mid_state_is_vacuum(Tsl,Tsr,sample_speed,ans,&l_or_r,
	                             spdnew,state_type) != YES)
	     {
	         screen("ERROR in sample_rsoln(), midstate should be a "
			"vacuum but is not!! : %d times\n", i+1);
		 screen("\tsample_speed = %"FFMT" uml = %"FFMT" umr = %"FFMT"\n"
			,sample_speed,uml,umr);

		 /* Aug 20 2002 : Myoung-Nyoun : use a finer tolerance : Start */
		 if(i < 1)
		 {
		     loop=1;
		     tolerance=1.e-10;
		     find_mid_state(Tsl,Tsr,0.0,&pml,&pmr,&uml,&umr,&ml,&mr,
				    &l_wave,&r_wave);
		     tolerance=1.e-6;
		     continue;
		 }
		 /* Aug 20 2002 : Myoung-Nyoun : use a finer tolerance : End */
		 
		 /* TMP */
		 screen("ERROR in sample_rsoln(), final error step\n");
		 {
		     double      ul = Vel(Tsl)[0];
		     double      ur = Vel(Tsr)[0];
		     double	pl = pressure(Tsl);
		     double	pr = pressure(Tsr);
		     double	il = acoustic_impedance(Tsl);
		     double	ir = acoustic_impedance(Tsr);
		     double	pbar = 0.5*(pl + pr);
		 
		     double	pstar = (il*pr + ir*pl + ir*il*(ul - ur))/(ir + il);
		     double	min_p = max(Min_pressure(Tsl),Min_pressure(Tsr));
		     double	vlmax = ul - riemann_wave_curve(Tsl,min_p);
		     double	vrmin = ur + riemann_wave_curve(Tsr,min_p);
		 
		     printf("either pstar = %"FFMT" > 0.1*pbar = %"FFMT"\n",pstar,0.1*pbar);
		     printf("or     vlmax = %"FFMT" > vrmin = %"FFMT"\n",vlmax,vrmin);
		 
		     printf("ul = %"FFMT", ur = %"FFMT"\n",ul,ur);
		     printf("pl = %"FFMT", pr = %"FFMT"\n",pl,pr);
		     printf("il = %"FFMT", ir = %"FFMT"\n",il,ir);
		     printf("pbar = %"FFMT", min_p = %"FFMT"\n",pbar,min_p);
		 }
		 printf("sample_speed = %"FFMT"\n",sample_speed);
		 printf("         uml = %"FFMT"\n",uml);
		 printf("         umr = %"FFMT"\n",umr);
		 (void) printf("sample_speed - uml = %"FFMT", "
			       "sample_speed - umr = %"FFMT"\n",
			       sample_speed-uml,sample_speed-umr);
		 verbose_print_state("Tsl : left state",Tsl);
		 verbose_print_state("Tsr : right state",Tsr);
		 
		 clean_up(ERROR);
	     }
	}

/* Aug 20 2002 : Myoung-Nyoun : loop tail : Start */
      }
/* Aug 20 2002 : Myoung-Nyoun : loop tail : End */

#if defined(DEBUG_GRIEMANN)
	if (debug_riem_sol == YES)
	{
	    (void) printf("Answer from sample_rsoln()\n");
	    (void) printf("l_or_r = %d (%s)\n",l_or_r,(l_or_r == LEFT_FAMILY) ?
	                "LEFT_FAMILY" : "RIGHT_FAMILY");
	    verbose_print_state("oned answer state",ans);
	}
	debug_print("riem_sol","Left sample_rsoln\n");
	if ((debug_riem_sol == YES) && (!debugging("riem_sol")))
	    (void) printf("Left sample_rsoln()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	return l_or_r;
}		/*end sample_rsoln*/



/*
*			find_mid_state():
*
*	Solves for the middle state, given the left and right states.
*	Uses the Van Leer interative scheme using a secant method
*	to solve for the mid state pressure.  See Colella and Glaz
*	"Efficient Solution Algorithms for the Riemann Problem
*	for Real Gases"  Jour. Comp. Phys. 59, No. 2, June 1985,
*	pp. 264-289.
*
*	For surface tension,  the pressure on the left side of a
*	contact is equal to pjump plus the pressure on the right side.
*/

#define Godunov_pressure(pl,vxl,ml,pr,vxr,mr)				   \
	((((vxl) - (vxr))*(mr)*(ml) + (pr)*(ml) + (pl)*(mr))/ ((ml) + (mr)))

LOCAL	int fms_algorithm = SECANT_RS;

EXPORT boolean find_mid_state(
	Locstate	         Tsl,
	Locstate	         Tsr,
	double		         pjump,
	double		         *ppstarl,
	double		         *ppstarr,
	double		         *pustarl,
	double		         *pustarr,
	double		         *pml,
	double		         *pmr,
	RIEMANN_SOLVER_WAVE_TYPE *pl_wave,
	RIEMANN_SOLVER_WAVE_TYPE *pr_wave)
{
	double	    vlmax, vrmin;
	double	    pr, pl;
	double	    p_start, p_min;
	double	    eps_u,eps_p; /* convergence limits for ustar, pstar */
	double       pmax, pmin;
	boolean	    status;
	const double SECANT_EPS = tolerance; /* TOLERANCE */

#if defined(DEBUG_GRIEMANN)
	debug_print("find_mid_state","Entered find_mid_state()\n");
	if (debugging("find_mid_state"))
	    debug_find_mid_state = YES;
	else if (debug_find_mid_state == YES)
	    (void) printf("Entered find_mid_state()\n");

	if (debug_find_mid_state == YES)
	{
	    (void) printf("States into find_mid_state()\n");
	    verbose_print_state("left state",Tsl);
	    verbose_print_state("right_state",Tsr);
	    (void) printf("Pressure jump due to surface tension is = %"FFMT"\n",
			  pjump);
	}
#endif /* defined(DEBUG_GRIEMANN) */

	if (debugging("bad_state"))
	{
	    boolean badl, badr;

	    badl = is_bad_state(Tsl,YES,"find_mid_state");
	    badr = is_bad_state(Tsr,YES,"find_mid_state");
	    if (badl || badr)
	    {
	        screen("ERROR in find_mid_state bad state detected as input\n");
	        if (badl)
	        {
	            fprint_raw_gas_data(stdout,Tsl,current_interface()->dim);
	            verbose_print_state("bad state Tsl",Tsl);
	        }
		else
	            verbose_print_state("Tsl",Tsl);
	        if (badr)
	        {
	            fprint_raw_gas_data(stdout,Tsr,current_interface()->dim);
	            verbose_print_state("bad state Tsr",Tsr);
	        }
		else
	            verbose_print_state("Tsr",Tsr);
		(void) printf("pjump = %"FFMT"\n",pjump);
	        clean_up(ERROR);
	    }
	}

	if (is_obstacle_state(Tsl) || is_obstacle_state(Tsr)) 
	{
	    screen("ERROR in find_mid_state obstacle state detected\n");
	    clean_up(ERROR);
	}

#if defined(COMBUSTION_CODE)
	if (Composition_type(Tsl) != PURE_NON_REACTIVE ||
	    Composition_type(Tsr) != PURE_NON_REACTIVE)
	{
	    if((Composition_type(Tsl) != THINFLAME) ||
	       (Composition_type(Tsr) != THINFLAME))
	    {
	        status = combust_find_mid_state(Tsl,Tsr,ppstarl,ppstarr,pustarl,
	                                        pustarr,pml,pmr,pl_wave,pr_wave);
#if defined(DEBUG_GRIEMANN)
	        left_find_mid_state("find_mid_state",*ppstarl,*ppstarr,
	                            *pustarl,*pustarr,*pml,*pmr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	        return status;
	    }
	}
#endif /* defined(COMBUSTION_CODE) */

	/* Initial guess using the linearized Godunov equation */

	initialize_riemann_solver(Tsl,Tsr,&p_start,&p_min,
	                          SECANT_EPS,&eps_u,&eps_p,find_mid_state);

	/*Check for vacuum in solution*/
	pr = p_min;
	pl = floor_pressure(p_min+pjump,p_min);
	vrmin = vel(0,Tsr) + riemann_wave_curve(Tsr,pr);
	vlmax = vel(0,Tsl) - riemann_wave_curve(Tsl,pl);
#if defined(DEBUG_GRIEMANN)
	if (debug_find_mid_state == YES)
	    (void) printf("Vacuum check, p_min = %"FFMT", pr = %"FFMT", pl = %"FFMT"\n"
			  "vrmin = %"FFMT", vlmax = %"FFMT"\n",p_min,pr,pl,vrmin,vlmax);
#endif /* defined(DEBUG_GRIEMANN) */
	if (vlmax < vrmin)
	{
#if defined(DEBUG_GRIEMANN)
	    if (debug_find_mid_state == YES)
	        (void) printf("Vacuum formed\n");
#endif /* defined(DEBUG_GRIEMANN) */
	    /* wave curves do not intersect,  a vacuum is formed */
	    *ppstarl = pl;
	    *pustarl = vlmax;
 	    *pml = mass_flux(pl,Tsl);
	    *pl_wave = (pl >= pressure(Tsl)) ? SHOCK : RAREFACTION;

	    *ppstarr = pr;
	    *pustarr = vrmin;
 	    *pmr = mass_flux(pr,Tsr);
	    *pr_wave = (pr >= pressure(Tsr)) ? SHOCK : RAREFACTION;

#if defined(DEBUG_GRIEMANN)
	    left_find_mid_state("find_mid_state",*ppstarl,*ppstarr,
	                        *pustarl,*pustarr,*pml,*pmr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	    return FUNCTION_SUCCEEDED;
	}

	pmax = HUGE_VAL;
	pmin = p_min;
	switch (fms_algorithm)
	{
	case GODUNOV_RS:
	    *ppstarr = p_start;
	    *ppstarl = floor_pressure(p_start+pjump,p_min);
	    status = godunov_find_mid_state(Tsl,Tsr,pjump,ppstarl,ppstarr,
	                                    p_min,pustarl,pustarr,pml,pmr,
	                                    pl_wave,pr_wave,pmin,pmax);
	    break;

	case SECANT_RS:
	default:
	    status = secant_find_mid_state(Tsl,Tsr,pjump,ppstarl,ppstarr,
	                                   pustarl,pustarr,pml,pmr,pl_wave,
	                                   pr_wave,p_start,p_min,eps_u,eps_p,
					   pmin,pmax);
	    break;
	}

#if defined(DEBUG_GRIEMANN)
	left_find_mid_state("find_mid_state",*ppstarl,*ppstarr,
	                    *pustarl,*pustarr,*pml,*pmr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	return status;
}		/*end find_mid_state*/

#define	vacuum_produced(pstar,ul_mid,ur_mid)				\
	((pstar <= 2.0*p_min) /*TOLERANCE*/ && (ul_mid < ur_mid))


LOCAL	boolean	secant_press_conv(
	double	pstar,
	double	p_old,
	double	eps_p,
	double	meps)
{
	double	eps = eps_p*pstar;

	eps = max(eps,meps);
	return (fabs(pstar-p_old) <= eps) ? YES : NO;
}		/*end secant_press_conv*/

/*ARGSUSED*/
LOCAL	boolean secant_vel_conv(
	double ulm,
	double ulm_old,
	double urm,
	double urm_old,
	double eps_u,
	double meps)
{
	double	eps, eps_l, eps_r;

	eps = 0.5*(fabs(ulm)+fabs(urm))*eps_u;	eps = max(eps,meps);
	if (fabs(ulm-urm) <= eps)
	    return YES;

	eps_l = fabs(ulm)*eps_u;	eps_l = max(eps_l,meps);
	eps_r = fabs(urm)*eps_u;	eps_r = max(eps_r,meps);

	return ((fabs(ulm-ulm_old) <= eps_l) && (fabs(urm-urm_old) <= eps_r)) ?
	    YES : NO;
}		/*end secant_vel_conv*/

LOCAL boolean secant_find_mid_state(
	Locstate	         Tsl,	/* left TGas or VGas states */
	Locstate	         Tsr,	/* right TGas or VGas states */
	double		         pjump,	/* Pressure jump across contact, pl-pr*/
	double		         *ppstarl,	/* left middle pressure */
	double		         *ppstarr,	/* right middle pressure */
	double		         *pustarl,	/* left middle velocity */
	double		         *pustarr,	/* right middle velocity */
	double		         *pml,
	double		         *pmr,
	RIEMANN_SOLVER_WAVE_TYPE *pl_wave,
	RIEMANN_SOLVER_WAVE_TYPE *pr_wave, /* wave families */
	double		         p_start,
	double		         p_min,
	double		         eps_u,	/*convergence limits for ustar, pstar*/
	double		         eps_p,	/*convergence limits for ustar, pstar*/
	double                    pmin,  /*solution lies in range */
	double                    pmax)  /*pmin <= p <= pmax      */
{
	const double meps = 10.0*MACH_EPS;/*TOLERANCE*/
	int	    i;
	double	    pstarr, pstarl;
	double	    ml,	mr;
	double	    ul_mid, ur_mid, ul_mid_old, ur_mid_old, p_old, p_old2;
	double	    gul_mid, gur_mid;
	double	    vul, vur;
	double	    numer, denom;
	double	    pl, pr, vxl, vxr;
	double	    poj, plj;
	double	    du_p_start, du_pstar, du_ppstar;
	boolean	    status;
	const double MIN_DENOM = 1.0e-10; /* TOLERANCE */
	const int   NUM_SEC_ITER = 20;/* Max # secant iterations */

#if defined(DEBUG_GRIEMANN)
	debug_print("find_mid_state","Entered secant_find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */

	p_old = p_start;

	poj = floor_pressure(p_old+pjump,p_min);
 	ml = mass_flux(poj,Tsl);
	pl = pressure(Tsl);	
	vxl = vel(0,Tsl);		

 	mr = mass_flux(p_old,Tsr);
	pr = pressure(Tsr);
	vxr = vel(0,Tsr);

	plj = floor_pressure(pl-pjump,p_min);
	pstarr = Godunov_pressure(plj,vxl,ml,pr,vxr,mr);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (pstarr <= p_min)
	    pstarr = p_min;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	pstarl = floor_pressure(pstarr+pjump,p_min);
	ul_mid = vxl;
	if (fabs(ml) > MIN_DENOM)
	    ul_mid += (plj - pstarl)/(ml);
	ur_mid = vxr;
	if (fabs(mr) > MIN_DENOM)
	    ur_mid += (pstarr - pr)/(mr);

#if defined(DEBUG_GRIEMANN)
	if (debug_find_mid_state == YES)
	{
	    (void) printf("eps_p = %"FFMT", eps_u = %"FFMT", p_min = %"FFMT", p_start = %"FFMT"\n",
	                  eps_p,eps_u,p_min,p_start);
	    (void) printf("\tGodunov iter: pstarl = %"FFMT", pstarr = %"FFMT" "
	                  "ul = %"FFMT", ur = %"FFMT"\n",pstarl,pstarr,ul_mid,ur_mid);
	    (void) printf("\t\tml = %"FFMT", mr = %"FFMT"\n\t\t u(sl) = %"FFMT", u(sr) = %"FFMT"\n",
	                  ml,mr,vxl,vxr);
	    (void) printf("\t\tP(sl) = %"FFMT", P(sr) = %"FFMT"\n",pl,pr);
	}
#endif /* defined(DEBUG_GRIEMANN) */

	                /*   secant method  */

	for (i = 0; i < NUM_SEC_ITER; ++i)
	{
	    ur_mid_old = ur_mid;
	    ul_mid_old = ul_mid;
	    p_old2 = p_old;
	    p_old = pstarr;
	    poj = floor_pressure(p_old+pjump,p_min);
	    ur_mid = vxr + riemann_wave_curve(Tsr,p_old);
	    ul_mid = vxl - riemann_wave_curve(Tsl,poj);
	    if (ur_mid <= ul_mid)
	    {
		if (pmin < p_old)
		    pmin = p_old;
	    }
	    else
	    {
		if (p_old < pmax)
		    pmax = p_old;
	    }
	    denom = (ur_mid - ur_mid_old) - (ul_mid - ul_mid_old);
	    numer = (ur_mid - ul_mid)*(p_old - p_old2);
	    if (fabs(denom) >= MIN_DENOM &&
	        fabs(numer) > 0.5*eps_p*p_old*fabs(denom))
	    {
	        pstarr = p_old - numer/denom;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	        if (pstarr <= p_min)
	            pstarr = p_min;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        pstarl = floor_pressure(pstarr+pjump,p_min);
	    }

#if defined(DEBUG_GRIEMANN)
	    if (debug_find_mid_state == YES)
	    {
	        (void) printf("after %d step: pstarl = %"FFMT" pstarr = %"FFMT", ul = %"FFMT", ur = %"FFMT"\n",
	                      i,pstarl,pstarr,ul_mid,ur_mid);
	    }
#endif /* defined(DEBUG_GRIEMANN) */

	    if (secant_press_conv(pstarr,p_old,eps_p,meps))
	    {
	        if (secant_vel_conv(ul_mid,ul_mid_old,ur_mid,
	                            ur_mid_old,eps_u,meps))
	        {
	            *ppstarl = pstarl;
	            *ppstarr = pstarr;
	            *pustarl = *pustarr = 0.5*(ur_mid + ul_mid);
 		    *pml = mass_flux(pstarl,Tsl);
 		    *pmr = mass_flux(pstarr,Tsr);
	            *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	            *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	            debug_print("find_mid_state","Left secant_find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            return FUNCTION_SUCCEEDED;
	        }
	        else if (vacuum_produced(pstarr,ul_mid,ur_mid))
	        {
	            /* Vacuum produced */
	            *ppstarr = pstarr = p_min;
	            *ppstarl = pstarl = floor_pressure(p_min+pjump,p_min);
	            *pustarl = ul_mid;
	            *pustarr = ur_mid;
 		    *pml = mass_flux(pstarl,Tsl);
 		    *pmr = mass_flux(pstarr,Tsr);
	            *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	            *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	            debug_print("find_mid_state","Left secant_find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            return FUNCTION_SUCCEEDED;
	        }
	        if (i > 0)
	            break;
	    }
	}

	/* Check for vacuum solution */
	poj = floor_pressure(p_min+pjump,p_min);
	vur = vxr + riemann_wave_curve(Tsr,p_min);
	vul = vxl - riemann_wave_curve(Tsl,poj);
	if (vul < vur)
	{
	    *ppstarr = pstarr = p_min;
	    *ppstarl = pstarl = floor_pressure(p_min+pjump,p_min);
	    *pustarl = ul_mid;
	    *pustarr = ur_mid;
 	    *pml = mass_flux(pstarl,Tsl);
 	    *pmr = mass_flux(pstarr,Tsr);
	    *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	    *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	    debug_print("find_mid_state","Left secant_find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	    return FUNCTION_SUCCEEDED;
	}

	poj = floor_pressure(pstarr+pjump,p_min);
	ur_mid = vxr + riemann_wave_curve(Tsr,pstarr);
	ul_mid = vxl - riemann_wave_curve(Tsl,poj);
	du_pstar = fabs(ur_mid - ul_mid);

	if (debugging("secant_fail"))
	{
	    double eps_l = fabs(ul_mid)*eps_u;
	    double eps_r = fabs(ur_mid)*eps_u;
	    double eps_m = 0.5*(fabs(ul_mid)+fabs(ur_mid))*eps_u;
	    eps_l = max(eps_l,meps);
	    eps_r = max(eps_r,meps);
	    eps_m = max(eps_m,meps);

	    (void) printf("WARNING in secant_find_mid_state(), "
	                  "secant method did not converge\n");
	    (void) printf("Number of iterations used = %d\n",i);
	    (void) printf("Convergence tolerances, eps_u = %"FFMT", meps = %"FFMT"\n",
	                  eps_u,meps);
	    (void) printf("p_old = %"FFMT", pstarl = %"FFMT" pstarr = %"FFMT", "
	                  "|p_old - pstarr| = %"FFMT"\n",
	                  p_old,pstarl,pstarr,fabs(p_old - pstarr));
	    (void) printf("pmin = %"FFMT", pmax = %"FFMT"\n\t",pmin,pmax);
	    (void) printf("ul_mid = %"FFMT", ul_mid_old = %"FFMT"\n\t",ul_mid,ul_mid_old);
	    (void) printf("|ul_mid - ul_mid_old| = %"FFMT", left tolerance = %"FFMT"\n",
	                   fabs(ul_mid - ul_mid_old),eps_l);
	    (void) printf("ur_mid = %"FFMT", ur_mid_old = %"FFMT"\n\t",ur_mid,ur_mid_old);
	    (void) printf("|ur_mid - ur_mid_old| = %"FFMT", right tolerance = %"FFMT"\n",
	                  fabs(ur_mid - ur_mid_old),eps_r);

	    (void) printf("pstarl = %"FFMT" pstarr = %"FFMT", "
	                  "ur_mid = %"FFMT", ul_mid = %"FFMT"\n\t",
	                  pstarl,pstarr,ur_mid,ul_mid);
	    (void) printf("|ur_mid - ul_mid| = %"FFMT", mid tolerance = %"FFMT"\n",
	                  du_pstar,eps_m);
	    print_Riemann_wave_curves(Tsl,Tsr,pjump,p_min,pstarr);
	}

	poj = floor_pressure(p_start+pjump,p_min);
	du_p_start = fabs(vxr - vxl + riemann_wave_curve(Tsr,p_start) +
	                              riemann_wave_curve(Tsl,poj));
	*ppstarr = (du_pstar < du_p_start) ? pstarr : p_start;
	*ppstarl = floor_pressure(*ppstarr+pjump,p_min);
	status = godunov_find_mid_state(Tsl,Tsr,pjump,ppstarl,ppstarr,
	                                p_min,pustarl,pustarr,pml,pmr,
	                                pl_wave,pr_wave,pmin,pmax);

	gur_mid = vxr + riemann_wave_curve(Tsr,*ppstarr);
	gul_mid = vxl - riemann_wave_curve(Tsl,*ppstarl);
	du_ppstar = fabs(gur_mid - gul_mid);

	if (debugging("secant_fail"))
	{
	    (void) printf("godunov_find_mid_state() started with %s = %"FFMT"\n",
	                  (du_pstar < du_p_start) ? "pstarr" : "p_start",
	                  (du_pstar < du_p_start) ? pstarr : p_start);

	    (void) printf("godunov_find_mid_state() returned\n");
	    (void) printf("*ppstarl = %"FFMT", *ppstarr = %"FFMT", "
	                  "gur_mid = %"FFMT", gul_mid = %"FFMT"\n\t",
	                  *ppstarr,*ppstarr,gur_mid,gul_mid);
	    (void) printf("|gur_mid - gul_mid| = %"FFMT"\n",du_ppstar);
	    (void) printf("pmin = %"FFMT", pmax = %"FFMT"\n",pmin,pmax);
	}

	/*
	*	In cases of convergence close to tolerance,
	*	godunov_find_mid_state() does not always improve the
	*	solution.  Use the better of the two results.
	*/

	if (du_pstar < du_ppstar)
	{
	    *ppstarr = pstarr;
	    *ppstarl = pstarl;
	    *pustarr = *pustarl = 0.5*(ur_mid + ul_mid);
 	    *pml = mass_flux(pstarl,Tsl);
 	    *pmr = mass_flux(pstarr,Tsr);
	    *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	    *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
	}

#if defined(DEBUG_GRIEMANN)
	debug_print("find_mid_state","Left secant_find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	return status;
}		/*end secant_find_mid_state*/


/*
*			godunov_find_mid_state():
*
*	Solves for the middle state, given the left and right states.
*	The solution uses Godunov's method to solve the nonlinear
*	functional equations.  The iteration method is due to Chorin,
*	J. Comp. Phys. 22 (1976) 517-533.
*/

LOCAL	boolean	Godunov_iter_conv(
	double	ml,
	double	mlo,
	double	mr,
	double	mro,
	double	eps_l,
	double	eps_r,
	double	meps)
{
	eps_l = max(ml*eps_l,meps);
	eps_r = max(mr*eps_r,meps);
	return ((fabs(ml-mlo) <= eps_l) && (fabs(mr-mro) <= eps_r)) ? YES : NO;
}		/*end Godunov_iter_conv*/


#define Godunov_velocity(pl,vxl,ml,pr,vxr,mr)				\
	(((pl) - (pr) + (mr)*(vxr) + (ml)*(vxl)) /  ((ml) + (mr)))

LOCAL boolean godunov_find_mid_state(
	Locstate	         Tsl,	/* left TGas or VGas states */
	Locstate	         Tsr,	/* right TGas or VGas states */
	double		         pjump,	/* Pressure jump across contact, pl-pr*/
	double		         *ppstarl,/* left middle pressure */
	double		         *ppstarr,/* right middle pressure */
	double		         p_min,
	double		         *pustarl,/* left middle velocity */
	double		         *pustarr,/* right middle velocity */
	double		         *pml,
	double		         *pmr,
	RIEMANN_SOLVER_WAVE_TYPE *pl_wave,		/* left wave family */
	RIEMANN_SOLVER_WAVE_TYPE *pr_wave,		/* right wave family */
	double                    pmin, /* solution lies in range */
	double                    pmax) /* pmin <= p <= pmax      */
{
	int	    i,j;
	double	    pstarr = *ppstarr;
	double	    pstarl = *ppstarl;
	double	    poj, ur_mid, ul_mid;
	double	    plj;
	double	    ml, mr;
	double	    ml_old, mr_old, p_tilde, alpha, p_old;
	double	    eps_l,eps_r;	 	/* Eos dependent convergence
	                                           limits for ml and mr */
	double	    pl = pressure(Tsl), pr = pressure(Tsr);
	double	    vxl = vel(0,Tsl), vxr = vel(0,Tsr);
	const double meps = 10.0*MACH_EPS;/*TOLERANCE*/
	const double R_EPS = tolerance; /*TOLERANCE, Governs convergence */
	const double R_EPS_FACTOR = 10.0; /* TOLERANCE, 
				                 * Controls oscillation in 
						 * ml and mr
						 */
	const int   NUM_GODUNOV_ITER = 20; /*TOLERANCE,
						    Max # Godunov iterations */

#if defined(DEBUG_GRIEMANN)
	debug_print("find_mid_state","Entered godunov_find_mid_state()\n");
	if ((debug_find_mid_state == YES) && (!debugging("find_mid_state")))
	    (void) printf("Entered godunov_find_mid_state()\n");
	if (debug_find_mid_state == YES)
	{
	    (void) printf("States into godunov_find_mid_state()\n");
	    verbose_print_state("left state",Tsl);
	    verbose_print_state("right_state",Tsr);
	}
#endif /* defined(DEBUG_GRIEMANN) */


/* Initial guess using the linearized Godunov equation */

	eps_l = eps_for_Godunov(Tsl,pstarl,R_EPS);
	eps_r = eps_for_Godunov(Tsr,pstarr,R_EPS);

 	ml = mass_flux(pstarl,Tsl);
 	mr = mass_flux(pstarr,Tsr);
#if defined(DEBUG_GRIEMANN)
	if (debug_find_mid_state == YES)
	{
	    (void) printf("\tpstarl = %"FFMT", pstarr = %"FFMT" ml = %"FFMT", mr = %"FFMT"\n",
	                  pstarl,pstarr,ml,mr);
	}
#endif /* defined(DEBUG_GRIEMANN) */

/*
*	Iterate at least twice to avoid spurious convergence when Press(Tsr)
*	= Press(Tsl).
*/

	plj = floor_pressure(pl-pjump,p_min);
	pstarr = Godunov_pressure(plj,vxl,ml,pr,vxr,mr);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (pstarr <= p_min)
	    pstarr = p_min;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	pstarl = floor_pressure(pstarr+pjump,p_min);
 	ml = mass_flux(pstarl,Tsl);
 	mr = mass_flux(pstarr,Tsr);

/* Iterative loop using Godunov's method, as adapted by Chorin */

	for (i = 0; i < NUM_GODUNOV_ITER; ++i)
	{
	    p_old = pstarr;
	    ml_old = ml;
	    mr_old = mr;
	    pstarr = Godunov_pressure(plj,vxl,ml,pr,vxr,mr);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (pstarr <= p_min)
	        pstarr = p_min;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    pstarl = floor_pressure(pstarr+pjump,p_min);
 	    ml = mass_flux(pstarl,Tsl);
 	    mr = mass_flux(pstarr,Tsr);
	    ul_mid = vxl - riemann_wave_curve(Tsl,pstarl);
	    ur_mid = vxr + riemann_wave_curve(Tsr,pstarr);
	    if (ur_mid < ul_mid)
	    {
		if (pmin < pstarr)
		    pmin = pstarr;
	    }
	    else if (ul_mid < ur_mid)
	    {
		if (pstarr < pmax)
		    pmax = pstarr;
	    }
	    else
	    {
	        *ppstarr = pstarr;
	        *ppstarl = pstarl;
	        *pustarl = *pustarr = 0.5*(ul_mid + ur_mid);
	        *pml = ml;
	        *pmr = mr;
	        *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	        *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	        left_find_mid_state("godunov_find_mid_state",*ppstarl,*ppstarr,
	                            *pustarl,*pustarr,ml,mr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	        return FUNCTION_SUCCEEDED;
	    }
#if defined(DEBUG_GRIEMANN)
	    if (debug_find_mid_state == YES)
	    {
	        (void) printf("\tpstarl = %"FFMT", pstarr = %"FFMT" ml = %"FFMT", mr = %"FFMT"\n",
	                      pstarl,pstarr,ml,mr);
	    }
#endif /* defined(DEBUG_GRIEMANN) */
	    if (Godunov_iter_conv(ml,ml_old,mr,mr_old,eps_l,eps_r,meps))
	    {
	        *ppstarr = pstarr;
	        *ppstarl = pstarl;
	        *pustarl = *pustarr = Godunov_velocity(plj,vxl,ml,pr,vxr,mr);
	        *pml = ml;
	        *pmr = mr;
	        *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	        *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	        left_find_mid_state("godunov_find_mid_state",*ppstarl,*ppstarr,
	                            *pustarl,*pustarr,ml,mr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	        return FUNCTION_SUCCEEDED;
	    }
	}

	/* Check for small oscillations in ml and mr */

	if ((fabs(pstarr - p_old) < max(R_EPS*p_old,meps)) && 
	    Godunov_iter_conv(ml,ml_old,mr,mr_old,R_EPS_FACTOR*eps_l,
	                      R_EPS_FACTOR*eps_r,meps))
	{
	    *ppstarr = pstarr;
	    *ppstarl = pstarl;
	    *pustarl = *pustarr = Godunov_velocity(plj,vxl,ml,pr,vxr,mr);
	    *pml = ml;
	    *pmr = mr;
	    *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	    *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	    left_find_mid_state("godunov_find_mid_state",*ppstarl,*ppstarr,
	                        *pustarl,*pustarr,ml,mr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	    return FUNCTION_SUCCEEDED;
	}

	/* Check for vacuum solution */
	poj = floor_pressure(p_min+pjump,p_min);
	ur_mid = vxr + riemann_wave_curve(Tsr,p_min);
	ul_mid = vxl - riemann_wave_curve(Tsl,poj);
	if (ul_mid < ur_mid)
	{
	    /* Mid state is a vacuum */
	    *ppstarr = pstarr = p_min;
	    *ppstarl = pstarl = floor_pressure(pstarr+pjump,p_min);
	    *pustarl = ul_mid;
	    *pustarr = ur_mid;
 	    *pml = mass_flux(pstarl,Tsl);
 	    *pmr = mass_flux(pstarr,Tsr);
	    *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	    *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	    debug_print("find_mid_state","Left godunov_find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	    return FUNCTION_SUCCEEDED;
	}

/* If first iteration loop fails, correct guess by method due to Chorin. */

	pstarl = floor_pressure(pstarr+pjump,p_min);
	eps_l = eps_for_Godunov(Tsl,pstarl,R_EPS);
	eps_r = eps_for_Godunov(Tsr,pstarr,R_EPS);

	alpha = 1.0;
	for (j = 0; j < 16; ++j)
	{
	    alpha *= 0.5;
#if defined(DEBUG_GRIEMANN)
            if (debug_find_mid_state == YES)
            {
                (void) printf("Correct guess by method due to Chorin\n");
		(void) printf("alpha = %g\n",alpha);
            }
#endif /* defined(DEBUG_GRIEMANN) */
	    for (i = 0; i < NUM_GODUNOV_ITER; ++i) 
	    {
	        p_old = pstarr;
	        ml_old = ml;
	        mr_old = mr;
	        p_tilde = Godunov_pressure(plj,vxl,ml,pr,vxr,mr);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	        if (p_tilde <= p_min)
	            p_tilde = p_min;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        pstarr = alpha*p_tilde + (1. - alpha)*(pstarr);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	        if (pstarr <= p_min)
	            pstarr = p_min;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        pstarl = floor_pressure(pstarr+pjump,p_min);
 	        ml = mass_flux(pstarl,Tsl);
 	        mr = mass_flux(pstarr,Tsr);
	        ul_mid = vxl - riemann_wave_curve(Tsl,pstarl);
	        ur_mid = vxr + riemann_wave_curve(Tsr,pstarr);
#if defined(DEBUG_GRIEMANN)
	        if (debug_find_mid_state == YES)
	        {
	            (void) printf("\tpstarl = %g, pstarr = %g\n",
		    		pstarl,pstarr);
		    (void) printf("\t    ml = %g,     mr = %g\n",ml,mr);
	            (void) printf("\tul_mid = %g, ur_mid = %g\n\n",
	                      	ul_mid,ur_mid);
	        }
#endif /* defined(DEBUG_GRIEMANN) */
	        if (ur_mid < ul_mid)
	        {
		    if (pmin < pstarr)
		        pmin = pstarr;
	        }
	        else if (ul_mid < ur_mid)
	        {
		    if (pstarr < pmax)
		        pmax = pstarr;
	        }
	        else
	        {
	            *ppstarr = pstarr;
	            *ppstarl = pstarl;
	            *pustarl = *pustarr = 0.5*(ul_mid + ur_mid);
	            *pml = ml;
	            *pmr = mr;
	            *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	            *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	            left_find_mid_state("godunov_find_mid_state",
			                *ppstarl,*ppstarr,*pustarl,*pustarr,
					ml,mr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	            return FUNCTION_SUCCEEDED;
	        }
	        if (Godunov_iter_conv(ml,ml_old,mr,mr_old,eps_l,eps_r,meps))
	        {
	            *ppstarr = pstarr;
	            *ppstarl = pstarl;
	            *pustarl = *pustarr=Godunov_velocity(plj,vxl,ml,pr,vxr,mr);
	            *pml = ml;
	            *pmr = mr;
	            *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	            *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	            if (debug_find_mid_state == YES)
	            {
	                (void) printf("alpha in find_mid_state = %"FFMT"\n",alpha);
	            }
	            left_find_mid_state("godunov_find_mid_state",
	                                *ppstarl,*ppstarr,*pustarl,*pustarr,
	                                ml,mr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	            return FUNCTION_SUCCEEDED;
	        }
	    }
	
	/* Check for small oscillations in ml and mr */

	    if ((fabs(pstarr - p_old) < max(R_EPS*p_old,meps)) && 
	        Godunov_iter_conv(ml,ml_old,mr,mr_old,R_EPS_FACTOR*eps_l,
	                          R_EPS_FACTOR*eps_r,meps))
	    {
	        *ppstarr = pstarr;
	        *ppstarl = pstarl;
	        *pustarl = *pustarr = Godunov_velocity(plj,vxl,ml,pr,vxr,mr);
	        *pml = ml;
	        *pmr = mr;
	        *pl_wave = (pstarl >= pl) ? SHOCK : RAREFACTION;
	        *pr_wave = (pstarr >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	        left_find_mid_state("godunov_find_mid_state",
	                           *ppstarl,*ppstarr,*pustarl,*pustarr,
	                            ml,mr,*pl_wave,*pr_wave);
#endif /* defined(DEBUG_GRIEMANN) */
	        return FUNCTION_SUCCEEDED;
	    }
	}
	if (debugging("secant_fail"))
	{
	    (void) printf("WARNING in godunov_find_mid_state(), "
	                  "Godunov method did not converge\n");
	    (void) printf("Convergence tolerances, eps_l = %"FFMT", eps_r = %"FFMT" "
			  "meps = %"FFMT"\n",
	                  eps_l,eps_r,meps);
	    (void) printf("p_old = %"FFMT", pstarl = %"FFMT" pstarr = %"FFMT", "
	                  "|p_old - pstarr| = %"FFMT"\n",
	                  p_old,pstarl,pstarr,fabs(p_old - pstarr));
	    (void) printf("pmin = %"FFMT", pmax = %"FFMT"\n",pmin,pmax);

	    (void) printf("ml = %"FFMT", ml_old = %"FFMT"\n",ml,ml_old);
	    (void) printf("|ml - ml_old| = %"FFMT", left tolerance = %"FFMT"\n",
	                   fabs(ml - ml_old),eps_l);

	    (void) printf("mr = %"FFMT", mr_old = %"FFMT"\n",mr,mr_old);
	    (void) printf("|mr - mr_old| = %"FFMT", right tolerance = %"FFMT"\n",
	                  fabs(mr - mr_old),eps_r);

	    (void) printf("pstarl = %"FFMT" pstarr = %"FFMT", "
	                  "ur_mid = %"FFMT", ul_mid = %"FFMT"\n\t",
	                  pstarl,pstarr,ur_mid,ul_mid);
	    (void) printf("|ur_mid - ul_mid| = %"FFMT"\n",fabs(ur_mid-ul_mid));
	    print_Riemann_wave_curves(Tsl,Tsr,pjump,p_min,pstarr);
	}
	return bisection_find_mid_state(Tsl,Tsr,pjump,ppstarl,ppstarr,
                                        p_min,pustarl,pustarr,pml,pmr,
					pl_wave,pr_wave,pmin,pmax);

}		/*end godunov_find_mid_state*/

LOCAL boolean bisection_find_mid_state(
        Locstate                 Tsl,   /* left TGas or VGas states */
        Locstate                 Tsr,   /* right TGas or VGas states */
        double                    pjump, /* Pressure jump across contact,pl-pr */
        double                    *ppstarl,      /* left middle pressure */
        double                    *ppstarr,      /* right middle pressure */
        double                    p_min,
        double                    *pustarl,      /* left middle velocity */
        double                    *pustarr,      /* right middle velocity */
        double                    *pml,
        double                    *pmr,
	RIEMANN_SOLVER_WAVE_TYPE *pl_wave,		/* left wave family */
	RIEMANN_SOLVER_WAVE_TYPE *pr_wave,		/* right wave family */
	double                    pmin,  /* solution lies in the range */
	double                    pmax)  /* pmin <= p <= pmax     */
{
        double           ml, mr, ml_old, mr_old;
        double           ur_mid,ul_mid;
        double           pl = pressure(Tsl), pr = pressure(Tsr);
        double           vxl = vel(0,Tsl), vxr = vel(0,Tsr);
        double           poj,p_mid;
        double           eps_l,eps_r;
        const double     R_EPS = tolerance;     /*TOLERANCE Governs convergence */
        const double     meps = 10.0*MACH_EPS;/*TOLERANCE*/
        int             i;

#if defined(DEBUG_GRIEMANN)
        debug_print("find_mid_state","Entered bisection_find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */

        eps_l = eps_for_Godunov(Tsl,*ppstarl,R_EPS);
        eps_r = eps_for_Godunov(Tsr,*ppstarr,R_EPS);
        poj = floor_pressure(p_min+pjump,p_min);

        vxl = vel(0,Tsl);
        vxr = vel(0,Tsr);
        pl = pressure(Tsl);
        pr = pressure(Tsr);

        mr = mass_flux(p_min,Tsr);
        ml = mass_flux(poj,Tsl);
        ul_mid = vxl - (poj - pl)/ml;
        ur_mid = vxr + (p_min - pr)/mr;
        if (ur_mid >= ul_mid)              /* vacuum case */
        {
            *ppstarl = p_min;
            *ppstarr = poj;
            *pustarl = ul_mid;
            *pustarr = ur_mid;
            *pml = ml;
            *pmr = mr;
            *pl_wave = (*ppstarl >= pl) ? SHOCK : RAREFACTION;
            *pr_wave = (*ppstarr >= pr) ? SHOCK : RAREFACTION;
            return YES;
        }
        else if (pmax == HUGE_VAL)
        {
            pmax = max(pl,pr);
            if (vxr < vxl)
            {
                pmax = max(pl,pr);
                for (i = 0; i < 10; ++i)
                {
                    pmin = pmax;
                    pmax *= 2.0;
                    poj = pmax + pjump;
                    mr = mass_flux(pmax,Tsr);
                    ml = mass_flux(poj,Tsl);
                    ur_mid = vxr + (pmax - pr)/mr;
                    ul_mid = vxl - (poj - pl)/ml;
                    if (ul_mid > ur_mid)
                        break;
                }
                if (ul_mid < ur_mid)
                {
	            (void) printf("WARNING in bisection_find_mid_state(), "
	                          "can't find upper bound for pressure\n");
	            *ppstarr = pmax;
	            *ppstarl = poj;
	            *pustarl = *pustarr = 0.5*(ul_mid + ur_mid);
	            *pml = ml;
	            *pmr = mr;
	            *pl_wave = (poj >= pl) ? SHOCK : RAREFACTION;
	            *pr_wave = (pmax >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	            (void) printf("States into find_mid_state()\n");
	            poj = floor_pressure(pmax+pjump,p_min);
	            verbose_print_state("left state",Tsl);
	            verbose_print_state("right_state",Tsr);
	            (void) printf("\nConvergence tolerances\n");
	            (void) printf("\tur_mid %"FFMT", ul_mid %"FFMT"\n",
	                          vxr + riemann_wave_curve(Tsr,pmax),
	                          vxl - riemann_wave_curve(Tsl,poj));
	            (void) printf("End Convergence toloerances\n");
	            (void) printf("\nReturned Data\n");
	            (void) printf("\tpstarl = %"FFMT" pstarr = %"FFMT" "
				  "ustarl = %"FFMT" ustarr = %"FFMT" ml = %"FFMT" mr = %"FFMT"\n",
	                            poj,pmax,*pustarl,*pustarr,ml,mr);
	            print_rsoln_wave("l_wave = ",*pl_wave,", ");
	            print_rsoln_wave("r_wave = ",*pr_wave,"\n");
	            debug_print("find_mid_state","Left find_mid_state()\n");
	            if ((debug_find_mid_state == YES) &&
			    (!debugging("find_mid_state")))
	                (void) printf("Left find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            return FUNCTION_FAILED;
                }
            }
        }
        mr = mass_flux(pmin,Tsr);
        ml = mass_flux(pmin+pjump,Tsl);
        i = 0;
	do
        {
            p_mid = 0.5*(pmin + pmax);
            poj = p_mid + pjump;
            ml_old = ml;
            mr_old = mr;
            mr = mass_flux(p_mid,Tsr);
            ml = mass_flux(poj,Tsl);
            ur_mid = vxr + (p_mid - pr)/mr;
            ul_mid = vxl - (poj - pl)/ml;
            if (ur_mid < ul_mid)
                pmin = p_mid;
            else if (ur_mid > ul_mid)
                pmax = p_mid;
            else
            {
                *ppstarr = p_mid;
                *ppstarl = poj;
                *pustarl = *pustarr = 0.5*(ur_mid + ul_mid);
                *pml = ml;
                *pmr = mr;
                *pl_wave = (poj >= pl) ? SHOCK : RAREFACTION;
                *pr_wave = (p_mid >= pr) ? SHOCK : RAREFACTION;
                return YES;
            }
            if (++i > 100)
	    {
	        (void) printf("WARNING in bisection_find_mid_state(), "
	                      "bisection method did not converge\n");
                *ppstarr = p_mid;
                *ppstarl = poj;
                *pustarl = *pustarr = 0.5*(ur_mid + ul_mid);
                *pml = ml;
                *pmr = mr;
                *pl_wave = (poj >= pl) ? SHOCK : RAREFACTION;
                *pr_wave = (p_mid >= pr) ? SHOCK : RAREFACTION;
#if defined(DEBUG_GRIEMANN)
	        (void) printf("States into find_mid_state()\n");
	        poj = floor_pressure(p_mid+pjump,p_min);
	        verbose_print_state("left state",Tsl);
	        verbose_print_state("right_state",Tsr);
	        (void) printf("\nConvergence tolerances\n");
	        (void) printf("\tfabs(ml - ml_old) %"FFMT", eps_l*ml %"FFMT"\n",
			      fabs(ml-ml_old),R_EPS*ml);
	        (void) printf("\tfabs(mr - mr_old) %"FFMT", eps_l*mr %"FFMT"\n",
			      fabs(mr-mr_old),R_EPS*mr);
	        (void) printf("\tur_mid %"FFMT", ul_mid %"FFMT"\n",
	                      vxr + riemann_wave_curve(Tsr,p_mid),
	                      vxl - riemann_wave_curve(Tsl,poj));
	        (void) printf("End Convergence toloerances\n");
	        (void) printf("\nReturned Data\n");
	        (void) printf("\tpstarl = %"FFMT" pstarr = %"FFMT" "
			      "ustarl = %"FFMT" ustarr = %"FFMT" ml = %"FFMT" mr = %"FFMT"\n",
	                      poj,p_mid,*pustarl,*pustarr,ml,mr);
	        print_rsoln_wave("l_wave = ",*pl_wave,", ");
	        print_rsoln_wave("r_wave = ",*pr_wave,"\n");
	        debug_print("find_mid_state","Left find_mid_state()\n");
	        if ((debug_find_mid_state == YES) &&
			(!debugging("find_mid_state")))
	            (void) printf("Left find_mid_state()\n");
#endif /* defined(DEBUG_GRIEMANN) */
		return FUNCTION_FAILED;
	    }
        }
        while (!Godunov_iter_conv(ml,ml_old,mr,mr_old,eps_l,eps_r,meps));
        *ppstarr = p_mid;
        *ppstarl = poj;
        *pustarl = *pustarr = 0.5*(ur_mid + ul_mid);
        *pml = ml;
        *pmr = mr;
        *pl_wave = (poj >= pl) ? SHOCK : RAREFACTION;
        *pr_wave = (p_mid >= pr) ? SHOCK : RAREFACTION;
        left_find_mid_state("bisection_find_mid_state",*ppstarl,*ppstarr,
                            *pustarl,*pustarr,ml,mr,*pl_wave,*pr_wave);
        return YES;
}       	/* end bisection_find_mid_state */


#if defined(DEBUG_GRIEMANN)
LOCAL	void left_find_mid_state(
	const char               *name,
	double		         pstarl,
	double		         pstarr,
	double		         ustarl,
	double		         ustarr,
	double		         ml,
	double		         mr,
	RIEMANN_SOLVER_WAVE_TYPE l_wave,
	RIEMANN_SOLVER_WAVE_TYPE r_wave)
{
	if (debug_find_mid_state == YES) 
	{
	    (void) printf("Answer from %s()\n",name);
	    if (pstarl != pstarr)
	        (void) printf("pstarl = %"FFMT" pstarr = %"FFMT" ",pstarl,pstarr);
	    else
	        (void) printf("pstar = %"FFMT" ",pstarl);
	    if (ustarl != ustarr)
	        (void) printf("ustarl = %"FFMT" ustarr = %"FFMT" ",ustarl,ustarr);
	    else
	        (void) printf("ustar = %"FFMT" ",ustarl);
	    (void) printf("ml = %"FFMT" mr = %"FFMT"\n",ml,mr);
	    print_rsoln_wave("l_wave = ",l_wave,", ");
	    print_rsoln_wave("r_wave = ",r_wave,"\n");
	}
	debug_print("find_mid_state","Left %s()\n",name);
	if ((debug_find_mid_state == YES) && (!debugging("find_mid_state")))
	    (void) printf("Left %s()\n",name);
}		/*end left_find_mid_state*/
#endif /* defined(DEBUG_GRIEMANN) */


/*
*			mid_state_is_vacuum():
*
*	Solves the riemann problem when a vacuum occurs
*/

LOCAL	boolean	mid_state_is_vacuum(
	Locstate	Tsl,
	Locstate	Tsr,
	double		sample_speed,
	Locstate	ans,
	int		*side,
	double		*spdnew,
	int		state_type)
{
	double		 min_fan_speed, max_fan_speed;
	double		 vlmax, vrmin;
	double		 c, W;
	double		 ul, ur, sl, sr, ml, mr;
	double		 pl, pr, pbar, il, ir;
	double		 pstar;
	double		 min_p;
	static	Locstate vac_state = NULL;

#if defined(DEBUG_GRIEMANN)
	debug_print("vacuum","Entered mid_state_is_vacuum()\n");
#endif /* defined(DEBUG_GRIEMANN) */

	if (vac_state == NULL)
	    (*Params(Tsl)->_alloc_state)(&vac_state,Params(Tsl)->sizest);

	/*
	* In order to avoid integrating the wave curve down to vacuum,
	* first solve a linearized Riemann problem with linearized wave curves
	*		p = pr + ir*(u - ur)
	*		p = pl - i-*(u - u-).
	* If the wave curves are convex,  the pressure of this linear
	* system is a lower bound for the pressure of the nonlinear
	* system.  We declare that no vacuum is produced if the linearized
	* solution has pressure greater than a tolerance fraction of the
	* average pressures of the left and right states.
	*/

	ul = Vel(Tsl)[0];		ur = Vel(Tsr)[0];
	pl = pressure(Tsl);		pr = pressure(Tsr);
	il = acoustic_impedance(Tsl);	ir = acoustic_impedance(Tsr);
	pbar = 0.5*(pl + pr);

	pstar = (il*pr + ir*pl + ir*il*(ul - ur))/(ir + il);
#if defined(DEBUG_GRIEMANN)
	if (debugging("vacuum"))
	    (void) printf("pstar = %"FFMT", pl = %"FFMT", pr = %"FFMT"\n"
			  "pbar = %"FFMT"\n",pstar,pl,pr,pbar);
#endif /* defined(DEBUG_GRIEMANN) */

	min_p = max(Min_pressure(Tsl),Min_pressure(Tsr));
	vlmax = ul - riemann_wave_curve(Tsl,min_p);
	vrmin = ur + riemann_wave_curve(Tsr,min_p);
	min_fan_speed = vel(0,Tsl) - sound_speed(Tsl);
	max_fan_speed = vel(0,Tsr) + sound_speed(Tsr); 
#if defined(DEBUG_GRIEMANN)
	if (debugging("vacuum"))
	    (void) printf("min_p = %"FFMT", vlmax = %"FFMT", vrmin = %"FFMT"\n"
			  "min_fan_speed = %"FFMT", max_fan_speed = %"FFMT"\n",
			  min_p,vlmax,vrmin,min_fan_speed,max_fan_speed);
#endif /* defined(DEBUG_GRIEMANN) */

	if (vlmax > vrmin)
	{
#if defined(DEBUG_GRIEMANN)
	    debug_print("vacuum","Left mid_state_is_vacuum(), vlmax > vrmin, "
			   "answer = NO\n");
#endif /* defined(DEBUG_GRIEMANN) */
	    return NO;
	}
	if (pl < min_p)	/*Left shock - right rarefaction*/
	{
	    ml = mass_flux(min_p,Tsl);
	    sl = vlmax - ml/Dens(Tsl);
	    if (sample_speed <= sl)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsl\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        set_state(ans,state_type,Tsl);
	        *spdnew = fabs(min_fan_speed);
	    }
	    else if (sample_speed <= vlmax)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State behind left shock\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        state_behind_sound_wave(Tsl,ans,&c,&W,0.0,ml,vlmax,
	                                min_p,state_type,BACKWARD_SHOCK_WAVE,
	                                SHOCK,LEFT_FAMILY);
	        *spdnew = fabs(vlmax) + c;
	    }
	    else if (sample_speed <= vrmin)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from right\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        state_on_adiabat_with_pr(Tsr,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vrmin;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vrmin);
	    }
	    else if (sample_speed <= max_fan_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State inside right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        state_on_adiabat_with_pr(Tsr,min_p,vac_state,TGAS_STATE);
	        (void) oned_state_in_rarefaction_fan(sample_speed,ur,Tsr,
	                                             vac_state,ans,state_type,
	                                             spdnew,RIGHT_FAMILY);
	    }
	    else
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsr\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        set_state(ans,state_type,Tsr);
	        *spdnew = fabs(max_fan_speed);
	    }
	}
	else if (pr < min_p) /*Right shock - left rarefaction*/
	{
	    mr = mass_flux(min_p,Tsr);
	    sr = vrmin + mr/Dens(Tsr);
	    if (sr <= sample_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsr\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        set_state(ans,state_type,Tsr);
	        *spdnew = fabs(max_fan_speed);
	    }
	    else if (vrmin <= sample_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State behind right shock\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        state_behind_sound_wave(Tsr,ans,&c,&W,0.0,mr,vrmin,
	                                min_p,state_type,FORWARD_SHOCK_WAVE,
	                                SHOCK,RIGHT_FAMILY);
	        *spdnew = fabs(vlmax) + c;
	    }
	    else if (vlmax <= sample_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from left\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        state_on_adiabat_with_pr(Tsl,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vlmax;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vlmax);
	    }
	    else if (min_fan_speed <= sample_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State inside left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        state_on_adiabat_with_pr(Tsl,min_p,vac_state,TGAS_STATE);
	        (void) oned_state_in_rarefaction_fan(sample_speed,ul,Tsl,
	                                             vac_state,ans,state_type,
	                                             spdnew,LEFT_FAMILY);
	    }
	    else
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsl\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        set_state(ans,state_type,Tsl);
	        *spdnew = fabs(min_fan_speed);
	    }
	}
	else if (min_p > Min_pressure(Tsl))
	{
	    /*Left rarefaction, right rarefaction to vacuum*/
	    if (sample_speed <= min_fan_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsl\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        set_state(ans,state_type,Tsl);
	        *spdnew = fabs(min_fan_speed);
	    }
	    else if (sample_speed >= max_fan_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsr\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        set_state(ans,state_type,Tsr);
	        *spdnew = fabs(max_fan_speed);
	    }
	    else if (vrmin <= sample_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State inside right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        state_on_adiabat_with_pr(Tsr,min_p,vac_state,TGAS_STATE);
	        (void) oned_state_in_rarefaction_fan(sample_speed,ur,Tsr,
	                                             vac_state,ans,state_type,
	                                             spdnew,RIGHT_FAMILY);
	    }
	    else if (vlmax <= sample_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from right\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        state_on_adiabat_with_pr(Tsr,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vrmin;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vrmin);
	    }
	    else
	    {
	        ml = mass_flux(min_p,Tsl);
	        state_behind_sound_wave(Tsl,ans,&c,&W,0.0,ml,vlmax,
	                                min_p,state_type,BACKWARD_SHOCK_WAVE,
	                                SHOCK,LEFT_FAMILY);
	        if (sample_speed < (vlmax - c))
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debugging("vacuum"))
	                (void) printf("State inside left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            *side = LEFT_FAMILY;
	            set_state(vac_state,TGAS_STATE,ans);
	            (void) oned_state_in_rarefaction_fan(sample_speed,ul,Tsl,
	                                                 vac_state,ans,
	                                                 state_type,spdnew,
	                                                 LEFT_FAMILY);
	        }
	        else
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debugging("vacuum"))
	                (void) printf("State behind left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            *spdnew = fabs(vlmax) + c;
	        }
	    }
	}
	else if (min_p > Min_pressure(Tsr))
	{
	    /*Right rarefaction, left rarefaction to vacuum*/
	    if (sample_speed <= min_fan_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsl\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        set_state(ans,state_type,Tsl);
	        *spdnew = fabs(min_fan_speed);
	    }
	    else if (sample_speed >= max_fan_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsr\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        set_state(ans,state_type,Tsr);
	        *spdnew = fabs(max_fan_speed);
	    }
	    else if (sample_speed <= vlmax)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State inside left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        state_on_adiabat_with_pr(Tsl,min_p,vac_state,TGAS_STATE);
	        (void) oned_state_in_rarefaction_fan(sample_speed,
	                                ul,Tsl,vac_state,
	                                ans,state_type,spdnew,LEFT_FAMILY);
	    }
	    else if (sample_speed <= vrmin)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from left\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        state_on_adiabat_with_pr(Tsl,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vlmax;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vlmax);
	    }
	    else
	    {
	        mr = mass_flux(min_p,Tsr);
	        state_behind_sound_wave(Tsr,ans,&c,&W,0.0,mr,vrmin,
	                                min_p,state_type,FORWARD_SHOCK_WAVE,
	                                SHOCK,RIGHT_FAMILY);
	        if (sample_speed > (vrmin + c))
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debugging("vacuum"))
	                (void) printf("State inside right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            *side = RIGHT_FAMILY;
	            set_state(vac_state,TGAS_STATE,ans);
	            (void) oned_state_in_rarefaction_fan(sample_speed,
	                                                 ur,Tsr,vac_state,
	                                                 ans,state_type,
	                                                 spdnew,RIGHT_FAMILY);
	        }
	        else
	        {
#if defined(DEBUG_GRIEMANN)
	            if (debugging("vacuum"))
	                (void) printf("State behind right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	            *spdnew = fabs(vrmin) + c;
	        }
	    }
	}
	else
	{
	    /* Rarefaction to vacuum on both sides */
	    if (sample_speed <= min_fan_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsl\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        set_state(ans,state_type,Tsl);
	        *spdnew = fabs(min_fan_speed);
	    }
	    else if (sample_speed >= max_fan_speed)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State copied from Tsr\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        set_state(ans,state_type,Tsr);
	        *spdnew = fabs(max_fan_speed);
	    }
	    else if (sample_speed <= vlmax)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State inside left fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        state_on_adiabat_with_pr(Tsl,min_p,vac_state,TGAS_STATE);
	        (void) oned_state_in_rarefaction_fan(sample_speed,ul,Tsl,
	                                             vac_state,ans,state_type,
	                                             spdnew,LEFT_FAMILY);
	    }
	    else if (sample_speed >= vrmin)
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State inside right fan\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        state_on_adiabat_with_pr(Tsr,min_p,vac_state,TGAS_STATE);
	        (void) oned_state_in_rarefaction_fan(sample_speed,ur,Tsr,
	                                             vac_state,ans,state_type,
	                                             spdnew,RIGHT_FAMILY);
	    }
	    else if (Dens(Tsl) > Dens(Tsr))
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from left\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        state_on_adiabat_with_pr(Tsl,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vlmax;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vlmax);
	    }
	    else if (Dens(Tsr) > Dens(Tsl))
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from right\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        state_on_adiabat_with_pr(Tsr,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vrmin;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vrmin);
	    }
	    else if (bernoulli_random_variable())
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from left\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = LEFT_FAMILY;
	        state_on_adiabat_with_pr(Tsl,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vlmax;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vlmax);
	    }
	    else
	    {
#if defined(DEBUG_GRIEMANN)
	        if (debugging("vacuum"))
	            (void) printf("State at vacuum from right\n");
#endif /* defined(DEBUG_GRIEMANN) */
	        *side = RIGHT_FAMILY;
	        state_on_adiabat_with_pr(Tsr,min_p,ans,TGAS_STATE);
	        Vel(ans)[0] = vrmin;
	        set_state(ans,state_type,ans);
	        *spdnew = fabs(vrmin);
	    }
	}
#if defined(DEBUG_GRIEMANN)
	debug_print("vacuum","Left mid_state_is_vacuum(), answer = YES\n");
#endif /* defined(DEBUG_GRIEMANN) */
	return YES;
}		/*end mid_state_is_vacuum*/


/*	
*		 oned_state_in_rarefaction_fan(): 	
*
*	This function computes the state along the space time line in a
*	rarefaction fan with slope dx/dt = speed.
*	U0 is the component of the velocity of state0 normal to the
*	fan.
*
*/

EXPORT	boolean oned_state_in_rarefaction_fan(
	double	 spl_sp,
	double	 u0,
	Locstate Ts0,
	Locstate mid,
	Locstate ans,
	int	 typeans,
	double	 *spdnew,
	int	 w_family)
{
	double		vans, cans;
	double		v;
	double		w;
	double		sign;
	boolean		vacuum_detected;
	int		tmp_type = TGAS_STATE;

	sign = (w_family == RIGHT_FAMILY) ? 1.0 : -1.0;
	w = sign*(spl_sp - u0) - sound_speed(Ts0);

	cans = oned_fan_state(w,Ts0,mid,ans,tmp_type,&vacuum_detected);

	vans = spl_sp - sign*cans;
	if (vacuum_detected)
	{
	    v = u0 + sign*riemann_wave_curve(Ts0,Press(ans));
	    vans = (w_family == LEFT_FAMILY) ? min(v,vans) : max(v,vans);
	}
#if defined(COMBUSTION_CODE)
	if (Composition_type(Ts0) == ZND) React(ans) = react(Ts0);
#endif /* defined(COMBUSTION_CODE) */

	Vel(ans)[0] = vans;
	set_type_of_state(ans,TGAS_STATE);
	set_state(ans,typeans,ans);
	*spdnew = fabs(vans) + cans;
	return vacuum_detected;
}		/*end oned_state_in_rarefaction_fan*/

LOCAL	boolean bernoulli_random_variable(void)
{
	static boolean	first = YES;

	if (first == YES)
	{
	    long seed = 0x3fd1;
	    first = NO;
	    srand48(seed);
	}

	return (lrand48() % 2) ? YES : NO;
}		/*end bernoulli_random_variable*/

#if defined(PHASE_CODE) && defined(TWOD)

/*
*		find_phase_mid_state():
*
*	Finds the middle states for phase transition problems in which shock
*	splitting and complex rarefaction curves occur.
*
*/

LOCAL void find_phase_mid_state(
	Locstate	         l_state,
	Locstate	         r_state,
	double		         *pstarl,
	double		         *pstarr,
	double		         *ustarl,
	double		         *ustarr,
	double		         *pml,
	double		         *pmr,
	double		         *rpi,
	double		         *rui,
	double		         *rri,
	int		         nri,
	double		         *lpi,
	double		         *lui,
	double		         *lri,
	int		         nli,
	RIEMANN_SOLVER_WAVE_TYPE *l_wave,
	RIEMANN_SOLVER_WAVE_TYPE *r_wave)
{
	RIEMANN_SOLVER_WAVE_TYPE l_wtype, r_wtype;
	int		         i, nr, nl;
	double		         ml, mr, pl, pr, ul, ur;
	double		         lu[2], lp[2], lr[2], rr[2], ru[2], rp[2];
	static	WAVE_CURVE       *left_wave = NULL, *right_wave = NULL;

	debug_print("ph_riem","Entered find_phase_mid_state\n");
	if (debugging("ph_riem"))
	    (void) printf("%d %d\n",nri,nli);

	if (right_wave == NULL)
	{
	    Gas_param	*params = Params(l_state);
	    size_t sizest = params->sizest;

	    scalar(&left_wave,sizeof(WAVE_CURVE));
	    scalar(&right_wave,sizeof(WAVE_CURVE));
	    for (i = 0; i < MAX_NUM_WAVES; ++i)
	    {
	        (*params->_alloc_state)(&left_wave->st[i],sizest);
	        (*params->_alloc_state)(&right_wave->st[i],sizest);
	    }
	}

	/* Check if no special behavior */

	Wave_curve(l_state) = Wave_curve(r_state) = NULL;
	if (*l_wave == SHOCK)
	{
	    for (i = 0; i < nli; ++i)
	    {
	        if (lpi[i] < *pstarl)
	        {
	            /* Determine if shock splitting */
	            if (!is_retrograde_bndry(lri[i],l_state))
	                continue;
	            Wave_curve(l_state) = left_wave;
	        }
	    }
	}
	if (*l_wave == RAREFACTION)
	{
	    for (i = 0; i < nli; ++i)
	    {
	        if (lpi[i] > *pstarl) 
	        {
	            debug_print("ph_riem","Split rarefaction\n");
                    Wave_curve(l_state) = left_wave;
	        }
	    }
	}

	if (*r_wave == SHOCK)
	{
	    for (i = 0; i < nri; ++i)
	    {
	        if (rpi[i] < *pstarr)
	        {
	            /* Determine if shock spritting */
	            if (!is_retrograde_bndry(rri[i],r_state))
	                continue;
	            Wave_curve(r_state) = right_wave;
	        }
	    }
	}
	if (*r_wave == RAREFACTION)
	{
	    for (i = 0; i < nri; ++i)
	    {
	        if (rpi[i] > *pstarr) 
	        {
	            debug_print("ph_riem","Split rarefaction\n");
	            Wave_curve(r_state) = right_wave;
	        }
	    }
	}

	if (Wave_curve(r_state) == NULL && Wave_curve(l_state) == NULL)
	    return;


	if ( nri < 2 || nli < 2)
	{
	    l_wtype = (*l_wave == SHOCK) ? RAREFACTION : SHOCK;
	    r_wtype = (*r_wave == SHOCK) ? RAREFACTION : SHOCK;
	    if (intrsct_wv_crv_wth_phs_bdry(l_state,r_state,rp,ru,rr,&nr,lp,lu,
	                                    lr,&nl,l_wtype,r_wtype))
	    {
	        for (i = nli; i < nl; ++i)
	        {
	            lpi[i] = lp[i];
	            lri[i] = lr[i];
	            lui[i] = lu[i];
	        }
	        for (i = nri; i < nr; ++i)
	        {
	            rpi[i] = rp[i];
	            rri[i] = rr[i];
	            rui[i] = ru[i];
	        }
	        nri = nr + nri;
	        nli = nl + nli;
	    }
	}
	make_wave_crv(LEFT_FAMILY,l_state,lpi,lri,lui,nli,left_wave);
	make_wave_crv(RIGHT_FAMILY,r_state,rpi,rri,rui,nri,right_wave);

	(void) find_mid_state(l_state,r_state,0.0,&pl,&pr,&ul,&ur,&ml,&mr,
	                      &l_wtype,&r_wtype);
	*pstarr = pr;
	*pstarl = pl;
	*ustarl = ul;
	*ustarr = ul;
	*pml = ml;
	*pmr = mr;
	*l_wave = l_wtype;
	*r_wave = r_wtype;
}		/*end find_phase_mid_state*/


LOCAL void get_sample_phase_state(
	double		         sample_speed,
	Locstate	         Ts,
	Locstate	         state,
	Locstate	         ans,
	int		         state_type,
	Locstate	         mid,
	size_t		         sizest,
	double		         m,
	double		         um,
	double		         pm,
	double		         *spdnew,
	int		         l_or_r,
	RIEMANN_SOLVER_WAVE_TYPE wtype)
{
	WAVE_CURVE     *wave_cur;
	double           u, p, p1;
	double           sign, us, ms, ps, rhos, shockspd1, shockspd;
	double           press, rho;
	double           rstart, rend, rhoc, ustart, uend, pstart, pend;
	int             num_waves;
	int             i;
	static Locstate st0 = NULL, st1 = NULL;

	if (st1 == NULL)
	{
	    (*Params(Ts)->_alloc_state)(&st0, Params(Ts)->sizest);
	    (*Params(Ts)->_alloc_state)(&st1, Params(Ts)->sizest);
	}

	wave_cur = Wave_curve(Ts);
	num_waves = wave_cur->num_waves;



	if (l_or_r == LEFT_FAMILY)
	    sign = -1.0;
	else
	    sign = 1.0;

	rho = Dens(Ts);
	p = pressure(Ts);
	u = vel(0, Ts);

	if (wtype == SHOCK)
	{
	    /* In shock branch */

	    /* Go from state closest to um on out */

	    for (i = 1; i < num_waves; ++i)
	    {
	        if (wave_cur->w_type[num_waves - i] != SHOCK)
	            continue;

	        press = Press(wave_cur->st[num_waves - i]);
	        if (press > pm)
	            continue;

	        if ((wave_cur->w_type[num_waves - i - 1] != SHOCK) ||
		    (wave_cur->w_type[num_waves - i - 1] == SHOCK &&
	             wave_cur->w_type[num_waves - i - 2] == SHOCK))
	        {
	            /* Single shock */
	            sample_shock_state(state,ans,state_type,Ts,um,pm,m,
	                               spdnew,sample_speed,sizest,l_or_r);
	            return;
	        }
	        else
	        {
	                /* Split shock */
	            us = vel(0,wave_cur->st[num_waves-1]);
	            ps = Press(wave_cur->st[num_waves - i]);
	            rhos = Dens(wave_cur->st[num_waves - 1]);
	            ms = sign * (ps - pm) / (us - um);
	            shockspd = us + (sign) * ms / rhos;
	            ms = sign * (ps - p) / (us - u);
	            shockspd1 = u + (sign) * ms / rho;
	            sample_split_shock(sample_speed,Ts,state,ans,
	                               wave_cur->st[num_waves - 1],
	                               spdnew,um,m,shockspd,shockspd1,
	                               sizest,l_or_r);
	            return;
	        }
	    }
	}
	else
	{
	    /* In rarefaction branch */
	    for (i = 0; i < num_waves; ++i)
	    {
	        if (wave_cur->w_type[i] == SHOCK)
	            continue;

	        press = Press(wave_cur->st[i]);
	        if (press < pm)
	            continue;
	        if (wave_cur->w_type[i] == COMPOSITE)
	        {
	            if (wave_cur->w_type[i + 1] == COMPOSITE)
	            {
	                /* Composite ends in rarefaction */
	                p1 = Press(wave_cur->st[i + 2]);
	                uend = wave_cur->uend[NO_PTS_ON_COMP - 1];
	                pend = wave_cur->pend[NO_PTS_ON_COMP - 1];
	                rend = wave_cur->rend[NO_PTS_ON_COMP - 1];
	                ustart = wave_cur->ustart[NO_PTS_ON_COMP - 1];
	                pstart = wave_cur->pstart[NO_PTS_ON_COMP - 1];
	                rstart = wave_cur->rstart[NO_PTS_ON_COMP - 1];
	                rhoc = wave_cur->rhoc[NO_PTS_ON_COMP - 1];
	                Dens(st0) = rend;
	                Press(st0) = pend;
	                Vel(st0)[0] = uend;
	                Set_params(st0,Ts);
	                set_type_of_state(st0,TGAS_STATE);
			reset_gamma(st0);
	                state_on_adiabat_with_pr(st0, pm, mid, state_type);
	                shockspd = ustart + (sign) * rhoc / rstart;
	                sample_comp_ending_in_raref(sample_speed,Ts,state,ans,
	                                            mid,st0,wave_cur->st[i+1],
	                                            spdnew,state_type,p1,
	                                            pstart,shockspd,
	                                            sizest,l_or_r);
	                return;
	            }
	            else
	            {
	                /* Composite ends in shock */
	                state_on_comp(pm,&rstart,&pstart,&ustart,&rhoc,
				      wave_cur);
	                p1 = Press(wave_cur->st[i + 1]);
	                shockspd = ustart + (sign) * rhoc / rstart;
	                sample_comp_ending_in_shock(sample_speed,Ts,state,ans,
	                                            mid,wave_cur->st[i+1],m,
	                                            spdnew,state_type,p1,
	                                            pstart,rstart,shockspd,
	                                            sizest,l_or_r);
	                return;
	            }
	        }
	        if (wave_cur->w_type[i] == RAREFACTION)
	        {
	            if (wave_cur->w_type[i + 1] != RAREFACTION)
	            {
	                /* Single rarefaction */
	                sample_raref_state(sample_speed,Ts,state,ans,sizest,
	                                   state_type,l_or_r,um,pm,spdnew);
	                return;
	            }
	            else
	            {
	                /* Split rarefaction */
	                sample_split_raref_state(sample_speed,Ts,state,
	                                         wave_cur->st[i],ans,sizest,
	                                         state_type,l_or_r,um,
	                                         pm,spdnew);
	                return;
	            }
	        }
	    }
	}
}		/*end get_sample_phase_state*/

LOCAL	void sample_shock_state(
	Locstate	state,
	Locstate	ans,
	int		state_type,
	Locstate	Ts,
	double		um,
	double		pm,
	double		m,
	double		*spdnew,
	double		sample_speed,
	size_t		sizest,
	int		l_or_r)
{
	double		shockspd, u, rho, sign;

	rho = Dens(Ts);
	u = vel(0,Ts);

	if (l_or_r == LEFT_FAMILY)
	    sign = -1.0;
	else
	    sign = 1.0;

	shockspd = u + (sign) * m / rho;

	switch (l_or_r)
	{
	case LEFT_FAMILY:
	    if (sample_speed <= shockspd)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + sound_speed(Ts);
	    }
	    else	/* Right of left shock */
	    {
	        Dens(ans) = m / (um - shockspd);
	        Set_params(ans,Ts);
	        set_vel_and_pr_across_shock(ans,state_type,Ts,rho,pm,um);
		reset_gamma(ans);
	        *spdnew = fabs(um) + sound_speed(ans);
	    }
	    break;
	case RIGHT_FAMILY:
	    if (sample_speed >= shockspd)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + sound_speed(Ts);
	    }
            else    /* Left of right shock */
	    {
	        Dens(ans) = -m / (um - shockspd);
	        Set_params(ans,Ts);
	        set_vel_and_pr_across_shock(ans,state_type,Ts,rho,pm,um);
		reset_gamma(ans);
	        *spdnew = fabs(um) + sound_speed(ans);
	    }
	}
}		/*end sample_shock_state*/

LOCAL	void sample_split_shock(
	double		sample_speed,
	Locstate	Ts,
	Locstate	state,
	Locstate	ans,
	Locstate	st0,
	double		*spdnew,
	double		um,
	double		m,
	double		shockspd,
	double		shockspd1,
	size_t		sizest,
	int		l_or_r)
{
	double		u, us;

	u = vel(0,Ts);
	us = vel(0,st0);

	switch (l_or_r)
	{
	case LEFT_FAMILY:
	    if (sample_speed >= shockspd)
	    {
	        Dens(ans) = m / (um - shockspd);
	        Set_params(ans,Ts);
	        Vel(ans)[0] = um;
	        set_type_of_state(ans,TGAS_STATE);
		reset_gamma(ans);
	        *spdnew = fabs(um) + sound_speed(ans);
	    }
	    else if (sample_speed <= shockspd1)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + sound_speed(Ts);
	    }
	    else
	    {
	         /* In split shock region */
	        ft_assign(ans,st0,sizest);
	        *spdnew = fabs(us) + sound_speed(st0);
	    }
	    break;
	case RIGHT_FAMILY:
	    if (sample_speed <= shockspd)
	    {
	        Dens(ans) = -m / (um - shockspd);
	        Set_params(ans,Ts);
	        Vel(ans)[0] = um;
	        set_type_of_state(ans,TGAS_STATE);
		reset_gamma(ans);
	        *spdnew = fabs(um) + sound_speed(ans);
	    }
	    else if (sample_speed >= shockspd1)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + sound_speed(Ts);
	    }
	    else
	    {
	         /* In split shock region */
	        
	        ft_assign(ans,st0,sizest);
	        *spdnew = fabs(us) + sound_speed(st0);
	    }
	}
}		/*end sample_split_shock*/


LOCAL	void sample_comp_ending_in_raref(
	double		sample_speed,
	Locstate	Ts,
	Locstate	state,
	Locstate	ans,
	Locstate	mid,
	Locstate	st0,
	Locstate	stph,
	double		*spdnew,
	int		state_type,
	double		p1,
	double		pstart,
	double		shockspd,
	size_t		sizest,
	int		l_or_r)
{
	double		u, us, uend, um, p, cm, c1, c0, cl, cr;

	uend = vel(0,st0);
	um = vel(0,mid);
	cm = sound_speed(mid);
	c1 = sound_speed(st0);
	c0 = sound_speed(Ts);
	p = pressure(Ts);
	u = vel(0,Ts);
        get_ph_sound_spd(&cl,&cr,stph);
	us = vel(0,stph);

	if (pstart <= p)
	{
	    /* Solution is single shock followed by rarefaction */
	    switch (l_or_r)
	    {
	    case LEFT_FAMILY:
	        if (sample_speed <= shockspd)
	        {
	            ft_assign(ans, state, sizest);
	            *spdnew = fabs(u) + sound_speed(Ts);
	        }
	        if (Between(sample_speed,shockspd,uend - c1))
	        {
	            ft_assign(ans, st0, sizest);
	            *spdnew = fabs(uend) + sound_speed(ans);
	        }
	        if (Between(sample_speed, uend - c1, um - cm))
	        {
	                            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,uend,st0,
	                                                 mid,ans,state_type,
	                                                 spdnew,LEFT_FAMILY);
	        }
	        else	/* Right of rarefaction */ 
	        {
                    ft_assign(ans, mid, sizest);
	            set_vel_and_pr_across_raref(ans,state_type,um);
	            *spdnew = fabs(um) + cm;
	        }
	        break;
	    case RIGHT_FAMILY:
	        if (sample_speed >= shockspd)
	        {
	            ft_assign(ans,state,sizest);
	            *spdnew = fabs(u) + sound_speed(Ts);
	        }
	        if (Between(sample_speed,shockspd,uend + c1))
	        {
	            ft_assign(ans,st0,sizest);
	            *spdnew = fabs(uend) + sound_speed(ans);
	        }
	        if (Between(sample_speed, uend + c1, um + cm))
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,uend,st0,
	                                                 mid,ans,state_type,
	                                                 spdnew,LEFT_FAMILY);
	        }
	        else
	        {
	            ft_assign(ans, mid, sizest);
	            set_vel_and_pr_across_raref(ans, state_type, um);
	            *spdnew = fabs(um) + cm;
	        }
	    }
	}
	else
	{
	    /*Solution is rarefaction followed by shock and then a rarefaction*/
	    switch (l_or_r)
	    {
	    case LEFT_FAMILY:
	        if (sample_speed <= u - c0)
	        {
	            ft_assign(ans,state,sizest);
	            *spdnew = fabs(u) + c0;
	        }
	        if (sample_speed >= um - cm)
	        {
	            ft_assign(ans,mid,sizest);
	            set_vel_and_pr_across_raref(ans,state_type,um);
	            *spdnew = fabs(um) + cm;
	        }
	        if (Between(pstart,p,p1))
	        {
	            if (Between(sample_speed,u - c0,shockspd))
	            {
	                /* In rarefaction fan */
	                (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                     stph,ans,
	                                                     state_type,spdnew,
	                                                     LEFT_FAMILY);
	            }
	        }
	        else
	        {
	            if (Between(sample_speed,us - cl,us - cr))
	            {
	                /* In constant state */
	                ft_assign(ans,stph,sizest);
	                *spdnew = fabs(us) + max(cr,cl);
	            }
	            if (Between(sample_speed,u-c0,us-cl))
	            {
	                (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                     mid,ans,
	                                                     state_type,spdnew,
	                                                     LEFT_FAMILY);
	            }
	            if (Between(sample_speed,us-cr,shockspd))
	            {
	                /* In rarefaction fan */
	                (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                     mid,ans,
	                                                     state_type,spdnew,
	                                                     LEFT_FAMILY);
	            }
	        }
	        if (Between(sample_speed,shockspd,uend - c1))
	        {
	            /* In constant state */
	            ft_assign(ans,stph,sizest);
	            *spdnew = fabs(us) + max(cr,cl);
	        }
	        if (Between(sample_speed, uend - c1, um - cm))
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,uend,stph,
	                                                 mid,ans,state_type,
	                                                 spdnew,LEFT_FAMILY);
	        }
	        else	/* Right of rarefaction */ 
	        {
                    ft_assign(ans, mid, sizest);
	            set_vel_and_pr_across_raref(ans,state_type,um);
	            *spdnew = fabs(um) + cm;
	        }
	        break;
	    case RIGHT_FAMILY:
	        if (sample_speed >= u + c0)
	        {
	            ft_assign(ans, state, sizest);
	            *spdnew = fabs(u) + c0;
	        }
	        if (sample_speed <= um + cm)
	        {
	            ft_assign(ans, mid, sizest);
	            set_vel_and_pr_across_raref(ans, state_type, um);
	            *spdnew = fabs(um) + cm;
	        }
	        if (Between(pstart,p,p1))
	        {
	            if (Between(sample_speed, u + c0,shockspd))
	            {
	                (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                     mid,ans,
	                                                     state_type,spdnew,
	                                                     RIGHT_FAMILY);
	            }
	        }
	        else
	        {
	            if (Between(sample_speed, us + cl,us + cr))
	            {
	                ft_assign(ans,stph,sizest);
	                *spdnew = fabs(us) + max(cl, cr);
	            }
	            if (Between(sample_speed,shockspd, us + cr))
	            {
	                 /* In rarefaction fan */
	                 (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                      mid,ans,
	                                                      state_type,spdnew,
	                                                      RIGHT_FAMILY);
	            }
	            if (Between(sample_speed,u + c0,us +cr))
	            {
                        (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                     mid,ans,
	                                                     state_type,spdnew,
	                                                     RIGHT_FAMILY);
	            }
	        }
	        if (Between(sample_speed, uend + c1, um + cm))
	        {
	            (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                 mid,ans,state_type,
	                                                 spdnew,RIGHT_FAMILY);
	        }
	        else
	        {
                    ft_assign(ans, mid, sizest);
	            set_vel_and_pr_across_raref(ans,state_type,um);
	            *spdnew = fabs(um) + cm;
	        }
	    }
	}
}		/*end sample_comp_ending_in_raref*/

LOCAL	void sample_comp_ending_in_shock(
	double		sample_speed,
	Locstate	Ts,
	Locstate	state,
	Locstate	ans,
	Locstate	mid,
	Locstate	stph,
	double		m,
	double		*spdnew,
	int		state_type,
	double		p1,
	double		pstart,
	double		rstart,
	double		shockspd,
	size_t		sizest,
	int		l_or_r)
{
	double		u, us, um, p, pm, c0, cl, cr;

	um = vel(0,mid);
	c0 = sound_speed(Ts);
	u = vel(0,Ts);
	pm = pressure(mid);
	p = pressure(Ts);
        get_ph_sound_spd(&cl,&cr,stph);
	us = vel(0,stph);

	if (pstart <= p)
	{
	    /* Solution is rarefying shock */
	    switch (l_or_r)
	    {
	    case LEFT_FAMILY:
	        if (sample_speed >= shockspd)
	        {
	            ft_assign(ans, state, sizest);
	            *spdnew = fabs(u) + c0;
	        }
	        else
	        {
	            Dens(ans) = m / (um - shockspd);
	            Set_params(ans,Ts);
	            set_vel_and_pr_across_shock(ans,state_type,Ts,rstart,pm,um);
		    reset_gamma(ans);
	            *spdnew = fabs(um) + sound_speed(ans);
	        }
	        break;
	    case RIGHT_FAMILY:
	        if (sample_speed <= shockspd)
	        {
	            ft_assign(ans,state,sizest);
	            *spdnew = fabs(u) + sound_speed(Ts);
	        }
	        else
	        {
	            Dens(ans) = -m / (um - shockspd);
	            Set_params(ans,Ts);
	            set_vel_and_pr_across_shock(ans,state_type,Ts,rstart,pm,um);
		    reset_gamma(ans);
	            *spdnew = fabs(um)+sound_speed(ans);
	        }
	    }
	}
	if (Between(p, pstart, p1))
	{
	    /* Solution is rarefaction followed by rarefying shock */
	    switch (l_or_r)
	    {
	    case LEFT_FAMILY:
	        if (sample_speed <= u - c0)
	        {
	            ft_assign(ans,state,sizest);
	            *spdnew = fabs(u) + c0;
	        }
	        if (sample_speed > shockspd)
	        {
	            Dens(ans) = m / (um - shockspd);
	            Set_params(ans,Ts);
	            set_vel_and_pr_across_shock(ans,state_type,Ts,rstart,pm,um);
		    reset_gamma(ans);
	            *spdnew = fabs(um) + sound_speed(ans);
	        }
	        else
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,u,
	                                                 Ts,mid,ans,state_type,
	                                                 spdnew,LEFT_FAMILY);
	        }
	        break;
	    case RIGHT_FAMILY:
	        if (sample_speed >= u + c0)
	        {
	            ft_assign(ans,state,sizest);
	            *spdnew = fabs(u) + c0;
	        }
	        if (sample_speed < shockspd)
	        {
	            Dens(ans) = - m / (um - shockspd);
	            Set_params(ans,Ts);
	            set_vel_and_pr_across_shock(ans,state_type,Ts,rstart,pm,um);
		    reset_gamma(ans);
	            *spdnew = fabs(um) + sound_speed(ans);
	        }
	        else
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,mid,
	                                                 ans,state_type,spdnew,
	                                                 RIGHT_FAMILY);
	        }
	    }
	}
	else
	{
	    /* Solution is split rarefaction followed by shock */
	    switch (l_or_r)
	    {
	    case LEFT_FAMILY:
	        if (sample_speed <= u - c0)
	        {
	            ft_assign(ans, state, sizest);
	            *spdnew = fabs(u) + c0;
	        }
	        if (Between(sample_speed, u - c0, us - cr))
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,stph,
	                                                 ans,state_type,
	                                                 spdnew,LEFT_FAMILY);
	        }
	        if (Between(sample_speed, us - cl,us - cr))
	        {
	            ft_assign(ans,stph,sizest);
	            *spdnew = fabs(us) + max(cl, cr);
	        }
	        if (Between(sample_speed, us - cl,shockspd))
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,us,stph,
	                                                 mid,ans,state_type,
	                                                 spdnew,LEFT_FAMILY);
	        }
	        else
	        {
	            /* In midstate */
	            Dens(ans) = m / (um - shockspd);
	            Set_params(ans,Ts);
	            set_vel_and_pr_across_shock(ans,state_type,Ts,rstart,pm,um);
		    reset_gamma(ans);
	            *spdnew = fabs(um) + sound_speed(ans);
	        }
	        break;
	    case RIGHT_FAMILY:
	        if (sample_speed >= u + c0)
	        {
	            ft_assign(ans,state,sizest);
	            *spdnew = fabs(u) + c0;
	        }
	        if (Between(sample_speed, u + c0, us + cr))
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,stph,
	                                                 ans,state_type,spdnew,
	                                                 RIGHT_FAMILY);
	        }
	        if (Between(sample_speed, us + cl,us + cr))
	        {
	            ft_assign(ans,stph,sizest);
	            *spdnew = fabs(us) + max(cl, cr);
	        }
	        if (Between(sample_speed, us + cl,shockspd))
	        {
	            /* In rarefaction fan */
	            (void) oned_state_in_rarefaction_fan(sample_speed,us,stph,
	                                                 mid,ans,state_type,
	                                                 spdnew,RIGHT_FAMILY);
	        }
	        else
	        {
	            /* In  midstate */
	            Dens(ans) = - m / (um - shockspd);
	            Set_params(ans,Ts);
	            set_vel_and_pr_across_shock(ans,state_type,Ts,rstart,pm,um);
		    reset_gamma(ans);
	            *spdnew = fabs(um) + sound_speed(ans);
	        }
	    }
	}
}		/*end sample_comp_ending_in_shock*/

LOCAL void sample_raref_state(
	double		sample_speed,
	Locstate	Ts,
	Locstate	state,
	Locstate	ans,
	size_t		sizest,
	int		state_type,
	int		l_or_r,
	double		um,
	double		pm,
	double		*spdnew)
{
	double		u, cm, c0; 
	static Locstate mid = NULL;

	c0 = sound_speed(Ts);
	u = vel(0,Ts);

	switch (l_or_r)
	{
	case LEFT_FAMILY:
	    if (sample_speed <= u - c0)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + c0;
	    }
	    else
	    {
	        state_on_adiabat_with_pr(Ts,pm,mid,state_type);
	        cm = sound_speed(mid);
	        if (sample_speed >= um - cm)
	        {
	            set_vel_and_pr_across_raref(ans,state_type,um);
	            ft_assign(ans,mid,sizest);
	            *spdnew = fabs(um) + cm;
	        }
	        else
	            (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,mid,
	                                                 ans,state_type,spdnew,
	                                                 LEFT_FAMILY);
	    }
	    break;
	case RIGHT_FAMILY:
	    if (sample_speed >= u + c0)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + c0;
	    }
	    else
	    {
	        state_on_adiabat_with_pr(Ts,pm,mid,state_type);
	        cm = sound_speed(mid);
	        if (sample_speed <= um + cm)
	        {
	            set_vel_and_pr_across_raref(ans,state_type,um);
	            ft_assign(ans,mid,sizest);
	            *spdnew = fabs(um) + cm;
	        }
	        else
	            (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,mid,
	                                                 ans,state_type,spdnew,
	                                                 LEFT_FAMILY);
	    }
	}
}		/*end sample_raref_state*/

LOCAL void sample_split_raref_state(
	double		sample_speed,
	Locstate	Ts,
	Locstate	state,
	Locstate	stph,
	Locstate	ans,
	size_t		sizest,
	int		state_type,
	int		l_or_r,
	double		um,
	double		pm,
	double		*spdnew)
{
	double		u, c0, cl, cr, cm, us; 
	static Locstate mid = NULL;

	if (mid == NULL)
	{
	    (*Params(Ts)->_alloc_state)(&mid, Params(Ts)->sizest);
	}

	c0 = sound_speed(Ts);
	u = vel(0,Ts);

	switch (l_or_r)
	{
	case LEFT_FAMILY:
	    if (sample_speed <= u - c0)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + c0;
	    }
	    else
	    {
	        state_on_adiabat_with_pr(Ts, pm, mid, state_type);
	        cm = sound_speed(mid);
	        us = vel(0,stph);
	        get_ph_sound_spd(&cl, &cr,stph);
	        if (sample_speed >= um - cm)
	        {
	            ft_assign(ans, mid, sizest);
	            set_vel_and_pr_across_raref(ans, state_type, um);
	            *spdnew = fabs(um) + cm;
	        }
	        else if (Between(sample_speed,us - cl,us - cr))
	        {
	                        /* In constant state */
	            debug_print("sample_phase","In split state\n");
	            set_state(ans,state_type,stph);
	            set_vel_and_pr_across_raref(ans,state_type,us);
	            *spdnew = fabs(us) + max(cr,cl);
	        }
	        else
	        {
	            debug_print("sample_phase","In rarefaction fan\n");
	            if (Between(sample_speed,us-cl,u - c0))
	                (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                     stph,ans,
	                                                     state_type,spdnew,
	                                                     LEFT_FAMILY);
	            else if (Between(sample_speed,us-cr,um - cm))
	                (void) oned_state_in_rarefaction_fan(sample_speed,us,
	                                                     stph,mid,ans,
	                                                     state_type,spdnew,
	                                                     LEFT_FAMILY);
	        }
	    }
	    break;
	case RIGHT_FAMILY:
	    if (sample_speed >= u + c0)
	    {
	        ft_assign(ans, state, sizest);
	        *spdnew = fabs(u) + c0;
	    }
	    else
	    {
	        state_on_adiabat_with_pr(Ts, pm, mid, state_type);
	        cm = sound_speed(mid);
	        us = vel(0,stph);
	        get_ph_sound_spd(&cl, &cr,stph);
	        if (sample_speed <= um + cm)
	        {
	            ft_assign(ans, mid, sizest);
	            set_vel_and_pr_across_raref(ans, state_type, um);
	            *spdnew = fabs(um) + cm;
	        }
	        else if (Between(sample_speed,us + cl,us + cr))
	        {
	                    /* In constant state */
	            debug_print("sample_phase","In split state\n");
	            set_state(ans,state_type,stph);
	            set_vel_and_pr_across_raref(ans, state_type, us);
	                    *spdnew = fabs(us) + max(cr, cl);
	        }
	        else
	        {
	            if (Between(sample_speed,us+cl,u + c0))
	                (void) oned_state_in_rarefaction_fan(sample_speed,u,Ts,
	                                                     stph,ans,
	                                                     state_type,spdnew,
	                                                     RIGHT_FAMILY);
	            if (Between(sample_speed,us+cr,um + cm))
	                (void) oned_state_in_rarefaction_fan(sample_speed,us,
	                                                     stph,mid,ans,
	                                                     state_type,spdnew,
	                                                     RIGHT_FAMILY);
	        }
	    }
	}
}		/*end sample_split_raref_state*/

LOCAL void set_vel_and_pr_across_shock(
	Locstate	ans,
	int		state_type,
	Locstate	Ts,
	double		rho,
	double		pm,
	double		um)
{
	set_type_of_state(ans,state_type);
	switch (state_type)
	{
	case TGAS_STATE:
	    Vel(ans)[0] = um;
	    Press(ans) = pm;
	    return;
	case GAS_STATE:
	    Mom(ans)[0] = Dens(ans)*um;
	    Energy(ans) = Dens(ans)*(0.5 * sqr(um) +
	                             specific_internal_energy(Ts) +
	                             0.5 * (pressure(Ts) + pm) *
	                             (1.0/rho - 1.0/Dens(ans)));
	    reset_gamma(ans);
	    return;
	case EGAS_STATE:
	    Vel(ans)[0] = um;
	    Energy(ans) = specific_internal_energy(Ts) +
	                  0.5*(pressure(Ts) + pm)*(1.0/rho - 1.0/Dens(ans));
	    reset_gamma(ans);
	    return;
	}
}		/*end set_vel_and_pr_across_shock*/

LOCAL void set_vel_and_pr_across_raref(
	Locstate	ans,
	int		state_type,
	double		um)
{
	set_type_of_state(ans,state_type);
	switch (state_type)
	{
	case TGAS_STATE:
	case EGAS_STATE:
	    Vel(ans)[0] = um;
	    break;
	case GAS_STATE:
	    Mom(ans)[0] = Dens(ans) * um;
	    Energy(ans) += 0.5 * um * Mom(ans)[0];
	    reset_gamma(ans);
	    break;
	}
}		/*end set_vel_and_pr_across_raref*/
#endif /* defined(PHASE_CODE) && defined(TWOD) */

LOCAL	void	print_Riemann_wave_curves(
	Locstate Tsl,
	Locstate Tsr,
	double	pjump,
	double	p_min,
	double	pstar)
{
	double	vxr = vel(0,Tsr);
	double	vxl = vel(0,Tsl);		
	double	poj, ur_mid, ul_mid;
	double	p_max;
	double	dp, dp1, p, pl, pr;
	int	i;

	verbose_print_state("Tsl",Tsl);
	verbose_print_state("Tsr",Tsr);
	(void) printf("pjump = %"FFMT", p_min = %"FFMT"\n",pjump,p_min);

	(void) printf("\n%-14s %-14s %-14s %-14s %-14s\n","pressure","poj",
	              "ur_mid","ul_mid","ur_mid-ul_mid");
	pl = pressure(Tsl);
	pr = pressure(Tsr);
	p_max = 2.0*max(pl,pr);
	dp = (p_max - p_min)/100.0;
	dp1 = (pstar-p_min)/10.0;
	dp = min(dp,dp1);
	for (i = 0; i <= 100; ++i)
	{
	    p = p_min + i*dp;
	    poj = floor_pressure(p+pjump,p_min);
	    ur_mid = vxr + riemann_wave_curve(Tsr,p);
	    ul_mid = vxl - riemann_wave_curve(Tsl,poj);
	    (void) printf("%-14g %-14g %-14g %-14g %-14g\n",p,poj,
	                  ur_mid,ul_mid,ur_mid-ul_mid);
	}
	(void) printf("\n");
}		/*end print_Riemann_wave_curve*/
