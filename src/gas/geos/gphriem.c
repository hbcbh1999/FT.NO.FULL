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
*				gphriem.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains phase dependent functions for the calculation of
*	Hugoniot curves, adiabatic curves and other support functions
*	needed for the solution of Riemann problems.
*
*/


#if defined(TWOD) && defined(PHASE_CODE)
#define	DEBUG_STRING 	"ph_riem"
#include <geos/sesame.h>

enum _LIQUID_BNDRY_TYPE {
	DENSITY_GIVEN = 1,
	ENERGY_GIVEN
};
typedef enum _LIQUID_BNDRY_TYPE LIQUID_BNDRY_TYPE;

	/* LOCAL Function Declarations */
LOCAL	boolean	dble_sonic_bisect(double,double*,POINTER);
LOCAL	boolean	is_liquid_bndry(double,LIQUID_BNDRY_TYPE,double*);
LOCAL	boolean	num_Hug_by_mass_flux_squared(double,Locstate,Locstate,
					     double,double,double);
LOCAL	boolean	ses_hug_bisect(double,double*,POINTER);
LOCAL	boolean	ses_ph_bisect(double,double*,POINTER);
LOCAL	boolean	sonic_shock(Locstate,Locstate,boolean,double,int);
LOCAL	void	get_Hug_intersection(Locstate,Locstate,int,WAVE_CURVE*);
LOCAL	void	get_intersection_state(Locstate,POINT*,double,double*,
				       double*,double*,int*);
LOCAL	void	get_phase_Smax_state(SESAME_EOS*,POINT**);
LOCAL	void	get_rare_intrst_state(Locstate,double,double*,
				      double*,double*,int*);
LOCAL	void	make_composite(Locstate,double,double,Locstate,Locstate,
			       WAVE_CURVE*);

/*
*			intrsct_wv_crv_wth_phs_bdry():
*
*	Do preliminary work to determine if special wave structure may occur
*	due to the inclusion of a phase boundary.
*
*	Returns NO if no intersections are found,  YES if at least one
*	wave curve crosses the phase boundary;
*/

EXPORT int intrsct_wv_crv_wth_phs_bdry(
	Locstate 	         Tsl,
	Locstate 	         Tsr,
	double		         *rpi,
	double		         *rui,
	double		         *rri,
	int		         *nri,
	double		         *lpi,
	double		         *lui,
	double		         *lri,
	int		         *nli,
	RIEMANN_SOLVER_WAVE_TYPE l_wave,
	RIEMANN_SOLVER_WAVE_TYPE r_wave)
{
	double		Sl, Sr;
	POINT		*pmax;

	DEBUG_ENTER(intrsct_wv_crv_wth_phs_bdry)
	*nri = *nli = 0;

	if ((multiphase_eos(Eos(Tsl)) == NO) || 
	    (multiphase_eos(Eos(Tsr)) == NO))
	{
		DEBUG_LEAVE(intrsct_wv_crv_wth_phs_bdry)
		return NO;
	}

	Sl = entropy(Tsl);	Sr = entropy(Tsr);

	if(Ses_params(Tsr)->Smax < Sl && Ses_params(Tsl)->Smax < Sr)
	{
		DEBUG_LEAVE(intrsct_wv_crv_wth_phs_bdry)
		return NO;
	}
	else
	{
		if (DEBUG)
		{
			(void) printf("Left state\n");
			(void) printf("%g %g\n",pressure(Tsl),vel(0,Tsl));
			(void) printf("Right state\n");
			(void) printf("%g %g\n",pressure(Tsr),vel(0,Tsr));
			(void) printf("Left wave and right wave\n");
			(void) printf("%d %d\n",l_wave,r_wave);
		}
		get_phase_Smax_state(Ses_params(Tsr),&pmax);
		if (Ses_params(Tsl)->Smax > Sl)
		{

			if (l_wave == RAREFACTION)
			{
				get_rare_intrst_state(Tsl,-1.0,
					lui,lpi,lri,nli);
			}
			else if (l_wave == SHOCK)
			{
				get_intersection_state(Tsl,pmax,-1.0,
					lui,lpi,lri,nli);
			}
		}
		if (Ses_params(Tsr)->Smax > Sr)
		{
			if (r_wave == RAREFACTION)
			{
				get_rare_intrst_state(Tsr,1.0,
					rui,rpi,rri,nri);
			}
			else if (r_wave == SHOCK)
			{
				get_intersection_state(Tsr,pmax,1.0,
					rui,rpi,rri,nri);
			}
		}
	}		
	if (DEBUG)
		(void) printf("%d %d\n",*nri,*nli);
        DEBUG_LEAVE(intrsct_wv_crv_wth_phs_bdry)
	return (*nri > 0 || *nli > 0) ? YES : NO;
}		/*end intrsct_wv_crv_wth_phs_bdry*/

LOCAL	void get_phase_Smax_state(
	SESAME_EOS 	*seos,
	POINT		**pmax)
{
	double		Smax = seos->Smax;
	INTERFACE 	*intfc = seos->fr[SESAME_RHO_ENTROPY]->interf;
	CURVE		**cur;
	BOND		*bond;

        DEBUG_ENTER(get_phase_Smax_state)
	for(cur = intfc->curves; cur && *cur;cur++)
	{
		if (wave_type(*cur) != PHASE_BOUNDARY) continue;

		for (bond = (*cur)->first; bond != NULL; bond = bond->next)
		{
			if(Coords(bond->start)[1] >= Smax)
			{
				*pmax = bond->start;
				DEBUG_LEAVE(get_phase_Smax_state)
				return;
			}
		}
	}		
        DEBUG_LEAVE(get_phase_Smax_state)
}		/*end get_phase_Smax_state*/

typedef struct {
	Front		*fr;
	Wave		*wv;
	double		slope;
	double		e0, p0, rho0;
	double		p, rho;
	double		x0, y0;
	CURVE		*phsbdry;
	SESAME_EOS	*seos;
} PH_HUG_PARAMS;

LOCAL	void get_intersection_state(
	Locstate 	Ts,
	POINT		*pmax,
	double		sign,
	double		*ui,
	double		*pi,
	double		*ri,
	int		*ni)
{
	PH_HUG_PARAMS 	ph_inv_params;
	double		rho0, p0, e0, u0, ymin;
	double		hug_s, hug_e;
	double		rho_s, rho_e, p_s, p_e, e_s, e_e;
	double		p, rho;
	double		epsilon, delta;
	SESAME_EOS	*seos = Ses_params(Ts);
	Front		*fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave		*wv = seos->wave[SESAME_RHO_ENTROPY];
	double		coords[2];
	Locstate 	state, state1;
	INTERFACE 	*intfc = fr->interf;
	CURVE		**cur, *phsbdry;
	BOND		*bond;
	static const double	HUGEPS = 1.0e-05; /*TOLERANCE*/

        DEBUG_ENTER(get_intersection_state)

	ph_inv_params.seos = seos;
	for(cur = intfc->curves;cur && *cur;cur++)
	{
		if(wave_type(*cur) == PHASE_BOUNDARY)
		{
			phsbdry = *cur;
			break;
		}
	}

	if (phsbdry == NULL)
	{
		screen("ERROR in get_intersection_state(), ");
		screen("can't find phase boundary\n");
		clean_up(ERROR);
	}

	rho0 = Dens(Ts);
	p0 = pressure(Ts);
	e0 = specific_internal_energy(Ts);
	u0 = vel(0,Ts);
	ymin = ses_rs_grid_from_entpy(entropy(Ts),seos);
	rho = ses_rs_rho_from_grid(Coords(pmax)[0],seos);
	for (bond = phsbdry->first; bond != NULL; bond = bond->next)
	{
		if ((Coords(bond->start)[1] < ymin) &&
		    (Coords(bond->end)[1] < ymin)) continue;

		state = left_state_at_point_on_curve(bond->start,bond,phsbdry);
		rho_s = ses_rs_rho_from_grid(Coords(bond->start)[0],seos);

		if (!is_retrograde_bndry(rho_s,Ts)) continue;

		e_s = ses_rs_colde(state) + ses_rs_rede(state)*
			(ses_rs_temp_from_var(ses_rs_Tvar(state),seos));
		p_s = ses_rs_coldp(state) + ses_rs_redp(state)*
			(ses_rs_temp_from_var(ses_rs_Tvar(state),seos))*rho_s;
		hug_s = (e_s - e0) + 0.5*(p_s + p0)*(1/rho_s - 1/rho0);
		state1 = left_state_at_point_on_curve(bond->end,bond,phsbdry);
		rho_e = ses_rs_rho_from_grid(Coords(bond->end)[0],seos);
		e_e = ses_rs_colde(state1) + ses_rs_rede(state1)*
			(ses_rs_temp_from_var(ses_rs_Tvar(state1),seos));
		p_e = ses_rs_coldp(state1) + ses_rs_redp(state1)*
			(ses_rs_temp_from_var(ses_rs_Tvar(state1),seos))*rho_e;
		hug_e = (e_e - e0) + 0.5*(p_e + p0)*(1/rho_e - 1/rho0);

		if (hug_e*hug_s > 0.0) continue;

		ph_inv_params.fr = fr;
		ph_inv_params.wv = wv;
		ph_inv_params.phsbdry = phsbdry;
		ph_inv_params.x0 = Coords(bond->start)[0];
		ph_inv_params.y0 = Coords(bond->start)[1];
		ph_inv_params.e0 = e0;
		ph_inv_params.p0 = p0;
		ph_inv_params.rho0 = rho0;
		ph_inv_params.slope =
			(Coords(bond->start)[0] - Coords(bond->end)[0]) /
			(Coords(bond->start)[1] - Coords(bond->end)[1]);
		epsilon = HUGEPS*e0;
		delta = HUGEPS*bond_length(bond);
		if (bisection_find_root(ses_ph_bisect,
				(POINTER)&ph_inv_params,0.0,&coords[1],
				Coords(bond->start)[1],Coords(bond->end)[1],
				epsilon,delta) == FUNCTION_FAILED)
		{
			screen("ERROR in get_intersection_state(), ");
			screen("unable to find root on left branch\n");
			clean_up(ERROR);
		}
		ri[*ni] = rho = ph_inv_params.rho;
        	pi[*ni] = p = ph_inv_params.p;
		ui[*ni] = u0 + sign*sqrt((1./rho0 - 1./rho)*(p - p0));
		(*ni)++;
	}
	if (DEBUG)
		(void) printf("%d\n",*ni);
        DEBUG_LEAVE(get_intersection_state)
}		/*end get_intersection_state*/

LOCAL boolean ses_ph_bisect(
	double		y,
	double		*hug,
	POINTER 	params)
{
	PH_HUG_PARAMS 	*prms = (PH_HUG_PARAMS *)params;
	CURVE		*phsbdry = prms->phsbdry;
	Front		*fr = prms->fr;
	Wave		*wv = prms->wv;
	double		rho0 = prms->rho0;
	double		e0 = prms->e0;
	double		p0 = prms->p0;
	double		slope = ((PH_HUG_PARAMS *) params)->slope;
	double		x0 = ((PH_HUG_PARAMS *) params)->x0;
	double		y0 = ((PH_HUG_PARAMS *) params)->y0;
	double		e, rho, p, T;
	double		coords[MAXD];
	double		var[NUM_SES_VAR];
	COMPONENT 	comp = COMP_PURE_PHASE;

        DEBUG_ENTER(ses_ph_bisect)
	coords[0] = x0 + slope*(y - y0);
	coords[1] = y;
	set_ses_intrp_flag(EVALUATE_PRESSURE|EVALUATE_ENERGY,
		           SESAME_RHO_ENTROPY);
	ses_solution(coords,comp,Hyper_surf(phsbdry),POSITIVE_SIDE,fr,
					(POINTER)wv,var);

	T = ses_rs_temp_from_var(var[RS_T],prms->seos);
	rho = prms->rho = ses_rs_rho_from_grid(coords[0],prms->seos);
        p = prms->p = rho*T*var[RS_RP] + var[RS_CP];
        e = T*var[RS_RE] + var[RS_CE];
	*hug = (e - e0) + 0.5*(p+p0)*(1.0/rho - 1.0/rho0);
	DEBUG_LEAVE(ses_ph_bisect)
	return FUNCTION_SUCCEEDED;
}		/*end ses_ph_bisect*/



LOCAL	void get_rare_intrst_state(
	Locstate 	Ts,
	double		sign,
	double		*ui,
	double		*pi,
	double		*ri,
	int		*ni)
{
	double		rho0_grid, u0, S0;
	double		tmp;
	SESAME_EOS	*seos = Ses_params(Ts);
	Front		*fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave		*wv = seos->wave[SESAME_RHO_ENTROPY];
	INTERFACE 	*intfc = fr->interf;
	CURVE		**cur, *phsbdry;
	BOND		*bond;
	COMPONENT 	comp;
	double		T, var[NUM_SES_VAR];
	double		slope;
	double		coords[2];
	static const double	TOLEPS = 1.0e-04; /*TOLERANCE*/

        DEBUG_ENTER(get_rare_intrst_state)

	*ni = 0;
	phsbdry = NULL;
	for(cur = intfc->curves; cur && *cur; cur++)
	{
		if (wave_type(*cur) == PHASE_BOUNDARY)
		{
			phsbdry = *cur;
			break;
		}
	}

	if (phsbdry == NULL)
	{
		screen("ERROR in get_rare_intrst_state(), ");
		screen("can't find phase boundary\n");
		clean_up(ERROR);
	}

	rho0_grid = ses_rs_grid_from_rho(Dens(Ts),seos);
	u0 = vel(0,Ts);
	S0 = entropy(Ts);
	coords[1] = ses_rs_grid_from_entpy(S0,seos);
	comp = COMP_PURE_PHASE;
	for (bond = phsbdry->first; bond != NULL; bond = bond->next)
	{
		if (!Between(coords[1],Coords(bond->start)[1],
				          Coords(bond->end)[1]))
			continue;

		slope = (Coords(bond->start)[0] - Coords(bond->end)[0])/
			(Coords(bond->start)[1] - Coords(bond->end)[1]);
		coords[0] = Coords(bond->start)[0] +
			(coords[1] - Coords(bond->start)[1])*slope;

		if (coords[0] > rho0_grid) continue;
		if (fabs(coords[0] - rho0_grid) < TOLEPS) continue;

		set_ses_intrp_flag(EVALUATE_PRESSURE,SESAME_RHO_ENTROPY);
		ses_solution(coords,comp,Hyper_surf(phsbdry),NEGATIVE_SIDE,
			     fr,(POINTER)wv,var);
		T = ses_rs_temp_from_var(var[RS_T],seos);
		ri[*ni] = ses_rs_rho_from_grid(coords[0],seos);
        	pi[*ni] = ri[*ni]*T*var[RS_RP] + var[RS_CP];
		ui[*ni] = u0 + sign*riemann_wave_curve(Ts,pi[*ni]);
		if (fabs(ui[*ni] - u0) < TOLEPS) continue;
		(*ni)++;
	}
	if (*ni == 2 && ri[1] < ri[0])
	{
		tmp = ri[1]; ri[1] = ri[0]; ri[0] = tmp;
		tmp = pi[1]; pi[1] = pi[0]; pi[0] = tmp;
		tmp = ui[1]; ui[1] = ui[0]; ui[0] = tmp;

	}
	if (DEBUG)
		(void) printf("%d\n",*ni);
        DEBUG_LEAVE(get_rare_intrst_state)
}		/*end get_rare_intrst_state*/


/*
*			is_retrograde_bndry()
*	Determines if the phase boundary is retrograde at the intersection
*	point. This tells if anomalous behavior will be expected for the
*	Riemann problem with phase transitions, and if so what type.
*	Retruns YES or NO.
*/


EXPORT	int is_retrograde_bndry(
	double		ri,
	Locstate 	state0)
{
	SESAME_EOS	*seos = Ses_params(state0);
	Front		*fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave		*wv = seos->wave[SESAME_RHO_ENTROPY];
	INTERFACE 	*intfc = fr->interf;
	CURVE		**cur, *phsbdry;
	BOND		*bond;
	Locstate 	state, state1;
	double		slope;
	double		press, press1, V, V1, dP, dV;
	double		S_grid, adbgam;
	double		coords[MAXD];
	double		var[NUM_SES_VAR];
	double		lri;


	DEBUG_ENTER(is_retrograde_bndry)
	if(is_liquid_bndry(ri,DENSITY_GIVEN,
			   Ses_params(state0)->sestab.tbls) == YES)
	{
		DEBUG_LEAVE(is_retrograde_bndry)
		return NO;
	}

	phsbdry = NULL;
	for(cur = intfc->curves; cur && *cur; cur++)
	{
		if (wave_type(*cur) == PHASE_BOUNDARY)
		{
			phsbdry = *cur;
			break;
		}
	}

	lri = ses_rs_grid_from_rho(ri,seos);
	for (bond = phsbdry->first; bond != NULL; bond = bond->next)
	{
		if (!Between(lri,Coords(bond->start)[0],
				    Coords(bond->end)[0]))
			continue;
		slope = (Coords(bond->end)[1] - Coords(bond->start)[1]) / 
			(Coords(bond->end)[0] - Coords(bond->start)[0]);
		S_grid = slope*(lri - Coords(bond->end)[0]) +
							Coords(bond->end)[1];
		state = left_state_at_point_on_curve(bond->start,bond,phsbdry);
		state1 = left_state_at_point_on_curve(bond->end,bond,phsbdry);
		press = ses_rs_redp(state)*
			ses_rs_rho_from_grid(Coords(bond->start)[0],seos)*
			ses_rs_temp_from_var(ses_rs_Tvar(state),seos) +
			ses_rs_coldp(state);
		press1 = ses_rs_redp(state1)*
			ses_rs_rho_from_grid(Coords(bond->end)[0],seos)*
			ses_rs_temp_from_var(ses_rs_Tvar(state1),seos) +
			ses_rs_coldp(state1);
		V  = 1.0 / ses_rs_rho_from_grid(Coords(bond->start)[0],seos);
		V1 = 1.0 / ses_rs_rho_from_grid(Coords(bond->end)[0],seos);
		dP = press1 - press;
		dV = V1 - V;
		coords[0] = lri;
		coords[1] = S_grid;
		set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);
		ses_solution(coords,COMP_PURE_PHASE,Hyper_surf(phsbdry),
			     NEGATIVE_SIDE,fr,(POINTER)wv,var);
		press = var[RS_RP]*ses_rs_rho_from_grid(lri,seos)*
			ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
		V = ses_rs_rho_from_grid(lri,seos);
		V = 1.0/V;
		adbgam = var[RS_AG];
	}

	/* See Menikoff and Plohr, The Riemann Problem for Fluid Flow
	*  of Real Materials, p. 101. 
	*/
	DEBUG_LEAVE(is_retrograde_bndry)
	return ((adbgam < -(V/press)*dP/dV) || (-(V/press)*dP/dV < 0.0)) ?
		YES : NO;

}		/*end is_retrograde_bndry*/


LOCAL boolean is_liquid_bndry(
	double		  x,
	LIQUID_BNDRY_TYPE type,
	double		  *tbls)
{
	boolean		is_lqd = NO;
	int		nr = (int)(tbls[2]);
	int		nt = (int)(tbls[3]);
	int		n4 = nr+nt+3*nr*nt+5;
	int		ntv = (int)(tbls[n4-1]);

        DEBUG_ENTER(is_liquid_bndry)

	switch (type)
	{
	case DENSITY_GIVEN:
		if(Between(x,tbls[n4 + 3*ntv],tbls[n4 + 4*ntv-1]))
			is_lqd = YES;
		break;

	case ENERGY_GIVEN:
		if(Between(x,tbls[n4 + 5*ntv],tbls[n4 +6*ntv-1]))
			is_lqd = YES;
		break;
	}

	DEBUG_LEAVE(is_liquid_bndry)
	return is_lqd;
}		/*end is_liquid_bndry*/

typedef struct {
	Locstate st00;
	double		r_retro;
} SONIC_INV_PARAMS;

LOCAL	void	make_composite(
	Locstate 	state,
	double		r_retro,
	double		p_phs,
	Locstate 	st0,
	Locstate 	st1,
	WAVE_CURVE 	*wave_cur)
{
	SONIC_INV_PARAMS sonic_inv_params;
	SESAME_EOS	*seos = Ses_params(state);
	Front		*fr = seos->fr[SESAME_RHO_ENTROPY];
	int		i;
	boolean 	on_pbdry;
	int		last;
	double		rho, drho, rho_grid, m2;
	double		rmin, rmax;
	double		a1s, a0s;
	double		epsilon, delta;
	static const double	SONEPS = 1.0e-05; /* TOLERANCE */
	static Locstate st00 = NULL;

	DEBUG_ENTER(make_composite)

	if (st00 == NULL)
	    (*Params(state)->_alloc_state)(&st00,Params(state)->sizest);

	drho = (fr->rect_grid->U[0] - fr->rect_grid->L[0])/NO_PTS_ON_COMP;

	rho_grid = ses_rs_grid_from_rho(r_retro,seos);

	on_pbdry = NO;
	last = NO;

	state_on_adiabat_with_pr(state,p_phs,st00,TGAS_STATE);
	for (i = 0; i < NO_PTS_ON_COMP; i++)
	{
	    rho_grid = rho_grid + i*drho;
	    rho = ses_rs_rho_from_grid(rho_grid,seos);
	    if (rho > Dens(st00))
	    {
	    	rho = Dens(st00);
	    	rho_grid = ses_rs_grid_from_rho(rho,seos);
	    	last = YES;
	    	if (Dens(st00) < Dens(state)) on_pbdry = YES;
	    }
 	    state_on_adiabat_with_dens(st00,rho,st0,TGAS_STATE);
	    if (sonic_shock(st0,st1,on_pbdry,r_retro,i) == NO)
	    {
	    	/* Have not found double sonic */
	    	wave_cur->ustart[i] = vel(0,st0);
	    	wave_cur->uend[i] = vel(0,st1);
	    	wave_cur->pstart[i] = pressure(st0);
	    	wave_cur->pend[i] = pressure(st1);
	    	wave_cur->rstart[i] = Dens(st0);
	    	wave_cur->rend[i] = Dens(st1);
	    	wave_cur->rhoc[i] = acoustic_impedance(st0);

	    }
	    else
	    {
		/* Passed double sonic */
	    	rmin = rho_grid - i*drho;
	    	rmin = ses_rs_rho_from_grid(rmin,seos);
	    	rmax = rho;
	    	last = YES;
			
	    	sonic_inv_params.st00  = st00;
	    	sonic_inv_params.r_retro  = r_retro;
	    	epsilon = SONEPS;
	    	delta = EPS*(rmax - rmin);
	    	if (bisection_find_root(dble_sonic_bisect,
	    		                (POINTER)&sonic_inv_params,0.0,&rho,
				        rmin,rmax,epsilon,delta) ==
		    FUNCTION_FAILED)
		{
		    (void) printf("WARNING in make_composite(), "
		                  "approximating double sonic\n");
		}
		else
		{
 		    state_on_adiabat_with_dens(st00,rho,st0,TGAS_STATE);
		    (void) sonic_shock(st0,st1,on_pbdry,r_retro,0);
		}
	    }
	    /* check to see if sonic on both sides */
	    a1s = acoustic_impedance_squared(st1);
	    a0s = acoustic_impedance_squared(st0);
	    /* This is the right way to check for sonic on both sides since
	    	the mass flux was set to rho0*c0 initially. */
	    m2 = a0s/a1s;
	    if ((fabs(m2 - 1.0) < SONEPS) || last)
	    {
	    	/* Sonic shock at both ends or last. You are done */
	    	if (DEBUG)
	    	    (void) printf("Double sonic\n");
	    	wave_cur->ustart[NO_PTS_ON_COMP-1] = vel(0,st0);
	    	wave_cur->uend[NO_PTS_ON_COMP-1] = vel(0,st1);
	    	wave_cur->pstart[NO_PTS_ON_COMP-1] = pressure(st0);
	    	wave_cur->pend[NO_PTS_ON_COMP-1] = pressure(st1);
	    	wave_cur->rstart[NO_PTS_ON_COMP-1] = Dens(st0);
	    	wave_cur->rend[NO_PTS_ON_COMP-1] = Dens(st1);
	       	wave_cur->rhoc[NO_PTS_ON_COMP-1] = acoustic_impedance(st0);
	    	break;
	    }
	}
	set_state(st0,TGAS_STATE,st0);
	set_state(st1,TGAS_STATE,st1);
	DEBUG_LEAVE(make_composite)
}		/*end make_composite*/

LOCAL boolean dble_sonic_bisect(
	double		rho,
	double		*m,
	POINTER 	params)
{
	SONIC_INV_PARAMS *prms = (SONIC_INV_PARAMS *)params;
	double		r_retro = prms->r_retro;
	Locstate 	st00 = prms->st00;
	double		a1s, a0s;
	boolean 	on_pbdry;
	static boolean 	first = YES;
	static Locstate st0 = NULL, st1 = NULL;

        DEBUG_ENTER(dble_sonic_bisect)
	if (first)
	{
		first = NO;
		(*Params(st00)->_alloc_state)(&st0,Params(st00)->sizest);
		(*Params(st00)->_alloc_state)(&st1,Params(st00)->sizest);
	}
	on_pbdry = NO;
 	state_on_adiabat_with_dens(st00,rho,st0,TGAS_STATE);
	if (fabs(rho - Dens(st00)) < EPS)
		on_pbdry = YES;
	if (sonic_shock(st0,st1,on_pbdry,r_retro,0) == NO)
	{
		a1s = acoustic_impedance_squared(st1);
		a0s = acoustic_impedance_squared(st0);
		*m = a0s/a1s - 1.0;
	}
	else
	{
		*m = -1.0;
	}
	DEBUG_LEAVE(dble_sonic_bisect)
	return FUNCTION_SUCCEEDED;
}		/*end dble_sonic_bisect*/


/*
* 
*			sonic_shock()
*
*	Determine a shock wave that is sonic to the given state0, so
*	m = rho0*c0. Returns state1.
*
*/

LOCAL boolean sonic_shock(
	Locstate 	state0,
	Locstate 	state1,
	boolean 	on_pbdry,
	double		r_retro,
	int		num)
{
	SESAME_EOS	*seos = Ses_params(state0);
	Front		*fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave		*wv = seos->wave[SESAME_RHO_ENTROPY];
	INTERFACE 	*intfc = fr->interf;
	CURVE		**cur, *phsbdry;
	COMPONENT 	comp;
	double		coords[MAXD];
	double		m, c, rho0, e0, p0, S0, var[NUM_SES_VAR];
	double		emin, emax, vmax;

        DEBUG_ENTER(sonic_shock)
	rho0 = Dens(state0);
	S0 = entropy(state0);
	coords[0] = ses_rs_grid_from_rho(rho0,seos);
	coords[1] = ses_rs_grid_from_entpy(S0,seos);

	set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);
	for(cur = intfc->curves;cur && *cur;cur++)
	{
	    if(wave_type(*cur) == PHASE_BOUNDARY)
	    {
	    	phsbdry = *cur;
	    	break;
	    }
	}
	if (on_pbdry == YES)
	{
	    comp = COMP_MIXED_PHASE;
	    ses_solution(coords,comp,Hyper_surf(phsbdry),POSITIVE_SIDE,
			 fr,(POINTER)wv,var);
	}
	else
	{
	    comp = nearest_interior_comp(YES,COMP_PURE_PHASE,
			                 coords,fr->interf);
	    ses_solution(coords,comp,NULL,POSITIVE_SIDE,fr,(POINTER)wv,var);
	}

	p0 = var[RS_RP]*rho0*
		ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
	c = var[RS_AG]*p0/ses_rs_rho_from_grid(rho0,seos);
	
	m = rho0*rho0*c;

	/* Must set limits on interval of solution first */
	e0 = specific_internal_energy(state0);
	if ((rho0 - r_retro) < EPS  || num == 0)
	{
	    coords[0] = ses_rs_grid_from_rho(r_retro,seos);
	    ses_solution(coords,COMP_PURE_PHASE,Hyper_surf(phsbdry),
			 NEGATIVE_SIDE,fr,(POINTER)wv,var);
	    emax = var[RS_RE]*ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CE];
	    vmax = 1.0/r_retro;
	    emin =  (2*m*e0 - p0*p0)/2*m;
	    Energy(state1) = emin -1.0;
	}
	else
	{
	    /* Set max energy to energy from previous try */
	    emax = Energy(state1);
	    vmax = 1.0/Dens(state1);
	}
	if (!num_Hug_by_mass_flux_squared(m,state0,state1,emin,emax,vmax))
	{
	    /* Found double sonic */
	    DEBUG_LEAVE(sonic_shock)
	    return YES;
	}
	else
	{
	    switch (state_type(state0))
	    {
	    case EGAS_STATE:
	    	break;
	    case TGAS_STATE:
	    	Press(state1) = pressure(state1);
	    	break;
	    case GAS_STATE:
	    	Energy(state1) *= Dens(state1);
	    	break;
	    case VGAS_STATE:
	    	set_state(state1,VGAS_STATE,state1);
	    	break;
	    default:
	    	screen("ERROR in sonic_shock(), "
	    	       "Unknown state type %d\n",state_type(state1));
	    	clean_up(ERROR);
	    }
	    DEBUG_LEAVE(sonic_shock)
	    return NO;
	}
}		/*end sonic_shock*/

EXPORT void make_wave_crv(
	int		l_or_r,
	Locstate 	state,
	double		*pi,
	double		*ri,
	double		*ui,
	int		ni,
	WAVE_CURVE 	*wave_cur)
{
	double		tmp;
	double		p_phs;
	int		i, j;
	static Locstate st0 = NULL, st1 = NULL;

        DEBUG_ENTER(make_wave_crv)
	if (st0 == NULL)
	{
	    (*Params(state)->_alloc_state)(&st0,Params(state)->sizest);
	    (*Params(state)->_alloc_state)(&st1,Params(state)->sizest);
	}
	
	/* Curves are just shocks and rarefactions */

	set_state(wave_cur->st[0],TGAS_STATE,state);
	wave_cur->w_type[0] = RAREFACTION;
	set_state(wave_cur->st[1],TGAS_STATE,state);
	wave_cur->w_type[1] = SHOCK;
	wave_cur->special = NO;
	
	/* No special behavior */
	if (ni == 0)
	{
	    DEBUG_LEAVE(make_wave_crv)
	    return;
	}

	if (DEBUG)
	    (void) printf("Special behavior\n");

	
	/* Initialize wave curve through given state */

	/* Order intersections by pressure */
	if (ni == 2)
	{
	    if ( pi[0] > pi[1])
	    {
	        tmp = pi[0];
	        pi[0] = pi[1];
	        pi[1] = tmp;
	        tmp = ui[0];
	        ui[0] = ui[1];
	        ui[1] = tmp;
	        tmp = ri[0];
	        ri[0] = ri[1];
	        ri[1] = tmp;
	    }
	}

	j = 0;
	wave_cur->num_waves = 0;
	for (i = 0; i < ni; i++)
	{
	    if (j > MAX_NUM_WAVES)
	    {
	    	screen("ERROR in make_wave_cur(), "
	    	       "too many wave curves\n");
	    	clean_up(ERROR);
	    }

	    if (pi[i] < pressure(state))
	    {
	    	/* On rarefaction branch */
	    	if (is_retrograde_bndry(ri[i],state))
	    	{ 
	    	    if (DEBUG)
	    		(void) printf("Retrograde behavior\n");
		    if (ni == 2 && pi[1] < pressure(state))
		    {
			p_phs = pi[1];
		    }
		    else
		    {
			p_phs = pressure(state);
		    }

		   /* Compute composite part of wave curve. This 
		    * consists of a shock that rarefies the fluid */

		    if (DEBUG)
			(void) printf("Making composite\n");
		    make_composite(state,ri[i],p_phs,st0,st1,wave_cur);
			    
		    wave_cur->special = YES;
		    set_state(wave_cur->st[j],TGAS_STATE,st1);
		    wave_cur->w_type[j] = COMPOSITE;
		    wave_cur->num_waves++;
		    j++;
		    Vel(wave_cur->st[j])[0] = ui[i];
		    Press(wave_cur->st[j]) = pi[i];
		    Dens(wave_cur->st[j]) = ri[i];
		    Set_params(wave_cur->st[j],state);
		    set_type_of_state(wave_cur->st[j],TGAS_STATE);
		    wave_cur->w_type[j] = COMPOSITE;
		    wave_cur->num_waves++;
		    j++;
		    set_state(wave_cur->st[j],TGAS_STATE,state);
		    wave_cur->w_type[j] = RAREFACTION;
		    wave_cur->num_waves++;
		    j++;
		}
		else
		{
		    if (DEBUG)
			(void) printf("Setting splti raref\n");
		    wave_cur->special = YES;
		    Vel(wave_cur->st[j])[0] = ui[i];
		    Press(wave_cur->st[j]) = pi[i];
		    Dens(wave_cur->st[j]) = ri[i];
		    Set_params(wave_cur->st[j],state);
		    set_type_of_state(wave_cur->st[j],TGAS_STATE);
		    wave_cur->w_type[j] = RAREFACTION;
		    wave_cur->num_waves++;
		    j++;
		    set_state(wave_cur->st[j],TGAS_STATE,state);
		    wave_cur->w_type[j] = RAREFACTION;
		    wave_cur->num_waves++;
		    j++;
		}
		set_state(wave_cur->st[j],TGAS_STATE,state);
		wave_cur->w_type[j] = SHOCK;
		wave_cur->num_waves++;
			
	    }
	    else
	    {
	    	/* On shock  branch */
	    	if (j < 1)
		    j = 1;
	    	if (is_retrograde_bndry(ri[i],state))
		{
		/* Shock splitting */
		    wave_cur->special = YES;
		    set_state(wave_cur->st[j],TGAS_STATE,state);
		    wave_cur->w_type[j] = SHOCK;
		    wave_cur->num_waves++;
		    j++;
		    Vel(wave_cur->st[j])[0] = ui[i];
		    Press(wave_cur->st[j]) = pi[i];
		    Dens(wave_cur->st[j]) = ri[i];
		    Set_params(wave_cur->st[j],state);
		    set_type_of_state(wave_cur->st[j],TGAS_STATE);
		    wave_cur->w_type[j] = SHOCK;
		    wave_cur->num_waves++;
		    j++;
		    get_Hug_intersection(wave_cur->st[j-1],
					 wave_cur->st[j],l_or_r,wave_cur);
		    j++;
		}	
		else
		{
			/* No splitting */
		    set_state(wave_cur->st[j],TGAS_STATE,state);
		    wave_cur->w_type[j] = SHOCK;
		    wave_cur->num_waves++;
		    j++;
		}
	    }
	}
	/* At the very least one will have a rarefaction and a shock */
	if (wave_cur->num_waves < 2)
	    wave_cur->num_waves = 2;
	if (DEBUG)
	{
	    for (j = 0; j < MAX_NUM_WAVES; j++)
	    {
		(void) printf("Wave type for wave curve is "
		              "%d\n",wave_cur->w_type[j]);
		(void) printf("%g\n",Vel(wave_cur->st[j])[0]);
		(void) printf("%d\n",wave_cur->num_waves);
	    }
	}
        DEBUG_LEAVE(make_wave_crv)
}		/*end make_wave_crv*/

typedef struct {
	double		p0, u0, V0, p1, u1, V1;
	int		l_or_r;
} SPLIT_INV_PARAMS;

/*
*			get_Hug_intersection()
*
*	Get the intersection of the Hugoniot H1 and H0 to determine region in 
*	which split shocks occur. 
*
*/

LOCAL	void get_Hug_intersection(
	Locstate 	st0,
	Locstate 	st1,
	int		l_or_r,
	WAVE_CURVE 	*wave_cur)
{
	SPLIT_INV_PARAMS split_inv_params;
	double		V1, V0, p1, p0; 
	double		u0, pmax, pmin, epsilon, delta;
	double		p, u, rho;
	int		num_waves;
	
        DEBUG_ENTER(get_Hug_intersection)
	V1 = Dens(st1);
	V1 = 1.0/V1;
	V0 = Dens(st0);
	V0 = 1.0/V0;
	p1 = pressure(st1);	
	p0 = pressure(st0);	

/* Set pmax == P  when V = 0 and curves cross in P-V plane */

	pmax = (p1 + p0)*(V1 - V0) + p0*V0 - p1*V1;
	pmax = pmax/(V1 - V0);

/* Set pmin == P when V = .9*V1 and curves corss in P-V plane to ensure
  it will not converge to p0 = p1 and V0 = V1 */

	pmin = (p1 + p0)*(V1 - V0) + p0*(V0 -.9*V1) -p1*.1*V1;

/* Iterate on values of p so u1 + (p - p1)/m1 = u0 + (p - p0)/m0 */

	split_inv_params.V1  = V1;
	split_inv_params.V0  = V0;
	split_inv_params.p1  = p1;
	split_inv_params.p0  = p0;
	split_inv_params.u1  = vel(0,st1);
	split_inv_params.l_or_r = l_or_r;
	u0 = split_inv_params.u0  = vel(0,st0);
	epsilon = EPS*u0;
	delta = EPS*(pmin - pmax);
	if (bisection_find_root(ses_hug_bisect,
				(POINTER)&split_inv_params,0.0,&p,
				pmin,pmax,epsilon,delta) == FUNCTION_FAILED)
	{
		screen("ERROR in get_Hug_intst(), ");
		screen("bisection does not converge\n");
		clean_up(ERROR);
	}
	else
	{
		
		num_waves = wave_cur->num_waves;
		
		rho = (p1 + p0)*(V0 - V1) - (p + p1)*V1 + (p + p0)*V0;
		rho = (p0 - p1)/rho;
		if (l_or_r == LEFT_FAMILY)
			u = u0 - (p - p0)/sqrt((p-p0)/V0 - 1.0/rho);
		else
			u = u0 + (p - p0)/sqrt((p-p0)/V0 - 1.0/rho);
	
		Vel(wave_cur->st[num_waves])[0] = u;
		Press(wave_cur->st[num_waves]) = p;
		Dens(wave_cur->st[num_waves]) = rho;
		Set_params(wave_cur->st[num_waves],st0);
		set_type_of_state(wave_cur->st[num_waves],TGAS_STATE);
		wave_cur->w_type[num_waves] = SHOCK;
		wave_cur->num_waves++;
	}
        DEBUG_LEAVE(get_Hug_intersection)
}		/*end get_Hug_intersection*/

LOCAL boolean ses_hug_bisect(
	double		p,
	double		*du,
	POINTER 	params)
{
	SPLIT_INV_PARAMS *prms = (SPLIT_INV_PARAMS *)params;
	double		u, V;
	double		V0 = prms->V0;
	double		V1 = prms->V1;
	double		p0 = prms->p0;
	double		p1 = prms->p1;
	double		u0 = prms->u0;
	double		u1 = prms->u1;
	int		l_or_r = prms->l_or_r;

        DEBUG_ENTER(ses_hug_bisect)
/* TODO: This may not be robust */
/* Get V from intersection in P_V plane */

	V = (p1 + p0)*(V0 - V1) - (p + p1)*V1 + (p + p0)*V0;
	V = V /(p0 - p1);
	
/* Check intersection in P-u plane */

	if (l_or_r == LEFT_FAMILY)
	{
		u = u0 - (p - p0)/sqrt((p-p0)/V0 - V);
		*du = u - u1 + (p - p1)/sqrt((p-p1)/(V1 - V));
	}
	else
	{
		u = u0 + (p - p0)/sqrt((p-p0)/V0 - V);
		*du = u - u1 - (p - p1)/sqrt((p-p1)/(V1 - V));
	}
	DEBUG_LEAVE(ses_hug_bisect)
	return FUNCTION_SUCCEEDED;
}		/*end ses_hug_bisect*/


EXPORT void state_on_comp(
	double		p,
	double		*rstart,
	double		*pstart,
	double		*ustart,
	double		*rhoc,
	WAVE_CURVE 	*wave_cur)
{
	int		i, num;
	double		slope;
	double		pend;
	double		pend1, pstart1, rstart1, ustart1, rhoc1;

        DEBUG_ENTER(state_on_comp)
	num = wave_cur->num_waves;
	num = num - 1;
	for (i = 0; i < NO_PTS_ON_COMP - 1; i++)
	{
		pend = wave_cur->pend[i];
		pend1 = wave_cur->pend[i+1];
		if (Between(p,pend1,pend))
		{
			*rhoc = wave_cur->rhoc[i];
			rhoc1 = wave_cur->rhoc[i+1];
			*pstart = wave_cur->pstart[i];
			pstart1 = wave_cur->pstart[i+1];
			*rstart = wave_cur->rstart[i];
			rstart1 = wave_cur->rstart[i+1];
			*ustart = wave_cur->ustart[i];
			ustart1 = wave_cur->ustart[i+1];
			slope = (p - pend)/(pend1 - pend);
			*rhoc = (rhoc1 - *rhoc)*slope + *rhoc;
			*pstart = (pstart1 - *pstart)*slope + *pstart;
			*ustart = (ustart1 - *ustart)*slope + *ustart;
			*rstart = (rstart1 - *rstart)*slope + *rstart;
			break;
		}
	}
        DEBUG_LEAVE(state_on_comp)
}		/*end state_on_comp*/

EXPORT void get_ph_sound_spd(
	double		*cl,
	double		*cr,
	Locstate 	state)
{
	SESAME_EOS	*seos = Ses_params(state);
	Front		*fr = seos->fr[SESAME_RHO_ENTROPY];
	Wave		*wv = seos->wave[SESAME_RHO_ENTROPY];
	INTERFACE 	*intfc = fr->interf;
	CURVE		**cur, *phsbdry;
	double		coords[MAXD];
	double		rho, S, p;
	double		var[NUM_SES_VAR];

        DEBUG_ENTER(get_ph_sound_spd)
	rho = Dens(state);
	S = entropy(state);

	coords[0] = ses_rs_grid_from_rho(rho,seos);
	coords[1] = ses_rs_grid_from_entpy(S,seos);

	for(cur = intfc->curves;cur && *cur;cur++)
	{
		if(wave_type(*cur) == PHASE_BOUNDARY)
		{
			phsbdry = *cur;
			break;
		}
	}

	set_ses_intrp_flag_all(SESAME_RHO_ENTROPY);
	ses_solution(coords,COMP_PURE_PHASE,Hyper_surf(phsbdry),NEGATIVE_SIDE,
		     fr,wv,var);
	p = var[RS_RP]*rho*
		ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
	*cl = sqrt(var[RS_AG]*p/rho);
	
	ses_solution(coords,COMP_MIXED_PHASE,Hyper_surf(phsbdry),
		     POSITIVE_SIDE,fr,wv,var);
	p = var[RS_RP]*rho*
		ses_rs_temp_from_var(var[RS_T],seos) + var[RS_CP];
	*cr = sqrt(var[RS_AG]*p/rho);
        DEBUG_LEAVE(get_ph_sound_spd)
}		/*end get_ph_sound_spd*/

LOCAL	boolean num_Hug_by_mass_flux_squared(
	double		mf_sqr,
	Locstate	st0,
	Locstate	ans,
	double		emin,
	double		emax,
	double		vmax0)
{
	int		n, unstable = 0;
	double		p0;
	double		e0, emax0;
	double		rho0, V0, V;
	double		fmin, fmax, fmax0, f, dfde;
	double		x, b, c, twomf_sqr;
	double		elast, flast;
	static const int MAX_NUM_ITER_NEWTON = 20;/*TOLERANCE*/
	static const int MAX_NUM_ITER_BISECTION = 40;/*TOLERANCE*/
	static const int UNSTABLE = 2;/*TOLERANCE*/

	DEBUG_ENTER(num_Hug_by_mass_flux_squared)
	Set_params(ans,st0);
	set_type_of_state(ans,EGAS_STATE);
	p0 = pressure(st0);
	e0 = specific_internal_energy(st0);
	rho0 = Dens(st0);
	V0 = 1.0/rho0;
	twomf_sqr = 2.0*mf_sqr;
	c = p0*p0 - twomf_sqr*e0;
	emax0 = emax + p0*vmax0 - 0.5*mf_sqr*vmax0*vmax0;
	b = p0 + mf_sqr*V0;
	fmin = mf_sqr - acoustic_impedance_squared(st0);
	fmax = fmax0 = -HUGE_VAL;
	flast = fmin;
	elast = emin;

	if ( Energy(ans) < emin || Energy(ans) > emax)
		Energy(ans) = 0.5*(emin+emax);
	if (DEBUG)
	{
		verbose_print_state("state 0",st0);
		(void) printf("Before loop\n");
		(void) printf("emin = %g, emax = %g, fmin = %g, fmax = %g\n",
			emin,emax,fmin,fmax);
		(void) printf("First value of e = %g\n",Energy(ans));
	}

	for ( n = 0; n < MAX_NUM_ITER_NEWTON; n++)
	{
		V = 2.0*(emax0 - Energy(ans))/
			(b + sqrt(c + twomf_sqr*Energy(ans)));
		Dens(ans) = 1.0/V;
		f = mf_sqr - (pressure(ans) - p0)/(V0 - V);

		if (DEBUG)
			(void) printf("rho = %g, e = %g, f = %g\n",Dens(ans),
			Energy(ans),f);

		if (fabs(f) < EPSILON*mf_sqr)
		{
			DEBUG_LEAVE(num_Hug_by_mass_flux_squared)
			return FUNCTION_SUCCEEDED;
		}

		if( f > 0.0 )
		{
			fmin = f;
			emin = Energy(ans);
		}
		else
		{
			fmax = f;
			emax = Energy(ans);
		}

		dfde = (f-flast)/(Energy(ans)-elast);
		flast = f;
		elast = Energy(ans);

		if (DEBUG)
			(void) printf("dfde = %g\n",dfde);
		Energy(ans) -= f/dfde;
		if( Energy(ans) < emin)
		{
			Energy(ans) = emin + 0.1*(emax-emin);
			if (unstable++ > UNSTABLE) break;
		}
		else if (Energy(ans) > emax)
		{
			Energy(ans) = emax - 0.1*(emax-emin);
			if (unstable++ > UNSTABLE) break;
		}
	}
		
	if (DEBUG)
		(void) printf("Using bisection, unstable = %s\n",
		(unstable > UNSTABLE) ? "YES" : "NO");

		/* Bisection method for robust alternative */

	for (n = 0; n < MAX_NUM_ITER_BISECTION; n++)
	{
		x = (fmax == fmax0) ? 0.5 : fmin/(fmin-fmax);
		if (x < 0.1) x = 0.2;
		else if(x > 0.9) x = 0.8;

		Energy(ans) = emin + x*(emax - emin);
		V = 2.0*(emax0 - Energy(ans))/
				(b + sqrt(c + twomf_sqr*Energy(ans)));
		Dens(ans) = 1.0/V;
		f = mf_sqr - (pressure(ans) - p0)/(V0 - V);
		if (DEBUG)
		{
			(void) printf("\nemin = %g, emax = %g, fmin = %g, fmax = %g\n",
				emin,emax,fmin,fmax);
			(void) printf("rho = %g, e = %g, f = %g\n",
				Dens(ans),Energy(ans),f);
		}
		if (fabs(f) < EPSILON*mf_sqr)
		{
			DEBUG_LEAVE(num_Hug_by_mass_flux_squared)
			return FUNCTION_SUCCEEDED;
		}

		if( f > 0.0 )
		{
			fmin = f;
			emin = Energy(ans);
		}
		else
		{
			fmax = f;
			emax = Energy(ans);
		}
	}
	(void) printf("WARNING in num_Hug_by_mass_flux(), ");
	(void) printf("No convergence to solution\n");
	DEBUG_LEAVE(num_Hug_by_mass_flux_squared)
	return FUNCTION_FAILED;
}		/*end num_Hug_by_mass_flux_squared*/

#endif /* defined(TWOD) && defined(PHASE_CODE) */

