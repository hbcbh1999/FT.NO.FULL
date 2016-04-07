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
*				sesinv.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if defined(SESAME_CODE) && defined(TWOD)
#define DEBUG_STRING    "ses_hyp"
#include <geos/sesame.h>

	/*LOCAL Function Prototypes*/
#if defined(__cplusplus)
extern "C" {    
#endif /* defined(__cplusplus) */
    LOCAL   LSODE_FUNC  dSdxi;
    LOCAL   LSODE_FUNC  dSdeta;
    LOCAL   LSODE_FUNC  integrate_entropy;
#if defined(__cplusplus)
}                          
#endif /* defined(__cplusplus) */
LOCAL	void	coldpe(double,SESAME_EOS*,double*,double*);
LOCAL	void	set_entropy_at_T_rho(double,double,Locstate,SESAME_EOS*);
LOCAL	void	ses_rt_set_entropy_on_isochore(CURVE*,int,int,
					       double,double,double**);
LOCAL	void	ses_rt_set_entropy_on_isotherm(CURVE*,int,int,
					       double,double,double**);

LOCAL	SESAME_EOS *lseos = NULL;

EXPORT	void	s2eos_lookup(
	double		rho,
	double		T,
	double		*p,
	double		*e,
	SESAME_EOS	*seos)
{
	double	*tbls = seos->sestab.tbls;
	static	int	ir = 1;

	FORTRAN_NAME(s2eos)(&ir,tbls,&rho,&T,p,e);

	if (T > 0.0)
	{
	    double	p2, e1;
	    p2 = p[2];
	    e1 = e[1];
	    p[2] = 0.5*(p2 + (p[0] - sqr(rho)*e1)/T);
	    e[1] = 0.5*(e1 + (p[0] - T*p2)/sqr(rho));
	}
}		/*end s2eos_lookup*/

EXPORT	void setrt(
	double		rho,
	double		T,
	Locstate	state,
	SESAME_EOS	*seos)
{
	double	ABS_SES_EPS = seos->ABS_SES_EPS;
	double	p[3], e[3], cold_p, cold_e;
	double	csq;
	double	x;


	/* Compute pressure, energy, adiabatic and Gruneisen gammas */

	s2eos_lookup(rho,T,p,e,seos);
	csq = sesame_rt_sound_speed_squared(rho,T,p,e);
	ses_rt_adb_gam(state) = rho*csq/p[0];
	ses_rt_gru_gam(state) = p[2]/(e[2]*rho);
	x = fabs(ses_rt_gru_gam(state)*
	         dlog_rho_drho_grid(rho,seos)*dT_grid_dlogT(T,seos));
	Gru_gam_abs_max(seos) = max(Gru_gam_abs_max(seos),x);

	coldpe(rho,seos,&cold_p,&cold_e);
	ses_rt_coldp(state) = cold_p;
	ses_rt_colde(state) = cold_e;
	if (T > ABS_SES_EPS)
	{
		ses_rt_redp(state) = (p[0] - cold_p)/(rho*T);
		ses_rt_rede(state) = (e[0] - cold_e)/T;
	}
	else
	{
		ses_rt_redp(state) = p[2]/rho;
		ses_rt_rede(state) = e[2];
	}
	Pressure_min(seos) = min(p[0],Pressure_min(seos));
	Pressure_max(seos) = max(p[0],Pressure_max(seos));
	Energy_min(seos) = min(e[0],Energy_min(seos));
	Energy_max(seos) = max(e[0],Energy_max(seos));

	ses_rt_S(state) = 0.0;		/*Not set yet*/
	ses_rt_riv(state) = 0.0;	/*Not set yet*/
}		/*end setrt*/

EXPORT	double	sesame_rt_sound_speed_squared(
	double	rho,
	double	T,
	double	*p,
	double	*e)
{
	double	csq;

	csq = p[1];
	if (fabs(e[2]) > fabs(p[2]) * 1e-4) /*TOLERANCE*/
		csq += (p[0] / (rho * rho) - e[1]) * p[2] / e[2];
	if (csq < 0.0)
	{
	    if (DEBUG)
	    {
		(void) printf("WARNING in sesame_rt_sound_speed_squared(), "
		              "negatve sound speed found at rho = %g, T = %g"
		              "\tp = %g, e = %g\n",rho,T,p[0],e[0]);
		(void) printf("\tdP/drho = %g, de/drho = %g\n",p[1],e[1]);
		(void) printf("\tdP/dT = %g, de/dT = %g\n",p[2],e[2]);
	    }
	    csq = 0.0;
	}
	return csq;
}		/*end sesame_rt_sound_speed_squared*/

LOCAL	void	coldpe(
	double		rho,
	SESAME_EOS	*seos,
	double		*cold_p,
	double		*cold_e)
{
	double	*tbls = seos->sestab.tbls;
	double	T;
	double	p[3], e[3];
	int	nr;

	nr = irint(tbls[2]);
	T = tbls[4+nr];
	s2eos_lookup(rho,T,p,e,seos);
	*cold_p = p[0];
	*cold_e = e[0];
}		/*end coldpe*/

/*
*			set_RT_entropy_from_cold_curve():
*
*	Computes the entropy values for the density-temperature sesame table.
*/

/*ARGSUSED*/
EXPORT	void	set_RT_entropy_from_cold_curve(
	Front		*fr,
	Wave		*wave,
	SESAME_EOS	*seos)
{
	BDRY_SIDE  side;
	RECT_GRID  *gr = fr->rect_grid;
	CURVE	   **c;
	INTERFACE  *intfc = fr->interf;
	Locstate   state;
	double	   xi;
	double	   T_min;
	double	   deta;/*Temperature spacing for entropy array*/
	double	   dxi; /*Density spacing for entropy array*/
	double	   **Sarray;
	double      S[2];
	double	   aerr = seos->abser0;
	double	   rerr = seos->reler0;
	double	   terr = seos->terr0;
	double	   eta0_T;
	int	   icoords[MAXD];
	int	   i, j;
	int	   ix, iy;
	int	   ieta0_T, neta_T;
	int	   xmax = wave->rect_grid->gmax[0];
	int	   ymax = wave->rect_grid->gmax[1];
	int	   neta;/*No. temp. points for entropy array*/
	int	   nxi; /*No. density points for entropy array*/
	int        istate;
	static int one = 1;

	lseos = seos;
	T_min = ses_rt_temp_from_grid(gr->L[1],seos);
	dxi = 0.5*gr->h[0];
	deta = 0.5*gr->h[1];
	nxi = 2*xmax + 1;
	neta = 2*ymax + 1;
	bi_array(&Sarray,nxi,neta,FLOAT);

	if (zero_temperature_cold_curve(seos) == YES)
	{
	    /*Minimum temperature is assumed to be zero*/
	    /*Set entropy on cold curve to zero*/

	    double	Sp;

	    for (i = 0; i < nxi; ++i)
	    	Sarray[i][0] = 0.0;

	    /*Use backward Euler to set entropy on second row*/

	    S[0] = S[1] = 0.0;
	    xi = gr->L[0];
	    eta0_T = gr->L[1] + deta;
	    for (i = 0; i < nxi; ++i)
	    {
	        seos->de_params.rho = ses_rt_rho_from_grid(xi,seos);
	        dSdeta(&one,&eta0_T,S,&Sp);
	        Sarray[i][1] = deta*Sp;
	        xi += dxi;
	    }
	    ieta0_T = 1;
	    neta_T = neta - 1;
	}
	else
	{
	    double	*S_T;
	    /*Integrate entropy along minimum temperature isotherm*/
	    uni_array(&S_T,nxi,FLOAT);
	    istate = 1;
	    S_T[0] = 0.0;
	    seos->de_params.temperature = T_min;
	    if (!ode_solver(dSdxi,S_T,gr->L[0],dxi,nxi,rerr,aerr,
			    terr,&istate,NULL))
	    {
		screen("ERROR in set_RT_entropy_from_cold_curve(), "
		       "ode_solver() failed for dSdxi\n");
		clean_up(ERROR);
	    }
	    for (i = 0; i < nxi; ++i)
	    	Sarray[i][0] = S_T[i];
	    ieta0_T = 0;
	    neta_T = neta;
	    eta0_T = gr->L[1];
	    free(S_T);
	}

	xi = gr->L[0];
	for (i = 0; i < nxi; ++i)
	{
	    istate = 1;
	    seos->de_params.rho = ses_rt_rho_from_grid(xi,seos);
	    if (!ode_solver(dSdeta,Sarray[i]+ieta0_T,eta0_T,deta,neta_T,
	    	            rerr,aerr,terr,&istate,NULL))
	    {
		screen("ERROR in set_RT_entropy_from_cold_curve(), "
		       "ode_solver() failed for dSdeta\n");
		clean_up(ERROR);
	    }
	    xi += dxi;
	}

	for (i = 0; i < nxi; ++i)
	for (j = 0; j < neta; ++j)
	{
	    Entropy_min(seos) = min(Sarray[i][j],Entropy_min(seos));
	    Entropy_max(seos) = max(Sarray[i][j],Entropy_max(seos));
	}

	/*Set interior states*/
	for (iy = 0; iy < ymax; ++iy)
	{
	    icoords[1] = iy;
	    for (ix = 0; ix < xmax; ++ix)
	    {
	    	icoords[0] = ix;
	    	state = Rect_state(icoords,wave);
	    	ses_rt_S(state) = Sarray[2*ix+1][2*iy+1];
	    }
	}

	/*Set front states*/
	for (c = intfc->curves; c && *c; ++c)
	{
	    side = rect_bdry_side_for_curve(NULL,NULL,*c,gr);
	    switch (side)
	    {
	    case LEFT_BDRY:
		ses_rt_set_entropy_on_isochore(*c,neta-2,0,deta,
					       gr->L[1],Sarray);
	    	break;
	    case RIGHT_BDRY:
		ses_rt_set_entropy_on_isochore(*c,neta-2,nxi-1,deta,
					       gr->L[1],Sarray);
	    	break;
	    case LOWER_BDRY:
		ses_rt_set_entropy_on_isotherm(*c,nxi-2,0,dxi,
					       gr->L[0],Sarray);
	    	break;
	    case UPPER_BDRY:
		ses_rt_set_entropy_on_isotherm(*c,nxi-2,neta-1,dxi,
					       gr->L[0],Sarray);
	    	break;
	    default:
		screen("ERROR in set_RT_entropy_from_cold_curve(), "
		       "unsupported boundary side %d\n",side);
		clean_up(ERROR);
	    }
	}
	free(Sarray);
}		/*end set_RT_entropy_from_cold_curve*/

EXPORT	boolean	zero_temperature_cold_curve(
	SESAME_EOS	*seos)
{
	double		T_min = Temp_min(seos);
	double		T_ref = Temp_ref(seos);

	return (T_min < 1.0e-6*T_ref) ? YES : NO; /*TOLERANCE*/
}		/*end zero_temperature_cold_curve*/

LOCAL	void	ses_rt_set_entropy_on_isotherm(
	CURVE		*c,
	int		imax,
	int		j,
	double		dxi,
	double		xi_min,
	double		**Sarray)
{
	BOND		*b;
	POINT		*p;
	double		xi, a;
	int		i;
	Locstate	sl, sr;

	b = c->first;
	p = b->start;
	xi = Coords(p)[0];
	i = (int)((xi - xi_min)/dxi);
	i = max(i,0);
	i = min(i,imax);
	a = (xi - xi_min)/dxi - i;
	a = max(a,0.0);
	a = min(a,1.0);
	sl = left_state_at_point_on_curve(p,b,c);
	sr = right_state_at_point_on_curve(p,b,c);
	ses_rt_S(sl) = ses_rt_S(sr) = (1.0 - a)*Sarray[i][j] + a*Sarray[i+1][j];
	for (b = c->first; b != NULL; b = b->next)
	{
	    p = b->end;
	    xi = Coords(p)[0];
	    i = (int)((xi - xi_min)/dxi);
	    i = max(i,0);
	    i = min(i,imax);
	    a = (xi - xi_min)/dxi - i;
	    a = max(a,0.0);
	    a = min(a,1.0);
	    sl = left_state_at_point_on_curve(p,b,c);
	    sr = right_state_at_point_on_curve(p,b,c);
	    ses_rt_S(sl) = ses_rt_S(sr) =
		        (1.0 - a)*Sarray[i][j] + a*Sarray[i+1][j];
	}
}		/*end ses_rt_set_entropy_on_isotherm*/

LOCAL	void	ses_rt_set_entropy_on_isochore(
	CURVE		*c,
	int		jmax,
	int		i,
	double		deta,
	double		eta_min,
	double		**Sarray)
{
	BOND		*b;
	POINT		*p;
	double		eta, a;
	int		j;
	Locstate	sl, sr;

	b = c->first;
	p = b->start;
	eta = Coords(p)[1];
	j = (int)((eta - eta_min)/deta);
	j = max(j,0);
	j = min(j,jmax);
	a = (eta - eta_min)/deta - j;
	a = max(a,0.0);
	a = min(a,1.0);
	sl = left_state_at_point_on_curve(p,b,c);
	sr = right_state_at_point_on_curve(p,b,c);
	ses_rt_S(sl) = ses_rt_S(sr) = (1.0 - a)*Sarray[i][j] + a*Sarray[i][j+1];
	for (b = c->first; b != NULL; b = b->next)
	{
	    p = b->end;
	    eta = Coords(p)[1];
	    j = (int)((eta - eta_min)/deta);
	    j = max(j,0);
	    j = min(j,jmax);
	    a = (eta - eta_min)/deta - i;
	    a = max(a,0.0);
	    a = min(a,1.0);
	    sl = left_state_at_point_on_curve(p,b,c);
	    sr = right_state_at_point_on_curve(p,b,c);
	    ses_rt_S(sl) = ses_rt_S(sr) =
		        (1.0 - a)*Sarray[i][j] + a*Sarray[i][j+1];
	}
}		/*end ses_rt_set_entropy_on_isochore*/


/*
*			set_RT_entropy_from_mid_point():
*
*	Computes the entropy values for the density-temperature sesame table.
*/

/*ARGSUSED*/
EXPORT	void	set_RT_entropy_from_mid_point(
	Front		*fr,
	Wave		*wave,
	SESAME_EOS	*seos)
{
	POINT			*pt;
	HYPER_SURF		*hs;
	HYPER_SURF_ELEMENT	*hse;
	INTERFACE		*intfc = fr->interf;
	Locstate		sl, sr, state;
	double			*coords;
	double			rho, T;
	int			icoords[MAXD];
	int			ix, iy;
	int			xmax = wave->rect_grid->gmax[0];
	int			ymax = wave->rect_grid->gmax[1];

	lseos = seos;

	/*Set front states*/
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&pt,&hse,&hs) != NO)
	{
	    slsr(pt,hse,hs,&sl,&sr);
	    rho = ses_rt_rho_from_grid(Coords(pt)[0],seos);
	    T = ses_rt_temp_from_grid(Coords(pt)[1],seos);
	    set_entropy_at_T_rho(rho,T,sl,seos);
	    ses_rt_S(sr) = ses_rt_S(sl);

	}

	/*Set interior states*/
	for (iy = 0; iy < ymax; ++iy)
	{
		icoords[1] = iy;
		for (ix = 0; ix < xmax; ++ix)
		{
			icoords[0] = ix;
			coords = Rect_coords(icoords,wave);
			state = Rect_state(icoords,wave);
	    		rho = ses_rt_rho_from_grid(coords[0],seos);
	    		T = ses_rt_temp_from_grid(coords[1],seos);
	    		set_entropy_at_T_rho(rho,T,state,seos);
		}
	}
}		/*end set_RT_entropy_from_mid_point*/

LOCAL	void	set_entropy_at_T_rho(
	double		rho,
	double		T,
	Locstate	state,
	SESAME_EOS	*seos)
{
	double	eta, xi, D;
	double	S[2];
	static const double DMIN = 1.0e-06;	/*TOLERANCE*/

	eta = ses_rt_grid_from_temp(T,seos) - seos->de_params.tlgmid;
	xi = ses_rt_grid_from_rho(rho,seos) - seos->de_params.rlgmid;
	D = hypot(eta,xi);
	S[0] = 0.0;
	if (D < DMIN)
	    S[1] = 0.0;
	else
	{
	    double	aerr = seos->abser0;
	    double	rerr = seos->reler0;
	    double	terr = seos->terr0;
	    int	istate = 1;

	    seos->de_params.eta = eta/D;
	    seos->de_params.xi = xi/D;
	    if (!ode_solver(integrate_entropy,S,0.0,D,2,
			    rerr,aerr,terr,&istate,NULL))
	    {
		screen("ERROR in set_entropy_at_T_rho(), "
		       "ode_solver() failed\n");
		clean_up(ERROR);
	    }
	}
	ses_rt_S(state) = S[1];
	Entropy_min(seos) = min(S[1],Entropy_min(seos));
	Entropy_max(seos) = max(S[1],Entropy_max(seos));
}		/*end set_entropy_at_T_rho*/

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */

/*
*			dSdxi():
*
*	Derivative of entropy along an isotherm.  Density is given
*	as a function of the variable xi.
*
*	This function is in a form suitable for use by the ode solver
*	lsode.  Since lsode is written in Fortran the arguments are
*	passed as pointers.
*
*	Input:
*		*neg = number of equations (unused, assumed to be 1)
*		*xi = ses_rt_grid_from_rho(rho,seos)
*		*y = entropy value, (unused)
*	Output:
*		*yp = dS/dxi
*/

/*ARGSUSED*/
LOCAL	void dSdxi(
	int		*neq,
	double		*xi,
	double		*y,
	double		*yp)
{
	double	T, rho;
	double	p[3], e[3];

	rho = ses_rt_rho_from_grid(*xi,lseos);
	T = lseos->de_params.temperature;
	s2eos_lookup(rho,T,p,e,lseos);
	yp[0] = (dses_rt_rho_from_grid(rho,lseos)/rho)*(rho*e[1] - p[0]/rho)/T;
}		/*end dSdxi*/

/*
*			dSdeta():
*
*	Derivative of entropy along an isochore.  Temperature is given
*	as a function of the variable eta.
*
*	This function is in a form suitable for use by the ode solver
*	lsode.  Since lsode is written in Fortran the arguments are
*	passed as pointers.
*
*	Input:
*		*neg = number of equations (unused, assumed to be 1)
*		*eta = ses_rt_grid_from_temp(T,seos)
*		*y = entropy value, (unused)
*	Output:
*		*yp = dS/deta
*/

/*ARGSUSED*/
LOCAL	void dSdeta(
	int		*neq,
	double		*eta,
	double		*y,
	double		*yp)
{
	double	T, rho;
	double	p[3], e[3];

	rho = lseos->de_params.rho;
	T = ses_rt_temp_from_grid(*eta,lseos);
	s2eos_lookup(rho,T,p,e,lseos);
	yp[0] = (dses_rt_temp_from_grid(T,lseos)/T)*e[2];
}		/*end dSdeta*/

/*
*			integrate_entropy():
*
*	Utility function for integrating entropy along the line from the mid
*	point in
*			ses_rt_grid_from_rho(density,seos) /
*			ses_rt_grid_from_temp(temperature,seos)
*	space to the point (rho, T).
*/

/*ARGSUSED*/
LOCAL	void integrate_entropy(
	int		*neq,
	double		*alpha,
	double		*y,
	double		*yp)
{
	double	rlgmid = lseos->de_params.rlgmid;
	double	tlgmid = lseos->de_params.tlgmid;
	double	xi = lseos->de_params.xi;
	double	eta = lseos->de_params.eta;
	double	T, rho;
	double	p[3], e[3];

	rho = ses_rt_rho_from_grid(*alpha*xi + rlgmid,lseos);
	T = ses_rt_temp_from_grid(*alpha*eta + tlgmid,lseos);
	s2eos_lookup(rho,T,p,e,lseos);
	eta *= dses_rt_temp_from_grid(T,lseos)/T;
	xi *= dses_rt_rho_from_grid(rho,lseos)/rho;
	yp[0] = eta*e[2] + xi*(rho*e[1] - p[0]/rho)/T;
}		/*end integrate_entropy*/

#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

#endif /* defined(SESAME_CODE) && defined(TWOD) */
