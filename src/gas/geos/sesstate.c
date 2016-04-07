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
*				sesstate.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if defined(SESAME_CODE) && defined(TWOD) && defined(PHASE_CODE)
#include <geos/sesame.h>

	/*LOCAL function prototypes*/
#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
    LOCAL   LSODE_FUNC      dsdeta1;
    LOCAL   LSODE_FUNC      dsdt;
    LOCAL   LSODE_FUNC      dsdxi;
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */
LOCAL	void	phasest(double,double*);

LOCAL	SESAME_EOS	*lseos = NULL;
LOCAL	double		rph[2], eph[2], pph[2];
LOCAL	int		istate;

/*
*				phbnd():
*
*	Sets up table window and computes the entropy on the reference curve.
*	S(rhomax, tmin) arbitrarily set to 0.
*
*	The rectangular (rho,T) window is first mapped onto the unit square.
*	This version calculates the entropy along the t2min isotherm so that
*	it can be used as an initial condition in calcualting the entropy
*	along the isopicts. Phase boundaries are allowed. The entropy is
*	set to 0 at the point (t2min, rmin).
*
*        input:
*		seos structure which includes
*			rmin, rmax, tpmin, tpmax = define rectangular region
*						   for which new inversions
*						   are to be done.
*			tbls = original SESAME tables.
*			abser0, reler0 = absolute and relative error
*					 tolerances for the ODE solver. 
*             		nrhyp, nthyp = number of rho and T points in the
*				       hyp structure.
*	output:
*		sref = Entropy values along the reference curve.
*/

EXPORT	void	phbnd(
	double		*sref,
	double		rcmin,
	double		rcmax,
	double		sinit,
	SESAME_EOS	*seos)
{
	boolean phase = multiphase_eos(seos);
	double	rmin = Rho_min(seos);
	double	rmax = Rho_max(seos);
	double	*tbls = seos->sestab.tbls;
	double	tmax = Temp_max(seos);
	double	tmin = Temp_min(seos);
	double	rerr = seos->reler0, aerr = seos->abser0;
	double	eta;
	int	nrhyp = Nrho_hyp(seos); 
	int	i;
	int	nr, nt;
	int	neta;


	seos->de_params.dxi = 0.0;
	/*
	*  Set tmin to be an isotherm given by Sesame so splining
	*  can be done and S calcualted.
	*/
	seos->de_params.rlgmin = log(rmin);
	nr = (int)(tbls[2]);
	nt = (int)(tbls[3]);
	for (i = 0; i <= nt; i++)
	{
	    if (tbls[nr + 4 + i] <= tmin && tbls[nr + 4 + i + 1] > tmin)
	        seos->de_params.t2min = tbls[nr + 4 + i];
	}
	seos->de_params.tlgmin = log(seos->de_params.t2min);
	seos->de_params.dlogr = log(rmax) - seos->de_params.rlgmin;
	seos->de_params.dlogt = log(tmax) - log(seos->de_params.t2min);
	neta = nrhyp;
	seos->de_params.deta = (rcmax-rcmin)/((neta-1.0)*seos->de_params.dlogr);
	if (phase == YES)
	    phasest(seos->de_params.t2min,tbls);

	/*            compute the entropy on the line xi = 0 */

	eta = (rcmin - seos->de_params.rlgmin)/seos->de_params.dlogr;

	lseos = seos;
	sref[0] = sinit;
	istate = 1;
	if (!ode_solver(dsdeta1,sref,eta,seos->de_params.deta,
		        neta,rerr,aerr,0.0,&istate,NULL))
	{
	    screen("ERROR in phbnd(), ode_solver() failed\n");
	    clean_up(ERROR);
	}
}		/*end phbnd*/

EXPORT	void	setspb(
	double		sinit,
	double		*s,
	int		npts,
	double		tph,
	SESAME_EOS	*seos)
{
	double	rerr = seos->reler0, aerr = seos->abser0;
	double	xi;
	int	nxi;

	nxi = npts + 1;
	seos->de_params.dxi = (tph - seos->de_params.tlgmin) / npts;
	lseos = seos;
	s[0] = sinit;
	xi = log(seos->de_params.t2min);
	istate = 1;
	if (!ode_solver(dsdt,s,xi,seos->de_params.dxi,nxi,
		        rerr,aerr,0.0,&istate,NULL))
	{
	    screen("ERROR in setspb(), ode_solver() failed\n");
	    clean_up(ERROR);
	}
}		/*end setspb*/

/*
*				sets():
*
*	Calculate the entropy along the isochores eta = const.
*/

EXPORT	void	sets(
	double		sinit,
	double		*s, 
	int		nht,
	int		*iplace,
	SESAME_EOS	*seos)
{
	double	*tt = seos->de_params.T;
	double	*tbls = seos->sestab.tbls;
	double	tmax = Temp_max(seos);
	double	tmin = Temp_min(seos);
	double	rerr = seos->reler0, aerr = seos->abser0;
	double	xi;
	int	i;
	int	nxi;

	xi = 0.0;
	if (*iplace != -1)
	    seos->de_params.t2min = tt[(*iplace)++];
	else
	{
	    int nr = (int)(tbls[2]), nt = (int)(tbls[3]);
	    for (i = 0; i < nt; i++)
	    {
	        if (tbls[nr + 4 + i] <= tmin && tbls[nr + 4 + i + 1] > tmin)
		    seos->de_params.t2min = tbls[nr + 4 + i];
	    }
	}
	seos->de_params.tlgmin = log(seos->de_params.t2min);
	seos->de_params.dlogt = log(tmax) - seos->de_params.tlgmin;
	nxi = nht + 1;
	seos->de_params.dxi = 1.0 / nht;

	/* set tbls,tt,ee,pp for dsdxi */

	lseos = seos;
	s[0] = sinit;
	istate = 1;
	if (!ode_solver(dsdxi,s,xi,seos->de_params.dxi,nxi,rerr,aerr,
		        0.0,&istate,NULL))
	{
	    screen("ERROR in sets(), ode_solver() failed\n");
	    clean_up(ERROR);
	}
}		/*end sets*/

LOCAL const double MIN_TEMP_DIFF = 0.01; /*TOLERANCE*/

LOCAL	void phasest(
	double		temp,
	double		*tbls)
{
	int		nt,nr,ntv,n4;
	int		i;
	double		el,eg,rl,rg;
	double		el1,eg1,rl1,rg1;
	double		pv1, pv;

	nr = (int)(tbls[2]);
	nt = (int)(tbls[3]);
	n4 = nr+nt+4+3*nr*nt;
	ntv = (int)(tbls[n4]);		 
	
	n4 = n4+1;

	for (i = 0; i < ntv-1; i++)
	{
	    if (fabs(temp - tbls[n4+ntv+i]) < MIN_TEMP_DIFF)
	    {
	    	rph[0] = tbls[n4+2*ntv+i];
	    	rph[1] = tbls[n4+3*ntv+i];
	    	pph[0] = pph[1] = tbls[n4+i];
	    	eph[0] = tbls[n4+4*ntv+i];
	    	eph[1] = tbls[n4+5*ntv+i];
	    	return;
	    }
	    else if (temp >= tbls[n4+ntv+i] && temp < tbls[n4+ntv+i+1])
	    {
	    	rg = tbls[n4+2*ntv+i];	rl = tbls[n4+3*ntv+i];
	    	pv = tbls[n4+i];	pv1 = tbls[n4+i+1];
	    	eg = tbls[n4+4*ntv+i];	el = tbls[n4+5*ntv+i];
	    	rg1 = tbls[n4+2*ntv+i+1];  rl1 = tbls[n4+3*ntv+i+1];
	    	eg1 = tbls[n4+4*ntv+i+1];  el1 = tbls[n4+5*ntv+i+1];
	    	rg = log(rg);
	    	rg1 = log(rg1);
	    	rl = log(rl);
	    	rl1 = log(rl1);
	    	pv = log(pv);
	    	pv1 = log(pv1);
	    	rg = (rg-rg1)/(log(1.+tbls[n4+ntv+i])-log(1.+tbls[n4+ntv+i+1]));
	    	rg = rg*(log(1. + temp) - log(1.+tbls[n4+ntv+i+1])) + rg1;
	    	rl = (rl-rl1)/(log(1.+tbls[n4+ntv+i])-log(1.+tbls[n4+ntv+i+1]));
	    	rl = rl*(log(1. + temp) - log(1.+tbls[n4+ntv+i+1])) + rl1;
	    	eg = (eg-eg1)/(log(1.+tbls[n4+ntv+i])-log(1.+tbls[n4+ntv+i+1]));
	    	eg = eg*(log(1. + temp) - log(1.+tbls[n4+ntv+i+1])) + eg1;
	    	el = (el-el1)/(log(1.+tbls[n4+ntv+i])-log(1.+tbls[n4+ntv+i+1]));
	    	el = el*(log(1. + temp) - log(1.+tbls[n4+ntv+i+1])) + el1;
	    	pv = (pv-pv1)/(log(1.+tbls[n4+ntv+i])-log(1.+tbls[n4+ntv+i+1]));
	    	pv = pv*(log(1. + temp) - log(1.+tbls[n4+ntv+i+1])) + pv1;
	    	rph[0] = exp(rg);
	    	rph[1] = exp(rl);
	    	eph[0] = eg;
	    	eph[1] = el;
	    	pph[0] = pph[1] = exp(pv);
	    }
	}
}		/*end phasest*/

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */

/*ARGSUSED*/
LOCAL	void dsdeta1(
	int		*neq,
	double		*peta,
	double		*y,
	double		*yp)
{
	double	*tbls = lseos->sestab.tbls;
	double	*rr = lseos->de_params.R;
	double	*ee = lseos->de_params.E;
	double	*pp = lseos->de_params.P;
	double	*slopee = lseos->de_params.slopeE;
	double	*slopep = lseos->de_params.slopeP;
	double	eta = *peta;
	double	x1, rho;
	double	p, e, dp;
	double	rlgmin, dlogr;
	double	t2min,de, r1, r2, s1, s2, p1, p2, e1, e2;
	int	nr, i, kmax, nrp3;
	static const double TOLEPS = 1.0e-05; /*TOLERANCE*/

	rlgmin = lseos->de_params.rlgmin;
	t2min = lseos->de_params.t2min;
	dlogr = lseos->de_params.dlogr;

	rho = rlgmin + eta*dlogr;
	nr = (int)(tbls[2]);

	if(rho <= rr[0])
	{
	    rho = exp(rho);
	    if (rph[0] <= rho && rho <= rph[1])
	    {
	    	/* linear approx in mixed phase */
	    	p = pph[0];
	    	de = rph[1]*rph[0]*(eph[1] - eph[0])/
	    		(rph[0] - rph[1]);
	    	de = -1.0*de/rho;
	    }
	    else
	    {
	    	p = exp(pp[0]);
	    	de = (exp(ee[0]))*slopee[0];
	    }
	    x1 = de - p/rho;
	    yp[0] = dlogr*x1/t2min;
	    return;
	}

	if (istate == 1)
	{
	    if (fabs(rph[0] - exp(rho)) < TOLEPS*dlogr)
	    {
	        /* linear approx in mixed phase */
	        rho = exp(rho);
	        p = pph[0];
	        de = rph[1]*(eph[1] - eph[0])/(rph[0] - rph[1]);
	        de = -1.0*de;
	        x1 = de - p/rho;
	        yp[0] = dlogr*x1/t2min;
	        return;
	    }
	    nrp3 = nr + 3;
	    for(i=0; i < nrp3; i++)
	    {

	        if (fabs(rph[1] - exp(rr[i])) < TOLEPS*dlogr
	           && max(fabs(rho - rr[i]),fabs(rho - rr[i+1]))< TOLEPS*dlogr)
	        {
	    	    p = exp(pp[i+1]);
	    	    de = (exp(ee[i+1]))*slopee[i+1];
	    	    rho = exp(rho);
	    	    x1 = de - p/rho;
	    	    yp[0] = dlogr*x1/t2min;
	    	    return;
	        }
	    }
	    for(i=0; i < nrp3; i++)
	    {
	        if (rho >= rr[i] && rho < rr[i+1])
	        {
	    	    r1 = rr[i];
	    	    r2 = rr[i+1];
	    	    p1 = pp[i];
	    	    p2 = pp[i+1];
	    	    s1 = slopep[i];
	    	    s2 = slopep[i+1];
	    	    kmax = i;
	    	    spline(r1,r2,p1,p2,s1,s2,rho,&p,&dp);
	    	    s1 = slopee[kmax];
	    	    s2 = slopee[kmax+1];
	    	    e1 = ee[kmax];
	    	    e2 = ee[kmax+1];
	    	    spline(r1,r2,e1,e2,s1,s2,rho,&e,&de);
	    	    p = exp(p);
	    	    rho = exp(rho);
	    	    de = (exp(e))*de;
	    	    x1 = de - p/rho;
	    	    yp[0] = dlogr*x1/t2min;
	    	    return;
	        }
	    }
	}
	else
	{
	    if (rph[0] < exp(rho) && exp(rho) <= rph[1])
	    {
	        /* linear approx in mixed phase */
	        rho = exp(rho);
	        p = pph[0];
	        de = rph[1]*rph[0]*(eph[1] - eph[0])/(rph[0] - rph[1]);
	        de = -1.0*de/rho;
	        x1 = de - p/rho;
	        yp[0] = dlogr*x1/t2min;
	        return;
	    }
	    nrp3 = nr + 3;
	    for(i=0; i < nrp3; i++)
	    {
	        if (rho > rr[i] && rho <= rr[i+1])
	        {
	    	    r1 = rr[i];
	    	    r2 = rr[i+1];
	    	    p1 = pp[i];
	    	    p2 = pp[i+1];
	    	    s1 = slopep[i];
	    	    s2 = slopep[i+1];
	    	    kmax = i;
	    	    spline(r1,r2,p1,p2,s1,s2,rho,&p,&dp);
	    	    s1 = slopee[kmax];
	    	    s2 = slopee[kmax+1];
	    	    e1 = ee[kmax];
	    	    e2 = ee[kmax+1];
	    	    spline(r1,r2,e1,e2,s1,s2,rho,&e,&de);
	    	    p = exp(p);
	    	    rho = exp(rho);
	    	    de = (exp(e))*de;
	    	    x1 = de - p/rho;
	    	    yp[0] = dlogr*x1/t2min;
	    	    return;
	        }
	    }
	}


	screen("ERROR in dsdeta1(), unable to locate density\n");
	clean_up(ERROR);
}		/*end dsdeta1*/

/*ARGSUSED*/
LOCAL	void	dsdt(
	int		*neq,
	double		*pxi,
	double		*y,
	double		*yp)
{
	double	*tbls = lseos->sestab.tbls;
	double	*e = lseos->de_params.pb_E;
	double	*p = lseos->de_params.pb_P;
	double	*t = lseos->de_params.pb_T;
	double	*r = lseos->de_params.pb_R;
	double	*slpe = lseos->de_params.pb_slopeE;
	double	*slpp = lseos->de_params.pb_slopeP;
	double	*slpr = lseos->de_params.pb_slopeR;
	double	xi = *pxi;
	double	temp;
	double	s1, s2, e1, e2, de;
	double	p1, p2, dp, r1, r2, dr;
	double	t1, t2;
	double	R, E, P;
	int	i, nt, nr, n4, ntv;

	temp = exp(xi);
	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);
	n4 = nr+nt+3*nr*nt+5;
	ntv = (int)(tbls[n4-1]);
	if (tbls[4] <= 0.) ntv = ntv-1;

	for (i = 0; i < ntv-1; i++)
	{
	    if (fabs(t[i+1] - xi) < MIN_TEMP_DIFF)
	    {
	    	de = exp(e[i+1])*slpe[i+1];
	    	P = exp(p[i+1]);
	    	R = exp(r[i+1]);
	    	dr = slpr[i+1];
	    	yp[0] = (de - P/R*dr)/temp;
	    	return;
	    }
	    if ( t[i] <= xi && t[i+1] >= xi)
	    {
	    	t1 = t[i];
	    	t2 = t[i+1];
	    	s1 = slpe[i];
	    	s2 = slpe[i+1];
	    	e1 = e[i];
	    	e2 = e[i+1];
	    	spline(t1,t2,e1,e2,s1,s2,xi,&E,&de);
	    	s1 = slpp[i];
	    	s2 = slpp[i+1];
	    	p1 = p[i];
	    	p2 = p[i+1];
	    	spline(t1,t2,p1,p2,s1,s2,xi,&P,&dp);
	    	s1 = slpr[i];
	    	s2 = slpr[i+1];
	    	r1 = r[i];
	    	r2 = r[i+1];
	    	spline(t1,t2,r1,r2,s1,s2,xi,&R,&dr);
	    	de = exp(E)*de;
	    	yp[0] = (de - exp(P)/exp(R)*dr)/temp;
	    	return;
	    }
	}
	screen("ERROR in dsdt(), unable to locate temperature\n");
	clean_up(ERROR);
}		/*end dsdt*/

/*ARGSUSED*/
LOCAL	void dsdxi(
	int		*neq,
	double		*pxi,
	double		*y,
	double		*yp)
{
	double	*tbls = lseos->sestab.tbls;
	double	*ee = lseos->de_params.E;
	double	*tt = lseos->de_params.T;
	double	*slope = lseos->de_params.slopeE;
	double	xi = *pxi;
	double	temp;
	double	e, de;
	double	dlogt;
	double	s1, s2, e1, e2;
	double	t1, t2, tlgmin;
	int	ntp1;
	int	i, nt;

	tlgmin = lseos->de_params.tlgmin;
	dlogt = lseos->de_params.dlogt;

	temp = exp(tlgmin + xi*dlogt);
	nt = (int)(tbls[3]);
	if (fabs(temp - tt[0]) < MIN_TEMP_DIFF)
	{
	    de = slope[0];
	    yp[0]= dlogt*(exp(ee[0])/(tt[0]+1.0))*de;
	    return;
	}	 
	ntp1 = nt + 1;


	for(i=0; i < ntp1; i++)
	{
	    if(fabs(temp - tt[i+1]) <= MIN_TEMP_DIFF)
	    {
	    	de = (exp(ee[i+1])/(tt[i+1]+1.0))*slope[i+1];
	    	yp[0] = dlogt*de;
	    	return;
	    }
	}
	for(i=0; i < ntp1; i++)
	{
	    if(temp > tt[i] && temp <= tt[i+1])
	    {    
	    	t1 = log(tt[i]+1.0);
	    	t2 = log(tt[i+1]+1.0);
	    	temp = log(temp+1.0);
	    	e1 = ee[i];
	    	e2 = ee[i+1];
	    	s1 = slope[i];
	    	s2 = slope[i+1];
	    	spline(t1,t2,e1,e2,s1,s2,temp,&e,&de);
	    	yp[0]= dlogt*(exp(e)/exp(temp))*de;
	    	return;
	    }
	}
	screen("ERROR in dsdxi(), ");
	screen("unable to locate temperature\n");
	clean_up(ERROR);
}		/*end dsdxi*/

#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

#endif /* defined(SESAME_CODE) && defined(TWOD) && defined(PHASE_CODE) */
