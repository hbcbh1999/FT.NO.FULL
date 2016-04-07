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
*                		grsolve.c
*
*    Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    Implementation of various approximate Riemann solvers to be used
*    in finite difference methods.
*
*/


#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

/* LOCAL Function Prototypes */
LOCAL	void	load_derived_fields(int,int,Vec_Gas*,Vec_Gas*);

/*ARGSUSED*/
EXPORT	void	g_muscl_exact_rsolver(
	int	   start,
	int	   end,
	double      **uL,
	Vec_Gas    *vlst,
	double      **uR,
	Vec_Gas    *vrst,
	double      **uM,
	Vec_Gas    *vmst,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	static Locstate sl = NULL, sr = NULL, smid = NULL;
	static double	dir[3] = { 1.0, 0.0, 0.0};
	Locstate	*state;
	double           *rhol, *rhor;
	double           *el, *er;
	double           *vl[3], *vr[3];
	int		dim;
	int		i,j;
	double           alpha = 0.0;

	debug_print("rsolve","Entered g_muscl_exact_rsolver()\n");

	dim = vmuscl->dim;
	state = vmuscl->vst->state + vmuscl->offset;

	if (sl == NULL) 
	{ 
	    size_t	sizest = vmuscl->sizest;

	    g_alloc_state(&sl,sizest);
	    g_alloc_state(&sr,sizest);
	    g_alloc_state(&smid,sizest);
	}
	rhol = uL[vmuscl->index.density];
	el = uL[vmuscl->index.energy];
	rhor = uR[vmuscl->index.density];
	er = uR[vmuscl->index.energy];
	for (i = 0; i < dim; ++i)
	{
	    vl[i] = uL[vmuscl->index.v[i]];
	    vr[i] = uR[vmuscl->index.v[i]];
	}
	for (j = start; j < end; ++j) 
	{
	    if (is_obstacle_state(state[j]))
	        continue;
	    Dens(sl) = (1-alpha)*rhol[j] + alpha*rhor[j];
	    Dens(sr) = (1-alpha)*rhor[j] + alpha*rhol[j];
	    Energy(sl) = el[j];
	    Energy(sr) = er[j];
	    for (i = 0; i < dim; ++i)
	    {
	    	Vel(sl)[i] = vl[i][j];
	    	Vel(sr)[i] = vr[i][j];
	    }
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                 /* Compute partial density for sl, sr */
                 Gas_param *params = Params(state[j]);
                 int    num_comps;
                 if((num_comps = params->n_comps) != 1)
                 {
                     for(i = 0; i < num_comps; i++)
                     {
                         pdens(sl)[i] = Dens(sl) * uL[vmuscl->index.prho[i]][j];                         pdens(sr)[i] = Dens(sr) * uR[vmuscl->index.prho[i]][j];
                     }
                 }
            }
	    Set_params(sl,state[j]);
	    set_type_of_state(sl,EGAS_STATE);
	    Set_params(sr,state[j]);
	    set_type_of_state(sr,EGAS_STATE);
	    reset_gamma(sl);
	    reset_gamma(sr);

	    riemann_solution(0.0,dir,sl,sr,smid,EGAS_STATE);

	    vmst->rho[j] = Dens(smid);
	    vmst->e[j] = Energy(smid);
	    for (i = 0; i < dim; ++i)
	    	vmst->v[i][j] = Vel(smid)[i];
	}
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Upwinding selection of the reconstructed mass fraction at zone interface */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
            {
                for (j = start; j < end; ++j)
                {
                    if (is_obstacle_state(state[j]))
                        continue;
                    if(vmst->v[0][j] > 0.0)
                    {
                        for(i = 0; i < num_comps; i++)
                            vmst->rho0[i][j] = uL[vmuscl->index.prho[i]][j];
                    }
                    else
                    {
                        for(i = 0; i < num_comps; i++)
                            vmst->rho0[i][j] = uR[vmuscl->index.prho[i]][j];
                    }
                }
                Vec_Gas_field_set(vmst,rho0) = YES;
            }
        }
	Vec_Gas_field_set(vmst,rho) = YES;
	Vec_Gas_field_set(vmst,e) = YES;
	Vec_Gas_field_set(vmst,v) = YES;
	if (Flux != NULL)
	    g_load_muscl_flux(start,end,uM,vmst,Flux,vmuscl);

	debug_print("rsolve","Left g_muscl_exact_rsolver()\n");
}		/*end g_muscl_exact_rsolver*/


EXPORT	void	g_exact_Riemann_midstate(
	int       start,
	int       end,
	Vec_Gas   *vlst,
	Vec_Gas   *vrst,
	double     *pm,
	double     *um,
	Vec_Muscl *vmuscl)
{
	Locstate	         *state;
	double                    pl, pr, ul, ur;
	double                    ml, mr;
	double                    *rhol, *rhor;
	double                    *el, *er;
	double                    *vl, *vr;
	int		         j;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	static Locstate          sl = NULL, sr = NULL;

	state = vmuscl->vst->state + vmuscl->offset;

	if (sl == NULL) 
	{ 
	    size_t	sizest = vmuscl->sizest;
	    if (sizest < sizeof(VGas))
	        sizest = sizeof(VGas);

	    g_alloc_state(&sl,sizest);
	    g_alloc_state(&sr,sizest);
	}
	rhol = vlst->rho;
	el = vlst->e;
	vl = vlst->v[0];
	rhor = vrst->rho;
	er = vrst->e;
	vr = vrst->v[0];
	for (j = start; j < end; ++j) 
	{
	    if (is_obstacle_state(state[j]))
	        continue;
	    Dens(sl) = rhol[j];
	    Dens(sr) = rhor[j];
	    Energy(sl) = el[j];
	    Energy(sr) = er[j];
	    Vel(sl)[0] = vl[j];
	    Vel(sr)[0] = vr[j];
            reset_gamma(sl);
            reset_gamma(sr);
	    Set_params(sl,state[j]);
	    set_type_of_state(sl,EGAS_STATE);
	    Set_params(sr,state[j]);
	    set_type_of_state(sr,EGAS_STATE);

	    set_state_for_find_mid_state(sl,sl);
	    set_state_for_find_mid_state(sr,sr);
	    find_mid_state(sl,sr,0.0,&pl,&pr,&ul,&ur,&ml,&mr,&l_wave,&r_wave);

	    pm[j] = 0.5*(pl+pr);
	    um[j] = 0.5*(ul+ur);
	}
}		/*end g_exact_Riemann_midstate*/

/*ARGSUSED*/
EXPORT	void g_gamma_law_fit_rsolver(
	int	   start,
	int	   end,
	double      **uL,
	Vec_Gas    *vlst,
	double      **uR,
	Vec_Gas    *vrst,
	double      **uM,
	Vec_Gas    *vmst,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	double     *rho_l = vlst->rho, *rho_r = vrst->rho, *rho_m = vmst->rho;
	double     *p_l = vlst->p, *p_r = vrst->p, *p_m = vmst->p;
	double     **v_l = vlst->v, **v_r = vrst->v, **v_m = vmst->v;
	double     *c_l = vlst->c, *c_r = vrst->c;
	double     *e_l = vlst->e, *e_r = vrst->e, *e_m = vmst->e;
	double     *FD_l = vlst->FD, *FD_r = vrst->FD;
	double     p0_l, p0_r, rho0_l, rho0_r, v0_l, v0_r;
	double     pi0_l, pi0_r, pinf_l, pinf_r;
	double     c0_l, c0_r, G0_l, G0_r;
	double     gam_l, gam_r, A0_l, A0_r, B0_l, B0_r, C0_l, C0_r;
	double     dp_l, dp_r;
	double     vm, pm, Sl, Sr;
	const double eps = 10.0*MACH_EPS;/*TOLERANCE*/
	boolean      left_shock, right_shock;
	int       dim = vmuscl->dim;
	int       j, k;

	load_derived_fields(start,end,vlst,vrst);

	g_gamma_law_fit_Riemann_midstate(start,end,vlst,vrst,p_m,v_m[0],vmuscl);

	for (j = start; j < end; ++j)
	{
	    rho0_l = rho_l[j];
	    c0_l = c_l[j];
	    G0_l = FD_l[j];
	    gam_l = 2.0*G0_l - 1.0;
	    A0_l = G0_l/gam_l;
	    B0_l = 1.0 - A0_l;
	    C0_l = 2.0*c0_l/(gam_l - 1.0);
	    p0_l = p_l[j];
	    pi0_l = rho0_l*c0_l*c0_l/gam_l;
	    pinf_l = pi0_l - p0_l;
	    v0_l = v_l[0][j];

	    rho0_r = rho_r[j];
	    c0_r = c_r[j];
	    G0_r = FD_r[j];
	    gam_r = 2.0*G0_r - 1.0;
	    A0_r = G0_r/gam_r;
	    B0_r = 1.0 - A0_r;
	    C0_r = 2.0*c0_r/(gam_r - 1.0);
	    p0_r = p_r[j];
	    pi0_r = rho0_r*c0_r*c0_r/gam_r;
	    pinf_r = pi0_r - p0_r;
	    v0_r = v_r[0][j];

	    pm = p_m[j];
	    vm = v_m[0][j];
	    dp_l = pm - p0_l;
	    left_shock = (dp_l/pi0_l >= -eps) ? YES : NO;
	    if (left_shock)
	        Sl = v0_l - c0_l*sqrt(1.0 + A0_l*dp_l/pi0_l);
	    else
	        Sl = v0_l - c0_l;
	    dp_r = pm - p0_r;
	    right_shock = (dp_r/pi0_r >= -eps) ? YES : NO;
	    if (right_shock)
	        Sr = v0_r + c0_r*sqrt(1.0 + A0_r*dp_r/pi0_r);
	    else
	        Sr = v0_r + c0_r;

	    if (Sl >= 0.0) /* Flux = Flux(left state) */
	    {
	        rho_m[j] = rho0_l;
	        e_m[j] = e_l[j];
	        p_m[j] = p0_l;
		for (k = 0; k < dim; ++k)
		    v_m[k][j] = v_l[k][j];
	    }
	    else if (Sr <= 0.0) /* Flux = Flux(right state) */
	    {
	        rho_m[j] = rho0_r;
	        e_m[j] = e_r[j];
	        p_m[j] = p0_r;
		for (k = 0; k < dim; ++k)
		    v_m[k][j] = v_r[k][j];
	    }
	    else if (0.0 <= vm) /* Flux from left wave */
	    {
	        if (left_shock) /* behind left shock */
		{
		    double mu2 = (gam_l-1.0)/(gam_l+1.0);
		    double pi_l = pm + pinf_l;
		    rho_m[j] = rho0_l*(pi_l + mu2*pi0_l)/(pi0_l + mu2*pi_l);
		    e_m[j] = e_l[j]+0.5*(pm*pm-p0_l*p0_l)*(1.0-mu2)/
		             (pi0_l + mu2*pi_l);
		    p_m[j] = pm;
	            v_m[0][j] = vm;
		}
		else /* rarefaction */
		{
		    double x, rhom_l, cm_l, vm_l;
		    x = 1.0 + dp_l/pi0_l;
		    if (x > 0.0)
		    {
		        rhom_l = rho0_l*pow(x,1.0/gam_l);
		        cm_l = sqrt(gam_l*(pm+pinf_l)/rhom_l);
		        vm_l = v0_l + C0_l*(1.0 - cm_l/c0_l);
		    }
		    else
		    {
		        x = 0.0;
		        rhom_l = Vacuum_dens(vlst->state[j]);
			cm_l = 0.0;
			vm_l = v0_l + C0_l;
		    }
		    if ((vm_l - cm_l) <= 0.0) /* behind left rarefaction */
		    {
		        rho_m[j] = rhom_l;
			e_m[j] = e_l[j] + p0_l/rho0_l - pm/rhom_l +
			        (c0_l*c0_l/(gam_l-1.0))*(pow(x,2.0*B0_l)-1.0);
		        p_m[j] = pm;
	                v_m[0][j] = vm;
		    }
		    else /* inside left rarefaction */
		    {
		        double mu2 = (gam_l-1.0)/(gam_l+1.0);
		        cm_l = mu2*(C0_l + v0_l);
			rho_m[j] = rho0_l*pow(cm_l/c0_l,2.0/(gam_l-1.0));
			p_m[j] = rho_m[j]*cm_l*cm_l/gam_l - pinf_l;
			x = 1.0 + (p_m[j]-p0_l)/pi0_l; x = max(x,0.0);
			e_m[j] = e_l[j] + p0_l/rho0_l - p_m[j]/rho_m[j] +
			        (c0_l*c0_l/(gam_l-1.0))*(pow(x,2.0*B0_l) - 1.0);
			v_m[0][j] = cm_l;
		    }
		}
		for (k = 1; k < dim; ++k)
		    v_m[k][j] = v_l[k][j];
	    }
	    else /* Flux from right wave */
	    {
	        if (right_shock) /* behind right shock */
	        {
		    double mu2 = (gam_r-1.0)/(gam_r+1.0);
		    double pi_r = pm + pinf_r;
		    rho_m[j] = rho0_r*(pi_r + mu2*pi0_r)/(pi0_r + mu2*pi_r);
		    e_m[j] = e_r[j]+0.5*(pm*pm-p0_r*p0_r)*(1.0-mu2)/
		             (pi0_r + mu2*pi_r);
		    p_m[j] = pm;
	            v_m[0][j] = vm;
	        }
	        else /* rarefaction */
	        {
		    double x, rhom_r, cm_r, vm_r;
		    x = 1.0 + dp_r/pi0_r;
		    if (x > 0.0)
		    {
		        rhom_r = rho0_r*pow(x,1.0/gam_r);
		        cm_r = sqrt(gam_r*(pm+pinf_r)/rhom_r);
		        vm_r = v0_r - C0_r*(1.0 - cm_r/c0_r);
		    }
		    else
		    {
		        x = 0.0;
		        rhom_r = Vacuum_dens(vrst->state[j]);
			cm_r = 0.0;
			vm_r = v0_r - C0_r;
		    }
		    if ((vm_r + cm_r) >= 0.0) /* behind right rarefaction */
		    {
		        rho_m[j] = rhom_r;
			e_m[j] = e_r[j] + p0_r/rho0_r - pm/rhom_r +
			    (c0_r*c0_r/(gam_r-1.0))*(pow(x,2.0*B0_r) - 1.0);
		        p_m[j] = pm;
	                v_m[0][j] = vm;
		    }
		    else /* inside right rarefaction */
		    {
		        double mu2 = (gam_r-1.0)/(gam_r+1.0);
		        cm_r = mu2*(C0_r - v0_r);
	 	        rho_m[j] = rho0_r*pow(cm_r/c0_r,2.0/(gam_r-1.0));
		        p_m[j] = rho_m[j]*cm_r*cm_r/gam_r - pinf_r;
			x = 1.0 + (p_m[j]-p0_r)/pi0_r; x = max(x,0.0);
			e_m[j] = e_r[j] + p0_r/rho0_r - p_m[j]/rho_m[j] +
			    (c0_r*c0_r/(gam_r-1.0))*(pow(x,2.0*B0_r) - 1.0);
		        v_m[0][j] = -cm_r;
		    }
	        }
	        for (k = 1; k < dim; ++k)
		    v_m[k][j] = v_r[k][j];
	    }
	}
	Vec_Gas_field_set(vmst,p) = YES;
	Vec_Gas_field_set(vmst,rho) = YES;
	Vec_Gas_field_set(vmst,e) = YES;
	Vec_Gas_field_set(vmst,v) = YES;
	g_load_muscl_flux(start,end,uM,vmst,Flux,vmuscl);
}		/*end g_gamma_law_fit_rsolver*/

/*ARGSUSED*/
EXPORT	void	g_gamma_law_fit_Riemann_midstate(
	int       start,
	int       end,
	Vec_Gas   *vlst,
	Vec_Gas   *vrst,
	double     *pm,
	double     *vm,
	Vec_Muscl *vmuscl)
{
	const double eps_u = 0.000001; /* TOLERANCE */
	const double eps_p = 0.000001; /* TOLERANCE */
	const double eps = 10.0*MACH_EPS;/*TOLERANCE*/
	double     *rho_l = vlst->rho, *rho_r = vrst->rho;
	double     *p_l = vlst->p, *p_r = vrst->p;
	double     *v_l = vlst->v[0], *v_r = vrst->v[0];
	double     *c_l = vlst->c, *c_r = vrst->c;
	double     *FD_l = vlst->FD, *FD_r = vrst->FD;
	double     m_l, m_r, i0_l, i0_r, c0_l, c0_r;
	double     gam_l, gam_r, dp_l, dp_r;
	double     p0_l, p0_r, v0_l, v0_r;
	double     pi0_l, pi0_r;
	double     A0_l, A0_r, B0_l, B0_r, C0_l, C0_r,  D_l, D_r;
	double     dv, vml, vmr;
	double     pnew;
	double     min_p, min_pl, min_pr;
	double     alpha;
	boolean      godunov, left_vac, right_vac;
	int       j, k;
	const int N = 4, M = 4;

	load_derived_fields(start,end,vlst,vrst);

	for (j = start; j < end; ++j)
	{
	    c0_l = c_l[j];
	    i0_l = rho_l[j]*c0_l;
	    gam_l = 2.0*FD_l[j] - 1.0;
	    pi0_l = rho_l[j]*c0_l*c0_l/gam_l;
	    A0_l = FD_l[j]/gam_l;
	    B0_l = 1.0 - A0_l;
	    C0_l = 2.0*c0_l/(gam_l - 1.0);
	    p0_l = p_l[j];
	    v0_l = v_l[j];
	    min_pl = Min_pressure(vlst->state[j]);

	    c0_r = c_r[j];
	    i0_r = rho_r[j]*c0_r;
	    gam_r = 2.0*FD_r[j] - 1.0;
	    pi0_r = rho_r[j]*c0_r*c0_r/gam_r;
	    A0_r = FD_r[j]/gam_r;
	    B0_r = 1.0 - A0_r;
	    C0_r = 2.0*c0_r/(gam_r - 1.0);
	    p0_r = p_r[j];
	    v0_r = v_r[j];
	    min_pr = Min_pressure(vrst->state[j]);

	    min_p = max(min_pl,min_pr);

	    dv = v0_l - v0_r;
	    pm[j] = 0.5*(p0_l + p0_r);
	    dp_l = pm[j] - p0_l;
	    if (dp_l/pi0_l >= -eps)
	        m_l = i0_l*sqrt(1.0 + A0_l*dp_l/pi0_l);
	    else
	    {
	        double x = 1.0 + dp_l/pi0_l; x = max(x,0.0);
	        m_l = i0_l*B0_l*(1.0 - x)/(1.0 - pow(x,B0_l));
	    }
	    dp_r = pm[j] - p0_r;
	    if (dp_r/pi0_r >= -eps)
	        m_r = i0_r*sqrt(1.0 + A0_r*dp_r/pi0_r);
	    else
	    {
	        double x = 1.0 + dp_r/pi0_r; x = max(x,0.0);
	        m_r = i0_r*B0_r*(1.0 - x)/(1.0 - pow(x,B0_r));
	    }
	    pm[j] = (m_l*m_r*dv + m_r*p0_l + m_l*p0_r)/(m_l+m_r);
	    if (pm[j] < min_p)
	        pm[j] = min_p;

	    alpha = 1.0;
	    for (k = 0; k < N*M; ++k)
	    {
	        godunov = NO;
	        left_vac = right_vac = NO;
	        dp_l = pm[j] - p0_l;
	        dp_r = pm[j] - p0_r;
		if (dp_l/pi0_l >= -eps)
		{
	            m_l = i0_l*sqrt(1.0 + A0_l*dp_l/pi0_l);
		    D_l = (1.0 - 0.5*A0_l*dp_l/(pi0_l + A0_l*dp_l))/m_l;
	            vml = v0_l - dp_l/m_l;
		}
		else
		{
		    double x = 1.0 + dp_l/pi0_l;
		    if (x > 0.0)
		    {
		        double y = pow(x,-A0_l);
		        m_l = i0_l*B0_l*(1.0 - x)/(1.0 - x*y);
		        D_l = y/i0_l;
	                vml = v0_l - C0_l*(x*y - 1.0);
		    }
		    else
		    {
		        left_vac = YES;
		        m_l = i0_l*B0_l;
	                vml = v0_l + C0_l;
			godunov = YES;
		    }
		}
		if (dp_r/pi0_r >= -eps)
		{
	            m_r = i0_r*sqrt(1.0 + A0_r*dp_r/pi0_r);
		    D_r = (1.0 - 0.5*A0_r*dp_r/(pi0_r + A0_r*dp_r))/m_r;
		    vmr = v0_r + dp_r/m_r;
		}
		else
		{
		    double x = 1.0 + dp_r/pi0_r;
		    if (x > 0.0)
		    {
		        double y = pow(x,-A0_r);
		        m_r = i0_r*B0_r*(1.0 - x)/(1.0 - x*y);
		        D_r = y/i0_r;
	                vmr = v0_r + C0_r*(x*y - 1.0);
		    }
		    else
		    {
		        right_vac = YES;
		        m_r = i0_r*B0_r;
	                vmr = v0_r - C0_r;
			godunov = YES;
		    }
		}
		if (left_vac && right_vac && (vml < vmr))
		{
		    pm[j] = min_p;
		    break;
		}
		if (!godunov)
		{
	            pnew = pm[j] + (dv - dp_l/m_l - dp_r/m_r)/(D_l+D_r);
		    if (pnew < min_p)
		        godunov = YES;
		}
		if (godunov || (k > 0 && k%M == 0))
		{
		    pnew = (m_l*m_r*dv + m_r*p0_l + m_l*p0_r)/(m_l+m_r);
		    if ((pnew < 0.0) || (k > 0 && k%M == 0))
		    {
		        alpha *= 0.5;
			pnew = alpha*pnew + (1.0 - alpha)*pm[j];
		    }
		}
	        if (pnew < min_p)
	            pnew = min_p;
		if ((fabs(pnew-pm[j])<pnew*eps_p || fabs(pnew-pm[j]) < eps)
		    && 
		    (fabs(vml - vmr) < 0.5*(fabs(vml)+fabs(vmr))*eps_u) ||
		    (fabs(vml - vmr) < eps))
		{
		    double pnnew = (m_l*m_r*dv + m_r*p0_l + m_l*p0_r)/(m_l+m_r);
		    if (fabs(pnnew-pm[j])<pnnew*eps_p || fabs(pnnew-pm[j])<eps)
		    {
	                pm[j] = pnew;
		        break;
		    }
		    break;
		}
	        pm[j] = pnew;
	    }
	    left_vac = right_vac = NO;
	    dp_l = pm[j] - p0_l;
	    if (dp_l/pi0_l >= -eps)
	    {
	        m_l = i0_l*sqrt(1.0 + A0_l*dp_l/pi0_l);
	        vml = v0_l - dp_l/m_l;
	    }
	    else
	    {
	        double x = 1.0 + dp_l/pi0_l;
	        if (x > 0.0)
		{
		    double y = pow(x,B0_l);
	            vml = v0_l - C0_l*(y - 1.0);
		}
	        else
	        {
	            left_vac = YES;
	            vml = v0_l + C0_l;
	        }
	    }
	    dp_r = pm[j] - p0_r;
	    if (dp_r/pi0_r >= -eps)
	    {
	        m_r = i0_r*sqrt(1.0 + A0_r*dp_r/pi0_r);
	        vmr = v0_r + dp_r/m_r;
	    }
	    else
	    {
	        double x = 1.0 + dp_r/pi0_r;
		if (x > 0.0)
		{
		    double y = pow(x,B0_r);
	            vmr = v0_r + C0_r*(y - 1.0);
		}
		else
		{
		    right_vac = YES;
	            vmr = v0_r - C0_r;
		}
	    }
	    vm[j] = 0.5*(vml + vmr);
	    if (((fabs(vml - vmr) > 0.5*(fabs(vml)+fabs(vmr))*eps_u) &&
		    (fabs(vml - vmr) > eps)) || left_vac || right_vac)
	        g_exact_Riemann_midstate(j,j+1,vlst,vrst,pm,vm,vmuscl);
	}
}		/*end g_gamma_law_fit_Riemann_midstate*/

/*ARGSUSED*/
EXPORT	void g_linear_us_up_rsolver(
	int	   start,
	int	   end,
	double      **uL,
	Vec_Gas    *vlst,
	double      **uR,
	Vec_Gas    *vrst,
	double      **uM,
	Vec_Gas    *vmst,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	double     *rho_l = vlst->rho, *rho_r = vrst->rho, *rho_m = vmst->rho;
	double     *p_l = vlst->p, *p_r = vrst->p, *p_m = vmst->p;
	double     **v_l = vlst->v, **v_r = vrst->v, **v_m = vmst->v;
	double     *c_l = vlst->c, *c_r = vrst->c;
	double     *e_l = vlst->e, *e_r = vrst->e, *e_m = vmst->e;
	double     *FD_l = vlst->FD, *FD_r = vrst->FD;
	double     p0_l, p0_r, rho0_l, rho0_r, v0_l, v0_r, i0_l, i0_r;
	double     c0_l, c0_r, e0_l, e0_r, rG0_l, G0_l, G0_r, rG0_r;
	double     dv_l, dv_r;
	double     vm, pm, Sl, Sr, m_l, m_r;
	double     pmin_l, pmax_l, pmin_r, pmax_r;
	double     vmin_r, vmax_r, vmin_l, vmax_l;
	int       dim = vmuscl->dim;
	int       j, k;
	const int N = 4;

	load_derived_fields(start,end,vlst,vrst);

	for (j = start; j < end; ++j)
	{
	    v0_l = v_l[0][j];
	    rho0_l = rho_l[j];
	    p0_l = p_l[j];
	    c0_l = c_l[j];
	    i0_l = rho_l[j]*c0_l;
	    G0_l = FD_l[j];
	    rG0_l = rho0_l*G0_l;

	    v0_r = v_r[0][j];
	    rho0_r = rho_r[j];
	    p0_r = p_r[j];
	    c0_r = c_r[j];
	    i0_r = rho0_r*c0_r;
	    G0_r = FD_r[j];
	    rG0_r = rho0_r*G0_r;

	    pmin_l = p_l[j] - 0.5*rho0_l*c0_l*c0_l/G0_l;
	    if (pmin_l < 0.0)
	    {
	        vmax_l = v0_l + (c0_l-sqrt(c0_l*c0_l-2.0*G0_l*p0_l/rho0_l));
	        pmin_l = 0.0;
	    }
	    else
	        vmax_l = v0_l + c0_l/G0_l;
	    if (G0_l > 2.0)
	    {
	        pmax_l =  HUGE_VAL;
		vmin_l = -HUGE_VAL;
	    }
	    else
	    {
		vmin_l = v0_l - c0_l/(1.0 - 0.5*G0_l);
	        pmax_l = p0_l + rho0_l*sqr(vmin_l - v0_l);
	    }
	    pmin_r = p_r[j] - 0.5*rho0_r*c0_r*c0_r/G0_r;
	    if (pmin_r < 0.0)
	    {
	        vmin_r = v0_r - (c0_r-sqrt(c0_r*c0_r-2.0*G0_r*p0_r/rho0_r));
	        pmin_r = 0.0;
	    }
	    else
	        vmin_r = v0_r - c0_r/G0_r;
	    if (G0_r > 2.0)
	    {
	        pmax_r = HUGE_VAL;
		vmax_r = HUGE_VAL;
	    }
	    else
	    {
		vmax_r = v0_r + c0_r/(1.0 - 0.5*G0_r);
	        pmax_r = p0_r + rho0_r*sqr(vmax_r - v0_r);
	    }
	    vm = 0.5*(v0_l+v0_r);
	    for (k = 0; k < N; ++k)
	    {
	        dv_l = (vm-v0_l);
	        dv_r = (vm-v0_r);
	        vm -= ( (p0_r+dv_r*(0.5*rG0_r*dv_r+i0_r)) -
		           (p0_l+dv_l*(0.5*rG0_l*dv_l-i0_l)) )/
			 ( (rG0_r*dv_r+i0_r) - (rG0_l*dv_l-i0_l) );
	    }
	    v_m[0][j] = vm;
	    dv_l = (vm-v0_l);
	    dv_r = (vm-v0_r);
	    pm = p_m[j] = 0.5*( (p0_r+dv_r*(0.5*rG0_r*dv_r+i0_r)) +
		              (p0_l+dv_l*(0.5*rG0_l*dv_l-i0_l)) );
	    if ((pm < pmin_l) || (pm < pmin_r) ||
	        (pm > pmax_l) || (pm > pmax_r) ||
	        (vm < vmin_l) || (vm < vmin_r) ||
	        (vm > vmax_l) || (vm > vmax_r))
	    {
		g_muscl_exact_rsolver(j,j+1,uL,vlst,uR,vrst,uM,vmst,
	                              NULL,vmuscl);
	        Vec_Gas_field_set(vmst,p) = NO;
		load_pressure(vmst,j,1);
	        Vec_Gas_field_set(vmst,p) = NO;
	    }
	    else
	    {
	        m_l = i0_l + 0.5*rho0_l*dv_l;
	        m_r = i0_r + 0.5*rho0_r*dv_r;
	        Sl = v0_l - ((m_l > rho0_l*c0_l) ? m_l/rho0_l : c0_l);
	        Sr = v0_r + ((m_r > rho0_r*c0_r) ? m_r/rho0_r : c0_r);
	        if (Sl > 0.0)
	        {
	            rho_m[j] = rho0_l;
	            e_m[j] = e_l[j];
	            p_m[j] = p0_l;
		    for (k = 1; k < dim; ++k)
		        v_m[k][j] = v_l[k][j];
	        }
	        else if (0.0 < vm)
	        {
	            e0_l = e_l[j];
		    dv_l /= c0_l;
		    rho_m[j] = rho0_l*(1.0 - 0.5*G0_l*dv_l)/
		                      (1.0 - (0.5*G0_l-1.0)*dv_l);
		    e_m[j] = e0_l + 0.5*(p0_l + pm)*(1.0/rho0_l - 1.0/rho_m[j]);
		    for (k = 1; k < dim; ++k)
		        v_m[k][j] = v_l[k][j];
	        }
	        else if (0.0 < Sr)
	        {
	            e0_r = e_r[j];
		    dv_r /= c0_r;
		    rho_m[j] = rho0_r*(1.0 + 0.5*G0_r*dv_r)/
		                      (1.0 + (0.5*G0_r-1.0)*dv_r);
		    e_m[j] = e0_r + 0.5*(p0_r + pm)*(1.0/rho0_r - 1.0/rho_m[j]);
		    for (k = 1; k < dim; ++k)
		        v_m[k][j] = v_r[k][j];
	        }
	        else
	        {
	            rho_m[j] = rho0_r;
	            e_m[j] = e_r[j];
	            p_m[j] = p0_r;
		    for (k = 1; k < dim; ++k)
		        v_m[k][j] = v_r[k][j];
	        }
	    }
	}
	Vec_Gas_field_set(vmst,p) = YES;
	Vec_Gas_field_set(vmst,rho) = YES;
	Vec_Gas_field_set(vmst,e) = YES;
	Vec_Gas_field_set(vmst,v) = YES;
	g_load_muscl_flux(start,end,uM,vmst,Flux,vmuscl);
}		/*end g_linear_us_up_rsolver*/

EXPORT	void	g_linear_us_up_Riemann_midstate(
	int       start,
	int       end,
	Vec_Gas   *vlst,
	Vec_Gas   *vrst,
	double     *pm,
	double     *vm,
	Vec_Muscl *vmuscl)
{
	double     *rho_l = vlst->rho, *rho_r = vrst->rho;
	double     *p_l = vlst->p, *p_r = vrst->p;
	double     *v_l = vlst->v[0], *v_r = vrst->v[0];
	double     *c_l = vlst->c, *c_r = vrst->c;
	double     *FD_l = vlst->FD, *FD_r = vrst->FD;
	double     p0_l, p0_r, rho0_l, rho0_r, v0_l, v0_r, i0_l, i0_r;
	double     c0_l, c0_r, rG0_l, G0_l, G0_r, rG0_r;
	double     pmin_l, pmax_l, pmin_r, pmax_r;
	double     vmin_r, vmax_r, vmin_l, vmax_l;
	double     dv_l, dv_r;
	int       j, k;
	const int N = 4;

	load_derived_fields(start,end,vlst,vrst);

	for (j = start; j < end; ++j)
	{
	    v0_l = v_l[j];
	    rho0_l = rho_l[j];
	    p0_l = p_l[j];
	    c0_l = c_l[j];
	    i0_l = rho_l[j]*c0_l;
	    G0_l = FD_l[j];
	    rG0_l = rho0_l*G0_l;

	    v0_r = v_r[j];
	    rho0_r = rho_r[j];
	    p0_r = p_r[j];
	    c0_r = c_r[j];
	    i0_r = rho0_r*c0_r;
	    G0_r = FD_r[j];
	    rG0_r = rho0_r*G0_r;

	    pmin_l = p_l[j] - 0.5*rho0_l*c0_l*c0_l/G0_l;
	    if (pmin_l < 0.0)
	    {
	        vmax_l = v0_l + (c0_l-sqrt(c0_l*c0_l-2.0*G0_l*p0_l/rho0_l));
	        pmin_l = 0.0;
	    }
	    else
	        vmax_l = v_l[j] + c0_l/G0_l;
	    if (G0_l > 2.0)
	    {
	        pmax_l =  HUGE_VAL;
		vmin_l = -HUGE_VAL;
	    }
	    else
	    {
		vmin_l = v0_l - c0_l/(1.0 - 0.5*G0_l);
	        pmax_l = p0_l + rho0_l*sqr(vmin_l - v0_l);
	    }
	    pmin_r = p_r[j] - 0.5*rho0_r*c0_r*c0_r/G0_r;
	    if (pmin_r < 0.0)
	    {
	        vmin_r = v0_r - (c0_r-sqrt(c0_r*c0_r-2.0*G0_r*p0_r/rho0_r));
	        pmin_r = 0.0;
	    }
	    else
	        vmin_r = v_r[j] - c0_r/G0_r;
	    if (G0_r > 2.0)
	    {
	        pmax_r = HUGE_VAL;
		vmax_r = HUGE_VAL;
	    }
	    else
	    {
		vmax_r = v0_r + c0_r/(1.0 - 0.5*G0_r);
	        pmax_r = p0_r + rho0_r*sqr(vmax_r - v0_r);
	    }
	    vm[j] = 0.5*(v_l[j]+v_r[j]);
	    for (k = 0; k < N; ++k)
	    {
	        dv_l = (vm[j]-v0_l);
	        dv_r = (vm[j]-v0_r);
	        vm[j] -= ( (p0_r+dv_r*(0.5*rG0_r*dv_r+i0_r)) -
		           (p0_l+dv_l*(0.5*rG0_l*dv_l-i0_l)) )/
			 ( (rG0_r*dv_r+i0_r) - (rG0_l*dv_l-i0_l) );
	    }
	    dv_l = (vm[j]-v0_l);
	    dv_r = (vm[j]-v0_r);
	    pm[j] = 0.5*( (p0_r+dv_r*(0.5*rG0_r*dv_r+i0_r)) +
		          (p0_l+dv_l*(0.5*rG0_l*dv_l-i0_l)) );
	    if ((pm[j] < pmin_l) || (pm[j] < pmin_r) ||
	        (pm[j] > pmax_l) || (pm[j] > pmax_r) ||
	        (vm[j] < vmin_l) || (vm[j] < vmin_r) ||
	        (vm[j] > vmax_l) || (vm[j] > vmax_r))
	    {
	        g_exact_Riemann_midstate(j,j+1,vlst,vrst,pm,vm,vmuscl);
	    }
	}
}		/*end g_linear_us_up_Riemann_midstate*/

LOCAL	void	load_derived_fields(
	int start,
	int end,
	Vec_Gas *vlst,
	Vec_Gas *vrst)
{
	double *c2, *c;
    	int i;

	load_pressure_and_gammas(vlst,start,end-start);
	if (!Vec_Gas_field_set(vlst,c))
	{
	    c2 = vlst->c2;
	    c = vlst->c;
	    for (i = start; i < end; ++i)
		c[i] = sqrt(c2[i]);
	    Vec_Gas_field_set(vlst,c) = YES;
	}
	load_pressure_and_gammas(vrst,start,end-start);
	if (!Vec_Gas_field_set(vrst,c))
	{
	    c2 = vrst->c2;
	    c = vrst->c;
	    for (i = start; i < end; ++i)
		c[i] = sqrt(c2[i]);
	    Vec_Gas_field_set(vrst,c) = YES;
	}
}		/*end load_derived_fields*/
