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
*				gmoc.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

	/* LOCAL Function Prototypes */
LOCAL	void adjust_for_strong_rarefaction(double,double,double,double,double,double,
                                           double*,double*,Locstate);
LOCAL	void first_order_moc(int,int,Vec_Muscl*);
LOCAL	void riemann_moc(int,int,Vec_Muscl*);

/*ARGSUSED*/
EXPORT	void	g_riemann_characteristic_solve(
	int       start,
	int       end,
	Vec_Muscl *vmuscl)
{
	riemann_moc(start,end,vmuscl);
	g_load_muscl_flux(start,end,vmuscl->uM1,&vmuscl->VM1st,
	                  &vmuscl->Flux1,vmuscl);
}		/*end g_riemann_characteristic_solve*/

LOCAL	void	riemann_moc(
	int       start,
	int       end,
	Vec_Muscl *vmuscl)
{
	RECT_GRID *gr = vmuscl->front->rect_grid;
	Vec_Gas   *vst = vmuscl->vst;
	Vec_Gas   *vlst = &vmuscl->VLst;
	Vec_Gas   *vmst = &vmuscl->VMst;
	Vec_Gas   *vrst = &vmuscl->VRst;
	Vec_Gas   *vans = &vmuscl->VM1st;
	double     *p_l = vlst->p, *p_m = vmst->p;
	double     *p_r = vrst->p, *p_a = vans->p;
	double     *c_a = vans->c;
	double     **v_l = vlst->v, **v_m = vmst->v;
	double     **v_r = vrst->v, **v_a = vans->v;
	double     *rho_l = vlst->rho, *rho_m = vmst->rho;
	double     *rho_r = vrst->rho, *rho_a = vans->rho;
	double     *g_l = vmuscl->uL[vmuscl->index.grav];
	double     *g_r = vmuscl->uR[vmuscl->index.grav];
	double     *g = vmuscl->uM1[vmuscl->index.grav];
	double     *c_l = vlst->c, *c_m = vmst->c, *c_r = vrst->c;
	double     *e_m = vmst->e, *e_a = vans->e;
	double     *FD_m, *FD_a;
	double     p0_l, p0_r, p0_m, v0_l, v0_r;
	double     c0_l, c0_m, c0_r, i0_l, i0_r;
	double     g0_l, g0_r;
	double     rho0_l, rho0_r, rho0_m;
	double     b0_m, c20_m;
	double     r, r0, r0_l, r0_r;
	double     W_l, W_r;
	double     dt = vmuscl->dt;
	double     p, dp, v, b, c, c2, rho;
	double     dv_l, dv_r;
	double     alpha = vmuscl->alpha;
	const double eps = MACH_EPS;
	const double MNlim = 1.5;
	int       dim = vmuscl->dim;
	int       j, k;

	if (vmuscl->index.FD >= 0)
	{
	    FD_m = vmst->FD;
	    FD_a = vans->FD;
	}
	else
	{
	    FD_m = FD_a = NULL;
	}
	rmidstate(start,end,vlst,vrst,p_a,v_a[0],vmuscl);
	for (j = start; j < end; ++j)
	{
	    p = p_a[j];
	    p0_l = p_l[j];     p0_m = p_m[j];     p0_r = p_r[j];
	    rho0_l = rho_l[j]; rho0_m = rho_m[j]; rho0_r = rho_r[j];

	    c0_m = c_m[j];
	    rho0_m = rho_m[j];
	    for (k = 1; k < dim; ++k)
	        v_a[k][j] = v_m[k][j];
	    dp = p - p0_m;
	    c20_m = c0_m*c0_m;
	    b0_m = rho0_m*c20_m;
	    if (FD_m != NULL)
	    {
	        double dV = -dp/(rho0_m*b0_m) + FD_m[j]*dp*dp/(rho0_m*b0_m*b0_m);

	        rho_a[j] = rho0_m/(1.0 + rho0_m*dV);
		e_a[j] = e_m[j] - 0.5*(p+p0_m)*dV;

		if ((dV < 0.0) && (MNlim*fabs(dp) < rho0_m*b0_m*fabs(dV)))
		{
		    double rho_tmp = rho0_m + dp/c20_m; 
		    if ((0.0 < rho_tmp) && (rho_tmp < rho_a[j]))
		    {
			rho_a[j] = rho_tmp;
	                e_a[j] = e_m[j] + dp*p0_m/(rho0_m*b0_m);
		    }
		}
	    }
	    else
	    {
	        rho_a[j] = rho0_m + dp/c20_m; 
	        e_a[j] = e_m[j] + dp*p0_m/(rho0_m*b0_m);
	    }
	    adjust_for_strong_rarefaction(rho0_m,e_m[j],p,p0_l,p0_m,p0_r,
	                                  rho_a+j,e_a+j,vst->state[j]);
	}
	Vec_Gas_field_set(vans,p) = YES;
	Vec_Gas_field_set(vans,rho) = YES;
	Vec_Gas_field_set(vans,e) = YES;
	Vec_Gas_field_set(vans,v) = YES;
	Vec_Gas_field_set(vans,c) = NO;
	load_sound_speed(vans,start,end-start);
	for (j = start; j < end; ++j)
	{
	      p0_l = p_l[j];     p0_m = p_m[j];     p0_r = p_r[j];
	      c0_l = c_l[j];     c0_m = c_m[j];     c0_r = c_r[j];
	    rho0_l = rho_l[j]; rho0_m = rho_m[j]; rho0_r = rho_r[j];

	    c20_m = c0_m*c0_m;
	    b0_m = rho0_m*c20_m;

	    v0_l = v_l[0][j];
	    v0_r = v_r[0][j];
	    g0_l = 0.5*(g[j] + g_l[j]);
	    g0_r = 0.5*(g[j] + g_r[j]);


	    rho = rho_a[j];
	    p = p_a[j];
	    v = v_a[0][j];
	    c = c_a[j];
	    c2 = c*c;
	    b = rho*c2;

	    dv_l = v - v0_l;
	    dv_r = v - v0_r;
	    i0_l = rho0_l*c0_l;
	    if (fabs(dv_l) > eps) /*TOLERANCE*/
	    {
	        W_l = fabs((p - p0_l)/dv_l);
		if (W_l/i0_l < eps)
		    W_l = i0_l;
	    }
	    else
	        W_l = i0_l;
	    i0_r = rho0_r*c0_r;
	    if (fabs(dv_r) > eps) /*TOLERANCE*/
	    {
	        W_r = fabs((p - p0_r)/dv_r);
		if (W_r/i0_r < eps)
		    W_r = i0_r;
	    }
	    else
	        W_r = i0_r;

	    v0_l = v0_l + g0_l*dt;
	    v0_r = v0_r + g0_r*dt;

	    if (alpha != 0.0)
	    {
	        r0_l = pos_radius(vlst->coords[j][0],gr);
	        r0_r = pos_radius(vrst->coords[j][0],gr);
		r0 = vans->coords[j][0];
	        r = pos_radius(0.0,gr);
		if (fabs(r0) < fabs(r)) /* At Origin */
		{
		    if (r0 > 0.0)
		    {
		        W_l    =    W_r;
		        v0_l   =  -v0_r;
		        p0_l   =   p0_r;
		        c0_l   =   c0_r;
		        rho0_l = rho0_r;
		    }
		    else
		    {
		        W_r    =    W_l;
		        v0_r   =  -v0_l;
		        p0_r   =   p0_l;
		        c0_r   =   c0_l;
		        rho0_r = rho0_l;
		    }
		}
		else
		{
	            v0_l -= 0.5*alpha*(v*c/r0 + v0_l*c0_l/r0_l)*dt;
	            v0_r += 0.5*alpha*(v*c/r0 + v0_r*c0_r/r0_r)*dt;
	        }
	    }

	    p = (p0_l*W_r+p0_r*W_l+W_l*W_r*(v0_l-v0_r))/(W_l+W_r);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (p < Min_pressure(vst->state[j]))
	        p = Min_pressure(vst->state[j]);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    p_a[j] = p;
	    v_a[0][j] = (p0_l - p0_r + W_l*v0_l + W_r*v0_r)/(W_l + W_r);
	    dp = p - p0_m;
	    if (FD_m != NULL)
	    {
	        double A, B, dV;
		A = -0.5*(1.0/(rho0_m*b0_m) + 1.0/(rho*b));
		B =  0.5*(FD_m[j]/(rho0_m*b0_m*b0_m) + FD_a[j]/(rho*b*b));
	        dV = A*dp + B*dp*dp;
	        rho_a[j] = rho0_m/(1.0 + rho0_m*dV);
		e_a[j] = e_m[j] - 0.5*(p+p0_m)*dV;

		if ((dV < 0.0) && (MNlim*fabs(dp) < rho0_m*b0_m*fabs(dV)))
		{
		    double rho_tmp = rho0_m + dp/c20_m; 
		    if ((0.0 < rho_tmp) && (rho_tmp < rho_a[j]))
		    {
			rho_a[j] = rho_tmp;
	                e_a[j] = e_m[j] + dp*p0_m/(rho0_m*b0_m);
		    }
		}
	    }
	    else
	    {
	        rho_a[j] = rho0_m + 2.0*dp/(c*c + c20_m); 
	        e_a[j] = e_m[j] + 0.5*(p/sqr(rho*c) + p0_m/sqr(rho0_m*c0_m))*dp;
	    }
	    adjust_for_strong_rarefaction(rho0_m,e_m[j],p_a[j],p0_l,p0_m,p0_r,
	                                  rho_a+j,e_a+j,vst->state[j]);
	}
	Vec_Gas_field_set(vans,p) = YES;
	Vec_Gas_field_set(vans,rho) = YES;
	Vec_Gas_field_set(vans,e) = YES;
	Vec_Gas_field_set(vans,v) = YES;
	Vec_Gas_field_set(vans,c) = NO;
	load_sound_speed(vans,start,end-start);
}		/*end riemann_moc*/

/*ARGSUSED*/
EXPORT	void	g_implicit_characteristic_solve(
	int       start,
	int       end,
	Vec_Muscl *vmuscl)
{
	RECT_GRID *gr = vmuscl->front->rect_grid;
	Vec_Gas   *vst = vmuscl->vst;
	Vec_Gas   *vlst = &vmuscl->VLst;
	Vec_Gas   *vmst = &vmuscl->VMst;
	Vec_Gas   *vrst = &vmuscl->VRst;
	Vec_Gas   *vans = &vmuscl->VM1st;
	double     *FD_m, *FD_a;
	double     *p_l = vlst->p, *p_m = vmst->p;
	double     *p_r = vrst->p, *p_a = vans->p;
	double     **v_l = vlst->v, **v_m = vmst->v;
	double     **v_r = vrst->v, **v_a = vans->v;
	double     *rho_l = vlst->rho, *rho_m = vmst->rho;
	double     *rho_r = vrst->rho, *rho_a = vans->rho;
	double     *c_l = vlst->c, *c_m = vmst->c, *c_r = vrst->c;
	double     *c_a = vans->c;
	double     *g_l = vmuscl->uL[vmuscl->index.grav];
	double     *g_r = vmuscl->uR[vmuscl->index.grav];
	double     *g = vmuscl->uM1[vmuscl->index.grav];
	double     *e_m = vmst->e, *e_a = vans->e;
	double     p0_l, p0_r, p0_m, v0_l, v0_m, v0_r;
	double     g0_l, g0_r;
	double     rho0_l, rho0_r, rho0_m;
	double     c0_l, c0_m, c0_r;
	double     c20_m, b0_m;
	double     r, r0, r0_l, r0_r;
	double     W_l, W_r;
	double     dt = vmuscl->dt;
	double     p, v, b, c, c2, rho;
	double     alpha;
	double     dp, dv;
	double     A, B, dV;
	int       j, n;
	const int N = 4;

	first_order_moc(start,end,vmuscl);

	if (vmuscl->index.FD >= 0)
	{
	    FD_m = vmst->FD;
	    FD_a = vans->FD;
	}
	else
	{
	    FD_m = FD_a = NULL;
	}

	for (j = start; j < end; ++j)
	{
	      p0_l = p_l[j];     p0_m = p_m[j];     p0_r = p_r[j];
	      c0_l = c_l[j];     c0_m = c_m[j];     c0_r = c_r[j];
	    rho0_l = rho_l[j]; rho0_m = rho_m[j]; rho0_r = rho_r[j];

	    c20_m = c0_m*c0_m;
	     b0_m = rho0_m*c20_m;

	    g0_l = 0.5*(g[j] + g_l[j]);
	    g0_r = 0.5*(g[j] + g_r[j]);
	    v0_l = v_l[0][j] + g0_l*dt;
	    v0_m = v_m[0][j];
	    v0_r = v_r[0][j] + g0_r*dt;
	    if (vmuscl->alpha != 0.0)
	    {
	        r0_l = pos_radius(vlst->coords[j][0],gr);
	        r0_r = pos_radius(vrst->coords[j][0],gr);
		r0 = vans->coords[j][0];
	        r = pos_radius(r0,gr);
		if (r0 < r) /* At origin */
		{
		      g0_l =   g0_r;
		      v0_l =  -v0_r;
		      p0_l =   p0_r;
		      c0_l =   c0_r;
		    rho0_l = rho0_r;
		     alpha = 0.0;
		}
		else
		    alpha = vmuscl->alpha;
	    }

	    for (n = 0; n < N; ++n)
	    {
		rho = rho_a[j];
		v = v_a[0][j];
		c = c_a[j];
		c2 = c*c;
		b = rho*c2;
		W_l = 0.5*(rho*c + rho0_l*c0_l);
		W_r = 0.5*(rho*c + rho0_r*c0_r);
		  
		if (alpha != 0.0)
		{
	            v0_l = v_l[0][j] + g0_l*dt -
		           0.5*alpha*(v*c/r + v_l[0][j]*c0_l/r0_l)*dt;
	            v0_r = v_r[0][j] + g0_r*dt +
		           0.5*alpha*(v*c/r + v_r[0][j]*c0_r/r0_r)*dt;
		}
	        p = (p0_l*W_r+p0_r*W_l+W_l*W_r*(v0_l-v0_r))/(W_l+W_r);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	        if (p < Min_pressure(vst->state[j]))
	            p = Min_pressure(vst->state[j]);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
		dp = p - p0_m;
	        v = (p0_l - p0_r + W_l*v0_l + W_r*v0_r)/(W_l + W_r);
		dv = v - v0_m;
		if (fabs(dv/c0_m) > 10.0)/*TOLERANCE*/
		{
		    riemann_moc(j,j+1,vmuscl);
		    break;
		}
		p_a[j] = p;
		v_a[0][j] = v;
	        if (FD_m != NULL)
	        {
		    A = -0.5*(1.0/(rho0_m*b0_m) + 1.0/(rho*b));
		    B =  0.5*(FD_m[j]/(rho0_m*b0_m*b0_m) + FD_a[j]/(rho*b*b));
	            dV = A*dp + B*dp*dp;
	            rho_a[j] = rho0_m/(1.0 + rho0_m*dV);
		    e_a[j] = e_m[j] - 0.5*(p+p0_m)*dV;
	        }
	        else
	        {
	            rho = rho0_m + 2.0*dp/(c2 + c20_m); 
	            e_a[j] = e_m[j] + 0.5*(p/(rho*b) + p0_m/(rho0_m*b0_m))*dp;
		}
	        adjust_for_strong_rarefaction(rho0_m,e_m[j],p,p0_l,p0_m,p0_r,
		                              &rho,e_a+j,vst->state[j]);
		rho_a[j] = rho;
	        Vec_Gas_field_set(vans,rho) = YES;
	        Vec_Gas_field_set(vans,e) = YES;
	        Vec_Gas_field_set(vans,p) = YES;
	        Vec_Gas_field_set(vans,v) = YES;
	        Vec_Gas_field_set(vans,c) = NO;
	        load_sound_speed(vans,j,1);
	    }
	}
	g_load_muscl_flux(start,end,vmuscl->uM1,vans,&vmuscl->Flux1,vmuscl);
}		/*end g_implicit_characteristic_solve*/

/*ARGSUSED*/
EXPORT	void	g_first_order_direct_characteristic_solve(
	int       start,
	int       end,
	Vec_Muscl *vmuscl)
{
	first_order_moc(start,end,vmuscl);
	g_load_muscl_flux(start,end,vmuscl->uM1,&vmuscl->VM1st,
	                  &vmuscl->Flux1,vmuscl);
}		/*end g_first_order_direct_characteristic_solve*/

/*ARGSUSED*/
LOCAL	void	first_order_moc(
	int       start,
	int       end,
	Vec_Muscl *vmuscl)
{
	RECT_GRID *gr = vmuscl->front->rect_grid;
	Vec_Gas   *vst = vmuscl->vst;
	Vec_Gas   *vlst = &vmuscl->VLst;
	Vec_Gas   *vmst = &vmuscl->VMst;
	Vec_Gas   *vrst = &vmuscl->VRst;
	Vec_Gas   *vans = &vmuscl->VM1st;
	double     *p_l = vlst->p, *p_m = vmst->p;
	double     *p_r = vrst->p, *p_a = vans->p;
	double     **v_l = vlst->v, **v_m = vmst->v;
	double     **v_r = vrst->v, **v_a = vans->v;
	double     *g_l = vmuscl->uL[vmuscl->index.grav];
	double     *g_r = vmuscl->uR[vmuscl->index.grav];
	double     *rho_l = vlst->rho, *rho_m = vmst->rho;
	double     *rho_r = vrst->rho, *rho_a = vans->rho;
	double     *c_l = vlst->c, *c_m = vmst->c, *c_r = vrst->c;
	double     *e_m = vmst->e, *e_a = vans->e;
	double     p0_l, p0_r, p0_m, v0_l, v0_r;
	double     c0_l, c0_m, c0_r;
	double     c20_m;
	double     g0_l, g0_r;
	double     rho0_l, rho0_r, rho0_m;
	double     r, r0, r0_l, r0_r;
	double     dt = vmuscl->dt;
	double     p, dp;
	double     alpha = vmuscl->alpha;
	int       dim = vmuscl->dim;
	int       j, k;

	for (j = start; j < end; ++j)
	{
	      p0_l = p_l[j];     p0_m = p_m[j];     p0_r = p_r[j];
	      c0_l = c_l[j];     c0_m = c_m[j];     c0_r = c_r[j];
	    rho0_l = rho_l[j]; rho0_m = rho_m[j]; rho0_r = rho_r[j];
	      g0_l = g_l[j];     g0_r = g_r[j];

	    c20_m = c0_m*c0_m;

	    v0_l = v_l[0][j] + g0_l*dt;
	    v0_r = v_r[0][j] + g0_r*dt;
	    if (alpha != 0.0)
	    {
	        r0_l = pos_radius(vlst->coords[j][0],gr);
	        r0_r = pos_radius(vlst->coords[j][0],gr);
		r0 = vans->coords[j][0];
	        r = pos_radius(r0,gr);
		if (r0 < r) /* At origin */
		{
		      v0_l =  -v0_r;
		      p0_l =   p0_r;
		      c0_l =   c0_r;
		    rho0_l = rho0_r;
		}
		else
		{
	            v0_l -= v_l[0][j]*alpha*c0_l*dt/r0_l;
	            v0_r += v_r[0][j]*alpha*c0_r*dt/r0_r;
		}
	    }

	    p = (p0_l*rho0_r*c0_r + p0_r*rho0_l*c0_l +
	                    rho0_r*c0_r*rho0_l*c0_l*(v0_l - v0_r))/
		         (rho0_l*c0_l + rho0_r*c0_r);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (p < Min_pressure(vst->state[j]))
	        p = Min_pressure(vst->state[j]);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    p_a[j] = p;
	    v_a[0][j] = (p0_l - p0_r + rho0_l*c0_l*v0_l + rho0_r*c0_r*v0_r)/
		            (rho0_l*c0_l + rho0_r*c0_r);
	    dp = p - p0_m;
	    rho_a[j] = rho0_m + dp/c20_m; 
	    e_a[j] = e_m[j] - 0.5*(p+p0_m)*(1.0/rho_a[j] - 1.0/rho0_m);
	    adjust_for_strong_rarefaction(rho0_m,e_m[j],p,p0_l,p0_m,p0_r,
	                                  rho_a+j,e_a+j,vst->state[j]);
	    for (k = 1; k < dim; ++k)
	        v_a[k][j] = v_m[k][j];
	}
	Vec_Gas_field_set(vans,p) = YES;
	Vec_Gas_field_set(vans,rho) = YES;
	Vec_Gas_field_set(vans,e) = YES;
	Vec_Gas_field_set(vans,v) = YES;
	Vec_Gas_field_set(vans,c) = NO;
	load_sound_speed(vans,start,end-start);
}		/*end first_order_moc*/


LOCAL	void adjust_for_strong_rarefaction(
	double     rho0,
	double     e0,
	double     p_a,
	double     p0_l,
	double     p0_m,
	double     p0_r,
	double     *rho_a,
	double     *e_a,
	Locstate  state)
{
	boolean strong_rarefaction = (*rho_a < 0.0) ? YES : NO;

	if (p_a < p0_m)
	{
	    double p_min = min(p0_l,p0_r);
	    p_min = min(p_min,p0_m);
	    if (fabs(p_a/p_min) > 10.0)/*TOLERANCE*/
	        strong_rarefaction = YES;
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (fabs(*rho_a)*(*e_a) < Min_energy(state))
	    strong_rarefaction = YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	if (strong_rarefaction)
	{
	    static	  Locstate st0 = NULL, sta = NULL;
	    if (st0 == NULL)
	    {
	        (*Params(state)->_alloc_state)(&st0,Params(state)->sizest);
	        (*Params(state)->_alloc_state)(&sta,Params(state)->sizest);
	    }
	    Dens(st0) = rho0;
	    Energy(st0) = e0;
	    Set_params(st0,state);
	    set_type_of_state(st0,EGAS_STATE);
	    reset_gamma(st0);
	    state_on_adiabat_with_pr(st0,p_a,sta,EGAS_STATE);
	    *rho_a = Dens(sta);
	    *e_a = Energy(sta);
	}
}		/*end adjust_for_strong_rarefaction*/
