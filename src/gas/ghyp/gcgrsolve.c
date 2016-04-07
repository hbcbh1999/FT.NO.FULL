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
*                	gcgrsolve.c
*
*    Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    An implementation of the Colella-Glaz approximate Riemann
*    solver.  See J. Comp. Physics 59, 264-289 (1985).
*
*/


#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

/*
*    The cg_data structure contains all the information
*    used by the Riemann solver
*/

typedef struct
{

	Vec_Muscl *vmuscl;
	Vec_Gas *vlst, *vrst, *vmst; /* left, right and mid states for */

	double *left_Gam, *right_Gam;/* left and right state Gammas */
	double *mid_Gam;             /* midstate Gamma */
	double *Gam_min, *Gam_max;   /* Extreme values of Gamma in the
	                             * neighborhood of a cell */
	double *pinfl, *pinfr; /* p_infinities for left and right states */
	double *Wl, *Wr;       /* Mass fluxes through left and right waves */
	double *pnew, *pold;   /* Last two pressures in secant iteration */
	double *vlnew, *vlold; /* Last two left velocities in secant iteration */
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double *min_delta_p;   /* minimum pressure tolerance */
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	double *vrnew, *vrold; /* Last two right velocities */
	double **worksp;
	MUSCL_FLUX *flux;	/* Conservative flux across cell boundaries */
	CG_PARAMS        Cg_params;
	int rsten;            /* Radius of finite difference stencil (2 for
	                       * MUSCL) */
	int start;            /* Starting index for calculations */
	int end;              /* Ending index for calculations */
	int max_vec_size;     /* The largest vector we've seen--space
	                       * allocated */
	int dim;              /* Spatial dimension of problem */
	size_t sizest;        /* Space to allocate for an EGAS structure */
} CG_data;


	/* LOCAL Function Declarations */
LOCAL    CG_data    *alloc_rsolve_data(int);
LOCAL    CG_data    *load_cg_rdata(int,int,Vec_Gas*,Vec_Gas*,Vec_Gas*,
                                   MUSCL_FLUX*,Vec_Muscl*);
LOCAL	double	Mfluxsqr(double,double,double,double,double,double,double,double,double);
LOCAL	double	DGam_star(double,double,double,double,double);
LOCAL	double	find_mid_c_rho_Gam(double,double*,double*,double,double,double,
				   Locstate,double,double,double,double,double,
				   double,double,double,double,double,CG_data*);
LOCAL    void    cell_bdy_in_rarefaction(double,double,double,double,double,double,
					 double*,double*,double*,double*);
LOCAL    void    correct_for_strong_waves(int,CG_data*);
LOCAL    void    find_cons_flux(int,int,CG_data*);
LOCAL    void    find_Gam_min_and_max(CG_data*);
LOCAL    void    find_mass_fluxes(CG_data*,double*);
LOCAL    void    find_midstates(CG_data*);
LOCAL    void    find_Gam_and_p_infinities(CG_data*);
LOCAL    void    free_rsolve_data(CG_data*);
LOCAL    void    get_first_guesses(CG_data*);
LOCAL    void    iterate_for_mid_pv(CG_data*);


/*ARGSUSED*/
EXPORT	void	g_cg_Riemann_midstate(
	int       start,
	int       end,
	Vec_Gas   *vlst,
	Vec_Gas   *vrst,
	double     *pm,
	double     *vm,
	Vec_Muscl *vmuscl)
{
	CG_data *rdata = NULL;
	Vec_Gas Vmst;

#if defined(DEBUG_MUSCL)
	debug_print("cg_rsolve","Entered cg_rsolve()\n");
#endif /* defined(DEBUG_MUSCL) */

	/*
	* The size of the current vector is end - start.  If this is bigger
	* than any of the uni_arrays we've allocated space for, or if we have never
	* allocated space, give up any current space and make some more.
	*/

	rdata = load_cg_rdata(start,end,vlst,vrst,NULL,NULL,vmuscl);
	zero_scalar(&Vmst,sizeof(Vec_Gas));
	Vmst.v = &vm;
	Vmst.p = pm;
	rdata->vmst = &Vmst;

	/*
	* This is where most of the work is.  Using the left and right
	* states use secant iteration to find the unique intersection
	* of wave curves corresponding to the midstate pressure and velocity
	*/

	iterate_for_mid_pv(rdata);
}		/*end g_cg_Riemann_midstate*/

/*ARGSUSED*/
EXPORT void cg_rsolve(
	int        start,
	int        end,
	double      **uL,
	Vec_Gas    *vlst,
	double      **uR,
	Vec_Gas    *vrst,
	double      **uM,
	Vec_Gas    *vmst,
	MUSCL_FLUX *flux,
	Vec_Muscl  *vmuscl)
{
	CG_data        *rdata = NULL;

#if defined(DEBUG_MUSCL)
	debug_print("cg_rsolve","Entered cg_rsolve()\n");
#endif /* defined(DEBUG_MUSCL) */

	/*
	* The size of the current vector is end - start.  If this is bigger
	* than any of the uni_arrays we've allocated space for, or if we have never
	* allocated space, give up any current space and make some more.
	*/

	rdata = load_cg_rdata(start,end,vlst,vrst,vmst,flux,vmuscl);

	/*
	* This is where most of the work is.  Using the left and right
	* states use secant iteration to find the unique intersection
	* of wave curves corresponding to the midstate pressure and velocity
	*/

	iterate_for_mid_pv(rdata);

	/*
	* We have the midstate pressure and velocity.  Now find the density,
	* energy, and momenta.
	*/

	find_midstates(rdata);

	/*
	* All right, now we know everything about the midstate, so find
	* the conserved-quantities flux.
	*/

	find_cons_flux(start,end,rdata);

#if defined(DEBUG_MUSCL)
	debug_print("cg_rsolve","Left cg_rsolve()\n");
#endif /* defined(DEBUG_MUSCL) */
}	        /*end cg_rsolve*/

LOCAL    CG_data *load_cg_rdata(
	int        start,
	int        end,
	Vec_Gas    *vlst,
	Vec_Gas    *vrst,
	Vec_Gas    *vmst,
	MUSCL_FLUX *flux,
	Vec_Muscl  *vmuscl)
{
	static CG_data    *rdata = NULL;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double *min_delta_p;
	double *min_pressure = vmuscl->vst->min_pressure + vmuscl->offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int   vsize = vmuscl->vsize;
	int   vectr_size = end - start;
	int   j;

	if ((rdata == NULL) || (vsize > rdata->max_vec_size))
	{
	    free_rsolve_data(rdata);
	    rdata = alloc_rsolve_data(vsize);
	}

	rdata->vmuscl = vmuscl;
	rdata->vlst = vlst;
	rdata->vrst = vrst;
	rdata->vmst = vmst;
	rdata->sizest = vmuscl->sizest;
	rdata->flux = flux;
	rdata->start = start;
	rdata->end = end;
	rdata->dim = vmuscl->dim;
	rdata->rsten = vmuscl->sten_rad;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	min_delta_p = rdata->min_delta_p;
	for (j = start; j < end; ++j)
	    min_delta_p[j] = 100.0*min_pressure[j];/*TOLERANCE*/
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

/*    
*    Fill in the pressure, sound speed, and Gruneisen coefficient
*       fields for the right and left states
*/

	load_pressure_and_gammas(rdata->vlst, start,vectr_size);
	load_pressure_and_gammas(rdata->vrst, start,vectr_size);

/*    Calculate the parameters p_infinity and Gamma    */

	find_Gam_and_p_infinities(rdata);


/*    Now find the values of Gam_min, Gam_max.  They are used to limit
	the Gamma calculated at a midstate */

	find_Gam_min_and_max(rdata);

	return rdata;
}	        /*end load_cg_rdata*/

/*
*            find_Gam_and_p_infinities():
*
*    Calculates p_infinity for each right and left state by the formula
*
*    p_infinity      = rho*c2/(GAM+1) - p
*    Gam        = (p + p_infinity)/(e*rho - p_infinity)
*/


LOCAL void find_Gam_and_p_infinities(
	CG_data        *rdata)
{
	int        j;
	int        start = rdata->start, end = rdata->end;
	double        *pinf_l = rdata->pinfl;
	double        *pinf_r = rdata->pinfr;
	double        *c2_l = rdata->vlst->c2;
	double        *c2_r = rdata->vrst->c2;
	double        *GAM_l = rdata->vlst->GAM;
	double        *GAM_r = rdata->vrst->GAM;
	double        *Gam_l = rdata->left_Gam;
	double        *Gam_r = rdata->right_Gam;
	double        *p_l = rdata->vlst->p;
	double        *rho_l = rdata->vlst->rho;
	double        *e_l = rdata->vlst->e;
	double        *p_r = rdata->vrst->p;
	double        *rho_r = rdata->vrst->rho;
	double        *e_r = rdata->vrst->e;

#define cg_pinf(p,c2,rho,GAM)    ((rho)*(c2)/(GAM+1.0) - (p))
	for (j = start; j < end; ++j)
	{
	    pinf_l[j] = cg_pinf(p_l[j],c2_l[j],rho_l[j],GAM_l[j]);
	    pinf_r[j] = cg_pinf(p_r[j],c2_r[j],rho_r[j],GAM_r[j]);
	}
#undef cg_pinf

	for (j = start; j < end; ++j)
	{
	    pinf_l[j] = max(pinf_l[j], 0.0);
	    pinf_r[j] = max(pinf_r[j], 0.0);
	}

#define cg_Gam(p,e,rho,pinf)    ((p) + (pinf))/((e)*(rho) - (pinf))
	for (j = start; j < end; ++j)
	{
	    Gam_l[j] = cg_Gam(p_l[j],e_l[j],rho_l[j],pinf_l[j]);
	    Gam_r[j] = cg_Gam(p_r[j],e_r[j],rho_r[j],pinf_r[j]);
	}
#undef cg_Gam
}	        /*end find_Gam_and_p_infinities*/


/*
*            find_Gam_min_and_max():
*
*    Finds extreme values of Gamma.  For each j it looks at a
*    neighborhood of cells with radius rsten centered at cell j.
*/


LOCAL void find_Gam_min_and_max(
	CG_data        *rdata)
{
	int        start = rdata->start, end = rdata->end;
	int        rsten = rdata->rsten;
	double      *Gam_max = rdata->Gam_max;
	double      *Gam_min = rdata->Gam_min;
	double      *Gaml = rdata->left_Gam;
	double      *Gamr = rdata->right_Gam;
	double      gam_min, gam_max;
	int        offset, j;

	for (j = start; j < end; ++j)
	    Gam_min[j] = Gam_max[j] = Gaml[j];

	for (offset = 1; offset <= rsten; ++offset)
	{
	    for (j = start+offset; j < end; ++j)
	    {
	        if (Gaml[j - offset] < Gamr[j - offset])
	        {
	            gam_min = Gaml[j - offset];
	            gam_max = Gamr[j - offset];
	        }
	        else
	        {
	            gam_min = Gamr[j - offset];
	            gam_max = Gaml[j - offset];
	        }
	        Gam_min[j] = min(Gam_min[j],gam_min);
	        Gam_max[j] = max(Gam_max[j],gam_max);
	    }
	    for (j = start; j < end-offset; ++j)
	    {
	        if (Gaml[j + offset] < Gamr[j + offset])
	        {
	            gam_min = Gaml[j + offset];
	            gam_max = Gamr[j + offset];
	        }
	        else
	        {
	            gam_min = Gamr[j + offset];
	            gam_max = Gaml[j + offset];
	        }
	        Gam_min[j] = min(Gam_min[j],gam_min);
	        Gam_max[j] = max(Gam_max[j],gam_max);
	    }
	}
}	        /*end find_Gam_min_and_max*/


/*
*            iterate_for_mid_pv():
*
*    Finds the intersection of the right and left state wave curves.
*    When finished it fills in the pressure and velocity uni_arrays.
*    It also returns with the uni_arrays Wl and Wr updated to the values based
*    on the most accurate pressures.
*/


LOCAL void iterate_for_mid_pv(
	CG_data        *rdata)
{
	Vec_Gas *vlst = rdata->vlst;
	Vec_Gas *vrst = rdata->vrst;
	Vec_Gas *vmst = rdata->vmst;
	int   start = rdata->start, end = rdata->end;
	int   pv_iterations = rdata->Cg_params.pv_iterations;
	int   j, n;
	double ptemp;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double *min_pressure = rdata->vmuscl->vst->min_pressure +
	                      rdata->vmuscl->offset;
	double *min_delta_p = rdata->min_delta_p;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	double *pnew = rdata->pnew;
	double *pold = rdata->pold;
	double *vlnew = rdata->vlnew;
	double *vlold = rdata->vlold;
	double *vrnew = rdata->vrnew;
	double *vrold = rdata->vrold;
	double *left_vn = vlst->v[0];
	double *left_p = vlst->p;
	double *right_vn = vrst->v[0];
	double *right_p = vrst->p;
	double *mid_vn = vmst->v[0];
	double *mid_p = vmst->p;
	double *Wl = rdata->Wl;
	double *Wr = rdata->Wr;
	double dp, dvl, dvr;
	double min_p_jump = rdata->Cg_params.min_p_jump;
	double min_v_jump = rdata->Cg_params.min_v_jump;
	double den;

#if defined(DEBUG_MUSCL)
	debug_print("cg_rsolve","Entered iterate_for_mid_pv()\n");
#endif /* defined(DEBUG_MUSCL) */
	get_first_guesses(rdata);

	for (n = 0; n < pv_iterations; ++n)
	{
	    find_mass_fluxes(rdata, pnew);
	    for (j = start; j < end; ++j)
	    {
	        dp = fabs(pnew[j] - pold[j]);
	        if ((dp < min_p_jump*mid_p[j])
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	                        || dp < min_delta_p[j]
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        )
	            continue;

	        vlold[j] = vlnew[j];
	        vrold[j] = vrnew[j];
	        ptemp = pnew[j];

	        vlnew[j] =  left_vn[j] - (ptemp -  left_p[j])/Wl[j];
	        vrnew[j] = right_vn[j] + (ptemp - right_p[j])/Wr[j];

	        dp = ptemp - pold[j];
	        dvl = vlnew[j] - vlold[j];
	        dvr = vrnew[j] - vrold[j];
	        den = fabs(dvl) + fabs(dvr);
	        if (den < fabs(min_v_jump*(left_p[j]/Wl[j] + 
	                       right_p[j]/Wr[j])))
	        {
	            vlnew[j] = vlold[j];
	            vrnew[j] = vrold[j];
	        }
	        else
	        {
	            pnew[j] = ptemp -
	                (vrnew[j] - vlnew[j])*fabs(dp)/den;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	            if (pnew[j] < min_pressure[j])
	                pnew[j] = min_pressure[j];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	            pold[j] = ptemp;
	        }
	    }

	}

	find_mass_fluxes(rdata, pnew);
	for (j = start; j < end; ++j)
	{
	    mid_p[j] = pnew[j];
	    mid_vn[j] = 0.5*(
	             left_vn[j] - (mid_p[j] -  left_p[j])/Wl[j] +
	            right_vn[j] + (mid_p[j] - right_p[j])/Wr[j]);
	}
	Vec_Gas_field_set(vmst,p) = YES;
#if defined(DEBUG_MUSCL)
	if (debugging("cg_rsolve"))
	{
	    (void) printf("Mid state pressures and normal velocities\n");
	    (void) printf("%-4s %-14s %-14s\n","n","pressure","velocity");
	    for (j = start; j < end; ++j)
	        (void) printf("%-4d %-14g %-14g\n",j,mid_p[j],mid_vn[j]);
	    (void) printf("\n");
	}
	debug_print("cg_rsolve","Left iterate_for_mid_pv()\n");
#endif /* defined(DEBUG_MUSCL) */
}	        /*end iterate_for_mid_pv*/




/*
*            get_first_guesses():
*
*    Find two guesses for the midstate pressure to start the
*    secant iterations.  The first guess is the average of the right
*    and left pressures and the second is found via a single Godunov
*    fixed point iteration step.
*/

LOCAL void get_first_guesses(
	CG_data        *rdata)
{
	int   j;
	int   start = rdata->start, end = rdata->end;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double *min_pressure = rdata->vmuscl->vst->min_pressure +
	                      rdata->vmuscl->offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	double *first_p = rdata->pold;
	double *scnd_p = rdata->pnew;
	double *scnd_vl = rdata->vlnew;
	double *scnd_vr = rdata->vrnew;
	double *Wl = rdata->Wl;
	double *Wr = rdata->Wr;
	double *left_vn = rdata->vlst->v[0];
	double *left_p = rdata->vlst->p;
	double *right_vn = rdata->vrst->v[0];
	double *right_p = rdata->vrst->p;
	double *mid_p = rdata->vmst->p;

	for (j = start; j < end; ++j)
	    mid_p[j] = first_p[j] = 0.5*(left_p[j] + right_p[j]);

	find_mass_fluxes(rdata, first_p);
	for (j = start; j < end; ++j)
	{
	    scnd_p[j] = ((left_vn[j] - right_vn[j])*Wr[j]*Wl[j] +
	            left_p[j]*Wr[j] + right_p[j]*Wl[j]) /
	            ( Wr[j] + Wl[j]);
	    scnd_vl[j] =  left_vn[j] - (first_p[j] -  left_p[j])/Wl[j];
	    scnd_vr[j] = right_vn[j] + (first_p[j] - right_p[j])/Wr[j];
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	for (j = start; j < end; ++j)
	    if (scnd_p[j] < min_pressure[j])
	        scnd_p[j] = min_pressure[j];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
}	        /*end get_first_guesses*/



/*
*            find_mass_fluxes():
*
*    Calculates the mass flux W for the forward and backward facing
*    waves.  This takes into account the modified equations that use
*    p_infinity.
*/


LOCAL void find_mass_fluxes(
	CG_data    *rdata,
	double      *mid_p)
{
	int      start = rdata->start;
	int      end = rdata->end;
	double    *Gam_min = rdata->Gam_min;
	double    *Gam_max = rdata->Gam_max;
	double    *pinfl = rdata->pinfl;
	double    *pinfr = rdata->pinfr;
	double    *Wl = rdata->Wl;
	double    *Wr = rdata->Wr;
	double    *left_Gam = rdata->left_Gam;
	double    *right_Gam = rdata->right_Gam;

	double    *left_c2 = rdata->vlst->c2;
	double    *left_p = rdata->vlst->p;
	double    *left_rho = rdata->vlst->rho;

	double    *right_c2 = rdata->vrst->c2;
	double    *right_p = rdata->vrst->p;
	double    *right_rho = rdata->vrst->rho;
	double    *Gam_starr = Wr;
	double    *Gam_starl = Wl;
	double    *Wsql = Wl;
	double    *Wsqr = Wr;
	double    Gam_bar, gam_bar_i;
	double    dp, pbar;
	int      j;
	double    min_p_jump = rdata->Cg_params.min_p_jump;
	double    sqr_min_mass_flux = rdata->Cg_params.sqr_min_mass_flux;

	for (j = start; j < end; ++j)
	{
	    /*
	    * Come up with a first guess for the midstate Gamma using jump
	    * equation (27).  We derive a Gamma for the right side of the
	    * contact and one for the left side.
	    */


	    Gam_bar = 0.5*(right_Gam[j] + left_Gam[j]);
	    gam_bar_i = 2.0*left_p[j]*right_p[j]/
			    (right_c2[j]*right_rho[j]* left_p[j] +
			      left_c2[j]* left_rho[j]*right_p[j]);

	    dp = mid_p[j] - right_p[j];
	    pbar = 0.5*(mid_p[j] + right_p[j]);
	    Gam_starr[j] = right_Gam[j] + DGam_star(Gam_bar,gam_bar_i,dp,
						    pbar,pinfr[j]);

	    dp = mid_p[j] - left_p[j];
	    pbar = 0.5*(mid_p[j] + left_p[j]);
	    Gam_starl[j] = left_Gam[j] + DGam_star(Gam_bar,gam_bar_i,dp,
						   pbar,pinfl[j]);
	}

	/*    Now limit the Gammas using Gam_min and Gam_max */

	for (j = start; j < end; ++j)
	{
	    Gam_starr[j] = max(Gam_min[j], min(Gam_starr[j], Gam_max[j]));
	    Gam_starl[j] = max(Gam_min[j], min(Gam_starl[j], Gam_max[j]));
	}
	/*    Find W^2 and then take square root */

	for (j = start; j < end; ++j)
	{
	    Wsqr[j] = Mfluxsqr(mid_p[j],right_p[j],right_rho[j],min_p_jump,
	                       right_c2[j],right_Gam[j],Gam_starr[j],pinfr[j],
			       sqr_min_mass_flux);

	    Wsql[j] = Mfluxsqr(mid_p[j],left_p[j],left_rho[j],min_p_jump,
	                       left_c2[j],left_Gam[j],Gam_starl[j],pinfl[j],
			       sqr_min_mass_flux);
	}


	for (j = start; j < end; ++j)
	{
	    Wl[j] = sqrt(Wsql[j]);
	    Wr[j] = sqrt(Wsqr[j]);
	}
}	        /*end find_mass_fluxes*/


LOCAL    double    Mfluxsqr(
	double    mp,
	double    p0,
	double    rho0,
	double    min_dp,
	double    c20,
	double    Gam0,
	double    Gam_star,
	double    pinf,
	double    sqr_min_mass_flux)
{
	double    dp;
	double    den;
	double    M2;

	dp = mp - p0;
	p0 += pinf;
	mp += pinf;

	den = dp + p0*(1.0 - Gam_star/Gam0);
	M2 = (fabs(den) <=  fabs(p0*min_dp)) ?
	     c20*sqr(rho0) : rho0*dp*(mp + Gam_star*0.5*(mp + p0))/den;

	return (M2 > sqr_min_mass_flux) ? M2 : sqr_min_mass_flux;
}	        /*end Mfluxsqr*/


/*
*			find_midstates():
*
*	By now we know the midstate pressures and velocities.  Now we need
*	to find out which state is along the cell boundary.  It could be
*	the midstate, the left state, the right state or we could be in
*	the middle of a rarefaction.  The first thing to do is use the
*	middle velocity to figure out which of the states S, right or left,
*	is closer to the cell boundary.  If u>0, the left is closer, otherwise
*	the right state is closer.  Then, based on whether p*>pS, i.e., whether
*	the separating wave is a shock or a rarefaction, we calculate the
*	wave speed.  If the wave is a shock everything is straightforward--
*	for positive wave speed we know one of the midstates is along the
*	cell boundary, otherwise it is state S.  If we decide that the
*	separating wave is a rarefaction things are more difficult.  Then
*	we find the speed of the midstate side, u* +/- c* and the speed
*	of the known state side of the wave, uS +/- cS.  If the cell
*	boundary is to one side, all is well and we proceed.  If the
*	boundary is in the middle of the rarefaction we use linear
*	interpolation to calculate the state at the boundary.
*/



LOCAL void find_midstates(
	CG_data        *rdata)
{
	Vec_Gas      *vlst = rdata->vlst;
	Vec_Gas      *vrst = rdata->vrst;
	Vec_Gas      *vmst = rdata->vmst;
	int          start = rdata->start, end = rdata->end;
	int          dim = rdata->dim;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double        *min_pressure = rdata->vmuscl->vst->min_pressure +
	                             rdata->vmuscl->offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	double        **left_v = vlst->v;
	double        *left_vn = left_v[0];
	double        *left_p = vlst->p;
	double        *left_rho = vlst->rho;
	double        *left_e = vlst->e;
	double        *left_Gam = rdata->left_Gam;
	double        *left_c2 = vlst->c2;
	double        *pinfl = rdata->pinfl;
	double        *Wl = rdata->Wl;

	double        **right_v = vrst->v;
	double        *right_vn = right_v[0];
	double        *right_p = vrst->p;
	double        *right_rho = vrst->rho;
	double        *right_Gam = rdata->right_Gam;
	double        *right_c2 = vrst->c2;
	double        *right_e = vrst->e;
	double        *pinfr = rdata->pinfr;
	double        *Wr = rdata->Wr;

	double        *mid_vn = vmst->v[0];
	double        *mid_p = vmst->p;
	double        *mid_rho = vmst->rho;
	double        *mid_Gam = rdata->mid_Gam;
	double        *mid_e = vmst->e;
	double        **mid_v = vmst->v;

	double        S_vn, S_vt[MAXD], S_W, S_p;
	double        S_rho, S_Gam, S_c2, S_e, S_vhat, S_pinf;
	double        Gam_bar, gam_bar_i, dp, pbar;
	double        mid_ws, S_ws;
	double        mid_vhat;
	double        sw_tol = rdata->Cg_params.sw_tol;
	double        min_p_jump = rdata->Cg_params.min_p_jump;

	int             i, j;
	static Locstate S = NULL;

#if defined(DEBUG_MUSCL)
	debug_print("cg_rsolve","Entered find_midstates()\n");
#endif /* defined(DEBUG_MUSCL) */

	/*    Allocate space for a Locstate if necessary    */

	if (S == NULL)
	    g_alloc_state(&S,rdata->sizest);

#define strong_wave(pm,pl,pr,rt)                    			\
	((fabs((pl)-(pm)) > (rt)*(pl)) || (fabs((pr)-(pm)) > (rt)*(pr)))

	for (j = start; j < end; ++j)
	{
	    if (
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	        (mid_p[j] < 2.0*min_pressure[j])/*TOLERANCE*/
	        ||
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        (strong_wave(mid_p[j],left_p[j],right_p[j],sw_tol))
	    )
	    {
	        correct_for_strong_waves(j,rdata);
	        continue;
	    }

	    /*
	     *  At this point I want to check for large relative pressure
	     *  jumps between the right and left states.  If a jump is too
	     *  large I'll send the states to the exact Riemann solver and
	     *  then continue this loop at the next cell.
	     *  To be implemented later.
	     */


	    /*
	     * Figure out which side of the middle contact wave is closest
	     * to the cell boundary and call the known state next to it S.
	     */

	    if (mid_vn[j] > 0.0)
	    {
	        S_vn = left_vn[j];
	        for (i = 1; i < dim; ++i)
	            S_vt[i] = left_v[i][j];
	        S_W = Wl[j];
	        S_c2 = left_c2[j];
	        S_p = left_p[j];
	        S_rho = left_rho[j];
	        S_Gam = left_Gam[j];
	        S_e = left_e[j];
	        S_vhat = -S_vn;
	        S_pinf = pinfl[j];
	        Set_params(S,vlst->state[j]);
	        mid_vhat = -mid_vn[j];
	    }
	    else
	    {
	        S_vn = right_vn[j];
	        for (i = 1; i < dim; ++i)
	            S_vt[i] = right_v[i][j];
	        S_W = Wr[j];
	        S_c2 = right_c2[j];
	        S_p = right_p[j];
	        S_rho = right_rho[j];
	        S_Gam = right_Gam[j];
	        S_e = right_e[j];
	        S_vhat = S_vn;
	        S_pinf = pinfr[j];
	        Set_params(S,vrst->state[j]);
	        mid_vhat = mid_vn[j];
	    }

	    /* Are the states connected by a rarefaction or a shock? Once
	     * we've answered that we can calculate the wave speed (or
	     * speeds for a rarefaction) and figure out what things are
	     * like along the cell boundary  */

	    dp = mid_p[j] - S_p;
	    if (fabs(dp) < S_p*min_p_jump)
	    {
	        mid_p[j] = S_p;
	        mid_rho[j] = S_rho;
	        mid_vn[j] = S_vn;
	        mid_Gam[j] = S_Gam;
	    }
	    else if (dp > 0.0)    /* A shock */
	    {

	        /* If the shock speed is negative, then the state along
	         * the boundary is S.  If it is positive, then the
	         * midstate is along the boundary and we don't change
	         * anything */
	        if ((S_vhat + S_W/S_rho) < 0.0)    
	        {    /* Shock speed negative, S along boundary */
	            mid_p[j] = S_p;
	            mid_rho[j] = S_rho;
	            mid_vn[j] = S_vn;
	            mid_Gam[j] = S_Gam;
	        }
	        else
	        {
	            mid_rho[j] = S_rho/(1.0  - dp*S_rho/sqr(S_W));
	            pbar = 0.5*(mid_p[j] + S_p);
	            Gam_bar = 0.5*(right_Gam[j] + left_Gam[j]);
	            gam_bar_i = 2.0*left_p[j]*right_p[j]/
			            (right_c2[j]*right_rho[j]* left_p[j] +
			              left_c2[j]* left_rho[j]*right_p[j]);
	            mid_Gam[j] = S_Gam + DGam_star(Gam_bar,gam_bar_i,
				                   dp,pbar,S_pinf);
	        }

	    }

	    else
	    {
	        /* A rarefaction */

	        S_ws = S_vhat + sqrt(S_c2);
	        if (S_ws < 0.0)
	        {
	            mid_p[j] = S_p;
	            mid_rho[j] = S_rho;
	            mid_vn[j] = S_vn;
	            mid_Gam[j] = S_Gam;
	        }
	        else 
	        {
	            mid_ws = mid_vhat + find_mid_c_rho_Gam(mid_p[j],mid_rho+j,
							   mid_Gam+j,S_p,
							   S_Gam,S_pinf,
							   S,S_rho,S_e,
	                                                   left_Gam[j],
							   right_Gam[j],
							   left_c2[j],
							   left_p[j],
							   left_rho[j],
	                				   right_c2[j],
							   right_p[j],
							   right_rho[j],
							   rdata);
	            if (mid_ws < 0.0)
	            {
	                cell_bdy_in_rarefaction(mid_ws,S_ws,S_p,S_rho,S_vn,
						S_Gam,mid_Gam+j,mid_p+j,
	                                        mid_rho+j,mid_vn+j);
	            }
	        }
	    }
	    mid_e[j] = (mid_p[j] + (mid_Gam[j] + 1.0)*S_pinf)/
		       (mid_Gam[j]*mid_rho[j]);

	    for (i = 1; i < dim; ++i)
	        mid_v[i][j] = S_vt[i];    /* Find tangential velocities */
	}
	Vec_Gas_field_set(vmst,rho) = YES;
	Vec_Gas_field_set(vmst,e) = YES;
	Vec_Gas_field_set(vmst,v) = YES;

#undef strong_wave

#if defined(DEBUG_MUSCL)
	if (debugging("cg_rsolve"))
	{
	    g_printout_vec_data("States from find_midstates()",
	                        mid_rho,mid_e,mid_v,dim,start,end,"muncons");
	    (void) printf("\nMid state pressures and Gammas\n");
	    (void) printf("%-4s %-14s %-14s\n","n","pressure","Gam");
	    for (j = start; j < end; ++j)
	        (void) printf("%-4d %-14g %-14g\n",j,mid_p[j],mid_Gam[j]);
	    (void) printf("\n");
	}
	debug_print("cg_rsolve","Left find_midstates()\n");
#endif /* defined(DEBUG_MUSCL) */
}	        /*end find_midstates*/

LOCAL void correct_for_strong_waves(
	int     j,
	CG_data *rdata)
{
	int             i, dim;
	Vec_Muscl       *vmuscl;
	Locstate        *state;
	Vec_Gas         *vrst = rdata->vrst, *vlst = rdata->vlst;
	double           *pm = rdata->vmst->p;
	double           *rm = rdata->vmst->rho;
	double           *em = rdata->vmst->e;
	double           **vm = rdata->vmst->v;
	static double    dir[3] = {1.0, 0.0, 0.0};
	static Locstate sl = NULL, sr = NULL, smid = NULL;

#if defined(DEBUG_MUSCL)
	debug_print("cg_rsolve","Entered correct_for_strong_waves(), j = %d\n",j);
#endif /* defined(DEBUG_MUSCL) */
	if (smid == NULL)
	{
	    g_alloc_state(&sl,rdata->sizest);
	    g_alloc_state(&sr,rdata->sizest);
	    g_alloc_state(&smid,rdata->sizest);
	    set_type_of_state(sl,EGAS_STATE);
	    set_type_of_state(sr,EGAS_STATE);
	}

	vmuscl = rdata->vmuscl;
	state = vmuscl->vst->state + vmuscl->offset;
	dim = vmuscl->dim;


	Dens(sl)   = vlst->rho[j];
	Energy(sl) = vlst->e[j];
	Dens(sr)   = vrst->rho[j];
	Energy(sr) = vrst->e[j];
	reset_gamma(sl);
	reset_gamma(sr);
	for (i = 0; i < dim; ++i)
	{
	    Vel(sl)[i] = vlst->v[i][j];
	    Vel(sr)[i] = vrst->v[i][j];
	}
	Set_params(sl,state[j]);
	Set_params(sr,state[j]);
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") &&
		(is_bad_state(sl,YES,"correct_for_strong_waves") ||
		 is_bad_state(sr,YES,"correct_for_strong_waves")))
	{
	    screen("ERROR in correct_for_strong_waves(), "
		   "bad state produced at index %d\n",j);
	    (void) printf("sl - ");
	    fprint_raw_gas_data(stdout,sl,current_interface()->dim);
	    (void) printf("sr - ");
	    fprint_raw_gas_data(stdout,sr,current_interface()->dim);
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
	riemann_solution(0.0,dir,sl,sr,smid,EGAS_STATE);
	pm[j] = pressure(smid);
	rm[j] = Dens(smid);
	em[j] = Energy(smid);
	for (i = 0; i < dim; ++i)
	    vm[i][j] = Vel(smid)[i];
#if defined(DEBUG_MUSCL)
	debug_print("cg_rsolve","Left correct_for_strong_waves()\n");
#endif /* defined(DEBUG_MUSCL) */
}	        /*end correct_for_strong_waves*/



/*
*            		find_cons_flux():
*
*    Calculates flux of conserved quantities based on midstate values.
*    Uses three different versions, one for each possible dimension.
*/


LOCAL void find_cons_flux(
	int        start,
	int        end,
	CG_data    *rdata)
{
	Vec_Muscl *vmuscl = rdata->vmuscl;
	double *u[5];
	double *p = rdata->vmst->p;
	int   i;

	u[vmuscl->index.density] = rdata->vmst->rho;
	u[vmuscl->index.energy] = rdata->vmst->e;
	for (i = 0; i < rdata->dim; ++i)
	    u[vmuscl->index.v[i]] = rdata->vmst->v[i];

	flux_vectors(start,end,u,p,rdata->flux,rdata->vmuscl);
}	        	/*end find_cons_flux*/

/*
*            	find_mid_c_rho_Gam():
*
*    Returns the midstate sound speed based on the midstate pressure and
*    the state across a rarefaction wave.  Computes the midstate density
*    and Gamma and places the results in mid_rho and mid_Gam.
*/

LOCAL double find_mid_c_rho_Gam(
	double        mid_p,
	double        *mid_rho,
	double        *mid_Gam,
	double        S_p,
	double        S_Gam,
	double        S_pinf,
	Locstate     S,
	double        S_rho,
	double        S_e,
	double        left_Gam,
	double        right_Gam,
	double        left_c2,
	double	     left_p,
	double        left_rho,
	double        right_c2,
	double	     right_p,
	double	     right_rho,
	CG_data      *rdata)
{
	double        dp, pbar;
	double        Gam_bar, gam_bar_i;
	double        c2;
	static Locstate    known_state = NULL;


	/*    Allocate space for a Locstate if necessary    */

	if (known_state == NULL)
	    g_alloc_state(&known_state, rdata->sizest);


	/*
	 * Use Gamma evolution equation (27) to get Gam on other side of
	 * rarefaction
	 */

	dp = mid_p - S_p;
	pbar = 0.5*(mid_p + S_p);
	Gam_bar = 0.5*(right_Gam + left_Gam);
	gam_bar_i = 2.0*left_p*right_p/(right_c2*right_rho* left_p +
			                 left_c2* left_rho*right_p);
	*mid_Gam = S_Gam + DGam_star(Gam_bar,gam_bar_i,dp,pbar,S_pinf);
	if ((1.0 + *mid_Gam) < 0.0)
	    *mid_Gam = S_Gam;

	/* Now get the midstate density */

	Dens(known_state) = S_rho;
	Energy(known_state) = S_e;
	reset_gamma(known_state);
	Set_params(known_state,S);
	set_type_of_state(known_state,EGAS_STATE);

	*mid_rho = dens_rarefaction(mid_p,known_state);

	/* Return the sound speed using Gam ~ GAM = c^2*rho/(p + pinf) - 1
	 * This is a good approximation since the evolution equation for Gam
	 * holds across a rarefaction */

	c2 = (*mid_Gam + 1.0)*(mid_p + S_pinf)/(*mid_rho);
	return (c2 > 0.0) ? sqrt(c2) : 0.0;
}	        /*end find_mid_c_rho_Gam*/


LOCAL	double	DGam_star(
	double	Gam_bar,
	double   gam_bar_i,
	double	dp,
	double	pbar,
	double	pinf)
{
	return Gam_bar*dp*(1.0/(pbar+pinf) - gam_bar_i*(Gam_bar+1.0)/pbar);
}		/*end DGam_star*/



/*
*			cell_bdy_in_rarefaction():
*
*    This routine calculates the state along the cell boundary when that
*    boundary is inside a rarefaction fan.  Eventually I want it to give 
*    exact answers in the polytropic and stiffened polytropic cases, but
*    for now I'll just use a linear interpolation of the states on either
*    side of the rarefaction as suggested by Colella and Glaz.
*/


LOCAL void cell_bdy_in_rarefaction(
	double        mid_ws,
	double        S_ws,
	double        S_p,
	double        S_rho,
	double        S_vn,
	double        S_Gam,
	double        *mid_Gam,
	double        *mid_p,
	double        *mid_rho,
	double        *mid_vn)
{

	double        a;

	a = mid_ws/(mid_ws - S_ws);
	*mid_p   += a*(S_p   - *mid_p);
	*mid_rho += a*(S_rho - *mid_rho);
	*mid_vn  += a*(S_vn  - *mid_vn);
	*mid_Gam += a*(S_Gam - *mid_Gam);
}	        /*end cell_bdy_in_rarefaction*/



/*
*            free_rsolve_data():
*
*    If we've entered cg_rsolve with more data than we have allocated
*    storage for we need to free up the old space before we can allocate
*    new space.
*/



LOCAL void free_rsolve_data(
	CG_data        *rdata)
{
	/* First make sure there is space to be freed */
	if (rdata == NULL)
	    return;

	free(rdata->worksp);

	zero_scalar(rdata,sizeof(CG_data));
}	        /*end free_rsolve_data*/

LOCAL CG_PARAMS Cg_params;

EXPORT    void    set_cg_params(
	CG_PARAMS	*cg_params)
{
	Cg_params = *cg_params;
}	        /*end set_cg_params*/

/*
*            alloc_rsolve_data():
*
*    We have either entered g_rsolve for the first time or with more data
*    than at any previous time.  Allocate new storage and point everything
*    to the right places.
*/

LOCAL CG_data *alloc_rsolve_data(
	int        mvsize)
{
	static CG_data    Rdata;

	zero_scalar(&Rdata,sizeof(CG_data));
	Rdata.max_vec_size = mvsize;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	bi_array(&Rdata.worksp,16,mvsize,FLOAT);
#else /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	bi_array(&Rdata.worksp,15,mvsize,FLOAT);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	Rdata.Cg_params = Cg_params;

	Rdata.left_Gam = Rdata.worksp[0];
	Rdata.right_Gam = Rdata.worksp[1];
	Rdata.mid_Gam = Rdata.worksp[2];
	Rdata.Gam_min = Rdata.worksp[3];
	Rdata.Gam_max = Rdata.worksp[4];
	Rdata.pinfl = Rdata.worksp[5];
	Rdata.pinfr = Rdata.worksp[6];
	Rdata.Wl = Rdata.worksp[7];
	Rdata.Wr = Rdata.worksp[8];
	Rdata.pnew = Rdata.worksp[9];
	Rdata.pold = Rdata.worksp[10];
	Rdata.vlnew = Rdata.worksp[11];
	Rdata.vlold = Rdata.worksp[12];
	Rdata.vrnew = Rdata.worksp[13];
	Rdata.vrold = Rdata.worksp[14];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Rdata.min_delta_p = Rdata.worksp[15];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	return &Rdata;
}	        /*end alloc_rsolve_data*/
