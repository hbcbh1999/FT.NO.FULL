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
*				ghypsub.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

	/* LOCAL Function Prototypes */

static	Muscl_Opts MOpts;

EXPORT	Muscl_Opts	*muscl_options(void)
{
	return &MOpts;
}		/*end muscl_options*/

/*ARGSUSED*/
EXPORT void oblique_FD(
	double		ds,
	double		dt,
	Tan_stencil	*tsten,
	Locstate	ans,
	Front		*fr)
{
	COMPONENT      comp = tsten->comp;
	const double    *dir = tsten->dir;
	const double    *crds0;
	Locstate       *states = tsten->states;
	Locstate       st, *state;
	double	       **coords;
	double	       *rho, *en_den;
	double	       *m[MAXD];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double	       *vacuum_dens;
	double	       *min_pressure;
	double	       *min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	RECT_GRID      *gr = fr->rect_grid;
	int	       i, j, dim = gr->dim;
	static int     npts, nrad;
	static Vec_Gas *vst = NULL;
	static Vec_Src *src = NULL;
	static double   **Q = NULL;
	static double   *coords_store;
	int	       bad_state_flag = debugging("bad_state") ? YES : NO;

	if (debugging("MUSCLob"))
	{
	    print_general_vector("point = ",Coords(tsten->p[0]),dim,", ");
	    print_general_vector("dir = ",dir,dim,", ");
	    (void) printf("ds = %g, dt = %g\n",ds,dt);
	}

	if (is_obstacle_state(states[0]))
	{
	    g_obstacle_state(ans,fr->sizest);
	    return;
	}
	if (Q == NULL) 
	{
	    npts = tsten->npts;
	    nrad = npts/2;
	    vst = g_alloc_vgas(vst,npts,dim);
	    src = muscl_alloc_vsrc(src,npts,wave_of_front(fr));
	    bi_array(&Q,3,3,FLOAT);
	    vst->Q = (const double* const*)Q;
	    uni_array(&coords_store,3*npts,FLOAT);
	    for (i = 0; i < npts; ++i)
	    	vst->coords[i] = coords_store + 3*i;
	}
	if (RegionIsFlowSpecified(ans,states[0],Coords(tsten->p[0]),
				  comp,comp,fr))
	    return;

	for (i = 0; i < dim; ++i)
	    m[i] = vst->m[i];
	set_rotation(Q,dir,dim);

	/* load the vector vst: */
	clear_Vec_Gas_set_flags(vst);
	rho = vst->rho;
	state = vst->state;
	coords = vst->coords;
	en_den = vst->en_den;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	vacuum_dens = vst->vacuum_dens;
	min_pressure = vst->min_pressure;
	min_energy = vst->min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	crds0 = Coords(tsten->p[0]);
	for (i = 0; i < npts; ++i)
	{
	    state[i] = st = states[i - nrad];
	    coords[i] = Coords(tsten->p[i - nrad]);
	    rho[i] = Dens(st);
	    en_den[i] = Energy(st);
	    for (j = 0; j < dim; ++j)
	    	m[j][i] = scalar_product(Q[j],Mom(st),dim);
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                /* Include partial density into vst */
                Gas_param *params = Params(st);
                int    num_comps;
                double  **rho0 = vst->rho0;
                if((num_comps = params->n_comps) != 1)
                {
                    for(j = 0; j < num_comps; j++)
                        rho0[j][i] = pdens(st)[j];
                }
            }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    vacuum_dens[i] = Vacuum_dens(st);
	    min_pressure[i] = Min_pressure(st);
	    min_energy[i] = Min_energy(st);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	}
	Vec_Gas_field_set(vst,state) = YES;
	Vec_Gas_field_set(vst,coords) = YES;
	Vec_Gas_field_set(vst,rho) = YES;
	Vec_Gas_field_set(vst,en_den) = YES;
	Vec_Gas_field_set(vst,m) = YES;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(states[0])->n_comps != 1)
                Vec_Gas_field_set(vst,rho0) = YES;
        }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Vec_Gas_field_set(vst,vacuum_dens) = YES;
	Vec_Gas_field_set(vst,min_pressure) = YES;
	Vec_Gas_field_set(vst,min_energy) = YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	set_params_jumps(vst,0,npts);

	oned_tangential_scheme(fr,tsten,0,npts,vst,src,dt,ds,dim);

	Dens(ans) = vst->rho[nrad];
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(states[0]);
            if(params->n_comps != 1)
            {
                for(j = 0; j < params->n_comps; j++)
                    pdens(ans)[j] = vst->rho0[j][nrad];
            }
        }
	Energy(ans) = vst->en_den[nrad];
	for (i = 0; i < dim; ++i)
	{
	    Mom(ans)[i] = 0.0;
	    for (j = 0; j < dim; ++j)
	    	Mom(ans)[i] += Q[j][i]*vst->m[j][nrad];
	}
	Set_params(ans,vst->state[nrad]);
	set_type_of_state(ans,GAS_STATE);
	reset_gamma(ans);

	if (is_bad_state(ans,bad_state_flag,"oblique_FD"))
	    LFoblique(ds,dt,tsten,ans,fr);

#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") && is_bad_state(ans,YES,"oblique_FD"))
	{
	    boolean first = YES;

	    if (first == YES)
	    {
	        int  i, nrad = tsten->npts/2;
	        char s[80];

		first = NO;
	        screen("ERROR in oblique_FD(), bad state produced\n");
		(void) printf("ans - ");
	        fprint_raw_gas_data(stdout,ans,current_interface()->dim);
	        print_general_vector("point = ",Coords(tsten->p[0]),dim,", ");
	        print_general_vector("dir = ",dir,dim,", ");
	        (void) printf("ds = %g, dt = %g\n",ds,dt);
	        print_Tan_stencil(fr,tsten);

	        for (i = -nrad; i <= nrad; ++i)
	        {
	            (void) sprintf(s,"tsten->leftst[%d]\n",i);
	            verbose_print_state(s,tsten->leftst[i]);
	            (void) sprintf(s,"tsten->rightst[%d]\n",i);
	            verbose_print_state(s,tsten->rightst[i]);
	        }
	        verbose_print_state("ans",ans);
	        (void) printf("comp = %d\n",comp);
	        (void) printf("Input hypersurface\n");
	        print_hypersurface(tsten->newhs);
		if (tsten->newhs != NULL)
	            print_interface(tsten->newhs->interface);
		add_to_debug("MUSCLob");
		add_to_debug("MUSCL");
		add_to_debug("MUSCL_SOLVER");
		add_to_debug("mtsten_all");
		add_to_debug("oned_MUSCL");
		/*oblique_FD(ds,dt,tsten,ans,fr); */
	        clean_up(ERROR);
	    }
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
	debug_print("MUSCLob","Left oblique_FD():\n");
}			/*end oblique_FD*/


EXPORT	Vec_Muscl *load_state_data(
	int		swp_num,
	int		*iperm,
	Front		*fr,
	Wave		*wave,
	Stencil         *sten,
	Tan_stencil	*tsten,
	Vec_Gas		*vst,
	Vec_Src		*src,
	int		offset,
	int		vsize,
	int		dim,
	double		dt,
	double		dn)
{
	return (*MOpts._load_state_data)(&MOpts,swp_num,iperm,fr,wave,sten,
	                                 tsten,vst,src,offset,vsize,dim,dt,dn);
}		/*end load_state_data*/


/*ARGSUSED*/
EXPORT	void g_load_muscl_flux(
	int	   start,
	int	   end,
	double      **uM,
	Vec_Gas	   *vmst,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	Locstate *state;
	double	 *pM;
	int	 nprms, *prms_jmp;
	int	 i, j, k;
	int	 dim;

	debug_print("mflux","Entered g_load_muscl_flux()\n");

	if (Flux == NULL)
	{
	    debug_print("mflux","Left g_load_muscl_flux()\n");
	    return;
	}

	if (vmuscl->avisc)
	    add_art_visc1(start,end,uM,vmuscl);

	dim = vmuscl->dim;
	pM = vmst->p;
	load_pressure(vmst,start,end-start);

	nprms = vmst->nprms;
	prms_jmp = vmst->prms_jmp;
	state = vmuscl->vst->state + vmuscl->offset;
	for (k = 0; k < nprms; ++k)
	{
	    if (is_obstacle_state(state[prms_jmp[k]]))
	    {
	    	for (j = prms_jmp[k]; j < prms_jmp[k+1]; ++j)
	    	{
	    	    for (i = 0; i < dim; ++i)
			uM[i+2][j] = 0.0;
	    	}
	    }
	}

	/* compute the flux  */
	flux_vectors(start,end,uM,pM,Flux,vmuscl);
	debug_print("mflux","Left g_load_muscl_flux()\n");
}		/*end g_load_muscl_flux*/

EXPORT	void copy_vec_state_params(
	Vec_Gas		*nvst,
	Vec_Gas		*ovst,
	int		offset,
	int		vsize)
{
	int		nprms, *oprms_jmp, *nprms_jmp;
	int		j;
	int		end = vsize+offset;

	nvst->state = ovst->state + offset;
	Vec_Gas_field_set(nvst,state) = YES;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	nvst->vacuum_dens = ovst->vacuum_dens + offset;
	Vec_Gas_field_set(nvst,vacuum_dens) = YES;
	nvst->min_pressure = ovst->min_pressure + offset;
	Vec_Gas_field_set(nvst,min_pressure) = YES;
	nvst->min_energy = ovst->min_energy + offset;
	Vec_Gas_field_set(nvst,min_energy) = YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	if (!Vec_Gas_field_set(ovst,prms_jmp))
	{
	    screen("ERROR in copy_vec_state_params(), "
	           "prms_jmp not set for ovst\n");
	    clean_up(ERROR);
	}
	oprms_jmp = ovst->prms_jmp;
	nprms_jmp = nvst->prms_jmp;
	for (j = 0; offset >= oprms_jmp[j] && j < end; ++j);
	nprms_jmp[0] = 0;
	for (nprms = 1; oprms_jmp[j] < end; ++j)
	{
	    nprms_jmp[nprms++] = oprms_jmp[j] - offset;
	}
	nprms_jmp[nprms] = vsize;
	nvst->nprms = nprms;
	Vec_Gas_field_set(nvst,prms_jmp) = YES;
}		/*end copy_vec_state_params*/


EXPORT	Vec_Gas* g_alloc_vgas(
	Vec_Gas* vst,
	int	 vctr_size,
	int	 dim)
{
	if (vst == NULL)
	{
	    scalar(&vst,sizeof(Vec_Gas));
	    vst->alloc.vst = YES;
	}

	/*MATRIX(vst,m,dim,vctr_size,FLOAT);  Changed on 2.10.2008 to    */
	/*MATRIX(vst,v,dim,vctr_size,FLOAT);  accomodate for TVD solver  */
	MATRIX(vst,m,MAXD,vctr_size,FLOAT);
	MATRIX(vst,v,MAXD,vctr_size,FLOAT);
	VECTOR(vst,rho,vctr_size,FLOAT);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /** According to g_compute_sizest() **/
            MATRIX(vst,rho0,g_nfloats()-(dim+2),vctr_size,FLOAT);
        }
	VECTOR(vst,en_den,vctr_size,FLOAT);
	VECTOR(vst,e,vctr_size,FLOAT);
	VECTOR(vst,re,vctr_size,FLOAT);
	VECTOR(vst,p,vctr_size,FLOAT);
	VECTOR(vst,c,vctr_size,FLOAT);
	VECTOR(vst,GAM,vctr_size,FLOAT);
	VECTOR(vst,c2,vctr_size,FLOAT);
	VECTOR(vst,FD,vctr_size,FLOAT);
	VECTOR(vst,prms_jmp,vctr_size+1,INT);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	VECTOR(vst,vacuum_dens,vctr_size,FLOAT);
	VECTOR(vst,min_pressure,vctr_size,FLOAT);
	VECTOR(vst,min_energy,vctr_size,FLOAT);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	VECTOR(vst,state,vctr_size,sizeof(Locstate));
	VECTOR(vst,coords,vctr_size,sizeof(double*));
	return vst;
}		/*end g_alloc_vgas*/


EXPORT	void	g_free_vgas(
	Vec_Gas* vst)
{
	if (vst == NULL)
	    return;

	Set_free(vst,m);
	Set_free(vst,v);
	Set_free(vst,rho);
	Set_free(vst,en_den);
	Set_free(vst,e);
	Set_free(vst,re);
	Set_free(vst,p);
	Set_free(vst,c);
	Set_free(vst,GAM);
	Set_free(vst,c2);
	Set_free(vst,FD);
	Set_free(vst,prms_jmp);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Set_free(vst,vacuum_dens);
	Set_free(vst,min_pressure);
	Set_free(vst,min_energy);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	Set_free(vst,state);
	Set_free(vst,coords);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /** According to g_compute_sizest() **/
            Set_free(vst,rho0);
        }
	if (vst->alloc.vst)
	    free(vst);
}		/*end g_free_vgas*/

EXPORT	void	g_free_vsrc(
	Vec_Src *src)
{
	if (src == NULL)
	    return;
	Set_free(src,mass);
	Set_free(src,energy);
	Set_free(src,mom);
	Set_free(src,radii);
	Set_free(src,rho0);
	if (src->alloc.src)
	    free(src);
}		/*end g_muscl_free_vsrc*/ 

EXPORT Vec_Muscl *g_muscl_free_wk_space(
	Vec_Muscl	*vmuscl)
{
	Vec_Eigen	*vegn;

	if (vmuscl == NULL)
	    return NULL;

	Set_free(vmuscl,u);
	Set_free(vmuscl,ucon);
	Set_free(vmuscl,source);
	Set_free(vmuscl,central);
	Set_free(vmuscl,forward);
	Set_free(vmuscl,backward);
	Set_free(vmuscl,dulim);
	Set_free(vmuscl,duf);
	Set_free(vmuscl,sgn);
	Set_free(vmuscl,du);
	Set_free(vmuscl,q);
	Set_free(vmuscl,dq);
	Set_free(vmuscl,uL);
	Set_free(vmuscl,pL);
	Set_free(vmuscl,uR);
	Set_free(vmuscl,pR);
	Set_free(vmuscl,right_wvs);
	Set_free(vmuscl,left_wvs);
	Set_free(vmuscl,source_sum);
	Set_free(vmuscl,uM);
	Set_free(vmuscl,pM);
	Set_free(vmuscl,Flux.F);
	Set_free(vmuscl,Flux.H);
	Set_free(vmuscl,uM0);
	Set_free(vmuscl,uM1);
	Set_free(vmuscl,Flux0.F);
	Set_free(vmuscl,Flux0.H);
	Set_free(vmuscl,Flux1.F);
	Set_free(vmuscl,Flux1.H);
	Set_free(vmuscl,awv);
	Set_free(vmuscl,grav);
	Set_free(vmuscl,worksp);

	vegn = &vmuscl->Vegn;
	Set_free(vegn,lambda);
	Set_free(vegn,sgnl);
	Set_free(vegn,sgnr);
	Set_free(vegn,sgnm);
	Set_free(vegn,l);
	Set_free(vegn,r);

	g_free_vgas(&vmuscl->VLst);
	g_free_vgas(&vmuscl->VMst);
	g_free_vgas(&vmuscl->VRst);
	g_free_vgas(&vmuscl->VM0st);
	g_free_vgas(&vmuscl->VM1st);

	if (vmuscl->avisc != NULL)
	{
	    Set_free(vmuscl->avisc,g);
	    Set_free(vmuscl->avisc,cs_ave);
	    Set_free(vmuscl->avisc,c_ave);
	    Set_free(vmuscl->avisc,vn_ave);
	    Set_free(vmuscl->avisc,b);
	    Set_free(vmuscl->avisc,visc);
	    Set_free(vmuscl->avisc,uconM);
	    Set_free(vmuscl->avisc,mdlambda);
	    if (vmuscl->avisc->alloc.avisc)
	        free(vmuscl->avisc);
	}
	if (vmuscl->msf != NULL)
	{
	    Set_free(vmuscl->msf,chi);
	    if (vmuscl->msf->alloc.msf)
	        free(vmuscl->msf);
	}
	Set_free(vmuscl,A);
	Set_free(vmuscl,dV);

	/* Erase references to freed storage */
	if (vmuscl->alloc.vmuscl)
	{
	    free(vmuscl);
	    vmuscl = NULL;
	}
	else
	{
	    zero_scalar(vmuscl,sizeof(Vec_Muscl));
	    set_no_alloc(vmuscl,vmuscl);
	}
	return vmuscl;
}		/*end g_muscl_free_wk_space*/

/*
*			compute_slope_limiting_coeffs():
*
*	Calculates limiting factors for the slopes from the linear
*	reconstruction.	 The factors, chi[j], are between zero and one -- 
*	zero for complete flattening at very strong shocks and one for
*	no flattening at very weak shocks or smooth regions.  Each
*	MUSCL slope slope[k][j] is replaced by slope[k][j]*chi[j].
*
*	Two parameters are used in calculating the chis.  eps is the minimum
*	relative pressure jump for a shock.
*	For eps > (p[j+1] - p[j-1])/min(p[j+1],p[j-1]) we have no limiting at
*	all, i.e. chi[j] = 1.
*
*	eta is a sensitivity factor.  If a zone is determined to be within
*	a shock, a small value of eta will result in strong limiting while
*	a large eta will result in weak limiting.  We actually use the
*	value ieta = 1.0/eta.
*
*	See P. Colella, 'A Direct Eulerian MUSCL Scheme for Gas Dynamics',
*	SIAM J. Sci. Stat. Comput., Vol. 6, No. 1 (1985).
*/

EXPORT	void	compute_slope_limiting_coeffs(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Locstate	*state;
	Gas_param	*prms;
	Vec_Gas		*vst;
	double		chibar;
	double		dlambda, lambdaR, lambdaL;
	double		s;
	double		m;
	double		dp, drho, du, dV;
	double		wv_sp;
	double		*p;
	double		*rho;
	double		msvj;
	double		*vn;
	double		Vjp1, Vjm1;
	double           *GAM, gam;
	int		offset;
	int		j, is;
	boolean		use_msf;
	double           Pjp1, Pjm1;
	double		eps;
	double		ieta;
	double		ch_spd_cutoff;
	double		*chi;
	double		*c, *c2;
	double		do_not_limit;
	double		D;
	double           K0;

	if (vmuscl->msf == NULL)
	    return;

	    /* Read arrays from Vec_Muscl structure */
	offset =	vmuscl->offset;
	vst =		vmuscl->vst;
	vn =		vst->v[0] + offset;
	rho =		vst->rho + offset;
	GAM =		vst->GAM + offset;
	p =		vst->p + offset;
	c =		vst->c + offset;
	c2 =		vst->c2 + offset;
	chi =		vmuscl->msf->chi;
	state =		vst->state + offset;

	prms =		NULL;
	eps =		0.0;
	msvj =		0.0;
	ieta =		0.0;
	ch_spd_cutoff =	0.0;
	for (j = start; j < end; ++j)
	{
	    chi[j] = 1.0;
	    if (prms != Params(state[j]))
	    {
	    	prms = Params(state[j]);
	    	use_msf = (prms != NULL) ?
	    		use_muscl_slope_flattening(prms->avisc) : NO;
	    	K0 = prms->avisc.contact_detector;
                if (use_msf)
	    	{
	    	    eps = prms->avisc.min_shock_jump;
	    	    msvj = prms->avisc.min_sp_vol_jump;
	    	    ieta = prms->avisc.msf_ieta;
	    	    ch_spd_cutoff = prms->avisc.char_speed_cutoff;
	    	}
		else
		{
	            eps =		0.0;
	            msvj =		0.0;
	            ieta =		0.0;
	            ch_spd_cutoff =	0.0;
		}
	    }
	    Pjp1 = c2[j+1]*rho[j+1]/(1.0+GAM[j+1]);
	    Pjm1 = c2[j-1]*rho[j-1]/(1.0+GAM[j-1]);
	    gam = 1.0+0.5*(GAM[j+1]+GAM[j-1]);
	    drho = fabs(rho[j+1]-rho[j-1])/min(rho[j+1],rho[j-1]);
	    dp = fabs(p[j+1] - p[j-1])/min(Pjp1,Pjm1);
	    if (K0*gam*drho < dp)
	    {
		/* Not a contact */
	        if ((p[j+1]-p[j])*(p[j]-p[j-1])<=0.0)/* non-monotone pressure */
		                                     /* use full limiting     */
	        {
		    chi[j] = 0.0;
		}
	        else if (use_msf)
	        {
	            /* No corrections except for sufficiently strong waves */
	            dp = p[j+1] - p[j-1];
		    Pjp1 = c2[j+1]*rho[j+1]/(1.0+GAM[j+1]);
		    Pjm1 = c2[j-1]*rho[j-1]/(1.0+GAM[j-1]);
	            if (fabs(dp) > eps*min(Pjp1,Pjm1))
		    {
	                du = vn[j+1] - vn[j-1];
	                s = dp*du;
	                if (s > 0.0)
	                {
	                    s = 1.0;
		            is = 1;
	                }
	                else if (s < 0.0)
	                {
	                    s = -1.0;
		            is = -1;
	                }
	                else
	                {
	                    s = 0.0;
		            is = 0;
	                }

	                lambdaR = vn[j+1] + s*c[j+1];
	                lambdaL = vn[j-1] + s*c[j-1];
	                do_not_limit = sqr(c[j+1] - c[j-1]);
	                if (lambdaR*lambdaL < ch_spd_cutoff*do_not_limit)
		        {
	                    if (du > 0.0)
	                    {
	    	                /* rarefaction */
	    	                dlambda = fabs(lambdaR - lambdaL);
	    	                lambdaR = fabs(lambdaR);
	    	                lambdaL = fabs(lambdaL);
		                D = lambdaR + lambdaL;
		                chibar = (D != 0.0) ? 1.0 - ieta*(dlambda/D) : 0.0;
	                    }
	                    else if (du < 0.0)
	                    {
	    	                /* transsonic shock */
	    	                Vjp1 = 1.0/rho[j+1];
		                Vjm1 = 1.0/rho[j-1];
	    	                dV = Vjp1 - Vjm1;
	    	                m = (fabs(dV)/min(Vjp1,Vjm1) < msvj) ?
			            rho[j]*c[j] : sqrt( fabs(dp/dV) );
	    	                lambdaR = fabs(lambdaR);
	    	                lambdaL = fabs(lambdaL);

	    	                wv_sp = fabs(m/rho[j+is] + s*vn[j+is]);
		                D = wv_sp + min(lambdaL,lambdaR);
	    	                chibar = (D != 0.0) ?
			            1.0 - ieta*(1.0-wv_sp/D) : 0.0;
	                    }
			    else
			        chibar = 1.0;
	                    chi[j] = max(0.0,chibar);
		        }
		    }
	        }
	    }
	}
	if ((vmuscl->tsten != NULL) &&
	    (vmuscl->tsten->hs != NULL) &&
	    (vmuscl->tsten->hs[0] != NULL) &&
	    (vmuscl->front->rect_grid->dim > 1) &&
	    is_scalar_wave(wave_type(vmuscl->tsten->hs[0])))
	{
	    double curvature = fabs(vmuscl->tsten->curvature);
	    double dh = mag_vector(vmuscl->front->rect_grid->h,
	                          vmuscl->front->rect_grid->dim);
	    double curvature_chi_mod =
	        1.0/(1.0 + curvature_factor(vmuscl->front)*curvature*dh);
	    for (j = start; j < end; ++j)
	        chi[j] *= curvature_chi_mod;
	}

}		/*end compute_slope_limiting_coeffs*/


EXPORT void g_muscl_alloc_phys_vecs(
	Wave		*wave,
	int		vs)
{
	RECT_GRID *gr = wave->rect_grid;
	int       dim = gr->dim;

	g_wave_vgas(wave) = g_alloc_vgas(g_wave_vgas(wave),vs,dim);
	g_wave_vsrc(wave) = muscl_alloc_vsrc(g_wave_vsrc(wave),vs,wave);

}		/*end g_muscl_alloc_phys_vecs*/

EXPORT	Vec_Src	*muscl_alloc_vsrc(
	Vec_Src *src,
	int     vsize,
	Wave    *wave)
{
	RECT_GRID *gr = wave->rect_grid;
	int		dim = gr->dim;
	static boolean	first = YES;
	static double	alpha;

	if (is_rotational_symmetry() && first)
	{
	    first = NO;
	    alpha = rotational_symmetry();
	}
	if (src == NULL && (source_terms_exist() == YES)) 
	{
	    scalar(&src,sizeof(Vec_Src));
	    set_alloc(src,src);
	}
	if (src == NULL)
	    return NULL;

	VECTOR(src,mass,vsize,FLOAT);
	VECTOR(src,energy,vsize,FLOAT);
	MATRIX(src,mom,dim,vsize,FLOAT);
	src->rmin = pos_radius(0.0,gr);
	if (is_rotational_symmetry() && alpha > 0.0)
	    VECTOR(src,radii,vsize,FLOAT);
	else
	    set_no_alloc(src,radii);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /** According to g_compute_sizest() **/
            MATRIX(src,rho0,g_nfloats()-(dim+2),vsize,FLOAT);
        }
        else
	set_no_alloc(src,rho0);
	return src;
}		/*end muscl_alloc_vsrc*/

EXPORT void g_muscl_free_phys_vecs(
	Wave		*wave)
{
	g_free_vgas(g_wave_vgas(wave));
	g_wave_vgas(wave) = NULL;
	g_free_vsrc(g_wave_vsrc(wave));
	g_wave_vsrc(wave) = NULL;
}		/*end g_muscl_free_phys_vecs*/

EXPORT	void g_load_VGas_state_vectors(
	int		offset,
	int		vsize,
	Vec_Gas		*vst,
	int		dim)
{
	int   j;
	double *rho = vst->rho + offset;
	double *en_den = vst->en_den + offset;
	double *re = vst->re + offset;
	double *e = vst->e + offset;
	double *c = vst->c + offset;
	double *c2 = vst->c2 + offset;
	double *mn = vst->m[0] + offset, *vn = vst->v[0] + offset;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double 	*min_energy = vst->min_energy + offset;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	switch (dim)
	{
#if defined(ONED)
	case 1:
	    for (j=0; j < vsize; ++j)
	    {
		vn[j] = mn[j]/rho[j];
		re[j] = en_den[j] - 0.5*mn[j]*vn[j];
		e[j] = re[j]/rho[j];
	    }
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                /* Turn partial density into mass fraction */
                Gas_param *params = Params(vst->state[offset]);
                double     *prho[MAX_NUM_GAS_COMPS];
                int       i, num_comps;
                if((num_comps = params->n_comps) != 1)
                {
                    for(i = 0; i < num_comps; i++)
                    {
                        prho[i] = vst->rho0[i] + offset;
                        for (j=0; j < vsize; ++j)
                            prho[i][j] = prho[i][j]/rho[j];
                    }
                }
            }
	    break;
#endif /* defined(ONED) */

#if defined(TWOD)
	case 2:
	    {
	        double *m2 = vst->m[1] + offset, *v2 = vst->v[1] + offset;
	        for (j=0; j < vsize; ++j)
	        {
		    vn[j] = mn[j]/rho[j];
		    v2[j] = m2[j]/rho[j];
		    re[j] = en_den[j] - 0.5*(mn[j]*vn[j] + m2[j]*v2[j]);
		    e[j] = re[j]/rho[j];
	        }
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    /* Turn partial density into mass fraction */
                    Gas_param *params = Params(vst->state[offset]);
                    double     *prho[MAX_NUM_GAS_COMPS];
                    int       i, num_comps;
                    if((num_comps = params->n_comps) != 1)
                    {
                        for(i = 0; i < num_comps; i++)
                        {
                            prho[i] = vst->rho0[i] + offset;
                            for (j=0; j < vsize; ++j)
                                prho[i][j] = prho[i][j]/rho[j];
                        }
                    }
                }
	    }
	    break;
#endif /* defined(TWOD) */

#if defined(THREED)
	case 3:
	    {
	        double *m2 = vst->m[1] + offset, *v2 = vst->v[1] + offset;
	        double *m3 = vst->m[2] + offset, *v3 = vst->v[2] + offset;
	        for (j=0; j < vsize; ++j)
	        {
		    vn[j] = mn[j]/rho[j];
		    v2[j] = m2[j]/rho[j];
		    v3[j] = m3[j]/rho[j];
		    re[j] = en_den[j] - 
				0.5*(mn[j]*vn[j] + m2[j]*v2[j] + m3[j]*v3[j]);
		    e[j] = re[j]/rho[j];
	        }
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    /* Turn partial density into mass fraction */
                    Gas_param *params = Params(vst->state[offset]);
                    double     *prho[MAX_NUM_GAS_COMPS];
                    int       i, num_comps;
                    if((num_comps = params->n_comps) != 1)
                    {
                        for(i = 0; i < num_comps; i++)
                        {
                            prho[i] = vst->rho0[i] + offset;
                            for (j=0; j < vsize; ++j)
                                prho[i][j] = prho[i][j]/rho[j];
                        }
                    }
                }
	    }
	    break;
#endif /* defined(THREED) */

	default: 
	    screen("ERROR in g_load_VGas_state_vectors(), "
	           "invalid dim = %d\n",dim);
	    clean_up(ERROR);
	    break;
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	for (j = 0; j < vsize; ++j)
	    if (re[j] < min_energy[j])
	    {
		re[j] = min_energy[j];
	        e[j] = re[j]/rho[j];
	    }
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	Vec_Gas_field_set(vst,re) = YES;
	Vec_Gas_field_set(vst,e) = YES;
	Vec_Gas_field_set(vst,v) = YES;

	load_pressure_and_gammas(vst,offset,vsize);

	    /* Sound speed */
	for (j = 0; j < vsize; ++j)
	    c[j] = sqrt(c2[j]);

	Vec_Gas_field_set(vst,c) = YES;

}	     /*end g_load_VGas_state_vectors*/

EXPORT	void	clear_Vec_Gas_set_flags(
	Vec_Gas *vst)
{
	Vec_Gas_field_set(vst,rho) = NO;
	Vec_Gas_field_set(vst,e) = NO; 
	Vec_Gas_field_set(vst,re) = NO; 
	Vec_Gas_field_set(vst,en_den) = NO;
	Vec_Gas_field_set(vst,p) = NO;
	Vec_Gas_field_set(vst,c) = NO;
	Vec_Gas_field_set(vst,GAM) = NO;
	Vec_Gas_field_set(vst,c2) = NO;
	Vec_Gas_field_set(vst,FD) = NO;
	Vec_Gas_field_set(vst,m) = NO;
	Vec_Gas_field_set(vst,v) = NO;
	Vec_Gas_field_set(vst,rho0) = NO;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Vec_Gas_field_set(vst,min_pressure) = NO;
	Vec_Gas_field_set(vst,min_energy) = NO;
	Vec_Gas_field_set(vst,vacuum_dens) = NO;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	Vec_Gas_field_set(vst,prms_jmp) = NO;
	Vec_Gas_field_set(vst,state) = NO;
	Vec_Gas_field_set(vst,coords) = NO;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            Vec_Gas_field_set(vst,rho0) = NO;
}		/*end clear_Vec_Gas_set_flags*/

EXPORT	Vec_Muscl *alloc_Vec_Muscl(
	Vec_Muscl *vmuscl)
{
	if (vmuscl == NULL)
	{
	    scalar(&vmuscl,sizeof(Vec_Muscl));
	    set_alloc(vmuscl,vmuscl);
	}
	else
	{
	    boolean alloc;
	    alloc = vmuscl->alloc.vmuscl;
	    zero_scalar(vmuscl,sizeof(Vec_Muscl));
	    vmuscl->alloc.vmuscl = alloc;
	}
	vmuscl->index.density     = INT_MIN;
	vmuscl->index.energy      = INT_MIN;
	vmuscl->index.v[0]        = INT_MIN;
	vmuscl->index.v[1]        = INT_MIN;
	vmuscl->index.v[2]        = INT_MIN;
	vmuscl->index.pressure    = INT_MIN;
	vmuscl->index.sound_speed = INT_MIN;
	vmuscl->index.GAM         = INT_MIN;
	vmuscl->index.FD          = INT_MIN;
	vmuscl->index.grav        = INT_MIN;
        {
            int i;
            for(i = 0; i < MAX_NUM_GAS_COMPS; i++)
                vmuscl->index.prho[i] = INT_MIN;
        }
	return vmuscl;
}		/*end alloc_Vec_Muscl*/


/*
*			g_assign_wave_state_vectors():
*
*	See comments preceding g_load_state_vectors.  This is the function that
*	loads the Vec_Gas back into the wave structure after the uni_array
*	solver.
*	Also, we need to make sure we don't reset obstacle states (see
*	g_init_obstacle_states() ).
*/

EXPORT	void g_assign_wave_state_vectors(
	int		swp_num,
	int		*iperm,
	Wave		*wv,
	Wave		*newwv,
	Vec_Gas		*vst,
	int		imin,
	int		imax,
	int		*icoords,
	int		pbuf)
{
	int		i, k;
	int		dim = wv->rect_grid->dim;
	int		idirs[MAXD];
	Locstate	state;
	double		speed;

	for (i = 0; i < dim; ++i)
	    idirs[i] = iperm[(i+swp_num)%dim];
	for (i = imin; i < imax; ++i)
	{
	    icoords[idirs[0]] = i + pbuf;
	    if (is_obstacle_state(Rect_state(icoords,wv)))
	        g_obstacle_state(Rect_state(icoords,newwv),wv->sizest);
	    else
	    {
	    	state = Rect_state(icoords,newwv);
	    	Dens(state) = vst->rho[i];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {                     
		    /* Compute differences for mass fraction */ 
                    Gas_param *params = Params(vst->state[i]);
                    int       j;
                    if(params->n_comps != 1)
                    {
                        for(j = 0; j < params->n_comps; j++)
                            pdens(state)[j] = vst->rho0[j][i];
                    }
                }
	    	for (k = 0; k < dim; ++k) 
	    	    Mom(state)[idirs[k]] = vst->m[k][i];
	    	Energy(state) = vst->en_den[i]; 
	    	Set_params(state,vst->state[i]);
	    	set_type_of_state(state,GAS_STATE);
		reset_gamma(state);
		speed = fabs(vel(idirs[swp_num],state)) + sound_speed(state);
		set_max_wave_speed(idirs[swp_num],speed,state,
				Rect_coords(icoords,wv),wv);
	    }
	}
}		/*end g_assign_wave_state_vectors*/

