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
*				gmuscl.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       New and improved version. Contains the drivers for the MUSCL scheme.
*/


#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

	/* LOCAL Function Declarations */
LOCAL	Vec_Muscl *g_load_muscl_state_data(Muscl_Opts*,int,int*,Front*,
                                                 Wave*,Stencil*,
						 Tan_stencil*,
						 Vec_Gas*,Vec_Src*,int,int,
						 int,double,double);
LOCAL	Vec_Muscl *g_muscl_alloc_wk_space(Muscl_Opts*,int,int,Vec_Muscl*);
LOCAL	boolean	g_muscl_load_source_vectors(int,int*,Tan_stencil*,double*,
                                            Vec_Gas*,Vec_Src*,int,int,
					    RECT_GRID*);
LOCAL	double	*g_set_vmuscl_vec_gas_storage(double**,Vec_Gas*,Vec_Muscl*);
LOCAL	void	apply_limiter(double**,int,int,Vec_Muscl*);
LOCAL	void	g_add_art_visc1(int,int,double**,Vec_Muscl*);
LOCAL	void	g_add_art_visc2(int,int,Vec_Muscl*);
LOCAL	void	g_bct_linear_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_compute_art_visc_coefs(int,int,struct _Vec_Muscl*);
LOCAL	void	g_compute_eigens(int,int,struct _Vec_Muscl*);
LOCAL	void	g_compute_lapidus_art_visc_coefs(int,int,Vec_Muscl*);
LOCAL	void	g_cons_src(int,int,int,int*,Tan_stencil*,struct _Vec_Muscl*);
LOCAL	void	g_cons_to_uncons(int,int,int,double**,double**);
LOCAL	void	g_eigenstate_linear_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_first_order_godunov_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_first_order_godunov_half_step(int,int,double,double,Vec_Muscl*);
LOCAL	void	g_half_step(int,int,double,double,Vec_Muscl*);
LOCAL	void	g_lin_rsoln(int,int,double**,double**,double**,Vec_Muscl*);
LOCAL	void	g_muscl_flux_vectors(int,int,double**,double*,MUSCL_FLUX*,
                                     Vec_Muscl*);
LOCAL	void	g_muscl_lin_approx_rsolver(int,int,double**,Vec_Gas*,double**,
                                           Vec_Gas*,double**,Vec_Gas*,
					   MUSCL_FLUX*,Vec_Muscl*);
LOCAL	void	g_state_linear_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_strong_wave_half_step(int,int,double,double,Vec_Muscl*);
LOCAL	void	g_uncons_to_cons(int,int,int,double**,double**);
LOCAL	void	left_multiply_state_by_matrix(double**,double***,double**,
					      int,int,int);
LOCAL	void    print_eigen(int,int,Vec_Eigen*,int);
LOCAL	void    print_linear_reconstruction(int,int,Vec_Muscl*);
LOCAL	void	limit_strong_wave_in_cell(Vec_Muscl*,int,int);


#if defined(DEBUG_MUSCL)
LOCAL	const char **set_muscl_debugging(int,int*,int*,Tan_stencil*,Front*,int);
LOCAL	const char **toggle_muscl_debugging(const char*,int,int*,int*,
					    Front*,int);
LOCAL	void	reset_muscl_debugging(const char**);
#endif /* defined(DEBUG_MUSCL) */


/*
*			oned_MUSCL():
*
*	One dimensional MUSCL code.   The input data are the conservative
*	state variables vst and the source terms src.   The MUSCL is a five
*	point scheme, therefore it updates the states from n1+2 to n2-3.
*	(n2 is the first index that does NOT get processed, thus at a
*	boundary n2-1 is the boundary state).
*/

/*ARGSUSED*/
EXPORT void oned_MUSCL(
	int		swp_num,
	int		*iperm,
	int		*icoords,
	Wave		*wave,
	Wave            *newwave,
	Front		*fr,
	Front           *newfr,
	Stencil         *sten,
	Tan_stencil	*tsten,
	int		offset,
	int		vsize,
	Vec_Gas		*vst,
	Vec_Src		*src,
	double		dt,
	double		dn,
	int		dim)
{
	Vec_Muscl  *vmuscl;
	double	   dtdni = dt/dn;
	double	   **ucon, **F, **source;
	double	   *uconk, *Fk, *sourcek;
	int        start, end;
	int        j, k, kmax;
	int        sten_rad;
/*#define	DEBUG_MUSCL */
#if defined(DEBUG_MUSCL)
	const char	**debug_strings;
#endif /* defined(DEBUG_MUSCL) */

	debug_print("MUSCL","Entered oned_MUSCL()\n");

#if defined(DEBUG_MUSCL)
	debug_strings = set_muscl_debugging(swp_num,iperm,
	        			    icoords,tsten,fr,vsize);
#endif /* defined(DEBUG_MUSCL) */

	vmuscl = load_state_data(swp_num,iperm,fr,wave,sten,tsten,
	        		 vst,src,offset,vsize,dim,dt,dn);
	sten_rad = vmuscl->sten_rad;
	
	/* compute the eigenvalues, eigenuni_arrays and max wavespeed*/
	compute_eigens(0,vsize,vmuscl);

	start = sten_rad-1;
	end = vsize-start;
	/* compute the coefficients for Lapidus and slope limiting viscosity */
	compute_art_visc_coefs(start,end,vmuscl);

	/* compute the linear reconstruction of the state variables */
	reconstructor(start,end,vmuscl);


	/* Evolve for half time step.  The returned data uL and uR
	 * are these evolved states on the two sides of each mesh edge. 
	 * They will be put in the Riemann solver. */
	start = sten_rad;
	end = vsize - sten_rad + 1;
	half_step(start,end,dt,dn,vmuscl);

	/* Solve the Riemann problem at each mesh edge with the states 
	 * uL and uR on the two sides.  The returned data uM is the solution
	 * state at the middle. */

#if defined(DEBUG_MUSCL_TIME)
	start_clock("Riemann_solver");
#endif /* defined(DEBUG_MUSCL_TIME) */

	rsolver(start,end,vmuscl->uL,&vmuscl->VLst,vmuscl->uR,&vmuscl->VRst,
	                  vmuscl->uM,&vmuscl->VMst,&vmuscl->Flux,vmuscl);
	if (vmuscl->avisc)
	    add_art_visc2(start,end,vmuscl);

#if defined(DEBUG_MUSCL)
	if (debugging("lriem"))
	{
	    int i;
	    double **uM = vmuscl->uM;
	    double *pM = vmuscl->pM;
	    g_printout_vec_data("Mid state on the edge, "
	        		"obtained by Riemann solver",
	        	        uM[0],uM[1],uM+2,dim,start,end,
	        		"muncons");
	    (void) printf("Mid state pressure\n"); 
	    (void) printf("%-4s %-14s\n","n","pM"); 
	    for (i = start; i < end; ++i)
	    	(void) printf("%-4d %-14g\n",i,pM[i]);
	    (void) printf("\n");
	}
	if (debugging("vflux")) 
	{
	    double **F = vmuscl->Flux.F;
	    g_printout_vec_data("Here are the flux uni_arrays",F[0],F[1],F+2,
	        		dim,start,end,"cons");
        }
#endif /* defined(DEBUG_MUSCL) */

#if defined(DEBUG_MUSCL_TIME)
	stop_clock("Riemann_solver");
#endif /* defined(DEBUG_MUSCL_TIME) */

	/* Set conservative source terms */
	sten_rad = vmuscl->sten_rad;
	start = sten_rad;
	end = vsize - sten_rad;
	cons_src(start,end,swp_num,iperm,tsten,vmuscl);

	/* compute the cell average of the approximate solution for the
	   next time step */
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
                kmax += num_comps;
        }
	ucon =		vmuscl->ucon;
	F =		vmuscl->Flux.F;
	source =	vmuscl->source;

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Convert mass fraction to partial density */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
            {
                for(k = 0; k < num_comps; k++)
                {
                    uconk = ucon[k+dim+2];
                    for (j = start; j < end; ++j)
                        uconk[j] *= ucon[0][j];
                }
            }
        }

	if (source != NULL)
	{
	    for (k = 0; k < kmax; ++k)
	    {
	    	uconk = ucon[k];
	    	Fk = F[k];
	    	sourcek = source[k];
	    	for (j = start; j < end; ++j)
		{
	    	    uconk[j] += dtdni*(Fk[j] - Fk[j+1]) + dt*sourcek[j];
		}
	    }
	}
	else
	{
	    for (k = 0; k < kmax; ++k)
	    {
	    	uconk = ucon[k];
	    	Fk = F[k];
	    	for (j = start; j < end; ++j)
		{
	    	    uconk[j] += dtdni*(Fk[j] - Fk[j+1]);
		}
	    }
	}

#if defined(DEBUG_MUSCL)
	if (debugging("vgdnv"))
	{
	    g_printout_vec_data("Here are the updated cell averages:",
	        	        ucon[0],ucon[1],ucon+2,dim,start,end,"cons");
	}
	if (debugging("vgdnv_ie"))
	{
	    g_print_internal_energy("Internal energy of updated cell averages:",
	        	            ucon,vmuscl,start,end);
	}
#endif /* defined(DEBUG_MUSCL) */

	debug_print("MUSCL","Left oned_MUSCL()\n");
#if defined(DEBUG_MUSCL)
	reset_muscl_debugging(debug_strings);
#endif /* defined(DEBUG_MUSCL) */
}		/*end oned_MUSCL*/

EXPORT	void	set_muscl_default_opts(
	Muscl_Opts  *mopts,
	Hyp_method  *method,
	int         dim)
{
	MUSCL_PromptType_Reconstructor *Sintrp;
	MUSCL_PromptType_Rsolver       *Rsolver;

	set_muscl_default_hooks(mopts,method,dim);

	uni_array(&Sintrp,5,sizeof(MUSCL_PromptType_Reconstructor));
	mopts->Sintrp = Sintrp;
	Sintrp[0].prompt = "Reconstruct density, energy, velocity";
	Sintrp[0].select = "d";
	Sintrp[0].reconstructor = g_state_linear_reconstructor;
	Sintrp[0].half_step = g_half_step;
	Sintrp[0].strong_wave_half_step = g_strong_wave_half_step;
	Sintrp[1].prompt = "Reconstruct eigen coordinates";
	Sintrp[1].select = "e";
	Sintrp[1].reconstructor = g_eigenstate_linear_reconstructor;
	Sintrp[1].half_step = g_half_step;
	Sintrp[1].strong_wave_half_step = g_strong_wave_half_step;
	Sintrp[2].prompt = "Bell-Colella-Trangenstein reconstruction";
	Sintrp[2].select = "b";
	Sintrp[2].reconstructor = g_bct_linear_reconstructor;
	Sintrp[2].half_step = g_half_step;
	Sintrp[2].strong_wave_half_step = g_strong_wave_half_step;
	Sintrp[3].prompt = "First order Godunov reconstruction (zero slopes)";
	Sintrp[3].select = "f";
	Sintrp[3].reconstructor = g_first_order_godunov_reconstructor;
	Sintrp[3].half_step = g_first_order_godunov_half_step;
	Sintrp[3].strong_wave_half_step = NULL;

	uni_array(&Rsolver,6,sizeof(MUSCL_PromptType_Rsolver));
	mopts->Rsolver = Rsolver;
	Rsolver[0].prompt = "Exact Riemann solver";
	Rsolver[0].select = "e";
	Rsolver[0].rsolver = g_muscl_exact_rsolver;
	Rsolver[0].rmidstate = g_exact_Riemann_midstate;
	Rsolver[1].prompt = "Linear approximate Riemann solver";
	Rsolver[1].select = "l";
	Rsolver[1].rsolver = g_muscl_lin_approx_rsolver;
	Rsolver[1].rmidstate = NULL;/*TODO*/
	Rsolver[2].prompt = "Colella-Glaz's approximate Riemann solver";
	Rsolver[2].select = "c";
	Rsolver[2].rsolver = cg_rsolve;
	Rsolver[2].rmidstate = g_cg_Riemann_midstate;
	Rsolver[3].prompt = "Linear US/UP fit (Dukowicz)";
	Rsolver[3].select = "d";
	Rsolver[3].rsolver = g_linear_us_up_rsolver;
	Rsolver[3].rmidstate = g_linear_us_up_Riemann_midstate;
	Rsolver[4].prompt = "Gamma Law fit";
	Rsolver[4].select = "g";
	Rsolver[4].rsolver = g_gamma_law_fit_rsolver;
	Rsolver[4].rmidstate = g_gamma_law_fit_Riemann_midstate;
}		/*end set_muscl_default_opts*/

EXPORT	void set_muscl_default_hooks(
	Muscl_Opts *mopts,
	Hyp_method  *method,
	int         dim)
{
	mopts->kmax = dim + 2;
	mopts->nfloats = dim+2;
	mopts->_npt_solver = point_FD;
	mopts->_one_side_npt_tang_solver = oblique_FD;
	mopts->_npt_tang_solver = g_two_side_npt_tang_solver;
        mopts->_npt_parab_tan_solver2d = parab_tan_solver2d;
        mopts->_npt_parab_tan_solver3d = parab_tan_solver3d;
	mopts->_print_internal_energy = g_print_internal_energy;
	mopts->_npts_tan_sten = method->npts_sten;
	mopts->_compute_art_visc_coefs = g_compute_art_visc_coefs;
	mopts->_load_state_vectors = g_load_state_vectors;
	mopts->_alloc_phys_vecs = g_muscl_alloc_phys_vecs;
	mopts->_free_phys_vecs = g_muscl_free_phys_vecs;
	mopts->_assign_wave_state_vectors = g_assign_wave_state_vectors;

	mopts->_reconstructor = g_state_linear_reconstructor;
	mopts->_load_state_data = g_load_muscl_state_data;
	mopts->_compute_eigens = g_compute_eigens;
	mopts->_flux_vectors = g_muscl_flux_vectors;
	mopts->_rsolver = g_muscl_exact_rsolver;
	mopts->_rmidstate = g_exact_Riemann_midstate;
	mopts->worksp_len = 5;
	mopts->_half_step = g_half_step;
	mopts->_strong_wave_half_step = g_strong_wave_half_step;
	mopts->_cons_src = g_cons_src;
	mopts->_add_art_visc1 = g_add_art_visc1;
	mopts->_add_art_visc2 = g_add_art_visc2;
	mopts->monotone_reconstruction = NO;
	mopts->link_reconstructions = NO;

	mopts->_Cg_params.pv_iterations = 4;
	mopts->_Cg_params.min_p_jump = 1e-06;          /*TOLERANCE*/
	mopts->_Cg_params.min_v_jump = 1e-06;          /*TOLERANCE*/
	mopts->_Cg_params.sw_tol = 100.0;              /*TOLERANCE*/
	mopts->_Cg_params.min_mass_flux = 1.0e-12;     /*TOLERANCE*/
	mopts->_Cg_params.sqr_min_mass_flux = 1.0e-24; /*TOLERANCE*/
}		/*end set_muscl_default_hooks*/
	

/*ARGSUSED*/
LOCAL	Vec_Muscl *g_load_muscl_state_data(
	Muscl_Opts	*mopts,
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
	static Vec_Muscl *vmuscl = NULL;
	double		 **u, **ucon, **source;
	const double* const *Q = vst->Q;
	boolean		 is_src;
	int		 i;

	if (vmuscl == NULL || vsize > vmuscl->max_vsize || dim > vmuscl->dim)
	{
	    vmuscl = g_muscl_free_wk_space(vmuscl);
	    vmuscl = g_muscl_alloc_wk_space(mopts,vsize,dim,vmuscl);
	}
	vmuscl->sizest = fr->sizest;
	vmuscl->idir = (iperm != NULL) ? iperm[swp_num] : dim;
	vmuscl->Q = Q;
	vmuscl->dt = dt;
	vmuscl->dn = dn;
	vmuscl->front = fr;
	vmuscl->wave = wave;
        vmuscl->sten = sten;
        vmuscl->tsten = tsten;
	if (wave != NULL)
	    vmuscl->sten_rad = wave->npts_sten/2;
	else if (tsten != NULL)
	    vmuscl->sten_rad = tsten->npts/2;
	else
	{
	    vmuscl->sten_rad = 2;
	    screen("ERROR in g_load_muscl_state_data(), can't determine "
	           "stencil size\n");
	    clean_up(ERROR);
	}
	vmuscl->vst = vst;
	vmuscl->src = src;
	vmuscl->offset = offset;	vmuscl->vsize = vsize;

	g_load_VGas_state_vectors(offset,vsize,vst,dim);

	clear_Vec_Gas_set_flags(&vmuscl->VMst);
	clear_Vec_Gas_set_flags(&vmuscl->VLst);
	clear_Vec_Gas_set_flags(&vmuscl->VRst);
	vmuscl->VMst.Q = vmuscl->VLst.Q = vmuscl->VRst.Q = Q;

	/* Set params jump points for vec gas states */
	copy_vec_state_params(&vmuscl->VMst,vst,offset,vsize);
	copy_vec_state_params(&vmuscl->VRst,vst,offset,vsize);
	copy_vec_state_params(&vmuscl->VLst,vst,offset-1,vsize);

	/*put the state variables and source terms in a form more
	  easily usable by the subroutines*/

	u = vmuscl->u;
	ucon = vmuscl->ucon;
	u[vmuscl->index.density] = ucon[0] = vst->rho + offset;
	u[vmuscl->index.energy] = vst->e + offset;
	ucon[1] = vst->en_den + offset;
	for (i = 0; i < dim; ++i)    
	{
	    u[vmuscl->index.v[i]] = vst->v[i] + offset;
	    ucon[i+2] = vst->m[i] + offset;
	}
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Include Mass fraction into vst */
            /* It is converted in g_load_VGas_state_vectors() */
            Gas_param *params = Params(vst->state[offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
            {
                for(i = 0; i < num_comps; i++)
                {
                    ucon[vmuscl->index.prho[i]] = u[vmuscl->index.prho[i]]
                               = vst->rho0[i] + offset;
                }
            }
        }
	if ((src != NULL) && ((source = vmuscl->source) != NULL))
	{
	    source[0] = src->mass + offset;
	    source[1] = src->energy + offset;
	    for (i=0; i < dim; ++i)	 
	    	source[i+2] = src->mom[i] + offset;
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                /* Include Mass fraction into vst */
                /* It is converted in g_load_VGas_state_vectors() */
                Gas_param *params = Params(vst->state[offset]);
                int    num_comps, j;
                if((num_comps = params->n_comps) != 1)
                {
                    for(i = 0; i < num_comps; i++)
                    {
                        source[vmuscl->index.prho[i]] = src->rho0[i] + offset;
                        /* New 051005 */
                        for(j = 0; j < vsize; j++)
                            source[vmuscl->index.prho[i]][j] = 0.0;
                        /* End of New 051005 */
                    }
                }
            }
	}

	if (vmuscl->grav != NULL)
	{
	    double time = fr->time + 0.5*dt;
	    double **coords = vst->coords + offset;
	    double *g = vmuscl->grav + offset;
	    int   j;

	    if (tsten == NULL)
	    {
	        for (j = 0; j < vsize; ++j)
	            g[j] = gravity(coords[j],time)[iperm[swp_num]];
	    }
	    else
	    {
	    	for (j = 0; j < vsize; ++j)
	    	    g[j] = scalar_product(Q[0],gravity(coords[j],time),dim);
	    }
	}

	is_src = g_muscl_load_source_vectors(swp_num,iperm,tsten,vmuscl->grav,
			                     vst,src,offset,vsize,
					     fr->rect_grid);

	vmuscl->is_src = is_src;

#if defined(DEBUG_MUSCL)
	if (debugging("load_Vgas")) 
	{
	    (void) printf("All quantities before MUSCL computation:\n");

	    g_printout_vec_data("The state variables",
				u[vmuscl->index.density],
				u[vmuscl->index.energy],
				u+vmuscl->index.v[0],
				dim,0,vsize,"muncons");
	    g_printout_vec_data("The conserved variables",
				ucon[0],ucon[1],ucon+2,dim,0,vsize,"cons");
	    if (is_src)
	    {
	    	g_printout_vec_data("The source terms",
	    			    source[0],source[1],source+2,dim,1,
	    			    vsize-1,"src");
	    }
 
	    (void) printf("Miscelaneous quantities:\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","pressure",
	    	          "sound speed","Grun. coef");
	    for (i = offset; i < vsize; ++i) 
	        (void) printf("%-4d %-14g %-14g %-14g\n",
			      i,vst->p[i],vst->c[i],vst->GAM[i]);
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_MUSCL) */
	return vmuscl;
}		/*end g_load_muscl_state_data*/


LOCAL boolean g_muscl_load_source_vectors(
	int		   swp_num,
	int		   *iperm,
	Tan_stencil	   *tsten,
	double		   *grav,
	Vec_Gas		   *vst,
	Vec_Src		   *src,
	int		   offset,
	int		   vs,
	RECT_GRID	   *gr)
{
	int		   i, j;
	int		   dim = gr->dim;
	const double* const *Q = vst->Q;
	double		   *rho_src;
	double		   *e_src;
	double		   *v_src[SMAXD];
	double		   *v_src0, *v_srci;
	static	boolean	   first = YES;
	double		*pr;
	double		*v0, *rho;
	static	double	alpha;

	if (src == NULL)
	    return NO;

	rho_src = src->mass + offset;
	e_src = src->energy + offset;

	for (i = 0; i < dim; ++i)
	    v_src[i] = src->mom[i] + offset;
	if (first == YES) 
	{
	    first = NO;
	    if (is_rotational_symmetry())
	    	alpha = rotational_symmetry();
	}

	for (j = 0; j < vs; ++j)
	{
	    rho_src[j] = 0.0;
	    e_src[j] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
	    v_srci = v_src[i];
	    for (j = 0; j < vs; ++j)
		v_srci[j] = 0.0;
	}
	if (grav != NULL)
	{
	    double	*g = grav + offset;
	    v_src0 = v_src[0];
	    for (j = 0; j < vs; ++j)
	    	v_src0[j] += g[j];
	}

	if (tsten == NULL)
	{
	    /* Rect sweep */

	    if (is_rotational_symmetry() && 
		(alpha > 0.0) && (iperm[swp_num]==0))
	    {
		/* Include cylindrical source terms */
 
		double *radii = src->radii + offset;
		double rmin = src->rmin;
		double *m;

		m = vst->m[0] + offset;
		pr = vst->p + offset;
		rho = vst->rho + offset;
		v0 = vst->v[0] + offset;
		for (j = 0; j < vs; ++j)
		{
		    if (fabs(radii[j]) > fabs(rmin))
		    {
		        rho_src[j] -= alpha*m[j]/radii[j];
		        e_src[j] -= alpha*pr[j]*v0[j]/(radii[j]*rho[j]);
		    }
		}
	    }
	}
	else if (is_rotational_symmetry())
	{
	    /* Tangential sweep */

	    /* For cylindrical coordinates the mass equation for the
	    *  tangential sweep becomes
	    *
	    *	drho/dt + d(rho*vtan)/dr = - rho*vtan/r * dr/dtan
	    *
	    *  where dr/dtan is the change in radius with respect to
	    *  changes in the tangential coordinate.  If we let
	    *  (r, z) = tT + nN, i.e. represent r and z by a combination
	    *  of uni_arrays T & N tangent and normal to the curve, dr/dtan
	    *  becomes d(tT[0])/dt = T[0] ( = Q[0][0] )
	    */

	    if (alpha > 0.0)
	    {
	        double rmin, rad;
		double *mom0;
	        POINT **pt;

	        mom0 = vst->m[0];
	        pt = tsten->p;
	        pr = vst->p;
	        rho = vst->rho;
	        v0 = vst->v[0];

		rmin = src->rmin;
	        for (j = 0; j < vs; ++j)
	        {
	            rad = pos_radius(Coords(pt[j-2])[0],gr);
		    if (fabs(rad) > fabs(rmin))
		    {
	                rho_src[j] -= alpha*mom0[j]*Q[0][0]/rad;
	                e_src[j] -= alpha*Q[0][0]*v0[j]*pr[j]/(rho[j]*rad);
		    }
	        }
	    }
	}

#if defined(DEBUG_MUSCL)
	if (debugging("vsrc"))
	{
	    g_printout_vec_data("The source uni_arrays from "
	        	        "g_muscl_load_source_vectors:",
	        	        src->mass,src->energy,src->mom,dim,
	        	        offset,vs,"src");
	}
#endif /* defined(DEBUG_MUSCL) */
	return YES;
}		/*end g_muscl_load_source_vectors*/


/*
*			g_muscl_alloc_wk_space():
*
*	Allocates work space arrays for the uni_arrayized MUSCL finite
*	difference method.
*
*	The memory map is as follows:
*
*	pM = worksp[0][0]
*	backward = worksp[0]
*	uL[i] = worksp[0][i]
*
*	central = worksp[1]
*	uR[i] = worksp[1][i]
*
*	forward = worksp[2]
*	du = worksp[2]
*	right/left_wvs[i] = uM[i] = worksp[2][i]
*
*	source_sum = worksp[3]
*	q = worksp[3]
*	dulim = worksp[3]
*	F[i] = worksp[3][i]
*
*	dq[i] = worksp[4][i]
*/

LOCAL	Vec_Muscl *g_muscl_alloc_wk_space(
	Muscl_Opts *mopts,
	int	   vsize,
	int	   dim,
	Vec_Muscl  *vmuscl)
{
	double		 ***worksp;
	AVISC		 Avisc;
	Vec_Eigen	 *vegn;
	int		 nfloats = mopts->nfloats;
	int		 negn = 3;
	int		 worksp_len = mopts->worksp_len;
	int              i, nvar_u;

	vmuscl = alloc_Vec_Muscl(vmuscl);
	vegn = &vmuscl->Vegn;

	vmuscl->Opts = *mopts;
	vmuscl->dim = dim;

	nvar_u = 0;
	vmuscl->index.density = nvar_u++;
	vmuscl->index.energy = nvar_u++;
	for (i = 0; i < dim; ++i)
	    vmuscl->index.v[i] = nvar_u++;
	vmuscl->nvar_u = nvar_u;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int num_comp;
            num_comp = mopts->kmax - (dim+2);
            for(i = 0; i < num_comp; ++i)
                vmuscl->index.prho[i] = nvar_u++;
            vmuscl->nvar_u = nvar_u;
        }

	zero_scalar(vegn,sizeof(Vec_Eigen));
	vegn->negn = negn;
	set_no_alloc(vegn,vegn);
	MATRIX(vegn,lambda,negn,vsize,FLOAT);
	MATRIX(vegn,sgnl,negn,vsize,FLOAT);
	MATRIX(vegn,sgnr,negn,vsize,FLOAT);
	MATRIX(vegn,sgnm,negn,vsize,FLOAT);
	TRI_ARRAY(vegn,l,negn,negn,vsize,FLOAT);
	TRI_ARRAY(vegn,r,negn,negn,vsize,FLOAT);

	if (is_gravity() == YES)
	    VECTOR(vmuscl,grav,vsize,FLOAT);

	TRI_ARRAY(vmuscl,worksp,worksp_len,nfloats,vsize,FLOAT);
	worksp = vmuscl->worksp;

	VECTOR(vmuscl,u,nfloats,sizeof(double*));
	VECTOR(vmuscl,ucon,nfloats,sizeof(double*));
	if (source_terms_exist() == YES)
	    VECTOR(vmuscl,source,nfloats,sizeof(double*));

	use_artificial_dissipation(&Avisc);
	if (use_lapidus_artificial_viscosity(Avisc))
	{
	    scalar(&vmuscl->avisc,sizeof(Vec_Avisc));
	    set_alloc(vmuscl->avisc,avisc);
	    MATRIX(vmuscl->avisc,g,3,vsize,FLOAT);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,cs_ave,worksp[0][0]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,c_ave,worksp[0][1]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,vn_ave,worksp[0][2]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,b,worksp[1]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,visc,worksp[2][0]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,uconM,worksp[1]);
	    set_no_alloc(vmuscl->avisc,mdlambda);
	}
	else if (use_upwind_artificial_viscosity(Avisc) ||
	         use_linear_artificial_viscosity(Avisc))
	{
	    scalar(&vmuscl->avisc,sizeof(Vec_Avisc));
	    set_alloc(vmuscl->avisc,avisc);

/*          041503, change to old code fashion
	    MATRIX(vmuscl->avisc,g,1,vsize,sizeof(double*));
*/ 
            MATRIX(vmuscl->avisc,g,1,vsize,sizeof(double));
/* end of change 041503 */

	    set_no_alloc(vmuscl->avisc,cs_ave);
	    set_no_alloc(vmuscl->avisc,c_ave);
	    set_no_alloc(vmuscl->avisc,vn_ave);
	    set_no_alloc(vmuscl->avisc,b);
	    set_no_alloc(vmuscl->avisc,visc);
	    set_no_alloc(vmuscl->avisc,uconM);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,mdlambda,worksp[0][0]);
	}

	if (vmuscl->avisc != NULL)
	{
	    vmuscl->avisc->use_lapidus=use_lapidus_artificial_viscosity(Avisc);
	    vmuscl->avisc->use_linear = use_linear_artificial_viscosity(Avisc);
	    vmuscl->avisc->use_upwind = use_upwind_artificial_viscosity(Avisc);
	}

	if (use_muscl_slope_flattening(Avisc))
	{
	    scalar(&vmuscl->msf,sizeof(Vec_MSF));
	    set_alloc(vmuscl->msf,msf);
	    VECTOR(vmuscl->msf,chi,vsize,FLOAT);
	}
	else
	    vmuscl->msf = NULL;
	vmuscl->monotone_reconstruction = mopts->monotone_reconstruction;
	vmuscl->link_reconstructions = mopts->link_reconstructions;

	vmuscl->max_vsize = vsize;

	        /* Set Linear reconstruction data structure */
	ASSIGN_ARRAY_POINTER(vmuscl,backward,worksp[0]);
	ASSIGN_ARRAY_POINTER(vmuscl,central,worksp[1]);
	ASSIGN_ARRAY_POINTER(vmuscl,forward,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,du,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,q,worksp[3]);

	        /* Set half step calculation data */

	ASSIGN_ARRAY_POINTER(vmuscl,uL,worksp[0]);
	ASSIGN_ARRAY_POINTER(vmuscl,uR,worksp[1]);    
	ASSIGN_ARRAY_POINTER(vmuscl,uM,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,right_wvs,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,left_wvs,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,Flux.F,worksp[3]);
	ASSIGN_ARRAY_POINTER(vmuscl,source_sum,worksp[3]);
	ASSIGN_ARRAY_POINTER(vmuscl,awv,worksp[3]);
	ASSIGN_ARRAY_POINTER(vmuscl,dq,worksp[4]);

	        /* Set Riemann solver date */
	vmuscl->pL = g_set_vmuscl_vec_gas_storage(vmuscl->uL,&vmuscl->VLst,
	        	                          vmuscl);
	vmuscl->pR = g_set_vmuscl_vec_gas_storage(vmuscl->uR,&vmuscl->VRst,
	        	                          vmuscl);
	vmuscl->pM = g_set_vmuscl_vec_gas_storage(vmuscl->uM,&vmuscl->VMst,
	        	                          vmuscl);
	set_no_alloc(vmuscl,A);
	set_no_alloc(vmuscl,dV);
	return vmuscl;
}		/*end g_muscl_alloc_wk_space*/


LOCAL	double *g_set_vmuscl_vec_gas_storage(
	double	  **u,
	Vec_Gas	  *vst,
	Vec_Muscl *vmuscl)
{
	int	  dim = vmuscl->dim;
	int	  vsize = vmuscl->max_vsize;
	int i;

	zero_scalar(vst,sizeof(Vec_Gas));
	ASSIGN_ARRAY_POINTER(vst,rho,u[vmuscl->index.density]);
	ASSIGN_ARRAY_POINTER(vst,e,u[vmuscl->index.energy]);
	VECTOR(vst,v,dim,sizeof(double*));
	for (i = 0; i < dim; ++i)
	    vst->v[i] = u[vmuscl->index.v[i]];
	VECTOR(vst,prms_jmp,vsize+1,INT);
	VECTOR(vst,p,vsize,FLOAT);
	VECTOR(vst,re,vsize,FLOAT);
	VECTOR(vst,GAM,vsize,FLOAT);
	VECTOR(vst,FD,vsize,FLOAT);
	VECTOR(vst,c2,vsize,FLOAT);
	VECTOR(vst,c,vsize,FLOAT);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int num_comp;
            num_comp = vmuscl->Opts.kmax - (dim+2);
            VECTOR(vst,rho0,num_comp,sizeof(double*));
            for(i = 0; i < num_comp; ++i)
                vst->rho0[i] = u[vmuscl->index.prho[i]];
        }
	return vst->p;
}		/*end g_set_vmuscl_vec_gas_storage*/


/*
*			g_muscl_flux_vectors();
*
*	This function computes the conservative fluxes for the system.
*/

LOCAL void g_muscl_flux_vectors(
	int	   start,
	int	   end,
	double	   **u,
	double	   *p,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	double *rho = u[vmuscl->index.density];
	double *e = u[vmuscl->index.energy];
	double *v0, *v1, *v2;
	double *F0, *F1, *F2, *F3, *F4;
	double **F = Flux->F;
	int   dim = vmuscl->dim;
	int	j;

	debug_print("mflux","Entered g_muscl_flux_vectors()\n");

	switch (dim)
	{
	case 1:
	    v0 = u[vmuscl->index.v[0]];
	    F0  = F[0]; F1 = F[1]; F2 = F[2];
	    for (j = start; j < end; ++j)
	    {
	    	/* mass */
	        F0[j] = rho[j]*v0[j];
	    	/* energy */
	        F1[j] = v0[j]*(rho[j]*(0.5*sqr(v0[j]) + e[j]) + p[j]);
	    	/* momentum */
	        F2[j] = F0[j]*v0[j] + p[j];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    /* Partial mass flux = total_mass_flux * mass_fraction*/
                    Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
                    int       i, num_comps;
                    if((num_comps = params->n_comps) != 1)
                    {
                        for(i = 0; i < num_comps; i++)
                            F[vmuscl->index.prho[i]][j] = F0[j]*u[vmuscl->index.prho[i]][j];
                    }
                }
	    }
	    break;

	case 2:
	    v0 = u[vmuscl->index.v[0]]; v1 = u[vmuscl->index.v[1]];
	    F0  = F[0]; 
	    F1 = F[1]; 
	    F2 = F[2];
	    F3 = F[3];
	    for (j = start-1; j < end; ++j)
	    {
	    	/* mass */
	        F0[j] = rho[j]*v0[j];
	    	/* energy */
	        F1[j] = v0[j]*
	            (rho[j]*(0.5*(sqr(v0[j])+sqr(v1[j]))+e[j])+p[j]);
	    	/* Sweep component of momentum */
	        F2[j] = F0[j]*v0[j] + p[j];
	    	/* Off sweep component of momentum */
	        F3[j] = F0[j]*v1[j];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    /* Partial mass flux = total_mass_flux * mass_fraction*/
                    Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
                    int       i, num_comps;
                    if((num_comps = params->n_comps) != 1)
                    {
                        for(i = 0; i < num_comps; i++)
                            F[vmuscl->index.prho[i]][j] = F0[j]*u[vmuscl->index.prho[i]][j];
                    }
                }
	    }
	    break;

	case 3:
	    v0 = u[vmuscl->index.v[0]];
	    v1 = u[vmuscl->index.v[1]];
	    v2 = u[vmuscl->index.v[2]];
	    F0  = F[0]; 
	    F1 = F[1]; 
	    F2 = F[2];
	    F3 = F[3];
	    F4 = F[4];
	    for (j = start; j < end; ++j)
	    {
	    	/* mass */
	        F0[j] = rho[j]*v0[j];
	    	/* energy */
	        F1[j] = v0[j]*
	    	    (rho[j]*(0.5*(sqr(v0[j])+sqr(v1[j])+sqr(v2[j]))+e[j])+p[j]);
	    	/* Sweep component of momentum */
	        F2[j] = F0[j]*v0[j] + p[j];
	    	/* Off sweep component of momentum */
	        F3[j] = F0[j]*v1[j];
	    	/* Off sweep component of momentum */
	        F4[j] = F0[j]*v2[j];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    /* Partial mass flux = total_mass_flux * mass_fraction*/
                    Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
                    int       i, num_comps;
                    if((num_comps = params->n_comps) != 1)
                    {
                        for(i = 0; i < num_comps; i++)
                            F[vmuscl->index.prho[i]][j] = F0[j]*u[vmuscl->index.prho[i]][j];
                    }
                }
	    }
	    break;

	default:
	    screen("ERROR in g_muscl_flux_vectors(), "
	           "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	    break;
	}
	debug_print("mflux","Left g_muscl_flux_vectors()\n");
}		/*end g_muscl_flux_vectors*/


/*
*			g_cons_src():
*
*	Transforms source terms from nonconservative to conservative form.
*	This function computes the conservative form of the source terms
*	for the Euler equations at the half-step level.  source*dt
*	is a approximation of the integral of the source vector over the
*	space/time cell [x[j]-dx/2,x[j]+dx/2]X[t,t+dt].
*
*	The macros below define the transformation from the nonconservation
*	representation of the Euler equations using density, specific
*	internal energy, and velocity into the conservative form 
*	using density, total energy, and momentum.
*/

LOCAL void g_cons_src(
	int		start,
	int		end,
	int		swp_num,
	int		*iperm,
	Tan_stencil	*tsten,
	Vec_Muscl	*vmuscl)
{
	double		*grav = vmuscl->grav;
	const double* const *Q = vmuscl->Q;
	int		i, j,dim;
	double		*rho_src, *E_src, *m_src[MAXD], *m_srci;
	double		**uM, **source;
	double		*rho;
	double		*v0;
	static boolean	first = YES;
	double		*vi;
	double		**F;
	double		*rho_flux, *E_flux;
	static	double	alpha;
        double           *prho0_src[MAX_NUM_GAS_COMPS], *prho0_flux[MAX_NUM_GAS_COMPS];
        int             n_comps = 1;

	if (((source = vmuscl->source) == NULL) || (vmuscl->is_src == NO)) 
	    return;

	dim = vmuscl->dim;
	uM = vmuscl->uM;
	rho = uM[0];
	v0 = uM[2];
	rho_src = source[0];
	E_src = source[1];
	for (i = 0; i < dim; ++i)
	    m_src[i] = source[i+2];
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* hook up partial density pointers */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if((n_comps = params->n_comps) != 1)
            {
                for(i = 0; i < params->n_comps; i++)
                    prho0_src[i] = source[vmuscl->index.prho[i]];
            }
        }

	if (first == YES)
	{
	    first = NO;
	    if (is_rotational_symmetry())
		alpha = rotational_symmetry();
	}
	for (j = start; j <= end; ++j)
	{
	    rho_src[j] = 0.0;
	    E_src[j] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
	    m_srci = m_src[i];
	    for (j = start; j <= end; ++j)
	    	m_srci[j] = 0.0;
	}
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* hook up partial density pointers */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if((n_comps = params->n_comps) != 1)
            {
                for (j = start; j <= end; ++j)
                {
                    for(i = 0; i < params->n_comps; i++)
                        prho0_src[i][j] = 0.0;
                }
            }
        }
	if (tsten == NULL)
	{
	    /* Rect sweep */

	    if (is_rotational_symmetry() &&
		(alpha > 0.0) && (iperm[swp_num]==0))
	    {
	        double	rmin, *radii;
	    	/* Include cylindrical source terms */

		rmin = vmuscl->src->rmin;
	    	radii = vmuscl->src->radii;
	    	F = vmuscl->Flux.F;
	    	rho_flux = F[0];
	    	E_flux = F[1];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                   n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                        prho0_flux[i] = F[vmuscl->index.prho[i]];
                }
	    	for (j = start; j < end; ++j)
	    	{
		    if (fabs(radii[j]) > fabs(rmin))
		    {
	    	        rho_src[j] -= alpha*0.5*(rho_flux[j]+rho_flux[j+1])/
				      radii[j];
	    	        E_src[j] -= alpha*0.5*(E_flux[j]+E_flux[j+1])/radii[j];
		    }
	        }
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                       n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                    {
                        for(j = start; j < end; ++j)
                        {
                            if (fabs(radii[j]) > fabs(rmin))
                            {
                                prho0_src[i][j] -= alpha*0.5*(prho0_flux[i][j]
					+prho0_flux[i][j+1])/radii[j];
                            }
                        }
                    }
                }
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    m_srci = m_src[i];
	    	    vi = uM[2+i];
	    	    for (j = start; j < end; ++j)
		    {
	    	    	m_srci[j] += rho_src[j]*0.5*(vi[j]+vi[j+1]);
		    }
	    	}
	    }
	    if (grav != NULL)
	    {
	    	m_srci = m_src[0];
	    	for (j = start; j < end; ++j)
	    	{
	    	    m_srci[j] += 0.5*(rho[j] + rho[j+1])*grav[j];
	    	    E_src[j] += 0.5*(rho[j]*v0[j] + rho[j+1]*v0[j+1])*grav[j];
	    	}
	    }
	}
	else
	{
	    /* Tangential sweep */

	    if (is_rotational_symmetry() &&
		alpha > 0.0)
	    {
	    	POINT	**pt;
	        RECT_GRID *gr = vmuscl->front->rect_grid;
	        double	rmin, rad;
	        double	a = 0.5*alpha*Q[0][0];

	    	/* Include cylindrical source terms */
	    	F = vmuscl->Flux.F;
	    	rho_flux = F[0];
	    	E_flux = F[1];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                   n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                        prho0_flux[i] = F[vmuscl->index.prho[i]];
                }
	    	pt = tsten->p;
		rmin = pos_radius(0.0,gr);
	    	for (j = start; j < end; ++j)
	    	{
	            rad = pos_radius(Coords(pt[j-2])[0],gr);
		    if (fabs(rad) > fabs(rmin))
		    {
	                rho_src[j] -= a*(rho_flux[j]+rho_flux[j+1])/rad;
	                E_src[j] -= a*(E_flux[j]+E_flux[j+1])/rad;
		    }
	        }
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                   n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                    {
                        for(j = start; j < end; ++j)
                        {
                            rad = pos_radius(Coords(pt[j-2])[0],gr);
                            if (fabs(rad) > fabs(rmin))
                                prho0_src[i][j] -= a*(prho0_flux[i][j]+prho0_flux[i][j+1])/rad;
                        }
                    }
                }
	        for (i = 0; i < dim; ++i)
	        {
	            m_srci = m_src[i];
	            vi = uM[2+i];
	            for (j = start; j < end; ++j)
                    {
	      	        m_srci[j] += rho_src[j]*0.5*(vi[j]+vi[j+1]);
                    }
	        }
	    }

	    if (grav != NULL)
	    {
	    	m_srci = m_src[0];
	    	for (j = start; j < end; ++j)
	        {
	            m_srci[j] += 0.5*(rho[j] + rho[j+1])*grav[j];
	            E_src[j] += 0.5*(rho[j]*v0[j] + rho[j+1]*v0[j+1])*grav[j];
	        }
	    }
	}
#if defined(DEBUG_MUSCL)
	if (debugging("csrc"))
	{
	    g_printout_vec_data("Here are the conservative source uni_arrays",
	    	                source[0],source[1],source+2,
	    	                vmuscl->dim,start,end,"src");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_cons_src*/


/*
*			g_state_linear_reconstructor():
*
*	This subroutine constructs piecewise linear functions for the
*	cell averages.	The construction is processed along the direction
*	of the right eigenuni_arrays.  A 'limiter' is applied to avoid
*	oscillation.  The slopes of the piecewise linear function in the
*	directions of the eigenuni_arrays are given by:
*
*	slope[m][j] = du[][j] * l[m][][j]
*
*	du[m][j] = min( 2 * |b[m][j]|, .5 * |c[m][j]|, 2 * |f[m][j]|) *
*					sgn(c[m][j])
*			if sgn(b[m][j]) == sgn(c[m][j]) == sgn(f[m][j])
*		   = 0		otherwise.
*
*	where
*	b[m][j] = (u[m][j]   - u[m][j-1])
*	c[m][j] = (u[m][j+1] - u[m][j-1])
*	f[m][j] = (u[m][j+1] - u[m][j]).
*
*	This limiter can be found (in scalar form) in van Leer, J. Comp.
*	Physics 23, p. 289 (1977).
*
*	Once again we are not introducing any notation for the last two
*	eigenvalues or eigenuni_arrays.
*
*	When using g_load_muscl_state_data() to load the state uni_arrays,
*	we have
*           u[vmuscl->index.density] = rho,
*           u[vmuscl->index.energy] = e,
*           u[vmuscl->index.v[i]] = velocity(i), 0 <= i < dim
*/

LOCAL void g_state_linear_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	double		**u, **dq, **backward, **central, **forward, **du;
	double		***l,***r;
	double		*uk, *bk, *ck, *fk;
	int		num_comps = 0;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_state_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	dim =		vmuscl->dim;
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if(params->n_comps != 1)
            {
                num_comps = params->n_comps;
                kmax += num_comps;
            }
        }
	u =		vmuscl->u;
	l =		vmuscl->Vegn.l;
	r =		vmuscl->Vegn.r;
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;
	du =		vmuscl->du;
	dq =		vmuscl->dq;

	for (k = 0; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    uk = u[k];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] =	      (uk[j]	 - uk[j - 1]);
	    	ck[j] = 0.5 * (uk[j + 1] - uk[j - 1]);
	    	fk[j] =	      (uk[j + 1] - uk[j]);
	    }
	}

	apply_limiter(du,start,end,vmuscl);
	
	/* New 051005 change argument "dim" to "dim+num_comps" */
	left_multiply_state_by_matrix(dq,l,du,dim+num_comps,start,end);
/*
	limit_strong_wave_in_cell(vmuscl,start,end);
*/
#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_state_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_state_linear_reconstructor*/


/*
*			g_bct_linear_reconstructor():
*
*
*    This subroutine constructs piecewise linear functions for the
*    cell averages.  The construction is processed along the direction
*    of the right eigenuni_arrays.	 A 'limiter' is applied to avoid
*    oscillation.  The slopes of the piecewise linear function in the
*    directions of the eigenuni_arrays are given by:
*
*	slope[m][j] = min( 2 * |b[m][j]|, .5 * |c[m][j]|, 2 * |f[m][j]|) *
*					sgn(c[m][j])
*			if sgn(b[m][j]) == sgn(c[m][j]) == sgn(f[m][j])
*		   = 0		otherwise.
*
*	where
*	b[m][j] = (u[][j]   - u[][j-1]) dot l[m][][j]
*	c[m][j] = (u[][j+1] - u[][j-1]) dot l[m][][j]
*	f[m][j] = (u[][j+1] - u[][j]  ) dot l[m][][j].
*
*	This limiter can be found (in scalar form) in van Leer, J. Comp.
*	Physics 23, p. 289 (1977).
*
*	Once again we are not introducing any notation for the last two
*	eigenvalues or eigenuni_arrays.
*/

LOCAL void g_bct_linear_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	double		**u, **backward, **central, **forward;
	double		***l;
	double		*bk, *ck, *fk;
	double		*u0, *u1, *u2, *uk;
	double		*lk0, *lk1, *lk2;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_bct_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	compute_slope_limiting_coeffs(start, end, vmuscl);

	dim =		vmuscl->dim;
	kmax =		dim+2;
	u =		vmuscl->u;
	l =		vmuscl->Vegn.l;
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;

	u0 = u[0];	u1 = u[1];	u2 = u[2];
	for (k = 0; k < 3; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] = ((u0[j] - u0[j - 1])	* lk0[j] +
		    	 (u1[j] - u1[j - 1])	* lk1[j] +
	    		 (u2[j] - u2[j - 1])	* lk2[j]);
		ck[j] = 0.5 * (
			(u0[j + 1] - u0[j - 1]) * lk0[j] +
			(u1[j + 1] - u1[j - 1]) * lk1[j] +
			(u2[j + 1] - u2[j - 1]) * lk2[j]);
		fk[j] = (
			(u0[j + 1] - u0[j])	* lk0[j] +
			(u1[j + 1] - u1[j])	* lk1[j] +
			(u2[j + 1] - u2[j])	* lk2[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    uk = u[k];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] =	      (uk[j]	 - uk[j - 1]);
	    	ck[j] = 0.5 * (uk[j + 1] - uk[j - 1]);
	    	fk[j] =	      (uk[j + 1] - uk[j]);
	    }
	}

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
            {
                for (k = vmuscl->index.prho[0]; k < vmuscl->index.prho[0]+num_comps; ++k)
                {
                    bk = backward[k]; ck = central[k]; fk = forward[k];
                    uk = u[k];
                    for (j = start; j < end; ++j)
                    {
                        bk[j] =       (uk[j]     - uk[j - 1]);
                        ck[j] = 0.5 * (uk[j + 1] - uk[j - 1]);
                        fk[j] =       (uk[j + 1] - uk[j]);
                    }
                }
            }
        }
	apply_limiter(vmuscl->dq,start,end,vmuscl);

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_bct_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_bct_linear_reconstructor*/

/*
*			g_eigenstate_linear_reconstructor():
*
*	This subroutine constructs piecewise linear functions for the
*	cell averages.	The construction is processed along the direction
*	of the right eigenuni_arrays.  A 'limiter' is applied to avoid
*	oscillation.  The slopes of the piecewise linear functions in the
*	directions of the eigenuni_arrays are given by:
*
*	slope[m][j] = min(2 * |b[m][j]|, .5 * |c[m][j]|, |f[m][j]|) *
*				sgn(c[m][j])
*			if sgn(b[m][j]) == sgn(c[m][j]) == sgn(f[m][j])
*		= 0	otherwise.
*
*	where
*	b[m][j] = (q[m][j]   - q[m][j-1])
*	c[m][j] = (q[m][j+1] - q[m][j-1]),
*	f[m][j] = (q[m][j+1] - q[m][j]	)
*
*	and q[m][j] = u[][j] dot l[m][][j]
*
*
*	This limiter can be found (in scalar form) in van Leer, J. Comp.
*	Physics 23, p. 289 (1977).
*
*	Once again we are not introducing any notation for the last two
*	eigenvalues or eigenuni_arrays.
*/

LOCAL void g_eigenstate_linear_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	double		**q = vmuscl->q, *qk;
	double		**u, **backward, **central, **forward;
	double		***l;
	double		*bk, *ck, *fk;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_eigenstate_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	dim =		vmuscl->dim;
	kmax =		dim+2;
	u =		vmuscl->u;
	l =		vmuscl->Vegn.l;
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;

	left_multiply_state_by_matrix(q,l,u,dim,start,end);
	for (k = 0; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    qk = q[k];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] =	      (qk[j]	 - qk[j - 1]);
	    	ck[j] = 0.5 * (qk[j + 1] - qk[j - 1]);
	    	fk[j] =	      (qk[j + 1] - qk[j]);
	    }
	}

	apply_limiter(vmuscl->dq,start,end,vmuscl);

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_eigenstate_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_eigenstate_linear_reconstructor*/

LOCAL	void	apply_limiter(
	double		**slope,
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;    /* j is the mesh index */
	int		kmax;
	double		**backward, **central, **forward;
	double		*bk, *ck, *fk;
	double		*slopek;
	double		sign, dl, abk, afk, ack;

	dim =		vmuscl->dim;
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
                kmax += num_comps;
        }
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	{
	    (void) printf("\n");
	    g_printout_vec_data("Backward differences",backward[0],backward[1],
	        		backward+2,dim,start,end,"muncons");
	    g_printout_vec_data("Central differences",central[0],central[1],
	        		central+2,dim,start,end,"muncons");
	    g_printout_vec_data("Forward differences",forward[0],forward[1],
	        		forward+2,dim,start,end,"muncons");
	}
#endif /* defined(DEBUG_MUSCL) */

	compute_slope_limiting_coeffs(start,end,vmuscl);

	/* Compute the standard Van Leer slope limiter */

	for (k = 0; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    slopek = slope[k];
	    for (j = start; j < end; ++j)
	    {
	    	if (bk[j]*fk[j] > 0.0)
	    	{
	    	    abk = fabs(bk[j]);
	    	    afk = fabs(fk[j]);
	    	    dl = 2.0*min(abk,afk);
	    	    sign = (ck[j] == 0.0) ? 0.0 : (ck[j] > 0.0) ? 1.0 : -1.0;
	    	    ack = fabs(ck[j]);
	    	    slopek[j] = sign*min(ack,dl);
	    	}
	    	else
	    	{
	    	    slopek[j] = 0.0;
	    	}
	    }
	}
	if (vmuscl->msf != NULL)
	{
	    double *chi = vmuscl->msf->chi;
	    for (k = 0; k < kmax; ++k)
	    {
	    	slopek = slope[k];
	    	for (j = start; j < end; ++j)
	    	    slopek[j] *= chi[j];
	    }
	}
}		/*end apply_limiter*/
	

LOCAL void g_first_order_godunov_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	double		**dq, *dqk;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_first_order_godunov_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	dim =		vmuscl->dim;
	kmax =		dim+2;
	dq =		vmuscl->dq;

	for (k = 0; k < kmax; ++k)
	{
	    dqk = dq[k];
	    for (j = start; j < end; ++j)
	        dqk[j] = 0.0;
	}

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_first_order_godunov_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_first_order_godunov_reconstructor*/


/*
*			 g_compute_eigens():
*
*	Calculate the eigenvalues and eigenuni_arrays of the linearized gas
*	dynamical equations in the normal/tangential direction.
*	The state vector is u = (rho,e,v0,v1,v2), which are to be
*	interpreted as cell averages and with all
*	uni_arrays written in these coordinates.  v0 is always the normal
*	velocity.
*
*
*	The eigenvalues are:
*
*		lambda[0] = v0 - c
*		lambda[1] = lambda[3] = lambda[4] = v0
*		lambda[2] = v0 + c
*
*	c2 = sound speed squared, GAM = Gruneisen exponent.
*	The right eigenuni_arrays are (column uni_arrays):
*
*	      r[][0]	       r[][1]	          r[][2]   r[][3] r[][4]
*
*	    ( -rho/c	         GAM/c       	   rho/c      0	     0	 )
*	    (-p/(rho*c) (GAM*p/rho-c2)/(rho*c)   p/(rho*c)    0	     0	 )
*	    (	 1	         0	            1	      0	     0	 )
*	    (	 0	         0	            0	      1	     0	 )
*	    (	 0	         0	            0	      0	     1	 )
*
*	The left eigenuni_arrays are (row uni_arrays):
*
*	l[0][] ( (GAM*p/rho-c2)(2*rho*c)  -GAM/(2c)  1/2      0	     0	 )
*	l[1][] (    p/(rho*c)	          -rho/c      0	      0	     0	 )
*	l[2][] (-(GAM*p/rho-c2)/(2*rho*c)  GAM/(2c)  1/2      0	     0	 )
*	l[3][] (       0		      0	      0	      1	     0	 )
*	l[4][] (       0		      0	      0	      0	     1	 )
*
*	Since l3 and l4 (resp. r3 and r4) are perpendicular to l0, l1, l2
*	(resp. r0, r1, r2), we shall not use the third or fourth 
*	components of l0, l1, l2 (resp. r0, r1, r2). We shall also not
*	introduce the notations l3, r3, lambda3, l4, r4, or lambda4.
*
*	Returns the maximum absolute value of the eigenvalues.
*/


LOCAL void g_compute_eigens(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Eigen  *vegn = &vmuscl->Vegn;
	Vec_Gas    *vst = vmuscl->vst;
	int        offset = vmuscl->offset;
	int        j, k;
	double      **lambda = vegn->lambda;
	double      ***l = vegn->l, ***r = vegn->r;
	double      *rho = vst->rho + offset;
	double      *v0 = vst->v[0] + offset;
	double      *c = vst->c + offset;
	double      *p = vst->p + offset;
	double      *GAM = vst->GAM + offset;
	double      *c2 = vst->c2 + offset;
	double      **sgnl = vegn->sgnl;
	double      **sgnr = vegn->sgnr;
	double      speed;
	double      *sgnlk, *sgnrk;
	double      *lambdak;
	double      *lambda0, *lambda1, *lambda2;
	double      *r00, *r01, *r02;
	double      *r10, *r11, *r12;
	double      *r20, *r21, *r22;
	double      *l00, *l01, *l02;
	double      *l10, *l11, *l12;
	double      *l20, *l21, *l22;
	Locstate   *state = vst->state + offset;

	lambda0 = lambda[0]; lambda1 = lambda[1]; lambda2 = lambda[2];
	r00 = r[0][0]; r01 = r[0][1]; r02 = r[0][2];
	r10 = r[1][0]; r11 = r[1][1]; r12 = r[1][2];
	r20 = r[2][0]; r21 = r[2][1]; r22 = r[2][2];

	l00 = l[0][0]; l01 = l[0][1]; l02 = l[0][2];
	l10 = l[1][0]; l11 = l[1][1]; l12 = l[1][2];
	l20 = l[2][0]; l21 = l[2][1]; l22 = l[2][2];
	
	for (j = start; j < end; ++j)
	{
	    lambda0[j] = v0[j] - c[j];
	    lambda1[j] = v0[j];
	    lambda2[j] = v0[j] + c[j];

	    r00[j] = -rho[j] / c[j];
	    r10[j] = -p[j] / (rho[j]*c[j]);
	    r20[j] = 1.0;

	    r01[j] = GAM[j] / c[j];
	    r11[j] = (GAM[j]*p[j]/rho[j]-c2[j]) / (rho[j]*c[j]);
	    r21[j] = 0.0;

	    r02[j] = rho[j] / c[j];
	    r12[j] = p[j] / (rho[j]*c[j]);
	    r22[j] = 1.0;

	    l00[j] = 0.5 * (GAM[j]*p[j]/rho[j]-c2[j]) / (rho[j]*c[j]);
	    l01[j] = -0.5 * GAM[j] / c[j];
	    l02[j] = 0.5;

	    l10[j] = p[j] / (rho[j]*c[j]);
	    l11[j] = -rho[j] / c[j];
	    l12[j] = 0.0;

	    l20[j] = 0.5 * (c2[j]-GAM[j]*p[j]/rho[j]) / (rho[j]*c[j]);
	    l21[j] = 0.5 * GAM[j] / c[j];
	    l22[j] = 0.5;
	}
	for (j = start; j < end; ++j)
	{
	    if (Params(state[j]) != NULL)
	    {
	        speed = Params(state[j])->avisc.sp_coef*(fabs(v0[j]) + c[j]);
	        if (vmuscl->wave != NULL)
	        {
	            set_max_wave_speed(vmuscl->idir,speed,vst->state[j+offset],
	                               vst->coords[j+offset],vmuscl->wave);
	        }
	        else
	        {
	            int    k, dim = vmuscl->dim;

	            for (k = 0; k < dim; ++k)
	            {
	                set_max_front_speed(k,vmuscl->Q[0][k]*speed,
	                                    vst->state[j+offset],
					    vst->coords[j+offset],
					    vmuscl->front);
	            }
	            set_max_front_speed(dim,speed/vmuscl->dn,
	                                vst->state[j+offset],
					vst->coords[j+offset],
					vmuscl->front);
	        }
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    sgnlk = sgnl[k]; sgnrk = sgnr[k]; lambdak = lambda[k];
	    for (j = start; j < end; ++j)
	    {
	        if (lambdak[j] > 0.0)
	        {
	            sgnlk[j] = 0.0;
	            sgnrk[j] = 1.0;
	        }
	        else
	        {
	            sgnlk[j] = 1.0;
	            sgnrk[j] = 0.0;
	        }
	    }
	}

#if defined(DEBUG_MUSCL)
	if (debugging("eigen")) print_eigen(start,end,vegn,1);
#endif /* defined(DEBUG_MUSCL) */

	return;
}        /*end g_compute_eigens*/


/*
*			    g_half_step():
*
*	This routine solves the linearized gas dynamics equations in each
*	mesh with the above piecewise linear data for half time step and
*	before interaction.  The resulting left/right states at the edge
*	of each mesh are:
*
*		   uL[][j] = u[][j-i] + duL[][j]
*		   uR[][j] = u[][j] + duR[][j]
*
*		   duL[][j] = (sum over lambda[k][j-1] >= 0)
*			     { 0.5*(1-lambda[k][j-1]*dt/dx)*dq[k][j-1]
*			     + 0.5*dt*sqrt(1+lambda[k][j-1]**2)}*r[][k][j-1]
*
*		   duR[][j] = (sum over lambda[k][j] < 0)
*			     { 0.5*(-1-lambda[k][j]*dt/dx)*dq[k][j]
*			     + 0.5*dt*dqrt(1+lambda[k][j]**2}*r[][k][j]
*
*	These data will be put in the Riemann solver.
*/
/* NOTE: The g_half_step is modified according to
 * "Multidimensioanl Upwind Methods for Hyperbolic Conservation Laws"
 *  By P. Colella, J.C.P. 87 171-200 (1990).
 *  to introduce a reference state.
 */

LOCAL void g_half_step(
	int		start,
	int		end,
	double		dt,
	double		dn,
	Vec_Muscl	*vmuscl)
{
	Vec_Eigen	*vegn;
	int		j, k;
	int		dim, kmax;
	double		**u, **source, **dq, **uL, **uR;
	double		**right_wvs, **left_wvs, **source_sum;
	double		***l, ***r;
	double		**lambda;
	double		**sgnl, **sgnr;
	double		*ssk;
	double		*rk0, *rk1, *rk2;
	double		*ws0, *ws1, *ws2, *wsk;
	double		*lambdak, *lambda1;
	double		*dqk;
	double		*sgnrk, *sgnr1, *sgnlk, *sgnl1;
	double		*uk, *uLk, *uRk;
	double		dtdni = dt/dn;
	static Locstate tmpst = NULL;
	double           lambdaM, lambdaL;
	int             num_comps = 0;

	dim =		vmuscl->dim;
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if(params->n_comps != 1)
            {
                num_comps = params->n_comps;
                kmax += num_comps;
            }
        }
	u =		vmuscl->u;
	source =	vmuscl->source;
	dq =		vmuscl->dq;
	uL =		vmuscl->uL;
	uR =		vmuscl->uR;
	right_wvs =	vmuscl->right_wvs;
	left_wvs =	vmuscl->left_wvs;
	source_sum =	vmuscl->source_sum;
	vegn =		&vmuscl->Vegn;
	l =		vegn->l;
	r =		vegn->r;
	lambda =	vegn->lambda;
	sgnl =		vegn->sgnl;
	sgnr =		vegn->sgnr;

	/* add up source term contributions for half time-step */

	if ((source != NULL) && (vmuscl->src != NULL))
	{
	    /* New 051005, change argument "source,dim" to "source,dim+num_comps" */
	    left_multiply_state_by_matrix(source_sum,l,source,dim+num_comps,start-1,end);
	    for (k = 0; k < kmax; ++k)
	    {
	    	ssk = source_sum[k];
		dqk = dq[k];
	    	for (j = start - 1; j < end; ++j)
		{
		    if(fabs(dqk[j]) > MACH_EPS)
                        ssk[j] *= 0.5 * dt;
                    else
                        ssk[j] = 0.0;
		}
	    }
	}
	else
	{
	    for (k = 0; k < kmax; ++k)
	    {
	    	ssk = source_sum[k];
	    	for (j = start - 1; j < end; ++j)
	            ssk[j] = 0.0;
	    }
	}

	/* add up contributions from right moving waves */
	for (k = 0; k < 3; ++k)
	{
	    wsk = right_wvs[k];
	    ssk = source_sum[k];
	    lambdak = lambda[k];
	    dqk = dq[k];
	    sgnrk = sgnr[k];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnrk[j-1]*(0.5*(1.0-dtdni*lambdak[j-1])*dqk[j-1] +
	    			ssk[j-1]);
	    }
	}

	for (k = 3; k < kmax; ++k)
	{
	    wsk = right_wvs[k];
	    ssk = source_sum[k];
	    lambda1 = lambda[1];
	    dqk = dq[k];
	    sgnr1 = sgnr[1];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnr1[j-1]*(0.5*(1.0-dtdni*lambda1[j-1])*dqk[j-1] +
					ssk[j-1]);
	    }
	}

	/* compute uL */
	for (k = 0; k < 3; ++k)
	{
	    rk0 = r[k][0];    rk1 = r[k][1];    rk2 = r[k][2];
	    ws0 = right_wvs[0]; ws1 = right_wvs[1]; ws2 = right_wvs[2];
	    uk = u[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	uLk[j] = uk[j-1] +
			 (ws0[j]*rk0[j-1] + ws1[j]*rk1[j-1] + ws2[j]*rk2[j-1]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    wsk = right_wvs[k];
	    uk = u[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	uLk[j] = uk[j-1] + wsk[j];
	    }
	}
	/*Check for negative densities */
	for (j = start; j < end; ++j)
	{
	    Locstate st0 = vmuscl->vst->state[j-1+vmuscl->offset];
	    if (vmuscl->VLst.rho[j] < 0.0)
	    {
		double    min_p = Min_pressure(st0);

		if (tmpst == NULL)
		    alloc_state(vmuscl->front->interf,&tmpst,
		                vmuscl->front->sizest);
		state_on_adiabat_with_pr(st0,min_p,tmpst,EGAS_STATE);
	        vmuscl->VLst.rho[j] = Dens(tmpst);
	        vmuscl->VLst.e[j] = Energy(tmpst);
		vmuscl->VLst.v[0][j] = u[2][j-1]+riemann_wave_curve(st0,min_p);
	    }
            /* If the reconstruction gives bad state, down to 1st order scheme */
            if(vmuscl->VLst.rho[j]*vmuscl->VLst.e[j] < Min_energy(st0))
            {
                vmuscl->VLst.rho[j] = u[0][j-1];
                vmuscl->VLst.e[j] = u[1][j-1];
                vmuscl->VLst.v[0][j] = u[2][j-1];
            }
	}
	Vec_Gas_field_set(&vmuscl->VLst,rho) = YES;
	Vec_Gas_field_set(&vmuscl->VLst,e) = YES;
	Vec_Gas_field_set(&vmuscl->VLst,v) = YES;

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Scaling the mass fraction */
            double       sum_f;
            if(num_comps != 1)
            {
                for (j = start; j < end; ++j)
                {
                    sum_f = 0.0;
                    for(k = 0; k < num_comps; k++)
                        sum_f += vmuscl->VLst.rho0[k][j];
                    for(k = 0; k < num_comps; k++)
                        vmuscl->VLst.rho0[k][j] *= 1.0/sum_f;
                }
                Vec_Gas_field_set(&vmuscl->VLst,rho0) = YES;
            }
        }

	/* add up contributions from left moving waves */
	for (k = 0; k < 3; ++k)
	{
	    wsk = left_wvs[k];
	    ssk = source_sum[k];
	    lambdak = lambda[k];
	    dqk = dq[k];
	    sgnlk = sgnl[k];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnlk[j]*(0.5*(-1.0 - dtdni*lambdak[j])*dqk[j] +
					ssk[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    wsk = left_wvs[k];
	    ssk = source_sum[k];
	    lambda1 = lambda[1];
	    dqk = dq[k];
	    sgnl1 = sgnl[1];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnl1[j]*(0.5*(-1.0 - dtdni*lambda1[j])*dqk[j] +
					ssk[j]);
	    }
	}

	/* compute uR */
	for (k = 0; k < 3; ++k)
	{
	    rk0 = r[k][0];    rk1 = r[k][1];    rk2 = r[k][2];
	    ws0 = left_wvs[0]; ws1 = left_wvs[1]; ws2 = left_wvs[2];
	    uk = u[k];
	    uRk = uR[k];
	    for (j = start; j < end; ++j)
	    {
	    	uRk[j] = uk[j] + (ws0[j]*rk0[j]+ws1[j]*rk1[j]+ws2[j]*rk2[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    wsk = left_wvs[k];
	    uk = u[k];
	    uRk = uR[k];
	    for (j = start; j < end; ++j)
	    {
	    	uRk[j] = uk[j] + wsk[j];
	    }
	}
	/*Check for negative densities */
	for (j = start; j < end; ++j)
	{
	    Locstate st0 = vmuscl->vst->state[j+vmuscl->offset];
	    if (vmuscl->VRst.rho[j] < 0.0)
	    {
		double    min_p = Min_pressure(st0);

		if (tmpst == NULL)
		    alloc_state(vmuscl->front->interf,&tmpst,
		                vmuscl->front->sizest);
		state_on_adiabat_with_pr(st0,min_p,tmpst,EGAS_STATE);
	        vmuscl->VRst.rho[j] = Dens(tmpst);
	        vmuscl->VRst.e[j] = Energy(tmpst);
		vmuscl->VRst.v[0][j] = u[2][j]-riemann_wave_curve(st0,min_p);
	    }
            /* If the reconstruction gives bad state, down to 1st order scheme */
            if(vmuscl->VRst.rho[j]*vmuscl->VRst.e[j] < Min_energy(st0))
            {
                vmuscl->VRst.rho[j] = u[0][j];
                vmuscl->VRst.e[j] = u[1][j];
                vmuscl->VRst.v[0][j] = u[2][j];
            }
	}
	Vec_Gas_field_set(&vmuscl->VRst,rho) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,e) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,v) = YES;

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Scaling the mass fraction */
            double       sum_f;
            if(num_comps != 1)
            {
                for (j = start; j < end; ++j)
                {
                    sum_f = 0.0;
                    for(k = 0; k < num_comps; k++)
                        sum_f += vmuscl->VRst.rho0[k][j];
                    for(k = 0; k < num_comps; k++)
                        vmuscl->VRst.rho0[k][j] *= 1.0/sum_f;
                }
                Vec_Gas_field_set(&vmuscl->VRst,rho0) = YES;
            }
        }
#if defined(DEBUG_MUSCL)
	if (debugging("half_step"))
	{
	    int i, j;

	    g_printout_vec_data("Left state on the edge",
				uL[vmuscl->index.density],
				uL[vmuscl->index.energy],
				uL+vmuscl->index.v[0],dim,start,end,"muncons");
	    g_printout_vec_data("Right state on the edge",
				uR[vmuscl->index.density],
				uR[vmuscl->index.energy],
				uR+vmuscl->index.v[0],dim,start,end,"muncons");

	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","density","int_en");
	    for (i = 0; i < dim; ++i)
		(void) printf(" %-4s%1d%-9s","vel[",i,"]");
	    (void) printf(" Solution before half step\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",u[vmuscl->index.density][j]);
		(void) printf(" %-14g",u[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",u[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");

	    (void) printf(" Solution after half step\n");
	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","l_density","l_int_en");
	    for (i = 0; i < dim; ++i)
		    (void) printf(" %-4s%1d%-9s","l_vel[",i,"]");
	    (void) printf(" %-14s %-14s","r_density","r_int_en");
	    for (i = 0; i < dim; ++i)
		    (void) printf(" %-4s%1d%-9s","r_vel[",i,"]");
	    (void) printf("\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",uL[vmuscl->index.density][j]);
		(void) printf(" %-14g",uL[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
			(void) printf(" %-14g",uL[vmuscl->index.v[k]][j]);
		(void) printf(" %-14g",uR[vmuscl->index.density][j]);
		(void) printf(" %-14g",uR[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
			(void) printf(" %-14g",uR[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_half_step*/


LOCAL void g_strong_wave_half_step(
	int		start,
	int		end,
	double		dt,
	double		dn,
	Vec_Muscl	*vmuscl)
{
	boolean  redo = NO;
	double *rhol, *rhor;
	double *el, *er;
	double **dq, **uL, **uR;
	int   j, k, kmax;

	g_half_step(start,end,dt,dn,vmuscl);

	dq =		vmuscl->dq;
	uL =		vmuscl->uL;
	uR =		vmuscl->uR;
	kmax =		vmuscl->dim+2;

	rhol = uL[vmuscl->index.density];
	el = uL[vmuscl->index.energy];
	rhor = uR[vmuscl->index.density];
	er = uR[vmuscl->index.energy];
	for (j = start; j < end; ++j)
	{
	    if ((rhol[j] < 0.0) || (el[j] < 0.0))
	    {
	    	for (k = 0; k < kmax; ++k)
	    	    dq[k][j-1] = 0.0;
	    	redo = YES;
	    }
	    if ((rhor[j] < 0.0) || (er[j] < 0.0))
	    {
	    	for (k = 0; k < kmax; ++k)
	    	    dq[k][j] = 0.0;
	    	redo = YES;
	    }
	}
	if (redo == YES)
	    g_half_step(start,end,dt,dn,vmuscl);
}		/*end g_strong_wave_half_step*/


/*ARGSUSED*/
LOCAL void g_first_order_godunov_half_step(
	int		start,
	int		end,
	double		dt,
	double		dn,
	Vec_Muscl	*vmuscl)
{
	int		j, k;
	int		dim, kmax;
	double		**u, **uL, **uR;
	double		*uk, *uLk, *uRk;

	dim =		vmuscl->dim;
	kmax =		dim+2;
	u =		vmuscl->u;
	uL =		vmuscl->uL;
	uR =		vmuscl->uR;

	/* set uL and uR */
	for (k = 0; k < kmax; ++k)
	{
	    uk = u[k]; uLk = uL[k]; uRk = uR[k];
	    for (j = start; j < end; ++j)
	    {
	    	uLk[j] = uk[j - 1];
	    	uRk[j] = uk[j];
	    }
	}

#if defined(DEBUG_MUSCL)
	if (debugging("half_step"))
	{
	    int i, j;

	    g_printout_vec_data("Left state on the edge",
				uL[vmuscl->index.density],
				uL[vmuscl->index.energy],
				uL+vmuscl->index.v[0],dim,start,end,"muncons");
	    g_printout_vec_data("Right state on the edge",
				uR[vmuscl->index.density],
				uR[vmuscl->index.energy],
				uR+vmuscl->index.v[0],dim,start,end,"muncons");

	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","density","int_en");
	    for (i = 0; i < dim; ++i)
	        (void) printf(" %-4s%1d%-9s","vel[",i,"]");
	    (void) printf(" Solution before half step\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",u[vmuscl->index.density][j]);
		(void) printf(" %-14g",u[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",u[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");

	    (void) printf(" Solution after half step\n");
	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","l_density","l_int_en");
	    for (i = 0; i < dim; ++i)
		(void) printf(" %-4s%1d%-9s","l_vel[",i,"]");
	    (void) printf(" %-14s %-14s","r_density","r_int_en");
	    for (i = 0; i < dim; ++i)
		(void) printf(" %-4s%1d%-9s","r_vel[",i,"]");
	    (void) printf("\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",uL[vmuscl->index.density][j]);
		(void) printf(" %-14g",uL[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",uL[vmuscl->index.v[k]][j]);
		(void) printf(" %-14g",uR[vmuscl->index.density][j]);
		(void) printf(" %-14g",uR[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",uR[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_first_order_godunov_half_step*/

LOCAL void print_eigen(
	int		n1,
	int		n2,
	Vec_Eigen	*vegn,
	int		dim)
{
	int		i, j, k, kmax = dim+2;
	double		**lambda = vegn->lambda;
	double		***l = vegn->l, ***r = vegn->r;

	(void) printf("PRINTOUT OF EIGENVALUES AND EIGENVECTORS\n");
	for (j = n1; j < n2; ++j)
	{
	    (void) printf("Eigenvalues and eigenuni_arrays at index %d\n",j);
	    (void) printf("\tEigenvalues - ");
	    for (k = 0; k < kmax; ++k)
	    	(void) printf("%-14g",lambda[k][j]);
	    (void) printf("\n\tLeft eigenuni_array bi_array, L\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	    (void) printf("%-14g ",l[i][k][j]);
	    	(void) printf("]\n");
	    }
	    (void) printf("\tRight eigenuni_array bi_array, R\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	    (void) printf("%-14g ",r[i][k][j]);
	    	(void) printf("]\n");
	    }
	    (void) printf("\tL*R\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	{
	    	    int	n;
	    	    double	tmp = 0.0;

	    	    for (n = 0; n < kmax; ++n)
	    	    	tmp += l[i][n][j]*r[n][k][j];
	    	    (void) printf("%-14g ",tmp);
	    	}
	    	(void) printf("]\n");
	    }
	    (void) printf("\tR*L\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	{
	    	    int	n;
	    	    double	tmp = 0.0;

	    	    for (n = 0; n < kmax; ++n)
	    	    	tmp += r[i][n][j]*l[n][k][j];
	    	    (void) printf("%-14g ",tmp);
	    	}
	    	(void) printf("]\n");
	    }
	    (void) printf("End Eigenvalues and eigenuni_arrays ");
	    (void) printf("at index %d\n\n",j);
	}
	(void) printf("END PRINTOUT OF EIGENVALUES AND EIGENVECTORS\n\n\n");
}		/*end print_eigen*/

LOCAL	void print_linear_reconstruction(
	int	  start,
	int	  end,
	Vec_Muscl *vmuscl)
{
	int	dim = vmuscl->dim, i, j, k, kmax;
	int	negn = vmuscl->Vegn.negn;
	double	**dq = vmuscl->dq, *dqk;
	double	*qk;
	double	*qLk, *qRk;
	double	**u = vmuscl->u;
	double	***l = vmuscl->Vegn.l, ***r = vmuscl->Vegn.r;
	static	double ***rm1 = NULL;
	static  double **q = NULL, **qL = NULL, **qR = NULL,
	              **uL = NULL, **uR = NULL;
	static  int   vsize = 0;

	kmax = dim+2;
	if (rm1 == NULL)
	    bi_array(&rm1,negn,negn,sizeof(double*));

	if (vsize < vmuscl->max_vsize)
	{
	    if (q != NULL)
	        free_these(5,q,qL,qR,uL,uR);
	    vsize = vmuscl->max_vsize;
	    bi_array(&q,kmax,vsize,FLOAT);
	    bi_array(&qL,kmax,vsize,FLOAT);
	    bi_array(&qR,kmax,vsize,FLOAT);
	    bi_array(&uL,kmax,vsize,FLOAT);
	    bi_array(&uR,kmax,vsize,FLOAT);
	}

	g_printout_vec_data("Monotonic MUSCL slopes",
			    dq[0],dq[1],dq+2,dim,start,end,"");

	left_multiply_state_by_matrix(q,l,u,dim,start,end);

	for (k = 0; k < kmax; ++k)
	{
	    qk = q[k]; qLk = qL[k]; qRk = qR[k];
	    dqk = dq[k];
	    for (j = start+1; j < end; ++j)
	    {
	    	qLk[j] = qk[j-1] + 0.5 *dqk[j-1];
	    	qRk[j] = qk[j]	 - 0.5 *dqk[j];
	    }
	}

	(void) output();
	(void) printf("%-4s %-14s %-14s","n","density","int_en");
	for (i = 0; i < dim; ++i)
	    (void) printf(" %-4s%1d%-9s","vel[",i,"]");
	(void) printf(" Input solution in physical state coordinates\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    (void) printf(" %-14g",u[vmuscl->index.density][j-1]);
	    (void) printf(" %-14g",u[vmuscl->index.energy][j-1]);
	    for (k = 0; k < dim; ++k)
	    	(void) printf(" %-14g",u[vmuscl->index.v[k]][j-1]);
	    (void) printf("\n");
	}
	(void) printf("\n");

	(void) output();
	(void) printf("%-4s","n");
	for (k = 0; k < kmax; ++k)
	    (void) printf(" u%-13d",k);
	(void) printf(" Input solution in eigen coordinates\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",q[k][j-1]);
	    (void) printf("\n");
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",q[k][j]);
	    (void) printf("\n");
	}
	(void) printf("\n");

	(void) output();
	(void) printf("%-4s","n");
	for (k = 0; k < kmax; ++k)
	    (void) printf(" u%-13d",k);
	(void) printf(" Reconstructed solution in eigen coordinates\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",qL[k][j]);
	    (void) printf("\n");
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",qR[k][j]);
	    (void) printf("\n");
	}
	(void) printf("\n");

	for (k = 0; k < negn; ++k)
	    for (j = 0; j < negn; ++j)
	    	rm1[k][j] = r[k][j]-1;
	left_multiply_state_by_matrix(uL,rm1,qL,dim,start+1,end);
	left_multiply_state_by_matrix(uR,r,qR,dim,start+1,end);

	(void) printf(" Reconstructed solution in "
	              "physical state coordinates\n");
	(void) output();
	(void) printf("%-4s %-14s %-14s","n","l_density","l_int_en");
	for (i = 0; i < dim; ++i)
	    (void) printf(" %-4s%1d%-9s","l_vel[",i,"]");
	(void) printf(" %-14s %-14s","r_density","r_int_en");
	for (i = 0; i < dim; ++i)
	    (void) printf(" %-4s%1d%-9s","r_vel[",i,"]");
	(void) printf("\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    (void) printf(" %-14g",uL[vmuscl->index.density][j]);
	    (void) printf(" %-14g",uL[vmuscl->index.energy][j]);
	    for (k = 0; k < dim; ++k)
	    	(void) printf(" %-14g",uL[vmuscl->index.v[k]][j]);
	    (void) printf(" %-14g",uR[vmuscl->index.density][j]);
	    (void) printf(" %-14g",uR[vmuscl->index.energy][j]);
	    for (k = 0; k < dim; ++k)
	    	(void) printf(" %-14g",uR[vmuscl->index.v[k]][j]);
	    (void) printf("\n");
	}
	(void) printf("\n");
}		/*end print_linear_reconstruction*/

LOCAL	void	limit_strong_wave_in_cell(
	Vec_Muscl	*vmuscl,
	int		start,
	int		end)
{
	int   dim = vmuscl->dim, i, j, k;
	double **du = vmuscl->du;
	double **dq = vmuscl->dq;
	double ***r = vmuscl->Vegn.r;
	double **lambda = vmuscl->Vegn.lambda;
	static double *dq_p,*dq_n;
	static double *du_p,*du_n;

	if (dq_p == NULL)
        {
	    uni_array(&dq_p,3,FLOAT);
	    uni_array(&dq_n,3,FLOAT);
	    uni_array(&du_p,3,FLOAT);
	    uni_array(&du_n,3,FLOAT);
	}
	for (j = start; j < end; ++j)
	{
	    if ((lambda[0][j] > 0.0 && 
	         lambda[1][j] > 0.0 &&
	         lambda[2][j] > 0.0) || 
		(lambda[0][j] <= 0.0 &&
		 lambda[1][j] <= 0.0 && 
		 lambda[2][j] <= 0.0))
		continue;
	    for (k = 0; k < 3; ++k)
	    {
	    	if (lambda[k][j] > 0.0)
		{
		    dq_p[k] = dq[k][j];
		    dq_n[k] = 0.0;
		}
		else
		{
		    dq_p[k] = 0.0;
		    dq_n[k] = dq[k][j];
		}
	    }
	    for (k = 0; k < 3; ++k)
	    {
	    	du_p[k] = dq_p[0]*r[k][0][j] + dq_p[1]*r[k][1][j]
			+ dq_p[2]*r[k][2][j];
	    	du_n[k] = dq_n[0]*r[k][0][j] + dq_n[1]*r[k][1][j]
			+ dq_n[2]*r[k][2][j];
	    }
	    for (k = 0; k < 3; ++k)
	    {
	    	if (fabs(du[k][j]) < fabs(du_p[k]) ||
		    fabs(du[k][j]) < fabs(du_n[k]))
		{
	    	    for (i = 0; i < 3; ++i)
		    	dq[i][j] = 0.0;
		    break;
	    	}
	    }
	}
}	/* end limit_strong_wave_in_cell */

LOCAL	void	left_multiply_state_by_matrix(
	double		**q,
	double		***l,
	double		**u,
	int		dim,
	int		start,
	int		end)
{
	double		*u0, *u1, *u2, *uk;
	double		*lk0, *lk1, *lk2;
	double		*qk;
	int		j, k;
	int		kmax = 2+dim;

	u0 = u[0]; u1 = u[1]; u2 = u[2];
	for (k = 0; k < 3; ++k)
	{
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    qk = q[k];
	    for (j = start; j < end; ++j)
	    	qk[j] = lk0[j]*u0[j] + lk1[j]*u1[j] + lk2[j]*u2[j];
	}
	for (k = 3; k < kmax; ++k)
	{
	    qk = q[k]; uk = u[k];
	    for (j = start; j < end; ++j)
	    	qk[j] = uk[j];
	}
}		/*end left_multiply_state_by_matrix*/




/*
*			g_compute_art_visc_coefs();
*/

LOCAL void g_compute_art_visc_coefs(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc = vmuscl->avisc;
	double		*g0;
	Locstate	*state;
	Gas_param	*prms;
	double		**lambda;
	double		dlambda, max_dlambda;
	int		i, j;
	int		negn;

	if (avisc == NULL)
	    return;
	if (avisc->use_lapidus)
	    g_compute_lapidus_art_visc_coefs(start,end,vmuscl);
	else
	{
	    g0 = avisc->g[0];
	    for (j = start; j < end; ++j)
	    	g0[j] = 0.0;
	}
	state = vmuscl->vst->state + vmuscl->offset;
	g0 = avisc->g[0];
	if (avisc->use_linear)
	{
	    double *c = vmuscl->vst->c;
	    double *vn = vmuscl->vst->v[0];
	    prms = NULL;
	    for (j = start; j < end; ++j)
	    {
	        if (Params(state[j]) != prms)
	    	    prms = Params(state[j]);
	        if (prms != NULL)
		{
		    double sp = (fabs(vn[j]) + c[j])*prms->avisc.sp_coef;;
	            g0[j] += 2.0*sp*prms->avisc.linear_visc_coef;
		}
	    }
	}
	if (avisc->use_upwind)
	{
	    lambda = vmuscl->Vegn.lambda;
	    negn = vmuscl->Vegn.negn;
	    prms = NULL;
	    for (j = start; j < end; ++j)
	    {
	        if (Params(state[j]) != prms)
	            prms = Params(state[j]);
	        if (prms == NULL)
	            continue;
	        max_dlambda = 0.0;
	        for (i = 0; i < negn; ++i)
	        {
	            dlambda = lambda[i][j] - lambda[i][j+1];
	            if (dlambda > max_dlambda)
	                max_dlambda = dlambda;
	        }
	        g0[j] += 2.0*prms->avisc.upwind_visc_coef*max_dlambda;
	    }
	}
}		/*end g_compute_art_visc_coefs*/

LOCAL void g_compute_lapidus_art_visc_coefs(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc = vmuscl->avisc;
	Vec_Eigen	*vegn;
	Vec_Gas		*vst;
	int		offset;
	double		**g;
	double		*cs_ave;
	double		*c_ave;
	double		*vn_ave;
	double		**b;
	double		*visc;
	double		*vn;
	double		*c2;
	double		*bi, *b0, *b1, *b2;
	double		*g0, *g1, *g2;
	Locstate	*state;
	Gas_param	*prms;
	double		**lambda, *lambdai;
	int		i, j, negn;

	        /* Read arrays from Vec_Muscl structure */
	offset = vmuscl->offset;
	vst = vmuscl->vst;
	vegn = &vmuscl->Vegn;
	g = avisc->g;
	cs_ave = avisc->cs_ave;
	c_ave = avisc->c_ave;
	vn_ave = avisc->vn_ave;
	b = avisc->b;
	visc = avisc->visc;
	vn = vst->v[0] + offset;
	c2 = vst->c2 + offset;
	state = vst->state + offset;
	lambda = vegn->lambda;
	negn = vegn->negn;

	prms = NULL;
	for (j = start; j < end; ++j)
	{
	    if (Params(state[j]) != prms)
	        prms = Params(state[j]);
	    visc[j] = (prms == NULL) ? 0.0 : prms->avisc.lapidus_visc_coef;
	}

	for (j = start; j < end; ++j)
	{
	    cs_ave[j] = 0.5*(c2[j+1] + c2[j]);
	    vn_ave[j] = 0.5*(vn[j+1] + vn[j]);
	}
	for (j = start; j < end; ++j)
	    c_ave[j] = sqrt(cs_ave[j]);
	
	for (i = 0; i < negn; ++i)
	{
	    bi = b[i];
	    lambdai = lambda[i];
	    for (j = start; j < end; ++j)
	        bi[j] = visc[j]*fabs(lambdai[j+1] - lambdai[j]);
	}

	b0 = b[0]; b1 = b[1]; b2 = b[2];
	g0 = g[0]; g1 = g[1]; g2 = g[2];
	for (j = start; j < end; ++j)
	{
	    g2[j] = 0.5*(b0[j] - 2.0*b1[j] + b2[j])/cs_ave[j];
	    g1[j] = -0.5/cs_ave[j]*(
	        		(2.0*vn_ave[j]+c_ave[j])*b0[j] -
	        		           4.0*vn_ave[j]*b1[j] +
	        		(2.0*vn_ave[j]-c_ave[j])*b2[j]);
	    g0[j] = 0.5/cs_ave[j]*(
	        		vn_ave[j]*(vn_ave[j]+c_ave[j])*b0[j] -
	        		2.0*(sqr(vn_ave[j])-cs_ave[j])*b1[j] +
	        		vn_ave[j]*(vn_ave[j]-c_ave[j])*b2[j]);
	}
#if defined(DEBUG_MUSCL)
	if (debugging("art_visc"))
	{
	    (void) printf("Artificial viscosity coefficients\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","g0","g1","g2");
	    for (j = start; j < end; ++j)
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,g0[j],g1[j],g2[j]);
	    (void) printf("\nEnd Artificial viscosity coefficients\n");
	    (void) printf("Wave speed gradients\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","b0","b1","b2");
	    for (j = start; j < end; ++j)
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,b0[j],b1[j],b2[j]);
	    (void) printf("End Wave speed gradients\n");
	    (void) printf("Mid wave speeds\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","v-c","v","v+c");
	    for (j = start; j < end; ++j)
	    {
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,vn_ave[j]-c_ave[j],
	        	      vn_ave[j],vn_ave[j]+c_ave[j]);
	    }
	    (void) printf("End Mid wave speeds\n");
	    (void) printf("Interpolation test\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n",
	        	  "n","G(v-c) - b0","G(v) - b1","G(v+c) - b2");
	    for (j = start; j < end; ++j)
	    {
	        double	vmc, v, vpc;
	        double	G0, G1, G2;

	        vmc = vn_ave[j] - c_ave[j];
	        v = vn_ave[j];
	        vpc = vn_ave[j] + c_ave[j];
	        G0 = g0[j] + g1[j]*vmc + g2[j]*vmc*vmc;
	        G1 = g0[j] + g1[j]*v + g2[j]*v*v;
	        G2 = g0[j] + g1[j]*vpc + g2[j]*vpc*vpc;
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,G0-b0[j],
	    		  G1-b1[j],G2-b2[j]);
	    }
	    (void) printf("End Interpolation test\n");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_compute_lapidus_art_visc_coefs*/


/*ARGSUSED*/
LOCAL	void g_muscl_lin_approx_rsolver(
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
	g_lin_rsoln(start,end,uL,uR,uM,vmuscl);

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps,i,j;
            if((num_comps = params->n_comps) != 1)
            {
                for (j = start; j < end; ++j)
                {
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
	g_load_muscl_flux(start,end,uM,vmst,Flux,vmuscl);
}		/*end g_muscl_lin_approx_rsolver*/


/*
*		      g_lin_rsoln():
*
*	This routine is the approximate Riemann solver. It solves the
*	linearized system. It does not handle sonic rarefactions properly.
*/

LOCAL void g_lin_rsoln(
	int	  start,
	int	  end,
	double     **uL,
	double     **uR,
	double     **uM,
	Vec_Muscl *vmuscl)
{
	double		**awv;
	Vec_Eigen	*vegn;
	int		j, k;
	int		kmax;
	double		***l, ***r;
	double		**lambda;
	double		*lambdak;
	double		**sgnm;
	double		*awv0, *awv1, *awv2, *awvk;
	double		*sgnmk, *sgnm1;
	double		*uR0, *uR1, *uR2, *uRk;
	double		*uL0, *uL1, *uL2, *uLk;
	double		*uMk;
	double		*lk0, *lk1, *lk2;
	double		*rk0, *rk1, *rk2;

	kmax =		vmuscl->dim+2;
	vegn =		&vmuscl->Vegn;
	l =		vegn->l;
	r =		vegn->r;
	sgnm =		vegn->sgnm;
	lambda =	vegn->lambda;
	awv =		vmuscl->awv;

	/* Do left contributions */
	for (k = 0; k < 3; ++k)
	{
	    sgnmk = sgnm[k];
	    lambdak = lambda[k];
	    for (j = start; j < end; ++j)
	    {
	    	sgnmk[j] = (lambdak[j-1] < 0.0) ? 1.0 : 0.0;
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    awvk = awv[k];
	    sgnmk = sgnm[k];
	    uR0 = uR[0]; uR1 = uR[1]; uR2 = uR[2];
	    uL0 = uL[0]; uL1 = uL[1]; uL2 = uL[2];
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnmk[j]*( (uR0[j] - uL0[j])*lk0[j-1] +
	    		             (uR1[j] - uL1[j])*lk1[j-1] +
	    		             (uR2[j] - uL2[j])*lk2[j-1]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    awvk = awv[k];
	    sgnm1 = sgnm[1];
	    uRk = uR[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnm1[j]*(uRk[j] - uLk[j]);
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    uMk = uM[k];
	    uLk = uL[k];
	    uRk = uR[k];
	    rk0 = r[k][0]; rk1 = r[k][1]; rk2 = r[k][2];
	    awv0 = awv[0]; awv1 = awv[1]; awv2 = awv[2];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] = 0.5*(uLk[j] + ( awv0[j]*rk0[j-1] +
	        			 awv1[j]*rk1[j-1] +
	        			 awv2[j]*rk2[j-1]));
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    uLk = uL[k];
	    uRk = uR[k];
	    uMk = uM[k];
	    awvk = awv[k];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] = 0.5*(uLk[j] + awvk[j]);
	    }
	}

	/* Do right contributions */
	for (k = 0; k < 3; ++k)
	{
	    sgnmk = sgnm[k]; lambdak = lambda[k];
	    for (j = start; j < end; ++j)
	    {
	    	sgnmk[j] = (lambdak[j] > 0.0) ? -1.0 : 0.0;
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    awvk = awv[k];
	    sgnmk = sgnm[k];
	    uR0 = uR[0]; uR1 = uR[1]; uR2 = uR[2];
	    uL0 = uL[0]; uL1 = uL[1]; uL2 = uL[2];
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnmk[j]*( (uR0[j] - uL0[j])*lk0[j] +
	        		     (uR1[j] - uL1[j])*lk1[j] +
	        		     (uR2[j] - uL2[j])*lk2[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    awvk = awv[k];
	    sgnm1 = sgnm[1];
	    uRk = uR[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnm1[j]*(uRk[j] - uLk[j]);
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    uMk = uM[k];
	    uLk = uL[k];
	    uRk = uR[k];
	    rk0 = r[k][0]; rk1 = r[k][1]; rk2 = r[k][2];
	    awv0 = awv[0]; awv1 = awv[1]; awv2 = awv[2];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] += 0.5*(uRk[j] + ( awv0[j]*rk0[j] +
	    			          awv1[j]*rk1[j] +
	    			          awv2[j]*rk2[j]));
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    uLk = uL[k];
	    uRk = uR[k];
	    uMk = uM[k];
	    awvk = awv[k];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] += 0.5*(uRk[j] + awvk[j]);
	    }
	}
}		/*end g_lin_rsoln*/

/*
*			g_add_art_visc1();
*
*	This function adds the first half of artificial viscosity.  It is
*	added to the mid states which will be used to compute conservative
*	fluxes and then in the conservative difference scheme.
*/

LOCAL void g_add_art_visc1(
	int		start,
	int		end,
	double           **uM,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc;
	int		dim;
	int		j, k;
	double		*g2, *g1;
	double		**u, **ucon;
	double		*p, **F, **g, **uconM;
	double		*uconMk, *uconk, *Fk;
	double		du;
	int		id = vmuscl->index.density;
	int		ie = vmuscl->index.energy;

	if ((avisc = vmuscl->avisc) == NULL)
	    return;
	g = avisc->g;
	g2 = g[2], g1 = g[1];
	if (g2 == NULL || g1 == NULL)
	    return;

	u =		vmuscl->u;
	dim =		vmuscl->dim;
	ucon =		vmuscl->ucon;
	F =		vmuscl->Flux.F;
	p =		vmuscl->vst->p + vmuscl->offset;
	uconM =		avisc->uconM;

	flux_vectors(start,end,u,p,&vmuscl->Flux,vmuscl);
	g_uncons_to_cons(start,end,dim,uM,uconM);

	for (k = 0; k < dim+2; ++k)
	{
	    uconk = ucon[k];
	    uconMk = uconM[k];
	    Fk = F[k];
	    for (j = start; j < end; ++j)
	    {
		du = 0.5*(g2[j-1]*(Fk[j] - Fk[j-1]) +
		                  g1[j-1]*(uconk[j] - uconk[j-1]));
	    	if ((k == id || k == ie) && du >= uconMk[j])
                {
                    const char *name = (k == id) ? "density" : "energy";
                    screen("WARNING in g_add_art_visc1():\n");
                    screen("nonlinead artificial viscosity causes middle ");
                    screen("state negative %s,\n",name);
                    screen("setting to use default (no viscosity)\n");
                }
		else
		    uconMk[j] -= du;
	    }
	}
	g_cons_to_uncons(start,end,dim,uconM,uM);
}		/*end g_add_art_visc1*/

/*
*			g_add_art_visc2();
*
*	This function adds the second half of the artificial viscosity.	 All
*	we need do here is to modify appropriately the flux uni_arrays already
*	computed.  Thus the conservative differencing algorithm does not 
*	need to be modified.
*/

LOCAL void g_add_art_visc2(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc;
	int		dim, j, k;
	double		*g0;
	double		*uconk, *Fk;
	double		**ucon, **F;

	if ((avisc = vmuscl->avisc) == NULL)
	    return;
	dim = vmuscl->dim;
	ucon = vmuscl->ucon;
	F = vmuscl->Flux.F;
	g0 = avisc->g[0];

	for (k = 0; k < dim+2; ++k)
	{
	    uconk = ucon[k];
	    Fk = F[k];
	    for (j = start; j < end; ++j)
	    	Fk[j] -= 0.5*g0[j-1]*(uconk[j] - uconk[j-1]);
	}
}		/*end g_add_art_visc2*/


LOCAL void g_uncons_to_cons(
	int		start,
	int		end,
	int		dim,
	double		**uncons,
	double		**cons)
{
	int		j;
	double		*rho = uncons[0],	*cons0 = cons[0];
	double		*e = uncons[1],		*en_den = cons[1];
	double		*v0 = uncons[2],	*m0 = cons[2];
	double		*v1, *v2, *m1, *m2, v_norm;

	switch (dim)
	{
	case 1:
	    for (j = start; j < end; ++j)
	    {
	    	m0[j] = v0[j]*rho[j];
	    	v_norm = sqr(v0[j]);
	    	en_den[j] = (e[j] + 0.5*v_norm)*rho[j];
	    	cons0[j] = rho[j];
	    }
	    break;
	case 2:
	    v1 = uncons[3]; m1 = cons[3];
	    for (j = start; j < end; ++j)
	    {
	    	m1[j] = v1[j]*rho[j];
	    	m0[j] = v0[j]*rho[j];
	    	v_norm = sqr(v1[j]) + sqr(v0[j]);
	    	en_den[j] = (e[j] + 0.5*v_norm)*rho[j];
	    	cons0[j] = rho[j];
	    }
	    break;
	case 3:
	    v1 = uncons[3]; m1 = cons[3];
	    v2 = uncons[4]; m2 = cons[4];
	    for (j = start; j < end; ++j)
	    {
	    	m2[j] = v2[j]*rho[j];
	    	m1[j] = v1[j]*rho[j];
	    	m0[j] = v0[j]*rho[j];
	    	v_norm = sqr(v2[j]) + sqr(v1[j]) + sqr(v0[j]);
	    	en_den[j] = (e[j] + 0.5*v_norm)*rho[j];
	        cons0[j] = rho[j];
	    }
	    break;
	}
}		/*end g_uncons_to_cons*/

LOCAL void g_cons_to_uncons(
	int		start,
	int		end,
	int		dim,
	double		**cons,
	double		**uncons)
{
	int		j;
	double		*rho = cons[0],		*uncons0 = uncons[0];
	double		*en_den = cons[1],	*e = uncons[1];
	double		*m0 = cons[2],		*v0 = uncons[2];
	double		*m1, *m2, *v1, *v2, v_norm;

	switch (dim)
	{
	case 1:
	    for (j = start; j < end; ++j)
	    {
	    	v0[j] = m0[j]/rho[j];
	    	v_norm = sqr(v0[j]);
	    	e[j] = en_den[j]/rho[j] - 0.5*v_norm;
	    	uncons0[j] = rho[j];
	    }
	    break;
	case 2:
	    m1 = cons[3]; v1 = uncons[3];
	    for (j = start; j < end; ++j)
	    {
	    	v1[j] = m1[j]/rho[j];
	    	v0[j] = m0[j]/rho[j];
	    	v_norm = sqr(v0[j]) + sqr(v1[j]);
	    	e[j] = en_den[j]/rho[j] - 0.5*v_norm;
	    	uncons0[j] = rho[j];
	    }
	    break;
	case 3:
	    m1 = cons[3]; v1 = uncons[3];
	    m2 = cons[4]; v2 = uncons[4];
	    for (j = start; j < end; ++j)
	    {
	    	v2[j] = m2[j]/rho[j];
	    	v1[j] = m1[j]/rho[j];
	    	v0[j] = m0[j]/rho[j];
	    	v_norm = sqr(v0[j]) + sqr(v1[j]) + sqr(v2[j]);
	    	e[j] = en_den[j]/rho[j] - 0.5*v_norm;
	    	uncons0[j] = rho[j];
	    }
	    break;
	}
}		/*end g_cons_to_uncons*/


#if defined(DEBUG_MUSCL)

/*
*		EXPLANATION OF DEBUGGING SWITCHES
*
*	The printout of the interior state debugging is so extensive
*	that it is impractical in most runs to printout unconditional
*	debug lines.  Therefore the following control strings for
*	debugging have been added to allow the restriction of debugging
*	to only certain areas of the program.  In particular it is possible
*	to restrict the printout as a function of the grid cell index
*	being computed.  The use of these debugging control strings is
*	enabled by turning on the debugging string MUSCL_SOLVER.
*
*		CONTROLLING DEBUGGING BY SWEEP DIRECTION
*
*	Key debug strings:
*		xonly	yonly	zonly
*
*	If any of these debug strings is turned on,  then the debugging
*	will be restricted to that corresponding sweep direction.  These
*	values are obviously mutually exclusive.
*
*		CONTROLLING DEBUGGING BY SWEEP TYPE
*
*	Key debug strings:
*		mtsten_all	mtsten_some	mtsten_none
*		mirreg_all	mirreg_some	mirreg_none
*		mreg_all	mreg_some	mreg_none
*
*	The debug strings are formed in two parts by combining the
*	prefixes mtsten,  mirreg,  and mreg,  with the suffixes
*	_all, _some, or _none.  The suffix arguments are mutually
*	exclusive with the least restrictive flags overriding the
*	more restrictive flags.  (For example requesting all and some 
*	yields all.)
*	
*	DEBUG STRING PREFIX	EXPLANATION
*
*	mtsten			Controls debugging for the tangential sweep.
*
*	mirreg			Controls debugging for the irregular stencil
*				sweeps.
*
*	mreg			Controls debugging for the regular stencil
*				sweep.
*
*	DEBUG STRING SUFFIX	EXPLANATION
*
*	_all			Turns on all debugging strings for MUSCL code.
*				The printout of these strings will be 
*				controled by the cell index and time step
*				control desribed below.
*
*	_some			Allows controled output of MUSCL code
*				debugging.  If this string is active you
*				can turn on specific debug strings
*				in the MUSCL code (eg cg_rsolve,  etc)
*				and the corresponding debug lines will
*				be controled by the cell index and
*				time step control described below.
*
*	_none			Turns off all debugging for MUSCL code
*				tangential sweep. In this case no debugging
*				output will be printed whether the individual
*				debug strings are turned on or not.
*
*
*		CONTROLLING DEBUGGING BY CELL INDEX
*
*	Key debug strings:
*		mrgrid
*		LX (followed by an integer not separated by a space, eg LX10)
*		LY (followed by an integer not separated by a space, eg LY10)
*		LZ (followed by an integer not separated by a space, eg LZ10)
*		UX (followed by an integer not separated by a space, eg UX10)
*		UY (followed by an integer not separated by a space, eg UY10)
*		UZ (followed by an integer not separated by a space, eg UZ10)
*
*	The grid cell index restriction is communicated by multiple debug 
*	strings.  If the debug string mrgrid is input this should be followed
*	by debug strings of the form LXnn, UXNN,  LYnn,  UYNN,  LZnn, UZNN,
*	where the nn and NN (which need not be all equal) are the indicies of
*	the cells for which debugging is desired to be turned on.  If any of
*	these are omited,  they default to the ends of the computational grid.
*
*			CONTROLLING DEBUGGING BY TIME STEP
*
*	Key debug strings:
*		mrtstep
*		mt (followed by an integer not separated by a space, eg mt10)
*		MT (followed by an integer not separated by a space, eg MT10)
*
*	The time step restriction is also communicated by multiple debug
*	strings.  If the debug string mrtstep is input this should be followed
*	by the debug strings mtnn and MTNN,  where nn is the first time to
*	print debuggin and NN is the last.
*
*	Example:	The following collection of debug strings
*	will print only output during the irregular sweep for cells with
*	indicies ix = 28 and iy = 50 for times 0 and 1.
*
*	: MUSCL_SOLVER
*	: mreg_none
*	: mtsten_none
*	: mirreg_all
*	: mrgrid
*	: LY28
*	: UY28
*	: LX50
*	: UX50
*	: mrtstep
*	: mt0
*	: MT1
*/


LOCAL	const char **set_muscl_debugging(
	int		sn,
	int		*ip,
	int		*ic,
	Tan_stencil	*tsten,
	Front		*fr,
	int		vsize)
{
	const char	**debug_strings = NULL;
	int		dim;

	if (!debugging("MUSCL_SOLVER"))
	    return NULL;
	if (tsten != NULL)
	{
	    if (debugging("mtsten_all"))
	    {
	    	debug_strings = toggle_muscl_debugging("ALL",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mtsten_some"))
	    {
	    	debug_strings = toggle_muscl_debugging("SOME",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mtsten_none"))
	    {
	    	debug_strings = toggle_muscl_debugging("NONE",sn,ip,ic,fr,NO);
	    }
	    if (debugging("oned_MUSCL"))
	    {
	    	(void) printf("oned_MUSCL called from oblique solver.\n");
	    	dim = fr->rect_grid->dim;
	    	print_general_vector("At coords ",Coords(tsten->p[0]),dim,
	        		     "\n\n");
	    }
	}
	else if (vsize == 5)
	{
	    if (debugging("mirreg_all"))
	    {
	    	debug_strings = toggle_muscl_debugging("ALL",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mirreg_some"))
	    {
	    	debug_strings = toggle_muscl_debugging("SOME",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mirreg_none"))
	    {
	    	debug_strings = toggle_muscl_debugging("NONE",sn,ip,ic,fr,NO);
	    }
	    if (debugging("oned_MUSCL"))
	    {
	    	dim = fr->rect_grid->dim;
	    	(void) printf("oned_MUSCL called from irregular solver.\n");
	    	(void) printf("Sweep direction = %d, ",ip[sn]);
	    	print_int_vector("at ic ",ic,dim,"\n");
	    	(void) printf("\n");
	    }
	}
	else 
	{
	    if (debugging("mreg_all"))
	    {
	    	debug_strings = toggle_muscl_debugging("ALL",sn,ip,ic,fr,YES);
	    }
	    else if (debugging("mreg_some"))
	    {
	    	debug_strings = toggle_muscl_debugging("SOME",sn,ip,ic,fr,YES);
	    }
	    else if (debugging("mreg_none"))
	    {
	    	debug_strings = toggle_muscl_debugging("NONE",sn,ip,ic,fr,YES);
	    }
	    if (debugging("oned_MUSCL"))
	    {
	    	dim = fr->rect_grid->dim;
	    	(void) printf("oned_MUSCL called from vector solver.\n");
	    	(void) printf("Sweep direction = %d, ",ip[sn]);
	    	print_int_vector("at ic ",ic,dim,"\n");
	    	(void) printf("\n");
	    }
	}
	return debug_strings;
}		/*end set_muscl_debugging*/

LOCAL	const char **toggle_muscl_debugging(
	const char	*howmuch,
	int		swp_num,
	int		*iperm,
	int		*icoords,
	Front		*fr,
	int		reg_sweep)
{
	int		  i, j, idir;
	int		  dim = fr->rect_grid->dim;
	int		  ndb;
	long		  int mt, MT;
	long		  int mg, MG;
	char		  **db;
	const char        **c, **s;
	static const char *debug_strings[100];
	static const char        *swpdb[3] = {"xonly", "yonly", "zonly"};
	static const char        *mnames[3] = {"LX", "LY", "LZ"};
	static const char        *Mnames[3] = {"UX", "UY", "UZ"};
	static const char        *muscl_debug_strings[] = { "cg_rsolve",
	        					    "MUSCL",
	        					    "oned_MUSCL",
	        					    "MUSCLob",
	        					    "lriem",
	        					    "vflux",
	        					    "load_Vgas",
	        					    "eigen",
	        					    "lcnst",
	        					    "half_step",
	        					    "csrc",
	        					    "vgdnv",
	        					    "vgdnv_ie",
	        					    "art_visc",
	        					    NULL};

	if (iperm != NULL)
	{
	    idir = iperm[swp_num];
	    if (!debugging(swpdb[idir]))
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    if (i == idir)
	        	continue;
	            if (debugging(swpdb[i]))
	        	howmuch = "NONE";
	        }
	    }
	    if (debugging("mrgrid"))
	    {
	    	db = debugging_names(&ndb);
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    if (reg_sweep == YES && i == idir)
	        	continue;
	    	    mg = 0;
	    	    MG = INT_MAX;
	    	    for (j = 0; j < ndb; ++j)
	    	    {
	    	        if (strncmp(db[j],mnames[i],2) == 0)
	    	            mg = atoi(db[j]+2);
	    	        else if (strncmp(db[j],Mnames[i],2) == 0)
	    	    	    MG = atoi(db[j]+2);
	    	    }
	    	    if (icoords[i] < mg || icoords[i] > MG)
	    	    	howmuch = "NONE";
	    	}
	    }
	}

	if (debugging("mrtstep"))
	{
	    mt = 0;
	    MT = INT_MAX;
	    db = debugging_names(&ndb);
	    for (i = 0; i < ndb; ++i)
	    {
	    	if (strncmp(db[i],"mt",2) == 0)
	    	    mt = atoi(db[i]+2);
	    	else if (strncmp(db[i],"MT",2) == 0)
	    	    MT = atoi(db[i]+2);
	    }
	    if (fr->step < mt || fr->step > MT)
	    	howmuch = "NONE";
	}

	c = debug_strings;
	if (strcmp(howmuch,"ALL") == 0)
	{
	    for (s = muscl_debug_strings; *s != NULL; ++s)
	    {
	    	if (debugging(*s))
	    	    continue;
	    	add_to_debug(*s);
	    	*c++ = "TURNED_ON";
	    	*c++ = *s;
	    }
	}
	if (strcmp(howmuch,"NONE") == 0)
	{
	    for (s = muscl_debug_strings; *s != NULL; ++s)
	    {
	    	if (!debugging(*s))
	            continue;
	    	remove_from_debug(*s);
	    	*c++ = "TURNED_OFF";
	    	*c++ = *s;
	    }
	}
	*c = NULL;
	return debug_strings;
}		/*end toggle_muscl_debugging*/

LOCAL	void reset_muscl_debugging(
	const char **debug_strings)
{
	const char **s;

	if (debug_strings == NULL)
	    return;

	for (s = debug_strings; *s != NULL; ++s)
	{
	    if (strcmp(*s++,"TURNED_OFF") == 0)
	    	add_to_debug(*s);
	    else
	    	remove_from_debug(*s);
	}
}		/*end reset_muscl_debugging*/


#endif /* defined(DEBUG_MUSCL) */
