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
*				gplm.c
*
*	Copyright 2000 by The University at Stony Brook, All rights reserved.
*
*	An implementation of the Colella Piecewise Linear Method.  See
*	"A Direct Eulerian MUSCL Scheme for Gas Dynamics"
*       P. Colella, SIAM J. Sci. Stat. Comput., Vol. 6, No. 1 (1985).
*
*
*	NOTE ON INDEXING:
*
*		Interior states are indexed in a row in the sweep direction.
*	If we let j denote the index in the sweep direction,  then cell
*	centered quantities (ie the input vst and src structures)
*	have index j.  In terms of the input arguments to oned_PLM() below,
*	the cell centered arrays have length vsize and start at an offset
*	offset from the arrays in the vst structure.  For example the density
*	array rho = vst->rho + offset. Cell edges are indexed so cell edge j
*	is the left cell edge of cell j.   Thus indices j in arrays associated
*	with the edge quantities (VMst, VRst, VLst, VM0st, VM1st) correspond
*	to the left cell edge of cell j.   Note that VRst,  and VLst contain
*	state information on the right and left of cell edge j respectively.
*/

#define DEBUG_PLM

#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>

/* LOCAL Function Prototypes */
LOCAL	Vec_Muscl	*g_load_plm_state_data(Muscl_Opts*,int,int*,Front*,
                                               Wave*,Stencil*,
					       Tan_stencil*,Vec_Gas*,Vec_Src*,
					       int,int,int,double,double);
LOCAL	Vec_Muscl *plm_alloc_wk_space(Muscl_Opts*,int,int,Vec_Muscl*);
LOCAL	double	  *set_plm_vec_gas_storage(double**,Vec_Gas*,double**,Vec_Muscl*);
LOCAL	void	  characteristic_traceback(int,int,Vec_Muscl*);
LOCAL	void	  g_plm_add_art_visc1(int,int,double**,Vec_Muscl*);
LOCAL	void	  g_plm_add_art_visc2(int,int,Vec_Muscl*);
LOCAL	void	  g_plm_compute_eigens(int,int,Vec_Muscl*);
LOCAL	void	  g_plm_flux_vectors(int,int,double**,double*,MUSCL_FLUX*,
                                     Vec_Muscl*);
LOCAL	void	  g_plm_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	  modify_flow_field_slope(double*,const double*,int,int);
LOCAL	void	  monotonize_reconstruction(const double*,double*,int,int);
LOCAL	void	  reconstruct_flow_field(const double*,double*,int,int);
LOCAL	void	  set_derived_field_slope_and_mid(const double*,const double*,
	                                          double*,double*,int,int);

/* Debugging Function Prototypes */
#if defined(DEBUG_PLM)
LOCAL	boolean	plm_interior_sten_reg(int,int,Stencil*);
LOCAL	void	print_data_header(int dim,const char*,const char*,const char*,
	                          const char*,const char*,const char*);
LOCAL	void	print_data_line(int,int,double**,double*,double*,double**,
                                double*,double*);
#endif /* defined(DEBUG_PLM) */

LOCAL  Locstate tmpst = NULL;
/*
*				oned_PLM():
*
*	One dimensional Colella (PLM) MUSCL code. The input data are the
*	conservative state variables vst and the source terms src.
*	The scheme uses a seven point therefore it updates the states from
*	n1+3 to n2-4.  (n2 is the first index that does NOT get processed,
*	thus at a boundary n2-1 is the boundary state).
*/

/*ARGSUSED*/
EXPORT void oned_PLM(
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
	Locstate   *state;
	double      **F, **H, **F0, **F1, **H0, **H1;
	double      *rho, **m, *E;
	double      *rhoF, *rhoH;
	double      **mF, **mH;
	double      *EF, *EH;
	double      *A, *dV;
	double      g, *grav, *gL, *gR;
	double      rho0, m0;
	double      speed;
	double      *c;
	double      *v0;
	double      **coords;
	int        start, end;
	int        j, k, nfloats;
	int        sten_rad;
	double      rmin;
#if defined(DEBUG_PLM)
	SWEEP_TYPE sweep_type = UNSET_SWEEP_TYPE;
#endif /* defined(DEBUG_PLM) */

#if defined(DEBUG_PLM)
	debug_print("plm","Entered oned_PLM()\n");
#endif /* defined(DEBUG_PLM) */
	if(is_rotational_symmetry())
	    rmin = pos_radius(0.0,fr->rect_grid);
	
#if defined(DEBUG_PLM)
	if (tsten)
	    sweep_type = TANGENTIAL_SWEEP;
	else if (sten)
	    sweep_type = IRREGULAR_INTERIOR_SWEEP;
	else
	    sweep_type = REGULAR_INTERIOR_SWEEP;
	
	if (debugging("plm"))
	{
	    (void) printf("dim = %d, offset = %d, vsize = %d, "
	                  "dt = %g, dn = %g, ",dim,offset,vsize,dt,dn);
	    switch (sweep_type)
	    {
	    case TANGENTIAL_SWEEP:
	        (void) printf("TANGENTIAL SWEEP\n");
		break;
	    case IRREGULAR_INTERIOR_SWEEP:
	        (void) printf("NEAR FRONT INTERIOR SWEEP\n");
		break;
	    case REGULAR_INTERIOR_SWEEP:
	        (void) printf("REGULAR INTERIOR SWEEP\n");
		break;
	    case UNSET_SWEEP_TYPE:
	        screen("ERROR in oned_PLM(), unset sweep type\n");
		clean_up(ERROR);
	    }
	}
#endif /* defined(DEBUG_PLM) */

	vmuscl = load_state_data(swp_num,iperm,fr,wave,sten,tsten,
	        		 vst,src,offset,vsize,dim,dt,dn);

	if (tmpst == NULL)
	    alloc_state(vmuscl->front->interf,&tmpst,vmuscl->front->sizest);

	compute_eigens(0,vsize,vmuscl);

#if defined(DEBUG_PLM)
	MusclSweepType(vmuscl) = sweep_type;
#endif /* defined(DEBUG_PLM) */
	sten_rad = vmuscl->sten_rad;
	start = sten_rad - 1;
	end = vsize - start;

	compute_art_visc_coefs(start,end,vmuscl);

	reconstructor(start,end,vmuscl);

	/* Compute states at cell edges at start of time step */
	++start;
	rsolver(start,end,vmuscl->uL,&vmuscl->VLst,vmuscl->uR,&vmuscl->VRst,
	                  vmuscl->uM0,&vmuscl->VM0st,&vmuscl->Flux0,vmuscl);
	load_pressure_and_gammas(&vmuscl->VM0st,start,end-start);

	characteristic_traceback(start,end,vmuscl);

	characteristic_solve(start,end,vmuscl);

	F = vmuscl->Flux.F;   H = vmuscl->Flux.H;
	F0 = vmuscl->Flux0.F; H0 = vmuscl->Flux0.H;
	F1 = vmuscl->Flux1.F; H1 = vmuscl->Flux1.H;
	nfloats = vmuscl->Opts.nfloats;
	for (k = 0; k < nfloats; ++k)
	{
	    for (j = start; j < end; ++j)
	    {
	        F[k][j] = 0.5*(F0[k][j] + F1[k][j]);
	        H[k][j] = 0.5*(H0[k][j] + H1[k][j]);
	    }
	}

	/* Conservative differencing */
	end--;
	rho = vmuscl->ucon[vmuscl->index.density];
	E = vmuscl->ucon[vmuscl->index.energy];
	m = vmuscl->ucon+vmuscl->index.v[0];
	A = vmuscl->A;
	dV = vmuscl->dV;
	rhoF = F[vmuscl->index.density];
	rhoH = H[vmuscl->index.density];
	mF = F + vmuscl->index.v[0];
	mH = H + vmuscl->index.v[0];
	EF = F[vmuscl->index.energy];
	EH = H[vmuscl->index.energy];

	grav = vmuscl->grav;              /* Gravity at time t  cell center */

	/* Gravities at time t+dt at cell edges */
	gL = vmuscl->uL[vmuscl->index.grav];
	gR = vmuscl->uR[vmuscl->index.grav];

	add_art_visc2(start,end,vmuscl);
	coords = vst->coords + offset;
	v0 = vst->v[0] + offset;
	for (j = start; j < end; ++j)
	{
	    rho0 = rho[j];
	    m0 = m[0][j];
	    if (is_rotational_symmetry() &&
		(vmuscl->alpha != 0.0) && (fabs(coords[j][0]) < rmin))
	    {
	        Locstate st0 = vst->state[j+offset];
		double    ke;
	        rho[j] -= (dt/dn)*0.5*(rho[j+1]+rho[j])*(v0[j+1] - v0[j]);
		state_on_adiabat_with_dens(st0,rho[j],tmpst,EGAS_STATE);
	        ke = m[0][j] = 0.0;
	        for (k = 1; k < dim; ++k)
		{
	            m[k][j] = rho[j]*vst->v[k][j];
		    ke += 0.5*m[k][j]*vst->v[k][j];
		}
		E[j] = rho[j]*Energy(tmpst)+ke;
	    }
	    else
	    {
	        rho[j] -= (dt/dV[j])*(A[j+1]*rhoF[j+1] - A[j]*rhoF[j]) +
	                  (dt/dn)*(rhoH[j+1] - rhoH[j]);
	        for (k = 0; k < dim; ++k)
	            m[k][j] -= (dt/dV[j])*(A[j+1]*mF[k][j+1] - A[j]*mF[k][j]) +
		               (dt/dn)*(mH[k][j+1] - mH[k][j]);

	        E[j] -= (dt/dV[j])*(A[j+1]*EF[j+1] - A[j]*EF[j]) +
		        (dt/dn)*(EH[j+1] - EH[j]);
	    }
	    g = 0.5*(grav[j] + 0.5*(gR[j]+gL[j+1]));
	    m[0][j] += dt*0.5*(rho0 + rho[j])*g;
	    E[j] += dt*0.5*(m0 + m[0][j])*g;
	}

	state = vst->state + offset;
	coords = vst->coords + offset;
	for (j = start; j < end; ++j)
	{
	    boolean use_godunov = NO;
	    if ((rho[j] < 0.0) || (E[j] < 0.0))
	        use_godunov = YES;
	    else
	    {
	        double e = E[j];
		for (k = 0; k < dim; ++k)
		    e -= 0.5*m[k][j]*m[k][j]/rho[j];
	        if (e < Min_energy(vst->state[j]))
		    use_godunov = YES;
	    }
	    if (use_godunov)
	    {
		if (tsten != NULL)
		{
		    if (debugging("godunov"))
		    {
	                double e = E[j];
			for (k = 0; k < dim; ++k)
		            e -= 0.5*m[k][j]*m[k][j]/rho[j];
			(void) printf("oned_PLM() calling godunovobl() for "
				      "tangential sweep\n");
			print_general_vector("coords = ",coords[j],dim,"\n");
			(void) printf("rho[%d] = %g, E[%d] = %g, e[%d] = %g, "
				      "Min_energy = %g\n",j,rho[j],j,E[j],j,e,
				      Min_energy(vst->state[j]));
		    }
		    godunovobl(dn,dt,tsten,tmpst,vmuscl->front);
		}
	        else if (sten != NULL)
		{
		    if (debugging("godunov"))
		    {
	                double e = E[j];
			for (k = 0; k < dim; ++k)
		            e -= 0.5*m[k][j]*m[k][j]/rho[j];
			(void) printf("oned_PLM() calling godunov() for "
				      "irregular cells sweep, j = %d, "
				      "idir = %d\n",j,iperm[swp_num]);
			print_general_vector("coords = ",coords[j],dim,"\n");
			(void) printf("rho[%d] = %g, E[%d] = %g, e[%d] = %g, "
				      "Min_energy = %g\n",j,rho[j],j,E[j],j,e,
				      Min_energy(vst->state[j]));
		    }
		    sten->reg_stencil = sten->prev_reg_stencil = NO;
		    godunov(dn,dt,tmpst,vmuscl->Q[0],swp_num,iperm,NULL,sten);
		}
		else
		{
		    static Stencil *tmpsten = NULL;
		    if (tmpsten == NULL)
		        tmpsten = alloc_stencil(3,vmuscl->front);
		    sten = tmpsten;
		    sten->fr = vmuscl->front;
		    sten->newfr = newfr;
		    sten->wave = vmuscl->wave;
		    sten->newwave = newwave;
		    for (k = 0; k < dim; ++k)
		    {
		        sten->icoords[-1][k] = icoords[k];
		        sten->icoords[ 0][k] = icoords[k];
		        sten->icoords[ 1][k] = icoords[k];
			Coords(sten->p[-1])[k] = coords[j-1][k];
			Coords(sten->p[ 0])[k] = coords[j][k];
			Coords(sten->p[ 1])[k] = coords[j+1][k];
		    }
		    sten->icoords[-1][vmuscl->idir] = j-offset-1;
		    sten->icoords[ 0][vmuscl->idir] = j-offset;
		    sten->icoords[ 1][vmuscl->idir] = j-offset+1;
		    sten->st[-1] = state[j-1];
		    sten->st[ 0] = state[j];
		    sten->st[ 1] = state[j+1];
		    sten->newcomp = component(coords[j],fr->interf);
		    sten->comp[-1] = sten->comp[0] = sten->comp[1] =
		        sten->newcomp;
		    sten->reg_stencil = sten->prev_reg_stencil = NO;
		    if (debugging("godunov"))
		    {
	                double e = E[j];
			for (k = 0; k < dim; ++k)
		            e -= 0.5*m[k][j]*m[k][j]/rho[j];
			(void) printf("oned_PLM() calling godunov() for "
				      "regular cells sweep, j = %d, "
				      "idir = %d\n",j,iperm[swp_num]);
			print_general_vector("coords = ",coords[j],dim,"\n");
			(void) printf("rho[%d] = %g, E[%d] = %g, e[%d] = %g, "
				      "Min_energy = %g\n",j,rho[j],j,E[j],j,e,
				      Min_energy(vst->state[j]));
#if defined(DEBUG_PLM)
			if (plm_interior_sten_reg(j,iperm[swp_num],sten))
			    (void) printf("Off front grid point\n");
			else
			    (void) printf("Near front grid point\n");
#endif /* defined(DEBUG_PLM) */
		    }
		    godunov(dn,dt,tmpst,vmuscl->Q[0],swp_num,iperm,NULL,sten);
		    sten = NULL;
		}
	        rho[j] = Dens(tmpst);
	        for (k = 0; k < dim; ++k)
		    m[k][j] = Mom(tmpst)[k];
		E[j] = Energy(tmpst);
	    }
	}

#if defined(DEBUG_PLM)
	if (debugging("plm"))
	{
	    static double rho_max, E_max, m_max[3], v_max[3], e_max;
	    static double rho_min, E_min, m_min[3], v_min[3], e_min;
	    static double ke_max, ke_min;
	    static int last_step = 0;
	    static FILE *maxfile = NULL;
	    static FILE *minfile = NULL;
	    double ke, e, v[3];

	    if (maxfile == NULL)
	    {
	        last_step = vmuscl->front->step;
		maxfile = fopen("PLM_MAX_VALUES","w");
		print_machine_parameters(maxfile);
		minfile = fopen("PLM_MIN_VALUES","w");
		print_machine_parameters(minfile);
	        (void) fprintf(maxfile,"%-6s %-14s %-14s %-14s %-14s",
		              "step","rho_max","E_max","ke_max","e_max");
	        (void) fprintf(minfile,"%-6s %-14s %-14s %-14s %-14s",
		              "step","rho_min","E_min","ke_min","e_min");
	        for (k = 0; k < dim; ++k)
		{
		    (void) fprintf(maxfile," m[%d]%-10s",k,"_max");
		    (void) fprintf(minfile," m[%d]%-10s",k,"_min");
		}
	        for (k = 0; k < dim; ++k)
		{
		    (void) fprintf(maxfile," v[%d]%-10s",k,"_max");
		    (void) fprintf(minfile," v[%d]%-10s",k,"_min");
		}
		(void) fprintf(maxfile,"\n");
		(void) fprintf(minfile,"\n");
	        rho_max = -HUGE_VAL; rho_min = HUGE_VAL;
		  E_max = -HUGE_VAL;   E_min = HUGE_VAL;
		  e_max = -HUGE_VAL;   e_min = HUGE_VAL;
		 ke_max = -HUGE_VAL;  ke_min = HUGE_VAL;
	        for (k = 0; k < dim; ++k)
		{
		    m_max[k] = v_max[k] = -HUGE_VAL;
		    m_min[k] = v_min[k] =  HUGE_VAL;
		}
	    }
	    for (j = start; j < end; ++j)
	    {
	        rho_max = max(rho[j],rho_max); rho_min = min(rho[j],rho_min);
	          E_max = max(E[j],E_max);       E_min = min(E[j],E_min);
	        for (ke = 0.0, k = 0; k < dim; ++k)
		{
		    m_max[k] = max(m[k][j],m_max[k]);
		    m_min[k] = min(m[k][j],m_min[k]);
		    v[k] = m[k][j]/rho[j];
		    v_max[k] = max(v[k],v_max[k]);
		    v_min[k] = min(v[k],v_min[k]);
		    ke += 0.5*m[k][j]*v[k];
		}
		e = (E[j]-ke)/rho[j];
		 e_max = max(e,e_max);    e_min = min(e,e_min);
		ke_max = max(ke,ke_max); ke_min = min(ke,ke_min);
	    }
	    if (vmuscl->front->step > last_step)
	    {
		(void) fprintf(maxfile,"%-6d %-14g %-14g %-14g %-14g",
		               last_step,rho_max,E_max,ke_max,e_max);
		(void) fprintf(minfile,"%-6d %-14g %-14g %-14g %-14g",
		               last_step,rho_min,E_min,ke_min,e_min);
	        for (k = 0; k < dim; ++k)
		{
		    (void) fprintf(maxfile," %-14g ",m_max[k]);
		    (void) fprintf(minfile," %-14g ",m_min[k]);
		}
	        for (k = 0; k < dim; ++k)
		{
		    (void) fprintf(maxfile," %-14g ",v_max[k]);
		    (void) fprintf(minfile," %-14g ",v_min[k]);
		}
		(void) fprintf(maxfile,"\n");
		(void) fprintf(minfile,"\n");
	        last_step = vmuscl->front->step;
	        rho_max = -HUGE_VAL; rho_min = HUGE_VAL;
		  E_max = -HUGE_VAL;   E_min = HUGE_VAL;
		  e_max = -HUGE_VAL;   e_min = HUGE_VAL;
	        for (k = 0; k < dim; ++k)
		{
		    m_max[k] = v_max[k] = -HUGE_VAL;
		    m_min[k] = v_min[k] =  HUGE_VAL;
		}
	    }
	}
#endif /* defined(DEBUG_PLM) */

	Vec_Gas_field_set(vst,c) = NO;
	Vec_Gas_field_set(vst,e) = NO;
	load_sound_speed(vst,start+offset,end-start);
	c = vst->c + offset;
	v0 = vst->v[0] + offset;
	if (vmuscl->wave != NULL)
	{
	    for (j = start; j < end; ++j)
	    {
	        v0[j] = m[0][j]/rho[j];
	        if (Params(state[j]) != NULL)
		{
		    /*speed = 2.0*Params(state[j])->avisc.sp_coef*
		            (fabs(v0[j])+c[j]);OLD VERSION*/
		    speed = Params(state[j])->avisc.sp_coef*
		            (fabs(v0[j])+c[j]);
		    set_max_wave_speed(vmuscl->idir,speed,state[j],
		                       coords[j],wave);
		}
	    }
	}
	else
	{
	    const double* const *Q = vst->Q;
	    for (j = start; j < end; ++j)
	    {
	        v0[j] = m[0][j]/rho[j];
	        if (Params(state[j]) != NULL)
		{
		    int k;
		    /*speed = 2.0*Params(state[j])->avisc.sp_coef*
		            (fabs(v0[j])+c[j]);OLD VERSION*/
		    speed = Params(state[j])->avisc.sp_coef*
		            (fabs(v0[j])+c[j]);
		    for (k = 0; k < dim; ++k)
		    {
		        set_max_front_speed(k,Q[0][k]*speed,
			                    state[j],coords[j],fr);
		    }
		    set_max_front_speed(dim,speed/dn,state[j],coords[j],fr);
		}
	    }
	}
#if defined(DEBUG_PLM)
	debug_print("plm","Left oned_PLM()\n");
#endif /* defined(DEBUG_PLM) */
}		/*end oned_PLM*/


EXPORT	void	set_plm_default_opts(
	Muscl_Opts  *mopts,
	Hyp_method  *method,
	int         dim)
{
	MUSCL_PromptType_Reconstructor        *Sintrp;
	MUSCL_PromptType_Rsolver              *Rsolver;
	MUSCL_PromptType_characteristic_solve *Moc;

	set_muscl_default_hooks(mopts,method,dim);

	mopts->_reconstructor = g_plm_reconstructor;
	mopts->_load_state_data = g_load_plm_state_data;
	mopts->_flux_vectors = g_plm_flux_vectors;
	mopts->_rsolver = g_muscl_exact_rsolver;
	mopts->_rmidstate = g_exact_Riemann_midstate;
	mopts->worksp_len = 0;
	mopts->_half_step = NULL;
	mopts->_strong_wave_half_step = NULL;
	mopts->_cons_src = NULL;
	mopts->_add_art_visc1 = g_plm_add_art_visc1;
	mopts->_add_art_visc2 = g_plm_add_art_visc2;

	mopts->_characteristic_solve = g_riemann_characteristic_solve;
	mopts->_compute_eigens = g_plm_compute_eigens;

	mopts->monotone_reconstruction = YES;
	mopts->link_reconstructions = YES;

	uni_array(&Sintrp,2,sizeof(MUSCL_PromptType_Reconstructor));
	mopts->Sintrp = Sintrp;
	Sintrp[0].prompt = "Colella Piecewise Linear Method";
	Sintrp[0].select = "cplm";
	Sintrp[0].reconstructor = mopts->_reconstructor;

	uni_array(&Rsolver,5,sizeof(MUSCL_PromptType_Rsolver));
	mopts->Rsolver = Rsolver;
	Rsolver[0].prompt = "Exact Riemann solver";
	Rsolver[0].select = "e";
	Rsolver[0].rsolver = g_muscl_exact_rsolver;
	Rsolver[0].rmidstate = g_exact_Riemann_midstate;
	Rsolver[1].prompt = "Colella-Glaz's approximate Riemann solver";
	Rsolver[1].select = "c";
	Rsolver[1].rsolver = cg_rsolve;
	Rsolver[1].rmidstate = g_cg_Riemann_midstate;
	Rsolver[2].prompt = "Linear US/UP fit (Dukowicz)";
	Rsolver[2].select = "d";
	Rsolver[2].rsolver = g_linear_us_up_rsolver;
	Rsolver[2].rmidstate = g_linear_us_up_Riemann_midstate;
	Rsolver[3].prompt = "Gamma Law fit";
	Rsolver[3].select = "g";
	Rsolver[3].rsolver = g_gamma_law_fit_rsolver;
	Rsolver[3].rmidstate = g_gamma_law_fit_Riemann_midstate;

	uni_array(&Moc,4,sizeof(MUSCL_PromptType_characteristic_solve));
	mopts->Moc = Moc;
	Moc[0].prompt = "Riemann problem based characteristic solve";
	Moc[0].select = "r";
	Moc[0].characteristic_solve = g_riemann_characteristic_solve;
	Moc[1].prompt = "Implicit in time characteristic solve";
	Moc[1].select = "i";
	Moc[1].characteristic_solve = g_implicit_characteristic_solve;
	Moc[2].prompt = "First order direct characteristic solve";
	Moc[2].select = "f";
	Moc[2].characteristic_solve = g_first_order_direct_characteristic_solve;

}		/*end set_plm_muscl_opts*/

/*ARGSUSED*/
LOCAL	Vec_Muscl *g_load_plm_state_data(
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
	RECT_GRID        *gr = fr->rect_grid;
	double            **CellEdge;
	double            **coords = vst->coords + offset;
	double            **u, **ucon;
	double            alpha = 0.0;
	const double* const *Q = vst->Q;
	int              i, j, iswp;
	static Vec_Muscl *vmuscl = NULL;

#if defined(DEBUG_PLM)
	debug_print("loadplm","Entered g_load_plm_state_data\n");
#endif /* defined(DEBUG_PLM) */

	if (vmuscl == NULL || vsize > vmuscl->max_vsize || dim > vmuscl->dim)
	{
	    vmuscl = g_muscl_free_wk_space(vmuscl);
	    vmuscl = plm_alloc_wk_space(mopts,vsize,dim,vmuscl);
	}

	/* Set geometry */
	vmuscl->coords = coords;
	CellEdge = vmuscl->CellEdge;
	if (tsten == NULL)
	{
	    iswp = iperm[swp_num];
	    for (j = 0; j < dim; ++j)
	    {
	        if (j == iswp)
		{
	            for (i = 0; i < vsize; ++i)
		        CellEdge[i][j] = coords[i][j] - 0.5*dn;
		    CellEdge[i][j] = coords[i-1][j] + 0.5*dn;
		}
		else
		{
	            for (i = 0; i < vsize; ++i)
		        CellEdge[i][j] = coords[i][j];
		    CellEdge[i][j] = coords[i-1][j];
	        }
	    }
	}
	else
	{
	    const double *dir = Q[0];
	    for (j = 0; j < dim; ++j)
	    {
	        for (i = 0; i < vsize; ++i)
		    CellEdge[i][j] = coords[i][j] - 0.5*dn*dir[j];
		CellEdge[i][j] = coords[i-1][j] + 0.5*dn*dir[j];
	    }
	}
	alpha = rotational_symmetry();
	if ((alpha == 0.0) || ((tsten == NULL) && iperm[swp_num] != 0))
	{
	    vmuscl->alpha = 0.0;
	    for (i = 0; i < vsize; ++i)
	    {
		vmuscl->A[i] = 1.0;
		vmuscl->dV[i] = dn;
	    }
	    vmuscl->A[i] = 1.0;
	}
	else
	{
	    if (alpha == 1.0)
	    {
	        vmuscl->alpha = (tsten == NULL) ? 1.0 : Q[0][0];
	        for (i = 0; i < vsize; ++i)
		{
		    vmuscl->A[i] = pos_radius(CellEdge[i][0],gr);
		    vmuscl->dV[i] = pos_radius(coords[i][0],gr)*dn;
		}
		vmuscl->A[i] = pos_radius(CellEdge[i][0],gr);
	    }
	    else if (alpha == 2.0)
	    {
	        double r;
	        vmuscl->alpha = 2.0;
	        for (i = 0; i < vsize; ++i)
		{
		    r = pos_radius(CellEdge[i][0],gr);
		    vmuscl->A[i] = r*r;
		    r = pos_radius(coords[i][0],gr);
		    vmuscl->dV[i] = (r*r+dn*dn/12)*dn;
		}
		r = pos_radius(CellEdge[i][0],gr);
		vmuscl->A[i] = r*r;
	    }
	    else
	    {
		screen("ERROR in g_load_plm_state_data(), "
		       "invalid value for alpha = %g\n",alpha);
		clean_up(ERROR);
	    }
	}

	vmuscl->sizest = fr->sizest;
	vmuscl->idir = (iperm != NULL) ? iperm[swp_num] : dim;
	vmuscl->Q = Q;
	vmuscl->dt = dt;
	vmuscl->dn = dn;
	vmuscl->sten = sten;
	vmuscl->tsten = tsten;
	vmuscl->front = fr;
	vmuscl->wave = wave;
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
	vmuscl->src = NULL;
	vmuscl->offset = offset;	vmuscl->vsize = vsize;

	clear_Vec_Gas_set_flags(&vmuscl->VMst);
	clear_Vec_Gas_set_flags(&vmuscl->VLst);
	clear_Vec_Gas_set_flags(&vmuscl->VRst);
	clear_Vec_Gas_set_flags(&vmuscl->VM0st);
	clear_Vec_Gas_set_flags(&vmuscl->VM1st);
	vmuscl->VMst.Q = vst->Q;
	vmuscl->VLst.Q = vst->Q;
	vmuscl->VRst.Q = vst->Q;
	vmuscl->VM0st.Q = vst->Q;
	vmuscl->VM1st.Q = vst->Q;

	g_load_VGas_state_vectors(offset,vsize,vst,dim);

	/* Set params jump points for vec gas states */
	copy_vec_state_params(&vmuscl->VMst,vst,offset,vsize);
	copy_vec_state_params(&vmuscl->VRst,vst,offset,vsize);
	copy_vec_state_params(&vmuscl->VLst,vst,offset-1,vsize);
	copy_vec_state_params(&vmuscl->VM0st,vst,offset,vsize);
	copy_vec_state_params(&vmuscl->VM1st,vst,offset,vsize);

	/*put the state variables and source terms in a form more
	  easily usable by the subroutines*/

	u = vmuscl->u;
	ucon = vmuscl->ucon;
	u[vmuscl->index.density] = vst->rho + offset;
	ucon[vmuscl->index.density] = vst->rho + offset;
	u[vmuscl->index.energy] = vst->e + offset;
	u[vmuscl->index.internal_energy_density] = vst->re + offset;
	ucon[vmuscl->index.energy] = vst->en_den + offset;
	for (i = 0; i < dim; ++i)    
	{
	    u[vmuscl->index.v[i]] = vst->v[i] + offset;
	    ucon[vmuscl->index.v[i]] = vst->m[i] + offset;
	}
	u[vmuscl->index.pressure] = vst->p + offset;
	u[vmuscl->index.sound_speed] = vst->c + offset;
	u[vmuscl->index.GAM] = vst->GAM + offset;
	if (vmuscl->index.FD >= 0)
	    u[vmuscl->index.FD] = vst->FD + offset;

	u[vmuscl->index.grav] = vmuscl->grav;

#if defined(DEBUG_PLM)
	if (debugging("loadplm"))
	{
	    double *A = vmuscl->A;
	    double *dV = vmuscl->dV;
	    (void) printf("\n%-5s %-14s %-14s\n","index","A","dV");
	    for (i = 0; i < vsize; ++i)
	        (void) printf("%-5d %-14g %-14g\n",i,A[i],dV[i]);
	    (void) printf("\n\n");
	}

	debug_print("loadplm","Left g_load_plm_state_data\n");
#endif /* defined(DEBUG_PLM) */
	return vmuscl;
}		/*end g_load_plm_state_data*/

/*
*			g_plm_reconstructor():
*
*	Computes the MUSCL slopes as defined in 
*	"A Direct Eulerian MUSCL Scheme for Gas Dynamics"
*	P. Colella, SIAM J. Sci. Stat. Comput., Vol. 6, No. 1 (1985).
*
*	The scheme uses a seven point stencil centered at each
*	grid point to compute the slopes.
*	
*/

LOCAL void g_plm_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int	        dim;
	int	        j, k;	/* j is the mesh index */
	int	        kmax;
	double	        *uk, *duk;
	double           **uL, **uR, *uLk, *uRk;
	double           **u, **du;
	Vec_Muscl_Index index = vmuscl->index;

#if defined(DEBUG_PLM)
	debug_print("plmrecon","Entered g_plm_reconstructor()\n");
#endif /* defined(DEBUG_PLM) */

	if (tmpst == NULL)
	    alloc_state(vmuscl->front->interf,&tmpst,vmuscl->front->sizest);

	dim =	   vmuscl->dim;
	kmax =	   dim+2;

	u  = vmuscl->u;
	du  = vmuscl->du;

#if defined(DEBUG_PLM)
	if (debugging("plmrecon"))
	    (void) printf("nvar_u = %d\n",vmuscl->nvar_u);
#endif /* defined(DEBUG_PLM) */

	reconstruct_flow_field(u[index.density],du[index.density],start,end);
	reconstruct_flow_field(u[index.energy],du[index.energy],start,end);
	reconstruct_flow_field(u[index.v[0]],du[index.v[0]],start,end);
	for (k = 0; k < dim; ++k)
	    reconstruct_flow_field(u[index.v[k]],du[index.v[k]],start,end);
	if (vmuscl->msf != NULL)
	{
	    double *chi = vmuscl->msf->chi;
	    compute_slope_limiting_coeffs(start,end,vmuscl);
	    modify_flow_field_slope(du[index.density],chi,start,end);
	    modify_flow_field_slope(du[index.energy],chi,start,end);
	    for (k = 0; k < dim; ++k)
	        modify_flow_field_slope(du[index.v[k]],chi,start,end);
#if defined(DEBUG_PLM)
	    if (debugging("plmrecon"))
	    {
	        (void) printf("Slope limiting coefficients\n");
	        (void) printf("%-5s %-14s\n","index","chi");
		for (j = start; j < end; ++j)
		    (void) printf("%-5d %-14g\n",j,chi[j]);
	    }
#endif /* defined(DEBUG_PLM) */
	}
	if (vmuscl->link_reconstructions)
	{
	    int l;
	    for (j = start; j < end; ++j)
	    {
	        for (k = 0; k < kmax; ++k)
		{
		    if (du[k][j] == 0.0)
		    {
	                for (l = 0; l < kmax; ++l)
			    du[l][j] = 0.0;
			break;
		    }
		}
	    }
	}
	if (vmuscl->monotone_reconstruction)
	{
	    monotonize_reconstruction(u[index.density],du[index.density],
	                              start,end-1);
	    monotonize_reconstruction( u[index.energy],du[index.energy],
	                              start,end-1);
	    for (k = 0; k < dim; ++k)
	        monotonize_reconstruction(u[index.v[k]],du[index.v[k]],
		                          start,end-1);
	}
	if ((vmuscl->alpha != 0) && (vmuscl->idir == 0))
	{
	    double **coords = vmuscl->coords;
	    int j0;

	    for (j0 = start; j0 < end; ++j0)
	        if (coords[j0][0] > 0.0)
		    break;
	    for (j = start; j < j0; ++j)
	    {
	        for (k = 0; k < kmax; ++k)
		{
		    du[k][j] = -du[k][2*j0 -1 - j];
		}
		du[index.v[0]][j] = du[index.v[0]][2*j0 -1 - j];
	    }
	}

	uL = vmuscl->uL;
	uR = vmuscl->uR;
	for (k = 0; k < kmax; ++k)
	{
	    uLk = uL[k];
	    uRk = uR[k];
	    uk = u[k];
	    duk = du[k];
	    for (j = start; j < end; ++j)
	    {
	        uRk[j]   = uk[j] - 0.5*duk[j];
	        uLk[j+1] = uk[j] + 0.5*duk[j];
	    }
	}
	for (j = start; j < end; ++j)
	{
	    if (vmuscl->VLst.rho[j] < 0.0)
	    {
	        Locstate st0 = vmuscl->vst->state[j-1+vmuscl->offset];
		double    min_p = Min_pressure(st0);

		state_on_adiabat_with_pr(st0,min_p,tmpst,EGAS_STATE);
	        vmuscl->VLst.rho[j] = Dens(tmpst);
	        vmuscl->VLst.e[j] = Energy(tmpst);
		vmuscl->VLst.v[0][j] = u[2][j-1]+riemann_wave_curve(st0,min_p);
	    }
	    if (vmuscl->VRst.rho[j] < 0.0)
	    {
	        Locstate st0 = vmuscl->vst->state[j+vmuscl->offset];
		double    min_p = Min_pressure(st0);

		state_on_adiabat_with_pr(st0,min_p,tmpst,EGAS_STATE);
	        vmuscl->VRst.rho[j] = Dens(tmpst);
	        vmuscl->VRst.e[j] = Energy(tmpst);
		vmuscl->VRst.v[0][j] = u[2][j]-riemann_wave_curve(st0,min_p);
	    }
	}
	if ((vmuscl->alpha != 0) && (vmuscl->idir == 0))
	{
	    double **coords = vmuscl->coords;
	    int j0;

	    for (j0 = start; j0 < end; ++j0)
	        if (coords[j0][0] > 0.0)
		    break;
	    if (j0 > start)
	    {
		vmuscl->VRst.v[0][j0] = 0.0;
		vmuscl->VLst.v[0][j0] = 0.0;
	    }
	}
	Vec_Gas_field_set(&vmuscl->VLst,rho) = YES;
	Vec_Gas_field_set(&vmuscl->VLst,e) = YES;
	Vec_Gas_field_set(&vmuscl->VLst,v) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,rho) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,e) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,v) = YES;
	load_internal_energy_density(&vmuscl->VLst,start+1,end-start);
	load_internal_energy_density(&vmuscl->VRst,start,end-start);

	/* Compute slopes for pressure and sound speed */
	load_pressure_and_gammas(&vmuscl->VLst,start+1,end-start);
	load_pressure_and_gammas(&vmuscl->VRst,start,end-start);
	for (j = start; j < end; ++j)
	{
	   vmuscl->VLst.c[j+1] = sqrt(vmuscl->VLst.c2[j+1]); 
	   vmuscl->VRst.c[j] = sqrt(vmuscl->VRst.c2[j]); 
	}
	Vec_Gas_field_set(&vmuscl->VLst,c) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,c) = YES;

	if (is_gravity())
	{
	    double t0 = vmuscl->front->time;
	    double t1 = vmuscl->front->time + vmuscl->dt;
	    double *gL = vmuscl->gL;
	    double *gR = vmuscl->gR;
	    double *gM0 = vmuscl->gM0;
	    double *gM1 = vmuscl->gM1;
	    double **CE = vmuscl->CellEdge;
	    int   j;

	    if (vmuscl->tsten == NULL)
	    {
	        int idir = vmuscl->idir;
	        for (j = start; j <= end; ++j)
		{
	            gM0[j] = gL[j] = gR[j] = gravity(CE[j],t0)[idir];
	            gM1[j] = gravity(CE[j],t1)[idir];
		}
	    }
	    else
	    {
	        const double *dir = vmuscl->Q[0];
	    	for (j = start; j <= end; ++j)
		{
	    	    gM0[j] = gL[j] = gR[j] =
		        scalar_product(dir,gravity(CE[j],t0),dim);
	    	    gM1[j] = scalar_product(dir,gravity(CE[j],t1),dim);
		}
	    }
	}
	else
	{
	    double *gL = vmuscl->gL;
	    double *gR = vmuscl->gR;
	    double *gM0 = vmuscl->gM0;
	    double *gM1 = vmuscl->gM1;
	    for (j = start; j <= end; ++j)
	    	gM0[j] = gM1[j] = gL[j] = gR[j] = 0.0;
	}

	set_derived_field_slope_and_mid(uL[index.internal_energy_density],
	                                uR[index.internal_energy_density],
					u[index.internal_energy_density],
					du[index.internal_energy_density],
					start,end);
	set_derived_field_slope_and_mid(uL[index.pressure],uR[index.pressure],
					u[index.pressure],du[index.pressure],
					start,end);
	set_derived_field_slope_and_mid(uL[index.sound_speed],
	                                uR[index.sound_speed],
					u[index.sound_speed],
					du[index.sound_speed],
					start,end);
	set_derived_field_slope_and_mid(uL[index.GAM],uR[index.GAM],
					u[index.GAM],du[index.GAM],
					start,end);
	if (index.FD >= 0)
	    set_derived_field_slope_and_mid(uL[index.FD],uR[index.FD],
					    u[index.FD],du[index.FD],
					    start,end);
	set_derived_field_slope_and_mid(uL[index.grav],uR[index.grav],
					u[index.grav],du[index.grav],
					start,end);
	Vec_Gas_field_set(vmuscl->vst,p) = YES;
	Vec_Gas_field_set(vmuscl->vst,c) = YES;
	Vec_Gas_field_set(vmuscl->vst,GAM) = YES;
	if (index.FD >= 0)
	    Vec_Gas_field_set(vmuscl->vst,FD) = YES;


#if defined(DEBUG_PLM)
	if (debugging("plmrecon"))
	{
	    double *rho, *e, **v, *p, *c;
	    double *drho, *de, **dv, *dp, *dc;
	    double *rhol, *el, **vl, *pl, *cl;
	    double *rhor, *er, **vr, *pr, *cr;
	    double **CE = vmuscl->CellEdge;
	    double **coords = vmuscl->coords;
	    double *gL = vmuscl->gL;
	    double *gR = vmuscl->gR;
	    double *gM0 = vmuscl->gM0;
	    double *gM1 = vmuscl->gM1;
	    const char *meanhdr, *reconhdr, *slphdr;

	    switch (MusclSweepType(vmuscl))
	    {
	    case TANGENTIAL_SWEEP:
	        meanhdr = "TANGENTIAL_SWEEP_MEAN_STATES";
	        reconhdr = "TANGENTIAL_SWEEP_RECONSTRUCTED_STATES";
	        slphdr = "TANGENTIAL_SWEEP_SLOPES";
		break;
	    case IRREGULAR_INTERIOR_SWEEP:
	        meanhdr = "IRREGULAR_INTERIOR_SWEEP_MEAN_STATES";
	        reconhdr = "IRREGULAR_INTERIOR_SWEEP_RECONSTRUCTED_STATES";
	        slphdr = "IRREGULAR_INTERIOR_SWEEP_SLOPES";
		break;
	    case REGULAR_INTERIOR_SWEEP:
	        meanhdr = "REGULAR_INTERIOR_SWEEP_MEAN_STATES";
	        reconhdr = "REGULAR_INTERIOR_SWEEP_RECONSTRUCTED_STATES";
	        slphdr = "REGULAR_INTERIOR_SWEEP_SLOPES";
		break;
	    case UNSET_SWEEP_TYPE:
	    default:
	        meanhdr = "UNSET_SWEEP_TYPE_MEAN_STATES";
	        reconhdr = "UNSET_SWEEP_TYPE_RECONSTRUCTED_STATES";
	        slphdr = "UNSET_SWEEP_TYPE_SLOPES";
	    }

	    rhol = uL[index.density];
	    el = uL[index.energy];
	    vl = uL + index.v[0];
	    pl = uL[index.pressure];
	    cl = uL[index.sound_speed];

	    rho  = u[index.density];
	    e  = u[index.energy];
	    v = u + index.v[0];
	    p  = u[index.pressure];
	    c  = u[index.sound_speed];

	    rhor = uR[index.density];
	    er = uR[index.energy];
	    vr = uR + index.v[0];
	    pr = uR[index.pressure];
	    cr = uR[index.sound_speed];

	    print_data_header(dim,"rho","e","v","p","c",meanhdr);
	    j = start;
	    for (j = start-2; j < (end+2); ++j)
	    {
	        print_data_line(dim,j,CE,rho,e,v,p,c);
		print_data_line(dim,j,CE+1,rho,e,v,p,c);
	    }
	    (void) printf("\n\n");

	    print_data_header(dim,"rho","e","v","p","c",reconhdr);
	    for (j = start; j < end; ++j)
	    {
		print_data_line(dim,j,CE,rhor,er,vr,pr,cr);
		print_data_line(dim,j+1,CE,rhol,el,vl,pl,cl);
	    }
	    (void) printf("\n\n");

	    drho  = du[index.density];
	    de  = du[index.energy];
	    dv = du + index.v[0];
	    dp  = du[index.pressure];
	    dc  = du[index.sound_speed];

	    print_data_header(dim,"drho","de","dv","dp","dc",slphdr);
	    for (j = start; j < end; ++j)
		print_data_line(dim,j,coords,drho,de,dv,dp,dc);
	    (void) printf("\n\n");

	    (void) output();
	    (void) printf("%-5s %-14s %-14s G0\n","index","x0","g");
	    for (j = start; j <= end; ++j)
	    {
	        (void) printf("%-5d %-14g %-14g\n",j,CE[j][0],gL[j]);
	        (void) printf("%-5d %-14g %-14g\n",j,CE[j][0],gR[j]);
	        (void) printf("%-5d %-14g %-14g\n",j,coords[j][0],gM0[j]);
	    }
	    (void) printf("\n\n");
	    (void) printf("%-5s %-14s %-14s G1\n","index","x0","g");
	    for (j = start; j <= end; ++j)
	        (void) printf("%-5d %-14g %-14g\n",j,coords[j][0],gM1[j]);
	    (void) printf("\n\n");
	}
	debug_print("plmrecon","Left g_plm_reconstructor()\n");
#endif /* defined(DEBUG_PLM) */
}		/*end g_plm_reconstructor*/

LOCAL	 void	reconstruct_flow_field(
	const double *u,
	double       *du,
	int         start,
	int         end)
{
	int    j;
	static double *b = NULL,     *bstore = NULL;
	static double *c = NULL,     *cstore = NULL;
	static double *f = NULL,     *fstore = NULL;
	static double *dulim = NULL, *dulimstore = NULL;
	static double *duf = NULL,   *dufstore = NULL;
	static double *sgn = NULL, *sgnstore = NULL;
	static int   vsize = 0;

	if (vsize < (end-start+2))
	{
	    if (bstore != NULL)
	    {
	        free_these(6,bstore,cstore,fstore,
		             dulimstore,dufstore,sgnstore);
	    }
	    vsize = end-start+2;
	    uni_array(&bstore,vsize,FLOAT);
	    uni_array(&cstore,vsize,FLOAT);
	    uni_array(&fstore,vsize,FLOAT);
	    uni_array(&dulimstore,vsize,FLOAT);
	    uni_array(&dufstore,vsize,FLOAT);
	    uni_array(&sgnstore,vsize,FLOAT);
	    b     = bstore     - start + 1;
	    c     = cstore     - start + 1;
	    f     = fstore     - start + 1;
	    dulim = dulimstore - start + 1;
	    duf   = dufstore - start + 1;
	    sgn   = sgnstore - start + 1;
	}

	for (j = start-1; j < end+1; ++j)
	{
	    b[j] = 2.0 * (u[j  ] - u[j-1]);
	    c[j] = 0.5 * (u[j+1] - u[j-1]);
	    f[j] = 2.0 * (u[j+1] - u[j  ]);
	}
	for (j = start-1; j < end+1; ++j)
	{
	    sgn[j] = (c[j] >= 0.0) ? 1.0 : -1.0;
	    if (b[j]*f[j] > 0.0)
		dulim[j] = (fabs(b[j])>fabs(f[j])) ? fabs(f[j]):fabs(b[j]);
	    else
		dulim[j] = 0.0;
	    duf[j] = (dulim[j] > fabs(c[j])) ? c[j] : dulim[j]*sgn[j];
	}
	for (j = start; j < end; ++j)
	{
	    du[j] = (2.0/3.0)*fabs((u[j+1] - 0.25*duf[j+1]) -
		                   (u[j-1] + 0.25*duf[j-1]));
	    if (du[j] > dulim[j])
		du[j] = dulim[j];
	    du[j] *= sgn[j];
	}
}		/*end reconstruct_flow_field*/

LOCAL	void	monotonize_reconstruction(
	const double *u,
	double       *du,
	int         start,
	int         end)
{
	double ul, ur, um, d;
	int   j;

	for (j = start; j < end; ++j)
	{
	    if (du[j]*du[j+1] > 0)
	    {
	        ul = u[j  ]+0.5*du[j  ];
		ur = u[j+1]-0.5*du[j+1];
		if ((ur-ul)*du[j] < 0.0)
		{
		    um = 0.5*(ul+ur);
		    d = 2.0*(um - u[j]);
		    if (d*du[j] < 0.0)
		        du[j] = 0.0;
		    else if (fabs(d) < fabs(du[j]))
		        du[j] = d;
		    d = 2.0*(u[j+1]-um);
		    if (d*du[j+1] < 0.0)
		        du[j+1] = 0.0;
		    else if (fabs(d) < fabs(du[j+1]))
		        du[j+1] = d;
		}
	    }
	}
}		/*end monotonize_reconstruction*/

LOCAL	void	modify_flow_field_slope(
	double       *duk,
	const double *chi,
	int         start,
	int         end)
{
	int j;

	for (j = start; j < end; ++j)
	    duk[j] *= chi[j];
}		/*end modify_flow_field_slope*/

LOCAL	void set_derived_field_slope_and_mid(
	const double *uLk,
	const double *uRk,
	double       *uk,
	double       *duk,
	int         start,
	int         end)
{
	int j;

	for (j = start; j < end; ++j)
	{
	    duk[j] =      uLk[j+1] - uRk[j];
	     uk[j] = 0.5*(uLk[j+1] + uRk[j]);
	}
}		/*end set_derived_field_slope_and_mid*/

/*ARGSUSED*/
LOCAL void g_plm_flux_vectors(
	int	   start,
	int	   end,
	double	   **u,
	double	   *p,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	double *rho, *e, **v, *v0;
	double *rhoF, *EF, **mF, *mF0;
	double *rhoH, *EH, **mH, *mH0;
	double **F, **H;
	double magv2;
	int   j, k, dim = vmuscl->dim;
	int   startj, endj;

#if defined(DEBUG_PLM)
	debug_print("mflux","Entered g_plm_flux_vectors()\n");
#endif /* defined(DEBUG_PLM) */

	if (Flux == NULL)
	{
#if defined(DEBUG_PLM)
	    debug_print("mflux","Flux == NULL\nLeft g_plm_flux_vectors()\n");
#endif /* defined(DEBUG_PLM) */
	    return;
	}

	rho = u[vmuscl->index.density];
	e = u[vmuscl->index.energy];
	v = u + vmuscl->index.v[0];
	v0 = v[0];

	F = Flux->F;
	H = Flux->H;
	rhoF = F[vmuscl->index.density];
	rhoH = H[vmuscl->index.density];
	EF = F[vmuscl->index.energy];
	EH = H[vmuscl->index.energy];
	mF = F + vmuscl->index.v[0];
	mF0 = mF[0];
	mH = H + vmuscl->index.v[0];
	mH0 = mH[0];

	for (j = start; j < end; ++j)
	{
	    rhoH[j] = 0.0;
	    mH0[j] = p[j];
	    for (k = 1; k < dim; ++k)
		mH[k][j] = 0.0;
	    EH[j] = 0.0;
	}
#if defined(DEBUG_PLM)
	if (debugging("mflux"))
	{
	    for (j = start; j < end; ++j)
		(void) printf("\tp[%d] = %g\n",j,p[j]);
	}
#endif /* defined(DEBUG_PLM) */
	startj = start;
	endj = end;
	if (vmuscl->alpha != 0.0)
	{
	    double **coords = vmuscl->coords;
	    double nu = sqrt(sqr(vmuscl->Q[0][1])+sqr(vmuscl->Q[0][2]));
	    double nv0;

	    /* For cells crossing the r = 0 axis,  enforce the
	     * boundary condition that the radial component of
	     * velocity at the cell edge r = 0 vanishes
	     */

	    for (j = start; coords[j-1][0]*coords[j][0] < 0.0 && j < end; ++j)
	    {
	        nv0 = nu*v0[j];
	        /* mass */
	        rhoF[j] = rho[j]*nv0;


	        /* Momentum */
	        mF0[j] = rhoF[j]*nv0;

	        magv2 = nv0*nv0;
	        for (k = 1; k < dim; ++k)
	        {
	            magv2 += v[k][j]*v[k][j];

		    mF[k][j] = rhoF[j]*v[k][j];
	        }

	        /* Energy */
	        EF[j] = rhoF[j]*(0.5*magv2 + e[j]) + nv0*p[j];
	    }
	    startj = j;
	    for (j = end-1; coords[j-1][0]*coords[j][0]<0.0 && j>startj; --j)
	    {
	        nv0 = nu*v0[j];
	        /* mass */
	        rhoF[j] = rho[j]*nv0;


	        /* Momentum */
	        mF0[j] = rhoF[j]*nv0;

	        magv2 = nv0*nv0;
	        for (k = 1; k < dim; ++k)
	        {
	            magv2 += v[k][j]*v[k][j];

		    mF[k][j] = rhoF[j]*v[k][j];
	        }

	        /* Energy */
	        EF[j] = rhoF[j]*(0.5*magv2 + e[j]) + nv0*p[j];
	    }
	    endj = j+1;
	}
	for (j = startj; j < endj; ++j)
	{
	    /* mass */
	    rhoF[j] = rho[j]*v0[j];


	    /* Momentum */
	    mF0[j] = rhoF[j]*v0[j];

	    magv2 = v0[j]*v0[j];
	    for (k = 1; k < dim; ++k)
	    {
	        magv2 += v[k][j]*v[k][j];

		mF[k][j] = rhoF[j]*v[k][j];
	    }

	    /* Energy */
	    EF[j] = rhoF[j]*(0.5*magv2 + e[j]) + v0[j]*p[j];
	}
#if defined(DEBUG_PLM)
	debug_print("mflux","Left g_plm_flux_vectors()\n");
#endif /* defined(DEBUG_PLM) */
}		/*end g_plm_flux_vectors*/

LOCAL	void	characteristic_traceback(
	int       start,
	int       end,
	Vec_Muscl *vmuscl)
{
	double      d;
	double      *vM0, *cM0;
	double      *v, *c, *dv, *dc;
	double      **uL, **uM, **uR, **u, **du, **uM0;
	double      **Lcrds, **Mcrds, **Rcrds;
	double      **CellEdge;
	double      lambda0, lambda, dlambda;
	double      dtdn;
	int        dim = vmuscl->dim;
	int        i, k, nvar_u;

	load_sound_speed(&vmuscl->VM0st,start,end-start);

	u = vmuscl->u;
	du = vmuscl->du;
	v = u[vmuscl->index.v[0]];
	c = vmuscl->vst->c;
	dv = du[vmuscl->index.v[0]];
	dc = vmuscl->du[vmuscl->index.sound_speed];

	uM0 = vmuscl->uM0;
	vM0 = uM0[vmuscl->index.v[0]];
	cM0 = vmuscl->VM0st.c;

	uL = vmuscl->uL;
	Lcrds = vmuscl->VLst.coords;
	uM = vmuscl->uM;
	Mcrds = vmuscl->VMst.coords;
	uR = vmuscl->uR;
	Rcrds = vmuscl->VRst.coords;
	nvar_u = vmuscl->nvar_u;
	dtdn = vmuscl->dt/vmuscl->dn;
	CellEdge = vmuscl->CellEdge;
	for (k = start; k < end; ++k)
	{
	    /* uL */
	    lambda0 = vM0[k]+cM0[k];
	    if (lambda0 > 0.0)
	    {
	        lambda = v[k-1]+c[k-1];	dlambda = dv[k-1]+dc[k-1];
	        d = (lambda + 0.5*dlambda)*dtdn/(1.0 + dlambda*dtdn);
	        if (d < -0.5) d = -0.5;
	        if (d >  0.5) d =  0.5;
		for (i = 0; i < dim; ++i)
		    Lcrds[k][i] = (1.0 - d)*CellEdge[k][i] + d*CellEdge[k-1][i];
		for (i = 0; i < nvar_u; ++i)
		    uL[i][k] = u[i][k-1] + du[i][k-1]*(0.5-d);
	    }
	    else if (lambda0 < 0.0)
	    {
	        lambda = v[k]+c[k];	dlambda = dv[k]+dc[k];
	        d = (0.5*dlambda - lambda)*dtdn/(1.0 + dlambda*dtdn);
	        if (d < -0.5) d = -0.5;
	        if (d >  0.5) d =  0.5;
		for (i = 0; i < dim; ++i)
		    Lcrds[k][i] = (1.0 - d)*CellEdge[k][i] + d*CellEdge[k+1][i];
		for (i = 0; i < nvar_u; ++i)
		    uL[i][k] = u[i][k] + du[i][k]*(d-0.5);
	    }
	    else
	    {
		for (i = 0; i < dim; ++i)
		    Lcrds[k][i] = CellEdge[k][i];
		for (i = 0; i < nvar_u; ++i)
		    uL[i][k] = uM0[i][k];
	    }
	    if (uL[0][k] < 0.0)
	    {
	        if (lambda0 > 0.0)
		{
		    for (i = 0; i < dim; ++i)
		        Lcrds[k][i] = 0.5*CellEdge[k][i] + 0.5*CellEdge[k-1][i];
		    for (i = 0; i < nvar_u; ++i)
		        uL[i][k] = u[i][k-1];
		}
	        else if (lambda0 < 0.0)
		{
		    for (i = 0; i < dim; ++i)
		        Lcrds[k][i] = 0.5*CellEdge[k][i] + 0.5*CellEdge[k+1][i];
		    for (i = 0; i < nvar_u; ++i)
		        uL[i][k] = u[i][k];
		}
	    }

	    /* uM */
	    lambda0 = vM0[k];
	    if (lambda0 > 0.0)
	    {
	        lambda = v[k-1];	dlambda = dv[k-1];
	        d = (lambda + 0.5*dlambda)*dtdn/(1.0 + dlambda*dtdn);
	        if (d < -0.5) d = -0.5;
	        if (d >  0.5) d =  0.5;
		for (i = 0; i < dim; ++i)
		    Mcrds[k][i] = (1.0 - d)*CellEdge[k][i] + d*CellEdge[k-1][i];
		for (i = 0; i < nvar_u; ++i)
		    uM[i][k] = u[i][k-1] + du[i][k-1]*(0.5-d);
	    }
	    else if (lambda0 < 0.0)
	    {
	        lambda = v[k];	dlambda = dv[k];
	        d = (0.5*dlambda - lambda)*dtdn/(1.0 + dlambda*dtdn);
	        if (d < -0.5) d = -0.5;
	        if (d >  0.5) d = 0.5;
		for (i = 0; i < dim; ++i)
		    Mcrds[k][i] = (1.0 - d)*CellEdge[k][i] + d*CellEdge[k+1][i];
		for (i = 0; i < nvar_u; ++i)
		    uM[i][k] = u[i][k] + du[i][k]*(d-0.5);
	    }
	    else
	    {
		for (i = 0; i < dim; ++i)
		    Mcrds[k][i] = CellEdge[k][i];
		for (i = 0; i < nvar_u; ++i)
		    uM[i][k] = uM0[i][k];
	    }
	    if (uM[0][k] < 0.0)
	    {
	        if (lambda0 > 0.0)
		{
		    for (i = 0; i < dim; ++i)
		        Mcrds[k][i] = 0.5*CellEdge[k][i] + 0.5*CellEdge[k-1][i];
		    for (i = 0; i < nvar_u; ++i)
		        uM[i][k] = u[i][k-1];
		}
	        else if (lambda0 < 0.0)
		{
		    for (i = 0; i < dim; ++i)
		        Mcrds[k][i] = 0.5*CellEdge[k][i] + 0.5*CellEdge[k+1][i];
		    for (i = 0; i < nvar_u; ++i)
		        uM[i][k] = u[i][k];
		}
	    }

	    /* uR */
	    lambda0 = vM0[k]-cM0[k];
	    if (lambda0 > 0.0)
	    {
	        lambda = v[k-1]-c[k-1];	dlambda = dv[k-1]-dc[k-1];
	        d = (lambda + 0.5*dlambda)*dtdn/(1.0 + dlambda*dtdn);
	        if (d < -0.5) d = -0.5;
	        if (d >  0.5) d = 0.5;
		for (i = 0; i < dim; ++i)
		    Rcrds[k][i] = (1.0 - d)*CellEdge[k][i] + d*CellEdge[k-1][i];
		for (i = 0; i < nvar_u; ++i)
		    uR[i][k] = u[i][k-1] + du[i][k-1]*(0.5-d);
	    }
	    else if (lambda0 < 0.0)
	    {
	        lambda = v[k]-c[k];	dlambda = dv[k]-dc[k];
	        d = (0.5*dlambda - lambda)*dtdn/(1.0 + dlambda*dtdn);
	        if (d < -0.5) d = -0.5;
	        if (d >  0.5) d = 0.5;
		for (i = 0; i < dim; ++i)
		    Rcrds[k][i] = (1.0 - d)*CellEdge[k][i] + d*CellEdge[k+1][i];
		for (i = 0; i < nvar_u; ++i)
		    uR[i][k] = u[i][k] + du[i][k]*(d-0.5);
	    }
	    else
	    {
		for (i = 0; i < dim; ++i)
		    Rcrds[k][i] = CellEdge[k][i];
		for (i = 0; i < nvar_u; ++i)
		    uR[i][k] = uM0[i][k];
	    }
	    if (uR[0][k] < 0.0)
	    {
	        if (lambda0 > 0.0)
		{
		    for (i = 0; i < dim; ++i)
		        Rcrds[k][i] = 0.5*CellEdge[k][i] + 0.5*CellEdge[k-1][i];
		    for (i = 0; i < nvar_u; ++i)
		        uR[i][k] = u[i][k-1];
		}
	        else if (lambda0 < 0.0)
		{
		    for (i = 0; i < dim; ++i)
		        Rcrds[k][i] = 0.5*CellEdge[k][i] + 0.5*CellEdge[k+1][i];
		    for (i = 0; i < nvar_u; ++i)
		        uR[i][k] = u[i][k];
		}
	    }
	}
	if ((vmuscl->alpha != 0) && (vmuscl->idir == 0))
	{
	    double **coords = vmuscl->coords;
	    int j0;

	    for (j0 = start; j0 < end; ++j0)
	        if (coords[j0][0] > 0.0)
		    break;
	    if (j0 > start)
	    {
	        for (i = 0; i < nvar_u; ++i)
		    uL[i][j0] = uR[i][j0];
	        uL[vmuscl->index.v[0]][j0] *= -1.0;
	        uM[vmuscl->index.v[0]][j0] = 0.0;
	    }
	}
}		/*end characteristic_traceback*/

LOCAL	Vec_Muscl *plm_alloc_wk_space(
	Muscl_Opts	*mopts,
	int		vsize,
	int		dim,
	Vec_Muscl       *vmuscl)
{
	AVISC		 Avisc;
	int		 nfloats = mopts->nfloats;
	int              i, nvar_u;

	vmuscl = alloc_Vec_Muscl(vmuscl);

	vmuscl->Opts = *mopts;
	vmuscl->dim = dim;

	nvar_u = 0;
	vmuscl->index.density = nvar_u++;
	vmuscl->index.energy = nvar_u++;
	for (i = 0; i < dim; ++i)
	    vmuscl->index.v[i] = nvar_u++;
	vmuscl->index.internal_energy_density = nvar_u++;
	vmuscl->index.pressure = nvar_u++;
	vmuscl->index.sound_speed = nvar_u++;
	vmuscl->index.GAM = nvar_u++;
#if DONT_COMPILE
	if ((mopts->_rmidstate == g_linear_us_up_Riemann_midstate) ||
	    (mopts->_rmidstate == g_gamma_law_fit_Riemann_midstate))
#endif /*DONT_COMPILE*/
	    vmuscl->index.FD = nvar_u++;

	vmuscl->index.grav = nvar_u++;

	vmuscl->nvar_u = nvar_u;
	MATRIX(vmuscl,CellEdge,vsize+1,dim,FLOAT);
	VECTOR(vmuscl,u,nvar_u,sizeof(double*));/* Allocate space for        */
	MATRIX(vmuscl,du,nvar_u,vsize,FLOAT);  /* pressure, sound speed,    */
					       /* Gruneisen exponent and    */
					       /* gravity                   */
	VECTOR(vmuscl,grav,vsize,FLOAT);

	VECTOR(vmuscl,ucon,nfloats,sizeof(double*));
	if (source_terms_exist() == YES)
	    VECTOR(vmuscl,source,nfloats,sizeof(double*));

	use_artificial_dissipation(&Avisc);
	if (use_lapidus_artificial_viscosity(Avisc))
	{
	    scalar(&vmuscl->avisc,sizeof(Vec_Avisc));
	    set_alloc(vmuscl->avisc,avisc);
	    MATRIX(vmuscl->avisc,g,3,vsize,FLOAT);
	    VECTOR(vmuscl->avisc,cs_ave,vsize,FLOAT);
	    VECTOR(vmuscl->avisc,c_ave,vsize,FLOAT);
	    VECTOR(vmuscl->avisc,vn_ave,vsize,FLOAT);
	    MATRIX(vmuscl->avisc,b,nfloats,vsize,FLOAT);
	    VECTOR(vmuscl->avisc,visc,vsize,FLOAT);
	    MATRIX(vmuscl->avisc,uconM,nfloats,vsize,FLOAT);
	    set_no_alloc(vmuscl->avisc,mdlambda);
	}
	else if (use_upwind_artificial_viscosity(Avisc) ||
	         use_linear_artificial_viscosity(Avisc))
	{
	    scalar(&vmuscl->avisc,sizeof(Vec_Avisc));
	    set_alloc(vmuscl->avisc,avisc);
	    MATRIX(vmuscl->avisc,g,1,vsize,FLOAT);
	    VECTOR(vmuscl->avisc,mdlambda,vsize,FLOAT);
	    set_no_alloc(vmuscl->avisc,cs_ave);
	    set_no_alloc(vmuscl->avisc,c_ave);
	    set_no_alloc(vmuscl->avisc,vn_ave);
	    set_no_alloc(vmuscl->avisc,b);
	    set_no_alloc(vmuscl->avisc,visc);
	    set_no_alloc(vmuscl->avisc,uconM);
	}

	if (vmuscl->avisc != NULL)
	{
	    vmuscl->avisc->use_lapidus=use_lapidus_artificial_viscosity(Avisc);
	    vmuscl->avisc->use_linear = use_linear_artificial_viscosity(Avisc);
	    vmuscl->avisc->use_upwind = use_upwind_artificial_viscosity(Avisc);
	    zero_scalar(&vmuscl->Vegn,sizeof(Vec_Eigen));
	    vmuscl->Vegn.negn = 3;
	    set_no_alloc(&vmuscl->Vegn,vegn);
	    MATRIX(&vmuscl->Vegn,lambda,vmuscl->Vegn.negn,vsize,FLOAT);
	}

	if (use_muscl_slope_flattening(Avisc))
	{
	    scalar(&vmuscl->msf,sizeof(Vec_MSF));
	    set_alloc(vmuscl->msf,msf);
	    VECTOR(vmuscl->msf,chi,vsize,FLOAT);
	}
	vmuscl->monotone_reconstruction = mopts->monotone_reconstruction;
	vmuscl->link_reconstructions = mopts->link_reconstructions;

	vmuscl->max_vsize = vsize;

	        /* Set half step calculation data */

	MATRIX(vmuscl,uL,nvar_u,vsize,FLOAT); /* Allocate space for       */
	MATRIX(vmuscl,uR,nvar_u,vsize,FLOAT); /* pressure and sound speed */
	MATRIX(vmuscl,uM,nvar_u,vsize,FLOAT); /* on cell edges            */

	MATRIX(vmuscl,Flux.F,nfloats,vsize,FLOAT);
	MATRIX(vmuscl,Flux.H,nfloats,vsize,FLOAT);

	MATRIX(vmuscl,uM0,nvar_u,vsize,FLOAT);
	MATRIX(vmuscl,Flux0.F,nfloats,vsize,FLOAT);
	MATRIX(vmuscl,Flux0.H,nfloats,vsize,FLOAT);

	MATRIX(vmuscl,uM1,nvar_u,vsize,FLOAT);
	MATRIX(vmuscl,Flux1.F,nfloats,vsize,FLOAT);
	MATRIX(vmuscl,Flux1.H,nfloats,vsize,FLOAT);

	        /* Set Riemann solver date */
	vmuscl->pL = set_plm_vec_gas_storage(vmuscl->uL,&vmuscl->VLst,NULL,
	                                     vmuscl);
	vmuscl->pM = set_plm_vec_gas_storage(vmuscl->uM,&vmuscl->VMst,NULL,
	                                     vmuscl);
	vmuscl->pR = set_plm_vec_gas_storage(vmuscl->uR,&vmuscl->VRst,NULL,
	                                     vmuscl);
	vmuscl->pM0=set_plm_vec_gas_storage(vmuscl->uM0,&vmuscl->VM0st,
					    vmuscl->CellEdge,vmuscl);
	vmuscl->pM1=set_plm_vec_gas_storage(vmuscl->uM1,&vmuscl->VM1st,
	                                    vmuscl->CellEdge,vmuscl);

	VECTOR(vmuscl,A,vsize+1,FLOAT);
	VECTOR(vmuscl,dV,vsize,FLOAT);

	vmuscl->gL = vmuscl->uL[vmuscl->index.grav];
	vmuscl->gM = vmuscl->uM[vmuscl->index.grav];
	vmuscl->gR = vmuscl->uR[vmuscl->index.grav];
	vmuscl->gM0 = vmuscl->uM0[vmuscl->index.grav];
	vmuscl->gM1 = vmuscl->uM1[vmuscl->index.grav];

	return vmuscl;
}		/*end plm_alloc_wk_space*/


LOCAL	double *set_plm_vec_gas_storage(
	double	  **u,
	Vec_Gas	  *vst,
	double     **coords,
	Vec_Muscl *vmuscl)
{
	int dim = vmuscl->dim;
	int	  vsize = vmuscl->max_vsize;
	int i;

	zero_scalar(vst,sizeof(Vec_Gas));
	ASSIGN_ARRAY_POINTER(vst,rho,u[vmuscl->index.density]);
	ASSIGN_ARRAY_POINTER(vst,e,u[vmuscl->index.energy]);
	ASSIGN_ARRAY_POINTER(vst,re,u[vmuscl->index.internal_energy_density]);
	VECTOR(vst,v,dim,sizeof(double*));
	for (i = 0; i < dim; ++i)
	    vst->v[i] = u[vmuscl->index.v[i]];
	ASSIGN_ARRAY_POINTER(vst,p,u[vmuscl->index.pressure]);
	ASSIGN_ARRAY_POINTER(vst,c,u[vmuscl->index.sound_speed]);
	ASSIGN_ARRAY_POINTER(vst,GAM,u[vmuscl->index.GAM]);
	if (vmuscl->index.FD >= 0)
	    ASSIGN_ARRAY_POINTER(vst,FD,u[vmuscl->index.FD]);
	else
	    vst->FD = NULL;
	VECTOR(vst,prms_jmp,vsize+1,INT);
	VECTOR(vst,c2,vsize,FLOAT);
	if (coords != NULL)
	    ASSIGN_ARRAY_POINTER(vst,coords,coords);
	else
	    MATRIX(vst,coords,vsize,dim,FLOAT);
	return vst->p;
}		/*end set_plm_vec_gas_storage*/

LOCAL void g_plm_compute_eigens(
        int             start,
	int             end,
	Vec_Muscl       *vmuscl)
{
	if (vmuscl->Vegn.lambda)
	{
	    Vec_Gas    *vst = vmuscl->vst;
	    int        offset = vmuscl->offset;
	    double      *v0 = vst->v[0] + offset;
	    double      *c = vst->c + offset;
	    int        j;

	    for (j = start; j < end; ++j)
	    {
	        vmuscl->Vegn.lambda[0][j] = v0[j] - c[j];
	        vmuscl->Vegn.lambda[1][j] = v0[j];
	        vmuscl->Vegn.lambda[2][j] = v0[j] + c[j];
	    }
#if defined(DEBUG_PLM)
	    if (debugging("plm"))
	    {
		(void) printf("%-14s %-14s %-14s %-14s\n",
			      "v - c","v","v + c","c");
	        for (j = start; j < end; ++j)
		    (void) printf("%-14g %-14g %-14g %-14g\n",
	                          vmuscl->Vegn.lambda[0][j],
	                          vmuscl->Vegn.lambda[1][j],
	                          vmuscl->Vegn.lambda[2][j],
				  c[j]);
		(void) printf("\n");
	    }
#endif /* defined(DEBUG_PLM) */
	}
}		/*end g_plm_compute_eigens*/

/*
*			g_plm_add_art_visc1();
*
*	This function adds the first half of artificial viscosity.  It is
*	added to the mid states which will be used to compute conservative
*	fluxes and then in the conservative difference scheme.
*/

LOCAL void g_plm_add_art_visc1(
	int		start,
	int		end,
	double           **uM,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc *avisc;
	int	  dim;
	int	  j, k;
	double	  *g2, *g1;
	double	  *p, **g;
	double     *rho,  **v, **m, *E, *v0, *m0;
	double     *rhoM, **vM, *vM0, *eM;
	double     ke;

	if ((avisc = vmuscl->avisc) == NULL)
	    return;
	g = avisc->g;
	g2 = g[2], g1 = g[1];
	if (g2 == NULL || g1 == NULL)
	    return;

	dim =		vmuscl->dim;
	rho = vmuscl->ucon[vmuscl->index.density];
	E   = vmuscl->ucon[vmuscl->index.energy];
	m   = vmuscl->ucon+vmuscl->index.v[0];
	v   = vmuscl->u + vmuscl->index.v[0];
	p   = vmuscl->u[vmuscl->index.pressure];
	m0 = m[0];
	v0 = v[0];

	rhoM = uM[vmuscl->index.density];
	vM   = uM + vmuscl->index.v[0];
	vM0  = vM[0];
	eM   = uM[vmuscl->index.energy];

	for (j = start; j < end; ++j)
	{
	    for (ke = 0.0, k = 0; k < dim; ++k)
	    {
	        ke += 0.5*vM[k][j]*vM[k][j];
	        vM[k][j] *= rhoM[j];
	    }
	    eM[j] = rhoM[j]*(ke + eM[j]);
	}
	for (j = start; j < end; ++j)
	{
	    rhoM[j] -= 0.5*(g2[j-1]*(m0[j] - m0[j-1]) +
	                    g1[j-1]*(rho[j]-rho[j-1]));
	    vM0[j] -= 0.5*(g2[j-1]*(m0[j]*v0[j]+p[j]-m0[j-1]*v0[j-1]-p[j-1]) +
	                   g1[j-1]*(m0[j]-m0[j-1]));
	    for (k = 1; k < dim; ++k)
	        vM[k][j] -= 0.5*(g2[j-1]*(m[k][j]*v0[j]-m[k][j-1]*v0[j-1]) +
	                         g1[j-1]*(m[k][j]-m[k][j-1]))/rhoM[j];
	    eM[j] -= 0.5*(g2[j-1]*((E[j]+p[j])*v0[j]-(E[j-1]+p[j-1])*v0[j-1]) +
	                    g1[j-1]*(E[j]-E[j-1]));
	}
	for (j = start; j < end; ++j)
	{
	    for (ke = 0.0, k = 0; k < dim; ++k)
	    {
	        vM[k][j] /= rhoM[j];
	        ke += 0.5*vM[k][j]*vM[k][j];
	    }
	    eM[j] = eM[j]/rhoM[j] - ke;
	    if (eM[j] < Min_energy(vmuscl->vst->state[j]))
	        eM[j] = Min_energy(vmuscl->vst->state[j]);
	}
}		/*end g_plm_add_art_visc1*/

LOCAL void g_plm_add_art_visc2(
        int             start,
	int             end,
	Vec_Muscl       *vmuscl)
{
	Vec_Avisc *avisc;
	double     **ucon, **H;
	double     *g0;
	int       k, j, dim;

	if ((avisc = vmuscl->avisc) == NULL)
	    return;
	dim = vmuscl->dim;
	g0 = avisc->g[0];
	H = vmuscl->Flux.H;
	ucon = vmuscl->ucon;
	for (j = start; j < end; ++j)
	    for (k = 0; k < dim+2; ++k)
	        H[k][j] -= 0.5*g0[j-1]*(ucon[k][j] - ucon[k][j-1]);
}		/*end g_plm_add_art_visc2*/

#if defined(DEBUG_PLM)
LOCAL	void	print_data_header(
	int dim,
	const char *rho,
	const char *e,
	const char *v,
	const char *p,
	const char *c,
	const char *mesg)
{
	int i;

	(void) output();
	(void) printf("%-5s ","index");
	for (i = 0; i < dim; ++i)
	    (void) printf("x%-13d ",i);
	(void) printf("%-14s ",rho);
	(void) printf("%-14s ",e);
	for (i = 0; i < dim; ++i)
	{
	    char vname[14];
	    (void) sprintf(vname,"%s%d",v,i);
	    (void) printf("%-14s ",vname);
	}
	(void) printf("%-14s ",p);
	(void) printf("%-14s ",c);
	(void) printf("%s\n",mesg);
}		/*end print_data_header*/

LOCAL	void	print_data_line(
	int   dim,
	int   j,
	double **crds,
	double *rho,
	double *e,
	double **v,
	double *p,
	double *c)
{
	int i;
	(void) printf(" %-5d ",j);
	for (i = 0; i < dim; ++i)
	    (void) printf("%-14g ",crds[j][i]);
	(void) printf("%-14g ",rho[j]);
	(void) printf("%-14g ",e[j]);
	for (i = 0; i < dim; ++i)
	    (void) printf("%-14g ",v[i][j]);
	(void) printf("%-14g ",p[j]);
	(void) printf("%-14g ",c[j]);
	(void) printf("\n");
}		/*end print_data_line*/

LOCAL	boolean plm_interior_sten_reg(
	int		is,
	int		idir,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	Front		*newfr = sten->newfr;
	int 		i, imax;
	int 		**icoords = sten->icoords, icrds[MAXD];
	int             endpt = stencil_radius(wave);
	int             vsten_rad = vsten_radius(wave);
	int             dim = newfr->interf->dim;
	int             gmax = newfr->rect_grid->gmax[idir];
	int             lbuf = newfr->rect_grid->lbuf[idir];
	int             ubuf = newfr->rect_grid->ubuf[idir];
	GRID_DIRECTION	prev_side, next_side;
	COMPONENT 	cmp;
	COMPONENT	new_comp = sten->newcomp, *comp = sten->comp;
	CRXING		*cross;

	for (i = -endpt; i <= endpt; ++i)
	    if ((equivalent_comps(comp[i],new_comp,newfr->interf) == NO) ||
		(sten->nc[i] != 0))
	    return NO;
	if (endpt >= vsten_rad)
	    return YES;

	for (i = 0; i < dim; ++i)
	    icrds[i] = icoords[0][i];

	imax = gmax + ubuf;
	switch (idir)
	{
	case 0:
	    prev_side = WEST;
	    next_side = EAST;
	    break;
	case 1:
	    prev_side = SOUTH;
	    next_side = NORTH;
	    break;
	case 2:
	    prev_side = LOWER;
	    next_side = UPPER;
	    break;
	}
	for (i = endpt+1; i <= vsten_rad; ++i)
	{
	    icrds[idir] = is - i;
	    cmp = (((is-i) < -lbuf) || ((is-i) >= imax)) ? 
	        exterior_component(newfr->interf) : Rect_comp(icrds,wave);
	    cross = Rect_crossing(icrds,next_side,wave);
	    if ((equivalent_comps(cmp,new_comp,newfr->interf) == NO) ||
		(cross != NULL))
	    	return NO;
		
	    icrds[idir] = is + i;
	    cmp = (((is+i) < -lbuf) || ((is+i) >= imax)) ? 
	        exterior_component(newfr->interf) : Rect_comp(icrds,wave);
	    cross = Rect_crossing(icrds,prev_side,wave);
	    if ((equivalent_comps(cmp,new_comp,newfr->interf) == NO) ||
		(cross != NULL))
		return NO;
	}
	return YES;
}		/*end plm_interior_sten_reg*/
#endif /* defined(DEBUG_PLM) */
