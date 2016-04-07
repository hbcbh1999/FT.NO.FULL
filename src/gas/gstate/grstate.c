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
*				grstate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Routines for the computation and evaluation of random state
*	fields.
*
*/


#include <gdecs/gdecs.h>

/* LOCAL Function prototypes */
LOCAL	void	independent_random_states(double,RANDOM_STATE*,INTERFACE*);
LOCAL	void	create_next_random_level(RANDOM_STATE*,INTERFACE*);

/*
*			random_velocity_inlet():
*	
*	Assigns the boundary state at the exterior point coords using
*	a random velocity perturbation of an ambient inlet state.
*
*/

/*ARGSUSED*/
EXPORT	void random_velocity_inlet(
	double		*coords,
	HYPER_SURF	*hs,
	Front		*front,
	POINTER		p2wave,
	Locstate	state)
{
	RANDOM_STATE    *rstate = (RANDOM_STATE*) boundary_state_data(hs);
	INTERFACE	*intfc = front->interf;
	RECT_GRID 	*rgr = &rstate->grid;
	RECT_GRID	*cgr = front->rect_grid;
	Locstate	cell_states[8];
	boolean		reflect_velocity[3];
	double		t[2], f[3];
	double		alpha[8];
	double		crds[3];
	int		*gmax = rgr->gmax;
	int		i, dim = cgr->dim, ic[3];
	int		ip, jp, kp;
	int		real_dim = dim;
	static Locstate state_at_level[2] = {NULL,NULL};

	if (state_at_level[0] == NULL)
	{
	    alloc_state(front->interf,&state_at_level[0],front->sizest);
	    alloc_state(front->interf,&state_at_level[1],front->sizest);
	}

	if (front->time > rstate->tlast)
	{
	    Locstate ***tmp_st;
	    tmp_st = rstate->old_st;
	    rstate->old_st = rstate->save_st;
	    rstate->save_st = tmp_st;
	    rstate->time_of_save = rstate->tlast - rstate->delta_t;
	    for (i = 0; i < 3; i++)
	    	rstate->seed_after_save[i] = rstate->seed_after_old[i];
	    while (front->time > rstate->tlast)
	    	create_next_random_level(rstate,intfc);
	}
	t[0] = (rstate->tlast - front->time)/rstate->delta_t;
	t[1] = 1.0 - (rstate->tlast - front->time)/rstate->delta_t;

	for (i = 0; i < 3; i++)
	{
	    ic[i] = 0;
	    crds[i] = 0.0;
	    reflect_velocity[i] = NO;
	}
	for (i = 0; i < dim; i++)
	{
	    double d;
	    crds[i] = coords[i];
	    d = crds[i] - cgr->GL[i];
	    if (d < 0.0)
	    {
	        if (rect_boundary_type(intfc,i,0) == REFLECTION_BOUNDARY)
	        {
	            crds[i] = cgr->GL[i] - d; 
	            reflect_velocity[i] = YES;
	        }
	        else
	            crds[i] = cgr->GL[i];
	    }
	    d = cgr->GU[i] - crds[i];
	    if (d < 0.0)
	    {
	        if (rect_boundary_type(intfc,i,1) == REFLECTION_BOUNDARY)
	        {
	            crds[i] = cgr->GU[i] + d; 
	            reflect_velocity[i] = YES;
	        }
	        else
	            crds[i] = cgr->GU[i];
	    }
	}

	if (rect_in_which(crds,ic,rgr) == FUNCTION_FAILED)
	{
	    screen("ERROR in random_velocity_inlet(), "
	           "coords outside of grid\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < dim; i++)
	{
	    if (cell_width(ic[i],i,rgr) == 0.0)
	    {
	    	f[i] = 0.0;
	    	real_dim--;
	    }
	    else
	    	f[i] = (crds[i]-cell_edge(ic[i],i,rgr))/cell_width(ic[i],i,rgr);
	}
	switch (real_dim)
	{
	case 0:
	    ft_assign(state_at_level[0],rstate->old_st[ic[0]][0][0],front->sizest);
	    ft_assign(state_at_level[1],rstate->new_st[ic[0]][0][0],front->sizest);
	    break;
	case 1:
	    alpha[0] = 1.0 - f[0];
	    alpha[1] = f[0];
	    ip = (ic[0] < gmax[0]) ? 1 : 0;
	    cell_states[0] = rstate->old_st[ic[0]][0][0];
	    cell_states[1] = rstate->old_st[ic[0]+ip][0][0];
	    g_linear_combination_of_states(alpha,cell_states,
	    		                   2,state_at_level[0]);
	    cell_states[0] = rstate->new_st[ic[0]][0][0];
	    cell_states[1] = rstate->new_st[ic[0]+ip][0][0];
	    g_linear_combination_of_states(alpha,cell_states,
	    		                   2,state_at_level[1]);
	    break;
	case 2:
	    alpha[0] = (1.0 - f[0])*(1.0 - f[1]);
	    alpha[1] = f[0]*(1.0 - f[1]);
	    alpha[2] = (1.0 - f[0])*f[1];
	    alpha[3] = f[0]*f[1];
	    for (i = 0; i < 4; i++)
	    {
	    	ip = (ic[0] < gmax[0]) ? i%2 : 0;
	    	jp = (ic[1] < gmax[1]) ? i/2 : 0;
	    	cell_states[i] = rstate->old_st[ic[0]+ip][ic[1]+jp][0];
	    }
	    g_linear_combination_of_states(alpha,cell_states,
	        	                   4,state_at_level[0]);
	    for (i = 0; i < 4; i++)
	    {
	    	ip = (ic[0] < gmax[0]) ? i%2 : 0;
	    	jp = (ic[1] < gmax[1]) ? i/2 : 0;
	    	cell_states[i] = rstate->new_st[ic[0]+ip][ic[1]+jp][0];
	    }
	    g_linear_combination_of_states(alpha,cell_states,
	        	                   4,state_at_level[1]);
	    break;
	case 3:
	    alpha[0] = (1.0 - f[0])*(1.0 - f[1])*(1.0 - f[2]);
	    alpha[1] = f[0]*(1.0 - f[1])*(1.0 - f[2]);
	    alpha[2] = (1.0 - f[0])*f[1]*(1.0 - f[2]);
	    alpha[3] = f[0]*f[1]*(1.0 - f[2]);
	    alpha[4] = (1.0 - f[0])*(1.0 - f[1])*f[2];
	    alpha[5] = f[0]*(1.0 - f[1])*f[2];
	    alpha[6] = (1.0 - f[0])*f[1]*f[2];
	    alpha[7] = f[0]*f[1]*f[2];
	    for (i = 0; i < 8; i++)
	    {
	    	ip = (ic[0] < gmax[0]) ? i%2 : 0;
	    	jp = (ic[1] < gmax[1]) ? (i%4)/2 : 0;
	    	kp = (ic[2] < gmax[2]) ? i%4 : 0;
	    	cell_states[i] = rstate->old_st[ic[0]+ip][ic[1]+jp][ic[2]+kp];
	    }
	    g_linear_combination_of_states(alpha,cell_states,
	        	                   8,state_at_level[0]);
	    for (i = 0; i < 8; i++)
	    {
	    	ip = (ic[0] < gmax[0]) ? i%2 : 0;
	    	jp = (ic[1] < gmax[1]) ? (i%4)/2 : 0;
	    	kp = (ic[2] < gmax[2]) ? i%4 : 0;
	    	cell_states[i] = rstate->new_st[ic[0]+ip][ic[1]+jp][ic[2]+kp];
	    }
	    g_linear_combination_of_states(alpha,cell_states,
	        	                   8,state_at_level[1]);
	    break;
	}
	g_linear_combination_of_states(t,state_at_level,2,state);
	if ((axisymmetric_random_region_about_origin(rstate,intfc) == YES) &&
		(RadialVelocityDecayScale(rstate) > 0.0))
	{
	    double decay_factor;
	    double r, a;
	    r = fabs(crds[0] - cgr->GL[0])/RadialVelocityDecayScale(rstate);
	    a = RadialVelocityDecayExponent(rstate);
	    decay_factor = 1.0 - exp(-pow(r,a));
	    Vel(state)[0] *= decay_factor;
	}
	for (i = 0; i < dim; i++)
	{
	    if (reflect_velocity[i] == YES)
	        Vel(state)[i] = -Vel(state)[i];
	}
	set_state(state,GAS_STATE,state);
	if (is_bad_state(state,YES,"random_velocity_inlet"))
	{
	    screen("ERROR in random_velocity_inlet(), bad state produced\n");
	    fprint_raw_gas_data(stdout,state,current_interface()->dim);
	    clean_up(ERROR);
	}
}		/*end random_velocity_inlet*/

EXPORT	void	generate_random_region(
	double		wgt,
	Locstate	***corr_states,
	RANDOM_STATE	*rstate,
	INTERFACE	*intfc)
{
	Locstate corr_state;
	int	*lbuf = rstate->grid.lbuf;
	int	*ubuf = rstate->grid.ubuf;
	size_t	sizest = Params(Mean(rstate))->sizest;
	int	gmax[3];
	int	l, n, m, ll, nn, mm, k;
	int	i;
	static	Locstate *state = NULL;
	static	double *alpha = NULL;
	static	int num_states = 0;

	if (state == NULL)
	{
	    num_states = rstate->M_ell;
	    uni_array(&state,num_states,sizeof(Locstate));
	    uni_array(&alpha,num_states,FLOAT);
	}
	if (num_states < rstate->M_ell)
	{
	    free_these(2,state,alpha);
	    num_states = rstate->M_ell;
	    uni_array(&state,num_states,sizeof(Locstate));
	    uni_array(&alpha,num_states,FLOAT);
	}
	for (i = 0; i < rstate->M_ell; i++)
	    alpha[i] = 1.0/rstate->M_ell;

	independent_random_states(wgt,rstate,intfc);
	for (i = 0; i < 3; i++)
	    gmax[i] = rstate->grid.gmax[i]+1;
	for (l = 0; l < gmax[0]; l++)
	for (m = 0; m < gmax[1]; m++)
	for (n = 0; n < gmax[2]; n++)
	{
	    corr_state = corr_states[l][m][n];
	    (*Params(Mean(rstate))->_clear_state)(state,sizest);
	    Set_params(state,Mean(rstate));
	    for (k = 0, ll = -lbuf[0]; ll <= ubuf[0]; ll++)
	    for (nn = -lbuf[1]; nn <= ubuf[1]; nn++)
	    for (mm = -lbuf[2]; mm <= ubuf[2]; mm++)
	    {
	    	if (! in_correlation_ellipse(ll,mm,nn,rstate))
	    		continue;
	    	state[k++] =
	            rstate->indep_st[l+ll+lbuf[0]][m+mm+lbuf[1]][n+nn+lbuf[2]];
	    }
	    g_linear_combination_of_states(alpha,state,rstate->M_ell,
	        	                   corr_state);
	}
}		/*end generate_random_region*/

EXPORT	void	relax_random_inlet_level(
	Locstate	***old_st,
	Locstate	***new_st,
	RANDOM_STATE    *rstate)
{
	Locstate state[2];
	double alpha[2];
	int gmax[3];
	int i, l, n, m;
	size_t sizest = Params(Mean(rstate))->sizest;
	static	Locstate new_state = NULL;

	if (new_state == NULL)
	    (*Params(Mean(rstate))->_alloc_state)(&new_state,sizest);

	alpha[1] = 1.0/rstate->N;
	alpha[0] = 1.0 - alpha[1];

	for (i = 0; i < 3; i++)
	    gmax[i] = rstate->grid.gmax[i]+1;
	for (l = 0; l < gmax[0]; l++)
	for (m = 0; m < gmax[1]; m++)
	for (n = 0; n < gmax[2]; n++)
	{
	    state[0] = old_st[l][m][n];
	    state[1] = new_st[l][m][n];
	    g_linear_combination_of_states(alpha,state,2,new_state);
	    ft_assign(state[1],new_state,sizest);
	    if (debugging("bad_state"))
	    {
		if (is_bad_state(state[1],YES,"relax_random_inlet_level"))
		{
		    screen("ERROR in relax_random_inlet_level(), "
			   "bad state detected\n");
	            fprint_raw_gas_data(stdout,state[1],current_interface()->dim);
		    clean_up(ERROR);
		}
	    }
	}
}		/*end relax_random_inlet_level*/

EXPORT	int	in_correlation_ellipse(
	int	l,
	int	n,
	int	m,
	RANDOM_STATE	*rstate)
{
	RECT_GRID *grid = &rstate->grid;
	double	*h = grid->h;
	double	x[3];
	double	**A = rstate->A;
	double	norm = 0.0;
	int i, j;

	x[0] = l*h[0]; x[1] = n*h[1]; x[2] = m*h[2];
	for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    norm += x[i]*A[i][j]*x[j];
	
	return (norm <= 1.0);
}		/*end in_correlation_ellipse*/

EXPORT	boolean axisymmetric_random_region_about_origin(
	RANDOM_STATE *rstate,
	INTERFACE    *intfc)
{
	RECT_GRID *cgr = computational_grid(intfc);
	RECT_GRID *rgr = &rstate->grid;
	int       dim = cgr->dim;

	if (dim != 2)
	    return NO;
	if (cgr->Remap.remap != CYLINDRICAL_REMAP)
	    return NO;
	if (cgr->L[0] < rgr->L[0])
	    return NO;
	return YES;
}		/*end axisymmetric_random_region_about_origin*/

LOCAL	void	create_next_random_level(
	RANDOM_STATE *rstate,
	INTERFACE    *intfc)
{
	Locstate ***tmp_st;
	int i;

	tmp_st = rstate->new_st;
	rstate->new_st = rstate->old_st;
	rstate->old_st = tmp_st;
	for (i = 0; i < 3; i++)
	    rstate->seed_after_old[i] = rstate->xsubi[i];
	generate_random_region(sqrt(2.0*rstate->N-1.0),
	    	               rstate->new_st,rstate,intfc);
	relax_random_inlet_level(rstate->old_st,rstate->new_st,rstate);
	rstate->tlast += rstate->delta_t;
}		/*end create_next_random_level*/

/*ARGSUSED*/
LOCAL	void independent_random_states(
	double        wgt,
	RANDOM_STATE *rstate,
	INTERFACE    *intfc)
{
	RECT_GRID	*grid = &rstate->grid;
	Locstate	***indep_st = rstate->indep_st;
	Locstate	state;
	Locstate	mu, sigma;
	boolean		do_radial_decay;
	double		vm, vs;
	int		gmax[3];
	int		l, n, m, i, dim = Params(Mean(rstate))->dim;
	unsigned short int *xsubi = rstate->xsubi;

	do_radial_decay = axisymmetric_random_region_about_origin(rstate,intfc);
	wgt *= sqrt((double) rstate->M_ell);
	mu = Mean(rstate);
	sigma = Sigma(rstate);
	for (i = 0; i < 3; i++)
	    gmax[i] = grid->gmax[i]+1+grid->lbuf[i]+grid->ubuf[i];
	for (l = 0; l < gmax[0]; l++)
	for (m = 0; m < gmax[1]; m++)
	for (n = 0; n < gmax[2]; n++)
	{
	    state = indep_st[l][m][n];
	    Set_params(state,mu);
	    set_type_of_state(state,state_type(mu));
	    do
	    {
	        Dens(state) = random_gaussian(Dens(mu),wgt*Dens(sigma),xsubi);
	        Energy(state) = random_gaussian(Energy(mu),
	        			        wgt*Energy(sigma),xsubi);
		vm = Mom(mu)[0];
		vs = Mom(sigma)[0];
		if ((do_radial_decay == YES) &&
		    (RadialVelocityDecayScale(rstate) > 0.0))
		{
		    double decay_factor;
		    double r, a;
		    r = (grid->GL[0] + (l+0.5)*grid->h[0]) /
			RadialVelocityDecayScale(rstate);
		    a = RadialVelocityDecayExponent(rstate);
		    decay_factor = 1.0 - exp(-pow(r,a));
		    vm *= decay_factor;
		    vs *= decay_factor;
		}
		Mom(state)[0] = random_gaussian(vm,wgt*vs,xsubi);
	        for (i = 1; i < dim; i++)
	        {
	            Mom(state)[i] = random_gaussian(Mom(mu)[i],
	        			            wgt*Mom(sigma)[i],xsubi);
	        }
		reset_gamma(state);
	    } while ((is_bad_state(state,NO,"independent_random_state")) ||
		     (scalar_product(Mom(state),Mom(mu),dim) < 0.0));
	}
}		/*end independent_random_state*/
