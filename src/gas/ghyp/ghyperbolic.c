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
*				ghyperbolic.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains g_tri_integral() and g_quad_integral, set to funciton
*		pointers in the wave structure.
*
*	Contains routines used in the hyperbolic step.
*
*/

#include <ghyp/ghyp.h>

#define VAPOR_COMP    4
#define AMBIENT_COMP  3

	/* LOCAL Function Declarations */
LOCAL	boolean	too_many_bad_gas_states(const char*,int,Front*);
LOCAL	void	g_init_obstacle_states(Vec_Gas*,int,int);
LOCAL	int	detect_ambient_vapor_mix(int,Stencil*);
LOCAL   int     local_LF_npt_tang_solver_switch(double,Tan_stencil*,Front*);
LOCAL   int     on_what = 0;


/*ARGSUSED*/
EXPORT	void point_FD(
	double		dh,
	double		dt,
	Locstate	ans,
	const double	*dir,
	int		swp_num,
	int		*iperm,
	int		*index,
	Stencil		*sten)
{
	Front		*fr = sten->fr;
	Front		*newfr = sten->newfr;
	RECT_GRID	*gr = fr->rect_grid;
	Locstate	st, *state;
	Wave 		*wave = sten->wave;
	Wave 		*newwave = sten->newwave;
	double		*rho, *en_den;
	double		*m[MAXD];
	double		**coords;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double		*vacuum_dens;
	double		*min_pressure;
	double		*min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	static double	alpha;
	int		i, j;
	int		dim = gr->dim;
	int		idirs[MAXD];
	static Vec_Gas	*vst = NULL;
	static Vec_Src *src = NULL;
	static int	endpt, npts;
	static double    **Q = NULL;

	if (is_obstacle_state(sten->st[0])) 
	{
	    g_obstacle_state(ans,sten->fr->sizest);
	    return;
	}
	if (vst == NULL) 
	{
	    /* Assumes stencil size never changes */

	    npts = sten->npts;
	    endpt = npts/2;
	    alloc_phys_vecs(wave,npts);
	    vst = g_wave_vgas(wave);
	    src = g_wave_vsrc(wave);
	    g_wave_vgas(wave) = NULL;
	    g_wave_vsrc(wave) = NULL;
	    bi_array(&Q,3,3,FLOAT);
	    vst->Q = (const double* const*)Q;
	    if (is_rotational_symmetry())
		alpha = rotational_symmetry();
	}

	if (RegionIsFlowSpecified(ans,sten->st[0],Coords(sten->p[0]),
				  sten->newcomp,sten->newcomp,fr))
	    return;

        /* Hard wired code */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(YES == detect_ambient_vapor_mix(iperm[swp_num],sten))
            {
                Gas_param   **prms = NULL, *vaporprms = NULL;
                uint64_t    igas;

                return_params_list(&prms);
                for (igas = 1; prms && *prms; ++igas, ++prms)
                {
                    if (igas == VAPOR_COMP-1)
                    {
                        vaporprms = *prms;
                        break;
                    }
                }
                return_params_list(&prms);
                for (igas = 1; prms && *prms; ++igas, ++prms)
                {
                    if (igas == AMBIENT_COMP-1)
                        break;
                }
                for(i = -endpt; i <= endpt; ++i)
                {
                    {
                        if(Params(sten->st[i]) == vaporprms)
                        {
                            pdens(sten->st[i])[0] = 0.0;
                            pdens(sten->st[i])[1] = Dens(sten->st[i]);
                            Params(sten->st[i]) = *prms;
                        }
                    }
                }
            }
        }

	clear_Vec_Gas_set_flags(vst);
	for (i = 0; i < 3; ++i)
	    for (j = 0; j < 3; ++j)
	        Q[i][j] = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    idirs[i] = iperm[(i+swp_num)%dim];
	    m[i] = vst->m[i];
	    Q[i][idirs[i]] = 1.0;
	}
	for (; i < 3; ++i)
	    Q[i][i] = 1.0;
	rho = vst->rho;
	state = vst->state;
	coords = vst->coords;
	en_den = vst->en_den;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	vacuum_dens = vst->vacuum_dens;
	min_pressure = vst->min_pressure;
	min_energy = vst->min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	for (i = 0; i < npts; ++i)
	{
	    state[i] = st = sten->ststore[i];
	    coords[i] = Coords(sten->pstore[i]);
	    rho[i] = Dens(st);
	    en_den[i] = Energy(st);
	    for (j = 0; j < dim; ++j)
	    	m[j][i] = Mom(st)[idirs[j]];
	    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                /* Include partial density into vst */
                Gas_param *params = Params(st);
                int    num_comps;
                double  *prho;
                double  **rho0 = vst->rho0;
                if((num_comps = params->n_comps) != 1)
                {
                    prho = pdens(st);
                    for(j = 0; j < num_comps; j++)
                        rho0[j][i] = prho[j];
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
            if(Params(sten->ststore[0])->n_comps != 1)
                Vec_Gas_field_set(vst,rho0) = YES;
        }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Vec_Gas_field_set(vst,vacuum_dens) = YES;
	Vec_Gas_field_set(vst,min_pressure) = YES;
	Vec_Gas_field_set(vst,min_energy) = YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	set_params_jumps(vst,0,npts);

	/* If using cylindrical coordinates and this is the
	 *  r-sweep include radius information for source
	 *  computation.  
	 */

	if (is_rotational_symmetry() && alpha > 0.0 && iperm[swp_num]==0)   
	{
	    double *radii = src->radii;

	    src->rmin = fabs(pos_radius(0.0,gr));
	    for (i = 0; i < npts; ++i)
	    	radii[i] = pos_radius(Coords(sten->p[0])[0]+(i-endpt)*dh,gr);
	}

	oned_interior_scheme(swp_num,iperm,sten->icoords[0],wave,newwave,
		             fr,newfr,sten,0,npts,vst,src,dt,dh,dim);

	Dens(ans) = vst->rho[endpt];
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vst->state[endpt]);
            if(params->n_comps != 1)
            {
                for(j = 0; j < params->n_comps; j++)
                    pdens(ans)[j] = vst->rho0[j][endpt];
            }
        }
	Energy(ans) = vst->en_den[endpt];
	for (i = 0; i < dim; ++i)
	    Mom(ans)[idirs[i]] = vst->m[i][endpt];
	Set_params(ans,vst->state[endpt]);
	set_type_of_state(ans,GAS_STATE);
	if (is_bad_state(ans,NO,"point_FD"))
	{
	    LF(dh,dt,ans,dir,swp_num,iperm,NULL,sten);
	}
	reset_gamma(ans);
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") && is_bad_state(ans,YES,"point_FD"))
	{
	    screen("ERROR in point_FD(), bad state found\n");
	    verbose_print_state("ans",ans);
	    print_Stencil(sten);
	    for (i = 0; i < npts; ++i)
	    {
		(void) printf("state[%d]\n",i);
		(void) printf("icoords %d %d\n",sten->icoords[i-endpt][0],sten->icoords[i-endpt][1]);
	        fprint_raw_gas_data(stdout,state[i],dim);
	    }
	    for (i = 0; i < npts; ++i)
	    {
		(void) printf("state[%d]",i);
	        verbose_print_state("",state[i]);
	    }
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */

}		/*end point_FD*/


/*
*			g_two_side_npt_tang():
*	
*	Calls one_side_npt_tang_solver() first on the left side and then on the
*	right side of the curve curve.
*/


EXPORT void g_two_side_npt_tang_solver(
	double		ds,
	double		dt,
	Tan_stencil	*sten,
	Locstate	ansl,
	Locstate	ansr,
	Front		*fr)
{
	INTERFACE	*intfc;
	COMPONENT	lcomp, rcomp;
	HYPER_SURF	*hs = sten->newhs;
	double		vn, v[MAXD], n[MAXD];
	int		i;
	static	size_t	sizest = 0;
	static	int	dim = 0;
        /* Test code */
        int             local_switch = NO;
        void (*one_side_npt_tang_solver)(double,double,Tan_stencil*,Locstate,
                        struct _Front*);
        static int      on_switch = 0, off_switch = 0;

	if (hs == NULL)
	    return;
	intfc = hs->interface;
	lcomp = negative_component(hs);
	rcomp = positive_component(hs);
	if (dim == 0)
	{
	    sizest = fr->sizest;
	    dim = intfc->dim;
	}

	if ((dim==2) && (sten->hs[0] != NULL))
	{
	    Locstate	sl, sr;

	    slsr(sten->p[0],sten->hse[0],sten->hs[0],&sl,&sr);
	    copy_state(sten->leftst[0],sl);
	    copy_state(sten->rightst[0],sr);
	}

        /* Test code */
        if(debugging("local_LF"))
        {
            local_switch = local_LF_npt_tang_solver_switch(ds,sten,fr);
            if(local_switch == YES)
            {
                one_side_npt_tang_solver =
                     fr->_one_side_npt_tang_solver;
                fr->_one_side_npt_tang_solver = LFoblique;
            }
        }
	switch (wave_type(hs))
	{
	case SUBDOMAIN_BOUNDARY:
	case PASSIVE_BOUNDARY:
	    g_obstacle_state(ansl,sizest);
	    g_obstacle_state(ansr,sizest);
	    break;

	/* Should this be applied in general? */
	case	NEUMANN_BOUNDARY:
	    dim = fr->rect_grid->dim;
	    normal(sten->p[0],sten->hse[0],sten->hs[0],n,fr);
	    if (is_excluded_comp(lcomp,intfc) == YES)
	    {
		sten->comp = rcomp;
		sten->states = sten->rightst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansr,fr);
		if (no_slip(hs))
		{
		    double alpha = 1.0 - adherence_coeff(hs);
		    alpha_state_velocity(alpha,ansr,dim);
		}
	    	zero_normal_velocity(ansr,n,dim);
	    	g_obstacle_state(ansl,sizest);
	    }
	    else
	    {
		sten->comp = lcomp;
		sten->states = sten->leftst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansl,fr);
		if (no_slip(hs))
		{
		    double alpha = 1.0 - adherence_coeff(hs);
		    alpha_state_velocity(alpha,ansl,dim);
		}
	    	zero_normal_velocity(ansl,n,dim);
	    	g_obstacle_state(ansr,sizest);
	    }
	    break;

	case	DIRICHLET_BOUNDARY:
	    if (is_excluded_comp(lcomp,intfc) == YES)
	    	g_obstacle_state(ansl,sizest);
	    else
	    {
		sten->comp = lcomp;
		sten->states = sten->leftst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansl,fr);
	    }
	    if (is_excluded_comp(rcomp,intfc) == YES)
	    	g_obstacle_state(ansr,sizest);
	    else
	    {
		sten->comp = rcomp;
		sten->states = sten->rightst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansr,fr);
	    }
	    break;

	default:
	    sten->comp = lcomp;
	    sten->states = sten->leftst;
	    one_side_npt_tang_solver(ds,dt,sten,ansl,fr);
	    sten->comp = rcomp;
	    sten->states = sten->rightst;
	    one_side_npt_tang_solver(ds,dt,sten,ansr,fr);
	    break;
	}
        /* Test code */
        if(debugging("local_LF"))
        {
            if(local_switch == YES)
                fr->_one_side_npt_tang_solver =
                    one_side_npt_tang_solver;
        }
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") && 
	    ((is_bad_state(ansl,YES,"g_two_side_npt_tang_solver")) ||
	     (is_bad_state(ansr,YES,"g_two_side_npt_tang_solver"))))
	{
	    int	    i, nrad = sten->npts/2;
	    char    s[80];
	    screen("ERROR in g_two_side_npt_tang_solver(), bad state found\n");
	    (void) printf("ansl - ");
	    fprint_raw_gas_data(stdout,ansl,current_interface()->dim);
	    (void) printf("ansr - ");
	    fprint_raw_gas_data(stdout,ansr,current_interface()->dim);

	    (void) printf("ds = %g, dt = %g\n",ds,dt);
	    print_general_vector("dir = ",sten->dir,dim,"\n");
	    print_Tan_stencil(fr,sten);

	    for (i = -nrad; i <= nrad; ++i)
	    {
	        (void) sprintf(s,"sten->leftst[%d]\n",i);
	        verbose_print_state(s,sten->leftst[i]);
	        (void) sprintf(s,"sten->rightst[%d]\n",i);
	        verbose_print_state(s,sten->rightst[i]);
	    }
	    verbose_print_state("ansl",ansl);
	    verbose_print_state("ansr",ansr);
	    (void) printf("Input hypersurface\n");
	    print_hypersurface(hs);
	    print_interface(intfc);
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
}		/*end g_two_side_npt_tang_solver*/


/*
*			check_ans():
*
*	Used in this file and ggodunov.c only.
*	Returns YES if states are physical,  NO otherwise.
*/

EXPORT	boolean check_ans(
	const char	*function,
	double		ds,
	double		dt,
	Locstate	ans,
	COMPONENT	comp,
	Stencil		*sten,
	int		increment_count)
{
	boolean		bad;
	int		nrad = sten->npts/2;
	int		*icoords = sten->icoords[0];
	int		i, dim = sten->fr->interf->dim;

	bad = is_bad_state(ans,YES,function);
	if (bad == NO)
	{
	    for (i = -nrad; (bad==NO) && i <= nrad; ++i)
	    	bad = is_bad_state(sten->st[i],YES,function);
	}
	if (bad == YES)
	{
	    if (debugging("bad_gas"))
	    {
 	        screen("WARNING - check_ans() detects bad gas state in %s\n",
		       function);
	        print_general_vector("Position = ",Coords(sten->p[0]),dim,"");
		print_int_vector(", icoords = ",icoords,dim,"\n");
		(void) printf("ds = %g, dt = %g, comp = %d\n",ds,dt,comp);
		print_Stencil(sten);
		for (i = -nrad; i <= nrad; ++i)
		{
		    (void) printf("state[%d]\n",i);
		    (*sten->fr->print_state)(sten->st[i]);
		}
		(void) printf("ANSWER\n");
		(*sten->fr->print_state)(ans);
		(void) printf("\n");
	    }

	    if (too_many_bad_gas_states("check_ans()",increment_count,sten->fr))
		clean_up(ERROR);
	    return NO;
	}
	return YES;
}		/*end check_ans*/

/*
*			check_gas():
*
*	Returns YES if states are physical,  NO otherwise.
*/

EXPORT boolean check_gas(
	const char      *function,
	Locstate	*sts,
	Locstate	ans,
	Tan_stencil	*sten,
	int		increment_count,
	Front		*fr)
{
	boolean		bad;
	int		i, j, nrad = sten->npts/2;
	int		dim = fr->rect_grid->dim;
	char		title[20];

	bad = is_bad_state(ans,YES,function);
	if (bad == NO)
	{
	    for (i = -nrad; (bad==NO) && i <= nrad; ++i)
	    	bad = is_bad_state(sts[i],YES,function);
	}
	if (bad == YES)
	{
	    if (debugging("bad_gas"))
	    {
	        (void) printf("WARNING in check_gas(),  bad gas state\n");
	        for (i = -nrad; i <= nrad; ++i)
	        {
	            (void) printf("curve[%d] = %llu, point[%d] = %llu",i,
	    		          curve_number(Curve_of_hs(sten->hs[i])),i,
	                          point_number(sten->p[i]));
	            if (sten->p[i] != NULL)
	            {
	                for (j = 0; j < dim; ++j)
	                    (void) printf(" %g",Coords(sten->p[i])[j]);
	            }
	            (void) printf("\n");
	            (void) sprintf(title,"state[%d]",i);
	            verbose_print_state(title,sts[i]);
	        }
	        verbose_print_state("NEW",ans);
	        (void) printf("\n");
	    }

	    if (too_many_bad_gas_states("check_gas()",increment_count,fr))
	    	clean_up(ERROR);
		
	    return NO;
	}
	return YES;
}		/*end check_gas*/


LOCAL	boolean too_many_bad_gas_states(
	const char	*mesg,
	int		increment_count,
	Front		*front)
{
	static const  int MAX_WARNINGS = 20;/*TOLERANCE*/
	static int total_warnings = 0;
	static int warnings_per_step = 0;
	static int timestep_of_last_warning = -1;

	if (!increment_count)
	    return NO;
	++total_warnings;
	if (front->step != timestep_of_last_warning)
	{
	    timestep_of_last_warning = front->step;
	    warnings_per_step = 0;
	}
	if (++warnings_per_step > MAX_WARNINGS)
	{
	    screen("Fatal ERROR in %s\nERROR - too many (%d) bad gas states "
	           "in a single time step\n",
		   mesg,warnings_per_step);
	    screen("ERROR - Total warnings = %d\n",total_warnings);
	    return YES;
	}
	return NO;
}		/*end too_many_bad_gas_states*/

/*
*			g_load_state_vectors():
*
*	This function loads the conservative state variables from the wave into
*	a Vec_Gas suitable for use by the vector solvers.  This function is
*	currently used only by the vector solvers, thus we make the assumption
*	that imin == 0.	 This means there are offset artificial states 
*	on either end of the Vec_Gas -- 0 to (offset-1) and (imax-offset) to
*	(imax-1).  
*
*	To simply load a Vec_Gas from arbitray imin to arbitray imax,
*	just pass offset = 0.
*/

EXPORT TIME_DYNAMICS g_load_state_vectors(
	int		swp_num,
	int		*iperm,
	Vec_Gas		*vst,
	int		imin,
	int		imax,
	Wave		*wv,
	Wave		*newwv,
	int		*icoords,
	int		pbuf)
{
	Locstate	rst;
	Locstate	*state = vst->state;
	double		**coords = vst->coords;
	double		*rho = vst->rho;
	double		*en_den = vst->en_den;
	double		*m[MAXD];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double		*vacuum_dens = vst->vacuum_dens;
	double		*min_pressure = vst->min_pressure;
	double		*min_energy = vst->min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int		i, j;
	int		dim = wv->rect_grid->dim;
	int		idirs[MAXD];
	static double    **Q = NULL;

	clear_Vec_Gas_set_flags(vst);
	for (j = 0; j < dim; ++j)
	{
	    idirs[j] = iperm[(j+swp_num)%dim];
	    m[j] = vst->m[j];
	}
	if (Q == NULL)
	    bi_array(&Q,3,3,FLOAT);
	Q[0][0] = Q[0][1] = Q[0][2] = 0.0;
	Q[1][0] = Q[1][1] = Q[1][2] = 0.0;
	Q[2][0] = Q[2][1] = Q[2][2] = 0.0;
	for (i = 0; i < dim; ++i)
	    Q[i][idirs[i]] = 1.0;
	for (; i < 3; ++i)
	    Q[i][i] = 1.0;
	vst->Q = (const double* const*)Q;
	for (i = imin; i < imax; ++i)
	{
	    icoords[idirs[0]] = i + pbuf;
	    state[i] = rst = Rect_state(icoords,wv);
	    coords[i] = Rect_coords(icoords,wv);
	    rho[i] = Dens(rst);
	    en_den[i] = Energy(rst);
	    for (j = 0; j < dim; ++j)
		m[j][i] = Mom(rst)[idirs[j]];
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                /* Include partial density into vst */
                Gas_param *params = Params(rst);
                int    num_comps;
                double  *prho;
                double  **rho0 = vst->rho0;
                if((params != NULL) &&
                   ((num_comps = params->n_comps) != 1))
                {
                    prho = pdens(rst);
                    for(j = 0; j < num_comps; j++)
                        rho0[j][i] = prho[j];
                }
            }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (Params(rst) != NULL)
	    {
	    	vacuum_dens[i] = Vacuum_dens(rst);
	    	min_pressure[i] = Min_pressure(rst);
	    	min_energy[i] = Min_energy(rst);
	    }
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	}
	Vec_Gas_field_set(vst,state) = YES;
	Vec_Gas_field_set(vst,coords) = YES;
	Vec_Gas_field_set(vst,rho) = YES;
	Vec_Gas_field_set(vst,en_den) = YES;
	Vec_Gas_field_set(vst,m) = YES;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            double  **rho0 = vst->rho0;
            if((Params(state[imin]) != NULL) &&
               (Params(state[imin])->n_comps != 1))
                Vec_Gas_field_set(vst,rho0) = YES;
        }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Vec_Gas_field_set(vst,min_pressure) = YES;
	Vec_Gas_field_set(vst,min_energy) = YES;
	Vec_Gas_field_set(vst,vacuum_dens) = YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	set_params_jumps(vst,imin,imax);
	if ((Params(state[0]) == NULL) && (vst->nprms <= 1))
	{
	    int	nrad = vsten_radius(wv);
	    copy_states_to_new_time_level(idirs,wv,newwv,imin+nrad,
				          imax-nrad,icoords,pbuf);
	    return CONSTANT_IN_TIME;
	}
	g_init_obstacle_states(vst,imin,dim);
	return DYNAMIC_IN_TIME;
}		/*end g_load_state_vectors*/

/*
*			g_init_obstacle_states():
*
*	This function copies a valid state into any index which is an
*	obstacle state.	 This allows the interior solver to ignore 
*	obstacle states completely.  See assign_wave_state_vectors().
*/

LOCAL void g_init_obstacle_states(
	Vec_Gas		*vst,
	int		imin,
	int		dim)
{
	int		start, end;     /*indices of obstacle state*/
	int		non_obst;     	/*index of state to be copied*/
	int		*prms_jmp = vst->prms_jmp;
	int		i, j, k;

	for (i = 0; i < vst->nprms; ++i)
	{
	    if (is_obstacle_state(vst->state[prms_jmp[i]]))
	    {
	        start = prms_jmp[i];
	        end = prms_jmp[i + 1];

	        non_obst = (start == imin) ? end : start - 1;

	        for (j = start; j < end; ++j)
	        {
	            vst->rho[j] = vst->rho[non_obst];
	            vst->en_den[j] = vst->en_den[non_obst];
	            for (k = 0; k < dim; ++k)
	            	vst->m[k][j] = vst->m[k][non_obst];
	            vst->state[j] = vst->state[non_obst];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	            vst->vacuum_dens[j] = vst->vacuum_dens[non_obst];
	            vst->min_pressure[j] = vst->min_pressure[non_obst];
	            vst->min_energy[j] = vst->min_energy[non_obst];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        }
	    }
	}
}		/*end g_init_obstacle_states*/

EXPORT	void	set_rotation(
	double	    **Q,
	const double *dir,
	int	    dim)
{
	double	dmin, mag;
	int	i, imin;

	for (i = 0; i < dim; ++i)
	    Q[0][i] = dir[i];
	switch (dim)
	{
	case 1:
	    break;
	case 2:
	    Q[1][0] = -Q[0][1];
	    Q[1][1] =  Q[0][0];
	    break;
	case 3:
	    dmin = fabs(dir[dim-1]);	imin = dim-1;
	    for (i = 0; i < (dim-1); ++i)
	    {
	    	if (fabs(dir[i]) < dmin)
	    	{
	    	    dmin = fabs(dir[i]);
	    	    imin = i;
	    	}
	    }
	    Q[1][imin] = 0.0;
	    Q[1][(imin+1)%3] = -Q[0][(imin+2)%3];
	    Q[1][(imin+2)%3] =  Q[0][(imin+1)%3];
	    mag = mag_vector(Q[1],dim);
	    Q[1][0] /= mag; Q[1][1] /= mag; Q[1][2] /= mag;
	    mag = vector_product(Q[0],Q[1],Q[2],dim);
	    Q[2][0] /= mag; Q[2][1] /= mag; Q[2][2] /= mag;
	    break;
	}
}		/*end set_rotation*/

EXPORT	boolean	g_detect_and_load_mix_state(
	int	idir,
	Stencil *sten,
	int	indx)
{
	COMPONENT  *comp = sten->comp;
	if(g_composition_type() != MULTI_COMP_NON_REACTIVE)
	    return NO;
	if (detect_ambient_vapor_mix(idir,sten) == YES)
	{
	    if(comp[indx] == VAPOR_COMP || comp[indx] == AMBIENT_COMP)
	    {
	        ft_assign(sten->worksp[indx],Rect_state(sten->icoords[indx],sten->wave),sten->fr->sizest);
	        sten->st[indx] = sten->worksp[indx];
	        return YES;
	    }
	}
	return NO;
}	/* end g_detect_and_load_mix_state */

LOCAL	int	detect_ambient_vapor_mix(
	int	idir,
	Stencil *sten)
{
        int     i;
        int     find_2nd = NO;
	int	endpt = stencil_radius(sten->wave);
        if(sten->newcomp == AMBIENT_COMP &&
           sten->comp[0] == VAPOR_COMP)
        {
            for(i = -endpt; i <= endpt; ++i)
            {
                if(sten->newcomp == AMBIENT_COMP && sten->comp[i] == VAPOR_COMP &&
                   Find_rect_comp(idir,sten->icoords[i],sten->newwave)
                   == AMBIENT_COMP)
                {
                    find_2nd = YES;
                    return YES;
                }
            }
        }
        else
            return NO;	
}	/* end detect_ambient_vapor_mix */

LOCAL int local_LF_npt_tang_solver_switch(
        double        ds,
        Tan_stencil  *sten,
        Front        *fr)
{
        COMPONENT  ncp, pcp;
        CURVE      *c = NULL;
        Locstate   *lsts = sten->leftst;
        Locstate   *rsts = sten->rightst;
        double      Tl, Tr;
        int        i, nrad = sten->npts/2;

        if(sten->newhs == NULL)
            return NO;
        if(wave_type(sten->newhs) <= FIRST_USER_BOUNDARY_TYPE)
            return NO;
        ncp = negative_component(sten->newhs);
        pcp = positive_component(sten->newhs);
        c = Curve_of_hs(sten->newhs);

        /* Size effact */
        if (c->num_points < 30)
        {
            on_what = 1;
            return YES;
        }

        if(wave_type(sten->newhs) == NEUMANN_BOUNDARY ||
           wave_type(sten->newhs) == DIRICHLET_BOUNDARY)
        {
            if (is_excluded_comp(ncp,fr->interf) == YES)    
            {
                /* Temperature effact */
                Tr = temperature(rsts[-nrad]);
                for (i = -nrad+1; i <= nrad; ++i)
                {
                    Tl = temperature(rsts[i]);
                    if ((Tr > Tl*1.2) ||
                        (Tl > Tr*1.2))
                    {
                        on_what = 2;
                        return YES;
                    }
                }
            }
            else
            {
                /* Temperature effact */
                Tl = temperature(lsts[-nrad]);
                for (i = -nrad+1; i <= nrad; ++i)
                {
                    Tr = temperature(lsts[i]);
                    if ((Tr > Tl*1.2) ||
                        (Tl > Tr*1.2))
                    {
                        on_what = 2;
                        return YES;
                    }
                }
            }
            return NO;
        }

        /* Size effact */
        if (c->num_points < 30)
        {
            on_what = 1;
            return YES;
        }
        /* Temperature effact */
        for (i = -nrad; i <= nrad; ++i)
        {
            Tl = temperature(lsts[i]);
            Tr = temperature(rsts[i]);
            if ((Tr > Tl*1.2) ||
                (Tl > Tr*1.2))
            {
                on_what = 2;
                return YES;
            }
        }
        /* Curvature effact */
        if(fabs(1.0/sten->curvature/ds) < 2.0)
        {
            on_what = 3;
            return YES;
        }
        return NO;
}	/* end local_LF_npt_tang_solver_switch */
