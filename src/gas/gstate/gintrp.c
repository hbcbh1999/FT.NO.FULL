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
*				gintrp.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines for the interpolation of state data
*
*/


#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
#if defined(TWOD)
LOCAL	boolean	correct_bad_params_in_quad(BILINEAR_ELEMENT*,TRI_SOLN*);
#endif /* defined(TWOD) */

LOCAL	void	cyl_lin_rgas_weighted_coefs(double*,double*,double,double,
					    RECT_GRID*);
LOCAL	void	cyl_lin_vol_rgas_weighted_coefs(double*,double*,double,double,
						RECT_GRID*);
LOCAL	void	cyl_lin_vol_weighted_coefs(double*,double*,double,double,
					   RECT_GRID*);
LOCAL	void	lin_unweighted_coefs(double*,double*,double,double,RECT_GRID*);
LOCAL	void	sph_lin_rgas_weighted_coefs(double*,double*,double,double,
					    RECT_GRID*);
LOCAL	void	sph_lin_vol_rgas_weighted_coefs(double*,double*,double,double,
						RECT_GRID*);
LOCAL	void	sph_lin_vol_weighted_coefs(double*,double*,double,double,
					   RECT_GRID*);


LOCAL	void	(*lin_weighted_coefs)(double*,double*,double,double,RECT_GRID*) =
						       lin_unweighted_coefs;
#if defined(TWOD) || defined(THREED)
LOCAL	void	cyl_tri_rgas_weighted_coefs(double*,double*,double*,double*,
					    double*,double*,RECT_GRID*);
LOCAL	void	cyl_tri_vol_weighted_coefs(double*,double*,double*,double*,
					   double*,double*,RECT_GRID*);
LOCAL	void	cyl_tri_vol_rgas_weighted_coefs(double*,double*,double*,double*,
						double*,double*,RECT_GRID*);
LOCAL	void	tri_unweighted_coefs(double*,double*,double*,double*,double*,double*,
				     RECT_GRID*);

LOCAL	void	(*tri_weighted_coefs)(double*,double*,double*,double*,
				      double*,double*,RECT_GRID*) =
				      tri_unweighted_coefs;
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
LOCAL	void	cyl_quad_rgas_weighted_coefs(double*,double*,double*,double*,
					     BILINEAR_ELEMENT*,RECT_GRID*);
LOCAL	void	cyl_quad_vol_rgas_weighted_coefs(double*,double*,double*,double*,
						 BILINEAR_ELEMENT*,RECT_GRID*);
LOCAL	void	cyl_quad_vol_weighted_coefs(double*,double*,double*,double*,
					    BILINEAR_ELEMENT*,RECT_GRID*);
LOCAL	void	quad_unweighted_coefs(double*,double*,double*,double*,
				      BILINEAR_ELEMENT*,RECT_GRID*);

LOCAL	void	(*quad_weighted_coefs)(double*,double*,double*,double*,
				       BILINEAR_ELEMENT*,RECT_GRID*) =
						       quad_unweighted_coefs;
#endif /* defined(TWOD) */
LOCAL  void 	least_sqr_fit_old(double*,double*,double*,int,int);
LOCAL  void 	least_sqr_fit(double*,double*,double*,int,int);

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
FORTRAN     void    FORTRAN_NAME(dqrdc)(double*,int*,int*,int*,double*,
				int*,double*,int*);
FORTRAN     void    FORTRAN_NAME(dqrsl)(double*,int*,int*,int*,double*,double*,
				double*,double*,double*,double*,double*,int*,int*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */


/*
*			INTERPOLATION ROUTINES

*
*		g_linear_combination_of_states()
*		g_lin_comb_states()
*		gt_lin_comb_states()
*		g_tri_lin_comb_states()
*		gt_tri_lin_comb_states()
*		g_tri_interpolator()
*		gt_tri_interpolator()
*		g_quad_interpolator()
*		gt_quad_interpolator()
*		g_grad_tri_interpolator()
*		g_grad_quad_interpolator()
*		g_tetra_interpolator()
*		g_cube_interpolator()
*
*	The following routines are accessed through function pointers
*	set in init_physics():
*	The lin_comb_state function is set to a function pointer in
*	the Front structure,  while the other functions are set to
*	function pointers in the Wave structure.  The difference between
*	the g_ and the gt_ versions is that the g_ do linear or bilinear
*	interpolation on the conserved quantities directly,  while the gt_
*	versions are linear in the log of the density and temperature
*	and the state velocity.
*
*	
*/

/*
*			g_linear_combination_of_states():
*
*	Simple unadorned function for forming a linear combination of
*	a group of states.  No checking is performed for vacumm states
*	and no corrections to the interpolation coefficients are
*	made to account for geometry.  This function is intended as a
*	low level utility.
*/

EXPORT void g_linear_combination_of_states(
	const double	*alpha,		/* coeffs in linear combination */
	const Locstate	*state,		/* array of states */
	const int	n,		/* number of states */
	Locstate	ans)
{
	int		dim;
	int		i, j;


	for (i = 1; i < n; i++)
	{
	    if (Different_params(state[i],state[0]))
	    {
		INTERFACE *intfc = current_interface();
		int	dim = intfc->dim;
	    	int	j;

	    	screen("ERROR in g_linear_combination_of_states(), "
	    	       "different params\n");
		(void) printf("number of states = %d\n",n);
		(void) printf("state = %p\n",state);
	    	for (j = 0; j < n; j++)
	    	{
		    (void) printf("state[%d] = %p\n",j,state[j]);
	    	    fprint_raw_gas_data(stdout,state[j],dim);
	    	}
	    	print_interface(intfc);
	    	clean_up(ERROR);
	    }
	    if (state_type(state[i]) != state_type(state[0]))
	    {
	    	char	s[80];
	    	int	j;

	    	screen("ERROR in g_linear_combination_of_states(), "
	    	       "different state types\n");
	    	for (j = 0; j < n; j++)
	    	{
	    	    (void)sprintf(s,"state[%d]",j);
	    	    verbose_print_state(s,state[j]);
	    	}
	    	print_interface(current_interface());
	    	clean_up(ERROR);
	    }
	}
	if (is_obstacle_state(state[0]))
	{
	    g_obstacle_state(ans,g_sizest());
	    return;
	}

	dim = Params(state[0])->dim;

	(*Params(state[0])->_clear_state)(ans,Params(state[0])->sizest);
	Set_params(ans,state[0]);
	set_type_of_state(ans,state_type(state[0]));
	for (i = 0; i < n; i++)
	{
	    Dens(ans) += alpha[i] * Dens(state[i]);
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                int   k;
                double sum = 0.0;
                /* pdens(ans) is consistent with Dens(ans) */
                if(Params(ans)->n_comps != 1)
                {
                    if(Params(ans)->n_comps != Params(state[i])->n_comps)
                    {
                        printf("ERROR: g_linear_combination_of_states\n");
                        printf("MULTI_COMP_NON_REACTIVE, sl and s2 have different # of comps.\n");
                        clean_up(ERROR);
                    }
                    for(k = 0; k < Params(ans)->n_comps; k++)
                    {
                        pdens(ans)[k] += alpha[i] * pdens(state[i])[k];
                        /* New 051005 */
                        if(fabs(pdens(ans)[k]) < 10.0*MACH_EPS && pdens(ans)[k] < 0.0)
                            pdens(ans)[k] = 0.0;
                        else if(fabs(pdens(ans)[k]) > 10.0*MACH_EPS && pdens(ans)[k] < 0.0)
                        {
                            printf("ERROR in g_linear_combination_of_states()\n");
                            printf("partial density < 0.0\n");
                            clean_up(ERROR);
                        }
                        /* End of New 051005 */
                    }
                    /* New 051005 */
                    for(k = 0; k < Params(ans)->n_comps; k++)
                        sum += pdens(ans)[k];
                    for(k = 0; k < Params(ans)->n_comps; k++)
                        pdens(ans)[k] = Dens(ans)/sum*pdens(ans)[k];
                    /* End of New 051005 */
                }
            }
	    Energy(ans) += alpha[i] * Energy(state[i]);
	    for (j = 0; j < dim; j++)
	    	Mom(ans)[j] += alpha[i] * Mom(state[i])[j];

#if defined(COMBUSTION_CODE)
	    switch(Composition_type(ans))
	    {
	    case ZND:
	    	Prod(ans) += alpha[i] * Prod(state[i]);
	    	break;
	    case TWO_CONSTITUENT_REACTIVE:
	    	Prod(ans) += alpha[i] * Prod(state[i]);
	    	Dens1(ans) += alpha[i] * Dens1(state[i]);
	    	break;
	    case PURE_NON_REACTIVE:
	    case PTFLAME:
	    case THINFLAME:
	    default:
	    	break;
	    }
#endif /* defined(COMBUSTION_CODE) */
	    reset_gamma(ans);
	}
}		/*end g_linear_combination_of_states*/


/*
*			g_lin_comb_states():
*			gt_lin_comb_states():
*
*	Forms a linear combination of two local states.  This function
*	is used primarily to calculate an interpolated state between the
*	two local states.  It can be (and is) used for general linear
*	combinations as well (e.g. finite difference approximations to
*	derivatives).  
*
*	With cylindrical coordinates,
*	the purpose is ASSUMED to be linear interpolation if the sum of the
*	coefficients is 1.  In this case, the interpolation coefficients
*	may be modified to reflect weights which are better for
*	cylindrical geometry.  (The modification depends on the
*	..._INTERP compilation flags.) 
*/


/*ARGSUSED*/
EXPORT void g_lin_comb_states(
	double		alpha,		/* coeffs in linear combination */
	double		beta,		/* coeffs in linear combination */
	double		*crds1,		/* location of s1		*/
	Locstate	s1,
	double		*crds2,		/* location of s1              */
	Locstate	s2,
	RECT_GRID	*gr,
	Locstate	ans)
{
	int		dim;
	int		i;

	switch (consistent_params_in_lin_comb(s1,s2))
	{
	case PARAMS_CONSISTENT:
	    break;
	case PARAMS_ALL_OBSTACLE_STATES:
	    g_obstacle_state(ans,g_sizest());
	    return;
	case PARAMS_INCONSISTENT:
	    inconsistent_params_in_lin_comb_message(alpha,beta,crds1,s1,
				                    crds2,s2,
						    "g_lin_comb_states",YES);
	    return;
	case OBSTACLE_STATE_FOUND:
	    obstacle_state_found_in_lin_comb_message(alpha,beta,
				                     crds1,s1,crds2,s2,
				                     "g_lin_comb_states",YES);
	    return;
	}

	if (alpha >= beta)
	    Set_params(ans,s1);
	else
	    Set_params(ans,s2);
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	if (is_rotational_symmetry())
	    (*lin_weighted_coefs)(&alpha,&beta,crds1[0],crds2[0],gr);

	Dens(ans) = alpha * Dens(s1) + beta * Dens(s2);
	for (i = 0; i < dim; i++)
	    Mom(ans)[i] = alpha * Mom(s1)[i] + beta * Mom(s2)[i];
	Energy(ans) = alpha * Energy(s1) + beta * Energy(s2);

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int   k;
            double sum = 0.0;
            /* pdens(ans) is consistent with Dens(ans) */
            if(Params(ans)->n_comps != 1)
            {
                if(Params(s1)->n_comps != Params(s2)->n_comps)
                {
                    printf("ERROR: g_lin_comb_states\n");
                    printf("MULTI_COMP_NON_REACTIVE, sl and s2 have different # of comps.\n");
                    clean_up(ERROR);
                }
                for(k = 0; k < Params(ans)->n_comps; k++)
                {
                    pdens(ans)[k] = alpha * pdens(s1)[k] + beta * pdens(s2)[k];
                    /* New 051005 */
                    if(fabs(pdens(ans)[k]) < 10.0*MACH_EPS && pdens(ans)[k] < 0.0)
                        pdens(ans)[k] = 0.0;
                    else if(fabs(pdens(ans)[k]) > 10.0*MACH_EPS && pdens(ans)[k] < 0.0)
                    {
                        printf("ERROR in g_lin_comb_states()\n");
                        printf("partial density %20.18g < 0.0\n",pdens(ans)[k]);
                        print_general_vector("crds1",crds1,dim,"\n");
                        verbose_print_state("s1", s1);
                        print_general_vector("crds2",crds2,dim,"\n");
                        verbose_print_state("s2", s2);
                        clean_up(ERROR);
                    }
                    /* End of New 051005 */
                }
                /* New 051005 */
                for(k = 0; k < Params(ans)->n_comps; k++)
                    sum += pdens(ans)[k];
                for(k = 0; k < Params(ans)->n_comps; k++)
                    pdens(ans)[k] = pdens(ans)[k]/sum*Dens(ans);
                /* End of New 051005 */
            }
        }
#if defined(COMBUSTION_CODE)
	switch(Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = alpha * Prod(s1) + beta * Prod(s2);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = alpha * Prod(s1) + beta * Prod(s2);
	    Dens1(ans) = alpha * Dens1(s1) + beta * Dens1(s2);
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
}		/*end g_lin_comb_states*/

/*ARGSUSED*/
EXPORT void gt_lin_comb_states(
	double		alpha,
	double		beta,
	double		*crds1,
	Locstate	s1,
	double		*crds2,
	Locstate	s2,
	RECT_GRID	*gr,
	Locstate	ans)
{
	int		i;
	int		dim;
	double		T1, T2;
	double		ra;
	static Locstate temp_ans = NULL;

	switch (consistent_params_in_lin_comb(s1,s2))
	{
	case PARAMS_CONSISTENT:
	    break;
	case PARAMS_ALL_OBSTACLE_STATES:
	    g_obstacle_state(ans,g_sizest());
	    return;
	case PARAMS_INCONSISTENT:
	    inconsistent_params_in_lin_comb_message(alpha,beta,crds1,s1,
				                    crds2,s2,
						    "gt_lin_comb_states",YES);
	    return;
	case OBSTACLE_STATE_FOUND:
	    obstacle_state_found_in_lin_comb_message(alpha,beta,
						     crds1,s1,crds2,s2,
				                     "gt_lin_comb_states",YES);
	    return;
	}

	if (temp_ans == NULL)
	{
	    (*Params(s1)->_alloc_state)(&temp_ans,Params(s1)->sizest);
	    set_type_of_state(temp_ans,FGAS_STATE);
	}

	if (alpha >= beta)
	    Set_params(ans,s1);
	else
	    Set_params(ans,s2);
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	if (is_rotational_symmetry())
	    (*lin_weighted_coefs)(&alpha,&beta,crds1[0],crds2[0],gr);

	Dens(ans) = ra = pow(Dens(s1),alpha)*pow(Dens(s2),beta);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            double sum;
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = pow(pdens(s1)[i],alpha)*pow(pdens(s2)[i],beta);
                sum = 0.0;
                for(i = 0; i < Params(ans)->n_comps; i++)
                    sum += pdens(ans)[i];
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = pdens(ans)[i]/sum*Dens(ans);
            }
        }
	for (i = 0; i < dim; i++)
	{
	    Mom(ans)[i] = ra*(alpha*vel(i,s1) + beta*vel(i,s2));
	}

	Dens(temp_ans) = Dens(ans);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(temp_ans)[i] = pdens(ans)[i];
            }
        }
	T1 = temperature(s1);	T2 = temperature(s2);
	Temperature(temp_ans) = pow(T1,alpha)*pow(T2,beta);
	Set_params(temp_ans,ans);
	reset_gamma(temp_ans);
	Energy(ans) = internal_energy(temp_ans) + kinetic_energy(ans);

#if defined(COMBUSTION_CODE)
	switch(Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = ra * (alpha*react(s1) + beta*react(s2));
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = ra * (alpha*react(s1) + beta*react(s2));
	    Dens1(ans) = pow(Dens1(s1),alpha)*pow(Dens1(s2),beta);
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
}		/*end gt_lin_comb_states*/

EXPORT CONSISTENT_PARAMS consistent_params_in_lin_comb(
	Locstate	s1,
	Locstate	s2)
{
	if (is_obstacle_state(s1) && is_obstacle_state(s2))
	    return PARAMS_ALL_OBSTACLE_STATES;
	if (is_obstacle_state(s1) || is_obstacle_state(s2))
	    return OBSTACLE_STATE_FOUND;

#if !defined(COMBUSTION_CODE)
	if (Different_params(s1,s2))
	    return PARAMS_INCONSISTENT;
#endif /* !defined(COMBUSTION_CODE) */
	return PARAMS_CONSISTENT;
}		/*end consistent_params_in_lin_comb*/

EXPORT	void obstacle_state_found_in_lin_comb_message(
	double      alpha,
	double      beta,
	double      *crds1,
	Locstate   s1,
	double      *crds2,
	Locstate   s2,
	const char *fname,
	boolean    fatal)
{
	static const char *mesg = "single obstacle state";
	int        dim;

	if (fatal == YES)
	    screen("ERROR in %s(), %s\n",fname,mesg);
	else
	    (void) printf("WARNING in %s(), %s\n",fname,mesg);
	dim = (is_obstacle_state(s1)) ? Params(s2)->dim :
				        Params(s1)->dim;
	(void) printf("alpha = %g, beta = %g\n",alpha,beta);
	print_general_vector("crds1 = ",crds1,dim,"");
	verbose_print_state("s1",s1);
	print_general_vector("crds2 = ",crds2,dim,"");
	verbose_print_state("s2",s2);
	if (fatal == YES)
	{
	    print_interface(current_interface());
	    show_intfc_states(current_interface());
	    clean_up(ERROR);
	}
}	/*end obstacle_state_found_in_lin_comb_message*/


EXPORT	void inconsistent_params_in_lin_comb_message(
	double      alpha,
	double      beta,
	double      *crds1,
	Locstate   s1,
	double      *crds2,
	Locstate   s2,
	const char *fname,
	boolean    fatal)
{
	static const char *mesg = "different params";
	int        dim;

	if (fatal == YES)
	    screen("ERROR in %s(), %s\n",fname,mesg);
	else
	    (void) printf("WARNING in %s(), %s\n",fname,mesg);
	(void) printf("params1 = %llu params2 = %llu\n",
		      gas_param_number(Params(s1)),
		      gas_param_number(Params(s2)));
	(void) printf("alpha = %g, beta = %g\n",alpha,beta);
	dim = (is_obstacle_state(s1)) ? Params(s1)->dim :
				        Params(s2)->dim;
	print_general_vector("crds1 = ",crds1,dim,"");
	verbose_print_state("s1",s1);
	print_general_vector("crds2 = ",crds2,dim,"");
	verbose_print_state("s2",s2);
	if (fatal == YES)
	{
	    print_interface(current_interface());
	    show_intfc_states(current_interface());
	    clean_up(ERROR);
	}
}	/*end inconsistent_params_in_lin_comb_message*/

#if defined(TWOD) || defined(THREED)

EXPORT CONSISTENT_PARAMS consistent_params_in_tri_lin_comb(
	Locstate s0,
	Locstate s1,
	Locstate s2)
{
	if (is_obstacle_state(s0) &&
	    is_obstacle_state(s1) &&
	    is_obstacle_state(s2)) 
	{
	    return PARAMS_ALL_OBSTACLE_STATES;
	}
#if !defined(COMBUSTION_CODE)
	if (Different_params(s0,s1) || Different_params(s0,s2))
	    return PARAMS_INCONSISTENT;
#endif /* !defined(COMBUSTION_CODE) */
	return PARAMS_CONSISTENT;
}		/*end consistent_params_in_tri_lin_comb*/

EXPORT void inconsistent_params_in_tri_lin_comb_mesg(
	double	 f0,
	double	 f1,
	double	 f2,
	double	 *crds0,
	Locstate s0,
	double	 *crds1,
	Locstate s1,
	double	 *crds2,
	Locstate s2,
	const char *fname,
	boolean    fatal)
{
	static const char *mesg = "different params";
	int dim;
	if (fatal == YES)
	    screen("ERROR in %s(), %s\n",fname,mesg);
	else
	    (void) printf("WARNING in %s(), %s\n",fname,mesg);
	(void) printf("params0 = %llu params1 = %llu params2 = %llu\n",
		      gas_param_number(Params(s0)),
		      gas_param_number(Params(s1)),
		      gas_param_number(Params(s2)));
	if (!is_obstacle_state(s0))
	    dim = Params(s0)->dim;
	else if (!is_obstacle_state(s1))
	    dim = Params(s1)->dim;
	else
	    dim = Params(s2)->dim;
	(void) printf("f0 = %g, f1 = %g, f2 = %g\n",f0,f1,f2);
	print_general_vector("crds0 = ",crds0,dim,"");
	verbose_print_state("s0",s0);
	print_general_vector("crds1 = ",crds1,dim,"");
	verbose_print_state("s1",s1);
	print_general_vector("crds2 = ",crds2,dim,"");
	verbose_print_state("s2",s2);
	if (fatal == YES)
	{
	    print_interface(current_interface());
	    show_intfc_states(current_interface());
	    clean_up(ERROR);
	}
}		/*end inconsistent_params_in_tri_lin_comb_mesg*/

/*
*			g_tri_interpolator():
*
*	Forms a linear combination of three local states.  This function
*	is generally used to calculate an interpolated state from the three
*	states at the vertices of a triangle.
*
*	With cylindrical coordinates,
*	the purpose is ASSUMED to be interpolation if the sum of the
*	coefficients is 1.  In this case, the interpolation coefficients
*	may be modified to reflect weights which are better for cylindrical
*	geometry.  (The modification depends on the ..._INTERP 
*	compilation flags.) 
*/


/*ARGSUSED*/
EXPORT boolean g_tri_interpolator(
	double		*f,
	LINEAR_ELEMENT	*t,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	double		f0 = f[0], f1 = f[1], f2 = f[2];
	Locstate	s0 = t->s[0], s1 = t->s[1], s2 = t->s[2];
					/* states in linear combination */

	return g_tri_lin_comb_states(f0,f1,f2,Coords(t->p[0]),s0,
		                     Coords(t->p[1]),s1,Coords(t->p[2]),s2,
		                     &soln->tri_grid->comp_grid,ans);
}		/*end g_tri_interpolator*/

/*ARGSUSED*/
EXPORT boolean gt_tri_interpolator(
	double		*f,
	LINEAR_ELEMENT	*t,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	double		f0 = f[0], f1 = f[1], f2 = f[2];
	Locstate	s0 = t->s[0], s1 = t->s[1], s2 = t->s[2];

	return gt_tri_lin_comb_states(f0,f1,f2,Coords(t->p[0]),s0,
		                      Coords(t->p[1]),s1,Coords(t->p[2]),
				      s2,&soln->tri_grid->comp_grid,ans);
}		/*end gt_tri_interpolator*/

/*
*			g_tri_lin_comb_states():
*			gt_tri_lin_comb_states():
*
*	Forms a linear combination of three states.  This function
*	is used primarily to calculate an interpolated state between the
*	three local states.  It can be (and is) used for general linear
*	combinations as well (e.g. finite difference approximations to
*	derivatives).  
*
*	With cylindrical coordinates,
*	the purpose is ASSUMED to be linear interpolation if the sum of
*	the coefficients is 1.  In this case, the interpolation coefficients
*	may be modified to reflect weights which are better for
*	cylindrical geometry.  (The modification depends on the
*	..._INTERP compilation flags.) 
*/

/*ARGSUSED*/
EXPORT boolean g_tri_lin_comb_states(
	double		f0,
	double		f1,
	double		f2,
	double		*crds0,
	Locstate	s0,
	double		*crds1,
	Locstate	s1,
	double		*crds2,
	Locstate	s2,
	RECT_GRID	*gr,
	Locstate	ans)
{
	int		i, dim;

	switch (consistent_params_in_tri_lin_comb(s0,s1,s2))
	{
	case PARAMS_CONSISTENT:
	    break;
	case PARAMS_ALL_OBSTACLE_STATES:
	    g_obstacle_state(ans,g_sizest());
	    return FUNCTION_SUCCEEDED;
	case PARAMS_INCONSISTENT:
	case OBSTACLE_STATE_FOUND:
	    inconsistent_params_in_tri_lin_comb_mesg(f0,f1,f2,crds0,s0,crds1,s1,
					             crds2,s2,
						     "g_tri_lin_comb_states",
						     NO);
	    return FUNCTION_FAILED;
	}
	Set_params(ans,s0);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	if (is_rotational_symmetry())
	    (*tri_weighted_coefs)(&f0,&f1,&f2,crds0,crds1,crds2,gr);

	Dens(ans) = f0 * Dens(s0) + f1 * Dens(s1) + f2 * Dens(s2);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ans)->n_comps != 1)
            {
                double sum;
                for(i = 0; i < Params(ans)->n_comps; i++)
                {
                    pdens(ans)[i] = f0*pdens(s0)[i] + f1*pdens(s1)[i]
                                  + f2*pdens(s2)[i];
                    /* New 051005 */
                    if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                        pdens(ans)[i] = 0.0;
                    else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                    {
                        printf("ERROR in g_tri_lin_comb_states()\n");
                        printf("partial density < 0.0\n");
                        clean_up(ERROR);
                    }
                    /* End of New 051005 */
                }
                /* TMP */
                sum = 0.0;
                for(i = 0; i < Params(ans)->n_comps; i++)
                    sum += pdens(ans)[i];
                /* New 051005 */
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = Dens(ans)/sum*pdens(ans)[i];
                /* End of New 051005 */

                if(fabs(sum - Dens(ans)) > 0.0000001)
                {
                    /* TMP */
                    int icrds[MAXD];
                    printf("ERROR g_tri_lin_comb_states()\n"); 
                    printf("partial density not consistent\n"); 
                    verbose_print_state("ans", ans); 
                    print_general_vector("crds0",crds0,dim,"\n");
                    verbose_print_state("s0", s0); 
                    print_general_vector("crds1",crds1,dim,"\n");
                    verbose_print_state("s1", s1); 
                    print_general_vector("crds2",crds2,dim,"\n");
                    verbose_print_state("s2", s2); 
                    rect_in_which(crds2,icrds,gr);
                    printf("Cell center for crds2[%d %d] h[%g %g]\n",
                      cell_center(icrds[0],0,gr), cell_center(icrds[1],1,gr),
                      gr->h[0], gr->h[1]);
                    clean_up(ERROR);
                }
            }
        }
	Energy(ans) = f0 * Energy(s0) + f1 * Energy(s1) + f2 * Energy(s2);
	for (i = 0; i < dim; i++)
	    Mom(ans)[i] = f0 * Mom(s0)[i] + f1 * Mom(s1)[i] + f2 * Mom(s2)[i];

#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = f0 * Prod(s0) + f1 * Prod(s1) + f2 * Prod(s2);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = f0 * Prod(s0) + f1 * Prod(s1) + f2 * Prod(s2);
	    Dens1(ans) = f0 * Dens1(s0) + f1 * Dens1(s1) + f2 * Dens1(s2);
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
	return FUNCTION_SUCCEEDED;
}		/*end g_tri_lin_comb_states*/

/*ARGSUSED*/
EXPORT boolean gt_tri_lin_comb_states(
	double		f0,
	double		f1,
	double		f2,
	double		*crds0,
	Locstate	s0,
	double		*crds1,
	Locstate	s1,
	double		*crds2,
	Locstate	s2,
	RECT_GRID	*gr,
	Locstate	ans)
{
	double		T0, T1, T2;
	int		i, dim;
	static Locstate temp_ans = NULL;

	switch (consistent_params_in_tri_lin_comb(s0,s1,s2))
	{
	case PARAMS_CONSISTENT:
	    break;
	case PARAMS_ALL_OBSTACLE_STATES:
	    g_obstacle_state(ans,g_sizest());
	    return FUNCTION_SUCCEEDED;
	case PARAMS_INCONSISTENT:
	case OBSTACLE_STATE_FOUND:
	    inconsistent_params_in_tri_lin_comb_mesg(f0,f1,f2,crds0,s0,crds1,s1,
					             crds2,s2,
						     "gt_tri_lin_comb_states",
						     NO);
	    return FUNCTION_FAILED;
	}

	if (temp_ans == NULL)
	{
	    (*Params(s1)->_alloc_state)(&temp_ans,Params(s1)->sizest);
	    set_type_of_state(temp_ans,FGAS_STATE);
	}

	Set_params(ans,s0);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	if (is_rotational_symmetry())
	    (*tri_weighted_coefs)(&f0,&f1,&f2,crds0,crds1,crds2,gr);

	Dens(ans) = pow(Dens(s0),f0)*pow(Dens(s1),f1)*pow(Dens(s2),f2);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            double sum;
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = pow(pdens(s0)[i],f0)*pow(pdens(s1)[i],f1)*
                                    pow(pdens(s2)[i],f2);
                sum = 0.0;
                for(i = 0; i < Params(ans)->n_comps; i++)
                    sum += pdens(ans)[i];
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = pdens(ans)[i]/sum*Dens(ans);
            }
        }
	for (i = 0; i < dim; i++)
	{
	    Mom(ans)[i] = Dens(ans)*(f0 * vel(i,s0) + f1 * vel(i,s1) + 
				     f2 * vel(i,s2));
	}

	Dens(temp_ans) = Dens(ans);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(temp_ans)[i] = pdens(ans)[i];
            }
        }
	T0 = temperature(s0);
	T1 = temperature(s1);
	T2 = temperature(s2);
	Temperature(temp_ans) = pow(T0,f0)*pow(T1,f1)*pow(T2,f2);
	Set_params(temp_ans,ans);
	reset_gamma(temp_ans);
	Energy(ans) = internal_energy(temp_ans) + kinetic_energy(ans);

#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = Dens(ans)*(f0*react(s0) + f1*react(s1) + f2*react(s2));
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = Dens(ans)*(f0*react(s0) + f1*react(s1) + f2*react(s2));
	    Dens1(ans) = pow(Dens1(s0),f0)*pow(Dens1(s1),f1)*pow(Dens1(s2),f2);
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
	return FUNCTION_SUCCEEDED;
}		/*end gt_tri_lin_comb_states*/

#endif /* defined(TWOD) || defined(THREED) */


#if defined(TWOD)

/*
*			g_quad_interpolator():
*
*	Forms a bilinear combination of four local states.  This function
*	is generally used to calculate an interpolated state from the four 
*	states at the vertices of a quadrangle.
*
*	With cylindrical coordinates,
*	the interpolation coefficients may be modified to 
*	reflect	weights which are better for cylindrical geometry.  (The 
*	modification depends on the ..._INTERP compilation flags.) 
*/


EXPORT void g_quad_interpolator(
	double		*f,
	BILINEAR_ELEMENT *q,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	Locstate	sts[(1<<MAXD)];
	Locstate	s0, s1, s2, s3;
	double		fx = f[0], fy = f[1];
	double		f0,f1,f2,f3;
	int		dim, i;

	states_on_bilinear_element(sts,q,soln->tri_grid);
	s0 = sts[0];	s1 = sts[1];	s2 = sts[2];	s3 = sts[3];

	if (is_obstacle_state(s0) && is_obstacle_state(s1)
		 && is_obstacle_state(s2) && is_obstacle_state(s3)) 
	{
	    g_obstacle_state(ans,g_sizest());
	    return;
	}
	if (	(	Different_params(s0,s1) ||
			Different_params(s0,s2) ||
			Different_params(s0,s3)	)
 	&&
		correct_bad_params_in_quad(q,soln) == NO)
	{
	    screen("ERROR in g_quad_interpolator(), different params\n");
	    (void) printf("params0 = %llu params1 = %llu ",
			  gas_param_number(Params(s0)),
			  gas_param_number(Params(s1)));
	    (void) printf("params2 = %llu params3 = %llu\n",
		          gas_param_number(Params(s2)),
			  gas_param_number(Params(s3)));
	    (void) printf("fx = %g, fy = %g\n",fx,fy);
	    print_BILINEAR_ELEMENT(q,soln->tri_grid);
	    verbose_print_state("s0",s0);
	    verbose_print_state("s1",s1);
	    verbose_print_state("s2",s2);
	    verbose_print_state("s3",s3);
	    print_interface(current_interface());
	    clean_up(ERROR);
	}
	Set_params(ans,s0);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;
	
	f1 = fx * (1.0 - fy);
	f2 = fx * fy;
	f3 = fy * (1.0 - fx);
	f0 = 1.0 - f1 - f2 - f3;

	if (is_rotational_symmetry())
	    (*quad_weighted_coefs)(&f0,&f1,&f2,&f3,q,&soln->tri_grid->comp_grid);

	Dens(ans) = f0*Dens(s0) + f1*Dens(s1) + f2*Dens(s2) + f3*Dens(s3);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ans)->n_comps != 1)
            {
                double sum;
                for(i = 0; i < Params(ans)->n_comps; i++)
                {
                    pdens(ans)[i] = f0*pdens(s0)[i] + f1*pdens(s1)[i] 
                                  + f2*pdens(s2)[i] + f3*pdens(s3)[i];
                    /* New 051005 */
                    if(fabs(pdens(ans)[i]) < 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                        pdens(ans)[i] = 0.0;
                    else if(fabs(pdens(ans)[i]) > 10.0*MACH_EPS && pdens(ans)[i] < 0.0)
                    {
                        printf("ERROR in g_quad_interpolator()\n");
                        printf("partial density < 0.0\n");
                        clean_up(ERROR);
                    }
                    /* End of New 051005 */
                }
                /* TMP */
                sum = 0.0;
                for(i = 0; i < Params(ans)->n_comps; i++)
                    sum += pdens(ans)[i];
                /* New 051005 */
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = Dens(ans)/sum*pdens(ans)[i];
                /* End of New 051005 */

                if(fabs(sum - Dens(ans)) > 0.0000001)
                {
                    printf("ERROR  g_quad_interpolator()\n");
                    printf("partial density not consistent\n");
                    verbose_print_state("ans", ans);
                    verbose_print_state("s0", s0);
                    verbose_print_state("s1", s1);
                    verbose_print_state("s2", s2);
                    verbose_print_state("s3", s3);
                    clean_up(ERROR);
                }
            }    
        }
	Energy(ans) = f0*Energy(s0) + f1*Energy(s1) + f2*Energy(s2) 
							+ f3*Energy(s3);
	for (i = 0; i < dim; i++)
	{
	    Mom(ans)[i] = f0*Mom(s0)[i] + f1*Mom(s1)[i] +
	    		  f2*Mom(s2)[i] + f3*Mom(s3)[i];
	}
#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = f0*Prod(s0) + f1*Prod(s1) +
	    		f2*Prod(s2) + f3*Prod(s3);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = f0*Prod(s0) + f1*Prod(s1) +
	    		f2*Prod(s2) + f3*Prod(s3);
	    Dens1(ans) = f0*Dens1(s0) + f1*Dens1(s1) +
	    		 f2*Dens1(s2) + f3*Dens1(s3);
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
}		/*end g_quad_interpolator*/


EXPORT void gt_quad_interpolator(
	double		*f,
	BILINEAR_ELEMENT *q,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	Locstate	sts[(1<<MAXD)];
	Locstate	s0, s1, s2, s3;
	double		fx = f[0], fy = f[1];
	double		f0,f1,f2,f3;
	double		T0, T1, T2, T3;
	int		dim, i;
	static Locstate temp_ans = NULL;

	states_on_bilinear_element(sts,q,soln->tri_grid);
	s0 = sts[0];	s1 = sts[1];	s2 = sts[2];	s3 = sts[3];

	if (is_obstacle_state(s0) && is_obstacle_state(s1)
		 && is_obstacle_state(s2) && is_obstacle_state(s3)) 
	{
		g_obstacle_state(ans,g_sizest());
		return;
	}
	if (	(	Different_params(s0,s1) ||
			Different_params(s0,s2) ||
			Different_params(s0,s3)	)
	&&
		correct_bad_params_in_quad(q,soln) == NO)
	{
		screen("ERROR in gt_quad_interpolator(), ");
		screen("different params\n");
		(void) printf("params0 = %llu params1 = %llu ",
			gas_param_number(Params(s0)),
			gas_param_number(Params(s1)));
		(void) printf("params2 = %llu params3 = %llu\n",
		        gas_param_number(Params(s2)),
			gas_param_number(Params(s3)));
		(void) printf("fx = %g, fy = %g\n",fx,fy);
		print_BILINEAR_ELEMENT(q,soln->tri_grid);
		verbose_print_state("s0",s0);
		verbose_print_state("s1",s1);
		verbose_print_state("s2",s2);
		verbose_print_state("s3",s3);
		print_interface(current_interface());
		clean_up(ERROR);
	}

	if (temp_ans == NULL)
	{
		(*Params(s0)->_alloc_state)(&temp_ans,Params(s0)->sizest);
		set_type_of_state(temp_ans,FGAS_STATE);
	}

	Set_params(ans,s0);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;
	
	f1 = fx * (1.0 - fy);
	f2 = fx * fy;
	f3 = fy * (1.0 - fx);
	f0 = 1.0 - f1 - f2 - f3;

	if (is_rotational_symmetry())
	    (*quad_weighted_coefs)(&f0,&f1,&f2,&f3,q,&soln->tri_grid->comp_grid);

	Dens(ans) = pow(Dens(s0),f0)*pow(Dens(s1),f1)*
			pow(Dens(s2),f2)*pow(Dens(s3),f3);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            double sum;
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = pow(pdens(s0)[i],f0)*pow(pdens(s1)[i],f1)*
                                    pow(pdens(s2)[i],f2)*pow(pdens(s3)[i],f3);
                sum = 0.0;
                for(i = 0; i < Params(ans)->n_comps; i++)
                    sum += pdens(ans)[i];
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = pdens(ans)[i]/sum*Dens(ans);
            }
        }
	for (i = 0; i < dim; i++)
	{
		Mom(ans)[i] = Dens(ans)*(f0*vel(i,s0) +
					f1*vel(i,s1) +
					f2*vel(i,s2) +
					f3*vel(i,s3));
	}

	Dens(temp_ans) = Dens(ans);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(temp_ans)[i] = pdens(ans)[i];
            }
        }
	T0 = temperature(s0);
	T1 = temperature(s1);
	T2 = temperature(s2);
	T3 = temperature(s3);
	Temperature(temp_ans) = pow(T0,f0)*pow(T1,f1)*pow(T2,f2)*pow(T3,f3);
	Set_params(temp_ans,ans);
	reset_gamma(temp_ans);
	Energy(ans) = internal_energy(temp_ans) + kinetic_energy(ans);
	
#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
		break;
	case ZND:
		Prod(ans) = Dens(ans)*(f0*react(s0) +
					f1*react(s1) +
					f2*react(s2) +
					f3*react(s3));
		break;
	case TWO_CONSTITUENT_REACTIVE:
		Prod(ans) = Dens(ans)*(f0*react(s0) +
					f1*react(s1) +
					f2*react(s2) +
					f3*react(s3));
		Dens1(ans) = pow(Dens1(s0),f0)*pow(Dens1(s1),f1)*
			pow(Dens1(s2),f2)*pow(Dens1(s3),f3);
		break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
}		/*end gt_quad_interpolator*/

LOCAL	boolean correct_bad_params_in_quad(
	BILINEAR_ELEMENT *q,
	TRI_SOLN	*soln)
{
	Locstate sts[(1<<MAXD)];
	HYPER_SURF	*hs;
	HYPER_SURF_ELEMENT *hse;
	double		ans[MAXD], t[MAXD];
	int		dim, i, j;
	int		n;
	int		bad;

#if !defined(COMBUSTION_CODE)
	if (g_composition_type() != PURE_NON_REACTIVE)
	    return YES;
#endif /* !defined(COMBUSTION_CODE) */

	dim = soln->intfc->dim;
	n = 1 << dim;
	states_on_bilinear_element(sts,q,soln->tri_grid);
	for (i = 1; i < n; i++)
	    if (Different_params(sts[0],sts[i]))
		break;
	if (i == n) return YES;
	/* Compare to 0 */
	for (j = 1; j < n; j++)
	{
	    if (j == i)
		continue;
	    if (Different_params(sts[0],sts[j]))
		break;
	}
	if (j == n) /* All other params agree with 0 */
	{
	    bad = i;
	}
	else
	{
	    for (j = 1; j < n; j++)
	    	if (Different_params(sts[i],sts[j]))
		    break;
	    if (j < n) return NO;
	    bad = 0;
	}

	/*	Three out four params agree on this element
	*	Use nearest interface point to substitute for
	*	the state in disagreement
	*/

	if (nearest_interface_point(Coords(q->p[bad]),q->comp,soln->intfc,
				    INCLUDE_BOUNDARIES,NULL,ans,t,&hse,&hs)
				    != YES)
	{
	    screen("ERROR in correct_bad_params_in_quad(), "
		   "nearest_interface_point() failed\n");
	    clean_up(ERROR);
	}
	state_along_hypersurface_element(q->comp,t,hse,hs,sts[bad]);
	for (i = 1; i < n; i++)
	    if (Different_params(sts[0],sts[i]))
		break;
	if (i < n)
	    return NO;
	return YES;
}		/*end correct_bad_params_in_quad*/


/*
*			g_grad_tri_interpolator():
*/

/* ARGSUSED */
EXPORT boolean g_grad_tri_interpolator(
	double		**df,
	LINEAR_ELEMENT	*t,
	TRI_SOLN	*soln,
	Locstate	*grad)
{
		/* TODO: code needed */

	screen("ERROR: grad_solution() not implemented\n");
	clean_up(ERROR);
	return FUNCTION_FAILED;
}		/*end g_grad_tri_interpolator*/



/*
*			g_grad_quad_interpolator():
*/

/* ARGSUSED */
EXPORT void g_grad_quad_interpolator(
	double		*f,
	double		*d,
	BILINEAR_ELEMENT *q,
	TRI_SOLN	*soln,
	Locstate	*grad)
{
		/* TODO: code needed */

	screen("ERROR: grad_solution() not implemented\n");
	clean_up(ERROR);
}		/*end g_grad_quad_interpolator*/
#endif /* defined(TWOD) */

/*
*	
*	There are currently several different interpolation options for
*	geometry.
*
*	   1) No geometry weighting applied to the interpolation
*	      coefficients (UNWEIGHTED_INTERP_COEFS).
*
*	   2) Linear interpolation using volume averaged states
*	      (RGAS_WEIGHTED_INTERP_COEFS),  in this case, r^a times the
*	      interpolated quantities are assumed to be linear,  where
*	      a = 0 for no rotational symmetry, a = 1 for cylindrical
*	      symmetry, and a = 2 for spherical symmetry.
*
*	   3) Linear interpolation using weights based
*	      on the volume separating states (VOLUME_WEIGHTED_INTERP_COEFS).
*
*	   4) Use the combination of 2 and 3.
*	      (RGAS_WEIGHTED_INTERP_COEFS | VOLUME_WEIGHTED_INTERP_COEFS).
*
*/

EXPORT	void	set_interpolation_weights(
	INTERPOLATION_WEIGHT		flag)
{
	double	alpha = rotational_symmetry();

	lin_weighted_coefs = lin_unweighted_coefs;
#if defined(TWOD) || defined(THREED)
	tri_weighted_coefs = tri_unweighted_coefs;
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
	quad_weighted_coefs = quad_unweighted_coefs;
#endif /* defined(TWOD) */
	if (alpha == 0.0)
	    return;
	if (alpha == 1.0)
	{
	    switch (flag)
	    {
	    case UNWEIGHTED_INTERP_COEFS:
	    	return;
	    case RGAS_WEIGHTED_INTERP_COEFS:
	    	lin_weighted_coefs = cyl_lin_rgas_weighted_coefs;
#if defined(TWOD) || defined(THREED)
	    	tri_weighted_coefs = cyl_tri_rgas_weighted_coefs;
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
	    	quad_weighted_coefs = cyl_quad_rgas_weighted_coefs;
#endif /* defined(TWOD) */
	    	return;
	    case RGAS_VOLUME_WEIGHTED_INTERP_COEFS:
	    	lin_weighted_coefs = cyl_lin_vol_rgas_weighted_coefs;
#if defined(TWOD) || defined(THREED)
	    	tri_weighted_coefs = cyl_tri_vol_rgas_weighted_coefs;
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
	    	quad_weighted_coefs = cyl_quad_vol_rgas_weighted_coefs;
#endif /* defined(TWOD) */
	    	return;
	    default:
	    case VOLUME_WEIGHTED_INTERP_COEFS:
	    	lin_weighted_coefs = cyl_lin_vol_weighted_coefs;
#if defined(TWOD) || defined(THREED)
	    	tri_weighted_coefs = cyl_tri_vol_weighted_coefs;
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
	    	quad_weighted_coefs = cyl_quad_vol_weighted_coefs;
#endif /* defined(TWOD) */
	    	return;
	    }
	}
	if (alpha == 2.0)
	{
	    switch (flag)
	    {
	    case UNWEIGHTED_INTERP_COEFS:
	    	return;
	    case RGAS_WEIGHTED_INTERP_COEFS:
	    	lin_weighted_coefs = sph_lin_rgas_weighted_coefs;
	    	return;
	    case RGAS_VOLUME_WEIGHTED_INTERP_COEFS:
	    	lin_weighted_coefs = sph_lin_vol_rgas_weighted_coefs;
	    	return;
	    default:
	    case VOLUME_WEIGHTED_INTERP_COEFS:
	    	lin_weighted_coefs = sph_lin_vol_weighted_coefs;
	    	return;
	    }
	}
}		/*end set_interpolation_weights*/



#define	lin_pos_radius(r1,r2,gr)					\
	if ((gr)->Remap.remap != IDENTITY_REMAP)			\
	{								\
	    double rmin = pos_radius(0.0,gr);				\
	    if ((r1 > 0.0 && r2 < 0.0) ||				\
	        (r1 < 0.0 && r2 > 0.0))					\
		return;							\
	    (r1) = pos_radius(r1,gr);					\
	    (r2) = pos_radius(r2,gr);					\
	    if ((fabs((r1)) <= rmin) || (fabs((r2)) <= rmin))		\
		return;							\
	}

#define	check_lin_coefs_for_modification(a,b,r1,r2,gr)			\
	lin_pos_radius(r1,r2,gr)					\
	if (fabs((a) + (b) - 1.0) > EPS)				\
	{ 								\
	    (void) printf("WARNING, Sum of coefficients in "		\
	                  "linear interpolator = %g\n",(a) + (b));	\
	    (void) printf("\tASSUME function is not used "		\
	                  "for interpolation, "				\
	                  "Coefficients will not be modified\n");	\
	    return;							\
	}								\

/*ARGSUSED*/
LOCAL	void	lin_unweighted_coefs(
	double		*alpha,
	double		*beta,
	double		r1,
	double		r2,
	RECT_GRID	*gr)
{
}		/*end lin_unweighted_coefs*/


/*ARGSUSED*/
LOCAL	void	cyl_lin_vol_weighted_coefs(
	double		*alpha,
	double		*beta,
	double		r1,
	double		r2,
	RECT_GRID	*gr)
{
	double		ra;

	check_lin_coefs_for_modification(*alpha,*beta,r1,r2,gr);
	ra = *alpha * r1 + *beta * r2;
	*alpha *= (r2 + ra) / (r1 + r2);
	*beta = 1.0 - *alpha;
}		/*end cyl_lin_vol_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	cyl_lin_rgas_weighted_coefs(
	double		*alpha,
	double		*beta,
	double		r1,
	double		r2,
	RECT_GRID	*gr)
{
	double		ra;

	check_lin_coefs_for_modification(*alpha,*beta,r1,r2,gr);
	ra = *alpha * r1 + *beta * r2;
	*alpha *= r1 / ra;
	*beta = 1.0 - *alpha;
}		/*end cyl_lin_rgas_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	cyl_lin_vol_rgas_weighted_coefs(
	double		*alpha,
	double		*beta,
	double		r1,
	double		r2,
	RECT_GRID	*gr)
{
	double		ra;

	check_lin_coefs_for_modification(*alpha,*beta,r1,r2,gr);
	ra = *alpha * r1 + *beta * r2;
	*alpha *= (r1* (r2 + ra)) / (ra * (r1 + r2));
	*beta = 1.0 - *alpha;
}		/*end cyl_lin_vol_rgas_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	sph_lin_vol_weighted_coefs(
	double		*alpha,
	double		*beta,
	double		r1,
	double		r2,
	RECT_GRID	*gr)
{
	double		ra;

	check_lin_coefs_for_modification(*alpha,*beta,r1,r2,gr);
	ra = *alpha * r1 + *beta * r2;
	*alpha *= (sqr(r2) + r2*ra + sqr(ra)) / (sqr(r1) + r1*r2 + sqr(r2));
	*beta = 1.0 - *alpha;
}		/*end sph_lin_vol_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	sph_lin_rgas_weighted_coefs(
	double		*alpha,
	double		*beta,
	double		r1,
	double		r2,
	RECT_GRID	*gr)
{
	double		ra;

	check_lin_coefs_for_modification(*alpha,*beta,r1,r2,gr);
	r1 = sqr(r1);	r2 = sqr(r2);
	ra = *alpha * r1 + *beta * r2;
	*alpha *= r1 / ra;
	*beta = 1.0 - *alpha;
}		/*end sph_lin_rgas_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	sph_lin_vol_rgas_weighted_coefs(
	double		*alpha,
	double		*beta,
	double		r1,
	double		r2,
	RECT_GRID	*gr)
{
	double		ra;

	check_lin_coefs_for_modification(*alpha,*beta,r1,r2,gr);
	ra = *alpha * r1 + *beta * r2;
	*alpha *= (sqr(r2) + r2*ra + sqr(ra)) / (sqr(r1) + r1*r2 + sqr(r2));
	r1 = sqr(r1);	r2 = sqr(r2);
	ra = *alpha * r1 + *beta * r2;
	*alpha *= r1 / ra;
	*beta = 1.0 - *alpha;
}		/*end sph_lin_vol_rgas_weighted_coefs*/

#if defined(TWOD) || defined(THREED)

#define	tri_pos_radius(r0,r1,r2,gr)					\
	if ((gr)->Remap.remap != IDENTITY_REMAP)			\
	{								\
	    double rmin = pos_radius(0.0,gr);				\
	    if (r0*r1 < 0.0 || r1*r2 < 0.0 || r0*r2 < 0.0)		\
	        return;							\
	    (r0) = pos_radius(r0,gr);					\
	    (r1) = pos_radius(r1,gr);					\
	    (r2) = pos_radius(r2,gr);					\
	    if ((fabs((r0)) <= rmin) || (fabs((r1)) <= rmin) || 	\
	    	(fabs((r2)) <= rmin))					\
		return;							\
	}

#define	check_tri_coefs_for_modification(f0,f1,f2,r0,r1,r2,gr)		\
	tri_pos_radius(r0,r1,r2,gr)					\
	if (fabs((f0) + (f1) + (f2) - 1.0) > EPS)			\
	{ 								\
	    (void) printf("WARNING Sum of coefficients "		\
	                  "in tri interpolator = %g\n",(f0)+(f1)+(f2));	\
	    (void) printf("\tASSUME function is not used "		\
	                  "for interpolation, "				\
	                  "Coefficients will not be modified\n");	\
	    return;							\
	}

/*ARGSUSED*/
LOCAL	void	tri_unweighted_coefs(
	double		*f0,
	double		*f1,
	double		*f2,
	double		*crds0,
	double		*crds1,
	double		*crds2,
	RECT_GRID	*gr)
{
}		/*end tri_unweighted_coefs*/


/*ARGSUSED*/
LOCAL	void	cyl_tri_vol_weighted_coefs(
	double		*f0,
	double		*f1,
	double		*f2,
	double		*crds0,
	double		*crds1,
	double		*crds2,
	RECT_GRID	*gr)
{
	double		r, r0,r1,r2;
	double		rsum;

	r0 = crds0[0];	r1 = crds1[0];	r2 = crds2[0];
	check_tri_coefs_for_modification(*f0,*f1,*f2,r0,r1,r2,gr);

	rsum = r0 + r1 + r2;
	r = *f0*r0 + *f1*r1 + *f2*r2;
	*f1 *= (r0 + r2 + r) / rsum;
	*f2 *= (r0 + r1 + r) / rsum;
	*f0 = 1.0 - *f1 - *f2;
}		/*end cyl_tri_vol_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	cyl_tri_rgas_weighted_coefs(
	double		*f0,
	double		*f1,
	double		*f2,
	double		*crds0,
	double		*crds1,
	double		*crds2,
	RECT_GRID	*gr)
{
	double		r, r0,r1,r2;

	r0 = crds0[0];	r1 = crds1[0];	r2 = crds2[0];
	check_tri_coefs_for_modification(*f0,*f1,*f2,r0,r1,r2,gr);

	r = *f0*r0 + *f1*r1 + *f2*r2;
	*f1 *= r1 / r;
	*f2 *= r2 / r;
	*f0 = 1.0 - *f1 - *f2;
}		/*end cyl_tri_rgas_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	cyl_tri_vol_rgas_weighted_coefs(
	double		*f0,
	double		*f1,
	double		*f2,
	double		*crds0,
	double		*crds1,
	double		*crds2,
	RECT_GRID	*gr)
{
	double		r, r0,r1,r2;
	double		rsum;

	r0 = crds0[0];
	r1 = crds1[0];
	r2 = crds2[0];
	check_tri_coefs_for_modification(*f0,*f1,*f2,r0,r1,r2,gr);

	r = *f0*r0 + *f1*r1 + *f2*r2;
	rsum = (r0 + r1 + r2)*r;
	*f1 *= r1 * (r0 + r2 + r) / rsum;
	*f2 *= r2 * (r0 + r1 + r) / rsum;
	*f0 = 1.0 - *f1 - *f2;
}		/*end cyl_tri_vol_rgas_weighted_coefs*/
#endif /* defined(TWOD) || defined(THREED) */

#if defined(TWOD)

#define	quad_pos_radius(r0,r1,r2,r3,gr)					\
	if ((gr)->Remap.remap != IDENTITY_REMAP)			\
	{								\
	    double rm = pos_radius(0.0,gr);				\
	    double rmax = max(max(r0,r1),max(r2,r3));			\
	    double rmin = min(min(r0,r1),min(r2,r3));			\
	    if (rmax*rmin < 0.0)					\
	    	return;							\
	    (r0) = pos_radius(r0,gr);					\
	    (r1) = pos_radius(r1,gr);					\
	    (r2) = pos_radius(r2,gr);					\
	    (r3) = pos_radius(r3,gr);					\
	    if ((fabs((r0))<=rm) || (fabs((r1))<=rm) || 		\
	    	(fabs((r2))<=rm) || (fabs((r3))<=rm))			\
		return;							\
	}

#define	check_quad_coefs_for_modification(f0,f1,f2,f3,r0,r1,r2,r3,gr)	\
	quad_pos_radius(r0,r1,r2,r3,gr)					\
	if (fabs((f0) + (f1) + (f2) +(f3) - 1.0) > EPS)			\
	{ 								\
	    (void) printf("WARNING: Sum of coefficients in quad "	\
	                  "interpolator = %g\n",(f0)+(f1)+(f2)+(f3));	\
	    (void) printf("\tASSUME function is not used "		\
	                  "for interpolation, "				\
	                  "Coefficients will not be modified\n");	\
	    return;							\
	}

/*ARGSUSED*/
LOCAL	void	quad_unweighted_coefs(
	double		 *f0,
	double		 *f1,
	double		 *f2,
	double		 *f3,
	BILINEAR_ELEMENT *q,
	RECT_GRID        *gr)
{
}		/*end quad_unweighted_coefs*/


/*ARGSUSED*/
LOCAL	void	cyl_quad_vol_weighted_coefs(
	double		 *f0,
	double		 *f1,
	double		 *f2,
	double		 *f3,
	BILINEAR_ELEMENT *q,
	RECT_GRID	 *gr)
{
	double		rsum;
	double		r, r0,r1,r2,r3;

	r0 = Coords(q->p[0])[0];
	r1 = Coords(q->p[1])[0];
	r2 = Coords(q->p[2])[0];
	r3 = Coords(q->p[3])[0];
	check_quad_coefs_for_modification(*f0,*f1,*f2,*f3,r0,r1,r2,r3,gr);

	r = *f0*r0 + *f1*r1 + *f2*r2 + *f3*r3;

	/* TODO: This may only be correct for rectangles */
	rsum = r0 + r1 + r2 + r3;
	*f0 *= (r1 + r2 + r3 + r) / rsum;
	*f1 *= (r0 + r2 + r3 + r) / rsum;
	*f2 *= (r0 + r1 + r3 + r) / rsum;
	*f3 = 1 - *f0 - *f1 - *f2;
}		/*end cyl_quad_vol_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	cyl_quad_rgas_weighted_coefs(
	double		 *f0,
	double		 *f1,
	double		 *f2,
	double		 *f3,
	BILINEAR_ELEMENT *q,
	RECT_GRID        *gr)
{
	double		r, r0,r1,r2,r3;

	r0 = Coords(q->p[0])[0];
	r1 = Coords(q->p[1])[0];
	r2 = Coords(q->p[2])[0];
	r3 = Coords(q->p[3])[0];
	check_quad_coefs_for_modification(*f0,*f1,*f2,*f3,r0,r1,r2,r3,gr);

	r = *f0*r0 + *f1*r1 + *f2*r2 + *f3*r3;
	*f0 *= r0/r;
	*f1 *= r1/r;
	*f2 *= r2/r;
	*f3 = 1 - *f0 - *f1 - *f2;
}		/*end cyl_quad_rgas_weighted_coefs*/

/*ARGSUSED*/
LOCAL	void	cyl_quad_vol_rgas_weighted_coefs(
	double		 *f0,
	double		 *f1,
	double		 *f2,
	double		 *f3,
	BILINEAR_ELEMENT *q,
	RECT_GRID	 *gr)
{
	double		rsum;
	double		r, r0,r1,r2,r3;

	r0 = Coords(q->p[0])[0];
	r1 = Coords(q->p[1])[0];
	r2 = Coords(q->p[2])[0];
	r3 = Coords(q->p[3])[0];
	check_quad_coefs_for_modification(*f0,*f1,*f2,*f3,r0,r1,r2,r3,gr);

	r = *f0*r0 + *f1*r1 + *f2*r2 + *f3*r3;

	/* TODO: This may only be correct for rectangles */
	rsum = (r0 + r1 + r2 + r3)*r;
	*f0 *= r0*(r1 + r2 + r3 + r) / rsum;
	*f1 *= r1*(r0 + r2 + r3 + r) / rsum;
	*f2 *= r2*(r0 + r1 + r3 + r) / rsum;
	*f3 = 1 - *f0 - *f1 - *f2;
}		/*end cyl_quad_vol_rgas_weighted_coefs*/

#endif /* defined(TWOD) */

/*
*
*			pos_radius():
*
*	This checks whether the radius is positive.  If it is, 
*	the function returns r; otherwise returns rmin = dr^3.
*
*	Nothing is done in cartesian geometry.
*
*	Input:
*	       	r               radius being checked
*	        gr              rect_grid used to set rmin
*	Output: 		corrected radius
*		If r <= rmin a message is written on the screen.
*/

EXPORT double pos_radius(
	double		r,
	RECT_GRID	*gr)
{
	double rmin = 0.25*gr->h[0];

	switch (gr->Remap.remap)
	{
	case IDENTITY_REMAP:
	    return r;
	case CYLINDRICAL_REMAP:
	case SPHERICAL_REMAP:
	    r = (r - gr->GL[0] > 0.0) ? r : (2.0*gr->GL[0] - r);
	    return (r > rmin) ? r : rmin;
	}
}		/*end pos_radius*/


#if defined(THREED)
/*
*			g_tetra_interpolator():
*
*	Forms a linear combination of three local states.  This function
*	is generally used to calculate an interpolated state from the four
*	states at the vertices of a tetrahedral.
*
*/

/*ARGSUSED*/
EXPORT boolean g_tetra_interpolator(
	double		*f,
	LINEAR_ELEMENT	*t,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	Locstate	*s = t->s;
					/* states in linear combination */
	int		dim;
	int		i,j;

	for (i = 1; i < 4; i++)
	{
	    if (Different_params(s[i],s[0]))
	    {
	    	int	j;
	    	screen("ERROR in g_tetra_interpolator(), ");
	    	screen("different params\n");
	    	for (j = 0; j < 4; j++)
	    	{
	    	    (void) printf("Params(s[%d]) = %llu\n",j,
		    		  gas_param_number(Params(s[j])));
	    	    (void) printf("t->p[%d] = ",j);
	    	    print_general_vector("",Coords(t->p[j]),3,"\n");
	    	}
	    	clean_up(ERROR);
	    }
	}

	if (is_obstacle_state(s[0]))
	{
	    g_obstacle_state(ans,g_sizest());
	    return FUNCTION_SUCCEEDED;
	}
	Set_params(ans,s[0]);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	Dens(ans) = 0.0;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = 0.0;     
            }
        }
	Energy(ans) = 0.0;
	for (j = 0; j < dim; j++)
	    Mom(ans)[j] = 0.0;
	for (i = 0; i < 4; i++)
	{
	    Dens(ans) += f[i]*Dens(s[i]);
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
	       int k;
               if(Params(ans)->n_comps != 1)
               {
                   for(k = 0; k < Params(ans)->n_comps; k++)
                   {
                       pdens(ans)[k] += f[i]*pdens(s[i])[k];
                    }
               }
            } 
	    Energy(ans) += f[i]*Energy(s[i]);
	    for (j = 0; j < dim; j++) 
	    	Mom(ans)[j] += f[i]*Mom(s[i])[j];
	}

#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = 0.0;
	    for (i = 0; i < 4; i++)
	    	Prod(ans) += f[i]*Prod(s[i]);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = 0.0;
	    Dens1(ans) = 0.0;
	    for (i = 0; i < 4; i++)
	    {
	    	Prod(ans) += f[i]*Prod(s[i]);
	    	Dens1(ans) += f[i]*Dens1(s[i]);
	    }
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
	return FUNCTION_SUCCEEDED;
}		/*end g_tetra_interpolator*/

/*
*			gt_tetra_interpolator():
*
*	Forms a linear combination of three local states.  This function
*	is generally used to calculate an interpolated state from the four
*	states at the vertices of a tetrahedral.
*
*/

/*ARGSUSED*/
EXPORT boolean gt_tetra_interpolator(
	double		*f,
	LINEAR_ELEMENT	*t,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	double		f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3];
	double		T0, T1, T2, T3;
	Locstate	s0 = t->s[0], s1 = t->s[1], s2 = t->s[2], s3 = t->s[3];
	Locstate	*s = t->s;
	static Locstate temp_ans = NULL;
					/* states in linear combination */
	int		dim;
	int		i;

	for (i = 1; i < 4; i++)
	{
	    if (Different_params(s[i],s[0]))
	    {
	    	screen("ERROR in gt_tetra_interpolator(), "
	    	       "different params\n");
	    	clean_up(ERROR);
	    }
	}

	if (is_obstacle_state(s[0]))
	{
	    g_obstacle_state(ans,g_sizest());
	    return FUNCTION_SUCCEEDED;
	}
	if (temp_ans == NULL)
	{
	    (*Params(s0)->_alloc_state)(&temp_ans,Params(s0)->sizest);
	    set_type_of_state(temp_ans,FGAS_STATE);
	}
	Set_params(ans,s0);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	Dens(ans) = pow(Dens(s0),f0)*pow(Dens(s1),f1)*
			pow(Dens(s2),f2)*pow(Dens(s3),f3);
	for (i = 0; i < dim; i++)
	{
	    Mom(ans)[i] = Dens(ans)*(f0 * vel(i,s0) +
				f1 * vel(i,s1) +
				f2 * vel(i,s2) +
				f3 * vel(i,s3));
	}

	Dens(temp_ans) = Dens(ans);
	T0 = temperature(s0);
	T1 = temperature(s1);
	T2 = temperature(s2);
	T3 = temperature(s3);
	Temperature(temp_ans) = pow(T0,f0)*pow(T1,f1)*pow(T2,f2)*pow(T3,f3);
	Set_params(temp_ans,ans);
	reset_gamma(temp_ans);
	Energy(ans) = internal_energy(temp_ans) + kinetic_energy(ans);

#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = Dens(ans)*(f0 * react(s0) + f1 * react(s1) +
				   f2 * react(s2) + f3 * react(s3));
		break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = Dens(ans)*(f0 * react(s0) + f1 * react(s1) +
	    		           f2 * react(s2) + f3 * react(s3));
	    Dens1(ans) = pow(Dens1(s0),f0)*pow(Dens1(s1),f1)*pow(Dens1(s2),f2)*
			 pow(Dens1(s3),f3);
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
	return FUNCTION_SUCCEEDED;
}		/*end gt_tetra_interpolator*/

/*
*			g_cube_interpolator():
*
*	Forms a bilinear combination of eight local states.  This function
*	is generally used to calculate an interpolated state from the eight 
*	states at the vertices of a cube.
*
*/


EXPORT void g_cube_interpolator(
	double		*ff,
	BILINEAR_ELEMENT *cb,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	double		Lx[2],Ly[2],Lz[2];
	double		f[8];
	int		dim, i,j,k;
	int		np = 0;
	Locstate	s[8];

	states_on_bilinear_element(s,cb,soln->tri_grid);

	for (i = 1; i < 8; i++)
	{
	    if (Different_params(s[i],s[0]))
	    {
	    	screen("ERROR in g_cube_interpolator(), "
	    	       "different params\n");
	    	print_BILINEAR_ELEMENT(cb,soln->tri_grid);
	    	for (j = 0; j < 8; j++)
	    	{
	    	    char str[5];

	    	    (void) sprintf(str,"[%d]",j);
	    	    verbose_print_state(str,s[j]);
	    	}
	    	print_interface(current_interface());
	    	clean_up(ERROR);
	    }
	}

	if (is_obstacle_state(s[0]))
	{
	    g_obstacle_state(ans,g_sizest());
	    return;
	}

	Set_params(ans,s[0]);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	Lx[0] = 1.0 - ff[0]; 		Lx[1] = ff[0];	
	Ly[0] = 1.0 - ff[1];		Ly[1] = ff[1];	
	Lz[0] = 1.0 - ff[2];		Lz[1] = ff[2];
	for (i = 0; i < 2; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
	    	for (k = 0; k < 2; k++)
	    	{
	    	    f[np++] = Lx[i]*Ly[j]*Lz[k];
	    	}
	    }
	}
	
	Dens(ans) = 0.0;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ans)->n_comps != 1)
            {
                for(i = 0; i < Params(ans)->n_comps; i++)
                    pdens(ans)[i] = 0.0;
            }
        }
	Energy(ans) = 0.0;
	for (j = 0; j < dim; j++) Mom(ans)[j] = 0.0;
	for (i = 0; i < 8; i++)
	{
	    Dens(ans) += f[i]*Dens(s[i]);
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
               if(Params(ans)->n_comps != 1)
               {
                   for(k = 0; k < Params(ans)->n_comps; k++)
                   {
                       pdens(ans)[k] += f[i]*pdens(s[i])[k];
                    }
               }
            }
	    Energy(ans) += f[i]*Energy(s[i]);
	    for (j = 0; j < dim; j++) 
	    	Mom(ans)[j] += f[i]*Mom(s[i])[j];
	}
	
#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    Prod(ans) = 0.0;
	    for (i = 0; i < 8; i++)
	    	Prod(ans) += f[i]*Prod(s[i]);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    Prod(ans) = 0.0;
	    Dens1(ans) = 0.0;
	    for (i = 0; i < 8; i++)
	    {
	    	Prod(ans) += f[i]*Prod(s[i]);
	    	Dens1(ans) += f[i]*Dens1(s[i]);
	    }
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
}		/*end g_cube_interpolator*/

/*
*			gt_cube_interpolator():
*
*	Forms a bilinear combination of eight local states.  This function
*	is generally used to calculate an interpolated state from the eight 
*	states at the vertices of a cube.
*
*/


EXPORT void gt_cube_interpolator(
	double		*ff,
	BILINEAR_ELEMENT *cb,
	TRI_SOLN	*soln,
	Locstate	ans)
{
	double		Lx[2],Ly[2],Lz[2];
	double		f[8];
	double		T;
	int		dim, i,j,k;
	int		np = 0;
	Locstate	s[8];
	static Locstate temp_ans = NULL;

	states_on_bilinear_element(s,cb,soln->tri_grid);

	for (i = 1; i < 7; i++)
	{
	    if (Different_params(s[i],s[0]))
	    {
	    	screen("ERROR in gt_cube_interpolator(), "
	    	       "different params\n");
	    	clean_up(ERROR);
	    }
	}


	if (is_obstacle_state(s[0]))
	{
	    g_obstacle_state(ans,g_sizest());
	    return;
	}
	if (temp_ans == NULL)
	{
	    (*Params(s[0])->_alloc_state)(&temp_ans,Params(s[0])->sizest);
	    set_type_of_state(temp_ans,FGAS_STATE);
	}

	Set_params(ans,s[0]);	/* TODO: FIX THIS */
	set_type_of_state(ans,GAS_STATE);
	dim = Params(ans)->dim;

	Lx[0] = 1.0 - ff[0]; 		Lx[1] = ff[0];	
	Ly[0] = 1.0 - ff[1];		Ly[1] = ff[1];	
	Lz[0] = 1.0 - ff[2];		Lz[1] = ff[2];
	for (i = 0; i < 2; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
	    	for (k = 0; k < 2; k++)
	    	{
	    	    f[np++] = Lx[i]*Ly[j]*Lz[k];
		}
	    }
	}
	
	Dens(ans) = 1.0;
	T = 1.0;
	for (j = 0; j < dim; j++) Vel(ans)[j] = 0.0;
	for (i = 0; i < 8; i++)
	{
	    Dens(ans) *= pow(Dens(s[i]),f[i]);
	    T *= pow(temperature(s[i]),f[i]);
	    for (j = 0; j < dim; j++) 
	    	Vel(ans)[j] += f[i]*Vel(s[i])[j];
	}
	for (j = 0; j < dim; j++) 
	    Mom(ans)[j] = Dens(ans)*Vel(ans)[j];
	Dens(temp_ans) = Dens(ans);
	Temperature(temp_ans) = T;
	Set_params(temp_ans,ans);
	reset_gamma(temp_ans);
	Energy(ans) = internal_energy(temp_ans) + kinetic_energy(ans);
	
#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	case PURE_NON_REACTIVE:
	case THINFLAME:
	case PTFLAME:
	    break;
	case ZND:
	    React(ans) = 0.0;
	    for (i = 0; i < 8; i++)
	    	React(ans) += f[i]*react(s[i]);
	    Prod(ans) = Dens(ans)*React(ans);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    React(ans) = 0.0;
	    Dens1(ans) = 1.0;
	    for (i = 0; i < 8; i++)
	    {
	    	React(ans) += f[i]*react(s[i]);
	    	Dens1(ans) *= pow(Dens1(s[i]),f[i]);
	    }
	    Prod(ans) = Dens(ans)*React(ans);
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);
}		/*end gt_cube_interpolator*/

#endif /* defined(THREED) */

EXPORT	boolean g_least_sqr_interpolator(
	LEAST_SQR_CLUSTER *cluster,
	TRI_SOLN        *soln,
	double		*crds,
	Locstate        ans)
{
	int dim = cluster->dim;
	int p = cluster->nc;
	int n = cluster->nr;
    	int i, j;
	Locstate *s = cluster->s;
	double *x = cluster->x;
    	static double *y, *b;
	static int cur_n,cur_p;

	if (n < p) 
	{
	    screen("ERROR: least square has too few points (%d)!\n",n);
	    clean_up(ERROR);
	}
	if (y == NULL)
	{
	    stat_vector(&y,MAX_LSQ_PTS,FLOAT);
	    stat_vector(&b,p,FLOAT);
	    cur_n = n;
	    cur_p = p;
	}
	else if (n > cur_n || p > cur_p)
	{
	    free_these(2,y,b);
	    stat_vector(&y,MAX_LSQ_PTS,FLOAT);
	    stat_vector(&b,p,FLOAT);
	    cur_n = n;
	    cur_p = p;
	}

	Set_params(ans,s[0]);
	set_type_of_state(ans,GAS_STATE);
        for (i = 0; i < n; ++i)
	    y[i] = Dens(s[i]);
        least_sqr_fit(x,y,b,n,p);
	Dens(ans) = b[0]*crds[0] + b[1]*crds[1] + b[2]*crds[2] + b[3];

        for (i = 0; i < n; ++i)
	    y[i] = Energy(s[i]);
        least_sqr_fit(x,y,b,n,p);
	Energy(ans) = b[0]*crds[0] + b[1]*crds[1] + b[2]*crds[2] + b[3];
	for (i = 0; i < dim; ++i)
	{
	    for (j = 0; j < n; ++j)
	    	y[j] = Mom(s[j])[i];
            least_sqr_fit(x,y,b,n,p);
	    for (j = 0; j < n; ++j)
	    	Mom(ans)[i] = b[0]*crds[0] + b[1]*crds[1] + 
				b[2]*crds[2] + b[3];
	}
#if defined(COMBUSTION_CODE)
	switch (Composition_type(ans))
	{
	    case PURE_NON_REACTIVE:
	    case THINFLAME:
            case PTFLAME:
	    	break;
	    case ZND:
	    	for (i = 0; i < n; ++i)
		    y[i] = Prod(s[i]);
                least_sqr_fit(x,y,b,n,p);
		Prod(ans) = b[0]*crds[0] + b[1]*crds[1] + b[2]*crds[2] + b[3];
	    case TWO_CONSTITUENT_REACTIVE:
	    	for (i = 0; i < n; ++i)
		    y[i] = Prod(s[i]);
                least_sqr_fit(x,y,b,n,p);
		Prod(ans) = b[0]*crds[0] + b[1]*crds[1] + b[2]*crds[2] + b[3];
	    	for (i = 0; i < n; ++i)
		    y[i] = Dens1(s[i]);
                least_sqr_fit(x,y,b,n,p);
		Dens1(ans) = b[0]*crds[0] + b[1]*crds[1] + b[2]*crds[2] + b[3];
	}
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(ans);

        return FUNCTION_SUCCEEDED;
}	/* end g_least_sqr_interpolator */


/*******************linpak routine*******************/

LOCAL	void least_sqr_fit_old(
	double *x, 
	double *y, 
	double *b, 
	int n, 
	int p)
{
	int i, ldx, job, info;
	static int *jpvt;
        static double *work,*qraux,*qy,*qty,*rsd,*xb,*bt,*x_copy;
	static int cur_n,cur_p;

	if (jpvt == NULL)
	{
	    stat_vector(&work,p,FLOAT);
	    stat_vector(&bt,p,FLOAT);
	    stat_vector(&jpvt,p,INT);

	    stat_vector(&qraux,n,FLOAT);
	    stat_vector(&qy,n,FLOAT);
	    stat_vector(&qty,n,FLOAT);
	    stat_vector(&rsd,n,FLOAT);
	    stat_vector(&xb,n,FLOAT);
	    stat_vector(&x_copy,n*p,FLOAT);
	    cur_n = n;
	    cur_p = p;
	}
	else if (n > cur_n || p > cur_p)
	{
	    free_these(9,work,bt,jpvt,qraux,qy,qty,rsd,xb,x_copy);

	    stat_vector(&work,p,FLOAT);
	    stat_vector(&bt,p,FLOAT);
	    stat_vector(&jpvt,p,INT);

	    stat_vector(&qraux,n,FLOAT);
	    stat_vector(&qy,n,FLOAT);
	    stat_vector(&qty,n,FLOAT);
	    stat_vector(&rsd,n,FLOAT);
	    stat_vector(&xb,n,FLOAT);
	    stat_vector(&x_copy,n*p,FLOAT);
	    cur_n = n;
	    cur_p = p;
	}

    	job = 1;    /*qr desposition with pivoting*/
    	ldx = n;    /*leading dimension*/

	for (i = 0 ; i < n*p ; i++) x_copy[i] = x[i];

        FORTRAN_NAME(dqrdc)(x_copy,&ldx,&n,&p,qraux,jpvt,work,&job);    

        job = 100;    

   	FORTRAN_NAME(dqrsl)(x_copy,&ldx,&n,&p,qraux,y,qy,qty,bt,rsd,
			xb,&job,&info);    
           
        /*restore the original sequence*/
        for (i = 0 ; i < p ; i++)
        {
	    b[jpvt[i]-1]=bt[i];
	}
        return;
}	/* end least_sqr_fit_old */

/*******************linpak routine******************/

LOCAL	void least_sqr_fit(
	double *x, 
	double *y, 
	double *b, 
	int n, 
	int p)
{
	int i, ldx, job, info;
	static int *jpvt;
        static double *work,*qraux,*qy,*qty,*rsd,*xb,*bt,*x_copy;
	static int cur_n,cur_p;

	if (jpvt == NULL)
	{
	    stat_vector(&work,p,FLOAT);
	    stat_vector(&bt,p,FLOAT);
	    stat_vector(&jpvt,p,INT);

	    stat_vector(&qraux,n,FLOAT);
	    stat_vector(&qy,n,FLOAT);
	    stat_vector(&qty,n,FLOAT);
	    stat_vector(&rsd,n,FLOAT);
	    stat_vector(&xb,n,FLOAT);
	    stat_vector(&x_copy,n*p,FLOAT);
	    cur_n = n;
	    cur_p = p;
	}
	else if (n > cur_n || p > cur_p)
	{
	    free_these(9,work,bt,jpvt,qraux,qy,qty,rsd,xb,x_copy);

	    stat_vector(&work,p,FLOAT);
	    stat_vector(&bt,p,FLOAT);
	    stat_vector(&jpvt,p,INT);

	    stat_vector(&qraux,n,FLOAT);
	    stat_vector(&qy,n,FLOAT);
	    stat_vector(&qty,n,FLOAT);
	    stat_vector(&rsd,n,FLOAT);
	    stat_vector(&xb,n,FLOAT);
	    stat_vector(&x_copy,n*p,FLOAT);
	    cur_n = n;
	    cur_p = p;
	}

    	job = 1;    /*qr desposition with pivoting*/
    	ldx = n;    /*leading dimension*/

	for (i = 0 ; i < n*p ; i++) x_copy[i] = x[i];

        FORTRAN_NAME(dqrdc)(x_copy,&ldx,&n,&p,qraux,jpvt,work,&job);    

        job = 100;    

   	FORTRAN_NAME(dqrsl)(x_copy,&ldx,&n,&p,qraux,y,qy,qty,bt,rsd,
			xb,&job,&info);    
           
        /*restore the original sequence*/
        for (i = 0 ; i < p ; i++)
        {
	    b[jpvt[i]-1]=bt[i];
	}
        return;
}	/* end least_sqr_fit */


/*Apply Laplacian smooth to a point. */
/*only apply it to interior points, make sure p is a point on tri. */
/*ref:  install_curve_points_state */
boolean    g_compute_point_smooth(
        SMOOTH_PARA     *smooth_para,
        SMOOTH_TOL      *stol,
        INTERFACE       *intfc)
{
	static Locstate	stl=NULL, str=NULL;
	Locstate	sl, sr, ansl, ansr;
	POINT		*p, *pt;
	TRI		*tri, **ptris;
	SURFACE		*surf;
	double		alpha, avep[3];
	int		i, j, nt;
	size_t		sizest;
	
	sizest = size_of_state(intfc);
	if(stl == NULL && sizest != 0)
	{
	    alloc_state(intfc, &stl, sizest);
	    alloc_state(intfc, &str, sizest);
	}
	tri = smooth_para->tri;
	p = smooth_para->pt;
	alpha = stol->alpha;

        for(i=0; i<3; i++)
            avep[i] = alpha*smooth_para->avep[i] + (1.0-alpha)*Coords(p)[i];

	ft_assign(Coords(p), avep, 3*FLOAT);
	nt = set_tri_list_around_point(p, tri, &ptris, intfc);
	for(i=0; i<nt; i++)
	    set_normal_of_tri(ptris[i]);

	if(Boundary_point(p))
	    return YES;
	
	if(sizest == 0)
	    return YES;

	if(g_composition_type() != PURE_NON_REACTIVE)
	{
	    printf("ERROR g_compute_point_smooth, average is not implemented");
	    clean_up(ERROR);
	}

	surf = tri->surf;
	slsr(p, Hyper_surf_element(tri), Hyper_surf(surf), &ansl, &ansr);
	
	/*for obstacle states, Dens, Mom, Energy are still valid. */
	for(i=0; i<nt; i++)
	{
	    /*on a surface, all states are GAS_STATE. */
	    /*set_type_of_state(ans,GAS_STATE); */
	    
	    j = Vertex_of_point(ptris[i], p);
	    pt = Point_of_tri(ptris[i])[Prev_m3(j)];
	    slsr(pt, Hyper_surf_element(ptris[i]), Hyper_surf(surf), 
		 &sl, &sr);
	    
	    if(i == 0)
	    {
		ft_assign(stl, sl, sizest);
		ft_assign(str, sr, sizest);
		continue;
	    }

	    Dens(stl) += Dens(sl);
	    Dens(str) += Dens(sr);
	    for(j=0; j<3; j++)
	    {
		Mom(stl)[j] += Mom(sl)[j];
		Mom(str)[j] += Mom(sr)[j];
	    }
	    Energy(stl) += Energy(sl);
	    Energy(str) += Energy(sr);
	}

	Dens(ansl) = alpha*Dens(stl)/nt + (1.0-alpha)*Dens(ansl);
	Dens(ansr) = alpha*Dens(str)/nt + (1.0-alpha)*Dens(ansr);
	for(j=0; j<3; j++)
	{
	    Mom(ansl)[j] = alpha*Mom(stl)[j]/nt + (1.0-alpha)*Mom(ansl)[j];
	    Mom(ansr)[j] = alpha*Mom(str)[j]/nt + (1.0-alpha)*Mom(ansr)[j];
	}
	Energy(ansl) = alpha*Energy(stl)/nt + (1.0-alpha)*Energy(ansl);
	Energy(ansr) = alpha*Energy(str)/nt + (1.0-alpha)*Energy(ansr);

	/*reset obstatcle states. */
	if(is_obstacle_state(ansl))
	    obstacle_state(intfc, ansl, sizest);
	if(is_obstacle_state(ansr))
	    obstacle_state(intfc, ansr, sizest);

	reset_gamma(ansl);
	reset_gamma(ansr);
	
	return YES;
}


