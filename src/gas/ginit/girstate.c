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
*				girstate.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains initialization routines for generating random state
*	fields.
*/

#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	void	allocate_correlated_state_array(RANDOM_STATE*);
LOCAL	void	allocate_independent_state_array(RANDOM_STATE*);
LOCAL	void	allocate_random_state_storage_arrays(RANDOM_STATE*);
LOCAL	void	copy_correlated_state_array(RANDOM_STATE*,RANDOM_STATE*);
LOCAL	void	copy_independent_state_array(RANDOM_STATE*,RANDOM_STATE*);
LOCAL	void	copy_random_levels(Locstate***,Locstate***,RANDOM_STATE*);
LOCAL	void	free_random_region_comp_type(COMP_TYPE*);
LOCAL	void	get_state_random_region(double*,Locstate,COMP_TYPE*,
	        			HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	gram_schmidt(double**,int);
LOCAL	void	init_random_correlation_ellipsoid(RANDOM_STATE*);
LOCAL	void	init_random_seed(unsigned short int*);
LOCAL	void	inlet_bounding_box(RANDOM_STATE*,HYPER_SURF*);
LOCAL	void	prompt_for_standard_deviation_state(Locstate,Gas_param*,int);
LOCAL	void	set_indep_to_correlation_weight(RANDOM_STATE*);
LOCAL	void	set_radial_velocity_decay_parameters(RANDOM_STATE*,INTERFACE*);
LOCAL	void	set_up_random_evaluation_grid(RANDOM_STATE*,RECT_GRID*);


EXPORT	void	set_random_region_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	RANDOM_STATE	*rstate;

	if (comp_type->type == RANDOM_REGION)
	    return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = RANDOM_REGION;
	rstate = allocate_random_state_structure(front->interf);
	comp_type->extra = (POINTER)rstate;
	comp_type->_get_state = get_state_random_region;
	comp_type->free_comp_type_extra = free_random_region_comp_type;
}		/*end set_random_region_comp_type*/


EXPORT	void	init_random_state_region(
	Locstate	mean,
	RANDOM_STATE	*rstate,
	Front		*front)
{
	IMPORT		boolean suppress_prompts;
	char		s[Gets_BUF_SIZE];
	double		*L, *U;
	int		i, stype;
	int		dim = front->rect_grid->dim;

	init_random_seed(rstate->xsubi);
	stype = state_type(mean);
	stype = init_state_type(stype);
	set_state(Mean(rstate),stype,mean);

	verbose_fprint_state(stdout,"Mean state",Mean(rstate));
	if (suppress_prompts == NO)
	    verbose_fprint_state(stderr,"Mean state",Mean(rstate));

	/* Set up the standard deviation state */
	prompt_for_standard_deviation_state(Sigma(rstate),
					    Params(Mean(rstate)),stype);

	/* Set up bounding box */
	rstate->grid.dim = dim;
	L = rstate->grid.L;
	U = rstate->grid.U;

	screen("Enter a bounding box for the random region (");
	for (i = 0; i < dim; ++i)
	    screen("L[%d] ",i);
	for (i = 0; i < dim-1; ++i)
	    screen("U[%d] ",i);
	screen("U[%d])\n",dim-1);
	screen("\tthis box will be extended periodically to encompass\n");
	screen("\tany actual coordinates contained in the random region\n");
	for (i = 0; i < dim; ++i)
	{
	    L[i] = front->rect_grid->L[i];
	    U[i] = front->rect_grid->U[i];
	}
	screen("Enter L, U (default =");
	for (i = 0; i < dim; ++i)
	    screen(" %g",L[i]);
	for (i = 0; i < dim; ++i)
	    screen(" %g",U[i]);
	screen("): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    char       *c;
	    static const char *seps = " \t";
	    c = strtok(s,seps);
	    for (i = 0; (c!=NULL) && (i<dim); c = strtok(NULL,seps), ++i)
	        (void) sscan_float(c,L+i);
	    if (i < dim)
	    {
	        screen("ERROR in init_random_state_region(), ");
	        screen("invalid input of bounding box\n");
	        clean_up(ERROR);
	    }
	    for (i = 0; (c!=NULL) && (i<dim); c = strtok(NULL,seps), ++i)
	        (void) sscan_float(c,U+i);
	    if (i < dim)
	    {
	        screen("ERROR in init_random_state_region(), ");
	        screen("invalid input of bounding box\n");
	        clean_up(ERROR);
	    }
	}

	init_random_correlation_ellipsoid(rstate);

	rstate->N = 10;
	screen("Enter the number of evalutions per lambda (dflt = %d): ",
		rstate->N);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&rstate->N);
	if (rstate->N <= 0)
	    rstate->N = 10;

	set_up_random_evaluation_grid(rstate,front->rect_grid);

	allocate_correlated_state_array(rstate);
	allocate_independent_state_array(rstate);

	set_indep_to_correlation_weight(rstate);

	generate_random_region(1.0,rstate->old_st,rstate,front->interf);
	free(rstate->indep_st);
	rstate->indep_st = NULL;
	free(rstate->indep_st_store);
	rstate->indep_st_store = NULL;
}		/*end init_random_state_region*/

/*ARGSUSED*/
LOCAL	void	get_state_random_region(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	RANDOM_STATE	*rstate;

	debug_print("init_states","Entered get_state_random_region()\n");
	if (ct->type != RANDOM_REGION)
	{
	    screen("ERROR in get_state_random_region(), "
	           "invalid comp_type\n");
	    clean_up(ERROR);
	}
	rstate = Random_state(ct);
	random_region_state(s,coords,rstate);
	set_state(s,stype,s);
}		/*end get_state_random_region*/

EXPORT	void	random_region_state(
	Locstate	s,
	double		*coords,
	RANDOM_STATE	*rstate)
{
	RECT_GRID	*gr = &rstate->grid;
	double		shft, crds[3];
	double		f[2][3];
	double		alpha[8];
	double		*L = gr->L, *U = gr->U;
	double		*h = gr->h;
	Locstate	states[8];
	int		ic[3];
	int		i, k, dim;
	int		l, m, n;
	static const double tol = 0.001;	/*TOLERANCE*/

	dim = gr->dim;

	/* Periodically translate coords so that crds is inside bbox*/
	for (i = 0; i < dim; ++i)
	{
	    crds[i] = coords[i];
	    if (crds[i] < L[i] - tol*h[i])
	    {
	        shft = ceil((L[i] - crds[i])/(U[i] - L[i]));
	        crds[i] += shft*(U[i] - L[i]);
	    }
	    else if (crds[i] < L[i])
	        crds[i] = L[i];
	    if (crds[i] > U[i] + tol*h[i])
	    {
	        shft = -floor((crds[i] - U[i])/(U[i] - L[i]));
	        crds[i] += shft*(U[i] - L[i]);
	    }
	    else if (crds[i] > U[i])
	        crds[i] = U[i];
	}
	for (i = 0; i < 3; ++i)
	    ic[i] = 0;
	if (rect_in_which(crds,ic,gr) == FUNCTION_FAILED)
	{
	    screen("ERROR in random_region_state(), "
	           "crds outside of grid\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < dim; ++i)
	{
	    if (ic[i] >= gr->gmax[i])
	        ic[i] = gr->gmax[i] - 1;
	    f[1][i] = (crds[i] - cell_edge(ic[i],i,gr))/cell_width(ic[i],i,gr);
	    f[0][i] = 1.0 - f[1][i];
	}
	for (i = dim; i < 3; ++i)
	{
	    f[1][i] = 0.0;
	    f[0][i] = 1.0 - f[1][i];
	}

	k = 1<<dim;
	for (i = 0; i < k; ++i)
	{
	    l = i%2;
	    m = (i/2)%2;
	    n = (i/4)%2;
	    states[i] = rstate->old_st[ic[0]+l][ic[1]+m][ic[2]+n];
	    alpha[i] = f[l][0]*f[m][1]*f[n][2];
	}
	g_linear_combination_of_states(alpha,states,k,s);
}		/*end random_region_state*/

LOCAL	void	free_random_region_comp_type(
	COMP_TYPE	*comp_type)
{
	RANDOM_STATE	*rstate;

	if (comp_type->type != RANDOM_REGION)
	    return;
	rstate = Random_state(comp_type);

	free_random_state_structure(rstate);

	comp_type->extra = NULL;
}		/*end free_random_region_comp_type*/

/*
*		prompt_for_random_flow_inlet():
*
*	Prompts for a random perturbation of a fixed state obtained
*	by an adjacent ambient region.  The requested data consists
*	of a dimensionless standard deviation for each state variable together
*	with an optional random number generator seed.
*
*	The perturbations of the state variables are of the form
*	y = Y,  where Y is N(mu,sigma^2).  Positivity constraints
*	on the state variables will cause the discarding of any data
*	that produces a negative density or pressure,  or that results
*	in a reversal of the velocity direction.  Thus the value of
*	sigma should generally be small relative to mu.
*/

EXPORT	int	prompt_for_random_flow_inlet(
	double		*coords,
	COMPONENT	comp,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	int		index,
	INIT_DATA	*init)
{
	BOUNDARY_STATE	Bstate;
	RANDOM_STATE	*rstate;
	char		s[Gets_BUF_SIZE];
	double		q, v[3], l;
	int		i,  dim = intfc->dim;
	int		stype = TGAS_STATE;

	rstate = allocate_random_state_structure(intfc);

	init_random_seed(rstate->xsubi);

	stype = init_state_type(stype);

	/* Set up mean state */
	set_ambient_boundary_state(Mean(rstate),coords,comp,intfc,init);
	set_state(Mean(rstate),stype,Mean(rstate));

	/* Set up the standard deviation state */
	prompt_for_standard_deviation_state(Sigma(rstate),
					    Params(Mean(rstate)),stype);

	inlet_bounding_box(rstate,hs);

	set_radial_velocity_decay_parameters(rstate,intfc);

	init_random_correlation_ellipsoid(rstate);

	for (i = 0; i < dim; ++i)
	    v[i] = vel(i,Mean(rstate));
	q = mag_vector(v,dim);
	if (q == 0.0)
	    q = 1.0;
	l = mag_vector(rstate->lambda,dim);
	rstate->tau = l/q;
	screen("Enter the correlation time tau (dflt = %g): ",rstate->tau);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&rstate->tau);

	rstate->N = 10;
	screen("Enter the number of evalutions per tau "
	       "and lambda (dflt = %d): ",rstate->N);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&rstate->N);
	if (rstate->N <= 0)
	    rstate->N = 10;

	set_up_random_evaluation_grid(rstate,computational_grid(intfc));

	allocate_random_state_storage_arrays(rstate);

	rstate->delta_t = rstate->tau/rstate->N;
	rstate->tlast = rstate->delta_t; /* Initial time + delta_t */

	set_indep_to_correlation_weight(rstate);

	generate_random_region(1.0,rstate->save_st,rstate,intfc);

	for (i = 0; i < 3; ++i)
	    rstate->seed_after_save[i] = rstate->xsubi[i];
	copy_random_levels(rstate->old_st,rstate->save_st,rstate);
	for (i = 0; i < 3; ++i)
	    rstate->seed_after_old[i] = rstate->xsubi[i];
	generate_random_region(sqrt(2.0*rstate->N-1.0),
	                       rstate->new_st,rstate,intfc);
	relax_random_inlet_level(rstate->old_st,rstate->new_st,rstate);

	Bstate._boundary_state_function = random_velocity_inlet;
	Bstate._boundary_state_function_name = strdup("random_velocity_inlet");
	Bstate._boundary_state = NULL;
	Bstate._boundary_state_data = (POINTER) rstate;
	Bstate._fprint_boundary_state_data =
	        		g_fprint_random_velocity_inlet_data;
	index = add_bstate_to_list(&Bstate,intfc,index);
	return index;
}		/*end prompt_for_random_flow_inlet*/

LOCAL	void set_radial_velocity_decay_parameters(
	RANDOM_STATE *rstate,
	INTERFACE    *intfc)
{
	RECT_GRID *rgr;
	char	  s[Gets_BUF_SIZE];

	if (axisymmetric_random_region_about_origin(rstate,intfc) == NO)
	    return;

	rgr = &rstate->grid;
	RadialVelocityDecayExponent(rstate) = 2.0;
	RadialVelocityDecayScale(rstate) = 0.1*(rgr->U[0] - rgr->L[0]);

	screen_print_long_string(
	    "The radial component of velocity near R = 0, is scaled to zero "
	    "by multipying it by a factor (1.0 - exp(-(r/R)^alpha)) to "
	    "ensure that axisymmetry is established near the origin.\n");

	screen("Enter the radial velocity decay exponent alpha (dflt = %g): ",
	       RadialVelocityDecayExponent(rstate));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&RadialVelocityDecayExponent(rstate));
	if (RadialVelocityDecayExponent(rstate) < 0.0)
	{
	    screen("ERROR in set_radial_velocity_decay_parameters(), "
		   "the the radial velocity decay exponent must "
		   "be non-negative\n");
	    clean_up(ERROR);
	}

	screen("Enter the radial velocity decay length scale R (dflt = %g): ",
	       RadialVelocityDecayScale(rstate));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&RadialVelocityDecayScale(rstate));
	if (RadialVelocityDecayScale(rstate) < 0.0)
	{
	    screen("ERROR in set_radial_velocity_decay_parameters(), "
		   "the the radial velocity decay scale must "
		   "be non-negative\n");
	    clean_up(ERROR);
	}
}		/*end set_radial_velocity_decay_parameters*/

LOCAL	void	init_random_seed(
	unsigned short int      *xsubi)
{
	char		s[Gets_BUF_SIZE];

	/* Set up random number seeds */
	xsubi[0] = 0x2ef2;
	xsubi[1] = 0x8705;
	xsubi[2] = 0x09c9;
	screen("Enter an optional set of three random seeds ");
	screen("for the random number generator\n");
	screen("\t(default = %d %d %d): ",xsubi[0],xsubi[1],xsubi[2]);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%hu %hu %hu",xsubi,xsubi+1,xsubi+2);
}		/*end init_random_seed*/


EXPORT	void	read_print_random_velocity_inlet_data(
	INIT_DATA      *init,
	const IO_TYPE   *io_type,
	BOUNDARY_STATE *bstate,
	INTERFACE      *intfc)
{
	FILE            *file = io_type->file;
	RANDOM_STATE	*rstate;
	Gas_param       **prms_lst;
	char		c, state_type_string[80];
	int		stype;
	int		nprms;
	int		i, j, k, dim = intfc->dim;
	int		gmax[3];
	size_t		sizest = size_of_state(intfc);
	uint64_t	*old_prms;

	bstate->_boundary_state_function = random_velocity_inlet;
	bstate->_fprint_boundary_state_data =
	        	g_fprint_random_velocity_inlet_data;
	rstate = allocate_random_state_structure(intfc);
	bstate->_boundary_state_data = (POINTER) rstate;

	nprms = read_Gas_param_list(init,io_type,&old_prms,&prms_lst,
	                            dim,sizest,NO);

	if (fgetstring(file,"State type = ") == FUNCTION_FAILED)
	{
	    screen("ERROR in read_print_random_velocity_inlet_data(), "
		   "invalid restart file, State type not found\n");
	    clean_up(ERROR);
	}
	(void) fscanf(file,"%s",state_type_string);
	stype = g_read_state_type_from_string(state_type_string);
	if (fgetstring(file,"Mean state = ") == FUNCTION_FAILED)
	{
	    screen("ERROR in read_print_random_velocity_inlet_data(), "
		   "invalid restart file, Mean state data not found\n");
	    clean_up(ERROR);
	}
	read_print_gas_data(io_type,&Mean(rstate),NO,
	                    stype,sizest,nprms,old_prms,prms_lst,dim);
	if (fgetstring(file,"Sigma state = ") == FUNCTION_FAILED)
	{
	    /* Try old style printout */
	    if (fgetstring(file,"Variance state = ") == FUNCTION_FAILED)
	    {
	        screen("ERROR in read_print_random_velocity_inlet_data(), "
		       "invalid restart file, standard deviation data "
		       "not found\n");
	        clean_up(ERROR);
	    }
	}
	read_print_gas_data(io_type,&Sigma(rstate),NO,
	                    stype,sizest,nprms,old_prms,prms_lst,dim);

	if (fgetstring(file,"Q = ") == FUNCTION_FAILED)
	{
	    screen("ERROR in read_print_random_velocity_inlet_data(), "
		   "invalid restart file, Q not found\n");
	    clean_up(ERROR);
	}
	if ((c = getc(file)) != '\f')   /* NOBINARY */
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < 3; ++i)
	    for (j = 0; j < 3; ++j)
	    	(void) fscan_float(file,&rstate->Q[i][j]);
	}
	else
	{
	    (void) getc(file);
	    for (i = 0; i < 3; ++i)
	        (void) read_binary_real_array(rstate->Q[i],3,io_type);
	}

	(void) fgetstring(file,"= ");
	if ((c = getc(file)) != '\f')   /* NOBINARY */
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < 3; ++i)
	    	(void) fscan_float(file,rstate->lambda+i);
	}
	else
	{
	    (void) getc(file);
	    (void) read_binary_real_array(rstate->lambda,3,io_type);
	}
	rstate->correlated = NO;
	for (i = 0; i < 3; ++i)
	{
	    if (rstate->lambda[i] != 0.0)
	    {
	    	rstate->correlated = YES;
	    	break;
	    }
	}


	(void) fgetstring(file,"= ");
	if ((c = getc(file)) != '\f')   /* NOBINARY */
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < 3; ++i)
	    	for (j = 0; j < 3; ++j)
	    	    (void) fscan_float(file,&rstate->A[i][j]);
	}
	else
	{
	    (void) getc(file);
	    for (i = 0; i < 3; ++i)
	        (void) read_binary_real_array(rstate->A[i],3,io_type);
	}


	(void) fgetstring(file,"= ");
	(void) fscanf(file,"%d",&rstate->N);

	rstate->tau = fread_float("= ",io_type);

	rstate->delta_t = rstate->tau/rstate->N;

	(void) fgetstring(file,"Evaluation grid =");
	read_rectangular_grid(io_type,&rstate->grid,YES,
			      &computational_grid(intfc)->Remap);

	set_indep_to_correlation_weight(rstate);

	(void) fgetstring(file,"= ");
	rstate->time_of_save = fread_float("= ",io_type);
	(void) fgetstring(file,"= ");
	(void) fscanf(file,"%hu %hu %hu",rstate->seed_after_save,
	                                 rstate->seed_after_save+1,
	                                 rstate->seed_after_save+2);

	for (i = 0; i < 3; ++i)
	    gmax[i] = rstate->grid.gmax[i]+1;

	allocate_random_state_storage_arrays(rstate);

	(void) fgetstring(file,"saved random states");
	for (i = 0; i < gmax[0]; ++i)
	for (j = 0; j < gmax[1]; ++j)
	for (k = 0; k < gmax[2]; ++k)
	    read_print_gas_data(io_type,&rstate->save_st[i][j][k],NO,
	        	        stype,sizest,nprms,old_prms,prms_lst,dim);

	if (axisymmetric_random_region_about_origin(rstate,intfc) == YES)
	{
	    if (fgetstring(file,"Radial velocity decay exponent = "))
	        RadialVelocityDecayExponent(rstate) = fread_float(NULL,io_type);
	    if (fgetstring(file,"Radial velocity decay length scale = "))
	        RadialVelocityDecayScale(rstate) = fread_float(NULL,io_type);
	}

	rstate->tlast = rstate->time_of_save + rstate->delta_t;

	for (i = 0; i < 3; ++i)
	    rstate->xsubi[i] = rstate->seed_after_save[i];
	if (rstate->time_of_save < rstate->delta_t)
	{
	    copy_random_levels(rstate->old_st,rstate->save_st,rstate);
	}
	else
	{
	    generate_random_region(sqrt(2.0*rstate->N-1.0),
	        	           rstate->old_st,rstate,intfc);
	    relax_random_inlet_level(rstate->save_st,rstate->old_st,rstate);
	}
	for (i = 0; i < 3; ++i)
	    rstate->seed_after_old[i] = rstate->xsubi[i];
	generate_random_region(sqrt(2.0*rstate->N-1.0),
			       rstate->new_st,rstate,intfc);
	relax_random_inlet_level(rstate->old_st,rstate->new_st,rstate);
}		/*end read_print_random_velocity_inlet_data*/

LOCAL	void	prompt_for_standard_deviation_state(
	Locstate	sigma,
	Gas_param	*params,
	int		stype)
{
	char       s[Gets_BUF_SIZE];
	int        i, dim = params->dim;
	static const char *xyz[] = { "x", "y", "z"};

	set_type_of_state(sigma,stype);
	Init_params(sigma,params);
	screen("Enter the density perturbation standard ");
	screen("deviation (dflt = %g): ",Dens(sigma));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Dens(sigma));
	switch (stype)
	{
	case GAS_STATE:
	    screen("Enter the total energy perturbation ");
	    screen("standard deviation (dflt = %g): ",Energy(sigma));
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscan_float(s,&Energy(sigma));
	    break;
	case TGAS_STATE:
	    screen("Enter the pressure perturbation ");
	    screen("standard deviation (dflt = %g): ",Press(sigma));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscan_float(s,&Press(sigma));
	    break;
	case EGAS_STATE:
	    screen("Enter the specific internal energy perturbation ");
	    screen("standard deviation (dflt = %g): ",Energy(sigma));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscan_float(s,&Energy(sigma));
	    break;
	default:
	    screen("ERROR in prompt_for_standard_deviation_state(), ");
	    fprint_state_type(stderr,
	        "unknown or unsupported state type = ",stype);
	    fprint_state_type(stdout,
	        "unknown or unsupported state type = ",stype);
	    clean_up(ERROR);
	}
	switch (stype)
	{
	case GAS_STATE:
	    for (i = 0; i < dim; ++i)
	    {
	        screen("Enter the %s momentum perturbation ",xyz[i]);
	        screen("standard deviation (dflt = %g): ",Mom(sigma)[i]);
	        (void) Gets(s);
	        if (s[0] != '\0')
	            (void) sscan_float(s,Mom(sigma)+i);
	    }
	    break;
	case TGAS_STATE:
	case EGAS_STATE:
	    for (i = 0; i < dim; ++i)
	    {
	        screen("Enter the %s velocity perturbation ",xyz[i]);
	        screen("standard deviation (dflt = %g): ",Vel(sigma)[i]);
	        (void) Gets(s);
	        if (s[0] != '\0')
	            (void) sscan_float(s,Vel(sigma)+i);
	    }
	    break;
	default:
	    break;
	}
}		/*end prompt_for_standard_deviation_state*/

LOCAL	void copy_random_levels(
	Locstate ***save_st,
	Locstate ***old_st,
	RANDOM_STATE *rstate)
{
	RECT_GRID	*grid = &rstate->grid;
	int		i, j, k;
	int		al[3];
	size_t		sizest = Params(Mean(rstate))->sizest;
	int		dim = grid->dim;

	rstate->time_of_save = rstate->tlast - rstate->delta_t;

	for (i = 0; i < dim; ++i)
	    al[i] = grid->gmax[i]+1;
	for (i = dim; i < 3; ++i)
	    al[i] = 1;

	for (i = 0; i < al[0]; ++i)
	for (j = 0; j < al[1]; ++j)
	for (k = 0; k < al[2]; ++k)
	    ft_assign(save_st[i][j][k],old_st[i][j][k],sizest);
}		/*end copy_random_levels*/

LOCAL	void	set_indep_to_correlation_weight(
	RANDOM_STATE	*rstate)
{
	int M_ell=0;
	int *lbuf = rstate->grid.lbuf;
	int *ubuf = rstate->grid.ubuf;
	int ll, nn, mm;

	for (ll = -lbuf[0]; ll <= ubuf[0]; ++ll)
	for (nn = -lbuf[1]; nn <= ubuf[1]; ++nn)
	for (mm = -lbuf[2]; mm <= ubuf[2]; ++mm)
	{
	    if (in_correlation_ellipse(ll,mm,nn,rstate))
	    	++M_ell;
	}
	rstate->M_ell = M_ell;
}		/*end set_indep_to_correlation_weight*/

LOCAL	void	init_random_correlation_ellipsoid(
	RANDOM_STATE *rstate)
{
	RECT_GRID  *grid = &rstate->grid;
	char	   s[Gets_BUF_SIZE];
	double	   rad;
	double	   *lambda = rstate->lambda;
	double	   **Q;
	int        i, j, dim = grid->dim;
	static const char *ordinal[] = {"first","second","third"};
	static const char *cardinal[] = {"one","two","three"};
	static const char *scan_3_floats = "%lf %lf %lf";

	rad = grid->U[0] - grid->L[0];
	for (i = 1; i < dim; ++i)
	    rad = min(rad,grid->U[i] - grid->L[i]);
	for (i = 0; i < dim; ++i)
	    lambda[i] = 0.1*rad;
	for (i = dim; i < 3; ++i)
	    lambda[i] = 0.0;
	Q = rstate->Q;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	    Q[i][j] = (i == j) ? 1.0 : 0.0;
	for (i = 0; i < dim; ++i)
	{
	    screen("Enter the %s "
	           "axis vector for the correlation ellipsoid\n\t"
	           "default = (%g",ordinal[i],Q[i][0]);
	    for (j = 1; j < dim; ++j)
	        screen(", %g",Q[i][j]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscanf(s,scan_3_floats,Q[i],Q[i]+1,Q[i]+2);
	    	lambda[i] = mag_vector(Q[i],3);
	    	gram_schmidt(Q,3);
	    }
	}
	screen("Enter the%s%s correlation length%s ",
	        (dim > 1) ? " " : "",
	        (dim > 1) ? cardinal[dim-1] : "",
	        (dim > 1) ? "s" : "");
	screen("(default = (%g",lambda[0]);
	for (j = 1; j < dim; ++j)
	    screen(", %g",lambda[j]);
	screen(")): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscanf(s,scan_3_floats,lambda,lambda+1,lambda+2);
	    for (i = 0; i < dim; ++i)
	    	lambda[i] = fabs(lambda[i]);
	}

}		/*end init_random_correlation_ellipsoid*/

LOCAL	void	allocate_random_state_storage_arrays(
	RANDOM_STATE *rstate)
{
	RECT_GRID	*grid = &rstate->grid;
	byte		*st;
	int		dim = grid->dim;
	int		al[3];
	size_t		sizest = Params(Mean(rstate))->sizest;
	int		i, j, k;

	for (i = 0; i < dim; ++i)
	    al[i] = grid->gmax[i]+1;
	for (i = dim; i < 3; ++i)
	    al[i] = 1;
	tri_array(&rstate->corr_st[0],al[0],al[1],al[2],sizeof(Locstate));
	tri_array(&rstate->corr_st[1],al[0],al[1],al[2],sizeof(Locstate));
	tri_array(&rstate->corr_st[2],al[0],al[1],al[2],sizeof(Locstate));
	uni_array(&rstate->corr_st_store,3*al[0]*al[1]*al[2],sizest);

	for (st = rstate->corr_st_store, i = 0; i < al[0]; ++i)
	for (j = 0; j < al[1]; ++j)
	for (k = 0; k < al[2]; ++k)
	{
	    rstate->corr_st[0][i][j][k] = (Locstate) st;
	    st += sizest;
	    rstate->corr_st[1][i][j][k] = (Locstate) st;
	    st += sizest;
	    rstate->corr_st[2][i][j][k] = (Locstate) st;
	    st += sizest;
	}
	rstate->old_st = rstate->corr_st[0];
	rstate->new_st = rstate->corr_st[1];
	rstate->save_st = rstate->corr_st[2];

	allocate_independent_state_array(rstate);
}		/*end allocate_random_state_storage_arrays*/

LOCAL	void	allocate_correlated_state_array(
	RANDOM_STATE *rstate)
{
	RECT_GRID	*grid = &rstate->grid;
	byte		*st;
	int		dim = grid->dim;
	int		al[3];
	size_t		sizest = Params(Mean(rstate))->sizest;
	int		i, j, k;

	for (i = 0; i < dim; ++i)
	    al[i] = grid->gmax[i]+1;
	for (i = dim; i < 3; ++i)
	    al[i] = 1;
	tri_array(&rstate->corr_st[0],al[0],al[1],al[2],sizeof(Locstate));
	uni_array(&rstate->corr_st_store,al[0]*al[1]*al[2],sizest);

	for (st = rstate->corr_st_store, i = 0; i < al[0]; ++i)
	for (j = 0; j < al[1]; ++j)
	for (k = 0; k < al[2]; ++k)
	{
	    rstate->corr_st[0][i][j][k] = (Locstate) st;
	    st += sizest;
	}
	rstate->old_st = rstate->corr_st[0];
}		/*end allocate_correlated_state_array*/

LOCAL	void	copy_correlated_state_array(
	RANDOM_STATE *copy,
	RANDOM_STATE *rstate)
{
	RECT_GRID	*grid = &rstate->grid;
	int		dim = grid->dim;
	int		al[3];
	size_t		sizest = Params(Mean(rstate))->sizest;
	int		i, j, k;

	for (i = 0; i < dim; ++i)
	    al[i] = grid->gmax[i]+1;
	for (i = dim; i < 3; ++i)
	    al[i] = 1;

	if (rstate->old_st != NULL)
	{
	    for (i = 0; i < al[0]; ++i)
	    for (j = 0; j < al[1]; ++j)
	    for (k = 0; k < al[2]; ++k)
	        ft_assign(copy->old_st[i][j][k],rstate->old_st[i][j][k],sizest);
	}
	if (rstate->new_st != NULL)
	{
	    for (i = 0; i < al[0]; ++i)
	    for (j = 0; j < al[1]; ++j)
	    for (k = 0; k < al[2]; ++k)
	        ft_assign(copy->new_st[i][j][k],rstate->new_st[i][j][k],sizest);
	}
	if (rstate->save_st != NULL)
	{
	    for (i = 0; i < al[0]; ++i)
	    for (j = 0; j < al[1]; ++j)
	    for (k = 0; k < al[2]; ++k)
	        ft_assign(copy->save_st[i][j][k],rstate->save_st[i][j][k],sizest);
	}
}		/*end allocate_correlated_state_array*/

LOCAL	void	allocate_independent_state_array(
	RANDOM_STATE *rstate)
{
	RECT_GRID *grid = &rstate->grid;
	byte		*st;
	int		dim = grid->dim;
	int		al[3];
	size_t		sizest = Params(Mean(rstate))->sizest;
	int		i, j, k;

	for (i = 0; i < dim; ++i)
	    al[i] = grid->gmax[i]+1+grid->lbuf[i]+grid->ubuf[i];
	for (i = dim; i < 3; ++i)
	    al[i] = 1;

	tri_array(&rstate->indep_st,al[0],al[1],al[2],sizeof(Locstate));
	uni_array(&rstate->indep_st_store,al[0]*al[1]*al[2],sizest);

	for (st = rstate->indep_st_store, i = 0; i < al[0]; ++i)
	for (j = 0; j < al[1]; ++j)
	for (k = 0; k < al[2]; ++k)
	{
	    rstate->indep_st[i][j][k] = (Locstate) st;
	    st += sizest;
	}
}		/*end allocate_independent_state_array*/

LOCAL	void	copy_independent_state_array(
	RANDOM_STATE *copy,
	RANDOM_STATE *rstate)
{
	RECT_GRID *grid = &rstate->grid;
	int		dim = grid->dim;
	int		al[3];
	size_t		sizest = Params(Mean(rstate))->sizest;
	int		i, j, k;

	for (i = 0; i < dim; ++i)
	    al[i] = grid->gmax[i]+1+grid->lbuf[i]+grid->ubuf[i];
	for (i = dim; i < 3; ++i)
	    al[i] = 1;

	for (i = 0; i < al[0]; ++i)
	for (j = 0; j < al[1]; ++j)
	for (k = 0; k < al[2]; ++k)
	{
	    ft_assign(copy->indep_st[i][j][k],rstate->indep_st[i][j][k],sizest);
	}
}		/*end copy_independent_state_array*/

LOCAL	void	set_up_random_evaluation_grid(
	RANDOM_STATE *rstate,
	RECT_GRID    *rgr)
{
	RECT_GRID *grid = &rstate->grid;
	double	*lambda = rstate->lambda;
	double	**Q = rstate->Q;
	double	**A = rstate->A;
	double	D[3];
	double	lambda_min = HUGE_VAL;
	double	*L = grid->L, *U = grid->U;
	double	dx;
	int	i, j, k, *gmax = grid->gmax;
	int	*lbuf = grid->lbuf;
	int	*ubuf = grid->ubuf;
	int	dim = grid->dim;

	rstate->correlated = YES;
	for (i = 0; i < dim; ++i)
	{
	    if (lambda[i] > 0.0 && lambda[i] < lambda_min)
	    	lambda_min = lambda[i];
	}
	if (lambda_min == HUGE_VAL)
	{
	    rstate->correlated = NO;
	    for (i = 0; i < dim; ++i)
	    {
	    	if ((U[i]-L[i]) > 0.0 && (U[i]-L[i]) < lambda_min)
	    	    lambda_min = U[i]-L[i];
	    }
	    if (lambda_min == HUGE_VAL)
	    {
	    	screen("ERROR in set_up_random_evaluation_grid(), "
	    	       "all scales are zero\n");
	    	clean_up(ERROR);
	    }
	}
	dx = lambda_min/rstate->N;
	for (i = 0; i < dim; ++i)
	{
	    double delta;
	    gmax[i] = (int) ceil((U[i]-L[i])/dx);
	    delta = 0.5*(gmax[i]*dx - (U[i]-L[i]));
	    U[i] += delta;
	    L[i] -= delta;
	}

	for (i = 0; i < dim; ++i)
	{
	    double buf_len;
	    for (buf_len = 0.0, j = 0; j < dim; ++j)
	    	buf_len += sqr(lambda[j])*sqr(Q[i][j]);
	    buf_len = sqrt(buf_len);
	    ubuf[i] = lbuf[i] = (int) ceil(buf_len/dx);
	}

	set_rect_grid(L,U,L,U,lbuf,ubuf,gmax,dim,&rgr->Remap,grid);

	for (i = 0; i < 3; ++i)
	    D[i] = (lambda[i] > 0.1*dx/rstate->N) ?
	    	   1.0/sqr(lambda[i]) : 100.0*sqr(rstate->N)/(dx*dx);

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	    for (A[i][j] = Q[0][i]*D[0]*Q[0][j], k = 1; k < 3; ++k)
	        A[i][j] += Q[k][i]*D[k]*Q[k][j];
}		/*end set_up_random_evaluation_grid*/

LOCAL	void inlet_bounding_box(
	RANDOM_STATE *rstate,
	HYPER_SURF *hs)
{
	RECT_GRID *gr = &rstate->grid;
	BOND *b;
	double *L = gr->L, *U = gr->U;
	int i, dim = gr->dim = hs->interface->dim;

	switch (dim)
	{
#if defined(ONED)
	case 1:
	    L[0] = Coords(Point_of_hs(hs))[0];
	    U[0] = Coords(Point_of_hs(hs))[0];
	    break;
#endif /* defined(ONED) */
	case 2:
	    b = Curve_of_hs(hs)->first;
	    for (i = 0; i < dim; ++i)
	    {
	    	L[i] = Coords(b->start)[i];
	    	U[i] = Coords(b->start)[i];
	    }

	    for (; b != NULL; b = b->next)
	    for (i = 0; i < dim; ++i)
	    {
	    	L[i] = min(L[i],Coords(b->end)[i]);
	    	U[i] = max(U[i],Coords(b->end)[i]);
	    }
	    break;
#if defined(THREED)
	case 3:
	    {
	        POINT *p;
	        SURFACE *s = Surface_of_hs(hs);
	        TRI *tri;
	        int j;
    
	        tri = first_tri(s);
	        p = Point_of_tri(tri)[0];
	        for (i = 0; i < dim; ++i)
	        {
	    	    L[i] = Coords(p)[i];
	    	    U[i] = Coords(p)[i];
	        }

	        for (; !at_end_of_tri_list(tri,s); tri = tri->next)
	        for (j = 0; j < 3; ++j)
	        for (p = Point_of_tri(tri)[j], i = 0; i < dim; ++i)
	        {
	    	    L[i] = min(L[i],Coords(p)[i]);
	    	    U[i] = max(U[i],Coords(p)[i]);
	        }
	    }
	    break;
#endif /* defined(THREED) */
	}
}		/*end inlet_bounding_box*/

LOCAL	void	gram_schmidt(
	double **Q,
	int dim)
{
	int i, j, k;
	double len;
	double sp;

	for (i = 0; i < dim; ++i)
	{
	    len = mag_vector(Q[i],dim);
	    for (j = 0; j < dim; ++j)
	        Q[i][j] /= len;
	    for (j = i+1; j < dim; ++j)
	    {
	    	sp = scalar_product(Q[i],Q[j],dim);
	    	for (k = 0; k < dim; ++k)
	    	    Q[j][k] -= sp*Q[i][k];
	    }
	}
}		/*end gram_schmidt*/


EXPORT	RANDOM_STATE	*allocate_random_state_structure(
	INTERFACE	*intfc)
{
	RANDOM_STATE	*rstate;
	ALIGN		*algn;
	size_t 		sizest = size_of_state(intfc);
	size_t		size, naRS, naS;

	naRS = num_aligns(sizeof(RANDOM_STATE));
	naS = num_aligns(sizest);

	size = sizeof(ALIGN)*(naRS+2*naS);
	scalar(&algn,size);
	rstate = (RANDOM_STATE*)algn;		algn += naRS;
	Mean(rstate) = (Locstate)algn;		algn += naS;
	Sigma(rstate) = (Locstate)algn;
	rstate->A[0] = &rstate->Astore[0][0];
	rstate->A[1] = &rstate->Astore[1][0];
	rstate->A[2] = &rstate->Astore[2][0];
	rstate->Q[0] = &rstate->Qstore[0][0];
	rstate->Q[1] = &rstate->Qstore[1][0];
	rstate->Q[2] = &rstate->Qstore[2][0];
	clear_state(intfc,Mean(rstate),sizest);
	clear_state(intfc,Sigma(rstate),sizest);
	RadialVelocityDecayExponent(rstate) = 0.0;
	RadialVelocityDecayScale(rstate) = 0.0;
	return rstate;
}		/*end allocate_random_state_structure*/

EXPORT	RANDOM_STATE	*copy_random_state_structure(
	RANDOM_STATE	*rstate,
	INTERFACE	*intfc)
{
	RANDOM_STATE	*copy;
	int             i, j;

	if (rstate == NULL)
	    return NULL;

	copy = allocate_random_state_structure(intfc);
	copy_state(Mean(copy),Mean(rstate));
	copy_state(Sigma(copy),Sigma(rstate));
	for (i = 0; i < 3; ++i)
	{
	    copy->lambda[i] = rstate->lambda[i];
	    for (j = 0; j < 3; ++j)
	    {
		copy->Q[i][j] = rstate->Q[i][j];
		copy->A[i][j] = rstate->A[i][j];
	    }
	    copy->xsubi[i] = rstate->xsubi[i];
	    copy->seed_after_save[i] = rstate->seed_after_save[i];
	    copy->seed_after_old[i] = rstate->seed_after_old[i];
	}
	copy->N = rstate->N;
	copy->tau = rstate->tau;
	copy->tlast = rstate->tlast;
	copy->delta_t = rstate->delta_t;
	copy->M_ell = rstate->M_ell;
	copy->correlated = rstate->correlated;
	copy_rect_grid(&copy->grid,&rstate->grid);
	copy->time_of_save = rstate->time_of_save;
	RadialVelocityDecayExponent(copy) = RadialVelocityDecayExponent(rstate);
	RadialVelocityDecayScale(copy) = RadialVelocityDecayScale(rstate);

	if (rstate->corr_st != NULL)
	{
	    allocate_correlated_state_array(copy);
	    copy_correlated_state_array(copy,rstate);
	}
	if (rstate->indep_st != NULL)
	{
	    allocate_independent_state_array(copy);
	    copy_independent_state_array(copy,rstate);
	}

	return copy;
}		/*end allocate_random_state_structure*/

EXPORT	void	free_random_state_structure(
	RANDOM_STATE	*rstate)
{
	int	i;
	if (rstate == NULL)
	    return;
	if (rstate->indep_st != NULL)
	    free(rstate->indep_st);
	if (rstate->indep_st_store != NULL)
	    free(rstate->indep_st_store);
	for (i = 0; i < 3; ++i)
	    if (rstate->corr_st[i] != NULL)
	    	free(rstate->corr_st[i]);
	if (rstate->corr_st_store != NULL)
	    free(rstate->corr_st_store);

	free(rstate);
}		/*end free_random_state_structure*/
