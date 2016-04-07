/************************************************************************************
frontier is a set of libraries that implements differnt types of front traking algorithms.
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
*				gipert.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Routines for initializing and analyzing Kelvin-Helmholtz
*	and Rayleigh-Taylor sinusoidal perturbation problems.
*
*	Multiple mode interfaces (random surfaces) are also treated.
*/

#if defined(TWOD) || defined(THREED)

#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL	boolean	mis_match(double,double*,POINTER);
LOCAL	double	RT_RS_fab(double,SINE_PERT*,SINE_PERT*,int,int,double);
LOCAL	double	RT_sigma(SINE_PERT*,SINE_PERT*,int,int,double);
LOCAL	int	get_next_mode(int,LIN_PERT*);
LOCAL	int	kh_single_mode(LIN_PERT*);
LOCAL	int	list_modes(NORMAL_MODE**,LIN_PERT*,int,int);
LOCAL	int	rt_single_mode(LIN_PERT*,INIT_DATA*);
LOCAL	void	F_bubble(int,double**,double*,double*,double,double);
LOCAL	void	bksub(int,int,int,int,int,double***);
LOCAL	void	difeq(int,int,double,double**,int*,LIN_PERT*);
LOCAL	void	free_kh_perturbed_comp_type(COMP_TYPE*);
LOCAL	void	free_kh_rt_sine_perturbed_comp_type(COMP_TYPE*);
LOCAL	void	free_rt_perturbed_comp_type(COMP_TYPE*);
LOCAL	void	get_KH_sine_perturbed_state(double*,Locstate,
					    COMP_TYPE*,HYPER_SURF*,INTERFACE*,
					    INIT_DATA*,int);
LOCAL	void	get_state_random_surface_perturbed(double*,Locstate,COMP_TYPE*,
						   HYPER_SURF*,INTERFACE*,
						   INIT_DATA*,int);
LOCAL	void	init_KH_sine_perturbation(SINE_PERT*,SINE_PERT*,int,double);
LOCAL	void	init_random_surface_perturbation(SINE_PERT*,SINE_PERT*,
						 int,double);
LOCAL	void	pinvs(int,int,int,int,int,int,double***,double**);
LOCAL	void	random_surface_perturbed_state(Locstate,double*,
					       double,COMPONENT,int);
LOCAL	void	red(int,int,int,int,int,int,int,int,int,int,int,
		    double***,double**);
LOCAL	void	set_kh_sine_perturbed_comp_type(COMP_TYPE*,Front*,
						int,PERT_TYPE);
LOCAL	void	set_pert_modes(SINE_PERT*,SINE_PERT*,int,int,
			       double**,double*,double*);
LOCAL	void	set_random_surface_perturbed_comp_type(COMP_TYPE*,Front*,
						       int,PERT_TYPE);
LOCAL	void	set_kh_perturbed_comp_type(COMP_TYPE*,Front*);
LOCAL	void	set_rt_perturbed_comp_type(COMP_TYPE*,Front*);
#if defined(TWOD)
LOCAL	void	make_random_curve(int,int,double,double,SINE_PERT*,
				  COMPONENT,COMPONENT,double,boolean,int);
#endif /* defined(TWOD) */


enum _BUBBLE_TYPE {
	MULTI_MODE      = 1,
	MULTI_BUBBLE    = 2,
	RANDOM_BUBBLE   = 3,
	FILE_AMPLITUDES = 4
};
typedef enum _BUBBLE_TYPE BUBBLE_TYPE;


LOCAL const int MAX_NUMBER_MODES = 40;  /* Maximum number of modes for */
					/* statistical regime. */

EXPORT	void init_random_surface(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	INTERFACE	*intfc = front->interf;
	char		s[Gets_BUF_SIZE];
	char		surf_type[Gets_BUF_SIZE];
	double		coords[3];
	double		mean_height;
	double		A, nu, wave_number, phase;
	double		*L = front->rect_grid->GL, *U = front->rect_grid->GU;
	double		rho_b, rho_a;
	double		amb_press;
	double		delta_v[MAXD-1];
	double		xi;
	double		surf_ten;
	int		shear;
	int		min_n, max_n;
	PERT_TYPE	init_type;
	int		number_modes;
	int		i, j, dim = front->interf->dim;
	BUBBLE_TYPE	bubble_type;
	double		**F_wave_number, *F_amplitude, *F_phase;
	double		g_z;
	Gas_param	*params, *params2;
	SINE_PERT	*rsp_a, *rsp_b;
	FOURIER_POLY	Fpoly;
#if defined(TWOD)
	CURVE		*cur;
	boolean            off_tracking;
	int		offset_N_bdry;
#if defined(USE_OVERTURE)
        int             num_points;  
        int             amr_refinecoeff = 1; 
#endif /* if defined(USE_OVERTURE) */ 
#endif /* defined(TWOD) */


	init_type = UNSET_PERT_TYPE;
	screen("\nLight fluid is above and heavy fluid below.\n"
	       "Four types of random surface problems are supported:\n"
	       "\tRayleigh-Taylor random surface with "
	       "linearized initial states (RT)\n"
	       "\tRayleigh-Taylor without linear analysis (RA)\n"
	       "\tKelvin-Helmholtz random surface (KH).\n"
	       "Enter choice: ");
	(void) Gets(s);
	if (strcmp(s,"RT") == 0 || strcmp(s,"rt") == 0)
	    init_type = RT;
	else if (strcmp(s,"RA") == 0 || strcmp(s,"ra") == 0)
	    init_type = RT_AMB;
	else if (strcmp(s,"KH") == 0 || strcmp(s,"kh") == 0)
	    init_type = KELHELM;
	else 
	{
	    screen("ERROR: unrecognized choice for random surface.\n");
	    clean_up(ERROR);
	}

	screen("Enter the choice of initial front description\n"
	       "Four types of descriptions are supported\n"
	       "       Multiple mode description (M, default),\n"
	       "       Multiple bubble description (B),\n"
	       "       Multiple random bubble description (R),\n"
	       "       Read initial amplitudes from a file (F).\n"
	       "       Enter choice: ");
	(void) Gets(surf_type);
	switch (surf_type[0]) 
	{
	case 'B':
	case 'b':
	    bubble_type = MULTI_BUBBLE;
	    number_modes = MAX_NUMBER_MODES;
	    break;
	case 'R':
	case 'r':
	    bubble_type = RANDOM_BUBBLE;
	    number_modes = random_bubble_num_modes(
			"in the initial interface",&min_n,&max_n,dim);
	    break;
	case 'F':
	case 'f':
	    bubble_type = FILE_AMPLITUDES;
	    number_modes = 1;  
	    break;
	case 'M':
	case 'm':
	default: 
	    bubble_type = MULTI_MODE;
	    screen("Enter the number of modes in the initial interface: ");
	    (void) Scanf("%d\n",&number_modes);
	    break;
	}

	if (init_type == KELHELM)
	{
	    set_kh_sine_perturbed_comp_type(comp_type(COMPB),front,
			                    number_modes,init_type);
	    set_kh_sine_perturbed_comp_type(comp_type(COMPA),front,
			                    number_modes,init_type);
	}
	else
	{
	    set_random_surface_perturbed_comp_type(comp_type(COMPB),front,
			                           number_modes,init_type);
	    set_random_surface_perturbed_comp_type(comp_type(COMPA),front,
			                           number_modes,init_type);
	}
	rsp_a = (SINE_PERT*)comp_type(COMPA)->extra;
	rsp_b = (SINE_PERT*)comp_type(COMPB)->extra;


/* For the interface as a whole: */

	screen("Enter the mean position of the front above L[%d]: ",dim-1);
	(void) Scanf("%f\n",&mean_height);
	mean_height += L[dim-1];
	rsp_a->z_intfc = mean_height;
	rsp_b->z_intfc = mean_height;
	rsp_a->z_bdry = U[dim-1];
	rsp_b->z_bdry = L[dim-1];


	/* Now enter the data needed to the SINE_PERT structures 
	   if amplitudes are being read in from a file. */

	rsp_a->read_amplitudes_from_file = 1;
	rsp_b->read_amplitudes_from_file = 1;




/* For multiple bubble description, make Fourier analysis of the front  */

	switch (bubble_type)
	{
	case MULTI_BUBBLE:
	    bi_array(&F_wave_number,number_modes,dim-1,FLOAT);
	    uni_array(&F_phase,number_modes,FLOAT);
	    uni_array(&F_amplitude,number_modes,FLOAT);
	    F_bubble(number_modes,F_wave_number,F_amplitude,F_phase,
		     U[0]-L[0],L[0]);
	    set_pert_modes(rsp_a,rsp_b,dim,number_modes,F_wave_number,
			   F_amplitude,F_phase);

	    free_these(3,F_wave_number,F_amplitude,F_phase);
	    break;

	case RANDOM_BUBBLE:
	    bi_array(&F_wave_number,number_modes,dim-1,FLOAT);
	    uni_array(&F_phase,number_modes,FLOAT);
	    uni_array(&F_amplitude,number_modes,FLOAT);
	    Fpoly.num_modes = number_modes;
	    Fpoly.A = F_amplitude;
	    Fpoly.phase = F_phase;
	    Fpoly.nu = F_wave_number;
	    Fpoly.dim = dim;
	    init_random_modes(0,min_n,max_n,number_modes,&Fpoly,L,U);
	    set_pert_modes(rsp_a,rsp_b,dim,number_modes,F_wave_number,
			   F_amplitude,F_phase);

	    break;

	case FILE_AMPLITUDES:
	    {
	        FILE *READAMPLITUDES;
		double delta_x, delta_y;  
		rsp_a->read_amplitudes_from_file = 0;
		rsp_b->read_amplitudes_from_file = 0;
		screen("Enter the number of initial interface" 
		        "amplitudes in the x direction: ");
		(void) Scanf("%d\n",&(rsp_a->Nx));
		rsp_b->Nx = rsp_a->Nx;
		screen("Enter the number of initial interface" 
		        "amplitudes in the y direction: ");
		(void) Scanf("%d\n",&(rsp_a->Ny));
		rsp_b->Ny = rsp_a->Ny;
		
		delta_x = (front->rect_grid->GU[0] - 
			   front->rect_grid->GL[0])/(rsp_b->Nx);
		delta_y = (front->rect_grid->GU[1] - 
			   front->rect_grid->GL[1])/(rsp_b->Ny);

		screen("Enter the name of the file containing the "
		       "initial interface amplitudes.\n "
		       "\tThe full directory path name of the file" 
		        "is preferred: ");
		(void) Gets(s);
		if((READAMPLITUDES = fopen(s,"r")) == NULL)
		{
		    screen("ERROR: the file %s does not exist in " 
		            "the indicated directory.\n",s); 
		    screen("Please give full directory path of" 
		            "the file.\n");
		    clean_up(ERROR);
		}

		/* It appears that this space need not be allocated for rsp_a 
		   also since it is rsp_b which is an argument to the 
		   make_random_surface() function call made later in this 
		   function, init_random_surface(). 24 Feb. 
		*/
		uni_array(&(rsp_b->x_coord),rsp_b->Nx,FLOAT);
		uni_array(&(rsp_b->y_coord),rsp_b->Ny,FLOAT);
		bi_array(&(rsp_b->amplitudes),(rsp_b->Nx),(rsp_b->Ny),FLOAT);
	    

		/* 
		   TODO: CHECK THAT THE FILE READAMPLITUDES CONTAINS THE 
		   STATED NUMBER OF ENTRIES. - BACKUP ERROR CHECK.
		*/
		for(i = 0; i < rsp_b->Nx; i++)
		    rsp_b->x_coord[i] = L[0] + (i + 0.5)*delta_x;
		for(j = 0; j < rsp_b->Ny; j++)
		    rsp_b->y_coord[j] = L[1] + (j + 0.5)*delta_y;

		for(i = 0; i < rsp_b->Nx; i++)
		{
		    for(j = 0; j < rsp_b->Ny; j++)
		    {
		        if ((fscanf(READAMPLITUDES,"%lf",
				    &(rsp_b->amplitudes[i][j]))) !=1)
			{
			    screen("ERROR: while reading data from file %s\n"
				   "Check if file %s contains %dx%d data "
				    "points.\n",
				   s,rsp_b->Nx,rsp_b->Ny);
			    clean_up(ERROR);
			}
		    }
		}
	    }
	    break;

	case MULTI_MODE:
	default:

	    /* For all modes: */

	    for (i = 0; i < number_modes; ++i) 
	    {
	    	static const char *dname[2] = {"x","y"};
	    	screen("Enter the amplitude of mode %d: ",i);
	    	(void) Scanf("%f\n",&A);
	    	screen("Enter the phase (in degrees) of mode %d: ",i);
	    	(void) Scanf("%f\n",&phase);
	    	phase = radians(phase);
	    	for (j = 0; j < dim-1; ++j)
	    	{
	    	    screen("Enter the number of periods in "
		           "the %s direction for mode %d: ",
			   dname[j],i);
		    (void) Scanf("%f\n",&nu);
		    wave_number = 2.*PI*nu/((U[j]-L[j]));
		    phase += L[j]*wave_number;
		    rsp_a->mode[i].wave_number[j] = wave_number;
		    rsp_b->mode[i].wave_number[j] = wave_number;
		}
		rsp_a->mode[i].phase = phase;
		rsp_b->mode[i].phase = phase;
		rsp_a->mode[i].amplitude = A;
		rsp_b->mode[i].amplitude = A;
	    }
	    screen("\n");
	    break;
	}

	screen("Choices for the perturbation boundary type are\n"
	       "\tPERIODIC (p)\n"
	       "\tSYMMETRIC (s)\n"
	       "\tUNMODIFIED (u)\n");
	for (j = 0; j < dim-1; ++j)
	{
	    screen("Enter the boundary type of perturbation "
		   "in direction %d (dflt = ",j);
	    rsp_b->pert_bdry_type[j] = rsp_a->pert_bdry_type[j];
	    switch (rsp_a->pert_bdry_type[j])
	    {
	    case PERIODIC:
		screen("p");
		break;
	    case SYMMETRIC:
		screen("s");
		break;
	    case UNMODIFIED:
	    default:
		screen("u");
		rsp_a->pert_bdry_type[j] = UNMODIFIED;
		rsp_b->pert_bdry_type[j] = UNMODIFIED;
		break;
	    }
	    screen("): ");
	    (void) Gets(s);
	    if (strncasecmp(s,"p",1) == 0)
	    {
		rsp_a->pert_bdry_type[j] = PERIODIC;
		rsp_b->pert_bdry_type[j] = PERIODIC;
	    }
	    else if (strncasecmp(s,"s",1) == 0)
	    {
		rsp_a->pert_bdry_type[j] = SYMMETRIC;
		rsp_b->pert_bdry_type[j] = SYMMETRIC;
	    }
	    else if (strncasecmp(s,"u",1) == 0)
	    {
		rsp_a->pert_bdry_type[j] = UNMODIFIED;
		rsp_b->pert_bdry_type[j] = UNMODIFIED;
	    }
	    switch (rsp_a->pert_bdry_type[j])
	    {
	    case PERIODIC:
		rect_boundary_type(intfc,j,0) = PERIODIC_BOUNDARY;
		rect_boundary_type(intfc,j,1) = PERIODIC_BOUNDARY;
		break;
	    case SYMMETRIC:
		rect_boundary_type(intfc,j,0) = REFLECTION_BOUNDARY;
		rect_boundary_type(intfc,j,1) = REFLECTION_BOUNDARY;
		break;
	    case UNMODIFIED:
	    default:
		rect_boundary_type(intfc,j,0) = UNKNOWN_BOUNDARY_TYPE;
		rect_boundary_type(intfc,j,1) = UNKNOWN_BOUNDARY_TYPE;
		break;
	    }
	}

	screen("Enter the density below, above: ");
	(void) Scanf("%f %f\n",&rho_b,&rho_a);

	off_tracking = NO;
	screen("Type y to turn off tracking for the contact: ");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	    off_tracking = YES;

	(void) prompt_for_eos_params(init,ip,YES,"");
	if (off_tracking)
	{
	    params = prompt_for_eos_params(init,ip,YES," for the whole domain");
	    params2 = params;
	}
	else
	{
	    params = prompt_for_eos_params(init,ip,YES," for the fluid below");
	    params2 = prompt_for_eos_params(init,ip,YES," for the fluid above");
	}

	screen("Enter the ambient pressure: ");
	(void) Scanf("%f\n",&amb_press);
	for (i = 0; i < dim - 1; ++i)
	    delta_v[i] = 0.0;
	if (init_type == KELHELM)
	    shear = YES;
	else
	{
	    screen("Add velocity shear across interface? (dflt = no): ");
	    (void) Gets(s);
	    shear = (s[0] == 'y' || s[0] == 'Y') ? YES : NO;
	}
	xi = 0.5;
	if (shear == YES)
	{
	    for (i = 0; i < dim - 1; ++i)
	    {
	    	screen("In coord dir %d, enter velocity jump"
	    	       " (above minus below): ",i);
	    	(void) Scanf("%f\n",&delta_v[i]);
	    	rsp_a->delta_v[i] = rsp_b->delta_v[i] = delta_v[i];
	    }
	    screen("Enter the velocity shear weight (default = %g): ",xi);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&xi);
	    }
	}
	surf_ten = prompt_for_surface_tension(CONTACT,"for the contact ");
	rsp_a->surf_ten = rsp_b->surf_ten = surf_ten;

#if defined(TWOD)	
	screen("There are two ways of implementing Neumann boundary\n"
	       "\tconditions.  Half grid offset boundaries (H) or\n"
	       "\treflecting boundary state (F, default).\n"
	       "Enter choice here: ");
	(void) Gets(s);
	if (s[0] == 'h' || s[0] == 'H')
	{
	    offset_N_bdry = YES;
	    front->neumann_bdry_state = NULL;
	    set_obstacle_comp_type(comp_type(COMPOBST),front);
	}
	else
	{
	    offset_N_bdry = NO;
	}
#endif /* defined(TWOD) */

	Init_params(rsp_a->amb_st,params2);
	Dens(rsp_a->amb_st) = rho_a;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            Gas_param *params = Params(rsp_a->amb_st);
            if(params->n_comps != 1)
            {
                    pdens(rsp_a->amb_st)[0] = rho_a;
                    pdens(rsp_a->amb_st)[1] = 0.0;
            }
        }
	Press(rsp_a->amb_st) = amb_press;
	zero_state_velocity(rsp_a->amb_st,dim);
	for (i = 0; i < dim - 1; ++i)
	    Vel(rsp_a->amb_st)[i] = (1.0 - xi)*delta_v[i];
	set_type_of_state(rsp_a->amb_st,TGAS_STATE);
	reset_gamma(rsp_a->amb_st);

	Init_params(rsp_b->amb_st,params);
	Dens(rsp_b->amb_st) = rho_b;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            Gas_param *params = Params(rsp_b->amb_st);
            if(params->n_comps != 1)
            {
                    pdens(rsp_b->amb_st)[0] = 0.0;
                    pdens(rsp_b->amb_st)[1] = rho_b;
            }
        }
	Press(rsp_b->amb_st) = amb_press;
	zero_state_velocity(rsp_b->amb_st,dim);
	for (i = 0; i < dim - 1; ++i)
	    Vel(rsp_b->amb_st)[i] = -xi*delta_v[i];
	set_type_of_state(rsp_b->amb_st,TGAS_STATE);
	reset_gamma(rsp_b->amb_st);

	for (i = 0; i < dim; ++i)
	    coords[i] = 0.5*(L[i]+U[i]);
	g_z = gravity(coords,initial_time(init))[dim-1];
	init_random_surface_perturbation(rsp_b,rsp_a,dim,g_z);

	switch (dim)
	{
#if defined(TWOD)
	case 2:

#if defined(USE_OVERTURE) 
            if(ip->root->overparam->numberOfRefinementLevels != 0 &&
               ip->root->overparam->refinementRatio != 0)
            {
                amr_refinecoeff = (int)pow(ip->root->overparam->refinementRatio,
                                  ip->root->overparam->numberOfRefinementLevels-1); 
                num_points = 2*front->rect_grid->gmax[0]*amr_refinecoeff; 
            }
            else
                num_points = 2*front->rect_grid->gmax[0];  
	    make_random_curve(CONTACT,num_points,L[0],U[0],
			      rsp_b,COMPB,COMPA,surf_ten,off_tracking,dim);
#else /* if defined(USE_OVERTURE) */
	    make_random_curve(CONTACT,2*front->rect_grid->gmax[0],L[0],U[0],
			      rsp_b,COMPB,COMPA,surf_ten,off_tracking,dim);
#endif /* if defined(USE_OVERTURE) */  

	    if (offset_N_bdry)
	    {
	    	RECT_GRID *gr = front->rect_grid;
	    	double y;

	    	exclude_comp(COMPOBST,intfc);
	    	y = U[1] - 0.51*cell_width(gr->gmax[1]-1,1,gr);
	    	make_neumann_curve(intfc,NULL,NULL,L[0],y,U[0],y,
				   YES,YES,COMPOBST,COMPA,&cur);
		y = L[1] + 0.51*cell_width(0,1,gr);
		make_neumann_curve(intfc,NULL,NULL,U[0],y,L[0],y,
				   YES,YES,COMPOBST,COMPB,&cur);
	    }
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:	
	    make_random_surface(CONTACT,front,rsp_b,COMPB,COMPA,surf_ten,dim);
	    /* Following added 22 March 2004   
	       if (rsp_b->read_amplitudes_from_file == 0)
	       {
	       free(&(rsp_b->x_coord));
	       free(&(rsp_b->y_coord));
	       free(&(rsp_b->amplitudes));
	       }
	    */
	    break;
#endif /* defined(THREED) */
	}	
}		/*end init_random_surface*/


/*
*			random_surface_perturbed_state():
*
*	Computes the fluid state s, belonging to component comp, at the
*	requested position coords at time t, according to the linear theory
*	for compressible Rayleigh-Taylor instability.  The state is given as
*	a sum of a 0th order equilibrium state and a 1st order perturbation,
*	where the perturbation is a superposition of the contributions due to
*	the independent Fourier modes of the perturbed material interface,
*	as computed in 	RT_single_mode_perturbation_state().  The equilibrium
*	state is assumed to be isothermal and hydrostatic.
*	
*	Prior to the sum of the perturbations (but after the computation of the
*	equilibrium state), the longitudinal (z) coordinate of the requested
*	position coords is offset by an amount equal to the displacement of
*	the interface position at the same transverse position.  This offset
*	produces a 2nd order change in the perturbation of the fluid state,
*	and is therefore permissible within the degree of approximation of the
*	linear theory.  The purpose of this offset is to reduce the size of
*	the difference between the left and right state pressures at points on
*	the material interface, which is ideally zero but has a truncation
*	error whose magnitude is 2nd order in the perturbation amplitude.
*/

LOCAL void random_surface_perturbed_state(
	Locstate	s,
	double		*coords,
	double		t,
	COMPONENT	comp,
	int		stype)
{
	SINE_PERT	*RS = (SINE_PERT*) comp_type(comp)->extra;
	Locstate	amb_st;
	int		i, m, dim, number_modes;
	double		g_z,z_intfc_pert;
	static Locstate tmp_st = NULL;
	double		crds[MAXD];

	number_modes = RS->number_modes;

	amb_st = RS->amb_st;
	dim = Params(amb_st)->dim;

	if (tmp_st == NULL)
	{
	    (*Params(amb_st)->_alloc_state)(&tmp_st,Params(amb_st)->sizest);
	    Set_params(tmp_st,amb_st);
	}

	z_intfc_pert = pert_interface(RS,coords,0.0,dim);
	g_z = gravity(coords,t)[dim - 1];
	for (i = 0; i < dim; ++i)
	    crds[i] = coords[i];
	crds[dim-1] = coords[dim-1] - (z_intfc_pert - RS->z_intfc);

	isothermal_stratified_state(s,coords[dim-1]-RS->z_intfc,g_z,amb_st);
	set_state(s,TGAS_STATE,s);
	for (i = 0; i < dim; ++i) Vel(s)[i] = Vel(amb_st)[i];
	if (RS->init_type == RT)
	{
	    for (m = 0; m < number_modes; ++m)
	    {
	    	RT_single_mode_perturbation_state(tmp_st,crds,t,amb_st,
						  RS->z_intfc,RS->z_bdry,
						  RS->mode+m,g_z);
	    	Dens(s) += Dens(tmp_st);
	    	Press(s) += Press(tmp_st);
	    	for (i = 0; i < dim; ++i)
		    Vel(s)[i] += vel(i,tmp_st);
	    }
	}
	set_state(s,stype,s);
}		/*end random_surface_perturbed_state*/

LOCAL	double	RT_sigma(
	SINE_PERT	*RS_b,
	SINE_PERT	*RS_a,
	int		i_mode,
	int		dim,
	double		g_z)
{
	boolean		   deb_RT_RS = (debugging("RT_RS")) ? YES : NO;
	double		   *k, k_sqr;
	double		   rho_a, rho_b;
	double		   sigma;
	int		   i;
	static const double epsilon = 0.000001; /*TOLERANCE*/

	double		   f[32+1], s[32+1]; /*Note the relation between the */
	static const int   max_steps = 32;   /*size of f and s and max_steps */

	debug_print("RT_RS","Entered RT_sigma()\n");

	if (deb_RT_RS == YES)
	{
	    (void) printf("Mode number %3d\n",i_mode);
	    (void) printf(" step\t\tsigma\t\t\tf\t\t\ts\n");
	}

	rho_b = Dens(RS_b->amb_st);
	rho_a = Dens(RS_a->amb_st);
	k = RS_b->mode[i_mode].wave_number;
	k_sqr = scalar_product(k,k,dim-1);

	i = 0;
	s[i] = 0.25*g_z*sqrt(k_sqr)*(rho_b - rho_a)/(rho_a + rho_b);
	sigma = sqrt(s[i]);
	f[i] = RT_RS_fab(s[i],RS_b,RS_a,i_mode,dim,g_z);
	if (deb_RT_RS == YES)
	    (void) printf("%5d%21g%24g%24g\n",i,sigma,f[i],s[i]);

	i = 1;
	s[i] = 4.0*g_z*sqrt(k_sqr)*(rho_b - rho_a)/(rho_a + rho_b);
	sigma = sqrt(s[i]);
	f[i] = RT_RS_fab(s[i],RS_b,RS_a,i_mode,dim,g_z);
	if (deb_RT_RS == YES)
	    (void) printf("%5d%21g%24g%24g\n",i,sigma,f[i],s[i]);

	while (fabs(f[i]-f[i-1]) > epsilon)
	{
	    if (i == max_steps)
	    {
		screen("ERROR in RT_sigma(), Secant method failed\n");
		clean_up(ERROR);
		break;
	    }

	    s[i+1] = -f[i-1]*(s[i]-s[i-1])/(f[i]-f[i-1])+s[i-1];
	    sigma = sqrt(s[i+1]);
	    f[i+1] = RT_RS_fab(s[i+1],RS_b,RS_a,i_mode,dim,g_z);
	    if (deb_RT_RS == YES)
	    {
	    	(void) printf("%5d%21g%24g%24g\n",
			      i+1,sigma,f[i+1],s[i+1]);
	    }
	    ++i;
	}
	debug_print("RT_RS","Left RT_sigma()\n");
	return sigma;
}		/*end RT_sigma*/

LOCAL double RT_RS_fab(
	double		s,
	SINE_PERT	*RS_b,
	SINE_PERT	*RS_a,
	int		i_mode,
	int		dim,
	double		g_z)
{
	double		f_a, f_b;
 	double		*k, k_sqr;

	k = RS_b->mode[i_mode].wave_number;
	k_sqr = scalar_product(k,k,dim-1);

	f_a = RT_RS_f(s,RS_a->amb_st,RS_a->z_intfc - RS_a->z_bdry,k_sqr,g_z);

	f_b = RT_RS_f(s,RS_b->amb_st,RS_b->z_intfc - RS_b->z_bdry,k_sqr,g_z);

	return f_b - f_a;
}		/*end RT_RS_fab*/

/*ARGSUSED*/
EXPORT	void	lin_pert_state(
	Locstate	state,
	double		*coords,
	Front		*fr,
	Wave		*wave,
	int		stype)
{
	COMPONENT	comp;
	SINE_PERT	*pert = Sine_pert(comp_type(COMPA));
	double		z;
	int		dim = fr->rect_grid->dim;

	z = pert_interface(pert,coords,fr->time,dim);
	comp = (z < coords[dim-1]) ? COMPA : COMPB;
	random_surface_perturbed_state(state,coords,fr->time,comp,stype);
}		/*end lin_pert_state*/

#if defined(TWOD)
/*
*		make_random_curve():
*
*	Make_random_surface_curve constructs the initial contact wave
*	for the "random surface" (statistics of fingers) problem.
*
*	TODO: merge this function with make_layer_surf() below.
*/

/*ARGSUSED*/
LOCAL void make_random_curve(
	int		w_type,
	int		num_points,
	double		xlower,
	double		xupper,
	SINE_PERT	*RS_b,
	COMPONENT	compb,
	COMPONENT	compa,
	double		surf_ten,
	boolean            off_tracking,
	int		dim)
{
	SINE_PERT	*pert = Sine_pert(comp_type(compb));
	int		i;
	double		dx,coords[MAXD];
	CURVE		*cur;
	NODE		*ns,*ne;

	coords[0] = xupper;
	coords[dim-1] = pert_interface(pert,coords,0.0,dim);
	ns = make_node(Point(coords));
	coords[0] = xlower;
	coords[dim-1] = pert_interface(pert,coords,0.0,dim);
	ne = make_node(Point(coords));
	cur = make_curve(compb,compa,ns,ne);
	dx = (xupper - xlower)/num_points;
	for (i = 1; i < num_points; ++i)
	{
	    coords[0] =  xlower + i*dx;
	    coords[dim-1] = pert_interface(pert,coords,0.0,dim);
	    if (insert_point_in_bond(Point(coords),cur->first,cur) !=
	        FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in make_random_curve(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	}

	wave_type(cur) = w_type;
	start_status(cur) = end_status(cur) = INCIDENT;
	surface_tension(cur) = surf_ten;

	if (off_tracking)
            untracked_hyper_surf(cur) = YES;
}		/*end make_random_curve*/
#endif /* defined(TWOD) */

LOCAL	void set_pert_modes(
	SINE_PERT	*rsp_a,
	SINE_PERT	*rsp_b,
	int		dim,
	int		number_modes,
	double		**F_wave_number,
	double		*F_amplitude,
	double		*F_phase)
{
	int		i, j;

	for (i = 0; i < number_modes; ++i) 
	{
	    for (j = 0; j < dim-1; ++j)
	    {
	    	rsp_a->mode[i].wave_number[j] = F_wave_number[i][j];
	    	rsp_b->mode[i].wave_number[j] = F_wave_number[i][j];
	    }
	    rsp_a->mode[i].amplitude = F_amplitude[i];
	    rsp_b->mode[i].amplitude = F_amplitude[i];
	    rsp_a->mode[i].phase = F_phase[i];
	    rsp_b->mode[i].phase = F_phase[i];
	}
}		/*end set_pert_modes*/



/*
*		init_random_surface_perturbation():
*/

LOCAL void init_random_surface_perturbation(
	SINE_PERT	*RS_b,
	SINE_PERT	*RS_a,
	int		dim,
	double		g_z)
{
	int		number_modes = RS_b->number_modes;
	int		i;

	debug_print("RS","Entered init_random_surface_perturbation()\n");
	if (RS_b->init_type == KELHELM)
	{
	    init_KH_sine_perturbation(RS_b,RS_a,dim,g_z);
	    debug_print("RS","Left init_random_surface_perturbation()\n");
	    return;
	}
	if (RS_b->init_type == RT) 
	{
	    /* Calculate sigma[i]. */
	    for (i = 0; i < number_modes; ++i)
	    	RS_b->mode[i].growth_rate = RS_a->mode[i].growth_rate =
				RT_sigma(RS_b,RS_a,i,dim,g_z);
	}
	debug_print("RS","Left init_random_surface_perturbation()\n");
}		/*end init_random_surface_perturbation*/


EXPORT	void init_kelvin_helmholtz(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	char		surf_type[5];
	double		mean_height;
	double		A, lambda[MAXD], wave_number[MAXD], phase;
	double		*L = front->rect_grid->GL, *U = front->rect_grid->GU;
	double		coords[3];
	double		rho_b, rho_a;
	double		amb_press, delta_v[MAXD-1];
	double		surf_ten;
	Gas_param	*params, *params2;
	SINE_PERT	*khsp_b,*khsp_a;
	int		number_modes, i, ms;
	int		j, dim = front->interf->dim;
	double		**F_wave_number, *F_amplitude, *F_phase;
	const double	*grav;
	double		g_z;

	screen("Enter the choice of initial front description\n"
	       "Two types of descriptions are supported\n"
	       "       Multiple mode description (MM),\n"
	       "       Multiple sine section description (MS).\n"
	       "       Enter choice: ");
	(void) Scanf("%s\n",surf_type);
	switch (surf_type[1]) 
	{
	case 'M':
	case 'm':
	    ms = 0;
	    break;
	case 'S':
	case 's':
	    ms = 1;
	    break;
	default: 
	    ms = -1;
	    screen("ERROR: unrecognized choice for initial "
	           "surface description.\n");
	    clean_up(ERROR);
	    break;
	}

	if (ms == 0) 
	{
	    screen("Enter the number of modes in the initial ");
	    screen("interface: ");
	    (void) Scanf("%d\n",&number_modes);
	}
	else 
	    number_modes = MAX_NUMBER_MODES;

	set_kh_sine_perturbed_comp_type(comp_type(COMPA),front,number_modes,
		                        KELHELM);
	khsp_a = (SINE_PERT*)comp_type(COMPA)->extra;
	set_kh_sine_perturbed_comp_type(comp_type(COMPB),front,number_modes,
		                        KELHELM);
	khsp_b = (SINE_PERT*)comp_type(COMPB)->extra;

	/* For the interface as a whole: */

	if (ms) 
	{
	    bi_array(&F_wave_number,number_modes,dim-1,FLOAT);
	    uni_array(&F_phase,number_modes,FLOAT);
	    uni_array(&F_amplitude,number_modes,FLOAT);

	    F_bubble(number_modes,F_wave_number,F_amplitude,F_phase,
		     U[0]-L[0],L[0]);
	    set_pert_modes(khsp_a,khsp_b,dim,number_modes,
			   F_wave_number,F_amplitude,F_phase);

	    free_these(3,F_wave_number,F_amplitude,F_phase);
	}
	else 
	{
	    for (i = 0; i < number_modes; ++i) 
	    {
	    	screen("Enter the amplitude of mode %d: ",i);
	    	(void) Scanf("%f\n",&A);
	    	khsp_a->mode[i].amplitude = A;
	    	khsp_b->mode[i].amplitude = A;
	    	screen("Enter the phase of mode %d: ",i);
	    	(void) Scanf("%f\n",&phase);
	    	phase = radians(phase);
	    	for (j = 0; j < dim-1; ++j)
	    	{
		    screen("For the coord direction %d ",j);
		    screen("Enter the wavelength of mode %d: ",i);
		    (void) Scanf("%f\n",&lambda[j]);
		    wave_number[j] = 2.*PI/(lambda[j]*(U[j]-L[j]));
		    khsp_a->mode[i].wave_number[j] = wave_number[j];
		    khsp_b->mode[i].wave_number[j] = wave_number[j];
		    phase += L[j]*wave_number[j];
		    khsp_a->mode[i].phase = phase;
		    khsp_b->mode[i].phase = phase;
		}
	    }
	}


	screen("Enter the mean height of the front above L[%d]: ",dim-1);
	(void) Scanf("%f\n",&mean_height);

	mean_height += L[dim-1];
	khsp_b->z_intfc = mean_height;
	khsp_b->z_bdry = L[dim-1];
	khsp_a->z_intfc = mean_height;
	khsp_a->z_bdry = U[dim-1];

	(void) prompt_for_eos_params(init,ip,YES,"");
	params = prompt_for_eos_params(init,ip,YES," for the fluid below");
	params2 = prompt_for_eos_params(init,ip,YES," for the fluid above");

	screen("Enter the density below, above: ");
	(void) Scanf("%f %f\n",&rho_b,&rho_a);
	screen("Enter the ambient pressure\n");
	screen("\tand velocity jump (above minus below): ");
	/* TODO threed: vel is a vector */
	(void) Scanf("%f %f\n",&amb_press,&delta_v[0]);
	for (i = 0; i < dim-1; ++i)
	{
	    khsp_b->delta_v[i] = khsp_a->delta_v[i] = delta_v[i];
	}
	surf_ten = prompt_for_surface_tension(CONTACT,"for the slip surface ");
	khsp_a->surf_ten = khsp_b->surf_ten = surf_ten;


	Init_params(khsp_a->amb_st,params2);
	Dens(khsp_a->amb_st) = rho_a;
	Press(khsp_a->amb_st) = amb_press;
	zero_state_velocity(khsp_a->amb_st,dim);
	set_type_of_state(khsp_a->amb_st,TGAS_STATE);
	reset_gamma(khsp_a->amb_st);

	Init_params(khsp_b->amb_st,params);
	Dens(khsp_b->amb_st) = rho_b;
	Press(khsp_b->amb_st) = amb_press;
	zero_state_velocity(khsp_b->amb_st,dim);
	set_type_of_state(khsp_b->amb_st,TGAS_STATE);
	reset_gamma(khsp_b->amb_st);

	for (i = 0; i < dim; ++i)
	    coords[i] = 0.5*(L[i]+U[i]);
	grav = gravity(coords,initial_time(init));
	g_z = grav[dim-1];
	init_KH_sine_perturbation(khsp_b,khsp_a,dim,g_z);

	switch (dim)
	{
#if defined(TWOD)
	case 2:
	    make_random_curve(CONTACT,2*front->rect_grid->gmax[0],L[0],U[0],
			      khsp_b,COMPB,COMPA,surf_ten,NO,dim);
		break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3: 
	    make_random_surface(CONTACT,front,khsp_b,
				COMPB,COMPA,surf_ten,dim);
	    break;
#endif /* defined(THREED) */
	}
}		/*end init_kelvin_helmholtz*/



/*
*			init_KH_sine_perturbation():
*
*	Initializes information needed by get_KH_sine_perturbed_state() and
*	print_sine_perturbed_front().  init_KH_sine_perturbation() calls the 
*	FORTRAN routine speed() to calculate the growth_rate and 
*	propagation_speed for the perturbed interface.
*/

LOCAL void init_KH_sine_perturbation(
	SINE_PERT	*sp_b,
	SINE_PERT	*sp_a,
	int		dim,
	double		g_z)
{
	double		csqr_a,csqr_b,vel_b[MAXD-1],vel_a[MAXD-1],dist_a,dist_b;
	double		g_rate,inc_g_rate,p_speed;
	double		rho_a, rho_b;
	int		i,number_modes = sp_b->number_modes;

	/*
	*	We work in the frame in which the propagation speed is zero.
	*	First vel_a and vel_b are set so this is TRUE in the
	*	incompressible limit.  Then the calculated propagation speed
	*	is subtracted from vel_a and vel_b to give the stream
	*	velocities above and below.
	*/

	rho_a = Dens(sp_a->amb_st);
	rho_b = Dens(sp_b->amb_st);
	dist_b = sp_b->z_intfc - sp_b->z_bdry;
	dist_a = sp_a->z_bdry - sp_a->z_intfc;
	for (i = 0; i < dim-1; ++i)
	{
		vel_b[i] = - rho_a * sp_b->delta_v[i] / (rho_b + rho_a);
		vel_a[i] =   rho_b * sp_a->delta_v[i] / (rho_b + rho_a);
	}
	csqr_b = sound_speed_squared(sp_b->amb_st);
	csqr_a = sound_speed_squared(sp_a->amb_st);

	for (i = 0; i < number_modes; ++i)
	{
		/*TODO 3D wavenumber */
		FORTRAN_NAME(speed)(&g_rate,&inc_g_rate,&p_speed,
			       &sp_b->mode[i].wave_number[0],
			       &rho_a,&rho_b,&csqr_a,&csqr_b,&dist_a,&dist_b,
			       &vel_a[0],&vel_b[0],&g_z,&sp_b->surf_ten);

		sp_b->mode[i].growth_rate = g_rate;
		sp_b->mode[i].propagation_speed = 0.;
	
		sp_a->mode[i].growth_rate = g_rate;
		sp_a->mode[i].propagation_speed = 0.;
	}
	for (i = 0; i < dim-1; ++i)
	{
		sp_a->stream_velocity[i] = vel_a[i];
		sp_b->stream_velocity[i] = vel_b[i];
	}
}		/*end init_KH_sine_perturbation*/

/*
*			get_KH_sine_perturbed_state():
*
*	Interfaces to the FORTRAN subroutine khstate() that
*	calculates a gas state obtained by perturbation
*	analysis of the Kelvin-Helmholtz instabilty.
*
*	See I. G. Currie, Fundamental Mechanics of Fluids, Chapter 6.,
*	or see Lamb's Hydrodynamics for the incompressible analysis;
*	for the compressible analysis, use Crocco's equation in place
*	of Bernoulli's equation and the wave equation in place of
*	Laplace's equation.
*/


/*ARGSUSED*/
LOCAL void get_KH_sine_perturbed_state(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*comp_type,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	SINE_PERT	*sp = (SINE_PERT*) comp_type->extra;
	Gas_param	*params = Params(sp->amb_st);
	int		dim;
	int		i,j;
	static Locstate tmps = NULL;

	if (tmps == NULL)
	{
	    (*Params(sp->amb_st)->_alloc_state)(&tmps,params->sizest);
	    Init_params(tmps,params);
	}

	dim = params->dim;

	set_state(s,TGAS_STATE,sp->amb_st);
	for (j = 0; j < dim-1; ++j)
		Vel(s)[j]   = sp->stream_velocity[j];
	Vel(s)[dim-1] = 0.0;

	for (i = 0; i < sp->number_modes; ++i)
	{
		KH_single_mode_state(tmps,coords,0.0,sp->amb_st,
			sp->stream_velocity[0],sp->z_intfc,
			sp->z_bdry,sp->mode+i);

		Dens(s) += Dens(tmps);
		Press(s) += Press(tmps);
		for (j = 0; j < dim; ++j)
			Vel(s)[j] += Vel(tmps)[j];
	}
	/*
	*  NOTE: The orignal code would have set
	*  Energy(s) = kinetic_energy(s) + Dens(s)*specific_enthalpy(s)
	*/
	set_state(s,stype,s);

#if defined(DEBUG_INIT_PERT)
	if (debugging("init"))
	{
	    for (j = 0; j < dim; ++j)
	    	(void) printf("coords[%d] = %g ",j,coords[j]);
	    (void) printf("t = %g\n",t);
	    (void) printf("z_intfc = %g z_bdry = %g\n",
	    	      sp->z_intfc,sp->z_bdry);
	    (void) printf("density = %g ambient pressure = %g\n",
	    	Dens(sp->amb_st),pressure(sp->amb_st));
	    for (j = 0; j < dim-1; ++j)
	    	(void) printf("stream vel[%d] = %g ",
	    		      j,sp->stream_velocity[j]);
	    (void) printf("\n");
	    g_print_state(s);
	    for (i = 0; i < sp->number_modes; ++i)
	    {
	        (void) printf("mode %d\n",i);
	        for (j = 0; j < dim-1; ++j)
	           (void) printf("wave_number[%d] = %g ",
				 j,sp->mode[i].wave_number[j]);
		(void) printf("\n");
		(void) printf("amplitude = %g ",sp->mode[i].amplitude);
		(void) printf("growth_rate = %g ",sp->mode[i].growth_rate);
		(void) printf("phase = %g\n",sp->mode[i].phase);
	    }
	}
#endif /* defined(DEBUG_INIT_PERT) */
}		/*end get_KH_sine_perturbed_state*/

/*
*			F_bubble():
*
*	This function calculates the Fourier coefficients
*	for the superposition of the individual bubbles
*	for the multiple bubble compressible Rayleigh-Taylor
*	simulation.   In this context,  a single "bubble"
*	with center c,  diameter d and amplitude A.  Is
*	the Fourier polynomial of degree nmode obtained
*	by truncating the Fourier series for the function
*
*		  _
*		 /  
*		 |  - 0.5 * A * (1 + cos(2*PI*(x-c)/d))  for |x - c| < d/2
*		/
*	f(x) = <
*		\
*		 |  0 					 otherwise
*		  \
*		   -
*
*	Each "bubble" may be modified by an additive constant
*	to adjust its mean position to the desired value.
*
*	TODO THREED
*/

LOCAL void F_bubble(
	int		nmode,
	double		**F_wave_number,
	double		*F_amplitude,
	double		*F_phase,
	double		X,
	double		XL)
{
	int		i, j, nb;
	double		*x, *a, *l;
	double		A, B, amp, c;

	screen("Enter the number of bubbles in the initial interface: ");
	(void) Scanf("%d\n",&nb);

	uni_array(&x,nb,FLOAT); uni_array(&a,nb,FLOAT); uni_array(&l,nb,FLOAT);
	for (i = 0; i < nb; ++i)
	{
		screen("Enter the center of bubble %d: ",i+1);
		(void) Scanf("%f\n",&x[i]);
		screen("Enter the amplitude of bubble %d: ",i+1);
		(void) Scanf("%f\n",&a[i]);
		screen("Enter the diameter of bubble %d: ",i+1);
		(void) Scanf("%f\n",&l[i]);
	}

	for (i = 1; i <= nmode; ++i) 
	{
		A = B = 0.0;
		for (j = 0; j < nb; ++j)
		{
			c = i * l[j] / X;
			if (fabs(c - 1.0) < EPS)
			{
				amp = -0.5*a[j]*l[j]/X;
			}
			else
				amp = a[j]*l[j]*sin(PI*c)/(PI*X*c*(sqr(c)-1.0));
			A += amp*cos(2.0*PI*i*(x[j] -XL)/X);
			B += amp*sin(2.0*PI*i*(x[j] -XL)/X);
		}
		F_amplitude[i-1]   = hypot(A,B);
		F_wave_number[i-1][0] = 2.0*PI*i/X;
		F_phase[i-1]       = 2.0*PI*i*XL/X + atan2(-A,B);
	}
	free_these(3,x,a,l);
}		/*end F_bubble*/


#define	Ct(i)	comp_type(layer[i]->comp)
#define	Rtp(i)	Rt_perturbed(Ct(i))

EXPORT	void	rt_ml_linear_pert(
	LAYER_SYS	*layer_sys,
	INIT_DATA	*init)
{
	Front		*front = layer_sys->front;
	int		i, j, k, nstep, num_modes, unstable_modes = 0;
	int	        num_layers = layer_sys->num_layers; 
	int		dim = front->rect_grid->dim; 
	COMP_TYPE_TYPE	prob_type = RT_PERTURBED;
	int		flag_unstable = YES;
	double		*velo;
	double		*amp, *phase, **wv_num, sig;
	FOURIER_POLY	*fpoly0,*fpoly;
	NORMAL_MODE	**normal_mode;
	LAYER		**layer = layer_sys->layer;
	LIN_PERT	*lin_pert;
	_RT_PERTURBED 	*rtp;

 	scalar(&lin_pert, sizeof(LIN_PERT));
	lin_pert->layer_sys = layer_sys;

	/* check validity, find problem type etc. */
	if (num_layers == 1)
	{
	    screen("\nERROR: in rt_ml_linear_pert().\n");
	    screen("The code assumes the perturbation is induced by\n");
	    screen("disturbances on interfaces, and there is none here.\n");
	    clean_up(ERROR);
	}

	for (i = 1; i <= num_layers; ++i)
	{
	    velo = Vel(Rt_kh(Ct(i))->ref_state);
	    if (velo[dim-1] != 0.0)
	    {
	    	screen("ERROR in rt_ml_linear_pert(), "
	    	       "v_z is nozero in the %d%s layer,"
	    	       " no linear perturbation theory applies.\n",
		       i,ordinal_suffix(i));
	    	clean_up(ERROR);
	    }
	    for (j = 0; j <= dim-2; ++j)
	        if (velo[j] != 0.0)
	    	    prob_type = KH_PERTURBED;
	}

	/* set up normal modes */
	fpoly0 = layer[1]->upper_surf->fpoly;
	num_modes = fpoly0->num_modes;
	uni_array(&normal_mode,num_modes,sizeof(NORMAL_MODE *));
	for (i = 0; i < num_modes; ++i)
	{
	    scalar(&normal_mode[i],sizeof(NORMAL_MODE));
	    normal_mode[i]->amp_max = fpoly0->A[i];
	    normal_mode[i]->phase = fpoly0->phase[i];
	    normal_mode[i]->wv_num = fpoly0->nu[i];
	    normal_mode[i]->ksq = 0.0;
	    for (k = 0; k <= dim-2; ++k)
	    	normal_mode[i]->ksq+=fpoly0->nu[i][k]*(fpoly0->nu[i][k]);
	}
	layer[1]->upper_surf->fpoly = NULL;

	/* changing comp_type */
	for (i = 1; i <= num_layers; ++i)
	{
	    if (prob_type == RT_PERTURBED)
	    {
	    	set_rt_perturbed_comp_type(Ct(i),front);
	    	rtp = Rtp(i);
	    	rtp->layer_label = layer[i]->layer_label;
	    	rtp->lower_surf = layer[i]->lower_surf;
	    	rtp->upper_surf = layer[i]->upper_surf;
	    	rtp->num_modes = -1;
	    	rtp->normal_mode = normal_mode;
	    }
	    else if (prob_type == KH_PERTURBED)
	    	set_kh_perturbed_comp_type(Ct(i),front);
	}

	/* solving for a single normal mode */
	while ((i = get_next_mode(num_modes, lin_pert)) < num_modes)
	{
	    if (i >= 0)
	    {
	    	lin_pert->normal_mode = normal_mode[i];
	    	switch (prob_type)
	    	{
	    	case RT_PERTURBED:
	    	    flag_unstable = (rt_single_mode(lin_pert,init) == YES) ?
				    YES : NO;
		    break;
		case KH_PERTURBED:
		    flag_unstable = (kh_single_mode(lin_pert) == YES) ?
				    YES : NO;
		    break;
		default:
		    screen("\nERROR: in rt_ml_linear_pert().\n");
		    screen("Unrecoganized problem.\n");
		}
	    }
	    unstable_modes = 
		   list_modes(normal_mode,lin_pert,num_modes,flag_unstable);
	}

	brdcst_info(lin_pert,normal_mode,num_modes,&unstable_modes);

	num_modes = unstable_modes;
	if (num_modes == 0)
	{
	    screen("\nERROR: in rt_ml_linear_pert().\n");
	    screen("All modes are stable!\n");
	    screen("There is no reason to apply ");
	    screen("linear perturbation theory in this case.\n");
	    clean_up(ERROR);
	}
	if (debugging("lin_pert"))
	{
	    (void) printf("The total number of unstable modes is %d\n",
	    	          num_modes);
	    (void) printf("The total number of interpolating points is %d.\n",
		          lin_pert->tot_pts);
	    for (i = 1; i <= num_layers; ++i)
	    {
	    	nstep = Rtp(i)->lin_pert_intvl;
	    	(void) printf("In the %d%s layer, ",i,ordinal_suffix(i));
	    	(void) printf("there are %d interpolating points.\n",nstep+1);
	    }
	    for (j = 0; j < num_modes; ++j)
	    {
	    	(void) printf("\nFor the %d%s mode :\n",j,ordinal_suffix(j));
	    	(void) printf("The wave number |k| is %g.\n",
	    		      sqrt(normal_mode[j]->ksq));
	    	sig = normal_mode[j]->sigma_r;
	    	(void) printf("The growth rate sigma is %g.\n",sig);
	    	for (i = 1; i < num_layers; ++i)
	    	{
	    	    nstep = Rtp(i)->lin_pert_intvl;
	    	    (void) printf("The amplitude of this mode on the ");
	    	    (void) printf("%d%s interface is %g.\n", 
	    			  i,ordinal_suffix(i),
	    			  normal_mode[j]->a[i][nstep]/sig);
	    	}
	    }
	}

	for (i = 1; i <= num_layers; ++i)
	    Rtp(i)->num_modes = num_modes;

	/* modify interfaces */
	for (i = 1; i < num_layers; ++i)
	{
	    fpoly = allocate_fourier_poly(num_modes,dim,NULL);
	    amp = fpoly->A;
	    phase = fpoly->phase;
	    wv_num = fpoly->nu;
	    nstep = Rtp(i)->lin_pert_intvl;
	    for (j = 0; j < num_modes; ++j)
	    {
	    	sig = normal_mode[j]->sigma_r;
	    	amp[j] = normal_mode[j]->a[i][nstep]/sig;
	    	phase[j] = normal_mode[j]->phase;
	    	for (k = 0; k <= dim-2; ++k)
	    	    wv_num[j][k] = normal_mode[j]->wv_num[k];
	    }
	    fpoly->z0 = layer[i]->upper_surf->pbar[dim-1];
	    layer[i]->upper_surf->fpoly = fpoly;
	}
	free(lin_pert);
	free(fpoly0);
}		/*end rt_ml_linear_pert*/


LOCAL	int	get_next_mode(
	int		n_modes,
	LIN_PERT	*lin_pert)
{
	static	int	ith_mode = -1;

	return (lin_pert->layer_sys->front->pp_grid->nn == 1) ?
		++ith_mode : pp_get_next_mode(&ith_mode,n_modes,lin_pert);
}		/*end get_next_mode*/

		

LOCAL	int	list_modes(
	NORMAL_MODE	**normal_mode,
	LIN_PERT	*lin_pert,
	int		n_modes,
	int		flag)
{
	static	int	unstable_modes = 0;

	if (lin_pert->layer_sys->front->pp_grid->nn > 1)
	{
		unstable_modes = pp_list_modes(normal_mode,lin_pert,
					n_modes,flag,unstable_modes);
		return unstable_modes;
	}
	if (flag == NO)
	{
		free(lin_pert->normal_mode);
		lin_pert->normal_mode = NULL;
		return	unstable_modes;
	}
	normal_mode[unstable_modes] = lin_pert->normal_mode; 
	return	++unstable_modes;
}		/*end list_modes*/




LOCAL	int	rt_single_mode(
	LIN_PERT	*lin_pert,
	INIT_DATA	*init)
{
	static const double XACC = 1.0e-12;		/*TOLERANCE*/
	static const double FACC = 1.0e-12;		/*TOLERANCE*/
	int 		i, j, k, nstep, ne, m, n_intvl, max_intvl; 
	int		*flag_root, *il, i_infc, num_unstable;
	int	        nl = lin_pert->layer_sys->num_layers;
	RECT_GRID 	*rect_grid = lin_pert->layer_sys->front->rect_grid;
	double		*L = rect_grid->GL, *U = rect_grid->GU;
	int 		dim = rect_grid->dim;
	const double	*grav;
	double		g_z;
	double		wv_num = sqrt(lin_pert->normal_mode->ksq);
	double		rho_m, rho_p, surf_ten, coords[MAXD], *r_c;
	double		tmp, ub, lb, factor;
	double		***yyy, ***yyy_tmp;
	double		lam_min, lam_max, dlam, lambda, **lam, mis_old, mis_new;
	double		max_a, *a, *b, **aa, **bb;
	Locstate	s;
	Front		*front = lin_pert->layer_sys->front;
	LAYER		**layer = lin_pert->layer_sys->layer;

	for (i = 0; i < dim; ++i)
	    coords[i] = 0.5*(L[i]+U[i]);
	grav = gravity(coords,initial_time(init));
	g_z = grav[dim-1];
	lin_pert->init = init;
	lin_pert->g_z = g_z;
	uni_array(&il,nl,INT);
	/* find how many unstable modes */
	/* 
	 * We have assumed that the only cause for an instability 
	 * is the density inversion across an interface.
	 */
	alloc_state(front->interf,&s,front->sizest);
	num_unstable = 0;
	for (i = 1; i < nl; ++i)
	{
	    r_c = Rtp(i)->rt_kh->ref_coords;
	    for (j = 0; j <= dim-2; ++j) 	coords[j] = r_c[j];
	    coords[dim-1] = layer[i]->upper_surf->pbar[dim-1];
	    Get_tgas_state(coords,s,Ct(i),front->interf,init);
	    rho_m = Dens(s);
	    Get_tgas_state(coords,s,Ct(i+1),front->interf,init);
	    rho_p = Dens(s);
	    surf_ten = layer[i]->upper_surf->surf_ten;
	    if (((rho_p-rho_m)*g_z+surf_ten*wv_num*wv_num) < 0.0)
	    {
	    	++num_unstable;
	    	il[num_unstable] = i;
	    }
	}
	free(s);
	if (num_unstable == 0)
	    return 	NO;

	/* find the eigenvalue lambda */
	ne = 2; 				       /* num of eqns */
	factor = 1.5;				       /* for discretization */
	max_intvl = num_unstable == 1 ? 4 : 2048;      /* for root searching */
	tmp = fabs(g_z*wv_num);
	uni_array(&flag_root,num_unstable+1,INT);
	for (i = 1; i <= num_unstable; ++i)	flag_root[i] = 1;
	bi_array(&lam,num_unstable+1,3,FLOAT);
	for (j = 1; j <= 2; ++j)	/* do it twice for better accuracy */
	{
	    m = 0;
	    for (i = 1; i <= nl; ++i)	
	    {
		ub = layer[i]->upper_surf->pbar[dim-1];
		lb = layer[i]->lower_surf->pbar[dim-1];
		nstep = (int) (factor*((ub-lb)/rect_grid->h[dim-1]+1));
		Rtp(i)->lin_pert_intvl = nstep;
		m += nstep+1;
	    }
	    lin_pert->tot_pts = m;
	    tri_array(&yyy_tmp,num_unstable,ne,m,DOUBLE);
	    yyy = yyy_tmp-1;
	    for (i = 1; i <= num_unstable; ++i) /* offset zero */
	    {
		--yyy[i];
		for (k = 1; k <= ne; ++k)
		    --yyy[i][k];
	    }

	    for (k = 1; k <= num_unstable; ++k) /* loop over unstable intfcs */
	    {
		lin_pert->intfc_label = il[k];
		lin_pert->pointer = (POINTER)yyy[k];

		n_intvl = num_unstable == 1 ? 1 : 8;
		while ((n_intvl *= 4) <= max_intvl)
		{
		    lam_min = tmp*XACC; 
		    lam_max = tmp*(1.0+XACC);
		    dlam = (lam_max-lam_min)/n_intvl;
		    (void) mis_match(lam_max,&mis_old,(POINTER)lin_pert);
		    for (i = 1; i <= n_intvl; ++i)
		    {
			if (mis_old == 0.0)
			{
			    lam[k][j] = lam_max;
			    break;
		        }
			lam_min = lam_max-dlam;
			(void) mis_match(lam_min,&mis_new,(POINTER)lin_pert);
			if (mis_old*mis_new < 0)
			{
	                    if (find_root(mis_match,(POINTER)
	                           lin_pert,0.0,&lam[k][j],lam_min,lam_max,
	                           FACC,lam_max*XACC) == FUNCTION_FAILED)
	                    {
	                        screen("ERROR in rt_single_mode(), ");
	                        screen("find_root for mis_match failed\n");
	                        clean_up(ERROR);
	                    }
	                    break;
		        }
			lam_max = lam_min;
			mis_old = mis_new;
		    }
		    if (i <= n_intvl)	break;
	        }
		if (n_intvl > max_intvl) 	flag_root[k] = 0;
	    }
	    if (j == 1)
	    {
		factor *= 2.0;
		free(yyy_tmp);
	    }
	}
	i = 0;
	for (k = 1; k <= num_unstable; ++k)		i += flag_root[k];
	if (i == 0)
	{
	    screen("\nERROR: in rt_single_mode().\n"
	           "There should be %d unstable mode(s) for ",num_unstable);
	    screen("the wave number |k| = %g, "
	           "but it cannot find any.\n"
	           "Possible cause---grid too coase.\n",wv_num);
	    clean_up(ERROR);
	}
	lambda = 0.0;
	i_infc = -1;
	for (k = 1; (k <= num_unstable) && (flag_root[k] > 0); ++k)
	{
	    tmp = (4.0*lam[k][2]-lam[k][1])/3.0; /* the error is (h_z)^4 */
	    if (tmp > lambda)
	    {
	    	lambda = tmp;
	    	i_infc = k;
	    }
	}
	if (i_infc == -1)
	{
	    screen("ERROR in rt_single_mode(), can't set i_infc\n");
	    clean_up(ERROR);
	}
	lin_pert->normal_mode->sigma_r = sqrt(lambda);
	lin_pert->normal_mode->sigma_i = 0.0;
	free(il);
	free(flag_root);
	free(lam);

	/* allocate storage and store eigenfunctions */
	uni_array(&aa,nl+1,sizeof(double *));
	uni_array(&bb,nl+1,sizeof(double *));
	for (k = 1, i = 1; i <= nl; ++i)
	{
		nstep = Rtp(i)->lin_pert_intvl;
		uni_array(&aa[i],nstep+1,FLOAT);
		uni_array(&bb[i],nstep+1,FLOAT);
		for (j = 0; j <= nstep; ++j, ++k)
		{
			aa[i][j] = yyy[i_infc][1][k];
			bb[i][j] = yyy[i_infc][2][k];
		}
	}
	lin_pert->normal_mode->a = aa;
	lin_pert->normal_mode->b = bb;
	free(yyy_tmp);

	if (debugging("lin_pert"))
	{
	    (void) printf("\nThe wave number is %g.\n",wv_num);
	    (void) printf("lambda is %g.\n",lambda);
	    for (i = 1; i < nl; ++i)
	    {
	    	nstep = Rtp(i)->lin_pert_intvl;
	    	a = lin_pert->normal_mode->a[i];
	    	(void) printf("On the %d%s interface, A is %g.\n",
			      i,ordinal_suffix(i),a[nstep]);
	    }
	    (void) printf("\n");
	}
		

	/* scale the eigenfunction */
	max_a = 0.0;
	for (i = 1; i < nl; ++i)
	{
		nstep = Rtp(i)->lin_pert_intvl;
		a = lin_pert->normal_mode->a[i];
		tmp = fabs(a[nstep]);
		if (tmp > max_a)	max_a = tmp;
	}
	if (max_a == 0.0)
	{
		screen("\nERROR: in rt_single_mode().\n");
		screen("Maximum perturbed normal velocity is zero.\n");
		clean_up(ERROR);
	}

	tmp = lin_pert->normal_mode->sigma_r*
	      (lin_pert->normal_mode->amp_max)/max_a;
	for (i = 1; i <= nl; ++i)
	{
	    nstep = Rtp(i)->lin_pert_intvl;
	    a = lin_pert->normal_mode->a[i];
	    b = lin_pert->normal_mode->b[i];
	    for (j = 0; j <= nstep; ++j)
	    	a[j] *= tmp, 	b[j] *= tmp;
	}
	return	YES;
}		/*end rt_single_mode*/



LOCAL	boolean	mis_match(
	double		lam,
	double		*mismatch,
	POINTER		pntr)
{
	int 		j, k, kp, ne, nb, m, il, kk;
	double 		**s, **s_tmp, ***c, ***c_tmp, **yy;
	LIN_PERT	*lin_pert = (LIN_PERT *)pntr;

	m = lin_pert->tot_pts;
	ne = 2;		nb = 1;
	bi_array(&s_tmp,ne,2*ne+1,DOUBLE);
	s = s_tmp-1;
	for (j = 1; j <= ne; ++j)
	    --s[j];
	tri_array(&c_tmp,ne,ne-nb+1,m+1,DOUBLE);
	c = c_tmp-1;
	for (j = 1; j <= ne; ++j)
	{
	    --c[j];
	    for (k = 1; k <= ne-nb+1; ++k)
		--c[j][k];
	}
	yy = (double **)lin_pert->pointer;
	il = lin_pert->intfc_label;

	difeq(1,il,lam,s,&kk,lin_pert);
	pinvs(2,2,3,5,1,1,c,s);
	for (k = 2; k <= m; ++k) 
	{
	    kp = k-1;
	    difeq(k,il,lam,s,&kk,lin_pert);
	    red(1,2,1,1,2,2,5,2,1,2,kp,c,s);
	    pinvs(1,2,2,5,1,k,c,s);
	}
	difeq(m+1,il,lam,s,&kk,lin_pert);
	red(1,1,3,3,4,4,5,2,1,2,m,c,s);
	pinvs(1,1,4,5,2,m+1,c,s);
	bksub(ne,nb,2,1,m,c);

	for (j = 1; j <= ne; ++j) for (k = 1; k <= m; ++k)
	    yy[j][k] = c[j][1][k];

	free(c_tmp);
	free(s_tmp);

	*mismatch = yy[1][kk]-yy[1][kk-1];
	return	FUNCTION_SUCCEEDED;
}		/*end mis_match*/


LOCAL	void	difeq(
	int		k,
	int		il,
	double		lambda,
	double		**s,
	int		*kk,
	LIN_PERT	*lin_pert)
{
	int		i, ll, rm, nstep;
	INIT_DATA	*init = lin_pert->init;
	int	        nl = lin_pert->layer_sys->num_layers; 
	int		dim = lin_pert->layer_sys->front->rect_grid->dim;
	double		ksq = lin_pert->normal_mode->ksq; 
	double		g_z = lin_pert->g_z;
	double		tmpa, tmpb, csq, rho, rho_prime, q;
	double		tmp, rho_p, rho_m, surf_ten;
	double		coords[MAXD], ub, lb, dz;
	Locstate	st; 
	Front		*front = lin_pert->layer_sys->front;
	LAYER		**layer = lin_pert->layer_sys->layer;

	if (k == 1)
	{
	    s[2][3] = 1.0;
	    s[2][4] = 0.0;
	    s[2][5] = 0.0;
	    return;
	}

	/* find the layer label for the point k and its remain */
	rm = k;
	for (ll = 1; ll <= nl; ++ll)
	{
	    i = Rtp(ll)->lin_pert_intvl + 1;
	    if (rm > i) 	rm -= i;
	    else		break;
	}

	if (ll > nl)
	{
	    s[1][3] = 1.0;
	    s[1][4] = 0.0;
	    s[1][5] = 0.0;
	    return;
	}

	/* find the physical quantities for the point k */
	alloc_state(front->interf,&st,front->sizest);
	ub = layer[ll]->upper_surf->pbar[dim-1];
	lb = layer[ll]->lower_surf->pbar[dim-1];
	nstep = Rtp(ll)->lin_pert_intvl;
	dz = (ub-lb)/nstep;
	for (i = 0; i <= dim-2; ++i)	
	    coords[i] = Rtp(ll)->rt_kh->ref_coords[i];
	if (rm > 1)
	{
	    coords[dim-1] = lb + (rm-1.5)*dz;
	    Get_tgas_state(coords,st,Ct(ll),front->interf,init);
	    rho = Dens(st);
	    csq = sound_speed_squared(st);
	    if (csq <= 0.0)
	    {
	    	screen("ERROR: in difeq().\n");
	    	print_general_vector("At point ",coords,dim,"");
	    	screen(", the sound speed squared is non-positive.\n");
		verbose_print_state("st",st);
	    	clean_up(ERROR);
	    }
	    rho_prime = get_rho_prime(st, Ct(ll),g_z);
	    q = g_z*(rho_prime/rho-g_z/csq);

	    tmpa = (ksq/lambda+1.0/csq)/rho;
	    tmpb = rho*(lambda+q);
	    s[1][1] = -1.0+dz*g_z/csq*0.5;
	    s[1][2] = dz*0.5*tmpa;
	    s[1][3] = 1.0+dz*g_z/csq*0.5;
	    s[1][4] = dz*0.5*tmpa;
	    s[1][5] = 0.0;
	    s[2][1] = dz*0.5*tmpb;
	    s[2][2] = -1.0-dz*g_z/csq*0.5;
	    s[2][3] = dz*0.5*tmpb;
	    s[2][4] = 1.0-dz*g_z/csq*0.5;
	    s[2][5] = 0.0;
	}
	else if ((rm == 1) && (ll > 1))
	{
	    surf_ten = layer[ll-1]->upper_surf->surf_ten;
	    coords[dim-1] = lb;
	    Get_tgas_state(coords,st,Ct(ll),front->interf,init);
	    rho_p = Dens(st);
	    Get_tgas_state(coords,st,Ct(ll-1),front->interf,init);
	    rho_m = Dens(st);

	    tmp = (rho_p-rho_m)*g_z+surf_ten*ksq;
	    if (ll == il + 1)
	    {
	    	s[1][1] = 0.5;
	    	s[1][2] = 0.0;
	    	s[1][3] = 0.5;
	    	s[1][4] = 0.0;
	    	s[1][5] = 1.0;
	    	*kk = k;
	    }
	    else
	    {
	    	s[1][1] = -1.0;
	    	s[1][2] = 0.0;
	    	s[1][3] = 1.0;
	    	s[1][4] = 0.0;
	    	s[1][5] = 0.0;
	    }
	    s[2][1] = 0.5*tmp;
	    s[2][2] = -1.0;
	    s[2][3] = 0.5*tmp;
	    s[2][4] = 1.0;
	    s[2][5] = 0.0;
	}
	free(st);

}		/*end difeq*/


LOCAL	void 	bksub(
	int		ne,
	int		nb,
	int		jf,
	int		k1,
	int		k2,
	double		***c)
{
	int 		nbf,im,kp,k,j,i;
	double 		xx;

	nbf = ne-nb;
	im = 1;
	for (k = k2; k >= k1; --k) 
	{
		if (k == k1) im = nbf+1;
		kp = k+1;
		for (j = 1; j <= nbf; ++j)
		{
			xx = c[j][jf][kp];
			for (i = im; i <= ne; ++i)
				c[i][jf][k] -= c[i][j][k]*xx;
		}
	}
	for (k = k1; k <= k2; ++k)
	{
		kp = k+1;
		for (i = 1; i <= nb; ++i)
			c[i][1][k] = c[i+nbf][jf][k];
		for (i = 1; i <= nbf; ++i)
			c[i+nb][1][k] = c[i][jf][kp];
	}
}		/*end bksub*/

LOCAL	void 	pinvs(
	int		ie1,
	int		ie2,
	int		je1,
	int		jsf,
	int		jc1,
	int		k,
	double		***c,
	double		**s)
{
	int 	js1, jpiv, jp, je2, jcoff, j, irow, ipiv, id, icoff, i; 
	int	*indxr, *indxr_tmp;
	double 	pivinv, piv, dum, big, *pscl, *pscl_tmp;

	uni_array(&indxr_tmp,ie2-ie1+1,INT);
	indxr = indxr_tmp-ie1;
	uni_array(&pscl_tmp,ie2-ie1+1,DOUBLE);
	pscl = pscl_tmp-ie1;

	je2 = je1+ie2-ie1;
	js1 = je2+1;
	for (i = ie1; i <= ie2; ++i) 
	{
	    big = 0.0;
	    for (j = je1; j <= je2; ++j)
	    	if (fabs(s[i][j]) > big) big = fabs(s[i][j]);
	    if (big == 0.0) 
	    {
	        pscl[i] = 0.0;
	        screen("ERROR in pinvs(), Singular matrix - row all 0.\n");
	        clean_up(ERROR);
	    }
	    else
	        pscl[i] = 1.0/big;
	    indxr[i] = 0;
	}
	for (id = ie1; id <= ie2; ++id) 
	{
	    piv = 0.0;
	    ipiv = jpiv = jp = -1;
	    for (i = ie1; i <= ie2; ++i)
	    {
	    	if (indxr[i] == 0) 
	    	{
	    	    big = 0.0;
	    	    for (j = je1; j <= je2; ++j) 
	    	    {
	    	    	if (fabs(s[i][j]) > big) 
	    	    	{
	    	    	    jp = j;
	    	    	    big = fabs(s[i][j]);
	    	    	}
	    	    }
	    	    if (big*pscl[i] > piv) 
	    	    {
	    	    	ipiv = i;
	    	    	jpiv = jp;
	    	    	piv = big*pscl[i];
	    	    }
	    	}
	    }
	    if ((ipiv == -1) || (jp == -1) || (jpiv == -1))
	    {
	    	screen("\nERROR in pinvs(), can set indicies.\n");
	    	clean_up(ERROR);
	    }
	    if (s[ipiv][jpiv] == 0.0) 
	    {
	    	screen("\nERROR in pinvs(), Singular bi_array.\n");
	    	clean_up(ERROR);
	    }
	    indxr[ipiv] = jpiv;
	    pivinv = 1.0/s[ipiv][jpiv];
	    for (j = je1; j <= jsf; ++j) s[ipiv][j] *= pivinv;
	    s[ipiv][jpiv] = 1.0;
	    for (i = ie1; i <= ie2; ++i) 
	    {
	    	if (indxr[i] != jpiv) 
	    	{
	    	    if (s[i][jpiv])
	    	    {
	    	    	dum = s[i][jpiv];
	    	    	for (j = je1; j <= jsf; ++j)
	    	    	    s[i][j] -= dum*s[ipiv][j];
	    	    	s[i][jpiv] = 0.0;
	    	    }
	    	}
	    }
	}
	jcoff = jc1-js1;
	icoff = ie1-je1;
	for (i = ie1; i <= ie2; ++i) 
	{
	    irow = indxr[i]+icoff;
	    for (j = js1; j <= jsf; ++j) 
	    	c[irow][j+jcoff][k] = s[i][j];
	}
	free(pscl_tmp);
	free(indxr_tmp);
}		/*end pinvs*/

LOCAL	void 	red(
	int		iz1,
	int		iz2,
	int		jz1,
	int		jz2,
	int		jm1,
	int		jm2,
	int		jmf,
	int		ic1,
	int		jc1,
	int		jcf,
	int		kc,
	double		***c,
	double		**s)
{
	int 		loff,l,j,ic,i;
	double 		vx;

	loff = jc1-jm1;
	ic = ic1;
	for (j = jz1; j <= jz2; ++j) 
	{
		for (l = jm1; l <= jm2; ++l) 
		{
			vx = c[ic][l+loff][kc];
			for (i = iz1; i <= iz2; ++i)
				s[i][l] -= s[i][j]*vx;
		}
		vx = c[ic][jcf][kc];
		for (i = iz1; i <= iz2; ++i)
			s[i][jmf] -= s[i][j]*vx;
		ic += 1;
	}
}		/*end red*/


/*ARGSUSED*/
LOCAL	int	kh_single_mode(
	LIN_PERT	*lin_pert)
{
	screen("\nERROR: in kh_single_mode().\n");
	screen("This function has not been written yet.\n");
	clean_up(ERROR);
	return	YES;
}		/*end kh_single_mode*/

EXPORT	double	get_rho_prime(
	Locstate	st,
	COMP_TYPE	*ct,
	double		g_z)
{
	STRATIFICATION_TYPE s_type;
	double		    csq;

	if (ct->type == RT_KH)
	    s_type = Rt_kh(ct)->stratification_type;
	else if (ct->type == RT_PERTURBED)
	    s_type = Rt_perturbed(ct)->rt_kh->stratification_type;
	else
	{
	    screen("ERROR in get_rho_prime(), invalid comp type "
		   "ct->type = %d\n",ct->type);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}

	switch (s_type)
	{
	case CONSTANT:
	case CONSTANT_DENSITY:
	case HYDRO_STATIC:
	    return 0.0;
	case ISOTHERMAL:
	    return sqr(Dens(st))*g_z*K_T(st);
	case ADIABATIC:
	    csq = sound_speed_squared(st);
	    return Dens(st)*g_z/csq;
	default:
	    screen("\nERROR: in get_rho_prime().\n");
	    screen("Unrecoganized thermodynamic type.");
	    clean_up(ERROR);
	}
	return	ERROR_FLOAT;
}		/*end get_rho_prime*/



LOCAL	void	set_kh_sine_perturbed_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front,
	int		num_modes,
	PERT_TYPE	init_type)
{
	SINE_PERT	*rsp;

	if (comp_type->type == KH_SINE_PERTURBED) /*ALREADY SET*/
		return;

	comp_type->type = KH_SINE_PERTURBED;
	rsp = alloc_sine_pert_structure(num_modes,front);
	comp_type->extra = (POINTER) rsp;
	rsp->init_type = init_type;
	rsp->pert_bdry_type[0] = UNMODIFIED;
	rsp->pert_bdry_type[1] = UNMODIFIED;
	rsp->pert_bdry_type[2] = UNMODIFIED;

	comp_type->_get_state = get_KH_sine_perturbed_state;
	comp_type->free_comp_type_extra = free_kh_rt_sine_perturbed_comp_type;
}		/*end set_kh_sine_perturbed_comp_type*/

LOCAL	void	set_random_surface_perturbed_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front,
	int		num_modes,
	PERT_TYPE	init_type)
{
	SINE_PERT	*rsp;

	if (comp_type->type == RANDOM_SURFACE_PERTURBED) /*ALREADY SET*/
		return;

	comp_type->type = RANDOM_SURFACE_PERTURBED;
	rsp = alloc_sine_pert_structure(num_modes,front);
	comp_type->extra = (POINTER)rsp;
	rsp->init_type = init_type;
	rsp->pert_bdry_type[0] = UNMODIFIED;
	rsp->pert_bdry_type[1] = UNMODIFIED;
	rsp->pert_bdry_type[2] = UNMODIFIED;

	comp_type->_get_state = get_state_random_surface_perturbed;
	comp_type->free_comp_type_extra = free_kh_rt_sine_perturbed_comp_type;
}		/*end set_random_surface_perturbed_comp_type*/

/*ARGSUSED*/
LOCAL	void	get_state_random_surface_perturbed(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*comp_type,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	debug_print("init_states","Entered get_state_random_surface_perturbed()\n");
	random_surface_perturbed_state(state,coords,initial_time(init),
				       comp_type->comp,stype);
}		/*end get_state_random_surface_perturbed*/

LOCAL	void	free_kh_rt_sine_perturbed_comp_type(
	COMP_TYPE	*comp_type)
{
	SINE_PERT	*rsp;

	if (comp_type->type != KH_SINE_PERTURBED)
		return;
	rsp = (SINE_PERT*)comp_type->extra;
	if (rsp != NULL)
		free(rsp);
	comp_type->extra = NULL;
}		/*end free_kh_rt_sine_perturbed_comp_type*/

LOCAL	void	set_rt_perturbed_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	_RT_KH		*rtkh = NULL;
	_RT_PERTURBED 	*rtp;

	if (comp_type->type == RT_PERTURBED) /*ALREADY SET*/
	    return;

	if (comp_type->type == RT_KH)
	    rtkh = Rt_kh(comp_type);

	if (rtkh == NULL)
	{
	    set_rt_kh_comp_type(comp_type,front);
	    rtkh = Rt_kh(comp_type);
	}
	comp_type->type = RT_PERTURBED;
	scalar(&rtp, sizeof(_RT_PERTURBED));
	rtp->rt_kh = rtkh;
	comp_type->extra = (POINTER)rtp;
	comp_type->_get_state = get_state_rt_kh_perturbed;
	comp_type->free_comp_type_extra = free_rt_perturbed_comp_type;

}		/*end set_rt_perturbed_comp_type*/

LOCAL	void	free_rt_perturbed_comp_type(
	COMP_TYPE	*comp_type)
{
	_RT_PERTURBED 	*rtp;
	if (comp_type->type != RT_PERTURBED)
		return;

	rtp = Rt_perturbed(comp_type);
	if (rtp == NULL)
		return;

	comp_type->type = RT_KH;
	comp_type->extra = (POINTER)rtp->rt_kh;
	free_rt_kh_comp_type(comp_type);
	free(rtp);
	comp_type->extra = NULL;
}		/*end free_rt_perturbed_comp_type*/

LOCAL	void	set_kh_perturbed_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	_RT_KH		*rtkh = NULL;

 	if (comp_type->type == KH_PERTURBED) /*ALREADY SET*/
	    return;

	if (comp_type->type == RT_KH)
	    rtkh = Rt_kh(comp_type);

	if (rtkh == NULL)
	{
	    set_rt_kh_comp_type(comp_type,front);
	    rtkh = Rt_kh(comp_type);
	}
	comp_type->type = KH_PERTURBED;
	comp_type->_get_state = get_state_rt_kh_perturbed;
	comp_type->free_comp_type_extra = free_kh_perturbed_comp_type;

}		/*end set_kh_perturbed_comp_type*/

LOCAL	void	free_kh_perturbed_comp_type(
	COMP_TYPE	*comp_type)
{
	if (comp_type->type != KH_PERTURBED)
	    return;

	comp_type->type = RT_KH;
	free_rt_kh_comp_type(comp_type);
	comp_type->extra = NULL;
}		/*end free_kh_perturbed_comp_type*/
#endif /* defined(TWOD) || defined(THREED) */
