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
*			spolars.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if defined(TWOD) || defined(THREED)
#include <gdecs/gdecs.h>
#include <ginit/ginit.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>

#if defined(__cplusplus)
LOCAL const int NUM_POINTS = 100;
LOCAL const int MAX_NUM_POINTS = 400;

LOCAL const unsigned int DONT_PLOT	  = 0x00;
LOCAL const unsigned int PLOT_SHOCK_POLAR = 0x01;
LOCAL const unsigned int PLOT_RAREF_POLAR = 0x02;
LOCAL const unsigned int PLOT_CPLUS_WAVE  = 0x04;
LOCAL const unsigned int PLOT_CMINUS_WAVE = 0x08;
LOCAL const unsigned int PLOT_SONIC_LOCUS = 0x10;
LOCAL const unsigned int PLOT_SUPERSONIC  = 0x20;
LOCAL const unsigned int PLOT_SUBSONIC	  = 0x40;

#else /* defined(__cplusplus) */
enum {
	NUM_POINTS = 100,
	MAX_NUM_POINTS = 400};

enum {
	DONT_PLOT	 = 0x00,
	PLOT_SHOCK_POLAR = 0x01,
	PLOT_RAREF_POLAR = 0x02,
	PLOT_CPLUS_WAVE	 = 0x04,
	PLOT_CMINUS_WAVE = 0x08,
	PLOT_SONIC_LOCUS = 0x10,
	PLOT_SUPERSONIC	 = 0x20,
	PLOT_SUBSONIC	 = 0x40};
#endif /* defined(__cplusplus) */

#define plot_shock_polar(plt_ctrl)		((plt_ctrl) & PLOT_SHOCK_POLAR)
#define plot_rarefaction_polar(plt_ctrl)	((plt_ctrl) & PLOT_RAREF_POLAR)
#define plot_cplus_wave(plt_ctrl)		((plt_ctrl) & PLOT_CPLUS_WAVE)
#define plot_cminus_wave(plt_ctrl)		((plt_ctrl) & PLOT_CMINUS_WAVE)
#define plot_sonic_locus(plt_ctrl)		((plt_ctrl) & PLOT_SONIC_LOCUS)
#define plot_supersonic(plt_ctrl)		((plt_ctrl) & PLOT_SUPERSONIC)
#define plot_subsonic(plt_ctrl)			((plt_ctrl) & PLOT_SUBSONIC)

enum {MAX_NUM_VAR = 15};

enum {
	TURN_ANG	= 0,
	PRESSURE	= 1,
	SPEC_VOL	= 2,
	DENSITY		= 3,
	VELX		= 4,
	VELY		= 5,
	SOUND_SPEED	= 6,
	FLOW_SPEED	= 7,
	MACH_NUMBER	= 8,
	MASS_FLUX	= 9,
	ADIABATIC_GAMMA = 10,
	GRUNEISEN_GAMMA = 11,
	TWO_D_STABILITY = 12,
	ENTROPY		= 13,
	ENTHALPY	= 14};

LOCAL	size_t sizest;

LOCAL	int 	use_logp_increments = YES;

	/* LOCAL Function Declarations */
LOCAL	double	turning_angle(double,double,Locstate);
LOCAL	double	prompt_for_rp_surface_tension(Locstate,Locstate);
LOCAL	int	compute_oned_rp(INIT_DATA*,INIT_PHYSICS*);
LOCAL	int	compute_refracted_vorticity(INIT_DATA*,INIT_PHYSICS*);
LOCAL	int	init_plot_shock_crossing(INIT_DATA*,INIT_PHYSICS*,double*,
					 Locstate*,int*,int*,double*,int*);
LOCAL	int	init_plot_shock_overtake(INIT_DATA*,INIT_PHYSICS*,double*,
					 Locstate*,int*,int*,double*,int*);
LOCAL	int	init_plot_shock_refraction(INIT_DATA*,INIT_PHYSICS*,
					   double*,Locstate*,int,int*,int*,
					   double*,int*);
LOCAL	int	init_plot_shock_transmission(INIT_DATA*,INIT_PHYSICS*,
					     Locstate*,int*,int*,double*,int*);
LOCAL	int	init_plot_sonic_locus(INIT_DATA*,INIT_PHYSICS*,Locstate*,
				      int*,int*,double*,int*);
LOCAL	int	init_plot_total_reflection(INIT_DATA*,INIT_PHYSICS*,double*,
					   Locstate*,int*,int*,double*,int*);
LOCAL	int	init_plot_wall_reflection(INIT_DATA*,INIT_PHYSICS*,double*,
					  Locstate*,int*,int*,double*,int*);
LOCAL	int	init_plot_shock_polars(INIT_DATA*,INIT_PHYSICS*,double*,
				       Locstate*,int*,int*,double*,int*);
LOCAL	double	min_vorticity_for_regular_refraction(RP_DATA*,double**,
                                                     double*,NODE_FLAG);
LOCAL	int	normal_shock_refraction(INIT_DATA*,INIT_PHYSICS*,Locstate*);
LOCAL	int	npt_w_speed_test(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	calculate_sonic_locus(Locstate,Locstate,double**,double*,double*,
				      int,int);
LOCAL	void	init_shock_polar_ahead_velocity(Locstate,double*,int*,
						double*,const char*);
LOCAL	void	init_state_across_contact(Locstate,Locstate,Gas_param*,double*,
					  const char*);
LOCAL	void	plot_polar(double,double,Locstate,double**,double,double,double*,
			   double*,int*,int);
LOCAL	void	plot_shock_refraction_vs_incident_angle(INIT_DATA*,
							INIT_PHYSICS*,int);
LOCAL	void	print_line(int,int,int,int,double**,double*,double*);
LOCAL	void	print_plot_header(const char**,size_t,double*,double*);
LOCAL	void	set_graph_line(double*,double*,double*,double*,Locstate,double,
			       double,double,double,double,double,double,
			       double,double,int);
LOCAL	void	set_plot_control(const char*,int*,int*,Locstate);
LOCAL	void	spolars_clean_up(void);
LOCAL	void	start_up(int,char**,INIT_DATA*);

LOCAL char *temporary_input_file = NULL;

#if !defined(__INTEL_COMPILER)
#pragma	noinline	turning_angle
#pragma	noinline	init_plot_shock_crossing
#pragma	noinline	init_plot_shock_overtake
#pragma	noinline	init_plot_shock_refraction
#pragma	noinline	init_plot_shock_transmission
#pragma	noinline	init_plot_sonic_locus
#pragma	noinline	init_plot_total_reflection
#pragma	noinline	init_plot_wall_reflection
#pragma	noinline	init_plot_shock_polars
#pragma	noinline	calculate_sonic_locus
#pragma	noinline	init_shock_polar_ahead_velocity
#pragma	noinline	init_state_across_contact
#pragma	noinline	normal_shock_refraction
#pragma	noinline	plot_polar
#pragma	noinline	plot_shock_refraction_vs_incident_angle
#pragma	noinline	print_line
#pragma	noinline	print_plot_header
#pragma	noinline	set_graph_line
#pragma	noinline	set_plot_control
#pragma	noinline	start_up
#endif /*!defined(__INTEL_COMPILER)*/

/*ARGSUSED*/
int sp_main(
	int		argc,
	char**		argv,
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	IMPORT	boolean	suppress_prompts;
	static		RECT_GRID	Rgr;
	char		s[Gets_BUF_SIZE];
	Locstate	state[3];
	double		theta0[3], Msq[3];
	double		pmax[3], pmin[3];
	double		**var[3];
	double		var_l[MAX_NUM_VAR], var_u[MAX_NUM_VAR];
	int		nfloats;
	int		max_pressure_given[3];
	int		plot_control[3];
	int		num_polars;
	int		npts[3];
	int		i;
	int		dim = 2;
	int		problem;
	int		SHOCK_POLAR;
	int		REG_REFL;
	int		MACH_REFL;
	int		TRANSMISSION;
	int		DIFFRACTION;
	int		SHOCK_OVERTAKE;
	int		SHOCK_CROSSING;
	int		SONIC_POINT_LOCUS;
	int		TOT_INT_REFL;
	int		ONED_RP;
	int		NPT_W_SPEED;
	int		VORTICITY;

	static const char *polar_header[] = { "TURNING_ANGLE",
					      "PRESSURE",
					      "SPECIFIC_VOLUME",
					      "DENSITY",
					      "VELX",
					      "VELY",
					      "SOUND_SPEED",
					      "FLOW_SPEED",
					      "MACH_NUMBER",
					      "MASS_FLUX",
					      "ADIABATIC_GAMMA",
					      "GRUNEISEN_GAMMA",
					      "TWO_D_STABILITY",
					      "ENTROPY",
					      "ENTHALPY",
					      NULL};

	suppress_prompts = NO;
	start_up(argc,argv,init);

	ip->root->grid->rect_grid = &Rgr;
	ip->root->front->rect_grid = &Rgr;
	ip->root->wave->rect_grid = &Rgr;
	Rgr.dim = dim;

	g_compute_sizest(PURE_NON_REACTIVE,&sizest,&nfloats,dim);
	g_set_sizeof_state(ip->root,sizest,nfloats);

	set_interface_hooks(dim,init);
	i_intfc(init) = make_interface(dim);

	for (i = 0; i < 3; ++i)
	{
	    max_pressure_given[i] = NO;
	    g_alloc_state(&state[i],sizeof(VGas));
	    set_type_of_state(state[i],VGAS_STATE);
	    npts[i] = NUM_POINTS;
	    theta0[i] = 0.0;
	    pmin[i] = HUGE_VAL;
	    pmax[i] = -HUGE_VAL;
	}

	for (i = 0; i < MAX_NUM_VAR; ++i)
	{
	    var_l[i] = HUGE_VAL;
	    var_u[i] = -HUGE_VAL;
	}

	set_binary_output(NO);
	i = 0;
	screen("The following shock polar configurations are available\n");
	screen("\tA single shock polar (S, or %d)\n",i);
	SHOCK_POLAR = i++;
	screen("\tRegular reflection (R, or %d)\n",i);
	REG_REFL = i++;
	screen("\tMach reflection (M, or %d)\n",i);
	MACH_REFL = i++;
	screen("\tTransmission node configuration (T, or %d)\n",i);
	TRANSMISSION = i++;
	screen("\tDiffraction node configuration (D, or %d)\n",i);
	DIFFRACTION = i++;
	screen("\tShock overtake configuration (O, or %d)\n",i);
	SHOCK_OVERTAKE = i++;
	screen("\tShock crossing configuration (C, or %d)\n",i);
	SHOCK_CROSSING = i++;
	screen("\tTotal internal reflection (TIR, or %d)\n",i);
	TOT_INT_REFL = i++;
	screen("\tSonic point locus (SPL, or %d)\n",i);
	SONIC_POINT_LOCUS = i++;
	screen("\tOne Dimensional Riemann Problem (RP, or %d)\n",i);
	ONED_RP = i++;
	screen("\tWave speed test(3PT, or %d)\n",i);
	NPT_W_SPEED = i++;
	screen("\tVorticity generated by refractions (V or %d)\n",i);
	VORTICITY = i++;
	screen("Enter the desired configuration [S or 0]: ");
	(void) Gets(s);
	problem = SHOCK_POLAR;
	if (s[0] != '\0')
	{
	    if (strcasecmp(s,"R") == 0)
	    	problem = REG_REFL;
	    else if (strcasecmp(s,"M") == 0)
	    	problem = MACH_REFL;
	    else if (strcasecmp(s,"T") == 0)
	    	problem = TRANSMISSION;
	    else if (strcasecmp(s,"D") == 0)
	    	problem = DIFFRACTION;
	    else if (strcasecmp(s,"O") == 0)
	    	problem = SHOCK_OVERTAKE;
	    else if (strcasecmp(s,"C") == 0)
	    	problem = SHOCK_CROSSING;
	    else if (strcasecmp(s,"TIR") == 0)
	    	problem = TOT_INT_REFL;
	    else if (strcasecmp(s,"SPL") == 0)
	    	problem = SONIC_POINT_LOCUS;
	    else if (strcasecmp(s,"RP") == 0)
	    	problem = ONED_RP;
	    else if (strcasecmp(s,"3PT") == 0)
	    	problem = NPT_W_SPEED;
	    else if (strcasecmp(s,"V") == 0)
	    	problem = VORTICITY;
	    else if (isdigit(s[0]))
	    	(void) sscanf(s,"%d",&problem);
	}
	if (problem == REG_REFL)
	    num_polars = init_plot_wall_reflection(init,ip,theta0,state,
				                   npts,plot_control,pmax,
						   max_pressure_given);
	else if (problem == MACH_REFL)
	    num_polars = init_plot_wall_reflection(init,ip,theta0,state,
				                   npts,plot_control,pmax,
						   max_pressure_given);
	else if (problem == TRANSMISSION)
	    num_polars = init_plot_shock_transmission(init,ip,state,
				                      npts,plot_control,
				                      pmax,max_pressure_given);
	else if (problem == DIFFRACTION)
	    num_polars = init_plot_shock_refraction(init,ip,theta0,state,dim,
						    npts,plot_control,pmax,
						    max_pressure_given);
	else if (problem == SHOCK_OVERTAKE)
	    num_polars = init_plot_shock_overtake(init,ip,theta0,state,
				                  npts,plot_control,pmax,
						  max_pressure_given);
	else if (problem == SHOCK_CROSSING)
	    num_polars = init_plot_shock_crossing(init,ip,theta0,state,
				                  npts,plot_control,pmax,
						  max_pressure_given);
	else if (problem == SONIC_POINT_LOCUS)
	    num_polars = init_plot_sonic_locus(init,ip,state,npts,plot_control,
					       pmax,max_pressure_given);
	else if (problem == TOT_INT_REFL)
	    num_polars = init_plot_total_reflection(init,ip,theta0,state,npts,
						    plot_control,pmax,
						    max_pressure_given);
	else if (problem == ONED_RP)
	    return compute_oned_rp(init,ip);
	else if (problem == NPT_W_SPEED)
	    return npt_w_speed_test(init,ip);
	else if (problem == VORTICITY)
	    compute_refracted_vorticity(init,ip);
	else
	    num_polars = init_plot_shock_polars(init,ip,theta0,state,npts,
						plot_control,pmax,
						max_pressure_given);

	if (num_polars == 0)
	{
	    clean_up(0);
	    return 0;
	}
	screen("Use log pressure increments [%s]?: ",
	    (use_logp_increments)?"yes":"no");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    use_logp_increments = YES;
	else if (s[0] == 'n' || s[0] == 'N')
	    use_logp_increments = NO;


	/* Set default pressure window */
	
	for (i = 0; i < num_polars; ++i)
	{
	    Msq[i] = mach_number_squared(state[i],NULL,NULL);
	    if (plot_shock_polar(plot_control[i]))
    	    {
	    	pmin[i] = min(pmin[i],pressure(state[i]));
	    	if (!max_pressure_given[i])
	    	{
	    	    pmax[i] = max(pmax[i],max_behind_shock_pr(Msq[i],state[i]));
	    	}
	    	var_l[PRESSURE] = min(var_l[PRESSURE],pmin[i]);
	    	var_u[PRESSURE] = max(var_u[PRESSURE],pmax[i]);
	    }
	    if (plot_rarefaction_polar(plot_control[i]))
	    {
	    	pmin[i] = min(pmin[i],0.0);
	    	pmax[i] = max(pmax[i],pressure(state[i]));
	    	var_l[PRESSURE] = min(var_l[PRESSURE],pmin[i]);
	    	var_u[PRESSURE] = max(var_u[PRESSURE],pmax[i]);
	    }
	}

	if (debugging("spolars"))
	{
	    for (i = 0; i < num_polars; ++i)
	    {
	    	char mesg[80];

	    	(void) printf("Base state %d, Msq[%d] = %g\n",i,i,Msq[i]);
	    	(void) sprintf(mesg,"state%d is %s, M%d = %g",i,
	    		       (Msq[i] >= 1.0) ? "supersonic" :
	    		       "subsonic",i,sqrt(Msq[i]));
	    	verbose_print_state(mesg,state[i]);
	    }
	}

	screen("To specify a lower bound for the pressure,\n");
	screen("\tenter the lower bound (current value = %g): ",
		var_l[PRESSURE]);
	(void) Gets(s);
	if (s[0] != '\0') (void) sscan_float(s,&var_l[PRESSURE]);

	screen("To specify a upper bound for the pressure,\n");
	screen("\tenter the upper bound (current value = %g): ",
		var_u[PRESSURE]);
	(void) Gets(s);
	if (s[0] != '\0') (void) sscan_float(s,&var_u[PRESSURE]);

	for (i = 0; i < num_polars; ++i)
	{
	    bi_array(&var[i],2*npts[i],MAX_NUM_VAR,FLOAT);
	    pmin[i] = max(pmin[i],var_l[PRESSURE]);
	    pmax[i] = min(pmax[i],var_u[PRESSURE]);
	    if (plot_shock_polar(plot_control[i])
	     || plot_rarefaction_polar(plot_control[i]))
	    {
	    	plot_polar(theta0[i],Msq[i],state[i],var[i],pmin[i],pmax[i],
			   var_l,var_u,npts+i,plot_control[i]);
	    }

	    if (plot_sonic_locus(plot_control[i]))
	    {
	    	calculate_sonic_locus(state[0],state[1],var[i],
				      var_l,var_u,npts[i],plot_control[i]);
	    }
	}


	/* Set default turning angle limits */

	screen("To specify a lower limit on the turning angle window,\n"
	       "\tenter the lower limit (current limit = %g): ",
		var_l[TURN_ANG]);
	(void) Gets(s);
	if (s[0] != '\0') (void) sscan_float(s,&var_l[TURN_ANG]);

	screen("To specify a upper limit on the turning angle window,\n"
	       "\tenter the upper limit (current limit = %g): ",
		var_u[TURN_ANG]);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&var_u[TURN_ANG]);
	
	for (i = 0; i < num_polars; ++i)
	{
	    int j;
	    if (plot_cminus_wave(plot_control[i]))
	    {
	    	print_plot_header(polar_header,15,var_l,var_u);
	    	for (j = 0; j < npts[i]; ++j) 
	    	{
	    	    print_line(j,npts[i],MAX_NUM_VAR,15,var[i],var_l,var_u);
	    	}
	    	(void) printf("\n");
	    }
	    if (plot_cplus_wave(plot_control[i]))
	    {
	    	print_plot_header(polar_header,15,var_l,var_u);
	    	for (j = 0; j < npts[i]; ++j) 
	    	{
	    	    print_line(j,npts[i],MAX_NUM_VAR,15,
			       var[i]+npts[i],var_l,var_u);
	    	}
	    	(void) printf("\n");
	    }
	}
	clean_up(0);
	return 0;
}		 /*end main*/

LOCAL int init_plot_wall_reflection(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	double		*theta0,
	Locstate	*state,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	char		s[Gets_BUF_SIZE];
	double		p1;
	double		shock_ang;
	int		is_turn_ang_positive;
	Gas_param	*params;
	static double	node_v[MAXD]; /*Statics initialized to zero*/

	params = init_eos_params(init,ip," ahead of the incident shock",NO);
	prompt_for_thermodynamics(state[0],params,
				  " of the gas ahead of the shock");
	screen("Enter the pressure behind the incident shock wave: ");
	(void) Scanf("%f\n",&p1);

	init_shock_polar_ahead_velocity(state[0],&p1,
		&max_pressure_given[0],&pmax[0],"for the incident gas ");
	screen("Enter the number of points on the incident shock ");
	screen("polar [%d]: ",NUM_POINTS);
	(void) Gets(s);
	if (s[0] != '\0') (void) sscanf(s,"%d",&npts[0]);
	set_plot_control("incident ",&plot_control[0],&npts[0],state[0]);
	screen("Enter p if the flow behind the first incident shock ");
	screen("is turned counter clockwise\n\t");
	screen("relative to the ahead flow, otherwise the flow ");
	screen("will be assumed\n\t");
	screen("to be turned in the clockwise direction.\n");
	screen("Enter choice here: ");
	(void) Gets(s);
	is_turn_ang_positive = (s[0] == 'p' || s[0] == 'P') ? YES : NO;
	Check_return(
	    s_polar_3(state[0],YES,p1,is_turn_ang_positive,NO,node_v,state[1],
		      &shock_ang,&theta0[1]),
	    init_plot_wall_reflection) 
	set_state(state[0],VGAS_STATE,state[0]);
	set_state(state[1],VGAS_STATE,state[1]);
	set_plot_control("reflected ",&plot_control[1],&npts[1],state[1]);
	return 2;
}		/*end init_plot_wall_reflection*/

LOCAL	int init_plot_shock_transmission(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	Locstate	*state,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	Gas_param	*params;
	static double	node_v[MAXD]; /*Statics initialized to zero*/

	params = init_eos_params(init,ip,
				 " ahead of the incident shock wave",NO);
	prompt_for_thermodynamics(state[0],params," of the incident gas ");
	init_shock_polar_ahead_velocity(state[0],
		(double *)NULL,&max_pressure_given[0],&pmax[0],
		"for the incident gas ");
	set_plot_control("incident ",&plot_control[0],&npts[0],state[0]);
	set_type_of_state(state[1],state_type(state[0]));
	params = init_eos_params(init,ip,
				 " ahead of the transmitted shock wave",NO);
	init_state_across_contact(state[1],state[0],params,node_v,
				" of the gas ahead of the transmitted shock");
	set_state(state[0],VGAS_STATE,state[0]);
	set_state(state[1],VGAS_STATE,state[1]);
	set_plot_control("transmitted ",&plot_control[1],&npts[1],state[1]);
	return 2;
}		/*end init_plot_shock_transmission*/

LOCAL	void plot_shock_refraction_vs_incident_angle(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	int		dim)
{
	char		         s[Gets_BUF_SIZE];
	char		         ffmt[20];
	int		         i, j;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	int		         nlines, npts;
	boolean		         is_reflected_shock;
	size_t		         width;
	double		         p0, p1;
	double		         pstarl, pstarr;
	double		         ustarl, ustarr;
	double		         ml, mr;
	double		         beta_l, beta_u;
	double		         beta, dbeta;
	double		         nv;
	double		         m, v0;
	double		         **vars;
	double                    v4, v5, a[2];
	Gas_param	         *params;
	NODE_FLAG	         flag;
	static RP_DATA	         *RP = NULL;
	static int	nvars;
	static double	**t = NULL;
	static double	node_v[MAXD]; /* Initialized to zero */
	static double	nor[3] = {0.0, -1.0, 0.0};
	static Locstate Tsl = NULL, Tsr = NULL;
	enum {
	    INCIDENT_ANGLE			=  0,
	    PSTAR,
	    SHEAR,
	    REFLECTED_ANGLE_1,
	    REFLECTED_ANGLE_2,
	    REFLECTED_ANGLE_3,
	    CONTACT_ANGLE,
	    TRANSMITTED_ANGLE,
	    INCIDENT_TURN_ANGLE,
	    REFLECTED_TURN_ANGLE,
	    TRANSMITTED_TURN_ANGLE,
	    THETA_ERR,
	    AHEAD_INCIDENT_MACH_NUMBER,
	    BEHIND_INCIDENT_MACH_NUMBER,
	    BEHIND_REFLECTED_MACH_NUMBER,
	    BEHIND_TRANSMITTED_MACH_NUMBER,
	    AHEAD_TRANSMITTED_MACH_NUMBER,
	    NUM_GRAPH_HEADERS
	};
	static const char *header[NUM_GRAPH_HEADERS+1];
	header[INCIDENT_ANGLE] = "INCIDENT_ANGLE";
	header[PSTAR] = "PSTAR";
	header[SHEAR] = "SHEAR";
	header[REFLECTED_ANGLE_1] = "REFLECTED_ANGLE_1";
	header[REFLECTED_ANGLE_2] = "REFLECTED_ANGLE_2";
	header[REFLECTED_ANGLE_3] = "REFLECTED_ANGLE_3";
	header[CONTACT_ANGLE] = "CONTACT_ANGLE";
	header[TRANSMITTED_ANGLE] = "TRANSMITTED_ANGLE";
	header[INCIDENT_TURN_ANGLE] = "INCIDENT_TURN_ANGLE";
	header[REFLECTED_TURN_ANGLE] = "REFLECTED_TURN_ANGLE";
	header[TRANSMITTED_TURN_ANGLE] = "TRANSMITTED_TURN_ANGLE";
	header[THETA_ERR] = "THETA_ERR";
	header[AHEAD_INCIDENT_MACH_NUMBER] = "AHEAD_INCIDENT_MACH_NUMBER";
	header[BEHIND_INCIDENT_MACH_NUMBER] = "BEHIND_INCIDENT_MACH_NUMBER";
	header[BEHIND_REFLECTED_MACH_NUMBER] = "BEHIND_REFLECTED_MACH_NUMBER";
	header[BEHIND_TRANSMITTED_MACH_NUMBER]="BEHIND_TRANSMITTED_MACH_NUMBER";
	header[AHEAD_TRANSMITTED_MACH_NUMBER] = "AHEAD_TRANSMITTED_MACH_NUMBER";
	header[NUM_GRAPH_HEADERS] = NULL;

	if (RP == NULL)
	{
	    g_alloc_state(&Tsl,sizeof(VGas));
	    g_alloc_state(&Tsr,sizeof(VGas));
	    for (nvars = 0; header[nvars] != NULL; ++nvars);
	    bi_array(&t,2,2,FLOAT);
	    t[1][0] = 1.0;
	    t[1][1] = 0.0;
	    RP = allocate_RP_DATA_structure(sizeof(VGas),NO,VGAS_STATE);
	}
	RP->ang_dir = COUNTER_CLOCK;

	clear_node_flag(flag);
	use_subsonic_state(flag) = YES;
	(void) prompt_for_eos_params(init,ip,NO,"");
	params = prompt_for_eos_params(init,ip,NO,"\n\tahead of the incident "
						  "shock wave");
	prompt_for_thermodynamics(RP->state[0], params," of the incident gas");
	set_state(RP->state[0],RP->stype,RP->state[0]);
	p0 = pressure(RP->state[0]);
	screen("Enter an optional normal component of velocity for the gas\n"
	       "\tahead of the incident shock [0]: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    double v[MAXD];
	    (void) sscan_float(s,&nv);
	    for (i = 0; i < dim; ++i)
	        v[i] = nv*nor[i];
	    add_velocity_to_state(RP->state[0],v);
	}
	Init_params(RP->state[6],
		    prompt_for_eos_params(init,ip,NO,"\n\tahead of the "
						     "transmitted shock wave"));
	set_type_of_state(RP->state[6],TGAS_STATE);
	for (i = 0; i < dim; ++i)
	    Vel(RP->state[6])[i] = Vel(RP->state[0])[i];
	Press(RP->state[6]) = p0;
        reset_gamma(RP->state[6]);
	screen("Enter the density of the gas ahead of the transmitted shock: ");
	(void) Scanf("%f\n",&Dens(RP->state[6]));
	screen("The slip across the upstream fluid interface is "
	       "defined as the difference\n"
	       "\tbetween the velocity on the transmitted side "
	       "of the interface and\n"
	       "\tthe the velocity on the incident shock side.\n"
	       "Enter an optional slip across the upstream "
	       "fluid interface [0]: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    double v;
	    (void) sscan_float(s,&v);
	    Vel(RP->state[6])[0] += v;
	}
	set_state(RP->state[6],RP->stype,RP->state[6]);

	prompt_for_behind_shock_state(RP->state[0],RP->state[1],NO,nor,
				      RP->stype,YES,init);
	p1 = pressure(RP->state[1]);
	m = mass_flux(p1,RP->state[0]);
	v0 = m/Dens(RP->state[0]);

	beta_l = 0.0;
	beta_u = 0.5*PI;
	screen("Enter optional limits on the incident angle, "
	       "between 0 and 90 degrees\n"
	       "\tEnter the lower limit [0 degrees]: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&beta_l);
	    if (beta_l < 0.0 || beta_l > 90.0)
		beta_l = 0.0;
	    beta_l = radians(beta_l);
	}
	screen("\tEnter the upper limit [90 degrees]: ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscan_float(s,&beta_u);
	    if (beta_u < 0.0 || beta_u > 90.0)
		beta_u = 90.0;
	    beta_u = radians(beta_u);
	}
	npts = NUM_POINTS;
	dbeta = (beta_u - beta_l)/npts;
	screen("The incident angle increment can be specified either "
	       "directly or indirectly\n\t"
	       "in terms of the number of intervals between the minimum "
	       "and maximum\n\t"
	       "values.\n"
	       "Enter choice (direct or indirect [indirect]: ");
	(void) Gets(s);
	if (s[0] == 'd' || s[0] == 'D')
	{
	    screen("Enter the incident angle increment "
	           "in degrees [%g]: ",degrees(dbeta));
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&dbeta);
	    	dbeta = radians(dbeta);
	    	npts = irint((beta_u - beta_l)/dbeta);
	    	dbeta = (beta_u - beta_l)/npts;
	    }
	}
	else
	{
	    screen("Enter the number of intervals to plot [%d]: ",npts);
	    (void) Gets(s);
	    if (s[0] != '\0')
		(void) sscanf(s,"%d",&npts);
	    dbeta = (beta_u - beta_l)/npts;
	}
	bi_array(&vars,npts+1,nvars,FLOAT);

	beta = beta_l;
	if (beta == 0.0)
	{
	    double     v[MAXD];
	    set_state(Tsl,VGAS_STATE,RP->state[1]);
	    zero_state_velocity(Tsl,dim);
	    for (i = 0; i < dim; ++i)
	    	v[i] = vel(i,RP->state[1]);
	    Vel(Tsl)[0] = scalar_product(v,nor,dim);
	    set_state(Tsr,VGAS_STATE,RP->state[6]);
	    zero_state_velocity(Tsr,dim);
	    for (i = 0; i < dim; ++i)
	    	v[i] = vel(i,RP->state[6]);
	    Vel(Tsr)[0] = scalar_product(v,nor,dim);
	    (void) find_mid_state(Tsl,Tsr,0.0,&pstarl,&pstarr,&ustarl,&ustarr,
				  &ml,&mr,&l_wave,&r_wave);
	    vars[0][INCIDENT_ANGLE]		= 0.0;
	    vars[0][PSTAR]			= 0.5*(pstarl+pstarr);
	    vars[0][SHEAR]                      =
	        Vel(RP->state[0])[0]-Vel(RP->state[6])[0];
	    vars[0][REFLECTED_ANGLE_1]		= 180.0;
	    vars[0][REFLECTED_ANGLE_2]		= 180.0;
	    vars[0][REFLECTED_ANGLE_3]		= 180.0;
	    vars[0][CONTACT_ANGLE]		= 180.0;
	    vars[0][TRANSMITTED_ANGLE]		= 180.0;
	    vars[0][INCIDENT_TURN_ANGLE]	= 0.0;
	    vars[0][REFLECTED_TURN_ANGLE]	= 0.0;
	    vars[0][TRANSMITTED_TURN_ANGLE]	= 0.0;
	    vars[0][THETA_ERR]			= 0.0;
	    beta += dbeta;
	    t[0][0] = cos(beta);
	    t[0][1] = sin(beta);
	    node_v[0] = v0/t[0][1];
	    if (is_regular_diffraction_node(NULL,node_v,NULL,t,RP,NULL,
				                &is_reflected_shock,NULL,
						DIFFRACTION_NODE,flag) 
						!= REGULAR_DIFFRACTION)
	    {
	    	screen("ERROR in "
	    	       "plot_shock_refraction_vs_incident_angle(), "
	    	       "is_regular_diffraction_node() failed for "
	    	       "small angle case\n");
	    	return;
	    }
	    vars[1][INCIDENT_ANGLE]			= degrees(beta);
	    vars[1][PSTAR] = pressure(RP->state[4]);
	    a[0] = cos(RP->ang[4]);
	    a[1] = sin(RP->ang[4]);
	    v4 = Vel(RP->state[4])[0]*a[0] + Vel(RP->state[4])[1]*a[1];
	    v5 = Vel(RP->state[5])[0]*a[0] + Vel(RP->state[5])[1]*a[1];
	    vars[1][SHEAR] = v5 - v4;
	    vars[1][REFLECTED_ANGLE_1] = degrees(RP->ang[1]);
	    vars[1][REFLECTED_ANGLE_2] = degrees(RP->ang[2]);
	    vars[1][REFLECTED_ANGLE_3] = degrees(RP->ang[3]);
	    vars[1][CONTACT_ANGLE] = degrees(RP->ang[4]);
	    vars[1][TRANSMITTED_ANGLE]		= RP->ang[5];
	    vars[1][INCIDENT_TURN_ANGLE] = degrees(RP->theta[0]);
	    vars[1][REFLECTED_TURN_ANGLE] = degrees(RP->theta[2]);
	    vars[1][TRANSMITTED_TURN_ANGLE]	= degrees(RP->theta[5]);
	    vars[1][THETA_ERR] = vars[1][INCIDENT_TURN_ANGLE] +
	    			 vars[1][REFLECTED_TURN_ANGLE] -
	    			 vars[1][TRANSMITTED_TURN_ANGLE];
	    vars[0][AHEAD_INCIDENT_MACH_NUMBER]	=
	    vars[1][AHEAD_INCIDENT_MACH_NUMBER]	= RP->M[0];
	    vars[0][BEHIND_INCIDENT_MACH_NUMBER]	=
	    vars[1][BEHIND_INCIDENT_MACH_NUMBER]	= RP->M[1];
	    vars[0][BEHIND_REFLECTED_MACH_NUMBER]	=
	    vars[1][BEHIND_REFLECTED_MACH_NUMBER]	= RP->M[4];
	    vars[0][BEHIND_TRANSMITTED_MACH_NUMBER]	=
	    vars[1][BEHIND_TRANSMITTED_MACH_NUMBER]	= RP->M[5];
	    vars[0][AHEAD_TRANSMITTED_MACH_NUMBER]	=
	    vars[1][AHEAD_TRANSMITTED_MACH_NUMBER]	= RP->M[6];
	    i = 2;
	    beta += dbeta;
	    nlines = 2;
	}
	else
	    nlines = 0;
	for (; nlines <= npts; ++nlines, beta += dbeta)
	{
	    t[0][0] = cos(beta);
	    t[0][1] = sin(beta);
	    node_v[0] = v0/t[0][1];
	    if (is_regular_diffraction_node(NULL,node_v,NULL,t,RP,NULL,
				            &is_reflected_shock,NULL,
					    DIFFRACTION_NODE,flag) 
					    != REGULAR_DIFFRACTION)
	    {
	    	break;
	    }
	    vars[nlines][INCIDENT_ANGLE] = degrees(beta);
	    vars[nlines][PSTAR] = pressure(RP->state[4]);
	    a[0] = cos(RP->ang[4]);
	    a[1] = sin(RP->ang[4]);
	    v4 = Vel(RP->state[4])[0]*a[0] + Vel(RP->state[4])[1]*a[1];
	    v5 = Vel(RP->state[5])[0]*a[0] + Vel(RP->state[5])[1]*a[1];
	    vars[nlines][SHEAR] = v5 - v4;
	    vars[nlines][REFLECTED_ANGLE_1] = degrees(RP->ang[1]);
	    vars[nlines][REFLECTED_ANGLE_2] = degrees(RP->ang[2]);
	    vars[nlines][REFLECTED_ANGLE_3] = degrees(RP->ang[3]);
	    vars[nlines][CONTACT_ANGLE] = degrees( RP->ang[4]);
	    vars[nlines][TRANSMITTED_ANGLE] = degrees(RP->ang[5]);
	    vars[nlines][INCIDENT_TURN_ANGLE] = degrees(RP->theta[0]);
	    vars[nlines][REFLECTED_TURN_ANGLE] = degrees(RP->theta[2]);
	    vars[nlines][TRANSMITTED_TURN_ANGLE] = degrees(RP->theta[5]);
	    vars[nlines][THETA_ERR] = vars[nlines][INCIDENT_TURN_ANGLE] +
	    			      vars[nlines][REFLECTED_TURN_ANGLE] -
	    			      vars[nlines][TRANSMITTED_TURN_ANGLE];
	    vars[nlines][AHEAD_INCIDENT_MACH_NUMBER]	 = RP->M[0];
	    vars[nlines][BEHIND_INCIDENT_MACH_NUMBER]	 = RP->M[1];
	    vars[nlines][BEHIND_REFLECTED_MACH_NUMBER]	 = RP->M[4];
	    vars[nlines][BEHIND_TRANSMITTED_MACH_NUMBER] = RP->M[5];
	    vars[nlines][AHEAD_TRANSMITTED_MACH_NUMBER]	 = RP->M[6];
	}
	width = strlen("BEHIND_TRANSMITTED_MACH_NUMBER");
	print_plot_header(header,width,NULL,NULL);
	(void) sprintf(ffmt,"%%-%dg ",(int)width);
	for (i = 0; i < nlines; ++i)
	{
	    for (j = 0; j < nvars; ++j)
	    	(void) printf(ffmt,vars[i][j]);
	    (void) printf("\n");
	}
	(void) printf("\n\n");

	free(vars);
}		/*end plot_shock_refraction_vs_incident_angle*/

LOCAL	int init_plot_shock_refraction(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	double		*theta0,
	Locstate	*state,
	int		dim,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	char		s[Gets_BUF_SIZE];
	double		p2, sin_beta;
	double		shock_ang;
	int		is_turn_ang_positive;
	Gas_param	*params;
	static double	node_v[3] = {0.0, 0.0, 0.0};
	static double	nor[] = {1.0, 0.0, 0.0};

	screen("Do you wish to solve a normal shock refraction [no]: ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    return normal_shock_refraction(init,ip,state);

	screen("Do you wish to plot the refraction parameters as a\n\t"
	       "function of incident angle [no]: ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    plot_shock_refraction_vs_incident_angle(init,ip,dim);
	    return 0;
	}

	params = init_eos_params(init,ip,
				 " ahead of the incident shock wave",NO);
	prompt_for_thermodynamics(state[0],params, " of the incident gas");

	prompt_for_behind_shock_state(state[0],state[2],NO,nor,
				      state_type(state[0]),YES,init);
	p2 = pressure(state[2]);
	init_shock_polar_ahead_velocity(state[0],&p2,&max_pressure_given[0],
	                                &pmax[0],"for the incident gas ");

	set_plot_control("incident ",&plot_control[0],&npts[0],state[0]);

	params = init_eos_params(init,ip,
				 " ahead of the transmitted shock wave",NO);
	init_state_across_contact(state[1],state[0],params,node_v,
				" of the gas ahead of the transmitted shock");
	set_plot_control("transmitted ",&plot_control[1],&npts[1],state[1]);

	screen("Enter p if the flow behind the incident \n"
	       "shock is turned counter clockwise relative to the \n"
	       "ahead flow, otherwise the flow will be assumed \n"
	       "to be turned in the clockwise direction.\n"
	       "\tEnter choice here: ");
	(void) Gets(s);
	is_turn_ang_positive = (s[0] == 'p' || s[0] == 'P') ? YES : NO;
	if (!s_polar_3(state[0],YES,p2,is_turn_ang_positive,NO,node_v,state[2],
		      &shock_ang,&theta0[2]))
	{
	    screen("ERROR in init_plot_shock_refraction(), "
	           "s_polar_3 failed\n");
	    clean_up(ERROR);
	}
	set_type_of_state(state[2],state_type(state[0]));
	set_state(state[0],VGAS_STATE,state[0]);
	set_state(state[1],VGAS_STATE,state[1]);
	set_state(state[2],VGAS_STATE,state[2]);
	set_plot_control("reflected ",&plot_control[2],&npts[2],state[2]);

	sin_beta = fabs(mass_flux(pressure(state[2]),state[0])/
	               (Dens(state[0])*mag_vector(Vel(state[0]),2)));
	if (fabs(sin_beta) < 1.0)
	{
	    NODE_FLAG      flag;
	    boolean is_reflected_shock;
	    static RP_DATA  *RP = NULL;
	    static double    **t;

	    if (RP == NULL)
	    {
	        RP = allocate_RP_DATA_structure(sizeof(VGas),NO,VGAS_STATE);
	        bi_array(&t,2,2,FLOAT);
	        t[1][0] = 1.0;
	        t[1][1] = 0.0;
	    }
	    clear_node_flag(flag);
	    use_subsonic_state(flag) = YES;
	    set_state(RP->state[0],VGAS_STATE,state[0]);
	    set_state(RP->state[1],VGAS_STATE,state[2]);
	    set_state(RP->state[6],VGAS_STATE,state[1]);
	    RP->ang_dir = (is_turn_ang_positive) ? COUNTER_CLOCK : CLOCKWISE;

	    t[0][0] = sqrt(1.0 - t[0][1]*t[0][1]);
	    t[0][1] = (RP->ang_dir == CLOCKWISE) ? -sin_beta : sin_beta;
	    if (is_regular_diffraction_node(NULL,node_v,NULL,t,RP,NULL,
	                                    &is_reflected_shock,NULL,
					    DIFFRACTION_NODE,flag)
					    == REGULAR_DIFFRACTION)
	    {
	        print_RP_DATA(RP,node_v);
	    }
	}

	return 3;
}		/*end init_plot_shock_refraction*/

LOCAL	int init_plot_total_reflection(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	double		*theta0,
	Locstate	*state,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	char		s[Gets_BUF_SIZE];
	double		p1;
	double		shock_ang;
	int		is_turn_ang_positive;
	Gas_param	*params;
	static double	node_v[MAXD]; /*Statics initialized to zero*/

	params = init_eos_params(init,ip,
				 " ahead of the incident shock wave",NO);
	prompt_for_thermodynamics(state[0],params, " of the incident gas");

	screen("Enter the pressure behind the incident shock wave: ");
	(void) Scanf("%f\n",&p1);
	init_shock_polar_ahead_velocity(state[0],
		    &p1,&max_pressure_given[0],&pmax[0],
		    "for the incident gas ");

	set_plot_control("incident ",&plot_control[0],&npts[0],state[0]);

	screen("Enter p if the flow behind the incident \n");
	screen("shock is turned counter clockwise relative to the \n");
	screen("ahead flow, otherwise the flow will be assumed \n");
	screen("to be turned in the clockwise direction.\n");
	screen("\tEnter choice here: ");
	(void) Gets(s);
	is_turn_ang_positive = (s[0] == 'p' || s[0] == 'P') ? YES : NO;
	Check_return(
	    s_polar_3(state[0],YES,p1,is_turn_ang_positive,NO,node_v,state[1],
		      &shock_ang,&theta0[1]),
	    init_plot_total_reflection) 
	set_type_of_state(state[1],state_type(state[0]));
	set_state(state[0],VGAS_STATE,state[0]);
	set_state(state[1],VGAS_STATE,state[1]);
	set_plot_control("reflected ",&plot_control[1],&npts[1],state[1]);
	return 2;
}		/*end init_plot_total_reflection*/

LOCAL	int init_plot_shock_overtake(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	double		*theta0,
	Locstate	*state,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	char		s[Gets_BUF_SIZE];
	double		p1, p2;
	double		shock_ang;
	int		is_turn_ang_positive;
	Gas_param	*params;
	static double	node_v[MAXD]; /*Statics initialized to zero*/

	params = init_eos_params(init,ip,"",NO);
	prompt_for_thermodynamics(state[0],params, " of the incident gas");
	screen("Enter the pressure behind the ");
	screen("shock being overtaken: ");
	(void) Scanf("%f\n",&p1);
	init_shock_polar_ahead_velocity(state[0],&p1,
		&max_pressure_given[0],&pmax[0],"for the incident gas ");

	set_plot_control("overtook ",&plot_control[0],&npts[0],state[0]);
	screen("Enter p if the flow behind the incident \n");
	screen("shock is turned counter clockwise relative to the \n");
	screen("ahead flow, otherwise the flow will be assumed \n");
	screen("to be turned in the clockwise direction.\n");
	screen("\tEnter choice here: ");
	(void) Gets(s);
	is_turn_ang_positive = (s[0] == 'p' || s[0] == 'P') ? YES : NO;
	Check_return(
	    s_polar_3(state[0],YES,p1,is_turn_ang_positive,NO,node_v,state[1],
		      &shock_ang,&theta0[1]),
	    init_plot_shock_overtake) 
	set_type_of_state(state[1],state_type(state[0]));
	set_plot_control("overtaking incident ",&plot_control[1],
			 &npts[1],state[1]);
	screen("Enter the pressure behind the ");
	screen("overtaking incident shock: ");
	(void) Scanf("%f\n",&p2);
	screen("Enter p if the flow behind the incident \n");
	screen("shock is turned counter clockwise relative to the \n");
	screen("ahead flow, otherwise the flow will be assumed \n");
	screen("to be turned in the clockwise direction.\n");
	screen("\tEnter choice here: ");
	(void) Gets(s);
	is_turn_ang_positive = (s[0] == 'p' || s[0] == 'P') ? YES : NO;
	Check_return(
	    s_polar_3(state[1],YES,p2,is_turn_ang_positive,NO,node_v,state[2],
		      &shock_ang,&theta0[2]),
	    init_plot_shock_overtake) 
	theta0[2] += theta0[1];
	set_type_of_state(state[2],state_type(state[1]));
	set_state(state[0],VGAS_STATE,state[0]);
	set_state(state[1],VGAS_STATE,state[1]);
	set_state(state[2],VGAS_STATE,state[2]);
	set_plot_control("reflected ",&plot_control[2],&npts[2],state[2]);
	return 3;
}		/*end init_plot_shock_overtake*/

LOCAL	int	compute_oned_rp(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Locstate	         left, right;
	Locstate	         left_mid, right_mid;
	double		         pjump;
	double		         pmidl, pmidr;
	double		         umidl, umidr;
	double		         ml, mr;
	double		         max_speed, speed[17];
	double                    posn;
	double		         ul, cl, uml, cml, umr, cmr, ur, cr;
	double		         x, dx, dt, cfl;
	int		         i;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	char		         s[Gets_BUF_SIZE];

	ip->root->front->rect_grid->dim = 1;
	g_alloc_state(&left,sizest);
	g_alloc_state(&left_mid,sizest);
	g_alloc_state(&right_mid,sizest);
	g_alloc_state(&right,sizest);

	(void) prompt_for_eos_params(init,ip,NO,"");

	Init_params(left,prompt_for_eos_params(init,ip,NO,
					       " for the left state"));
	prompt_for_ref_state(" on the left",left,TGAS_STATE,Params(left),init);

	Init_params(right,prompt_for_eos_params(init,ip,NO,
						" for the right state"));
	prompt_for_ref_state(" on the right",right,TGAS_STATE,
			     Params(right),init);

	pjump = prompt_for_rp_surface_tension(left,right);


	(void) find_mid_state(left,right,pjump,&pmidl,&pmidr,&umidl,&umidr,&ml,
			      &mr,&l_wave,&r_wave);

	(void) printf("umidl - (Vel(left)[0] - "
	              "riemann_wave_curve(left,pmidl)) = %18.16g\n",
		      umidl-(Vel(left)[0]-riemann_wave_curve(left,pmidl)));
	(void) printf("umidr - (Vel(right)[0] + "
	              "riemann_wave_curve(right,pmidr)) = %18.16g\n",
		       umidr-(Vel(right)[0]+riemann_wave_curve(right,pmidr)));
	
	Dens(left_mid) = (l_wave == SHOCK) ?
				dens_Hugoniot(pmidl,left) :
				dens_rarefaction(pmidl,left);
	Vel(left_mid)[0] = umidl;
	Press(left_mid) = pmidl;
	Set_params(left_mid,left);
	set_type_of_state(left_mid,TGAS_STATE);
	reset_gamma(left_mid);

	Dens(right_mid) = (r_wave == SHOCK) ?
				dens_Hugoniot(pmidr,right) :
				dens_rarefaction(pmidr,right);
	Vel(right_mid)[0] = umidr;
	Press(right_mid) = pmidr;
	Set_params(right_mid,right);
	set_type_of_state(right_mid,TGAS_STATE);
	reset_gamma(right_mid);

	ul = vel(0,left);		cl = sound_speed(left);
	uml = vel(0,left_mid);		cml = sound_speed(left_mid);
	umr = vel(0,right_mid);		cmr = sound_speed(right_mid);
	ur = vel(0,right);		cr = sound_speed(right);

	speed[0] = ul - cl; speed[1] = ul; speed[2] = ul + cl;
	if (l_wave == SHOCK)
	{
	    speed[3] = speed[4] = ul - ml/Dens(left);
	}
	else
	{
	    speed[3] = speed[0];
	    speed[4] = uml - cml;
	}
	speed[5] = uml - cml; speed[6] = uml; speed[7] = uml + cml;
	speed[8] = 0.5*(uml+umr);
	speed[9] = umr - cmr; speed[10] = umr; speed[11] = umr + cmr;
	if (r_wave == SHOCK)
	{
	    speed[12] = speed[13] = ur + mr/Dens(right);
	}
	else
	{
	    speed[12] = speed[11];
	    speed[13] = ur + cr;
	}
	speed[14] = ur - cr; speed[15] = ur; speed[16] = ur + cr;

	verbose_print_state("left",left);
	verbose_print_state("left_mid",left_mid);
	verbose_print_state("right_mid",right_mid);
	verbose_print_state("right",right);

	(void) printf("ml = %18.16g, mr = %18.16g\n",ml,mr);
	(void) printf("l_wave = %s, r_wave = %s\n",
		(l_wave == SHOCK) ? "SHOCK" : "RAREFACTION",
		(r_wave == SHOCK) ? "SHOCK" : "RAREFACTION");

	if (l_wave == SHOCK)
	    (void) printf("left shock speed = %18.16g\n",speed[3]);
	else
	{
	    (void) printf("left leading edge rarefaction speed = %18.16g\n",
			  speed[3]);
	    (void) printf("left trailing edge rarefaction speed = %18.16g\n",
			  speed[4]);
	}
	(void) printf("contact speed = %18.16g\n",speed[8]);
	if (r_wave == SHOCK)
	    (void) printf("right shock speed = %18.16g\n",speed[13]);
	else
	{
	   (void) printf("right trailing edge rarefaction speed = %18.16g\n",
			 speed[12]);
	   (void) printf("right leading edge rarefaction speed = %18.16g\n",
			 speed[13]);
	}
	max_speed = 0.0;
	for (i = 0; i < 17; ++i)
	    max_speed = max(max_speed,fabs(speed[i]));
	(void) printf("maximum wave speed = %18.16g\n",max_speed);

	x = 0.0;
	screen("Enter a coordinate position (dflt = %18.16g): ",x);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&x);

	dx = 1.0; /*TOLERANCE*/
	screen("Enter a grid spacing distance (dflt = %18.16g): ",dx);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&dx);

	cfl = 0.75;
	screen("Enter a CFL factor (dflt = %18.16g): ",cfl);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&cfl);

	dt = cfl * dx / max_speed;
	screen("Enter a positive time (dflt = %18.16g): ",dt);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&dt);

	(void) printf("Position of waves at time %18.16g\n",dt);
	if (l_wave == SHOCK)
	{
	    posn = x + dt*speed[3];
	    (void) printf("left shock position = %18.16g\n",posn);
	    if (x != 0.0)
	        (void) printf("              ratio = %18.16g\n",posn/x);
	}
	else
	{
	    posn = x + dt*speed[3];
	    (void) printf("left leading edge rarefaction position = %18.16g\n",
			  posn);
	    if (x != 0.0)
	        (void) printf("                                 "
			      "ratio = %18.16g\n",posn/x);
	    posn = x + dt*speed[4];
	    (void) printf("left trailing edge rarefaction position = %18.16g\n",
			  posn);
	    if (x != 0.0)
	        (void) printf("                                  "
			      "ratio = %18.16g\n",posn/x);
	}
	posn = x + dt*speed[8];
	(void) printf("contact position = %18.16g\n",posn);
	if (x != 0.0)
	    (void) printf("           ratio = %18.16g\n",posn/x);
	if (r_wave == SHOCK)
	{
	    posn = x + dt*speed[13];
	    (void) printf("right shock position = %18.16g\n",posn);
	    if (x != 0.0)
	        (void) printf("               ratio = %18.16g\n",posn/x);
	}
	else
	{
	    posn = x + dt*speed[12];
	    (void) printf("right trailing edge rarefaction position = "
			  "%18.16g\n",posn);
	    if (x != 0.0)
	        (void) printf("                                   "
			      "ratio = %18.16g\n",posn/x);
	    posn = x + dt*speed[13];
	    (void) printf("right leading edge rarefaction position = %18.16g\n",
			  posn);
	    if (x != 0.0)
	        (void) printf("                                  "
			      "ratio = %18.16g\n",posn/x);
	}
	clean_up(0);
	return 0;
}		/*end compute_oned_rp*/

/*ARGSUSED*/
LOCAL	int	compute_refracted_vorticity(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Gas_param      *params;
	double          beta, p1, m;
	double          nor[2], node_v[2];
	double          w, w0, p0, dw, w_max, w_min;
	double          win_bar, wout_bar, A, B, residual;
	double          v4, v5, v6, a[2];
	boolean           is_reflected_shock;
	int            n, Npts;
	const int      num_points = 100;
	char	       s[Gets_BUF_SIZE];
	NODE_FLAG      flag;
	static RP_DATA *RP = NULL;
	static double   *win, *wout, *pout;
	static double   **t;
	const char *sd = "(v ahead incident shock - v ahead transmitted shock)";

	if (RP == NULL)
	{
	    RP = allocate_RP_DATA_structure(sizeof(VGas),NO,VGAS_STATE);
	    RP->ang_dir = COUNTER_CLOCK;
	    bi_array(&t,2,2,FLOAT);
	    uni_array(&win,num_points+2,FLOAT);
	    uni_array(&wout,num_points+2,FLOAT);
	    uni_array(&pout,num_points+2,FLOAT);
	    t[1][0] = 1.0;
	    t[1][1] = 0.0;
	}
	clear_node_flag(flag);
	use_subsonic_state(flag) = YES;

	(void) prompt_for_eos_params(init,ip,NO,"");
	params = prompt_for_eos_params(init,ip,NO,"\n\tahead of the incident "
						  "shock wave");
	prompt_for_thermodynamics(RP->state[0],params," of the incident gas");
	set_state(RP->state[0],VGAS_STATE,RP->state[0]);
	screen("Enter the incident angle in degrees: ");
	(void) Scanf("%f\n",&beta);
	beta = fabs(radians(beta));
	t[0][0] = cos(beta);    t[0][1] = sin(beta);
	nor[0] = sin(beta);	nor[1] = -cos(beta);
	prompt_for_behind_shock_state(RP->state[0],RP->state[1],NO,nor,
				      state_type(RP->state[0]),YES,init);
	p1 = pressure(RP->state[1]);
	m = mass_flux(p1,RP->state[0]);
	node_v[0] = m/(Dens(RP->state[0])*nor[0]);
	node_v[1] = 0.0;
	params = prompt_for_eos_params(init,ip,NO,"\n\tahead of the "
						     "transmitted shock wave");
	init_state_across_contact(RP->state[6],RP->state[0],params,node_v,
				" of the gas ahead of the transmitted shock");
	set_state(RP->state[6],VGAS_STATE,RP->state[6]);

	w_max = 100.0*sound_speed(RP->state[0]);
	screen("Enter the maximum vorticity to be plotted (dflt = %g): ",
	       w_max);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&w_max);

	w_min = min_vorticity_for_regular_refraction(RP,t,node_v,flag);
	v6 = Vel(RP->state[6])[0];
	Vel(RP->state[6])[0] = Vel(RP->state[0])[0];
	if (is_regular_diffraction_node(NULL,node_v,NULL,t,RP,NULL,
	                                &is_reflected_shock,NULL,
					DIFFRACTION_NODE,flag)
					== REGULAR_DIFFRACTION)
	{
	    a[0] = cos(RP->ang[4]);
	    a[1] = sin(RP->ang[4]);
	    v4 = Vel(RP->state[4])[0]*a[0] + Vel(RP->state[4])[1]*a[1];
	    v5 = Vel(RP->state[5])[0]*a[0] + Vel(RP->state[5])[1]*a[1];
	    w0 = v5 - v4;
	    p0 = 0.5*(pressure(RP->state[4])+pressure(RP->state[5]));
	    (void) printf("\nIncoming shear %s = 0.0, "
	                  "Outgoing shear %s = %g\n",sd,sd,w0);
	    print_RP_DATA(RP,node_v);
	}
	dw = (w_max - w_min)/num_points;
	for (w=w_min, n=0; (n <= num_points+1) && (w <= w_max); ++n, w += dw)
	{
	    win[n] = w - v6;
	    Vel(RP->state[6])[0] = -win[n];
	    if (is_regular_diffraction_node(NULL,node_v,NULL,t,RP,NULL,
	                                    &is_reflected_shock,NULL,
					    DIFFRACTION_NODE,flag)
					    == REGULAR_DIFFRACTION)
	    {
	        a[0] = cos(RP->ang[4]);
	        a[1] = sin(RP->ang[4]);
	        v4 = Vel(RP->state[4])[0]*a[0] + Vel(RP->state[4])[1]*a[1];
	        v5 = Vel(RP->state[5])[0]*a[0] + Vel(RP->state[5])[1]*a[1];
		wout[n] = v5 - v4;
		pout[n] = 0.5*(pressure(RP->state[4])+pressure(RP->state[5]));
	        (void) printf("\nIncoming shear %s = %g, "
		              "Outgoing shear %s = %g\n",
			      sd,win[n],sd,wout[n]);
		print_RP_DATA(RP,node_v);
	    }
	    if ((w < 0.0) && (w+dw > 0.0))
	    {
	        win[++n] = -v6;
		wout[n] = w0;
		pout[n] = p0;
	    }
	}
	Npts = n;
	Vel(RP->state[6])[0] = v6;


	win_bar = wout_bar = 0.0;
	for (n = 0; n < Npts; ++n)
	{
	    win_bar += win[n];
	    wout_bar += wout[n];
	}
	win_bar /= Npts;
	wout_bar /= Npts;
	A = B = 0.0;
	for (n = 0; n < Npts; ++n)
	{
	    B += (win[n]-win_bar)*(wout[n]-wout_bar);
	    A += (win[n]-win_bar)*(win[n]-win_bar);
	}
	B /= A;
	A = wout_bar - B*win_bar;
	residual = 0.0;
	for (n = 0; n < Npts; ++n)
	    residual += (wout[n]-A-B*win[n])*(wout[n]-A-B*win[n]);
	residual = sqrt(residual/n);
	(void) printf("Least squares fit\n");
	(void) printf("win = %g + %g * win, residual = %g\n",A,B,residual);
	(void) printf("vorticity for zero initial shear = %g\n",w0);

	(void) output();
	(void) printf("%-18s %-18s %-18s %-18s %-18s\n",
	              "INCOMING_VORTICITY","OUTGOING_VORTICITY",
		      "DELTA_VORTICITY","LINEAR_RESIDUAL","PRESSURE");
	for (n = 0; n < Npts; ++n)
	    (void) printf("%-18g %-18g %-18g %-18g %-18g\n",
	                  win[n],wout[n],wout[n]-win[n],wout[n]-A-B*win[n],
			  pout[n]);
	(void) printf("\n");

	clean_up(0);
	return 0;
}		/*end compute_refracted_vorticity*/

LOCAL	double	min_vorticity_for_regular_refraction(
	RP_DATA   *RP,
	double     **t,
	double     *node_v,
	NODE_FLAG flag)
{
	double w, w_regular, w_irregular;
	double c6;
	double v0, v6;
	boolean  is_reflected_shock;
	int   n;
	const int N = 10;

	v0 = Vel(RP->state[0])[0];
	v6 = Vel(RP->state[6])[0];
	if ((v0 - node_v[0])*(v6 - node_v[0]) <= 0.0)
	    return 0.0;
	if (is_regular_diffraction_node(NULL,node_v,NULL,t,RP,NULL,
	                                 &is_reflected_shock,NULL,
					 DIFFRACTION_NODE,flag)
					 != REGULAR_DIFFRACTION)
	{
	    Vel(RP->state[6])[0] = v6;
	    return 0.0;
	}
	w_regular = 0.0;

	c6 = sound_speed(RP->state[6]);
	w_irregular = c6 + (v6 - node_v[0]);
	if (w_irregular >= 0.0)
	    return 0.0;

	for (n = 0; n < N; ++n)
	{
	    w = 0.5*(w_regular+w_irregular);
	    Vel(RP->state[6])[0] = v6 - w;
	    if (is_regular_diffraction_node(NULL,node_v,NULL,t,RP,NULL,
	                                    &is_reflected_shock,NULL,
					    DIFFRACTION_NODE,flag)
					    == REGULAR_DIFFRACTION)
	        w_regular = w;
	    else
	        w_irregular = w;
	}
	Vel(RP->state[6])[0] = v6;
	return w_regular;
}		/*end min_vorticity_for_regular_refraction*/

LOCAL	int	npt_w_speed_test(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front		*front = ip->root->front;
	Grid		*grid = ip->root->grid;
	Locstate	sl, sr;
	Locstate	ansl, ansr;
	char		s[Gets_BUF_SIZE];
	double		*n, *p;
	double		V[3];
	int		i, j, dim = 1;
	static WSSten   *sten = NULL;

	screen("Enter the spatial dimension (default = %d): ",dim);
	(void) Gets(s);
	if (s[0] != 0)
	    (void) sscanf(s,"%d",&dim);
	if ((dim <= 0) || (dim > MAXD))
	{
	    screen("ERROR in npt_w_speed_test(), "
		   "Invalid dimension\n");
	    clean_up(ERROR);
	    return ERROR;
	}
	front->rect_grid->dim = dim;
	set_interface_hooks(dim,init);
	init_remap_and_rect_grid(grid->rect_grid,ip);

	g_init_physics(init,ip);
	sizest = front->sizest;

	if (sten == NULL)
	    sten = AllocDefaultWSSten(front);
	else
	    ClearWSStenData(sten);

	sl = sten->sl[0];
	sr = sten->sr[0];
	sten->front = front;
	sten->wave = ip->root->wave;
	sten->ncomp = FIRST_DYNAMIC_COMPONENT;
	sten->pcomp = FIRST_DYNAMIC_COMPONENT+1;
	sten->hs = NULL;
	n = sten->nor;
	p = sten->coords;
	g_alloc_state(&ansl,sizest);
	g_alloc_state(&ansr,sizest);

	n[0] = 1.0;	n[1] = 0.0;	n[2] = 0.0;
	p[0] = 0.0;	p[1] = 0.0;	p[2] = 0.0;
	front->interf = make_interface(dim);

	if (dim > 1)
	{
	    screen("Enter the wave normal vector (default = ");
	    if (dim == 2)
		screen("<%g, %g>",n[0],n[1]);
	    else if (dim == 3)
		screen("<%g, %g, %g>",n[0],n[1],n[2]);
	    screen("): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	if (dim == 2)
	    	    (void) Scanf("%f %f\n",n,n+1);
	    	else if (dim == 3)
	    	    (void) Scanf("%f %f %f\n",n,n+1,n+2);
	    }
	}
	screen("Enter the position of the front (default = ");
	if (dim == 2)
	    screen("(%g, %g)",p[0],p[1]);
	else if (dim == 3)
	    screen("(%g, %g, %g(",p[0],p[1],p[2]);
	screen("): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if (dim == 2)
	    	(void) Scanf("%f %f\n",p,p+1);
	    else if (dim == 3)
	    	(void) Scanf("%f %f %f\n",p,p+1,p+2);
	}

	sten->dn = grid_size_in_direction(n,front->rect_grid->h,dim);
	for (i = 0; i < sten->nsts; ++i)
	{
	    for (j = 0; j < dim; ++j)
	    {
	        sten->lcrds[i][j] = p[j] - i*sten->dn*n[j];
	        sten->rcrds[i][j] = p[j] + i*sten->dn*n[j];
	    }
	}
	screen("Enter the time step: ");
	(void) Scanf("%f\n",&sten->dt);

	(void) prompt_for_eos_params(init,ip,NO,"");

	Init_params(sl,prompt_for_eos_params(init,ip,NO," for the left side"));
	for (i = sten->nsts-1; i >= 0; --i)
	{
	    (void) sprintf(s," sl[%d]",i);
	    prompt_for_ref_state(s,sten->sl[i],GAS_STATE,Params(sl),init);
	}

	Init_params(sten->sr[0],
		    prompt_for_eos_params(init,ip,NO," for the right side"));
	for (i = 0; i < sten->nsts; ++i)
	{
	    (void) sprintf(s," sr[%d]",i);
	    prompt_for_ref_state(s,sten->sr[i],GAS_STATE,Params(sten->sr[0]),
				 init);
	}

	sten->w_type = prompt_for_wave_type("",front->interf,ip);
	if (sten->w_type == CONTACT)
	    sten->pjump = prompt_for_rp_surface_tension(sl,sr);

	npt_w_speed(sten,ansl,ansr,V);

	for (i = sten->nsts-1; i >= 0; --i)
	{
	    (void) sprintf(s," sl[%d]",i);
	    verbose_print_state(s,sten->sl[i]);
	}
	for (i = 0; i < sten->nsts; ++i)
	{
	    (void) sprintf(s," sr[%d]",i);
	    verbose_print_state(s,sten->sr[i]);
	}
	verbose_print_state("ansl",ansl);
	verbose_print_state("ansr",ansr);

	clean_up(0);
	return 0;
}		/*end npt_w_speed_test*/

LOCAL	double	prompt_for_rp_surface_tension(
	Locstate	left,
	Locstate	right)
{
	char		s[Gets_BUF_SIZE];
	double pjump = 0.0;

	screen("Do you want to use surface tension (default = no): ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
		double		tension = 0.0, r = 0.0, M = 0.0;

		screen("Enter the coefficient of ");
		screen("surface tension (default = %g): ",tension);
		(void) Gets(s);
		if (s[0] != '\0') (void) sscan_float(s,&tension);
		screen("Enter the optional coefficient of ");
		screen("dynamic surface tension (default = %g): ",M);
		(void) Gets(s);
		if (s[0] != '\0')
			(void) sscan_float(s,&M);
		if (M != 0.0)
		{
		       double		U = 0.0, dh = 1.0;

		       screen("Enter the velocity shear across ");
		       screen("the interface (default = %g): ",U);
		       (void) Gets(s);
		       if (s[0] != '\0') (void) sscan_float(s,&U);
		       screen("Enter the lenght scale dh (default = %g): ",dh);
		       (void) Gets(s);
		       if (s[0] != '\0') (void) sscan_float(s,&dh);

		       tension += M*Dens(left)*Dens(right)*sqr(U)*dh / 
			       (2.0*PI*(Dens(left)+Dens(right)));
		}

		screen("Enter the radius of mean curvature (default = %g): ",r);
		(void) Gets(s);
		if (s[0] != '\0') (void) sscan_float(s,&r);
		pjump = -tension*r;
	}
	return pjump;
}		/*end prompt_for_rp_surface_tension*/

LOCAL	int init_plot_shock_crossing(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	double		*theta0,
	Locstate	*state,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	char		s[Gets_BUF_SIZE];
	double		p1, p2;
	double		shock_ang;
	int		is_turn_ang_positive;
	Gas_param	*params;
	static double	node_v[MAXD]; /*Statics initialized to zero*/

	params = init_eos_params(init,ip,"",NO);
	prompt_for_thermodynamics(state[0],params, " of the incident gas");
	init_shock_polar_ahead_velocity(state[0],
		(double *)NULL,&max_pressure_given[0],&pmax[0],
		"for the incident gas ");

	set_plot_control("first incident ",&plot_control[0],&npts[0],state[0]);

	screen("Enter the pressure behind the ");
	screen("first incident shock: ");
	(void) Scanf("%f\n",&p1);
	screen("Enter p if the flow behind the incident \n");
	screen("shock is turned counter clockwise relative to the \n");
	screen("ahead flow, otherwise the flow will be assumed \n");
	screen("to be turned in the clockwise direction.\n");
	screen("\tEnter choice here: ");
	(void) Gets(s);
	is_turn_ang_positive = (s[0] == 'p' || s[0] == 'P') ? YES : NO;
	Check_return(
	    s_polar_3(state[0],YES,p1,is_turn_ang_positive,NO,node_v,state[1],
		      &shock_ang,&theta0[1]),
	    init_plot_shock_crossing) 
	set_type_of_state(state[1],state_type(state[0]));
	set_plot_control("first reflected ",&plot_control[1],&npts[1],state[1]);

	screen("Enter the pressure behind the ");
	screen("second incident shock: ");
	(void) Scanf("%f\n",&p2);
	screen("Enter p if the flow behind the incident \n");
	screen("shock is turned counter clockwise relative to the \n");
	screen("ahead flow, otherwise the flow will be assumed \n");
	screen("to be turned in the clockwise direction.\n");
	screen("\tEnter choice here: ");
	(void) Gets(s);
	is_turn_ang_positive = (s[0] == 'p' || s[0] == 'P') ? YES : NO;
	Check_return(
	    s_polar_3(state[0],YES,p2,is_turn_ang_positive,NO,node_v,state[2],
		      &shock_ang,&theta0[2]),
	    init_plot_shock_crossing) 
	set_type_of_state(state[2],state_type(state[0]));
	set_state(state[0],VGAS_STATE,state[0]);
	set_state(state[1],VGAS_STATE,state[1]);
	set_state(state[2],VGAS_STATE,state[2]);
	set_plot_control("second reflected ",&plot_control[2],&npts[2],state[2]);
	return 3;
}		/*end init_plot_shock_crossing*/

LOCAL	int init_plot_sonic_locus(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	Locstate	*state,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	char		s[Gets_BUF_SIZE];
	Gas_param	*params;

	params = init_eos_params(init,ip,"",NO);
	prompt_for_thermodynamics(state[0],params, " of the incident gas");
	init_shock_polar_ahead_velocity(state[0],
		(double *)NULL,&max_pressure_given[0],&pmax[0],
		"for state0 ");
	set_plot_control("first incident ",&plot_control[0],&npts[0],state[0]);
	copy_state(state[1],state[0]);
	init_shock_polar_ahead_velocity(state[1],
		(double *)NULL,&max_pressure_given[1],&pmax[1],"for state1 ");
	set_plot_control("second incident ",&plot_control[1],&npts[1],state[1]);
	plot_control[2] = PLOT_SONIC_LOCUS | PLOT_CPLUS_WAVE;
	screen("Enter the number of points on the sonic locus ");
	screen("[%d]: ",npts[2]);
	(void) Gets(s);
	if (s[0] != '\0')
	{
		(void) sscanf(s,"%d",&npts[2]);
		npts[2] = min(npts[2],MAX_NUM_POINTS);
	}
	return 3;
}		/*end init_plot_sonic_locus*/

LOCAL	int init_plot_shock_polars(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	double		*theta0,
	Locstate	*state,
	int		*npts,
	int		*plot_control,
	double		*pmax,
	int		*max_pressure_given)
{
	char		s[Gets_BUF_SIZE];
	int		i, num_polars;
	double		v0x, v0y, vx, vy;
	Gas_param	*params;

	screen("Enter the number of polars to be plotted (<3) [1]: ");
	(void) Gets(s);
	if (s[0] != 0) (void) sscanf(s,"%d",&num_polars);
	else
		num_polars = 1;
	for (i = 0; i < num_polars; ++i)
	{
		(void) sprintf(s," for shock polar %d",i);
		params = init_eos_params(init,ip,"",NO);
		prompt_for_thermodynamics(state[i], params,s);
		init_shock_polar_ahead_velocity(state[i],
			    (double *)NULL,&max_pressure_given[i],
			   &pmax[i],(char *) NULL);
		set_state(state[i],VGAS_STATE,state[i]);
		set_plot_control("",&plot_control[i],&npts[i],state[i]);
	}
	v0x = vel(0,state[0]);
	v0y = vel(1,state[0]);
	for (i = 1; i < num_polars; ++i)
	{
		vx = vel(0,state[i]);
		vy = vel(1,state[i]);
		theta0[i] = atan2(v0x*vy-v0y*vx,v0x*vx+v0y*vy);
	}
	return num_polars;
}		/*end init_plot_shock_polars*/

LOCAL	void init_shock_polar_ahead_velocity(
	Locstate	state,
	double		*p2,
	int		*max_pressure_given,
	double		*pmax,
	const char	*mesg)
{
	double		M0, c0, p0, V0, Vmax;
	double		theta;
	double		inc_ang;
	char		s[Gets_BUF_SIZE];
	static	char	lmesg = '\0';

	if (mesg == NULL)
	    mesg = &lmesg;

	screen("The velocity data %smay be entered in one of\n"
	       "the following methods.\n"
	       "\tVelocity of the ahead state (V)\n"
	       "\tMach number of the ahead state (M)\n",mesg);
	if (p2 != NULL)
	    screen("\tIncident angle of shock (I)\n");
	if (max_pressure_given != NULL)
	    screen("\tPressure behind normal shock ");
	screen("(maximum pressure) (P)\n"
	       "Enter choice: ");
	(void) Gets(s);
	switch (s[0])
	{
	case 'V':
	case 'v':
	    screen("Enter x velocity and y velocity of the incident gas: ");
	    (void) Scanf("%f %f\n",&Vel(state)[0],&Vel(state)[1]);
	    screen("Enter an optional velocity to translate the "
	           "above velocities into\n\tthe steady state frame: ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        double	nvx, nvy;
	    	(void) sscanf(s,"%lf %lf",&nvx,&nvy);
	    	Vel(state)[0] -= nvx;
	    	Vel(state)[1] -= nvy;
	    }
	    break;
	case 'M':
	case 'm':
	    screen("Enter the Mach number of the ahead state: ");
	    (void) Scanf("%f\n",&M0);
	    c0 = sound_speed(state);
	    screen("Enter an optional turn velocity [0]: ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&theta);
	    	Vel(state)[0] = -c0*M0*cos(theta);
	    	Vel(state)[1] = -c0*M0*sin(theta);
	    }
	    else
	    {
	    	Vel(state)[0] = -c0*M0;
	    	Vel(state)[1] = 0.0;
	    }
	    break;
	case 'I':
	case 'i':
	    if (p2 == NULL)
	    {
	    	screen("Invalid choice\n");
	    	clean_up(ERROR);
	    }
	    screen("Enter the angle (in degrees) between the\n\t");
	    screen("incident shock and incoming streamline: ");
	    (void) Scanf("%f\n",&inc_ang);
	    inc_ang = radians(inc_ang);

	    Vel(state)[0] = -mass_flux(*p2,state)/
				fabs(Dens(state) * sin(inc_ang));
	    Vel(state)[1] = 0.0;
	    break;
	case 'P':
	case 'p':
	    if (max_pressure_given == NULL)
	    {
	    	screen("Invalid choice\n");
	    	clean_up(ERROR);
	    }
	    *max_pressure_given = YES;
	    screen("Enter the maximum pressure on the shock polar: ");
	    (void) Scanf("%f\n",pmax);
	    Vmax =1.0/dens_Hugoniot(*pmax,state);
	    p0 = pressure(state);
	    V0 = 1.0/Dens(state);
	    Vel(state)[0] = V0*sqrt((*pmax - p0)/(V0 - Vmax));
	    Vel(state)[1] = 0.0;
	    break;
	default:
	    screen("NO SUCH CHOICE\n");
	    clean_up(ERROR);
	}
	if (state_type(state) == GAS_STATE)
	{
	    Mom(state)[0] = Dens(state)*Vel(state)[0];
	    Mom(state)[1] = Dens(state)*Vel(state)[1];
	    Energy(state) += kinetic_energy(state);
            reset_gamma(state);
		    
	}
	if (debugging("sonic_point"))
	{
	    static	double node_v[MAXD]; /*Statics initialized to zero*/
	    double s_ang, t_ang;
	    double *av = NULL, *rv = NULL;
	    double M0sq;
	    static Locstate sonic_state = NULL;

	    if (sonic_state == NULL)
	    	g_alloc_state(&sonic_state,sizest);

	    M0sq = mach_number_squared(state,av,rv);
	    Check_return(
		s_polar_3(state,YES,pressure_at_sonic_point(M0sq,state),
			  NO,NO,node_v,sonic_state,&s_ang,&t_ang),
		init_shock_polar_ahead_velocity)
	    verbose_print_state("Sonic Point State",sonic_state);
	    (void) fprintf(stderr,"Sonic point: density = %g, pressure = %g, ",
			   Dens(sonic_state),pressure(sonic_state));
	    (void) fprintf(stderr,"internal energy = %g\n",
			   specific_internal_energy(sonic_state));
	    screen("Shock angle = %g, turning angle = %g\n\n",s_ang,t_ang);
	}
}		/*end init_shock_polar_ahead_velocity*/

LOCAL	void print_plot_header(
	const char	**header,
	size_t		width,
	double		*var_l,
	double		*var_u)
{
	char		cfmt[20], ffmt[20];
	double		*vl, *vu;
	const char	**s;

	(void) sprintf(cfmt,"%%-%ds ",(int)width);
	(void) sprintf(ffmt,"%%-%dg ",(int)width);
	(void) output();
	for (s = header; *s != NULL; ++s)
	    (void) printf(cfmt,*s);
	(void) printf("\n");
	if (var_l != NULL && var_u != NULL)
	{
	    for (vl = var_l, s = header; *s != NULL; ++vl, ++s)
	    	(void) printf(ffmt,*vl);
	    (void) printf("\n");
	    for (vu = var_u, s = header; *s != NULL; ++vu, ++s)
	    	(void) printf(ffmt,*vu);
	    (void) printf("\n");
	}
}		/*end print_plot_header*/

LOCAL	void print_line(
	int		j,
	int		npts,
	int		num_var,
	int		width,
	double		**var,
	double		*var_l,
	double		*var_u)
{
	char		ffmt[20];
	double		dvar, varmid;
	double		theta_ratio;
	double		theta_l = var_l[TURN_ANG];
	double		theta_u = var_u[TURN_ANG];
	int		k;

	(void) sprintf(ffmt,"%%-%dg ",width);
	if (var[j][TURN_ANG] < theta_l)
	{
		if ((j > 0) && (var[j-1][TURN_ANG] > theta_l))
		{
			(void) printf(ffmt,theta_l);
			theta_ratio = (theta_l - var[j-1][TURN_ANG]) /
					(var[j][TURN_ANG] - var[j-1][TURN_ANG]);
			for (k = 1; k < num_var; ++k)
			{
				dvar = var[j][k] - var[j-1][k];
				varmid = var[j-1][k] + theta_ratio*dvar;
				(void) printf(ffmt,varmid);
			}
			(void) printf("\n");
		}
		(void) printf(ffmt,theta_l);
		for (k = 1; k < num_var; ++k)
			(void) printf(ffmt,var[j][k]);
		(void) printf("\n");
		if (j < npts-1 && var[j+1][TURN_ANG] > theta_l)
		{
			(void) printf(ffmt,theta_l);
			theta_ratio = (theta_l - var[j+1][TURN_ANG]) /
					(var[j+1][TURN_ANG] - var[j][TURN_ANG]);
			for (k = 1; k < num_var; ++k)
			{
				dvar = var[j+1][k] - var[j][k];
				varmid = var[j+1][k] + theta_ratio*dvar;
				(void) printf(ffmt,varmid);
			}
			(void) printf("\n");
		}
	}
	else if (var[j][TURN_ANG] > theta_u)
	{
		if (j > 0 && var[j-1][TURN_ANG] < theta_u)
		{
			(void) printf(ffmt,theta_u);
			theta_ratio = (theta_u - var[j-1][TURN_ANG]) /
					(var[j][TURN_ANG] - var[j-1][TURN_ANG]);
			for (k = 1; k < num_var; ++k)
			{
				dvar = var[j][k] - var[j-1][k];
				varmid = var[j-1][k] + theta_ratio*dvar;
				(void) printf(ffmt,varmid);
			}
			(void) printf("\n");
		}
		(void) printf(ffmt,theta_u);
		for (k = 1; k < num_var; ++k)
			(void) printf(ffmt,var[j][k]);
		(void) printf("\n");
		if (j < npts-1 && var[j+1][TURN_ANG] < theta_u)
		{
			(void) printf(ffmt,theta_u);
			theta_ratio = (theta_u - var[j+1][TURN_ANG]) /
					(var[j+1][TURN_ANG] - var[j][TURN_ANG]);
			for (k = 1; k < num_var; ++k)
			{
				dvar = var[j+1][k] - var[j][k];
				varmid = var[j+1][k] + theta_ratio*dvar;
				(void) printf(ffmt,varmid);
			}
			(void) printf("\n");
		}
	}
	else
	{
		for (k = 0; k < num_var; ++k)
			(void) printf(ffmt,var[j][k]);
		(void) printf("\n");
	}
}		/*end print_line*/


LOCAL void set_plot_control(
	const char	*message,
	int		*plot_control,
	int		*num_points,
	Locstate	state)
{
	double		Msq;
	char		s[Gets_BUF_SIZE];

	
	if (message == NULL) message = "";

	Msq = mach_number_squared(state,(double *)NULL,(double *)NULL);
	if (Msq < 1.0)
	{
		screen("The state ahead of the %sshock and rarefaction polar ",
			message);
		screen("is subsonic\n");
		*plot_control = DONT_PLOT;
		return;
	}

	screen("The %sshock and rarefaction polars are ",message);
	screen("divided into four parts.\n\t");
	screen("You will now be prompted or your choices on plotting ");
	screen("portions of the ");
	screen("%sshock and rarefaction polars.\n",message);
	*plot_control = DONT_PLOT;
	screen("Plot the right (C+) branch of the shock polar [yes]: ");
	(void) Gets(s);
	if (s[0] != 'n' && s[0] != 'N')
		*plot_control |= PLOT_SHOCK_POLAR | PLOT_CPLUS_WAVE;
	screen("Plot the left (C-) branch of the shock polar [yes]: ");
	(void) Gets(s);
	if (s[0] != 'n' && s[0] != 'N')
		*plot_control |= PLOT_SHOCK_POLAR | PLOT_CMINUS_WAVE;
	screen("Plot the left (C+) branch of the rarefaction polar [yes]: ");
	(void) Gets(s);
	if (s[0] != 'n' && s[0] != 'N')
		*plot_control |= PLOT_RAREF_POLAR | PLOT_CPLUS_WAVE;
	screen("Plot the right (C-) branch of the rarefaction polar [yes]: ");
	(void) Gets(s);
	if (s[0] != 'n' && s[0] != 'N')
		*plot_control |= PLOT_RAREF_POLAR | PLOT_CMINUS_WAVE;

	*plot_control |= PLOT_SUPERSONIC;
	*plot_control |= PLOT_SUBSONIC;
	if (plot_shock_polar(*plot_control))
	{
		screen("Plot supersonic portion of shock polar [yes]: ");
		(void) Gets(s);
		if (s[0] == 'n' || s[0] == 'N')
			*plot_control &= ~PLOT_SUPERSONIC;
		screen("Plot subsonic portion of shock polar [yes]: ");
		(void) Gets(s);
		if (s[0] == 'n' || s[0] == 'N')
			*plot_control &= ~PLOT_SUBSONIC;
	}

	screen("Enter the number of points on the shock polar ");
	screen("[%d]: ",*num_points);
	(void) Gets(s);
	if (s[0] != '\0')
	{
		(void) sscanf(s,"%d",num_points);
		*num_points = min(*num_points,MAX_NUM_POINTS);
	}
}		/*end set_plot_control*/

LOCAL	void plot_polar(
	double		theta0,
	double		M0sq,
	Locstate	state0,
	double		**var,
	double		pmin,
	double		pmax,
	double		*var_l,
	double		*var_u,
	int		*pnpts,
	int		plt_ctrl)
{
	int		npts = *pnpts;
	double		dp;
	double		p, p0, rho0, q0sq, H0, ais0;
	double		t_ang;
	double		rhodiff_min;
	int		i, k = 2*npts-1;
	static Locstate st1 = NULL;

	if (st1 == NULL) g_alloc_state(&st1,sizeof(VGas));

	Set_params(st1,state0);
	set_type_of_state(st1,state_type(state0));
	H0 = specific_enthalpy(state0);
	ais0 = acoustic_impedance_squared(state0);
	q0sq = M0sq * sound_speed_squared(state0);
	rho0 = Dens(state0);
	p0 = pressure(state0);
	if (!plot_subsonic(plt_ctrl))
	{
		pmax = min(pmax,pressure_at_sonic_point(M0sq,state0));
	}
	if (!plot_supersonic(plt_ctrl))
	{
		pmin = max(pmin,pressure_at_sonic_point(M0sq,state0));
	}
	rhodiff_min = EPSILON*rho0;
	if (use_logp_increments)
		dp = pow((1.0 + pmax)/(1.0 + pmin),1.0/((double) (npts - 1)));
	else
		dp = (pmax - pmin)/((double) (npts - 1));
	for (i = 0; i < npts; ++i)
	{
		if (use_logp_increments)
			p = (1.0 + pmin)*pow(dp,(double) i) - 1.0;
		else
			p = pmin + ((double) i)*dp;
		t_ang = turning_angle(p,M0sq,state0);

		if (p > p0)
		{
			state_w_pr_on_Hugoniot(state0,p,st1,state_type(st1));
		}
		else
		{
			state_on_adiabat_with_pr(state0,p,st1,state_type(st1));
		}
		set_graph_line(var[i],var[k-i],var_l,var_u,
			st1,p,t_ang,theta0,p0,rho0,H0,ais0,q0sq,
			rhodiff_min,plt_ctrl);

	}
}		/*end plot_polar*/

LOCAL	void calculate_sonic_locus(
	Locstate	state0,
	Locstate	state1,
	double		**var,
	double		*var_l,
	double		*var_u,
	int		npts,
	int		plt_ctrl)
{
	double		qbase_sqr;
	double		s;
	double		t_ang, s_ang;
	double		vx, vy, vx0, vy0, vx1, vy1;
	double		p0, rho0, H0, ais0, rhodiff_min;
	double		nptsm1 = (double) (npts - 1);
	int		i, k = 2*npts-1;
	static Locstate base_state = NULL, sonic_state = NULL;
	static double	node_v[MAXD]; /*Statics initialized to zero*/

	if (base_state == NULL)
	{
	    g_alloc_state(&base_state,sizeof(VGas));
	    g_alloc_state(&sonic_state,sizeof(VGas));
	}
	copy_state(base_state,state0);
	p0 = pressure(state0);
	rho0 = Dens(state0);
	H0 = specific_enthalpy(state0);
	ais0 = acoustic_impedance_squared(state0);
	rhodiff_min = EPSILON*rho0;
	vx0 = vel(0,state0);
	vy0 = vel(1,state0);
	vx1 = vel(0,state1);
	vy1 = vel(1,state1);
	for (i = 0; i < npts; ++i)
	{
	    s = ((double) i)/nptsm1;
	    vx = (1.0 - s)*vx0 + s*vx1;	vy = (1.0 - s)*vy0 + s*vy1;
	    if (state_type(base_state) == GAS_STATE)
	    {
	    	Mom(base_state)[0] = Dens(base_state) * vx;
	    	Mom(base_state)[1] = Dens(base_state) * vy;
	    }
	    else
	    {
	    	Vel(base_state)[0] = vx;
	    	Vel(base_state)[1] = vy;
	    }
	    qbase_sqr = sqr(vx) + sqr(vy);
	    Check_return(
		s_polar_3(base_state,YES,
			  pressure_at_sonic_point(0.0,base_state),
			  NO,NO,node_v,sonic_state,&s_ang,&t_ang),
		calculate_sonic_locus)
	    set_graph_line(var[i],var[k-i],var_l,var_u,sonic_state,
			   pressure(sonic_state),t_ang,0.0,p0,rho0,H0,
			   ais0,qbase_sqr,rhodiff_min,plt_ctrl);
	}
}		/*end calculate_sonic_locus*/

LOCAL	void set_graph_line(
	double		*var,
	double		*ovar,
	double		*var_l,
	double		*var_u,
	Locstate	st1,
	double		p,
	double		t_ang,
	double		theta0,
	double		p0,
	double		rho0,
	double		H0,
	double		ais0,
	double		q0sq,
	double		rhodiff_min,
	int		plt_ctrl)
{
	double		rho, V, S, q, H, c, M, mf, vx, vy, gam[2];
	double		t_plus, t_minus, stab;

	var[PRESSURE] = ovar[PRESSURE] = p;
	var_u[PRESSURE] = max(p,var_u[PRESSURE]);
	var_l[PRESSURE] = min(p,var_l[PRESSURE]);

	ovar[TURN_ANG] = t_plus = theta0 + t_ang;
	if (plot_cplus_wave(plt_ctrl))
	{
		var_u[TURN_ANG] = max(t_plus,var_u[TURN_ANG]);
		var_l[TURN_ANG] = min(t_plus,var_l[TURN_ANG]);
	}

	var[TURN_ANG] = t_minus = theta0 - t_ang;
	if (plot_cminus_wave(plt_ctrl))
	{
		var_u[TURN_ANG] = max(t_minus,var_u[TURN_ANG]);
		var_l[TURN_ANG] = min(t_minus,var_l[TURN_ANG]);
	}

	rho = Dens(st1);
	var[DENSITY] = ovar[DENSITY] = rho;
	var_u[DENSITY] = max(rho,var_u[DENSITY]);
	var_l[DENSITY] = min(rho,var_l[DENSITY]);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	V = (rho > Vacuum_dens(st1)) ? 1.0/rho : HUGE_VAL;
#else /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	V = 1.0/rho;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	var[SPEC_VOL] = ovar[SPEC_VOL] = V;
	var_u[SPEC_VOL] = max(V,var_u[SPEC_VOL]);
	var_l[SPEC_VOL] = min(V,var_l[SPEC_VOL]);

	c = sound_speed(st1);
	var[SOUND_SPEED] = ovar[SOUND_SPEED] = c;
	var_u[SOUND_SPEED] = max(c,var_u[SOUND_SPEED]);
	var_l[SOUND_SPEED] = min(c,var_l[SOUND_SPEED]);

	S = entropy(st1);
	var[ENTROPY] = ovar[ENTROPY] = S;
	var_u[ENTROPY] = max(S,var_u[ENTROPY]);
	var_l[ENTROPY] = min(S,var_l[ENTROPY]);

	gam[0] = adiabatic_gamma(st1);
	var[ADIABATIC_GAMMA] = ovar[ADIABATIC_GAMMA] = gam[0];
	var_u[ADIABATIC_GAMMA] = max(gam[0],var_u[ADIABATIC_GAMMA]);
	var_l[ADIABATIC_GAMMA] = min(gam[0],var_l[ADIABATIC_GAMMA]);

	gam[1] = gruneisen_gamma(st1);
	var[GRUNEISEN_GAMMA] = ovar[GRUNEISEN_GAMMA] = gam[1];
	var_u[GRUNEISEN_GAMMA] = max(gam[1],var_u[GRUNEISEN_GAMMA]);
	var_l[GRUNEISEN_GAMMA] = min(gam[1],var_l[GRUNEISEN_GAMMA]);

	stab = gam[0] - 1.0 - gam[1];
	var[TWO_D_STABILITY] = ovar[TWO_D_STABILITY] = stab;
	var_u[TWO_D_STABILITY] = max(stab,var_u[TWO_D_STABILITY]);
	var_l[TWO_D_STABILITY] = min(stab,var_l[TWO_D_STABILITY]);

	mf = (fabs(rho - rho0) < rhodiff_min) ? 1.0 :
		rho*rho0*(p-p0)/(ais0*(rho - rho0));
	var[MASS_FLUX] = ovar[MASS_FLUX] = mf;
	var_u[MASS_FLUX] = max(mf,var_u[MASS_FLUX]);
	var_l[MASS_FLUX] = min(mf,var_l[MASS_FLUX]);

	H = specific_enthalpy(st1);
	var[ENTHALPY] = ovar[ENTHALPY] = H;
	var_u[ENTHALPY] = max(H,var_u[ENTHALPY]);
	var_l[ENTHALPY] = min(H,var_l[ENTHALPY]);

	q = sqrt(q0sq + 2.0*(H0 - H));
	var[FLOW_SPEED] = ovar[FLOW_SPEED] = q;
	var_u[FLOW_SPEED] = max(q,var_u[FLOW_SPEED]);
	var_l[FLOW_SPEED] = min(q,var_l[FLOW_SPEED]);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	M = (rho > Vacuum_dens(st1)) ? q/c : HUGE_VAL;
#else /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	M = q/c;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	var[MACH_NUMBER] = ovar[MACH_NUMBER] = M;
	var_u[MACH_NUMBER] = max(M,var_u[MACH_NUMBER]);
	var_l[MACH_NUMBER] = min(M,var_l[MACH_NUMBER]);

	var[VELX] = vx = q*cos(var[TURN_ANG]);
	var[VELY] = vy = q*sin(var[TURN_ANG]);

	if (plot_cminus_wave(plt_ctrl))
	{
		var_u[VELX] = max(vx,var_u[VELX]);
		var_l[VELX] = min(vx,var_l[VELX]);

		var_u[VELY] = max(vy,var_u[VELY]);
		var_l[VELY] = min(vy,var_l[VELY]);
	}

	ovar[VELX] = vx = q*cos(ovar[TURN_ANG]);
	ovar[VELY] = vy = q*sin(ovar[TURN_ANG]);


	if (plot_cplus_wave(plt_ctrl))
	{
		var_u[VELX] = max(vx,var_u[VELX]);
		var_l[VELX] = min(vx,var_l[VELX]);

		var_u[VELY] = max(vy,var_u[VELY]);
		var_l[VELY] = min(vy,var_l[VELY]);
	}
}		/*end set_graph_line*/

LOCAL	double turning_angle(
	double		p,
	double		M0sq,
	Locstate	state0)
{
	double		theta;

	if (steady_state_wave_curve(p,M0sq,&theta,state0) == FUNCTION_FAILED)
	{
	    screen("ERROR in turning_angle(), "
	           "steady_state_wave_curve() failed\n");
	    clean_up(ERROR);
	}
	return theta;
}		/*end turning_angle*/

LOCAL	int normal_shock_refraction(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	Locstate	*state)
{
	Gas_param	         *params;
	double		         pml, pmr, uml, umr, mr, ml; /* midstate */
	double		         s;
	double		         m;
	double		         cfac;
	double		         Am, Ap;
	RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;		/* SHOCK or RAREFACTION */
	static double	         node_v[MAXD]; /*Statics initialized to zero*/
	static double	         nor[] = {1.0, 0.0, 0.0};
	static Locstate          sl = NULL, sml = NULL, smr = NULL, sr = NULL;

	if (sl == NULL)
	{
	    g_alloc_state(&sl,sizeof(VGas));
	    g_alloc_state(&sml,sizeof(VGas));
	    g_alloc_state(&smr,sizeof(VGas));
	    g_alloc_state(&sr,sizeof(VGas));
	}
	params = init_eos_params(init,ip,
				 "\n\tahead of the incident shock wave",NO);
	prompt_for_state(state[0],UNKNOWN_STATE,params," of the incident gas");

	prompt_for_behind_shock_state(state[0],state[2],NO,nor,
				      state_type(state[0]),YES,init);
	params = init_eos_params(init,ip,
				 " ahead of the transmitted shock wave",NO);
	init_state_across_contact(state[1],state[0],params,node_v,
				" of the gas ahead of the transmitted shock");
	set_state(sl,VGAS_STATE,state[2]);
	set_state(sr,VGAS_STATE,state[1]);
	(void) find_mid_state(sl,sr,0.0,&pml,&pmr,&uml,&umr,&ml,&mr,
			      &l_wave,&r_wave);
	if (l_wave == SHOCK)
	    state_w_pr_on_Hugoniot(sl,pml,sml,TGAS_STATE);
	else
	    state_on_adiabat_with_pr(sl,pml,sml,TGAS_STATE);
	if (r_wave == SHOCK)
	    state_w_pr_on_Hugoniot(sr,pmr,smr,TGAS_STATE);
	else
	    state_on_adiabat_with_pr(sr,pmr,smr,TGAS_STATE);
	Vel(sml)[0] = uml;
	Vel(smr)[0] = umr;

	(void) printf("\n\nRIEMANN SOLUTION FOR SHOCK REFRACTION\n");
	verbose_print_state("State ahead of incident shock",state[0]);
	verbose_print_state("\nState behind incident shock",sl);
	verbose_print_state("\nState behind reflected wave",sml);
	verbose_print_state("\nState behind transmitted shock",smr);
	verbose_print_state("\nState ahead of transmitted shock",sr);
	m = mass_flux(pressure(sl),state[0]);
	(void) printf("\nIncident shock mass flux = %g\n",m);
	(void) printf("Incident wave pressure ratio = %g\n",
		pressure(sl)/pressure(state[0]));
	(void) printf("Incident wave density ratio = %g\n",
		      Dens(sl)/Dens(state[0]));
	s = vel(0,state[0]) + fabs(m)/Dens(state[0]);
	(void) printf("Incident shock velocity = %g\n",s);
	(void) printf("Incident shock Mach number = %g\n",
		fabs(m)/acoustic_impedance(state[0]));
	(void) printf("Incident shock behind Mach number = %g\n",
		fabs(m)/acoustic_impedance(sl));
	Am = (Dens(state[1]) - Dens(state[0])) /
	     (Dens(state[0]) + Dens(state[1]));
	Ap = (Dens(smr) - Dens(sml))/(Dens(smr) + Dens(sml));
	cfac = 1.0 - 0.5*(uml+umr)/s;
	(void) printf("Compression factor = %g\n",cfac);
	(void) printf("Preshock Atwood number = %g\n",Am);
	(void) printf("Postshock Atwood number = %g\n",Ap);
		
	(void) printf("Reflected wave mass flux = %g\n",ml);
	(void) printf("Reflected wave pressure ratio = %g\n",pml/pressure(sl));
	(void) printf("Reflected wave density ratio = %g\n",Dens(sml)/Dens(sl));
	if (l_wave == SHOCK)
	{
	    (void) printf("Reflected shock mass flux = %g\n",ml);
	    s = vel(0,sl) - fabs(ml)/Dens(sl);
	    (void) printf("Reflected shock velocity = %g\n",s);
	    (void) printf("Reflected shock ahead Mach number = %g\n",
	    	fabs(ml)/acoustic_impedance(sl));
	    (void) printf("Reflected shock behind Mach number = %g\n",
	    	fabs(ml)/acoustic_impedance(sml));
	}
	else
	{
	    (void) printf("Reflected rarefaction ");
	    (void) printf("leading edge impedance = %g\n",
			  acoustic_impedance(sl));
	    s = vel(0,sl) - sound_speed(sl);
	    (void) printf("Reflected rarefaction ");
	    (void) printf("leading edge velocity = %g\n",s);

	    (void) printf("Reflected rarefaction ");
	    (void) printf("trailing edge impedance = %g\n",
	    	acoustic_impedance(sml));
	    s = vel(0,sml) - sound_speed(sml);
	    (void) printf("Reflected rarefaction ");
	    (void) printf("trailing edge velocity = %g\n",s);
	}
	(void) printf("Transmitted shock mass flux = %g\n",mr);
	(void) printf("Transmitted wave pressure ratio = %g\n",
		      pmr/pressure(sr));
	(void) printf("Transmitted wave density ratio = %g\n",
		      Dens(smr)/Dens(sr));
	s = vel(0,sr) + fabs(mr)/Dens(sr);
	(void) printf("Transmitted shock velocity = %g\n",s);
	(void) printf("Transmitted shock ahead Mach number = %g\n",
		fabs(mr)/acoustic_impedance(sr));
	(void) printf("Transmitted shock behind Mach number = %g\n",
		fabs(mr)/acoustic_impedance(smr));
	clean_up(0);
	return 0;
}		/*end normal_shock_refraction*/


/*
*				start_up():
*
*	Calls Routines to handle system error messages, and prints
*	Current messages at top of output or on the screen.
*
*/

#define SHIFT (--argc,++argv)

LOCAL	void start_up(
	int	     argc,
	char	     **argv,
	INIT_DATA    *init)
{
	IMPORT	boolean	suppress_prompts;

	init_clean_up(spolars_clean_up,NULL);

	SHIFT;
	while (argc && argv[0][0]=='-')
	{
	    if (strncmp(argv[0],"-i",2) == 0)
	    {
	    	static char infile[1024];
	    	SHIFT;
		temporary_input_file = infile;
	    	stripcomm(infile,argv[0]);
	    	if (freopen(infile,"r",stdin) == NULL)
	    	{
	    	    screen("ERROR in start_up(), "
	    	           "can't reopen %s to %s\n","stdin",infile);
	    	    clean_up(ERROR);
	    	}
	    	suppress_prompts = YES;
	    }
	    else if (strncmp(argv[0],"-o",2) == 0)
	    {
	    	SHIFT;
	    	if (freopen(argv[0],"w",stdout) == NULL)
	    	{
	    	    screen("ERROR in start_up(), "
	    	           "can't reopen %s to %s\n","stdout",argv[0]);
	    	    clean_up(ERROR);
	    	}
	    }
	    else if (strncmp(argv[0],"-e",2) == 0)
	    {
	    	SHIFT;
	    	if (freopen(argv[0],"w",stderr) == NULL)
	    	{
	    	    screen("ERROR in start_up(), ");
	    	    screen("can't reopen %s to %s\n","stderr",argv[0]);
	    	    clean_up(ERROR);
	    	}
	    }
	    SHIFT;
	}

	set_error_immediate(stdout);
	record_print_version(stdout);
	print_title(stdout,title(init));
	init_prompting_and_debugging(init);
}		/*end start_up*/

LOCAL	void spolars_clean_up(void)
{
	if ((temporary_input_file != NULL) && (is_io_node(pp_mynode())))
	    (void) unlink(temporary_input_file);
}		/*end spolars_clean_up*/


LOCAL	void init_state_across_contact(
	Locstate	state1,
	Locstate	state0,
	Gas_param	*params,
	double		*node_v,
	const char	*mesg)
{
	int		i, dim;
	char		s[Gets_BUF_SIZE];
	double		slip, alpha;
	double		v0[MAXD];
	double		M1, q1, c1;
	static	char	lmesg = '\0';

	if (mesg == NULL)
	    mesg = &lmesg;
	Init_params(state1,params);
	set_type_of_state(state1,TGAS_STATE);
	dim = params->dim;

	screen("Enter the density%s: ",mesg);
	(void) Scanf("%f\n",&Dens(state1));
	Press(state1) = pressure(state0);
	for (i = 0; i < dim; ++i)
	    Vel(state1)[i] = vel(i,state0);
        reset_gamma(state1);
		

	if (dim == 1)
	    return;
	screen("The velocity data may be entered in "
	       "one of the following ways\n"
	       "\tDirect entry of velocity (V)\n"
	       "\tSteady Mach number (M)\n"
	       "\tSteady flow speed (Q)\n"
	       "\tSlip across contact (S)\n"
	       "\tNo slip (N, default)\n"
	       "Enter choice: ");
	(void) Gets(s);
	switch (s[0])
	{
	case 'V':
	case 'v':
	    screen("Enter the %d components of the state velocity: ",dim);
	    for (i = 0; i < dim; ++i)
	    	(void) Scanf("%f",&Vel(state1)[i]);
	    screen("Enter an optional velocity to translate the "
	           "above velocities into\n\tthe steady state frame: ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        double	nvx, nvy;
	    	(void) sscanf(s,"%lf %lf",&nvx,&nvy);

	    	Vel(state1)[0] -= nvx;
	    	Vel(state1)[1] -= nvy;
	    }
	    break;
	case 'M':
	case 'm':
	    screen("Enter the steady state Mach number of the state: ");
	    (void) Scanf("%f\n",&M1);
	    c1 = sound_speed(state1);
	    for (i = 0; i < dim; ++i)
	    	v0[i] = vel(i,state0) - node_v[i];
	    alpha = c1*M1/mag_vector(v0,dim);
	    for (i = 0; i < dim; ++i)
	    	Vel(state1)[i] = alpha*v0[i] + node_v[i];
	    break;
	case 'Q':
	case 'q':
	    screen("Enter the steady state flow speed of the state: ");
	    (void) Scanf("%f\n",&q1);
	    for (i = 0; i < dim; ++i)
	    	v0[i] = vel(i,state0) - node_v[i];
	    alpha = q1/mag_vector(v0,dim);
	    for (i = 0; i < dim; ++i)
	    	Vel(state1)[i] = alpha*v0[i] + node_v[i];
	    break;
	case 'S':
	case 's':
	    slip = 0.0;
	    screen("Enter the the slip across the contact (dflt = 0): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&slip);
	    }
	    for (i = 0; i < dim; ++i)
	    	v0[i] = vel(i,state0) - node_v[i];
	    alpha = 1.0 + slip/mag_vector(v0,dim);
	    for (i = 0; i < dim; ++i)
	    	Vel(state1)[i] = alpha*v0[i] + node_v[i];
	    break;
	case 'N':
	case 'n':
	default:
	    break;
	}
	return;
}		/*end init_state_across_contact*/

#endif /* defined(TWOD) || defined(THREED) */
