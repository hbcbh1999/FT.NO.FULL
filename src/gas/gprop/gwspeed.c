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
*
*				gwspeed.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	g_w_speed() and g_npt_w_speed() calculate the wave speeds and 
*	states for a tracked wave.  The first-order computation is done
*	in g_w_speed().  The function g_npt_w_speed() calculates the
*	second-order corrections.
*
*	Source-term corrections to front states due to gravitation or 
*	cylindrical symmetry are included in g_npt_w_speed() through
*	the method of characteristics or throught the function include_source().
*	Also the characteristic equations in
*	g_npt_w_speed() are modified due to cylindrical source terms,
*	which become important near the axis of symmetry.
*/

#include <gdecs/gdecs.h>

	/* LOCAL Function Prototypes */
LOCAL	boolean    invalid_shock(double,double,double,double);
LOCAL	boolean	strong_interaction(Locstate,Locstate,int,
	                           double,double,double,double);
LOCAL	boolean	strong_wave(Locstate,double,double);
LOCAL	boolean	fC(double,double*,POINTER);
LOCAL   boolean    pressure_switch_on(Locstate,Locstate,Locstate);
LOCAL	double	wall_limiter(WSSten*,SIDE,NptWSpeedOpts*);
LOCAL	void	apply_bc(Locstate,Wave*,HYPER_SURF*);
LOCAL	void	IncomingStatesAtContact(WSSten*,double,double,
	                                Locstate,Locstate,double*,double*,
	                                NptWSpeedOpts*);
LOCAL   void    bdry_npt_w_speed(WSSten*,Locstate,Locstate,NptWSpeedOpts*);
LOCAL	void	bdry_wspeed(int,double*,Locstate,Locstate,Locstate,Locstate,
	                    int,Front*);
LOCAL	void	contact_cheap_moc(double*,Locstate,Locstate,Locstate,Locstate,
				  Locstate,Locstate,double,double,double,double,
				  double,double*,double*,Front*);
LOCAL	void	contact_filter_outgoing_wave(Locstate,Locstate,Locstate,double,
	                                     WSSten*,SIDE,NptWSpeedOpts*);
LOCAL	void	contact_moc_plus_rh(double*,double,Locstate,Locstate,Locstate,
			            Locstate,Locstate,Locstate,double,double,
				    double,double*,double*,Front*);
LOCAL	void	contact_riemann(Locstate,Locstate,Locstate,Locstate,
	                     Locstate,double,Locstate,double,Locstate,Locstate,
	                     WSSten*,double*,NptWSpeedOpts*,double*,double,double);
LOCAL	void	g_Left_w_speed(const char*,Locstate,Locstate,double*,int);
LOCAL	void	midstate(Locstate,Locstate,double,double,double,int,
                         RIEMANN_SOLVER_WAVE_TYPE,int);
LOCAL	void	neumann_cheap_moc(double*,Locstate,double,double,Locstate,
	                          SIDE,Locstate,double,double*,Front*);
LOCAL	void	onedw_speed(Locstate,Locstate,int,Locstate,Locstate,Locstate,
	                    Locstate,double*,int,double,Front*);
LOCAL	void	shock_ahead_state_cheap_moc(double*,Locstate,Locstate,Locstate,
	                                    Locstate,Locstate,double,double,double,
					    double,double*,double*,int,double,
					    Front*);
LOCAL	void	wspeed_neumann_riem_inv_moc(double*,Locstate,double,double,
	                                    Locstate,SIDE,Locstate,double,double*,
	                                    Front*);
LOCAL	void	wspeed_shock_ahead_state_riem_inv_moc(double*,Locstate,Locstate,
	                                              Locstate,Locstate,
						      Locstate,double,double,
						      double,double,double*,
						      double*,int,double,Front*);

#if defined(COMBUSTION_CODE)
LOCAL	void	CJ_state_behind_curved_wave(Locstate,Locstate,double,double*,
	                                    int,int);
LOCAL   void flame_moc_plus_rh(double*,Locstate,Locstate,Locstate,
                                    Locstate,Locstate,Locstate,Locstate,
                                    double,double,double,double,double,double*,
                                    double*,double,double,double*,double,double,
                                    Front*);
LOCAL  boolean flame_moc_plus_rh_signed(Locstate,Locstate,Locstate,
				    Locstate,Locstate,boolean,
				    double,double,double,double*,
				    double*,double*,int);
LOCAL  double flame_moc_plus_rh_guessp(Locstate,Locstate,Locstate,
				    Locstate,Locstate,double,
				    double,double,double,double*,
				    double*,double*,int);
LOCAL  double flame_moc_plus_rh_deltap(Locstate,Locstate,Locstate,
				    Locstate,Locstate,double,
				    double,double,double,double*,
				    double*,double*,int);
LOCAL  double flame_moc_plus_rh_findstates(double,double,
				    Locstate,Locstate,Locstate,
				    Locstate,Locstate,
				    double,double,double,double*,
				    double*,double*,int);
#endif /* defined(COMBUSTION_CODE) */

#if !defined(__INTEL_COMPILER)
#pragma	noinline	strong_interaction
#pragma	noinline	fC
#pragma noinline        bdry_npt_w_speed
#pragma	noinline	bdry_wspeed
#pragma	noinline	contact_cheap_moc
#pragma	noinline	contact_moc_plus_rh
#pragma	noinline	midstate
#pragma	noinline	neumann_cheap_moc
#pragma	noinline	contact_riemann
#pragma	noinline	shock_ahead_state_cheap_moc
#if defined(COMBUSTION_CODE)
#pragma	noinline	CJ_state_behind_curved_wave
#endif /* defined(COMBUSTION_CODE) */
#endif /*!defined(__INTEL_COMPILER)*/

LOCAL	NptWSpeedOpts	G3ptOpts = {
	-1.0,		/*Mach_tol*/
	-1.0,		/*A_tol*/
	-1.0,		/*Wall_limiter*/
	MOC_TYPE_UNSET,	/*vector_moc*/
	MOC_TYPE_UNSET,	/*scalar_moc*/
	NO,		/*_scalar_filter_outgoing_waves_at_contact*/
	NO,		/*use_cheap_ahead_state_moc*/
	NO,		/*use_neumann_cheap_moc*/
	NULL,
	NULL,
	IDENTITY_REMAP,
	MODIFIED_EULER
};

EXPORT	void	set_npt_wspeed_options(
	NptWSpeedOpts	*opts)
{
	if (opts->vector_ahead_state_moc == NULL)
	{
	    opts->vector_ahead_state_moc =
		(opts->use_cheap_ahead_state_moc == YES) ?
		shock_ahead_state_cheap_moc :
		wspeed_shock_ahead_state_riem_inv_moc;
	}
	if (opts->neumann_moc == NULL)
	{
	    opts->neumann_moc = (opts->use_neumann_cheap_moc == YES) ?
		                neumann_cheap_moc :
		                wspeed_neumann_riem_inv_moc;
	}
	G3ptOpts = *opts;
}		/*end set_npt_wspeed_options*/

#if !defined(UNRESTRICTED_THERMODYNAMICS)
#define pressure_is_near_vacuum(p,st) ((p) < 2.0*Min_pressure(st)) /*TOLERANCE*/
#else /* !defined(UNRESTRICTED_THERMODYNAMICS) */
#define pressure_is_near_vacuum(p,st)	(NO)
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */



#if defined(DEBUG_W_SPEED)

LOCAL	boolean debug_w_speed = NO;

#define	Left_w_speed(f,ansl,ansr,w,dim) g_Left_w_speed(f,ansl,ansr,w,dim)

EXPORT void g_set_debug_w_speed(boolean y_or_n)
{
	debug_w_speed = y_or_n;
}	/*end g_set_debug_w_speed*/

LOCAL	void	g_Left_w_speed(
	const char	*func,
	Locstate	ansl,
	Locstate	ansr,
	double		*w,
	int		dim)
{
	double		len;
	int		i;

	if (debug_w_speed == YES) 
	{
	    (void) printf("Final left, right states:\n");
	    verbose_print_state("ansl",ansl);
	    verbose_print_state("ansr",ansr);
	}
	if (debugging("wave_speeds")) 
	{
	    if (w != NULL)
	    {
	        (void) printf("\nFinal wave velocity in %s() = ",func);
	        print_general_vector("",w,dim,"");
	        for (len = 0.0, i = 0; i < dim; ++i)
	            len += sqr(w[i]);
	        len = sqrt(len);
	        (void) printf("speed = %g\n",len);
	    }
	}
	if ((debug_w_speed == YES) && (!debugging("w_speed")))
	    (void) printf("Left %s()\n",func);
	debug_print("w_speed","Left %s()\n",func);
}		/*end g_Left_w_speed*/


#else /* defined(DEBUG_W_SPEED) */

	/*   Make Left_w_speed a macro call so that  */
	/* there is no function call when not debugging */

#define	Left_w_speed(func,ansl,ansr,w,dim)

#endif /* defined(DEBUG_W_SPEED) */

EXPORT	void	print_g_npt_w_speed_opts(
	NptWSpeedOpts	*opts)
{
	printf("Current values for options for g_npt_w_speed\n");
	printf("\tA wave is defined to be strong if "
	       "|1 - (1/(rho*c)*|dp/du|| > Mach_tol or\n"
	       "\t\t|rhol - rhor|/(rhol+rhor) > A_tol\n");
	printf("\tMach_tol = %g\n",opts->Mach_tol);
	printf("\tA_tol = %g\n",opts->A_tol);
	printf("\tNeumann boundary states are computed by an average "
	       "of a reflection\n"
	       "\t\tsymmetry contact propagation and a method of characterics\n"
	       "\t\tcalculation.  The weight of the symmetry contact result\n"
	       "\t\tis proportional to the flow gradient.  The wall limiter\n"
	       "\t\tvalue gives this proportionality constant.\n");
	printf("\tWall_limiter = %g\n",opts->Wall_limiter);
	printf("\tvector_moc = ");
	switch (opts->vector_moc)
	{
	case MOC_PLUS_RIEMANN:
	    printf("MOC_PLUS_RIEMANN\n");
	    break;
	case MOC_PLUS_RH:
	    printf("MOC_PLUS_RH\n");
	    break;
	case RIEMANN:
	    printf("RIEMANN\n");
	    break;
	default:
	    printf("ERROR in print_g_npt_w_speed_opts(), "
	           "Invalid value %d for vector_moc\n",opts->vector_moc);
	    clean_up(ERROR);
	}
	printf("\tscalar_moc = ");
	switch (opts->scalar_moc)
	{
	case MOC_PLUS_RIEMANN:
	    printf("MOC_PLUS_RIEMANN\n");
	    break;
	case MOC_PLUS_RH:
	    printf("MOC_PLUS_RH\n");
	    break;
	case RIEMANN:
	    printf("RIEMANN "
	           "%s filtering of outgoing waves\n",
	           (opts->_scalar_filter_outgoing_waves_at_contact == YES) ?
	           "with" : "without");
	    break;
	default:
	    printf("ERROR in print_g_npt_w_speed_opts(), "
	           "Invalid value %d for scalar_moc\n",opts->scalar_moc);
	    clean_up(ERROR);
	}
	/*
	if (opts->vector_ahead_state_moc == NULL)
	{
	    opts->vector_ahead_state_moc =
		(opts->use_cheap_ahead_state_moc == YES) ?
		shock_ahead_state_cheap_moc :
		wspeed_shock_ahead_state_riem_inv_moc;
	}
	*/
	printf("\tvector_ahead_state_moc = %p ",opts->vector_ahead_state_moc);
	if (opts->vector_ahead_state_moc ==
	                                wspeed_shock_ahead_state_riem_inv_moc)
	    printf("shock_ahead_state_riem_inv_moc\n");
	else if (opts->vector_ahead_state_moc ==
	                                shock_ahead_state_cheap_moc)
	    printf("shock_ahead_state_cheap_moc\n");
	else
	    printf("user defined\n");

	/*
	if (opts->neumann_moc == NULL)
	{
	    opts->neumann_moc = (opts->use_neumann_cheap_moc == YES) ?
		                neumann_cheap_moc :
		                wspeed_neumann_riem_inv_moc;
	}
	*/
	printf("\tneumann_moc = %p ",opts->neumann_moc);
	if (opts->neumann_moc == wspeed_neumann_riem_inv_moc)
	    printf("neumann_riem_inv_moc\n");
	else if (opts->neumann_moc == neumann_cheap_moc)
	    printf("neumann_cheap_moc\n");
	else
	    printf("user defined\n");
	if (opts->remap != IDENTITY_REMAP)
	{
	    printf("\tgeom_source_method = ");
	    switch (opts->geom_source_method)
	    {
	    case MODIFIED_EULER:
	        printf("MODIFIED_EULER");
	        break;
	    case ANALYTIC_SOLUTION:
	        printf("ANALYTIC_SOLUTION");
	        break;
	    case BACKWARD_EULER:
	        printf("BACKWARD_EULER");
	        break;
	    }
	    printf("\n");
	}
	printf("End current values for options for g_npt_w_speed\n");
}		/*end print_g_npt_w_speed_opts*/



/*
*			g_w_speed():
*
*	Computes wave speeds and states for a tracked wave to first order.
*
*	Input:	sl,sr   - left and right states close to oriented front
*		nor     - unit normal to front (left->right)
*		w_type  - type of tracked wave
*		dt	- time step, used for the inclusion of source terms.
*		pjump   - For surface tension,  the pressure on the negative
*                         side of a contact is equal to the pressure on the
*                         positive side plus pjump.
*
*	Output:	ansl,ansr - updated left and right states.
*		W         - calculated wave velocity
*
*/

/*ARGSUSED*/
EXPORT	void g_w_speed(
	double		*pt,
	Locstate	sl,
	Locstate	sr,
	Locstate	ansl,
	Locstate	ansr,
	double		*W,
	double		pjump,
	double		*nor,
	int		w_type,
	Front		*fr)
{
	int		dim = fr->rect_grid->dim;
	static size_t	sizest = 0;
	static Locstate Tsl = NULL, Tsr = NULL, left = NULL, right = NULL;

#if defined(COMBUSTION_CODE)
	static Locstate CJ = NULL;
#endif /* defined(COMBUSTION_CODE) */

	double s;		/* normal wave speed */
	double vtanl[MAXD], vtanr[MAXD];	/* tangential velocity */
	double vnorl, vnorr;     /* normal velocity              */
	int i;

	if (is_obstacle_state(sl) && is_obstacle_state(sr)) 
	{
	    obstacle_state(fr->interf,ansl,fr->sizest);
	    obstacle_state(fr->interf,ansr,fr->sizest);
	    for (i = 0; i < dim; ++i)
	        W[i] = 0.0;
	    return;
	}
	if (right == NULL)
	{
	    Gas_param *params=(!is_obstacle_state(sl)) ? Params(sl):Params(sr);
	    sizest = params->sizest;
	    (*params->_alloc_state)(&Tsl,sizest);
	    (*params->_alloc_state)(&Tsr,sizest);
	    (*params->_alloc_state)(&left,max(sizeof(VGas),sizest));
	    (*params->_alloc_state)(&right,max(sizeof(VGas),sizest));
#if defined(COMBUSTION_CODE)
	    (*params->_alloc_state)(&CJ,sizest);
#endif /* defined(COMBUSTION_CODE) */
	}

#if defined(DEBUG_W_SPEED)
	debug_print("w_speed","Entered g_w_speed()\n");
	if (debugging("w_speed"))
	{
	    debug_w_speed = YES;
	}
	else if (debug_w_speed == YES)
	    (void) printf("Entered g_w_speed()\n");
	if (debug_w_speed == YES) 
	{
	    (void) printf("INPUT DATA TO W_SPEED()\n");
	    fprint_wave_type(stdout,"wave type = ",w_type,"\n",fr->interf);
	    print_general_vector("nor = ",nor,dim,"\n");
	    verbose_print_state("sl",sl);
	    verbose_print_state("sr",sr);
	    (void) printf("END INPUT DATA TO W_SPEED()\n");
	}
#endif /* defined(DEBUG_W_SPEED) */

	if (w_type < FIRST_PHYSICS_WAVE_TYPE)
	{

	    for (i = 0; i < dim; ++i)
	        W[i] = 0.0;
	    bdry_wspeed(w_type,nor,sl,sr,ansl,ansr,GAS_STATE,fr);
	    Left_w_speed("g_w_speed",ansl,ansr,W,dim);
	    return;
	}


	set_state(Tsl,TGAS_STATE,sl);
	set_state(Tsr,TGAS_STATE,sr);
	vnorl = vnorr = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    vnorl += nor[i] * Vel(Tsl)[i];
	    vnorr += nor[i] * Vel(Tsr)[i];
	}
	for (i = 0; i < dim; ++i)
	{
	    vtanl[i] = Vel(Tsl)[i] - nor[i] * vnorl;
	    vtanr[i] = Vel(Tsr)[i] - nor[i] * vnorr;
	}
	Vel(Tsl)[0] = vnorl;
	Vel(Tsr)[0] = vnorr;
	for (i = 1; i < dim; ++i)
	    Vel(Tsl)[i] =  Vel(Tsr)[i] = 0.0;

	onedw_speed(Tsl,Tsr,w_type,left,right,ansl,ansr,&s,TGAS_STATE,pjump,fr);
	
	for (i = dim - 1; i >= 0; --i)
	{
	    W[i] = nor[i] * s;
	    Vel(ansl)[i] =  vtanl[i] + nor[i] * Vel(ansl)[0];
	    Vel(ansr)[i] =  vtanr[i] + nor[i] * Vel(ansr)[0];
	}

	set_state(ansl,GAS_STATE,ansl);
	if (ansr != ansl)
	    set_state(ansr,GAS_STATE,ansr);

#if defined(DEBUG_W_SPEED)
	Left_w_speed("g_w_speed",ansl,ansr,W,dim);
#endif /* defined(DEBUG_W_SPEED) */
}	/*end g_w_speed*/



LOCAL	void bdry_wspeed(
	int		w_type,
	double		*nor,
	Locstate	sl,
	Locstate	sr,
	Locstate	ansl,
	Locstate	ansr,
	int		state_type,
	Front		*fr)
{
	size_t		sizest = fr->sizest;
	double		mnor;
	int		i, dim;

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	    fprint_wave_type(stdout,"wave_type = ",w_type,"\n",fr->interf);
#endif /* defined(DEBUG_W_SPEED) */
	set_type_of_state(ansl,GAS_STATE);
	set_type_of_state(ansr,GAS_STATE);
	switch (w_type)
	{
	case PASSIVE_BOUNDARY:
	case SUBDOMAIN_BOUNDARY:
	    /* PASSIVE and SUBDOMAIN states not set or used */
	    obstacle_state(fr->interf,ansl,sizest);
	    obstacle_state(fr->interf,ansr,sizest);
	    break;

	case DIRICHLET_BOUNDARY:
	case TIME_DIRICHLET_BOUNDARY:
	    copy_state(ansl,sl);
	    copy_state(ansr,sr);
	    break;

	case NEUMANN_BOUNDARY:
	    if (is_obstacle_state(sl))
	        obstacle_state(fr->interf,ansl,sizest);
	    else
	    {
	        Dens(ansl) = Dens(sl);
	        Energy(ansl) = energy(sl);
	        mnor = 0.0;
	        dim = Params(sl)->dim;
	        for (i = 0; i < dim; ++i)
	            mnor +=  nor[i] * mom(i,sl);
	        for (i = 0; i < dim; ++i)
	            Mom(ansl)[i] = mom(i,sl) - nor[i] * mnor;
	        /* QUESTION?  Should the Energy be adjusted to
	        *  offset the change in kinetic energy? */
	        /* ANSWER: Yes, if any thermodynamic call is
	        *  made, even through a subroutine, before
	        *  tangential momentum is readded to ansl.
	        *  This change has to be done everywhere, in a
	        *  consistent fashion, when it is done */

#if defined(COMBUSTION_CODE)
	        if (Composition_type(sl) == ZND)
	            Prod(ansl) = Prod(sl);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ansl);

	        Set_params(ansl,sl);
	        if (state_type != state_type(ansl)) 
	            set_state(ansl,state_type,ansl);
	    }

	    if (is_obstacle_state(sr))
	        obstacle_state(fr->interf,ansr,sizest);
	    else
	    {
	        Dens(ansr) = Dens(sr);
	        Energy(ansr) = energy(sr);
	        mnor = 0.0;
	        dim = Params(sr)->dim;
	        for (i = 0; i < dim; ++i)
	            mnor +=  nor[i] * mom(i,sr);
	        for (i = 0; i < dim; ++i)
	            Mom(ansr)[i] = mom(i,sr) - nor[i] * mnor;
	        /* QUESTION?  Should the Energy be adjusted to
	        *  offset the change in kinetic energy? */

#if defined(COMBUSTION_CODE)
	        if (Composition_type(sr) == ZND)
	            Prod(ansr) = Prod(sr); 
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ansr);

	        Set_params(ansr,sr);
	        if (state_type != state_type(ansr)) 
	            set_state(ansr,state_type,ansr);
	    }
	    break;
	
	default:
	    printf("ERROR in  bdry_wspeed(), ");
	    fprint_wave_type(stdout,"unknown wave_type = ",w_type,"\n",
	                     fr->interf);
	    clean_up(ERROR);
	}
}	/*end bdry_wspeed*/


/*
*			onedw_speed():
*
*	One dimensional wave speed operator.
*
*	Input	sl, sr		- States at time t near curve
*		w_type		- type of tracked wave
*		pwall		- Wall information (Unused if nonwall)
*
*	Output	left/right	- States corrected for wall interactions
*		Wx		- pointer to wave velocity
*		ansl, ansr	- States at time t + dt near curve
*		state_type	- type of state representation
*		pjump		- Pressure jump across contact
*
*/

LOCAL	void onedw_speed(
	Locstate	sl,
	Locstate	sr,
	int		w_type,
	Locstate	left,
	Locstate	right,
	Locstate	ansl,
	Locstate	ansr,
	double		*Wx,
	int		state_type,
	double		pjump,
	Front		*fr)
{
	RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
#if defined(ONED)
	int		dim;
#endif /* defined(ONED) */
	double		pml, pmr, uml, umr ,mr,ml; /* midstate quantities */
	double		cl, cr;
#if defined(COMBUSTION_CODE)
	static		Locstate CJ = NULL;
#endif /* defined(COMBUSTION_CODE) */
	static	size_t sizest = 0;
	
	if (is_obstacle_state(sl) && is_obstacle_state(sr)) 
	{
	    obstacle_state(fr->interf,ansl,fr->sizest);
	    obstacle_state(fr->interf,left,fr->sizest);
	    obstacle_state(fr->interf,ansr,fr->sizest);
	    obstacle_state(fr->interf,right,fr->sizest);
	    *Wx = 0.0;
	    return;
	}
	if (sizest == 0)
	{
	    Gas_param *params = (!is_obstacle_state(sl))?Params(sl):Params(sr);

	    sizest = params->sizest;
#if defined(COMBUSTION_CODE)
	    (*params->_alloc_state)(&CJ,sizest);
#endif /* defined(COMBUSTION_CODE) */

	}

#if defined(DEBUG_W_SPEED)
	debug_print("w_speed","Entered onedw_speed()\n");
	if (debugging("w_speed"))
	{
	    debug_w_speed = YES;
	}
	else if (debug_w_speed == YES)
	    (void) printf("Entered onedw_speed()\n");

	if (debug_w_speed == YES)
	{
	    fprint_wave_type(stdout,"wave_type = ",w_type,"\n",fr->interf);
	    verbose_print_state("sl",sl);
	    verbose_print_state("sr",sr);
	}
#endif /* defined(DEBUG_W_SPEED) */

	
#if defined(ONED)
	dim = (is_obstacle_state(sl)) ? Params(sr)->dim : Params(sl)->dim;
	if ((dim == 1) && ((w_type < FIRST_PHYSICS_WAVE_TYPE) || 
	                    (w_type == TIME_DIRICHLET_BOUNDARY)) )
	{
	    static double nor[3] = {1.0, 0.0, 0.0};
	    *Wx = 0.0;

	    bdry_wspeed(w_type,nor,sl,sr,ansl,ansr,state_type,fr);
	    copy_state(left,sl);
	    copy_state(right,sr);
	    Left_w_speed("onedw_speed",ansl,ansr,Wx,dim);
	    return;
	}
#endif /* defined(ONED) */


	set_type_of_state(ansl,GAS_STATE);
	set_type_of_state(ansr,GAS_STATE);
	set_state_for_find_mid_state(left,sl);
	set_state_for_find_mid_state(right,sr);


	if (find_mid_state(left,right,pjump,&pml,&pmr,&uml,&umr,&ml,&mr,
	                       &l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	{
	    (void) printf("WARNING in onedw_speed(), "
	                  "find_mid_state() did not converge\n");
	    verbose_print_state("left",left);
	    verbose_print_state("right",right);
	    (void) printf("pjump = %g\n",pjump);
	    (void) printf("pml = %g, pmr = %g\n",pml,pmr);
	    (void) printf("uml = %g, umr = %g\n",uml,umr);
	    (void) printf("ml = %g, mr = %g\n",ml,mr);
	    print_rsoln_wave("l_wave = ",l_wave,", ");
	    print_rsoln_wave("r_wave = ",r_wave,"\n");
	}

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	{
	    (void) printf("l_wave = %s, r_wave = %s\n",rsoln_wave_name(l_wave),
	                  rsoln_wave_name(r_wave));
	    (void) printf("midstate pressures = %g %g\n",pml,pmr);
	    (void) printf("midstate velocity = %g %g\n",uml,umr);
	}
#endif /* defined(DEBUG_W_SPEED) */

	if (is_backward_wave(w_type))
	{
	    state_behind_sound_wave(left,ansr,NULL,Wx,pjump,ml,uml,
	                            pml,state_type,w_type,l_wave,LEFT_FAMILY);
#if defined(COMBUSTION_CODE)
	    if (r_wave == SHOCK || r_wave == RAREFACTION)
	        Set_params(ansr,left);
#endif /* defined(COMBUSTION_CODE) */
#if DONT_COMPILE /*OLD VERSION*/
	    if ((l_wave==RAREFACTION) && is_rarefaction_trailing_edge(w_type))
#endif /*DONT_COMPILE*/ /*OLD VERSION*/
	    if (is_rarefaction_trailing_edge(w_type))
	        copy_state(ansl,ansr);
	    else
	        set_state(ansl,state_type,sl);
	}

	else if (is_scalar_wave(w_type)) 
	{
	    midstate(left,ansl,ml,uml,pml,state_type,l_wave,LEFT_FAMILY);
	    midstate(right,ansr,mr,umr,pmr,state_type,r_wave,RIGHT_FAMILY);
	    cl = sound_speed(ansl);
	    cr = sound_speed(ansr);
	    *Wx = ((cl/(cl+cr))*uml + (cr/(cl+cr))*umr);
/*Previous Code:    *Wx = 0.5*(uml+umr); */
/* Evaluate cl = sound_speed(ansl); cr = sound_speed(ansr); */
/*TM - Try *Wx = ((cl/(cl+cr))*uml + (cr/(cl+cr))*umr); */
	}

	else if (is_forward_wave(w_type)) 
	{
	    state_behind_sound_wave(right,ansl,NULL,Wx,pjump,mr,umr,
	                            pmr,state_type,w_type,r_wave,RIGHT_FAMILY);
#if defined(COMBUSTION_CODE)
	    if (r_wave == SHOCK || r_wave == RAREFACTION)
	        Set_params(ansl,right);
#endif /* defined(COMBUSTION_CODE) */

#if DONT_COMPILE /*OLD VERSION*/
	    if ((r_wave==RAREFACTION) && is_rarefaction_trailing_edge(w_type))
#endif /*DONT_COMPILE*/ /*OLD VERSION*/
	    if (is_rarefaction_trailing_edge(w_type))
	        copy_state(ansr,ansl);
	    else
	        set_state(ansr,state_type,sr);
	}
	else
	{
	    screen("ERROR: unknown wave_type %d in w_speed\n",w_type);
	    clean_up(ERROR);
	}

#if defined(ONED)
	if ((dim == 1) && (state_type != TGAS_STATE))
	{
	    set_state(left,state_type,left);
	    set_state(right,state_type,right);
	}
#endif /* defined(ONED) */

	Left_w_speed("onedw_speed",ansl,ansr,Wx,1);
	return;
}	/*end onedw_speed*/

/*ARGSUSED*/
EXPORT	void state_behind_sound_wave(
	Locstate	         ahead,
	Locstate	         ans,
	double		         *cans,
	double		         *Wx,
	double		         pjump,
	double		         m,
	double		         um,
	double		         pm,
	int		         st_type_ans,
	int		         w_type,
	RIEMANN_SOLVER_WAVE_TYPE wave,
	int		         family)
{
	double	s, c;
	double	sgn = (family == LEFT_FAMILY) ? -1.0 : 1.0;
	double   pa;
	int	i, dim = Params(ahead)->dim;

#if defined(DEBUG_W_SPEED)
	debug_print("w_speed","Entered state_behind_sound_wave()\n");
	if (debug_w_speed == YES)
	{
	    if (!debugging("w_speed"))
	        (void) printf("Entered state_behind_sound_wave()\n");
	    switch (wave) 
	    {
	    case SHOCK:
	        (void) printf("shock\n");
	        break;

	    case RAREFACTION:
	        (void) printf("rarefaction\n");
	        break;

#if defined(COMBUSTION_CODE)
	    case STRONG_DET:
	        (void) printf("strong detonation\n");
	        break;

	    case CJ_DET:	
	        (void) printf("CJ-detonation\n");
	        (void) printf("CJ pressure = %g \n",pressure(ans));
	        (void) printf("CJ density = %g \n",Dens(ans));
	        (void) printf("CJ velocity = %g \n",um);
	        break;
#endif /* defined(COMBUSTION_CODE) */
	    
	    default:
	        screen("ERROR - unknown wave %d "
	               "in state_behind_sound_wave()\n",wave);
	        clean_up(ERROR);
	        break;
	    }
	}
#endif /* defined(DEBUG_W_SPEED) */
#if defined(COMBUSTION_CODE)
	if ((wave == CJ_DET) || (wave == STRONG_DET))
	{
	    CJ_state_behind_curved_wave(ahead,ans,pjump,Wx,st_type_ans,family);
	    return;
	}
#endif /* defined(COMBUSTION_CODE) */

	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = 0.0;
	switch (wave) 
	{
	case SHOCK:
	    s = Vel(ahead)[0] + sgn*m/Dens(ahead);
	    Dens(ans) = -sgn*m/(um - s);	
	    pa = pressure(ahead);
	    if (invalid_shock(pa,Dens(ahead),pm,Dens(ans)))
		Dens(ans) = dens_Hugoniot(pm,ahead);
            /* Mass fraction does not change across shock waves */
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                int    k;
                if(Params(ahead)->n_comps != 1)
                {
                    for(k = 0; k < Params(ahead)->n_comps; k++)
                        pdens(ans)[k] = (pdens(ahead)[k]/Dens(ahead))*Dens(ans);                }
            }
	    Set_params(ans,ahead);
	    switch (st_type_ans)
	    {
	    case EGAS_STATE:
	        set_type_of_state(ans,EGAS_STATE);
	        Vel(ans)[0] = um;
	        Energy(ans) = specific_internal_energy(ahead) +
	                      0.5*(pa+pm)*(1.0/Dens(ahead) - 1.0/Dens(ans));
#if defined(COMBUSTION_CODE)
	        if (Composition_type(ans) == ZND)
	            React(ans) = React(ahead);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	        break;
	    case GAS_STATE:
	        set_type_of_state(ans,GAS_STATE);
	        Mom(ans)[0] = Dens(ans)*um;
	        Energy(ans) = Dens(ans)*(0.5*sqr(um) +
	                      specific_internal_energy(ahead) + 
	                      0.5*(pa+pm)*(1.0/Dens(ahead) - 1.0/Dens(ans)));
#if defined(COMBUSTION_CODE)
	        if (Composition_type(ans) == ZND)
	            Prod(ans) = Prod(ahead);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	        break;
	    case TGAS_STATE:
	        set_type_of_state(ans,TGAS_STATE);
	        Vel(ans)[0] = um;
	        Press(ans) = pm;
#if defined(COMBUSTION_CODE)
	        if (Composition_type(ans) == ZND)
	            React(ans) = React(ahead);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	        break;
	    default:
	        set_type_of_state(ans,TGAS_STATE);
	        Vel(ans)[0] = um;
	        Press(ans) = pm;
#if defined(COMBUSTION_CODE)
	        if (Composition_type(ans) == ZND)
	            React(ans) = React(ahead);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	        set_state(ans,st_type_ans,ans);
	        break;
	    }
	    if (cans != NULL)
	        c = sound_speed(ans);
	    break;

	case RAREFACTION:
	    if (is_rarefaction_trailing_edge(w_type))
	    {
	        switch (st_type_ans)
	        {
	        case EGAS_STATE:
	            state_on_adiabat_with_pr(ahead,pm,ans,EGAS_STATE);
	            Vel(ans)[0] = um;
	            break;
	        case GAS_STATE:
	            state_on_adiabat_with_pr(ahead,pm,ans,GAS_STATE);
	            Mom(ans)[0] = Dens(ans)*um;
	            Energy(ans) += 0.5*Mom(ans)[0]*um;
		    break;
	        case TGAS_STATE:
	            state_on_adiabat_with_pr(ahead,pm,ans,TGAS_STATE);
	            Vel(ans)[0] = um;
		    break;
		default:
	            state_on_adiabat_with_pr(ahead,pm,ans,TGAS_STATE);
	            Vel(ans)[0] = um;
	            set_state(ans,st_type_ans,ans);
		    break;
	        }
		reset_gamma(ans);
	        c = sound_speed(ans);
	        s = um + sgn*c;
	    }
	    else
	    {
	        Dens(ans) = Dens(ahead);
                /* Mass fraction does not change across rarefaction waves */
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    int    k;
                    if(Params(ahead)->n_comps != 1)
                    {
                        for(k = 0; k < Params(ahead)->n_comps; k++)
                            pdens(ans)[k] = (pdens(ahead)[k]/Dens(ahead))*Dens(ans);
                    }
                }
	        Set_params(ans,ahead);
	        switch (st_type_ans)
	        {
	        case EGAS_STATE:
	            set_type_of_state(ans,EGAS_STATE);
	            Vel(ans)[0] = vel(0,ahead);
	            Energy(ans) = specific_internal_energy(ahead);
	            break;
	        case GAS_STATE:
	            set_type_of_state(ans,GAS_STATE);
	            Mom(ans)[0] = mom(0,ahead);
	            Energy(ans) = energy(ahead);
	            break;
	        case TGAS_STATE:
	            set_type_of_state(ans,TGAS_STATE);
	            Vel(ans)[0] = vel(0,ahead);
	            Press(ans) = pressure(ahead);
	            break;
		default:
	            set_type_of_state(ans,TGAS_STATE);
	            Vel(ans)[0] = vel(0,ahead);
	            Press(ans) = pressure(ahead);
	            set_state(ans,st_type_ans,ans);
	            break;
	        }
		reset_gamma(ans);
	        c = sound_speed(ahead);
	        s = vel(0,ahead) + sgn*c;
	    }
	    break;

#if defined(COMBUSTION_CODE)
	case STRONG_DET:
	    Set_other_params(ans,ahead);
	    s = Vel(ahead)[0] + m/Dens(ahead);
	    Dens(ans) = -sgn*m/(um - s);
	    switch (st_type_ans)
	    {
	    case EGAS_STATE:
	        set_type_of_state(ans,EGAS_STATE);
	        Vel(ans)[0] = um;
	        Energy(ans) = specific_internal_energy(ahead) +
	                      0.5*(pressure(ahead)+pm)*
	                      (1.0/Dens(ahead) - 1.0/Dens(ans));
	        break;
	    case GAS_STATE:
	        set_type_of_state(ans,GAS_STATE);
	        Mom(ans)[0] = Dens(ans)*um;
	        /* Is this correct for Stong Det?*/
	        Energy(ans) = Dens(ans)*(0.5*sqr(um) +
	                      specific_internal_energy(ahead) + 
	                      0.5*(pressure(ahead)+pm)*
	                      (1.0/Dens(ahead) - 1.0/Dens(ans)));
	    case TGAS_STATE:
	        set_type_of_state(ans,TGAS_STATE);
	        Vel(ans)[0] = um;
	        Press(ans) = pm;
	        break;
	    default:
	        set_type_of_state(ans,TGAS_STATE);
	        Vel(ans)[0] = um;
	        Press(ans) = pm;
	        set_state(ans,st_type_ans,ans);
	        break;
	    }
		reset_gamma(ans);
	    if (cans != NULL)
	        c = sound_speed(ans);
	    break;

	case CJ_DET:	
	    /*Is this right? */
	    s = CJ_det(ans,st_type_ans,ahead,family);
	    Set_other_params(ans,ahead);
	    /*
	    c = sound_speed(ans);
	    s = vel(0,ans) + sgn*c;
	    */
	    break;
#endif /* defined(COMBUSTION_CODE) */
	
	default:
	    screen("ERROR - unknown wave %d "
	           "in state_behind_sound_wave()\n",wave);
	    clean_up(ERROR);
	    break;
	}
	reset_gamma(ans);

	if (cans != NULL)
	    *cans = c;
	*Wx = s;

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	{
	    verbose_print_state("Answer from state_behind_sound_wave",ans);
	    (void) printf("Wave speed = %g, sound speed = %g\n",
	                  s,sound_speed(ans));

	}
	if ((debug_w_speed == YES) && (!debugging("w_speed")))
	    (void) printf("Left state_behind_sound_wave()\n");
	debug_print("w_speed","Left state_behind_sound_wave()\n");
#endif /* defined(DEBUG_W_SPEED) */
	return;
}	/*end state_behind_sound_wave*/

LOCAL	void midstate(
	Locstate	         ahead,
	Locstate	         ans,
	double		         m,
	double		         um,
	double		         pm,
	int		         st_type_ans,
	RIEMANN_SOLVER_WAVE_TYPE wave,
	int		         family)
{
	double sgn, pa;
#if defined(COMBUSTION_CODE)
	static Locstate CJ = NULL;
#endif /* defined(COMBUSTION_CODE) */
	

	sgn = (family == LEFT_FAMILY) ? 1.0 : -1.0;
	switch (wave) 
	{
	case SHOCK:
	    pa = pressure(ahead);
	    Set_params(ans,ahead);
	    set_type_of_state(ans,st_type_ans);
	    Dens(ans) = sgn*m/(um - vel(0,ahead) + sgn*m/Dens(ahead));
	    if (invalid_shock(pa,Dens(ahead),pm,Dens(ans)))
		Dens(ans) = dens_Hugoniot(pm,ahead);
            /* scaling */
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                int  i;
                if(Params(ahead)->n_comps != 1)
                {
                    for(i = 0; i < Params(ahead)->n_comps; i++)
                        pdens(ans)[i] = (pdens(ahead)[i]/Dens(ahead))*Dens(ans);                }
            }
	    switch (st_type_ans)
	    {
	    case TGAS_STATE:
	        Vel(ans)[0] = um;
	        Press(ans) = pm;
#if defined(COMBUSTION_CODE)
	        if (Composition_type(ans) == ZND)
	            React(ans) = React(ahead);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	        break;
	    case EGAS_STATE:
	        Vel(ans)[0] = um;
	        Energy(ans) = specific_internal_energy(ahead) +
	                      0.5*(pa+pm)*(1.0/Dens(ahead) - 1.0/Dens(ans));
#if defined(COMBUSTION_CODE)
	        if (Composition_type(ans) == ZND)
	            React(ans) = React(ahead);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	        break;
	    case GAS_STATE:
	        Mom(ans)[0] = Dens(ans)*um;
	        Energy(ans) = Dens(ans)*(0.5*sqr(um) +
	                      specific_internal_energy(ahead) + 
	                      0.5*(pa+pm)*(1.0/Dens(ahead) - 1.0/Dens(ans)));
#if defined(COMBUSTION_CODE)
	        if (Composition_type(ans) == ZND)
	            Prod(ans) = Prod(ahead);
#endif /* defined(COMBUSTION_CODE) */
		reset_gamma(ans);
	    }
	    break;

	case RAREFACTION:
	    state_on_adiabat_with_pr(ahead,pm,ans,st_type_ans);
	    if (pressure_is_near_vacuum(pm,ahead))
	    {
	        um = vel(0,ahead) - sgn*riemann_wave_curve(ahead,pm);
	    }
	    switch (st_type_ans)
	    {
	    case TGAS_STATE:
	    case EGAS_STATE:
	        Vel(ans)[0] = um;
	        break;
	    case GAS_STATE:
	        Mom(ans)[0] = Dens(ans)*um;
	        Energy(ans) += 0.5*Mom(ans)[0]*um;
		reset_gamma(ans);
	    }
	    break;

#if defined(COMBUSTION_CODE)
	case STRONG_DET:
	    pa = pressure(ahead);
	    Set_other_params(ans,ahead);
	    Dens(ans) = sgn*m/(um - Vel(ahead)[0] + sgn*m/Dens(ahead));
	    set_type_of_state(ans,st_type_ans);
	    switch (st_type_ans)
	    {
	    case TGAS_STATE:
	        Vel(ans)[0] = um;
	        Press(ans) = pm;
		reset_gamma(ans);
	        break;
	    case EGAS_STATE:
	        Vel(ans)[0] = um;
                      Energy(ans) = specific_internal_energy(ahead) +
	                          0.5*(pa+pm)*(1.0/Dens(ahead) - 1.0/Dens(ans));
		reset_gamma(ans);
	        break;
	    case GAS_STATE:
	        Mom(ans)[0] = Dens(ans)*um;
	        Energy(ans) = Dens(ans)*(0.5*sqr(um) +
	                      specific_internal_energy(ahead) + 
	                      0.5*(pa+pm)*(1.0/Dens(ahead) - 1.0/Dens(ans)));
		reset_gamma(ans);
	    }
	    break;

	case CJ_DET:
	    if (CJ == NULL)
	    {
	        (*Params(ahead)->_alloc_state)(&CJ,Params(ahead)->sizest);
	    }
	    CJ_det(CJ,st_type_ans,ahead,family);
	    Set_other_params(ans,ahead);
	    Set_params(CJ,ans);
#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES)
	    {
	        (void) printf("CJ pressure = %g \n",Press(CJ));
	        (void) printf("CJ density = %g \n",Dens(CJ));
	        (void) printf("CJ velocity = %g \n",Vel(CJ)[0]);
	    }
#endif /* defined(DEBUG_W_SPEED) */
	    state_on_adiabat_with_pr(CJ,pm,ans,st_type_ans);
	    if (st_type_ans == TGAS_STATE || st_type_ans == EGAS_STATE)
	        Vel(ans)[0] = um;
	    else if (st_type_ans == GAS_STATE)
	    {
	        Mom(ans)[0] = Dens(ans)*um;
	        Energy(ans) += 0.5*Mom(ans)[0]*um;
	    }
	    reset_gamma(ans);
	    break;
#endif /* defined(COMBUSTION_CODE) */

	default:
	    screen("ERROR: unknown wave %d in midstate()\n",wave);
	    clean_up(ERROR);
	    break;
	}
}	/*end midstate*/





/*
*			g_npt_w_speed():
*
*	Finds the states and wave speed for the solution of a non-local Riemann
*	problem.  Source-term corrections to front states and (in the case of 
*	cylindrical symmetry) to the method of characteristics are included.
*
*
*	Input
*	sten		- WSSten structure containing all input data
*
*	Output
*	ansl, ansr	- left/right propagated states
*	W		- wave velocity
*	
*/

EXPORT	void g_npt_w_speed(
	WSSten		*sten,
	Locstate	ansl,
	Locstate	ansr,
	double		*W)
{
	double		*pt = sten->coords;
	Locstate	 sll;
	Locstate	  sl;
	Locstate	  sr;
	Locstate	 srr;
	double		*nor = sten->nor;
	double		dn = sten->dn;
	double		dt = sten->dt;
	double		pjump = sten->pjump;
	int		w_type = sten->w_type;
	Front		*front = sten->front;
	RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
	RECT_GRID	*gr = front->rect_grid;
	INTERFACE	*intfc = front->interf;
	Gas_param	**params;
	NptWSpeedOpts	Opts;
	int		i, dim = front->interf->dim;
	size_t		sizest = front->sizest;
	double		rr,rl,r;		/* densities */
	double		pml, pmr, pl, pr;	/* pressures */
	double		vtanr[MAXD], vtanl[MAXD], V[MAXD];     /* velocities */
	double		ul,ur,uml,umr,vl,vr;	/* normal velocities */
	double		s;			/* wave speed	*/
	double		cl,cr,cm;		/* sound speeds */
	double		dn1,dn2,dn3;		/* origins of characteristics
	                                           in mesh units */
	double		mr,ml;			/* miscellaneous variables */
	double           dnl, dnr;
	double           difful,diffur,dterm; /*diffusion term*/
	static double	g[MAXD];
	static Locstate left = NULL, right = NULL, mid = NULL,
	                mid_l = NULL, mid_r = NULL, Tsl = NULL, Tsr = NULL;
	static Locstate	st_l1 = NULL, st_l2 = NULL, st_l3 = NULL,
	                st_r1 = NULL, st_r2 = NULL, st_r3 = NULL;
#if defined(COMBUSTION_CODE)
	static Locstate CJ = NULL;
#endif /* defined(COMBUSTION_CODE) */

	Opts = G3ptOpts;

#if defined(DEBUG_W_SPEED)
	debug_print("w_speed","Entered g_npt_w_speed()\n");
	if (debugging("w_speed"))
	    debug_w_speed = YES;
	else if (debug_w_speed == YES)
	    (void) printf("Entered g_npt_w_speed()\n");
#endif /* defined(DEBUG_W_SPEED) */

	if (Tsr == NULL) 
	{
	    alloc_state(front->interf,&st_l1,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&st_l2,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&st_l3,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&st_r1,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&st_r2,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&st_r3,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&mid,sizest);
	    alloc_state(front->interf,&mid_l,sizest);
	    alloc_state(front->interf,&mid_r,sizest);
	    alloc_state(front->interf,&left,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&right,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&Tsl,max(sizeof(VGas),sizest));
	    alloc_state(front->interf,&Tsr,max(sizeof(VGas),sizest));

#if defined(COMBUSTION_CODE)
	    alloc_state(front->interf,&CJ,sizest);
#endif /* defined(COMBUSTION_CODE) */

	}
	eval_gravity(pt,front->time,g);

	set_ws_slopes(sten);

	sll  =  sten->sl[1];
	sl   =  sten->sl[0];
	sr   =  sten->sr[0];
	srr  =  sten->sr[1];

	/* Get NS mass diffusion term */
	params = gas_params_list(intfc);
	dterm = params[0]->avisc.diffusivity_coef;
	
#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	{
	    (void) printf("INPUT DATA\n");
	    PrintWSSten(sten);
	    (void) printf("END INPUT DATA\n");
            print_g_npt_w_speed_opts(&Opts);
	}
#endif /* defined(DEBUG_W_SPEED) */

	if (w_type < FIRST_PHYSICS_WAVE_TYPE)
	{
	    bdry_npt_w_speed(sten,ansl,ansr,&Opts);
	    if (w_type == MOVABLE_BODY_BOUNDARY)
	    {
		HYPER_SURF *hs = sten->hs;
		double omega_dt,crds_com[MAXD];
		omega_dt = angular_velo(hs)*dt;
	    	for (i = 0; i < dim; ++i)
		{
		    W[i] = center_of_mass_velo(hs)[i];
		    crds_com[i] = sten->coords[i] - center_of_mass(hs)[i];
		}
		W[0] += -angular_velo(hs)*crds_com[1]*cos(omega_dt) -
			 angular_velo(hs)*crds_com[0]*sin(omega_dt);
		W[1] +=  angular_velo(hs)*crds_com[0]*cos(omega_dt) -
			 angular_velo(hs)*crds_com[1]*sin(omega_dt);
	    }
	    else
	    {
	    	for (i = 0; i < dim; ++i)
	            W[i] = 0.0;
	    }
	    Left_w_speed("g_npt_w_speed",ansl,ansr,W,dim);
	    return;
	}

	/* Calculate diffusion source terms */
	
	if (w_type == CONTACT)
	{
            difful = 3*dterm*(Dens(srr) - 2*Dens(sl) + 
	    		Dens(sll))/(dn*dn*Dens(sl));
            diffur = 3*dterm*(Dens(srr) - 2*Dens(sr) + 
	    		Dens(sll))/(dn*dn*Dens(sr));
	}
	else
	{
	    difful = diffur = 0.0;
	}

	    /* Set left,right state variables */

	set_type_of_state(ansl,GAS_STATE);
	set_type_of_state(ansr,GAS_STATE);
	set_state_for_find_mid_state(right,sr);
	set_state_for_find_mid_state(left,sl);
	copy_state(Tsr,right);
	copy_state(Tsl,left);

	rr  = Dens(right);		rl  = Dens(left);
	ur = ul = 0.0;
	for ( i = 0; i < dim; ++i)
	{
	    ur += nor[i] * Vel(right)[i];
	    ul += nor[i] * Vel(left)[i];
	}
	for ( i = 0; i < dim; ++i)
	{
	    vtanr[i] = Vel(right)[i] - nor[i] * ur;
	    vtanl[i] = Vel(left)[i]  - nor[i] * ul;
	}

	Vel(right)[0] = ur;	Vel(left)[0] = ul;
	for ( i = 1; i < dim; ++i)
	    Vel(right)[i] = Vel(left)[i] = Vel(mid)[i] = 0.0;

	if (find_mid_state(left,right,pjump,&pml,&pmr,&uml,&umr,&ml,&mr,
	                       &l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	{
	    (void) printf("WARNING in g_npt_w_speed(), find_mid_state() "
	                  "did not converge for left/right RP\n");
	}

	if(debugging("tst_nan"))
	{
	    printf("#tst_nan\n");
	    printf("%24.16e  %24.16e\n", pml, pmr);
	    printf("%24.16e  %24.16e\n", uml, umr);
	    printf("%24.16e  %24.16e\n", ml, mr);
	    printf("%d  %d\n", l_wave, r_wave);
	}
	
	if (strong_interaction(left,right,w_type,pml,pmr,uml,umr) == YES)
	{
#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES)
	    {
	        (void) printf("WARNING in g_npt_w_speed(), "
	                      "strong interaction detected\n");
	    }
#endif /* defined(DEBUG_W_SPEED) */
	    Opts.vector_moc = RIEMANN;
	    Opts.scalar_moc = RIEMANN;
	}
#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	{
	    (void) printf("Answer from find_mid_state()\n"
	                  "pml = %g, pmr = %g\n"
	                  "uml = %g, umr = %g\n"
	                  "ml = %g, mr = %g\n",
	                  pml,pmr,uml,umr,ml,mr);
	    print_rsoln_wave("l_wave = ",l_wave,", ");
	    print_rsoln_wave("r_wave = ",r_wave,"\n");
	}
#endif /* defined(DEBUG_W_SPEED) */

	if (is_backward_wave(w_type))
	{
#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES)
	        (void) printf("Backward wave\n");
#endif /* defined(DEBUG_W_SPEED) */
	    state_behind_sound_wave(left,mid,&cm,&s,pjump,ml,uml,pml,
	                TGAS_STATE,w_type,l_wave,LEFT_FAMILY);

	    uml = Vel(mid)[0];
	    cl = sound_speed(Tsl);
	    for (i = 0; i < dim; ++i)
	        W[i] = nor[i] * s;

	    /* Find feet of characteristics from wave location */

	    dn2 = 1.0 + (s - ul) * dt / dn;
	    dn3 = 1.0 + (s - ul - cl) * dt / dn;

	    dn1 = (s + sound_speed(sll));
	    for (i = 0; i < dim; ++i)
	        dn1 -= (nor[i] * vel(i,sll));
	    dn1 = 1.0 + dn1 * dt / dn;
	    dn1 = min(dn1,1.0 + ((s - uml + cm) * dt / dn));
	    dn1 = max(dn1,dn2);
	    dn1 = min(dn1,1.0);

	    /* Find states ahead of shock at feet of characteristics */

	    ws_interpolate(st_l3,dn3,NEGATIVE_SIDE,TGAS_STATE,sten);
	    ws_interpolate(st_l2,dn2,NEGATIVE_SIDE,TGAS_STATE,sten);
	    ws_interpolate(st_l1,dn1,NEGATIVE_SIDE,TGAS_STATE,sten);

#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES) 
	    {
	        (void) printf("Interpolated states\n");
	        (void) printf("dn1 = %g, dn2 = %g, dn3 = %g\n",dn1,dn2,dn3);
	        verbose_print_state("st_l1",st_l1);
	        verbose_print_state("st_l2",st_l2);
	        verbose_print_state("st_l3",st_l3);
	    }
#endif /* defined(DEBUG_W_SPEED) */

	    if (l_wave == RAREFACTION)
	    {
	        if (w_type==BACKWARD_SOUND_WAVE_TE)
	        {
	            (*Opts.vector_ahead_state_moc)(pt,Tsl,st_l1,st_l2,st_l3,
						   left,-dn,1.0-dn1,1.0-dn2,
						   1.0-dn3,nor,W,YES,dt,front);
	            w_speed(pt,sr,srr,right,Tsr,V,pjump,nor,w_type,front);
	            set_state(left,GAS_STATE,left);
	            w_speed(pt,left,right,ansl,ansr,V,pjump,nor,w_type,front);
	            for (i = 0; i < dim; ++i)
	                W[i] = 0.5*(W[i] + V[i]);
	            copy_state(ansl,ansr);
	        }
		else
	        {
	            (*Opts.vector_ahead_state_moc)(pt,Tsl,st_l1,st_l2,st_l3,
						   left,-dn,1.0-dn1,1.0-dn2,
						   1.0-dn3,nor,W,YES,dt,front);


	            s = - sound_speed(left);
	            for (i = 0; i < dim; ++i)
	                s += nor[i]*Vel(left)[i];
	            for (i = 0; i < dim; ++i)
	            {
	                V[i] = nor[i] * s;
	                W[i] = 0.5*(W[i] + V[i]);
	            }

	            set_state(ansl,GAS_STATE,left);
	            copy_state(ansr,ansl);
	        }
	    }
	    else
	    {

	        /* Find state behind shock at foot of characteristic */
	        /*   Find foot of characteristic from wave location  */

	        dn1 = (s + sound_speed(srr));
	        for (i = 0; i < dim; ++i)
	            dn1 -= (nor[i] * vel(i,srr));
	        dn1 *= dt / dn;
	        dn1 = max(dn1,((s - umr + cm) * dt / dn));
	        dn1 = min(dn1,1.0);
	        dn1 = max(dn1,0.0);


	        ws_interpolate(st_r1,dn1,POSITIVE_SIDE,TGAS_STATE,sten);
#if defined(DEBUG_W_SPEED)
	        if (debug_w_speed == YES) 
	        {
	            (void) printf("Interpolated right state\n");
	            (void) printf("dn1 = %g\n",dn1);
	            verbose_print_state("st_r1",st_r1);
	        }
#endif /* defined(DEBUG_W_SPEED) */

	        switch (Opts.vector_moc)
	        {
	        case RIEMANN:
	            /* Find left state after time dt */
	            /* Solve characteristic PDEs to get left state */

	            (*Opts.vector_ahead_state_moc)(pt,Tsl,st_l1,st_l2,st_l3,
						   left,
						   -dn,1.0-dn1,1.0-dn2,1.0-dn3,
						   nor,W,NO,0.0,front);

	            set_state(ansl,GAS_STATE,left);

#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("left state after dt\n");
	                verbose_print_state("ansl",ansl);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	        /* Separate out back mode at foot of characteristic */
	        /* by solving Riemann problem between sr, s1 */

	            set_state_for_find_mid_state(st_r1,st_r1);
	            Vel(st_r1)[0] = scalar_product(Vel(st_r1),nor,dim);
	            Set_params(st_r1,sr);

	            if (find_mid_state(right,st_r1,0.0,&pl,&pr,&vl,&vr,&ml,&mr,
				       &l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	            {
	                (void) printf("WARNING in g_npt_w_speed(), "
	                              "find_mid_state() did not converge for "
	                              "backward wave right/st_r1 RP\n");
	            }

	            switch (l_wave) 
	            {
	            case SHOCK:
	                Set_params(st_r1,sr);
	                r = ml/(vl - ur + ml/rr);
	                if (invalid_shock(pressure(Tsr),Dens(Tsr),pl,r))
		            r = dens_Hugoniot(pl,Tsr);
	                break;

	            case RAREFACTION:
	                Set_params(st_r1,sr);
	                r = dens_rarefaction(pl,Tsr);
	                break;

#if defined(COMBUSTION_CODE)
	            case STRONG_DET:
	                Set_other_params(st_r1,sr);
	                r = ml/(vl - ur + ml/rr);
	                break;

	            case CJ_DET:
	                CJ_det(CJ,TGAS_STATE,right,LEFT_FAMILY);
	                Set_other_params(st_r1,sr);
	                r = Dens(CJ);
	                vl = Vel(CJ)[0];
	                break;
#endif /* defined(COMBUSTION_CODE) */
	        
	            default:
	                screen("ERROR: unknown wave ");
	                screen("%d in g_npt_w_speed\n",l_wave);
	                clean_up(ERROR);
	            }

	            Dens(st_r1) = r;
	            Press(st_r1) = pl;
	            set_type_of_state(st_r1,TGAS_STATE);
	            for (i = 0; i < dim; ++i)
	                Vel(st_r1)[i] = vtanr[i] + nor[i] * vl;
#if defined(COMBUSTION_CODE)
	            if (Composition_type(sr) == ZND)
	                React(st_r1) = React(Tsr);
#endif /* defined(COMBUSTION_CODE) */
		    reset_gamma(st_r1);
	        
	        /* Solve Riemann problem to find right state */

	            set_state(st_r1,GAS_STATE,st_r1);
	            copy_state(left,ansl);
	            w_speed(pt,left,st_r1,ansl,ansr,V,pjump,nor,w_type,front);
#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("States after w_speed\n");
	                verbose_print_state("left",left);
	                verbose_print_state("st_r1",st_r1);
	                verbose_print_state("ansl",ansl);
	                verbose_print_state("ansr",ansr);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	        /* Use centered difference in time for wave speed */

	            for (i = 0; i < dim; ++i)
	                W[i] = 0.5*(W[i] + V[i]);
                    include_source(pt,ansl,ul,dt,nor,W,g,gr,difful,w_type);
                    include_source(pt,ansr,uml,dt,nor,W,g,gr,diffur,w_type);
	            break;

	        case MOC_PLUS_RIEMANN:
	        /* Find left state after time dt */
	        /* Solve characteristic PDEs to get left state */

	            (*Opts.vector_ahead_state_moc)(pt,Tsl,st_l1,st_l2,st_l3,
						   left,
						   -dn,1.0-dn1,1.0-dn2,1.0-dn3,
						   nor,W,NO,0.0,front);

	            set_state(ansl,GAS_STATE,left);
	            set_state(st_r1,GAS_STATE,st_r1);

#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("left state after dt\n");
	                verbose_print_state("ansl",ansl);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	            /* Solve Riemann problem to find right state */

	            copy_state(left,ansl);
	            w_speed(pt,left,st_r1,ansl,ansr,V,pjump,nor,w_type,front);
#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("States after w_speed\n");
	                verbose_print_state("left ",left);
	                verbose_print_state("st_r1",st_r1);
	                verbose_print_state("ansl ",ansl);
	                verbose_print_state("ansr ",ansr);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	        /* Use centered difference in time for wave speed */

	            for (i = 0; i < dim; ++i)
	                W[i] = 0.5*(W[i] + V[i]);

	            include_source(pt,ansl,ul,dt,nor,W,g,gr,difful,w_type);
	            include_source(pt,ansr,uml,dt,nor,W,g,gr,diffur,w_type);
	            break;
	        
	        case MOC_PLUS_RH:

	        /* Find left state after time dt */
	        /* Solve characteristic PDEs to get left state */

	            (*Opts.vector_ahead_state_moc)(pt,Tsl,st_l1,st_l2,st_l3,
						   left,
						   -dn,1.0-dn1,1.0-dn2,1.0-dn3,
						   nor,W,YES,dt,front);

	            set_state(ansl,GAS_STATE,left);

#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("left state after dt\n");
	                verbose_print_state("ansl",ansl);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	    /* Use characteristic equations and two Riemann invariants */

	            set_state(st_r1,TGAS_STATE,st_r1);
	            if (!shock_moc_plus_rh(pt,left,mid,st_r1,ansr,
	                                   -dn1*dn,nor,W,w_type,front))
		    {
			MOC_TYPE vector_moc = G3ptOpts.vector_moc;
			G3ptOpts.vector_moc = RIEMANN;
			g_npt_w_speed(sten,ansl,ansr,W);
			G3ptOpts.vector_moc = vector_moc;
			return;
		    }
	        }
	    }    
	}
	else if (is_scalar_wave(w_type))
	{
#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES)
	        (void) printf("Contact\n");
#endif /* defined(DEBUG_W_SPEED) */
	    set_type_of_state(mid_l,TGAS_STATE);
	    /* Use mid as temporary storage to compute sound speed on left */
	    Press(mid_l) = pml;
	    Vel(mid_l)[0] = uml;
	    switch (l_wave) 
	    {
	    case SHOCK:
	        Set_params(mid_l,sl);
	        Dens(mid_l) = ml/(uml - ul + ml/rl);
	        if (invalid_shock(pressure(Tsl),Dens(Tsl),pml,Dens(mid_l)))
		    Dens(mid_l) = dens_Hugoniot(pml,Tsl);
                /* Mass fraction does not change if not accross contact */
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    if(Params(sl)->n_comps != 1)
                    {
                        for(i = 0; i < Params(sl)->n_comps; i++)
                            pdens(mid_l)[i] = (pdens(sl)[i]/Dens(sl))*Dens(mid_l);
                    }
                }
		reset_gamma(mid_l);
	        break;

	    case RAREFACTION:
	        Set_params(mid_l,sl);
	        Dens(mid_l)= dens_rarefaction(Press(mid_l),Tsl);
		if(debugging("tst_nan"))
			printf("mid_l %24.16e\n", Dens(mid_l));

                /* Mass fraction does not change if not accross contact */
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    if(Params(Tsl)->n_comps != 1)
                    {
                        for(i = 0; i < Params(Tsl)->n_comps; i++)
                            pdens(mid_l)[i] = (pdens(Tsl)[i]/Dens(Tsl))*Dens(mid_l);
                    }
                }
	        if (pressure_is_near_vacuum(pml,left))
	            Vel(mid_l)[0] = ul - riemann_wave_curve(left,Press(mid_l));
		reset_gamma(mid_l);
	        break;

#if defined(COMBUSTION_CODE)
	    case STRONG_DET:
	        Set_other_params(mid_l,sl);
	        Dens(mid_l) = ml/(uml - ul + ml/rl);
		reset_gamma(mid_l);
	        break;

	    case CJ_DET:
	        CJ_det(CJ,TGAS_STATE,left,LEFT_FAMILY);
	        Set_other_params(mid_l,sl);
	        Set_params(CJ,mid_l);
	        Dens(mid_l) = dens_rarefaction(Press(mid_l),CJ);
		reset_gamma(mid_l);
	        break;
#endif /* defined(COMBUSTION_CODE) */
	    
	    default:
	        screen("ERROR: unknown wave %d in w_speed\n",l_wave);
	        clean_up(ERROR);
	        break;
	    }
#if defined(COMBUSTION_CODE)
	    if (Composition_type(sl) == ZND)
	        React(mid_l) = React(Tsl); 
#endif /* defined(COMBUSTION_CODE) */
	    cl = sound_speed(mid_l);

	    set_type_of_state(mid_r,TGAS_STATE);
	    Press(mid_r) = pmr;	/* Store pm in mid */
	    Vel(mid_r)[0] = umr;
	    switch (r_wave) 
	    {
	    case SHOCK:
	        Set_params(mid_r,sr);
	        Dens(mid_r) = -mr/(umr - ur - mr/rr);
	        if (invalid_shock(pressure(Tsr),Dens(Tsr),pmr,Dens(mid_r)))
		    Dens(mid_r) = dens_Hugoniot(pmr,Tsr);
		
		if(debugging("tst_nan"))
			printf("mid_r %24.16e\n", Dens(mid_r));

                /* Mass fraction does not change if not accross contact */
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    if(Params(sr)->n_comps != 1)
                    {
                        for(i = 0; i < Params(sr)->n_comps; i++)
                            pdens(mid_r)[i] = (pdens(sr)[i]/Dens(sr))*Dens(mid_r);
                    }
                }
                reset_gamma(mid_r);
	        break;

	    case RAREFACTION:
	        Set_params(mid_r,sr);
	        Dens(mid_r)= dens_rarefaction(pmr,Tsr);

                /* Mass fraction does not change if not accross contact */
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    if(Params(Tsr)->n_comps != 1)
                    {
                        for(i = 0; i < Params(Tsr)->n_comps; i++)
                            pdens(mid_r)[i] = (pdens(Tsr)[i]/Dens(Tsr))*Dens(mid_r);
                    }
                }
	        if (pressure_is_near_vacuum(pmr,right))
	            Vel(mid_r)[0] = ur + riemann_wave_curve(right,Press(mid_r));
                reset_gamma(mid_r);
	        break;

#if defined(COMBUSTION_CODE)
	    case STRONG_DET:
	        Set_other_params(mid_r,sr);
	        Dens(mid_r) = -mr/(umr - ur - mr/rr);
                reset_gamma(mid_r);
	        break;

	    case CJ_DET:
	        CJ_det(CJ,TGAS_STATE,right,RIGHT_FAMILY);
	        Set_other_params(mid_r,sr);
	        Set_params(CJ,mid_r);
	        Dens(mid_r) = dens_rarefaction(Press(mid_r),CJ);
                reset_gamma(mid_r);
	        break;
#endif /* defined(COMBUSTION_CODE) */

	    default:
	        screen("ERROR: unknown wave %d in w_speed\n",r_wave);
	        clean_up(ERROR);
	        break;
	    }
	    Set_params(mid_r,sr);
#if defined(COMBUSTION_CODE)
	    if (Composition_type(sl) == ZND)
	        React(mid_r) = React(Tsr); 
#endif /* defined(COMBUSTION_CODE) */
	    cr = sound_speed(mid_r);

	    for (i = 0; i < dim; ++i)
		 W[i] = nor[i]*((cl/(cl+cr))*uml + (cr/(cl+cr))*umr);
/*Previous Code: W[i] = nor[i] * 0.5 * (uml + umr);  -TOM */
/*TM - Try W[i] = nor[i]*((cl/(cl+cr))*uml + (cr/(cl+cr))*umr); */
	    
	    /* Find right and left states after time dt */

	    if(debugging("tst_nan"))
		    printf("#cl cr  %24.16e  %24.16e\n", cl, cr);

	    /* Find far left and right states */
	    IncomingStatesAtContact(sten,cl,cr,st_l3,st_r1,&dnl,&dnr,&Opts);

#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES) 
	    {
	        (void) printf("Interpolated states\n");
	        verbose_print_state("st_r1",st_r1);
	        verbose_print_state("st_l3",st_l3);
	    }
#endif /* defined(DEBUG_W_SPEED) */

#if !defined(COMBUSTION_CODE)
	    Set_params(st_r1,sr);    Set_params(st_l3,sl);
#endif /* !defined(COMBUSTION_CODE) */

	    switch (Opts.scalar_moc)
	    {
	    case RIEMANN:
                contact_riemann(left,Tsl,right,Tsr,st_l3,dnl,st_r1,dnr,ansl,
                                ansr,sten,W,&Opts,g,difful,diffur);
	        break;

	    case MOC_PLUS_RIEMANN:
	        set_state(st_r1,TGAS_STATE,st_r1);
	        set_state(st_l3,TGAS_STATE,st_l3);
	        /* Find corrected left and right states */
	        /* by solving characteristic equations */
	        contact_cheap_moc(pt,Tsl,Tsr,st_l3,st_r1,ansl,ansr,
				  cl,cr,pjump,dnl,dnr,nor,W,front);
	        break;


	    case MOC_PLUS_RH:
	        set_state(st_r1,TGAS_STATE,st_r1);
	        set_state(st_l3,TGAS_STATE,st_l3);
	        /* Use characteristic equations and */
	        /* two Riemann invariants */
	        contact_moc_plus_rh(pt,pmr,Tsl,Tsr,st_l3,st_r1,
	                            ansl,ansr,pjump,dnl,dnr,nor,W,front);
	        break;
	    }
	}
	else if (is_forward_wave(w_type)) 
	{
#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES)
	        (void) printf("Forward Wave\n");
#endif /* defined(DEBUG_W_SPEED) */
	    state_behind_sound_wave(right,mid,&cm,&s,pjump,mr,umr,pmr,
	                TGAS_STATE,w_type,r_wave,RIGHT_FAMILY);

	    umr = Vel(mid)[0];
	    cr = sound_speed(Tsr);

	    for (i = 0; i < dim; ++i)
	        W[i] = nor[i] * s;

	        /* Find feet of characteristics ahead */
	        /*  of the shock from wave location   */

	    dn1 = (s - ur + cr) * dt / dn;
	    dn2 = (s - ur) * dt / dn;

	    dn3 = s - sound_speed(srr);
	    for (i = 0; i < dim; ++i)
	        dn3 -= nor[i] * vel(i,srr);
	    dn3 *= dt / dn;
	    dn3 = max(dn3,((s - umr - cm) * dt / dn));
	    dn3 = min(dn3,dn2);
	    dn3 = max(dn3,0.0);

	    /* Find states ahead of the shock at feet of characteristics */

	    ws_interpolate(st_r1,dn1,POSITIVE_SIDE,TGAS_STATE,sten);
	    ws_interpolate(st_r2,dn2,POSITIVE_SIDE,TGAS_STATE,sten);
	    ws_interpolate(st_r3,dn3,POSITIVE_SIDE,TGAS_STATE,sten);

#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES) 
	    {
	        (void) printf("Interpolated states\n");
	        (void) printf("dn1 = %g, dn2 = %g, dn3 = %g\n",dn1,dn2,dn3);
	        verbose_print_state("st_r1",st_r1);
	        verbose_print_state("st_r2",st_r2);
	        verbose_print_state("st_r3",st_r3);
	    }
#endif /* defined(DEBUG_W_SPEED) */

	    if (r_wave==RAREFACTION)
	    {
	        if (w_type==FORWARD_SOUND_WAVE_TE)
	        {
	            (*Opts.vector_ahead_state_moc)(pt,Tsr,st_r1,st_r2,st_r3,
						   right,dn,dn1,dn2,dn3,
					           nor,W,YES,dt,front);

	            w_speed(pt,sll,sr,Tsl,left,V,pjump,nor,w_type,front);
	            set_state(right,GAS_STATE,right);
	            w_speed(pt,left,right,ansl,ansr,V,pjump,nor,w_type,front);
	            for (i = 0; i < dim; ++i)
	                W[i] = 0.5*(W[i] + V[i]);
	            copy_state(ansr,ansl);
	        }
		else
	        {
	            (*Opts.vector_ahead_state_moc)(pt,Tsr,st_r1,st_r2,st_r3,
						   right,dn,dn1,dn2,dn3,
					           nor,W,YES,dt,front);


	            s = sound_speed(right);
	            for (i = 0; i < dim; ++i)
	                s += nor[i] * Vel(right)[i];
	            for (i = 0; i < dim; ++i)
	            {
	                V[i] = nor[i] * s;
	                W[i] = 0.5*(W[i] + V[i]);
	            }

	            set_state(ansr,GAS_STATE,right);
	            copy_state(ansl,ansr);
	        }
	    }
	    else
	    {
	        /* Find foot of characteristic behind */
	        /*    shock from wave location      */

	        dn3 = s - sound_speed(sll);
	        for (i = 0; i < dim; ++i)
	            dn3 -= nor[i] * vel(i,sll);
	        dn3 = 1.0 + dn3 * dt / dn;
	        dn3 = min(dn3,(1.0 + (s - uml - cm) * dt / dn));
	        dn3 = min(dn3,1.0);
	        dn3 = max(dn3,0.0);

	        /* Find state behind shock at foot of characteristic */

	        ws_interpolate(st_l3,dn3,NEGATIVE_SIDE,TGAS_STATE,sten);
#if defined(DEBUG_W_SPEED)
	        if (debug_w_speed == YES) 
	        {
	            (void) printf("Left interpolated state\n""dn3 = %g\n",dn3);
	            verbose_print_state("st_l3",st_l3);
	        }
#endif /* defined(DEBUG_W_SPEED) */


	        switch (Opts.vector_moc)
	        {
	        case RIEMANN:

	            /* Find right state after time dt */
	            /* Solve characteristic PDEs for right state */

	            (*Opts.vector_ahead_state_moc)(pt,Tsr,st_r1,st_r2,st_r3,
	                                           right,dn,dn1,dn2,dn3,
						   nor,W,NO,0.0,front);
	            
	            set_state(ansr,GAS_STATE,right);

#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("right state after dt\n");
	                verbose_print_state("ansr",ansr);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	                /*   Separate out forward mode at foot  */
	                /* of characteristic by solving Riemann */
	                /*    problem between st_l3, sl       */

	            set_state_for_find_mid_state(st_l3,st_l3);
	            Vel(st_l3)[0] = scalar_product(Vel(st_l3),nor,dim);
	            Set_params(st_l3,sl);

	            if (find_mid_state(st_l3,left,0.0,&pl,&pr,&vl,&vr,&ml,&mr,
				       &l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	            {
	                (void) printf("WARNING in g_npt_w_speed(), "
	                              "find_mid_state() did not converge for "
	                              "forward wave st_l3/left RP\n");
	            }

	            switch (r_wave) 
	            {
	            case SHOCK:
	                r = -mr/(vr - ul - mr/rl);
	                if (invalid_shock(pressure(Tsl),Dens(Tsl),pr,r))
		            r = dens_Hugoniot(pr,Tsl);
	                Set_params(st_l3,sl);
	                break;

	            case RAREFACTION:
	                r = dens_rarefaction(pr,Tsl);
	                Set_params(st_l3,sl);
	                break;

#if defined(COMBUSTION_CODE)
	            case STRONG_DET:
	                Set_other_params(st_l3,sl);
	                r = -mr/(vr - ul - mr/rl);
	                break;

	            case CJ_DET:
	                CJ_det(CJ,TGAS_STATE,left,RIGHT_FAMILY);
	                Set_other_params(st_l3,sl);
	                r = Dens(CJ);
	                vr = Vel(CJ)[0];
	                break;
#endif /* defined(COMBUSTION_CODE) */
	            
	            default:
	                screen("ERROR: unknown wave ");
	                screen("%d in g_npt_w_speed\n",r_wave);
	                clean_up(ERROR);
	                break;
	            }

	            set_type_of_state(st_l3,TGAS_STATE);
	            Dens(st_l3) = r;
	            Press(st_l3) = pr;
	            for ( i = 0; i < dim; ++i)
	                Vel(st_l3)[i] = vtanl[i] + nor[i] * vr;
                    /* Mass fraction of st_l3 is preserved, scaling is performed */
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                        int    k;
                        double  sum;
                        if(Params(st_l3)->n_comps != 1)
                        {
                            sum = 0.0;
                            for(k = 0; k < Params(st_l3)->n_comps; k++)
                                sum += pdens(st_l3)[k];
                            for(k = 0; k < Params(st_l3)->n_comps; k++)
                                pdens(st_l3)[k] *= Dens(st_l3)/sum;
                        }
                    }
#if defined(COMBUSTION_CODE)
	            if (Composition_type(sl) == ZND)
	                React(st_l3) = React(Tsl);
#endif /* defined(COMBUSTION_CODE) */
		    reset_gamma(st_l3);
	            
	            /* Solve Riemann problem to find right state */

	            set_state(st_l3,GAS_STATE,st_l3);
	            copy_state(right,ansr);
	            w_speed(pt,st_l3,right,ansl,ansr,V,pjump,nor,w_type,front);
#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("States after w_speed\n");
	                verbose_print_state("st_l3",st_l3);
	                verbose_print_state("right",right);
	                verbose_print_state("ansl ",ansl);
	                verbose_print_state("ansr ",ansr);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	            /* Use centered in time difference for wave speed */

	            for ( i= 0; i < dim; ++i)
	                W[i] = 0.5*(W[i] + V[i]);
	            include_source(pt,ansl,umr,dt,nor,W,g,gr,difful,w_type);
	            include_source(pt,ansr,ur,dt,nor,W,g,gr,diffur,w_type);
	            break;

	        case MOC_PLUS_RIEMANN:

	            /* Find right state after time dt */
	            /* Solve characteristic PDEs for right state */

	            (*Opts.vector_ahead_state_moc)(pt,Tsr,st_r1,st_r2,st_r3,
						   right,dn,dn1,dn2,dn3,
						   nor,W,NO,0.0,front);
	            
	            set_state(ansr,GAS_STATE,right);
	            set_state(st_l3,GAS_STATE,st_l3);

#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("right state after dt\n");
	                verbose_print_state("ansr",ansr);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	            /* Solve Riemann problem to find left state */

	            copy_state(right,ansr);
	            w_speed(pt,st_l3,right,ansl,ansr,V,pjump,nor,w_type,front);

#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("States after w_speed\n");
	                verbose_print_state("st_l3",st_l3);
	                verbose_print_state("right",right);
	                verbose_print_state("ansl",ansl);
	                verbose_print_state("ansr",ansr);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	            /* Use centered in time difference for wave speed */

	            for (i = 0; i < dim; ++i)
	                W[i] = 0.5*(W[i] + V[i]);

	            include_source(pt,ansl,umr,dt,nor,W,g,gr,difful,w_type);
	            include_source(pt,ansr,ur,dt,nor,W,g,gr,diffur,w_type);
	            break;

	        case MOC_PLUS_RH:

	            /* Find right state after time dt */
	            /* Solve characteristic PDEs for right state */

	            (*Opts.vector_ahead_state_moc)(pt,Tsr,st_r1,st_r2,st_r3,
						   right,dn,dn1,dn2,dn3,
						   nor,W,YES,dt,front);
	            
	            set_state(ansr,GAS_STATE,right);

#if defined(DEBUG_W_SPEED)
	            if (debug_w_speed == YES) 
	            {
	                (void) printf("right state after dt\n");
	                verbose_print_state("ansr",ansr);
	            }
#endif /* defined(DEBUG_W_SPEED) */

	                /* Use characteristic equations and */
	                /*    two Riemann invariants        */

	            set_state(st_l3,TGAS_STATE,st_l3);
	            if (!shock_moc_plus_rh(pt,right,mid,st_l3,ansl,
				           (1.0 - dn3)*dn,nor,W,w_type,front))
		    {
			MOC_TYPE vector_moc = G3ptOpts.vector_moc;
			G3ptOpts.vector_moc = RIEMANN;
			g_npt_w_speed(sten,ansl,ansr,W);
			G3ptOpts.vector_moc = vector_moc;
			return;
		    }
	        }
	    }
	}
#if defined(COMBUSTION_CODE)
	else if (is_thinflame_wave(w_type))
    	{
            double cl, cr, w, ql, qr;
 		
            cr = sound_speed(sr);
	    cl = sound_speed(sl);
	       
	    if(Unburned(sr))
	    {
		w = flame_velocity(sr,sr);
			
		/* To satisfy the exothermic conditions ql-qr < 0 
		 * it should be ql = 0 and qr = Heat_release(st_deflg) 
		 * since in this situation right side is unburned. 
		 * But it doesn't matter to which side one adds the 
		 * heat as long as the sign in front of (ql-qr) is 
		 * chosen properly.
		*/
		ql = Heat_release(sl);
		qr = 0; 			

		s = ur + w;			/* Absolute flame velocity */
		for(i=0; i<dim; i++)
			W[i] = nor[i]*s;

		/* Find feet of characteristics. There are three 
		 * characteristics, two on positive side and one on 
		 * the negative side. The characteristic on the 
		 * positive side are u and u-c; on the negative side
		 * it is u+c. States on the feet of the characteristics
		 * are found through interpolation. 
                */
		dn1 = 1.0+(s-ul-cl)*dt/dn;
		dn2 = (s - ur)*dt/dn;
		dn3 = (s - ur + cr)*dt/dn;
		
		ws_interpolate(st_l1,dn1,NEGATIVE_SIDE,TGAS_STATE,sten);
		ws_interpolate(st_r2,dn2,POSITIVE_SIDE,TGAS_STATE,sten);
		ws_interpolate(st_r3,dn3,POSITIVE_SIDE,TGAS_STATE,sten);

#if defined(DEBUG_W_SPEED)
	    	if (debug_w_speed == YES) 
	    	{
		    (void) printf("Interpolated states\n");
		    (void) printf("dn1 = %g, dn2 = %g, dn3 = %g\n",dn1,dn2,dn3);
		    verbose_print_state("st_l1",st_l1);
		    verbose_print_state("st_r2",st_r2);
		    verbose_print_state("st_r3",st_r3);
	    	}
#endif /* defined(DEBUG_W_SPEED) */
	        switch (Opts.vector_moc)
	        {
		    /* Only flame_cheap_moc_plus_rh method is
		    *  implmented. This method finds the states
		    *  using characteristic equations and 
		    *  Rankine-Hugoniot conditions.
		    */
		    case RIEMANN:
			set_state(st_l1,TGAS_STATE,st_l1);
			set_state(st_r2,TGAS_STATE,st_r2);
			set_state(st_r3,TGAS_STATE,st_r3);
			flame_moc_plus_rh(pt,Tsl,Tsr,st_l1,st_r2,st_r3,
					          ansl,ansr,pjump,dn1,dn2,dn3,
                                                  dn,nor,W,cl,cr,g,ql,qr,front); 
			break;

		    case MOC_PLUS_RIEMANN:
			set_state(st_l1,TGAS_STATE,st_l1);
			set_state(st_r2,TGAS_STATE,st_r2);
			set_state(st_r3,TGAS_STATE,st_r3);
			flame_moc_plus_rh(pt,Tsl,Tsr,st_l1,st_r2,st_r3,
					          ansl,ansr,pjump,dn1,dn2,dn3,
                                                  dn,nor,W,cl,cr,g,ql,qr,front); 
			break;

		    case MOC_PLUS_RH:
			set_state(st_l1,TGAS_STATE,st_l1);
			set_state(st_r2,TGAS_STATE,st_r2);
			set_state(st_r3,TGAS_STATE,st_r3);
			flame_moc_plus_rh(pt,Tsl,Tsr,st_l1,st_r2,st_r3,
					          ansl,ansr,pjump,dn1,dn2,dn3,
                                                  dn,nor,W,cl,cr,g,ql,qr,front); 
			break;
		}
	    }	 
	    else
	    {
		w = flame_velocity(sl,sl);
			
		/* To satisfy the exothermic conditions ql-qr < 0 
		 * it should be ql = 0 and qr = Heat_release(st_deflg) 
		 * since in this situation right side is unburned. 
		 * But it doesn't matter to which side one adds the 
		 * heat as long as the sign in front of (ql-qr) is 
		 * chosen properly.
		*/
		qr = Heat_release(sr);
		ql = 0; 			

		s = ul - w;			/* Absolute flame velocity */
		for(i=0; i<dim; i++)
			W[i] = nor[i]*s;

		/* Find feet of characteristics. There are three 
		 * characteristics, two on positive side and one on 
		 * the negative side. The characteristic on the 
		 * positive side are u and u-c; on the negative side
		 * it is u+c. States on the feet of the characteristics
		 * are found through interpolation. 
                */
		dn1 = 1.0+(s-ul-cl)*dt/dn;
		dn2 = 1.0-w*dt/dn;
		dn3 = (s - ur + cr)*dt/dn;
		
		ws_interpolate(st_l1,dn1,NEGATIVE_SIDE,TGAS_STATE,sten);
		ws_interpolate(st_l2,dn2,NEGATIVE_SIDE,TGAS_STATE,sten);
		ws_interpolate(st_r3,dn3,POSITIVE_SIDE,TGAS_STATE,sten);


	        switch (Opts.vector_moc)
	        {
		    /* Only flame_cheap_moc_plus_rh method is
		    *  implmented. This method finds the states
		    *  using characteristic equations and 
		    *  Rankine-Hugoniot conditions.
		    */
		    case RIEMANN:
			set_state(st_l1,TGAS_STATE,st_l1);
			set_state(st_l2,TGAS_STATE,st_l2);
			set_state(st_r3,TGAS_STATE,st_r3);
			flame_moc_plus_rh(pt,Tsl,Tsr,st_l1,st_l2,st_r3,
					          ansl,ansr,pjump,dn1,dn2,dn3,
                                                  dn,nor,W,cl,cr,g,ql,qr,front); 
			break;

		    case MOC_PLUS_RIEMANN:
			set_state(st_l1,TGAS_STATE,st_l1);
			set_state(st_l2,TGAS_STATE,st_l2);
			set_state(st_r3,TGAS_STATE,st_r3);
			flame_moc_plus_rh(pt,Tsl,Tsr,st_l1,st_l2,st_r3,
					          ansl,ansr,pjump,dn1,dn2,dn3,
                                                  dn,nor,W,cl,cr,g,ql,qr,front); 
			break;

		    case MOC_PLUS_RH:
			set_state(st_l1,TGAS_STATE,st_l1);
			set_state(st_l2,TGAS_STATE,st_l2);
			set_state(st_r3,TGAS_STATE,st_r3);
			flame_moc_plus_rh(pt,Tsl,Tsr,st_l1,st_l2,st_r3,
					          ansl,ansr,pjump,dn1,dn2,dn3,
                                                  dn,nor,W,cl,cr,g,ql,qr,front); 
			break;
		}
	    }	 

	}
#endif /* defined(COMBUSTION_CODE) */
	else
	{
	    screen("ERROR: unknown wave_type %d in g_npt_w_speed\n",w_type);
	    clean_up(ERROR);
	}

	Left_w_speed("g_npt_w_speed",ansl,ansr,W,dim);
}	/*end g_npt_w_speed*/

/*#bjet2 */
EXPORT	void	tst_fprint_Tan_stencil(
	FILE		*file,
	Front		*front,
	Tan_stencil	*sten)
{
	int		i;
	int		dim = front->rect_grid->dim;
	int		nrad = sten->npts/2;
	char		s[80];

	(void) fprintf(file,"Data for Tan_stencil %p\n",sten);
	(void) fprintf(file,"npts = %d\n",sten->npts);
	for (i = -nrad; i <= nrad; ++i)
	{
	    (void) sprintf(s,"Coords(sten->p[%d]) = ",i);
	    fprint_general_vector(file,s,Coords(sten->p[i]),dim,"\n");
	    (void) fprintf(file,"sten->t[%d] = %g\n",i,sten->t[i]);
	    (void) fprintf(file,"sten->leftst[%d]\n",i);
	    verbose_print_state("leftst ", sten->leftst[i]);
	    /*fprint_state_data(file,sten->leftst[i],front->interf); */
	    (void) fprintf(file,"sten->rightst[%d]\n",i);
	    verbose_print_state("rightst ", sten->rightst[i]);
	    /*fprint_state_data(file,sten->rightst[i],front->interf); */
	}
}		/*end fprint_Tan_stencil*/

LOCAL	void bdry_npt_w_speed(
	WSSten		*wssten,
	Locstate	ansl,
	Locstate	ansr,
	NptWSpeedOpts   *opts)
{
	Locstate	 sll =  wssten->sl[1];
	Locstate	  sl =  wssten->sl[0];
	Locstate	  sr =  wssten->sr[0];
	Locstate	 srr =  wssten->sr[1];
	COMPONENT	l_comp = wssten->ncomp;
	COMPONENT	r_comp = wssten->pcomp;
	COMPONENT       comp;
	double		*nor = wssten->nor;
	double		dn = wssten->dn;
	double		dt = wssten->dt;
	double           W[3];
	int		w_type = wssten->w_type;
	Front		*front = wssten->front;
	Wave		*wave = wssten->wave;
	SIDE            side;
	HYPER_SURF	*hs = wssten->hs;
	double		*pt = wssten->coords;
	double           alpha;
	static	int	nrad;
	static	Tan_stencil *sten = NULL;
	static	Locstate left = NULL, right = NULL, Tsl = NULL, Tsr = NULL;
	static	Locstate st_l3 = NULL, st_r1 = NULL;
	static	Locstate stemp = NULL;
	static  Locstate ansl_1, ansl_2, ansr_1, ansr_2;

	double		ul, ur;		/* normal velocities */
	double		cl, cr;		/* sound speeds */
	double		dn1, dn3;	/* origins of characteristics
	                                   in mesh units */
	size_t		sizest = front->sizest;
	INTERFACE       *intfc = front->interf;
	int             i, j, k, dim = intfc->dim;

#if defined(DEBUG_W_SPEED)
	debug_print("w_speed","Entered bdry_npt_w_speed()\n");
	if (debugging("w_speed"))
	    debug_w_speed = YES;
	else if (debug_w_speed == YES)
	    (void) printf("Entered bdry_npt_w_speed()\n");
	if (debug_w_speed == YES)
	    print_wave_type("wave type = ",w_type,"\n",intfc);
#endif /* defined(DEBUG_W_SPEED) */

	if (Tsr == NULL) 
	{
	    nrad = front->npts_tan_sten/2;
	    sten = alloc_tan_stencil(front,nrad);
	    set_default_tan_stencil(sten);
	    alloc_state(intfc,&st_l3,max(sizeof(VGas),sizest));
	    alloc_state(intfc,&st_r1,max(sizeof(VGas),sizest));
	    alloc_state(intfc,&left,max(sizeof(VGas),sizest));
	    alloc_state(intfc,&right,max(sizeof(VGas),sizest));
	    alloc_state(intfc,&Tsl,sizest);
	    alloc_state(intfc,&Tsr,sizest);
	    alloc_state(intfc,&stemp,sizest);
	    alloc_state(intfc,&ansl_1,sizest);
	    alloc_state(intfc,&ansl_2,sizest);
	    alloc_state(intfc,&ansr_1,sizest);
	    alloc_state(intfc,&ansr_2,sizest);
	}

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES) 
	{
	    print_general_vector("nor = ",nor,dim,"");
	    (void) printf(", dn = %g, dt = %g\n",dn,dt);
	    (void) printf("Initial data -- sll, sl, sr, srr:\n");
	    verbose_print_state("sll",sll);
	    verbose_print_state("sl ",sl );
	    verbose_print_state("sr ",sr );
	    verbose_print_state("srr",srr);
	}
#endif /* defined(DEBUG_W_SPEED) */

	set_type_of_state(ansl,GAS_STATE);
	set_type_of_state(ansr,GAS_STATE);
	switch (w_type)
	{
	case PASSIVE_BOUNDARY: 
	case SUBDOMAIN_BOUNDARY:
	    /* PASSIVE and SUBDOMAIN states not set or used */
	    obstacle_state(intfc,ansl,sizest);
	    obstacle_state(intfc,ansr,sizest);
	    break;

	case DIRICHLET_BOUNDARY:
	    k = min(wssten->nsts-1,nrad);
	    for (i = 0; i < dim; ++i)
	        Coords(sten->p[0])[i] = pt[i];
	    
	    for (j = 1; j <= nrad; ++j)
	    {
	        for (i = 0; i < dim; ++i)
	        {
	            Coords(sten->p[j])[i]  = pt[i] + j*dn*nor[i];
	            Coords(sten->p[-j])[i] = pt[i] - j*dn*nor[i];
	        }
	    }
	    if (is_obstacle_state(sl)) /*left side is exterior*/
	    {
		comp = r_comp;

		for (j = 0; j <= nrad; ++j)
	            copy_state(sten->rightst[j],sr);
	        for (j = 1; j <= nrad; ++j)
	            hyp_solution(Coords(sten->p[-j]),l_comp,hs,UNKNOWN_SIDE,
			         front,wave,sten->rightst[-j],
				 sten->rightst[-j+1]);
		for (j = -nrad; j <= nrad; ++j)
	            Set_params(sten->rightst[j],sr);
	    }
	    else
	    {
		comp = l_comp;
		
		for (j = 0; j <= k; ++j)
		    copy_state(sten->rightst[-j],wssten->sl[j]);
		for (; j <= nrad; ++j)
		    hyp_solution(Coords(sten->p[-j]),comp,hs,UNKNOWN_SIDE,
				 front,wave,sten->rightst[-j],
				 sten->rightst[-j+1]);
		for (j = 1; j <= nrad; ++j)
		    copy_state(sten->rightst[j],sl);
		for (j = -nrad; j <= nrad; ++j)
	            Set_params(sten->rightst[j],sl);
	    }
	    sten->newhs = NULL;
	    sten->comp = comp;
	    sten->states = sten->rightst;
	    sten->dir = nor;
	    
	    one_side_npt_tang_solver(dn,dt,sten,ansl,front);

	    if (is_obstacle_state(sr)) /*right side is exterior*/
	    {
		comp = l_comp;
	        for (j = 0; j <= nrad; ++j)
	            copy_state(sten->rightst[-j],sl);
	        for (j = 1; j <= nrad; ++j)
	            hyp_solution(Coords(sten->p[j]),r_comp,hs,UNKNOWN_SIDE,
			         front,wave,sten->rightst[j],
				 sten->rightst[j-1]);
		for (j = -nrad; j <= nrad; ++j)
	            Set_params(sten->rightst[j],sl);
	    }
	    else
	    {
		comp = r_comp;
	        for (j = 0; j <= k; ++j)
	            copy_state(sten->rightst[j],wssten->sr[j]);
	        for (; j <= nrad; ++j)
	            hyp_solution(Coords(sten->p[j]),comp,hs,UNKNOWN_SIDE,
			         front,wave,sten->rightst[j],
				 sten->rightst[j-1]);
	        for (j = 1; j <= nrad; ++j)
	            copy_state(sten->rightst[-j],sr);
		for (j = -nrad; j <= nrad; ++j)
	            Set_params(sten->rightst[j],sr);
	    }
	    sten->newhs = NULL;
	    sten->comp = comp;
	    sten->states = sten->rightst;
	    sten->dir = nor;
	    
	    one_side_npt_tang_solver(dn,dt,sten,ansr,front);

#if defined(DEBUG_W_SPEED)
	    if (debug_w_speed == YES) 
	    {
	        (void) printf("Intermediate left, right states:\n");
	        verbose_print_state("ansl",ansl);
	        verbose_print_state("ansr",ansr);
	    }
#endif /* defined(DEBUG_W_SPEED) */

	        /*  Source terms are included above in the two  */
	        /* calls to one_side_npt_tang_solver(). */

	    if (!is_obstacle_state(sl))
	    {
	        evaluate_dirichlet_boundary_state(Coords(sten->p[0]),hs,
			front,wave,ansl);
		riemann_solution(0.0,nor,ansl,ansr,ansl,GAS_STATE);
#if !defined(COMBUSTION_CODE)
	        Set_params(ansl,sl);
#endif /* !defined(COMBUSTION_CODE) */
	    }
	    
	    if (!is_obstacle_state(sr))
	    {
		evaluate_dirichlet_boundary_state(Coords(sten->p[0]),hs,
			front,wave,ansl);
	        riemann_solution(0.0,nor,ansl,ansr,ansr,GAS_STATE);
#if !defined(COMBUSTION_CODE)
	        Set_params(ansr,sr);
#endif /* !defined(COMBUSTION_CODE) */
	    }

	    if (is_obstacle_state(sl))
	        obstacle_state(intfc,ansl,sizest);
	    if (is_obstacle_state(sr))
	        obstacle_state(intfc,ansr,sizest);

	    break;

	case NEUMANN_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
	    /* Use the Method of Characteristics in the normal direction */


	    if (is_obstacle_state(sr)) 
		side = NEGATIVE_SIDE;
	    else if (is_obstacle_state(sl)) 
		side = POSITIVE_SIDE;
	    else
	    {
	        side = UNKNOWN_SIDE;
		screen("ERROR in bdry_npt_w_speed(), invalid state "
		       "on NEUMANN_BOUNDARY,  no obstacle state found\n");
		clean_up(ERROR);
	    }
	    reflect_wssten(wssten,side,front);
	    wssten->w_type = CONTACT;
	    npt_w_speed(wssten,ansl_1,ansr_1,W);
	    alpha = wall_limiter(wssten,side,opts);
	    wssten->w_type = NEUMANN_BOUNDARY;

	    if (side == POSITIVE_SIDE) /* Right Side */
	    {
	        obstacle_state(intfc,ansl,sizest);
		for (i = 0; i < wssten->nsts; ++i)
		{
	            obstacle_state(intfc,wssten->sl[i],sizest);
	            obstacle_state(intfc,wssten->tsl[i],sizest);
	            obstacle_state(intfc,wssten->dsl[i],sizest);
		}

	        set_state(Tsr,TGAS_STATE,sr);
	        ur = 0.0;
	        if (include_wall_normal_velocity(front) == YES)
	            ur += scalar_product(nor,Vel(Tsr),dim);
	        cr = sound_speed(Tsr);
	        dn1 = cr * dt / dn;
	        ws_interpolate(st_r1,dn1,POSITIVE_SIDE,TGAS_STATE,wssten);
#if defined(DEBUG_W_SPEED)
	        if (debug_w_speed == YES)
	        {
	            (void) printf("Interpolated state, dn1 = %g\n",dn1);
	            verbose_print_state("st_r1",st_r1);
		    (void) printf("Calling neumann_moc()\n");
	        }
#endif /* defined(DEBUG_W_SPEED) */
	        (*opts->neumann_moc)(pt,Tsr,ur,cr,st_r1,side,
				     right,cr*dt,nor,front);
	        set_state(ansr_2,GAS_STATE,right);

	        if (Params(Tsr)->avisc.heat_cond != 0.0)
	        {
	            double heat_cond = Params(Tsr)->avisc.heat_cond;
		    double Tw = temperature(Tsr);
		    double Gam = gruneisen_gamma(Tsr);
		    double dS = heat_cond * (temperature(st_r1)/Tw - 1.0);
		    Dens(ansr_1) *= exp(-Gam*Tw*dS/sqr(cr));
		    reset_gamma(ansr_1);
		}
		interpolate_states(front,alpha,1.0-alpha,pt,
		                   ansr_1,pt,ansr_2,ansr);
		if (no_slip(hs))
		{
		    double alpha = 1.0 - adherence_coeff(hs);
		    alpha_state_velocity(alpha,ansr,dim);
		}
		zero_normal_velocity(ansr,nor,dim);
	    }
	    else if (side == NEGATIVE_SIDE) /* Left Side */
	    {
	        obstacle_state(intfc,ansr,sizest);
		for (i = 0; i < wssten->nsts; ++i)
		{
	            obstacle_state(intfc,wssten->sr[i],sizest);
	            obstacle_state(intfc,wssten->tsr[i],sizest);
	            obstacle_state(intfc,wssten->dsr[i],sizest);
		}

	        set_state(Tsl,TGAS_STATE,sl);
	        ul = 0.0;
	        if (include_wall_normal_velocity(front) == YES)
	            ul += scalar_product(nor,Vel(Tsl),dim);
	        cl = sound_speed(Tsl);
	        dn3 = 1.0 - cl * dt / dn;
	        ws_interpolate(st_l3,dn3,NEGATIVE_SIDE,TGAS_STATE,wssten);
#if defined(DEBUG_W_SPEED)
	        if (debug_w_speed == YES)
	        {
	            (void) printf("Interpolated state, dn3 = %g\n",dn3);
	            verbose_print_state("st_l3",st_l3);
		    (void) printf("Calling neumann_moc()\n");
	        }
#endif /* defined(DEBUG_W_SPEED) */
	        (*opts->neumann_moc)(pt,Tsl,ul,cl,st_l3,side,
				     left,-cl*dt,nor,front);
	        set_state(ansl_2,GAS_STATE,left);

	        if (Params(Tsl)->avisc.heat_cond != 0.0)
	        {
	            double heat_cond = Params(Tsl)->avisc.heat_cond;
		    double Tw = temperature(Tsl);
		    double Gam = gruneisen_gamma(Tsl);
		    double dS = heat_cond * (temperature(st_l3)/Tw - 1.0);
		    Dens(ansl_1) *= exp(-Gam*Tw*dS/sqr(cl));
    		    reset_gamma(ansr_1);
	        }
		interpolate_states(front,alpha,1.0-alpha,pt,
		                   ansl_1,pt,ansl_2,ansl);
		if (no_slip(hs))
		{
		    double alpha = 1.0 - adherence_coeff(hs);
		    alpha_state_velocity(alpha,ansl,dim);
		}
		zero_normal_velocity(ansl,nor,dim);
	    }
	    break;

	default:
	    screen("ERROR in bdry_npt_w_speed(), ");
	    print_wave_type("unknown wave_type = ",w_type,"\n",intfc);
	    clean_up(ERROR);
	}
	Left_w_speed("bdry_npt_w_speed",ansl,ansr,NULL,dim);
	return;
}		/*end bdry_npt_w_speed*/


LOCAL	void	apply_bc(
	Locstate   st,
	Wave       *wv,
	HYPER_SURF *hs)
{
	if (boundary_state_function(hs) ==
		constant_pressure_flow_through_boundary_state)
	{
	    Locstate bst = boundary_state(hs);
	    int      st_type = state_type(st);
	    set_state(st,TGAS_STATE,st);
	    Dens(st) = Dens(bst);
	    Press(st) = pressure(bst);
	    set_state(st,st_type,st);
	}
	/* For time dependent pressure Dirichlet boundary */
	else if(boundary_state_function(hs) ==
		time_dep_pressure_flow_through_boundary_state)
	{
	    Locstate    bst = boundary_state(hs);
	    int         st_type = state_type(st);
	    FD_DATA     *fd_data = (FD_DATA*)boundary_state_data(hs);
	    double	tr = fd_data->tr;
	    double	tp = fd_data->tp;
	    double	ts = fd_data->ts;
	    double	pr_a = fd_data->pr_a;
	    double	pr_p = fd_data->pr_p;
	    double       time = wv-> time;

	    set_state(bst,TGAS_STATE,bst);
	    set_state(st,TGAS_STATE,st);

	    if (time < tr) Press(st) = (pr_p-pr_a)/tr*time + pr_a;
	    else if (time <= tr+tp) Press(st) = pr_p;
	    else if (time <= tr+tp+ts) Press(st) = -(pr_p-pr_a)/ts
	    			*(time-tr-tp) + pr_p;
	    else Press(st) = pr_a;

	    Dens(bst) = Dens(st) = density(st);
	    Press(st) = Press(bst);
	    set_state(st,st_type,st);
	}
}		/*end apply_bc*/

/*
*			include_source():
*
*	This function incorporates the source terms introduced by gravity
*	and 3-D cylindrical symmetry.  Operator splitting (in time) is
*	used.
*
*	Input:	pt      - starting location of point being updated
*               ans 	- state to be updated
*		vnor0   - normal component of velocity at start of step
*		dt	- time step size
*		nor	- unit normal to front
*		W       - front velocity
*		g       - gravity uni_array
*		gr      - computational grid
*		w_type  - wave type for hypersurface being propagated
*
*	Output:	ans	- updated state
*			  obstacle states are not updated.	
*
*	The operator split algorithm updates the equation
*
*	rho_t + alpha*rho*nor[0]*(u,nor)/r(t) = 0
*	u_t = (g,nor)*nor
*	S_t = 0
*
*	where nor is the normal vector to the surface, S is the entropy,
*	u is the vector velocity, rho is the mass density, and g is the body
*	source. The geometry factor, alpha,  is zero for rectangular geometry,
*	one for cylindrical, and two for spherical geometry.  The position
*	r moves with velocity W[0],  r(t) = r0 + W[0]*t.
*/		

EXPORT	void include_source(
	double		*pt,
	Locstate	ans,
	double           vnor0,
	double		dt,
	double		*nor,
	double		*W,
	double		*g,
	RECT_GRID	*gr,
	double		diffu,
	int             w_type)
{
	static boolean	first = YES;
	static boolean	is_grav;
	static boolean	add_source_terms;
	double		g_factor;
	int		dim = gr->dim;
	int		stype;
	int		i;
	boolean		diffu_term = NO;


	static GEOMETRY	geom;
	static double	alpha;
	GEOM_SOURCE_METHOD geom_source_method = G3ptOpts.geom_source_method;

	if (first == YES) 
	{
	    first = NO;
	    is_grav = is_gravity();
	    add_source_terms = is_grav;

	    if (is_rotational_symmetry())
            {
	        geom = Geometry(&alpha);
	        if (geom != RECTANGULAR)
	            add_source_terms = YES;
            }
	}

	if (!add_source_terms )
	    return;

	if (is_obstacle_state(ans))
	    return;

	stype = state_type(ans);
	set_state(ans,TGAS_STATE,ans);
	if ((is_grav) && (g != NULL))
	{
	    if (!diffu_term)
            {
	        /*Dens(ans) += diffu*Dens(ans)*dt; */
	        g_factor = scalar_product(nor,g,dim)*dt;
	        for (i = 0; i < dim; ++i)
	            Vel(ans)[i] += nor[i]*g_factor;
	    }
	    else
            {
	        /*Dens(ans) += diffu*Dens(ans)*dt; */
		g_factor = scalar_product(nor,g,dim)*dt;
		diffu = max(diffu,0.0);
	        for (i = 0; i < dim; ++i)
	            Vel(ans)[i] = exp(-diffu*dt)*Vel(ans)[i] + 
			               nor[i]*g_factor;
	    }
	}
	else
	{
	    diffu = max(diffu,0.0);
	    Dens(ans) += diffu*Dens(ans)*dt;
	}
	if (is_rotational_symmetry() && geom != RECTANGULAR)
	{
	    double    rmin, r0, r1;
	    double    v[3];

	    rmin = pos_radius(0.0,gr);
	    r0 = pos_radius(pt[0],gr);
	    r1 = pos_radius(pt[0] + W[0]*dt,gr);

	    if ((fabs(r0) > rmin) && (fabs(r1) > rmin))
	    {
	        double rho, c2, pm, p;
		double rho_fac;
		double pa, pb, p0;
	        double vn;
		double rm = (geom == SPHERICAL) ?
	            2.0*(r0*r0+r0*r1+r1*r1)/(3.0*(r0+r1)) : 0.5*(r0+r1);
                int number_of_interations;
		const int MAX_NUM_ITERATIONS = 20;

		switch (geom_source_method)
		{
		case ANALYTIC_SOLUTION:
		    if (geom == SPHERICAL)
			rho_fac = (r0*r0)/(r1*r1);
	            else if (geom == CYLINDRICAL)
			rho_fac = r0/r1;
		    else
		    {
			rho_fac = 1.0;
			screen("ERROR in include_source(), unknown geometry\n");
			clean_up(ERROR);
		    }
		    if (!is_scalar_wave(w_type))
		    {
	                vn = nor[0]*0.5*(vnor0+scalar_product(Vel(ans),nor,dim));
		        rho_fac *= exp(alpha*(W[0]-vn)*dt/rm);
		    }
		    rho = Dens(ans)*rho_fac;
		    Press(ans) = pressure_rarefaction(rho,ans);
		    Dens(ans) = rho;
		    reset_gamma(ans);
		    break;

		case  BACKWARD_EULER:
	            vn = nor[0]*0.5*(vnor0+scalar_product(Vel(ans),nor,dim));
                    pa = p0 = Press(ans);
                    number_of_interations = 0;
	            for (i = 0; i < dim; ++i)
		        v[i] = Vel(ans)[i];
                    do
                    {
                        pb = pa;
                        c2 = sound_speed_squared(ans);
                        pa = p0 - alpha*Dens(ans)*c2*vn*dt/rm;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
                        pa = max(Min_pressure(ans),pa);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
                        state_on_adiabat_with_pr(ans,pa,ans,state_type(ans));
                    } while ((fabs(pb - pa) > EPSILON*pa) &&
                         (++number_of_interations < MAX_NUM_ITERATIONS));
	            for (i = 0; i < dim; ++i)
		        Vel(ans)[i] = v[i];
		    break;

		case MODIFIED_EULER:
		default:
	            vn = nor[0]*0.5*(vnor0+scalar_product(Vel(ans),nor,dim));
	            p = Press(ans);
	            c2 = sound_speed_squared(ans);
	            rho = Dens(ans);
	            pm = p - 0.5*alpha*rho*c2*vn*dt/r0;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	            pm = max(Min_pressure(ans),pm);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	            Dens(ans) = dens_rarefaction(pm,ans);
	            Press(ans) = pm;
		    reset_gamma(ans);
	            c2 = sound_speed_squared(ans);

	            p -= alpha*Dens(ans)*c2*vn*dt/rm;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	            p = max(Min_pressure(ans),p);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	            Dens(ans) = dens_rarefaction(p,ans);
	            Press(ans) = p;
		    reset_gamma(ans);
		    break;
		}
                /* Mass fraction is preserved, scaling is performed */
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    int    k;
                    double  sum;
                    if(Params(ans)->n_comps != 1)
                    {
                        sum = 0.0;
                        for(k = 0; k < Params(ans)->n_comps; k++)
                            sum += pdens(ans)[k];
                        for(k = 0; k < Params(ans)->n_comps; k++)
                            pdens(ans)[k] *= Dens(ans)/sum;
                    }
                }
	    }
	}

	set_state(ans,stype,ans);
}	/*end include_source*/

LOCAL	boolean	strong_wave(
	Locstate	s0,
	double		pm,
	double		um)
{
	double	Mach_tol = G3ptOpts.Mach_tol;
	double	i0 = acoustic_impedance(s0);
	double   p0 = pressure(s0);
	double	idu = i0*fabs(um - vel(0,s0));
	double	dp = fabs(pm - p0);

	if (Mach_tol == HUGE_VAL)
	    return NO;
	
	if ((dp + idu) < EPSILON*p0) /*TOLERANCE*/
	    return NO;
	if (fabs(dp - idu) <= Mach_tol*max(idu,dp))
	    return NO;
	return YES;
}		/*end strong_wave*/

LOCAL	boolean	strong_interaction(
	Locstate	sl,
	Locstate	sr,
	int		w_type,
	double		pml,
	double		pmr,
	double		uml,
	double		umr)
{
	double	rhol, rhor;
	double	A_tol = G3ptOpts.A_tol;

	if (debugging("no_si"))
	    return NO;

	if (is_forward_wave(w_type))
	{
	    if (strong_wave(sl,pml,uml))
	        return YES;
	    rhol = dens_Hugoniot(pml,sl);
	    rhor = dens_Hugoniot(pmr,sr);
	    if (fabs(rhol - rhor)/(rhol+rhor) > A_tol)
	    {
#if defined(DEBUG_W_SPEED)
	        if (debug_w_speed == YES)
	            (void) printf("strong density jump for forward wave "
				  "rhol = %g, rhor = %g, "
				  "A = %g, A_tol = %g\n",rhol,rhor,
				  (rhol - rhor)/(rhol+rhor),A_tol);
#endif /* defined(DEBUG_W_SPEED) */
	        return YES;
	    }
	    if (is_rarefaction_wave(w_type))
	    {
	        if (strong_wave(sr,pmr,umr))
	            return YES;
	    }
	}
	else if (is_scalar_wave(w_type))
	{
	    if ((strong_wave(sl,pml,uml)) || (strong_wave(sr,pmr,umr)))
	        return YES;
	}
	else if (is_backward_wave(w_type))
	{
	    if (strong_wave(sr,pmr,umr))
	        return YES;
	    rhol = dens_Hugoniot(pml,sl);
	    rhor = dens_Hugoniot(pmr,sr);
	    if (fabs(rhol - rhor)/(rhol+rhor) > A_tol)
	    {
#if defined(DEBUG_W_SPEED)
	        if (debug_w_speed == YES)
	            (void) printf("strong density jump for backward wave "
				  "rhol = %g, rhor = %g, "
				  "A = %g, A_tol = %g\n",rhol,rhor,
				  (rhol - rhor)/(rhol+rhor),A_tol);
#endif /* defined(DEBUG_W_SPEED) */
	        return YES;
	    }
	    if (is_rarefaction_wave(w_type))
	    {
	        if (strong_wave(sl,pml,uml))
	            return YES;
	    }
	}
	return NO;
}		/*end strong_interaction*/

/*
*			cheap method of characteristics
*
*	The following functions solve the Euler equations in characteristic
*	form.  The particular equations solved are
*
*		d P           d N
*              ----   - rho*c ---  = 0		where l- = N - c
*               d l-          d l-
*
*		d T
*              ----   = 0			where l0 = N
*	        d l0
*
*		d P          d rho
*              ----  -  c^2 ------   = 0	where l0 = N
*		d l0         d l0
*
*		d P           d N
*              ----   + rho*c ---  = 0		where l+ = N + c
*               d l+          d l+
*
*	Here P is the pressure, rho the density, and N and T are the
*	normal and tangential components of velocity respectively.
*/


/*
*			neumann_cheap_moc():
*
*	Performs a simple method of characteristic update of the
*	state near a neumann boundary.  The method is to use
*	a finite difference approximation to the characteristic
*	equations for the pressure and velocity.
*
*	Assumes that all of the states are in TGas format.
*/

/*ARGSUSED*/
LOCAL	void neumann_cheap_moc(
	double    *pt,		/* Wall point position           */
	Locstate sw,		/* State at the wall             */
	double    uw,		/* Wall velocity                 */
	double    cw,		/* Sound speed at the wall       */
	Locstate sfoot,		/* Incoming characteristic state */
	SIDE     side,		/* Interior side of the wall     */
	Locstate ans,		/* Output answer state           */
	double    dn,		/* Normal direction mesh spacing */
	double    *nor,		/* Wall normal                   */
	Front    *front)	/* Front structure               */
{
	RECT_GRID *gr = front->rect_grid;
	double     ufoot, p, rw, pw, du, vnor;
	double     dp, heat_cond;
	double     dS, Tw, Gam;
	double     time = front->time;
	double     dt = front->dt;
	double     g1[3], g;
	double     pt1[3];
	int       i, dim;

	static	double	alpha;
	static	boolean	first = YES;

	if (first == YES)
	{
	    first = NO;
	    alpha = rotational_symmetry();
	}

	debug_print("w_speed","Entered neumann_cheap_moc()\n");

	dim = gr->dim;
	for (i = 0; i < dim; ++i)
	    pt1[i] = pt[i]+dn*nor[i];
	eval_gravity(pt1,time,g1);
	g = scalar_product(g1,nor,dim);
	ufoot = scalar_product(nor,Vel(sfoot),dim);
	rw = Dens(sw);
	pw = Press(sw);
	du = uw - ufoot;

	    /* source term correction to MOC: */

	if (is_gravity() == YES)
	    du -= g*dt;
	if (side == NEGATIVE_SIDE)
	    du = -du;

	if (is_rotational_symmetry() && alpha > 0.0)
	{
	    double rmin = pos_radius(0.0,gr);
	    double rad = pos_radius(pt1[0],gr);
	    if (fabs(rad) > fabs(rmin))
	        du -= alpha*nor[0]*cw*ufoot*dt/rad;
	}

	Press(ans) = p = Press(sfoot) + rw * cw * du;
	heat_cond = Params(sw)->avisc.heat_cond;
	if (heat_cond != 0.0)
	{
	    Tw = temperature(sw);
	    Gam = gruneisen_gamma(sw);
	    dS = heat_cond * (temperature(sfoot)/Tw - 1.0);
	    dp = p - pw - rw*Tw*Gam*dS;
	}
	else
	    dp = p - pw;
	Dens(ans) = rw + dp / sqr(cw);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(i = 0; i < Params(sw)->n_comps; i++)
                pdens(ans)[i] = (pdens(sw)[i]/Dens(sw))*Dens(ans);
        }
	vnor = scalar_product(nor,Vel(sw),dim);
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = Vel(sw)[i] - nor[i] * vnor;
	Set_params(ans,sw);
	set_type_of_state(ans,TGAS_STATE);
	reset_gamma(ans);
	debug_print("w_speed","Left neumann_cheap_moc()\n");
}	/*end neumann_cheap_moc*/

LOCAL	void	wspeed_neumann_riem_inv_moc(
	double    *pt,
	Locstate sw,
	double	 uw,
	double	 cw,
	Locstate sfoot,
	SIDE     side,
	Locstate ans,
	double    dn,
	double	 *nor,
	Front    *fr)
{
	debug_print("w_speed","Entered wspeed_neumann_riem_inv_moc()\n");
	neumann_riem_inv_moc(pt,sw,uw,cw,sfoot,side,ans,dn,nor,fr);
	debug_print("w_speed","Left wspeed_neumann_riem_inv_moc()\n");
}	/*end wspeed_neumann_riem_inv_moc*/




/*
*		shock_ahead_state_cheap_moc():
*
*	Performs a simple method of characteristics update of the state
*	ahead of a shock wave, by using a simple finite difference version
*	of the characteristic equations for the pressure and velocity.
*
*		Assumes that all state are in TGas format.
*	
*/

LOCAL	void shock_ahead_state_cheap_moc(
	double     *pt,
	Locstate  st0,
	Locstate  st1,
	Locstate  st2,
	Locstate  st3,
	Locstate  ans,
	double     dn,
	double     f1,
	double     f2,
	double     f3,
	double     *nor,
	double     *W,
	int       add_source,
	double     dt,
	Front     *front)
{
	RECT_GRID *gr = front->rect_grid;
	double	  time = front->time;
	double	  u1, u2, u3;
	double	  u, vtan[MAXD];
	double	  r2, c2, p1, p2, p3;
	double	  p;
	double     g, pt2[3];
	int	  i, dim;
	double	  alpha = rotational_symmetry();

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	    (void) printf("Entered shock_ahead_state_cheap_moc()\n");
#endif /* defined(DEBUG_W_SPEED) */

	dim = Params(st1)->dim;
	u1 = u2 = u3 = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    u1 += nor[i]*Vel(st1)[i];
	    u2 += nor[i]*Vel(st2)[i];
	    u3 += nor[i]*Vel(st3)[i];
	    pt2[i] = pt[i] + f2*dn*nor[i];
	}
	if (add_source)
	{
	    g = scalar_product(gravity(pt2,time),nor,dim);
	}
	else
	    g = 0.0;
#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	    (void) printf("u1 = %g, u2 = %g, u3 = %g\n",u1,u2,u3);
#endif /* defined(DEBUG_W_SPEED) */
	p1 = Press(st1);
	r2 = Dens(st2);
	p2 = Press(st2);
	c2 = sound_speed(st2);
	p3 = Press(st3);
	Set_params(ans,st1);
	set_type_of_state(ans,TGAS_STATE);
	p = 0.5 * ( p3 + p1  + r2*c2*(u3 - u1) );

	if (is_rotational_symmetry() && add_source && alpha > 0.0)
	{
	    double rmin = pos_radius(0.0,gr);
	    double rad = pos_radius(pt[0],gr);
	    if (fabs(rad) > fabs(rmin))
	        p -= alpha*nor[0]*r2*sqr(c2)*u2*dt/rad;
	}

	Press(ans) = p;
	Dens(ans) = r2 + (p - p2) / sqr(c2);
	u = 0.5 * ( (u3 + u1) + (p3 - p1) / (r2 * c2));
	if (add_source && (is_gravity() == YES))
	    u += g*dt;
	for (i = 0; i < dim; ++i)
	{
	    vtan[i] = Vel(st2)[i] - nor[i]*u2;
	    /*Vel(ans)[i] = vtan[i]*r2/Dens(ans) + nor[i]*u;*/
	    Vel(ans)[i] = vtan[i] + nor[i]*u;
	}
	reset_gamma(ans);
	
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(ans) < Vacuum_dens(ans))
	{
	    shock_ahead_state_riem_inv_moc(pt,st0,st1,st2,st3,ans,dn,f1,f2,f3,
					   nor,W,add_source,dt,front);
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	{
	    verbose_print_state("Answer from shock_ahead_state_cheap_moc()",
	                        ans);
	    (void) printf("Left shock_ahead_state_cheap_moc()\n");
	}
#endif /* defined(DEBUG_W_SPEED) */

}	/*end shock_ahead_state_cheap_moc*/

LOCAL	void wspeed_shock_ahead_state_riem_inv_moc(
	double     *pt,
	Locstate  st0,
	Locstate  st1,
	Locstate  st2,
	Locstate  st3,
	Locstate  ans,
	double     dn,
	double     f1,
	double     f2,
	double     f3,
	double     *nor,
	double     *W,
	int       add_source,
	double     dt,
	Front     *front)
{
	shock_ahead_state_riem_inv_moc(pt,st0,st1,st2,st3,ans,dn,f1,f2,f3,
				       nor,W,add_source,dt,front);
}	/*end wspeed_shock_ahead_state_riem_inv_moc*/


LOCAL	void	contact_riemann(
	Locstate	left,
	Locstate	Tsl,
	Locstate	right,
	Locstate	Tsr,
	Locstate	st_l3,
	double           dnl,
	Locstate	st_r1,
	double           dnr,
	Locstate	ansl,
	Locstate	ansr,
	WSSten          *sten,
	double		*W,
	NptWSpeedOpts	*opts,
	double		*grav,
	double		difful,
	double		diffur)
{
	double	*pt = sten->coords;
	double	V[3];
	double   pjump = sten->pjump;
	double	*nor = sten->nor;
	double	dt = sten->dt;
	double   vnor0;
	int	w_type = sten->w_type;
	Front	*front = sten->front;
	int	i, dim = front->rect_grid->dim;

	contact_filter_outgoing_wave(right,Tsr,st_r1,dnr,
		                     sten,POSITIVE_SIDE,opts);
#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES) 
	{
	    verbose_print_state("st_r1 after separation of back "
	                        "mode at foot of characteristic",st_r1);
	}
#endif /* defined(DEBUG_W_SPEED) */

	contact_filter_outgoing_wave(left,Tsl,st_l3,dnl,sten,NEGATIVE_SIDE,opts);

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES) 
	{
	    verbose_print_state("st_l3 after separation of back "
	                        "mode at foot of characteristic",st_l3);
	}
#endif /* defined(DEBUG_W_SPEED) */

	/* Solve Riemann problem to find left and right states */

	set_state(st_r1,GAS_STATE,st_r1);
	set_state(st_l3,GAS_STATE,st_l3);
	w_speed(pt,st_l3,st_r1,ansl,ansr,V,pjump,nor,w_type,front);

/*#if defined(DEBUG_W_SPEED) */
/*	if (debug_w_speed == YES)  */
	if(debugging("tst_nan"))
	{
	    (void) printf("States after w_speed\n");
	    verbose_print_state("st_l3",st_l3);
	    verbose_print_state("st_r1",st_r1);
	    verbose_print_state("ansl ",ansl);
	    verbose_print_state("ansr ",ansr);
	}
/*#endif*/  /*defined(DEBUG_W_SPEED) */ 

	if (is_gravity() == YES)
	{
	    double g = scalar_product(grav,nor,dim);

	    /*
	    *  source terms:  include_source() takes care of the left
	    *  and right states, but we also need to accelerate the
	    *  contact if there's an external force (which for now is
	    *  only gravity).
	    */

	    for (i = 0; i < dim; ++i)
	        V[i] += nor[i] * g * dt;
	}

	/* Use centered in time difference for wave speed */

	vnor0 = scalar_product(W,nor,dim);
	for (i = 0; i < dim; ++i)
	    W[i] = 0.5*(W[i] + V[i]);

	include_source(pt,ansl,vnor0,dt,nor,W,grav,front->rect_grid,
			difful,w_type);
	include_source(pt,ansr,vnor0,dt,nor,W,grav,front->rect_grid,
			diffur,w_type);
}		/*end contact_riemann*/

#if DONT_COMPILE
LOCAL	void	filter_wave(Locstate,Locstate,Locstate,double*,SIDE);
LOCAL	void	filter_wave(
	Locstate	sl,
	Locstate	sr,
	Locstate	sfilter,
	double		*nor,
	SIDE		side)
{
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave, wave;
	double	 *vahead, vl[3], vr[3];
	double	 pl, pr, ul, ur, ml, mr, p, u, m;
	int	 state_type = state_type(sfilter);
	int	 i, dim = Params(sl)->dim;
	int      family;
	Locstate ahead;
	static	Locstate	left = NULL, right = NULL;

	if (left == NULL)
	{
	    (*Params(sl)->_alloc_state)(&left,max(sizeof(VGas),
	                                Params(sl)->sizest));
	    (*Params(sr)->_alloc_state)(&right,max(sizeof(VGas),
	                                Params(sr)->sizest));
	}
	set_state_for_find_mid_state(left,sl);
	set_state_for_find_mid_state(right,sr);
	Vel(left)[0] = scalar_product(VelocityVector(left,vl),nor,dim);
	Vel(right)[0] = scalar_product(VelocityVector(right,vr),nor,dim);

	if (find_mid_state(left,right,0.0,&pl,&pr,&ul,&ur,&ml,&mr,
	                   &l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	{
	    (void) printf("WARNING in filter_wave(), "
	                  "find_mid_state() did not converge for ");
	                  "contact wave RP\n");
	}
	if (side == POSITIVE_SIDE)
	{
	    ahead = right;
	    p = pr;
	    u = ur;
	    m = mr;
	    wave = r_wave;
	    vahead = vr;
	    family = RIGHT_FAMILY;
	}
	else if (side == NEGATIVE_SIDE)
	{
	    ahead = left;
	    p = pl;
	    u = ul;
	    m = ml;
	    wave = l_wave;
	    vahead = vl;
	    family = LEFT_FAMILY;
	}
	else
	{
	    screen("ERROR in filter_wave, unknown side %d\n",side);
	    clean_up(ERROR);
	}
	midstate(ahead,sfilter,m,u,p,TGAS_STATE,wave,family);
	for (i = 0; i < dim; ++i)
	    Vel(sfilter)[i] = vahead[i] + (u - Vel(ahead)[0])*nor[i];
	set_state(sfilter,state_type,sfilter);
}		/*end filter_wave*/
#endif /* DONT_COMPILE */


/*ARGSUSED*/
LOCAL	void	contact_filter_outgoing_wave(
	Locstate	scontact,
	Locstate	Tscontact,
	Locstate	sfilter,
	double           dsfilter,
	WSSten		*sten,
	SIDE		side,
	NptWSpeedOpts	*opts)
{
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave, wave;
	Locstate	sl, sr;
	double	*nor = sten->nor;
	double	pl, pr, ul, ur, ml, mr, p, u, m;
	double	r;
	double	rc, uc = Vel(scontact)[0];
	double	sgn;
	double   vnor_sfilter, vtan[3];
	Front   *front = sten->front;
	int	dim = front->rect_grid->dim;
	int	i;
	static	Locstate tmpst = NULL;
#if defined(COMBUSTION_CODE)
	int	family;
#endif /* defined(COMBUSTION_CODE) */
        double   r_sfilter;

	if (tmpst == NULL)
	    alloc_state(front->interf,&tmpst,front->sizest);

	set_state(sfilter,TGAS_STATE,sfilter);
	vnor_sfilter = scalar_product(Vel(sfilter),nor,dim);
	for (i = 0; i < dim; ++i)
	    vtan[i] = Vel(Tscontact)[i] - nor[i]*uc;
	if (opts->_scalar_filter_outgoing_waves_at_contact == NO)
	{
	    for (i = 0; i < dim; ++i)
	        Vel(sfilter)[i] = vtan[i] + vnor_sfilter*nor[i];
	    return;
	}
        r_sfilter = Dens(sfilter);

	rc = Dens(scontact);

	/* Load normal component of velociy of sfilter in FRONT local coords */
	/* The other components are not needed here. */
	Vel(sfilter)[0] = vnor_sfilter;

	if (side == POSITIVE_SIDE)
	{
	    sl = scontact;
	    sr = sfilter;
	    sgn = -1.0;
#if defined(COMBUSTION_CODE)
	    family = LEFT_FAMILY;
#endif /* defined(COMBUSTION_CODE) */
	}
	else
	{
	    sl = sfilter;
	    sr = scontact;
	    sgn = 1.0;
#if defined(COMBUSTION_CODE)
	    family = RIGHT_FAMILY;
#endif /* defined(COMBUSTION_CODE) */
	}

	/* Separate out back mode at foot of characteristic */
	/* by solving Riemann problem between scontact, sfilter */

	set_state_for_find_mid_state(sfilter,sfilter);

	if (!find_mid_state(sl,sr,0.0,&pl,&pr,&ul,&ur,&ml,&mr,&l_wave,&r_wave))
	{
	    (void) printf("WARNING in contact_filter_outgoing_wave(), "
	                  "find_mid_state() did not converge for "
	                  "contact wave RP\n");
	}

	if (side == POSITIVE_SIDE)
	{
	    p = pl;
	    u = ul;
	    m = ml;
	    wave = l_wave;
	}
	else
	{
	    p = pr;
	    u = ur;
	    m = -mr;
	    wave = r_wave;
	}

	switch (wave) 
	{
	case SHOCK:
	    r = m/(u - uc + m/rc);
	    if (invalid_shock(pressure(Tscontact),Dens(Tscontact),p,r))
		r = dens_Hugoniot(p,Tscontact);
	    Set_params(sfilter,Tscontact);
	    break;

	case RAREFACTION:
	    r = dens_rarefaction(p,Tscontact);
	    Set_params(sfilter,Tscontact);
	    if (pressure_is_near_vacuum(p,Tscontact))
	        u = uc + sgn*riemann_wave_curve(scontact,p);
	    break;

#if defined(COMBUSTION_CODE)
	case STRONG_DET:
	    Set_other_params(sfilter,Tscontact);
	    r = m/(u - uc + m/rc);
	    break;

	case CJ_DET:
	    CJ_det(tmpst,TGAS_STATE,scontact,family);
	    Set_other_params(sfilter,Tscontact);
	    r = Dens(tmpst);
	    u = Vel(tmpst)[0];
	    break;
#endif /* defined(COMBUSTION_CODE) */

	default:
	    screen("ERROR in contact_filter_outgoing_wave(), "
	           "unknown wave %d\n",wave);
	    clean_up(ERROR);
	    break;
	}

	Dens(sfilter) = r;	Press(sfilter) = p;
        /* scaling */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(sfilter)->n_comps != 1)
            {
                for(i = 0; i < Params(sfilter)->n_comps; i++)
                    pdens(sfilter)[i] = (pdens(sfilter)[i]/r_sfilter)*Dens(sfilter);
            }
        }
	set_type_of_state(sfilter,TGAS_STATE);
	for (i = 0; i < dim; ++i)
	    Vel(sfilter)[i] = vtan[i] + u*nor[i];

#if defined(COMBUSTION_CODE)
	if (Composition_type(scontact) == ZND)
	    React(sfilter) = React(Tscontact);
#endif /* defined(COMBUSTION_CODE) */
	reset_gamma(sfilter);
}		/*end contact_filter_outgoing_wave*/

/*
*		contact_cheap_moc():
*
*	Updates the states on a contact discontinuity using
*	a simple finite difference approximation to the
*	characteristic equations for the pressure and normal 
*	velocity.
*
*	Assumes the input states are in TGas format.
*/

LOCAL	void contact_cheap_moc(
	double	  *pt,
	Locstate  sl,
	Locstate  sr,
	Locstate  sl3,
	Locstate  sr1,
	Locstate  ansl,
	Locstate  ansr,
	double	  cl,
	double	  cr,
	double	  pjump,
	double	  dnl,
	double     dnr,
	double	  *nor,
	double	  *W,
	Front     *front)
{
	RECT_GRID *gr = front->rect_grid;
	double     u1, u3, p1, p3;		/* 1, 3 state nor vel, press */
	double     pansl, pansr;		/* answer press */
	double     uansl, uansr, u;	/* answer normal velocity */
	double     ul, ur;
	double     rr, pr, rl, pl;		/* left, right density, press */
	double     den, rcr, rcl;
	double     pt3[3], pt1[3];
	double     dgdt, gv3[3], gv1[3], g1dt, g3dt;
	double     dt = front->dt;
	double     time = front->time;
	double     pfac = 0.0, ulfac = 0.0, urfac = 0.0, ufac = 0.0;
	int       i, dim = gr->dim;

	double	alpha = rotational_symmetry();

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	    (void) printf("Entered contact_cheap_moc()\n");
#endif /* defined(DEBUG_W_SPEED) */

	rr = Dens(sr);		pr = Press(sr);
	rl = Dens(sl);		pl = Press(sl);
	p1 = Press(sr1);	p3 = Press(sl3);
	u1 = u3 = ul = ur = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    u1 += nor[i] * Vel(sr1)[i];
	    u3 += nor[i] * Vel(sl3)[i];
	    ul += nor[i] * Vel(sl)[i];
	    ur += nor[i] * Vel(sr)[i];
	    pt3[i] = pt[i] + dnl*nor[i];
	    pt1[i] = pt[i] + dnr*nor[i];
	}
	eval_gravity(pt3,time,gv3);
	eval_gravity(pt1,time,gv1);
	g3dt = scalar_product(gv3,nor,dim)*dt;
	g1dt = scalar_product(gv1,nor,dim)*dt;
	dgdt = g1dt - g3dt;
	/*ur = ul = 0.5 * (ul + ur);*/

	rcr = rr*cr;	rcl = rl*cl;	den = rcr + rcl;
	if (is_rotational_symmetry() && alpha > 0.0) 
	{
	    double rmin = pos_radius(0.0,gr);
	    double radr = pos_radius(pt1[0],gr);
	    double radl = pos_radius(pt3[0],gr);
	    if ((fabs(radr) > fabs(rmin)) && (fabs(radl) > fabs(rmin)))
	    {
	        pfac = -nor[0]*dt*rcl*rcr*alpha*(u1*cr/radr + u3*cl/radl)/den;
		ulfac = (nor[0]*dt*alpha*cl*rcl*u3/radl)/den;
		urfac = (nor[0]*dt*alpha*cr*rcr*u1/radr)/den;
		ufac = urfac - ulfac;
	    }
	}

	pansr = (rcr*p3+rcl*p1 + rcl*rcr*(u3-u1-dgdt) - rcr*pjump)/den + pfac;
	pansl = (rcr*p3+rcl*p1 + rcl*rcr*(u3-u1-dgdt) + rcl*pjump)/den + pfac;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pansr = max(Min_pressure(sr),pansr);
	pansl = max(Min_pressure(sl),pansl);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	u = (rcl*(u3+g3dt) + rcr*(u1+g1dt) + p3 - p1 - pjump)/den + ufac;

	uansl = (pressure_is_near_vacuum(pansl,sl)) ? u3+g3dt+p3/rcl-ulfac : u;
	uansr = (pressure_is_near_vacuum(pansr,sr)) ? u1+g1dt-p1/rcr+urfac : u;

	for (i = 0; i < dim; ++i)
	{
	    Vel(ansl)[i] = Vel(sl)[i]  + (uansl - ul) * nor[i];
	    Vel(ansr)[i] = Vel(sr)[i]  + (uansr - ur) * nor[i];
	}
	Dens(ansr) = rr + (pansr - pr) / (cr * cr);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* mass fraction of sr1 is preserved */
            if(Params(sr1)->n_comps != 1)
            {
                for(i = 0; i < Params(sr1)->n_comps; i++)
                    pdens(ansr)[i] = (pdens(sr1)[i]/Dens(sr1))*Dens(ansr);
            }
        }
	Press(ansr) = pansr;
	Set_params(ansr,sr);
	set_type_of_state(ansr,TGAS_STATE);
        reset_gamma(ansr);

	set_state(ansr,GAS_STATE,ansr);

	Dens(ansl) = rl + (pansl - pl) / (cl * cl);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* mass fraction of sl3 is preserved */
            if(Params(sr1)->n_comps != 1)
            {
                for(i = 0; i < Params(sl3)->n_comps; i++)
                    pdens(ansl)[i] = (pdens(sl3)[i]/Dens(sl3))*Dens(ansl);
            }
        }
	Press(ansl) = pansl;
	Set_params(ansl,sl);
	set_type_of_state(ansl,TGAS_STATE);
        reset_gamma(ansl);

	set_state(ansl,GAS_STATE,ansl);

	    /* Use centered in time difference for wave speed */

	for (i = 0; i < dim; ++i)
	    W[i] = 0.5*(W[i] + nor[i]*u);

#if defined(DEBUG_W_SPEED)
	if (debug_w_speed == YES)
	{
	    verbose_print_state("Answer right from contact_cheap_moc()",ansr);
	    verbose_print_state("Answer left from contact_cheap_moc()",ansl);
	    (void) printf("Left contact_cheap_moc()\n");
	}
#endif /* defined(DEBUG_W_SPEED) */
}	/*end contact_cheap_moc*/

/*
*			contact_moc_plus_rh():
*
*	The states on opposite sides of a contact discontinuity
*	are updated using the method of characteristics and
*	the Rankine-Hugoniot relations across the contact.
*
*	All input states are assumed to be in TGas format.
*
*	The answer states ansl and ansr are returned in Gas format.
*
*/

typedef struct {
	Locstate st3l, st1r;
	Locstate st2l, st2r;
	Locstate ansl, ansr;
	double rc3l, rc1r;
	double u3l, u1r;
	double p3l, p1r;
	double pjump;
	double p_min;
	double alpha;
	double rfac;
} C_MPRH;


LOCAL void contact_moc_plus_rh(
	double	 *pt,
	double	 pm,
	Locstate st2l,
	Locstate st2r,
	Locstate st3l,
	Locstate st1r,
	Locstate ansl,
	Locstate ansr,
	double	 pjump,
	double	 dnl,
	double    dnr,
	double	 *nor,
	double	 *W,
	Front	*front)
{
	RECT_GRID *gr = front->rect_grid;
	const double meps = MACH_EPS;/*TOLERANCE*/
	double     time = front->time;
	double     dt = front->dt;
	double     pt1[3], pt3[3], pta[3];
	double     g1, g3, gv1[3], gv3[3], gva[3];
	double	  u1r, u2r, u2l, u3l;
	double	  ul, ur;
	double	  pi, u;
	double	  x0, a, b;
	double	  delta, epsilon;
	double	  p_min = max(Min_pressure(st3l),Min_pressure(st1r));
	int	  i, dim = gr->dim;
	static const int mnth = 10; /*TOLERANCE*/
	C_MPRH	  Cmprh;
	double     alpha = rotational_symmetry();
	Cmprh.alpha = alpha;

#if defined(DEBUG_W_SPEED)
	debug_print("cmph","Entered contact_moc_plus_rh()\n");
#endif /* defined(DEBUG_W_SPEED) */

	Cmprh.pjump = pjump;
	Cmprh.p_min = p_min;

	Cmprh.ansr = ansr;
	Cmprh.st2r = st2r;
	Cmprh.st1r = st1r;
	Cmprh.rc1r = acoustic_impedance(st1r);
	Cmprh.p1r = pressure(st1r);

	Cmprh.ansl = ansl;
	Cmprh.st2l = st2l;
	Cmprh.st3l = st3l;
	Cmprh.rc3l = acoustic_impedance(st3l);
	Cmprh.p3l = pressure(st3l);

	u1r = u2r = u2l = u3l = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    u1r += nor[i] * Vel(st1r)[i];
	    u2r += nor[i] * Vel(st2r)[i];
	    u2l += nor[i] * Vel(st2l)[i];
	    u3l += nor[i] * Vel(st3l)[i];
	    pt1[i] = pt[i] + dnr*nor[i];
	    pt3[i] = pt[i] + dnl*nor[i];
	    pta[i] = pt[i] + W[i]*dt;
	}
	eval_gravity(pt3,time,gv3);
	eval_gravity(pt1,time,gv1);
	eval_gravity(pta,time+dt,gva);
	for (g1 = 0.0, g3 = 0.0, i = 0; i < dim; ++i)
	{
	    g1 += 0.5*(gva[i]+gv1[i])*nor[i];
	    g3 += 0.5*(gva[i]+gv3[i])*nor[i];
	}
	u1r += g1*dt;
	u3l += g3*dt;

#if defined(DEBUG_W_SPEED)
	if (debugging("cmph"))
	{
	    (void) printf("Input data into contact_moc_plus_rh()\n");
	    verbose_print_state("st3l",st3l);
	    verbose_print_state("st2l",st2l);
	    verbose_print_state("st2r",st2r);
	    verbose_print_state("st1r",st1r);
	    (void) printf("pm = %g, pjump = %g\n",pm,pjump);
	    (void) printf("dnl = %g, dnr = %g, dt = %g\n",dnl,dnr,dt);
	    (void) printf("g1 = %g, g3 = %g\n",g1,g3);
	    print_general_vector("nor = ",nor,dim,"\n");
	}
#endif /* defined(DEBUG_W_SPEED) */


	Cmprh.rfac = 0.0;
	if (is_rotational_symmetry() && alpha > 0.0)
	{
	    double c3l, c1r;
	    double rmin, rad, radl, radr;

	    rmin = pos_radius(0.0,gr);
	    radl = pos_radius(pt3[0],gr);
	    radr = pos_radius(pt1[0],gr);
	    rad = pos_radius(pta[0],gr);

	    c3l = sound_speed(st3l);
	    c1r = sound_speed(st1r);

	    rmin = fabs(rmin);
	    if ((fabs(rad) > rmin) && (fabs(radl) > rmin) && 
	    	(fabs(radr) > rmin))
	    {
	        u3l *= (1.0 - 0.5*nor[0]*dt*c3l*alpha/radl );
	        u1r *= (1.0 + 0.5*nor[0]*dt*c1r*alpha/radr );
	        Cmprh.rfac = 0.5*nor[0]*dt*alpha/rad;
	    }
	}
	Cmprh.u1r = u1r;
	Cmprh.u3l = u3l;

	x0 = pm;
	delta = max(meps, x0*EPS);
	(void) fC(x0, &epsilon, (POINTER) &Cmprh);
	epsilon = fabs(epsilon)*EPS;
	epsilon = max(epsilon, meps);
	a = 0.5*x0;/*TOLERANCE*/
	b = 1.5*x0;/*TOLERANCE*/
#if defined(DEBUG_W_SPEED)
	if (debugging("cmph"))
	    (void) printf("Initial search interval = (%g, %g)\n",a,b);
#endif /* defined(DEBUG_W_SPEED) */
	if (find_root(fC,(POINTER) &Cmprh,0.0,&pi,a,b,epsilon,delta) ==
	                                        FUNCTION_FAILED)
	{
	    double du;
	    /* Check for vacuum production */

#if defined(DEBUG_W_SPEED)
	    if (debugging("cmph"))
	        (void) printf("First find_root() failed\n");
#endif /* defined(DEBUG_W_SPEED) */

	    (void) fC(p_min,&du,(POINTER) &Cmprh);

	    a = p_min;
	    b = 10.0*max(Cmprh.p3l,Cmprh.p1r);
	    if (du <= 0.0)
	    {
	        ul = vel(0,ansl);
	        ur = vel(0,ansr);
	        u = (Dens(ansr) > Dens(ansl)) ? ur : ul;
	    }
	    else if ( (find_root(fC,(POINTER) &Cmprh,0.0,&pi,
	                         a,b,epsilon,delta) == FUNCTION_SUCCEEDED)
	        ||
	              (search_harder_for_root(fC,(POINTER) &Cmprh,0.0,&pi,a,b,
					      &a,&b,p_min,HUGE_VAL,mnth,epsilon,
					      delta) == FUNCTION_SUCCEEDED))
	    {
#if defined(DEBUG_W_SPEED)
	        if (debugging("cmph"))
	            (void) printf("Second find_root() succeeded\n");
#endif /* defined(DEBUG_W_SPEED) */
	        ul = vel(0,ansl);
	        ur = vel(0,ansr);
	        if (pressure(ansl) <= p_min || pressure(ansr) <= p_min)
	        {
	            u = (Dens(ansr) > Dens(ansl)) ? ur : ul;
	        }
	        else
	        {
	            u = ul = ur = 0.5*(ul + ur);
	        }
	    }
	    else
	    {
	        ul = vel(0,ansl);
	        ur = vel(0,ansr);
	        u = ul = ur = 0.5*(ul + ur);
	        if (debugging("cmph"))
	        {
	            (void) printf("WARNING in contact_moc_plus_rh(), "
	                          "Failure in fC() for contact\n");
	            (void) printf("x0 = %g, a = %g, b = %g\n",x0,a,b);
	            (void) fC(pi,&du,(POINTER) &Cmprh);
	            (void) printf("pi = %g, du(pi) = %g\n",pi,du);
	            (void) printf("Velocities at vacuum, ");
	            (void) printf("ul = %g, ur = %g, ur - ul = %g\n",
	                          ul,ur,ur-ul);
	            print_function_values(fC,(POINTER) &Cmprh,0.0,a,b,100,
	                                  "fC",stdout);
	        }
	    }
	} 	
	else
	{
	    ul = vel(0,ansl);
	    ur = vel(0,ansr);
	    if (pressure(ansl) <= p_min || pressure(ansr) <= p_min)
	    {
	        u = (Dens(ansr) > Dens(ansl)) ? ur : ul;
	    }
	    else
	    {
	        u = ul = ur = 0.5*(ul + ur);
	    }
	}

	/* This assumes state is such as velocities are stored */
	for (i = 0; i < dim; ++i)
	{
	    Vel(ansl)[i] = Vel(st2l)[i] + (ul - u2l)*nor[i];
	    Vel(ansr)[i] = Vel(st2r)[i] + (ur - u2r)*nor[i];
	}
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* mass fraction of st3l is preserved */
            if(Params(st3l)->n_comps != 1)
            {
                for(i = 0; i < Params(st3l)->n_comps; i++)
                    pdens(ansl)[i] = (pdens(st3l)[i]/Dens(st3l))*Dens(ansl);
            }
        }
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* mass fraction of st1r is preserved */
            if(Params(st1r)->n_comps != 1)
            {
                for(i = 0; i < Params(st3l)->n_comps; i++)
                    pdens(ansr)[i] = (pdens(st1r)[i]/Dens(st1r))*Dens(ansr);
            }
        }
	set_state(ansl,GAS_STATE,ansl);
	set_state(ansr,GAS_STATE,ansr);

	/* Use centered in time difference for wave speed */

	for (i = 0; i < dim; ++i)
	    W[i] = 0.5*(W[i] + nor[i] * u);

#if defined(DEBUG_W_SPEED)
	if (debugging("cmph"))
	{
	    verbose_print_state("ansl ",ansl);
	    verbose_print_state("ansr ",ansr);
	    print_general_vector("W = ",W,dim,"\n");
	}
	debug_print("cmph","Left contact_moc_plus_rh()\n");
#endif /* defined(DEBUG_W_SPEED) */
}	/*end contact_moc_plus_rh*/

/* function for contact */

LOCAL boolean fC(
	double		pr,
	double		*fans,
	POINTER		prm)
{
	C_MPRH		*cmprh = (C_MPRH *) prm;
	double		rcr, rcl;
	double		pl, ul, ur, kl, kr;
	Locstate	ansl = cmprh->ansl, ansr = cmprh->ansr;

	state_on_adiabat_with_pr(cmprh->st2r,pr,ansr,TGAS_STATE);
	rcr = 0.5 * (cmprh->rc1r + acoustic_impedance(ansr));

	pl = floor_pressure(pr + cmprh->pjump,cmprh->p_min);
	state_on_adiabat_with_pr(cmprh->st2l,pl,ansl,TGAS_STATE);
	rcl = 0.5 * (cmprh->rc3l + acoustic_impedance(ansl));
	
	kl = kr = 1.0;
	if (is_rotational_symmetry() && cmprh->alpha > 0.0)
	{
	    double cl, cr;
	    double rfac = cmprh->rfac;

	    cl = sound_speed(ansl);
	    kl = 1.0/(1.0 + cl*rfac);
	    cr = sound_speed(ansr);
	    kr = 1.0/(1.0 - cr*rfac);
	}

	Vel(ansl)[0] = ul = kl*(cmprh->u3l - (pl - cmprh->p3l)/rcl);
	Vel(ansr)[0] = ur = kr*(cmprh->u1r + (pr - cmprh->p1r)/rcr);

	*fans = ul - ur;

#if defined(DEBUG_W_SPEED)
	if (debugging("cmph"))
	{
	    (void) printf("In fC(), pl = %g, pr = %g, ",pl,pr);
	    (void) printf("ul = %g, ur = %g, ul - ur = %g\n",ul,ur,ul-ur);
	}
#endif /* defined(DEBUG_W_SPEED) */

	return FUNCTION_SUCCEEDED;
}	/*end fC*/


LOCAL	void	IncomingStatesAtContact(
	WSSten	      *sten,
	double	      cl,
	double	      cr,
	Locstate      st_l3,
	Locstate      st_r1,
	double	      *dnl,
	double         *dnr,
	NptWSpeedOpts *opts)
{
	double		*pt = sten->coords;
	Locstate	 sll =  sten->sl[1];
	Locstate	  sl =  sten->sl[0];
	Locstate	  sr =  sten->sr[0];
	Locstate	 srr =  sten->sr[1];
	double		dn = sten->dn;
	double		dt = sten->dt;
	double		dn1, dn3;
	Locstate	slll, srrr;

	dn1 = cr * dt / dn;
	*dnr = cr * dt;
	dn3 = 1.0 - cl * dt / dn;
	*dnl = -cl*dt;
	if (sten->nsts > 2)
	{
	    /* Find states at feet of characteristic */

	    slll = sten->sl[2];
	    srrr = sten->sr[2];
	    ws_interpolate(st_r1,dn1,POSITIVE_SIDE,TGAS_STATE,sten);
	    ws_interpolate(st_l3,dn3,NEGATIVE_SIDE,TGAS_STATE,sten);

	    if(debugging("tst_nan"))
	    {
		printf("#IncomingStatesAtContact %d\n", sten->nsts);
	        verbose_print_state("st_l3",st_l3);
	        verbose_print_state("st_r1",st_r1);
	    }
	    if (opts->scalar_moc != RIEMANN)
	    {
	        if ((pressure_switch_on(sr,srr,srrr) == NO) ||
	            (pressure_switch_on(sl,sll,slll) == NO))
	            opts->scalar_moc = RIEMANN;
	    }
	}
	else
	{
	    static Locstate wk_slll = NULL, wk_srrr = NULL;
	    Front	    *front = sten->front;
	    Wave	    *wave = sten->wave;
	    double	    *crdsll = sten->lcrds[1];
	    double	    *crdsrr = sten->rcrds[1];
	    double           crdslll[3], crdsrrr[3];
	    double	    *nor = sten->nor;
	    COMPONENT	    pcomp = sten->pcomp;
	    COMPONENT	    ncomp = sten->ncomp;
	    HYPER_SURF	    *hs = sten->hs;
	    size_t	    sizest = front->sizest;
	    int             i, dim = front->rect_grid->dim;

	    if (wk_slll == NULL)
	    {
		scalar(&wk_slll,sizest);
		scalar(&wk_srrr,sizest);
	    }
	    slll = wk_slll;
	    srrr = wk_srrr;
	    for (i = 0; i < dim; ++i)
	    {
		crdslll[i] = pt[i] - nor[i] * 2.0 * dn;
		crdsrrr[i] = pt[i] + nor[i] * 2.0 * dn;
	    }
	    hyp_solution(crdslll,ncomp,hs,NEGATIVE_SIDE,front,wave,slll,sl);
	    hyp_solution(crdsrrr,pcomp,hs,POSITIVE_SIDE,front,wave,srrr,sr);
	    if (pressure_switch_on(sr,srr,srrr))
		interpolate_states(front,1.0-dn1,dn1,pt,sr,crdsrr,srr,st_r1);
	    else
	    {
	        copy_state(st_r1,srr);
		opts->scalar_moc = RIEMANN;
	    }
	    set_state(st_r1,TGAS_STATE,st_r1);
	    if (pressure_switch_on(sl,sll,slll))
		interpolate_states(front,1.0-dn3,dn3,crdsll,sll,pt,sl,st_l3);
	    else
	    {
		copy_state(st_l3,sll);
		opts->scalar_moc = RIEMANN;
	    }
	    set_state(st_l3,TGAS_STATE,st_l3);
	}
}		/*end IncomingStatesAtContact*/

LOCAL   boolean pressure_switch_on(
	Locstate s1,
	Locstate s2,
	Locstate s3)
{
	double p1,p2,p3;
	p1 = pressure(s1);
	p2 = pressure(s2);
	p3 = pressure(s3);

	return (fabs(p2-p1) > fabs(p3-p2) || (p1-p2)*(p2-p3) < 0.0) ? NO : YES;
}		/* end pressure_switch_on */


LOCAL	double	wall_limiter(
	WSSten        *wssten,
	SIDE          side,
	NptWSpeedOpts *opts)
{
	double    alpha;
	double    *nor = wssten->nor;
	const double    eps = 10.0*MACH_EPS;
	double    wl = opts->Wall_limiter;
	int      dim = wssten->front->rect_grid->dim;
	int      i, j, nsts = wssten->nsts;
	Locstate *st = (side == POSITIVE_SIDE) ? wssten->sr : wssten->sl;
	static   double *p, *rho, *e, *c, **v, *vn, *s;
	static   int   N = 0;

	if (wl == HUGE_VAL)
	    alpha = 1.0;
	else if (wl <= 0.0)
	    alpha = 0.0;
	else
	{
	    if ((p == NULL) || (N < nsts))
	    {
	        if (p != NULL)
	            free_these(7,p,rho,e,c,v,vn,s);
	        N = nsts;
	        uni_array(&p,N,FLOAT);
	        uni_array(&rho,N,FLOAT);
	        uni_array(&e,N,FLOAT);
	        uni_array(&c,N,FLOAT);
	        bi_array(&v,N,3,FLOAT);
	        uni_array(&vn,N,FLOAT);
	        uni_array(&s,N,FLOAT);
	    }
    
	    for (i = 0; i < nsts; ++i)
	    {
	        p[i] = pressure(st[i]);
	        rho[i] = Dens(st[i]);
	        e[i] = specific_internal_energy(st[i]);
	        c[i] = sound_speed(st[i]);
	        (void) VelocityVector(st[i],v[i]); 
	        s[i] = mag_vector(v[i],dim);
	        vn[i] = scalar_product(v[i],nor,dim);
	    }
	    alpha = 0.0;
	    for (i = 0; i < (nsts-1); ++i)
	    {
	        for (j = i+1; j < nsts; ++j)
	        {
	            double x;
    
		    x = fabs(p[j] - p[i])/((j-i)*(p[i] + eps));
		    alpha = max(alpha,x);
		    x = fabs(rho[j] - rho[i])/((j-i)*(rho[i] + eps));
		    alpha = max(alpha,x);
		    x = fabs(e[j] - e[i])/((j-i)*(e[i] + eps));
		    alpha = max(alpha,x);
		    x = fabs(vn[j] - vn[i])/((j-i)*(c[i] + eps));
		    alpha = max(alpha,x);
		    x = fabs(s[j] - s[i])/((j-i)*(c[i] + eps));
		    alpha = max(alpha,x);
		    x = vector_product(v[i],v[j],NULL,dim)/(s[i]*s[j] + eps);
		    alpha = max(alpha,x);
	        }
	    }
	    alpha *= wl;
	    alpha = min(alpha,1.0);
	}
	return alpha;
}		/*end wall_limiter */

LOCAL	boolean invalid_shock(
	double pa,
	double rhoa,
	double pb,
	double rhob)
{
	return (((pb - pa)*(rhob - rhoa) < 0.0) || (rhob < 0.0)) ? YES : NO;
}

#if defined(COMBUSTION_CODE)

LOCAL	void CJ_state_behind_curved_wave(
	Locstate	ahead,
	Locstate	ans,
	double		pjump,
	double		*Wx,
	int		st_type_ans,
	int		family)
{
	double	    pCJ, uCJ, rhoCJ, speedCJ;
	double	    corr_fac=1.0;
	static const double BFF	= 0.3;	/* TOLERANCE */

	debug_print("w_speed","Entered CJ_state_behind_curved_wave()\n");
	speedCJ = CJ_det(ans,st_type_ans,ahead,family);
	pCJ = pressure(ans);
	rhoCJ = Dens(ans);
	uCJ = vel(0,ans);

	if ((Params(ans)->speed_corr * fabs(pjump)) > (BFF * speedCJ))
	{
	    corr_fac = BFF * speedCJ / (Params(ans)->speed_corr * pjump);
	}
	pjump *= corr_fac;
	pCJ -= Params(ans)->p_corr * pjump;
	uCJ -= Params(ans)->u_corr * pjump;
	rhoCJ -= Params(ans)->rho_corr * pjump;
	speedCJ -= Params(ans)->speed_corr * pjump;

	Dens(ans) = rhoCJ;
	Vel(ans)[0] = uCJ;
	Press(ans) = pCJ;
	set_type_of_state(ans,TGAS_STATE);

	set_state(ans,st_type_ans,ans);

	Set_other_params(ans,ahead);
	*Wx = speedCJ;
	debug_print("w_speed","Left CJ_state_behind_curved_wave()\n");
	return;
}	/*end CJ_state_behind_curved_wave*/


LOCAL   void flame_moc_plus_rh(
        double     *pt,
        Locstate  sl,
        Locstate  sr,
        Locstate  st1,
        Locstate  st2,
        Locstate  st3,
        Locstate  ansl,
        Locstate  ansr,
        double     pjump,
        double     dn1,
        double     dn2,
        double     dn3,
        double     dn,
        double     *nor,
        double     *W,
        double     cl,
        double     cr,
	double     *g,
	double     ql,
	double     qr,
        Front     *front)
{
        RECT_GRID *gr = front->rect_grid;
        double     u1, u2, u3, usl, usr, ull, urr;
	double     c1, c3;
	double 	  g1dt,g2dt,g3dt;
       	double     pt1[3], pt2[3], pt3[3];
	const double 	  *gv1;
	const double 	  *gv2;
	const double 	  *gv3;
	double     dt = front->dt;
        double     time = front->time;
        int       i, dim = gr->dim;
	
	double     ul, ur;
	double     rl, rr;
	double     pl, pr;
	double     s, LHS, gam, V[MAXD];
	double     newpl,newpr;
	unsigned int iter=0;
	unsigned int maxiter=100;
			
	double rmin, radl, radr2, radr3;
	double ulfac = 0.0;
	double urfac = 0.0;
	double vflame;
	Locstate unburned_state = NULL, burned_state = NULL;	

	/* Debugging variables */
	double mu2,heat,tmp1,tmp2,w1,w2,s1,M,M1,M2,M3;
	FILE *fpt;
	/* Debugging variables */

	double	  alpha = rotational_symmetry();

	/* Changed code here, Dec 9, 2003 */
	u1 = 0.0;
	u2 = 0.0;
	u3 = 0.0;
	usl = 0.0;
	usr = 0.0;

	for(i=0; i<dim; i++)
	{
	    u1 += nor[i]*Vel(st1)[i];
	    u2 += nor[i]*Vel(st2)[i];
	    u3 += nor[i]*Vel(st3)[i];
	    usl += nor[i]*Vel(sl)[i];
	    usr += nor[i]*Vel(sr)[i];
	    pt1[i] = pt[i] - (1.- dn1)*dn*nor[i];
	    if(Burned(sl))
		pt2[i] = pt[i] + dn2*dn*nor[i];		
	    else
		pt2[i] = pt[i] - (1. - dn2)*dn*nor[i];
	    pt3[i] = pt[i] + dn3*dn*nor[i];
	}

	if(is_gravity() == YES) 
	{
	    gv1 = gravity(pt1,time);
	    gv2 = gravity(pt2,time);			
	    gv3 = gravity(pt3,time);
	    g1dt = scalar_product(gv1,nor,dim)*dt;
	    g3dt = scalar_product(gv3,nor,dim)*dt;
	}
	else
	{
	    g1dt = 0.0;
	    g3dt = 0.0;
	}
	
	ulfac = -g1dt;
	urfac = +g3dt;	
	c1 = sound_speed(st1);
 	c3 = sound_speed(st3);

	Set_params(ansl,sl);
	set_type_of_state(ansl,TGAS_STATE);
	Set_params(ansr,sr);
	set_type_of_state(ansr,TGAS_STATE);

	if(is_rotational_symmetry() && alpha > 0.0)
	{
	    rmin = pos_radius(0.0,gr);		
	    radl = pos_radius(pt1[0],gr);		 
	    radr3 = pos_radius(pt3[0],gr);
			
	    if(fabs(radl) > fabs(rmin))
		ulfac += nor[0]*alpha*c1*u1*dt/radl;
	    if(fabs(radr3) > fabs(rmin))
		urfac += nor[0]*alpha*c3*u3*dt/radr3;
	}


	if(flame_moc_plus_rh_signed(st1,st2,st3,ansl,ansr,YES,ulfac,
				urfac,ql-qr,nor,&ull,&urr,dim)==NO)
        {   
	    if(flame_moc_plus_rh_signed(st1,st2,st3,ansl,ansr,NO,ulfac,
				    urfac,ql-qr,nor,&ull,&urr,dim)==NO)
	    {	
	       	printf("ERROR: No root found in flame_moc_plus_rh\n");
	       	clean_up(ERROR);
	    }
	}


	/* if ansr is unburned, set unburned_state = ansr, and 
	 * burned_state = ansl; else vice_versa. This is done 
	 * because in combustion, burned and unburned states have 
	 * differente equations and some equations are solved with 
	 * respect to unburned state and some with respect to burned.
	*/	
	if((Burned(ansl) == YES))
	{
	    /* Forward-moving flame */
	    unburned_state = ansr;
	    burned_state = ansl;
	}
	else
	{
	    /* Backward-moving flame */
	    unburned_state = ansl;
	    burned_state = ansr;
	}

	for(i=0; i<dim; i++)
	{
		 Vel(ansl)[i] = Vel(sl)[i] + (ull - usl)*nor[i];
		 Vel(ansr)[i] = Vel(sr)[i] + (urr - usr)*nor[i];
	}

	Set_params(ansl,sl);
	Set_params(ansr,sr);
	set_type_of_state(ansl,TGAS_STATE);
	set_type_of_state(ansr,TGAS_STATE);
	set_state(ansl,GAS_STATE,ansl);
	set_state(ansr,GAS_STATE,ansr);
	
}       /*end flame_moc_plus_rh*/ 


LOCAL boolean flame_moc_plus_rh_signed(
	Locstate	st1,
	Locstate	st2,
	Locstate	st3,
	Locstate	ansl,
	Locstate	ansr,
	boolean	sign,
	double	ulfac,
	double	urfac,
	double	heat,
	double	*nor,
	double*	pull,
	double*       purr,
	int	        dim)
{
        double pneg, ppos, delta;
        double gneg, gpos;
	double beta;
	double alpha = Press(st2);
	double pl,pr,rhol,rhor;
 
	/*if (debugging("my_point"))
	{
	    int i;
	    double pl,p,g;
	    FILE *file;
	    char fname[100];

	    strcpy(fname,"root");
	    sprintf(fname, "%s%i", fname,fc);
	    printf("fc %d fname %s\n",fc,fname);
	    fc++;
	    pl = 0;
	    delta = max(0.01*fabs(alpha), EPS);
	    file = fopen(fname,"w");
	    for (i = 0; i < 1000; i++)
	    {
		p = pl+i*delta;
		g = flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,p,ulfac,
		                        urfac,heat,nor,pull,purr,dim);
		fprintf(file,"%f	%f\n",p,g);
	    }
	    fclose(file);
	    
	}*/
        pneg = alpha;
        gneg = flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,pneg,ulfac,
			urfac,heat,nor,pull,purr,dim);
	/* Changed code here - 9th Dec, 2003 */
	if (fabs(gneg)<EPS)
	    return YES;
        
        delta = max(0.01*fabs(alpha), EPS); 
        do {
            ppos = pneg + ((sign == YES)?+delta:-delta);
            if(ppos<0)
            	ppos = EPS;
                      
            gpos = flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,ppos,
	    		ulfac,urfac,heat,nor,pull,purr,dim);

            if ((gneg*gpos)>0.0)
 	    {
		/*TMP */
		if (fabs(gpos) > fabs(gneg))
			return NO;
		else
                pneg = ppos;
                gneg = gpos;
                if ((ppos<=0.1*alpha)||(ppos>=10.*alpha))
                { 
		    if (sign == NO)     
		    {
		        printf("Warning: iteration out of range!\n");
		        printf("alpha %f ppos = %f\n",alpha,ppos);
			clean_up(ERROR);
			flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,
				alpha,ulfac,urfac,heat,nor,pull,purr,dim); 
			return YES;
		    }
		    else
			return NO;
		    }                             
                }
 
        } while((fabs(gpos)>EPS)&&((gneg*gpos)>0.0));
	if (fabs(gpos)<EPS)
	    return YES;

	/* Find ansl and ansr */
	beta = flame_moc_plus_rh_findstates(pneg,ppos,st1,st2,st3,ansl,
			ansr,ulfac,urfac,heat,nor,pull,purr,dim);
	/* Find ansl and ansr */
	/*TMP check if ans if physical */
        pl = Press(ansl);
	pr = Press(ansr);
	rhol = Dens(ansl);
	rhor = Dens(ansr);
 	if (pl < 0 || pr < 0 || rhol < 0 || rhor < 0)
	{
	    printf("Warning: unphysical state, turn around\n");
	    printf("pl %f pr %f rhol %f rhor %f gpos %f gneg %f\n",pl,pr,rhol,rhor,gpos,gneg);
	    if (sign == NO)
		clean_up(ERROR);
	    else
	        return NO;
	}

        return YES;
}

LOCAL double flame_moc_plus_rh_guessp(
	Locstate	st1,
	Locstate	st2,
	Locstate	st3,
	Locstate	ansl,
	Locstate	ansr,
	double	p,
	double	ulfac,
	double	urfac,
	double	heat,
	double	*nor,
	double*	pull,
	double* 	purr,
	int		dim)
{
	double pl,pr,rl,rr,vl,vr,g;
	double u1, u2, u3;
	int i;	
	double vflame;
	double gam = adiabatic_gamma(st1);
	double mu2 = (gam - 1.)/(gam + 1.);

	double Pr_ansl;

	u1 = u2 = u3 = 0.0;
	for (i=0; i<dim; i++)
	{
	    u1  += nor[i]*Vel(st1)[i];
	    u2  += nor[i]*Vel(st2)[i];
	    u3  += nor[i]*Vel(st3)[i];
	}
	/* if ansr is unburned, set unburned_state = ansr, and 
	 * burned_state = ansl; else vice_versa. This is done 
	 * because in combustion, burned and unburned states have 
	 * differente equations and some equations are solved with 
	 * respect to unburned state and some with respect to burned.
	*/	
	if ((Burned(ansl)))
	{
	    /* Forward-moving flame */
	    /* Changed code here  - 9th December, 2003*/
	    Press(ansr)     = p;
	    *purr = u3 + urfac + (Press(ansr) - 
	    		Press(st3))/(Dens(st3)*sound_speed(st3));
	    Dens(ansr) = Dens(st2) + (Press(ansr) - 
	    		Press(st2))/(sound_speed(st2)*sound_speed(st2));
		
	    vflame = flame_velocity(ansr,ansr);
		
	    Press(ansl) = vflame*(u1 - ulfac - *purr) + 
	    		vflame*Press(st1)/(Dens(st1)*sound_speed(st1)) 
			+ Press(ansr)/(Dens(ansr));
	    Press(ansl)     = Press(ansl)/(vflame/(Dens(st1)*sound_speed(st1)) 
	    		+ 1/(Dens(ansr)));
		
	    *pull = u1 - (Press(ansl) - Press(st1))/(Dens(st1)*
	    		sound_speed(st1)) - ulfac;
	    Dens(ansl) = 1/(1/Dens(ansr) - (Press(ansl) - 
	    		Press(ansr))/(vflame*vflame*Dens(ansr)*Dens(ansr)));

	    g = Press(ansr)*(1/Dens(ansr) - mu2/Dens(ansl)) -
			Press(ansl)*(1/Dens(ansl) - mu2/Dens(ansr)) +
		    	2*mu2*heat;
	    reset_gamma(ansl);
	    reset_gamma(ansr);
	    return g;
	}
	else
	{
	    /* Backward-moving flame */
	    Press(ansl)     = p;
	    *pull = u1 - ulfac - (Press(ansl) - 
	    			Press(st1))/(Dens(st1)*sound_speed(st1));
	    Dens(ansl) = Dens(st2) + (Press(ansl) - 
	    			Press(st2))/(sound_speed(st2)*sound_speed(st2));
		
	    vflame = flame_velocity(ansl,ansl);
		
	    Press(ansr) = *pull - u3 - urfac + 
	    			Press(st3)/(Dens(st3)*sound_speed(st3)) + 
				Press(ansl)/(Dens(ansl)*vflame);
	    Press(ansr) = Press(ansr)/(1/(Dens(st3)*sound_speed(st3)) + 
	    			1/(Dens(ansl)*vflame));
		
	    *purr = u3 + (Press(ansr) - Press(st3))/(Dens(st3)*
	    			sound_speed(st3)) + urfac;
	    Dens(ansr) = 1/(1/Dens(ansl) + (Press(ansl) - Press(ansr))/
	    			(vflame*vflame*Dens(ansl)*Dens(ansl)));

	    g = Press(ansr)*(1/Dens(ansr) - mu2/Dens(ansl)) -
		    Press(ansl)*(1/Dens(ansl) - mu2/Dens(ansr)) +
		    2*mu2*heat;
	    reset_gamma(ansl);
	    reset_gamma(ansr);
	    return g;
	}
}

LOCAL double flame_moc_plus_rh_deltap(
	Locstate	st1,
	Locstate	st2,
	Locstate	st3,
	Locstate	ansl,
	Locstate	ansr,
	double	p,
	double	ulfac,
	double	urfac,
	double	heat,
	double	*nor,
	double*	pull,
	double*	purr,
	int	dim)
{
	double DP = max(0.0001*p, EPS);
	double dp = DP; 

        double gneg = flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,p,ulfac,
			urfac,heat,nor,pull,purr,dim);
        double gpos = flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,p+DP,
			ulfac,urfac,heat,nor,pull,purr,dim);
        double gdif = max((gpos-gneg)/DP, EPS);
            
              dp = max(fabs(gneg/gdif), EPS);
                        
        return dp;
}

LOCAL double flame_moc_plus_rh_findstates(
	double       pneg,
	double       ppos,
	Locstate	st1,
	Locstate	st2,
	Locstate	st3,
	Locstate	ansl,
	Locstate	ansr,
	double	ulfac,
	double	urfac,
	double	heat,
	double	*nor,
	double*	pull,
	double* 	purr,
	int	dim)
{
	double pmid, gmid;
 	double gneg = flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,pneg,
			ulfac,urfac,heat,nor,pull,purr,dim);
	do
	{
	     pmid = (pneg+ppos)/2;
             gmid = flame_moc_plus_rh_guessp(st1,st2,st3,ansl,ansr,pmid,
	   		ulfac,urfac,heat,nor,pull,purr,dim);

             if(gneg > 0)
	     {   
        	if(gmid > 0)
                     pneg = pmid;
        	else
                     ppos = pmid;
             } 
             else
	     {   
        	if(gmid > 0)
                     ppos = pmid;
        	else
                     pneg = pmid;
             } 

	} while((fabs(gmid)>EPS)&&(fabs(pneg-ppos)>1e-6));

	return pmid;
}
#endif /* defined(COMBUSTION_CODE) */
