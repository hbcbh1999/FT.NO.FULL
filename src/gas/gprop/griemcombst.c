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
*				griemcombst.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*    	Contains the Riemann solvers for detonation gas dynamics.
*
*	CJ_det() calculates the state behind a CJ detonation.
*/

#include <gdecs/gdecs.h>

#if defined(COMBUSTION_CODE)

enum {
	FLOW          = 0,
	IGNITE        = 1,
	ABOVE_IGNIT   = 1,
	BELOW_IGNIT   = 0,
	ERR_CHECK     = 1,
	NO_CHECK      = 0,
	DO_NEWTON     = 1,
	NO_NEWTON     = 0,
	NOT_AVAIL     = 0,
	AVAILABLE     = 1
};


LOCAL int l_CJ = NOT_AVAIL,r_CJ = NOT_AVAIL;

	/* LOCAL Function Declarations */
LOCAL	double	shock_raref_p_guess(Locstate,int,Locstate,int,
				    Locstate,Locstate,double,double,int);
LOCAL	double	shock_shock_p_guess(Locstate,Locstate);
LOCAL	double	shock_shock_state_guess(Locstate);
LOCAL	int	burn_curve(double*,double*,double,Locstate,int);
LOCAL	int	curve(double*,double*,double,Locstate,int );
LOCAL	int	detonation(double*,double*,double,Locstate,int);
LOCAL	int	mid_temp(Locstate,double,double,int,int,int);
LOCAL	int	rarefaction(double*,double*,double,Locstate,int);
LOCAL	int	shock(double*,double*,double,Locstate,int);
LOCAL	int	vacuum(double*,double*,Locstate,Locstate);
LOCAL	void	intersection(double*,double*,RIEMANN_SOLVER_WAVE_TYPE*,
                             RIEMANN_SOLVER_WAVE_TYPE*,Locstate,Locstate,
			     int,int);
LOCAL	void	get_p_guess(double*,int*,Locstate,Locstate,Locstate,Locstate,
			    int,int);
LOCAL	void	newton(double*,double*,int*,int*,double,Locstate,Locstate,
		       int (*)(double*,double*,double,Locstate,int),
		       int (*)(double*,double*,double,Locstate,int));
LOCAL	void	reference_values(Locstate,Locstate,int,int);


/*
*			combust_onedrsoln():
*
*	Determines the solution of a Riemann problem along a given
*  ray, specified by the slope "sample_speed" in the space-time
*  plane.
*/



EXPORT int combust_onedrsoln(
	double		sample_speed,
	Locstate	sl,
	Locstate	sr,
	Locstate	ans,
	double		*spdnew,
	int		state_type)
{
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave, s_or_r;
	int		l_or_r = ERROR;         /* left or right of contact */
	double		side;			/* left = -1.0 right = 1.0 */
	double		vr,vl;			/* normal velocities */
	double		vstart,vxstart;
	double		vxl,vxr;		/* x-velocity on each side */
	double		mr,ml;   		/* midstate fluxes */
	double		m_start;
	double		shock_speed;
	double		cside, cmid;
	double		cCJ, MCJ;
	static boolean	first = YES;
	static size_t	sizest;

	static Locstate s_start = NULL,Tsl = NULL,
	                Tsr = NULL,mid,CJ = NULL;	/* TGas states */

	debug_print("riem_sol","Entering combust_onedrsoln\n");

	if ((is_obstacle_state(sl)) && (is_obstacle_state(sr)))
	{
		g_obstacle_state(ans,g_sizest());
		*spdnew = 0.0;
		return l_or_r;
	}

	if (first)
	{
		first = NO;
		Gas_param *params = (is_obstacle_state(sl) != NULL) ?
			Params(sl) : Params(sr);
		sizest = params->sizest;
		(*params->_alloc_state)(&Tsl,max(sizest,sizeof(VGas)));
		(*params->_alloc_state)(&Tsr,max(sizest,sizeof(VGas)));
		(*params->_alloc_state)(&mid,sizest);
		(*params->_alloc_state)(&CJ,sizest);
		(*params->_alloc_state)(&s_start,sizest);
	}


	set_state(Tsl,TGAS_STATE,sl);
	set_state(Tsr,TGAS_STATE,sr);

	vxr = Vel(Tsr)[0];
	vxl = Vel(Tsl)[0];


	vl = Vel(Tsl)[0];
	vr = Vel(Tsr)[0];

#if defined(DEBUG_GRIEMANN)
	if(debugging("riem_sol"))
	{
		(void) printf("in combust_onedrsoln() ");
		(void) printf("before combust_find_mid_state\n");
		verbose_print_state("sl",sl);
		verbose_print_state("sr",sr);
	}
#endif /* defined(DEBUG_GRIEMANN) */

	(void) combust_find_mid_state(Tsl,Tsr,&Press(mid),&Press(mid),
		                      &Vel(mid)[0],&Vel(mid)[0],
				      &ml,&mr,&l_wave,&r_wave);
	reset_gamma(mid);

#if defined(DEBUG_GRIEMANN)
	if (debugging("riem_sol"))
	{
		(void) printf("combust_find_mid_state returns ");
		(void) printf("to riemann_solution\n");
		(void) printf("pm=%g vm=%g ml=%g mr=%g\n",
			      Press(mid),Vel(mid)[0],ml,mr);
		(void) printf("lwave=%d rwave=%d \n",l_wave,r_wave);
	}
#endif /* defined(DEBUG_GRIEMANN) */
	
	if (sample_speed < Vel(mid)[0])
	{
		m_start      	= ml;
		vstart		= vl;
		ft_assign(s_start,Tsl,sizest);
		s_or_r       	= l_wave;
		l_or_r       	= LEFT_FAMILY;
		side		= -1.0;
		vxstart		= vxl;
	}
	
	else
	{
		m_start      	= mr;
		vstart		= vr;
		ft_assign(s_start,Tsr,sizest);
		s_or_r       	= r_wave;
		l_or_r       	= RIGHT_FAMILY;
		side		= 1.0;
		vxstart		= vxr;
	}

#if defined(DEBUG_GRIEMANN)
	debug_print("riem_sol","In riem_sol wave=%d side=%d \n",s_or_r,l_or_r);
#endif /* defined(DEBUG_GRIEMANN) */

	if ((s_or_r == SHOCK) || (s_or_r == STRONG_DET))
	{
			/* shock or strong detonation */

		shock_speed = vstart + (side * m_start/Dens(s_start));
#if defined(DEBUG_GRIEMANN)
		if (debugging("znd_speed"))
		{
			if (Press(mid) > (1.5 * pressure(s_start)))
			{
				(void) printf("shock speed = %g\n",shock_speed);
			}
		}
#endif /* defined(DEBUG_GRIEMANN) */

		if ((side * (sample_speed - shock_speed)) >= 0.)
		{
			*spdnew = fabs(Vel(ans)[0]) + sound_speed(s_start);

			set_state(ans,state_type,s_start);
			return l_or_r;
		}
		else
		{
			Dens(ans) = side * m_start/(shock_speed - Vel(mid)[0]);
			Vel(ans)[0] = Vel(mid)[0];
			Set_params(ans,s_start);
			set_type_of_state(ans,TGAS_STATE);
			if (s_or_r == STRONG_DET)
				Set_other_params(ans,s_start);
			Press(ans) = Press(mid);	
		        reset_gamma(ans);
				
			if (Composition_type(ans) == ZND)
				React(ans) = React(s_start);
			*spdnew = fabs(Vel(ans)[0]) + sound_speed(ans);
			if (state_type != state_type(ans)) 
				set_state(ans,state_type,ans);
			return l_or_r;
		}
	}
	else if (s_or_r == CJ_DET)  /* CJ-detonation */
	{
		Vel(s_start)[0] = vstart;
		shock_speed = CJ_det (CJ,TGAS_STATE,s_start,l_or_r);
		cCJ = sound_speed(CJ);
		MCJ = Dens(CJ) * cCJ;
		shock_speed = vstart + (side * MCJ / Dens(s_start));
		Vel(s_start)[0] = vxstart;

		if ((side * (sample_speed - shock_speed)) >= 0.)
		{
			*spdnew = fabs(Vel(ans)[0]) + sound_speed(s_start);

			set_state(ans,state_type,s_start);
			return l_or_r;
		}
		else	
		{
			Press(s_start) = Press(CJ);
			Vel(s_start)[0] = Vel(CJ)[0];
			vstart = Vel(CJ)[0];
			Dens(s_start) = Dens(CJ);
			cside = cCJ;
		        reset_gamma(s_start);
				
			Set_params(s_start,CJ);
			set_type_of_state(s_start,TGAS_STATE);
		}
	}

	else if (s_or_r == RAREFACTION)
	{
		cside = sound_speed(s_start);
		if ((side * (sample_speed - vstart)) >= cside)
		{
			*spdnew = fabs(Vel(ans)[0]) + sound_speed(s_start);

			set_state(ans,state_type,s_start);
			return l_or_r;
		}
	}

	/*   between rarefaction and slip line */

	state_on_adiabat_with_pr(s_start,Press(mid),mid,TGAS_STATE);
	cmid = sound_speed(mid);

	if ((side * (sample_speed - Vel(mid)[0])) <= cmid)
	{
		*spdnew = fabs(Vel(mid)[0]) + sound_speed(mid);
		set_state(ans,state_type,mid);
	}
	else
	{

	/*   last case : point is in rarefaction fan  */

		(void) oned_state_in_rarefaction_fan(sample_speed,
			        Vel(s_start)[0],s_start,mid,ans,
				state_type,spdnew,l_or_r);
	}
	Set_params(ans,s_start);

	return l_or_r;
}		/*end combust_onedrsoln*/


/*
*			combust_find_mid_state():
*
* 	This routine does the logic and runs the procedure to find the mid
*  state.  It finds the intersection of the v(p) - left and v(p) - right
*  curves for burned and unburned v(p) curves.
*  We only burn if forced to burn. 
*
* 	The output is pmid, vmid, l_wave and r_wave (mid_state pressure
*  and velocity and the type of wave on left and right, shock, 
*  rarefaction, strong detonation or CJ-detonation).   
*  Also, found are the fluxes across the wave on each side.
*
* 	The input is the TGas states l_state and r_state.
*/

EXPORT boolean combust_find_mid_state(
	Locstate	         l_state,
	Locstate	         r_state,
	double		         *pmidl,
	double		         *pmidr,
	double		         *vmidl,
	double		         *vmidr,
	double		         *left_flux,
	double		         *right_flux,
	RIEMANN_SOLVER_WAVE_TYPE *l_wave,
	RIEMANN_SOLVER_WAVE_TYPE *r_wave)
{
	int        burn_left, burn_right;
	static boolean    first = YES;
	double        cl, cr;
	double        vright_min,vleft_max;
	double        pmid, vmid;
	static Locstate CJ = NULL;

	if (is_obstacle_state(l_state) || is_obstacle_state(r_state))
	{
	    screen("ERROR: obstacle state in combust_find_mid_state\n");
	    clean_up(ERROR);
	}
	if (first)
	{
	    first = NO;
	    (*Params(l_state)->_alloc_state)(&CJ,Params(l_state)->sizest);
	}

#if defined(DEBUG_GRIEMANN)
	if(debugging("find_mid_state"))
	{
	    verbose_print_state("In combust_find_mid_state  l_state",l_state);
	    verbose_print_state("r_state",r_state);
	}
#endif /* defined(DEBUG_GRIEMANN) */

	l_CJ = NOT_AVAIL;
	r_CJ = NOT_AVAIL;

	if (Composition_type(l_state) == PTFLAME || 
	   (Composition_type(l_state) == THINFLAME))
	{
	    burn_left = (Burned(l_state)) ? BURNED : UNBURNED;
	    burn_right = (Burned(r_state)) ? BURNED : UNBURNED;
	}

	cl = sound_speed(l_state);
	cr = sound_speed(r_state);


	if ((fabs((Press(l_state) - Press(r_state)) / (Press(r_state) + 
	    Press(l_state))) < Tol_pressure(l_state))
	                    &&
	    (fabs((Vel(l_state)[0] - Vel(r_state)[0]) /
	      (0.000001 + fabs(Vel(l_state)[0]) + fabs(Vel(r_state)[0]))) < /*TOLERANCE*/
	                             Tol_pressure(r_state)))
	{
	    *pmidl = *pmidr = .5 * (Press(l_state) + Press(r_state));
	    *vmidl = *vmidr = .5 * (Vel(l_state)[0] + Vel(r_state)[0]);

	    *l_wave = RAREFACTION;
	    *r_wave = RAREFACTION;

	    *left_flux  = Dens(l_state) * cl; 
	    *right_flux = Dens(r_state) * cr;
	    return FUNCTION_SUCCEEDED;
	}


	if (vacuum(&vright_min, &vleft_max, l_state, r_state))
	{
	    *pmidl = *pmidr = Min_pressure(r_state);
	    *vmidl = *vmidr = .5 * (vleft_max + vright_min); 

	    *l_wave = RAREFACTION;
	    *r_wave = RAREFACTION;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if ((fabs(*vmidl-Vel(l_state)[0])) <
	            (Press(l_state) * Tol_alpha(l_state)))
	        *left_flux  = Dens(l_state) * cl;
	    else
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        *left_flux = Press(l_state)/(*vmidl - Vel(l_state)[0]);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if ((fabs(*vmidr-Vel(r_state)[0])) <
	            (Tol_alpha(r_state) * Press(r_state)))
	        *right_flux = Dens(r_state) * cr;
	    else
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        *right_flux = Press(r_state)/(Vel(r_state)[0] - *vmidr);
	    return FUNCTION_SUCCEEDED;
	}

	if (Composition_type(l_state) != PTFLAME ||
	   (Composition_type(l_state) != THINFLAME))
	{
	    intersection(&pmid,&vmid,l_wave,r_wave,l_state,r_state,FLOW,FLOW);
	    *pmidl = *pmidr = pmid;
	    *vmidl = *vmidr = vmid;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if ((fabs(vmid - Vel(l_state)[0])) <= 
	        (Tol_alpha(l_state) * fabs(pmid - Press(l_state))))
	        *left_flux = Dens(l_state) * cl;
	    else
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        *left_flux = (pmid - Press(l_state))/(Vel(l_state)[0] - vmid);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if ((fabs(vmid - Vel(r_state)[0])) <= 
	        (Tol_alpha(r_state) * fabs(pmid - Press(r_state))))
	        *right_flux = Dens(r_state) * cr;
	    else
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        *right_flux = (Press(r_state) - pmid) /
	                (Vel(r_state)[0] - vmid);

	    return FUNCTION_SUCCEEDED;
	}

	if ((burn_left == BURNED) && (burn_right == BURNED))
	{
	    /* if both sides are burned just find the intersection */

	    intersection(&pmid,&vmid,l_wave,r_wave,l_state,r_state,FLOW,FLOW);
	    *pmidl = *pmidr = pmid;
	    *vmidl = *vmidr = vmid;
	}

	if ((burn_left == UNBURNED) && (burn_right == BURNED))
	{

	/* If left is unburned and right is burned find 
	*  intersection, check mid-left temperature to see
	*  if we need a detonation. If we do, find the new
	*  intersection and check the mid-left temp for error.
	*/

	    intersection(&pmid,&vmid,l_wave,r_wave,l_state,r_state,FLOW,FLOW);
	    *pmidl = *pmidr = pmid;
	    *vmidl = *vmidr = vmid;
	    if (mid_temp(l_state,pmid,vmid,LEFT_FAMILY,*l_wave,NO_CHECK) == 
	        ABOVE_IGNIT)
	    {
	        intersection(&pmid,&vmid,l_wave,r_wave,l_state,
	                     r_state,IGNITE,FLOW);
	        *pmidl = *pmidr = pmid;
	        *vmidl = *vmidr = vmid;
	        mid_temp(l_state,pmid,vmid,LEFT_FAMILY,*l_wave,ERR_CHECK);
	    }
	}

	if ((burn_left == BURNED) && (burn_right == UNBURNED))
	{

	/* If left is burned and right is unburned find 
	*  intersection, check mid-right temperature to see
	*  if we need a detonation. If we do, find the new
	*  intersection and check the mid-right temp for error.
	*/

	    intersection(&pmid,&vmid,l_wave,r_wave,l_state,r_state,FLOW,FLOW);
	    *pmidl = *pmidr = pmid;
	    *vmidl = *vmidr = vmid;
	    if (mid_temp(r_state,pmid,vmid,RIGHT_FAMILY,*r_wave,NO_CHECK) ==
	        ABOVE_IGNIT)
	    {
	        intersection(&pmid,&vmid,l_wave,r_wave,l_state,r_state,
	                     FLOW,IGNITE);
	        *pmidl = *pmidr = pmid;
	        *vmidl = *vmidr = vmid;
	        mid_temp(r_state,pmid,vmid,RIGHT_FAMILY,*r_wave,ERR_CHECK);
	    }
	}

	if ((burn_left == UNBURNED) && (burn_right == UNBURNED))
	{

	/* If left is unburned and right is unburned find 
	*  intersection, check mid-right temperature to see
	*  if we need a detonation on right. If we do, find the new
	*  intersection and check the mid-right temp for error.
	*  Then check to see if there is also need for a 
	*  detonation on the left. If so find the new intersection
	*  (of burn_curve and burn_curve) and check for error
	*  i.e. that mid_state temperatures are high enough.
	*/

	    intersection(&pmid,&vmid,l_wave,r_wave,l_state,r_state,FLOW,FLOW);
	    *pmidl = *pmidr = pmid;
	    *vmidl = *vmidr = vmid;
	    if (mid_temp(r_state, pmid, vmid, RIGHT_FAMILY,*r_wave,NO_CHECK) ==
	        ABOVE_IGNIT)
	    {
	        intersection(&pmid,&vmid,l_wave,r_wave,l_state,
			     r_state,FLOW,IGNITE);
	        *pmidl = *pmidr = pmid;
	        *vmidl = *vmidr = vmid;
	        mid_temp(r_state,pmid,vmid,RIGHT_FAMILY,*r_wave,ERR_CHECK);

	        if (mid_temp(l_state, pmid, vmid, LEFT_FAMILY,
		             *l_wave,NO_CHECK) == ABOVE_IGNIT)
	        {
	            intersection(&pmid,&vmid,l_wave,r_wave,l_state,
	                         r_state,IGNITE,IGNITE);
	            *pmidl = *pmidr = pmid;
	            *vmidl = *vmidr = vmid;
	            mid_temp(r_state,pmid,vmid,RIGHT_FAMILY,*r_wave,ERR_CHECK);
	            mid_temp(l_state, pmid, vmid, LEFT_FAMILY,
	                     *l_wave, ERR_CHECK);
	        }

	    }

	    else if (mid_temp(l_state, pmid, vmid,LEFT_FAMILY,
	        *l_wave, NO_CHECK) == ABOVE_IGNIT)
	    {

	        /* If there was no need for a detonation on the
	        *  right, check the left. If need be find the
	        *  intersection (burn_curve and curve) for a
	        *  left detonation. Error check the temp, test if
	        *  a right detonation is also needed. Get the new
	        *  intersection if so and check the new mid_state
	        *  temperatures support detonation
	        */

	        intersection(&pmid,&vmid,l_wave,r_wave,l_state,
			     r_state,IGNITE,FLOW);
	        *pmidl = *pmidr = pmid;
	        *vmidl = *vmidr = vmid;

	        mid_temp(l_state, pmid, vmid, LEFT_FAMILY,*l_wave, ERR_CHECK);

	        if (mid_temp(r_state, pmid, vmid,RIGHT_FAMILY,
	            *r_wave, NO_CHECK) == ABOVE_IGNIT)
	        {
	            intersection(&pmid,&vmid,l_wave,r_wave,l_state,
	                         r_state,IGNITE,IGNITE);
	            *pmidl = *pmidr = pmid;
	            *vmidl = *vmidr = vmid;
	            mid_temp(r_state,pmid,vmid,RIGHT_FAMILY,*r_wave,ERR_CHECK);
	            mid_temp(l_state,pmid,vmid,LEFT_FAMILY,*l_wave,ERR_CHECK);
	        }

	    }
	}

	if (*l_wave == CJ_DET)
	{
	    CJ_det(CJ,TGAS_STATE,l_state,LEFT_FAMILY);
	    Set_other_params(CJ,l_state);
	    *left_flux = Dens(CJ) * sound_speed(CJ);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	else if ((fabs(vmid - Vel(l_state)[0])) <=
	    (Tol_alpha(l_state) * fabs(pmid - Press(l_state))))
	{
	    *left_flux = Dens(l_state) * cl;
	}
	 else if (((fabs(vmid - Vel(l_state)[0])) <
	              (Tol_pressure(l_state) * Vel(l_state)[0]))
	     && ((fabs(pmid - Press(l_state))) < 
	     (Tol_pressure(l_state) * Press(l_state))))
	{
	    *left_flux = Dens(l_state) * cl;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	else
	{
	    *left_flux = (pmid - Press(l_state))/(Vel(l_state)[0] - vmid);
	}

	if (*r_wave == CJ_DET)
	{
	    CJ_det(CJ,TGAS_STATE,r_state,RIGHT_FAMILY);
	    Set_other_params(CJ,r_state);
	    *right_flux = Dens(CJ) * sound_speed(CJ);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	else if ((fabs(vmid - Vel(r_state)[0])) <=
	    (Tol_alpha(r_state) * fabs(pmid - Press(r_state))))
	{
	    *right_flux = Dens(r_state) * cr;
	}
	 else if (((fabs(vmid - Vel(r_state)[0])) <
	          (Tol_pressure(r_state) * Vel(r_state)[0])) &&
	     ((fabs(pmid - Press(r_state))) < 
	          (Tol_pressure(r_state) * Press(r_state))))
	{
	    *right_flux = Dens(r_state) * cr;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	else
	{ 
	    *right_flux = (Press(r_state) - pmid) / (Vel(r_state)[0] - vmid);
	}
	return FUNCTION_SUCCEEDED;
}	     /*end combust_find_mid_state*/

/*
*			mid_temp():
*
*  	This routine finds the mid-state temperature (for CJ-det finds it
*  just inside the detonation and before the rarefaction).  
*  If error_check = ERR_CHECK it makes sure the mid_temp is above
*  t_crit (if the temperature is too low the program stops and
*  an error message is printed). 
*  In any case it returns whether mid_temp is above or below ignition
*  temperature.  
* 
* 	There is no output except possibly an error message and whether
*  mid-temp is above or below ignition temperature.
*  If error_check = NO_CHECK then no error check is made. 
*/

LOCAL int mid_temp(
	Locstate	start,
	double		pmid,
	double		vmid,
	int		l_or_r,
	int		s_or_r,
	int		error_check)
{
	double		mside,shock_speed,rhomid;
	double		cstart;
	double		temp;
	double		t_crit;
	int		side;
	static Locstate CJ = NULL;
	static Locstate mid = NULL;
	static boolean	first = YES;

	if (first)
	{
		first = NO;
		(*Params(start)->_alloc_state)(&CJ,Params(start)->sizest);
		(*Params(start)->_alloc_state)(&mid,Params(start)->sizest);
	}
	if (l_or_r == LEFT_FAMILY) side = -1;
	else side = 1;
	Press(mid) = pmid;
	Vel(mid)[0]  = vmid;
        reset_gamma(mid);
		
	if ((s_or_r == SHOCK) || (s_or_r == STRONG_DET))
	{

	/* inside shock or strong detonation */

		cstart = sound_speed(start); 
#if !defined(UNRESTRICTED_THERMODYNAMICS)
		if ((fabs(vmid - Vel(start)[0])) <= 
			(Tol_alpha(start) * fabs(pmid - Press(start))))
		{
			mside = Dens(start) * cstart;
		}
 		else if (((fabs(vmid - Vel(start)[0])) < 
			       (Tol_pressure(start) * fabs(Vel(start)[0]))) &&
			 ((fabs(pmid - Press(start))) < 
			       (Tol_pressure(start) * Press(start))))
		{

			mside = Dens(start) * cstart;
		}
		else
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
			mside = side * (pmid - Press(start)) / 
					(vmid - Vel(start)[0]);

		shock_speed = Vel(start)[0] + (side * mside / Dens(start));
		if (mside < Tol_pressure(start)) rhomid = Dens(start);
		else	rhomid = side * mside / (shock_speed-vmid);
		temp = pmid / rhomid;
	}

	else if (s_or_r == RAREFACTION)
	{

	/* inside rarefaction */

		state_on_adiabat_with_pr(start,Press(mid),mid,TGAS_STATE);
		rhomid = Dens(mid);
		temp =  pmid / rhomid ;
	}

	else
	{
	/* inside CJ-detonation */

		shock_speed = CJ_det(CJ,TGAS_STATE,start,l_or_r);
		temp =  Press(CJ) / Dens(CJ);
		debug_print("riem_sol","temp = %g \n",temp);
	}
	t_crit = Params(start)->critical_temperature;

	if ((error_check == ERR_CHECK) && (temp < t_crit))
	{
		screen("ERROR: detonation but temperature too low \n");
		if (l_or_r == LEFT_FAMILY)
			(void) printf("LEFT_FAMILY");
		else
			(void) printf("RIGHT_FAMILY");
		if (s_or_r == SHOCK)
			(void) printf(" SHOCK\n");
		else
			(void) printf(" RAREFACTION\n");
		(void) printf("temp = %g pmid = %g vmid = %g\n",temp,pmid,vmid);
		clean_up(ERROR);
	}
	if (temp >=  t_crit) return(ABOVE_IGNIT);
	else return(BELOW_IGNIT);

}		/*end mid_temp*/

/*
*			intersection():
*
*  	This routine organizes the finding of the intersection of the
*  velocity-left curve and velocity-right curve given whether we
*  are using the non-igniting curve or the combustion curve.
*  This is passed in by l_ig_or_flow and r_ig_or_flow.
*
*  	The left and right states are input here.
*  The output are the midstate values of pressure and velocity
*  and the types of wave on left and right.
*/

LOCAL void intersection(
	double		         *pmid,
	double		         *vmid,
	RIEMANN_SOLVER_WAVE_TYPE *s_or_rl,
	RIEMANN_SOLVER_WAVE_TYPE *s_or_rr,
	Locstate	         l_state,
	Locstate	         r_state,
	int		         l_ig_or_flow,
	int		         r_ig_or_flow)
{
	int		do_newton;
	double  		pguess;
	double  		dvmdp;
	static boolean	first = YES;
	static Locstate ref_l = NULL,ref_r = NULL;

	if (first)
	{
	    first = NO;
	    (*Params(l_state)->_alloc_state)(&ref_l,Params(l_state)->sizest);
	    (*Params(r_state)->_alloc_state)(&ref_r,Params(r_state)->sizest);
	}

	debug_print("riem_sol","Entering intersection\n");


/* Get reference values for each curve:
*  that is, the CJ state for detonations, and the input
*  state for no ignition.
*/

	reference_values(ref_l,l_state,l_ig_or_flow,LEFT_FAMILY);

	reference_values(ref_r,r_state,r_ig_or_flow,RIGHT_FAMILY);


/* Get a good starting value for pressure for the newton
*  iteration. (In the case or rarefactions or CJ-detonations
*  as the only waves, the "exact" value of pmid is found.
*  No newton iteration is done in that case.)
*/

	get_p_guess(&pguess,&do_newton,l_state,ref_l,r_state,
		ref_r,l_ig_or_flow,r_ig_or_flow);

#if defined(DEBUG_GRIEMANN)
	debug_print("find_mid_state"," In intersection() pguess = %g \n",pguess);
#endif /* defined(DEBUG_GRIEMANN) */
	
	if (pguess < 0.)
	{
		(void) printf("ERROR: pguess is NEGATIVE\n");
		(void) printf("pguess = %g in intersection\n",pguess);
		clean_up(ERROR);
	}

	/* perform the newton iteration if necessary */

	if(do_newton == DO_NEWTON)
	{
		if ((l_ig_or_flow == FLOW) && (r_ig_or_flow == FLOW))
		{
			newton(pmid,vmid,(int*)s_or_rl,(int*)s_or_rr,pguess,l_state,
				r_state,curve,curve);
		}

		if ((l_ig_or_flow == FLOW) && (r_ig_or_flow == IGNITE))
		{
			newton(pmid,vmid,(int*)s_or_rl,(int*)s_or_rr,pguess,l_state,
				r_state,curve,burn_curve);
		}

		if ((l_ig_or_flow == IGNITE) && (r_ig_or_flow == FLOW))
		{
			newton(pmid,vmid,(int*)s_or_rl,(int*)s_or_rr,pguess,l_state,
				r_state,burn_curve,curve);
		}

		if ((l_ig_or_flow == IGNITE) && (r_ig_or_flow == IGNITE))
		{
			newton(pmid,vmid,(int*)s_or_rl,(int*)s_or_rr,pguess,l_state,
				r_state,burn_curve,burn_curve);
		}
	}

	/* if no newton iteration was needed find the midstate velocity  */

	else
	{
		if (l_ig_or_flow == FLOW)
		{
			curve(vmid,&dvmdp,pguess,l_state,LEFT_FAMILY);
			*s_or_rl = RAREFACTION;
		}
		else
		{
			burn_curve(vmid,&dvmdp,pguess,l_state,LEFT_FAMILY);
			*s_or_rl = CJ_DET;
		}
		if (r_ig_or_flow == FLOW) *s_or_rr = RAREFACTION;
		else *s_or_rr = CJ_DET;

		*pmid = max(Tol_pressure(l_state),pguess);
	}
#if defined(DEBUG_GRIEMANN)
	debug_print("find_mid_state","End intersection() pmid = %g vmid = %g\n ",
		*pmid, *vmid);
#endif /* defined(DEBUG_GRIEMANN) */
} 		/*end intersection*/

/*
*			newton():
*
* 	Newton finds the pressure and velocity in the middle
*  state by Newton iteration on vl(P)-vr(P)
*
*  	Vmid and pmid are found,  as well as the resulting wave on
*  each side. s_or_r = 1 for rarefaction, 2 for a shock,
*  3 for a strong detonation and 4 for a CJ_detonation.
*
*  	The left and right states are the inputs. 
*/

LOCAL void newton(
	double		*pmid,
	double		*vmid,
	int		*s_or_rl,
	int		*s_or_rr,
	double		pguess,
	Locstate	l_state,
	Locstate	r_state,
	int		(*curve_left) (double*,double*,double,Locstate,int),
	int		(*curve_right) (double*,double*,double,Locstate,int))
{
	int		count;     
	double		vmidl,vmidr;
	double		dvmdpl,dvmdpr;
	double		tol_rel;
	double		correction;

	tol_rel  = .334*(Press(l_state) + Press(r_state) +
			 pguess)*Tol_pressure(l_state);
	*pmid = max(pguess,tol_rel);
	correction = tol_rel + 1.;
	count = 0;
	while ((fabs(correction) > tol_rel) && (count < 20)) 
	{
#if defined(DEBUG_GRIEMANN)
	    /* lw and rw stand for left and right wave */
	    debug_print("find_mid_state","In newton() pmid = %g lw = %d rw = %d\n",
			           *pmid, *s_or_rl, *s_or_rr);
#endif /* defined(DEBUG_GRIEMANN) */
	    count++; 
	    *s_or_rl = (*curve_left)(&vmidl,&dvmdpl,*pmid,l_state,LEFT_FAMILY); 
	    *s_or_rr = (*curve_right)(&vmidr,&dvmdpr,*pmid,r_state,
				      RIGHT_FAMILY);

	    correction = ((vmidr - vmidl) / (dvmdpl-dvmdpr));
	    if (debugging("riem_sol"))
	    {
	        (void) printf("vmidl= %g vmidr= %g dvmdpl= %g dvmdpr= %g\n",
	    	             vmidl,vmidr,dvmdpl,dvmdpr);
		(void) printf("l_wave = %d, r_wave = %d correction = %g\n",
			      *s_or_rl,*s_or_rr,correction);
	    }
	    *pmid = *pmid + correction;
	    tol_rel = .334 * (Press(l_state) + Press(r_state) + 
			      *pmid) * Tol_pressure(l_state);
	    if (*pmid < tol_rel)
		*pmid = tol_rel;
	} 
	if (count >= 20)
	{
	    (void) printf(" There have been too many newton iterations!!! \n");
	    (void) printf(" pleft  vleft  rholeft");
	    (void) printf("       pright  vright  rhoright\n");
	    (void) printf("%g %g %g %g %g %g\n", 
	    	          Press(l_state),Vel(l_state)[0], Dens(l_state), 
	    	          Press(r_state), Vel(r_state)[0], Dens(r_state));
	    (void) printf(" latest correction = %g \n", correction);
	}
	*vmid = .5 * (vmidl+vmidr);

#if defined(DEBUG_GRIEMANN)
	if (debugging("find_mid_state"))
	{
	    (void) printf("End newton() pmid = %g vmid = %g\n",*pmid, *vmid);
	    (void) printf("s_or_rl = %d s_or_rr = %d\n",*s_or_rl, *s_or_rr);
	}
#endif /* defined(DEBUG_GRIEMANN) */

}		/*end newton*/

/*
*			reference_values():
*
* 	Given the side state and whether or not combustion is desired,
*  return a state on the v(p) curve to be used for iteration.
*
*  	For no combustion just use the original state. This divides
*  the v(p) curve into two parts: a rarefaction piece and a shock piece.
*
*  	For combustion,  the CJ-state is returned. This divides the
*  v(p) burn-curve into two pieces: the CJ-detonation followed by
*  a rarefaction piece and the strong detonation piece.
*/

LOCAL void reference_values(
	Locstate	ref,	    		/* TGas states */
	Locstate	start,
	int		ig_or_flow,
	int		l_or_r)
{
	static	size_t	sizest;
	static	boolean	first = YES;
	
	if (first)
	{
		first = NO;
		sizest = Params(start)->sizest;
	}

	if (ig_or_flow == FLOW) 	/* no combustion */
	{
		ft_assign(ref,start,sizest);	
	}

	else	   			/* combustion */
	{
		CJ_det(ref,TGAS_STATE,start,l_or_r);
		Set_other_params(ref,start);
	}
}		/*end reference_values*/

/*
*			get_p_guess():
*
* 	This routine finds a good starting pressure to iterate on to
*  find pmid.  In the rarefaction (or CJ-detonation)-rarefaction
*  (or CJ-detonation) case pmid is computed.
*
*  	Input are the left and right original states and the left and
*  right reference states and whether or not combustion is taking
*  place on each side.
*
*  	The value of pguess is output. Also output is whether or not
*  this is actually pmid or whether newton iteration must be 
*  performed next.
*  If pguess is pmid then we set do_newton to NO_NEWTON,  i.e.
*  not additional newton iterations need be performed. If
*  pguess is just a starting guess then do_newton is DO_NEWTON.
*/

LOCAL void get_p_guess(
	double		*pguess,
	int		*do_newton,
	Locstate	l_state,
	Locstate	ref_l,
	Locstate	r_state,
	Locstate	ref_r,
	int		l_ig_or_flow,
	int		r_ig_or_flow)
{
	int 		left_wave, right_wave;
	double		vr_at_pref_l, vl_at_pref_r;
	double		dvrdpl, dvldpr;
	double		p_min, eps_u, eps_p;


	/*  Find the values of vl(p) at the right reference pressure */
	/*  and of vr(p) at the left reference pressure              */

	if (l_ig_or_flow == FLOW)
		curve(&vl_at_pref_r,&dvldpr,Press(ref_r),l_state,LEFT_FAMILY);

	else
		burn_curve(&vl_at_pref_r,&dvldpr,Press(ref_r),l_state,
			LEFT_FAMILY);

	if (r_ig_or_flow == FLOW)
		curve(&vr_at_pref_l,&dvrdpl,Press(ref_l),r_state,RIGHT_FAMILY);

	else
		burn_curve(&vr_at_pref_l,&dvrdpl,Press(ref_l),r_state,
			RIGHT_FAMILY);

/* Decide what type of wave we have on each side
*  RAREFACTION means rarefaction if no combustion but CJ-detonation
*  if there is combustion on that side.
*  SHOCK means shock if there's no combustion but STRONG DETONATION
*  if there is combustion on that side.
*/

	if (vr_at_pref_l >= Vel(ref_l)[0]) left_wave = RAREFACTION;
	else left_wave = SHOCK;

	if (vl_at_pref_r <= Vel(ref_r)[0]) right_wave = RAREFACTION;
	else right_wave = SHOCK;

/* Find a good initial guess for pmid in the SHOCK (or strong 
*  detonation)-SHOCK (or strong detonation) case.    
*/

	if ((left_wave == SHOCK) && (right_wave == SHOCK))
	{
		*pguess = shock_shock_p_guess(l_state,r_state);
	}

/* Find a good initial guess for pmid in the case of a shock
*  (or strong detonation) on one side and a rarefaction (or
*  CJ-detonation) on the other side.
*/

	if (left_wave != right_wave)
	{
		*pguess = shock_raref_p_guess(l_state,l_ig_or_flow,
			r_state,r_ig_or_flow,ref_l,ref_r,
			vl_at_pref_r,vr_at_pref_l,left_wave);
	}

/* Find pmid for the rarefaction(or CJ-detonation)-rarefaction
*  (or CJ-detonation) case.
*/

	if ((left_wave == RAREFACTION) && (right_wave == RAREFACTION))
	{
	    initialize_riemann_solver(ref_l,ref_r,pguess,&p_min,EPS,
				      &eps_u,&eps_p,find_mid_state);
	}
	*do_newton = NO_NEWTON;
}		/*end get_p_guess*/



/*
*			shock_shock_p_guess():
*
*  	Input the original left and right states and return a good initial
*  guess for pmid.  The value returned is less than pmid so that
*  the newton's method guarantees convergence.
*/

LOCAL double shock_shock_p_guess(
	Locstate	l_state,
	Locstate	r_state)
{
	double guess, pguess;
	double right_const, left_const;

	/* Increase vr(p) curve and decrease vl(p) curve to get
	*  an estimate on pmid that is better for bigger
	*  velocity differences . (see notes).
	*/


	right_const = shock_shock_state_guess(r_state);
	left_const = shock_shock_state_guess(l_state);

	guess  = (Vel(l_state)[0] - Vel(r_state)[0]) /
		(left_const + right_const);
	pguess = guess * guess;

	/* Choose largest of Press(l_state), Press(r_state) and pguess */

	pguess = max(pguess,Press(l_state));
	pguess = max(pguess,Press(r_state));

	return(pguess);
}		/*end shock_shock_p_guess*/



/*
*			shock_raref_p_guess():
*
* 	This routine find a good initial guess for pmid when there is
*  a shock ( or strong detonation) on one side and a rarefaction
*  (or CJ-detonation) on the other side.
*  Input the original left and right states, the left and right
*  reference states and whether there is ignition on each side.
*
*  	A guess for pmid that is less than pmid ( to guarantee newton
*  method convergence) is returned.
*/

LOCAL double shock_raref_p_guess(
	Locstate	l_state,
	int		l_ig_or_flow,
	Locstate	r_state,
	int		r_ig_or_flow,
	Locstate	ref_l,
	Locstate	ref_r,
	double		vl_at_pref_r,
	double		vr_at_pref_l,
	int		left_wave)
{
	double		pguess;
	double		vmidl,vmidr,dvmdpr,dvmdpl;

	/* Find the intersection of the lines connecting vl(pref_l)
	*  to vl(pref_r) and vr(pref_l) to vr(pref_r) and use that
	*  value of pressure (which is greater than pmid) for one
	*  newton iteration. This new pressure is less than pmid.
	*/

	pguess = Press(ref_l) + ((Press(ref_l) - Press(ref_r)) *
		(Vel(ref_l)[0] - vr_at_pref_l) /
		(vl_at_pref_r + vr_at_pref_l - Vel(ref_l)[0] - Vel(ref_r)[0]));

	if (l_ig_or_flow == FLOW)
		curve(&vmidl,&dvmdpl,pguess,l_state,LEFT_FAMILY);
	else burn_curve(&vmidl,&dvmdpl,pguess,l_state,LEFT_FAMILY);

	if (r_ig_or_flow == FLOW)
		curve(&vmidr,&dvmdpr,pguess,r_state,RIGHT_FAMILY);
	else burn_curve(&vmidr,&dvmdpr,pguess,r_state,RIGHT_FAMILY);

	pguess = pguess + ((vmidr - vmidl) / (dvmdpl - dvmdpr));

	/* Decide which value of pressure to use left refernce pressure, */
	/* right reference pressure or the newton iterate                */

	/* Left reference pressure is better than */
	/* the iterate:				  */

	if ((left_wave == SHOCK) && (pguess < Press(ref_l))) 
		pguess = Press(ref_l);

	/* Right_wave must be a shock, and right reference pressure */
	/* is a better guess than the iterate: 			   */

	if ((left_wave == RAREFACTION) && (pguess < Press(ref_r))) 
		pguess = Press(ref_r);

	return(pguess);
}		/*end shock_raref_p_guess*/



/*	
*			 shock_shock_state_guess():
*
*	This function is used in shock_shock_p_guess
*	See this function in griemann.c
*
*/

LOCAL double shock_shock_state_guess(
	Locstate	state0)
{
	double		state_const;
	double		GAM;

	GAM = gruneisen_gamma(state0);
	state_const = sqrt(2. / (Dens(state0) * GAM));

	return state_const;
}		/*end shock_shock_state_guess*/


/*
*			curve():
*
*   	Given the side(l_or_r), an initial state "start" and pmid, 
*  curve finds vmid and dvmdp (derivative of v(p)) and whether
*  there is a shock or a rarefaction.
*  The value of curve is 1 for a rarefaction
*  and 2 for a shock. Note in this routine there is no combustion.
*
*  	To connect starting regions in which Press(l_state) = Press(r_state),
*  and Vel(l_state)[0] = Vel(r_state)[0],
*  rarefaction calculations are used.  ( In reality, the only "waves"
*  are contacts in this case). If we use shocks, later on in
*  calculations in the mid_temp routine and the wave curve routine
*  we end up with quotients 0/0.
*/

LOCAL int curve(
	double		*vmid,
	double		*dvmdp,
	double		pmid,
	Locstate	start,
	int		l_or_r)
{
	if (pmid > Press(start))
		return shock(vmid,dvmdp,pmid,start,l_or_r);
	else
		return rarefaction(vmid,dvmdp,pmid,start,l_or_r);
}		/*end curve*/

/*
*			burn_curve():
*
*  	Given the side, the initial state "start" and pmid,
*  burn_curve finds vmid and dvmdp and whether there is a strong detonation
*  or a CJ-det followed by a rarefaction.  The value of burn_curve
*  is 3 for a strong detonation and 4 for a CJ_detonation followed
*  by a rarefaction. 
*
*	To connect starting regions in which Press(l_state) = Press(r_state), 
*  and Vel(l_state)[0] = Vel(r_state)[0],
*  rarefaction calculations are used.  ( In reality, the only "waves"
*  are contacts in this case). If we use shocks, later on in
*  calculations in the mid_temp routine and the wave curve routine
*  we end up with quotients 0/0.
*/

LOCAL int burn_curve(
	double		*vmid,
	double		*dvmdp,
	double		pmid,
	Locstate	start,
	int		l_or_r)
{
	static Locstate CJ = NULL;

	if (CJ == NULL)
	{
		(*Params(start)->_alloc_state)(&CJ,Params(start)->sizest);
	}

	CJ_det(CJ,TGAS_STATE,start,l_or_r);

	if (pmid > Press(CJ))
	{
		return detonation(vmid,dvmdp,pmid,start,l_or_r);
	}
	else
	{
		Set_other_params(CJ,start);
		rarefaction(vmid,dvmdp,pmid,CJ,l_or_r);
		return CJ_DET;
	}
}		/*end burn_curve*/

/*
*			CJ_det():
*
* 	This routine finds the state behind a CJ_detonation.
*  The inputs are the initial state "start"
*  and the side (l_or_r, LEFT_FAMILY or RIGHT_FAMILY) we are on.
*/

EXPORT double CJ_det(
	Locstate	CJ,
	int		st_type_CJ,
	Locstate	start,
	int		l_or_r)
{
	double		CJ_speed;

	int avail;

	avail = (l_or_r == LEFT_FAMILY) ? l_CJ : r_CJ;

	CJ_speed = CJ_state(CJ,st_type_CJ,start,l_or_r,avail);

	if (l_or_r == LEFT_FAMILY) l_CJ = AVAILABLE;
	else r_CJ = AVAILABLE;
	return CJ_speed;
}		/*end CJ_det*/




/*
*			shock():
*
*	This routine, shock,  calculates vmid and dv/dp at a given
*  	input value of p (pressure) for a shock.
*  	Given are the initial state "start" and pmid,  as well as 
*  	the side l_or_r ( 1 for right and -1 for left).
*/

LOCAL int shock(
	double		*vmid,
	double		*dvmdp,
	double		p1,
	Locstate	start,
	int		l_or_r)
{
	double		side;
	double		pmid, pstart, dp, twom, m2;
	double		Vmid, Vstart, dV;
	double		gam, GAM;
	static Locstate mid = NULL;

	if (mid == NULL)
		(*Params(start)->_alloc_state)(&mid,Params(start)->sizest);
	
	side = (l_or_r == LEFT_FAMILY) ? -1. : 1.;
	state_w_pr_on_Hugoniot(start,p1,mid,TGAS_STATE);
	*vmid = Vel(start)[0] + side * riemann_wave_curve(start,p1);

	pstart = pressure(start);		Vstart = 1.0/Dens(start);
	pmid = pressure(mid);		Vmid = 1.0/Dens(mid);
	dp = pmid - pstart;			dV = Vmid - Vstart;
	m2 = -dp/dV;
	twom = 2.0*sqrt(m2);
	gam = adiabatic_gamma(start);
	GAM = gruneisen_gamma(start);
	*dvmdp = side*(1. + m2*(Vmid+.5*GAM*dV)/(gam*p1-.5*GAM*dp))/twom;;

	return  SHOCK;
}		/*end shock*/


/*
*			detonation():
*
*	This routine  calculates vmid and dv/dp at a given
*  input value of p (pressure) for a detonation.
*  Given are the initial state "start" and pmid,  as well as 
*  the side l_or_r ( 1 for right and -1 for left).
*
*/

LOCAL int detonation(
	double		*vmid,
	double		*dvmdp,
	double		p1,
	Locstate	start,
	int		l_or_r)
{
	double		side;
	double		pmid, pstart, dp, twom, m2;
	double		rmid, Vmid, Vstart, dV;
	double		c2mid, GAM;
	double		m;
	static Locstate mid = NULL, ubmid = NULL;

	if (mid == NULL)
	{
		(*Params(start)->_alloc_state)(&ubmid,Params(start)->sizest);
		(*Params(start)->_alloc_state)(&mid,Params(start)->sizest);
	}

	side = (l_or_r == LEFT_FAMILY) ? -1. : 1.;
	Set_other_params(mid,start);
	state_w_pr_on_Hugoniot(start,p1,mid,TGAS_STATE);
	pstart = pressure(start);
	Vstart = 1.0/Dens(start);
	pmid = pressure(mid);
	rmid = Dens(mid);
	Vmid = 1.0/rmid;
	dp = pmid - pstart;			dV = Vmid - Vstart;
	m2 = -dp/dV;
	m = sqrt(m2);
	*vmid = Vel(start)[0] + side * m * (Vstart - Vmid);

	twom = 2.0*sqrt(m2);
	c2mid = sound_speed_squared(mid);
	GAM = gruneisen_gamma(mid);
	Set_params(ubmid,start);
	state_w_pr_on_Hugoniot(start,pmid,ubmid,TGAS_STATE);
	if (Dens(ubmid) < Dens(mid))
	{
		screen("ERROR in detonation(), ");
		screen("inconsistent strong detonation\n");
		clean_up(ERROR);
	}

	*dvmdp = side*(1. + m2*(Vmid+.5*GAM*dV)/(c2mid*rmid-.5*GAM*dp))/twom;;

	return STRONG_DET;
}		/*end detonation*/

/*
*			rarefaction():
*
* 	This routine ,rarefaction, calculates vmiddle and dv/dp
*	at a given input value of p for a rarefaction.
*	Given are the initial state "start" and pmid,  as well as 
*	the side l_or_r ( 1 for right and -1 for left).
*/

LOCAL int rarefaction(
	double		*vmid, 		/* velocity in middle	*/
	double		*dvmdp,		/* dv/dp at pmid	*/
	double		pmid,		/* pressure in middle	*/
	Locstate	start,		/* initial state 	*/
	int		l_or_r)		/* left or right wave	*/
{
	static Locstate mid = NULL;
	double		side;

	if (mid == NULL)
	{
		(*Params(start)->_alloc_state)(&mid,Params(start)->sizest);
	}

	side = (l_or_r == LEFT_FAMILY) ? -1. : 1.;

	state_on_adiabat_with_pr(start,pmid,mid,TGAS_STATE);
	*vmid = Vel(start)[0] + side * riemann_wave_curve(start,pmid);
	*dvmdp = side / (Dens(mid)*sound_speed(mid));
	return  RAREFACTION;
}		/*end rarefaction*/






/*
*			vacuum():
*
*	This routine tests whether the mid_state is a vacuum 
*  The minimum velocity of any state to which the right state can
*  be connected ( vr_min) and the maximum velocity of any state
*  to which the left side may be connected ( vl_min) are
*  computed and are the output.
* 	The input is the left and right states.
* 	The value of vacuum is YES in the case of a midstate vacuum
*  and NO in the case of a middle state with positive pressure.
*/

LOCAL int vacuum(
	double		*vr_min,
	double		*vl_max,
	Locstate	l_state,
	Locstate	r_state)
{
	double		dvmdp;

	curve(vr_min,&dvmdp,Min_pressure(r_state),r_state,RIGHT_FAMILY);
	curve(vl_max ,&dvmdp,Min_pressure(l_state),l_state,LEFT_FAMILY);

	if ((*vr_min) >= (*vl_max)) return (YES) ;
	else return (NO) ;
}		/*end vacuum*/


#endif /* defined(COMBUSTION_CODE) */
