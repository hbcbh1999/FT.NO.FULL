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
*                               s2phase-eos.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <geos/s2phase.h>

        /* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL   double   S2PHASE_internal_energy(Locstate);
LOCAL   double   S2PHASE_specific_internal_energy(Locstate);
LOCAL   double   S2PHASE_pressure(Locstate);
LOCAL   double   S2PHASE_sound_speed_squared(Locstate);
LOCAL   double   S2PHASE_acoustic_impedance_squared(Locstate);

        /* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL   double   S2PHASE_entropy(Locstate);
LOCAL   double   S2PHASE_gruneisen_gamma(Locstate);
LOCAL   double   S2PHASE_temperature(Locstate);

	/* MATERIAL PROPERTY FUNCTIONS */
LOCAL	double	S2PHASE_bulk_viscosity(Locstate);
LOCAL	double	S2PHASE_shear_viscosity(Locstate);
LOCAL   double   S2PHASE_density(Locstate);

        /* VECTORIZED THERMODYNAMIC FUNCTIONS */
LOCAL   void    S2PHASE_single_eos_load_pressure(Vec_Gas*,int,int);

        /* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
        /* Purely Thermodynamic Hugoniot Functions */

        /* Velocity Related Hugoniot Functions */

        /* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL   double   S2PHASE_dens_rarefaction(double,Locstate);
LOCAL   double   S2PHASE_pressure_rarefaction(double,Locstate);

        /* General Wave Curve Functions */
LOCAL   double   S2PHASE_mass_flux(double,Locstate);
LOCAL   double   S2PHASE_mass_flux_squared(double,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
LOCAL   double   S2PHASE_oned_fan_state(double,Locstate,Locstate,Locstate,
					int,boolean*);

        /* Functions to Compute Riemann Solutions */
LOCAL   double   S2PHASE_riemann_wave_curve(Locstate,double);
LOCAL   void    S2PHASE_set_state_for_find_mid_state(Locstate,Locstate);
LOCAL   double   S2PHASE_eps_for_Godunov(Locstate,double,double);
LOCAL	void	S2PHASE_initialize_riemann_solver(Locstate,Locstate,double*,
		                                  double*,double,double*,double*,
					boolean(*)(Locstate,Locstate,
		                                double,double*,double*,double*,
						double*,double*,double*,
						RIEMANN_SOLVER_WAVE_TYPE*,
						RIEMANN_SOLVER_WAVE_TYPE*));


        /* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */

#if defined(COMBUSTION_CODE)
        /* DETONATION SPECIFIC UTILITY FUNCTIONS */
#endif /* defined(COMBUSTION_CODE) */

        /* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */

        /* INITIALIZATION UTILITY FUNCTIONS */
LOCAL   EOS*    S2PHASE_free_EOS_params(EOS*);
LOCAL   void    S2PHASE_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                              Gas_param*);
LOCAL   void    S2PHASE_fprint_EOS_params(FILE*,Gas_param*);
LOCAL   void    S2PHASE_prompt_for_EOS_params(INIT_DATA*,Gas_param*,const char*,
	                                      const char*);

LOCAL   void    set_eos_function_hooks(EOS*);

LOCAL   void    calculate_S2PHASE_coefs(S2PHASE_EOS*);
LOCAL   boolean    eq_for_gamma_l(double,double*,POINTER);

/*??*/
LOCAL   boolean    eq_for_rho(double,double*,POINTER);

LOCAL   boolean    oned_fan_aux(double,double*,POINTER);
LOCAL   double   int_dp_over_rho_c(double,Locstate,Locstate);
LOCAL   double   exact_int_dp_over_rho_c(double,Locstate,Locstate);

/* eneregy integral E = A (p-p_v) - p_vl(1/rho-1/rho_v) - (p/rho-p_v/rho_v) */
LOCAL   double   integrate_for_energy(double,S2PHASE_EOS*);

LOCAL   double   int_c_drho_over_rho(double,double,Locstate);

/* root finder for an equation gv=a_v*rho_v/p_v with respect to a_v */
LOCAL   boolean    root_for_gamma(double,double,double *,double,double,double);
LOCAL   double   pvl(double,double,double,double);
LOCAL   double   pv(double,double,double,double,double);
LOCAL   double   gv(double,double,double,double,double);
/* ??*/


EXPORT  EOS     *set_S2PHASE_eos(
        EOS     *eos)
{
        if (eos == NULL)
            scalar(&eos,sizeof(S2PHASE_EOS));
        (void) set_GENERIC_eos(eos);
        set_eos_function_hooks(eos);
        return eos;
}       /*end set_S2PHASE_eos*/


/*
*                       S2PHASE_free_EOS_params():
*
*       Frees the storage allocated for an equation of state parameter
*       function.
*/

LOCAL   EOS*    S2PHASE_free_EOS_params(
        EOS *eos)
{
        free(eos);
        return NULL;
}               /*end S2PHASE_free_EOS_params*/

/*
*                       S2PHASE_prompt_for_EOS_params():
*
*       Prompts for equation of state parameters.
*/

/*ARGSUSED*/
LOCAL   void    S2PHASE_prompt_for_EOS_params(
	INIT_DATA  *init,
        Gas_param  *params,
        const char *message1,
        const char *message2)
{
        S2PHASE_EOS *s2pheos = (S2PHASE_EOS*) params->eos;
        char s[121];
        double ssp;

	
        screen("Enter the isentropic two phase eos parameters.\n"
               "Enter sound speed for the saturated liquid: ");
        (void) Gets(s);
        (void) sscan_float(s,&ssp);
        s2pheos->a_sat_l = ssp*ssp;

        screen("Enter the adiabatic exponent gamma "
	       "(=ratio of specific heats) of the saturated vapor: ");
        (void) Gets(s);
        (void) sscan_float(s,&s2pheos->gamma_v);

        screen("Enter density of the saturated liquid: ");
        (void) Gets(s);
        (void) sscan_float(s,&s2pheos->rho_sat_l);

        screen("Enter density of the saturated vapor: ");
        (void) Gets(s);
        (void) sscan_float(s,&s2pheos->rho_sat_v);

        screen("Enter pressure of the saturated liquid: ");
        (void) Gets(s);
        (void) sscan_float(s,&s2pheos->p_sat_l);

        screen("Enter the temperature of the liquid: ");
        (void) Gets(s);
        (void) sscan_float(s,&s2pheos->t_sat_l);

        screen("Enter the specific heat at constant volume of the liquid: ");
        (void) Gets(s);
        (void) sscan_float(s,&s2pheos->cv_l);

	/* Moved to parabolic step check */
	/*
        screen("Type 'y' to have the Navier-Stokes terms "
	       "(physical viscosity)\n\tcomputed for this eos model "
	       "(y, n(dflt)): ");
        (void) Gets(s);
	if (strncmp(s,"y",1) == 0)
	*/

	{        

	    /* Moved to parabolic step check */
	    /*
 	    s2pheos->eos._compute_ns_terms = YES;
	    */

	    screen("Enter dynamic shear viscosity of the saturated liquid: ");
	    (void) Gets(s);
	    (void) sscan_float(s,&s2pheos->shear_visc_l);
	    
	    screen("Enter dynamic shear viscosity of the saturated vapor: ");
	    (void) Gets(s);
	    (void) sscan_float(s,&s2pheos->shear_visc_v);
	    
	    screen("Enter bulk viscosity of the saturated liquid: ");
	    (void) Gets(s);
	    (void) sscan_float(s,&s2pheos->bulk_visc_l);
	    
	    screen("Enter bulk viscosity of the saturated vapor: ");
	    (void) Gets(s);
	    (void) sscan_float(s,&s2pheos->bulk_visc_v);
	}

        calculate_S2PHASE_coefs(s2pheos);
}	/*end S2PHASE_prompt_for_EOS_params*/


LOCAL   void    calculate_S2PHASE_coefs(
        S2PHASE_EOS *s2pheos)
{
        double a_l = s2pheos->a_sat_l;
	double gamma_v = s2pheos->gamma_v;
        double p_l = s2pheos->p_sat_l;
        double rho_l = s2pheos->rho_sat_l;
        double rho_v = s2pheos->rho_sat_v;
        double delta_e;
        double t_l = s2pheos->t_sat_l;
        double cv_l = s2pheos->cv_l;
	
	double a_v;
	if(root_for_gamma(gamma_v,a_l,&a_v,rho_l,rho_v,p_l)
	   == FUNCTION_FAILED)
	{
	    screen("ERROR in calculate_S2PHASE_coefs(), "
		   "root_for_gamma() failed\n");
	    clean_up(ERROR);
	}
	s2pheos->a_sat_v = a_v;
	s2pheos->p_vl = pvl(a_l,a_v,rho_l,rho_v);
	s2pheos->p_sat_v = pv(a_l,a_v,rho_l,rho_v,p_l);
	s2pheos->eta_v = s2pheos->p_sat_v/pow(rho_v,gamma_v);
        s2pheos->t_sat_v = t_l;
        s2pheos->R_v = s2pheos->p_sat_v/s2pheos->t_sat_v/rho_v;
        s2pheos->S_0 = log(s2pheos->eta_v)*s2pheos->R_v/(gamma_v-1);
        s2pheos->e_sat_v = s2pheos->p_sat_v/(gamma_v-1)/rho_v;
        delta_e = integrate_for_energy(rho_l,s2pheos);
        s2pheos->delta_e = delta_e;
        s2pheos->e_sat_l = s2pheos->e_sat_v + delta_e; 
        s2pheos->gamma_l = sqrt(0.25 + a_l/(t_l*cv_l)) + 0.5;
        s2pheos->R_l = cv_l*(s2pheos->gamma_l - 1);
        s2pheos->p_inf = s2pheos->R_l*rho_l*t_l-p_l;
        s2pheos->eta_l = (s2pheos->p_inf + p_l)/pow(rho_l,s2pheos->gamma_l);
        s2pheos->e_inf = (s2pheos->gamma_l*s2pheos->p_inf+p_l)
		/(s2pheos->gamma_l-1)/rho_l-s2pheos->e_sat_l;

}
LOCAL   boolean  root_for_gamma(
	double gamma_v,
	double a_sat_l,
	double *a_v,
	double rho_sat_l,
	double rho_sat_v,
	double p_sat_l)
{
        double a_0,a_1,g_0,g_1; /* initial guess, next guess, previous value, next value */
	double step;            /* initial progression step */
	int   i,j;             /* dummy indices */
	int   tag=1;           /* 1 : no success, 0 : success */
	int   ITER1=10;        /* number of iteration for reducing step size */
	int   ITER2=10000;     /* number of iteration for finding higher gamma */
	int   ITER3=1000;       /* number of iteration for mid-point method */
	double tol=1.e-14;      /* tolerance for mid-point method */
	double err;             /* size of interval containing root */

	const char *fname = "root_for_gamma()";

	a_0=0;
	step=1;
	for(j=0;j<ITER1 && tag;j++)
	{
	  for(i=0;i<ITER2;i++)
	  {
	    double a=a_0+step;
	    double g=gv(a_sat_l,a,rho_sat_l,rho_sat_v,p_sat_l);
	    if(g<0)
	    {
	      step=step/2;
	      break;
	    }
	    if(g>gamma_v)
	    {
	      tag=0;
	      a_1=a;
	      g_1=g;
	      break;
	    }
	    a_0=a;
	    g_0=g;
	  }
	}

	if(tag)
	{
	  screen("ERROR in %s, gamma values are lower than "
		 "specified value %g in input file\n",fname,gamma_v);
	  return FUNCTION_FAILED;
	}

	for(i=0;i<ITER3;i++)
	{
	  double a_m=(a_0+a_1)/2;
	  double g_m=gv(a_sat_l,a_m,rho_sat_l,rho_sat_v,p_sat_l);

	  err=(a_1-a_0)/2;
	  if( (err < tol)  || (a_m+err == a_m) )
          {
	    *a_v=a_m;
	    return FUNCTION_SUCCEEDED;
	  }

	  if(g_m>=gamma_v)
	  {
	    a_1=a_m;
	    g_1=g_m;
	  }
	  else
	  {
	    a_0=a_m;
	    g_0=g_m;
	  }
	}

	screen("ERROR in %s, failed to find root for "
	       "specified value %g in input file within %d iteration\n",fname,gamma_v,ITER3);
	return FUNCTION_FAILED;

}/* end root_for_gamma */

/*
 *                    rho_sat_v*a_sat_v*rho_sat_l*a_sat_l*(rho_sat_v - rho_sat_l)
 *               pvl=-------------------------------------------------------------
 *                     rho_sat_v*rho_sat_v*a_sat_v - rho_sat_l*rho_sat_l*a_sat_l
 */
LOCAL double pvl(
      double a_sat_l,
      double a_sat_v,
      double rho_sat_l,
      double rho_sat_v)
{
      double TINY = 1.e-48;    /* if small number < TINY occurs, then clean up */
      double tmp=rho_sat_v*rho_sat_v*a_sat_v-rho_sat_l*rho_sat_l*a_sat_l;

      const char *fname="pvl()";

      if((tmp > 0 && tmp < TINY) || (tmp<0 && -tmp<TINY))
      {
	screen("ERROR in %s, too small number %g in denominator\n",fname,tmp);
	clean_up(ERROR);
      }

      return rho_sat_v*a_sat_v*rho_sat_l*a_sat_l*(rho_sat_v-rho_sat_l)/tmp;
}/* end pvl */

/*
 *                                      rho_sat_v*rho_sat_v*a_sat_v
 *               pv=p_sat_l + pvl * log ---------------------------
 *                                      rho_sat_l*rho_sat_l*a_sat_l
 */
LOCAL double pv(
      double a_sat_l,
      double a_sat_v,
      double rho_sat_l,
      double rho_sat_v,
      double p_sat_l)
{
      double TINY = 1.e-48;    /* if small number < TINY occurs, then clean up */
      double tmp=rho_sat_v*rho_sat_v*a_sat_v/(rho_sat_l*rho_sat_l*a_sat_l);

      const char *fname="pv()";

      if(tmp < TINY)
      {
	screen("ERROR in %s, too small number %g in argument of log\n",fname,tmp);
	clean_up(ERROR);
      }

      return p_sat_l + pvl(a_sat_l,a_sat_v,rho_sat_l,rho_sat_v)*log(tmp);
}/* end pv */

/*
 *                   a_sat_v*rho_sat_v
 *               gv=-------------------
 *                         pv
 */
LOCAL double gv(
      double a_sat_l,
      double a_sat_v,
      double rho_sat_l,
      double rho_sat_v,
      double p_sat_l)
{
      double TINY = 1.e-48;    /* if small number < TINY occurs, then clean up */
      double tmp=pv(a_sat_l,a_sat_v,rho_sat_l,rho_sat_v,p_sat_l);

      const char *fname="gv()";

      if((tmp > 0 && tmp < TINY) || (tmp < 0 && -tmp < TINY))
      {
	screen("ERROR in %s, too small number %g in denominator\n",fname,tmp);
	clean_up(ERROR);
      }

      return a_sat_v*rho_sat_v/tmp;
}/* end gv*/

LOCAL   void    set_eos_function_hooks(
        EOS *eos)
{
        /* PRIMARY THERMODYNAMIC FUNCTIONS */
        eos->_internal_energy = S2PHASE_internal_energy;
        eos->_specific_internal_energy = S2PHASE_specific_internal_energy;
        eos->_pressure = S2PHASE_pressure;
        eos->_sound_speed_squared = S2PHASE_sound_speed_squared;
        eos->_acoustic_impedance_squared = S2PHASE_acoustic_impedance_squared;

        /* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
        eos->_entropy = S2PHASE_entropy;
	eos->_gruneisen_gamma = S2PHASE_gruneisen_gamma;
        eos->_temperature = S2PHASE_temperature;

	/* MATERIAL PROPERTY FUNCTIONS */
	eos->_bulk_viscosity = S2PHASE_bulk_viscosity;
	eos->_shear_viscosity = S2PHASE_shear_viscosity;
        eos->_density = S2PHASE_density;

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
        eos->_single_eos_load_pressure =
            S2PHASE_single_eos_load_pressure;

	/* Purely Thermodynamic Hugoniot Functions */
        eos->_dens_Hugoniot = S2PHASE_dens_rarefaction;

        /* Purely Thermodynamic Adiabatic Wave Curve Functions */
        eos->_dens_rarefaction = S2PHASE_dens_rarefaction;
        eos->_pressure_rarefaction = S2PHASE_pressure_rarefaction;

        /* General Wave Curve Functions */
        eos->_mass_flux = S2PHASE_mass_flux;
        eos->_mass_flux_squared = S2PHASE_mass_flux_squared;
	eos->_oned_fan_state = S2PHASE_oned_fan_state;
       
        /* Functions for the Evaluation of Riemann Solutions */

        /* Functions to Compute Riemann Solutions */
        eos->_riemann_wave_curve = S2PHASE_riemann_wave_curve;
        eos->_set_state_for_find_mid_state =
	    S2PHASE_set_state_for_find_mid_state;
        eos->_eps_for_Godunov = S2PHASE_eps_for_Godunov;
	eos->_initialize_riemann_solver = S2PHASE_initialize_riemann_solver;

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */

        eos->_prompt_for_EOS_params = S2PHASE_prompt_for_EOS_params;
        eos->_fprint_EOS_params = S2PHASE_fprint_EOS_params;
        eos->_read_print_EOS_params = S2PHASE_read_print_EOS_params;
        eos->_free_EOS_params = S2PHASE_free_EOS_params;

}


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*                       S2PHASE_internal_energy():
*
*       Returns the internal energy per unit volume of a state.
*/

LOCAL   double   S2PHASE_internal_energy(
        Locstate state)
{
        double   rho_sat_v = Rho_sat_v(state);
        double   rho_sat_l = Rho_sat_l(state);
        double   e_sat_v = E_sat_v(state);
        double   rho = Dens(state);
        double   alpha = (rho - rho_sat_l)/(rho_sat_v - rho_sat_l);
        
        if (rho < 0)
        {
            screen("ERROR in S2PHASE_internal_energy(), ");
            screen("density is negative : rho = %g\n",rho);
	    screen("state :\n");
	    screen("        specific entropy = %g\n",    entropy(state));
	    fprint_state_type(stdout,"State type = ",state_type(state));
	    screen("Params state = %lu\n",gas_param_number(Params(state)));
            clean_up(ERROR);
        }

        if (alpha >= 1) /* vapor phase */
            return Eta_v(state)*pow(rho,Gamma_v(state))/(Gamma_v(state)-1);

        if (alpha < 1 && alpha > 0) /* mixed phase */ 
            return (e_sat_v + integrate_for_energy(rho,
				S2PHASE_Eos(state)))*rho;

        if (alpha <= 0)  /* liquid phase */
            return Eta_l(state)*pow(rho,Gamma_l(state))/(Gamma_l(state)-1) 
			+ P_inf(state) - E_inf(state)*rho;

}               /*end S2PHASE_internal_energy*/


/*
*                       S2PHASE_specific_internal_energy():
*
*       Returns the internal energy per unit volume of a state.
*/
LOCAL   double   S2PHASE_specific_internal_energy(
        Locstate state)
{
        double   rho_sat_v = Rho_sat_v(state);
        double   rho_sat_l = Rho_sat_l(state);
        double   e_sat_v = E_sat_v(state);
        double   rho = Dens(state);
        double   alpha = (rho - rho_sat_l)/(rho_sat_v - rho_sat_l);

        if (rho < 0)
        {
       	    screen("ERROR in S2PHASE_specific_internal_energy(), ");
            screen("density is negative : rho = %g\n",rho);
	    screen("state :\n");
	    screen("        specific entropy = %g\n",    entropy(state));
	    fprint_state_type(stdout,"State type = ",state_type(state));
	    screen("Params state = %lu\n",gas_param_number(Params(state)));
       	    clean_up(ERROR); 
        }

        if (alpha >= 1.0) /* vapor phase */
          return Eta_v(state)*pow(rho,Gamma_v(state)-1.0)/(Gamma_v(state)-1.0);

        else if (alpha < 1.0 && alpha > 0.0) /* mixed phase */ 
          return e_sat_v + integrate_for_energy(rho,S2PHASE_Eos(state));
        
        else /* alpha <= 0  liquid phase */
          return Eta_l(state)*pow(rho,Gamma_l(state)-1.0)/(Gamma_l(state)-1.0) +
                 P_inf(state)/rho - E_inf(state);
}               /*end S2PHASE_specific_internal_energy*/


/*
*                       S2PHASE_pressure():
*
*       Returns the thermodynamic pressure of a state.
*
*                                    dE  |
*                            P = -  ---- |
*                                    dV  |S
*
*       Where E = specific internal energy,  V = specific volume,  and
*       S = specific entropy.
*/

LOCAL   double   S2PHASE_pressure(
        Locstate state)
{
	double p_l,p_vl,a_l,a_v,rho_l,rho_v,rho,alpha;

	if (is_obstacle_state(state))
            return HUGE_VAL;

        p_l = P_sat_l(state);
        p_vl = P_vl(state);
        a_l = A_sat_l(state);
        a_v = A_sat_v(state);
        rho_l = Rho_sat_l(state);
        rho_v = Rho_sat_v(state);
        rho = Dens(state);
        alpha = (rho - rho_l)/(rho_v - rho_l);

        if (rho < 0.0)
        {
            screen("ERROR in S2PHASE_pressure(), ");
            screen("density is negative : rho = %g\n",rho);
	    screen("state :\n");
	    screen("        specific entropy = %g\n",    entropy(state));
	    fprint_state_type(stdout,"State type = ",state_type(state));
	    screen("Params state = %lu\n",gas_param_number(Params(state)));
	    clean_up(ERROR); 
        }
 
        if (alpha >= 1.0) /* vapor phase */
            return Eta_v(state)*pow(rho,Gamma_v(state));

	else if (alpha < 1.0 && alpha > 0.0) /* mixed phase */ 
            return p_l + p_vl*log((rho_v*a_v*(rho_l+alpha*(rho_v-rho_l)))/
                        (rho_l*(rho_v*a_v-alpha*(rho_v*a_v-rho_l*a_l))));
        
	else  /* liquid phase */
            return Eta_l(state)*pow(rho,Gamma_l(state)) - P_inf(state);
}               /*end S2PHASE_pressure*/


/*
*                       S2PHASE_sound_speed_squared():
*
*       Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*                       c = ---- |
*                           drho |S
*/

LOCAL   double   S2PHASE_sound_speed_squared(
        Locstate state)
{
        double a_l = A_sat_l(state);
        double a_v = A_sat_v(state);
        double rho_l = Rho_sat_l(state);
        double rho_v = Rho_sat_v(state);
        double rho = Dens(state);
        double alpha = (rho - rho_l)/(rho_v - rho_l);
 
        if (rho < 0)
        {
            screen("ERROR in S2PHASE_sound_speed_squared(), ");
            screen("density is negative : rho = %g\n",rho);
	    screen("state :\n");
	    screen("        specific entropy = %g\n",    entropy(state));
	    fprint_state_type(stdout,"State type = ",state_type(state));
	    screen("Params state = %lu\n",gas_param_number(Params(state)));
            clean_up(ERROR); 
          }

        if (alpha >= 1) /* vapor phase */
          return Gamma_v(state)*Eta_v(state)*pow(rho,Gamma_v(state)-1);

        else if (alpha < 1 && alpha > 0) /* mixed phase */ 
          return 1/((alpha*rho_v+(1-alpha)*rho_l)*
                    (alpha/(rho_v*a_v)+(1-alpha)/(rho_l*a_l)));
        
        else /* alpha <= 0  liquid phase */
          return Gamma_l(state)*Eta_l(state)*pow(rho,Gamma_l(state)-1);

}               /*end S2PHASE_sound_speed_squared*/


/*
*               S2PHASE_acoustic_impedance_squared():
*
*       Returns the square of the local acoustic impedance of the state.
*
*                        2     dP  |
*                       i = - ---- |
*                              dV  |S
*/

LOCAL   double   S2PHASE_acoustic_impedance_squared(
        Locstate state)
{
  return Dens(state)*Dens(state)*sound_speed_squared(state);

}               /*end S2PHASE_acoustic_impedance_squared*/



/*********** SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS*********/
/*
*                       S2PHASE_entropy():
*/

LOCAL   double   S2PHASE_entropy(
        Locstate state)
{
	return S_0(state);
} /* end S2PHASE_entropy() */


/*
*                       S2PHASE_temperature():
*/

LOCAL   double   S2PHASE_temperature(
        Locstate state)
{
        double alpha = (Dens(state) - Rho_sat_l(state))/
                      (Rho_sat_v(state) - Rho_sat_l(state));
 
        if (Dens(state) < 0)
          {
            printf("ERROR in S2PHASE_temperature(), ");
            screen("density is negative : rho = %g\n",Dens(state));
	    screen("state :\n");
	    screen("        specific entropy = %g\n",    entropy(state));
	    fprint_state_type(stdout,"State type = ",state_type(state));
	    screen("Params state = %lu\n",gas_param_number(Params(state)));
	    clean_up(ERROR); 
          }

        if (alpha >= 1) /* vapor phase */
          return Eta_v(state)*pow(Dens(state), Gamma_v(state)-1)/R_v(state);

        if (alpha < 1 && alpha > 0) /* mixed phase */ /*wrong but does not
affect computations*/ 
          return T_sat_l(state) + alpha*(T_sat_v(state) - T_sat_l(state));
        
        if (alpha <= 0)  /* liquid phase */
          return Eta_l(state)*pow(Dens(state),Gamma_l(state)-1)/R_l(state);


} /* end S2PHASE_temperarure() */

/*
*			S2PHASE_gruneisen_gamma():
*
*	Returns the dimensionless Gruneisen exponent
*
*
*                                                 dP/dE |
*		GAMMA = - d(log T)/d(log V) |  =  -----  V
*                                            S     rho
*
*	As usual P = thermodynamic pressure,  V = specific volume
*	rho = density, E = specific internal energy,
*	and  S = specific entropy.
*
*
*/
/*ARGSUSED*/
LOCAL	double	S2PHASE_gruneisen_gamma(
	Locstate state)
{
	return 0.0;
}		/*end S2PHASE_gruneisen_gamma*/

/***************END SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS *******/

/************** MATERIAL PROPERTY FUNCTIONS ********************************/

LOCAL	double	S2PHASE_bulk_viscosity(
	Locstate state)
{
        double rho_l = Rho_sat_l(state);
        double rho_v = Rho_sat_v(state);
        double rho = Dens(state);
        double alpha = (rho - rho_l)/(rho_v - rho_l);
        double mu_l = Bulk_visc_l(state);
	double mu_v = Bulk_visc_v(state);

        if (alpha >= 1) /* mostly vapor */
          return mu_v;

	else if (alpha > 0) /* mixed phase */ 
          return alpha*mu_v+(1-alpha)*mu_l; /* linear from mu to zero */
        
	else  /* liquid phase */
          return mu_l;
}	/*end S2PHASE_bulk_viscosity */

LOCAL	double	S2PHASE_shear_viscosity(
	Locstate state)
{
        double rho_l = Rho_sat_l(state);
        double rho_v = Rho_sat_v(state);
        double rho = Dens(state);
        double alpha = (rho - rho_l)/(rho_v - rho_l);
        double mu_l = Shear_visc_l(state);
	double mu_v = Shear_visc_v(state);

        if (rho < 0)
        {
            printf("ERROR in S2PHASE_shear_viscosity(), ");
            screen("density is negative : rho = %g\n",rho);
	    screen("state :\n");
	    screen("        specific entropy = %g\n",    entropy(state));
	    fprint_state_type(stdout,"State type = ",state_type(state));
	    screen("Params state = %lu\n",gas_param_number(Params(state)));
	    clean_up(ERROR); 
	    return ERROR_FLOAT;
	}
	
	if(debugging("schmidt_viscosity"))
	  /* geometric linear from mu_l,
	     Kubota Model (1989) */
	  mu_v= mu_l*rho/rho_l;

        if (alpha > 1) /* mostly vapor */
          return mu_v;

	else if (alpha > 0) /* mixed phase */ 
	{
	    double x = alpha * rho_v/rho;

	    if(debugging("schmidt_viscosity"))
	    /* geometric linear from mu_l,
	       Kubota Model (1989) */
	      return mu_l*rho/rho_l;
	    if(debugging("dukler_viscosity"))
	    /* linear from mu_l to mu_v, 
	       Dukler Model (1964), 
	       laminar flow */
	      return alpha*mu_v+(1-alpha)*mu_l;

	    /* weighted linear from mu_l to mu_v, 
	       Cicchitti Model (1960), 
	       laminar flow */
	    return x*mu_v+(1-x)*mu_l;
	}
	else  /* liquid phase */
          return mu_l;
}	/*end S2PHASE_shear_viscosity */


/*
*               S2PHASE_density():
*
*   Returns density after specification of pressure for TGAS_STATE
*
*      gamma_l    P + P_inf
*   rho       = ----------- 	: liquid phase
*                   eta_l
*
*      gamma_v      P 
*   rho       = ------- 	: vapor phase
*                 eta_v
*
*                       rho_l ( rho_v^2 a_v - rho_l^2 a_l)
*   rho= ------------------------------------------------------------------ 
*        rho_l(rho_v a_v-rho_l a_l)+exp((P_l-P)/P_vl)rho_v a_v(rho_v-rho_l)
*				: mixed phase
*/
LOCAL   double   S2PHASE_density(
        Locstate state)
{
        double rho;
	switch(state_type(state))
	{
	case    EGAS_STATE:
	case	FGAS_STATE:
	case	GAS_STATE:
	case	VGAS_STATE:
	    rho=Dens(state);
	    break;
	case	TGAS_STATE:
	    if(Press(state) >= P_sat_l(state))
		rho=pow((Press(state)+P_inf(state))/Eta_l(state),
				1/Gamma_l(state));
	    else if(Press(state) <= P_sat_v(state))
	        rho=pow(Press(state)/Eta_v(state),1/Gamma_v(state));
	    else
	        {
		  double r_l=Rho_sat_l(state);
		  double r_v=Rho_sat_v(state);
		  double a_l=A_sat_l(state);
		  double a_v=A_sat_v(state);
		  double p_l=P_sat_l(state);
		  double p_vl=P_vl(state);
		  double ratio=exp( -(Press(state)-p_l)/p_vl );
		  rho=r_l*(r_v*r_v*a_v-r_l*r_l*a_l)/(r_l*(r_v*a_v-
		  		r_l*a_l)+ratio*r_v*a_v*(r_v-r_l));
		}
	    break;
	default:
	    screen("ERROR in S2PHASE_density(), no such state type\n");
	    clean_up(ERROR);
	}
	return rho;
}               /*end S2PHASE_density*/
/************** END MATERIAL PROPERTY FUNCTIONS ****************************/

/***************VECTORIZED THERMODYNAMIC FUNCTIONS *************************/
/***************END VECTORIZED THERMODYNAMIC FUNCTIONS *********************/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/

/***************Purely Thermodynamic Hugoniot Functions*********************/
/***************End Purely Thermodynamic Hugoniot Functions*****************/


/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/
/*      
*                       S2PHASE_dens_rarefaction():
*
*       Given the state st0 and the pressure on the other side of
*       a simple wave in steady irrotational flow, this
*       function returns the density on the other side.
*
*       The answer is give by the solution of the ordinary differential
*       equation
*
*               dh/dP = V,  h(p0) = h0;
*
*       where h is the specific enthalpy,  and the derivatives are taken
*       at constant entropy.
*/

/* actually, the density found is not related at all to st0.  
   We use st0 just to get the s2phase constants p_sat_l, p_sat_v, etc. 
   The density and pressure are directly related, so we are simply 
   finding the density corresponding to the given pressure */

LOCAL   double   S2PHASE_dens_rarefaction(
        double    p1,
        Locstate st0)
{
        double p_l = P_sat_l(st0);
        double p_v = P_sat_v(st0);
        double rho_l = Rho_sat_l(st0);
        double rho_v = Rho_sat_v(st0);
	double p_vl, a_l, a_v;
        double rho1;

	if (rho_v <= rho_l)
	{
            if (p1 <= p_v) /* vapor phase */
                return pow(p1/Eta_v(st0),1.0/Gamma_v(st0));
	    else  if (p_l <= p1) /*liquid phase*/
                return pow((p1+P_inf(st0))/Eta_l(st0),1.0/Gamma_l(st0));
	    else /* mixed phase */ 
            {
                p_vl = P_vl(st0);
                a_l = A_sat_l(st0);
                a_v = A_sat_v(st0);
		rho1 = rho_l*(rho_v*rho_v*a_v - rho_l*rho_l*a_l)/
		       (rho_v*a_v*(rho_v-rho_l)*exp(-(p1 - p_l)/p_vl) -
			 rho_l*(rho_l*a_l - rho_v*a_v));
                return rho1;
            }
	}
	else
	{
	    screen("\nERROR in S2PHASE_dens_rarefaction: rho_v = %g > rho_l = %g\n"
		   ,rho_v,rho_l);
	    clean_up(ERROR);
	    return ERROR_FLOAT;
	}

}               /*end S2PHASE_dens_rarefaction*/

LOCAL   boolean eq_for_rho(
        double   rho,
        double   *function_for_p,
        POINTER state)
{
        double p_l = P_sat_l(state);
        double p_vl = P_vl(state);
        double a_l = A_sat_l(state);
        double a_v = A_sat_v(state);
        double rho_l = Rho_sat_l(state);
        double rho_v = Rho_sat_v(state);
        double alpha = (rho - rho_l)/(rho_v - rho_l);

        *function_for_p = p_l +
		p_vl*log((rho_v*a_v*(rho_l+alpha*(rho_v-rho_l)))/ 
		(rho_l*(rho_v*a_v-alpha*(rho_v*a_v-rho_l*a_l))));

        return FUNCTION_SUCCEEDED;
}
/* end eq_for_rho() */
/*      
*                       S2PHASE_pressure_rarefaction():
*
*       Given the state st0 and the density on the other side of
*       a simple wave in steady irrotational flow, this
*       function returns the pressure on the other side.
*
*       The answer is give by the solution of the ordinary differential
*       equation
*
*               de/dV = -P,  e(V0) = e0;
*
*       where e is the specific internal energy,  and the derivatives are taken
*       at constant entropy.
*/

LOCAL   double   S2PHASE_pressure_rarefaction(
        double    rho1,
        Locstate st0)
{
        static Locstate st1;
        static size_t   sizest;
        sizest = Params(st0)->sizest;

	if (st1 == NULL)
            (*Params(st0)->_alloc_state)(&st1,sizest);
        set_type_of_state(st1,GAS_STATE);
        Set_params(st1,st0);
        Dens(st1) = rho1;
        return pressure(st1);
}               /*end S2PHASE_pressure_rarefaction*/

/*
*               S2PHASE_set_state_for_find_mid_state():
*
*       Copies the Gas state st into the thermodynamic
*       version Tst, for some EOSs a VGas state is set.
*
*       Technical function added for enhanced performance.
*/

LOCAL   void    S2PHASE_set_state_for_find_mid_state(
        Locstate Tst,
        Locstate st)
{
        set_state(Tst,TGAS_STATE,st);
}               /*end S2PHASE_set_state_for_find_mid_state*/

/*
*                       S2PHASE_eps_for_Godunov():
*
*       Returns a tolerance to be used to determine convergence of the
*       of the Riemann solver.
*
*       Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL   double   S2PHASE_eps_for_Godunov(
        Locstate state,
        double pstar,
        double r_eps)
{
        return r_eps;
}               /*end S2PHASE_eps_for_Godunov*/

/*
*			S2PHASE_initialize_riemann_solver()
*
*	Computes the epsilons and the initial guess for the pressure
*	in the secant iteration of find_mid_state.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	void	S2PHASE_initialize_riemann_solver(
	Locstate Tsl,
	Locstate Tsr,
	double *pstar,
	double *p_min,
	double eps,
	double *eps_u,
	double *eps_p,
	boolean (*fd_md_st)(Locstate,Locstate,double,double*,double*,
			double*,double*,double*,double*,
			RIEMANN_SOLVER_WAVE_TYPE*,RIEMANN_SOLVER_WAVE_TYPE*))
{
	double	pl = pressure(Tsl),		pr = pressure(Tsr);

	if ((pl < 0.0) || (pr < 0.0))
	{
	    screen("ERROR in S2PHASE_initialize_riemann_solver(),\n "
		   "before computation pl = %lf, pr = %lf ",pl,pr);
	    clean_up(ERROR);
	}

	*pstar = 0.5*(pl + pr);
	*p_min = max(Min_pressure(Tsr),Min_pressure(Tsl));
	*pstar = max(*pstar,*p_min);
	*eps_u = *eps_p = eps;
}		/*end S2PHASE_initialize_riemann_solver*/


/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/

/*
*                       S2PHASE_riemann_wave_curve():
*
*       Evalutes the forward wave family wave curve defined by
*
*                _
*               |
*               |
*               |                                1/2
*               |   [ (Pstar  -  P0) * ( V0 - V) ]     if Pstar > P0
*               |
*               |
*               / 
*              /
*              \
*               \               
*               |
*               |        / Pstar     |
*               |       /            |
*               |       \      dP    |
*               |        \   ------  |                 if Pstar < P0
*               |         \   rho c  |
*               |         /          |
*               |        / P0        | S
*               |_
*
*/

LOCAL   double   S2PHASE_riemann_wave_curve(
        Locstate st0,
        double    pstar)
{
	double p0 = pressure(st0);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
        if (pstar < Min_pressure(st0))
            pstar = Min_pressure(st0);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

        if (pstar >= p0)
        {
	    double V0 = 1.0/Dens(st0);
	    double Vstar = 1.0/dens_Hugoniot(pstar,st0);
	    if (Vstar > V0)
	    	return 0.0;
	    else
	    	return sqrt((pstar-p0)*(V0 - Vstar));
	}
	else
	    return int_dp_over_rho_c(pstar,st0,NULL);
}               /*end S2PHASE_riemann_wave_curve*/

LOCAL   double  dRF(double,POINTER); /* what is this line for? */

LOCAL   double    int_c_drho_over_rho(
        double    rho0,
        double    rho1,
        Locstate state)
{
	Spline            *spl;
        double             Iv, I0, I1, Il;
        double             rho_v = Rho_sat_v(state);
        double             rho_l = Rho_sat_l(state);
        double             p_v, p_l;
	double             alpha0, alpha1;
	double             gam_v, Gam_v, gam_l, Gam_l, p_inf;
	int               ier;
	int               one;

	alpha0 = (rho0 - rho_l)/(rho_v - rho_l);
	alpha1 = (rho1 - rho_l)/(rho_v - rho_l);

	if (alpha0 >= 1.0 && alpha1 >= 1.0)/*both vapor*/
	{
            p_v = P_sat_v(state);
	    gam_v = Gamma_v(state);
	    Gam_v = gam_v - 1.0;
	    return (2.0/Gam_v)*sqrt(gam_v*p_v/rho_v)*
		    (pow(rho1/rho_v,0.5*Gam_v) - pow(rho0/rho_v,0.5*Gam_v));
	}
	else if (alpha0 <= 0.0 && alpha1 <= 0.0)/*both liquid*/
	{
            p_l = P_sat_l(state);
	    gam_l = Gamma_l(state);
	    Gam_l = gam_l - 1.0;
            p_inf = P_inf(state);
	    return (2.0/Gam_l)*sqrt(gam_l*(p_l+p_inf)/rho_l)*
		    (pow(rho1/rho_l,0.5*Gam_l) - pow(rho0/rho_l,0.5*Gam_l));
	}
	if (Int_c_drho_over_rho_mix_spline(state) == NULL) /*Unset*/
	{
            double  a_l = A_sat_l(state);
            double  a_v = A_sat_v(state);
	    double  *w, *wrk;
	    double  *x, *y, dx;
	    double  *t, *c;
	    double  rtol, atol;
	    double  s, fp;
	    double  xb, xe;
	    double  abserr;
	    int    lwrk, *iwrk;
	    int    i, m, n, k, iopt, nest;
	    int    neval;
	    QUADRATURE_STATUS ierq;
            Locstate st1;

	    scalar(&spl,sizeof(Spline));
	    Int_c_drho_over_rho_mix_spline(state) = spl;
	    spl->m = m = 1001; /*TOLERANCE*/
	    iopt = 0;/* Smoothing spline */
	    uni_array(&spl->x,m,FLOAT);
	    uni_array(&spl->y,m,FLOAT);
	    spl->k = k = 3; /* Cubic spline approximates T_ref */
	    spl->nest = nest = m + k + 1;
	    uni_array(&spl->t,nest,FLOAT);
	    uni_array(&spl->c,nest,FLOAT);

	    /* Work space storage */
	    uni_array(&w,m,FLOAT);
	    for (i = 0; i < m; ++i)
	        w[i] = 1.0;
	    lwrk = m*(k+1)+nest*(7+3*k);
	    uni_array(&wrk,lwrk,FLOAT);
	    uni_array(&iwrk,nest,FLOAT);

	    x = spl->x;
	    y = spl->y;
	    t = spl->t;
	    c = spl->c;

	    dx = fabs(rho_l - rho_v)/(m - 1);
	    x[0] = min(rho_l,rho_v);
	    for (i = 1; i < m; ++i)
		x[i] = x[0] + i*dx;
	    x[m-1] = max(rho_l,rho_v);
	    xb = x[0]; xe = x[m-1];

	    for (i = 0; i < k; ++i)
	        t[i] = x[0];
	    for (i = k; i < (nest-k); ++i)
	        t[i] = x[i-k];
	    for (i = nest-k; i < nest; ++i)
	        t[i] = x[m-1];

            (*Params(state)->_alloc_state)(&st1,Params(state)->sizest);
            set_type_of_state(st1,GAS_STATE);
            Set_params(st1,state);

	    spl->rtol = rtol = 1.0e-6;/*TOLERANCE*/
	    spl->atol = 0.5*(a_l+a_v)*rtol;
	    y[0] = 0.0;
	    for (i = 1; i < m; ++i)
	    {
		double I;
	        I = SimpRule(dRF,(POINTER)st1,x[i-1],x[i],atol,rtol,
	                     &abserr,&neval,&ierq);
	        switch (ierq)
	        {
	        case INVALID_EPSILON:
	            screen("ERROR in int_c_drho_over_rho(), invalid epsilons\n"
		           "atol = %"FFMT", rtol = %"FFMT"\n",atol,rtol);
	            clean_up(ERROR);
	            break;
	        case INACCURATE_INTEGRAL:
	            if ((fabs(abserr) > 20.0*atol) &&
		        (fabs(abserr) > 20.0*rtol*I))
	            {
		    /*
		     * Don't print a warning if we are close to satifying the
		     * error requirements
		     */
	                (void) printf("WARNING in int_c_drho_over_rho(), "
			              "inaccurate result\n \tneval = %d, "
			              "abserr = %"FFMT" result = %"FFMT"\n",
			              neval,abserr,I);
	            }
	            break;
	        case ACCURATE_INTEGRAL:
	        default:
	            break;
	        }
		y[i] = y[i-1] + I;
	    }

	    s = 0.0; /* s = 0 implies interpolation */
	    FORTRAN_NAME(curfit)(&iopt,&m,x,y,w,&xb,&xe,&k,&s,&nest,&n,t,c,&fp,
		                 wrk,&lwrk,iwrk,&ier);
	    switch (ier)
	    {
	    case 0:
	        break;
	    case -1:
	        break;
	    case -2:
	        break;
	    case 1:
	        screen("ERROR in int_c_drho_over_rho() spline fit failed\n"
	               "error flag = %d, not enough storage\n",ier);
	        clean_up(ERROR);
	        break;
	    case 2:
	        screen("ERROR in int_c_drho_over_rho() spline fit failed\n"
	               "error flag = %d, s too small?\n",ier);
	        clean_up(ERROR);
	        break;
	    case 3:
	        screen("ERROR in int_c_drho_over_rho() spline fit failed\n"
	               "error flag = %d, max interations exceeded "
		       "s too small?\n",ier);
	        clean_up(ERROR);
	        break;
	    case 10:
	        screen("ERROR in int_c_drho_over_rho() spline fit failed\n"
	               "error flag = %d, invalid input data\n");
	        clean_up(ERROR);
	        break;
	    default:
	        screen("ERROR in int_c_drho_over_rho() spline fit failed\n"
	               "error flag = %d, unknown error\n",ier);
	        clean_up(ERROR);
	        break;
	    }
	    spl->n = n;
	    free_these(4,w,wrk,iwrk,st1);
	}
	Iv = Il = 0.0;
	if (alpha0 >= 1.0) /*State 0 is vapor*/
	{
            p_v = P_sat_v(state);
	    gam_v = Gamma_v(state);
	    Gam_v = gam_v - 1.0;
	    Iv += (2.0/Gam_v)*sqrt(gam_v*p_v/rho_v)*
		    (1.0 - pow(rho0/rho_v,0.5*Gam_v));
	    rho0 = rho_v;
	}
	else if (alpha0 <= 0.0) /*State 0 is liquid*/
	{
            p_l = P_sat_l(state);
	    gam_l = Gamma_l(state);
	    Gam_l = gam_l - 1.0;
            p_inf = P_inf(state);
	    Il += (2.0/Gam_l)*sqrt(gam_l*(p_l+p_inf)/rho_l)*
		    (1.0 - pow(rho0/rho_l,0.5*Gam_l));
	    rho0 = rho_l;
	}
	if (alpha1 >= 1.0) /*State 1 is vapor*/
	{
            p_v = P_sat_v(state);
	    gam_v = Gamma_v(state);
	    Gam_v = gam_v - 1.0;
	    Iv += (2.0/Gam_v)*sqrt(gam_v*p_v/rho_v)*
		    (pow(rho1/rho_v,0.5*Gam_v) - 1.0);
	    rho1 = rho_v;
	}
	else if (alpha1 <= 0.0)
	{
            p_l = P_sat_l(state);
	    gam_l = Gamma_l(state);
	    Gam_l = gam_l - 1.0;
            p_inf = P_inf(state);
	    Il += (2.0/Gam_l)*sqrt(gam_l*(p_l+p_inf)/rho_l)*
		    (pow(rho1/rho_l,0.5*Gam_l) - 1.0);
	    rho1 = rho_l;
	}
	spl = Int_c_drho_over_rho_mix_spline(state);
	one = 1;
	FORTRAN_NAME(splev)(spl->t,&spl->n,spl->c,&spl->k,&rho1,&I1,&one,&ier);
	FORTRAN_NAME(splev)(spl->t,&spl->n,spl->c,&spl->k,&rho0,&I0,&one,&ier);

        return Iv + I1 - I0 + Il;
}		/*end int_c_drho_over_rho */


/*ARGSUSED*/
LOCAL	double dRF(
	double	rho,
	POINTER prms)
{
    	Locstate state = (Locstate) prms;
        double    p_vl = P_vl(state);
        double    a_l = A_sat_l(state);
        double    a_v = A_sat_v(state);
        double    rho_l = Rho_sat_l(state);
        double    rho_v = Rho_sat_v(state);
	double    c2;

	if (rho <= 0.0)
	    return 0.0;
	c2 = (p_vl/rho)*(rho_v*rho_v*a_v - rho_l*rho_l*a_l)/
	               ((rho_v-rho)*rho_v*a_v+(rho-rho_l)*rho_l*a_l);

	return (c2 > 0.0) ? sqrt(c2)/rho : 0.0;
}		/*end dRF*/


/*
*                       S2PHASE_mass_flux_squared():
*
*       Returns square of the mass flux across a wave.
*
*                                2
*                2   | (P - P0) |
*               m  = | -------  |
*                    | (U - U0) |
*
*       Where 
*               P0 = pressure ahead of the shock
*               U0 = velocity ahead of the shock
*               P = pressure behind the shock
*               U = velocity behind the shock
*
*/

LOCAL   double   S2PHASE_mass_flux_squared(
        double p,
        Locstate st0)
{
        double p0 = pressure(st0);

        if (fabs(p - p0) < p0*EPSILON)
           return acoustic_impedance_squared(st0);

	{
            double du = riemann_wave_curve(st0,p);
            
            return sqr((p-p0)/du);
	}

}               /*end S2PHASE_mass_flux_squared*/


LOCAL   double   S2PHASE_mass_flux(
        double p,
        Locstate st0)
{
        return sqrt(S2PHASE_mass_flux_squared(p,st0));
}               /*end S2PHASE_mass_flux*/


/*
*                       S2PHASE_fprint_EOS_params():
*
*       Prints the parameters that define the given equation of state.
*       NOTE:  This is not really an initialization function,  but it is
*       convenient to locate it next the the corresponding read function.
*/

LOCAL   void    S2PHASE_fprint_EOS_params(
        FILE *file,
        Gas_param *params)
{
        S2PHASE_EOS *s2phase = (S2PHASE_EOS *)params->eos;

        (void) fprintf(file,"\tEquation of state = %d ISENTROPIC_TWO_PHASE\n",
                       ISENTROPIC_TWO_PHASE);
        fprint_float(file,"\ta_sat_l = ",s2phase->a_sat_l,"\n");
        fprint_float(file,"\ta_sat_v  = ",s2phase->a_sat_v,"\n");
        fprint_float(file,"\tp_sat_l  = ",s2phase->p_sat_l,"\n");
        fprint_float(file,"\tp_inf  = ",s2phase->p_inf,"\n");
        fprint_float(file,"\tp_sat_v  = ",s2phase->p_sat_v,"\n");
        fprint_float(file,"\tp_vl  = ",s2phase->p_vl,"\n");
        fprint_float(file,"\trho_sat_l  = ",s2phase->rho_sat_l,"\n");
        fprint_float(file,"\trho_sat_v  = ",s2phase->rho_sat_v,"\n");
        fprint_float(file,"\tcv_l  = ",s2phase->cv_l,"\n");
        fprint_float(file,"\te_sat_l  = ",s2phase->e_sat_l,"\n");
        fprint_float(file,"\te_sat_v  = ",s2phase->e_sat_v,"\n");
        fprint_float(file,"\tdelta_e  = ",s2phase->delta_e,"\n");
        fprint_float(file,"\te_inf  = ",s2phase->e_inf,"\n");
        fprint_float(file,"\tS_0  = ",s2phase->S_0,"\n");
        fprint_float(file,"\tR_l  = ",s2phase->R_l,"\n");
        fprint_float(file,"\tR_v  = ",s2phase->R_v,"\n");
        fprint_float(file,"\tt_sat_l  = ",s2phase->t_sat_l,"\n");
        fprint_float(file,"\tgamma_l  = ",s2phase->gamma_l,"\n");
        fprint_float(file,"\tgamma_v  = ",s2phase->gamma_v,"\n");
        fprint_float(file,"\teta_l  = ",s2phase->eta_l,"\n");
        fprint_float(file,"\teta_v  = ",s2phase->eta_v,"\n");
	fprint_float(file,"\tliquid shear viscosity  = ",
		     s2phase->shear_visc_l,"\n");
	fprint_float(file,"\tvapor shear viscosity  = ",
		     s2phase->shear_visc_v,"\n");
	fprint_float(file,"\tliquid bulk viscosity  = ",
		     s2phase->bulk_visc_l,"\n");
	fprint_float(file,"\tvapor bulk viscosity  = ",
		     s2phase->bulk_visc_v,"\n");
}               /*end S2PHASE_fprint_EOS_params */

/*
*                       S2PHASE_read_print_EOS_params():
*
*       Reads the equation of state data as printed by 
*	S2PHASE_fprint_EOS_params. This is restart function.
*/

/*ARGSUSED*/
LOCAL   void    S2PHASE_read_print_EOS_params(
        INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	FILE        *file = io_type->file;
        S2PHASE_EOS *s2phase = (S2PHASE_EOS *)params->eos;

	s2phase->int_c_drho_over_rho_mix_spline = NULL;
	s2phase->a_sat_l = fread_float("a_sat_l = ",io_type);
        s2phase->a_sat_v = fread_float("a_sat_v  = ",io_type);
        s2phase->p_sat_l = fread_float("p_sat_l  = ",io_type);
        s2phase->p_inf = fread_float("p_inf  = ",io_type);
        s2phase->p_sat_v = fread_float("p_sat_v  = ",io_type);
        s2phase->p_vl = fread_float("p_vl  = ",io_type);
        s2phase->rho_sat_l = fread_float("rho_sat_l  = ",io_type);
        s2phase->rho_sat_v = fread_float("rho_sat_v  = ",io_type);
        s2phase->cv_l = fread_float("cv_l  = ",io_type);
        s2phase->e_sat_l = fread_float("e_sat_l  = ",io_type);
        s2phase->e_sat_v = fread_float("e_sat_v  = ",io_type);
        s2phase->delta_e = fread_float("delta_e  = ",io_type);
        s2phase->e_inf = fread_float("e_inf  = ",io_type);
        s2phase->S_0 = fread_float("S_0  = ",io_type);
        s2phase->R_l = fread_float("R_l  = ",io_type);
        s2phase->R_v = fread_float("R_v  = ",io_type);
        s2phase->t_sat_l = fread_float("t_sat_l  = ",io_type);
        s2phase->gamma_l = fread_float("gamma_l  = ",io_type);
        s2phase->gamma_v = fread_float("gamma_v  = ",io_type);
        s2phase->eta_l = fread_float("eta_l  = ",io_type);
        s2phase->eta_v = fread_float("eta_v  = ",io_type);
	if (fgetstring(file,"d_visc_l  = ")) /*Old style*/
	{
            s2phase->shear_visc_l = fread_float(NULL,io_type);
            s2phase->shear_visc_v = fread_float("d_visc_v  = ",io_type);
            s2phase->bulk_visc_l = fread_float("b_visc_l  = ",io_type);
            s2phase->bulk_visc_v = fread_float("b_visc_v  = ",io_type);
	}
	else if (fgetstring(file,"liquid shear viscosity  = "))
	{
            s2phase->shear_visc_l = fread_float(NULL,io_type);
            s2phase->shear_visc_v =
		fread_float("vapor shear viscosity  = ",io_type);
            s2phase->bulk_visc_l =
		fread_float("liquid bulk viscosity  = ",io_type);
            s2phase->bulk_visc_v =
		fread_float("vapor bulk viscosity  = ",io_type);
	}
}		/*end S2PHASE_read_print_EOS_params */

/* integration along isentrope from rho_v to rho *
 * eneregy integral E = A (p-p_v) - p_vl(1/rho-1/rho_v) - (p/rho-p_v/rho_v)
 * where A = (rho_v*a_v - rho_l*a_l)/(rho_v^2*a_v - rho_l^2*a_l)
 */
LOCAL   double integrate_for_energy(
        double       rho,
        S2PHASE_EOS *s2eos)
{
        double rho_v = s2eos->rho_sat_v;
        double rho_l = s2eos->rho_sat_l;
        double p_v  = s2eos->p_sat_v;
        double p_l  = s2eos->p_sat_l;
        double p_vl = s2eos->p_vl;
        double a_v  = s2eos->a_sat_v;
        double a_l  = s2eos->a_sat_l;
	double A,p;
        
	A=rho_v*rho_v*a_v - rho_l*rho_l*a_l;
	if( A == 0.0 || rho < rho_v || rho > rho_l)
	{
	    screen("\nError in integrate_for_energy(): A=%g, rho=%g, rho_v=%g, rho_l=%g\n",A, rho,rho_v,rho_l);
	    clean_up(ERROR);
	    return ERROR_FLOAT;
	}
	A=(rho_v*a_v - rho_l*a_l)/A;
	p=p_l+p_vl*log(rho*rho_v*a_v*(rho_v-rho_l)/rho_l
		       /((rho_v*rho_v*a_v-rho_l*rho_l*a_l) - rho*(rho_v*a_v-rho_l*a_l)));
	return A*(p-p_v) - p_vl*(1/rho-1/rho_v) - (p/rho-p_v/rho_v);

} /* end integrate_for_energy() */

LOCAL   double old_integrate_for_energy(
        double       rho,
        S2PHASE_EOS *s2eos)
{
        double rho_v = s2eos->rho_sat_v;
        double rho_l = s2eos->rho_sat_l;
        double p_v  = s2eos->p_sat_v;
        double p_l  = s2eos->p_sat_l;
        double p_vl = s2eos->p_vl;
        double a_v  = s2eos->a_sat_v;
        double a_l  = s2eos->a_sat_l;
        double e, x;
	double A, B, dp, V_v, V;
        
	A = (rho_l*a_l - rho_v*a_v)*rho_l/(rho_v*a_v*(rho_v-rho_l));
	B = (rho_v*rho_v*a_v - rho_l*rho_l*a_l)/(rho_v*a_v*(rho_v-rho_l));

	x = rho/(A*rho + B*rho_l);
	dp = (p_v-p_l)/p_vl;
	V_v = 1.0/rho_v;
	V = 1.0/rho;
        e = p_l*(V_v-V) + (p_vl/(B*rho_l))*(exp(-dp)*(1.0-dp)-(1.0-log(x))/x);
        
        return e;
}


/*
*			S2PHASE_oned_fan_state():
*
*	This is a utility function provided for the evaluation of states
*	in a simple wave.   Given sta, it solves for stm using the
*	equation:
*
*	                 / p_m        |            / p_m        |
*	                /             |           /             |
*	                \       dP    |           \       G dP  |
*	w = c_m - c_a +  \    -----   |         =  \     ------ |
*	                  \   rho c   |             \    rho c  |
*	                  /           |             /           |
*	                 /p_a         | S = S_a    / p_a        | S = S_a
*
*                                            / c_m        |
*                                           /             |
*                                           \        dc   |
*                                         =  \     ------ |
*                                             \     mu^2  |
*                                             /           |
*                                            / c_a        | S = S_a
*
*
*	here c is the sound speed,  rho the density,  S the specific entropy,
*	p the pressure,  and mu^2 = (G - 1)/G,  where G is the fundamental
*	derivative of gas dynamics.  The returned state1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	state1 can be used to evaluate the pressure, density, and sound speed
*	of the state inside the rarefaction fan.
*	
*	Input data:
*		w = value of w as defined above
*		sta = state ahead of fan
*		stb = state behind fan
*
*	Output data:
*		stm = state inside fan
*		vacuum = 1 if stm is a vacuum,  0 otherwise
*
*	Returns the sound speed of the answer state stm.
*/

struct FAN_AUX	{
	Locstate sta;
	Locstate stm;
	double	ca, cm;
};

LOCAL	double	S2PHASE_oned_fan_state(
	double    w,
	Locstate sta,
	Locstate stb,
	Locstate stm,
	int      stype_m,
	boolean     *vacuum)
{
	double	pa = pressure(sta);
	double	pb = pressure(stb);
	double	pm, ca;
	double   wb;
	struct FAN_AUX	Fprms;

	*vacuum = NO;
	Fprms.sta = sta;
	Fprms.stm = stm;
	Fprms.ca = ca = sound_speed(sta);
	Set_params(stm,sta);

	if (!oned_fan_aux(pb,&wb,(POINTER)&Fprms))
	{
	    screen("ERROR in S2PHASE_oned_fan_state(), can't evalute "
		   "state behind rarefaction fan\n");
	    screen("pa = %g, pb = %g\n",pa,pb);
	    screen("wb = %g\n",wb);
	    clean_up(ERROR);
	    return Fprms.cm;
	}
	if (((w <= wb) && (wb <= 0.0)) || ((0.0 <= wb) && (wb <= w)))
	{
	    set_state(stm,stype_m,stb);
            Fprms.cm = sound_speed(stm); /* New Code */
	}
	else if (((wb <= 0.0) && (0.0 <= w)) || ((w <= 0.0) && (0.0 <= wb)))
	{
	    set_state(stm,stype_m,sta);
            Fprms.cm = sound_speed(stm); /* New Code */
	}
	else if (find_root(oned_fan_aux,(POINTER)&Fprms,w,&pm,pa,pb,
		      EPSILON*ca,EPSILON*(pb-pa))) 
	{
	    set_state(stm,stype_m,stm);

	    if(debugging("check_find_root"))
	    {
		screen("find_root finds a root : pm = %g\n",pm);
		screen("\tpressure(stm)=%g\n",pressure(stm));
		screen("\tpa=%g, pb=%g\n",pa,pb); 
	    }
	}
	else
	{
	    *vacuum = YES;
	    screen("ERROR in S2PHASE_oned_fan_state(), can't find root\n");
	    clean_up(ERROR);
	    return Fprms.cm;
	}
	return Fprms.cm;
}		/* end S2PHASE_oned_fan_state*/

LOCAL	boolean	oned_fan_aux(
	double	p,
	double	*w,
	POINTER	prms)
{
	struct FAN_AUX	*fprms = (struct FAN_AUX  *)prms;
	Locstate sta = fprms->sta;
	Locstate stm = fprms->stm;
	double    cm;

	*w = int_dp_over_rho_c(p,sta,stm);
	fprms->cm = cm = sound_speed(stm);

	*w += cm - fprms->ca;
	return FUNCTION_SUCCEEDED;
}		/*end oned_fan_aux*/

/*
*
*			int_dp_over_rho_c():
*
*
*	Returns the value of the integral:
*
*
*            / p1        |
*           /            |
*           \      dP    |            2
*            \   ------  |    =  ----------- ( c1 - c0 )  for vapor phase
*             \   rho c  |       gamma_v - 1
*             /          |
*            / p0        | S
*
*
*        =  2 SQRT(P_vl)( SQRT(V0 - V) - SQRT(V1 - V) ) for mixed phase
*
*                                 2
*                         =  ----------- ( c1 - c0 )  for liquid phase
*                            gamma_l - 1
*
*             rho_l a_l - rho_v a_v           1           1
*       V = -------------------------   V0 = ----   V1 = ----
*           rho_l^2 a_l - rho_v^2 a_v        rho0        rho1
*
*	Also returns the state on the adiabat containing state0 with
*	pressure P.
*
*/
LOCAL	double	int_dp_over_rho_c(
	double		p1,
	Locstate	state0,
	Locstate	state1)
{
	double	ans;
        double   p0 = pressure(state0);
	double   dp;

	DEBUG_ENTER(int_dp_over_rho_c)
	
	if (p1 == -HUGE_VAL) /*Reduce to finite interval*/
	    p1 = pressure_rarefaction(0.001*Dens(state0),state0);/*TOLERANCE*/

	if (state1 == NULL)
	{
	    static Locstate tmpst = NULL;

	    if (tmpst == NULL)
	    	(*Params(state0)->_alloc_state)(&tmpst,Params(state0)->sizest);
	    state1 = tmpst;
	}
	set_state(state1,EGAS_STATE,state0);

	if (p0 != 0.0) 
	    dp = (p1 - p0)/p0;
	else if (p1 != 0.0) 
	    dp = (p1 - p0)/p1;
	else 
	{
	    set_state(state1,EGAS_STATE,state0); /* Need this ? */
	    return 0.0;   
	}

	if(p0 == p1)                             /* Need this ? */
	{
	    set_state(state1,EGAS_STATE,state0); /* Need this ? */
	    return 0.0;
	}

	ans = exact_int_dp_over_rho_c(p1,state0,state1); 

	state_on_adiabat_with_pr(state1,p1,state1,EGAS_STATE); /* Set to P1 */

	DEBUG_LEAVE(int_dp_over_rho_c)
	return ans;
}		/*end int_dp_over_rho_c*/

LOCAL   double   exact_int_dp_over_rho_c(
	double           p1,
	Locstate        state0,
	Locstate	state1)
{
	double	ans;
	double   rho0 = Dens(state0);
	double   rho1;
	double   c0 = sound_speed(state0);
	double   c1;
        double   p0 = pressure(state0);
	double	sgn=1;
        double   p_vl = P_vl(state0);
        double   a_l = A_sat_l(state0);
        double   a_v = A_sat_v(state0);
        double   rho_l = Rho_sat_l(state0);
        double   rho_v = Rho_sat_v(state0);
	double   g_v = Gamma_v(state0);
	double   g_l = Gamma_l(state0);
	double	V = (rho_l*a_l - rho_v*a_v)/(rho_l*rho_l*a_l - rho_v*rho_v*a_v);

	state_on_adiabat_with_pr(state1,p1,state1,EGAS_STATE); /* Set to P1 */
	rho1 = Dens(state1);
	c1 = sound_speed(state1);

	if(rho1 < rho0)
	{
	    sgn = rho0;
	    rho0 = rho1;
	    rho1 = sgn;
	    sgn = c0;
	    c0 = c1;
	    c1 = sgn;
	    sgn = -1;
	}

	if(rho1 <= rho_v)         /* Only vapor */
	{
	    ans=2./(g_v-1.)*(c1-c0);
	}
	else if(rho0 >= rho_l)    /* Only liquid */
	{
	    ans=2./(g_l-1.)*(c1-c0);
	}
	else if(rho0 < rho_v)     /* Contains vapor */
	{
	    ans = 2./(g_v-1.)*(sqrt(a_v)-c0);      /* vapor part */
	    if(rho1 > rho_l)                       /* Contains liquid */
	    {
	        ans+= 2.*sqrt(p_vl)*(sqrt(1./rho_v - V)-sqrt(1./rho_l - V));   
						   /* mixed part */
		ans+= 2./(g_l-1.)*(c1-sqrt(a_l));                              
						   /* liquid part */
	    }
	    else
	        ans+= 2.*sqrt(p_vl)*(sqrt(1./rho_v - V)-sqrt(1./rho1 - V));    
						   /* mixed part */
	}
	else                      /* No vapor */
	{
	    ans=0;
	    if(rho1 > rho_l)      /* Contains liquid */
	    {
	        ans+= 2./(g_l-1.)*(c1-sqrt(a_l));   /* liquid part */
		rho1=rho_l;
	    }
	    ans+= 2.*sqrt(p_vl)*(sqrt(1./rho0 - V)-sqrt(1./rho1 - V));         
	    					    /* mixed part */
	}
	return ans*sgn;
}	/* end exact_int_dp_over_rho_c */

LOCAL	void	S2PHASE_single_eos_load_pressure(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	static Locstate tmpst;
	Locstate *state = vst->state + offset;
	double *rho = vst->rho + offset;
	double *p = vst->p + offset;
	int k, j, dim;

	dim = Params(state[0])->dim;
	if (tmpst == NULL)
	{
	    (*Params(state[0])->_alloc_state)(&tmpst,Params(state[0])->sizest);
	    zero_scalar(tmpst,Params(state[0])->sizest);
	    set_type_of_state(tmpst,GAS_STATE);
	}
	for (k = 0; k < vsize; k++)
	{
	    Dens(tmpst) = rho[k];
	    Set_params(tmpst,state[k]);
	    p[k] = pressure(tmpst);
	}
}		/*end S2PHASE_single_eos_load_pressure*/
