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
*                               girm_linear.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Richtmyer linearized theroy for shock-contact interaction.
*
*	The following are some suggestions that might help in running this
*	solver.
*
*	1.  Use units in time that are O(1), in particular it is desirable
*	    that dt be O(1).
*
*	2.  It seems that the code works best with a rather small CFL factor.
*           such as 0.1.
*
*	3. Use a fairly fine mesh for example 300 zones for the self-similar
*	   region seems to be suitable for many cases.
*
*	4.  The code seems to have a tendency to blow-up if too large an
*           initial amplitude is selected.  Of course this is out of the
*	    linear regime anyway.
*
*	5.  Choose the endpoint values of the region where log(time)
*	    is used in the solution to be as small a possible.
*
*/
 
#if defined(TWOD)
#include <ginit/ginit.h>


typedef	struct {
	double	ksq, a0m0;
	double	a0p0, a1p0, a2p0, a3p0;
	double  *c, *rho, *g;   /* for rarefaction wave */
	double	w1, w2, w3;
	double	v0, v1, v2, v3, u0, u3, c1sq, c2sq, K1, K2, U0, deno1, deno2;
	Locstate	sta, st0, st1, st2, st3;
}	PHY_PARAM;

typedef	struct {
	int	n1, n2, n3;
	double	dx1, dx2, dxr;
	double	dels, s_start, s_end, t_end, smooth;
	double	dt_prt_vel, dt_prt_prf;
	int     tri_prt; /* Print trigrid? */
	int     nx, ny; /* Grid zones in each direction for trigrid printing */
	double	lx, ux, ly, uy, y0; /* Corners and contact posn for trigrid
				     *printing */
	double   a0m0phys; /* Physical perturbation for scaling of
			   * pressure perturbations */
	double   t_offset; /* Time to start printing.  Used to synchronize
		           * linear and nonlinear printouts */
}	NUM_PARAM;


	 /* LOCAL Function Declarations */
LOCAL	NUM_PARAM*	setup_num_param(INIT_DATA*,INIT_PHYSICS*,LAYER_SYS*,
					PHY_PARAM*);
LOCAL	PHY_PARAM*	setup_phy_param(LAYER_SYS*,Locstate,Locstate,Locstate,
					double);
LOCAL	double	f_A0(double,PHY_PARAM*);
LOCAL	double	find_dp(double,double,double*,PHY_PARAM*,NUM_PARAM*);
LOCAL	double	p21_at_w2(PHY_PARAM*,NUM_PARAM*);
LOCAL	double	refracted_circulation_per_unit_len_per_sin_inc_ang(Locstate,
								   Locstate,
								   Locstate,
								   Locstate,
								   Locstate);
LOCAL	void	empirical_growth_rates(PHY_PARAM*);
LOCAL	void	info_inside_rarefaction_wave(PHY_PARAM*,NUM_PARAM*);
LOCAL	void	init_conditions(PHY_PARAM*,NUM_PARAM*);
LOCAL	void	march_forward(INTERFACE*,PHY_PARAM*,NUM_PARAM*);
LOCAL	void	rm_linear_print_interface(INTERFACE*,double,double,double,
					  double,double,double,NUM_PARAM*,
					  PHY_PARAM*,int);
LOCAL	void	rm_linear_print_tri_soln(INTERFACE*,double,double,
					 PHY_PARAM*,NUM_PARAM*);
LOCAL	void	small_t_solution(PHY_PARAM*,NUM_PARAM*);
LOCAL	void	solve_1d_rp(LAYER_SYS*,Locstate,Locstate,Locstate,double*);
LOCAL	void	update(double,double,int,PHY_PARAM*,NUM_PARAM*);
LOCAL	void	update_inside_rr(double,double,int,PHY_PARAM*,NUM_PARAM*);


	/* LOCAL global varibles */
LOCAL	int	n1, n2, n3, rflct_wv_type;
LOCAL	double  *pm = NULL, *p = NULL, *pp = NULL, *qm = NULL, *q = NULL,
               *qp = NULL, *Am = NULL, *A = NULL, *Ap = NULL, *Bm = NULL,
	       *B = NULL, *Bp = NULL, *Vm = NULL, *V = NULL, *Vp = NULL,
		a1m, a1, a1p, a2m, a2, a2p, a3;

EXPORT	void	rm_linear(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	LAYER_SYS	*layer_sys)
{
	size_t		sizest = layer_sys->front->sizest;
	INTERFACE	*intfc = layer_sys->front->interf;
	double		pstar;
	Locstate	sta, st3, st0;
	NUM_PARAM       *num_param;
	PHY_PARAM       *phy_param;
	int		dim = intfc->dim;

	if (dim != 2)
	{
	    screen("ERROR in rm_linear(), dim = %d != 2 not supported\n",dim);
	    clean_up(ERROR);
	}

	alloc_state(intfc,&sta,sizest);
	alloc_state(intfc,&st3,sizest);
	alloc_state(intfc,&st0,sizest);

	solve_1d_rp(layer_sys,sta,st0,st3,&pstar);

	phy_param = setup_phy_param(layer_sys,sta,st0,st3,pstar);

	empirical_growth_rates(phy_param);

	num_param = setup_num_param(init,ip,layer_sys,phy_param);

	if ((rflct_wv_type == RAREFACTION) && (n3 > 0))
	    info_inside_rarefaction_wave(phy_param,num_param);

	init_conditions(phy_param,num_param);

	march_forward(intfc,phy_param,num_param);

	clean_up(0);
}		/*end rm_linear*/


LOCAL	void	solve_1d_rp(
	LAYER_SYS	*layer_sys,
	Locstate	sta,
	Locstate	st0,
	Locstate	st3,
	double		*ppstar)
{
	int		         dim = layer_sys->front->rect_grid->dim;
	size_t		         sizest = layer_sys->front->sizest;
	int	  	         num_layers = layer_sys->num_layers;
	int                      i, j;
	RIEMANN_SOLVER_WAVE_TYPE l_wv, r_wv;
	double		         ustarl, ustarr;
	double		         pstarl, pstarr;
	double		         ml, mr;
	LAYER		         **layer = layer_sys->layer;
	COMP_TYPE	         *ct;

	if (layer[num_layers]->lower_surf->wv_type != BACKWARD_SHOCK_WAVE)
	{
	    screen("ERROR in solve_1d_rp(), upper most intfc should be a "
		   "backward shock wave!\n");
	    clean_up(ERROR);
	}
	
	/* assuming the left state is always ahead */
	i = num_layers - 1;
	ct = comp_type(layer[i]->comp);
	set_state(sta,TGAS_STATE,ct->extra);
	Vel(sta)[0] = Vel(sta)[dim-1]; 
	Vel(sta)[dim-1] = 0;

	if (is_scalar_wave(layer[i]->lower_surf->wv_type) == 0)
	{
	    screen("ERROR in solve_1d_rp(), the first intfc is not contact!\n");
	    clean_up(ERROR);
	}

	/* find left (ahead) state */
	j = i-1;
	ct = comp_type(layer[j]->comp);
	if (ct->type != AMBIENT)
	{
	    screen("ERROR in solve_1d_rp(), Ahead state is not ambient!\n");
	    clean_up(ERROR);
	}
	set_state(st0,TGAS_STATE,ct->extra);
	Vel(st0)[0] = Vel(st0)[dim-1];
	Vel(st0)[dim-1] = 0;

	/* find right (behind) state */
	j = i+1;
	ct = comp_type(layer[j]->comp);
	if (ct->type != AMBIENT)
	{
	    screen("ERROR in solve_1d_rp(), Behind state is not ambient!\n");
	    clean_up(ERROR);
	}
	ft_assign(st3,ct->extra,sizest);
	set_state(st3,TGAS_STATE,st3);
	Vel(st3)[0] = Vel(st3)[dim-1];
	Vel(st3)[dim-1] = 0;

	/* solving the Riemann problem */
	(void) find_mid_state(st0,st3,0.0,&pstarl,&pstarr,&ustarl,&ustarr,
			      &ml,&mr,&l_wv,&r_wv);

	*ppstar = 0.5*(pstarl + pstarr);

	if (l_wv != SHOCK)
	{
	    screen("ERROR in sovle_1d_rp(), "
	           "The transmitted wave is not a shock!\n");
	    clean_up(ERROR);
	}

	rflct_wv_type = r_wv;
	(void) printf("\nThe reflected wave is a %s wave.\n",
		      rflct_wv_type == SHOCK ? "shock" : "rarefaction");

	Vel(st0)[0] -= ustarl;
	Vel(st3)[0] -= ustarr;
}		/*end solve_1d_rp*/


LOCAL	PHY_PARAM	*setup_phy_param(
	LAYER_SYS	*layer_sys,
	Locstate	sta,
	Locstate	st0,
	Locstate	st3,
	double		pstar)
{ 	
	Front	  *front = layer_sys->front;
	INTERFACE *intfc = front->interf;
	size_t	  sizest = front->sizest;
	double	  nor[MAXD];	/* have to be double! */
	double	  pa, p0, p1, p2, p3, c0, c1, c2, c3;
	double	  massflux, strength, mach, temp;
	double	  w1, w2, w3, i12, i22, g_gamma;
	double	  v0, v1, v2, v3, u0, c1sq, c2sq, K1, K2, U0, a0m0;
	Locstate  st1, st2;
	PHY_PARAM *phy_param;

	alloc_state(intfc,&st1,sizest);
	alloc_state(intfc,&st2,sizest); 	
	scalar(&phy_param,sizeof(PHY_PARAM));

	phy_param->ksq = 1.0;
	phy_param->a0m0 = a0m0 = layer_sys->layer[1]->upper_surf->fpoly->A[0];
	if (a0m0 == 0.0)
	    phy_param->a0m0 = a0m0 = 1.0;

	/* set up states in regions 1 & 2 */
	nor[0] = -1.0;
	(void) s_polar_4(BEHIND_PRESSURE,pstar,&w1,nor,st0,st1,TGAS_STATE);
	if (rflct_wv_type == SHOCK)
	{
	    nor[0] = 1.0;
	    (void) s_polar_4(BEHIND_PRESSURE,pstar,&w2,nor,st3,st2,TGAS_STATE);
	    w3 = w2;
	}
	else if (rflct_wv_type == RAREFACTION)
	{
	    w3 = Vel(st3)[0] + sound_speed(st3);
	    state_on_adiabat_with_pr(st3,pstar,st2,TGAS_STATE);
	    w2 = sound_speed(st2);
	    if (w2 > w3)
	    {
	        screen("ERROR in setup_phy_param(), w2 > w3 !\n");
	        clean_up(ERROR);
	    }
	}
	else
	{
	    screen("ERROR in setup_phy_param(), "
	           "Reflected wave is not a shock or a rarefaction!\n");
	    clean_up(ERROR);
	}

	pa = Press(sta); 
	p0 = Press(st0); 	p1 = Press(st1); 
	p2 = Press(st2); 	p3 = Press(st3);

	massflux =  mass_flux((double)p3,sta);
	(void) printf("\nmass flux of the incident shock = %g\n",massflux);

	phy_param->U0 = U0 = massflux/Dens(sta);

	mach = U0/sound_speed(sta);
	(void) printf("\nMach number of the incident shock = %g\n",mach);

	strength =  (p3-pa)/p3;
	(void) printf("\nIncident shock strength = %g\n",(p3-pa)/p3);
	(void) printf("Transmitted shock strength = %g\n",(p1-p0)/p1);
	strength =  (p2-p3)/p2;
	if (rflct_wv_type == RAREFACTION)
	    (void) printf("Reflected rarefaction strength = %g\n",strength);
	else
	    (void) printf("Reflected shock strength = %g\n",strength);

	phy_param->w1 = w1;
	phy_param->w2 = w2;
	phy_param->w3 = w3;
	if (rflct_wv_type == RAREFACTION)
	    (void) printf("\nw1 = %g,   w2 = %g,   w3 = %g\n",w1,w2,w3);
	else
	    (void) printf("\nw1 = %g,   w2 = %g\n",w1,w2);

	phy_param->sta = sta;
	phy_param->st0 = st0;
	phy_param->st1 = st1;
	phy_param->st2 = st2;
	phy_param->st3 = st3;

	c0 = sound_speed(st0); c1 = sound_speed(st1);
	c2 = sound_speed(st2); c3 = sound_speed(st3);
	phy_param->c1sq = c1sq = c1*c1;
	phy_param->c2sq = c2sq = c2*c2;
	(void) printf("c0 = %g,   c1 = %g,   c2 = %g,   c3 = %g\n",c0,c1,c2,c3);
	(void) printf("p0 = %g,   p1 = %g,   p2 = %g,   p3 = %g\n",
		      pressure(st0),pressure(st1),pressure(st2),pressure(st3));
	(void) printf("rho0 = %g,   rho1 = %g,   rho2 = %g,   rho3 = %g\n",
		      Dens(st0),Dens(st1),Dens(st2),Dens(st3));
	(void) printf("s0 = %g,   s1 = %g,   s2 = %g,   s3 = %g\n",
		      entropy(st0),entropy(st1),entropy(st2),entropy(st3));
	(void) printf("u0 = %g,   u1 = %g,   u2 = %g,   u3 = %g\n",
		      Vel(st0)[0],Vel(st1)[0],Vel(st2)[0],Vel(st3)[0]);

	temp = 1.0/w1 + 1.0/c1 + 1.0/w2 + 1.0/c2;
	(void) printf("\noptimal layer thickness---L1 = %g*L, L2 = %g*L\n",
		      (1.0/w2+1.0/c2)/temp,(1.0/w1 + 1.0/c1)/temp);
	(void) printf("approximate signal contamination time---%g*L\n\n",
		      (1.0/w2+1.0/c2)*(1.0/w1 + 1.0/c1)/temp);

	phy_param->v0 = v0 = 1.0/Dens(st0);	
	phy_param->v1 = v1 = 1.0/Dens(st1);
	phy_param->v2 = v2 = 1.0/Dens(st2);	
	phy_param->v3 = v3 = 1.0/Dens(st3);
	phy_param->u0 = u0 = Vel(st0)[0];		
	phy_param->u3 = Vel(st3)[0];
	i12 = acoustic_impedance_squared(st1);
	i22 = acoustic_impedance_squared(st2);
	g_gamma = gruneisen_gamma(st1);
	/*
	*  Ki = rho^2*c^2*dV/dP|h, see Menikoff-Plohr 
	*  ``The Riemann Problem for Fluid Flow of Real Materials'',
	*  Rev. Mod. Phys., V 61, 1989 pp. 75--130
	*  Note that our Ki here is acutally 1/Ki, using the notion
	*  of Ki in Richtmyer's original paper.
	*/
	phy_param->K1 = K1 = (1 + 0.5*g_gamma*(v1 - v0)/v1)*i12/
		             (i12 - 0.5*g_gamma*(p1 - p0)*Dens(st1));
	g_gamma = gruneisen_gamma(st2);
	phy_param->K2 = K2 = (1 + 0.5*g_gamma*(v2 - v3)/v2)*i22/
		             (i22 - 0.5*g_gamma*(p2 - p3)*Dens(st2));

	phy_param->deno1 = w1 + 0.5*c1sq/w1 + 0.5*w1*K1;
	phy_param->deno2 = w2 + 0.5*c2sq/w2 + 0.5*w2*K2;

	phy_param->a0p0 = (1.0 - u0/U0)*a0m0;
	phy_param->a1p0 = (1.0 - (u0 + w1)/U0)*a0m0;
	phy_param->a2p0 = (1.0 - (u0 - w2)/U0)*a0m0;
	if (rflct_wv_type == RAREFACTION)
	    phy_param->a3p0 = (1.0 - (u0 - w3)/U0)*a0m0;
	else
	    phy_param->a3p0 = 0.0;
	(void) printf("\n");
	(void) printf("a0(+0)/a0(-0) = %g\n",phy_param->a0p0/a0m0);
	(void) printf("a1(+0)/a0(-0) = %g\n",phy_param->a1p0/a0m0);
	(void) printf("a2(+0)/a0(-0) = %g\n",phy_param->a2p0/a0m0);
	if (rflct_wv_type == RAREFACTION)
	    (void) printf("a3(+0)/a0(-0) = %g\n",phy_param->a3p0/a0m0);

	return	phy_param;
}		/*end setup_phy_param*/


LOCAL	void	empirical_growth_rates(PHY_PARAM *phy_param)
{
	double		u0, k, a0, rho1, rho2, tmp;
	double		g1, g2, a1, a2, s_t2, s_r2;
	double		a0m0 = phy_param->a0m0;
	Locstate	sta = phy_param->sta;
	Locstate	st0 = phy_param->st0;
	Locstate	st1 = phy_param->st1;
	Locstate	st2 = phy_param->st2;
	Locstate	st3 = phy_param->st3;

	u0 = Vel(st0)[0];		
	k = sqrt(phy_param->ksq);

	(void) printf("\n--------------------------------------------------\n");
	(void) printf("Various empirical formulas for the growth rate:\n");
	(void) printf("\nI. Richtmyer's impulsive model:\n");

	rho1 = Dens(st0);		
	rho2 = Dens(sta);
	a0 = 1.0;
	(void) printf("\npre-rho, pre-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);
	a0 = phy_param->a0p0/a0m0;
	(void) printf("\npre-rho, post-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);
	a0 = 0.5*(a0m0 + phy_param->a0p0)/a0m0;
	(void) printf("\npre-rho, avg-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);

	rho1 = 0.5*(Dens(st0)+Dens(st1));		
	rho2 = 0.5*(Dens(sta)+Dens(st2));
	a0 = 1.0;
	(void) printf("\navg-rho, pre-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);
	a0 = phy_param->a0p0/a0m0;
	(void) printf("\navg-rho, post-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);
	a0 = 0.5*(a0m0 + phy_param->a0p0)/a0m0;
	(void) printf("\navg-rho, avg-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);

	rho1 = Dens(st1);		
	rho2 = Dens(st2);
	a0 = 1.0;
	(void) printf("\npost-rho, pre-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);
	a0 = phy_param->a0p0/a0m0;
	(void) printf("\npost-rho, post-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);
	a0 = 0.5*(a0m0 + phy_param->a0p0)/a0m0;
	(void) printf("\npost-rho, avg-a0, \t a0_dot = %g\n",
		      u0*(rho1-rho2)/(rho1+rho2)*k*a0);

	(void) printf("\nII. Fraley's perturbation model ");
	(void) printf("(ideal gas EOS only):\n");
	g1 = gruneisen_gamma(st1) + 1.0;
	g2 = gruneisen_gamma(st2) + 1.0;
	a1 = phy_param->a1p0;	 
	a2 = phy_param->a2p0;
	s_t2 = 1.0 - Press(st0)/Press(st1);	s_t2 *= s_t2;
	s_r2 = 1.0 - Press(st2)/Press(st3);	s_r2 *= s_r2;
	tmp = 1 + (g1+1)*s_t2/4/(g1*g1*(2-s_t2)-s_t2);
	a1 /= tmp;
	tmp = 1 + (g2+1)*s_r2/4/(g2*g2*(2-s_r2)-s_r2);
	a2 /= tmp;
	a2 *= -Vel(st3)[0]/Vel(st0)[0];
	rho1 = Dens(st1);		
	rho2 = Dens(st2);
	tmp = (rho1*a1 - rho2*a2)/(rho1+rho2)/a0m0;
	(void) printf("\t \t a0_dot = %g\n",u0*k*tmp);
	
	(void) printf("\nIII. vortex sheet model:\n");
	a0 = 1.0;
	tmp = refracted_circulation_per_unit_len_per_sin_inc_ang(sta,st0,
		                                                 st1,st2,st3);

	(void) printf("\t \t a0_dot = %g\n",0.5*k*a0*tmp);
	(void) printf("----------------------------------------------------\n");
}		/*end empirical_growth_rates*/

LOCAL	double refracted_circulation_per_unit_len_per_sin_inc_ang(
	Locstate	sta,
	Locstate	st0,
	Locstate	st1,
	Locstate	st2,
	Locstate	st3)
{
	double		Ha, H0, H1, H2;
	double		rhoa, ma;

	rhoa = Dens(sta);
	ma =  mass_flux(pressure(st3),sta);
	Ha = specific_enthalpy(sta);
	H0 = specific_enthalpy(st0);
	H1 = specific_enthalpy(st1);
	H2 = specific_enthalpy(st2);

	return rhoa*(Ha - H2 + H1 - H0)/ma;
}		/*end refracted_circulation_per_unit_len_per_sin_inc_ang*/


LOCAL	NUM_PARAM	*setup_num_param(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	LAYER_SYS	*layer_sys,
	PHY_PARAM	*phy_param)
{
	Front	    *front = layer_sys->front;
	RECT_GRID   *gr = front->rect_grid;
	STOP	    *stop = initial_stop(init);
	int	    dim = gr->dim;
	int	    n;
	char	    s[Gets_BUF_SIZE];
	double	    w1 = phy_param->w1;
	double	    w2 = phy_param->w2;
	double	    w3 = phy_param->w3;
	double	    cfl, delx, cmax, t_end, dx1, dx2;
	double	    dels;
	NUM_PARAM   *num_param;
	static const char  *dname[] = {"x", "y", "z"};
	static const double THIN = 1.0e-5; /*TOLERANCE*/
	static const int   MIN_PTS = 4;   /*TOLERANCE*/

	scalar(&num_param,sizeof(NUM_PARAM));

	n = gr->gmax[dim-1];
	screen("\nEnter the number of points in the self-similar (%s/t) "
	       "interval between the\n\ttransmitted and reflected waves "
	       "(dflt = %d): ",dname[dim-1],n);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&n);
	check_int_input("",n,1,1,GE_);
	n1 = (int) (n*w1/(w1+w3));
	n1 = num_param->n1 = max(n1,MIN_PTS);
	dx1 = num_param->dx1 = w1/n1;	
	n2 = (int) (n*w2/(w1+w3));
	n2 = num_param->n2 = max(n2,MIN_PTS);
	dx2 = num_param->dx2 = w2/n2;	
	if ((rflct_wv_type == RAREFACTION) && (fabs(w3-w2) > THIN*fabs(w3)))
	{
	    n3 = (int) (n*(w3-w2)/(w1+w3));
	    n3 = num_param->n3 = max(n3,MIN_PTS);
	    num_param->dxr = (w3-w2)/n3;
	}
	else
	{
	    n3 = num_param->n3 = 0;
	    num_param->dxr = 0.0;
	}

	cfl = min(Time_step_factor(ip->root->front),0.1); /*TOLERANCE*/
	screen("Input the cfl factor (dflt = %g): ",cfl);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&cfl);
	check_float_input("",cfl,0.0,1.0,GE_AND_LE);

	delx = min(dx1,dx2);
	cmax = max(phy_param->c1sq,phy_param->c2sq);
	cmax = sqrt(cmax);
	num_param->dels = dels = min(cfl*delx/cmax,1.0);
	(void) printf("dels = %g\n\n",dels);

	num_param->s_start = 0.5*exp(dels); /*TOLERANCE*/
	screen("The unperturbed wave interaction is assumed "
	       "to occur at time 0,\n\t"
	       "the self-similar solver starts at a positive time t_start,  "
	       "enter the\n\t"
	       "positive time at which the linearized solution "
	       "should be initialized\n\t"
	       "(dflt = %g)\n",num_param->s_start);
	screen("Input t_start: ");
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&num_param->s_start);
	check_float_input("for t_start",num_param->s_start,0.0,0.0,GE_);

	num_param->s_end = num_param->s_start + 2.0*exp(dels); /*TOLERANCE*/
	screen("For small times,  the solver uses log(time) as the "
	       "timelike variable,\n\t"
	       "enter the time t_end at which the solver should "
	       "switch from solving in\n\t"
	       "log(time) to solving in real time (dflt = %g)\n",
	       num_param->s_end);
	screen("Input t_end: ");
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&num_param->s_end);
	check_float_input("for t_end",num_param->s_end,num_param->s_start,
			  num_param->s_start,GE_);

	num_param->s_start = log(num_param->s_start);
	num_param->s_end = log(num_param->s_end);

	t_end = stop->_stop_time;
	screen("Enter the termination time of the simulation, (dflt = %g): ",
	       t_end);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&t_end);
	num_param->t_end = t_end;
	check_float_input("",t_end,0.0,0.0,GE_);

	num_param->dt_prt_vel = t_end/100.0;
	num_param->dt_prt_prf = t_end/100.0;
	screen("Input time interval to print out velocity of a0 (dflt = %g): ",
	       num_param->dt_prt_vel);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&num_param->dt_prt_vel);

	screen("Input time interval to print out "
	       "p and w profiles (dflt = %g): ",num_param->dt_prt_prf);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&num_param->dt_prt_prf);

	num_param->smooth = 0.1; /*TOLERANCE*/
	screen("Input the smoothing factor (dflt = %g): ",
	       num_param->smooth);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&num_param->smooth);
	check_float_input("",num_param->smooth,0.0,0.5,GE_AND_LE);

	if (Prt_mode(prt_opts(init)[TRI_STATES]) != PRINTING_OFF)
	{
	    LAYER_SURF *contact = layer_sys->layer[1]->upper_surf;
	    num_param->tri_prt = YES;
	    num_param->nx = gr->gmax[0];
	    num_param->ny = gr->gmax[1];
	    num_param->lx = gr->GL[0];
	    num_param->ux = gr->GU[0];
	    num_param->ly = gr->GL[1];
	    num_param->uy = gr->GU[1];

	    num_param->y0 = contact->pbar[dim-1];
	    num_param->a0m0phys = contact->fpoly->A[0];

	    num_param->t_offset = 0.0;
	    screen("Enter time offset for synchronization of tri plots with\n"
	           "\ta nonlinear simulation (dflt = %g): ",
		   num_param->t_offset);
	    (void) Gets(s);
	    if (s[0] != '\0')
	        (void) sscan_float(s,&num_param->t_offset);
	}
	else
	{
	    num_param->tri_prt = NO;
	    num_param->t_offset = 0.0;
	}

	return	num_param;
}		/*end setup_num_param*/


LOCAL   void    info_inside_rarefaction_wave(
	PHY_PARAM	*phy_param,
	NUM_PARAM	*num_param)
{
	boolean  vacuum;
        int      i;
        double	 *g, *c, *rho;
        double	 speed;
	double	 spdnew;
        Locstate st, ahead = phy_param->st3, behind = phy_param->st2;
 
        (*Params(ahead)->_alloc_state)(&st,Params(ahead)->sizest);
        uni_array(&c,n3+1,sizeof(double));
        uni_array(&rho,n3+1,sizeof(double));
        uni_array(&g,n3+1,sizeof(double));
 
        for (i = 0; i <= n3; i++)
        {
            speed = phy_param->w2 + i*num_param->dxr;
	    vacuum = oned_state_in_rarefaction_fan(speed,Vel(ahead)[0],
				                   ahead,behind,st,TGAS_STATE,
						   &spdnew,RIGHT_FAMILY);
            if (vacuum == YES)
            {
                screen("ERROR in info_inside_rarefaction_wave(), "
                       "vacuum detected inside rarefaction wave!\n");
                clean_up(ERROR);
            }
            c[i] = sound_speed(st);
	    rho[i] = Dens(st);
            g[i] = fundamental_derivative(st);
        }
 
        phy_param->c = c;
        phy_param->rho = rho;
        phy_param->g = g;
 
        free(st);
}		/*end info_inside_rarefaction_wave*/


LOCAL	void	init_conditions(
	PHY_PARAM	*phy_param,
	NUM_PARAM	*num_param)
{
	int		i, n = n1 + n2 + 1, use_random = 0;
	double		t;
	char		choice[Gets_BUF_SIZE];

	uni_array(&pm,n,sizeof(double));
	uni_array(&p, n,sizeof(double));
	uni_array(&pp,n,sizeof(double));
	uni_array(&qm,n,sizeof(double));
	uni_array(&q, n,sizeof(double));
	uni_array(&qp,n,sizeof(double));
	/* making offsets */
	pm += n1;	p += n1;	pp += n1;
	qm += n1;	q += n1;	qp += n1;

	if ((rflct_wv_type == RAREFACTION) && (n3 > 0))
	{
	    uni_array(&Am,n3+1,sizeof(double));
	    uni_array(&A, n3+1,sizeof(double));
	    uni_array(&Ap,n3+1,sizeof(double));
	    uni_array(&Bm,n3+1,sizeof(double));
	    uni_array(&B, n3+1,sizeof(double));
	    uni_array(&Bp,n3+1,sizeof(double));
	    uni_array(&Vm,n3+1,sizeof(double));
	    uni_array(&V, n3+1,sizeof(double));
	    uni_array(&Vp,n3+1,sizeof(double));
	}

	screen("Type y to use random initial conditions for the pressure: ");
	(void) Gets(choice);
	if ((choice[0] == 'y') || (choice[0] == 'Y'))
	{
	    long	seed;

	    use_random = 1;
	    screen("Input the seed for the random number generator: ");
	    (void) Scanf("%d\n",&seed);
	    srand48(seed);
	}

	/* zeroth order amplitudes */ 
	if (use_random)		/* amp +- 25% */
	{
	    a1m = a1 = (1.0 + 0.25*(2.0*drand48()-1))*phy_param->a1p0;
	    a2m = a2 = (1.0 + 0.25*(2.0*drand48()-1))*phy_param->a2p0;
	    a3 = (1.0 + 0.25*(2.0*drand48()-1))*phy_param->a3p0;
	}
	else
	{
	    a1m = a1 = phy_param->a1p0;
	    a2m = a2 = phy_param->a2p0;
	    a3 = phy_param->a3p0;
	}

	/* zeroth order solution inside rarefaction */ 
	if ((rflct_wv_type == RAREFACTION) && (n3 > 0))
	{
	    double w2 = phy_param->w2,
	    	  ksq = phy_param->ksq, dxr = num_param->dxr, xi;

	    for (i = n3; i >= 0; i--)
	    {
	    	xi = w2 + i*dxr;
	    	Am[i] = A[i] = f_A0(xi,phy_param);
	    	Bm[i] = B[i] = 0.0;
	    	xi = w2 + (i+0.5)*dxr;
	    	if (i == n3)	
	    	    Vm[i] = V[i] = 0.0;
	    	else		
	    	    Vm[i] = V[i] = V[i+1] + ksq*f_A0(xi,phy_param)*dxr;
		}
	}
	if (use_random && (rflct_wv_type == RAREFACTION) && (n3 > 0))
	{
	    for (i = 1; i < n3; i++)
	    {
	    	Am[i] *= drand48(); A[i] *= drand48();
	    	Bm[i] *= drand48(); B[i] *= drand48();
	    	Vm[i] *= drand48(); V[i] *= drand48();
	    }
	}

	/* first order solution in regions 1 & 2 */ 
	if (use_random)
	{
	    for (i = -n1; i <= n2; i++)
	        pp[i] = drand48();
	    for (i = 0; n3 && i <= n3; i++) 	
	    {
	    	Ap[i] = drand48();
	    	Bp[i] = drand48();
	    	Vp[i] = drand48();
	    }
	    a2p = drand48();
	}
	else
	    small_t_solution(phy_param,num_param);

	/* values for the first step */
	t = exp(num_param->s_start);
	for (i = -n1; i <= n2; i++)	
	{
	    pm[i] = t*pp[i];
	    qm[i] = t*pp[i];
	}
	if ((rflct_wv_type == RAREFACTION) && (n3 > 0))
	{
	    for (i = 0; i <= n3; i++)
	    {
	    	Am[i] += t*t*Ap[i];
	    	Bm[i] += t*t*Bp[i];
	    	Vm[i] += t*t*Vp[i];
	    }
	    a2m += t*t*a2p;
	}
		

	/* values for the second step */
	t = exp(num_param->s_start + num_param->dels);
	for (i = -n1; i <= n2; i++)	
	{
	    p[i] = t*pp[i];
	    q[i] = t*pp[i];
	}
	if ((rflct_wv_type == RAREFACTION) && (n3 > 0))
	{
	    for (i = 0; i <= n3; i++)
	    {
	    	A[i] += t*t*Ap[i];
	    	B[i] += t*t*Bp[i];
	    	V[i] += t*t*Vp[i];
	    }
	    a2 += t*t*a2p;
	}
}		/*end init_conditions*/


LOCAL	double	f_A0(
	double xi,
	PHY_PARAM *phy_param)
{
	static	boolean	first = YES;
	static	double	slp, inc;

	if (first == YES)
	{
	    double w2, w3, A0w2, A0w3;

	    first = NO;
	    w2 = phy_param->w2;
	    w3 = phy_param->w3;
	    A0w2 = -a2/phy_param->g[0];
	    A0w3 = -a3/phy_param->g[n3];
	    slp = (A0w3-A0w2)/(w3-w2);
	    inc = -slp*w2 + A0w2;
	}

	return	slp*xi + inc;
}		/*end f_A0*/


LOCAL	void	small_t_solution(
	PHY_PARAM	*phy_param,
	NUM_PARAM	*num_param)
{
	int	i;
	double	dx1 = num_param->dx1, dx2 = num_param->dx2; 
	double	w1 = phy_param->w1, w2 = phy_param->w2, U0 = phy_param->U0,
		ksq = phy_param->ksq, a0m0 = phy_param->a0m0,
		c1sq = phy_param->c1sq, c2sq = phy_param->c2sq,
		u0 = phy_param->u0, u3 = phy_param->u3,
		v1 = phy_param->v1, v2 = phy_param->v2,
		deno1 = phy_param->deno1, deno2 = phy_param->deno2; 
	double	pn2, ca1, sigma, tau, tmp;
	double	beta1, beta2, p0;

	tmp = -deno1*(v2*w1+v1*w2) + (w1*w1-c1sq)*v2;
	sigma = -deno1*v1/tmp;
	ca1 = ksq*c1sq*u0*w1/v1*(1.0 - (u0 + w1)/U0)*a0m0;
	tau = v1*ca1/tmp;

	pn2 = 0.0;		
	if (rflct_wv_type == SHOCK)
	{
	    double ca2 = ksq*c2sq*u3*w2/v2*(1.0 - (u0 - w2)/U0)*a0m0;

	    tmp = deno2 - (w2*w2-c2sq)*sigma;
	    pn2 = ((w2*w2-c2sq)*tau + ca2)/tmp;
	}
	else if ((rflct_wv_type == RAREFACTION) && (n3 > 0))
	    pn2 = p21_at_w2(phy_param,num_param); 

	beta2 = sigma*pn2 + tau;
	beta1 = v2*beta2/v1;
	p0 = pn2 - beta2*w2;

	(void) printf("\nSmall t pressure amplitude = %g\n",p0/a0m0);
	(void) printf("Small t normal velocity amplitude = %g\n",v1*beta1/a0m0);

	for (i = -n1; i <= 0; i++)
	    pp[i] = beta1*i*dx1 + p0;

	for (i = 1; i <= n2; i++)
	    pp[i] = beta2*i*dx2 + p0;
}		/*end small_t_solution*/


LOCAL	double	p21_at_w2(
	PHY_PARAM	*phy_param,
	NUM_PARAM	*num_param)
{
	int		i;
	double		*g = phy_param->g, *c = phy_param->c;
	double		*rho = phy_param->rho;
	double		ksq = phy_param->ksq, dxr = num_param->dxr;

	/* temporary use Ap, Bp and Vp as A^(2), B^(2) and V^(2), etc. */
	Bp[n3] = Bp[n3-1] = 0.0;
	Ap[n3] = 0.0;
	Ap[n3-1] = 0.5*(1.0-2.0/g[n3-1])*Bp[n3-1] - 0.25*c[n3-1]*V[n3-1];
	Vp[n3] = Vp[n3-1] = 0.0;
	for (i = n3-1; i > 0; i--)
	{
	    Bp[i-1] = Bp[i+1] - 2*dxr*(0.5/c[i]*Bp[i] - 0.25*V[i]);
	    Ap[i-1] = 0.5*(1.0-2.0/g[i-1])*Bp[i-1] - 0.25*c[i-1]*V[i-1];
	    Vp[i-1] = Vp[i+1] - 2*dxr*(2.0/c[i]*Vp[i] - ksq*(Ap[i]-Bp[i]));
	}
	a2p = g[0]*Ap[0] + (2.0 - g[0])*Bp[0];

	return	-rho[0]*c[0]*(Bp[0] + 0.5*c[0]*V[0]);
}		/*end p21_at_w2*/


#define	Smoothing(a,b,c)						\
    (num_param->smooth*(a)+num_param->smooth*(c)+(1-2*num_param->smooth)*(b))

LOCAL	void	march_forward(
	INTERFACE *intfc,
	PHY_PARAM *phy_param,
	NUM_PARAM *num_param)
{
	int	i, is_s_part = 1, j = 0, n;
	double	dels = num_param->dels, t_end = num_param->t_end;
	double	dt_prt_prf = num_param->dt_prt_prf, t_prt_prf = 0;
	double	dt_prt_vel = num_param->dt_prt_vel, t_prt_vel = 0;
	double	*tp, expdelsm1, t, s, delt, del; 
	double	a0m0 = phy_param->a0m0, amp0 = phy_param->a0p0;
	double	*rtime, *ramp0,*ramp1,*ramp2, *radot0,*racc0, *p_w1,*p_w2,*p_0;
	double	wm0, w0, wp0, alpha1 = phy_param->v1/num_param->dx1/2; 
	double	t_offset = num_param->t_offset;

	n = ((int)(num_param->t_end/dt_prt_vel)) + 2;
	uni_array(&rtime, n,sizeof(double));
	uni_array(&ramp0, n,sizeof(double));
	uni_array(&ramp1, n,sizeof(double));
	uni_array(&ramp2, n,sizeof(double));
	uni_array(&radot0,n,sizeof(double));
	uni_array(&racc0, n,sizeof(double));
	uni_array(&p_w1,  n,sizeof(double));
	uni_array(&p_w2,  n,sizeof(double));
	uni_array(&p_0,   n,sizeof(double)); 

	s = num_param->s_start;		t = exp(s);
	wm0 = alpha1*(3*pm[0] - 4*pm[-1] + pm[-2]);
	s += dels;			t = exp(s);
	w0 = alpha1*(3*p[0] - 4*p[-1] + p[-2]);

	expdelsm1 = -1.0 + exp(dels);
	delt = t*expdelsm1;
	while (t <= t_end)
	{
	    if (s < num_param->s_end)
	    {
	    	s += dels;
	    	delt = t*expdelsm1;
	    }
	    else if (is_s_part == 1)
	    {
	    	is_s_part = 0;
	    	for (i = -n1; i <= n2; i++)
	    	    q[i] /= t,	qm[i] /= t - delt;
	    }

	    del = is_s_part ? dels : delt;
	    update(t,del,is_s_part,phy_param,num_param);
	    wp0 = -(is_s_part ? 1.0 : 1.0/t)*
		  alpha1*(3*p[0]-4*p[-1]+p[-2])*del*2 + wm0;

	    /* temporal smoothing operation */
	    for (i = -n1; i <= n2; i++)
	    {
	    	p[i] = Smoothing(pm[i],p[i],pp[i]);
	    	q[i] = Smoothing(qm[i],q[i],qp[i]);
	    }
	    a1 = Smoothing(a1m,a1,a1p);
	    a2 = Smoothing(a2m,a2,a2p);

	    /* calculate amplitude */
	    amp0 += 0.5*(wm0 + w0)*delt;

	    /* printing */
	    if (t > t_prt_vel + t_offset)
	    {
	    	rtime[j] = t;
	    	ramp0[j] = amp0/a0m0;
	    	ramp1[j] = a1/a0m0;
	    	ramp2[j] = a2/a0m0;
	    	radot0[j] = w0/a0m0;
	    	racc0[j] = (wp0 - wm0)/2.0/del/a0m0;
	    	p_0[j] = p[0]/a0m0;
	    	p_w1[j] = p[-n1]/a0m0;
	    	p_w2[j] = p[n2]/a0m0;
	    	if (is_s_part)
		    racc0[j] /= t;
	    	j++;
		t_prt_vel += dt_prt_vel;
	    }
	    if (t > t_prt_prf + t_offset)
	    {
	    	double	x = -phy_param->w1, dxr = num_param->dxr;
	    	double	*c = phy_param->c, *rho = phy_param->rho;

	    	(void) printf("\n\ntime = %20.8g\n",t);
	    	(void) output();
	    	(void) printf("%20s %20s\n","xi","del_p");
	    	for (i = -n1; i <= n2; i++)
	    	{
	    	    if (t_end < 1.0)
			(void) printf("%20.8g %20.8g\n",x,p[i]/t/a0m0);
		    else
			(void) printf("%20.8g %20.8g\n",x,p[i]/a0m0);
		    x += (i < 0) ? num_param->dx1 : num_param->dx2;
		}
		x -= num_param->dx2;
		for (i = 0; n3 && i <= n3; i++)
		{
		    (void) printf("%20.8g %20.8g\n",
				  x,rho[i]*c[i]*(A[i] - B[i])/t/a0m0);
		    x += dxr;
		}
	    	if (n3 > 0)	/* solution inside rarefaction */
	    	{
	    	    (void) output();
	    	    (void) printf("%20s %20s\n","xi","del_A");
	    	    for (i = 0; i <= n3; i++)
	    	    	(void) printf("%20.8g %20.8g\n",i*dxr,A[i]/a0m0);
	    	    (void) output();
	    	    (void) printf("%20s %20s\n","xi","del_B");
	    	    for (i = 0; i <= n3; i++)
	    	    	(void) printf("%20.8g %20.8g\n",i*dxr,B[i]/a0m0);
	    	    (void) output();
	    	    (void) printf("%20s %20s\n","xi","del_V");
	    	    for (i = 0; i <= n3; i++)
	    	    	(void) printf("%20.8g %20.8g\n",i*dxr,V[i]/a0m0);
	    	}

	    	/* Print tri grid if requested */
	    	if (num_param->tri_prt)
	    	    rm_linear_print_tri_soln(intfc,t,amp0,phy_param,num_param);

	    	t_prt_prf += dt_prt_prf;
	    }

	    if (is_s_part == 1)
		t = exp(s);
	    else
		t += delt;

	    tp  = pm, pm = p,	p = pp,	pp = tp;	
	    tp  = qm, qm = q,	q = qp,	qp = tp;	
	    wm0 = w0, w0 = wp0;
	    a1m = a1, a1 = a1p;
	    a2m = a2, a2 = a2p;
	    if ((rflct_wv_type == RAREFACTION) && (n3 > 0))
	    {
	    	tp = Am, Am = A, A = Ap, Ap = tp;
	    	tp = Bm, Bm = B, B = Bp, Bp = tp;
	    	tp = Vm, Vm = V, V = Vp, Vp = tp;
	    }
	}

	(void) output();
	(void) printf("%14s%14s%14s%14s%14s%14s%14s%14s%14s\n","time",
		      "a0","a0_dot","a0_ddot","a1","a2","p_w1","p_w2","p_0");
	for (i = 0; i < j; i++)
	    (void) printf("%14g%14g%14g%14g%14g%14g%14g%14g%14g\n",
			  rtime[i],ramp0[i],radot0[i],racc0[i],
			  ramp1[i],ramp2[i],p_w1[i],p_w2[i],p_0[i]);
}		/*end march_forward*/


#define	ldf(f,i,dx)		(-0.5*(3*f[i] - 4*f[i+1] + f[i+2])/(dx))
#define	dfdx(f,i,dx)		((f[i+1] - f[i-1]) / (dx) /2.0)
#define	d2fdx2(f,i,dxsq)	((f[i+1] + f[i-1] - 2*f[i]) / (dxsq))

LOCAL	void	update(
	double		t,
	double		del,
	int		is_s_part,
	PHY_PARAM	*phy_param,
	NUM_PARAM	*num_param)
{
	int		i;
	double		ft2, ft, o_or_z;
	double		x, tmp;
	static	boolean	first = YES;
	static	double	dx1, dx2, dx1sq, dx2sq;
	static	double	c1sq, c2sq, v1, v2, ksq, alpha1, alpha2; 
	static	double	cts1, cts2, cts3, crs1 = 0.0, crs2 = 0.0, crs3 = 0.0; 

	if (first)
	{
	    double w1 = phy_param->w1, w2 = phy_param->w2; 
	    double K1 = phy_param->K1, K2 = phy_param->K2,
		  u0 = phy_param->u0, u3 = phy_param->u3,
		  v0 = phy_param->v0, v3 = phy_param->v3,
		  deno1 = phy_param->deno1, deno2 = phy_param->deno2; 

	    first = NO;
	    dx1 = num_param->dx1; 
	    dx2 = num_param->dx2; 
	    dx1sq = dx1*dx1, dx2sq = dx2*dx2;
	    ksq = phy_param->ksq;
	    c1sq = phy_param->c1sq, 	c2sq = phy_param->c2sq;
	    v1 = phy_param->v1,		v2 = phy_param->v2; 
	    alpha1 = v1/dx1/2;		alpha2 = v2/dx2/2;
	    cts1 = -0.5*v0*v1/(v0-v1)*(1/w1-w1*K1/c1sq);
	    cts2 = (w1*w1 - c1sq)/deno1;
	    cts3 = ksq*c1sq*u0*w1/v1/deno1;
	    if ((rflct_wv_type == SHOCK) && (v2 != v3))
	    {
	    	crs1 = 0.5*v3*v2/(v3-v2)*(1/w2-w2*K2/c2sq);
	    	crs2 = (w2*w2 - c2sq)/deno2;
	    	crs3 = ksq*c2sq*u3*w2/v2/deno2;
	    }

	}


	if (is_s_part)		
	{
	    ft2 = 1.0;	ft = 1.0;	o_or_z = 1.0;
	}
	else 		
	{
	    ft2 = 1.0/t/t;	ft = 1.0/t;	o_or_z = 0.0;
	}


	/* between the transmitted shock and the contact surface */
	for (i = -n1+1; i < 0; i++)
	{
	    x = i*dx1;
	    qp[i] = (ft2*((c1sq-x*x)*d2fdx2(p,i,dx1sq) - 2*x*dfdx(p,i,dx1) 
		     - t*t*ksq*c1sq*p[i]) + ft*2*x*dfdx(q,i,dx1) + o_or_z*q[i])*
		    2*del + qm[i];
	    pp[i] = q[i]*2*del + pm[i];
	}


	/* between cont. surf. and reflected shock or TE of rarefaction wave */
	for (i = 1; i < n2; i++)
	{
	    x = i*dx2;
	    qp[i] = (ft2*((c2sq-x*x)*d2fdx2(p,i,dx2sq) - 2*x*dfdx(p,i,dx2) 
		     - t*t*ksq*c2sq*p[i]) + ft*2*x*dfdx(q,i,dx2) + o_or_z*q[i])*
		    2*del + qm[i];
	    pp[i] = q[i]*2*del + pm[i];
	}


	/* on the contact surface */
	pp[0] = -(alpha2*(-4*pp[1] + pp[2]) + alpha1*(-4*pp[-1] + pp[-2]))/
		 (alpha1 + alpha2)/3;
	qp[0] = (pp[0] - pm[0])/del - qm[0];


	/* on the transmitted shock */
	a1p = t*ft*cts1*p[-n1]*2*del + a1m;
	pp[-n1] = q[-n1]*2*del + pm[-n1];
	tmp = 0.5*(3*pp[-n1] - 4*pp[-n1+1] + pp[-n1+2])/dx1;
	qp[-n1] = ft*(cts2*tmp + t*cts3*a1p);


	/* on the reflected shock or trailing edge of RR */
	if ((rflct_wv_type == SHOCK) || (n3 == 0)) /* thin rarefaction fan */
	{
	    a2p = t*ft*crs1*p[n2]*2*del + a2m;
	    pp[n2] = q[n2]*2*del + pm[n2];
	    tmp = 0.5*(3*pp[n2] - 4*pp[n2-1] + pp[n2-2])/dx2;
	    qp[n2] = ft*(crs2*tmp + t*crs3*a2p);
	}
	else 
	{
	    double g0 = phy_param->g[0], c0 = phy_param->c[0],
	    	  rho0 = phy_param->rho[0];

	    update_inside_rr(t,del,is_s_part,phy_param,num_param);
	    a2p = ft*(a2 + g0*A[0] + (2-g0)*B[0])*2*del + a2m;
	    pp[n2] = c0*rho0*(Ap[0] - Bp[0] + 1.0/g0*a2p)/t;
	    qp[n2] = (3*pp[n2] -4*p[n2] + pm[n2])/2.0/del;
	    /*
	    tmp = -0.5*(3*Bp[0] - 4*Bp[1] + Bp[2])/dxr;
	    qp[n2] = -ft*rho0*c0*c0*(2.0*tmp/t + Vp[0]*t);
	    */
	}
}		/*end update*/
	

#define	Interpolate(f,i,frac)		((1-(frac))*f[i] + (frac)*f[i+1])

LOCAL	void	update_inside_rr(
	double		t,
	double		del,
	int		is_s_part,
	PHY_PARAM	*phy_param,
	NUM_PARAM	*num_param)
{
	int		i, j;
	double		ft = is_s_part ? 1.0 : 1.0/t;
	double		cs, xi, A0, B0, V0, c0, dot, h, frac, tmp;
	double		*g = phy_param->g, *c = phy_param->c;
	double		ksq = phy_param->ksq, dxr = num_param->dxr;
	double		end = phy_param->w3 - phy_param->w2;

	/* update A */
	Ap[n3] = -1.0/g[n3]*a3;
	for (i = 0; i < n3; i++)
	    Ap[i] = Am[i] + 2*del*ft*((1.0-2.0/g[i])*B[i] - 0.5*c[i]*t*t*V[i]);	

	/* update B */
	Bp[n3] = 0.0;
	for (i = 0; i < n3; i++)
	{
	    Bp[i] = 0.0;
	    xi = i*dxr;
	    cs = ft*2*c[i];
	    xi += cs*del;
	    if (xi < end)
	    {
	    	tmp = xi/dxr; j = (int)(floor(tmp)); frac = tmp - j;
	    	c0 = Interpolate(c,j,frac);
	    	B0 = Interpolate(B,j,frac);
	    	V0 = Interpolate(V,j,frac);
	    	dot = ft*(B0 + 0.5*c0*t*t*V0);
	    	cs = ft*2*c0;
	    	xi += cs*del;
	    	if (xi < end)
	    	{
	    	    tmp = xi/dxr; j = (int)(floor(tmp)); frac = tmp - j;
	    	    B0 = Interpolate(Bm,j,frac);
	    	    Bp[i] = B0 + 2*del*dot;
	    	}
	    }
	}

	/* update V */
	Vp[n3] = 0.0;
	for (i = 0; i < n3; i++)
	{
	    xi = i*dxr;
	    cs = ft*c[i];
	    xi += cs*del;
	    if (xi < end)
	    {
	    	tmp = xi/dxr; j = (int)(floor(tmp)); frac = tmp - j;
	    	c0 = Interpolate(c,j,frac);
	    	A0 = Interpolate(A,j,frac);
	    	B0 = Interpolate(B,j,frac);
	    	dot = ft*ksq*c0*(A0 - B0);
	    	cs = ft*c0;
	    	xi += cs*del;
	    	if (xi < end)
	    	{
	    	    tmp = xi/dxr; j = (int)(floor(tmp)); frac = tmp - j;
		    V0 = Interpolate(Vm,j,frac);
		    Vp[i] = V0 + 2*del*dot;
		}
		else 
		{
		    h = del - (xi - end)/cs;
		    dot = ft*ksq*c[n3]*(-a3/g[n3]);
		    Vp[i] = (h+del)*dot;
		}
	    }
	    else
	    {
	    	h = del - (xi - end)/cs;
	    	dot = ft*ksq*c[n3]*(-a3/g[n3]);
	    	Vp[i] = h*dot;
	    }
	}
}		/*end update_inside_rr*/

/*
* 			rm_linear_print_tri_soln():
*
*	Prints pressure using TRI_SOLN format for plotting with tri.
*
*	Note: The routine tests to make sure that the perturbations
*	on the transmitted and reflected waves don't intersect with the
*	perturbed contact.  With sufficiently large perturbations this
*	can happen at early enough
*	times.  If this is the case there will be a warning and no data
*	will be printed.
*
*	RLH 1-18-94
*/

LOCAL void rm_linear_print_tri_soln(
	INTERFACE *intfc,
	double	  t,
	double	  amp0,
	PHY_PARAM *phy_param,
	NUM_PARAM *num_param)
{
	double lx = num_param->lx, /* Domain size */
	      ly = num_param->ly;
	double ux = num_param->ux,
	      uy = num_param->uy;
	double y0 = num_param->y0; /* y0--avg. contact posn. */
	double dx, dy;
	double w1 = phy_param->w1, /* w1, w2, w3--wave speeds */
	      w2 = phy_param->w2,
	      w3 = phy_param->w3;
	double xi;
	double x, y, cos_kx;
	double p0 = Press(phy_param->st0), /* p0, p1, p3--ahead, middle and */
	      p1 = Press(phy_param->st1), /* behind pressures from         */
	      p3 = Press(phy_param->st3); /* Riemann problem               */
	double a0m0phys = num_param->a0m0phys;
	double amp1 = a0m0phys*a1, /* amp1, amp2, amp3, a0 = scaled amplitudes */
	      amp2 = a0m0phys*a2,
	      amp3 = a0m0phys*a3;
	double a0 = a0m0phys*amp0;
	int   i, j;
	int   nx = num_param->nx, ny = num_param->ny; /* Number of zones */
	int	num_nodes, offset;

	if (debugging("tri_lin"))
	{
	    (void) printf("\na0 = %g\namp1 = %g\tamp2 = %g\tamp3 = %g\n",
			  a0,amp1,amp2,amp3);
	}

	(void) output();
	(void) printf("TRI_SOLN:  PRESSURE\n");
	(void) printf("#Point Source Data\n");
	(void) printf("#0\n");
	(void) printf("#Clipping Boundary\tCXL = %g CXU = %g "
		      "CYL = %g CYU = %g\n",lx,ux,ly,uy);

	/* Self intersections?  If so don't print any data */
	if ((fabs(amp1 - a0) > w1*t) ||
	    (fabs(amp2 - a0) > w2*t) ||
           ((rflct_wv_type==RAREFACTION) && (fabs(amp3-amp2) > (w3-w2)*t)))
	{
	    (void) printf("NODES 0 1 0 1\n");
	    (void) printf("RECTANGLES 0 0 0 1\n");
	    (void) printf("TRIANGLES 0 0 0 1\n");
	    (void) printf("END\n");
	    return;
	}

	dx = (ux - lx)/(nx - 1);
	dy = (uy - ly)/(ny - 1);
	num_nodes = nx * ny;

	(void) printf("NODES %d 1 0 1\n",num_nodes);

	x = lx;

	/* Print node data from bottom to top, left to right */
	for (i = 0; i < nx; i++, x += dx)
	{
	    y = ly;

	    /* cos_kx is the periodic term in the perurbed quantities.
	    *  That is, q(x,y,t) = q(y,t) + dq(y,t)*cos_kx
	    */
	    cos_kx = cos(2.0*PI/(ux - lx) * (x - lx));

	    /* Print pressure ahead of transmitted shock */
	    for (j = 0; j < ny && y < y0 - w1*t - amp1*cos_kx; j++, y += dy)
	    	(void) printf("%g %g %g\n",x,y,p0);

	    /* Print pressure between transmitted shock and contact */
	    for (; j < ny && y < y0 - a0*cos_kx; j++, y += dy)
	    {
	    	xi = (y - y0 + a0*cos_kx)*w1 / (w1*t - (a0 - amp1)*cos_kx);
	    	(void) printf("%g %g %g\n",x,y,
			      p1 - cos_kx*find_dp(xi,t,p,phy_param,num_param));
	    }

	    /* Print pressure between contact and reflected shock 
							or trailing edge */
	    for (; j < ny && y < y0 + w2*t - amp2*cos_kx; j++, y += dy)
	    {
	    	xi = (y - y0 + a0*cos_kx)*w2 / (w2*t - (amp2 - a0)*cos_kx);
	    	(void) printf("%g %g %g\n",x,y,
			      p1 - cos_kx*find_dp(xi,t,p,phy_param,num_param));
	    }
		
	    if (rflct_wv_type == RAREFACTION)
	    {
	    	/* Print pressure between rarefaction edges */
	    	for (; j < ny &&  y < y0 + w3*t - amp3*cos_kx; j++, y += dy)
	    	{
	    	    xi = w2 + (w3 - w2)*(y - (y0 + w2*t - amp2*cos_kx))/
	    	    	((w3 - w2)*t + (amp2 - amp3)*cos_kx);
	    	    (void) printf("%g %g %g\n",x,y,
		        p1 - cos_kx*find_dp(xi,t,p,phy_param,num_param));
	    	}
	    }

	    /* Finally, print pressure ahead of reflected wave */
	    for (; j < ny; j++, y += dy)
	    	(void) printf("%g %g %g\n",x,y,p3);
	}

	/* List nodes making up rectangles */
	(void) printf("RECTANGLES %d 0 0 1\n",(nx-1)*(ny-1));

	for (j = 0; j < ny-1; j++)
	{
	    for (i = 0; i < nx-1; i++)
	    {
	        offset = i*ny;
	        (void) printf("%d %d %d %d\n",j + offset,j + offset + ny,
		    	      j + offset +  1,j + offset + ny + 1);
	    }
	}

	/* We only have rectangles--no triangles */
	(void) printf("TRIANGLES 0 0 0 1\n");
	(void) printf("END\n");


	/* Print interface for merging with the tri solution plot */

	rm_linear_print_interface(intfc,y0,t,a0,amp1,amp2,amp3,
				  num_param,phy_param,rflct_wv_type);
}		/*end rm_linear_print_tri_soln*/

/*
*			find_dp():
*
*	For rm_linear_print_tri_soln().  Given xi (which is y/t) find
*	an spatially-interpolated pressure perturbation at time t.
*	For the region between the transmitted shock and either a reflected
*	shock or a rarefaction trailing edge the perturbation is found by
*	looking in the array p[].  Within a rarefaction the formula
*   	dp = rho*c*(A-B)/t		is used.
*
*/

LOCAL double find_dp(
	double		xi,
	double		t,
	double		*p,
	PHY_PARAM	*phy_param,
	NUM_PARAM	*num_param)
{
	double		dx, alpha, dp1, dp2;
	int		m,n;
	double		*c = phy_param->c, *rho = phy_param->rho;

	if ( (rflct_wv_type != RAREFACTION) ||
	    ((rflct_wv_type == RAREFACTION) && (xi < phy_param->w2)))
	{
	    /* Above (dx>0) || below (dx<0) contact? */
	    dx = (xi<0) ? num_param->dx1 : num_param->dx2;
		
	    /* Find perturbations for interpolation */
	    n = (int)(floor(xi/dx));
	    alpha = xi/dx - n;

	    dp1 = p[n];	dp2 = p[n+1];

	    /* Return linear interpolation between dp1 and dp2 */
	    return ((dp2 - dp1)*alpha*dx + dp1)*num_param->a0m0phys;
	}
	else 
	{
	    xi -= phy_param->w2;
	    dx = num_param->dxr;

	    n = (int)(floor(xi/dx));
	    alpha = xi/dx - n;

	    m = n + 1;
	    dp1 = rho[n]*c[n]*(A[n] - B[n]);
	    dp2 = rho[m]*c[m]*(A[m] - B[m]);

	    return ((dp2 - dp1)*alpha*dx + dp1)/t*num_param->a0m0phys;
	}
}		/*end find_dp*/

LOCAL void rm_linear_print_interface(
	INTERFACE *phys_intfc,
	double	  y0,
	double	  t,
	double	  a0,
	double	  amp1,
	double	  amp2,
	double	  amp3,
	NUM_PARAM *num_param,
	PHY_PARAM *phy_param,
	int	  r_wv_type)
{
	INTERFACE	*intfc;
	NODE		*nstart, *n1, *n2;
	double		lx = num_param->lx, ux = num_param->ux;
	double		ly = num_param->ly, uy = num_param->uy;
	double		coords1[2], coords2[2];
	static	double	phase, A;
	static	double	**nu = NULL;
	static	int	num_points = 0;
	static	boolean	first = YES;
	static	FOURIER_POLY	fpoly;
	static	RECT_GRID	top_grid, comp_grid;

	if (first)
	{
	    REMAP *remap = &topological_grid(phys_intfc).Remap;

	    int	gmax[2];
	    first = NO;
	    bi_array(&nu,2,1,FLOAT);

	    nu[0][0] = 2.0*PI/(ux - lx);
	    nu[1][0] = 1.0;
	    phase = 2.0*PI*lx/(ux - lx) - PI/2.0;

	    fpoly.num_modes = 1;
	    fpoly.dim = 2;
	    fpoly.nu = nu;
	    fpoly.phase = &phase;
	    fpoly.A = &A;

	    num_points = 2*num_param->nx;

	    coords1[0] = lx; coords1[1] = ly;
	    coords2[0] = ux; coords2[1] = uy;
	    gmax[0] = num_param->nx;
	    gmax[1] = num_param->ny;
	
	    set_rect_grid(coords1,coords2,coords1,coords2,
	    		  NOBUF,NOBUF,gmax,2,remap,&top_grid);
	    set_rect_grid(coords1,coords2,coords1,coords2,
	    		  NOBUF,NOBUF,gmax,2,remap,&comp_grid);

	    set_binary_output(NO);
	    (void) printf("\n\n\n");
	    (void) output();
	    (void) printf("\t\t\tINITIAL DATA:\n\n\n");
	    print_rectangular_grid(&comp_grid);

	    (void) printf("\n       stop_time = %-10g               "
			  "stop_step = %-10d\n",0.0,0);
	    (void) printf("\n\t\tComputational Area = %g\n",
			  (ux - lx)*(uy - ly));
	    (void) printf("\n\t\tRemap Geometry:  IDENTITY_REMAP\n");

	    (void) printf("\n\t\tPrinting Interval:  %g real time units",0.0);
	    (void) printf("\n\n\n\n\n");
	}

	intfc = make_interface(2);
	set_current_interface(intfc);

	set_topological_grid(intfc,&top_grid);
	set_computational_grid(intfc,&comp_grid);

	/* Create contact */
	fpoly.z0 = y0;
	A = a0;

	g_make_fourier_curve(CONTACT,num_points,lx,ux,&fpoly,
			     EXTERIOR_COMP,EXTERIOR_COMP,0.0);

	/* Create transmitted shock */
	fpoly.z0 = y0 - phy_param->w1*t;
	A = amp1;

	g_make_fourier_curve(FORWARD_SHOCK_WAVE,num_points,lx,ux,&fpoly,
			     EXTERIOR_COMP,EXTERIOR_COMP,0.0);

	/* Create reflected wave(s) */

	fpoly.z0 = y0 + phy_param->w2*t;
	A = amp2;

	if (r_wv_type != RAREFACTION)
	    g_make_fourier_curve(BACKWARD_SHOCK_WAVE,num_points,
			lx,ux,&fpoly,EXTERIOR_COMP,EXTERIOR_COMP,0.0);

	else
	{
	    g_make_fourier_curve(BACKWARD_SOUND_WAVE_TE,num_points,
			         lx,ux,&fpoly,EXTERIOR_COMP,EXTERIOR_COMP,0.0);
		
	    fpoly.z0 = y0 + phy_param->w3*t;
	    A = amp3;
		
	    g_make_fourier_curve(BACKWARD_SOUND_WAVE_LE,num_points,
			         lx,ux,&fpoly,EXTERIOR_COMP,EXTERIOR_COMP,0.0);
	}

	/* Make boundary curves */

	coords1[0] = lx; coords1[1] = ly;
	coords2[0] = ux; coords2[1] = ly;

	nstart = n1 = make_node(Point(coords1));
	n2 = make_node(Point(coords2));

	/* bottom */
	(void) make_curve(EXTERIOR_COMP,EXTERIOR_COMP,n1,n2);

	coords1[0] = ux; coords1[1] = uy;
	n1 = make_node(Point(coords1));

	/* right */
	(void) make_curve(EXTERIOR_COMP,EXTERIOR_COMP,n2,n1);

	coords2[0] = lx; coords2[1] = uy;
	n2 = make_node(Point(coords2));

	/* top */
	(void) make_curve(EXTERIOR_COMP,EXTERIOR_COMP,n1,n2);

	/* left */
	(void) make_curve(EXTERIOR_COMP,EXTERIOR_COMP,n2,nstart);

	(void) output();
	(void) printf("\tTIME DATA:  t = %g\n",t);

	(void) output();
	(void) printf("\t\t\tFRONT DATA:\n");

	print_interface(intfc);
	(void) delete_interface(intfc);
}		/*end rm_linear_print_interface*/
#endif /* defined(TWOD) */
