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
*				gsesphase.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Uses the SESAME and GRIZZLY tabular equation of state to initialize
*	hyp solution functions for the evaluation of the equation
*	of state for gas dynamics with multiple phase materials.
*/

#if defined(TWOD) && defined(PHASE_CODE)

#define DEBUG_STRING    "ses_hyp"

#include <geos/sesame.h>

	/* LOCAL Function Declarations */
LOCAL	double	set_S_state(double*,double,double,int,SESAME_EOS*);
LOCAL	void	boundary_states(int,double,double*,Front*,SESAME_EOS*,COLD_CURVE*,
				PHASE_BDRY*,int);
LOCAL	void	clear_utility_vectors(DEPARAMS*);
LOCAL	void	get_cold_data(SESAME_EOS*,COLD_CURVE*,double,double*,double*);
LOCAL	void	get_rho_at_crit_pt(PHASE_BDRY*,int,double*);
LOCAL	void	init_PE_cross_grid(Wave*,Front*,SESAME_EOS*,int,double*,double*);
LOCAL	void	insert_density_crossings(Front*,Wave*,SESAME_EOS*,COLD_CURVE*,
					 PHASE_BDRY*);
LOCAL	void	lookspl2(double,SESAME_EOS*,double*,double*);
LOCAL	void	rho_T_phase_initializer(Front*,double*,COMPONENT,Locstate,
					SESAME_EOS*,COLD_CURVE*,
					double,double*,double*);
LOCAL	void	search(double,SESAME_EOS*,int,int,int*);
LOCAL	void	set_S_on_pb(double,double,SESAME_EOS*,PHASE_BDRY*,double*,int);
LOCAL	void	set_entropy_on_ph_bnd(SESAME_EOS*,PHASE_BDRY*,Front*);
LOCAL	void	set_rho_T_phase_state(Locstate,double*,SESAME_EOS*,
				      COLD_CURVE*,double,double*,double*);
LOCAL	void	set_sinit(double*,double,SESAME_EOS*,PHASE_BDRY*);
LOCAL	void	set_vapor_dome(int,int,double,double,double,double,double,double,
			       double,double,NODE*,NODE*,double*,double*,double*,
			       double*,double*,double*,Front*,SESAME_EOS*);

LOCAL const double TOLEPS = 1.0e-05; /*TOLERANCE*/


LOCAL void rho_T_phase_initializer(
	Front		*front,
	double		*coords,
	COMPONENT	comp,
	Locstate	state,
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	double		sinit,
	double		*pp,
	double		*EE)
{
	size_t		sizest = sizeof(SES_RT_STATE);

	if (is_exterior_comp(comp,front->interf))
	{
		clear_state(front->interf,state,sizest);
	}
	else
		set_rho_T_phase_state(state,coords,seos,cold_curve,
				      sinit,pp,EE);
}		/*end rho_T_phase_initializer*/


LOCAL void set_rho_T_phase_state(
	Locstate	state,
	double		*coords,
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	double		sinit,
	double		*pp,
	double		*EE)
{

	double		rho_grid = coords[0];
	double		T_grid = coords[1];
	double		rho, T;
	double		p, e, cold_p, cold_e;

	rho = ses_rt_rho_from_grid(rho_grid,seos);
	T = ses_rt_temp_from_grid(T_grid,seos);
	if ((fabs(rho - Rho_min(seos)) < TOLEPS) &&
	    (fabs(T - Temp_min(seos))  < TOLEPS*(Temp_min(seos))))
	{
		seos->pmin = pp[0];
	}
	get_cold_data(seos,cold_curve,rho,&cold_e, &cold_p);
	p = pp[0];	
	e = EE[0];
	ses_rt_coldp(state) = cold_p;
	ses_rt_colde(state) = cold_e;
	ses_rt_redp(state) = (p - cold_p)/(rho*T);
	ses_rt_rede(state) = (e - cold_e)/T;
	ses_rt_adb_gam(state) = pp[1]/EE[1];
	ses_rt_gru_gam(state) = pp[1]/(EE[1]*rho);
	ses_rt_S(state) = sinit;
	Entropy_min(seos) = min(sinit,Entropy_min(seos));
	Entropy_max(seos) = max(sinit,Entropy_max(seos));
	Pressure_min(seos) = min(p,Pressure_min(seos));
	Pressure_max(seos) = max(p,Pressure_max(seos));
	Energy_min(seos) = min(e,Energy_min(seos));
	Energy_max(seos) = max(e,Energy_max(seos));
	if (DEBUG)
	{
		(void) printf("EOS - rho %g log(rho) %g, T %g log(T) %g\n",
			rho,log(rho),T,log(T));
		(void) printf("\tp %g\t\tlog(p) %g\n",p,log(p));
		(void) printf("\tcold p %g\t\tred. p %g\n",
			      cold_p,ses_rt_redp(state));
		(void) printf("\te %g\t\tlog(e+1) %g\n",e,log(e+1));
		(void) printf("\tcold e %g\t\tred. e %g\n",
			      cold_e,ses_rt_rede(state));
		(void) printf("\tsinit %g\n",sinit);
	}

}		/*end set_rho_T_phase_state*/



/*
*			init_RT_interior_states():
*
*	Initializes the states in a wave structure by calling
*
*		rho_T_phase_initializer(coords,comp,state,params)
*
*	at the centers of the grid blocks of wave->rect_grid.
*/

EXPORT void init_RT_interior_states(
	Wave		*wave,
	Front		*front,
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	PHASE_BDRY	*phase_bound)
{
	INTERFACE	*intfc;
	double		*tbls = seos->sestab.tbls;
	int		nt;
	int		nthyp = Ntemp_hyp(seos);
	int		ix,iy;
	size_t		sizest = front->sizest;
	int		place;
	int		ntv, nr, n4, flag;
	double		*s;
	double		tmin = Temp_min(seos);
	double		x, y;
	double		*T = seos->de_params.T;
	double		EE[2], pp[2];
	double		sinit,*et,*pt;
	double		r[2],rede[2],redp[2];
	double		*coords;
	double		svec[2], t2min;
	double		tph, rc;
	int		hT = nthyp - 1;
	int		n_pts_ds = hT + 2;
	int		icoords[MAXD];
	Locstate	state;
	COMPONENT	comp;
	int		status;
	int		i;

	DEBUG_ENTER(init_RT_interior_states)

	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);
	n4 = nr+nt+3*nr*nt+5;
	ntv = (int)(tbls[n4-1]);
	if (tbls[4] <= 0.) ntv = ntv-1;
	t2min = tmin;
	for(i = 0; i < nt-1; ++i)
	{
		if (tbls[4+nr+i]<= tmin && tbls[4+nr+i+1] >tmin)
			t2min = tbls[4+nr+i];
	}

	if (wave->sizest == 0)
	{
		DEBUG_LEAVE(init_RT_interior_states)
		return;
	}

	set_ses_intrp_flag_all(SESAME_RHO_TEMP);
	status = init_ses_hyp_soln_func(wave,front);
	if (DEBUG) (void) printf("Ended init_ses_hyp_soln_funct\n");

	if (status != GOOD_STEP)
	{
		screen("ERROR: init_ses_hyp_soln_func() failed\n");
		clean_up(ERROR);
	}

	intfc = front->interf;

	clear_utility_vectors(&seos->de_params);

	uni_array(&s,n_pts_ds,FLOAT);
	for (i = 0; i < n_pts_ds; ++i) s[i] = 0.0;

	for (ix = 0; ix < wave->rect_grid->gmax[0]; ++ix)
	{
		icoords[0] = ix;
		icoords[1] = 0;
		coords = Rect_coords(icoords,wave);
		x = coords[0];
		place = -1;
		init_PE_phase_grid(seos,x,phase_bound,&place);
		if (place == -1)
			set_sinit(&sinit,x,seos,phase_bound);
		else
		{
			tph = ses_rt_grid_from_temp(T[place],seos);
			set_S_on_pb(tph,t2min,seos,phase_bound,svec,ntv);
			get_rho_at_crit_pt(phase_bound,ntv,&rc);
			sinit = (ses_rt_rho_from_grid(x,seos) > rc) ? svec[1] :
								 svec[0];
		}
		sets(sinit,s,hT,&place,seos);

		for (iy = 0; iy < wave->rect_grid->gmax[1]; ++iy)
		{
			icoords[0] = ix;
			icoords[1] = iy;
			coords = Rect_coords(icoords,wave);
			x = coords[0];
			y = coords[1];
			comp = Rect_comp(icoords,wave);
			state = Rect_state(icoords,wave);

			if (is_exterior_comp(comp,intfc))
			{
				clear_state(intfc,state,sizest);
			}
			else
			{
			    sinit = set_S_state(s,y,T[place],place,seos);
			    lookspl2(y,seos,pp,EE);
			    if (comp == COMP_MIXED_PHASE)
			    {
			        get_temp_state(y,ntv,seos,phase_bound,r,
			        	rede,redp,&flag);
			        set_S_on_pb(y,t2min,seos,phase_bound,svec,ntv);
			        pp[0] = redp[0];
			        EE[0] = rede[1] +
			               (r[1]*r[0]*(rede[1] - rede[0])*
			               ((1.0/ses_rt_rho_from_grid(x,seos)) -
			               (1.0/r[1]))) / (r[0] - r[1]);
			        sinit = svec[1] +
			               (r[1]*r[0]*(svec[1] - svec[0])*
			               ((1.0/ses_rt_rho_from_grid(x,seos)) -
			               (1.0/r[1]))) / (r[0] - r[1]);
			    }
			    rho_T_phase_initializer(front,coords,comp,
			    	  state,seos,cold_curve,sinit,pp,EE);
			}
		}
	}
	free(s);

	uni_array(&et,wave->rect_grid->gmax[0]+8,FLOAT);
	uni_array(&pt,wave->rect_grid->gmax[0]+8,FLOAT);
	for (i = 0; i < wave->rect_grid->gmax[0]+8; ++i)
	{
		et[i] = 0.0;
		pt[i] = 0.0;
	}

	for(iy = 0; iy < wave->rect_grid->gmax[1]; ++iy)
	{
		init_PE_cross_grid(wave,front,seos,iy,pt,et);
	}
	free_these(2,et,pt);

	DEBUG_LEAVE(init_RT_interior_states)
}		/*end init_RT_interior_states*/

LOCAL void lookspl2(
	double		y,
	SESAME_EOS	*seos,
	double		*pp,
	double		*EE)
{
	int		phase = multiphase_eos(seos);
	double		*T = seos->de_params.T;
	double		*p = seos->de_params.P;
	double		*E = seos->de_params.E;
	double		*slopeP = seos->de_params.slopeP;
	double		*slopeE = seos->de_params.slopeE;
	double		ly = y;
	double		*tbls = seos->sestab.tbls;
	double		tmin = Temp_min(seos);
	double		tmax = Temp_max(seos );
	double		t1, t2, s1, s2, e1, e2, p1, p2;
	double		emin;
	int		nt,nr;
	int		i;

	DEBUG_ENTER(lookspl2)
	ly = ses_rt_temp_from_grid(ly,seos);
	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);

	emin = tbls[nr+nt + nr*nt +4];
	if (emin <= 0.0) emin = fabs(emin) + .01;
	else emin = 0.0;
	if (ly <= T[0])
	{
		EE[0] = (exp(E[0])) - emin;
		pp[0] = exp(p[0]);
		EE[1] = (EE[0] - exp(E[1]) +emin) / (T[0] - T[1]);
		pp[1] = (pp[0] - exp(p[1])) / (T[0] - T[1]);
		DEBUG_LEAVE(lookspl2)
		return;
	}


	for(i = 0; i < nt; ++i)
	{
		if ((phase == 1) &&
		(max(fabs(ly-T[i]),fabs(ly - T[i+1])) < TOLEPS*log(tmax/tmin)))
		{
			EE[0] = (exp(E[i])) - emin;
			pp[0] = exp(p[i]);
			EE[1] = (EE[0] - exp(E[i-1]) +emin) / (T[i] - T[i-1]);
			pp[1] = (pp[0] - exp(p[i-1])) / (T[i] - T[i-1]);
			DEBUG_LEAVE(lookspl2)
			return;
		}
		if (ly > T[i] && ly <= T[i+1] && T[i+1] > 0.)
		{
			t1 = ses_rt_grid_from_temp(T[i],seos);
			t2 = ses_rt_grid_from_temp(T[i+1],seos);
			e1 = E[i];
			e2 = E[i+1];
			s1 = slopeE[i];
			s2 = slopeE[i+1];
			ly = ses_rt_grid_from_temp(ly,seos);
			spline(t1,t2,e1,e2,s1,s2,ly,&EE[0],&EE[1]);
			EE[0] = (exp(EE[0])) - emin;
			EE[1] = ((EE[0] + emin) / exp(ly))*EE[1];
			p1 = p[i];
			p2 = p[i+1];
			s1 = slopeP[i];
			s2 = slopeP[i+1];
			spline(t1,t2,p1,p2,s1,s2,ly,&pp[0],&pp[1]);
			pp[0] = exp(pp[0]);
			pp[1] = (pp[0]/exp(ly))*pp[1];
			DEBUG_LEAVE(lookspl2)
			return;
		}
	}
	DEBUG_LEAVE(lookspl2)
}		/*end lookspl2*/


EXPORT void init_new_phase_bound(
	Wave		*wave,
	Front		*fr,
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	PHASE_BDRY	*phase_bound)
{
	INTERFACE	*intfc = fr->interf;
	NODE		*ns, *ne;
	double		*cold_Pph = phase_bound->cp;
	double		*cold_Eph = phase_bound->ce;
	double		*rho = phase_bound->rvar;
	double		*T = phase_bound->Tvar;
	double		*Pph = phase_bound->rp;
	double		*Eph = phase_bound->re;
	double		*R = seos->de_params.R;
	double		rho_grid, T_grid;
	double		*tbls = seos->sestab.tbls;
	double		rmin = Rho_min(seos);
	double		tmin = Temp_min(seos);
	double		rmax = Rho_max(seos);
	double		tmax = Temp_max(seos);
	double		*red_T, *red_R, *red_E, *red_P, *cld_P, *cld_E;
	double		m, r, p1, e1, es, ee, ps, pe, cpe, cps, cee, ces;
	double		temp;
	double		rhomax_grid, rhomin_grid, Tmax_grid, Tmin_grid;
	double		dT;
	double		rede[2], redp[2], rph[2];
	double		coords[MAXD];
	int		ntemp = (int)(tbls[3]);
	int		nr = (int)(tbls[2]);
	int		n4 = 5+nr+ntemp+3*nr*ntemp; /* Address of beginning 
						     * of phase curves */
	int		nt = (int)(tbls[n4-1]);
	int		nrhyp = fr->rect_grid->gmax[0];
	int		nthyp = fr->rect_grid->gmax[1];
	int		numph, i, j, ktot, ntot, n, i1, i2, iflag, flag;
	boolean		sav_intrp = interpolate_intfc_states(intfc);
	int		num_temp_pts = 2*(nrhyp+nthyp);
	int		max_num_pts = 4*nt*num_temp_pts;

	DEBUG_ENTER(init_new_phase_bound)
	rhomin_grid = ses_rt_grid_from_rho(rmin,seos);
	rhomax_grid = ses_rt_grid_from_rho(rmax,seos);
	Tmin_grid = ses_rt_grid_from_temp(tmin,seos);
	Tmax_grid = ses_rt_grid_from_temp(tmax,seos);

	uni_array(&red_T,max_num_pts,FLOAT);
	uni_array(&red_P,max_num_pts,FLOAT);
	uni_array(&cld_P,max_num_pts,FLOAT);
	uni_array(&red_E,max_num_pts,FLOAT);
	uni_array(&cld_E,max_num_pts,FLOAT);
	uni_array(&red_R,max_num_pts,FLOAT);
	for (i = 0; i < max_num_pts; ++i)
	{
		red_T[i] = 0.0;
		red_P[i] = 0.0;
		cld_P[i] = 0.0;
		red_E[i] = 0.0;
		cld_E[i] = 0.0;
		red_R[i] = 0.0;
	}

	dT = (Tmax_grid - Tmin_grid) / ((double) num_temp_pts);

	ktot = 0;
	i1 = 0;
	numph = nt;
	if (tbls[n4+2*nt] <= 0.0)
		numph = nt-1;
	if (numph == nt)
		ntot = 2*nt-1;
	else 
	{
		ntot = 2*nt - 3;
		i1 = 1;
	}


	n = -1;
	for(i = 0; i < ntot; ++i)
	{
		if (T[i] <= 0.0)
		{
			p1 = 0.0;
			e1 = 0.0;
		}
		else
		{
			get_cold_data(seos,cold_curve,rho[i],&cee, &cpe);
			p1 = (Pph[i]-cpe) / (rho[i]*T[i]);
			e1 = (Eph[i]-cee) / (T[i]);
		}
		red_P[ktot] = p1;
		red_E[ktot] = e1;
		red_R[ktot] = rho[i];
		red_T[ktot] = T[i];
		cld_E[ktot] = cee;
		cld_P[ktot] = cpe;
		++ktot;
		if (i == numph-1)
		{
			n = ktot;
			continue;
		}
		for (j = 0; j < num_temp_pts; ++j)
		{
			T_grid = (i < numph -1) ? Tmin_grid + j*dT :
					Tmin_grid +(num_temp_pts-j-1)*dT;
			temp = ses_rt_temp_from_grid(T_grid,seos);
			if ((temp > T[i] && temp < T[i+1] && i < numph -1)
			    || (temp < T[i] && temp > T[i+1] && i >= numph-1))
			{
				red_T[ktot] = temp;
				get_temp_state(T_grid,numph,seos,phase_bound,
					       rph,rede,redp,&flag);
				r = rph[0];
				if (i >= numph-1)
				{
				    red_R[ktot] = rph[1];
				    get_cold_data(seos,cold_curve,rph[1],
						  &cee, &cpe);
				    p1 = (redp[1]-cpe) /
					(red_R[ktot]*red_T[ktot]);
				    e1 = (rede[1]-cee) / (red_T[ktot]);
				    red_E[ktot] = e1;
				    red_P[ktot] = p1;
				}
				else
				{
				    red_R[ktot] = rph[0];
				    get_cold_data(seos,cold_curve,rph[0],
						  &cee, &cpe);
				    p1 = (redp[0]-cpe) /
					(red_R[ktot]*red_T[ktot]);
				    e1 = (rede[0]-cee) / (red_T[ktot]);
				    red_E[ktot] = e1;
				    red_P[ktot] = p1;
				}
				cld_E[ktot] = cee;
				cld_P[ktot] = cpe;
				++ktot;
			}
		}
	}

	if (n == -1)
	{
		screen("ERROR in init_new_phase_bound(), ");
		screen("n not set\n");
		(void) printf("ntot = %d, numph = %d\n",ntot,numph);
		clean_up(ERROR);
	}

	red_T[ktot] =T[ntot];
	red_R[ktot] = rho[ntot];
	cld_P[ktot] = cold_Pph[numph];
	cld_E[ktot] = cold_Eph[numph];
	if (T[ntot] <= 0.0)
	{
		red_P[ktot] = 0.0;
		red_E[ktot] = 0.0;
	}
	else
	{
		red_P[ktot] = (Pph[ntot]-cld_P[ktot]) / (rho[ntot]*T[ntot]);
		red_E[ktot] = (Eph[ntot]-cld_E[ktot]) / T[ntot];
	}
	if (DEBUG) (void) printf("Adding phase boundary to interface\n");

	/* Add phase boundary to interface*/





	if (DEBUG) (void) printf("Start and end states set\n");
	i1 = -1;
	i2 = -1;
	iflag = 0;

	for (i = 0; i < ktot; ++i)
	{
		if ((red_R[i] < rmin) && (red_R[i+1] >= rmin) &&
			(red_T[i] <= tmax) && (red_T[i] >= tmin))
		{
			m = (red_T[i+1] - red_T[i]) / (red_R[i+1]-red_R[i]);
			temp = m*(rmin - red_R[i+1]) + red_T[i+1];
			r = rmin;
			T_grid = ses_rt_grid_from_temp(temp,seos);
			m = (Eph[i+1] - Eph[i]) / (red_T[i+1]- red_T[i]);
			es = m*(temp - red_T[i+1]) +Eph[i+1];
			m = (Pph[i+1] - Pph[i]) / (red_T[i+1]- red_T[i]);
			ps = m*(temp - red_T[i+1]) +Pph[i+1];
			get_cold_data(seos,cold_curve,r,&ces,&cps);
			es = (es - ces)/temp;
			ps = (ps - cps)/(rmin*temp);
			coords[0] = rhomin_grid;
			coords[1] = T_grid;
			ns = make_node(Point(coords));
			node_type(ns) = PHASE_BDRY_NODE;
			i1 = i+1;
			if (red_R[i1] > rmax) i1 = 0;
		}


		/* Set states at pt */


		if ((red_R[i] <= rmax) && (red_R[i+1] > rmax) &&
			(red_T[i] <= tmax) && (red_T[i] >= tmin))
		{
			m = (red_T[i+1] - red_T[i]) / (red_R[i+1]-red_R[i]);
			temp = m*(rmax - red_R[i+1]) + red_T[i+1];
			m = (Eph[i+1] - Eph[i]) / (red_T[i+1]- red_T[i]);
			ee = m*(temp - red_T[i+1]) +Eph[i+1];
			m = (Pph[i+1] - Pph[i]) / (red_T[i+1]- red_T[i]);
			pe = m*(temp - red_T[i+1]) +Pph[i+1];
			get_cold_data(seos,cold_curve,r,&cee,&cpe);
			ee = (ee - cee)/temp;
			pe = (pe - cpe)/(rmax*temp);
			r = rmax;
			T_grid = ses_rt_grid_from_temp(temp,seos);
			coords[0] = rhomax_grid;
			coords[1] = T_grid;
			ne = make_node(Point(coords));
			node_type(ne) = PHASE_BDRY_NODE;
			i2 = i;
			if (red_R[i2] < rmin) i2 = 0;
		}

	}

	for (i = 0; i < n; ++i)
	{
		if ((red_T[i]<= tmin) && (red_T[i+1] > tmin) &&
			(red_R[i+1]<= rmax) && (red_R[i+1] >= rmin) &&
			(i1 == -1))
		{
			get_temp_state(Tmin_grid,numph,seos,phase_bound,rph,
				rede,redp,&flag);
			r = rph[0];
			es = rede[0];
			ps = redp[0];
			get_cold_data(seos,cold_curve,r,&ces,&cps);
			temp = tmin;
			es = (es - ces)/temp;
			ps = (ps - cps)/(r*temp);
			rho_grid = ses_rt_grid_from_rho(r,seos);
			coords[0] = rho_grid;
			coords[1] = Tmin_grid;
			ns = make_node(Point(coords));
			node_type(ns) = PHASE_BDRY_NODE;
			i1 = i+1;
		}

		if ((red_T[i]<= tmax) && (red_T[i+1] > tmax) &&
			(red_R[i]<= rmax) && (red_R[i] >= rmin))
		{
			m = (red_T[i+1] - red_T[i]) / (red_R[i+1]-red_R[i]);
			r = (tmax-red_T[i+1])/m + red_R[i+1];
			m = (red_T[i+1] - red_T[i]) / (red_E[i+1]-red_E[i]);
			ee = (tmax-red_T[i+1])/m + red_E[i+1];
			m = (red_T[i+1] - red_T[i]) / (red_P[i+1]-red_P[i]);
			pe = (tmax-red_T[i+1])/m + red_P[i+1];
			m = (cld_E[i+1] - cld_E[i]) / (red_T[i+1]-red_T[i]);
			cee = (tmax-red_T[i+1])*m + cld_E[i+1];
			m = (cld_P[i+1] - cld_P[i]) / (red_T[i+1]-red_T[i]);
			cpe = (tmax-red_T[i+1])*m + cld_P[i+1];
			temp = tmax;
			rho_grid = ses_rt_grid_from_rho(r,seos);
			coords[0] = rho_grid;
			coords[1] = Tmax_grid;
			ne = make_node(Point(coords));
			node_type(ne) = PHASE_BDRY_NODE;
			i2 = i;
		}
	}
	if (DEBUG) (void) printf( "setting vapor\n");
	if ((i1 != -1) && (i2 != -1))
	{
		set_vapor_dome(i1,i2,ee,es,pe,ps,cps,ces,cpe,cee,ns,ne,
			       red_R,red_P,red_E,red_T,cld_P,cld_E,fr,seos);
		i1 = -1;
		i2 = -1;
		iflag = 1;
	}
	for (i = n; i < ktot; ++i)
	{
		if ((red_T[i] >= tmin) && (red_T[i+1] < tmin) &&
			(red_R[i]<= rmax) && (red_R[i] >= rmin))
		{
			get_temp_state(Tmin_grid,numph,seos,phase_bound,
					rph,rede,redp,&flag);
			r = rph[1];
			ee = rede[1];
			pe = redp[1];
			get_cold_data(seos,cold_curve,r,&cee,&cpe);
			temp = tmin;
			rho_grid = ses_rt_grid_from_rho(r,seos);
			ee = (ee - cee)/temp;
			pe = (pe - cpe)/(r*temp);
			coords[0] = rho_grid;
			coords[1] = Tmin_grid;
			ne = make_node(Point(coords));
			node_type(ne) = PHASE_BDRY_NODE;
			i2 = i;
			if (fabs(red_T[i]- tmin) < TOLEPS) i2 = i-1;
		}

		if ((red_T[i]>= tmax) && (red_T[i+1] < tmax) &&
			(red_R[i+1]<= rmax) && (red_R[i+1] >= rmin))
		{
			m = (red_T[i+1] - red_T[i]) / (red_R[i+1]-red_R[i]);
			r = (tmax-red_T[i+1])/m + red_R[i+1];
			temp = tmax;
			m = (Eph[i+1] - Eph[i]) / (red_T[i+1]- red_T[i]);
			es = m*(temp - red_T[i+1]) +Eph[i+1];
			m = (Pph[i+1] - Pph[i]) / (red_T[i+1]- red_T[i]);
			ps = m*(temp - red_T[i+1]) +Pph[i+1];
			get_cold_data(seos,cold_curve,r,&ces,&cps);
			es = (es - ces)/temp;
			ps = (ps - cps)/(r*temp);
			rho_grid = ses_rt_grid_from_rho(r,seos);
			coords[0] = rho_grid;
			coords[1] = Tmax_grid;
			ns = make_node(Point(coords));
			node_type(ns) = PHASE_BDRY_NODE;
			i1 = i+1;
		}
	}

	if (DEBUG) (void) printf( "setting vapor\n");
	if ((i1 != -1) && (i2 != -1))
	{
		if (iflag == 0)
		{
			set_vapor_dome(i1,i2,ee,es,pe,ps,cps,ces,cpe,cee,ns,ne,
				red_R,red_P,red_E,red_T,cld_P,cld_E,fr,seos);
		}
		if (iflag == 1)
		{
			set_vapor_dome(i1,i2,ee,es,pe,ps,cps,ces,cpe,cee,
				ns,ne,red_R,red_P,red_E,red_T,cld_P,cld_E,
				fr,seos);
		}

	}

	insert_density_crossings(fr,wave,seos,cold_curve,phase_bound);

	if (multiphase_eos(seos) == YES)
	{
	    /* LISA,  make this a subroutine */
	    double		rhog, rhol;
	    double		t2min, svec[2];
	    double		*S;
	    double		lgrmin, lgrmax;
	    double		tph;

	    t2min = tmin;
	    for(i = 0; i < ntemp-1; ++i)
	    {
	    	if (tbls[4+nr+i]<= tmin && tbls[4+nr+i+1] >tmin)
	    		t2min = tbls[4+nr+i];
	    }

	    clear_utility_vectors(&seos->de_params);
	    uni_array(&S,nrhyp,FLOAT);

	    init_PE_spline(seos,phase_bound);
	    for(i = 0; i < seos->n_pts_ref_curve; ++i)
	    	seos->S_on_ref_curve[i] = 0.0;
	    if (Temp_min(seos) > 0.0 && phase_bound->place[0] > 0 )
	    {
	    	rhog = R[phase_bound->place[0]];
	    	lgrmin = log(rmin);
	    	phbnd(S,lgrmin,rhog,0.0,seos);
	    	for(i = 0; i < nrhyp; ++i)
	    		seos->S_on_ref_curve[i] = S[i];
	    }
	    set_entropy_on_ph_bnd(seos,phase_bound,fr);
	    if (phase_bound->place[1] > 0)
	    {
	    	rhol = R[phase_bound->place[1]];
	    	tph = ses_rt_grid_from_temp(t2min,seos);
	    	lgrmax = log(rmax);
	    	set_S_on_pb(tph,t2min,seos,phase_bound,svec,numph);
	    	phbnd(S,rhol,lgrmax,svec[1],seos);
	    	for(i = 0; i < nrhyp; ++i)
	    		seos->S_on_ref_curve[i+nrhyp] = S[i];
	    }
	    if (DEBUG)
	    {
	    	(void) printf("S on ref is\n");
	    	for(i = 0; i < 2*(Nrho_hyp(seos))+2; ++i)
	    	{
	    		(void) printf("seos->S_on_ref_curve[%d] = %g\n",
	    			i,seos->S_on_ref_curve[i]);
	    	}
	    }
	    free(S);
	}
	interpolate_intfc_states(intfc) = sav_intrp;

	free_these(6,red_E,red_P,red_R,red_T,cld_P,cld_E);
	DEBUG_LEAVE(init_new_phase_bound)
}		/*end init_new_phase_bound*/

/*
*			ses_phase_states()
*
*	Set the entropy, adb_gam, and gru_gam on the phase boundary.
*
*/

EXPORT	void ses_phase_states(
	Front		*fr,
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound)
{
	CURVE		**cur;
	BOND		*b;
	Locstate	state;
	INTERFACE	*intfc = fr->interf;
	double		*E, *P, *slopeP, *slopeE;
	int		nthyp = Ntemp_hyp(seos);
	int		hT = nthyp - 1;
	int		n_pts_ds = hT + 2;
	double		r, temp;
	double		*sinit;
	double		pp[2], EE[2];
	int		place,i;

	DEBUG_ENTER(ses_phase_states)

	clear_utility_vectors(&seos->de_params);
	P = seos->de_params.P;
	E = seos->de_params.E;
	slopeP = seos->de_params.slopeP;
	slopeE = seos->de_params.slopeE;

	uni_array(&sinit,n_pts_ds,FLOAT);
	for(i = 0; i < n_pts_ds; ++i)
		sinit[i] = 0.;

	for(cur = intfc->curves; cur && *cur; ++cur)
	{
	    if (wave_type(*cur) == PHASE_BOUNDARY)
	    {
		r = Coords((*cur)->start->posn)[0];

		/* Set grugam on phase boundary */

		place = -1;
		init_PE_phase_grid(seos,r,phase_bound,&place);
		temp = Coords((*cur)->start->posn)[1];
		if (place != -1)
		{
		    state = right_start_state(*cur);
		    ses_rt_adb_gam(state) = (exp(P[place])*slopeP[place]) /
			(exp(E[place])*slopeE[place]);
		    ses_rt_gru_gam(state) = ses_rt_adb_gam(state)*(1/exp(r));
		    state = left_start_state(*cur);
		    ses_rt_adb_gam(state) = (exp(P[place+1])*slopeP[place+1]) /
			(exp(E[place+1])*slopeE[place+1]);
		    ses_rt_gru_gam(state) = ses_rt_adb_gam(state)*(1/exp(r));
		}
		else
		{
		    state = right_start_state(*cur);
		    lookspl2(temp,seos,pp,EE);
		    ses_rt_adb_gam(state) = pp[1]/EE[1];
		    ses_rt_gru_gam(state) = (1.0/exp(r))*ses_rt_adb_gam(state);
		    state = left_start_state(*cur);
		    ses_rt_adb_gam(state) = pp[1]/EE[1];
		    ses_rt_gru_gam(state) = (1.0/exp(r))*ses_rt_adb_gam(state);
		}
		for (b = (*cur)->first; b != NULL; b = b->next)
		{
		    r = Coords(b->end)[0];

		    place = -1;
		    init_PE_phase_grid(seos,r,phase_bound,&place);
		    if (place != -1)
		    {
			state = right_state_at_point_on_curve(b->end,b,*cur);
			ses_rt_adb_gam(state) = (exp(P[place])*slopeP[place]) /
			    (exp(E[place])*slopeE[place]);
			ses_rt_gru_gam(state) =ses_rt_adb_gam(state)*(1/exp(r));
			state = left_state_at_point_on_curve(b->end,b,*cur);
			ses_rt_adb_gam(state) = (exp(P[place+1])*
			    slopeP[place+1]) / (exp(E[place+1])*
			    slopeE[place+1]);
			ses_rt_gru_gam(state) =ses_rt_adb_gam(state)*(1/exp(r));
		    }
		    else
		    {
			lookspl2(temp,seos,pp,EE);
			state = right_state_at_point_on_curve(b->end,b,*cur);
			ses_rt_adb_gam(state) = pp[1]/EE[1];
			ses_rt_gru_gam(state) = pp[1]/(EE[1]*exp(r));
			state = left_state_at_point_on_curve(b->end,b,*cur);
			ses_rt_adb_gam(state) = pp[1]/EE[1];
			ses_rt_gru_gam(state) = pp[1]/(EE[1]*exp(r));
		    }
		}
	    }
	}
	if (debugging("ses_print_hyp"))
		verbose_ses_show_intfc_states(fr->interf,SESAME_RHO_TEMP,seos);
	interpolate_intfc_states(intfc) = NO;
	free(sinit);
	DEBUG_LEAVE(ses_phase_states)
}		/*end set_phase_states*/



EXPORT	void set_boundary_states(
	Front		*fr,
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound,
	COLD_CURVE	*cold_curve)
{
	CURVE		**cur;
	size_t		sizest = fr->sizest;
	INTERFACE	*intfc = fr->interf;
	double		rmin = Rho_min(seos);
	double		tmin = Temp_min(seos);
	double		tmin_grid;
	double		*tbls = seos->sestab.tbls;
	double		rho, drho, T_grid;
	double		*Tvec, *Evec, *Pvec, *slopPvec, *slopEvec;
	double		*sinit;
	double           s;
	double		*T, *P, *E, *slopeE, *slopeP;
	double		r[8], rede[8], redp[8];
	double		pp[2], EE[2], coords[MAXD];
	double		t2min, svec[2], tph, rc;
	int		nrhyp = Nrho_hyp(seos);
	int		nthyp = Ntemp_hyp(seos);
	int		i, j, place, ix, flag;
	int		nt, nr, n4, ntv;
	int		hT = nthyp - 1;
	int		n_pts_ds = hT + 2;
	int		nrhypp1;
	static const int	NUM_EXTRA_PTS = 3; /*TOLERANCE*/

	DEBUG_ENTER(set_boundary_states)
	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);
	n4 = nr+nt+3*nr*nt+5;
	ntv = (int)(tbls[n4-1]);
	if (tbls[4] <= 0.0)
	    ntv = ntv-1;
	t2min = tmin;
	for(i = 0; i < nt-1; ++i)
	{
	    if (tbls[4+nr+i]<= tmin && tbls[4+nr+i+1] >tmin)
	    	t2min = tbls[4+nr+i];
	}


	tmin_grid = ses_rt_grid_from_temp(tmin,seos);


	clear_utility_vectors(&seos->de_params);
	E = seos->de_params.E;
	P = seos->de_params.P;
	T = seos->de_params.T;
	slopeE = seos->de_params.slopeE;
	slopeP = seos->de_params.slopeP;

	uni_array(&sinit,n_pts_ds,FLOAT);
	for (i = 0; i < n_pts_ds; ++i)
	    sinit[i] = 0.0;

	for (i = 0; i < 8; ++i)
	{
	    r[i] = 0.0;
	    rede[i] = 0.0;
	    redp[i] = 0.0;
	}

	for(cur = intfc->curves; cur && *cur; ++cur)
	{
	    if (!is_exterior_comp(negative_component((*cur)),intfc))
	    	continue;

		
	    if (node_type((*cur)->start) == ERROR)
	    	node_type((*cur)->start) = FIXED_NODE;
	    if (node_type((*cur)->end) == ERROR)
	    	node_type((*cur)->end) = FIXED_NODE;
	    wave_type(*cur) = DIRICHLET_BOUNDARY;
	    start_status(*cur) = FIXED;
	    end_status(*cur) = FIXED;
	    rho = Coords((*cur)->start->posn)[0];
	    T_grid = Coords((*cur)->start->posn)[1];
	    place = -1;
	    init_PE_phase_grid(seos,rho,phase_bound,&place);
	    if (place == -1)
	    	set_sinit(&s,rho,seos,phase_bound);
	    else
	    {
	    	tph = ses_rt_grid_from_temp(T[place],seos);
	    	set_S_on_pb(tph,t2min,seos,phase_bound,svec,ntv);
	    	get_rho_at_crit_pt(phase_bound,ntv,&rc);
	    	s = (ses_rt_rho_from_grid(rho,seos) > rc) ? svec[1] : svec[0];
	    }
	    sets(s,sinit,hT,&place,seos);
	    if (positive_component((*cur)) == COMP_MIXED_PHASE)
	    {
	    	get_temp_state(T_grid,ntv,seos,phase_bound,r,rede,redp,&flag);
		set_S_on_pb(T_grid,t2min,seos,phase_bound,svec,ntv);
		s = svec[1] + (r[1]*r[0]*(svec[1] - svec[0])*
			      ((1.0/ses_rt_rho_from_grid(rho,seos)) -
			      (1.0/r[1]))) / (r[0] - r[1]);
	    }
	    else
	    	s = set_S_state(sinit,T_grid,T[place],place,seos);
	    lookspl2(T_grid,seos,pp,EE);
	    coords[0] = Coords((*cur)->start->posn)[0];
	    coords[1] = Coords((*cur)->start->posn)[1];
	    set_rho_T_phase_state(left_start_state(*cur),coords,
	    		          seos,cold_curve,s,pp,EE);
	    ft_assign(right_start_state(*cur),left_start_state(*cur),sizest);
	    rho = Coords((*cur)->end->posn)[0];
	    T_grid = Coords((*cur)->end->posn)[1];
	    place = -1;
	    init_PE_phase_grid(seos,rho,phase_bound,&place);
	    if (place == -1)
	    	set_sinit(&s,rho,seos,phase_bound);
	    else
	    {
	    	tph = ses_rt_grid_from_temp(T[place],seos);
	    	set_S_on_pb(tph,t2min,seos,phase_bound,svec,ntv);
	    	get_rho_at_crit_pt(phase_bound,ntv,&rc);
	    	s = (ses_rt_rho_from_grid(rho,seos) > rc) ? svec[1] : svec[0];
	    }
	    sets(s,sinit,hT,&place,seos);
	    if (positive_component((*cur)) == COMP_MIXED_PHASE)
	    {
	    	get_temp_state(T_grid,ntv,seos,phase_bound,r,rede,redp,&flag);
		set_S_on_pb(T_grid,t2min,seos,phase_bound,svec,ntv);
		s = svec[1] + (r[1]*r[0]*(svec[1] - svec[0])*
			      ((1.0/ses_rt_rho_from_grid(rho,seos)) -
			      (1.0/r[1]))) / (r[0] - r[1]);
	    }
	    else
	    	s = set_S_state(sinit,T_grid,T[place],place,seos);
	    lookspl2(T_grid,seos,pp,EE);
	    coords[0] = Coords((*cur)->end->posn)[0];
	    coords[1] = Coords((*cur)->end->posn)[1];
	    set_rho_T_phase_state(left_end_state(*cur),coords,seos,
				  cold_curve,s,pp,EE);
	    ft_assign(right_end_state(*cur),left_end_state(*cur),sizest);
	}
	get_phase_temp_state(tmin_grid,fr,r,rede,redp,&flag,seos);

	/* Allocate temporary storage */

	nrhypp1 = nrhyp +1;

	uni_array(&Tvec,nrhypp1*(nt +2),FLOAT);
	uni_array(&Evec,nrhypp1*(nt +2),FLOAT);
	uni_array(&Pvec,nrhypp1*(nt +2),FLOAT);
	uni_array(&slopPvec,nrhypp1*(nt+2),FLOAT);
	uni_array(&slopEvec,nrhypp1*(nt +2),FLOAT);
	for (i = 0; i < nrhypp1*(nt+2); ++i)
	{
	    Tvec[i] = 0.0;
	    Evec[i] = 0.0;
	    Pvec[i] = 0.0;
	    slopPvec[i] = 0.0;
	    slopEvec[i] = 0.0;
	}

	/*Set states at boundary of interface */


	drho = fr->rect_grid->h[0];
	for (i = 0; i < nrhyp+1; ++i)
	{
	    rho = cell_edge(i,0,fr->rect_grid);
	    place = -1;
	    init_PE_phase_grid(seos,rho,phase_bound,&place);
	    if (place == -1)
	    	set_sinit(&s,rho,seos,phase_bound);
	    else
	    {
	    	tph = ses_rt_grid_from_temp(T[place],seos);
	    	set_S_on_pb(tph,t2min,seos,phase_bound,svec,ntv);
	    	get_rho_at_crit_pt(phase_bound,ntv,&rc);
	    	s = (ses_rt_rho_from_grid(rho,seos) > rc) ? svec[1] : svec[0];
	    }
	    sets(s,sinit,hT,&place,seos);
	    boundary_states(place,rho,sinit,
	    		fr,seos,cold_curve,phase_bound,i);
	    for(j = 0; j < nt+2; ++j)
	    {
	    	Tvec[i*(nt+2) + j] = T[j];
	    	Evec[i*(nt+2) + j] = E[j];
	    	Pvec[i*(nt+2) + j] = P[j];
	    	slopEvec[i*(nt+2) + j] = slopeE[j];
	    	slopPvec[i*(nt+2) + j] = slopeP[j];
	    }
	    if (flag == 2 && i < nrhyp)
	    {
	    	for(j = 1; j <= NUM_EXTRA_PTS; ++j)
	    	{
	    	    rho = ses_rt_grid_from_rho(rmin,seos) + i*drho +
	    	    	j*drho/(NUM_EXTRA_PTS + 1.0);
	    	    if (rho > r[1])
	    	    {
	    	    	place = -1;
	    	    	init_PE_phase_grid(seos,rho,phase_bound,&place);
	    	    	if (place == -1)
	    	    	{
	    	    	    set_sinit(&s,rho,seos,phase_bound);
	    	    	}
	    	    	else
	    	    	{
	    	    	    tph = ses_rt_grid_from_temp(T[place],seos);
	    	    	    set_S_on_pb(tph,t2min,seos,phase_bound,svec,ntv);
	    	    	    get_rho_at_crit_pt(phase_bound,ntv,&rc);
	    	    	    s = (ses_rt_rho_from_grid(rho,seos) > rc) ?
				svec[1] : svec[0];
	    	    	}
	    	    	sets(s,sinit,hT,&place,seos);
	    	    	ix = 5;
	    	    	boundary_states(place,rho,sinit,fr,seos,
	    	    		cold_curve,phase_bound,ix);
	    	    }
	    	}
	    }
	}
	if (debugging("ses_print_hyp"))
		verbose_ses_show_intfc_states(fr->interf,SESAME_RHO_TEMP,seos);
	interpolate_intfc_states(intfc) = NO;
	set_cross_states(fr,seos,phase_bound,Tvec,Evec,
			 Pvec,slopPvec,slopEvec);
	if (debugging("ses_print_hyp"))
		verbose_ses_show_intfc_states(fr->interf,SESAME_RHO_TEMP,seos);

	free(sinit);
	free_these(5,Tvec,Evec,Pvec,slopPvec,slopEvec);
	DEBUG_LEAVE(set_boundary_states)
}		/*end set_boundary_states*/



EXPORT	void cold_PE_spline(
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve)
{
	double		*tbls = seos->sestab.tbls;
	double		*x, *y, *slp, *z, *c, *d, *e, *dx, *dy;
	double		temp, rhol, rhog, el, eg, pv;
	double		rmax = Rho_max(seos);
	int		phase = multiphase_eos(seos);
	double		*cold_rho = cold_curve->rvar;
	double		*cold_P = cold_curve->cp;
	double		*cold_E = cold_curve->ce;
	double		*cold_Pslp = cold_curve->slcp;
	double		*cold_Eslp = cold_curve->slce;
	int		nr, nt, n4, ntv;
	int		i, j, k, jold, i1, jmax;

	DEBUG_ENTER(cold_PE_spline)
	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);

	if (phase == 1)
	{
		n4 = 4+nr+nt+3*nr*nt;
		ntv = (int)(tbls[n4]);
		rhog = tbls[n4+2*ntv+1];
		rhol = tbls[n4+3*ntv+1];
		pv = tbls[n4+1];
		temp = tbls[n4+ntv+1];
		eg = tbls[n4+4*ntv+1];
		el = tbls[n4+5*ntv +1];
	}

	for (i = 0; i < nr; ++i)
	{
		if (tbls[4+i] >= rmax)
		{
			jmax = i;
			break;
		}
	}


	uni_array(&slp,nr+4,FLOAT);
	uni_array(&x,nr+4,FLOAT);
	uni_array(&y,nr+4,FLOAT);
	uni_array(&z,nr+4,FLOAT);
	uni_array(&dx,nr+4,FLOAT);
	uni_array(&dy,nr+4,FLOAT);
	uni_array(&c,nr+4,FLOAT);
	uni_array(&d,nr+4,FLOAT);
	uni_array(&e,nr+4,FLOAT);
	for (i = 0; i < nr+4; ++i)
	{
		slp[i] = 0.0;
		x[i] = 0.0;
		y[i] = 0.0;
		z[i] = 0.0;
		dx[i] = 0.0;
		dy[i] = 0.0;
		c[i] = 0.0;
		d[i] = 0.0;
		e[i] = 0.0;
	}

	jold = 0;
	j = 0;
	i1 = 0;
	if (tbls[4] <= 0.0) i1 = 1;
	for(i = 0; i <= jmax; ++i)
	{
	    x[j] = ses_rt_grid_from_rho(tbls[4+i+i1],seos);
	    y[j] = tbls[4+nr+nt+i+i1];
	    z[j] = tbls[4+nr+nt+nr*nt+i+i1];
	    ++j;
	    if ((phase == 1) && (fabs(tbls[nr+4] - temp) < TOLEPS))
	    {
		if (i != nr-1-i1)
		{
		    if ((tbls[4+i+i1] <= rhol) && 
			(tbls[4+i+1+i1] > rhol))
		    {
			x[j] = ses_rt_grid_from_rho(rhol,seos);
			z[j] = el;
			y[j] = pv;
			++j;
			if (splcomp(x,y,slp,j,c,d,e,dx,dy) >= 0)
			{
			    screen("ERROR in cold_PE_spline(), ");
			    screen("splcomp failed\n");
			    clean_up (ERROR);
			}
			for(k = 0; k < j; ++k)
			{
			    cold_rho[k+jold] =
				ses_rt_rho_from_grid(x[k],seos);
			    cold_P[k+jold] = y[k];
			    cold_Pslp[k+jold] = slp[k];
			}
			if (splcomp(x,z,slp,j,c,d,e,dx,dy) >= 0)
			{
			    screen("ERROR in cold_PE_spline(), ");
			    screen("splcomp failed\n");
			    clean_up (ERROR);
			}
			for(k = 0; k < j; ++k)
			{
			    cold_E[k+jold] = z[k];
			    cold_Eslp[k+jold] = slp[k];
			}
			jold = j +jold;
			x[0] = ses_rt_grid_from_rho(rhol,seos);
			/* y[0] = log(pv +1.0); */
			z[0] = el;
			j = 1;
		    }
		    else if ((tbls[4+i+i1] < rhog) && 
			(tbls[4+i+1+i1] >= rhog) && (rhog > 0))
		    {
			x[j] = ses_rt_grid_from_rho(rhog,seos);
			z[j] = eg;
			y[j] = pv;
			++j;
			if (splcomp(x,y,slp,j,c,d,e,dx,dy) >= 0)
			{
			    screen("ERROR in cold_PE_spline(), ");
			    screen("splcomp failed\n");
			    clean_up (ERROR);
			}
			for(k = 0; k < j; ++k)
			{
			    cold_rho[k+jold] = 
				ses_rt_rho_from_grid(x[k],seos);
			    cold_P[k+jold] = y[k];
			    cold_Pslp[k+jold] = slp[k];
			}
			if (splcomp(x,z,slp,j,c,d,e,dx,dy) >= 0)
			{
			    screen("ERROR in cold_PE_spline(), ");
			    screen("splcomp failed\n");
			    clean_up (ERROR);
			}
			for(k = 0; k < j; ++k)
			{
			    cold_E[k+jold] = z[k];
			    cold_Eslp[k+jold] = slp[k];
			}
			jold = j +jold;
			x[0] = ses_rt_grid_from_rho(rhog,seos);
			y[0] = pv;
			z[0] = eg;
			j = 1;
		    }
		}
	    }
	}
	if (j >1)
	{

		if (splcomp(x,y,slp,j,c,d,e,dx,dy) >= 0)
		{
			screen("ERROR in cold_PE_spline(), ");
			screen("splcomp failed\n");
			clean_up (ERROR);
		}
		for(k = 0; k < j; ++k)
		{
			cold_rho[k+jold] = ses_rt_rho_from_grid(x[k],seos);
			cold_P[k+jold] = y[k];
			cold_Pslp[k+jold] = slp[k];
		}
		if (splcomp(x,z,slp,j,c,d,e,dx,dy) >= 0)
		{
			screen("ERROR in cold_PE_spline(), ");
			screen("splcomp failed\n");
			clean_up (ERROR);
		}
		for(k = 0; k < j; ++k)
		{
			cold_E[k+jold] = z[k];
			cold_Eslp[k+jold] = slp[k];
		}
		jold = j +jold;
	}
	if (DEBUG)
	{
		(void) printf("At end of cold_PE_spline(), nr = %d\n",nr);
		(void) printf("%-14s %-14s %-14s\n",
			      "COLD DENSITY","COLD PRESSURE","COLD PSLOPE");
		for (i = 0; i < nr; ++i)
			(void) printf("%-14g %-14g %-14g\n",
				cold_rho[i],cold_P[i],cold_Pslp[i]);
		(void) printf("\n");
		(void) printf("End of printout of cold quantities\n");
	}
	free_these(9,slp,x,y,z,dx,dy,c,d,e);
	DEBUG_LEAVE(cold_PE_spline)
}		/*end cold_PE_spline*/

EXPORT	void phase_spline(
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	PHASE_BDRY	*phase_bound)
{
	double		*tbls = seos->sestab.tbls;
	int		nr = (int)(tbls[2]);
	int		nt = (int)(tbls[3]);
	int		n4 = nr+nt+4+3*nr*nt;
	double		*cold_rho = cold_curve->rvar;
	double		*cold_P = cold_curve->cp;
	double		*cold_E = cold_curve->ce;
	double		*cold_Pslp = cold_curve->slcp;
	double		*cold_Eslp = cold_curve->slce;
	double		*cold_Pph = phase_bound->cp;
	double		*cold_Eph = phase_bound->ce;
	double		*rho = phase_bound->rvar;
	double		*T = phase_bound->Tvar;
	double		*red_Pph = phase_bound->rp;
	double		*red_Eph = phase_bound->re;
	int		numph;
	int		ntv = (int)(tbls[n4]);
	int		i, j, i1;
	double		r1, r2, y1, y2, s1, s2, r, y0, yp;

	DEBUG_ENTER(phase_spline)

	n4 = n4+1;
	numph = ntv;
	phase_bound->dim = 2*ntv;
	if (tbls[n4+2*ntv] <= 0.0)
	{
		numph = ntv - 1;
		phase_bound->dim = 2*ntv - 1;
	}
	for (i = 0; i < numph; ++i)
	{
		if (tbls[n4+2*ntv] <= 0.0)
		{
		/* Exclude the point rho = 0.0 from the phase boundary */

			i1 = i+1;
			rho[i] = tbls[n4+2*ntv+i1];
			red_Pph[i] = tbls[n4+i1];
			T[i] = tbls[n4+ntv+i1];
			red_Eph[i] = tbls[n4+4*ntv+i1];
		}
		else
		{
			rho[i] = tbls[n4+2*ntv+i];	
			rho[ntv+i] = tbls[n4+4*ntv-i-1];
			red_Pph[i] = tbls[n4+i];
			red_Pph[ntv+i] = tbls[n4+ntv-i-1];
			T[i] = tbls[n4+ntv+i];		
			T[ntv+i] = tbls[n4+2*ntv-i-1];
			red_Eph[i] = tbls[n4+4*ntv+i];
			red_Eph[ntv+i] = tbls[n4+6*ntv-i-1];
		}
	}

	if (tbls[n4+2*ntv] <= 0.0)
		for(i = 0; i < ntv; ++i)
		{
			rho[numph+i] = tbls[n4+4*ntv-i-1];
			red_Pph[numph+i] = tbls[n4+ntv-i-1];
			red_Eph[numph+i] = tbls[n4+6*ntv-i-1];
			T[numph+i] = tbls[n4+2*ntv-i-1];
		}

	if ((fabs(tbls[nr+4]-tbls[n4+ntv]) < TOLEPS) && (numph == ntv))
	{
		cold_Pph[0] = tbls[n4];
		cold_Pph[2*ntv-1] = tbls[n4];
		cold_Eph[0] = tbls[n4+4*ntv];
		cold_Eph[2*ntv-1] = tbls[n4+6*ntv-1];
	}


	for(i = 0; i < ntv; ++i)
	{
		if (tbls[n4+2*ntv+i] <= cold_rho[0])
		{
			cold_Pph[i] = 0.0;
			cold_Eph[i] = tbls[n4+4*ntv];
		}
		else for(j = 0; j < nr+3; ++j)
		{
			if ((tbls[n4+2*ntv+i] > cold_rho[j]) && 
				(tbls[n4+2*ntv+i] <= cold_rho[j+1]))
			{
				r1 = ses_rt_grid_from_rho(cold_rho[j],seos);
				r2 = ses_rt_grid_from_rho(cold_rho[j+1],seos);
				y1 = cold_P[j];
				y2 = cold_P[j+1];
				s1 = cold_Pslp[j];
				s2 = cold_Pslp[j+1];
				r = ses_rt_grid_from_rho(tbls[n4+2*ntv+i],seos);
				spline(r1,r2,y1,y2,s1,s2,r,&y0,&yp);
				cold_Pph[i] = y0;
				y1 = cold_E[j];
				y2 = cold_E[j+1];
				s1 = cold_Eslp[j];
				s2 = cold_Eslp[j+1];
				spline(r1,r2,y1,y2,s1,s2,r,&y0,&yp);
				cold_Eph[i] = y0;
			}
		}
	}
	for(i = 0; i < ntv; ++i)
	{
		for(j = 0; j < nr+3; ++j)
		{
			if ((tbls[n4+3*ntv+i] >= cold_rho[j]) && 
				(tbls[n4+3*ntv+i] < cold_rho[j+1]))
			{
				r1 = ses_rt_grid_from_rho(cold_rho[j],seos);
				r2 = ses_rt_grid_from_rho(cold_rho[j+1],seos);
				y1 = cold_E[j];
				y2 = cold_E[j+1];
				s1 = cold_Eslp[j];
				s2 = cold_Eslp[j+1];
				r = ses_rt_grid_from_rho(tbls[n4+3*ntv+i],seos);
				spline(r1,r2,y1,y2,s1,s2,r,&y0,&yp);
				cold_Eph[i+ntv] = y0;
				y1 = cold_P[j];
				y2 = cold_P[j+1];
				s1 = cold_Pslp[j];
				s2 = cold_Pslp[j+1];
				spline(r1,r2,y1,y2,s1,s2,r,&y0,&yp);
				cold_Pph[i+ntv] = y0;
			}
		}
	}
	DEBUG_LEAVE(phase_spline)
}		/*end phase_spline*/


/*
*				lookspl():
*
*	Determines y,dy/dx,z,dz/dx,by looking up points on an isotherm spline.
*	Note output is given in log variables.
*/

EXPORT	void lookspl(
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound,
	double		x,
	int		i,
	int		ii,
	int		itotal,
	int		jtotal,
	double		*y,
	double		*z,
	double		*dy,
	double		*dz)
{
	double		lx = x;
	double		*tmp1 = seos->slpp;
	double		*tmp2 = seos->slpe;
	double		*tbls = seos->sestab.tbls;
	double		rmax = Rho_max(seos);
	double		rmin = Rho_min(seos);
	int		phase = multiphase_eos(seos);
	int		k,j,kmax,nr,nt;
	double		s3,s4,y1,y2,s1,s2,t1,t2,z2,z1;
	double		dlogr,T;
	double		er,pr;
	double		temp,dtemp,rho;
	double		emin;

	DEBUG_ENTER(lookspl)

	y1 = 0.0;
	y2 = 0.0;
	z1 = 0.0;
	z2 = 0.0;
	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);
	emin = tbls[nr+nt + nr*nt +4];
	emin = (emin <= 0.0) ? fabs(emin) + .01 : 0.0;
	dlogr = log(rmax/rmin);
	rho = ses_rt_rho_from_grid(lx,seos);
	if (lx <= tmp1[ii*jtotal])
	{
		temp = tmp1[itotal*jtotal];
		*dy = temp;
		temp = tmp2[itotal*jtotal + ii*jtotal];
		*dz = temp;
		search(rho,seos,nr,4,&j);
		j = j+4;
		temp = log(tbls[nr+nt+j+(i-4)*nr]);
		*y = temp;
		temp = log(tbls[nr+nt+nr*nt+j+(i-4)*nr]+emin);
		*z = temp;
		DEBUG_LEAVE(lookspl)
		return;
	}
	for(k = 1; k < jtotal; ++k)
	{
		if ((tmp1[ii*jtotal+k] > lx) && (tmp1[ii*jtotal+k-1] <= lx))
		{
			t1 = tmp1[ii*jtotal+k-1];
			s1 = tmp1[itotal*jtotal + ii*jtotal+k-1];
			s3 = tmp2[itotal*jtotal + ii*jtotal+k-1];
			t2 = tmp1[ii*jtotal+k];
			s2 = tmp1[itotal*jtotal +ii*jtotal+k];
			s4 = tmp2[itotal*jtotal + ii*jtotal+k];
			kmax = k;
		}
	}
	if ((fabs(tmp1[ii*jtotal+jtotal-1] - lx) <= TOLEPS*dlogr) && 
		(fabs(tmp1[ii*jtotal+jtotal-1] - tmp1[ii*jtotal+jtotal-2]) >
								TOLEPS*dlogr))
	{
		kmax = jtotal -1;
		t1 = tmp1[ii*jtotal+jtotal-1];
		s1 = tmp1[itotal*jtotal + ii*jtotal+jtotal-1];
		search(rho,seos,nr,4,&j);
		j = j+4;
		temp = log(tbls[nr+nt+j+(i-4)*nr]);
		*y = temp;
		temp = log(tbls[nr+nt+nr*nt +j+(i-4)*nr]+emin);
		*z = temp;
		temp = s1;
		*dy = temp;
		temp = s2;
		*dz = temp;
		DEBUG_LEAVE(lookspl)
		return;
	}
	if (fabs(tmp1[ii*jtotal+kmax -1] -lx) <= TOLEPS*dlogr)
	{
		*dy = tmp1[itotal*jtotal + ii*jtotal+k-1];
		*dz = tmp2[itotal*jtotal + ii*jtotal+k-1];
		search(rho,seos,nr,4,&j);
		j = j+4;
		temp = log(tbls[nr+nt+j+(i-4)*nr]);
		*y = temp;
		temp = log(tbls[nr+nt+nr*nt+j+(i-4)*nr]+emin);
		*z = temp;
		DEBUG_LEAVE(lookspl)
		return;
	}

	if (kmax != jtotal -1)
	{
		if (phase == YES)
		{
			if (fabs(tmp1[ii*jtotal+kmax] - tmp1[ii*jtotal+kmax+1])
								<= TOLEPS )
			{
				if (!get_phase_state(tmp1[ii*jtotal+kmax+1],
						seos,phase_bound,&T,&er,&pr))
				{
					screen("ERROR in lookspl(), ");
					screen("get_phase_state() ");
					screen("returns NO\n");
					clean_up(ERROR);
				}
				T = exp(T);
				z2 = er + emin;
				y2 = pr;
				z2 = log(z2);
				y2 = log(y2);
			}
			if (kmax > 1)
			{
				if (fabs(tmp1[ii*jtotal+kmax-1] -
					tmp1[ii*jtotal+kmax-2]) <= TOLEPS)
				{
					if (!get_phase_state(
						tmp1[ii*jtotal+kmax-1],
						seos,phase_bound,&T,&er,&pr))
					{
						screen("ERROR in lookspl(), ");
						screen("get_phase_state() ");
						screen("returns NO\n");
						clean_up(ERROR);
					}
					T = exp(T);
					z1 = er + emin;
					y1 = pr;
					z1 = log(z1);
					y1 = log(y1);
				}
			}
		}
		if ((phase ==NO) ||
			(fabs(tmp1[ii*jtotal+kmax-1] -
					tmp1[ii*jtotal+kmax-2]) > TOLEPS) ||
			(kmax <= 1))
		{
			search(rho,seos,nr,4,&j);
			j = j+4;
			y1 = log(tbls[nr+nt+j+(i-4)*nr]);
			z1 = log(tbls[nr+nt+nr*nt+j+(i-4)*nr] + emin);
		}
		if ((phase == NO) ||
			(fabs(tmp1[ii*jtotal+kmax] - tmp1[ii*jtotal+kmax+1]) >
				TOLEPS))
		{
			search(rho,seos,nr,4,&j);
			j = j+4;
			y2 = log(tbls[nr+nt+j+1+(i-4)*nr]);
			z2 = log(tbls[nr+nt+nr*nt +j+1+(i-4)*nr] + emin);
		}
		spline(t1,t2,y1,y2,s1,s2,lx,&temp,&dtemp);
		*y = temp;
		*dy = dtemp;
		spline(t1,t2,z1,z2,s3,s4,lx,&temp,&dtemp);
		*z = temp;
		*dz = dtemp;
	}
	DEBUG_LEAVE(lookspl)
}		/*end lookspl*/

LOCAL	void get_cold_data(
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	double		rho,
	double		*cold_e,
	double		*cold_p)
{
	double		*tbls = seos->sestab.tbls;
	double		*cold_P = cold_curve->cp;
	double		*cold_rho = cold_curve->rvar;
	double		*cold_E = cold_curve->ce;
	double		*cold_Pslp = cold_curve->slcp;
	double		*cold_Eslp = cold_curve->slce;
	double		lrho = rho;
	double		tmp;
	double		r1, r2, s1, s2, y1, y2, dE, dP;
	int		j;
	int		nr = (int)(tbls[2]);
	int		nrp3;

	DEBUG_ENTER(get_cold_data)
	nrp3 = 3 + nr;
	if (lrho <= cold_rho[0] )
	{
		tmp = cold_E[0];
		*cold_e = tmp;
		tmp = cold_P[0];
		*cold_p = tmp;
		DEBUG_LEAVE(get_cold_data)
		return;
	}
	for(j = 0; j < nrp3; ++j)
	{
		if (fabs(lrho - cold_rho[j+1]) < TOLEPS)
		{
			*cold_e = cold_E[j+1];
			*cold_p = cold_P[j+1];
			DEBUG_LEAVE(get_cold_data)
			return;
		}
		if ((lrho >= cold_rho[j]) && (lrho <= cold_rho[j+1]))
		{
			r1 = ses_rt_grid_from_rho(cold_rho[j],seos);
			r2 = ses_rt_grid_from_rho(cold_rho[j+1],seos);
			y1 = cold_E[j];
			y2 = cold_E[j+1];
			s1 = cold_Eslp[j];
			s2 = cold_Eslp[j+1];
			lrho = ses_rt_grid_from_rho(lrho,seos);
			spline(r1,r2,y1,y2,s1,s2,lrho,&tmp,&dE);
			*cold_e = tmp;
			tmp = 0.0;
			y1 = cold_P[j];
			y2 = cold_P[j+1];
			s1 = cold_Pslp[j];
			s2 = cold_Pslp[j+1];
			spline(r1,r2,y1,y2,s1,s2,lrho,&tmp,&dP);
			*cold_p = tmp;
			DEBUG_LEAVE(get_cold_data)
			return;
		}
	}
	if (DEBUG)
	{
		(void) printf("Warning in get_cold_data. cold data not set\n");
		(void) printf("lrho = %g\n",lrho);
		for(j = 0; j < nrp3; ++j)
			(void) printf("cold_rho[%d] = %g\n",j,cold_rho[j]);
		(void) printf("cold_rho[%d] = %g\n",nrp3,cold_rho[nrp3]);
		(void) printf("\n\n");
	}
	DEBUG_LEAVE(get_cold_data)

}		/*end get_cold_data*/

LOCAL void init_PE_cross_grid(
	Wave		*wave,
	Front		*fr,
	SESAME_EOS	*seos,
	int		iy,
	double		*pt,
	double		*et)
{
	Locstate    	state;
	COMPONENT    	comp;
	int        	ix, i, j, k, itotal, jold, flag;
	int        	jmin, nr, nt;
	double        	*tbls = seos->sestab.tbls;
	double		**fstore;
	double        	E, re, ce, p, rp, cp, csq, adbgam;
	double        	x, y;
	double        	*slp, *xx, *yy, *zz, *c, *d, *e, *dx, *dy;
	double        	rho[8], rede[8], redp[8];
	double        	emin, a, b;
	double        	*coords;
	int        	icoords[MAXD];
	int        	ksav;

	DEBUG_ENTER(init_PE_cross_grid)
	itotal = wave->rect_grid->gmax[0] + 8;

	bi_array(&fstore,9,itotal,FLOAT);
	c   = fstore[0];
	d   = fstore[1];
	e   = fstore[2];
	yy  = fstore[3];
	zz  = fstore[4];
	dy  = fstore[5];
	xx  = fstore[6];
	dx  = fstore[7];
	slp = fstore[8];
	for (i = 0; i < itotal; ++i)
	{
	    c[i] = 0.0;
	    d[i] = 0.0;
	    e[i] = 0.0;
	    yy[i] = 0.0;
	    zz[i] = 0.0;
	    dy[i] = 0.0;
	    xx[i] = 0.0;
	    dx[i] = 0.0;
	    slp[i] = 0.0;
	}

	for (i = 0; i < 8; ++i)
	{
	    rho[i] = 0.0;
	    rede[i] = 0.0;
	    redp[i] = 0.0;
	}

	jold = 0;
	i = 0; /* LISA:rename this variable to distinguish from simple index*/
	icoords[0] = 0;
	icoords[1] = iy;
	coords = Rect_coords(icoords,wave);
	y = coords[1];
	get_phase_temp_state(y,fr,rho,rede,redp,&flag,seos);
	y = ses_rt_temp_from_grid(y,seos);
	nr = (int)(tbls[2]);
	nt = (int)(tbls[3]);
	emin = tbls[nr+nt + nr*nt +4];
	emin = (emin <= 0.0) ? fabs(emin) + 0.01 : 0.0;
	ksav = 0;
	for(ix = 0; ix < wave->rect_grid->gmax[0]; ++ix)
	{
	    icoords[0] = ix;
	    coords = Rect_coords(icoords,wave);
	    x = coords[0];
	    comp = Rect_comp(icoords,wave);
	    state = Rect_state(icoords,wave);

	    re = ses_rt_rede(state);
	    ce = ses_rt_colde(state);
	    rp = ses_rt_redp(state);
	    cp = ses_rt_coldp(state);
	    xx[i] = x;
	    x = ses_rt_rho_from_grid(x,seos);
	    E = y*re +ce + emin;
	    p = y*x*rp+cp;
	    zz[i] = log(E);
	    yy[i] = log(p);
	    x = ses_rt_grid_from_rho(x,seos);
	    icoords[0] = ix+1;
	    icoords[1] = iy;
	    coords = Rect_coords(icoords,wave);
	    ++i;
	    if (flag > 0)
	    {
	        for(j = 0; j < flag; ++j)
	        {
	            if ((ix < wave->rect_grid->gmax[0]-1) &&
	                (x < rho[j]) &&
	                (coords[0] >= rho[j]))
	            {
	                xx[i] = rho[j];
	                re = rede[j];
	                rp = redp[j];
	                x = ses_rt_rho_from_grid(rho[j],seos);
	                E = re + emin;
	                p = rp;
	                yy[i] = log(p);
	                zz[i] = log(E);
	                ++i;
	                a = (yy[1] - yy[0]) / (xx[1] - xx[0]);
	                b = (yy[i-1] - yy[i-2]) / (xx[i-1] - xx[i-2]);
	                if( splcomp2(xx,yy,slp,i,c,d,e,dx,dy,a,b) >= 0)
	                {
	                    screen("ERROR in init_PE_cross_grid(), ");
	                    screen("splcomp2 failed\n");
	                    clean_up (ERROR);
	                }
	                for(k = 0; k < i; ++k)
	                {
	                    pt[k] = slp[k];
	                }
	                a = (zz[1] - zz[0]) / (xx[1] - xx[0]);
	                b = (zz[i-1] - zz[i-2]) / (xx[i-1] - xx[i-2]);
	                if (splcomp2(xx,zz,slp,i,c,d,e,dx,dy,a,b) >= 0)
	                {
	                    screen("ERROR in init_PE_cross_grid(), ");
	                    screen("splcomp2 failed\n");
	                    clean_up (ERROR);
	                }
	                for(k = 0; k < i; ++k)
	                {
	                    et[k] = slp[k];
	                }
	                if (j>0)
	                {
	                    for(k = 1; k < i-1; ++k)
	                    {
	                        icoords[0] = ksav;
	                        icoords[1] = iy;
	                        coords = Rect_coords(icoords,wave);
	                        x = coords[0];
	                        y = coords[1];
	                        x = ses_rt_rho_from_grid(x,seos);
	                        y = ses_rt_temp_from_grid(y,seos);
	                        comp = Rect_comp(icoords,wave);
	                        state = Rect_state(icoords,wave);
	                        rp = ses_rt_redp(state);
	                        cp = ses_rt_coldp(state);
	                        re = ses_rt_rede(state);
	                        ce = ses_rt_colde(state);
	                        p = y*x*rp+cp;
	                        E = y*re + ce;
	                        adbgam = ses_rt_adb_gam(state);
	                        if (comp == COMP_MIXED_PHASE)
	                        {
	                            /* linear approx in mixed phase */
	                            pt[k] = 0;
	                            rho[0] = ses_rt_rho_from_grid(rho[0],seos);
	                            rho[1] = ses_rt_rho_from_grid(rho[1],seos);
	                            et[k] = rho[1]*rho[0]*(rede[1] - rede[0]) /
	                                (rho[0] - rho[1]);
	                            et[k] = -1.0*et[k]/(x*x);
	                            csq = (p/(x*x) - et[k])*adbgam;
	                            rho[0] = ses_rt_grid_from_rho(rho[0],seos);
	                            rho[1] = ses_rt_grid_from_rho(rho[1],seos);
	                        }
	                        else
	                            csq = (p/x)*pt[k] + (p/(x*x) -(E/x)*
	                                    et[k])*adbgam;

	                        if (csq < 0.0) 
	                        {
	                            csq = TOLEPS;
	                            if (DEBUG)
	                                (void) printf("Warning csq < 0\n");
	                        }
	                        ses_rt_adb_gam(state) = x*csq/p;
	                        ++ksav;
	                    }
	                }
	                else if (j == 0)
	                {
	                    for(k = 0; k < i-1; ++k)
	                    {
	                        icoords[0] = k;
	                        icoords[1] = iy;
	                        coords = Rect_coords(icoords,wave);
	                        x = coords[0];
	                        y = coords[1];
	                        x = ses_rt_rho_from_grid(x,seos);
	                        y = ses_rt_temp_from_grid(y,seos);
	                        comp = Rect_comp(icoords,wave);
	                        state = Rect_state(icoords,wave);
	                        rp = ses_rt_redp(state);
	                        cp = ses_rt_coldp(state);
	                        re = ses_rt_rede(state);
	                        ce = ses_rt_colde(state);
	                        p = y*x*rp+cp;
	                        E = y*re + ce;
	                        adbgam = ses_rt_adb_gam(state);
	                        if (comp == COMP_MIXED_PHASE)
	                        {
	                            pt[k] = 0;
	                            rho[0] = ses_rt_rho_from_grid(rho[0],seos);
	                            rho[1] = ses_rt_rho_from_grid(rho[1],seos);
	                            et[k] = rho[1]*rho[0]*(rede[1] - rede[0]) /
	                                    (rho[0] - rho[1]);
	                            et[k] = -1.0*et[k]/(x*x);
	                            csq = (p/(x*x) - et[k])*adbgam;
	                            rho[0] = ses_rt_grid_from_rho(rho[0],seos);
	                            rho[1] = ses_rt_grid_from_rho(rho[1],seos);
	                        }
	                        else
	                            csq = (p/x)*pt[k] + (p/(x*x) -(E/x)*
	                        et[k])*adbgam;
	                        if (csq < 0.0) 
	                        {
	                            csq = TOLEPS;
	                            if (DEBUG)
	                                (void) printf("Warning csq < 0\n");
	                        }
	                        ses_rt_adb_gam(state) = x*csq/p;
	                        ++ksav;
	                    }
	                }
	                xx[0] = rho[j];
	                yy[0] = log(p);
	                zz[0] = log(E+emin);
	                jold = jold +i;
	                i = 1;
	            }
	        }
	    }
	}
	if (i > 1)
	{
	    a = (yy[1] - yy[0]) / (xx[1] - xx[0]);
	    b = (yy[i-1] - yy[i-2]) / (xx[i-1] - xx[i-2]);
	    if (splcomp2(xx,yy,slp,i,c,d,e,dx,dy,a,b) >= 0)
	    {
	        screen("ERROR in init_PE_cross_grid(), ");
	        screen("splcomp2 failed\n");
	        clean_up (ERROR);
	    }
	    for(k = 0; k < i; ++k)
	    {
	        pt[k] = slp[k];
	    }
	    a = (zz[1] - zz[0]) / (xx[1] - xx[0]);
	    b = (zz[i-1] - zz[i-2]) / (xx[i-1] - xx[i-2]);
	    if (splcomp2(xx,zz,slp,i,c,d,e,dx,dy,a,b) >= 0)
	    {
	        screen("ERROR in init_PE_cross_grid(), ");
	        screen("splcomp2 failed\n");
	        clean_up (ERROR);
	    }
	    for(k = 0; k < i; ++k)
	    {
	        et[k] = slp[k];
	    }
	    if (ksav == 0)
	    {
	        jmin = 0;
	    }
	    else
	    {
	        jmin = 1;
	    }
	    for(k = jmin; k < i; ++k)
	    {
	        icoords[0] = ksav;
	        icoords[1] = iy;
	        coords = Rect_coords(icoords,wave);
	        x = coords[0];
	        y = coords[1];
	        x = ses_rt_rho_from_grid(x,seos);
	        y = ses_rt_temp_from_grid(y,seos);
	        comp = Rect_comp(icoords,wave);
	        state = Rect_state(icoords,wave);
	        rp = ses_rt_redp(state);
	        cp = ses_rt_coldp(state);
	        re = ses_rt_rede(state);
	        ce = ses_rt_colde(state);
	        p = y*x*rp+cp;
	        E = y*re + ce;
	        adbgam = ses_rt_adb_gam(state);
	        if (comp == COMP_MIXED_PHASE)
	        {
	            pt[k] = 0;
	            rho[0] = ses_rt_rho_from_grid(rho[0],seos);
	            rho[1] = ses_rt_rho_from_grid(rho[1],seos);
	            et[k] = rho[1]*rho[0]*(rede[1] - rede[0]) /
	                (rho[0] - rho[1]);
	            et[k] = -1.0*et[k]/(x*x);
	            csq = (p/(x*x) - et[k])*adbgam;
	            rho[0] = ses_rt_grid_from_rho(rho[0],seos);
	            rho[1] = ses_rt_grid_from_rho(rho[1],seos);
	        }
	        else
	            csq = (p/x)*pt[k] + (p/(x*x) -(E/x)*et[k])*adbgam;
	        if (csq < 0.0) 
	        {
	            csq = TOLEPS;
	            if (DEBUG)
	                (void) printf("Warning csq < 0\n");
	        }
	        ses_rt_adb_gam(state) = x*csq/p;
	        ++ksav;
	    }
	}
	free(fstore);
	free_these(3,rho,rede,redp);
	DEBUG_LEAVE(init_PE_cross_grid)
}	    /*end init_PE_cross_grid*/


LOCAL	void set_entropy_on_ph_bnd(
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound,
	Front		*fr)
{
	INTERFACE	*intfc = fr->interf;
	CURVE		**cur, *phsbdry;
	BOND		*b, *b1;
	Locstate	state;
	double		*tbls = seos->sestab.tbls;
	double		tmin = Temp_min(seos);
	double		t2min;
	double		*T, *P, *slopeP;
	double		*rho, *E, *slopeE, *slopeR;
	double		*c, *d, *e, *dx, *dy;
	double		sinit;
	double		svec[2];
	double		temp;
	double		r[2], er[2], pr[2];
	int		flag;
	int		n_pts_dome = phase_bound->n_pts_dome;
	int		i, j;
	int		nr = (int)(tbls[2]);
	int		nt = (int)(tbls[3]);
	int		n4 = nr+nt+3*nr*nt+5;
	int		ntv = (int)(tbls[n4-1]);

	DEBUG_ENTER(set_entropy_on_ph_bnd)
	if (tbls[4] <= 0.0)
		ntv = ntv-1;

	seos->de_params.pb_len = ntv;

	uni_array(&seos->de_params.pb_T,seos->de_params.pb_len,FLOAT);
	T = seos->de_params.pb_T;
	uni_array(&seos->de_params.pb_P,seos->de_params.pb_len,FLOAT);
	P = seos->de_params.pb_P;
	uni_array(&seos->de_params.pb_E,seos->de_params.pb_len,FLOAT);
	E = seos->de_params.pb_E;
	uni_array(&seos->de_params.pb_R,seos->de_params.pb_len,FLOAT);
	rho = seos->de_params.pb_R;
	uni_array(&seos->de_params.pb_slopeP,seos->de_params.pb_len,FLOAT);
	slopeP = seos->de_params.pb_slopeP;
	uni_array(&seos->de_params.pb_slopeE,seos->de_params.pb_len,FLOAT);
	slopeE = seos->de_params.pb_slopeE;
	uni_array(&seos->de_params.pb_slopeR,seos->de_params.pb_len,FLOAT);
	slopeR = seos->de_params.pb_slopeR;
	uni_array(&dx,ntv,FLOAT);
	uni_array(&dy,ntv,FLOAT);
	uni_array(&c,ntv,FLOAT);
	uni_array(&d,ntv,FLOAT);
	uni_array(&e,ntv,FLOAT);
	for (i = 0; i < ntv; ++i)
	{
		T[i] = 0.0;
		P[i] = 0.0;
		slopeP[i] = 0.0;
		E[i] = 0.0;
		rho[i] = 0.0;
		slopeE[i] = 0.0;
		slopeR[i] = 0.0;
		dx[i] = 0.0;
		dy[i] = 0.0;
		c[i] = 0.0;
		d[i] = 0.0;
		e[i] = 0.0;
	}

	for(cur = intfc->curves; cur && *cur; ++cur)
	{
		if (wave_type(*cur) == PHASE_BOUNDARY)
		{
			phsbdry = *cur;
			if (DEBUG)
			    print_curve(*cur);
		}
	}

	t2min = tmin;
	for(i = 0; i < nt-1; ++i)
	{
		if (tbls[4+nr+i]<= tmin && tbls[4+nr+i+1] >tmin)
			t2min = tbls[4+nr+i];
	}

	temp = ses_rt_grid_from_temp(t2min,seos);
	get_temp_state(temp,ntv,seos,phase_bound,r,er,pr,&flag);
	T[0] = log(t2min);
	P[0] = log(pr[0]);
	E[0] = log(er[0]);
	rho[0] = log(r[0]);
	
	j = 1;
	for (i = 0; i < ntv; ++i)
	{
		if (phase_bound->Tvar[i] > t2min)
		{
			T[j] = phase_bound->Tvar[i];
			T[j] = log(T[j]);
			P[j] = phase_bound->rp[i];
			E[j] = phase_bound->re[i];
			E[j] = log(E[j]);
			P[j] = log(P[j]);
			rho[j] = phase_bound->rvar[i];
			rho[j] = log(rho[j]);
			++j;
		}
	}

	if ((splcomp(T,E,slopeE,j,c,d,e,dx,dy) >= 0) ||
	    (splcomp(T,rho,slopeR,j,c,d,e,dx,dy) >= 0) ||
	    (splcomp(T,P,slopeP,j,c,d,e,dx,dy) >= 0))
	{
		screen("ERROR in set_entropy_on_ph_bnd(), ");
		screen("splcomp failed\n");
		clean_up(ERROR);
	}



	/* Compute S on vapor side of vapor dome */
	set_sinit(&sinit,rho[0],seos,phase_bound);
	setspb(sinit,phase_bound->S,n_pts_dome,T[j-1],seos);
	temp = Coords(phsbdry->first->start)[1];
	set_S_on_pb(temp,t2min,seos,phase_bound,svec,ntv);
	state = right_start_state(phsbdry);
	ses_rt_S(state) = svec[0];
	state = left_start_state(phsbdry);
	ses_rt_S(state) = svec[0];
	state = right_end_state(phsbdry);
	ses_rt_S(state) = svec[1];
	state = left_end_state(phsbdry);
	ses_rt_S(state) = svec[1];
	Entropy_min(seos) = min(svec[0],Entropy_min(seos));
	Entropy_max(seos) = max(svec[0],Entropy_max(seos));
	Entropy_min(seos) = min(svec[1],Entropy_min(seos));
	Entropy_max(seos) = max(svec[1],Entropy_max(seos));
	
	for (b = (phsbdry)->first; b != NULL; b = b->next)
	{
		if (Coords(b->end)[1] < Coords(b->start)[1]) continue;
		temp = Coords(b->end)[1];
		set_S_on_pb(temp,t2min,seos,phase_bound,svec,ntv);
		state = left_state_at_point_on_curve(b->end,b,phsbdry);
		ses_rt_S(state) = svec[0];
		state = right_state_at_point_on_curve(b->end,b,phsbdry);
		ses_rt_S(state) = svec[0];
		for (b1 = (phsbdry)->first; b1 != NULL; b1 = b1->next)
		{
			if ((fabs(temp - Coords(b1->end)[1]) < TOLEPS) &&
			    (b->end != b1->end))
			{
			    state = left_state_at_point_on_curve(b1->end,
								 b1,phsbdry);
			    ses_rt_S(state) = svec[1];
			    state = right_state_at_point_on_curve(b1->end,
								  b1,phsbdry);
			    ses_rt_S(state) = svec[1];
			}
		}
		Entropy_min(seos) = min(svec[0],Entropy_min(seos));
		Entropy_max(seos) = max(svec[0],Entropy_max(seos));
		Entropy_min(seos) = min(svec[1],Entropy_min(seos));
		Entropy_max(seos) = max(svec[1],Entropy_max(seos));
	}
	free_these(12,
		seos->de_params.pb_T,
		seos->de_params.pb_P,
		seos->de_params.pb_R,
		seos->de_params.pb_E,
		seos->de_params.pb_slopeP,
		seos->de_params.pb_slopeE,
		seos->de_params.pb_slopeR,
		dx,dy,c,d,e);
	seos->de_params.pb_T = NULL;
	seos->de_params.pb_P = NULL;
	seos->de_params.pb_R = NULL;
	seos->de_params.pb_E = NULL;
	seos->de_params.pb_slopeP = NULL;
	seos->de_params.pb_slopeE = NULL;
	seos->de_params.pb_slopeR = NULL;
	seos->de_params.pb_len = 0;
	DEBUG_LEAVE(set_entropy_on_ph_bnd)
}		/*end get_entropy_on_ph_bnd*/

LOCAL	void get_rho_at_crit_pt(
	PHASE_BDRY	*phase_bound,
	int		ntv,
	double		*rc)
{
	*rc = phase_bound->rvar[ntv];
}		/*end get_rho_at_crit_pt*/

LOCAL	void set_S_on_pb(
	double		temp,
	double		tmin,
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound,
	double		*svec,
	int		ntv)
{
	int		i;
	int		npts = phase_bound->n_pts_dome;
	double		ltemp = temp;
	double		dT;
	double		*Tvar = phase_bound->Tvar;
	double		*S = phase_bound->S;
	double		r[2], er[2], pr[2];
	double		t1, t2;
	int		flag;


	DEBUG_ENTER(set_S_on_pb)
	get_temp_state(ltemp,ntv,seos,phase_bound,r,er,pr,&flag);
	ltemp = ses_rt_temp_from_grid(ltemp,seos);
	ltemp = log(ltemp);
	dT = (log(phase_bound->Tvar[ntv-1]) - log(tmin))/npts;
	if (ltemp <= log(tmin))
		svec[0] = S[0];

	for (i = 0; i < npts; ++i)
	{
		t1 = log(tmin) +i*dT;
		t2 = log(tmin) + (i+1)*dT;
		if (ltemp >= t1 && ltemp <= t2)
		{
			svec[0] = (S[i+1] - S[i])*(ltemp - t1)/(t2 - t1)
				+S[i];
			break;
		}
	}

	ltemp = exp(ltemp);

	if (ltemp < Tvar[ntv])
	svec[1] = svec[0] + ((er[1] - er[0]) +pr[0]*(1/r[1] - 1/r[0]))/ltemp;
	else
	svec[1] = svec[0];
	DEBUG_LEAVE(set_S_on_pb)
}		/*end get_S_on_pb*/


LOCAL	void insert_density_crossings(
	Front		*fr,
	Wave		*wv,
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	PHASE_BDRY	*phase_bound)
{
	INTERFACE	*intfc = fr->interf;
	CURVE		**cur, *phsbdry;
	BOND		*b;
	POINT		*pt;
	RECT_GRID	*gr = wv->rect_grid;
	double		*tbls = seos->sestab.tbls;
	double		T, pr[2], r[2], er[2], E, S, ce, cp;
	double		rho, rhost, rhoend;
	double		coords[MAXD];
	int		ix, xmax = gr->gmax[0];
	int		numph, nr, nt, n4;
	int		j, flag;
	size_t		 sizest = fr->sizest;

	DEBUG_ENTER(insert_density_crossings)

	for(cur = intfc->curves; cur && *cur; ++cur)
	{
		if (wave_type(*cur) == PHASE_BOUNDARY)
		{
			phsbdry = *cur;
			if (DEBUG)
			{
				(void) printf("Phase boundary before ");
				(void) printf("insert_density_crossings()\n");
				print_curve(*cur);
			}
		}
	}
	if (phsbdry == NULL)
	{
		clean_up(ERROR);
	}

	nr = (int)(tbls[2]);
	nt = (int)(tbls[3]);
	n4 = nr+nt+3*nr*nt+5;
	numph = (int)(tbls[n4-1]);
	if (tbls[n4+2*numph] <= 0.0)
		numph = numph - 1;

	for (ix = 0; ix < xmax; ++ix)
	{
		rho = cell_center(ix,0,gr);

		/* Check for density crossing on phase boundary */
		if (!get_phase_hyp_state(rho,seos,&T,&E,&S))
			continue;

		get_temp_state(T,numph,seos,phase_bound,r,er,pr,&flag);
		for (j = 0; j < flag; ++j)
		{
		    for (b = (phsbdry)->first; b != NULL; b = b->next)
		    {
			rhost = Coords(b->start)[0];
			rhoend = Coords(b->end)[0];
			if (rhost < ses_rt_grid_from_rho(r[j],seos) && 
				ses_rt_grid_from_rho(r[j],seos) < rhoend )
			{
			    coords[0] = ses_rt_grid_from_rho(r[j],seos);
			    coords[1] = T;
			    pt = Point(coords);
			    get_cold_data(seos,cold_curve,r[j],&ce,&cp);
			    ses_rt_coldp(left_state(pt)) = cp;
			    ses_rt_colde(left_state(pt)) = ce;
			    ses_rt_redp(left_state(pt)) = (pr[j] - cp) /
				    (r[j]*ses_rt_temp_from_grid(T,seos));
			    ses_rt_rede(left_state(pt)) = (er[j] - ce) /
				    ses_rt_temp_from_grid(T,seos);
			    Pressure_min(seos) = min(pr[j],Pressure_min(seos));
			    Pressure_max(seos) = max(pr[j],Pressure_max(seos));
			    Energy_min(seos) = min(er[j],Energy_min(seos));
			    Energy_max(seos) = max(er[j],Energy_max(seos));
			    ft_assign(right_state(pt),left_state(pt),sizest);
			    if (insert_point_in_bond(pt,b,phsbdry) !=
				FUNCTION_SUCCEEDED)
	                    {
	                        screen("ERROR in insert_density_crossings(), "
		                       "insert_point_in_bond() failed\n");
	                        clean_up(ERROR);
	                    }
			    if (b == phsbdry->last) continue;
			    b = b->next;
			}
		    }
		}
	}
	if (DEBUG)
	{
	    (void) printf("Phase boundary after "
	                  "insert_density_crossings()\n");
	    print_curve(phsbdry);
	}

	DEBUG_LEAVE(insert_density_crossings)
}		/*end insert_density_crossings*/

LOCAL	void set_vapor_dome(
	int		i1,
	int		i2,
	double		ee,
	double		es,
	double		pe,
	double		ps,
	double		cps,
	double		ces,
	double		cpe,
	double		cee,
	NODE		*ns,
	NODE		*ne,
	double		*red_R,
	double		*red_P,
	double		*red_E,
	double		*red_T,
	double		*cld_P,
	double		*cld_E,
	Front		*fr,
	SESAME_EOS	*seos)
{
	INTERFACE	*intfc = fr->interf;
	CURVE		*curve;
	size_t		sizest = fr->sizest;
	POINT		*pt;
	int		i;
	double		coords[MAXD];
	double		tmin = Temp_min(seos);
	double		tmax = Temp_max(seos);

	DEBUG_ENTER(set_vapor_dome)
	if (i1 > i2)
	{
	    i = i1;
	    i1 = i2;
	    i2 = i;
	    curve = make_curve(COMP_PURE_PHASE,COMP_MIXED_PHASE,ne,ns);
	    wave_type(curve) = PHASE_BOUNDARY;
	    start_status(curve) = FIXED;
	    end_status(curve) = FIXED;
	    if (DEBUG)
		(void) printf( "Vapor dome set\n");
	    ses_rt_coldp(left_start_state(curve)) = cps;
	    ses_rt_colde(left_start_state(curve)) = ces;
	    ses_rt_redp(left_start_state(curve)) = ps;
	    ses_rt_rede(left_start_state(curve)) = es;
	    ft_assign(right_start_state(curve),left_start_state(curve),sizest);
	    ses_rt_coldp(left_end_state(curve)) = cpe;
	    ses_rt_colde(left_end_state(curve)) = cee;
	    ses_rt_redp(left_end_state(curve)) = pe;
	    ses_rt_rede(left_end_state(curve)) = ee;
	    ft_assign(right_end_state(curve),left_end_state(curve),sizest);
	    interpolate_intfc_states(intfc) = NO;
	}
	else 
	{
	    curve = make_curve(COMP_PURE_PHASE,COMP_MIXED_PHASE,ns,ne);
	    wave_type(curve) = PHASE_BOUNDARY;
	    start_status(curve) = INCIDENT;
	    end_status(curve) = INCIDENT;
	    if (DEBUG)
		(void) printf( "Vapor dome set\n");
	    ses_rt_coldp(left_start_state(curve)) = cps;
	    ses_rt_colde(left_start_state(curve)) = ces;
	    ses_rt_redp(left_start_state(curve)) = ps;
	    ses_rt_rede(left_start_state(curve)) = es;
	    ft_assign(right_start_state(curve),left_start_state(curve),sizest);
	    ses_rt_coldp(left_end_state(curve)) = cpe;
	    ses_rt_colde(left_end_state(curve)) = cee;
	    ses_rt_redp(left_end_state(curve)) = pe;
	    ses_rt_rede(left_end_state(curve)) = ee;
	    ft_assign(right_end_state(curve),left_end_state(curve),sizest);
	    interpolate_intfc_states(intfc) = NO;
	}


	if (i1 != i2)
	{
	    for(i = i1; i <= i2; ++i)
	    {
	        if ((fabs(red_T[i] - red_T[i-1]) > TOLEPS) &&
		     (fabs(red_T[i] - tmin) > TOLEPS*(tmax - tmin)))
	        {
	    	    coords[0] = ses_rt_grid_from_rho(red_R[i],seos);
	    	    coords[1] = ses_rt_grid_from_temp(red_T[i],seos);
	    	    pt = Point(coords);

	    	    /* Set states at pt */

	    	    ses_rt_coldp(left_state(pt)) = cld_P[i];
	    	    ses_rt_colde(left_state(pt)) = cld_E[i];
	    	    ses_rt_redp(left_state(pt)) = red_P[i];
	    	    ses_rt_rede(left_state(pt)) = red_E[i];
	    	    ft_assign(right_state(pt),left_state(pt),sizest);
	    	    if (insert_point_in_bond(pt,curve->last,curve) !=
			FUNCTION_SUCCEEDED)
	            {
	                screen("ERROR in set_vapor_dome(), "
		               "insert_point_in_bond() failed\n");
	                clean_up(ERROR);
	            }
		}
	    }
	}
	DEBUG_LEAVE(set_vapor_dome)
}		/*end set_vapor_dome*/

/*ARGSUSED*/
EXPORT	int get_phase_state(
	double		rho,
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound,
	double		*T,
	double		*er,
	double		*pr)
{
	double		*red_E = phase_bound->re;
	double		*red_P = phase_bound->rp;
	double		*r = phase_bound->rvar;
	double		*Temp = phase_bound->Tvar;
	int		dim = phase_bound->dim;
	double		slope;
	double		lrho;
	int		i;

	DEBUG_ENTER(get_phase_state)

	lrho = ses_rt_rho_from_grid(rho,seos);
	for(i = 0; i < dim; ++i)
	{
		if (Between(lrho,r[i],r[i+1]))
		{
			if (Temp[i+1] <= 0.0) break;
			lrho = ses_rt_grid_from_rho(lrho,seos);
			if (fabs(r[i] - r[i+1]) < TOLEPS)
			{
			    *T = ses_rt_grid_from_temp(Temp[i],seos);
			    *er = red_E[i];
			    *pr = red_P[i];
			}
			else
			{
			    slope = (lrho - ses_rt_grid_from_rho(r[i],seos)) /
			    	(ses_rt_grid_from_rho(r[i+1],seos) - 
					ses_rt_grid_from_rho(r[i],seos));
			    *T = (ses_rt_grid_from_temp(Temp[i+1],seos) -
			          ses_rt_grid_from_temp(Temp[i],seos))*slope + 
			    	ses_rt_grid_from_temp(Temp[i],seos);
			    slope = (*T - ses_rt_grid_from_temp(Temp[i],seos)) /
			    	(ses_rt_grid_from_temp(Temp[i+1],seos) -
			    	ses_rt_grid_from_temp(Temp[i],seos));
			    *pr = (log(red_P[i+1]) - log(red_P[i]))*slope +
			    	log(red_P[i]);
			    *pr = exp(*pr);
			    *er = (red_E[i+1] - red_E[i])*slope + red_E[i];

			}
			DEBUG_LEAVE(get_phase_state)
			return YES;
		}
	}
	*T = ERROR_FLOAT;
	*pr = ERROR_FLOAT;
	*er = ERROR_FLOAT;
	DEBUG_LEAVE(get_phase_state)
	return NO;
}		/*end get_phase_state*/

EXPORT	int get_phase_hyp_state(
	double		rho,
	SESAME_EOS	*seos,
	double		*T,
	double		*E,
	double		*S)
{
	CURVE		**cur;
	BOND		*b;
	Locstate	state1,state2;
	Front		*fr = seos->fr[SESAME_RHO_TEMP];
	INTERFACE	*intfc = fr->interf;
	double		temp1, colde, rede, tmp;
	double		ec, er;

	for(cur = intfc->curves; cur && *cur; ++cur)
	{
	    if (wave_type(*cur) == PHASE_BOUNDARY)
	    {
		for (b = (*cur)->first; b != NULL; b = b->next)
		{
		    if (Between(rho,Coords(b->start)[0],Coords(b->end)[0]))
		    {
			tmp = (rho-Coords(b->start)[0]) /
			      (Coords(b->start)[0]-Coords(b->end)[0]);
			*T = (Coords(b->start)[1] - Coords(b->end)[1])*tmp + 
				Coords(b->start)[1];
			temp1 = ses_rt_temp_from_grid(*T,seos);

			state1 = right_state_at_point_on_curve(b->start,b,*cur);
			state2 = right_state_at_point_on_curve(b->end,b,*cur);
			tmp = (*T - Coords(b->start)[1]) / 
			    (Coords(b->start)[1] - Coords(b->end)[1]);
			colde = ses_rt_colde(state1);
			rede = ses_rt_rede(state1);
			er = colde +
			   rede*ses_rt_temp_from_grid(Coords(b->start)[1],seos);
			colde = ses_rt_colde(state2);
			rede = ses_rt_rede(state2);
			ec = colde +
			   rede*ses_rt_temp_from_grid(Coords(b->end)[1],seos);
			*E = tmp*(ec -er) +er;
			rede = ses_rt_S(state1)-ses_rt_S(state2);
			*S = rede*tmp + ses_rt_S(state1);
			return YES;
		    }
		    else if (fabs(Coords(b->end)[0] -rho) < TOLEPS)
		    {
			*T = Coords(b->end)[1];
			temp1 = ses_rt_temp_from_grid(*T,seos);
			state1= right_state_at_point_on_curve(b->end,b,*cur);
			ec = ses_rt_colde(state1);
			er = ses_rt_rede(state1);
			*E = ec + (temp1)*er;
			*S = ses_rt_S(state1);
			return YES;
		    }
		}
	    }
	}
	return NO;
}		/*end get_phase_hyp_state*/

EXPORT	void get_phase_temp_state(
	double		T,
	Front		*fr,
	double		*r,
	double		*er,
	double		*pr,
	int		*flag,
	SESAME_EOS	*seos)
{
	CURVE		**cur;
	BOND		*b;
	Locstate	state1, state2;
	INTERFACE	*intfc = fr->interf;
	double		tmp;
	double		press, press1, press2;
	double		enrgy, enrgy1, enrgy2;
	int		j;

	DEBUG_ENTER(get_phase_temp_state)

	*flag = 0;
	for(cur = intfc->curves; cur && *cur; ++cur)
	{
	    if (wave_type(*cur) == PHASE_BOUNDARY)
	    {
	    	j = 0;
	    	for (b = (*cur)->first; b != NULL; b = b->next)
	    	{
	    	    if (!Between(T,Coords(b->start)[1],
	    			      Coords(b->end)[1]))
	    		    continue;

	       	    /*
	    	     *	Due to small pressures at the base
	    	     *	of the vapor domes,	interpolation
	    	     *	on the log of the pressure is preferred.
	    	     */
	    	    tmp = (T-Coords(b->start)[1]) /
	    		    (Coords(b->start)[1]-Coords(b->end)[1]);
	    	    r[j] = (Coords(b->start)[0] -
	    		    Coords(b->end)[0]) * tmp;
	    	    r[j] = r[j] + Coords(b->start)[0];
	    	    tmp = (T-Coords(b->start)[1]) /
	    		    (Coords(b->start)[1]-Coords(b->end)[1]);
	    	    state1 = right_state_at_point_on_curve(b->start,
	    						   b,*cur);
	    	    state2 = right_state_at_point_on_curve(b->end,
	    						   b,*cur);
	    	    press = ses_rt_redp(state1)*
	    		    ses_rt_rho_from_grid(Coords(b->start)[0],seos)*
	    		    ses_rt_temp_from_grid(Coords(b->start)[1],seos)
	    				    + ses_rt_coldp(state1);
	    	    press1 = ses_rt_redp(state2)*
	    		    ses_rt_rho_from_grid(Coords(b->end)[0],seos)*
	    		    ses_rt_temp_from_grid(Coords(b->end)[1],seos) +
	    				    ses_rt_coldp(state2);
	    	    press2 = (log(press) - log(press1)) /
	    		     (Coords(b->start)[1] - Coords(b->end)[1])*
	    		     (T - Coords(b->end)[1]) + log(press1);
	    	    pr[j] = exp(press2);

	    	    enrgy = ses_rt_rede(state1)*
	    		   ses_rt_temp_from_grid(Coords(b->start)[1],seos) +
	    			   ses_rt_colde(state1);
	    	    enrgy1 = ses_rt_rede(state2)*
	    		     ses_rt_temp_from_grid(Coords(b->end)[1],seos) +
	    			    ses_rt_colde(state2);
	    	    enrgy2 = (enrgy - enrgy1) /
	    		     (Coords(b->start)[1] - Coords(b->end)[1])*
	    		     (T - Coords(b->end)[1]) + enrgy1;
	    	    er[j] = enrgy2;

	    	    if (j > 0 && fabs(r[j] - r[j-1]) < TOLEPS) 
	    		    continue;
	    	    else
	    	    {
	    		++j;
	    		*flag = j;
	    		if (*flag >= 2)
	    		{
	    		    DEBUG_LEAVE(get_phase_temp_state)
	    		    return;
	    		}
	    	    }
	    	}
	    }
	}
	DEBUG_LEAVE(get_phase_temp_state)
}		/*end get_phase_temp_state*/

LOCAL	void search(
	double		x,
	SESAME_EOS	*seos,
	int		nsch,
	int		nstart,
	int		*j)
{
	double		*tbls = seos->sestab.tbls;
	int		i,tmp;

	DEBUG_ENTER(search)

	if (x <= tbls[nstart])
	{
		tmp = 0;
		*j = tmp;
		DEBUG_LEAVE(search)
		return;
	}
	if (x >= tbls[nstart+nsch-1])
	{
		tmp = nsch - 1;
		*j = tmp;
		DEBUG_LEAVE(search)
		return;
	}

	for(i = 0; i < nsch-1; ++i)
	{
		if (fabs(tbls[i+nstart]-x) <= TOLEPS*x)
		{
			tmp = i;
			*j = i;
			DEBUG_LEAVE(search)
			return;
		}
		if (tbls[i+nstart] <= x && tbls[i+1+nstart] > x)
		{
			tmp = i;
			*j = i;
			DEBUG_LEAVE(search)
			return;
		}
	}
	DEBUG_LEAVE(search)
}		/*end search*/

LOCAL	void set_sinit(
	double		*sinit,
	double		rho,
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound)
{
	double		lrho = rho;
	double		*sref = seos->S_on_ref_curve;
	double		*slpp = seos->slpp;
	double		rmin = Rho_min(seos);
	double		rmax = Rho_max(seos);
	int		nrhyp = Nrho_hyp(seos);
	int		*place = phase_bound->place;
	int		i;
	double		drho, rmin_grid, rph1, rph2;
	double		dlogr;
	double		r1, r2;


	DEBUG_ENTER(set_sinit)

	dlogr = log(rmax/rmin);
	rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	lrho = (rho - rmin_grid)/dlogr;
	if (place[1] != -1) 
	{
		rph2 = slpp[place[1]];
		drho = (log(rmax) - rph2)/(dlogr*(nrhyp - 1.0));
		if (fabs(rph2 - rho) < TOLEPS)
		{
			*sinit = sref[nrhyp];
			DEBUG_LEAVE(set_sinit)
			return;
		}
		if (fabs(log(rmax) - rho) < TOLEPS)
		{
			*sinit = sref[2*nrhyp-1];
			DEBUG_LEAVE(set_sinit)
			return;
		}
		for (i = 0; i < nrhyp-1; ++i)
		{
			r1 = i*drho + 
				(rph2 - log(rmin))/dlogr;
			r2 = (i+1)*drho +
				(rph2 - log(rmin))/dlogr;
			if (Between(lrho,r1,r2))
			{
				*sinit = ((sref[i+nrhyp+1] - sref[i+nrhyp]) /
					(r2- r1))*
					(lrho- r1) + sref[i+nrhyp];
				DEBUG_LEAVE(set_sinit)
				return;
			}
		}
	}
	if (place[0] != -1) 
	{
		rph1 = slpp[place[0]];
		if (fabs(rho - rph1) < TOLEPS)
		{
			*sinit = sref[nrhyp - 1];
			DEBUG_LEAVE(set_sinit)
			return;
		}
		drho = (rph1 - log(rmin))/(dlogr*(nrhyp - 1.0));
		for (i = 0; i < nrhyp-1; ++i)
		{
			r1 = i*drho; 
			r2 = (i+1)*drho;
			if (Between(lrho,r1,r2))
			{
			    *sinit = ((sref[i+1] - sref[i])/(r2- r1))*
					(lrho- r1) + sref[i];
			    DEBUG_LEAVE(set_sinit)
			    return;
			}
		}
	}
	if (DEBUG)
	{
		if (Between(rho,rph1,rph2))
		{
			(void) printf("Warning in set_sinit\n");
			(void) printf("rho out of range\n");
		}
		(void) printf("WARNING in set_sinit(), ");
		(void) printf("S not set\n");
		(void) printf("rho = %g, rph1 = %g, rph2 = %g\n",rho,rph1,rph2);
	}
	DEBUG_LEAVE(set_sinit)
}		/*end ses_sinit*/

LOCAL	double set_S_state(
	double		*Svec,
	double		temp,
	double		T,
	int		place,
	SESAME_EOS	*seos)
{
	double		S;
	double		tmin = Temp_min(seos);
	double		tmax = Temp_max(seos);
	double		*tbls = seos->sestab.tbls;
	double		*Tvec;
	double		ltemp = temp;
	double		t2min, dlogt, dT;
	int		nthyp = Ntemp_hyp(seos);
	int		i, ii;
	int		nr, nt;
	int		hT = nthyp - 1;
	int		n_pts_ds = hT + 2;

	DEBUG_ENTER(set_S_state)

	nr = (int)(tbls[2]);
	nt = (int)(tbls[3]);

	uni_array(&Tvec,n_pts_ds,FLOAT);
	for (i = 0; i < n_pts_ds; ++i) Tvec[i] = 0.0;

	ltemp = ses_rt_temp_from_grid(ltemp,seos);
	dT = 1.0/hT;

	t2min = tmin;
	if (place == -1)
	{
	    for(i = 0; i < nt-1; ++i)
	    {
	    	if (tbls[4+nr+i]<= tmin && tbls[4+nr+i+1] >tmin)
	    	    t2min = tbls[4+nr+i];
	    }
	}
	else
	    t2min = T;

	dlogt = log(tmax) - log(t2min);
	ltemp = (log(ltemp) - log(t2min))/dlogt;

	ii = 0;
	for(i = 0; i < n_pts_ds - 1; ++i)
	{
	    Tvec[ii] = i*dT;
	    ++ii;
	}

	if (ltemp <= 0.0)
	{
	    S = Svec[0];
	    free(Tvec);
	    DEBUG_LEAVE(set_S_state)
	    return S;
	}
	if (ltemp >= 1.0)
	{
	    S = Svec[ii-1];
	    free(Tvec);
	    DEBUG_LEAVE(set_S_state)
	    return S;
	}

	for(i = 0; i < ii-1; ++i)
	{
	    if ((Tvec[i] < ltemp) && (Tvec[i+1] >= ltemp))
	    {
	    	S = (Svec[i+1] - Svec[i])*(ltemp - Tvec[i]) /
	    	     (Tvec[i+1] - Tvec[i]) + Svec[i];
	    	free(Tvec);
	    	DEBUG_LEAVE(set_S_state)
	    	return S;
	    }
	}
	free(Tvec);
	DEBUG_LEAVE(set_S_state)
	return -HUGE_VAL;
}		/*end set_S_state*/


LOCAL	void boundary_states(
	int		place,
	double		rho,
	double		*sinit,
	Front		*fr,
	SESAME_EOS	*seos,
	COLD_CURVE	*cold_curve,
	PHASE_BDRY	*phase_bound,
	int		call)
{
	CURVE		**cur,*C;
	POINT		*pt;
	INTERFACE	*intfc = fr->interf;
	size_t		sizest = fr->sizest;
	double		*T = seos->de_params.T;
	double		coords[MAXD];
	double		rmin = Rho_min(seos);
	double		tmin = Temp_min(seos);
	double		rmax = Rho_max(seos);
	double		tmax = Temp_max(seos);
	double		rmax_grid,tmax_grid;
	double		rmin_grid,tmin_grid;
	double		*tbls = seos->sestab.tbls;
	double		dT, T_grid;
	double           s;
	double		tph;
	double		pp[2], EE[2];
	double		r[2], rede[2], redp[2];
	double		svec[2], t2min;
	int		flag;
	int		ntv, n4, nr, nt;
	int		i;
	int		nrhyp = Nrho_hyp(seos);
	int		nthyp = Ntemp_hyp(seos);

	DEBUG_ENTER(boundary_states)

	nr = (int)(tbls[2]);
	nt = (int)(tbls[3]);
	n4 = nr+nt+3*nr*nt+5;
	ntv = (int)(tbls[n4-1]);
	if (tbls[4] <= 0.) ntv = ntv-1;
	t2min = tmin;
	for(i = 0; i < nt-1; ++i)
	{
		if (tbls[4+nr+i]<= tmin && tbls[4+nr+i+1] >tmin)
			t2min = tbls[4+nr+i];
	}

	rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	tmin_grid = ses_rt_grid_from_temp(tmin,seos);
	dT = fr->rect_grid->h[1];

	flag = 0;
	if (call == 0)
	{
	    for(cur = intfc->curves; cur && *cur; ++cur)
	    {
	    	if (is_exterior_comp((negative_component(*cur)),intfc))
	    	{
	    	    if ((Coords((*cur)->start->posn)[0] <= rmax_grid) &&
	    		(Coords((*cur)->start->posn)[1] <= tmax_grid))
	    	    {
	    	        rmax_grid = Coords((*cur)->start->posn)[0];
	    	        tmax_grid = Coords((*cur)->start->posn)[1];
	    	        tph = Coords((*cur)->end->posn)[1];
	    	        C = *cur;
	    	    }
	    	}
	    }

	    dT = fr->rect_grid->h[1];


	    for (i = 1, T_grid = tmin_grid + dT; i < nthyp; ++i, T_grid += dT)
	    {
	        if (T_grid <= tph)
	        {
	            lookspl2(T_grid,seos,pp,EE);
	            if (positive_component(C) == COMP_PURE_PHASE)
	            {
	                s = set_S_state(sinit,T_grid,T[place],place,seos);
	            }
	            else
	            {
	                get_temp_state(T_grid,ntv,seos,phase_bound,r,
	            		       rede,redp,&flag);
	                pp[0] = redp[0];
	                EE[0] = rede[1] + (r[1]*r[0]*(rede[1] - rede[0])*
		                ((1.0/ses_rt_rho_from_grid(rho,seos)) -
		                (1.0/r[1])))/(r[0] - r[1]);
	                set_S_on_pb(T_grid,t2min,seos,
	                	    phase_bound,svec,ntv);
	                s = svec[1] + (r[1]*r[0]*(svec[1] - svec[0])*
		                ((1.0/ses_rt_rho_from_grid(rho,seos)) - 
		                (1.0/r[1])))/(r[0] - r[1]);
	                flag = 0;
	            }
	            coords[0] = rmin_grid;
	            coords[1] = T_grid;
	            pt = Point(coords);
	            set_rho_T_phase_state(left_state(pt),coords,
	    	    		      seos,cold_curve,s,pp,EE);
	            ft_assign(right_state(pt),left_state(pt),sizest);
	            if (!insert_point_in_bond(pt,C->last,C))
	            {
	                screen("ERROR in boundary_states(), "
	                       "insert_point_in_bond() failed\n");
	                clean_up(ERROR);
	            }
		}
	    }
	    rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	    tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	    rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	    tmin_grid = ses_rt_grid_from_temp(tmin,seos);
	    if ( tph < tmax_grid)
	    {
	    	for(cur = intfc->curves; cur && *cur; ++cur)
	    	{
	    	    if (!is_exterior_comp(negative_component((*cur)),intfc))
			continue;
	            if ((Coords((*cur)->start->posn)[0] <= rmin_grid)
					 &&
			(Coords((*cur)->start->posn)[1] <= tph)
					 &&
			(Coords((*cur)->start->posn)[1] > tmin_grid))
		    {
			C = *cur;
		    }
		}
		for (i = 1, T_grid=tmin_grid+dT; i < nthyp; ++i, T_grid += dT)
		{
		    if (T_grid < tph) continue;
		    s = set_S_state(sinit,T_grid,T[place],place,seos);
		    lookspl2(T_grid,seos,pp,EE);
		    coords[0] = rmin_grid;
		    coords[1] = T_grid;
		    pt = Point(coords);
		    set_rho_T_phase_state(left_state(pt),coords,
					  seos,cold_curve,s,pp,EE);
		    ft_assign(right_state(pt),left_state(pt),sizest);
		    if (!insert_point_in_bond(pt,C->last,C))
	            {
	                screen("ERROR in boundary_states(), "
		               "insert_point_in_bond() failed\n");
	                clean_up(ERROR);
	            }
		}
	    }
	}
	rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	tmin_grid = ses_rt_grid_from_temp(tmin,seos);
	if (call == (nrhyp))
	{
	    rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	    tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	    rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	    tmin_grid = ses_rt_grid_from_temp(tmin,seos);

	    for(cur = intfc->curves; cur && *cur; ++cur)
	    {
	    	if (!is_exterior_comp(negative_component((*cur)),intfc))
		    continue;

		if ((Coords((*cur)->end->posn)[0] >= rmin_grid) &&
		    (Coords((*cur)->end->posn)[1] <= tmax_grid))
		{
		    rmin_grid = Coords((*cur)->end->posn)[0];
		    tmax_grid = Coords((*cur)->end->posn)[1];
		    tph = Coords((*cur)->start->posn)[1];
		    C = *cur;
		}
	    }
	    for (i = 1, T_grid = tmin_grid + dT; i < nthyp; ++i, T_grid += dT)
	    {
	        if (T_grid < tph)
	        {
	    	    lookspl2(T_grid,seos,pp,EE);
	    	    if (positive_component(C) == COMP_PURE_PHASE)
		    {
			s = set_S_state(sinit,T_grid,T[place],place,seos);
		    }
		    else
		    {
		        get_temp_state(T_grid,ntv,seos,phase_bound,r,
		    		rede,redp,&flag);
		        pp[0] = redp[0];
		        EE[0] = rede[1] + (r[1]*r[0]*(rede[1] - rede[0])*
				    ((1.0/ses_rt_rho_from_grid(rho,seos)) -
				    (1.0/r[1])))/(r[0] - r[1]);
			set_S_on_pb(T_grid,t2min,seos,phase_bound,svec,ntv);
			s = svec[1] + (r[1]*r[0]*(svec[1] - svec[0])*
					((1.0/ses_rt_rho_from_grid(rho,seos)) -
					(1.0/r[1])))/(r[0] - r[1]);
			flag = 0;
		    }
		    coords[0] = rho;	coords[1] = T_grid;
		    pt = Point(coords);
		    set_rho_T_phase_state(left_state(pt),coords,
		    		      seos,cold_curve,s,pp,EE);
		    ft_assign(right_state(pt),left_state(pt),sizest);
		    if (!insert_point_in_bond(pt,C->first,C))
	            {
	                screen("ERROR in boundary_states(), "
		               "insert_point_in_bond() failed\n");
	                clean_up(ERROR);
	            }
		 }
	    }
	    rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	    tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	    rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	    tmin_grid = ses_rt_grid_from_temp(tmin,seos);

	    if ( tph < tmax_grid)
	    {
	        for(cur = intfc->curves; cur && *cur; ++cur)
	        {
	    	    if (!is_exterior_comp(negative_component((*cur)),intfc))
			continue;
		    if ((Coords((*cur)->start->posn)[0] >= rmax_grid) &&
		        (Coords((*cur)->end->posn)[0] >= rmax_grid) &&
		        (Coords((*cur)->start->posn)[1] >= tph) &&
			(Coords((*cur)->start->posn)[1] >= tmax_grid))
		    {
		        C = *cur;
		    }
	        }
		for (i = 1, T_grid=tmin_grid+dT; i < nthyp; ++i, T_grid += dT)
	        {
		    if (T_grid >= tph)
		    {
			s = set_S_state(sinit,T_grid,T[place],place,seos);
			lookspl2(T_grid,seos,pp,EE);
			coords[0] = rmax_grid;
			coords[1] = T_grid;
			pt = Point(coords);
			set_rho_T_phase_state(left_state(pt),
				       coords,seos,cold_curve,s,pp,EE);
			ft_assign(right_state(pt),left_state(pt),sizest);
			if (!insert_point_in_bond(pt,C->first,C))
	                {
	                    screen("ERROR in boundary_states(), "
		                   "insert_point_in_bond() failed\n");
	                    clean_up(ERROR);
	                }
		    }
		}
	    }
	}

	rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	tmin_grid = ses_rt_grid_from_temp(tmin,seos);

	if ((call > 0) && (call < (nrhyp)))
	{
	    /*Set states along tmin isotherm */
	    for(cur = intfc->curves; cur && *cur; ++cur)
	    {
	    	if (!is_exterior_comp(negative_component((*cur)),intfc))
	 	    continue;
	    	if ((Coords((*cur)->start->posn)[1] <= tmin_grid) &&
	    	    (Coords((*cur)->end->posn)[1] <= tmin_grid) &&
	    	    (Coords((*cur)->end->posn)[0] <= rho) &&
	    	    (Coords((*cur)->start->posn)[0] >= rho))
	    	{
	    	    C = *cur;
	    	}
	    }
	    lookspl2(tmin_grid,seos,pp,EE);
	    if (positive_component(C) == COMP_PURE_PHASE)
	    {
	        s = set_S_state(sinit,tmin_grid,T[place],place,seos);
	    }
	    else
	    {
	        get_temp_state(tmin_grid,ntv,seos,phase_bound,r,
				rede,redp,&flag);
	        pp[0] = redp[0];
	        EE[0] = rede[1] + (r[1]*r[0]*(rede[1] - rede[0])*
	    		((1.0/ses_rt_rho_from_grid(rho,seos)) -
	    		(1.0/r[1])))/(r[0] - r[1]);
	        set_S_on_pb(tmin_grid,t2min,seos,phase_bound,svec,ntv);
		s = svec[1] + (r[1]*r[0]*(svec[1] - svec[0])*
				((1.0/ses_rt_rho_from_grid(rho,seos)) -
				(1.0/r[1])))/(r[0] - r[1]);
	        flag = 0;
	    }
	    coords[0] = rho;
	    coords[1] = tmin_grid;
	    pt = Point(coords);
	    set_rho_T_phase_state(left_state(pt),coords,seos,cold_curve,
	    	          s,pp,EE);
	    ft_assign(right_state(pt),left_state(pt),sizest);
	    if (!insert_point_in_bond(pt,C->first,C))
	    {
	        screen("ERROR in boundary_states(), "
	               "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }

	    /*Set states along tmax isotherm */

	    for(cur = intfc->curves; cur && *cur; ++cur)
	    {
	    	if (!is_exterior_comp(negative_component((*cur)),intfc))
	    	    continue;
	    	if ((Coords((*cur)->start->posn)[1] >= tmax_grid) &&
	    	    (Coords((*cur)->end->posn)[1] >= tmax_grid) &&
	    	    Between(rho,Coords((*cur)->start->posn)[0],
	    	    	Coords((*cur)->end->posn)[0]))
	    	{
	    	    C = *cur;
	    	}
	    }
	    tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	    lookspl2(tmax_grid,seos,pp,EE);
	    if (positive_component(C) == COMP_PURE_PHASE)
	    {
	    	s = set_S_state(sinit,tmax_grid,T[place],place,seos);
	    }
	    else
	    {
	    	lookspl2(tmax_grid,seos,pp,EE);
	    	get_temp_state(tmax_grid,ntv,seos,phase_bound,r,
	    	           rede,redp,&flag);
	    	pp[0] = redp[0];
	    	EE[0] = rede[1] + (r[1]*r[0]*(rede[1] - rede[0])*
				((1.0/ses_rt_rho_from_grid(rho,seos)) -
				(1.0/r[1])))/(r[0] - r[1]);
	    	set_S_on_pb(tmax_grid,t2min,seos,phase_bound,
				svec,ntv);
	    	s = svec[1] + (r[1]*r[0]*(svec[1] - svec[0])*
				((1.0/ses_rt_rho_from_grid(rho,seos)) -
				(1.0/r[1])))/(r[0] - r[1]);
	    	flag = 0;
	    }
	    coords[0] = rho;
	    coords[1] = tmax_grid;
	    pt = Point(coords);
	    set_rho_T_phase_state(left_state(pt),coords,seos,cold_curve,
				      s,pp,EE);
	    ft_assign(right_state(pt),left_state(pt),sizest);
	    if (!insert_point_in_bond(pt,C->last,C))
	    {
	        screen("ERROR in boundary_states(), "
	               "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	}
	DEBUG_LEAVE(boundary_states)
}		/*end boundary_states*/

LOCAL	void	clear_utility_vectors(
	DEPARAMS	*de_params)
{
	int i;
	for (i = 0; i < de_params->len; ++i)
	{
		de_params->R[i] = 0.0;
		de_params->E[i] = 0.0;
		de_params->P[i] = 0.0;
		de_params->T[i] = 0.0;
		de_params->slopeE[i] = 0.0;
		de_params->slopeP[i] = 0.0;
	}
}		/*end clear_utility_vectors*/

#endif /* defined(TWOD) && defined(PHASE_CODE) */
