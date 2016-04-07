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
*			gsesspline.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Routines for splining for the SESASME EOS
*/

#if defined(TWOD) && defined(PHASE_CODE)
#define DEBUG_STRING    "ses_hyp"
#include <geos/sesame.h>

LOCAL const double TOLEPS = 1.0e-05; /*TOLERANCE*/

	/* LOCAL Function Declarations */
LOCAL	int	is_node_of_phase_boundary(NODE*);
LOCAL	void	adb_gam_state(Front*,SESAME_EOS*,double*,double*,double*,
			      double,double*,double*);
LOCAL	void	get_de_and_dp(double,Front*,double*,double*,double,double,double,
			      double*,double*,double*,double*,double*,int,double*,
			      SESAME_EOS*);
LOCAL	void	set_remaining_states(NODE*,Locstate,Locstate,Front*);
LOCAL	void	set_vectors(Front*,double*,double*,double*,SESAME_EOS*,
			    double,double*,double*,double*,double*,double*);

EXPORT	void init_PE_spline(
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound)
{
	double		*R = seos->de_params.R;
	double		*P = seos->de_params.P;
	double		*E = seos->de_params.E;
	double		*slopeP = seos->de_params.slopeP;
	double		*slopeE = seos->de_params.slopeE;
	int		phase = multiphase_eos(seos);
	double		*tbls = seos->sestab.tbls;
	int		*place = phase_bound->place;
	int		dim;
	double		*slpp = seos->slpp;
	double		*slpe = seos->slpe;
	double		rmin = Rho_min(seos);
	double		rmax = Rho_max(seos);
	double		tmin = Temp_min(seos);
	double		tmax = Temp_max(seos);
	double		*slp, *x, *y, *z, *dx, *dy, *c, *d, *e;
	int		nt, i, j, nr, ntv, imin, imax, jmin, jmax, jj, k, ii;
	int		n4, jold, itotal, jtotal;
	int		flag;
	double		rl, rg, pv, temp;
	double		el, eg, tv;
	double		rede[2], redp[2], rho[2];
	double		emin, a, b;

	DEBUG_ENTER(init_PE_spline)

	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);
	emin = tbls[nr+nt + nr*nt +4];
	emin = (emin <= 0.0) ? fabs(emin) + .01 : 0.0;
	n4 = nr+nt+3*nr*nt+5;
	ntv = (int)(tbls[n4-1]);
	dim = (tbls[4] <= 0.0) ? ntv-1 : ntv;
	place[1] = -1;
	place[0] = -1;

	for (i = 0; i < nt-1; i++)
	{
		if (tbls[4+nr+i] <= tmin && tbls[nr+4+i+1] > tmin)
		{
			imin = i;
		}
		if (tbls[4+nr+i] < tmax && tbls[4+nr+i+1] >= tmax)
		{
			imax = i+1;
		}
	}

	itotal = imax - imin + 9;
	for (i = 0; i < nr-1; i++)
	{
		if (tbls[4+i] <= rmin && tbls[4+i+1] > rmin)
		{
			jmin = i;
		}
		if (tbls[4+i] < rmax && tbls[4+i+1] >= rmax)
		{
			jmax = i+1;
		}
	}
	for (i = 0; i < nr; i++)
	{
		if (tbls[4+i] >= rmax)
		{
			jmax = i;
			break;
		}
	}

	jtotal = jmax - jmin + 9;

	uni_array(&slp,jtotal,FLOAT);
	uni_array(&x,jtotal,FLOAT);
	uni_array(&y,jtotal,FLOAT);
	uni_array(&z,jtotal,FLOAT);
	uni_array(&dx,jtotal,FLOAT);
	uni_array(&dy,jtotal,FLOAT);
	uni_array(&c,jtotal,FLOAT);
	uni_array(&d,jtotal,FLOAT);
	uni_array(&e,jtotal,FLOAT);
	for (i = 0; i < jtotal; i++)
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

	/* Set up coarse grid isotherm splining.*/

	jj = 0;
	if (tbls[imin+4+nr] <= 0.0)imin = imin +1;
	tv = tbls[n4 + 2*ntv - 1];

	for(i = imin+4; i <= imax+4; i++)
	{
	    temp = tbls[i+nr];

	    jold = 0;
	    temp = ses_rt_grid_from_temp(temp,seos);
	    flag = 0;
	    get_temp_state(temp,dim,seos,phase_bound,rho,rede,redp,&flag);
	    temp = ses_rt_temp_from_grid(temp,seos);
	    if (flag == 2)
	    {
		rg = rho[0];
		rl = rho[1];
		eg = rede[0];
		el = rede[1];
		pv = redp[0];
	    }
	    j = jmin + 4;
	    ii = i-imin-4;
	    for(j = jmin+4; j <= jmax +4; j++)
	    {
		x[jj] = tbls[j];
		y[jj] = tbls[nr+nt+j+(i-4)*nr];
		z[jj] = tbls[nr+nt+nr*nt +j+(i-4)*nr];
		x[jj] = log(x[jj]);
		y[jj] = log(y[jj]);
		z[jj] = log(z[jj] + emin);
		if (fabs(x[1] - x[0]) < TOLEPS) jj = 0;
		jj++;
		if (!phase || temp > tv) continue;
		if (tbls[j+1] >= rg && tbls[j] < rg)
		{
		    x[jj] = log(rg);
		    y[jj] = log(pv);
		    z[jj] = log(eg + emin);
		    a = (y[1] - y[0])/(x[1] - x[0]);
		    b = (y[jj] - y[jj -1])/(x[jj] - x[jj-1]);
		    jj++;
		    if (splcomp2(x,y,slp,jj,c,d,e,dx,dy,a,b) >= 0)
		    {
			screen("ERROR in init_PE_spline(), ");
			screen("splcomp2 failed\n");
			clean_up(ERROR);
		    }
		    for(k = 0; k < jj; k++)
		    {
			slpp[itotal*jtotal+ii*jtotal+k] = slp[k];
			slpp[ii*jtotal + k] = x[k];
			if (ii == 0)
			{
			    R[k] = x[k];
			    P[k] = y[k];
			    slopeP[k] = slp[k];
			}
		    }
		    jj--;
		    a = (z[1] - z[0])/(x[1] - x[0]);
		    b = (z[jj] - z[jj -1])/(x[jj] - x[jj-1]);
		    jj++;
		    if (splcomp2(x,z,slp,jj,c,d,e,dx,dy,a,b) >= 0)
		    {
			screen("ERROR in init_PE_spline(), ");
			screen("splcomp2 failed\n");
			clean_up(ERROR);
		    }
		    for(k = 0; k < jj; k++)
		    {
			slpe[itotal*jtotal+ii*jtotal+k] = slp[k];
			slpe[ii*jtotal + k] = x[k];
			if (ii == 0)
			{
			    E[k] = z[k];
			    slopeE[k] = slp[k];
			}
		    }
		    jold = jj;
		    if (ii == 0) place[0] = jold;
		    jj = 0;
		    x[0] = log(rg);
		    y[0] = log(pv);
		    z[0] = log(eg+emin);
		    jj++;
		}
		if ((tbls[j+1] >= rl) && (tbls[j] < rl) && (j < (jmax+4)))
		{
		    x[jj] = log(rl);
		    y[jj] = log(pv);
		    z[jj] = log(el +emin);
		    a = (y[1] - y[0])/(x[1] - x[0]);
		    b = (y[jj] - y[jj -1])/(x[jj] - x[jj-1]);
		    jj++;
		    if (splcomp2(x,y,slp,jj,c,d,e,dx,dy,a,b) >= 0)
		    {
			screen("ERROR in init_PE_spline(), ");
			screen("splcomp2 failed\n");
			clean_up(ERROR);
		    }
		    for(k = 0; k < jj; k++)
		    {
			slpp[ii*jtotal + (k+jold)] = x[k];
			slpp[itotal*jtotal+ii*jtotal+(k+jold)] = slp[k];
			if (ii == 0)
			{
			    R[k+jold] = x[k];
			    P[k+jold] = y[k];
			    slopeP[k+jold] = slp[k];
			}
		    }
		    jj--;
		    a = (z[1] - z[0])/(x[1] - x[0]);
		    b = (z[jj] - z[jj -1])/(x[jj] - x[jj-1]);
		    jj++;
		    if (splcomp2(x,z,slp,jj,c,d,e,dx,dy,a,b) >= 0)
		    {
			screen("ERROR in init_PE_spline(), ");
			screen("splcomp2 failed\n");
			clean_up(ERROR);
		    }
		    for(k = 0; k < jj; k++)
		    {
			slpe[ii*jtotal + (k+jold)] = x[k];
			slpe[itotal*jtotal+ii*jtotal+(k+jold)] = slp[k];
			if (ii == 0)
			{
			    E[k+jold] = z[k];
			    slopeE[k+jold] = slp[k];
			}
		    }
		    jold = jj+jold;
		    if (ii == 0) place[1] = jold;
		    jj = 0;
		    x[0] = log(rl);
		    y[0] = log(pv);
		    z[0] = log(el + emin);
		    jj++;
		}
	    }

	    if (jj > 1)
	    {
		jj--;
		a = (y[1] - y[0])/(x[1] - x[0]);
		b = (y[jj] - y[jj -1])/(x[jj] - x[jj-1]);
		jj++;
		if (splcomp2(x,y,slp,jj,c,d,e,dx,dy,a,b) >= 0)
		{
		    screen("ERROR in init_PE_spline(), ");
		    screen("splcomp2 failed\n");
		    clean_up(ERROR);
		}
		for(k = 0; k < jj; k++)
		{
		    slpp[ii*jtotal + k+jold] = x[k];
		    slpp[itotal*jtotal+ii*jtotal + k+jold] = slp[k];
		    if (ii == 0)
		    {
			R[k+jold] = x[k];
			P[k+jold] = y[k];
			slopeP[k+jold] = slp[k];
		    }
		}
		jj--;
		a = (z[1] - z[0])/(x[1] - x[0]);
		b = (z[jj] - z[jj -1])/(x[jj] - x[jj-1]);
		jj++;
		if (splcomp2(x,z,slp,jj,c,d,e,dx,dy,a,b) >= 0)
		{
		    screen("ERROR in init_PE_spline(), ");
		    screen("splcomp2 failed\n");
		    clean_up(ERROR);
		}
		for(k = 0; k < jj; k++)
		{
		    slpe[ii*jtotal + (k+jold)] = x[k];
		    slpe[itotal*jtotal+ii*jtotal + (k+jold)] = slp[k];
		    if (ii == 0)
		    {
			E[k+jold] = z[k];
			slopeE[k+jold] = slp[k];
		    }
		}
		jold = 0;
		jj = 0;
	    }
	}
	if (DEBUG)
	{
		(void) printf("In init_PE_spline()\n");
		(void) printf("itotal = %d, jtotal = %d\n",itotal,jtotal);
		(void) printf("%-14s %-14s %-14s %-14s %-14s\n",
			"k","R","E","P","slopeP");
		(void) printf("place[0] = %d, place[1] = %d\n",
			      place[0],place[1]);
		for(k = 0; k < jtotal; k++)
		{
			(void) printf("%-14d %-14g %-14g %-14g %-14g\n",
				k,R[k],E[k],P[k],slopeP[k]);
		}
		(void) printf("\n\n");
	}
	free_these(9,slp,x,y,z,dx,dy,c,d,e);
	DEBUG_LEAVE(init_PE_spline)
}		/*end init_PE_spline*/


EXPORT	void init_PE_phase_grid(
	SESAME_EOS	*seos,
	double		rho,
	PHASE_BDRY	*phase_bound,
	int		*place)
{
	double        *T = seos->de_params.T;
	double        *p = seos->de_params.P;
	double        *E = seos->de_params.E;
	double        *slopeP = seos->de_params.slopeP;
	double        *slopeE = seos->de_params.slopeE;
	double        *tbls = seos->sestab.tbls;
	double        rmin = Rho_min(seos);
	double        rmax = Rho_max(seos);
	double        tmin = Temp_min(seos);
	double        tmax = Temp_max(seos);
	double        Tb;
	double        rho_grid;
	double        ys, dz, er, pr;
	double        *xx, *yy, *zz, *slp;
	double        x1;
	double        *c, *d, *e, *dx, *dy;
	double        emin;
	double        y, z;
	double        rede[2], redp[2], r[2];
	double        V;
	int        flag;
	int        k, nr, nt, itotal, jtotal, jtot;
	int        i, j, imax, imin, jmax, jmin;
	int        jold, ii;
	int        n4, ntv;

	DEBUG_ENTER(init_PE_phase_grid)

	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);
	n4 = nr+nt+3*nr*nt+5;
	ntv = (int)(tbls[n4-1]);
	if (tbls[4] <= 0.) ntv = ntv-1;
	emin = tbls[nr+nt + nr*nt +4];
	emin = (emin <= 0.0) ? fabs(emin) + .01 : 0.0;
	for(i = 0; i < nr-1; i++)
	{
	    if (tbls[4+i]<= rmin && tbls[4+i+1]> rmin)  jmin = i;
	    if (tbls[4+i]< rmax && tbls[4+i+1]>= rmax) jmax = i+1;
	}
	for(i = 0; i < nt-1; i++)
	{
	    if (tbls[4+nr+i]<= tmin && tbls[4+nr+i+1] >tmin) imin = i;
	    if (tbls[4+nr+i]< tmax && tbls[4+nr+i+1]>= tmax) imax = i+1;
	}
	jtot = jmax - jmin +1;
	itotal = imax- imin + 9;
	jtotal = jmax - jmin + 9;

	uni_array(&c,nt+2,FLOAT);
	uni_array(&d,nt+2,FLOAT);
	uni_array(&e,nt+2,FLOAT);
	uni_array(&yy,nt+2,FLOAT);
	uni_array(&dy,nt+2,FLOAT);
	uni_array(&xx,nt+2,FLOAT);
	uni_array(&dx,nt+2,FLOAT);
	uni_array(&zz,nt+2,FLOAT);
	uni_array(&slp,nt+2,FLOAT);
	for (i = 0; i < nt+2; i++)
	{
	    c[i] = 0.0;
	    d[i] = 0.0;
	    e[i] = 0.0;
	    yy[i] = 0.0;
	    dy[i] = 0.0;
	    xx[i] = 0.0;
	    dx[i] = 0.0;
	    zz[i] = 0.0;
	    slp[i] = 0.0;
	}

	j = 0;
	jtot = 0;

	Tb = 0.0;

	rho_grid = rho;
	rho = ses_rt_rho_from_grid(rho,seos);
	V = 1.0/rho;
	if (get_phase_state(rho_grid,seos,phase_bound,&Tb,&er,&pr))
	{
	    x1 = Tb;
	    jtot++;
	}
	jold = 0;
	j = 0;

	for(i = imin+4; i <= imax+4; i++)
	{
	    ii = i-imin-4;
	    xx[j] = ses_rt_grid_from_temp(tbls[i+nr],seos);
	    if (jtot != 0 && j == 1 && 
	    fabs(xx[0]- xx[1])< TOLEPS && jold == 0)
	    {
	        j = 1;
	    }
	    else if ((jtot != 0) && (fabs(xx[0]-Tb) < TOLEPS) &&
	        (j == 0) && (jold == 0))
	    {
	        get_temp_state(xx[j],ntv,seos,phase_bound,r,rede,redp,&flag);
	        yy[j] = redp[0];
	        zz[j] = rede[1] + (r[0]*(rede[1] - rede[0])*
	                (r[1]*V - 1.0))/(r[0] - r[1]);
	        yy[j] = log(yy[j]);
	        zz[j] = log(zz[j] + emin);
	        jtot = 0;
	        j++;
	    }
	    else 
	    {
	        get_temp_state(xx[j],ntv,seos,phase_bound,r,rede,redp,&flag);
	        if (flag == 2 && Between(rho,r[0],r[1]))
	        {
	            yy[j] = redp[0];
	            zz[j] = rede[1] + (r[0]*(rede[1] - rede[0])*
	                (r[1]*V - 1.0))/(r[0] - r[1]);
	            yy[j] = log(yy[j]);
	            zz[j] = log(zz[j] + emin);
	        }
	        else 
	        {
	            lookspl(seos,phase_bound,rho_grid,
	                i,ii,itotal,jtotal,&y,&z,&ys,&dz);
	            zz[j] = z;
	            yy[j] = y;
	        }
	        j++;
	    }
	    if ((fabs(ses_rt_grid_from_temp(tbls[i+nr],seos) - x1) < TOLEPS) &&
	        (j > 1))
	    {
	        j--;
	        get_temp_state(xx[j],ntv,seos,phase_bound,r,rede,redp,&flag);
	        yy[j] = redp[0];
	        zz[j] = rede[1] + (r[0]*(rede[1] - rede[0])*
	                (r[1]*V - 1.0))/(r[0] - r[1]);
	        yy[j] = log(yy[j]);
	        zz[j] = log(zz[j] + emin);
	        j++;
	        if (splcomp(xx,yy,slp,j,c,d,e,dx,dy) >= 0)
	        {
	            screen("ERROR in init_PE_phase_grid(), ");
	            screen("splcomp failed\n");
	            clean_up (ERROR);
	        }
	        for(k = 0; k < j; k++)
	        {
	            p[k+jold] = yy[k];
	            slopeP[k+jold] = slp[k];
	        }
	        if (splcomp(xx,zz,slp,j,c,d,e,dx,dy) >= 0)
	        {
		screen("ERROR in init_PE_phase_grid(), ");
		screen("splcomp failed\n");
	            clean_up (ERROR);
	        }
	        for(k = 0; k < j; k++)
	        {
	            T[k+jold] = ses_rt_temp_from_grid(xx[k],seos);
	            slopeE[k+jold] = slp[k];
	            E[k+jold] = zz[k];
	        }
	        jold = j;
	        *place = jold-1;
	        xx[0] = x1;
	        yy[0] = yy[jold - 1];
	        zz[0] = zz[jold - 1];
	        jtot = 0;
	        j = 1;
	    }
	    else if ((ses_rt_grid_from_temp(tbls[i+nr],seos) < x1) &&
	        (ses_rt_grid_from_temp(tbls[i+nr+1],seos) > x1) &&
	        (i != imax+4) &&
	        (fabs(ses_rt_grid_from_temp(tbls[i+nr],seos) - x1) > TOLEPS) &&
	        (fabs(ses_rt_grid_from_temp(tbls[i+nr+1],seos) - x1) > TOLEPS))
	    {
	        xx[j] = x1;
	        get_temp_state(xx[j],ntv,seos,phase_bound,r,rede,redp,&flag);
	        yy[j] = redp[0];
	        if (r[0] < r[1])
	            zz[j] = rede[1] + (r[0]*(rede[1] - rede[0])*
	                (r[1]*V - 1.0))/(r[0] - r[1]);
	        else
	            zz[j] = rede[1];
	        yy[j] = log(yy[j]);
	        zz[j] = log(zz[j] + emin);
	        j++;
	        if (splcomp(xx,yy,slp,j,c,d,e,dx,dy) >= 0)
	        {
		screen("ERROR in init_PE_phase_grid(), ");
		screen("splcomp failed\n");
	            clean_up (ERROR);
	        }
	        for(k = 0; k < j; k++)
	        {
	            slopeP[k+jold] = slp[k];
	            p[k+jold] = yy[k];
	        }
	        if (splcomp(xx,zz,slp,j,c,d,e,dx,dy) >= 0)
	        if (j == 0)
	        {
		screen("ERROR in init_PE_phase_grid(), ");
		screen("splcomp failed\n");
	            clean_up (ERROR);
	        }
	        for(k = 0; k < j; k++)
	        {
	            T[k+jold] = ses_rt_temp_from_grid(xx[k],seos);
	            slopeE[k+jold] = slp[k];
	            E[k+jold] = zz[k];
	        }
	        jold = j;
	        *place = jold-1;
	        xx[0] = x1;
	        yy[0] = yy[jold-1];
	        zz[0] = zz[jold - 1];
	        jtot = 0;
	        j = 1;
	    }
	}
	if (j > 1)
	{
	    if (splcomp(xx,yy,slp,j,c,d,e,dx,dy) >= 0)
	    {
	        screen("ERROR in init_PE_phase_grid(), ");
	        screen("splcomp failed\n");
	        clean_up (ERROR);
	    }
	    for(k = 0; k < j; k++)
	    {
	        slopeP[k+jold] = slp[k];
	        p[k+jold] = yy[k];
	    }
	    if (splcomp(xx,zz,slp,j,c,d,e,dx,dy) >= 0)
	    {
	        screen("ERROR in init_PE_phase_grid(), ");
	        screen("splcomp failed\n");
	        clean_up (ERROR);
	    }
	    for(k = 0; k < j; k++)
	    {
	        slopeE[k+jold] = slp[k];
	        E[k+jold] = zz[k];
	        T[k+jold] = ses_rt_temp_from_grid(xx[k],seos);
	    }
	}
	for(j = k+jold; j < nt+2; j++)
	{
	    T[j] = -1;
	    E[j] = -1;
	    slopeE[j] = -1;
	    p[j] = -1;
	    slopeP[j] = -1;
	}
	if (DEBUG)
	{
	    for(j = 0; j < nt+2; j++)
	    {
	        (void) printf("In init_PE_phase_grid()\n");
	        (void) printf("%-14s %-14s %-14s %-14s %-14s %-14s\n",
	              "j","T","EXP(E)-EMIN","SLOPEE","P","SLOPEP");
	        (void) printf("%-14d %-14g %-14g %-14g %-14g %-14g\n",
	              j,T[j],exp(E[j])- emin,
	              slopeE[j],p[j],slopeP[j]);
	        (void) printf("\n\n");
	    }
	}
	if (*place == 0) *place = -1;
	free_these(9,dy,xx,dx,d,e,zz,yy,c,slp);
	DEBUG_LEAVE(init_PE_phase_grid)
}        /*end init_PE_phase_grid*/

EXPORT    void set_cross_states(
	Front        *fr,
	SESAME_EOS    *seos,
	PHASE_BDRY    *phase_bound,
	double        *Tvec,
	double        *Evec,
	double        *Pvec,
	double        *slopPvec,
	double        *slopEvec)
{
	int        j, flag;
	ORIENTATION    orient;
	CURVE        **cur, *C;
	NODE        **n, *node;
	BOND        *b;
	Locstate    state1, state2;
	INTERFACE    *intfc = fr->interf;
	double        *tbls = seos->sestab.tbls;
	int        nrhyp = Nrho_hyp(seos);
	int        nthyp = Ntemp_hyp(seos);
	int        nr, nt, ntv;
	int        i;
	double        tmin = Temp_min(seos);
	double        tmax = Temp_max(seos);
	double        rmin = Rho_min(seos);
	double        rmax = Rho_max(seos);
	double        rho, temp, *r, *E, *P, *et, *pt;
	double        adbgam, e, p, rp, csq;
	double        rmin_grid, rmax_grid, tmin_grid, tmax_grid;
	double        *rede, *redp, *colde, *coldp, *rph;
	double        emin, dT, de;

	DEBUG_ENTER(set_cross_states)

	nt = (int)(tbls[3]);
	nr = (int)(tbls[2]);
	emin = tbls[nr+nt + nr*nt +4];
	emin = (emin <= 0.0) ? fabs(emin) + 0.01 : 0.0;
	ntv = (int)(tbls[nr+nt+3*nr*nt+4]);
	if (tbls[4] <= 0.) ntv = ntv-1;

	uni_array(&r,nrhyp+2,FLOAT);
	uni_array(&E,nrhyp+2,FLOAT);
	uni_array(&P,nrhyp+2,FLOAT);
	for (i = 0; i < nrhyp+2; i++)
	{
	    r[i] = 0.0;
	    E[i] = 0.0;
	    P[i] = 0.0;
	}

	uni_array(&et,nrhyp+8,FLOAT);
	uni_array(&pt,nrhyp+8,FLOAT);
	for (i = 0; i < nrhyp+8; i++)
	{
	    et[i] = 0.0;
	    pt[i] = 0.0;
	}

	/*
	* LISA,     it is questionable to allocate uni_arrays of known length
	* dynamically either declare them explictly (double r[8],...)
	* or this is an error if 8 might be too small.
	* These arrays were reference up to value ntv.
	*/
	uni_array(&rph,ntv+8,FLOAT);
	uni_array(&rede,ntv+8,FLOAT);
	uni_array(&redp,ntv+8,FLOAT);
	uni_array(&colde,ntv+8,FLOAT);
	uni_array(&coldp,ntv+8,FLOAT);
	if (DEBUG)
	    (void) printf("In set_cross_states(), ntv = %d\n",ntv);
	for (i = 0; i < ntv+8; i++)
	{
	    rph[i] = 0.0;
	    rede[i] = 0.0;
	    redp[i] = 0.0;
	    colde[i] = 0.0;
	    coldp[i] = 0.0;
	}

	rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	tmin_grid = ses_rt_grid_from_temp(tmin,seos);


	/* Set adbgam on the phase boundary */

	for(cur = intfc->curves; cur && *cur; cur++)
	{
	    if (wave_type(*cur) != PHASE_BOUNDARY) continue;

	    temp = Coords((*cur)->start->posn)[1];
	    rho = Coords((*cur)->start->posn)[0];
	    get_temp_state(temp,ntv,seos,phase_bound,rph,rede,redp,&flag);
	    temp = ses_rt_temp_from_grid(temp,seos);
	    set_vectors(fr,r,E,P,seos,temp,Tvec,
	            Evec,Pvec,slopPvec,slopEvec);
	    adb_gam_state(fr,seos,r,P,E,temp,et,pt);
	    rho = ses_rt_rho_from_grid(rho,seos);
	    state1= right_start_state(*cur);
	    state2= left_start_state(*cur);
	    rp = ses_rt_redp(state1);
	    p = rp*rho*temp + ses_rt_coldp(state1);
	    e = (ses_rt_rede(state1))*temp + ses_rt_colde(state1);
	    adbgam = ses_rt_adb_gam(state1);

	    /* linear approx in mixed phase */
	    de = rph[1]*rph[0]*(rede[1] - rede[0])/(rph[0] - rph[1]);
	    de = -1.0*de/(rho*rho);
	    csq = (p/(rho*rho) - de)*adbgam;
	    ses_rt_adb_gam(state1) = rho*csq/p;
	    adbgam = ses_rt_adb_gam(state2);
	    get_de_and_dp(rho,fr,E,P,e,p,adbgam,
	        rede,redp,rph,et,pt,flag,&csq,seos);
	    if (csq <= 0.0)
	    {
	         csq = ses_rt_adb_gam(state1) + TOLEPS;
	        if (DEBUG)
	            (void) printf("Warning csq set positive on phbd\n");
	    }
	    ses_rt_adb_gam(state2) = csq;
	    node = (*cur)->start;
	    set_remaining_states(node,state1,state2,fr);
	    for (b = (*cur)->first; b != NULL; b = b->next)
	    {
	        rho = Coords(b->end)[0];
	        temp = Coords(b->end)[1];
	        get_temp_state(temp,ntv,seos,phase_bound,rph,
	                rede,redp,&flag);
	        temp = ses_rt_temp_from_grid(temp,seos);
	        set_vectors(fr,r,E,P,seos,temp,Tvec,
	                Evec,Pvec,slopPvec,slopEvec);
	        adb_gam_state(fr,seos,r,P,E,temp,et,pt);
	        rho = ses_rt_rho_from_grid(rho,seos);
	        state1 = right_state_at_point_on_curve(b->end,b,*cur);
	        state2= left_state_at_point_on_curve(b->end,b,*cur);
	        rp = ses_rt_redp(state1);
	        p = rp*rho*temp + ses_rt_coldp(state1);
	        e = (ses_rt_rede(state1))*temp + ses_rt_colde(state1);
	        adbgam = ses_rt_adb_gam(state1);

	        /* linear approx in mixed phase */

	        /* Guard against degneracy at top of phase bound */
	        if ( rede[1] < rede[0])
	        {
	            de = rph[1]*rph[0]*(rede[1] - rede[0])/
	                (rph[0] - rph[1]);
	            de = -1.0*de/(rho*rho);
	        }
	        else
	            de = 0.0;

	        csq = (p/(rho*rho) - de)*adbgam;
	        ses_rt_adb_gam(state1) = rho*csq/p;

	        adbgam = ses_rt_adb_gam(state2);
	        get_de_and_dp(rho,fr,E,P,e,p,adbgam,
	                rede,redp,rph,et,pt,flag,&csq,seos);
	        if (csq <= 0.0) 
	        {
	            csq = ses_rt_adb_gam(state1) + TOLEPS;
	            if (DEBUG)
	            (void) printf("Warning csq set positive on phbnd\n");
	        }
	        ses_rt_adb_gam(state2) = csq;
	        if (b == (*cur)->last)
	        {
	            node = (*cur)->end;
	            set_remaining_states(node,state1,state2,fr);
	        }
	    }
	}

	/* Set boundary adb gam at boundary nodes */

	for (n = intfc->nodes; n && *n; n++)
	{
	    if (is_node_of_phase_boundary(*n)) continue;

	    rho = Coords((*n)->posn)[0];
	    temp = Coords((*n)->posn)[1];
	    get_temp_state(temp,ntv,seos,phase_bound,rph,rede,redp,&flag);
	    temp = ses_rt_temp_from_grid(temp,seos);
	    set_vectors(fr,r,E,P,seos,temp,Tvec,Evec,Pvec,slopPvec,slopEvec);
	    adb_gam_state(fr,seos,r,P,E,temp,et,pt);
	    rho = ses_rt_rho_from_grid(rho,seos);

	    C = NULL;
	    for (cur = (*n)->in_curves; cur && *cur; cur++)
	    {
	        C = *cur;
	        orient = NEGATIVE_ORIENTATION;
	    }
	    if (C == NULL)
	    {
	    for (cur = (*n)->in_curves; cur && *cur; cur++)
	    {
	        C = *cur;
	        orient = POSITIVE_ORIENTATION;
	    }
	    }
	    if (C == NULL)
	    {
	    screen("ERROR in set_cross_states(), ");
	    screen("Curveless node found\n");
	    clean_up(ERROR);
	    }

	    state1 = Left_state_at_node(C,orient);
	    state2 = Right_state_at_node(C,orient);
	    rp = ses_rt_redp(state1);
	    p = rp*rho*temp + ses_rt_coldp(state1);
	    e = (ses_rt_rede(state1))*temp + ses_rt_colde(state1);
	    adbgam = ses_rt_adb_gam(state1);
	    get_de_and_dp(rho,fr,E,P,e,p,adbgam,rede,redp,rph,
	        et,pt,flag,&csq,seos);
	    if (csq <= 0.0) 
	    {
	    csq = TOLEPS;
	    if (DEBUG)
	        (void) printf("Warning csq set positive on curve\n");
	    }

	    for (cur = (*n)->in_curves; cur && *cur; cur++)
	    {
	    state1 = left_end_state(*cur);
	    state2 = right_end_state(*cur);
	    ses_rt_adb_gam(state1) = csq;
	    ses_rt_adb_gam(state2) = csq;
	    }
	    for (cur = (*n)->out_curves; cur && *cur; cur++)
	    {
	    state1 = left_start_state(*cur);
	    state2 = right_start_state(*cur);
	    ses_rt_adb_gam(state1) = csq;
	    ses_rt_adb_gam(state2) = csq;
	    }
	}

	/*Set states at boundary points */

	rmax_grid = ses_rt_grid_from_rho(rmax,seos);
	tmax_grid = ses_rt_grid_from_temp(tmax,seos);
	rmin_grid = ses_rt_grid_from_rho(rmin,seos);
	tmin_grid = ses_rt_grid_from_temp(tmin,seos);

	if (DEBUG) (void) printf("Setting boundary states\n");

	dT = fr->rect_grid->h[1];

	for(j = 1; j < nthyp; j++)
	{
	    temp = tmin_grid + j*dT;
	    rho = rmin_grid;
	    temp = ses_rt_temp_from_grid(temp,seos);
	    set_vectors(fr,r,E,P,seos,temp,
	            Tvec,Evec,Pvec,slopPvec,slopEvec);
	    adb_gam_state(fr,seos,r,P,E,temp,et,pt);
	    temp = ses_rt_grid_from_temp(temp,seos);
	    get_temp_state(temp,ntv,seos,phase_bound,rph,rede,redp,&flag);
	    rho = ses_rt_rho_from_grid(rho,seos);
	    for(cur = intfc->curves; cur && *cur; cur++)
	    {
	        if (is_exterior_comp(negative_component((*cur)),intfc))
	        {
	            if ((Coords((*cur)->end->posn)[0] <= rmin_grid) &&
	            (Coords((*cur)->start->posn)[1] <= temp) &&
	            (Coords((*cur)->start->posn)[0] <= rmin_grid))
	            {
	                C = *cur;
	            }
	        }
	    }
	    for (b = C->first; b != NULL; b = b->next)
	    {

	        if (!Between(temp,Coords(b->start)[1],
	                     Coords(b->end)[1]) ||
	            (temp >= Coords(b->start)[1])) continue;
	        if (b == (C)->last) continue;

	        state1= right_state_at_point_on_curve(b->end,b,C);
	        state2= left_state_at_point_on_curve(b->end,b,C);
	        rp = ses_rt_redp(state2);
	        p = rp*rho*ses_rt_temp_from_grid(temp,seos) +
	            ses_rt_coldp(state2);
	        e = (ses_rt_rede(state2))*ses_rt_temp_from_grid(temp,seos)
	            + ses_rt_colde(state2);
	        adbgam = ses_rt_adb_gam(state2);
	        get_de_and_dp(rho,fr,E,P,e,p,adbgam,rede,redp,
	            rph,et,pt,flag,&csq,seos);
	        if (csq <= 0.0) 
	        {
	            csq = TOLEPS;
	            if (DEBUG)
	            (void) printf("Warning csq set positive on curve\n");
	        }
	        ses_rt_adb_gam(state1) = csq;
	        ses_rt_adb_gam(state2) = csq;
	    }

	    rho = rmax_grid;
	    rho = ses_rt_rho_from_grid(rho,seos);
	    for(cur = intfc->curves; cur && *cur; cur++)
	    {
	        if (is_exterior_comp(negative_component((*cur)),intfc))
	        {
	            if ((Coords((*cur)->end->posn)[0] >= rmax_grid) &&
	            (Coords((*cur)->start->posn)[1] >= temp) &&
	            (Coords((*cur)->start->posn)[0] >= rmax_grid))
	            {
	                C = *cur;
	            }
	        }
	    }
	    for (b = (C)->first; b != NULL; b = b->next)
	    {
	        if (!Between(temp,Coords(b->start)[1],
	                     Coords(b->end)[1]) ||
	            (temp >= Coords(b->start)[1])) continue;
	        if (b == C->last) continue;

	        state1= right_state_at_point_on_curve(b->end,b,C);
	        state2= left_state_at_point_on_curve(b->end,b,C);
	        rp = ses_rt_redp(state2);
	        p = rp*rho*ses_rt_temp_from_grid(temp,seos) +
	                ses_rt_coldp(state2);
	        e = (ses_rt_rede(state2))*ses_rt_temp_from_grid(temp,seos)
	            + ses_rt_colde(state2);
	        adbgam = ses_rt_adb_gam(state2);
	        get_de_and_dp(rho,fr,E,P,e,p,adbgam,rede,redp,
	            rph,et,pt,flag,&csq,seos);
	        if (csq <= 0.0) 
	        {
	            csq = TOLEPS;
	            if (DEBUG)
	            (void) printf("Warning csq set positive on curve\n");
	        }
	        ses_rt_adb_gam(state1) = csq;
	        ses_rt_adb_gam(state2) = csq;
	    }
	}




	temp = tmin_grid;
	temp = ses_rt_temp_from_grid(temp,seos);
	set_vectors(fr,r,E,P,seos,temp,Tvec,Evec,Pvec,slopPvec,slopEvec);
	adb_gam_state(fr,seos,r,P,E,temp,et,pt);
	temp = ses_rt_grid_from_temp(temp,seos);
	get_temp_state(temp,ntv,seos,phase_bound,rph,rede,redp,&flag);
	for(cur = intfc->curves; cur && *cur; cur++)
	{
	    if (!is_exterior_comp(negative_component((*cur)),intfc))
	        continue;
	    if ((Coords((*cur)->end->posn)[1] <= temp) && 
	        (Coords((*cur)->start->posn)[1] <= temp))
	    {
	        C = *cur;
	        for (b = (C)->first; b != NULL; b = b->next)
	        {
	        if (b == (C)->last) continue;
	        rho = Coords(b->end)[0];
	        state1= right_state_at_point_on_curve(b->end,b,C);
	        state2= left_state_at_point_on_curve(b->end,b,C);
	        rho = ses_rt_rho_from_grid(rho,seos);
	        rp = ses_rt_redp(state2);
	        p = rp*rho*ses_rt_temp_from_grid(temp,seos) +
	            ses_rt_coldp(state2);
	        e = (ses_rt_rede(state2))*ses_rt_temp_from_grid(temp,seos) +
	                ses_rt_colde(state2);
	        adbgam = ses_rt_adb_gam(state2);
	        if (positive_component(C) == COMP_MIXED_PHASE)
	        {
	            /* linear approx in mixed phase */

	            /* Guard against degneracy at top of 
	                    phase bound */
	            if ( rede[1] < rede[0])
	            {
	                de = rph[1]*rph[0]*(rede[1] - rede[0])/
	                (rph[0] - rph[1]);
	                de = -1.0*de/(rho*rho);
	            }
	            else
	                de = 0.0;
	            csq = (p/(rho*rho) - de)*adbgam;
	            if (csq <= 0.0) 
	            {
	                csq = TOLEPS;
	                if (DEBUG)
	                (void) printf("csq < 0 in mixed phase\n");
	            }
	            ses_rt_adb_gam(state1) = rho*csq/p;
	            ses_rt_adb_gam(state2) = rho*csq/p;

	        }
	        else
	        {
	            get_de_and_dp(rho,fr,E,P,e,p,adbgam,rede,redp,
	                rph,et,pt,flag,&csq,seos);
	            if (csq <= 0.0)
	            {
	                csq = TOLEPS;
	                if (DEBUG)
	                (void) printf("csq < 0 om curve\n");
	            }
	            ses_rt_adb_gam(state1) = csq;
	            ses_rt_adb_gam(state2) = csq;
	        }
	        }
	    }
	}

	temp = tmax_grid;
	temp = ses_rt_temp_from_grid(temp,seos);
	set_vectors(fr,r,E,P,seos,temp,Tvec,Evec,Pvec,slopPvec,slopEvec);
	adb_gam_state(fr,seos,r,P,E,temp,et,pt);
	temp = ses_rt_grid_from_temp(temp,seos);
	get_temp_state(temp,ntv,seos,phase_bound,rph,rede,redp,&flag);
	for(cur = intfc->curves; cur && *cur; cur++)
	{
	    if (!is_exterior_comp(negative_component((*cur)),intfc))
	        continue;
	    if ((Coords((*cur)->end->posn)[1] >= temp) && 
	        (Coords((*cur)->start->posn)[1] >= temp))
	    {
	        C = *cur;
	        for (b = (C)->first; b != NULL; b = b->next)
	        {
	        /* End state of curve already set */
	        if (b == (C)->last) continue;
	        rho = Coords(b->end)[0];
	        state1= right_state_at_point_on_curve(b->end,b,C);
	        state2= left_state_at_point_on_curve(b->end,b,C);
	        rho = ses_rt_rho_from_grid(rho,seos);
	        rp = ses_rt_redp(state2);
	        p = rp*rho*ses_rt_temp_from_grid(temp,seos) +
	            ses_rt_coldp(state2);
	        e = (ses_rt_rede(state2))*ses_rt_temp_from_grid(temp,seos) +
	                ses_rt_colde(state2);
	        adbgam = ses_rt_adb_gam(state2);
	        if (positive_component(C) == COMP_MIXED_PHASE)
	        {
	            /* linear approx in mixed phase */

	            /* Guard against degneracy at top of 
	                    phase bound */
	            if ( rede[1] < rede[0])
	            {
	                de = rph[1]*rph[0]*(rede[1] - rede[0])/
	                (rph[0] - rph[1]);
	                de = -1.0*de/(rho*rho);
	            }
	            else
	                de = 0.0;
	            csq = (p/(rho*rho) - de)*adbgam;
	            ses_rt_adb_gam(state1) = rho*csq/p;
	            ses_rt_adb_gam(state2) = rho*csq/p;

	        }
	        else
	        {
	            get_de_and_dp(rho,fr,E,P,e,p,adbgam,rede,redp,
	                rph,et,pt,flag,&csq,seos);
	            if (csq <= 0.0)
	            {
	                csq = TOLEPS;
	                if (DEBUG)
	                (void) printf("csq < 0 om curve\n");
	            }
	            ses_rt_adb_gam(state1) = csq;
	            ses_rt_adb_gam(state2) = csq;
	        }
	        }
	    }
	}

	interpolate_intfc_states(intfc) = NO;

	free_these(5,r,P,E,et,pt);
	free_these(5,rph,rede,redp,colde,coldp);

	DEBUG_LEAVE(set_cross_states)
}        /*end set_cross_states*/

LOCAL    int is_node_of_phase_boundary(
	NODE        *n)
{
	CURVE        **cur;

	for (cur = n->in_curves; cur && *cur; cur++)
	{
	    if (wave_type(*cur) == PHASE_BOUNDARY)
	        return YES;
	}
	for (cur = n->out_curves; cur && *cur; cur++)
	{
	    if (wave_type(*cur) == PHASE_BOUNDARY)
	        return YES;
	}
	return NO;
}        /*end is_node_of_phase_boundary*/

LOCAL void get_de_and_dp(
	double		rho,
	Front		*fr,
	double		*Evec,
	double		*Pvec,
	double		e,
	double		p,
	double		adbgam,
	double		*eph,
	double		*pph,
	double		*rph,
	double		*et,
	double		*pt,
	int		flag,
	double		*csq,
	SESAME_EOS	*seos)
{
	double        rmin_grid = fr->rect_grid->L[0];
	double        rmax_grid = fr->rect_grid->U[0];
	double        drho;
	double        lrho = rho;
	double        e1, e2, p1, p2, r1, r2;
	double        slpe1, slpe2, slpp1, slpp2, val, de, dp;
	double        lrph[2], leph[2], lpph[2];
	int        nrhyp = fr->rect_grid->gmax[0];
	int        nintst, numph;
	int        i, i1,    j;

	DEBUG_ENTER(get_de_and_dp)

	lrho = ses_rt_grid_from_rho(lrho,seos);
	drho = (rmax_grid - rmin_grid)/(nrhyp);
	nintst = 0;

	for (i= 0; i < flag; i++)
	{
	    if (Between(ses_rt_grid_from_rho(rph[i],seos),rmin_grid,rmax_grid))
	    {
	        lrph[nintst] = ses_rt_grid_from_rho(rph[i],seos);
	        leph[nintst] = eph[i];
	        lpph[nintst] = pph[i];
	        nintst++;
	    }
	}

	numph = nintst;
	for (i = 0; i < nrhyp; i++)
	{
	    for (j = 0; j < nintst; j++)
	    {
	        if (fabs(lrph[j] - rmin_grid+i*drho) < TOLEPS*drho)
	        {
	        numph = nintst - j - 1;
	        }
	    }
	}

	if (lrho <=  rmin_grid)
	{
	    de = et[0];
	    dp = pt[0];
	    lrho = ses_rt_rho_from_grid(lrho,seos);
	    *csq = (p/lrho)*dp + adbgam*(p/(lrho*lrho) - 
	        (e/lrho)*de);
	    *csq = *csq*lrho/p;
	    if (*csq <=0.0) *csq = TOLEPS;
	    DEBUG_LEAVE(get_de_and_dp)
	    return;
	}

	if (lrho >=  rmax_grid)
	{
	    de = et[nrhyp+nintst+numph];
	    dp = pt[nrhyp+nintst+numph];
	    lrho = ses_rt_rho_from_grid(lrho,seos);

	    *csq = (p/lrho)*dp + adbgam*(p/(lrho*lrho) - 
	        (e/lrho)*de);
	    *csq = *csq*lrho/p;
	    if (*csq <=0.0) *csq = TOLEPS;
	    DEBUG_LEAVE(get_de_and_dp)
	    return;
	}
	/* Get the derivative */

	for (i = 0; i < nrhyp; i++)
	{
	    if (Between(lrho,rmin_grid+i*drho,rmin_grid+(i+1)*drho))
	    {
	        r1 = rmin_grid+i*drho;
	        e1 = Evec[i];
	        p1 = Pvec[i];
	        slpe1 = et[i];
	        slpp1 = pt[i];
	        r2 = rmin_grid + (i+1)*drho;
	        e2 = Evec[i+1];
	        p2 = Pvec[i+1];
	        slpe2 = et[i+1];
	        slpp2 = pt[i+1];
	        i1 = i;
	        break;
	    }
	}

	for (j = 0; j < nintst; j++)
	{
	    if (fabs(lrph[j] - lrho) < TOLEPS*drho && j == 0)
	    {
	        if (nintst == 2)
	        {
	            de = (lrph[0] <= r1) ? et[i1] : et[i1+1];
	            dp = (lrph[0] <= r1) ? pt[i1] : pt[i1+1];
	        }
	        else
	        {
	            de = et[i1+numph+1];
	            dp = pt[i1+numph+1];
	        }
	        lrho = ses_rt_rho_from_grid(lrho,seos);
	        *csq = (p/lrho)*dp + adbgam*(p/(lrho*lrho) - 
	            (e/lrho)*de);
	        *csq = *csq*lrho/p;
		DEBUG_LEAVE(get_de_and_dp)
	        return;
	    }
	    if (fabs(lrph[j] - lrho) < TOLEPS*drho && j == 1)
	    {
	        de = et[i1+nintst+numph];
	        dp = pt[i1+nintst+numph];
	        lrho = ses_rt_rho_from_grid(lrho,seos);
	        *csq = (p/lrho)*dp + adbgam*(p/(lrho*lrho) - 
	            (e/lrho)*de);
	        *csq = *csq*lrho/p;
		DEBUG_LEAVE(get_de_and_dp)
	        return;
	    }
	}

	for (j = 0; j < nintst; j++)
	{
	    if (r1 < lrph[j] && lrph[j] < lrho) 
	    {
	        r1 = lrph[j];
	        e1 = leph[j];
	        p1 = lpph[j];
	    }
	    if (r2 > lrph[j] && lrph[j] > lrho)
	    {
	        r2 = lrph[j];
	        e2 = leph[j];
	        p2 = lpph[j];
	    }
	}

	if ((nintst == 2 && r1 >= lrph[1]) || (nintst == 1 && r1 >= lrph[0]))
	{
	    slpe1 = et[i1+nintst+numph];
	    slpp1 = pt[i1+nintst+numph];
	    slpe2 = et[i1+nintst+numph+1];
	    slpp2 = pt[i1+nintst+numph+1];
	}

	spline(r1,r2,p1,p2,slpp1,slpp2,lrho,&val,&dp);
	spline(r1,r2,e1,e2,slpe1,slpe2,lrho,&val,&de);
	lrho = ses_rt_rho_from_grid(lrho,seos);

	*csq = (p/lrho)*dp + adbgam*(p/(lrho*lrho) - (e/lrho)*de);
	*csq = *csq*lrho/p;
	if (*csq <=0.0) *csq = TOLEPS;

	DEBUG_LEAVE(get_de_and_dp)
	return;
}        /*end get_de_and_dp*/

/*
*            adb_gam_state():
*
*    Set dE/drph and dP/drph for adiabatic gamma calcultion. If rph is
*    on the phase boundary then the appropriate values are determined.
*/

LOCAL    void adb_gam_state(
	Front        *fr,
	SESAME_EOS    *seos,
	double        *r,
	double        *P,
	double        *E,
	double        temp,
	double        *et,
	double        *pt)
{
	double        *rho, *rede, *redp, *coldp, *colde;
	double        y;
	double        *dx, *dy, *xx, *yy, *zz, *c, *d, *e, *slp;
	double        *tbls = seos->sestab.tbls;
	double        emin;
	double        a, b;
	int        nrhyp = Nrho_hyp(seos);
	int        i, j, ii, k, jold, flag;
	int        nr, nt;

	DEBUG_ENTER(adb_gam_state)

	/*
	* LISA,     it is questionable to allocate uni_arrays of known length
	* dynamically. Either declare them explictly (double r[8],...)
	* or this is an error if 8 might be too small.
	*/
	uni_array(&rho,8,FLOAT);
	uni_array(&rede,8,FLOAT);
	uni_array(&redp,8,FLOAT);
	uni_array(&colde,8,FLOAT);
	uni_array(&coldp,8,FLOAT);
	for (i = 0; i < 8; i++)
	{
	    rho[i] = 0.0;
	    rede[i] = 0.0;
	    redp[i] = 0.0;
	    colde[i] = 0.0;
	    coldp[i] = 0.0;
	}

	uni_array(&xx,nrhyp+8,FLOAT);
	uni_array(&yy,nrhyp+8,FLOAT);
	uni_array(&zz,nrhyp+8,FLOAT);
	uni_array(&c,nrhyp+8,FLOAT);
	uni_array(&d,nrhyp+8,FLOAT);
	uni_array(&dx,nrhyp+8,FLOAT);
	uni_array(&dy,nrhyp+8,FLOAT);
	uni_array(&e,nrhyp+8,FLOAT);
	uni_array(&slp,nrhyp+8,FLOAT);
	for (i = 0; i < nrhyp+8; i++)
	{
	    xx[i] = 0.0;
	    yy[i] = 0.0;
	    zz[i] = 0.0;
	    c[i] = 0.0;
	    d[i] = 0.0;
	    dx[i] = 0.0;
	    dy[i] = 0.0;
	    e[i] = 0.0;
	    slp[i] = 0.0;
	}

	y = ses_rt_grid_from_temp(temp,seos);
	ii = 0;
	jold = 0;
	nr = (int)(tbls[2]);
	nt = (int)(tbls[3]);
	emin = tbls[nr+nt + nr*nt +4];
	emin = (emin <= 0.0) ? fabs(emin) + 0.01 : 0.0;
	get_phase_temp_state(y,fr,rho,rede,redp,&flag,seos);
	for(i = 0; i < nrhyp+1; i++)
	{
	    if (flag >0)
	    {
	        for(j = 0; j < flag; j++)
	        {
	        if ((ii == 1) &&
	        (fabs(rho[j] - ses_rt_grid_from_rho(r[i],seos)) <= TOLEPS))
	            i++;
	        }
	    }
	    xx[ii] = ses_rt_grid_from_rho(r[i],seos);
	    if ((ii == 1) && (fabs(xx[0] - xx[1]) < TOLEPS)) ii --;
	    yy[ii] = E[i];
	    zz[ii] = P[i];
	    ii++;
	    if (flag > 0)
	    {
	        for(j = 0; j < flag; j++)
	        {
	            if ((i != nrhyp-1) && (i != 0) &&
	                (ses_rt_rho_from_grid(rho[j],seos)<= r[i+1])
	                && 
	                (ses_rt_rho_from_grid(rho[j],seos) > r[i]))
	            {
	                xx[ii] = rho[j];
	                yy[ii] = rede[j] +emin;
	                yy[ii] = log(yy[ii]);
	                zz[ii] = redp[j];
	                zz[ii] = log(zz[ii]);
	                if (fabs(xx[0] - xx[1]) < TOLEPS)
	                    break;
	                a = (zz[1] - zz[0])/(xx[1] - xx[0]);
	                b = (zz[ii] - zz[ii -1])/
	                    (xx[ii] - xx[ii-1]);
	                ii++;
	                if (splcomp2(xx,zz,slp,ii,c,d,e,dx,dy,a,b) >= 0)
	                {
	                    screen("ERROR in adb_gam_state(), ");
	                    screen("splcomp2 failed\n");
	                    clean_up (ERROR);
	                }
	                for(k = 0; k < ii; k++)
	                {
	                    pt[k+jold] = slp[k];
	                }
	                ii--;
	                a = (yy[1] - yy[0])/(xx[1] - xx[0]);
	                b = (yy[ii] - yy[ii -1])/(xx[ii] - xx[ii-1]);
	                ii++;
	                if (splcomp2(xx,yy,slp,ii,c,d,e,dx,dy,a,b) >= 0)
	                {
	                    screen("ERROR in adb_gam_state(), ");
	                    screen("splcomp2 failed\n");
	                    clean_up (ERROR);
	                }
	                for(k = 0; k < ii; k++)
	                {
	                    et[k+jold] = slp[k];
	                }
	                jold = jold + ii;
	                ii = 0;
	                xx[ii] = rho[j];
	                yy[ii] = rede[j] + emin;
	                zz[ii]    = redp[j];
	                yy[ii] = log(yy[ii]);
	                zz[ii] = log(zz[ii]);
	                ii = 1;
	            }
	        }
	    }
	}
	if (ii > 1)
	{
	    ii--;
	    a = (zz[1] - zz[0])/(xx[1] - xx[0]);
	    b = (zz[ii] - zz[ii -1])/(xx[ii] - xx[ii-1]);
	    ii++;
	    if (splcomp2(xx,zz,slp,ii,c,d,e,dx,dy,a,b) >= 0)
	    {
	        screen("ERROR in adb_gam_state(), ");
	        screen("splcomp2 failed\n");
	        clean_up (ERROR);
	    }
	    for(k = 0; k < ii; k++)
	    {
	        pt[k+jold] = slp[k];
	    }
	    ii--;
	    a = (yy[1] - yy[0])/(xx[1] - xx[0]);
	    b = (yy[ii] - yy[ii -1])/(xx[ii] - xx[ii-1]);
	    ii++;
	    if (splcomp2(xx,yy,slp,ii,c,d,e,dx,dy,a,b) >= 0)
	    {
	    screen("ERROR in adb_gam_state(), ");
	    screen("splcomp2 failed\n");
	        clean_up (ERROR);
	    }
	    for(k = 0; k < ii; k++)
	    {
	        et[k+jold] = slp[k];
	    }
	}
	if (DEBUG)
	{
	    (void) printf("In adb_gam_state, pt and et are\n");
	    for (i = 0; i < nrhyp+4; i++)
	    {
	        (void) printf("pt[%d] = %g, et[%d] = %g\n",
	                  i,pt[i],i,et[i]);
	    }
	    (void) printf("\n");
	}
	free_these(5,rho,rede,redp,colde,coldp);
	free_these(9,xx,yy,zz,slp,c,d,e,dx,dy);
	DEBUG_LEAVE(adb_gam_state)
}        /*end adb_gam_state*/

LOCAL    void set_vectors(
	Front        *fr,
	double        *r,
	double        *E,
	double        *P,
	SESAME_EOS    *seos,
	double        temp,
	double        *Tvec,
	double        *Evec,
	double        *Pvec,
	double        *slopPvec,
	double        *slopEvec)
{
	int        nrhyp = Nrho_hyp(seos);
	double        tmin = Temp_min(seos);
	double        tmax = Temp_max(seos);
	double        *tbls = seos->sestab.tbls;
	double        rmin_grid = fr->rect_grid->L[0];
	double        t1,t2,e1,e2,p1,p2,s1,s2,dp;
	double        val;
	double        ltemp = temp;
	double        drho;
	int        i,j;
	int        nt;

	DEBUG_ENTER(set_vectors)

	nt = (int)(tbls[3]);
	ltemp = ses_rt_grid_from_temp(ltemp,seos);
	drho = fr->rect_grid->h[0];
	for(j = 0; j < nrhyp+1; j++)
	{
	    r[j] = ses_rt_rho_from_grid(rmin_grid + j*drho,seos);
	    for(i = 0; i < nt+1; i++)
	    {
	        if ((temp <= Tvec[j*(nt+2)+i+1]) &&
	            (temp > Tvec[j*(nt+2)+i]) &&
	            (Tvec[j*(nt+2)+i+1] > 0.0))
	        {
	            t1 = ses_rt_grid_from_temp(Tvec[j*(nt+2)+i],seos);
	            t2 = ses_rt_grid_from_temp(Tvec[j*(nt+2)+i+1],seos);
	            e1 = Evec[j*(nt+2)+i];
	            e2 = Evec[j*(nt+2)+i+1];
	            p1 = Pvec[j*(nt+2)+i];
	            p2 = Pvec[j*(nt+2)+i+1];
	            s1 = slopEvec[j*(nt+2)+i];
	            s2 = slopEvec[j*(nt+2)+i+1];
	            spline(t1,t2,e1,e2,s1,s2,ltemp,&val,&dp);
	            E[j] = val;
	            s1 = slopPvec[j*(nt+2)+i];
	            s2 = slopPvec[j*(nt+2)+i+1];
	            spline(t1,t2,p1,p2,s1,s2,ltemp,&val,&dp);
	            P[j] = val;
	            break;
	        }
	        else if (temp <= Tvec[j*(nt+2)])
	        {
	            E[j] = Evec[j*(nt+2)];
	            P[j] = Pvec[j*(nt+2)];
	            break;
	        }
	        else if (max(fabs(Tvec[j*(nt+2)+i] - temp),
	        fabs(Tvec[j*(nt+2)+i+1]- temp)) <= TOLEPS*(tmax-tmin))
	        {
	            E[j] = Evec[j*(nt+2)+i];
	            P[j] = Pvec[j*(nt+2)+i];
	            break;
	        }
	    }
	}
	DEBUG_LEAVE(set_vectors)
}        /*end set_vectors*/

/*
*            set_remaining_states()
*
*    Set states on nodes of curves adjacent to phase boundary 
*    The phase boundary  is oriented with COMP_MIXED_PHASE on the right.
*/

LOCAL    void set_remaining_states(
	NODE        *n,
	Locstate    str,
	Locstate    stl,
	Front        *fr)
{
	INTERFACE *intfc = fr->interf;
	CURVE        **cur;
	Locstate state, state1;

	DEBUG_ENTER(set_remaining_states)

	for (cur = n->in_curves; cur && *cur; cur++)
	{
	    if (!is_exterior_comp(negative_component((*cur)),intfc))
	        continue;
	    state = right_end_state(*cur);
	    state1 = left_end_state(*cur);
	    if (positive_component((*cur)) == COMP_MIXED_PHASE)
	    {
	        ses_rt_rede(state) = ses_rt_rede(str);
	        ses_rt_redp(state) = ses_rt_redp(str);
	        ses_rt_colde(state) = ses_rt_colde(str);
	        ses_rt_coldp(state) = ses_rt_coldp(str);
	        ses_rt_adb_gam(state) = ses_rt_adb_gam(str);
	        ses_rt_gru_gam(state) = ses_rt_gru_gam(str);
	        ses_rt_S(state) = ses_rt_S(str);
	        ses_rt_rede(state1) = ses_rt_rede(str);
	        ses_rt_redp(state1) = ses_rt_redp(str);
	        ses_rt_colde(state1) = ses_rt_colde(str);
	        ses_rt_coldp(state1) = ses_rt_coldp(str);
	        ses_rt_S(state1) = ses_rt_S(str);
	        ses_rt_adb_gam(state1) = ses_rt_adb_gam(str);
	        ses_rt_gru_gam(state1) = ses_rt_gru_gam(str);
	    }
	    else
	    {
	        ses_rt_rede(state) = ses_rt_rede(stl);
	        ses_rt_redp(state) = ses_rt_redp(stl);
	        ses_rt_colde(state) = ses_rt_colde(stl);
	        ses_rt_coldp(state) = ses_rt_coldp(stl);
	        ses_rt_adb_gam(state) = ses_rt_adb_gam(stl);
	        ses_rt_gru_gam(state) = ses_rt_gru_gam(stl);
	        ses_rt_S(state) = ses_rt_S(stl);
	        ses_rt_rede(state1) = ses_rt_rede(stl);
	        ses_rt_redp(state1) = ses_rt_redp(stl);
	        ses_rt_colde(state1) = ses_rt_colde(stl);
	        ses_rt_coldp(state1) = ses_rt_coldp(stl);
	        ses_rt_S(state1) = ses_rt_S(stl);
	        ses_rt_adb_gam(state1) = ses_rt_adb_gam(stl);
	        ses_rt_gru_gam(state1) = ses_rt_gru_gam(stl);
	    }
	}
	for (cur = n->out_curves; cur && *cur; cur++)
	{
	    if (!is_exterior_comp(negative_component((*cur)),intfc))
	        continue;
	    state = right_start_state(*cur);
	    state1 = left_start_state(*cur);
	    if (positive_component((*cur)) == COMP_MIXED_PHASE)
	    {
	        ses_rt_rede(state) = ses_rt_rede(str);
	        ses_rt_redp(state) = ses_rt_redp(str);
	        ses_rt_colde(state) = ses_rt_colde(str);
	        ses_rt_coldp(state) = ses_rt_coldp(str);
	        ses_rt_S(state) = ses_rt_S(str);
	        ses_rt_adb_gam(state) = ses_rt_adb_gam(str);
	        ses_rt_gru_gam(state) = ses_rt_gru_gam(str);
	        ses_rt_rede(state1) = ses_rt_rede(str);
	        ses_rt_redp(state1) = ses_rt_redp(str);
	        ses_rt_colde(state1) = ses_rt_colde(str);
	        ses_rt_coldp(state1) = ses_rt_coldp(str);
	        ses_rt_S(state1) = ses_rt_S(str);
	        ses_rt_adb_gam(state1) = ses_rt_adb_gam(str);
	        ses_rt_gru_gam(state1) = ses_rt_gru_gam(str);
	    }
	    else
	    {
	        ses_rt_rede(state) = ses_rt_rede(stl);
	        ses_rt_redp(state) = ses_rt_redp(stl);
	        ses_rt_colde(state) = ses_rt_colde(stl);
	        ses_rt_coldp(state) = ses_rt_coldp(stl);
	        ses_rt_S(state) = ses_rt_S(stl);
	        ses_rt_adb_gam(state) = ses_rt_adb_gam(stl);
	        ses_rt_gru_gam(state) = ses_rt_gru_gam(stl);
	        ses_rt_rede(state1) = ses_rt_rede(stl);
	        ses_rt_redp(state1) = ses_rt_redp(stl);
	        ses_rt_colde(state1) = ses_rt_colde(stl);
	        ses_rt_coldp(state1) = ses_rt_coldp(stl);
	        ses_rt_S(state1) = ses_rt_S(stl);
	        ses_rt_adb_gam(state1) = ses_rt_adb_gam(stl);
	        ses_rt_gru_gam(state1) = ses_rt_gru_gam(stl);
	    }
	}
	interpolate_intfc_states(intfc) = NO;
	DEBUG_LEAVE(set_remaining_states)
}        /*end set_remaining_states*/

EXPORT    void get_temp_state(
	double		temp,
	int		n,
	SESAME_EOS	*seos,
	PHASE_BDRY	*phase_bound,
	double		*r,
	double		*rede,
	double		*redp,
	int		*flag)
{
	double       *rho = phase_bound->rvar;
	double       *T = phase_bound->Tvar;
	double       *Pph = phase_bound->rp;
	double       *Eph = phase_bound->re;
	int         i;
	static const double MIN_TEMP_DIFF = 1.0e-2; /*TOLERANCE*/

	DEBUG_ENTER(get_temp_state)
	*flag = 0;
	if (ses_rt_temp_from_grid(temp,seos) > T[n])
	{
	    if (DEBUG)
	    {
	        (void) printf("Temperature above max, no vars set\n");
	        (void) printf("temp = %g, ",temp);
		(void) printf("ses_rt_temp_from_grid(temp,seos) = %g, ",
	            	      ses_rt_temp_from_grid(temp,seos));
	        (void) printf("T[%d] = %g\n",n,T[n]);
	    }
	    DEBUG_LEAVE(get_temp_state)
	    return;
	}

	for(i = 0; i < n; i++)
	{
	    if (fabs(ses_rt_temp_from_grid(temp,seos) - T[i])< MIN_TEMP_DIFF)
	    {
	        r[0] = rho[i];
	        redp[0] = Pph[i];
	        rede[0] = Eph[i];
	        r[1] = rho[2*n-1-i];
	        redp[1] = Pph[2*n-1-i];
	        rede[1] = Eph[2*n-1-i];
	        *flag = 2;
	        if (DEBUG)
	        {
	            (void) printf("Temperature near grid value\n");
	        }
		DEBUG_LEAVE(get_temp_state)
	        return;
	    }
	}

	for(i = 0; i < n-1; i++)
	{
	    if ((T[i] < ses_rt_temp_from_grid(temp,seos)) &&
	        (T[i+1] > ses_rt_temp_from_grid(temp,seos)))
	    {
	        r[0] = ((ses_rt_grid_from_rho(rho[i+1],seos) -
	                ses_rt_grid_from_rho(rho[i],seos))/
	                (ses_rt_grid_from_temp(T[i+1],seos) -
	                ses_rt_grid_from_temp(T[i],seos)))*
	                (temp - ses_rt_grid_from_temp(T[i],seos)) +
	                ses_rt_grid_from_rho(rho[i],seos);
	        r[0] = ses_rt_rho_from_grid(r[0],seos);
	        redp[0] = ((log(Pph[i+1]) - log(Pph[i]))/
	            (ses_rt_grid_from_temp(T[i+1],seos) -
	            ses_rt_grid_from_temp(T[i],seos)))*
	            (temp - ses_rt_grid_from_temp(T[i],seos)) +
	            log(Pph[i]);
	        redp[0] = exp(redp[0]);
	        rede[0] = ((Eph[i+1] - Eph[i])/
	                (ses_rt_grid_from_temp(T[i+1],seos) - 
	                ses_rt_grid_from_temp(T[i],seos)))*
	                (temp - ses_rt_grid_from_temp(T[i],seos)) +
	                Eph[i];
	        *flag = 2;
	        if (DEBUG)
	        {
	            (void) printf("Temperature between T[%d] and T[%d]\n",
	                i,i+1);

	        }
	    }
	    if ((T[n+i] > ses_rt_temp_from_grid(temp,seos)) &&
	        (T[n+i+1] < ses_rt_temp_from_grid(temp,seos)))
	    {
	        r[1] = ((ses_rt_grid_from_rho(rho[n+i+1],seos) -
	                ses_rt_grid_from_rho(rho[n+i],seos))/
	                (ses_rt_grid_from_temp(T[n+i+1],seos) -
	                ses_rt_grid_from_temp(T[n+i],seos)))*
	                (temp - ses_rt_grid_from_temp(T[n+i],seos)) +
	                ses_rt_grid_from_rho(rho[n+i],seos);
	        r[1] = ses_rt_rho_from_grid(r[1],seos);
	        redp[1] = ((log(Pph[n+i+1]) - log(Pph[n+i]))/
	                (ses_rt_grid_from_temp(T[n+i+1],seos) -
	                ses_rt_grid_from_temp(T[n+i],seos)))*
	                (temp - ses_rt_grid_from_temp(T[n+i],seos)) +
	                log(Pph[n+i]);
	        redp[1] = exp(redp[1]);
	        rede[1] = ((Eph[n+i+1] - Eph[n+i])/
	                (ses_rt_grid_from_temp(T[n+i+1],seos) - 
	                ses_rt_grid_from_temp(T[n+i],seos)))*
	                (temp - ses_rt_grid_from_temp(T[n+i],seos)) +
	                Eph[n+i];
	        if (DEBUG)
	        {
	            (void) printf("Temperature between T[%d] and T[%d]\n",
	                n+i,n+i+1);

	        }
	    }
	}
	DEBUG_LEAVE(get_temp_state)

}        /*end get_temp_state*/
#endif /* defined(TWOD) && defined(PHASE_CODE) */
