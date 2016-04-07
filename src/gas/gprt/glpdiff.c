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
*			glpdiff.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if defined(TWOD) || defined(THREED)

#include <gdecs/gdecs.h>

	/* LOCAL Function Declarations */
LOCAL	boolean	is_block_on_front(int*,Wave*,Front*,double*,int);
LOCAL	void	L1_diff(Grid*,Wave*,Front*,OUTPUT_DATA*);
LOCAL	void	L2_diff3d(Grid*,Wave*,Front*,OUTPUT_DATA*);
LOCAL	void	L2_diff2d(Grid*,Wave*,Front*,OUTPUT_DATA*);
LOCAL	void	L2_regular_integral(Locstate,Locstate,Locstate,double*,double*,
				    double*,double*,double*,double*,int**,Wave*,
				    Front*,Locstate*);
LOCAL	void	subdomain_integral(Locstate,Locstate,Locstate,double*,double*,
				   double*,double*,double*,double*,double*,double*,
				   Wave*,Front*,OUTPUT_DATA*);
LOCAL	void	zero_states(Locstate,Locstate,Locstate,double*,double*,double*,
			    double*,double*,double*,Front*);
LOCAL	void	add_product_to_total(Locstate,Locstate,Locstate,double*,double*,
				     double*,double*,double*,double*,int,Locstate,
				     Locstate,Locstate,Locstate,Locstate,
				     Locstate,double,double,double,double,double,
				     double,double*,double*,double*,double*,
				     double*,double*);
LOCAL	void	add_multiple_of_product_to_total(Locstate,Locstate,Locstate,
						 double*,double*,double*,double*,
						 double*,double*,Locstate,
						 Locstate,Locstate,Locstate,
						 Locstate,Locstate,double,
						 double,double,double,double,
						 double,double*,double*,double*,
						 double*,double*,double*,
						 Front*,double);
LOCAL	void	state_sqrt(Locstate,Locstate,Locstate,double*,double*,double*,
			   double*,double*,double*,int);
LOCAL	void	absolute_value(Locstate,Locstate,Locstate,double*,double*,
			       double*,double*,double*,double*,int);
LOCAL	void	multiplied_by_factor(Locstate,Locstate,Locstate,double*,double*,
				     double*,double*,double*,double*,int,double);
LOCAL	void	add_to_total(Locstate,Locstate,Locstate,double*,double*,double*,
			     double*,double*,double*,Locstate,Locstate,Locstate,
			     double,double,double,double*,double*,double*,int);
LOCAL	void	find_diff_and_other_states(Locstate,double*,double*,double*,
					   double*,double*,double*,Locstate,
					   Locstate,int);
LOCAL	void	Lp_printout(OUTPUT_DATA*,Grid*,Locstate,Locstate,Locstate,
			    double,double,double,double*,double*,double*,int);

#if !defined(__INTEL_COMPILER)
#pragma noinline	L1_diff
#pragma	noinline	L2_diff3d
#pragma	noinline	L2_diff2d
#pragma	noinline	L2_regular_integral
#pragma	noinline	subdomain_integral
#pragma	noinline	zero_states
#pragma	noinline	add_product_to_total
#pragma	noinline	add_multiple_of_product_to_total
#pragma	noinline	state_sqrt
#pragma	noinline	absolute_value
#pragma	noinline	multiplied_by_factor
#pragma	noinline	add_to_total
#pragma	noinline	find_diff_and_other_states
#pragma	noinline	Lp_printout
#pragma noinline        is_block_on_front
#endif /*!defined(__INTEL_COMPILER)*/

/*
*			Lp_diff():
*
*       Computes the Lp norm of the difference between the numerical solution
*       computed using the full Euler equations and the approximate solution
*       given by the linearized equations of motion.
*/

/*ARGSUSED*/
EXPORT void Lp_diff(
	Grid		*grid,
	Wave		*wave,
	Front		*fr,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	Lp_Diff_data    *lp_data = Lp_data(data);
	int             dim = fr->rect_grid->dim;
	int             p = lp_data->p;

	if (p == 1)
	    L1_diff(grid,wave,fr,data);
	else if (dim == 3)
	    L2_diff3d(grid,wave,fr,data);
	else if (dim == 2)
	    L2_diff2d(grid,wave,fr,data);
	else
	{
	    screen("L2_diff is not defined for dim = %d\n",dim);
	    clean_up(ERROR);
	}
}		/*end Lp_diff*/

LOCAL  void L1_diff(
	Grid            *grid,
	Wave            *wave,
	Front           *fr,
	OUTPUT_DATA     *data)
{
	COMPONENT       comp;
	INTERFACE       *intfc = fr->interf;
	Lp_Diff_data    *lp_data = Lp_data(data);
	int             iz, iy, ix, i;
	int             dim = fr->rect_grid->dim;
	int             ratio = lp_data->ratio;
	int             gmax[3];
	int		stype = lp_data->stype;
	double           L[3], h[3], dv;
	double           coords[3];
	double           nonlin_pr, lin_pr, diff_pr;
	double           nonlin_vel[3], lin_vel[3], diff_vel[3];
	double           nonlin_vel_total[3], lin_vel_total[3];
	double		diff_vel_total[3];
	double           nonlin_pr_total, lin_pr_total, diff_pr_total;
	static Locstate nonlin_total = NULL, lin_total = NULL,
	                diff_total = NULL;
	static Locstate nonlin = NULL, lin = NULL, diff = NULL;

	for (i = 0; i < dim; i++)
	{
	    L[i] = fr->rect_grid->L[i];
	    h[i] = fr->rect_grid->h[i];
	    h[i] /= ratio;
	    gmax[i] = fr->rect_grid->gmax[i] * ratio;
	}
	for (i = dim; i < 3; i++)
	{
	    L[i] = h[i] = 0.0;
	    gmax[i] = 1;
	}
	if (nonlin == NULL)
	{
	    alloc_state(fr->interf,&nonlin,fr->sizest);
	    alloc_state(fr->interf,&lin,fr->sizest);
	    alloc_state(fr->interf,&diff,fr->sizest);
	    alloc_state(fr->interf,&nonlin_total,fr->sizest);
	    alloc_state(fr->interf,&lin_total,fr->sizest);
	    alloc_state(fr->interf,&diff_total,fr->sizest);
	}

	zero_states(nonlin_total,lin_total,diff_total,&nonlin_pr_total,
	    &lin_pr_total,&diff_pr_total,nonlin_vel_total,lin_vel_total,
	    diff_vel_total,fr);

	for (iz = 0; iz < gmax[2]; iz++)
	{
	    coords[2] = L[2] + (iz + 0.5) * h[2];
	    for (iy = 0; iy < gmax[1]; iy++)
	    {
		coords[1] = L[1] + (iy + 0.5) * h[1];
		for (ix = 0; ix < gmax[0]; ix++)
		{
		    coords[0] = L[0] + (ix + 0.5) * h[0];

		    comp = component(coords,intfc);
		    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				 fr,wave,nonlin,NULL);
		    (*lp_data->_alternate_state)(lin,coords,fr,wave,stype);

		    find_diff_and_other_states(diff,&nonlin_pr,&lin_pr,&diff_pr,
			                       nonlin_vel,lin_vel,diff_vel,
					       nonlin,lin,dim);
		    absolute_value(nonlin,lin,diff,&nonlin_pr,&lin_pr,&diff_pr,
			           nonlin_vel,lin_vel,diff_vel,dim);
		    add_to_total(nonlin_total,lin_total,diff_total,
			         &nonlin_pr_total,&lin_pr_total,&diff_pr_total,
			         nonlin_vel_total,lin_vel_total,diff_vel_total,
			         nonlin,lin,diff,nonlin_pr,lin_pr,diff_pr,
			         nonlin_vel,lin_vel,diff_vel,dim);

		}
	    }
	}

	dv = h[0];
	for (i = 1; i < dim; i++)
	    dv *= h[i];
	multiplied_by_factor(nonlin_total,lin_total,diff_total,
	    &nonlin_pr_total,&lin_pr_total,&diff_pr_total,
	    nonlin_vel_total,lin_vel_total,diff_vel_total,dim,dv);

	pp_global_sum((double*)diff_total,2+dim);
	pp_global_sum(&diff_pr_total,1);
	pp_global_sum(diff_vel_total,dim);
	pp_global_sum((double*)lin_total,2+dim);
	pp_global_sum(&lin_pr_total,1);
	pp_global_sum(lin_vel_total,dim);
	pp_global_sum((double*)nonlin_total,2+dim);
	pp_global_sum(&nonlin_pr_total,1);
	pp_global_sum(nonlin_vel_total,dim);

	Lp_printout(data,grid,diff_total,lin_total,nonlin_total,
	    diff_pr_total,lin_pr_total,nonlin_pr_total,
	    diff_vel_total,lin_vel_total,nonlin_vel_total,dim);
}               /*end L1_diff*/

LOCAL  void L2_diff3d(
	Grid            *grid,
	Wave            *wave,
	Front           *fr,
	OUTPUT_DATA     *data)
{
	Lp_Diff_data	*lp_data = Lp_data(data);
	COMPONENT	comp;
	INTERFACE	*intfc = fr->interf;
	int		stype = lp_data->stype;
	int             ix, iy, iz, i, j;
	int             *gmax = fr->rect_grid->gmax;
	int             dim = 3;
	int             ratio = lp_data->ratio;
	int             npts, **icoords, indx;
	int             i1, i2, i3, i4, j1, j2, j3, j4;
	SINE_PERT	*pert = lp_data->pert;
	double           m1, m2, m3, m4, m5, dv;
	double           L[3], U[3], h[3], coords[3], h_sub[3];
	double           tolerance;
	double           nonlin_pr_total, lin_pr_total, diff_pr_total;
	double           nonlin_pr, lin_pr, diff_pr;
	double           nonlin_vel[3], lin_vel[3], diff_vel[3];
	double           nonlin_vel_total[3], lin_vel_total[3],
	                diff_vel_total[3];
	static double    *nonlin_prb = NULL, *lin_prb = NULL, *diff_prb = NULL;
	static double    **nonlin_velb = NULL, **lin_velb = NULL,
	                **diff_velb = NULL;
	static Locstate nonlin_total = NULL, lin_total = NULL,
	                diff_total = NULL;
	static Locstate *lin_regular_grid_state = NULL;
	static Locstate nonlin = NULL, lin = NULL, diff = NULL;
	static Locstate *nonlinb = NULL, *linb = NULL, *diffb = NULL;
	static double    *z_lin_intfc = NULL;

	bi_array(&icoords,8,3,INT);
	for (i = 0; i < dim; i++)
	{
	    L[i] = fr->rect_grid->L[i];
	    U[i] = fr->rect_grid->U[i];
	    h[i] = fr->rect_grid->h[i];
	    h_sub[i] = h[i] / ratio;
	}
	tolerance = h[2] * 1.0e-6;

	if (nonlin == NULL)
	{
	    alloc_state(fr->interf,&nonlin,fr->sizest);
	    alloc_state(fr->interf,&lin,fr->sizest);
	    alloc_state(fr->interf,&diff,fr->sizest);
	    alloc_state(fr->interf,&nonlin_total,fr->sizest);
	    alloc_state(fr->interf,&lin_total,fr->sizest);
	    alloc_state(fr->interf,&diff_total,fr->sizest);
	    npts = (gmax[0] + 1) * (gmax[1] + 1) * gmax[2];
	    uni_array(&lin_regular_grid_state,npts,sizeof(Locstate));
	    for (i = 0; i < npts; i++)
	        alloc_state(fr->interf,&lin_regular_grid_state[i],fr->sizest);
	    npts = (gmax[0] + 1) * (gmax[1] + 1) * 2;
	    uni_array(&nonlinb,npts,sizeof(Locstate));
	    uni_array(&linb,npts,sizeof(Locstate));
	    uni_array(&diffb,npts,sizeof(Locstate));
	    uni_array(&nonlin_prb,npts,FLOAT);
	    uni_array(&lin_prb,npts,FLOAT);
	    uni_array(&diff_prb,npts,FLOAT);
	    bi_array(&nonlin_velb,npts,3,FLOAT);
	    bi_array(&lin_velb,npts,3,FLOAT);
	    bi_array(&diff_velb,npts,3,FLOAT);
	    for (i = 0; i < npts; i++)
		alloc_state(fr->interf,&diffb[i],fr->sizest);
	    uni_array(&z_lin_intfc,npts,FLOAT);
	    npts /= 2;
	    for (i = 0; i < npts; i++)
	    {
		alloc_state(fr->interf,&nonlinb[i],fr->sizest);
		alloc_state(fr->interf,&linb[i],fr->sizest);
	    }
	}
	zero_states(nonlin_total,lin_total,diff_total,&nonlin_pr_total,
	    &lin_pr_total,&diff_pr_total,nonlin_vel_total,lin_vel_total,
	    diff_vel_total,fr);
	for (iz = 0, i = 0; iz < gmax[2]; iz++)
	{
	    coords[2] = L[2] + (iz + 0.5) * h[2];
	    for (iy = -1; iy < gmax[1]; iy++)
	    {
		coords[1] = L[1] + (iy + 0.5) * h[1];
		for (ix = -1; ix < gmax[0]; ix++)
		{
		    coords[0] = L[0] + (ix + 0.5) * h[0];
		    (*lp_data->_alternate_state)(lin_regular_grid_state[i++],
			                         coords,fr,wave,stype);
		}
	    }
	}
	coords[2] = L[2];
	for (iy = -1, i = 0; iy < gmax[1]; iy++)
	{
	    coords[1] = L[1] + (iy + 0.5) * h[1];
	    for (ix = -1; ix < gmax[0]; ix++)
	    {
		coords[0] = L[0] + (ix + 0.5) * h[0];
		z_lin_intfc[i++] = pert_interface(pert,coords,fr->time,dim);
	    }
	    coords[1] += 0.5 * h[1];
	    for (ix = -1; ix < gmax[0]; ix++)
	    {
		coords[0] = L[0] + (ix + 1) * h[0];
		z_lin_intfc[i++] = pert_interface(pert,coords,fr->time,dim);
	    }
	}
	for (iz = 1; iz < gmax[2]; iz++)
	{
	    icoords[0][2] = icoords[1][2] = iz;
	    icoords[2][2] = icoords[3][2] = iz;
	    icoords[4][2] = icoords[5][2] = iz - 1;
	    icoords[6][2] = icoords[7][2] = iz - 1;
	    coords[2] = L[2] + (iz - 0.5) * h[2];
	    for (iy = 0; iy < gmax[1]; iy++)
	    {
		icoords[0][1] = icoords[1][1] = iy;
		icoords[4][1] = icoords[5][1] = iy;
		icoords[2][1] = icoords[3][1] = iy - 1;
		icoords[6][1] = icoords[7][1] = iy - 1;
		coords[1] = L[1] + (iy - 0.5) * h[1];
		indx = (gmax[0] + 1) * 2 * iy;
		for (ix = 0; ix < gmax[0]; ix++, indx++)
		{
		    icoords[0][0] = icoords[3][0] = ix;
		    icoords[4][0] = icoords[7][0] = ix;
		    icoords[1][0] = icoords[2][0] = ix - 1;
		    icoords[5][0] = icoords[6][0] = ix - 1;
		    coords[0] = L[0] + (ix - 0.5) * h[0];
		    if (is_block_on_front(icoords[0],wave,fr,z_lin_intfc,indx))
		    {
			subdomain_integral(diff,nonlin,lin,&diff_pr,&nonlin_pr,
			                   &lin_pr,diff_vel,nonlin_vel,
					   lin_vel,coords,h_sub,wave,fr,data);
		    }
		    else
		    {
			L2_regular_integral(diff,nonlin,lin,&diff_pr,&nonlin_pr,
			                    &lin_pr,diff_vel,nonlin_vel,
					    lin_vel,icoords,wave,fr,
			                    lin_regular_grid_state);
		    }

		    add_to_total(nonlin_total,lin_total,diff_total,
			         &nonlin_pr_total,&lin_pr_total,&diff_pr_total,
			         nonlin_vel_total,lin_vel_total,diff_vel_total,
			         nonlin,lin,diff,nonlin_pr,lin_pr,diff_pr,
			         nonlin_vel,lin_vel,diff_vel,dim);
		}
	    }
	}

	npts = (gmax[0] + 1) * (gmax[1] + 1);
	m1 = 2.0 / 27.0;        m2 = 1.0 / 27.0;
	m3 = m2 * 0.25;         m4 = m2 * 0.5;
	m5 = m2 * 0.125;
	for (iz = 0; iz < gmax[2]; iz += gmax[2]-1)
	{
	    icoords[0][2] = iz;
	    coords[2] = (iz == 0) ? L[2] + tolerance : U[2] - tolerance;

	    for (iy = -1, i = 0; iy < gmax[1]; iy++)
	    {
		icoords[0][1] = iy;
		coords[1] = L[1] + (iy + 0.5) * h[1];
		for (ix = -1; ix < gmax[0]; ix++, i++)
		{
		    icoords[0][0] = ix;
		    coords[0] = L[0] + (ix + 0.5) * h[0];
		    comp = component(coords,intfc);
		    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				 fr,wave,nonlinb[i],NULL);
		    (*lp_data->_alternate_state)(linb[i],coords,fr,wave,stype);
		    find_diff_and_other_states(diffb[i],nonlin_prb+i,
			lin_prb+i,diff_prb+i,nonlin_velb[i],
			lin_velb[i],diff_velb[i],nonlinb[i],linb[i],dim);
		    j = i + npts;
		    nonlinb[j] = Rect_state(icoords[0],wave);
		    indx = (icoords[0][2] * (gmax[1] + 1) + icoords[0][1] + 1)
			* (gmax[0] + 1) + icoords[0][0] + 1;
		    linb[j] = lin_regular_grid_state[indx];
		    find_diff_and_other_states(diffb[j],nonlin_prb+j,
			                       lin_prb+j,diff_prb+j,
					       nonlin_velb[j],lin_velb[j],
					       diff_velb[j],nonlinb[j],
					       linb[j],dim);
		}
	    }

	    for (iy = 0; iy < gmax[1]; iy++)
	    {
		i = (iy + 1) * (gmax[0] + 1) + 1;
		for (ix = 0; ix < gmax[0]; ix++, i++)
		{
		    j = i + npts;
		    i1 = i - 1;         j1 = j - 1;
		    i2 = i - gmax[0] -1;        j2 = j - gmax[0] -1;
		    i3 = i2 - 1;                j3 = j2 - 1;
		    i4 = i2 + 1;                j4 = j2 + 1;

		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[i],
		                                     linb[i],linb[i],
						     diffb[i],diffb[i],
						     nonlin_prb[i],
		                                     nonlin_prb[i],lin_prb[i],
						     lin_prb[i],diff_prb[i],
		                                     diff_prb[i],nonlin_velb[i],
						     nonlin_velb[i],
		                                     lin_velb[i],lin_velb[i],
						     diff_velb[i],diff_velb[i],
						     fr,m1);
		     add_multiple_of_product_to_total(nonlin_total,lin_total,
						      diff_total,
						      &nonlin_pr_total,
						      &lin_pr_total,
		                                      &diff_pr_total,
						      nonlin_vel_total,
						      lin_vel_total,
		                                      diff_vel_total,
						      nonlinb[j],nonlinb[j],
		                                      linb[j],linb[j],
						      diffb[j],diffb[j],
						      nonlin_prb[j],
		                                      nonlin_prb[j],lin_prb[j],
						      lin_prb[j],diff_prb[j],
		                                      diff_prb[j],
						      nonlin_velb[j],
						      nonlin_velb[j],
		                                      lin_velb[j],lin_velb[j],
						      diff_velb[j],
						      diff_velb[j],fr,m1);

		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[i1],
		                                     linb[i],linb[i1],diffb[i],
						     diffb[i1],nonlin_prb[i],
		                                     nonlin_prb[i1],lin_prb[i],
						     lin_prb[i1],diff_prb[i],
		                                     diff_prb[i1],
						     nonlin_velb[i],
						     nonlin_velb[i1],
		                                     lin_velb[i],lin_velb[i1],
						     diff_velb[i],
						     diff_velb[i1],fr,m2);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[j],nonlinb[j1],
		                                     linb[j],linb[j1],diffb[j],
						     diffb[j1],nonlin_prb[j],
		                                     nonlin_prb[j1],lin_prb[j],
						     lin_prb[j1],diff_prb[j],
		                                     diff_prb[j1],
						     nonlin_velb[j],
						     nonlin_velb[j1],
		                                     lin_velb[j],lin_velb[j1],
						     diff_velb[j],
						     diff_velb[j1],fr,m2);

		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[i2],
		                                     linb[i],linb[i2],diffb[i],
						     diffb[i2],nonlin_prb[i],
		                                     nonlin_prb[i2],lin_prb[i],
						     lin_prb[i2],diff_prb[i],
		                                     diff_prb[i2],
						     nonlin_velb[i],
						     nonlin_velb[i2],
		                                     lin_velb[i],lin_velb[i2],
						     diff_velb[i],
						     diff_velb[i2],fr,m2);
	            add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[j],nonlinb[j2],
		                                     linb[j],linb[j2],diffb[j],
						     diffb[j2],nonlin_prb[j],
		                                     nonlin_prb[j2],lin_prb[j],
						     lin_prb[j2],diff_prb[j],
		                                     diff_prb[j2],
						     nonlin_velb[j],
						     nonlin_velb[j2],
		                                     lin_velb[j],lin_velb[j2],
						     diff_velb[j],diff_velb[j2],
						     fr,m2);

		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[j],
		                                     linb[i],linb[j],
						     diffb[i],diffb[j],
						     nonlin_prb[i],
		                                     nonlin_prb[j],
						     lin_prb[i],lin_prb[j],
						     diff_prb[i],
		                                     diff_prb[j],
						     nonlin_velb[i],nonlin_velb[j],
		                                     lin_velb[i],lin_velb[j],
						     diff_velb[i],diff_velb[j],
						     fr,m1);
		
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[i3],
		                                     linb[i],linb[i3],diffb[i],
						     diffb[i3],nonlin_prb[i],
		                                     nonlin_prb[i3],lin_prb[i],
						     lin_prb[i3],diff_prb[i],
		                                     diff_prb[i3],
						     nonlin_velb[i],
						     nonlin_velb[i3],
		                                     lin_velb[i],lin_velb[i3],
						     diff_velb[i],
						     diff_velb[i3],fr,m3);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,nonlinb[j],
						     nonlinb[j3],
		                                     linb[j],linb[j3],diffb[j],
						     diffb[j3],nonlin_prb[j],
		                                     nonlin_prb[j3],lin_prb[j],
						     lin_prb[j3],diff_prb[j],
		                                     diff_prb[j3],
						     nonlin_velb[j],
						     nonlin_velb[j3],
		                                     lin_velb[j],lin_velb[j3],
						     diff_velb[j],
						     diff_velb[j3],fr,m3);

		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[i4],
		                                     linb[i],linb[i4],diffb[i],
						     diffb[i4],nonlin_prb[i],
		                                     nonlin_prb[i4],lin_prb[i],
						     lin_prb[i4],diff_prb[i],
		                                     diff_prb[i4],
						     nonlin_velb[i],
						     nonlin_velb[i4],
		                                     lin_velb[i],lin_velb[i4],
						     diff_velb[i],
						     diff_velb[i4],fr,m3);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,nonlinb[j],
						     nonlinb[j4],
		                                     linb[j],linb[j4],diffb[j],
						     diffb[j4],nonlin_prb[j],
		                                     nonlin_prb[j4],lin_prb[j],
						     lin_prb[j4],diff_prb[j],
		                                     diff_prb[j4],
						     nonlin_velb[j],
						     nonlin_velb[j4],
		                                     lin_velb[j],lin_velb[j4],
						     diff_velb[j],
						     diff_velb[j4],fr,m3);

		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i1],nonlinb[j],
		                                     linb[i1],linb[j],
						     diffb[i1],diffb[j],
						     nonlin_prb[i1],
		                                     nonlin_prb[j],lin_prb[i1],
						     lin_prb[j],diff_prb[i1],
		                                     diff_prb[j],
						     nonlin_velb[i1],
						     nonlin_velb[j],
		                                     lin_velb[i1],lin_velb[j],
						     diff_velb[i1],
						     diff_velb[j],fr,m4);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[j1],
		                                     linb[i],linb[j1],diffb[i],
						     diffb[j1],nonlin_prb[i],
		                                     nonlin_prb[j1],lin_prb[i],
						     lin_prb[j1],diff_prb[i],
		                                     diff_prb[j1],
						     nonlin_velb[i],
						     nonlin_velb[j1],
		                                     lin_velb[i],lin_velb[j1],
						     diff_velb[i],
						     diff_velb[j1],fr,m4);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i2],nonlinb[j],
		                                     linb[i2],linb[j],diffb[i2],
						     diffb[j],nonlin_prb[i2],
		                                     nonlin_prb[j],lin_prb[i2],
						     lin_prb[j],diff_prb[i2],
		                                     diff_prb[j],
						     nonlin_velb[i2],
						     nonlin_velb[j],
		                                     lin_velb[i2],lin_velb[j],
						     diff_velb[i2],
						     diff_velb[j],fr,m4);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,nonlinb[i],
						     nonlinb[j2],
		                                     linb[i],linb[j2],diffb[i],
						     diffb[j2],nonlin_prb[i],
		                                     nonlin_prb[j2],lin_prb[i],
						     lin_prb[j2],diff_prb[i],
		                                     diff_prb[j2],
						     nonlin_velb[i],
						     nonlin_velb[j2],
		                                     lin_velb[i],lin_velb[j2],
						     diff_velb[i],
						     diff_velb[j2],fr,m4);

		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i3],nonlinb[j],
		                                     linb[i3],linb[j],diffb[i3],
						     diffb[j],nonlin_prb[i3],
		                                     nonlin_prb[j],lin_prb[i3],
						     lin_prb[j],diff_prb[i3],
		                                     diff_prb[j],
						     nonlin_velb[i3],
						     nonlin_velb[j],
		                                     lin_velb[i3],lin_velb[j],
						     diff_velb[i3],
						     diff_velb[j],fr,m5);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[j3],
		                                     linb[i],linb[j3],diffb[i],
						     diffb[j3],nonlin_prb[i],
		                                     nonlin_prb[j3],lin_prb[i],
						     lin_prb[j3],diff_prb[i],
		                                     diff_prb[j3],
						     nonlin_velb[i],
						     nonlin_velb[j3],
		                                     lin_velb[i],lin_velb[j3],
						     diff_velb[i],
						     diff_velb[j3],fr,m5);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,nonlinb[i4],
						     nonlinb[j],linb[i4],
						     linb[j],diffb[i4],
						     diffb[j],nonlin_prb[i4],
		                                     nonlin_prb[j],lin_prb[i4],
						     lin_prb[j],diff_prb[i4],
		                                     diff_prb[j],
						     nonlin_velb[i4],
						     nonlin_velb[j],
		                                     lin_velb[i4],lin_velb[j],
						     diff_velb[i4],
						     diff_velb[j],fr,m5);
		    add_multiple_of_product_to_total(nonlin_total,lin_total,
						     diff_total,
						     &nonlin_pr_total,
						     &lin_pr_total,
		                                     &diff_pr_total,
						     nonlin_vel_total,
						     lin_vel_total,
		                                     diff_vel_total,
						     nonlinb[i],nonlinb[j4],
		                                     linb[i],linb[j4],diffb[i],
						     diffb[j4],nonlin_prb[i],
		                                     nonlin_prb[j4],lin_prb[i],
						     lin_prb[j4],diff_prb[i],
		                                     diff_prb[j4],
						     nonlin_velb[i],
						     nonlin_velb[j4],
		                                     lin_velb[i],lin_velb[j4],
						     diff_velb[i],
						     diff_velb[j4],fr,m5);
	        }
	    }
	}
	dv = h[0] * h[1] * h[2];
	multiplied_by_factor(nonlin_total,lin_total,diff_total,&nonlin_pr_total,
		             &lin_pr_total,&diff_pr_total,nonlin_vel_total,
			     lin_vel_total,diff_vel_total,dim,dv);

	pp_global_sum((double*)diff_total,2+dim);
	pp_global_sum(&diff_pr_total,1);
	pp_global_sum(diff_vel_total,dim);
	pp_global_sum((double*)lin_total,2+dim);
	pp_global_sum(&lin_pr_total,1);
	pp_global_sum(lin_vel_total,dim);
	pp_global_sum((double*)nonlin_total,2+dim);
	pp_global_sum(&nonlin_pr_total,1);
	pp_global_sum(nonlin_vel_total,dim);

	state_sqrt(nonlin_total,lin_total,diff_total,&nonlin_pr_total,
	    &lin_pr_total,&diff_pr_total,nonlin_vel_total,
	    lin_vel_total,diff_vel_total,dim);

	Lp_printout(data,grid,diff_total,lin_total,nonlin_total,
		diff_pr_total,lin_pr_total,nonlin_pr_total,
		diff_vel_total,lin_vel_total,nonlin_vel_total,dim);
}		/*end L2_diff3d*/

LOCAL  void L2_diff2d(
	Grid            *grid,
	Wave            *wave,
	Front           *fr,
	OUTPUT_DATA     *data)
{
	Lp_Diff_data	*lp_data = Lp_data(data);
	COMPONENT	comp;
	INTERFACE	*intfc = fr->interf;
	SINE_PERT	*pert = lp_data->pert;
	int		stype = lp_data->stype;
	int		ix, iy, i, j; 
	int		*gmax = fr->rect_grid->gmax;
	int		dim = 2;
	int		ratio = lp_data->ratio;
	int             npts, **icoords, indx;
	int		i1, j1;
	double		m1, m2, m3, dv;
	double		L[2], U[2], h[2], coords[3], h_sub[3];
	double		tolerance;
	double		lin_vel_total[2];
	double		nonlin_vel_total[2];
	double		diff_vel_total[2];
	double		nonlin_pr, lin_pr, diff_pr;
	double		nonlin_vel[2], lin_vel[2], diff_vel[2];
	double		nonlin_pr_total, lin_pr_total, diff_pr_total;
	static double	*nonlin_prb = NULL, *lin_prb = NULL, *diff_prb = NULL;
	static double	**nonlin_velb = NULL, **lin_velb = NULL,
	                **diff_velb = NULL;
	static Locstate *lin_regular_grid_state = NULL;
	static Locstate nonlin = NULL, lin = NULL, diff = NULL;
	static Locstate	nonlin_total = NULL, lin_total = NULL,
	                diff_total = NULL;
	static Locstate *nonlinb = NULL, *linb = NULL, *diffb = NULL;
	static double    *y_lin_intfc = NULL;

	bi_array(&icoords,4,2,INT);
	for (i = 0; i < dim; i++)
	{
	    L[i] = fr->rect_grid->L[i];
	    U[i] = fr->rect_grid->U[i];
	    h[i] = fr->rect_grid->h[i];
	    h_sub[i] = h[i] / ratio;
	}
	coords[2] = 0.0;
	h_sub[2] = 0.0;
	tolerance = h[1] * 1.0e-6;

	if (nonlin == NULL)
	{
	    alloc_state(fr->interf,&nonlin,fr->sizest);
	    alloc_state(fr->interf,&lin,fr->sizest);
	    alloc_state(fr->interf,&diff,fr->sizest);
	    alloc_state(fr->interf,&nonlin_total,fr->sizest);
	    alloc_state(fr->interf,&lin_total,fr->sizest);
	    alloc_state(fr->interf,&diff_total,fr->sizest);
	    npts = (gmax[0] + 1) * gmax[1];
	    uni_array(&lin_regular_grid_state,npts,sizeof(Locstate));
	    for (i = 0; i < npts; i++)
	    	alloc_state(fr->interf,&lin_regular_grid_state[i],fr->sizest);
	    npts = (gmax[0] + 1) * 2;
	    uni_array(&nonlinb,npts,sizeof(Locstate));
	    uni_array(&linb,npts,sizeof(Locstate));
	    uni_array(&diffb,npts,sizeof(Locstate));
	    uni_array(&nonlin_prb,npts,FLOAT);
	    uni_array(&lin_prb,npts,FLOAT);
	    uni_array(&diff_prb,npts,FLOAT);
	    bi_array(&nonlin_velb,npts,2,FLOAT);
	    bi_array(&lin_velb,npts,2,FLOAT);
	    bi_array(&diff_velb,npts,2,FLOAT);
	    for (i = 0; i < npts; i++)
	    {
		alloc_state(fr->interf,&diffb[i],fr->sizest);
	    }
	    uni_array(&y_lin_intfc,npts,FLOAT);
	    npts = gmax[0] + 1;
	    for (i = 0; i < npts; i++)
	    {
		alloc_state(fr->interf,&nonlinb[i],fr->sizest);
		alloc_state(fr->interf,&linb[i],fr->sizest);
	    }
	}
	zero_states(nonlin_total,lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		&diff_pr_total,nonlin_vel_total,lin_vel_total,diff_vel_total,fr);

	for (iy = 0, i = 0; iy < gmax[1]; iy++) 
	{
	    coords[1] = L[1] + (iy + 0.5) * h[1];
	    for (ix = -1; ix < gmax[0]; ix++)
	    {
	        coords[0] = L[0] + (ix + 0.5) * h[0];
		(*lp_data->_alternate_state)(lin_regular_grid_state[i++],
					     coords,fr,wave,stype);
	    }
	}
	coords[1] = L[1];
	for (ix = -1, i = 0; ix < gmax[0]; ix++)
	{
	    coords[0] = L[0] + (ix + 0.5) * h[0];
	    y_lin_intfc[i++] = pert_interface(pert,coords,fr->time,dim);
	    coords[0] = L[0] + (ix + 1) * h[0];
	    y_lin_intfc[i++] = pert_interface(pert,coords,fr->time,dim);
	}
	for (iy = 1; iy < gmax[1]; iy++) 
	{
	    icoords[0][1] = icoords[1][1] = iy;
	    icoords[2][1] = icoords[3][1] = iy - 1;
	    coords[1] = L[1] + (iy - 0.5) * h[1];
	    for (ix = 0; ix < gmax[0]; ix++) 
	    {
	        icoords[0][0] = icoords[3][0] = ix;
		icoords[1][0] = icoords[2][0] = ix - 1;
		coords[0] = L[0] + (ix - 0.5) * h[0];

		if (is_block_on_front(icoords[0],wave,fr,y_lin_intfc,ix*2))
		    subdomain_integral(diff,nonlin,lin,&diff_pr,&nonlin_pr,
				       &lin_pr,diff_vel,nonlin_vel,lin_vel,
				       coords,h_sub,wave,fr,data);
		else
		    L2_regular_integral(diff,nonlin,lin,&diff_pr,&nonlin_pr,
					&lin_pr,diff_vel,nonlin_vel,lin_vel,
					icoords,wave,fr,lin_regular_grid_state);

		add_to_total(nonlin_total,lin_total,diff_total,
			     &nonlin_pr_total,&lin_pr_total,&diff_pr_total,
			     nonlin_vel_total,lin_vel_total,diff_vel_total,
			     nonlin,lin,diff,nonlin_pr,lin_pr,diff_pr,
			     nonlin_vel,lin_vel,diff_vel,dim);
	    }
	}

	npts = gmax[0] + 1;
	m1 = 1.0 / 9.0;		m2 = m1 * 0.5;		m3 = m2 * 0.5;

	for (iy = 0; iy < gmax[1]; iy += gmax[1] - 1)
	{
	    icoords[0][1] = iy;
	    coords[1] = (iy == 0) ? L[1] + tolerance : U[1] - tolerance;

	    for (ix = -1, i = 0; ix < gmax[0]; ix++, i++)
	    {
		icoords[0][0] = ix;
		coords[0] = L[0] + (ix + 0.5) * h[0];
		comp = component(coords,intfc);
		hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
			     fr,wave,nonlinb[i],NULL);
		(*lp_data->_alternate_state)(linb[i],coords,fr,wave,stype);
		find_diff_and_other_states(diffb[i],nonlin_prb+i,
		    lin_prb+i,diff_prb+i,nonlin_velb[i],
		    lin_velb[i],diff_velb[i],nonlinb[i],linb[i],dim);
		j = i + npts;
		nonlinb[j] = Rect_state(icoords[0],wave);
		indx = icoords[0][1] * (gmax[0] + 1) + icoords[0][0] + 1;
		linb[j] = lin_regular_grid_state[indx];
		find_diff_and_other_states(diffb[j],nonlin_prb+j,
		    lin_prb+j,diff_prb+j,nonlin_velb[j],
		    lin_velb[j],diff_velb[j],nonlinb[j],linb[j],dim);
	    }

	    for (ix = 0, i = 1; ix < gmax[0]; ix++, i++)
	    {
		j = i + npts;
		i1 = i - 1;		j1 = j - 1;

		add_multiple_of_product_to_total(nonlin_total,
		    lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		    &diff_pr_total,nonlin_vel_total,lin_vel_total,
		    diff_vel_total,nonlinb[i],nonlinb[i],
		    linb[i],linb[i],diffb[i],diffb[i],nonlin_prb[i],
		    nonlin_prb[i],lin_prb[i],lin_prb[i],diff_prb[i],
		    diff_prb[i],nonlin_velb[i],nonlin_velb[i],
		    lin_velb[i],lin_velb[i],diff_velb[i],diff_velb[i],fr,m1);
		add_multiple_of_product_to_total(nonlin_total,
		    lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		    &diff_pr_total,nonlin_vel_total,lin_vel_total,
		    diff_vel_total,nonlinb[j],nonlinb[j],
		    linb[j],linb[j],diffb[j],diffb[j],nonlin_prb[j],
		    nonlin_prb[j],lin_prb[j],lin_prb[j],diff_prb[j],
		    diff_prb[j],nonlin_velb[j],nonlin_velb[j],
		    lin_velb[j],lin_velb[j],diff_velb[j],diff_velb[j],fr,m1);

		add_multiple_of_product_to_total(nonlin_total,
		    lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		    &diff_pr_total,nonlin_vel_total,lin_vel_total,
		    diff_vel_total,nonlinb[i],nonlinb[i1],
		    linb[i],linb[i1],diffb[i],diffb[i1],nonlin_prb[i],
		    nonlin_prb[i1],lin_prb[i],lin_prb[i1],diff_prb[i],
		    diff_prb[i1],nonlin_velb[i],nonlin_velb[i1],
		    lin_velb[i],lin_velb[i1],diff_velb[i],diff_velb[i1],fr,m2);
		add_multiple_of_product_to_total(nonlin_total,
		    lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		    &diff_pr_total,nonlin_vel_total,lin_vel_total,
		    diff_vel_total,nonlinb[j],nonlinb[j1],
		    linb[j],linb[j1],diffb[j],diffb[j1],nonlin_prb[j],
		    nonlin_prb[j1],lin_prb[j],lin_prb[j1],diff_prb[j],
		    diff_prb[j1],nonlin_velb[j],nonlin_velb[j1],
		    lin_velb[j],lin_velb[j1],diff_velb[j],diff_velb[j1],fr,m2);

		add_multiple_of_product_to_total(nonlin_total,
		    lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		    &diff_pr_total,nonlin_vel_total,lin_vel_total,
		    diff_vel_total,nonlinb[i],nonlinb[j],
		    linb[i],linb[j],diffb[i],diffb[j],nonlin_prb[i],
		    nonlin_prb[j],lin_prb[i],lin_prb[j],diff_prb[i],
		    diff_prb[j],nonlin_velb[i],nonlin_velb[j],
		    lin_velb[i],lin_velb[j],diff_velb[i],diff_velb[j],fr,m1);

		add_multiple_of_product_to_total(nonlin_total,
		    lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		    &diff_pr_total,nonlin_vel_total,lin_vel_total,
		    diff_vel_total,nonlinb[i1],nonlinb[j],
		    linb[i1],linb[j],diffb[i1],diffb[j],nonlin_prb[i1],
		    nonlin_prb[j],lin_prb[i1],lin_prb[j],diff_prb[i1],
		    diff_prb[j],nonlin_velb[i1],nonlin_velb[j],
		    lin_velb[i1],lin_velb[j],diff_velb[i1],diff_velb[j],fr,m3);

		add_multiple_of_product_to_total(nonlin_total,
		    lin_total,diff_total,&nonlin_pr_total,&lin_pr_total,
		    &diff_pr_total,nonlin_vel_total,lin_vel_total,
		    diff_vel_total,nonlinb[i],nonlinb[j1],
		    linb[i],linb[j1],diffb[i],diffb[j1],nonlin_prb[i],
		    nonlin_prb[j1],lin_prb[i],lin_prb[j1],diff_prb[i],
		    diff_prb[j1],nonlin_velb[i],nonlin_velb[j1],
		    lin_velb[i],lin_velb[j1],diff_velb[i],diff_velb[j1],fr,m3);
	    }	
	}
	dv = h[0] * h[1];
	multiplied_by_factor(nonlin_total,lin_total,diff_total,&nonlin_pr_total,
		&lin_pr_total,&diff_pr_total,nonlin_vel_total,lin_vel_total,
		diff_vel_total,dim,dv);

	pp_global_sum((double*)diff_total,2+dim);
	pp_global_sum(&diff_pr_total,1);
	pp_global_sum(diff_vel_total,dim);
	pp_global_sum((double*)lin_total,2+dim);
	pp_global_sum(&lin_pr_total,1);
	pp_global_sum(lin_vel_total,dim);
	pp_global_sum((double*)nonlin_total,2+dim);
	pp_global_sum(&nonlin_pr_total,1);
	pp_global_sum(nonlin_vel_total,dim);

	state_sqrt(nonlin_total,lin_total,diff_total,&nonlin_pr_total,
		   &lin_pr_total,&diff_pr_total,nonlin_vel_total,
		   lin_vel_total,diff_vel_total,dim);

	Lp_printout(data,grid,diff_total,lin_total,nonlin_total,
		    diff_pr_total,lin_pr_total,nonlin_pr_total,
		    diff_vel_total,lin_vel_total,nonlin_vel_total,dim);
}		/*end L2_diff2d*/

LOCAL  void L2_regular_integral(
        Locstate	diff_total,
        Locstate	nonlin_total,
        Locstate	lin_total,
        double		*diff_pr_total,
        double 		*nonlin_pr_total,
        double		*lin_pr_total,
        double		*diff_vel_total,
        double		*nonlin_vel_total,
        double		*lin_vel_total,
	int		**icoords,
        Wave 		*wave,
	Front          	*fr,
	Locstate	*lin_regular_grid_state)
{
        int	       	i, j, indx, npts;
	int             *gmax = fr->rect_grid->gmax;
	int            	dim = fr->rect_grid->dim;
	double		factor;
	double		nonlin_pr[8], lin_pr[8], diff_pr[8];
	static double	*nonlin_vel[8], *lin_vel[8], *diff_vel[8];
        static Locstate	*nonlin = NULL, *lin = NULL, *diff = NULL;
	static Locstate nonlin_tmp = NULL, lin_tmp = NULL, diff_tmp = NULL;
	double		nonlin_pr_tmp, lin_pr_tmp, diff_pr_tmp;
	double		nonlin_vel_tmp[3], lin_vel_tmp[3], diff_vel_tmp[3];

	npts = (dim == 3) ? 8 : 4;
        if (nonlin == NULL)
        {
	    uni_array(&nonlin,npts,sizeof(Locstate));
	    uni_array(&lin,npts,sizeof(Locstate));
	    uni_array(&diff,npts,sizeof(Locstate));
	    for (i = 0; i < npts; i++)
	    {
		alloc_state(fr->interf,&diff[i],fr->sizest);
		uni_array(&nonlin_vel[i],3,FLOAT);
		uni_array(&lin_vel[i],3,FLOAT);
		uni_array(&diff_vel[i],3,FLOAT);
	    }
	    alloc_state(fr->interf,&nonlin_tmp,fr->sizest);
	    alloc_state(fr->interf,&lin_tmp,fr->sizest);
	    alloc_state(fr->interf,&diff_tmp,fr->sizest);
        }

	zero_states(nonlin_total,lin_total,diff_total,nonlin_pr_total,lin_pr_total,
		diff_pr_total,nonlin_vel_total,lin_vel_total,diff_vel_total,fr);
	for (i = 0; i < npts; i++)
	{
	    nonlin[i] = Rect_state(icoords[i],wave);
	    if (dim == 3)
	       indx = (icoords[i][2] * (gmax[1] + 1) + icoords[i][1] + 1)
			* (gmax[0] + 1) + icoords[i][0] + 1;
	    else
	       indx = icoords[i][1] * (gmax[0] + 1) + icoords[i][0] + 1;
	    lin[i] = lin_regular_grid_state[indx];
	    find_diff_and_other_states(diff[i],nonlin_pr+i,lin_pr+i,diff_pr+i,
	       nonlin_vel[i],lin_vel[i],diff_vel[i],nonlin[i],lin[i],dim);
	}

	for (i = 0; i < npts; i++)
	    add_product_to_total(nonlin_total,lin_total,diff_total,
		nonlin_pr_total,lin_pr_total,diff_pr_total,
		nonlin_vel_total,lin_vel_total,diff_vel_total,
		dim,nonlin[i],nonlin[i],lin[i],lin[i],diff[i],diff[i],
		nonlin_pr[i],nonlin_pr[i],lin_pr[i],lin_pr[i],
		diff_pr[i],diff_pr[i],nonlin_vel[i],nonlin_vel[i],
		lin_vel[i],lin_vel[i],diff_vel[i],diff_vel[i]);
	for (i = 0; i < 4; i++)
	{
	    j = (i + 1) % 4;
	    add_product_to_total(nonlin_total,lin_total,diff_total,
		nonlin_pr_total,lin_pr_total,diff_pr_total,
		nonlin_vel_total,lin_vel_total,diff_vel_total,
		dim,nonlin[i],nonlin[j],lin[i],lin[j],diff[i],diff[j],
		nonlin_pr[i],nonlin_pr[j],lin_pr[i],lin_pr[j],
		diff_pr[i],diff_pr[j],nonlin_vel[i],nonlin_vel[j],
		lin_vel[i],lin_vel[j],diff_vel[i],diff_vel[j]);
	     if (dim == 3)
	     {
		add_product_to_total(nonlin_total,lin_total,diff_total,
			nonlin_pr_total,lin_pr_total,diff_pr_total,
			nonlin_vel_total,lin_vel_total,diff_vel_total,
			dim,nonlin[i+4],nonlin[j+4],lin[i+4],lin[j+4],diff[i+4],diff[j+4],
			nonlin_pr[i+4],nonlin_pr[j+4],lin_pr[i+4],lin_pr[j+4],
			diff_pr[i+4],diff_pr[j+4],nonlin_vel[i+4],nonlin_vel[j+4],
			lin_vel[i+4],lin_vel[j+4],diff_vel[i+4],diff_vel[j+4]);
		add_product_to_total(nonlin_total,lin_total,diff_total,
			nonlin_pr_total,lin_pr_total,diff_pr_total,
			nonlin_vel_total,lin_vel_total,diff_vel_total,
			dim,nonlin[i],nonlin[i+4],lin[i],lin[i+4],diff[i],diff[i+4],
			nonlin_pr[i],nonlin_pr[i+4],lin_pr[i],lin_pr[i+4],
			diff_pr[i],diff_pr[i+4],nonlin_vel[i],nonlin_vel[i+4],
			lin_vel[i],lin_vel[i+4],diff_vel[i],diff_vel[i+4]);
	     }
	}
	zero_states(nonlin_tmp,lin_tmp,diff_tmp,&nonlin_pr_tmp,&lin_pr_tmp,
		&diff_pr_tmp,nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,fr);
	for (i = 0; i < 2; i++)
	{
	    add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		&nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		dim,nonlin[i],nonlin[i+2],lin[i],lin[i+2],diff[i],diff[i+2],
		nonlin_pr[i],nonlin_pr[i+2],lin_pr[i],lin_pr[i+2],
		diff_pr[i],diff_pr[i+2],nonlin_vel[i],nonlin_vel[i+2],
		lin_vel[i],lin_vel[i+2],diff_vel[i],diff_vel[i+2]);
	    if (dim == 3)
	    {
		j = 2 * i;
		add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		    &nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		    dim,nonlin[i+4],nonlin[i+6],lin[i+4],lin[i+6],diff[i+4],diff[i+6],
		    nonlin_pr[i+4],nonlin_pr[i+6],lin_pr[i+4],lin_pr[i+6],
		    diff_pr[i+4],diff_pr[i+6],nonlin_vel[i+4],nonlin_vel[i+6],
		    lin_vel[i+4],lin_vel[i+6],diff_vel[i+4],diff_vel[i+6]);
		add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		    &nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		    dim,nonlin[i],nonlin[i+5],lin[i],lin[i+5],diff[i],diff[i+5],
		    nonlin_pr[i],nonlin_pr[i+5],lin_pr[i],lin_pr[i+5],
		    diff_pr[i],diff_pr[i+5],nonlin_vel[i],nonlin_vel[i+5],
		    lin_vel[i],lin_vel[i+5],diff_vel[i],diff_vel[i+5]);
		add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		    &nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		    dim,nonlin[i+1],nonlin[i+4],lin[i+1],lin[i+4],diff[i+1],diff[i+4],
		    nonlin_pr[i+1],nonlin_pr[i+4],lin_pr[i+1],lin_pr[i+4],
		    diff_pr[i+1],diff_pr[i+4],nonlin_vel[i+1],nonlin_vel[i+4],
		    lin_vel[i+1],lin_vel[i+4],diff_vel[i+1],diff_vel[i+4]);
		add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		    &nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		    dim,nonlin[j],nonlin[7],lin[j],lin[7],diff[j],diff[7],
		    nonlin_pr[j],nonlin_pr[7],lin_pr[j],lin_pr[7],
		    diff_pr[j],diff_pr[7],nonlin_vel[j],nonlin_vel[7],
		    lin_vel[j],lin_vel[7],diff_vel[j],diff_vel[7]);
		add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		    &nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		    dim,nonlin[j+4],nonlin[3],lin[j+4],lin[3],diff[j+4],diff[3],
		    nonlin_pr[j+4],nonlin_pr[3],lin_pr[j+4],lin_pr[3],
		    diff_pr[j+4],diff_pr[3],nonlin_vel[j+4],nonlin_vel[3],
		    lin_vel[j+4],lin_vel[3],diff_vel[j+4],diff_vel[3]);
	    }
	}
	multiplied_by_factor(nonlin_tmp,lin_tmp,diff_tmp,&nonlin_pr_tmp,&lin_pr_tmp,
		&diff_pr_tmp,nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,dim,0.5);
	add_to_total(nonlin_total,lin_total,diff_total,
		nonlin_pr_total,lin_pr_total,diff_pr_total,
		nonlin_vel_total,lin_vel_total,diff_vel_total,
		nonlin_tmp,lin_tmp,diff_tmp,nonlin_pr_tmp,lin_pr_tmp,diff_pr_tmp,
		nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,dim);
	if (dim == 3)
	{
	    zero_states(nonlin_tmp,lin_tmp,diff_tmp,&nonlin_pr_tmp,&lin_pr_tmp,
	        &diff_pr_tmp,nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,fr);
	    for (i = 0; i < 2; i++)
	    {
		add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		    &nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		    dim,nonlin[i],nonlin[i+6],lin[i],lin[i+6],diff[i],diff[i+6],
		    nonlin_pr[i],nonlin_pr[i+6],lin_pr[i],lin_pr[i+6],
		    diff_pr[i],diff_pr[i+6],nonlin_vel[i],nonlin_vel[i+6],
		    lin_vel[i],lin_vel[i+6],diff_vel[i],diff_vel[i+6]);
		add_product_to_total(nonlin_tmp,lin_tmp,diff_tmp,
		    &nonlin_pr_tmp,&lin_pr_tmp,&diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,
		    dim,nonlin[i+2],nonlin[i+4],lin[i+2],lin[i+4],diff[i+2],diff[i+4],
		    nonlin_pr[i+2],nonlin_pr[i+4],lin_pr[i+2],lin_pr[i+4],
		    diff_pr[i+2],diff_pr[i+4],nonlin_vel[i+2],nonlin_vel[i+4],
		    lin_vel[i+2],lin_vel[i+4],diff_vel[i+2],diff_vel[i+4]);
	    }
	    multiplied_by_factor(nonlin_tmp,lin_tmp,diff_tmp,&nonlin_pr_tmp,&lin_pr_tmp,
		    &diff_pr_tmp,nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,dim,0.25);
	    add_to_total(nonlin_total,lin_total,diff_total,
	    	    nonlin_pr_total,lin_pr_total,diff_pr_total,
		    nonlin_vel_total,lin_vel_total,diff_vel_total,
		    nonlin_tmp,lin_tmp,diff_tmp,nonlin_pr_tmp,lin_pr_tmp,diff_pr_tmp,
		    nonlin_vel_tmp,lin_vel_tmp,diff_vel_tmp,dim);
	}
	factor = (dim == 3) ? 1.0 / 27.0 : 1.0 / 9.0;
	multiplied_by_factor(nonlin_total,lin_total,diff_total,
	    nonlin_pr_total,lin_pr_total,diff_pr_total,
	    nonlin_vel_total,lin_vel_total,diff_vel_total,dim,factor);
}		/*end L2_regular_integral*/

LOCAL  void subdomain_integral(
	Locstate        diff_total,
	Locstate        nonlin_total,
	Locstate        lin_total,
	double           *diff_pr_total,
	double           *nonlin_pr_total,
	double           *lin_pr_total,
	double           *diff_vel_total,
	double           *nonlin_vel_total,
	double           *lin_vel_total,
	double           *L,
	double		*h,
	Wave            *wave,
	Front           *fr,
	OUTPUT_DATA     *data)
{
	Lp_Diff_data    *lp_data = Lp_data(data);
	COMPONENT       comp;
	INTERFACE       *intfc = fr->interf;
	int		stype = lp_data->stype;
	int		iz, iy, ix;
	int		dim = fr->rect_grid->dim;
	int             p = lp_data->p;
	int             ratio = lp_data->ratio;
	int		zmax;
	double		coords[3];
	double		factor;
	double		nonlin_pr, lin_pr, diff_pr;
	double		nonlin_vel[3], lin_vel[3], diff_vel[3];
	static Locstate nonlin = NULL, lin = NULL, diff = NULL;

	if (nonlin == NULL)
	{
		alloc_state(fr->interf,&nonlin,fr->sizest);
		alloc_state(fr->interf,&lin,fr->sizest);
		alloc_state(fr->interf,&diff,fr->sizest);
	}

	zero_states(nonlin_total,lin_total,diff_total,nonlin_pr_total,
		lin_pr_total,diff_pr_total,nonlin_vel_total,lin_vel_total,
		diff_vel_total,fr);

	zmax = (dim ==3) ? ratio : 1;
	for (iz = 0; iz < zmax; iz++)
	{
	    coords[2] = L[2] + (iz + 0.5) * h[2];
	    for (iy = 0; iy < ratio; iy++)
	    {
		coords[1] = L[1] + (iy + 0.5) * h[1];
		for (ix = 0; ix < ratio; ix++)
		{
		    coords[0] = L[0] + (ix + 0.5) * h[0];
		    comp = component(coords,intfc);
		    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				 fr,wave,nonlin,NULL);
		    (*lp_data->_alternate_state)(lin,coords,fr,wave,stype);

		    find_diff_and_other_states(diff,&nonlin_pr,&lin_pr,&diff_pr,
			     nonlin_vel,lin_vel,diff_vel,nonlin,lin,dim);
		    if (p == 1)
		    {
			 absolute_value(nonlin,lin,diff,
					&nonlin_pr,&lin_pr,&diff_pr,
			                nonlin_vel,lin_vel,diff_vel,dim);
			 add_to_total(nonlin_total,lin_total,diff_total,
			     nonlin_pr_total,lin_pr_total,diff_pr_total,
			     nonlin_vel_total,lin_vel_total,diff_vel_total,
			     nonlin,lin,diff,nonlin_pr,lin_pr,diff_pr,
			     nonlin_vel,lin_vel,diff_vel,dim);
		     }
		     else
			 add_product_to_total(nonlin_total,lin_total,diff_total,
			     nonlin_pr_total,lin_pr_total,diff_pr_total,
			     nonlin_vel_total,lin_vel_total,diff_vel_total,
			     dim,nonlin,nonlin,lin,lin,diff,diff,nonlin_pr,nonlin_pr,
			     lin_pr,lin_pr,diff_pr,diff_pr,nonlin_vel,nonlin_vel,
			     lin_vel,lin_vel,diff_vel,diff_vel);
		}
	    }
	}

	factor = 1.0 / ratio / ratio;
	if (dim == 3) factor /= ratio;
	multiplied_by_factor(nonlin_total,lin_total,diff_total,
	    nonlin_pr_total,lin_pr_total,diff_pr_total,
	    nonlin_vel_total,lin_vel_total,diff_vel_total,dim,factor);
}		/*end subdomain_integral*/

LOCAL  void zero_states(
	Locstate        nonlin,
	Locstate        lin,
	Locstate        diff,
	double           *nonlin_pr,
	double           *lin_pr,
	double           *diff_pr,
	double           *nonlin_vel,
	double           *lin_vel,
	double           *diff_vel,
	Front		*fr)
{
	int		i, dim = fr->rect_grid->dim;

	clear_state(fr->interf,nonlin,fr->sizest);
	clear_state(fr->interf,lin,fr->sizest);
	clear_state(fr->interf,diff,fr->sizest);
	*nonlin_pr = *lin_pr = *diff_pr = 0.0;
	for (i = 0; i < dim; i++)
	{
		nonlin_vel[i] = 0.0;
		lin_vel[i] = 0.0;
		diff_vel[i] = 0.0;
	}
}		/*end zero_states*/

LOCAL  void add_product_to_total(
        Locstate        nonlin_total,
        Locstate        lin_total,
        Locstate        diff_total,
        double           *nonlin_pr_total,
        double           *lin_pr_total,
        double           *diff_pr_total,
        double           *nonlin_vel_total,
        double           *lin_vel_total,
        double           *diff_vel_total,
	int		dim,
        Locstate        nonlin1,
        Locstate        nonlin2,
        Locstate        lin1,
        Locstate        lin2,
	Locstate	diff1,
	Locstate	diff2,
	double		nonlin_pr1,
	double		nonlin_pr2,
	double           lin_pr1,
	double           lin_pr2,
	double           diff_pr1,
	double           diff_pr2,
	double           *nonlin_vel1,
	double           *nonlin_vel2,
	double           *lin_vel1,
	double           *lin_vel2,
	double           *diff_vel1,
	double           *diff_vel2)
{
	int		i;

	Dens(nonlin_total) += Dens(nonlin1) * Dens(nonlin2);
	Energy(nonlin_total) += Energy(nonlin1) * Energy(nonlin2);
	Dens(lin_total) += Dens(lin1) * Dens(lin2);
	Energy(lin_total) += Energy(lin1) * Energy(lin2);
	Dens(diff_total) += Dens(diff1) * Dens(diff2);
	Energy(diff_total) += Energy(diff1) * Energy(diff2);
	*diff_pr_total += diff_pr1 * diff_pr2;
	*lin_pr_total += lin_pr1 * lin_pr2;
	*nonlin_pr_total += nonlin_pr1 * nonlin_pr2;
	for (i = 0; i < dim; i++)
	{
		Mom(nonlin_total)[i] += Mom(nonlin1)[i] * Mom(nonlin2)[i];
		Mom(lin_total)[i] += Mom(lin1)[i] * Mom(lin2)[i];
		Mom(diff_total)[i] += Mom(diff1)[i] * Mom(diff2)[i];
		diff_vel_total[i] += diff_vel1[i] * diff_vel2[i];
		lin_vel_total[i] += lin_vel1[i] * lin_vel2[i];
		nonlin_vel_total[i] += nonlin_vel1[i] * nonlin_vel2[i];
	}
	reset_gamma(lin_total);
	reset_gamma(nonlin_total);
	reset_gamma(diff_total);
}		/*end add_product_to_total*/

LOCAL  void add_multiple_of_product_to_total(
        Locstate        nonlin_total,
        Locstate        lin_total,
        Locstate        diff_total,
        double           *nonlin_pr_total,
        double           *lin_pr_total,
        double           *diff_pr_total,
        double           *nonlin_vel_total,
        double           *lin_vel_total,
        double           *diff_vel_total,
        Locstate        nonlin1,
        Locstate        nonlin2,
        Locstate        lin1,
        Locstate        lin2,
	Locstate	diff1,
	Locstate	diff2,
	double		nonlin_pr1,
	double		nonlin_pr2,
	double           lin_pr1,
	double           lin_pr2,
	double           diff_pr1,
	double           diff_pr2,
	double           *nonlin_vel1,
	double           *nonlin_vel2,
	double           *lin_vel1,
	double           *lin_vel2,
	double           *diff_vel1,
	double           *diff_vel2,
	Front           *fr,
	double		factor)
{
	int		dim = fr->rect_grid->dim;
	double           nonlin_pr, lin_pr, diff_pr;
	double           nonlin_vel[3], lin_vel[3], diff_vel[3];
	static Locstate nonlin = NULL, lin = NULL, diff = NULL;

	if (nonlin == NULL)
	{
	    alloc_state(fr->interf,&nonlin,fr->sizest);
	    alloc_state(fr->interf,&lin,fr->sizest);
	    alloc_state(fr->interf,&diff,fr->sizest);
	}
	zero_states(nonlin,lin,diff,&nonlin_pr,&lin_pr,&diff_pr,nonlin_vel,
		lin_vel,diff_vel,fr);
	add_product_to_total(nonlin,lin,diff,&nonlin_pr,&lin_pr,&diff_pr,
		nonlin_vel,lin_vel,diff_vel,dim,nonlin1,nonlin2,lin1,lin2,
		diff1,diff2,nonlin_pr1,nonlin_pr2,lin_pr1,lin_pr2,
		diff_pr1,diff_pr2,nonlin_vel1,nonlin_vel2,lin_vel1,lin_vel2,
		diff_vel1,diff_vel2);
	multiplied_by_factor(nonlin,lin,diff,&nonlin_pr,&lin_pr,&diff_pr,
		nonlin_vel,lin_vel,diff_vel,dim,factor);
	add_to_total(nonlin_total,lin_total,diff_total,nonlin_pr_total,
		lin_pr_total,diff_pr_total,nonlin_vel_total,lin_vel_total,
		diff_vel_total,nonlin,lin,diff,nonlin_pr,lin_pr,diff_pr,
		nonlin_vel,lin_vel,diff_vel,dim);
}		/*end add_multiple_of_product_to_total*/

LOCAL  void state_sqrt(
        Locstate        nonlin_total,
        Locstate        lin_total,
        Locstate        diff_total,
        double           *nonlin_pr_total,
        double           *lin_pr_total,
        double           *diff_pr_total,
        double           *nonlin_vel_total,
        double           *lin_vel_total,
        double           *diff_vel_total,
	int		dim)
{
	int		i;

	Dens(nonlin_total) = sqrt(Dens(nonlin_total));
	Energy(nonlin_total) = sqrt(Energy(nonlin_total));
	Dens(lin_total) = sqrt(Dens(lin_total));
	Energy(lin_total) = sqrt(Energy(lin_total));
	Dens(diff_total) = sqrt(Dens(diff_total));
	Energy(diff_total) = sqrt(Energy(diff_total));
	*diff_pr_total = sqrt(*diff_pr_total);
	*lin_pr_total = sqrt(*lin_pr_total);
	*nonlin_pr_total = sqrt(*nonlin_pr_total);

	for (i = 0; i < dim; i++)
	{
		Mom(nonlin_total)[i] = sqrt(Mom(nonlin_total)[i]);
		Mom(lin_total)[i] = sqrt(Mom(lin_total)[i]);
		Mom(diff_total)[i] = sqrt(Mom(diff_total)[i]);
		diff_vel_total[i] = sqrt(diff_vel_total[i]);
		lin_vel_total[i] = sqrt(lin_vel_total[i]);
		nonlin_vel_total[i] = sqrt(nonlin_vel_total[i]);
	}
	reset_gamma(lin_total);
	reset_gamma(nonlin_total);
	reset_gamma(diff_total);
}		/*end state_sqrt*/

LOCAL  void absolute_value(
	Locstate        nonlin,
	Locstate        lin,
	Locstate        diff,
	double           *nonlin_pr,
	double           *lin_pr,
	double           *diff_pr,
	double           *nonlin_vel,
	double           *lin_vel,
	double           *diff_vel,
	int             dim)
{
	int 		i;

	Dens(diff) = fabs(Dens(diff));
	Energy(diff) = fabs(Energy(diff));
	*nonlin_pr = fabs(*nonlin_pr);
	*lin_pr = fabs(*lin_pr);
	*diff_pr = fabs(*diff_pr);
	for (i = 0; i < dim; i++)
	{
		Mom(nonlin)[i] = fabs(Mom(nonlin)[i]);
		Mom(lin)[i] = fabs(Mom(lin)[i]);
		Mom(diff)[i] = fabs(Mom(diff)[i]);
		nonlin_vel[i] = fabs(nonlin_vel[i]);
		lin_vel[i] = fabs(lin_vel[i]);
		diff_vel[i] = fabs(diff_vel[i]);
	}
	reset_gamma(diff);
}		/*end absolute_value*/

LOCAL  void multiplied_by_factor(
	Locstate        nonlin_total,
	Locstate        lin_total,
	Locstate        diff_total,
	double           *nonlin_pr_total,
	double           *lin_pr_total,
	double           *diff_pr_total,
	double           *nonlin_vel_total,
	double           *lin_vel_total,
	double           *diff_vel_total,
	int             dim,
	double		factor)
{
	int i;

	Dens(nonlin_total) *= factor;
	Energy(nonlin_total) *= factor;
	Dens(lin_total) *= factor;
	Energy(lin_total) *= factor;
	Dens(diff_total) *= factor;
	Energy(diff_total) *= factor;
	*diff_pr_total *= factor;
	*lin_pr_total *= factor;
	*nonlin_pr_total *= factor;
	for (i = 0; i < dim; i++)
	{
		Mom(nonlin_total)[i] *= factor;
		Mom(lin_total)[i] *= factor;
		Mom(diff_total)[i] *= factor;
		diff_vel_total[i] *= factor;
		lin_vel_total[i] *= factor;
		nonlin_vel_total[i] *= factor;
	}
	reset_gamma(nonlin_total);
	reset_gamma(lin_total);
	reset_gamma(diff_total);
}		/*end multiplied_by_factor*/

LOCAL  void add_to_total(
	Locstate        nonlin_total,
	Locstate        lin_total,
	Locstate        diff_total,
	double           *nonlin_pr_total,
	double           *lin_pr_total,
	double           *diff_pr_total,
	double           *nonlin_vel_total,
	double           *lin_vel_total,
	double           *diff_vel_total,
	Locstate        nonlin,
	Locstate        lin,
	Locstate        diff,
	double           nonlin_pr,
	double           lin_pr,
	double           diff_pr,
	double           *nonlin_vel,
	double           *lin_vel,
	double           *diff_vel,
	int             dim)
{
	int i;

	Dens(nonlin_total) += Dens(nonlin);
        Energy(nonlin_total) += Energy(nonlin);
	Dens(lin_total) += Dens(lin);
	Energy(lin_total) += Energy(lin);
	Dens(diff_total) += Dens(diff);
	Energy(diff_total) += Energy(diff);
	*diff_pr_total += diff_pr;
	*lin_pr_total += lin_pr;
	*nonlin_pr_total += nonlin_pr;
	for (i = 0; i < dim; i++)
	{
		Mom(nonlin_total)[i] += Mom(nonlin)[i];
		Mom(lin_total)[i] += Mom(lin)[i];
		Mom(diff_total)[i] += Mom(diff)[i];
		diff_vel_total[i] += diff_vel[i];
		lin_vel_total[i] += lin_vel[i];
		nonlin_vel_total[i] += nonlin_vel[i];
	}
	reset_gamma(nonlin_total);
	reset_gamma(lin_total);
	reset_gamma(diff_total);

}		/*end add_to_total*/

LOCAL  void find_diff_and_other_states(
	Locstate        diff,
	double		*nonlin_pr,
	double		*lin_pr,
	double		*diff_pr,
	double		*nonlin_vel,
	double		*lin_vel,
	double		*diff_vel,
	Locstate	nonlin,
	Locstate	lin,
	int		dim)
{
	int		i;

	Dens(diff) = Dens(nonlin) - Dens(lin);
	Energy(diff) = Energy(nonlin) - Energy(lin);
	*nonlin_pr = pressure(nonlin);
	*lin_pr = pressure(lin);
	*diff_pr = *nonlin_pr - *lin_pr;
	for (i = 0; i < dim; i++)
	{
		Mom(diff)[i] = Mom(nonlin)[i] - Mom(lin)[i];
		lin_vel[i] = vel(i,lin);
		nonlin_vel[i] = vel(i,nonlin);
		diff_vel[i] = nonlin_vel[i] - lin_vel[i];
	}
	reset_gamma(diff);
}		/*end find_diff_and_other_states*/

LOCAL  void Lp_printout(
	OUTPUT_DATA     *data,
	Grid            *grid,
	Locstate	diff_total,
	Locstate	lin_total,
	Locstate	nonlin_total,
	double		diff_pr_total,
	double		lin_pr_total,
	double		nonlin_pr_total,
	double		*diff_vel_total,
	double		*lin_vel_total,
	double		*nonlin_vel_total,
	int		dim)
{
	FILE       *file = Output_file(data);
	int	   i;
	int        p = Lp_data(data)->p;
	static const char *MOM[] = {"X-MOMENTUM", "Y-MOMENTUM", "Z-MOMENTUM"};
	static const char *VEL[] = {"X-VELOCITY", "Y-VELOCITY", "Z-VELOCITY"};
	static boolean	first = YES;

	if (is_io_node(pp_mynode()))
	{
	    if (first == YES)
	    {
		first = NO;
	        if (file == NULL)
	        {
	            Output_file(data) = file = fopen(Output_filename(data),"w");
		    if (debugging("nobuf"))
			setbuf(file,NULL);
		}
		print_machine_parameters(file);

	        (void) foutput(file);
	        (void) fprintf(file,"%-21s","TIME");
		if (p == 1)
		{
	            (void) fprintf(file,"%-21s%-21s%-21s",
		       "L1_DIFF_DENSITY","L1_DIFF_ENERGY","L1_DIFF_PRESSURE");
		    for (i = 0; i < dim; i++)
		    {
		    	(void) fprintf(file,"L1_DIFF_%-13s",MOM[i]);
		    	(void) fprintf(file,"L1_DIFF_%-13s",VEL[i]);
		    }
	            (void) fprintf(file,"%-21s%-21s%-21s",
			"L1_LIN_DENSITY","L1_LIN_ENERGY","L1_LIN_PRESSURE");
		    for (i = 0; i < dim; i++)
		    {
		    	(void) fprintf(file,"L1_LIN_%-14s",MOM[i]);
		    	(void) fprintf(file,"L1_LIN_%-14s",VEL[i]);
		    }
	            (void) fprintf(file,"%-21s%-21s%-21s",
				   "L1_NONLIN_DENSITY","L1_NONLIN_ENERGY",
				   "L1_NONLIN_PRESSURE");
		    for (i = 0; i < dim; i++)
		    {
		    	(void) fprintf(file,"L1_NONLIN_%-11s",MOM[i]);
		    	(void) fprintf(file,"L1_NONLIN_%-11s",VEL[i]);
		    }
		}
		else
		{
	            (void) fprintf(file,"%-21s%-21s%-21s",
		        "L2_DIFF_DENSITY","L2_DIFF_ENERGY","L2_DIFF_PRESSURE");
		    for (i = 0; i < dim; i++)
		    {
		    	(void) fprintf(file,"L2_DIFF_%-13s",MOM[i]);
		    	(void) fprintf(file,"L2_DIFF_%-13s",VEL[i]);
		    }
	            (void) fprintf(file,"%-21s%-21s%-21s",
			"L2_LIN_DENSITY","L2_LIN_ENERGY","L2_LIN_PRESSURE");
		    for (i = 0; i < dim; i++)
		    {
		    	(void) fprintf(file,"L2_LIN_%-14s",MOM[i]);
		    	(void) fprintf(file,"L2_LIN_%-14s",VEL[i]);
		    }
	            (void) fprintf(file,"%-21s%-21s%-21s",
			    "L2_NONLIN_DENSITY","L2_NONLIN_ENERGY",
			    "L2_NONLIN_PRESSURE");
		    for (i = 0; i < dim; i++)
		    {
		    	(void) fprintf(file,"L2_NONLIN_%-11s",MOM[i]);
		    	(void) fprintf(file,"L2_NONLIN_%-11s",VEL[i]);
		    }
		 }
		 (void) fprintf(file,"\n");

	    } /* if (first == YES) */

	    (void) fprintf(file,"%-21g",grid->time);
	    (void) fprintf(file,"%-21g%-21g%-21g",
		    Dens(diff_total),Energy(diff_total),diff_pr_total);
	    for (i = 0; i < dim; i++)
	    {
		    (void) fprintf(file,"%-21g",Mom(diff_total)[i]);
		    (void) fprintf(file,"%-21g",diff_vel_total[i]);
	    }
	    (void) fprintf(file,"%-21g%-21g%-21g",
		    Dens(lin_total),Energy(lin_total),lin_pr_total);
	    for (i = 0; i < dim; i++)
	    {
		    (void) fprintf(file,"%-21g",Mom(lin_total)[i]);
		    (void) fprintf(file,"%-21g",lin_vel_total[i]);
	    }
	    (void) fprintf(file,"%-21g%-21g%-21g",
		    Dens(nonlin_total),Energy(nonlin_total),nonlin_pr_total);
	    for (i = 0; i < dim; i++)
	    {
	        (void) fprintf(file,"%-21g",Mom(nonlin_total)[i]);
	        (void) fprintf(file,"%-21g",nonlin_vel_total[i]);
	    }
	    (void) fprintf(file,"\n");
	    (void) fflush(file);
	} /* if (is_io_node(pp_mynode())) */
}		/*end Lp_printout*/

LOCAL  boolean is_block_on_front(
	int             *icoords,
	Wave            *wave,
	Front           *fr,
	double           *lin_intfc,
	int             indx)
{
	TRI_GRID        *tri_grid = wave_tri_soln(wave)->tri_grid;
	BLK_EL0         *blk_el0;
	int             dim = fr->rect_grid->dim;
	TG_PT           **p;
	double           zl, zu;

	blk_el0 = &Blk_el0(icoords,tri_grid);
	if (!blk_el0_is_bilinear(blk_el0))
	    return YES;

	p = blk_el0_bilinear_el(blk_el0)->p;
	switch (dim)
	{
#if defined(ONED)
	case 1:
	    break;
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	    zl = Coords(p[0])[1];
	    zu = Coords(p[(1<<dim)-1])[1];
	    if (((zl < lin_intfc[indx])   && (lin_intfc[indx]   < zu)) ||
		((zl < lin_intfc[indx+1]) && (lin_intfc[indx+1] < zu)) ||
		((zl < lin_intfc[indx+2]) && (lin_intfc[indx+2] < zu)))
		    return YES;
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    {
		int i;
	        int *gmax = fr->rect_grid->gmax;
	        zl = Coords(p[0])[2];
	        zu = Coords(p[(1<<dim)-1])[2];
	        i = indx;
	        if ((zl < lin_intfc[i]) && (lin_intfc[i] < zu))
		    return YES;
	        i = indx + 1;
	        if ((zl < lin_intfc[i]) && (lin_intfc[i] < zu))
		    return YES;
	        i += gmax[0];
	        if ((zl < lin_intfc[i]) && (lin_intfc[i] < zu))
		    return YES;
	        i += gmax[0] + 1;
	        if ((zl < lin_intfc[i]) && (lin_intfc[i] < zu))
		    return YES;
	        i += 1;
	        if ((zl < lin_intfc[i]) && (lin_intfc[i] < zu))
		    return YES;
	    }
	    break;
#endif /* defined(THREED) */
	}

	return NO;
}		/*end is_block_on_front*/
#endif /* defined(TWOD) || defined(THREED) */
