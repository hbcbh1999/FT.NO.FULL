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
*				girpregion.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <ginit/ginit.h>

/* LOCAL Function prototypes */
LOCAL	COMP_TYPE *clone_elliptical_comp_type(ELLIPSOID*,COMP_TYPE*,
                                                    Front*,INIT_PHYSICS*);
LOCAL	ELLIPSOID *clone_ellipsoid(LAYER*,ELLIPSOID*);
LOCAL	LAYER	  *add_new_layer(LAYER_SYS*,COMPONENT);
LOCAL	void	  set_riemann_problem_ellipsoid(Front*,LAYER*,ELLIPSOID*,
                                                RIEMANN_SOLVER_WAVE_TYPE,
						RIEMANN_SOLVER_WAVE_TYPE,double,
						double,Locstate,Locstate,
						Locstate,Locstate,
						INIT_PHYSICS*);

/*ARGSUSED*/
EXPORT	void	set_up_riemann_problem_region(
	int		layer_label,
	int		region_label,
	int		surf_label,
	int		ellip_label,
	double		*coords,
	double		*nor,
	SIDE            ahead_side,
	Locstate	ahead_st,
	Locstate	st,
	LAYER_SYS	*layer_sys,
	INIT_PHYSICS	*ip,
	INIT_DATA	*init)
{
	COMP_TYPE		 *ct;
	_RAREFACTION_WAVE_1D	 *rw1d;
	LAYER                    *lyr, *llyr, *ulyr;
	Front	                 *front = layer_sys->front;
	INTERFACE                *intfc = front->interf;
	ELLIPSOID                *ellip;
	LAYER_SURF               *lsurf;
	Locstate		 sl, sr;
	Locstate                 left, right;
	Locstate                 sml, smr;
	double                    vl, vr;
	double                    pjump;
	double                    pml, pmr, uml, umr, ml, mr;
	double                    cl, cr, cml, cmr, Wl, Wr;
	double                    W, V;
	double                    dt, dh;
	RIEMANN_SOLVER_WAVE_TYPE l_wave, r_wave;
	int                      i, dim = front->rect_grid->dim;
	int                      w_type;
	size_t                   sizest = front->sizest;

	debug_print("layer","Entered set_up_riemann_problem_region()\n");

	alloc_state(intfc,&left,max(sizeof(VGas),sizest));
	alloc_state(intfc,&right,max(sizeof(VGas),sizest));
	alloc_state(intfc,&sml,sizest);
	alloc_state(intfc,&smr,sizest);

	if (ahead_side == POSITIVE_SIDE)
	{
	    sl = st;
	    sr = ahead_st;
	}
	else
	{
	    sl = ahead_st;
	    sr = st;
	}
	set_state(right,TGAS_STATE,sr);
	set_state(left,TGAS_STATE,sl);
	lyr = layer_sys->layer[layer_label];
	if (ellip_label > 0)
	{
	    if (ellip_label <= lyr->num_ellips)
	        ellip = lyr->ellip[ellip_label];
	    else
	    {
	        screen("ERROR in set_up_riemann_problem_region(), "
		       "invalid ellip_label %d > num_ellips %d\n",
		       ellip_label,lyr->num_ellips);
		clean_up(ERROR);
	    }
	    lsurf = NULL;
	}
	else
	{
	    ellip = NULL;
	    if (surf_label == layer_label)
	    {
	        lsurf = lyr->upper_surf;
		llyr = lyr;
		ulyr = layer_sys->layer[layer_label+1];
	    }
	    else if (surf_label == layer_label-1)
	    {
	        lsurf = lyr->lower_surf;
		ulyr = lyr;
		llyr = layer_sys->layer[layer_label-1];
	    }
	    else
	    {
	        screen("ERROR in set_up_riemann_problem_region(), "
		       "invalid surf_label %d layer_label %d\n",
		       surf_label,layer_label);
		clean_up(ERROR);
	    }
	}

	if (ellip != NULL)
	{
	    vr = Vel(right)[0];
	    vl = Vel(left)[0];
	    pjump = -ellip->surf_tension/
	             distance_between_positions(coords,ellip->cen,dim);
	}
	else
	{
	    vr = scalar_product(Vel(right),nor,dim);
	    vl = scalar_product(Vel(left),nor,dim);
	    pjump = 0.0;
	}
	zero_state_velocity(right,dim);
	Vel(right)[0] = vr;
	zero_state_velocity(left,dim);
	Vel(right)[0] = vl;
	set_state_for_find_mid_state(right,right);
	set_state_for_find_mid_state(left,left);
	if (find_mid_state(left,right,pjump,&pml,&pmr,&uml,&umr,&ml,&mr,
	                       &l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	{
	    screen("ERROR in set_up_riemann_problem_region(), "
	           "find_mid_state() did not converge\n");
	    verbose_print_state("left",left);
	    verbose_print_state("right",right);
	    (void) printf("pjump = %g\n"
	                  "pml = %g, pmr = %g\n"
	                  "uml = %g, umr = %g\n"
	                  "ml = %g, mr = %g\n",
	                  pjump,pml,pmr,uml,umr,ml,mr);
	    (void) printf("l_wave = %s, r_wave = %s\n",
	                  rsoln_wave_name(l_wave),rsoln_wave_name(r_wave));
	    clean_up(ERROR);
	}

	w_type = (l_wave == RAREFACTION) ?
	    BACKWARD_SOUND_WAVE_TE : BACKWARD_SHOCK_WAVE;
	state_behind_sound_wave(left,sml,&cml,&Wl,0.0,ml,uml,pml,TGAS_STATE,
	                        w_type,l_wave,LEFT_FAMILY);
	w_type = (r_wave == RAREFACTION) ?
	    FORWARD_SOUND_WAVE_TE : FORWARD_SHOCK_WAVE;
	state_behind_sound_wave(right,smr,&cmr,&Wr,0.0,mr,umr,pmr,TGAS_STATE,
	                        w_type,r_wave,RIGHT_FAMILY);

	cl = sound_speed(left);
	cr = sound_speed(right);
	W = max(fabs(Wl),fabs(Wr));
	V = fabs(vl) + cl;
	W = max(W,V);
	V = fabs(vr) + cr;
	W = max(W,V);
	V = fabs(uml) + cml;
	W = max(W,V);
	V = fabs(umr) + cmr;
	W = max(W,V);
	for (dh = HUGE_VAL, i = 0; i < dim; ++i)
	    dh = min(dh,front->rect_grid->h[i]);
	dt = 0.1*dh/W;/*TOLERANCE*/
	layer_sys->dt = min(layer_sys->dt,dt);

	if (debugging("layer"))
	{
	    (void) printf("States from Riemann solution\n");
	    verbose_print_state("left state",left);
	    verbose_print_state("left mid state",sml);
	    verbose_print_state("right mid state",smr);
	    verbose_print_state("right state",right);
	    (void) printf("l_wave = %s, r_wave = %s\n",
	                  rsoln_wave_name(l_wave),rsoln_wave_name(r_wave));
	    (void) printf("Wave speeds\n");
	    if (l_wave == RAREFACTION)
	    {
	        (void) printf("Left rarefaction leading edge speed = %g\n",
		              vl-cl);
	        (void) printf("Left rarefaction trailing edge speed = %g\n",
		              uml-cml);
	    }
	    else if (l_wave == SHOCK)
	        (void) printf("Left shock speed = %g\n",Wl);
	    (void) printf("Contact speed = %g (uml = %g, umr = %g)\n",
	                  0.5*(uml+umr),uml,umr);
	    if (r_wave == RAREFACTION)
	    {
	        (void) printf("Right rarefaction trailing edge speed = %g\n",
		              umr+cmr);
	        (void) printf("Right rarefaction leading edge speed = %g\n",
		              vr+cr);
	    }
	    else if (r_wave == SHOCK)
	        (void) printf("Right shock speed = %g\n",Wr);

	}

	if (ellip == NULL)
	{
	    LAYER      *rlyr, *mlyr;
	    LAYER_SURF *rlyr_le, *rlyr_te, *shock;
	    double      *nor = lsurf->nor;
	    double      vml, vmr;

	    vml = Vel(sml)[0];
	    vmr = Vel(smr)[0];
	    for (i = 0; i < dim; ++i)
	    {
		Vel(sml)[i] = vel(i,sl) + (vml-vl)*nor[i];
		Vel(smr)[i] = vel(i,sr) + (vmr-vr)*nor[i];
	    }
	    if (l_wave == RAREFACTION)
	    {
	        rlyr = add_new_layer(layer_sys,new_component(NEW_COMP));
	        rlyr->lower_surf = alloc_layer_surf();
	        rlyr->upper_surf = alloc_layer_surf();
	        *rlyr->upper_surf = *lsurf;
	        *rlyr->lower_surf = *lsurf;
		rlyr->lower_surf->reset_position = YES;
		rlyr->upper_surf->reset_position = YES;
	        rlyr->lower_surf->surf_ten = 0.0;
	        rlyr->upper_surf->surf_ten = 0.0;

		ct = comp_type(rlyr->comp);
	        set_rarefaction_wave_1d_comp_type(ct,front);
	        rw1d = Rarefaction_wave_1d(ct);
		rw1d->l_or_r = LEFT_FAMILY;
		rw1d->zbar = lsurf->pbar[dim-1];
		rw1d->tbar = -HUGE_VAL; /*To be set later */
		rw1d->el_lead = rw1d->el_trail = NULL;
		set_state(rw1d->stl,TGAS_STATE,sl);
		set_state(rw1d->stt,TGAS_STATE,sml);
		rw1d->spl = vl-cl;
		rw1d->spt = vml-cml;

	        mlyr = add_new_layer(layer_sys,new_component(NEW_COMP));

	        if (nor[dim-1] < 0.0)
	        {
	    	    rlyr_le = rlyr->upper_surf;
	    	    rlyr_te = rlyr->lower_surf;
	    	    mlyr->upper_surf = rlyr_te;
	    	    mlyr->lower_surf = lsurf;
		    rw1d->zt = rw1d->zmin = -HUGE_VAL;/*To be set later*/
	            rw1d->stmin = rw1d->stt;
		    rw1d->zl = rw1d->zmax =  HUGE_VAL;/*To be set later*/
	            rw1d->stmax = rw1d->stl;
		    ulyr->lower_surf = rlyr_le;
	        }
	        else
	        {
	    	    rlyr_le = rlyr->lower_surf;
	    	    rlyr_te = rlyr->upper_surf;
	    	    mlyr->lower_surf = rlyr_te;
	    	    mlyr->upper_surf = lsurf;
		    rw1d->zl = rw1d->zmin = -HUGE_VAL;/*To be set later*/
	            rw1d->stmin = rw1d->stl;
		    rw1d->zt = rw1d->zmax =  HUGE_VAL;/*To be set later*/
	            rw1d->stmax = rw1d->stt;
		    llyr->upper_surf = rlyr_le;
	        }
		rw1d->lead = rlyr_le;
		rw1d->trail = rlyr_te;
	        rlyr_le->l_comp = lsurf->l_comp;
	        rlyr_le->r_comp = rlyr_te->l_comp = rlyr->comp;
	        rlyr_le->wv_type = BACKWARD_SOUND_WAVE_LE;
	        rlyr_te->wv_type = BACKWARD_SOUND_WAVE_TE;
	        rlyr_te->r_comp = lsurf->l_comp = mlyr->comp;
	        for (i = 0; i < dim; ++i)
		{
		    rlyr_le->velocity[i] = rw1d->spl*nor[i];
		    rlyr_te->velocity[i] = rw1d->spt*nor[i];
		}
	    }
	    else
	    {
	        mlyr = add_new_layer(layer_sys,new_component(NEW_COMP));
	        shock = alloc_layer_surf();
	        *shock = *lsurf;
	        if (nor[dim-1] < 0.0)
	        {
	    	    mlyr->upper_surf = shock;
	    	    mlyr->lower_surf = lsurf;
		    ulyr->lower_surf = shock;
	        }
	        else
	        {
	    	    mlyr->lower_surf = shock;
	    	    mlyr->upper_surf = lsurf;
		    llyr->upper_surf = shock;
	        }
	        shock->l_comp = lsurf->l_comp;
	        shock->wv_type = BACKWARD_SHOCK_WAVE;
	        shock->r_comp = lsurf->l_comp = mlyr->comp;
		shock->reset_position = YES;
	        for (i = 0; i < dim; ++i)
		    shock->velocity[i] = Wl*nor[i];
	    }
	    ct = comp_type(mlyr->comp);
	    set_ambient_comp_type(ct,front);
	    set_state(Ambient(ct),GAS_STATE,sml);
	    if (r_wave == RAREFACTION)
	    {
	        rlyr = add_new_layer(layer_sys,new_component(NEW_COMP));
	        rlyr->lower_surf = alloc_layer_surf();
	        rlyr->upper_surf = alloc_layer_surf();
	        *rlyr->upper_surf = *lsurf;
	        *rlyr->lower_surf = *lsurf;
		rlyr->lower_surf->reset_position = YES;
		rlyr->upper_surf->reset_position = YES;
	        rlyr->lower_surf->surf_ten = 0.0;
	        rlyr->upper_surf->surf_ten = 0.0;

		ct = comp_type(rlyr->comp);
	        set_rarefaction_wave_1d_comp_type(ct,front);
	        rw1d = Rarefaction_wave_1d(ct);
		rw1d->l_or_r = RIGHT_FAMILY;
		rw1d->zbar = lsurf->pbar[dim-1];
		rw1d->tbar = -HUGE_VAL; /*To be set later */
		rw1d->el_lead = rw1d->el_trail = NULL;
		set_state(rw1d->stl,TGAS_STATE,sr);
		set_state(rw1d->stt,TGAS_STATE,smr);
		rw1d->spl = vr+cr;
		rw1d->spt = vmr+cmr;
		rw1d->lead = rlyr_le;
		rw1d->trail = rlyr_te;

	        mlyr = add_new_layer(layer_sys,new_component(NEW_COMP));

	        if (nor[dim-1] < 0.0)
	        {
	    	    rlyr_le = rlyr->lower_surf;
	    	    rlyr_te = rlyr->upper_surf;
	    	    mlyr->lower_surf = rlyr_te;
	    	    mlyr->upper_surf = lsurf;
		    rw1d->zl = rw1d->zmin = -HUGE_VAL;/*To be set later*/
	            rw1d->stmin = rw1d->stl;
		    rw1d->zt = rw1d->zmax =  HUGE_VAL;/*To be set later*/
	            rw1d->stmax = rw1d->stt;
		    llyr->upper_surf = rlyr_le;
	        }
	        else
	        {
	    	    rlyr_le = rlyr->upper_surf;
	    	    rlyr_te = rlyr->lower_surf;
	    	    mlyr->upper_surf = rlyr_te;
	    	    mlyr->lower_surf = lsurf;
		    rw1d->zt = rw1d->zmin = -HUGE_VAL;/*To be set later*/
	            rw1d->stmin = rw1d->stt;
		    rw1d->zl = rw1d->zmax =  HUGE_VAL;/*To be set later*/
	            rw1d->stmax = rw1d->stl;
		    ulyr->lower_surf = rlyr_le;
	        }
	        rlyr_le->r_comp = lsurf->r_comp;
	        rlyr_le->l_comp = rlyr_te->r_comp = rlyr->comp;
	        rlyr_le->wv_type = FORWARD_SOUND_WAVE_LE;
	        rlyr_te->wv_type = FORWARD_SOUND_WAVE_TE;
	        rlyr_te->l_comp = lsurf->r_comp = mlyr->comp;
	        for (i = 0; i < dim; ++i)
	        {
	    	    rlyr_le->velocity[i] = rw1d->spl*nor[i];
	    	    rlyr_te->velocity[i] = rw1d->spt*nor[i];
	        }
	    }
	    else
	    {
	        mlyr = add_new_layer(layer_sys,new_component(NEW_COMP));
	        shock = alloc_layer_surf();
	        *shock = *lsurf;
	        if (nor[dim-1] < 0.0)
	        {
	    	    mlyr->lower_surf = shock;
	    	    mlyr->upper_surf = lsurf;
		    llyr->upper_surf = shock;
	        }
	        else
	        {
	    	    mlyr->upper_surf = shock;
	    	    mlyr->lower_surf = lsurf;
		    ulyr->lower_surf = shock;
	        }
	        shock->r_comp = lsurf->r_comp;
	        shock->wv_type = FORWARD_SHOCK_WAVE;
	        shock->l_comp = lsurf->r_comp = mlyr->comp;
		shock->reset_position = YES;
	        for (i = 0; i < dim; ++i)
	    	    shock->velocity[i] = Wr*nor[i];
	    }
	    ct = comp_type(mlyr->comp);
	    set_ambient_comp_type(ct,front);
	    set_state(Ambient(ct),GAS_STATE,smr);

	    lsurf->wv_type = CONTACT;
	    lsurf->reset_position = YES;
	    for (i = 0; i < dim; ++i)
	        lsurf->velocity[i] = 0.5*(vml+vmr)*nor[i];
	}
	else
	{
	    set_riemann_problem_ellipsoid(front,lyr,ellip,l_wave,r_wave,
	                                  Wl,Wr,sl,sml,smr,sr,ip);
	}


	free_these(4,left,right,sml,smr);
	debug_print("layer","Left set_up_riemann_problem_region()\n");
}		/*end set_up_riemann_problem_region*/

LOCAL	void	set_riemann_problem_ellipsoid(
	Front		         *front,
	LAYER                    *lyr,
	ELLIPSOID                *ellip,
	RIEMANN_SOLVER_WAVE_TYPE l_wave,
	RIEMANN_SOLVER_WAVE_TYPE r_wave,
	double                    Wl,
	double                    Wr,
	Locstate                 sl,
	Locstate                 sml,
	Locstate                 smr,
	Locstate                 sr,
	INIT_PHYSICS		 *ip)
{
	_ELLIPTICAL          *rel, *mel, *el;
	_RAREFACTION_WAVE_1D *rw1d;
	COMP_TYPE            *rct, *mct, *ct;
	ELLIPSOID            *le_ellip, *te_ellip, *shock;
	INTERFACE            *intfc = front->interf;
	double                rmax;
	int                  i, dim = intfc->dim;

	ct = comp_type(ellip->compin);
	rmax = max_radii(ellip);
	if (l_wave == RAREFACTION)
	{
	    le_ellip = clone_ellipsoid(lyr,ellip);
	    le_ellip->wv_type = BACKWARD_SOUND_WAVE_LE;
	    te_ellip = clone_ellipsoid(lyr,ellip);
	    te_ellip->wv_type = BACKWARD_SOUND_WAVE_TE;
	    if (ellip->nor_orient == POSITIVE_ORIENTATION)
	    {
	        mct = clone_elliptical_comp_type(ellip,ct,front,ip);
	        rct = clone_elliptical_comp_type(te_ellip,ct,front,ip);
		if (ct->type == ELLIPTICAL)
		{
		    el = Elliptical(ct);
		    el->ellipsoid = le_ellip;
		}
	        le_ellip->compin = ellip->compin;
	        le_ellip->compout = rct->comp;
	        te_ellip->compin = le_ellip->compout;
	        te_ellip->compout = mct->comp;
	        ellip->compin = te_ellip->compout;
	        inner_ellipsoid(le_ellip) = inner_ellipsoid(ellip);
	        outer_ellipsoid(le_ellip) = te_ellip;
	        inner_ellipsoid(te_ellip) = le_ellip;
	        outer_ellipsoid(te_ellip) = ellip;
	        inner_ellipsoid(ellip) = te_ellip;
	    }
	    else if (ellip->nor_orient == NEGATIVE_ORIENTATION)
	    {
	        rct = clone_elliptical_comp_type(le_ellip,ct,front,ip);
	        mct = clone_elliptical_comp_type(te_ellip,ct,front,ip);
	        le_ellip->compout = ellip->compout;
	        le_ellip->compin = rct->comp;
	        te_ellip->compout = le_ellip->compin;
	        te_ellip->compin = mct->comp;
	        ellip->compout = te_ellip->compin;
	        outer_ellipsoid(le_ellip) = outer_ellipsoid(ellip);
	        inner_ellipsoid(le_ellip) = te_ellip;
	        outer_ellipsoid(te_ellip) = le_ellip;
	        inner_ellipsoid(te_ellip) = ellip;
	        outer_ellipsoid(ellip) = te_ellip;
	    }
	    else
	    {
	        screen("ERROR in set_up_riemann_problem_region(), "
	               "ellip->nor_orient not set\n");
	        clean_up(ERROR);
	    }
	    rel = Elliptical(rct);
	    rel->rw1d = rw1d = allocate_RAREFACTION_WAVE_1D(front);
	    rw1d->l_or_r = LEFT_FAMILY;
	    rw1d->zbar = rmax;
	    rw1d->tbar = -HUGE_VAL; /*To be set later */
	    rw1d->el_lead = le_ellip;
	    rw1d->el_trail = te_ellip;
	    rw1d->lead = rw1d->trail = NULL;
	    set_state(rw1d->stl,TGAS_STATE,sl);
	    set_state(rw1d->stt,TGAS_STATE,sml);
	    le_ellip->reset_position = te_ellip->reset_position = YES;
	    rw1d->spl = vel(0,sl) - sound_speed(sl);
	    rw1d->spt = vel(0,sml) - sound_speed(sml);
	    for (i = 0; i < dim; ++i)
	    {
	        le_ellip->vr[i] = rw1d->spl*ellip->rad[i]/rmax;
	        te_ellip->vr[i] = rw1d->spt*ellip->rad[i]/rmax;
	    }
	    if (ellip->nor_orient == POSITIVE_ORIENTATION)
	    {
	        rw1d->zl = rw1d->zmin = -HUGE_VAL;
		rw1d->zt = rw1d->zmax = HUGE_VAL;
		rw1d->stmin = rw1d->stl;
		rw1d->stmax = rw1d->stt;
	    }
	    else
	    {
	        rw1d->zt = rw1d->zmin = -HUGE_VAL;
		rw1d->zl = rw1d->zmax = HUGE_VAL;
		rw1d->stmin = rw1d->stt;
		rw1d->stmax = rw1d->stl;
	    }
	}
	else
	{
	    shock = clone_ellipsoid(lyr,ellip);
	    shock->wv_type = BACKWARD_SHOCK_WAVE;
	    if (ellip->nor_orient == POSITIVE_ORIENTATION)
	    {
	        mct = clone_elliptical_comp_type(ellip,ct,front,ip);
		if (ct->type == ELLIPTICAL)
		{
		    el = Elliptical(ct);
		    el->ellipsoid = shock;
		}
	        shock->compin = ellip->compin;
	        shock->compout = mct->comp;
	        ellip->compin = shock->compout;
	        inner_ellipsoid(shock) = inner_ellipsoid(ellip);
	        outer_ellipsoid(shock) = ellip;
	        inner_ellipsoid(ellip) = shock;
	    }
	    else if (ellip->nor_orient == NEGATIVE_ORIENTATION)
	    {
	        mct = clone_elliptical_comp_type(shock,ct,front,ip);
	        shock->compout = ellip->compout;
	        shock->compin = mct->comp;
	        ellip->compout = shock->compin;
	        outer_ellipsoid(shock) = outer_ellipsoid(ellip);
	        inner_ellipsoid(shock) = ellip;
	        outer_ellipsoid(ellip) = shock;
	    }
	    else
	    {
	        screen("ERROR in set_up_riemann_problem_region(), "
	               "ellip->nor_orient not set\n");
	        clean_up(ERROR);
	    }
	    shock->reset_position = YES;
	    for (i = 0; i < dim; ++i)
	        shock->vr[i] = Wl*ellip->rad[i]/rmax;
	}
	mel = Elliptical(mct);
	set_state(mel->state,TGAS_STATE,sml);
	if (r_wave == RAREFACTION)
	{
	    le_ellip = clone_ellipsoid(lyr,ellip);
	    le_ellip->wv_type = FORWARD_SOUND_WAVE_LE;
	    te_ellip = clone_ellipsoid(lyr,ellip);
	    te_ellip->wv_type = FORWARD_SOUND_WAVE_TE;
	    if (ellip->nor_orient == NEGATIVE_ORIENTATION)
	    {
	        mct = clone_elliptical_comp_type(ellip,ct,front,ip);
	        rct = clone_elliptical_comp_type(te_ellip,ct,front,ip);
		if (ct->type == ELLIPTICAL)
		{
		    el = Elliptical(ct);
		    el->ellipsoid = le_ellip;
		}
	        le_ellip->compin = ellip->compin;
	        le_ellip->compout = rct->comp;
	        te_ellip->compin = le_ellip->compout;
	        te_ellip->compout = mct->comp;
	        ellip->compin = te_ellip->compout;
	        inner_ellipsoid(le_ellip) = inner_ellipsoid(ellip);
	        outer_ellipsoid(le_ellip) = te_ellip;
	        inner_ellipsoid(te_ellip) = le_ellip;
	        outer_ellipsoid(te_ellip) = ellip;
	        inner_ellipsoid(ellip) = te_ellip;
	    }
	    else if (ellip->nor_orient == POSITIVE_ORIENTATION)
	    {
	        rct = clone_elliptical_comp_type(le_ellip,ct,front,ip);
	        mct = clone_elliptical_comp_type(te_ellip,ct,front,ip);
	        le_ellip->compout = ellip->compout;
	        le_ellip->compin = rct->comp;
	        te_ellip->compout = le_ellip->compin;
	        te_ellip->compin = mct->comp;
	        ellip->compout = te_ellip->compin;
	        outer_ellipsoid(le_ellip) = outer_ellipsoid(ellip);
	        inner_ellipsoid(le_ellip) = te_ellip;
	        outer_ellipsoid(te_ellip) = le_ellip;
	        inner_ellipsoid(te_ellip) = ellip;
	        outer_ellipsoid(ellip) = te_ellip;
	    }
	    else
	    {
	        screen("ERROR in set_up_riemann_problem_region(), "
	               "ellip->nor_orient not set\n");
	        clean_up(ERROR);
	    }
	    rel = Elliptical(rct);
	    rel->rw1d = rw1d = allocate_RAREFACTION_WAVE_1D(front);
	    rw1d->l_or_r = RIGHT_FAMILY;
	    rw1d->zbar = rmax;
	    rw1d->tbar = -HUGE_VAL; /*To be set later */
	    rw1d->el_lead = le_ellip;
	    rw1d->el_trail = te_ellip;
	    rw1d->lead = rw1d->trail = NULL;
	    set_state(rw1d->stl,TGAS_STATE,sr);
	    set_state(rw1d->stt,TGAS_STATE,smr);
	    le_ellip->reset_position = te_ellip->reset_position = YES;
	    rw1d->spl = vel(0,sr) + sound_speed(sr);
	    rw1d->spt = vel(0,smr) + sound_speed(smr);
	    for (i = 0; i < dim; ++i)
	    {
	        le_ellip->vr[i] = rw1d->spl*ellip->rad[i]/rmax;
	        te_ellip->vr[i] = rw1d->spt*ellip->rad[i]/rmax;
	    }
	    if (ellip->nor_orient == NEGATIVE_ORIENTATION)
	    {
	        rw1d->zl = rw1d->zmin = -HUGE_VAL;
		rw1d->zt = rw1d->zmax = HUGE_VAL;
		rw1d->stmin = rw1d->stl;
		rw1d->stmax = rw1d->stt;
	    }
	    else
	    {
	        rw1d->zt = rw1d->zmin = -HUGE_VAL;
		rw1d->zl = rw1d->zmax = HUGE_VAL;
		rw1d->stmin = rw1d->stt;
		rw1d->stmax = rw1d->stl;
	    }
	}
	else
	{
	    shock = clone_ellipsoid(lyr,ellip);
	    shock->wv_type = FORWARD_SHOCK_WAVE;
	    if (ellip->nor_orient == NEGATIVE_ORIENTATION)
	    {
	        mct = clone_elliptical_comp_type(ellip,ct,front,ip);
		if (ct->type == ELLIPTICAL)
		{
		    el = Elliptical(ct);
		    el->ellipsoid = shock;
		}
	        shock->compin = ellip->compin;
	        shock->compout = mct->comp;
	        ellip->compin = shock->compout;
	        inner_ellipsoid(shock) = inner_ellipsoid(ellip);
	        outer_ellipsoid(shock) = ellip;
	        inner_ellipsoid(ellip) = shock;
	    }
	    else if (ellip->nor_orient == POSITIVE_ORIENTATION)
	    {
	        mct = clone_elliptical_comp_type(shock,ct,front,ip);
	        shock->compout = ellip->compout;
	        shock->compin = mct->comp;
	        ellip->compout = shock->compin;
	        outer_ellipsoid(shock) = outer_ellipsoid(ellip);
	        inner_ellipsoid(shock) = ellip;
	        outer_ellipsoid(ellip) = shock;
	    }
	    else
	    {
	        screen("ERROR in set_up_riemann_problem_region(), "
	               "ellip->nor_orient not set\n");
	        clean_up(ERROR);
	    }
	    shock->reset_position = YES;
	    for (i = 0; i < dim; ++i)
	        shock->vr[i] = Wr*ellip->rad[i]/rmax;
	}
	mel = Elliptical(mct);
	set_state(mel->state,TGAS_STATE,smr);

	ellip->wv_type = CONTACT;
	ellip->reset_position = YES;
	for (i = 0; i < dim; ++i)
	    ellip->vr[i] = 0.5*(vel(0,sml)+vel(0,smr))*ellip->rad[i]/rmax;
}		/*end set_riemann_problem_ellipsoid*/

LOCAL	COMP_TYPE *clone_elliptical_comp_type(
	ELLIPSOID    *ellip,
	COMP_TYPE    *ct,
	Front        *front,
	INIT_PHYSICS *ip)
{
	_ELLIPTICAL          *el, *nel;
	COMP_TYPE            *nct;
	INTERFACE            *intfc = front->interf;
	int                  i, dim;

	nct = comp_type(new_component(NEW_COMP));
	set_elliptical_comp_type(nct,ip);
	nel = Elliptical(nct);
	nel->ellipsoid = ellip;
	switch (ct->type)
	{
	case ELLIPTICAL:
	    dim = intfc->dim;
	    el = Elliptical(ct);
	    nel->rstate = copy_random_state_structure(el->rstate,intfc);
	    ft_assign(nel->state,el->state,front->sizest);
	    ft_assign(nel->wkstate[0],el->wkstate[0],front->sizest);
	    ft_assign(nel->wkstate[1],el->wkstate[1],front->sizest);
	    for (i = 0; i < dim; ++i)
	        nel->weight[i] = el->weight[i];
	    nel->r0 = el->r0;
	    nel->rw1d = el->rw1d;
	    nel->stratification_type = el->stratification_type;
	    break;
	case AMBIENT:
	    nel->rstate = NULL;
	    set_state(nel->state,TGAS_STATE,Ambient(ct));
	    break;
	default:
	    screen("ERROR in copy_elliptical_comp_type(), "
	           "comp_type %s not supported\n",comp_type_name(ct->type));
	    clean_up(ERROR);
	    break;
	}
	return nct;
}		/*end clone_elliptical_comp_type*/

LOCAL	ELLIPSOID	*clone_ellipsoid(
	LAYER     *lyr,
	ELLIPSOID *template_ellip)
{
	ELLIPSOID *ellip;

	ellip = allocate_ellipsoid(template_ellip,0);
	ellip->fpoly = template_ellip->fpoly;
	ellip->lpoly = template_ellip->lpoly;
	if (lyr != NULL)
	    lyr->ellip[++lyr->num_ellips] = ellip;

	return ellip;
}		/*end clone_ellipsoid*/

LOCAL	LAYER *add_new_layer(
	LAYER_SYS *layer_sys,
	COMPONENT comp)
{
	LAYER *lyr;
	scalar(&lyr,sizeof(LAYER));
	lyr->prev = layer_sys->layer[layer_sys->num_layers];
	layer_sys->layer[layer_sys->num_layers]->next = lyr;
	layer_sys->layer[++layer_sys->num_layers] = lyr;
	lyr->layer_label = layer_sys->num_layers;
	lyr->comp = comp;
	lyr->num_ellips = 0;
	lyr->ellip = NULL;
	lyr->lower_surf = lyr->upper_surf = NULL;
	return lyr;
}		/*end add_new_layer*/

