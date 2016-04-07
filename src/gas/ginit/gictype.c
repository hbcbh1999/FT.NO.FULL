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
*				gictype.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <ginit/ginit.h>

/* LOCAL Function prototypes */
LOCAL   boolean	perturbed_layer_surf(LAYER_SURF*);
LOCAL	boolean	perturbed_rarefaction_layer(_RAREFACTION_WAVE_1D*);
LOCAL	void	get_state_trans_layer(double*,Locstate,COMP_TYPE*,HYPER_SURF*,
	                              INTERFACE*,INIT_DATA*,int);
LOCAL	void	free_rarefaction_wave_1d_comp_type(COMP_TYPE*);
LOCAL	void	free_stretching_comp_type(COMP_TYPE*);
LOCAL	void	free_trans_layer_comp_type(COMP_TYPE*);
LOCAL   void    free_2d_riemann_comp_type(COMP_TYPE*);
LOCAL	void	get_state_rarefaction_wave_1d(double*,Locstate,
					      COMP_TYPE*,HYPER_SURF*,
					      INTERFACE*,INIT_DATA*,int);
LOCAL	void	get_state_stretching(double*,Locstate,COMP_TYPE*,
				     HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	set_trans_layer_comp_type(COMP_TYPE*,Front*);


/*
*	The program get_comp_type_type() is trying to be smart.
*	An alternative (cleaner) design is to make the code prompt for
*	comp_type->type for each layer. 
*	ymy 7/5/93
*/

/*ARGSUSED*/
EXPORT	COMP_TYPE	*g_get_comp_type_type(
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_PHYSICS	*ip,
	INIT_DATA	*init)
{
	Front	  *front = ip->root->front;
	LAYER	  *lyr = layer_sys->layer[layer_label];
	char	  choice[Gets_BUF_SIZE];
	COMP_TYPE *ct;

	ct = (region_label > 0) ? comp_type(lyr->ellip[region_label]->compin) :
		        	  comp_type(lyr->comp);

#if defined(TWOD) || defined(THREED)
	if ((ALL_CONTACT(layer_sys->flag) == YES) && (is_gravity() == YES))
	    set_rt_kh_comp_type(ct,front);
	else
#endif /* defined(TWOD) || defined(THREED) */
	if ((region_label == 0) &&
	    is_rarefaction_wave(lyr->lower_surf->wv_type) &&
	    is_rarefaction_wave(lyr->upper_surf->wv_type))
	{
	    set_rarefaction_wave_1d_comp_type(ct,front);
	}
#if defined(TWOD) || defined(THREED)
	else if ((lyr->num_ellips == 0) &&
	     is_scalar_wave(lyr->lower_surf->wv_type) &&
	     is_scalar_wave(lyr->upper_surf->wv_type))
	{
	    screen("Is this (the %d%s layer)",
	    	   layer_label,ordinal_suffix(layer_label));
	    screen(" a transition layer? (n, dflt): ");
	    (void) Gets(choice);
	    if ((choice[0] == 'y') || (choice[0] == 'Y'))
	    {
	    	set_trans_layer_comp_type(ct,front);
	    }
	}
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
	else if ((region_label > 0) &&
		 (ELLIPTICAL_REGION(layer_sys->flag) == YES) &&
		 (region_label == lyr->num_ellips))
	{
	    if (lyr->ellip[region_label]->wv_type == NEUMANN_BOUNDARY || 
		lyr->ellip[region_label]->wv_type == MOVABLE_BODY_BOUNDARY)
	    {
	    	set_obstacle_comp_type(ct,front);
	    	return ct;
	    }
	    else
	    	set_elliptical_comp_type(ct,ip);
	}
#endif /* defined(TWOD) */
	if (region_label > 0)
	{
	    (void) sprintf(choice,
	                   "\nEnter the component type for the %d%s region "
	                   "inside the %d%s layer\n",
	    	           region_label,ordinal_suffix(region_label),
		           layer_label,ordinal_suffix(layer_label));
	}
	else
	{
	    (void) sprintf(choice,
	                   "\nEnter the component type for the %d%s layer\n",
	    	           layer_label,ordinal_suffix(layer_label));
	}
	prompt_for_comp_type_type(ct,choice,ip);
	return ct;
}		/*end g_get_comp_type_type*/

EXPORT	void	prompt_for_comp_type_type(
	COMP_TYPE    *ct,
	const char   *mesg,
	INIT_PHYSICS *ip)
{
	Front	  *front = ip->root->front;
	char	   choice[Gets_BUF_SIZE];
	boolean	   default_comp_type_set = NO;
	const char *default_choice;
	int	   i, len;

	screen("%sChoices are\n"
	       "\tAmbient region, (AR)\n"
	       "\tRandom perturbation region, (RA)\n"
	       "\tRayleigh-Taylor, (RT)\n"
	       "\tKelvin-Helmholtz, (KH)\n"
	       "\tOne dimensional rarefaction, (RW1D)\n"
	       "\tOne dimensional overlay, (1DO)\n"
	       "\tTransitional Layer, (TL)\n"
	       "\tElliptical region, (EL)\n"
	       "\tStretching, (ST)\n"
	       "\tTabulated region, (TR)\n"
	       "\tObstacle, (obstacle)\n",mesg);
	default_comp_type_set = YES;
	switch (ct->type)
	{
	case AMBIENT:
	    default_choice = "ar";
	    break;
	case RANDOM_REGION:
	    default_choice = "ra";
	    break;
	case RT_KH:
	    default_choice = "rt";
	    break;
	case RAREFACTION_WAVE_1D:
	    default_choice = "rw1d";
	    break;
	case TRANS_LAYER:
	    default_choice = "tl";
	    break;
	case ELLIPTICAL:
	    default_choice = "el";
	    break;
	case STRETCHING:
	    default_choice = "st";
	    break;
	case ONE_DIMENSIONAL_OVERLAY:
	    default_choice = "1do";
	    break;
	case TABULATED_REGION:
	    default_choice = "tr";
	    break;
	case OBSTACLE:
	    default_choice = "obstacle";
	    break;
	default:
	    default_comp_type_set = NO;
	    default_choice = "none";
	    screen("No default available,  component type unset\n");
	    break;
	}
	screen("Default = %s\n",default_choice);
	screen("Enter choice: ");
	(void) Gets(choice);
	len = (int) strlen(choice);
	for (i = 0; i < len; ++i)
	    choice[i] = tolower(choice[i]);
	if ((choice[0] == '\0') && (default_comp_type_set==NO))
	{
	    screen("No default component type available, enter choice: ");
	    (void) Gets(choice);
	    if (choice[0] == '\0')
	    {
	        screen("ERROR in g_get_comp_type_type(), "
		       "no valid default choice for component type exists\n");
	        clean_up(ERROR);
	    }
	}
	if ((choice[0] != '\0') && (strcmp(choice,default_choice) != 0))
	{
	    if (strcmp(choice,"ar")==0)
	    	set_ambient_comp_type(ct,front);
	    else if (strcmp(choice,"rw1d")==0)
	    	set_rarefaction_wave_1d_comp_type(ct,front);
	    else if (strcmp(choice,"tl")==0)
	    	set_trans_layer_comp_type(ct,front);
	    else if (strcmp(choice,"st")==0)
	    	set_stretching_comp_type(ct,front);
	    else if (strcmp(choice,"obstacle")==0)
	    	set_obstacle_comp_type(ct,front);
	    else if (strcmp(choice,"ra")==0)
	    	set_random_region_comp_type(ct,front);
	    else if (strcmp(choice,"tr")==0)
	    	set_tabulated_region_comp_type(ct,front);
#if defined(ONED)
	    else if (strcmp(choice,"1do")==0)
	    	set_1d_overlay_comp_type(ct);
#endif /* defined(ONED) */
#if defined(TWOD) || defined(THREED)
	    else if ((strcmp(choice,"rt")==0) || (strcmp(choice,"kh")==0))
	    	set_rt_kh_comp_type(ct,front);
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
	    else if (strcmp(choice,"el")==0)
	    	set_elliptical_comp_type(ct,ip);
#endif /* defined(TWOD) */
	    else
	    {
	    	screen("ERROR in g_get_comp_type_type(), "
	    	       "invalid choice (%s) for component type\n",
	    		choice);
	    	clean_up(ERROR);
	    }
	}
}		/*end prompt_for_comp_type_type*/

LOCAL	void	set_trans_layer_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	_TRANS_LAYER	*extra;

	if (comp_type->type == TRANS_LAYER) /*ALREADY SET*/
		return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = TRANS_LAYER;
	scalar(&extra,sizeof(_TRANS_LAYER));
	comp_type->extra = (POINTER)extra;
	alloc_state(front->interf,&extra->lower_st,front->sizest);
	alloc_state(front->interf,&extra->upper_st,front->sizest);

	comp_type->_get_state = get_state_trans_layer;
	comp_type->free_comp_type_extra = free_trans_layer_comp_type;

}		/*end set_trans_layer_comp_type*/


/*ARGSUSED*/
LOCAL	void	get_state_trans_layer(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	int		i, dim = ct->params->dim;
	double		z0, z1, z, q0, q1, frac;
	_TRANS_LAYER	*t_l = Trans_layer(ct);

#define Interpolate(q0, q1, frac)         ((frac)*q0 + (1-(frac))*q1)

	debug_print("init_states","Entered get_state_trans_layer()\n");
	set_type_of_state(s,TGAS_STATE);
	z = coords[dim-1];
	z0 = get_surf_height(coords, t_l->lower_surf);
	z1 = get_surf_height(coords, t_l->upper_surf);
	frac = 1.0 - (z-z0)/(z1-z0);
	frac = min(frac,1.0);
	frac = max(frac,0.0);

	q0 = Dens(t_l->lower_st);	 
	q1 = Dens(t_l->upper_st);
	Dens(s) = Interpolate(q0, q1, frac);

	q0 = pressure(t_l->lower_st);
	q1 = pressure(t_l->upper_st);
	Press(s) = Interpolate(q0, q1, frac);

	for (i = 0; i < dim; ++i)
	{
	    q0 = vel(i, t_l->lower_st);
	    q1 = vel(i, t_l->upper_st);
	    Vel(s)[i] = Interpolate(q0, q1, frac);
	}

	if ((ct->params != Params(t_l->lower_st)) ||
	    (ct->params != Params(t_l->upper_st)))
	{
	    screen("ERROR in get_state_trans_layer(): "
	           "inconsistent params!\n");
	    clean_up(ERROR);
	}
	Init_params(s,ct->params);

	set_state(s,stype,s);
#undef Interpolate
}		/*end get_state_trans_layer*/

LOCAL	void	free_trans_layer_comp_type(
	COMP_TYPE	*comp_type)
{
	_TRANS_LAYER *extra;
	if (comp_type->type != TRANS_LAYER)
		return;

	extra = Trans_layer(comp_type);
	if (extra == NULL)
		return;

	if (extra->lower_st != NULL)
		free(extra->lower_st);
	if (extra->upper_st != NULL)
		free(extra->upper_st);
	free(extra);
	comp_type->extra = NULL;
}		/*end free_trans_layer_comp_type*/


EXPORT	void	set_stretching_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	_STRETCHING	*extra;

	if (comp_type->type == STRETCHING) /*ALREADY SET*/
	    return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = STRETCHING;
	scalar(&extra,sizeof(_STRETCHING));
	comp_type->extra = (POINTER)extra;
	alloc_state(front->interf,&extra->ambient,front->sizest);
	comp_type->_get_state = get_state_stretching;
	comp_type->free_comp_type_extra = free_stretching_comp_type;
}		/*end set_stretching_comp_type*/

EXPORT	void	set_2d_riemann_comp_type(
	COMP_TYPE	*comp_type,
	Gas_param       *params,
	INIT_DATA       *init)
{
	_RIEMANN_2D	*extra;
	int		i,num_sections;
	Locstate	*states;
	double		*start_angle;
	char            message[100];

	screen("Enter number of sections for the 2D Riemann problems: ");
	Scanf("%d\n",&num_sections);

	uni_array(&states,num_sections,sizeof(Locstate));
	uni_array(&start_angle,num_sections,FLOAT);

	for (i = 0; i < num_sections; ++i)
        {
	    (*params->_alloc_state)(&states[i],params->sizest);
	    screen("Enter start angle (in degree between 0 to 360) "
		   "for section %d: ");
	    Scanf("%f\n",&start_angle[i]);
	    start_angle[i] *= PI/180.0;
	    sprintf(message," section %d",i+1);
	    prompt_for_ref_state(message,states[i],TGAS_STATE,params,init);
	    set_state(states[i],GAS_STATE,states[i]);
	}

	scalar(&extra,sizeof(_RIEMANN_2D));
	extra->num_sections = num_sections;
	extra->start_angle = start_angle;
	extra->states = states;
	comp_type->extra = (POINTER)extra;
	comp_type->_get_state = get_state_2d_riemann;
	comp_type->free_comp_type_extra = free_2d_riemann_comp_type;
}		/*end set_2d_riemann_comp_type */

LOCAL   void    free_2d_riemann_comp_type(
        COMP_TYPE       *comp_type)
{
	int i,num_sections;
        Locstate *states;
	double *start_angle;

        if (comp_type->extra != NULL)
        {
	    num_sections = Riemann_2d(comp_type)->num_sections;
	    start_angle = Riemann_2d(comp_type)->start_angle;
	    states = Riemann_2d(comp_type)->states;
	    for (i = 0; i < num_sections; ++i)
		free(states[i]);
	    free(start_angle);
	    free(comp_type->extra);
        }
        comp_type->extra = NULL;
}               /*end free_2d_riemann_comp_type */


/*ARGSUSED*/
LOCAL	void get_state_stretching(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		stype)
{
	_STRETCHING	*str = Stretching(ct);
	double		v[3], f[3], w;
	int		i, j, dim;

	debug_print("init_states","Entered get_state_stretching()\n");
	set_state(state,stype,str->ambient);
	dim = Params(state)->dim;
	for (i = 0; i < dim; ++i)
	{
	    f[i] = (coords[i] - str->L[i])/(str->U[i] - str->L[i]);
	    v[i] = 0.0;
	}
	for (i = 0; i < 1<<dim; ++i)
	{
	    w = 1.0;
	    for (j = 0; j < dim; ++j)
	    	w *= ((i>>j)%2)?f[j]:1.0-f[j];
	    for (j = 0; j < dim; ++j)
	    	v[j] += w*str->v[i][j];
	}
	add_velocity_to_state(state,v);
}		/*end get_state_stretching*/

LOCAL	void	free_stretching_comp_type(
	COMP_TYPE	*comp_type)
{
	_STRETCHING *extra;
	if (comp_type->type != STRETCHING)
		return;

	extra = Stretching(comp_type);
	if (extra == NULL)
	    return;

	if (extra->ambient != NULL)
	    free(extra->ambient);
	free(extra);
	comp_type->extra = NULL;
}		/*end free_stretching_comp_type*/

EXPORT	void	set_rarefaction_wave_1d_comp_type(
	COMP_TYPE	*comp_type,
	Front		*front)
{
	_RAREFACTION_WAVE_1D	*extra;

	if (comp_type->type == RAREFACTION_WAVE_1D) /*ALREADY SET*/
		return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = RAREFACTION_WAVE_1D;
	scalar(&extra,sizeof(_RAREFACTION_WAVE_1D));
	comp_type->extra = (POINTER)extra;
	alloc_state(front->interf,&extra->stl,front->sizest);
	alloc_state(front->interf,&extra->stt,front->sizest);

	comp_type->_get_state = get_state_rarefaction_wave_1d;
	comp_type->free_comp_type_extra = free_rarefaction_wave_1d_comp_type;

}		/*end set_rarefaction_wave_1d_comp_type*/

LOCAL	void	free_rarefaction_wave_1d_comp_type(
	COMP_TYPE	*comp_type)
{
	_RAREFACTION_WAVE_1D *rw1d;
	if (comp_type->type != RAREFACTION_WAVE_1D)
		return;

	rw1d = Rarefaction_wave_1d(comp_type);
	if (rw1d == NULL)
		return;

	if (rw1d->stl != NULL)
		free(rw1d->stl);
	if (rw1d->stt != NULL)
		free(rw1d->stt);
	free(rw1d);
	comp_type->extra = NULL;
}		/*end free_rarefaction_wave_1d_comp_type*/

/*ARGSUSED*/
LOCAL	void	get_state_rarefaction_wave_1d(
	double		*coords,
	Locstate	s,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int             stype)
{
	int		dim = ct->params->dim;
	int		i;
	double		spdnew, mnew;
	Locstate	stl, stt;
	_RAREFACTION_WAVE_1D	*rw1d = Rarefaction_wave_1d(ct);

	debug_print("init_states","Entered get_state_rarefaction_wave_1d()\n");
	if (ct->type != RAREFACTION_WAVE_1D)
	{
	    screen("ERROR in get_state_rarefaction_wave_1d(), "
	           "inconsistent comp_type->type\n");
	    clean_up(ERROR);
	}

	stl = rw1d->stl;
	stt = rw1d->stt;

	if (rarefaction_edge_at_coords(coords,hs,LEADING_EDGE) == YES)
	    set_state(s,TGAS_STATE,stl);
	else if (rarefaction_edge_at_coords(coords,hs,TRAILING_EDGE) == YES)
	    set_state(s,TGAS_STATE,stt);
	else
	{
	    double z, zbar, tbar, zmin, zmax;
	    if (perturbed_rarefaction_layer(rw1d) == YES)
	    {
		double zl, zt;
		zl = get_surf_height(coords,rw1d->lead);
		zt = get_surf_height(coords,rw1d->trail);
		if (zl < zt)
		{
		    zmin = zl;
		    zmax = zt;
		}
		else
		{
		    zmin = zt;
		    zmax = zl;
		}
	        tbar = rw1d->tbar;
		zbar = 0.5*(zl + zt - (rw1d->spl + rw1d->spt)*tbar);
	    }
	    else
	    {
	        zmin = rw1d->zmin;
	        zmax = rw1d->zmax;
	        tbar = rw1d->tbar;
		zbar = rw1d->zbar;
	    }

	    z = coords[dim-1];
	    if (z <= zmin)
	        set_state(s,TGAS_STATE,rw1d->stmin);
	    else if (z >= zmax)
	        set_state(s,TGAS_STATE,rw1d->stmax);
	    else
	    {
	        double speed, vl;
		vl = Vel(stl)[dim-1];
	        speed = (z - zbar)/tbar;
	        (void) oned_state_in_rarefaction_fan(speed,vl,stl,stt,s,
						     TGAS_STATE,&spdnew,
						     rw1d->l_or_r);
	    }
	}
	mnew = Vel(s)[0];
	for (i = 0; i < dim-1; ++i)
	    Vel(s)[i] = Vel(stl)[i];
	Vel(s)[dim-1] = mnew;
	set_state(s,stype,s);
}		/*end get_state_rarefaction_wave_1d*/

LOCAL	boolean	perturbed_rarefaction_layer(
	_RAREFACTION_WAVE_1D *rw1d)
{
	if (rw1d->lead == NULL || rw1d->trail == NULL)
	    return NO;
	return ((perturbed_layer_surf(rw1d->lead) == YES) ||
			(perturbed_layer_surf(rw1d->trail) == YES)) ? YES : NO;
}		/*end perturbed_rarefaction_layer*/

LOCAL   boolean perturbed_layer_surf(
	LAYER_SURF *l)
{
	return (l->fpoly != NULL) ? YES : NO;
}               /*end perturbed_ellipsoid*/

