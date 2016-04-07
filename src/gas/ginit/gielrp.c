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
*				gielrp.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains initialization routines for elliptical Riemann
*	problems in gas dynamics.
*
*/

#include <ginit/ginit.h>

enum _ELLIPTICAL_REGION_TYPE {
    AMBIENT_ELLIPTICAL_REGION                           = 0,
    RANDOM_AMBIENT_ELLIPTICAL_REGION,
    ELLIPTICAL_VELOCITY_REGION,
    RANDOM_ELLIPTICAL_VELOCITY_REGION,
    RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION,
    RANDOM_RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION,
    ELLIPTICAL_RAREFACTION_REGION,
    ONE_DIMENSIONAL_ELLIPTICAL_OVERLAY,
    TABULATED_ELLIPTICAL_REGION,
    READ_TGAS_STATE_FROM_1D_FILE_REGION
};
typedef enum _ELLIPTICAL_REGION_TYPE ELLIPTICAL_REGION_TYPE;

	/* LOCAL Function Prototypes */

LOCAL	COMPONENT prompt_for_elliptical_boundary_comps(ELLIPSOID*,COMPONENT,
						       Front*,const char*,
						       const char*,const char*);
LOCAL	COMPONENT prompt_for_elliptical_component(COMPONENT,const char*,
						  const char*,const char*,
						  const char*);
LOCAL	boolean	add_ray_to_boundary(double*,double*,int,const char*,const char*,
	                            const char*);
LOCAL	boolean	inside_box(const double*,const double*,const double*,int);
LOCAL	boolean    inside_ellipsoid(double*,ELLIPSOID*);
LOCAL	boolean    nested_ellipsoids(ELLIPSOID*,ELLIPSOID*);
LOCAL	boolean	prompt_for_elliptic_region_type(const INIT_PHYSICS*,
						const INIT_DATA*,
						ELLIPTICAL_REGION_TYPE*,
						ELLIPSOID*,ELLIPSOID*);
LOCAL	boolean    perturbed_ellipsoid(ELLIPSOID*);
LOCAL	boolean    perturbed_elliptical_rarefaction_region(_RAREFACTION_WAVE_1D*);
LOCAL	boolean	rarefaction_region(ELLIPSOID*,ELLIPSOID*);
LOCAL	double	mean_pressure_for_region(COMP_TYPE*,ELLIPSOID*,ELLIPSOID*,
                                         double,double,RECT_GRID*);
LOCAL	double	prompt_for_wave_number(double,int,double,double,const char*);
LOCAL	int	read_vector(double*,double*,int);
LOCAL	void	Gas_from_axis_elliptical_state(Locstate,COMP_TYPE*);
LOCAL	void	axis_elliptical_state_from_Gas(COMP_TYPE*,Locstate);
LOCAL	void	connect_ellipsoid_to_bdry(ELLIPSOID*,Front*,const char*,
					  INIT_PHYSICS*);
LOCAL	void	free_elliptical_comp_type(COMP_TYPE*);
LOCAL	void	get_state_elliptical(double*,Locstate,COMP_TYPE*,
				     HYPER_SURF*,INTERFACE*,INIT_DATA*,int);
LOCAL	void	init_central_gravity_states(INIT_PHYSICS*,INIT_DATA*,double*,
					   COMPONENT*,Gas_param**,
					   ELLIPSOID**,int,double,double);
LOCAL	void	init_elrp_in_to_out_states(INIT_PHYSICS*,INIT_DATA*,double*,
					   COMPONENT*,Gas_param**,
					   ELLIPSOID**,int,double,double);
LOCAL	void	init_elrp_out_to_in_states(INIT_PHYSICS*,INIT_DATA*,double*,
					   COMPONENT*,Gas_param**,
					   ELLIPSOID**,int,double,double);
LOCAL	void	init_elrp_params(INIT_DATA*,INIT_PHYSICS*,
				 Gas_param**,const char*,int*,int);
LOCAL	void	init_elrp_states(INIT_PHYSICS*,INIT_DATA*,double*,COMPONENT*,
				 int*,Gas_param**,ELLIPSOID**,int,double,double);
LOCAL	void	init_input_modes(int,int,int,int,FOURIER_POLY*,double*,double*);
LOCAL	void	init_state_behind_elliptic_wave(COMPONENT,COMPONENT,
						Gas_param*,Gas_param*,
						ELLIPSOID*,ELLIPSOID*,
						ELLIPSOID*,double*,
						INIT_PHYSICS*,INIT_DATA*,
						const char*,const char*);
LOCAL	void	make_plane_wave(double,COMPONENT,COMPONENT,double,double,
				int,RECT_GRID*);
LOCAL	void	prompt_for_elliptical_state(INIT_PHYSICS*,COMPONENT,
					    ELLIPSOID*,ELLIPSOID*,
					    ELLIPSOID*,Gas_param*,
					    const char*,INIT_DATA*);
LOCAL	double	prompt_for_elrp_plane_wave(const char*,double*,COMP_TYPE*,
					   COMP_TYPE*,Front*,Gas_param*,
					   int,INIT_DATA*,INIT_PHYSICS*,
					   double,double);
LOCAL	void	prompt_for_rarefaction_edge(_RAREFACTION_WAVE_1D*,
	                                    RAREFACTION_EDGE_TYPE,
					    ELLIPSOID*,ELLIPSOID*,const char*);
LOCAL	void	prompt_for_rotation(const char*,double**,int);
LOCAL	void	set_consistent_ellip_orient(ELLIPSOID*,double,double);
LOCAL	void	set_up_ambient_elliptic_region(boolean,COMP_TYPE*,Gas_param*,
					    const char*,Front*,INIT_DATA*);
LOCAL   void    set_up_read_state_from_1d_file_region(COMP_TYPE*,Gas_param*,
                                            const char*,Front*,INIT_DATA*);
LOCAL	void set_up_elliptical_overlay_region(COMP_TYPE*,Gas_param*,const char*,
	                                      INIT_PHYSICS*,INIT_DATA*);
LOCAL	void	input_spherical_modes(int,int,int,int,FOURIER_POLY*,
				      double*,double*);
LOCAL   void    get_state_1dtable(double*,Locstate,COMP_TYPE*,HYPER_SURF*,
					INTERFACE*,INIT_DATA*,int);
LOCAL   void    free_1dtable_comp_type(COMP_TYPE*);
LOCAL 	void 	init_SN_layers(Front*,double,int,int,double*,double*,
				double*,double*,Locstate);


/*
*	The various regions in this initialization are referred to
*	by indices as indicated in the table below.
*
*	Inside the inner ellipsoid		0
*	Between inner and outer ellipsoids	1
*	Outside the outer ellipsoid		2
*	Above the upper plane wave		3
*	Below the lower plane wave		4
*	Behind the ellipsoidal wall		5
*/

EXPORT	void init_el_riem_prob(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	COMPONENT   min_comp = FIRST_DYNAMIC_COMPONENT;
	Front	   *front = ip->root->front;
	INTERFACE   *intfc = front->interf;
	char	   s[Gets_BUF_SIZE];
	char	   mesg[256];
	const char *elname;
	char	   where[256];
	ELLIPSOID  **ellip;
	double	   upper_height = HUGE_VAL, lower_height = -HUGE_VAL;
	double	   *L = front->rect_grid->L;
	double	   *h = front->rect_grid->h;
	double	   dflt_rad[3];
	double	   *pr_region;
	int	   i, k, dim = front->interf->dim;
	int	   num_ellip;
	int	   *wv_type;
	COMPONENT  *comps;
	Gas_param  **prms;
	LAYER_FLAG Flag;

	CLEAR_LAYER_FLAG(Flag);
	ALLOW_SECTORED_ELLIPSOIDS(Flag) = YES;
	if (dim < 2)
	{
	    screen("ERROR in init_el_riem_prob(), dim < 2 not supported\n");
	    clean_up(ERROR);
	    return;
	}

	elname = (dim == 3) ? "ellipsoid" : "ellipse";
	num_ellip = 1;
	screen("Enter the number (>= 1) of %s waves (dflt = %d): ",
		(dim==2)?"elliptical":"ellipsoidal",num_ellip);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&num_ellip);
	if (num_ellip < 1)
	{
	    screen("ERROR in init_el_riem_prob(), no waves specified\n");
	    clean_up(ERROR);
	    return;
	}
	uni_array(&comps,num_ellip+3,sizeof(COMPONENT));
	uni_array(&prms,num_ellip+3,sizeof(Gas_param*));
	uni_array(&ellip,num_ellip,sizeof(ELLIPSOID*));
	uni_array(&wv_type,num_ellip+3,INT);
	uni_array(&pr_region,num_ellip+3,FLOAT);

	for (i = 0; i < num_ellip+3; ++i)
	{
	    wv_type[i] = UNKNOWN_WAVE_TYPE;
	    pr_region[i] = -HUGE_VAL;
	    comps[i] = NO_COMP;
	}

	    /* Choose and locate tracked waves */

	(void) sprintf(where,"%s",(num_ellip > 1) ? "inner" : "");
	INNERMOST_ELLIPSOID(Flag) = YES;
	ellip[0] = prompt_for_ellipsoid(front,where,NULL,NULL,&Flag,ip);
	INNERMOST_ELLIPSOID(Flag) = NO;
	wv_type[0] = ellip[0]->wv_type;
	if ((wv_type[0] < FIRST_PHYSICS_WAVE_TYPE) && (num_ellip > 1))
	{
	    ellip[0]->compin = comps[0] = COMPOBST;
	    (void)SetConstantFlowRegion(comps[0],return_obst_state(),intfc);
	    set_obstacle_comp_type(comp_type(comps[0]),front);
	}
	else
	{
	    ellip[0]->compin = comps[0] =
	        prompt_for_elliptical_component(min_comp,"inside",
						where,elname,"");
	    (void) sprintf(mesg,"inside the %s",where);
	    comps[0] = prompt_for_elliptical_boundary_comps(ellip[0],comps[0],
						       front,mesg,elname,"");
	}
	for (i = 1; i < num_ellip; ++i)
	{
	    (void) sprintf(mesg,"%d-%s",i,ordinal_suffix(i));
	    for (k = 0; k < dim; ++k)
	    	dflt_rad[k] = ellip[i-1]->rad[k] + 0.5*h[k];
	    ellip[i] = prompt_for_ellipsoid(front,mesg,ellip[i-1]->cen,
					    dflt_rad,&Flag,ip);
	    (void) sprintf(mesg,"%d-%s and the %d-%s",
			   i-1,ordinal_suffix(i-1),i,ordinal_suffix(i));
	    if (comps[i] == NO_COMP)
		comps[i] = comps[i-1]+1;
	    ellip[i-1]->compout = ellip[i]->compin = comps[i] =
	        prompt_for_elliptical_component(comps[i],"between",
						mesg,elname,"s");
	    (void) sprintf(mesg,"%d-%s and the %d-%s",
			   i-1,ordinal_suffix(i-1),i,ordinal_suffix(i));
	    (void) sprintf(mesg,"between the %d-%s and the %d-%s",
				i-1,ordinal_suffix(i-1),i,ordinal_suffix(i));
	    comps[i] = prompt_for_elliptical_boundary_comps(ellip[i],comps[i],
						       front,mesg,elname,"s");
	    wv_type[i] = ellip[i]->wv_type;
	    if (nested_ellipsoids(ellip[i-1],ellip[i]) == YES)
	    {
	        outer_ellipsoid(ellip[i-1]) = ellip[i];
	        inner_ellipsoid(ellip[i]) = ellip[i-1];
	    }
	    if (wv_type[i] < FIRST_PHYSICS_WAVE_TYPE)
	    {
	    	num_ellip = i+1;
	    	break;
	    }
	}
	(void) sprintf(mesg,"%d-%s",num_ellip-1,ordinal_suffix(num_ellip-1));
	if (wv_type[num_ellip-1] < FIRST_PHYSICS_WAVE_TYPE)
	{
	    ellip[num_ellip-1]->compout = comps[num_ellip] = COMPOBST;
	    (void)SetConstantFlowRegion(comps[num_ellip],return_obst_state(),
					intfc);
	    set_obstacle_comp_type(comp_type(comps[num_ellip]),front);
	}
	else
	{
	    ellip[num_ellip-1]->compout = comps[num_ellip] =
	        prompt_for_elliptical_component(comps[num_ellip-1]+1,"outside",
					    mesg,elname,"");
	    comps[num_ellip+1] = comps[num_ellip+2] = comps[num_ellip];
	    if (dim != 3)
	    	connect_ellipsoid_to_bdry(ellip[num_ellip-1],front,"outer ",ip);
	    screen("To have upper plane wave, enter the height above L[%d]: ",
		   dim-1);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&upper_height);
	        upper_height += L[dim-1];
	    	wv_type[num_ellip] = prompt_for_wave_type(
						"for the upper plane wave",
						front->interf,ip);
		comps[num_ellip+1] =
		    prompt_for_elliptical_component(comps[num_ellip]+1,
						    "above","upper",
						    "plane wave","");
	    }
	    else
		comps[num_ellip+1] = comps[num_ellip];
	    screen("To have lower plane wave, enter the height above L[%d]: ",
		   dim-1);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	    	(void) sscan_float(s,&lower_height);
	        lower_height += L[dim-1];
	    	wv_type[num_ellip+1] = prompt_for_wave_type(
						"for the lower plane wave",
						front->interf,ip);
		comps[num_ellip+2] =
		    prompt_for_elliptical_component(comps[num_ellip+1]+1,
						    "below","lower",
						    "plane wave","");
	    }
	    else
		comps[num_ellip+2] = comps[num_ellip];
	}



		/* EOS for elliptical regions */

	init_elrp_params(init,ip,prms,elname,wv_type,num_ellip);

	init_elrp_states(ip,init,pr_region,comps,wv_type,prms,ellip,
	                 num_ellip,lower_height,upper_height);
	screen("\n");

	/* NOTE:  assumes explosion or implosion into vr = 0 */
	/* TODO: revise radial_analysis so shock is not assumed
	*  to propagate first. */

	for (i = 0; i < num_ellip; ++i)
	{
	    HYPER_SURF *hs;

	    if (is_scalar_wave(ellip[i]->wv_type))
		ellip[i]->layer_index = ++num_layers(front->interf);
	    set_consistent_ellip_orient(ellip[i],pr_region[i],pr_region[i+1]);
	    hs = make_ellipsoid(ellip[i],comps[i],comps[i+1],front);
	    if (wv_type[i] < FIRST_PHYSICS_WAVE_TYPE)
	    {	
	    	(void) sprintf(s,"%s %d",elname,i);
		bstate_index(hs) = prompt_for_boundary_state(wave_type(hs),
                                       s,ellip[i]->cen,negative_component(hs),
                                       -1,hs,init,ip);
	    }
	    make_ellip_region_boundaries(ellip[i],front);
	}

	if (wv_type[num_ellip-1] >= FIRST_PHYSICS_WAVE_TYPE)
	{
	    if (wv_type[num_ellip] != UNKNOWN_WAVE_TYPE) 
	    {
	    	make_plane_wave(upper_height,comps[num_ellip+1],
	    		        comps[num_ellip],pr_region[num_ellip+1],
	    		        pr_region[num_ellip],
	    		        wv_type[num_ellip],front->rect_grid);
	    	screen("Is the flow above the upper plane wave "
	    	       "constant (dflt = no): ");
	    	(void) Gets(s);
	    	if (s[0] == 'y' || s[0] == 'Y')
	    	    (void)SetConstantFlowRegion(comps[num_ellip+1],
	    			Ambient(comp_type(comps[num_ellip+1])),intfc);
	    }
	    if (wv_type[num_ellip+1] != UNKNOWN_WAVE_TYPE) 
	    {
	    	make_plane_wave(lower_height,comps[num_ellip],
	    		        comps[num_ellip+2],pr_region[num_ellip],
			        pr_region[num_ellip+2],
			        wv_type[num_ellip+1],front->rect_grid);
		screen("Is the flow below the lower plane wave "
		       "constant (dflt = no): ");
		(void) Gets(s);
		if (s[0] == 'y' || s[0] == 'Y')
		    (void)SetConstantFlowRegion(comps[num_ellip+2],
		    		Ambient(comp_type(comps[num_ellip+2])),intfc);
	    }
	}
	free_these(5,comps,ellip,wv_type,prms,pr_region);
}		/*end init_el_riem_prob*/

LOCAL	COMPONENT prompt_for_elliptical_component(
	COMPONENT  comp,
	const char *side,
	const char *mesg,
	const char *elname,
	const char *plural)
{
	char s[Gets_BUF_SIZE];
	screen("Enter the component number for the region\n\t%s "
	       "the %s %s%s (default = %d): ",side,mesg,elname,plural,comp);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%d",&comp);
	return comp;
}		/*end prompt_for_elliptical_component*/

LOCAL	COMPONENT prompt_for_elliptical_boundary_comps(
	ELLIPSOID  *ellip,
	COMPONENT  comp,
	Front      *front,
	const char *mesg,
	const char *elname,
	const char *plural)
{
	INTERFACE  *intfc = front->interf;
	char	   mesg1[256];
	const char *sn[] = {"lower","upper"};
	const char *an[] = {"angular","azmuth","polar"};
	int        k, dim = front->rect_grid->dim;

	for (k = 0; k < 4; ++k)
	{
	    if (ellip->btype[k] == UNKNOWN_BOUNDARY_TYPE)
		continue;
	    if (ellip->btype[k] < FIRST_PHYSICS_WAVE_TYPE)
	    {
		ellip->bcomp[k] = COMPOBST;
	        (void)SetConstantFlowRegion(COMPOBST,return_obst_state(),intfc);
	        set_obstacle_comp_type(comp_type(COMPOBST),front);
	    }
	    else
	    {
	        (void) sprintf(mesg1,"%s %s side of the region %s",
			       sn[k%2],an[(dim-2)*(1+k/2)],mesg);
		ellip->bcomp[k] = comp =
			prompt_for_elliptical_component(comp+1,"bounding",
						        mesg1,elname,plural);
	    }
	}
	return comp;
}		/*end prompt_for_elliptical_boundary_comps*/

LOCAL	void init_elrp_states(
	INIT_PHYSICS *ip,
	INIT_DATA    *init,
	double	     *pr_region,
	COMPONENT    *comps,
	int	     *wv_type,
	Gas_param    **prms,
	ELLIPSOID    **ellip,
	int	     num_ellip,
	double        lower_height,
	double        upper_height)
{
	Front	  *front = ip->root->front;
	RECT_GRID *gr = front->rect_grid;
	int	  i, dim = front->rect_grid->dim;
	double	  nor[MAXD];
	char	  s[Gets_BUF_SIZE];

	screen("Choose order of initialization, possible types are\n"
		"\tCentral gravity (from inner to outer) (g), \n"
		"\tFrom inner to outer (i), \n"
		"\tFrom outer to inner (o), \n"
		"Enter choice: ");
	(void) Gets(s);

	switch (s[0])
	{
	case 'g':
	case 'G':
	    init_central_gravity_states(ip,init,pr_region,comps,
			               prms,ellip,num_ellip,
				       lower_height,upper_height);
	    break;
	case 'i':
	case 'I':
	    init_elrp_in_to_out_states(ip,init,pr_region,comps,
			               prms,ellip,num_ellip,
				       lower_height,upper_height);
	    break;
	case 'o':
	case 'O':
	    init_elrp_out_to_in_states(ip,init,pr_region,comps,
			               prms,ellip,num_ellip,
				       lower_height,upper_height);
	    break;
	default: 
	    screen("Error: no such choice!\n");
	    clean_up(ERROR);
	}


		/* Prompt for states along plane waves */

	if (wv_type[num_ellip] != UNKNOWN_WAVE_TYPE) 
	{
	    for (i = 1; i < dim; ++i)
		nor[i-1] = 0.0;
	    nor[dim-1] = -1.0;
	    pr_region[num_ellip+1] =
	        prompt_for_elrp_plane_wave("upper",nor,
		                           comp_type(comps[num_ellip]),
					   comp_type(comps[num_ellip+1]),
					   front,prms[num_ellip+1],
					   wv_type[num_ellip],init,ip,
					   upper_height,gr->U[dim-1]);
	}
	else
	    pr_region[num_ellip+1] = pr_region[num_ellip];

	if (wv_type[num_ellip+1] != UNKNOWN_WAVE_TYPE) 
	{
	    for (i = 1; i < dim; ++i)
		nor[i-1] = 0.0;
	    nor[dim-1] = 1.0;
	    pr_region[num_ellip+2] =
	        prompt_for_elrp_plane_wave("lower",nor,
				           comp_type(comps[num_ellip]),
				           comp_type(comps[num_ellip+2]),
				           front,prms[num_ellip+2],
				           wv_type[num_ellip+1],init,ip,
					   gr->L[dim-1],lower_height);
	}
	else
	    pr_region[num_ellip+2] = pr_region[num_ellip];
}		/*end init_elrp_states*/

LOCAL	void init_elrp_out_to_in_states(
	INIT_PHYSICS *ip,
	INIT_DATA    *init,
	double	     *pr_region,
	COMPONENT    *comps,
	Gas_param    **prms,
	ELLIPSOID    **ellip,
	int	     num_ellip,
	double        lower_height,
	double        upper_height)
{
	Front	   *front = ip->root->front;
	RECT_GRID  *gr = front->rect_grid;
	INTERFACE  *intfc = front->interf;
	COMP_TYPE  *ct;
	int	   i, dim = gr->dim;
	const char *elname = (dim == 3) ? "ellipsoid" : "ellipse";
	char	   where[256];
	char	   s[Gets_BUF_SIZE];

	if (prms[num_ellip] == Params(return_obst_state()))
	{
	    pr_region[num_ellip] = 0.0;
	    (void)SetConstantFlowRegion(comps[num_ellip],
					return_obst_state(),intfc);
	    set_obstacle_comp_type(comp_type(comps[num_ellip]),front);
	    (void) sprintf(where," inside the %s%s",
			   (num_ellip > 1) ? "outer " : "",elname);
	    prompt_for_elliptical_state(ip,comps[num_ellip-1],
					ellip[num_ellip-1],
					(num_ellip > 1) ?
					    ellip[num_ellip-2] : NULL,
					ellip[num_ellip-1],
				        prms[num_ellip-1],where,init);
	    if (num_ellip == 1)
		return;
	    num_ellip--;
	}
	else
	{
	    (void) sprintf(where," outside the %s%s",
			   (num_ellip > 1) ? "outer " : "",elname);
	    prompt_for_elliptical_state(ip,comps[num_ellip],
					ellip[num_ellip-1],
					ellip[num_ellip-1],NULL,
				        prms[num_ellip],where,init);
	    ct = comp_type(comps[num_ellip]);
	    pr_region[num_ellip] =
	        mean_pressure_for_region(ct,ellip[num_ellip-1],
		                         NULL,lower_height,upper_height,gr);
	}

	for (i = num_ellip-1; i > 0; i--)
	{
	    (void) sprintf(where," between the %d-%s and the %d-%s %ss",
		           i,ordinal_suffix(i),i-1,ordinal_suffix(i-1),elname);
	    if (ellip[i]->wv_type == BACKWARD_SHOCK_WAVE)
		(void) strcpy(s,"y");
	    else
	    {
	        screen("Type 'y' to enter an arbitrary state "
	               "for the region\n\t%s: ",where);
	        (void) Gets(s);
	    }
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	        prompt_for_elliptical_state(ip,comps[i],
					     ellip[i],ellip[i-1],ellip[i],
					     prms[i],where,init);
		pr_region[i] =
		    mean_pressure_for_region(comp_type(comps[i]),
		                             ellip[i-1],ellip[i],
					     lower_height,upper_height,gr);
	    }
	    else
	        init_state_behind_elliptic_wave(comps[i+1],comps[i],
						prms[i+1],prms[i],
						ellip[i],ellip[i-1],ellip[i],
						pr_region+i,ip,init,
						"inside",where);
	}
	if (prms[0] == Params(return_obst_state()))
	{
	    pr_region[0] = 0.0;
	    (void)SetConstantFlowRegion(comps[0],return_obst_state(),intfc);
	    set_obstacle_comp_type(comp_type(comps[0]),front);
	}
	else
	{
	    (void) sprintf(where," inside the inner %s%s",
			   (num_ellip > 1) ? "inner " : "",elname);
	    if (ellip[i]->wv_type == BACKWARD_SHOCK_WAVE)
		(void) strcpy(s,"y");
	    else
	    {
	        screen("Type 'y' to enter an arbitrary state "
	               "for the region\n\t%s: ",where);
	        (void) Gets(s);
	    }
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	        prompt_for_elliptical_state(ip,comps[0],
					    ellip[0],
					    NULL,ellip[0],
					    prms[0],where,init);
		pr_region[0] =
		    mean_pressure_for_region(comp_type(comps[0]),
		                             NULL,ellip[0],
					     lower_height,upper_height,gr);
	    }
	    else
	        init_state_behind_elliptic_wave(comps[1],comps[0],prms[1],
						prms[0],ellip[0],NULL,ellip[0],
						pr_region+0,ip,init,
						"inside",where);
	}
}		/*end init_elrp_out_to_in_states*/

LOCAL	void init_central_gravity_states(
	INIT_PHYSICS *ip,
	INIT_DATA    *init,
	double	     *pr_region,
	COMPONENT    *comps,
	Gas_param    **prms,
	ELLIPSOID    **ellip,
	int	     num_ellip,
	double        lower_height,
	double        upper_height)
{
	Front	   *front = ip->root->front;
	RECT_GRID  *gr = front->rect_grid;
	int	   i,imin,imax,num_intervals;
	char	   s[Gets_BUF_SIZE];
	double	   rho_c,pres_c,temp_c;
	double	   range,r0,radius,dr;
	double	   *den,*pre,*grv,*pos;
	Locstate   st_start;
	double	   *GL = gr->GL;
	double	   *GU = gr->GU;
	int	   dim = gr->dim;
	COMP_TYPE  *ctype;
	ONED_TABLE *table;

	alloc_state(front->interf,&st_start,front->sizest);
	Params(st_start) = prms[0];

	screen("There are two types of input for the core states\n"
		"\tCore density and pressure (p),\n"
		"\tCore density and temperature (t),\n"
		"Enter choice: ");
	Gets(s);
	switch (s[0])
	{
	case 'p':
	case 'P':
	    screen("Enter the core density and pressure: ");
	    (void) Scanf("%f %f\n",&rho_c,&pres_c);
	    Dens(st_start) = rho_c;
	    Press(st_start) = pres_c;
	    set_type_of_state(st_start,TGAS_STATE);
	    printf("In init_central_gravity_states() pressure = %f\n",
			    pressure(st_start));
	    reset_gamma(st_start);
	    break;
	case 't':
	case 'T':
	    screen("Enter the core density and temperature: ");
	    (void) Scanf("%f %f\n",&rho_c,&temp_c);
	    Dens(st_start) = rho_c;
	    Temperature(st_start) = temp_c;
	    set_type_of_state(st_start,FGAS_STATE);
	    reset_gamma(st_start);
	    Press(st_start) = pressure(st_start);
	    set_type_of_state(st_start,TGAS_STATE);
	    reset_gamma(st_start);
	    break;
	default:
	    screen("Error: unknown choice!\n");
	    clean_up(ERROR);
	}

	range = 0.0;
	for (i = 0; i < dim; ++i)
	    range += sqr(GU[i] - GL[i]);
	range = sqrt(range);
	screen("The range of data is from 0 to %g\n",range);
	screen("Enter the number of intervals for the data: ");
	(void) Scanf("%d\n",&num_intervals);

	scalar(&table,sizeof(ONED_TABLE));
	table->nx = num_intervals+1;
	uni_array(&table->d,num_intervals+1,FLOAT);
	uni_array(&table->p,num_intervals+1,FLOAT);
	uni_array(&table->v,num_intervals+1,FLOAT);
	uni_array(&table->pt,num_intervals+1,FLOAT);

	uni_array(&grv,num_intervals+1,FLOAT);
	pos = table->pt;
	den = table->d;
	pre = table->p;
	for (i = 0; i < num_intervals; ++i)
	    table->v[i] = 0.0;

	pos[0] = grv[0] = 0.0;
	pre[0] = Press(st_start);
	den[0] = Dens(st_start);
	imin = 0;
	dr = range/num_intervals;

	for (i = 0; i < num_ellip; ++i)
	{
#if defined(COMBUSTION_CODE)
	    screen("Type 'y' if the gas in %d-th layer has burned: ",i);
	    Gets(s);
	    if (s[0] == 'y')
	    	Params(st_start) = prms[i]->other_params;
	    else
#endif /* defined(COMBUSTION_CODE) */
	    	Params(st_start) = prms[i];
	    radius = ellip[i]->rad[0];
	    imax = (int)(radius/dr);
	    init_SN_layers(front,dr,imin,imax,pos,den,pre,grv,st_start);
	    r0 = pos[imax];
	    imin = imax;
	    Dens(st_start) = den[imax];
	    Press(st_start) = pre[imax];
	    ctype = comp_type(comps[i]);
	    ctype->type = ONE_DIMENSIONAL_TABLE;
	    ctype->extra = (POINTER) table;
	    ctype->params = Params(st_start);
            ctype->_get_state = get_state_1dtable;
	    if (i == 0)
            	ctype->free_comp_type_extra = free_1dtable_comp_type;
	    else
            	ctype->free_comp_type_extra = NULL;
	}
#if defined(COMBUSTION_CODE)
	screen("Type 'y' if the gas in %d-th layer has burned: ",num_ellip);
	Gets(s);
	if (s[0] == 'y')
	    Params(st_start) = prms[num_ellip]->other_params;
	else
#endif /* defined(COMBUSTION_CODE) */
	    Params(st_start) = prms[num_ellip];
	imax = num_intervals;
	init_SN_layers(front,dr,imin,imax,pos,den,pre,grv,st_start);
	ctype = comp_type(comps[num_ellip]);
	ctype->type = ONE_DIMENSIONAL_TABLE;
	ctype->extra = (POINTER) table;
	ctype->params = Params(st_start);
        ctype->_get_state = get_state_1dtable;
        ctype->free_comp_type_extra = NULL;
	free_these(2,grv,st_start);
}		/*end init_central_gravity_states*/

LOCAL	void init_elrp_in_to_out_states(
	INIT_PHYSICS *ip,
	INIT_DATA    *init,
	double	     *pr_region,
	COMPONENT    *comps,
	Gas_param    **prms,
	ELLIPSOID    **ellip,
	int	     num_ellip,
	double        lower_height,
	double        upper_height)
{
	Front	   *front = ip->root->front;
	RECT_GRID  *gr = front->rect_grid;
	INTERFACE  *intfc = front->interf;
	int	   i, dim = gr->dim;
	int	   istart;
	const char *elname = (dim == 3) ? "ellipsoid" : "ellipse";
	char	   where[256];
	char	   s[Gets_BUF_SIZE];

	if (prms[0] == Params(return_obst_state()))
	{
	    pr_region[0] = 0.0;
	    (void)SetConstantFlowRegion(comps[0],return_obst_state(),intfc);
	    set_obstacle_comp_type(comp_type(comps[0]),front);
	    (void) sprintf(where," outside the %s%s",
	                   (num_ellip > 1) ? "inner " : "",elname);
	    prompt_for_elliptical_state(ip,comps[1],
					ellip[0],
					ellip[0],(num_ellip > 1) ?
					    ellip[1] : NULL,
				        prms[1],where,init);
	    if (num_ellip == 1)
		return;
	    istart = 2;
	}
	else
	{
	    (void) sprintf(where," inside the %s%s",
			   (num_ellip > 1) ? "inner " : "",elname);
	    prompt_for_elliptical_state(ip,comps[0],ellip[0],NULL,ellip[0],
					prms[0],where,init);
	    pr_region[0] =
	        mean_pressure_for_region(comp_type(comps[0]),NULL,ellip[0],
		                         lower_height,upper_height,gr);
	    istart = 1;
	}
		
	for (i = istart; i < num_ellip; ++i)
	{
	    (void) sprintf(where," between the %d-%s and the %d-%s %ss",
		           i-1,ordinal_suffix(i-1),i,ordinal_suffix(i),elname);
	    if (ellip[i-1]->wv_type == FORWARD_SHOCK_WAVE)
		(void) strcpy(s,"y");
	    else
	    {
	        screen("Type 'y' to enter an arbitrary state "
	               "for the region\n\t%s: ",where);
	        (void) Gets(s);
	    }
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	        prompt_for_elliptical_state(ip,comps[i],
					    ellip[i],ellip[i-1],ellip[i],
					    prms[i],where,init);
		pr_region[i] =
		    mean_pressure_for_region(comp_type(comps[i]),
		                             ellip[i-1],ellip[i],
					     lower_height,upper_height,gr);
	    }
	    else
	        init_state_behind_elliptic_wave(comps[i-1],comps[i],prms[i-1],
						prms[i],
						ellip[i-1],ellip[i-1],ellip[i],
						pr_region+i,ip,init,
						"outside",where);
	}
	if (prms[num_ellip] == Params(return_obst_state()))
	{
	    pr_region[0] = 0.0;
	    (void)SetConstantFlowRegion(comps[num_ellip],
					return_obst_state(),intfc);
	    set_obstacle_comp_type(comp_type(comps[num_ellip]),front);
	}
	else
	{
	    i = num_ellip;
	    (void) sprintf(where," outside the %s%s",
			   (num_ellip > 1) ? "outer " : "",elname);
	    if (ellip[i-1]->wv_type == FORWARD_SHOCK_WAVE)
		(void) strcpy(s,"y");
	    else
	    {
	        screen("Type 'y' to enter an arbitrary state "
	               "for the region\n\t%s: ",where);
	        (void) Gets(s);
	    }
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	        prompt_for_elliptical_state(ip,comps[i],
					    ellip[i-1],ellip[i-1],NULL,
					    prms[i],where,init);
		pr_region[i] =
		    mean_pressure_for_region(comp_type(comps[i]),
		                             ellip[i-1],NULL,
					     lower_height,upper_height,gr);
	    }
	    else
	        init_state_behind_elliptic_wave(comps[i-1],comps[i],prms[i-1],
						prms[i],ellip[i-1],ellip[i-1],
						NULL,pr_region+i,ip,init,
						"outside",where);
	}
}		/*end init_elrp_in_to_out_states*/

LOCAL	double	mean_pressure_for_region(
	COMP_TYPE *ct,
	ELLIPSOID *inner,
	ELLIPSOID *outer,
	double     l,
	double     u,
	RECT_GRID *gr)
{
	switch (ct->type)
	{
	case AMBIENT:
	    return pressure(Ambient(ct));
	case RANDOM_REGION:
	    return pressure(Mean(Random_state(ct)));
	case ELLIPTICAL:
	    return pressure(Elliptical(ct)->state);
	case ONE_DIMENSIONAL_OVERLAY:
	    return -HUGE_VAL;/*TODO Compute this value?*/
	case ONE_DIMENSIONAL_TABLE:
	    return HUGE_VAL;
	case TABULATED_REGION:
	    {
	        TABULATED_REGION_DATA *table = Tabulated_region_data(ct);
		double p;
		double L[3], U[3];
		double *h;
		double x[3];
		int   *n;
	        int   i, j, N, indx;
		int   dim = gr->dim;

		for (i = 0; i < dim; ++i)
		{
		    L[i] = gr->L[i];
		    U[i] = gr->U[i];
		}
		L[dim-1] = max(L[dim-1],l);
		U[dim-1] = min(U[dim-1],u);
		h = table->h;
		n = table->n;
	        for (p = 0.0, indx = 0, N = 0; indx < table->len; ++indx)
		{
		    for (i = indx, j = 0; j < dim; ++j)
		    {
			x[j] = table->GL[j] + h[j]*(i%n[j]);
			i /= n[j];
		    }
		    if (((inner != NULL) &&  inside_ellipsoid(x,inner)) ||
		        ((outer != NULL) && !inside_ellipsoid(x,outer)) ||
			(!inside_box(x,L,U,dim)))
		        continue;
	            p += table->p[indx];
		    ++N;
		}
		if (N > 0)
		    p /= N;
		return p;
	    }
	default:
	    screen("ERROR in mean_pressure_for_region(), "
		   "invalid component type\n");
	    clean_up(ERROR);
	}
	return -HUGE_VAL;
}		/*end mean_pressure_for_region*/

LOCAL	boolean	inside_box(
	const double *x,
	const double *l,
	const double *u,
	int   dim)
{
	int i;

	for (i = 0; i < dim; ++i)
	{
	    if (x[i] < l[i])
	        return NO;
	    if (x[i] > u[i])
	        return NO;
	}
	return YES;
}		/*end inside_box*/

/*ARGSUSED*/
LOCAL	void	init_state_behind_elliptic_wave(
	COMPONENT	acomp,
	COMPONENT	bcomp,
	Gas_param	*aprms,
	Gas_param	*bprms,
	ELLIPSOID	*ellipsoid,
	ELLIPSOID	*inner,
	ELLIPSOID	*outer,
	double		*pbehind,
	INIT_PHYSICS	*ip,
	INIT_DATA	*init,
	const char	*mesg,
	const char      *where)
{
	_ELLIPTICAL	       *ellip;
	_RAREFACTION_WAVE_1D   *rw1d = NULL;
	ELLIPTICAL_REGION_TYPE ert;
	GRAVITY                *grav_data = gravity_data(init);
	COMP_TYPE	       *ctype;
	Front		       *front = ip->root->front;
	INTERFACE	       *intfc = front->interf;
	boolean		       random;
	double		       nor[MAXD];
	int		       wv_type = ellipsoid->wv_type;
	int		       i, dim = front->rect_grid->dim;
	Locstate	       ahead, behind;
	size_t		       sizest = front->sizest;
	char		       s[Gets_BUF_SIZE];
#if defined(COMBUSTION_CODE)
	const char	       *elname = (dim == 3) ? "ellipsoid" : "ellipse";
#endif /* defined(COMBUSTION_CODE) */

	if (comp_type(acomp)->type == ONE_DIMENSIONAL_OVERLAY)
	{
	    prompt_for_elliptical_state(ip,bcomp,ellipsoid,inner,outer,bprms,
					where,init);
	    return;
	}

	alloc_state(intfc,&ahead,sizest);
	alloc_state(intfc,&behind,sizest);
	Gas_from_axis_elliptical_state(ahead,comp_type(acomp));
	ellip = Elliptical(comp_type(acomp));
	if (ellip->stratification_type != CONSTANT)
	{
	    double dr, r = 0.0;
	    for (i = 0; i < dim; ++i)
	        r = max(r,ellipsoid->rad[i]);
	    dr = r - ellip->r0;
	    get_state_in_stratified_region(ellip->stratification_type,behind,
					   dr,grav_data->G,ahead);
	    set_state(ahead,state_type(ahead),behind);
	    clear_state(intfc,behind,sizest);
	}

	nor[0] = (is_backward_wave(wv_type)) ? -1.0 : 1.0;
	for (i = 1; i < dim; ++i)
	    nor[i] = 0.0;

	random = prompt_for_elliptic_region_type(ip,init,&ert,inner,outer);
	printf("ert = %d\n",ert);
	if (ert == ONE_DIMENSIONAL_ELLIPTICAL_OVERLAY)
	{
	    set_up_elliptical_overlay_region(comp_type(bcomp),bprms,where,
					     ip,init);
	}
	else if (ert == TABULATED_ELLIPTICAL_REGION)
	{
	    prompt_for_elliptical_state(ip,bcomp,ellipsoid,inner,outer,
	    				bprms,s,init);
	}
	else
	{
	    if (is_shock_wave(wv_type))
	    {
#if defined(COMBUSTION_CODE)
	        if (aprms->composition_type != PURE_NON_REACTIVE)
	        {
	            (void) sprintf(s," %s the %s",mesg,elname);
	            prompt_for_elliptical_state(ip,bcomp,ellipsoid,inner,outer,
	    				    aprms,s,init);
	    	    Gas_from_axis_elliptical_state(behind,comp_type(bcomp));
	        }
	        else
#endif /* defined(COMBUSTION_CODE) */
	    	prompt_for_behind_shock_state(ahead,behind,YES,nor,
	    				      GAS_STATE,YES,init);
	    }
	    else if (is_scalar_wave(wv_type))
	    {
#if defined(COMBUSTION_CODE)
	        if (bprms->composition_type != PURE_NON_REACTIVE)
	        {
	            (void) sprintf(s," %s the %s",mesg,elname);
	            prompt_for_elliptical_state(ip,bcomp,ellipsoid,inner,outer,
	    				        bprms,s,init);
	    	    Gas_from_axis_elliptical_state(behind,comp_type(bcomp));
	        }
	        else
#endif /* defined(COMBUSTION_CODE) */
	            prompt_for_behind_contact_state(ahead,behind,
	    					bprms,YES,nor,GAS_STATE,init);
	    }
	    else if ((ert == ELLIPTICAL_RAREFACTION_REGION) &&
	             (rarefaction_region(inner,outer) == YES))
	    {
	        rw1d = allocate_RAREFACTION_WAVE_1D(front);
	        if (is_rarefaction_leading_edge(wv_type))
	        {
	    	    copy_state(rw1d->stl,ahead);
	            prompt_for_rarefaction_edge(rw1d,TRAILING_EDGE,
	    				        inner,outer,mesg);
	            set_state(behind,GAS_STATE,rw1d->stt);
	        }
	        else
	        {
	    	    copy_state(rw1d->stt,ahead);
	            prompt_for_rarefaction_edge(rw1d,LEADING_EDGE,
	    				        inner,outer,mesg);
	            set_state(behind,GAS_STATE,rw1d->stl);
	        }
	    }
	    else
	    {
	        set_state(behind,GAS_STATE,ahead);
	    }

	    *pbehind = pressure(behind);
	    ctype = comp_type(bcomp);
	    if (ert == AMBIENT_ELLIPTICAL_REGION)
	    {
	        set_ambient_comp_type(ctype,front);
	        ft_assign(Ambient(ctype),behind,sizest);
	        screen("Is the flow inside the wave constant (dflt = no): ");
	        (void) Gets(s);
	        if (s[0] == 'y' || s[0] == 'Y')
	        	(void)SetConstantFlowRegion(bcomp,Ambient(ctype),intfc);
	    }
	    else if(ert == READ_TGAS_STATE_FROM_1D_FILE_REGION)
            {
                set_up_read_state_from_1d_file_region(ctype,bprms,mesg,
				front,init);
            }
	    else if (ert == RANDOM_AMBIENT_ELLIPTICAL_REGION)
	    {
	        set_random_region_comp_type(ctype,front);
	        init_random_state_region(behind,Random_state(ctype),front);
	    }
	    else
	    {
	        double maxr;

	        set_elliptical_comp_type(ctype,ip);
	        ctype->params = Params(behind);
	        axis_elliptical_state_from_Gas(ctype,behind);
	        ellip = Elliptical(ctype);
	        ellip->ellipsoid = ellipsoid;
	        ellip->rw1d = rw1d;
	        for (maxr = 0.0, i = 0; i < dim; ++i)
	            maxr = max(maxr,ellipsoid->rad[i]);
	        for (i = 0; i < dim; ++i)
	            ellip->weight[i] = ellipsoid->rad[i] / maxr;
	        if ((ert==RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION) ||
	            (ert==RANDOM_RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION))
	        {
	            ellip->stratification_type = prompt_for_stratification("");
	            if (ellip->stratification_type != CONSTANT)
	            {
	                ellip->r0 = maxr;
	    	    screen("Enter the reference radius for the stratification "
	    	           "(dflt = %g): ",ellip->r0);
	    	    (void) Gets(s);
	    	    if (s[0] != '\0')
	    	        (void) sscan_float(s,&ellip->r0);
	            }
	        }
	        if (random == YES)
	        {
	            ellip->rstate = allocate_random_state_structure(intfc);
	            alloc_state(intfc,ellip->wkstate,sizest);
	            alloc_state(intfc,ellip->wkstate+1,sizest);
	            init_random_state_region(ellip->state,ellip->rstate,front);
	        }
	        if (ellip->rw1d != NULL)
	        {
	            verbose_print_state("Initialized elliptical state "
	    			    "for rarefaction region",ellip->state);
	            print_rarefaction_wave_1d(ellip->rw1d);
	        }
	        else
	            verbose_print_state("Initialized elliptical state",
	    			    ellip->state);
	    }
	}
	free_these(2,ahead,behind);
}		/* end init_state_behind_elliptic_wave*/

LOCAL	void init_elrp_params(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	Gas_param	**prms,
	const char	*elname,
	int		*wv_type,
	int		num_ellip)
{
	int	j;
	char	mesg[80];
	char	s[Gets_BUF_SIZE];

	(void) prompt_for_eos_params(init,ip,YES,"");
	if (wv_type[0] < FIRST_PHYSICS_WAVE_TYPE)
	{
	    if (num_ellip > 1)
	        prms[0] = Params(return_obst_state());
	    else
	    {
		screen("Is the region inside the %s active (dflt = y): ",
		       elname);
		(void) Gets(s);
		if ((s[0] == 'N') || (s[0] == 'n'))
		{
	            (void) sprintf(mesg,"\n\tfor the region outside the %s",
			           elname);
		    prms[0] = Params(return_obst_state());
	            prms[1] = prompt_for_eos_params(init,ip,YES,mesg);
		}
		else
		{
	            (void) sprintf(mesg,"\n\tfor the region inside the %s",
			           elname);
	            prms[0] = prompt_for_eos_params(init,ip,YES,mesg);
		    prms[1] = Params(return_obst_state());
		}
	    }
	}
	else
	{
	    (void) sprintf(mesg,"\n\tfor the region inside the %s%s",
	                   (num_ellip > 1) ? "inner " : "",elname);
	    prms[0] = prompt_for_eos_params(init,ip,YES,mesg);
	}
	for (j = 1; j < num_ellip; ++j)
	{
	    if (is_scalar_wave(wv_type[j-1]) ||
		(wv_type[j-1] < FIRST_PHYSICS_WAVE_TYPE))
	    {
	        (void) sprintf(mesg,"\n\tfor the region between the "
				    "%d-%s and %d-%s %ss",
				    j-1,ordinal_suffix(j-1),
				    j,ordinal_suffix(j),elname);
	        prms[j] = prompt_for_eos_params(init,ip,YES,mesg);
	    }
	    else
	        prms[j] = prms[j-1];
	}
	if ((wv_type[num_ellip-1]<FIRST_PHYSICS_WAVE_TYPE) && (num_ellip>1))
	    prms[num_ellip] = Params(return_obst_state());
	else if (is_scalar_wave(wv_type[num_ellip-1]))
	{
	    (void) sprintf(mesg,"\n\tfor the region outside the %s%s",
			   (num_ellip > 1) ? "outer " : "",elname);
	    prms[num_ellip] = prompt_for_eos_params(init,ip,YES,mesg);
	}
	else
	    prms[num_ellip] = prms[num_ellip-1];
	screen("\n");
}		/*end init_elrp_params*/

LOCAL	void	make_plane_wave(
	double		height,
	COMPONENT	comp0,
	COMPONENT	comp1,
	double		p0,
	double		p1,
	int		wv_type,
	RECT_GRID	*rgr)
{
	HYPER_SURF	*hs = NULL;
	CURVE		*cur;
	double		coords1[MAXD], coords2[MAXD];
	double		*L = rgr->L, *U = rgr->U;
	int		i, dim = rgr->dim;

	for (i = 0; i < dim-1; ++i)
	{
	    coords1[i] = L[i];
	    coords2[i] = U[i];
	}
	coords1[dim-1] = height;
	coords2[dim-1] = height;
	switch (dim)
	{
	case 2:
	    cur = make_curve(comp0,comp1,make_node(Point(coords1)),
			     make_node(Point(coords2)));
	    hs = Hyper_surf(cur);
	    start_status(cur) = INCIDENT;
	    end_status(cur) = INCIDENT;
	    wave_type(cur) = wv_type;
	    if ((wv_type == FORWARD_SHOCK_WAVE && p0 < p1) ||
		(wv_type == BACKWARD_SHOCK_WAVE && p0 > p1))
	    {
		invert_curve(cur);
		wave_type(cur) = wv_type;
	    }
	    break;
	case 3:
	    /* TODO */
	    screen("ERROR in make_plane_wave(), 3D code needed\n");
	    clean_up(ERROR);
	}
	if (is_scalar_wave(wv_type))
	    layer_index(hs) = ++num_layers(hs->interface);
}		/*end make_plane_wave*/

LOCAL	void	set_consistent_ellip_orient(
	ELLIPSOID	*ellip,
	double		pin,
	double		pout)
{
	double		tmp;
	int		i, dim = ellip->dim, wv_type = ellip->wv_type;

	if (pin == pout)
	    return;
	if ((wv_type == FORWARD_SHOCK_WAVE && pin < pout) ||
	    (wv_type == BACKWARD_SHOCK_WAVE && pin > pout))
	{
	    if (ellip->nor_orient != NEGATIVE_ORIENTATION)
	    {
	    	ellip->nor_orient = NEGATIVE_ORIENTATION;
	    	for (i = 0; i < dim-1; ++i)
	    	{
	    	    tmp = ellip->ThetaS[i];
	    	    ellip->ThetaS[i] = ellip->ThetaE[i];
	    	    ellip->ThetaE[i] = tmp;
	    	}
	    }
	}
	if ((wv_type == FORWARD_SHOCK_WAVE && pout < pin) ||
			(wv_type == BACKWARD_SHOCK_WAVE && pout > pin))
	{
	    if (ellip->nor_orient != POSITIVE_ORIENTATION)
	    {
	    	ellip->nor_orient = POSITIVE_ORIENTATION;
	    	for (i = 0; i < dim-1; ++i)
	    	{
	    	    tmp = ellip->ThetaS[i];
	    	    ellip->ThetaS[i] = ellip->ThetaE[i];
	    	    ellip->ThetaE[i] = tmp;
	    	}
	    }
	}
}		/*end set_consistent_ellip_orient*/

LOCAL	double prompt_for_elrp_plane_wave(
	const char	*mesg,
	double		*nor,
	COMP_TYPE	*actype,
	COMP_TYPE	*ctype,
	Front		*front,
	Gas_param	*prms,
	int		wv_type,
	INIT_DATA	*init,
	INIT_PHYSICS    *ip,
	double           lower_height,
	double           upper_height)
{
	RECT_GRID *gr = front->rect_grid;
	char	  s[120];
	Locstate  ahead;

	if (mesg != NULL)
	    (void) sprintf(s," for the %s wave",mesg);
	else
	    (void) sprintf(s," for the wave");

	alloc_state(front->interf,&ahead,front->sizest);
	prompt_for_comp_type_type(ctype,s,ip);
	ctype->params = prms;
	if (ctype->type == TABULATED_REGION)
	{
	    set_up_read_state_from_file_region(ctype,s,front);
	    return mean_pressure_for_region(ctype,NULL,NULL,
	                                    lower_height,upper_height,gr);
	}

	Gas_from_axis_elliptical_state(ahead,actype);

	switch(wv_type) 
	{
	case FORWARD_SHOCK_WAVE:
	case BACKWARD_SHOCK_WAVE:
	    prompt_for_behind_shock_state(ahead,Ambient(ctype),YES,
	                                  nor,GAS_STATE,YES,init);
	    break;
	case CONTACT:
	    prompt_for_behind_contact_state(ahead,Ambient(ctype),prms,YES,nor,
					    GAS_STATE,init);
	    break;
	default:
	    screen("No such wave type in init type%s\n",mesg);
	    clean_up(ERROR);
	    break;
	}
	free(ahead);
	return mean_pressure_for_region(ctype,NULL,NULL,
	                                lower_height,upper_height,gr);
}		/*end prompt_for_elrp_plane_wave*/

EXPORT  ELLIPSOID *prompt_for_ellipsoid(
	Front	         *front,
	const char       *mesg,
	double	         *dflt_center,
	double	         *dflt_radii,
	const LAYER_FLAG *flag,
	INIT_PHYSICS     *ip)
{
	IMPORT	boolean suppress_prompts;
	ELLIPSOID  *ellipsoid = NULL;
	ELLIPSOID  Ellip;
	INTERFACE  *intfc = front->interf;
	RECT_GRID  *gr = front->rect_grid;
	boolean    prompt_for_boundaries;
	double      ang;
	char	   s[Gets_BUF_SIZE];
	int	   i, dim = gr->dim;
	int        btype;
	const char *elname, *in_out;

	if (flag == NULL)
	{
	    static LAYER_FLAG Noflag;
	    flag = &Noflag;
	}
	if (dim < 2)
	{
	    screen("ERROR in prompt_for_ellipsoid(), dim = %d < 2 "
		   "not supported\n");
	    clean_up(ERROR);
	}
	elname = (dim == 3) ? "ellipsoid" : "ellipse";

	set_default_ellipsoid_structure(&Ellip,dim);
	Ellip.nor_orient = POSITIVE_ORIENTATION;

	prompt_for_boundaries = NO;
	screen("The %s %s can be chosen to as\n"
	       "\tclosed (full) (C, F, default),\n",mesg,elname);
	if (allow_sectored_ellipsoids(flag) == YES)
	{
	    screen("\thalf (H),\n");
	    screen("\tquarter (Q),");
	    if (dim == 3)
	        screen("\n\teighth (E),");
	    screen(" or\n\tangle sector specified (A).\n");
	}
	else
	{
	    if (dim == 2)
	    {
	        screen("\thalf (H), or\n\tquarter (Q).\n");
	    }
	    else
	    {
	        screen("\thalf (H),\n"
	               "\tquarter (Q), or\n"
	               "\teighth (E).\n");
	    }
	}
	screen("Enter choice: ");
	(void) Gets(s);
	switch (tolower(s[0]))
	{
	case 'a':
	    if (allow_sectored_ellipsoids(flag) != YES)
	    {
		screen("ERROR in prompt_for_ellipsoid(),  sectored "
		       "ellipsoids not allowed for this initialization\n");
		clean_up(ERROR);
	    }
	    Ellip.closed = NO;
	    prompt_for_boundaries = YES;
	    if (dim == 2)
		screen("Enter the low angle for the sector: ");
	    else
	        screen("Enter the lower angular value of "
		       "the azimuth in degrees: ");
	    (void) Scanf("%f\n",Ellip.ThetaS);
	    if (dim == 2)
		screen("Enter the high angle for the sector: ");
	    else
	        screen("Enter the upper angular value of the "
		       "azimuth in degrees: ");
	    (void) Scanf("%f\n",Ellip.ThetaE);
	    if (dim == 3)
	    {
	        screen("Enter the lower polar angle in degrees: ");
	        (void) Scanf("%f\n",Ellip.ThetaS+1);
	        screen("Enter the upper polar angle in degrees: ");
	        (void) Scanf("%f\n",Ellip.ThetaE+1);
	    }

	    Ellip.ThetaS[0] = radians(Ellip.ThetaS[0]);
	    Ellip.ThetaE[0] = radians(Ellip.ThetaE[0]);
	    if (Ellip.nor_orient == POSITIVE_ORIENTATION)
	    {
		if (Ellip.ThetaE[0] < Ellip.ThetaS[0])
		{
		    ang = Ellip.ThetaS[0];
		    Ellip.ThetaS[0] = Ellip.ThetaE[0];
		    Ellip.ThetaE[0] = ang;

		    btype = Ellip.btype[0];
		    Ellip.btype[0] = Ellip.btype[1];
		    Ellip.btype[1] = btype;
		}
		Ellip.ThetaS[0] = normalized_angle(Ellip.ThetaS[0]);
		Ellip.ThetaE[0] = normalized_angle(Ellip.ThetaE[0]);
		while (Ellip.ThetaE[0] < Ellip.ThetaS[0])
		    Ellip.ThetaE[0] += 2.0*PI;
	    }
	    else if (Ellip.nor_orient == NEGATIVE_ORIENTATION)
	    {
		if (Ellip.ThetaS[0] < Ellip.ThetaE[0])
		{
		    ang = Ellip.ThetaS[0];
		    Ellip.ThetaS[0] = Ellip.ThetaE[0];
		    Ellip.ThetaE[0] = ang;

		    btype = Ellip.btype[0];
		    Ellip.btype[0] = Ellip.btype[1];
		    Ellip.btype[1] = btype;
		}
		Ellip.ThetaS[0] = normalized_angle(Ellip.ThetaS[0]);
		Ellip.ThetaE[0] = normalized_angle(Ellip.ThetaE[0]);
		while (Ellip.ThetaE[0] > Ellip.ThetaS[0])
		    Ellip.ThetaE[0] -= 2.0*PI;
	    }
	    if (dim == 3)
	    {
		while (Ellip.ThetaS[1] < 0.0)
		    Ellip.ThetaS[1] += 360.0;
		while (Ellip.ThetaS[1] > 360.0)
		    Ellip.ThetaS[1] -= 360.0;
		if (Ellip.ThetaS[1] > 180.0)
		    Ellip.ThetaS[1] = 360.0 - Ellip.ThetaS[1];
	        Ellip.ThetaS[1] = radians(Ellip.ThetaS[1]);
		while (Ellip.ThetaE[1] < 0.0)
		    Ellip.ThetaE[1] += 360.0;
		while (Ellip.ThetaE[1] > 360.0)
		    Ellip.ThetaE[1] -= 360.0;
		if (Ellip.ThetaE[1] > 180.0)
		    Ellip.ThetaE[1] = 360.0 - Ellip.ThetaE[1];
	        Ellip.ThetaE[1] = radians(Ellip.ThetaE[1]);
		if (Ellip.ThetaE[1] < Ellip.ThetaS[1])
		{
		    ang = Ellip.ThetaS[1];
		    Ellip.ThetaS[1] = Ellip.ThetaE[1];
		    Ellip.ThetaE[1] = ang;
		}
	    }
	    for (i = 0; i < 4; ++i)
		Ellip.rbdry[i] = NO;
	    break;
	case 'h':
	    Ellip.ThetaS[0] = -0.5*PI;
	    Ellip.ThetaE[0] = 0.5*PI;
	    if (dim == 3)
	    {
	    	Ellip.ThetaS[1] = 0.0;
	    	Ellip.ThetaE[1] = PI;
	    }
	    Ellip.closed = NO;
	    for (i = 0; i < 4; ++i)
		Ellip.rbdry[i] = YES;
	    break;
	case 'q':
	    Ellip.closed = NO;
	    Ellip.ThetaS[0] = 0.0;
	    Ellip.ThetaE[0] = 0.5*PI;
	    if (dim == 3)
	    {
	    	Ellip.ThetaS[1] = 0.0;
	    	Ellip.ThetaE[1] = PI;
	    }
	    for (i = 0; i < 4; ++i)
		Ellip.rbdry[i] = YES;
	    break;
	case 'e':
	    Ellip.closed = NO;
	    Ellip.ThetaS[0] = 0.0;
	    Ellip.ThetaE[0] = 0.5*PI;
	    if (dim == 3)
	    {
	    	Ellip.ThetaS[1] = 0.0;
	    	Ellip.ThetaE[1] = 0.5*PI;
	    }
	    for (i = 0; i < 4; ++i)
		Ellip.rbdry[i] = YES;
	    break;
	case 'f':
	case 'c':
	default:
	    Ellip.closed = YES;
	    Ellip.ThetaS[0] = -PI;
	    Ellip.ThetaE[0] =  PI;
	    if (dim == 3)
	    {
	    	Ellip.ThetaS[1] = 0.0;
	    	Ellip.ThetaE[1] = PI;
	    }
	    for (i = 0; i < 4; ++i)
		Ellip.rbdry[i] = NO;
	    break;
	}
	if (dim == 2)
	    screen("This %s has been oriented in the counter clockwise direction\n"
	             "so that the forward moving normal directions point towards\n"
	             "the outside of the %s.\n",elname,elname);

	screen("Enter the %d coordinates of the center of the %s %s",
	       dim,mesg,elname);
	if (dflt_center != NULL)
	{
	    static const char *dflt = "\n\t(dflt = ";
	    print_general_vector(dflt,dflt_center,dim,"): ");
	    if (suppress_prompts == NO)
	        fprint_general_vector(stderr,dflt,dflt_center,dim,"): ");
	}
	else
	    screen(": ");
	if (read_vector(Ellip.cen,dflt_center,dim) != dim)
	{
	    screen("ERROR in prompt_for_ellipsoid(), insufficient "
		   "number of coordinates for %s center\n",elname);
	    clean_up(ERROR);
	}

	screen("Enter the common radius or the %d mean radii of the %s %s",
	       dim,mesg,elname);
	if (dflt_radii != NULL)
	{
	    static const char *dflt = "\n\t(dflt = ";
	    print_general_vector(dflt,dflt_radii,dim,"): ");
	    if (suppress_prompts == NO)
	        fprint_general_vector(stderr,dflt,dflt_radii,dim,"): ");
	}
	else
	    screen(": ");
	if (read_vector(Ellip.rad,dflt_radii,dim) == 0)
	{
	    screen("ERROR in prompt_for_ellipsoid(), insufficient "
		   "number of coordinates for %s radii\n",elname);
	    clean_up(ERROR);
	}

	screen("Enter a scaling factor for the radius (dflt = %g): ",
	       Ellip.scale);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&Ellip.scale);

	prompt_for_rotation(mesg,Ellip.Q,dim);

	switch (dim)
	{
	case 2:
	    screen("Type 'y' if the %s %s is to be fourier polynomial "
	           "perturbed: ",mesg,elname);
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
	    {
	        FOURIER_POLY *fpoly;
	        int          num_modes_i = 0, num_modes_r = 0;
	        int          min_n_i, max_n_i;
	        int          min_n_r, max_n_r;

	        screen("Should the perturbation be ");
	        screen("(R)andom (default), (I)nput, or\n\t");
	        screen("(M)ixed input and random?: ");
	        (void) Gets(s);
	        switch (s[0])
	        {
	        case 'm':
	        case 'M':
	    	    (void) sprintf(s,"for the user input modes on this %s",
				   elname);
	    	    num_modes_i = random_bubble_num_modes(s,&min_n_i,&max_n_i,
							  dim);
	    	    (void) sprintf(s,"for the random modes on this %s",elname);
	    	    num_modes_r = random_bubble_num_modes(s,&min_n_r,&max_n_r,
							  dim);
	            ellipsoid = allocate_ellipsoid(&Ellip,
						   num_modes_i+num_modes_r);
		    fpoly = ellipsoid->fpoly;
	    	    init_input_modes(0,min_n_i,max_n_i,num_modes_i,fpoly,
				     ellipsoid->ThetaS,ellipsoid->ThetaE);
	    	    init_random_modes(num_modes_i,min_n_r,max_n_r,num_modes_r,
	        	              fpoly,ellipsoid->ThetaS,
				      ellipsoid->ThetaE);
	    	    break;
	        case 'i':
	        case 'I':
	    	    (void) sprintf(s,"for this %s",elname);
	    	    num_modes_i = random_bubble_num_modes(s,&min_n_i,
							  &max_n_i,dim);
	            ellipsoid = allocate_ellipsoid(&Ellip,num_modes_i);
		    fpoly = ellipsoid->fpoly;
	    	    init_input_modes(0,min_n_i,max_n_i,num_modes_i,fpoly,
	    	                     ellipsoid->ThetaS, ellipsoid->ThetaE);
	    	    break;
	        case 'r':
	        case 'R':
	        default:
	    	    (void) sprintf(s,"for this random %s",elname);
	    	    num_modes_r = random_bubble_num_modes(s,&min_n_r,
							  &max_n_r,dim);
	            ellipsoid = allocate_ellipsoid(&Ellip,num_modes_r);
		    fpoly = ellipsoid->fpoly;
	    	    init_random_modes(0,min_n_r,max_n_r,num_modes_r,
	    	                      fpoly,ellipsoid->ThetaS,
				      ellipsoid->ThetaE);
	    	    break;
	        }
	    }
	    break;
	case 3:
	    screen("Type 'y' if the %s %s is to be perturbed by spherical "
	           "harmonics: ",mesg,elname);
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
	    {
	        int num_modes_i = 0;
	        int min_n_i, max_n_i;

	    	(void) sprintf(s,"for this %s",elname);
	    	num_modes_i = spherical_num_modes(s,&min_n_i,&max_n_i,
						  Ellip.ThetaE);
	    	(void) sprintf(s,"for the user input modes on this %s",elname);
	        ellipsoid = allocate_ellipsoid(&Ellip,num_modes_i);
	    	input_spherical_modes(0,min_n_i,max_n_i,num_modes_i,
		                      ellipsoid->fpoly,ellipsoid->ThetaS,
				      ellipsoid->ThetaE);
	    }
	    break;
	}

	if (ellipsoid == NULL)
	    ellipsoid = allocate_ellipsoid(&Ellip,0);

	if (dim == 2)
	{
	    screen("Type 'y' if the %s %s is to be Legendre polynomial "
	           "perturbed: ",mesg,elname);
	    (void) Gets(s);
	    if (s[0] == 'Y' || s[0] == 'y')
	    {
	        (void) sprintf(s,"for this Legendre perturbed %s",elname);
	        ellipsoid->lpoly = get_legendre_coeffs(0.0,s);
	    }
	}

	if (all_contact(flag) == YES)
	    ellipsoid->wv_type = CONTACT;
	else
	{
	    in_out = (ellipsoid->nor_orient == POSITIVE_ORIENTATION) ?
		     "outward" : "inward";
	    (void) sprintf(s, "for the %s %s (normal/forward direction = %s)",
	    	       mesg, elname,in_out);
	    ellipsoid->wv_type = prompt_for_wave_type(s,front->interf,ip);
	}

	if (dim == 3)
	    (void) sprintf(s,"for the %s ellipsoidal surface ",mesg);
	else
	    (void) sprintf(s,"for the %s elliptical curve ",mesg);

	ellipsoid->surf_tension =
	    prompt_for_surface_tension(ellipsoid->wv_type,s);
	if (ellipsoid->wv_type == MOVABLE_BODY_BOUNDARY)
	    prompt_for_rigid_body_params(&ellipsoid->rgb_params,s,dim);

	ellipsoid->untracked = NO;
	if (ellipsoid->wv_type >= FIRST_PHYSICS_WAVE_TYPE)
	{
	    screen("Type y to turn off tracking for the %s %s: ", mesg, elname);
	    (void) Gets(s);
	    if ( s[0] == 'y' || s[0] == 'Y' )
	        ellipsoid->untracked = YES;
	}
	if ((prompt_for_boundaries == YES) &&
	    (ellipsoid->untracked != YES) &&
	    (    (ellipsoid->wv_type >= FIRST_PHYSICS_WAVE_TYPE) ||
	         (innermost_ellipsoid(flag) == NO)    ))
	{
	    int        n = 1<<(dim-1), k = 2*(dim-2);
	    const char *bn[] = {"the lower angular curve",
		                "the upper angular curve",
		                "the lower azmuth surface",
		                "the upper azmuth surface",
		                "the lower polar surface",
		                "the upper polar surface"};

	    for (i = 0; i < n; ++i)
	        ellipsoid->btype[i] = prompt_for_wave_type(bn[i+k],intfc,ip);
	}
	 
#if DONT_COMPILE
	if ((ellipsoid->untracked == NO) && (ellipsoid->fpoly == NULL) &&
	    (ellipsoid->lpoly == NULL))
	{
	    double r = ellipsoid->rad[0];
	    int i;
	    for (i = 0; i < dim; ++i)
		if ((ellipsoid->cen[i] != 0.0) || (ellipsoid->rad[i] != r))
		    break;
	    if (i == dim)
	    {
		screen("Do you wish to enforce origin pointing "
		       "normal uni_arrays (dflt = no): ");
		(void) Gets(s);
		if (strcasecmp(s,"yes") == 0)
		    ellipsoid->n2o_enforced = YES;
	    }
	}
#endif /* DONT_COMPILE */
	return ellipsoid;
}		/*end prompt_for_ellipsoid*/

LOCAL	void	connect_ellipsoid_to_bdry(
	ELLIPSOID    *ellip,
	Front        *front,
	const char   *mesg,
	INIT_PHYSICS *ip)
{
	BDRY_SIDE  bside;
	INTERFACE  *intfc = front->interf;
	RECT_GRID  *gr = front->rect_grid;
	boolean    add_ray[4];
	double      p[3];
	double      theta[4][2];
	double      tol = grid_tolerance(gr);/*TOLERANCE*/
	int        i, j, n, dim = gr->dim;
	int        sd;
	const char *elname;
	static const char *nn[] = { "start",
			            "end",
			            "lower azmuth lower polar",
			            "upper azmuth lower polar",
			            "lower azmuth upper polar",
			            "upper azmuth upper polar"};
	static const char *sn[] = { "for the lower angle curve",
			            "for the upper angle curve",
			            "for the lower azmuth surface",
			            "for the upper azmuth surface",
			            "for the lower polar surface",
			            "for the upper polar surface"};
	static const int side[4][4][2] = { 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 1, 1, 0, 0, 0, 0,
	                                   2, 0, 2, 1, 0, 3, 1, 3};
        static const int av[4][4][2] = {0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 1, 1, 0, 0, 0, 0,
                                        1, 2, 0, 3, 0, 3, 1, 2};
	

	if (ellip->closed == YES)
	    return;

	elname = (dim == 3) ? "ellipsoid" : "ellipse";

	for (i = 0; i < 2; ++i)
	{
	    theta[2*i][0] = ellip->ThetaS[0];
	    theta[2*i+1][0] = ellip->ThetaE[0];

	    theta[i][1] = ellip->ThetaS[1];
	    theta[i+2][1] = ellip->ThetaE[1];
	    add_ray[2*i] = add_ray[2*i+1] = NO;
	}

	n = 1<<(dim-1);

	for (i = 0; i < n; ++i)
	{
	    coords_on_pert_ellipsoid(theta[i],p,ellip);
	    bside = boundary_side(p,gr,tol);
	    if (bside != NOT_A_BDRY)
		continue;
	    add_ray[i] = add_ray_to_boundary(ellip->odir[i],theta[i],
					     dim,nn[2*(dim-2)+i],mesg,elname);
	    if (add_ray[i] != YES)
		continue;
	    for (j = 0; j < 2; ++j)
	    {
		sd = side[dim][i][j];
	        if ((add_ray[av[dim][i][j]] != YES) ||
		    (ellip->obtype[sd] != UNKNOWN_BOUNDARY_TYPE))
		    continue;
	        ellip->obtype[sd] = prompt_for_wave_type(sn[2*(dim-2)+sd],
					                     intfc,ip);
		if (ellip->obtype[sd] == UNKNOWN_BOUNDARY_TYPE)
		    continue;
		if (ellip->obtype[sd] < FIRST_PHYSICS_WAVE_TYPE)
		{
		    ellip->obcomp[sd] = COMPOBST;
	            (void)SetConstantFlowRegion(ellip->obcomp[sd],
						return_obst_state(),intfc);
	            set_obstacle_comp_type(comp_type(ellip->obcomp[sd]),front);
		}
		else
		{
		    ellip->obcomp[sd] =
			prompt_for_elliptical_component(ellip->compout,
						        "outside",mesg,elname,
							"");
	        }
	    }
	}
}		/*end connect_ellipsoid_to_bdry*/

LOCAL	boolean	add_ray_to_boundary(
	double	   *odir,
	double      *theta,
	int        dim,
	const char *sd,
	const char *mesg,
	const char *elname)
{
	char  s[Gets_BUF_SIZE];
	double ang[2];
	const char *aname[] = {"","","azmuth and polar "};
	const char *plural[] = {"","","s"};

	ang[0] = theta[0];
	ang[1] = theta[1];
	screen("To install a front from the %s node of the %s%s "
	       "to the\n\tcomputational boundary,  enter the %sdirection "
	       "angle%s (degrees) of the ray\n\tto the boundary "
	       "(type c to use %g",sd,mesg,elname,aname[dim],plural[dim],
	       degrees(ang[0]));
	if (dim == 3)
	    screen(" %g",degrees(ang[1]));
	screen("): ");
	(void) Gets(s);
	if (s[0] == '\0')
	    return NO;
	switch (dim)
	{
	case 2:
	    if (tolower(s[0]) != 'c')
	        (void) sscan_float(s,ang);
	    odir[0] = cos(ang[0]);
	    odir[1] = sin(ang[1]);
	    break;
	case 3:
	    if (tolower(s[0]) != 'c')
	    {
	        static const char *fmt = "%lf %lf";
	        if (sscanf(s,fmt,theta,theta+1) != 2)
		    return NO;
	    }
	    odir[0] = cos(ang[0])*sin(ang[1]);
	    odir[1] = sin(ang[0])*sin(ang[1]);
	    odir[2] = cos(ang[1]);
	}
	return YES;
}		/*end add_ray_to_boundary*/


LOCAL	int	read_vector(
	double	*outvec,
	double	*dfltvec,
	int	dim)
{
	char	   s[Gets_BUF_SIZE];
	char	   *st;
	static const char *separators = " \t";
	int        n = 0;
	int	   i;

	if (dim < 1)
	    return 0;
	(void) Gets(s);
	if (s[0] == '\0')
	{
	    if (dfltvec != NULL)
	    {
	    	for (i = 0; i < dim; ++i)
	    	    outvec[i] = dfltvec[i];
	    	return dim;
	    }
	    else
	    	return 0;
	}

	st = strtok(s,separators);
	if (sscan_float(st,outvec) != 1)
	{
	    screen("ERROR in read_vector(), no input uni_array\n");
	    clean_up(ERROR);
	}
	++n;
	for (i = 1; i < dim; ++i)
	    outvec[i] = outvec[0];
	for (i = 1; i < dim; ++i)
	{
	    outvec[i] = outvec[i-1];
	    if ((st = strtok(NULL,separators)) == NULL)
	    	break;
	    (void) sscan_float(st,outvec+i);
	    ++n;
	}
	return n;
}		/*read_vector*/

LOCAL	void	prompt_for_rotation(
	const char *mesg,
	double	   **Q,
	int	   dim)
{
	double	   **T, **S;
	double	   alpha;
	double	   cosa, sina;
	char	   s[Gets_BUF_SIZE];
	int	   i, j, k, l;
	const char *elname = (dim == 3) ? "ellipsoid" : "ellipse";

	screen("Type 'y' to rotate the %s %s: ",mesg,elname);
	(void) Gets(s);
	if (s[0] != 'y' && s[0] != 'Y') return;

	T = NULL;
	if (dim == 3)
	{
	    int	n;
	    double    cost, sint, cosp, sinp;
	    double    magr, r[MAXD];
	    double    phi, theta;

	    screen_print_long_string(
		"The axis of rotation of the rotation can be entered in two "
		"different ways. Either enter the two spherical coordinates "
		"(theta, phi) in degrees of the vector or the three "
		"rectangular coordinates of the uni_array. The default is the "
		"z-axis. Type '?' for an explanation of theta and phi.\n");
	    screen("Enter coordinates: ");
	    (void) Gets(s);
	    if (s[0] == '?')
	    {
	        screen_print_long_string(
	            "The angle phi is the counter-clockwise angle of the "
	    	    "axis vector with respect to the positive z axis, "
	    	    "and the angle theta is the counter-clockwise angle of "
	    	    "the projection of the axis vector in the x-y plane "
	            "with the positive x-axis.\n");
	        screen("Enter coordinates: ");
	        (void) Gets(s);
	    }
	    if (s[0] != '\0')
	    {
	        bi_array(&T,dim,dim,FLOAT);
	        n = sscanf(s,"%lf %lf %lf",r,r+1,r+2);
	        if (n == 2)
	        {
	            theta = radians(r[0]);
		    cost = cos(theta);
		    sint = sin(theta);
		    phi = radians(r[1]);
		    cosp = cos(phi);
		    sinp = sin(phi);
		}
		else if (n == 3)
		{
		    magr = mag_vector(r,dim);
		    for (i = 0; i < dim; ++i)
		        r[i] /= magr;
		    cosp = r[2];
		    sinp = hypot(r[0],r[1]);
		    if (fabs(sinp) < EPSILON)
		    {
		        sinp = 0.0;
		        cosp = 1.0;
		        cost = 1.0;
		        sint = 0.0;
		    }
		    else
		    {
		        cost = r[0]/sinp;
		        sint = r[1]/sinp;
		    }
		}
		else
		{
		    sinp = ERROR_FLOAT;
		    cosp = ERROR_FLOAT;
		    cost = ERROR_FLOAT;
		    sint = ERROR_FLOAT;
		    screen("ERROR in prompt_for_rotation(), "
		           "Invalid uni_array\n");
		    clean_up(ERROR);
		}
		T[0][0] = cost*cosp; T[0][1] = -sint; T[0][2] = cost*sinp;
		T[1][0] = sint*cosp; T[1][1] = cost;  T[1][2] = sint*sinp;
		T[2][0] = -sinp;     T[2][1] = 0.0;   T[2][2] = cosp;
	    }
	}
	bi_array(&S,dim,dim,FLOAT);
	for (i = 0; i < dim; ++i)
	    for (j = 0; j < dim; ++j)
	    	S[i][j] = (i == j) ? 1.0 : 0.0;
	screen("Enter the rotation angle (in degrees) "
	       "about the axis of rotation: ");
	(void) Scanf("%f\n",&alpha);
	alpha = radians(alpha);
	cosa = cos(alpha);	sina = sin(alpha);
	S[0][0] = cosa;	S[0][1] = -sina;
	S[1][0] = sina;	S[1][1] =  cosa;
	if (T != NULL)
	{
	    for (i = 0; i < dim; ++i)
	    for (j = 0; j < dim; ++j)
	    {
	        Q[i][j] = 0.0;
	        for (k = 0; k < dim; ++k)
	        for (l = 0; l < dim; ++l)
	            Q[i][j] += T[i][k]*S[k][l]*T[j][l];
	    }
	}
	else
	{
	    for (i = 0; i < dim; ++i)
	    for (j = 0; j < dim; ++j)
	        Q[i][j] = S[i][j];
	}

	if (T != NULL)
	    free(T);
	free(S);
}		/*end prompt_for_rotation*/

LOCAL void axis_elliptical_state_from_Gas(
	COMP_TYPE	*comp_type,
	Locstate	ans)
{
	_ELLIPTICAL	*ellip;

	ellip = Elliptical(comp_type);
	set_state(ellip->state,TGAS_STATE,ans);
	if (debugging("elrp"))
	{
	    verbose_print_state("ellip->state",ellip->state);
	}
}		/*end axis_elliptical_state_from_Gas*/


LOCAL void Gas_from_axis_elliptical_state(
	Locstate	ans,
	COMP_TYPE	*comp_type)
{
	_ELLIPTICAL	*ellip;

	if (comp_type->type == AMBIENT)
	{
	    set_state(ans,GAS_STATE,Ambient(comp_type));
	    return;
	}
	if (comp_type->type == RANDOM_REGION)
	{
	    set_state(ans,GAS_STATE,Mean(Random_state(comp_type)));
	    return;
	}

	ellip = Elliptical(comp_type);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(ellip->state)->n_comps != 1)
                set_state(ans,TGAS_STATE,ellip->state);
            else
                set_state(ans,GAS_STATE,ellip->state);
            return;
        }
	set_state(ans,GAS_STATE,ellip->state);
}		/*end Gas_from_axis_elliptical_state*/

LOCAL void prompt_for_elliptical_state(
	INIT_PHYSICS	*ip,
	COMPONENT	comp,
	ELLIPSOID	*ellipsoid,
	ELLIPSOID	*inner,
	ELLIPSOID	*outer,
	Gas_param	*params,
	const char	*message,
	INIT_DATA	*init)
{
	ELLIPTICAL_REGION_TYPE ert;
	boolean		       random;
	char		       s[Gets_BUF_SIZE];
	Front		       *front = ip->root->front;
	COMP_TYPE	       *ctype = comp_type(comp);
	_ELLIPTICAL	       *ellpt;
	double		       rho,pr,vr, maxr;
	int		       i, dim = ellipsoid->dim;

	screen("Initialize the elliptical state%s -\n",message);

	random = prompt_for_elliptic_region_type(ip,init,&ert,inner,outer);

	if ((ert == AMBIENT_ELLIPTICAL_REGION) ||
	    (ert == RANDOM_AMBIENT_ELLIPTICAL_REGION))
	{
	    set_up_ambient_elliptic_region(random,ctype,params,
					   message,front,init);
	    return;
	}
	else if(ert == READ_TGAS_STATE_FROM_1D_FILE_REGION)
        {
            set_up_read_state_from_1d_file_region(ctype,params,message,
			    front,init);
            return;
        }
	else if (ert == TABULATED_ELLIPTICAL_REGION)
	{
	    set_tabulated_region_comp_type(ctype,front);
	    ctype->params = params;
	    set_up_read_state_from_file_region(ctype,message,front);
	    return;
	}
	else if (ert == ONE_DIMENSIONAL_ELLIPTICAL_OVERLAY)
	{
	    set_up_elliptical_overlay_region(ctype,params,message,ip,init);
	    return;
	}

	set_elliptical_comp_type(ctype,ip);
	ellpt = Elliptical(ctype);
	ellpt->ellipsoid = ellipsoid;
	for (maxr = 0.0, i = 0; i < dim; ++i)
	    maxr = max(maxr,ellipsoid->rad[i]);
	for (i = 0; i < dim; ++i)
	    ellpt->weight[i] = ellipsoid->rad[i]/maxr;

	if ((ert == ELLIPTICAL_RAREFACTION_REGION) &&
	    (rarefaction_region(inner,outer) == YES))
	{
	    ellpt->rw1d = allocate_RAREFACTION_WAVE_1D(front);
	    screen("\tEnter the gas state at the rarefaction leading "
		   "edge%s\n\t\trho, pr, radial velocity: ",message);
	    (void) Scanf("%f %f %f\n",&rho,&pr,&vr);
	    Dens(ellpt->rw1d->stl) = rho;
	    Press(ellpt->rw1d->stl) = pr;
	    Vel(ellpt->rw1d->stl)[0] = vr;
	    Init_params(ellpt->rw1d->stl,params);
	    set_type_of_state(ellpt->rw1d->stl,TGAS_STATE);
	    reset_gamma(ellpt->rw1d->stl);
	    prompt_for_rarefaction_edge(ellpt->rw1d,TRAILING_EDGE,
					inner,outer,message);
	    verbose_print_state("Initialized elliptical state "
				"for rarefaction region",ellpt->state);
	    print_rarefaction_wave_1d(ellpt->rw1d);
	}
	else
	{
	    screen("\tEnter the gas state%s\n\t\trho, pr, radial velocity: ",
		   message);
	    (void) Scanf("%f %f %f\n",&rho,&pr,&vr);

#if defined(COMBUSTION_CODE)
	    prompt_for_burning(&params,message);
#endif /* defined(COMBUSTION_CODE) */
	    Dens(ellpt->state) = rho;
	    Press(ellpt->state) = pr;
	    RadialVelocity(ellpt) = vr;
	    Init_params(ellpt->state,params);
	    set_type_of_state(ellpt->state,TGAS_STATE);
	    ctype->params = params;
	    reset_gamma(ellpt->state);

	    verbose_print_state("Initialized elliptical state",ellpt->state);
	}

	if ((ert == RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION) ||
	    (ert == RANDOM_RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION))
	{
	    ellpt->stratification_type = prompt_for_stratification("");
	    if (ellpt->stratification_type != CONSTANT)
	    {
	        ellpt->r0 = maxr;
		screen("Enter the reference radius for the stratification "
		       "(dflt = %g): ",ellpt->r0);
		(void) Gets(s);
		if (s[0] != '\0')
		    (void) sscan_float(s,&ellpt->r0);
	    }
	}

	if (random == YES)
	{
	    ellpt->rstate = allocate_random_state_structure(front->interf);
	    alloc_state(front->interf,ellpt->wkstate,front->sizest);
	    alloc_state(front->interf,ellpt->wkstate+1,front->sizest);
	    init_random_state_region(ellpt->state,ellpt->rstate,front);
	}
}		/*end prompt_for_elliptical_state*/

LOCAL	void set_up_elliptical_overlay_region(
	COMP_TYPE    *ctype,
	Gas_param    *params,
	const char   *message,
	INIT_PHYSICS *ip,
	INIT_DATA    *init)
{
	set_1d_overlay_comp_type(ctype);
	ctype->params = params;
	prompt_for_1d_overlay(init,ctype,ip,message,RADIAL_OVERLAY);
}		/*end set_up_elliptical_overlay_region*/

LOCAL	boolean	prompt_for_elliptic_region_type(
	const INIT_PHYSICS     *ip,
	const INIT_DATA        *init,
	ELLIPTICAL_REGION_TYPE *ert,
	ELLIPSOID              *inner,
	ELLIPSOID              *outer)
{
	GRAVITY *grav_data = gravity_data(init);
	boolean random = NO;
	boolean allow_radial_stratified;
	char    s[Gets_BUF_SIZE];

	*ert = ELLIPTICAL_VELOCITY_REGION;
	allow_radial_stratified =
	    ((problem_type(ip) == RADIAL_RAYLEIGH_TAYLOR) &&
	     (grav_data->type == RADIAL_GRAVITY)) ? YES : NO;

	screen("The choices for flow initialization types are\n");
	screen("\tAmbient region (constant initial conditions) (A%s)\n",
	       (*ert == AMBIENT_ELLIPTICAL_REGION) ? ", default" : "");
	screen("\tTabulated region (initial conditions read from a file) "
	       "(TR%s)\n",
	       (*ert == TABULATED_ELLIPTICAL_REGION) ? ", default" : "");
	screen("\tRandom perturbation of an ambient region (RA%s)\n",
	       (*ert == RANDOM_AMBIENT_ELLIPTICAL_REGION) ? ", default" : "");
	screen("\tElliptical region "
	       "(constant thermodynamics, radial velocity) (E%s)\n",
	       (*ert == ELLIPTICAL_VELOCITY_REGION) ? ", default" : "");
	if (rarefaction_region(inner,outer) == YES)
	    screen("\tElliptical rarefaction region "
	           "(1d rarefaction, radial velocity) (ER%s)\n",
	           (*ert == ELLIPTICAL_RAREFACTION_REGION) ? ", default" : "");
	screen("\tRandom perturbation of elliptical region (RE%s)\n",
	       (*ert == RANDOM_ELLIPTICAL_VELOCITY_REGION) ? ", default" : "");
	screen("\tRead states from 1d file region "
               "(density, velicity, pressure (RS1%s)\n",
               (*ert == READ_TGAS_STATE_FROM_1D_FILE_REGION) ? ", default" : "");
	if (allow_radial_stratified == YES)
	{
	    screen("\tRadial gravity stratified region (RGSR%s)\n",
		   (*ert == RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION) ?
		       ", default" : "");
	    screen("\tRandom perturbation of radial gravity "
		   "stratified region (RRGSR%s)\n",
		   (*ert==RANDOM_RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION) ?
		   ", default" : "");
	}
	screen("\tOverlay of one dimensional flow (O%s)\n",
	       (*ert == ONE_DIMENSIONAL_ELLIPTICAL_OVERLAY) ? ", default" : "");
	screen("Enter choice: ");
	(void) Gets(s);
	if (strcasecmp(s,"A") == 0)
	{
	    *ert = AMBIENT_ELLIPTICAL_REGION;
	    random = NO;
	}
	else if (strcasecmp(s,"TR") == 0)
	{
	    *ert = TABULATED_ELLIPTICAL_REGION;
	    random = NO;
	}
	else if (strcasecmp(s,"RA") == 0)
	{
	    *ert = RANDOM_AMBIENT_ELLIPTICAL_REGION;
	    random = YES;
	}
	else if ((rarefaction_region(inner,outer) == YES) &&
		 (strncasecmp(s,"ER",2) == 0))
	{
	    *ert = ELLIPTICAL_RAREFACTION_REGION;
	    outer_ellipsoid(inner) = outer;
	    inner_ellipsoid(outer) = inner;
	    random = NO;
	}
	else if (strcasecmp(s,"E") == 0)
	{
	    *ert = ELLIPTICAL_VELOCITY_REGION;
	    random = NO;
	}
	else if (strcasecmp(s,"RE") == 0)
	{
	    *ert = RANDOM_ELLIPTICAL_VELOCITY_REGION;
	    random = YES;
	}
	else if (strcasecmp(s,"RS1") == 0)
        {
            *ert = READ_TGAS_STATE_FROM_1D_FILE_REGION;
            random = NO;
        }
	else if (strcasecmp(s,"O") == 0)
	{
	    *ert = ONE_DIMENSIONAL_ELLIPTICAL_OVERLAY;
	    random = NO;
	}
	else if ((allow_radial_stratified == YES) &&
	         (strcasecmp(s,"RGSR") == 0))
	{
	    *ert = RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION;
	    random = NO;
	}
	else if ((allow_radial_stratified == YES) &&
	         (strcasecmp(s,"RRGSR") == 0))
	{
	    *ert = RANDOM_RADIAL_STRATIFIED_ELLIPTICAL_VELOCITY_REGION;
	    random = YES;
	}
	return random;
}		/*end prompt_for_elliptic_region_type*/

LOCAL	boolean	rarefaction_region(
	ELLIPSOID *inner,
	ELLIPSOID *outer)
{
	int wt_in, wt_out;

	if ((inner == NULL) || (outer == NULL))
	    return NO;

	wt_in = inner->wv_type;
	wt_out = outer->wv_type;
	if (!(is_rarefaction_wave(wt_in) && is_rarefaction_wave(wt_out)))
	    return NO;

	if (is_forward_wave(wt_in) && !is_forward_wave(wt_out))
	    return NO;

	if (is_backward_wave(wt_in) && !is_backward_wave(wt_out))
	    return NO;

	if (is_rarefaction_leading_edge(wt_in) &&
	    !is_rarefaction_trailing_edge(wt_out))
	    return NO;

	if (is_rarefaction_trailing_edge(wt_in) &&
	    !is_rarefaction_leading_edge(wt_out))
	    return NO;
	return YES;
}		/*end rarefaction_region*/

LOCAL	void	set_up_ambient_elliptic_region(
	boolean	   random,
	COMP_TYPE  *ctype,
	Gas_param  *params,
	const char *message,
	Front	   *front,
	INIT_DATA  *init)
{
	char		s[Gets_BUF_SIZE];
	Locstate	mean, st;

	prompt_for_ambient_state(ctype,params,message,front,init);
	st = Ambient(ctype);
	if (random == NO)
	{
	    screen("Is the flow for this region constant (dflt = no): ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	       (void)SetConstantFlowRegion(ctype->comp,st,front->interf);
	}
	else
	{
	    alloc_state(front->interf,&mean,front->sizest);
	    copy_state(mean,st);
	    (*ctype->free_comp_type_extra)(ctype);
	    st = NULL;
	    set_random_region_comp_type(ctype,front);
	    init_random_state_region(mean,Random_state(ctype),front);
	    free(mean);
	}
}		/*end set_up_ambient_region*/

/*
*			get_state_elliptical():
*
*	Sets state with given density and pressure and radial elliptical
*	velocity.  The velocity is given by
*
*                                    T
*	v = vr * Q * diag(weight) * Q * dr
*
*	where dr is the displacement vector of coords with respect
*	to the center of the ellipse,  Q is the rotation matrix of
*	the ellipse, weight is vector of weights between 0 and 1,
*	and diag(weight) is the diagonal matrix for the	weight vector weight.
*	The weights are determined by the ratio of the axis lengths of the
*	ellipse to the length of the maximum axis in the rotated coordinates.
*/

/*ARGSUSED*/
LOCAL	void get_state_elliptical(
	double		*coords,
	Locstate	state,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int             stype)
{
	_ELLIPTICAL *ellip = Elliptical(ct);
	GRAVITY     *grav_data = gravity_data(init);
	double 	    length;
	double	    vr;
	double 	    *cen = ellip->ellipsoid->cen;
	double 	    **Q = ellip->ellipsoid->Q;
	double 	    *weight = ellip->weight;
	double 	    dr[MAXD], drp[MAXD], v[MAXD], vp[MAXD];
	int 	    i, j, dim = ellip->ellipsoid->dim;
	
	debug_print("init_states","Entered g_state_elliptical()\n");
	if (ct->type != ELLIPTICAL)
	{
	    screen("ERROR in get_state_elliptical(), "
	           "inconsistent comp_type->type.\n");
	    clean_up(ERROR);
	}
	if ((ellip->stratification_type != CONSTANT) &&
	    (grav_data->type == RADIAL_GRAVITY))
	{
	    double r, dr;
	    int   dim = intfc->dim;
	    r = distance_between_positions(coords,grav_data->center,dim);
	    dr = r - ellip->r0;
	    get_state_in_stratified_region(ellip->stratification_type,state,
					   dr,grav_data->G,ellip->state);
	    vr = RadialVelocity(ellip);
	}
	else if (ellip->rw1d != NULL)
	{
	    _RAREFACTION_WAVE_1D *rw1d = ellip->rw1d;
	    Locstate stl, stt;

	    stl = rw1d->stl;
	    stt = rw1d->stt;
	    if (rarefaction_edge_at_coords(coords,hs,LEADING_EDGE) == YES)
	        set_state(state,TGAS_STATE,stl);
	    else if (rarefaction_edge_at_coords(coords,hs,TRAILING_EDGE) == YES)
	        set_state(state,TGAS_STATE,stt);
	    else
	    {
	        double r = distance_between_positions(coords,cen,dim);
	        double rmin, rmax, zbar, tbar;
	        if (perturbed_elliptical_rarefaction_region(rw1d) == YES)
	        {
		    double rl, rt;
		    double theta[3], dir[3];
		    double crdsl[3], crdst[3];
		    for (i = 0; i < dim; ++i)
			dir[i] = (coords[i] - cen[i])/r;
		    theta[0] = angle(dir[0],dir[1]);
		    if (dim == 3)
			theta[1] = angle(dir[2],hypot(dir[0],dir[1]));
		    coords_on_pert_ellipsoid(theta,crdsl,rw1d->el_lead);
		    rl = distance_between_positions(crdsl,cen,dim);
		    coords_on_pert_ellipsoid(theta,crdst,rw1d->el_trail);
		    rt = distance_between_positions(crdst,cen,dim);
		    if (rl < rt)
		    {
			rmin = rl;
			rmax = rt;
		    }
		    else
		    {
			rmin = rt;
			rmax = rl;
		    }
		    tbar = rw1d->tbar;
	            zbar = 0.5*(rl + rt - (rw1d->spl + rw1d->spt)*tbar);
	        }
		else
		{
	            rmin = rw1d->zmin;
	            rmax = rw1d->zmax;
		    zbar = rw1d->zbar;
		    tbar = rw1d->tbar;
		}
	        if (r <= rmin)
		    set_state(state,TGAS_STATE,rw1d->stmin);
	        else if (r >= rmax)
		    set_state(state,TGAS_STATE,rw1d->stmax);
	        else
	        {
	            double speed, vl, spdnew;
	            speed = (r - zbar)/tbar;
	            vl = vel(0,stl);
	            (void) oned_state_in_rarefaction_fan(speed,vl,stl,stt,state,
							 TGAS_STATE,&spdnew,
						         rw1d->l_or_r);
	        }
	    }
	    vr = vel(0,state);
	}
	else
	{
	    set_state(state,TGAS_STATE,ellip->state);
	    vr = RadialVelocity(ellip);
	}
	for (i = 0; i < dim; ++i)
	    dr[i] = coords[i] - cen[i];
	length = mag_vector(dr,dim);
	if (length > EPSILON) 
	{
	    for (i = 0; i < dim; ++i)
	    	dr[i] /= length;
	    for (i = 0; i < dim; ++i)
	    {
	    	for (drp[i] = 0.0, j = 0; j < dim; ++j)
		    drp[i] += Q[j][i]*dr[j];
	    	vp[i] = vr * weight[i] * drp[i];
	    }
	    for (i = 0; i < dim; ++i)
	    {
	    	for (v[i] = 0.0, j = 0; j < dim; ++j)
		    v[i] += Q[i][j]*vp[i];
	    	Vel(state)[i] = v[i];
	    }
	}
	else 
	{
	    for (i = 0; i < dim; ++i)
	    	Vel(state)[i] = 0.0;
	}
	if (ellip->rstate != NULL)
	{
	    Locstate	states[3];
	    static const double alpha[3] = {1, 1, -1};

	    states[0] = ellip->wkstate[0];
	    states[1] = ellip->wkstate[1];
	    states[2] = Mean(ellip->rstate);
	    set_state(states[0],state_type(Mean(ellip->rstate)),state);
	    random_region_state(states[1],coords,ellip->rstate);
	    g_linear_combination_of_states(alpha,states,3,state);
	}
        /* This is for the oned sod-tube test problem */
        /* init partial density */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(debugging("implosion_multi"))
            {
                int    num_comps;
                double   *rho;
                if((num_comps = Params(state)->n_comps) != 1)
                {
                    rho = pdens(state);
                    if(coords[0] <= 0.05)
                        rho[0] = 0.8*Dens(state);
                    else if (coords[0] > 0.05 && coords[0] <= 0.1)
                        rho[0] = 0.3*Dens(state);
                    else
                        rho[0] = 0.1*Dens(state);
                    rho[1] = Dens(state) - rho[0];
                }
            }
        }
	set_state(state,stype,state);
}		/*end get_state_elliptical*/

LOCAL	boolean perturbed_elliptical_rarefaction_region(
	_RAREFACTION_WAVE_1D *rw1d)
{
	if (rw1d->el_lead == NULL || rw1d->el_trail == NULL)
	    return NO;
	return ((perturbed_ellipsoid(rw1d->el_lead) == YES) ||
	        (perturbed_ellipsoid(rw1d->el_trail) == YES)) ? YES : NO;
}		/*end perturbed_elliptical_rarefaction_region*/

LOCAL	boolean perturbed_ellipsoid(
	ELLIPSOID *el)
{
	return ((el->fpoly != NULL) || (el->lpoly != NULL)) ? YES : NO;
}		/*end perturbed_ellipsoid*/

LOCAL   void init_input_modes(
	int		offset,
	int		min_n,
	int		max_n,
	int		nmodes,
	FOURIER_POLY	*fpoly,
	double		*L,
	double		*U)
{
	double		**wv_num = fpoly->nu+offset, *A = fpoly->A + offset;
	double		*phase = fpoly->phase + offset;
        int		dim = fpoly->dim;
        int		i, j, n, m[2];
        char		s[Gets_BUF_SIZE];
	static const char *fmt = "%lf %lf";
	static const char *direction[2] = {"in the azmuthal direction",
					   "in the polar direction"};

	(void) printf("\tnumber of modes::%d\n",nmodes);
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= max_n - min_n; ++i)
	    {
	        screen("Enter amplitude and phase (deg.) for mode %d: ", i);
	        (void) Gets(s);
	        (void) sscanf(s,fmt, &A[i], &phase[i]);
	        phase[i] = radians(phase[i]);
		wv_num[i][0] = prompt_for_wave_number(min_n+i,i,
						      L[0],U[0],"");
	        phase[i] -= L[0]*wv_num[i][0];
	    }
            break;
        case 3:
            i = 0;
            for (n = min_n; n <= max_n; ++n)
            {
                for (m[0] = 0; m[0] <= n; ++m[0])
                {
		    m[1] = n - m[0];
		    screen("Enter amplitude and phase (deg.) for mode %d: ", i);
		    (void) Gets(s);
	            (void) sscanf(s,fmt, &A[i], &phase[i]);
		    phase[i] = radians(phase[i]);

		    for (j = 0; j < 2; ++j)
                    {
                        wv_num[i][j] = prompt_for_wave_number(m[j],i,
						              L[j],U[j],
							      direction[j]);
                        phase[i] += L[j]*wv_num[i][j];
                    }
                    ++i;
                }
            }
            break;
        default:
	    screen("ERROR in init_input_modes(), invalid dim = %d\n",dim);
            clean_up(ERROR);
            break;
        }
}		/*end init_input_modes*/

LOCAL	double	prompt_for_wave_number(
	double nu,
	int   i,
	double L,
	double U,
	const char *direction)
{
	double  wv_num;
	double  period;
	size_t j, len;
        char   *c, s[Gets_BUF_SIZE];

	wv_num = 2.0*PI*nu/(U-L);
	period = (U-L)/nu;
	screen("The frequency and wave number for mode %d%s ",i,direction);
	screen("over the angular interval\n\t"
	       "%g PI -> %g PI (%g -> %g degrees), can be entered in\n\t"
	       "three different ways.\n\t",L/PI,U/PI,degrees(L),degrees(U));
	screen("1. Enter the frequency by simply typing the value\n\t"
	       "   or the string \"frequency =\" followed by the desired "
	           "frequency\n\t"
	       "2. Enter the wave number by typing \"wave number =\" "
                   "followed by the\n\t"
               "   desired wave number.\n\t"
	       "3. Enter the angular period (in degrees) by typing "
                   "\"period =\"\n\t"
	       "   followed by the desired angular period.\n");
	screen("Enter the desired data (dflt frequency   = %g or\n"
	       "                             wave number = %g or\n"
	       "                             period      = %g): ",
	       nu,wv_num,degrees(period));
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((c = strrchr(s,'=')) != NULL)
	    {
		++c;
		if ((strchr(s,'f') != NULL) || (strchr(s,'F') != NULL))
		{
	            (void) sscan_float(c,&nu);
	            wv_num = 2.0*PI*nu/(U-L);
		    period = (U-L)/nu;
		}
		else if ((strchr(s,'w') != NULL) || (strchr(s,'W') != NULL))
	        {
	            (void) sscan_float(c,&wv_num);
		    nu = wv_num*(U-L)/(2.0*PI);
		    period = (U-L)/nu;
	        }
		else if ((strchr(s,'p') != NULL) || (strchr(s,'P') != NULL))
		{
	            (void) sscan_float(c,&period);
		    period = radians(period);
		    nu = (U-L)/period;
	            wv_num = 2.0*PI*nu/(U-L);
		}
	    }
	    else if (sscan_float(s,&nu) == 1)
	    {
	        wv_num = 2.0*PI*nu/(U-L);
		period = (U-L)/nu;
	    }
	    else
	    {
	        screen("ERROR in prompt_for_wave_number(), "
		       "invalid input %s\n",s);
		clean_up(ERROR);
	    }
	}
	(void) sprintf(s,"\tfor mode %d%s",i,direction);
	(void) printf("%s::frequency = %g\n",s,nu);
	len = strlen(s);
	for (j = 1; j < len; ++j)
	    s[j] = ' ';
	(void) printf("%s::wave number = %g\n",s,wv_num);
	(void) printf("%s::period = %g\n",s,degrees(period));
	return wv_num;
}		/*end prompt_for_wave_number*/

                                 
EXPORT	ELLIPSOID  *allocate_ellipsoid(
	ELLIPSOID *template_ellipsoid,
	int	  nmodes)
{
	ELLIPSOID	*ellipsoid;
	size_t		size, naE = 0;
	int		i, j, dim = template_ellipsoid->dim;

	if (nmodes > 0)
	{
	    size_t	naFP, naA, naP, naNUP, naNUS;
	    naE = num_aligns(sizeof(ELLIPSOID));
	    naFP = num_aligns(sizeof(FOURIER_POLY));
	    naA = num_aligns(nmodes*FLOAT);
	    naP = num_aligns(nmodes*FLOAT);
	    naNUP = num_aligns(nmodes*sizeof(double*));
	    naNUS = num_aligns(nmodes*(dim-1)*FLOAT);
	    size = sizeof(ALIGN)*(naE + naFP + naA + naP + naNUP + naNUS);
	}
	else
	    size = sizeof(ELLIPSOID);

	scalar(&ellipsoid,size);

	set_default_ellipsoid_structure(ellipsoid,dim);

	for (i = 0; i < 3; ++i)
	{
	    ellipsoid->cen[i] = template_ellipsoid->cen[i];
	    ellipsoid->rad[i] = template_ellipsoid->rad[i];
	    ellipsoid->vr[i] = template_ellipsoid->vr[i];
	    for (j = 0; j < 3; ++j)
	    	ellipsoid->Q[i][j] = template_ellipsoid->Q[i][j];
	}
	for (j = 0; j < 2; ++j)
	{
	    ellipsoid->ThetaS[j] = template_ellipsoid->ThetaS[j];
	    ellipsoid->ThetaE[j] = template_ellipsoid->ThetaE[j];
	}
	ellipsoid->closed = template_ellipsoid->closed;
	ellipsoid->nor_orient = template_ellipsoid->nor_orient;
	ellipsoid->compin = template_ellipsoid->compin;
	ellipsoid->compout = template_ellipsoid->compout;
	ellipsoid->fpoly = template_ellipsoid->fpoly;
	ellipsoid->hs = template_ellipsoid->hs;
	ellipsoid->surf_tension = template_ellipsoid->surf_tension;
	ellipsoid->wv_type = template_ellipsoid->wv_type;
	ellipsoid->dim = template_ellipsoid->dim;
	ellipsoid->untracked = template_ellipsoid->untracked;
	ellipsoid->reset_position = template_ellipsoid->reset_position;
	ellipsoid->layer_index = 0;
	ellipsoid->scale = template_ellipsoid->scale;
	ellipsoid->_make_ellipsoid = template_ellipsoid->_make_ellipsoid;
	inner_ellipsoid(ellipsoid) = outer_ellipsoid(ellipsoid) = NULL;
	for (i = 0; i < 4; ++i)
	{
	    ellipsoid->btype[i] = template_ellipsoid->btype[i];
	    ellipsoid->bcomp[i] = template_ellipsoid->bcomp[i];
	    ellipsoid->rbdry[i] = template_ellipsoid->rbdry[i];
	    ellipsoid->obtype[i] = template_ellipsoid->obtype[i];
	    ellipsoid->obcomp[i] = template_ellipsoid->obcomp[i];
	}

	if (nmodes > 0)
	{
	    ALIGN *fpstore = ((ALIGN*)ellipsoid) + naE;
	    ellipsoid->fpoly = allocate_fourier_poly(nmodes,dim,fpstore);
	}

	return ellipsoid;
}		/*end allocate_ellipsoid*/


EXPORT	void	set_default_ellipsoid_structure(
	ELLIPSOID	*ellipsoid,
	int		dim)
{
	int		i, j;

	zero_scalar(ellipsoid,sizeof(ELLIPSOID));
	ellipsoid->Q = ellipsoid->Qrows;
	for (i = 0; i < 3; ++i)
	    ellipsoid->Q[i] = &ellipsoid->Qstore[i][0];

	for (i = 0; i < 3; ++i)
	{
	    ellipsoid->cen[i] = 0.0;
	    ellipsoid->rad[i] = 1.0;
	    ellipsoid->vr[i] = -HUGE_VAL;
	    for (j = 0; j < 3; ++j)
	    	ellipsoid->Q[i][j] = (i == j) ? 1.0 : 0.0;
	}
	for (i = 0; i < 4; ++i)
	{
	    ellipsoid->btype[i] = UNKNOWN_BOUNDARY_TYPE;
	    ellipsoid->bcomp[i] = NO_COMP;
	    ellipsoid->obtype[i] = UNKNOWN_BOUNDARY_TYPE;
	    ellipsoid->rbdry[i] = NO;
	    ellipsoid->obcomp[i] = NO_COMP;
	}
	ellipsoid->scale = 1.0;
	ellipsoid->dim = dim;
	ellipsoid->n2o_enforced = NO;
	ellipsoid->compin = ellipsoid->compout = NO_COMP;
	ellipsoid->reset_position = NO;
	switch (dim)
	{
	case 2:
	    ellipsoid->_make_ellipsoid = g_make_ellipse;
	    ellipsoid->_make_ellip_region_boundaries =
		g_make_ellip_region_boundaries2d;
	    break;
	case 3:
	    ellipsoid->_make_ellipsoid = g_make_ellipsoid;
	    ellipsoid->_make_ellip_region_boundaries =
		g_make_ellip_region_boundaries3d;
	    break;
	default:
	    screen("ERROR in set_default_ellipsoid_structure(), "
		   "dim = %d not supported\n",dim);
	    clean_up(ERROR);
	}
}		/*end set_default_ellipsoid_structure*/

EXPORT	void	set_elliptical_comp_type(
	COMP_TYPE	*comp_type,
	INIT_PHYSICS	*ip)
{
	_ELLIPTICAL	*ellip;

	if (comp_type->type == ELLIPTICAL) /*ALREADY SET*/
	    return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = ELLIPTICAL;
	scalar(&ellip,sizeof(_ELLIPTICAL));
	ellip->stratification_type = CONSTANT;
	ellip->r0 = -HUGE_VAL;
	ellip->rw1d = NULL;
	alloc_state(ip->root->front->interf,&ellip->state,
		    ip->root->front->sizest);
	comp_type->extra = (POINTER) ellip;
	comp_type->_get_state = get_state_elliptical;
	comp_type->free_comp_type_extra = free_elliptical_comp_type;
}		/*end set_elliptical_comp_type*/

LOCAL	void	free_elliptical_comp_type(
	COMP_TYPE	*comp_type)
{
	_ELLIPTICAL	*ellip;
	if (comp_type->type != ELLIPTICAL)
	    return;
	ellip = Elliptical(comp_type);

	if (ellip == NULL)
	    return;
	if (ellip->state != NULL)
	    free(ellip->state);
	if (ellip->wkstate[0] != NULL)
	    free(ellip->wkstate[0]);
	if (ellip->wkstate[1] != NULL)
	    free(ellip->wkstate[1]);
	if (ellip->rstate != NULL)
	    free_random_state_structure(ellip->rstate);
	if (ellip->rw1d != NULL)
	{
	    free(ellip->rw1d->stl);
	    free(ellip->rw1d->stt);
	    free(ellip->rw1d);
	}
	free(ellip);
	comp_type->extra = NULL;
}		/*end free_elliptical_comp_type*/

EXPORT	void	print_ellipsoid(
	ELLIPSOID	*ellip,
	INTERFACE	*intfc)
{
	int i, j;
	int dim = ellip->dim;

	(void) printf("ELLIPSOID structure %p\n",ellip);
	(void) printf("dim = %d\n",ellip->dim);
	print_general_vector("center = ",ellip->cen,dim,"\n");
	print_general_vector("radii = ",ellip->rad,dim,"\n");
	if (ellip->Q == NULL)
	    (void) printf("Rotation matrix = NULL\n");
	else
	{
	    (void) printf("Q = ");
	    for (i = 0; i < dim; ++i)
	    {
		for (j = 0; j < dim; ++j)
		    (void) printf("\t%-14g ",ellip->Q[i][j]);
		(void) printf("\n");
	    }
	}
	(void) printf("ThetaS[0] = %g (%g degrees)\n",
		      ellip->ThetaS[0],degrees(ellip->ThetaS[0]));
	(void) printf("ThetaS[1] = %g (%g degrees)\n",
		      ellip->ThetaS[1],degrees(ellip->ThetaS[1]));
	(void) printf("ThetaE[0] = %g (%g degrees)\n",
		      ellip->ThetaE[0],degrees(ellip->ThetaE[0]));
	(void) printf("ThetaE[1] = %g (%g degrees)\n",
		      ellip->ThetaE[1],degrees(ellip->ThetaE[1]));
	(void) printf("closed = %s\n",y_or_n(ellip->closed));
	print_orientation("nor_orient = ",ellip->nor_orient,"\n");
	(void) printf("compin = %d, compout = %d\n",
		      ellip->compin,ellip->compout);
	if (ellip->fpoly == NULL)
	    (void) printf("ellip->fpoly = NULL\n");
	else
	{
	    (void) printf("ellip->fpoly = %p\n",ellip->fpoly);
	    (void) printf("\tnum_modes = %d\n",ellip->fpoly->num_modes);
	    (void) printf("\tdim = %d\n",ellip->fpoly->dim);
	    for (i = 0; i < ellip->fpoly->num_modes; ++i)
	    {
		(void) printf("\tnu[%d] = ",i);
		print_general_vector("",ellip->fpoly->nu[i],
				     ellip->fpoly->dim,"\n");
	    }
	    (void) printf("\tz0 = %g\n",ellip->fpoly->z0);
	    for (i = 0; i < ellip->fpoly->num_modes; ++i)
		(void) printf("\tA[%d] = %g\n",i,ellip->fpoly->A[i]);
	    for (i = 0; i < ellip->fpoly->num_modes; ++i)
		(void) printf("\tphase[%d] = %g (%g degrees)\n",
			      i,ellip->fpoly->phase[i],
			      degrees(ellip->fpoly->phase[i]));
	}
	if (ellip->hs == NULL)
	    (void) printf("hs = NULL\n");
	else
	{
	    (void) printf("hs =\n");
	    print_hypersurface(ellip->hs);
	}
	(void) printf("surf_tension = %g\n",ellip->surf_tension);
	print_wave_type("wv_type = ",ellip->wv_type,"\n",intfc);
	(void) printf("untracked = %s\n",y_or_n(ellip->untracked));
	(void) printf("layer_index = %d\n",ellip->layer_index);
	(void) printf("End ELLIPSOID structure %p\n",ellip);
}		/*end print_ellipsoid*/

LOCAL void set_up_read_state_from_1d_file_region(
	COMP_TYPE	*ctype,
      	Gas_param	*params,
      	const char	*message,
      	Front		*front,
      	INIT_DATA	*init)
{       
	ONED_TABLE     *table;
        FILE            *fd, *fp, *fvr, *fpt;
        char            den[Gets_BUF_SIZE];
        char            pre[Gets_BUF_SIZE];
        char            rvel[Gets_BUF_SIZE]; 
        char            pt[Gets_BUF_SIZE];
 	char		s[Gets_BUF_SIZE];
        char            *data_list[4];   
        int             rmax, ir;
	char            tname[Gets_BUF_SIZE],*gasdir,oname[1024];


        screen("Initialize the elliptical state%s -\n",message);

        screen("\tEnter the integer for the number of data: ");
        (void) Scanf("%d\n", &rmax);

        screen("\nEnter the file name of density: ");
        (void) Gets(den);
        
        screen("\nEnter the file name of pressure: ");
        (void) Gets(pre);
        
        screen("\nEnter the file name of velocity: ");
        (void) Gets(rvel);
        
        screen("\nEnter the file name of position: ");
        (void) Gets(pt);

	screen("\nDo the above files exist (dflt = y): ");
	(void) Gets(s);

	strcpy(oname,output_filename(init));
	gasdir = get_dirname(dirname(oname));
	sprintf(tname,"%s/%s",gasdir,den);
        fd = fopen(tname, "r");
	sprintf(tname,"%s/%s",gasdir,pre);
        fp = fopen(tname, "r");
	sprintf(tname,"%s/%s",gasdir,rvel);
        fvr = fopen(tname, "r");
	sprintf(tname,"%s/%s",gasdir,pt);
        fpt = fopen(tname, "r");

        if (ctype->free_comp_type_extra != NULL)
	    (*ctype->free_comp_type_extra)(ctype); 
       
        ctype->type = ONE_DIMENSIONAL_TABLE;
        scalar(&table,sizeof(ONED_TABLE));
    
         table->nx = rmax;
  
        uni_array(&(table->d), table->nx, FLOAT);
        uni_array(&(table->p), table->nx, FLOAT);
        uni_array(&(table->v), table->nx, FLOAT);
        uni_array(&(table->pt), table->nx, FLOAT);

        ctype->extra = (POINTER) table;
        ctype->_get_state = get_state_1dtable;
        ctype->free_comp_type_extra = free_1dtable_comp_type;
#if defined(COMBUSTION_CODE)
	prompt_for_burning(&params,message);
#endif /* defined(COMBUSTION_CODE) */	
	ctype->params = params;

        for(ir = 0; ir < rmax; ir++)
            { 
              if(fscan_float(fd, &(table->d)[ir]) != 1)
                { 
                  screen("ERROR in set_up_read_state_from_1d_file_region(), "
                         "density file reading failed\n");
                  clean_up(ERROR);
                }
            
              if(fscan_float(fp, &(table->p)[ir]) != 1)
                { 
                  screen("ERROR in set_up_read_state_from_1d_file_region(), "
                         "pressure file reading failed\n");
                  clean_up(ERROR);
                }

              if(fscan_float(fvr, &(table->v)[ir]) != 1)
                { 
                  screen("ERROR in set_up_read_state_from_1d_file_region(), "
                         " velocity file reading failed\n");
                  clean_up(ERROR);
                }
              
              if(fscan_float(fpt, &(table->pt)[ir]) != 1)
                { 
                  screen("ERROR in set_up_read_state_from_1d_file_region(), "
                         " position file reading failed\n");
                  clean_up(ERROR);
                }
            }           

   
	fclose(fd);
        fclose(fp);
        fclose(fvr);
        fclose(fpt);

}	/*end set_up_read_state_from_1d_file_region*/

LOCAL void init_SN_layers(
	Front		*front,
	double		dr,
	int		imin,
	int		imax,
	double		*pos,
	double		*den,
	double		*pre,
	double		*grv,
	Locstate	st_start)
{
	double	G = astrophys_grav_constant();
	double	g0,g1;
	double	r0,r1;
	Locstate st0,st1;
	int	i,i_break = 0;
	boolean	break_here = NO;
	static	int indx = 0;
	char	dname[100],pname[100],gname[100];
	FILE	*dfile,*gfile,*pfile;

	debug_print("init_SN","Entered init_SN_layers()\n");

	alloc_state(front->interf,&st0,front->sizest);
	alloc_state(front->interf,&st1,front->sizest);

	/* Initialize center state */

	Params(st0) = Params(st1) = Params(st_start);
	Dens(st0) = Dens(st1) = Dens(st_start);
	Press(st0) = Press(st1) = pressure(st_start);
	Vel(st0)[0] = Vel(st0)[1] = 0.0;
	Vel(st1)[0] = Vel(st1)[1] = 0.0;
	set_type_of_state(st0,TGAS_STATE);
	set_type_of_state(st1,TGAS_STATE);
	reset_gamma(st0);
	reset_gamma(st1);

	g0 = g1 = grv[imin];
	r0 = r1 = pos[imin];

	if (debugging("init_SN"))
	{
	    sprintf(dname,"SN_density-%d",indx);
	    sprintf(pname,"SN_pressure-%d",indx);
	    sprintf(gname,"SN_gravity-%d",indx);
	    dfile = fopen(dname,"w");
	    gfile = fopen(gname,"w");
	    pfile = fopen(pname,"w");
	    fprintf(dfile,"\"density\"\n");
	    fprintf(gfile,"\"gravity\"\n");
	    fprintf(pfile,"\"pressure\"\n");
	    fprintf(dfile,"%f %g\n",pos[imin],den[imin]);
	    fprintf(gfile,"%f %g\n",pos[imin],grv[imin]);
	    fprintf(pfile,"%f %g\n",pos[imin],pre[imin]);
	}

	for (i = imin+1; i <= imax; ++i)
	{
	    r1 = r0 + dr;
	    Dens(st1) = Dens(st0)*(1.0 - dr*g0/sound_speed_squared(st0)); 
	    Press(st1) = Press(st0) - dr*Dens(st0)*g0;
	    g1 = (g0*sqr(r0) + dr*4.0*PI*G*Dens(st0)*sqr(r0))/sqr(r1);
	    reset_gamma(st1);

	    Dens(st0) = Dens(st0) - 0.5*dr*
	    		(Dens(st0)*g0/sound_speed_squared(st0) + 
			 Dens(st1)*g1/sound_speed_squared(st1));
	    if (Dens(st0) < 1e-12)
	    {
		Dens(st0) = 1e-12;
		i_break = i;
		break_here = YES;
	    }
	    Press(st0) = Press(st0) - 0.5*dr*(Dens(st0)*g0 + Dens(st1)*g1);
	    g0 = (g0*sqr(r0) + dr*2.0*PI*G*(Dens(st0)*sqr(r0) + 
				            Dens(st1)*sqr(r1)))/sqr(r1);
	    reset_gamma(st0);

	    pos[i] = r1;
	    den[i] = max(Dens(st0),Vacuum_dens(st0));
	    pre[i] = max(Press(st0),Min_pressure(st0));
	    grv[i] = g0;
	    r0 = r1;
	    if (debugging("init_SN"))
	    {
	    	fprintf(dfile,"%f %g\n",pos[i],den[i]);
	    	fprintf(gfile,"%f %g\n",pos[i],grv[i]);
	    	fprintf(pfile,"%f %g\n",pos[i],pre[i]);
	    }
	    if (break_here)
	    	break;
	}
	if (break_here)
	{
	    for (i = i_break+1; i <= imax; ++i)
	    {
	        r1 = r0 + dr;
	    	pos[i] = r1;
	    	den[i] = den[i-1];
	    	pre[i] = pre[i-1];
	    	grv[i] = grv[i-1];
	    	r0 = r1;
	    	if (debugging("init_SN"))
	    	{
	    	    fprintf(dfile,"%f %g\n",pos[i],den[i]);
	    	    fprintf(gfile,"%f %g\n",pos[i],grv[i]);
	    	    fprintf(pfile,"%f %g\n",pos[i],pre[i]);
	    	}
	    }
	}
	indx++;
	if (debugging("init_SN"))
	{
	    fclose(dfile);
	    fclose(gfile);
	    fclose(pfile);
	}
	debug_print("init_SN","Left init_SN_layers()\n");
} /* end init_SN_layers */

LOCAL	void prompt_for_rarefaction_edge(
	_RAREFACTION_WAVE_1D  *rw1d,
	RAREFACTION_EDGE_TYPE type,
	ELLIPSOID             *inner,
	ELLIPSOID             *outer,
	const char            *message)
{
	char     s[Gets_BUF_SIZE];
	double    sgn, pr, rl, rt, vl, vt, spl, spt;
	Locstate stl = rw1d->stl;
	Locstate stt = rw1d->stt;

	screen("Enter the pressure at the %s rarefaction %s edge: ",
	       message,
	       (type == TRAILING_EDGE) ? "trailing" : "leading");
	(void) Scanf("%f\n",&pr);
	if (is_forward_wave(inner->wv_type))
	{
	    rw1d->l_or_r = RIGHT_FAMILY;
	    sgn = 1.0;
	}
	else
	{
	    rw1d->l_or_r = LEFT_FAMILY;
	    sgn = -1.0;
	}
	if (type == TRAILING_EDGE)
	{
	    state_on_adiabat_with_pr(stl,pr,stt,TGAS_STATE);
	    vl = vel(0,stl);
	    Vel(stt)[0] = vt = vl + sgn*riemann_wave_curve(stl,pr);
	}
	else
	{
	    state_on_adiabat_with_pr(stt,pr,stl,TGAS_STATE);
	    pr = pressure(stt);
	    vt = vel(0,stt);
	    Vel(stl)[0] = vl = vt - sgn*riemann_wave_curve(stl,pr);
	}
	rw1d->spl = spl = vl + sgn*sound_speed(stl);
	rw1d->spt = spt = vt + sgn*sound_speed(stt);
	if (is_rarefaction_leading_edge(inner->wv_type))
	{
	    rl = max_radii(inner);
	    rt = max_radii(outer);
	    rw1d->el_lead = inner;
	    rw1d->el_trail = outer;
	}
	else
	{
	    rl = max_radii(outer);
	    rt = max_radii(inner);
	    rw1d->el_lead = outer;
	    rw1d->el_trail = inner;
	}
	screen("Enter the leading edge position (dflt = %g): ",rl);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&rl);
	screen("Enter the trailing edge position (dflt = %g): ",rt);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&rt);
	rw1d->zl = rl;
	rw1d->zt = rt;
	if (rl < rt)
	{
	    rw1d->zmin = rl;
	    rw1d->stmin = stl;
	    rw1d->zmax = rt;
	    rw1d->stmax = stt;
	}
	else
	{
	    rw1d->zmin = rt;
	    rw1d->stmin = stt;
	    rw1d->zmax = rl;
	    rw1d->stmax = stl;
	}
	rw1d->tbar = (rl - rt)/(spl - spt);
	rw1d->zbar = 0.5*(rl + rt - (spl + spt)*rw1d->tbar);
	screen("Enter the wave space-time center r and the time elapsed since "
	       "\n\tthe origination of the wave (dflt = %g %g): ",rw1d->zbar,
						                  rw1d->tbar);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    static const char *fmt = "%lf %lf";
	    if (sscanf(s,fmt,&rw1d->zbar,&rw1d->tbar) != 2)
	    {
		screen("ERROR in prompt_for_rarefaction_edge(), "
		       "invalid input of wave space-time center\n");
		clean_up(ERROR);
	    }
	}
}		/*end prompt_for_rarefaction_edge*/

EXPORT	double	max_radii(
	ELLIPSOID *ellip)
{
	double r;
	int   i, dim = ellip->dim;
	for (r = 0.0, i = 0; i < dim; ++i)
	    r = max(r,ellip->rad[i]);
	return ellip->scale*r;
}		/*end max_radii*/

LOCAL	boolean nested_ellipsoids(
	ELLIPSOID *inner,
	ELLIPSOID *outer)
{
	double *cin;
	double *rin;
	double **Qin;
	double p[3];
	int   i, j, dim;

	if ((inner == NULL) || (outer == NULL))
	    return NO;

	dim = inner->dim;
	cin = inner->cen;
	rin = inner->rad;
	Qin = inner->Q;

	if (inside_ellipsoid(cin,outer) == NO)
	    return NO;

	for (i = 0; i < dim; ++i)
	{
	    for (j = 0; j < dim; ++j)
		p[j] = cin[j] + rin[i]*Qin[j][i];
	    if (inside_ellipsoid(p,outer) == NO)
	        return NO;
	    for (j = 0; j < dim; ++j)
		p[j] = cin[j] - rin[i]*Qin[j][i];
	    if (inside_ellipsoid(p,outer) == NO)
	        return NO;
	}

	/* The convex hull of the extreme points of the inner ellipsoid
	 * are inside the outer ellipsoid,  this is enough for now.
	 */
	return YES;
}		/*end nested_ellipsoids*/

LOCAL	boolean inside_ellipsoid(
	double     *p,
	ELLIPSOID *ellip)
{
	double *cen = ellip->cen;
	double *rad = ellip->rad;
	double **Q = ellip->Q;
	double l, dp[3];
	int   i, j, dim;

	dim = ellip->dim;
	for (i = 0; i < dim; ++i)
	{
	    for (dp[i] = 0.0, j = 0; j < dim; ++j)
	        dp[i] += (p[j] - cen[j])*Q[j][i];
	}

	for (l = 0.0, i = 0; i < dim; ++i)
	    l += sqr(dp[i]/rad[i]);
	return (l < 1.0) ? YES : NO;
}		/*end nested_ellipsoids*/

LOCAL	void input_spherical_modes(
	int		offset,
	int		min_n,
	int		max_n,
	int		nmodes,
	FOURIER_POLY	*fpoly,
	double		*L,
	double		*U)
{
	double		**wv_num = fpoly->nu+offset, *A = fpoly->A + offset;
	double		*phase = fpoly->phase + offset;
	int		l,m,i,pl,pm;

	pm = irint(2.0*PI/(U[0] - L[0]));
	pl = irint(PI/(U[1] - L[1]));

	(void) printf("\tnumber of modes::%d\n",nmodes);
	i = 0;
	for (l = min_n; l <= max_n; l = l+pl)
	{
	    for (m = 0; m <= l; m = m+pm)
            {
	        (void) printf("Prompt for spherical mode %d\n",i);
	        (void) printf("\tfrequency for the azimuth mode %d ::%d\n",i,l);
		(void) printf("\tfrequency for the polar angle mode %d ::%d\n",
			      i,m);

		screen("Enter amplitude and phase (deg.) for mode %d: ", i);
		(void) Scanf("%f %f\n",&A[i],&phase[i]);
		phase[i] = radians(phase[i]);
		wv_num[i][0] = l;
		wv_num[i][1] = m;
             	++i;
            }
	}
}		/*end input_spherical_modes*/


LOCAL void get_state_1dtable(
      	double		*coords,
	Locstate	state,
	COMP_TYPE	*ct,
	HYPER_SURF	*hs,
	INTERFACE	*intfc,
	INIT_DATA	*init,
	int		type)
{
	ONED_TABLE    	*table = One_d_table(ct);
	int	    	nr = table->nx;
	int	    	i;
	double	    	r;
	double	    	f0, f1;
	double         	vr;
	double	    	n[2];

	debug_print("init_states","Entered get_state_1dtable()\n");
	if (ct->type != ONE_DIMENSIONAL_TABLE)
        {
	    screen("ERROR in get_state_1dtable(), "
	           "inconsistent comp_type->type.\n");
	    clean_up(ERROR);
	}
      
      	set_type_of_state(state, TGAS_STATE);

      	r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]);
      	r = max(r, 0.001);
 
      	n[0] = coords[0]/r;
      	n[1] = coords[1]/r;


      	if (r <= (table->pt)[0])
        {  
            Press(state) = (table->p)[0];
            Dens(state)  = (table->d)[0];
	    vr = (table->v)[0];
	    Vel(state)[0] = vr*n[0];
	    Vel(state)[1] = vr*n[1];
	}
      	else if (r >= (table->pt)[nr-1])
       	{
            Press(state) = (table->p)[nr-1];
            Dens(state)  = (table->d)[nr-1];
	    vr = (table->v)[nr-1];
	    Vel(state)[0] = vr*n[0];
	    Vel(state)[1] = vr*n[1];  
	}
      	else
       	{
            i = 1;
	    while ((r > (table->pt)[i]) & (i < nr))
		 i++;
            f0 = ((table->pt)[i] - r)/((table->pt)[i] - (table->pt)[i-1]);
	    f1 = 1. - f0;
            Press(state) = f0*(table->p)[i-1] + f1*(table->p)[i];
            Dens(state) = f0*(table->d)[i-1] + f1*(table->d)[i];
	    vr = f0*(table->v)[i-1] + f1*(table->v)[i];
	    Vel(state)[0] = vr*n[0];
	    Vel(state)[1] = vr*n[1];  
	}
      	Params(state) = ct->params;
      	set_state(state,GAS_STATE,state);

}     	/*end get_state_1dtable*/

LOCAL	void	free_1dtable_comp_type(
	COMP_TYPE	*comp_type)
{
	ONED_TABLE *table;
	if (comp_type->type != ONE_DIMENSIONAL_TABLE)
		return;

	table = One_d_table(comp_type);
	if (table == NULL)
		return;
        if (table->p != NULL)
         free(table->p);
        if (table->d != NULL)
         free(table->d);
        if (table->v != NULL)
         free(table->v);
        if (table->pt != NULL)
        free(table->pt);
	
	free(table);
	comp_type->extra = NULL;
}		/*end free_1dtable_comp_type*/


EXPORT	void prompt_for_rigid_body_params(
	RIGID_BODY_PARAMS *rgb_params,
	const char *mesg,
	int dim)
{
	int 	i;

	screen("Enter the total mass for %s: ",
	       (mesg != NULL) ? mesg : "\0");
	Scanf("%f \n",&rgb_params->total_mass);
	screen("Enter the center of mass for %s: ",
	       (mesg != NULL) ? mesg : "\0");
	for (i = 0; i < dim; ++i) 
	    Scanf("%f ",&rgb_params->center_of_mass[i]);
	Scanf("\n");
	screen("Enter the moment of inertial for %s: ",
	       (mesg != NULL) ? mesg : "\0");
	Scanf("%f \n",&rgb_params->mom_of_inertial);
}	/* end prompt_for_rigid_body_params */
