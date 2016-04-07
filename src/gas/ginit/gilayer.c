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
*				gilayer.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <ginit/ginit.h>

 
	/* LOCAL Function Declarations */
LOCAL	LAYER_SYS	*prompt_for_layer_sys(INIT_PHYSICS*);
LOCAL	LAYER_FLAG	set_layer_sys_flag(int);
LOCAL	LAYER_SURF	*bndry_surf(double,RECT_GRID*);
LOCAL	LAYER_SURF	*prompt_for_layer_surf(Front*,LAYER*,const LAYER_FLAG*,
					       INIT_PHYSICS*);
LOCAL	Locstate	st_info_acrx_intfc(int,int,int,int,double*,double*,SIDE,
					   Locstate,int,Locstate,LAYER_SYS*,
					   INIT_PHYSICS*,INIT_DATA*);
LOCAL	UT_SHOCK	*alloc_untrack_shock_structure(Front*);
LOCAL	boolean	comp_type_processed(LAYER_SYS*,int,int);
LOCAL	double	moving_frame_speed(LAYER_SYS*);
LOCAL	int	get_wv_type(int,int,LAYER*);
LOCAL	int	get_ahead_state_coords_etc(LAYER_SYS*,int,int,double*,double*,
                                           SIDE*,Locstate,int*,int*,int,
					   INIT_DATA*);
LOCAL	void	change_ref_frame(LAYER_SYS*,double);
LOCAL	void	check_for_consistent_layer(char*,LAYER*,Front*);
LOCAL	void	get_comp_type_extra(INIT_DATA*,INIT_PHYSICS*,
				    LAYER_SYS*,int,int);
LOCAL	void	get_comp_type_params(INIT_DATA*,INIT_PHYSICS*,LAYER_SYS*,
				     int,int);
LOCAL	void	get_info_on_ellip(RECT_GRID*,char*,
				  ELLIPSOID*,double*,double*,INIT_DATA*);
LOCAL	void	get_info_on_surf(RECT_GRID*,LAYER_SURF*,double*,double*,int);
LOCAL	void	goto_moving_frame(LAYER_SYS*);
LOCAL	void	make_all_intfc(LAYER_SYS*,INIT_DATA*);
LOCAL	void	make_untracked_layer_shock(Front*,LAYER_SURF*,INIT_DATA*);
LOCAL	void	print_rm_scale_factors(LAYER_SYS*);
LOCAL	void	prompt_for_leading_edge_state(Locstate,Locstate,int);
LOCAL	void	prompt_for_trailing_edge_state(Locstate,Locstate,int);
LOCAL	void	region_name(char*,int,int,int,int);
LOCAL	void	set_extra_ambient(COMP_TYPE*,LAYER_SYS*,int,int,
				  INIT_PHYSICS*,INIT_DATA*);
LOCAL	void	set_extra_elliptical(COMP_TYPE*,LAYER_SYS*,int,int,
				     INIT_PHYSICS*,INIT_DATA*);
LOCAL	void	set_extra_random_region(COMP_TYPE*,LAYER_SYS*,int,int,
					INIT_PHYSICS*,INIT_DATA*);
LOCAL	void	set_extra_rarefaction_wave_1d(COMP_TYPE*,LAYER_SYS*,
					      int,int,INIT_PHYSICS*,INIT_DATA*);
LOCAL	void	set_extra_rt_kh(COMP_TYPE*,LAYER_SYS*,int,int,
				INIT_PHYSICS*,INIT_DATA*);
LOCAL	void	set_extra_stretching(COMP_TYPE*,LAYER_SYS*,int,int,INIT_DATA*);
LOCAL	void	set_extra_trans_layer(COMP_TYPE*,LAYER_SYS*,int,int,
				      INIT_PHYSICS*,INIT_DATA*);
LOCAL	void	set_untracked_surface(LAYER_SURF*,Front*);
LOCAL	boolean	input_coordinates(double*,int);
LOCAL	void	free_1d_overlay_comp_type(COMP_TYPE*);
LOCAL	void	set_extra_1d_overlay(INIT_DATA*,COMP_TYPE*,
				     INIT_PHYSICS*,LAYER_SYS*);


/*		LAYER AND REGION LABELING CONVENTION
*
*	Layer_label goes from 1 to num_layers from bottom to top.
*	Surface_label goes from 1 to num_layers-1 from bottom to top.
*	The ith surface separates the ith layer and the (i+1)th layer.
*
*	The notion of ellipsoids and regions is only valid within each layer.
*	Region_label starts from 0.
*	Ellipsoid_label starts from 1.
*	The region immediately within ellipsoid i (i = 1,...,) is labeled  i;
*	The region out of EVERY ellipsoid is labeled 0;
*
*	An ellipsoid may NOT contain any other ellipsoids whose
*	ellipsoid_labels are smaller than the label of itself!!!
*/

EXPORT	void 	init_multi_layer(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	LAYER_SYS	*layer_sys;

	debug_print("layer","Entered init_multi_layer()\n");

	supports_riemann_problem_waves(ip) = YES;
	layer_sys = prompt_for_layer_sys(ip);

	init_comp_type(init,ip,layer_sys);
	
	if (MOVING_FRAME(layer_sys->flag) == YES)
	    goto_moving_frame(layer_sys);

	switch (problem_type(ip))
	{
	case RICHTMYER_LINEAR_THEORY:
	    rm_linear(init,ip,layer_sys);           /* never returns */
	    break;

	case RAYLEIGH_TAYLOR:
	    if (RT_LINEAR_THEORY(layer_sys->flag) == YES)
	        rt_ml_linear_pert(layer_sys,init);
	    break;

	default:
	    break;
	}

	make_all_intfc(layer_sys,init);

	if (debugging("layer"))
	{
	    (void) printf("Initial interface\n");
	    print_interface(ip->root->front->interf);
	}

	if (problem_type(ip) == RICHTMYER_MESHKOV)
	    print_rm_scale_factors(layer_sys);

	debug_print("layer","Left init_multi_layer()\n");
}		/*end init_multi_layer*/


LOCAL	LAYER_SYS	*prompt_for_layer_sys(
	INIT_PHYSICS	*ip)
{
	Front	          *front = ip->root->front;
	int 	          i, num_layers;
	LAYER_FLAG        flag;
	RECT_GRID         *rect_grid = front->rect_grid;
	int	          dim = rect_grid->dim;
	char 	          mesg[100]; 
	LAYER	          **layer;
	LAYER_SYS         *layer_sys;
	static const char *xyz[] = { "x", "y", "z"};

	flag = set_layer_sys_flag(problem_type(ip));

	screen("\nEnter number of layers: ");
	(void) Scanf("%d\n",&num_layers);
	if (num_layers > 1)
	{
	    screen("\nThe system consists of %d layer(s) "
	           "and %d interface(s), both the layer \n"
	           "label and the interface label increase in the "
	           "positive %s direction.\n\n",
	           num_layers,num_layers-1,xyz[dim-1]);
	}

	uni_array(&layer,5*num_layers+1,sizeof(LAYER *));  /* layer[0] unused */
							/* Each layer can
							 * be expanded to as
							 * many as five 
							 * layers to
							 * accomodate Riemann
							 * problem solutions
							 */
	for (i = 1; i <= num_layers; ++i)
	{
	    layer[i] = prompt_for_each_layer(i,num_layers,&flag,front,layer,ip);

	    if (i == num_layers)
	        break;

	    /* check avg_position */
	    (void) sprintf(mesg,"for the average position of the %d%s surface",
			   i, ordinal_suffix(i));
	    check_for_consistent_layer(mesg,layer[i],front);
	}

	/* Set components of layer surfs */
	for (i = 1; i <= num_layers; ++i)
	{
	    layer[i]->upper_surf->l_comp = layer[i]->comp;
	    layer[i]->lower_surf->r_comp = layer[i]->comp;
	}

	scalar(&layer_sys,sizeof(LAYER_SYS));
	layer_sys->dt = HUGE_VAL;
	layer_sys->num_layers = num_layers;
	layer_sys->layer = layer;
	layer_sys->front = front;
	layer_sys->flag = flag;
	return	layer_sys;
}		/*end prompt_for_layer_sys*/

LOCAL	void	check_for_consistent_layer(
	char	*mesg,
	LAYER	*layer,
	Front	*front)
{
	LAYER_SURF	*lower_surf = layer->lower_surf;
	LAYER_SURF	*upper_surf = layer->upper_surf;
	double	*L = front->rect_grid->L, *U = front->rect_grid->U;
	double	cp, vp[3];
	int	dim = front->rect_grid->dim;
	int	k = dim-1;

	cp = vector_product(upper_surf->nor,lower_surf->nor,vp,dim);
	if (cp <= EPSILON*mag_vector(front->rect_grid->h,dim))
	{
	    /*parallel curves or surfaces*/

	    check_float_input(mesg,upper_surf->pbar[k],lower_surf->pbar[k],
			      U[k],GE_AND_LE);
	    return;
	}
	if (dim == 2)
	{
	    double	x, y;
	    double	l, u;

	    u = scalar_product(upper_surf->nor,upper_surf->pbar,dim);
	    l = scalar_product(lower_surf->nor,lower_surf->pbar,dim);
	    x = (u*lower_surf->nor[1] - l*upper_surf->nor[1])/vp[0];
	    y = (l*upper_surf->nor[0] - u*lower_surf->nor[0])/vp[0];
	    check_float_input(mesg,x,L[0],U[0],LE_OR_GE);
	    check_float_input(mesg,y,L[1],U[1],LE_OR_GE);
	    y = get_surf_height(lower_surf->pbar,upper_surf);
	    check_float_input(mesg,y,lower_surf->pbar[1],U[1],GE_);
	}
	/*TODO: add checking for 3D*/
}		/*end check_for_consistent_layer*/

/*ARGSUSED*/
EXPORT 	void 	make_layer_surf(
	Front		*front,
	LAYER_SURF	*surf,
	INIT_DATA	*init)
{
	RECT_GRID	*rect_grid = front->rect_grid;
	HYPER_SURF	*hs = NULL;
	int	 	dim = rect_grid->dim;

	if ((surf == NULL) || (surf->created))
	    return;

	surf->created = YES;
	if (is_shock_wave(surf->wv_type) && (surf->untracked != NULL))
	    make_untracked_layer_shock(front,surf,init);

	switch (dim)
	{
	case 1:
	{
	    POINT	*pt;

	    pt = make_point(surf->pbar,surf->l_comp,surf->r_comp);

	    hs = Hyper_surf(pt);
	    wave_type(pt) = surf->wv_type;
	    if (surf->untracked != NULL)
	    	untracked_hyper_surf(pt) = YES;
	}
	    break;
	case 2:
	{
	    CURVE 	*cur;
	    NODE 	*ns, *ne;
	    double	dx, x0, x1;
	    double	*L = rect_grid->L, *U = rect_grid->U;
	    double	*pbar = surf->pbar, *nor = surf->nor;
	    double 	coords[MAXD];
	    int		i;
	    int 	num_pts = 2*rect_grid->gmax[0];

	    x0 = L[0]; 	x1 = U[0];

	    coords[0] = x1;	
	    coords[1] = get_surf_height(coords,surf);
	    if (coords[1] < L[1])
	    	x1 = pbar[0] - (nor[1]/nor[0])*(L[1] - pbar[1]);
	    else if (coords[1] > U[1])
	    	x1 = pbar[0] - (nor[1]/nor[0])*(U[1] - pbar[1]);
	    ns = make_node(Point(coords));

	    coords[0] = x0;	
	    coords[1] = get_surf_height(coords,surf);
	    if (coords[1] < L[1])
	    	x0 = pbar[0] - (nor[1]/nor[0])*(L[1] - pbar[1]);
	    else if (coords[1] > U[1])
	    	x0 = pbar[0] - (nor[1]/nor[0])*(U[1] - pbar[1]);
	    ne = make_node(Point(coords));

	    cur = make_curve(surf->l_comp,surf->r_comp,ns,ne);
	    hs = Hyper_surf(cur);

	    if (surf->fpoly != NULL)
	    {
	    	dx = (x1 - x0)/num_pts;
	    	for (i = 1; i < num_pts; ++i)
	    	{
	    	    coords[0] = x0 + i*dx;
	    	    coords[1] = get_surf_height(coords,surf);
	    	    if (insert_point_in_bond(Point(coords),cur->first,cur) !=
			FUNCTION_SUCCEEDED)
	            {
	                screen("ERROR in make_layer_surf(), "
		               "insert_point_in_bond() failed\n");
	                clean_up(ERROR);
	            }
	    	}
		do_not_redistribute(cur) = YES;
	    }
	
	    set_is_bdry(ns);	set_is_bdry(ne);
	    start_status(cur) = end_status(cur) = INCIDENT;
	    wave_type(cur) = surf->wv_type;
	    surface_tension(cur) = surf->surf_ten;
	    if (surf->untracked != NULL)
	    	untracked_hyper_surf(cur) = YES;
	}
	    break;
	case 3:
	{
	    BOND		*b;
	    CURVE		**c, *curves[4];
	    HYPER_SURF		*hstmp;
	    HYPER_SURF_ELEMENT	*hse;
	    INTERFACE		*intfc = front->interf;
	    NODE		**n, *nodes[4];
	    POINT		*p;
	    SURFACE		*s;
	    double		*L = computational_grid(intfc)->VL;
	    double		*U = computational_grid(intfc)->VU;
	    double		coords[MAXD];
	    int			i;

	    if (surf->untracked != NULL)
	        return;

	    /*TODO: code needed for 3D rotated surface*/
	    s = make_surface(surf->l_comp,surf->r_comp,NULL,NULL);
	    hs = Hyper_surf(s);
	    wave_type(s) = surf->wv_type;
	    coords[0] = L[0];
	    coords[1] = L[1];
	    coords[2] = surf->pbar[2];
	    nodes[0] = make_node(Point(coords));
	    coords[0] = U[0];
	    nodes[1] = make_node(Point(coords));
	    coords[1] = U[1];
	    nodes[2] = make_node(Point(coords));
	    coords[0] = L[0];
	    nodes[3] = make_node(Point(coords));
	    for (i = 0; i < 4; ++i)    
	    {
	        curves[i] = make_curve(NO_COMP,NO_COMP,nodes[i],nodes[(i+1)%4]);
		install_curve_in_surface_bdry(s,curves[i],POSITIVE_ORIENTATION);
	    }
	    planar_surface_triangulation(s,rect_grid,YES);

	    /* Perturb points on the interface */

	    if (surf->fpoly != NULL)
	    {
	        for (n = intfc->nodes; n && *n; ++n)
	        {
	            p = (*n)->posn;
	            Coords(p)[dim-1] = get_surf_height(Coords(p),surf);
	        }
		for (c = intfc->curves; c && *c; ++c)
		{
		    for (b = (*c)->first; b != (*c)->last; b = b->next)
		    {
		        p = b->end;
		        Coords(p)[dim-1] = get_surf_height(Coords(p),surf);
		    }
		}
		(void) next_point(intfc,NULL,NULL,NULL);
		while(next_point(intfc,&p,&hse,&hstmp))
		{
		    if (hs != hstmp)
			continue;
		    if (!Boundary_point(p))
		    	Coords(p)[dim-1] = get_surf_height(Coords(p),surf);
		}
	    }
	    reset_normal_on_intfc(intfc);
	    surface_tension(s) = surf->surf_ten;
	}
	    break;
	}
	if (is_scalar_wave(surf->wv_type))
	    layer_index(hs) = ++num_layers(hs->interface);
}		/*end make_layer_surf*/

EXPORT 	void 	make_untracked_layer_shock(
	Front		*front,
	LAYER_SURF	*surf,
	INIT_DATA	*init)
{
	UT_SHOCK	*utsw = surf->untracked;
	COMP_TYPE	*ctypel, *ctyper;
	RECT_GRID	*gr = front->rect_grid;
	int		i, dim = gr->dim;

	if (utsw->width < front->rect_grid->h[dim-1])
	    return;

	ctypel = comp_type(surf->l_comp);
	ctyper = comp_type(surf->r_comp);

	for (i = 0; i < dim; ++i)
	    utsw->posn[i] = surf->pbar[i];
	scalar(&utsw->ctype0,sizeof(COMP_TYPE));
	scalar(&utsw->ctype1,sizeof(COMP_TYPE));
	if (is_forward_wave(surf->wv_type))
	{
	    Get_state(utsw->posn,utsw->state0,ctyper,NULL,front->interf,
		      init,GAS_STATE);
	    Get_state(utsw->posn,utsw->state1,ctypel,NULL,front->interf,
		      init,GAS_STATE);
	    *utsw->ctype0 = *ctyper;
	    *utsw->ctype1 = *ctypel;
	    for (i = 0; i < dim; ++i)
	    	utsw->nor[i] = surf->nor[i];
	}
	else
	{
	    Get_state(utsw->posn,utsw->state0,ctypel,NULL,front->interf,
		      init,GAS_STATE);
	    Get_state(utsw->posn,utsw->state1,ctyper,NULL,front->interf,
		      init,GAS_STATE);
	    *utsw->ctype0 = *ctypel;
	    *utsw->ctype1 = *ctyper;
	    for (i = 0; i < dim; ++i)
	    	utsw->nor[i] = -surf->nor[i];
	}
	utsw->_wave_type = surf->wv_type;
	set_untracked_shock_wave_comp_type(ctypel,utsw,front);
	set_untracked_shock_wave_comp_type(ctyper,utsw,front);
}		/*end make_untracked_layer_shock*/

LOCAL  void    make_all_intfc(
	LAYER_SYS	*layer_sys,
	INIT_DATA	*init)
{
        int             i;
        int     	num_layers = layer_sys->num_layers;
        Front   	*front = layer_sys->front;
        LAYER   	**layer = layer_sys->layer;
 
        for (i = 1; i <= num_layers; ++i)
        {
            ELLIPSOID       *el;
            int             j, num_ellips = layer[i]->num_ellips;
            for (j = 1; j <= num_ellips; ++j)
            {
                el = layer[i]->ellip[j];
	        if (is_scalar_wave(el->wv_type))
		    el->layer_index = ++num_layers(front->interf);
		(void) make_ellipsoid(el,el->compin,el->compout,front);
            }
            make_layer_surf(front,layer[i]->lower_surf,init);
            make_layer_surf(front,layer[i]->upper_surf,init);
        }
}		/*end make_all_intfc*/
 

LOCAL	LAYER_FLAG set_layer_sys_flag(
	int        problem)
{
	LAYER_FLAG Flag;
	char		choice[Gets_BUF_SIZE];

	CLEAR_LAYER_FLAG(Flag);
	switch (problem)
	{
	case	RAYLEIGH_TAYLOR:
	    ALL_CONTACT(Flag) = YES;
	    screen("Use RT linearized perturbation theory? (y, dlft): ");
	    (void) Gets(choice);
	    if ((choice[0] != 'n') && (choice[0] != 'N'))
		RT_LINEAR_THEORY(Flag) = YES;
	    break;
	case	BUBBLES_DROPS:
	case	SHOCKED_THERMAL_LAYER:
	case	EXPANDING_SHELLS:
	    HAS_ELLIPSOID(Flag) = ELLIPTICAL_REGION(Flag) = YES;
	    break;
	case	SHOCK_JET:
	    HAS_ELLIPSOID(Flag) = MOVING_FRAME(Flag) = YES;
	    break;
	case	RICHTMYER_MESHKOV:
	    MOVING_FRAME(Flag) = YES;
	    break;
	default:
	    break;
	}
	return Flag;
}		/*end set_layer_sys_flag*/


EXPORT	LAYER	*prompt_for_each_layer(
	int		 layer_label,
	int		 num_layers,
	const LAYER_FLAG *flag,
	Front		 *front,
	LAYER		 **layer,
	INIT_PHYSICS	 *ip)
{
	RECT_GRID	*rect_grid = front->rect_grid;
	int		dim = rect_grid->dim;
	char		s[Gets_BUF_SIZE];
	LAYER		*lyr;
	char 		mesg[100]; 
	const char 	*elname = (dim == 3) ? "ellipsoid" : "ellipse";
	int		j, num_ellips;
	ELLIPSOID	**ellip;

	scalar(&lyr, sizeof(LAYER));

	lyr->layer_label = layer_label;
	if (layer_label == 1)
	    lyr->comp = FIRST_DYNAMIC_COMPONENT;
	else
	{
	    lyr->comp = layer[layer_label-1]->comp +
			layer[layer_label-1]->num_ellips + 1;
	    lyr->prev = layer[layer_label-1];
	    layer[layer_label-1]->next = lyr;
	}

	screen("Enter the component label for layer %d (default = %d): ",
		lyr->layer_label,lyr->comp);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    (void) sscanf(s,"%d",&lyr->comp);
	    if (lyr->comp < FIRST_DYNAMIC_COMPONENT)
	    {
	    	(void) printf("WARNING in prompt_for_each_layer(), "
	    	              "component entered in less than "
	    	              "FIRST_DYNAMIC_COMPONENT\n"
	    	              "Is this correct? (n): ");
	    	(void) Gets(s);
	    	if (s[0] != 'y' || s[0] != 'Y')
	    	{
	    	    screen("ERROR in prompt_for_each_layer(), "
	    		   "invalid component\n");
	    	    clean_up(ERROR);
	    	}
	    }
	}
	screen("Layer %d has component label %d\n",lyr->layer_label,lyr->comp);
	new_component(lyr->comp);

	lyr->num_ellips = 0;
	num_ellips = lyr->num_ellips;
	if (has_ellipsoid(flag) == YES)
	{
	    (void) sprintf(mesg,"in the %d%s layer",
	    	           layer_label,ordinal_suffix(layer_label));
	    screen("Enter the number of %ss %s (default = %d): ",
	    	   elname, mesg, num_ellips);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d",&num_ellips);
	    lyr->num_ellips = num_ellips;
	}
	if (num_ellips > 0)
	{
	    uni_array(&ellip,5*num_ellips+1,sizeof(ELLIPSOID *));
	    lyr->ellip = ellip;
	}
	else
	    lyr->ellip = NULL;
	for (j = 1; j <= num_ellips ; ++j)
	{
	    (void) sprintf(mesg,"%d%s",j,ordinal_suffix(j));
	    ellip[j] = prompt_for_ellipsoid(front,mesg,NULL,NULL,flag,ip);
	    ellip[j]->compin = lyr->comp+j;
	    screen("Enter the component number for the region inside\n"
	           "\t%s number %d (default = %d): ",elname,j,ellip[j]->compin);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    	(void) sscanf(s,"%d",&ellip[j]->compin);
	    new_component(ellip[j]->compin);
	    ellip[j]->compout = lyr->comp;
	    if (j > 1)
	    {
	    	screen("Enter the component label of the region that"
	    	       " immediately contains the %s %s"
	    	       " (default = %d): ",
	    	       mesg,elname,ellip[j]->compout);
	    	(void) Gets(s);
	    	if (s[0] != '\0')
	    	    (void) sscanf(s,"%d",&ellip[j]->compout);
	    }
	}

	if (layer_label == 1)
	    lyr->lower_surf = bndry_surf(rect_grid->L[dim-1],rect_grid);
	else
	    lyr->lower_surf = layer[layer_label-1]->upper_surf;

	if (layer_label == num_layers)
	    lyr->upper_surf = bndry_surf(rect_grid->U[dim-1],rect_grid);
	else
	    lyr->upper_surf = prompt_for_layer_surf(front,lyr,flag,ip);

	return	lyr;
}		/*end prompt_for_each_layer*/


LOCAL	LAYER_SURF	*bndry_surf(
	double		height,
	RECT_GRID	*gr)
{
	LAYER_SURF	*surf;
	double		*L = gr->L, *U = gr->U;
	int		i, dim = gr->dim;
	int		k = dim - 1;

	surf = alloc_layer_surf();
	surf->dim = gr->dim;
	for (i = 0; i < k; ++i)
	{
	    surf->pbar[i] = 0.5*(L[i] + U[i]);
	    surf->nor[i] = 0.0;
	}
	surf->pbar[k] = height;
	surf->nor[k] = 1.0;
	surf->created = YES;

	return	surf;
}		/*end bndry_surf*/


EXPORT	LAYER_SURF	*alloc_layer_surf(void)
{
	LAYER_SURF	*surf;
	int		i;

	scalar(&surf, sizeof(LAYER_SURF)); 
	surf->l_comp = surf->r_comp = ERROR;
	for (i = 0; i < 3; ++i)
	{
	    surf->pbar[i] = HUGE_VAL;
	    surf->velocity[i] = -HUGE_VAL;
	}
	surf->s_max = -HUGE_VAL;
	surf->s_min = HUGE_VAL;
	for (i = 0; i < 3; ++i)
	    surf->nor[i] = HUGE_VAL;
	surf->fpoly = NULL;
	surf->wv_type = UNKNOWN_WAVE_TYPE;
	surf->dim = -1;
	surf->surf_ten = 0.0;
	surf->untracked = NULL;
	surf->reset_position = NO;
	surf->layer_index = 0;
	surf->created = NO;

	return	surf;
}		/*end alloc_layer_surf*/


LOCAL   LAYER_SURF *prompt_for_layer_surf(
	Front		 *front,
	LAYER		 *layer,
	const LAYER_FLAG *flag,
	INIT_PHYSICS	 *ip)
{
	IMPORT boolean	suppress_prompts;
	int		il = layer->layer_label;
	int		dim = front->rect_grid->dim;
	int		j, wv_type;
	boolean		flag_add_perturbation, rotated;
	char		mesg[100], mesg2[100], choice[Gets_BUF_SIZE];
	static		double	dflt_nor[4][3] = { {0.0, 0.0, 0.0},
						   {1.0, 0.0, 0.0},
						   {0.0, 1.0, 0.0},
						   {0.0, 0.0, 1.0} };
	LAYER_SURF	*surf;

	surf = alloc_layer_surf();
	surf->dim = dim;

	(void) sprintf(mesg," for the %d%s surface",il,ordinal_suffix(il));

	flag_add_perturbation =
		((rt_linear_theory(flag) == YES) && (il == 1)) ? YES : NO;

	screen("Input either the average height %s\n\t"
	       "or the coordinates of a point on this surface: ",mesg);
	(void) Gets(choice);
	if (choice[0] != '\0')
	{
	    double	*U = front->rect_grid->U, *L = front->rect_grid->L;
	    int	i;
	    char	*c;

	    for (j = 0, c = strtok(choice," \t"); j < dim && c != NULL;
						++j, c = strtok(NULL," \t"))
		(void) sscan_float(c,surf->pbar+j);
	    if (j == 1)
	    {
	    	surf->pbar[dim-1] = surf->pbar[0];
	    	for (i = 0; i < dim-1; ++i)
	    	    surf->pbar[i] = 0.5*(L[i] + U[i]);
	    }
	    else if (j < dim)
	    {
	    	screen("ERROR in prompt_for_layer_surf(), "
	    	       "improper input of surf coords\n");
	    	clean_up(ERROR);
	    }
	}
	else
	{
	    screen("ERROR in prompt_for_layer_surf(), "
	           "improper input of surf coords\n");
	    clean_up(ERROR);
	}

	surf->l_comp = NO_COMP;
	surf->r_comp = NO_COMP;

	for (j = 0; j < dim; ++j)
	    surf->nor[j] = dflt_nor[dim][j];

	rotated = NO;
	if (dim > 1)
	{
	    screen("Enter ");
	    if (dim == 2)
	    {
	    	screen("either the angle (in degrees) that the\n\t");
	    	screen("line makes with the x-axis (dflt = 0) or\n\t");
	    }
	    screen("the average normal to the surface ");
	    print_general_vector("[default = ",surf->nor,dim,"]: ");
	    if (suppress_prompts == NO)
	    	fprint_general_vector(stderr,"[default = ",surf->nor,dim,"]: ");
	    (void) Gets(choice);
	    if (choice[0] != '\0')
	    {
	    	double	mag;
	    	char	*c;

		for (j = 0, c = strtok(choice," \t"); j < dim && c != NULL;
    						++j, c = strtok(NULL," \t"))
    		    (void) sscan_float(c,surf->nor+j);
		if ((j == 1) && (dim == 2))
		{
		    double theta = radians(surf->nor[0]) + 0.5*PI;
		    if (surf->nor[0] == 0)
		    {
		    	surf->nor[1] = 1.0;
		    }
		    else
		    {
		    	surf->nor[0] = cos(theta);
		    	surf->nor[1] = sin(theta);
		    }
		}
    		else if (j < dim)
    		{
    	            screen("ERROR in prompt_for_layer_surf(), "
    	                   "improper normal uni_array\n");
    	            clean_up(ERROR);
	    	}
	    	mag = mag_vector(surf->nor,dim);
	    	for (j = 0; j < dim; ++j)
	    		surf->nor[j] /= mag;
	    }
	}

	/*currently you can not both rotate and perturb a surface*/
	for (j = 0; j < dim; ++j)
	    if (surf->nor[j] != dflt_nor[dim][j])
	    	rotated = YES;
	if (rotated == YES)
	    flag_add_perturbation = NO;

	(void) sprintf(mesg2,"%s with normal",mesg);
	sprint_general_vector(mesg2+strlen(mesg2)," ",surf->nor,dim,". ");
	if (all_contact(flag) == YES)
	    surf->wv_type = wv_type = CONTACT;
	else if (surf->wv_type == UNKNOWN_WAVE_TYPE)
	    surf->wv_type = wv_type =
	        prompt_for_wave_type(mesg2,front->interf,ip);
	else
	    wv_type = surf->wv_type;
	surf->surf_ten = prompt_for_surface_tension(wv_type, mesg);

	surf->fpoly = NULL;

	if ((rt_linear_theory(flag) != YES) && is_scalar_wave(wv_type) &&
	    (dim > 1) && (rotated == NO))
	{
	    screen("Type y if you want to add small perturbations to the "
	           "%d%s interface (dflt=%s): ",il,ordinal_suffix(il),
	    	   (flag_add_perturbation==YES)?"yes":"no");
	    (void) Gets(choice);
	    if ((choice[0] == 'y') || (choice[0] == 'Y'))
	    	flag_add_perturbation = YES;
	    else if ((choice[0] == 'n') || (choice[0] == 'N'))
	    	flag_add_perturbation = NO;
	}

	screen("Type y to turn off tracking %s: ", mesg);
	(void) Gets(choice);
	if ((choice[0] == 'y') || (choice[0] == 'Y'))
	    set_untracked_surface(surf,front);
	if (flag_add_perturbation == YES)
	{
	    surf->fpoly = get_fourier_coeffs(front->rect_grid->GL,
	    				     front->rect_grid->GU,
					     front->rect_grid->dim,mesg);
	    surf->fpoly->z0 = surf->pbar[dim-1];
	}
	screen("\n");

	return	surf;
}		/*end prompt_for_layer_surf*/

LOCAL	UT_SHOCK	*alloc_untrack_shock_structure(
	Front	*front)
{
	UT_SHOCK *utsw;
	scalar(&utsw,sizeof(UT_SHOCK));
	alloc_state(front->interf,&utsw->state0,front->sizest);
	alloc_state(front->interf,&utsw->state1,front->sizest);
	utsw->free_with_comp_type = NO;
	return utsw;
}		/*end alloc_untrack_shock_structure*/

LOCAL	void	set_untracked_surface(
	LAYER_SURF      *surf,
	Front		*front)
{
	RECT_GRID	*gr = front->rect_grid;
	char		choice[80];
	int		n;
	int		dim = gr->dim;

	surf->untracked = alloc_untrack_shock_structure(front);

	surf->untracked->width = 0;
	n = 0;
	screen("Enter the width of the untracked wave ");
	screen("in mesh units (dflt = %d): ",n);
	(void) Gets(choice);
	if (choice[0] != '\0')
		(void) sscanf(choice,"%d",&n);
	surf->untracked->width = n*fabs(scalar_product(gr->h,surf->nor,dim));
}		/*set_untracked_surface*/

EXPORT 	void    init_comp_type(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	LAYER_SYS	*layer_sys)
{
	int		i, dim, layer_label, ellip_label;
	int		num_layers = layer_sys->num_layers;
	int		num_ellips;
	double           dt;
	COMP_TYPE       *ct;
	LAYER		**layer = layer_sys->layer;

	/* prompt for all EOS models */
	(void) prompt_for_eos_params(init,ip,YES,"");

	for (i = 1; i <= num_layers; ++i)
	{
	    screen("\nPlease enter a layer label for prompting: ");
	    (void) Scanf("%d\n",&layer_label);
	    check_int_input("for the layer label",layer_label,1,num_layers,
			    GE_AND_LE);
	    num_ellips = layer[layer_label]->num_ellips;
	    for (ellip_label = 0; ellip_label <= num_ellips; ++ellip_label)
	    {
		if (comp_type_processed(layer_sys,layer_label,ellip_label)==YES)
	            continue;
		(void) get_comp_type_type(layer_sys,layer_label,
					  ellip_label,ip,init);
		get_comp_type_params(init,ip,layer_sys,layer_label,ellip_label);
		get_comp_type_extra(init,ip,layer_sys,layer_label,ellip_label);
	    }
	}

	/* If Riemann problem waves are set addition action is
	 * required to set up rarefaction layers
	 */

	dt = layer_sys->dt;
	num_layers = layer_sys->num_layers;
	dim = ip->root->front->rect_grid->dim;
	for (layer_label = 1; layer_label <= num_layers; ++layer_label)
	{
	    _RAREFACTION_WAVE_1D *rw1d;
	    LAYER                *lyr = layer_sys->layer[layer_label];
	    ELLIPSOID            *ellip;

	    num_ellips = lyr->num_ellips;
	    for (ellip_label = 1; ellip_label <= num_ellips; ++ellip_label)
	    {
	        ellip = lyr->ellip[ellip_label];
	        ct = comp_type(ellip->compin);
		if (ct->type == ELLIPTICAL)
		{
	            _ELLIPTICAL *el = Elliptical(ct);
		    if ((rw1d = el->rw1d) != NULL)
		    {
		        rw1d->tbar = -dt;
		        rw1d->zl = rw1d->zbar + dt*rw1d->spl;
		        rw1d->zt = rw1d->zbar + dt*rw1d->spt;
		        if (rw1d->zl < rw1d->zt)
		        {
		            rw1d->zmin = rw1d->zl;
		            rw1d->zmax = rw1d->zt;
		        }
		        else
		        {
		            rw1d->zmin = rw1d->zt;
		            rw1d->zmax = rw1d->zl;
		        }
		    }
		}
		if (ellip->reset_position)
		{
		    for (i = 0; i < dim; ++i)
		        ellip->rad[i] += ellip->vr[i]*dt;
		    ellip->reset_position = NO;
		}
	    }
	    ct = comp_type(lyr->comp);
	    if (ct->type == RAREFACTION_WAVE_1D)
	    {
	        rw1d = Rarefaction_wave_1d(ct);
		rw1d->tbar = -dt;
		rw1d->zl = rw1d->zbar + dt*rw1d->spl;
		rw1d->zt = rw1d->zbar + dt*rw1d->spt;
		if (rw1d->zl < rw1d->zt)
		{
		    rw1d->zmin = rw1d->zl;
		    rw1d->zmax = rw1d->zt;
		}
		else
		{
		    rw1d->zmin = rw1d->zt;
		    rw1d->zmax = rw1d->zl;
		}
	    }
	    if (lyr->upper_surf->reset_position)
	    {
	        for (i = 0; i < dim; ++i)
		    lyr->upper_surf->pbar[i] += lyr->upper_surf->velocity[i]*dt;
	        lyr->upper_surf->reset_position = NO;
	    }
	    if (lyr->lower_surf->reset_position)
	    {
	        for (i = 0; i < dim; ++i)
		    lyr->lower_surf->pbar[i] += lyr->lower_surf->velocity[i]*dt;
	        lyr->lower_surf->reset_position = NO;
	    }
	}
}		/*end init_comp_type*/


LOCAL	boolean	comp_type_processed(
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label)
{
	LAYER		*lyr = layer_sys->layer[layer_label];
	COMPONENT	compin;
	int		i;
	static COMPONENT	*comps = NULL;
	static int		comp_list_length = 0, num_comps = 0;

	compin = (region_label > 0) ? lyr->ellip[region_label]->compin
				    : lyr->comp;

	if (comps == NULL)
	{
	    comp_list_length = max_num_comps();
	    uni_array(&comps,comp_list_length,sizeof(COMPONENT));
	    num_comps = 0;
	}
	for (i = 0; i < num_comps; ++i)
	    if (comps[i] == compin)
	    	return YES;

	/*new component being processed,  add it to the list */

	if (num_comps >= comp_list_length)
	{
	    COMPONENT	*new_comps;
	    int		new_comp_list_length = 2*comp_list_length;
	    uni_array(&new_comps,comp_list_length,sizeof(COMPONENT));
	    for (i = 0; i < comp_list_length; ++i)
	    	new_comps[i] = comps[i];
	    free(comps);
	    comps = new_comps;
	    comp_list_length = new_comp_list_length;
	}
	comps[num_comps++] = compin;
	return NO;
}		/*end comp_type_processed*/




LOCAL 	void	get_comp_type_params(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label)
{
	LAYER		**layer = layer_sys->layer, *lyr = layer[layer_label];
	COMP_TYPE	*ct;
	COMPONENT	compin;
	int		dim = layer_sys->front->rect_grid->dim;
	int		j, num_ellips = lyr->num_ellips;
	char		name[150];
	Gas_param	*prm;

	compin = (region_label>0) ? lyr->ellip[region_label]->compin:lyr->comp;
	ct = comp_type(compin);

	if (ct->params != NULL)
	    return;

	if (ct->type == OBSTACLE)
	    return;

	/* check nearby layers */
	if ((region_label == 0) &&
	    is_vector_wave(lyr->lower_surf->wv_type) &&
	    (prm = comp_type(layer[layer_label-1]->comp)->params) != NULL)
	{
	    ct->params = prm;
	    return;
	}
	else if ((region_label == 0) &&
		 is_vector_wave(lyr->upper_surf->wv_type) &&
		 (prm = comp_type(layer[layer_label+1]->comp)->params) != NULL)
	{
	    ct->params = prm;
	    return;
	}

	/* check if contains some ellipsoids */
	for (j = 1; j <= num_ellips && lyr->ellip; ++j)
	{
	    if ((lyr->ellip[j]->compout == compin) &&
	        is_vector_wave(lyr->ellip[j]->wv_type) &&
	        ((prm = comp_type(lyr->ellip[j]->compin)->params) != NULL))
	    {
	    	ct->params = prm;
	    	return;
	    }
	}

	/* check outside region of an ellipsoid */
	if ((region_label > 0) &&
	    is_vector_wave(lyr->ellip[region_label]->wv_type) &&
	    ((prm = comp_type(lyr->ellip[region_label]->compout)->params)
	     							!= NULL))
	{
	    ct->params = prm;
	    return;
	}

	region_name(name, layer_label, region_label, dim, num_ellips);
	ct->params = prompt_for_eos_params(init,ip, YES, name); 

	/* check params consitency of an untracked contact */
	if ((region_label == 0) &&
	    is_scalar_wave(lyr->lower_surf->wv_type) &&
	    (lyr->lower_surf->untracked != NULL) &&
	    ((prm = comp_type(layer[layer_label-1]->comp)->params) != NULL) &&
	    (ct->params != prm))
	{
	    screen("\nERROR in get_comp_type_params(), ");
	    screen("The %d%s surface is an untracked contact, ",
		   layer_label-1, ordinal_suffix(layer_label-1));
	    screen("the EOS model on its two sides cannot be different!\n");
	    clean_up(ERROR);
	}
	else if ((region_label == 0) &&
		 is_scalar_wave(lyr->upper_surf->wv_type) &&
		 (lyr->upper_surf->untracked != NULL) &&
		 ((prm = comp_type(layer[layer_label+1]->comp)->params)
		  						!= NULL) &&
		 (ct->params != prm))
	{
	    screen("\nERROR in get_comp_type_params(), ");
	    screen("The %d%s surface is an untracked contact, ",
	           layer_label, ordinal_suffix(layer_label));
	    screen("the EOS model on its two sides cannot be different!\n");
	    clean_up(ERROR);
	}
}		/*end get_comp_type_params*/

/*ARGSUSED*/
LOCAL 	void	get_comp_type_extra(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip,
	LAYER_SYS	*layer_sys,
	int		ll,
	int		rl)
{
	LAYER		*lyr = layer_sys->layer[ll];
	COMP_TYPE	*ct;
	COMPONENT	compin;

	compin = (rl>0) ? lyr->ellip[rl]->compin : lyr->comp;
	ct = comp_type(compin);

	switch (ct->type)
	{
	case	AMBIENT:
	    set_extra_ambient(ct,layer_sys,ll,rl,ip,init);
	    break;
	case	RANDOM_REGION:
	    set_extra_random_region(ct,layer_sys,ll,rl,ip,init);
	    break;
	case 	RT_KH:
	    set_extra_rt_kh(ct,layer_sys,ll,rl,ip,init);
	    break;
	case 	RAREFACTION_WAVE_1D:
	    set_extra_rarefaction_wave_1d(ct,layer_sys,ll,rl,ip,init);
	    break;
	case 	TRANS_LAYER:
	    set_extra_trans_layer(ct,layer_sys,ll,rl,ip,init);
	    break;
	case	ELLIPTICAL:
	    set_extra_elliptical(ct,layer_sys,ll,rl,ip,init);
	    break;
	case	STRETCHING:
	    set_extra_stretching(ct,layer_sys,ll,rl,init);
	    break;
	case    TABULATED_REGION:
	    set_up_read_state_from_file_region(ct," for the layer",
	                                       ip->root->front);
	    break;
	case	ONE_DIMENSIONAL_OVERLAY:
	    set_extra_1d_overlay(init,ct,ip,layer_sys);
	    break;
	case	OBSTACLE:
	    break;
	default:
	    screen("ERROR in get_comp_type_extra(), "
	           "This type of comp_type is not surported.\n");
	    clean_up(ERROR);
	    break;
	}
}		/*end get_comp_type_extra*/


LOCAL	void	set_extra_ambient(
	COMP_TYPE	*ct,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_PHYSICS    *ip,
	INIT_DATA	*init)
{
	Locstate	st = Ambient(ct);
	SIDE            ahead_side;
	int		surf_label, el;
	int		w_type;
	double		coords[MAXD], nor[MAXD];
	Locstate	ahead_st;
	Front		*front = layer_sys->front;

	alloc_state(front->interf,&ahead_st,front->sizest);
			   
	w_type = get_ahead_state_coords_etc(layer_sys,layer_label,region_label,
				            coords,nor,&ahead_side,ahead_st,
					    &surf_label,&el,NO,init);
			    
	st = st_info_acrx_intfc(layer_label,region_label,surf_label,el,
				coords,nor,ahead_side,ahead_st,w_type,
				st,layer_sys,ip,init);
	ct->extra = (POINTER) st;
	if (is_obstacle_state(st))
	{
	    screen("\nERROR in set_extra_ambient(), Obstacle state!\n");
	    clean_up(ERROR);
	}
	set_state(st, GAS_STATE, st);

	free(ahead_st);
}		/*end set_extra_ambient*/

LOCAL	void	set_extra_random_region(
	COMP_TYPE	*ct,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_PHYSICS    *ip,
	INIT_DATA	*init)
{
	Locstate	st = Mean(Random_state(ct));
	SIDE            ahead_side;
	int		surf_label, el;
	int		w_type;
	double		coords[MAXD], nor[MAXD];
	Locstate	ahead_st;
	Front		*front = layer_sys->front;

	alloc_state(front->interf,&ahead_st,front->sizest);
			   
	w_type = get_ahead_state_coords_etc(layer_sys,layer_label,region_label,
				            coords,nor,&ahead_side,ahead_st,
					    &surf_label,
					    &el,NO,init);
			    
	st = st_info_acrx_intfc(layer_label,region_label,surf_label,el,
				  coords,nor,ahead_side,
				  ahead_st,w_type,st,layer_sys,ip,init);
	if (is_obstacle_state(st))
	{
		screen("\nERROR in set_extra_random_region(), ");
		screen("Obstacle state!\n");
		clean_up(ERROR);
	}
	init_random_state_region(st,Random_state(ct),front);
	free(ahead_st);
}		/*end set_extra_random_region*/


LOCAL	void	set_extra_rt_kh(
	COMP_TYPE	*ct,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_PHYSICS    *ip,
	INIT_DATA	*init)
{
	Front		*front = layer_sys->front;
	SIDE            ahead_side;
	int		dim = front->rect_grid->dim;
	int		j, surf_label, el;
	int		w_type;
	char		name[150];
	double		coords[MAXD], nor[MAXD];
	Locstate	ahead_st;
	_RT_KH		*extra = Rt_kh(ct);

	alloc_state(front->interf,&ahead_st,front->sizest);
			   
	w_type = get_ahead_state_coords_etc(layer_sys,layer_label,region_label,
		   	                    coords,nor,&ahead_side,ahead_st,
					    &surf_label,&el,YES,init);
			    
	region_name(name,layer_label,region_label,dim,
		    layer_sys->layer[layer_label]->num_ellips);

	extra->stratification_type = prompt_for_stratification(name);
	for (j = 0; j < dim; ++j)
	    extra->ref_coords[j] = coords[j];

	extra->ref_state = st_info_acrx_intfc(layer_label,region_label,
					      surf_label,el,coords,nor,
					      ahead_side,ahead_st,w_type,
					      extra->ref_state,layer_sys,
					      ip,init);
	free(ahead_st);
}		/*end set_extra_rt_kh*/

LOCAL	void	set_extra_rarefaction_wave_1d(
	COMP_TYPE	*ct,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_PHYSICS    *ip,
	INIT_DATA	*init)
{
	LAYER		*lyr = layer_sys->layer[layer_label];
	Front		*front = layer_sys->front;
	int		dim = front->rect_grid->dim;
	int		surf_label, el, sign, wv_type;
	int		w_type;
	double		coords[MAXD], nor[MAXD];
	double		zl, zt, spl, spt;
	LAYER_SURF 	*ls = lyr->lower_surf, *us = lyr->upper_surf;
	SIDE            ahead_side;
	Locstate	ahead_st, st0;
	_RAREFACTION_WAVE_1D	*rw1d = Rarefaction_wave_1d(ct);

	alloc_state(front->interf,&ahead_st,front->sizest);
			   
	w_type = get_ahead_state_coords_etc(layer_sys,layer_label,region_label,
			   	            coords,nor,&ahead_side,ahead_st,
					    &surf_label,&el,YES,init);

	wv_type = get_wv_type(surf_label,el,lyr);
			    
	/* checking for consistency */
	if (ls->nor[dim-1] != us->nor[dim-1])
	{
	    screen("ERROR in set_extra_rarefaction_wave_1d(), "
	           "Inconsistent normal directions.\n");
	    clean_up(ERROR);
	}

	/* The information on nor is not needed, we have implicitly assumed 
	 * that nor[dim-1] is in the z direction.
	 */

	if (!((is_rarefaction_leading_edge(ls->wv_type) &&
		 is_rarefaction_trailing_edge(us->wv_type))
		 		 ||
		 (is_rarefaction_leading_edge(us->wv_type) &&
		  is_rarefaction_trailing_edge(ls->wv_type))))
	{
	    screen("ERROR in set_extra_rarefaction_wave_1d(), "
	           "Inconsistent 1d rarefaction wave.\n");
	    clean_up(ERROR);
	}

	st0 = st_info_acrx_intfc(layer_label,region_label,surf_label,el,
				 coords,nor,ahead_side,ahead_st,w_type,NULL,
				 layer_sys,ip,init);
	if (is_rarefaction_leading_edge(ls->wv_type))
	{
	    zl = ls->pbar[dim-1];
	    zt = us->pbar[dim-1];
	    rw1d->l_or_r =  LEFT_FAMILY;
	    sign = 1;
	    if (ls->wv_type == wv_type)
	    {
	    	copy_state(rw1d->stl,st0);
	    	prompt_for_trailing_edge_state(rw1d->stt,st0,sign);
	    }
	    else
	    {
	    	prompt_for_leading_edge_state(rw1d->stl,st0,sign);
	    	copy_state(rw1d->stt,st0);
	    }
	    rw1d->lead = ls;
	    rw1d->trail = us;
	}
	else
	{
	    zl = us->pbar[dim-1];
	    zt = ls->pbar[dim-1];
	    rw1d->l_or_r = RIGHT_FAMILY;
	    sign = -1;
	    if (us->wv_type == wv_type)
	    {
	    	copy_state(rw1d->stl,st0);
	    	prompt_for_trailing_edge_state(rw1d->stt,st0,sign);
	    }
	    else
	    {
	    	prompt_for_leading_edge_state(rw1d->stl,st0,sign);
	    	copy_state(rw1d->stt,st0);
	    }
	    rw1d->lead = us;
	    rw1d->trail = ls;
	}

	rw1d->spl = spl = vel(dim-1,rw1d->stl) - sign*sound_speed(rw1d->stl);
	rw1d->spt = spt = vel(dim-1,rw1d->stt) - sign*sound_speed(rw1d->stt);

	rw1d->zl = zl;
	rw1d->zt = zt;
	if (zl < zt)
	{
	    rw1d->zmin = zl;
	    rw1d->stmin = rw1d->stl;
	    rw1d->zmax = zt;
	    rw1d->stmax = rw1d->stt;
	}
	else
	{
	    rw1d->zmin = zt;
	    rw1d->stmin = rw1d->stt;
	    rw1d->zmax = zl;
	    rw1d->stmax = rw1d->stl;
	}
	rw1d->tbar = (zl - zt)/(spl - spt);
	rw1d->zbar = 0.5*(zl + zt - (spl + spt)*rw1d->tbar);

	free(st0);
	free(ahead_st);
}		/*end set_extra_rarefaction_wave_1d*/


LOCAL	void	prompt_for_trailing_edge_state(
	Locstate	st,
	Locstate	stl,
	int		sign)
{
	int		dim = Params(stl)->dim;
	int		i;
	double		ptmp;

	for (i = 0; i < dim-1; ++i)
	    Vel(st)[i] = Vel(stl)[i];
	screen("Enter the pressure at the trailing edge\n\t"
	       "of the rarefaction wave (< %g): ", Press(stl));
	(void) Scanf("%f\n",&ptmp);
	check_float_input("for the pressure",ptmp,0.0,Press(stl),GE_AND_LE);
	state_on_adiabat_with_pr(stl,ptmp,st,TGAS_STATE);
	Vel(st)[dim-1] = Vel(stl)[dim-1] - sign*riemann_wave_curve(stl,ptmp);
}		/*end prompt_for_trailing_edge_state*/


LOCAL	void	prompt_for_leading_edge_state(
	Locstate	st,
	Locstate	stt,
	int		sign)
{
	int		dim = Params(stt)->dim;
	int		i;
	double		ptmp;

	for (i = 0; i < dim-1; ++i)
	    Vel(st)[i] = Vel(stt)[i];
	screen("Enter the pressure at the leading edge\n\t"
	       "of the rarefaction wave (> %g): ",Press(stt));
	(void) Scanf("%f\n",&ptmp);
	check_float_input("for the pressure",ptmp,Press(stt),0.0,GE_);
	state_on_adiabat_with_pr(stt,ptmp,st,TGAS_STATE);
	Vel(st)[dim-1] = Vel(stt)[dim-1] - sign*riemann_wave_curve(stt,ptmp);
}		/*end prompt_for_leading_edge_state*/


LOCAL	void	set_extra_trans_layer(
	COMP_TYPE	*ct,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_PHYSICS    *ip,
	INIT_DATA	*init)
{
	LAYER		*lyr = layer_sys->layer[layer_label];
	Front		*front = layer_sys->front;
	int		dim = front->rect_grid->dim;
	int		surf_label, el, is_upper_surf;
	int		w_type;
	char		mesg[120], name[150];
	double		coords[MAXD], nor[MAXD]; 
	Locstate	ahead_st, st0, st1;
	LAYER_SURF	*surf;
	SIDE            ahead_side;
	_TRANS_LAYER	*t_l = Trans_layer(ct);

	alloc_state(front->interf,&ahead_st,front->sizest);
			   
	w_type = get_ahead_state_coords_etc(layer_sys,layer_label,region_label,
				            coords,nor,&ahead_side,ahead_st,
					    &surf_label,&el,NO,init);
	is_upper_surf = (surf_label == layer_label) ? 1 : 0;

	if (is_upper_surf)
	{
		st0 = t_l->upper_st;
		st1 = t_l->lower_st;
	}
	else
	{
		st1 = t_l->upper_st;
		st0 = t_l->lower_st;
	}
	(void) st_info_acrx_intfc(layer_label,region_label,surf_label,
				 el,coords,nor,ahead_side,ahead_st,
				 w_type,st0,layer_sys,ip,init);

	t_l->lower_surf = lyr->lower_surf;
	t_l->upper_surf = lyr->upper_surf;

	region_name(name,layer_label,region_label,
		    dim,layer_sys->layer[layer_label]->num_ellips);
	(void) sprintf(mesg,"Prompting for the state %s at point ",name);
	surf = is_upper_surf ? lyr->lower_surf : lyr->upper_surf;
	coords[dim-1] = get_surf_height(coords,surf);
	print_general_vector(mesg,coords,dim,"\n");
	prompt_for_ref_state("",st1,TGAS_STATE,ct->params,init);
	
	free(ahead_st);
}		/*end set_extra_trans_layer*/


LOCAL	void	set_extra_elliptical(
	COMP_TYPE	*ct,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_PHYSICS    *ip,
	INIT_DATA	*init)
{
	LAYER		*lyr = layer_sys->layer[layer_label];
	char		s[Gets_BUF_SIZE];
	Front		*front = layer_sys->front;
	int		dim = front->rect_grid->dim;
	int		i, surf_label, el;
	int		w_type;
	double		coords[MAXD], nor[MAXD], maxr = 0.0;
	Locstate	ahead_st, st0;
	SIDE            ahead_side;
	_ELLIPTICAL	*ellpt;

	alloc_state(front->interf,&ahead_st,front->sizest);
	ellpt = Elliptical(ct);
			   
	w_type = get_ahead_state_coords_etc(layer_sys,layer_label,region_label,
				            coords,nor,&ahead_side,ahead_st,
					    &surf_label,&el,YES,init);
			    
	st0 = st_info_acrx_intfc(layer_label,region_label,surf_label,
				 el,coords,nor,ahead_side,ahead_st,w_type,
				 NULL,layer_sys,ip,init);

	set_state(ellpt->state,TGAS_STATE,st0);
	RadialVelocity(ellpt) = Vel(st0)[dim-1];
	ellpt->ellipsoid = lyr->ellip[el];
	for (i = 0; i < dim; ++i)
	    maxr = max(maxr,ellpt->ellipsoid->rad[i]);
	for (i = 0; i < dim; ++i)
	    ellpt->weight[i] =  ellpt->ellipsoid->rad[i] / maxr;
	
	screen("Type 'y' to initialize random perturbations in this ");
	screen("elliptical region: ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	{
	    ellpt->rstate = allocate_random_state_structure(front->interf);
	    alloc_state(front->interf,ellpt->wkstate,front->sizest);
	    alloc_state(front->interf,ellpt->wkstate+1,front->sizest);
	    init_random_state_region(ellpt->state,ellpt->rstate,front);
	}
	free(ahead_st);
}		/*end set_extra_elliptical*/

/*ARGSUSED*/
LOCAL	void	set_extra_stretching(
	COMP_TYPE	*ct,
	LAYER_SYS	*layer_sys,
	int		layer_label,
	int		region_label,
	INIT_DATA	*init)
{
	LAYER	    *lyr = layer_sys->layer[layer_label];
	Front	    *front = layer_sys->front;
	RECT_GRID   *rgr = front->rect_grid;
	_STRETCHING *str = Stretching(ct);
	double	    *L = rgr->L, *U = rgr->U;
	int	    i, j, k, dim = rgr->dim;
	char	    mesg[256];
	static const char  *dname[] = {"x","y","z"};

	(void) sprintf(mesg," for all regions with component number %d",
		       ct->comp);
	prompt_for_thermodynamics(str->ambient,ct->params,mesg);

	for (i = 0; i < dim-1; ++i)
	{
		str->L[i] = L[i];
		str->U[i] = U[i];
	}
	str->L[dim-1] = lyr->lower_surf->pbar[dim-1];
	str->U[dim-1] = lyr->upper_surf->pbar[dim-1];

	if (region_label > 0) /*Assumes layer containing ellipse is already initialized*/
	{
	    _STRETCHING *lstr = Stretching(comp_type(lyr->comp));
	    for (i = 0; i < 1<<dim; ++i)
	    	for (j = 0; j < dim; ++j)
	    	    str->v[i][j] = lstr->v[i][j];
	}
	for (i = 0; i < 1<<dim; ++i)
	{
	    screen("Enter the velocity on the");
	    for (k = 0; k < dim; ++k)
	    	screen(" %s %s",((i>>k)%2)?"upper":"lower",dname[k]);
	    screen("\n\tcorner of the layer: ");
	    for (j = 0; j < dim; ++j)
	    	(void) Scanf("%f",&str->v[i][j]);
	    (void) Scanf("\n");
	}
}		/*end set_extra_stretching*/

LOCAL	int	get_ahead_state_coords_etc(
	LAYER_SYS	*layer_sys,
	int		ll,
	int		rl,
	double		*coords,
	double		*nor,
	SIDE            *ahead_side,
	Locstate	ahead_st,	/* TGAS state */
	int		*surf_label,
	int		*ellip_label,
	int		ignr_pert,	/* ignore perturbed intfc shape? */
	INIT_DATA	*init)
{
	LAYER		**layer = layer_sys->layer, *lyr = layer[ll];
	LAYER_SURF      *lyr_surf;
	ELLIPSOID	**ellip = lyr->ellip;
	COMP_TYPE	*ct;
	COMPONENT       comp;
	Front		*front = layer_sys->front;
	RECT_GRID	*r_g = front->rect_grid;
	int		num_layers = layer_sys->num_layers;
	int		num_ellips = lyr->num_ellips;
	int		w_type;
	char		name[150];
	int		dim = r_g->dim;

	obstacle_state(front->interf,ahead_st,front->sizest);
	w_type = UNKNOWN_WAVE_TYPE;

	*ellip_label = *surf_label = 0;
	*ahead_side = UNKNOWN_SIDE;
	if (num_ellips == 0)
	{
	    if (ll > 1)
	    {
	        get_info_on_surf(r_g,lyr->lower_surf,coords,nor,ignr_pert);
		*surf_label = ll - 1;
	    }
	    else if (ll < num_layers)
	    {
	        get_info_on_surf(r_g,lyr->upper_surf,coords,nor,ignr_pert);
		*surf_label = ll;
	    }
	    else			/* num_layers == 1 */
		get_info_on_surf(r_g,lyr->lower_surf,coords,nor,ignr_pert);
	}
	else if (rl == 0)
	{
	    region_name(name,ll,1,dim,num_ellips);
	    get_info_on_ellip(r_g,name,ellip[1],coords,nor,init);
	    *ellip_label = 1;
	}

	/* check nearby layers */
	if ((rl == 0) && (ll > 1) &&
	    ((ct = comp_type((comp=layer[ll-1]->comp)))->extra != NULL))
	{
	    lyr_surf = lyr->lower_surf;
	    get_info_on_surf(r_g,lyr_surf,coords,nor,ignr_pert);
	    Get_tgas_state(coords,ahead_st,ct,front->interf,init);
	    w_type = lyr->lower_surf->wv_type;
	    *surf_label = ll - 1;
	    *ahead_side = (comp == lyr_surf->l_comp) ?
	        NEGATIVE_SIDE : POSITIVE_SIDE;
	}
	else if ((rl == 0) && (ll < num_layers) &&
		 ((ct = comp_type((comp=layer[ll+1]->comp)))->extra != NULL))
	{
	    lyr_surf = lyr->upper_surf;
	    get_info_on_surf(r_g,lyr_surf,coords,nor,ignr_pert);
	    Get_tgas_state(coords,ahead_st,ct,front->interf,init);
	    w_type = lyr->upper_surf->wv_type;
	    *surf_label = ll;
	    *ahead_side = (comp == lyr_surf->l_comp) ?
	        NEGATIVE_SIDE : POSITIVE_SIDE;
	}
	/* check nearby layers */

	if (rl > 0) /* assuming prompting from outside in */
	{
	    region_name(name,ll,1,dim,num_ellips);
	    get_info_on_ellip(r_g,name,ellip[rl],coords,nor,init);
	    *ellip_label = rl;
	    if ((ct = comp_type(ellip[rl]->compout))->extra != NULL)
	    {
	    	Get_tgas_state(coords,ahead_st,ct,front->interf,init);
	    	w_type = ellip[rl]->wv_type;
		*ahead_side = (ellip[rl]->nor_orient == POSITIVE_ORIENTATION) ?
		    POSITIVE_SIDE : NEGATIVE_SIDE;
	    }
	}
	return w_type;
}		/*end get_ahead_state_coords_etc*/


LOCAL	void	region_name(
	char		*mesg,
	int		layer_label,
	int		region_label,
	int		dim,
	int		num_ellips)
{
	char	   mesg2[100]; 
	const char *elname = (dim == 3) ? "ellipsoid" : "ellipse";

	if (region_label == 0)
		(void) sprintf(mesg2, "%s", (num_ellips == 0) ?
			       "" : " in the out most region");
	else
		(void) sprintf(mesg2, " immediately inside the %d%s %s", 
			     region_label,ordinal_suffix(region_label),elname);
	(void) sprintf(mesg, "for the material%s of the %d%s layer", 
		       mesg2,layer_label,ordinal_suffix(layer_label));
}		/*end region_name*/


LOCAL	int	get_wv_type(
	int		surf_label,
	int		ellip_label,
	LAYER		*lyr)
{
	int layer_label = lyr->layer_label;
	int wv_type = UNKNOWN_WAVE_TYPE;

	if (ellip_label == 0)
	{
	    if (surf_label == layer_label)
	        wv_type = lyr->upper_surf->wv_type;
	    else if (surf_label == layer_label-1)
	        wv_type = lyr->lower_surf->wv_type;
	}
	else if ((ellip_label > 0) && (ellip_label <= lyr->num_ellips))
	    wv_type = lyr->ellip[ellip_label]->wv_type;
	return	wv_type;
}		/*end get_wv_type*/


LOCAL	Locstate	st_info_acrx_intfc(
	int		layer_label,
	int		region_label,
	int		surf_label,
	int		ellip_label,
	double		*coords,
	double		*nor,
	SIDE            ahead_side,
	Locstate	ahead_st,
	int		wv_type,
	Locstate	st,
	LAYER_SYS	*layer_sys,
	INIT_PHYSICS    *ip,
	INIT_DATA	*init)
{
	LAYER		*lyr = layer_sys->layer[layer_label];
	Front		*front = layer_sys->front;
	COMP_TYPE	*ct;
	COMPONENT	compin;
	int		dim;
	size_t		sizest;
	int		num_ellips = lyr->num_ellips;
	int		num_layers = layer_sys->num_layers;
	char		mesg[120], name[150];

	compin = (region_label > 0) ? lyr->ellip[region_label]->compin 
				    : lyr->comp;
	ct = comp_type(compin);
	if (ct->type == OBSTACLE)
	{
	    if (st != NULL)
    		obstacle_state(front->interf,st,front->sizest);
	    else
	    	st = return_obst_state();
	    return st;
	}
	dim = ct->params->dim;
	sizest = ct->params->sizest;

	region_name(name,layer_label,region_label,dim,num_ellips);
	if (wv_type == UNKNOWN_WAVE_TYPE)
	    wv_type = get_wv_type(surf_label,ellip_label,lyr); 

	if (st == NULL)
	    alloc_state(front->interf,&st,sizest);

	if (is_obstacle_state(ahead_st))
	{
	    (void) sprintf(mesg,"Prompting for the state %s at point ",name);
	    print_general_vector(mesg,coords,dim,"\n");
	    prompt_for_ref_state("",st,TGAS_STATE,ct->params,init);
	}
	else if (wv_type == RIEMANN_PROBLEM_WAVE)
	{
	    (void) sprintf(mesg,"Prompting for the state %s at point ",name);
	    print_general_vector(mesg,coords,dim,"\n");
	    prompt_for_ref_state("",st,TGAS_STATE,ct->params,init);
	    set_up_riemann_problem_region(layer_label,region_label,surf_label,
	                                  ellip_label,coords,nor,ahead_side,
					  ahead_st,st,layer_sys,ip,init);
	}
	else if (is_scalar_wave(wv_type))
	{
	    if (((surf_label == layer_label) && (layer_label < num_layers) &&
		 (comp_type(layer_sys->layer[layer_label+1]->comp)->type 
							== TRANS_LAYER))
		    	 ||
		((surf_label == layer_label-1) && (layer_label > 1) &&
		 (comp_type(layer_sys->layer[layer_label-1]->comp)->type 
							== TRANS_LAYER))
			 ||
		(ct->type == TRANS_LAYER))
		{
		    set_state(st,TGAS_STATE,ahead_st);
		    Init_params(st,ct->params);
		}
		else
		    prompt_for_behind_contact_state(ahead_st,
					st,ct->params,YES,nor,TGAS_STATE,init);
	}
	else if (is_shock_wave(wv_type))
	{
	    boolean isforward = is_forward_wave(wv_type) ? YES : NO;
	    prompt_for_behind_shock_state(ahead_st,st,YES,nor,TGAS_STATE,
					  isforward,init);
	}
	else if (is_rarefaction_wave(wv_type))
	{
	    set_state(st,TGAS_STATE,ahead_st);
	}
	else
	{
	    screen("ERROR in st_info_acrx_intfc(), "
	           "Unsuported wave type (%d) %s\n",wv_type,
		   wave_type_as_string(wv_type,front->interf));
	    clean_up(ERROR);
	}
	prompt_for_constant_flow_region(ct->comp,st,front->interf,init);
	return	st;		/* TGAS_STATE */
}		/*end st_info_acrx_intfc*/

EXPORT	void g_prompt_for_constant_flow_region(
	COMPONENT  comp,
	Locstate   st,
	INTERFACE  *intfc)
{
	char	s[Gets_BUF_SIZE];
	
	screen("Is the flow in this region constant in time? "
	       "(default = no): ");
	(void) Gets(s);
	if ((s[0] == 'y') || (s[0] == 'Y'))
	{
	    int stype = state_type(st);
	    set_state(st,GAS_STATE,st);
	    (void)SetConstantFlowRegion(comp,st,intfc);
	    set_state(st,stype,st);
	}
}	/* end g_prompt_for_constant_flow_region */ 

LOCAL	void 	get_info_on_surf(
	RECT_GRID	*rect_grid,
	LAYER_SURF	*surf,
	double		*coords,
	double		*nor,
	int		ignr_pert)
{
	int		i, dim = rect_grid->dim;

	for (i = 0; i < dim; ++i)
	    coords[i] = surf->pbar[i];
	if ((ignr_pert == NO) && (surf->fpoly != NULL))
	    coords[dim-1] = get_surf_height(coords,surf);

	for (i = 0; i < dim; ++i)
	    nor[i] = surf->nor[i];
	if ((ignr_pert == NO) && (surf->fpoly != NULL))
	{
	    for (i = 0; i < dim; ++i)
	          nor[i] += 0.0; /* well you figure out what to add :-) */
	}
}		/*end get_info_on_surf*/


LOCAL	void 	get_info_on_ellip(
	RECT_GRID	*rect_grid,
	char		*mesg,
	ELLIPSOID	*ellip,
	double		*coords,
	double		*nor,
	INIT_DATA	*init)
{
	int		i, j, k, dim = rect_grid->dim;
	double		g[MAXD], temp = 0.0;
	double		M[MAXD][MAXD];
	double		**Q = ellip->Q;
	double		alpha;
	char		choice[Gets_BUF_SIZE];

	/* find the ray direction */
	if (is_gravity() == YES)
	{
	    static	int	sign = 0;
	    const double	*grav = gravity(ellip->cen,initial_time(init));

	    for (i = 0; i < dim; ++i)
	    	g[i] = grav[i];
	    temp = mag_vector(g,dim);
	    for (i = 0; i < dim; ++i)
		g[i] /= temp;
	    if (is_scalar_wave(ellip->wv_type) && (sign == 0))
	    {
	    	screen("\n%s\n"
	    	       "is it heavier or lighter than "
	    	       "the surounding material?\n"
	    	       "\theavier (h)\n"
	    	       "\tlighter (l, dlft)\n"
	    	       "Enter choice: ",mesg);
	    	(void) Gets(choice);
	    	sign = ((choice[0] == 'h') || (choice[0] == 'H')) ? 1 : -1;
	    }
	    for (i = 0; i < dim; ++i)
	    	g[i] *= sign;
	}
	else
	{
	    for (i = 0; i < dim; ++i)
	    	g[i] = 0.0;
	    g[dim-1] = 1.0;
	}

	/* Coords is the intersection of the ray g and the ellipsoid.
	 * Nor is the normal of the ellipsoid at coords.
	 */

	for (i = 0; i < dim; ++i)
	for (j = 0; j < dim; ++j)
	    for (M[i][j] = 0.0, k = 0; k < dim; ++k)
	    	M[i][j] += Q[i][k]*Q[j][k]/sqr(ellip->rad[k]);

	for (alpha = 0, i = 0; i < dim; ++i)
	for (j = 0; j < dim; ++j)
	    alpha += g[i]*M[i][j]*g[j];
	alpha = 1.0/sqrt(alpha);

	for (i = 0; i < dim; ++i)
	    coords[i] = ellip->cen[i] + alpha*g[i];

	for (i = 0; i < dim; ++i)
	    for (nor[i] = 0.0, j = 0; j < dim; ++j)
	    	nor[i] += M[i][j]*g[j];
	temp = mag_vector(nor,dim);
	for (i = 0; i < dim; ++i)
	    nor[i] /= temp;

	/* be sure nor[] is consistent with the nor_orient field in ellip */
	if (ellip->nor_orient == NEGATIVE_ORIENTATION)
	{
	    for (i = 0; i < dim; ++i)
	    	nor[i] *= -1.0;
	}
}		/*end get_info_on_ellip*/


LOCAL	void	goto_moving_frame(
	LAYER_SYS	*layer_sys)
{
	char		choice[Gets_BUF_SIZE];
	double		speed;

	speed = moving_frame_speed(layer_sys);
	screen("\nThe default option for the moving frame velocity is such "
	       "that the average speed\n of the contact surfaces is zero "
	       "after the shock interactions.\n"
	       "Enter the vertical velocity of the moving frame"
	       " (hit return, dflt): ");
	(void) Gets(choice);
	if (choice[0] != '\0')
	   (void) sscan_float(choice,&speed);

	if (speed != 0.0)
	    change_ref_frame(layer_sys, speed);
}		/*end goto_moving_frame*/

LOCAL	double	du;

LOCAL	double	moving_frame_speed(
	LAYER_SYS	*layer_sys)
{
	COMP_TYPE	*ct;
	LAYER		**layer = layer_sys->layer;
	LAYER_SURF	*surf;
	Locstate	rst, lst, mst;
	double		tot_speed = 0;
	double		pstarl = 0.0, pstarr = 0.0;
	double		ustarl = 0.0, ustarr = 0.0;
	double		ml, mr;
	double		nor[MAXD];
	double		rho_rstar, rho_lstar;
	double		ua = 0.0;
	double		n;
	int		wv_type;
	int		dim = layer_sys->front->rect_grid->dim;
	int     	num_layers = layer_sys->num_layers;
	int             i, j, k, si = 0, num_contact, up = YES;
	RIEMANN_SOLVER_WAVE_TYPE l_wv, r_wv;
	size_t		sizest = layer_sys->front->sizest;

	(void) printf("\nOne dimensional data for "
		      "Richtmyer-Meshkov interaction\n");

	/* find out the shock interface and its direction */
	for (i = 1; i < num_layers; ++i)
	{
	    surf = layer[i]->upper_surf;
	    wv_type = surf->wv_type;
	    n = surf->nor[dim-1];
	    if (((wv_type == FORWARD_SHOCK_WAVE)  && (n > 0)) ||
	        ((wv_type == BACKWARD_SHOCK_WAVE) && (n < 0)))
	    {
	    	up = YES;
	    	si = i;
	    	break;
	    }
	    else if (((wv_type == FORWARD_SHOCK_WAVE)  && (n < 0)) ||
		     ((wv_type == BACKWARD_SHOCK_WAVE) && (n > 0)))
	    {
	    	up = NO;
	    	si = i;
	    	break;
	    }
	}
	for (i = 1; i < dim; ++i)
	    nor[i] = 0.0;
	nor[0] = (up == YES) ? 1.0 : -1.0;
	
	/* assuming the right state is always ahead */
	alloc_state(layer_sys->front->interf,&rst, sizest);
	alloc_state(layer_sys->front->interf,&lst, sizest);
	alloc_state(layer_sys->front->interf,&mst, sizest);
	i = si;
	while ((i >= 1) && (i < num_layers))
	{
	    surf = (up == YES) ? layer[i+1]->upper_surf : layer[i]->lower_surf; 
	    if (!is_scalar_wave(surf->wv_type))
	    	break;

	    /* find right (ahead) state */
	    j = (up == YES) ? i+2 : i-1;
	    ct = comp_type(layer[j]->comp);
	    if (ct->type == AMBIENT)
	    	set_state(rst, TGAS_STATE, Ambient(ct));
	    else if (ct->type == RANDOM_REGION)
	    	set_state(rst, TGAS_STATE,Mean(Random_state(ct)));
	    else
	    {
	    	(void) printf("WARNING: in moving_frame_speed(), "
	    	              "Ahead state is not ambient!\n");
	    	break;
	    }
	    Vel(rst)[0] = nor[0]*Vel(rst)[dim-1];
	    for (k = 1; k < dim; ++k) Vel(rst)[k] = 0.0;
	    ua += Vel(rst)[0];

	    /* find left (behind) state */
	    if ((pstarl == 0.0) && (ustarl == 0.0))
	    {
	        j = (up == YES) ? i : i+1;
	        ct = comp_type(layer[j]->comp);
	        if (ct->type == AMBIENT)
	            set_state(lst,TGAS_STATE,Ambient(ct));
	        else if (ct->type == RANDOM_REGION)
	            set_state(lst,TGAS_STATE,Mean(Random_state(ct)));
	        else
	        {
	            (void) printf("WARNING: in moving_frame_speed(), "
	                          "Behind state is not ambient!\n");
	            break;
	        }
	        verbose_print_state("Behind shock state",lst);

	        j = (up == YES) ? i+1 : i;
	        ct = comp_type(layer[j]->comp);
	        if (ct->type == AMBIENT)
	            set_state(mst,TGAS_STATE,Ambient(ct));
	        else if (ct->type == RANDOM_REGION)
	            set_state(mst,TGAS_STATE,Mean(Random_state(ct)));
	        else
	        {
	            (void) printf("WARNING: in moving_frame_speed(), "
		    	          "Ahead state is not ambient!\n");
		    break;
		}
		verbose_print_state("Ahead shock state",mst);
		ua += nor[0]*Vel(Ambient(ct))[dim-1];
	    }
	    else
	    {
	        j = (up == YES ? i+1 : i);
	        ct = comp_type(layer[j]->comp);
	        (void) s_polar_4(BEHIND_PRESSURE,pstarr,&ustarr,nor,
	        		 Ambient(ct),lst,TGAS_STATE);
	    }
	    Vel(lst)[0] = nor[0]*Vel(lst)[dim-1];
	    for (k = 1; k < dim; ++k)
	        Vel(lst)[k] = 0.0;
	    verbose_print_state("Ahead contact state",rst);

	    (void) find_mid_state(lst,rst,0.0,&pstarl,&pstarr,
	    		          &ustarl,&ustarr,&ml,&mr,&l_wv,&r_wv);

	    (void) printf("\nPreshock Atwood number = %g\n",
	    	          (Dens(mst)-Dens(rst))/(Dens(mst)+Dens(rst)));

	    rho_rstar = dens_Hugoniot(pstarr,rst);
	    rho_lstar = (pstarl < pressure(lst)) ?
		dens_rarefaction(pstarl,lst) : dens_Hugoniot(pstarl,lst);
	    (void) printf("Postshock Atwood number = %g\n",
			  (rho_lstar-rho_rstar)/(rho_lstar+rho_rstar));
	    (void) printf("Postshock pressure = %g\n",pstarr);

	    ustarr *= nor[0];
	    ustarl *= nor[0];

	    if (r_wv != SHOCK)
	    {
	        screen("ERROR in moving_frame_speed(), "
	               "How can a right wave not be a shock here!\n");
	        clean_up(ERROR);
	    }
	    tot_speed += 0.5*(ustarr + ustarl);
	    up == YES ? i++ : i--;
	}
	num_contact = (up == YES) ? i - si : si - i;
	ua /= num_contact + 1;
	if (num_contact > 0)
	{
	    tot_speed /= num_contact;
	    du = fabs(ua - tot_speed);
	}
	else
	{
	    tot_speed = 0.0;
	    du = 0.0;
	}
	(void) printf("Translation speed to 1D steady contact = %g\n",
		      tot_speed);
	(void) printf("End One dimensional data for ");
	(void) printf("Richtmyer-Meshkov interaction\n");

	free_these(3,lst,rst,mst);
	return	tot_speed;
}		/*end moving_frame_speed*/

LOCAL	void	print_rm_scale_factors(
	LAYER_SYS	*layer_sys)
{
	LAYER		**layer = layer_sys->layer;
	LAYER_SURF	*surf;
	double	a = 0.0;
	int	i;
	int     num_layers = layer_sys->num_layers;
	int	num_contact = 0;

	if (du == 0.0)
	    return;

	for (i = 1; i < num_layers; ++i)
	{
	    surf = layer[i]->upper_surf;
	    if (! is_scalar_wave(surf->wv_type))
	    	continue;
	    ++num_contact;
	    a += fabs(surf->s_max - surf->s_min);
	}
	if (num_contact == 0)
	    return;
	
	a /= 2*num_contact;
	(void) printf("\nAverage change in interface velocity, du = %g\n",du);
	(void) printf("Average preshocked interface amplitude, a = %g\n",a);
	(void) printf("Time scale, a/du = %g\n\n",a/du);
}		/*end print_rm_scale_factors*/

LOCAL	void	change_ref_frame(
	LAYER_SYS	*layer_sys,
	double		speed)
{
	int		dim = layer_sys->front->rect_grid->dim;
	int		i, j;
	int		l, m, n, gmax[3];
	double		v[3];
	LAYER		*lyr;
	COMP_TYPE	*ct;
	COMPONENT	compin;
	RANDOM_STATE	*rstate;
	Locstate	st;

	for (i = 0; i < dim; ++i)
		v[i] = 0.0;
	v[dim-1] = -speed;
	for (i = 1; i <= layer_sys->num_layers; ++i)
	{
	    lyr = layer_sys->layer[i];
	    for (j = 0; j <= lyr->num_ellips; ++j)
	    {
	        compin = (j == 0) ? lyr->comp : lyr->ellip[j]->compin;
	        switch ((ct = comp_type(compin))->type)
	        {
	        case	AMBIENT:
	            st = Ambient(ct);
		    add_velocity_to_state(st,v);
	            break;
	        case	RANDOM_REGION:
		    rstate = Random_state(ct);
	            st = Mean(rstate);
		    add_velocity_to_state(st,v);
		    for (l = 0; l < 3; ++l)
			gmax[l] = 1;
		    for (l = 0; l < dim; ++l)
			gmax[l] = rstate->grid.gmax[l]+1;
		    for (l = 0; l < gmax[0]; ++l)
		    for (m = 0; m < gmax[1]; ++m)
		    for (n = 0; n < gmax[2]; ++n)
			add_velocity_to_state(rstate->old_st[l][m][n],v);
	            break;
	        case	RAREFACTION_WAVE_1D:
	            st = Rarefaction_wave_1d(ct)->stl;
		    add_velocity_to_state(st,v);
	            st = Rarefaction_wave_1d(ct)->stt;
		    add_velocity_to_state(st,v);
	            break;
	        case	TRANS_LAYER:
	            st = Trans_layer(ct)->lower_st;
		    add_velocity_to_state(st,v);
	            st = Trans_layer(ct)->upper_st;
		    add_velocity_to_state(st,v);
	            break;
	        /* add whatever you need here */
	        default:
	            screen("ERROR in change_ref_frame(), ");
	            screen("Unsupported comp_type->type.\n");
	            clean_up(ERROR);
	        }
	    }
	}
}		/*end change_ref_frame*/


/*---------------------------------------------------------------------------*/
/*		
*	Functions ft_assigning a state for a given coordinate are listed below.
*	All states returned by get_state_whatever() functions are in GAS_STATE.
*/

EXPORT	boolean	rarefaction_edge_at_coords(
	double                 *coords,
	HYPER_SURF            *hs,
	RAREFACTION_EDGE_TYPE type)
{
	int wv_type;
	if (hs == NULL)
	    return NO;
	wv_type = wave_type(hs);
	if ((type == LEADING_EDGE) && is_rarefaction_leading_edge(wv_type))
	    return YES;
	if ((type == TRAILING_EDGE) && is_rarefaction_trailing_edge(wv_type))
	    return YES;
	if (is_rarefaction_wave(wv_type))
	    return NO;
	if (hs->interface->dim == 2)
	{
	    CURVE **c, *cur = Curve_of_hs(hs);
	    NODE  *n = NULL;

	    if (coords == Coords(cur->start->posn))
		n = cur->start;
	    else if (coords == Coords(cur->end->posn))
		n = cur->end;
	    else
		return NO;
	    for (c = n->in_curves; c && *c; ++c)
	    {
		wv_type = wave_type(*c);
	        if ((type == LEADING_EDGE) &&
		    is_rarefaction_leading_edge(wv_type))
	            return YES;
	        if ((type == TRAILING_EDGE) &&
		    is_rarefaction_trailing_edge(wv_type))
	            return YES;
	    }
	    for (c = n->out_curves; c && *c; ++c)
	    {
		wv_type = wave_type(*c);
	        if ((type == LEADING_EDGE) &&
		    is_rarefaction_leading_edge(wv_type))
	            return YES;
	        if ((type == TRAILING_EDGE) &&
		    is_rarefaction_trailing_edge(wv_type))
	            return YES;
	    }
	}
	return NO;
}		/*end rarefaction_edge_at_coords*/

EXPORT	double	get_surf_height(
	double		*coords,
	LAYER_SURF	*surf)
{
	double	height;

	if (surf->fpoly != NULL)
	    height = fourier_poly(coords,surf->fpoly);
	else
	{
	    double *pbar = surf->pbar, *nor = surf->nor;
	    int	  i, k = surf->dim-1;;

	    height = surf->pbar[k];;
	    for (i = 0; i < k; ++i)
	    	height -= (nor[i]/nor[k])*(coords[i] - pbar[i]);
	}

	surf->s_min = min(height,surf->s_min);
	surf->s_max = max(height,surf->s_max);
	return height;
}		/*end get_surf_height*/


/*ARGSUSED*/
LOCAL	void	set_extra_1d_overlay(
	INIT_DATA	*init,
	COMP_TYPE	*ctin,
	INIT_PHYSICS	*ip,
	LAYER_SYS	*layer_sys)
{
	MOVING_FRAME(layer_sys->flag) = NO;
	prompt_for_1d_overlay(init,ctin,ip,NULL,OVERLAY_TYPE_UNSET);
}		/*end set_extra_1d_overlay*/

EXPORT	void	prompt_for_1d_overlay(
	INIT_DATA	*init,
	COMP_TYPE	*ctin,
	INIT_PHYSICS	*ip,
	const char      *mesg,
	OVERLAY_TYPE    overlay_type)
{
	IO_TYPE         IO_type;
	COMPONENT	compin;
	COMPONENT	comp;
	COMP_TYPE	*ct;
	FILE		*file = NULL;
	Front		*front = ip->root->front;
	INTERFACE	*cur_intfc;
	ONED_OVERLAY	*olay = One_d_overlay(ctin);
	RECT_GRID	*gr = ip->root->front->rect_grid;
	char		s[Gets_BUF_SIZE];
	double		magv;
	int		nvars;
	int		i;
	int		step;
	int		dim = gr->dim;
	int		grids_set;

	compin = ctin->comp;

	if ((mesg != NULL) && (strlen(mesg) != 0))
	{
	    screen("Initializing a one dimensional overlay for the region\n"
		   "\t%s\n",mesg);
	}
	if (overlays_exist() == YES)
	{
	    screen("Enter the component number of an existing region "
		   "whose one dimensional overlay\n\t"
	           "you wish to use,  or enter the current region's "
		   "component number (%d)\n\t"
	           "or RETURN to initiate a new overlay structure.\n",compin);
	    screen("Enter component number (default creates a new overlay): ");
	    (void) Gets(s);
	    if (s[0] != '0')
	    {
	    	(void) sscanf(s,"%d",&comp);
	    	if ((comp<FIRST_DYNAMIC_COMPONENT) || (comp>max_num_comps()))
		{
		    screen("ERROR in set_extra_1d_overlay(), "
		           "invalid component %d\n",comp);
		    clean_up(ERROR);
		}
		ct = comp_type(comp);
		if (ct->type != ONE_DIMENSIONAL_OVERLAY)
		{
		    screen("ERROR in set_extra_1d_overlay(), component %d "
			   "is not a 1D overlay region\n",comp);
		    clean_up(ERROR);
		}
		free_1d_overlay_comp_type(ctin);
		ctin->extra = ct->extra;
		ctin->free_comp_type_extra = NULL;
		return;
	    }
	}
	screen("Enter the file name of the overlay data: ");
	(void) Gets(s);
	if ((s[0] == '\0') || ((file = fopen(s,"r")) == NULL))
	{
	    screen("ERROR in set_extra_1d_overlay(), input file ");
	    if (s[0] == '\0')
	    	screen("not supplied\n");
	    else
	    	screen("%s does not exist or is unreadable\n",s);
	    clean_up(ERROR);
	}
	determine_io_type(file,&IO_type);
	screen("Enter the time step at which to find the overlay data: ");
	(void) Scanf("%d\n",&step);
	cur_intfc = current_interface();
	if (position_file_at_read_time(file,step,NULL,NULL) == FUNCTION_FAILED)
	{
	    screen("ERROR in set_extra_1d_overlay(), input file "
	           "does not contain time step data\n");
	    clean_up(ERROR);
	}
	olay->prt = ip->prt;
	olay->front = front;
	nvars = olay->prt->n_restart_vars;
	olay->is = g_set_input_solution(olay->prt);
	set_interface_hooks(1,init);
	set_size_of_intfc_state(StateSize(init));
	if ((olay->intfc1d=read_print_interface(init,&IO_type,YES,&grids_set))==NULL)
	{
	    screen("ERROR in set_extra_1d_overlay(), "
	           "can not read interface\n");
	    clean_up(ERROR);
	}
	if (read_state_variables(&IO_type,nvars,olay->intfc1d,olay->is,1) ==
							FUNCTION_FAILED)
	{
	    screen("ERROR in set_extra_1d_overlay(), "
	           "can not read state variables\n");
	    clean_up(ERROR);
	}
	(void) Fclose(file);
	set_current_interface(cur_intfc);
	if (overlay_type != OVERLAY_TYPE_UNSET)
	    olay->overlay_type = overlay_type;
	else
	{
	    olay->overlay_type = RECTANGULAR_OVERLAY;
	    screen("Enter the remap geometry of the one dimensional data\n\t");
	    screen("Choices are\n"
		   "\tRECTANGULAR (default)\n"
		   "\tRADIAL\n"
		   "\tCYLINDRICAL\n"
		   "Enter choice: ");
	(void) Gets(s);
	    if (strcasecmp(s,"RADIAL") == 0)
	        olay->overlay_type = RADIAL_OVERLAY;
	    else if (strcasecmp(s,"RECTANGULAR") == 0)
	        olay->overlay_type = RECTANGULAR_OVERLAY;
	    else if (strcasecmp(s,"CYLINDRICAL") == 0)
	        olay->overlay_type = CYLINDRICAL_OVERLAY;
	}

	olay->origin[0] = olay->origin[1] = olay->origin[2] = 0.0;
	switch (olay->overlay_type)
	{
	case RADIAL_OVERLAY:
	    olay->direction[0] = 0.0;
	    olay->direction[1] = 0.0;
	    olay->direction[2] = 0.0;
	    break;
	case CYLINDRICAL_OVERLAY:
	    olay->direction[0] = 0.0;
	    olay->direction[1] = 0.0;
	    olay->direction[2] = 1.0;
	    break;
	case RECTANGULAR_OVERLAY:
	    olay->direction[0] = 1.0;
	    olay->direction[1] = 0.0;
	    olay->direction[2] = 0.0;
	    break;
	default:
	    screen("ERROR in prompt_for_1d_overlay(), invalid "
		   "overlay type %d\n",olay->overlay_type);
	    clean_up(ERROR);
	}
	screen("Enter the coordinates of the origin of "
	       "the one dimensional solution\n\t",
	       olay->origin[0],olay->origin[1],olay->origin[2]);
	if (input_coordinates(olay->origin,3) == FUNCTION_FAILED)
	{
	    screen("ERROR in prompt_for_1d_overlay(), "
	           "improper origin input\n");
	    clean_up(ERROR);
	}
	switch (olay->overlay_type)
	{
	case RADIAL_OVERLAY:
	    break;
	case CYLINDRICAL_OVERLAY:
	    screen("Enter the direction of the cylindrical symmetry axis\n\t");
	    if (input_coordinates(olay->direction,3) == FUNCTION_FAILED)
	    {
		screen("ERROR in prompt_for_1d_overlay(), "
		       "improper direction input\n");
		clean_up(ERROR);
	    }
	    magv = mag_vector(olay->direction,dim);
	    for (i = 0; i < dim; ++i)
	    	olay->direction[i] /= magv;
	    break;
	case RECTANGULAR_OVERLAY:
	    screen("Enter the direction of the one dimensional axis\n\t");
	    if (input_coordinates(olay->direction,3) == FUNCTION_FAILED)
	    {
		screen("ERROR in set_extra_1d_overlay(), "
		       "improper direction input\n");
		clean_up(ERROR);
	    }
	    magv = mag_vector(olay->direction,dim);
	    for (i = 0; i < dim; ++i)
	    	olay->direction[i] /= magv;
	    break;
	case OVERLAY_TYPE_UNSET:
	default:
	    screen("ERROR in prompt_for_1d_overlay(), unset overlay type\n");
	    clean_up(ERROR);
	}
}		/*end prompt_for_1d_overlay*/

LOCAL	boolean	input_coordinates(
	double	*coords,
	int	dim)
{
	char	s[Gets_BUF_SIZE];
	int	i, n;

	screen("(default = %g",coords[0]);
	for (i = 1; i < dim; ++i)
	    screen(", %g",coords[i]);
	screen("): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    double crds[3];
	    static const char *fmt = "%lf %lf %lf";

	    n = sscanf(s,fmt,crds,crds+1,crds+2);
	    if (n != dim)
	    	return FUNCTION_FAILED;
	    for (i = 0; i < dim; ++i)
		coords[i] = crds[i];
	}
	return FUNCTION_SUCCEEDED;
}		/*end input_coordinates*/


EXPORT	void	set_1d_overlay_comp_type(
	COMP_TYPE	*comp_type)
{
	ONED_OVERLAY	*olay;

	if (comp_type->type == ONE_DIMENSIONAL_OVERLAY) /*ALREADY SET*/
		return;

	if (comp_type->free_comp_type_extra != NULL)
	    (*comp_type->free_comp_type_extra)(comp_type);

	comp_type->type = ONE_DIMENSIONAL_OVERLAY;
	scalar(&olay,sizeof(ONED_OVERLAY));
	olay->overlay_type = OVERLAY_TYPE_UNSET;
	comp_type->extra = (POINTER)olay;

	comp_type->_get_state = get_state_1d_overlay;
	comp_type->free_comp_type_extra = free_1d_overlay_comp_type;

}		/*end set_1d_overlay_comp_type*/


LOCAL	void	free_1d_overlay_comp_type(
	COMP_TYPE	*comp_type)
{
	ONED_OVERLAY *olay;
	if (comp_type->type != ONE_DIMENSIONAL_OVERLAY)
		return;

	olay = One_d_overlay(comp_type);
	if (olay == NULL)
		return;

	if (olay->is != NULL)
		free(olay->is);
	free(olay);
	comp_type->extra = NULL;
}		/*end free_1d_overlay_comp_type*/
