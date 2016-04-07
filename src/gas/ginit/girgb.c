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
*				girgb.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <ginit/ginit.h>

 
EXPORT	void 	init_fluid_rigid_body(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	LAYER_FLAG      flag;
	LAYER_SYS       *layer_sys;
	LAYER           **layer;
	Front           *front = ip->root->front;
	RECT_GRID	*gr = front->rect_grid;
	INTERFACE	*intfc = front->interf;
	COMPONENT	left_c,right_c;
	RIGID_BODY_PARAMS rgb_params;
	CURVE		**curves;
	POINTER		func_params;
	double		(*func)(POINTER,double*);
	int 		i, j, k, dim = gr->dim;
	int		num_rigid_bodys;	
	char		mesg[100];
	char		s[Gets_BUF_SIZE];
	int		is_uniform;
	int		num_segs;

	debug_print("rigid_body","Entered init_fluid_rigid_body()\n");

	scalar(&layer_sys,sizeof(LAYER_SYS));

	screen("Enter yes to initialize a shock in the region: ");
	(void) Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    layer_sys->num_layers = 2;
	else
	    layer_sys->num_layers = 1;
	uni_array(&layer,5*layer_sys->num_layers+1,sizeof(LAYER *));
	layer_sys->layer = layer;
	layer_sys->front = front;
	CLEAR_LAYER_FLAG(layer_sys->flag);

	for (i = 1; i <= layer_sys->num_layers; ++i)
	{
	    layer[i] = prompt_for_each_layer(i,layer_sys->num_layers,
	    			&layer_sys->flag,front,layer,ip);
	}
	if (layer_sys->num_layers == 2)
	{
	    layer[1]->upper_surf->r_comp = layer[1]->comp;
	    layer[1]->upper_surf->l_comp = layer[2]->comp;
	    make_layer_surf(front,layer[1]->upper_surf,init);
	    right_c = layer[2]->comp;
	}
	else
	    right_c = layer[1]->comp;
	left_c = COMPOBST;
	set_obstacle_comp_type(comp_type(COMPOBST),front);

	/** Start specification of rigid bodies **/

	screen("Enter number of rigid bodys: ");
	Scanf("%d\n",&num_rigid_bodys);
	screen("Available geometrical shapes of rigid body are\n");
	screen("\tEllipse (E)\n");
	screen("\tTilt Ellipse (N)\n"); /*JD Kim */
        screen("\tTriangle (T)\n");     /*JD Kim */
        screen("\tRectangle (R)\n");    /*JD Kim */
	screen("\tCosmos (C)\n");
	screen("\tTaegeuk (G)\n");	/*JD Kim */
	screen("\tWing (W)\n");		/*JD Kim */
	screen("\tPropeller (P)\n"); 	/*JD Kim */

	for (i = 1; i <= num_rigid_bodys; ++i)
	{
	    screen("Enter the shape of rigid body %d: ");
	    (void) Gets(s);
	    screen("Is the density of rigid body uniform?(1:yes/0:no): ");
	    Scanf("%d \n",&is_uniform);
	    switch (s[0])
	    {
	    case 'e':
	    case 'E':
	        func_params = init_ellipse_params(gr);
		func = ellipse_func;
		if (is_uniform)
		{
		    ELLIP2D_PARAMS *eparams = (ELLIP2D_PARAMS*)func_params;
	    	    screen("Enter the total mass for %d-th rigid body: ",i+1);
	    	    Scanf("%lf \n",&rgb_params.total_mass);
                    rgb_params.center_of_mass[0] = eparams->x0;
		    rgb_params.center_of_mass[1] = eparams->y0;
        	    rgb_params.mom_of_inertial = (rgb_params.total_mass)*
		    		(sqr(eparams->a)+sqr(eparams->b))/5.0;
		    printf("moment of inertial of uniform ellipseis %lf\n",
		    		rgb_params.mom_of_inertial);
		}
		else
		{
	    	    sprintf(mesg,"the %d-th rigid body",i);
	    	    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
		}
	        break;
	    case 'n':   /*JD Kim */
            case 'N':   /*JD Kim */
                func_params = init_ellipse_tilt_params(gr);     /*JD Kim */
                func = ellipse_tilt_func;       /*JD Kim */
		if (is_uniform)
		{
		    ELLIP2D_TILT_PARAMS *etparams = 
		    		(ELLIP2D_TILT_PARAMS*)func_params;
	    	    screen("Enter the total mass for %d-th rigid body: ",i+1);
	    	    Scanf("%lf \n",&rgb_params.total_mass);
                    rgb_params.center_of_mass[0] = etparams->x0;
                    rgb_params.center_of_mass[1] = etparams->y0;
                    rgb_params.mom_of_inertial = (rgb_params.total_mass)*
		    		(sqr(etparams->a)+sqr(etparams->b))/5.0;
		    printf("moment of inertial of uniform ellipse is %lf\n",
		    		rgb_params.mom_of_inertial);
		}
		else
		{
	    	    sprintf(mesg,"the %d-th rigid body",i);
	    	    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
		}
                break;
            case 't':   /*JD Kim */
            case 'T':   /*JD Kim */
                func_params = init_triangle_params(gr); /*JD */
                func = triangle_func;   /*JD */
		if (is_uniform)
		{
		    double a,b,h;
		    TRIANGLE_PARAMS *tparams = (TRIANGLE_PARAMS*)func_params;
	    	    screen("Enter the total mass for %d-th rigid body: ",i+1);
	    	    Scanf("%lf \n",&rgb_params.total_mass);
                    rgb_params.center_of_mass[0] = 1.0/3.0*(tparams->x[0] + 
					tparams->x[1] + tparams->x[2]);
                    rgb_params.center_of_mass[1] = 1.0/3.0*(tparams->y[0] + 
					tparams->y[1] + tparams->y[2]);
		    a = tparams->x[2] - tparams->x[1];
		    b = tparams->x[1] - tparams->x[0];
		    h = tparams->y[2] - tparams->y[1];
                    rgb_params.mom_of_inertial = (b*b*b*h - b*b*h*a + b*h*a*a 
		    			+ b*h*h*h)/36.0;
		    printf("moment of inertial of uniform triangle is %lf\n",
		    			rgb_params.mom_of_inertial);
		}
		else
		{
	    	    sprintf(mesg,"the %d-th rigid body",i);
	    	    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
		}
                break;
            case 'r':  /*JD Kim */
            case 'R':  /*JD Kim */
                func_params = init_rectangle_params(gr);        /*JD*/
                func = rectangle_func;  /*JD */
		if (is_uniform)
		{
		    RECTANGLE_PARAMS *rparams = (RECTANGLE_PARAMS*)func_params;
	    	    screen("Enter the total mass for %d-th rigid body: ",i+1);
	    	    Scanf("%lf \n",&rgb_params.total_mass);
                    rgb_params.center_of_mass[0] = 0.5*(2.0*rparams->x0 + 
					rparams->a);
                    rgb_params.center_of_mass[1] = 0.5*(2.0*rparams->y0 + 
					rparams->b);
                    rgb_params.mom_of_inertial = rgb_params.total_mass*
		    		(sqr(rparams->a)+sqr(rparams->b))/4.0;
		    printf("moment of inertial of uniform rectangle is %lf\n",
		    		rgb_params.mom_of_inertial);
		}
		else
		{
	    	    sprintf(mesg,"the %d-th rigid body",i);
	    	    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
		}
                break;
	    case 'c':
	    case 'C':
		func_params = init_cosmos_params(gr);
		func = cosmos_func;
		if (is_uniform)
		{
		    COSMOS_PARAMS *cparams = (COSMOS_PARAMS*)func_params;
		    screen("Enter the total mass for %d-th rigid body: ",i+1);
		    Scanf("%lf \n",&rgb_params.total_mass);
		    rgb_params.center_of_mass[0] = 1.0;
		    rgb_params.center_of_mass[1] = 2.0;
		    rgb_params.mom_of_inertial = 1.0;
		    printf("moment of inertial of uniform comma is %lf\n",
		    		rgb_params.mom_of_inertial);
		}
		else
		{
		    sprintf(mesg,"the %d-th rigid body",i);
                    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
                }
		break;
	    case 'g':  /*JD Kim */
            case 'G':  /*JD Kim */
                func_params = init_taegeuk_params(gr);
                func = taegeuk_func;
                if (is_uniform)
                {
                    TAEGEUK_PARAMS *tgparams = (TAEGEUK_PARAMS*)func_params;
                    screen("Enter the total mass for %d-th rigid body: ",i+1);
                    Scanf("%lf \n",&rgb_params.total_mass);
                    rgb_params.center_of_mass[0] = tgparams->x2;
                    rgb_params.center_of_mass[1] = tgparams->y;
                    rgb_params.mom_of_inertial = 0.001;
                    printf("moment of inertial of uniform taegeuk is %lf\n",
		    		rgb_params.mom_of_inertial);
                }
                else
                {
                    sprintf(mesg,"the %d-th rigid body",i);
                    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
                }
		break;

	    case 'w':  /*JD Kim */
            case 'W':  /*JD Kim */
                func_params = init_wing_params(gr);
                func = wing_func;
                if (is_uniform)
                {
                    WING_PARAMS *wparams = (WING_PARAMS*)func_params;
                    screen("Enter the total mass for %d-th rigid body: ",i+1);
                    Scanf("%lf \n",&rgb_params.total_mass);
                    rgb_params.center_of_mass[0] = wparams->x1;
                    rgb_params.center_of_mass[1] = wparams->y1;
                    rgb_params.mom_of_inertial = 0.1;
                    printf("moment of inertial of uniform wing is %lf\n",
                                rgb_params.mom_of_inertial);
                }
                else
                {
                    sprintf(mesg,"the %d-th rigid body",i);
                    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
                }
		break;

	    case 'p':  /*JD Kim */
	    case 'P':  /*JD Kim */
	        func_params = init_propeller_params(gr);
		func = propeller_func;
		if (is_uniform)
		{
		    PROPELLER_PARAMS *pparams = (PROPELLER_PARAMS*)func_params;
		    screen("Enter the total mass for %d-th rigid body: ",i+1);
		    Scanf("%lf \n",&rgb_params.total_mass);
		    rgb_params.center_of_mass[0] = pparams->x[0];
		    rgb_params.center_of_mass[1] = pparams->y[0];
		    printf("COM of uniform propeller is %lf %lf\n",
		    	rgb_params.center_of_mass[0],rgb_params.center_of_mass[1]);
		    rgb_params.mom_of_inertial = 0.1;
		    printf("moment of inertial of uniform propeller is %lf\n",
		    		rgb_params.mom_of_inertial);
		}
		else
		{
		    sprintf(mesg,"the %d-th rigid body",i);
		    prompt_for_rigid_body_params(&rgb_params,mesg,dim);
		}
		break;
	    }	
	    screen("Type yes if rigid body will only rotate around COM: ");
	    (void) Gets(s);
	    if (s[0] == 'y' || s[0] == 'Y')
	    {
	    	rgb_params.rotation_only = YES;
		rgb_params.vertical_motion_only = NO;
	    }
	    else
	    {
	    	rgb_params.rotation_only = NO;
		screen("Type yes if you want vertical motion only?: ");
                (void) Gets(s);
                if (s[0] == 'y' || s[0] == 'Y')
                    rgb_params.vertical_motion_only = YES;
                else
                    rgb_params.vertical_motion_only = NO;
	    }

	    curves = make_level_curves(gr,intfc,left_c,right_c,func,
				func_params,NO,&num_segs);
	    if (curves == NULL)
	    {
		screen("Cannot make rigid body %d\n",i);
		clean_up(ERROR);
	    }
	    for (k = 0; k < num_segs; ++k)
	    {
	    	if (curves[k]->start == curves[k]->end)
		    node_type(curves[k]->start) = CLOSED_NODE;
	    	start_status(curves[k]) = end_status(curves[k]) = INCIDENT;
	    	wave_type(curves[k]) = MOVABLE_BODY_BOUNDARY;
	    	total_mass(curves[k]) = rgb_params.total_mass;
	    	mom_inertial(curves[k]) = rgb_params.mom_of_inertial;
	    	angular_velo(curves[k]) = 0.0;
	    	surface_tension(curves[k]) = 0.0;
	    	body_index(curves[k]) = i-1;
	    	if (rgb_params.rotation_only == YES)
	    	    motion_type(curves[k]) = ROTATION;
	    	else if (rgb_params.vertical_motion_only == YES)
	    	    motion_type(curves[k]) = VERTICAL_MOTION;
	    	else
		    motion_type(curves[k]) = FREE_MOTION;
	    	for (j = 0; j < dim; ++j)
	    	{
	    	    center_of_mass(curves[k])[j] = rgb_params.center_of_mass[j];
		    center_of_mass_velo(curves[k])[j] = 0.0;
	    	}
	    }
	    printf("i = %d num_rigid_bodys = %d\n",i,num_rigid_bodys);
	}

	/** End specification of rigid bodies **/

	init_comp_type(init,ip,layer_sys);
	if (debugging("rigid_body"))
	{
	    gview_plot_interface("out",intfc);
	    make_interface_topology_lists(intfc);
	    show_COMP(stdout,intfc);
	}
	debug_print("rigid_body","Left init_fluid_rigid_body()\n");
}	/* end init_fluid_rigid_body */

