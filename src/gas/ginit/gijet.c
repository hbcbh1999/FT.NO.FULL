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
*				gijet.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Initialized an injection jet.
*/

#if defined(TWOD) || defined(THREED)
#define DEBUG_STRING    "init_jet"
#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
#if defined(TWOD)
LOCAL	double eval_poly(double,const double*,int);
LOCAL	void	init_injection_inlet_jet2d(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_fuel_injection_jet2d(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	init_fuel_injection_jet3d(INIT_DATA*,INIT_PHYSICS*);
LOCAL	void	insert_elliptical_curve_section(BOND*,CURVE*,double,double,double*,
	                                        double*,double,Front*);
LOCAL	void	insert_linear_section(BOND*,CURVE*,double*,double*,Front*);
LOCAL	void	insert_polynomial_curve_section(BOND*,CURVE*,double,double,
	                                        double*,int,Front*);
#endif /* defined(TWOD) */
#if defined(THREED)
LOCAL	void	init_injection_inlet_jet3d(INIT_DATA*,INIT_PHYSICS*);
LOCAL   void 	map_rectangle_to_circle(double*,double*,double,double,double,double);
LOCAL	void	make_jet_surface(Front*,double,double,int,double,
				 COMPONENT,COMPONENT,COMPONENT);
LOCAL 	void  	reset_surface_sort_status(SURFACE*);
LOCAL 	void  	spherical_surface(SURFACE*,double,double,double*,RECT_GRID*);
LOCAL 	void 	surfaces_on_bonds_of_curve(INTERFACE*);
LOCAL  	double   jet_inner_func(POINTER,double*);
LOCAL  	double   jet_outer_func(POINTER,double*);
LOCAL  	double   jet_nozzle_func(POINTER,double*);
LOCAL  	double   jet_plane_func(POINTER,double*);
LOCAL  	double   jet_top_bdry_func(POINTER,double*);
LOCAL  	double   jet_bottom_bdry_func(POINTER,double*);
LOCAL  	boolean    make_jet_surfaces(RECT_GRID*,COMPONENT,COMPONENT,COMPONENT,
                        double (*func1)(POINTER,double*),POINTER,
                        double (*func2)(POINTER,double*),POINTER,
			double (*func3)(POINTER,double*),POINTER,
		        double (*func4)(POINTER,double*),POINTER,
		        SURFACE**,CURVE**,int*,int*);
/*TMP*/
#define RANGE_GRID 1.0
#define NUMBER_GRID 40.0
#define epsilon RANGE_GRID/NUMBER_GRID*0.001

LOCAL   void    assign_exterior_comp(COMPONENT***, RECT_GRID);
LOCAL   boolean 	exist_func_root(double,double);
LOCAL   int     install_grid_crx2(double (*func)(POINTER,double*),POINTER,
                        EG_CRX*,RECT_GRID,COMPONENT,COMPONENT,int);

LOCAL   int     install_grid_topbdry_crx(double (*func)(POINTER,double*),POINTER,
		        EG_CRX*,RECT_GRID,COMPONENT,COMPONENT,int);
LOCAL   int     install_grid_bottombdry_crx(double (*func)(POINTER,double*),POINTER,
                        EG_CRX*,RECT_GRID,COMPONENT,COMPONENT,int);
LOCAL   int     install_curve_crx(double (*func_1)(POINTER,double*),POINTER,
                        double (*func_2)(POINTER,double*),POINTER,EG_CRX*,
			RECT_GRID,int*,POINT***,COMPONENT,COMPONENT,COMPONENT,
			int**,int*);
LOCAL   int 	face_crx_in_dir(double (*func1)(POINTER,double*),POINTER, 
			double (*func2)(POINTER,double*),POINTER,double***,double*,int);
LOCAL	double	return_sign(double f);
LOCAL	boolean	find_func_root(double (*func)(POINTER,double*),
                        POINTER,double*,double*,double*);
LOCAL   double 	partial_func_derivative(double(*func)(POINTER,double*),
                        POINTER,double*,int);
LOCAL   double   partial_func_derivative2(double(*func)(POINTER,double*),
                        POINTER,double*,int,int);
LOCAL   int     which_3comp(COMPONENT,COMPONENT,COMPONENT,COMPONENT,COMPONENT,
		             COMPONENT,COMPONENT); 
/*TMP END*/
#endif /* defined(THREED) */

typedef struct {
        double *cen;             /* Center of buttom of the jet */
        double r1;               /* Lengths of radii in upper section */
        double r2;               /* Lengths of radii in lower section */
        double *h;               /* height of different sections */
        double ratio;            /* ratio of height from contact 
	                           surface to top of the nozzle */
} JET_PARAMS;

EXPORT  void init_fuel_injection_jet(
        INIT_DATA       *init,
        INIT_PHYSICS    *ip)
{
	switch (ip->root->front->rect_grid->dim)
	{
#if defined(TWOD)
	case 2:
	    init_fuel_injection_jet2d(init,ip);
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
        case 3:
            init_3comp_jet3d(init,ip);
            /*init_fuel_injection_jet3d(init,ip); */
            break;
#endif /* defined(THREED) */
	default:
	    screen("ERROR in init_fuel_injection_jet(), "
		   "unsupported or invalid dimension %d\n",
		   ip->root->front->rect_grid->dim);
	    clean_up(ERROR);
	    break;
	}
}	/* end init_fuel_injection_jet */

/*
*		init_fuel_injection_jet2d():
*
*	Initializes a high pressure chamber for the injection jet.
*
*
*       |           
*       |           
*       |<-n_wu------------------------------>2_________________________
*	|                  _                  /#########################|
*	|                  |                 /##########################|
*	|                 n_hu              /###########################|
*       |                  |               /############################|
*       |                  -              /#############################|
*       |<-n_wm------------------------->/##############################|
*       |                  _             \##############################|
*       |                  |              \#############################|
*	|                 nh_l             \############################|
*       |                  |                \###########################|
*       |                  -                 \##########################|
*	|<-nw_l------------------------------>3          |##############|
*	|                                            _   |##############|
*	|                                            |   |##############|
*	|                                            |   |##############|1
*	|                                           c_h  |##############|  _
*	|                                            |   |##############|  |
*	|                                            |   |##############|  |
*	|                                            -  4|##############|  |
*	|                                           _    |##############|  |
*	|                                           |    |##############|  |
*	|                                           |    |##############|  |
*	|                                           |    |##############|  |
*	|                                          h_i   |##############|  h_o
*	|                                           |    |##############|  |
*	|                                           |    |##############|  |
*	|                                           |    |##############|  |
*	|<-------------------------r_i--------------|--->|##############|  |
*	|                                           |    |##############|  |
*	|                                           |    |##############|  |
*	|                                           |    |<-----w_w---->|  |
*	|                                           |    |##############|  |
*	|___________________________________________-____|##############|  -
*
*	Points 1 and 2, and/or points 3 and 4 may be optionally connected
*	by an elliptical section to create a curved outer (respectively inner)
*	boundary cap.
*/

LOCAL	void init_fuel_injection_jet2d(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front      *front = ip->root->front;
	INTERFACE  *intfc = front->interf;
	RECT_GRID  *gr = front->rect_grid;
	BOND       *b;
	CURVE      *wall, *inwall, *jet;
	CURVE      **curves;
	POINT      *p;
	NODE       *ns, *ne;
	Gas_param  *paramsa, *paramsb, *paramsc;
	double      cen[2], rad[2];
	double      alpha;
	double      *L = gr->GL;
	double      *U = gr->GU;
	double	   *h = gr->h;
	double	   coords[MAXD];
	double      nor[3];
	double      y0, y1, y0p, y1p, x0, x1, c[4];
	double      theta;
	char 	   s[1024];
	int        i, dim = gr->dim;
	double      r_i, w_w, h_i, h_o, nw_l, c_h, nw_m, nh_l, nw_u, nh_u;
	double      nh_max;
	double      tol = MIN_SC_SEP(intfc);
	boolean       crop_x_bdry = NO;

	debug_print("jet","Entered init_fuel_injection_jet2d()\n");

	set_obstacle_comp_type(comp_type(COMPOBST),front);

	screen("\nFuel injection jet enters from below.\n");
	(void) sprintf(s,
            "The jet is initialized as a half jet centered "
            "about r = %g of the computational domain. ",L[0]);
	screen_print_long_string(s);

	i = 0;
	screen("The fuel chamber and nozzle may be entered in "
	       "one of three ways\n"
	       "\t(1) A preset geometry consisting of a fuel chamber,\n"
	       "\t    an elliptical cross section cap, an a Lavel nozzle\n"
	       "\t(2) Direct entry of the chamber wall and nozzle shape\n"
	       "\t(3) Read the coordinates of the chamber wall and nozzle\n"
	       "\t    from a file\n");
	screen("Enter choice (default = %d): ",i);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscanf(s,"%d",&i) != 1) || (i < 1) || (i > 3))
	    {
	        screen("ERROR in init_fuel_injection_jet2d(), "
		       "invalid choice = %d\n",i);
	        clean_up(ERROR);
	    }
	}
	    
	if ((i == 2) || (i == 3))
	{
	    double *x, *y;
	    int   N;

	    screen("The chamber wall is assumed to be oriented so that the\n"
		   "\tright side of the curve bounding the chamber wall is\n"
		   "\tadjacent to the interior of the computational domain.\n");
	    if (i == 2)
	    {
		screen("Enter the number of points on the curve: ");
		if ((Scanf("%d",&N) != 1) || (N <= 0))
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid number of points\n");
	            clean_up(ERROR);
		}
		uni_array(&x,N,FLOAT);
		uni_array(&y,N,FLOAT);
		screen("\nEnter %d ordered pairs for the coordinates of the "
		       "chamber wall.\n",N);
		for (i = 0; i < N; ++i)
		{
		    screen("\t: ");
		    if ((Scanf("%f %f\n",x+i,y+i) != 2) ||
			    (x[i] < L[0]) || (y[i] < L[1]) ||
			    (U[0] < x[i]) || (U[1] < y[i]))
		    {
	                screen("ERROR in init_fuel_injection_jet2d(), "
		               "invalid input of coordinate pair\n");
	                clean_up(ERROR);
		    }
		}
	    }
	    else
	    {
	        FILE *file;
		double tmp;
		screen("The point data for the chamber wall should consist of\n"
		       "\tprecisely N pairs of ASCII numbers. The number of\n"
		       "\tpoints in the wall is obtained counting the number\n"
		       "\tpairs in the file\n");
		screen("Enter the file name for the chamber wall data: ");
		(void) Gets(s);
		if ((s[0] == '\0') || ((file = fopen(s,"r")) == NULL))
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "can't open %s\n",s);
	            clean_up(ERROR);
		}
		for (N = 0; (fscan_float(file,&tmp) == 1); ++N);
		if (N%2)
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid (odd) number of points in %s\n",s);
	            clean_up(ERROR);
		}
		N /= 2;
		uni_array(&x,N,FLOAT);
		uni_array(&y,N,FLOAT);
		rewind(file);
		for (i = 0; i < N; ++i)
		{
		    (void) fscan_float(file,x+i);
		    (void) fscan_float(file,y+i);
		}
		(void) fclose(file);
	    }
	    coords[0] = x[0];
	    coords[1] = y[0];
	    ns = make_node(Point(coords));
	    set_is_bdry(ns);
	    node_type(ns) = FIXED_NODE;
	    coords[0] = x[N-1];
	    coords[1] = y[N-1];
	    ne = make_node(Point(coords));
	    set_is_bdry(ne);
	    node_type(ne) = FIXED_NODE;
	    wall = make_curve(COMPOBST,COMPA,ns,ne);
	    wave_type(wall) = NEUMANN_BOUNDARY;
	    start_status(wall) = end_status(wall) = FIXED;
	    for (i = 1; i < (N-1); ++i)
	    {
	        coords[0] = x[i];
	        coords[1] = y[i];
		insert_point_in_bond(Point(coords),wall->last,wall);
	    }
	    free_these(2,x,y);

	    screen("Redistribute the input curve (default = yes): ");
	    (void) Gets(s);
	    if (strncasecmp(s,"n",1) != 0)
	    {
		if (!equi_curve_redistribute(front,wall,YES))
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
			   "redistribution of wall boundary failed\n");
		    clean_up(ERROR);
		}
	    }
	}
	else
	{
	    screen("Enter the inner radius of the fuel container: ");
	    (void) Scanf("%f\n",&r_i);
	    if ((r_i <= 0.0) || (U[0] < (L[0] + r_i)))
	    {
	        screen("ERROR in init_fuel_injection_jet2d(), "
		       "invalid inner radius r_i = %"FFMT"\n",r_i);
	        (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	        clean_up(ERROR);
	    }

	    screen("Enter the width of the fuel container wall: ");
            (void) Scanf("%f\n",&w_w);
	    if (w_w <= 0.0)
	    {
	        screen("ERROR in init_fuel_injection_jet2d(), "
		       "invalid container width = %"FFMT"\n",w_w);
	        clean_up(ERROR);
	    }
   
	    h_i = 0.0;
	    screen("Enter the height of the interior vertical wall\n"
		   "\tof the fuel chamber (default = %g): ",h_i);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        sscan_float(s,&h_i);
	        if ((h_i < 0.0) || (U[1] < (L[1] + h_i)))
	        {
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid interior wall height h_i = %"FFMT"\n",h_i);
	            (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	            clean_up(ERROR);
	        }
	    }

	    h_o = h_i;
	    screen("Enter the height of the exterior vertical wall\n"
		   "\tof the fuel chamber (default = %g): ",h_o);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        sscan_float(s,&h_o);
	        if ((h_o < 0.0) || (U[1] < (L[1] + h_o)))
	        {
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid exterior wall height h_o = %"FFMT"\n",h_o);
	            (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	            clean_up(ERROR);
	        }
	    }

	    c_h = 0.0;
	    screen("Enter the height of the cap on the fuel chamber\n"
	           "\t(default or zero = no cap): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        if ((sscan_float(s,&c_h) != 1) || (c_h < 0.0) ||
		        (U[1] < (L[1] + h_i + c_h) ))
	        {
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid cap height c_h = %"FFMT"\n",c_h);
	            (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	            clean_up(ERROR);
	        }
	    }
	    screen("Enter the width of the lower (inlet) "
	           "section of the nozzle: ");
	    (void) Scanf("%f\n",&nw_l);
	    if ((nw_l <= 0.0) || (r_i < nw_l) || (U[0] < (L[0] + nw_l)))
	    {
	        screen("ERROR in init_fuel_injection_jet2d(), "
		       "invalid lower inlet width nw_l = %"FFMT"\n",nw_l);
	        (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	        clean_up(ERROR);
	    }
	    nw_u = nw_m = nw_l;
    
	    nh_l = 0.0;
	    screen("Enter the height of the lower section of the nozzle\n"
	           "\t(default or zero = no lower section): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        if ((sscan_float(s,&nh_l) != 1) || (nh_l < 0.0) ||
		        (U[1] < (L[1] + h_i + c_h + nh_l) ))
	        {
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid lower section nozzle "
			   "height nh_l = %"FFMT"\n",nh_l);
	            (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	            clean_up(ERROR);
	        }
	    }
	    if (nh_l > 0.0)
	    {
	        screen("Enter the width of the middle inlet of the nozzle\n"
	               "\t(default = %g): ",nw_m);
	        (void) Gets(s);
	        if (s[0] != '\0')
	        {
	            if ((sscan_float(s,&nw_m) != 1) || (nw_m < 0.0) ||
		            ((r_i + w_w) < nw_m))
	            {
	                screen("ERROR in init_fuel_injection_jet2d(), invalid "
			        "middle inlet width nw_m = %"FFMT"\n",
		               nw_m);
	                clean_up(ERROR);
	            }
	        }
	    }

	    nh_u = 0.0;
	    screen("Enter the height of the upper section of the nozzle\n"
	           "\t(default or zero = no upper section): ");
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
	        if ((sscan_float(s,&nh_u) != 1) || (nh_u < 0.0) ||
		        (U[1] < (L[1] + h_i + c_h + nh_l + nh_u) ))
	        {
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid upper section nozzle "
			   "height nh_u = %"FFMT"\n",
		           nh_u);
	            (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	            clean_up(ERROR);
	        }
	    }
	    if (nh_u > 0.0)
	    {
	        screen("Enter the width of the upper outlet of the nozzle\n"
	               "\t(default = %g): ",nw_u);
	        (void) Gets(s);
	        if (s[0] != '\0')
	        {
	            if ((sscan_float(s,&nw_u) != 1) || (nw_u < 0.0) ||
		            ((r_i + w_w) < nw_u))
	            {
	                screen("ERROR in init_fuel_injection_jet2d(), "
		               "invalid oulet section "
			       "nozzle width nw_u = %"FFMT"\n",
		               nw_u);
	                clean_up(ERROR);
	            }
	        }
	    }
	    coords[0] = L[0] + r_i + w_w;
	    coords[1] = L[1];
	    ns = make_node(Point(coords));
	    set_is_bdry(ns);
	    node_type(ns) = FIXED_NODE;
	    coords[0] = L[0] + r_i;
	    coords[1] = L[1];
	    ne = make_node(Point(coords));
	    set_is_bdry(ne);
	    node_type(ne) = FIXED_NODE;
	    wall = make_curve(COMPOBST,COMPA,ns,ne);
	    wave_type(wall) = NEUMANN_BOUNDARY;
	    start_status(wall) = end_status(wall) = FIXED;
	    if (h_o > 0.0)
	    {
	        coords[0] = L[0] + r_i + w_w;
	        coords[1] = L[1] + h_o;
		insert_linear_section(wall->last,wall,Coords(ns->posn),
			              coords,front);
	    }
	    if (h_o < (h_i + c_h + nh_l + nh_u))
	    {
		cen[0] = L[0] + nw_u;
		cen[1] = L[1] + h_o;
		rad[0] = w_w + r_i - nw_u;
		rad[1] = h_i + c_h + nh_l + nh_u - h_o;
                insert_elliptical_curve_section(wall->last,wall,0.0,0.5*PI,
			                        cen,rad,0.0,front);
	    }
	    else
	    {
	        coords[0] = L[0] + nw_u;
	        coords[1] = L[1] + h_i + c_h + nh_l + nh_u;
		insert_linear_section(wall->last,wall,Coords(wall->last->start),
			              coords,front);
	    }
	    if (nh_u > 0.0)
	    {
	        coords[0] = L[0] + nw_m;
	        coords[1] = L[1] + h_i + c_h + nh_l;
		insert_linear_section(wall->last,wall,Coords(wall->last->start),
			              coords,front);
		if (nh_l > 0.0)
		{
	            coords[0] = L[0] + nw_l;
	            coords[1] = L[1] + h_i + c_h;
		    insert_linear_section(wall->last,wall,
			                  Coords(wall->last->start),
					  coords,front);
		}
	    }
	    else
	    {
	        coords[0] = L[0] + nw_l;
	        coords[1] = L[1] + h_i + c_h;
		insert_linear_section(wall->last,wall,Coords(wall->last->start),
			              coords,front);
	    }
	    if (c_h > 0.0)
	    {
		cen[0] = L[0] + nw_l;
		cen[1] = L[1] + h_i;
		rad[0] = r_i - nw_l;
		rad[1] = c_h;
                insert_elliptical_curve_section(wall->last,wall,0.5*PI,0.0,
			                        cen,rad,0.0,front);
	    }
	    insert_linear_section(wall->last,wall,Coords(wall->last->start),
			          Coords(wall->last->end),front);
	}
	never_redistribute(wall) = YES;

	/* Clip the nozzle wall at the upper x-boundary if necessary */

	if (U[0] < Coords(wall->start->posn)[0])
	{
	    ORIENTATION orient = POSITIVE_ORIENTATION;
	    
	    for (b = wall->first; b != NULL; b = b->next)
	    {
		if (Coords(b->end)[0] < U[0])
		    break;
	    }
	    if (b == NULL)
	    {
	        screen("ERROR in init_fuel_injection_jet2d(), "
		       "invalid nozzle, all points are out of bounds\n");
	        clean_up(ERROR);
	    }
	    alpha = (U[0] - Coords(b->start)[0])/
		    (Coords(b->end)[0] - Coords(b->start)[0]);
	    coords[0] = U[0];
	    coords[1] = (1.0 - alpha)*Coords(b->start)[1] +
		        alpha*Coords(b->end)[1];
	    ns = wall->start;
	    Coords(ns->posn)[0] = coords[0];
	    Coords(ns->posn)[1] = coords[1];
	    set_bond_length(wall->first,dim);
	    p = b->end;
	    while(Point_adjacent_to_node(wall,orient) != p)
		(void) delete_point_adjacent_to_node(front,wall,orient);
	    crop_x_bdry = YES;
	}

	nh_max = Coords(wall->start->posn)[1];
	for (b = wall->first; b != NULL; b = b->next)
	    if (nh_max < Coords(b->end)[1])
		nh_max = Coords(b->end)[1];

	coords[1] = nh_max;
	screen("Enter the height at which the jet connects "
	       "to the nozzle wall\n\t(default = %g): ",nh_max-L[1]);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscan_float(s,coords+1) != 1) || (nh_max < coords[1]))
	    {
	        screen("ERROR in init_fuel_injection_jet2d(), "
		       "invalid jet height, height %g too high\n",coords[1]);
	        clean_up(ERROR);
	    }
	    coords[1] += L[1];
	}

	for (b = wall->last; b != NULL; b = b->prev)
	{
	    if (Coords(b->start)[1] > nh_max - tol*h[1])
		break;
	}
	alpha = (coords[1] - Coords(b->start)[1])/
	        (Coords(b->end)[1] - Coords(b->start)[1]);
	coords[0] = (1.0 - alpha)*Coords(b->start)[0] + alpha*Coords(b->end)[0];
	if (_scaled_separation(coords,Coords(b->start),h,dim) < tol)
	    p = b->start;
	else if (_scaled_separation(coords,Coords(b->end),h,dim) < tol)
	{
	    p = b->end;
	    if (b->next == NULL)
		p = b->start;
	    else
	        b = b->next;
	}
	else
	{
	    p = Point(coords);
	    insert_point_in_bond(p,b,wall);
	    b = b->next;
	}
	normal(p,Hyper_surf_element(b),Hyper_surf(wall),nor,front);
	curves = split_curve(p,b,wall,COMPOBST,COMPB,COMPOBST,COMPA);
	inwall = curves[1];
	no_slip(Hyper_surf(inwall)) = NO; /*default*/
	screen("The inner injector wall is set as a slip Nuemann ");
	screen("boundary.\nFor no slip Neumann boundary, type 'ns'.\n");
	screen("Enter choice: ");
	(void) Gets(s);
	if (s[0] == 'n' || s[0] == 'N')
	{
	    double adherence_coeff;
	    no_slip(Hyper_surf(inwall)) = YES; 
	    adherence_coeff = 1.0;
	    screen("Enter the adherence coefficient\n");
	    screen("\tadherence_coeff = 0.0 for fully slippery\n");
	    screen("\tadherence_coeff = 1.0 for fully non-slippery\n");
	    screen("\tdefault = %g (non-slippery)\n",adherence_coeff);
	    screen("Enter the value of adherence coefficient: ");
	    (void) Gets(s);
	    sscan_float(s,&adherence_coeff);
	    if (adherence_coeff < 0.0 || adherence_coeff > 1.0)
	    {
		screen("ERROR in init_fuel_injection_jet2d(), "
                       "invalid input of partial slip parameter\n");
                clean_up(ERROR);
	    }
	    if(adherence_coeff == 0.0)
	        no_slip(Hyper_surf(inwall)) = NO; /* turn back to slip */
	    else
	        adherence_coeff(inwall) = adherence_coeff;
	}
	wall = NULL;

	ns = inwall->start;
	node_type(ns) = NEUMANN_NODE;
	screen("Select nozzle-jet node type, choices are\n"
	       "\tNeumann node (n, default)\n"
	       "\tAttached B node (a)\n"
	       "Enter choice: ");
	Gets(s);
	if (s[0] == 'n' || s[0] == 'N')
	    node_type(ns) = NEUMANN_NODE;
	else if (s[0] == 'a' || s[0] == 'A')
	    node_type(ns) = ATTACHED_B_NODE;


	if (nor[0] != 0.0)
	{
	    x1 = Coords(ns->posn)[0];
	    y1 = Coords(ns->posn)[1];
	    y1p = nor[1]/nor[0];

	    screen("Enter the tangent angle of jet surface\n"
		   "\tat the nozzle wall (default = %g): ",
		   degrees(atan2(-nor[1],-nor[0])));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		if (sscan_float(s,&theta) != 1)
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid input of jet tangent at nozzle\n");
		    clean_up(ERROR);
		}
		y1p = tan(radians(theta));
	    }

	    x0 = L[0];
	    y0 = y1 - 0.5*y1p*(x1 - x0);
	    y0p = 0.0;

	    screen("Enter the change in height of the jet from the nozzle\n"
		   "\tto the interior wall or symmetry axis (default = %g): ",
		   y0-y1);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		if (sscan_float(s,&y0) != 1)
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid input of jet height at interior boundary\n");
		    clean_up(ERROR);
		}
		y0 += y1;
	    }
	    screen("Enter the tangent angle of jet surface at the "
		   "interior boundary\n\t(default = %g): ",atan(y0p));
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		if (sscan_float(s,&theta) != 1)
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid input of jet tangent at nozzle\n");
		    clean_up(ERROR);
		}
		y0p = tan(radians(theta));
	    }
	    coords[0] = x0;
	    coords[1] = y0;
	    ne = make_node(Point(coords));
	    set_is_bdry(ne);
	    jet = make_curve(COMPA,COMPB,ns,ne);
	    c[0] = ( (x0*x0*x0*y1 - x1*x1*x1*y0) +
		     3.0*x0*x1*(x1*y0-x0*y1) +
		     x0*x1*(x1*x1*y0p - x0*x0*y1p) +
		     x0*x0*x1*x1*(y1p - y0p) ) / ((x0-x1)*(x0-x1)*(x0-x1));
	    c[1] = ( (x0*x0*x0*y1p - x1*x1*x1*y0p) +
		     2.0*x0*x1*(x0*y0p - x1*y1p) +
		     x0*x1*(x0*y1p - x1*y0p)) / ((x0-x1)*(x0-x1)*(x0-x1));
	    c[2] = ( (x1*x1*y1p - x0*x0*y0p) +
		     2.0*(x1*x1*y0p - x0*x0*y1p) +
		     3.0*(x0*y0 - x1*y1) +
		     x0*x1*(y1p - y0p) +
		     3.0*(x1*y0 - x0*y1) ) / ((x0-x1)*(x0-x1)*(x0-x1));
	    c[3] = ( 2.0*(y1 - y0) + (x0 - x1)*(y0p + y1p) ) /
		     ((x0-x1)*(x0-x1)*(x0-x1));
            insert_polynomial_curve_section(jet->first,jet,x1,x0,c,3,front);
	}
	else
	{
	    cen[0] = L[0];
	    cen[1] = Coords(ns->posn)[1];
	    rad[0] = Coords(ns->posn)[0];
	    rad[1] = 2.0*h[1];
	    screen("Enter the change in height of the jet from the nozzle "
		   "to the interior wall or symmetry axis (default = %g): ",
		   rad[1]);
	    (void) Gets(s);
	    if (s[0] != '\0')
	    {
		if ((sscan_float(s,rad+1) != 1) || (rad[1] < 0.0))
		{
	            screen("ERROR in init_fuel_injection_jet2d(), "
		           "invalid input of jet height at interior boundary\n");
		    clean_up(ERROR);
		}
	    }
	    coords[0] = cen[0];
	    coords[1] = cen[1] + rad[1];
	    ne = make_node(Point(coords));
	    set_is_bdry(ne);
	    jet = make_curve(COMPA,COMPB,ns,ne);
            insert_elliptical_curve_section(jet->first,jet,0.0,0.5*PI,
			                    cen,rad,0.0,front);
	}
	start_status(jet) = INCIDENT;
	end_status(jet) = INCIDENT;
	wave_type(jet) = CONTACT;

	(void) prompt_for_eos_params(init,ip,YES,"");
	paramsa = prompt_for_eos_params(init,ip,YES,"\n\tfor the injected gas");
	paramsb = prompt_for_eos_params(init,ip,YES,"\n\tfor the ambient gas");
	paramsc = prompt_for_eos_params(init,ip,YES,"\n\tfor the bubbles in the injected gas");
	prompt_for_ambient_state(comp_type(COMPA),paramsa,
                                 " for the injected gas",front,init);
	prompt_for_ambient_state(comp_type(COMPB),paramsb,
                                 " for the ambient gas",front,init);
	prompt_for_ambient_state(comp_type(COMPC),paramsc,
	                         " for the bubbles in the injected gas",front,init);
	surface_tension(jet) =  prompt_for_surface_tension(CONTACT,
						           "for the jet ");

	rect_boundary_type(intfc,1,0) = MIXED_TYPE_BOUNDARY;
	if (crop_x_bdry)
	  rect_boundary_type(intfc,0,1) = MIXED_TYPE_BOUNDARY;
	

	debug_print("jet","Left init_fuel_injection_jet2d()\n");
}	/* end init_fuel_injection_jet2d */

#if defined(THREED)
/*
*		init_fuel_injection_jet3d():
*
*	Initializes a gas chamber for the injection jet.
*
*
*       |           
*       |                                      
*       |-------------- r_n ---->|______________________________________ 
*	|    -                   |######################################|
*	|    |                   |######################################|
*	|    |          COMPB    |######################################|
*       |    |                   |######################################|
*       |   h_u                  |######################################|
*       |____|___contact_________|##################### COMPOBST #######|
*       |    |                   |######################################|
*       |    |h2                 |######################################|
*	|    |                   |######################################|
*       |    -                -         \###############################|
*       |                     |               \#########################|  
*	|                     |                   \#####################|  
*	|                     |                     \###################|  
*	|                    h_l                      \#################|  
*	|                     |                        |################|  
*	|                     |                         |###############|  
*	|           COMPA     |                          \##############|  
*	|                     -                          |##############|  -
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  h_i
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|<-------------------------r_i------------------>|##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |<-----w_w---->|  |
*	|                                                |##############|  |
*	|________________________________________________|##############|  -
*
*/
LOCAL	void init_fuel_injection_jet3d(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front           *front = ip->root->front;
	INTERFACE       *intfc = front->interf;
	RECT_GRID       *gr = front->rect_grid;
	CURVE           *wall, *inwall, **jet;
	CURVE           **curves;
	SURFACE	        *surf_annulus, *surf_orifice, *surf_cap, *surf_bdry1,
		        *surf_bdry2, *surf_bdry3, *surf_bdry4;
	SURFACE	        **surfs ;
	JET_PARAMS      jet_paramsa, jet_paramsb;
	BDRY_BOX_PARAMS bp;
	Gas_param       *gas_paramsa, *gas_paramsb,  *gas_paramsc;
	double           *L = gr->GL;
	double           *U = gr->GU;
	double	        *h = gr->h;
	double	        coords[MAXD];
	char 	        s[1024];
	int             i, dim = gr->dim;
	double           r_i, w_w, h_i, r_n, h_l, h_u;
	double           cen[3], h_in[3], h_out[3], ratio;
	double           tol = MIN_SC_SEP(intfc);
        int             is,ic;
	int             me[3], *G;
	double           surf_tension = 0.0;

	debug_print("jet","Entered init_fuel_injection_jet3d()\n");

	set_obstacle_comp_type(comp_type(COMPOBST),front);

	/* read geometrics of the jet */
	screen("\nFuel injection jet enters from below.\n");
	screen("Please enter the geometric parameters of the fuel jet.\n");
	
	screen("Enter the center values of the floor of fuel container: ");
	(void) Scanf("%f %f %f\n",&cen[0],&cen[1],&cen[2]);
	if ((L[0] > cen[0]) || (U[0] < cen[0]) ||
	    (L[1] > cen[1]) || (U[1] < cen[1]) ||
	    (L[2] > cen[2]) || (U[2] < cen[2]))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid center of the floor x0 = %"FFMT", "
	       "y0 = %"FFMT", z0 = %"FFMT"\n",cen[0],cen[1],cen[2]);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	    (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	    clean_up(ERROR);
	}

	screen("Enter the inner radius of the fuel container: ");
	(void) Scanf("%f\n",&r_i);
	if ((r_i <= 0.0) || (U[0] < (L[0] + r_i)))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid inner radius r_i = %"FFMT"\n",r_i);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    clean_up(ERROR);
	}

	screen("Enter the width of the fuel container wall: ");
        (void) Scanf("%f\n",&w_w);
	if (w_w <= 0.0)
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
		       "invalid container width = %"FFMT"\n",w_w);
	    clean_up(ERROR);
	}
   
	screen("Enter the radium of the nozzle: ");
	(void) Scanf("%f\n",&r_n);
	if ((r_n <= 0.0) || (r_i < r_n) || (U[0] < (L[0] + r_n)))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid nozzle radium r_n = %"FFMT"\n",r_n);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    clean_up(ERROR);
	}
   
	h_i = 0.0;
	screen("Enter the height of the interior vertical wall\n"
	   "\tof the fuel chamber (default = %g): ",h_i);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    sscan_float(s,&h_i);
	    if ((h_i < 0.0) || (U[2] < (L[2] + h_i)))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
		           "invalid interior wall height h_i = %"FFMT"\n",h_i);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        clean_up(ERROR);
	    }
	}

	h_l = 0.0;
	screen("Enter the height of the curve section of the nozzle\n"
	       "\t(default or zero = no lower section): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscan_float(s,&h_l) != 1) || (h_l < 0.0) ||
	        (U[2] < (L[2] + h_i + h_l) ))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
	           "invalid curve section nozzle "
		   "height nh_l = %"FFMT"\n",h_l);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        clean_up(ERROR);
	    }
	}

	h_u = 0.0;
	screen("Enter the height of the cylinder section of the nozzle\n"
	           "\t(default or zero = no cylinder section): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscan_float(s,&h_u) != 1) || (h_u < 0.0) ||
	        (U[2] < (L[2] + h_i + h_l + h_u) ))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
	           "invalid cylinder section nozzle "
		   "height nh_u = %"FFMT"\n",
	           h_u);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        clean_up(ERROR);
	    }
	}
	
	ratio = 0.5;
	screen("Enter the ratio of the height of the contace section\n"
	       "\t to the top of the nozzle (default = %f): ",ratio);
	(void) Scanf("%f\n",&ratio);
	if ((ratio > 1.0))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid ratio of contace in nozzle = %"FFMT"\n",ratio);
	    clean_up(ERROR);
	}
	
	(void) prompt_for_eos_params(init,ip,YES,"");
	gas_paramsb = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the ambient gas");
	gas_paramsa = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the injected gas");
	gas_paramsc = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the bubbles in the injected gas");

	prompt_for_ambient_state(comp_type(COMPB),gas_paramsb,
                                 " for the ambient gas",front,init);
	prompt_for_ambient_state(comp_type(COMPA),gas_paramsa,
                                 " for the injected gas",front,init);
	prompt_for_ambient_state(comp_type(COMPC),gas_paramsc,
		                 " for the bubbles in the injected gas",front,init);

	/* start make surface of jet */
	h_in[0] = h_i;
	h_in[1] = h_l;
	h_in[2] = h_u;

	h_out[0] = h_i;
	h_out[1] = h_l;
	h_out[2] = h_u;
	
	jet_paramsa.cen = cen;
	jet_paramsa.r1 = r_i;
	jet_paramsa.r2 = r_n;
	jet_paramsa.h = h_in;
	jet_paramsa.ratio = ratio;
	
	jet_paramsb.cen = cen;
	jet_paramsb.r1 = r_i + w_w;
	jet_paramsb.r2 = r_n;
	jet_paramsb.h = h_out;
	jet_paramsb.ratio = ratio;

	bp.L = gr->GL;         bp.U = gr->GU;
	uni_array(&surfs,7,sizeof(SURFACE*));
	uni_array(&jet,3,sizeof(CURVE*));
         
	is = ic = 0;
       	make_jet_surfaces(gr,COMPA,COMPOBST,COMPB,
                        jet_inner_func,(POINTER)&jet_paramsa,
		        jet_outer_func,(POINTER)&jet_paramsb,
			jet_top_bdry_func,(POINTER)&bp,
			jet_bottom_bdry_func,(POINTER)&bp,
			surfs,jet,&is,&ic);
         
	printf("End make_jet_surfaces,is = %d,ic = %d\n",is,ic);
	
	
        for (i = 0; i < ic; i++)
	{
	    switch (jet[i]->number)
	    {
	        case 1:    
		    curve_type(jet[i]) = ATTACHED_B_CURVE;	
		    break;
		case 2:
		    curve_type(jet[i]) = FIXED_CURVE;
		    break;
		case 3:
		    curve_type(jet[i]) = FIXED_CURVE;
		    break;
		default:
		    screen("UNKNOWN CURVE number,code needed!\n");
		    clean_up(ERROR);
	    }
	}

	/* Every subdomain read this info */
        surf_tension =  prompt_for_surface_tension(CONTACT,
                                "for the jet ");
	for(i = 0; i < is; i++)
	{
 	    switch (surfs[i]->number)
	    {
 	        case 1:
		    wave_type(surfs[i]) = NEUMANN_BOUNDARY;	
	            break;
		case 2:
		    wave_type(surfs[i]) = NEUMANN_BOUNDARY;
		    break;
		case 3:
                    surface_tension(surfs[i]) =  surf_tension;
	            wave_type(surfs[i]) = CONTACT;
		    break;
		case 4:
		    wave_type(surfs[i]) = DIRICHLET_BOUNDARY;
		    set_is_bdry(surfs[i]);
		    break;
		case 5:
		    wave_type(surfs[i]) = PASSIVE_BOUNDARY;
		    set_is_bdry(surfs[i]);
	            break;		    
		case 6:
		    wave_type(surfs[i]) = DIRICHLET_BOUNDARY;
		    set_is_bdry(surfs[i]);
		    break;
		case 7:
                    wave_type(surfs[i]) = DIRICHLET_BOUNDARY;
                    set_is_bdry(surfs[i]);
                    break;
                default:
                    screen("UNKNOWN SURFACE number,code needed!\n");
                    clean_up(ERROR);
	    }
	}

	/* Set the top boundary in Z dirction as  MIXED_TYPE.
	 * The inlet of the container is there.
	 */
	/*rect_boundary_type(intfc,2,1) = MIXED_TYPE_BOUNDARY; */

	debug_print("jet","Left init_fuel_injection_jet3d()\n");
}	/* end init_fuel_injection_jet3d */

#endif /* defined(THREED) */

EXPORT	void init_injection_inlet_jet(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	switch (ip->root->front->rect_grid->dim)
	{
#if defined(TWOD)
	case 2:
	    init_injection_inlet_jet2d(init,ip);
	    break;
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	    init_injection_inlet_jet3d(init,ip);
	    break;
#endif /* defined(THREED) */
	default:
	    screen("ERROR in init_injection_inlet_jet(), "
		   "unsupported or invalid dimension %d\n",
		   ip->root->front->rect_grid->dim);
	    clean_up(ERROR);
	    break;
	}
}		/*end init_injection_inlet_jet*/

#if defined(TWOD)

		/* possible types of jet curves */

enum _JET_TYPE {
	FULL_JET = 1,
	HALF_JET = 2
};
typedef enum _JET_TYPE JET_TYPE;


LOCAL	void init_injection_inlet_jet2d(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	COMPONENT	ecomp;
	CURVE		*cur, *inlet;
	ELLIPSOID	Ellip;
	HYPER_SURF	*jet_hs;
	Front		*front = ip->root->front;
	INTERFACE	*intfc = front->interf;
	Locstate	si;
	NODE		*ns, *ne;
	NODE		*ins, *ine;
	RECT_GRID	*gr = front->rect_grid;
	char		s[1024];
	double		A, beta, beta_max, v[MAXD];
	double		theta0;
	double		iradius, imid, ilower, iupper;
	double		jradius, jmid, jlower, jupper;
	double		*L = front->rect_grid->GL, *U = front->rect_grid->GU;
	double		coords[3];
	double		ci;
	JET_TYPE	jet;
	boolean		attached_jet;
	int		i, dim = intfc->dim;
	Gas_param	*params, *params2;

	debug_print("jet","Entered init_injection_inlet_jet2d()\n");

	if (dim != 2)
	{
	    screen("ERROR in init_injection_inlet_jet2d(), "
		   "unsupported dim %d\n",dim);
	    clean_up(ERROR);
	    debug_print("jet","Left init_injection_inlet_jet2d()\n");
	    return;
	}

	iradius = 0.1*(U[0] - L[0]);
	imid = 0.5*(U[0] + L[0]);
	switch (gr->Remap.remap)
	{
	case IDENTITY_REMAP:
	    jet = FULL_JET;
	    ilower = imid - iradius;
	    iupper = imid + iradius;
	    break;
	case CYLINDRICAL_REMAP:
	    jet = HALF_JET;
	    ilower = RL_eff(init);
	    iupper = ilower + 2.0*iradius;
	    break;
	default:
	    screen("ERROR in init_injection_inlet_jet2d(), "
		   "unsupported remap\n");
	    clean_up(ERROR);
	    debug_print("jet","Left init_injection_inlet_jet2d()\n");
	    return;
	}
	imid = 0.5*(ilower + iupper);

	screen("\nJet enters from below.\n");
	(void) sprintf(s,
	    "The jet may be initialized as either a half jet centered "
	    "about %s = %g, or a full inlet along the lower boundary "
	    "of the computational domain\n",
	    (gr->Remap.remap == CYLINDRICAL_REMAP) ? "r" : "x",L[0]);
	screen_print_long_string(s);
	screen("Enter the jet type (dflt = %s): ",
	       (jet == FULL_JET) ? "full" : "half");
	(void) Gets(s);
	if ((s[0] == 'F') || (s[0] == 'f'))
	    jet = FULL_JET;
	else if ((s[0] == 'H') || (s[0] == 'h'))
	    jet = HALF_JET;
	switch (jet)
	{
	case FULL_JET:
	    screen("Enter the lower coordinate of the inlet (dflt = %g): ",
		   ilower);
	    (void) Gets(s);
	    if (s[0] != '\0')
		(void) sscan_float(s,&ilower);
	    screen("Enter the upper coordinate of the inlet (dflt = %g): ",
		   iupper);
	    (void) Gets(s);
	    if (s[0] != '\0')
		(void) sscan_float(s,&iupper);
	    if (iupper < ilower)
	    {
	        screen("ERROR in init_injection_inlet_jet2d(), invalid "
		       "input, upper coords less than lower coords\n");
	        clean_up(ERROR);
	        debug_print("jet","Left init_injection_inlet_jet2d()\n");
	        return;
	    }
	    imid   = 0.5*(ilower + iupper);
	    iradius = 0.5*(iupper - ilower);
	    break;

	case HALF_JET:
	    screen("Enter the radius of the inlet (dflt = %g): ",iradius);
	    (void) Gets(s);
	    if (s[0] != '\0')
		(void) sscan_float(s,&iradius);
	    iupper = ilower + iradius;
	    imid = ilower;
	    break;
	}

	attached_jet = YES;
	screen("Do you wish the jet to remain attached to the inlet edge? "
	       "(dflt = %s): ",y_or_n(attached_jet));
	Gets(s);
	if (s[0] == 'y' || s[0] == 'Y')
	    attached_jet = YES;
	if (s[0] == 'n' || s[0] == 'N')
	    attached_jet = NO;

	jradius = iradius;
	jlower = ilower;
	jupper = iupper;
	jmid = imid;
	if (attached_jet == NO)
	{
	    switch (jet)
	    {
	    case FULL_JET:
	        screen("Enter the lower coordinate (<= %g) of the "
		       "jet (dflt = %g): ",ilower,jlower);
	        (void) Gets(s);
	        if (s[0] != '\0')
		    (void) sscan_float(s,&jlower);
		if (jlower > ilower)
		{
	            screen("ERROR in init_injection_inlet_jet2d(), invalid "
		           "input, left jet edge greater "
			   "than left inlet edge\n");
	            clean_up(ERROR);
		}
	        screen("Enter the upper coordinate (>= %g) of the "
		       "jet (dflt = %g): ",iupper,jupper);
	        (void) Gets(s);
	        if (s[0] != '\0')
		    (void) sscan_float(s,&jupper);
	        if (jupper < iupper)
	        {
	            screen("ERROR in init_injection_inlet_jet2d(), invalid "
		           "input, right jet edge less "
			   "than right inlet edge\n");
	            clean_up(ERROR);
	            debug_print("jet","Left init_injection_inlet_jet2d()\n");
	            return;
	        }
	        jmid   = 0.5*(jlower + jupper);
	        jradius = 0.5*(jupper - jlower);
	        break;

	    case HALF_JET:
	        screen("Enter the radius the jet (>= %g) (dflt = %g): ",
		       iradius,jradius);
	        (void) Gets(s);
	        if (s[0] != '\0')
		    (void) sscan_float(s,&jradius);
	        if (jradius < iradius)
	        {
	            screen("ERROR in init_injection_inlet_jet2d(), invalid "
		           "input, jet radius less "
			   "than inlet radius\n");
	            clean_up(ERROR);
	            debug_print("jet","Left init_injection_inlet_jet2d()\n");
	            return;
	        }
	        jupper = jlower + jradius;
	        jmid = jlower;
	        break;
	    }
	}

	set_default_ellipsoid_structure(&Ellip,dim);
	Ellip.cen[0] = jmid;
	Ellip.cen[1] = L[1];
	Ellip.rad[0] = jradius;
	Ellip.rad[1] = A = jradius;;
	Ellip.ThetaS[0] = 0.0;
	Ellip.ThetaE[0] = (jet == FULL_JET) ? PI : 0.5*PI;
	Ellip.closed = NO;
	Ellip.nor_orient = POSITIVE_ORIENTATION;
	Ellip.compin = COMPB;
	Ellip.compout = COMPA;
	Ellip.fpoly = NULL;
	Ellip.hs = NULL;
	Ellip.surf_tension = 0.0;
	Ellip.wv_type = CONTACT;
	Ellip.dim = dim;
	Ellip.untracked = NO;
	Ellip.layer_index = ++num_layers(intfc);

	screen("Enter the initial height of the jet (dflt = %g): ",A);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&A);
	beta_max = atan(jradius/A);
	beta = 0.0;
	screen("Enter the angle in degrees that the jet makes with the "
	       "vertical on\n\t"
	       "its right hand side, the absolute value of "
	       "this angle must be less than\n\t"
	       "%g degrees, (dflt = %g): ",degrees(beta_max),degrees(beta));
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscan_float(s,&beta);
	beta = radians(beta);
	if (beta > beta_max)
	{
	    screen("ERROR in init_injection_inlet_jet2d(), invalid "
		   "input, angle of jet edge with vertical too large\n");
	    clean_up(ERROR);
	    debug_print("jet","Left init_injection_inlet_jet2d()\n");
	    return;
	}
	theta0 = asin((A/jradius)*tan(beta));
	Ellip.rad[0] = jradius/cos(theta0);
	Ellip.rad[1] = A/(1.0 - sin(theta0));
	Ellip.cen[1] = L[1] - Ellip.rad[1]*sin(theta0);
	Ellip.ThetaS[0] = theta0;
	Ellip.ThetaE[0] = (jet == FULL_JET) ? PI - theta0 : 0.5*PI;

	(void) prompt_for_eos_params(init,ip,YES,"");
	params = prompt_for_eos_params(init,ip,YES,"\n\tfor the ambient gas");
	params2 = prompt_for_eos_params(init,ip,YES,"\n\tfor the injected gas");
	prompt_for_ambient_state(comp_type(COMPA),params,
		                 " for the ambient gas",front,init);
	prompt_for_ambient_state(comp_type(COMPB),params2,
		                 " for the injected gas",front,init);
	si = Ambient(comp_type(COMPB));
	for (i = 0; i < dim; ++i)
	    v[i] = vel(i,si);
	if ((gr->Remap.remap == CYLINDRICAL_REMAP) && (jet == HALF_JET))
	{
	    if (v[0] != 0.0)
	    {
		screen("ERROR in init_injection_inlet_jet2d(), "
		       "a cylindrical half jet with vr != 0\n");
		clean_up(ERROR);
	    }
	}
	if (v[1] <= 0.0)
	{
	    (void) printf("WARNING in init_injection_inlet_jet2d(), "
			  "negative injection velocity\n");
	}
	ci = sound_speed(si);
	if (ci > fabs(v[1]))
	{
	    (void) printf("WARNING in init_injection_inlet_jet2d(), "
			  "subsonic flow at inlet\n");
	}

	Ellip.surf_tension = prompt_for_surface_tension(CONTACT,"for the jet ");

	if (debugging("jet"))
	{
	    (void) printf("Ellip\n");
	    print_ellipsoid(&Ellip,intfc);
	}
	jet_hs = make_ellipsoid(&Ellip,COMPB,COMPA,front);
	set_is_bdry(Curve_of_hs(jet_hs)->start);
	set_is_bdry(Curve_of_hs(jet_hs)->end);
	ecomp = exterior_component(intfc);
	if (attached_jet == YES)
	{
	    if (jet == FULL_JET)
	    {
	        ins = Curve_of_hs(jet_hs)->start;
		ine = Curve_of_hs(jet_hs)->end;
		node_type(ins) = node_type(ine) = ATTACHED_B_NODE;
	    }
	    else
	    {
		node_type(Curve_of_hs(jet_hs)->start) = ATTACHED_B_NODE;
		node_type(Curve_of_hs(jet_hs)->end) = NEUMANN_NODE;
	        ins = Curve_of_hs(jet_hs)->start;
	        coords[0] = ilower;
	        coords[1] = L[1];
	        ine = make_node(Point(coords));
	        node_type(ine) = FIXED_NODE;
	        set_is_bdry(ine);
	    }
	}
	else
	{
	    node_type(Curve_of_hs(jet_hs)->start) = NEUMANN_NODE;
	    node_type(Curve_of_hs(jet_hs)->end) = NEUMANN_NODE;

	    coords[0] = iupper;
	    coords[1] = L[1];
	    ins = make_node(Point(coords));
	    set_is_bdry(ins);
	    node_type(ins) = FIXED_NODE;

	    coords[0] = ilower;
	    coords[1] = L[1];
	    ine = make_node(Point(coords));
	    set_is_bdry(ine);
	    node_type(ine) = FIXED_NODE;

	    ns = Curve_of_hs(jet_hs)->start;
	    ne = ins;
	    cur = make_curve(ecomp,Ellip.compin,ns,ne);
	    set_is_bdry(cur);
	    start_status(cur) = end_status(cur) = FIXED;
	    wave_type(cur) = NEUMANN_BOUNDARY;

	    if (jet == FULL_JET)
	    {
	        ns = ine;
	        ne = Curve_of_hs(jet_hs)->end;
	        cur = make_curve(ecomp,Ellip.compin,ns,ne);
	        set_is_bdry(cur);
	        start_status(cur) = end_status(cur) = FIXED;
	        wave_type(cur) = NEUMANN_BOUNDARY;

	        ns = Curve_of_hs(jet_hs)->end;
	        coords[0] = RL_eff(init);
	        coords[1] = L[1];
	        ne = make_node(Point(coords));
	        set_is_bdry(ne);
	        node_type(ns) = FIXED_NODE;
	        cur = make_curve(ecomp,Ellip.compout,ns,ne);
	        set_is_bdry(cur);
	        start_status(cur) = end_status(cur) = FIXED;
	        wave_type(cur) = NEUMANN_BOUNDARY;
	    }
	}

	coords[0] = U[0];
	coords[1] = L[1];
	ns = make_node(Point(coords));
	set_is_bdry(ns);
	node_type(ns) = FIXED_NODE;
	ne = Curve_of_hs(jet_hs)->start;
	cur = make_curve(ecomp,Ellip.compout,ns,ne);
	set_is_bdry(cur);
	start_status(cur) = end_status(cur) = FIXED;
	wave_type(cur) = NEUMANN_BOUNDARY;

	inlet = make_curve(ecomp,Ellip.compin,ins,ine);
	set_is_bdry(inlet);
	wave_type(inlet) = DIRICHLET_BOUNDARY;
	bstate_index(inlet) = prompt_for_boundary_state(wave_type(inlet),
							"inlet",NULL,
						        Ellip.compin,
						        -1,Hyper_surf(inlet),
						        init,ip);

	rect_boundary_type(intfc,1,0) = MIXED_TYPE_BOUNDARY;
	if (debugging("jet"))
	{
	    (void) printf("Interface after init_injection_inlet_jet2d()\n");
	    print_interface(intfc);
	}
	debug_print("jet","Left init_injection_inlet_jet2d()\n");
}		/*end init_injection_inlet_jet2d*/
#endif /* defined(TWOD) */

#if defined(THREED)


LOCAL	void init_injection_inlet_jet3d(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front	  *front = ip->root->front;
	double	  radius, height;
	double	  surf_ten;
	Gas_param *params, *params2;
	COMPONENT JET_EXTERIOR = exterior_component(front->interf);
	COMPONENT JET_INLET    = FIRST_DYNAMIC_COMPONENT;
	COMPONENT JET_AMBIENT  = JET_INLET + 1;
	DEBUG_ENTER(init_injection_inlet_jet3d)

	screen("\nJet enters from below.\n");
	screen("\tInitial jet surface is semi-ellipsoid\n"
	       "\tEnter the height of the initial jet perturbation: ");
	(void) Scanf("%f\n",&height);
	screen("\tEnter the radius of the jet inlet orifice: ");
	(void) Scanf("%f\n",&radius);  

	(void) prompt_for_eos_params(init,ip,YES,"");
	params = prompt_for_eos_params(init,ip,YES,"\n\tfor the ambient gas");
	params2 = prompt_for_eos_params(init,ip,YES,"\n\tfor the injected gas");
	prompt_for_ambient_state(comp_type(JET_AMBIENT),params,
				 " for the ambient gas",front,init);
	prompt_for_ambient_state(comp_type(JET_INLET),params2,
				 " for the jet",front,init);

	surf_ten = prompt_for_surface_tension(CONTACT,"for the jet ");
	make_jet_surface(front,height,radius,BOTTOM,surf_ten,
			 JET_EXTERIOR,JET_INLET,JET_AMBIENT);

	DEBUG_LEAVE(init_injection_inlet_jet3d)
}		/*end init_injection_inlet_jet3d*/

/*
*		make_jet_surface():
*
*	make_jet_surface constructs the initial contact wave
*	for an injection jet.
*
*       Author: J.D. Pinezich 9/3/98
*	Rewritten by: X. L. Li 4/4/03
*
*
*/

LOCAL 	void 	make_jet_surface( 
     	Front	     *front,
	double	     height,
	double 	     radius,
	int	     top_or_bottom,
	double 	     surf_ten,
	COMPONENT    JET_EXTERIOR,
	COMPONENT    JET_INLET,
	COMPONENT    JET_AMBIENT)
{
	BOND	   *b;
	INTERFACE  *intfc = front->interf;
	NODE	   *nodes[8];
	CURVE 	   *curves[8];
	CURVE      **inlet_boundary;
	POINT	   *pmin,*pmel,*pmeu,*pmax;
	POINT 	   *p;
	SURFACE	   *s[3], *surf_annulus, *surf_orifice, *surf_cap;
	SURFACE	   **surfs;
	double	   coords[MAXD], midpoint[MAXD];
	double	   x_o[2],y_o[2],x_i[2],y_i[2];
	double	   edge;
	double      X0, Y0, X, Y, l0, lmax;
	int	   i, j, dim = intfc->dim;
	const int  index[] = { 0, 1, 3, 2 };
	RECT_GRID  *comp_grid = front->rect_grid;
	RECT_GRID  dual_grid;
	double	   *dual_VL, *dual_VU;
	double      *L = comp_grid->GL, *U = comp_grid->GU;
	double      *h = comp_grid->h, hx = h[0], hy = h[1];
	int	   *dual_gmax,i_edge[2];
	boolean	   odd_gmax[2] = {YES,YES};
	PLANE_PARAMS plane_params;
	ELLIP_PARAMS ellipsoid_params;
	BDRY_BOX_PARAMS bp;
	double N[3],P[3],*cen,*rad;
	int         is,ic;

	DEBUG_ENTER(make_jet_surface)

	if (top_or_bottom == BOTTOM)
	{
	    midpoint[2] = comp_grid->L[2];
	    height = fabs(height);
	}
	else
	{
	    midpoint[2]  = comp_grid->U[2];
	    height = -fabs(height);
	}

	/* make x and y coordinates */

	/*
	*  	y_o[1] N3------------c2------------N2    
	*		|                          |    
	*               |                          |    
	*       y_i[1] 	|       N7---c6----N6      |    
	*		|        |         |       |    
	*              c3       c7         c5      c1    
	*               |        |         |       |    
	*       y_i[0]  |       N4---c4----N5      |    
	*               |                          |    
	*               |                          |    
	*       y_o[0] N0------------c0------------N1   
	*                                               
	*            x_o[0]   x_i[0]     x_i[1]  x_o[1] 
	*/

	if (debugging("mk_jet_surf"))
	{
	    (void) printf("comp_grid:\n");
	    print_RECT_GRID_structure(comp_grid);
	}


	/* make surfaces using components :    		  *
	 *                                                *
	 *                                                *
	 *                   JET_AMBIENT                  *
	 *                                                *
	 *                        /\                      *
	 *                       /||\                     *
	 *                        ||                      *
	 *                        ||                      *
	 *                                                *
	 *                     -------                    *
	 *                   /         \                  *
	 *                  /           \                 *
	 *                 /             \                *
	 *                |   JET_INLET   |               *
	 *                |               |               *
	 *   -------------------------------------------- *
	 *                                                *
	 *            
	 *            JET_EXTERIOR                  *
	 *                                                */


	cen = ellipsoid_params.cen;
	rad = ellipsoid_params.rad;
	N[0] = N[1] = 0,0; 	N[2] = 1.0;
	rad[0] = rad[1] = radius;
	rad[2] = height;
	for (i = 0; i < 3; ++i)	
	{
	    P[i] = comp_grid->L[i];
	    cen[i] = comp_grid->L[i] + 0.5*(comp_grid->U[i] 
	    			- comp_grid->L[i]);
	}
	for (i = 0; i < 3; ++i)
	{
            plane_params.N[i] = N[i];
	    plane_params.P[i] = P[i];
	}    
	bp.L = comp_grid->VL;         bp.U = comp_grid->VU;
        for (i = 0; i < 3; ++i)
        {
            if (comp_grid->lbuf[i] != 0) bp.L[i] -= comp_grid->h[i];
            if (comp_grid->ubuf[i] != 0) bp.U[i] += comp_grid->h[i];
        }

	uni_array(&surfs,3,sizeof(SURFACE*));

	make_jet_surfaces(comp_grid,JET_EXTERIOR,JET_INLET,JET_AMBIENT,
			plane_func,(POINTER)&plane_params,
			ellipsoid_func,(POINTER)&ellipsoid_params,
			jet_top_bdry_func,(POINTER)&ellipsoid_params,
			jet_bottom_bdry_func,(POINTER)&bp,
			surfs,inlet_boundary,&is,&ic);

	set_is_bdry(surf_orifice);

	for (i = 4; i < 8; ++i)
	    install_curve_in_surface_bdry(surf_orifice,curves[i],
					  POSITIVE_ORIENTATION);

	oblique_planar_surface_triangulation(surf_orifice,front->rect_grid);

	if (debugging("mk_jet_surf"))
	    summarize_interface("mk_jet_surf","square_hole",intfc,
				XY_PLANE,"make_jet_surface","square_hole");

	/*
	*  Remap square regions to circular regions correspond to the
	*  desired geometry.
	*/

	for (i = 5; i < 8; ++i)
	{
	    CURVE *curve1, *curve2;

	    curve1 = nodes[i]->in_curves[0];
	    curve2 = nodes[i]->out_curves[0];
	    inlet_boundary[0] = join_curves(curve1,curve2,NO_COMP,NO_COMP,NULL);
	    if (inlet_boundary[0] == NULL)
	    {
		screen("ERROR in make_jet_surface(), join_curves failed\n");
		clean_up(ERROR);
	    }
	    (void) delete_node(nodes[i]);
	    nodes[i] = NULL;
	}
	for (i = 4; i < 8; ++i)
	    curves[i] = NULL;

	s[0] = surf_annulus;
	s[1] = surf_orifice;
	s[2] = surf_cap;
	for (i = 0; i < 3; ++i)
	    reset_surface_sort_status(s[i]);

	p = inlet_boundary[0]->first->start;
	map_rectangle_to_circle(Coords(p),midpoint,X0,Y0,radius,lmax);
	sorted(p) = YES;
	for (b = inlet_boundary[0]->first; b != inlet_boundary[0]->last; 
			b = b->next)
	{
	    p = b->end;
	    map_rectangle_to_circle(Coords(p),midpoint,X0,Y0,radius,lmax);
	    set_bond_length(b,dim);
	    sorted(p) = YES;
	}
	set_bond_length(inlet_boundary[0]->last,dim);

	for (i = 0; i < 3; ++i)
	{
	    TRI *t;
	    for (t = first_tri(s[i]); !at_end_of_tri_list(t,s[i]);
		 t = t->next)
	    {
	        for (j = 0; j < 3; ++j)
	        {
		    p = Point_of_tri(t)[j];
	            if (sorted(p) == NO)
		    {
	                map_rectangle_to_circle(Coords(p),midpoint,X0,Y0,
						radius,lmax);
	                sorted(p) = YES;
		    }
	        }
	    }
	}
        reset_surface_sort_status(surf_cap);
	spherical_surface(surf_cap,height,radius,midpoint,comp_grid);
	reset_normal_on_intfc(intfc);

	interface_reconstructed(intfc) = NO;

	if (debugging("mk_jet_surf"))
	    summarize_interface("mk_jet_surf","round_hole",intfc,
				XY_PLANE,"make_jet_surface","round_hole");

	if (debugging("mk_jet_surf"))
	    surfaces_on_bonds_of_curve(intfc);

	if (top_or_bottom == BOTTOM)
            rect_boundary_type(intfc,2,0) = MIXED_TYPE_BOUNDARY;
	else
            rect_boundary_type(intfc,2,1) = MIXED_TYPE_BOUNDARY;

	if (!i_consistent_interface(intfc))
	{
	    screen("ERROR in make_jet_surface(), inconsistent interface\n");
	    clean_up(ERROR);
	}

	if (debugging("mk_jet_surf"))
	    points_of_interface(intfc);

	DEBUG_LEAVE(make_jet_surface)
}		/*end make_jet_surface*/

LOCAL   void 	spherical_surface(
	SURFACE   *s,
	double 	  height,
	double 	  r_circle,
	double 	  *m,
	RECT_GRID *comp_grid)
{
	TRI         *t;
	POINT       *p;
	int         j;
	double       x, y, r;
	double       *h = comp_grid->h;
	double       rtol, ztol;
	static const double tolerance = 0.01; /*TOLERANCE*/

	rtol = tolerance*hypot(h[0],h[1]);
	ztol = tolerance*h[2];
	for (t = first_tri(s); !at_end_of_tri_list(t,s); t = t->next)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(t)[j];
	        if (sorted(p) == NO)
		{
		    x = Coords(p)[0];
		    y = Coords(p)[1];
		    r = hypot(x-m[0],y-m[1]);
		    if (fabs(r - r_circle) < rtol)
		        Coords(p)[2] = m[2];
		    else
		    {
			double ratio = r/r_circle;
			Coords(p)[2] = (ratio < 1.0) ?
		                       height*sqrt(1.0 - sqr(ratio)) : m[2];
			if (fabs(Coords(p)[2] - m[2]) < ztol)
			    Coords(p)[2] = m[2];
		    }
		    sorted(p) = YES;
		}
	    }
	}
} 		/*end spherical_surface*/


LOCAL   void 	reset_surface_sort_status(
	SURFACE		*s)
{
	TRI *t;

	for (t = first_tri(s); !at_end_of_tri_list(t,s); t = t->next)
	{
	    sorted(Point_of_tri(t)[0]) = NO;
	    sorted(Point_of_tri(t)[1]) = NO;
	    sorted(Point_of_tri(t)[2]) = NO;
	}
}       	/*end reset_surface_sort_status*/

/*                                                                      *
 *	int   surfaces_on_bonds_of_curve()                              *
 *                                                                      *
 *      loop through all bonds in interface                             *
 *      printout surfaces of curve_of_bond and surfaces of bond         * 
 *      return 0 if order of surface storage for curves agrees with     *
 *      order for bonds.                                                *
 *                                                                      *
 *      TODO: fully implement                                           */

LOCAL   void	surfaces_on_bonds_of_curve(
	INTERFACE 	*intfc)
{
  	BOND 		*b;
	CURVE 		*c;
	SURFACE 	*s;
	BOND_TRI	*t;
	int 		j, i = 0;

	(void) next_bond(intfc,NULL,NULL);
	(void) printf(" bond#   bondX  curve#   curveX   curve surfaces     "
		      "::      bond surfaces\n");
	while (next_bond(intfc,&b,&c))
	{
	    (void) printf("%d  %llu %5d  %llu   : ",i,bond_number(b,intfc),
			  index_of_pointer((POINTER*)intfc->curves,c),
			  curve_number(c));
	    j = 0; 
	    while ((c)->pos_surfaces != NULL &&
		   (s = (c)->pos_surfaces[j++]) != NULL)
	        (void) printf(" [+]%llu ",surface_number(s));
	    j = 0;
	    while ((c)->neg_surfaces != NULL &&
		   (s = (c)->neg_surfaces[j++]) != NULL)
	        (void) printf(" [-]%llu ",surface_number(s));
	    (void) printf(" :: ");
	    j = 0;
	    while (Btris(b) != NULL && (t = Btris(b)[j++]) != NULL)
	        (void) printf(" %llu ",
			      surface_number(Surface_of_tri(t->tri)));
	    (void) printf("\n");
	    ++i;
	}
} 		/*end surfaces_on_bonds_of_curve*/

/*
*			map_rectangle_to_circle():
*
*	Defines a planar transformation with the following properties.
*
*	1. Maps concentric rectangles inside the rectangle R0 of center m and
*          half side lengths X0 and Y0 to concentric circles centered at m
*          inside the circle C0 of radius r0.  R0 is mapped to C0.
*
*	2.  For region outside R0 the mapping transforms continuously from
*           the above mapping to the identity mapping. This is accomplished by
*           linearly interpolating between the map in (1) and the identity
*           mapping in the region outside R0 and inside Rmax,  where Rmax
*           is the rectangle concentric and similar to R0 with half diagonal
*	    length lmax.
*
*	3. Outside of Rmax the mapping is an identity.
*
*	Input:
*	        c    = coordinates to be transformed.
*               m    = center of the transformation.
*               X0   = half X width of the rectangle R0.
*               Y0   = half Y height of the rectangle R0.
*	        r0   = radius of image of the region defined in 1.
*	        lmax = half diagonal length of the rectangle at which the
*                      tranformation becomes the identity. 
*
*	Output:
*	        c   = transformed coordinates.
*
*/

LOCAL   void 	map_rectangle_to_circle( 
	double *c,
	double *m,
	double X0,
	double Y0,
	double r0,
	double lmax)
{
	double        r, g, alpha;
	double        x, y;
	double        ax, ay;
	double        l, l0;

	l0 = hypot(X0,Y0);

  	x = c[0] - m[0];
	y = c[1] - m[1];
	r = hypot(x,y);

	ax = fabs(x)/X0;
	ay = fabs(y)/Y0;
	l = l0*max(ax,ay);

	g = l*r0/l0;
	if (l <= l0)
	    alpha = (r > 0.0) ? g/r : 1.0;
	else if ((l0 < l) && (l <= lmax))
	{
	    double t = (l - l0)/(lmax - l0);
	    alpha = (1.0 - t)*g/r + t;
	}
	else
	    alpha = 1.0;

	c[0] = alpha*x + m[0];
	c[1] = alpha*y + m[1];

  	return;
}    		/*end map_rectangle_to_circle*/

LOCAL   double jet_inner_func(
        POINTER func_params,
        double *coords)
{
        JET_PARAMS *params;
        double *cen,r0,r1,r2,h1,h2,h3,ratio,arg;
        double tmp;

        params = (JET_PARAMS *)func_params;
        cen = params->cen;
        r1 = params->r1;
        r2 = params->r2;
        ratio = params->ratio;

        h1 = params->cen[2] - params->h[0];
        h2 = h1 - params->h[1];
        h3 = h2 - params->h[2]*(1.0-ratio);
        if (ratio < 0)
            h3 = h2 - params->h[2];
	
        arg = 1.0;
        if (coords[2] >= h1)
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r1);
        else if (coords[2] >= h3)
        {
            r0 = sqr(r1)-sqr(coords[2]-h1);
            if (coords[2] >= h2)
                arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-r0;
            else
                arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r2);
        }
        if ((ratio < 0.0) && (coords[2] < h3))
        {
            double a2,b2,ytmp;
            a2 = sqr(1.2*r2);
            ytmp = h3 - 1.*r2;
            b2 = sqr(h3-ytmp)/(a2-sqr(r2-cen[0]));
            r0 = a2 - sqr(coords[2]-ytmp)/b2;
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-r0;
        }
                                                                                                    
        return arg;
}       /* end jet_inner_func */

LOCAL   double jet_outer_func(
        POINTER func_params,
        double *coords)
{
        JET_PARAMS *params;
        double *cen,r0,r1,r2,h1,h2,h3,ratio,arg;
	
        params = (JET_PARAMS *)func_params;
        cen = params->cen;
        r1 = params->r1;
        r2 = params->r2;
        ratio = params->ratio;
        h1 = params->cen[2] - params->h[0];
        h2 = h1 - params->h[1] - params->h[2]*(1.0-ratio);
        h3 = h2 - params->h[2]*ratio;
        if (ratio < 0)
            h2 = h3 = h1 - params->h[1] - params->h[2];

        arg = 1.0;
        if (coords[2] >= h1)
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r1);
        else if (coords[2] >= h3)
        {
            r0 = sqr(r1) - sqr(coords[2]-h1);
            if (ratio > 0.0)
            {
                if (coords[2] > h2)
                    arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-r0;
                else if ((sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])) <= r0)
                    arg = -sqr(coords[0]-cen[0])-sqr(coords[1]-cen[1])+sqr(r2);
            }
            else
                arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-r0;
        }
        return arg;
}       /* end jet_outer_func */

LOCAL   double jet_nozzle_func(
        POINTER func_params,
        double *coords)
{
        JET_PARAMS *params;
        double *cen,r1,r2,h1,h2,ratio,arg;
                                                                                                    
        params = (JET_PARAMS *)func_params;
        cen = params->cen;
        r1 = params->r1;
        r2 = params->r2;
        ratio = params->ratio;
                                                                                                    
        h1 = params->h[0] + params->h[1];
        h2 = params->h[0] + params->h[1]*ratio;
                                                                                                    
        arg = 1.0;
                                                                                                    
        if (coords[2] >= h1)
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r1);
        else if (coords[2] >= h2)
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r2);
        if ((ratio == 0) && (coords[2] < h2))
        {
            double a2,b2,r0,ytmp;
            a2 = sqr(1.5*r2);
            ytmp = h2 - 1.5*r2;
            b2 = sqr(h2-ytmp)/(a2-sqr(r2-cen[0]));
            r0 = a2 - sqr(coords[2]-ytmp)/b2;
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-r0;
        }
                                                                                                    
        return arg;
}       /* end jet_nozzle_func */
                                                                                                    
LOCAL   double jet_plane_func(
        POINTER func_params,
        double *coords)
{
        JET_PARAMS *params;
        double *cen,r1,r2,h1,ratio,arg;
                                                                                                    
        params = (JET_PARAMS *)func_params;
        cen = params->cen;
        r1 = params->r1;
        r2 = params->r2;
        ratio = params->ratio;
                                                                                                    
        h1 = params->h[0];
                                                                                                    
        arg = 1.0;
        if ((coords[2] > h1) && ((sqr(coords[0]-cen[0])+
                sqr(coords[1]-cen[1])) >= sqr(r2)))
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r1);
                                                                                                    
        if ((coords[2] + coords[0] <= 0.595) &&
            (coords[2] + coords[0] >= 0.545) &&
            (coords[0] < (0.4 + 0.1*cos(PI/4.0))))
            arg = -1.0;
                                                                                                    
        return arg;
}       /* end jet_plane_func */


LOCAL   double jet_top_bdry_func(
        POINTER func_params,
        double *coords)
{
	BDRY_BOX_PARAMS *params;
        double *U;
	        
	params = (BDRY_BOX_PARAMS *)func_params;
	U = params->U;
        if (coords[2] > U[2])
            return 1.0;
	else
	    return -1.0;	
}       /* end jet_top_bdry_func */

LOCAL   double jet_bottom_bdry_func(
        POINTER func_params,
        double *coords)
{
        BDRY_BOX_PARAMS *params;
        double *L;
		
	params = (BDRY_BOX_PARAMS *)func_params;
	L = params->L;
	
        if (coords[2] < L[2])
            return -1.0;
        else
            return 1.0;
}       /* end jet_bottom_bdry_func */


EXPORT boolean make_jet_surfaces(
		RECT_GRID   *rgr,
		COMPONENT   comp0,
		COMPONENT   comp1,
		COMPONENT   comp2,
		double       (*func1)(POINTER,double*),
		POINTER     func1_params,
		double       (*func2)(POINTER,double*),
		POINTER     func2_params,
		double       (*func3)(POINTER,double*),
		POINTER     func3_params,
		double       (*func4)(POINTER,double*),
		POINTER     func4_params,
		SURFACE     **s,
		CURVE       **c,
		int         *num_surfs,
		int         *num_curves)
{
	/* num_crx        all crx in computational domain
	 * num_crx0       crx between comp0 and comp1
	 * num_crx1       crx between comp0 and comp2 plus num_crx0
	 * num_crx2       crx between comp1 and comp2 plus num_crx1
	 * num_crx3       crx between comp0 and NO_COMP plus num_crx2
	 * num_curve_crx  crx on the grid face
	 */

	int		i,j,num_crx,num_crx0,num_crx1,num_crx2,
	                num_crx3,num_crx4,num_crx5,
	                num_crx6,*gmax,num_curve_crx,num_crx7;
	RECT_GRID	dual_gr;
	SURFACE		**surfs;
	EG_CRX		Eg_crx;
	BLK_INFO	blk_info;
	static BLK_CRX  *blk_crx;
	CURVE           **curve;
	NODE		*ns,*ne;
	POINT           ***p;
	int             **curve_crx;
	int             curve_n[3];
	int		is = 0;
	int             ic = 0;

	if (blk_crx == NULL) 
	{
	    blk_crx = alloc_blk_crx(NO);
	}
	for(i=0;i<3;i++)
	    curve_n[i]=0;
        
	zero_scalar(&Eg_crx,sizeof(EG_CRX));
	set_grid_for_surface_construction(&dual_gr,rgr);
	uni_array(&surfs,7,sizeof(SURFACE*));
	uni_array(&curve,3,sizeof(CURVE*));
	bi_array(&p,3,2,sizeof(POINT*));
	gmax = dual_gr.gmax;
	tri_array(&Eg_crx.comp,gmax[0]+1,gmax[1]+1,gmax[2]+1,
			sizeof(COMPONENT));
        
	reset_domain_comp(Eg_crx.comp,dual_gr);

	assign_negative_comp(func1,func1_params,Eg_crx.comp,
			dual_gr,comp0);
       
	assign_intersection_comp(func1,func1_params,func2,func2_params,
			Eg_crx.comp,dual_gr,comp1,
			POSITIVE_SIDE,NEGATIVE_SIDE);
	assign_positive_comp(func2,func2_params,Eg_crx.comp,
			dual_gr,comp2);
	
	    /* set exterior comp for make bdry surface */
	assign_exterior_comp(Eg_crx.comp,dual_gr);

	num_crx = count_crx_through_comp(gmax,Eg_crx.comp);
	printf("num_crx = %d\n",num_crx);
	alloc_grid_crx_mem(&Eg_crx,gmax,num_crx,YES);
        
	num_crx0 = install_grid_crx(func1,func1_params,&Eg_crx,dual_gr,
				comp0,comp1);
	num_crx1 = install_grid_crx(func2,func2_params,&Eg_crx,dual_gr,
				comp1,comp2);
	num_crx2 = install_grid_crx2(func1,func1_params,&Eg_crx,dual_gr,
				comp2,comp0,num_crx1);
        num_crx3 = install_grid_topbdry_crx(func3,func3_params,&Eg_crx,dual_gr,
			        comp0,NO_COMP,num_crx2);
        num_crx4 = install_grid_topbdry_crx(func3,func3_params,&Eg_crx,dual_gr,
                                comp1,NO_COMP,num_crx3);
	num_crx5 = install_grid_topbdry_crx(func3,func3_params,&Eg_crx,dual_gr,
			        comp2,NO_COMP,num_crx4);
	num_crx6 = install_grid_bottombdry_crx(func4,func4_params,&Eg_crx,
				dual_gr,comp2,NO_COMP,num_crx5);
        bi_array(&curve_crx,3,num_crx-num_crx6,sizeof(int));
	num_crx7 = num_crx6;
	num_curve_crx = install_curve_crx(func1,func1_params,func2,func2_params,
		               	&Eg_crx, dual_gr,&num_crx7,p,comp0,comp1,comp2,
				curve_crx,curve_n);

	/* make surfaces and make curves*/
        
	/* make the inwall and the contact curve 
           and the curve between inwall and the top bdry */
	if (num_crx0 != 0)
	{
	    surfs[is] = make_surface(comp1,comp0,(CURVE**)NULL,(CURVE**)NULL);
            for (i = 0; i < num_crx0; ++i)
	        Eg_crx.crx_store[i].s = surfs[is];
	    surfs[is]->number = 1;
	    if(curve_n[0] != 0)
	    {
		if(pp_numnodes() == 0)
		{
		    ns = ne = make_node(p[0][0]);	
		    curve[ic] = make_curve(NO_COMP,NO_COMP,ns,ne);
		    curve[ic]->start->posn = p[0][0];
		    curve[ic]->end->posn = p[0][0];
		}
		else
		{
                    ns =  make_node(p[0][0]);
                    ne = make_node(p[0][1]);
                    curve[ic] = make_curve(NO_COMP,NO_COMP,ns,ne);
                    curve[ic]->start->posn = p[0][0];
                    curve[ic]->end->posn = p[0][1];
		}
                curve[ic]->first = curve[ic]->last = NULL;
                add_to_pointers((POINTER)curve[ic],
                               (POINTER**)&surfs[is]->pos_curves);
                add_to_pointers((POINTER)surfs[is],
                               (POINTER**)&curve[ic]->pos_surfaces);
		curve[ic]->number = 1;
                ic++;
	    }
            if(curve_n[1] != 0)
            {
                if(pp_numnodes() == 0)
                {
                    ns = ne = make_node(p[1][0]);
                    curve[ic] = make_curve(NO_COMP,NO_COMP,ns,ne);
                    curve[ic]->start->posn = p[1][0];
                    curve[ic]->end->posn = p[1][0];
                }
		else
                {    
                    ns =  make_node(p[1][0]);
                    ne = make_node(p[1][1]);
                    curve[ic] = make_curve(NO_COMP,NO_COMP,ns,ne);
                    curve[ic]->start->posn = p[1][0];
                    curve[ic]->end->posn = p[1][1];
                } 
                curve[ic]->first = curve[ic]->last = NULL;
                add_to_pointers((POINTER)curve[ic],
                               (POINTER**)&surfs[is]->pos_curves);
                add_to_pointers((POINTER)surfs[is],
                               (POINTER**)&curve[ic]->pos_surfaces);
		curve[ic]->number = 2;
                ic++;
            }
            is++;
	}
	/* make the outwall and the curve between outwall and top bdry */
	if ((num_crx1-num_crx0) != 0)
	{
	    surfs[is] = make_surface(comp2,comp1,(CURVE**)NULL,(CURVE**)NULL);
	    for (i = num_crx0; i < num_crx1; ++i)
	        Eg_crx.crx_store[i].s = surfs[is];
	    surfs[is]->number = 2;
            if(curve_n[0] != 0)
            {
                add_to_pointers((POINTER)curve[0],
                               (POINTER**)&surfs[is]->pos_curves);
                add_to_pointers((POINTER)surfs[is],
                               (POINTER**)&curve[0]->pos_surfaces);
            }
            if(curve_n[2] != 0)
            {
                if(pp_numnodes() == 0)
                {
                    ns = ne = make_node(p[2][0]);
                    curve[ic] = make_curve(NO_COMP,NO_COMP,ns,ne);
                    curve[ic]->start->posn = p[2][0];
                    curve[ic]->end->posn = p[2][0];
                }
                else
                {
                    ns =  make_node(p[2][0]);
                    ne = make_node(p[2][1]);
                    curve[ic] = make_curve(NO_COMP,NO_COMP,ns,ne);
                    curve[ic]->start->posn = p[2][0];
                    curve[ic]->end->posn = p[2][1];
                }
                curve[ic]->first = curve[ic]->last = NULL;
                add_to_pointers((POINTER)curve[ic],
                               (POINTER**)&surfs[is]->pos_curves);
                add_to_pointers((POINTER)surfs[is],
                               (POINTER**)&curve[ic]->pos_surfaces);
		curve[ic]->number = 3;
                ic++;
            }
            is++;
	}
	/* make the Contact surface */
	if ((num_crx2-num_crx1) != 0)
	{
	    surfs[is] = make_surface(comp0,comp2,(CURVE**)NULL,(CURVE**)NULL);
	    for (i = num_crx1; i < num_crx2; ++i)
	        Eg_crx.crx_store[i].s = surfs[is];
	    surfs[is]->number = 3;
            if(curve_n[0] != 0)
            {
                add_to_pointers((POINTER)curve[0],
                               (POINTER**)&surfs[is]->pos_curves);
                add_to_pointers((POINTER)surfs[is],
                               (POINTER**)&curve[0]->pos_surfaces);
            }
            is++;       
	}
        /* make the surface between the fuel container and the top bdry */ 
	if ((num_crx3-num_crx2) != 0)
        {
            surfs[is] = make_surface(NO_COMP,comp0,(CURVE**)NULL,(CURVE**)NULL);
            for (i = num_crx2; i < num_crx3; ++i)
                Eg_crx.crx_store[i].s = surfs[is];
	    surfs[is]->number = 4;
            if(curve_n[1] != 0)
            {
                if(curve_n[0] != 0)
                {
                    add_to_pointers((POINTER)curve[1],
                                   (POINTER**)&surfs[is]->neg_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[1]->neg_surfaces);
                }
                else
                {
                    add_to_pointers((POINTER)curve[0],
                                   (POINTER**)&surfs[is]->neg_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[0]->neg_surfaces);
                }
            }
            is++;
        }
        /* make the surface between the ambient and the top bdry */
        if ((num_crx4-num_crx3) != 0)
        {
            surfs[is] = make_surface(NO_COMP,comp1,(CURVE**)NULL,(CURVE**)NULL); 
	    for (i = num_crx3; i < num_crx4; ++i)
                Eg_crx.crx_store[i].s = surfs[is];
	    surfs[is]->number = 5;
            if(curve_n[1] != 0)
            {
                if(curve_n[0] != 0)
                {
                    add_to_pointers((POINTER)curve[1],
                                   (POINTER**)&surfs[is]->pos_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[1]->pos_surfaces);
                }
                else
                {
                    add_to_pointers((POINTER)curve[0],
                                   (POINTER**)&surfs[is]->pos_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[0]->pos_surfaces);
                }
            }
            if(curve_n[2] != 0)
            {
                if(ic == 3)
                {
                    add_to_pointers((POINTER)curve[2],
                                   (POINTER**)&surfs[is]->neg_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[2]->neg_surfaces);
                }
                else if(ic == 2)
                {
                    add_to_pointers((POINTER)curve[1],
                                   (POINTER**)&surfs[is]->neg_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[1]->neg_surfaces);
                }
                else
                {
                    add_to_pointers((POINTER)curve[0],
                                   (POINTER**)&surfs[is]->neg_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[0]->neg_surfaces);
                }
            }
            is++;
        }
        /* make the bdry surface  on top bdry */
	if ((num_crx5-num_crx4) != 0)
        {
            surfs[is] = make_surface(NO_COMP,comp2,(CURVE**)NULL,(CURVE**)NULL);
	    for (i = num_crx4; i < num_crx5; ++i)
                Eg_crx.crx_store[i].s = surfs[is];
	    surfs[is]->number = 6;
            if(curve_n[2] != 0)
            {
                if(ic == 3)
                { 
                    add_to_pointers((POINTER)curve[2],
                                   (POINTER**)&surfs[is]->pos_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[2]->pos_surfaces);
                }
                else if(ic == 2)
                {
                    add_to_pointers((POINTER)curve[1],
                                   (POINTER**)&surfs[is]->pos_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[1]->pos_surfaces);
                }
                else
                {
                    add_to_pointers((POINTER)curve[0],
                                   (POINTER**)&surfs[is]->pos_curves);
                    add_to_pointers((POINTER)surfs[is],
                                   (POINTER**)&curve[0]->pos_surfaces);
                }
            }
            is++;
        }
        /* make the bdry surface on bottom bdry */
	if ((num_crx6-num_crx5) != 0)
	{
            surfs[is] = make_surface(NO_COMP,comp2,(CURVE**)NULL,(CURVE**)NULL);    
	    for (i = num_crx5; i < num_crx6; ++i)
	        Eg_crx.crx_store[i].s = surfs[is];
	    surfs[is]->number = 7;
	    is++;
	}
       	if (num_curve_crx != 0) 
	{
            Eg_crx.num_curves = blk_info.num_curves = ic;
	    uni_array(&Eg_crx.curves,Eg_crx.num_curves,sizeof(CURVE*));
            uni_array(&blk_info.curves,blk_info.num_curves,sizeof(CURVE*));
	    
	    for(i = 0; i < blk_info.num_curves; i++)
	    {	    
	        Eg_crx.curves[i] = curve[i];
                blk_info.curves[i] = curve[i];
	    }

	    for(i = 0; i < curve_n[0]; i++)
	        Eg_crx.crx_store[curve_crx[0][i]].c = curve[0];
	    for(i = 0; i < curve_n[1]; i++)
            {
                if(curve_n[0] != 0)
	            Eg_crx.crx_store[curve_crx[1][i]].c = curve[1];
                else
                    Eg_crx.crx_store[curve_crx[1][i]].c = curve[0];
            }
	    for(i = 0; i < curve_n[2]; i++)
            {
                if(ic == 3)
	            Eg_crx.crx_store[curve_crx[2][i]].c = curve[2];
                else if(ic == 2)
                    Eg_crx.crx_store[curve_crx[2][i]].c = curve[1];
                else
                    Eg_crx.crx_store[curve_crx[2][i]].c = curve[0];
            } 
	}
	blk_info.num_surfs = is;
	uni_array(&blk_info.surfs,blk_info.num_surfs,sizeof(SURFACE*));
	uni_array(&blk_info.cur_tris,blk_info.num_surfs,sizeof(TRI*));

	for (i = 0; i < blk_info.num_surfs; i++)
	{
	    first_tri(surfs[i]) = last_tri(surfs[i]) = NULL; 
	    surfs[i]->num_tri = 0;
	    blk_info.cur_tris[i] = NULL;
	    blk_info.surfs[i] = surfs[i];
	}
	
	blk_crx->comps[0] = comp0;
	blk_crx->comps[1] = comp1;
	blk_crx->comps[2] = comp2;
	blk_crx->comps[3] = NO_COMP;
	blk_crx->blk_info = &blk_info; 

	make_grid_surfaces(blk_crx,&Eg_crx,gmax,YES);

        for (i = 0; i < blk_info.num_surfs; ++i)
	{
	    last_tri(blk_info.surfs[i]) = blk_info.cur_tris[i];
	    last_tri(blk_info.surfs[i])->next = 
	                tail_of_tri_list(blk_info.surfs[i]);
	    first_tri(blk_info.surfs[i])->prev = 
	                head_of_tri_list(blk_info.surfs[i]);
	}

	if (num_curve_crx != 0)
        {
            for(i=0; i< blk_info.num_curves; i++)	
            {
	        curve[i]->last->next = curve[i]->first->prev = NULL;
	    }
        }
	for (i = 0; i < blk_info.num_curves; i++) 
	    reorder_curve_link_list(Eg_crx.curves[i]);
	for (i = 0; i < blk_info.num_surfs; ++i)
	{
	    reset_intfc_num_points(surfs[i]->interface);
	}

	free_grid_crx_mem(&Eg_crx,YES);
	free_these(3,Eg_crx.comp,blk_info.surfs,blk_info.cur_tris);
	if (num_curve_crx != 0) 
	    free_these(1,blk_info.curves);

	*num_surfs = is;
	*num_curves = ic;
	for (i = 0; i < blk_info.num_surfs; i++)
	    s[i] = surfs[i];
	for (i = 0; i < blk_info.num_curves; i++)
        {
	    c[i] = curve[i];
        }

	free_these(3,surfs,curve,p);

	return YES;
}	/* end make_jet_surfaces */
#endif /* defined(THREED) */

#if defined(TWOD)
/*
*		insert_linear_section():
*
*	Inserts a linear section running from start to end into the bond b.
*/

LOCAL	void	insert_linear_section(
	BOND  *b,
	CURVE *c,
	double *start,
	double *end,
	Front *fr)
{
	RECT_GRID *gr = fr->rect_grid;
	double     *h = gr->h;
	double     L, ds;
	double     coords[2];
	double     space;
	int       i, N;
	int       dim = gr->dim;
	double     tol = MIN_SC_SEP(fr->interf);

	space = Front_spacing(fr,
	    	    (wave_type(c) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE) ?
	    		VECTOR_WAVE : GENERAL_WAVE);

	L = _scaled_separation(start,end,h,dim);
	N = irint(L/space);
	ds = L/N;

	if (_scaled_separation(start,Coords(b->start),h,dim) > tol)
	{
	    insert_point_in_bond(Point(start),b,c);
	    b = b->next;
	}
	for (i = 1; i < N; ++i)
	{
	    double alpha;
	    alpha = i*ds/L;
	    coords[0] = (1.0-alpha)*start[0] + alpha*end[0];
	    coords[1] = (1.0-alpha)*start[1] + alpha*end[1];
	    insert_point_in_bond(Point(coords),b,c);
	    b = b->next;
	}
	if (_scaled_separation(end,Coords(b->end),h,dim) > tol)
	{
	    insert_point_in_bond(Point(end),b,c);
	    b = b->next;
	}
}		/*end insert_linear_section*/

/*
*		insert_elliptical_curve_section():
*
*	Inserts an elliptical section with center cen and axes with
*	radii rad rotated by angle phi, into the bond b of curve c.
*	The inserted points are chosen to be separated by appropriate
*	front spacing factor in the scaled metric.
*/

LOCAL	const double *el_rad, *el_h;
LOCAL	double el_phi;
LOCAL	double el_len(double,POINTER);
#if defined(__cplusplus)
extern "C" {    
#endif /* defined(__cplusplus) */
    LOCAL   LSODE_FUNC  iel_len;
#if defined(__cplusplus)
}                          
#endif /* defined(__cplusplus) */

LOCAL	void	insert_elliptical_curve_section(
	BOND  *b,
	CURVE *c,
	double ThetaS,
	double ThetaE,
	double *cen,
	double *rad,
	double phi,
	Front *fr)
{
	RECT_GRID         *gr = fr->rect_grid;

#if defined(USE_OVERTURE)
        double             h[MAXD];  
        int               amr_refinecoeff = 1;
#else /* if defined(USE_OVERTURE) */ 
	double             *h = gr->h;
#endif /* if defined(USE_OVERTURE) */

	double             L;
	double             epsabs, epsrel;
	double             abserr;
	double             space;
	double             *theta;
	double             ds;
	double             cp, sp, coords[2];
	double             x, y;
	int               i, N, dim = gr->dim;
	int               neval;
	int               istate;
	QUADRATURE_STATUS ier;
	double             tol = MIN_SC_SEP(fr->interf);

#if defined(USE_OVERTURE)
        if(fr->NumberOfLevels != 0)
            amr_refinecoeff = (int)pow(2, fr->NumberOfLevels-1);
        for(i = 0; i < dim; i++)
            h[i] = gr->h[i]/amr_refinecoeff;
#endif /* if defined(USE_OVERTURE) */

	space = Front_spacing(fr,
	    	    (wave_type(c) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE) ?
	    		VECTOR_WAVE : GENERAL_WAVE);

	el_rad = rad; el_h = h;
	el_phi = phi;

	/* Compute the arc length in the scaled metric of the curve
	 * section to be inserted.
	 */
	epsrel = 1.0e-6;
	epsabs = sqrt((rad[0]/h[0])*(rad[1]/h[1]))*epsrel;
	L = dqng(el_len,NULL,ThetaS,ThetaE,
		      epsabs,epsrel,&abserr,&neval,&ier);
	switch (ier)
	{
	case INVALID_EPSILON:
	    screen("ERROR in insert_elliptical_curve_section(), "
		   "invalid epsilons\n"
		   "epsabs = %"FFMT", epsrel = %"FFMT"\n",epsabs,epsrel);
	    clean_up(ERROR);
	    break;
	case INACCURATE_INTEGRAL:
	    if ((fabs(abserr) > 20.0*epsabs) && (fabs(abserr) > 20.0*epsrel*L))
	    {
		/*
		 * Don't print a warning is we are close to satifying the
		 * error requirements
		 */
	        (void) printf("WARNING in insert_elliptical_curve_section(), "
			      "inaccurate result\n \tneval = %d, "
			      "abserr = %"FFMT" result = %"FFMT"\n",
			      neval,abserr,L);
	    }
	    break;
	case ACCURATE_INTEGRAL:
	default:
	    break;
	}
	N = irint(fabs(L)/space);
	if (N < 1)	/*Nothing to do, the elliptical section would be */
	    return;	/* too short */

	uni_array(&theta,N+1,FLOAT);

	ds = L/N;

	/*
	*  Integrate with respect to arc length to find angles with
	*  eqi-distance separation in the scaled metric along the ellipse
	*/
	theta[0] = ThetaS;
	istate = 1;
	if ((!ode_solver(iel_len,theta,0.0,ds,N,epsrel,
			 epsrel*fabs(ThetaE-ThetaS),
		         abserr,&istate,NULL)) || (istate != 2))
	{
	    screen("ERROR in insert_elliptical_curve_section() can't solve ODE\n"
	           "error flag = %d\n",istate);
	    clean_up(ERROR);
	}
	theta[N] = ThetaE;


	cp = cos(phi); sp = sin(phi);
	x = rad[0]*cos(theta[0]-phi);
	y = rad[1]*sin(theta[0]-phi);
	coords[0] = cen[0] + x*cp + y*sp;
	coords[1] = cen[1] - x*sp + y*cp;
	if (_scaled_separation(coords,Coords(b->start),h,dim) > tol)
	{
	    insert_point_in_bond(Point(coords),b,c);
	    b = b->next;
	}
	for (i = 1; i < N; ++i)
	{
	    x = rad[0]*cos(theta[i]-phi);
	    y = rad[1]*sin(theta[i]-phi);
	    coords[0] = cen[0] + x*cp + y*sp;
	    coords[1] = cen[1] - x*sp + y*cp;
	    insert_point_in_bond(Point(coords),b,c);
	    b = b->next;
	}
	x = rad[0]*cos(theta[N]-phi);
	y = rad[1]*sin(theta[N]-phi);
	coords[0] = cen[0] + x*cp + y*sp;
	coords[1] = cen[1] - x*sp + y*cp;
	if (_scaled_separation(coords,Coords(b->end),h,dim) > tol)
	{
	    insert_point_in_bond(Point(coords),b,c);
	    b = b->next;
	}

	free(theta);
}		/*end insert_elliptical_curve_section*/

/*ARGSUSED*/
LOCAL	double el_len(
	double theta,
	POINTER prms)
{
	IMPORT const double *el_rad, *el_h;
	IMPORT double el_phi;
	double dx, dy;
	double dx1, dy1;
	double c, s;

	dx1 = -el_rad[0]*sin(theta-el_phi);
	dy1 =  el_rad[1]*cos(theta-el_phi);
	c = cos(el_phi);
	s = sin(el_phi);
	dx =  dx1*c + dy1*s;
	dy = -dx1*s + dy1*c;
	return hypot(dx/el_h[0],dy/el_h[1]);
}		/*end el_len*/

#if defined(__cplusplus)
extern "C" {    
#endif /* defined(__cplusplus) */

/*ARGSUSED*/
LOCAL	void iel_len(
	int		*neq,
	double		*x,
	double		*y,
	double		*yp)
{
	*yp = 1.0/el_len(*y,NULL);
}		/*end iel_len*/

#if defined(__cplusplus)
}                          
#endif /* defined(__cplusplus) */

/*
*		insert_polynomial_curve_section():
*
*	Inserts a polynomial curve y = p[x] into the curve c inside bond b.
*/

LOCAL	const double *pc, *poly_h;
LOCAL	int   d;

LOCAL	double poly_len(double,POINTER);
#if defined(__cplusplus)
extern "C" {    
#endif /* defined(__cplusplus) */
    LOCAL   LSODE_FUNC  ipoly_len;
#if defined(__cplusplus)
}                          
#endif /* defined(__cplusplus) */

LOCAL	void	insert_polynomial_curve_section(
	BOND  *b,
	CURVE *c,
	double xs,
	double xe,
	double *coef,
	int   n,
	Front *fr)
{
	RECT_GRID         *gr = fr->rect_grid;

#if defined(USE_OVERTURE)
        double             h[MAXD];
        int               amr_refinecoeff = 1;
#else /* if defined(USE_OVERTURE) */
        double             *h = gr->h;
#endif /* if defined(USE_OVERTURE) */

	double             L;
	double             epsabs, epsrel;
	double             abserr;
	double             space;
	double             *x;
	double             ds;
	double             coords[2];
	int               i, N, dim = gr->dim;
	int               neval;
	int               istate;
	QUADRATURE_STATUS ier;
	double             tol = MIN_SC_SEP(fr->interf);

#if defined(USE_OVERTURE)
        if(fr->NumberOfLevels != 0)
            amr_refinecoeff = (int)pow(2, fr->NumberOfLevels-1);
        for(i = 0; i < dim; i++)
            h[i] = gr->h[i]/amr_refinecoeff;
#endif /* if defined(USE_OVERTURE) */

	space = Front_spacing(fr,
	    	    (wave_type(c) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE) ?
	    		VECTOR_WAVE : GENERAL_WAVE);

	pc = coef;
	poly_h = h;
	d = n;

	/* Compute the arc length in the scaled metric of the curve
	 * section to be inserted.
	 */
	epsrel = 1.0e-6;
	epsabs = fabs(xs - xe)*epsrel;
	L = dqng(poly_len,NULL,xs,xe,
		      epsabs,epsrel,&abserr,&neval,&ier);
	switch (ier)
	{
	case INVALID_EPSILON:
	    screen("ERROR in insert_polynomial_curve_section(), "
		   "invalid epsilons\n"
		   "epsabs = %"FFMT", epsrel = %"FFMT"\n",epsabs,epsrel);
	    clean_up(ERROR);
	    break;
	case INACCURATE_INTEGRAL:
	    if ((fabs(abserr) > 20.0*epsabs) && (fabs(abserr) > 20.0*epsrel*L))
	    {
		/*
		 * Don't print a warning is we are close to satifying the
		 * error requirements
		 */
	        (void) printf("WARNING in insert_polynomial_curve_section(), "
			      "inaccurate result\n \tneval = %d, "
			      "abserr = %"FFMT" result = %"FFMT"\n",
			      neval,abserr,L);
	    }
	    break;
	case ACCURATE_INTEGRAL:
	default:
	    break;
	}
	N = irint(fabs(L)/space);
	if (N < 1)	/*Nothing to do, the elliptical section would be */
	    return;	/* too short */

	uni_array(&x,N+1,FLOAT);

	ds = L/N;

	/*
	*  Integrate with respect to arc length to find angles with
	*  eqi-distance separation in the scaled metric along the ellipse
	*/
	x[0] = xs;
	istate = 1;
	if ((!ode_solver(ipoly_len,x,xs,ds,N,epsrel,
			 epsrel*fabs(xs-xe),
		         abserr,&istate,NULL)) || (istate != 2))
	{
	    screen("ERROR in insert_polynomial_curve_section() can't solve ODE\n"
	           "error flag = %d\n",istate);
	    clean_up(ERROR);
	}
	x[N] = xe;


	coords[0] = xs;
	coords[1] = eval_poly(coords[0],coef,n);
	if (_scaled_separation(coords,Coords(b->start),h,dim) > tol)
	{
	    insert_point_in_bond(Point(coords),b,c);
	    b = b->next;
	}
	for (i = 1; i < N; ++i)
	{
	    coords[0] = x[i];
	    coords[1] = eval_poly(coords[0],coef,n);
	    insert_point_in_bond(Point(coords),b,c);
	    b = b->next;
	}
	coords[0] = x[N];
	coords[1] = eval_poly(coords[0],coef,n);
	if (_scaled_separation(coords,Coords(b->end),h,dim) > tol)
	{
	    insert_point_in_bond(Point(coords),b,c);
	    b = b->next;
	}

	free(x);
}		/*end insert_polynomial_curve_section*/

/*ARGSUSED*/
LOCAL	double poly_len(
	double x,
	POINTER prms)
{
	IMPORT	const double *pc, *poly_h;
	double dp;
	int   i;

	dp = 0.0;
	for (i=d-1; i >= 0; --i)
	    dp = (i+1)*pc[i+1] + x*dp;

	return hypot(1.0/poly_h[0],dp/poly_h[1]);
}		/*end poly_len*/

#if defined(__cplusplus)
extern "C" {    
#endif /* defined(__cplusplus) */

/*ARGSUSED*/
LOCAL	void ipoly_len(
	int		*neq,
	double		*x,
	double		*y,
	double		*yp)
{
	*yp = 1.0/poly_len(*y,NULL);
}		/*end ipoly_len*/
#if defined(__cplusplus)
}                          
#endif /* defined(__cplusplus) */

LOCAL	double eval_poly(
	double       x,
	const double *coef,
	int         n)
{
    	double p;
	int   i;

	p = 0.0;
	for (i=n; i >= 0; --i)
	    p = coef[i] + x*p;
	return p;
}

#endif /* defined(TWOD) */

#endif /* defined(TWOD) || defined(THREED) */


#if defined(THREED)

LOCAL   int install_grid_crx2(
        double (*func)(POINTER,double*),
        POINTER func_params,
        EG_CRX *eg_crx,
        RECT_GRID grid,
        COMPONENT comp0,
        COMPONENT comp1,
	int       num_crx)
{
        double coords1[3];
        double coords2[3];
        double crds_crx[3];
        double *L = grid.L;
        double *h = grid.h;
        int *gmax = grid.gmax;
        int dim = grid.dim;
        int i,j,k;
        int n_crx = num_crx;
        BBI_POINT ****x_crx = eg_crx->x_crx;
        BBI_POINT ****y_crx = eg_crx->y_crx;
        BBI_POINT ****z_crx = eg_crx->z_crx;
        BBI_POINT *crx_store = eg_crx->crx_store;
        COMPONENT ***comp = eg_crx->comp;
        /*TMP */
        JET_PARAMS *params;
        double *cen,h1,h2,h3,ratio;

	int debug_flag = NO;

	/* install x-crossings */

	for (j = 0; j <= gmax[1]; ++j)
	{
	    coords1[1] = coords2[1] = L[1] + j*h[1];
	    for (k = 0; k <= gmax[2]; ++k)
	    {
		coords1[2] = coords2[2] = L[2] + k*h[2];
		for (i = 0; i < gmax[0]; ++i)
		{

		    if (((comp[i][j][k] == comp0) && (comp[i+1][j][k] == comp1)) ||
		        ((comp[i][j][k] == comp1) && (comp[i+1][j][k] == comp0)))
		    {
			coords1[0] = L[0] + i*h[0];
		        coords2[0] = L[0] + (i+1)*h[0];

		 	if (! grid_line_crx_in_dir(func,func_params,
			    	dim,coords1,coords2,crds_crx,0))
			{
			    screen("ERROR: in install_grid_crx(), no x-crxing!");
			    clean_up(ERROR);
			}
			if (crds_crx[0] - coords1[0] < 0.004*h[0])
			    crds_crx[0] = coords1[0] + 0.004*h[0];
			if (coords2[0] - crds_crx[0] < 0.004*h[0])
			    crds_crx[0] = coords2[0] - 0.004*h[0];
			crx_store[n_crx].p = Point(crds_crx);
			x_crx[i][j][k] = &crx_store[n_crx++];
		    }
		}
	    }
	}

	/* install y-crossings */

	for (i = 0; i <= gmax[0]; ++i)
	{
	    coords1[0] = coords2[0] = L[0] + i*h[0];
	    for (k = 0; k <= gmax[2]; ++k)
	    {
		coords1[2] = coords2[2] = L[2] + k*h[2];
		for (j = 0; j < gmax[1]; ++j)
		{
	            /* TMP		 */
		    /*
		    if(i == 0 && j == 13 && k == 53)
		    {
		        debug_flag == YES;
	                printf("test y-crossings comp0 %d, comp1 %d, i, j, k =[%d %d %d]\n", 
					comp0, comp1, i, j, k);		
		    }	    
		    if(i == 0 && j == 13 && k == 54)
		    {
		        debug_flag == YES;
	                printf("test y-crossings comp0 %d, comp1 %d, i, j, k =[%d %d %d]\n", 
					comp0, comp1, i, j, k);		
			printf("comp[i][j][k]= %d, comp[i][j+1][k] %d\n",
					comp[i][j][k], comp[i][j+1][k]);
		        coords1[1] = L[1] + j*h[1];
		        coords2[1] = L[1] + (j+1)*h[1];
			print_general_vector("coords1",coords1,3,"\n");
			print_general_vector("coords2",coords2,3,"\n");
		    }	    
		    */
		    /* END TMP */
		    if (((comp[i][j][k] == comp0) && (comp[i][j+1][k] == comp1)) ||
		        ((comp[i][j][k] == comp1) && (comp[i][j+1][k] == comp0)))
		    {
		        coords1[1] = L[1] + j*h[1];
		        coords2[1] = L[1] + (j+1)*h[1];

			/* TMP     */
			/*
			if(debug_flag == YES)
			{
			    printf("before grid_line_crx_in_dir\n");	
			    print_general_vector("coords1",coords1,3,"\n");
			    print_general_vector("coords2",coords2,3,"\n");
			}	
                        */
			if (! grid_line_crx_in_dir(func,func_params,
			    dim,coords1,coords2,crds_crx,1))
			{
			    screen("ERROR: in install_grid_crx(), no y-crxing!");
			    clean_up(ERROR);
			}
			{
			double cc1=crds_crx[1];
			if (crds_crx[1] - coords1[1] < 0.004*h[1])
			    crds_crx[1] = coords1[1] + 0.004*h[1];
			if (coords2[1] - crds_crx[1] < 0.004*h[1])
			    crds_crx[1] = coords2[1] - 0.004*h[1];
			}
			crx_store[n_crx].p = Point(crds_crx);
                        y_crx[i][j][k] = &crx_store[n_crx++];
		    }
		}
	    }
	}
	/* END NEW ADDED 050106 */

        /* install z-crossings */

        for (i = 0; i <= gmax[0]; ++i)
        {
            coords1[0] = coords2[0] = L[0] + i*h[0];
            for (j = 0; j <= gmax[1]; ++j)
            {
                coords1[1] = coords2[1] = L[1] + j*h[1];
                for (k = 0; k < gmax[2]; ++k)
                {
                    if (((comp[i][j][k] == comp0) && (comp[i][j][k+1] == comp1)) ||
                        ((comp[i][j][k] == comp1) && (comp[i][j][k+1] == comp0)))
                    {
                        coords1[2] = L[2] + k*h[2];
                        coords2[2] = L[2] + (k+1)*h[2];
                        if (! grid_line_crx_in_dir(func,func_params,
                            dim,coords1,coords2,crds_crx,2))
                        {
                           /*TMP */
                            printf("in install_grid_crx,i = %d,j =%d,j =%d\n",i,j,k);
                            screen("ERROR: in install_grid_crx(), no z-crxing!");
                            clean_up(ERROR);
                        }
                        if (crds_crx[2] - coords1[2] < 0.004*h[2])
                            crds_crx[2] = coords1[2] + 0.004*h[2];
                        if (coords2[2] - crds_crx[2] < 0.004*h[2])
                            crds_crx[2] = coords2[2] - 0.004*h[2];

                        crx_store[n_crx].p = Point(crds_crx);
                        z_crx[i][j][k] = &crx_store[n_crx++];
                    }
                }
            }
        }
        return n_crx;
}       /* end install_grid_crx */


LOCAL   int install_grid_topbdry_crx(
	double (*func)(POINTER,double*),
        POINTER func_params,
        EG_CRX *eg_crx,
        RECT_GRID grid,
        COMPONENT comp0,
        COMPONENT comp1,
	int  num_crx)
{   
        double coords1[3];
        double coords2[3];
        double crds_crx[3];
        double f1,f2;
        double *L = grid.L;
        double *U = grid.U;
        double *h = grid.h;
        int *gmax = grid.gmax;
        int dim = grid.dim;
        int i,j,k,m;
        int n_crx = num_crx;
	BBI_POINT ****x_crx = eg_crx->x_crx;
        BBI_POINT ****y_crx = eg_crx->y_crx;
        BBI_POINT ****z_crx = eg_crx->z_crx;
        BBI_POINT *crx_store = eg_crx->crx_store;
        COMPONENT ***comp = eg_crx->comp;

	/* install z-crossings */
	                                                                                
        for (i = 0; i <= gmax[0]; ++i)
        {
            coords1[0] = coords2[0] = L[0] + i*h[0];
            for (j = 0; j <= gmax[1]; ++j)
            {
                coords1[1] = coords2[1] = L[1] + j*h[1];
	        for (k = (int)(0.5*gmax[2]); k <= gmax[2]-1; ++k)
		{	
                    if (((comp[i][j][k] == comp0) && (comp[i][j][k+1] == comp1)) ||
                       ((comp[i][j][k] == comp1) && (comp[i][j][k+1] == comp0)))
	            {
		        coords1[2] = L[2] + k*h[2];
		        coords2[2] = L[2] + (k+1)*h[2];

                        f1 = (*func)(func_params,coords1);
                        f2 = (*func)(func_params,coords2);
                        if((f1 < 0 && f2 > 0) || (f1 > 0 && f2 < 0))
                        {
			    for(m =0; m <3; m++)
			        crds_crx[m] = (coords1[m] + coords2[m])/2.0;
                        }
                        else
                        {
                            screen("ERROR: in install_grid_topbdry_crx(), no z-crxing!");
                            (void)printf("comp0 =%d,comp1 = %d,i = %d,j =%d,k =%d\n",comp0,comp1,i,j,k);
                            clean_up(ERROR);
                        }
	                crx_store[n_crx].p = Point(crds_crx);
                        z_crx[i][j][k] = &crx_store[n_crx++];
		    }
		}
	    }
	}
	return n_crx;
}   /*end install_grid_topbdry_crx*/

LOCAL   int install_grid_bottombdry_crx(
        double (*func)(POINTER,double*),
        POINTER func_params,
        EG_CRX *eg_crx,
        RECT_GRID grid,
        COMPONENT comp0,
        COMPONENT comp1,
	int  num_crx)
{
	double coords1[3];
        double coords2[3];
        double crds_crx[3];
        double f1,f2;
        double *L = grid.L;
        double *U = grid.U;
        double *h = grid.h;
        int *gmax = grid.gmax;
        int dim = grid.dim;
        int i,j,k,m;
        int n_crx = num_crx;

	BBI_POINT ****x_crx = eg_crx->x_crx;
	BBI_POINT ****y_crx = eg_crx->y_crx;
	BBI_POINT ****z_crx = eg_crx->z_crx;
	BBI_POINT *crx_store = eg_crx->crx_store;
	COMPONENT ***comp = eg_crx->comp;
	
	
	        /* install z-crossings */
	                                                                           
	for (i = 0; i <= gmax[0]; ++i)
        {
            coords1[0] = coords2[0] = L[0] + i*h[0];
            for (j = 0; j <= gmax[1]; ++j)
            {
                coords1[1] = coords2[1] = L[1] + j*h[1];
                for (k = 0; k <= (int)(0.5*gmax[2]); ++k)
                {
                    if (((comp[i][j][k] == comp0) && (comp[i][j][k+1] == comp1)) ||
                        ((comp[i][j][k] == comp1) && (comp[i][j][k+1] == comp0)))
                    {
                        coords1[2] = L[2] + k*h[2];
                        coords2[2] = L[2] + (k+1)*h[2];

                        f1 = (*func)(func_params,coords1);
                        f2 = (*func)(func_params,coords2);
                        if((f1 < 0 && f2 > 0) || (f1 > 0 && f2 < 0))
                        {
                            for(m =0; m <3; m++)
                                crds_crx[m] = (coords1[m] + coords2[m])/2.0;
                        }
                        else
                        {
                            screen("ERROR: in install_grid_bottombdry_crx(), no z-crxing!");
                            (void)printf("i = %d,j =%d,k =%d\n",i,j,k);
                            clean_up(ERROR);
                        }

		        crx_store[n_crx].p = Point(crds_crx);
		        z_crx[i][j][k] = &crx_store[n_crx++];
                    }
		}
	    }
	}
	return n_crx;
}    /*end install_grid_bottombdry_crx*/

/**********************************************************************
 *     This function sets NO_COMP outside the compution domain 
 *     in z direction in order to make boundary surface for jet3d.
 *********************************************************************/

LOCAL   void assign_exterior_comp(
	COMPONENT ***comp,
	RECT_GRID gr)
{
        int i,j,k;
        int *gmax = gr.gmax;
        double *GL = gr.GL;
        double *GU = gr.GU;
        double *L = gr.L;
        double *h = gr.h;
        double coords[3];

        for (i = 0; i <= gmax[0]; ++i)
        {
            coords[0] = L[0] + i*h[0];
            for (j = 0; j <= gmax[1]; ++j)
            {
                coords[1] = L[1] + j*h[1];
                for (k = 0; k <= gmax[2]; ++k)
                {
                    coords[2] = L[2] + k*h[2];
                    if(coords[2] <= GU[2] && coords[2] >= GL[2] )
                        continue;
                    comp[i][j][k] = NO_COMP;
	        }
            }
        }

}

/* Finding curve_crx acording to the following rule:
*    (a) There are exactly 3 component on the four vertex.
*    (b) components at diagonal position are distinct. c00 != c11, c01 != c11
*
*/

LOCAL   int install_curve_crx(
        double (*func_1)(POINTER,double*),
	POINTER func_1_params,
	double (*func_2)(POINTER,double*),
	POINTER func_2_params,
	EG_CRX *eg_crx,
	RECT_GRID grid,
	int *n_crx,
	POINT  ***p,
	COMPONENT comp0,
	COMPONENT comp1,
	COMPONENT comp2,
	int    **curve_crx,
	int    curve_n[])
{
	double      ***face_coords, face_crds_crx[3];
	double      *L = grid.L;
	double      *h = grid.h;
	int        *gmax = grid.gmax;
	int        i,j,k;
	static int n_curve_crx;
	int        count_p0, count_p1, count_p2;
        	
	BBI_POINT  ****x_curve_crx = eg_crx->x_curve_crx;
	BBI_POINT  ****y_curve_crx = eg_crx->y_curve_crx;
	BBI_POINT  ****z_curve_crx = eg_crx->z_curve_crx;
	BBI_POINT  *crx_store = eg_crx->crx_store;
	COMPONENT  ***comp = eg_crx->comp;
	/*TMP */
        double      *cen,r2,h1,h2,h3,ratio,crds[2];
        JET_PARAMS *params;
        int        num_crx_curve1;
        
        params = (JET_PARAMS *)func_1_params;
        cen = params->cen;
        r2 = params->r2;
        ratio = params->ratio;
        h1 = params->cen[2] - params->h[0];
        h2 = h1 - params->h[1];
        h3 = h2 - params->h[2]*(1.0-ratio);

	/*End */
	count_p0 = count_p1 = count_p2 = 0;
	tri_array(&face_coords,3,2,2,sizeof(double));

	/* install x-crossings */
        for (k = 0; k < gmax[2]; ++k)
        {
            face_coords[2][0][0] = L[2] + k*h[2];
            face_coords[2][0][1] = L[2] + (k+1)*h[2];
            face_coords[2][1][0] = L[2] + k*h[2];
            face_coords[2][1][1] = L[2] + (k+1)*h[2];
				    
            for (j = 0; j < gmax[1]; ++j)
            {
	        face_coords[1][0][0] = L[1] + j*h[1];
	        face_coords[1][0][1] = L[1] + j*h[1];
	        face_coords[1][1][0] = L[1] + (j+1)*h[1];
	        face_coords[1][1][1] = L[1] + (j+1)*h[1];
		for (i = 0; i <= gmax[0]; ++i)
		{
		    x_curve_crx[i][j][k] = NULL;
		    if (is_curve_crx(comp[i][j][k],comp[i][j+1][k],
		            comp[i][j][k+1],comp[i][j+1][k+1]))
		    {
		        face_coords[0][0][0] = L[0] + i*h[0];
		        face_coords[0][0][1] = L[0] + i*h[0];
		        face_coords[0][1][0] = L[0] + i*h[0];
		        face_coords[0][1][1] = L[0] + i*h[0];
			if (! face_crx_in_dir(func_1,func_1_params,
				func_2,func_2_params,face_coords,face_crds_crx,0))
			{   
			    screen("ERROR: in install_curve_crx(), no x-curve_crxing!");
			    clean_up(ERROR);
			}
			
                        /* TMP */
			/*
			if(face_crds_crx[2] < 0.2 && face_crds_crx[2] > 0.19)
			{	
			    print_general_vector("x-face_crds_crx",face_crds_crx,3,"\n");
			    printf("comp: [%d] [%d] \n", comp[i][j][k+1], comp[i][j+1][k+1]); 
			    printf("comp: [%d] [%d] \n", comp[i][j][k], comp[i][j+1][k]); 
			    printf("which_3comp = %d\n",
				    which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                    comp[i][j+1][k],comp[i][j][k+1],
                                    comp[i][j+1][k+1]));
                        }
			*/

			if (face_crds_crx[1] - face_coords[1][0][0] < 0.004*h[1])
			    face_crds_crx[1] = face_coords[1][0][0] + 0.004*h[1];
			if (face_coords[1][1][1] - face_crds_crx[1] < 0.004*h[1])
			    face_crds_crx[1] = face_coords[1][1][1] - 0.004*h[1];
			if (face_crds_crx[2] - face_coords[2][0][0] < 0.004*h[2])
			    face_crds_crx[2] = face_coords[2][0][0] + 0.004*h[2];
			if (face_coords[2][1][1] - face_crds_crx[2] < 0.004*h[2])
			    face_crds_crx[2] = face_coords[2][1][1] - 0.004*h[2];
			    
			crx_store[*n_crx].p = Point(face_crds_crx);
			if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                      comp[i][j+1][k],comp[i][j][k+1],
                                      comp[i][j+1][k+1]) == 1)
			{
			    curve_crx[0][curve_n[0]++]=*n_crx;
			    if(count_p0 < 2)
			    {
				p[0][count_p0] = crx_store[*n_crx].p;
				++count_p0;
				
			    }
			}
			if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                      comp[i][j+1][k],comp[i][j][k+1],
                                      comp[i][j+1][k+1]) == 2)
                        {
                            curve_crx[1][curve_n[1]++]=*n_crx;
                            if(count_p1 < 2)
                            {
                                 p[1][count_p1] = crx_store[*n_crx].p;
				 ++count_p1;
                            }
                        }
			if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                      comp[i][j+1][k],comp[i][j][k+1],
                                      comp[i][j+1][k+1]) == 3)
                        {
                             curve_crx[2][curve_n[2]++]=*n_crx;
			     if(count_p2 < 2) 
			     {
				 p[2][count_p2] = crx_store[*n_crx].p;
				 ++count_p2;
			     }
			} 
     		        x_curve_crx[i][j][k] = &crx_store[(*n_crx)++];
			n_curve_crx++;
		    }
		}
	    }
	}

	/* install y-crossings */
        for (k = 0; k < gmax[2]; ++k)
        {
            face_coords[2][0][0] = L[2] + k*h[2];
            face_coords[2][0][1] = L[2] + k*h[2];
            face_coords[2][1][0] = L[2] + (k+1)*h[2];
            face_coords[2][1][1] = L[2] + (k+1)*h[2];
						    
 	    for (i = 0; i < gmax[0]; ++i)
            {
	        face_coords[0][0][0] = L[0] + i*h[0];
	        face_coords[0][0][1] = L[0] + (i+1)*h[0];
	        face_coords[0][1][0] = L[0] + i*h[0];
	        face_coords[0][1][1] = L[0] + (i+1)*h[0];
		for (j = 0; j <= gmax[1]; ++j)
		{
		    y_curve_crx[i][j][k] = NULL;
		    if (is_curve_crx(comp[i][j][k],comp[i][j][k+1],
		            comp[i+1][j][k],comp[i+1][j][k+1]))
		    {
		        face_coords[1][0][0] = L[1] + j*h[1];
		        face_coords[1][0][1] = L[1] + j*h[1];
		        face_coords[1][1][0] = L[1] + j*h[1];
		        face_coords[1][1][1] = L[1] + j*h[1];
			if (! face_crx_in_dir(func_1,func_1_params,
				func_2,func_2_params,face_coords,face_crds_crx,1))
			{
			    screen("ERROR: in install_curve_crx(), no y-curve_crxing!");
			    clean_up(ERROR);
			}
			
			if (face_crds_crx[0] - face_coords[0][0][0] < 0.004*h[0])
			    face_crds_crx[0] = face_coords[0][0][0] + 0.004*h[0];
			if (face_coords[0][1][1] - face_crds_crx[0] < 0.004*h[0])
			    face_crds_crx[0] = face_coords[0][1][1] - 0.004*h[0];
			if (face_crds_crx[2] - face_coords[2][0][0] < 0.004*h[2])
			    face_crds_crx[2] = face_coords[2][0][0] + 0.004*h[2];
			if (face_coords[2][1][1] - face_crds_crx[2] < 0.004*h[2])
			    face_crds_crx[2] = face_coords[2][1][1] - 0.004*h[2];
			    
			crx_store[*n_crx].p = Point(face_crds_crx);

			if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
			             comp[i][j][k+1],comp[i+1][j][k],
		                     comp[i+1][j][k+1]) == 1)
		        {
		            curve_crx[0][curve_n[0]++]=*n_crx;
		            if(count_p0 < 2)
		            {
		                p[0][count_p0] = crx_store[*n_crx].p;
				++count_p0;
		            }
		        }
                        if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                      comp[i][j][k+1],comp[i+1][j][k],
                                      comp[i+1][j][k+1]) == 2)
                        {
                            curve_crx[1][curve_n[1]++]=*n_crx;
                            if(count_p1 < 2)
                            {
                                p[1][count_p1] = crx_store[*n_crx].p;
				++count_p1;
                            }
                        }
			if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                      comp[i][j][k+1],comp[i+1][j][k],
                                      comp[i+1][j][k+1]) == 3)
                        {
                            curve_crx[2][curve_n[2]++]=*n_crx;
                            if(count_p2 < 2)
                            {
                                p[2][count_p2] = crx_store[*n_crx].p;
				++count_p2;
                            }
                        }
			
    	 	        y_curve_crx[i][j][k] = &crx_store[(*n_crx)++];
			n_curve_crx++;
		    }
		}
	    }
	}

	/* install z-crossings */

	for (i = 0; i < gmax[0]; ++i)
	{
	    face_coords[0][0][0] = L[0] + i*h[0];
	    face_coords[0][0][1] = L[0] + i*h[0];
	    face_coords[0][1][0] = L[0] + (i+1)*h[0];
	    face_coords[0][1][1] = L[0] + (i+1)*h[0];
	    for (j = 0; j < gmax[1]; ++j)
	    {
		face_coords[1][0][0] = L[1] + j*h[1];
		face_coords[1][0][1] = L[1] + (j+1)*h[1];
		face_coords[1][1][0] = L[1] + j*h[1];
		face_coords[1][1][1] = L[1] + (j+1)*h[1];
		for (k = 0; k <= gmax[2]; ++k)
		{
		    z_curve_crx[i][j][k] = NULL;
		    if (is_curve_crx(comp[i][j][k],comp[i+1][j][k],
		            comp[i][j+1][k],comp[i+1][j+1][k]))
		    {   
		        face_coords[2][0][0] = L[2] + k*h[2];
		        face_coords[2][0][1] = L[2] + k*h[2];
		        face_coords[2][1][0] = L[2] + k*h[2];
		        face_coords[2][1][1] = L[2] + k*h[2];
			if (! face_crx_in_dir(func_1,func_1_params,
				func_2,func_2_params,face_coords,face_crds_crx,2))
			{
			    screen("ERROR: in install_curve_crx(), no z-curve_crxing!");
			    clean_up(ERROR);
			}
			
			if (face_crds_crx[0] - face_coords[0][0][0] < 0.004*h[0])
			    face_crds_crx[0] = face_coords[0][0][0] + 0.004*h[0];
			if (face_coords[0][1][1] - face_crds_crx[0] < 0.004*h[0])
			    face_crds_crx[0] = face_coords[0][1][1] - 0.004*h[0];
			if (face_crds_crx[1] - face_coords[1][0][0] < 0.004*h[1])
			    face_crds_crx[1] = face_coords[1][0][0] + 0.004*h[1];
			if (face_coords[1][1][1] - face_crds_crx[1] < 0.004*h[1])
			    face_crds_crx[1] = face_coords[1][1][1] - 0.004*h[1];
			
			crx_store[*n_crx].p = Point(face_crds_crx);
                        
			if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                     comp[i+1][j][k],comp[i][j+1][k],
                                     comp[i+1][j+1][k]) == 1)
                        {
                            curve_crx[0][curve_n[0]++]=*n_crx;
                            if(count_p0 < 2)
                            {
                                p[0][count_p0] = crx_store[*n_crx].p;
				++count_p0;
                            }
                        }
                        if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                      comp[i+1][j][k],comp[i][j+1][k],
                                      comp[i+1][j+1][k]) == 2)
                        {
                            curve_crx[1][curve_n[1]++]=*n_crx;
                            if(count_p1 < 2)
                            {
                                p[1][count_p1] = crx_store[*n_crx].p;
				++count_p1;
                            }
                        }
                        if(which_3comp(comp0,comp1,comp2,comp[i][j][k],
                                      comp[i+1][j][k],comp[i][j+1][k],
                                      comp[i+1][j+1][k]) == 3)
                        {
                            curve_crx[2][curve_n[2]++]=*n_crx;
                            if(count_p2 < 2)
                            {
                                p[2][count_p2] = crx_store[*n_crx].p;
				++count_p2;
                            }
                        }
			z_curve_crx[i][j][k] = &crx_store[(*n_crx)++];
			n_curve_crx++;
		    }
		}
	    }
	}
	free_these(1,face_coords);
	return n_curve_crx;
}	/* end install_curve_crx */


LOCAL	boolean exist_func_root(
	double x1,
	double x2)
{
	if (((x1 > 0) && (x2 > 0)) || ((x1 < 0) && (x2 < 0)))
	    return NO;
	return YES;
}

/* Modified Newton's Iteration: to find the root of two equations system:
*  f(x,y)=0  g(x,y)=0 in a given domian.. x0 < x < x1, y0 < y < y1..
***case-01
*                       /|\
*                        |y-axis
*                        |  
*                        |    b2
*                     c10 -----\--------------- c11
*		         |      \g(x,y)        |
*			 |       \           __|a2
*			 |        \       __/  |
*			 |         \   __/     |
*			 |         _o_/        |
*		         |f(x,y) _/  \         |
*			 |     _/     \        |
*			 |   _/        \       |
*		      c00 --/-----------\------ -------> x-axis
*			   a1           b1     c01 
*
***case-02
*                       /|\
*                        |y-axis
*                        |  
*                        |    
*                     c10 -\--------------- c11
*		         |  \              |
*			 |   \f(x,y)       |
*			 |    \            | / g(x,y)
*			 |     \           |/
*			 |      \          /
*		         |       \        /|
*			 |        \      / |
*			 |         \    /  |
*		      c00 ----------\--/--- -------> x-axis
*			           a1  b1     c01 
*
*
*            |df/dx  df/dy|         | x |          |f(x,y)|
* Jacobian = |            |     X = |   |   F(X) = |      |
*            |dg/dx  dg/dy|         | y |          |g(x,y)|
*
*           n+1   n     -1      0    n
*          X   = X  - Jacobian(X )F(X )
*/

LOCAL	int face_crx_in_dir(
	double (*func_1)(POINTER,double*),
	POINTER func_1_params,
	double (*func_2)(POINTER,double*),
	POINTER func_2_params,
	double ***face_coords,
	double *face_crds_crx,
	int dir)
{
	int i,j,k,num_iter;
	double c00[3],c01[3],c10[3],c11[3],a1[3],b1[3];
	double Jacobian[2][2],DetJ,InvJ[2][2];
	double f0,g0,func_f[2][2],func_g[2][2];

	double epsilon_1 = 1.0e-10;
	double coeff = 0.2;
	num_iter = 100;
	
	for (i = 0; i < 3; i++)
	{
	    c00[i] = face_coords[i][0][0];
	    c01[i] = face_coords[i][0][1];
	    c10[i] = face_coords[i][1][0];
	    c11[i] = face_coords[i][1][1];
	}
	func_f[0][0] = (*func_1)(func_1_params,c00);
	func_f[0][1] = (*func_1)(func_1_params,c01);
	func_f[1][0] = (*func_1)(func_1_params,c10);
	func_f[1][1] = (*func_1)(func_1_params,c11);
	func_g[0][0] = (*func_2)(func_2_params,c00);
	func_g[0][1] = (*func_2)(func_2_params,c01);
	func_g[1][0] = (*func_2)(func_2_params,c10);
	func_g[1][1] = (*func_2)(func_2_params,c11);
	
	for (i = 0; i < 2; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
		if ((fabs(func_f[i][j]) < epsilon) && (fabs(func_g[i][j]) < epsilon))
	        {
		    for (k = 0; k < 3; k++)
		        face_crds_crx[k] = face_coords[k][i][j];
	            return YES;
		}
	    }
	}
	
	/* set the central point of the rectangular grid as the initial trial point */
	for (i = 0; i < 3; i++)
	    face_crds_crx[i] = (face_coords[i][0][0] + face_coords[i][1][1])*0.5;

	/*TMP */
	{
	    double funcf,funcg,*cen,r1,r2,r3,h1,h2,h3,ratio;
	    JET_PARAMS *params,*params2;
	    
	    params = (JET_PARAMS *)func_1_params;
	    params2= (JET_PARAMS *)func_2_params;
	    cen = params->cen;
	    r1 = params->r1;
	    r2 = params->r2;
	    r3 = params2->r1;
	    ratio = params->ratio;
	    h1 = params->cen[2] - params->h[0];
	    h2 = h1 - params->h[1];
	    h3 = h2 - params->h[2]*(1.0-ratio);
	    if (ratio <= 0)
	        h3 = h2 - params->h[2];
	    
	    funcf = (*func_1)(func_1_params,face_crds_crx);
	    funcg = (*func_2)(func_2_params,face_crds_crx);
	    if (fabs(face_crds_crx[2] - h3) < 0.008)
	    {
		face_crds_crx[2] =  h3;    
	        if ((fabs(funcf) < epsilon))
		    return YES;
	        if(dir == 0)
	        {   
		    if(face_crds_crx[1] < 0)	
		        face_crds_crx[1] = -sqrt(sqr(r2)-sqr(face_crds_crx[0]));
		    else
			face_crds_crx[1] = sqrt(sqr(r2)-sqr(face_crds_crx[0]));    
	        }
	        if(dir == 1)
	        {
	            if(face_crds_crx[0] < 0)		
		        face_crds_crx[0] = -sqrt(sqr(r2)-sqr(face_crds_crx[1]));
		    else
			face_crds_crx[0] = sqrt(sqr(r2)-sqr(face_crds_crx[1]));    
	        }
		return YES;
	    }
	    if (fabs(face_crds_crx[2] - cen[2]) < epsilon && 
	         fabs(sqr(face_crds_crx[0])+sqr(face_crds_crx[1])-sqr(r1)) < 0.01 )
	    {
          	if ((fabs(funcf) < epsilon))
		    return YES;
		if(dir == 0)
		{
		    if(face_crds_crx[1] < 0)
		        face_crds_crx[1] = -sqrt(sqr(r1)-sqr(face_crds_crx[0]));
		    else
		        face_crds_crx[1] = sqrt(sqr(r1)-sqr(face_crds_crx[0]));
		}
		if(dir == 1)
		{
		    if(face_crds_crx[0] < 0)
		        face_crds_crx[0] = -sqrt(sqr(r1)-sqr(face_crds_crx[1]));
		    else
		        face_crds_crx[0] = sqrt(sqr(r1)-sqr(face_crds_crx[1]));
	        }
		face_crds_crx[2] = cen[2];
	        return YES;
	    }
	    if (fabs(face_crds_crx[2] - cen[2]) < epsilon &&
	        fabs(sqr(face_crds_crx[0])+sqr(face_crds_crx[1])-sqr(r3)) < 0.01 )
	    {
	        if ((fabs(funcf) < epsilon))
	            return YES;
	        if(dir == 0)
	        {
	            if(face_crds_crx[1] < 0)
	                face_crds_crx[1] = -sqrt(sqr(r3)-sqr(face_crds_crx[0]));
	            else
	                face_crds_crx[1] = sqrt(sqr(r3)-sqr(face_crds_crx[0]));
	        }
	        if(dir == 1)
	        {
	            if(face_crds_crx[0] < 0)
	                face_crds_crx[0] = -sqrt(sqr(r3)-sqr(face_crds_crx[1]));
	            else
	                face_crds_crx[0] = sqrt(sqr(r3)-sqr(face_crds_crx[1]));
	        }
		face_crds_crx[2]=cen[2];
	        return YES;
	   }
	}
	/*END */
	/* calculate Jacobian */
        Jacobian[0][0] = partial_func_derivative(func_1,
                         func_1_params,face_crds_crx,(dir+1)%3);
        Jacobian[0][1] = partial_func_derivative(func_1,
                         func_1_params,face_crds_crx,(dir+2)%3);
        Jacobian[1][0] = partial_func_derivative(func_2,
                         func_2_params,face_crds_crx,(dir+1)%3);
        Jacobian[1][1] = partial_func_derivative(func_2,
                         func_2_params,face_crds_crx,(dir+2)%3);


	/* calculate determent of Jacobian */
	DetJ = Jacobian[0][0]*Jacobian[1][1] - Jacobian[0][1]*Jacobian[1][0];
        if(DetJ != 0)
        {	
	    /* calcunate inverse of Jacobian */
	    InvJ[0][0] =  Jacobian[1][1]/DetJ;
	    InvJ[0][1] = -Jacobian[0][1]/DetJ;
	    InvJ[1][0] = -Jacobian[1][0]/DetJ;
	    InvJ[1][1] =  Jacobian[0][0]/DetJ;

	    
	    f0 = (*func_1)(func_1_params,face_crds_crx);
	    g0 = (*func_2)(func_2_params,face_crds_crx);
	
	    /* calculate root of g(x,y)=0 and f(x,y)=0 with Modified Newton Method */
	    for (i = 0; i < num_iter; ++i)
	    {
	        face_crds_crx[(dir+1)%3] += -InvJ[0][0]*f0 - InvJ[0][1]*g0;
	        face_crds_crx[(dir+2)%3] += -InvJ[1][0]*f0 - InvJ[1][1]*g0;
	        f0 = (*func_1)(func_1_params,face_crds_crx);
	        g0 = (*func_2)(func_2_params,face_crds_crx);
	    }
        
                                                
	    /* assign curve_crx on the grid cell */
	    if ((fabs(f0) < epsilon_1) && (fabs(g0) < epsilon_1))
	    {
	        return YES;
	    }
	    else
	    {
	        if (exist_func_root(func_f[0][0],func_f[0][1]) && 
	            exist_func_root(func_g[0][0],func_g[0][1]))
	        {
	            if (! find_func_root(func_1,func_1_params,a1,c00,c01))
		    {
		        screen("ERROR: in face_crx_in_dir, No a1 root is "
		                   "found on 00 - 01 for f(x,y)\n");
		        clean_up(ERROR);
		    }
		    if (! find_func_root(func_2,func_2_params,b1,c00,c01))
		    {
		        screen("ERROR: in face_crx_in_dir, No b1 root is "
		                   "found on 00 - 01 for g(x,y)\n");
		        clean_up(ERROR);
		    }
		    for (i = 0; i < 3; i++)
		        face_crds_crx[i] = 0.5*(b1[i]+a1[i]);
		    face_crds_crx[(dir+1)%3] += coeff*fabs(b1[(dir+2)%3]-a1[(dir+2)%3]);
		    return YES;
	        }
	        else if (exist_func_root(func_f[0][0],func_f[1][0]) && 
	                 exist_func_root(func_g[0][0],func_g[1][0]))
	        {
	            if (! find_func_root(func_1,func_1_params,a1,c00,c10))
		    {
		        screen("ERROR: in face_crx_in_dir, No a1 root is "
		                   "found on 00 - 10 for f(x,y)\n");
		        clean_up(ERROR);
		    }
		    if (! find_func_root(func_2,func_2_params,b1,c00,c10))
		    {
		        screen("ERROR: in face_crx_in_dir, No b1 root is "
		                   "found on 00 - 10 for g(x,y)\n");
		        clean_up(ERROR);
		    }
		    for (i = 0; i < 3; i++)
		        face_crds_crx[i] = 0.5*(b1[i]+a1[i]);
		    face_crds_crx[(dir+2)%3] += coeff*fabs(b1[(dir+1)%3]-a1[(dir+1)%3]);
		    return YES;
	        }
	        else if (exist_func_root(func_f[1][0],func_f[1][1]) && 
	                 exist_func_root(func_g[1][0],func_g[1][1]))
	        {
	            if (! find_func_root(func_1,func_1_params,a1,c10,c11))
		    {
		        screen("ERROR: in face_crx_in_dir, No a1 root is "
		                   "found on 10 - 11 for f(x,y)\n");
		        clean_up(ERROR);
		    }
		    if (! find_func_root(func_2,func_2_params,b1,c10,c11))
		    {
		        screen("ERROR: in face_crx_in_dir, No b1 root is "
		                   "found on 10 - 11 for g(x,y)\n");
		        clean_up(ERROR);
		    }
		    for (i = 0; i < 3; i++)
		        face_crds_crx[i] = 0.5*(b1[i]+a1[i]);
		    face_crds_crx[(dir+1)%3] -= coeff*fabs(b1[(dir+2)%3]-a1[(dir+2)%3]);
		    return YES;
	        }
	        else if (exist_func_root(func_f[0][1],func_f[1][1]) && 
	                 exist_func_root(func_g[0][1],func_g[1][1]))
	        {
	            if (! find_func_root(func_1,func_1_params,a1,c01,c11))
		    {
		        screen("ERROR: in face_crx_in_dir, No a1 root is "
		                   "found on 01 - 11 for f(x,y)\n");
		        clean_up(ERROR);
		    }
		    if (! find_func_root(func_2,func_2_params,b1,c01,c11))
		    {
		        screen("ERROR: in face_crx_in_dir, No b1 root is "
		                   "found on 01 - 11 for g(x,y)\n");
		        clean_up(ERROR);
		    }
		    for (i = 0; i < 3; i++)
		        face_crds_crx[i] = 0.5*(b1[i]+a1[i]);
		    face_crds_crx[(dir+2)%3] -= coeff*fabs(b1[(dir+1)%3]-a1[(dir+1)%3]);
		    return YES;
	        }
                else
	        {
	           /*   screen("WARNING: in face_crx_in_dir, Cannt find a curve crx\n"); */
		    for (i = 0; i < 3; i++)
		        face_crds_crx[i] = 0.5*(face_coords[i][1][1]+face_coords[i][0][0]);
		    return YES;
	        }
            } 
	}
        return YES;
}	/* end face_crx_in_dir */

LOCAL	double partial_func_derivative(
	double (*func)(POINTER,double*),
	POINTER func_params,
	double *coords,
	int dir)
{
	if (func == ellipsoid_func)
	{
	    ELLIP_PARAMS *params;
	    const double *cen, *rad;
	    params = (ELLIP_PARAMS *)func_params;
	    cen = params->cen;
	    rad = params->rad;
	    return 2.0*(coords[dir]-cen[dir])/sqr(rad[dir]);
	}
	else if (func == bdry_box_func)
	{
	    return 0;
	}
	else if (func == plane_func)
	{
	    PLANE_PARAMS *params;
	    const double *N;
	    params = (PLANE_PARAMS*)func_params;
	    N = params->N;
	    return N[dir];
	}
	else if (func == jet_nozzle_func)
	{
	    JET_PARAMS *params;
	    double *cen,h1,h2,ratio;
	    params = (JET_PARAMS *)func_params;
	    cen = params->cen;
	    ratio = params->ratio;
	    h1 = params->h[0] + params->h[1];
	    h2 = params->h[0] + params->h[1]*ratio;
	    if (dir == 2) 
	        return 0;
	    else if (coords[2] >= h1)
	        return 2.0*(coords[dir]-cen[dir]);
	    else if (coords[2] >= h2)
	        return 2.0*(coords[dir]-cen[dir]);
	    else 
	        return 0;
	}
	else if (func == jet_plane_func)
	{
	    JET_PARAMS *params;
	    double *cen,h1,h2,ratio;
	    params = (JET_PARAMS *)func_params;
	    cen = params->cen;
	    ratio = params->ratio;
	    h1 = params->h[0];
	    if (dir == 2) 
	        return 0;
	    else if (coords[2] >= h1)
	        return 2.0*(coords[dir]-cen[dir]);
	    else if ((coords[2]-cen[2]) == (cen[0]-coords[0]))
	    {
	        if (dir == 0) return 1.0;
		if (dir == 2) return -1.0;
		else return 0.0;
	    }
	    else 
	        return 0;
	}

        if (func == jet_inner_func)
        {
            JET_PARAMS *params;
            double *cen,r0,r1,r2,h1,h2,h3,ratio;
            params = (JET_PARAMS *)func_params;
            cen = params->cen;
            r1 = params->r1;
            r2 = params->r2;
            ratio = params->ratio;
            h1 = params->cen[2] - params->h[0];
            h2 = h1 - params->h[1];
            h3 = h2 - params->h[2]*(1.0-ratio);
            
            if (dir == 2)
                return 0;
            else if (coords[2] >= h3)
                return 2.0*(coords[dir]-cen[dir]);
            else
                return 0;
           
        }
        if (func == jet_outer_func)
        {
            JET_PARAMS *params;
            double *cen,r0,r1,r2,h1,h2,h3,ratio;
            params = (JET_PARAMS *)func_params;
            cen = params->cen;
            r1 = params->r1;
            r2 = params->r2;
            ratio = params->ratio;
            h1 = params->cen[2] - params->h[0];
            h2 = h1 - params->h[1] - params->h[2]*(1.0-ratio);
            h3 = h2 - params->h[2]*ratio;
     
            if (dir ==2)
                return 0;
            else if (coords[2] >= h3)
            {
                return 2.0*(coords[dir]-cen[dir]); 
            }
            else 
               return 0;
        }

	else
	{
	    screen("WARNING: in partial_func_derivative, No function "
	                   "is matched\n"); 
	    clean_up(ERROR);
	    return 0;
	}
}	/* end partial_func_derivative */

LOCAL   double partial_func_derivative2(
        double (*func)(POINTER,double*),
        POINTER func_params,
        double *coords,
        int dir,
        int nd )
{
        if (func == jet_inner_func)
        {
            JET_PARAMS *params;
            double *cen,r0,r1,r2,h1,h2,h3,ratio;
            params = (JET_PARAMS *)func_params;
            cen = params->cen;
            r1 = params->r1;
            r2 = params->r2;
            ratio = params->ratio;
            h1 = params->cen[2] - params->h[0];
            h2 = h1 - params->h[1];
            h3 = h2 - params->h[2]*(1.0-ratio);

            if (ratio <= 0)
                h3 = h2 - params->h[2];

            if (ratio >= 0 && coords[2] < h3)
                return 0;
            else
                return 2.0*(coords[nd]-cen[nd]);
        }
        else if (func == jet_outer_func)
        {
            JET_PARAMS *params;
            double *cen,r0,r1,r2,h1,h2,h3,ratio;
            params = (JET_PARAMS *)func_params;
            cen = params->cen;
            r1 = params->r1;
            r2 = params->r2;
            ratio = params->ratio;
            h1 = params->cen[2] - params->h[0];
            h2 = h1 - params->h[1] - params->h[2]*(1.0-ratio);
            h3 = h2 - params->h[2]*ratio;
            if (ratio < 0)
                h2 = h3 = h1 - params->h[1] - params->h[2];
            if (coords[2] >= h1)
            {
                if ((dir == 0 && nd == 0) ||(dir == 1 && nd == 1))
                    return 0;
                else
                    return 2.0*(coords[nd]-cen[nd]);
            }
            else if (coords[2] > h3)
            {
                if (ratio < 0)
                {
                    if ((dir == 0 && nd == 0) ||(dir == 1 && nd == 1))
                        return 0;
                    else
                        return 2.0*(coords[nd]-cen[nd]);
                }   
                else
                {
                     r0 = sqr(r1) - sqr(coords[2]-h1);
                     if (coords[2] > h2)
                     {
                         if ((dir == 0 && nd == 0) ||(dir == 1 && nd == 1))
                             return 0;
                         else
                             return 2.0*(coords[nd]-cen[nd]);
                     }
                     else if((sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])) <= r0)
                     {
                         if ((dir == 0 && nd == 0) ||(dir == 1 && nd == 1))
                             return 0;
                         else
                             return -2.0*(coords[nd]-cen[nd]);
                     }
                 }
            }
            return 0; 
        }

}

/* end partial_func_derivative2 */

LOCAL	double return_sign(double f)
{
	if (f >= 0) return 1.0;
	return -1.0;
}	/* end return_sign */

/* find function root with bisection method
*  at interval [a,b], f(a)*f(b)<0, 
*  Bisection method:
*  1. s = sign(f(a));
*  2. m = (a+b)/2;
*  3. if sign(f(m)) == s, then a = m;
*     else b = m;
*  4. return to 1. unless |b-a| is sufficiently small
*/

LOCAL	boolean find_func_root(
        double (*func)(POINTER,double*), 
	POINTER func_params,
	double *c0,
	double *c1,
	double *c2)
{               
	int i,ii, n_crx, num_iter = 0;
	
	double fa,fb,fm,s,ca[3],cb[3];

	fa = (*func)(func_params,c1);
	fb = (*func)(func_params,c2);
	if (!exist_func_root(fa,fb)) return NO;
	
	if (fabs(fa) <= epsilon) 
	{
	    for (i = 0; i < MAXD; i++)
                c0[i] = c1[i];
	    return YES;
	}
	if (fabs(fb) <= epsilon) 
	{
	    
	    for (i = 0; i < MAXD; i++)
                c0[i] = c2[i];
	    return YES;
	}

	n_crx = 0;

	if (c1[0] != c2[0])
	{
	    ii = 0;
	    ++n_crx;
	}
	if (c1[1] != c2[1])
	{
	    ii = 1;
	    ++n_crx;
	}
	if (c1[2] != c2[2])
	{
	    ii = 2;
	    ++n_crx;
	}
	if (n_crx == 0)
	{
	    screen("ERROR: in find_func_root\n");
	    clean_up(ERROR);
	}
	if (n_crx > 1)
	{
	    screen("ERROR: in find_func_root\n");
	    clean_up(ERROR);
	}
	
	for (i = 0; i < 3; i++)
	    c0[i] = ca[i] = cb[i] = c1[i];
	ca[ii] = c1[ii];    
	cb[ii] = c2[ii];
	
	while ((fabs((cb[ii]-ca[ii])/ca[ii]) > epsilon) && (++num_iter<100))
	{
	    c0[ii] = (ca[ii]+cb[ii])*0.5;
	    fa = (*func)(func_params,ca);
	    s = return_sign(fa);
	    fm = (*func)(func_params,c0); 
	    if (return_sign(fm) == s) ca[ii] = c0[ii];
	    else cb[ii] = c0[ii];
	}
	return YES;
}       /*end find_func_root*/

LOCAL	void return_original_coords(
	BLK_CRX *blk_crx,
	double *h)
{
	int i;
	for (i = 0; i < 2; i++)
	{
	    if (blk_crx->crx[1][i][0] != NULL)
	    {
	        h[0] = Coords(blk_crx->crx[1][i][0]->p)[0];
	        break;
	    }
	    if (blk_crx->crx[2][0][i] != NULL)
	    {
	        h[0] = Coords(blk_crx->crx[2][0][i]->p)[0];
	        break;
	    }
	    if (blk_crx->curve_crx[0][0] != NULL)
	    {
	        h[0] = Coords(blk_crx->curve_crx[0][0]->p)[0];
	        break;
	    }
	    if (blk_crx->crx[1][i][1] != NULL)
	    {
	        h[0] = Coords(blk_crx->crx[1][i][1]->p)[0]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    if (blk_crx->crx[2][1][i] != NULL)
	    {
	        h[0] = Coords(blk_crx->crx[2][1][i]->p)[0]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    if (blk_crx->curve_crx[0][1] != NULL)
	    {
	        h[0] = Coords(blk_crx->curve_crx[0][1]->p)[0]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    h[0] = 1000;
	}
	for (i = 0; i < 2; i++)
	{
	    if (blk_crx->crx[2][i][0] != NULL)
	    {
	        h[1] = Coords(blk_crx->crx[2][i][0]->p)[1];
	        break;
	    }
	    if (blk_crx->crx[0][0][i] != NULL)
	    {
	        h[1] = Coords(blk_crx->crx[0][0][i]->p)[1];
	        break;
	    }
	    if (blk_crx->curve_crx[1][0] != NULL)
	    {
	        h[1] = Coords(blk_crx->curve_crx[1][0]->p)[1];
	        break;
	    }
	    if (blk_crx->crx[2][i][1] != NULL)
	    {
	        h[1] = Coords(blk_crx->crx[2][i][1]->p)[1]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    if (blk_crx->crx[0][1][i] != NULL)
	    {
	        h[1] = Coords(blk_crx->crx[0][1][i]->p)[1]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    if (blk_crx->curve_crx[1][1] != NULL)
	    {
	        h[1] = Coords(blk_crx->curve_crx[1][1]->p)[1]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    h[1] = 1000;
	}
	for (i = 0; i < 2; i++)
	{
	    if (blk_crx->crx[0][i][0] != NULL)
	    {
	        h[2] = Coords(blk_crx->crx[0][i][0]->p)[2];
	        break;
	    }
	    if (blk_crx->crx[1][0][i] != NULL)
	    {
	        h[2] = Coords(blk_crx->crx[1][0][i]->p)[2];
	        break;
	    }
	    if (blk_crx->curve_crx[2][0] != NULL)
	    {
	        h[2] = Coords(blk_crx->curve_crx[2][0]->p)[2];
	        break;
	    }
	    if (blk_crx->crx[0][i][1] != NULL)
	    {
	        h[2] = Coords(blk_crx->crx[0][i][1]->p)[2]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    if (blk_crx->crx[1][1][i] != NULL)
	    {
	        h[2] = Coords(blk_crx->crx[1][1][i]->p)[2]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    if (blk_crx->curve_crx[2][1] != NULL)
	    {
	        h[2] = Coords(blk_crx->curve_crx[2][1]->p)[2]-RANGE_GRID/NUMBER_GRID;
	        break;
	    }
	    h[2] = 1000;
	}
}	/* end return_original_coords */

LOCAL	void	check_consistency_of_crx(
	BLK_CRX	*blk_crx,
	boolean	include_curve_crx,
	int	num_case)
{
	int dir,i,j,k;
	double h[3];

	return_original_coords(blk_crx,h);
	
	for (dir = 0; dir < 3; dir++)
	{
	    for (i = 0; i < 2; ++i)
	    {
	        for (j = 0; j < 2; ++j)
	        {
	            if (blk_crx->crx[dir][i][j] != NULL)
	            {
		        for (k = 0; k < 3; k++)
			{
	                    if ((Coords(blk_crx->crx[dir][i][j]->p)[k] < h[k]-epsilon) ||
		                (Coords(blk_crx->crx[dir][i][j]->p)[k] > 
				(h[k]+RANGE_GRID/NUMBER_GRID+epsilon)))
		            {
			        screen("ERROR in crx position at blk_crx->crx[%d][%d][%d] in " 
			        "%d direction for case-%d\n",dir,i,j,k,num_case);
				return;
			    }
			}
		    }
		}
		if ((include_curve_crx) && (blk_crx->curve_crx[dir][i] != NULL))
		{
		    for (k = 0; k < 3; k++)
		    {
		        if ((Coords(blk_crx->curve_crx[dir][i]->p)[k] < h[k]-epsilon) ||
		          (Coords(blk_crx->curve_crx[dir][i]->p)[k] > 
			  (h[k]+RANGE_GRID/NUMBER_GRID+epsilon)))
		        {
			    screen("ERROR in curve_crx position at blk_crx->curve_crx[%d][%d] in "
		            "%d direction  for case-%d\n",dir,i,k,num_case);
		            return;
			}
		    }
	        }
	    }
	}
}	/* end check_consistency_of_crx */


/*              
 *    for Z and X cut   Y (Z)      
 *                      ^ 
 *    s3   s4           |
 *    s1   s2           |------> X (Y)
 *                     /
 *                    / 
 *                   v
 *                   Z (X)
 *    for Y cut
 *                     Z
 *    s2  s4           ^   ^ Y
 *    s1  s3           |  /
 *                     | /
 *                     |/------> X 
 *                   
 */
LOCAL   int  which_3comp(
	COMPONENT   comp0,
	COMPONENT   comp1,
	COMPONENT   comp2,
	COMPONENT   comps1,
	COMPONENT   comps2,
	COMPONENT   comps3,
	COMPONENT   comps4)
{       
	/* NEW. Exclude diagonal case */
	/*
	if(comps1 != comps2 &&
	   comps3 != comps4 &&
           comps1 != comps3 &&
           comps2 != comps4)	   
	    return 0;
	*/    
	/* END NEW */

	if(comps1==comps2)
	{
	    if(!(comps3==comps4))
            {		    
                if(comps1==comp2 && ((comps3==comp0 && comps4==comp1) 
			||(comps3==comp1 && comps4==comp0)))
		    return 1;
		/* NEW */
                if(comps1==comp0 && ((comps3==comp2 && comps4==comp1)
		        ||(comps3==comp1 && comps4==comp2)))
	            return 1;
		/* END NEW */
		if(comps1==NO_COMP && ((comps3==comp0 && comps4==comp1)
					||(comps3==comp1 && comps4==comp0)))
		    return 2;
		if(comps1==NO_COMP && ((comps3==comp1 && comps4==comp2)
					||(comps3==comp2 && comps4==comp1)))
	            return 3;
	    }
	}
	if(comps1==comps3)
        {
            if(!(comps2==comps4))
            {
                if(comps1==comp2 && ((comps2==comp0 && comps4==comp1) 
			||(comps2==comp1 && comps4==comp0))) 
		    return 1;
		/* NEW */
                if(comps1==comp0 && ((comps2==comp2 && comps4==comp1)
		        ||(comps2==comp1 && comps4==comp2)))
	            return 1;
		/* END NEW */
		if(comps1==NO_COMP && ((comps2==comp0 && comps4==comp1)
			||(comps2==comp1 && comps4==comp0)))
                    return 2;
                if(comps1==NO_COMP && ((comps2==comp1 && comps4==comp2)
			||(comps2==comp2 && comps4==comp1)))
	            return 3;
	     }
	}
	/*  comps2==comps3 and comps1==comps4 are diagonal case */
	if(comps2==comps3)
        {
            if(!(comps1==comps4))
            {
                if(comps2==comp2 && ((comps1==comp0 && comps4==comp1)
			||(comps1==comp1 && comps4==comp0)))
		    return 1;
                if(comps2==NO_COMP && ((comps1==comp0 && comps4==comp1)
			||(comps1==comp1 && comps4==comp0)))
                    return 2;
                if(comps2==NO_COMP && ((comps1==comp1 && comps4==comp2) 
			||(comps1==comp2 && comps4==comp1)))
                    return 3;
            }
        }
	if(comps1==comps4)
        {
            if(!(comps2==comps3))
            {
                if(comps1==comp2 && ((comps2==comp0 && comps3==comp1)
			||(comps2==comp1 && comps3==comp0)))
		    return 1;
                if(comps1==NO_COMP && ((comps2==comp0 && comps3==comp1)
			||(comps2==comp1 && comps3==comp0)))
                    return 2;
                if(comps1==NO_COMP && ((comps2==comp1 && comps3==comp2)
			||(comps2==comp2 && comps3==comp1)))
                    return 3;
             }
         }
	if(comps2==comps4)
        {
            if(!(comps1==comps3))
            {
                if(comps2==comp2 && ((comps1==comp0 && comps3==comp1)
		        ||(comps1==comp1 && comps3==comp0)))
	            return 1;
		/* NEW */
                if(comps2==comp0 && ((comps1==comp2 && comps3==comp1)
		        ||(comps1==comp1 && comps3==comp2)))
	            return 1;
		/* END NEW */
	        if(comps2==NO_COMP && ((comps1==comp0 && comps3==comp1) 
			||(comps1==comp1 && comps3==comp0)))
                    return 2;
                if(comps2==NO_COMP && ((comps1==comp1 && comps3==comp2)
        		||(comps1==comp2 && comps3==comp1)))
                    return 3;
            }
         }
         if(comps4==comps3)
         {
             if(!(comps1==comps2))
             {
                if(comps3==comp2 && ((comps1==comp0 && comps2==comp1)
                        ||(comps1==comp1 && comps2==comp0)))
                    return 1;
		/* NEW */
                if(comps3==comp0 && ((comps1==comp2 && comps2==comp1)
		        ||(comps1==comp1 && comps2==comp2)))
	            return 1;
		/* END NEW */
                if(comps3==NO_COMP && ((comps1==comp0 && comps2==comp1)
                        ||(comps1==comp1 && comps2==comp0)))
                    return 2;
                if(comps3==NO_COMP && ((comps1==comp1 && comps2==comp2)
                        ||(comps1==comp2 && comps2==comp1)))
                    return 3;
              }
          }
          return 0;		 
}	
		    
/* This function is different from reverse_bond, the tris on 
 * the bond do not adjust.
 */
LOCAL  void reverse_bond_only(
        BOND *b)
{
        INTERFACE *intfc = current_interface();
                                                                                                                                     
        if (intfc == NULL)
            return;
	(*i_user_interface(intfc)._reverse_bond)(b);
        return;
}       /*end reverse_bond_only */
                                                                                                                                     
#endif /* defined(THREED) */
