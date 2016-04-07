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
*				grectstat.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Rect state statistics are based on the rectangular grid and ignore
*	the presence of fronts.  This type of statistic is easy to implement
*	and is more appropriate than the other types of statistics when there
*	are singularities in the quantity being computed, such as divergence
*	of velocity across a shock.
*/

#if defined(TWOD) || defined(THREED)

#define DEBUG_STRING "rsstat"

#include <gdecs/gdecs.h>

LOCAL const char *HEADERFIELD = "%-22s";
LOCAL const char *FIELD       = "%-22g";
LOCAL const char *INT_FIELD   = "%-22d";

	/* Maximum number of rect state statistics */

enum {MAX_NUM_RECT_STATS = 20};

	/* Possible quantities for rect state based statistics */

enum _RECT_STATISTIC {
	NO_RECT_STATISTIC = -1,
	AVG_RHO_CURLV_SQR =  1,
	R_XX,
	R_YY,
	MIXED_ZONE_FRAC,
	POT_ENG,
	RHO_SQR,
	KE_TOTAL,
	E_TOTAL,
	KE_X,
	KE_Y,
	RHOX,
	RHOY,
	RHOVX,
	RHOVY
};
typedef enum _RECT_STATISTIC RECT_STATISTIC;


typedef struct {
	OUTPUT_DATA	odata;

		/* Indices for rect state stats controls */
	int		num_rect_stats;
	RECT_STATISTIC	rect_stats_selected[MAX_NUM_RECT_STATS];
	int		num_rect_stat_vars;
	int		refinement_ratio;  /* fineness for stats relative
					    * to computational grid*/
	double		*rect_height_data;

	int		RXX_selected;
	int		RYY_selected;

	FILE		*rss_height_file;
	FILE		*rss_time_file;

	char		*rss_height_dname,	*rss_height_filename;
	char		*rss_time_dname,	*rss_time_filename;

} Rect_state_stat_data;

#define RSS_data(data)		((Rect_state_stat_data *)data)

	/* LOCAL Function Declarations */
LOCAL	void	compute_rect_state_stats_at_layer(Rect_state_stat_data*,
						  double,int,Front*,Wave*);
LOCAL	void	print_rss_height_header(Rect_state_stat_data*,Front*);
LOCAL	void	print_rss_time_header(Rect_state_stat_data*);
LOCAL	void	printout_rect_state_stats(Rect_state_stat_data*,Front*);
LOCAL	void	record_rect_state_stats(Grid*,Wave*,Front*,Printplot*,
					OUTPUT_DATA*,boolean);
LOCAL	double	sum_stat_vector(double*,int,int,int);


/*
*			init_rect_state_stats():
*/

EXPORT	void init_rect_state_stats(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Printplot	*prt)
{
	Rect_state_stat_data *rssdata;
	int		i, entry;
	char		s[Gets_BUF_SIZE];
	
	screen("Type 'y' to request rect state statistics: ");
	(void) Gets(s);
	if ((s[0] != 'Y') && (s[0] != 'y')) return;

	scalar(&rssdata,sizeof(Rect_state_stat_data));

		/* Default settings for Rect_state_stat_data structure */
	rssdata->num_rect_stats = 0;
	rssdata->num_rect_stat_vars = 0;
	rssdata->refinement_ratio = 1;
	rssdata->RXX_selected = rssdata->RYY_selected = NO;
	rssdata->num_rect_stats = 0;
	for (i = 0; i < MAX_NUM_RECT_STATS; i++)
	    rssdata->rect_stats_selected[i] = NO_RECT_STATISTIC;
	rssdata->rss_height_file = NULL;
	rssdata->rss_time_file = NULL;

	screen("Enter selections, one to a line, (end with -1)\n\n"
	       "\t 1) Average squared mass weighted vorticity\n"
	       "\t 2) Average R_xx, KE_x, m_x\n"
	       "\t 3) Average R_yy, KE_y, m_y\n"
	       "\t 4) Fraction of mixed zones\n"
	       "\t 5) Potential Energy\n"
	       "\t 6) Density squared\n"
	       "\t 7) Total Kinetic Energy\n"
	       "\t 8) Total Energy\n"
	       "\t 9) X contribution to Kinetic Energy\n"
	       "\t 10) Y contribution to Kinetic Energy\n"
	       "\t 11) Total Mass (X)\n"
	       "\t 12) Total Mass (Y)\n"
	       "\t 13) X component Momentum\n"
	       "\t 14) Y component Momentum\n");

	for (i = 0; i <= MAX_NUM_RECT_STATS; i++)
	{
		screen("\t: ");
		(void) Gets(s);
		(void) sscanf(s,"%d",&entry);
		if (entry <= 0)
			break;

		if (i == MAX_NUM_RECT_STATS)
		{
	    	screen("ERROR in init_rect_state_stats(), "
	    	       "too many stats requested.\n");
			clean_up(ERROR);
		}

	    switch (entry)
	    {
	    case AVG_RHO_CURLV_SQR:
	        rssdata->rect_stats_selected[i] = AVG_RHO_CURLV_SQR;
		break;
	    case R_XX:
	        rssdata->rect_stats_selected[i] = R_XX;
		break;
	    case R_YY:
	        rssdata->rect_stats_selected[i] = R_YY;
		break;
	    case MIXED_ZONE_FRAC:
	        rssdata->rect_stats_selected[i] = MIXED_ZONE_FRAC;
		break;
	    case POT_ENG:
	        rssdata->rect_stats_selected[i] = POT_ENG;
		break;
	    case RHO_SQR:
	        rssdata->rect_stats_selected[i] = RHO_SQR;
		break;
	    case KE_TOTAL:
	        rssdata->rect_stats_selected[i] = KE_TOTAL;
		break;
	    case E_TOTAL:
	        rssdata->rect_stats_selected[i] = E_TOTAL;
		break;
	    case KE_X:
	        rssdata->rect_stats_selected[i] = KE_X;
		break;
	    case KE_Y:
	        rssdata->rect_stats_selected[i] = KE_Y;
		break;
	    case RHOX:
	        rssdata->rect_stats_selected[i] = RHOX;
		break;
	    case RHOY:
	        rssdata->rect_stats_selected[i] = RHOY;
		break;
	    case RHOVX:
	        rssdata->rect_stats_selected[i] = RHOVX;
		break;
	    case RHOVY:
	        rssdata->rect_stats_selected[i] = RHOVY;
		break;
	    }

		if (entry == R_XX)
		{
			rssdata->RXX_selected = YES;
			rssdata->rect_stats_selected[++i] = RHOX;
			rssdata->rect_stats_selected[++i] = KE_X;
			rssdata->rect_stats_selected[++i] = RHOVX;
		}
		else if (entry == R_YY)
		{
			rssdata->RYY_selected = YES;
			rssdata->rect_stats_selected[++i] = RHOY;
			rssdata->rect_stats_selected[++i] = KE_Y;
			rssdata->rect_stats_selected[++i] = RHOVY;
		}
	}

	rssdata->num_rect_stats = i;
	if (rssdata->num_rect_stats <= 0)
	{
	    screen("Error in init_rect_state_stats(), "
		   "stats requested, but no variables specified.\n");
		clean_up(ERROR);
	}

	screen("\nEnter refinement ratio (default = %d): ",
	       rssdata->refinement_ratio);
	(void) Gets(s);
	if (s[0] != '\0')
		(void) sscanf(s,"%d",&rssdata->refinement_ratio);

	/* Assumes upper global boundaries are not reflecting or periodic. */
	rssdata->num_rect_stat_vars =
	    (front->pp_grid->Global_grid.gmax[1]+1)*rssdata->refinement_ratio*
		(rssdata->num_rect_stats + 1 - rssdata->RXX_selected -
		 rssdata->RYY_selected);

	uni_array(&rssdata->rect_height_data, 
	       rssdata->num_rect_stat_vars,sizeof(double));

	init_output_data(init,&rssdata->odata,grid,prt,NULL,NO,NO,NO);
	add_user_output_function(record_rect_state_stats,&rssdata->odata,prt);

	rssdata->rss_height_file =
	    open_data_file(front,"rect state height data",YES,NO,NULL,
			   &rssdata->rss_height_dname,NULL,
			   &rssdata->rss_height_filename);

	rssdata->rss_time_file =
		open_data_file(front,"rect state time data",
			       YES,NO,NULL,&rssdata->rss_time_dname,
			              NULL,&rssdata->rss_time_filename);

	screen("\n");

} 		/*end init_rect_state_stats*/


/*ARGSUSED*/
LOCAL void record_rect_state_stats(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*odata,
	boolean		about_to_stop)
{
	Rect_state_stat_data	*rssdata = RSS_data(odata);
	RECT_GRID	*rgr = wave->rect_grid;
	boolean		first_layer = NO;
	double		z;
	int		iy, iz;
	int		fiy, fimin, fimax;
	int 		*gmax = rgr->gmax;
	int		refinement_ratio = rssdata->refinement_ratio;
	int		dim = rgr->dim;
	static boolean	first = YES;
	static double	dh = 0.0;
	static int	iymin, iymax;
	static int	izmin, izmax;
	static int	top_layer;	/* Index of top mz layer */
	static int	bottom_layer;	/* Index of bottom mz layer */

	DEBUG_ENTER(record_rect_state_stats)

	if (first == YES)
	{
	    first = NO;

	    iymin = 0;	iymax = gmax[1]+1;
	    izmin = 0;	izmax = 1;

	    if (dim == 2)
	    {
	    	top_layer    = iymax-1;
	    	bottom_layer = iymin;
	    }
	    else if (dim == 3)
	    {
	    	top_layer    = izmax-1;
	    	bottom_layer = izmin;
	    }
	    dh = rgr->h[dim-1]/refinement_ratio;

	    print_rss_time_header(rssdata);
	}

	print_rss_height_header(rssdata,front);

	for (iz = izmin; iz < izmax; iz++)
	{
	    first_layer = YES;

	    for (iy = iymin; iy < iymax; iy++)
	    {
		fiy = refinement_ratio/2;
	        if (iy == bottom_layer)
		    fimin = 0;
		else 
	        {
		    fimin = (refinement_ratio%2 == 1) ?
			    -fiy : (-fiy) + 1;
	        }
		fimax = (iy == top_layer) ? 0 : fiy;

		/* This loop breaks up each layer based on refinement_ratio */
		for (fiy = fimin; fiy <= fimax; fiy++)
		{
		    z = cell_edge(iy,dim-1,rgr) + fiy*dh;
		    if (z == rgr->GL[1])
			z += EPSILON; /*TOLERANCE*/
		    else if (z == rgr->GU[1])
			z -= EPSILON; /*TOLERANCE*/

		    compute_rect_state_stats_at_layer(rssdata,z,first_layer,
						      front,wave);
		    first_layer = NO;
	        }
	    }
        }

	printout_rect_state_stats(rssdata,front);

	print_graph_footer(rssdata->rss_height_file,"RECT STATE HEIGHT",YES);
	if (about_to_stop == YES)
		print_graph_footer(rssdata->rss_time_file,
				   "RECT STATE TIME",YES);

	DEBUG_LEAVE(record_rect_state_stats)

}		/*end record_rect_state_stats*/


LOCAL void compute_rect_state_stats_at_layer(
	Rect_state_stat_data *rssdata,
	double 		z,
	int		first_layer,
	Front 		*front,
	Wave 		*wave)
{
	COMPONENT	comp;
	INTERFACE	*intfc = front->interf;
	double		rhocurlvsqr = 0.0;
	double		ke_x = 0.0, ke_y = 0.0, ke_tot = 0.0, e_tot = 0.0;
	double		rhox = 0.0, rhoy = 0.0;
	double		rhovx = 0.0, rhovy = 0.0;
	double		rho_sqr = 0.0, pot_en = 0.0;
	double		contribution;
	double		coords[MAXD], x, dx, dy;
	double		grav;
	int 		i, j, k, ii, jj, *gmax = wave->rect_grid->gmax;
	static Locstate	N = NULL, S = NULL, E = NULL, W = NULL, here = NULL;
	static boolean	first = YES;
	static double	XL, XU, YU, YL, hx, hy, h, mix_frac;
	static double	*rect_height_data = NULL;
	static int	fine, num_rect_stats, stat_posn;
	static int	offset;

	if (first == YES)
	{
	    size_t 	sizest = front->sizest;
	    RECT_GRID	*rgr = wave->rect_grid;
	    PP_GRID	*pp_grid = front->pp_grid;
	    int		icoords[MAXD];
	    int		my_id = pp_mynode();

	    first = NO;

	    alloc_state(front->interf,&N, sizest);
	    alloc_state(front->interf,&W, sizest);
	    alloc_state(front->interf,&E, sizest);
	    alloc_state(front->interf,&S, sizest);

	    alloc_state(front->interf,&here, sizest); 

	    XL = rgr->L[0];
	    YL = rgr->L[1];
	    YU = rgr->U[1];
	    XU = rgr->U[0];

	    fine = rssdata->refinement_ratio;
	    h = rgr->h[0];
	    hx = 0.5*h/fine;
	    hy = 0.5*rgr->h[1]/fine;

	    num_rect_stats = rssdata->num_rect_stats;


	    find_Cartesian_coordinates(my_id,pp_grid,icoords);
	    offset = icoords[1]*fine*rgr->gmax[1]* (num_rect_stats + 1 -
			 rssdata->RXX_selected - rssdata->RYY_selected);
	    rect_height_data = rssdata->rect_height_data + offset;
	}

	if (first_layer)
	{
	    int		mixed = 0, icoords[MAXD];
     	    TRI_GRID	*grid = wave_tri_soln(wave)->tri_grid;
	    BLK_EL0		*blk_el0;

	    for (k = 0; k < num_rect_stats; k++)
	    {
	    	switch (rssdata->rect_stats_selected[k])
	    	{
	    	case MIXED_ZONE_FRAC:

		    /* Loop over tri_grid looking for
		     * mixed zones.  The first and last rows/columns
		     * would be counted as mixed if they have boundary
		     * curves, so we'll just skip them to make life
		     * easy.  This should have very little effect on
		     * the calculation.
		     */

	    	    for (jj = 1; jj < gmax[1]; jj++)
	    	    {
	    	        icoords[1] = jj;

	    	        for (ii = 1; ii < gmax[0]; ii++)
	    	        {
	    	    	    icoords[0] = ii;
	    	            blk_el0 = &Regular_blk_el0(icoords,grid);

			    if (!blk_el0_is_bilinear(blk_el0))
				mixed++;
			}
		    }
		    mix_frac = mixed/((gmax[0]-1.0)*(gmax[1]-1.0));
		    break;

		default:
		    break;
		}
	    }

	    zero_scalar(rssdata->rect_height_data, 
			rssdata->num_rect_stat_vars*FLOAT);
	    stat_posn = 0;
	}
	

	/* TODO:This section of code is very wasteful.  It calls hyp_solution()
	 * for every single state even though most of the time the state
	 * has already been obtained, either during this call to the function
	 * or during the previous one.  Make this more efficient!
	 */	
	for (i = 0; (i < gmax[0]) && (z < YU) && (z >= YL); i++)
	{
	    x = XL + i*h;
		
	    /* Keep us off the very edge of the domain */
	    if (i == 0)	x += EPSILON*h;	/* TOLERANCE */

	    for (j = 0; j < fine; j++)
	    {
	    	dx = 2.0*hx;	/* Distance between N & S, E & W 
	    			 * for derivatives*/
	    	dy = 2.0*hy;
	    	coords[0] = x + j*hx;
	    	coords[1] = z;

	    	comp = component(coords,intfc);
	        grav = gravity(coords,front->time)[1];
	    	hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
			     front,wave,here,NULL);

		    coords[1] += hy;

		    comp = component(coords,intfc);
		hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,front,
				wave,N,NULL);

		    if (is_obstacle_state(N))
		    {
			ft_assign(N,here,front->sizest);
			dy /= 2.0;
		    }

		    coords[1] -= 2.0*hy;

		    comp = component(coords,intfc);
		hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,front,
				wave,S,NULL);

		    if (is_obstacle_state(S))
		    {
			ft_assign(S,here,front->sizest);
			/*
			 * Note: by dividing again we're assuming
			 * that only one of S & N is an obstacle
			 */
			dy /= 2.0;
		    }

		    coords[0] += hx;
		    coords[1] += hy;

		    comp = component(coords,intfc);
		hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,front,
				wave,E,NULL);

		    if (is_obstacle_state(E))
		    {
			ft_assign(E,here,front->sizest);
			dx /= 2.0;
		    }

		    coords[0] -= 2.0*hx;

		    comp = component(coords,intfc);
		hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,front,
				wave,W,NULL);

		    if (is_obstacle_state(W))
		    {
			ft_assign(W,here,front->sizest);
			dx /= 2.0;
		    }
		
		for (k = 0; k < num_rect_stats; k++)
		{

		    switch (rssdata->rect_stats_selected[k])
		    {
		    case AVG_RHO_CURLV_SQR:
			contribution  = (Mom(E)[1] - Mom(W)[1])/dx;
			contribution -= (Mom(N)[0] - Mom(S)[0])/dy;
			rhocurlvsqr += sqr(contribution);
			break;

		    case R_XX:
		    case R_YY:
			break; /* These are functions of other averages */

		    case MIXED_ZONE_FRAC:
			break; /* Already calculated */

		    case POT_ENG:
			/*
			*  Potential energy increases against gravity,
			*  so we subtract
			*/	
			pot_en -= grav*Dens(here)*z;
			break;

		    case RHO_SQR:
			rho_sqr += sqr(Dens(here));
			break;

		    case KE_TOTAL:
		    	ke_tot += 0.5*Dens(here)*
				  (sqr(vel(0,here)) + sqr(vel(1,here)));
		    	break;

		    case E_TOTAL:
		    	e_tot += energy(here);
			break;

		    case KE_X:
		    	ke_x += Dens(here)*sqr(vel(0,here))/2.0;
		    	break;

		    case KE_Y:
		    	ke_y += Dens(here)*sqr(vel(1,here))/2.0;
		    	break;

		    case RHOX:
		    	rhox += Dens(here);
			break;

		    case RHOY:
		    	rhoy += Dens(here);
			break;

		    case RHOVX:
			rhovx += Mom(here)[0];
			break;

		    case RHOVY:
			rhovy += Mom(here)[1];
			break;

		    default:
			screen("ERROR: bad statistic index in "
			       "record_rect_state_stats()\n");
			clean_up(ERROR);
		    }
		}
	    }
	}

	rect_height_data[stat_posn++] = (z >= YL && z < YU) ? z : 0.0;
			
	for (k = 0; k < num_rect_stats; k++)
	{
	    switch (rssdata->rect_stats_selected[k])
	    {
	    case AVG_RHO_CURLV_SQR:
	    	rect_height_data[stat_posn++] = rhocurlvsqr;
	    	break;
	    case R_XX:
	    case R_YY:
	    	break;	/* These wait until we know <rho>, etc. */
	    case MIXED_ZONE_FRAC:
	    	rect_height_data[stat_posn++] = mix_frac;
	    	break;
	    case POT_ENG:
	    	rect_height_data[stat_posn++] = pot_en;
	    	break;
	    case RHO_SQR:
	    	rect_height_data[stat_posn++] = rho_sqr;
	    	break;
	    case KE_TOTAL:
	    	rect_height_data[stat_posn++] = ke_tot;
	    	break;
	    case E_TOTAL:
	    	rect_height_data[stat_posn++] = e_tot;
	    	break;	
	    case KE_X:
	    	rect_height_data[stat_posn++] = ke_x;
	    	break;
	    case KE_Y:
	    	rect_height_data[stat_posn++] = ke_y;
	    	break;
	    case RHOX:
	    	rect_height_data[stat_posn++] = rhox;
	    	break;
	    case RHOY:
	    	rect_height_data[stat_posn++] = rhoy;
	    	break;	
	    case RHOVX:
	    	rect_height_data[stat_posn++] = rhovx;
	    	break;		
	    case RHOVY:
	    	rect_height_data[stat_posn++] = rhovy;
	    	break;		
	    }
	}

}		/*end compute_rect_state_stats_at_layer*/


LOCAL	void	printout_rect_state_stats(
	Rect_state_stat_data	*rssdata,
	Front			*front)
{
	double		rho, rhovxsqr, rhovx, rhovysqr, rhovy, value;
	double		*rect_height_data;
	int		i, k, stat_posn, first_layer = YES;
	FILE		*height_file, *time_file;
	static 	boolean	first = YES;
	static	double	GXL, GXU, GYL, GYU, hx;
	static	int	stat_vector_stride;
	static	int	num_rect_stats, num_stats_up, num_nodes_across;
	static	long	num_stat_vars;

	if (first == YES)
	{
		first = NO;

		num_stat_vars    = rssdata->num_rect_stat_vars;
		num_rect_stats   = rssdata->num_rect_stats;
		num_nodes_across = front->pp_grid->gmax[0];
		/* Note: assumes top boundary not reflecting or periodic */
		num_stats_up = (front->pp_grid->Global_grid.gmax[1] + 1)*
				rssdata->refinement_ratio;
		GXL = front->rect_grid->GL[0];
		GXU = front->rect_grid->GU[0];
		GYL = front->rect_grid->GL[1];
		GYU = front->rect_grid->GU[1];
		hx  = front->rect_grid->h[0]/rssdata->refinement_ratio;
		/* The number of stored quantities in the uni_array
		 * rect_height_data[] is the number of statistics
		 * plus 1 for z minus 2 for R_xx and R_yy (these
		 * two are computed from other stats).
		 */
		stat_vector_stride = num_rect_stats + 1 - rssdata->RXX_selected
			                                - rssdata->RYY_selected;
	}

	rect_height_data = rssdata->rect_height_data;
	pp_global_sum(rect_height_data, num_stat_vars);

	height_file = rssdata->rss_height_file;
	time_file   = rssdata->rss_time_file;

	(void) fprintf(time_file,INT_FIELD,front->step);
	(void) fprintf(time_file,    FIELD,front->time);

	stat_posn = 0;

	for (i = 0; i < num_stats_up; i++)
	{
	    /* print z */
	    (void) fprintf(height_file,FIELD,
			   rect_height_data[stat_posn++]/num_nodes_across);

	    for (k = 0; k < num_rect_stats; k++)
	    {
		switch (rssdata->rect_stats_selected[k])
		{
		case AVG_RHO_CURLV_SQR:
		case POT_ENG:
		case RHO_SQR:
		case KE_TOTAL:
		case E_TOTAL:
		case KE_X:
		case KE_Y:
		case RHOX:
		case RHOY:
		case RHOVX:
		case RHOVY:
			if (first_layer)
			{
			    /* Compute total value of a variable */
			    value = sum_stat_vector(rect_height_data,
						    stat_vector_stride,
						    stat_posn,num_stats_up);
			    value *= hx*(GYU - GYL)/num_stats_up;
			    (void) fprintf(time_file, FIELD, value);
			}
			/* Compute row average of a variable */
			(void) fprintf(height_file, FIELD, 
				rect_height_data[stat_posn++]*hx/(GXU - GXL));
			break;

		case R_XX:
			rho = rect_height_data[stat_posn];
			rhovxsqr = 2.0*rect_height_data[stat_posn+1];
			rhovx = rect_height_data[stat_posn+2];

			(void) fprintf(height_file, FIELD,
				(rhovxsqr - sqr(rhovx)/rho)*hx/(GXU - GXL));

			break;

		case R_YY:
			rho = rect_height_data[stat_posn];
			rhovysqr = 2.0*rect_height_data[stat_posn+1];
			rhovy = rect_height_data[stat_posn+2];

			(void) fprintf(height_file, FIELD,
				(rhovysqr - sqr(rhovy)/rho)*hx/(GXU - GXL));

			break;

		case MIXED_ZONE_FRAC:
			if (first_layer)
			{
				value = sum_stat_vector(rect_height_data,
						  stat_vector_stride,stat_posn,
						  num_stats_up);
				value /= num_stats_up*num_nodes_across;
				(void) fprintf(time_file, FIELD, value);
			}
			(void) fprintf(height_file, FIELD, 
				rect_height_data[stat_posn++]*hx/(GXU - GXL));
			break;
		}
	    }

	    (void) fprintf(height_file,"\n");

	    if (first_layer)
		    (void) fprintf(time_file, "\n");

	    first_layer = NO;
	}

}		/*end printout_rect_state_stats*/

LOCAL double sum_stat_vector(
	double *uni_array,
	int   stride,
	int   offset,
	int   entries)
{
	int i;
	double total = 0.0;

	uni_array += offset;

	for (i = 0; i < entries; i++)
		total += uni_array[i*stride];

	return total;
}		/*end sum_stat_vector*/


LOCAL void print_rss_height_header(
	Rect_state_stat_data	*rssdata,
	Front			*front)
{
	FILE		*height_file = rssdata->rss_height_file;
	int		i;

	(void) foutput(height_file);
	(void) fprintf(height_file,
		  "\t\tBEGIN RECT STATE HEIGHT DEPENDENT STATISTICAL DATA\n\n");
	(void) fprintf(height_file, "\n\n\tt = %g,\tj = %d\n\n",
		       front->time,front->step);

	(void) foutput(height_file);
	(void) fprintf(height_file,HEADERFIELD,"Y");
	for (i = 0; i < rssdata->num_rect_stats; i++)
	{
	    switch (rssdata->rect_stats_selected[i])
	    {
	    case AVG_RHO_CURLV_SQR:
		(void) fprintf(height_file,HEADERFIELD,"H");
		break;

	    case R_XX:
	    case R_YY:
		break;

	    case MIXED_ZONE_FRAC:
		(void) fprintf(height_file,HEADERFIELD,"MZ_FRAC");
		break;

	    case POT_ENG:
		(void) fprintf(height_file,HEADERFIELD,"POT_ENG");
		break;

	    case RHO_SQR:
		(void) fprintf(height_file,HEADERFIELD,"RHO_SQR");
		break;

	    case KE_TOTAL:
		(void) fprintf(height_file,HEADERFIELD,"KE_TOTAL");
		break;

	    case E_TOTAL:
		(void) fprintf(height_file,HEADERFIELD,"E_TOTAL");
		break;

	    case KE_X:
		(void) fprintf(height_file,HEADERFIELD,"KE_X");
		break;

	    case KE_Y:
		(void) fprintf(height_file,HEADERFIELD,"KE_Y");
		break;

	    case RHOX:
		(void) fprintf(height_file,HEADERFIELD,"RHOX");
		break;

	    case RHOY:
		(void) fprintf(height_file,HEADERFIELD,"RHOY");
		break;

	    case RHOVX:
		(void) fprintf(height_file,HEADERFIELD,"RHOVX");
		break;

	    case RHOVY:
		(void) fprintf(height_file,HEADERFIELD,"RHOVY");
		break;

	    default:
		screen("ERROR in print_rss_height_header(), "
		       "bad statatistic chosen\n");
		clean_up(ERROR);
	    }
	}
	(void) fprintf(height_file,"\n");

}		/*end print_rss_height_header*/


LOCAL void print_rss_time_header(
	Rect_state_stat_data	*rssdata)
{
	FILE		*time_file = rssdata->rss_time_file;
	int		i;

	(void) foutput(time_file);
	(void) fprintf(time_file,
		  "\t\tBEGIN RECT STATE TIME DEPENDENT STATISTICAL DATA\n\n");

	(void) foutput(time_file);
	(void) fprintf(time_file,HEADERFIELD,"STEP");
	(void) fprintf(time_file,HEADERFIELD,"TIME");
	for (i = 0; i < rssdata->num_rect_stats; i++)
	{
		switch (rssdata->rect_stats_selected[i])
		{
		case AVG_RHO_CURLV_SQR:
			(void) fprintf(time_file,HEADERFIELD,"H");
			break;

		case R_XX:
		case R_YY:
			break;

		case MIXED_ZONE_FRAC:
			(void) fprintf(time_file,HEADERFIELD,"MZ_FRAC");
			break;

		case POT_ENG:
			(void) fprintf(time_file,HEADERFIELD,"POT_ENG");
			break;

		case RHO_SQR:
			(void) fprintf(time_file,HEADERFIELD,"RHO_SQR");
			break;

		case KE_TOTAL:
			(void) fprintf(time_file,HEADERFIELD,"KE_TOTAL");
			break;

		case E_TOTAL:
			(void) fprintf(time_file,HEADERFIELD,"E_TOTAL");
			break;

		case KE_X:
			(void) fprintf(time_file,HEADERFIELD,"KE_X");
			break;

		case KE_Y:
			(void) fprintf(time_file,HEADERFIELD,"KE_Y");
			break;

		case RHOX:
			(void) fprintf(time_file,HEADERFIELD,"RHOX");
			break;

		case RHOY:
			(void) fprintf(time_file,HEADERFIELD,"RHOY");
			break;

		case RHOVX:
			(void) fprintf(time_file,HEADERFIELD,"RHOVX");
			break;

		case RHOVY:
			(void) fprintf(time_file,HEADERFIELD,"RHOVY");
			break;
		}
	}
	(void) fprintf(time_file,"\n");
}		/*end print_rss_time_header*/

#endif /* defined TWOD || defined THREED */
