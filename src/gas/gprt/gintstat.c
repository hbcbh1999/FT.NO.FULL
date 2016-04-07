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
*				gintstat.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
* 	Routines to sum fluid quantities at material interfaces within layers.
*/


#include <gprt/glayer.h>

	/* LOCAL Function Prototypes */

LOCAL	void record_intfc_stats(Grid*,Wave*,Front*,Printplot*,
				OUTPUT_DATA*,boolean);
LOCAL	void print_intfc_stats(const Grid*,const Front*,OUTPUT_DATA*,boolean);
LOCAL	void accumulate_intfc_stats_planar2d(Grid*,Wave*,Front*,
					     const Intfc_stats*);
LOCAL	void accumulate_intfc_stats_planar3d(Grid*,Wave*,Front*,
					     const Intfc_stats*);
LOCAL	void accumulate_intfc_stats_spherical2d(Grid*,Wave*,Front*,
						const Intfc_stats*);
LOCAL	void accumulate_intfc_stats_spherical3d(Grid*,Wave*,Front*,
						const Intfc_stats*);
LOCAL	void accumulate_intfc_data(const Intfc_data*,Intfc_data*);
LOCAL 	int  compute_crossings(double,double,double,double,double,int*,double*);


/*
*			init_intfc_stats():
*
*	Initializer for computation of interface averages on planar layers.
*/

/*ARGSUSED*/
EXPORT	void init_intfc_stats(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Printplot	*prt)
{
        const RECT_GRID	*rgrid = front->rect_grid;
        const int	dim = rgrid->dim;
	Intfc_stats	*istats;
	char		s[Gets_BUF_SIZE];
	int		i;

	if (dim == 1)
	    return;

	screen("Type 'y' to request interface statistics: ");
	(void) Gets(s);
	if ((s[0] != 'Y') && (s[0] != 'y'))
	    return;

	scalar(&istats,sizeof(Intfc_stats));
	zero_scalar(istats,sizeof(Intfc_stats));

	Output_mode(&istats->odata) = EXACT_TIME;
	Output_time_freq(&istats->odata) = 1.0;
	Output_start_time(&istats->odata) = 0;
	Output_in_binary(&istats->odata) = NO;

	init_output_data(init,&istats->odata,grid,prt,
			 "interface statistics",YES,NO,YES);
	add_user_output_function(record_intfc_stats,&istats->odata,prt);

	screen("Compute interface statistics for\n");
	screen("\t\tplanar ('p', default) or radial geometry ('r'): ");
	(void) Gets(s);
	if (s[0] == 'R' || s[0] == 'r')
	{
	    double dh = (dim == 2 ? min(rgrid->h[0],rgrid->h[1])
			: min(min(rgrid->h[0],rgrid->h[1]),
			      rgrid->h[2]));
	    int result;

	    istats->geom = spherical;
	    screen("Enter the coordinates of the origin: ");
	    if (dim == 2) 
	        result = Scanf("%f %f\n",&istats->origin[0],
			       &istats->origin[1]);
	    else 
	        result = Scanf("%f %f %f\n",&istats->origin[0],
			       &istats->origin[1],&istats->origin[2]);
	    if (result != dim)
	    {
	        screen("ERROR in init_intfc_stats(), "
		       "Insufficient number of coordinates.\n");
		clean_up(ERROR);
	    }

	    /* default radius limits for spherical geometry */

	    istats->h_min = dh;

	    /* As a cheap method for guessing default h_max, use the distance
	       from the origin to the farthest domain boundary */

	    istats->h_max = -HUGE_VAL;
	    for (i = 0; i < dim; i++)
	    {
		double l = max(fabs(rgrid->L[i]-istats->origin[i]),
			      fabs(rgrid->U[i]-istats->origin[i]));
	        istats->h_max = max(l,istats->h_max);
	    }
	    istats->n_layers = (int)((istats->h_max - istats->h_min)/dh) + 1;
	}
        else
	{
            istats->geom = planar;
	    
	    /* defaults for planar geometry */

	    istats->h_min = rgrid->GL[dim-1] + rgrid->h[dim-1];
	    istats->h_max = rgrid->GU[dim-1] - rgrid->h[dim-1];
	    istats->n_layers = (int)rgrid->gmax[dim-1];
	}

	screen("Specify the range of %s and number of layers in the ",
	       (istats->geom == spherical ? "radius" : "height"));
	screen("computation of the\n\tinterface averages ");
	screen("(default = %g %g %d): ", 
	       istats->h_min, istats->h_max, istats->n_layers);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%lf %lf %d", 
			  &istats->h_min, &istats->h_max, &istats->n_layers);
	screen("\n");

	istats->dh = (istats->h_max - istats->h_min)/(istats->n_layers - 1);

	/* TODO: If there is surface tension or mass transfer across the
           interface, then there are interface averages for each distinct
           material.  For now, we assume that pressure and normal velocity are
           continuous across interfaces. */

	uni_array(&istats->data,istats->n_layers,sizeof(Intfc_data));

	/* Create the column header */

	if (dim == 2)
	    (void) sprintf(istats->column_header,
		           "%14s %9s %14s %14s %14s %14s %14s %14s\n",
		           "HEIGHT", "COUNT", "NOR_X", "NOR_Z", "V_DOT_NOR", 
		           "PRE_NOR_X", "PRE_NOR_Z", "PV_DOT_NOR");
	else
	    (void) sprintf(istats->column_header,
		           "%14s %9s %14s %14s %14s %14s "
			   "%14s %14s %14s %14s %14s\n",
		           "HEIGHT", "COUNT", "LENGTH",
			   "NOR_X", "NOR_Y", "NOR_Z",
		           "V_DOT_NOR",
			   "PRE_NOR_X", "PRE_NOR_Y", "PRE_NOR_Z",
		           "PV_DOT_NOR");
}		/*end init_intfc_stats*/

/*
*			record_intfc_stats():
*
*	Main routine for computation of interface averages on planar or
*	spherical layers.  Note that the argument prt is never used, but it
*	should not be removed because this function has to conform to the
*	expected prototype in add_user_output_function() (see
*	init_intfc_stats()).  */

/*ARGSUSED*/
LOCAL	void record_intfc_stats(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*out,
	boolean		about_to_stop)
{
        const Intfc_stats *istats = (Intfc_stats*)out;
	const int	dim = front->rect_grid->dim;
	const int	myid = pp_mynode();
	const int	nn = pp_numnodes();
	const int	n_layers = istats->n_layers;
	const size_t	size_data = sizeof(Intfc_data);

	register int	i, id;

	static Intfc_data*	tmp_intfc_data = NULL;

	start_clock("record_intfc_stats");

	for (i = 0; i < n_layers; i++)
	    zero_scalar(&istats->data[i],size_data);

	/* collect interface totals */

	start_clock("accumulate_intfc_stats");

	switch(dim) {
	case 2:
	    if (istats->geom == spherical)
	      accumulate_intfc_stats_spherical2d(grid,wave,front,istats);
	    else	
		accumulate_intfc_stats_planar2d(grid,wave,front,istats);
	    break;
        case 3:
	    if (istats->geom == spherical)
	        accumulate_intfc_stats_spherical3d(grid,wave,front,istats);

	    else
	        accumulate_intfc_stats_planar3d(grid,wave,front,istats);
	    break;
	}

	if (nn > 1)
	{
	    if (tmp_intfc_data == NULL)
		uni_array(&tmp_intfc_data,n_layers,size_data);

	    if (!is_io_node(myid))
		pp_send(INTFC_SUM_ID+myid,(POINTER)istats->data,
			n_layers*size_data,IO_NODE_ID);
	    else
	    {
		for (id = 0; id < nn; id++)
		{
		    if (id == myid) continue;
		    pp_recv(INTFC_SUM_ID+id,id,(POINTER)tmp_intfc_data,
			    n_layers*size_data);

		    /* add the contribution from node id to totals */

		    for (i = 0; i < n_layers; i++)
		        accumulate_intfc_data(&tmp_intfc_data[i],
					      &istats->data[i]);
		}
	    }
	}  /* if (nn > 1) */

	stop_clock("accumulate_intfc_stats");

	/* All that's left to do is output, which is handled by io_node. */

	if (is_io_node(myid))
	    print_intfc_stats(grid,front,out,about_to_stop);

	stop_clock("record_intfc_stats");
}		/*end record_intfc_stats*/


LOCAL	void print_intfc_stats(
	const Grid	*grid,
	const Front	*front,
	OUTPUT_DATA	*out,
	boolean		about_to_stop)
{
        const Intfc_stats *istats = (Intfc_stats*)out;

	const int	dim = front->rect_grid->dim;
        int 		i, j;
        char 		fname[Gets_BUF_SIZE];

	start_clock("print_intfc_stats");

        (void) sprintf(fname,"%s.ts%s",Output_filename(out),
		       right_flush(grid->step,TSTEP_FIELD_WIDTH));
	Output_file(out) = fopen(fname,"w");
	if (Output_file(out) == NULL)
	{
	    screen("ERROR in print_intfc_stats(), "
		   "Unable to open file %s.\n",fname);
	    clean_up(ERROR);
	}
	print_machine_parameters(Output_file(out));

	print_graph_header(Output_file(out),"INTERFACE",
			   istats->column_header,YES,front);
	    
	for (i = 0; i < istats->n_layers; i++)
	{
	    (void) fprintf(Output_file(out),"%15.5e",
			   istats->h_min+i*istats->dh);
	    (void) fprintf(Output_file(out),"%10d",
			   istats->data[i].count);
	    if (dim == 3)
	      (void) fprintf(Output_file(out),"%15.5e",
			     istats->data[i].nor_fac);
	    for (j = 0; j < dim; j++)
	      (void) fprintf(Output_file(out),"%15.5e",
			     istats->data[i].nor[j]);
	    (void) fprintf(Output_file(out),"%15.5e",
			   istats->data[i].vdotn);
	    for (j = 0; j < dim; j++)
	      (void) fprintf(Output_file(out),"%15.5e",
			     istats->data[i].pre[j]);
	    (void) fprintf(Output_file(out),"%15.5e",
			   istats->data[i].pvdotn);
	    (void) fprintf(Output_file(out),"\n");
	}
	    
	print_graph_footer(Output_file(out),"INTERFACE",about_to_stop);
	fclose(Output_file(out));
	stop_clock("print_intfc_stats");
}		/*end print_intfc_stats*/


/*
*			accumulate_intfc_stats_planar2d():
*
* 	Step through every bond on every contact surface.  For each bond, see 
*	where it crosses on the vertical grid defined in the Intfc_stats 
*	structure.  For each crossing, add to the accumulation for fluid 
*	quantities at the particular z by interpolation on the states at the 
*	start and end points of the bond, and increment the counter for that z.
*
*/

/*ARGSUSED*/
LOCAL	void accumulate_intfc_stats_planar2d(
	Grid		  *grid,
	Wave		  *wave,
	Front		  *front,
	const Intfc_stats *istats)
{
	const int	zdir = front->rect_grid->dim - 1;
	const double 	h_min = istats->h_min;
	const double 	h_max = istats->h_max;
	const double	dh = istats->dh;

	register int	i, j;
	int 		n_crx, tmp_n_crx;
	double		nor[MAXD], v[MAXD], vdotn, pre, nor_x;

	INTERFACE	*intfc = front->interf;
	CURVE		**c;
	BOND		*b;
	double	        *pcrds1, *pcrds2, wgt;

	static int	max_n_crx, *crx_index = NULL;
	static double	*crx_z = NULL;
	static boolean	first = YES;
	static Locstate tmpstl = NULL, tmpstr = NULL, tmpst = NULL;
	static POINT	*tmppt = NULL;
 
	if (first == YES)
	{
		first = NO;

		max_n_crx = 1;
		uni_array(&crx_index,max_n_crx,sizeof(int));
		uni_array(&crx_z,max_n_crx,sizeof(double));
		
		/* Static_point creates storage for 2 states in addition to
		   the coordinates. */

		tmppt = Static_point(intfc);
		tmpstl = left_state(tmppt);
		tmpstr = right_state(tmppt);
	}

	for (c = intfc->curves; c && *c; c++)
	{
	    if ((!is_scalar_wave(wave_type(*c))) && 
		(!is_thinflame_wave(wave_type(*c))))
		continue;
	    
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		if (b->length <= 0.0) continue;

		pcrds1 = Coords(b->start);
		pcrds2 = Coords(b->end);
 
		/* Make sure the bond b is on the actual grid.  By convention,
		   continue if the starting point of b is not on. */

		for (i = 0; i <= zdir; i++)
		    if (pcrds1[i] < front->rect_grid->L[i] ||
		        pcrds1[i] >= front->rect_grid->U[i])
			break;
		if (i <= zdir) continue;

		/* See where the bond b crosses z values on the statistical
		   grid in Intfc_stats.  Then allocate (if necessary) and
		   load an array of z-index values where this bond crosses. */

		tmp_n_crx = (int)(fabs(pcrds2[zdir]-pcrds1[zdir])/dh) + 2;

		if (tmp_n_crx > max_n_crx)
		{
		    max_n_crx = tmp_n_crx;
		    free(crx_index);
		    free(crx_z);
		    uni_array(&crx_index,max_n_crx,sizeof(int));
		    uni_array(&crx_z,max_n_crx,sizeof(double));
		}
			
		n_crx = compute_crossings(pcrds1[zdir],pcrds2[zdir],
					  h_min,h_max,dh,crx_index,crx_z);

		if (n_crx > max_n_crx)
		{
		    screen("ERROR in accumulate_intfc_stats(), "
		           "number of crossings, %d, exceeds "
		           "size of crossing array, %d.\n",n_crx,max_n_crx);
		    clean_up(ERROR);
		}

		for (i = 0; i < n_crx; i++)
		{
		    /* find the relative coordinate of this crossing 
		       along the bond */
		    wgt = (crx_z[i] - pcrds1[zdir])
			/ (pcrds2[zdir] - pcrds1[zdir]);
		  
		    /* get the absolute coordinates of and normal uni_array
		       at this crossing */

		    Coords(tmppt)[zdir] = crx_z[i];
		    for (j = 0; j < zdir; j++)
			Coords(tmppt)[j] = pcrds1[j] + 
			    wgt * (pcrds2[j] - pcrds1[j]);
		    normal(tmppt,Hyper_surf_element(b),Hyper_surf(*c),
			   nor,front);
		    nor_x = fabs(nor[0]);

		    /* get the states at this crossing */
		    left_state_along_bond(wgt,b,*c,tmpstl);
		    right_state_along_bond(wgt,b,*c,tmpstr);

		    /* NOTE: Use the convention that the interface states are
		       sampled from the material corresponding to params index
		       0.  Thus we need to make sure that nor points INTO this
		       fluid (it points from left state to right state by
		       default). */
		
		    if (index_of_Gas_param(Params(tmpstl)) == 0)
		    {
		        tmpst = tmpstl;
			for (j = 0; j <= zdir; j++)
			    nor[j] /= -nor_x;
		    }
		    else
		    {	
		        tmpst = tmpstr;
			for (j = 0; j <= zdir; j++)
			    nor[j] /= nor_x;
		        if (index_of_Gas_param(Params(tmpstr)) != 0)
			{
			    screen("ERROR in accumulate_intfc_stats(), "
				   "Can't match state at bond.\n");
			    clean_up(ERROR);
			}
		    }   
		    for (j = 0; j <= zdir; j++)
		        v[j] = vel(j,tmpst);
		    pre = pressure(tmpst);
		    vdotn = scalar_product(v,nor,zdir+1);

		    /* add in the fluid quantities at this crossing
		       and increment count */

		    istats->data[crx_index[i]].count++;
		    istats->data[crx_index[i]].nor_fac += 1.0;
		    istats->data[crx_index[i]].vdotn += vdotn;
		    istats->data[crx_index[i]].pvdotn += pre*vdotn;
		    for (j = 0; j <= zdir; j++)
		    {
			istats->data[crx_index[i]].nor[j] += nor[j];
			istats->data[crx_index[i]].pre[j] += pre*nor[j];
		    }
		    
		}  /* for (i = 0; i < n_crx; i++) */

	    }  /* for (b = (*c)->first; ...) */

	}  /* for (c = intfc->curves; ...) */

}		/* end accumulate_intfc_stats_planar2d() */


/*
*		accumulate_intfc_stats_planar3d():
*
* 	Step through every triangle on every contact surface.  For each 
*	triangle, see where it crosses on the vertical grid defined in the 
*	Intfc_stats structure.  For each crossing, find the line segment of
*       intersection of the triangle with the (x,y) plane, get the state at
*       the midpoint of this segment by interpolation on the states at the
*       vertices of the triangle, add to the accumulation for fluid 
*       quantities at the particular z, weighting the contribution by the
*       length of the segment, and add the length to the normalization factor.
*/

/*ARGSUSED*/
LOCAL	void accumulate_intfc_stats_planar3d(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	const Intfc_stats *istats)
{
	RECT_GRID   *rgrid = front->rect_grid;
	double       *L = rgrid->U, *U = rgrid->U;
	const int   zdir = rgrid->dim - 1;
	const double h_min = istats->h_min;
	const double h_max = istats->h_max;
	const double dh = istats->dh;

	register int	i, j;
	int 		n_crx, tmp_n_crx;
	double		nor[MAXD], v[MAXD], vdotn, pre, nor_xy;

	INTERFACE	*intfc = front->interf;
	SURFACE		**s;
	TRI		*t;
	double	        *pcrds[3], pcrds1[MAXD], pcrds2[MAXD];
	double 		wgt1, wgt2, nor_len, length;
	int		toppt, botpt, midpt;
	Locstate 	topstl, topstr, botstl, botstr, midstl, midstr;

	static int	max_n_crx, *crx_index = NULL;
	static double	*crx_z = NULL;
	static boolean	first = YES;
	static Locstate tmpstl = NULL, tmpstr = NULL, tmpst = NULL;
	static POINT	*tmppt = NULL;
 
	if (first == YES)
	{
		first = NO;

		max_n_crx = 1;
		uni_array(&crx_index,max_n_crx,sizeof(int));
		uni_array(&crx_z,max_n_crx,sizeof(double));
		
		/* Static_point creates storage for 2 states in addition to
		   the coordinates. */

		tmppt = Static_point(intfc);
		tmpstl = left_state(tmppt);
		tmpstr = right_state(tmppt);
	}

	for (s = intfc->surfaces; s && *s; s++)
	{
	    if ((!is_scalar_wave(wave_type(*s))) &&
		(!is_thinflame_wave(wave_type(*s))))
		continue;
	
	    for (t = first_tri(*s); t != last_tri(*s); t = t->next)
	    {
		for (i = 0; i < 3; i++)
		    pcrds[i] = Coords(Point_of_tri(t)[i]);
 
		/* Make sure the triangle t is on the ACTUAL grid.  By
		   convention, continue if the first vertex of t is not on. */

		for (i = 0; i <= zdir; i++)
		    if (pcrds[0][i] < L[i] || pcrds[0][i] >= U[i])
			break;
		if (i <= zdir) continue;

		/* find the indices of the bottom (smallest z coord), top
		   (largest z coord), and middle vertices along the z dir */

		for (i = 1, toppt = 0, botpt = 0; i < 3; i++)
		{
		    if (pcrds[i][zdir] >= pcrds[toppt][zdir])
			toppt = i;
		    if (pcrds[i][zdir] < pcrds[botpt][zdir])
			botpt = i;
		}
		midpt = (botpt + 1) % 3;
		if (midpt == toppt)
		    midpt = (midpt + 1) % 3;

		if (botpt == midpt || botpt == toppt || midpt == toppt)
		{
		    screen("ERROR in accumulate_intfc_stats(), "
		           "bottom (%d), middle (%d) and top (%d) vertices "
		           "of triangle are not distinct.\n",
			   botpt,midpt,toppt);
		    clean_up(ERROR);
		}

		/* See where the triangle t crosses z values on the statistical
		   grid in Intfc_stats.  Then allocate (if necessary) and load
		   an array of z-index values where this triangle crosses. */

		tmp_n_crx = (int)(fabs(pcrds[toppt][zdir]
				       - pcrds[botpt][zdir]) / dh) + 2;
		if (tmp_n_crx > max_n_crx)
		{
		    max_n_crx = tmp_n_crx;
		    free(crx_index);
		    free(crx_z);
		    uni_array(&crx_index,max_n_crx,sizeof(int));
		    uni_array(&crx_z,max_n_crx,sizeof(double));
		}
			
		n_crx = compute_crossings(pcrds[botpt][zdir],pcrds[toppt][zdir],
					  h_min,h_max,dh,crx_index,crx_z);

		if (n_crx > max_n_crx)
		{
		    screen("ERROR in accumulate_intfc_stats(), "
		           "number of crossings, %d, exceeds "
		           "size of crossing array, %d.\n",n_crx,max_n_crx);
		    clean_up(ERROR);
		}

		if (n_crx > 0)
		{
		    /* we'll need the left and right states on each vertex
		       of the triangle t */

		    slsr(Point_of_tri(t)[botpt],Hyper_surf_element(t),
			 Hyper_surf(*s),&botstl,&botstr);
		    slsr(Point_of_tri(t)[midpt],Hyper_surf_element(t),
			 Hyper_surf(*s),&midstl,&midstr);
		    slsr(Point_of_tri(t)[toppt],Hyper_surf_element(t),
			 Hyper_surf(*s),&topstl,&topstr);
		}
		else
		    continue;

		for (i = 0; i < n_crx; i++)
		{
		    const double *tnor;
		    /* Find the line segment formed by the intersection of this
		       triangle with the (x,y) plane at z = crx_z[i].  One
		       endpoint of this segment is always on the (botpt,toppt)
		       side; the other endpoint is on either the (botpt,midpt)
		       side or the (midpt,toppt) side depending on the z
		       coordinate of the midpt vertex. */

		    wgt1 = (crx_z[i] - pcrds[botpt][zdir]) 
			/ (pcrds[toppt][zdir] - pcrds[botpt][zdir]);

		    if (pcrds[midpt][zdir] > crx_z[i])
			wgt2 = (crx_z[i] - pcrds[botpt][zdir])
			    / (pcrds[midpt][zdir] - pcrds[botpt][zdir]);
		    else
			wgt2 = (crx_z[i] - pcrds[midpt][zdir])
			    / (pcrds[toppt][zdir] - pcrds[midpt][zdir]);
		  
		    /* get the absolute coordinates of the endpoints of
		       the line segment and its length */

		    pcrds1[zdir] = pcrds2[zdir] = crx_z[i];
		    for (j = 0; j < zdir; j++)
		    {
			pcrds1[j] = pcrds[botpt][j] 
			    + wgt1 * (pcrds[toppt][j] - pcrds[botpt][j]);
			if (pcrds[midpt][zdir] > crx_z[i])
			    pcrds2[j] = pcrds[botpt][j] 
				+ wgt2 * (pcrds[midpt][j] - pcrds[botpt][j]);
			else
			    pcrds2[j] = pcrds[midpt][j]
				+ wgt2 * (pcrds[toppt][j] - pcrds[midpt][j]);
		    }
		    length = sqrt(sqr(pcrds2[0]-pcrds1[0])
				  + sqr(pcrds2[1]-pcrds1[1]));

		    /* if the triangle is just touching the (x,y) plane at
		       z = crx_z[i], then the length of the intersection is
		       zero, hence no contribution to the interfacial avg. */

		    if (length <= 0.0)
		        continue;
		        
		    /* get the state and normal vector at the midpoint 
		       of the line segment */

		    Coords(tmppt)[zdir] = crx_z[i];
		    for (j = 0; j < zdir; j++)
			Coords(tmppt)[j] = 0.5 * (pcrds1[j] + pcrds2[j]);
		    tnor = Tri_normal(t);
		    nor_len = Mag3d(tnor);
		    for (j = 0; j <= zdir; j++)
		      nor[j] = tnor[j]/nor_len;
		    nor_xy = sqrt(nor[0]*nor[0] + nor[1]*nor[1]);

		    if (nor_xy <= 0.0)
		    {
		        (void) printf("WARNING in accumulate_intfc_stats(), "
			              "normal to triangle at z = %g "
			              "is aligned along the z direction\n",
				      crx_z[i]);
			(void) printf("\tnor[0] = %g, nor[1] = %g, "
				      "nor[2] = %g\n",nor[0],nor[1],nor[2]);
			(void) printf("\ttop: pcrds[%d] = (%g %g %g)\n",toppt,
			              pcrds[toppt][0],pcrds[toppt][1],
				      pcrds[toppt][2]);
			(void) printf("\tmid: pcrds[%d] = (%g %g %g)\n",midpt,
			              pcrds[midpt][0],pcrds[midpt][1],
				      pcrds[midpt][2]);
			(void) printf("\tbot: pcrds[%d] = (%g %g %g)\n",botpt,
			              pcrds[botpt][0],pcrds[botpt][1],
				      pcrds[botpt][2]);
		    }

		    if (pcrds[midpt][zdir] > crx_z[i])
		    {
			if (tri_interpolate_states(front,1.0-0.5*(wgt1+wgt2),
					           0.5*wgt2,0.5*wgt1,NULL,
						   botstl,NULL,midstl,NULL,
						   topstl,tmpstl)
						       != FUNCTION_SUCCEEDED)
			{
			    screen("ERROR in "
				   "accumulate_intfc_stats_planar3d(), "
				   "tri_interpolate_states() failed on "
				   "the negative side for "
				   "pcrds[%d][%d] > crx_z[%d]\n",midpt,zdir,i);
			    clean_up(ERROR);
			}
			if (tri_interpolate_states(front,1.0-0.5*(wgt1+wgt2),
						   0.5*wgt2,0.5*wgt1,NULL,
						   botstr,NULL,midstr,NULL,
						   topstr,tmpstr)
						       != FUNCTION_SUCCEEDED)
			{
			    screen("ERROR in "
				   "accumulate_intfc_stats_planar3d(), "
				   "tri_interpolate_states() failed on "
				   "the positive side for "
				   "pcrds[%d][%d] > crx_z[%d]\n",midpt,zdir,i);
			    clean_up(ERROR);
			}
		    }
		    else
		    {
			if (tri_interpolate_states(front,0.5*(1.0-wgt1),
						   0.5*(1.0-wgt2),
						   0.5*(wgt1+wgt2),NULL,
						   botstl,NULL,midstl,NULL,
					           topstl,tmpstl)
						       != FUNCTION_SUCCEEDED)
			{
			    screen("ERROR in "
				   "accumulate_intfc_stats_planar3d(), "
				   "tri_interpolate_states() failed on "
				   "the negative side for "
				   "pcrds[%d][%d] <= crx_z[%d]\n",midpt,zdir,i);
			    clean_up(ERROR);
			}
			if (tri_interpolate_states(front,0.5*(1.0-wgt1),
					           0.5*(1.0-wgt2),
						   0.5*(wgt1+wgt2),
					           NULL,botstr,NULL,midstr,NULL,
					           topstr,tmpstr)
						       != FUNCTION_SUCCEEDED)
			{
			    screen("ERROR in "
				   "accumulate_intfc_stats_planar3d(), "
				   "tri_interpolate_states() failed on "
				   "the positive side for "
				   "pcrds[%d][%d] <= crx_z[%d]\n",midpt,zdir,i);
			    clean_up(ERROR);
			}
		    }

		    /* NOTE: Use the convention that the interface states are
		       sampled from the material corresponding to params index
		       0.  Thus we need to make sure that nor points INTO this
		       fluid (it points from left state to right state by
		       default). */

		    if (index_of_Gas_param(Params(midstl)) == 0)
		    { 	
		        tmpst = tmpstl;
			for (j = 0; j <= zdir; j++)
			    nor[j] /= -nor_xy;
		    }
		    else
		    {
		        tmpst = tmpstr;
			for (j = 0; j <= zdir; j++)
			    nor[j] /= nor_xy;
		        if (index_of_Gas_param(Params(midstr)) != 0)
			{
			    screen("ERROR in accumulate_intfc_stats(), "
				   "Can't match state at bond.\n");
			    clean_up(ERROR);
			}
		    }
		    for (j = 0; j <= zdir; j++)
		        v[j] = vel(j,tmpst);
		    pre = pressure(tmpst);
		    vdotn = scalar_product(v,nor,zdir+1);

		    /* add in the fluid quantities at this crossing and
		       increment count and nor_fac */

		    istats->data[crx_index[i]].count++;
		    istats->data[crx_index[i]].nor_fac += length;
		    istats->data[crx_index[i]].vdotn += length*vdotn;
		    istats->data[crx_index[i]].pvdotn += length*pre*vdotn;
		    for (j = 0; j <= zdir; j++)
		    {
			istats->data[crx_index[i]].nor[j] += length*nor[j];
			istats->data[crx_index[i]].pre[j] 
			    += length*pre*nor[j];
		    }

		}  /* for (i = 0; i < n_crx; i++) */

	    }  /* for (t = first_tri(s); ...) */

	}  /* for (s = intfc->surfaces; ...) */

}		/* end accumulate_intfc_stats_planar3d() */


/*ARGSUSED*/
LOCAL	void accumulate_intfc_stats_spherical2d(
	Grid		  *grid,
	Wave		  *wave,
	Front		  *front,
	const Intfc_stats *istats)
{
	screen("ERROR in accumulate_intfc_stats_spherical2d(), "
	       "Interface stats are not yet available "
	       "for 3D spherical layers.\n");
	clean_up(ERROR);
}		/* end accumulate_intfc_stats_spherical2d() */


/*ARGSUSED*/
LOCAL	void accumulate_intfc_stats_spherical3d(
	Grid		  *grid,
	Wave		  *wave,
	Front		  *front,
	const Intfc_stats *istats)
{
	screen("ERROR in accumulate_intfc_stats_spherical3d(), "
	       "Interface stats are not yet available "
	       "for 3D spherical layers.\n");
	clean_up(ERROR);
}		/* end accumulate_intfc_stats_spherical3d() */


LOCAL	void	accumulate_intfc_data(
	const Intfc_data* idata1,
	Intfc_data* idata2)
{
        int j;
	for (j = 0; j < MAXD; j++)
	{
	    idata2->nor[j] += idata1->nor[j];
	    idata2->pre[j] += idata1->pre[j];
	}
	idata2->vdotn += idata1->vdotn;
	idata2->pvdotn += idata1->pvdotn;
	idata2->nor_fac += idata1->nor_fac;
	idata2->count += idata1->count;
}		/*end accumulate_intfc_data*/


/*
*			compute_crossings():
*
*       Find all of the crossings, including endpoints, of the segment [z1,z2] 
*	in the 1D grid defined on the closed interval [z_min,z_max] with spacing
*	dz.  For each crossing, the index of the grid point and its associated
*	z value are stored in the arrays index and z_value, respectively. 
*/

/*ARGSUSED*/
LOCAL 	int compute_crossings(
	double 	z1,
	double 	z2,
        double	z_min,
        double 	z_max,
        double 	dz,
        int 	*index,
        double	*z_value)
{
        double z;
	int n, n1, n2, count = 0;
	
	if (z1 > z2)
	{
	    z = z2; z2 = z1; z1 = z;
	}
	if (z_min > z_max)
	{
	    z = z_max; z_max = z_min; z_min = z;
	}
	if (z2 < z_min || z1 > z_max)
	    return 0;
	
	if (z1 < z_min)
	    z1 = z_min;
	if (z2 > z_max)
	    z2 = z_max;
	
	n1 = irint((z1 - z_min)/dz);
	n2 = irint((z2 - z_min)/dz) + 1;
	
	for (n = n1; n <= n2; n++)
	{
	    z = z_min + n*dz;
	    if ((z >= z1) && (z <= z2))
	    {
		index[count] = n;
		z_value[count] = z;
		count++;
	    }
	}
	return count;
}		/*end compute_crossings*/
