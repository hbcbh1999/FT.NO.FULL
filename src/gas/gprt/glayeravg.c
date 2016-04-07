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
*				glayeravg.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains routines to average gas states in various layers
*
*/

#if defined(TWOD) || defined(THREED)

#include <gprt/glayer.h>

LOCAL	void init_column_headers(Layer_stats*,int);
LOCAL	void record_layer_stats(Grid*,Wave*,Front*,Printplot*,
				OUTPUT_DATA*,boolean);
LOCAL	void print_layer_stats(const Grid*,const Front*,OUTPUT_DATA*,boolean);


/*			init_layer_stats()
*
*	Initializer for computation of layer averages.
*/

/*ARGSUSED*/
EXPORT	void init_layer_stats(
	INIT_DATA	*init,
	Front		*front,
	Grid		*grid,
	Printplot	*prt)
{
	INTERFACE	*intfc = front->interf;
	const RECT_GRID *rgrid = front->rect_grid;
	const size_t	size_bst = sizeof(Big_State);
	const size_t	size_turb = sizeof(Turbulence_moments);
        const int	dim = rgrid->dim;
	const int	n_params = num_gas_params(intfc);
	Layer_stats	*lstats;
	char		s[Gets_BUF_SIZE];
	int		i;
	
	if (dim == 1)
	    return;

	screen("Type 'y' to request layer statistics: ");
	(void) Gets(s);
	if ((s[0] != 'Y') && (s[0] != 'y'))
	    return;

	scalar(&lstats,sizeof(Layer_stats));
	zero_scalar(lstats,sizeof(Intfc_stats));

	uni_array(&lstats->bst_out,n_params,sizeof(FILE*));
	uni_array(&lstats->bst_fname,n_params,sizeof(char*));
	uni_array(&lstats->bst_dname,n_params,sizeof(char*));

	/* defaults */
	Output_mode(&lstats->odata) = EXACT_TIME;
	Output_time_freq(&lstats->odata) = 1.0;
	Output_start_time(&lstats->odata) = 0;
	Output_in_binary(&lstats->odata) = NO;

	init_output_data(init,&lstats->odata,grid,prt,NULL,NO,NO,YES);
	add_user_output_function(record_layer_stats,&lstats->odata,prt);

	lstats->n_params = n_params;    /* number of distinct materials */
        lstats->rfactor = 2.0;		/* default sub-grid refinement */

        screen("Enter a sub-grid refinement factor for the\n"
               "\tpointwise averaging (default = %g): ",lstats->rfactor);
        (void) Gets(s);
        if (strlen(s) != 0)
	    (void) sscan_float(s,&lstats->rfactor);

	screen("Compute layer statistics for\n");
	screen("\t\tplanar ('p', default) or radial geometry ('r'): ");
	(void) Gets(s);
	if (s[0] == 'R' || s[0] == 'r')
	{
	    double dh = (dim == 2 ? min(rgrid->h[0],rgrid->h[1])/lstats->rfactor
			: min(min(rgrid->h[0],rgrid->h[1]),
			      rgrid->h[2])/lstats->rfactor);
	    int result;

	    lstats->geom = spherical;
	    screen("Enter the coordinates of the origin: ");
	    if (dim == 2) 
	        result = Scanf("%f %f\n",&lstats->origin[0],
			       &lstats->origin[1]);
	    else 
	        result = Scanf("%f %f %f\n",&lstats->origin[0],
			       &lstats->origin[1],&lstats->origin[2]);
	    if (result != dim)
	    {
		screen("ERROR in init_layer_stats(), "
		       "Insufficient number of coordinates.\n");
		clean_up(ERROR);
	    }

	    /* default radius limits for spherical geometry */

	    lstats->h_min = dh;

	    /* As a cheap method for guessing default h_max, use the distance
	       from the origin to the farthest domain boundary */

	    lstats->h_max = -HUGE_VAL;
	    for (i = 0; i < dim; ++i)
	    {
		double l = max(fabs(rgrid->GL[i]-lstats->origin[i]),
			      fabs(rgrid->GU[i]-lstats->origin[i]));
		if (l > lstats->h_max)
		    lstats->h_max = l;
	    }
	    lstats->n_layers = (int)((lstats->h_max - lstats->h_min)/dh) + 1;
	}
        else
	{
            lstats->geom = planar;
	    
	    /* defaults for planar geometry */

	    lstats->h_min = rgrid->GL[dim-1] + rgrid->h[dim-1]/lstats->rfactor;
	    lstats->h_max = rgrid->GU[dim-1] - rgrid->h[dim-1]/lstats->rfactor;
	    lstats->n_layers = (int)(rgrid->gmax[dim-1]*lstats->rfactor);
	}

	screen("Specify the range of %s and number of layers in the ",
	       (lstats->geom == spherical ? "radius" : "height"));
	screen("computation of the\n\tlayer averages ");
	screen("(default = %g %g %d): ", 
	       lstats->h_min, lstats->h_max, lstats->n_layers);
	(void) Gets(s);
	if (s[0] != '\0')
	    (void) sscanf(s,"%lf %lf %d", 
			  &lstats->h_min, &lstats->h_max, &lstats->n_layers);
	screen("\n");

        uni_array(&lstats->bst_data,lstats->n_layers*n_params,size_bst);

	lstats->dh = (lstats->h_max - lstats->h_min)/(lstats->n_layers - 1);

	for (i = 0; i < n_params; ++i)
	{
	    (void) sprintf(s,"fluid %d",i);
	    lstats->bst_out[i] = open_data_file(front,s,NO,YES,
						NULL,&lstats->bst_dname[i],
						NULL,&lstats->bst_fname[i]);
	}

	screen("Type 'y' to include turbulence moments in the");
	screen(" layer averaging: ");
	(void) Gets(s);
	lstats->include_turb = (s[0] == 'Y' || s[0] == 'y' ? YES : NO);

	if (lstats->include_turb == YES)
	{
	    uni_array(&lstats->turb_data,lstats->n_layers*n_params,size_turb);
	    uni_array(&lstats->turb_out,n_params,sizeof(FILE*));
	    uni_array(&lstats->turb_fname,n_params,sizeof(char*));
	    uni_array(&lstats->turb_dname,n_params,sizeof(char*));
	    for (i = 0; i < n_params; ++i)
	    {
		(void) sprintf(s,"fluid %d",i);
		lstats->turb_out[i] = open_data_file(front,s,NO,YES,
						  NULL,&lstats->turb_dname[i],
						  NULL,&lstats->turb_fname[i]);
	    }
	}

	init_column_headers(lstats,dim);
}		/*end init_layer_stats*/


LOCAL	void init_column_headers(
	Layer_stats* lstats,
	int dim)
{
	char	xyz[3] = {'X', 'Y', 'Z'};
	char	rtp[3] = {'R', 'T', 'P'};
        char 	tmpstr[Gets_BUF_SIZE];
	int i, j;

	/* First the first-order moments: */

	(void) sprintf(lstats->bst_col_header,"%14s %14s %14s",
		       (lstats->geom == planar ? "HEIGHT" : "RADIUS"), 
		       "FRAC", "DEN");
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(tmpstr," %12s_%c","VEL",
		           (lstats->geom == planar ? xyz[i] : rtp[i]));
	    (void) strcat(lstats->bst_col_header,tmpstr);
	}
	(void) sprintf(tmpstr," %14s %14s %14s\n", "PRE", "KE", "IE");
	(void) strcat(lstats->bst_col_header,tmpstr);

	/* Now the second-order (turbulent) moments: */

	(void) sprintf(lstats->turb_col_header,"%14s",
		(lstats->geom == planar ? "HEIGHT" : "RADIUS"));
	(void) sprintf(tmpstr," %14s %14s %14s","DEN*DEN","DEN*KE","DEN*IE");
	(void) strcat(lstats->turb_col_header,tmpstr);
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(tmpstr," %12s_%c","DEN*VEL",
		    (lstats->geom == planar ? xyz[i] : rtp[i]));
	    (void) strcat(lstats->turb_col_header,tmpstr);
	}
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(tmpstr," %12s_%c","DKV",
		    (lstats->geom == planar ? xyz[i] : rtp[i]));
	    (void) strcat(lstats->turb_col_header,tmpstr);
	}
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(tmpstr," %12s_%c","DEV",
		    (lstats->geom == planar ? xyz[i] : rtp[i]));
	    (void) strcat(lstats->turb_col_header,tmpstr);
	}
	for (i = 0; i < dim; ++i)
	{
	    (void) sprintf(tmpstr," %12s_%c","PRE*VEL",
		    (lstats->geom == planar ? xyz[i] : rtp[i]));
	    (void) strcat(lstats->turb_col_header,tmpstr);
	}
	for (i = 0; i < dim; ++i)
	{
	    for (j = i; j < dim; ++j)
	    {
	        (void) sprintf(tmpstr," %11s_%c%c","DVV",
			(lstats->geom == planar ? xyz[i] : rtp[i]),
			(lstats->geom == planar ? xyz[j] : rtp[j]));
		(void) strcat(lstats->turb_col_header,tmpstr);
	    }
	}
	for (i = 0; i < dim; ++i)
	{
	    for (j = i; j < dim; ++j)
	    {
	        (void) sprintf(tmpstr," %11s_%c%c","VV",
			(lstats->geom == planar ? xyz[i] : rtp[i]),
			(lstats->geom == planar ? xyz[j] : rtp[j]));
		(void) strcat(lstats->turb_col_header,tmpstr);
	    }
	}
	(void) sprintf(tmpstr,"\n");
	(void) strcat(lstats->turb_col_header,tmpstr);
}		/*end init_column_headers*/

/*			record_layer_stats()
*
*	Main routine for computation of averages in planar or spherical layers.
*	Note that the argument prt is never used, but it should not be removed
*	because this function has to conform to the expected prototype in
*	add_user_output_function() (see init_intfc_stats()).
*/

/*ARGSUSED*/
LOCAL	void record_layer_stats(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*out,
	boolean		about_to_stop)
{
        const Layer_stats *lstats = (Layer_stats*)out;

	const double	*origin = (lstats->geom == spherical ? 
				   lstats->origin : NULL);
	const boolean	include_turb = lstats->include_turb;
	const int	n_layers = lstats->n_layers;
	const int	n_params = lstats->n_params;
	const int	nn = pp_numnodes();
	const int	myid = pp_mynode();
	const size_t	size_bst = sizeof(Big_State);
	const size_t	size_turb = sizeof(Turbulence_moments);

	Big_State* 		bst_data = lstats->bst_data;
	Turbulence_moments* 	turb_data = lstats->turb_data;

	double 		h;
	register int	i;
	int		n, p;

	static Big_State* 		tmp_bst = NULL;
	static Turbulence_moments* 	tmp_turb = NULL;

	debug_print("glayer","Entered record_layer_stats()\n");
        start_clock("record_layer_stats");

	/* Accumulation and normalization of layer totals are done by the
	   I/O node, which therefore needs additional storage for the data
	   which it receives from the other nodes. */

	if (is_io_node(myid) && nn > 1 && tmp_bst == NULL)
	{
	    uni_array(&tmp_bst,n_layers*n_params,size_bst);
	    if (include_turb == YES)
	        uni_array(&tmp_turb,n_layers*n_params,size_turb);
	}

	/* Accumulate layer totals local to this node */

	start_clock("accumulate_layer_stats");

	switch(include_turb)
	{
	case YES:
	    for (i = 0; i < n_layers; ++i)
	    {
	        h = lstats->h_min + i*lstats->dh;
		accumulate_state_in_layer(wave,front,h,
					  &bst_data[n_params*i],
					  &turb_data[n_params*i],
					  lstats->rfactor,origin,
					  NO,NO,NO,NO);
	    }
	    break;
	case NO:
	    for (i = 0; i < n_layers; ++i)
	    {
	        h = lstats->h_min + i*lstats->dh;
		accumulate_state_in_layer(wave,front,h,
					  &bst_data[n_params*i],NULL,
					  lstats->rfactor,origin,
					  NO,NO,NO,NO);
	    }
	    break;
	}

	/* Send layer totals to the I/O node */

	if (nn > 1)
	{
	    if (is_io_node(myid))
	    {
	        for (n = 0; n < nn; ++n)
		{
		    if (n != myid)
		    {
		        pp_recv(LAYER_SUM_ID,n,(POINTER)tmp_bst,
				n_params*n_layers*size_bst);
			if (debugging("glayer"))
			{
			    (void) printf("Node %d:\n",n);
			    (void) printf("\tlayer 0 count = %d and %d\n",
				          tmp_bst[0].count, tmp_bst[1].count);
			    (void) printf("\tlayer %d count = %d and %d\n",
				          n_layers/2,
					  tmp_bst[2*(n_layers/2)].count,
				          tmp_bst[2*(n_layers/2)+1].count);
			    (void) printf("\tlayer %d count = %d and %d\n",
				          n_layers-1,
					  tmp_bst[2*(n_layers-1)].count,
				          tmp_bst[2*(n_layers-1)+1].count);
			}
			for (i = 0; i < n_layers; ++i)
			    for (p = 0; p < n_params; ++p)
			        accumulate_state_totals(&tmp_bst[n_params*i+p],
							NULL,
							&bst_data[n_params*i+p],
							NULL);
			if (include_turb == YES)
			{
			    pp_recv(LAYER_SUM_ID+1,n,(POINTER)tmp_turb,
				    n_params*n_layers*size_turb);
			    if (debugging("glayer"))
			    {
				(void) printf("\tlayer 0 dd sum = %g and %g\n",
					      tmp_turb[0].dd,tmp_turb[1].dd);
				(void) printf("\tlayer %d dd sum = %g and %g\n",
					      n_layers/2,
					      tmp_turb[2*(n_layers/2)].dd,
					      tmp_turb[2*(n_layers/2)+1].dd);
				(void) printf("\tlayer %d dd sum = %g and %g\n",
					      n_layers-1,
					      tmp_turb[2*(n_layers-1)].dd,
					      tmp_turb[2*(n_layers-1)+1].dd);
			    }
			    for (i = 0; i < n_layers; ++i)
			        for (p = 0; p < n_params; ++p)
				    accumulate_state_totals(NULL,
						&tmp_turb[n_params*i+p],
						NULL,
						&turb_data[n_params*i+p]);
			}
		    }
		}
	    }
	    else
	    {
	        pp_send(LAYER_SUM_ID,(POINTER)bst_data,
			n_params*n_layers*size_bst,IO_NODE_ID);
		if (include_turb == YES)
		    pp_send(LAYER_SUM_ID+1,(POINTER)turb_data,
			    n_params*n_layers*size_turb,IO_NODE_ID);
	    }
	}

	if (debugging("glayer"))
	    (void) printf("Global accumulation of layer totals is finished.\n");

	stop_clock("accumulate_layer_stats");

	/* Normalize totals and print from the I/O node. */

	if (is_io_node(myid))
	{
	    int n_points;
	    for (i = 0; i < n_layers; ++i)
	    {
	        n_points = 0;
		for (p = 0; p < n_params; ++p)
		    n_points += bst_data[n_params*i+p].count;
		if (debugging("glayer"))
		    (void) printf("In record_layer_stats(), n_points = %d\n",
				  n_points);
		if (include_turb == YES)
		    for (p = 0; p < n_params; ++p)
		        normalize_state_totals(&bst_data[n_params*i+p],
					       &turb_data[n_params*i+p],
					       n_points);
		else
		    for (p = 0; p < n_params; ++p)
		        normalize_state_totals(&bst_data[n_params*i+p],
					       NULL,n_points);
	    }
	    print_layer_stats(grid,front,out,about_to_stop);
	}

        stop_clock("record_layer_stats");
	debug_print("glayer","Left record_layer_stats()\n");
}		/*end record_layer_stats*/


LOCAL	void print_layer_stats(
	const Grid	*grid,
	const Front	*front,
	OUTPUT_DATA	*out,
	boolean		about_to_stop)
{
        const Layer_stats *lstats = (Layer_stats*)out;

	const int	dim = front->rect_grid->dim;
        const int	n_params = lstats->n_params;
	const int	n_layers = lstats->n_layers;

	const Big_State* 	  bst_data = lstats->bst_data;
	const Turbulence_moments* turb_data = lstats->turb_data;

        int 		i, j, k, p;
	double 		h, dk;
        char 		tmpstr[Gets_BUF_SIZE], fname[Gets_BUF_SIZE];

	start_clock("print_layer_stats");
  
	for (p = 0; p < n_params; ++p)
	{
	    (void) sprintf(fname,"%s.ts%s",lstats->bst_fname[p],
		    right_flush(grid->step,TSTEP_FIELD_WIDTH));
	    lstats->bst_out[p] = fopen(fname,"w");
	    if (lstats->bst_out[p] == NULL)
	    {
	        screen("\nERROR in print_layer_stats():\n");
		screen("Unable to open file %s.\n",fname);
		clean_up(ERROR);
	    }
	    print_machine_parameters(lstats->bst_out[p]);
						
	    (void) sprintf(tmpstr,"FLUID %d LAYER",p);
	    print_graph_header(lstats->bst_out[p],tmpstr,
			       lstats->bst_col_header,YES,front);

	    for (i = 0; i < n_layers; ++i)
	    {
		h = lstats->h_min + i*lstats->dh;
		(void) fprintf(lstats->bst_out[p],"%15.5e%15.5e%15.5e",h,
			       bst_data[n_params*i+p].frac,
			       bst_data[n_params*i+p].d);
		for (j = 0; j < dim; ++j)
		    (void) fprintf(lstats->bst_out[p],"%15.5e",
				   bst_data[n_params*i+p].v[j]);
		(void) fprintf(lstats->bst_out[p],"%15.5e%15.5e%15.5e",
			       bst_data[n_params*i+p].p,
			       bst_data[n_params*i+p].k,
			       bst_data[n_params*i+p].e);
		if (debugging("glayer"))
		    (void) fprintf(lstats->bst_out[p],"%10d",
				   bst_data[n_params*i+p].count);
		(void) fprintf(lstats->bst_out[p],"\n");
	    }
	    print_graph_footer(lstats->bst_out[p],tmpstr,about_to_stop);
	    fclose(lstats->bst_out[p]);

	    if (lstats->include_turb == YES)
	    {
	        (void) sprintf(fname,"%s.ts%s",lstats->turb_fname[p],
			right_flush(grid->step,TSTEP_FIELD_WIDTH));
		lstats->turb_out[p] = fopen(fname,"w");
		if (lstats->turb_out[p] == NULL)
		{
		    screen("\nERROR in print_layer_stats():\n");
		    screen("Unable to open file %s.\n",fname);
		    clean_up(ERROR);
		}
		print_machine_parameters(lstats->turb_out[p]);
						
		(void) sprintf(tmpstr,"FLUID %d TURBULENT LAYER",p);
		print_graph_header(lstats->turb_out[p],tmpstr,
				   lstats->turb_col_header,YES,front);
		for (i = 0; i < n_layers; ++i)
		{
		    h = lstats->h_min + i*lstats->dh;
		    (void) fprintf(lstats->turb_out[p],"%15.5e",h);
		    (void) fprintf(lstats->turb_out[p],"%15.5e",
				   turb_data[n_params*i+p].dd);
		    for (j = 0, dk = 0.0; j < dim; ++j)
		        dk += 0.5*turb_data[n_params*i+p].dvv[sym_index(j,j)];
		    (void) fprintf(lstats->turb_out[p],"%15.5e",dk);
		    (void) fprintf(lstats->turb_out[p],"%15.5e",
				   turb_data[n_params*i+p].de);
		    for (j = 0; j < dim; ++j)
		        (void) fprintf(lstats->turb_out[p],"%15.5e",
				       turb_data[n_params*i+p].dv[j]);
		    for (j = 0; j < dim; ++j)
		        (void) fprintf(lstats->turb_out[p],"%15.5e",
				       turb_data[n_params*i+p].dkv[j]);
		    for (j = 0; j < dim; ++j)
		        (void) fprintf(lstats->turb_out[p],"%15.5e",
				       turb_data[n_params*i+p].dev[j]);
		    for (j = 0; j < dim; ++j)
		        (void) fprintf(lstats->turb_out[p],"%15.5e",
				       turb_data[n_params*i+p].pv[j]);
		    for (j = 0; j < dim; ++j)
		        for (k = j; k < dim; ++k)
			    (void) fprintf(lstats->turb_out[p],"%15.5e",
			        turb_data[n_params*i+p].dvv[sym_index(j,k)]);
		    for (j = 0; j < dim; ++j)
		        for (k = j; k < dim; ++k)
			    (void) fprintf(lstats->turb_out[p],"%15.5e",
			        turb_data[n_params*i+p].vv[sym_index(j,k)]);
		    (void) fprintf(lstats->turb_out[p],"\n");
		}
		print_graph_footer(lstats->turb_out[p],tmpstr,about_to_stop);
		fclose(lstats->turb_out[p]);
	    }
	}

	stop_clock("print_layer_stats");
}		/*end record_layer_stats*/
	
#endif  /* defined(TWOD) || defined(THREED) */
