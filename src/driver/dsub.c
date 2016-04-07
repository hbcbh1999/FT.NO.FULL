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
*				dsub.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains subsidiary routines for the driver directory.
*/


#include <driver/ddecs.h>


	/* LOCAL Function Declarations */
LOCAL	Wave	*d_copy_wave(Wave*);
LOCAL	Front	*d_copy_front(Front*);
#if defined(USE_OVERTURE)
LOCAL   Front   *d_deep_copy_front(Front*);
#endif /* if defined(USE_OVERTURE) */

LOCAL	void	print_Grid_Stats_data_structure(Printplot*);
LOCAL	void	print_PRINT_CONTROL_structure(PRINT_CONTROL*);
LOCAL	void	print_OUTPUT_DATA_structure(OUTPUT_DATA*);
LOCAL	void	print_statistical_variables(int,double*,const char*);
LOCAL	int	do_cross_section(CROSS_SECTION*,Front*);
LOCAL	void	print_cross_section(FILE*,CROSS_SECTION*,Front*,Wave*);
LOCAL	boolean 	merge_equivalent_comp_list(INTERFACE*);
LOCAL 	void 	scatter_merged_components(int,COMPONENT*,Front*);
LOCAL 	void 	set_intfc_comp_to_newcomp(INTERFACE*,COMPONENT,COMPONENT);

#if defined(USE_HDF)
LOCAL	void	print_HDF_plot_data_structure(HDF_plot_data*);
#endif /* defined(USE_HDF) */


LOCAL	void print_statistical_variables(
	int	   n_var,
	double	   *stat_var,
	const char *name)
{
	int		var;

	if ((stat_var == NULL) || (n_var <= 0))
	{
	    (void) printf("statistical variable %s unft_assigned\n",name);
	    return;
	}
	(void) printf("\t\t\tstatistical variable %s\n\tvar\tvalue\n",name);
	for (var = 0;  var < n_var;  var++)
	    (void) printf("\t%d\t%g\n",var,stat_var[var]);
}		/*end print_statistical_variables*/






/****************************************************************************
*									    *
*			ddecs.h structures				    *
*									    *
****************************************************************************/

/*
*			print_PRINT_CONTROL_structure():
*/

LOCAL void print_PRINT_CONTROL_structure(
	PRINT_CONTROL	*ptc)
{
	PRINT_OPTIONS	*pto = &print_options(ptc);
	if (ptc == NULL)
	{
		(void) printf("PRINT_CONTROL structure NULL\n");
		return;
	}

	(void) printf("Begin PRINT_CONTROL structure %p\n",(POINTER)ptc);
	if (mesh_time_output(prt_mode(pto)))
	{
		(void) printf("prt_mode = %s\n","MESH_TIME");

		(void) printf("print_step_interval(mesh time mode) %d\n",
			      print_step_interval(pto));
		(void) printf("print_start_step (mesh time mode) %d\n",
			      print_start_step(pto));
		(void) printf("next printout time step %d\n",
			      next_print_step(ptc));
	}
	else
	{
		if (exact_time_output(prt_mode(pto)))
			(void) printf("prt_mode = %s\n","EXACT_TIME");
		else if (constant_time_output(prt_mode(pto)))
			(void) printf("prt_mode = %s\n","CONSTANT_TIME");

		(void) printf("print_interval(real time mode) %g\n",
			      print_time_interval(pto));
		(void) printf("print_start_time (real time mode) %g\n",
			      print_start_time(pto));
		(void) printf("next printout time %g\n",
			      next_print_time(ptc));
		(void) printf("print time tolerance %g\n",
			      print_time_tolerance(ptc));
	}
	(void) printf("End PRINT_CONTROL_structure %p\n",(POINTER)ptc);
}		/*end print_PRINT_CONTROL_structure*/

/*
*			print_OUTPUT_DATA_structure():
*/

LOCAL	void print_OUTPUT_DATA_structure(
	OUTPUT_DATA	*odata)
{
	if (odata == NULL)
	{
		(void) printf("OUTPUT_DATA structure NULL\n");
		return;
	}

	(void) printf("Begin OUTPUT_DATA structure %p\n",(POINTER)odata);
	(void) printf("prt = %p\n",(POINTER)Output_printplot(odata));
	print_PRINT_CONTROL_structure(Print_control(odata));
	(void) printf("directory name = %s\n",(Output_dir(odata) != NULL) ?
		      Output_dir(odata) : "NULL");
	(void) printf("file name = %s\n",(Output_filename(odata) != NULL) ?
		      Output_filename(odata) : "NULL");
	(void) printf("logfile = %p\n",(POINTER)Output_file(odata));
	(void) printf("binary output = %s\n",(Output_in_binary(odata) == YES) ?
		      "YES" : "NO");
	(void) printf("End OUTPUT_DATA structure %p\n",(POINTER)odata);
}		/*end print_OUTPUT_DATA_structure*/

/*
*			print_Printplot_structure()
*/

EXPORT	void print_Printplot_structure(
	Printplot	*prt)
{
	int		i;

	if (prt == NULL)
	{
	    (void) printf("\t\tPrintplot structure NULL\n");
	    return;
	}

	(void) printf("\n\n\n\t\tPrintplot structure %p\n\n",(POINTER)prt);
	for (i = 0; i < NUM_OUTPUT_FORMATS; i++)
	{
	    if (prt->main_output_data[i] != NULL)
	    {
	    	(void) printf("Printing controls for %s \n",
			      output_format_name(i));
	    	print_OUTPUT_DATA_structure(prt->main_output_data[i]);
	    }
	}

	print_Grid_Stats_data_structure(prt);

	(void) printf("\nnum rect state variables %d",prt->n_rect_state_vars);
	(void) printf("\toutput_soln %p\n",(POINTER)prt->output_soln);
	if ((prt->n_rect_state_vars > 0) && (prt->output_soln != NULL))
	{
	    for (i = 0;  i < prt->n_rect_state_vars;  i++)
	    	print_OUTPUT_SOLN_structure(prt->output_soln[i]);
	}
	(void) printf("\nread variables n_restart_vars %d",prt->n_restart_vars);
	(void) printf("\trestart_soln %p\n",(POINTER)prt->restart_soln);
	if ((prt->n_restart_vars > 0) && (prt->restart_soln != NULL))
	{
	    for (i = 0;  i < prt->n_restart_vars;  i++)
	    	print_INPUT_SOLN_structure(prt->restart_soln[i]);
	}

	(void) printf("\n\t\tnum tri plot variables %d\n",(int)prt->n_tri_vars);
	for (i = 0;  i < prt->n_tri_vars;  i++)
	    (void) printf("plot name: %-14.14s,  plot function %-14p\n",
	    	          prt->tri_plot_name[i],prt->tri_plot_function[i]);

#if defined(USE_HDF)
	(void) printf("\n\t\tnum HDF plot variables %d\n",prt->n_HDF_vars);
	for (i = 0;  i < prt->n_HDF_vars;  i++)
	    print_HDF_plot_data_structure(prt->HDF_data+i);
	(void) printf("\n\t\tnum SDS plot variables %d\n",prt->n_SDS_vars);
	for (i = 0;  i < prt->n_SDS_vars;  i++)
	    print_HDF_plot_data_structure(prt->SDS_data+i);
#endif /* defined(USE_HDF) */

	(void) printf("\n\t\tEnd Printplot structure %p\n\n",(POINTER)prt);
}		/*end print_Printplot_structure*/


#if defined(USE_HDF)
LOCAL void print_HDF_plot_data_structure(
	HDF_plot_data	*hdf_data)
{
	HDF_PRINT_OPTIONS	*opts = &HDF_print_opts(hdf_data);
	int var;
	(void) printf("\n\tHDF_plot_data structure %p\n",hdf_data);

	(void) printf("Number of plotting variables %d\n",hdf_num_vars(opts));
	(void) printf("dim %d\n",hdf_data->dim);
	print_int_vector("pixels ",hdf_pixels(opts),hdf_data->dim,"\n");
	print_general_vector("L0 ",hdf_L0(opts),hdf_data->dim,"\n");
	print_general_vector("U0 ",hdf_U0(opts),hdf_data->dim,"\n");
	print_general_vector("V ",hdf_V(opts),hdf_data->dim,"\n");
	for (var = 0; var < hdf_num_vars(opts); var++)
	{
	    HDF_FRAME_OPTS *fopts = &HDF_frame_opts(hdf_data->frame_data[var]);
	    (void) printf("plot_name[%d] %s, ",var,hdf_plot_name(fopts));
	    (void) printf("plot_function[%d] %p\n",var,
			  hdf_plot_function(fopts));
	    (void) printf("plot_filter[%d] %p\n",var,
			  hdf_plot_filter(fopts));
	    if (hdf_dflt_scale(fopts) == YES)
	    {
	    	(void) printf("Default scaling is on for %s\n",
			      hdf_plot_name(fopts));
	    }
	    else
	    {
	    	(void) printf("Default scaling is off for %s, ",
			      hdf_plot_name(fopts));
	    	(void) printf("min = %g, max = %g\n",
			      hdf_scale_min(fopts),hdf_scale_max(fopts));
	    }
	}

	(void) printf("\n\tEnd HDF_plot_data structure %p\n",
		      (POINTER)hdf_data);
}		/*end print_HDF_plot_data_structure*/
#endif /* defined(USE_HDF) */

LOCAL void print_Grid_Stats_data_structure(
	Printplot	*prt)
{
	Grid_Stats_data *gs_data = prt->gs_data;
	int		i;

	if (gs_data == NULL) return;

	(void) printf("\n\t\t\tStatistical Variables\n");
	(void) printf("nfloats %d\n",gs_data->nfloats);
	(void) printf("\tinitial_present %p\t\tpresent %p\n",
		      (POINTER)gs_data->initial_present,
		      (POINTER)gs_data->present);
	(void) printf("\tincremented %p\t\tsource_incremented %p\n",
		      (POINTER)gs_data->incremented,
		      (POINTER)gs_data->source_incremented);
	(void) printf("\tinhom_source_present %p\n",
		      (POINTER)gs_data->inhom_source_present);
	print_statistical_variables(gs_data->nfloats,
				    gs_data->initial_present,"initial_present");
	print_statistical_variables(gs_data->nfloats,
				    gs_data->present,"present");
	print_statistical_variables(gs_data->nfloats,
				    gs_data->incremented,"incremented");
	(void) printf("\tnum_point_sources %d\n",gs_data->num_point_sources);
	for( i = 0;  i < gs_data->num_point_sources;  i++ )
	{
		char	ss[35];

		(void) sprintf(ss,"for source[%d]",i);
		print_statistical_variables(gs_data->nfloats,
					    gs_data->source_incremented[i],ss);
	}
	print_statistical_variables(gs_data->nfloats,
				    gs_data->inhom_source_present,
				    "inhom_source");

	(void) printf("stat_var() %p\tstat_point_flux() %p\n",
		      gs_data->stat_var,
		      gs_data->stat_point_flux);
	(void) printf("stat_var() %p\tstat_point_flux() %p\n",
		      gs_data->stat_var,
		      gs_data->stat_point_flux);
	(void) printf("stat_flux() %p\n",gs_data->stat_flux);
	(void) printf("inhom_source() %p\n",
		      gs_data->stat_inhom_source);
}		/*end print_Grid_Stats_data_structure*/



/*
*			d_print_Grid_structure()
*/

EXPORT	void d_print_Grid_structure(
	Grid	*grid)
{
	(void) printf("\n\n\n\t\tGrid %p structure\n",(POINTER)grid);
	if( grid == NULL )
	{
	    (void) printf("\t\tstructure not yet allocated\n");
	    (void) printf("\n\t\tEnd Grid %p structure\n\n",(POINTER)grid);
	    return;
	}

	print_RECT_GRID_structure(grid->rect_grid);

	(void) printf("\t\t\tCurrent Time Grid\n");
	(void) printf("dt %-14g time %-14g step %-14d dt_lim %-14g\n",grid->dt,
		      grid->time,grid->step,dt_lim(grid));

	print_stopping_criteria(grid);

	(void) printf("\t\t\tPause Criteria\n");
	(void) printf("pause_mode %s\tpause_time %g\tpause_step %d\n",
		      mesh_time_output(pause_mode(grid))     ? "MESH_TIME"     :
		      exact_time_output(pause_mode(grid))    ? "EXACT_TIME"    :
		      constant_time_output(pause_mode(grid)) ? "CONSTANT_TIME" :
							       "INVALID VALUE",
		      pause_time(grid),pause_step(grid));
	
	(void) printf("is initialization complete %d\n",
		      grid->initialization_complete);

	(void) printf("\n\t\tEnd Grid %p structure\n\n",(POINTER)grid);
}		/*end d_print_Grid_structure*/


/*
*			output_format_name():
*
*	Converts an integer value for an output format into a descriptive
*	character string.
*/

EXPORT const char *output_format_name(
	int	format)
{
	const char *name;

	switch (format)
	{
	case PRT_FRONTS:
	    name = "front_plots";
	    break;
	case RECT_STATES:
	    name = "interior_states";
	    break;
	case TRI_STATES:
	    name = "tri_plots";
	    break;
	case HDF_STATES:
	    name = "HDF_plots";
	    break;
	case SDS_STATES:
	    name = "SDS_plots";
	    break;
	case VTK_STATES:
	    name = "VTK_plots";
	    break;
	case GD_MOVIE:
	    name = "GD_movie_plots";
	    break;
	default:
	    name = "unknown format";
	    break;
	}
	return name;
}		/*end output_format_name*/



/*ARGSUSED*/
EXPORT void print_cross_sectional_graphs(
	Grid		*grid,
	Wave		*wave,
	Front		*front,
	Printplot	*prt,
	OUTPUT_DATA	*data,
	boolean		about_to_stop)
{
	Cross_Sections_data *cr_data = Cr_Data(data);
	CROSS_SECTION	    *cr_sec;
	FILE		    *file = Output_file(data);

	if ((cr_data->cross_sections == NULL) ||
	    (front->_fprint_header_for_graph_curve_states == NULL) ||
	    (front->_fgraph_curve_states == NULL))
		return;
	(void) foutput(file);
	(void) fprintf(file,"\t\tCROSS SECTIONAL PLOTS OF STATES\n\n");
	for (cr_sec=cr_data->cross_sections; cr_sec!=NULL; cr_sec=cr_sec->next)
	    print_cross_section(file,cr_sec,front,wave);
	(void) fprintf(file,
		       "\n\n\t\tEND OF CROSS SECTIONAL PLOTS OF STATES\n\n");
}		/*end print_cross_sectional_graphs*/


LOCAL	void print_cross_section(
	FILE		*file,
	CROSS_SECTION	*cr_sec,
	Front		*front,
	Wave		*wave)
{
	BOND      *b;
	COMPONENT comp;
	CROSS	  *cross, *cr;
	CURVE	  *cur, *ctmp, *cc;
	INTERFACE *temp_intfc, *sav_intfc, *intfc;
	Locstate  sl, sr;
	NODE	  *ns, *ne;
	POINT	  *p;
	double	  length, s, cp;
	int	  i;
	int	  N = cr_sec->num_pts;

	if (!do_cross_section(cr_sec,front)) return;
	sav_intfc = current_interface();
	set_size_of_intfc_state(0);
	set_copy_intfc_states(NO);
	intfc = front->interf;
	temp_intfc = copy_interface(intfc);
	interpolate_intfc_states(temp_intfc) = NO;

	/* Copy cross section curve into temp_intfc */
	/* Only this curve will have states */

	size_of_state(temp_intfc) = front->sizest;
	ns = make_node(Point(cr_sec->pcoords[0]));
	ne = make_node(Point(cr_sec->pcoords[N-1]));
	cur = make_curve(NO_COMP,NO_COMP,ns,ne);
	for (i = N - 2; i > 0; i--)
	{
	    p = Point(cr_sec->pcoords[i]);
	    if (insert_point_in_bond(p,cur->first,cur) != FUNCTION_SUCCEEDED)
	    {
		screen("ERROR in print_cross_section(), "
		       "insert_point_in_bond() failed\n");
		clean_up(ERROR);
	    }
	}

	/* Assign states on cross section curve */

	b = cur->first;
	slsr(b->start,Hyper_surf_element(b),Hyper_surf(cur),&sl,&sr);
	comp = component(Coords(b->start),intfc);
	hyp_solution(Coords(b->start),comp,NULL,UNKNOWN_SIDE,
		     front,wave,sl,NULL);
	ft_assign(sr,sl,front->sizest);

	while (b != NULL)
	{
	    slsr(b->end,Hyper_surf_element(b),Hyper_surf(cur),&sl,&sr);
	    comp = component(Coords(b->end),intfc);
	    hyp_solution(Coords(b->end),comp,NULL,UNKNOWN_SIDE,
			 front,wave,sl,NULL);
	    ft_assign(sr,sl,front->sizest);
	    b = b->next;
	}
	if (intersections(temp_intfc,&cross,NO) == FUNCTION_FAILED) 
	{
	    (void) foutput(file);
	    screen("ERROR in print_cross_section(), "
		   "intersections() failed\n");
	    fprint_interface(file,intfc);
	    clean_up(ERROR);
	}
	for (cr = cross; cr != NULL; cr = cr->next)
	{
	    BOND *b, *btmp, *b2, *bb;
	    if (cr->c1 != cur && cr->c2 != cur) continue;
	    if (cr->c1 != cur)
	    {
	    	ctmp = cr->c1;
	    	btmp = cr->b1;
	    	cr->c1 = cr->c2;
	    	cr->b1 = cr->b2;
	    	cr->c2 = ctmp;
	    	cr->b2 = btmp;
	    }
	    if (cr->c2 == cur) continue;
	    vector_product_on_bonds(cr->b1,cr->b2,2,&cp);
	    b = cr->b1;
	    p = cr->p;
	    if (insert_point_in_bond(p,b,cur) != FUNCTION_SUCCEEDED)
	    {
		screen("ERROR in print_cross_section(), "
		       "insert_point_in_bond() failed\n");
		clean_up(ERROR);
	    }
	    rcl_after_insert_point(cr,p,b);
	    p = copy_point(cr->p);
	    if (insert_point_in_bond(p,b,cur) != FUNCTION_SUCCEEDED)
	    {
		screen("ERROR in print_cross_section(), "
		       "insert_point_in_bond() failed\n");
		clean_up(ERROR);
	    }
	    rcl_after_insert_point(cr,p,b);
	    s = separation(p,cr->b2->start,front->interf->dim) /
	    		   bond_length(cr->b2);
	    cc = correspond_curve(cr->c2);
	    for (b2 = cr->c2->first, bb = cc->first; b2 != cr->b2;
	    				b2 = b2->next, bb = bb->next)
	    	;
	    slsr(b->end,Hyper_surf_element(b),Hyper_surf(cur),&sl,&sr);
	    if (cp > 0.0)
	    	left_state_along_bond(s,bb,cc,sl);
	    else
	    	right_state_along_bond(s,bb,cc,sl);
	    ft_assign(sr,sl,front->sizest);
	    slsr(b->next->end,Hyper_surf_element(b),Hyper_surf(cur),&sl,&sr);
	    if (cp > 0.0)
	    	right_state_along_bond(s,bb,cc,sl);
	    else
	    	left_state_along_bond(s,bb,cc,sl);
	    ft_assign(sr,sl,front->sizest);
	}

	length = 0.0;
	fprint_header_for_graph_curve_states(file,front,cr_sec->message);
	fgraph_curve_states(file,cur,front,&length);
	(void) fprintf(file,"\nEnd %s\n\n",cr_sec->message);

	set_current_interface(sav_intfc);
	(void) delete_interface(temp_intfc);
}		/*end print_cross_section*/

/*ARGSUSED*/
LOCAL	int do_cross_section(
	CROSS_SECTION	*cr_sec,
	Front		*front)
{
	if (cr_sec == NULL) return NO;
	if (front->pp_grid->nn > 1)
	{
	    /* Make sure that this cross section lies solely within
	     *  this subdomain */
					 
	    double   *VL = front->rect_grid->VL;
	    double   *VU = front->rect_grid->VU;
	    int	N = cr_sec->num_pts;
	    int     i, j, dim = front->interf->dim;
	    for (i = 0; i < dim; i++)
	    {
		for (j = 0; j < N; j++)
		    if (!Between(cr_sec->pcoords[j][i],VL[i],VU[i]))
		        return NO;
	    }
	}
	return YES;
}		/*end do_cross_section*/

/*ARGSUSED*/
EXPORT	boolean regrid_front_and_wave(
	Front		*front,
	Wave		*wave,
	RECT_GRID	*fgr,
	RECT_GRID	*tgr,
	boolean		restart_init)
{
	COMPONENT	comp;
	Front		*tempfront;
	Locstate	state;
	RECT_GRID	*sav_gr = front->rect_grid;
	Wave		Tempwave;
	double		*coords;
	boolean		sav_interp = interpolate_intfc_states(front->interf);
	int		icoords[3];
	int		ix, xmax = fgr->gmax[0];
	double		dt_frac = 1.0;
	int		iy, ymax = (fgr->dim > 1) ? fgr->gmax[2] : 1;
	int		iz, zmax = (fgr->dim == 3) ? fgr->gmax[2] : 1;

	tempfront = copy_front(front);
	tempfront->rect_grid = fgr;
	set_size_of_intfc_state(front->sizest);
	set_copy_intfc_states(YES);
	set_add_to_correspond_list(YES);
	tempfront->interf = copy_interface(front->interf);
	if (tempfront->interf == NULL)
	{
	    screen("ERROR in regrid_front_and_wave(), "
	           "unable to copy interface\n");
	    free_front(tempfront);
	    return NO;
	}
	set_topological_grid(tempfront->interf,tgr);
	set_computational_grid(tempfront->interf,fgr);
	interpolate_intfc_states(tempfront->interf) = YES;

	Redistribution_count(tempfront) = 0;
	tempfront->dt_frac = &dt_frac;
	if (redistribute(tempfront,YES,restart_init) != GOOD_REDISTRIBUTION)
	{
	    screen("ERROR in regrid_front_and_wave(), "
	           "redistribute() failed\n");
	    free_front(tempfront);
	    return NO;
	}

	assign_wave_parameters(&Tempwave,wave);
	Tempwave.rect_grid = fgr;
	start_clock("init_hyp_solution");
	if (init_hyp_solution_function(&Tempwave,tempfront) != GOOD_STEP)
	{
	    screen("ERROR in regrid_front_and_wave(), "
	           "init_hyp_solution_function() failed\n");
	    free_front(tempfront);
	    return NO;
	}
	stop_clock("init_hyp_solution");

	for (iz = 0; iz < zmax; iz++)
	{
	    icoords[2] = iz;
	    for (iy = 0; iy < ymax; iy++)
	    {
	    	icoords[1] = iy;
	    	for (ix = 0; ix < xmax; ix++) 
	    	{
	    	    icoords[0] = ix;
	    	    coords = Rect_coords(icoords,&Tempwave);
	    	    comp = Rect_comp(icoords,&Tempwave);
	    	    state = Rect_state(icoords,&Tempwave);
	    	    hyp_solution(coords,comp,NULL,UNKNOWN_SIDE,
				 front,wave,state,NULL);
	    	}
	    }
	}

	*sav_gr = *fgr; /* Copy new grid info to original rect grid */
	assign_wave_pointers(wave,&Tempwave);
	assign_interface_and_free_front(front,tempfront);
	interpolate_intfc_states(front->interf) = sav_interp;

	/* Reset rect grid to original structures */

	front->rect_grid = wave->rect_grid = sav_gr;
	set_computational_grid(front->interf,front->rect_grid);

	return YES;
}		/*end regrid_front_and_wave*/


#define		MAX_MERGE_COMP		20

EXPORT	void delete_untracked_hyper_surfaces(
	Front		*front,
	Wave		*wave)
{
	HYPER_SURF	**hs;
	Front		*tempfront;
	INTERFACE	*intfc = front->interf;
	size_t		sizest;
	boolean		sav_interp;
	boolean		untrack = NO;
	int		num_mc;	
	COMPONENT	merge_comps[MAX_MERGE_COMP*3];

	debug_print("untrack","Entered delete_untracked_hyper_surfaces()\n");

	for (hs = intfc->hss; hs && *hs; hs++)
	{
	    if (untracked_hyper_surf(*hs) == YES)
	    {
	    	untrack = YES;
	    	break;
	    }
	}

	if (pp_max_status(untrack) == NO)
	{
	    if (merge_equivalent_comp_list(intfc))
	    {
		reinit_hyp_solution_function(wave,front);
	    }
	    debug_print("untrack","Left delete_untracked_hyper_surfaces()\n");
	    return;
	}

	tempfront = copy_front(front);
	sizest = front->sizest;
	set_size_of_intfc_state(sizest);
	set_copy_intfc_states(YES);
	set_add_to_correspond_list(YES);
	tempfront->interf = copy_interface(front->interf);
	if (tempfront->interf == NULL)
	{
	    screen("ERROR in delete_untracked_hypersurfaces(), "
	           "unable to copy interface\n");
	    free_front(tempfront);
	    clean_up(ERROR);
	}
	sav_interp = interpolate_intfc_states(front->interf);
	interpolate_intfc_states(tempfront->interf) = YES;

	intfc = tempfront->interf;
	num_mc = 0;
	switch (intfc->dim)
	{
	case 1:
	    {
	        POINT	**p;

	        if (debugging("untrack"))
		{
		    (void) printf("Deleting selected points from intfc %p\n",
				  intfc);
		    (void) printf("Points to be deleted\n");
	            for (p = intfc->points; p && *p; p++)
		    {
	                if (untracked_hyper_surf(*p) == YES)
			    print_point(*p);
		    }
		    print_interface(intfc);
		}
	        for (p = intfc->points; p && *p; p++)
	        {
	            if (untracked_hyper_surf(*p) == YES)
	            {
			merge_comps[3*num_mc] = negative_component(*p);
			merge_comps[3*num_mc+1] = positive_component(*p);
			merge_comps[3*num_mc+2] = negative_component(*p);
			num_mc++;
	                if (untrack_point(*p,negative_component(*p),
					  tempfront) == YES)
	            	    p = intfc->points - 1;
	            }
	        }
	    }
	    break;
	case 2:
	    {
	    	CURVE		**c;
	    	O_CURVE		oc;
	    	UNTRACK_FLAG	flag;

	    	flag.start_states_set = YES;
	    	flag.end_states_set = YES;

	    	for (c = intfc->curves; c && *c; c++)
	    	{
	    	    if (untracked_hyper_surf(*c) == YES)
	    	    {
			merge_comps[3*num_mc] = negative_component(*c);
			merge_comps[3*num_mc+1] = positive_component(*c);
			merge_comps[3*num_mc+2] = negative_component(*c);
			num_mc++;
	    	    	oc.curve = *c;
	    	    	oc.orient = POSITIVE_ORIENTATION;
	    	    	if (untrack_curve(&oc,NULL,negative_component(*c),0.0,
	    	    		          tempfront,NULL,NULL,flag) == YES)
	    	    	    c = intfc->curves - 1;
	    	    }
	    	}
	    }
	    break;
	case 3:
	    {
		SURFACE  **s;
		for (s = intfc->surfaces; s && *s; s++)
		{
	    	    if (untracked_hyper_surf(*s) == YES)
	    	    {
			merge_comps[3*num_mc] = negative_component(*s);
			merge_comps[3*num_mc+1] = positive_component(*s);
			merge_comps[3*num_mc+2] = negative_component(*s);
			num_mc++;
	    	    	if (untrack_surface(*s,negative_component(*s),
					    tempfront) == YES)
	    	    	    s = intfc->surfaces - 1;
	    	    }
		}
	    }
	    break;
	}
	scatter_merged_components(num_mc,merge_comps,tempfront);

	assign_interface_and_free_front(front,tempfront);
	interpolate_intfc_states(front->interf) = sav_interp;

	reinit_hyp_solution_function(wave,front);

	debug_print("untrack","Left delete_untracked_hyper_surfaces()\n");
}		/*end delete_untracked_hyper_surfaces*/

LOCAL	boolean merge_equivalent_comp_list(
	INTERFACE *intfc)
{
	EQUIV_COMPS     *e_comps;
	int i;
	COMPONENT	newcomp;
	POINT **pp;
	CURVE **cc;
	SURFACE **ss;

	if (E_comps(intfc) == NULL) return NO;

	for (e_comps = E_comps(intfc); e_comps; e_comps = e_comps->next)
	{
	    newcomp = e_comps->comp[0];
	    switch(intfc->dim)
	    {
	    case 1:
	    	for (pp = intfc->points; pp && *pp; ++pp)
		{
		    for (i = 1; i < e_comps->n_equiv; ++i)
		    {
		    	if (positive_component(*pp) == e_comps->comp[i])
			    positive_component(*pp) = e_comps->comp[0];
		    	if (negative_component(*pp) == e_comps->comp[i])
			    negative_component(*pp) = e_comps->comp[0];
		    }
		}
	    	break;
	    case 2:
	    	for (cc = intfc->curves; cc && *cc; ++cc)
		{
		    for (i = 1; i < e_comps->n_equiv; ++i)
		    {
		    	if (positive_component(*cc) == e_comps->comp[i])
			    positive_component(*cc) = e_comps->comp[0];
		    	if (negative_component(*cc) == e_comps->comp[i])
			    negative_component(*cc) = e_comps->comp[0];
		    }
		}
	    	break;
	    case 3:
	    	for (ss = intfc->surfaces; ss && *ss; ++ss)
		{
		    for (i = 1; i < e_comps->n_equiv; ++i)
		    {
		    	if (positive_component(*ss) == e_comps->comp[i])
			    positive_component(*ss) = e_comps->comp[0];
		    	if (negative_component(*ss) == e_comps->comp[i])
			    negative_component(*ss) = e_comps->comp[0];
		    }
		}
	    }
	}
	return YES;
}	/* end merge_equivalent_comp_list */

EXPORT	void	d_set_default_front_parameters(
	INIT_DATA	*init,
	Front		*fr)
{
	h_set_default_front_parameters(init,fr);
	fr->_copy_front =			d_copy_front;
	fr->_copy_into_front =			d_copy_into_front;
	fr->_print_Front_structure =		d_print_D_Front_structure;

#if defined(USE_OVERTURE)
        fr->_deep_copy_front =                  d_deep_copy_front;
#endif /*if defined(USE_OVERTURE) */
}		/*end d_set_default_front_parameters*/

/*
*			d_copy_front():
*
*	Basic default function for copying a front structure.
*	Allocates storage for the new front and copies the
*	argument into the new structure.
*/

LOCAL	Front *d_copy_front(
	Front		*fr)
{
	Front		*newfr;

	scalar(&newfr,sizeof(D_Front));
	copy_into_front(newfr,fr);
	return newfr;
}		/*end d_copy_front*/


#if defined(USE_OVERTURE)
/*
*                       d_deep_copy_front():
*
*       Basic default function for copying a front structure.
*       Allocates storage for the new front and copies the
*       argument into the new structure.
*/
LOCAL   Front *d_deep_copy_front(
        Front           *fr)
{
        Front           *newfr;

        scalar(&newfr,sizeof(D_Front));
        copy_into_front(newfr,fr);
        scalar(&(newfr->rect_grid), sizeof(RECT_GRID));
        *(newfr->rect_grid) = *(fr->rect_grid);
        scalar(&(newfr->pd_flag),sizeof(Patch_bdry_flag));
        return newfr;
}               /*end d_deep_copy_front*/
#endif /* if defined(USE_OVERTURE) */

/*
*			d_copy_into_front():
*
*	Copies fr into newfr.  Assumes newfr is already allocated.
*/

EXPORT	void d_copy_into_front(
	Front		*newfr,
	Front		*fr)
{
	h_copy_into_front(newfr,fr);
	DriverFrontExtension(newfr) = DriverFrontExtension(fr);
}		/*end d_copy_into_front*/


/* 
*			d_set_default_wave_parameters():
*	
*	Initialize function pointers for D_ extensions to Wave structure.
*/

EXPORT	void	d_set_default_wave_parameters(
	INIT_DATA *init,
	Wave	  *wave)
{
	h_set_default_wave_parameters(init,wave);

	switch (wave->rect_grid->dim)
	{
	case 1:
		plot_states_function(wave) = plot_states1d;
		break;
	case 2:
		plot_states_function(wave) = plot_states2d;
		break;
	case 3:
		plot_states_function(wave) = plot_states3d;
		break;
	}

	/* Currently only one driver extension. */
	wave->_copy_wave = d_copy_wave;
	wave->_copy_into_wave = d_copy_into_wave;

}		/*end d_set_default_wave_parameters*/

LOCAL	Wave	*d_copy_wave(
	Wave		*wave)
{
	Wave		*newwave;

	scalar(&newwave,sizeof(D_Wave));
	copy_into_wave(newwave,wave);
	return newwave;
}		/*end d_copy_wave*/

EXPORT	void	d_copy_into_wave(
	Wave	*newwave,
	Wave	*wave)
{
	h_copy_into_wave(newwave,wave);
	DriverWaveExtension(newwave) = DriverWaveExtension(wave);
}		/*end d_copy_into_wave*/

EXPORT	void	d_assign_wave_parameters(
	Wave	*newwave,
	Wave	*wave)
{
	h_assign_wave_parameters(newwave,wave);
	plot_states_function(newwave) = plot_states_function(wave);
	chart_of_wave(newwave) = chart_of_wave(wave);
}		/*end d_assign_wave_parameters*/


LOCAL void scatter_merged_components(
	int num_mc,
	COMPONENT *merge_comps,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	int myid = pp_mynode();
	int i,j,nc;
	COMPONENT comp1,comp2,new_comp;

	if (pp_numnodes() == 1)
	    return;
	for (i = 0; i < pp_numnodes(); ++i)
	{
	    if (i == myid) continue;
	    pp_send_all(myid,&num_mc,INT);
	    pp_send_all(myid,merge_comps,num_mc*3*sizeof(COMPONENT));
	}
	for (i = 0; i < pp_numnodes(); ++i)
	{
	    if (i == myid) continue;
	    pp_recv(i,i,&nc,INT);
	    pp_recv(i,i,merge_comps,nc*3*sizeof(COMPONENT));
	    for (j = 0; j < nc; ++j)
	    {
	    	comp1 = merge_comps[3*j];
		comp2 = merge_comps[3*j+1];
		new_comp = merge_comps[3*j+2];
		set_equivalent_comps(new_comp,comp1,intfc);
		set_equivalent_comps(new_comp,comp2,intfc);
	    }
	}
	merge_equivalent_comp_list(intfc);
}	/* end scatter_merged_components */

